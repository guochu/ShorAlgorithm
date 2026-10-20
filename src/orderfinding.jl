# =============================================================================
# orderfinding.jl — 量子求阶（Shor 的量子核）
#
# 移植自 VQC `src/utility/shor/orderfinding.jl` 与 `orderfinding_sqc.jl`。
#
# 四种算法（`QuantumOrderFindingAlgorithm` 的子类型单例，n = bits(N)）：
#   * `SemiclassicalQFT`：单控制比特 + 半经典（测量反馈）QFT，2n+3 个
#     qubit（默认，Beauregard 布局）；打包线路含测量与经典反馈（IfOp）；
#   * `SemiclassicalLNN`：单控制比特 + 半经典 QFT，Fowler LNN 布局，
#     2n+4 个 qubit；打包线路同上；
#   * `FullQFT`：完整 QPE，Beauregard 电路，4n+2 个 qubit；
#   * `LNNQFT`：完整 QPE，Fowler LNN 电路，4n+3 个 qubit。
#
#   注：两个半经典版的 IPE 控制位（Beauregard 的 t / LNN 的 ki）都是
#   干净的控制位——LNN 的 ki 虽在乘法块内部被路由，但每轮结束都会
#   复原（实证：ki=0 分支恒等、ki=1 分支恰为 ×a）。
#
# 两条使用路径：
#   * `order_finding(a, N; alg)`：驱动函数，内部模拟 + 测量，返回 y ≈ s/r；
#   * `order_finding_circuit(a, N; alg)`：把整个计算打包成一条线路
#     （半经典版含测量 / 经典反馈 / 输出寄存器），作用到量子态上即得解。
#
# 返回值 / 读出值 `y ∈ [0,1)` 为相位估计读数：y ≈ s/r（r 为 a 模 N 的阶）。
# =============================================================================

# ── 算法类型 ─────────────────────────────────────────────────────────────────

"""
    QuantumOrderFindingAlgorithm

求阶算法的抽象基类。具体算法为携带配置（`verbosity`）的单例子类型，
经 `order_finding` / `order_finding_circuit` 的 `alg` 关键字派发；
旧接口的 `less_qubits` / `algorithm` 组合已折叠进算法类型本身：

* `SemiclassicalQFT()`：`less_qubits=true,  algorithm=:qft`（默认）；
* `SemiclassicalLNN()`：`less_qubits=true,  algorithm=:lnn`；
* `FullQFT()`：`less_qubits=false, algorithm=:qft`；
* `LNNQFT()`：`less_qubits=false, algorithm=:lnn`。
"""
abstract type QuantumOrderFindingAlgorithm end

"""
    SemiclassicalQFT(; verbosity=0)

单控制比特 + 半经典 QFT（Griffiths–Niu 测量反馈，Beauregard 单控制
技巧）：仅 `2n+3` 个 qubit，**默认算法**。
"""
struct SemiclassicalQFT <: QuantumOrderFindingAlgorithm
    verbosity::Int
    SemiclassicalQFT(; verbosity::Int = 0) = new(verbosity)
end

"""
    FullQFT(; verbosity=0)

完整 QPE + Beauregard 全连接电路：`4n+2` 个 qubit。
"""
struct FullQFT <: QuantumOrderFindingAlgorithm
    verbosity::Int
    FullQFT(; verbosity::Int = 0) = new(verbosity)
end

"""
    LNNQFT(; verbosity=0)

完整 QPE + Fowler–Devitt–Hollenberg 线性近邻电路：`4n+3` 个 qubit。
"""
struct LNNQFT <: QuantumOrderFindingAlgorithm
    verbosity::Int
    LNNQFT(; verbosity::Int = 0) = new(verbosity)
end

"""
    SemiclassicalLNN(; verbosity=0)

单控制比特 + 半经典 QFT（Griffiths–Niu 测量反馈，Fowler LNN 布局）：
仅 `2n+4` 个 qubit。IPE 控制位为 LNN 布局的 `ki`（qubit `n+2`），
win（累积寄存器）初值 `|1⟩` 置于 qubit `n+1`（LSB），b 工作寄存器
以 Fourier 零态 `|0̃⟩` 起步。
"""
struct SemiclassicalLNN <: QuantumOrderFindingAlgorithm
    verbosity::Int
    SemiclassicalLNN(; verbosity::Int = 0) = new(verbosity)
end

# ── FullQFT：Beauregard 完整 QPE ─────────────────────────────────────────────

"""
    _beauregard_qpe_circuit(a, N) -> Circuit

完整相位估计求阶线路（Beauregard，4n+2 个 qubit，n = bits(N)，
寄存器位宽 L = n+1）：

    win  = 2n+2        累计寄存器 LSB，初值 |1⟩
    计数  = 2n+3..4n+2  m = 2n 位计数寄存器，初态 |+…+⟩

线路 = `X(2n+2)`、计数寄存器 H、`modular_exponentiator_circuit(N, a, 2n)`、
计数寄存器逆 QFT。测量计数寄存器得 `y ≈ s/r`。
对应旧 `order_finding_qft_circuit`（旧顶层 `shift(+1)` 已内建）。
"""
function _beauregard_qpe_circuit(a::Integer, N::Integer)
    _check_register_args(N, a)
    n = ndigits(N, base=2)
    m = 2n
    c = Circuit(4n + 2)
    push!(c, X(2n + 2))                        # win LSB = 1
    for q in (2n+3):(4n+2)
        push!(c, H(q))
    end
    append!(c, modular_exponentiator_circuit(N, a, m))
    push!(c, block(sdagger(qft_circuit(collect((2n+3):(4n+2)))); name = :IQFT))
    return c
end

# ── LNNQFT：Fowler 完整 QPE ──────────────────────────────────────────────────

"""
    _lnn_qpe_circuit(a, N) -> Circuit

完整相位估计求阶线路的 LNN（线性近邻）版本（Fowler，4n+3 个 qubit，
n = bits(N)）：

    ki   = 3n+1          初值 |1⟩（经首次 SWAP 链注入 kx）
    计数  = 1..2n        m = 2n 位计数寄存器，初态 |+…+⟩

线路 = `X(3n+1)`、计数寄存器 H、`lnn_modular_exponentiator_circuit(N, a, 2n)`、
计数寄存器逆 QFT。对应旧 `order_finding_lnn_circuit`。
"""
function _lnn_qpe_circuit(a::Integer, N::Integer)
    _check_register_args(N, a)
    n = ndigits(N, base=2)
    m = 2n
    c = Circuit(m + 2n + 3; cregs = CReg[])
    push!(c, X(m + n + 1))                     # ki = |1⟩（旧 gate(3n, X)，m = 2n）
    for q in 1:m
        push!(c, H(q))
    end
    append!(c, lnn_modular_exponentiator_circuit(N, a, m))
    push!(c, block(sdagger(qft_circuit(m)); name = :IQFT))   # 计数寄存器 = qubit 1..m
    return c
end

# ── SemiclassicalQFT：单控制比特 + 测量反馈 ──────────────────────────────────

# ── SemiclassicalLNN：LNN 布局的测量反馈 ─────────────────────────────────────

# ── 打包：把整个求阶计算包装成一条线路 ───────────────────────────────────────

"""a^(2^(k-1)) mod N 及其逆的前 m 个幂。"""
function _mod_powers(a::Int, N::Int, m::Int)
    astore = Vector{Int}(undef, m)
    rastore = Vector{Int}(undef, m)
    a_cur = a
    ra_cur = invmod(a, N)
    for k in 1:m
        astore[k] = a_cur
        rastore[k] = ra_cur
        a_cur = mod(a_cur * a_cur, N)
        ra_cur = mod(ra_cur * ra_cur, N)
    end
    return astore, rastore
end

"""
    _semiclassical_qpe_circuit(a, N) -> Circuit

半经典求阶的**单线路打包**：把“H → 受控 U → 经典反馈相位 → H → 测量 →
复位”的整个迭代显式编码进 IR——测量写回经典寄存器 `r`（`MeasOp`），
反馈相位按 Griffiths–Niu **逐比特线性展开**为 `IfOp` 条件相位：第 k
轮对每个已测比特 r_j 施加条件 `PHASE(-π/2^(k-j))`（最新比特权重 1/2，
与 `order_finding_scq` 的 `PHASE(-π·0.r_{k-1}…r_1)` 等价），共
m(m-1)/2 个 IfOp；顺带的条件 `X` 把刚学到的结果拷贝到**输出寄存器**。

重复子线路封装为命名块（`draw(c; unroll=false)` 以块为单元）：
b 寄存器 `QFT_b`、每轮受控 U 幂 `U(a·)`。

布局（4n+3 个 qubit，m = 2n，n = bits(N)）：t=1，b=2..n+2，
win=n+3..2n+2，ctl=2n+3，out=2n+4..4n+3。out 的第 k 位存第 k 轮
测量比特（测量自 LSB 向 MSB 进行，out[1] 即读出值的 LSB）；以
out[1] 为 LSB 读出整数值 / 2^m ≈ s/r。实际运行请优先使用
`order_finding` 驱动版本（免于逐分支求值的经典开销）。
"""
function _semiclassical_qpe_circuit(a::Integer, N::Integer)
    _check_register_args(N, a)
    n = ndigits(N, base=2)
    m = 2n
    ctl = 2n + 3
    out0 = 2n + 3                       # out 寄存器起始：qubit out0+k 存第 k 轮结果
    rreg = CReg(:r, m)
    c = Circuit(4n + 3; cregs=[rreg])
    qs_b = collect(2:(n+2))
    astore, rastore = _mod_powers(a, N, m)
    push!(c, X(2n + 2))                 # win LSB = 1
    push!(c, block(qft_circuit(qs_b); name = :QFT_b))   # b 进入 Fourier 基
    for k in 1:m
        push!(c, H(ctl))
        push!(c, block(modular_u_circuit(N, astore[m-k+1], rastore[m-k+1]);
                       name = Symbol("U(", astore[m-k+1], ")")))
        if k >= 2
            # 反馈角 = -π·0.r_{k-1}…r_1（最新比特权重 1/2），逐比特分解
            # 为条件相位之和；同一轮顺带把第 k-1 轮结果拷到输出寄存器。
            for j in 1:(k-1)
                if_then(c, rreg[m-j+1] == 1,
                        Circuit([PHASE(-π / 2.0^(k-j), ctl)]))
            end
            if_then(c, rreg[m-k+2] == 1, Circuit([X(out0 + k - 1)]))
        end
        push!(c, H(ctl))
        push!(c, measure(ctl, rreg[m-k+1]))
        push!(c, reset(ctl))
    end
    # 最后一轮 r_m 无后继反馈分支，单独按取值拷贝（r_m = rreg[1]）
    if_then(c, rreg[1] == 1, Circuit([X(out0 + m)]))
    return c
end

"""
    _semiclassical_lnn_qpe_circuit(a, N) -> Circuit

`SemiclassicalLNN` 的**单线路打包**：结构与 `_semiclassical_qpe_circuit`
完全同构，只是把 Beauregard 布局换成 Fowler LNN 布局——

    win = 1..n+1（LSB = n+1，初值 |1⟩），ki = n+2（IPE 控制位），
    kx = n+3，b = n+4..2n+4（Fourier 工作寄存器，初态 |0̃⟩ 经 H 制备），
    out = 2n+5..2n+4+m。

每轮：H(ki) → `U(a·)` 块（受控模乘，内部自带 b 的 QFT 括号）→
逐比特反馈 IfOp（条件相位 + 条件 `X` 拷贝到 out）→ H(ki) →
测量写 `r` → 复位 ki。以 out[1] 为 LSB 读出整数值 / 2^m ≈ s/r。
"""
function _semiclassical_lnn_qpe_circuit(a::Integer, N::Integer)
    _check_register_args(N, a)
    n = ndigits(N, base = 2)
    m = 2n
    ki = n + 2
    wlsb = n + 1
    out0 = 2n + 5                       # out 寄存器：qubit out0+k-1 存第 k 轮结果
    rreg = CReg(:r, m)
    c = Circuit(2n + 4 + m; cregs = [rreg])
    astore, rastore = _mod_powers(a, N, m)
    push!(c, X(wlsb))                   # win = |1⟩
    for q in (n+4):(2n+4)               # b = Fourier 零态 |0̃⟩（kx 置 |0⟩）
        push!(c, H(q))
    end
    for k in 1:m
        push!(c, H(ki))
        push!(c, block(lnn_modular_u_circuit(N, astore[m-k+1], rastore[m-k+1]);
                       name = Symbol("U(", astore[m-k+1], ")")))
        if k >= 2
            # 反馈角 = -π·0.r_{k-1}…r_1（最新比特权重 1/2），逐比特分解；
            # 同一轮顺带把第 k-1 轮结果拷到输出寄存器。
            for j in 1:(k-1)
                if_then(c, rreg[m-j+1] == 1,
                        Circuit([PHASE(-π / 2.0^(k-j), ki)]))
            end
            if_then(c, rreg[m-k+2] == 1, Circuit([X(out0 + k - 2)]))
        end
        push!(c, H(ki))
        push!(c, measure(ki, rreg[m-k+1]))
        push!(c, reset(ki))
    end
    # 最后一轮 r_m 无后继反馈分支，单独按取值拷贝（r_m = rreg[1]）
    if_then(c, rreg[1] == 1, Circuit([X(out0 + m - 1)]))
    return c
end

# ── 统一入口（后端注入）─────────────────────────────────────────────────────

"相位读数寄存器（按测量顺序，`qs[1]` 权重 1/2）对应的 qubit 列表。"
_readout_qubits(alg::QuantumOrderFindingAlgorithm, n::Int, m::Int) =
    _readout_qubits(typeof(alg), n, m)
_readout_qubits(::Type{SemiclassicalQFT}, n::Int, m::Int) = collect((2n+4):(4n+3))
_readout_qubits(::Type{SemiclassicalLNN}, n::Int, m::Int) = collect((2n+5):(4n+4))
_readout_qubits(::Type{FullQFT},          n::Int, m::Int) = collect((2n+3):(4n+2))
_readout_qubits(::Type{LNNQFT},           n::Int, m::Int) = collect(1:m)

"""给线路追加相位读出测量：寄存器 `qs` 逐位测入经典寄存器 `:ph`。"""
function _with_phase_readout!(c::Circuit, qs::Vector{Int})
    ph = CReg(:ph, length(qs))
    push!(c.cregs, ph)
    for (k, q) in enumerate(qs)
        push!(c, measure(q, ph[k]))
    end
    return c
end

"""从计数键提取 `:ph` 相位读数样本：`y = 整数值 / 2^m`。"""
function _phase_samples(counts::Dict{String,Int}, m::Int)
    samples = Tuple{Float64,Int}[]                  # (y, 出现次数)
    for (key, n) in counts
        for part in split(key, ',')
            name, bits = split(part, ':')
            if name == "ph"
                push!(samples, (parse(Int, bits, base = 2) / 2.0^m, n))
            end
        end
    end
    isempty(samples) && throw(ArgumentError("phase register :ph not found in counts"))
    return samples
end

"""
    order_finding(a, N; alg=SemiclassicalQFT(), backend, shots=1) -> Float64

求阶唯一入口，按 `alg` 派发到具体实现：

* `SemiclassicalQFT()`：单控制比特半经典 QFT（2n+3 qubit，默认）；
* `SemiclassicalLNN()`：单控制比特半经典 QFT 的 LNN 版（2n+4 qubit）；
* `FullQFT()`：完整 QPE，Beauregard 电路（4n+2 qubit）；
* `LNNQFT()`：完整 QPE，Fowler 近邻电路（4n+3 qubit）。

执行经 `QuantumCircuits.simulate(c, backend; shots)` 完成——**backend
为必填的 `QuantumCircuits.Interface.Backend` 实例**（算法包与具体
模拟器解耦；如 VQC 的 `StateVectorBackend()`）。运行 `shots` 个样本，
返回出现最频繁的相位读数 `y ≈ s/r`（r 为 a 模 N 的阶）；用
`continued_fraction` 从 `y` 提取阶 `r`。

半经典算法依赖中途测量与经典反馈，需后端支持 `:mid_measure` 能力。
"""
function order_finding(a::Int, N::Int;
                       alg::QuantumOrderFindingAlgorithm = SemiclassicalQFT(),
                       backend::QuantumCircuits.Interface.Backend,
                       shots::Int = 1)
    _check_order_finding_input(a, N)
    shots >= 1 || throw(ArgumentError("shots must be >= 1"))
    if alg isa Union{SemiclassicalQFT, SemiclassicalLNN}
        supports(backend, :mid_measure) ||
            throw(ArgumentError("$(typeof(backend)) does not support :mid_measure, " *
                                "required by $alg"))
    end

    n = ndigits(N, base = 2)
    m = 2n
    c = _with_phase_readout!(order_finding_circuit(a, N; alg = alg),
                             _readout_qubits(alg, n, m))
    r = simulate(c, backend; shots = shots)
    samples = _phase_samples(r.counts, m)
    y, _ = argmax(s -> s[2], samples)               # 出现最频繁的读数
    (alg.verbosity >= 1) && println("order finding ($alg): y = $y")
    return y
end

"以新的 `verbosity` 复制算法配置。"
_with_verbosity(alg::QuantumOrderFindingAlgorithm, verbosity::Int) =
    typeof(alg)(; verbosity=verbosity)

"""
    order_finding_circuit(a, N; alg=SemiclassicalQFT()) -> Circuit

把整个求阶计算打包成**一条线路**，作用到 `|0…0⟩` 态上并（对半经典版）
用 `simulate` 执行即可得到解：

* `SemiclassicalQFT()`：含测量 / 经典反馈（IfOp）/ 复位的动态线路；
  解编码在输出寄存器 out = 2n+4..4n+3（读值 / 2^m ≈ s/r）；
* `SemiclassicalLNN()`：同上的 LNN 版（2n+4+m qubit），out = 2n+5..2n+4+m；
* `FullQFT()`：纯酉线路，解在计数寄存器 2n+3..4n+2 的测量分布中；
* `LNNQFT()`：同上的 LNN 版，计数寄存器 1..2n。
"""
function order_finding_circuit(a::Integer, N::Integer;
                               alg::QuantumOrderFindingAlgorithm = SemiclassicalQFT())
    return _order_finding_circuit(alg, a, N)
end

_order_finding_circuit(::SemiclassicalQFT, a::Integer, N::Integer) =
    _semiclassical_qpe_circuit(a, N)
_order_finding_circuit(::SemiclassicalLNN, a::Integer, N::Integer) =
    _semiclassical_lnn_qpe_circuit(a, N)
_order_finding_circuit(::FullQFT, a::Integer, N::Integer) =
    _beauregard_qpe_circuit(a, N)
_order_finding_circuit(::LNNQFT, a::Integer, N::Integer) =
    _lnn_qpe_circuit(a, N)
