# =============================================================================
# fowler.jl — Fowler–Devitt–Hollenberg 型 LNN（线性近邻）模算术电路
#
# 移植自 VQC `src/utility/shor/Fowler.jl`（0-based → 1-based；
# 原代码中的 `shift(·, k)` 在本实现中保留为同值 `shift`，因为各积木
# 已按旧书写位置 +1 输出）。
# 寄存器位宽统一取 `L = bits(N) + 1`（`_register_bits`）。
#
# 寄存器布局（1-based，canonical 位置）：
#
#   LNN 模加法器 `lnn_modular_adder_circuit`（n+6 个 qubit，n = L-1；
#   对应旧 x1/MS/x2/ki/kx/b 布局）：
#       x1 = 1, MS = 2, x2 = 3, ki = 4, kx = 5, b = 6..n+6
#
#   LNN 受控模乘法器 / 模乘 U 块（2n+4 个 qubit）：
#       x = 1..n（其中 x[n-1] = n+1，位 n-1 为加法器 MS 辅助位）
#       ki = n+2, kx = n+3, b = n+4..2n+4
#
#   LNN 模指电路（m+2n+3 个 qubit）：
#       x = 1..m（计数/指数寄存器），ki = m+n+1，kx = m+n+2，
#       b = m+n+3..m+2n+3
#
# `PLUS` / `MINUS` 是 Beauregard/Fowler 技巧中的 ±π/4 全相位旋转
# （e^{±iπ/4}·RX(±π/2)）：受控使用时全局相位变为相对相位，故必须以
# ConstGate 保存精确矩阵，而不能用 RX 近似。
# =============================================================================

const PLUS = ConstGate(:PLUS, ComplexF64[0.5+0.5im 0.5-0.5im; 0.5-0.5im 0.5+0.5im])
const MINUS = ConstGate(:MINUS, ComplexF64[0.5-0.5im 0.5+0.5im; 0.5+0.5im 0.5-0.5im])

"""
    lnn_fourier_adder_circuit(a, L; base=1) -> Circuit

LNN Fourier 域整数加法的相位阶梯（不含 QFT 括号、不含受控化）：
`PHASE(compute_phase(a, i))` 作用于 `base+i`（i = 0..L-1）。
对应旧 `LNNAdderCircuit`。
"""
function lnn_fourier_adder_circuit(a::Integer, L::Integer; base::Int = 1)
    (0 <= a < 1 << L) ||
        throw(ArgumentError("a = $a does not fit into a $L-bit register."))
    c = Circuit(base + L - 1)
    for i in 0:(L-1)
        push!(c, PHASE(compute_phase(a, i), base + i))
    end
    return c
end

"""
    lnn_controlled_adder_circuit(a, L; base=1, inverse=false) -> Circuit

LNN 受控 Fourier 域加法：`CPHASE–SWAP` 阶梯，作用在
`(base, base+1, …, base+L)` 上（受控化由阶梯与 `ki/kx` 协同完成）。
`inverse=true` 为旧 `LNNCAdderCircuit(a, inverse=true)` 的逐门对应
（门序反转、CPHASE 控制目标互换，相位不变）。
"""
function lnn_controlled_adder_circuit(a::Integer, L::Integer; base::Int = 1,
                                      inverse::Bool = false)
    (0 <= a < 1 << L) ||
        throw(ArgumentError("a = $a does not fit into a $L-bit register."))
    c = Circuit(base + L)
    if !inverse
        for i in 0:(L-1)
            push!(c, CPHASE(compute_phase(a, i), base + i, base + i + 1))
            push!(c, SWAP(base + i, base + i + 1))
        end
    else
        for i in (L-1):-1:0
            push!(c, CPHASE(compute_phase(a, i), base + i + 1, base + i))
            push!(c, SWAP(base + i + 1, base + i))
        end
    end
    return c
end

# 内部便利：受控 PLUS / MINUS（控制位前置，与旧 CONTROL_PLUS/MINUS 一致）
_ctrl_plus(q::Int, c::Int) = ctrl(PLUS(q), c)
_ctrl_minus(q::Int, c::Int) = ctrl(MINUS(q), c)

"""
    lnn_modular_adder_circuit(N, a; base=1) -> Circuit

LNN 双受控模加法器（ki、kx 协同的近邻实现）。
canonical 布局：x1=1, MS=2, x2=3, ki=4, kx=5, b=6..n+6。
对应旧 `LNNModAdderCircuit`。
"""
function lnn_modular_adder_circuit(N::Integer, a::Integer; base::Int = 1)
    _check_register_args(N, a)
    L = _register_bits(N)
    n = L - 1
    c = Circuit(n + 6)
    push!(c, _ctrl_plus(5, 4))                 # gate((3,4), CONTROL_PLUS)
    push!(c, CX(3, 4))                         # gate((2,3), CNOT)
    push!(c, _ctrl_minus(5, 4))                # gate((3,4), CONTROL_MINUS)
    push!(c, CX(3, 4))                         # gate((2,3), SWAP*CNOT)：先 CNOT 后 SWAP
    push!(c, SWAP(3, 4))
    push!(c, _ctrl_plus(5, 4))                 # gate((3,4), CONTROL_PLUS)
    append!(c, shift(lnn_controlled_adder_circuit(a, L), 4))
    append!(c, dagger(shift(lnn_fourier_adder_circuit(N, L), 4)))
    append!(c, dagger(shift(qft_circuit(L), 4)))
    push!(c, SWAP(2, 3))                       # gate((1,2), SWAP)
    push!(c, SWAP(3, 4))                       # gate((2,3), SWAP)
    push!(c, CX(5, 4))                         # gate((4,3), SWAP*CNOT)：先 CNOT 后 SWAP
    push!(c, SWAP(5, 4))
    append!(c, move_site(5, n + 5))            # move_site(4, n+4)
    append!(c, shift(qft_circuit(L), 3))       # shift(_QFT(L), 3)
    append!(c, shift(lnn_controlled_adder_circuit(N, L; inverse=true), 3))
    append!(c, dagger(shift(lnn_controlled_adder_circuit(a, L), 4)))
    append!(c, dagger(shift(qft_circuit(L), 5)))
    append!(c, move_site(5, n + 6))            # move_site(4, n+5)
    push!(c, X(5))                             # gate(4, X)
    push!(c, CX(5, 4))                         # gate((4,3), CNOT)
    push!(c, X(5))                             # gate(4, X)
    push!(c, SWAP(3, 4))                       # gate((2,3), SWAP)
    push!(c, SWAP(2, 3))                       # gate((1,2), SWAP)
    push!(c, SWAP(1, 2))                       # gate((0,1), SWAP)
    append!(c, shift(qft_circuit(L), 4))       # shift(_QFT(L), 4)
    append!(c, shift(lnn_controlled_adder_circuit(a, L; inverse=true), 4))
    push!(c, _ctrl_minus(5, 4))                # gate((3,4), CONTROL_MINUS)
    push!(c, SWAP(3, 4))                       # gate((2,3), CNOT*SWAP)：先 SWAP 后 CNOT
    push!(c, CX(3, 4))
    push!(c, _ctrl_plus(5, 4))                 # gate((3,4), CONTROL_PLUS)
    push!(c, CX(3, 4))                         # gate((2,3), CNOT)
    push!(c, SWAP(2, 3))                       # gate((1,2), SWAP)
    push!(c, SWAP(1, 2))                       # gate((0,1), SWAP)
    push!(c, _ctrl_minus(5, 4))                # gate((3,4), CONTROL_MINUS)
    return base == 1 ? c : shift(c, base - 1)
end

"""
    lnn_controlled_modular_multiplier_circuit(N, a; base=1) -> Circuit

LNN 受控模乘-累加：`|ki⟩|x⟩|b⟩ → |ki⟩|x⟩|b + ki·a·x (mod N)⟩`。
canonical 布局：x=1..n（x[n-1] = n+1，位 n-1 为加法器 MS 辅助位），
ki = n+2, kx = n+3, b = n+4..2n+4。对应旧 `LNNCModMultiplierCircuit`。
"""
function lnn_controlled_modular_multiplier_circuit(N::Integer, a::Integer;
                                                   base::Int = 1)
    _check_register_args(N, a)
    L = _register_bits(N)
    n = L - 1
    c = Circuit(2n + 4)
    mv = move_site(n - 1, 1)                   # move_site(n-2, 0)
    for i in 0:(n-1)
        apow = mod(a * 2^i, N)
        iter = Circuit(2n + 4)
        append!(iter, shift(lnn_modular_adder_circuit(N, apow), n - 2))
        append!(iter, mv)
        push!(c, block(iter; name = Symbol("LNNADD(", apow, ")")))
    end
    return base == 1 ? c : shift(c, base - 1)
end

# 受控 SWAP（FREDKIN 的 LNN 分解），作用于 canonical (1,2,3)，
# 交换 1↔3、以 2 为控制。对应旧 `_cswap()`。
function _lnn_cswap()
    c = Circuit(3)
    push!(c, CX(3, 2))                         # CNOTGate(2, 1)
    push!(c, _ctrl_plus(3, 2))                 # gate((1,2), CONTROL_PLUS)
    push!(c, CX(1, 2))                         # CNOTGate(0, 1)
    push!(c, _ctrl_minus(3, 2))                # gate((1,2), CONTROL_MINUS)
    push!(c, CX(1, 2))                         # gate((0,1), SWAP*CNOT)：先 CNOT 后 SWAP
    push!(c, SWAP(1, 2))
    push!(c, _ctrl_plus(3, 2))                 # gate((1,2), SWAP*CONTROL_PLUS)：先 CPLUS 后 SWAP
    push!(c, SWAP(2, 3))
    push!(c, CX(2, 1))                         # CNOTGate(1, 0)
    return c
end

# LNN 交织置换。对应旧 `mesh(n)`。
_lnn_mesh(n::Int) = begin
    c = Circuit(2n)
    for i in 0:(n-2)
        append!(c, move_site(n + i + 1, 2i + 2))   # move_site(n+i, 2i+1)
    end
    c
end

# 双 n 位寄存器的 LNN 受控交换，作用于 canonical (2, …, 2n+1)。
# 对应旧 `control_swap(n)`。
function _lnn_control_swap(n::Int)
    c = Circuit(2n + 1)
    append!(c, shift(_lnn_mesh(n), 1))
    for i in 0:2:(2n-1)
        append!(c, shift(_lnn_cswap(), i))
    end
    append!(c, dagger(_lnn_mesh(n)))
    return c
end

"""
    lnn_modular_u_circuit(N, a, ra; base=1) -> Circuit

LNN 受控模乘酉算子 U_a：`|ki⟩|x⟩|b=0⟩ → |ki⟩|a^ki·x (mod N)⟩|b=0⟩`，
`ra` 必须是 `a` 模 `N` 的逆（`invmod(a, N)`）。
canonical 布局同 `lnn_controlled_modular_multiplier_circuit`。
对应旧 `LNNModUCircuit`。
"""
function lnn_modular_u_circuit(N::Integer, a::Integer, ra::Integer; base::Int = 1)
    _check_register_args(N, a)
    (mod(a * ra, N) == 1) ||
        throw(ArgumentError("ra = $ra is not the inverse of a = $a modulo N = $N."))
    L = _register_bits(N)
    n = L - 1
    c = Circuit(2n + 4)
    push!(c, block(lnn_controlled_modular_multiplier_circuit(N, a);
                     name = Symbol("LNN-CMUL(", a, ")")))
    mid = Circuit(2n + 4)
    append!(mid, sdagger(shift(qft_circuit(L), n + 3)))  # shift(_QFT(L))', n+3
    append!(mid, move_site(n, 1))                        # move_site(n-1, 0)
    append!(mid, move_site(n + 3, 2))                    # move_site(n+2, 1)
    append!(mid, move_site(n + 3, 3))                    # move_site(n+2, 2)
    append!(mid, shift(move_site(1, n + 1), n + 3))      # shift(move_site(0, n), n+3)
    append!(mid, shift(_lnn_control_swap(n), 2))         # shift(control_swap(n), 2)
    append!(mid, move_site(2, n + 2))                    # move_site(1, n+1)
    append!(mid, move_site(1, n))                        # move_site(0, n-1)
    append!(mid, move_site(2n + 3, n + 2))               # move_site(2n+2, n+1)
    append!(mid, move_site(2n + 4, n + 4))               # move_site(2n+3, n+3)
    append!(mid, shift(qft_circuit(L), n + 3))           # shift(_QFT(L), n+3)
    push!(c, block(mid; name = :LNNUNCOMP))
    push!(c, block(sdagger(lnn_controlled_modular_multiplier_circuit(N, ra));
                     name = Symbol("LNN-CMUL(", ra, ")†")))
    return base == 1 ? c : shift(c, base - 1)
end

"""
    lnn_modular_exponentiator_circuit(N, a, m; base=1) -> Circuit

LNN 受控模指：`|x⟩|ki=1⟩ → |x⟩|aˣ (mod N)⟩`（调用前须先把 ki 置 |1⟩——
首次 SWAP 链会把它注入 kx，见 `order_finding_circuit(; alg=LNNQFT())` 的用法）。
布局（m+2n+3 个 qubit，L = bits(N)+1，n = L-1）：x = 1..m（计数/指数
寄存器），ki = m+n+1，kx = m+n+2，b = m+n+3..m+2n+3。
对应旧 `LNNModExponentiatorCircuit`。
"""
function lnn_modular_exponentiator_circuit(N::Integer, a::Integer, m::Int;
                                           base::Int = 1)
    _check_register_args(N, a)
    (m >= 1) || throw(ArgumentError("m must be positive"))
    L = _register_bits(N)
    n = L - 1
    ra = invmod(a, N)
    qft_blk = shift(qft_circuit(L), m + n + 2)           # shift(_QFT(L), m+n+2)
    c = Circuit(m + 2n + 3)
    push!(c, block(qft_blk; name = :QFT_b))
    for i in (m-1):-1:0
        mv = move_site(i + 1, m + n + 1)                 # move_site(i, m+n)
        rnd = Circuit(m + 2n + 3)
        append!(rnd, mv)
        append!(rnd, shift(lnn_modular_u_circuit(N, a, ra), m - 1))
        append!(rnd, sdagger(mv))
        push!(c, block(rnd; name = Symbol("U(", a, ")")))
        a = mod(a * a, N)
        ra = mod(ra * ra, N)
    end
    push!(c, block(sdagger(qft_blk); name = Symbol("QFT_b†")))
    return base == 1 ? c : shift(c, base - 1)
end
