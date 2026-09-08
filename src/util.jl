# =============================================================================
# util.jl — 电路工具：SWAP 链 / 平移 / 相位 / QFT / 受控交换块 / 校验
#
# 移植自 VQC `src/utility/shor/util.jl`（0-based 位置 → 1-based）。
# 模算术直接使用 Base 的 `mod` / `invmod` / `powermod`，不再需要
# 自定义定长整数类型；寄存器位宽由 `_register_bits(N)` 统一推导。
# 约定：本包各算术寄存器**首个（号码最小）qubit 为最高有效位**（MSB），
# `qft_circuit(qs)` 同样以 `qs[1]` 为 MSB，与旧 `_QFT` 的约定一致。
# =============================================================================

"""SWAP 链：把 qubit `i` 的内容搬到 qubit `j`（1-based；经相邻 SWAP 实现）。"""
function move_site(i::Int, j::Int)
    c = Circuit(max(i, j))
    (i == j) && return c
    if i < j
        for k in i:(j-1)
            push!(c, SWAP(k, k + 1))
        end
    else
        for k in i:-1:(j+1)
            push!(c, SWAP(k, k - 1))
        end
    end
    return c
end

"""把整条线路的所有量子比特索引平移 `k`（寄存器搬移/复用积木时使用）。

支持 `BlockOp`（递归平移块体；带 `mapping` 的块平移映射表）。"""
function shift(c::Circuit, k::Int)
    k == 0 && return copy(c)
    ops = Vector{Operation}(undef, length(c.ops))
    for (i, op) in enumerate(c.ops)
        ops[i] = if op isa GateOp
            GateOp(op.gate, op.qubits .+ k, op.params)
        elseif op isa ReinitOp
            ReinitOp(op.qubits .+ k)
        elseif op isa BarrierOp
            BarrierOp(op.qubits .+ k)
        elseif op isa BlockOp
            if op.mapping === nothing
                BlockOp(op.name, shift(op.body, k), op.n, nothing)
            else
                BlockOp(op.name, op.body, op.n, op.mapping .+ k)
            end
        else
            throw(ArgumentError("cannot shift operation of type $(typeof(op))"))
        end
    end
    return Circuit(ops)
end

"""
    compute_phase(a, j) -> Float64

Fourier 域加法的旋转角：`L` 位寄存器中权重为 `2^(L-1-j)`（`j` 从 0 计，
`qs[1]` = MSB）的比特施加 `PHASE(compute_phase(a, j))` 即在 Fourier 域
加上整数 `a`。闭式解 `2π·(a mod 2^(j+1)) / 2^(j+1)`，与旧实现的
逐位求和公式等价。
"""
compute_phase(a::Integer, j::Integer) = 2π * mod(a, 1 << (j + 1)) / (1 << (j + 1))

"""
    qft_circuit(qs) -> Circuit
    qft_circuit(L::Int)

量子 Fourier 变换（教科书线路，含末尾位序翻转 SWAP）。

`qs` 为寄存器的 1-based qubit 列表，**`qs[1]` 为最高有效位**；
`qft_circuit(L)` 等价于 `qft_circuit(collect(1:L))`。

所实现的酉算子：`|x⟩ → 2^(-L/2) Σ_y e^(2πi·x·y/2^L) |y⟩`，其中 `x`、`y`
按 `qs[1]` = MSB 的位权读值（与旧 VQC `_QFT` 逐矩阵等价）。
"""
function qft_circuit(qs::AbstractVector{Int})
    L = length(qs)
    L >= 1 || throw(ArgumentError("qubit list must be non-empty"))
    allunique(qs) || throw(ArgumentError("duplicate qubits in $qs"))
    c = Circuit(maximum(qs))
    for j in 1:L
        push!(c, H(qs[j]))
        for k in (j+1):L
            # 受控相位角 2π/2^(k-j+1)，控制 qs[k]、目标 qs[j]（CPHASE 对称）
            push!(c, CPHASE(π / 2.0^(k-j), qs[k], qs[j]))
        end
    end
    for j in 1:(L >> 1)
        push!(c, SWAP(qs[j], qs[L+1-j]))
    end
    return c
end

qft_circuit(L::Int) = qft_circuit(collect(1:L))

"""
    controlled_swap_block(c, n, i0, j0) -> Circuit

受控 SWAP 块：在控制比特 `c` 下交换两组寄存器
`(i0, i0+1, …, i0+n-1) ↔ (j0, j0+1, …, j0+n-1)`。
每组对用 `CNOT – CCX – CNOT` 分解（旧 `control_swap_block`）。
"""
function controlled_swap_block(c::Int, n::Int, i0::Int, j0::Int)
    (i0 == j0) && throw(ArgumentError("position not allowed."))
    blk = Circuit(max(c, i0 + n - 1, j0 + n - 1))
    for k in 0:(n-1)
        push!(blk, CX(j0+k, i0+k))
        push!(blk, CCX(c, i0+k, j0+k))
        push!(blk, CX(j0+k, i0+k))
    end
    return blk
end

# ── 寄存器位宽与参数校验 ─────────────────────────────────────────────────────

"""
模算术寄存器位宽：比 `N` 的比特数多一位（`L = bits(N) + 1`），
保证 `2N ≤ 2^L`——这是 Beauregard 模加法器溢出检测正确性的充要条件。
所有电路积木统一由 `N` 推导位宽，组合时天然一致。
"""
_register_bits(N::Integer) = ndigits(N, base=2) + 1

function _check_register_args(N::Integer, a::Integer)
    (N >= 2) || throw(ArgumentError("N must be >= 2, got $N."))
    (0 <= a < N) || throw(ArgumentError("a must satisfy 0 <= a < N, got a = $a, N = $N."))
    return nothing
end

# ── 测量比特串 → 相位浮点 ────────────────────────────────────────────────────

"""把测量比特串 `bits`（`bits[1]` 权重 1/2）读成 [0,1) 内的相位小数。"""
function phase_to_float(bits::AbstractVector{Int})
    y = 0.0
    for (i, b) in enumerate(bits)
        y += b * 2.0^(-i)
    end
    return y
end

"""
    sdagger(c) -> Circuit

块感知的 `dagger`：逆序 + 逐操作求逆；`BlockOp` 递归求逆并在块名后
附加 `†`，使 `draw` 中的正/逆块可区分。对不含块的线路与 `dagger`
等价。
"""
function sdagger(c::Circuit)
    ops = Vector{Operation}(undef, length(c.ops))
    for (i, op) in enumerate(reverse(c.ops))
        ops[i] = if op isa BlockOp
            BlockOp(Symbol(op.name, "†"), sdagger(op.body), op.n, op.mapping)
        else
            inv(op)
        end
    end
    return Circuit(ops; cregs = copy(c.cregs))
end

# ── 求阶输入校验 ─────────────────────────────────────────────────────────────

function _check_order_finding_input(a::Int, N::Int)
    (a == 1) && throw(ArgumentError("the input a can not be 1."))
    (N == 1) && throw(ArgumentError("the input N can not be 1."))
    (gcd(a, N) == 1) || throw(ArgumentError("a and N are not co-prime to each other."))
    return nothing
end
