# =============================================================================
# beauregard.jl — Beauregard 型 QFT 模算术电路（全连接）
#
# 移植自 VQC `src/utility/shor/Beauregard.jl`（0-based → 1-based）。
# 寄存器位宽统一取 `L = bits(N) + 1`（`_register_bits`），保证 2N ≤ 2^L。
#
# 寄存器布局（1-based，寄存器内**首 qubit = MSB**）：
#
#   Fourier 加法器（作用于 L 位寄存器 qs）：
#       qs[1..L]，qs[1] = MSB（权重 2^(L-1)）
#
#   模加法器 `modular_adder_circuit`（L+3 个 qubit）：
#       t  = 1           溢出辅助位，|0⟩ 进 |0⟩ 出
#       b  = 2..L+1      被加数寄存器，b[0]（MSB）= 2
#       c2 = L+2         控制位 2
#       c1 = L+3         控制位 1
#
#   受控模乘法器 / 模乘 U 块（2n+3 个 qubit，n = L-1）：
#       t   = 1
#       b   = 2..L+1
#       win = L+2..2n+2  乘数窗口寄存器（n 位）
#       c   = 2n+3       总控制位
#
#   模指电路（2n+2+m 个 qubit）：
#       t   = 1
#       b   = 2..L+1     恒为 |0…0⟩
#       win = L+2..2n+2  累计寄存器，初值 |1⟩（LSB = 2n+2），输出 a^x mod N
#       x   = 2n+3..2n+2+m  指数寄存器
# =============================================================================

"""
    fourier_adder_circuit(a, qs) -> Circuit

Fourier 域整数加法（不含 QFT 括号）：寄存器 `qs`（`qs[1]` = MSB，宽度
`L = length(qs)`）须已处于 Fourier 基；对第 `j` 位施加
`PHASE(compute_phase(a, j))`。要求 `0 ≤ a < 2^L`。
对应旧 `QFTAdderCircuit`。
"""
function fourier_adder_circuit(a::Integer, qs::AbstractVector{Int})
    (0 <= a < 1 << length(qs)) ||
        throw(ArgumentError("a = $a does not fit into a $(length(qs))-bit register."))
    c = Circuit(maximum(qs))
    for j in 0:(length(qs)-1)
        push!(c, PHASE(compute_phase(a, j), qs[j+1]))
    end
    return c
end

"""
    controlled_fourier_adder_circuit(a, qs, ctrl) -> Circuit

控制位 `ctrl` 下的 Fourier 域加法（`qs` 已在 Fourier 基）。
对应旧 `QFTCAdderCircuit`。
"""
function controlled_fourier_adder_circuit(a::Integer, qs::AbstractVector{Int},
                                          ctl::Int)
    (0 <= a < 1 << length(qs)) ||
        throw(ArgumentError("a = $a does not fit into a $(length(qs))-bit register."))
    c = Circuit(maximum(qs))
    for j in 0:(length(qs)-1)
        push!(c, ctrl(PHASE(compute_phase(a, j), qs[j+1]), ctl))
    end
    return c
end

"""
    double_controlled_fourier_adder_circuit(a, qs, c1, c2) -> Circuit

双控制 Fourier 域加法（`c1`、`c2` 同时为 1 才加）。对应旧 `QFTCCAdderCircuit`。
"""
function double_controlled_fourier_adder_circuit(a::Integer, qs::AbstractVector{Int},
                                                 c1::Int, c2::Int)
    (0 <= a < 1 << length(qs)) ||
        throw(ArgumentError("a = $a does not fit into a $(length(qs))-bit register."))
    c = Circuit(maximum(qs))
    for j in 0:(length(qs)-1)
        push!(c, ctrl(PHASE(compute_phase(a, j), qs[j+1]), c1, c2))
    end
    return c
end

"""
    modular_adder_circuit(N, a) -> Circuit

双受控模加法器：`|c1 c2⟩|b⟩ → |c1 c2⟩|b + c1c2·a (mod N)⟩`，
辅助位 `t`（qubit 1）保持 |0⟩；b 寄存器以 Fourier 基进、Fourier 基出。
要求 `0 ≤ a < N`。对应旧 `QFTCCModAdderCircuit`（含溢出回加与 t 位复位）。
"""
function modular_adder_circuit(N::Integer, a::Integer)
    _check_register_args(N, a)
    L = _register_bits(N)
    qs_b = collect(2:(L+1))                    # b 寄存器，qs_b[1] = MSB
    qft = qft_circuit(qs_b)
    c1, c2 = L + 3, L + 2
    c = Circuit(L + 3)
    # 1. CC-add(a)
    append!(c, double_controlled_fourier_adder_circuit(a, qs_b, c1, c2))
    # 2. add(-N)
    append!(c, dagger(fourier_adder_circuit(N, qs_b)))
    # 3-4. 回到计算基，把溢出位（b 的 MSB）异或到 t
    append!(c, dagger(qft))
    push!(c, CX(2, 1))
    # 5-6. 若 t=1（未溢出）则加回 N
    append!(c, qft)
    append!(c, controlled_fourier_adder_circuit(N, qs_b, 1))
    # 7-8. 减回 a
    append!(c, dagger(double_controlled_fourier_adder_circuit(a, qs_b, c1, c2)))
    append!(c, dagger(qft))
    # 9-11. 复位 t
    push!(c, X(2))
    push!(c, CX(2, 1))
    push!(c, X(2))
    # 12-13. 回 Fourier 基，再加 a
    append!(c, qft)
    append!(c, double_controlled_fourier_adder_circuit(a, qs_b, c1, c2))
    return c
end

"""
    controlled_modular_multiplier_circuit(N, a) -> Circuit

受控模乘-累加：`|c⟩|x⟩|b⟩ → |c⟩|x⟩|b + c·a·x (mod N)⟩`。
b 以 Fourier 基进、Fourier 基出。对应旧 `QFTCModMultiplierCircuit`。
"""
function controlled_modular_multiplier_circuit(N::Integer, a::Integer)
    _check_register_args(N, a)
    L = _register_bits(N)
    n = L - 1
    c = Circuit(2n + 3)
    mv1 = move_site(2n + 3, n + 4)             # 把 c 搬到模加法器 c1 位
    for i in 0:(n-1)
        mv2 = move_site(2n + 2 - i, n + 3)     # 把 win 的第 n-1-i 位搬到 c2 位
        apow = mod(a * 2^i, N)
        iter = Circuit(2n + 3)
        append!(iter, mv2)
        append!(iter, mv1)
        append!(iter, modular_adder_circuit(N, apow))
        append!(iter, sdagger(mv1))
        append!(iter, sdagger(mv2))
        push!(c, block(iter; name = Symbol("CMADD(", apow, ")")))
    end
    return c
end

"""
    modular_u_circuit(N, a, ra) -> Circuit

受控模乘酉算子 U_a：`|c⟩|win⟩|b=0⟩ → |c⟩|a^c·win (mod N)⟩|b=0⟩`，
其中 `ra` 必须是 `a` 模 `N` 的逆（`invmod(a, N)`）。
结构：受控乘 → 逆 QFT → 受控交换 b/win → QFT → 逆受控乘（b 位清零技巧）。
对应旧 `QFTModUCircuit`。
"""
function modular_u_circuit(N::Integer, a::Integer, ra::Integer)
    _check_register_args(N, a)
    (mod(a * ra, N) == 1) ||
        throw(ArgumentError("ra = $ra is not the inverse of a = $a modulo N = $N."))
    L = _register_bits(N)
    n = L - 1
    qs_b = collect(2:(L+1))
    qft = qft_circuit(qs_b)
    c = Circuit(2n + 3)
    push!(c, block(controlled_modular_multiplier_circuit(N, a);
                     name = Symbol("CMUL(", a, ")")))
    mid = Circuit(2n + 3)
    append!(mid, sdagger(qft))
    append!(mid, controlled_swap_block(2n + 3, n, 3, n + 3))   # b[1..n] ↔ win
    append!(mid, qft)
    push!(c, block(mid; name = Symbol("CSWAP+", L)))
    push!(c, block(sdagger(controlled_modular_multiplier_circuit(N, ra));
                     name = Symbol("CMUL(", ra, ")†")))
    return c
end

"""
    modular_exponentiator_circuit(N, a, m) -> Circuit

受控模指：`|x⟩|win=1⟩ → |x⟩|aˣ (mod N)⟩`，b 寄存器（qubit 2..L+1）恒为 0。

布局（2n+2+m 个 qubit，L = bits(N)+1，n = L-1）：t=1，b=2..L+1，
win=L+2..2n+2（初值 |1⟩，LSB = qubit `2n+2`），x=2n+3..2n+2+m（m 位指数，
高位在前处理：第 k 块施加 U_{a^(2^k)}）。对应旧 `QFTCModExponentiatorCircuit`。
"""
function modular_exponentiator_circuit(N::Integer, a::Integer, m::Int)
    _check_register_args(N, a)
    (m >= 1) || throw(ArgumentError("m must be positive"))
    L = _register_bits(N)
    n = L - 1
    ra = invmod(a, N)
    qs_b = collect(2:(L+1))
    qft = qft_circuit(qs_b)
    c = Circuit(2n + 2 + m)
    push!(c, block(qft; name = :QFT_b))      # b 进入 Fourier 基
    for i in (m-1):-1:0
        mv = shift(move_site(i + 1, 1), 2n + 2)   # 把 x[i] 搬到乘法器控制位 c
        rnd = Circuit(2n + 2 + m)
        append!(rnd, mv)
        append!(rnd, controlled_modular_multiplier_circuit(N, a))
        append!(rnd, sdagger(qft))
        append!(rnd, controlled_swap_block(2n + 3, n, 3, n + 3))  # b[1..n] ↔ win
        append!(rnd, qft)
        append!(rnd, sdagger(controlled_modular_multiplier_circuit(N, ra)))
        append!(rnd, sdagger(mv))
        push!(c, block(rnd; name = Symbol("U(", a, ")")))   # 受控 U_{a^(2^i)}
        a = mod(a * a, N)
        ra = mod(ra * ra, N)
    end
    push!(c, block(sdagger(qft); name = Symbol("QFT_b†")))  # b 回计算基
    return c
end
