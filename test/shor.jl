using Test
using Random
using ShorAlgorithm
using QuantumCircuits
using VQC

# ── 测试辅助 ──────────────────────────────────────────────────────────────────

"""在寄存器 `qs`（`qs[1]` = MSB）上装载整数值 `value`（基态制备）。"""
function load_value!(state, qs::Vector{Int}, value::Int)
    L = length(qs)
    for (j, q) in enumerate(qs)
        ((value >> (L - j)) & 1) == 1 && apply!(state, X(q))
    end
    return state
end

"""读出处于计算基态的寄存器（`qs[1]` = MSB）：返回 (值, 峰值概率)。"""
function read_value(state, qs::Vector{Int})
    ps = probabilities(state)
    v = argmax(ps) - 1                     # 小端序基矢索引
    val = 0
    for q in qs                            # qs[1] = MSB
        val = (val << 1) | ((v >> (q - 1)) & 1)
    end
    return val, ps[v+1]
end

"""`L` 位小端序基矢索引 ↔ MSB-first 位权值（位反转）。"""
bitrev(j::Int, L::Int) = sum(((j >> (L - 1 - k)) & 1) << k for k in 0:(L-1); init=0)

@testset "ShorAlgorithm" verbose = true begin
    @testset "compute_phase" begin
        @test compute_phase(13, 0) ≈ π                 # 2π·1/2
        @test compute_phase(13, 1) ≈ π / 2             # 2π·1/4
        @test compute_phase(13, 3) ≈ 13π / 8           # 2π·13/16
        @test compute_phase(0, 5) == 0.0
        @test compute_phase(1, 5) ≈ 2π / 64
        @test compute_phase(63, 5) ≈ 2π * 63 / 64
    end

    @testset "QFT 电路 vs DFT" begin
        # 小端序基矢索引 j 处的振幅 = e^(2πi·x·y/d)/√d，y = bitrev(j)
        for L in 1:4
            d = 1 << L
            circ = qft_circuit(L)
            for x in 0:d-1
                state = load_value!(zero_state(L), collect(1:L), x)
                simulate!(circ, state)
                for j in 0:d-1
                    y = bitrev(j, L)
                    @test storage(state)[j+1] ≈ cis(2π * x * y / d) / sqrt(d) atol = 1e-12
                end
            end
        end
        # 逆 QFT 可逆性
        state = load_value!(zero_state(4), collect(1:4), 11)
        simulate!(qft_circuit(4), state)
        simulate!(dagger(qft_circuit(4)), state)
        v, p = read_value(state, collect(1:4))
        @test v == 11 && p > 1 - 1e-10
    end

    @testset "Fourier 加法器" begin
        L = 3
        qs = collect(1:L)
        for a in (1, 5)
            for b in 0:(1 << L)-1
                state = load_value!(zero_state(L), qs, b)
                simulate!(qft_circuit(qs), state)
                simulate!(fourier_adder_circuit(a, qs), state)
                simulate!(dagger(qft_circuit(qs)), state)
                v, p = read_value(state, qs)
                @test v == (b + a) % (1 << L)
                @test p > 1 - 1e-10
            end
        end
    end

    @testset "模加法器（Beauregard）" begin
        # 寄存器位宽 L 须满足 2N ≤ 2^L（Beauregard 条件），取 L = bits(7)+1 = 4
        N = 7
        L = 4
        circ = modular_adder_circuit(N, 5)
        # 布局：t=1, b=2..5（2 为 MSB），c2=6, c1=7
        # 契约：b 寄存器 Fourier 域进、Fourier 域出 → 测试时加 QFT 括号
        qft = qft_circuit(collect(2:(L+1)))
        for b in 0:N-1
            for (c1, c2) in ((1, 1), (1, 0), (0, 1), (0, 0))
                state = zero_state(7)
                load_value!(state, [2, 3, 4, 5], b)
                c1 == 1 && apply!(state, X(7))
                c2 == 1 && apply!(state, X(6))
                simulate!(qft, state)
                simulate!(circ, state)
                simulate!(dagger(qft), state)
                v, p = read_value(state, [2, 3, 4, 5])
                tval, _ = read_value(state, [1])
                expected = (c1 == 1 && c2 == 1) ? (b + 5) % N : b
                @test v == expected
                @test tval == 0                            # t 位复位
                @test p > 1 - 1e-10
            end
        end
    end

    @testset "受控模乘法器（Beauregard）" begin
        N = 7
        L = 4
        n = L - 1
        circ = controlled_modular_multiplier_circuit(N, 3)
        # 布局：t=1, b=2..5, win=6..8（6 为 MSB），c=9
        # 契约：|c⟩|win⟩|b⟩ → |c⟩|win⟩|b + c·a·win (mod N)⟩（Fourier 域进、出）
        qft = qft_circuit(collect(2:(L+1)))
        for c in (0, 1)
            for win in 0:(1 << n)-1
                state = zero_state(9)
                load_value!(state, [6, 7, 8], win)
                c == 1 && apply!(state, X(9))
                simulate!(qft, state)
                simulate!(circ, state)
                simulate!(dagger(qft), state)
                vwin, p = read_value(state, [6, 7, 8])
                vb, _ = read_value(state, [2, 3, 4, 5])
                vt, _ = read_value(state, [1])
                @test vwin == win                            # win 寄存器不变
                @test vb == (c == 1 ? (3 * win) % N : 0)     # b 累加 c·a·win
                @test vt == 0
                @test p > 1 - 1e-10
            end
        end
    end

    @testset "模指电路（Beauregard）" begin
        # N=15 需 5 位寄存器（2N ≤ 2^L）
        N = 15
        L = 5
        n = L - 1
        m = 2
        circ = modular_exponentiator_circuit(N, 7, m)
        # 布局：t=1, b=2..6, win=7..10（7 为 MSB），x=11..12（11 为 MSB）
        for e in 0:(1 << m)-1
            state = zero_state(2n + 2 + m)
            apply!(state, X(2n + 2))                       # win LSB = 1
            load_value!(state, collect((2n+3):(2n+2+m)), e)
            simulate!(circ, state)
            vwin, p = read_value(state, collect((L+2):(2n+2)))
            vb, _ = read_value(state, collect(2:(L+1)))
            @test vwin == powermod(7, e, N)
            @test vb == 0
            @test p > 1 - 1e-10
        end
    end

    @testset "连分数 / 素性 / 完全幂" begin
        @test continued_fraction(0.5, 10) == 1 // 2
        @test continued_fraction(1 / 3, 10) == 1 // 3
        @test continued_fraction(0.75, 7) == 3 // 4
        @test continued_fraction(0.0, 5) == 0 // 1
        # 相位估计读数 y = round(s·2^m/r)/2^m 的还原
        s, r, m = 3, 11, 8
        y = round(Int, s * (1 << m) / r) / (1 << m)
        @test continued_fraction(y, 20) == 3 // 11
        @test is_prime(97) && !is_prime(91) && is_prime(7919) && !is_prime(1)
        @test perfect_power(27) == (3, 3)
        @test perfect_power(64) == (8, 2)
        @test perfect_power(15) === nothing
    end

    @testset "求阶：完整 QFT（N=15）" begin
        # 完整线路为 18 量子比特，代价较高，仅抽查有代表性的阶
        orders = Dict(2 => 4, 4 => 2, 7 => 4, 14 => 2)
        for (a, r) in orders
            y = order_finding(a, 15; alg=FullQFT())
            k = round(Int, y * r)
            @test abs(y - k / r) < 0.01
        end
    end

    @testset "求阶：单控制比特半经典" begin
        # 修复 modular_u_circuit 缺失 sdagger 后（见 beauregard.jl），
        # 每次测量都应落在 s/r 格点上。
        for (a, r) in ((2, 6), (4, 3), (8, 2), (10, 6), (17, 6))
            for _ in 1:3
                y = order_finding(a, 21)              # 默认 alg = SemiclassicalQFT()
                k = round(Int, y * r)
                @test abs(y - k / r) < 0.05
            end
        end
        @test_throws ArgumentError order_finding(1, 15)
        @test_throws ArgumentError order_finding(3, 15)   # gcd(3,15) ≠ 1
        @test_throws TypeError order_finding(2, 15; alg="bad")
    end

    @testset "求阶：单控制比特半经典 LNN" begin
        # Fowler 布局上的 Griffiths–Niu 反馈（SemiclassicalLNN）：
        # 每次测量都应落在 s/r 格点上。
        for (a, N, r) in ((7, 15, 4), (2, 21, 6), (4, 21, 3), (8, 21, 2), (10, 21, 6))
            for _ in 1:3
                y = order_finding(a, N; alg=SemiclassicalLNN())
                k = round(Int, y * r)
                @test abs(y - k / r) < 0.05
            end
        end
    end

    @testset "打包线路：order_finding_circuit" begin
        # FullQFT / LNNQFT：纯酉线路（位宽冒烟；求解分布已由求阶 testset 覆盖）
        @test QuantumCircuits.nqubits(order_finding_circuit(7, 15; alg=FullQFT())) == 18
        @test QuantumCircuits.nqubits(order_finding_circuit(7, 15; alg=LNNQFT())) == 19
        # SemiclassicalQFT：测量 + 经典反馈（IfOp）单线路，out 寄存器给出解
        c = order_finding_circuit(7, 15; alg=SemiclassicalQFT())
        @test QuantumCircuits.nqubits(c) == 19
        # 块级视图：重复子线路以命名块为单元（QFT_b / U(a·)）
        s_view = QuantumCircuits.draw(c; unroll=false)
        @test occursin("QFT_b", s_view) && occursin("U(7)", s_view)
        for _ in 1:5
            out = simulate(c, zero_state(19))       # 驱动测量 / IfOp / 复位
            pm = marginal_probabilities(out, collect(12:19)) # out 寄存器边缘分布（12 为 LSB）
            v = argmax(pm) - 1
            @test v % 64 == 0 && v ÷ 64 in (0, 1, 2, 3)   # y = v/2^8 ≈ s/4
            @test maximum(pm) > 1 - 1e-10             # out 寄存器确定（辅助位不承诺干净）
        end
        # SemiclassicalLNN：LNN 布局的半经典打包线路（win/ki/kx/b + out）
        cl = order_finding_circuit(7, 15; alg=SemiclassicalLNN())
        @test QuantumCircuits.nqubits(cl) == 20      # 2n+4+m = 4n+4
        sl = QuantumCircuits.draw(cl; unroll=false)
        @test occursin("U(7)", sl)
        for _ in 1:3
            out = simulate(cl, zero_state(20))
            pm = marginal_probabilities(out, collect(13:20))  # out 寄存器边缘分布（13 为 LSB）
            v = argmax(pm) - 1
            @test v % 64 == 0 && v ÷ 64 in (0, 1, 2, 3)   # y = v/2^8 ≈ s/4
            @test maximum(pm) > 1 - 1e-10
        end
    end

    @testset "shor 分解" begin
        rng = MersenneTwister(42)
        # small_primes 置空：绕过试除捷径，真正走量子求阶路径
        kw = (small_primes=Int[], rng=rng)
        p, q = shor(15; kw...)
        @test p * q == 15 && 1 < p < 15
        p, q = shor(21; kw...)
        @test p * q == 21 && Set((p, q)) == Set((3, 7))
        p, q = shor(33; kw...)
        @test p * q == 33 && Set((p, q)) == Set((3, 11))
        p, q = shor(35; kw...)
        @test p * q == 35 && Set((p, q)) == Set((5, 7))
        # 完整 QPE 线路驱动的分解（计算量更大，仅验一个小例）
        p, q = shor(15; kw..., alg=FullQFT())
        @test p * q == 15 && Set((p, q)) == Set((3, 5))
        # 试除捷径
        p, q = shor(3 * 5 * 47 * 97)
        @test p * q == 3 * 5 * 47 * 97
        @test_throws ArgumentError shor(13)                   # 素数
    end

    @testset "求阶：LNN（Fowler）" begin
        # 注：LNN 电路为旧实现的逐门移植，端到端以相位读数落格校验
        y = order_finding(7, 15; alg=LNNQFT())
        k = round(Int, y * 4)
        @test abs(y - k / 4) < 0.05
    end
end
