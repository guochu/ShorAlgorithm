# ============================================================================
# mps_lnn_benchmark.jl
#
# 复现 Dang, Hill & Hollenberg (2018),
# "Optimising matrix product state simulations of Shor's algorithm"
# (arXiv:1712.07311) 中的计算（轻量版，单机可跑）：
#
#   Part 1  参数验证：论文 Table I / Table II / Fig. 2 全部 (l, N, a)
#           案例的周期 r、尾零数 α、奇因子 η = r/2^α（纯经典，秒级）。
#   Part 2  Fig. 1(c) 分布复现：l = 5, N = 21, a = 2（r = 6）——用
#           `SemiclassicalLNN` 打包线路（Fowler LNN 布局 + IPE 逐轮测量，
#           与论文的半经典 LNN 方案同族）在 MPSSimulator 的
#           MatrixProductBackend 上做 MPS 采样，统计 S 读数分布并验证
#           连分数恢复 r。
#
# 用法（在 MPSSimulator 仓库目录下运行，LOAD_PATH 方式解析本地包）：
#   julia --startup-file=no mps_lnn_benchmark.jl [shots]
#
# 论文的大规模算例（l = 16/17/20，分别用到 2/24/216 个超算节点）不在
# 本脚本默认范围内；l ≥ 11 的采样需要放宽 MPS 键维截断并预留更长运行
# 时间，可在工作站上以更大 shots 运行本脚本扩展。
# ============================================================================
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "QuantumCircuits"))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "MPSSimulator"))
push!(LOAD_PATH, joinpath(@__DIR__, "..", "..", "VQC"))
push!(LOAD_PATH, joinpath(@__DIR__, "..", ".."))

using QuantumCircuits
using MPSSimulator
using ShorAlgorithm

shots = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 50
seed = 20240703                            # 论文发表日，固定可复现

# ── 经典求阶（暴力幂模；论文案例 r ≤ 5×10^5，瞬间完成）────────────────────
function find_order(a::Int, N::Int)
    x = a % N
    r = 1
    while x != 1
        x = x * a % N
        r += 1
    end
    return r
end

# ────────────────────────────────────────────────────────────────────────────
# Part 1: 参数验证（对照论文 Table I / Table II / Fig. 2）
# ────────────────────────────────────────────────────────────────────────────
# (l, N, a, 论文 r, 论文 α, 论文 η)
paper_cases = [
    (11, 1943,   2,    924, 2,    231),   # Fig. 2
    (13, 8189,  10,   3870, 1,   1935),   # Table I
    (14, 16351,  2,   8036, 2,   2009),   # Table I
    (15, 32663,  6,  16104, 3,   2013),   # Table I
    (16, 56759,  2,  28140, 2,   7035),   # Table II
    (17, 124631, 2,  57516, 2,  14379),   # Table II
    (20, 961307, 5, 479568, 4,  29973),   # Table II (η = 29973)
]

println("="^78)
println("Part 1: 参数验证（l, N, a) → (r, α, η)，对照论文 Table I / II / Fig. 2")
println("="^78)
println("l    N       a    r       α   η      论文一致")

function verify_params(paper_cases)
    all_ok = true
    for (l, N, a, r_paper, α_paper, η_paper) in paper_cases
        r = find_order(a, N)
        α = trailing_zeros(r)
        η = r >> α
        ok = (r == r_paper && α == α_paper && η == η_paper)
        all_ok &= ok
        println(l, "   ", N, "   ", a, "   ", r, "   ", α, "   ", η, "   ", ok ? "✓" : "✗ (论文: $r_paper, $α_paper, $η_paper)")
    end
    return all_ok
end
all_ok = verify_params(paper_cases)
println(all_ok ? "全部一致 ✓" : "存在不一致 ✗")
println()

# ────────────────────────────────────────────────────────────────────────────
# Part 2: Fig. 1(c) 分布复现（l = 5, N = 21, a = 2 → r = 6, m = 2l = 10）
# ────────────────────────────────────────────────────────────────────────────
println("="^78)
println("Part 2: Fig. 1(c) 分布复现（l=5, N=21, a=2, r=6; SemiclassicalLNN + MPS 采样 $shots shots）")
println("="^78)

const N21 = 21
const a21 = 2
r21 = find_order(a21, N21)
n = ndigits(N21, base = 2)      # 5
m = 2n                          # 10 计数比特
println("qubits = 2n+4+m = ", 2n + 4 + m, "   r = ", r21, "   理想峰位 s ≈ k·2^m/r = ",
        join([round(Int, k * 2^m / r21) for k in 0:(r21-1)], ", "))

c = order_finding_circuit(a21, N21; alg = SemiclassicalLNN())
backend = MatrixProductBackend()
t0 = time()
res = QuantumCircuits.simulate(c, backend; shots = shots, seed = seed)
t1 = time()
println("采样耗时（含首次编译）: ", round(t1 - t0, digits = 1), " s")
println()

# ── 连分数后处理：遍历 y 的连分数收敛子分母 q，验证 a^q ≡ 1 (mod N) ────────
function recover_r(y::Float64, a::Int, N::Int, qmax::Int)
    frac = rationalize(y)                       # Float64 → 精确二进制分数
    num, den = numerator(frac), denominator(frac)
    q2, q1 = 1, 0                               # 收敛子分母递推 q_{k-2}=1, q_{k-1}=0
    while den != 0
        ak = num ÷ den
        num, den = den, num % den
        q = ak * q1 + q2
        q > qmax && break
        q2, q1 = q1, q
        q > 0 && powermod(a, q, N) == 1 && return q
    end
    return 0
end

# 解码 counts（键格式 "r:<比特串>"，比特串按测量顺序 r_1..r_m，LSB 位先测）
hist = Dict{Int,Int}()          # S 读数（s ≈ y·2^m）
valid = 0                       # 连分数恢复出合法阶的 shot 数
informative = 0                 # s ≠ 0 的 shot 数（s = 0 不含阶信息）
for (k, v) in res.counts
    dig = filter(ch -> ch in ('0', '1'), k)
    bits = [parse(Int, ch) for ch in dig]
    y = ShorAlgorithm.phase_to_float(reverse(bits))
    s = round(Int, y * (1 << m))
    hist[s] = get(hist, s, 0) + v
    # 连分数恢复 r：分母 ≤ 2^m 的最佳逼近给出 r 的因子；6 ∤ 2^m 时
    # 逐倍数验证（标准 Shor 经典后处理）
    s == 0 && continue
    global informative += v
    r_est = recover_r(y, a21, N21, 1 << m)
    if r_est > 0
        global valid += v
        println("  s=", lpad(s, 4), " ×", v, " → r = ", r_est, "  ✓ (a^r ≡ 1 mod ", N21, ")")
    else
        println("  s=", lpad(s, 4), " ×", v, " → r 恢复失败")
    end
end

println("S 读数分布（s: count）与峰位（理想峰间距 ≈ ", round(2^m / r21, digits = 1), "）:")
maxc = maximum(values(hist))
for s in sort(collect(keys(hist)))
    cnt = hist[s]
    bar = repeat("█", round(Int, 40 * cnt / maxc))
    println("  s=", lpad(s, 4), "  ", lpad(cnt, 3), "  ", bar)
end
println()
println("连分数恢复：", valid, "/", informative, " 个含信息 shots（s ≠ 0）得到合法阶 a^r ≡ 1 mod ", N21)
factor1 = gcd(powermod(a21, r21 ÷ 2, N21) + 1, N21)
println("由 r = ", r21, " 经典后处理得因子：", factor1, " × ", N21 ÷ factor1)
