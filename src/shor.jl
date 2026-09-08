# =============================================================================
# shor.jl — 完整 Shor 因数分解（经典预处理 + 量子求阶 + 经典后处理）
#
# 量子部分只承担"求阶"一步；其余（随机选基、gcd 捷径、连分数展开、
# 奇偶性/平方根校验）为经典计算。对应旧实现中注释掉的
# `continuous_fraction` 后处理，按教科书算法补全。
# =============================================================================

"""
    continued_fraction(y, qmax) -> Rational{Int}

返回分母不超过 `qmax` 且最接近 `y ∈ [0,1)` 的分数（连分数渐近分数）。
用于把相位估计读数 `y ≈ s/r` 还原为 `s//r`。
"""
function continued_fraction(y::Real, qmax::Int)
    (0 <= y <= 1) || throw(ArgumentError("y must lie in [0, 1], got $y"))
    qmax >= 1 || throw(ArgumentError("qmax must be positive"))
    x = Float64(y)
    # y == 1 与 y == 0 的平凡情形
    x ≥ 1 && return 1 // 1
    num_prev, den_prev = 0, 1     # 渐近分数 k-2
    num, den = 1, 0               # 渐近分数 k-1
    best = 0 // 1
    best_err = abs(y - 0.0)
    for _ in 1:64
        ai = floor(Int, x)
        num_next = ai * num + num_prev
        den_next = ai * den + den_prev
        (den_next > qmax) && break
        num_prev, num = num, num_next
        den_prev, den = den, den_next
        cand = Rational{Int}(num, den)
        err = abs(Float64(cand) - y)
        if err < best_err
            best = cand
            best_err = err
        end
        frac = x - ai
        (frac <= 0.0) && break          # 精确终止
        x = 1 / frac
    end
    return best
end

"""
    is_prime(n) -> Bool

Miller–Rabin 素性测试；对 `n < 3.3×10²⁴` 使用固定基组为确定性判定。
"""
function is_prime(n::Integer)
    n = Int(n)
    (n < 2) && return false
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
        (n % p == 0) && return n == p
    end
    d = n - 1
    s = 0
    while d % 2 == 0
        d = d >> 1
        s += 1
    end
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
        x = powermod(a, d, n)
        (x == 1 || x == n - 1) && continue
        composite = true
        for _ in 1:(s-1)
            x = (x * x) % n
            (x == n - 1) && (composite = false; break)
        end
        composite && return false
    end
    return true
end

"""
    perfect_power(n) -> Union{Tuple{Int,Int},Nothing}

若 `n = p^k`（k ≥ 2）则返回 `(p, k)`，否则返回 `nothing`。
"""
function perfect_power(n::Integer)
    n = Int(n)
    (n >= 4) || return nothing
    for k in 2:ndigits(n, base=2)
        r = round(Int, n^(1.0 / k))
        for c in (r-1, r, r+1)
            (c >= 2 && c^k == n) && return (c, k)
        end
    end
    return nothing
end

"""
    order_from_phase(y, a, N) -> Union{Int,Nothing}

从相位读数 `y ≈ s/r` 中提取阶 `r`（连分数 + `a^r ≡ 1 (mod N)` 校验）。
无效时返回 `nothing`。
"""
function order_from_phase(y::Real, a::Int, N::Int)
    φ = continued_fraction(y, N)
    denominator(φ) == 0 && return nothing
    r = denominator(φ)
    (r >= 1 && powermod(a, r, N) == 1) ? r : nothing
end

"""
    shor(n; max_attempts=100, shots=3, alg=SemiclassicalQFT(),
         small_primes=[3, 5, …, 47], verbosity=0, rng=Random.default_rng())
         -> Tuple{Int,Int}

完整 Shor 因数分解：返回 `n` 的一对非平凡因子 `(p, q)`。

流程：偶数 / 小素数 / 完全幂的 classical 捷径 → 随机选 `a`，
`gcd(a, n) > 1` 直接得因子 → 量子求阶（`order_finding`）→ 连分数提取
`r` → `r` 偶且 `a^(r/2) ≢ -1 (mod n)` 时由
`gcd(a^(r/2) ∓ 1, n)` 得因子。

参数：

* `max_attempts`：随机基的最大尝试次数；
* `shots`：每次尝试的量子测量次数；
* `alg`：量子求阶算法（`QuantumOrderFindingAlgorithm` 子类型单例，
  见 `order_finding`）；
* `small_primes`：试除捷径的素数表（置空可强制走量子路径，测试用）；
* `rng`：随机源。

`n` 为素数时抛 `ArgumentError`。
"""
function shor(n::Integer;
              max_attempts::Int = 100,
              shots::Int = 3,
              alg::QuantumOrderFindingAlgorithm = SemiclassicalQFT(),
              small_primes::Vector{Int} = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
              verbosity::Int = 0,
              rng::AbstractRNG = Random.default_rng())
    n = Int(n)
    (n > 1) || throw(ArgumentError("n must be > 1."))
    is_prime(n) && throw(ArgumentError("$n is prime."))

    # ── classical 捷径 ──
    iseven(n) && return (2, n >> 1)
    for p in small_primes
        (n % p == 0) && return (p, n ÷ p)
    end
    pp = perfect_power(n)
    if pp !== nothing
        p, _ = pp
        (verbosity >= 1) && println("$n = $(pp[1])^$(pp[2]); splitting off $(pp[1]).")
        return (p, n ÷ p)          # k ≥ 2 保证非平凡
    end

    # ── 量子循环 ──
    for attempt in 1:max_attempts
        a = rand(rng, 2:(n-1))
        g = gcd(a, n)
        if g > 1
            (verbosity >= 1) && println("attempt $attempt: lucky gcd($a, $n) = $g.")
            return (g, n ÷ g)
        end
        for _ in 1:shots
            y = order_finding(a, n; alg=_with_verbosity(alg, verbosity))
            r = order_from_phase(y, a, n)
            (r === nothing || isodd(r)) && continue
            x = powermod(a, r >> 1, n)
            (x == n - 1) && continue        # a^(r/2) ≡ -1：退化，换基重试
            f1 = gcd(x - 1, n)
            f2 = gcd(x + 1, n)
            if 1 < f1 < n
                return (f1, n ÷ f1)
            elseif 1 < f2 < n
                return (f2, n ÷ f2)
            end
        end
        (verbosity >= 1) && println("attempt $attempt with a = $a failed; retrying.")
    end
    throw(ErrorException("shor: failed to factor $n within $max_attempts attempts."))
end
