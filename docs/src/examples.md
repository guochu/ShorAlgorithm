# 示例：N = 15, a = 7

本文以最小的非平凡例子 **N = 15 = 3 × 5**、**a = 7**（`gcd(7, 15) = 1`）
完整走一遍 Shor 算法的量子核——量子求阶，并画出其中使用的 QFT 线路。

```@example ex
using ShorAlgorithm
using QuantumCircuits
using VQC   # 态矢量模拟后端

N, a = 15, 7
n = ndigits(N, base=2)      # 4：寄存器位宽 L = n+1，计数位 m = 2n
```

## 1. 量子求阶

`a = 7` 模 `N = 15` 的阶为 `r = 4`（7, 4, 13, 1）。默认算法
`SemiclassicalQFT`（单控制比特 + 半经典 QFT，仅 2n+3 = 11 个 qubit，
n = bits(N) = 4）返回相位估计读数 `y ≈ s/r`：

```@example ex
y = order_finding(a, N)          # 2n+3 = 11 qubit
```

用连分数把 `y ≈ s/r` 还原为阶：

```@example ex
φ = continued_fraction(y, N)
```

```@example ex
r = denominator(φ)               # → 4
powermod(a, r, N)                # 校验：a^r mod N == 1
```

## 2. 求阶中使用的 QFT 线路

两种完整 QPE 算法（`FullQFT` / `LNNQFT`）以及半经典版的 b 寄存器初始
QFT、计数寄存器读出，都使用同一个 QFT 构造器
[`qft_circuit`](@ref)。它实现

`|x⟩ → 2^(-L/2) Σ_y e^(2πi·x·y/2^L) |y⟩`（含末尾位序翻转 SWAP）。

以 3 位为例画出其线路（为清晰起见去掉默认经典寄存器）：

```@example ex
nocreg(c::Circuit) = Circuit(c.ops; cregs = CReg[])   # 去掉默认经典线
qft3 = nocreg(qft_circuit(3))
nothing
```

```@example ex
print(draw(qft3))
nothing
```

结构即教科书 QFT：对第 `j` 位施加 `H`，再与所有更高位做
受控相位 `CPHASE(π/2^(k-j))`，最后翻转位序。N = 15 求阶时计数寄存器
为 `m = 2n = 8` 位，使用的是同一构造按 8 位展开的同构线路
（`draw(dagger(qft_circuit(8)))` 即读出用的逆变换）。

4 位版本（对应把 15 装入 L = 5 位寄存器的另一常用规模）：

```@example ex
qft4 = nocreg(qft_circuit(4))
nothing
```

```@example ex
print(draw(qft4))
nothing
```

## 3. 完整 QPE 线路的块级视图

`order_finding_circuit` / `order_finding` 门数庞大（N = 15, a = 7 时
FullQFT 版含数千个门），因此构造时把重复子线路封装成了命名 `block`
（`BlockOp`）：受控模加 `CMADD(a·2^k)`、受控模乘 `CMUL(a)`、
交换-消去段 `CSWAP+`、每轮受控 U 幂 `U(a^(2^k))` 以及 b 寄存器的
`QFT_b` / `QFT_b†`。`draw(c; unroll = false)` 把每个块渲染为
单个命名盒，顶层结构一目了然：

```@example ex
cfull = order_finding_circuit(a, N; alg = FullQFT())
nothing
```

```@example ex
print(draw(cfull; unroll = false))
nothing
```

`unroll` / `unroll!` 可手动把块展开为平面门序列（模拟后端会自动
展开，语义与平面版完全一致）；`draw` 默认 `unroll = true`，画的
正是展开后的平面图。

## 4. 打包成一条线路

`order_finding_circuit` 的半经典版（`SemiclassicalQFT`，19 qubit）
进一步内含测量、`IfOp` 经典反馈与输出寄存器——作用到 `|0…0⟩` 上经
`simulate` 执行后，解可直接从输出寄存器（qubit 12..19）读出：

```julia
c = order_finding_circuit(a, N)              # SemiclassicalQFT 打包
out = simulate(c, zero_state(QuantumCircuits.nqubits(c)))
pm = marginal_probabilities(out, collect(2n+4:4n+3))  # 输出寄存器边缘分布
v = argmax(pm) - 1
# v / 2^m ≈ s/r；建议多次采样验证
```

## 5. 换用其他算法

```@example ex
y_slnn = order_finding(a, N; alg = SemiclassicalLNN()) # 2n+4 qubit
y_lnn  = order_finding(a, N; alg = LNNQFT())           # 4n+3 qubit
y_full = order_finding(a, N; alg = FullQFT())          # 4n+2 qubit
```

## 6. 完整因数分解

```@example ex
shor(N)
```

`shor` 内部自动完成：偶数 / 小素数 / 完全幂捷径 → 随机选基与 gcd
捷径 → 量子求阶 → 连分数提取阶 → `gcd(a^(r/2) ∓ 1, N)` 提取因子。
