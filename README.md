# ShorAlgorithm.jl

基于 [QuantumCircuits.jl](../QuantumCircuits)（量子线路 IR）与
[VQC.jl](../VQC)（态矢量模拟后端）的 **Shor 算法应用层包**。
把旧 VQC 中的 Shor 算法电路移植为 1-based、小端序的
`QuantumCircuits.Circuit`，提供四种量子求阶算法与完整的经典 + 量子
混合因数分解。

## 特性

* **四种求阶算法**，经统一入口 `order_finding` 的 `alg` 关键字派发：

  | 算法类型 | 电路族 | qubit 数 | 说明 |
  |---|---|---|---|
  | `SemiclassicalQFT()`（默认） | Beauregard | 2n+3 | 单控制比特 + 半经典 QFT（测量反馈） |
  | `SemiclassicalLNN()` | Fowler | 2n+4 | 半经典 QFT 的线性近邻（LNN）版 |
  | `FullQFT()` | Beauregard | 4n+2 | 完整 QPE |
  | `LNNQFT()` | Fowler | 4n+3 | 完整 QPE 的 LNN 版 |

  其中 `n = bits(N)`。日志级别由算法类型携带
  （如 `SemiclassicalQFT(verbosity=1)`）。

* **线路打包**：`order_finding_circuit(a, N; alg)` 把整个求阶计算
  （半经典版含测量、`IfOp` 经典反馈与输出寄存器）包装成一条
  `Circuit`，作用到量子态上即得解。
* **电路积木**：QFT、Fourier 域加法器、（双）受控模加法器、受控模
  乘法器、模乘 U 块、模指电路——均有全连接（Beauregard）与 LNN
  （Fowler）两族，可独立复用。
* **完整因数分解** `shor`：试除 / 完全幂 / gcd 捷径 + 量子求阶 +
  连分数后处理。

## 安装

本包尚未注册；在项目环境中以路径安装本包及其依赖：

```julia
using Pkg
Pkg.develop([
    PackageSpec(path = "../QuantumCircuits"),
    PackageSpec(path = "../VQC"),
    PackageSpec(path = "."),
])
```

## 快速上手

```julia
using ShorAlgorithm
using VQC   # 态矢量模拟后端（加载后自动接入 QuantumCircuits IR）

# 量子求阶：a = 7 模 N = 15 的阶为 r = 4
y = order_finding(7, 15)                     # 相位读数 y ≈ s/r
r = denominator(continued_fraction(y, 15))   # → 4

# 切换算法
y = order_finding(7, 15; alg = FullQFT())
y = order_finding(7, 15; alg = LNNQFT())

# 完整因数分解
shor(15)   # → (3, 5)
```

把整个求阶计算打包成一条线路（半经典版含测量 / 经典反馈）：

```julia
using QuantumCircuits, VQC

c = order_finding_circuit(7, 15)      # SemiclassicalQFT 打包：19 qubit
out = simulate(c, zero_state(19))
# 解编码在输出寄存器 out = 2n+4..4n+3（此处 12..19）：
pm = marginal_probabilities(out, collect(12:19))
```

## 文档

```julia
julia --project=docs docs/make.jl
```

构建后的文档位于 `docs/build/`；`docs/src/examples.md` 含
`N = 15, a = 7` 的完整走查（含 QFT 线路图）。

## 测试

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

## 参考

* P. W. Shor, *Polynomial-Time Algorithms for Prime Factorization and
  Discrete Logarithms on a Quantum Computer*, SIAM J. Comput. 26, 1484 (1997).
* S. Beauregard, *Circuit for Shor's algorithm using 2n+3 qubits*,
  Quant. Inf. Comput. 3, 175 (2003).
* R. R. Griffiths, C.-S. Niu, *Semiclassical Fourier Transform for
  Quantum Computation*, PRL 76, 3228 (1996).
* A. G. Fowler, S. J. Devitt, L. C. L. Hollenberg,
  *Implementation of Shor's Algorithm on a Linear Nearest Neighbour
  Qubit Array*, Quant. Inf. Comput. 4, 237 (2004).
