# ShorAlgorithm.jl

基于 [QuantumCircuits.jl](https://github.com)（量子线路 IR）与
[VQC.jl](https://github.com)（态矢量模拟后端）的 **Shor 算法应用层包**。

把旧 VQC 中的 Shor 算法电路移植为 1-based、小端序的
`QuantumCircuits.Circuit` 中间表示，提供**四种量子求阶算法**与完整的
经典 + 量子混合因数分解。

## 快速上手

```@example index
using ShorAlgorithm
using VQC   # 态矢量模拟后端（加载后自动接入 QuantumCircuits IR）

# 量子求阶：a = 7 模 N = 15 的阶为 r = 4
y = order_finding(7, 15)                     # 相位读数 y ≈ s/r
r = denominator(continued_fraction(y, 15))

r
```

切换算法只需更换 `alg`（派发见 [`order_finding`](@ref)）：

```@example index
order_finding(7, 15; alg = FullQFT())
```

完整因数分解：

```@example index
shor(15)
```

## 求阶算法一览

| 算法类型 | 电路族 | qubit 数 | 说明 |
|---|---|---|---|
| `SemiclassicalQFT()`（默认） | Beauregard | 2n+3 | 单控制比特 + 半经典 QFT（测量反馈） |
| `SemiclassicalLNN()` | Fowler | 2n+4 | 单控制比特 + 半经典 QFT 的 LNN 版 |
| `FullQFT()` | Beauregard | 4n+2 | 完整 QPE |
| `LNNQFT()` | Fowler | 4n+3 | 完整 QPE 的 LNN 版 |

两个半经典版的 IPE 控制位（Beauregard 的 t / LNN 的 ki）都是干净的
控制位：LNN 的 ki 虽在乘法块内部被路由，但每轮结束都会复原
（ki = 0 分支恒等，ki = 1 分支恰为 ×a）。

其中 `n = bits(N)`。日志级别由算法类型携带
（如 `SemiclassicalQFT(verbosity = 1)`）。

## 目录

```@contents
Pages = ["index.md", "examples.md", "api.md"]
Depth = 1
```
