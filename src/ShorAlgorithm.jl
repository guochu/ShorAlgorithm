"""
    ShorAlgorithm

基于 `QuantumCircuits` IR 的 Shor 算法应用层包：把 VQC（旧 0-based
QCircuit 实现）中的 Shor 算法电路移植为 1-based、小端序的
`QuantumCircuits.Circuit`，并用 `VQC` 态矢量后端执行。

分层：

* **电路工具**：`move_site`、`shift`、`compute_phase`、`qft_circuit`、
  `controlled_swap_block`；模算术直接使用 Base 的 `mod` / `invmod` /
  `powermod`，寄存器位宽统一取 `L = bits(N) + 1`（保证 `2N ≤ 2^L`）；
* **电路积木（Beauregard，全连接）**：Fourier 加法器、（双）受控模加法器、
  受控模乘法器、模乘 U 块、模指电路（`modular_exponentiator_circuit`）；
* **电路积木（Fowler，LNN 近邻）**：`lnn_*` 系列对应电路；
* **算法**：`order_finding`（唯一入口，`alg` 关键字按
  `QuantumOrderFindingAlgorithm` 的子类型单例派发：
  `SemiclassicalQFT`（半经典 Beauregard，默认）/ `SemiclassicalLNN`
  （半经典 Fowler 近邻）/ `FullQFT`（Beauregard 完整 QPE）/
  `LNNQFT`（Fowler 近邻完整 QPE），日志级别由 `alg` 携带）、
  `order_finding_circuit`（把整个求阶计算打包成一条线路——半经典版
  含测量 / 经典反馈 IfOp / 输出寄存器）、`shor`（完整因数分解）。

约定（与 QuantumCircuits / VQC 一致）：

* 比特索引 **1-based**，状态小端序（qubit 1 = 最低有效位）；
* 门矩阵：`GateOp.qubits` 列表第一个比特为矩阵最高位；
* 各算术寄存器沿用 Beauregard/Fowler 论文的布局——**寄存器内部
  首个（号码最小的）qubit 为最高有效位**，见各电路 docstring 的布局图。
"""
module ShorAlgorithm

using Random
using QuantumCircuits
using VQC

export compute_phase,
       continued_fraction, is_prime, perfect_power,
       # 电路工具
       move_site, shift, sdagger, qft_circuit, controlled_swap_block,
       # Beauregard（全连接）
       fourier_adder_circuit,
       controlled_fourier_adder_circuit,
       double_controlled_fourier_adder_circuit,
       modular_adder_circuit,
       controlled_modular_multiplier_circuit,
       modular_u_circuit,
       modular_exponentiator_circuit,
       order_finding_circuit,
       # Fowler（LNN 近邻）
       lnn_fourier_adder_circuit,
       lnn_controlled_adder_circuit,
       lnn_modular_adder_circuit,
       lnn_controlled_modular_multiplier_circuit,
       lnn_modular_u_circuit,
       lnn_modular_exponentiator_circuit,
       # 算法
       QuantumOrderFindingAlgorithm, SemiclassicalQFT, SemiclassicalLNN, FullQFT, LNNQFT,
       order_finding,
       shor

include("util.jl")
include("beauregard.jl")
include("fowler.jl")
include("orderfinding.jl")
include("shor.jl")

end # module
