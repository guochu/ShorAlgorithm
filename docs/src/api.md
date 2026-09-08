# API 参考

## 求阶

```@docs
order_finding
QuantumOrderFindingAlgorithm
SemiclassicalQFT
SemiclassicalLNN
FullQFT
LNNQFT
order_finding_circuit
continued_fraction
```

## 求阶驱动实现

```@docs
ShorAlgorithm
ShorAlgorithm.order_finding_by_qft
ShorAlgorithm.order_finding_by_lnn
ShorAlgorithm.order_finding_scq
ShorAlgorithm.order_finding_scq_lnn
ShorAlgorithm.order_from_phase
ShorAlgorithm.phase_to_float
ShorAlgorithm._register_bits
ShorAlgorithm._mod_powers
ShorAlgorithm._with_verbosity
```

## 打包线路（内部构造）

```@docs
ShorAlgorithm._semiclassical_qpe_circuit
ShorAlgorithm._semiclassical_lnn_qpe_circuit
ShorAlgorithm._beauregard_qpe_circuit
ShorAlgorithm._lnn_qpe_circuit
```

## 电路工具

```@docs
qft_circuit
sdagger
compute_phase
move_site
shift
controlled_swap_block
```

## 模算术电路积木（Beauregard，全连接）

```@docs
fourier_adder_circuit
controlled_fourier_adder_circuit
double_controlled_fourier_adder_circuit
modular_adder_circuit
controlled_modular_multiplier_circuit
modular_u_circuit
modular_exponentiator_circuit
```

## 模算术电路积木（Fowler，LNN 近邻）

```@docs
lnn_fourier_adder_circuit
lnn_controlled_adder_circuit
lnn_modular_adder_circuit
lnn_controlled_modular_multiplier_circuit
lnn_modular_u_circuit
lnn_modular_exponentiator_circuit
```

## 因数分解

```@docs
shor
is_prime
perfect_power
```
