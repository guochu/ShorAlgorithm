using Test
using QuantumCircuits
using QuantumCircuits.Interface: simulate, simulate!
using VQC
using MPSSimulator

# 两个后端包都导出 apply! / probabilities / marginal_probabilities（实现各异），
# 显式消歧：纯线路测试组固定用 VQC 的态矢量版本；MPS 组仅经 Interface 驱动。
using VQC: apply!, probabilities, marginal_probabilities

# 测试用后端（算法包本身不依赖任何模拟器）：
#   * VQC 的态矢量实现；
#   * MPSSimulator 的矩阵乘积态实现（跨后端验证 Interface 契约）。
const backend = VQC.StateVectorBackend()
const mps_backend = MPSSimulator.MatrixProductBackend()

# ShorAlgorithm 测试入口
include("shor.jl")
