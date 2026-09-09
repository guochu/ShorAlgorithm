using Test
using QuantumCircuits
using QuantumCircuits.Interface: simulate, simulate!
using VQC

# 测试用后端：VQC 的态矢量实现（算法包本身不依赖 VQC）
const backend = VQC.StateVectorBackend()

# ShorAlgorithm 测试入口
include("shor.jl")
