import Pkg

const DOCS = @__DIR__                    # <pkg>/docs
const PKG = dirname(DOCS)                # <pkg>

Pkg.activate(DOCS; io = devnull)
# 本包与工作区依赖尚未注册，按相对路径安装（幂等）
Pkg.develop([
    Pkg.PackageSpec(path = joinpath(dirname(PKG), "QuantumCircuits")),
    Pkg.PackageSpec(path = joinpath(dirname(PKG), "VQC")),
    Pkg.PackageSpec(path = PKG),
]; io = devnull)
Pkg.instantiate(; io = devnull)

using ShorAlgorithm
using Documenter

DocMeta.setdocmeta!(ShorAlgorithm, :DocTestSetup, :(using ShorAlgorithm); recursive = true)

makedocs(;
    modules = [ShorAlgorithm],
    sitename = "ShorAlgorithm.jl",
    authors = "Guo Chu <guochu604b@gmail.com>",
    pages = [
        "主页" => "index.md",
        "示例：N = 15, a = 7" => "examples.md",
        "API 参考" => "api.md",
    ],
)

# 部署到远端仓库时启用（需要 Git 远程信息）：
# deploydocs(;
#     repo = "github.com/<user>/ShorAlgorithm.jl",
#     devbranch = "main",
# )
