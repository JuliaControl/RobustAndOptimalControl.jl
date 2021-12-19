ENV["GKSwstype"] = 322 # workaround for gr segfault on GH actions
# ENV["GKS_WSTYPE"]=100 # try this if above does not work
using Documenter, RobustAndOptimalControl

using Plots
gr()


makedocs(
      sitename = "RobustAndOptimalControl Documentation",
      doctest = true,
      modules = [RobustAndOptimalControl],
      pages = [
            "Home" => "index.md",
      ],
      format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
)

deploydocs(
      deps = Deps.pip("pygments", "mkdocs", "python-markdown-math", "mkdocs-cinder"),
      repo = "github.com/JuliaControl/RobustAndOptimalControl.jl.git",
)
