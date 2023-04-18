ENV["GKSwstype"] = 322 # workaround for gr segfault on GH actions
# ENV["GKS_WSTYPE"]=100 # try this if above does not work
using Documenter, RobustAndOptimalControl, ControlSystemsBase

# Make it possible to load the CollapsedExample Documenter extension
let exts = joinpath(@__DIR__, "ext")
      if !(exts in LOAD_PATH)
            pushfirst!(LOAD_PATH, exts)
      end
end
import DocumenterCollapsedExample

using Plots
gr()


makedocs(
      sitename = "RobustAndOptimalControl Documentation",
      doctest = true,
      modules = [RobustAndOptimalControl, ControlSystemsBase],
      pagesonly = false,
      pages = [
            "Home" => "index.md",
            "Uncertainty modeling" => "uncertainty.md",
            "Examples" => [
                  "Simple mixed-sensitivity ``H_\\infty`` design" => "hinf_DC.md",
                  "General ``H_\\infty`` design" => "hinf_connection.md",
                  "LQG control with integral action" => "lqg_disturbance.md",
                  "Control design for a pendulum on a cart" => "cartpole.md",
                  "``H_2`` design with a known disturbance" => "measurable_disturbance.md",
                  "When are systems similar?" => "similarity.md",
            ],
            "API" => "api.md",
      ],
      format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
)

deploydocs(
      deps = Deps.pip("pygments", "mkdocs", "python-markdown-math", "mkdocs-cinder"),
      repo = "github.com/JuliaControl/RobustAndOptimalControl.jl.git",
)
