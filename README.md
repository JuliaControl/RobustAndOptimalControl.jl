# RobustAndOptimalControl

[![Build Status](https://github.com/JuliaControl/RobustAndOptimalControl.jl/workflows/CI/badge.svg)](https://github.com/JuliaControl/RobustAndOptimalControl.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaControl/RobustAndOptimalControl.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaControl/RobustAndOptimalControl.jl)



**Work in progress** Expect the functionality and API of this package to be immature and break often.


This package aims to be en experimental testbed for APIs and algorithms which may eventually make their way into ControlSystems.jl

All examples in the folder `examples` are tested in the CI tests, they currently serve as the only available documentation. The tests may further serve as documentation for functionality not yet covered by examples.


## Installation
```julia
pkg> add RobustAndOptimalControl
```
## Acknowledgement
The code for the H∞ design was originally written by Marcus Greiff @mgreiff in the PR https://github.com/JuliaControl/ControlSystems.jl/pull/192 