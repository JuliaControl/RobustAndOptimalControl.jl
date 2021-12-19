# RobustAndOptimalControl

[![Build Status](https://github.com/JuliaControl/RobustAndOptimalControl.jl/workflows/CI/badge.svg)](https://github.com/JuliaControl/RobustAndOptimalControl.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaControl/RobustAndOptimalControl.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaControl/RobustAndOptimalControl.jl)



**Work in progress** Expect the functionality and API of this package to be immature and break often.


This package aims to be en experimental testbed for APIs and algorithms which may eventually make their way into ControlSystems.jl

All examples in the folder `examples` are tested in the CI tests, they currently serve as the only available documentation. The tests may further serve as documentation for functionality not yet covered by examples.

Summary of additional functionality:
- `LQGProblem`:
  + Return an LQG object that describes the closed control loop around the process `sys=ss(A,B,C,D)` where the controller is of LQG-type.
  + If only an `ExtendedStateSpace` system is provided, the system `P` is assumed to correspond to the H₂ optimal control problem with ``` C1'C1    = Q1 D12'D12  = Q2 B1*B1'   = R1 D21*D21' = R2 ``` and an `LQGProblem` with the above covariance matrices is returned.
- `NamedStateSpace`: See `named_ss` for a convenient constructor.
- `add_disturbance`: See CCS pp.
- `add_input_integrator`: Augment the output of `sys` with the integral of input at index `ui`, i.e.,  `y_aug = [y; ∫u[ui]]` 
- `add_input_differentiator`: Augment the output of sys with the difference u(k+1)-u(k)
- `add_output_differentiator`: Augment the output of `sys` with the numerical difference (discrete-time derivative) of output, i.e., `y_aug = [y; (y-y_prev)/sys.Ts]` To add both an integrator and a differentiator to a SISO system, use ``` ``` 
- `add_output_integrator`: add_output_integrator(sys::AbstractStateSpace{<:Discrete}, ind = 1; ϵ = 0) Augment the output of `sys` with the integral of output at index `ind`, i.e.,  `y_aug = [y; ∫y[ind]]` To add both an integrator and a differentiator to a SISO system, use ``` Gd = add_output_integrator(add_output_differentiator(G), 1) ``` Note: numerical integration is subject to numerical drift.
- `bilinearc2d`:
  + Balanced Bilinear transformation in State-Space.
  + Applies a Balanced Bilinear transformation to a discrete-time extended statespace object 
  + Applies a Balanced Bilinear transformation to a discrete-time statespace object 
- `bilineard2c`:
  + Balanced Bilinear transformation in State-Space.
  + Applies a Balanced Bilinear transformation to continuous-time extended statespace object 
  + Applies a Balanced Bilinear transformation to continuous-time statespace object 
- `closedloop`: Closed-loop system as defined in Glad and Ljung eq.
- `connect`: Create complicated feedback interconnection. See also [complicated_feedback.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/examples/complicated_feedback.jl).
- `controller_reduction`: Minimize    ||(K-Kᵣ) W||∞ if out=false See Robust and Optimal Control Ch 19.1 out indicates if the weight will be applied as output or input weight.
- `expand_symbol`: Takes a symbol and an integer and returns a vector of symbols with increasing numbers appended to the end.
- `extended_controller`:
  + Returns an expression for the controller that is obtained when state-feedback `u = -L(xᵣ-x̂)` is combined with a Kalman filter with gain `K` that produces state estimates x̂.
  + Takes a controller and returns an `ExtendedStateSpace` version which has augmented input `[r; y]` and output `y` (`z` output is 0-dim).
- `find_lft`: Given an systems `sys` with uncertain coefficients in the form of `StaticParticles`, find a lower linear fractional transformation `M` such that `lft(M, δ) ≈ sys`.
- `frequency_weighted_reduction`: Find Gr such that ||Wₒ(G-Gr)Wᵢ||∞ is minimized.
- `h2synthesize`: Synthesize H₂-optimal controller K and calculate the closed-loop transfer function from `w` to `z`.
- `hinfassumptions`: Check the assumptions for using the γ-iteration synthesis in Theorem 1.
- `hinfpartition`: Transform a SISO or MIMO system G, with weighting functions WS, WU, WT into and LFT with an isolated controller, and write the resulting system, P(s), on a state-space form.
- `hinfsignals`: Use the extended state-space model, a plant and the found controller to extract the closed loop transfer functions operating solely on the state-space.
- `hinfsynthesize`: Computes an H-infinity optimal controller K for an extended plant P such that ||F_l(P, K)||∞ < γ for the largest possible γ given P.
- `hsvd`: Return the Hankel singular values of `sys`, computed as the eigenvalues of `QP` Where `Q` and `P` are the Gramians of `sys`.
- `makeweight`: Create a weighting function that goes from gain `low` at zero frequency, through gain `mid` to gain `high` at ∞ # Arguments: - `low`: A number specifying the DC gain  - `mid`: A number specifying the frequency at which the gain is 1, or a tuple `(freq, gain)`.
- `measure`: Return a system with specified states as measurement outputs.
- `named_ss`: Create a `NamedStateSpace` system. StateSpace systems with named inputs, outputs and states. See also [complicated_feedback.jl](https://github.com/JuliaControl/RobustAndOptimalControl.jl/blob/master/examples/complicated_feedback.jl).
- `specificationplot`: This function visualizes the control synthesis using the hInf_synthesize with the three weighting functions {WS(jω), WU(jω), WT(jω)} inverted and scaled by γ, against the corresponding transfer fucntions {S(jω), C(jω)S(jω), T(jω)}, to verify visually that the specifications are met.
- `sumblock`: Create a summation node that sums (or subtracts) vectors of length `n`.
- `δ`: Create an uncertain element of `N` uniformly distributed samples ∈ [-1, 1] 
- `ControlSystems.feedback`: Feedback between two named systems.
- `ControlSystems.observer_controller`: Returns an expression for the feedback controller `u = Cy` that is obtained when state-feedback `u = -Lx̂` is combined with a Kalman filter with gain `K` that produces state estimates x̂.
- `RobustAndOptimalControl._assertrealandpsd`: Check that a matrix is real and PSD - throw an error otherwise.
- `RobustAndOptimalControl._checkfeasibility`: Check the feasibility of the computed solutions Xinf, Yinf and the algebraic Riccatti equations, return true if the solution is valid, and false otherwise.
- `RobustAndOptimalControl._coordinatetransformqr`: Use the QR decomposition to find a transformaiton [Tl, Tr] such that Tl*A*Tr becomes [0;I], [0 I] or I depending on the dimensionality of A.
- `RobustAndOptimalControl._coordinatetransformsvd`: Use the SVD to find a transformaiton [Tl, Tr] such that Tl*A*Tr becomes [0;I], [0 I] or I depending on the dimensionality of A.
- `RobustAndOptimalControl._detectable`: Applies the Hautus lemma to check if the pair is detectable 
- `RobustAndOptimalControl._input2ss`: Helper function used for type conversion in hinfpartition() 
- `RobustAndOptimalControl._scalematrix`: Find a left and right transform of A such that Tl*A*Tr = [I, 0], or Tl*A*Tr = [I; 0], depending on the dimensionality of A.
- `RobustAndOptimalControl._solvehamiltonianare`: Solves a hamiltonian Alebraic Riccati equation using the Schur-decomposition, for additional details, see   @article{laub1979schur,   } 
- `RobustAndOptimalControl._solvematrixequations`: Solves the dual matrix equations in the γ-iterations (equations 7-12 in Doyle).
- `RobustAndOptimalControl._stabilizable`: Applies the Hautus lemma to check if the pair is stabilizable 
- `RobustAndOptimalControl._synthesizecontroller`: Syntheize a controller by operating on the scaled state-space description of the system (i.e., the state-space realization of `P̄`) using the solutions from the γ-iterations.
- `RobustAndOptimalControl._transformp2pbar`: Transform the original system P to a new system Pbar, in which D12bar = [0; I] and D21bar = [0 I] in order to satisfy the feasibility assumption A3 (see Doyle) 
- `RobustAndOptimalControl._γiterations`: Rune the complete set of γ-iterations over a specified search interval with a set number of iterations.
- `RobustAndOptimalControl.controller_reduction_weight`: Lemma 19.1 See Robust and Optimal Control Ch 19.1 


## Installation
```julia
pkg> add RobustAndOptimalControl
```
## Acknowledgement
The code for the H∞ design was originally written by Marcus Greiff @mgreiff in the PR https://github.com/JuliaControl/ControlSystems.jl/pull/192 