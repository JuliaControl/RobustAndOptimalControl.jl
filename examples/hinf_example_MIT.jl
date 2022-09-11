using ControlSystemsBase, RobustAndOptimalControl
using Plots

"""
This is a simple SISO example which was used for debugging the implementation,
as this exact example was use in the lecture notes of the "Principles of Optimal
Control" cours of the MIT OpenCourseWare [1], where the choice of weighting
functions and dynamics gave rise to an H-infinity optimal cotnroller with a
γ of approximately 1.36, where, in our case, we get a controller at 0.93

[1] https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-323-principles-of-optimal-control-spring-2008/lecture-notes/lec15.pdf
"""

# Define the process
G   = tf([200], [0.025,1.0025,10.1,1])

# Sensitivity weight function
M, wB, A = 1.5, 10, 1e-4
WS = tf([1/M, wB],[1, wB*A])

# Output sensitivity weight function
WU = ss(0.1)

# Complementary sensitivity weight function
WT = []

# Form augmented P dynamics in state-space
P = hinfpartition(G, WS, WU, WT)

# Check that the assumptions are satisfied
flag = hinfassumptions(P)

# Synthesize the H-infinity optimal controller
C, γ = hinfsynthesize(P, γrel=1)

# Extract the transfer functions defining some signals of interest
Pcl, S, CS, T = hinfsignals(P, G, C)

## Plot the specifications
plot(
  specificationplot([S, CS, T], [ss(WS), WU, WT], γ),
## Plot the closed loop gain from w to z
  specificationplot(Pcl, γ; s_labels=["\$\\sigma(P_{cl}(j\\omega))\$"], w_labels=["\$\\gamma\$"])
) |> display
