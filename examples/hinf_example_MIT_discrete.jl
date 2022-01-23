using ControlSystems, RobustAndOptimalControl
using Plots

"""
This is a simple SISO example which was used for debugging the implementation,
as this exact example was use in the lecture notes of the "Principles of Optimal
Control" cours of the MIT OpenCourseWare [1], however, this example is in discrete time

[1] https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-323-principles-of-optimal-control-spring-2008/lecture-notes/lec15.pdf

The example can be set to visualize and save plots using the two variables
  makeplots - true/false (true if plots are to be generated, false for testing)
"""
makeplots = true

# Define the process
ts = 0.005
Gd = bilinearc2d(ss(tf([200], [0.025,1.0025,10.1,1])),ts)

# Sensitivity weight function
M, wB, A = 1.5, 10, 1e-4
WS = tf([1/M, wB],[1, wB*A])

# Output sensitivity weight function
WU = ss(0.1)

# Complementary sensitivity weight function
WT = []

# Create continuous time approximation of the process
Gc = bilineard2c(ss(Gd))

# Form augmented P dynamics in state-space
Pc = hinfpartition(Gc, WS, WU, WT)

# Check that the assumptions are satisfied
flag = hinfassumptions(Pc)

# Synthesize the H-infinity optimal controller
Cc, γ = hinfsynthesize(Pc)

# Extract the transfer functions defining some signals of interest, but do so
# using discrete equivalent of the continuous time objects Pc, Cc and Gc
Pcl, S, CS, T = hinfsignals(
  bilinearc2d(Pc, ts),
  bilinearc2d(Gc, ts),
  bilinearc2d(Cc, ts)
)

# Visualize results
if makeplots
  specificationplot([S, CS, T], [WS, WU, WT], γ)
  specificationplot(Pcl, γ; s_labels=["\$\\sigma(P_{cl}(j\\omega))\$"], w_labels=["\$\\gamma\$"])
end
