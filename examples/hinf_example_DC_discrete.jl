using Plots
using ControlSystems, RobustAndOptimalControl
"""
This is a simple SISO example with a pole in the origin, corresponding to the
DC servos used in the Lund laboratories. It serves to exeplify how the syntheis
can be done for simple SISO systems, and also demonstrates how we chan verify
if the problem is feasible to solve using the ARE method.

The example can be set to visualize and save plots using the variables
  makeplots - true/false (true if plots are to be generated, false for testing)
  SavePlots - true/false (true if plots are to be saved, false for testing)
"""
makeplots = true

# Define the process
ts = 0.01
epsilon = 1e-5
Gd = ss(c2d(tf([11.2], [1, 0.12]) * tf([1], [1, epsilon]), ts))

# Sensitivity weight function
M, wB, A = 1.5, 20.0, 1e-8
WS = tf([1/M, wB],[1, wB*A])

# Output sensitivity weight function
WU = ss(1.0)

# Complementary sensitivity weight function
WT = []

# Create continuous time approximation of the process
Gc = bilineard2c(ss(Gd))

# Form the P in the LFT Fl(P,C) as a partitioned state-space object
Pc = hinfpartition(Gc, WS, WU, WT)

# Check if the problem is feasible
flag = hinfassumptions(Pc)

# Synthesize the H-infinity optimal controller
flag, Cc, γ = hinfsynthesize(Pc)

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
