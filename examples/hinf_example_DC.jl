using Plots
using ControlSystems
"""
This is a simple SISO example with a pole in the origin, corresponding to the
DC servos used in the Lund laboratories. It serves to exeplify how the syntheis
can be done for simple SISO systems, and also demonstrates how we chan verify
if the problem is feasible to solve using the ARE method.

The example can be set to visualize and save plots using the variables
  MakePlots - true/false (true if plots are to be generated, false for testing)
  SavePlots - true/false (true if plots are to be saved, false for testing)
"""
MakePlots, SavePlots = true, false

# Define the process
Gtrue   = tf([11.2], [1, 0.12,0])

# Sensitivity weight function
M, wB, A = 1.5, 20.0, 1e-8
WS = tf([1/M, wB],[1, wB*A])

# Output sensitivity weight function
WU = ss(1)

# Complementary sensitivity weight function
WT = []

# Form the P in the LFT F_l(P,C) as a partitioned state-space object
P = hinfpartition(Gtrue, WS, WU, WT)

# Check if the system is feasible for synythesis
flag = hinfassumptions(P)

# Since it is not, modify the plant desciption
epsilon = 1e-5
G = tf([11.2], [1, 0.12]) * tf([1], [1, epsilon])

# Form the P in the LFT Fl(P,C) as a partitioned state-space object
P = hinfpartition(G, WS, WU, WT)

# Check if the problem is feasible
flag = hinfassumptions(P)

# Synthesize the H-infinity optimal controller
flag, C, γ = hinfsynthesize(P)

# Extract the transfer functions defining some signals of interest
Pcl, S, CS, T = hinfsignals(P, G, C)

## Plot the specifications
if MakePlots
  specificationplot([S, CS, T], [WS, WU, WT], γ)
  if SavePlots
    savefig("example_DC_specifications.pdf")
  end
end

## Plot the closed loop gain from w to z
if MakePlots
  specificationplot(Pcl, γ; s_labels=["\$\\sigma(P_{cl}(j\\omega))\$"], w_labels=["\$\\gamma\$"])
  if SavePlots
    savefig("example_DC_clgain.pdf")
  end
end
