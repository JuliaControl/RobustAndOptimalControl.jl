using ControlSystems, RobustAndOptimalControl
using Plots
using LinearAlgebra
"""
This is a simple SISO example with integrator dynamics corresponding to the
quad tank process in the lab.

The example can be set to visualize and save plots using the variables
  makeplots - true/false (true if plots are to be generated, false for testing)
  SavePlots - true/false (true if plots are to be saved, false for testing)
"""
makeplots = true

# Define the proces parameters
k1, k2, kc, g = 3.33, 3.35, 0.5, 981
A1, A3, A2, A4 = 28, 28, 32, 32
a1, a3, a2, a4= 0.071, 0.071, 0.057, 0.057
h01, h02, h03, h04 = 12.4, 12.7, 1.8, 1.4
T1, T2 = (A1/a1)*sqrt(2*h01/g), (A2/a2)*sqrt(2*h02/g)
T3, T4 = (A3/a3)*sqrt(2*h03/g), (A4/a4)*sqrt(2*h04/g)
c1, c2 = (T1*k1*kc/A1), (T2*k2*kc/A2)
γ1, γ2 = 0.7, 0.6

# Define the process dynamics
A = [-1/T1     0 A3/(A1*T3)          0;
     0     -1/T2          0 A4/(A2*T4);
     0         0      -1/T3          0;
     0         0          0      -1/T4];
B = [γ1*k1/A1     0;
     0                γ2*k2/A2;
     0                (1-γ2)*k2/A3;
     (1-γ1)*k1/A4 0              ];
C = [kc 0 0 0;
     0 kc 0 0];
D = zeros(2,2)
G = ss(A,B,C,D);

# Sensitivity weight function
WS = makeweight(100, (0.1, 1), 1/10) * I(2)

# Output sensitivity weight function
# WUelement = 5*tf([1,1],[0.1,1])
WU = tf(0.01) .* I(2)

# Complementary sensitivity weight function
WT = makeweight(1/10, (0.1, 1), 10) * I(2)

# Form augmented P dynamics in state-space
P = hinfpartition(G, WS, WU, WT)

# Check that the assumptions are satisfied
flag = hinfassumptions(P)

# Synthesize the H-infinity optimal controller
C, γ = hinfsynthesize(P)

Pcl, S, CS, T = hinfsignals(P, G, C)

## Plot the specifications
# TODO figure out why I get segmentation errors when using ss instead of tf for
# the weighting functions, makes no sense at all
if makeplots
  specificationplot([S, CS, T], [WS[1,1], 0.01, WT[1,1]], γ)

## Plot the closed loop gain from w to z
# TODO figure out why the legends don't seem to work in this case

  specificationplot(Pcl, γ; s_labels=["\$\\sigma(P_{cl}(j\\omega))\$"], w_labels=["\$\\gamma\$"])

  times = [i for i in range(0, stop=300, length=10000)]
  plot(step(T, times))
end
