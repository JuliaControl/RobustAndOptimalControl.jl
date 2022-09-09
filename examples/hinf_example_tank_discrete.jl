using ControlSystemsBase, RobustAndOptimalControl
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
ts = 0.01
Gd = bilinearc2d(ss(A,B,C,D), ts)

# Sensitivity weight function
WSelement = 100*tf([0.1,1],[1000,1])
WS = [WSelement 0; 0 WSelement]
iWSelement = 1/WSelement
iWS = [iWSelement 0; 0 iWSelement]

# Output sensitivity weight function
WUelement = 100*tf([1,1],[0.1,1]) ##
WUelement = ss(0.1)
WU = [WUelement 0; 0 WUelement]
iWUelement = 1/WUelement
iWU = [iWUelement 0; 0 iWUelement]

# Complementary sensitivity weight function
WTelement = tf([10,0.1],[1,1])
WT  = [WTelement 0; 0 WTelement]
iWTelement = 1/WTelement
iWT = [iWTelement 0; 0 iWTelement]

# Create continuous time approximation of the process
Gc = bilineard2c(Gd)

# Form the P in the LFT Fl(P,C) as a partitioned state-space object
Pc = hinfpartition(Gc, WS, WU, WT)

# Check if the problem is feasible
flag = hinfassumptions(Pc)

# Synthesize the H-infinity optimal controller
Cc, γ = hinfsynthesize(Pc)

# Extract the transfer functions defining some signals of interest
Pcl, S, CS, T = hinfsignals(Pc, Gc, Cc)

Pcl = bilinearc2d(Pcl, ts)
S   = bilinearc2d(S,  ts)
CS  = bilinearc2d(CS, ts)
T   = bilinearc2d(T,  ts)

if makeplots
    # Specifications
    specificationplot([S, CS, T], [WSelement, 0.1, WTelement], γ)

    # Closed-loop H-infinity norm
    specificationplot(Pcl, γ; s_labels=["\$\\sigma(P_{cl}(j\\omega))\$"], w_labels=["\$\\gamma\$"])

    # Stepresponse
    times = 0:ts:300
    plot(step(T, times))
end
