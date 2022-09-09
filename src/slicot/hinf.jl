using SLICOT_jll
slicot = SLICOT_jll.libslicot_path

G = ssrand(2,3,4);
P = partition(G, 1, 1)


# SUBROUTINE SB10AD( JOB, N, M, NP, NCON, NMEAS, GAMMA, A, LDA,
#                         B, LDB, C, LDC, D, LDD, AK, LDAK, BK, LDBK, CK,
#                         LDCK, DK, LDDK, AC, LDAC, BC, LDBC, CC, LDCC,
#                         DC, LDDC, RCOND, GTOL, ACTOL, IWORK, LIWORK,
#                         DWORK, LDWORK, BWORK, LBWORK, INFO )


# JOB     (input) INTEGER
# Indicates the strategy for reducing the GAMMA value, as
# follows:
# = 1: Use bisection method for decreasing GAMMA from GAMMA
#      to GAMMAMIN until the closed-loop system leaves
#      stability.
# = 2: Scan from GAMMA to 0 trying to find the minimal GAMMA
#      for which the closed-loop system retains stability.
# = 3: First bisection, then scanning.
# = 4: Find suboptimal controller only.

job = 1

# N       (input) INTEGER
# The order of the system.  N >= 0.
N = P.nx

# M       (input) INTEGER
# The column size of the matrix B.  M >= 0.
M = P.nu + P.nw

# NP      (input) INTEGER
# The row size of the matrix C.  NP >= 0.
NP = P.ny + P.nz

# NCON    (input) INTEGER
# The number of control inputs (M2).  M >= NCON >= 0,
# NP-NMEAS >= NCON.
NCON = P.nu
M2 = NCON

# NMEAS   (input) INTEGER
# The number of measurements (NP2).  NP >= NMEAS >= 0,
# M-NCON >= NMEAS.
NMEAS = P.ny
NP2 = NMEAS

# GAMMA   (input/output) DOUBLE PRECISION
# The initial value of gamma on input. It is assumed that
# gamma is sufficiently large so that the controller is
# admissible. GAMMA >= 0.
# On output it contains the minimal estimated gamma.
GAMMA = 20

# A       (input) DOUBLE PRECISION array, dimension (LDA,N)
# The leading N-by-N part of this array must contain the
# system state matrix A.
A = P.A

# LDA     INTEGER
# The leading dimension of the array A.  LDA >= max(1,N).
LDA = N

# B       (input) DOUBLE PRECISION array, dimension (LDB,M)
# The leading N-by-M part of this array must contain the
# system input matrix B.
B = P.B

# LDB     INTEGER
# The leading dimension of the array B.  LDB >= max(1,N).
LDB = N

# C       (input) DOUBLE PRECISION array, dimension (LDC,N)
# The leading NP-by-N part of this array must contain the
# system output matrix C.
C = P.C

# LDC     INTEGER
# The leading dimension of the array C.  LDC >= max(1,NP).
LDC = NP

# D       (input) DOUBLE PRECISION array, dimension (LDD,M)
# The leading NP-by-M part of this array must contain the
# system input/output matrix D.
D = P.D

# LDD     INTEGER
# The leading dimension of the array D.  LDD >= max(1,NP).
LDD = NP

# AK      (output) DOUBLE PRECISION array, dimension (LDAK,N)
# The leading N-by-N part of this array contains the
# controller state matrix AK.
AK = zeros(N,N)

# LDAK    INTEGER
# The leading dimension of the array AK.  LDAK >= max(1,N).
LDAK = N

# BK      (output) DOUBLE PRECISION array, dimension (LDBK,NMEAS)
# The leading N-by-NMEAS part of this array contains the
# controller input matrix BK.
BK = zeros(N,NMEAS)

# LDBK    INTEGER
# The leading dimension of the array BK.  LDBK >= max(1,N).
LDBK = N

# CK      (output) DOUBLE PRECISION array, dimension (LDCK,N)
# The leading NCON-by-N part of this array contains the
# controller output matrix CK.
CK = zeros(NCON,N)

# LDCK    INTEGER
# The leading dimension of the array CK.
# LDCK >= max(1,NCON).
LDCK = NCON

# DK      (output) DOUBLE PRECISION array, dimension (LDDK,NMEAS)
# The leading NCON-by-NMEAS part of this array contains the
# controller input/output matrix DK.
DK = zeros(NCON,NMEAS)

# LDDK    INTEGER
# The leading dimension of the array DK.
# LDDK >= max(1,NCON).
LDDK = NCON

# AC      (output) DOUBLE PRECISION array, dimension (LDAC,2*N)
# The leading 2*N-by-2*N part of this array contains the
# closed-loop system state matrix AC.
AC = zeros(2N, 2N)

# LDAC    INTEGER
# The leading dimension of the array AC.
# LDAC >= max(1,2*N).
LDAC = 2N

# BC      (output) DOUBLE PRECISION array, dimension (LDBC,M-NCON)
# The leading 2*N-by-(M-NCON) part of this array contains
# the closed-loop system input matrix BC.
BC = zeros(2N, M-NCON)

# LDBC    INTEGER
# The leading dimension of the array BC.
# LDBC >= max(1,2*N).
LDBC = 2N

# CC      (output) DOUBLE PRECISION array, dimension (LDCC,2*N)
# The leading (NP-NMEAS)-by-2*N part of this array contains
# the closed-loop system output matrix CC.
CC = zeros(NP-NMEAS, 2N)

# LDCC    INTEGER
# The leading dimension of the array CC.
# LDCC >= max(1,NP-NMEAS).
LDCC = NP-NMEAS

# DC      (output) DOUBLE PRECISION array, dimension (LDDC,M-NCON)
# The leading (NP-NMEAS)-by-(M-NCON) part of this array
# contains the closed-loop system input/output matrix DC.
DC = zeros(NP-NMEAS, M-NCON)

# LDDC    INTEGER
# The leading dimension of the array DC.
# LDDC >= max(1,NP-NMEAS).
LDDC = NP-NMEAS

# RCOND   (output) DOUBLE PRECISION array, dimension (4)
#          For the last successful step:
# RCOND(1) contains the reciprocal condition number of the
#          control transformation matrix;
# RCOND(2) contains the reciprocal condition number of the
#          measurement transformation matrix;
# RCOND(3) contains an estimate of the reciprocal condition
#          number of the X-Riccati equation;
# RCOND(4) contains an estimate of the reciprocal condition
#          number of the Y-Riccati equation.

RCOND = zeros(4)

# Tolerances
# GTOL    DOUBLE PRECISION
# Tolerance used for controlling the accuracy of GAMMA
# and its distance to the estimated minimal possible
# value of GAMMA.
# If GTOL <= 0, then a default value equal to sqrt(EPS)
# is used, where EPS is the relative machine precision.
GTOL = 1e-3

# ACTOL   DOUBLE PRECISION
# Upper bound for the poles of the closed-loop system
# used for determining if it is stable.
# ACTOL <= 0 for stable systems.
ACTOL = 0

# Workspace
# IWORK   INTEGER array, dimension (LIWORK)

# LIWORK  INTEGER
# The dimension of the array IWORK.
LIWORK = max(2*max(N,M-NCON,NP-NMEAS,NCON,NMEAS),N*N)
IWORK = zeros(Int, LIWORK)

# DWORK   DOUBLE PRECISION array, dimension (LDWORK)
# On exit, if INFO = 0, DWORK(1) contains the optimal
# value of LDWORK.

# LDWORK  INTEGER

# where
M1  = M   - M2
NP1 = NP - NP2
ND1 = NP1 - M2
ND2 = M1 - NP2
LW1 = N*M + NP*N + NP*M + M2*M2 + NP2*NP2;
LW2 = max( ( N + NP1 + 1 )*( N + M2 ) +
             max( 3*( N + M2 ) + N + NP1, 5*( N + M2 ) ),
           ( N + NP2 )*( N + M1 + 1 ) +
             max( 3*( N + NP2 ) + N + M1, 5*( N + NP2 ) ),
           M2 + NP1*NP1 + max( NP1*max( N, M1 ),
                               3*M2 + NP1, 5*M2 ),
           NP2 + M1*M1 +  max( max( N, NP1 )*M1,
                               3*NP2 + M1, 5*NP2 ) );
LW3 = max( ND1*M1 + max( 4*min( ND1, M1 ) + max( ND1,M1 ),
                         6*min( ND1, M1 ) ),
           NP1*ND2 + max( 4*min( NP1, ND2 ) +
                                           max( NP1,ND2 ),
                          6*min( NP1, ND2 ) ) );
LW4 = 2*M*M + NP*NP + 2*M*N + M*NP + 2*N*NP;
LW5 = 2*N*N + M*N + N*NP;
LW6 = max( M*M   + max( 2*M1, 3*N*N +
                        max( N*M, 10*N*N + 12*N + 5 ) ),
           NP*NP + max( 2*NP1, 3*N*N +
                        max( N*NP, 10*N*N + 12*N + 5 ) ));
LW7 = M2*NP2 + NP2*NP2 + M2*M2 +
      max( ND1*ND1 + max( 2*ND1, ( ND1 + ND2 )*NP2 ),
           ND2*ND2 + max( 2*ND2, ND2*M2 ), 3*N,
           N*( 2*NP2 + M2 ) +
           max( 2*N*M2, M2*NP2 +
                        max( M2*M2 + 3*M2, NP2*( 2*NP2 +
                             M2 + max( NP2, N ) ) ) ) );



# The dimension of the array DWORK.
LDWORK = LW1 + max(1,LW2,LW3,LW4,LW5 + max(LW6,LW7))
LDWORK *= 2 # For good performance, LDWORK must generally be larger.
DWORK = zeros(LDWORK)

# BWORK   LOGICAL array, dimension (LBWORK)

# LBWORK  INTEGER
# The dimension of the array BWORK.  LBWORK >= 2*N.
LBWORK = 2N
BWORK = falses(LBWORK)

INFO = Ref(0)

AF64 = Vector{Float64}
ABool = Vector{Bool}

@ccall SB10AD( JOB::Int, N::Int, M::Int, NP::Int, NCON::Int, NMEAS::Int, GAMMA::Float64, A::AF64, LDA::Int,
                            B::AF64, LDB::Int, C::AF64, LDC::Int, D::AF64, LDD::Int, AK::AF64, LDAK::Int, BK::AF64, LDBK::Int, CK::AF64,
                            LDCK::Int, DK::AF64, LDDK::Int, AC::AF64, LDAC::Int, BC::AF64, LDBC::Int, CC::AF64, LDCC::Int,
                            DC::AF64, LDDC::Int, RCOND::AF64, GTOL::Float64, ACTOL::Float64, IWORK, LIWORK::Int,
                            DWORK::AF64, LDWORK::Int, BWORK::ABool, LBWORK::Int, INFO::Int )