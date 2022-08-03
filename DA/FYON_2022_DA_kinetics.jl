#=
This file contains all DA model gating functions
=#

# Gating functions
boltz(V, A, B) = 1 / (1 + exp(-(V-A)/B))
tauX(V, A, B, D, E) = A - B / (1 + exp((V+D)/E))

# Na-current (m = activation variable, h = inactivation variable)
m_inf(V) = boltz(V, -30.0907, 9.7264)
tau_m(V) = 0.01 + 1.0 / ((-(15.6504 + 0.4043*V)/(exp(-19.565 -0.5052*V)-1.0)) + 3.0212*exp(-7.4630e-3*V))
h_inf(V) = boltz(V, -54.0289, -10.7665)
tau_h(V) = 0.4 + 1.0 / ((5.0754e-4*exp(-6.3213e-2*V)) + 9.7529*exp(0.13442*V))

# Kd-current (n = activation variable)
n_inf(V) = boltz(V, -25., 12.)
tau_n(V) = tauX(V, 20., 18., 38., -10.)

# L-type Ca-current (m = activation variable)
mCaL_inf(V) = boltz(V, -50., 2.)
tau_mCaL(V) = tauX(V, 30., 28., 45., -3.)

# N-type Ca-current (m = activation variable)
mCaN_inf(V) = boltz(V, -30., 7.)
tau_mCaN(V) = tauX(V, 30., 25., 55., -6.)

# ERG K-current (o = activation variable, i = intermediate variable)
a0ERG(V) = 0.0036 * exp(0.0759*V)
b0ERG(V) = 1.2523e-5 * exp(-0.0671*V)
aiERG(V) = 0.1 * exp(0.1189*V)
biERG(V) = 0.003 * exp(-0.0733*V)
o_inf(V) = a0ERG(V)*biERG(V) / (a0ERG(V)*(aiERG(V)+biERG(V)) + b0ERG(V)*biERG(V))
i_inf(V) = a0ERG(V)*aiERG(V) / (a0ERG(V)*(aiERG(V)+biERG(V)) + b0ERG(V)*biERG(V))

# NMDA current
NMDA_inf(V) = 1 / (1 + Mg*exp(-0.08*V)/10.)
