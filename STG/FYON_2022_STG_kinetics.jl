#=
This file contains all STG model gating functions
=#

# Gating functions
boltz(V, A, B) = 1 / (1 + exp((V+A)/B))
tauX(V, A, B, D, E) = A - B / (1+exp((V+D)/E))

# Na-current (m = activation variable, h = inactivation variable)
mNa_inf(V) = boltz(V, 25.5, -5.29)
tau_mNa(V) = tauX(V, 1.32, 1.26, 120., -25.)
hNa_inf(V) = boltz(V, 48.9, 5.18)
tau_hNa(V) = (0.67 / (1 + exp((V+62.9)/-10.0))) * (1.5 + 1 / (1+exp((V+34.9)/3.6)))

# T-type Ca-current (m = activation variable, h = inactivation variable)
mCaT_inf(V) = boltz(V, 27.1, -7.2)
tau_mCaT(V) = tauX(V, 21.7, 21.3, 68.1, -20.5)
hCaT_inf(V) = boltz(V, 32.1, 5.5)
tau_hCaT(V) = tauX(V, 105., 89.8, 55., -16.9)

# Slow Ca-current (m = activation variable, h = inactivation variable)
mCaS_inf(V) = boltz(V, 33., -8.1)
tau_mCaS(V) = 1.4 + (7 / ((exp((V+27)/10)) + (exp((V+70)/-13))))
hCaS_inf(V) = boltz(V, 60., 6.2)
tau_hCaS(V) = 60 + (150 / ((exp((V+55)/9)) + (exp((V+65)/-16))))

# A-current (m = activation variable, h = inactivation variable)
mA_inf(V) = boltz(V, 27.2, -8.7)
tau_mA(V) = tauX(V, 11.6, 10.4, 32.9, -15.2)
hA_inf(V) = boltz(V, 56.9, 4.9)
tau_hA(V) = tauX(V, 38.6, 29.2, 38.9, -26.5)

# Ca K-current (m = activation variable)
mKCa_inf(V, Ca) = (Ca / (Ca+3)) * (1 / (1+exp((V+28.3)/-12.6)))
tau_mKCa(V) = tauX(V, 90.3, 75.1, 46., -22.7)

# Kd-current (m = activation variable)
mKd_inf(V) = boltz(V, 12.3, -11.8)
tau_mKd(V) = tauX(V, 7.2, 6.4, 28.3, -19.2)

# H-current (m = activation variable)
mH_inf(V) = boltz(V, 70., 6.)
tau_mH(V) = tauX(V, 272., -1499., 42.2, -8.73)
