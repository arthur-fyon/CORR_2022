#=
This file contains all derivatives of STG model gating functions
=#

# Derivative of the gating function
dboltz(V, A, B) = -exp((A + V)/B) / (B * (exp((A + V)/B) + 1)^2) #boltz(V, A, B) = 1 / (1 + exp((V+A)/B))

# Na-current (m = activation variable, h = inactivation variable)
dmNa(V) = dboltz(V, 25.5, -5.29)
dhNa(V) = dboltz(V, 48.9, 5.18)

# T-type Ca-current (m = activation variable, h = inactivation variable)
dmCaT(V) = dboltz(V, 27.1, -7.2)
dhCaT(V) = dboltz(V, 32.1, 5.5)

# Slow Ca-current (m = activation variable, h = inactivation variable)
dmCaS(V) = dboltz(V, 33., -8.1)
dhCaS(V) = dboltz(V, 60., 6.2)

# A-current (m = activation variable, h = inactivation variable)
dmA(V) = dboltz(V, 27.2, -8.7)
dhA(V) = dboltz(V, 56.9, 4.9)

# Ca K-current (m = activation variable)
dmKCa(V, Ca) = (5*Ca) / (63*exp((5*V)/63 + 283/126) * (Ca + 3) * (1/exp((5*V)/63 + 283/126) + 1)^2) # dmKCa/dV
dmKCa2(V, Ca) = 1 / ((Ca + 3) * (1/exp((5*V)/63 + 283/126) + 1)) - Ca / ((Ca + 3)^2 * (1/exp((5*V)/63 + 283/126) + 1)) # dmKCa/dCa

# Kd-current (m = activation variable)
dmKd(V) = dboltz(V, 12.3, -11.8)

# H-current (m = activation variable)
dmH(V) = dboltz(V, 70., 6.)

# Ca dynamics (dCa_dot/dV)
dCa(V, Ca, gCaT, gCaS) = (47*gCaS*exp((5*V)/31 + 300/31)*(V - VCa))/(310*(exp(- (10*V)/81 - 110/27) + 1)^3*(exp((5*V)/31 + 300/31) + 1)^2) -
    (47*gCaT)/(50*(exp(- (5*V)/36 - 271/72) + 1)^3*(exp((2*V)/11 + 321/55) + 1)) -
    (47*gCaS*exp(- (10*V)/81 - 110/27)*(V - VCa))/(135*(exp(- (10*V)/81 - 110/27) + 1)^4*(exp((5*V)/31 + 300/31) + 1)) -
    (47*gCaS)/(50*(exp(- (10*V)/81 - 110/27) + 1)^3*(exp((5*V)/31 + 300/31) + 1)) -
    (47*gCaT*exp(- (5*V)/36 - 271/72)*(V - VCa))/(120*(exp(- (5*V)/36 - 271/72) + 1)^4*(exp((2*V)/11 + 321/55) + 1)) +
    (47*gCaT*exp((2*V)/11 + 321/55)*(V - VCa))/(275*(exp(- (5*V)/36 - 271/72) + 1)^3*(exp((2*V)/11 + 321/55) + 1)^2)
