#=
This file contains differential equations describing the STG model
=#

include("FYON_2022_DA_kinetics.jl") # Include DA model gating functions

## DA model from Drion 2011 - current-clamp mode
function DA_ODE(du, u, p, t)
    # Parameters
    Iapp  = p[1] # Amplitude of constant applied current
    gNa   = p[2] # Sodium current maximal conductance
    gKd   = p[3] # Delayed-rectifier potassium current maximal conductance
    gCaL  = p[4] # L-type calcium current maximal conductance
    gCaN  = p[5] # N-type calcium current maximal conductance
    gERG  = p[6] # ERG potassium current maximal conductance
    gNMDA = p[7] # NMDA current maximal conductance
    gleak = p[8] # Leak current maximal conductance
    C     = p[9] # Membrane capacitance

    # Variables
    V    = u[1] # Membrane potential
    m    = u[2] # Sodium current activation
    h    = u[3] # Sodium current inactivation
    n    = u[4] # Delayed-rectifier potassium current activation
    mCaL = u[5] # L-type calcium current activation
    mCaN = u[6] # N-type calcium current activation
    oERG = u[7] # ERG potassium current activation
    iERG = u[8] # ERG current intermediate

    # ODEs
    du[1]= 1/C*(- gNa*m^3*h*(V-VNa) - gKd*n^3*(V-VK) - gCaL*mCaL^2*(V-VCa) -
                  gCaN*mCaN*(V-VCa) - gERG*oERG*(V-VK) - gleak*(V-Vleak) -
                  gNMDA*(V-VNMDA) / (1 + Mg*exp(-0.08*V)/10.) + Iapp)

    du[2] = (1/tau_m(V)) * (m_inf(V) - m)
    du[3] = (1/tau_h(V)) * (h_inf(V) - h)
    du[4] = (1/tau_n(V)) * (n_inf(V) - n)
    du[5] = (1/tau_mCaL(V)) * (mCaL_inf(V) - mCaL)
    du[6] = (1/tau_mCaN(V)) * (mCaN_inf(V) - mCaN)
    du[7] = a0ERG(V) * (1-oERG-iERG) + biERG(V) * iERG - oERG * (aiERG(V)+b0ERG(V))
    du[8] = aiERG(V) * oERG - biERG(V) * iERG
end

## DA model from Drion 2011 - current-clamp mode with hyperpolarizing step
function DA_hyper_ODE(du, u, p, t)
    # Parameters
    Iapp  = p[1]  # Amplitude of constant applied current
    gNa   = p[2]  # Sodium current maximal conductance
    gKd   = p[3]  # Delayed-rectifier potassium current maximal conductance
    gCaL  = p[4]  # L-type calcium current maximal conductance
    gCaN  = p[5]  # N-type calcium current maximal conductance
    gERG  = p[6]  # ERG potassium current maximal conductance
    gNMDA = p[7]  # NMDA current maximal conductance
    gleak = p[8]  # Leak current maximal conductance
    C     = p[9]  # Membrane capacitance
    I1    = p[10] # Amplitude of the first step input
    I2    = p[11] # Amplitude of the second step input
    ti1   = p[12] # Starting time of the first step input
    tf1   = p[13] # Ending time of the first step input
    ti2   = p[14] # Starting time of the second step input
    tf2   = p[15] # Ending time of the second step input

    # Variables
    V    = u[1] # Membrane potential
    m    = u[2] # Sodium current activation
    h    = u[3] # Sodium current inactivation
    n    = u[4] # Delayed-rectifier potassium current activation
    mCaL = u[5] # L-type calcium current activation
    mCaN = u[6] # N-type calcium current activation
    oERG = u[7] # ERG potassium current activation
    iERG = u[8] # ERG current intermediate

    # ODEs
    du[1]= 1/C*(- gNa*m^3*h*(V-VNa) - gKd*n^3*(V-VK) - gCaL*mCaL^2*(V-VCa) -
                  gCaN*mCaN*(V-VCa) - gERG*oERG*(V-VK) - gleak*(V-Vleak) -
                  gNMDA*(V-VNMDA) / (1 + Mg*exp(-0.08*V)/10.) + Iapp +
                  I1*pulse(t, ti1, tf1) + I2*pulse(t, ti2, tf2))

    du[2] = (1/tau_m(V)) * (m_inf(V) - m)
    du[3] = (1/tau_h(V)) * (h_inf(V) - h)
    du[4] = (1/tau_n(V)) * (n_inf(V) - n)
    du[5] = (1/tau_mCaL(V)) * (mCaL_inf(V) - mCaL)
    du[6] = (1/tau_mCaN(V)) * (mCaN_inf(V) - mCaN)
    du[7] = a0ERG(V) * (1-oERG-iERG) + biERG(V) * iERG - oERG * (aiERG(V)+b0ERG(V))
    du[8] = aiERG(V) * oERG - biERG(V) * iERG
end

## DA model from Drion 2011 - current-clamp mode with ion channel blockade
function DA_block_ODE(du, u, p, t)
    # Parameters
    Iapp  = p[1]  # Amplitude of constant applied current
    gNa   = p[2]  # Sodium current maximal conductance
    gKd   = p[3]  # Delayed-rectifier potassium current maximal conductance
    gCaL  = p[4]  # L-type calcium current maximal conductance
    gCaN  = p[5]  # N-type calcium current maximal conductance
    gERG  = p[6]  # ERG potassium current maximal conductance
    gNMDA = p[7]  # NMDA current maximal conductance
    gleak = p[8]  # Leak current maximal conductance
    C     = p[9]  # Membrane capacitance
    ti    = p[10] # Starting time of the ion channel blockade
    tf    = p[11] # Ending time of the ion channel blockade

    # Ion channel blockade
    if t > ti && t < tf
        gCaL = 0
    end

    # Variables
    V    = u[1] # Membrane potential
    m    = u[2] # Sodium current activation
    h    = u[3] # Sodium current inactivation
    n    = u[4] # Delayed-rectifier potassium current activation
    mCaL = u[5] # L-type calcium current activation
    mCaN = u[6] # N-type calcium current activation
    oERG = u[7] # ERG potassium current activation
    iERG = u[8] # ERG current intermediate

    # ODEs
    du[1]= 1/C*(- gNa*m^3*h*(V-VNa) - gKd*n^3*(V-VK) - gCaL*mCaL^2*(V-VCa) -
                  gCaN*mCaN*(V-VCa) - gERG*oERG*(V-VK) - gleak*(V-Vleak) -
                  gNMDA*(V-VNMDA) / (1 + Mg*exp(-0.08*V)/10.) + Iapp)

    du[2] = (1/tau_m(V)) * (m_inf(V) - m)
    du[3] = (1/tau_h(V)) * (h_inf(V) - h)
    du[4] = (1/tau_n(V)) * (n_inf(V) - n)
    du[5] = (1/tau_mCaL(V)) * (mCaL_inf(V) - mCaL)
    du[6] = (1/tau_mCaN(V)) * (mCaN_inf(V) - mCaN)
    du[7] = a0ERG(V) * (1-oERG-iERG) + biERG(V) * iERG - oERG * (aiERG(V)+b0ERG(V))
    du[8] = aiERG(V) * oERG - biERG(V) * iERG
end

## DA model from Drion 2011 - current-clamp mode along neuromodulation path
function DA_ODE_neuromod(du, u, p, t)
    # Parameters
    Iapp     = p[1]  # Amplitude of constant applied current
    gNa      = p[2]  # Sodium current maximal conductance
    gKd      = p[3]  # Delayed-rectifier potassium current maximal conductance
    gCaL_min = p[4]  # L-type calcium current maximal conductance minimal value
    gCaL_max = p[5]  # L-type calcium current maximal conductance maximal value
    gERG     = p[6]  # ERG potassium current maximal conductance
    gNMDA    = p[7]  # NMDA current maximal conductance
    gleak    = p[8]  # Leak current maximal conductance
    C        = p[9]  # Membrane capacitance
    a0       = p[10] # First coefficient of the regression line
    a1       = p[11] # Second coefficient of the regression line

    # Computing gCaL
    gCaL = gCaL_max - (gCaL_max - gCaL_min) * t/Tfinal

    # Computing gCaN
    gCaN = a0 + a1 * gCaL

    # Variables
    V    = u[1] # Membrane potential
    m    = u[2] # Sodium current activation
    h    = u[3] # Sodium current inactivation
    n    = u[4] # Delayed-rectifier potassium current activation
    mCaL = u[5] # L-type calcium current activation
    mCaN = u[6] # N-type calcium current activation
    oERG = u[7] # ERG potassium current activation
    iERG = u[8] # ERG current intermediate

    # ODEs
    du[1]= 1/C*(- gNa*m^3*h*(V-VNa) - gKd*n^3*(V-VK) - gCaL*mCaL^2*(V-VCa) -
                  gCaN*mCaN*(V-VCa) - gERG*oERG*(V-VK) - gleak*(V-Vleak) -
                  gNMDA*(V-VNMDA) / (1 + Mg*exp(-0.08*V)/10.) + Iapp)

    du[2] = (1/tau_m(V)) * (m_inf(V) - m)
    du[3] = (1/tau_h(V)) * (h_inf(V) - h)
    du[4] = (1/tau_n(V)) * (n_inf(V) - n)
    du[5] = (1/tau_mCaL(V)) * (mCaL_inf(V) - mCaL)
    du[6] = (1/tau_mCaN(V)) * (mCaN_inf(V) - mCaN)
    du[7] = a0ERG(V) * (1-oERG-iERG) + biERG(V) * iERG - oERG * (aiERG(V)+b0ERG(V))
    du[8] = aiERG(V) * oERG - biERG(V) * iERG
end

## Stimulation function
heaviside(t) = (1 + sign(t)) / 2
pulse(t, ti, tf) = heaviside(t-ti) - heaviside(t-tf)
