#=
This file contains differential equations describing the STG model
=#

include("FYON_2022_STG_kinetics.jl") # Include STG model gating functions

## STG model from Liu 1998 - current-clamp mode
function STG_ODE(du, u, p, t)
    # Parameters
    Iapp  = p[1]  # Amplitude of constant applied current
    gNa   = p[2]  # Sodium current maximal conductance
    gCaT  = p[3]  # T-type calcium current maximal conductance
    gCaS  = p[4]  # Slow calcium current maximal conductance
    gA    = p[5]  # A-type potassium current maximal conductance
    gKCa  = p[6]  # Calcium controlled potassium current maximal conductance
    gKd   = p[7]  # Delayed-rectifier potassium current maximal conductance
    gH    = p[8]  # H current maximal conductance
    gleak = p[9]  # Leak current maximal conductance
    C     = p[10] # Membrane capacitance

    # Variables
    V    = u[1]  # Membrane potential
    mNa  = u[2]  # Sodium current activation
    hNa  = u[3]  # Sodium current inactivation
    mCaT = u[4]  # T-type calcium current activation
    hCaT = u[5]  # T-type calcium current inactivation
    mCaS = u[6]  # Slow calcium current activation
    hCaS = u[7]  # Slow calcium current inactivation
    mA   = u[8]  # A-type potassium current activation
    hA   = u[9]  # A-type potassium current inactivation
    mKCa = u[10] # Calcium controlled potassium current activation
    mKd  = u[11] # Delayed-rectifier potassium current activation
    mH   = u[12] # H current activation
    Ca   = u[13] # Calcium concentration

    # ODEs
    du[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
                   gKd*mKd^4*(V-VK) - gH*mH*(V-VH) - gleak*(V-Vleak) + Iapp)

    du[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
    du[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
    du[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
    du[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
    du[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
    du[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
    du[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
    du[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
    du[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
    du[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
    du[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)

    du[13] = (-0.94*(gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
end

function STG_ODE_HS(du, u, p, t)
    # Parameters
    Iapp  = p[1](t)  # Amplitude of constant applied current
    gNa   = p[2]     # Sodium current maximal conductance
    gCaT  = p[3]     # T-type calcium current maximal conductance
    gCaS  = p[4]     # Slow calcium current maximal conductance
    gA    = p[5]     # A-type potassium current maximal conductance
    gKCa  = p[6]     # Calcium controlled potassium current maximal conductance
    gKd   = p[7]     # Delayed-rectifier potassium current maximal conductance
    gH    = p[8]     # H current maximal conductance
    gleak = p[9]     # Leak current maximal conductance
    C     = p[10]    # Membrane capacitance

    # Variables
    V    = u[1]  # Membrane potential
    mNa  = u[2]  # Sodium current activation
    hNa  = u[3]  # Sodium current inactivation
    mCaT = u[4]  # T-type calcium current activation
    hCaT = u[5]  # T-type calcium current inactivation
    mCaS = u[6]  # Slow calcium current activation
    hCaS = u[7]  # Slow calcium current inactivation
    mA   = u[8]  # A-type potassium current activation
    hA   = u[9]  # A-type potassium current inactivation
    mKCa = u[10] # Calcium controlled potassium current activation
    mKd  = u[11] # Delayed-rectifier potassium current activation
    mH   = u[12] # H current activation
    Ca   = u[13] # Calcium concentration

    # ODEs
    du[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
                   gKd*mKd^4*(V-VK) - gH*mH*(V-VH) - gleak*(V-Vleak) + Iapp)

    du[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
    du[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
    du[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
    du[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
    du[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
    du[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
    du[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
    du[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
    du[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
    du[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
    du[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)

    du[13] = (-0.94*(gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
end

function STG_ODE_CR(du, u, p, t)
    # Parameters
    Iapp  = p[1]     # Amplitude of constant applied current
    gNa   = p[2]     # Sodium current maximal conductance
    gCaT  = p[3]     # T-type calcium current maximal conductance
    gCaS  = p[4]     # Slow calcium current maximal conductance
    gA    = p[5](t)  # A-type potassium current maximal conductance
    gKCa  = p[6]     # Calcium controlled potassium current maximal conductance
    gKd   = p[7]     # Delayed-rectifier potassium current maximal conductance
    gH    = p[8]     # H current maximal conductance
    gleak = p[9]     # Leak current maximal conductance
    C     = p[10]    # Membrane capacitance

    # Variables
    V    = u[1]  # Membrane potential
    mNa  = u[2]  # Sodium current activation
    hNa  = u[3]  # Sodium current inactivation
    mCaT = u[4]  # T-type calcium current activation
    hCaT = u[5]  # T-type calcium current inactivation
    mCaS = u[6]  # Slow calcium current activation
    hCaS = u[7]  # Slow calcium current inactivation
    mA   = u[8]  # A-type potassium current activation
    hA   = u[9]  # A-type potassium current inactivation
    mKCa = u[10] # Calcium controlled potassium current activation
    mKd  = u[11] # Delayed-rectifier potassium current activation
    mH   = u[12] # H current activation
    Ca   = u[13] # Calcium concentration

    # ODEs
    du[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
                   gKd*mKd^4*(V-VK) - gH*mH*(V-VH) - gleak*(V-Vleak) + Iapp)

    du[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
    du[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
    du[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
    du[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
    du[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
    du[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
    du[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
    du[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
    du[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
    du[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
    du[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)

    du[13] = (-0.94*(gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
end

## STG model from Liu 1998 - current-clamp mode with hyperpolarizing step
function STG_hyper_ODE(du, u, p, t)
    # Parameters
    Iapp  = p[1]  # Amplitude of constant applied current
    gNa   = p[2]  # Sodium current maximal conductance
    gCaT  = p[3]  # T-type calcium current maximal conductance
    gCaS  = p[4]  # Slow calcium current maximal conductance
    gA    = p[5]  # A-type potassium current maximal conductance
    gKCa  = p[6]  # Calcium controlled potassium current maximal conductance
    gKd   = p[7]  # Delayed-rectifier potassium current maximal conductance
    gH    = p[8]  # H current maximal conductance
    gleak = p[9]  # Leak current maximal conductance
    C     = p[10] # Membrane capacitance
    I1    = p[11] # Amplitude of the first step input
    I2    = p[12] # Amplitude of the second step input
    ti1   = p[13] # Starting time of the first step input
    tf1   = p[14] # Ending time of the first step input
    ti2   = p[15] # Starting time of the second step input
    tf2   = p[16] # Ending time of the second step input

    # Variables
    V    = u[1]  # Membrane potential
    mNa  = u[2]  # Sodium current activation
    hNa  = u[3]  # Sodium current inactivation
    mCaT = u[4]  # T-type calcium current activation
    hCaT = u[5]  # T-type calcium current inactivation
    mCaS = u[6]  # Slow calcium current activation
    hCaS = u[7]  # Slow calcium current inactivation
    mA   = u[8]  # A-type potassium current activation
    hA   = u[9]  # A-type potassium current inactivation
    mKCa = u[10] # Calcium controlled potassium current activation
    mKd  = u[11] # Delayed-rectifier potassium current activation
    mH   = u[12] # H current activation
    Ca   = u[13] # Calcium concentration

    # ODEs
    du[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
                   gKd*mKd^4*(V-VK) - gH*mH*(V-VH) - gleak*(V-Vleak) + Iapp +
                   I1*pulse(t, ti1, tf1) + I2*pulse(t, ti2, tf2))

    du[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
    du[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
    du[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
    du[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
    du[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
    du[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
    du[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
    du[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
    du[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
    du[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
    du[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)

    du[13] = (-0.94*(gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
end

## STG model from Liu 1998 - current-clamp mode with ion channel blockade
function STG_block_ODE(du, u, p, t)
    # Parameters
    Iapp  = p[1]  # Amplitude of constant applied current
    gNa   = p[2]  # Sodium current maximal conductance
    gCaT  = p[3]  # T-type calcium current maximal conductance
    gCaS  = p[4]  # Slow calcium current maximal conductance
    gA    = p[5]  # A-type potassium current maximal conductance
    gKCa  = p[6]  # Calcium controlled potassium current maximal conductance
    gKd   = p[7]  # Delayed-rectifier potassium current maximal conductance
    gH    = p[8]  # H current maximal conductance
    gleak = p[9]  # Leak current maximal conductance
    C     = p[10] # Membrane capacitance
    ti    = p[11] # Starting time of the ion channel blockade
    tf    = p[12] # Ending time of the ion channel blockade

    # Ion channel blockade
    if t > ti && t < tf
        gCaS = 0
    end

    # Variables
    V    = u[1]  # Membrane potential
    mNa  = u[2]  # Sodium current activation
    hNa  = u[3]  # Sodium current inactivation
    mCaT = u[4]  # T-type calcium current activation
    hCaT = u[5]  # T-type calcium current inactivation
    mCaS = u[6]  # Slow calcium current activation
    hCaS = u[7]  # Slow calcium current inactivation
    mA   = u[8]  # A-type potassium current activation
    hA   = u[9]  # A-type potassium current inactivation
    mKCa = u[10] # Calcium controlled potassium current activation
    mKd  = u[11] # Delayed-rectifier potassium current activation
    mH   = u[12] # H current activation
    Ca   = u[13] # Calcium concentration

    # ODEs
    du[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
                   gKd*mKd^4*(V-VK) - gH*mH*(V-VH) - gleak*(V-Vleak) + Iapp)

    du[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
    du[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
    du[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
    du[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
    du[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
    du[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
    du[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
    du[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
    du[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
    du[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
    du[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)

    du[13] = (-0.94*(gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
end

## STG model from Liu 1998 - current-clamp mode along neuromodulation path
function STG_ODE_neuromod(du, u, p, t)
    # Parameters
    Iapp   = p[1]  # Amplitude of constant applied current
    gNa    = p[2]  # Sodium current maximal conductance
    gCaT   = p[3]  # T-type calcium current maximal conductance
    gA_min = p[4]  # A-type potassium current maximal conductance minimum value
    gA_max = p[5]  # A-type potassium current maximal conductance maximum value
    gKCa   = p[6]  # Calcium controlled potassium current maximal conductance
    gKd    = p[7]  # Delayed-rectifier potassium current maximal conductance
    gH     = p[8]  # H current maximal conductance
    gleak  = p[9]  # Leak current maximal conductance
    C      = p[10] # Membrane capacitance
    a0     = p[11] # First coefficient of the regression line
    a1     = p[12] # Second coefficient of the regression line

    # Computing gA
    gA = gA_max - (gA_max - gA_min) * t/Tfinal

    # Computing gCaS
    gCaS = a0 + a1 * gA

    # Variables
    V    = u[1]  # Membrane potential
    mNa  = u[2]  # Sodium current activation
    hNa  = u[3]  # Sodium current inactivation
    mCaT = u[4]  # T-type calcium current activation
    hCaT = u[5]  # T-type calcium current inactivation
    mCaS = u[6]  # Slow calcium current activation
    hCaS = u[7]  # Slow calcium current inactivation
    mA   = u[8]  # A-type potassium current activation
    hA   = u[9]  # A-type potassium current inactivation
    mKCa = u[10] # Calcium controlled potassium current activation
    mKd  = u[11] # Delayed-rectifier potassium current activation
    mH   = u[12] # H current activation
    Ca   = u[13] # Calcium concentration

    # ODEs
    du[1] = 1/C*(- gNa*mNa^3*hNa*(V-VNa) - gCaT*mCaT^3*hCaT*(V-VCa) -
                   gCaS*mCaS^3*hCaS*(V-VCa) - gA*mA^3*hA*(V-VK) - gKCa*mKCa^4*(V-VK) -
                   gKd*mKd^4*(V-VK) - gH*mH*(V-VH) - gleak*(V-Vleak) + Iapp)

    du[2] = (1/tau_mNa(V)) * (mNa_inf(V) - mNa)
    du[3] = (1/tau_hNa(V)) * (hNa_inf(V) - hNa)
    du[4] = (1/tau_mCaT(V)) * (mCaT_inf(V) - mCaT)
    du[5] = (1/tau_hCaT(V)) * (hCaT_inf(V) - hCaT)
    du[6] = (1/tau_mCaS(V)) * (mCaS_inf(V) - mCaS)
    du[7] = (1/tau_hCaS(V)) * (hCaS_inf(V) - hCaS)
    du[8] = (1/tau_mA(V)) * (mA_inf(V) - mA)
    du[9] = (1/tau_hA(V)) * (hA_inf(V) - hA)
    du[10] = (1/tau_mKCa(V)) * (mKCa_inf(V, Ca) - mKCa)
    du[11] = (1/tau_mKd(V)) * (mKd_inf(V) - mKd)
    du[12] = (1/tau_mH(V)) * (mH_inf(V) - mH)

    du[13] = (-0.94*(gCaT*mCaT^3*hCaT*(V-VCa) + gCaS*mCaS^3*hCaS*(V-VCa)) - Ca + 0.05)/20
end

## Stimulation function
heaviside(t) = (1 + sign(t)) / 2
pulse(t, ti, tf) = heaviside(t-ti) - heaviside(t-tf)
