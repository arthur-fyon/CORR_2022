#=
This file computes the dynamic input conductances of the STG model as well as
applying the compensation algorithm for a set of maximal conductances (inputs of the algorithm)
=#

include("FYON_2022_STG_kinetics.jl") # Include model gating functions
include("FYON_2022_STG_gs_derivatives.jl") # Include derivatives of model gating functions (dX/dV)

## DIC computation functions
# These 2 functions compute variable contributions in the fast, slow and ultraslow timescales (wfs(V) and wsu(V))
function dist(tauX, tauA, tauB)
   (log(tauB) - log(tauX))/(log(tauB) - log(tauA))
end

function var_contribution(tauX, tauf, taus, tauu)
  wfs = 1.
  wsu = 1.
  if tauX < tauf
   wfs = 1.
   wsu = 1.
  elseif tauf <= tauX < taus
   wfs = dist(tauX, tauf, taus)
   wsu = 1.
  elseif taus <= tauX < tauu
   wfs = 0.
   wsu = dist(tauX, taus, tauu)
  else
   wfs = 0.
   wsu = 0.
  end
  return wfs, wsu
end

# This function computes the dynamic input conductances and Istatic
function DICs(V, Iapp, gNa, gCaT, gCaS, gA, gKCa, gKd, gH, gleak)
  # Initializing some variables, indices and flag
  gf = zeros(length(V))
  gs = zeros(length(V))
  gu = zeros(length(V))
  gin = zeros(length(V))
  Istatic = zeros(length(V))
  ith = 1
  iosc = 1
  flag_th = 0

  # loop over all acceptable values of the membrane voltage
  for i = 2 : length(V)
    # Defining the 3 time constants of the 3 timescales
    tau_fast = tau_mNa(V[i])
    tau_slow = tau_mKd(V[i])
    tau_uslow = tau_mH(V[i])

    # Computing the wfs(V) and wsu(V) for all gating variables
    (wfs_mNa, wsu_mNa) = var_contribution(tau_mNa(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_hNa, wsu_hNa) = var_contribution(tau_hNa(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_mCaT, wsu_mCaT) = var_contribution(tau_mCaT(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_hCaT, wsu_hCaT) = var_contribution(tau_hCaT(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_mCaS, wsu_mCaS) = var_contribution(tau_mCaS(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_hCaS, wsu_hCaS) = var_contribution(tau_hCaS(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_mA, wsu_mA) = var_contribution(tau_mA(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_hA, wsu_hA) = var_contribution(tau_hA(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_mKd, wsu_mKd) = var_contribution(tau_mKd(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_mKCa, wsu_mKCa) = var_contribution(tau_mKCa(V[i]), tau_fast, tau_slow, tau_uslow)
    (wfs_mH, wsu_mH) = var_contribution(tau_mH(V[i]), tau_fast, tau_slow, tau_uslow)

    # Calcium dynamic in the ultraslow timescale
    wfs_mKCa2 = 0.
    wsu_mKCa2 = 0.

    # Computing equilibrium value of [Ca]
    Ca = -0.94*(gCaT*mCaT_inf(V[i])^3*hCaT_inf(V[i])*(V[i]-VCa) + gCaS*mCaS_inf(V[i])^3*hCaS_inf(V[i])*(V[i]-VCa)) + 0.05

    # Computing dV_dot/dV(V) appearing in the fast timescale (Ohm's law)
    dvdot_dv = - gNa*mNa_inf(V[i])^3*hNa_inf(V[i]) - gCaT*mCaT_inf(V[i])^3*hCaT_inf(V[i]) - gCaS*mCaS_inf(V[i])^3*hCaS_inf(V[i]) -
                 gA*mA_inf(V[i])^3*hA_inf(V[i]) - gKCa*mKCa_inf(V[i], Ca)^4 - gKd*mKd_inf(V[i])^4 - gH*mH_inf(V[i]) - gleak

    # Computing fast input conductance (g_fast(V))
    gf[i] = - dvdot_dv + wfs_mNa*3*gNa*mNa_inf(V[i])^2*hNa_inf(V[i])*(V[i]-VNa)*dmNa(V[i]) + wfs_hNa*gNa*mNa_inf(V[i])^3*(V[i]-VNa)*dhNa(V[i]) +
              wfs_mCaT*3*gCaT*mCaT_inf(V[i])^2*hCaT_inf(V[i])*(V[i]-VCa)*dmCaT(V[i]) + wfs_hCaT*gCaT*mCaT_inf(V[i])^3*(V[i]-VCa)*dhCaT(V[i]) +
              wfs_mCaS*3*gCaS*mCaS_inf(V[i])^2*hCaS_inf(V[i])*(V[i]-VCa)*dmCaS(V[i]) + wfs_hCaS*gCaS*mCaS_inf(V[i])^3*(V[i]-VCa)*dhCaS(V[i]) +
              wfs_mA*3*gA*mA_inf(V[i])^2*hA_inf(V[i])*(V[i]-VK)*dmA(V[i]) + wfs_hA*gA*mA_inf(V[i])^3*(V[i]-VK)*dhA(V[i]) +
              wfs_mKd*gKd*4*mKd_inf(V[i])^3*(V[i]-VK)*dmKd(V[i]) + wfs_mKCa*gKCa*4*mKCa_inf(V[i], Ca)^3*(V[i]-VK)*dmKCa(V[i], Ca) +
              wfs_mKCa2*gKCa*4*mKCa_inf(V[i], Ca)^3*(V[i]-VK).*dmKCa2(V[i], Ca).*dCa(V[i], Ca, gCaT, gCaS) + wfs_mH*gH*(V[i]-VH)*dmH(V[i])

    # Computing slow input conductance (g_slow(V))
    gs[i] = (wsu_mNa - wfs_mNa)*3*gNa*mNa_inf(V[i])^2*hNa_inf(V[i])*(V[i]-VNa)*dmNa(V[i]) + (wsu_hNa - wfs_hNa)*gNa*mNa_inf(V[i])^3*(V[i]-VNa)*dhNa(V[i]) +
            (wsu_mCaT - wfs_mCaT)*3*gCaT*mCaT_inf(V[i])^2*hCaT_inf(V[i])*(V[i]-VCa)*dmCaT(V[i]) + (wsu_hCaT - wfs_hCaT)*gCaT*mCaT_inf(V[i])^3*(V[i]-VCa)*dhCaT(V[i]) +
            (wsu_mCaS - wfs_mCaS)*3*gCaS*mCaS_inf(V[i])^2*hCaS_inf(V[i])*(V[i]-VCa)*dmCaS(V[i]) + (wsu_hCaS - wfs_hCaS)*gCaS*mCaS_inf(V[i])^3*(V[i]-VCa)*dhCaS(V[i]) +
            (wsu_mA - wfs_mA)*3*gA*mA_inf(V[i])^2*hA_inf(V[i])*(V[i]-VK)*dmA(V[i]) + (wsu_hA - wfs_hA)*gA*mA_inf(V[i])^3*(V[i]-VK)*dhA(V[i]) +
            (wsu_mKd - wfs_mKd)*gKd*4*mKd_inf(V[i])^3*(V[i]-VK)*dmKd(V[i]) + (wsu_mKCa - wfs_mKCa)*gKCa*4*mKCa_inf(V[i], Ca)^3*(V[i]-VK)*dmKCa(V[i], Ca) +
            (wsu_mKCa2 - wfs_mKCa2)*gKCa*4*mKCa_inf(V[i], Ca)^3*(V[i]-VK).*dmKCa2(V[i], Ca).*dCa(V[i], Ca, gCaT, gCaS) + (wsu_mH - wfs_mH)*gH*(V[i]-VH)*dmH(V[i])

    # Computing ultraslow input conductance (g_uslow(V))
    gu[i] = (1 - wsu_mNa)*3*gNa*mNa_inf(V[i])^2*hNa_inf(V[i])*(V[i]-VNa)*dmNa(V[i]) + (1 - wsu_hNa)*gNa*mNa_inf(V[i])^3*(V[i]-VNa)*dhNa(V[i]) +
            (1 - wsu_mCaT)*3*gCaT*mCaT_inf(V[i])^2*hCaT_inf(V[i])*(V[i]-VCa)*dmCaT(V[i]) + (1 - wsu_hCaT)*gCaT*mCaT_inf(V[i])^3*(V[i]-VCa)*dhCaT(V[i]) +
            (1 - wsu_mCaS)*3*gCaS*mCaS_inf(V[i])^2*hCaS_inf(V[i])*(V[i]-VCa)*dmCaS(V[i]) + (1 - wsu_hCaS)*gCaS*mCaS_inf(V[i])^3*(V[i]-VCa)*dhCaS(V[i]) +
            (1 - wsu_mA)*3*gA*mA_inf(V[i])^2*hA_inf(V[i])*(V[i]-VK)*dmA(V[i]) + (1 - wsu_hA)*gA*mA_inf(V[i])^3*(V[i]-VK)*dhA(V[i]) +
            (1 - wsu_mKd)*gKd*4*mKd_inf(V[i])^3*(V[i]-VK)*dmKd(V[i]) + (1 - wsu_mKCa)*gKCa*4*mKCa_inf(V[i], Ca)^3*(V[i]-VK)*dmKCa(V[i], Ca) +
            (1 - wsu_mKCa2)*gKCa*4*mKCa_inf(V[i], Ca)^3*(V[i]-VK).*dmKCa2(V[i], Ca).*dCa(V[i], Ca, gCaT, gCaS) + (1 - wsu_mH)*gH*(V[i]-VH)*dmH(V[i])

    # Detecting the threshold voltage
    if flag_th == 0 && gf[i-1] + gs[i-1] + gu[i-1] > 0 && gf[i] + gs[i] + gu[i] <= 0
      ith = i
      flag_th = 1
    end

    # Detecting the up state voltage
    if Istatic[i-1] > 0 && Istatic[i] <= 0
      iosc = i
    end

    # Computing Istatic
    Istatic[i] = - gNa*mNa_inf(V[i])^3*hNa_inf(V[i])*(V[i]-VNa) - gCaT*mCaT_inf(V[i])^3*hCaT_inf(V[i])*(V[i]-VCa) -
                   gCaS*mCaS_inf(V[i])^3*hCaS_inf(V[i])*(V[i]-VCa) - gA*mA_inf(V[i])^3*hA_inf(V[i])*(V[i]-VK) -
                   gKCa*mKCa_inf(V[i],Ca)^4*(V[i]-VK) - gKd*mKd_inf(V[i])^4*(V[i]-VK) -
                   gH*mH_inf(V[i])*(V[i]-VH) - gleak*(V[i]-Vleak) + Iapp
  end

  # Normalizing everything by g_leakage
  gf = gf ./ gleak
  gs = gs ./ gleak
  gu = gu ./ gleak
  gin = gf + gs + gu
  Istatic = Istatic ./ gleak

  return ith, iosc, gf, gs, gu, gin, Istatic
end

## Compensation algorithm functions
# This function computes the values of (gNa, gA, gH) that generate the values of the DICs
# at threshold voltage for the given set of conductances (gCaT, gCaS, gKd, gKCa, gleak)
function DICs_gmax_initNaAH(gCaT, gCaS, gKd, gKCa, gleak, gf_th, gs_th, gu_th, Vth)
  # Defining the 3 time constants of the 3 timescales
  tau_fast = tau_mNa(Vth)
  tau_slow = tau_mKd(Vth)
  tau_uslow = tau_mH(Vth)

  # Computing the wfs(V) and wsu(V) for all gating variables
  (wfs_mNa, wsu_mNa) = var_contribution(tau_mNa(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_hNa, wsu_hNa) = var_contribution(tau_hNa(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mCaT, wsu_mCaT) = var_contribution(tau_mCaT(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_hCaT, wsu_hCaT) = var_contribution(tau_hCaT(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mCaS, wsu_mCaS) = var_contribution(tau_mCaS(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_hCaS, wsu_hCaS) = var_contribution(tau_hCaS(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mA, wsu_mA) = var_contribution(tau_mA(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_hA, wsu_hA) = var_contribution(tau_hA(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mKd, wsu_mKd) = var_contribution(tau_mKd(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mKCa, wsu_mKCa) = var_contribution(tau_mKCa(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mH, wsu_mH) = var_contribution(tau_mH(Vth), tau_fast, tau_slow, tau_uslow)

  # Calcium dynamic in the ultraslow timescale
  wfs_mKCa2 = 0.
  wsu_mKCa2 = 0.

  # Computing equilibrium value of [Ca]
  Ca_th = -0.94*(gCaT*mCaT_inf(Vth)^3*hCaT_inf(Vth)*(Vth-VCa) + gCaS*mCaS_inf(Vth)^3*hCaS_inf(Vth)*(Vth-VCa)) + 0.05

  # Computing dV_dot/dV(V) appearing in the fast timescale (Ohm's law)
  dvdot_dv =  - gCaT*mCaT_inf(Vth)^3*hCaT_inf(Vth) - gCaS*mCaS_inf(Vth)^3*hCaS_inf(Vth) -
                gKCa*mKCa_inf(Vth, Ca_th)^4 - gKd*mKd_inf(Vth)^4 - gleak

  # Initializing the linear system of the compensation algorithm
  A = zeros(3, 3)
  B = zeros(3, 1)

  # Filling the A matrix of the compensation algorithm (= contribution of the unknowns to the DICs)
  A[1, 1] = (1/gleak) * (mNa_inf(Vth)^3*hNa_inf(Vth) + wfs_mNa*3*mNa_inf(Vth)^2*hNa_inf(Vth)*(Vth-VNa)*dmNa(Vth) + wfs_hNa*mNa_inf(Vth)^3*(Vth-VNa)*dhNa(Vth))
  A[1, 2] = (1/gleak) * (mA_inf(Vth)^3*hA_inf(Vth) + wfs_mA*3*mA_inf(Vth)^2*hA_inf(Vth)*(Vth-VK)*dmA(Vth) + wfs_hA*mA_inf(Vth)^3*(Vth-VK)*dhA(Vth))
  A[1, 3] = (1/gleak) * (mH_inf(Vth) + wfs_mH*(Vth-VH)*dmH(Vth))
  A[2, 1] = (1/gleak) * ((wsu_mNa - wfs_mNa)*3*mNa_inf(Vth)^2*hNa_inf(Vth)*(Vth-VNa)*dmNa(Vth) + (wsu_hNa - wfs_hNa)*mNa_inf(Vth)^3*(Vth-VNa)*dhNa(Vth))
  A[2, 2] = (1/gleak) * ((wsu_mA - wfs_mA)*3*mA_inf(Vth)^2*hA_inf(Vth)*(Vth-VK)*dmA(Vth) + (wsu_hA - wfs_hA)*mA_inf(Vth)^3*(Vth-VK)*dhA(Vth))
  A[2, 3] = (1/gleak) * ((wsu_mH - wfs_mH)*(Vth-VH)*dmH(Vth))
  A[3, 1] = (1/gleak) * ((1 - wsu_mNa)*3*mNa_inf(Vth)^2*hNa_inf(Vth)*(Vth-VNa)*dmNa(Vth) + (1 - wsu_hNa)*mNa_inf(Vth)^3*(Vth-VNa)*dhNa(Vth))
  A[3, 2] = (1/gleak) * ((1 - wsu_mA)*3*mA_inf(Vth)^2*hA_inf(Vth)*(Vth-VK)*dmA(Vth) + (1 - wsu_hA)*mA_inf(Vth)^3*(Vth-VK)*dhA(Vth))
  A[3, 3] = (1/gleak) * ((1 - wsu_mH)*(Vth-VH)*dmH(Vth))

  # Filling the B matrix of the compensation algorithm (= contribution of the non-unknowns to the DICs)
  B[1, 1] = gf_th - (1/gleak) * (- dvdot_dv + wfs_mCaT*3*gCaT*mCaT_inf(Vth)^2*hCaT_inf(Vth)*(Vth-VCa)*dmCaT(Vth) + wfs_hCaT*gCaT*mCaT_inf(Vth)^3*(Vth-VCa)*dhCaT(Vth) +
                                              wfs_mCaS*3*gCaS*mCaS_inf(Vth)^2*hCaS_inf(Vth)*(Vth-VCa)*dmCaS(Vth) + wfs_hCaS*gCaS*mCaS_inf(Vth)^3*(Vth-VCa)*dhCaS(Vth) +
                                              wfs_mKd*gKd*4*mKd_inf(Vth)^3*(Vth-VK)*dmKd(Vth) + wfs_mKCa*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK)*dmKCa(Vth,Ca_th) +
                                              wfs_mKCa2*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK).*dmKCa2(Vth, Ca_th).*dCa(Vth, Ca_th, gCaT, gCaS))
  B[2, 1] = gs_th - (1/gleak) * ((wsu_mCaT - wfs_mCaT)*3*gCaT*mCaT_inf(Vth)^2*hCaT_inf(Vth)*(Vth-VCa)*dmCaT(Vth) + (wsu_hCaT - wfs_hCaT)*gCaT*mCaT_inf(Vth)^3*(Vth-VCa)*dhCaT(Vth) +
                                 (wsu_mCaS - wfs_mCaS)*3*gCaS*mCaS_inf(Vth)^2*hCaS_inf(Vth)*(Vth-VCa)*dmCaS(Vth) + (wsu_hCaS - wfs_hCaS)*gCaS*mCaS_inf(Vth)^3*(Vth-VCa)*dhCaS(Vth) +
                                 (wsu_mKd - wfs_mKd)*gKd*4*mKd_inf(Vth)^3*(Vth-VK)*dmKd(Vth) + (wsu_mKCa - wfs_mKCa)*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK)*dmKCa(Vth, Ca_th) +
                                 (wsu_mKCa2 - wfs_mKCa2)*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK).*dmKCa2(Vth, Ca_th).*dCa(Vth, Ca_th, gCaT, gCaS))
  B[3, 1] = gu_th - (1/gleak) * ((1 - wsu_mCaT)*3*gCaT*mCaT_inf(Vth)^2*hCaT_inf(Vth)*(Vth-VCa)*dmCaT(Vth) + (1 - wsu_hCaT)*gCaT*mCaT_inf(Vth)^3*(Vth-VCa)*dhCaT(Vth) +
                                 (1 - wsu_mCaS)*3*gCaS*mCaS_inf(Vth)^2*hCaS_inf(Vth)*(Vth-VCa)*dmCaS(Vth) + (1 - wsu_hCaS)*gCaS*mCaS_inf(Vth)^3*(Vth-VCa)*dhCaS(Vth) +
                                 (1 - wsu_mKd)*gKd*4*mKd_inf(Vth)^3*(Vth-VK)*dmKd(Vth) + (1 - wsu_mKCa)*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK)*dmKCa(Vth, Ca_th) +
                                 (1 - wsu_mKCa2)*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK).*dmKCa2(Vth, Ca_th).*dCa(Vth, Ca_th, gCaT, gCaS))

  # Solving the linear system
  g_sol = \(A, B)
  return (g_sol[1], g_sol[2], g_sol[3])
end

# This function computes the values of (gCaS, gA) that generate the values of the DICs (slow and uslow)
# at threshold voltage for the given set of conductances (gNa, gCaT, gKd, gKCa, gH, gleak)
function DICs_gmax_neuromodCaSA(gNa, gCaT, gKd, gKCa, gH, gleak, gs_th, gu_th, Vth)
  # Defining the 3 time constants of the 3 timescales
  tau_fast = tau_mNa(Vth)
  tau_slow = tau_mKd(Vth)
  tau_uslow = tau_mH(Vth)

  # Computing the wfs(V) and wsu(V) for all gating variables
  (wfs_mNa, wsu_mNa) = var_contribution(tau_mNa(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_hNa, wsu_hNa) = var_contribution(tau_hNa(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mCaT, wsu_mCaT) = var_contribution(tau_mCaT(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_hCaT, wsu_hCaT) = var_contribution(tau_hCaT(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mCaS, wsu_mCaS) = var_contribution(tau_mCaS(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_hCaS, wsu_hCaS) = var_contribution(tau_hCaS(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mA, wsu_mA) = var_contribution(tau_mA(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_hA, wsu_hA) = var_contribution(tau_hA(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mKd, wsu_mKd) = var_contribution(tau_mKd(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mKCa, wsu_mKCa) = var_contribution(tau_mKCa(Vth), tau_fast, tau_slow, tau_uslow)
  (wfs_mH, wsu_mH) = var_contribution(tau_mH(Vth), tau_fast, tau_slow, tau_uslow)

  # Calcium dynamic in the ultraslow timescale
  wfs_mKCa2 = 0.
  wsu_mKCa2 = 0.

  # Computing equilibrium value of [Ca] (have to approximate gCaS)
  Ca_th = -0.94*(gCaT*mCaT_inf(Vth)^3*hCaT_inf(Vth)*(Vth-VCa) + 10*mCaS_inf(Vth)^3*hCaS_inf(Vth)*(Vth-VCa)) + 0.05

  # Initializing the linear system of the compensation algorithm
  A = zeros(2, 2)
  B = zeros(2, 1)

  # Filling the A matrix of the compensation algorithm (= contribution of the unknowns to the DICs)
  A[1, 1] = (1/gleak) * ((wsu_mCaS - wfs_mCaS)*3*mCaS_inf(Vth)^2*hCaS_inf(Vth)*(Vth-VCa)*dmCaS(Vth) + (wsu_hCaS - wfs_hCaS)*mCaS_inf(Vth)^3*(Vth-VCa)*dhCaS(Vth))
  A[1, 2] = (1/gleak) * ((wsu_mA - wfs_mA)*3*mA_inf(Vth)^2*hA_inf(Vth)*(Vth-VK)*dmA(Vth) + (wsu_hA - wfs_hA)*mA_inf(Vth)^3*(Vth-VK)*dhA(Vth))
  A[2, 1] = (1/gleak) * ((1 - wsu_mCaS)*3*mCaS_inf(Vth)^2*hCaS_inf(Vth)*(Vth-VCa)*dmCaS(Vth) + (1 - wsu_hCaS)*mCaS_inf(Vth)^3*(Vth-VCa)*dhCaS(Vth))
  A[2, 2] = (1/gleak) * ((1 - wsu_mA)*3*mA_inf(Vth)^2*hA_inf(Vth)*(Vth-VK)*dmA(Vth) + (1 - wsu_hA)*mA_inf(Vth)^3*(Vth-VK)*dhA(Vth))

  # Filling the B matrix of the compensation algorithm (= contribution of the non-unknowns to the DICs)
  B[1, 1] = gs_th - (1/gleak) * ((wsu_mNa - wfs_mNa)*gNa*3*mNa_inf(Vth)^2*hNa_inf(Vth)*(Vth-VNa)*dmNa(Vth) + (wsu_hNa - wfs_hNa)*gNa*mNa_inf(Vth)^3*(Vth-VNa)*dhNa(Vth) +
                                (wsu_mCaT - wfs_mCaT)*3*gCaT*mCaT_inf(Vth)^2*hCaT_inf(Vth)*(Vth-VCa)*dmCaT(Vth) + (wsu_hCaT- wfs_hCaT)*gCaT*mCaT_inf(Vth)^3*(Vth-VCa)*dhCaT(Vth) +
                                (wsu_mKd - wfs_mKd)*gKd*4*mKd_inf(Vth)^3*(Vth-VK)*dmKd(Vth) + (wsu_mKCa - wfs_mKCa)*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK)*dmKCa(Vth, Ca_th) +
                                (wsu_mKCa2 - wfs_mKCa2)*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK).*dmKCa2(Vth, Ca_th).*dCa(Vth, Ca_th, gCaT, 10.) + (wsu_mH - wfs_mH)*gH*(Vth-VH)*dmH(Vth))
  B[2, 1] = gu_th - (1/gleak) * ((1 - wsu_mNa)*gNa*3*mNa_inf(Vth)^2*hNa_inf(Vth)*(Vth-VNa)*dmNa(Vth) + (1 - wsu_hNa)*gNa*mNa_inf(Vth)^3*(Vth-VNa)*dhNa(Vth) +
                                (1 - wsu_mCaT)*3*gCaT*mCaT_inf(Vth)^2*hCaT_inf(Vth)*(Vth-VCa)*dmCaT(Vth) + (1 - wsu_hCaT)*gCaT*mCaT_inf(Vth)^3*(Vth-VCa)*dhCaT(Vth) +
                                (1 - wsu_mKd)*gKd*4*mKd_inf(Vth)^3*(Vth-VK)*dmKd(Vth) + (1 - wsu_mKCa)*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK)*dmKCa(Vth, Ca_th) +
                                (1 - wsu_mKCa2)*gKCa*4*mKCa_inf(Vth, Ca_th)^3*(Vth-VK).*dmKCa2(Vth, Ca_th).*dCa(Vth, Ca_th, gCaT, 10.) + (1 - wsu_mH)*gH*(Vth-VH)*dmH(Vth))

  # Solving the linear system
  g_sol = \(A, B)
  return (g_sol[1], g_sol[2])
end
