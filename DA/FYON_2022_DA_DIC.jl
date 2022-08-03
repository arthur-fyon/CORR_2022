#=
This file computes the dynamic input conductances of the DA model as well as
applying the compensation algorithm for a set of maximal conductances (inputs of the algorithm)
=#

include("FYON_2022_DA_kinetics.jl") # Include model gating functions
include("FYON_2022_DA_gs_derivatives.jl") # Include derivatives of model gating functions (dX/dV)

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
function DICs(V, Iapp, gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak)
  # Initializing some variables, indices and flag
  gf = zeros(length(V))
  gs = zeros(length(V))
  gu = zeros(length(V))
  gin = zeros(length(V))
  Istatic = zeros(length(V))
  ith = 1
  iosc = 1
  flag_th = 0
  flag_start = 0

  # loop over all acceptable values of the membrane voltage
  for i = 2 : length(V)
    # Defining the 3 time constants of the 3 timescales
    tau_f = tau_m(V[i])
    tau_s = tau_n(V[i])
    tau_u = 100. # Fixed because no tau for the ERG current

    # Computing the wfs(V) and wsu(V) for all gating variables
    (wfs_m, wsu_m) = var_contribution(tau_m(V[i]), tau_f, tau_s, tau_u)
    (wfs_h, wsu_h) = var_contribution(tau_h(V[i]), tau_f, tau_s, tau_u)
    (wfs_n, wsu_n) = var_contribution(tau_n(V[i]), tau_f, tau_s, tau_u)
    (wfs_mCaL, wsu_mCaL) = var_contribution(tau_mCaL(V[i]), tau_f, tau_s, tau_u)
    (wfs_mCaN, wsu_mCaN) = var_contribution(tau_mCaN(V[i]), tau_f, tau_s, tau_u)

    # ERG current dynamic in the ultraslow timescale
    wfs_oERG = 0.
    wsu_oERG = 0.

    # NMDA current dynamic in the fast timescale
    wfs_NMDA = 1.
    wsu_NMDA = 1.

    # Removing NMDA current because single cells are considered (still possible to include it)
    gNMDA = 0.

    # Computing dV_dot/dV(V) appearing in the fast timescale (Ohm's law)
    dvdot_dv = - gNa*m_inf(V[i])^3*h_inf(V[i]) - gKd*n_inf(V[i])^3 - gCaL*mCaL_inf(V[i])^2 -
                 gCaN*mCaN_inf(V[i]) - gERG*o_inf(V[i]) - gleak - gNMDA*NMDA_inf(V[i])

    # Computing fast input conductance (g_fast(V))
    gf[i] = - dvdot_dv + wfs_m*3*gNa*m_inf(V[i])^2*h_inf(V[i])*(V[i]-VNa)*dm(V[i]) + wfs_h*gNa*m_inf(V[i])^3*(V[i]-VNa)*dh(V[i]) +
              wfs_n*3*gKd*n_inf(V[i])^2*(V[i]-VK)*dn(V[i]) + wfs_mCaL*2*gCaL*mCaL_inf(V[i])*(V[i]-VCa)*dmCaL(V[i]) +
              wfs_mCaN*gCaN*(V[i]-VCa)*dmCaN(V[i]) + wfs_oERG*gERG*(V[i]-VK)*doERG(V[i]) +
              wfs_NMDA*gNMDA*(V[i]-VNMDA)*dNMDA(V[i])

    # Computing slow input conductance (g_slow(V))
    gs[i] = (wsu_m - wfs_m)*3*gNa*m_inf(V[i])^2*h_inf(V[i])*(V[i]-VNa)*dm(V[i]) + (wsu_h - wfs_h)*gNa*m_inf(V[i])^3*(V[i]-VNa)*dh(V[i]) +
            (wsu_n - wfs_n)*3*gKd*n_inf(V[i])^2*(V[i]-VK)*dn(V[i]) + (wsu_mCaL - wfs_mCaL)*2*gCaL*mCaL_inf(V[i])*(V[i]-VCa)*dmCaL(V[i]) +
            (wsu_mCaN - wfs_mCaN)*gCaN*(V[i]-VCa)*dmCaN(V[i]) + (wsu_oERG - wfs_oERG)*gERG*(V[i]-VK)*doERG(V[i]) +
            (wsu_NMDA - wfs_NMDA)*gNMDA*(V[i]-VNMDA)*dNMDA(V[i])

    # Computing ultraslow input conductance (g_uslow(V))
    gu[i] = (1 - wsu_m)*3*gNa*m_inf(V[i])^2*h_inf(V[i])*(V[i]-VNa)*dm(V[i]) + (1 - wsu_h)*gNa*m_inf(V[i])^3*(V[i]-VNa)*dh(V[i]) +
            (1 - wsu_n)*3*gKd*n_inf(V[i])^2*(V[i]-VK)*dn(V[i]) + (1 - wsu_mCaL)*2*gCaL*mCaL_inf(V[i])*(V[i]-VCa)*dmCaL(V[i]) +
            (1 - wsu_mCaN)*gCaN*(V[i]-VCa)*dmCaN(V[i]) + (1 - wsu_oERG)*gERG*(V[i]-VK)*doERG(V[i]) +
            (1 - wsu_NMDA)*gNMDA*(V[i]-VNMDA)*dNMDA(V[i])

    # Detecting the threshold voltage
    if flag_th == 0 && gf[i-1] + gs[i-1] + gu[i-1] > 0 && gf[i] + gs[i] + gu[i] <= 0
      ith = i
      flag_th = 1
    end

    # Detecting the up state voltage
    if Istatic[i-1] > 0 && Istatic[i] <= 0
      iosc = i
    end

    # Computes Istatic
    Istatic[i] = - gNa*m_inf(V[i])^3*h_inf(V[i])*(V[i]-VNa) - gKd*n_inf(V[i])^3*(V[i]-VK) -
                   gCaL*mCaL_inf(V[i])^2*(V[i]-VCa) - gCaN*mCaN_inf(V[i])*(V[i]-VCa) -
                   gERG*o_inf(V[i])*(V[i]-VK) - gleak*(V[i]-Vleak) + Iapp - gNMDA*NMDA_inf(V[i])*(V[i]-VNMDA)
  end

  # If Vth had not been found using the algorithm, assign it to default value
  if ith < 10
    ith = findall(V .== -55.5)[1][1]
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
# This function computes the values of (gNa, gCaN, gERG) that generate the values of the DICs
# at threshold voltage for the given set of conductances (gKd, gCaL, gNMDA, gleak)
function DICs_gmax_initNaCaNERG(gKd, gCaL, gNMDA, gleak, gf_th, gs_th, gu_th, Vth)
  # Defining the 3 time constants of the 3 timescales
  tau_f = tau_m(Vth)
  tau_s = tau_n(Vth)
  tau_u = 100.

  # Computing the wfs(V) and wsu(V) for all gating variables
  (wfs_m, wsu_m) = var_contribution(tau_m(Vth), tau_f, tau_s, tau_u)
  (wfs_h, wsu_h) = var_contribution(tau_h(Vth), tau_f, tau_s, tau_u)
  (wfs_n, wsu_n) = var_contribution(tau_n(Vth), tau_f, tau_s, tau_u)
  (wfs_mCaL, wsu_mCaL) = var_contribution(tau_mCaL(Vth), tau_f, tau_s, tau_u)
  (wfs_mCaN, wsu_mCaN) = var_contribution(tau_mCaN(Vth), tau_f, tau_s, tau_u)

  # ERG current dynamic in the ultraslow timescale
  wfs_oERG = 0.
  wsu_oERG = 0.

  # NMDA current dynamic in the fast timescale
  wfs_NMDA = 1.
  wsu_NMDA = 1.

  # Removing NMDA current because single cells are considered (still possible to include it)
  gNMDA = 0.

  # Computing dV_dot/dV(V) appearing in the fast timescale (Ohm's law)
  dvdot_dv = - gKd*n_inf(Vth)^3 - gCaL*mCaL_inf(Vth)^2 - gleak - gNMDA*NMDA_inf(Vth)

  # Initializing the linear system of the compensation algorithm
  A = zeros(3, 3)
  B = zeros(3, 1)

  # Filling the A matrix of the compensation algorithm (= contribution of the unknowns to the DICs)
  A[1, 1] = (1/gleak) * (m_inf(Vth)^3*h_inf(Vth) + wfs_m*3*m_inf(Vth)^2*h_inf(Vth)*(Vth-VNa)*dm(Vth) + wfs_h*m_inf(Vth)^3*(Vth-VNa)*dh(Vth))
  A[1, 2] = (1/gleak) * (mCaN_inf(Vth) + wfs_mCaN*(Vth-VCa)*dmCaN(Vth))
  A[1, 3] = (1/gleak) * (o_inf(Vth) + wfs_oERG*(Vth-VK)*doERG(Vth))
  A[2, 1] = (1/gleak) * ((wsu_m - wfs_m)*3*m_inf(Vth)^2*h_inf(Vth)*(Vth-VNa)*dm(Vth) + (wsu_h - wfs_h)*m_inf(Vth)^3*(Vth-VNa)*dh(Vth))
  A[2, 2] = (1/gleak) * ((wsu_mCaN - wfs_mCaN)*(Vth-VCa)*dmCaN(Vth))
  A[2, 3] = (1/gleak) * ((wsu_oERG - wfs_oERG)*(Vth-VK)*doERG(Vth))
  A[3, 1] = (1/gleak) * ((1 - wsu_m)*3*m_inf(Vth)^2*h_inf(Vth)*(Vth-VNa)*dm(Vth) + (1 - wsu_h)*m_inf(Vth)^3*(Vth-VNa)*dh(Vth))
  A[3, 2] = (1/gleak) * ((1 - wsu_mCaN)*(Vth-VCa)*dmCaN(Vth))
  A[3, 3] = (1/gleak) * ((1 - wsu_oERG)*(Vth-VK)*doERG(Vth))

  # Filling the B matrix of the compensation algorithm (= contribution of the non-unknowns to the DICs)
  B[1, 1] = gf_th - (1/gleak) * (- dvdot_dv + wfs_n*3*gKd*n_inf(Vth)^2*(Vth-VK)*dn(Vth) +
                                   wfs_mCaL*2*gCaL*mCaL_inf(Vth)*(Vth-VCa)*dmCaL(Vth) +
                                   wfs_NMDA*gNMDA*(Vth-VNMDA)*dNMDA(Vth))
  B[2, 1] = gs_th - (1/gleak) * ((wsu_n - wfs_n)*3*gKd*n_inf(Vth)^2*(Vth-VK)*dn(Vth) +
                                 (wsu_mCaL - wfs_mCaL)*2*gCaL*mCaL_inf(Vth)*(Vth-VCa)*dmCaL(Vth) +
                                 (wsu_NMDA - wfs_NMDA)*gNMDA*(Vth-VNMDA)*dNMDA(Vth))
  B[3, 1] = gu_th - (1/gleak) * ((1 - wsu_n)*3*gKd*n_inf(Vth)^2*(Vth-VK)*dn(Vth) +
                                 (1 - wsu_mCaL)*2*gCaL*mCaL_inf(Vth)*(Vth-VCa)*dmCaL(Vth) +
                                 (1 - wsu_NMDA)*gNMDA*(Vth-VNMDA)*dNMDA(Vth))

  # Solving the linear system
  g_sol = \(A, B)
  return (g_sol[1], g_sol[2], g_sol[3])
end

# This function computes the values of (gCaL, gCaN) that generate the values of the DICs (slow and uslow)
# at threshold voltage for the given set of conductances (gNa, gKd, gERG, gNMDA, gleak)
function DICs_gmax_neuromodCaLCaN(gNa, gKd, gERG, gNMDA, gleak, gs_th, gu_th, Vth)
  # Defining the 3 time constants of the 3 timescales
  tau_f = tau_m(Vth)
  tau_s = tau_n(Vth)
  tau_u = 100.

  # Computing the wfs(V) and wsu(V) for all gating variables
  (wfs_m, wsu_m) = var_contribution(tau_m(Vth), tau_f, tau_s, tau_u)
  (wfs_h, wsu_h) = var_contribution(tau_h(Vth), tau_f, tau_s, tau_u)
  (wfs_n, wsu_n) = var_contribution(tau_n(Vth), tau_f, tau_s, tau_u)
  (wfs_mCaL, wsu_mCaL) = var_contribution(tau_mCaL(Vth), tau_f, tau_s, tau_u)
  (wfs_mCaN, wsu_mCaN) = var_contribution(tau_mCaN(Vth), tau_f, tau_s, tau_u)

  # ERG current dynamic in the ultraslow timescale
  wfs_oERG = 0.
  wsu_oERG = 0.

  # NMDA current dynamic in the fast timescale
  wfs_NMDA = 1.
  wsu_NMDA = 1.

  # Removing NMDA current because single cells are considered (still possible to include it)
  gNMDA = 0.

  # Initializing the linear system of the compensation algorithm
  A = zeros(2, 2)
  B = zeros(2, 1)

  # Filling the A matrix of the compensation algorithm (= contribution of the unknowns to the DICs)
  A[1, 1] = (1/gleak) * ((wsu_mCaL - wfs_mCaL)*2*mCaL_inf(Vth)*(Vth-VCa)*dmCaL(Vth))
  A[1, 2] = (1/gleak) * ((wsu_mCaN - wfs_mCaN)*(Vth-VCa)*dmCaN(Vth))
  A[2, 1] = (1/gleak) * ((1 - wsu_mCaL)*2*mCaL_inf(Vth)*(Vth-VCa)*dmCaL(Vth))
  A[2, 2] = (1/gleak) * ((1 - wsu_mCaN)*(Vth-VCa)*dmCaN(Vth))

  # Filling the B matrix of the compensation algorithm (= contribution of the non-unknowns to the DICs)
  B[1, 1] = gs_th - (1/gleak) * ((wsu_m - wfs_m)*3*gNa*m_inf(Vth)^2*h_inf(Vth)*(Vth-VNa)*dm(Vth) + (wsu_h - wfs_h)*gNa*m_inf(Vth)^3*(Vth-VNa)*dh(Vth) +
                                 (wsu_n - wfs_n)*3*gKd*n_inf(Vth)^2*(Vth-VK)*dn(Vth) + (wsu_oERG - wfs_oERG)*gERG*(Vth-VK)*doERG(Vth) +
                                 (wsu_NMDA - wfs_NMDA)*gNMDA*(Vth-VNMDA)*dNMDA(Vth))
  B[2, 1] = gu_th - (1/gleak) * ((1 - wsu_m)*3*gNa*m_inf(Vth)^2*h_inf(Vth)*(Vth-VNa)*dm(Vth) + (1 - wsu_h)*gNa*m_inf(Vth)^3*(Vth-VNa)*dh(Vth) +
                                 (1 - wsu_n)*3*gKd*n_inf(Vth)^2*(Vth-VK)*dn(Vth) + (1 - wsu_oERG)*gERG*(Vth-VK)*doERG(Vth) +
                                 (1 - wsu_NMDA)*gNMDA*(Vth-VNMDA)*dNMDA(Vth))

  # Solving the linear system
  g_sol = \(A, B)
  return (g_sol[1], g_sol[2])
end
