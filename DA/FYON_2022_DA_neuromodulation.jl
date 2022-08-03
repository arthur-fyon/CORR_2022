#=
This file generates an initial set of DA neurons using the DA model as well as
neuromodulating this set to tune its firing pattern
=#

# Include DICs computation
include("FYON_2022_DA_DIC.jl") # Include the computation of DICs and compensation algorithms

## Functions to generate the initial set of DA neurons
# This function create 'ncells' DA neurons with a maximal variability in terms of
# maximal conductances while maintaining a stable firing pattern (here spiking)
function degeneracy_fixDICs_neuromod(ncells, gs_th, gu_th, Vth)
  # Initialiizing some variables
  g_all_init = zeros(ncells, 7)
  ICs_th_init = zeros(ncells, 8)

  # Fitted linear relation
  gf_th = -3.89388 * gs_th - 11.05758

  # Loop over all DA neurons
  for n = 1 : ncells
    # Generating random gleak and gCal and gKd (proportional to gleak)
    gleak = 0.008667 + 0.008667 * rand(1, 1)[1]
    gCaL = gleak * (0.015 + 0.06 * rand(1, 1)[1]) / 0.013
    gKd = gleak * (6. + 4. * rand(1, 1)[1]) / 0.013
    gNMDA = gleak * (0.12) / 0.013

    # Compute gNa, gCaN and gERG using the compensation algorithm to observe stable spiking
    (gNa, gCaN, gERG) = DICs_gmax_initNaCaNERG(gKd, gCaL, gNMDA, gleak, gf_th, gs_th, gu_th, Vth)

    # Computing DICs & co of the current DA neuron
    (ith, iosc, gf, gs, gu, gin, Istatic) = DICs(V, 0., gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak)

    # Saving DICs & co at threshold and up state voltage and all maximal conductances
    ICs_th_init[n, :] = [V[ith], V[iosc], gf[ith], gs[ith], gs[iosc], gu[ith], gin[ith], Istatic[ith]]
    g_all_init[n, :] = [gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak]
  end

  return g_all_init, ICs_th_init
end

# This function create 'ncells' DA neurons with a maximal variability in terms of
# maximal conductances while maintaining a stable firing pattern (here spiking)
# and varying only active properties of the neuron (maintaining passive properties)
function degeneracy_fixDICs_neuromodDIC(ncells, gs_th, gu_th, Vth)
  # Initialiizing some variables
  g_all_init = zeros(ncells, 7)
  ICs_th_init = zeros(ncells, 8)

  # Fitted linear relation
  gf_th = -3.89388 * gs_th - 11.05758

  # Loop over all DA neurons
  for n = 1 : ncells
    # Generating random gCal and gKd (proportional to gleak)
    gleak = 0.013
    gCaL = gleak * (0.015 + 0.06 * rand(1, 1)[1]) / 0.013
    gKd = gleak * (6. + 4. * rand(1, 1)[1]) / 0.013
    gNMDA = gleak * (0.12) / 0.013

    # Compute gNa, gCaN and gERG using the compensation algorithm to observe stable spiking
    (gNa, gCaN, gERG) = DICs_gmax_initNaCaNERG(gKd, gCaL, gNMDA, gleak, gf_th, gs_th, gu_th, Vth)

    # Computing DICs & co of the current DA neuron
    (ith, iosc, gf, gs, gu, gin, Istatic) = DICs(V, 0., gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak)

    # Saving DICs & co at threshold and up state voltage and all maximal conductances
    ICs_th_init[n, :] = [V[ith], V[iosc], gf[ith], gs[ith], gs[iosc], gu[ith], gin[ith], Istatic[ith]]
    g_all_init[n, :] = [gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak]
  end

  return g_all_init, ICs_th_init
end

# This function create 'ncells' DA neurons with a maximal variability in terms of
# maximal conductances while maintaining a stable firing pattern (here spiking)
# and varying only passive properties of the neuron (maintaining active properties)
function degeneracy_fixDICs_neuromodleak(ncells, gs_th, gu_th, Vth)
  # Initialiizing some variables
  g_all_init = zeros(ncells, 7)
  ICs_th_init = zeros(ncells, 8)

  # Fitted linear relation
  gf_th = -3.89388 * gs_th - 11.05758

  # Loop over all DA neurons
  for n = 1:ncells
    # Generating random gleak and non random gCal and gKd (proportional to gleak)
    gleak = 0.008667 + 0.011 * rand(1, 1)[1]
    gCaL = gleak * (0.03) / 0.013
    gKd = gleak * (6.) / 0.013
    gNMDA = gleak * (0.12) / 0.013

    # Compute gNa, gCaN and gERG using the compensation algorithm to observe stable spiking
    (gNa, gCaN, gERG) = DICs_gmax_initNaCaNERG(gKd, gCaL, gNMDA, gleak, gf_th, gs_th, gu_th, Vth)

    # Computing DICs & co of the current DA neuron
    (ith, iosc, gf, gs, gu, gin, Istatic) = DICs(V, 0., gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak)

    # Saving DICs & co at threshold and up state voltage and all maximal conductances
    ICs_th_init[n, :] = [V[ith], V[iosc], gf[ith], gs[ith], gs[iosc], gu[ith], gin[ith], Istatic[ith]]
    g_all_init[n, :] = [gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak]
  end

  return g_all_init, ICs_th_init
end

## This function neuromodulate a set of DA neurons to switch from a certain firing pattern
# to another in a robust way using the compensation algorithm with gCaL and gCaN
function neuromodCaLCaN(ncells, g_all, ICs, gs_th, gu_th)
  # Initialiizing some variables
  g_all_neuromod = zeros(ncells, 7)
  ICs_th_neuromod = zeros(ncells, 8)

  # Loop over all DA neurons
  for n = 1:ncells
    # Extracting fixed maximal conductances as well as threshold voltage
    gNa = g_all[n, 1]
    gKd = g_all[n, 2]
    gERG = g_all[n, 5]
    gNMDA = g_all[n, 6]
    gleak = g_all[n, 7]
    Vth = ICs[n, 1]

    # Compute gCaL and gCaN using the compensation algorithm to neuromodulate firing pattern
    (gCaL, gCaN) = DICs_gmax_neuromodCaLCaN(gNa, gKd, gERG, gNMDA, gleak, gs_th, gu_th, Vth)

    # Computing DICs & co of the current DA neuron
    (ith, iosc, gf, gs, gu, gin, Istatic) = DICs(V, 0., gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak)

    # Saving DICs & co at threshold and up state voltage and all maximal conductances
    ICs_th_neuromod[n, :] = [V[ith], V[iosc], gf[ith], gs[ith], gs[iosc], gu[ith], gin[ith], Istatic[ith]]
    g_all_neuromod[n, :] = [gNa, gKd, gCaL, gCaN, gERG, gNMDA, gleak]
  end

  return g_all_neuromod, ICs_th_neuromod
end
