#=
This file contains functions to extract characteristics of the firing pattern
as well as some functions to plot complicated graphs
=#

using Statistics, Plots, StatsPlots, LaTeXStrings, Printf

## Functions extracting characteristics of the firing pattern

# This function extracts the spiking frequency of a spiking firing pattern
function extract_frequency(V, t)
    # Defining thresholds
    spike_up_threshold = 10.
    spike_down_threshold = 0.

    # Detecting spikes
    spike_detected = 0
    spike_times = []
    for i in 1:length(V)
        if V[i] > spike_up_threshold && spike_detected == 0 # Start of spike
            append!(spike_times, t[i])
            spike_detected = 1
        end
        if V[i] < spike_down_threshold && spike_detected == 1 # End of spike
            spike_detected = 0
        end
    end

    # If the neuron is silent
    if length(spike_times) < 2
        return NaN
    end

    # Calculating all interspike intervals
    ISI=[]
    for i in 2 : length(spike_times)
        append!(ISI, spike_times[i] - spike_times[i-1])
    end

    # If the neuron is silent
    if length(ISI) < 2
        return NaN
    end

    # Computing the spiking frequency
    T = mean(ISI) / 1000 # in seconds
    f = 1 / T # in Hz

    return f
end

# This function extracts characteristics of a bursting firing pattern
function extract_burstiness(V, t)
    # Defining thresholds
    spike_up_threshold = 10.
    spike_down_threshold = 0.

    # Detecting spikes
    spike_detected = 0
    spike_times = []
    for i in 1 : length(V)
        if V[i] > spike_up_threshold && spike_detected == 0 # Start of spike
            append!(spike_times, t[i])
            spike_detected = 1
        end
        if V[i] < spike_down_threshold && spike_detected == 1 # End of spike
            spike_detected = 0
        end
    end

    # If the neuron is silent
    if length(spike_times) < 3
        return NaN, NaN, NaN, NaN
    end

    # Calculating all interspike intervals
    ISI = []
    for i in 2 : length(spike_times)
        append!(ISI, spike_times[i] - spike_times[i-1])
    end

    # Defining a threshold to separate intraburst from interburst ISI
    max_ISI = maximum(ISI)
    min_ISI = minimum(ISI)
    half_ISI = (max_ISI+min_ISI)/2

    # If ISI too constant, neuron is spiking
    if max_ISI - min_ISI < 25
        return NaN, NaN, NaN, NaN
    end

    # Detecting the first spike of a burst
    first_spike_burst = findall(x -> x > half_ISI, ISI)

    # Computing the interburst frequency
    Ts = ISI[first_spike_burst]
    interburst_T = mean(Ts) / 1000 # in seconds
    interburst_f = 1 / interburst_T # in Hz

    # Computing the number of spikes per burst
    nb_spike_burst = []
    for i in 2 : length(first_spike_burst)
        append!(nb_spike_burst, first_spike_burst[i] - first_spike_burst[i-1])
    end

    # If spiking
    if length(nb_spike_burst) < 2
        return NaN, NaN, NaN, NaN
    end
    nb_spike_per_burst = round(mean(nb_spike_burst))

    # If no bursting
    if nb_spike_per_burst < 1.5 || nb_spike_per_burst > 500
        burstiness = NaN
        intraburst_f = NaN
        nb_spike_per_burst = NaN
        interburst_f = NaN
    else # Else, bursting: computing the intraburst frequency
        intra_spike_burst = findall(x -> x < half_ISI, ISI)
        Ts_intraburst = ISI[intra_spike_burst]
        T_intraburst = mean(Ts_intraburst) / 1000 # in seconds
        intraburst_f = 1 / T_intraburst # in Hz

        burstiness = (nb_spike_per_burst * intraburst_f) / interburst_T
    end

    return burstiness, nb_spike_per_burst, intraburst_f, interburst_f
end

## Functions to plot complicated graphs

# This function plots main direction of any dimensionality reduction techniques
# in a heatmap way
function heatmap_dir(dir_val, nb_channels)
    # Creating the first bin for the highest variance direction
    bin = 1
    val = abs(dir_val.values[nb_channels - bin + 1])
    str_val = @sprintf "%.2e" val
    p1 = heatmap(1:1, 1:nb_channels, reshape(dir_val.vectors[:, nb_channels - bin + 1], nb_channels, 1),
                 aspect_ratio=1, size=(200, 500), grid=false, xlabel=L"$\lambda = %$str_val$",
                 axis=false, ticks=false, c=cgrad([:dodgerblue3, :gray93, :orangered3]),
                 colorbar=false, clim=(-1, 1),
                 yticks=(1:nb_channels, ["gNa", "gCaT", "gCaS", "gA", "gKCa", "gKd", "gH", "gleak"]))

    # Second bin
    bin = bin + 1
    val = abs(dir_val.values[nb_channels - bin + 1])
    str_val = @sprintf "%.2e" val
    p2 = heatmap(1:1, 1:nb_channels, reshape(dir_val.vectors[:, nb_channels - bin + 1], nb_channels, 1),
                 aspect_ratio=1, size=(200, 500), grid=false, xlabel=L"$\lambda = %$str_val$",
                 axis=false, ticks=false, c=cgrad([:dodgerblue3, :gray93, :orangered3]),
                 colorbar=false, clim=(-1, 1))

    # Third bin
    bin = bin + 1
    val = abs(dir_val.values[nb_channels - bin + 1])
    str_val = @sprintf "%.2e" val
    p3 = heatmap(1:1, 1:nb_channels, reshape(dir_val.vectors[:, nb_channels - bin + 1], nb_channels, 1),
                 aspect_ratio=1, size=(200, 500), grid=false, xlabel=L"$\lambda = %$str_val$",
                 axis=false, ticks=false, c=cgrad([:dodgerblue3, :gray93, :orangered3]),
                 colorbar=false, clim=(-1, 1))

    # Fourth bin
    bin = bin + 1
    val = abs(dir_val.values[nb_channels - bin + 1])
    str_val = @sprintf "%.2e" val
    p4 = heatmap(1:1, 1:nb_channels, reshape(dir_val.vectors[:, nb_channels - bin + 1], nb_channels, 1),
                 aspect_ratio=1, size=(200, 500), grid=false, xlabel=L"$\lambda = %$str_val$",
                 axis=false, ticks=false, c=cgrad([:dodgerblue3, :gray93, :orangered3]),
                 colorbar=false, clim=(-1, 1))

    # Fifth bin
    bin = bin + 1
    val = abs(dir_val.values[nb_channels - bin + 1])
    str_val = @sprintf "%.2e" val
    p5 = heatmap(1:1, 1:nb_channels, reshape(dir_val.vectors[:, nb_channels - bin + 1], nb_channels, 1),
                 aspect_ratio=1, size=(200, 500), grid=false, xlabel=L"$\lambda = %$str_val$",
                 axis=false, ticks=false, c=cgrad([:dodgerblue3, :gray93, :orangered3]),
                 colorbar=false, clim=(-1, 1))

    # Sixth bin
    bin = bin + 1
    val = abs(dir_val.values[nb_channels - bin + 1])
    str_val = @sprintf "%.2e" val
    p6 = heatmap(1:1, 1:nb_channels, reshape(dir_val.vectors[:, nb_channels - bin + 1], nb_channels, 1),
                 aspect_ratio=1, size=(200, 500), grid=false, xlabel=L"$\lambda = %$str_val$",
                 axis=false, ticks=false, c=cgrad([:dodgerblue3, :gray93, :orangered3]),
                 colorbar=false, clim=(-1, 1))

    # Seventh bin
    bin = bin + 1
    val = abs(dir_val.values[nb_channels - bin + 1])
    str_val = @sprintf "%.2e" val
    p7 = heatmap(1:1, 1:nb_channels, reshape(dir_val.vectors[:, nb_channels - bin + 1], nb_channels, 1),
                 aspect_ratio=1, size=(200, 500), grid=false, xlabel=L"$\lambda = %$str_val$",
                 axis=false, ticks=false, c=cgrad([:dodgerblue3, :gray93, :orangered3]),
                 colorbar=false, clim=(-1, 1))

    # Eigth bin
    bin = bin + 1
    val = abs(dir_val.values[nb_channels - bin + 1])
    str_val = @sprintf "%.2e" val
    p8 = heatmap(1:1, 1:nb_channels, reshape(dir_val.vectors[:, nb_channels - bin + 1], nb_channels, 1),
                 aspect_ratio=1, size=(200, 500), grid=false, xlabel=L"$\lambda = %$str_val$",
                 axis=false, ticks=false, c=cgrad([:dodgerblue3, :gray93, :orangered3]),
                 colorbar=false, clim=(-1, 1))

    # Creating the colorbar bin
    colors =  -1. : 0.002 : 1.
    p90 = heatmap(1:1, colors, reshape(colors, length(colors), 1), size=(200, 500),
                 grid=false, axis=false, xticks=false, colorbar=false,
                 c=cgrad([:dodgerblue3, :gray93, :orangered3]), clim=(-1, 1), ymirror=true,
                 yticks=(-1. : 0.2 : 1.), ylims=(-1, 1), yaxis=true)

    # Arranging everything
    CC = plot(p1, p2, p3, p4, p5, p6, p7, p8, p90, layout=(1, nb_channels + 1),
              size=(1200, 500), margin=5Plots.mm)

    return CC
end

# This function plots a scatter matrix for all dimensions of the STG model
# with the two first main directions of the dimensionality reduction technique
function scatter_matrix(dir_val, nb_channels, g_all, mean_vec, Na_max, CaT_max, CaS_max,
                        A_max, KCa_max, Kd_max, H_max, leak_max, s1, s2, color_p, m_shape)

    pNaCaT = scatter(g_all[:, 1], g_all[:, 2], label="", markerstrokewidth=0.1,
                     color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_{CaT}")
    xlims!((0, Na_max))
    ylims!((0, CaT_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaCaS = scatter(g_all[:, 1], g_all[:, 3], label="", markerstrokewidth=0.1,
                     color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_{CaS}")
    xlims!((0, Na_max))
    ylims!((0, CaS_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaA = scatter(g_all[:, 1], g_all[:, 4], label="", markerstrokewidth=0.1,
                   color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_{A}")
    xlims!((0, Na_max))
    ylims!((0, A_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaKCa = scatter(g_all[:, 1], g_all[:, 5], label="", markerstrokewidth=0.1,
                     color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_{KCa}")
    xlims!((0, Na_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaKd = scatter(g_all[:, 1], g_all[:, 6], label="", markerstrokewidth=0.1,
                    color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_{Kd}")
    xlims!((0, Na_max))
    ylims!((0, Kd_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaH = scatter(g_all[:, 1], g_all[:, 7], label="", markerstrokewidth=0.1,
                   color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_{H}")
    xlims!((0, Na_max))
    ylims!((0, H_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaleak = scatter(g_all[:, 1], g_all[:, 8], label="", markerstrokewidth=0.1,
                      color=color_p, grid=false, ticks=false,
                      xrotation=45.0, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"g_{leak}")
    xlabel!(L"\bar{g}_{Na}")
    xlims!((0, Na_max))
    ylims!((0, leak_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTCaS = scatter(g_all[:, 2], g_all[:, 3], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                      color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, CaS_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTA = scatter(g_all[:, 2], g_all[:, 4], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, A_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTKCa = scatter(g_all[:, 2], g_all[:, 5], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                      color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdCaT = scatter(g_all[:, 6], g_all[:, 2], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaT_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTH = scatter(g_all[:, 2], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, H_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTleak = scatter(g_all[:, 2], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                       color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, CaT_max))
    xlabel!(L"\bar{g}_{CaT}")
    ylims!((0, leak_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pACaS = scatter(g_all[:, 4], g_all[:, 3], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaS_max))
    xlims!((0, A_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSKCa = scatter(g_all[:, 3], g_all[:, 5], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                      color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaS_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdCaS = scatter(g_all[:, 6], g_all[:, 3], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaS_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSH = scatter(g_all[:, 3], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaS_max))
    ylims!((0, H_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSleak = scatter(g_all[:, 3], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                       color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, CaS_max))
    xlabel!(L"\bar{g}_{CaS}")
    ylims!((0, leak_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAKCa = scatter(g_all[:, 4], g_all[:, 5], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, A_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdA = scatter(g_all[:, 6], g_all[:, 4], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                   color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, A_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAH = scatter(g_all[:, 4], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                  color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, A_max))
    ylims!((0, H_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAleak = scatter(g_all[:, 4], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, A_max))
    xlabel!(L"\bar{g}_{A}")
    ylims!((0, leak_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdKCa = scatter(g_all[:, 6], g_all[:, 5], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, KCa_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKCaH = scatter(g_all[:, 5], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, KCa_max))
    ylims!((0, H_max))

    plot!([mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKCaleak = scatter(g_all[:, 5], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                       color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, KCa_max))
    xlabel!(L"\bar{g}_{KCa}")
    ylims!((0, leak_max))

    plot!([mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdH = scatter(g_all[:, 6], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                   color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, Kd_max))
    ylims!((0, H_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdleak = scatter(g_all[:, 6], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                      color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, Kd_max))
    xlabel!(L"\bar{g}_{Kd}")
    ylims!((0, leak_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pHleak = scatter(g_all[:, 7], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, H_max))
    xlabel!(L"\bar{g}_{H}")
    ylims!((0, leak_max))

    plot!([mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    CC = plot(pNaKd, pNaCaT, pKdCaT, pNaA, pKdA, pCaTA, pNaCaS, pKdCaS, pCaTCaS, pACaS,
              pNaKCa, pKdKCa, pCaTKCa, pAKCa, pCaSKCa, pNaH, pKdH, pCaTH, pAH, pCaSH, pKCaH,
              pNaleak, pKdleak, pCaTleak, pAleak, pCaSleak, pKCaleak, pHleak, size =(900, 900),
              layout = @layout([° _ _ _ _ _ _; ° ° _ _ _ _ _; ° ° ° _ _ _ _;
                ° ° ° ° _ _ _; ° ° ° ° ° _ _; ° ° ° ° ° ° _; ° ° ° ° ° ° °]), margin=3Plots.mm)

    return CC
end

# This function plots a scatter matrix for all dimensions of the STG model
# with the two first main directions of the dimensionality reduction technique
function scatter_matrix_zcolor(dir_val, nb_channels, g_all, mean_vec, Na_max, CaT_max, CaS_max,
                        A_max, KCa_max, Kd_max, H_max, leak_max, s1, s2, color_p, m_shape)

    pNaCaT = scatter(g_all[:, 1], g_all[:, 2], label="", markerstrokewidth=0.1,
                     zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_{CaT}")
    xlims!((0, Na_max))
    ylims!((0, CaT_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaCaS = scatter(g_all[:, 1], g_all[:, 3], label="", markerstrokewidth=0.1,
                     zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_{CaS}")
    xlims!((0, Na_max))
    ylims!((0, CaS_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaA = scatter(g_all[:, 1], g_all[:, 4], label="", markerstrokewidth=0.1,
                   zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_{A}")
    xlims!((0, Na_max))
    ylims!((0, A_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaKCa = scatter(g_all[:, 1], g_all[:, 5], label="", markerstrokewidth=0.1,
                     zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_{KCa}")
    xlims!((0, Na_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaKd = scatter(g_all[:, 1], g_all[:, 6], label="", markerstrokewidth=0.1,
                    zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_{Kd}")
    xlims!((0, Na_max))
    ylims!((0, Kd_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaH = scatter(g_all[:, 1], g_all[:, 7], label="", markerstrokewidth=0.1,
                   zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_{H}")
    xlims!((0, Na_max))
    ylims!((0, H_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaleak = scatter(g_all[:, 1], g_all[:, 8], label="", markerstrokewidth=0.1,
                      zcolor=color_p, grid=false, ticks=false,
                      xrotation=45.0, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"g_{leak}")
    xlabel!(L"\bar{g}_{Na}")
    xlims!((0, Na_max))
    ylims!((0, leak_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTCaS = scatter(g_all[:, 2], g_all[:, 3], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                      zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, CaS_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTA = scatter(g_all[:, 2], g_all[:, 4], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, A_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTKCa = scatter(g_all[:, 2], g_all[:, 5], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                      zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdCaT = scatter(g_all[:, 6], g_all[:, 2], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaT_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTH = scatter(g_all[:, 2], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, H_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTleak = scatter(g_all[:, 2], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                       zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, CaT_max))
    xlabel!(L"\bar{g}_{CaT}")
    ylims!((0, leak_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pACaS = scatter(g_all[:, 4], g_all[:, 3], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaS_max))
    xlims!((0, A_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSKCa = scatter(g_all[:, 3], g_all[:, 5], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                      zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaS_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdCaS = scatter(g_all[:, 6], g_all[:, 3], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaS_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSH = scatter(g_all[:, 3], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaS_max))
    ylims!((0, H_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSleak = scatter(g_all[:, 3], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                       zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, CaS_max))
    xlabel!(L"\bar{g}_{CaS}")
    ylims!((0, leak_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAKCa = scatter(g_all[:, 4], g_all[:, 5], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, A_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdA = scatter(g_all[:, 6], g_all[:, 4], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                   zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, A_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAH = scatter(g_all[:, 4], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                  zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, A_max))
    ylims!((0, H_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAleak = scatter(g_all[:, 4], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, A_max))
    xlabel!(L"\bar{g}_{A}")
    ylims!((0, leak_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdKCa = scatter(g_all[:, 6], g_all[:, 5], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, KCa_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKCaH = scatter(g_all[:, 5], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, KCa_max))
    ylims!((0, H_max))

    plot!([mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKCaleak = scatter(g_all[:, 5], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                       zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, KCa_max))
    xlabel!(L"\bar{g}_{KCa}")
    ylims!((0, leak_max))

    plot!([mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdH = scatter(g_all[:, 6], g_all[:, 7], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                   zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, Kd_max))
    ylims!((0, H_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdleak = scatter(g_all[:, 6], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                      zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, Kd_max))
    xlabel!(L"\bar{g}_{Kd}")
    ylims!((0, leak_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pHleak = scatter(g_all[:, 7], g_all[:, 8], label="", markerstrokewidth=0.1, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, H_max))
    xlabel!(L"\bar{g}_{H}")
    ylims!((0, leak_max))

    plot!([mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    CC = plot(pNaKd, pNaCaT, pKdCaT, pNaA, pKdA, pCaTA, pNaCaS, pKdCaS, pCaTCaS, pACaS,
              pNaKCa, pKdKCa, pCaTKCa, pAKCa, pCaSKCa, pNaH, pKdH, pCaTH, pAH, pCaSH, pKCaH,
              pNaleak, pKdleak, pCaTleak, pAleak, pCaSleak, pKCaleak, pHleak, size =(900, 900),
              layout = @layout([° _ _ _ _ _ _; ° ° _ _ _ _ _; ° ° ° _ _ _ _;
                ° ° ° ° _ _ _; ° ° ° ° ° _ _; ° ° ° ° ° ° _; ° ° ° ° ° ° °]), margin=3Plots.mm)

    return CC
end
