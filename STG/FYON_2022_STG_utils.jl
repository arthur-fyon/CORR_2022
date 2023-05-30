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

    total_var = sum(dir_val.values)
    eig_val_decreasing = reverse(dir_val.values) ./ total_var
    # Creating the first bin for the highest variance direction
    bin = 1
    val = eig_val_decreasing[1]*100
    str_val = @sprintf "%d" val
    p1 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]),
                 colorbar=false, clim=(0, 1), tickfontsize=18, guidefontsize=15,
                 yticks=(1:nb_channels, reverse([L"\bar{g}_\mathrm{Na}", L"\bar{g}_\mathrm{CaT}", L"\bar{g}_\mathrm{CaS}",
                 L"\bar{g}_\mathrm{A}", L"\bar{g}_\mathrm{KCa}", L"\bar{g}_\mathrm{Kd}", L"\bar{g}_\mathrm{H}", L"g_\mathrm{leak}"])))

    # Second bin
    bin = bin + 1
    val = eig_val_decreasing[2]*100
    str_val = @sprintf "%d" val
    p2 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Third bin
    bin = bin + 1
    val = eig_val_decreasing[3]*100
    str_val = @sprintf "%d" val
    p3 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Fourth bin
    bin = bin + 1
    val = eig_val_decreasing[4]*100
    str_val = @sprintf "%d" val
    p4 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Fifth bin
    bin = bin + 1
    val = eig_val_decreasing[5]*100
    str_val = @sprintf "%d" val
    p5 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Sixth bin
    bin = bin + 1
    val = eig_val_decreasing[6]*100
    str_val = @sprintf "%d" val
    p6 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Seventh bin
    bin = bin + 1
    val = eig_val_decreasing[7]*100
    str_val = @sprintf "%d" val
    p7 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Eigth bin
    bin = bin + 1
    val = eig_val_decreasing[8]*100
    str_val = @sprintf "%d" val
    p8 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%", right_margin=5Plots.mm,
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Creating the colorbar bin
    colors =  -1. : 0.002 : 1.
    p90 = heatmap(1:1, colors, reshape(colors, length(colors), 1), size=(200, 500),
                 grid=false, axis=false, xticks=false, colorbar=false,
                 c=cgrad([:gray93, :orangered3]), clim=(0, 1), ymirror=true,
                 yticks=(0 : 0.2 : 1, [L"0\%", L"20\%", L"40\%", L"60\%", L"80\%", L"100\%"]), ylims=(0, 1), yaxis=true)

    # Arranging everything
    CC = plot(p1, p2, p3, p4, p5, p6, p7, p8, layout=(1, nb_channels),
              size=(600, 500))

    return CC
end

# This function plots main direction of any dimensionality reduction techniques
# in a heatmap way
function heatmap_dir4(dir_val, nb_channels)

    total_var = sum(dir_val.values)
    eig_val_decreasing = reverse(dir_val.values) ./ total_var
    # Creating the first bin for the highest variance direction
    bin = 1
    val = eig_val_decreasing[1]*100
    str_val = @sprintf "%d" val
    p1 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"$\mathrm{PC}1 (%$str_val\%)$",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]),
                 colorbar=false, clim=(0, 1), tickfontsize=18, guidefontsize=15,
                 yticks=(1:nb_channels, reverse([L"\bar{g}_\mathrm{Na}", L"\bar{g}_\mathrm{CaT}", L"\bar{g}_\mathrm{CaS}",
                 L"\bar{g}_\mathrm{A}", L"\bar{g}_\mathrm{KCa}", L"\bar{g}_\mathrm{Kd}", L"\bar{g}_\mathrm{H}", L"g_\mathrm{leak}"])))

    # Second bin
    bin = bin + 1
    val = eig_val_decreasing[2]*100
    str_val = @sprintf "%d" val
    p2 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"$\mathrm{PC}2 (%$str_val\%)$",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Third bin
    bin = bin + 1
    val = eig_val_decreasing[3]*100
    str_val = @sprintf "%d" val
    p3 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"$\mathrm{PC}3 (%$str_val\%)$",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Fourth bin
    bin = bin + 1
    val = eig_val_decreasing[4]*100
    str_val = @sprintf "%d" val
    p4 = heatmap(1:1, 1:nb_channels, reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"$\mathrm{PC}4 (%$str_val\%)$",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Creating the colorbar bin
    colors =  -1. : 0.002 : 1.
    p90 = heatmap(1:1, colors, reshape(colors, length(colors), 1), size=(200, 500),
                 grid=false, axis=false, xticks=false, colorbar=false,
                 c=cgrad([:gray93, :orangered3]), clim=(0, 1), ymirror=true,
                 yticks=(0 : 0.2 : 1, [L"0\%", L"20\%", L"40\%", L"60\%", L"80\%", L"100\%"]), ylims=(0, 1), yaxis=true)

    # Arranging everything
    CC = plot(p1, p2, p3, p4, layout=(1, 4),
              size=(600, 500), margin=5Plots.mm)

    return CC
end

# This function plots a scatter matrix for all dimensions of the STG model
# with the two first main directions of the dimensionality reduction technique
function scatter_matrix3x3(g_all, maxs, color_p, m_shape, names; flag=0, dir_val=Nothing, mean_vec=Nothing, s1=Nothing, flag2=0, g_all2=Nothing, m_shape2=Nothing, s2=Nothing) # flag = 0 --> correlation, elseif 1 PC1 #flag2 = 1 -> 2 g_all

    cors = NaN * ones(3, 3)

    p12 = scatter(g_all[:, 1], g_all[:, 2], label="", markerstrokewidth=0., color=color_p, top_margin=12Plots.mm,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 1], g_all2[:, 2], label="", color=color_p,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    annotate!(maxs[1]/2, maxs[2]*1.3, Plots.text(names[1], :black, :center, 18))

    xlims!((0, maxs[1]))
    ylims!((0, maxs[2]))

    if flag == 0
        line_12 = fit(g_all[:, 1], g_all[:, 2], 1)
        s0 = minimum(g_all[:, 1])
        sn = maximum(g_all[:, 1])
        plot!([s0, sn], [line_12(s0), line_12(sn)], linewidth=2, label="", linecolor="black")

        cors[1, 1] = cor(g_all[:, 1], g_all[:, 2])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values, mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values],
              [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values, mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end




    p13 = scatter(g_all[:, 1], g_all[:, 3], label="", markerstrokewidth=0., color=color_p,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 1], g_all2[:, 3], label="", color=color_p,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    xlims!((0, maxs[1]))
    ylims!((0, maxs[3]))

    if flag == 0
        line_13 = fit(g_all[:, 1], g_all[:, 3], 1)
        s0 = minimum(g_all[:, 1])
        sn = maximum(g_all[:, 1])
        plot!([s0, sn], [line_13(s0), line_13(sn)], linewidth=2, label="", linecolor="black")

        cors[2, 1] = cor(g_all[:, 1], g_all[:, 3])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values, mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values],
              [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values, mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end



    p14 = scatter(g_all[:, 1], g_all[:, 4], label="", markerstrokewidth=0., color=color_p,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 1], g_all2[:, 4], label="", color=color_p,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    xlims!((0, maxs[1]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_14 = fit(g_all[:, 1], g_all[:, 4], 1)
        s0 = minimum(g_all[:, 1])
        sn = maximum(g_all[:, 1])
        plot!([s0, sn], [line_14(s0), line_14(sn)], linewidth=2, label="", linecolor="black")

        cors[3, 1] = cor(g_all[:, 1], g_all[:, 4])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values, mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values],
              [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values, mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end



    p23 = scatter(g_all[:, 2], g_all[:, 3], label="", markerstrokewidth=0., color=color_p,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 2], g_all2[:, 3], label="", color=color_p,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    xlims!((0, maxs[2]))
    ylims!((0, maxs[3]))

    if flag == 0
        line_23 = fit(g_all[:, 2], g_all[:, 3], 1)
        s0 = minimum(g_all[:, 2])
        sn = maximum(g_all[:, 2])
        plot!([s0, sn], [line_23(s0), line_23(sn)], linewidth=2, label="", linecolor="black")

        cors[2, 2] = cor(g_all[:, 2], g_all[:, 3])
    elseif flag == 1
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values, mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values],
              [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values, mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end



    p24 = scatter(g_all[:, 2], g_all[:, 4], label="", markerstrokewidth=0., color=color_p,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 2], g_all2[:, 4], label="", color=color_p,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    xlims!((0, maxs[2]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_24 = fit(g_all[:, 2], g_all[:, 4], 1)
        s0 = minimum(g_all[:, 2])
        sn = maximum(g_all[:, 2])
        plot!([s0, sn], [line_24(s0), line_24(sn)], linewidth=2, label="", linecolor="black")

        cors[3, 2] = cor(g_all[:, 2], g_all[:, 4])
    elseif flag == 1
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values, mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values],
              [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values, mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end



    p34 = scatter(g_all[:, 3], g_all[:, 4], label="", markerstrokewidth=0., color=color_p, right_margin=15Plots.mm,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 3], g_all2[:, 4], label="", color=color_p,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    annotate!(maxs[3]*1.3, maxs[4]/2, Plots.text(names[4], :black, :center, 18))
    xlims!((0, maxs[3]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_34 = fit(g_all[:, 3], g_all[:, 4], 1)
        s0 = minimum(g_all[:, 3])
        sn = maximum(g_all[:, 3])
        plot!([s0, sn], [line_34(s0), line_34(sn)], linewidth=2, label="", linecolor="black")

        cors[3, 3] = cor(g_all[:, 3], g_all[:, 4])
        display(cors)
    elseif flag == 1
        plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values, mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values],
              [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values, mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end

    p21 = plot(axis=false, ticks=false, labels=false)
    xlims!((-1, 1))
    ylims!((-1, 1))
    annotate!(0, 0, Plots.text(names[2], :black, :center, 18))

    p32 = plot(axis=false, ticks=false, labels=false)
    xlims!((-1, 1))
    ylims!((-1, 1))
    annotate!(0, 0, Plots.text(names[3], :black, :center, 18))




    CC = plot(p12, p21, p13, p23, p32, p14, p24, p34, size =(500, 500),
              layout = @layout([° ° _; ° ° °; ° ° °]), margin=3Plots.mm)

    return CC
end


# This function plots a scatter matrix for all dimensions of the STG model
# with the two first main directions of the dimensionality reduction technique
function scatter_matrix(dir_val, nb_channels, g_all, mean_vec, Na_max, CaT_max, CaS_max,
                        A_max, KCa_max, Kd_max, H_max, leak_max, s1, s2, color_p, m_shape)

    pNaCaT = scatter(g_all[:, 1], g_all[:, 2], label="", markerstrokewidth=0.,
                     color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{CaT}")
    xlims!((0, Na_max))
    ylims!((0, CaT_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaCaS = scatter(g_all[:, 1], g_all[:, 3], label="", markerstrokewidth=0.,
                     color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{CaS}")
    xlims!((0, Na_max))
    ylims!((0, CaS_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaA = scatter(g_all[:, 1], g_all[:, 4], label="", markerstrokewidth=0.,
                   color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{A}")
    xlims!((0, Na_max))
    ylims!((0, A_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaKCa = scatter(g_all[:, 1], g_all[:, 5], label="", markerstrokewidth=0.,
                     color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{KCa}")
    xlims!((0, Na_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaKd = scatter(g_all[:, 1], g_all[:, 6], label="", markerstrokewidth=0.,
                    color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{Kd}")
    xlims!((0, Na_max))
    ylims!((0, Kd_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaH = scatter(g_all[:, 1], g_all[:, 7], label="", markerstrokewidth=0.,
                   color=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{H}")
    xlims!((0, Na_max))
    ylims!((0, H_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaleak = scatter(g_all[:, 1], g_all[:, 8], label="", markerstrokewidth=0.,
                      color=color_p, grid=false, ticks=false,
                      xrotation=45.0, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    ylabel!(L"g_\mathrm{leak}")
    xlabel!(L"\bar{g}_\mathrm{Na}")
    xlims!((0, Na_max))
    ylims!((0, leak_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTCaS = scatter(g_all[:, 2], g_all[:, 3], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                      color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, CaS_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTA = scatter(g_all[:, 2], g_all[:, 4], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, A_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTKCa = scatter(g_all[:, 2], g_all[:, 5], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                      color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdCaT = scatter(g_all[:, 6], g_all[:, 2], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaT_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTH = scatter(g_all[:, 2], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, H_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTleak = scatter(g_all[:, 2], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                       color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, CaT_max))
    xlabel!(L"\bar{g}_\mathrm{CaT}")
    ylims!((0, leak_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pACaS = scatter(g_all[:, 4], g_all[:, 3], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaS_max))
    xlims!((0, A_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSKCa = scatter(g_all[:, 3], g_all[:, 5], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                      color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaS_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdCaS = scatter(g_all[:, 6], g_all[:, 3], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaS_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSH = scatter(g_all[:, 3], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaS_max))
    ylims!((0, H_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSleak = scatter(g_all[:, 3], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                       color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, CaS_max))
    xlabel!(L"\bar{g}_\mathrm{CaS}")
    ylims!((0, leak_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAKCa = scatter(g_all[:, 4], g_all[:, 5], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, A_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdA = scatter(g_all[:, 6], g_all[:, 4], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                   color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, A_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAH = scatter(g_all[:, 4], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                  color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, A_max))
    ylims!((0, H_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAleak = scatter(g_all[:, 4], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, A_max))
    xlabel!(L"\bar{g}_\mathrm{A}")
    ylims!((0, leak_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdKCa = scatter(g_all[:, 6], g_all[:, 5], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, KCa_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKCaH = scatter(g_all[:, 5], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                    color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, KCa_max))
    ylims!((0, H_max))

    plot!([mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKCaleak = scatter(g_all[:, 5], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                       color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, KCa_max))
    xlabel!(L"\bar{g}_\mathrm{KCa}")
    ylims!((0, leak_max))

    plot!([mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdH = scatter(g_all[:, 6], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                   color=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, Kd_max))
    ylims!((0, H_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdleak = scatter(g_all[:, 6], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                      color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, Kd_max))
    xlabel!(L"\bar{g}_\mathrm{Kd}")
    ylims!((0, leak_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pHleak = scatter(g_all[:, 7], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape,
                     color=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, H_max))
    xlabel!(L"\bar{g}_\mathrm{H}")
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

    pNaCaT = scatter(g_all[:, 1], g_all[:, 2], label="", markerstrokewidth=0.,
                     zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{CaT}")
    xlims!((0, Na_max))
    ylims!((0, CaT_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaCaS = scatter(g_all[:, 1], g_all[:, 3], label="", markerstrokewidth=0.,
                     zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{CaS}")
    xlims!((0, Na_max))
    ylims!((0, CaS_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaA = scatter(g_all[:, 1], g_all[:, 4], label="", markerstrokewidth=0.,
                   zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{A}")
    xlims!((0, Na_max))
    ylims!((0, A_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaKCa = scatter(g_all[:, 1], g_all[:, 5], label="", markerstrokewidth=0.,
                     zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{KCa}")
    xlims!((0, Na_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaKd = scatter(g_all[:, 1], g_all[:, 6], label="", markerstrokewidth=0.,
                    zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{Kd}")
    xlims!((0, Na_max))
    ylims!((0, Kd_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaH = scatter(g_all[:, 1], g_all[:, 7], label="", markerstrokewidth=0.,
                   zcolor=color_p, grid=false, ticks=false, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"\bar{g}_\mathrm{H}")
    xlims!((0, Na_max))
    ylims!((0, H_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pNaleak = scatter(g_all[:, 1], g_all[:, 8], label="", markerstrokewidth=0.,
                      zcolor=color_p, grid=false, ticks=false,
                      xrotation=45.0, tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false), guidefontsize=18)
    ylabel!(L"g_\mathrm{leak}")
    xlabel!(L"\bar{g}_\mathrm{Na}")
    xlims!((0, Na_max))
    ylims!((0, leak_max))

    plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels], mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1], mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTCaS = scatter(g_all[:, 2], g_all[:, 3], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                      zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, CaS_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTA = scatter(g_all[:, 2], g_all[:, 4], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, A_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTKCa = scatter(g_all[:, 2], g_all[:, 5], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                      zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdCaT = scatter(g_all[:, 6], g_all[:, 2], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaT_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTH = scatter(g_all[:, 2], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaT_max))
    ylims!((0, H_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaTleak = scatter(g_all[:, 2], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                       zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, CaT_max))
    xlabel!(L"\bar{g}_\mathrm{CaT}")
    ylims!((0, leak_max))

    plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels], mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1], mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pACaS = scatter(g_all[:, 4], g_all[:, 3], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaS_max))
    xlims!((0, A_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSKCa = scatter(g_all[:, 3], g_all[:, 5], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                      zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaS_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdCaS = scatter(g_all[:, 6], g_all[:, 3], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, CaS_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSH = scatter(g_all[:, 3], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, CaS_max))
    ylims!((0, H_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pCaSleak = scatter(g_all[:, 3], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                       zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, CaS_max))
    xlabel!(L"\bar{g}_\mathrm{CaS}")
    ylims!((0, leak_max))

    plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels], mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1], mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAKCa = scatter(g_all[:, 4], g_all[:, 5], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, A_max))
    ylims!((0, KCa_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdA = scatter(g_all[:, 6], g_all[:, 4], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                   zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, A_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAH = scatter(g_all[:, 4], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                  zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, A_max))
    ylims!((0, H_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pAleak = scatter(g_all[:, 4], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, A_max))
    xlabel!(L"\bar{g}_\mathrm{A}")
    ylims!((0, leak_max))

    plot!([mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels], mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1], mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdKCa = scatter(g_all[:, 6], g_all[:, 5], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    ylims!((0, KCa_max))
    xlims!((0, Kd_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKCaH = scatter(g_all[:, 5], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                    zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, KCa_max))
    ylims!((0, H_max))

    plot!([mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKCaleak = scatter(g_all[:, 5], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                       zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, KCa_max))
    xlabel!(L"\bar{g}_\mathrm{KCa}")
    ylims!((0, leak_max))

    plot!([mean_vec[5] - s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels], mean_vec[5] + s1*dir_val.vectors[:, nb_channels][5]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[5] - s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1], mean_vec[5] + s2*dir_val.vectors[:, nb_channels-1][5]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdH = scatter(g_all[:, 6], g_all[:, 7], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                   zcolor=color_p, grid=false, ticks=false, guidefontsize=18)
    xlims!((0, Kd_max))
    ylims!((0, H_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[7] - s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels], mean_vec[7] + s1*dir_val.vectors[:, nb_channels][7]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[7] - s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1], mean_vec[7] + s2*dir_val.vectors[:, nb_channels-1][7]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pKdleak = scatter(g_all[:, 6], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                      zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, Kd_max))
    xlabel!(L"\bar{g}_\mathrm{Kd}")
    ylims!((0, leak_max))

    plot!([mean_vec[6] - s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels], mean_vec[6] + s1*dir_val.vectors[:, nb_channels][6]*dir_val.values[nb_channels]],
          [mean_vec[8] - s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels], mean_vec[8] + s1*dir_val.vectors[:, nb_channels][8]*dir_val.values[nb_channels]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)

    plot!([mean_vec[6] - s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1], mean_vec[6] + s2*dir_val.vectors[:, nb_channels-1][6]*dir_val.values[nb_channels-1]],
          [mean_vec[8] - s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1], mean_vec[8] + s2*dir_val.vectors[:, nb_channels-1][8]*dir_val.values[nb_channels-1]],
          arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)



    pHleak = scatter(g_all[:, 7], g_all[:, 8], label="", markerstrokewidth=0., tickfontsize=10, markershape=m_shape, legend=false, c=cgrad(:thermal, rev=false),
                     zcolor=color_p, grid=false, ticks=false, xrotation=45.0, guidefontsize=18)
    xlims!((0, H_max))
    xlabel!(L"\bar{g}_\mathrm{H}")
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
