using LaTeXStrings
using Plots
using Plots.PlotMeasures
using StatsPlots

function sim_plots(path, label, N, u, s, p, V, q_pred, q, h, PF, hh, min_ut, max_ut, min_V, max_V)
    font = 10
    xfont = 8

    M = 1e6

    ## Plot 1: Generation Outflow
    plot1 = plot(1:N, u, label="Outflow", lw=2, legend=false) #, legend = :bottomright)
    # Maximum Release 
    hline!(plot1, [max_ut], color=:red, linestyle=:dash) #, label="Max")
    # Minimum Release
    hline!(plot1, [min_ut], color=:red, linestyle=:dash) #, label="Min")
    xlabel!(plot1, "Hour",xguidefontsize=xfont)
    ylabel!(plot1, "Flow(m3/hr)")
    title!(plot1, "Generation Outflow", titlefontsize=font)

    ## Plot 2: Hydropower generation 
    plot2 = plot(1:N, p, label="Generation", lw=2, ylim=(0, PF), legend=false) #, legend = :bottomright)
    # Add feeder capacity
    hline!(plot2, [PF], color=:red, linestyle=:dash, legend=false)
    # Add hydraulic head
    # plot!(plot2, 1:N, hh, label="H. Head", lw=2,legend = :bottomright, color=:red, linestyle=:dash)
    xlabel!(plot2, "Hour",xguidefontsize=xfont)
    ylabel!(plot2, "MWh")
    title!(plot2, "Hydropower Generation", titlefontsize=font)

    ## Plot 3: Volume 
    plot3 = plot(1:N, V, lw=2, legend=false)
    xlabel!(plot3, "Hour",xguidefontsize=xfont)
    ylabel!(plot3, "Volume (m3)")
    title!(plot3, "Volume", titlefontsize=font)

    ## Plot 4: Inflow
    plot4 = plot(1:N, q_pred, lw=2, legend=false)
    plot!(plot4, q, color=:green, linestyle=:dash, legend=false)
    xlabel!(plot4, "Hour",xguidefontsize=xfont)
    ylabel!(plot4, "Flow (m3/hr)")
    title!(plot4, "Inflow", titlefontsize=font)

    ## Plot 5: Spill Outflow
    plot5 = plot(1:N, s, lw=2, legend=false)
    xlabel!(plot5, "Hour",xguidefontsize=xfont)
    ylabel!(plot5, "Flow (m3/hr)")
    title!(plot5, "Spill Outflow", titlefontsize=font)

    ## Plot 6: Hydraulic Head
    plot6 = plot(1:N, h, lw=2, legend=false)
    # Maximum Release 
    hline!(plot6, [max_V], color=:red, linestyle=:dash) #, label="Max")
    # Minimum Release
    hline!(plot6, [min_V], color=:red, linestyle=:dash) #, label="Min")
    xlabel!(plot6, "Hour",xguidefontsize=xfont)
    ylabel!(plot6, "Elevation (m)")
    title!(plot6, "Forebay Elevation", titlefontsize=font)
  
    
    fig = plot(plot1, plot2, plot5, plot3, plot4, plot6, layout=(3, 2))
    savefig(fig,  path * "/" * label * ".png");

end


function ramp_rate(u)
    plot(u[2:end]-u[1:end-1], label = "Simulated Ramp Rate", legend = :outertopright)   
    hline!([RR_up], color=:red, linestyle=:dash, label="Max Ramp Rate")
    hline!([RR_dn], color=:red, linestyle=:dash, label="Min Ramp Rate")
    title!("Relaxed Ramp Rate for Simulated Policy")
end 


function hhead_plots(path, N, hh, hh_d)
    
    ## Plot 1: Hydraulic Head
    plot1 = plot(1:N, hh, label="Hydraulic Head", lw=2, legend = :topright)
    xlabel!(plot1, "Hour")
    ylabel!(plot1, "Hydraulic Head [m]")
    title!(plot1, "Available Hydraulic Head [m]")

    ## Plot 2: Hydraulic Head Derivative
    plot2 = plot(1:N, hh_d, label="Derivative of Hydraulic Head", lw=2, legend = :topright)
    xlabel!(plot2, "Hour")
    ylabel!(plot2, "Derivative of Hydraulic Head")
    title!(plot2, "Derivative of Hydraulic Head")

end


function monte_carlo_plot(path, label, T, q1_M, q1, scale)
    
    ##Monte Carlo Inflows
    #=
    plot4 = plot(1:T, q1_M, lw=2, legend=false)
    plot!(plot4, q1, color=:green, linestyle=:dash, legend=false)
    xlabel!(plot4, "Hour")
    ylabel!(plot4, "Flow (m3/hr)")
    title!(plot4, "Inflow")
    =#

    # Calculate average predicted inflow over M simulations
    mean_soc = [mean(q1_M[t, :]) for t in 1:T]
    # Calculate onfidence spread over M simulations
    p25_soc = [quantile(q1_M[t, :], 0.25) for t in 1:T]
    p75_soc = [quantile(q1_M[t, :], 0.75) for t in 1:T]

    time_hours = 1:T
    p = plot(time_hours, p25_soc, fillrange=p75_soc, fillcolor=:green, linecolor=:lightblue, linewidth=1, label="25th-75th percentile", size=(1200, 300))
    plot!(p, time_hours, mean_soc, linewidth=2, linecolor=:blue, label="Predicted (Avg)", xlabel="Time", ylabel="Value", size=(800, 400)) 
    plot!(p, time_hours, q1, color=:red, linewidth=1.5, linestyle=:dash, legend=true, label="Observed")
    xlabel!(p, "Hour")
    ylabel!(p, "Flow (m3/hr)")
    title!(p, "Inflow Spread over M Simulations, Upstream Scale x" * string(scale))

    savefig(p, path * "/" * label * "MC.png")

end