using LaTeXStrings
using Plots

function sim_plots(path, label, T, u, p, V, q, PF, hh)
    font = 10
    xfont = 8

    ## Plot 1: Water release
    plot1 = plot(1:T, u, label="Outflow", lw=2, legend = :topright)
    # Add maximum th 
    hline!(plot1, [max_ut], color=:red, linestyle=:dash, label="Max")
    # Add min th
    hline!(plot1, [min_ut], color=:red, linestyle=:dash, label="Min")
    xlabel!(plot1, "Hour",xguidefontsize=xfont)
    ylabel!(plot1, "Water Release (m3)")
    title!(plot1, "Hourly Water Release", titlefontsize=font)

    ## Plot 2: Hydropower generation 
    plot2 = plot(1:T, p, label="Generation", lw=2,legend = :topright, ylim=(0, PF))
    # Add feeder capacity
    hline!(plot2, [PF], color=:red, linestyle=:dash, label="Capacity")
    # Add hydraulic head
    plot!(1:T, hh, label="H. Head", lw=2,legend = :topright, color=:red, linestyle=:dash,)
    xlabel!(plot2, "Hour",xguidefontsize=xfont)
    ylabel!(plot2, "MWh")
    title!(plot2, "Hydropower Generation", titlefontsize=font)

    ## Plot 3: Volume 
    plot3 = plot(1:T, V, lw=2, legend=false)
    xlabel!(plot3, "Hour",xguidefontsize=xfont)
    ylabel!(plot3, "Reservoir Volume (m3/hr)")
    title!(plot3, "Volume", titlefontsize=font)

    ## Plot 4: Inflow
    plot4 = plot(1:T, q, lw=2, legend=false)
    xlabel!(plot4, "Hour",xguidefontsize=xfont)
    ylabel!(plot4, "System Inflow (m3/hr)")
    title!(plot4, "Inflow", titlefontsize=font)
    
    fig = plot(plot1, plot2, plot3, plot4, layout=(2, 2))
    savefig(fig,  path * "/" * label * ".png");

end


function ramp_rate(path, T, u)
    plot(u[2:end]-u[1:end-1], label = "Simulated Ramp Rate", legend = :outertopright)   
    hline!([RR_up], color=:red, linestyle=:dash, label="Max Ramp Rate")
    hline!([RR_dn], color=:red, linestyle=:dash, label="Min Ramp Rate")
    title!("Relaxed Ramp Rate for Simulated Policy")
    savefig(path * "/ramp_rate.png");
end 


function hhead_plots(path, T, hh, hh_d)
    
    ## Plot 1: Hydraulic Head
    plot1 = plot(1:T, hh, label="Hydraulic Head", lw=2, legend = :topright)
    xlabel!(plot1, "Hour")
    ylabel!(plot1, "Hydraulic Head [m]")
    title!(plot1, "Available Hydraulic Head [m]")
    savefig(plot1, path * "/hh.png");

    ## Plot 2: Hydraulic Head Derivative
    plot2 = plot(1:T, hh_d, label="Derivative of Hydraulic Head", lw=2, legend = :topright)
    xlabel!(plot2, "Hour")
    ylabel!(plot2, "Derivative of Hydraulic Head")
    title!(plot2, "Derivative of Hydraulic Head")
    savefig(plot2, path * "/hh_d.png");

end

