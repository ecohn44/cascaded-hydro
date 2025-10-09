using JuMP
using Gurobi
using MAT
using Printf
using CSV
using DataFrames
using Dates
using Plots
using Base.Filesystem
using Ipopt
using LaTeXStrings
using BenchmarkTools
include("/Users/elizacohn/Desktop/cascaded-hydro/dataload.jl")
include("/Users/elizacohn/Desktop/cascaded-hydro/plots.jl")
include("/Users/elizacohn/Desktop/cascaded-hydro/methods.jl")

# -----------------  STATIC PARAMETERS  ----------------- #

global s2hr = 3600       # seconds in an hour (delta t)
global eta = .9          # efficiency of release-energy conversion
global rho_w = 1000      # density of water [kg/m^3]
global g = 9.8           # acceleration due to gravity [m/s^2]

# ----------------- UNIT 1 SITE PARAMETERS (BONNEVILE)  ----------------- #

global a1 = 0.0026;                         # hydraulic head parameter 1 
global b1 = 0.4404;                         # hydraulic head parameter 2 
global min_ut1 = s2hr*cfs_to_m3s(30600)     # min water release rate  [m3/hr] 
global max_ut1 = s2hr*cfs_to_m3s(134100)    # max water release rate [m3/hr] 
global min_h1 = ft_to_m(70)                 # min forebay elevation levels [m]
global max_h1 = ft_to_m(77)                 # max forebay elevation levels [m]
global min_V1 = (min_h1/a1)^(1/b1)          # min volume [m3]
global max_V1 = (max_h1/a1)^(1/b1)          # max volume [m3]
global RR_dn1 = s2hr*cfs_to_m3s(-6300)      # down ramp rate limit [m3/hr] 
global RR_up1 = s2hr*cfs_to_m3s(6000)       # up ramp rate limit [m3/hr]
global F1 = 1154                            # nameplate capacity [MW]

# ----------------- UNIT 2 SITE PARAMETERS (DALLES)  ----------------- #

global a2 = 4.1;                            # hydraulic head parameter 1 
global b2 = 0.134;                          # hydraulic head parameter 2 
global min_ut2 = s2hr*cfs_to_m3s(51700)     # min water release rate  [m3/hr] 
global max_ut2 = s2hr*cfs_to_m3s(161600)    # max water release rate [m3/hr] 
global min_h2 = ft_to_m(155)                # min forebay elevation levels [m]
global max_h2 = ft_to_m(165)                # max forebay elevation levels [m]
global min_V2 = (min_h2/a2)^(1/b2)          # min volume [m3]
global max_V2 = (max_h2/a2)^(1/b2)          # max volume [m3]
global RR_dn2 = s2hr*cfs_to_m3s(-4900)      # down ramp rate limit [m3/hr] 
global RR_up2 = s2hr*cfs_to_m3s(4700)       # up ramp rate limit [m3/hr]
global F2 = 1780                            # nameplate capacity [MW]

## ----------- DATA SAMPLING FOR CONVEX HULL APPROXIMATION ----------- ##

global M = 7 
global N = 7

# -----------------  DATA LOAD  ----------------- #
println("--- DATA LOAD BEGIN ---")

gage, inflow, storage = fullsim_dataload();

# Filter Dataset
start_date = DateTime("2023-06-01T00:00:00")
end_date   = DateTime("2023-06-01T06:00:00") # Run for T = 6 hrs
inflow_s = inflow[(inflow.datetime .>= start_date) .&& (inflow.datetime .<= end_date), :]
global T = nrow(inflow_s);

# Storage Levels [m3]
SOC_01 = 0.5
SOC_02 = 0.5 

# Initial Conditions
global V0_01 = SOC_01*(max_V1 - min_V1) + min_V1  
global V0_02 = SOC_02*(max_V2 - min_V2) + min_V2 

# Inflow [m3/hr] 
q1 = s2hr*inflow_s.bon_inflow_m3s         # historic downstream inflow to 01 Bonneville [m3/s] --> [m3/hr]
q2 = s2hr*inflow_s.tda_inflow_m3s         # historic upstream inflow to 02 Dalles Dam [m3/s] --> [m3/hr]

println("--- DATA LOAD COMPLETE ---")


## ----------- SIMULATIONS ----------- ##
 
method = "MINLP"
# method = "CHA"

framework = "DET"
# framework = "DUI"

println("--- SIMULATION BEGIN ---")

if method == "CHA"
    model, obj, s1, lam1, V1, u1, p1, s2, lam2, V2, u2, p2 = convex_hull_approx()
elseif method == "MINLP"
    # model, obj, s1, V1, u1, p1, s2, V2, u2, p2 = MINLP()
    model, obj, V1, p1, u1, s1, q1_pred, V2, p2, u2, s2 = MINLP_loop(q1, q2) #, framework)
end

println("Objective: " * string(obj))

# model_report(model)
variable_report(method, obj, p1, u1, s1, p2, u2, s2)

println("--- SIMULATION COMPLETE ---")

# -----------------  PLOTS  ----------------- #

printplot = false 

head1 = a1 .* (V1.^b1)
p1_max =  (eta * g * rho_w * u1 .* head1)/(3.6e9)

#head2 = a2 .* (V2.^b2)
#p2_max =  (eta * g * rho_w * u2 .* head2)/(3.6e9)

if printplot
    # Create directory for this run 
    dir = "./plots/" ;
    stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS ") * method;
    path = dir * stamp;
    mkdir(path)

    sim_plots(path, "Unit1", T, u1, s1, p1, V1, q1, head1, F1, p1_max, min_ut1, max_ut1, min_h1, max_h1)
    sim_plots(path, "Unit2", T, u2, s2, p2, V2, q2, head2, F1, p2_max, min_ut2, max_ut2, min_h2, max_h2)
end 

