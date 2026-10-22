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
include("/Users/elizacohn/Desktop/cascaded-hydro/julia/dataload.jl")
include("/Users/elizacohn/Desktop/cascaded-hydro/julia/plots.jl")
include("/Users/elizacohn/Desktop/cascaded-hydro/julia/methods.jl")

# -----------------  STATIC PARAMETERS  ----------------- #

global s2hr = 3600       # seconds in an hour (delta t)
global eta = .9          # efficiency of release-energy conversion
global rho_w = 1000      # density of water [kg/m^3]
global g = 9.8           # acceleration due to gravity [m/s^2]

# ----------------- UNIT 1 SITE PARAMETERS (BONNEVILE)  ----------------- #

scale1 = 1
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

scale2 = 1  # scale resecale the reservoir
global a2 = 4.1;                            # hydraulic head parameter 1 
global b2 = 0.134;                          # hydraulic head parameter 2 
global min_ut2 = s2hr*cfs_to_m3s(51700)     # min water release rate  [m3/hr] 
global max_ut2 = scale2*s2hr*cfs_to_m3s(161600)    # max water release rate [m3/hr] 
global min_h2 = ft_to_m(155)                # min forebay elevation levels [m]
global max_h2 = ft_to_m(165)                # max forebay elevation levels [m]
global min_V2 = (min_h2/a2)^(1/b2)          # min volume [m3]
global max_V2 = (max_h2/a2)^(1/b2)          # max volume [m3]
global RR_dn2 = s2hr*cfs_to_m3s(-4900)      # down ramp rate limit [m3/hr] 
global RR_up2 = s2hr*cfs_to_m3s(4700)       # up ramp rate limit [m3/hr]
global F2 = 1780                            # nameplate capacity [MW]

## --------------- DDU SEASONAL OLS PARAMETERS --------------- ##

wet_params = (
    center_day = 157,    # DoY
    inflow_mean =  2.58e7,
    inflow_std = 9.21e6,
    outflow_mean = 2.37e7,
    outflow_std = 9.12e6,
    # -- DDU Params -- #
    constant = -0.0018,
    coef1 = 0.836,   # inflow_lag1 [m3/hr]
    coef2 = 0.165,   # outflow_lag1 [m3/hr]
    resid_var = 0.0959,
    # -- DIU Params -- #
    AR_const = -0.0004,
    AR_coef = 0.995,    # inflow_lag1 [m3/hr]
    AR_resid_var = 0.1042,
    # -- ARCH-X Params -- #
    omega = 0.172,   
    alpha = 0.815,   
    gamma = 0.00832,  
    error_mean = 0.069,  
    error_std = 0.0364
)

dry_params = (
    center_day = 275,    # DoY
    inflow_mean =  1.06e7,
    inflow_std = 2.47e6,
    outflow_mean = 9.69e6,
    outflow_std = 2.79e6,
    # -- DDU Params -- #
    constant = -0.0013,
    coef1 = 0.830,   # inflow_lag1 [m3/hr]
    coef2 = 0.169,   # outflow_lag1 [m3/hr]
    resid_var = 0.295,
    # -- DIU Params -- #
    AR_const = 0.0020,
    AR_coef = 0.950,    # inflow_lag1 [m3/hr]
    AR_resid_var = 0.3188,
    # -- ARCH-X Params -- #
    omega = 0.179,   
    alpha = 0.796,   
    gamma = 0.0404,   
    error_mean = 0.21,  
    error_std = 0.11
)

## ----------- HYDRAULIC HEAD SUB-INTERVALS ----------- ##

# Function to create N sub-intervals and their midpoint references
function create_intervals_and_references(min_h, max_h, N)
    # Create N+1 breakpoints to define N intervals
    breakpoints = range(min_h, max_h, length = N + 1)
    breakpoints = collect(breakpoints)
    
    # Create interval bounds
    left_bounds = breakpoints[1:end-1]   # First N points
    right_bounds = breakpoints[2:end]    # Last N points
    
    # Create reference vector with midpoints
    reference_values = [right_bounds[i] for i in 1:N]

    return left_bounds, right_bounds, reference_values
end


# Create intervals and references
N = 20
global h1_lbounds, h1_rbounds, h1_refvals = create_intervals_and_references(min_h1, max_h1, N)
global h2_lbounds, h2_rbounds, h2_refvals = create_intervals_and_references(min_h2, max_h2, N)


## ----------- SIMULATION SETTINGS ----------- ##

## Seasonality
season = "DRY"
# season = "WET"

## Solution Method
# method = "MINLP"
method = "PWL"

## Uncertainty Framework
framework = "DET"
# framework = "DIU"
# framework = "DDU"

# -----------------  DATA LOAD  ----------------- #

println("--- DATA LOAD BEGIN ---")

if season == "DRY"
    params = dry_params
end

if season == "WET"
    params = wet_params
end

_, inflow, _ = fullsim_dataload();

# Filter Dataset
D = 0    
global T = 12 + 24*D   # Number of simulation hours
global lag = 1  # Number of lag terms in OLS model
year = 2022     # Simulation year 

#sim_center_date = DateTime("2023-06-01T00:00:00")
dt = Date(year) + Day(params.center_day - 1)      
sim_center_date = DateTime(string(dt, "T00:00:00"))
start_date = sim_center_date - Hour(T/2) - Hour(lag)
end_date   = sim_center_date + Hour(T/2 - 1)
inflow_s = inflow[(inflow.datetime .>= start_date) .&& (inflow.datetime .<= end_date), :]

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


# -----------------  PLOT DIRECTORY ----------------- #

make_dir = true 

if make_dir
    dir = "./plots/" ;
    stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS ") * " " * method * " " * season * " " * framework * " T = " * string(T);
    path = dir * stamp;
    mkdir(path)
end

## ----------- SIMULATIONS ----------- ##

monte_carlo = false 
M = 20
q1_M = zeros(Float64, T, M)

if monte_carlo
    for m = 1:M

        model, obj, V1, p1, u1, s1, q1_pred, std_hat, V2, p2, u2, s2 = simulation_loop(q1, q2, framework, method, params)
        variable_report(method, framework, season, obj, p1, u1, s1, p2, u2, s2, N)

        # Save simulation behavior
        q1_M[:, m] = q1_pred

    end 

    # Observe how the generation profile / volume changes as we scale up or down the upstream Unit
    # Test number 1: increasing the maximum outflow in upstream unit 
    # how does it change during the wet or dry season?

    ## Plot: Monte Carlo Inflow
    monte_carlo_plot(path, "Unit1", T, q1_M, q1[2:end], scale2)

else
    println("--- SIMULATION BEGIN ---")

    model, obj, V1, p1, u1, s1, q1_pred, std_hat, V2, p2, u2, s2 = simulation_loop(q1, q2, framework, method, params)

    println("Objective: " * string(obj))

    println("--- SIMULATION COMPLETE ---")

    # Multi-period Nonlinear Framework 
    # model, obj, s1, V1, u1, p1, s2, V2, u2, p2 = MINLP()

    # model_report(model)
    variable_report(method, framework, season, obj, p1, u1, s1, p2, u2, s2, N)

    printplot = true

    if printplot
        head1 = a1 .* (V1.^b1)
        p1_max =  (eta * g * rho_w * u1 .* head1)/(3.6e9)

        head2 = a2 .* (V2.^b2)
        p2_max =  (eta * g * rho_w * u2 .* head2)/(3.6e9)

        sim_plots(path, "Unit1", T, u1, s1, p1, V1, q1_pred, q1[(1 + lag):end], head1, F1, p1_max, min_ut1, max_ut1, min_h1, max_h1)
        sim_plots(path, "Unit2", T, u2, s2, p2, V2, q2[(1 + lag):end], q2[(1 + lag):end], head2, F2, p2_max, min_ut2, max_ut2, min_h2, max_h2)
    end 

end



