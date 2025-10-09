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

global M = 30 
global N = 30

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
V0_01 = SOC_01*(max_V1 - min_V1) + min_V1  
V0_02 = SOC_02*(max_V2 - min_V2) + min_V2 

# Inflow [m3/hr] 
q1 = s2hr*inflow_s.bon_inflow_m3s         # historic downstream inflow to 01 Bonneville [m3/s] --> [m3/hr]
q2 = s2hr*inflow_s.tda_inflow_m3s         # historic upstream inflow to 02 Dalles Dam [m3/s] --> [m3/hr]

println("--- DATA LOAD COMPLETE ---")

### ----------- SIMULATIONS ----------- ###
 
## ----------- DATA SAMPLING FOR CONVEX HULL APPROXIMATION ----------- ##

# Independent Sample Vectors
V1_sample = range(min_V1, max_V1, length = M)
u1_sample = range(min_ut1, max_ut1, length = N)

V2_sample = range(min_V2, max_V2, length = M)
u2_sample = range(min_ut2, max_ut2, length = N)

# Precomputed Matrix of Volume x Water Release Products
Y1 = [u * (V^b1) for V in V1_sample, u in u1_sample] 
Y2 = [u * (V^b2) for V in V2_sample, u in u2_sample]

# Create the optimization model
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "NumericFocus", 3)

## ----------- SIMULATIONS ----------- ##

### Unit 01: Bonneville Dam (Downstream)
## Define variables
@variable(model, s1[1:T] >= 0) # Spill Outflow [m3/hr]
@variable(model, p1[1:T] >= 0) # Power Generation
@variable(model, lam1[1:M, 1:N, 1:T] >= 0)  # lambda matrix


# Enhanced binary variables for 2D SOS2 constraints
@variable(model, z_block[1:M-1, 1:N-1, 1:T], Bin)  # Binary for each 2x2 block selection

# 2D SOS2 constraints - only one 2x2 block can be active at a time
for t in 1:T
    # Only one block can be selected
    @constraint(model, sum(z_block[m, n, t] for m in 1:M-1, n in 1:N-1) <= 1)
    
    # Lambda values can only be nonzero within the selected block
    for m in 1:M, n in 1:N
        # Determine which blocks contain this (m,n) point
        blocks_containing_mn = []
        
        # Check all possible blocks that could contain point (m,n)
        for bm in max(1, m-1):min(M-1, m), bn in max(1, n-1):min(N-1, n)
            if bm <= M-1 && bn <= N-1 && m <= bm+1 && n <= bn+1
                push!(blocks_containing_mn, (bm, bn))
            end
        end
        
        # Lambda can only be active if at least one containing block is active
        if !isempty(blocks_containing_mn)
            @constraint(model, lam1[m, n, t] <= sum(z_block[bm, bn, t] for (bm, bn) in blocks_containing_mn))
        else
            # This shouldn't happen with proper indexing, but as safety
            @constraint(model, lam1[m, n, t] == 0)
        end
    end
end


## Expressions
@expression(model, V1[t=1:T], sum(lam1[m,n,t] * V1_sample[m] for m=1:M, n=1:N))
@expression(model, u1[t=1:T], sum(lam1[m,n,t] * u1_sample[n] for m=1:M, n=1:N))
@expression(model, y1[t=1:T], sum(lam1[m,n,t] * Y1[m,n]      for m=1:M, n=1:N))

## Initial conditions (t = 1)
@constraint(model, MassBalInit1, V1[1] == V0_01)
@constraint(model, RampRateInit1, u1[1] == min_ut1)

## Constraints
@constraint(model, MassBal1[t in 2:T], V1[t] == V1[t-1] + q1[t] - u1[t] - s1[t])
@constraint(model, Release1[t in 2:T], min_ut1 <= u1[t] <= max_ut1)
@constraint(model, ReleaseEnergy1[t in 1:T], p1[t] == (eta * g * rho_w * a1 * y1[t])/(3.6e9))
@constraint(model, RampRate1[t in 2:T], RR_dn1 <= u1[t] - u1[t-1] <= RR_up1)
@constraint(model, FeederCap1[t in 1:T], 0 <= p1[t] <= F1)
@constraint(model, Volume1[t in 1:T], min_V1 <= V1[t] <= max_V1)
@constraint(model, lambda1[t in 1:T], sum(lam1[:,:,t]) == 1)

### Unit 02: Dalles Dam (Downstream)
## Define variables
@variable(model, s2[1:T] >= 0) # Spill Outflow [m3/hr]
@variable(model, p2[1:T] >= 0)
@variable(model, lam2[1:M, 1:N, 1:T] >= 0)  # lambda matrix

## Expressions
@expression(model, V2[t=1:T], sum(lam2[m,n,t] * V2_sample[m] for m=1:M, n=1:N))
@expression(model, u2[t=1:T], sum(lam2[m,n,t] * u2_sample[n] for m=1:M, n=1:N))
@expression(model, y2[t=1:T], sum(lam2[m,n,t] * Y2[m,n]      for m=1:M, n=1:N))

## Initial conditions
@constraint(model, MassBalInit2, V2[1] == V0_02)
@constraint(model, RampRateInit2, u2[1] == min_ut2)

## Constraints
@constraint(model, MassBal2[t in 2:T], V2[t] == V2[t-1] + q2[t] - u2[t] - s2[t])
@constraint(model, Release2[t in 2:T], min_ut2 <= u2[t] <= max_ut2)
@constraint(model, ReleaseEnergy2[t in 1:T], p2[t] == (eta * g * rho_w * a2 * y2[t])/(3.6e9))
@constraint(model, RampRate2[t in 2:T], RR_dn2 <= u2[t] - u2[t-1] <= RR_up2)
@constraint(model, FeederCap2[t in 1:T], 0 <= p2[t] <= F2)
@constraint(model, Volume2[t in 1:T], min_V2 <= V2[t] <= max_V2)
@constraint(model, lambda2[t in 1:T], sum(lam2[:,:,t]) == 1)

# Objective function
@objective(model, Max, sum(p1 + p2)) 

# Solve the optimization problem
optimize!(model)

# Total Generation
obj = objective_value(model);

# Testbed Variable Save (replaces return)
s1 = value.(s1)
lam1 = value.(lam1)
V1 = value.(V1)
u1 = value.(u1)
p1 = value.(p1)
s2 = value.(s2)
lam2 = value.(lam2)
V2 = value.(V2) 
u2 = value.(u2)
p2 = value.(p2)

println("Objective: " * string(obj))

# model_report(model)
variable_report("CHA", obj, p1, u1, s1, p2, u2, s2)

println("--- SIMULATION COMPLETE ---")

# -----------------  PLOTS  ----------------- #

printplot = false 

head1 = a1 .* (V1.^b1)
p1_max =  (eta * g * rho_w * u1 .* head1)/(3.6e9)

head2 = a2 .* (V2.^b2)
p2_max =  (eta * g * rho_w * u2 .* head2)/(3.6e9)

if printplot
    # Create directory for this run 
    dir = "./plots/" ;
    stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS ") * method;
    path = dir * stamp;
    mkdir(path)

    sim_plots(path, "Unit1", T, u1, s1, p1, V1, q1, head1, F1, p1_max, min_ut1, max_ut1, min_h1, max_h1)
    sim_plots(path, "Unit2", T, u2, s2, p2, V2, q2, head2, F1, p2_max, min_ut2, max_ut2, min_h2, max_h2)
end 

# -----------------  CLAUDE DIAGNOSTICS  ----------------- #

# Check the approximation error for Unit 01
approximation_error(u1, V1, y1, b1)

# Check if decision variables are in bounds and how close they are to sample grid points
# check_bounds(V1, u1, V1_sample, u1_sample)

# Print out lambda distribution 
lambda_distribution(lam1, V1_sample, u1_sample)