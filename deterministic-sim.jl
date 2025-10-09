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

# -----------------  STATIC PARAMETERS  ----------------- #

global s2hr = 3600       # seconds in an hour (delta t)
global eta = .9          # efficiency of release-energy conversion
global rho_w = 1000      # density of water [kg/m^3]
global g = 9.8           # acceleration due to gravity [m/s^2]

# ----------------- UNIT 1 SITE PARAMETERS (BONNEVILE)  ----------------- #

global min_ut1 = s2hr*cfs_to_m3s(30600)     # min water release rate  [m3/hr] 
global max_ut1 = s2hr*cfs_to_m3s(134100)    # max water release rate [m3/hr] 
global RR_dn1 = s2hr*cfs_to_m3s(-6300)      # down ramp rate limit [m3/hr] 
global RR_up1 = s2hr*cfs_to_m3s(6000)       # up ramp rate limit [m3/hr]
global a1 = 1;                        # hydraulic head parameter 1 
global b1 = 0.2;                       # hydraulic head parameter 2 
global F1 = 1154                       # nameplate capacity [MW]

# ----------------- UNIT 2 SITE PARAMETERS (DALLES)  ----------------- #

global min_ut2 = s2hr*cfs_to_m3s(51700)     # min water release rate  [m3/hr] 
global max_ut2 = s2hr*cfs_to_m3s(161600)    # max water release rate [m3/hr] 
global RR_dn2 = s2hr*cfs_to_m3s(-4900)      # down ramp rate limit [m3/hr] 
global RR_up2 = s2hr*cfs_to_m3s(4700)       # up ramp rate limit [m3/hr]
global a2 = 1;                        # hydraulic head parameter 1 
global b2 = 0.2;                       # hydraulic head parameter 2 
global F2 = 1780                       # nameplate capacity [MW]

# -----------------  DATA LOAD  ----------------- #
println("--- DATA LOAD BEGIN ---")

flow, inflow, bon, tda = fullsim_dataload();

# Filter Dataset
start_date = DateTime("2023-06-01T00:00:00")
end_date   = DateTime("2023-06-01T23:59:59") # 1 Day
flow_s = flow[(flow.datetime .>= start_date) .&& (flow.datetime .<= end_date), :]
inflow_s = inflow[(inflow.datetime .>= start_date) .&& (inflow.datetime .<= end_date), :]
N = nrow(flow_s);

# Storage Levels [m3]
V0_01 = bon.S_m3[end]
V0_02 = tda.S_m3[end]

# Historic Inflow [m3/s] (~check units)
# q1 = flow_s.down_inflow_m         # downstream inflow to 01 Bonneville
q2 = s2hr*inflow_s.inflow_m3s       # upstream inflow to 02 Dalles Dam [m3/s] --> [m3/hr]
q1 = copy(q2)                       # temp until streamflow rating curve compelte

println("--- DATA LOAD COMPLETE ---")

## ----------- RUN BASELINE ----------- ##

println("--- SIMULATION BEGIN ---")

# Create the optimization model
model = Model(Gurobi.Optimizer)
# model = Model(Ipopt.Optimizer)
# set_silent(model) 

## Define variables
# Unit 01: Bonneville Dam (Downstream)
@variable(model, V1[1:N] >= 0) # Reservoir Volume [m3]
@variable(model, p1[1:N] >= 0) # Generation [MWh]
@variable(model, u1[1:N])      # Generation Outflow [m3/hr] (~check units)
set_upper_bound.(u1, max_ut1)
set_lower_bound.(u1, min_ut1)
# To Do: Upper/Lower bounds on Volume1

# Unit 02: The Dalles Dam (Upstream)
@variable(model, V2[1:N] >= 0)
@variable(model, p2[1:N] >= 0)
@variable(model, u2[1:N])
set_upper_bound.(u2, max_ut2)
set_lower_bound.(u2, min_ut2)
# To Do: Upper/Lower bounds on Volume2

## Initial conditions
# Unit 1
@constraint(model, MassBalInit1, V1[1] == V0_01)
@constraint(model, RampRateInit1, u1[1] == min_ut1)
# Unit 2
@constraint(model, MassBalInit2, V2[1] == V0_02)
@constraint(model, RampRateInit2, u2[1] == min_ut2)

# Objective function
@objective(model, Max, sum(p1 + p2)) 

## Constraints
# Unit 1
@constraint(model, MassBal1[t in 2:N], V1[t] == V1[t-1] + q1[t] - u1[t])
# @constraint(model, ReleaseEnergy1[t in 1:N], p1[t] <= (eta * g * rho_w * u1[t] * a1 * (V1[t]^b1))/(3.6e9))
@constraint(model, ReleaseEnergy1[t in 1:N], p1[t] <= (eta * g * rho_w * u1[t] * a1 * (V0_01^b1))/(3.6e9))
@constraint(model, Release1[t in 2:N], min_ut1 <= u1[t] <= max_ut1)
@constraint(model, RampRate1[t in 2:N], RR_dn1 <= u1[t] - u1[t-1] <= RR_up1)
@constraint(model, FeederCap1[t in 1:N], 0 <= p1[t] <= F1)

# Unit 2
@constraint(model, MassBal2[t in 2:N], V2[t] == V2[t-1] + q2[t] - u2[t])
# @constraint(model, ReleaseEnergy2[t in 1:N], p2[t] <= (eta * g * rho_w * u2[t] * a2 * (V2[t]^b2))/(3.6e9))
@constraint(model, ReleaseEnergy2[t in 1:N], p2[t] <= (eta * g * rho_w * u2[t] * a2 * (V0_02^b2))/(3.6e9))
@constraint(model, Release2[t in 2:N], min_ut2 <= u2[t] <= max_ut2)
@constraint(model, RampRate2[t in 2:N], RR_dn2 <= u2[t] - u2[t-1] <= RR_up2)
@constraint(model, FeederCap2[t in 1:N], 0 <= p2[t] <= F2)

# Solve the optimization problem
optimize!(model)

# Revenue
obj = objective_value(model);
println(obj)

# -----------------  PLOTS  ----------------- #

# Create directory for this run 
stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS") ;
dir = "./plots/" ;
path = dir * stamp;
# path = "./plots/06-23-2025 14.25.56"
mkdir(path)

hh1 =  (eta * g * rho_w * value.(u1) * a1 .* (value.(V1).^b1))/(3.6e9)
sim_plots(path, "Unit1", N, value.(u1), value.(p1), value.(V1), q1, F1, hh1, min_ut1, max_ut1)

hh2 =  (eta * g * rho_w * value.(u2) * a2 .* (value.(V2).^b2))/(3.6e9)
sim_plots(path, "Unit2", N, value.(u2), value.(p2), value.(V2), q2, F2, hh2, min_ut2, max_ut2)

println("--- SIMULATION COMPLETE ---")
