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
include("/Users/elizacohn/Desktop/fpv-hydro-dispatch/dataload.jl")
include("/Users/elizacohn/Desktop/cascaded-hydro/plots.jl")

global s2hr = 3600  # seconds in an hour (delta t)
global min_ut = s2hr*cfs_to_m3s(5000)    # min daily release limit [m3/s] to [m3/hr] 
global max_ut = s2hr*cfs_to_m3s(25000)   # max daily release limit [m3/s] to [m3/hr] 
global RR_dn = s2hr*cfs_to_m3s(-2500) # down ramp rate limit [m3/s] to [m3/hr] 
global RR_up = s2hr*cfs_to_m3s(4000)  # up ramp rate limit [m3/s] to [m3/hr] 
global F = 1000   # max feeder capacity [MW]  
global eta = .775     # efficiency of release-energy conversion
global rho_w = 1000 # density of water [kg/m^3]
global g = 9.8      # acceleration due to gravity [m/s^2]
global a = 15;      # hydraulic head parameter 1 
global b = 0.2;     # hydraulic head parameter 2 


# -----------------  DATA LOAD  ----------------- #
println("--- SIMULATION BEGIN ---")

y = "22"
m = 1
T = 24 
N = 7
del = 10 # constant time delay between units

daily, _, _ = fullsim_dataload();
daily_s = filter(row -> row[:year] == y && parse(Int, row[:month]) == m, daily)

V0 = 1e-3 * daily_s.storage[1] # initial storage conditions for month 
U = sum(daily_s.release[1:N]) # monthly water contract
q1 = 1e4 * dailyflow_to_hourly(daily_s.inflow[1:N], T) # inflow
dq = zeros(T*N) # intermediary flow

## ----------- RUN BASELINE ----------- ##

# Create the optimization model
# model = Model(Gurobi.Optimizer)
model = Model(Ipopt.Optimizer)
set_silent(model) 

## Define variables
# Unit 1
@variable(model, V1[1:T*N] >= 0)
@variable(model, p1[1:T*N] >= 0)
@variable(model, u1[1:T*N])
set_upper_bound.(u1, max_ut)
set_lower_bound.(u1, min_ut)

# Unit 2
@variable(model, V2[1:T*N] >= 0)
@variable(model, p2[1:T*N] >= 0)
@variable(model, u2[1:T*N])
@variable(model, q2[1:T*N])
set_upper_bound.(u2, max_ut)
set_lower_bound.(u2, min_ut)

## Initial conditions
# Unit 1
@constraint(model, MassBalInit, V1[1] == V0)
@constraint(model, RampRateInit, u1[1] == min_ut)
# Unit 2
@constraint(model, MassBalInit2, V2[1] == V0)
@constraint(model, RampRateInit2, u2[1] == min_ut)

# Objective function
@objective(model, Max, sum(p1 + p2))

## Constraints
# Unit 1
@constraint(model, MassBal[t in 2:T*N], V1[t] == V1[t-1] + q1[t] - u1[t])
@constraint(model, ReleaseEnergy[t in 1:T*N], p1[t] <= (eta * g * rho_w * u1[t] * a * (V1[t]^b))/(3.6e9))
@constraint(model, Release[t in 2:T*N], min_ut <= u1[t] <= max_ut)
@constraint(model, RampRate[t in 2:T*N], RR_dn <= u1[t] - u1[t-1] <= RR_up)
@constraint(model, FeederCap[t in 1:T*N], 0 <= p1[t] <= F)
@constraint(model, WaterContract, sum(u1) == U)

# Unit 2
@constraint(model, Inflow2A[t in 1:del], q2[t] == dq[t])
@constraint(model, Inflow2B[t in (del+1):T*N],q2[t] == u1[t-del] + dq[t])
@constraint(model, MassBal2[t in 2:T*N], V2[t] == V2[t-1] + q2[t] - u2[t])
@constraint(model, ReleaseEnergy2[t in 1:T*N], p2[t] <= (eta * g * rho_w * u2[t] * a * (V2[t]^b))/(3.6e9))
@constraint(model, Release2[t in 2:T*N], min_ut <= u2[t] <= max_ut)
@constraint(model, RampRate2[t in 2:T*N], RR_dn <= u2[t] - u2[t-1] <= RR_up)
@constraint(model, FeederCap2[t in 1:T*N], 0 <= p2[t] <= F)
@constraint(model, WaterContract2, sum(u2) == U) 

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
# path = "./plots/04-16-2025 baseline tuning"
mkdir(path)

hh1 =  (eta * g * rho_w * value.(u1) * a .* (value.(V1).^b))/(3.6e9)
sim_plots(path, "Unit1", T, N, value.(u1), value.(p1), value.(V1), q1, F, hh1)

hh2 =  (eta * g * rho_w * value.(u2) * a .* (value.(V2).^b))/(3.6e9)
sim_plots(path, "Unit2", T, N, value.(u2), value.(p2), value.(V2), value.(q2), F, hh2)

println("--- SIMULATION COMPLETE ---")
