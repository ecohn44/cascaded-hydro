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

a1 = 0.0026;                         # hydraulic head parameter 1 
b1 = 0.4404;                         # hydraulic head parameter 2 
min_ut1 = s2hr*cfs_to_m3s(30600)     # min water release rate  [m3/hr] 
max_ut1 = s2hr*cfs_to_m3s(134100)    # max water release rate [m3/hr] 
min_h1 = ft_to_m(70)                 # min forebay elevation levels [m]
max_h1 = ft_to_m(77)                 # max forebay elevation levels [m]
min_V1 = (min_h1/a1)^(1/b1)          # min volume [m3]
max_V1 = (max_h1/a1)^(1/b1)          # max volume [m3]
RR_dn1 = s2hr*cfs_to_m3s(-6300)      # down ramp rate limit [m3/hr] 
RR_up1 = s2hr*cfs_to_m3s(6000)       # up ramp rate limit [m3/hr]
F1 = 1154                            # nameplate capacity [MW]

# ----------------- UNIT 2 SITE PARAMETERS (DALLES)  ----------------- #

a2 = 4.1;                            # hydraulic head parameter 1 
b2 = 0.134;                          # hydraulic head parameter 2 
min_ut2 = s2hr*cfs_to_m3s(51700)     # min water release rate  [m3/hr] 
max_ut2 = s2hr*cfs_to_m3s(161600)    # max water release rate [m3/hr] 
min_h2 = ft_to_m(155)                # min forebay elevation levels [m]
max_h2 = ft_to_m(165)                # max forebay elevation levels [m]
min_V2 = (min_h2/a2)^(1/b2)          # min volume [m3]
max_V2 = (max_h2/a2)^(1/b2)          # max volume [m3]
RR_dn2 = s2hr*cfs_to_m3s(-4900)      # down ramp rate limit [m3/hr] 
RR_up2 = s2hr*cfs_to_m3s(4700)       # up ramp rate limit [m3/hr]
F2 = 1780                            # nameplate capacity [MW]

# -----------------  DATA LOAD  ----------------- #
println("--- DATA LOAD BEGIN ---")

gage, inflow, storage = fullsim_dataload();

# Filter Dataset
start_date = DateTime("2023-06-01T00:00:00")
end_date   = DateTime("2023-06-01T06:00:00") # Run for T = 6 hrs
inflow_s = inflow[(inflow.datetime .>= start_date) .&& (inflow.datetime .<= end_date), :]
T = nrow(inflow_s);

# Storage Levels [m3]
V0_01 = af_to_m3(697000) #storage.bon_S_m3[end]
V0_02 = storage.tda_S_m3[end]

# Inflow [m3/hr] 
q1 = s2hr*inflow_s.bon_inflow_m3s         # historic downstream inflow to 01 Bonneville [m3/s] --> [m3/hr]
q2 = s2hr*inflow_s.tda_inflow_m3s         # historic upstream inflow to 02 Dalles Dam [m3/s] --> [m3/hr]

println("--- DATA LOAD COMPLETE ---")

## ----------- DATA SAMPLING FOR CONVEX HULL APPROXIMATION ----------- ##

M = 10 
N = 15

# Independent Sample Vectors
V1_sample = range(min_V1, max_V1, length = M)
u1_sample = range(min_ut1, max_ut1, length = N)

V2_sample = range(min_V2, max_V2, length = M)
u2_sample = range(min_ut2, max_ut2, length = N)

# Define Hydraulic Head Function
f1(V) = a1 * (V^b1)
f2(V) = a2 * (V^b2)

# Precomputed Matrix of Power Outputs
P1 = [eta * g * rho_w * f1(V) * u for u in u1_sample, V in V1_sample] / (3.6e9)
P2 = [eta * g * rho_w * f2(V) * u for u in u2_sample, V in V2_sample] / (3.6e9)

## ----------- RUN BASELINE ----------- ##

println("--- SIMULATION BEGIN ---")

# Create the optimization model
model = Model(Gurobi.Optimizer)

### Unit 01: Bonneville Dam (Downstream)
## Define variables
@variable(model, s1[1:T] >= 0) # Spill Outflow [m3/hr]
@variable(model, lam1[1:N, 1:M, 1:T] >= 0)  # lambda matrix

## Initial conditions
@constraint(model, MassBalInit1, sum(lam1[i,j,1] * V1_sample[j] for i=1:N, j=1:M) == V0_01)
@constraint(model, RampRateInit1, sum(lam1[i,j,1] * u1_sample[i] for i=1:N, j=1:M) == min_ut1)

## Expressions
@expression(model, V1[t=1:T], sum(lam1[i,j,t] * V1_sample[j] for i=1:N, j=1:M))
@expression(model, u1[t=1:T], sum(lam1[i,j,t] * u1_sample[i] for i=1:N, j=1:M))
@expression(model, p1[t=1:T], sum(lam1[i,j,t] * P1[i,j]      for i=1:N, j=1:M))

## Constraints
@constraint(model, MassBal1[t in 2:T], V1[t] == V1[t-1] + q1[t] - u1[t] - s1[t])
@constraint(model, Release1[t in 2:T], min_ut1 <= u1[t] <= max_ut1)
@constraint(model, RampRate1[t in 2:T], RR_dn1 <= u1[t] - u1[t-1] <= RR_up1)
@constraint(model, FeederCap1[t in 1:T], 0 <= p1[t] <= F1)
@constraint(model, Volume1[t in 1:T], min_V1 <= V1[t] <= max_V1)
@constraint(model, lambda1[t in 1:T], sum(lam1[:,:,t]) == 1)

### Unit 02: Dalles Dam (Downstream)
## Define variables
@variable(model, s2[1:T] >= 0) # Spill Outflow [m3/hr]
@variable(model, lam2[1:N, 1:M, 1:T] >= 0)  # lambda matrix

## Initial conditions
@constraint(model, MassBalInit2, sum(lam2[i,j,1] * V2_sample[j] for i=1:N, j=1:M) == V0_02)
@constraint(model, RampRateInit2, sum(lam2[i,j,1] * u2_sample[i] for i=1:N, j=1:M) == min_ut2)

## Expressions
@expression(model, V2[t=1:T], sum(lam2[i,j,t] * V2_sample[j] for i=1:N, j=1:M))
@expression(model, u2[t=1:T], sum(lam2[i,j,t] * u2_sample[i] for i=1:N, j=1:M))
@expression(model, p2[t=1:T], sum(lam2[i,j,t] * P2[i,j]      for i=1:N, j=1:M))

## Constraints
@constraint(model, MassBal2[t in 2:T], V2[t] == V2[t-1] + q2[t] - u2[t] - s2[t])
@constraint(model, Release2[t in 2:T], min_ut2 <= u2[t] <= max_ut2)
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
println(obj)

# -----------------  MODEL REPORT  ----------------- #

# Check the solution method used
method_used = MOI.get(model, Gurobi.ModelAttribute("ConcurrentWinMethod"))
status = MOI.get(model, MOI.TerminationStatus())

# Check iteration counts
iter_count = MOI.get(model, Gurobi.ModelAttribute("IterCount"))
bar_iter_count = MOI.get(model, Gurobi.ModelAttribute("BarIterCount"))

# Check if there are any general constraints (which could make it behave like MIP)
num_gen_constrs = MOI.get(model, Gurobi.ModelAttribute("NumGenConstrs"))
num_sos = MOI.get(model, Gurobi.ModelAttribute("NumSOS"))

println("Method used: ", method_used)  # 0=primal simplex, 1=dual simplex, 2=barrier
println("Status: ", status)
println("Simplex iterations: ", iter_count)
println("Barrier iterations: ", bar_iter_count)
println("General constraints: ", num_gen_constrs)
println("SOS constraints: ", num_sos)

# -----------------  PLOTS  ----------------- #

# Create directory for this run 
stamp = Dates.format(now(), "mm-dd-yyyy HH.MM.SS") * " CONVEX HULL";
dir = "./plots/" ;
path = dir * stamp;
mkdir(path)

p1_max =  (eta * g * rho_w * value.(u1) * a1 .* (value.(V1).^b1))/(3.6e9)
head1 = a1 .* (value.(V1).^b1)
sim_plots(path, "Unit1", T, value.(u1), value.(s1), value.(p1), value.(V1), q1, head1, F1, p1_max, min_ut1, max_ut1, min_h1, max_h1)

p2_max =  (eta * g * rho_w * value.(u2) * a2 .* (value.(V2).^b2))/(3.6e9)
head2 = a2 .* (value.(V2).^b2)
sim_plots(path, "Unit2", T, value.(u2), value.(s2), value.(p2), value.(V2), q2, head2, F1, p2_max, min_ut2, max_ut2, min_h2, max_h2)


println("--- SIMULATION COMPLETE ---")
