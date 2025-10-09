using MAT
using Printf
using CSV
using DataFrames
using Dates
using Base.Filesystem
using LaTeXStrings

function cfs_to_m3s(cfs)
    # Conversion factor: 1 cfs is approximately 0.0283168 m³/s
    conversion_factor = 0.0283168
    return cfs * conversion_factor
end

function min_max_normalize(column)
    col_min = minimum(column)
    col_max = maximum(column)
    return (column .- col_min) ./ (col_max - col_min)
end

function inverse_minmax(norm_values, min_val, max_val)
    return norm_values .* (max_val - min_val) .+ min_val
end

function af_to_m3(acft)
    # Conversion factor: 1 acre-foot is approximately 1233.48 cubic meters
    conversion_factor = 1233.48
    return acft * conversion_factor
end

function ft_to_m(ft)
    # Conversion factor: Divide the length value by 3.281
    conversion_factor = 3.281
    return ft / conversion_factor
end

function cfs_to_af(ft)
    # Conversion factor: 1 CFS sustained for 1 day is 1.9835 acre-feet
    conversion_factor = 1.9835
    return ft / conversion_factor
end


function fullsim_dataload()

    # LMP_path = string("");
    flow_path = string("/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/flowdata.csv");
    bon_storage_path = string("/Users/elizacohn/Desktop/cascaded-hydro/streamflow-data/bonneville/BON6S_daily.csv");
    tda_storage_path = string("/Users/elizacohn/Desktop/cascaded-hydro/streamflow-data/dalles/TDA6S_daily.csv");

    ## Data Set 1: Local Marginal Price
    # RTP = DataFrame(CSV.File(LMP_path)); 

    ## Data Set 2: Flow Data
    flow = DataFrame(CSV.File(flow_path));
    select!(flow, Not(:Column1))
    # (To Do): convert ft to cfs & convert to m3/s from cfs
    flow.down_inflow_m = ft_to_m(flow.down_inflow);
    flow.up_outflow_m = ft_to_m(flow.up_outflow);
    flow.datetime = DateTime.(flow.datetime, dateformat"yyyy-mm-dd HH:MM:SS")

    ## Dataset 3: Storage Levels
    bon = DataFrame(CSV.File(bon_storage_path));
    rename!(bon, Symbol("S (unit:cfs)") => :S_cfs);
    bon.S_cumcfs = cumsum(bon.S_cfs)
    bon.S_m3 = af_to_m3(cfs_to_af(bon.S_cumcfs))

    tda = DataFrame(CSV.File(tda_storage_path));
    rename!(tda, Symbol("S (unit:cfs)") => :S_cfs);
    tda.S_cumcfs = cumsum(tda.S_cfs)
    tda.S_m3 = af_to_m3(cfs_to_af(tda.S_cumcfs))

    return flow, bon, tda
end