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
    gage_path = string("/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/USGS-gage-data.csv");
    tda_inflow_path = string("/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/tda-inflow-USACE.csv");
    bon_inflow_path = string("/Users/elizacohn/Desktop/cascaded-hydro/simulation-data/bon-inflow-USACE.csv");
    bon_storage_path = string("/Users/elizacohn/Desktop/cascaded-hydro/streamflow-data/bonneville/BON6S_daily.csv");
    tda_storage_path = string("/Users/elizacohn/Desktop/cascaded-hydro/streamflow-data/dalles/TDA6S_daily.csv");

    ## Data Set 1: Local Marginal Price
    # RTP = DataFrame(CSV.File(LMP_path)); 

    ## Data Set 2: USGS Gauge Data
    gage = DataFrame(CSV.File(gage_path));
    select!(gage, Not(:Column1))
    gage.down_inflow_m = ft_to_m(gage.down_inflow);
    gage.up_outflow_m = ft_to_m(gage.up_outflow);
    gage.datetime = DateTime.(gage.datetime, dateformat"yyyy-mm-dd HH:MM:SS")
    gage = filter(row -> minute(row.datetime) == 0, gage)

    ## Data Set C: USACE Inflow Data
    inflow = DataFrame(CSV.File(tda_inflow_path));
    inflow.datetime = DateTime.(inflow.datetime, dateformat"yyyy-mm-dd HH:MM:SS")
    inflow.tda_inflow_m3s = 1000*cfs_to_m3s(inflow.inflow_kcfs) # convert kcfs --> cfs --> m3/s 
    bon_inflow = DataFrame(CSV.File(bon_inflow_path));
    inflow.bon_inflow_m3s = bon_inflow.inflow_m3s

    ## Dataset 4: BPA Storage Levels
    storage = DataFrame(CSV.File(bon_storage_path));
    rename!(storage, Symbol("S (unit:cfs)") => :bon_S_cfs);
    storage.bon_S_cumcfs = cumsum(storage.bon_S_cfs)
    storage.bon_S_m3 = af_to_m3(cfs_to_af(storage.bon_S_cumcfs))

    tda_S = DataFrame(CSV.File(tda_storage_path));
    rename!(tda_S, Symbol("S (unit:cfs)") => :S_cfs);
    tda_S.S_cumcfs = cumsum(tda_S.S_cfs)
    storage.tda_S_m3 = af_to_m3(cfs_to_af(tda_S.S_cumcfs))

    return flow, inflow, storage
end