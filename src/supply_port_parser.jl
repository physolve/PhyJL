using CSV
using DataFrames
using Plots

file="SupplyPort_2_5.txt"
path = joinpath(@__DIR__, "data", "resultSupply", file);
result_port_raw = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true))
plt2 = plot(result_port_raw.Elapsed, result_port_raw.Pressure, seriestype=:scatter, label="data")