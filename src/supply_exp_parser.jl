using DelimitedFiles
using DataFrames
using Plots
file="2025-04-09_fastResult_R0_6.csv"
dir = "resultLeakage"
path = joinpath(@__DIR__, "data", dir, file);
data, header = readdlm(path, ',', header=true);
resultSupply_raw = DataFrame(data, vec(header))
x_sec= resultSupply_raw[:,1]./1024
resultSupply_raw[:,1] = x_sec
plt1 = plot(resultSupply_raw.time, resultSupply_raw.value)
xlims!(8,12.0)
