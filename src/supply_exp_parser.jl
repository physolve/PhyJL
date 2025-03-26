using DelimitedFiles
using DataFrames
using Plots
file="fastResult_26_4.csv"
path = joinpath(@__DIR__, "data", "resultSupply", file);
data, header = readdlm(path, ',', header=true);
resultSupply_raw = DataFrame(data, vec(header))
x_sec= resultSupply_raw[:,1]./1024
resultSupply_raw[:,1] = x_sec
plt1 = plot(resultSupply_raw.time, resultSupply_raw.value)
xlims!(10.0,19.0)
