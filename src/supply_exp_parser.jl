using DelimitedFiles
using DataFrames
using Plots
file="2025-04-10_fastResult_R1_6.csv"
dir = "resultLeakage"
path = joinpath(@__DIR__, "data", dir, file);
data, header = readdlm(path, ',', header=true);
resultSupply_raw = DataFrame(data, vec(header))
x_sec = resultSupply_raw[:,1]./1024
resultSupply_raw[:,1] = x_sec

# total_seconds = 3

# index_first = findfirst(==(6), x_sec)
# index_last = findfirst(==(9), x_sec)
# total_points = index_last - index_first
# points_per_second = Int(total_points/total_seconds)
# theme(:wong2)
# plt1 = plot(
#     resultSupply_raw.time, resultSupply_raw.value;
#     lw=2, ls=:dot, label="Давление в SQ", framestyle = :box,
# )
# xlims!(5.5,9.5)
# # annotate!(6, 0, text("\$\\bar \\delta \$", :left))
# vline!([6],label = "Открыт AR2")
# vline!([9],label = "Закрыт AR2")
# title!("Подача водорода в SQ из генератора")
# xlabel!("Время, с")
# ylabel!("Давление, бар")
# quiver!([9], [0.4], quiver=([-3], [0]), linecolor=:darkolivegreen)
# quiver!([6], [0.4], quiver=([3], [0]), linecolor=:darkolivegreen)
# annotate!([7.5], [0.8], ("$points_per_second точек в секунду", 12))

plt1 = plot(resultSupply_raw.time, resultSupply_raw.value)
# xlims!(8,12.0)
