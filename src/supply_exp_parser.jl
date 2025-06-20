using DelimitedFiles
using DataFrames
using Plots
file="fastResult_3.csv"
dir = "resultSupply"
path = joinpath(@__DIR__, "data", dir, file);
data, header = readdlm(path, ',', header=true);
resultSupply_raw = DataFrame(data, vec(header))
x_sec = resultSupply_raw[:,1]./1280
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

plt2 = plot(plt1,result_port_raw.Elapsed,result_port_raw.Pressure)



open(path,"r") do f
    s = readline(f)
    infos = split(s,"\t")
    global pressure_inlet = parse(Float64, split(infos[3]," ")[3])
end

result_port_raw = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true))

# пересчет в г/с




prev_vol = (result_port_raw.Pressure[1]*(v_B)/t_B)*n_t/n_p #+v_C2
rate_B = []
real_vol_B = []

model_flow_rate = []
model_vol_B = []
prev_model_vol_B = prev_vol
dp = 0
dt = 0
flow_coefficient = 0.0005
dp_test_plot = []
initial_rate = 0
prev_model_vol_B += initial_rate
for i in eachindex(result_port_raw.Elapsed)
    if(i==length(result_port_raw.Elapsed))
        push!(real_vol_B, real_vol_B[end])
        push!(dp_test_plot, dp_test_plot[end])
        push!(rate_B, rate_B[end])
        push!(model_flow_rate, model_flow_rate[end])
        push!(model_vol_B, model_vol_B[end])
        break
    end
    dp = pressure_inlet-result_port_raw.Pressure[i+1]
    if(dp < 0.33) dp = 0 end
    push!(dp_test_plot, dp)
    dt = result_port_raw.Elapsed[i+1]-prev_time
    cur_vol = (result_port_raw.Pressure[i+1]*(v_B)/t_B)*n_t/n_p #+v_C2
    if(i > 1)
        push!(real_vol_B, cur_vol)
        push!(rate_B, (cur_vol-prev_vol)/(dt))
        if(result_port_raw.Pressure[i+1] < 1/2*pressure_inlet)
            push!(model_flow_rate, 0.471*6950*flow_coefficient*pressure_inlet*sqrt(1/(0.07*t_B))*16.6667)
        else
            push!(model_flow_rate, 6950*flow_coefficient*pressure_inlet*(1-2*dp/(3*pressure_inlet))*sqrt(dp/(pressure_inlet*0.07*t_B))*16.6667)
        end
    else
        push!(real_vol_B, prev_vol)
        push!(rate_B, 0)
        push!(model_flow_rate, 0)
    end
    global prev_model_vol_B += last(model_flow_rate)*dt
    push!(model_vol_B, prev_model_vol_B)
    global prev_vol = cur_vol
    global prev_time = result_port_raw.Elapsed[i+1]
end

theme(:wong2)

# plot!(twinx(),result_port_raw.Elapsed, dp_test_plot,legend=:topleft)
# plot_orig = plot(result_port_raw.Elapsed,[real_vol_B,model_vol_B,result_port_raw."Modelled cm3 H2"], label = ["H2 sensor in B" "H2 model julia" "H2 model GRAM"])
plot_orig = plot(result_port_raw.Elapsed,[real_vol_B,model_vol_B];
    lw=2, label = ["Датчик (B+C2)" "Модель Julia"],framestyle = :box,
)
xlabel!("Время, с")
ylabel!("Объем \$H_2,\\ см^3 \$")
# plot value of deviation in bar

file_model="SupplyPort_4_28_0.txt"
path_model = joinpath(@__DIR__, "data", "resultSupply", file_model);
result_port_model = DataFrame(CSV.File(path_model, delim = '\t', header=2, skipto=3, silencewarnings = true))
plt4 = plot(result_port_raw.Elapsed,result_port_raw.Pressure;
    label = "Давление (B+C2)", lw=2, framestyle = :box,
)
plot_presure = plot(plt4,result_port_model.Elapsed,result_port_model.Pressure;
    label = "Модель GRAM", lw=2, legend=:topleft,
)
xlabel!("Время, с")
ylabel!("Давление, бар")
title!("Расчет давления")
xlims!(0,10)
xticks!([0:1:10;])
vline!([7.78],primary=false,lw=2,linecolor=:blue,linestyle = :dash)
annotate!(7.9, 3, text("7.78 с", 10, :left))

plt6 = plot(result_port_raw.Elapsed,[rate_B,model_flow_rate];
    lw=2, label = ["Датчик (B+C2)" "Модель Julia"], framestyle = :box, legend=:bottomleft,
)
plot(plt6, result_port_model.Elapsed,result_port_model."Flow rate"; 
    lw=2, label = "Пре-расчет GRAM",
)
xlabel!("Время, с")
ylabel!("Поток \$H_2,\\ см^3/с \$")
xlims!(0,10)
xticks!([0:1:10;])
vline!([7.78],primary=false,lw=2,linecolor=:blue,linestyle = :dash)
annotate!(8, 140, text("7.78 с", 10, :left))

plot_model = plot(plot_orig, result_port_model.Elapsed,result_port_model."Modelled cm3 H2", label = "Пре-расчет GRAM", color = 4)
title!("Сравнение расчетов")
plot!(legend=:top, legendcolumns=3)
xlims!(0,8)
xticks!([0:1:8;])
