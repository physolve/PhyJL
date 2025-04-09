using DelimitedFiles
using CSV
using DataFrames
using Plots

file="2025-04-08_R0_5.txt"
dir = "resultLeakage"
path = joinpath(@__DIR__, "data", dir, file);

pressure_storage = 50
pressure_reaction = 0
open(path,"r") do f
    s = readline(f)
    infos = split(s,"\t")
    global pressure_storage = parse(Float64, split(infos[3]," ")[3])
    global pressure_reaction = parse(Float64, split(infos[4]," ")[3])
end
result_port_raw = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true))

# calculate rate with current value of volume
v_B = 55.5544 # cm3
v_C2 = 51.689 # cm3
t_B = 273+26
n_p = 1.01 #bar
n_t = 273

v_E = 25.7941
v_D2 = 13.4308
t_E = 273+26

prev_time = result_port_raw.Elapsed[1]
prev_vol_S = (result_port_raw.Storage[1]*(v_B)/t_B)*n_t/n_p # +v_C2 optional
rate_S = []
real_vol_S = []

prev_vol_R = (result_port_raw.Reaction[1]*(v_E)/t_E)*n_t/n_p # +v_D2 optional
rate_R = []
real_vol_R = []

model_flow_rate = []
model_vol_R = []
prev_model_vol_R = prev_vol_R
dp = 0
dt = 0

# turn = 7 # from file
# flow_coefficient = 0.00054*turn-0.0014 #
flow_coefficient = 0.00085
dp_test_plot = []

initial_rate = 0
prev_model_vol_R += initial_rate


R = 8.314           # Universal gas constant (J/mol·K)
M_H2 = 2.016e-3     # Molar mass of H2 (kg/mol)
T = 298             # Temperature (K)
R_H2 = R / M_H2     # Specific gas constant (J/kg·K)

for i in eachindex(result_port_raw.Elapsed)
    if(i==length(result_port_raw.Elapsed))
        push!(real_vol_R, real_vol_R[end])
        push!(dp_test_plot, dp_test_plot[end])
        push!(rate_R, rate_R[end])
        push!(model_flow_rate, model_flow_rate[end])
        push!(model_vol_R, model_vol_R[end])
        break
    end
    dp = result_port_raw.Storage[i+1] - result_port_raw.Reaction[i+1]
    if(dp < 0.33) dp = 0 end
    push!(dp_test_plot, dp)
    global dt = result_port_raw.Elapsed[i+1]-prev_time
    cur_vol_R = (result_port_raw.Reaction[i+1]*(v_E)/t_E)*n_t/n_p

    if(i > 1)
        push!(real_vol_R, cur_vol_R)
        push!(rate_R, (cur_vol_R-prev_vol_R)/dt)
        p1 = result_port_raw.Storage[i+1]
        p2 = result_port_raw.Reaction[i+1]
        # if(p2/p1 < 0.528) 
        #     push!(model_flow_rate, 0.471*6950*flow_coefficient*p1*sqrt(1/(0.07*t_B))*16.6667)
        # else
        #     push!(model_flow_rate, 6950*flow_coefficient*p1*(1-2*dp/(3*p1))*sqrt(dp/(p1*0.07*t_B))*16.6667) 
        # end
        #
        # p1 = result_port_raw.Storage[i+1]
        # p2 = result_port_raw.Reaction[i+1] 
        # if(p2/p1 < 0.528) 
        #     push!(model_flow_rate, 152.98*flow_coefficient*p1*sqrt(1/(0.07*t_E))*16.6667)
        # else
        #     push!(model_flow_rate, 306.91*flow_coefficient*sqrt((p1^2-p2^2)/(0.07*t_E))*16.6667)
        # end
        # if(p2/p1 < 0.528) 
        #     push!(model_flow_rate, 0.021*1360*flow_coefficient*p1/sqrt(0.0696)*sqrt(0.53*p1/t_E)*16.6667)
        # else
        #     push!(model_flow_rate, 0.014*1360*flow_coefficient*p1/sqrt(0.0696)*sqrt(dp/t_E)*16.6667)
        # end
        # p_avg = (result_port_raw.Storage[i+1] + result_port_raw.Reaction[i+1])/2
        # push!(model_flow_rate, 1360*flow_coefficient*sqrt((dp*p_avg)/(0.0696*273.15)))
        # p1 = result_port_raw.Storage[i+1]
        # p2 = result_port_raw.Reaction[i+1]
        gamma = 1.41
        M = 2.016e-3
        if(p2/p1 < 0.528) # sonic (chocked)
            push!(model_flow_rate, flow_coefficient*1.7e-5*p1*1e5
            *sqrt(gamma*M/(1*R*t_E)*(2/(gamma+1))^((gamma+1)/(gamma-1)))
            /M*22400/t_E*n_t
            )
        else
            term = (p2/p1)^(2/gamma)-(p2/p1)^((gamma+1)/gamma)
            push!(model_flow_rate, flow_coefficient*1.7e-5*p1*1e5*
            sqrt((2*gamma*M)/((gamma-1)R*t_E)*term)
            /M*22400/t_E*n_t
            )
            
        end
        # p1 = result_port_raw.Storage[i+1]
        # p2 = result_port_raw.Reaction[i+1] 
        # if(p2/p1 < 0.528) 
        #     push!(model_flow_rate, 26.9*flow_coefficient*p1*sqrt(2.016/t_E)*R*t_E*10/(2.016*p1))
        # else
        #     push!(model_flow_rate, 22.67*flow_coefficient*p1*sqrt((2.016/t_E)*(1-(p2/p1)^2))*R*t_E*10/(2.016*p1))
        # end
    else
        push!(real_vol_R, prev_vol_R)
        push!(rate_R, 0)
        push!(model_flow_rate, 0)
    end
    global prev_model_vol_R += last(model_flow_rate)*dt #/t_E*n_t/n_p
    push!(model_vol_R, prev_model_vol_R)
    global prev_vol_R = cur_vol_R
    global prev_time = result_port_raw.Elapsed[i+1]
end
real_vol_R
plt6 = plot(result_port_raw.Elapsed,[model_flow_rate,rate_R], label = ["H2 rate cm3/s model" "H2 rate cm3/s sensor"], title = "Rate of $file", titlefont=font(12))
plt6 = plot(plt6,result_port_raw.Elapsed,result_port_raw."Flow rate", label = "H2 rate cm3/s model GRAM")
# plot!(twinx(),result_port_raw.Elapsed, dp_test_plot,legend=:topleft)
plot_orig = plot(result_port_raw.Elapsed,[real_vol_R,model_vol_R], label = ["H2 sensor in R" "H2 model julia"], title = "Volume of $file", titlefont=font(12))
plot_orig = plot(plot_orig, result_port_raw.Elapsed, result_port_raw."Modelled cm3 H2", label = "H2 model GRAM")
# plot value of deviation in bar

# file_model="2025-04-04_R0_1_model.txt"
# path_model = joinpath(@__DIR__, "data", dir, file_model);
# result_port_model = DataFrame(CSV.File(path_model, delim = '\t', header=2, skipto=3, silencewarnings = true))
# plt4 = plot(result_port_raw.Elapsed,result_port_raw.Reaction)
# plt5 = plot(plt4,result_port_model.Elapsed,result_port_model.Reaction)
# plot(plt6, result_port_model.Elapsed,result_port_model."Flow rate", label = "H2 rate cm3/s pre model GRAM")
# plot_model = plot(plot_orig, result_port_model.Elapsed,result_port_model."Modelled cm3 H2", label = "H2 pre model GRAM", color = 4)
