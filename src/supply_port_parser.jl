using CSV
using DataFrames
using Plots

file="SupplyPort_2_3.txt"
path = joinpath(@__DIR__, "data", "resultSupply", file);
result_port_raw = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true))
plt2 = plot(result_port_raw.Elapsed, result_port_raw.Pressure, seriestype=:scatter, label="data")

result_port_raw."Flow rate".=result_port_raw."Flow rate".*16.6667

# calculate rate with current value of volume
v_B = 55.5544 # cm3
v_C2 = 51.689 # cm3
t_B = 273+26
n_p = 1 #bar
n_t = 273
prev_time = 0
prev_vol = (result_port_raw.Pressure[1]*(v_B+v_C2)/t_B)*n_t/n_p
rate_B = []
real_vol_B = []
pressure_inlet = 32.042
model_flow_rate = []
model_vol_B = []
prev_model_vol_B = prev_vol
dp = 0
dt = 0
flow_coefficient = 0.001225

for i in eachindex(result_port_raw.Elapsed)
    cur_vol = (result_port_raw.Pressure[i]*(v_B+v_C2)/t_B)*n_t/n_p
    push!(real_vol_B, cur_vol)
    dp = pressure_inlet-result_port_raw.Pressure[i]
    dt = result_port_raw.Elapsed[i]-prev_time
    if(i > 2)
        push!(rate_B, (cur_vol-prev_vol)/(dt))
        if(result_port_raw.Pressure[i] < 1/2*pressure_inlet)
            push!(model_flow_rate, 0.471*6950*flow_coefficient*pressure_inlet*sqrt(1/(0.07*t_B))*16.6667)
        else
            push!(model_flow_rate, 6950*flow_coefficient*pressure_inlet*(1-2*dp/(3*pressure_inlet))*sqrt(dp/(pressure_inlet*0.07*t_B))*16.6667)
        end
    else
        push!(rate_B, 0)
        push!(model_flow_rate, 0)
    end
    global prev_model_vol_B += last(model_flow_rate)*dt
    push!(model_vol_B, prev_model_vol_B)
    global prev_vol = cur_vol
    global prev_time = result_port_raw.Elapsed[i]
end
plot(result_port_raw.Elapsed,[model_flow_rate,rate_B,result_port_raw."Flow rate"], label = ["H2 rate cm3/s model" "H2 rate cm3/s sensor" "H2 rate cm3/s model GRAM"])
plot(result_port_raw.Elapsed,[real_vol_B, model_vol_B,result_port_raw."Modelled cm3 H2"], label = ["H2 sensor in B" "H2 model julia" "H2 model GRAM"])

# rate_B.= / result_port_raw.Elapsed.
