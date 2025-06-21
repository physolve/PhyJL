using DelimitedFiles
using CSV
using DataFrames
using Plots
using Glob

# calculate rate with current value of volume
v_B = 55.5544 # cm3
v_C2 = 51.689 # cm3
t_B = 273+26
n_p = 1 #bar
n_t = 273
prev_time = 0

M = 2.016e-3  # kg/mol

Cv = 0.0005 * 2.2 * 1.7e-5 # Convert Cv to SI units (m³/s·Pa^0.5)

gamma = 1.41
R = 8.314
T = 298

dir = "resultSupply"
dir_glob = joinpath(@__DIR__, "data", dir);
files_AR0 = glob("2025-06-*AR0*.txt",dir_glob)
files_AR1 = glob("2025-06-*AR1*.txt",dir_glob)
files_AR0_model = glob("2025-06-*AR0*model.txt",dir_glob)
files_AR1_model = glob("2025-06-*AR1*model.txt",dir_glob)
files_real_AR0 = setdiff(files_AR0, files_AR0_model)
files_real_AR1 = setdiff(files_AR1, files_AR1_model)

fast_AR1 = glob("2025-06-*AR1*.csv",dir_glob)
# path = joinpath(@__DIR__, "data", "resultSupply", file);

os_dependent_slash = '/'
if Sys.iswindows()
    os_dependent_slash = raw"\\"
end
# parsing TEST WITH LOW PRESSURE 4 bar and high pressure 20 bar 
compare_plots = []

files_AR1 = zip(files_real_AR1,files_AR1_model,fast_AR1)

for (file_real, file_model, file_fast) in files_AR1
    local path = file_real
    local plot_title = split(path, os_dependent_slash)[end]
    local port_pres
    open(path,"r") do f
        s = readline(f)
        infos = split(s,"\t")
        port_pres = parse(Float64, split(infos[2]," ")[3])
    end
    local result_port_real = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true)) 
    path = file_model
    plot_title = split(path, os_dependent_slash)[end]
    port_pres
    open(path,"r") do f
        s = readline(f)
        infos = split(s,"\t")
        port_pres = parse(Float64, split(infos[2]," ")[3])
    end
    local result_port_model = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true)) 
    # plot!(result_port_raw.Elapsed, result_port_raw.Pressure, label = "model")
    push!(compare_plots, plot([result_port_real.Elapsed,result_port_model.Elapsed], [result_port_real.Pressure,result_port_model.Pressure], label = ["real" "model"], title = "GRAM"))
    # push!(compare_plots, plot([result_port_real.Elapsed,result_port_model.Elapsed], [result_port_real.Rate,result_port_model.Rate], label = ["real" "model"]))

    local fast_rate = []
    local fast_time = []
    local model_rate = []
    path = file_fast
    local data, header = readdlm(path, ',', header=true);
    local resultSupply_raw = DataFrame(data, vec(header))
    local x_sec = resultSupply_raw[:,1]
    local rate = 0

    local last_time = resultSupply_raw.time[1]
    local last_pressure = resultSupply_raw.value[1]

    local pressure_inlet = resultSupply_raw.value[end] # rewrite for non end variant

    local v_S = result_port_real.Volume[end]

    for i in eachindex(resultSupply_raw.time)
        # calculate real rate using dP/dt
        if(i==1)
            push!(fast_rate, 0)
            push!(fast_time, 0)
            continue
        end
        if(resultSupply_raw.time[i]-last_time < 0.2)
            continue
        end
        rate = (resultSupply_raw.value[i]-last_pressure)*v_S/(R*T*10)/(resultSupply_raw.time[i]-last_time)
        push!(fast_rate, rate)
        push!(fast_time, last_time)
        last_time = resultSupply_raw.time[i]
        last_pressure = resultSupply_raw.value[i]
    end
    local model_pressure = []
    local fast_pressure = resultSupply_raw.value[1]
    for i in eachindex(resultSupply_raw.time)
        if(i==1)
            push!(model_pressure, fast_pressure)
            push!(model_rate, 0)
            continue
        end
        if(fast_pressure > 0.528*pressure_inlet)
            if(abs(fast_pressure - pressure_inlet) < 0.1)
                rate = 0
            else
                rate = 0.3*sqrt(pressure_inlet-fast_pressure)*Cv*pressure_inlet*1.01*1e5*sqrt(2*gamma/(R/M*(gamma-1)*T))*sqrt((fast_pressure/pressure_inlet)^(2/gamma)-(fast_pressure/pressure_inlet)^((gamma+1)/gamma))
            end
        else
            rate = Cv*pressure_inlet*1.01*1e5*sqrt(gamma/(R/M*T))*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)))
        end
        push!(model_rate, rate/M)
        fast_pressure += rate*(resultSupply_raw.time[i]-resultSupply_raw.time[i-1])/M*R*T/v_S*10
        push!(model_pressure, fast_pressure)
    end
    push!(compare_plots, plot(resultSupply_raw.time, [resultSupply_raw.value, model_pressure], label = ["real" "model"], title = "Julia"))
    # push!(compare_plots, plot([fast_time,resultSupply_raw.time], [fast_rate, model_rate], label = ["real" "model"]))
end

# xlabel!("Время, с")
# title = ["($i)" for j in 1:1, i in 1:length(compare_plots)]
plot(compare_plots..., layout = length(compare_plots), framestyle = :box, size=(1000,1000))

# write reverse calc for pressure_inlet
