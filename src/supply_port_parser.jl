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

Cv = 0.0005* 1.7e-5  # Convert Cv to SI units (m³/s·Pa^0.5)

gamma = 1.41
R = 8.314
T = 298

pressure_inlet = [32,45.1,53] #[55,55,55]

dir = "resultSupply"
dir_glob = joinpath(@__DIR__, "data", dir);
files_AR0 = glob("2025-04-*AR0*.txt",dir_glob)
files_AR1 = glob("2025-04-*AR1*.txt",dir_glob)
files_AR0_model = glob("2025-04-*AR0*model.txt",dir_glob)
files_AR1_model = glob("2025-04-*AR1*model.txt",dir_glob)
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
plot(0,0)
file_n = 0
for file in fast_AR1
    global file_n += 1
    local fast_rate = []
    # local path = joinpath(@__DIR__, "data", dir, file);
    local path = file
    local data, header = readdlm(path, ',', header=true);
    local resultSupply_raw = DataFrame(data, vec(header))
    local x_sec = resultSupply_raw[:,1]./1280
    resultSupply_raw[:,1] = x_sec
    # local resultSupply_raw = DataFrame(CSV.File(path, delim = ',', header=0, skipto=1, silencewarnings = true))
    # local plot_title = split(path, os_dependent_slash)[end]
    
    local rate = 0
    for i in eachindex(resultSupply_raw.time)
        if(resultSupply_raw.value[i] > 0.528*pressure_inlet[file_n])
            rate = Cv*pressure_inlet[file_n]*1.01*1e5*sqrt(2*gamma/(R/M*(gamma-1)*T))*sqrt((resultSupply_raw.value[i]/pressure_inlet[file_n])^(2/gamma)-(resultSupply_raw.value[i]/pressure_inlet[file_n])^((gamma+1)/gamma))
        else
            rate = Cv*pressure_inlet[file_n]*1.01*1e5*sqrt(gamma/(R/M*T))*(2/(gamma+1))^((gamma+1)/(2*(gamma-1)))
        end
        push!(fast_rate, rate)
    end
    
    # for i in eachindex(resultSupply_raw.time)
    #     if(i>1 && i%1000==0)
    #         local bar_rate = (resultSupply_raw.value[i]-resultSupply_raw.value[i-999])/(resultSupply_raw.time[i]-resultSupply_raw.time[i-999])
    #         rate = bar_rate*v_B*mass_H2/(10*R*T)
    #     end
    #     push!(fast_rate, rate)
    # end
    push!(compare_plots, plot(resultSupply_raw.time, fast_rate))
    # plot!(resultSupply_raw.time, resultSupply_raw.value)
end
# xlims!(0,10)
# xlabel!("Время, с")
# ylabel!("Давление, бар")
# title!("Подача через AR1")
plot(compare_plots..., layout = length(compare_plots), framestyle = :box, size=(1000,1000))