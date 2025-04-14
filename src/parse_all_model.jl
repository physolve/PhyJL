using DelimitedFiles
using CSV
using DataFrames
using Plots
using Glob
# pack
dir = "resultLeakage"
dir_glob = joinpath(@__DIR__, "data", dir);
files = glob("2025-04-*.txt",dir_glob)
files_model = glob("2025-04-*model.txt",dir_glob)
files_real = setdiff(files, files_model)

if(length(files_model)!=length(files_real))
    if length(files_model)<length(files_real)
        println("Missing model files!!!")
    else
        println("Missing real files!!!")
    end
    return
end

files_pair = tuple.(files_real, files_model)
os_dependent_slash = '/'
if Sys.iswindows()
    os_dependent_slash = raw"\\"
end

all_plots = []

for pair in files_pair
    local path = pair[1]
    local plot_title = split(path, os_dependent_slash)[end]
    open(path,"r") do f
        s = readline(f)
        infos = split(s,"\t")
    end
    local result_port_raw = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true))

    local result_model = DataFrame(CSV.File(pair[2], delim = '\t', header=2, skipto=3, silencewarnings = true))

    local residual_pressure = [] # Elapsed less than model!
    for i in eachindex(result_port_raw.Elapsed)
        # residual for reaction 
        local index_close = findmin(x->abs(x-result_port_raw.Elapsed[i]), result_model.Elapsed)[2]
        if !isnothing(index_close)
            push!(residual_pressure,result_port_raw.Reaction[i] - result_model.Reaction[index_close])
        else
            push!(residual_pressure,0)
        end
    end
    
    local plot_pressure_real = plot(result_port_raw.Elapsed,[result_port_raw.Storage,result_port_raw.Reaction], label = ["Pressure S real" "Pressure R real"],
    title = "Pressure of $(plot_title)", titlefont=font(10))
    local pressure_real_and_model = plot(plot_pressure_real,result_model.Elapsed,[result_model.Storage,result_model.Reaction],
    label = ["Pressure S model GRAM" "Pressure R model GRAM"], legend=:bottomright)
    local residual_plot = plot(result_port_raw.Elapsed, residual_pressure, linecolor=:red, legend=false)
    hline!([0],primary=false,lw=1, linecolor=:black, linestyle = :dash)
    ylims!(-1, 1)
    push!(all_plots,plot(pressure_real_and_model, residual_plot, layout = (2, 1)))
end

plot(all_plots..., layout = length(all_plots), size=(1200,1200))
