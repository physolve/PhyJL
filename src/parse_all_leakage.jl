# parse all adata
# R1 case

using DelimitedFiles
using CSV
using DataFrames
using Plots
using Glob

# parameters
v_B = 55.5544 # cm3
v_C2 = 51.689 # cm3
t_B = 273+26
n_p = 1.01 #bar
n_t = 273

v_E = 25.7941
v_D2 = 13.4308
t_E = 273+26

# Constants
gamma = 1.41
M = 2.016e-3  # kg/mol
R = 8.314
T = 298
V1 = 55.5544e-6  # m³
V2 = 25.79e-6    # m³

# needle dependent
Cv_R0 = 0.00085 * 0.8 * 1.7e-5  # Convert Cv to SI units (m³/s·Pa^0.5)
Cv_R1 = 0.00017 * 1.7e-5  # Convert Cv to SI units (m³/s·Pa^0.5)
# pack
dir = "resultLeakage"
dir_glob = joinpath(@__DIR__, "data", dir);
files_R0 = glob("2025-04-*R0*.txt",dir_glob)
files_R1 = glob("2025-04-*R1*.txt",dir_glob)
files_R0_model = glob("2025-04-*R0*model.txt",dir_glob)
files_R1_model = glob("2025-04-*R1*model.txt",dir_glob)
files_real_R0 = setdiff(files_R0, files_R0_model)
files_real_R1 = setdiff(files_R1, files_R1_model)

all_plots_R1 = []

os_dependent_slash = '/'
if Sys.iswindows()
    os_dependent_slash = raw"\\"
end
# parsing TEST WITH LOW PRESSURE 4 bar and high pressure 20 bar 
for file in files_real_R1
    local path = file
    local pressure_storage = 50
    local pressure_reaction = 0

    local plot_title = split(path, os_dependent_slash)[end]
    open(path,"r") do f
        s = readline(f)
        infos = split(s,"\t")
        pressure_storage = parse(Float64, split(infos[3]," ")[3])
        pressure_reaction = parse(Float64, split(infos[4]," ")[3])
    end
    local result_port_raw = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true))
    
    # Initial conditions
    local Cv = Cv_R1
    local P1 = pressure_storage*1e5  # Pa
    local P2 = pressure_reaction*1e5  # Pa
    local t_final = result_port_raw.Elapsed[end]  # s ? from file last 
    local dt = 0.01  # Time step
    # Store results
    local time = []
    local p1 = []
    local p2 = []
    local flow_rate = []
    local ch_cf_R1 = 1.35*log10(pressure_storage-pressure_reaction)
    local subs_cf_R1 = 0.75*log10(pressure_storage-pressure_reaction)
    for t in 0:dt:t_final
        local delta_P = P1 - P2
        if delta_P > 0
            # Flow from R1 to R2
            Pup, Pdown = P1, P2
            # Check if flow is choked
            if Pdown / Pup < (2 / (gamma + 1))^(gamma / (gamma - 1))
                # Choked flow 
                m_dot = Cv * ch_cf_R1 * Pup * sqrt((gamma * M) / (R * T) * (2 / (gamma + 1)) ^ ((gamma + 1) / (gamma - 1)))
            else
                # Subsonic flow
                term = (Pdown / Pup) ^ (2 / gamma) - (Pdown / Pup) ^ ((gamma + 1) / gamma)
                m_dot = Cv * Pup * subs_cf_R1 * sqrt((2 * gamma * M) / ((gamma - 1) * R * T) * term)
            end
        else
            # No flow (delta_P <= 0)
            m_dot = 0
        end
        # Update masses and pressures
        n1 = (P1 * V1) / (R * T)
        n2 = (P2 * V2) / (R * T)
        
        dn = m_dot * dt / M  # Moles transferred
        n1 -= dn
        n2 += dn
        
        P1 = n1 * R * T / V1
        P2 = n2 * R * T / V2
        
        # Store results
        push!(time,t)
        push!(p1,(P1 / 1e5))  # Convert to bar
        push!(p2,(P2 / 1e5))
        push!(flow_rate,(m_dot * 1e3))  # mg/s ???
    end
    real_m_dot_rate_R = []
    local current_mole = 0
    local current_time = 0
    local residual_pressure = [] # Elapsed less than model!
    for i in eachindex(result_port_raw.Elapsed)
        if(i==1)
            current_mole = result_port_raw.Reaction[1]*1e5*V2/(R * T)
             current_time = result_port_raw.Elapsed[1]
            push!(real_m_dot_rate_R,0)
        else
            dNu = result_port_raw.Reaction[i]*1e5*V2/(R * T)-current_mole
            local dTime = result_port_raw.Elapsed[i]-current_time 
            if(dTime > 1e-3)
                push!(real_m_dot_rate_R, dNu/dTime*M*1e3)
                current_mole = result_port_raw.Reaction[i]*1e5*V2/(R * T)
                current_time = result_port_raw.Elapsed[i]
            else
                push!(real_m_dot_rate_R, real_m_dot_rate_R[end]) 
            end
        end
            # residual for reaction 
        local index_close = findmin(x->abs(x-result_port_raw.Elapsed[i]), time)[2]
        if !isnothing(index_close)
            push!(residual_pressure,result_port_raw.Reaction[i] - p2[index_close])
        else
            push!(residual_pressure,0)
        end
    end
    local plot_pressure_real = plot(result_port_raw.Elapsed,[result_port_raw.Storage,result_port_raw.Reaction], label = ["Pressure S real" "Pressure R real"],
    title = "Pressure of $(plot_title)", titlefont=font(10))
    local pressure_real_and_model = plot(plot_pressure_real,time,[p1,p2],
    label = ["Pressure S model Julia" "Pressure R model Julia"], legend=:bottomright)
    local residual_plot = plot(result_port_raw.Elapsed,residual_pressure, linecolor=:red, legend=false)
    hline!([0],primary=false,lw=1, linecolor=:black, linestyle = :dash)
    ylims!(-1, 1)
    push!(all_plots_R1,plot(pressure_real_and_model, residual_plot, layout = (2, 1)))
end

plot(all_plots_R1..., layout = length(all_plots_R1), size=(1200,1200))


#FOR R0#
all_plots_R0 = []
for file in files_real_R0
    local path = file
    local pressure_storage = 50
    local pressure_reaction = 0
    local plot_title = split(file,os_dependent_slash)[end]
    open(path,"r") do f
        s = readline(f)
        infos = split(s,"\t")
        pressure_storage = parse(Float64, split(infos[3]," ")[3])
        pressure_reaction = parse(Float64, split(infos[4]," ")[3])
    end
    local result_port_raw = DataFrame(CSV.File(path, delim = '\t', header=2, skipto=3, silencewarnings = true))
    
    # Initial conditions
    local Cv = Cv_R0
    local P1 = pressure_storage*1e5  # Pa
    local P2 = pressure_reaction*1e5  # Pa
    local t_final = result_port_raw.Elapsed[end]  # s ? from file last 
    local dt = 0.01  # Time step
    # Store results
    local time = []
    local p1 = []
    local p2 = []
    local flow_rate = []
    local ch_cf_R0 = 1.35*log10(pressure_storage-pressure_reaction)
    local subs_cf_R0 = 0.75*log10(pressure_storage-pressure_reaction)
    for t in 0:dt:t_final
        local delta_P = P1 - P2
        if delta_P > 0
            # Flow from R1 to R2
            Pup, Pdown = P1, P2
            # Check if flow is choked
            if Pdown / Pup < (2 / (gamma + 1))^(gamma / (gamma - 1))
                # Choked flow
                m_dot = Cv * Pup * ch_cf_R0 * sqrt((gamma * M) / (R * T) * (2 / (gamma + 1)) ^ ((gamma + 1) / (gamma - 1)))
            else
                # Subsonic flow
                # dp_coef = 0.00455*delta_P*1e-5+0.66451
                term = (Pdown / Pup) ^ (2 / gamma) - (Pdown / Pup) ^ ((gamma + 1) / gamma)
                m_dot = Cv * Pup * subs_cf_R0 * sqrt((2 * gamma * M) / ((gamma - 1) * R * T) * term)
            end
        else
            # No flow (delta_P <= 0)
            m_dot = 0
        end
        # Update masses and pressures
        n1 = (P1 * V1) / (R * T)
        n2 = (P2 * V2) / (R * T)
        
        dn = m_dot * dt / M  # Moles transferred
        n1 -= dn
        n2 += dn
        
        P1 = n1 * R * T / V1
        P2 = n2 * R * T / V2
        
        # Store results
        push!(time,t)
        push!(p1,(P1 / 1e5))  # Convert to bar
        push!(p2,(P2 / 1e5))
        push!(flow_rate,(m_dot * 1e3))  # mg/s ???
    end
    real_m_dot_rate_R = []
    local current_mole = 0
    local current_time = 0
    local residual_pressure = [] # Elapsed less than model!
    for i in eachindex(result_port_raw.Elapsed)
        if(i==1)
            current_mole = result_port_raw.Reaction[1]*1e5*V2/(R * T)
             current_time = result_port_raw.Elapsed[1]
            push!(real_m_dot_rate_R,0)
        else
            dNu = result_port_raw.Reaction[i]*1e5*V2/(R * T)-current_mole
            local dTime = result_port_raw.Elapsed[i]-current_time 
            if(dTime > 1e-3)
                push!(real_m_dot_rate_R, dNu/dTime*M*1e3)
                current_mole = result_port_raw.Reaction[i]*1e5*V2/(R * T)
                current_time = result_port_raw.Elapsed[i]
            else
                push!(real_m_dot_rate_R, real_m_dot_rate_R[end]) 
            end
        end
            # residual for reaction 
        local index_close = findmin(x->abs(x-result_port_raw.Elapsed[i]), time)[2]
        if !isnothing(index_close)
            push!(residual_pressure,result_port_raw.Reaction[i] - p2[index_close])
        else
            push!(residual_pressure,0)
        end
    end
    local plot_pressure_real = plot(result_port_raw.Elapsed,[result_port_raw.Storage,result_port_raw.Reaction], label = ["Pressure S real" "Pressure R real"],
    title = "Pressure of $(plot_title)", titlefont=font(10))
    local pressure_real_and_model = plot(plot_pressure_real,time,[p1,p2],
    label = ["Pressure S model Julia" "Pressure R model Julia"], legend=:bottomright)
    local residual_plot = plot(result_port_raw.Elapsed,residual_pressure, linecolor=:red, legend=false)
    hline!([0],primary=false,lw=1, linecolor=:black, linestyle = :dash)
    ylims!(-1, 1)
    push!(all_plots_R0,plot(pressure_real_and_model, residual_plot, layout = (2, 1)))
end

plot(all_plots_R0..., layout = length(all_plots_R0), size=(1200,1200))

all_plots = []
append!(all_plots, all_plots_R0, all_plots_R1)
plot(all_plots..., layout = length(all_plots), size=(1200,1200))