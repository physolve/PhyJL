using DelimitedFiles
using CSV
using DataFrames
using Plots
runT = 5
file="2025-04-10_R1_$(runT).txt"
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

# Constants
Cv = 0.00017 * 1.7e-5  # Convert Cv to SI units (m³/s·Pa^0.5)
gamma = 1.41
M = 2.016e-3  # kg/mol
R = 8.314
T = 298
V1 = 55.5544e-6  # m³
V2 = 25.79e-6    # m³

# Initial conditions
P1 = result_port_raw.Storage[1]*1e5  # Pa
P2 = result_port_raw.Reaction[1]*1e5  # Pa
t_final = result_port_raw.Elapsed[end]  # s ? from file last 
dt = 0.01  # Time step

# Store results
time = []
p1 = []
p2 = []
flow_rate = []
delta_P = P1 - P2
dp_coef_ch = 0.0375*delta_P*1e-5+1
(2 / (gamma + 1))^(gamma / (gamma - 1))*result_port_raw.Storage[1] # pressure of P2 for subsonic
dp_coef_subs = 1 #0.016875*delta_P*1e-5+0.525

for t in 0:dt:t_final
    global delta_P = P1 - P2
    if delta_P > 0
        # Flow from R1 to R2
        Pup, Pdown = P1, P2
        # Check if flow is choked
        if Pdown / Pup < (2 / (gamma + 1))^(gamma / (gamma - 1))
            # Choked flow
            m_dot = Cv * dp_coef_ch * Pup * sqrt((gamma * M) / (R * T) * (2 / (gamma + 1)) ^ ((gamma + 1) / (gamma - 1)))
        else
            # Subsonic flow
            # dp_coef = 0.00455*delta_P*1e-5+0.66451
            term = (Pdown / Pup) ^ (2 / gamma) - (Pdown / Pup) ^ ((gamma + 1) / gamma)
            m_dot = Cv * dp_coef_subs * Pup * sqrt((2 * gamma * M) / ((gamma - 1) * R * T) * term)
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
    
    global P1 = n1 * R * T / V1
    global P2 = n2 * R * T / V2
    
    # Store results
    push!(time,t)
    push!(p1,(P1 / 1e5))  # Convert to bar
    push!(p2,(P2 / 1e5))
    push!(flow_rate,(m_dot * 1e3))  # mg/s ???
end

real_m_dot_rate_R = []
current_mole = 0
current_time = 0
for i in eachindex(result_port_raw.Elapsed)
    if(i==1)
        global current_mole = result_port_raw.Reaction[1]*1e5*V2/(R * T)
        global current_time = result_port_raw.Elapsed[1]
        push!(real_m_dot_rate_R,0)
    else
        dNu = result_port_raw.Reaction[i]*1e5*V2/(R * T)-current_mole
        dt = result_port_raw.Elapsed[i]-current_time 
        if(dt > 1e-3)
            push!(real_m_dot_rate_R, dNu/dt*M*1e3)
            global current_mole = result_port_raw.Reaction[i]*1e5*V2/(R * T)
            global current_time = result_port_raw.Elapsed[i]
        else
            push!(real_m_dot_rate_R, real_m_dot_rate_R[end]) 
        end
    end
end
plot_pressure_real = plot(result_port_raw.Elapsed,[result_port_raw.Storage,result_port_raw.Reaction], label = ["Pressure S real" "Pressure R real"],
title = "Pressure of $file", titlefont=font(12))
pressure_real_and_model = plot(plot_pressure_real,time,[p1,p2], label = ["Pressure S model Julia" "Pressure R model Julia"])
plot_rate_real = plot(result_port_raw.Elapsed, real_m_dot_rate_R, label = "Rate real",
title = "Rate of $file", titlefont=font(12))
rate_real_and_model = plot(plot_rate_real, time, flow_rate, label = "Rate model Julia")
plot(pressure_real_and_model, rate_real_and_model)

file_model="2025-04-10_R1_$(runT)_model.txt"
path_model = joinpath(@__DIR__, "data", dir, file_model);
result_port_model = DataFrame(CSV.File(path_model, delim = '\t', header=2, skipto=3, silencewarnings = true))
pressure_model = plot(pressure_real_and_model,result_port_model.Elapsed,[result_port_model.Storage,result_port_model.Reaction], 
label = ["Pressure S model GRAM" "Pressure R model GRAM"])
rate_model = plot(rate_real_and_model,result_port_model.Elapsed, result_port_model."Flow rate", label = "Rate model GRAM")
plot(pressure_model, rate_model)

# plot(plt6, result_port_model.Elapsed,result_port_model."Flow rate", label = "H2 rate cm3/s pre model GRAM")
# plot_model = plot(plot_orig, result_port_model.Elapsed,result_port_model."Modelled cm3 H2", label = "H2 pre model GRAM", color = 4)
# plot(rate_model, result_port_raw.Elapsed, result_port_raw."Flow rate", label = "Rate GRAM on spot")