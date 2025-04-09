import numpy as np
import matplotlib.pyplot as plt
import time

# -- 1. Constants and Initial Conditions --

# Gas Properties (Hydrogen)
M = 2.016e-3  # Molar mass (kg/mol)
k = 1.41       # Specific heat ratio (dimensionless)

# Thermodynamic Constants
R = 8.314      # Universal gas constant (J/(mol*K))

# System Parameters
V1_cm3 = 55.0
V2_cm3 = 25.0
V1 = V1_cm3 * 1e-6 # Reservoir 1 volume (m^3)
V2 = V2_cm3 * 1e-6 # Reservoir 2 volume (m^3)
Cv = 0.00013     # Valve flow coefficient (dimensionless, standard US units)
T_C = 20.0       # Room temperature (°C)
T = T_C + 273.15 # Temperature (K)

# Initial Conditions
P1_0_bar = 7.0
P2_0_bar = 0.5
P1_0 = P1_0_bar * 1e5 # Initial pressure in Reservoir 1 (Pa)
P2_0 = P2_0_bar * 1e5 # Initial pressure in Reservoir 2 (Pa)

# -- 2. Calculate Derived Constants --

# Critical pressure ratio for choked flow
r_c = (2 / (k + 1))**(k / (k - 1))

# Estimate mass flow constants (needs careful derivation based on Cv standard)
# We use the initial choked flow calculation result (approx. 1.97e-6 kg/s at 7 bar)
# Q_m0 = 1.97e-6 # kg/s (approximate mass flow at t=0)
# C_choked = Q_m0 / P1_0 # kg/(s*Pa)
# This method requires knowing Q_m0 accurately.
# Let's re-calculate Q_m0 using a standard Cv formula interpretation for choked flow.
# Q (SCFH) = 1360 * Cv * P1_psia * sqrt(1 / (G * T_R))
# G_H2 = M / (28.97e-3) = 2.016 / 28.97 approx 0.0696
# T_R = (T_C * 9/5 + 32) + 459.67 # Temperature in Rankine
T_R = T * 9/5 # Temperature in Rankine is T(K) * 1.8
P1_0_psia = P1_0_bar * 14.5038
Q_scfh = 1360 * Cv * P1_0_psia * np.sqrt(1 / (0.0696 * T_R))
# Convert SCFH to kg/s
# Standard conditions for SCFH often 60°F (15.6°C = 288.7K) and 1 atm (101325 Pa)
# Density at standard conditions: rho_std = P_std * M / (R * T_std)
rho_std = 101325 * M / (R * 288.7)
# Conversion: 1 ft^3 = 0.028317 m^3, 1 hr = 3600 s
Q_m0 = Q_scfh * 0.028317 * rho_std / 3600

# Now calculate choked constant C_choked based on this Q_m0 and initial P1_0 (Pa)
C_choked = Q_m0 / P1_0 # kg/(s*Pa)

# Calculate non-choked constant ensuring continuity at r_c
# Assuming Q_m_non_choked = C_non_choked * sqrt(P1^2 - P2^2)
# At transition: C_choked * P1 = C_non_choked * sqrt(P1^2 - (r_c*P1)^2)
# C_non_choked = C_choked / sqrt(1 - r_c**2)
C_non_choked = C_choked / np.sqrt(1 - r_c**2) # kg/(s*Pa)

# Pressure change constants
K1 = (R * T) / (M * V1) # J/(kg*m^3) or m^2/s^2
K2 = (R * T) / (M * V2) # J/(kg*m^3) or m^2/s^2

# -- 3. Simulation Parameters --
dt = 0.1          # Time step (seconds) - adjust if needed
# Estimate time scale - very rough estimate: time ~ V*dP / Q_avg
# Total moles transferred ~ P_final*(V1+V2)/(RT) - P1_0*V1/(RT) - P2_0*V2/(RT)
# Average flow rate is hard to guess, maybe half initial?
# time_est = abs(P1_0-P_final)*V1 / (R*T/M * Q_m0 * 0.5) # Very rough
# Let's run for a fixed long duration first
t_max = 7200       # Maximum simulation time (seconds) - e.g., 2 hours
pressure_tolerance = 100 # Stop when pressure difference is below this (Pa) (e.g. 0.001 bar)

# -- 4. Initialize Data Storage --
times = [0.0]
P1_vals = [P1_0]
P2_vals = [P2_0]
Qm_vals = [Q_m0] # Store initial flow rate

# -- 5. Simulation Loop (Euler Method) --
t_start = time.time()
t = 0.0
P1 = P1_0
P2 = P2_0
Q_m = Q_m0 # Initial flow rate for first calculation step

while t < t_max:
    # Check for equilibrium
    if abs(P1 - P2) < pressure_tolerance:
        print(f"\nEquilibrium reached near t = {t:.2f} s")
        break

    # Calculate current pressure ratio (avoid division by zero)
    if P1 <= P2 or P1 < 1e-3: # Added P1 check for safety
        Q_m = 0.0
        r = 1.0 # Effectively non-choked, flow is zero
    else:
        r = P2 / P1
        # Determine Flow Regime and Calculate Q_m
        if r <= r_c:
            # Choked flow
            Q_m = C_choked * P1
        else:
            # Non-choked flow
            # Use max(0,...) for numerical stability near equilibrium
            pressure_term_sq = max(0, P1**2 - P2**2)
            Q_m = C_non_choked * np.sqrt(pressure_term_sq)

    # Calculate Pressure Changes
    dP1 = -K1 * Q_m * dt
    dP2 = +K2 * Q_m * dt

    # Update Pressures
    P1_new = P1 + dP1
    P2_new = P2 + dP2

    # Prevent pressures from going below zero in case of large dt
    P1 = max(1e-3, P1_new) # Avoid zero or negative pressure
    P2 = max(1e-3, P2_new)

    # Update Time
    t = t + dt

    # Store results
    times.append(t)
    P1_vals.append(P1)
    P2_vals.append(P2)
    Qm_vals.append(Q_m) # Store the flow rate *used* for this step's calculation

# Ensure last calculated flow rate is added if loop terminated by t_max
if t >= t_max:
     print(f"\nSimulation stopped at t_max = {t_max:.2f} s")
# Calculate final Q_m if needed (will be near zero)
# Qm_vals.append(Q_m) # Q_m used in last step is already appended

t_end = time.time()
print(f"Simulation duration: {t_end - t_start:.3f} seconds")

# -- 6. Post-Processing and Output --
P1_bar = np.array(P1_vals) / 1e5
P2_bar = np.array(P2_vals) / 1e5
Qm_kgs = np.array(Qm_vals)

# Convert mass flow to SLPM (Standard Liters Per Minute)
# Standard conditions: 0°C (273.15 K) and 1 atm (101325 Pa) often used for SLPM
rho_std_slpm = 101325 * M / (R * 273.15) # kg/m^3 at STP for H2
# Conversion Q_slpm = Q_m (kg/s) / rho_std (kg/m^3) * (1 m^3 / 1000 L) * (60 s / 1 min)
Q_slpm = Qm_kgs / rho_std_slpm * 60 * 1000

print(f"\nInitial Conditions:")
print(f"  P1 = {P1_0_bar:.2f} bar, P2 = {P2_0_bar:.2f} bar")
print(f"  Initial Flow Rate = {Q_slpm[0]:.4f} SLPM ({Qm_kgs[0]*1e6:.2f} mg/s)")
print(f"\nFinal State (at t = {times[-1]:.2f} s):")
print(f"  P1 = {P1_bar[-1]:.4f} bar, P2 = {P2_bar[-1]:.4f} bar")
# Calculate theoretical P_final
P_final_bar = (P1_0_bar * V1_cm3 + P2_0_bar * V2_cm3) / (V1_cm3 + V2_cm3)
print(f"  Theoretical P_final = {P_final_bar:.4f} bar")
print(f"  Final Flow Rate = {Q_slpm[-1]:.4f} SLPM ({Qm_kgs[-1]*1e6:.2f} mg/s)")
print(f"Critical pressure ratio r_c = {r_c:.4f}")

# Find transition time (approximate)
transition_index = -1
for i in range(len(P1_vals)):
    if P1_vals[i] > 1e-3 and P2_vals[i] / P1_vals[i] > r_c:
        transition_index = i
        break

if transition_index > 0 :
    print(f"Approximate transition from choked to non-choked flow around t = {times[transition_index]:.2f} s")
    print(f"  Conditions at transition: P1={P1_bar[transition_index]:.3f} bar, P2={P2_bar[transition_index]:.3f} bar, Ratio={P2_bar[transition_index]/P1_bar[transition_index]:.3f}")
else:
     print("Flow remained choked or simulation ended before transition.")


# -- 7. Plotting --
plt.figure(figsize=(12, 8))

# Plot Pressures
plt.subplot(2, 1, 1)
plt.plot(times, P1_bar, label=f'P1 (V={V1_cm3} cm³)')
plt.plot(times, P2_bar, label=f'P2 (V={V2_cm3} cm³)')
plt.axhline(P_final_bar, color='gray', linestyle='--', label=f'P_final (Theoretical) = {P_final_bar:.2f} bar')
if transition_index > 0:
    plt.axvline(times[transition_index], color='red', linestyle=':', label=f'Choked/Non-choked transition ≈ {times[transition_index]:.1f} s')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (bar)')
plt.title('Reservoir Pressures vs Time')
plt.legend()
plt.grid(True)

# Plot Flow Rate
plt.subplot(2, 1, 2)
plt.plot(times, Q_slpm, label='Flow Rate (SLPM)')
# Optional: Plot mass flow rate as well
# plt.plot(times, Qm_kgs * 1e6, label='Flow Rate (mg/s)') # mg/s for visibility
if transition_index > 0:
    plt.axvline(times[transition_index], color='red', linestyle=':', label=f'Choked/Non-choked transition ≈ {times[transition_index]:.1f} s')
plt.xlabel('Time (s)')
plt.ylabel('Flow Rate (SLPM)')
plt.title('Hydrogen Flow Rate vs Time (Cv = {})'.format(Cv))
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()