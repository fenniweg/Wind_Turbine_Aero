''''
Question 2 Run the wind turbine at max C_p all the way to rated Power. Find rated wind speed V_= rated and omega_max- PLot omega(V_0)

'''

### IMPORT LIBRARIES AND DATA ###
import numpy as np
import matplotlib.pyplot as plt

from load_data import v_min,v_max,P_rated,R,rho
from old.Q1 import cp_max_polynomial, optimum_lambda_polynomial, optimum_theta_polynomial

print("Compute rated wind speed and maximum rotational speed at rated wind speed")
C_p_max = cp_max_polynomial #from Q1 results
optimum_lambda = optimum_lambda_polynomial #from Q1 results
optimum_theta = optimum_theta_polynomial #from Q1 results


#Calculate rotor area
A = np.pi*R**2

#Calculate rated wind speed
v_rated = (P_rated/(0.5*rho*A*C_p_max))**(1/3)

#Calculate maximum rotational speed at rated wind speed
omega_max = (optimum_lambda*v_rated)/R #[rad/s]
rpm_max = omega_max*60/(2*np.pi) #[rpm]

#Output results
print(f"Rated wind speed: {v_rated:.2f} m/s")
print(f"Maximum rotational speed at rated wind speed: {omega_max:.2f} rad/s ({rpm_max:.2f} rpm)")

print('Compute omega as a function of wind speed from cut-in to max speed')
print('Compute power as a function of wind speed from cut-in to max speed ')

v_sweep = np.linspace(v_min, v_rated, 100) #Sweep wind speed from cut-in to rated
omega_sweep = (optimum_lambda*v_sweep)/R #Calculate omega for each wind speed
rpm_sweep = omega_sweep*60/(2*np.pi) #Convert omega to rpm
P_sweep = 0.5*rho*A*C_p_max*v_sweep**3 #Calculate power for each wind speed
v_sweep = np.append(v_sweep, np.linspace(v_rated, v_max, 100)) #Sweep wind speed from rated to max
rpm_sweep = np.append(rpm_sweep, np.full(100, rpm_max)) #Omega is constant at omega_max above rated wind speed
P_sweep = np.append(P_sweep, np.full(100, P_rated)) #Power is constant at rated power above rated wind speed

#Plot 1, omega vs wind speed from cut-in to max speed
plt.figure(figsize=(10,6))
plt.plot(v_sweep, rpm_sweep, 'b-', linewidth=2.5, label='Rotational Speed $\\omega(V_0)$')
# Add a vertical dashed line to highlight the rated wind speed point
plt.axvline(v_rated, color='red', linestyle='--', 
            label=f'Rated Wind Speed = {v_rated:.2f} m/s')

plt.axhline(rpm_max, color='green', linestyle='--',
            label=f'Maximum Rotational Speed = {rpm_max:.2f} rpm')

# Format the plot with titles, labels, and grid
plt.title('DTU 10MW: Rotational Speed vs Wind Speed', fontsize=14)
plt.xlabel('Wind Speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Rotational Speed [RPM]', fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig('Figures/omega_vs_wind_speed.png', dpi=300)


#Plot power against wind speed from cut-in to max speed
plt.figure(figsize=(10,6))
plt.plot(v_sweep, P_sweep, 'g-', linewidth=2.5, label='Power $P(V_0)$')
plt.axvline(v_rated, color='red', linestyle='--', 
            label=f'Rated Wind Speed = {v_rated:.2f} m/s')
plt.axhline(P_rated, color='blue', linestyle='--',
            label=f'Rated Power = {P_rated:.2f} W')
plt.grid(True, linestyle=':', alpha=0.7)
plt.title('DTU 10MW: Power vs Wind Speed', fontsize=14)
plt.xlabel('Wind Speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Power [W]', fontsize=12)
plt.legend(fontsize=11)
plt.savefig('Figures/power_vs_wind_speed.png', dpi=300)

print('Plots saved in Figures folder')
