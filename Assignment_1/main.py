'''Calculate DTU Wind Energy 10MW wind turbine performance using BEM theory.
This code calculates the performance of the DTU Wind Energy 10MW wind turbine using Blade Element Momentum (BEM) theory. 
The code computes the power coefficient (Cp) as a function of tip speed ratio (lambda) and pitch angle (theta_p) using two methods: Polynomial and Interpolation. 
The rated wind speed and maximum rotational speed at rated wind speed are also calculated.
Stall ad feather pitching are implemented to limit the mechanical power to P_rated and the rotational speed to omega_max.
The results are plotted to visualize the performance characteristics of the wind turbine.'''

### IMPORT LIBRARIES AND DATA ###

import glob
import os

import numpy as np
from bem import BEM_algorithm, solve_bem, solve_pitch
from load_data import A, P_rated, R, rho, theta_p, tip_speed_ratio, v_max, v_min
from scipy.integrate import quad

#lookup global file path for results folder if exxistent otherwise create it
#check if results directory exists, if not create it
results_path = glob.glob('**/results', recursive=True)[0]



Cp = np.zeros((len(tip_speed_ratio), len(theta_p), 2)) #Initialize array to store Cp values for each method
Ct = np.zeros((len(tip_speed_ratio), len(theta_p), 2)) #Initialize array to store Ct values for each method
print("Starting BEM_algorithm computations for all combinations of tip speed ratios and pitch angles...")
v = 10
#Outer loop over methods, inner loop over tip speed ratios and pitch angles
#Loop over all combinations of tip speed ratios and pitch angles to compute thrust, torque, and power using BEM_algorithm
for method in ['Polynomial','Madsen']:
    print(f"Running BEM_algorithm for method: {method}")

    for s in tip_speed_ratio:
        for theta in theta_p:
            Cp_value,Ct_value = BEM_algorithm(s,theta,v,method) #Call BEM_algorithm function to compute Cp for combinaion 
            Cp[np.where(tip_speed_ratio==s)[0][0], np.where(theta_p==theta)[0][0], 0 if method=='Polynomial' else 1] = Cp_value #Store Cp value in array for correct method
            Ct[np.where(tip_speed_ratio==s)[0][0], np.where(theta_p==theta)[0][0], 0 if method=='Polynomial' else 1] = Ct_value #Store Ct value in array for correct method
    print(f"Completed BEM_algorithm for method: {method}")

# Find the maximum Cp value and its corresponding tip speed ratio and pitch angle for each method

optimum_index_polynomial = np.unravel_index(np.argmax(Cp[:,:,0]), Cp[:,:,0].shape)
optimum_index_madsen = np.unravel_index(np.argmax(Cp[:,:,1]), Cp[:,:,1].shape)

cp_max_polynomial = Cp[optimum_index_polynomial[    0], optimum_index_polynomial[1], 0]
cp_max_madsen = Cp[optimum_index_madsen[0], optimum_index_madsen[1], 1]

optimum_lambda_polynomial = tip_speed_ratio[optimum_index_polynomial[0]]
optimum_theta_polynomial = theta_p[optimum_index_polynomial[1]]

optimum_lambda_madsen = tip_speed_ratio[optimum_index_madsen[0]]
optimum_theta_madsen = theta_p[optimum_index_madsen[1]]
print("Search for maximum Cp values completed.")    
#results
print(f"Polynomial Method: Cp_max = {cp_max_polynomial:.4f}, λ_max = {optimum_lambda_polynomial:.2f}, θp_max = {optimum_theta_polynomial:.2f},Ct_max = {Ct[optimum_index_polynomial[0], optimum_index_polynomial[1], 0]:.4f}")
print(f"Madsen Method: Cp_max = {cp_max_madsen:.4f}, λ_max = {optimum_lambda_madsen:.2f}, θp_max = {optimum_theta_madsen:.2f},Ct_max = {Ct[optimum_index_madsen[0], optimum_index_madsen[1], 1]:.4f}") 

# Contour plot of Cp as a function of tip speed ratio and pitch angle for both methods
THETA_GRID, LAMBDA_GRID = np.meshgrid(theta_p, tip_speed_ratio)


# Save all data needed for plotting in a separate script with global path name
np.savez(os.path.join(results_path, 'Cp_plotting_data.npz'),
    theta_p=theta_p,
    tip_speed_ratio=tip_speed_ratio,
    THETA_GRID=THETA_GRID,
    LAMBDA_GRID=LAMBDA_GRID,
    Cp_polynomial=Cp[:, :, 0],
    Cp_madsen=Cp[:, :, 1],
    Ct_polynomial=Ct[:, :, 0],
    Ct_madsen=Ct[:, :, 1],
    cp_max_polynomial=cp_max_polynomial,
    cp_max_madsen=cp_max_madsen,
    optimum_lambda_polynomial=optimum_lambda_polynomial,
    optimum_theta_polynomial=optimum_theta_polynomial,
    optimum_lambda_madsen=optimum_lambda_madsen,
    optimum_theta_madsen=optimum_theta_madsen,
)
print("Saved plotting data to 'results/Cp_plotting_data.npz'.")

#Question 2: Compute rated wind speed and maximum rotational speed at rated wind speed

print("Compute rated wind speed and maximum rotational speed at rated wind speed")
C_p_max = cp_max_madsen #from Q1 results
optimum_lambda = optimum_lambda_madsen #from Q1 results
optimum_theta = optimum_theta_madsen #from Q1 results
#Calculate rated wind speed
v_rated = (P_rated/(0.5*rho*A*C_p_max))**(1/3)
#Calculate maximum rotational speed at rated wind speed
omega_max = (optimum_lambda*v_rated)/R #[rad/s]
rpm_max = omega_max*60/(2*np.pi) #[rpm]

print('Finding rated wind speed and maximum rotational speed at rated wind speed')
print('Madsen method: rated wind speed =', v_rated, 'm/s', 'maximum rotational speed at rated wind speed =', omega_max, 'rad/s (', rpm_max, 'rpm)')

C_p_max = cp_max_polynomial #from Q1 results
optimum_lambda = optimum_lambda_polynomial #from Q1 results
optimum_theta = optimum_theta_polynomial #from Q1 results

#Calculate rated wind speed
v_rated = (P_rated/(0.5*rho*A*C_p_max))**(1/3)

#Calculate maximum rotational speed at rated wind speed
omega_max = (optimum_lambda*v_rated)/R #[rad/s]
rpm_max = omega_max*60/(2*np.pi) #[rpm]

#Output results
print('Polynomial method: rated wind speed =', v_rated, 'm/s', 'maximum rotational speed at rated wind speed =', omega_max, 'rad/s (', rpm_max, 'rpm)')
v_sweep = np.linspace(v_min, v_rated, 20) #Sweep wind speed from cut-in to rated
omega_sweep = (optimum_lambda*v_sweep)/R #Calculate omega for each wind speed
rpm_sweep = omega_sweep*60/(2*np.pi) #Convert omega to rpm
P_sweep = 0.5*rho*A*C_p_max*v_sweep**3 #Calculate power for each wind speed
v_sweep = np.append(v_sweep, np.linspace(v_rated, v_max, 100)) #Sweep wind speed from rated to max
rpm_sweep = np.append(rpm_sweep, np.full(100, rpm_max)) #Omega is constant at omega_max above rated wind speed
P_sweep = np.append(P_sweep, np.full(100, P_rated)) #Power is constant at rated power above rated wind speed
omega_sweep = np.append(omega_sweep, np.full(100, omega_max)) #Omega is constant at omega_max above rated wind speed
# Save sweep data for use in other scripts/plotting
np.savez(os.path.join(results_path, 'sweep_plotting_data.npz'),
         v_sweep=v_sweep,
         omega_sweep=omega_sweep,
         P_sweep=P_sweep,
         v_rated=v_rated,
         omega_max=omega_max)
print("Saved plotting data to 'results/sweep_plotting_data.npz'.")


### QUESTION 3: Implement stall and feather pitching to limit the mechanical power to P_rated and the rotational speed to omega_max.

# ---------------------------------------------------------------------------
# Find the pitch setting that makes the power equal to the rated power.
# For femather pitching, a larger pitch angle reduces Cp.
# For stall pitching, a smaller pitch angle reduces Cp.
# --------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Sweep between cut-in and cut-out wind speed, leave out rated wind speed since it is already known.
# ---------------------------------------------------------------------------
print("Compute pitch angle, Cp, Ct, T and P as a function of wind speed from rated to max speed for feather and stall pitching")
v_sweep = np.linspace(v_rated+0.1, v_max, 10)
theta_p_sweep_feather = np.zeros_like(v_sweep)
theta_p_sweep_stall = np.zeros_like(v_sweep)
Cp_sweep_feather = np.zeros_like(v_sweep)
Ct_sweep_feather = np.zeros_like(v_sweep)
Cp_sweep_stall = np.zeros_like(v_sweep)
Ct_sweep_stall = np.zeros_like(v_sweep)
P_sweep_feather = np.zeros_like(v_sweep)
P_sweep_stall = np.zeros_like(v_sweep)
T_sweep_feather = np.zeros_like(v_sweep)
T_sweep_stall = np.zeros_like(v_sweep)

for i, v in enumerate(v_sweep):
    lambda_i = omega_max * R / v

    cp_target = P_rated / (0.5 * rho * A * v**3)

    # Feather: increase pitch angle above the optimal setting to reduce Cp.
    # Stall: decrease pitch angle below the optimal setting to reduce Cp.
    theta_f = solve_pitch(0.0, 40.0, lambda_i, cp_target)
    theta_s = solve_pitch(-40.0, 0.0, lambda_i, cp_target)
    theta_p_sweep_feather[i] = theta_f
    theta_p_sweep_stall[i] = theta_s

    #compute final Cp, Ct, T and P for feather and stall pitching
    Cp_sweep_feather[i], Ct_sweep_feather[i] = solve_bem(lambda_i, theta_f, v)
    Cp_sweep_stall[i],Ct_sweep_stall[i] = solve_bem(lambda_i, theta_s,v)
    P_sweep_feather[i] = 0.5*rho*A*Cp_sweep_feather[i]*v**3
    P_sweep_stall[i] = 0.5*rho*A*Cp_sweep_stall[i]*v**3
    T_sweep_feather[i] = Ct_sweep_feather[i]*0.5*rho*A*v**2
    T_sweep_stall[i] = Ct_sweep_stall[i]*0.5*rho*A*v**2
    
#get cp,ct,power,thrust,pitch from 0 to rated wind speed for plotting
v_sweep_below_rated = np.linspace(v_min, v_rated, 15)
Cp_sweep_below_rated = np.zeros_like(v_sweep_below_rated)
Ct_sweep_below_rated = np.zeros_like(v_sweep_below_rated)
T_sweep_below_rated = np.zeros_like(v_sweep_below_rated)
P_sweep_below_rated = np.zeros_like(v_sweep_below_rated)

for i, v in enumerate(v_sweep_below_rated):
    omega_i = (optimum_lambda*v)/R
    lambda_i = omega_i * R / v
    Cp_sweep_below_rated[i], Ct_sweep_below_rated[i] = solve_bem(lambda_i, optimum_theta,v)
    P_sweep_below_rated[i] = 0.5*rho*A*Cp_sweep_below_rated[i]*v**3
    T_sweep_below_rated[i] = Ct_sweep_below_rated[i]*0.5*rho*A*v**2

#put arrays together for plotting
v_sweep = np.concatenate((v_sweep_below_rated, v_sweep))
theta_p_sweep_feather = np.concatenate((np.full_like(v_sweep_below_rated, optimum_theta), theta_p_sweep_feather))
theta_p_sweep_stall = np.concatenate((np.full_like(v_sweep_below_rated,optimum_theta),theta_p_sweep_stall))
Cp_sweep_feather = np.concatenate((Cp_sweep_below_rated, Cp_sweep_feather))
Ct_sweep_feather = np.concatenate((Ct_sweep_below_rated, Ct_sweep_feather))
Cp_sweep_stall = np.concatenate((Cp_sweep_below_rated, Cp_sweep_stall))
Ct_sweep_stall = np.concatenate((Ct_sweep_below_rated, Ct_sweep_stall))
P_sweep_feather = np.concatenate((P_sweep_below_rated, P_sweep_feather))
P_sweep_stall = np.concatenate((P_sweep_below_rated, P_sweep_stall))
T_sweep_feather = np.concatenate((T_sweep_below_rated, T_sweep_feather))
T_sweep_stall = np.concatenate((T_sweep_below_rated, T_sweep_stall))

np.savez(os.path.join(results_path, 'Q3_pitch_control_plotting_data.npz'),
         v_sweep=v_sweep,
            theta_p_sweep_feather=theta_p_sweep_feather,
            theta_p_sweep_stall=theta_p_sweep_stall,
            Cp_sweep_feather=Cp_sweep_feather,
            Ct_sweep_feather=Ct_sweep_feather,
            Cp_sweep_stall=Cp_sweep_stall,
            Ct_sweep_stall=Ct_sweep_stall,
            P_sweep_feather=P_sweep_feather,
            P_sweep_stall=P_sweep_stall,
            T_sweep_feather=T_sweep_feather,
            T_sweep_stall=T_sweep_stall,
            v_rated=v_rated,
            rpm_max=rpm_max,
            P_rated=P_rated)
print("Saved plotting data to 'results/Q3_pitch_control_plotting_data.npz")
print('Simulation done :)')

print('Calculating anual energy output AEP with cut-out at 25 m/s and 20 m/s and energy lost between the two cut-out speeds')

A_weibull = 9.0
k_weibull = 1.9

V_cut_in = 4.0
V_cut_out_25 = 25.0
V_cut_out_20 = 20.0

area = np.pi * R**2

V_rated = (P_rated/(0.5 * rho * area * C_p_max ))**(1/3)


# Weibull probability density function
def weibull(V):
    return (k_weibull / A_weibull* (V / A_weibull)**(k_weibull - 1)* np.exp(-(V / A_weibull)**k_weibull))


# Wind turbine power curve
def power(V, V_cut_out):
    if V < V_cut_in:
        return 0.0
    elif V < V_rated:
        return 0.5*rho*area*C_p_max *V**3
    elif V < V_cut_out:
        return P_rated
    else:
        return 0.0


# AEP with 25 m/s 
AEP_25 = 8760 * quad(lambda V: power(V, V_cut_out_25) * weibull(V),
    0,
    40
)[0]


# AEP with 20 m/s 
AEP_20 = 8760 * quad(
    lambda V: power(V, V_cut_out_20) * weibull(V),
    0,
    40
)[0]


# Energy lost between 20 and 25 m/s
energy_lost = AEP_25 - AEP_20


AEP_25_GWh = AEP_25 / 1e9
AEP_20_GWh = AEP_20 / 1e9
energy_lost_GWh = energy_lost / 1e9



print(f"AEP with cut-out at 25 m/s = {AEP_25_GWh:.2f} GWh/year")
print(f"AEP with cut-out at 20 m/s = {AEP_20_GWh:.2f} GWh/year")
print(f"Energy lost = {energy_lost_GWh:.2f} GWh/year")


