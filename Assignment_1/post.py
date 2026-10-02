''' This file is used for post pprocessing of the results from main. It generates plots needed for the report and saves them in the Figures folder. It also saves the results in a text file for further analysis.
The plots generated are:
1. Cp as a function of tip speed ratio and pitch angle for both Polynomial and Madsen methods.
2. Power, pitch angle, thrust and Cp/Ct as a function of wind speed for feather and stall pitching.'''



import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from load_data import P_rated

#find global path
if not os.path.exists('Figures'):
    os.makedirs('Figures')

results_path = glob.glob('**/results', recursive=True)[0]
figures_path = glob.glob('**/Figures', recursive=True)[0]

#load npz files with results from main.py
npzfile = np.load(os.path.join(results_path, 'Cp_plotting_data.npz'))

theta_p = npzfile['theta_p']
tip_speed_ratio = npzfile['tip_speed_ratio']
THETA_GRID = npzfile['THETA_GRID']
LAMBDA_GRID = npzfile['LAMBDA_GRID']
Cp_polynomial = npzfile['Cp_polynomial']
Cp_madsen = npzfile['Cp_madsen']
Ct_polynomial = npzfile['Ct_polynomial']
Ct_madsen = npzfile['Ct_madsen']
cp_max_polynomial = npzfile['cp_max_polynomial']
cp_max_madsen = npzfile['cp_max_madsen']
optimum_lambda_polynomial = npzfile['optimum_lambda_polynomial']
optimum_theta_polynomial = npzfile['optimum_theta_polynomial']
optimum_lambda_madsen = npzfile['optimum_lambda_madsen']
optimum_theta_madsen = npzfile['optimum_theta_madsen']

npzfile = np.load(os.path.join(results_path, 'sweep_plotting_data.npz'))

v_sweep1 = npzfile['v_sweep']
omega_sweep = npzfile['omega_sweep']
P_sweep = npzfile['P_sweep']
omega_max = npzfile['omega_max']

npzfile = np.load(os.path.join(results_path, 'Q3_pitch_control_plotting_data.npz'))
v_sweep = npzfile['v_sweep']
theta_p_sweep_feather = npzfile['theta_p_sweep_feather']
theta_p_sweep_stall = npzfile['theta_p_sweep_stall']
Cp_sweep_feather = npzfile['Cp_sweep_feather']
Ct_sweep_feather = npzfile['Ct_sweep_feather']
Cp_sweep_stall = npzfile['Cp_sweep_stall']
Ct_sweep_stall = npzfile['Ct_sweep_stall']
P_sweep_feather = npzfile['P_sweep_feather']
P_sweep_stall = npzfile['P_sweep_stall']
T_sweep_feather = npzfile['T_sweep_feather']
T_sweep_stall = npzfile['T_sweep_stall']
v_rated = npzfile['v_rated']

print("Imported plotting data from 'results' folder for Q1, Q2 and Q3.")

'Q1: Contour plot of Cp and Ct as a function of tip speed ratio and pitch angle for both Madsen and Polynomial methods'

fig, axs = plt.subplots(2, 2, figsize=(15, 12), sharex=True, sharey=True)

for row, (method, cp_data, ct_data, theta_opt, lambda_opt) in enumerate([
    ('Madsen', Cp_madsen, Ct_madsen, optimum_theta_madsen, optimum_lambda_madsen),
    ('Polynomial', Cp_polynomial, Ct_polynomial,
     optimum_theta_polynomial, optimum_lambda_polynomial),
]):
    for ax, data, coefficient,cmap_c in [
        (axs[row, 0], cp_data, 'C_P','plasma'),
        (axs[row, 1], ct_data, 'C_T','viridis'),
    ]:
        
        contour = ax.contourf(THETA_GRID, LAMBDA_GRID, data, levels=20,
                              cmap=cmap_c, rasterized=True)
        
       
        if coefficient == 'C_P':
            ax.scatter(theta_opt, lambda_opt, color='b', s=80, marker='*',
                   label='Maximum $C_P$ ')
            ax.legend()
        if coefficient == 'C_T':
            ax.scatter(theta_opt, lambda_opt, color='r', s=80, marker='o',
                   label='$C_T$ ')
            ax.legend()
        ax.set_title(f'{method} Method - ${coefficient}(\\lambda, \\theta_p)$ Contour')
        ax.set_xlabel('Pitch Angle $\\theta_p$ [deg]')
        ax.set_ylabel('Tip Speed Ratio $\\lambda [-]$')
        ax.tick_params(labelsize=12)
        #add colorbar to each subplot
        cbar = fig.colorbar(contour, ax=ax)
        cbar.set_label(f'${coefficient}$ value [-]', fontsize=13)
        ax.grid(True)
        
#all fonts and labels are set to size 12 for consistency   
plt.rcParams.update({'font.size': 13})
plt.tight_layout()

fig.savefig(os.path.join(figures_path,'Cp_contour_comparison_raster.pdf'), dpi=300)


'Q2: Plot omega and P vs wind speed up to max speed'
# Plot 1, omega vs wind speed from cut-in to max speed
plt.figure(figsize=(10, 6))
plt.plot(v_sweep1, omega_sweep, 'b-', linewidth=2.8, label='Rotational Speed $\\omega(V_0)$')
# Add a vertical dashed line to highlight the rated wind speed point
plt.axvline(v_rated, color='red', linestyle='--', linewidth=2.0,
            label='Rated Wind Speed $V_{o,rated}$')

plt.axhline(omega_max, color='green', linestyle='--', linewidth=2.0,
            label='Maximum Rotational Speed $\\omega_{max}$')

# Format the plot with titles, labels, and grid
plt.title('DTU 10MW: Rotational Speed vs Wind Speed', fontsize=13)
plt.xlabel('Wind Speed $V_0$ [m/s]', fontsize=13)
plt.ylabel('Rotational Speed [rad/s]', fontsize=13)
plt.grid(True, linestyle=':', alpha=0.7, linewidth=1.0)
plt.legend(fontsize=13)
plt.tight_layout()
plt.savefig(os.path.join(figures_path,'omega_vs_wind_speed.pdf'), dpi=300,format='pdf', bbox_inches='tight')



# # Plot power against wind speed from cut-in to max speed
# plt.figure(figsize=(10, 6))
# plt.plot(v_sweep1, P_sweep, 'g-', linewidth=2.8, label='Power $P(V_0)$')
# plt.axvline(v_rated, color='red', linestyle='--', linewidth=2.0,
#             label=f'Rated Wind Speed = {v_rated:.2f} m/s')
# plt.axhline(P_rated, color='blue', linestyle='--', linewidth=2.0,
#             label=f'Rated Power = {P_rated:.2f} W')
# plt.grid(True, linestyle=':', alpha=0.7, linewidth=1.0)
# plt.title('DTU 10MW: Power vs Wind Speed', fontsize=13)
# plt.xlabel('Wind Speed $V_0$ [m/s]', fontsize=13)
# plt.ylabel('Power [W]', fontsize=13)
# plt.legend(fontsize=13)
# plt.tight_layout()
# plt.savefig('Figures/power_vs_wind_speed.pdf', dpi=300, format='pdf', bbox_inches='tight')

print('Plots for Q2 saved in Figures folder')


# ---------------------------------------------------------------------------
# Plot the required quantities.
# ---------------------------------------------------------------------------

# 4 seperate plots: 1) P vs V0, 2) theta_p vs V0, 3) T vs V0, 4) Cp and Ct vs V0 with each feather and stall pitch



plt.figure(figsize=(10, 6))
plt.plot(v_sweep, P_sweep_feather, 'b-', lw=2.8, label='Feather')
plt.plot(v_sweep, P_sweep_stall, 'r--', lw=2.8, label='Stall')
plt.axhline(P_rated, color='k', ls='--', lw=1.8, label='Rated power')
plt.xlabel('Wind speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Mechanical power $P$ [W]', fontsize=12)
plt.title('Power regulation by pitching', fontsize=12)
plt.legend(fontsize=13)
plt.tight_layout()
plt.grid(True, alpha=0.3)
plt.savefig(os.path.join(figures_path,'Q3_power_control.pdf'), dpi=300,format='pdf', bbox_inches='tight')


plt.figure(figsize=(10, 6))
plt.plot(v_sweep, theta_p_sweep_feather, 'b-', lw=2.8, label='Feather')
plt.plot(v_sweep, theta_p_sweep_stall, 'r--', lw=2.8, label='Stall')
plt.xlabel('Wind speed $V_0$ [m/s]', fontsize=12)
plt.ylabel(r'Pitch angle $\theta_p$ [deg]', fontsize=12)
plt.title(r'Pitch angle required for rated power at $\omega_{max}$', fontsize=12)
plt.legend(fontsize=13)
plt.tight_layout()
plt.grid(True, alpha=0.3)
plt.savefig(os.path.join(figures_path,'Q3_pitch_angle.pdf'), dpi=300,format='pdf', bbox_inches='tight')


plt.figure(figsize=(10, 6))
plt.plot(v_sweep, T_sweep_feather, 'b-', lw=2.8, label='Feather')
plt.plot(v_sweep, T_sweep_stall, 'r--', lw=2.8, label='Stall')
plt.xlabel('Wind speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Thrust $T$ [N]', fontsize=12)
plt.title('Thrust at rated-speed power limit', fontsize=12)
plt.legend(fontsize=13)
plt.tight_layout()
plt.grid(True, alpha=0.3)
plt.savefig(os.path.join(figures_path,'Q3_thrust.pdf'), dpi=300,format='pdf', bbox_inches='tight')


plt.figure(figsize=(10, 6))
plt.plot(v_sweep, Cp_sweep_feather, 'b-', lw=2.8, label=r'$C_p$ feather')
plt.plot(v_sweep, Cp_sweep_stall, 'r--', lw=2.8, label=r'$C_p$ stall')
plt.plot(v_sweep, Ct_sweep_feather, 'b:', lw=2.4, label=r'$C_T$ feather')
plt.plot(v_sweep, Ct_sweep_stall, 'r:', lw=2.4, label=r'$C_T$ stall')
plt.xlabel('Wind speed $V_0$ [m/s]', fontsize=12)
plt.ylabel('Coefficient value [-]', fontsize=12)
plt.title('Dimensionless coefficients for each pitch strategy', fontsize=12)
plt.legend(fontsize=13)
plt.tight_layout()
plt.grid(True, alpha=0.3)

plt.savefig(os.path.join(figures_path,'Q3_coefficients.pdf'), dpi=300,format='pdf', bbox_inches='tight')
print("Plots for Q3 saved in Figures folder")


