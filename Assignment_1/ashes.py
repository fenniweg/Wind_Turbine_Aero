''' Use this file to call the BEM algorithm and compute the aerodynamics loads with your own BEM
code for Vo=5, 9, 11, and 20 m/s. Try and explain the source of any differences you may
see'''
import matplotlib.pyplot as plt
import numpy as np

import bem
from load_data import (
    A,
    P_rated,
    R,
    blade_dat,
    r_ashes,
    rho,
    v_5_ashes,
    v_9_ashes,
    v_11_ashes,
    v_20_ashes,
)

r_list = blade_dat['r'].values
r_compare = r_list / R
theta_opt = -0.11
omega_max = 0.98
cp_max = 0.47
v_rated = 11.19
optimum_lambda = omega_max * R / v_rated

v_compare = [5,9,11,20] #m/s, windpseeds to compare with ashe
p_n_compare_feather = np.zeros((len(v_compare), len(r_list)))
p_n_compare_stall = np.zeros((len(v_compare), len(r_list)))
p_t_compare_stall = np.zeros((len(v_compare), len(r_list)))
p_t_compare_feather = np.zeros((len(v_compare), len(r_list)))
p_n_ashes = [v_5_ashes[0:40], v_9_ashes[0:40], v_11_ashes[0:40], v_20_ashes[0:40]]
p_t_ashes = [v_5_ashes[40:80], v_9_ashes[40:80], v_11_ashes[40:80], v_20_ashes[40:80]]

for i, v in enumerate(v_compare):
    omega_max = (optimum_lambda*v)/R
    if v>= v_rated:
        omega_max = (optimum_lambda*v_rated)/R
    lambda_i = omega_max * R / v
    rpm = omega_max * 60 / (2*np.pi)
    print('Computing loads for V_0 =', v, 'm/s')
    print('Rotor speed =', rpm, 'rpm')
    print('Tip speed ratio =', lambda_i)
    cp_target = P_rated / (0.5 * rho * A * v**3)
    # Feather: increase pitch angle above the optimal setting to reduce Cp.
    # Stall: decrease pitch angle below the optimal setting to reduce Cp.
    if v >= v_rated:
        theta_f = bem.solve_pitch(0, 40.0, lambda_i, cp_target)
    else:
        theta_f = 0 #below rated wind speed, pitch angle is 0
    #compute final Cp, Ct, T and P for feather and stall pitching
    p_n_compare_feather[i,:],p_t_compare_feather[i,:]= bem.BEM_algorithm(lambda_i, theta_f,return_loads = True,V_0 = v)
    print('Feather pitch angle =', theta_f, 'deg')
    



# Separate plots for normal and tangential loads at each wind speed and compare with ashes data and save as pdf
for i, v in enumerate(v_compare):
    fig,axs = plt.subplots(1,2,figsize=(12,5))
    ax_n, ax_t = axs
    ax_n.plot(r_compare, p_n_compare_feather[i, :], 'b-', lw=2.5, label='Feather')
    ax_n.plot(r_ashes, p_n_ashes[i], 'k--', lw=2.5, label='ASHES')

    ax_n.set_xlabel('Blade radius $r/R$')
    ax_n.set_ylabel(r'$p_n$ [N/m]')
    ax_n.set_title(f'Normal load at $V_0$ = {v} m/s')
    ax_n.grid(True, alpha=0.3)

    ax_t.plot(r_compare, p_t_compare_feather[i, :], 'r--', lw=2.5, label='Feather')
    
    ax_t.plot(r_ashes, p_t_ashes[i], 'k--', lw=2.5, label='ASHES')
    ax_t.set_xlabel('Blade radius $r/R$')
    ax_t.set_ylabel(r'$p_t$ [N/m]')
    ax_t.set_title(f'Tangential load at $V_0$ = {v} m/s')
    ax_t.grid(True, alpha=0.3)

    ax_n.set_xlabel('Blade radius $r/R$')
    ax_t.set_xlabel('Blade radius $r/R$')
    ax_n.legend()
    ax_t.legend()
    fig.tight_layout()
    fig.savefig(rf'Figures\loads_{v}ms.pdf')
    plt.close(fig)


