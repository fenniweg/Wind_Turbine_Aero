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
r_compare = r_list/R
theta_opt = 0.00
omega_max = 1.01
cp_max = 0.4661
v_rated = 11.43
optimum_lambda = 7.86

v_compare = [5,9,11,20] #m/s, windpseeds to compare with ashes
p_n_compare = np.zeros((len(v_compare), len(r_list)))
p_t_compare = np.zeros((len(v_compare), len(r_list)))
p_n_ashes = [v_5_ashes[0:40], v_9_ashes[0:40], v_11_ashes[0:40], v_20_ashes[0:40]]
p_t_ashes = [v_5_ashes[40:80], v_9_ashes[40:80], v_11_ashes[40:80], v_20_ashes[40:80]]

for i, v in enumerate(v_compare):
    if v >= v_rated:
            #if wind speed is above rated, use rated omega to compute loads
            omega_i = omega_max
            lambda_i = omega_max * R / v
            cp_target = P_rated / (0.5 * rho * A * v**3)
            theta_i = bem.solve_pitch(0, 40.0, lambda_i, cp_target)
            
    else:
         
        omega_i = (optimum_lambda*v)/R
        lambda_i = omega_i * R / v
        theta_i = theta_opt #below rated wind speed, pitch angle is opt

    print('Computing loads for V_0 =', v, 'm/s with theta_i =', theta_i, 'deg')

    p_n_compare[i,:],p_t_compare[i,:]= bem.BEM_algorithm(lambda_i, theta_i,v,return_loads = True)
    #tip speed
    rpm = omega_i * 60 / (2*np.pi)
    
    print('Omega =', omega_i, 'rad/s')
    print('Rotor speed =', rpm, 'rpm')
    print('Tip speed ratio =', lambda_i)
    print('Pitch angle =', theta_i, 'deg')
    cp_i,_ = bem.BEM_algorithm(lambda_i, theta_i,v)
    print('Cp =', cp_i)

    
    

    



# Separate plots for normal and tangential loads at each wind speed and compare with ashes data and save as pdf
for i, v in enumerate(v_compare):
    fig,axs = plt.subplots(1,2,figsize=(12,5))
    ax_n, ax_t = axs
    ax_n.plot(r_compare, p_n_compare[i, :], 'b-', lw=2.5, label='BEM Code')
    ax_n.plot(r_ashes, p_n_ashes[i], 'k--', lw=2.5, label='ASHES')

    ax_n.set_xlabel('Blade radius $r/R$')
    ax_n.set_ylabel(r'$p_n$ [N/m]')
    ax_n.set_title(f'Normal load at $V_0$ = {v} m/s')
    ax_n.grid(True, alpha=0.3)

    ax_t.plot(r_compare, p_t_compare[i, :], 'r--', lw=2.5, label='BEM Code')
    
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


print('Plots saved in Figures folder')


