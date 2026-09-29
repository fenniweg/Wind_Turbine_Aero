# Q3: limit power at high wind speed by pitching the blades.
# Compute the pitch setting theta_p(V0) between cut-in and cut-out wind speed when
# the mechanical power is limited to P_rated and the rotational speed to omega_max.
# This is done for both feather and stall pitching.

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import brentq
from load_data import v_max, P_rated, R, rho,A
from old.Q2 import omega_max, rpm_max, v_rated,optimum_lambda, optimum_theta
from bem import BEM_algorithm

def function_to_solve(theta_p, lambda_i, cp_target):
    Cp  = BEM_algorithm(lambda_i, theta_p)
    return Cp - cp_target


def solve_pitch(theta_p_low, theta_p_high, lambda_i, cp_target):
    """Solve Cp(lambda_i, theta_p) = cp_target using Brent's method."""
    f_low = function_to_solve(theta_p_low, lambda_i, cp_target)
    f_high = function_to_solve(theta_p_high, lambda_i, cp_target)

    # At rated wind speed the target is exactly the maximum Cp, so the root is at the optimum pitch.
    if np.isclose(cp_target, BEM_algorithm(lambda_i, optimum_theta), rtol=1e-8, atol=1e-8):
        return optimum_theta
    if np.isclose(f_low, 0.0, atol=1e-10):
        return theta_p_low
    if np.isclose(f_high, 0.0, atol=1e-10):
        return theta_p_high

    return brentq(function_to_solve, theta_p_low, theta_p_high, args=(lambda_i, cp_target))

# ---------------------------------------------------------------------------
# Find the pitch setting that makes the power equal to the rated power.
# For feather pitching, a larger pitch angle reduces Cp.
# For stall pitching, a smaller pitch angle reduces Cp.
# --------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Sweep between cut-in and cut-out wind speed, leave out rated wind speed since it is already known.
# ---------------------------------------------------------------------------
v_sweep = np.linspace(v_min+0.1, v_max, 10)
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
    Cp_sweep_feather[i], Ct_sweep_feather[i],T_sweep_feather[i],P_sweep_feather[i] = BEM_algorithm(lambda_i, theta_f,Loads = True)
    Cp_sweep_stall[i], Ct_sweep_stall[i],T_sweep_stall[i],P_sweep_stall[i] = BEM_algorithm(lambda_i, theta_s,Loads = True)

# ---------------------------------------------------------------------------
# Plot the required quantities.
# ---------------------------------------------------------------------------

#4 plots: 1) P vs V0, 2) theta_p vs V0, 3) T vs V0, 4) Cp and Ct vs V0 with each feather and stall pitch

fig,axs = plt.subplots(2,2,figsize=(12,10))
axs[0,0].plot(v_sweep, P_sweep_feather, 'b-', lw=2.5, label='Feather')
axs[0,0].plot(v_sweep, P_sweep_stall, 'r--', lw=2.5, label='Stall')
axs[0,0].axhline(P_rated, color='k', ls='--', lw=1.5, label='Rated power')
axs[0,0].set_xlabel('Wind speed $V_0$ [m/s]')
axs[0,0].set_ylabel('Mechanical power $P$ [W]')
axs[0,0].set_title('Power regulation by pitching')
axs[0,0].grid(True, alpha=0.3)


axs[0,1].plot(v_sweep, theta_p_sweep_feather, 'b-', lw=2.5, label='Feather')
axs[0,1].plot(v_sweep, theta_p_sweep_stall, 'r--', lw=2.5, label='Stall')
axs[0,1].set_xlabel('Wind speed $V_0$ [m/s]')
axs[0,1].set_ylabel(r'Pitch angle $\theta_p$ [deg]')
axs[0,1].set_title(r'Pitch angle required for rated power at $\omega_{max}$')
axs[0,1].grid(True, alpha=0.3)

axs[1,0].plot(v_sweep, T_sweep_feather, 'b-', lw=2.5, label='Feather')
axs[1,0].plot(v_sweep, T_sweep_stall, 'r--', lw=2.5, label='Stall')
axs[1,0].set_xlabel('Wind speed $V_0$ [m/s]')
axs[1,0].set_ylabel('Thrust $T$ [N]')
axs[1,0].set_title('Thrust at rated-speed power limit')
axs[1,0].grid(True, alpha=0.3)

axs[1,1].plot(v_sweep, Cp_sweep_feather, 'b-', lw=2.5, label=r'$C_p$ feather')
axs[1,1].plot(v_sweep, Cp_sweep_stall, 'r--', lw=2.5, label=r'$C_p$ stall')
axs[1,1].plot(v_sweep, Ct_sweep_feather, 'b:', lw=2.0, label=r'$C_T$ feather')
axs[1,1].plot(v_sweep, Ct_sweep_stall, 'r:', lw=2.0, label=r'$C_T$ stall')
axs[1,1].set_xlabel('Wind speed $V_0$ [m/s]')
axs[1,1].set_ylabel('Coefficient value')    
axs[1,1].set_title('Dimensionless coefficients for each pitch strategy')
axs[1,1].grid(True, alpha=0.3)
for ax in axs.flat:
    ax.legend()

plt.tight_layout()

plt.savefig('Figures/Q3_pitch_control.png', dpi=300)
