"""
Question 1: Compute the highest obtainable power coefficient, Cp,max (λmax, θp,max).
Using BEM two methods the polynomial method and the Madsen et al. method. 
Compare the results of the two methods by plotting the power coefficient, Cp, as a function of the tip speed ratio, λ and pitch angles, θp. 

"""

### IMPORT LIBRARIES AND DATA ###

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bem import BEM_algorithm

#Define range of tip speed ratios (lambda) and pitch angles (theta_p) to test
tip_speed_ratio= np.linspace(5.0, 10.0, 10)      # Tip speed ratio range 
theta_p = np.linspace(-4.0, 3.0, 10)       # Pitch angle range [deg] 

#Initialize array to store Cp values for each method, with dimensions (len(tip_speed_ratio), len(theta_p), 2) for two methods
#Initialize tip speed ratio and pitch angle values in the first row and column of the array for each method to assign values to the correct indices later

Cp = np.zeros((len(tip_speed_ratio), len(theta_p), 2)) #Initialize array to store Cp values for each method


print("Starting BEM_algorithm computations for all combinations of tip speed ratios and pitch angles...")

#Outer loop over methods, inner loop over tip speed ratios and pitch angles
#Loop over all combinations of tip speed ratios and pitch angles to compute thrust, torque, and power using BEM_algorithm
for method in ['Polynomial','Madsen']:
    print(f"Running BEM_algorithm for method: {method}")

    for s in tip_speed_ratio:
        for theta in theta_p:
            Cp_value = BEM_algorithm(s,theta,method) #Call BEM_algorithm function to compute Cp for combinaion 
            Cp[np.where(tip_speed_ratio==s)[0][0], np.where(theta_p==theta)[0][0], 0 if method=='Polynomial' else 1] = Cp_value #Store Cp value in array for correct method
    print(f"Completed BEM_algorithm for method: {method}")

df = pd.DataFrame(Cp[:,:,0], index=tip_speed_ratio, columns=theta_p) #Create DataFrame for Polynomial method
df2 = pd.DataFrame(Cp[:,:,1], index=tip_speed_ratio, columns=theta_p) #Create DataFrame for Madsen method


#extract optimum
# Find the maximum Cp value and its corresponding tip speed ratio and pitch angle for each method
optimum_index_polynomial = np.unravel_index(np.argmax(Cp[:,:,0]), Cp[:,:,0].shape)
optimum_index_madsen = np.unravel_index(np.argmax(Cp[:,:,1]), Cp[:,:,1].shape)

cp_max_polynomial = Cp[optimum_index_polynomial[0], optimum_index_polynomial[1], 0]
cp_max_madsen = Cp[optimum_index_madsen[0], optimum_index_madsen[1], 1]

optimum_lambda_polynomial = tip_speed_ratio[optimum_index_polynomial[0]]
optimum_theta_polynomial = theta_p[optimum_index_polynomial[1]]

optimum_lambda_madsen = tip_speed_ratio[optimum_index_madsen[0]]
optimum_theta_madsen = theta_p[optimum_index_madsen[1]]
print("Search for maximum Cp values completed.")    
print("\n=================== RESULTS ===================")
#results
print(f"Polynomial Method: Cp_max = {cp_max_polynomial:.4f}, λ_max = {optimum_lambda_polynomial:.2f}, θp_max = {optimum_theta_polynomial:.2f}")
print(f"Madsen Method: Cp_max = {cp_max_madsen:.4f}, λ_max = {optimum_lambda_madsen:.2f}, θp_max = {optimum_theta_madsen:.2f}") 
#Save results to text file to use in Q2.py

print("\n=================== PLOTTING ===================")
# Contour plot of Cp as a function of tip speed ratio and pitch angle for both methods
THETA_GRID, LAMBDA_GRID = np.meshgrid(theta_p, tip_speed_ratio)

fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for ax, cp_data, title, theta_opt, lambda_opt in [
    (axs[0], Cp[:, :, 0], 'Polynomial Method - $C_p(\\lambda, \\theta_p)$ Contour',
     optimum_theta_polynomial, optimum_lambda_polynomial),
    (axs[1], Cp[:, :, 1], 'Madsen Method - $C_p(\\lambda, \\theta_p)$ Contour',
     optimum_theta_madsen, optimum_lambda_madsen),
]:
    contour = ax.contourf(THETA_GRID, LAMBDA_GRID, cp_data, levels=20, cmap='viridis')
    fig.colorbar(contour, ax=ax, label='$C_p$')
    ax.scatter(theta_opt, lambda_opt, color='r', s=80, marker='*', label='Max $C_p$')
    ax.set_title(title)
    ax.set_xlabel('Pitch Angle $\\theta_p$ [deg]')
    ax.set_ylabel('Tip Speed Ratio $\\lambda$')
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.savefig('Figures/Cp_contour_comparison.png', dpi=300)
print("Contour plots saved as 'Figures/Cp_contour_comparison.png'.")


            





