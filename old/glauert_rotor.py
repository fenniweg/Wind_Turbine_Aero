##Exercsie 1
# calculate Glauert optimum rotor

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from load_data import blade_dat, airfoil_data


#function to compute effiecency for one tip speed integrating over r/R whole rotor blade
def C_p_lambda(tip_speed, R, r, B, F):
    C_p = 0
    # radial discretization
    
    


#function to compute efficiency for given values of lambda  integrating over r/R whole rotor blade

def c_p_lambda(list_lambda, R, B, C_l, C_d, F):
    list_c_p = []

    # radial discretization
    mu = np.linspace(0.01, 0.99, 200)  # mu = r/R

    for l in list_lambda:

        integrand = []

        for m in mu:

            x = l * m

            # Newton iteration for a
            a = 0.33
            tol = 1e-6
            max_iter = 100

            for _ in range(max_iter):
                f_a = 16*a**3 - 24*a**2 + a*(9 - 3*x**2) - 1 + x**2
                df_da = 48*a**2 - 48*a + 9 - 3*x**2

                a_new = a - f_a/df_da

                if abs(a_new - a) < tol:
                    break

                a = a_new

            a_prime = (1 - 3*a)/(4*a - 1)

            psi = np.arctan((1-a)/((1+a_prime)*x))

            velocity_ratio = np.sin(psi)/(1-a)

            C_n = C_l*np.cos(psi) + C_d*np.sin(psi)
            C_t = C_l*np.sin(psi) - C_d*np.cos(psi)

            chord = (
                8*a*(1-a)*np.pi*F*m
                * velocity_ratio**2
                /(B*C_n) * R
            )

            # local integrand
            f = (
                m
                * (chord/R)
                * (1-a)**2
                / np.sin(psi)**2
                * C_t
            )

            integrand.append(f)

        # numerical integration over r/R
        integral = np.trapz(integrand, mu)

        CP = l * B / np.pi * integral

        list_c_p.append(CP)

    return list_c_p


# #a
# R = 10 #m
# B = 3
# alpha_design = 4 #degree

# C_l = 0.8
# C_d = 0
# F = 1
# list_lambda = np.arange(2,20)  
# list_c_p = c_p_lambda(list_lambda,R,B,C_l,C_d,F)

# #b
# C_d = 0.006
# list_lambda_b = np.arange(2,16)
# list_c_p_b = c_p_lambda(list_lambda_b,R,B,C_l,C_d,F)

# #plot C_P from a and b vs lambda
# plt.figure(figsize=(8, 6))
# plt.subplot(1, 2, 1)
# plt.plot(list_lambda, list_c_p)
# plt.xlabel('Tip speed ratio (lambda)')
# plt.ylabel(' (C_P)')    
# plt.grid()
# plt.subplot(1, 2, 2)
# plt.plot(list_lambda_b, list_c_p_b, 'r--')
# plt.xlabel('Tip speed ratio (lambda)')
# plt.ylabel(' (C_P)')
# plt.grid()
# plt.show()


#c
l = 8
R = 10
B = 3
alpha_design = 4
C_l = 0.8
C_d = 0.006
F =1

a_list = []
chord_list = []
twist_list = []
r_list = []
r_over_R_list = np.arange(0.001, 1, 0.001)

for r_over_R in r_over_R_list:
    r = r_over_R * R
    r_list.append(r_over_R)
    a = 0.33
    tol = 1e-6
    max_iter = 100
    x = l*r_over_R
    for i in range(max_iter):
            f_a = 16*a**3 - 24*a**2 + a*(9 - 3*x**2) - 1 + x**2
            df_da = 48*a**2 - 48*a + 9 - 3*x**2
            a_new = a - f_a/df_da
            if abs(a_new - a) < tol:
                break
            a = a_new

    a_list.append(a)
    a_prime = (1-3*a)/(4*a-1)
    psi = np.arctan((1-a)/((1+a_prime)*x))
    print(f"a': {psi}")
    
    velocity_ratio = np.sin(psi)/(1-a)
    
            
    C_n = C_l*np.cos(psi) + C_d*np.sin(psi)
    C_t = C_l*np.sin(psi) - C_d*np.cos(psi)
    twist = psi*180/np.pi - alpha_design
    twist_list.append(twist)
    
    chord = R*8*a*(1-a)*np.pi*F*(r_over_R)*velocity_ratio**2/(B*C_n)
    chord_list.append(chord)


#plot a, chord and twist vs r/R
plt.figure(figsize=(12, 5)) 
plt.subplot(1, 3, 1)
plt.plot(r_list, a_list)
plt.xlabel('r/R')
plt.ylabel('Axial Induction Factor (a)')
plt.grid()
plt.subplot(1, 3, 2)
plt.plot(r_list, chord_list)
plt.xlabel('r/R')
plt.ylabel('Chord Length (m)')
plt.grid()
plt.subplot(1, 3, 3)
plt.plot(r_list, twist_list)
plt.xlabel('r/R')
plt.ylabel('Twist Angle (degrees)')
plt.tight_layout()
plt.grid()

plt.show()

    


    


