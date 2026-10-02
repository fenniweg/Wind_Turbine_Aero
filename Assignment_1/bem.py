'''This file contains functions to be used for Assignemnt 1 of Wind Energy course. 
The functions include BEM_algorithm, double_interpolation, function_to_solve, and solve_pitch.
'''

import numpy as np
from load_data import A, R, V_rated, airfoil_data, blade_dat, n_blades, rho
from scipy.optimize import brentq

omega_max = 1.01


def double_interpolation(alpha,t_over_c):
    """Double interpolation to find C_l and C_d for a given angle of attack (alpha) and thickness/chord ratio (t_over_c)."""
    C_l_thickness = np.zeros(6)
    C_d_thickness = np.zeros(6)
    
    for k in range (6):
        C_l_thickness[k] = np.interp(alpha,airfoil_data[k][:,0],airfoil_data[k][:,1])
        C_d_thickness[k] = np.interp(alpha,airfoil_data[k][:,0],airfoil_data[k][:,2])

    #Interpolate to actual thickness
    thickness_profile = [100,60,48,36,30.1,24.1] #thickness profile of the 6 airfoils in percentage
    #Use np.argsort to sort the thickness profile and corresponding C_l and C_d values
    sort_indices = np.argsort(thickness_profile)
    x_p_sorted = np.array(thickness_profile)[sort_indices]
    y_p_sorted_C_l = np.array(C_l_thickness)[sort_indices]
    y_p_sorted_C_d = np.array(C_d_thickness)[sort_indices]

    C_l = np.interp(t_over_c,x_p_sorted,y_p_sorted_C_l)
    C_d = np.interp(t_over_c,x_p_sorted,y_p_sorted_C_d)
    return C_l, C_d



def function_to_solve(theta_p, lambda_i, cp_target):
    """Function to solve for the pitch angle (theta_p) that achieves a target power coefficient (cp_target) at a given tip speed ratio (lambda_i)."""
    Cp,_  = BEM_algorithm(lambda_i, theta_p,V_rated)
    return Cp - cp_target


def solve_pitch(theta_p_low, theta_p_high, lambda_i, cp_target):
    """Solve Cp(lambda_i, theta_p) = cp_target using Brent's method."""
    optimum_theta = -0.11 # Initialize optimum_theta to a default value
    f_low = function_to_solve(theta_p_low, lambda_i, cp_target)
    f_high = function_to_solve(theta_p_high, lambda_i, cp_target)

    # At rated wind speed the target is exactly the maximum Cp, so the root is at the optimum pitch.
    cp,_ = BEM_algorithm(lambda_i, optimum_theta,V_rated)
    if np.isclose(cp_target, cp, rtol=1e-8, atol=1e-8):
        return optimum_theta
    if np.isclose(f_low, 0.0, atol=1e-10):
        return theta_p_low
    if np.isclose(f_high, 0.0, atol=1e-10):
        return theta_p_high

    return brentq(function_to_solve, theta_p_low, theta_p_high, args=(lambda_i, cp_target))

def solve_bem(lambda_tip, theta_p_deg,V_0, method='Polynomial', return_loads=False):
    if V_0 >= V_rated:
        omega = omega_max
    else:
        omega = lambda_tip * V_0 / R
   #Load blade data
    N_elem = len(blade_dat)
    r_list = blade_dat['r'].values
    chord_list = blade_dat['c'].values
    beta_list = blade_dat['beta'].values
    t_over_c_list = blade_dat['t/c'].values
   
    p_n_list = np.zeros(len(r_list))
    p_t_list = np.zeros(len(r_list))
    
    # Loop over radial elements up to penultimate node to prevent Prandtl tip divergence (Assignment #1 tip)
    for i in range(N_elem - 1): #parte dal primo elemento 2.8 e arriva al penultimo 88.45
        r_i = r_list[i]
        c_i = chord_list[i]
        beta_i = beta_list[i]
        tc_i = t_over_c_list[i]
        
        sigma_i = (c_i * n_blades) / (2 * np.pi * r_i)  # Solidity (Slide 15)
        
        # BEM Iteration Setup
        a = 0.0
        a_prime = 0.0
        f_relax = 0.1
        
        for icount in range(1000):  # Safety iteration limit (Slide 5)
            # Flow Angle phi
            tan_phi = ((1.0 - a) * V_0) / ((1.0 + a_prime) * omega * r_i)
            phi = np.arctan(tan_phi)
            
            # Local Angle of Attack alpha (Slide 13)
            theta_loc = beta_i + theta_p_deg
            alpha_deg = np.degrees(phi) - theta_loc
            
            # Double Interpolation for Lift and Drag
            Cl, Cd = double_interpolation(alpha_deg,tc_i)
            
            # Load Coefficients Cn and Ct (Slide 15/16)
            Cn = Cl * np.cos(phi) + Cd * np.sin(phi)
            Ct = Cl * np.sin(phi) - Cd * np.cos(phi)
            
            # Prandtl Tip Loss Correction F (Slide 14)
            sin_phi_abs = np.abs(np.sin(phi))
            if sin_phi_abs < 1e-6:
                F = 1.0
            else:
                exp_arg = -(n_blades/ 2.0) * (R - r_i) / (r_i * sin_phi_abs)
                F = (2.0 / np.pi) * np.arccos(np.exp(exp_arg))
            F = max(F, 1e-4)  # Avoid division by zero
            
            # Calculate CT local
            dCT = ((1.0 - a)**2 * Cn * sigma_i) / (np.sin(phi)**2)
            
            # Axial Induction Update (Slide 1 & Assignment #1)
            if method == 'Polynomial':  # Equation (1) Classical Glauert
                if a <= 0.33:
                    a_star = (sigma_i * Cn) / (4.0 * F * np.sin(phi)**2) * (1.0 - a)
                else:
                    a_star = dCT / (4.0 * F * (1.0 - 0.25 * (5.0 - 3.0 * a) * a))
            else:  # Equation (2) Madsen et al.
                CT_F = dCT / F
                a_star = 0.246 * CT_F + 0.0586 * (CT_F**2) + 0.0883 * (CT_F**3)
                
            # Tangential Induction Update
            a_prime_star = (sigma_i * Ct) / (4.0 * F * np.sin(phi) * np.cos(phi)) * (1.0 + a_prime)
            
            # Apply Underrelaxation (Slide 1)
            a_new = f_relax * a_star + (1.0 - f_relax) * a
            a_prime_new = f_relax * a_prime_star + (1.0 - f_relax) * a_prime
            
            # Convergence check
            if abs(a_new - a) < 1e-6 and abs(a_prime_new - a_prime) < 1e-6:
                a, a_prime = a_new, a_prime_new
                break
                
            a, a_prime = a_new, a_prime_new
            
        # Element loads [N/m] (Slide 17)
        Vrel = np.sqrt((V_0 * (1.0 - a))**2 + (omega * r_i * (1.0 + a_prime))**2)
        p_n_list[i] = 0.5 * rho * (Vrel**2) * c_i * Cn
        p_t_list[i] = 0.5 * rho * (Vrel**2) * c_i * Ct

    # Force last element at the tip r=R to 0 N/m (Assignment #1 tip)
    p_n_list[-1] = 0.0
    p_t_list[-1] = 0.0

    if return_loads == True:
            return  p_n_list, p_t_list
    
    # Integrate loads across the blade span (Slide 18)
    Thrust = n_blades * np.trapz(p_n_list, r_list)
    Torque = n_blades * np.trapz(r_list * p_t_list, r_list)
    Power = omega * Torque
    
    # Dimensionless Coefficients
    A = np.pi * R**2
    Cp = Power / (0.5 * rho * A * V_0**3)
    CT = Thrust / (0.5 * rho * A * V_0**2)
    
    

    return Cp, CT

def BEM_algorithm (s,theta_p,V_0,method = 'Polynomial',return_loads = False):
    '''
    BEM_algorithm computes the power coefficient (Cp) and thrust coefficient (CT) for a given tip speed ratio (s), 
    pitch angle (theta_p), and method ('Polynomial' or 'Madsen'). Polynomial is standard.'''
    if V_0 > V_rated:
            omega = omega_max
    else:
            omega = s * V_0 / R
     #Load blade data
    r_list = blade_dat['r'].values
    chord_list = blade_dat['c'].values
    beta_list = blade_dat['beta'].values
    t_over_c_list = blade_dat['t/c'].values

    p_n_list = np.zeros(len(r_list))
    p_t_list = np.zeros(len(r_list))

    #loop through each blade element
    for i in range(len(r_list)-1): #last element is tip, skip to avoid numerical issues
        r = r_list[i]
        chord = chord_list[i]
        beta = beta_list[i]
        t_over_c = t_over_c_list[i]

        sigma = chord*n_blades/(2*np.pi*r) 

    #Initialze a and a_prime, convergence tolerance

        a = 0
        a_prime = 0
        epsilon = 1e-6

        count = 0

        while True:
            count += 1
            #Calculate flowangle
            psi = np.arctan((1-a)*V_0/((1+a_prime)*omega*r)) #radians
    
            #Compute local angle of attack alpha
            alpha = np.degrees(psi)-(beta+theta_p) #degrees
    
            #Lookup C_l and C_d from airfoil data based on alpha with double interpolation
            C_l, C_d = double_interpolation(alpha,t_over_c)

            #Compute normal and tangential force coefficients
            C_n = C_l*np.cos(psi) + C_d*np.sin(psi)
            C_t = C_l*np.sin(psi) - C_d*np.cos(psi)
    
            
            dC_T = (1-a)**2*C_n*sigma/(np.sin(psi)**2)
    
            # Prandtl Tip Loss Correction F 
            if np.abs(np.sin(psi)) < 1e-6:
                F = 1.0
            else:
                exp_arg = -(n_blades / 2.0) * (R - r) / (r * np.abs(np.sin(psi)))
                F = (2.0 / np.pi) * np.arccos(np.exp(exp_arg))
                F = max(F, 1e-4)  # Avoid division by zero

            #Axial Induction Factor a and Tangential Induction Factor a_prime for chosen method
            if method  == 'Polynomial':
                if a <= 0.33: #no correction for a_star applied
                    a_star = sigma*C_n/(4*F*np.sin(psi)**2)*(1-a)
                else: # a_star with correction
                    a_star = dC_T / (4.0 * F * (1.0 - 0.25 * (5.0 - 3.0 * a) * a))
            elif method =='Madsen': #Madsen et al
                    CT_F = dC_T / F
                    a_star = 0.246 * CT_F + 0.0586 * (CT_F**2) + 0.0883 * (CT_F**3)

            else:
                raise ValueError("Invalid method. Choose 'Polynomial' or 'Madsen'.")

            #Update a and a_prime with relaxation factor 0.1
            a_new = 0.1*a_star+(1-0.1)*a #Update a with relaxation factor 0.1

            a_prime_star = sigma*C_t/(4*F*np.sin(psi)*np.cos(psi))*(1+a_prime)
            a_prime_new = 0.1*a_prime_star+(1-0.1)*a_prime

            #Check for convergence
            if abs(a_new - a) < epsilon and abs(a_prime_new - a_prime) < epsilon:
                # print(f"Converged after {i} iterations")
                a, a_prime = a_new, a_prime_new
                break
            else:
                a = a_new
                a_prime = a_prime_new

            if count > 1000:
                print(f"Warning: BEM did not converge after 1000 iterations in element {i}")
                #break and go to next element
        
                break
            


        #if converged, calculate p_n and p_t for this blade element
            

        V_rel = np.sqrt((V_0 * (1.0 - a))**2 + (omega * r * (1.0 + a_prime))**2)

        p_n = 0.5*rho*(V_rel)**2*chord*C_n
        p_t = 0.5*rho*(V_rel)**2*chord*C_t

        #Save p_n and p_t for each blade element
        p_n_list[i] = (p_n)
        p_t_list[i]=(p_t)
        

    #Force last element to be zero to avoid numerical issues at the tip
    p_n_list[-1] = 0.0
    p_t_list[-1] = 0.0
    if return_loads == True:
        return p_n_list,p_t_list
    # Integrate over the blade elements to get total thrust and torque
    thrust = n_blades*np.trapz(p_n_list,r_list)
    torque = n_blades*np.trapz(r_list*p_t_list,r_list)

    P = torque*omega

    #Dimensionless Coefficient
    Cp = P/(0.5*rho*A*V_0**3)
    CT = thrust/(0.5*rho*A*V_0**2)

    return Cp,CT


