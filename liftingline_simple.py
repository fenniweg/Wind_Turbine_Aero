
import numpy as np
import matplotlib.pyplot as plt
import bem



# ============================================================
# Wing geometry
# ============================================================

b = 20
chord = 1.2


# ============================================================
# Basic lifting line script that can compute the induced wind
# for a known circulation distribution
# ============================================================


# ============================================================
# Setting up the influence matrix
#
# Defining wing and placing points on the blade using a
# cosine distribution
# ============================================================

npoint = 50

delta_theta = np.pi / (npoint - 1)

theta = np.zeros(npoint)
yp = np.zeros(npoint)
xp = np.zeros(npoint)
zp = np.zeros(npoint)


for i in range(npoint):
    theta[i] = i * delta_theta
    yp[i] = -0.5 * b * np.cos(theta[i])
    xp[i] = 0.0
    zp[i] = 0.0



# ============================================================
# Calculating pivot points, where the induction is calculated
# ============================================================

ypiv = np.zeros(npoint - 1)
xpiv = np.zeros(npoint - 1)
zpiv = np.zeros(npoint - 1)

for i in range(npoint - 1):
    theta_mean = 0.5 * (theta[i] + theta[i + 1])
    ypiv[i] = -0.5 * b * np.cos(theta_mean)
    xpiv[i] = 0.0
    zpiv[i] = 0.0



# Building up influence matrix by taking one pivot point
# at a time
# ============================================================

influence = np.zeros((npoint - 1, npoint - 1))

for i in range(npoint - 1):

    # Pivot point
    Px = xpiv[i]
    Py = ypiv[i]
    Pz = zpiv[i]

    # Finding the influence from all horseshoe vortices
    # at this pivot point for a unit circulation in each
    # horseshoe vortex

    for k in range(npoint - 1):

        # A unit circulation
        gamma_unit = 1.0

        # Analytical Biot-Savart for a half-infinite
        # straight vortex line

        indu1 = 0.0

        h = yp[k + 1] - Py
        indu1 = gamma_unit / (4.0 * np.pi * h)

        h = Py - yp[k]
        indu2 = gamma_unit / (4.0 * np.pi * h)

        w = indu1 + indu2

        influence[i, k] = w


# ============================================================
# Now the influence matrix has been created and can be used
# ============================================================

alpha_g = 2.0          # deg


# Speed and density
vo = 100.0             # m/s
rho = 1.225            # kg/m^3


# ============================================================
# Initialization (first iteration)
#
# Putting the induced wind (downwash) to zero in every
# pivot point
# ============================================================

w = np.zeros(npoint - 1)
alpha_eff = np.zeros(npoint - 1)
alpha_i = np.zeros(npoint - 1)
cl = np.zeros(npoint - 1)
v_rel = np.zeros(npoint - 1)
gamma = np.zeros(npoint - 1)
w_star = np.zeros(npoint - 1)
R = np.zeros(npoint - 1)
l = np.zeros(npoint - 1)
d = np.zeros(npoint - 1)



# ============================================================
# Iteration starts
# (beta is a relaxation factor)
# ============================================================

beta = 0.01


for n in range(1000):
    # ============================================================
    # For every section on the wing
    # ============================================================
    for k in range(npoint - 1):
        
        v_rel[k] = np.sqrt(vo**2 + w[k]**2) # Relative velocity at the pivot point
        alpha_i[k] = np.arctan(w[k]/ vo)  # Induced angle of attack at the pivot point
        alpha_eff[k] = alpha_g -np.degrees( alpha_i[k] )# Effective angle of attack at the pivot point

        cl[k] = 0.11*alpha_eff[k]# Lift coefficient at the pivot point
        R[k] = 0.5*rho*v_rel[k]**2 * cl[k] # Dynamic pressure at the pivot point
        l = R[k] * np.cos(alpha_i[k]) # Lift per unit span at the pivot point
        d[k] = np.tan(alpha_i[k]) * l # Induced drag per unit span at the pivot point
        

        gamma[k] = 0.5*cl[k] * v_rel[k] * chord  # Circulation at the pivot point
     


    w_star= influence @ gamma #Induced wind at every section due to the circulation distribution

    w_new = beta*w_star+(1-beta)*w # Relaxation of the induced wind at the pivot point

    #check for convergence
    if n % 100 == 0:
        print('Iteration', n, 'Max change in induced wind =', np.max(np.abs(w-w_new)))

    if (np.max(np.abs(w-w_new)))< 1e-6:
        print('Converged after', n, 'iterations')
        break
    w = w_new

#analytical solution for an elliptical circulation distribution
#plot lift, induced drag, induced wind(downwash) and effective angle of attack distribution along the span of the wing
#4 subplots

# plt.figure(figsize=(12, 8))
# plt.subplot(2, 2, 1)
# plt.plot(ypiv, cl, 'b-', lw=2.5, label='Lift coefficient distribution')
# plt.xlabel('Spanwise location $y$ [m]')
# plt.ylabel('Lift coefficient $C_l$')
# plt.title('Lift coefficient distribution along the span')
# plt.grid(True, alpha=0.3)   

# plt.subplot(2, 2, 3)
# plt.plot(ypiv, w, 'r-', lw=2.5, label='Induced wind (downwash) distribution')
# plt.xlabel('Spanwise location $y$ [m]')
# plt.ylabel('Induced wind $w$ [m/s]')
# plt.title('Induced wind (downwash) distribution along the span')
# plt.grid(True, alpha=0.3)   

# plt.subplot(2, 2, 4)
# plt.plot(ypiv, alpha_eff, 'g-', lw=2.5, label='Effective angle of attack distribution')
# plt.xlabel('Spanwise location $y$ [m]')
# plt.ylabel('Effective angle of attack $\\alpha_{eff}$ [deg]')
# plt.title('Effective angle of attack distribution along the span')
# plt.grid(True, alpha=0.3)

# plt.subplot(2, 2, 2)
# plt.plot(ypiv,d , 'm-', lw=2.5, label=' induced drag')
# plt.xlabel('Spanwise location $y$ [m]')
# plt.ylabel('induced drag $d$ [N/m]')
# plt.title('Induced drag distribution along the span')
# plt.grid(True, alpha=0.3)   

# plt.show()  
