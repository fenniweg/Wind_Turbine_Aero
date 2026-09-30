'''
This file loads the necessary data for the wind turbine analysis, including airfoil data, blade data, and turbine parameters. 
The data is used in subsequent modules for calculations related to the performance of the wind turbine.
Data for the wind tubrine is defined here, including rotor radius, number of blades, rated power, cut-in and cut-out wind speeds, and air density.
'''

import glob
import os

import numpy as np
import pandas as pd

#lookup global file path for data folder
data_path = glob.glob('**/data', recursive=True)[0]


#Load blade data and name columns for easier access using the global data path
blade_dat = pd.read_fwf(os.path.join(data_path, 'bladedat.txt'), header=None)
blade_dat.columns =['r','c','beta','t/c']



#load airfoil data and store in 3D array
airfoil_data = np.zeros((6, 105,4))  # 6 airfoils,105 data points each,4 columns (alpha, C_l(alpha), C_d(alpha), C_m(alpha))

#each airfoil has one thickness/chord ratio t/c  not imolemned her but directly in bem code, can be added 

airfoil_data[0] = pd.read_csv(os.path.join(data_path, 'cylinder.txt'), header=None, sep=None, engine='python')
airfoil_data[1] = pd.read_csv(os.path.join(data_path, 'FFA-W3-600.txt'), header=None, sep=None, engine='python')
airfoil_data[2] = pd.read_csv(os.path.join(data_path, 'FFA-W3-480.txt'), header=None, sep=None, engine='python')
airfoil_data[3] = pd.read_csv(os.path.join(data_path, 'FFA-W3-360.txt'), header=None, sep=None, engine='python')
airfoil_data[4] = pd.read_csv(os.path.join(data_path, 'FFA-W3-301.txt'), header=None, sep=None, engine='python')
airfoil_data[5] = pd.read_csv(os.path.join(data_path, 'FFA-W3-241.txt'), header=None, sep=None, engine='python')


test = pd.read_fwf(os.path.join(data_path, 'bladedat.txt'), header=None)
#Data for 10 MW Turbine

R = 89.17 #Rotor radius in meters
n_blades = 3 #number of blades
P_rated = 10.64e6 #rated power in watts
v_min = 4 #cut in wind speed
v_max = 25 #cut out wind speed
rho = 1.225 #air density in kg/m^3

A = np.pi*R**2 #rotor area in m^2

V_rated= 11.19

#Define range of tip speed ratios (lambda) and pitch angles (theta_p) to test
tip_speed_ratio= np.linspace(5.0, 10.0, 15)      # Tip speed ratio range 
theta_p = np.linspace(-4.0, 3.0, 15)       # Pitch angle range [deg] 


#import data from ashes
v_5_ashes = pd.read_csv(os.path.join(data_path, 'ashes_5.txt'), header=None, sep=None, engine='python')
v_9_ashes = pd.read_csv(os.path.join(data_path, 'ashes_9.txt'), header=None, sep=None, engine='python')
v_11_ashes = pd.read_csv(os.path.join(data_path, 'ashes_11.txt'), header=None, sep=None, engine='python')
v_20_ashes = pd.read_csv(os.path.join(data_path, 'ashes_20.txt'), header=None, sep=None, engine='python')
r_ashes = pd.read_csv(os.path.join(data_path, 'blade_ashes.txt'), header=None, sep=None, engine='python')

v_5_ashes = v_5_ashes.to_numpy().flatten()
v_9_ashes = v_9_ashes.to_numpy().flatten()
v_11_ashes = v_11_ashes.to_numpy().flatten()
v_20_ashes = v_20_ashes.to_numpy().flatten()
r_ashes = r_ashes.to_numpy().flatten()  

r_ashes = r_ashes/ 89.17 #normalize radius to rotor radius

