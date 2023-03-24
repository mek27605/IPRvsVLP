# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 12:47:43 2023

@author: emink
"""

import fluidproperties as fpo
import BeggsandBrill as BBe
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define input parameters

oil_rate = 100  # Oil Flow Rate (STB/D)
water_rate = 50  # Water Flow Rate (STB/D)
WC = water_rate / (oil_rate + water_rate)  # Water Cut (Unitless)
GOR = 300  # Gas-Oil Ratio (SCF/STB)
gas_grav = 0.65  # Specific Gravity of Gas (Unitless)
oil_grav = 35    # API Gravity of Oil (API)
wtr_grav = 1.07  # Water Specific Gravity (Unitless)
diameter = 2.441  # Tubing Inside Diameter (inches)
angle = 90.0   # Well Deviation Angle (Degrees)
thp = 150.0  # Tubing Head Pressure (psia)
tht = 100.0  # Tubing Head Temperature (Fahrenheit)
twf = 150.0  # Bottomhole Temperature (Fahrenheit)
depth = 5000  # Well Depth (feet)
roughness = 0.0006 # Tubing roughness (inches)
Psep = 114.7  # Separator pressure (psia)
Tsep = 50  # Separator Temperature (Fahrenheit)


# Define a function to calculate temperature gradient
def temperature_gradient(T1, T2, Depth):
    if Depth == 0:
        return 0
    else:
        return abs(T1 - T2) / Depth
    

    
# Calculate temperature gradient and create temperature array
t_grad = temperature_gradient(tht, twf, depth)  # Gradient in °F/ft
depths = np.linspace(0, depth, 50)  # Depth array with 50 members
temps = tht + t_grad * depths  # Temperature array with 50 members


# Define a function to calculate pressure traverse
def pressure_traverse(oil_rate):
    """This function calculates the pressure traverse which represents the pressure distiribution      
       throughout the wells."""
    pressure_list = []
    dpdz_list = []

    for i in range(len(depths)):
        if i == 0:
            pressure_list.append(thp)
        else:
            dz = (depths[i] - depths[i-1])
            pressure = pressure_list[i-1] + dz * dpdz_list[i-1]
            pressure_list.append(round(pressure, 3))

        dpdz_step = BBe.beggsandbrill(pressure_list[i], temps[i], oil_rate, WC, GOR, gas_grav, oil_grav, wtr_grav, diameter, angle, roughness, Psep,Tsep)
        dpdz_list.append(round(dpdz_step, 5))

    return pressure_list, dpdz_list

p, dpdz =pressure_traverse(oil_rate)

def VLP(rates):
    """This function calculates the bottomhole pressure for each flow rate."""
    
    bhp_list =[]
    for q in rates:
        p, dpdz= pressure_traverse(q)
        bhp_step = p[-1]
        bhp_list.append(round(bhp_step,3))
    return bhp_list

rate_array = np.linspace(100, 5000, 50)  #Creates a rate array with 50 members btw 100 to 5000 STB/D
bhp_array= VLP(rate_array)    #Runs the VLP function

ratevspress=pd.DataFrame(data={"Flow Rate":rate_array,"BHP":bhp_array})
ratevspress.head(10)

plt.plot(rate_array, bhp_array)
plt.xlabel('Liquid Rate')
plt.ylabel('Bottomhole Pressure')
plt.title('VLP Curve')
plt.show()
