# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
import psapy.Vogel as IPR
import psapy.FluidProps as FluidProps
import psapy.FluidProps
import psapy.BeggsandBrill as BB
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
sns.set(color_codes=True)
from ipywidgets import *


# ko = float(input("Enter perm(Darcy)"))
# h = float(input("Enter pay zone thickness(ft)"))
# muo = float(input("Enter the viscosity(cP):"))
# Bo = float(input("Enter formation volume factor:"))
# s = float(input("Enter the skin factor:"))
# Pr = float(input("Enter the reservoir pressure(psi):"))
# Pb = float(input("Enter the buble point pressure(psi):"))
# ct = float(input("Enter the total compressibility(psi^-1):"))
# A = float(input("Enter the drainage area(acres):"))
# re=np.sqrt(43560*A/np.pi)
# rw = float(input("Enter the wellbore radius(ft):"))

# Psat = FluidProps.Pbub(Temp,75,100,GasGrav, API, GOR)
# Psat


# J= (0.00708*ko*h)/(muo*Bo*(np.log(re/rw)-0.75+s))

# qmax= J*Pr/1.8


def ipr(Pr,Pb,ptest,qtest,re=1000,rw=0.25,h=100,ko=100,mu=100, Bo=1):  
    # if ptest>=Pb:
    #     j=qtest/(Pr-ptest)    
    # else:
    #     j=qtest/(Pr-Pb+Pb/1.8*(1-0.2*ptest/Pr-0.8*(ptest/Pr)**2))    
    # print("Productivity index is {}".format(j))
    
    if Pb>=Pr:
       list_pwf2=[]
       ylim=Pr  
       list_q=[]
       list_pwf2=[]
       j_darcy=0.00708*ko*h/(mu*Bo*np.log(0.472*re/rw))
       print("Productivity index(J) is {} STB/D-psi".format(round(j_darcy,4)))
       qob=0
       for i in range(round(Pr,-2),-100,-100):
           qomax=qtest/(1-0.2*ptest/Pr-0.8*(ptest/Pr)**2)
           qo=qomax*(1-0.2*(i/Pr)-0.8*(i/Pr)**2)
           list_pwf2.append(i)
           list_q.append(round(qo))
       print(list_q,list_pwf2)    
    else:   
        #The case for undersaturated reservoirs where Pb is less than Pr
        ylim=Pr    
        list_pwf2=[]
        list_q=[]
        for i in range(round(Pr,-1),-100,-100):
            list_pwf2.append(i)
            
            if ptest>=Pb:
                #When the entered test pressure is below bubble point
                j1=qtest/(Pr-ptest)       
                qob=j1*(Pr-Pb)
                
                if i>=Pb:
                    #When the loop pressure value of i is greater than bubble point                    
                    qo=j1*(Pr-i)
                    list_q.append(round(qo))
                else:
                    #When the loop pressure value of i is less than bubble point
                    qo=qob+(j1*Pb*(1-0.2*i/Pb-0.8*(i/Pb)**2)/1.8)
                    list_q.append(round(qo))
                    
            else:
                #When the entered test pressure is above bubble point
                j1=((qtest/((Pr-Pb)+Pb/1.8*(1-0.2*ptest/Pr-0.8*(ptest/Pr)**2))))  
                qob=((qtest/((Pr-Pb)+Pb/1.8*(1-0.2*ptest/Pr-0.8*(ptest/Pr)**2))*(Pr-Pb)))
                if i>=Pb:
                    qo=j1*(Pr-i)
                    list_q.append(round(qo))
                else:
                    qo=qob+(j1*Pb/1.8*(1-0.2*i/Pb-0.8*(i/Pb)**2))
                    list_q.append(round(qo))  
    return list_q,list_pwf2










