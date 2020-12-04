# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 20:13:34

@author: johna
"""
# MECE 6397, SciComp, HW8 Question 3
#Solve two different schemes for a 1-D wave equation, c=1
#Backward and Forward
#Up to me to decide how far to extend x
#\https://github.com/jeander5/-MECE_6397_HW8_COMP/

#imports
import math
import numpy as np
import matplotlib.pyplot as plt 

def DIF(L ,N):
    """discretizes an interval with N interior points"""
#Inputs are interval length number of interior points  
#Returns the discretize domain and the interval spacing
#Includes endpoints    
    h = L/(N+1)
    x = np.linspace(0, L, N+2)
    return(x, h)
    
def IC(x):
    """Returns the initial condition for this Homework Problem"""
#note this also includes the boundary condition u(x=0,t)=0    
    lenx=len(x)
    func_vals=np.zeros(lenx)
    k=0
    while x[k]<1:
        func_vals[k]=x[k]*(1-x[k])
        k=k+1
    return (func_vals)
#that really didnt need to be a function

def u_exact_func(x,t):
    """Returns the Analytical solution for this Assignment"""
#Inputs are the whole x domain and a single time t.

    lenx=len(x)
    func_vals=np.zeros(lenx)
    k=0
    while x[k]-t<1:
        func_vals[k] = - x[k]*x[k] + x[k] + 2*x[k]*t - t*t - t
        k=k+1
    return (func_vals)

def Backward_Scheme_func(u_n,C):
    """Returns the u values at the next time step using the backward scheme """
#Inputs are the u values at the previous, nth, time step, for a previously discretized domain           
#Define constant phi outside of the loop.
    phi=1-C
    N=len(u_n)
    func_vals=np.ones(N)
#this line is for a boundary condition of u(0,t0)    
    func_vals[0]=u_n[0]
    for j in range(1,N):
        func_vals[j]=u_n[j]*phi+C*u_n[j-1]
    return (func_vals)
    

#Advance Solution Forward_Backward Scheme
def ASF_BS_func(T,dt, u_n):   
    """Moves the Solution forward T/dt number of times using the backward scheme """
#inputs for now are total Time, dt, and most up to date u vals
    N=round(T/dt)
#this still uses the dx assigned outside of the function    
    C=c*dt/dx
    for j in range(N):
        u_n=Backward_Scheme_func(u_n,C)
    return(u_n)    
            
def Forward_Scheme_func(u_n,C,BC):
    """Returns the u values at the next time step using the forward scheme """
#Inputs are the u values at the previous, nth, time step, for a previously discretized domain
#And also the right boundary condition, u(x=L,t), from the exact function             
#Define constant phi outside of the loop.
    phi=1+C
    N=len(u_n)
    func_vals=np.ones(N)
#this line is for a boundary condition of u(0,t0)    
    func_vals[0]=u_n[0]
    for j in range(1,N-1):
        func_vals[j]=u_n[j]*phi-C*u_n[j+1]    
    func_vals[-1]=u_n[j+1]*phi-C*BC
    # Last equation different, for the right side boundary condition   
    return (func_vals)

def ASF_FS_func(T,dt, u_n):   
    """Moves the Solution forward T/dt number of times using the forward scheme """
#inputs for now are total Time, dt, and most up to date u vals
    N=round(T/dt)
#this still uses the dx assigned outside of the function      
    C=c*dt/dx
#getting Boundary condition from the exact function   
    BC=u_exact_func(x,T)[-1]
    for j in range(N):
        u_n=Forward_Scheme_func(u_n,C,BC)
    return(u_n)    
    
#Call the functions as needed for the report
    
#let length of the interval be 4 pi     
L=4*math.pi
#Number of internal points
N=600
#wave speed, is simply 1 here
c=1
#calling the discretizing the interval function
x, dx=DIF(L, N)
#I will define a dt here as well, 
dt = dx*4
C = c*dt/dx
Uo=IC(x)    
      
#Plotting, uncomment and modify as needed    

#u_exact=u_exact_func(x,0)
#fig2, ax2 = plt.subplots()
#plt.grid(1)
#plt.plot(x,u_exact,'b')
#u_exact=u_exact_func(x,5.368)
##u_7=ASF_BS_func(3,dt,Uo)
#u_6=ASF_FS_func(5.368,dt,Uo)
#plt.plot(x,u_exact,'r')
#plt.plot(x,u_6,'-.c')
#ax2.set(ylim=(0, 0.5))
#ax2.set(xlim=(0, 12))
##u_7=ASF_BS_func(5,dt,Uo)
##plt.plot(x,u_7,'m--')
##I will just append a legend string as I go along
#legend_string=([])
#legend_string=(['Exact Solution t= 0s'])
#legend_string.append('Exact solution t=5.368 s')
#legend_string.append('Forward Scheme t=5.368s')
##legend_string.append('Backward Scheme t=5.75s')
##legend_string.append('Backward Scheme t=5s')
##this legend string reports the external courant number, 
##so if I input a new delta t to the advancing function the legend string will be wrong
#ax2.title.set_text('Courant Number = %s \n $\Delta$ t= %s s \n $\Delta$ x= %s units'
#                   %(round(C,4),round(dt,4),round(dx,4)))
#ax2.legend(legend_string)
#plt.xlabel('x')
#plt.ylabel('u(x,t)')    
