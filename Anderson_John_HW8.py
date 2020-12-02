# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 20:13:34

@author: johna
"""
# MECE 6397, SciComp, HW8 Question 3
#Solve two different schemes for a 1-D wave equation, c=1
#Backward and Forward
#Up to me to decide how far to extend x
#\https://github.com/jeander5/-MECE_6397_HW8_COMP/upload

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
#note this also includes the boundary condtion u(x=0,t)=0    
    lenx=len(x)
    func_vals=np.zeros(lenx)
    k=0
    while x[k]<1:
        func_vals[k]=x[k]*(1-x[k])
        k=k+1
    return (func_vals)
#I dont really like that function. The IC has a conditional statement 
#so I have one too.
#alternatively I could locate index of x that corresponds to x[i]<1
#and then apply the Initial condition up to that index.
#if i can do that withou calling any other functions I should    
#dont worry about it

#Note right biundary cndition is just 0
    
#we will just let length of interval be pi for now.    
L=4*math.pi
N=1000
#wave speed, is implied is simply 1 here
c=1
#starting N kinda big just to get nicer plots

x, dx=DIF(L, N)
#I will define a dt here as well, 
dt = 0.01
C = c*dt/dx
Uo=IC(x)


#
#first lets get the analytic solution,
#simplky from the general solution and the initail condition...
#expanding out the binomials and rearranging...
#-x*x +x+2*x*t-t*t-t for 0<x-t<1
#and 0 for 1<= x-t
#hmm how to code this? Well Im gonna need some conditional statemetns
#but thats okay, if the solution has conditional statements the code does too.
#I will make a function, input will be x domain and a single time t.

def u_exact_func(x,t):
    """Returns the Analytical solution for this Assignment"""
#this is just a right traveling wave, 
#note that its a parabola though not a sine 
#Inputs are a the whole x domain and a single time t.
#so I will just have to call this badboy multiple times, which I am fine with.
#its better this way actually, rather than inputting multiple times.
#I can already see how I will do the graphs this gonna work quite nicely.
#    
    lenx=len(x)
    func_vals=np.zeros(lenx)
    k=0
    while x[k]-t<1:
        func_vals[k] = - x[k]*x[k] + x[k] + 2*x[k]*t - t*t - t
        k=k+1
    return (func_vals)

#u_exact=u_exact_func(x,0)
#fig, ax = plt.subplots()
#plt.grid(1)
#plt.plot(x,u_exact,'r:')
#u_exact=u_exact_func(x,1)
#plt.plot(x,u_exact,'b')
#u_exact=u_exact_func(x,1.5)
#plt.plot(x,u_exact,'g.-')
#u_exact=u_exact_func(x,2)
#plt.plot(x,u_exact,'c-.')
#ax.set(ylim=(0, 1))
##I will append the legend string as needed, for example
#legend_string=(['t= 0s'])
#legend_string.append('t=1s')
#legend_string.append('t=1.5s')
#legend_string.append('t=2s')
#ax.legend(legend_string)
#Very Nice! right traveling wave.
#okay im getting some wierd error that shouldnt be happening if t is not an integer....what is going on....
#oh no its if L-t<1 the wave is travelling outside of my discretized domain.
#I like this

#Okay, scheme 1. BACKWARD
#Im gonna define a constant here
phi=1-C
#and I really so no need to make a matrix here. 
#just a good ol fashioned for loop will suffice
#actually I should make a function.
def Backward_Scheme_func(u_n,C):
    """Returns the u values at the next time step """
#Inputs are teh u values at the previous, nth, time step, for a previously discretized domain
#        
#And the just the Courant number        
#Im gonna define a constant here
    phi=1-C
    N=len(u_n)
#len(u_n) -1 because im not solving for tjat given boundary condition
#okay I will just solve for it
    func_vals=np.ones(N)
#this line is for a boundary condition of u(0,t0)    
    func_vals[0]=u_n[0]
    for j in range(1,N):
        func_vals[j]=u_n[j]*phi+C*u_n[j-1]
    return (func_vals)
    
#ok now its just a matter of how long I wanna carry out the solution for, just call it that many times
#I dont wanna go to function crazy but I wanna make a function like carry out the solution for N times
#Advance Solution Forward_Backward Scheme
def ASF_BS_func(T,dt, u_n):   
    """Moves the Solution forward T/dt number of times """
#inputs for now are total Time, dt, and most up to date u vals
#
#I kinda wanna be able to change everything on the fly, grid, Courant number, ect
#It will just mean more inputs to the function. Thats okay. Im fine with that. 
#Im going for controllablity for this assignment.
#Ill just keep it simple for now  
    N=round(T/dt)
    real_time=N*dt
    #Need to have this real time in here incase user inputs a T that isnt divisible by delta    
    n=len(u_n)
#    func_vals[0]=np.ones(n)
#    okay so this calculates a new Courant number, but still uses the dx from outside the function.
# So the way it is written now i can still modify the courant number    
    C=c*dt/dx
    for j in range(N):
        u_n=Backward_Scheme_func(u_n,C)
    return(u_n)    
            
#example lets compare at 7 seconds        
u_exact=u_exact_func(x,0)
fig, ax = plt.subplots()
plt.grid(1)
plt.plot(x,u_exact,'b')
u_exact=u_exact_func(x,3)
u_7=ASF_BS_func(3,dt,Uo)
plt.plot(x,u_exact,'r')
plt.plot(x,u_7,':g')
ax.set(ylim=(0, 1))
legend_string=(['Exact Solution t= 0s'])
legend_string.append('Exact solution t=3s')
legend_string.append('Backward Scheme t=3s')
#this legend string reports the external courant number, so if I input a new delta t to the adnacing function the legend string will be wrong
ax.title.set_text('Courant Number = %s \n $\Delta$ t= %s s \n $\Delta$ x= %s units'%(round(C,4),round(dt,4),round(dx,4)))
ax.legend(legend_string) 
    
#Okay I see I gotta do something about that left point, no thats just how its gonna be....    
#i dont like these functions they are taking a lot of time to run many points
#oh my dt was just really small. so it just did that inner loop like 12 billions times lol

#I will do the forward scheme tomorrow