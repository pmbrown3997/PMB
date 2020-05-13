# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt

#%% Problem 1

X = []

for ii in range(0,101):
    x = ii
    X.append(x)
    
    
F = []

for jj in range(0,101):
    y = 836773.401*m.exp(-574.528*X[jj])
    F.append(y)
    
plt.figure(1)
plt.grid()
plt.plot(X,F) 
plt.title('Instentaneous Deposition Profile Problem 1')
plt.xlabel('Distance within material[m]')
plt.ylabel('Fluence[J/m^2]')
    
    
#%% Problem 2

hv = [1, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 80, 100]

sigma = [42507, 13463, 5817.5, 1736.9, 726.66, 367.79,79210.71, 88.196, 45.827, 15.693, 8.6561, 5.0678, 4.1258, 3.7238, 3.4928, 3.2093, 3.018]

print('Enter Number between 1 and 100')
y = input()

x = float(y)

xo = 0
yo = 0
x1 = 0
y1 = 0

Sigma = 0

for i in range(0,17):
    s = hv[i]
    if x == s:
        Sigma = sigma[i]
        break
    else:
        Sigma == 0
        
if Sigma == 0:
    for ii in range(0,17):
        s = hv[ii]
        if x > s:
            xo = s
            yo = sigma[ii]
        else:
            jj = ii-1
            xo = hv[jj]
            yo = sigma[jj]
            x1 = s
            y1 = sigma[ii]
            
            Sigma = (yo*(x1-x)+y1*(x-xo))/(x1-xo)
            break
else:
    Sigma = Sigma
        
print('Cross Section = ', Sigma, 'Barns')        
        
#%% Problem 3

a = 7.56*10**(-16)
k = 1.38*10**(-23)
Na = 6.022*10**(-23)

def BE(Z_given):
    BE = (2.517*10**(-18))*Z_given**(7/3)
    return(BE)
    
row = [2710, 7300, 11343] #Al, Fe, Pb

A = [0.027, 0.065, 0.208] #Al, Fe, Pb

Z = [13, 26, 82] #Al, Fe, Pb

def Yield(row_given, A_given, Z_given, m_given):
    Y = row_given*(Na/A_given)*BE(Z_given)*(2+Z_given) + a*(m_given/row_given)*((2/(3*k))*BE(Z_given))**4
    return(Y)
    
    
M = []

for ii in range(0,101):
    M.append(ii)


Yal =[]
Yfe =[]
Ypb =[]

for jj in range(0,101):
    Yal.append(Yield(row[0], A[0], Z[0], M[jj]))    
    Yfe.append(Yield(row[1], A[1], Z[1], M[jj])) 
    Ypb.append(Yield(row[2], A[2], Z[2], M[jj])) 
    
    
plt.figure(2)
plt.grid()
plt.loglog(M,Yal, 'k', label='Al') 
plt.loglog(M,Yfe, 'b', label='Fe') 
plt.loglog(M,Ypb, 'r', label='Pb') 
plt.legend(loc='upper right')
plt.title('Burnout Yield for Cases of Al, Fe, and Pb Problem 3')
plt.xlabel('Mass of Case[Kg]')
plt.ylabel('Yield[J]')    
    
    
YAL = Yield(row[0], A[0], Z[0], M[1])    
YFE = Yield(row[1], A[1], Z[1], M[1])    
YPB = Yield(row[2], A[2], Z[2], M[1])   

YAL = (10/7)*YAL
YFE = (10/7)*YFE
YPB = (10/7)*YPB

Conversion = 2.39*10**(-13)

YAL = Conversion*YAL
YFE = Conversion*YFE
YPB = Conversion*YPB
        
print('The Yield for Al is', YAL, 'kT')
print('The Yield for Fe is', YFE, 'kT')
print('The Yield for Pb is', YPB, 'kT')
 
 
#%% Problem 4 

MI = 0.08999 #For part b

def BuildUp(A1_given, A2_given, c1_given, c2_given, mu_given):
    MFP = mu_given*MI
    BUF = A1_given*m.exp(c1_given*MFP) + A2_given*m.exp(c2_given*MFP)
    return(BUF)
    
print('BUF is =', BuildUp(-114.1, 115.1, 0.14, 0.16, 0.01749))


BUF = [1, 1.000000806, 1.00000521, 1.00000444]
P = [0.055064, 0.250586, 0.376856, 0.245272]

R = 7.07107*10**3 #For Part b

Y = 2.34304*10**14

mu = [-2.21, -0.2392, -0.03782, -0.01749]

Sum = 0

for ii in range(0,4):
    Sum = Sum + P[ii]*BUF[ii]*m.exp(mu[ii]*MI)
    
x = (4*m.pi*R**2)   
        
F = (Y/x)*Sum

print('Fluence is =', F)


















