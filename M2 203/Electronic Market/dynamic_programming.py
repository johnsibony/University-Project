#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 22:16:24 2019

@author: v
"""

import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt

# Utilities
def h(x, eta):
    return eta * x

def g(x, psi):
    return psi * x 

def C(x,n,lbda,psi,eta,sigma=0.3,tau=0.5):
    res = lbda*n*g(n/tau,psi) + lbda*(x-n)*tau*h(n/tau,eta) + 0.5*(lbda**2)*(sigma**2)*tau*((x-n)**2)
    return res



def dynamic(nb_T, X, lbda, psi, eta, plot='True'):


    ### Initialization
    u = np.zeros(shape=(nb_T, X+1), dtype="float64")      
    b = np.zeros(shape=(nb_T, X+1), dtype="int")          
    inventoryforX = np.zeros(shape=(nb_T,1), dtype="int")     
    inventoryforX[0] = X
    N = []                                                  
    tau = 0.5
    
    for x in range(X+1):
        u[nb_T - 1, x] = np.exp(x * h(x/tau, eta))
        b[nb_T - 1, x] = x
    
    ### Backwards Step
    for t in range(nb_T-2, -1, -1):
        for x in range(X+1):
            
            best_value = u[t+1,0] * np.exp(C(x, x, lbda, psi, eta))
            best_n = x
            
            for n in range(x):
                
                current_value = u[t+1,x-n] * np.exp(C(x, n, lbda, psi, eta))
                
                if current_value < best_value:
                    best_value = current_value
                    best_n = n   
               
            u[t,x] = best_value
            b[t,x] = best_n
    
  
    for t in range(1, nb_T):
        inventoryforX[t] = inventoryforX[t-1] - b[t,inventoryforX[t-1]]
        N.append(b[t,inventoryforX[t-1]])
    
    N = np.asarray(N)
    
    return u, b, inventoryforX, N


'''
First experiment with eta = 0.01, psi = 0.01 for different values of lambda
'''

u1, b1, y1, N1 = dynamic(nb_T=100, X=1000, lbda=0.001, psi=0.01, eta=0.01, plot='False')
u2, b2, y2, N2 = dynamic(nb_T=100, X=1000, lbda=0.01, psi=0.01, eta=0.01, plot='False')
u3, b3, y3, N3 = dynamic(nb_T=100, X=1000, lbda=0.02, psi=0.01, eta=0.01, plot='False')
u4, b4, y4, N4 = dynamic(nb_T=100, X=1000, lbda=0.05, psi=0.01, eta=0.01, plot='False')
u5, b5, y5, N5 = dynamic(nb_T=100, X=1000, lbda=0.1, psi=0.01, eta=0.01, plot='False')



from matplotlib.ticker import MaxNLocator

ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

pas = 8/100
x = [ 9+i* pas for i in range(100)]

plt.plot(x,y1, color='purple', lw=2, label='$\lambda=0.001$')
plt.plot(x,y2, color='blueviolet', lw=2, label='$\lambda=0.01$')
plt.plot(x,y3, color='magenta', lw=2, label='$\lambda=0.02$')
plt.plot(x,y4, color='hotpink', lw=2, label='$\lambda=0.05$')
plt.plot(x,y5, color='pink', lw=2, label='$\lambda=0.1$')

plt.xlabel('Time')
plt.ylabel('Shares remaining in the portfolio')
plt.legend(loc='best')
plt.title('Trading with linear market impact, $\eta=0.01$ and $\psi=0.01$')
plt.savefig('1.png')
plt.show()



'''
First experiment with eta = 0.05, psi = 0.05 for different values of lambda
'''


u1, b1, y1, N1 = dynamic(nb_T=100, X=1000, lbda=0.001, psi=0.05, eta=0.05, plot='False')
u2, b2, y2, N2 = dynamic(nb_T=100, X=1000, lbda=0.01, psi=0.05, eta=0.05, plot='False')
u3, b3, y3, N3 = dynamic(nb_T=100, X=1000, lbda=0.02, psi=0.05, eta=0.05, plot='False')
u4, b4, y4, N4 = dynamic(nb_T=100, X=1000, lbda=0.05, psi=0.05, eta=0.05, plot='False')
u5, b5, y5, N5 = dynamic(nb_T=100, X=1000, lbda=0.1, psi=0.05, eta=0.05, plot='False')




ax = plt.figure().gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

pas = 8/100
x = [ 9+i* pas for i in range(100)]

plt.plot(x,y1, color='purple', lw=2, label='$\lambda=0.001$')
plt.plot(x,y2, color='blueviolet', lw=2, label='$\lambda=0.01$')
plt.plot(x,y3, color='magenta', lw=2, label='$\lambda=0.02$')
plt.plot(x,y4, color='hotpink', lw=2, label='$\lambda=0.05$')
plt.plot(x,y5, color='pink', lw=2, label='$\lambda=0.1$')

plt.xlabel('Time')
plt.ylabel('Shares remaining in the portfolio')
plt.legend(loc='best')
plt.title('Trading with linear market impact, $\eta=0.05$ and $\psi=0.05$')
plt.savefig('2.png')
plt.show()

















