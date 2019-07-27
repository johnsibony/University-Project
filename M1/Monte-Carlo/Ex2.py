#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 12:25:37 2016

@author: root
"""

import numpy as np
#import scipy.stats as sps
import matplotlib.pyplot as plt
import scipy.misc as msc

n=10
T=1
N=1000000
nbsim=1000000 # Pour l'Allocation Optimale
beta=0.1
sigma=0.5
K=1

def randGrid(): #simulation de la grille aléatoire
    t=[0]       #initialisation T_0=0
    i=0
    while(t[len(t)-1]<T):
        t.append(np.minimum(t[len(t)-1]-np.log(np.random.rand())/beta,T))
        i=i+1
    #plt.hist(NT, normed=1,histtype='bar', bins=np.int(N**(1/3)))
    #plt.show()
    return t

def brownianForward(t): #simulation forward du mouvement brownien sur une grille t donnée
    W=np.zeros(len(t))
    X=np.random.randn(len(t)-1)*np.sqrt(accr(t))
    W[1:]=np.cumsum(X)
    return W
    #t=np.linspace(0,1,n+1)
    #plt.plot(t,W)
    #plt.show()

def accr(W):    #calcul des accroissements 
    X=np.zeros(len(W)-1)
    for i in range(len(X)):
        X[i]=W[i+1]-W[i]
    return X

def mu(x):  #fonction mu
    return (0.1*(np.sqrt(np.exp(x))-1)-1/8)
    
def rec(t,W):   #definition récursive des Xt
    k=len(t)
    X=np.zeros(k)
    accr_t,accr_W=accr(t),accr(W)
    for i in range(k-1):
        X[i+1]=X[i]+mu(X[i])*accr_t[i]+sigma*accr_W[i]
    return X
    
def pricing(X): #fonction g de pricing de l'option
    return(np.maximum((np.exp(X)-K),0)) 
    
def funcPsi(t,W): #fonction Psi
    p=1
    NT=np.int(len(t)-2)
    X=rec(t,W)
    accr_t,accr_W=accr(t),accr(W) #calcul des accroissements du mouvements brownien et de la grille 
    for i in range(NT):
        p=p*(mu(X[i+1])-mu(X[i]))*accr_W[i+1]/(sigma*beta*accr_t[i+1]) #calcul à part du produit
    return(np.exp(beta*T)*(pricing(X[len(X)-1])-pricing(X[len(X)-2]))*p)

def q1(): #Methode de monte carlo classique sur la grille aléatoire
    psi=np.zeros(N)
    for i in range (N):
        t=randGrid()    #simulation de la grille aléatoire
        W=brownianForward(t)
        psi[i]=funcPsi(t,W)
    I=np.mean(psi) #estimateur
    s=np.std(psi,ddof=1) #variance de l'estimateur
    print("Par la méthode de Monte Carlo classique:")
    print("Estimateur :",I)
    print("Variance empirique :",s)
    print("Intervalle de confiance a 95% :",[I-1.96*s/np.sqrt(N),I+1.96*s/np.sqrt(N)])
    print("Erreur relative maximale a 95% :",100*1.96*s/(I*np.sqrt(N)),"%")
    t=range(1,N+1)
    meanList=np.cumsum(psi)/t
    CI_high=I+1.96*s/np.sqrt(N) #borne sup de l'intervalle de confiance
    CI_low=I-1.96*s/np.sqrt(N)   #borne inf de l'intervalle de confiance
    plt.plot(t, meanList, "red", label="Empirical Mean")
    plt.axhline(y=CI_high, color='green', hold=True, linewidth=0.5, label="CI higer bound")
    plt.axhline(y=CI_low, color='green', hold=True, linewidth=0.5, label="CI lower bound") 
    plt.xlabel("Number of simulations")
    plt.legend(loc="best")
    plt.grid(color='gray', linestyle='dashed')
    axes = plt.gca()
    axes.set_ylim(0.18, 0.22)
    plt.show()
    
def probaNT(k): #probabilité que NT = k
    return np.exp(-beta*T)*(beta*T)**k/msc.factorial(k)

def strates(): #definition du nombre des strates nécessaire par allocation proportionnelle 
    s=[np.int(np.round(probaNT(0)*N))]
    i=1
    k=np.round(probaNT(i)*N)
    while(k>0): #la strate NT=k est admissible tant que la partie entiere du nombre de tirage dans la strat est non nulle
        s.append(np.int(k))
        i=i+1
        k=np.round(probaNT(i)*N)
    print('nombres de strates: ',len(s))
    return s
    #print(s)
    #print(np.sum(s))

def temps(k): #simulation des Ti sachant NT=k 
    t=np.zeros(k+2)
    t[1:k+1]=np.sort(np.random.rand(k)) #réarangement croissant de k uniformes 
    t[k+1]=1
    #print(t)
    return(t)

def q2_1(): #Methode de Monte Carlo par stratification avec Allocation proportionnelle
    S=strates() #definition du nombre de strates 
    psi=np.zeros(len(S))
    vstrat=np.zeros(len(S)) 
    pk=np.asarray(range(len(S)))
    pk=probaNT(pk) #definition des probabilités de tomber dans chaque strate
    W=np.random.randn(S[0])
    for i in range(len(S)): #calcul de l esperance conditionnnelle par strate
        psi_k=np.zeros(S[i])
        for j in range(S[i]):
            t=temps(i)
            W=brownianForward(t)
            psi_k[j]=funcPsi(t,W) #simulation de psi sachant NT=i
        vstrat[i]=np.std(psi_k, ddof=1)**2 #calcul de la variance empirique par strate
        psi[i]=np.mean(psi_k) #calcul de l esperance conditionnelle empirique dans la strat i 
    Is=np.sum(psi*pk)
    stotal=np.sqrt(np.mean(vstrat)) #calcul pondéré par pk de la variance totale
    print("Par la méthode de Monte Carlo par stratification proportionnelle:")
    print("Estimateur :",Is)
    print("Variance total :",stotal)
    print("Intervalle de confiance a 95% :",[Is-1.96*stotal/np.sqrt(N),Is+1.96*stotal/np.sqrt(N)])
    print("Erreur relative maximale a 95% :",100*1.96*stotal/(Is*np.sqrt(N)),"%")
        
def q2_2():  #Methode de Monte Carlo par stratification avec Allocation Optimale
    psi=[]                                    
    N=[]    
    k=0
    a=np.exp(-beta*T)
    p=[]
    while a*nbsim >= 1 : #definition des strates
       p.append(a)
       N.append(k)
       k=k+1
       a=(((beta*T)**k)/msc.factorial(k))*np.exp(-beta*T)
    vstrat=[]        #vecteur de variance empirique par strate
    n=10000
    for i in range(k):  #première estimation pour obtenir une variance empirique dans chaque strate et définir l'allocation optimale
       ksi=[]
       for l in range(n):
           U=np.random.rand(N[i])*T
           for j in range(N[i]):
               for r in range(j,N[i]):
                   if U[r]<U[j] : 
                       tmp=U[r]
                       U[r]=U[j]
                       U[j]=tmp
           U=np.append(U,T)
           U=np.insert(U, 0, 0)
           X=[0]
           Z=np.random.randn(N[i]+1)
           for q in range(0,N[i]+1):
               Z[q] = (Z[q]*np.sqrt(U[q+1]-U[q]))
           W=np.zeros(N[i]+2)
           W[1:]=np.cumsum(Z)
           for j in range(1,N[i]+2) :
               X.append(X[j-1] + mu(X[j-1])*(U[j]-U[j-1])+sigma*(W[j]-W[j-1]))
           P=1
           for j in range(1,N[i]+1) :
               P=P*(((mu(X[j])-mu(X[j-1]))*(W[j+1]-W[j]))/(sigma*beta*(U[j+1]-U[j])))
           if N[i]>0:
               tmp = 1
           else:
               tmp = 0
           ksi.append(np.exp(beta*T)*(pricing(X[N[i]+1])-pricing(X[N[i]])*tmp)*P)
       V=np.mean(ksi)
       S=[]
       for l in range(n):
           S.append((ksi[l]-V)**2)
       Var=(float(1)/float(n-1))*np.sum(S)
       vstrat.append(Var)       #variance empirique par strate
    Vartot = []
    pp = []    #allocation optimale
    for i in range(len(p)):                 
       pp.append(p[i]*np.sqrt(vstrat[i]))
    s=sum(pp)
    for i in range(len(p)):
       pp[i]=pp[i]/s
    for i in range(N[k-1]):
       ksi=[]
       for l in range(int(nbsim*pp[i])):   #allocation optimale
           U=np.random.rand(N[i])*T    
           for j in range(N[i]):     #tri des uniformes pour obtenir la subdivision
               for r in range(j,N[i]):
                   if U[r]<U[j] : 
                       tmp = U[r]
                       U[r] = U[j]
                       U[j] = tmp
           U=np.append(U, T)
           U=np.insert(U, 0, 0)
           X=[0]
           Z=np.random.randn(N[i]+1)
           for q in range(0,N[i]+1):
               Z[q]=(Z[q]*np.sqrt(U[q+1] - U[q]))
           W=np.zeros(N[i]+2)
           W[1:]=np.cumsum(Z)
           for j in range(1,N[i]+2) :
               X.append(X[j-1]+mu(X[j-1])*(U[j]-U[j-1])+sigma*(W[j]-W[j-1]))
           P=1
           for j in range(1,N[i]+1) :
               P=P*(((mu(X[j])-mu(X[j-1]))*(W[j+1]-W[j]))/(sigma*beta*(U[j+1]-U[j])))
           if N[i]>0:
               tmp = 1
           else:
               tmp = 0
           ksi.append(np.exp(beta*T)*(pricing(X[N[i]+1])-pricing(X[N[i]])*tmp)*P)
       V=np.mean(ksi)    #moyenne dans chaque strate
       S=[]
       for l in range(int(nbsim*pp[i])):
           S.append((ksi[l] - V)**2)
       Var=(float(1)/float((nbsim*pp[i])-1))*np.sum(S)    #variance dans chaque strate
       Vartot.append(Var)
       psi.append(V)
    M=0   #estimateur
    for i in range (N[k-1]):
       M=M+(psi[i]*p[i])
    VT=0
    for i in range(N[k-1]):
       VT=VT+p[i]*np.sqrt(Vartot[i])
    VT=(VT)**2      #variance totale
    ICmoins=M-((1.96*VT)/np.sqrt(nbsim))
    ICplus=M+((1.96*VT)/np.sqrt(nbsim))
    print('Par la méthode de Monte Carlo par stratification avec Allocation Optimale:')
    print('Estimateur :', M)
    print('Variance :', VT)
    print('Les bornes de l intervalle de confiance sont :', ICmoins)
    print('et ', ICplus)  
    print('La variance totale est ',Vartot)
    print('La probabilité de tomber dans chaque strates est : ', p)