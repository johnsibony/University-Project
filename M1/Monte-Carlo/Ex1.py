#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 16:05:53 2016

@author: root
"""

import numpy as np 
import scipy.stats as sps
import matplotlib.pyplot as plt

#N=1000000
n=10
ntotal=1000000
sigma=0.5
K=1
T=1
delta=T/n

def brownian(): #simulation forward du mouvement brownien
    W=np.zeros(n+1)
    X=np.random.randn(n)/np.sqrt(n)
    W[1:]=np.cumsum(X)
    #t=np.linspace(0,1,n+1)
    #plt.plot(t,W)
    #plt.show()
    return W

def mu(x): #fonction mu
    return (0.1*(np.sqrt(np.exp(x))-1)-1/8)
    
def rec(W): #definition récurrente des Xt
    X=np.zeros(n+1)
    for i in range(n):
        X[i+1]=X[i]+mu(X[i])*delta + sigma*(W[i+1]-W[i])
    #t=np.linspace(0,1,n+1)
    #plt.plot(t,X)
    #plt.show()
    return(X[10])
    
def pricing(x): #fonction g de pricing de l option
    return(np.maximum((np.exp(x)-K),0))
    
def q1(): #Methode de Monte Carlo classique pour le Schéma d Euler
    XT=np.zeros(ntotal)
    for i in range(ntotal):
        W=brownian()
        XT[i]=rec(W)
    F=pricing(XT)
    I=np.mean(F) #estimateur
    s=np.std(F,ddof=1) #variance
    print("Par la méthode de Monte Carlo classique:")
    print("Estimateur :",I)
    print("Variance empirique:",s)
    print("Intervalle de confiance a 95% :",[I-1.96*s/np.sqrt(ntotal),I+1.96*s/np.sqrt(ntotal)])
    print("Erreur relative maximale a 95% :",100*1.96*s/(I*np.sqrt(ntotal)),"%")
    t=range(1,ntotal+1)
    meanList=np.cumsum(F)/t
    CI_high=I-1.96*s/np.sqrt(ntotal) #borne sup de l'intervalle de confiance
    CI_low=I+1.96*s/np.sqrt(ntotal)   #borne inf de l'intervalle de confiance
    plt.plot(t, meanList, "red", label="Empirical Mean")
    plt.axhline(y=CI_high, color='green', hold=True, linewidth=0.5, label="CI higer bound")
    plt.axhline(y=CI_low, color='green', hold=True, linewidth=0.5, label="CI lower bound") 
    plt.xlabel("Number of simulations")
    axes = plt.gca()
    axes.set_ylim(0.19, 0.21)
    plt.legend(loc="best")
    plt.grid(color='gray', linestyle='dashed')
    plt.show()
    
def q2(): #Methode de Monte Carlo avec Variable antithetique
    XT=np.zeros(ntotal)
    XTa=np.zeros(ntotal)
    for i in range(ntotal):
        W=brownian()
        XT[i]=rec(W) #variable standard
        XTa[i]=rec(-W) #variable antithetique
    Fa=pricing(XT)+pricing(XTa) 
    Ia=np.mean(Fa)/2 #estimateur
    sa=np.std(Fa,ddof=1)/2  #variance
    print("Par la méthode de Monte Carlo Antithetique:")
    print("Estimateur :",Ia)
    print("Variance empirique:",sa)
    print("Intervalle de confiance a 95% :",[Ia-1.96*sa/np.sqrt(ntotal),Ia+1.96*sa/np.sqrt(ntotal)])
    print("Erreur relative maximale a 95% :",100*1.96*sa/(Ia*np.sqrt(ntotal)),"%")
    t=range(1,ntotal+1)
    meanList=np.cumsum(Fa/2)/t
    CI_high=Ia-1.96*sa/np.sqrt(ntotal) #borne sup de l'intervalle de confiance
    CI_low=Ia+1.96*sa/np.sqrt(ntotal)   #borne inf de l'intervalle de confiance
    plt.plot(t, meanList, "red", label="Empirical Mean")
    plt.axhline(y=CI_high, color='green', hold=True, linewidth=0.5, label="CI higer bound")
    plt.axhline(y=CI_low, color='green', hold=True, linewidth=0.5, label="CI lower bound")
    plt.xlabel("Number of simulations")
    axes = plt.gca()
    axes.set_ylim(0.19, 0.21)
    plt.legend(loc="best")
    plt.grid(color='gray', linestyle='dashed') 
    plt.show()    

def xTilde(W): #definition de xtilde
    return(sigma*W[10])

def bEst(): #calcul du coefficient b optimal
    k=np.int(ntotal/10)
    X=np.zeros(k)
    Y=np.zeros(k)
    for i in range(k):
        W=brownian()
        X[i]=rec(W)
        Y[i]=xTilde(W)
    X=pricing(X)
    Y=pricing(Y)
    return np.cov(X,Y)[0][1]/np.var(Y) #b optimal
    
def q3(): #Methode de Monte Carlo par Variable de Controle 
    m=np.exp(sigma**2/2)*sps.norm.cdf(-(np.log(K)-sigma**2)/sigma)-K*sps.norm.cdf(-(np.log(K))/sigma) #Calcul de la moyenne de Xtilde
    b=bEst()
    print("b estim",b)
    X=np.zeros(ntotal)
    Xtilde=np.zeros(ntotal)
    for i in range(ntotal): #calcul de Xtilde
        W=brownian()
        X[i]=rec(W)
        Xtilde[i]=xTilde(W)
    X=pricing(X) #calcul de Xt
    Xtilde=pricing(Xtilde)
    Fc=X-b*(Xtilde-m)
    Ic=np.mean(Fc) #estimateur
    sc=np.std(Fc,ddof=1) #variance
    print("Par la méthode de Monte Carlo par Variable de Controle:")
    print("Estimateur :",Ic)
    print("Variance empirique:",sc)
    print("Intervalle de confiance a 95% :",[Ic-1.96*sc/np.sqrt(ntotal),Ic+1.96*sc/np.sqrt(ntotal)])
    print("Erreur relative maximale a 95% :",100*1.96*sc/(Ic*np.sqrt(ntotal)),"%")
    CI_high=Ic-1.96*sc/np.sqrt(ntotal) #borne sup de l'intervalle de confiance
    CI_low=Ic+1.96*sc/np.sqrt(ntotal) #borne inf de l'intervalle de confiance
    t=range(1,ntotal+1)
    meanList=np.cumsum(Fc)/t
    CI_high=Ic-1.96*sc/np.sqrt(ntotal) #borne sup de l'intervalle de confiance
    CI_low=Ic+1.96*sc/np.sqrt(ntotal) #borne inf de l'intervalle de confiance
    plt.plot(t, meanList, "red", label="Empirical Mean")
    plt.axhline(y=CI_high, color='green', hold=True, linewidth=0.5, label="CI higer bound")
    plt.axhline(y=CI_low, color='green', hold=True, linewidth=0.5, label="CI lower bound")
    plt.xlabel("Number of simulations")
    plt.legend(loc="best")
    plt.grid(color='gray', linestyle='dashed')
    axes = plt.gca()
    axes.set_ylim(0.19, 0.21)
    plt.show()

def brownianAccr(W): #calcul des accroissements du mouvement brownien
    X=np.zeros(n)
    for i in range(n):
        X[i]=W[i+1]-W[i]
    return X
    
def fDerivee_i(i,X,teta): #calcul de la derivee de la fonction d'importane
    s=0
    for j in range(n):
        s=s+2*X[j]*teta-teta**2
    return (-1/delta)*(X[i]-teta)*np.exp((-1/(2*delta))*s)*(pricing(recAccr(X)))**2
    #return (-1/delta)*(X[i]-teta[i])*np.exp((-1/(2*delta))*n*np.mean(2*teta*X-teta**2)) #deuxieme methode de calcul
    
def tetaRec_i(i,teta_i,gamma,X): #descente de gradient stochastique pour la ieme coordonnée du parametre
    return teta_i-gamma*fDerivee_i(i,X,teta_i)

def recAccr(X): #definition des accroissements de Xt
    Z=np.zeros(n+1)
    for i in range(n):
        Z[i+1]=Z[i]+mu(Z[i])*delta + sigma*(X[i])
    #t=np.linspace(0,1,n+1)
    #plt.plot(t,X)
    #plt.show()
    return(Z[10])
    
def q4(): #Méthode de Monte Carlo par Fonction d'importance  
    Ii=0
    si=0
    teta=np.zeros(n) #parametre de la fonction d'importance
    for l in range(ntotal):
        W=brownian()
        X=brownianAccr(W)
        Ii= l/(l+1)*Ii + 1/(l+1)*pricing(recAccr(X+teta))*np.exp((-1/(2*delta))*n*np.mean(2*teta*X+teta**2)) #estimateur actualisé
        si= l/(l+1)*si + 1/(l+1)*pricing(recAccr(X+teta))*np.exp((-1/(delta))*n*np.mean(2*teta*X+teta**2)) #variance empirique actualisée
        for i in range(n): #descente de gradient sur teta coordonnée par coordonnée
            teta[i]=tetaRec_i(i,teta[i],1/(100*(l+1)),X)
    print("Par la méthode de Monte Carlo par Fonction d'importance:")
    print("Estimateur :",Ii)
    print("Parametre minimiseur de la famille de densité ",teta)
    print("Variance empirique :",si)
    print("Intervalle de confiance a 95% :",[Ii-1.96*si/np.sqrt(ntotal),Ii+1.96*si/np.sqrt(ntotal)])
    print("Erreur relative maximale a 95% :",100*1.96*si/(Ii*np.sqrt(ntotal)),"%")
    
def q4_graph():
    Ii=np.zeros(ntotal)
    si=np.zeros(ntotal)
    teta=np.zeros(n) #parametre de la fonction d'importance
    for l in range(ntotal-1):
        W=brownian()
        X=brownianAccr(W)
        Ii[l+1]= l/(l+1)*Ii[l] + 1/(l+1)*pricing(recAccr(X+teta))*np.exp((-1/(2*delta))*n*np.mean(2*teta*X+teta**2)) #estimateur actualisé
        si[l+1]= l/(l+1)*si[l] + 1/(l+1)*pricing(recAccr(X+teta))*np.exp((-1/(delta))*n*np.mean(2*teta*X+teta**2)) #variance empirique actualisée
        for i in range(n): #descente de gradient sur teta coordonnée par coordonnée
            teta[i]=tetaRec_i(i,teta[i],1/(100*(l+1)),X)
    print("Par la méthode de Monte Carlo par Fonction d'importance:")
    print("Estimateur :",Ii[ntotal-1])
    print("Variance empirique :",si[ntotal-1])
    print("Intervalle de confiance a 95% :",[Ii-1.96*si/np.sqrt(ntotal),Ii+1.96*si/np.sqrt(ntotal)])
    print("Erreur relative maximale a 95% :",100*1.96*si/(Ii*np.sqrt(ntotal)),"%")
    t=range(1,ntotal+1)
    CI_high=Ii[ntotal-1]-1.96*si[ntotal-1]/np.sqrt(ntotal) #borne sup de l'intervalle de confiance
    CI_low=Ii[ntotal-1]+1.96*si[ntotal-1]/np.sqrt(ntotal)
    plt.plot(t, Ii, "red", label="Empirical Mean")
    plt.axhline(y=CI_high, color='green', hold=True, linewidth=0.5, label="CI higer bound")
    plt.axhline(y=CI_low, color='green', hold=True, linewidth=0.5, label="CI lower bound")
    plt.xlabel("Number of simulations")
    plt.legend(loc="best")
    plt.grid(color='gray', linestyle='dashed')
    axes = plt.gca()
    axes.set_ylim(0.19, 0.21)
    #plt.figure(2)   #graph pour la variance empirique
    #plt.plot(t, si, "red", label="Empirical Variance")
    #plt.xlabel("Number of simulations")
    #plt.legend(loc="best")
    #plt.grid(color='gray', linestyle='dashed')
    plt.show()