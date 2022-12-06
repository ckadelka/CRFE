#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 15:32:28 2021

@author: ckadelka
"""

import numpy as np
import scipy 
import matplotlib.pyplot as plt



threshold = 1.5
alpha=0.1
beta=0.25
belief = 10
sigma=1

k=5

distribution  = scipy.stats.norm


mu_alpha = scipy.optimize.bisect(lambda x: distribution.cdf(threshold,x,sigma)-(1-alpha),-5*sigma+threshold,5*sigma+threshold)
betas=np.zeros(k)
dummy_sum = sum([k-j+(j-1)*belief for j in range(1,k+1)])
mus = np.zeros(k)
for i in range(k):
    betas[i] = beta*(k*(k-(i+1)+i*belief))/dummy_sum
    mus[i] = scipy.optimize.bisect(lambda x: distribution.cdf(threshold,x,sigma)-betas[i],-5*sigma+threshold,5*sigma+threshold)
    
    
    
#PLOTTING
x = np.linspace(mu_alpha-sigma*3,max(mus)+sigma*3,200)
  
f,ax=plt.subplots()  
ax.plot(x,distribution.pdf(x,mu_alpha,sigma))
for i in range(k):
    ax.plot(x,distribution.pdf(x,mus[i],sigma))
ax.plot([threshold,threshold],ax.get_ylim(),'k-',lw=2)

    
f,ax=plt.subplots()  
ax.plot(x,distribution.cdf(x,mu_alpha,sigma))
for i in range(k):
    ax.plot(x,distribution.cdf(x,mus[i],sigma))
ax.plot([threshold,threshold],ax.get_ylim(),'k-',lw=2)    




all_mus = np.append(mus,mu_alpha)
n_genes_annotated = np.append([100]*len(mus),8000)

f,ax=plt.subplots()  
ax.plot(x,distribution.pdf(x,mu_alpha,sigma))
ax.plot(x, np.dot(n_genes_annotated[:k],[distribution.pdf(x,all_mus[i],sigma) for i in range(k)])/sum(n_genes_annotated[:k]) )
ax.plot([threshold,threshold],ax.get_ylim(),'k-',lw=2)

