#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:27:00 2020

@author: JaylenJames

Script will take the recreation of CO's inputs from paper and create 
    a 3D histogram of the data.


This file includes changes to histogram4.py after plotting CO's input correctly
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from sklearn import preprocessing
import pandas as pd
from IPython import get_ipython
from stereogram_function import stereogram
import seaborn as sns
from scipy import stats
from scipy import integrate
from sklearn.mixture import GaussianMixture
from matplotlib.colors import LogNorm

#Setting to show plots in a separte window. Set to 'inline' to show in console
#get_ipython().run_line_magic('matplotlib', 'qt')


#Recreating CO's inputs


double_int_results = np.zeros([10,10])

alpha = range(1,11)
beta = range(1,11)
B = 1
s=10
k=10
delta = B/s
index_row = 0

for beta_iteration in range(10):
    
    index_col = 0
    for alpha_iteration in range(10):
        
        f = lambda y, x: (90/19)*x*np.sqrt(1-x**2)*(1-y)*((1-x**2)+(1-y))
        
        y_low = (beta[beta_iteration]-1)/k
        y_high = beta[beta_iteration]/k
        
        x_low = (alpha[alpha_iteration]-1)*delta
        x_high = alpha[alpha_iteration]*delta
        
        result, error = integrate.dblquad(f, y_low , y_high , lambda x: x_low, lambda x: x_high );
        
        double_int_results[index_row][index_col] = result
        
        index_col += 1
    
    index_row +=1 








xedges = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
yedges = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]



if __name__ == '__main__':
    # Create a figure for plotting the data as a 3D histogram.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    #
    #hist, xedges, yedges = np.histogram2d(min_semi_axis_norm, y_square_norm, bins=(xedges, yedges))
    
    
    
    # Construct arrays for the anchor positions of the 16 bars.
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1], indexing="ij")
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    zpos = 0
    
    # Construct arrays with the dimensions for the 16 bars.
    dx = 0.1 * np.ones_like(zpos)
    dy = 0.1 * np.ones_like(zpos)
    dz = double_int_results.ravel()
    
    
    bar_color = 'r' #np.array([1,1,2,2,3,3])
    bar_color_set = ['r','r','b','b','g','g']
    bar_color_array = np.asarray(bar_color_set)
    
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color = None ,zsort='average')
    
    plt.show()

P, Q = stereogram() #Obtain values for 3D spheroids histogram



####### Begin Calculation of 3D speroid size and shape frequencies ##########

#Set constants to be used in calculation
#H_bar = np.mean(min_semi_axis_norm)#1           #replace with correct value
H_bar = 298/(65*np.pi) #Estimate through trial and error in this code
#H_bar = 304/(135*np.pi) #value from CO's Part 1 paper

#delta = 1/s

g = np.zeros([10,10])           #initialize blank g matrix
val_mat = np.zeros([10,10])
    
i_index_row = 0

for i in range(1, 11):  #cycle through row elements of g matrix
    j_index_col = 0
    
    for j in range(1, 11):   #cycle through column elements of g matrix
        alpha_index_row = 0
        
        for alpha in range(1,11): 
            beta_index_col = 0
            
            for beta in range(1,11):
                
                val_mat[alpha_index_row][beta_index_col] = P[i_index_row][alpha_index_row]*double_int_results[alpha_index_row][beta_index_col]*Q[beta_index_col][j_index_col];
            
                beta_index_col += 1
            
            alpha_index_row += 1
            
        val_mat_sum = np.sum(val_mat)
        mltd_val_mat = (H_bar/delta)*val_mat_sum
        g[i_index_row][j_index_col] = mltd_val_mat
        
        j_index_col += 1
        
    i_index_row += 1


g_ij_sum = np.sum(g)


if __name__ == '__main__':
    # Create a figure for plotting the data as a 3D histogram.
    ######Add axis lables to the graphs!!!#######
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    #xpos = [range(g.shape[0])]
    #ypos = [range(g.shape[1])]
    #xpos, ypos = np.meshgrid(xpos, ypos)
    #xpos = xpos.flatten('F')
    #ypos = ypos.flatten('F')
    #zpos = np.zeros_like(xpos)
    
    
    xpos, ypos = np.meshgrid(xedges[:-1], yedges[:-1], indexing="ij")
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    zpos = 0
    
    
    dx = 0.1 * np.ones_like(zpos)
    dy = 0.1 * np.ones_like(zpos)
    dz = g.ravel()
    
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color= None, zsort='average')
    
    plt.show()









