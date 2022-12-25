#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:06:07 2020

@author: JaylenJames
"""

import matplotlib.pyplot as plt
import numpy as np
from IPython import get_ipython
from stereogram_function import stereogram


def plot_discrete_3d_likelihood(frequencies, 
                                color = None, title = None, fontsize = 12, xlabel = '$b$', ylabel = '$x^2$',
                                xedges = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
                                yedges = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
                                Light_Source = None):
    """
    Return a plot of the 3D distribution of freuquencies, given frequency values at each location.
        The plot will be a 10 by 10 plot from 0 to 1 in both randomm variable
        spaces if edges are not provided.
        
        The xedges and yedges must be equal in symmetry in this version.
    
    Parameters
    ----------
    section_measurement_frequencies : 2 by 2 matrix
        The frequencies associated with m and y^2 measurments from a section.
    dtype : float
    default: Data from CO's paper
    
    
        
    Returns
    -------
    out : ndarray
        Frequencies of 2D measurements with EM applied
        
    
    Examples
    --------
   
    """
    
    frequencies = frequencies
    
    
    
    #xedges = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    #yedges = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.xaxis.set_rotate_label(False) 
    ax.yaxis.set_rotate_label(False)    # disable auto rotation
    
    ax.set_xlabel(xlabel, fontsize = fontsize, rotation = 0)
    ax.set_ylabel(ylabel, fontsize = fontsize, rotation = 0)
    
    
    ax.view_init(28, 55)
    plt.title(title, fontsize = fontsize)
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
    
    
    dx = xedges[1] * np.ones_like(zpos)
    dy = yedges[1] * np.ones_like(zpos)
    dz = frequencies.ravel()
    
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color= color, zsort='average', lightsource = Light_Source)
    
    
    
    
    
    #plt.show()
    return fig, ax