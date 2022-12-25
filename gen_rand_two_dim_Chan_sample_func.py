# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 09:32:01 2021

@author: Jaylen

This function is to be used with the output of the 3D frequencies predicted
    using Chan's EM method.
    

"""

import numpy as np
from scipy import stats

#Begin function here. Return the b and x2 values

def gen_rand_two_dim_Chan_sample_func(probabilities, m_min_semi_axis_values, y2_eccentricities, size_of_drawing = 100):

    identifier = np.arange(len(probabilities))
    
    probabilities_Chan_EM = probabilities 
    min_semi_axis_norm = m_min_semi_axis_values[:len(probabilities)]
    y_square = y2_eccentricities[:len(probabilities)]
    
    
    probs_and_measures = np.column_stack((identifier, probabilities_Chan_EM, min_semi_axis_norm, y_square))
    
    
    xk = probs_and_measures[:,0]
    pk = probs_and_measures[:,1]
    
    
    EM_probabilities_class_Chan = stats.rv_discrete(name='EM_probabilites_class_Chan', values=(xk, pk))
    
    size_of_drawing = size_of_drawing
    sampled_identifier_values = EM_probabilities_class_Chan.rvs(size = size_of_drawing)
    
    
    spheroid_sizes_via_Chan_EM = probs_and_measures[sampled_identifier_values]
    
    EM_b = spheroid_sizes_via_Chan_EM[:,2]
    EM_x2 = spheroid_sizes_via_Chan_EM[:, 3]
    
    return EM_b, EM_x2




if __name__ == "__main__":
    
    # Use this to run the above function for testing and debugging.
    
    ##################################
    # Import EM via Chan Probabilities
    ##################################
    
    #Mac
    try:
        probabilities_Chan_EM = np.genfromtxt('/Users/JaylenJames/Stereology Project/EM_probabilities_Chan.csv', delimiter=',')
    
    #PC    
    except:
        probabilities_Chan_EM = np.genfromtxt('C:/Users/Jaylen/Stereology Work/Dist_of_Micro/EM_probabilities_Chan.csv', delimiter=',')
        
    
    # Full Data - Microstructure measurement data into Numpy arrays
    maj_axis = np.loadtxt("180806_MTR_Metrics_allMajAx.txt")
    min_axis = np.loadtxt("180806_MTR_Metrics_allMinAx.txt")
    
    min_semi_axis = min_axis/2
    maj_semi_axis = maj_axis/2
    y_square = 1 - (min_axis/maj_axis)**2
    
    min_semi_axis_norm = min_semi_axis/np.max(min_semi_axis)
    
    
    probabilities_dataset_size = len(probabilities_Chan_EM)
    
    #Taking norm of minor semi axis as if there were 'probabilities_dataset_size' number of measurements
    min_semi_axis_norm = min_semi_axis[:probabilities_dataset_size]/np.max(min_semi_axis[:probabilities_dataset_size])
    
    EM_b, EM_x2 = gen_rand_two_dim_Chan_sample_func(probabilities_Chan_EM, min_semi_axis_norm, y_square)


