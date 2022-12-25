# -*- coding: utf-8 -*-

"""

Created on Mon Oct 18 19:36:00 2021

@author: Jaylen

Script uses memmaping efficiently, in effort to be run on VM 1.

"""

 

import numpy as np

import pandas as pd
import scipy.integrate as integrate
import time
#import tables as tb
#import h5py

from tempfile import TemporaryFile
import os
import csv
import sys

from plot_discrete_3d_likelihood_function import plot_discrete_3d_likelihood

from plot_COs_inputs_and_unfolded_distribution import *
#from gen_rand_two_dim_sample_func import gen_rand_two_dim_sample_func
from gen_rand_two_dim_Chan_sample_func import gen_rand_two_dim_Chan_sample_func


t_start_script = time.time()


#os.chdir('/home/jupyter/MTR_Project_Folder')


number_of_samples = 200 
multiple = 20 #109
number_of_bins_selection = 10


###############################################vvvvvvvvvvvvvvvvvvvvvvvvvvv
#Experimental Data, Filtered Data - Txt read
filt_min_semi_axis = np.loadtxt(r"filtered_minor_semi_axis_values.txt")
filt_maj_semi_axis = np.loadtxt(r"filtered_major_semi_axis_values.txt")
filt_y_square = np.loadtxt(r"filtered_y_square_values.txt")

#Down select from experimental data using a uniform random distribution
filt_min_semi_axis = np.random.choice(filt_min_semi_axis, number_of_samples)
filt_maj_semi_axis = np.random.choice(filt_maj_semi_axis, number_of_samples)
filt_y_square = np.random.choice(filt_y_square, number_of_samples)

filt_min_semi_axis_norm = filt_min_semi_axis/np.max(filt_min_semi_axis) # would multiply by largest minor semi axis measurement length in data set. For CO's synthetic data, there is no such value. 
filt_maj_semi_axis = filt_min_semi_axis_norm/(np.sqrt(1 - filt_y_square))

# np.savez("./outfile_two_dim_sample_vals_selected", filt_min_semi_axis_norm = filt_min_semi_axis_norm, 
#                                                      filt_maj_semi_axis = filt_maj_semi_axis,
#                                                     filt_y_square = filt_y_square)



a_hist, x_edges_, y_edges_ = np.histogram2d(filt_min_semi_axis_norm, filt_y_square, bins=number_of_bins_selection, range=[[0, 1], [0, 1]], density = True)
a_hist_norm = a_hist/np.sum(a_hist)
#plot_discrete_3d_likelihood(a_hist_norm, xlabel = '$m$', ylabel = '$y^2$')


min_semi_axis_norm = pd.DataFrame(filt_min_semi_axis_norm)
maj_semi_axis = pd.DataFrame(filt_maj_semi_axis)
ellipse_eccentricity = pd.DataFrame(filt_y_square)

###############################################^^^^^^^^^^^^^^^^^^^^^^^^^^^^



###############################################vvvvvvvvvvvvvvv
#Sythetic Data
# CO_two_dim_freq = double_int_results
# CO_unfolded_freq = g
# CO_unfolded_freq_normalized = CO_unfolded_freq/np.sum(CO_unfolded_freq)


# min_semi_axis_norm, filt_y_square = gen_rand_two_dim_sample_func(CO_two_dim_freq, number_of_samples, bins_per_var = number_of_bins_selection)


# np.savez("./outfile_synth_two_dim_sample_vals_selected", min_semi_axis_norm = min_semi_axis_norm, 
#                                                    filt_y_square = filt_y_square)



# a_hist, x_edges_, y_edges_ = np.histogram2d(min_semi_axis_norm.iloc[:][0], filt_y_square.iloc[:][0], bins=number_of_bins_selection, range=[[0, 1], [0, 1]], density = True)
# a_hist_norm = a_hist/np.sum(a_hist)
# #plot_discrete_3d_likelihood(a_hist_norm, xlabel = '$m$', ylabel = '$y^2$')



# filt_min_semi_axis = min_semi_axis_norm
# #filt_maj_semi_axis = filt_maj_semi_axis
# filt_y_square = filt_y_square


# filt_min_semi_axis_norm = filt_min_semi_axis #filt_min_semi_axis/np.max(filt_min_semi_axis) # would multiply by largest minor semi axis measurement length in data set. For CO's synthetic data, there is no such value. 
# filt_maj_semi_axis = filt_min_semi_axis_norm/(np.sqrt(1 - filt_y_square))



# # OPTIONAL - Addition of Noise to sample
# # rand_vals = np.random.uniform(-0.01, 0.01, number_of_samples)

# # filt_min_semi_axis_norm = filt_min_semi_axis_norm + rand_vals[:,None] #Addition of noise to the sample
# # filt_maj_semi_axis = filt_maj_semi_axis + rand_vals[:,None] #Addition of noise to the sample
# # filt_y_square = filt_y_square + rand_vals[:,None] #Addition of noise to the sample

# np.random.seed(10)

# filt_min_semi_axis_norm = np.random.normal(0 + filt_min_semi_axis_norm, 0.01) #Addition of noise to the sample
# filt_maj_semi_axis = np.random.normal(0 + filt_maj_semi_axis, 0.01) #Addition of noise to the sample
# filt_y_square = np.random.normal(0 + filt_y_square, 0.01) #Addition of noise to the sample


# min_semi_axis_norm = pd.DataFrame(filt_min_semi_axis_norm)
# maj_semi_axis = pd.DataFrame(filt_maj_semi_axis)
# ellipse_eccentricity = pd.DataFrame(filt_y_square)


# np.savez("./outfile_noise_for_synth_two_dim_sample_vals_selected", min_semi_axis_norm = min_semi_axis_norm, 
#                                                                maj_semi_axis = maj_semi_axis,
#                                                                ellipse_eccentricity = ellipse_eccentricity, 
#                                                                filt_min_semi_axis_norm = filt_min_semi_axis_norm,
#                                                                filt_maj_semi_axis = filt_maj_semi_axis,
#                                                                filt_y_square = filt_y_square)
#################################################^^^^^^^^^^^^^^^^^^^^^^^








###########################################################################################################################################
# End of selecting input data
###########################################################################################################################################




#min_semi_axis_norm = pd.DataFrame(min_semi_axis/np.max(min_semi_axis))
del maj_semi_axis # min_semi_axis

ellipse_dataset_size = number_of_samples #len(maj_semi_axis)


# define Phi function
def phi_func(min_semi_axis_len, eccentricity_2D):

    int_function = lambda theta: ((1.0 - eccentricity_2D * (np.sin(theta))**2)**0.5 * np.sin(theta))
    integral_component, _ = integrate.quad(int_function, 0.0, np.pi/2.0 )

    phi = min_semi_axis_len * (1.0 - eccentricity_2D)**-0.5 * integral_component

    return phi


t_ini = time.time()

# Initialize probabilities prior to first iteration of EM
p_new_init = np.squeeze(np.ones([ellipse_dataset_size, 1], dtype=np.float32) / ellipse_dataset_size)


p_new = np.memmap("./p_new_file.dat", dtype = "float32", mode = "w+", shape=(ellipse_dataset_size,))
p_new[:] = p_new_init[:]

del p_new_init

t_fin_calculating_theta_components = time.time()

print("Done Initializing p_new. Done with first cell. Time for completion: ", t_fin_calculating_theta_components - t_start_script)






t_init_1_3 = time.time()
i_p_array_common_init = np.tile(np.arange(0,ellipse_dataset_size), [int(ellipse_dataset_size/multiple), 1]).squeeze()

i_p_array_common = np.memmap("./i_p_array_common_file.dat", dtype = "int32", mode = "w+", shape=(int(ellipse_dataset_size/multiple), ellipse_dataset_size))
i_p_array_common[:] = i_p_array_common_init[:]

del i_p_array_common_init
#j_array_1_3 = np.tile(np.atleast_2d(np.arange(0,int(ellipse_dataset_size/17))).T, [1, ellipse_dataset_size]).squeeze()
#
##first half of ind val test
#indicator_value_test = np.squeeze((np.int32(min_semi_axis_norm.values[i_p_array_common] > min_semi_axis_norm.values[j_array_1_3]) *\
#                                   np.int32(ellipse_eccentricity.values[i_p_array_common] > ellipse_eccentricity.values[j_array_1_3]) *\
#                                   np.int32(ellipse_eccentricity.values[i_p_array_common] != 1))) > 0
#
#
#theta_components = np.zeros(np.shape(indicator_value_test), dtype=np.float32)
#theta_components[indicator_value_test] = (min_semi_axis_norm.values[i_p_array_common[indicator_value_test], 0]**-1.0 *\
#                                           ((min_semi_axis_norm.values[i_p_array_common[indicator_value_test], 0]**2.0 - min_semi_axis_norm.values[j_array_1_3[indicator_value_test], 0]**2.0) *\
#                                           (1.0 - ellipse_eccentricity.values[i_p_array_common[indicator_value_test], 0]) * \
#                                           (ellipse_eccentricity.values[i_p_array_common[indicator_value_test],0] - ellipse_eccentricity.values[j_array_1_3[indicator_value_test], 0]) )**-0.5)




#del j_array_1_3, indicator_value_test

t_fin_1_3 = time.time()
print("Time for cell completion:", t_fin_1_3 - t_init_1_3)

theta_components_2 = np.memmap("./theta_components_file.dat", dtype = "float32", mode = "w+", shape=(ellipse_dataset_size, ellipse_dataset_size))

start_row = 0
j_start_row = 1
j_array_2_3 = np.memmap("./j_array_2_3_file.dat", dtype = "int32", mode = "w+", shape=(int((ellipse_dataset_size/multiple)), ellipse_dataset_size))
for rep in np.arange(1,multiple):
    
    j_end_row = int(ellipse_dataset_size/(multiple/(rep)))
    j_array_2_3_init = np.tile(np.atleast_2d(np.arange(j_start_row, j_end_row)).T, [1, ellipse_dataset_size]).squeeze()
    j_start_row = int(ellipse_dataset_size/(multiple/rep))
    
    j_array_2_3[0:np.size(j_array_2_3_init,axis=0),:] = j_array_2_3_init[:]
    
    del j_array_2_3_init
    
    #first half of ind val test
    indicator_value_test = np.squeeze((np.int32(min_semi_axis_norm.values[i_p_array_common] > min_semi_axis_norm.values[j_array_2_3]) *\
                                           np.int32(ellipse_eccentricity.values[i_p_array_common] > ellipse_eccentricity.values[j_array_2_3]) *\
                                           np.int32(ellipse_eccentricity.values[i_p_array_common] != 1))) > 0

    theta_components_part = np.zeros(np.shape(indicator_value_test), dtype=np.float32)
    theta_components_part[indicator_value_test] = (min_semi_axis_norm.values[i_p_array_common[indicator_value_test], 0]**-1.0 *\
                                               ((min_semi_axis_norm.values[i_p_array_common[indicator_value_test], 0]**2.0 - min_semi_axis_norm.values[j_array_2_3[indicator_value_test], 0]**2.0) *\
                                               (1.0 - ellipse_eccentricity.values[i_p_array_common[indicator_value_test], 0]) * \
                                               (ellipse_eccentricity.values[i_p_array_common[indicator_value_test],0] - ellipse_eccentricity.values[j_array_2_3[indicator_value_test], 0]) )**-0.5)
    del indicator_value_test
    #theta_components = np.append(theta_components, theta_components_part, axis = 0)
    #indicator_value_test = np.append(indicator_value_test, indicator_value_test_2_3, axis = 0)
    
    
    end_row = np.int(ellipse_dataset_size/(multiple/rep))
    #end_row = np.size(theta_components_part, axis = 0)
    theta_components_2[start_row:end_row,:] = theta_components_part[:]
    
    start_row = np.int(ellipse_dataset_size/(multiple/rep))
    #start_row = int(end_row +1)
    #del j_array_2_3
    
    
    print("A rep is complete", rep)
    
    del theta_components_part#
                          
del i_p_array_common, j_array_2_3 #, indicator_value_test




t_init_mem = time.time()


# Save theta components to disk, remove this variable, and memmap it
#t_c_map = np.memmap("./theta_components_file.dat", dtype = "float32", mode = "w+", shape=(ellipse_dataset_size, ellipse_dataset_size))

#t_c_map.flush()
#theta_components_2 = np.memmap("./theta_components_file.dat", dtype='float32', mode='r', shape=(ellipse_dataset_size, ellipse_dataset_size))
#del theta_components

t_fin_mem = time.time()

print("Done memmapping theta components. Time for completion:", t_fin_mem - t_init_mem)





#t_fin_mem = time.time()



phi = np.squeeze(np.zeros([1, ellipse_dataset_size], dtype=np.float32))
for i in np.arange(0, ellipse_dataset_size): 
    phi[i] = phi_func(min_semi_axis_norm.values[i], ellipse_eccentricity.values[i])
phi_m1 = 1.0 / phi
del phi, ellipse_eccentricity, min_semi_axis_norm
phi_m1[phi_m1 == np.inf] = np.spacing(np.float32(1.0))

#numerator_components = np.zeros(np.shape(theta_components_2), dtype=np.float32)

print("Done Calculating Phi. Starting EM Loop.")

while_loop_counter = 0

t_init_while_loop = time.time()

while True:

    # E-step
    numerator_components = p_new * theta_components_2
    
    
    #Apply Memmap for numerator_components
    #np.save("../numerator_components_file.dat", numerator_components)
    numerator_components_lite = np.memmap("./numerator_components_file.dat", dtype = "float32", mode = "w+", shape=(ellipse_dataset_size, ellipse_dataset_size))
    numerator_components_lite[:] = numerator_components[:]
    #n_c_map.flush()
    #numerator_components_lite = np.memmap("./numerator_components_file.dat", dtype='float32', mode='r', shape=(ellipse_dataset_size, ellipse_dataset_size))
    del numerator_components
    
    
    denominator = np.sum(numerator_components_lite, axis=1, dtype=np.float32)
    
    
    division_for_eta = numerator_components_lite/denominator[:,None]
    division_for_eta[np.isnan(division_for_eta)] = np.spacing(np.float32(1.0))
    
    #Apply Memmap for division_for_eta
    #np.save("../division_for_eta_file.dat", division_for_eta)
    division_for_eta_lite = np.memmap("./division_for_eta_file.dat", dtype = "float32", mode = "w+", shape=(ellipse_dataset_size, ellipse_dataset_size))
    division_for_eta_lite[:] = division_for_eta[:]
    #d_e_map.flush()
    #division_for_eta_lite = np.memmap("./division_for_eta_file.dat", dtype='float32', mode='r', shape=(ellipse_dataset_size, ellipse_dataset_size))
    del division_for_eta
    
    eta_i = np.sum(division_for_eta_lite, axis=0)
    eta_i[eta_i == np.inf] = np.spacing(np.float32(1.0))
    
    # M-step
    p_old = p_new
    tmp = eta_i * phi_m1
    
    
    p_new = tmp / np.sum(tmp)
    
    convergence_metric = np.sum((p_new - p_old)**2)
    
    while_loop_counter += 1
    print("While loops completed:", while_loop_counter)
    if convergence_metric < 0.05:
        break

t_fin = time.time()

print("Done with EM")

#%%
  
"""

////////////////////////////////////////////////////

"""




#Reassign the variables since they were deleted
min_semi_axis = pd.DataFrame(filt_min_semi_axis)
maj_semi_axis = pd.DataFrame(filt_maj_semi_axis)
ellipse_eccentricity = pd.DataFrame(filt_y_square)
min_semi_axis_norm = pd.DataFrame(min_semi_axis/np.max(min_semi_axis))

EM_b, EM_x2 = gen_rand_two_dim_Chan_sample_func(p_new, min_semi_axis_norm, filt_y_square, size_of_drawing = number_of_samples)

a_Chan_hist, x_edges_Ch, y_edges_Ch = np.histogram2d(EM_b, EM_x2, bins=number_of_bins_selection, range=[[0, 1], [0, 1]], density = True)
a_Chan_hist_norm = a_Chan_hist/np.sum(a_Chan_hist)
#plot_discrete_3d_likelihood(a_Chan_hist_norm)

t_end_script = time.time()

print("Time for full Script:", t_end_script - t_start_script)
print("Time for full EM Process:", t_fin - t_ini)
print("Time for EM While Loop:", t_fin - t_init_while_loop)


#print("New mem map method time:", t_fin_mem - t_init_mem)

full_script_time = t_end_script - t_start_script
full_EM_process_time = t_fin - t_ini
full_EM_while_loop_time = t_fin - t_init_while_loop

#Saving variablez
# np.savez("./outfile_important_data", p_new = p_new, 
#                       a_hist_norm = a_hist_norm,
#                       a_Chan_hist_norm = a_Chan_hist_norm,
#                       #CO_unfolded_freq = CO_unfolded_freq,
#                       full_script_time = full_script_time,
#                       full_EM_process_time = full_EM_process_time,
#                       full_EM_while_loop_time = full_EM_while_loop_time)

#loading variablez
# important_output_arrays = np.load('./outfile_important_data_100K_Exp_Run_3.npz')
# sorted(important_output_arrays.files)
# p_new = important_output_arrays['p_new'] #etc.
# a_hist_norm = important_output_arrays['a_hist_norm']
# a_Chan_hist_norm = important_output_arrays['a_Chan_hist_norm']

# plot_discrete_3d_likelihood(a_hist_norm, xlabel = '$m$', ylabel = '$y^2$')
# plot_discrete_3d_likelihood(a_Chan_hist_norm, color = 'olivedrab')




#%%
#Only use this section when calculating probabilities for the synthetic data.

# #Plotting heatmap
# import matplotlib
# import mpl_toolkits
# import seaborn as sns



# full_subsampled_diff = a_Chan_hist_norm - CO_unfolded_freq
# fig = plt.figure()

# ax = sns.heatmap(full_subsampled_diff, vmin = -0.03, vmax = 0.03)

# ax.set_title(number_of_samples, fontsize=18)
# ax.set_xlabel('$x^2$', fontsize=18)
# ax.set_ylabel('$b$', fontsize=18)

# results_RSS = np.sum((full_subsampled_diff)**2)
# results_RS = np.sum(np.abs(full_subsampled_diff))
# print("Residual Sum of Squares:", results_RSS)
# print("Residual Sum (More Accurate):", results_RS)


