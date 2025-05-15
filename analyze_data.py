from profile_class import profile_group
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import ticker
import matplotlib as mpl
import pickle


def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)




path_1 = "remote_results/t0_005/rot/profiles"
path_2 = "remote_results/t0_005/no_rot/profiles"
path_3 = "remote_results/t0_005/constant_flux/profiles"
path_4 = "remote_results/t0_005/constant_rot/profiles"
path_5 = "remote_results/t0_005/F_ratio_5/profiles"


num_files_1 = 80
num_files_2 = 80
num_files_3 = 43
num_files_4 = 55
num_files_5 = 70

test_1_switch = True
test_2_switch = True
test_3_switch = True
test_4_switch = True
test_5_switch = True


if test_1_switch == True:
    test_1 = profile_group()
    test_1.define_constant(fratio = 10, t_0 = 0.005, omega_num =  0.1 / 3e-6, pr= 0.1, ra_number = 2e9, lamb = 1e-7, tau = 0.1)
    test_1.load_profiles(path_1, num_files_1, box_dimension =3)
    test_1.measure_velocity(path_1, num_files_1, box_dimension =3)
    test_1.average_profile_time_based(num_to_combine = 100, average_velocity = True)
    # test_1.average_profile_grid_based(display_resolution = 256)



if test_2_switch == True:
    test_2 = profile_group()
    test_2.define_constant(fratio = 10, t_0 = 0.005, pr= 0.1, ra_number = 2e9, lamb = 1e-7, tau = 0.1)
    test_2.load_profiles(path_2, num_files_2, box_dimension =3)
    test_2.measure_velocity(path_2, num_files_2, box_dimension =3)
    test_2.average_profile_time_based(num_to_combine = 100, average_velocity = True)
    # test_2.average_profile_grid_based(display_resolution = 256)


if test_3_switch == True:
    test_3 = profile_group()
    test_3.define_constant(fratio = 10, pr= 0.1, ra_number = 2e9, lamb = 1e-7, tau = 0.1)
    test_3.load_profiles(path_3, num_files_3, box_dimension =3)
    test_3.measure_velocity(path_3, num_files_3, box_dimension =3)
    test_3.average_profile_time_based(num_to_combine = 100, average_velocity = True)
    # test_3.average_profile_grid_based(display_resolution = 256)



if test_4_switch == True:
    test_4 = profile_group()
    test_4.define_constant(fratio = 10, omega_num =  0.1 / 3e-6, pr= 0.1, ra_number = 2e9, tau = 0.1)
    test_4.load_profiles(path_4, num_files_4, box_dimension =3)
    test_4.measure_velocity(path_4, num_files_4, box_dimension =3)
    test_4.average_profile_time_based(num_to_combine = 100, average_velocity = True)
    # test_4.average_profile_grid_based(display_resolution = 256)


if test_5_switch == True:
    test_5 = profile_group()
    test_5.define_constant(fratio = 5, omega_num =  0.1 / 3e-6, pr= 0.1, ra_number = 2e9, tau = 0.1)
    test_5.load_profiles(path_5, num_files_5, box_dimension =3)
    test_5.measure_velocity(path_5, num_files_5, box_dimension =3)
    test_5.average_profile_time_based(num_to_combine = 100, average_velocity = True)
    # test_5.average_profile_grid_based(display_resolution = 256)





# method: "heating" or "cooling"
# boundary: "relative" or "absolute"
# relative/absolute_percentra = float number 
if test_1_switch == True:
    test_1.analyze_c(method = "cooling", relative_percenta = 0.05)
if test_2_switch == True:
    test_2.analyze_c(method = "cooling", relative_percenta = 0.05)
if test_3_switch == True:
    test_3.analyze_c(method = "cooling", relative_percenta = 0.03)
if test_4_switch == True:
    test_4.analyze_c(method = "cooling", relative_percenta = 0.05)
if test_5_switch == True:
    test_5.analyze_c(method = "cooling", relative_percenta = 0.05)



if test_1_switch == True:
    test_1.analyze_cflux(method = "cooling")
if test_2_switch == True:
    test_2.analyze_cflux(method = "cooling")
if test_3_switch == True:
    test_3.analyze_cflux(method = "cooling")
if test_4_switch == True:
    test_4.analyze_cflux(method = "cooling")
if test_5_switch == True:
    test_5.analyze_cflux(method = "cooling")



if test_1_switch == True:
    test_1.analyze_rho(method = "cooling", relative_percenta = 0.01)
if test_2_switch == True:
    test_2.analyze_rho(method = "cooling", relative_percenta = 0.01)
if test_3_switch == True:
    test_3.analyze_rho(method = "cooling", relative_percenta = 0.01)
if test_4_switch == True:
    test_4.analyze_rho(method = "cooling", relative_percenta = 0.01)
if test_5_switch == True:
    test_5.analyze_rho(method = "cooling", relative_percenta = 0.01)




# function_to_fit: 1. t^0.5 2. power law 
# outcome domain: array, specify the range where the function is evaluated
# c, c_flux, rho, dsdz: False/True
# Stacking results: allows not to erase the previous results
t_i = 0.005
time = np.logspace(-6, 1, 1000)
time_short = np.logspace(np.log10(t_i), 1, 1000)



if test_1_switch == True:
    fit_2 = test_1.fit_for_h(function_to_fit = 54, outcome_domain = time_short, c_flux = True, t_i  = t_i)

if test_2_switch == True:
    fit_3 = test_2.fit_for_h(function_to_fit = 53, outcome_domain = time_short, c_flux = True, t_i  = t_i)

if test_3_switch == True:
    fit_1 = test_3.fit_for_h(function_to_fit = 51, outcome_domain = time_short, c_flux = True, t_i  = t_i)


if test_5_switch == True:
    fit_5 = test_5.fit_for_h(function_to_fit = 2, outcome_domain = time_short, c = True, t_i  = t_i)





velo_range =  "all"
avg_met = "avg"
test_1.analyze_velocity(test_1.h_rho, box_dimension = 3, velocity_range = velo_range, average_method = avg_met)
test_2.analyze_velocity(test_2.h_rho, box_dimension = 3, velocity_range = velo_range, average_method = avg_met)
test_3.analyze_velocity(test_3.h_rho, box_dimension = 3, velocity_range = velo_range, average_method = avg_met)
test_4.analyze_velocity(test_4.h_rho, box_dimension = 3, velocity_range = velo_range, average_method = avg_met)
test_5.analyze_velocity(test_5.h_rho, box_dimension = 3, velocity_range = velo_range, average_method = avg_met)



test_1.average_all_h()
test_2.average_all_h()
test_3.average_all_h()
test_4.average_all_h()
test_5.average_all_h()



test_1.compute_epsilon(h = test_1.h_rho)
test_2.compute_epsilon(h = test_2.h_rho)
test_3.compute_epsilon(h = test_3.h_rho)
test_4.compute_epsilon(h = test_4.h_rho)
test_5.compute_epsilon(h = test_5.h_rho)



test_1.compute_gamma(h = test_1.h_rho)
test_2.compute_gamma(h = test_2.h_rho)
test_3.compute_gamma(h = test_3.h_rho)
test_4.compute_gamma(h = test_4.h_rho)
test_5.compute_gamma(h = test_5.h_rho)


test_package = [test_1, test_2, test_3, test_4, test_5]

save_object(test_package, 'test_analysis.pkl')


