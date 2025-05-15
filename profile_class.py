import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import h5py
import scipy 
import os



# fitting function 
# t,h are the data, time domiain sepcify the desired t range where the fitting returns h(t)
def t_05(t, h, time_domain, t_i, h_i):

    func_05 = lambda x, k: (h_i**2 + k * (x- t_i) ) **0.5
    fit, _ = scipy.optimize.curve_fit(func_05, t, h)

    return fit, func_05(time_domain, fit[0])


def t_06(t, h, time_domain, t_i, h_i):

    func_05 = lambda x, k: h_i + k * ( x **(0.6) - t_i **0.6)
    fit, _ = scipy.optimize.curve_fit(func_05, t, h)

    return fit, func_05(time_domain, fit[0])


def t_power(t, h, time_domain, t_i, h_i):
    func_power = lambda x, k, b: (h_i**(1/b) + k * (x- t_i) ) **b
    fit, _ = scipy.optimize.curve_fit(func_power, t, h)

    return fit, func_power(time_domain, fit[0], fit[1])


def h_nr_fitk(t, h, time_domain, t_0, f_ratio, t_i, h_i):

    # k =  2* (1 - epsilon + 2 * gamma) * |ds_0/dz|
    func_hnr = lambda t, k: np.sqrt(h_i**2  + k * f_ratio * t_0 * np.log (abs( (t + t_0) / t_i / 2)) )


    fit, _ = scipy.optimize.curve_fit(func_hnr, t, h)

    return fit, func_hnr(time_domain, fit[0])




def h_wr_fitk(t, h, time_domain, t_0, f_ratio, ra_num, pr_num, omega_num, t_i, h_i):

    # k =  24 * gamma /gamma'
    func_hwr = lambda t, k : np.power( h_i**(12/5) + k * f_ratio**(6/5) * (ra_num*pr_num/omega_num**3)**(1/5)
     * t_0 * ( (t_0 / (t_0 + t_i))**0.2 - (t_0 / (t_0 + t))**0.2 ), 5/12)

    fit, _ = scipy.optimize.curve_fit(func_hwr,t , h)
    print(t_i, h_i)
    return fit, func_hwr(time_domain,fit[0] )




def numer_fit_func(t, h, time_domain, t_0, f_ratio, ra_num, pr_num, omega_num, t_i, h_i, fit = "constant_rot", method = "midpoint"):
    

    if fit == "constant_rot":
        def dhEdt(t, hE, ep, gamma):
            dEdt = 1
            dhdt =  gamma * f_ratio **(6/5) * (pr_num*ra_num/ omega_num**3)**(1/5) * hE[0]**(3/5) / ( hE[0]**2 /2 - f_ratio * (1- ep) * hE[1] )
            return np.array([dhdt, dEdt])
    if fit == "constant_nr":
        def dhEdt(t, hE, ep, gamma):
            dEdt = 1
            dhdt = gamma *  f_ratio  * hE[0] / ( hE[0]**2 /2 - f_ratio * (1- ep) * hE[1] )
            return np.array([dhdt, dEdt])
    if fit == "deflux_rot": 
        def dhEdt(t, hE, ep, gamma):
            dEdt = t_0 / (t_0 + t)
            g_th  = hE[0]**2 /2 - f_ratio * (1- ep) * (t_0 * hE[1] )
            fh =  gamma * f_ratio **(6/5) * (pr_num*ra_num/ omega_num**3)**(1/5) * hE[0]**(3/5) * (t_0 / (t_0 + t))**(6/5)
            dhdt = fh / g_th
            return np.array([dhdt, dEdt])
    if fit == "deflux_nr":
        def dhEdt(t, hE, ep, gamma):
            dEdt = t_0 / (t_0 + t)
            g_th  = hE[0]**2 /2 - f_ratio * (1- ep) * (t_0 * hE[1] )
            fh =   gamma *  f_ratio * (t_0/ (t_0 + t)) * hE[0]
            dhdt = fh / g_th
            return np.array([dhdt, dEdt])



    # evaluation through different methods
    if method == "rk4":
        def func_eval(t, ep, gamma):

            dt = np.diff(t)
            h = np.zeros(len(t))
            h[0] = h_i

            for i in range(len(t)-1):    
                f = dhdt(t[i], h[i], ep, gamma)
                f1 = dhdt(t[i] + dt/2, h[i] + f*dt[i]/2, ep, gamma)
                f2 = dhdt(t[i] + dt/2, h[i] + f1*dt[i]/2, ep, gamma)
                f3 = dhdt(t[i], h[i] + f2*dt[i], ep, gamma)
                x[i] = h[i] + dt[i]*(f + 2*f1 + 2*f2 + f3)/6

            return h


    if method == "midpoint":
        def func_eval(t, ep, gamma):

            dt = np.diff(t)

            for i in range(len(t)-1):
                h1 = hE[i] + 0.5*dt[i] * dhEdt(t[i], hE[i], ep, gamma)
                f1 = dhEdt(t[i]+ dt[i]/2, h1, ep, gamma)
                hE[i+1] = hE[i] + f1*dt[i]

            return hE[:,0]

    # initial values for fit
    hE = np.zeros((len(t), 2))
    hE[0, 0] = h_i
    print(fit)
    if fit == "constant_rot" or "constant_nr":
        hE[0, 1] = 0.005

    if fit == "deflux_rot" or "deflux_nr":
        hE[0, 1] = 0.005 * np.log((0.005 + t_i) / 0.005)


    # fit
    fit, _ = scipy.optimize.curve_fit(func_eval, t, h)


    # initial values for evaluation
    hE = np.zeros((len(time_domain), 2))
    hE[0, 0] = h_i
    if fit.all() == "constant_rot" or "constant_nr":
        hE[0, 1] = 0.005

    if fit.all() == "deflux_rot" or "deflux_nr":
        hE[0, 1] = 0.005 * np.log((0.005 + t_i) / 0.005)


    return fit, func_eval(time_domain, 0.4, 0.9)

    # return fit, func_eval(time_domain, fit[0], fit[1])








def find_highest_peak(array):
    n = len(array)
    highest_peak_index = 0
    highest_peak_value = float('-inf')
    
    for i in range(n):
        # Check if it's a peak
        if (i == 0 or array[i] >= array[i-1]) and (i == n-1 or array[i] >= array[i+1]):
            # Update highest peak if current is larger
            if array[i] > highest_peak_value:
                highest_peak_value = array[i]
                highest_peak_index = i
    
    return highest_peak_index, highest_peak_value





class profile_group:
    def _init_(self):

        # when initilize the objects the constant the trial 
        # should be input manually for the setup, which contains
        # all the parameters: fratio, ra, tau, prandtl etc.
        self.num_of_files = 0


    def define_constant(self, fratio = 0, ra_number = 0, pr = 0, tau = 0, t_0 = 0, omega_num = 0, lamb = 0):
        self.F_ratio = fratio
        self.Ra_value = ra_number
        self.prandtl = pr
        self.tau = tau
        self.t_0 =  t_0
        self.omega_num = omega_num
        self.lamb = lamb 


    def exmaine_data(self, path):
        file_path = path + str("/profiles_s1.h5")
        f = h5py.File(file_path, 'r')
        sca =  f["scales"]
        tsk = f['tasks']
  
        self.sca_keys = sca.keys()
        self.tsk_keys = tsk.keys()
        # print("scale keys:", self.sca_keys)
        # print("taks keys:", self.tsk_keys)

        print(tsk["KE_hor"])


    # PART 1: loading and averaging the data

    # Load the profile to initial setup 
    # if one wants to restore to the original resolution ex. 1024
    # set reset it with load_profile again
    def load_profiles(self, path, num_of_profiles, box_dimension = 2, max_write = 100, method = "cooling"):

        # path need to be the string of the folder contains the files
        # num_of_profiles determins the total loaded files, int
        # max_write matches the the number of profile in each pack of the profile
        file_path_pool = path + str("/profiles_s")
        num_file_read = int(num_of_profiles)
        self.dimension = box_dimension


        if self.dimension == 2:
            # for 2D box, there's only z coordinates
            # load in the file once to construct the z coordinate
            path =  f"{file_path_pool}{1}.h5"
            f = h5py.File(path, 'r')
            a_group_key  = list(f.keys())[0]
            s = list(f[a_group_key])[6]
            # store the z coordinates, which is untreated, before averaging
            self.z_coord = f['/scales/'+s][:]

            if self.z_coord.ndim + 1 !=  box_dimension:
                raise ValueError("detected dimension and input dimension unmatched")
                

            # load all the data that will be potentially needed
            self.time = np.zeros([num_file_read * max_write])
            self.c_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.c_flux_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.rho_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.t_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.t_flux_data = np.zeros([num_file_read * max_write, len(self.z_coord)])

            for ind in range(num_file_read):

                path =  f"{file_path_pool}{ind+1}.h5"
                f = h5py.File(path, 'r')
                tsk = f["tasks"]
                sca = f["scales"]

                t_relative = sca["sim_time"][:]
                c_data = np.squeeze( tsk['C'][:], axis = 1)
                c_flux_data = np.squeeze( tsk["C_conv_flux"][:], axis = 1)
                t_data = np.squeeze( tsk["T"][:], axis = 1)

                if method == "heating":
                    rho_data = c_data - t_data
                if method == "cooling":
                    rho_data = c_data - t_data + 1 

                # storing
                self.time[ind * max_write : (ind + 1)* max_write ] = t_relative
                self.c_data[ind * max_write : (ind + 1)* max_write ] = c_data
                self.c_flux_data[ind * max_write : (ind + 1)* max_write] = c_flux_data
                self.rho_data[ind * max_write : (ind + 1)* max_write] = rho_data
                self.t_data[ind * max_write : (ind + 1)* max_write] = t_data


        if self.dimension == 3:

            path =  f"{file_path_pool}{1}.h5"
            f = h5py.File(path, 'r')
            a_group_key  = list(f.keys())[0]
            s = list(f[a_group_key])[6]

            self.z_coord = f['/scales/'+s][:]
            self.time = np.zeros([num_file_read * max_write])
            self.c_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.c_flux_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.rho_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.t_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.t_flux_data = np.zeros([num_file_read * max_write, len(self.z_coord)])


            for ind in range(num_file_read):

                path =  f"{file_path_pool}{ind+1}.h5"
                f = h5py.File(path, 'r')
                tsk = f["tasks"]
                sca = f["scales"]

                t_relative = sca["sim_time"]
                c_data = tsk['C'][:]
                c_flux_data = tsk["C_conv_flux"][:] - self.tau * tsk["C_grad"][:]
                t_data = tsk["T"][:]
                t_flux_data = (tsk["T_conv_flux"][:] - tsk["T_grad"][:]) / self.F_ratio


                # T0 from setup, assuming linear
                T0 = 1
                if method == "heating":
                    rho_data = c_data - t_data
                if method == "cooling":
                    rho_data = c_data - t_data + T0 


                # storing
                self.time[ind * max_write : (ind + 1)* max_write ] = t_relative
                self.c_data[ind * max_write : (ind + 1)* max_write ] = np.squeeze(c_data)
                self.c_flux_data[ind * max_write : (ind + 1)* max_write ] = np.squeeze(c_flux_data)
                self.rho_data[ind * max_write : (ind + 1)* max_write ] = np.squeeze(rho_data)
                self.t_data[ind * max_write : (ind + 1)* max_write ] = np.squeeze(t_data)
                self.t_flux_data[ind * max_write : (ind + 1)* max_write ] = np.squeeze(t_flux_data)


        # count and report the number of file
        print("Profiles start time:", self.time[0])
        print("End time:", self.time[-1])
        print("total Profiles loaded:", len(self.time))




    # grid based averaging to a smaller size, ex. if display_res = 256
    # where original coord 1024, then average over 4 grid point.
    def average_profile_grid_based(self, display_resolution = 256):

        if self.dimension == 2:
            if len(self.z_coord) % display_resolution != 0:
                raise ValueError("Unequal division, choose a number would split the box evenly")
            
            division = len(self.z_coord) // display_resolution
                
            # smooth the grid to match the resolution
            reshaped_zcoord = self.z_coord.reshape(display_resolution, division)
            self.z_coord = np.mean(reshaped_zcoord, axis= -1)

            # take average over some grid points 
            reshaped_cdata = self.c_data.reshape(len(self.time), display_resolution, division)
            reshaped_cfluxdata = self.c_data.reshape(len(self.time), display_resolution, division)
            reshaped_rhodata = self.c_data.reshape(len(self.time), display_resolution, division)

            self.c_data = np.mean(reshaped_cdata, axis= -1)
            self.c_flux_data = np.mean(reshaped_cfluxdata, axis= -1)
            self.rho_data = np.mean(reshaped_rhodata, axis= -1)



        if self.dimension == 3:
            if len(self.z_coord) % display_resolution != 0:
                raise ValueError("Unequal division, choose a number would split the box evenly", f"the box is {len(self.z_coord)}")

            division = len(self.z_coord) // display_resolution
                
            # smooth the grid to match the resolution
            reshaped_zcoord = self.z_coord.reshape(display_resolution, division)
            self.z_coord = np.mean(reshaped_zcoord, axis= -1)

            # take average over some grid points 
            reshaped_cdata = self.c_data.reshape(len(self.time), display_resolution, division)
            reshaped_cfluxdata = self.c_data.reshape(len(self.time), display_resolution, division)
            reshaped_rhodata = self.c_data.reshape(len(self.time), display_resolution, division)

            self.c_data = np.mean(reshaped_cdata, axis= -1)
            self.c_flux_data = np.mean(reshaped_cfluxdata, axis= -1)
            self.rho_data = np.mean(reshaped_rhodata, axis= -1)

        print("... ...")
        print("Coordinates averaging complete")
        print("Now the resolution is:", len(self.z_coord))




    # another way to average is to time_based average 
    def average_profile_time_based(self, num_to_combine = 100, average_velocity = False):

        # One can choose a number of average that not divide the entire pool evenly
        # However, the remainder of the profile will be dicarded,
        # hence the end time maybe different from before averaging the profiles 

        if self.dimension == 2:
            # discard the extra profile at first
            rest_structure = len(self.time)//num_to_combine
            remainder = len(self.time) % num_to_combine

            if remainder != 0:
                rest_time = self.time[: - remainder]
                rest_c = self.c_data[: -remainder]
                rest_cflux = self.c_flux_data[: -remainder]
                rest_rho = self.rho_data[: -remainder]
                rest_temp = self.t_data[: -remainder]
            else:
                rest_time = self.time
                rest_c = self.c_data
                rest_cflux = self.c_flux_data
                rest_rho = self.rho_data
                rest_temp = self.t_data
            
            reshaped_time = rest_time.reshape(rest_structure, num_to_combine)
            reshaped_cdata = rest_c.reshape(rest_structure, num_to_combine, len(self.z_coord))
            reshaped_cfluxdata = rest_cflux.reshape(rest_structure, num_to_combine, len(self.z_coord))
            reshaped_rhodata = rest_rho.reshape(rest_structure, num_to_combine, len(self.z_coord))
            reshaped_tdata = rest_temp.reshape(rest_structure, num_to_combine, len(self.z_coord))

            self.time = np.mean(reshaped_time, axis= -1)
            self.c_data = np.mean(reshaped_cdata, axis = 1)
            self.c_flux_data = np.mean(reshaped_cfluxdata, axis = 1)
            self.rho_data = np.mean(reshaped_rhodata, axis = 1)
            self.t_data = np.mean(reshaped_tdata, axis = 1)

            if average_velocity:
                if remainder != 0:
                    rest_ux = self.u_x_data[: -remainder]
                    rest_uz = self.u_z_data[: -remainder]
                else:
                    rest_ux = self.u_x_data
                    rest_uz = self.u_z_data
                
                reshaped_ux = rest_ux.reshape(rest_structure, num_to_combine, len(self.z_coord))
                reshaped_uz = rest_uz.reshape(rest_structure, num_to_combine, len(self.z_coord))

                self.u_x_data = np.mean(reshaped_ux, axis = 1)
                self.u_z_data = np.mean(reshaped_uz, axis = 1)
            else:
                None

            print("... ...")
            print("Time averaging complete")
            print("Profiles start time:", self.time[0])
            print("End time:", self.time[-1])
            print("total Profiles reduced to:", len(self.time))
            print("dt =", self.time[1] - self.time[0])


        if self.dimension == 3:
            # discard the extra profile at first
            rest_structure = len(self.time)//num_to_combine
            remainder = len(self.time) % num_to_combine

            if remainder != 0:
                rest_time = self.time[: - remainder]
                rest_c = self.c_data[: -remainder]
                rest_cflux = self.c_flux_data[: -remainder]
                rest_rho = self.rho_data[: -remainder]
                rest_temp = self.t_data[: -remainder]
                rest_tempflux = self.t_flux_data[: -remainder]
            else:
                rest_time = self.time
                rest_c = self.c_data
                rest_cflux = self.c_flux_data
                rest_rho = self.rho_data
                rest_temp = self.t_data
                rest_tempflux = self.t_flux_data

            
            reshaped_time = rest_time.reshape(rest_structure, num_to_combine)
            reshaped_cdata = rest_c.reshape(rest_structure, num_to_combine, len(self.z_coord))
            reshaped_cfluxdata = rest_cflux.reshape(rest_structure, num_to_combine, len(self.z_coord))
            reshaped_rhodata = rest_rho.reshape(rest_structure, num_to_combine, len(self.z_coord))
            reshaped_tdata = rest_temp.reshape(rest_structure, num_to_combine, len(self.z_coord))
            reshaped_tempfluxdata = rest_tempflux.reshape(rest_structure, num_to_combine, len(self.z_coord))

            self.time = np.mean(reshaped_time, axis= -1)
            self.c_data = np.mean(reshaped_cdata, axis = 1)
            self.c_flux_data = np.mean(reshaped_cfluxdata, axis = 1)
            self.rho_data = np.mean(reshaped_rhodata, axis = 1)
            self.t_data = np.mean(reshaped_tdata, axis = 1)
            self.t_flux_data = np.mean(reshaped_tempfluxdata, axis = 1)

            if average_velocity:
                if remainder != 0:
                    rest_ux = self.u_x_data[: -remainder]
                    rest_uz = self.u_z_data[: -remainder]
                    rest_uy = self.u_y_data[: -remainder]
                else:
                    rest_ux = self.u_x_data
                    rest_uz = self.u_z_data
                    rest_uy = self.u_y_data
                
                reshaped_ux = rest_ux.reshape(rest_structure, num_to_combine, len(self.z_coord))
                reshaped_uz = rest_uz.reshape(rest_structure, num_to_combine, len(self.z_coord))
                reshaped_uy = rest_uy.reshape(rest_structure, num_to_combine, len(self.z_coord))

                self.u_x_data = np.mean(reshaped_ux, axis = 1)
                self.u_z_data = np.mean(reshaped_uz, axis = 1)
                self.u_y_data = np.mean(reshaped_uy, axis = 1)

            else:
                None

            print("... ...")
            print("Time averaging complete")
            print("Profiles start time:", self.time[0])
            print("End time:", self.time[-1])
            print("total Profiles reduced to:", len(self.time))
            print("dt =", self.time[1] - self.time[0])




    # PART 1 extra: to measure some extra parameters
    # load this before any averaging
    def measure_velocity(self, path,  num_of_profiles, box_dimension = 3, max_write = 100):
        
        file_path_pool = path + str("/profiles_s")
        num_file_read = int(num_of_profiles)
        path =  f"{file_path_pool}{1}.h5"
        f = h5py.File(path, 'r')
        a_group_key  = list(f.keys())[0]
        s = list(f[a_group_key])[6]


        if box_dimension == 2:

            self.u_z_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.u_x_data = np.zeros([num_file_read * max_write, len(self.z_coord)])

            for ind in range(num_file_read):

                path =  f"{file_path_pool}{ind+1}.h5"
                f = h5py.File(path, 'r')
                tsk = f["tasks"]
                sca = f["scales"]

                uz_data = np.squeeze( tsk['KE_vert'][:], axis = 1)
                ux_data = np.squeeze( tsk["KE_hor"][:], axis = 1)

                # storing
                self.u_x_data[ind * max_write : (ind + 1)* max_write] = ux_data
                self.u_z_data[ind * max_write : (ind + 1)* max_write] = uz_data


        if box_dimension == 3:

            self.u_z_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.u_x_data = np.zeros([num_file_read * max_write, len(self.z_coord)])
            self.u_y_data = np.zeros([num_file_read * max_write, len(self.z_coord)])

            for ind in range(num_file_read):

                path =  f"{file_path_pool}{ind+1}.h5"
                f = h5py.File(path, 'r')
                tsk = f["tasks"]
                sca = f["scales"]

                uz_data = np.squeeze( tsk['KE_vert'][:])
                ux_data = np.squeeze( tsk["KE_hor_x"][:])
                uy_data = np.squeeze( tsk["KE_hor_y"][:])


                # storing
                self.u_x_data[ind * max_write : (ind + 1)* max_write] = ux_data
                self.u_z_data[ind * max_write : (ind + 1)* max_write] = uz_data
                self.u_y_data[ind * max_write : (ind + 1)* max_write] = uy_data




    def analyze_velocity(self, h_size, box_dimension = 3, velocity_range = "all", average_method = "max"):
        self.u_rms_all = np.zeros(np.shape(h_size))
        self.u_rms_x = np.zeros(np.shape(h_size))
        self.u_rms_y = np.zeros(np.shape(h_size))
        self.u_rms_z = np.zeros(np.shape(h_size))
        self.tflux_avg = np.zeros(np.shape(h_size))

        for i in range(len(h_size)):
            ind = int( (1-h_size[i]) * len(self.z_coord)) 
            if h_size[i] == 0:
                ind = 0

            # cubic = scipy.interpolate.CubicSpline(self.z_coord, self.t_flux_data[i])
            # fine_grid =  np.linspace(self.z_coord[0], self.z_coord[-1], 3000)
            # fine_profile = cubic(fine_grid)
            # ind_fine =  np.where(fine_grid < (1- h_size[i]))[-1] + 1


            if box_dimension == 2:
                u_mag_squared = self.u_x_data[i][ind:] + self.u_z_data[i][ind:]
            if box_dimension == 3:
                
                u_mag_squared_all = self.u_x_data[i][ind:] + self.u_y_data[i][ind:] + self.u_z_data[i][ind:]
                u_mag_squared_z = self.u_z_data[i][ind:]
                u_mag_squared_x = self.u_x_data[i][ind:]
                u_mag_squared_y = self.u_y_data[i][ind:]

            if average_method == "avg":
                u_rms_all = scipy.integrate.simpson( u_mag_squared_all, x = self.z_coord[ind:])/ h_size[i]
                u_rms_x = scipy.integrate.simpson( u_mag_squared_x, x = self.z_coord[ind:])/ h_size[i]
                u_rms_y = scipy.integrate.simpson( u_mag_squared_y, x = self.z_coord[ind:])/ h_size[i]
                u_rms_z = scipy.integrate.simpson( u_mag_squared_z, x = self.z_coord[ind:])/ h_size[i]


                flux_integral =  scipy.integrate.simpson( self.t_flux_data[i][ind:], x = self.z_coord[ind:])
                tflux_avg = flux_integral / h_size[i]

            if average_method == "max":
                u_rms_all = np.max(u_mag_squared_all)
                u_rms_x = np.max(u_mag_squared_x)
                u_rms_y = np.max(u_mag_squared_y)
                u_rms_z = np.max(u_mag_squared_z)
                tflux_avg = np.max(self.t_flux_data[i][ind:])


            self.u_rms_all[i] = np.sqrt(u_rms_all)
            self.u_rms_x[i] = np.sqrt(u_rms_x)
            self.u_rms_y[i] = np.sqrt(u_rms_y)
            self.u_rms_z[i] = np.sqrt(u_rms_z)
            self.tflux_avg[i] = tflux_avg

        print("root mean square velocity computed")


    # PART 2: analyze data

    # Analyze is a series of function to compute and measure the size of convection zone h(t)
    # There're 4 ways to do it: 
    # 1. looking at concertration, where it drops, i.e by 3%, is th h(t)
    # 2. the peak in the conpositional flux is roughly the size of h(t)
    # 3. similar for density with 1., where it drops is the h(t)
    # 4. solute conservation, ds/dz(method see function) is a way, but terribly behaves when later time
    def analyze_c(self, method = "heating", boundary = "relative", 
                  relative_percenta = 0.05, absolute_percenta = 0.001):

        self.h_c = np.zeros(len(self.time))

        # mode 01: analyze composition
        # h is the first point, which difference with h[0] exceed tol 
        if method == "heating":
            for i in range(len(self.time)):

                if boundary == "relative":
                    bound = self.c_data[i][0]
                    tol = bound * relative_percenta
                elif boundary == "absolute":
                    tol = absolute_percenta 

                # interpolation
                cubic = scipy.interpolate.CubicSpline(self.z_coord, self.c_data[i])
                fine_grid =  np.linspace(self.z_coord[0], slef.z_coord[-1], 5000)
                fine_profile = cubic(fine_grid)

                ind = np.where(abs(fine_profile - bound) > tol)
                record = fine_grid[ind]
                try:
                    self.h_c[i] = record[0]
                except IndexError:
                    print("convection zone reached boudary")

        elif method == "cooling":
            for i in range(len(self.time)):
                if boundary == "relative":
                    bound = self.c_data[i][-1]
                    tol = bound * relative_percenta
                elif boundary == "absolute":
                    tol = absolute_percenta 
                
                # interpolation
                cubic = scipy.interpolate.CubicSpline(self.z_coord, self.c_data[i])
                fine_grid =  np.linspace(self.z_coord[0], self.z_coord[-1], 5000)
                fine_profile = cubic(fine_grid)

                ind = np.where(abs(fine_profile - bound) > tol)
                record = fine_grid[ind]
                try:
                    self.h_c[i] =  1 - record[-1]
                except IndexError:
                    print("convection zone reached boudary")

        print("... ...")
        print ("Analyze with composition completed")
                    

                    
    def analyze_t(self, method = "heating", boundary = "relative", 
                  relative_percenta = 0.03, absolute_percenta = 0.001):

        self.h_t = np.zeros(len(self.time))

        # mode 01: analyze composition
        # h is the first point, which difference with h[0] exceed tol 
        if method == "heating":
            for i in range(len(self.time)):

                if boundary == "relative":
                    bound = self.t_data[i][0]
                    tol = bound * relative_percenta
                elif boundary == "absolute":
                    tol = absolute_percenta 
                
                ind = np.where(abs(self.t_data[i] - bound) > tol)
                record = self.z_coord[ind]
                try:
                    self.h_t[i] = record[0]
                except IndexError:
                    print("convection zone reached boudary")

        elif method == "cooling":
            for i in range(len(self.time)):
                if boundary == "relative":
                    bound = self.t_data[i][-1]
                    tol = bound * relative_percenta
                elif boundary == "absolute":
                    tol = absolute_percenta 
                
                ind = np.where(abs(self.t_data[i] - bound) > tol)
                record = self.z_coord[ind]
                try:
                    self.h_t[i] =  1 - record[-1]
                except IndexError:
                    print("convection zone reached boudary")

        print("... ...")
        print ("Analyze with temp completed")
                    

    def analyze_cflux(self, method = "heating"):

        self.h_cflux = np.zeros(len(self.time))
        
        # mode 02: compare composition flux
        # Composition flux in a gaussian profile, 
        # with highest peak located at h 

        if method == "heating":
            for i in range(len(self.time)):
                peak, _ = find_highest_peak(self.c_flux_data[i])
                self.h_cflux[i]=  self.z_coord[peak]

        elif method == "cooling":
            for i in range(len(self.time)):
                peak, _ = find_highest_peak(self.c_flux_data[i])
                self.h_cflux[i]=  1  - self.z_coord[peak]

        print ("Analyze with compositional flux completed")


    
    def analyze_rho(self, method = "heating", boundary = "relative", 
                    relative_percenta = 0.03, absolute_percenta = 0.001):

        self.h_rho = np.zeros(len(self.time))

        # mode 03: analyze density
        # since rho = rho* alpha (C - T), T in unit of beta/alpha * S
        # h is the first point, which difference with h[0] exceed tol 
        # where tol is 10% of the plateau in the profile 
        if method == "heating":
            
            for i in range(len(self.time)):
                bound = np.max(self.rho_data[i])

                if boundary == "relative":
                    tol = bound * relative_percenta
                elif boundary == "absolute":
                    tol = absolute_percenta 
                
                max_ind = np.where(self.rho_data[i] == bound)[0]

                # to avoid the dip in profile at boundary, flattern them with max value
                # to avoid changing th data, use a copy to flattern
                rho_copy = self.rho_data.copy()
                rho_list_copy = rho_copy[i]
                rho_list_copy[:int(max_ind)] = bound

                # interpolation
                cubic = scipy.interpolate.CubicSpline(self.z_coord, rho_list_copy)
                fine_grid =  np.linspace(self.z_coord[0], self.z_coord[-1], 5000)
                fine_profile = cubic(fine_grid)

                ind = np.where(abs(fine_profile - bound) > tol)
                record = fine_grid[ind] 
                try:
                    self.h_rho[i] = record[0]
                except IndexError:
                    print("convection zone reached boudary")


        if method == "cooling":
            
            for i in range(len(self.time)):
                bound = np.min(self.rho_data[i])

                if boundary == "relative":
                    tol = bound * relative_percenta
                elif boundary == "absolute":
                    tol = absolute_percenta 
                
                min_ind = np.where(self.rho_data[i] == bound)[0]

                # to avoid the dip in profile at boundary, flattern them with max value
                # to avoid changing th data, use a copy to flattern
                rho_copy = self.rho_data.copy()
                rho_list_copy = rho_copy[i]
                rho_list_copy[int(min_ind):] = bound

                # interpolation
                cubic = scipy.interpolate.CubicSpline(self.z_coord, rho_list_copy)
                fine_grid =  np.linspace(self.z_coord[0], self.z_coord[-1], 5000)
                fine_profile = cubic(fine_grid)

                ind = np.where(abs(fine_profile - bound) > tol)
                record = fine_grid[ind] 
                try:
                    self.h_rho[i] = 1 - record[-1]
                except IndexError:
                    print("convection zone reached boudary")


        print("... ...")
        print ("Analyze with density completed")



    def analyze_dsdz(self, method = "heating", relative_percenta = 0.05, selective_percenta = 0.5):
            
        self.h_dsdz =  np.zeros(len(self.time))

        # mode 04: use the relation that: Delta_C = 0.5 dC_0/dz h
        # for linear initial profile, where we have dC_0/dz=1, then 
        # relation is that: h = 2 * Delta_C, where delta_C is the c of the convection zone
        if method == "heating":
            for i in range(len(self.time)):
                bound = self.c_data[i][0]
                tol = bound * relative_percenta
                ind = np.where(abs(self.c_data[i] - bound) > tol)
                record = self.z_coord[ind]
                try:
                    ds_ave = self.c_data[i][: int(record[0]* selective_percenta)]
                    self.h_dsdz[i] = 2*(1 - np.mean(ds_ave))
                except IndexError:
                    print("convection zone reached boudary")

        if method == "cooling":
            for i in range(len(self.time)):
                bound = self.c_data[i][0]
                tol = bound * relative_percenta
                ind = np.where(abs(self.c_data[i] - bound) > tol)
                record = self.z_coord[ind]
                try:
                    ds_ave = self.c_data[i][: int(record[0]* selective_percenta)]
                    self.h_dsdz[i] = 2* np.mean(ds_ave)
                except IndexError:
                    print("convection zone reached boudary")

        print("... ...")
        print ("Analyze with solute conservation completed")




    def compute_grav_eng(self, method = "cooling", region = "full", h_profile = 0):
        # region = "full"  or "h_only"
        # gravitational energy: E/g = \int_0^h \delta_rho * z * dz
        # dimension of the quanity: g * \beta * S_0 * H^3 * \rho_0
        # prefix: Ra * Pr 
        self.grav_eng = np.zeros(len(self.time))

        if region == "full":
            for i in range(len(self.time)):
                if method == "cooling":
                    delta_rho = self.rho_data[i]  - (1 - self.z_coord)
                if method == "heating":
                    delta_rho = self.rho_data[i] - (1 - self.z_coord)

                integrand =  delta_rho * self.z_coord
                E_g = scipy.integrate.simpson(y = integrand, x = self.z_coord)

                self.grav_eng[i] = E_g * self.Ra_value * self.prandtl


        if region == "h_only":
            for i in range(len(self.time)):
                if method == "cooling":
                    delta_rho = self.rho_data[i] - (1 - self.z_coord)
                if method == "heating":
                    delta_rho = self.rho_data[i] - (1 - self.z_coord)

                ind = np.where( h_profile[i] > self.z_coord )

                integrand =  delta_rho[ind] * self.z_coord[ind]
                E_g = scipy.integrate.simpson(y = integrand, x = self.z_coord[ind])

                self.grav_eng[i] = E_g * self.Ra_value * self.prandtl


        return self.grav_eng





    def compute_thermal_eng(self, method = "cooling", region = "full", h_profile = 0):
        # thermal energy: E/ (c_p * \rho_0) = \int_0^h \Delta_t * dz
        # dimension: \rho_0 * c_P * H * T_0 = F_c * k_T / H^2
        # prefix = Ra * Pr / lambda
        self.therm_eng = np.zeros(len(self.time))
        prefix =  self.Ra_value * self.prandtl / self.lamb

        if region == "full":
            for i in range(len(self.time)):
                if method == "cooling":
                    integrand =  1 - self.t_data[i]
                    E_t = scipy.integrate.simpson(y = integrand, x = self.z_coord)
                self.therm_eng[i] = E_t

        if region == "h_only":
            for i in range(len(self.time)):
                ind = np.where( h_profile[i] > self.z_coord )
                if method == "cooling":
                    integrand = (1 - self.t_data[i]) * self.rho_data[i]
                    E_t = scipy.integrate.simpson(y = integrand[ind], x = self.z_coord[ind])

                self.therm_eng[i] = E_t

        return self.therm_eng



    def compute_kinetic_eng(self, method = "cooling", region = "full", h_profile = 0):

        self.kin_eng = np.zeros(len(self.time))
        prefix =  1

        if region == "full":
            for i in range(len(self.time)):
                if method == "cooling":
                    integrand =  self.u_z_data[i] +  self.u_x_data[i]  +  self.u_y_data[i] 
                    E_k = scipy.integrate.simpson(y = integrand, x = self.z_coord)

                self.kin_eng[i] = E_k /2

        if region == "h_only":
            for i in range(len(self.time)):
                ind = np.where( h_profile[i] > self.z_coord )
                if method == "cooling":
                    integrand = self.u_z_data[i] +  self.u_x_data[i]  +  self.u_y_data[i] 
                    E_k = scipy.integrate.simpson(y = integrand[ind], x = self.z_coord[ind])

                self.kin_eng[i] = E_k /2 

        return self.kin_eng





    def compute_surface_flux(self, function="deflux"):
        prefix = self.F_ratio * self.Ra_value / self.lamb
        self.surf_eng = np.zeros(len(self.time))

        if function == "deflux":
            # func = self.t_0 / (self.t_0 + self.time) 
            func = self.t_0 * np.log((self.t_0 + self.time)/ self.t_0 )
        
        if function == "constant":
            # func = np.ones(len(self.time))
            func = self.time
        
        for i in range(len(self.time)):
            E_surf = scipy.integrate.simpson(y = func[:i+1], x = self.time[:i+1])
            self.surf_eng[i] = func[i] * self.F_ratio

        return self.surf_eng




    # computation of some constant
    # epsilon is the thermal coeficient 

    def average_all_h(self):
        self.h_avg = (self.h_c + self.h_cflux + self.h_rho) / 3

        return self.h_avg


    def compute_epsilon(self, h):

        self.epsilon_measured = np.zeros(len(self.time))


        for i in range(len(self.time)):
            cubic = scipy.interpolate.CubicSpline(self.z_coord, self.t_data[i])
            fine_grid =  np.linspace(self.z_coord[0], self.z_coord[-1], 2000)
            fine_profile = cubic(fine_grid)

            ind = np.where(fine_grid > 1 - h[i])
            if h[i] == 0.0: 
                ind = np.where(fine_grid > 0)


            tot_fluxeng = scipy.integrate.simpson(y = 1 - fine_profile, x = fine_grid)
            convect_fluxeng = scipy.integrate.simpson(y = 1 - fine_profile[ind], x = fine_grid[ind])

            epi_ratio = (tot_fluxeng - convect_fluxeng) / tot_fluxeng

            self.epsilon_measured[i] = epi_ratio


        return self.epsilon_measured



    def compute_gamma(self, h):
        adv = np.roll(h, -1)
        back = np.roll(h, 1)
        adv[-1] = adv[-2]
        back[0] = 0 

        self.dhdt = (adv - back) / (2 * (self.time[1] - self.time[0]))
        self.gamma_measured = np.zeros(len(self.time))


        for i in range(len(self.time)):

            cubic = scipy.interpolate.CubicSpline(self.z_coord, 1 - self.t_data[i])
            fine_grid =  np.linspace(self.z_coord[0], self.z_coord[-1], 2000)
            fine_profile = cubic(fine_grid)

            ind = np.where(fine_grid > 1 - h[i])
            if h[i] == 0.0: 
                ind = np.where(fine_grid > 0)

            convect_fluxeng = scipy.integrate.simpson(y = fine_profile[ind], x = fine_grid[ind])

            if h[i] == 0.0:
                del_t = convect_fluxeng 
                del_s = 1/2

            else:
                del_t = convect_fluxeng / h[i]
                print(del_t)
                del_s = h[i]/2



            gamma = (1 / self.F_ratio ) * ( del_s - del_t) * (1/ self.t_flux_data[i][-1]) * self.dhdt[i]



            self.gamma_measured[i] = gamma


        return self.gamma_measured





    # PART 3: fitting functions 
    def fit_for_h (self, function_to_fit, outcome_domain, c = False, c_flux = False,
                    rho =  False, dsdz = False, stacking_results = False, t_i = 0):

        # for each measured h(t) in different ways, function is to fit sperately,
        # False indicates not to be included in the fitting
        # fitting is done with all H(t) contributed, not an average of individual contribution
        # specify t_i, h_i will be computed from h profile


        # initialize or stacking depends on the saetting 
        if stacking_results:
            print("Notice, results are stacking")
            print("Already fitted for:", self.fit_type)
        else:
            self.fit_type = []
            self.fit_param = []

        fit_data_pool = []
        fit_time_pool = []

        if c:
            fit_data_pool = fit_data_pool + list(self.h_c)
            fit_time_pool = fit_time_pool + list(self.time)
        else:
            None
        if c_flux:
            fit_data_pool = fit_data_pool + list(self.h_cflux)
            fit_time_pool = fit_time_pool + list(self.time)
        else:
            None
        if rho:
            fit_data_pool = fit_data_pool + list(self.h_rho)
            fit_time_pool = fit_time_pool + list(self.time)
        else:
            None


        # find h_i with t_i specified:
        time_pool = np.array(fit_time_pool)
        data_pool = np.array(fit_data_pool)
        cubic = scipy.interpolate.CubicSpline(time_pool, data_pool)

        fine_t = np.linspace(t_i, self.time[-1], 100)
        fine_h = cubic(fine_t)
        h_i = fine_h[0]
        ind = np.where(data_pool> h_i)

        # function to fit is an interger number
        # 1: fit through square root function: h~t^1/2
        # 2: fit through power law: h~t^a, reporting a
        # 3: fit log function, without roataion, fit for k
        # 4: 
        # 51: numerically fitting for constant flux, no rotaion
        # 52: numerically fitting for constant flux, rotaion
        # 53: numerically fitting for decrease flux, no rotaion
        # 54: numerically fitting for decrease flux, rotaion


        if function_to_fit == 1:
            fit_param, outcome = t_06(time_pool[ind], data_pool[ind], outcome_domain, t_i, h_i)
            self.fit_param.append(fit_param)
            self.fit_type.append("t^0.5")

        if function_to_fit == 2:
            fit_param, outcome = t_power(time_pool[ind], data_pool[ind],outcome_domain, t_i, h_i)
            self.fit_param.append(fit_param)
            self.fit_type.append("power law")

        if function_to_fit == 3:
            fit_param, outcome = h_nr_fitk(time_pool[ind], data_pool[ind], outcome_domain, self.t_0, self.F_ratio, t_i, h_i)
            self.fit_param.append(fit_param)
            self.fit_type.append("log func_nr")

        if function_to_fit == 4:
            fit_param, outcome = h_wr_fitk(time_pool[ind], data_pool[ind], outcome_domain, 
                                        self.t_0, self.F_ratio, self.Ra_value, self.prandtl, self.omega_num, t_i, h_i)
            self.fit_param.append(fit_param)
            self.fit_type.append("log func_wr")

        if function_to_fit == 4:
            fit_param, outcome = h_wr_fitk(time_pool[ind], data_pool[ind], outcome_domain, 
                                        self.t_0, self.F_ratio, self.Ra_value, self.prandtl, self.omega_num, t_i, h_i)
            self.fit_param.append(fit_param)
            self.fit_type.append("log func_wr")

            
        if function_to_fit == 51:
            fit_param, outcome = numer_fit_func(time_pool[ind], data_pool[ind], outcome_domain, 
                                        self.t_0, self.F_ratio, self.Ra_value, self.prandtl, self.omega_num, t_i, h_i,
                                        fit = "constant_nr")
            self.fit_param.append(fit_param)
            self.fit_type.append("numer_cnr")

        if function_to_fit == 52:
            fit_param, outcome = numer_fit_func(time_pool[ind], data_pool[ind], outcome_domain, 
                                        self.t_0, self.F_ratio, self.Ra_value, self.prandtl, self.omega_num, t_i, h_i,
                                        fit = "constant_rot")
            self.fit_param.append(fit_param)
            self.fit_type.append("numer_crot")
            
        if function_to_fit == 53:
            fit_param, outcome = numer_fit_func(time_pool[ind], data_pool[ind], outcome_domain, 
                                        self.t_0, self.F_ratio, self.Ra_value, self.prandtl, self.omega_num, t_i, h_i,
                                        fit = "deflux_nr")
            self.fit_param.append(fit_param)
            self.fit_type.append("numer_denr")

        if function_to_fit == 54:
            fit_param, outcome = numer_fit_func(time_pool[ind], data_pool[ind], outcome_domain, 
                                        self.t_0, self.F_ratio, self.Ra_value, self.prandtl, self.omega_num, t_i, h_i,
                                        fit = "deflux_rot")
            self.fit_param.append(fit_param)
            self.fit_type.append("numer_derot")

        print("... ...")
        print("Current fit pool:", self.fit_type, self.fit_param)
        return outcome






    def fit_param_epi_gamma(self, hdata, h0,  t_domain, function_type = 1):

        F_ratio = self.F_ratio
        ds0dz = 1
        t_0 = t_domain[0]

        def f_th_deflux(t, h, epsilon, gamma):
            g_th = 0.5 * ds0dz * h - (1-epsilon) * F_ratio * t_0 * np.log(abs ((t_0 + t) / t_0) ) / h
            f = gamma * F_ratio * t_0 / (t_0 + t) 
            return f / g_th

        def f_th_constant_flux(t, h, epsilon, gamma):
            g_th = 0.5 * ds0dz * h - (1-epsilon) * F_ratio * t / h
            f = gamma * F_ratio 
            return f / g_th

        
        if function_type == 1: 
                
            def func_numer(t, epsilon, gamma):

                dt = np.diff(t)
                h = np.zeros(len(t))
                h[0] = h0

                for i in range(len(t)-1):
                    h1 = h[i] + 0.5*dt[i] * f_th_constant_flux(t[i], h[i], 0.3, gamma)
                    f1 = f_th_constant_flux(t[i]+ dt[i]/2, h1, 0.3, gamma)

                    h[i+1] = h[i] + f1*dt[i]
                return h

        if function_type == 2:
                
            def func_numer(t, epsilon, gamma):

                dt = np.diff(t)
                h = np.zeros(len(t))
                h[0] = h0

                for i in range(len(t)-1):
                    h1 = h[i] + 0.5*dt[i] * f_th_deflux(t[i], h[i], epsilon, gamma)
                    f1 = f_th_deflux(t[i]+ dt[i]/2, h1, epsilon, gamma)

                    h[i+1] = h[i] + f1*dt[i]
                return h

        if function_type == 3:
            def func_numer(t, epsilon, gamma):
                c = 1-  epsilon  +2*gamma
                prefix = np.sqrt(2*c)
                h = prefix * F_ratio**0.5 * np.sqrt(t)
                return h


        fit, _ = scipy.optimize.curve_fit(func_numer, self.time, hdata)

        print("epsilon:", fit[0])
        print("gamma:", fit[1])

        # return func_numer(t_domain, fit[0], fit[1])

        return func_numer(t_domain, 0.8, 0.4)


