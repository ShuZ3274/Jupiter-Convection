import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import ticker
import matplotlib as mpl
import pickle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch



with open('test_analysis.pkl', 'rb') as inp:
    test_1, test_2, test_3, test_4, test_5 = pickle.load(inp)



# modes are 
# "single_profile": plot for the single test, the profiles and h by differnt measurement
# "profiles": plot profiles like: c, cfluc, t, t_flux, etc.
# "analysis": different atype of analysis, data to look at
# "h": size of the convection zone for all chosen test case
mode =  input("mode:")
tests_include = [1, 2, 3, 4, 5]
save_fig = False
log_scale = False



# single profile
if mode == "single_profile":
    test_case =  test_3
    down_factor = 8
    pick = 2


    fig, ax= plt.subplots(1,2, figsize=(15, 6))
    shrinked_num = len(test_case.time) // down_factor
    color = mpl.colormaps["viridis"](np.linspace(0,1, shrinked_num))

    for i in range(shrinked_num):
        velo = np.sqrt(np.abs(test_case.u_z_data[i * down_factor]))
        ax[0].plot(test_case.z_coord, velo, color = color[i],  alpha = 0.10 + i* (0.7/ (shrinked_num - 1)))

        ax[0].vlines(x = 1 -  test_case.h_cflux[i * down_factor], ymin = 0, ymax = 0.5 * np.max(velo)
                , color = 'red',  alpha = 0.20 + i * (0.7/ (shrinked_num - 1)), label = "flux_peak" )
        ax[0].vlines(x = 1 -  test_case.h_rho[i * down_factor], ymin = 0, ymax = 0.5 * np.max(velo)
                , color = 'darkgrey',  alpha = 0.20 + i * (0.7/ (shrinked_num - 1)), label = "flux_peak" )

    ax[0].set_xlabel("z_coord")
    ax[0].set_ylabel("U_z")

    ax[1].plot(test_case.z_coord, test_case.c_data[pick * down_factor], linestyle = "--", color = color[pick], label = 'c' )
    ax[1].plot(test_case.z_coord, test_case.rho_data[pick * down_factor], color = color[pick], label = 'rho' )
    ax[1].plot(test_case.z_coord, test_case.c_flux_data[pick * down_factor] / np.max(test_case.c_flux_data[pick * down_factor]),
            linestyle = 'dotted', color = color[pick], label = 'c_flux' )
    ax[1].plot(test_case.z_coord, test_case.t_flux_data[pick * down_factor],linestyle = 'dashdot', color = color[pick], label = 't_flux' )

    ax[1].vlines(x = 1 -  test_case.h_cflux[pick * down_factor], ymin = 0, ymax = 1
            , color = 'red',  alpha = 0.20 + pick * (0.7/ (shrinked_num - 1)), label = "flux_peak" )
    ax[1].vlines(x = 1 -  test_case.h_c[pick * down_factor], ymin = 0, ymax = 1
            , color = 'blue',  alpha = 0.20 + pick * (0.7/ (shrinked_num - 1)), label = r"$5\%c$" )
    ax[1].vlines(x = 1 -  test_case.h_rho[pick * down_factor], ymin = 0, ymax = 1
            , color = 'darkgrey',  alpha = 0.20 + pick * (0.7/ (shrinked_num - 1)), label = "rho" )


    ax[1].set_xlabel("z_coord")
    ax[1].set_ylabel("normalized_data")
    plt.legend()
    plt.show()

    if save_fig == True:
        fig.savefig('single_scale_test3.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)








# profiles

if mode == "profiles":
    

    fig, ax= plt.subplots(2,2, figsize=(15, 10))

    tot = np.array([test_1, test_2, test_3, test_5])
    tot_name = np.array(["deflux_rot", "deflux_nr", "con_nr", "con_rot"])

    
    for j , test_case in enumerate(tot):

        color = mpl.colormaps["viridis"](np.linspace(0,1, len(test_case.time)))
        object_int = test_case.t_data

        for i in range(len(test_case.time)):
            ax[j//2, j%2].plot(test_1.z_coord, object_int[i], color = color[i],  alpha = 0.10 + i* (0.7/ (len(test_case.time) - 2)))

            ax[j//2, j%2].vlines(x =1 -  test_case.h_rho[i], ymin = 0.7 * np.max(object_int), ymax = np.max(object_int)
                    , color = 'red',  alpha = 0.10 + i* (0.7/ (len(test_case.time) - 2)) )

            ax[j//2, j%2].set_xlabel("z_coord")
            ax[j//2, j%2].set_ylabel(rf"$F_H$")
            ax[j//2, j%2].set_title(tot_name[j])

            ax[j//2, j%2].plot(1- test_case.h_rho, test_case.tflux_avg, color = "black")

            if log_scale ==  True:
                ax[j//2, j%2].set_yscale('log')
                ax[j//2, j%2].set_xscale('log')
    plt.show()

    if save_fig == True:
        fig.savefig('Flux_th.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)




# temp001

if mode == "temp001":
    
    fig, ax= plt.subplots(3,1, figsize=(5, 10))

    test_int = test_1


    ax[0].plot( test_int.time, test_int.u_rms_all,  "-o", markersize = 4, color = 'deeppink', label = rf"rms")
    ax[0].plot( test_int.time , test_int.u_rms_x, "-s", markersize = 4, color = 'darkgreen', label = rf"x")
    ax[0].plot( test_int.time, test_int.u_rms_y, "-o", markersize = 4, color = 'darkgrey', alpha = 0.7, label = rf"y")
    ax[0].plot( test_int.time, test_int.u_rms_z, "-o", markersize = 4, color = 'black', label = rf"z")
    ax[0].set_ylabel(r"velo", fontsize=12)


    ax[1].plot( test_int.time, test_int.h_c,  "-o", markersize = 4, color = 'deeppink', label = rf"c")
    ax[1].plot( test_int.time , test_int.h_cflux, "-s", markersize = 4, color = 'darkgreen', label = rf"cflux")
    ax[1].plot( test_int.time , test_int.tflux_avg, "--", markersize = 4, color = 'grey', label = rf"tflux")
    ax[1].set_ylabel(r"$F_H$", fontsize=12)

    fh = test_int.tflux_avg **2

    fh_dummy =  np.logspace(-5, 1, 2000)
    u_scale_line_1 =  test_int.u_rms_all[-10] * fh_dummy ** (1/3)* 5
    u_scale_line_2 =  test_int.u_rms_all[-10] * fh_dummy ** (1/5) 


    ax[2].plot( test_int.h_rho[:-2] * fh[:-2] , test_int.u_rms_z[:-2]  ,  "-o", markersize = 4, color = 'deeppink', label = rf"uz")
    ax[2].plot( test_int.h_rho[:-2] * fh[:-2] , test_int.u_rms_all[:-2] , "-s", markersize = 4, color = 'darkgreen', label = rf"rms")
    ax[2].plot( fh_dummy, u_scale_line_1, "--", color = "darkgrey", alpha = 0.7, label = "1/3")
    ax[2].plot( fh_dummy, u_scale_line_2, color = "darkgrey", alpha = 0.7, label = "1/5")

    ax[2].set_xlabel(r"hF_H", fontsize=12)
    ax[2].set_ylabel(r"$velo$", fontsize=12)

    ax[2].set_xscale("log")
    ax[2].set_yscale("log")
    ax[2].set_ylim(ymax = 1e3, ymin =1e2)
    ax[2].set_xlim(xmax = 1e-2, xmin =1e-4)

    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.grid()
    plt.show()
    if save_fig ==  True:
        fig.savefig('velo_h_ross_test1.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)



"hfu relation"
if mode == "hfu":

    fig, ax= plt.subplots(2,1, figsize=(8, 10))


    fh_dummy =  np.logspace(-5, 1, 2000)
    u_scale_line_1 =  test_1.u_rms_all[-10] * fh_dummy ** (1/3)
    u_scale_line_2 =  test_1.u_rms_all[-10] * fh_dummy ** (1/5) 

    ax[0].plot(test_1.h_rho * test_1.tflux_avg **2, test_1.u_rms_z, label = "dwr")
    ax[0].plot(test_2.h_rho * test_2.tflux_avg, test_2.u_rms_z, label = "dnr")
    ax[0].plot(test_3.h_rho * test_3.tflux_avg, test_3.u_rms_z, label = "cnr")
    ax[0].plot(test_4.h_rho * test_4.tflux_avg **2, test_4.u_rms_z, label = "cwr10")
    ax[0].plot(test_5.h_rho * test_5.tflux_avg **2, test_5.u_rms_z, label = "cwr5")
    ax[0].set_xlabel(r"hF_H", fontsize=12)
    ax[0].set_ylabel(r"$U_z$", fontsize=12)
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_ylim(ymax = 1e3, ymin =1e2)
    ax[0].set_xlim(xmax = 1, xmin =1e-3)

    ax[1].plot(test_1.h_rho * test_1.tflux_avg **2, test_1.u_rms_all, label = "dwr")
    ax[1].plot(test_2.h_rho * test_2.tflux_avg, test_2.u_rms_all, label = "dwr")
    ax[1].plot(test_3.h_rho * test_3.tflux_avg, test_3.u_rms_all, label = "dwr")
    ax[1].plot(test_4.h_rho * test_4.tflux_avg **2, test_4.u_rms_all, label = "dwr")
    ax[1].plot(test_5.h_rho * test_5.tflux_avg **2, test_5.u_rms_all, label = "dwr")
    ax[1].set_xlabel(r"hF_H", fontsize=12)
    ax[1].set_ylabel(r"$U_{rms}$", fontsize=12)
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    ax[1].set_ylim(ymax = 1e3, ymin =1e2)
    ax[1].set_xlim(xmax = 1, xmin =1e-3)

    ax[0].plot( fh_dummy, u_scale_line_1*2, "--", color = "darkgrey", alpha = 0.7, label = "1/3")
    ax[0].plot( fh_dummy, u_scale_line_2, color = "darkgrey", alpha = 0.7, label = "1/5")
    ax[1].plot( fh_dummy, u_scale_line_1*2, "--", color = "darkgrey", alpha = 0.9, label = "1/3")
    ax[1].plot( fh_dummy, u_scale_line_2*1.5, color = "darkgrey", alpha = 0.7, label = "1/5")


    ax[0].legend()

    plt.grid()
    plt.show()
    if save_fig ==  True:
        fig.savefig('hfu.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)






# "analysis"
if mode == "analysis":
    mpl.style.use('seaborn-v0_8')
    fig2, ax = plt.subplots(figsize=(8, 6))
    interest_data = "rossby"

    if interest_data == "rms":
        ax.plot( test_1.h_c, test_1.u_rms,  "-o", markersize = 4, color = 'deeppink', label = rf"deflux_rotation")
        ax.plot( test_2.h_c, test_2.u_rms, "-s", markersize = 4, color = 'darkblue', label = rf"deflux_no rot")
        ax.plot( test_3.h_c, test_3.u_rms, "-s", markersize = 4, color = 'darkgreen', label = rf"constant")
        ax.plot( test_4.time, test_4.u_rms, "-o", markersize = 4, color = 'darkgrey', label = rf"constant_rot_10")
        ax.plot( test_5.time, test_5.u_rms, "-o", markersize = 4, color = 'black', label = rf"constant_rot_5")
        if log_scale == True:
            plt.xscale("log")
            plt.yscale("log")
        ax.set_xlabel(r"time", fontsize=12)
        ax.set_ylabel(r"$rms$", fontsize=12)
        # plt.ylim(ymax = 1e4, ymin = 1)
        # plt.xlim(xmin = 1e-2, xmax = 1)
        ax.grid(which='major', alpha=0.85)
        ax.grid(which='minor', alpha=0.5)



    if interest_data == "rossby":

        ax.plot( test_1.time, test_1.u_rms_z * 3e-6 / 0.1 / test_1.h_rho,  "-o", markersize = 4, color = 'deeppink', label = rf"deflux_rotation")
        ax.plot( test_4.time, test_4.u_rms_z * 3e-6 / 0.1 / test_4.h_rho, "-o", markersize = 4, color = 'darkgrey', label = rf"constant_rot_10")
        ax.plot( test_5.time, test_5.u_rms_z * 3e-6 / 0.1 / test_5.h_rho, "-o", markersize = 4, color = 'black', label = rf"constant_rot_5")
        if log_scale == True:
            plt.xscale("log")
            plt.yscale("log")
        ax.set_xlabel(r"time", fontsize=12)
        ax.set_ylabel(r"$Ro$", fontsize=12)
        ax.grid(which='major', alpha=0.85)
        ax.grid(which='minor', alpha=0.5)

    plt.legend()
    plt.grid()
    plt.show()
    if save_fig ==  True:
        fig2.savefig('rossby_uz.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)







"h"
if mode == "h":
    t_i = 0.005
    time = np.logspace(-6, 1, 1000)
    time_short = np.logspace(np.log10(t_i), 1, 1000)


    fit_2 = test_1.fit_for_h(function_to_fit = 54, outcome_domain = time_short, c_flux = True, t_i  = t_i)
    fit_3 = test_2.fit_for_h(function_to_fit = 53, outcome_domain = time_short, c_flux = True, t_i  = t_i)
    fit_1 = test_3.fit_for_h(function_to_fit = 51, outcome_domain = time_short, c_flux = True, t_i  = t_i)
    fit_5 = test_5.fit_for_h(function_to_fit = 2, outcome_domain = time_short, c = True, t_i  = t_i)


    mpl.style.use('seaborn-v0_8')
    fig2, ax = plt.subplots(figsize=(8.5, 6))


    ax.vlines(x = test_1.t_0, ymin = 1e-6, ymax = 1, color = 'plum', linewidth = 0.9, alpha = 0.8)
    ax.text(test_1.t_0 * 1.01, 0.5, rf'$t_0$ ={test_1.t_0}', fontsize = 9, color = 'plum', alpha = 1, rotation = - 90)
    ax.vlines(x = test_2.t_0, ymin = 1e-6, ymax = 1, color = 'blue', linewidth = 0.9, alpha = 0.3)


    ax.plot(time_short, fit_2, alpha = 0.6, color = 'deeppink', linestyle = "--")
    ax.plot(time_short, fit_1, alpha = 0.2, color = 'darkgreen', linestyle = "--")
    ax.plot(time_short, fit_3, alpha = 0.6, color = 'darkblue', linestyle = "--")


    ax.plot(test_1.time, test_1.h_rho, "-o", markersize = 4, color = 'deeppink', label = rf"deflux_rotation")
    ax.plot(test_2.time, test_2.h_rho, "-s", markersize = 4, color = 'darkblue', label = rf"deflux_no rot")
    ax.plot(test_3.time, test_3.h_rho, "-s", markersize = 4, color = 'darkgreen', label = rf"constant")
    ax.plot(test_4.time, test_4.h_rho, "-s", markersize = 4, color = 'darkgrey', label = rf"constant_rot_10")
    ax.plot(test_5.time, test_5.h_rho, "-s", markersize = 4, color = 'black', label = rf"constant_rot_5")



    ax.tick_params(axis='y', direction = "in", right = True, length = 10, which = "major")
    ax.tick_params(axis='x',direction='in',top=True, length = 10, which = "major")
    ax.tick_params(axis='y', direction = "in", right = True, length = 4, which = "minor")
    ax.tick_params(axis='x',direction='in', length = 4, which = "minor")



    # ax.text(0.052, 0.47, f'k= {round(test_1.fit_param[0][0], 2)}', fontsize = 11, color = 'deeppink')
    # ax.text(0.052, 0.22, f'k= {round(test_2.fit_param[0][0], 2)}', fontsize = 11, color = 'darkblue')

    # ax.text(0.005, 0.045, r'Fit: $h^2(t) = k \left( \frac{F_0}{F_c} t_0 \right) ln \left( \frac{t + t_0}{t_0} \right)$', fontsize=14,color='black',rotation = 0)
    # ax.text(0.01, 0.032, r'$\tau = 0.07$', fontsize=12,color='black',rotation = 0, alpha = 0.8)
    # ax.text(0.01, 0.025, r'$Pr = 0.5$', fontsize=12,color='black',rotation = 0, alpha = 0.8)
    # ax.text(0.01, 0.02, r'$R = 10^{10}$', fontsize=12,color='black',rotation = 0, alpha = 0.8)
    # ax.text(0.01, 0.016, r'$F_0 / F_c = 10$', fontsize=12,color='black',rotation = 0, alpha = 0.8)


    plt.xscale("log")
    plt.yscale("log")
    ax.set_xlabel(r"Time", fontsize=12)
    ax.set_ylabel(r"$h(t)$", fontsize=12)
    plt.ylim(ymax = 1, ymin = 0.01)
    plt.xlim(xmin = 1e-4, xmax = 1e-1)
    ax.grid(which='major', alpha=0.85)
    ax.grid(which='minor', alpha=0.5)

    plt.legend()
    plt.grid()
    plt.show()
    if save_fig == True:
        fig2.savefig('deflux_scale_flux.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)






"hsin"
if mode == "hsin":
    t_i = 0.005
    time = np.logspace(-6, 1, 1000)
    time_short = np.logspace(np.log10(t_i), 1, 1000)

    
    test_int = test_5

    # fit_c = test_int.fit_for_h(function_to_fit = 1, outcome_domain = time_short, c_flux= True, t_i  = t_i)

    mpl.style.use('seaborn-v0_8')
    fig2, ax = plt.subplots(figsize=(8.5, 6))


    ax.vlines(x = test_1.t_0, ymin = 1e-6, ymax = 1, color = 'plum', linewidth = 0.9, alpha = 0.8)
    ax.text(test_1.t_0 * 1.01, 0.5, rf'$t_0$ ={test_1.t_0}', fontsize = 9, color = 'plum', alpha = 1, rotation = - 90)
    ax.vlines(x = test_2.t_0, ymin = 1e-6, ymax = 1, color = 'blue', linewidth = 0.9, alpha = 0.3)


    ax.plot(test_int.time, test_int.h_c, markersize = 4, color = 'darkgreen', label = rf"c")
    ax.plot(test_int.time, test_int.h_rho, "--", markersize = 4, color = 'darkgreen', label = rf"rho")
    ax.plot(test_int.time, test_int.h_cflux, "-s", markersize = 4, color = 'darkgreen', alpha = 0.8, label = rf"clux")


    # ax.plot(time_short, fit_c, alpha = 0.9, color = 'plum', linestyle = "-.")

    ax.tick_params(axis='y', direction = "in", right = True, length = 10, which = "major")
    ax.tick_params(axis='x',direction='in',top=True, length = 10, which = "major")
    ax.tick_params(axis='y', direction = "in", right = True, length = 4, which = "minor")
    ax.tick_params(axis='x',direction='in', length = 4, which = "minor")



    # ax.text(0.052, 0.47, f'k= {round(test_int.fit_param[0][0] / np.sqrt(test_int.F_ratio), 2)}', fontsize = 11, color = 'deeppink')
    # ax.text(0.052, 0.22, f'k= {round(test_2.fit_param[0][0], 2)}', fontsize = 11, color = 'darkblue')

    # ax.text(0.005, 0.045, r'Fit: $h^2(t) = k \left( \frac{F_0}{F_c} t_0 \right) ln \left( \frac{t + t_0}{t_0} \right)$', fontsize=14,color='black',rotation = 0)
    # ax.text(0.01, 0.032, r'$\tau = 0.07$', fontsize=12,color='black',rotation = 0, alpha = 0.8)
    # ax.text(0.01, 0.025, r'$Pr = 0.5$', fontsize=12,color='black',rotation = 0, alpha = 0.8)
    # ax.text(0.01, 0.02, r'$R = 10^{10}$', fontsize=12,color='black',rotation = 0, alpha = 0.8)
    # ax.text(0.01, 0.016, r'$F_0 / F_c = 10$', fontsize=12,color='black',rotation = 0, alpha = 0.8)


    plt.xscale("log")
    plt.yscale("log")
    ax.set_xlabel(r"Time", fontsize=12)
    ax.set_ylabel(r"$h(t)$", fontsize=12)
    plt.ylim(ymax = 1, ymin = 0.01)
    plt.xlim(xmin = 1e-4, xmax = 1e-1)
    ax.grid(which='major', alpha=0.85)
    ax.grid(which='minor', alpha=0.5)

    plt.legend()
    plt.grid()
    plt.show()
    if save_fig == True:
        fig2.savefig('test5_fitall.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)








# "veloh"
if mode == "velohf":

    h_dmy = np.linspace(0, 1, 100)
    h_nr = h_dmy**(1/3)
    h_r = h_dmy**(1/5)

    mpl.style.use('seaborn-v0_8')
    fig2, ax = plt.subplots(figsize=(8.5, 6))

    plot_target = "hsin"

    if plot_target == "h":
        ax.plot(h_dmy, h_nr, linestyle="--", color = 'deeppink', alpha = 0.4)
        ax.plot( test_1.h_rho, test_1.u_rms_all,  "-o", markersize = 4, color = 'deeppink', label = rf"deflux_rotation")
        ax.plot( test_2.h_rho , test_2.u_rms, "-s", markersize = 4, color = 'deeppink', label = rf"deflux_no rot")
        

        ax.plot(h_dmy, h_r, linestyle = "--", color = 'darkgrey', alpha = 0.4)

        ax.plot( test_3.h_rho[:-2] , test_3.u_rms[:-2], "-s", markersize = 4, color = 'darkgrey', label = rf"constant")
        ax.plot( test_4.h_rho, test_4.u_rms, "-o", markersize = 4, color = 'darkgrey', alpha = 0.7, label = rf"constant_rot_10")
        ax.plot( test_5.h_rho, test_5.u_rms, "-o", markersize = 4, color = 'black', label = rf"constant_rot_5")
        ax.set_xlabel(r"h-size", fontsize=12)


    if plot_target == "hsin":
        test_int = test_4

        ax.plot( test_int.time, test_int.u_rms_all,  "-o", markersize = 4, color = 'deeppink', label = rf"rms")
        ax.plot( test_int.time , test_int.u_rms_x, "-s", markersize = 4, color = 'darkgreen', label = rf"x")
        
        ax.plot( test_int.time, test_int.u_rms_y, "-o", markersize = 4, color = 'darkgrey', alpha = 0.7, label = rf"y")
        ax.plot( test_int.time, test_int.u_rms_z, "-o", markersize = 4, color = 'black', label = rf"z")
        ax.set_xlabel(r"time", fontsize=12)




    if plot_target == "f":
        ax.plot(h_dmy, h_nr, linestyle="--", color = 'deeppink', alpha = 0.4)
        ax.plot(test_1.tflux_avg, test_1.u_rms,  "-o", markersize = 4, color = 'deeppink', label = rf"deflux_rotation")
        ax.plot(test_2.tflux_avg , test_2.u_rms, "-s", markersize = 4, color = 'deeppink', label = rf"deflux_no rot")
        

        ax.plot(h_dmy, h_r, linestyle = "--", color = 'darkgrey', alpha = 0.4)

        ax.plot(  test_3.tflux_avg[:-2], test_3.u_rms[:-2], "-s", markersize = 4, color = 'darkgrey', label = rf"constant")
        ax.plot( test_4.tflux_avg , test_4.u_rms, "-o", markersize = 4, color = 'darkgrey', alpha = 0.7, label = rf"constant_rot_10")
        ax.plot( test_5.tflux_avg, test_5.u_rms, "-o", markersize = 4, color = 'black', label = rf"constant_rot_5")
        ax.set_xlabel(r"F_H", fontsize=12)



    if plot_target == "hf":
        ax.plot(h_dmy, h_nr, linestyle="--", color = 'deeppink', alpha = 0.4)
        ax.plot( test_1.h_rho * test_1.tflux_avg ** 2, test_1.u_rms,  "-o", markersize = 4, color = 'deeppink', label = rf"deflux_rotation")
        ax.plot( test_2.h_rho * test_2.tflux_avg , test_2.u_rms, "-s", markersize = 4, color = 'deeppink', label = rf"deflux_no rot")
        

        ax.plot(h_dmy, h_r, linestyle = "--", color = 'darkgrey', alpha = 0.4)

        ax.plot( test_3.h_rho[:-2] * test_3.tflux_avg[:-2], test_3.u_rms[:-2], "-s", markersize = 4, color = 'darkgrey', label = rf"constant")
        ax.plot( test_4.h_rho* test_4.tflux_avg ** 2, test_4.u_rms, "-o", markersize = 4, color = 'darkgrey', alpha = 0.7, label = rf"constant_rot_10")
        ax.plot( test_5.h_rho* test_5.tflux_avg ** 2, test_5.u_rms, "-o", markersize = 4, color = 'black', label = rf"constant_rot_5")
        ax.set_xlabel(r"hf", fontsize=12)

    
    if log_scale == True:
        plt.xscale("log")
        plt.yscale("log")

    ax.set_ylabel(r"$velo$", fontsize=12)
    # plt.ylim(ymax = 1e4, ymin = 1e2)
    # plt.xlim(xmin = 1e-2, xmax = 1)
    ax.grid(which='major', alpha=0.85)
    ax.grid(which='minor', alpha=0.5)


    plt.legend()
    plt.grid()
    plt.show()
    if save_fig ==  True:
        fig2.savefig('velo_time_test4.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)






"temp section"
if mode == "temp1":
    t_i = 0.005
    time = np.logspace(-6, 1, 1000)
    time_short = np.logspace(np.log10(t_i), 1, 1000)

    fit_2 = test_1.fit_for_h(function_to_fit = 54, outcome_domain = time_short, c_flux = True, t_i  = t_i)
    fit_3 = test_2.fit_for_h(function_to_fit = 53, outcome_domain = time_short, c_flux = True, t_i  = t_i)
    fit_1 = test_3.fit_for_h(function_to_fit = 51, outcome_domain = time_short, c_flux = True, t_i  = t_i)
    fit_5 = test_5.fit_for_h(function_to_fit = 2, outcome_domain = time_short, c = True, t_i  = t_i)


    mpl.style.use('seaborn-v0_8')
    fig2, ax = plt.subplots(figsize=(8.5, 6))


    ax.vlines(x = test_1.t_0, ymin = 1e-6, ymax = 1, color = 'plum', linewidth = 0.9, alpha = 0.8)
    ax.text(test_1.t_0 * 1.01, 0.5, rf'$t_0$ ={test_1.t_0}', fontsize = 9, color = 'plum', alpha = 1, rotation = - 90)
    ax.vlines(x = test_2.t_0, ymin = 1e-6, ymax = 1, color = 'blue', linewidth = 0.9, alpha = 0.3)


    ax.plot(time_short, fit_2, alpha = 0.6, color = 'deeppink', linestyle = "--")
    ax.plot(time_short, fit_1, alpha = 0.2, color = 'darkgreen', linestyle = "--")
    ax.plot(time_short, fit_3, alpha = 0.6, color = 'darkblue', linestyle = "--")


    ax.plot(test_1.time, test_1.h_cflux, "-o", markersize = 4, color = 'deeppink', label = rf"deflux_rotation")
    ax.plot(test_2.time, test_2.h_cflux, "-s", markersize = 4, color = 'darkblue', label = rf"deflux_no rot")
    ax.plot(test_3.time, test_3.h_cflux, "-s", markersize = 4, color = 'darkgreen', label = rf"constant")
    ax.plot(test_4.time, test_4.h_c, "-s", markersize = 4, color = 'darkgrey', label = rf"constant_rot_10")
    ax.plot(test_5.time, test_5.h_cflux, "-s", markersize = 4, color = 'black', label = rf"constant_rot_5")



    ax.tick_params(axis='y', direction = "in", right = True, length = 10, which = "major")
    ax.tick_params(axis='x',direction='in',top=True, length = 10, which = "major")
    ax.tick_params(axis='y', direction = "in", right = True, length = 4, which = "minor")
    ax.tick_params(axis='x',direction='in', length = 4, which = "minor")


    plt.xscale("log")
    plt.yscale("log")
    ax.set_xlabel(r"Time", fontsize=12)
    ax.set_ylabel(r"$h(t)$", fontsize=12)
    plt.ylim(ymax = 1, ymin = 0.01)
    plt.xlim(xmin = 1e-4, xmax = 1e-1)
    ax.grid(which='major', alpha=0.85)
    ax.grid(which='minor', alpha=0.5)

    plt.legend()
    plt.grid()
    plt.show()


    if save_fig == True:
        fig2.savefig('temp1.pdf', format="pdf",dpi=300, bbox_inches='tight', pad_inches = 0.05)




if mode == "plt01":
    plt.style.use('apj.mplstyle')

    fig, ax= plt.subplots(figsize=(9, 6))

    tot = np.array([test_1, test_2, test_3, test_4, test_5])
    tot_name = np.array(["deflux_rot", "deflux_nr", "con_nr", "Rotating Constant Flux", "con_rot_5"])


    chose = 3

    test_case = tot[chose]
    color = mpl.colormaps["cividis"](np.linspace(0,1, len(test_case.time)))
    color_2 = mpl.colormaps["RdPu"](np.linspace(0,1, len(test_case.time)))
    object_int = test_case.rho_data

    for i in range(len(test_case.time)):
        ax.plot(test_1.z_coord, object_int[i], color = color[i], alpha = 0.30 + i* (0.65/ (len(test_case.time) - 2)))

        ax.vlines(x =1 -  test_case.h_rho[i], ymax = 0.2* np.sqrt(test_case.h_rho[i]), ymin = 0, color = color[i])



    profiles_color =  color[len(test_case.time) // 2]
    flux_color = "black"
    # ax1 = ax.twinx()
    # ax1.plot(1- test_case.h_rho, test_case.tflux_avg, color = flux_color, label = r"<$\tilde{F_H}$> over h")
    # ax1.plot( [1-test_case.h_rho[20], 1-test_case.h_rho[20]] , [0, 0.1] , color = color_2[20], label = 'h')


    minor_legnth = 4
    major_length = 10
    axistitle_fontsize =14
    label_fontsize = 17
    
    

    tic_dir = 'in'
    ax.tick_params(axis='y', direction = tic_dir, right = False, length = major_length, which = "major")
    ax.tick_params(axis='y', direction = tic_dir, right = False, length = minor_legnth, which = "minor")
    # ax1.tick_params(axis='y', direction = tic_dir, right = True, length = major_length, which = "major", labelcolor = flux_color)
    # ax1.tick_params(axis='y', direction = tic_dir, right = True, length = minor_legnth, which = "minor", labelcolor = flux_color)

    ax.tick_params(axis='x',direction= tic_dir,top=False, length = major_length, which = "major")
    ax.tick_params(axis='x',direction= tic_dir, length = minor_legnth, which = "minor")



    ax.set_xlabel(r"$\tilde{z}$ ", fontsize = label_fontsize)

    ax.set_ylabel(r"$\rho/ \rho_0$", fontsize = label_fontsize)
    # ax1.set_ylabel(r"$\tilde{F_H} / (F_0 / F_c) $", color = flux_color, fontsize = label_fontsize)

    ax.set_title(tot_name[chose])


    
    # ins_ax = fig.add_axes([0.14, 0.3, .13, .25])  # [x, y, width, height] w.r.t. fig

    # pick_ins =  20
    # x = np.linspace(1-test_case.h_rho[pick_ins], 1, 20)
    # y1 = np.ones(20) 
    # y2 = np.zeros(20)

    # ins_ax.fill_between(x, y1, y2, where=(y1 > y2), color=color_2[pick_ins], alpha=0.7, label='y1 > y2')


    # ins_ax.plot(test_case.z_coord, test_case.t_flux_data[pick_ins], color="black")
    # ins_ax.tick_params(axis='y', right = True, left = False, labelleft = False, labelright = True )


    # cmap = plt.get_cmap("viridis")
    # cmap22 = plt.get_cmap("RdPu")
    # norm = plt.Normalize(vmin= test_case.time[0], vmax=test_case.time[-1])
    # sm1 = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    # sm2 = plt.cm.ScalarMappable(norm=norm, cmap=cmap22)


    # cbar1 = fig.colorbar(sm1, ax=ax, orientation='vertical', shrink = 0.8, pad = -0.11)
    # cbar1.set_label(r'Time($\tilde{t}$)', labelpad=0.05)

    # cbar2 = fig.colorbar(sm2, ax=ax, orientation='vertical', shrink = 0.8, pad = 0.09)

    # # cbar1.ax.set_yscale('symlog', linthresh=1e-6)

    # cbar2.ax.yaxis.set_ticks([])
    # # cbar2.set_label(right = False)
    # cbar2.ax.yaxis.label.set_visible(False)


    ax.text(0.04, 0.03, f'Convection \n Boundary \n Location', fontsize = 12, color = color_2[-1])
    ax.text(1 - test_case.h_rho[len(test_case.h_rho)//2], 0.22, r'$\leftarrow$ Time', fontsize = 18)


    # ax.text(0.5, 1.1, r'Average $\tilde{F_H}$ over h', fontsize = 12, color = "black")



    # # Create connection from main figure to inset
    # con = ConnectionPatch(
    #     xyA=(1 - test_case.h_rho[pick_ins], test_case.tflux_avg[pick_ins]),  # Point in main axes
    #     xyB=(1, 1),  # Same point in inset axes
    #     coordsA="data",          # Data coordinates for A
    #     coordsB="data",          # Data coordinates for B
    #     axesA=ax1,                # Main axes
    #     axesB=ins_ax,          # Inset axes
    #     color="red",             # Line color
    #     linestyle="dashed",      # Line style
    #     linewidth=1              # Line width
    # )

    # # Add the connection to the figure (not to either axes)
    # fig.add_artist(con)



    # ax1.legend(loc = [0.05, 0.3])
    plt.show()



    if save_fig == True:
        fig.savefig('poster_plot01.png', format="png",dpi=300, bbox_inches='tight', pad_inches = 0.05)




if mode == "plt03":
    plt.style.use('apj.mplstyle')

    fig, ax= plt.subplots(figsize=(8, 8))


    t_i = 0.005
    time = np.logspace(-6, 1, 1000)
    time_short = np.logspace(np.log10(t_i), 1, 1000)

    color = mpl.colormaps["viridis"](np.linspace(0,1, len(test_3.time)))
    color_2 = mpl.colormaps["RdPu"](np.linspace(0,1, len(test_3.time)))
    
    non_rot_color = color[0]
    Rot_color = color_2[-15]


    fit_1 = test_1.fit_for_h(function_to_fit = 54, outcome_domain = time_short, rho = True, t_i  = t_i)
    fit_2 = test_2.fit_for_h(function_to_fit = 53, outcome_domain = time_short, c_flux = True, t_i  = t_i)

    fit_3 = test_3.fit_for_h(function_to_fit = 1, outcome_domain = time_short, rho = True, t_i  = t_i)
    fit_4 = test_4.fit_for_h(function_to_fit = 52, outcome_domain = time_short, rho = True, t_i  = t_i)



    ax.vlines(x = test_1.t_0, ymin = 1e-6, ymax = 1, color = color[20], linewidth = 1, alpha = 0.6)
    ax.text(test_1.t_0 * 1.01, 0.5, rf'$t_0$ ={test_1.t_0}', fontsize = 13, color = color[20], alpha = 0.6, rotation = - 90)
    ax.vlines(x = test_2.t_0, ymin = 1e-6, ymax = 1, color = 'blue', linewidth = 0.9, alpha = 0.3)

    ax.plot(time_short, fit_2 * 0.8, alpha = 0.6, color = non_rot_color, linestyle = "--", label = 'Prediction(Non-Rotating)')
    ax.plot(time_short, fit_1, alpha = 0.6, color = Rot_color, linestyle = "--",  label = 'Prediction(Rotating)')

    ax.plot(time_short, fit_3, alpha = 0.6, color = non_rot_color, linestyle = "--")
    ax.plot(time_short, fit_4 * 0.83, alpha = 0.6, color = Rot_color, linestyle = "--")

    ax.plot(test_3.time[:-3], test_3.h_rho[:-3], "-o", markersize = 4, color = non_rot_color, label = rf"Non-rotating, constant flux")
    ax.plot(test_2.time, test_2.h_rho, "-s", markersize = 4, color = non_rot_color, label = rf"Non-rotating, decreasing flux")
    
    ax.plot(test_4.time, test_4.h_rho, "-o", markersize = 4, color = Rot_color, label = rf"Rotating, constant flux")
    ax.plot(test_1.time, test_1.h_rho, "-s", markersize = 4, color = Rot_color, label = rf"Rotating, decreasing flux")


    # ax.plot(test_5.time, test_5.h_rho, "-s", markersize = 4, color = 'black', label = rf"constant_rot_5")


    ax.tick_params(axis='y', direction = "in", right = True, length = 10, which = "major")
    ax.tick_params(axis='x',direction='in',top=True, length = 10, which = "major")
    ax.tick_params(axis='y', direction = "in", right = True, length = 4, which = "minor")
    ax.tick_params(axis='x',direction='in', length = 4, which = "minor")

    plt.xscale("log")
    plt.yscale("log")
    ax.set_xlabel(r"Time ($\tilde{t}$)", fontsize=16)
    ax.set_ylabel(r"Convection Zone Size$(h)$", fontsize=16)

    ax.text(0.013, 0.75, f'Reaching \n Fully Mixed \n State', fontsize = 14, color = non_rot_color)
    ax.text(0.05, 0.25, f'Reaching \n Saturation', fontsize = 14, color = Rot_color)

    plt.ylim(ymax = 1, ymin = 0.09)
    plt.xlim(xmin = 2e-3, xmax = 1e-1)
    ax.grid(which='major', alpha=0.85)
    ax.grid(which='minor', alpha=0.5)


    plt.legend(loc = 'lower right')
    plt.grid()
    plt.show()
    if save_fig == True:
        fig.savefig('poster_plot03.png', format="png",dpi=300, bbox_inches='tight', pad_inches = 0.05)





if mode == "plt02":

    fig, ax= plt.subplots(figsize=(7.5, 5.5))


    color = mpl.colormaps["viridis"](np.linspace(0,1, len(test_3.time)))
    color_2 = mpl.colormaps["RdPu"](np.linspace(0,1, len(test_3.time)))

    fh_dummy =  np.logspace(-5, 1, 2000)
    u_scale_line_1 =  test_1.u_rms_all[-10] * fh_dummy ** (1/3)
    u_scale_line_2 =  test_1.u_rms_all[-10] * fh_dummy ** (1/5) 


    non_rot_color = color[0]
    Rot_color = color_2[-15]


    ax.plot( fh_dummy, u_scale_line_1 * 1.5, "--", color = "darkgrey", alpha = 0.9)
    ax.plot( fh_dummy, u_scale_line_2, color = "darkgrey", alpha = 0.9)

    ax.plot(test_3.h_rho[:-3] * test_3.tflux_avg[:-3] , test_3.u_rms_z[:-3] ,'.', color = non_rot_color, label = "Non-rotating, constant flux")
    ax.plot(test_2.h_rho * test_2.tflux_avg, test_2.u_rms_z,'2', color = non_rot_color, label = "Non-rotating, decreasing flux")
    
    
    ax.plot(test_4.h_rho * test_4.tflux_avg **2, test_4.u_rms_z,'.', color = Rot_color, label = "Rotating, constant flux")
    ax.plot(test_1.h_rho * test_1.tflux_avg **2, test_1.u_rms_z, "2", color= Rot_color,  label =  rf"Rotating, decreasing flux")



    ax.text(0.08, 750, r"$U_z \sim (hF)^{1/3}$", color = "darkgrey", fontsize = 16)
    ax.text(0.2, 370, r"$U_z \sim (hF^2)^{1/5}$", color = "darkgrey", fontsize = 16)

    ax.set_xlabel(r"f (h, $F_H$)", fontsize=16)
    ax.set_ylabel(r"<$U_z$>", fontsize=16)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(ymax = 1e3, ymin =0.9e2)
    ax.set_xlim(xmax = 1, xmin =1e-3)


    ax.legend(loc = 'lower right')



    ins_ax = fig.add_axes([0.16, 0.6, .22, .26])  # [x, y, width, height] w.r.t. fig
    
    test_case =  test_4
    pick_ins =  10
    x = np.linspace(1-test_case.h_rho[pick_ins], 1, 20)
    y1 = np.ones(20) * np.max(test_case.u_z_data[pick_ins] ** 0.5)
    y2 = np.zeros(20)

    ins_ax.fill_between(x, y1, y2, where=(y1 > y2), color=color_2[pick_ins], alpha=0.7, label='y1 > y2')


    ins_ax.plot(test_case.z_coord, test_case.u_z_data[pick_ins] ** 0.5, linewidth = 1.2, color="black", alpha = 0.8)
    ins_ax.tick_params(axis='y', right = True, left = False, labelleft = False, labelright = True )
    ins_ax.text(0.62, 10, r"<$U_z$>", fontsize = 13, color = "black")


    # Create connection from main figure to inset
    con = ConnectionPatch(
        xyA=(test_case.h_rho[pick_ins] * test_case.tflux_avg[pick_ins] **2, test_case.u_rms_z[pick_ins]),  # Point in main axes
        xyB=(1, 0),  # Same point in inset axes
        coordsA="data",          # Data coordinates for A
        coordsB="data",          # Data coordinates for B
        axesA=ax,                # Main axes
        axesB=ins_ax,          # Inset axes
 color="red",             # Line color
        linestyle="dashed",      # Line style
        linewidth=1              # Line width
    )

    # Add the connection to the figure (not to either axes)
    fig.add_artist(con)




    plt.grid()
    plt.show()
    if save_fig ==  True:
        fig.savefig('poster_plot02.png', format="png",dpi=300, bbox_inches='tight', pad_inches = 0.05)





if mode == "epi":
    plt.plot(test_1.time, test_1.epsilon_measured, label =  "dwr")
    plt.plot(test_2.time, test_2.epsilon_measured, label =  "dnr")
    plt.plot(test_3.time, test_3.epsilon_measured, label =  "cnr")
    plt.plot(test_4.time, test_4.epsilon_measured, label =  "cwr10")
    plt.plot(test_5.time, test_5.epsilon_measured, label =  "cwr5")


    plt.xlabel("time")
    plt.ylabel(r"$\epsilon$")
    plt.legend()
    plt.show()


if mode == "gamma":
    plt.plot(test_1.time, test_1.gamma_measured, label =  "dwr")
    plt.plot(test_2.time, test_2.gamma_measured, label =  "dnr")
    plt.plot(test_3.time[:-3], test_3.gamma_measured[:-3], label =  "cnr")
    plt.plot(test_4.time, test_4.gamma_measured, label =  "cwr10")
    plt.plot(test_5.time, test_5.gamma_measured, label =  "cwr5")

    # plt.xscale("log")
    # plt.yscale("log")


    plt.xlabel("time")
    plt.ylabel(r"$\gamma$")
    plt.legend()
    plt.show()



if mode == "hcom":
    test_case = test_3
    plt.plot(test_case.time, test_case.h_c, label = "c")
    plt.plot(test_case.time, test_case.h_cflux, label = "cflux")
    plt.plot(test_case.time, test_case.h_rho, label = "rho")
    plt.plot(test_case.time, test_case.h_avg, label = "avg")


    plt.xlabel("time")
    plt.ylabel("h")
    plt.legend()
    plt.show()