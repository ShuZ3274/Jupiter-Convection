"""
Dedalus script simulating 2D foir fully compressible convections
with compositionally startified, isothermal initial condition
cooled by top boundary of the domain with constant heat flux
An extension to: https://iopscience.iop.org/article/10.3847/2041-8213/ae0459
"""



# calculation specific dependency
import numpy as np
import dedalus.public as d3
from   dedalus.extras import flow_tools
from   dedalus.tools  import post
# ststus reporting
import logging
logger = logging.getLogger(__name__)
# Multi-processing for cluster usage
from   mpi4py import MPI
from   collections import OrderedDict
# Visulization
import matplotlib.pyplot as plt






# ----------------------------------- Problem parameters -------------------------------------------

# independent fluid properties
Ra       = 2e8   # g*L^3 / kappa_T*nu
Pr       = 0.1   # Prandtl number, nu / kappa_T
gamma    = 5/3   # monatomic ideal gas
n_rho    = 4     
Le       = 10    # Lewis number, kappa_T / kappa_X


# independent simulation parameters
aspect   = 4     # ratio between vertical and horizontal direction
Nz       = 1024  # resolution
dealias  = 3/2


# dependent parameters
ln16     = np.log(16) 


# box dimensions
Lz       = 1
Nx, Lx   = Nz*2, Lz * aspect
Ny, Ly   = Nz*2, Lz * aspect
 





# ----------------------------------- Dedalus domain setup -------------------------------------------

dtype        = np.float64
coords       = d3.CartesianCoordinates('x', 'z')
dist         = d3.Distributor(coords, dtype=dtype)
xbasis       = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)   # periodic
zbasis       = d3.ChebyshevT(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)    # bounded
bases        = (xbasis, zbasis)
h_bases      = (xbasis)  # horizontal bases collection
x, z         = dist.local_grids(xbasis, zbasis)
ex, ez       = coords.unit_vector_fields(dist)

# variables
u            = dist.VectorField(coords, name='u', bases=bases)
T1           = dist.Field(name='T', bases=bases)
ln_rho1      = dist.Field(name='ln_rho1', bases=bases)
x1           = dist.Field(name='x1', bases=bases)

# Taus
tau_T1       = dist.Field(name='tau_T1', bases=h_bases)
tau_T2       = dist.Field(name='tau_T2', bases=h_bases)
tau_u1       = dist.VectorField(coords, name='tau_u1', bases=h_bases)
tau_u2       = dist.VectorField(coords, name='tau_u2', bases=h_bases)
tau_x1       = dist.Field(name='tau_x1', bases=h_bases)
tau_x2       = dist.Field(name='tau_x2', bases=h_bases)
lift         = lambda A: d3.Lift(A, zbasis.derivative_basis(1), -1)

# Background fields (only depend on z)
T0           = dist.Field(name='T0', bases=(zbasis,))
ln_rho0      = dist.Field(name='ln_rho0', bases=(zbasis))
x0           = dist.Field(name='x0', bases = zbasis)
grad_T0      = dist.VectorField(coords, name='grad_T0', bases=(zbasis,))
grad_ln_T0   = dist.VectorField(coords, name='grad_ln_T0', bases=(zbasis,))
grad_ln_rho0 = dist.VectorField(coords, name='grad_ln_rho0', bases=(zbasis,))
grad_x0      = dist.VectorField(coords, name='grad_x0', bases=(zbasis,))

nu           = dist.Field(name='nu', bases=(zbasis,))
grad_nu      = dist.VectorField(coords, name='grad_nu', bases=(zbasis,))
grad_kappa   = dist.VectorField(coords, name='grad_kappa', bases=(zbasis,))
kappa_T0     = dist.Field(name='kappa', bases=(zbasis,))
kappa_X0     = dist.Field(name='kappa', bases=(zbasis,))






# ----------------------------------- Initial conditions -------------------------------------------
# subscript 0 denotes initial background, 1 denotes the perturbation to the background


# linear composition isothermal
L_o_H    = 15/16*(n_rho/ln16 - 1)
T_func       = lambda z: 1
grad_T_func  = lambda z: 0
x0_func      = lambda z: 1 - z
grad_x0_func = lambda z: -1
rho_func     = lambda z: (16/(1 + 15*z))**(n_rho/ln16)
ln_rho_func  = lambda z: n_rho/ln16*np.log(16/(1 + 15*z))
heat_flux    = (15/16*(2-1/gamma) + (1-1/gamma)*L_o_H) * (5/2)



T0['g']         = T_func(z)
grad_T0['g'][0] = 0
grad_T0['g'][1] = grad_T_func(z)
T1['g']         = 0

x0['g']         = x0_func(z)
grad_x0['g'][0] = 0
grad_x0['g'][1] = grad_x0_func(z)
x1['g']         = 0

ln_rho0['g']   = np.log(rho_func(z))
grad_ln_rho0   = d3.grad(ln_rho0)
ln_rho1['g']   = 0

nu['g']        = 1 / rho_func(z)
kappa_T0['g']  = 1 / rho_func(z)
kappa_X0['g']  = 1 / rho_func(z)
grad_nu        = d3.grad(nu)
grad_kappa_T   = d3.grad(kappa_T0)







# ----------------------------------- Equations -------------------------------------------

# Substitutions
# Taus
grad_u        = d3.grad(u) - ez*lift(tau_u1) # First-order reduction
grad_T1       = d3.grad(T1) - ez*lift(tau_T1)
grad_x1       = d3.grad(x1) - ez*lift(tau_x1)
Lap_T1        = d3.div(grad_T1)
grad_ln_rho1  = d3.grad(ln_rho1)
div_u         = d3.trace(grad_u)

#viscosity
eye = dist.TensorField((coords, coords), name='I',bases=(zbasis))
eye['g']      = 0
eye['g'][0,0] = 1
eye['g'][1,1] = 1

#viscosity
E             = 0.5*(grad_u + d3.trans(grad_u))
sigma         = 2*(E - (1/3)*div_u*eye)
div_sigma     = d3.div(sigma)
VH            = 2 * nu * (d3.trace(d3.dot(E,E)) - (1/3)*div_u*div_u)
viscous_diffusion_L = nu*(div_sigma + sigma@grad_ln_rho0) + sigma@grad_nu
viscous_diffusion_R = nu*sigma@grad_ln_rho1

# Thermal diffusion
mu_inv_term         = 1 - 15/16*(x0 + x1)
thermal_diffusion_L =  kappa_T0*(Lap_T1 + d3.dot(grad_T1, grad_ln_rho0) + d3.dot(grad_T0, grad_ln_rho1)) + d3.dot(grad_kappa,grad_T1)
# thermal_diffusion_R =  kappa_T0*(grad_T1 @ grad_ln_rho1 + grad_T0 @ grad_ln_rho0 + d3.div(grad_T0) + d3.grad(np.log(mu_inv_term)) @  (grad_T0+grad_T1) ) + grad_kappa @ grad_T0 
thermal_diffusion_R =  kappa_T0*(grad_T1 @ grad_ln_rho1 + d3.dot(d3.grad(mu_inv_term)/mu_inv_term , grad_T1))


#Coefficients after nondimensionalization
# dimensionless heat capacity
cV         = 1/(gamma - 1) * (1 - 15/16*(x0 + x1))
# Buoyancy
C_B        = Ra * Pr / L_o_H
# viscous diffusion
C_V_diff   = Pr
C_V_Diss   = L_o_H / Ra / cV



# momentum terms
grad_p1_L = grad_T1 - 15/16*( x0*grad_T1 + x1*grad_T0 + T0*grad_x1 + T1*grad_x0 + T0*x0*grad_ln_rho1 + T0*x1*grad_ln_rho0 + T1*x0*grad_ln_rho0) \
            + T0*grad_ln_rho1 + T1*grad_ln_rho0
grad_p1_R = T1*grad_ln_rho1 - 15/16*( x1*grad_T1 + T1*grad_x1 + T0*x1*grad_ln_rho1 + T1*x0*grad_ln_rho1 +T1*x1*grad_ln_rho0 + T1*x1*grad_ln_rho1 )
# chemical diffusion:
dx2dt_L   = kappa_X0 * (grad_ln_rho0@grad_x1 + grad_ln_rho1@grad_x0 + d3.div(grad_x1)) + d3.grad(kappa_X0)@grad_x1
dx2dt_R   = kappa_X0 * (grad_ln_rho1@grad_x1 + grad_ln_rho0@grad_x0 + d3.div(grad_x0)) + d3.grad(kappa_X0)@grad_x0


# define the problem
variables = [ln_rho1, T1, u, x1, tau_u1, tau_u2, tau_T1, tau_T2, tau_x1, tau_x2]
problem = d3.IVP(variables, namespace=locals())


# mass continuity
problem.add_equation("dt(ln_rho1) + u@grad_ln_rho0 + div_u = - u@grad_ln_rho1") 
# momentum equation
problem.add_equation("dt(u) + C_B*grad_p1_L - C_V_diff*viscous_diffusion_L + lift(tau_u2) =\
                    - u@grad(u) - C_B*grad_p1_R + C_V_diff*viscous_diffusion_R") 
# energy equation
problem.add_equation("dt(T1) + u@grad_T0 + (gamma-1)*T0*div_u + lift(tau_T2) - gamma*thermal_diffusion_L = \
                    -u@grad_T1 - (gamma-1)*T1*div_u + gamma*thermal_diffusion_R+ C_V_Diss*VH")
# chemical transposrt
problem.add_equation("dt(x1) + u@grad_x0 - 1/Le*dx2dt_L  + lift(tau_x2) = - u@grad_x1 + 1/Le*dx2dt_R")





# Boundary Conditions
heat_flux  = (15/16*(2-1/gamma) + (1-1/gamma)*L_o_H) * (5/2)  # rho_0 * kappa_t0 * cp = 5/2
rho_mu_inv = np.exp(ln_rho0 + ln_rho1)*(1 - 15/16*(x0 + x1))
print("thershold heatflux:", heat_flux)
problem.add_equation("ez@(grad(T1)(z=Lz)) = -2* heat_flux / (5/2*rho_mu_inv(z=Lz))")
problem.add_equation("ez@(grad(T1)(z=0))  = 0")


# conserved mass fraction (forced 0 on perturbation)
problem.add_equation("ez@(grad_x1(z=0)) = -grad_x0_func(0)")
problem.add_equation("ez@(grad_x1(z=1)) = -grad_x0_func(1)")

# Inpenetrable
problem.add_equation("ez@(u(z=1)) = 0")
problem.add_equation("ez@(u(z=0)) = 0")
# Stress free
strain_rate = d3.grad(u) + d3.trans(d3.grad(u))
shear_stress_top = ez@(strain_rate(z=Lz))
shear_stress_bot = ez@(strain_rate(z=0))
problem.add_equation("ex@shear_stress_top = 0")
problem.add_equation("ex@shear_stress_bot = 0")










# ----------------------------------- Solver construction -------------------------------------------

# solver property
timestepper   = d3.RK443
cfl_safety    = 0.1
max_timestep  = 1e-3
max_sim_time  = 10
timestep      = 1e-6 * max_timestep # first time step


solver = problem.build_solver(timestepper, ncc_cutoff=1e-8)
solver.stop_sim_time = max_sim_time

# CFL
CFL = d3.CFL(solver, 
             initial_dt=timestep, 
             cadence   =1, 
             safety    =cfl_safety, 
             threshold =0.1,
             max_change=2, 
             min_change=0.1, 
             max_dt    =max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property(np.sqrt(u@u)/Pr, name='Re')
flow.add_property(np.sqrt(L_o_H / Ra / Pr *u@u/ (1 - 15/16*(x0 + x1))), name='Mach')




# ----------------------------------- Data storage -------------------------------------------
# saved files are hp5 binary files for optimization of space allocation
# snapshots contain full data (2D or 3D), saved on a less frequent cadence
# profiles are horizontally averaged data (only contains z-direction), saved more frequently

# check saved with wall time
checkpoint = solver.evaluator.add_file_handler('checkpoint', 
                                               wall_dt=3600, 
                                               max_writes=1, 
                                               parallel='gather')
checkpoint.add_tasks(solver.state, layout='g')



# snapshots will save all raw variables
snapshots = solver.evaluator.add_file_handler('snapshots', 
                                              sim_dt    =0.2*max_timestep, 
                                              max_writes=100)
snapshots.add_task(ln_rho1, name='r1')
snapshots.add_task(T1, name='T1')
snapshots.add_task(x1, name='x1')
snapshots.add_task(u@ex, name='ux')
snapshots.add_task(u@ez, name='uz')



# profiles contain plane average values
profiles = solver.evaluator.add_file_handler('profiles', 
                                             sim_dt    =0.01*max_timestep, 
                                             max_writes=100)

plane_avg  = lambda A: d3.Integrate(A, coords['x'])/(Lx)
# volume_avg = lambda A: d3.Integrate(A, (coords['x'], coords['z']))/(aspect*Lz*Lz)

x2    = x0 + x1
x2    = x0 + x1
Tfull = T0 + T1
lnr2  = ln_rho1 + ln_rho0
mu_inv = 1 - 15/16*x2
ux_mag = ex@u
uz_mag = ez@u
dT2    = ez@d3.grad(Tfull)
dx2    = ez@d3.grad(x2)
rhou2  = np.exp(lnr2) * u@u
rhoTmu = np.exp(lnr2) * Tfull * mu_inv
rhomu  = np.exp(lnr2) * mu_inv
rhomudT= np.exp(lnr2) * mu_inv*dT2

profiles.add_task(plane_avg(rhou2), name = "rhou2")
profiles.add_task(plane_avg(rhoTmu), name = "rhoTmu")
profiles.add_task(plane_avg(rhomu), name = "rhomu")
profiles.add_task(plane_avg(x2), name = "x2")
profiles.add_task(plane_avg(lnr2), name = "lnr2")
profiles.add_task(plane_avg(Tfull), name = "T2")
profiles.add_task(plane_avg(dT2), name = "dT2")
profiles.add_task(plane_avg(ux_mag), name = "ux")
profiles.add_task(plane_avg(uz_mag), name = "uz")
profiles.add_task(plane_avg(dx2), name = "dx2")
profiles.add_task(plane_avg(rhomudT), name = "rhomudT")


# charateristic values
Ma2  = L_o_H / Ra / Pr * (u@u / (1 - 15/16*(x0 + x1)))       # Mach number squared
profiles.add_task((plane_avg(Ma2))**0.5, name = "Mach")











# ----------------------------------- Main Loop -------------------------------------------

restart = False
if not restart:
    T1.fill_random('g', seed=42, distribution='normal', scale=1e-6) # Random noise
    T1.low_pass_filter(scales=0.25)
    file_handler_mode = 'overwrite'
else:
    write, timestep = solver.load_state('checkpoints/checkpoints_s20.h5')
    file_handler_mode = 'append'


# Main loop
startup_iter = 10
try:
    logger.info('Starting main loop')
    while solver.proceed:

        timestep = CFL.compute_timestep()

        # break loop if CFL give too small timestep
        if timestep < 1e-12:
            break
      
        solver.step(timestep)


        if (solver.iteration-1) % 10 == 0:
            avg_Re = flow.volume_integral('Re') / (Lx * Lz)
            avg_Ma = flow.volume_integral('Mach') / (Lx * Lz)
            logger.info("Iter=%i, Time=%e, dt=%e, avg(Re)=%f, avg(Ma)=%f"
                        %(solver.iteration, solver.sim_time, timestep, avg_Re, avg_Ma))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()