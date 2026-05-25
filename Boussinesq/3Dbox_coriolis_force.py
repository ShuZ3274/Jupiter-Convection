"""
Dedalus script simulating 3D rotationally-constrainted, compositionally-stabilized, entraining convection.
For more details, please refer to: https://iopscience.iop.org/article/10.3847/2041-8213/ae0459
"""

import numpy as np
import dedalus.public as d3
import logging
from mpi4py import MPI
logger = logging.getLogger(__name__)



# Parameters - load in from parameter file (Nx, Nz, Ra, Pr, etc)
from control_parameters import parameters
locals().update(parameters)
# Additional Parameters
timestepper = d3.SBDF2
# timestepper = d3.RK443

cfl_safety = 0.2
max_timestep = 1e-5
dtype = np.float64
ncpu = MPI.COMM_WORLD.size
log2 = np.log2(ncpu)
if log2 == int(log2):
    mesh = [int(2**np.ceil(log2/2)),int(2**np.floor(log2/2))]
logger.info("running on processor mesh={}".format(mesh))




# Bases
coords = d3.CartesianCoordinates('x','y', 'z')
dist = d3.Distributor(coords, dtype=dtype,mesh=mesh)
xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, aspect*Lz), dealias=dealias)
ybasis = d3.RealFourier(coords['y'], size=Ny, bounds=(0, aspect*Lz), dealias=dealias)
zbasis = d3.ChebyshevT(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)


# Fields
bases = (xbasis, ybasis, zbasis)
horiz_bases = (xbasis, ybasis)

p = dist.Field(name='p', bases=bases)
T = dist.Field(name='T', bases=bases)
C = dist.Field(name='C', bases=bases)
u = dist.VectorField(coords, name='u', bases=bases)

tau_p = dist.Field(name='tau_p')
tau_T1 = dist.Field(name='tau_T1', bases=horiz_bases)
tau_T2 = dist.Field(name='tau_T2', bases=horiz_bases)
tau_C1 = dist.Field(name='tau_C1', bases=horiz_bases)
tau_C2 = dist.Field(name='tau_C2', bases=horiz_bases)
tau_u1 = dist.VectorField(coords, name='tau_u1', bases=horiz_bases)
tau_u2 = dist.VectorField(coords, name='tau_u2', bases=horiz_bases)


x, y, z= dist.local_grids(xbasis, ybasis, zbasis)
x_de, y_de, z_de = dist.local_grids(xbasis, ybasis, zbasis, scales=(dealias, dealias, dealias))
ex, ey, ez = coords.unit_vector_fields(dist)
lift_basis = zbasis.derivative_basis(1)
lift = lambda A: d3.Lift(A, lift_basis, -1)

grad_u = d3.grad(u) + ez*lift(tau_u1) # First-order reduction
grad_T = d3.grad(T) + ez*lift(tau_T1) # First-order reduction
grad_C = d3.grad(C) + ez*lift(tau_C1) # First-order reduction


# operators
grad = d3.grad
dot = d3.DotProduct
vorticity = d3.curl(u)

FK_vert = dot(ez,u)**2
FK_hor_x = dot(ex,u)**2
FK_hor_y = dot(ey,u)**2


#for stress-free BCs
strain_rate = d3.grad(u) + d3.trans(d3.grad(u))
shear_stress_top = ez@(strain_rate(z=Lz))
shear_stress_bot = ez@(strain_rate(z=0))


#Initial / Background profiles
T['g'] = 1
C['g'] = 1 - z


# Problem
problem = d3.IVP([p, T, C, u, tau_p, tau_T1, tau_T2, tau_C1, tau_C2, tau_u1, tau_u2], namespace=locals())
problem.add_equation("trace(grad_u) + tau_p = 0")
problem.add_equation("dt(T) - div(grad_T) + lift(tau_T2) = - (u@grad(T))")
problem.add_equation("dt(C) - tau*div(grad_C)  + lift(tau_C2) = - (u@grad(C))")
# problem.add_equation("dt(u) - Prandtl*div(grad_u) + grad(p) - Rayleigh*Prandtl*(T - C)*ez  + lift(tau_u2) = - cross(vorticity,u)")


# with rotation
problem.add_equation("dt(u) - Prandtl*div(grad_u) + grad(p) - Rayleigh*Prandtl*(T - C)*ez + lift(tau_u2) + omega * cross(ez, u) = - cross(vorticity,u)")



# function for decreasing flux 
# typical convection start time
t_0 = t0
temp = dist.Field()
t = dist.Field()
# initial time
t["g"]  = max_timestep * 0.01
def flux_function(*args):
    t = args[0].data
    # function that contorl boundary condition
    temp["g"] = t_0 /(t_0 + t)
    return temp["g"]

def flux_t(*args, domain=temp.domain, F= flux_function):
    return d3.GeneralFunction(dist, domain, layout='g', tensorsig=(), dtype=np.float64, func=F, args=args)


# Boundary Conditions on T and C
problem.add_equation("ez@(grad(C)(z=0)) = 0")
problem.add_equation("ez@(grad(C)(z=Lz)) = 0")
problem.add_equation("ez@(grad(T)(z=Lz)) = - F_ratio * flux_t(t)")
problem.add_equation("ez@(grad(T)(z=0)) = 0")



# Inpenetrable
problem.add_equation("ez@(u(z=Lz)) = 0")
problem.add_equation("ez@(u(z=0)) = 0")



# Stress free
problem.add_equation("ex@shear_stress_top = 0")
problem.add_equation("ex@shear_stress_bot = 0")

problem.add_equation("ey@shear_stress_top = 0")
problem.add_equation("ey@shear_stress_bot = 0")


# Pressure gauge
problem.add_equation("integ(p) = 0")

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time



# To restart:
#  (1) change use_checkpoint to True
#  (2) change name of checkpoint below to the one you wish to use
#use_checkpoint = True
use_checkpoint = False
timestep = max_timestep * 0.01
if use_checkpoint:
    write, timestep = solver.load_state('./checkpoint/checkpoint_s2.h5', -1)
else:
    # Initial conditions
    noise = dist.Field(name='noise', bases=bases)
    noise.fill_random('g', seed=42, distribution='normal', scale=1e-5) # Random noise
    noise.low_pass_filter(scales=0.25)

    noise.change_scales(dealias)
    T.change_scales(dealias)
    T['g'] += noise['g']


# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=1e-3, max_writes=10)
snapshots.add_task(T, name='T')
snapshots.add_task(C, name='C')


plane_avg = lambda A: d3.Integrate(A, (coords['x'], coords['y']))/((aspect*Lz)**2)
profiles = solver.evaluator.add_file_handler('profiles', sim_dt=1e-5, max_writes=100)
profiles.add_task(plane_avg(T), name='T')
profiles.add_task(plane_avg(dot(ez, u*T)), name='T_conv_flux')
profiles.add_task(plane_avg(dot(ez, grad(T))), name='T_grad')

profiles.add_task(plane_avg(C), name='C')
profiles.add_task(plane_avg(dot(ez, u*C)), name='C_conv_flux')
profiles.add_task(plane_avg(dot(ez, grad(C))), name='C_grad')

profiles.add_task(plane_avg(dot(ez, 0.5*u*dot(u,u))), name='KE_flux')
profiles.add_task(plane_avg(FK_vert), name='KE_vert')
profiles.add_task(plane_avg(FK_hor_x), name='KE_hor_x')
profiles.add_task(plane_avg(FK_hor_y), name='KE_hor_y')


# Checkpoint
checkpoint = solver.evaluator.add_file_handler('checkpoint', wall_dt=3600, max_writes=1, parallel='gather')
checkpoint.add_tasks(solver.state, layout='g')


# CFL
CFL = d3.CFL(solver, initial_dt=timestep, cadence=1, safety=cfl_safety, threshold=0.05,
             max_dt=max_timestep)
CFL.add_velocity(u)

# Flow properties
flow = d3.GlobalFlowProperty(solver, cadence=10)
flow.add_property(np.sqrt(u@u), name='Pe')

# Main loop
startup_iter = 10
try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()

      
        if timestep < 1e-11:
            break

        t["g"] = solver.sim_time          
        solver.step(timestep)

        
        if (solver.iteration-1) % 10 == 0:
            max_Pe = flow.max('Pe')
            logger.info('Iteration=%i, Time=%e, dt=%e, max(Pe)=%f' %(solver.iteration, solver.sim_time, timestep, max_Pe))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()


