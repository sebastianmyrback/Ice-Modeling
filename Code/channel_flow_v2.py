"""
This script tries to reproduce numerical simulations/results from Adrian Hirn's dissertation "Finite Element Approximation of Problemsin Non-Newtonian Fluid Mechanics". In particular, Example 3 from Chapter 6 is considered.

The example uses an analytical solution of the p-Stokes equations on a unit square.

CH COMMENTS:

General
-------

I am assuming that you are using Python3
> from __future__ import print_function
is unnecessary, since this was one of the (major) differences when
Python changed major version from 2 to 3 (print is a real function in
Python3)

Also, if you look at the `fenics` module (i.e.
/path_to_fenics/fenics/__init__.py) that you are importing, it is only
a wrapper to import dolfin, so this is also unnecessary.

Boundary condition specific
---------------------------
Boundary condition works with supplying the kwarg `method =
"pointwise"`. However, it is then better to call/instantiate the
Corner class directly, and not marking points.

Parameters
----------
Regarding the analytical solution, remember that the parameters used herein should also be set in the pstokes module. Otherwise the results will depend on ice paramters.

So, the `p`-parameter must be set. This can be done with the
pstokes.pstokes_parameters (I hope). See below.

Similarly, the solution should be dependent on the choice of the
epsilon parameter in Hirn. This is implemented in the pstokes module
as part of the parameters with key 'eps_reg' (for epsilon
regularization). Also, you should check what measure Hirn uses for
`h`, but I suspect he used the the maximum cell diameter of the mesh.
This is a detail, but good to know.

The viscosity must be set as well, which is unity herein. In glaciology, the viscosity paramter is scaled by some 'ice softness' thingy, usually callet A0. I think that you can set this to unity, and it should be equivalent to Hirn's mu_0 = 1.

Also, due to the weird mixture of small an large numbers in the physical parameters in ice-sheet modeling, the equations are usually scale to be in MPa (pressure) and m/year (velocity). This scaling is the default in the pstokes module, but should probably *not* be used in this case. So turn the scaling of in the paramters (kwarg "scale_units").

To see what parameters are used in the pstokes module, do
> pstokes.pstokes_paramters.parameter_info()

Finally, you have to check what the body force of the system is. This
is chosen in a way so that it fulfills Hirn's example.
"""

import numpy as np
from dolfin import *
import pstokes as ps
import matplotlib.pyplot as plt

title_font = {'weight': 'bold', 'size': 20} # Set fonts for the plots
x_label_font = {'size': 14}
y_label_font = {'size': 14}
x_tick_font = {'size': 12}
y_tick_font = {'size': 12}


# SET UP PROBLEM PARAMETERS AND CONSTANTS

# CH: you don't seem to use these....
mu = 1             # kinematic viscosity
rho = 1            # density

b = -1.45
mu0 = 1.0
eps0 = 1.0
p = 1.3


# Create mesh and define function spaces
nx = 30
ny = 30

# Build mesh
mesh = RectangleMesh(
    Point(-0.5, -0.5), Point(0.5, 0.5),
    nx, ny)

# CH: determine mesh size by diameter of cell (max)
h = mesh.hmax()
eps = eps0*h**(2/p)

# CH: set the paramters for the p-Stokes module:
# Carreau model epsilon, this is just set a constant in the pstokes
# module, but is squared in Hirn's definition, so do that.
ps.pstokes_parameters['eps_reg'] = eps**2
ps.pstokes_parameters['p'] = p  # set the power-law paramter
ps.pstokes_parameters['A0'] = mu0  # unity pre-viscosity parameter, we believe this corresponds to Hirn's mu_0.
ps.pstokes_parameters['rho'] = rho
ps.pstokes_parameters['scale_units'] = False  # physical units m/s and Pa
# CH: etc....

# Expressions for velocity and pressure from A.Hirn
boundary_exp_v = Expression(('pow(pow(pow(x[0], 2) + pow(x[1], 2), 0.5), 7) * x[0]',
                             'pow(pow(pow(x[0], 2) + pow(x[1], 2), 0.5), 7) * (-x[1])'), degree=5)
boundary_exp_p = Expression('pow(pow(pow(x[0], 2) + pow(x[1], 2), 0.5), b) * x[0] * x[1]', b=b, degree=5)


# Boundary Conditions
class Boundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


class Corner(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], -0.5) and near(x[1], -0.5)


# no slip at boundary
boundary = Boundary()
sub_domains_v = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
sub_domains_v.set_all(0)
boundary.mark(sub_domains_v, 1)

# Define function spaces
V = VectorElement("Lagrange", mesh.ufl_cell(), 1)
P = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, V*P)
w = Function(W)

# set the actual Dirichlet conditions
noslip_bc_v = DirichletBC(W.sub(0), boundary_exp_v, sub_domains_v, 1)
boundary_p = DirichletBC(W.sub(1), boundary_exp_p, Corner(),
                         method='pointwise')
bcs = [noslip_bc_v, boundary_p]

# Solve the problem
# linear_solver = 'gmres'
# preconditioner = 'ilu'


linear_solver = 'mumps'
preconditioner = 'default'

newton_parameters = Parameters('newton_solver')
newton_parameters.add('linear_solver', linear_solver)  # the linear solver
newton_parameters.add('preconditioner', preconditioner)   # the preconditioner
newton_parameters.add('maximum_iterations',  50)  #maximum number of allowed iterations for Newton solver
newton_parameters.add('relaxation_parameter', 1.0)  # If 1: a normal Newton solver. If between 0 and 1, the new solution will be a mix of the old and new one
newton_parameters.add('relative_tolerance',  1e-7)  #Newton solver will stop when the residual r divided by right hand side is smaller than this number



# CH: note that this solver takes a kwarg "body_force". This will most
# likely have to be specified, together with other relevant paremters
# for the problem.
# You will have to calculate the specific bodyforce that fulfills the equations.
# body_force = ..... ????

# Set body force
(u, p) = w.split(True)

body_force = -div(ps.tau(u)) + grad(p)

# CH: final comment: Seems that IP works better than GLS, but not sure.
ps.viscosity_newton(w, W, bcs=bcs, body_force=body_force,
                    newton_solver_params=newton_parameters,
                    stabilization='IP',
                    symmetric_stokes=False)


# Split the mixed solution using deepcopy
# (needed for further computation on coefficient vector)
(u, p) = w.split(True)
u.rename("velocity", "Velocity in [{}]".format(
    ps.pstokes_parameters.working_units['velocity']))
p.rename("pressure", "Pressure in [{}]".format(
    ps.pstokes_parameters.working_units['stress']))


plt.figure(figsize=(12, 6),  dpi=120)
plt.title('Velocity Profile', **title_font)
cax = plot(u, cmap=plt.cm.jet, antialiased=True)
cbar = plt.colorbar(cax)
cbar.ax.set_title('Velocity profile')
plt.axis('tight')
plt.show()

plt.figure(figsize=(12, 6),  dpi=120)
plt.title('Pressure Profile', **title_font)
cax = plot(p, cmap=plt.cm.jet, antialiased=True)
cbar = plt.colorbar(cax)
plt.axis('tight')
plt.show()

plt.figure(figsize=(12, 6),  dpi=120)
# Horizontal component
u_hor = plt.subplot(2, 2, 1)
cax = plot(u[0], cmap=plt.cm.jet, antialiased=True)
cbar = plt.colorbar(cax)
cbar.ax.set_title('m/yr')
u_hor.set_title('$u_x$', **x_label_font)
plt.axis('tight')
plt.ylabel('z (m)', **y_label_font)

# Vertical component
u_ver = plt.subplot(2, 2, 2)
cax = plot(u[1], cmap=plt.cm.jet, antialiased=True)
cbar = plt.colorbar(cax)
u_ver.set_title('$u_y$', **x_label_font)
plt.axis('tight')

# Pressure
u_ver = plt.subplot(2, 2, 3)
cax = plot(p, cmap=plt.cm.jet, antialiased=True)
cbar = plt.colorbar(cax)
u_ver.set_title('$p$', **x_label_font)
plt.axis('tight')


plt.ylabel('y (m)', **y_label_font)
plt.xlabel('x (m)', **x_label_font)

plt.show()

