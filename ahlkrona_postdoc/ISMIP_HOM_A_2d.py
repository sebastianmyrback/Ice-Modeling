"""Running the basic ISMIP-HOM examples to test for convergence.

"""
__author__ = "Christian Helanow"
__mail__ = " helanow@iastate.edu"
__copyright__ = "Copyright (c) 2020 {}".format(__author__)

import os
from dolfin import *
import pstokes as ps
import matplotlib.pyplot as plt


savedir = "ISMIP_HOM_A_2d"
if not os.path.exists(savedir):
    os.mkdir(savedir)

# SET UP BASIC 2D MESH
nx = 40
ny = 20
mesh_length = 80000
height = 1000
base_mesh = RectangleMesh(
    Point(0, 0), Point(mesh_length, height),
    nx, ny)

# ISMIP-HOM A GEOMETRY
# inclination of surface, in radians
alpha = .5 * pi / 180
bump_amp = 500

# expressions for bed and surface
surface_str = '- x[0] * tan(alpha)'
bed_str = '{ss} + -height + amp * sin(2*pi*x[0]/L)'.format(
    ss=surface_str)
surface = Expression(surface_str, alpha=alpha, degree=2)
bed = Expression(bed_str, alpha=alpha, L=mesh_length, amp=bump_amp,
                        height=height, degree=2)

ps.pstokes_parameters['eps_reg'] = 1e-7
ps.pstokes_parameters['p'] = 4/3


####### SET UP STOKES PROBLEM
# BCS
class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) and on_boundary

# no slip at bed
bottomboundary = BottomBoundary()
sub_domains = MeshFunction("size_t", base_mesh, base_mesh.topology().dim() - 1)
sub_domains.set_all(0)
bottomboundary.mark(sub_domains, 1)

# periodic bc
PBC = ps.PeriodicBoundaryX(base_mesh, map_tol=1e-7)

# FUNCTION SPACES
V = VectorElement("Lagrange", base_mesh.ufl_cell(), 1)
P = FiniteElement("Lagrange", base_mesh.ufl_cell(), 1)
W = FunctionSpace(base_mesh, V * P, constrained_domain=PBC)
w = Function(W)                 # solution vector

# set the actual Dirichlet conditions
noslip_bc = DirichletBC(W.sub(0), Constant((0, 0)), sub_domains, 1)
bcs = [noslip_bc]

# DEFORM MESH
ps.deform_mesh_to_geometry(base_mesh, surface, bed)

# SOLVE THE PROBLEM

newton_parameters = Parameters('newton_solver') 
newton_parameters.add('linear_solver', 'gmres')  # the linear solver
newton_parameters.add('preconditioner', 'ilu')   # the preconditioner
newton_parameters.add('maximum_iterations',  50)  #maximum number of allowed iterations for Newton solver 
newton_parameters.add('relaxation_parameter', 1.0)  # If 1: a normal Newton solver. If between 0 and 1, the new solution will be a mix of the old and new one
newton_parameters.add('relative_tolerance',  1e-7)  #Newton solver will stop when the residual r divided by right hand side is smaller than this number


# solve equation
ps.viscosity_newton(w, W, bcs=bcs,
                    newton_solver_params=newton_parameters,
                    stabilization='GLS',
                    symmetric_stokes=False)

# Split the mixed solution using deepcopy
# (needed for further computation on coefficient vector)
(u, p) = w.split(True)
u.rename("velocity", "Velocity in [{}]".format(
    ps.pstokes_parameters.working_units['velocity']))
p.rename("pressure", "Pressure in [{}]".format(
    ps.pstokes_parameters.working_units['stress']))

# write all the variables
xdmf_file = XDMFFile("{}/results.xdmf".format(savedir))
xdmf_file.parameters["flush_output"] = True
xdmf_file.parameters["functions_share_mesh"] = True
xdmf_file.parameters["rewrite_function_mesh"] = True
xdmf_file.write(u, 0)
xdmf_file.write(p, 0)

# python-plot

title_font = {'weight': 'bold', 'size': 20} # Set fonts for the plots
x_label_font = {'size': 14}
y_label_font = {'size': 14}
x_tick_font = {'size': 12}
y_tick_font = {'size': 12}

plt.figure(figsize=(12, 6),  dpi=200)

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
u_ver.set_title('$u_z$', **x_label_font)
plt.axis('tight')

# Pressure
u_ver = plt.subplot(2, 2, 3)
cax = plot(p, cmap=plt.cm.jet, antialiased=True)
cbar = plt.colorbar(cax)
u_ver.set_title('$p$', **x_label_font)
plt.axis('tight')


plt.ylabel('z (m)', **y_label_font)
plt.xlabel('x (m)', **x_label_font)

plt.show()
