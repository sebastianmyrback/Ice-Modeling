"""Running the basic ISMIP-HOM examples to test for convergence.

"""
__author__ = "Christian Helanow"
__mail__ = " helanow@iastate.edu"
__copyright__ = "Copyright (c) 2020 {}".format(__author__)

import os
from dolfin import *
import pstokes as ps
import matplotlib.pyplot as plt
import numpy as np

# python-plot

title_font = {'weight': 'bold', 'size': 20} # Set fonts for the plots
x_label_font = {'size': 14}
y_label_font = {'size': 14}
x_tick_font = {'size': 12}
y_tick_font = {'size': 12}

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

class TopBoundary(SubDomain):
      def inside(self, x, on_boundary):
          return near(x[1], height) and on_boundary

topboundary = TopBoundary()

# no slip at bed
bottomboundary = BottomBoundary()
sub_domains = MeshFunction("size_t", base_mesh, base_mesh.topology().dim() - 1)
sub_domains.set_all(0)
bottomboundary.mark(sub_domains, 1)

# GREJ SOM ÄNDRADES
topboundary.mark(sub_domains, 2)

# FIXME: hack for surface velocity
top_boundary = MeshFunction("size_t", base_mesh, 0)
top_boundary.set_all(0)
topboundary.mark(top_boundary, 1)

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

linear_solver = 'gmres'
preconditioner = 'none'

newton_parameters = Parameters('newton_solver') 
newton_parameters.add('linear_solver', linear_solver)  # the linear solver
newton_parameters.add('preconditioner', preconditioner)   # the preconditioner
newton_parameters.add('maximum_iterations',  50)  #maximum number of allowed iterations for Newton solver 
newton_parameters.add('relaxation_parameter', 1.0)  # If 1: a normal Newton solver. If between 0 and 1, the new solution will be a mix of the old and new one
newton_parameters.add('relative_tolerance',  1e-7)  #Newton solver will stop when the residual r divided by right hand side is smaller than this number


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


# GET SURFACE VELOCITY
# FIXME: this only works for the bottom for some reason....
#us = ps.get_boundary_function(u, sub_domains, 2)


# hack for surface vel
surf_indexes = top_boundary.where_equal(1)  # the MeshFunction marked with 1
surf_coors = base_mesh.coordinates()[surf_indexes]  # surface coordinates
dofs = W.dofmap().entity_dofs(base_mesh, 0)  # get dofs for points in W
coor_output = w.vector()[dofs].reshape(-1, 3)[surf_indexes]  # all values on coors
surf_vel = coor_output[:, 0:2]
speed = np.hypot(*surf_vel.T)
# plot with matplotlib
u_x = surf_vel[:,0].T       # Horisontal Surface Velocity
u_x_max = np.amax(u_x)      # Maximal Hor. Surf. Velocity
x_coord = surf_coors[:,0].T # Horisontal Coordinates

# 'True' for plot
if False:
    plt.quiver(*surf_coors.T, *surf_vel.T, speed)
plt.show()

# Horisontal Bed Coordinates
y_coord = -x_coord*np.tan(alpha) - height + bump_amp*np.sin(2*np.pi*x_coord/mesh_length)    # Equation for bed coordinates (QUICK FIX)


# EXPORT MAX VELOCITIES

# data = np.array([[mesh_length/1000], [u_x_max]])
# data = data.T
# datafile_path = "/Users/sebastianmyrback/Documents/KTH/Kandidatexamensarbete/Ice-Modeling/ahlkrona_postdoc/max_v.txt"
# with open(datafile_path, 'a+') as datafile_id:
#      # here you open the ascii file
#
#      np.savetxt(datafile_id, data, fmt=['%.4f', '%.4f'])
#      # here the ascii file is written.
#
# datafile_id.close()


# Maximal Surface Velocities for full-Stokes solution
max_u_arr = np.array([[5.0000, 11.7169],
[10.0000, 22.5055],
[20.0000, 46.4237],
[40.0000, 73.3268],
[80.0000, 94.6126],
[160.0000, 107.7985]])


# IMPORT BENCHMARKS

# Surface Velocities

b005 = np.array([np.loadtxt("aas1b005_surf.txt")[:, :2], np.loadtxt("ssu1b005.txt")[:, :2]])
b010 = np.array([np.loadtxt("aas1b010_surf.txt")[:, :2], np.loadtxt("ssu1b010.txt")[:, :2], np.loadtxt("spr1b010.txt")[:, :2]])
b020 = np.array([np.loadtxt("aas1b020_surf.txt")[:, :2], np.loadtxt("ssu1b020.txt")[:, :2], np.loadtxt("spr1b020.txt")[:, :2]])
b040 = np.array([np.loadtxt("aas1b040_surf.txt")[:, :2], np.loadtxt("ssu1b040.txt")[:, :2], np.loadtxt("spr1b040.txt")[:, :2]])
b080 = np.array([np.loadtxt("aas1b080_surf.txt")[:, :2], np.loadtxt("ssu1b080.txt")[:, :2], np.loadtxt("spr1b080.txt")[:, :2]])
b160 = np.array([np.loadtxt("aas1b160_surf.txt")[:, :2], np.loadtxt("ssu1b160.txt")[:, :2], np.loadtxt("spr1b160.txt")[:, :2]])

# Maximum Velocities

aas1_maxima = np.array([[5, np.amax(b005[0][:, 1])], [10, np.amax(b010[0][:, 1])], [20, np.amax(b020[0][:, 1])],
                             [40, np.amax(b040[0][:, 1])], [80, np.amax(b080[0][:, 1])], [160, np.amax(b160[0][:, 1])]])
ssu1_maxima = np.array([[5, np.amax(b005[1][:, 1])], [10, np.amax(b010[1][:, 1])], [20, np.amax(b020[1][:, 1])],
                             [40, np.amax(b040[1][:, 1])], [80, np.amax(b080[1][:, 1])], [160, np.amax(b160[1][:, 1])]])
spr1_maxima = np.array([[10, np.amax(b010[2][:, 1])], [20, np.amax(b020[2][:, 1])],
                             [40, np.amax(b040[2][:, 1])], [80, np.amax(b080[2][:, 1])], [160, np.amax(b160[0][:, 1])]])

# Plot Surface Velocities
# True for plot
if False:
    plt.figure(figsize=(12, 6),  dpi=150)
    plt.plot(x_coord/mesh_length, u_x, label='full-Stokes solution')
    plt.plot(b005[0][:, 0], b005[0][:, 1], label='aas1b model')
    plt.plot(b005[1][:, 0], b005[1][:, 1], label='ssu1b model')
    #plt.plot(b005[2][:, 0], b005[2][:, 1], label='spr1b model (not FEM)')
    plt.title('Length Scale: ' + str(int(mesh_length/1000)) + ' km', **title_font)
    plt.ylabel('Horisontal Surface Velocity [m/yr]', **y_label_font)
    plt.xlabel('Normalized Longitudal Dimension')
    plt.legend(prop={'size': 12})
    plt.show()


    # PLot maximum surface velocities

    plt.figure(figsize=(12, 6),  dpi=150)
    plt.plot(max_u_arr[:,0], max_u_arr[:,1], label='full-Stokes solution')
    plt.plot(aas1_maxima[:,0], aas1_maxima[:,1], label='aas1b model')
    plt.plot(ssu1_maxima[:,0], ssu1_maxima[:,1], label='ssu1b model')
    plt.plot(spr1_maxima[:,0], spr1_maxima[:,1], label='spr1b model')
    plt.title('Velocity Maxima at Surface', **title_font)
    plt.ylabel('Maximal Horisontal Surface Velocity [m/yr]', **y_label_font)
    plt.xlabel('Length of Domain [km]')
    plt.legend(prop={'size': 12})
    plt.show()

# write all the variables
xdmf_file = XDMFFile("{}/results.xdmf".format(savedir))
xdmf_file.parameters["flush_output"] = True
xdmf_file.parameters["functions_share_mesh"] = True
xdmf_file.parameters["rewrite_function_mesh"] = True
xdmf_file.write(u, 0)
xdmf_file.write(p, 0)


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
