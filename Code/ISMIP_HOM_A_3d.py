"""Running the basic ISMIP-HOM examples to test for convergence.

"""
__author__ = "Christian Helanow"
__mail__ = " helanow@iastate.edu"
__copyright__ = "Copyright (c) 2020 {}".format(__author__)

import os
from dolfin import *
import pstokes as ps

# choose simulations type
mesh_type = "extruded"
# mesh_type = "unstructured"

savedir = "ISMIP_HOM_A_3d"
if not os.path.exists(savedir):
    os.mkdir(savedir)

# import mesh, either `extruded` or `unstructured`
xdmf_mesh = "meshes/{}.xdmf".format(mesh_type)
with XDMFFile(MPI.comm_world,
                     xdmf_mesh) as xdmf_infile:
    base_mesh = Mesh(MPI.comm_world)
    xdmf_infile.read(base_mesh)

# ISMIP-HOM A GEOMETRY
# set the gemometry parameters for the domain
# domain x,y length [m] (this is also hard coded in the mesh, so must
# conform to that for accurate domain representation)
mesh_length = 80000
# Height
height = 1000

# inclination of surface, in radians
alpha = .5 * pi / 180
bump_amp = 500

# expressions for bed and surface
surface_str = '- x[0] * tan(alpha)'
bed_str = '{ss} + -height + amp * sin(2*pi*x[0]/L) * sin(2*pi*x[1]/L)'.format(
    ss=surface_str)
surface = Expression(surface_str, alpha=alpha, degree=2)
bed = Expression(bed_str, alpha=alpha, L=mesh_length, amp=bump_amp,
                        height=height, degree=2)

####### SET UP STOKES PROBLEM
# BCS
class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 0.0) and on_boundary

# no slip at bed
bottomboundary = BottomBoundary()
sub_domains = MeshFunction("size_t", base_mesh, base_mesh.topology().dim() - 1)
sub_domains.set_all(0)
bottomboundary.mark(sub_domains, 1)

PBC = ps.PeriodicBoundaryXY(base_mesh, map_tol=1e-7)

# FUNCTION SPACES
V = VectorElement("Lagrange", base_mesh.ufl_cell(), 1)
P = FiniteElement("Lagrange", base_mesh.ufl_cell(), 1)
W = FunctionSpace(base_mesh, V * P, constrained_domain=PBC)
w = Function(W)                 # solution vector

# set the actual Dirichlet conditions
noslip_bc = DirichletBC(W.sub(0), Constant((0, 0, 0)), sub_domains, 1)
bcs = [noslip_bc]

# DEFORM MESH
ps.deform_mesh_to_geometry(base_mesh, surface, bed)

# SOLVE THE PROBLEM
newton_parameters = Parameters('newton_solver')
newton_parameters.add('linear_solver', 'mumps')
newton_parameters.add('preconditioner', 'default')
newton_parameters.add('maximum_iterations',  50)
newton_parameters.add('relaxation_parameter', 1.0)
newton_parameters.add('relative_tolerance',  1e-7)

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


