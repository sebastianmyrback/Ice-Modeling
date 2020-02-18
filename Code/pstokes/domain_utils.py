# -*- coding: utf-8 -*-
"""
This module contains some meshing and boundary condition utils for the
pstokes module.

"""
import numpy as np
from dolfin import *


class PeriodicBoundaryBase(SubDomain):
    def __init__(self, mesh, map_tol=1e-10):
        """The extent of the domain. Must be rectangular for double
        periodic bcs.

        Parameters
        ----------
        mesh : dolfin Mesh
            The mesh for which the computations are to be periodic.
        """
        # FIXME: Due to pybind11 (intrinsics between Python MRO and
        # C++), the normal super inheritance can cause issues. Use
        # traditional explicit style instead.
        # super(PeriodicBoundaryBase, self).__init__(map_tol=map_tol)
        SubDomain.__init__(self, map_tol=map_tol)
        coors = mesh.coordinates()
        self.xmin = MPI.min(MPI.comm_world, coors[:, 0].min())
        self.xmax = MPI.max(MPI.comm_world, coors[:, 0].max())
        self.ymin = MPI.min(MPI.comm_world, coors[:, 1].min())
        self.ymax = MPI.max(MPI.comm_world, coors[:, 1].max())
        self.dim = mesh.geometric_dimension()
        self.map_tol = map_tol


class PeriodicBoundaryX(PeriodicBoundaryBase):
    """Maps the the x-coordinate boundaries to the same nodes.

    Parameters
    ----------
    mesh : dolfin Mesh
        The mesh for which the computations are to be periodic.

    """
    def __init__(self, mesh, map_tol=1e-7):
        # FIXME: Due to pybind11 (intrinsics between Python MRO and
        # C++), the normal super inheritance can cause issues. Use
        # traditional explicit style instead.
        # super(PeriodicBoundaryX, self).__init__(mesh, map_tol=map_tol)
        PeriodicBoundaryBase.__init__(self, mesh, map_tol=map_tol)
    # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return bool(near(x[0], self.xmin) and on_boundary)

    # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        y[0] = x[0] - self.xmax
        y[1] = x[1]
        if self.dim == 3:
            y[2] = x[2]


class PeriodicBoundaryXY(PeriodicBoundaryBase):
    def __init__(self, mesh, map_tol=1e-7):
        """The extent of the domain. Must be rectangular for double
        periodic bcs."""
        # FIXME: Due to pybind11 (intrinsics between Python MRO and
        # C++), the normal super inheritance can cause issues. Use
        # traditional explicit style instead.
        # super(PeriodicBoundaryXY, self).__init__(mesh, map_tol=map_tol)
        PeriodicBoundaryBase.__init__(self, mesh, map_tol=map_tol)

        if self.dim != 3:
            raise ValueError(
                "Dimension of mesh must be 3 for XY periodic BCs")

    def inside(self, x, on_boundary):
        """
        Return True if on left or bottom boundary AND NOT on one
        of the two corners (xmin, 1) and (1, 0).
        """
        return bool((near(x[0], self.xmin) or near(x[1], self.ymin)) and
                    (not ((near(x[0], self.xmin) and near(x[1], self.ymax))
                          or (near(x[0], self.xmax) and near(x[1], self.ymin))))
                    and on_boundary)

    def map(self, x, y):
        """
        Remap the values on the top and right sides to the bottom and left
        sides.
        """
        if near(x[0], self.xmax) and near(x[1], self.ymax):
            y[0] = x[0] - self.xmax
            y[1] = x[1] - self.ymax
            y[2] = x[2]
        elif near(x[0], self.xmax):
            y[0] = x[0] - self.xmax
            y[1] = x[1]
            y[2] = x[2]
        elif near(x[1], self.ymax):
            y[0] = x[0]
            y[1] = x[1] - self.ymax
            y[2] = x[2]
        else:
            y[0] = x[0]
            y[1] = x[1]
            y[2] = x[2]


def deform_mesh_to_geometry(mesh, surf_exp, bed_exp,
                            allow_extrapolation=True):
    """
    Vertically deforms the mesh to the surface and bed geometry.

    .. note:: Assumes that the base is at z=0

    Parameters
    ----------
    mesh : dolfin.cpp.mesh.Mesh
        The base mesh to be deformed.
    surf_exp : dolfin.function.expression.Expression
        The expression/function for the surface.
    bed_exp : dolfin.function.expression.Expression
        The expression/function for the bed.
    allow_extrapolation : bool
        If True (default), allow extrapolation for `surf_exp` and
        `bed_exp`.
    """
    dim = mesh.geometric_dimension()
    if (dim != 2) and dim != 3:
        raise ValueError("Can only deform 2- and 3-dimensional meshes.")
    # FIXME: make a local (deep) copy of the mesh coordinates. This
    # seems to prevent certain issues when deforming the mesh
    # (compared to when the mesh coordinates are deformed directly)
    # which gives the error:
    # *** Error:   Unable to create mesh entity.
    # *** Reason:  Mesh entity index -1 out of range.....
    # *** Where:   This error was encountered inside MeshEntity.cpp.
    coors = mesh.coordinates().copy()
    # dummy function space
    V = FunctionSpace(mesh, "CG", 1)
    surf = interpolate(surf_exp, V)
    bed = interpolate(bed_exp, V)
    if allow_extrapolation:
        surf.set_allow_extrapolation(True)
        bed.set_allow_extrapolation(True)
    min_height = MPI.min(MPI.comm_world, coors[:, dim-1].min())
    max_height = MPI.max(MPI.comm_world, coors[:, dim-1].max())
    mesh_height = max_height - min_height
    for x in coors:
        x[dim-1] = ((x[dim-1] / mesh_height) * (surf(*x) - bed(*x)))
        x[dim-1] = x[dim-1] + bed(x[0], x[1])
    # assign modified coordinates to mesh coordinates
    mesh.coordinates()[:] = coors


def is_periodic(coors, axis, tol):
    min_coor = MPI.min(MPI.comm_world, coors[:, axis].min())
    max_coor = MPI.max(MPI.comm_world, coors[:, axis].max())
    min_ind = np.where(coors[:, axis] == min_coor)
    max_ind = np.where(coors[:, axis] == max_coor)
    min_points = coors[min_ind]
    max_points = coors[max_ind]
    min_points = np.delete(min_points, axis, axis=1)
    max_points = np.delete(max_points, axis, axis=1)
    min_sort_ind = np.lexsort((min_points[:, 1], min_points[:, 0]))
    max_sort_ind = np.lexsort((max_points[:, 1], max_points[:, 0]))
    min_sorted = min_points[min_sort_ind]
    max_sorted = max_points[max_sort_ind]
    periodic = np.allclose(min_sorted, max_sorted, atol=tol, rtol=0)
    return periodic


def check_mesh_periodicity(mesh, axis=0, tol=1e-8):
    coors = mesh.coordinates()
    names = ['x', 'y', 'z']
    periodic = is_periodic(coors, axis, tol)
    if periodic:
        print("The {}-coordinates are periodic...".format(names[axis]))
    else:
        print("FAIL: The {}-coordinates are NOT periodic".format(names[axis]))
        
