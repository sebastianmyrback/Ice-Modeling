# -*- coding: utf-8 -*-
"""
This module contains some helper functions for the pstokes module.

"""
import dolfin


class pStokesParameters(dolfin.Parameters):
    """Documentation for StokesParameters

    """
    def __init__(self):
        dolfin.Parameters.__init__(self)
        self.rename("p-Stokes_parameters")
        # power of p-Stokes: p = 4/3 => Glen's n = 3 for ice
        self.add('p', 4/3)
        # density of ice, kg/m^3
        self.add('rho', 910)
        # gravity, m/s^2
        self.add('gravity', 9.81)
        # ice softness Pa^3/s (\approx 1e-16 Pa^3/a)
        self.add('A0', 3.1688e-24)
        # secs per year
        self.add('spy', int(60*60*24*365.25))
        # epsilon to prevent infinite viscosity
        self.add('eps_reg', 10**(-12))
        # Whether to use MPa and year, instead of Pa and sec
        self.add('scale_units', True)
        self.working_units = {'stress': 'MPa', 'velocity': 'm/a'}
        self._meta_data = {
            'p': 'p-Stokes exponent',
            'rho': 'ice density [kg/m^3]',
            'gravity': 'gravity [m/s^2]',
            'A0': 'ice-softness parameter [Pa^-3/s]',
            'spy': 'second per year',
            'eps_reg': 'strain-rate regulation parameter (infinite viscosity)',
            'scale_units': 'Use [MPa] and [m/a] if true, otherwise [Pa] and [m/s]'}
        self.init_global_parameters()
        self.parameter_info()

    def init_global_parameters(self):
        pstokes_info("Initiating module pStokes")
        pstokes_info("Setting global parameter['allow_extrapolation'] = True",
                     indent=True)
        pstokes_info("If you want to change this back to FEniCS default, use:",
                     "parameter['allow_extrapolation'] = False", indent=True)
        dolfin.parameters['allow_extrapolation'] = True

    def parameter_info(self, indent=2):
        pstokes_info("\np-Stokes Parameters", indent=indent)
        pstokes_info("Working units are {} (stress) and {} (velocity)".format(
            self.working_units['stress'], self.working_units['velocity']),
                     indent=indent)
        print_strings = []
        fmt_string = "{{:<12}} | {{:<10{}}} | {{}}"
        for k in self.keys():
            val = self.__getitem__(k)
            fmt = ".3g" if isinstance(val, float) else ""
            p_string = fmt_string.format(fmt)
            s = p_string.format(k + ':', val, self._meta_data[k])
            print_strings.append(s)
        line_str = "-" * max(map(len, print_strings))
        print_strings = [line_str] + print_strings + [line_str]
        pstokes_info(*print_strings, indent=indent)

    def __setitem__(self, item, value):
        dolfin.Parameters.__setitem__(self, item, value)
        if item == 'scale_units' and value:
            self.working_units['stress'] = 'MPa'
            self.working_units['velocity'] = 'm/a'
            pstokes_info(
                "Changing working units to {} (stress) and {} (velocity)".format(
                    self.working_units['stress'], self.working_units['velocity']))
        elif item == 'scale_units' and not value:
            self.working_units['stress'] = 'Pa'
            self.working_units['velocity'] = 'm/s'
            pstokes_info(
                "Changing working units to {} (stress) and {} (velocity)".format(
                    self.working_units['stress'], self.working_units['velocity']))
            pstokes_info(
                "WARNING! Using SI unit values can lead to very poor results",
                indent=True)


def mpi_print(*args, **kwargs):
    """Print string with rank 0 using when MPI.

    Parameters
    ----------
    string : str
        String to be printed.
    """
    if dolfin.MPI.comm_world.Get_rank() == 0:
        print(*args, **kwargs)


def pstokes_info(*string_args, indent=False, **kwargs):
    """Appends the string to information string.

    Prints 'pStokes info:::' + string
    """
    si = "pStokes info :::: "
    if indent is True:
        ind_str = " " * len(si)
    elif isinstance(indent, int) and indent >= 0:
        ind_str = " " * indent
    else:
        indent = False
    string_args = [s.replace('\n', '\n ' + ind_str) for s in string_args]
    if len(string_args) == 1 and not indent:
        mpi_print(si, *string_args, **kwargs)
    elif len(string_args) == 1:
        mpi_print(ind_str, *string_args, **kwargs)
    elif indent:
        for s in string_args:
            mpi_print(ind_str, s, **kwargs)
    else:
        mpi_print(si, string_args[0], **kwargs)
        for s in string_args[1::]:
            mpi_print(ind_str, s, **kwargs)


def print_summary(w, params=None):
    """Prints velocity and pressure summary from Mixed Function.

    Parameters
    ----------
    w : dolfin MixedElement Function
        Containing velocity and pressure on mixed function space, W.
    """
    (u, p) = w.split(True)
    dim = w.geometric_dimension()
    vel = u.vector()[:].reshape((-1, dim))
    coor_labels = ['X', 'Y', 'Z']
    def var_range(x):
        xmin = dolfin.MPI.min(dolfin.MPI.comm_world, x.min())
        xmax = dolfin.MPI.max(dolfin.MPI.comm_world, x.max())
        return (xmin, xmax)
    vel_unit = "" if not params else params.working_units['velocity']
    p_unit = "" if not params else params.working_units['stress']
    pstokes_info("Summary output:")
    for i in range(dim):
        pstokes_info(
            "Velocity {}: ({:<.4g}, {:<.4g}) {}".format(
                coor_labels[i], *var_range(vel[:, i]), vel_unit),
            indent=True)
    pstokes_info(
        "Pressure: ({:<.4g}, {:<.4g}) {}".format(
            *var_range(p.vector()[:]), p_unit),
        indent=True)


# for projecting the stress field
def local_project(v, V, u=None):
    """Element-wise projection using LocalSolver

    Parameters
    ----------
    v : dolfin Function
        Function to be interpolated element wise.
    V : dolfin FunctionSpace
        The funciton space to be interpolated onto, most likely DG0.
    u : dolfin Function, optional
        Resulting function with interpolated values.
    """
    dv = dolfin.TrialFunction(V)
    v_ = dolfin.TestFunction(V)
    a_proj = dolfin.inner(dv, v_)*dolfin.dx
    b_proj = dolfin.inner(v, v_)*dolfin.dx
    solver = dolfin.LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u is None:
        u = dolfin.Function(V)
        solver.solve_local_rhs(u)
        return u
    else:
        solver.solve_local_rhs(u)
        return

