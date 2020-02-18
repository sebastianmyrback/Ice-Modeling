"""
Module containing the helper functions and physics for the p-Stokes
equations.
"""
__author__ = "Christian Helanow"
__mail__ = " helanow@iastate.edu"
__copyright__ = "Copyright (c) 2020 {}".format(__author__)

import numpy as np
from dolfin import *
from slepc4py import SLEPc

from .helpers import mpi_print, pstokes_info, pStokesParameters, \
    print_summary

# instantiate global and local (p-Stokes) parameters
pstokes_parameters = pStokesParameters()


def get_gravitational_force(dim):
    """Returns the gravitational force.

    Calculates the gravitational force of dimension `dim`. This is
    expressed in units [m/a^2], which is typical for glaciology.

    Parameters
    ----------
    dim : int
        The geometrical dimension.
    """
    # seconds per year squared
    spy2 = pstokes_parameters['spy']**2
    if pstokes_parameters['scale_units']:
        # acceleration in m/a^2 (downward direction)
        g_acc = -pstokes_parameters['gravity'] * spy2
        # scaled density 
        density = pstokes_parameters['rho'] / (10**6 * spy2)
    else:
        g_acc = -pstokes_parameters['gravity']
        density = pstokes_parameters['rho']
    force = g_acc * density
    bf_list = [0, ] * dim
    bf_list[-1] = force
    body_force = Constant(bf_list)
    return body_force


def mu_eff(u):
    """ Constitutive eq. p-Stokes
    
    Returns the effective viscosity of ice at a given velocity `u`.

    Parameters
    ----------
    u : dolfin Function
        Function containing the velocity.
    """
    n = 1 / (pstokes_parameters['p'] - 1)  # Glen's n
    if pstokes_parameters['scale_units']:
        # change working units to MPa and year (instead of Pa and sec).
        # FIXME: this is shady, since this parameter depends on the power.
        # Just changing the power is probably incorrect.
        A0 = pstokes_parameters['A0'] * 10**(6*n) * pstokes_parameters['spy']
    else:
        A0 = pstokes_parameters['A0']
    eps_reg = pstokes_parameters['eps_reg']
    eff_strain_squared = (0.5 * inner(sym(grad(u)), sym(grad(u)))
                          + Constant(eps_reg))
    _mu = 0.5 * (A0**(-1/n) * eff_strain_squared**((1-n)/(2*n)))
    return _mu


def tau(u):
    """Deviatoric stress tensor in terms of strain rates."""
    return 2 * mu_eff(u) * sym(grad(u))


def stokes_variation(w, w_test, body_force=None, symmetric=False):
    (u, p) = split(w)
    (v, q) = split(w_test)
    if not body_force:
        body_force = get_gravitational_force(w.geometric_dimension())
    # if symmetric, use sym-grad
    if symmetric:
        grad_test = sym(grad(v))
    else:
        grad_test = grad(v)
    # symmetric in div for stabilization (i.e. -/-)
    form = inner(tau(u), grad_test) * dx \
        - div(v) * p * dx \
        - div(u) * q * dx \
        - inner(body_force, v) * dx
    return form


# Compute solution
def viscosity_newton(w, W, bcs=None, body_force=None,
                     newton_solver_params=None,
                     stabilization='IP',
                     symmetric_stokes=False):
    """
    This is a manual version of the Newton solver. It uses the
    variable viscosity in the stabilization parameter, as it is
    supposed to be.
    INPUT:
    rel_tol: residual/relative tolerance.
    abs_tol: absolute tolerance. l2 norm of residual.
    iter_max: maximum Newton iterations. Also set if not specified
            directly.
    relax: relaxation parameter in Newton solver.
    If not set, the default parameters from IcePhysics::prm will
    be set.
    """
    mpi_print("****************************************")
    mpi_print("**** Solving the p-Stokes equations ****")
    mpi_print("****************************************")
    pstokes_parameters.parameter_info()
    if pstokes_parameters['scale_units']:
        pstokes_info("Using scaled units MPa and years.")
    else:
        pstokes_info("Using units Pa and seconds.")
    prm = NewtonSolver().parameters

    if newton_solver_params:
        prm.update(newton_solver_params)
    if not bcs:
        pass
    else:
        try:
            # copy bcs, locally: otherwise homogenize etc. affects subsequent
            # runs
            bcs = [DirichletBC(bc) for bc in bcs]
            for bc in bcs:
                bc.apply(w.vector())
                bc.homogenize()
        except TypeError:
            bcs = DirichletBC(bcs)
            bcs.apply(w.vector())
            bcs.homogenize()
    # set tolerance, max iterations and relaxation
    # parameter
    rel_tol = prm['relative_tolerance']
    abs_tol = prm['absolute_tolerance']
    iter_max = prm['maximum_iterations']
    relax = 1.0 if not prm['relaxation_parameter'] else prm['relaxation_parameter']
    # string to be used when printing out residual info
    res_string = ("Newton it: {} | " +
                  "r (abs): {:.3e} (tol: {:.3e}) | " +
                  "r (rel): {:.3e} (tol: {:.3e})")
    # ::::::::NEWTON SOLVER STARTS HERE:::::::::::
    w_trial = TrialFunction(W)
    w_test = TestFunction(W)
    (u, p) = split(w_trial)
    (v, q) = split(w_test)

    # specify the viscosity dependent stabilization
    visc_approx = mu_eff(w.sub(0))
    n_facet = FacetNormal(W.mesh())
    # SET GLS OR IP STABILIZATION
    if not body_force:
        body_force = get_gravitational_force(w_trial.geometric_dimension())
    if stabilization == 'GLS':
        pstokes_info("Using GLS stabilization")
        stab_visc = MinCellEdgeLength(W.mesh())  # Min or Max?
        GLS_visc_stab = -1 / (24*visc_approx) * dot(
            stab_visc*(grad(p) - body_force),
            stab_visc*grad(q))*dx(W.mesh())
        Complete_form = (
            stokes_variation(w_trial, w_test,
                             body_force=body_force,
                             symmetric=symmetric_stokes)
            + GLS_visc_stab)
    elif stabilization == 'IP':
        pstokes_info("Using IP stabilization")
        s = 2
        dim = W.mesh().geometric_dimension()
        if dim == 3:
            h_stab = MaxFacetEdgeLength(W.mesh())
        else:
            h_stab = MaxCellEdgeLength(W.mesh())
        IP_visc_stab = (
            1 / 24 *
            (avg(h_stab)**(s+1)/avg(visc_approx)*jump(grad(p), n_facet) *
             jump(grad(q), n_facet)*dS))
        Complete_form = (
            stokes_variation(w_trial, w_test,
                             body_force=body_force,
                             symmetric=symmetric_stokes)

            - IP_visc_stab)
    elif not stabilization:
        pstokes_info("Running unstabilized pStokes")
        Complete_form = stokes_variation(w_trial, w_test,
                                         body_force=body_force,
                                         symmetric=symmetric_stokes)
    else:
        msg = (
            "{} is not implemented as a stabilization type with"
            " the viscosity_newton method.")
        msg.format(stabilization)
        raise TypeError(msg)

    Complete_form = action(Complete_form, w)

    # calculate the jacobian for this case
    # is this really necessary??? Yes!
    J = derivative(Complete_form, w, w_trial)
    # define vector for incrementing the solution
    U_inc = Function(W)
    # dummy value to start iteration
    rel_res = 1.1 * rel_tol
    abs_res = 1.1 * abs_tol
    nIter = 0
    while (rel_res > rel_tol and abs_res > abs_tol) and nIter < iter_max:
        nIter += 1
        A, b = assemble_system(J, -Complete_form,
                               bcs, keep_diagonal=False)
        # specify the type of linear solver
        linear_solver = prm['linear_solver']

        preconditioner = prm['preconditioner']

        print('linear solver:' + str(linear_solver) + '. Preconditioner: ' + str(preconditioner))

        # Check symmetry
        Jac = PETScMatrix()
        assemble(J, Jac)
        #print('Is symmetric', np.linalg.norm(Jac.array() - Jac.array().T) < 1E-10)
        
        #print("Condition number is:",np.linalg.cond(Jac.array()))


        if False:

            # Checking condition number
            print("Assembling svd and eigenvalue problem")

            JacCopy = PETScMatrix()
            assemble(J, tensor=JacCopy)
            print("Computing singular values...")

            Svdsolver = SLEPc.SVD()
            Svdsolver.create()
            Svdsolver.setOperator(JacCopy.mat())
            Svdsolver.setType(Svdsolver.Type.LAPACK)
            Svdsolver.solve()

            nconv = Svdsolver.getConverged()
            if nconv > 0:
                vsvd, usvd = JacCopy.mat().getVecs()
                print(range(nconv))
                for i in range(nconv):
                    sigma = Svdsolver.getSingularTriplet(i, usvd, vsvd)
                    print("Singular value #", i, ": ", sigma)
            else:
                print("svd solver dit not converge")

            print("Computing eigenvalues...")

            eigensolver = SLEPcEigenSolver(JacCopy)
            eigensolver.parameters['solver'] = 'lapack'
            eigensolver.solve()

            # Printing eigenvalues

            matrixsize = np.shape(JacCopy.array())[0]
            for i in range(matrixsize):
                r, c, rx, cx = eigensolver.get_eigenpair(i)
                print("Eigenvalue #", i, ": ", r, "+i", c)


        """ -------- DIRECT/KRYLOV SOLVER --------------"""
        if linear_solver == 'mumps':
            solve(A, U_inc.vector(), b, linear_solver)

        else:
            MAX_ITERS = 10000
            solver = PETScKrylovSolver(linear_solver, preconditioner)
            solver.ksp().setGMRESRestart(MAX_ITERS)
            #solver.ksp().setMINRESRestart(MAX_ITERS)
            solver.parameters['relative_tolerance'] = 1E-2  # iterations stop if relative residual is smaller than this
            solver.parameters['absolute_tolerance'] = 1E-10  # iterations stop if residual is smaller than this
            solver.parameters['maximum_iterations'] = MAX_ITERS  # maximum allowed iterations of linear solver
            solver.parameters['monitor_convergence'] = True  # just to see what's going on
            solver.parameters['report'] = True
            # info(solver.parameters, 1)  # if you want to see what parameters there are to set
            solver.solve(A, U_inc.vector(), b)

        # assemble abs residual after first solve
        if nIter == 1:
            residual0 = b.norm('l2')
        # calculate residuals
        abs_res = b.norm('l2')
        rel_res = abs_res/residual0
        # update solution
        w.vector()[:] += relax*U_inc.vector()
        if MPI.comm_world.Get_rank() == 0:
            mpi_print(res_string.format(nIter, abs_res, abs_tol, rel_res, rel_tol))
        no_newton_iter = nIter
    # check for convergence
    if nIter == iter_max and rel_res > rel_tol:
        newton_convergence = False
    else:
        newton_convergence = True
    print_summary(w, params=pstokes_parameters)


def normal_newton(w, W, bcs=None, body_force=None,
                  newton_params=None,
                  adaptive_error=False,
                  adaptive_tol=1e-5):
    # trial and test
    dim = W.mesh().geometric_dimension()
    w_trial = TrialFunction(W)
    w_test = TestFunction(W)
    F = stokes_variation(w, w_test, body_force=body_force)
    J = derivative(F, w, w_trial)  # jacobian
    problem = NonlinearVariationalProblem(F, w, bcs=bcs, J=J)
    if adaptive_error:
        # M = inner(tau(w.sub(0))-w.sub(1)*Identity(dim),
        #           tau(w.sub(0))-w.sub(1)*Identity(dim))*dx()
        M = w[2]* dx()
        solver = AdaptiveNonlinearVariationalSolver(problem, M)
    else:
        solver = NonlinearVariationalSolver(problem)
    # parameters
    prm = solver.parameters

    if not newton_params:
        newton_params = Parameters('newton_solver')
        newton_params.add('linear_solver', 'mumps')
        newton_params.add('preconditioner', 'default')
        newton_params.add('maximum_iterations', 25)
        newton_params.add('relaxation_parameter', 0.7)
        newton_params.add('relative_tolerance', 1e-7)
    if adaptive_error:
        prm['nonlinear_variational_solver']['newton_solver'].update(newton_params)
        prm["error_control"][
            "dual_variational_solver"]["linear_solver"] = newton_params['linear_solver']
        solver.solve(adaptive_tol)
    else:
        prm['newton_solver'].update(newton_params)
        solver.solve()

