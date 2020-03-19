"""
This script tries to reproduce numerical simulations/results from Adrian Hirn's dissertation "Finite Element Approximation
of Problemsin Non-Newtonian Fluid Mechanics". In particular, Example 3 from Chapter 6 is considered.
The example uses an analytical solution of the p-Stokes equations on a unit square.
"""

import numpy as np
from dolfin import *
import pstokes_v2 as ps
import matplotlib.pyplot as plt
import sympy as sp

x1, x2 = sp.symbols('x[0] x[1]')

title_font = {'weight': 'bold', 'size': 20} # Set fonts for the plots
x_label_font = {'size': 14}
y_label_font = {'size': 14}
x_tick_font = {'size': 12}
y_tick_font = {'size': 12}


# Set up problem parameters and constants

mu = 1             # kinematic viscosity
rho = 1            # density

b = -1.45
mu0 = 1.0          # kinematic viscosity


# Solve the problem one time
def solve(nx, ny, eps0, linear_solver, preconditioner):

    # Booleans for different methods
    sympy_method = False
    Josefin_method = True

    simple_solution = False
    sinus_solution = True
    Hirn_solution = False

    if simple_solution or sinus_solution:
        p_param = 2.0
    if Hirn_solution:
        p_param = 1.3

    # Build mesh
    mesh = RectangleMesh(
        Point(-0.5, -0.5), Point(0.5, 0.5),
        nx, ny)

    # Determine mesh size by diameter of cell (max)
    h_alt = mesh.hmax()
    print(nx, 'nx')
    h = np.maximum(1/nx, 1/ny)
    print(h_alt, 'h_alt')
    print(h, 'h')
    eps = eps0*h_alt**(2/p_param)
    print(eps, 'current eps')

    # Carreau model epsilon, this is just set a constant in the pstokes
    # module, but is squared in Hirn's definition, so do that.
    ps.pstokes_parameters['eps_reg'] = eps**2
    ps.pstokes_parameters['p'] = p_param  # set the power-law paramter
    ps.pstokes_parameters['A0'] = mu0  # unity pre-viscosity parameter, we believe this corresponds to Hirn's mu_0.
    ps.pstokes_parameters['rho'] = rho
    ps.pstokes_parameters['scale_units'] = False  # physical units m/s and Pa

    # Expressions for velocity and pressure from A.Hirn
    if simple_solution:
        # v(x) = (y, -x)
        # p(x) = 0
        v_e = Expression(('x[1]', '-x[0]'), degree=5)
        p_e = Expression('0', degree=5)

    if sinus_solution:
        v_e = Expression(('sin(pi*x[0])*sin(pi*x[1])', 'sin(pi*x[0])*sin(pi*x[1])'), degree=5)
        # p(x) = 0
        p_e = Expression('0', degree=5)

    if Hirn_solution:
        # v(x) = |x| ^ 7 * (x1, -x2)
        v_e = Expression(('pow(hypot(x[0], x[1]), 7) * x[0]',
                                 'pow(hypot(x[0], x[1]), 7) * (-x[1])'), degree=5)
        # p(x) = |x|^b*x1*x2
        p_e = Expression('(x[0] == 0 && x[1] == 0) ?0.0 :pow(hypot(x[0], x[1]), b) * x[0] * x[1]')


    # Boundary Conditions
    class Boundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary

    class Corner(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0], -0.5) and near(x[1], -0.5)

    # No-slip at boundary
    boundary_v = Boundary()
    sub_domains_v = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    sub_domains_v.set_all(0)
    boundary_v.mark(sub_domains_v, 1)

    # Define function spaces
    V = VectorElement("Lagrange", mesh.ufl_cell(), 1)
    P = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
    W = FunctionSpace(mesh, V*P)
    w = Function(W)

    # Set the actual Dirichlet conditions
    noslip_bc_v = DirichletBC(W.sub(0), v_e, sub_domains_v, 1)
    boundary_p = DirichletBC(W.sub(1), p_e, Corner(), method='pointwise')
    bcs = [noslip_bc_v, boundary_p]

    # Solve the problem
    newton_parameters = Parameters('newton_solver')
    newton_parameters.add('linear_solver', linear_solver)  # the linear solver
    newton_parameters.add('preconditioner', preconditioner)   # the preconditioner
    newton_parameters.add('maximum_iterations',  50)  #maximum number of allowed iterations for Newton solver
    newton_parameters.add('relaxation_parameter', 0.6)  # If 1: a normal Newton solver. If between 0 and 1, the new solution will be a mix of the old and new one   Hirn puts 0.75
    newton_parameters.add('relative_tolerance',  1e-7)  #Newton solver will stop when the residual r divided by right hand side is smaller than this number

    # Set body force explicitly: (f = - div(S) + grad(p) in Hirn example 3, chapter 6)

    # Pressure: (Note the conditional statement to avoid singularities in origo)

    # Simple test problem
    p_x_cond = 0
    p_y_cond = 0


    """ Method 1 (Projection); not very good """
    p_anal = project(p_e, W.sub(1).collapse())
    u_anal = project(v_e, W.sub(0).collapse())
    #body_force = project(-div(ps.tau(u_anal)) + grad(p_anal), W.sub(0).collapse())

    """ Method 3: This one uses sympy to extract the explicit expressions.
        Currently the most reliable one, although it still produces strange
        results regarding the pressure. """

    if sympy_method:
        if Hirn_solution:
            v_symb_x = (x1**2 + x2**2)**7 * x1
            v_symb_y = -(x1**2 + x2**2)**7 * x2
        if sinus_solution:
            v_symb_x = sin(pi * x1) * sin(pi * x2)
            v_symb_y = sin(pi * x1) * sin(pi * x2)
        if simple_solution:
            v_symb_x = x2
            v_symb_y = -x1
        Dv11 = v_symb_x.diff(x1, 1)
        Dv12 = 0.5*(v_symb_x.diff(x2, 1) + v_symb_y.diff(x1, 1))
        Dv21 = 0.5*(v_symb_y.diff(x1, 1) + v_symb_x.diff(x2, 1))
        Dv22 = v_symb_y.diff(x2, 1)

        Dv_norm_sq = Dv11*Dv11 + Dv12*Dv12 + Dv21*Dv21 + Dv22*Dv22      # Inner product for tensors

        # Components of extra stress tensor with Carreau model (the factor 0.5 here
        # I cannot motivate), but it's there in Christians p-Stokes module)

        S11 = mu0*(eps**2 + 0.5*Dv_norm_sq) ** ((p_param-2)/2) * Dv11
        S12 = mu0*(eps**2 + 0.5*Dv_norm_sq) ** ((p_param-2)/2) * Dv12
        S21 = mu0*(eps**2 + 0.5*Dv_norm_sq) ** ((p_param-2)/2) * Dv21
        S22 = mu0*(eps**2 + 0.5*Dv_norm_sq) ** ((p_param-2)/2) * Dv22

        # Tensor divergence:
        dS1 = S11.diff(x1, 1) + S12.diff(x2, 1)
        dS2 = S21.diff(x1, 1) + S22.diff(x2, 1)

        # Body force: f = - div(S) + grad(p) ( = 0 for this problem)
        body_force_x = -dS1 + p_x_cond
        body_force_y = -dS2 + p_y_cond

    """ Method 4: Josefin's Matlab script"""
    if Josefin_method:
        # Enforce p(x=0,y=0) = 0. Might not be well-posed.
        #p_x_cond = sp.Piecewise((0, sp.Eq(x1, 0) & sp.Eq(x2, 0)), (p_x, True))
        #p_y_cond = sp.Piecewise((0, sp.Eq(x1, 0) & sp.Eq(x2, 0)), (p_y, True))


        # Hirn problem
        if Hirn_solution:
            body_force_x = '(x[0] == 0 && x[1] == 0) ? 0.0 :(29*pow(x[0],2)*x[1])/(20*pow((pow(x[0],2) + pow(x[1],2)),(69/40))) - x[1]/pow((pow(x[0],2) + pow(x[1],2)),(29/40)) - 7*pow(2,(2 - p_param/2))*x[0]*pow((pow(x[0],2) + pow(x[1],2)),(3/2))*(8*pow(x[0],2) +3*pow(x[1],2))*pow((357*pow(x[0],2)*pow(x[1],12) + 875*pow(x[0],4)*pow(x[1],10) + 1295*pow(x[0],6)*pow(x[1],8) + 1295*pow(x[0],8)*pow(x[1],6) + 875*pow(x[0],10)*pow(x[1],4) + 357*pow(x[0],12)*pow(x[1],2) + 2*pow(eps,2) + 65*pow(x[0],14) + 65*pow(x[1],14)),(p_param/2 - 1))-7*pow(2,(3 - p_param/2))*x[0]*pow((pow(x[0],2) + pow(x[1],2)),(13/2))*(p_param/2 - 1)*(8*pow(x[0],2) + pow(x[1],2))*(46*pow(x[0],2)*pow(x[1],2) + 65*pow(x[0],4) + 51*pow(x[1],4))*pow((357*pow(x[0],2)*pow(x[1],12) + 875*pow(x[0],4)*pow(x[1],10) + 1295*pow(x[0],6)*pow(x[1],8) + 1295*pow(x[0],8)*pow(x[1],6) + 875*pow(x[0],10)*pow(x[1],4) + 357*pow(x[0],12)*pow(x[1],2) + 2*pow(eps,2) + 65*pow(x[0],14) + 65*pow(x[1],14)),(p_param/2 - 2))'
            body_force_y = '(x[0] == 0 && x[1] == 0) ? 0.0 :(29*pow(x[1],2)*x[0])/(20*pow((pow(x[0],2) + pow(x[1],2)),(69/40))) - x[0]/pow((pow(x[0],2) + pow(x[1],2)),(29/40)) + 7*pow(2,(2 - p_param/2))*x[1]*pow((pow(x[0],2) + pow(x[1],2)),(3/2))*(8*pow(x[1],2) +3*pow(x[0],2))*pow((357*pow(x[0],2)*pow(x[1],12) + 875*pow(x[0],4)*pow(x[1],10) + 1295*pow(x[0],6)*pow(x[1],8) + 1295*pow(x[0],8)*pow(x[1],6) + 875*pow(x[0],10)*pow(x[1],4) + 357*pow(x[0],12)*pow(x[1],2) + 2*pow(eps,2) + 65*pow(x[0],14) + 65*pow(x[1],14)),(p_param/2 - 1))+7*pow(2,(3 - p_param/2))*x[1]*pow((pow(x[0],2) + pow(x[1],2)),(13/2))*(p_param/2 - 1)*(8*pow(x[1],2) + pow(x[0],2))*(46*pow(x[0],2)*pow(x[1],2) + 65*pow(x[1],4) + 51*pow(x[0],4))*pow((357*pow(x[0],2)*pow(x[1],12) + 875*pow(x[0],4)*pow(x[1],10) + 1295*pow(x[0],6)*pow(x[1],8) + 1295*pow(x[0],8)*pow(x[1],6) + 875*pow(x[0],10)*pow(x[1],4) + 357*pow(x[0],12)*pow(x[1],2) + 2*pow(eps,2) + 65*pow(x[0],14) + 65*pow(x[1],14)),(p_param/2 - 2))'

        # Sinus problem
        if sinus_solution:
            body_force_x = '(16*pow(pi,2)*pow((3*pow(pi,2) - 2*pow(pi,2)*cos(2*pi*(x[0] + x[1])) + 8*pow(eps,2) - pow(pi,2)*cos(2*pi*(x[0] - x[1]))),(p_param/2))*((pow(pi,2)*cos(pi*(3*x[0] - x[1])))/2 - pow(pi,2)*cos(3*pi*(x[0] + x[1])) - 8*pow(eps,2)*cos(pi*(x[0] + x[1])) + (3*pow(pi,2)*cos(pi*(x[0] - x[1])))/4 + (pow(pi,2)*cos(pi*(x[0] - 3*x[1])))/2 - (3*pow(pi,2)*cos(pi*(x[0] + 3*x[1])))/2 +   (pow(pi,2)*cos(pi*(3*x[0] + x[1])))/2 + (pow(pi,2)*cos(3*pi*(x[0] - x[1])))/4 + 4*pow(eps,2)*cos(pi*(x[0] - x[1])) + (p_param*pow(pi,2)*cos(pi*(x[0] - x[1])))/4 + (p_param*pow(pi,2)*cos(pi*(x[0] + 3*x[1])))/2 - (p_param*pow(pi,2)*cos(pi*(3*x[0] + x[1])))/2 - (p_param*pow(pi,2)*cos(3*pi*(x[0] - x[1])))/4 - p_param*pow(pi,2)*cos(pi*(x[0] + x[1])) + p_param*pow(pi,2)*cos(3*pi*(x[0] + x[1]))))/(pow(2,((3*p_param)/2))*pow((pow(pi,2)*cos(2*pi*(x[0] - x[1])) + 2*pow(pi,2)*cos(2*pi*(x[0] + x[1])) - 8*pow(eps,2) - 3*pow(pi,2)),2))'
            body_force_y = '(16*pow(pi,2)*pow((3*pow(pi,2) - 2*pow(pi,2)*cos(2*pi*(x[0] + x[1])) + 8*pow(eps,2) - pow(pi,2)*cos(2*pi*(x[0] - x[1]))),(p_param/2))*((pow(pi,2)*cos(pi*(3*x[0] - x[1])))/2 - pow(pi,2)*cos(3*pi*(x[0] + x[1])) - 8*pow(eps,2)*cos(pi*(x[0] + x[1])) + (3*pow(pi,2)*cos(pi*(x[0] - x[1])))/4 + (pow(pi,2)*cos(pi*(x[0] - 3*x[1])))/2 +   (pow(pi,2)*cos(pi*(x[0] + 3*x[1])))/2 - (3*pow(pi,2)*cos(pi*(3*x[0] + x[1])))/2 + (pow(pi,2)*cos(3*pi*(x[0] - x[1])))/4 + 4*pow(eps,2)*cos(pi*(x[0] - x[1])) + (p_param*pow(pi,2)*cos(pi*(x[0] - x[1])))/4 - (p_param*pow(pi,2)*cos(pi*(x[0] + 3*x[1])))/2 + (p_param*pow(pi,2)*cos(pi*(3*x[0] + x[1])))/2 - (p_param*pow(pi,2)*cos(3*pi*(x[0] - x[1])))/4 - p_param*pow(pi,2)*cos(pi*(x[0] + x[1])) + p_param*pow(pi,2)*cos(3*pi*(x[0] + x[1]))))/(pow(2,((3*p_param)/2))*pow((pow(pi,2)*cos(2*pi*(x[0] - x[1])) + 2*pow(pi,2)*cos(2*pi*(x[0] + x[1])) - 8*pow(eps,2) - 3*pow(pi,2)),2))'

        # Simple problem
        if simple_solution:
            body_force_x = '0'
            body_force_y = '0'

        body_force = Expression((body_force_x, body_force_y), p_param=p_param, eps=eps, degree=5)

    # Write body force as an Expression (careful with the mixing of variable names here)
    if sympy_method:
        body_force_x_code = sp.printing.ccode(body_force_x)
        body_force_y_code = sp.printing.ccode(body_force_y)
        body_force = Expression((body_force_x_code, body_force_y_code), degree=5)



    # Solve problem
    ps.viscosity_newton(w, W, bcs=bcs,
                        body_force=body_force,
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

    return [u, p, v_e, p_e, W, h]


# Compute the errors in velocity
def compute_errors_vel(u_e, u):

    # H1 seminorm
    E = errornorm(u_e, u, norm_type='H10')

    return E


# Compute the errors in the pressure
def compute_errors_p(p_e, p):

    # L2 norm
    E = errornorm(p_e, p, norm_type='L2')

    return E


# Compute convergence rates
def compute_convergence(u_error, p_error, h_arr):
    u_error_i = u_error[:-1]
    u_error_im1 = u_error[1:]
    p_error_i = p_error[:-1]
    p_error_im1 = p_error[1:]

    h_i = h_arr[:-1]
    h_im1 = h_arr[1:]
    print(u_error, 'u_error')
    print(u_error_i, 'u_error_i')
    print(u_error_im1, 'u_error_im1')
    print(h_i, 'h_i')
    print(h_arr, 'h_arr')
    print(h_im1, 'h_im1')

    # Convergence rate (see fenics book)
    r_vel = np.log(u_error_i/u_error_im1)/np.log(h_i/h_im1)
    r_p = np.log(p_error_i/p_error_im1)/np.log(h_i/h_im1)

    print(r_vel, 'convergence rate vel')
    print(r_p, 'convergence rate p')


# Plot velocity & pressure profiles
def plot_profiles(u, p, u_e, p_e, W):

    plt.figure(figsize=(12, 6),  dpi=120)
    plt.title('Velocity Profile', **title_font)
    cax = plot(u, cmap=plt.cm.jet, antialiased=True)
    cbar = plt.colorbar(cax)
    cbar.ax.set_title('$v$')
    plt.axis('tight')
    plt.show()

    plt.figure(figsize=(12, 6),  dpi=120)
    plt.title('Pressure Profile', **title_font)
    cax = plot(p, cmap=plt.cm.jet, antialiased=True)
    cbar = plt.colorbar(cax)
    cbar.ax.set_title('$p$')
    plt.axis('tight')
    plt.show()

    # Plot differences between analytical and simulated solutions.
    u_diff = Function(W.sub(0).collapse())
    u_e = interpolate(u_e, W.sub(0).collapse())
    u_diff.vector()[:] = u.vector()[:] - u_e.vector()[:]
    plt.figure(figsize=(12, 6),  dpi=120)
    plt.title('Difference between analytical and computed velocity', **title_font)
    cax = plot(u_diff, cmap=plt.cm.jet, antialiased=True)
    cbar = plt.colorbar(cax)
    cbar.ax.set_title('$\mathbf{u}_h - \mathbf{u}$')
    plt.axis('tight')
    plt.show()

    p_diff = Function(W.sub(1).collapse())
    p_e = interpolate(p_e, W.sub(1).collapse())
    p_diff.vector()[:] = p.vector()[:] - p_e.vector()[:]
    plt.figure(figsize=(12, 6),  dpi=120)
    plt.title('Difference between analytical and computed pressure', **title_font)
    cax = plot(p_diff, cmap=plt.cm.jet, antialiased=True)
    cbar = plt.colorbar(cax)
    cbar.ax.set_title('$p_h - p$')
    plt.axis('tight')
    plt.show()


# Plot how the error propagates towards the mesh size
def error_vs_grid(u_error, p_error, h_arr, eps0, lin_sol, precond, nx, ny):
    #hirn_p_eps1 = [2.66e-02, 1.14e-02, 5.47e-03, 2.67e-03, 1.31e-03, 6.46e-04]
    #hirn_v_eps1 = [1.47e-02, 6.47e-03, 2.99e-03, 1.43e-03, 7.00e-04, 3.46e-04]
    hirn_p_eps1 = [2.66e-02, 1.14e-02, 5.47e-03]    # If you want to compare our error to Hirns, use same
    hirn_v_eps1 = [1.47e-02, 6.47e-03, 2.99e-03]    # number of elements in these arrays as in 'nx_arr' in 'run()'

    plt.figure(figsize=(12, 6), dpi=120)
    plt.title(r'$ \Vert \nabla(v - v_h^{\epsilon})\Vert$ vs h for $\epsilon_0 = $' + str(eps0) + '. Maximum mesh size = ' + str((nx, ny)), **title_font)
    plt.plot(h_arr, u_error, label=(lin_sol, precond))
    #plt.plot(h_arr, hirn_v_eps1, label='Hirn´s velocity error')
    plt.xlabel('h', **x_label_font)
    plt.ylabel(r'$\Vert \nabla(v - v_h^{\epsilon})\Vert$', **y_label_font)
    plt.legend()
    plt.axis('tight')
    plt.show()

    plt.figure(figsize=(12, 6), dpi=120)
    plt.title(r'$ \Vert p - p_h^{\epsilon}\Vert$ vs h for $\epsilon_0 = $' + str(eps0) + '. Maximum mesh size = ' + str((nx, ny)), **title_font)
    plt.plot(h_arr, p_error, label=(lin_sol, precond))
    #plt.plot(h_arr, hirn_p_eps1, label='Hirn´s pressure error')
    plt.xlabel('h', **x_label_font)
    plt.ylabel('$\Vert p - p_h^{\epsilon}\Vert$', **y_label_font)
    plt.legend()
    plt.axis('tight')
    plt.show()


# main
def run():

    linear_solver = 'gmres'
    preconditioner = 'ilu'
    #linear_solver = 'mumps'
    #preconditioner = 'default'

    eps0 = 1.0  # constant epsilon parameter
    p_param = 2.0  # the 'p' in p-Stokes

    # Create arrays of mesh sizes
    # n = 16, 32, 64, 128, 256, 512

    nx_arr = np.array([16, 32, 64])
    ny_arr = nx_arr

    # Allocate memory
    u_error = np.zeros(nx_arr.shape)
    p_error = np.zeros(ny_arr.shape)
    h_arr = np.zeros(nx_arr.shape)
    #h_arr = 1/nx_arr

    i = 0

    for nx, ny in zip(nx_arr, ny_arr):
        print(nx_arr[i], ' current grid size')
        [u, p, u_e, p_e, W, h] = solve(nx, ny, eps0, linear_solver, preconditioner)
        u_error[i] = compute_errors_vel(u_e, u)
        p_error[i] = compute_errors_p(p_e, p)
        h_arr[i] = h
        i += 1

    if True:
        plot_profiles(u, p, u_e, p_e, W)

    print(u_error, 'u_error')
    print(p_error, 'p_error')

    # Plot error vs mesh size
    error_vs_grid(u_error, p_error, h_arr, eps0, linear_solver, preconditioner, nx_arr[-1], ny_arr[-1])

    # Compute convergence rates
    compute_convergence(u_error, p_error, h_arr)

    """ Ignore This """
    v_mumps_default = [0.024681763923475285, 0.02569289421021222, 0.02604939999265724, 0.02614152162869911,
                       0.026167591727440046, 0.026174692164698612]
    v_mumps_default = [0.026174692164698612, 0.026167591727440046, 0.02614152162869911, 0.02604939999265724, 0.02569289421021222, 0.024681763923475285]
    p_mumps_default = [0.06587862829124322, 0.04442976012421918, 0.03260378887200377, 0.03072719087671432,
                       0.03050679234384097, 0.030482045237350305]

    """--------------------------"""

    # plt.figure(figsize=(12, 6), dpi=120)
    # plt.title(r'$ \Vert \nabla(v - v_h^{\epsilon})\Vert$ vs h for $\epsilon_0 = $' + str(
    #     eps0) + '. Maximum mesh size = ' + str((int(nx_arr[-1]), int(ny_arr[-1]))), **title_font)
    # plt.plot(h_arr_, u_error_direct_eps1_reg_h, label=(linear_solver, preconditioner))
    # plt.plot(h_arr_, hirn_v_eps1, label='Hirn´s velocity error')
    # plt.xlabel('h', **x_label_font)
    # plt.ylabel(r'$\Vert \nabla(v - v_h^{\epsilon})\Vert$', **y_label_font)
    # plt.axis('tight')
    # plt.legend()
    # plt.show()
    #
    # plt.figure(figsize=(12, 6), dpi=120)
    # plt.title(r'$ \Vert p - p_h^{\epsilon}\Vert$ vs h for $\epsilon_0 = $' + str(eps0) + '. Maximum mesh size = ' + str(
    #     (int(nx_arr[-1]), int(ny_arr[-1]))), **title_font)
    # plt.plot(h_arr_, p_error_direct_eps1_reg_h, label=(linear_solver, preconditioner))
    # plt.plot(h_arr_, hirn_p_eps1, label='Hirn´s pressure error')
    # plt.xlabel('h', **x_label_font)
    # plt.ylabel('$\Vert p - p_h^{\epsilon}\Vert$', **y_label_font)
    # plt.axis('tight')
    # plt.legend()
    # plt.show()


run()