## Preliminary observations of the linear solvers of FEniCS applied to artificial glacier model

1. Switching *p* towards 2 (linear) gave extremely low velocity profiles.
2. Velocity profile for eps = 1e-2 was higher than for eps = 1e-7. 

## Tables of working linear solvers and preconditioners applied to this problem

Here, an uppper limit of 50 000 linear solver iterations was set. Some combinations would possible converge for a greater number of allowed iterations, but we had to set some limit. Listed is all available linear solvers, however we are only interested in the Krylov solvers which are

* bicgstab        
* cg                           
* default                              
* gmres                 
* minres                             
* richardson                      
* tfqmr


### Preconditioners:

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Generalized minimal residual method | Algebraic multigrid 'amg'| No |
| Generalized minimal residual method | Default preconditioner  'default'| Yes |
| Generalized minimal residual method | Hypre algebraic multigrid 'hypre_amg' | No |
| Generalized minimal residual method | Hypre parallel incomplete LU factorization 'hypre_euclid'| No|
| Generalized minimal residual method | Hypre parallel sparse approximate inverse 'hypre_parasails'| No |
| Generalized minimal residual method | Incomplete Cholesky factorization 'icc' | No |
| Generalized minimal residual method | Incomplete LU factorization 'ilu'| Yes |
| Generalized minimal residual method | Jacobi iteration 'jacobi'| No |
| Generalized minimal residual method | No preconditioner 'none'| No |
| Generalized minimal residual method | PETSc algebraic multigrid 'petsc_amg'| No |
| Generalized minimal residual method | Successive over-relaxation 'sor'| No |

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Biconjugate gradient stabilized method 'bicgstab' | Algebraic multigrid 'amg'| Yes but incorrect solution (might also converge only for small grid sizes)|
| Biconjugate gradient stabilized method 'bicgstab' | Default preconditioner  'default'| Yes |
| Biconjugate gradient stabilized method 'bicgstab' | Hypre algebraic multigrid 'hypre_amg' | No |
| Biconjugate gradient stabilized method 'bicgstab' | Hypre parallel incomplete LU factorization 'hypre_euclid'| No|
| Biconjugate gradient stabilized method 'bicgstab' | Hypre parallel sparse approximate inverse 'hypre_parasails'| No |
| Biconjugate gradient stabilized method 'bicgstab' | Incomplete Cholesky factorization 'icc' | No |
| Biconjugate gradient stabilized method 'bicgstab' | Incomplete LU factorization 'ilu'| Yes |
| Biconjugate gradient stabilized method 'bicgstab' | Jacobi iteration 'jacobi'| No |
| Biconjugate gradient stabilized method 'bicgstab' | No preconditioner 'none'| No |
| Biconjugate gradient stabilized method 'bicgstab' | PETSc algebraic multigrid 'petsc_amg'| No, but does only ~2 Krylov iters. for each Newton iteration|
| Biconjugate gradient stabilized method 'bicgstab' | Successive over-relaxation 'sor'| No, but does only ~2 Krylov iters. for each Newton iteration |

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Conjugate gradient method 'cg' | All | No|

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Minimal residual method 'minres' | None | Yes|

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| MUltifrontal Massively Parallel Sparse direct Solver 'petsc' | amg | Yes|

### Linear solvers 

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Biconjugate gradient stabilized method 'bicgstab'| Incomplete LU factorization | Yes |
| Conjugate gradient method 'cg'| Incomplete LU factorization | No |
| Default | Incomplete LU factorization | Yes |
| Generalized minimal residual method 'gmres'| Incomplete LU factorization | Yes|
| Minimal residual method 'minres'| Incomplete LU factorization | No |
| MUltifrontal Massively Parallel Sparse direct Solver 'mumps'| Incomplete LU factorization | No ('Unknown Krylov Method')|
| PETSc built in LU solver 'petsc'| Incomplete LU factorization | No ('Unknown Krylov Method')|
| Richardson method 'richardson'| Incomplete LU factorization | No |
| Transpose-free quasi-minimal residual method 'tfqmr'| Incomplete LU factorization | Yes |
| Unsymmetric MultiFrontal sparse LU factorization 'umfpack'| Incomplete LU factorization | No ('Unknown Krylov Method')|


