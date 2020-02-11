# Preliminary observations of the linear solvers of FEniCS applied to artificial glacier model

1. Switching *p* towards 2 (linear) gave extremely low velocity profiles.
2. Velocity profile for eps = 1e-2 was higher than for eps = 1e-7. 

## Tables of working linear solvers and preconditioners applied to this problem

Here, an uppper limit of 50 000 linear solver iterations was set. Some combinations would possible converge for a greater number of allowed iterations, but we had to set some limit. 

We will investigate two types of iterative solvers: Krylov solvers and multigrid solvers.

## Krylov Solvers

The, in FEniCS, available Krylov *linear solvers* are:

* bicgstab       |  Biconjugate gradient stabilized method      
* cg             |  Conjugate gradient method                   
* default        |  default Krylov method                       
* gmres          |  Generalized minimal residual method         
* minres         |  Minimal residual method                     
* richardson     |  Richardson method                           
* tfqmr          |  Transpose-free quasi-minimal residual method

The, in FEniCS, available Krylov *preconditioners* are:

* amg              |  Algebraic multigrid                       
* default          |  default preconditioner                    
* hypre_amg        |  Hypre algebraic multigrid (BoomerAMG)     
* hypre_euclid     |  Hypre parallel incomplete LU factorization
* hypre_parasails  |  Hypre parallel sparse approximate inverse 
* icc              |  Incomplete Cholesky factorization         
* ilu              |  Incomplete LU factorization               
* jacobi           |  Jacobi iteration                          
* none             |  No preconditioner                         
* petsc_amg        |  PETSc algebraic multigrid                 
* sor              |  Successive over-relaxation  


| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Generalized minimal residual method | Algebraic multigrid 'amg'| No |
| **Generalized minimal residual method** | **Default preconditioner  'default'**| **Yes**|
| Generalized minimal residual method | Hypre algebraic multigrid 'hypre_amg' | No |
| Generalized minimal residual method | Hypre parallel incomplete LU factorization 'hypre_euclid'| No|
| Generalized minimal residual method | Hypre parallel sparse approximate inverse 'hypre_parasails'| No |
| Generalized minimal residual method | Incomplete Cholesky factorization 'icc' | No |
| **Generalized minimal residual method** | **Incomplete LU factorization 'ilu'**| **Yes** |
| Generalized minimal residual method | Jacobi iteration 'jacobi'| No |
| Generalized minimal residual method | No preconditioner 'none'| No |
| Generalized minimal residual method | PETSc algebraic multigrid 'petsc_amg'| No |
| Generalized minimal residual method | Successive over-relaxation 'sor'| No |

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Biconjugate gradient stabilized method 'bicgstab' | Algebraic multigrid 'amg'| No|
| **Biconjugate gradient stabilized method 'bicgstab'** | **Default preconditioner  'default'**| **Yes** |
| Biconjugate gradient stabilized method 'bicgstab' | Hypre algebraic multigrid 'hypre_amg' | No |
| Biconjugate gradient stabilized method 'bicgstab' | Hypre parallel incomplete LU factorization 'hypre_euclid'| No|
| Biconjugate gradient stabilized method 'bicgstab' | Hypre parallel sparse approximate inverse 'hypre_parasails'| No |
| Biconjugate gradient stabilized method 'bicgstab' | Incomplete Cholesky factorization 'icc' | No |
| **Biconjugate gradient stabilized method 'bicgstab'** | **Incomplete LU factorization 'ilu'**| **Yes** |
| Biconjugate gradient stabilized method 'bicgstab' | Jacobi iteration 'jacobi'| No |
| Biconjugate gradient stabilized method 'bicgstab' | No preconditioner 'none'| No |
| Biconjugate gradient stabilized method 'bicgstab' | PETSc algebraic multigrid 'petsc_amg'| No, but does only ~2 Krylov iters. for each Newton iteration|
| Biconjugate gradient stabilized method 'bicgstab' | Successive over-relaxation 'sor'| No, but does only ~2 Krylov iters. for each Newton iteration |

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Conjugate gradient method 'cg' | All | No|

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Minimal residual method 'minres' | Hypre algebraic multigrid 'hypre_amg' | No|
| Minimal residual method 'minres' | Default preconditioner  'default' | No|
| Minimal residual method 'minres' | Hypre algebraic multigrid 'hypre_amg' | No|
| Minimal residual method 'minres' | Hypre parallel incomplete LU factorization 'hypre_euclid' | No|
| Minimal residual method 'minres' | Hypre parallel sparse approximate inverse 'hypre_parasails'| No|
| **Minimal residual method 'minres'** | **Incomplete Cholesky factorization 'icc'** | **Yes up to ~ 20x10 grid**|
| Minimal residual method 'minres' | Incomplete LU factorization 'ilu'| No |
| Minimal residual method 'minres' | Jacobi iteration 'jacobi' | No|
| **Minimal residual method 'minres'** | **No preconditioner 'none'** | **Yes** |
| Minimal residual method 'minres' | PETSc algebraic multigrid 'petsc_amg' | No|
| Minimal residual method 'minres' | Successive over-relaxation 'sor' | No|

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Richardson method 'richardson' | All | No|

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| Transpose-free quasi-minimal residual method 'tfqmr'| Hypre algebraic multigrid 'hypre_amg' | No|
| **Transpose-free quasi-minimal residual method 'tfqmr'** | **Default preconditioner  'default'** | **Yes**|
| Transpose-free quasi-minimal residual method 'tfqmr' | Hypre algebraic multigrid 'hypre_amg' | No|
| Transpose-free quasi-minimal residual method 'tfqmr' | Hypre parallel incomplete LU factorization 'hypre_euclid' | No|
| Transpose-free quasi-minimal residual method 'tfqmr' | Hypre parallel sparse approximate inverse 'hypre_parasails'| No|
| Transpose-free quasi-minimal residual method 'tfqmr' | Incomplete Cholesky factorization 'icc' | No |
| **Transpose-free quasi-minimal residual method 'tfqmr'** | **Incomplete LU factorization 'ilu'**| **Yes** |
| **Transpose-free quasi-minimal residual method 'tfqmr'** | **Jacobi iteration 'jacobi'** | **Yes**|
| **Transpose-free quasi-minimal residual method 'tfqmr'** | **No preconditioner 'none'** | **Yes up to ~30x15 grid** |
| Transpose-free quasi-minimal residual method 'tfqmr' | PETSc algebraic multigrid 'petsc_amg' | No|
| Transpose-free quasi-minimal residual method 'tfqmr' | Successive over-relaxation 'sor' | No|

### Working Combinations

| Linear Solver | Preconditioner| Converging Solution|
| ------------- |:-------------:| -----:|
| **Generalized minimal residual method** | **Incomplete LU factorization 'ilu'**| **Yes** |
| **Biconjugate gradient stabilized method 'bicgstab'** | **Incomplete LU factorization 'ilu'**| **Yes** |
| **Minimal residual method 'minres'** | **No preconditioner 'none'** | **Yes** |
| **Transpose-free quasi-minimal residual method 'tfqmr'** | **Incomplete LU factorization 'ilu'**| **Yes** |
| **Transpose-free quasi-minimal residual method 'tfqmr'** | **Jacobi iteration 'jacobi'** | **Yes**|
| **Transpose-free quasi-minimal residual method 'tfqmr'** | **No preconditioner 'none'** | **Yes up to ~30x15 grid** |



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


