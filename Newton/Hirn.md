## Hirn

Hirn föreslår en förbättringar av Newtons metod för att lösa p-Stokes ekvationer. 

### Steg

1. Stabilisering av det instabila FEM problemet. Här används LPS (vad används i FEniCS lösningen?)
2. De algebraiska problemet löses mha Newtons metod.
3. De linjära delproblemen löses med GMRES och multigrid som prekonditionerare. 

### Diskretisering

Diskretiseringen av problemet görs 


### Funderingar

* Hirn får p-Stokes att konvergera mha GMRES samt multigrid prekonditionerare, det lyckas inte vi med. Beror det på att han använder LPS medan vi använder GLS som stabilisering?


## Kvadratproblemet

Vi har nu lyckats implementera kvadratproblemet. Det fungerar bra, om än långsamt för direkta lösare. För vissa linjära lösare fungerar det snabbt och bra upp till en viss mesh storlek där de **linjära lösarna inte konvergerar** ("KSP_DIVERGED_BREAKDOWN
A breakdown in the Krylov method was detected so the method could not continue to enlarge the Krylov space. Could be due to a singlular matrix or preconditioner").

| Linear Solver + Preconditioner| 16x16 | 32x32 | 64x64 | 128x128 |
| --------------------:| -----:|------:| -----:| -----:|
| (gmres, ilu) | yes | yes | yes | no |
| (gmres, jacobi) | yes | yes | yes | no |
| (gmres, none) | yes | yes | yes | no |
| (bicgstab, ilu) | yes | yes | yes | no |
| (minres, none) | yes | yes | yes | yes |

* Är Newtons metod implementerad på exakt samma sätt som Hirn? I Hirn på sida. 169 skriver han

"The number
within the brackets represents the total number of iterations performed by the step-size
control and it equals the number l^star that appears in Algorithm 3.1 with lambda = 3/4. As initial
guess for Newton’s method on level l, we chose the FE solution corresponding to level l−1."

Gör verkligen vi detta? För han får konvergens till TOL = 10exp-11 på cirka 7-8 iterationer. För oss krävs det mer kring 20 stycken för TOL = 1exp-7.

Vidare, gällande Newtons metod implementerad av Hirn så skriver han såhär:

"Step 3 of Algorithm 3.1 includes the step-size control which is crucial when highly nonlinear
p-structure problems are solved via Newton’s method. In general, the Newton update zk
h
is weighted by the relaxation parameter lambda. The step-size control enables the globalization
of Newton’s method, i.e., the independence of the convergence with respect to the choice
of u0
h. If l^star = 0, then Algorithm 3.1 performs one full Newton cycle."

### Konvergensordning

För eps0 = 1, h = 1/nx, lambda = 0.6, och maximal gridsize 256:

[0.27125261 0.08736796 0.02327888 0.00582187] convergence rate vel
[0.63997459 0.22650198 0.02950423 0.00180077] convergence rate p

Detta säger absolut ingenting.. Vi får helt enkelt hoppas på att felet ligger i hur vi mäter normen.
