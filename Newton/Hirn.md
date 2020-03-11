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


### Kvadratproblemet

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
