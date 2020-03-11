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

Vi har nu lyckats implementera kvadratproblemet. Det fungerar bra, om än långsamt för direkta lösare. För linjära lösare fungerar det snabbt och bra upp till en viss mesh storlek där de linjära lösarna inte konvergerar.

