## Hirn

Hirn föreslår en förbättringar av Newtons metod för att lösa p-Stokes ekvationer. 

### Steg

1. Stabilisering av det instabila FEM problemet. Här används LPS (vad används i FEniCS lösningen?)
2. De algebraiska problemet löses mha Newtons metod.
3. De linjära delproblemen löses med GMRES och multigrid som prekonditionerare. 

### Diskretisering

Diskretiseringen av problemet görs 
