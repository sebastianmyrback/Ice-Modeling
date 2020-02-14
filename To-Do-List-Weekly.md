### Vecka 6
**Jämförelse med ISMIP-HOP experiment B:**

Main focus will be to set up the code in a way that we can compare relevant parameters measured in the ISMIP-HOP experiment B. 
These can be shown in the *ISMIP-HOP* folder. The following results are benchmarked.
Huvudfokus kommer vara att modifiera koden för att kunna ta fram de relevanta resultat som går att jämföra med benchmarken. 
Dessa hittas i *ISMIP-HOP* mappen. Följande resultat är benchmarkade:

* **Surface velocity in the x direction for each domain length.** Gjort
* Dp for each domain length.
* **Maximum surface velocity value in x direction for each domain length.** Gjort
* Minimum surface velocity value in x direction for each domain length.
* Maximum dp value

Det vi behöver göra denna vecka är alltså att fundera på hur vi skall använda dessa benchmarks exakt för att jämföra med våra
lösningar. En idé kanske är att göra en medelvärdesbildning av dessa och jämföra med vårt resultat? 

Vi behöver också såklart ta fram alla ovanstående parametrar. 

### v.7-8

* Vi har nu lyckats ta fram plottar för konvergens av krylovlösare med prekonditionerare både mot antal newton iterationer samt antal krylov iterationer. De ser hyfsade ut.

De mer praktiska sakerna som behövs göras framöver är

* Behöver klura lite på hur man tar fram dessa påstådda multigrid *linjära lösare* som Josefin nämnt. Jag hittar bara multigrid prekonditionerare, som inte konvergerar här (t.ex. 'amg'). 
* Behöver lyckas hämta ut värdena för trycket vid glaciärbotten, för att kunna ta fram Dp = p - rho*g*h. 

De mer teoretiska sakerna som behövs förstås är

* Hur beror konvergens av linjära lösare/prekonditionerare på konditionstalet?
* Hur beror konvergens av lösningen på utseendet/egenskaperna hos vår systemmatris?
* Exakt hur interagerar Newton iterationerna med de linjära iterationerna?
