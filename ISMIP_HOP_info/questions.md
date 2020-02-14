#### Frågor till Josefin ####

### v.7
* Målet till nästa vecka (v.7) är att jämföra linjära lösare och prekonditionerare. Vi undrar hur det är tänkt att göras. 
Vår plan är att mäta konvergenshastigheten för de olika metoderna, genom att kanske använda benchmark-värdena? 
Eller behöver vi bara beräkna residualen för att mäta konvergenshastigheten? 

* Hur gör vi med benchmark-värdena? Ska vi medelvärdensbilda alla olika full-stokes-resultat? 

* Funkar vissa prekonditionerare endast med vissa linjära lösare?

* I pstokes.py används bara Krylov lösare vad det ser ut som. Skall vi testa andra typer av iterative lösare än Krylov?

* Försöker förstå hastighetsvektorn, för att kunna ta fram endast de värden som är på ytan.

### v.8-9

* Vi har nu lyckats ta fram plottar för konvergens av krylovlösare med prekonditionerare både mot antal newton iterationer samt antal krylov iterationer. De ser hyfsade ut.

De mer praktiska sakerna som behövs göras framöver är

* Behöver klura lite på hur man tar fram dessa påstådda multigrid *linjära lösare* som Josefin nämnt. Jag hittar bara multigrid prekonditionerare, som inte konvergerar här (t.ex. 'amg'). 
* Behöver lyckas hämta ut värdena för trycket vid glaciärbotten, för att kunna ta fram Dp = p - rho*g*h. 

De mer teoretiska sakerna som behövs förstås är

* Hur beror konvergens av linjära lösare/prekonditionerare på konditionstalet?
* Hur beror konvergens av lösningen på utseendet/egenskaperna hos vår systemmatris?
* Exakt hur interagerar Newton iterationerna med de linjära iterationerna?
