#### Frågor till Josefin ####
* Målet till nästa vecka (v.7) är att jämföra linjära lösare och prekonditionerare. Vi undrar hur det är tänkt att göras. 
Vår plan är att mäta konvergenshastigheten för de olika metoderna, genom att kanske använda benchmark-värdena? 
Eller behöver vi bara beräkna residualen för att mäta konvergenshastigheten? 

* Hur gör vi med benchmark-värdena? Ska vi medelvärdensbilda alla olika full-stokes-resultat? 

* Funkar vissa prekonditionerare endast med vissa linjära lösare?

* I pstokes.py används bara Krylov lösare vad det ser ut som. Skall vi testa andra typer av iterative lösare än Krylov?

* Försöker förstå hastighetsvektorn, för att kunna ta fram endast de värden som är på ytan.

