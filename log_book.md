## Loggbok Kandidatexamensarbete

### 27/1 - Möte med Josefin
Vi hade möte Med vår handledare Josefin om hur projektet skulle avgränsas och utformas. I samband med detta skissades en projektplan med preliminära deadlines.

### 28/1 - Projektplan skrevs
Projektplanen skrevs och skickades till Josefin för vidare revidering.

### 31/1 - Fippla med kod
Idag började vi kolla på iskoden som modellerar ett teoretiskt istäcke som rör sig utefter en sinusformad berggrund. Vi ändrade på ett par parameterar och försöker få oss en uppfattning om hur koden funkar.
En jämförelse för vilka linjära lösare/prekonditionerare som funkar på detta problem gjordes.

### 4/2 - Översikt

Vi har börjat sätta oss in i ISMIP-HOP B experimentet, och med detta har frågor uppstått. Dessa har samanställts i dokumentet *questions*, och handlar om själva felanalysen. Tagit fram några preliminära residual-norm plottar.

### 5/2 - Möte med Josefin 

Jonatan och Sebastian hade några frågor specifikt angående ISMIP-HOP experimentet samt detaljer i koden. Vi fick dessa klargjorda. Efteråt satt vi och började lyckas få fram plottar där vi undersöker residualnormen mot antalet iterationer (både Newton och Krylov iterationer). 

### 13/2 - Möte med Josefin
Under dagens möte diskuterades 'Krylov sub-spaces', josefin bad oss undersöka beroendet av grid size för de olika lösarna. Undersökningen av lösningsmetoderna gav överensstämmande resultat med benchmark-datan. Josefin bad oss förklara olika lösare och varför de olika fungerar olika bra. Josefin uppmanade oss att läsa senaste IPCC-rapporten om kryosfären.

### INSERT SPONTANMÖTE MED JOSEFIN OCH CHRISTIAN

### 19/2 - Framtagning av fler Newton plottar
Jag har idag tagit fram en undersökning för hur konditionstalet beror på newton iterationerna. Ett minimum av konditionstalet finnes för lägre antal iterationer, medan peakar kan ses för konditionstal mellan 6-10 framför allt. 
