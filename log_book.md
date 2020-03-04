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
Under dagens möte diskuterades 'Krylov sub-spaces', Josefin bad oss undersöka beroendet av grid size för de olika lösarna. Undersökningen av lösningsmetoderna gav överensstämmande resultat med benchmark-datan. Josefin bad oss förklara olika lösare och varför de olika fungerar olika bra. Josefin uppmanade oss att läsa senaste IPCC-rapporten om kryosfären.

### 18/2 SPONTANMÖTE MED JOSEFIN OCH CHRISTIAN
Jonatan och Sebastian hade ett spontanmöte med Josefin och Christian då alla råkade befinna sig i Kräftriket samtidigt och flera frågor hade börjat ta form. Under mötet diskuterades vilken avgränsning som skulle göra gällandes linjära lösare då Sebastian upptäckt ett sätt att få ännu fler att konvergera. Vi berättade om hur vissa lösare inte konvergerade för specifika storlekar på elementen. I punktform kom vi fram till detta
* Behåll proportionerna på elementen.
* Begränsa studien till specifik storlek på tvärsnittet!
* Kolla hur konditionstalet beror på h(ena sidlängden på elementen)
* Jämför konditionstalen med p och om möjligt även epsilon
* Plocka ut egenvärden från enbart A-matrisen 
* Plotta konditionstal mot newton-iterationer

### 20/2 - Undersökning av hur antalet Krylov-iterationer beror på mesh
Jonatan har idag tagit fram plottar för hur antalet krylov-iterationer beror på storleken på elementen. Reslutaten visade att GMRES och bicgstab med ilu som prekonditionerare gav minst krylov-iterationer innan konvergens.

### 19/2 - Framtagning av fler Newton plottar
Sebastian har idag tagit fram en undersökning för hur konditionstalet beror på newton iterationerna. Ett minimum av konditionstalet finnes för lägre antal iterationer, medan peakar kan ses för konditionstal mellan 6-10 framför allt.

### 27/2 Möte med forskningsgruppen + möte med Josefin och Christian
Idag var var vi och presenterade vårt projekt för Josfin, Christian och Josefins två doktorander. Vi fick också höra vad de andra i gruppen forskar kring. Efter mötet pratade vi med Josefin och Christian om projektet och pratade om hur vi skulle ta fram konditionstalet för A-matrisen när Christian har skrivit kod som tar ut den. Vi tog också fram ett alternativt sätt att beräkna konditionstalet, baserat på singulärvärdena, som gav i allmänhet lägre konditionstal och ett annat beteende hos plottarna. 

### 4/3 Möte med Josefin och Christian

* Jämföra med exakt lösning, dvs litet epsilon.
* Använd hastigheten för randen enligt Hirn.
* Sätt en punkt i hörnet av trycket enligt Hirn.
