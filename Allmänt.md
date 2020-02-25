# Översikt

Tack för mötet idag! Här kommer att-göra-listan med tillhörande
länkar/instruktioner. Ibland är instruktionerna på engelska, eftersom
jag ursprungligen skickat dem till min engelsktalande doktorand.

* Installera fenics, https://fenicsproject.org/download/	Jag tycker
det är lättast att göra det via anaconda

-----

* I denna länk, i sektionen "Examples" så finns följande exempel som ni
kan prova

ft01_poisson.py
ft02_poisson_membrane.py
ft03_heat.py
ft04_heat_gaussian.py
ft05_poisson_nonlinear.py
ft06_elasticity.py
ft07_navier_stokes_channel.py
ft08_navier_stokes_cylinder.py

but there are a couple of things to fix as they are restructuring fenics:

- The command "interactive()" has been taken away. I imported matplotlib
and used plt.show() instead
- wherever it says vector().array(), it should instead be
vector().get_local()
- In tutorial 6, the nabla_div operator has not been imported. It can be
imported like this: from ufl import nabla_div
- In tutorial 6 I also didnt manage to make the plot work at all, but it
saves some vtu-files, which you can look at using paraview.
https://www.paraview.org/	Tomas can install paraview for you. Also
ask him about matlab.
- In tutorial 8, comment out everything with "progress". It takes a lot
of time to run, so I changed from range(num_steps) to range(500)

The text that goes with the tutorials is this one:
https://fenicsproject.org/pub/tutorial/pdf/fenics-tutorial-vol1.pdf

-----

* När ni fått lite feeling för fenics så kan ni testa att köra iskoden.
Jag kommer mejla den nu i veckan, jag ska städa den lite först.

-----

* För att läsa på är denna bok av Larson och Bengzon bra:
http://matematicas.unex.es/~coco/Modelizacion/FEMLarson_Bengzon.pdf

I kapitel 6 skriver de om metoder för att lösa linjära system, och i
kapitel 12 skriver de om hur man använder FEM för fluidmekanik-problem

-----

* Kolla in de linjära lösare och prekonditionerare som finns i fenics:
https://fenicsproject.org/pub/tutorial/html/._ftut1017.html

-----

* Ta en titt på datan (bedmap och bedmachine) och fundera på hur ni kan
titta på dem. Själv har jag använt bedmap men inte bedmachine.
Bedmachine är i formatet "netcdf" som är ett format som ofta används i
klimat och i biologi för stor data

-----

* Installera gmsh, i min erfarenhet verkar det som att man gör bäst i
att inte välja senaste versionen, utan en gammal. Jag kör med den som är
längst ner på denna sida: http://gmsh.info/bin/Linux/older/

Just download it, you will get a folder with a bin and a
share folder, and in the bin folder there is an executable.

When you open it you see the gmsh GUI. They have a
graphical interface for creating meshes where you can draw points and
lines, which will result in a file describing the domain, that has the
extension .geo.  The GUI used to drive
me a little bit crazy but it can be nice to try to get a feeling for how
it works. While you click around there will be a .geo file automatically
created, saving the info on how you want your domain.

What I used to do in the end was to write the .geo files directly
without clicking around in the GUI. I think I have a file that creates a
domain from data-points, I will see if I can find it, somewhere (den här
ska jag fortfarande leta upp).

Once you have the geo file you can use gmsh to mesh it, which will
create a .msh file containing information about where your vertices are,
and how they are sorted into elements.

That .msh file you can convert to fenics format using the tool
dolfin-convert by writing this in the terminal:

dolfin-convert filename.msh filename.xml

And then you can load it into fenics

-----

* För att visualisera simuleringar brukar jag använda ParaView:
https://www.paraview.org/ Men jag tror det räcker om ni visualiserar mha
python/fenics

-----

Jag tänker att det som vore jättebra att ha i rapporten är 1) en tabell
som liknar den i pappret "Optimal numerical solvers...", i Fig 2, för en
viss glaciär, 2) en text/teori som förklarar varför vissa lösare är
bättre än andra (baserat på egenskaperna hos systemet Ax=b som kommer
från p-Stokes) 3)  en studie om hur antal iterationer och konditionstal
varierar om ni byter geometri eller randvillkor (t ex om ni gör så att
isen glider mer mot underlaget), i den studien kan ni kanske välja den
lösaren+prekonditioneraren som visade sig vara bäst i uppgift 1.

Det blev ett väldigt lång lista, men ni behöver absolut inte vara
stressade över att ni måste göra allt detta redan nu!

Hälsningar,

Josefin

--

## Möte med Josefin 27/1

Testa olika linjära lösare på teoretisk sinusvåg (olika grid-size, cluster 3D) och ge en översikt om fördelar och nackdelar (nerskrivet och klart 19/2)
Pröva linjära lösare på Arolla(schweizisk glaciär, botten är frusen, en del börjar slidea) Gör samma analys som i 1) Klart 26/2
Kolla om matriserna är symmetrisk för systemet ovan. Undersök konditionstalet. Positivt definit?
Undersöka Newton-lösaren (Ax = b, A ej linjär utan beror på x) för teoretiska problem och övergå(möjligtvis, förhoppningsvis!!) till riktig data. Läs i Larsson & Bengzon om icke-linjära lösare för info. Be Josefin skicka killen med kodens kod(klart 25/3).
Kolla möjligtvis på riktig geometri(26-31/3)


Josefins artificiella iskod: 

Hej!
 
Tack! Insåg när jag skrev detta mejl att det blev hemskt långt, så jag delar upp det i punktform:

Planen ser bra ut! Några småsaker:


Titlen är jättebra!
Ser att det står "måna" istället för "mån" på ett ställe
I första meningen bör det nog vara applicerade istället för applicerad, eftersom ekvationer är i plural
Förslag ista meningen: byt lösningar till linjära lösare. Och om ni vill fokusera mkt på Newton-lösaren kan ni kanske ta bort ordet primärt :)
 
*Jag kom på att om ni vill komma in på Newton-metoden snabbt så kan ni ju alltid börja med den efter ni gjort ISMIP-HOM experimentet, så kan ni göra Arolla-glaciären efter Newton om ni vill.*
 
*Alla benchmark experiment som heter något med ISMIP-HOM finns beskrivna här: http://homepages.ulb.ac.be/~fpattyn/ismip/*
 
Den nya koden finns här: https://bitbucket.org/christianhelanow/ahlkrona_postdoc/src/master/ Kommer ni åt detta repositorium? Annars ska jag be Christian ge er rättigheter. Jag håller på att lägga till lite grejer så ni lätt kan byta linjära lösare och prekonditionerare, jag mejlar igen när jag är klar med det.  Filen ni ska köra är den som heter ismip_hom_a_2d.py
 
 
* Angående metoden för att få Newtons metod att fungera bättre, så är den del av en doktorsavhandling av en tysk kille som heter Adrian Hirn. Sektion 6.2 i http://archiv.ub.uni-heidelberg.de/volltextserver/13173/1/dissertation_hirn.pdf handlar om konvergens av Newtons metod för p-Stokes. Den här avhandlingen är mycket teknisk men ni kan kolla lite på de numeriska experimenten, och hur han formulerar själva metoden i praktiken. Saker att göra här är:
 
 
Göra om hans experiment (kanske inte alla, vi kan prata om det när det närmarr sig) för att se att ni får samma resultat. Christian Helanow som börjar jobba i min grupp den 17:e februari har implementerat detta i koden jag skickade, men tagit bort det igen. När han kommer kan han lägga tillbaka det.
Kolla experimentellt om hans upper bound för konditionstalet (eq. 6.14) stämmer
Gradvis ändra på hans experiment så det med och mer liknar ISMIP-HOM B (isen som åker över en sinusformad mark) och se om metoden funkar lika bra på detta exempel.
Om den inte fungerar lika bra - försök komma på en modifikation av metoden som gör att den funkar för is.
 	   
Om denna metod går att få att fungera bra på is vore det ett riktigt bra bidrag till ismodellerings-communityn!
 
Som ni ser i ekvation 9.20 i bengzon coh larson http://matematicas.unex.es/~coco/Modelizacion/FEMLarson_Bengzon.pdf så får man om man använder Newtons metod inte systemet Ax=b utan istället Jd=r. I slutändan är det egentligen alltså matrisen J som man behöver kolla konditionstal och symmetri hos. Jag mejlar er koden för hur ni kan kolla konditionstal och symmetri i nästa mejl.
 
För att kolla vilka grejer ni kan ändra på för den linjära lösaren kan ni använda detta kommando:
			info(LinearVariationalSolver.default_parameters(), 1)
 
Hälsningar,
 
Josefin
 

## Fortsättning: 

Här kommer en liten förklaring av det jag lagt till, ni kan börja med
att testa att ändra på parametrar och få en känsla för vad som händer.
Det blev ett väldigt långt mejl igen, men jag har känslan att jag inte
förklarat så bra, så har ni några frågor tveka inte att mejla!


### Ändra på material-parametrar:
-------------------------------
För att ändra på parametrarna epsilon och p har jag lagt till dessa
rader i ismip_hom_a_2d.py:

	ps.pstokes_parameters['eps_reg']=1e-12
	ps.pstokes_parameters['p']=1.33

Se slidsen jag skickade innan jul för att kolla vad de är, i
spänningtensorn finns dessa parametrar. p bestämmer hur olinjärt
materialet är, teoretiskt sätt kan p ligga mellan 1 (maximalt olinjärt,
viskositeten blir väldigt hög i vissa punkter), och 2 (helt linjärt,
används för att modellera vanliga Newtonska fluider). För is sätter man
oftast p=4/3 men det är egentligen en stor diskussion ifall det är rätt.
Det kan hursomhelst vara intressant att titta hur p-värdet påverkar
konvergensen av den linjära lösaren och Newtons lösare.
Epsilon är en parameter som gör så att viskositeten inte blir oändlig.
Man tror att epsilon är noll eller vädligt litet i verkligheten, men i
praktiken blir det svårt att få lösarna att
fungera när epsilon är litet. Jag har sett Grönlands-simuleringar där
folk sätter epsilon till 1e-2..

### Ändra på parametrar rörande Newton-lösaren:
--------------------------------------------

I ismip_hom_a_2d.py finns även dessa rader, som sätter parametrar
gällande Newton-lösaren (eller de första två raderna gäller den linjära
lösaren):
```python
newton_parameters = Parameters('newton_solver')
newton_parameters.add('linear_solver', 'gmres')  # the linear solver
newton_parameters.add('preconditioner', 'ilu')   # the preconditioner
newton_parameters.add('maximum_iterations',  50)  #maximum number of
allowed iterations for Newton solver
newton_parameters.add('relaxation_parameter', 1.0)  # If 1: a normal
Newton solver. If between 0 and 1, the new solution will be a mix of the
old and new one
newton_parameters.add('relative_tolerance',  1e-7)  #Newton solver will
stop when the residual r divided by right hand side is smaller than this
number
```
### Ändra på parametrar rörande linjära-lösaren:
--------------------------------------------

i pstokes.py har jag lagt till dessa rader
```python
		solver = KrylovSolver(linear_solver, preconditioner)
     		solver.parameters['relative_tolerance']=1E-2	# iterations
stop if relative residual is smaller than this
        	solver.parameters['absolute_tolerance']=1E-10   # iterations
stop if residual is smaller than this
        	solver.parameters['maximum_iterations']=10000   # maximum
allowed iterations of linear solver
        	solver.parameters['monitor_convergence'] = True # just to see
what's going on
        	solver.parameters['report'] = True

        	#info(solver.parameters, 1)
```

Så att ni kan ändra på parametrar angående den linjära lösaren. För att
det sak funka löses systemet såhär (i Christians kod står det inte
solver. framför)
```python
        	solver.solve(A, U_inc.vector(), b)
		
```

### Kolla symmetri och konditionstal
--------------------------------------------

Jag har även lagt till dessa rader i pstokes.py för att kolla om
Jacobianen J är symmetrisk:

```python
Check symmetry
        	Jac = PETScMatrix()
        	assemble(J, Jac)
        	print('Is symmetric', np.linalg.norm(Jac.array() -
Jac.array().T) < 1E-10)
```

För att kolla konditionstalet finns denna rad:

	print("Condition number is:",np.linalg.cond(Jac.array()))

Dock ser jag att konditionstalet alltid ser väldigt stort ut lite
oavsett vad man sätter p och epsilon till, det är ett litet mysterium
just nu, jag ska dubbelkolla det.

Hälsningar,

Josefin

### Sammanfattning v.8
Hej,

Ursäkta mitt sena svar! Jag svarar "inline" i ert mejl:

>
>
> Det vi gjort det senaste är:
>
>
> 1.  Undersökt mer om konditionstalets relation till meshstorlek och
> kryloviterationer.


Här ser man ju iaf att konditionstal och krylov-iterationer ökar med
minskande meshstorlek. Två grejer att fundera på här är:

- Öker det som h^2 (detta är vad teorin säger)?
- Varför är det så stort även för grova nät? Detta kan vi diskutera
tillsammans på gruppmötet

>
> 2.  Hur konditionstalet beror på antalet Newton iterationer.
>

Det ser ut att vara någorlunda stabilt som ni säger, och alltid av
storleksordningen 1e18. Så det var ju skönt, det hade blivit rörigt att
tänka om det varierade mycket :)

> 3.  Hur konditionstalet beror på parametern p.
>

Detta var ett knivigt problem. Dock ser konditionstalet även här
relativt stabilt ut, dvs det är alltid av storleksordningen 1e18. Man
kan inte tydligt se att förhållandet som finns i ekvation 6.14 i Hirns
avhandling gäller. Detta kan vi diskutera på gruppmötet, jag ska fundera
lite.

Christian håller på att fixa så ni kan ta ut bara A-delen av matrisen,
man kan ju hoppas att det blir tydligare då.

> 4.  Hur konditionstalet beror på parametern epsilon.
>
Jag hittar inte dom här plottarna, visa dom hemskt gärna på gruppmötet!
Det vore intressant att se om ekvation 6.14 verkar gälla här.
>
>
>
> Sedan undrar jag också lite över relevansen av att plotta mot alla p i
> intervallet (1,2). I din artikel " Equal-Order Stabilized Finite Element
> Approximation of the p-Stokes Equations on Anisotropic Cartesian Meshes"
> så skriver du att p \approx 1,33 = 4/3 för is. Undrar därmed hur mkt
> fokus vi ska lägga på andra p:n än just detta, men det är ju något vi
> kan diskutera mer när vi ses.
>
Det är mest för att kunna kolla att det stämmer med Hirns formel 6.14,
så vi förstår om vår ismodell beter sig enligt teori och om inte, varför.
>
> Det som kommer vara vårt fokus framöver är väl egentligen att komma på
> ett sätt att extrahera A-matrisen, samt ta ut beloppet av
> töjningstensorn som också var en del i hans uppskattning, kanske med
> hjälp av er, då det verkade lite klurigt att göra utan erfarenhet av
> FEniCS. Hittade någon grej här men har inte grävt överdrivet mycket
> (https://bitbucket.org/fenics-apps/cbc.block/src/master/block/block_assemble.py).
> Ur denna hoppas vi kunna göra en jämförelse med Hirns resultat och inte
> långt efter förhoppningsvis kunna dra några slutsatser om de linjära
> lösarna baserad på all information vi framtagit.
>
Christian håller på!
>
> Dock lär vi ju nästa vecka börja kika mer på Hirns experiment  för
> Newton's metod, och jag personligen känner att en del tid där också
> kommer behöva läggas på att läsa in sig lite på teorin.
>
Ja om ni inte riktigt vet hur ni ska gå vidare innan vi diskuterar på
gruppmötet så lägg tiden på att läsa på på detta, och på att skriva på
rapporten.
>
>
>
En annan sak jag tänkte på är att det vore snyggt att ha plottarna med
residualen vs antal newton iterationer så att det är log på y-axeln men
inte x-axeln, men om det är tradigt att ändra så strunt i  det.

Bra jobbat, det är inte ett lätt problem ni har men jag tycker det är
sjukt intressant!
