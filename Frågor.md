### pstokes
* Varför en faktor 0.5 i uttrycket för strain-rate tensorn(se uttrycket för mu_eff i pstokes.py)? Borde vara 1.
* Hur kan man ändra startgissningen i fenics? Det borde minska antalet newtoniterationer innan konvergens.
* Hur ska undersökningen av konditionstalet för arolla gå till?
* Varför konvergerar ISMIP-HOM så mycket snabbare med SI_units=True?
* Hur kan Hirn få konvergens i sina experiment för eps_0 = 0? Anledningen till att vi har epsilon i Carreu-modellen är väl för att undvika singulariteter på grund av den icke-linjära viskositeten?
* Funkar inte att beräkna konditionstal för arolla genom singulärvärden
* Hur pass försiktiga ska vi vara med att ändra "konstanterna"(normen av DV, c) i Hirns övre begränsning av konditionstalet? Borde dessa vara samma när vi undersöker epsilon och p, eller är det rimligt att skruva på värdena för att få bättre resultat?
