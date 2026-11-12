# MAPK Research Lab: Streamlit aplikacija, verzija 5

## Pokretanje na računaru

Raspakujte ceo ZIP u jedan folder. Sačuvajte zajedno app.py, matematičke
module i foldere sa rezultatima. U Anaconda Prompt-u ili terminalu pređite
u folder MAPK_Python_v5, pa pokrenite:

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Aplikacija se otvara u browseru, obično na http://localhost:8501.
Ovo je lokalna aplikacija: nije potrebno objavljivanje na internetu.
Za zaustavljanje pritisnite Ctrl+C u terminalu. Pri sledećem pokretanju
nije potrebno ponavljati instalaciju. Na Windows-u možete pokrenuti i
run_app.bat, posle instalacije zavisnosti u aktivnom Python okruženju.

Streamlit aplikacija se pokreće iz terminala navedenom komandom;
F5 u Spyder-u ostaje način pokretanja prethodnih analitičkih skripti.
Možete i dalje koristiti sve skripte iz ranijih verzija bez interfejsa.

## Šta je uključeno

1. **Odzivi:** izbor početnog scenarija, V1, k11,k12 i koncentracija obe
   fosfataze; interaktivni vremenski i fazni grafici svih osam koncentracija.
   Prikaz celog ili samo kasnog dela odziva i preuzimanje PNG grafika.
2. **Ravnoteža i stabilnost:** ravnotežne koncentracije, sopstvene vrednosti
   i Jakobijan u pet nezavisnih koordinata; konzervacija i rezidual.
3. **Hopf i Ljapunov:** na zahtev se računaju pragovi i ravnotežna grana
   za primenjene parametre, a zatim, posebnim dugmetom, prvi Ljapunovljev
   koeficijent i lokalna klasifikacija. Prikazani su smer i stabilnost
   malih ciklusa u odnosu na V1.
4. **Osetljivost:** pregled svih prethodno proverenih scenarija verzije 4,
   pripadajućih pragova i simulacionih metrika, uz grafike i CSV tabele.
   Ovaj tab prikazuje sačuvanu studiju: nije nova analiza trenutnih parametara.
5. **Model i preuzimanje:** opis pretpostavki, svi primenjeni parametri,
   vremenski podaci kao CSV, parametri i sažetak kao JSON i Jakobijan kao CSV.

## Parametri i rad sa rezultatima

- Izaberi početni scenario ili ručno promeni parametre u bočnoj traci.
- Pritisni **Primeni i simuliraj**. Promena unapred pripremljenog scenarija
  samo popunjava polja; rezultati se menjaju tek posle primene.
- Plava poruka iznad tabova uvek prikazuje parametre izračunatog rezultata.
  Tako ne treba pomešati rezultate sa još neprimenjenim unosom.
- Posle primene novih parametara brišu se stari Hopfovi i Ljapunovljevi
  rezultati. Izračunaj ih ponovo na odgovarajućem tabu.
- Promena gornje granice pretrage zahteva ponovni klik na računanje pragova.
- Preuzete datoteke se snimaju kroz browser na tvoju lokaciju za preuzimanje.

Početno se prikazuje oscilatorna radna tačka V1=5, k11=k12=0.0025 1/s,
fosfataze po 1 nM. Dostupni su i referentni model, ravnoteža pri V1=0.25
i primeri smanjenih fosfataza.

U interfejsu su omogućeni V1=0.02-20 nM/s, k11,k12=0-0.025 1/s i
fosfataze=0.25-2 nM. To su istraživačke granice unosa, ne eksperimentalno
potvrđeni fiziološki opsezi. Nije svaki unos unapred analiziran; ako
numerički postupak ne uspe, prikazuje se poruka bez izmišljanja rezultata.

## Matematički status

Matematičko jezgro mapk_model.py nije menjano. Aplikacija koristi postojeće
funkcije za dinamiku, ravnotežu, pragove i koeficijente. Hillov koeficijent
ostaje n=1, a Va,Vd,kp su nula. Ukupne količine su 100,300,300 nM i
početni uslovi su isti kao u prethodnim verzijama. Sve detaljne pretpostavke
ostaju navedene u README_SR.md i ostalim analitičkim pregledima.

Ustaljivanje se proverava istim kriterijumima kao u analizi osetljivosti:
simulacija se po potrebi produžava za još tri bloka po 24000 s, nakon
izabranog početnog trajanja. Ako kriterijumi nisu ispunjeni, prikazuje se
da rezultat nije potvrđen; period i amplituda se ne predstavljaju kao
ustaljeni. Duži prelazni proces nije dokaz periodičnog rešenja.

Hopfova pretraga prati jednu ravnotežnu granu od V1=0.02 do izabrane
gornje granice. Veoma uski intervali i druge grane mogu biti propušteni.
Prvi Ljapunovljev koeficijent koristi Kuznetsovljevu konvenciju i ||q||=1
u koordinatama nM. Koeficijenti se sada mogu računati i za nov primenjeni
skup parametara. Za skoro nulti koeficijent potrebna je dodatna provera;
pri |l1|<1e-9 interfejs traži dodatnu proveru umesto klasifikacije. Ovaj prag je praktična numerička zaštita, ne matematički dokaz degeneracije. Globalna grana
periodičnih orbita i Floquetovi multiplikatori nisu uključeni.

Ovaj MAPK model još nije model Alzheimerove bolesti. Vezu sa tom bolešću,
potrebne module, ulaze i izlaze, kao i validaciju treba odrediti zasebnim
pregledom literature i eksperimentalnih podataka. U aplikaciji nisu dodate
pretpostavljene bolesti, terapijski efekti ili biološki zaključci.

## Provere interfejsa

test_streamlit_app.py proverava početni oscilatorni scenario, računanje
Hopfovih pragova i negativnih Ljapunovljevih koeficijenata, promenu scenarija,
prelazak na ravnotežu pri smanjenoj fosfatazi, uklanjanje zastarelih
rezultata i referentnu ravnotežu iz rada. Pokretanje:

```bash
python test_streamlit_app.py
```

Glavni tokovi provereni su preko Streamlit AppTest-a. Grafički prikaz i klik na preuzimanje u stvarnom browseru nisu provereni u ovom okruženju jer pokretanje browsera nije dostupno. Testirano sa Streamlit 1.66.1 i Plotly 7.1.0. Minimalne verzije iz
requirements.txt nisu zasebno testirane. Playwright nije potreban za
pokretanje aplikacije i nije dodat korisničkim zavisnostima.

Za narednu fazu korisno je prvo dovršiti rekonstrukcije NFκB i Wnt,
zatim proveriti literaturu i izabrati zajedničke merljive veličine sa
eksperimentalnim saradnicima. Žurnalsku novinu treba proceniti tek nakon
te provere, odvojeno od uspešne reprodukcije starog modela.
