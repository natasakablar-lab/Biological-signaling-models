# IL1 NFκB Streamlit v3 — ravnoteže i stabilnost

Zadržani su puni ODE (22), redukovani ODE (16), DAE (16+6) i poređenje svih
22 koncentracije. Dodata je kartica „Ravnoteže i bifurkacije“.

## Pokretanje

Raspakovati ceo ZIP i iz direktorijuma aplikacije pokrenuti:

    python -m pip install -r requirements.txt
    python -m streamlit run app.py

Za DAE i poređenje sva tri zapisa:

    python -m pip install -r requirements-dae.txt
    python -m streamlit run app.py

Ako je prethodna verzija već radila, obično je dovoljno pokrenuti novu app.py
iz novog direktorijuma sa istim Python okruženjem. START_WINDOWS.bat služi
pokretanju nakon instalacije.

## Nova analiza

1. U bočnom meniju postaviti parametre i početne uslove, pa izračunati odzive.
2. Otvoriti „Ravnoteže i bifurkacije“ i kliknuti „Izračunaj ravnotežu“.
3. Pregledati svih 22 ravnotežnih koncentracija i 16 svojstvenih vrednosti.
4. Izabrati parametar i opseg faktora, pa „Prati ravnotežnu granu“.
5. Birati koncentraciju za dijagram; preuzeti ravnotežu i granu kao CSV.

Analiza uvek koristi parametre i šest konstanti poslednjeg uspešnog proračuna.
Rezultati prethodne postavke ne prikazuju se kao aktuelni. Posle promene
postavke ponoviti računanje. Poređenje tri zapisa i dalje se pokreće samim
uključivanjem opcije, kao u verziji 2.1.

Metoda ravnoteže trenutno zahteva sve stope i svih šest konstanti strogo
pozitivne. Granični slučajevi poput ukupnog IL1=0 ili k4=0 nisu izolovane
ravnoteže iste strukture i jasno se prijavljuju kao posebna analiza.

## Nalaz za početni model

Fizička ravnoteža je lokalno asimptotski stabilna. Zbog razgradnje i odsustva
obnavljanja IKK, sve IKK vrste imaju nultu ravnotežnu koncentraciju.

Devet pretraga (k11,k7,k8,k6,a10,k13,a11,k19,k4), sa faktorima 0.01–100 i
81 tačkom po parametru, ne daje Hopf kandidata na pronađenoj grani.
To nije globalni dokaz odsustva bifurkacija i ne uvodi izmene biologije.

## Fajlovi i ponavljanje

- NALAZI_RAVNOTEZE.md: detaljan nalaz, izvođenje ravnoteže, dokaz granice IKK,
  svi brojevi, opsezi i ograničenja pretrage.
- equilibria.py: skalarna konstrukcija fizičke ravnoteže, 16-D Jakobijan,
  praćenje parametarske grane i lokalizacija potencijalnih Hopf kandidata.
- FORMULACIJE.md: puni, redukovani i DAE zapis iz prethodne verzije.
- results_equilibria/: CSV, JSON, NPZ spektri i devet grafika.

    python run_equilibrium_analysis.py
    python verify_equilibria.py
    python verify_equilibrium_app.py

Pretraga je praćenje po prirodnom parametru sa skalarnim rešavanjem ravnoteže,
a ne pseudoarclength nastavak. Skriveni/tangencijalni koreni i uske grane mogu
zahtevati posebne metode. Svaka tačka proverava rezidual svih 22 ODE,
šest bilansa i fizičku prihvatljivost. Prikaz Hopf kandidata nije automatska
potvrda nedegenerisane Hopf bifurkacije ili postojanja stabilnog ciklusa.

Kinetički parametri su i dalje predmet provere literature i eksperimentalnih
uslova. Autorka je potvrdila izvođenje modela, jedinice k6/k19 i aktivni x1.
