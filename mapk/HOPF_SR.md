# Provera Hopfovih pragova, verzija 2

Napomena za paket verzije 3: ovaj fajl opisuje prethodnu fazu. Prvi
Ljapunovljev koeficijent je sada izračunat; videti LYAPUNOV_SR.md.

Pokrenite `hopf_analysis.py` u istom folderu u kome je `mapk_model.py`.
U Spyder-u otvorite hopf_analysis.py i pritisnite F5. Instalacija zavisnosti
je ista kao u README_SR.md. Prethodna skripta mapk_model.py ostaje dostupna.
Ova analiza ne pokreće ponovo ceo skup prethodnih simulacija.

## Rezultati

V1 je u nM/s, periodi u sekundama. Period u pragu računat je kao
T = 2*pi/omega, gde je kritični par sopstvenih vrednosti +/- i*omega.

| Scenario | Prag | V1 Python | V1 u radu | Period Python | Period u radu |
|---|---|---:|---:|---:|---:|
| Referentni, k11=k12=0 | H1 | 0.3842452864 | 0.3842 | 637.736630 | 637.7 |
| Referentni, k11=k12=0 | H2 | 6.2480479397 | 6.248 | 1227.011026 | 1227 |
| Procenjeni, k11=k12=0.0025 1/s | H1 | 0.3862273464 | nije navedeno | 635.347710 | nije navedeno |
| Procenjeni, k11=k12=0.0025 1/s | H2 | 6.1485359755 | nije navedeno | 1234.033107 | nije navedeno |

Referentne vrednosti pragova i perioda slažu se sa radom do objavljenog
zaokruživanja. Parametri nisu podešavani da bi se pragovi pogodili: korišćene
su iste pretpostavke kao u verziji 1 (fosfataze po 1 nM; dodatne brzine
isključene u referentnom scenariju).

Na praćenoj fizički dozvoljenoj ravnotežnoj grani, u pregledanom rasponu
0.1 <= V1 <= 9, pronađena su po dva prelaza stabilnosti za oba scenarija.
To nije iscrpna pretraga svih mogućih grana ravnoteže ili globalnih bifurkacija.

Za referentni scenario ravnoteža je lokalno asimptotski stabilna ispod H1,
nestabilna između H1 i H2 i ponovo lokalno asimptotski stabilna iznad H2.
Isto važi za procenjeni scenario, sa njegovim pomerenim pragovima.
Ovo se odnosi na petodimenzionalnu konzervacionu mnogostrukost.

## Numerički postupak i provere

1. Nastavak fizički dozvoljene ravnotežne grane uz prethodno rešenje kao
   početnu aproksimaciju. Rešavaju se ravnotežne jednačine, ne koristi se
   završno stanje vremenske simulacije kao zamena za nestabilnu ravnotežu.
2. Spektar 5D Jakobijana, bez tri nulte sopstvene vrednosti konzervacije.
3. Uokviravanje promene znaka najvećeg realnog dela i Brentova metoda.
4. Na svakom pragu postoji jedan prost kompleksni par sa realnim delom
   numerički nula i nenultom frekvencijom. Preostale tri sopstvene vrednosti
   imaju strogo negativne realne delove.
5. Prelaz je transverzalan: d(Re lambda)/dV1 nije nula. U referentnom
   slučaju nagibi su +0.00752067 i -0.00155457; u procenjenom slučaju
   +0.00737168 i -0.00145235. Nagibi su ponovo računati sa prepolovljenim
   korakom parametra.
6. Obrnuti nastavak na drugačijoj mreži reprodukuje pragove uz odstupanje
   manje od 4e-13 nM/s. To je unutrašnja numerička konzistentnost, a ne
   biološka preciznost parametara.
7. Kompleksni korak za Jakobijan nezavisno je upoređen sa centralnim
   razlikama na dve veličine koraka; odstupanja kritične sopstvene vrednosti
   manja su od 1e-12 1/s u ovom izračunavanju.
8. Simulacije malih perturbacija duž kritičnog moda pri V1=0.98*H i
   V1=1.02*H potvrđuju očekivani rast/opadanje. Amplituda je projektovana
   na kompleksni mod pomoću inverzne matrice sopstvenih vektora.
   Provera traje tri lokalna perioda i testira lokalnu linearnu dinamiku.

## Granica zaključka

Ovo su numerički potvrđeni linearni Hopfovi uslovi i pragovi gubitka/povratka
stabilnosti. Prvi Ljapunovljev koeficijent nije računat, pa ovom verzijom
ne određujemo superkritičnost/subkritičnost, niti isključujemo degenerisani
Hopfov slučaj. Nisu računati Floquetovi multiplikatori ili nastavak
periodičnih orbita. Za punu lokalnu klasifikaciju sledeći korak je prvi
Ljapunovljev koeficijent; za granu ciklusa kontinuacija periodičnih rešenja.

Ne treba izjednačavati period 2*pi/omega u pragu sa periodom ciklusa pri
V1=5. To su različite tačke u prostoru parametara.

## Izlazni fajlovi

U results_hopf:

- hopf_thresholds.png: ravnotežna grana KPP i najveći realni deo spektra;
- hopf_thresholds.csv: pragovi, frekvencije, periodi i nagibi;
- reference_equilibrium_branch.csv i estimated_equilibrium_branch.csv:
  sve ravnotežne koncentracije, stabilnost i reziduali duž grane;
- hopf_summary.json: svi parametri, kritični spektri i numeričke provere.

Grafik prikazuje samo ravnoteže. Isprekidana crvena linija predstavlja
nestabilne ravnoteže, ne maksimum ili minimum oscilacija.

Izvor poređenja: dostavljeni rad N. A. Kablar, MTNS 2010, str. 2367.
