# Osetljivost MAPK modela: verzija 4

## Pokretanje

U Spyder-u otvorite sensitivity_analysis.py i pokrenite F5. Skripta koristi iste
zavisnosti i ne menja matematičko jezgro. Rezultati su već uključeni u paket.
Ponovno pokretanje zamenjuje izlaze istih naziva u results_sensitivity.
Broj scenarija i adaptivno produžavanje simulacija znače da ova skripta može
trajati duže od prethodnih. Napredak se ispisuje u konzoli.

## Opseg i značenje analize

Radna tačka: k11=k12=0.0025 1/s, obe fosfataze po 1 nM.
Ostali parametri su isti kao u prethodnoj verziji. Opsezi su scenariji:
nisu eksperimentalni intervali poverenja ili potvrđeni fiziološki rasponi.

- k11,k12: 0, 0.00125, 0.0025, 0.005, 0.01 1/s.
- Fosfataze drugog (phosphatase_KK) i trećeg (phosphatase_K) nivoa:
  0.5, 0.75, 1, 1.25, 1.5 nM, nezavisno.
- Zajednička mreža k11 x k12: 25 kombinacija; fosfataze su po 1 nM.
- Ukupno 34 imenovana scenarija, uključujući referentni; referentni je
  numerički isti kao zajednička kombinacija (0,0), pa ima 33 različita skupa.
- Ravnoteže i linearni Hopfovi prelazi: V1 od 0.02 do 20 nM/s.
- Vremenski odzivi: V1=0.25 i V1=5 nM/s, isti početni uslovi kao ranije.
- Lokalna osetljivost pragova i perioda u pragu: koraci +/-1% i +/-0.5%.
- Lokalna osetljivost KPP ravnoteže: implicitno diferenciranje 5D sistema,
  nezavisno provereno centralnim razlikama. Pri V1=5 to je nestabilna
  ravnoteža; njena koncentracija nije srednja vrednost oscilatornog odziva.

## Glavni nalazi

Oko radne tačke k11 i k12 deluju suprotno na H1: povećanje k11 podiže,
a povećanje k12 spušta prag. Na H2 oba lokalno deluju u smeru snižavanja.
Zajednička promena zato može prikriti različite pojedinačne efekte.
Nelinearni zajednički efekat izračunat je direktno u svakoj tački mreže,
bez pretpostavke da se pojedinačni efekti prosto sabiraju.

Fosfataze imaju veću lokalnu relativnu osetljivost pragova od dodatnih
katalitičkih konstanti u ovoj radnoj tački. Promene nisu globalno monotone:
na primer H1 za phosphatase_KK prvo opada pa raste u ispitanim tačkama.
To nije opšte pravilo o povećanju/smanjenju fosfataza u ćeliji.

## Lokalna normalizovana osetljivost pragova

S=(p/H)*dH/dp. Na primer S=0.295 znači da lokalno povećanje parametra
za 1% pomera prag približno +0.295%. Prikazan je korak +/-0.5%.

| Parametar | S za H1 | S za H2 |
|---|---:|---:|
| k11 | +0.10180 | -0.00604 |
| k12 | -0.08975 | -0.00898 |
| phosphatase_KK | +0.29464 | +0.33149 |
| phosphatase_K | -0.17048 | +0.08674 |

## Izabrani rezultati pri V1=5

A=(KPP_max-KPP_min)/2 je poluamplituda, ne razlika maksimum-minimum.
Period je iz ustaljenih vrhova odziva, a ne 2*pi/omega u Hopfovoj tački.

| Scenario | H1 | H2 | A (nM) | T (s) |
|---|---:|---:|---:|---:|
| baseline | 0.386227 | 6.148536 | 79.149 | 1291.04 |
| k11_0.01 | 0.516436 | 6.111321 | 77.152 | 1267.79 |
| k12_0.01 | 0.355038 | 6.014631 | 68.928 | 1249.19 |
| phosphatase_KK_0.5 | 0.401542 | 4.333015 | 0.000 | nema oscilacija |
| phosphatase_KK_1.5 | 0.469900 | 6.854816 | 96.286 | 1288.93 |
| phosphatase_K_0.5 | 0.515152 | 4.148622 | 0.000 | nema oscilacija |
| phosphatase_K_1.5 | 0.380041 | 6.298625 | 97.052 | 1383.89 |
| joint_0.01_0.01 | 0.423156 | 5.831611 | 64.032 | 1192.63 |

U svim ispitanim scenarijima pri V1=0.25 odziv je dostigao ravnotežu.
Pri V1=5 smanjenje pojedinačne fosfataze na 0.5 nM pomera H2 ispod 5,
pa radna vrednost V1 izlazi iz intervala nestabilne ravnoteže. Time se
objašnjava prestanak oscilacija bez potrebe da opšti efekat fosfataze
proglasimo stabilizujućim ili destabilizujućim za sve uslove.

## Provere i ograničenja

- Pragovi: potvrđen imaginarni par, transverzalni prelaz i tri preostala
  stabilna moda. Dva prelaza nađena su u svakom scenariju na praćenoj grani.
- Periodičnost: najmanje sedam vrhova i dolina u kasnoj polovini odziva;
  promena perioda i visine ekstrema mora biti manja od 0.5%.
- Simulacija počinje sa 24000 s i po potrebi se produžava do 96000 s.
  Svi priloženi scenariji su klasifikovani; nema unresolved slučajeva.
- Ravnoteža u simulaciji: odstupanje svih stanja <1e-4 nM i brzine <1e-6 nM/s.
- Provereni su nenegativnost koncentracija i tri očuvana ukupna zbira.
- Na šest reprezentativnih slučajeva ponovljena je pretraga pragova
  dvostruko gušćom mrežom i odziv pri V1=5 nezavisnim Radau rešavačem.
- Lokalni izvodi pragova/perioda provereni su na dva koraka parametra.
- Skeniranje mreže ne isključuje vrlo uske dodatne intervale ili druge grane.
- Ustaljene oscilacije u konačnoj simulaciji nisu Floquetov dokaz stabilnosti.
  Prvi Ljapunovljev koeficijent nije ponovo računat za svih 34 scenarija;
  prethodna superkritična klasifikacija ne prenosi se automatski na njih.
- Nije rađena zajednička mreža obe fosfataze niti globalna statistička
  analiza (npr. Sobol); ovo su OAT, lokalna i mrežna analiza k11,k12.

## Fajlovi

- sensitivity_analysis.py: cela analiza i generisanje ovog pregleda.
- results_sensitivity/scenarios.json: svi parametri i metrički rezultati.
- threshold_scenarios.csv: pragovi i periodi u pragovima.
- dynamic_scenarios.csv: ravnoteže, minimumi, maksimumi, poluamplitude, periodi.
- local_sensitivities.csv: izvodi i normalizovane osetljivosti H1,H2 i njihovih perioda.
- local_equilibrium_sensitivities.csv: lokalna osetljivost KPP ravnoteže.
- validation.json: detalji provera gušće mreže i rešavača.
- one_parameter_sensitivity.png: pragovi, poluamplitude i periodi po parametru.
- joint_k11_k12_thresholds.png: pragovi po zajedničkoj mreži; ćelije predstavljaju
  diskretne scenarije, ne jednako razmaknutu numeričku osu ili interpolaciju.
- phosphatase_time_responses.png i *_V5_trajectory.csv: tri primera odziva.

Matematički model i njegove prethodno navedene pretpostavke ostaju isti.
Ovaj paket je osnova za sledeću fazu: Streamlit interfejs.
