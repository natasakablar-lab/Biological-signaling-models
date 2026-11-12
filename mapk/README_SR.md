# MAPK: rekonstrukcija u Pythonu, verzija 1

Izvor: Nataša A. Kablar, "MAPK Module: Biological Basis, Structure,
Mathematical Model and Dynamical Analyse", MTNS 2010, str. 2363-2368.
Rekonstrukcija na osnovu dostavljenog PDF-a; originalni MATLAB fajlovi nisu korišćeni.

## Pokretanje

Raspakujte ZIP u sopstveni folder. Python 3.10 ili noviji:

```bash
python -m pip install -r requirements.txt
python mapk_model.py
```

U Spyder-u otvorite mapk_model.py i pokrenite ceo fajl (F5), uz instalirane
pakete u istom Python okruženju. Kod je nezavisan od MATLAB-a i Streamlit-a.
Slike i numerički rezultati snimaju se u podfolder results pored skripte.
Ponovno pokretanje zamenjuje rezultate istih naziva. Grafici se čuvaju kao PNG,
bez potrebe za grafičkim interfejsom. Priloženi rezultati su već izračunati.

## Model i pretpostavke

Implementirane su diferencijalne jednačine (1)-(8) i izrazi za brzine sa
str. 2365. Oznake x u izrazima za brzine i y u ODE sistemu tumače se kao
iste koncentracije. Stanja su redom:
KKK, KKKP, KK, KKP, KKPP, K, KP, KPP. Vreme je u sekundama,
koncentracije u nM. Na grafiku je vreme prikazano u minutima.

| Parametar | Vrednost | Poreklo / status |
|---|---:|---|
| V1 | 0.25 ili 5 nM/s | Osnovni i oscilatorni scenario iz teksta rada |
| V2 | 0.25 nM/s | Tabela 3 |
| n, KI | 1; 9 nM | Tabela 3 |
| k3,k4,k7,k8 | 0.025 1/s | Tabela 3 |
| k5,k6 | 0.75 1/s | Tabela 3 |
| k9,k10 | 0.5 1/s | Tabela 3 |
| KM1,KM2 | 10; 8 nM | Tabela 3 |
| KM3 do KM10 | 15 nM | Tabela 3 |
| Ukupna fosfataza KK i K | po 1 nM | Pretpostavka koja reprodukuje navedenu osnovnu ravnotežu |
| k11,k12: referentni scenario | 0 | Dodatne reakcije isključene |
| k11,k12: procenjeni scenario | po 0.0025 1/s | Radna pretpostavka: 10% k7,k8; nije eksperimentalno određeno |
| KM11,KM12 | po 15 nM | Pretpostavka po analogiji sa KM7,KM8 |
| Va,Vd,kp | 0 | Cross-talk i pozitivna sprega isključeni u ovoj verziji |
| KMa,KMd | po 10 nM | Tabela 3; neaktivni pri Va=Vd=0 |
| KMp | 10 nM | Neaktivna pretpostavka pri kp=0 |

Va,Vd u kodu su efektivne maksimalne brzine, u nM/s. To nisu katalitičke
konstante ka1,kd1 iz tabele, u 1/s: za prelaz bi bila potrebna i odgovarajuća
koncentracija aktivatora/deaktivatora. Obe brzine su zasad nula.

Početno: KKK=100, KK=300, K=300 nM; ostala stanja su nula.

Za procenjeni scenario koristimo:

v11 = k11 * KKP * K / (KM11 + K)

v12 = k12 * KKP * KP / (KM12 + KP)

Ove brzine predstavljaju KKP kao katalizator reakcija K -> KP i KP -> KPP;
KKP se u tim reakcijama ne troši, u skladu sa objavljenim ODE sistemom.
U tekstualnim hemijskim šemama rada postoje proizvodi koji sugerišu promenu
stanja katalizatora; takvu drugačiju stehiometriju nismo uvodili u ODE.
ATP/ADP nisu eksplicitna stanja; model implicitno pretpostavlja raspoloživ
energetski supstrat i nije zatvoren termodinamički model ćelije.

## Šta znači fizička prihvatljivost ove verzije

Za nenegativne koncentracije i pozitivne konstante, brzine su nenegativne,
imaju odgovarajuće dimenzije i saturaciju po supstratu. Na granici gde je
neka koncentracija nula, njena brzina promene ne vodi van nenegativne oblasti.
Tri zbira su očuvana:

- KKK + KKKP = 100 nM;
- KK + KKP + KKPP = 300 nM;
- K + KP + KPP = 300 nM.

To je osobina objavljenog osmostanog efektivnog ODE modela. Kompleksi iz
Tabele 1 nisu posebna stanja: ne tvrdimo da je suma slobodnih oblika tačna
ukupna koncentracija u eksplicitnom modelu sa kompleksima. Za taj model
bile bi potrebne dodatne reakcije i konstante.

Nenegativnost i konzervacija omogućavaju fizički prihvatljivu simulaciju,
ali same po sebi ne dokazuju eksperimentalnu realističnost izabranih k11,k12.
Njihova vrednost je ilustrativna, a ne statistička procena ili fit podataka.
Kao narednu proveru ima smisla varirati ih nezavisno, na primer
0, 0.00025, 0.0025 i 0.025 1/s, kao scenarije, bez tvrdnje da je to potvrđen
biološki opseg. KKP može biti mnogo zastupljeniji od KKPP, pa desetina
katalitičke konstante ne znači nužno desetinu ukupnog protoka.

## Redukcija, ravnoteža i stabilnost

Pet nezavisnih stanja su z=(KKK,KK,KKPP,K,KPP).
Ostala se računaju kao KKKP=100-KKK, KKP=300-KK-KKPP i KP=300-K-KPP.
Jakobijan se računa kompleksnim korakom u ovih pet koordinata; time se
izbegavaju tri nulte sopstvene vrednosti nastale konzervacijom u 8D modelu.
Jedinica Jakobijana i sopstvenih vrednosti je 1/s.

Za osnovni scenario, V1=0.25, k11=k12=0:

| Koncentracija | U radu (nM) | Python (nM) |
|---|---:|---:|
| KKK | 81.938 | 81.937969 |
| KK | 272.22 | 272.215222 |
| KKPP | 7.8503 | 7.850335 |
| K | 288.5 | 288.498325 |
| KPP | 2.5738 | 2.573760 |

Vrednosti se slažu do štampanog zaokruživanja. To je snažna provera
rekonstrukcije osnovnog scenarija, ali ne dokazuje jedinstvenost parametara.
Osnovni spektar 5D Jakobijana:

- -0.00211770 +/- 0.00882047 i;
- -0.01417130;
- -0.02856631;
- -0.03515554.

Ravnoteža je lokalno asimptotski stabilna na odgovarajućoj konzervacionoj
mnogostrukosti, sa prigušenim oscilatornim modom. Kompleksni deo daje
lokalni period oko 712.3 s, ne period održavanog ciklusa.

| Scenario | Najveći realni deo spektra (1/s) | Zaključak o ravnoteži | Odziv |
|---|---:|---|---|
| Referentni V1=0.25 | -0.00211770 | Lokalno asimptotski stabilna | Prigušivanje |
| Referentni V1=5 | +0.00096343 | Nestabilna | Održavane oscilacije u simulaciji |
| Procenjeni V1=0.25 | -0.00211124 | Lokalno asimptotski stabilna | Prigušivanje |
| Procenjeni V1=5 | +0.00087529 | Nestabilna | Održavane oscilacije u simulaciji |

Procenjeni slučaj pri V1=0.25 daje KPP=2.836057 nM.
Za V1=5 kasni period vrhova KPP je približno 1312 s u referentnom i 1291 s
u procenjenom scenariju. To su deskriptivne procene pri rezoluciji uzorkovanja
2 s i nisu periodi baš u Hopfovim tačkama. Kontinuacija i Floquet analiza nisu
urađene: održavane oscilacije u konačnoj simulaciji nisu dokaz stabilnosti
graničnog ciklusa. Nestabilna ravnoteža ne treba da se naziva marginalno
stabilnom samo zato što se pojavljuje periodično kretanje.

## Provere i fajlovi

Skripta proverava sve četiri konfiguracije:

- uspeh integracije i nenegativnost bez odsecanja koncentracija;
- očuvanje tri ukupne količine (najveće odstupanje oko 1.1e-12 nM);
- slaganje DOP853 i Radau (najveća razlika manja od 1.4e-6 nM);
- slaganje 8D i 5D integracije;
- rezidual pronađene fizički dozvoljene ravnoteže;
- Jakobijan poređen sa centralnim konačnim razlikama;
- tri približno nulte sopstvene vrednosti punog 8D Jakobijana.

Detaljne vrednosti, svi parametri i sopstvene vrednosti su u results/summary.json.
results/equilibria.csv sadrži sve ravnoteže; *_trajectory.csv sve vremenske
odzive; *_jacobian_5D.csv matrice u navedenom redosledu pet stanja.
time_responses.png i late_phase_portraits.png su grafici.

Testirano sa Python paketima NumPy 2.3.5 i SciPy 1.17.0 u ovom okruženju.
Minimalne verzije u requirements.txt nisu zasebno testirane.

## Sledeća faza

Numerička kontinuacija ravnoteže po V1, provera Hopfovih pragova iz rada,
zatim nezavisna osetljivost na k11,k12 i fosfataze. Posle provere modela
može se dodati Streamlit interfejs bez promene matematičkog jezgra.
