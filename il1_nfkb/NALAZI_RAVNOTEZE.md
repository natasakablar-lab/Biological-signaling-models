# IL1 NFκB — ravnoteže i prva bifurkaciona analiza

## Glavni nalaz

Podrazumevani potvrđeni model ima fizički prihvatljivu, lokalno asimptotski
stabilnu ravnotežu na fiksnoj klasi šest bilansnih konstanti. U devet
jednoparametarskih pretraga pronađena grana ostaje stabilna na ispitanim
mrežama. Nije detektovan Hopf kandidat niti prelazak svojstvene vrednosti
kroz nulu. To nije globalni dokaz odsustva bifurkacija ili drugih grana.

Nijedna reakcija ni vrednost osnovnog modela nije promenjena radi dobijanja
bifurkacije. Svi scan opsezi služe matematičkom istraživanju, bez tvrdnje da
su eksperimentalno potvrđeni. Poreklo kinetičkih parametara ostaje otvoreno.

## Bilans IKK određuje ravnotežnu granicu

Za ukupni IKK važi:

$$T_{IKK}=x_1+x_3+x_5+x_7+x_{21}+x_{22},\qquad \dot T_{IKK}=-k_4x_1.$$

Ako je k4>0, u nenegativnoj ravnoteži mora biti x1=0. Tada jednačina za x1
ima samo nenegativne članove:

$$0=(d_1+k_5)x_3+k_2x_5+d_4x_7+k_{19}x_{22}.$$

Za pozitivne stope sledi x3=x5=x7=x22=0. Iz jednačine za x22 zatim sledi:

$$a_{18}x_{21}x_{20}=0.$$

Pozitivne bilansne konstante IL1, ILR, MyD88, IRAK i TRAF6 i reverzibilna
vezivanja daju x20>0 u ravnoteži. Zato je i x21=0. Dakle:

$$x_1^*=x_3^*=x_5^*=x_7^*=x_{21}^*=x_{22}^*=0.$$

To je granica pozitivnog ortanta za IKK vrste, uz nenulte NFκB i IL1
koncentracije. Vrednosti k4 i k19 ne menjaju ovu ravnotežu dok su pozitivne,
ali menjaju linearizaciju i prelazne odzive.

Za periodičnu orbitu perioda P:

$$0=T_{IKK}(P)-T_{IKK}(0)=-k_4\int_0^P x_1(t)\,dt.$$

Nenegativnost i kontinuitet daju x1(t)=0 duž cele periodične orbite. Jednačina
za x1 potom daje i x3=x5=x7=x22=0. Trajno aktivan IKK zato nije moguć u
samoodrživom periodičnom režimu zatvorenog modela. Ovaj argument ne isključuje
svaki mogući ciklus NFκB podsistema na grani sa ugašenim IKK.

## Konstrukcija ravnoteže bez dugog vremenskog čekanja

Neka je n=x10 i b=x2. Definišemo:

$$L(n)=\frac{a_{10}a_{12}n}{(d_{12}+k_{13})d_{10}+k_{13}a_{12}n},$$
$$A_2(n)=\frac{a_3k_{13}L(n)}{a_{11}(d_3+k_3)},\qquad
A_1(n)=\frac{a_3d_{11}n/a_{11}+k_{13}L(n)}{d_3+k_3},$$
$$U(n)=\frac{k_8}{k_7}(k_6+k_{11}n^2).$$

Ravnotežni bilans IκB daje kvadratnu jednačinu:

$$k_3A_2b^2+(k_9+k_3A_1)b-U=0.$$

Za pozitivan koren koristimo numerički stabilan oblik:

$$b=\frac{2U}{k_9+k_3A_1+\sqrt{(k_9+k_3A_1)^2+4k_3A_2U}}.$$

Zatim:

$$x_{11}=L(n)b,\quad x_4=\frac{d_{11}n+k_{13}x_{11}}{a_{11}},\quad
x_6=A_2b^2+A_1b,$$
$$x_8=\frac{k_6+k_{11}n^2}{k_7},\qquad
x_9=\frac{a_{10}b}{d_{10}+k_{13}a_{12}n/(d_{12}+k_{13})}.$$

Preostaje jedna skalarna jednačina NFκB bilansa:

$$F(n)=x_4(n)+x_6(n)+n+x_{11}(n)-c_1=0,\qquad 0<n<c_1.$$

Važi F(0)=−c1<0 i F(c1)>0. Postoji bar jedan koren za pozitivne stope i c1.
Program proverava promene znaka na kombinovanoj linearnoj i logaritamskoj
mreži i koristi Brentov metod. Detektovane grane i njihovi reziduali se
proveravaju. Broj pronađenih korena na mreži nije dokaz globalne jedinstvenosti.

IL1 vezivanja rešavaju se nezavisno pomoću četiri ravnotežna odnosa:

$$x_{14}=\frac{a_{14}}{d_{14}}x_{12}x_{13},\quad
x_{16}=\frac{a_{15}}{d_{15}}x_{15}x_{14},$$
$$x_{18}=\frac{a_{16}}{d_{16}}x_{17}x_{16},\quad
x_{20}=\frac{a_{17}}{d_{17}}x_{19}x_{18},$$

uz pet fiksnih bilansa. Slobodne koncentracije računaju se u logaritamskim
koordinatama radi njihove pozitivnosti.

## Podrazumevana ravnoteža

Šest konstanti je (0.1,0.15,0.1,0.05,0.05,0.05) μM, istim redosledom kao
u FORMULACIJE.md. Sve ravnotežne koncentracije su ispod.

| Stanje | Vrsta | Ravnoteža μM |
|---|---|---:|
| x1 | IKK_P | 0 |
| x2 | IkB | 0.1794667796 |
| x3 | IKK_P:IkB | 0 |
| x4 | NFkB | 0.0002835956694 |
| x5 | IKK_P:IkB:NFkB | 0 |
| x6 | IkB:NFkB | 0.09749624571 |
| x7 | IKK_P:IkB:NFkB (x7) | 0 |
| x8 | IkBt | 0.005508189123 |
| x9 | IkBn | 0.1417312312 |
| x10 | NFkBn | 0.0003727827009 |
| x11 | IkBn:NFkBn | 0.001847375915 |
| x12 | IL-1 | 0.05096312499 |
| x13 | IL-R | 0.0009631249918 |
| x14 | C14 (IL-1:IL-R) | 0.04908385934 |
| x15 | MyD88 | 4.698432916e-05 |
| x16 | C16 (MyD88:C14) | 0.002306172204 |
| x17 | IRAK | 0.002353156533 |
| x18 | C18 (IRAK:C16) | 0.005426784187 |
| x19 | TRAF6 | 0.00777994072 |
| x20 | C20 (TRAF6:C18) | 0.04222005928 |
| x21 | IKK | 0 |
| x22 | C22 (IKK:C20) | 0 |

Najveći realni deo svojstvene vrednosti redukovanog Jakobijana je
**-7.3796454140e-06 s⁻¹**. Svih 16 svojstvenih vrednosti je realno
i negativno u ovoj početnoj parametarskoj tački. Lokalna ravnoteža je
asimptotski stabilna. Najsporija linearna vremenska skala je približno
135508 s (37.64 h).

Rezidual svih 22 ravnotežnih jednačina iznosi 1.016e-20 μM/s,
a odstupanje od konstantnih bilansa 9.714e-16 μM.
Udeo indukovane sinteze k11·n² u ukupnoj ravnotežnoj sintezi transkripta je
0.1487%. To opisuje ovu tačku i ne
predstavlja opšti dokaz o mogućnosti oscilacija.

## Praćenje grana i pretraga kandidata

Za svaki parametar koristi se 81 logaritamski raspoređena vrednost,
od 0.01 do 100 puta njegove podrazumevane vrednosti. Ostali parametri i
svih šest bilansnih konstanti ostaju fiksni: ukupno 729 izračunatih tačaka.
U svakoj se računaju svih 22 ravnotežnih koncentracija i svih 16 svojstvenih
vrednosti analitičkog redukovanog Jakobijana sa aktuelnim parametrima.

| Parametar | Najmanja vrednost | Najveća vrednost | Najveći max Re λ u pretrazi s⁻¹ | Hopf kandidati |
|---|---:|---:|---:|---:|
| k11 | 0.000165 | 1.65 | -7.3346266e-06 | 0 |
| k7 | 2.8e-06 | 0.028 | -1.2035965e-06 | 0 |
| k8 | 4.08e-05 | 0.408 | -1.2035965e-06 | 0 |
| k6 | 1.54e-08 | 0.000154 | -1.2035965e-06 | 0 |
| a10 | 3e-06 | 0.03 | -7.4471923e-07 | 0 |
| k13 | 0.000138 | 1.38 | -7.2600536e-06 | 0 |
| a11 | 0.0009 | 9 | -7.2622437e-06 | 0 |
| k19 | 1.2e-06 | 0.012 | -1.1722336e-06 | 0 |
| k4 | 1.2e-06 | 0.012 | -7.4188072e-08 | 0 |

U svim ispitanim tačkama pronađen je jedan fizički koren skalarne ravnotežne
jednačine. Svi vodeći realni delovi su negativni; skalarni izvod F'(n*) je
pozitivan na pronađenim korenima. Nije nađen kandidat za Hopf ili fold na
ispitanim mrežama. Numerički podaci i svih 16 spektara čuvaju se uz svaki parametar.

Kompleksni modovi se između susednih tačaka povezuju uparivanjem svojstvenih
vrednosti. Promena znaka realnog dela kompleksnog moda daje kandidata koji
se dodatno lokalizuje i proverava transverzalnost. Ni takav kandidat bez
provere nedegenerisanosti i odgovarajuće analize periodičnih orbita ne bi bio
predstavljen kao potpuno potvrđena Hopf bifurkacija.

Metoda je praćenje po prirodnom parametru sa rešavanjem skalarne ravnoteže;
nije pseudoarclength nastavak. Uske oblasti nestabilnosti, tangencijalni koreni,
druge grane i zajedničko menjanje više parametara mogu izmaći ovoj pretrazi.
Ne tvrdimo da bifurkacije ne postoje za sve pozitivne parametre.

## Nezavisne provere

- Puni model integrisan je do 5000000 s isključivo kao numerička
  provera asimptotskog stanja, bez tumačenja tog trajanja kao eksperimentalne prognoze.
  Najveće odstupanje završnog stanja od izračunate ravnoteže je
  1.329e-14 μM.
- Newtonov metod primenjen je na nezavisni sistem od 16 ravnotežnih jednačina.
  Najveća razlika je 1.693e-15 μM.
- Testovi proveravaju reziduale, bilanse, fizičku prihvatljivost,
  aktuelne parametre u Jakobijanu, nezavisno diferenciranje i granične uslove.

## Granice trenutne analize i nastavak

Metoda u aplikaciji trenutno zahteva strogo pozitivne kinetičke parametre i
svih šest bilansnih konstanti. k4=0, k19=0 i ukupni IL1=0 zahtevaju posebnu
granicu sa dodatnim očuvanjima ili familijama ravnoteža; program ih jasno odbija
umesto prikazivanja pogrešne stabilnosti izolovane ravnoteže.

Sledeći koraci su provera literature i porekla parametara, eventualna
višeparametarska pretraga i zasebna analiza graničnih slučajeva. Ako je cilj
trajno aktivan IKK ili trajni IL1 odziv, reakcije obnavljanja, recikliranja ili
spoljnog ulaza moraju se posebno biološki opravdati. Takve izmene nisu uvedene
u ovaj potvrđeni model.

## Metodski izvori

Brentov metod: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html

Kuznetsov i saradnici, MatCont metodologija nastavka i detekcije bifurkacija:
https://biblio.ugent.be/publication/2129354

Govaerts i saradnici, definisanje Hopf tačaka za parametarske Jakobijane:
https://doi.org/10.1137/S0036142995282182

Ovi izvori služe metodama analize; ne navode se kao izvori bioloških parametara.
