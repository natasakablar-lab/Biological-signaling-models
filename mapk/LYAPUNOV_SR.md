# Prvi Ljapunovljev koeficijent: MAPK, verzija 3

## Pokretanje

Raspakujte ceo paket u jedan folder. U istom Python okruženju instalirajte
zavisnosti:

```bash
python -m pip install -r requirements.txt
python lyapunov_analysis.py
```

Ako već imate ostale pakete iz verzije 2, nova zavisnost je SymPy:

```bash
python -m pip install sympy
```

U Spyder-u otvorite lyapunov_analysis.py i pokrenite ga pomoću F5. Skripta
učitava rezultate pragova iz results_hopf/hopf_summary.json, koji su
uključeni u paket. Ako promenite model ili parametre, prvo ponovo pokrenite
hopf_analysis.py da osvežite te podatke, a zatim lyapunov_analysis.py.
Rezultati se čuvaju u results_lyapunov; postojeći fajlovi istih naziva
zamenjuju se pri ponovnom pokretanju. Prethodne skripte ostaju dostupne.

## Formula i normalizacija

Koristimo pet nezavisnih koncentracija z=(KKK,KK,KKPP,K,KPP), izraženih
u nM, sa vremenom u sekundama. Ne vršimo dodatno skaliranje koncentracija.
U svakoj Hopfovoj tački razvijamo desnu stranu oko ravnoteže:

f(z*+u)=A u + B(u,u)/2 + C(u,u,u)/6 + O(||u||^4).

A je Jakobijan, B je sirovi drugi, a C sirovi treći izvod (faktorijeli
nisu ugrađeni u B i C). Desni i adjungovani sopstveni vektor zadovoljavaju:

A q = i*omega*q; A^T p = -i*omega*p;
conj(p)^T q = 1; conj(q)^T q = 1.

Uzimamo:

h11 = -A^(-1) B(q,conj(q));
h20 = (2*i*omega*I-A)^(-1) B(q,q);
G21 = <p,C(q,q,conj(q)) + 2*B(q,h11) + B(conj(q),h20)>;
l1 = Re(G21)/(2*omega).

Inverzne matrice se u kodu ne formiraju eksplicitno: rešavaju se linearni
sistemi. Ovaj 5D sistem je ključan: puni 8D Jakobijan ima tri nulte
sopstvene vrednosti zbog konzervacije i ne može se invertovati u formuli.

Ovo je Kuznetsovljeva konvencija. Brojna vrednost l1 zavisi od skaliranja
koordinata i q, dok je znak relevantan za klasifikaciju. Izabrana
normalizacija ||q||_2=1 uklanja proizvoljno skaliranje q u ovim koordinatama.
Koeficijenti ovde ne treba direktno da se porede po veličini sa softverom
koji koristi drugačiju normalizaciju.

U originalnom vremenu, odgovarajući realni kubni koeficijent u kompleksnoj
normalnoj formi je Re(G21)/2 = omega*l1. Posle pozitivnog skaliranja vremena
tau=omega*t, radijalni kubni koeficijent je l1. I ovaj originalni vremenski
koeficijent je snimljen, kako bi konvencija bila potpuno proverljiva.

Izvor formule: Yu. A. Kuznetsov, Elements of Applied Bifurcation Theory,
Second Edition, Springer, 1998, jednačina (5.62), štampana strana 178
(PDF strana 199). Primer 5.4, strane 178-181, korišćen je i za proveru.

https://www.ma.ic.ac.uk/~dturaev/kuznetsov.pdf

Sažeto objašnjenje autora, uključujući normalizaciju i znak koeficijenta:
https://var.scholarpedia.org/article/Andronov-Hopf

## Rezultati

| Scenario | Tačka | V1 (nM/s) | l1 (navedena konvencija) | Tip |
|---|---|---:|---:|---|
| Referentni: k11=k12=0 | H1 | 0.3842452864 | -8.73913664e-5 | Superkritični |
| Referentni: k11=k12=0 | H2 | 6.2480479397 | -1.30598229e-4 | Superkritični |
| Procenjeni: k11=k12=0.0025 1/s | H1 | 0.3862273464 | -1.11566831e-4 | Superkritični |
| Procenjeni: k11=k12=0.0025 1/s | H2 | 6.1485359755 | -1.24576135e-4 | Superkritični |

Sva četiri koeficijenta su negativna i numerički jasno različita od nule.
Zajedno sa prethodno potvrđenom transverzalnošću i preostalim stabilnim
modovima, to daje numeričku klasifikaciju sva četiri Hopfova slučaja kao
nedegenerisana i superkritična na odgovarajućoj konzervacionoj mnogostrukosti.
Ovo je prvi Ljapunovljev koeficijent normalne forme u Hopfovoj tački,
a ne Ljapunovljev eksponent vremenske trajektorije.

## Na kojoj strani praga su stabilni mali ciklusi?

Neka je mu=Re(lambda), a s=dmu/dV1 u pragu. Lokalno:

dr/dt = mu*r + (omega*l1)*r^3 + ...,
mu približno s*(V1-V1_H).

Za l1<0 stabilni mali ciklus postoji na strani mu>0, sa vodećom procenom
r^2 približno -s*(V1-V1_H)/(omega*l1). To važi dovoljno blizu praga.

- H1: s>0. Stabilni mali ciklusi su za V1 malo IZNAD H1. Rastom V1
  ravnoteža gubi stabilnost i nastaju stabilne oscilacije.
- H2: s<0. Stabilni mali ciklusi su za V1 malo ISPOD H2. Rastom V1 do H2
  mali ciklus se skuplja u ravnotežu, a ravnoteža postaje stabilna.

Oba slučaja su superkritična. Obrnut smer promene stabilnosti u H2 nije
razlog da se H2 proglasi subkritičnim.

Lokalni rezultat ne dokazuje da se ciklusi iz H1 i H2 globalno spajaju u
jednu granu, da postoji jedinstven ciklus kroz ceo interval, ili da nema
drugih globalnih bifurkacija. Za to su potrebni nastavak periodičnih
orbita i, za njihovu numeričku stabilnost, Floquetovi multiplikatori.

## Provere

1. Drugi i treći izvodi dobijeni su simboličkom diferencijacijom ISTE
   funkcije koja se koristi za simulacije. Nije zasebno prepisan model.
2. Simbolički Jakobijan poklapa se sa Jakobijanom kompleksnog koraka.
3. Nezavisno su računati B i C preko centralnih razlika numeričkog
   Jakobijana za korake 0.04,0.02,0.01 i 0.005 nM. Relativna odstupanja
   l1 pri poslednjem koraku manja su od 9e-8 u sva četiri slučaja.
4. Provereni su reziduali sopstvenih vektora, <p,q>=1 i ||q||=1.
5. Promena faze q za 0.731 rad ne menja l1 (uz ponovno normiranje p).
6. Za tačan kubni radijalni oscilator sa a=+/-1 i omega=1.7 dobijeno je
   l1=2*a/omega, što proverava znak i faktor normalizacije.
7. Za Kuznetsovljev 3D primer 5.4 pri omega=1, sa kvadratnom
   nelinearnošću, dobijeno je l1=-0.3 uz ||q||=1. Objavljena vrednost je
   -0.9 uz q=(1,i,-1); podela sa ||q||^2=3 daje naš rezultat. Ovaj primer
   posebno proverava članove sa B koji ne postoje u kubnom benchmarku.

Sve provere prolaze. Priložene rezultate generisala je SymPy 1.14.0;
minimalna verzija u requirements.txt nije zasebno testirana.

## Izlazni fajlovi

- results_lyapunov/lyapunov_coefficients.csv: glavna tabela.
- results_lyapunov/lyapunov_summary.json: koeficijenti, G21, tri doprinosa
  formuli, normalizacije i rezultati nezavisnih provera.
- results_lyapunov/*_derivatives.npz: A,B,C, ravnoteža u 5D i V1;
  može se otvoriti sa numpy.load radi dalje analize.

Parametri dodatnih reakcija ostaju radne procene iz verzije 1. Negativni
l1 i matematička klasifikacija sami po sebi nisu eksperimentalna potvrda
njihove biološke vrednosti.
