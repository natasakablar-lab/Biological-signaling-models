# IL1 NFκB — puni ODE, redukovani ODE i DAE

Autorka je 5. novembra 2026. potvrdila izvedene jednačine, bazalni izvor k6,
ispravljene jedinice k6 i k19 i x1 kao aktivni IKK_P. Specifična komponenta
je označena kao TRAF6 (u ranijim tabelama TRAF). Poreklo i eksperimentalna
primenljivost kinetičkih parametara ostaju za proveru literature.

## Puni ODE model

$$\dot x=f(x;p)=S v(x;p),\quad x\in\mathbb R^{22}.$$

S ima 22 reda, 35 kolona i rang 16. Sve reakcije i svih 22 eksplicitnih
jednačina nalaze se u model.py. Puni model ostaje referenca za koncentracije.
Periodični reseti iz istorijskog MATLAB Main fajla nisu uključeni.

## Šest konstantnih bilansa

$$c_1=x_4+x_5+x_6+x_7+x_{10}+x_{11},$$
$$c_2=x_{12}+x_{14}+x_{16}+x_{18}+x_{20}+x_{22},$$
$$c_3=x_{13}+x_{14}+x_{16}+x_{18}+x_{20}+x_{22},$$
$$c_4=x_{15}+x_{16}+x_{18}+x_{20}+x_{22},$$
$$c_5=x_{17}+x_{18}+x_{20}+x_{22},$$
$$c_6=x_{19}+x_{20}+x_{22}.$$

Redosled c odgovara NFκB, IL1, ILR, MyD88, IRAK, TRAF6. Važi c=Cx0 i CS=0.
Podrazumevano c=(0.1,0.15,0.1,0.05,0.05,0.05) μM. Konstante se računaju
iz unetih početnih uslova i ostaju fiksne tokom pojedinačne simulacije.

## Redukovani ODE

$$q=(x_1,x_2,x_3,x_5,x_6,x_7,x_8,x_9,x_{10},x_{11},x_{14},x_{16},x_{18},x_{20},x_{21},x_{22})^T.$$

Šest rekonstruisanih koncentracija:

$$x_4=c_1-x_5-x_6-x_7-x_{10}-x_{11},$$
$$x_{12}=c_2-x_{14}-x_{16}-x_{18}-x_{20}-x_{22},$$
$$x_{13}=c_3-x_{14}-x_{16}-x_{18}-x_{20}-x_{22},$$
$$x_{15}=c_4-x_{16}-x_{18}-x_{20}-x_{22},$$
$$x_{17}=c_5-x_{18}-x_{20}-x_{22},$$
$$x_{19}=c_6-x_{20}-x_{22}.$$

A=(4,12,13,15,17,19) označava ove indekse (brojanje od 1), D ostalih 16.
Matrica B ima jedinične redove na D, a na A redove −C_D. Vektor b ima c na
pozicijama A i nule na D. U originalnom redosledu koncentracija:

$$x=Bq+b(c),\qquad \dot q=f_D(Bq+b(c);p).$$
$$CB=0,\quad Cb=c,\qquad J_{red}=f_{D,x}(Bq+b(c);p)B.$$

Ovo je tačna redukcija pomoću bilansa. Nema kvazistacionarne aproksimacije,
pretpostavke brzih reakcija niti odsecanja koncentracija. Analitički
reduced_jacobian je dostupan za kasniju analizu ravnoteža i bifurkacija na
fiksnoj klasi očuvanja c.

## DAE i descriptor zapis

U originalnom redosledu svih 22 koncentracije:

$$F_D=\dot x_D-f_D(x;p)=0,\qquad F_A=Cx-c=0.$$

Šest odabranih ODE zamenjeno je bilansima, umesto dodavanja bilansa svim ODE.
Ekvivalentno:

$$E\dot x=H(x;c,p),\qquad H_D=f_D(x;p),\quad H_A=c-Cx.$$

E je dijagonalna matrica 22×22: Eii=1 za i u D, Eii=0 za i u A. Rang E je 16.
Pošto je C_A=I6, algebarske koncentracije se jednoznačno određuju iz preostalih,
a DAE je diferencijalnog indeksa 1. To je ekvivalentna DAE postavka potvrđenog
modela pomoću zakona očuvanja.

Saglasni početni uslovi zadovoljavaju:

$$Cx_0=c,\qquad \dot x_0=f(x_0;p),\qquad C\dot x_0=0.$$

Program proverava saglasnost ako se eksplicitno zadaju konstante c.
DAE koristi pravi implicitni rešavač SUNDIALS IDA preko scikit-sundae.
Ne postoji tihi prelazak na ODE integraciju. Bez tog paketa puni i redukovani
ODE ostaju dostupni, a izbor DAE daje jasno uputstvo za instalaciju.

Dokumentacija rešavača:
https://scikit-sundae.readthedocs.io/en/stable/api/sksundae/ida/
https://sundials.readthedocs.io/en/latest/ida/Mathematics_link.html

## Sedmi bilans

$$T_{IKK}=x_1+x_3+x_5+x_7+x_{21}+x_{22},\qquad \dot T_{IKK}=-k_4x_1.$$

To nije konstanta i ne daje dodatnu eliminaciju bez nove dinamičke promenljive.
Aktivni i neaktivni IKK ostaju različite koncentracije.

## Provere

verify_formulations.py proverava četiri scenarija: uslove iz tabela,
bez početnog IL1, bez početnog aktivnog IKK i pozitivne početne komplekse
uz izmenjene kinetičke parametre. Sva 22 odziva porede se na istoj mreži.

Provereni su CS=0, CB=0, rangovi S i E, početna konzistentnost, šest bilansa,
DAE reziduali i nenegativnost uz numeričku toleranciju. Analitički Jakobijan
proveren je nezavisnim complex-step diferenciranjem.

Linearizovani descriptor par (H_x,E) ima 16 konačnih i šest beskonačnih modova.
Konačni spektar poklapa se sa spektrom J_red. Test je izvršen u početnoj tački,
koja nije ravnoteža; rezultat se zato ne tumači kao analiza stabilnosti ravnoteže.

Stvarne greške su u results/verification_formulations.json. Matrice C, B, E,
mapa koordinata i svi odzivi dostupni su i kao CSV/JSON. Slaganje zapisa
potvrđuje matematičku i numeričku ekvivalentnost, a ne eksperimentalnu
validaciju biološkog modela. Bifurkaciona i sensitivity analiza su naredna faza.
