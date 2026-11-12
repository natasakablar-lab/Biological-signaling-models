"""MAPK Research Lab. Start with: python -m streamlit run app.py"""
from pathlib import Path
from dataclasses import asdict,replace
from io import BytesIO
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import streamlit as st
from mapk_model import Parameters,NAMES,KEEP
from sensitivity_analysis import dynamic_metrics,thresholds

ROOT=Path(__file__).resolve().parent
st.set_page_config(page_title='MAPK Research Lab',page_icon='🧬',layout='wide')
PRESETS={
    'Radna tačka · oscilacije':dict(V1=5.,k11=.0025,k12=.0025,phosphatase_KK=1.,phosphatase_K=1.),
    'Radna tačka · ravnoteža':dict(V1=.25,k11=.0025,k12=.0025,phosphatase_KK=1.,phosphatase_K=1.),
    'Referentni model · oscilacije':dict(V1=5.,k11=0.,k12=0.,phosphatase_KK=1.,phosphatase_K=1.),
    'Referentni model · ravnoteža':dict(V1=.25,k11=0.,k12=0.,phosphatase_KK=1.,phosphatase_K=1.),
    'Fosfataza KK · 0.5 nM':dict(V1=5.,k11=.0025,k12=.0025,phosphatase_KK=.5,phosphatase_K=1.),
    'Fosfataza K · 0.5 nM':dict(V1=5.,k11=.0025,k12=.0025,phosphatase_KK=1.,phosphatase_K=.5),
}


def preset_changed():
    for key,value in PRESETS[st.session_state.preset].items():
        st.session_state['input_'+key]=value


@st.cache_data(show_spinner=False,max_entries=32)
def run_simulation(settings,horizon,method):
    p=Parameters(**settings)
    metrics,t,y=dynamic_metrics(p,horizon=horizon,method=method,keep=True)
    # Solve equilibrium by continuation at the SAME active parameter values.
    from sensitivity_analysis import safe_equilibrium
    guess=np.array([90.,290.,2.,295.,.5])
    for v in np.linspace(.02,p.V1,40):
        eq,J,eigs,res=safe_equilibrium(replace(p,V1=float(v)),guess);guess=eq[KEEP]
    return dict(parameters=settings,metrics=metrics,t=t,y=y,equilibrium=eq,
                jacobian=J,eigenvalues=eigs,equilibrium_residual=res,method=method)


@st.cache_data(show_spinner=False,max_entries=32)
def run_hopf(settings,upper):
    p=Parameters(**settings)
    grid=np.unique(np.r_[np.geomspace(.02,.2,25),np.linspace(.2,1.,61),np.linspace(1.,upper,220)])
    points,branch=thresholds(p,grid)
    return dict(parameters=settings,upper=upper,points=points,
                branch=np.array([[v,*y,a] for v,y,a in branch]))


@st.cache_data(show_spinner=False,max_entries=16)
def run_lyapunov(settings,points):
    from lyapunov_analysis import derivative_builder,coefficient
    from sensitivity_analysis import safe_equilibrium
    p=Parameters(**settings)
    derivatives=derivative_builder(p)
    guess=np.array([90.,290.,2.,295.,.5]);rows=[]
    for i,h in enumerate(points,1):
        for v in np.linspace(.02,h['V1'],60):
            y,J,eigs,res=safe_equilibrium(replace(p,V1=float(v)),guess);guess=y[KEEP]
        A,B,C=derivatives(guess,h['V1'])
        l1,G21,terms,checks=coefficient(A,B,C)
        near_zero=abs(l1)<1e-9
        rows.append(dict(tacka=f'H{i}',V1=h['V1'],l1=l1,
                         tip='blizu nule · dodatna provera' if near_zero else 'superkritična' if l1<0 else 'subkritična',
                         ReG21_over_2=float(G21.real/2),
                         strana_malog_ciklusa='neodređeno' if near_zero else ('iznad praga' if h['slope']>0 else 'ispod praga')
                         if l1<0 else ('ispod praga' if h['slope']>0 else 'iznad praga'),
                         stabilnost_malog_ciklusa='neodređeno' if near_zero else 'stabilan' if l1<0 else 'nestabilan'))
    return rows


def csv_bytes(frame):return frame.to_csv(index=False).encode('utf-8-sig')


def png_bytes(t,y,selected):
    fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
    for label in selected:ax.plot(t/60,y[NAMES.index(label)],label=label,lw=1.2)
    ax.set(xlabel='Time (min)',ylabel='Concentration (nM)');ax.grid(alpha=.2)
    if selected:ax.legend()
    buf=BytesIO();fig.savefig(buf,format='png',dpi=200);plt.close(fig)
    return buf.getvalue()


st.title('MAPK Research Lab')
st.caption('Rekonstrukcija modela N. A. Kablar · MTNS 2010 · simulacija i dinamička analiza')
with st.sidebar:
    st.subheader('Parametri modela')
    st.selectbox('Početni scenario',list(PRESETS),key='preset',on_change=preset_changed)
    if 'input_V1' not in st.session_state:preset_changed()
    with st.form('simulation_settings'):
        st.number_input('V1 · nM/s',min_value=.02,max_value=20.,step=.05,format='%.6f',key='input_V1')
        st.number_input('k11 · 1/s',min_value=0.,max_value=.025,step=.00025,format='%.6f',key='input_k11')
        st.number_input('k12 · 1/s',min_value=0.,max_value=.025,step=.00025,format='%.6f',key='input_k12')
        st.number_input('Fosfataza KK · nM',min_value=.25,max_value=2.,step=.05,key='input_phosphatase_KK')
        st.number_input('Fosfataza K · nM',min_value=.25,max_value=2.,step=.05,key='input_phosphatase_K')
        horizon=st.selectbox('Početno trajanje · s',[12000,24000,48000],index=1,key='horizon')
        method=st.selectbox('Numerički rešavač',['DOP853','Radau'],key='solver')
        submitted=st.form_submit_button('Primeni i simuliraj',type='primary')
    st.caption('Simulacija se po potrebi produžava. Posle promene parametara pritisni „Primeni i simuliraj”.')
    st.caption('k11 i k12 su radne procene. n=1; ostali parametri ostaju kao u rekonstrukciji.')

if submitted or 'simulation' not in st.session_state:
    settings=asdict(replace(Parameters(),**{k:st.session_state['input_'+k] for k in PRESETS[next(iter(PRESETS))]}))
    try:
        with st.spinner('Računam odziv i proveravam ustaljivanje…'):
            result=run_simulation(settings,float(horizon),method)
        st.session_state.simulation=result
        st.session_state.pop('hopf',None);st.session_state.pop('lyapunov',None)
    except Exception as exc:
        st.error(f'Izračunavanje nije završeno: {exc}')
        if 'simulation' not in st.session_state:st.stop()
        st.warning('Ispod su prethodni uspešno izračunati rezultati; novi parametri nisu primenjeni.')

r=st.session_state.simulation;m=r['metrics'];active=r['parameters']
st.info(f"Primenjeno: V1={active['V1']:g} nM/s · k11={active['k11']:g}, k12={active['k12']:g} 1/s · "
        f"fosfataze KK/K={active['phosphatase_KK']:g}/{active['phosphatase_K']:g} nM")
status={'equilibrium':'Dostignuta ravnoteža','settled_oscillations':'Ustaljene oscilacije','unresolved':'Ustaljivanje nije potvrđeno'}
cols=st.columns(4)
cols[0].metric('Odziv',status[m['status']])
cols[1].metric('Poluamplituda KPP', '—' if m['KPP_half_amplitude_nM'] is None else f"{m['KPP_half_amplitude_nM']:.3f} nM")
cols[2].metric('Period', '—' if m['period_s'] is None else f"{m['period_s']:.2f} s")
cols[3].metric('Simulirano vreme',f"{m['horizon_s']/60:.0f} min")
if m['status']=='unresolved':st.warning('Prikazani odziv još nije zadovoljio kriterijume ustaljivanja. Amplituda i period nisu potvrđeni.')

tabs=st.tabs(['Odzivi','Ravnoteža i stabilnost','Hopf i Ljapunov','Osetljivost','Model i preuzimanje'])
with tabs[0]:
    selected=st.multiselect('Koncentracije',list(NAMES),default=['KKKP','KKPP','KPP'])
    late=st.checkbox('Prikaži samo poslednju polovinu simulacije',value=False)
    mask=r['t']>=r['t'][-1]/2 if late else np.ones(len(r['t']),dtype=bool)
    fig=go.Figure()
    for label in selected:
        fig.add_trace(go.Scatter(x=r['t'][mask]/60,y=r['y'][NAMES.index(label),mask],name=label,mode='lines'))
    fig.update_layout(xaxis_title='Vreme (min)',yaxis_title='Koncentracija (nM)',height=430,
                      margin=dict(l=20,r=20,t=20,b=20),legend=dict(orientation='h'))
    st.plotly_chart(fig,width='stretch')
    c1,c2=st.columns(2)
    with c1:
        fx=st.selectbox('Fazni prikaz · horizontalna osa',list(NAMES),index=4)
    with c2:
        fy=st.selectbox('Fazni prikaz · vertikalna osa',list(NAMES),index=7)
    phase=go.Figure(go.Scatter(x=r['y'][NAMES.index(fx),mask],y=r['y'][NAMES.index(fy),mask],mode='lines'))
    phase.update_layout(xaxis_title=fx+' (nM)',yaxis_title=fy+' (nM)',height=350,margin=dict(t=10,b=20))
    st.plotly_chart(phase,width='stretch')
    st.caption('Ustaljene oscilacije u konačnoj simulaciji nisu zamena za Floquetovu analizu. Poluamplituda je (maksimum−minimum)/2.')
    st.download_button('Preuzmi grafik · PNG',png_bytes(r['t'][mask],r['y'][:,mask],selected),file_name='MAPK_odziv.png',mime='image/png')

with tabs[1]:
    eig=r['eigenvalues'];alpha=float(max(eig.real))
    stability='lokalno asimptotski stabilna' if alpha < -1e-8 else 'nestabilna' if alpha>1e-8 else 'blizu kritične tačke'
    st.subheader('Ravnoteža je '+stability)
    st.caption('Stabilnost se odnosi na konzervacionu mnogostrukost sa pet nezavisnih stanja. Ravnoteža može biti nestabilna iako odziv ima ustaljene oscilacije.')
    c1,c2=st.columns(2)
    c1.dataframe(pd.DataFrame({'Stanje':NAMES,'Ravnoteža (nM)':r['equilibrium']}),hide_index=True,width='stretch')
    c2.dataframe(pd.DataFrame({'Re λ (1/s)':eig.real,'Im λ (1/s)':eig.imag}),hide_index=True,width='stretch')
    st.write('Jakobijan u koordinatama: KKK, KK, KKPP, K, KPP')
    st.dataframe(pd.DataFrame(r['jacobian'],index=[NAMES[i] for i in KEEP],columns=[NAMES[i] for i in KEEP]),width='stretch')
    st.caption(f"Rezidual ravnoteže: {r['equilibrium_residual']:.2e} nM/s · najveće odstupanje ukupnih količina: {m['total_drift_nM']:.2e} nM")

with tabs[2]:
    st.subheader('Ravnotežna grana i linearni Hopfovi pragovi')
    st.caption('Menja se V1; ostali primenjeni parametri ostaju fiksni. Ova pretraga ne nastavlja periodične orbite.')
    upper=st.selectbox('Gornja granica pretrage V1 · nM/s',[10.,20.],index=1)
    if st.button('Izračunaj Hopfove pragove'):
        try:
            with st.spinner('Pratim ravnotežnu granu…'):
                st.session_state.hopf=run_hopf(active,upper)
            st.session_state.pop('lyapunov',None)
        except Exception as exc:
            st.error(f'Pretraga nije završena: {exc}')
    if 'hopf' in st.session_state:
        h=st.session_state.hopf
        st.caption(f"Izračunati opseg: 0.02–{h['upper']:g} nM/s.")
        if not h['points']:st.warning('Na praćenoj grani nije pronađen prelaz u zadatom opsegu. To ne isključuje druge grane ili prelaze van opsega.')
        else:st.dataframe(pd.DataFrame(h['points']).rename(columns={'V1':'V1 (nM/s)','T_s':'Period u pragu (s)','slope':'dReλ/dV1','KPP_eq':'KPP u pragu (nM)'}),hide_index=True,width='stretch')
        b=h['branch'];fig=go.Figure()
        for stable,label,dash,color in ((True,'Stabilna ravnoteža','solid','#17859b'),(False,'Nestabilna ravnoteža','dash','#c45c38')):
            fig.add_trace(go.Scatter(x=b[:,0],y=np.where((b[:,-1]<0)==stable,b[:,8],np.nan),name=label,line=dict(dash=dash,color=color)))
        for point in h['points']:fig.add_vline(x=point['V1'],line_dash='dot',line_color='#996529')
        fig.update_layout(xaxis_title='V1 (nM/s)',yaxis_title='Ravnotežna KPP (nM)',height=380)
        st.plotly_chart(fig,width='stretch')
        st.download_button('Preuzmi pragove · CSV',csv_bytes(pd.DataFrame(h['points'])),file_name='MAPK_Hopf.csv',mime='text/csv')
        if h['points'] and st.button('Izračunaj prvi Ljapunovljev koeficijent'):
            try:
                with st.spinner('Računam izvode do trećeg reda i koeficijente…'):
                    st.session_state.lyapunov=run_lyapunov(active,h['points'])
            except Exception as exc:st.error(f'Koeficijent nije izračunat: {exc}')
        if 'lyapunov' in st.session_state:
            st.dataframe(pd.DataFrame(st.session_state.lyapunov),hide_index=True,width='stretch')
            st.caption('Kuznetsovljeva konvencija: l1=Re(G21)/(2ω), ||q||₂=1, <p,q>=1; koncentracije u nM. Klasifikacija je lokalna.')
            st.download_button('Preuzmi koeficijente · CSV',csv_bytes(pd.DataFrame(st.session_state.lyapunov)),file_name='MAPK_Lyapunov.csv',mime='text/csv')

with tabs[3]:
    st.subheader('Provereni scenariji osetljivosti')
    st.info('Ovaj pregled prikazuje prethodno izračunate scenarije verzije 4. Podaci ne zavise od trenutnih podešavanja u bočnoj traci.')
    path=ROOT/'results_sensitivity'/'scenarios.json'
    if path.exists():
        saved=json.loads(path.read_text())
        name=st.selectbox('Scenario iz analize osetljivosti',list(saved))
        s=saved[name];params=s['parameters']
        st.write({key:params[key] for key in ('k11','k12','phosphatase_KK','phosphatase_K')})
        st.dataframe(pd.DataFrame(s['thresholds']),hide_index=True,width='stretch')
        st.dataframe(pd.DataFrame(s['dynamics']).T, width='stretch')
        for fname,label in [('one_parameter_sensitivity.png','Pojedinačni uticaji'),('joint_k11_k12_thresholds.png','Zajednička promena k11 i k12')]:
            with st.expander(label):st.image(str(ROOT/'results_sensitivity'/fname),width='stretch')
        st.download_button('Preuzmi kompletnu tabelu scenarija · CSV',(ROOT/'results_sensitivity'/'dynamic_scenarios.csv').read_bytes(),file_name='MAPK_sensitivity.csv',mime='text/csv')
    else:st.warning('Nedostaje results_sensitivity. Raspakuj ceo paket ili prethodno pokreni sensitivity_analysis.py.')

with tabs[4]:
    st.subheader('Matematički model')
    st.markdown('Tri nivoa: **KKK → KK → K**. Izlaz **KPP** negativnom povratnom spregom utiče na aktivaciju prvog nivoa. Reakcije sa **k11** i **k12** omogućavaju dodatnu aktivaciju trećeg nivoa preko KKP.')
    st.latex(r'v_{11}=\frac{k_{11}[KKP][K]}{K_{M11}+[K]},\qquad v_{12}=\frac{k_{12}[KKP][KP]}{K_{M12}+[KP]}')
    st.markdown('Ukupne koncentracije su **100 / 300 / 300 nM**. Početno su aktivirani oblici nula. Cross-talk i pozitivna povratna sprega su isključeni; Hillov koeficijent je n=1. Dodatne konstante i koncentracije fosfataza imaju status modelskih pretpostavki.')
    st.caption('Aplikacija služi istraživanju dinamike ovog modela. Nije model Alzheimerove bolesti; tu vezu tek treba zasnovati na literaturi i eksperimentalnim podacima.')
    st.dataframe(pd.DataFrame(active.items(),columns=['Parametar','Primenjena vrednost']),hide_index=True,width='stretch')
    frame=pd.DataFrame(r['y'].T,columns=NAMES);frame.insert(0,'time_s',r['t'])
    st.download_button('Preuzmi vremenske podatke · CSV',csv_bytes(frame),file_name='MAPK_trajectory.csv',mime='text/csv')
    payload={'parameters':active,'metrics':m,'solver':r['method'],'model':'Kablar MTNS 2010 reconstruction v5'}
    st.download_button('Preuzmi parametre i sažetak · JSON',json.dumps(payload,indent=2,ensure_ascii=False).encode(),file_name='MAPK_run.json',mime='application/json')
    st.download_button('Preuzmi Jakobijan · CSV',csv_bytes(pd.DataFrame(r['jacobian'],columns=[NAMES[i] for i in KEEP])),file_name='MAPK_Jacobian.csv',mime='text/csv')
    st.markdown('**Dokumentacija:** README_SR.md, HOPF_SR.md, LYAPUNOV_SR.md i SENSITIVITY_SR.md u preuzetom paketu.')
