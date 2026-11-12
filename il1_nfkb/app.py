"""Start with: streamlit run app.py"""
from dataclasses import asdict
import io,json
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from model import *
from formulations import *
from equilibria import equilibrium_roots, assess_equilibrium, branch_scan, SCAN_PARAMETERS

st.set_page_config(page_title='IL-1 · NF-κB laboratorija',layout='wide')
st.title('IL-1 · NF-κB biološki modul')
st.caption(TITLE+' | ODE · DAE · Ravnoteže i stabilnost | Verzija 3')
st.info('Autorka je potvrdila izvedene jednačine, bazalni izvor k6, jedinice k6 i k19 i x1 = aktivni IKK_P. Sva tri zapisa prikazuju svih 22 koncentracije. Poreklo i eksperimentalna primenljivost parametara ostaju za proveru literature.')

presets={'Tabele objavljenog rada':'published','Bez IL-1, sa početnim aktivnim IKK':'without_IL1','IL-1 uz početni IKK_P = 0':'IL1_only'}
with st.sidebar:
    st.header('Simulacija')
    preset_label=st.selectbox('Početni scenario',list(presets))
    preset=presets[preset_label]
    t_end=st.number_input('Kraj simulacije (s)',min_value=1.,max_value=1000000.,value=50000.,step=1000.)
    points=st.number_input('Broj izlaznih tačaka',min_value=101,max_value=10001,value=2501,step=100)
    mode_labels={'Puni ODE · 22':'full','Redukovani ODE · 16':'reduced','DAE · 16 + 6':'dae'}
    formulation=mode_labels[st.selectbox('Matematički zapis',list(mode_labels))]
    method=st.selectbox('ODE rešavač',['Radau','BDF','LSODA'])
    st.caption('DAE uvek koristi SUNDIALS IDA. Izbor iznad važi za ODE zapise.')
    checkall=st.checkbox('Uporedi sva tri zapisa',value=False,disabled=not ida_available())
    if not ida_available(): st.caption('Za DAE i poređenje sva tri zapisa: pip install -r requirements-dae.txt')
    compare=st.checkbox('Uporedi sa scenarijem bez IL-1',value=True)
    st.caption('Poređenje: postavlja samo x12(0)=0, uz ostale unete parametre i početne uslove. Za prilagođene početne komplekse x14, x16, x18, x20 ili x22 ovaj scenario i dalje sadrži vezani IL-1.')

params0=asdict(Parameters())
param_df=pd.DataFrame({'Parametar':list(params0),'Vrednost':list(params0.values()),
                      'Jedinica':[parameter_unit(k) for k in params0]})
ydef=initial_conditions(preset)
y_df=pd.DataFrame({'Stanje':[f'x{i+1}' for i in range(22)],'Vrsta':NAMES,'Vrednost (μM)':ydef})
with st.form('inputs'):
    col1,col2=st.columns(2)
    with col1:
        st.subheader('Parametri')
        pedit=st.data_editor(param_df,disabled=['Parametar','Jedinica'],hide_index=True,
            key='params',height=360,column_config={'Vrednost':st.column_config.NumberColumn(format='%.8g',min_value=0.)})
    with col2:
        st.subheader('Početne koncentracije')
        yedit=st.data_editor(y_df,disabled=['Stanje','Vrsta'],hide_index=True,
            key='initial_'+preset,height=360,column_config={'Vrednost (μM)':st.column_config.NumberColumn(format='%.8g',min_value=0.)})
    submitted=st.form_submit_button('Izračunaj odzive',type='primary')

@st.cache_data(show_spinner=False)
def compute(pdct,y0,tend,n,solver,formulation):
    pp=Parameters(**pdct)
    s=simulate_formulation(formulation,pp,y0,tend,n,solver)
    return s,diagnostics(s,pp)

# The comparison checkbox is an action: calculate immediately on first enable.
comparison_requested = checkall and ('result' not in st.session_state or not st.session_state.result[-1])
if submitted or 'result' not in st.session_state or comparison_requested:
    try:
        pdct=dict(zip(pedit['Parametar'],pedit['Vrednost'].astype(float)))
        y0=yedit['Vrednost (μM)'].to_numpy(dtype=float)
        with st.spinner('Računam odzive...'):
            sol,diag=compute(pdct,y0,float(t_end),int(points),method,formulation)
            base=None
            if compare:
                ybase=y0.copy(); ybase[11]=0
                base,_=compute(pdct,ybase,float(t_end),int(points),method,formulation)
            all_solutions={}
            if checkall:
                for mode in ['full','reduced','dae']:
                    all_solutions[mode]=sol if mode==formulation else compute(pdct,y0,float(t_end),int(points),method,mode)[0]
        st.session_state.result=(sol,diag,base,pdct,y0.tolist(),method,all_solutions)
    except Exception as exc: st.error(str(exc)); st.stop()

sol,diag,base,pdct,y0,used_method,all_solutions=st.session_state.result
st.caption(f'Poslednji rezultat: 0–{sol.t[-1]:g} s, {len(sol.t)} tačaka, {sol.solver_name}, zapis {sol.formulation}. Posle izmene ulaza kliknite „Izračunaj odzive“.')
tabs=st.tabs(['Odzivi','Fazni portreti','Jednačine i izvori','Bilansi i izvoz','Poređenje zapisa','Ravnoteže i bifurkacije'])
with tabs[0]:
    selected=st.multiselect('Prikazana stanja',list(range(22)),default=DEFAULT_SELECTION,
        format_func=lambda i:f'x{i+1} · {NAMES[i]}')
    fig=go.Figure()
    for i in selected:
        fig.add_trace(go.Scatter(x=sol.t,y=sol.y[i],name=f'x{i+1} {NAMES[i]}'))
        if base is not None:
            fig.add_trace(go.Scatter(x=base.t,y=base.y[i],name=f'x{i+1} bez slobodnog IL-1',line={'dash':'dot'}))
    fig.update_layout(xaxis_title='Vreme (s)',yaxis_title='Koncentracija (μM)',height=510,legend={'orientation':'h'})
    st.plotly_chart(fig,width='stretch')
    rows=[]
    for i in selected:
        k=int(np.argmax(sol.y[i])); rows.append({'Stanje':f'x{i+1} {NAMES[i]}','Maksimum (μM)':sol.y[i,k],
          'Vreme maksimuma na mreži (s)':sol.t[k],'Konačna vrednost (μM)':sol.y[i,-1]})
    st.dataframe(pd.DataFrame(rows),hide_index=True)
    st.caption('Maksimumi se računaju na izlaznoj mreži. Dugački prelazni odzivi sami po sebi ne potvrđuju granični ciklus ili bifurkaciju.')
with tabs[1]:
    cols=st.columns(3)
    ids=[]
    for col,label,default in zip(cols,['X osa','Y osa','Z osa'],[9,1,0]):
        with col: ids.append(st.selectbox(label,list(range(22)),index=default,format_func=lambda i:f'x{i+1} {NAMES[i]}'))
    a,b,c=ids
    fig=go.Figure(go.Scatter3d(x=sol.y[a],y=sol.y[b],z=sol.y[c],mode='lines',line={'color':sol.t,'colorscale':'Viridis','width':4,'colorbar':{'title':'s'}}))
    fig.add_trace(go.Scatter3d(x=[sol.y[a,0],sol.y[a,-1]],y=[sol.y[b,0],sol.y[b,-1]],z=[sol.y[c,0],sol.y[c,-1]],mode='markers',marker={'size':5,'color':['green','red']},name='Početak / kraj'))
    fig.update_layout(height=600,scene={'xaxis_title':NAMES[a]+' (μM)','yaxis_title':NAMES[b]+' (μM)','zaxis_title':NAMES[c]+' (μM)'})
    st.plotly_chart(fig,width='stretch')
with tabs[2]:
    st.write('Radni izvor: Paper With IL-1.pdf, 31. jul 2007, str. 3–7. Parametri i početni uslovi prepisani su i provereni prema četiri dostavljene tabele iz konačnog rada.')
    st.markdown('[Bibliografski zapis koji je dostavila autorka](https://www.semanticscholar.org/paper/Mathematical-Model-of-IL-1-NFkB-Biological-Module-Kablar-Kvrgic/a56c5cd1783d27cf526395abb7e5090c717aa651)')
    st.success('Potvrđeno: reakcija (15) je bazalni izvor ∅ → x8 sa brzinom k6. k6 ima jedinicu μM/s, k19 jedinicu s⁻¹. Specifična komponenta označena je kao TRAF6.')
    st.write('Tabela stanja označava x1 kao IKK_P, a tabela početnih uslova kao IKK. Ovde x1(0)=0,1 pripada aktivnom IKK_P; x21(0)=0,1 neaktivnom IKK. Početni ukupni IKK je 0,2 μM.')
    import sympy as sp
    st.subheader('Redukcija i DAE zapis')
    st.write('Algebarski određene koncentracije: x4, x12, x13, x15, x17, x19. Diferencijalne koordinate: '+', '.join(f'x{i+1}' for i in DIFFERENTIAL)+'.')
    st.latex(r'x=Bq+b(c),\quad \dot q=f_D(Bq+b(c))')
    st.latex(r'\dot x_D=f_D(x),\quad 0=Cx-c,\quad C_A=I_6')
    st.write('DAE je indeksa 1. Šest ODE zamenjeno je bilansima; ukupan IKK nije konstanta i ne daje sedmu eliminaciju.')
    xs,eqs=symbolic_equations()
    for i,expr in enumerate(eqs):
        with st.expander(f'x{i+1} · {NAMES[i]}'):
            st.latex(r'\frac{dx_{'+str(i+1)+r'}}{dt}='+sp.latex(expr))
    st.caption('Svi članovi nastaju iz stehiometrijske matrice. Transkript se ne troši pri translaciji; NF-κB se ne troši pri katalitičkoj transkripciji.')
with tabs[3]:
    st.subheader('Numeričke provere trenutnog odziva')
    st.json(diag)
    st.write('Čuva se ukupni NF-κB i ukupne količine IL-1, IL-R, MyD88, IRAK i TRAF6 u svim kompleksima. Ukupni IKK opada samo razgradnjom slobodnog IKK_P: dT_IKK/dt = −k4·x1.')
    if diag['minimum_concentration'] < -1e-8: st.warning('Uočene su negativne koncentracije iznad numeričke tolerancije.')
    if max(diag['conservation_max_abs_drift'].values())>1e-7: st.warning('Bilansi odstupaju više od 1e-7 μM; proverite tolerancije i ulaze.')
    df=pd.DataFrame(sol.y.T,columns=[f'x{i+1}_{n}' for i,n in enumerate(NAMES)]); df.insert(0,'time_s',sol.t)
    st.download_button('Preuzmi sve odzive CSV',df.to_csv(index=False).encode(),file_name='IL1_NFkB_responses.csv',mime='text/csv')
    if sol.formulation=='dae':
        residuals=np.column_stack([dae_residual(t,x,xp,sol.totals,Parameters(**pdct)) for t,x,xp in zip(sol.t,sol.y.T,sol.yp.T)])
        st.metric('Najveći DAE rezidual',f'{np.max(np.abs(residuals)):.3e}')
    cfg={'formulation':sol.formulation,'solver_actual':sol.solver_name,'parameters':pdct,'initial_conditions':y0,'t_end_s':float(sol.t[-1]),'points':len(sol.t),'ode_solver':used_method,'diagnostics':diag}
    st.download_button('Preuzmi postavku JSON',json.dumps(cfg,indent=2).encode(),file_name='IL1_NFkB_configuration.json',mime='application/json')

with tabs[4]:
    if not all_solutions:
        st.write('U bočnom meniju uključite „Uporedi sva tri zapisa“. Poređenje se računa odmah.')
    else:
        comparison=pd.DataFrame(comparison_table(all_solutions))
        st.write('Isti parametri, početni uslovi i izlazna mreža. Odstupanja su prema punom ODE modelu.')
        state_index=st.selectbox('Koncentracija za poređenje',list(range(22)),index=9,
            format_func=lambda i:f'x{i+1} {NAMES[i]}')
        trajectories=go.Figure()
        for mode,label,dash in [('full','Puni ODE','solid'),('reduced','Redukovani ODE','dash'),('dae','DAE · IDA','dot')]:
            solution=all_solutions[mode]
            trajectories.add_trace(go.Scatter(x=solution.t,y=solution.y[state_index],name=label,line={'dash':dash}))
        trajectories.update_layout(xaxis_title='Vreme (s)',yaxis_title='Koncentracija (μM)',height=400)
        st.plotly_chart(trajectories,width='stretch')
        st.caption('Tri krive se gotovo preklapaju; numeričke razlike su prikazane ispod.')
        st.dataframe(comparison,hide_index=True)
        errfig=go.Figure()
        for mode in ['reduced','dae']:
            delta=np.max(np.abs(all_solutions[mode].y-all_solutions['full'].y),axis=0)
            errfig.add_trace(go.Scatter(x=sol.t,y=delta,name=mode))
        errfig.update_layout(xaxis_title='Vreme (s)',yaxis_title='Najveće odstupanje 22 stanja (μM)',height=380)
        st.plotly_chart(errfig,width='stretch')
        st.download_button('Preuzmi poređenje CSV',comparison.to_csv(index=False).encode(),file_name='comparison_formulations.csv',mime='text/csv')

with tabs[5]:
    st.subheader('Fizičke ravnoteže i lokalna stabilnost')
    st.write('Analiza koristi parametre i šest bilansnih konstanti poslednjeg uspešnog proračuna simulacije. Stabilnost se računa iz Jakobijana redukovanog modela sa 16 koordinata.')
    st.info('Za pozitivne stope k4 i k19 i pozitivne bilansne konstante fizička ravnoteža ima x1=x3=x5=x7=x21=x22=0. Ukupni IKK se u zatvorenom modelu troši. k19 i k4 utiču na prolazni odziv i brzinu relaksacije, ali ne menjaju ovu ravnotežu dok su pozitivni.')
    c_analysis=totals(np.array(y0))
    analysis_signature=json.dumps({'parameters':pdct,'totals':c_analysis.tolist()},sort_keys=True)
    if st.button('Izračunaj ravnotežu',key='equilibrium_button'):
        try:
            eqs=equilibrium_roots(Parameters(**pdct),c_analysis)
            eq=eqs[0]
            metrics=assess_equilibrium(eq,Parameters(**pdct),c_analysis)
            st.session_state.equilibrium_result=(analysis_signature,eq,metrics,len(eqs))
        except Exception as exc: st.error(str(exc))
    er=st.session_state.get('equilibrium_result')
    if er and er[0]==analysis_signature:
        _,eq,metrics,number=er
        st.write(f'Pronađeno ravnoteža na skalarnoj mreži: {number}. Prikazana je prva.')
        c1,c2,c3=st.columns(3)
        c1.metric('Stabilnost',metrics['stability'])
        c2.metric('Najveći Re λ (s⁻¹)',f"{metrics['spectral_abscissa']:.6e}")
        c3.metric('Rezidual svih 22 ODE (μM/s)',f"{metrics['full_residual_inf']:.3e}")
        eqframe=pd.DataFrame({'Stanje':[f'x{i+1}' for i in range(22)],'Vrsta':NAMES,'Ravnoteža (μM)':eq})
        st.dataframe(eqframe,hide_index=True)
        eigenframe=pd.DataFrame({'Re λ (s⁻¹)':metrics['eigen_real'],'Im λ (s⁻¹)':metrics['eigen_imag']})
        st.dataframe(eigenframe,hide_index=True)
        st.download_button('Preuzmi ravnotežu CSV',eqframe.to_csv(index=False).encode(),file_name='equilibrium.csv',mime='text/csv')
    elif er: st.caption('Postavka se promenila. Ponovo izračunajte ravnotežu za aktuelne parametre i bilanse.')
    st.subheader('Praćenje ravnotežne grane po parametru')
    with st.form('branch_inputs'):
        parameter=st.selectbox('Parametar za analizu',SCAN_PARAMETERS,key='branch_parameter')
        bcols=st.columns(3)
        low=bcols[0].number_input('Najmanji faktor početne vrednosti',min_value=.000001,value=.01,format='%.6g')
        high=bcols[1].number_input('Najveći faktor početne vrednosti',min_value=.00001,value=100.,format='%.6g')
        nodes=bcols[2].number_input('Tačke na logaritamskoj mreži',min_value=5,max_value=401,value=81,step=10)
        scan_submitted=st.form_submit_button('Prati ravnotežnu granu')
    if scan_submitted:
        try:
            with st.spinner('Računam ravnotežnu granu i svih 16 svojstvenih vrednosti...'):
                branch=branch_scan(parameter,Parameters(**pdct),c_analysis,float(low),float(high),int(nodes))
            st.session_state.branch_result=(analysis_signature,branch)
        except Exception as exc: st.error(str(exc))
    br=st.session_state.get('branch_result')
    if br and br[0]==analysis_signature:
        branch=br[1]; frame=pd.DataFrame(branch['rows']); key=branch['parameter']
        st.caption(f'Prikazano prethodno izračunato praćenje parametra {key}.')
        chosen=st.selectbox('Ravnotežna koncentracija za dijagram',list(range(22)),index=9,format_func=lambda i:f'x{i+1} {NAMES[i]}')
        branchplot=go.Figure(go.Scatter(x=frame['value'],y=frame[f'x{chosen+1}'],mode='lines+markers',marker={'color':frame['spectral_abscissa'],'colorscale':'RdBu','showscale':True,'colorbar':{'title':'max Re λ'}}))
        branchplot.update_layout(xaxis_title=key+' ('+parameter_unit(key)+')',xaxis_type='log',yaxis_title='Ravnotežna koncentracija (μM)',height=380)
        st.plotly_chart(branchplot,width='stretch')
        spectral=go.Figure(go.Scatter(x=frame['value'],y=frame['spectral_abscissa'],name='max Re λ'))
        spectral.add_hline(y=0,line_color='black')
        spectral.update_layout(xaxis_title=key,xaxis_type='log',yaxis_title='Najveći Re λ (s⁻¹)',height=340)
        st.plotly_chart(spectral,width='stretch')
        if branch['hopf_candidates']: st.warning('Pronađeni su kandidati za Hopf. Potrebna je dodatna provera nedegenerisanosti.'); st.json(branch['hopf_candidates'])
        else: st.success('Na ovoj grani i u ovoj mreži nije detektovan kandidat za Hopf bifurkaciju.')
        st.caption('Praćenje koristi pozitivan parametar i skalarno rešavanje ravnoteže. Ovo nije pseudoarclength nastavak ni globalni dokaz odsustva bifurkacija. Uži prelazi, tangencijalni koreni i druge grane mogu zahtevati dodatnu analizu.')
        st.dataframe(frame[['value','spectral_abscissa','stability','full_residual_inf','scalar_slope_at_root','detected_roots']],hide_index=True)
        st.download_button('Preuzmi ravnotežnu granu CSV',frame.to_csv(index=False).encode(),file_name=f'branch_{key}.csv',mime='text/csv')
    elif br: st.caption('Parametri ili bilansi su promenjeni. Ponovite praćenje grane.')
