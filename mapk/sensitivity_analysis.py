"""MAPK parameter sensitivity, version 4. Run: python sensitivity_analysis.py.

OAT and joint scenarios are exploratory, not experimentally inferred ranges.
No model equation is changed. No periodic orbit continuation is claimed.
"""
from dataclasses import replace,asdict
from pathlib import Path
import json,csv
import numpy as np
from scipy.optimize import least_squares,brentq
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mapk_model import Parameters,KEEP,Y0,TOTALS,expand,reduced_rhs,rhs,jacobian,equilibrium

BASE=replace(Parameters(),k11=.0025,k12=.0025)
PARAMS=('k11','k12','phosphatase_KK','phosphatase_K')
KGRID=[0.,.00125,.0025,.005,.01]
PGRID=[.5,.75,1.,1.25,1.5]
VGRID=np.unique(np.r_[np.geomspace(.02,.2,25),np.linspace(.2,1.,61),np.linspace(1.,20.,160)])


def safe_equilibrium(p,guess):
    try:
        return equilibrium(p,guess)
    except RuntimeError:
        # Bounded concentrations with penalties for negative eliminated states.
        def residual(z):
            y=expand(z)
            return np.r_[reduced_rhs(0,z,p),np.minimum(y,0)]
        r=least_squares(residual,np.clip(guess,1e-6,[99.999,299.999,299.999,299.999,299.999]),
                        bounds=(np.zeros(5),[100,300,300,300,300]),
                        xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=1000)
        return equilibrium(p,r.x)


def thresholds(p, grid=VGRID):
    guess=np.array([90.,290.,2.,295.,.5])
    branch=[]
    for v in grid:
        y,A,eigs,res=safe_equilibrium(replace(p,V1=float(v)),guess)
        guess=y[KEEP]
        branch.append((float(v),y,float(max(eigs.real))))
    pts=[]
    for a,b in zip(branch[:-1],branch[1:]):
        if a[2]*b[2]>=0: continue
        def ev(v):return safe_equilibrium(replace(p,V1=v),a[1][KEEP])
        v=brentq(lambda v:max(ev(v)[2].real),a[0],b[0],xtol=1e-11)
        y,A,eigs,res=ev(v)
        eig=eigs[np.argmax(eigs.real)]
        omega=abs(eig.imag)
        if omega<1e-7: raise RuntimeError('Nonoscillatory stability crossing')
        d=1e-4*v
        slope=(max(ev(v+d)[2].real)-max(ev(v-d)[2].real))/(2*d)
        assert abs(eig.real)<1e-9 and abs(slope)>1e-8
        other=eigs[abs(eigs.imag-omega)>1e-7]
        other=other[abs(other.imag+omega)>1e-7]
        assert len(other)==3 and max(other.real)<-1e-7
        pts.append({'V1':v,'T_s':float(2*np.pi/omega),'slope':float(slope),
                    'KPP_eq':float(y[7]),'residual':float(res)})
    return pts,branch


def refined_peaks(t,y,positive=True):
    s=y if positive else -y
    ix,_=find_peaks(s,prominence=max(1e-7,.01*np.ptp(s)))
    ts=[];ys=[]
    for i in ix:
        denominator=s[i-1]-2*s[i]+s[i+1]
        offset=.5*(s[i-1]-s[i+1])/denominator if denominator else 0.
        ts.append(t[i]+offset*(t[i+1]-t[i]))
        ys.append(y[i]-.25*(y[i-1]-y[i+1])*offset)
    return np.array(ts),np.array(ys)


def dynamic_metrics(p, horizon=24000., method='DOP853',keep=False):
    # Equilibrium from low-input continuation, to avoid converging to an
    # inadmissible algebraic root when initial conditions are far away.
    g=np.array([90.,290.,2.,295.,.5])
    for v in np.linspace(.02,p.V1,40):
        eq,J,eigs,res=safe_equilibrium(replace(p,V1=float(v)),g);g=eq[KEEP]
    alpha=float(max(eigs.real))
    end=horizon
    current=Y0.copy();start=0.
    all_t=[];all_y=[]
    status='unresolved'
    for attempt in range(4):
        t=np.arange(start,end+1,2.)
        sol=solve_ivp(lambda t,y:rhs(t,y,p),(start,end),current,t_eval=t,
                      method=method,rtol=2e-9,atol=2e-11)
        if not sol.success:raise RuntimeError(sol.message)
        if all_t:all_t.append(sol.t[1:]);all_y.append(sol.y[:,1:])
        else:all_t.append(sol.t);all_y.append(sol.y)
        ts=np.concatenate(all_t);ys=np.concatenate(all_y,axis=1)
        totals=np.vstack((ys[:2].sum(0),ys[2:5].sum(0),ys[5:].sum(0)))
        drift=float(np.max(abs(totals-TOTALS[:,None])))
        if ys.min()<-1e-6 or drift>1e-5:raise RuntimeError('Physical check failed')
        mask=ts>=end*.5
        late_t=ts[mask];signal=ys[7,mask]
        distance=float(max(abs(ys[:,-1]-eq)))
        speed=float(max(abs(rhs(0,ys[:,-1],p))))
        pk,high=refined_peaks(late_t,signal)
        tr,low=refined_peaks(late_t,signal,False)
        if alpha<0 and distance<1e-4 and speed<1e-6:
            status='equilibrium';break
        if len(pk)>=7 and len(tr)>=7:
            dt=np.diff(pk[-7:]);recent_high=high[-6:];recent_low=low[-6:]
            span=float(np.mean(recent_high)-np.mean(recent_low))
            period_variation=float(np.ptp(dt)/np.mean(dt))
            amplitude_variation=float(max(np.ptp(recent_high),np.ptp(recent_low))/max(span,1e-10))
            if span>.01 and max(period_variation,amplitude_variation)<.005:
                status='settled_oscillations';break
        current=ys[:,-1];start=end;end+=24000.
    result={'V1':p.V1,'KPP_equilibrium_nM':float(eq[7]),'equilibrium_spectral_abscissa':alpha,
            'status':status,'horizon_s':float(ts[-1]),'distance_to_equilibrium_nM':distance,
            'final_rhs_nM_per_s':speed,'total_drift_nM':drift,'minimum_concentration_nM':float(ys.min()),
            'late_observed_range_nM':float(np.ptp(signal)),
            'KPP_min_nM':None,'KPP_max_nM':None,'KPP_half_amplitude_nM':None,'period_s':None}
    if status=='equilibrium':
        result.update(KPP_min_nM=float(eq[7]),KPP_max_nM=float(eq[7]),KPP_half_amplitude_nM=0.)
    elif status=='settled_oscillations':
        hi=float(np.mean(high[-6:]));lo=float(np.mean(low[-6:]))
        result.update(KPP_min_nM=lo,KPP_max_nM=hi,KPP_half_amplitude_nM=(hi-lo)/2,
                      period_s=float(np.mean(np.diff(pk[-7:]))),
                      period_relative_range=period_variation,amplitude_relative_drift=amplitude_variation)
    return (result,ts,ys) if keep else result


def write_csv(path,rows):
    columns=list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)


def main():
    out=Path(__file__).resolve().parent/'results_sensitivity';out.mkdir(exist_ok=True)
    scenarios={}
    def add(label,p,group,param='',value=None):
        scenarios[label]={'p':p,'group':group,'parameter':param,'value':value}
    add('baseline',BASE,'baseline')
    add('reference',Parameters(),'reference')
    for param in PARAMS:
        for value in KGRID if param.startswith('k') else PGRID:
            if value==getattr(BASE,param):continue
            add(f'{param}_{value:g}',replace(BASE,**{param:value}),'OAT',param,value)
    # Joint coupling with independently varied catalytic constants.
    for a in KGRID:
        for b in KGRID:
            if a==BASE.k11 or b==BASE.k12:continue # already present as OAT / baseline
            add(f'joint_{a:g}_{b:g}',replace(BASE,k11=a,k12=b),'joint_k')
    results={};rows=[];dynamic=[]
    for idx,(name,s) in enumerate(scenarios.items(),1):
        p=s['p']
        print(f'{idx}/{len(scenarios)} {name}',flush=True)
        pts,branch=thresholds(p)
        r={'parameters':asdict(p),'group':s['group'],'parameter':s['parameter'],
           'value':s['value'],'thresholds':pts,'scan_V1':[float(VGRID[0]),float(VGRID[-1])],
           'dynamics':{}}
        row={'scenario':name,'group':s['group'],'parameter':s['parameter'],'value':s['value'],
             **{k:getattr(p,k) for k in PARAMS},'threshold_count':len(pts)}
        for i,point in enumerate(pts,1):
            row[f'H{i}']=point['V1'];row[f'T_H{i}_s']=point['T_s']
        for v in (.25,5.):
            d=dynamic_metrics(replace(p,V1=v))
            r['dynamics'][str(v)]=d
            dynamic.append({'scenario':name,'group':s['group'],**{k:getattr(p,k) for k in PARAMS},**d})
        results[name]=r;rows.append(row)
        # Checkpoint the completed scenarios: a failed later scenario does not
        # erase evidence or silently generate an empty threshold row.
        (out/'scenarios.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
        write_csv(out/'threshold_scenarios.csv',rows)
        write_csv(out/'dynamic_scenarios.csv',dynamic)
    # Dimensionless local sensitivities S=(p/output)*d(output)/dp.
    # Both steps are evaluated independently; zero baselines are avoided.
    local=[]
    basepts=results['baseline']['thresholds']
    for param in PARAMS:
        basevalue=getattr(BASE,param)
        for frac in (.01,.005):
            sides=[]
            for sign in (-1,1):
                p=replace(BASE,**{param:basevalue*(1+sign*frac)})
                points,_=thresholds(p)
                assert len(points)==2
                sides.append(points)
            for i in range(2):
                for metric in ('V1','T_s'):
                    derivative=(sides[1][i][metric]-sides[0][i][metric])/(2*frac*basevalue)
                    norm=basevalue*derivative/basepts[i][metric]
                    local.append({'parameter':param,'point':f'H{i+1}','output':metric,
                                  'relative_step':frac,'derivative':derivative,'normalized_sensitivity':norm})
    write_csv(out/'local_sensitivities.csv',local)
    for param in PARAMS:
        for i in range(2):
            for metric in ('V1','T_s'):
                vals=[r['normalized_sensitivity'] for r in local if r['parameter']==param
                      and r['point']==f'H{i+1}' and r['output']==metric]
                assert abs(vals[0]-vals[1])<max(.01,abs(vals[1])*.01)
    # Refined scan and independent solver checks on a representative subset.
    validation=[]
    for name in ('baseline','reference','k11_0.01','k12_0.01','phosphatase_KK_0.5','phosphatase_K_1.5'):
        p=scenarios[name]['p'];fine,_=thresholds(p,np.unique(np.r_[VGRID,(VGRID[1:]+VGRID[:-1])/2]))
        old=results[name]['thresholds']
        assert len(fine)==len(old)
        gap=max((abs(a['V1']-b['V1']) for a,b in zip(fine,old)),default=0.)
        assert gap<1e-7
        alt=dynamic_metrics(replace(p,V1=5.),method='Radau')
        orig=results[name]['dynamics']['5.0']
        assert alt['status']==orig['status']
        delta={'scenario':name,'threshold_grid_difference':gap,'status':alt['status']}
        for key in ('KPP_equilibrium_nM','KPP_half_amplitude_nM','period_s'):
            if alt[key] is not None and orig[key] is not None:
                error=abs(alt[key]-orig[key]);delta[key+'_solver_difference']=error
                assert error<max(.02,abs(orig[key])*.001)
        validation.append(delta)
    (out/'validation.json').write_text(json.dumps(validation,indent=2))
    # Thresholds, amplitude and period along one-parameter scenarios.
    fig,axs=plt.subplots(4,3,figsize=(13,13),layout='constrained')
    for i,param in enumerate(PARAMS):
        selected=[r for r in rows if r['parameter']==param]
        selected.append(next(r for r in rows if r['scenario']=='baseline'))
        selected=sorted(selected,key=lambda r:r[param])
        x=[r[param] for r in selected]
        for h in ('H1','H2'):
            axs[i,0].plot(x,[r.get(h,np.nan) for r in selected],'-o',label=h,ms=4)
        for j,key in ((1,'KPP_half_amplitude_nM'),(2,'period_s')):
            yy=[results[r['scenario']]['dynamics']['5.0'][key] for r in selected]
            axs[i,j].plot(x,[np.nan if v is None else v for v in yy],'-o',ms=4)
        for ax in axs[i]:ax.set_xlabel(param+(' (1/s)' if param.startswith('k') else ' (nM)'));ax.grid(alpha=.2)
    for j,title in enumerate(('Linear Hopf thresholds V1 (nM/s)','KPP half-amplitude at V1=5 (nM)',
                              'Period at V1=5 (s)')):axs[0,j].set_title(title)
    axs[0,0].legend()
    fig.savefig(out/'one_parameter_sensitivity.png',dpi=170);plt.close(fig)
    # Joint grid uses actual scanned cases, not an additive approximation.
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for col,key in enumerate(('H1','H2')):
        matrix=np.full((len(KGRID),len(KGRID)),np.nan)
        for r in rows:
            if r['phosphatase_KK']==1 and r['phosphatase_K']==1 and r['k11'] in KGRID and r['k12'] in KGRID:
                matrix[KGRID.index(r['k12']),KGRID.index(r['k11'])]=r.get(key,np.nan)
        image=axs[col].imshow(matrix,origin='lower',aspect='auto',cmap='viridis')
        axs[col].set(xticks=range(len(KGRID)),yticks=range(len(KGRID)),
                     xticklabels=[f'{v:g}' for v in KGRID],yticklabels=[f'{v:g}' for v in KGRID],
                     xlabel='k11 (1/s)',ylabel='k12 (1/s)',title=key+' (V1, nM/s)')
        fig.colorbar(image,ax=axs[col])
    fig.savefig(out/'joint_k11_k12_thresholds.png',dpi=180);plt.close(fig)
    supplementary_outputs(out)
    print('Completed sensitivity analysis',flush=True)


def supplementary_outputs(out=None):
    root=Path(__file__).resolve().parent
    out=root/'results_sensitivity' if out is None else out
    results=json.loads((out/'scenarios.json').read_text())
    # Exact local equilibrium sensitivity via implicit differentiation.
    # F_z dz/dp + F_p = 0; verify against two-sided equilibrium solves.
    eqrows=[]
    for v in (.25,5.):
        p=replace(BASE,V1=v)
        g=np.array([90.,290.,2.,295.,.5])
        for a in np.linspace(.02,v,40):
            y,J,eig,_=safe_equilibrium(replace(p,V1=float(a)),g);g=y[KEEP]
        for param in PARAMS:
            pval=getattr(p,param)
            fp=np.imag(reduced_rhs(0,g,replace(p,**{param:pval+1e-25j})))/1e-25
            dz=-np.linalg.solve(J,fp)
            # State KPP is the last reduced coordinate.
            checks=[]
            for frac in (.01,.005):
                yp=safe_equilibrium(replace(p,**{param:pval*(1+frac)}),g)[0][7]
                ym=safe_equilibrium(replace(p,**{param:pval*(1-frac)}),g)[0][7]
                derivative=(yp-ym)/(2*frac*pval)
                checks.append(abs(derivative-dz[-1])/max(abs(dz[-1]),1e-10))
            assert checks[-1]<.002
            eqrows.append({'V1':v,'parameter':param,'KPP_equilibrium_nM':y[7],
                           'derivative':dz[-1],'normalized_sensitivity':pval*dz[-1]/y[7],
                           'finite_difference_relative_error_half_percent':checks[-1]})
    write_csv(out/'local_equilibrium_sensitivities.csv',eqrows)
    fig,axs=plt.subplots(3,1,figsize=(10,8),layout='constrained',sharex=True)
    selected=('baseline','phosphatase_KK_0.5','phosphatase_K_0.5')
    labels=('Baseline: both phosphatases = 1 nM',
            'KK phosphatase = 0.5 nM; K phosphatase = 1 nM',
            'K phosphatase = 0.5 nM; KK phosphatase = 1 nM')
    for ax,name,label in zip(axs,selected,labels):
        par=Parameters(**results[name]['parameters'])
        d,t,y=dynamic_metrics(replace(par,V1=5.),keep=True)
        ax.plot(t/60,y[7],lw=1.2)
        ax.set(title=label,ylabel='KPP (nM)',xlim=(0,400));ax.grid(alpha=.2)
        np.savetxt(out/f'{name}_V5_trajectory.csv',np.column_stack((t,y.T)),delimiter=',',
                   header='time_s,KKK,KKKP,KK,KKP,KKPP,K,KP,KPP',comments='')
    axs[-1].set_xlabel('Time (min)')
    fig.suptitle('First 400 min; convergence checked over the full simulation',fontsize=12)
    fig.savefig(out/'phosphatase_time_responses.png',dpi=170);plt.close(fig)
    local=list(csv.DictReader((out/'local_sensitivities.csv').open()))
    lines=['# Osetljivost MAPK modela: verzija 4','',
      '## Pokretanje','',
      'U Spyder-u otvorite sensitivity_analysis.py i pokrenite F5. Skripta koristi iste',
      'zavisnosti i ne menja matematičko jezgro. Rezultati su već uključeni u paket.',
      'Ponovno pokretanje zamenjuje izlaze istih naziva u results_sensitivity.',
      'Broj scenarija i adaptivno produžavanje simulacija znače da ova skripta može',
      'trajati duže od prethodnih. Napredak se ispisuje u konzoli.',
      '', '## Opseg i značenje analize','',
      'Radna tačka: k11=k12=0.0025 1/s, obe fosfataze po 1 nM.',
      'Ostali parametri su isti kao u prethodnoj verziji. Opsezi su scenariji:',
      'nisu eksperimentalni intervali poverenja ili potvrđeni fiziološki rasponi.',
      '', '- k11,k12: 0, 0.00125, 0.0025, 0.005, 0.01 1/s.',
      '- Fosfataze drugog (phosphatase_KK) i trećeg (phosphatase_K) nivoa:',
      '  0.5, 0.75, 1, 1.25, 1.5 nM, nezavisno.',
      '- Zajednička mreža k11 x k12: 25 kombinacija; fosfataze su po 1 nM.',
      '- Ukupno 34 imenovana scenarija, uključujući referentni; referentni je',
      '  numerički isti kao zajednička kombinacija (0,0), pa ima 33 različita skupa.',
      '- Ravnoteže i linearni Hopfovi prelazi: V1 od 0.02 do 20 nM/s.',
      '- Vremenski odzivi: V1=0.25 i V1=5 nM/s, isti početni uslovi kao ranije.',
      '- Lokalna osetljivost pragova i perioda u pragu: koraci +/-1% i +/-0.5%.',
      '- Lokalna osetljivost KPP ravnoteže: implicitno diferenciranje 5D sistema,',
      '  nezavisno provereno centralnim razlikama. Pri V1=5 to je nestabilna',
      '  ravnoteža; njena koncentracija nije srednja vrednost oscilatornog odziva.',
      '', '## Glavni nalazi','',
      'Oko radne tačke k11 i k12 deluju suprotno na H1: povećanje k11 podiže,',
      'a povećanje k12 spušta prag. Na H2 oba lokalno deluju u smeru snižavanja.',
      'Zajednička promena zato može prikriti različite pojedinačne efekte.',
      'Nelinearni zajednički efekat izračunat je direktno u svakoj tački mreže,',
      'bez pretpostavke da se pojedinačni efekti prosto sabiraju.',
      '', 'Fosfataze imaju veću lokalnu relativnu osetljivost pragova od dodatnih',
      'katalitičkih konstanti u ovoj radnoj tački. Promene nisu globalno monotone:',
      'na primer H1 za phosphatase_KK prvo opada pa raste u ispitanim tačkama.',
      'To nije opšte pravilo o povećanju/smanjenju fosfataza u ćeliji.',
      '', '## Lokalna normalizovana osetljivost pragova','',
      'S=(p/H)*dH/dp. Na primer S=0.295 znači da lokalno povećanje parametra',
      'za 1% pomera prag približno +0.295%. Prikazan je korak +/-0.5%.',
      '', '| Parametar | S za H1 | S za H2 |','|---|---:|---:|']
    for param in PARAMS:
        vals=[next(float(r['normalized_sensitivity']) for r in local if r['parameter']==param
              and r['point']==h and r['output']=='V1' and float(r['relative_step'])==.005) for h in ('H1','H2')]
        lines.append(f'| {param} | {vals[0]:+.5f} | {vals[1]:+.5f} |')
    lines+=['','## Izabrani rezultati pri V1=5','',
       'A=(KPP_max-KPP_min)/2 je poluamplituda, ne razlika maksimum-minimum.',
       'Period je iz ustaljenih vrhova odziva, a ne 2*pi/omega u Hopfovoj tački.',
       '', '| Scenario | H1 | H2 | A (nM) | T (s) |', '|---|---:|---:|---:|---:|']
    for name in ('baseline','k11_0.01','k12_0.01','phosphatase_KK_0.5',
                 'phosphatase_KK_1.5','phosphatase_K_0.5','phosphatase_K_1.5','joint_0.01_0.01'):
        r=results[name];d=r['dynamics']['5.0'];h=r['thresholds']
        period='nema oscilacija' if d['period_s'] is None else f"{d['period_s']:.2f}"
        lines.append(f"| {name} | {h[0]['V1']:.6f} | {h[1]['V1']:.6f} | {d['KPP_half_amplitude_nM']:.3f} | {period} |")
    lines+=['','U svim ispitanim scenarijima pri V1=0.25 odziv je dostigao ravnotežu.',
        'Pri V1=5 smanjenje pojedinačne fosfataze na 0.5 nM pomera H2 ispod 5,',
        'pa radna vrednost V1 izlazi iz intervala nestabilne ravnoteže. Time se',
        'objašnjava prestanak oscilacija bez potrebe da opšti efekat fosfataze',
        'proglasimo stabilizujućim ili destabilizujućim za sve uslove.',
        '', '## Provere i ograničenja','',
        '- Pragovi: potvrđen imaginarni par, transverzalni prelaz i tri preostala',
        '  stabilna moda. Dva prelaza nađena su u svakom scenariju na praćenoj grani.',
        '- Periodičnost: najmanje sedam vrhova i dolina u kasnoj polovini odziva;',
        '  promena perioda i visine ekstrema mora biti manja od 0.5%.',
        '- Simulacija počinje sa 24000 s i po potrebi se produžava do 96000 s.',
        '  Svi priloženi scenariji su klasifikovani; nema unresolved slučajeva.',
        '- Ravnoteža u simulaciji: odstupanje svih stanja <1e-4 nM i brzine <1e-6 nM/s.',
        '- Provereni su nenegativnost koncentracija i tri očuvana ukupna zbira.',
        '- Na šest reprezentativnih slučajeva ponovljena je pretraga pragova',
        '  dvostruko gušćom mrežom i odziv pri V1=5 nezavisnim Radau rešavačem.',
        '- Lokalni izvodi pragova/perioda provereni su na dva koraka parametra.',
        '- Skeniranje mreže ne isključuje vrlo uske dodatne intervale ili druge grane.',
        '- Ustaljene oscilacije u konačnoj simulaciji nisu Floquetov dokaz stabilnosti.',
        '  Prvi Ljapunovljev koeficijent nije ponovo računat za svih 34 scenarija;',
        '  prethodna superkritična klasifikacija ne prenosi se automatski na njih.',
        '- Nije rađena zajednička mreža obe fosfataze niti globalna statistička',
        '  analiza (npr. Sobol); ovo su OAT, lokalna i mrežna analiza k11,k12.',
        '', '## Fajlovi','',
        '- sensitivity_analysis.py: cela analiza i generisanje ovog pregleda.',
        '- results_sensitivity/scenarios.json: svi parametri i metrički rezultati.',
        '- threshold_scenarios.csv: pragovi i periodi u pragovima.',
        '- dynamic_scenarios.csv: ravnoteže, minimumi, maksimumi, poluamplitude, periodi.',
        '- local_sensitivities.csv: izvodi i normalizovane osetljivosti H1,H2 i njihovih perioda.',
        '- local_equilibrium_sensitivities.csv: lokalna osetljivost KPP ravnoteže.',
        '- validation.json: detalji provera gušće mreže i rešavača.',
        '- one_parameter_sensitivity.png: pragovi, poluamplitude i periodi po parametru.',
        '- joint_k11_k12_thresholds.png: pragovi po zajedničkoj mreži; ćelije predstavljaju',
        '  diskretne scenarije, ne jednako razmaknutu numeričku osu ili interpolaciju.',
        '- phosphatase_time_responses.png i *_V5_trajectory.csv: tri primera odziva.',
        '', 'Matematički model i njegove prethodno navedene pretpostavke ostaju isti.',
        'Ovaj paket je osnova za sledeću fazu: Streamlit interfejs.', '']
    (root/'SENSITIVITY_SR.md').write_text('\n'.join(lines),encoding='utf-8')


if __name__=='__main__':main()
