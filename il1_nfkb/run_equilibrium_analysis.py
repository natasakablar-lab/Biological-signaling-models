"""Reproduce equilibrium, stability and one-parameter branch scans."""
from pathlib import Path
from dataclasses import asdict
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import root
from model import *
from formulations import *
from equilibria import *


def run():
    out=Path(__file__).parent/'results_equilibria'; out.mkdir(exist_ok=True)
    p=Parameters(); c=totals(initial_conditions()); eq=equilibrium_roots(p,c)[0]
    info=assess_equilibrium(eq,p,c)
    report={'parameters':asdict(p),'fixed_totals':c.tolist(),'equilibrium':eq.tolist(),'baseline_metrics':info,
            'scans':{},'scope':'Natural-parameter scans of detected physical branches; no global completeness claim.'}
    with (out/'equilibrium_22_states.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['state','species','concentration_uM'])
        w.writerows((f'x{i+1}',NAMES[i],float(v)) for i,v in enumerate(eq))
    np.savetxt(out/'equilibrium_spectrum.csv',np.column_stack([info['eigen_real'],info['eigen_imag']]),delimiter=',',header='real_per_s,imag_per_s',comments='')
    # Independent long transient, followed by a 16-D Newton solve.
    long=simulate_formulation('full',p,t_end=5e6,points=1001,rtol=1e-10,atol=1e-13)
    nr=root(lambda q:reduced_rhs(0,q,c,p),long.y[DIFFERENTIAL,-1],jac=lambda q:reduced_jacobian(q,c,p),tol=1e-11)
    newton_x=reconstruct(nr.x,c)
    report['independent_checks']={'long_time_s':float(long.t[-1]),
      'long_transient_max_difference_uM':float(np.max(np.abs(long.y[:,-1]-eq))),
      'newton_success':bool(nr.success),'newton_max_difference_uM':float(np.max(np.abs(newton_x-eq))),
      'newton_full_residual_inf':float(np.max(np.abs(rhs(0,newton_x,p))))}
    np.savetxt(out/'long_transient_22_states.csv',np.column_stack([long.t,long.y.T]),delimiter=',',header=','.join(['time_s']+[f'x{i+1}' for i in range(22)]),comments='')
    for parameter in SCAN_PARAMETERS:
        result=branch_scan(parameter,p,c,points=81)
        rows=result['rows']
        with (out/f'branch_{parameter}.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        spectra=result['spectra']; np.savez_compressed(out/f'spectrum_{parameter}.npz',values=[r['value'] for r in rows],eigenvalues=spectra)
        report['scans'][parameter]={'factor_range':[.01,100.],'points':81,
          'max_spectral_abscissa':max(r['spectral_abscissa'] for r in rows),
          'min_spectral_abscissa':min(r['spectral_abscissa'] for r in rows),
          'maximum_detected_roots':max(r['detected_roots'] for r in rows),
          'max_equilibrium_residual':max(r['full_residual_inf'] for r in rows),
          'min_scalar_slope_at_root':min(r['scalar_slope_at_root'] for r in rows),
          'hopf_candidates':result['hopf_candidates']}
        print(parameter,report['scans'][parameter],flush=True)
        fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
        values=np.array([r['value'] for r in rows])
        axs[0,0].plot(values,[r['x10'] for r in rows],color='seagreen'); axs[0,0].set_ylabel('Equilibrium NF-κB$_n$ (μM)')
        axs[0,1].plot(values,[r['x2'] for r in rows],color='seagreen'); axs[0,1].set_ylabel('Equilibrium IκB (μM)')
        axs[1,0].plot(values,[r['spectral_abscissa'] for r in rows]); axs[1,0].axhline(0,color='black',lw=.8); axs[1,0].set_ylabel('Max Re λ (s⁻¹)')
        for j in range(16): axs[1,1].plot(values,spectra[:,j].real,lw=.8)
        axs[1,1].axhline(0,color='black',lw=.8); axs[1,1].set_ylabel('All 16 Re λ (s⁻¹)')
        for ax in axs.flat: ax.set_xscale('log');ax.set_xlabel(parameter+' ('+parameter_unit(parameter)+')');ax.grid(alpha=.2)
        fig.suptitle('Physical equilibrium branch | '+parameter+' | fixed six conserved totals')
        fig.savefig(out/f'branch_{parameter}.png',dpi=170);plt.close(fig)
    (out/'analysis_summary.json').write_text(json.dumps(report,indent=2))
    print('INDEPENDENT',report['independent_checks'],flush=True)
    return report

if __name__=='__main__': run()
