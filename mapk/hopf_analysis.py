"""Equilibrium continuation and linear Hopf thresholds for MAPK.

Run beside mapk_model.py: python hopf_analysis.py
This computes equilibrium branches, NOT periodic-orbit continuation.
"""
from dataclasses import replace, asdict
from pathlib import Path
import json
import csv
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mapk_model import Parameters, KEEP, Y0, equilibrium, reduced_rhs, expand


def trace_branch(p, grid, initial=None):
    guess = np.array([81.94,272.22,7.85,288.5,2.574]) if initial is None else initial
    rows = []
    for value in grid:
        y,J,eig,residual = equilibrium(replace(p,V1=float(value)),guess)
        guess = y[KEEP]
        rows.append((float(value), y, eig, residual))
    return rows


def leading_pair(eig):
    positive = eig[eig.imag > 1e-7]
    if not len(positive):
        raise RuntimeError('No complex eigenvalue pair')
    return positive[np.argmax(positive.real)]


def locate(p, grid, initial=None):
    branch = trace_branch(p,grid,initial)
    alpha = np.array([max(row[2].real) for row in branch])
    crossings = np.flatnonzero(alpha[:-1]*alpha[1:] < 0)
    points = []
    for idx in crossings:
        guess = branch[idx][1][KEEP]
        def evaluate(v):
            return equilibrium(replace(p,V1=float(v)),guess)
        value = brentq(lambda v: max(evaluate(v)[2].real),
                       grid[idx],grid[idx+1],xtol=2e-12,rtol=1e-13)
        y,J,eig,residual = evaluate(value)
        pair = leading_pair(eig)
        delta = 1e-4*value
        slope = (leading_pair(evaluate(value+delta)[2]).real
                 -leading_pair(evaluate(value-delta)[2]).real)/(2*delta)
        delta2 = delta/2
        slope2 = (leading_pair(evaluate(value+delta2)[2]).real
                  -leading_pair(evaluate(value-delta2)[2]).real)/(2*delta2)
        others = eig[np.abs(eig.imag-abs(pair.imag))>1e-7]
        others = others[np.abs(others.imag+abs(pair.imag))>1e-7]
        # Independent Jacobian: central differences with two step sizes.
        def finite_eigen(h):
            x=y[KEEP]
            par=replace(p,V1=value)
            E=np.eye(5)
            finite=np.column_stack([(reduced_rhs(0,x+h*E[j],par)
                                     -reduced_rhs(0,x-h*E[j],par))/(2*h)
                                    for j in range(5)])
            return leading_pair(np.linalg.eigvals(finite))
        fd1,fd2 = finite_eigen(1e-4),finite_eigen(5e-5)
        assert abs(pair.real)<1e-10
        assert abs(slope)>1e-7 and np.all(others.real < -1e-6)
        assert abs(fd1-pair)<1e-7 and abs(fd2-pair)<1e-7
        assert abs(slope-slope2)<1e-7
        points.append(dict(V1_nM_per_s=float(value),period_s=float(2*np.pi/pair.imag),
             omega_per_s=float(pair.imag),critical_pair_real_per_s=float(pair.real),
             crossing_slope=float(slope),crossing_slope_half_step=float(slope2),
             other_eigenvalues_per_s=[[float(v.real),float(v.imag)] for v in others],
             all_eigenvalues_per_s=[[float(v.real),float(v.imag)] for v in eig],
             equilibrium_nM=y.tolist(),rhs_residual=float(residual),
             finite_difference_eigenvalue_error=float(max(abs(fd1-pair),abs(fd2-pair))),
             bracket=[float(grid[idx]),float(grid[idx+1])]))
    return branch,points


def local_check(p, value, eq_guess):
    """Tiny modal perturbation; fit oscillatory-mode growth near equilibrium.

    Initial point is physical. Short, very small perturbations check local
    linear stability without confusing critical slowing with a stable cycle.
    """
    p=replace(p,V1=value)
    y,J,eig,_=equilibrium(p,eq_guess)
    vals,vecs=np.linalg.eig(J)
    idx=np.argmax(np.where(vals.imag>1e-7,vals.real,-np.inf))
    lam=vals[idx]
    eigenvector=vecs[:,idx]
    perturb=1e-3*eigenvector.real/np.max(abs(eigenvector.real))
    x0=y[KEEP]+perturb
    assert expand(x0).min()>0
    duration=3*2*np.pi/lam.imag
    t=np.linspace(0,duration,1201)
    sol=solve_ivp(lambda t,x:reduced_rhs(t,x,p),(0,duration),x0,t_eval=t,
                  method='DOP853',rtol=2e-12,atol=2e-13)
    if not sol.success: raise RuntimeError(sol.message)
    # Inverse eigenvector matrix projects onto the complex mode.
    projection=np.linalg.inv(vecs)[idx]
    amplitude=np.abs(projection@(sol.y-y[KEEP,None]))
    fitted=float(np.polyfit(t,np.log(amplitude),1)[0])
    expected=float(lam.real)
    assert abs(fitted-expected)<max(5e-8,abs(expected)*.02)
    return dict(V1_nM_per_s=value,predicted_growth_per_s=expected,
                simulated_growth_per_s=fitted,duration_s=float(duration))


def main():
    out=Path(__file__).resolve().parent/'results_hopf'
    out.mkdir(exist_ok=True)
    grid=np.unique(np.r_[np.linspace(.1,1.,181),np.linspace(1.,9.,321)])
    scenarios={'reference':Parameters(),
               'estimated':replace(Parameters(),k11=.0025,k12=.0025)}
    report={}
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for col,(name,p) in enumerate(scenarios.items()):
        branch,points=locate(p,grid)
        assert len(points)==2, (name,len(points))
        # Repeat continuation from high V1 to low V1, on a shifted grid.
        _,reverse=locate(p,np.linspace(9.,.1,401),branch[-1][1][KEEP])
        reverse=sorted(reverse,key=lambda v:v['V1_nM_per_s'])
        gap=max(abs(a['V1_nM_per_s']-b['V1_nM_per_s']) for a,b in zip(points,reverse))
        assert gap<1e-7
        checks=[]
        for point in points:
            for factor in (.98,1.02):
                checks.append(local_check(p,point['V1_nM_per_s']*factor,
                                          np.array(point['equilibrium_nM'])[KEEP]))
        report[name]={'parameters':asdict(p),'thresholds':points,
                      'reverse_continuation_max_threshold_difference':gap,
                      'local_simulation_checks':checks,
                      'note':'Linear Hopf criteria verified. First Lyapunov coefficient and Floquet multipliers not computed.'}
        a=np.array([max(v[2].real) for v in branch])
        x=np.array([v[0] for v in branch])
        y=np.array([v[1][7] for v in branch])
        for stable,style,label,color in ((True,'-','Stable equilibrium','#167d9a'),
                                          (False,'--','Unstable equilibrium','#c45332')):
            axs[0,col].plot(x,np.where((a<0)==stable,y,np.nan),style,color=color,label=label,lw=2)
        axs[1,col].plot(x,a,color='#344c7b',lw=1.8)
        axs[1,col].axhline(0,color='gray',lw=.8)
        for k,point in enumerate(points,1):
            v=point['V1_nM_per_s']
            for ax in axs[:,col]: ax.axvline(v,color='#996529',ls=':',lw=1)
            axs[0,col].scatter(v,point['equilibrium_nM'][7],color='#996529',zorder=4)
            axs[0,col].annotate(f'H{k}: {v:.6f}',(v,point['equilibrium_nM'][7]),
                               xytext=(7,14 if k==1 else 12),textcoords='offset points',fontsize=9)
        axs[0,col].set(title='Reference (k11=k12=0)' if name=='reference' else 'Assumed k11=k12=0.0025 /s',
                       ylabel='Equilibrium KPP (nM)')
        axs[0,col].legend(fontsize=9,loc='lower right')
        axs[1,col].set(ylabel='Largest real eigenvalue part (1/s)')
        for ax in axs[:,col]:
            ax.set_xlabel('V1 (nM/s)')
            ax.grid(alpha=.18)
        with (out/f'{name}_equilibrium_branch.csv').open('w',newline='') as f:
            w=csv.writer(f)
            w.writerow(['V1_nM_per_s','KKK','KKKP','KK','KKP','KKPP','K','KP','KPP','max_real_eigenvalue_per_s','rhs_residual'])
            for v,y,eig,res in branch: w.writerow([v,*y,max(eig.real),res])
    fig.savefig(out/'hopf_thresholds.png',dpi=180)
    plt.close(fig)
    paper_v=[.3842,6.248]
    paper_t=[637.7,1227.]
    for i,point in enumerate(report['reference']['thresholds']):
        point['paper_V1']=paper_v[i]
        point['paper_period_s']=paper_t[i]
        point['V1_difference_from_paper']=point['V1_nM_per_s']-paper_v[i]
        point['period_difference_from_paper_s']=point['period_s']-paper_t[i]
    (out/'hopf_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    with (out/'hopf_thresholds.csv').open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['scenario','point','V1_nM_per_s','omega_per_s','period_s','crossing_slope'])
        for name,r in report.items():
            for i,h in enumerate(r['thresholds'],1):
                w.writerow([name,f'H{i}',h['V1_nM_per_s'],h['omega_per_s'],h['period_s'],h['crossing_slope']])
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
