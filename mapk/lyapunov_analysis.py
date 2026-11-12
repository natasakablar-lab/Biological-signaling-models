"""First Lyapunov coefficient, Kuznetsov convention, for the 5D MAPK model.

Run: python lyapunov_analysis.py (requires SymPy in addition to version 2).
Coordinates in nM, time in seconds, ||q||_2=1 and <p,q>=1.
Formula: Kuznetsov, Elements of Applied Bifurcation Theory, 2nd ed., (5.62).
"""
from pathlib import Path
from dataclasses import asdict, replace
from types import SimpleNamespace
from itertools import combinations_with_replacement, permutations
import json
import csv
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from mapk_model import Parameters, KEEP, reduced_rhs, jacobian, equilibrium


def derivative_builder(parameters):
    """Differentiate the SAME RHS used in simulations, to order three."""
    z=sp.symbols('z0:5',real=True)
    v=sp.Symbol('V1',real=True)
    p=SimpleNamespace(**{k:sp.Rational(str(val)) for k,val in asdict(parameters).items()})
    p.V1=v
    f=list(reduced_rhs(0,np.array(z,dtype=object),p))
    n=len(z)
    keys=[]
    expressions=[]
    for order in (1,2,3):
        for i in range(n):
            for inds in combinations_with_replacement(range(n),order):
                keys.append((order,i,inds))
                expressions.append(sp.diff(f[i],*[z[j] for j in inds]))
    numeric=sp.lambdify((*z,v),expressions,'numpy',cse=True)
    def evaluate(x,V1):
        tensors={d:np.zeros((n,)*(d+1)) for d in (1,2,3)}
        for (order,i,inds),value in zip(keys,numeric(*x,V1)):
            for perm in set(permutations(inds)):
                tensors[order][(i,*perm)]=value
        return tensors[1],tensors[2],tensors[3]
    return evaluate


def coefficient(A,B,C, phase=0.):
    """B=D2f and C=D3f, with no factorials inside the tensors.

    f(x0+u)=A u + B(u,u)/2 + C(u,u,u)/6 + ... .
    Complex scalar product is np.vdot(p,q)=conj(p).T @ q.
    """
    vals,vecs=np.linalg.eig(A)
    choices=np.where(vals.imag>1e-7)[0]
    idx=choices[np.argmin(abs(vals[choices].real))]
    omega=float(vals[idx].imag)
    q=vecs[:,idx]/np.linalg.norm(vecs[:,idx])*np.exp(1j*phase)
    left,lvec=np.linalg.eig(A.T)
    p=lvec[:,np.argmin(abs(left+1j*omega))]
    p=p/np.conj(np.vdot(p,q))
    b=lambda u,v: np.einsum('ijk,j,k->i',B,u,v)
    c=lambda u,v,w:np.einsum('ijkl,j,k,l->i',C,u,v,w)
    h11=-np.linalg.solve(A,b(q,q.conj()))
    h20=np.linalg.solve(2j*omega*np.eye(len(A))-A,b(q,q))
    terms=np.array([np.vdot(p,c(q,q,q.conj())),
                    2*np.vdot(p,b(q,h11)), np.vdot(p,b(q.conj(),h20))])
    G21=sum(terms)
    l1=float(G21.real/(2*omega))
    checks={'q_norm':float(np.linalg.norm(q)),
            'p_q_error':float(abs(np.vdot(p,q)-1)),
            'right_eigenvector_residual':float(np.linalg.norm(A@q-1j*omega*q)),
            'left_eigenvector_residual':float(np.linalg.norm(A.T@p+1j*omega*p))}
    assert checks['p_q_error']<1e-10
    assert checks['right_eigenvector_residual']<1e-9
    assert checks['left_eigenvector_residual']<1e-9
    return l1,G21,terms,checks


def finite_tensors(x,p,step):
    """Independent derivatives of numerical RHS via complex-step Jacobian.

    B: centered first differences of J. C: centered mixed second differences
    of J, including the repeated-index stencil at offsets +/-2*step.
    """
    n=len(x)
    E=np.eye(n)
    fun=lambda u:reduced_rhs(0,u,p)
    J=lambda u:jacobian(fun,u)
    A=J(x)
    B=np.empty((n,n,n));C=np.empty((n,n,n,n))
    for k in range(n):
        B[:,:,k]=(J(x+step*E[k])-J(x-step*E[k]))/(2*step)
        for l in range(n):
            C[:,:,k,l]=(J(x+step*E[k]+step*E[l])-J(x+step*E[k]-step*E[l])
                       -J(x-step*E[k]+step*E[l])+J(x-step*E[k]-step*E[l]))/(4*step**2)
    return A,B,C


def benchmark():
    """Exact cubic radial oscillator: expected l1=2*a/omega for ||q||=1."""
    x,y=sp.symbols('x y',real=True)
    tests=[]
    for a in (-1.,1.):
        omega=1.7
        f=sp.Matrix([-omega*y+a*x*(x*x+y*y),omega*x+a*y*(x*x+y*y)])
        zero={x:0,y:0};z=(x,y)
        A=np.array(f.jacobian(z).subs(zero),float)
        B=np.array([[[sp.diff(fi,z[j],z[k]).subs(zero) for k in range(2)]
                     for j in range(2)] for fi in f],float)
        C=np.array([[[[sp.diff(fi,z[j],z[k],z[l]).subs(zero) for l in range(2)]
                      for k in range(2)] for j in range(2)] for fi in f],float)
        l1,*_=coefficient(A,B,C)
        assert abs(l1-2*a/omega)<1e-12
        tests.append({'a':a,'omega':omega,'expected':2*a/omega,'computed':l1})
    # A 3D example checks quadratic terms absent from the radial benchmark.
    # Kuznetsov 2nd ed. Example 5.4, omega=1, q=(1,i,-1).
    # With unit-norm q, published l1=-0.9 must be divided by ||q||^2=3.
    A=np.array([[0.,1.,0.],[0.,0.,1.],[-1.,-1.,-1.]])
    B=np.zeros((3,3,3));B[2,0,0]=2.
    C=np.zeros((3,3,3,3))
    l1,*_=coefficient(A,B,C)
    assert abs(l1-(-.3))<1e-12
    tests.append({'example':'Kuznetsov Example 5.4, omega=1, unit q',
                  'expected':-.3,'computed':l1})
    return tests


def main():
    root=Path(__file__).resolve().parent
    source=root/'results_hopf'/'hopf_summary.json'
    if not source.exists():
        raise FileNotFoundError('Run hopf_analysis.py first, then lyapunov_analysis.py')
    points=json.loads(source.read_text())
    out=root/'results_lyapunov';out.mkdir(exist_ok=True)
    report={'convention':'Kuznetsov l1=Re(G21)/(2*omega), raw nM coordinates, unit Euclidean q, <p,q>=1',
            'formula_source':'https://www.ma.ic.ac.uk/~dturaev/kuznetsov.pdf, Eq. (5.62), Example 5.4',
            'benchmark':benchmark(),'cases':{}}
    rows=[]
    for name,case in points.items():
        base=Parameters(**case['parameters'])
        print(f'Building symbolic derivatives: {name}',flush=True)
        derivatives=derivative_builder(base)
        for index,point in enumerate(case['thresholds'],1):
            guess=np.array(point['equilibrium_nM'])[KEEP]
            # Refine the saved Hopf point before evaluating higher derivatives.
            v0=point['V1_nM_per_s']
            value=brentq(lambda v:max(equilibrium(replace(base,V1=v),guess)[2].real),
                          v0*.999,v0*1.001,xtol=1e-13)
            p=replace(base,V1=value)
            y,J,eig,residual=equilibrium(p,guess)
            x=y[KEEP]
            A,B,C=derivatives(x,value)
            assert np.max(abs(A-J))<1e-11
            l1,G21,terms,checks=coefficient(A,B,C)
            phase=coefficient(A,B,C,phase=.731)[0]
            assert abs(phase-l1)<1e-12
            finite=[]
            for h in (.04,.02,.01,.005):
                approx=coefficient(*finite_tensors(x,p,h))[0]
                finite.append({'step_nM':h,'l1':approx,'relative_error':abs(approx-l1)/abs(l1)})
            assert max(v['relative_error'] for v in finite[-2:])<1e-4
            key=f'{name}_H{index}'
            answer={'V1':value,'l1':l1,'G21':[float(G21.real),float(G21.imag)],
                    'cubic_coefficient_original_time_ReG21_over_2':float(G21.real/2),
                    'classification':'supercritical' if l1<0 else 'subcritical',
                    'stable_small_cycle_side_of_V1':('above' if point['crossing_slope']>0 else 'below') if l1<0 else None,
                    'crossing_slope':point['crossing_slope'],
                    'equilibrium_residual':residual,'normalization_checks':checks,
                    'phase_invariance_error':abs(phase-l1),'finite_difference_checks':finite,
                    'G21_terms':[[float(t.real),float(t.imag)] for t in terms]}
            report['cases'][key]=answer
            rows.append([name,f'H{index}',value,l1,G21.real/2,answer['classification'],
                         answer['stable_small_cycle_side_of_V1']])
            np.savez(out/f'{key}_derivatives.npz',A=A,B=B,C=C,x=x,V1=value)
            print(key,l1,answer['classification'],flush=True)
    (out/'lyapunov_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    with (out/'lyapunov_coefficients.csv').open('w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['scenario','point','V1','l1_Kuznetsov','ReG21_over_2_original_time',
                    'classification','stable_small_cycle_side_of_V1'])
        w.writerows(rows)


if __name__=='__main__':main()
