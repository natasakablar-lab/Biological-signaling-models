"""Tests for analytical steady-state construction and spectral calculations."""
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.optimize import root
from model import *
from formulations import *
from equilibria import *


def main():
    p=Parameters();c=totals(initial_conditions());x=equilibrium_roots(p,c)[0]
    assert max(abs(rhs(0,x,p)))<1e-12
    assert max(abs(totals(x)-c))<1e-12
    assert np.all(x>=0)
    assert np.all(x[[0,2,4,6,20,21]]==0)
    ev=np.linalg.eigvals(reduced_jacobian(x[DIFFERENTIAL],c,p))
    assert ev.real.max()<0
    # Independent 16-D root solver from a perturbed equilibrium.
    q=x[DIFFERENTIAL]; seed=q.copy();seed[1]*=1.1
    nr=root(lambda q:reduced_rhs(0,q,c,p),seed,jac=lambda q:reduced_jacobian(q,c,p),tol=1e-11)
    assert np.max(np.abs(reconstruct(nr.x,c)-x))<1e-9
    # Complex-step scalar derivative checked against a central difference.
    h=x[9]*1e-4
    finite=(scalar_balance(x[9]+h,p,c[0])-scalar_balance(x[9]-h,p,c[0]))/(2*h)
    analytic=assess_equilibrium(x,p,c)['scalar_slope_at_root']
    assert np.isclose(finite,analytic,rtol=1e-6)
    # The equilibrium is independent of strictly positive k4 and k19.
    for key in ['k4','k19']:
        for factor in [.01,100.]:
            pp=replace(p,**{key:getattr(p,key)*factor})
            assert np.max(np.abs(equilibrium_roots(pp,c)[0]-x))<1e-12
    # Parameter-dependent Jacobian must be used, not the default Jacobian.
    pp=replace(p,k7=p.k7*10)
    xx=equilibrium_roots(pp,c)[0]
    good=reduced_jacobian(xx[DIFFERENTIAL],c,pp)
    eps=1e-25
    numerical=np.column_stack([np.imag(reduced_rhs(0,xx[DIFFERENTIAL].astype(complex)+1j*eps*np.eye(16)[i],c,pp))/eps for i in range(16)])
    assert np.max(np.abs(good-numerical))<1e-12
    cases=[(replace(p,k4=0),c),(p,np.array([.1,0,.1,.05,.05,.05]))]
    for pp,cc in cases:
        try: equilibrium_roots(pp,cc)
        except ValueError: pass
        else: raise AssertionError('Unsupported boundary case was accepted')
    scan=branch_scan('k11',points=11)
    assert len(scan['rows'])==11 and not scan['hopf_candidates']
    report={'baseline_residual_inf':float(max(abs(rhs(0,x,p)))),
      'baseline_max_real_eigenvalue':float(ev.real.max()),
      'independent_Newton_difference_uM':float(np.max(np.abs(reconstruct(nr.x,c)-x))),
      'variable_parameter_Jacobian_check':'passed','boundary_case_guards':'passed','small_branch_scan':'passed'}
    out=Path(__file__).parent/'results_equilibria';out.mkdir(exist_ok=True)
    (out/'verification_equilibria.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
