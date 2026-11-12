"""Equivalent full ODE, reduced ODE and genuine index-1 DAE formulations.

All public solutions return all 22 concentrations in the original order.
The DAE is solved by SUNDIALS IDA (scikit-sundae); there is no ODE fallback.
"""
from types import SimpleNamespace
from dataclasses import asdict
import numpy as np
from scipy.integrate import solve_ivp
from model import Parameters, initial_conditions, rhs, simulate, CONSERVATION, REACTIONS, S

ALGEBRAIC=np.array([3,11,12,14,16,18],dtype=int)  # x4,x12,x13,x15,x17,x19
DIFFERENTIAL=np.array([i for i in range(22) if i not in ALGEBRAIC],dtype=int)
C=np.zeros((6,22))
for row,idx in enumerate(CONSERVATION.values()): C[row,idx]=1.
assert np.array_equal(C[:,ALGEBRAIC],np.eye(6))
# x = B*q + b(c), where q contains the 16 differential coordinates.
B=np.zeros((22,16)); B[DIFFERENTIAL]=np.eye(16); B[ALGEBRAIC]=-C[:,DIFFERENTIAL]
E=np.diag([0. if i in ALGEBRAIC else 1. for i in range(22)])


def totals(x): return C@np.asarray(x,dtype=float)


def reconstruct(q,c):
    """Reconstruct 22 states from 16 coordinates and six conserved totals.

    Accepts q with shape (16,) or (16,N); no clipping of concentrations.
    """
    q=np.asarray(q); c=np.asarray(c)
    if q.shape[0]!=16 or c.shape!=(6,): raise ValueError('Expected q[16,...], c[6].')
    out=B@q
    if q.ndim==1: out[ALGEBRAIC]+=c
    elif q.ndim==2: out[ALGEBRAIC]+=c[:,None]
    else: raise ValueError('q must have one or two dimensions.')
    return out


def reduced_rhs(t,q,c,p=Parameters()):
    return rhs(t,reconstruct(q,c),p)[DIFFERENTIAL]


def full_jacobian(x,p=Parameters()):
    """Exact mass-action Jacobian, valid even at zero concentrations."""
    dv=np.zeros((len(REACTIONS),22))
    for j,(_,left,_,key) in enumerate(REACTIONS):
        for i,n in left.items():
            term=getattr(p,key)*n
            for h,m in left.items(): term*=x[h]**(m-1 if h==i else m)
            dv[j,i]=term
    return S@dv


def reduced_jacobian(q,c,p=Parameters()):
    return full_jacobian(reconstruct(q,c),p)[DIFFERENTIAL]@B


def dae_residual(t,x,xp,c,p=Parameters()):
    """F=0 in original variable order: 16 dynamics and six constraints."""
    res=np.empty(22)
    res[DIFFERENTIAL]=np.asarray(xp)[DIFFERENTIAL]-rhs(t,x,p)[DIFFERENTIAL]
    res[ALGEBRAIC]=totals(x)-c
    return res


def consistent_initial_derivative(x0,p=Parameters()):
    """Full ODE derivative also satisfies differentiated conservation laws."""
    return rhs(0,np.asarray(x0),p)


def ida_available():
    try:
        from sksundae.ida import IDA
        return True
    except ImportError: return False


def validate(p,x0,t_end,points):
    x0=np.asarray(x0,dtype=float)
    if x0.shape!=(22,) or not np.isfinite(x0).all() or (x0<0).any():
        raise ValueError('Potrebna su 22 konačna nenegativna početna uslova.')
    if any(not np.isfinite(v) or v<0 for v in asdict(p).values()):
        raise ValueError('Parametri moraju biti konačni i nenegativni.')
    if not np.isfinite(t_end) or t_end<=0 or not isinstance(points,int) or not 2<=points<=20001:
        raise ValueError('Neispravan vremenski interval ili broj tačaka.')
    return x0


def simulate_formulation(formulation='full',p=Parameters(),x0=None,t_end=50000.,points=2501,
                         method='Radau',rtol=1e-9,atol=1e-12,conserved_totals=None):
    x0=initial_conditions() if x0 is None else x0
    x0=validate(p,x0,t_end,points)
    c=totals(x0) if conserved_totals is None else np.asarray(conserved_totals,dtype=float)
    if c.shape!=(6,) or not np.isfinite(c).all(): raise ValueError('Neispravni konstantni bilansi.')
    if np.max(np.abs(totals(x0)-c))>1e-12:
        raise ValueError('Početni uslovi nisu saglasni sa zadatim konstantnim bilansima.')
    ts=np.linspace(0,t_end,points)
    if formulation=='full':
        sol=simulate(p,x0,t_end,points,method,rtol,atol)
        sol.formulation='full'; sol.coordinates=22; sol.solver_name=method; sol.totals=c
        return sol
    if formulation=='reduced':
        opts={}
        if method in ['Radau','BDF','LSODA']:
            opts['jac']=lambda t,q:reduced_jacobian(q,c,p)
        sol=solve_ivp(lambda t,q:reduced_rhs(t,q,c,p),(0,t_end),x0[DIFFERENTIAL],
            t_eval=ts,method=method,rtol=rtol,atol=atol,**opts)
        if not sol.success: raise RuntimeError(sol.message)
        sol.q=sol.y.copy(); sol.y=reconstruct(sol.q,c)
        sol.formulation='reduced'; sol.coordinates=16; sol.solver_name=method; sol.totals=c
        return sol
    if formulation!='dae': raise ValueError('Unknown formulation')
    try:
        from sksundae.ida import IDA
    except ImportError as exc:
        raise ImportError('Za pravi DAE rešavač instalirajte: pip install scikit-sundae') from exc
    def resfn(t,x,xp,res): res[:]=dae_residual(t,x,xp,c,p)
    solver=IDA(resfn,algebraic_idx=ALGEBRAIC.tolist(),rtol=rtol,atol=atol)
    out=solver.solve(ts,x0,consistent_initial_derivative(x0,p))
    if not out.success: raise RuntimeError(out.message)
    return SimpleNamespace(t=out.t,y=out.y.T,yp=out.yp.T,success=True,
        message=out.message,formulation='dae',coordinates=16,solver_name='SUNDIALS IDA',totals=c)


def comparison_table(solutions):
    """Max discrepancy against full ODE, on an identical time grid."""
    reference=solutions['full']
    rows=[]
    for name,sol in solutions.items():
        if sol.y.shape!=reference.y.shape or not np.array_equal(sol.t,reference.t):
            raise ValueError('Comparison requires identical time grids')
        delta=np.abs(sol.y-reference.y)
        for i in range(22):
            rows.append({'formulation':name,'state':f'x{i+1}',
                'max_abs_error_uM':float(delta[i].max()),
                'rms_error_uM':float(np.sqrt(np.mean(delta[i]**2)))})
    return rows
