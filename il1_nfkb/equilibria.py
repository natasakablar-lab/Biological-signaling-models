"""Physical equilibria and natural-parameter branch tracking.

Scope: strictly positive rates and six strictly positive conserved totals.
No parameters or reactions are changed. No transient extrema are labeled
as equilibrium or periodic-orbit bifurcations.
"""
from dataclasses import asdict,replace
from functools import lru_cache
import numpy as np
from scipy.optimize import brentq,least_squares,linear_sum_assignment
from model import Parameters,initial_conditions,rhs,NAMES
from formulations import totals,DIFFERENTIAL,reduced_jacobian

SCAN_PARAMETERS=['k11','k7','k8','k6','a10','k13','a11','k19','k4']


def validate_analysis(p,c):
    if any(not np.isfinite(v) or v<=0 for v in asdict(p).values()):
        raise ValueError('Ova analiza ravnoteže zahteva strogo pozitivne kinetičke parametre. Nulte stope menjaju granične slučajeve i zakone očuvanja.')
    c=np.asarray(c,dtype=float)
    if c.shape!=(6,) or not np.isfinite(c).all() or (c<=0).any():
        raise ValueError('Potrebno je šest strogo pozitivnih bilansnih konstanti; slučaj bez IL-1 zahteva posebnu analizu granice.')
    return c


def downstream_at_n(n,p=Parameters()):
    """Closed-form stationary downstream states at given nuclear NF-kB n.

    Solve the IkB production/degradation quadratic using its positive root.
    Form is stable when the quadratic coefficient is small or zero.
    """
    L=p.a10*p.a12*n/((p.d12+p.k13)*p.d10+p.k13*p.a12*n)
    A2=p.a3*p.k13*L/(p.a11*(p.d3+p.k3))
    A1=(p.a3*p.d11*n/p.a11+p.k13*L)/(p.d3+p.k3)
    U=p.k8*(p.k6+p.k11*n*n)/p.k7
    aa=p.k3*A2; bb=p.k9+p.k3*A1
    b=2*U/(bb+np.sqrt(bb*bb+4*aa*U))
    j=L*b; free_nf=(p.d11*n+p.k13*j)/p.a11
    bound_nf=A2*b*b+A1*b
    x=np.zeros(22,dtype=np.result_type(n,float))
    x[[1,3,5,7,8,9,10]]=[b,free_nf,bound_nf,(p.k6+p.k11*n*n)/p.k7,
        p.a10*b/(p.d10+p.k13*p.a12*n/(p.d12+p.k13)),n,j]
    return x


def scalar_balance(n,p,c_nf):
    x=downstream_at_n(n,p)
    return x[[3,5,9,10]].sum()-c_nf


@lru_cache(maxsize=128)
def _upstream_cached(rates,conserved):
    c=np.array(conserved)
    K=np.array(rates)
    def build(logf):
        l,r,m,i,t=np.exp(logf)
        C14=K[0]*l*r; C16=K[1]*m*C14
        C18=K[2]*i*C16; C20=K[3]*t*C18
        return np.array([l,r,C14,m,C16,i,C18,t,C20])
    def residual(logf):
        l,r,A,m,B,i,D,t,F=build(logf)
        return (np.array([l+A+B+D+F,r+A+B+D+F,m+B+D+F,i+D+F,t+F])-c)/c
    sol=least_squares(residual,np.log(c/10),xtol=1e-13,ftol=1e-13,gtol=1e-13,max_nfev=800)
    if not sol.success or np.max(np.abs(residual(sol.x)))>1e-9:
        raise RuntimeError('Nije pronađena dovoljno tačna ravnoteža IL-1 vezivanja.')
    return build(sol.x)


def equilibrium_roots(p=Parameters(),c=None,grid_points=401):
    """Find sign-changing roots of the scalar balance; not a global root proof.

    The positive bracket guarantees at least one root. Additional roots and
    narrow/tangent branches can require dedicated continuation tools.
    """
    c=totals(initial_conditions()) if c is None else c
    c=validate_analysis(p,c)
    grid=np.unique(np.r_[0.,np.geomspace(c[0]*1e-12,c[0],grid_points),np.linspace(0,c[0],grid_points)])
    values=np.array([scalar_balance(n,p,c[0]) for n in grid])
    roots=[]
    for lo,hi,flo,fhi in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
        if flo*fhi<0:
            root=brentq(lambda n:scalar_balance(n,p,c[0]),lo,hi,xtol=5e-16,rtol=1e-13)
            if not roots or abs(root-roots[-1])>1e-12*c[0]: roots.append(root)
    if values[-1]==0: roots.append(grid[-1])
    if not roots: raise RuntimeError('Skalarna ravnoteža nije pronađena.')
    upstream=_upstream_cached(tuple([p.a14/p.d14,p.a15/p.d15,p.a16/p.d16,p.a17/p.d17]),tuple(c[1:]))
    equilibria=[]
    for n in roots:
        x=downstream_at_n(n,p); x[11:20]=upstream
        metrics=assess_equilibrium(x,p,c)
        if not metrics['physical'] or metrics['full_residual_inf']>1e-9 or metrics['conservation_error_inf']>1e-9:
            raise RuntimeError('Pronađeni kandidat nije prihvatljiva ravnoteža.')
        equilibria.append(x)
    return equilibria


def assess_equilibrium(x,p,c):
    eig=np.linalg.eigvals(reduced_jacobian(x[DIFFERENTIAL],c,p))
    leading=eig[np.argmax(eig.real)]
    alpha=float(leading.real)
    # Absolute guard for near-zero modes; near-neutral is not called Hopf.
    status='stabilna' if alpha < -1e-10 else 'nestabilna' if alpha > 1e-10 else 'blizu neutralnosti'
    slope=float(np.imag(scalar_balance(x[9]+1j*1e-25,p,c[0]))/1e-25)
    return {'physical':bool(np.min(x)>=-1e-12),
        'full_residual_inf':float(np.max(np.abs(rhs(0,x,p)))),
        'conservation_error_inf':float(np.max(np.abs(totals(x)-c))),
        'spectral_abscissa':alpha,'leading_real':float(leading.real),'leading_imag_abs':float(abs(leading.imag)),
        'stability':status,'scalar_slope_at_root':slope,
        'eigen_real':eig.real.tolist(),'eigen_imag':eig.imag.tolist()}


def branch_scan(parameter='k11',p=Parameters(),c=None,factor_min=.01,factor_max=100.,points=81):
    """Natural-parameter tracking of a detected physical equilibrium branch.

    Uses a scalar balance corrector at each positive parameter. Not pseudo-
    arclength continuation; no claim to exhaust all folds or hidden branches.
    """
    c=totals(initial_conditions()) if c is None else validate_analysis(p,c)
    validate_analysis(p,c)
    if parameter not in asdict(p): raise ValueError('Nepoznat parametar.')
    if factor_min<=0 or factor_max<=factor_min or not 5<=points<=401: raise ValueError('Neispravan opseg pretrage.')
    factors=np.geomspace(factor_min,factor_max,points)
    rows=[]; spectra=[]; states=[]; previous=None; previous_eig=None
    for factor in factors:
        value=getattr(p,parameter)*factor; pp=replace(p,**{parameter:float(value)})
        roots=equilibrium_roots(pp,c)
        x=min(roots,key=lambda xx:abs(xx[9]-(previous[9] if previous is not None else xx[9])))
        info=assess_equilibrium(x,pp,c)
        eig=np.array(info['eigen_real'])+1j*np.array(info['eigen_imag'])
        if previous_eig is not None:
            r,col=linear_sum_assignment(abs(previous_eig[:,None]-eig[None,:]))
            eig=eig[col]
        rows.append({'parameter':parameter,'value':float(value),'factor':float(factor),
          **{k:info[k] for k in ['spectral_abscissa','leading_real','leading_imag_abs','stability','full_residual_inf','conservation_error_inf','scalar_slope_at_root']},
          'detected_roots':len(roots),**{f'x{i+1}':float(v) for i,v in enumerate(x)}})
        previous=x; previous_eig=eig; states.append(x); spectra.append(eig)
    candidates=[]
    # Look for matched complex modes crossing the imaginary axis.
    for j in range(16):
        for i in range(points-1):
            a,b=spectra[i][j],spectra[i+1][j]
            if a.imag<=1e-8 or b.imag<=1e-8 or a.real*b.real>=0: continue
            target=.5*(a+b)
            def mode_at(logfactor):
                pp=replace(p,**{parameter:getattr(p,parameter)*np.exp(logfactor)})
                roots=equilibrium_roots(pp,c)
                xx=min(roots,key=lambda v:abs(v[9]-states[i][9]))
                ee=np.linalg.eigvals(reduced_jacobian(xx[DIFFERENTIAL],c,pp))
                e=min(ee,key=lambda v:abs(v-target))
                return e,xx,pp
            logroot=brentq(lambda z:mode_at(z)[0].real,np.log(factors[i]),np.log(factors[i+1]),xtol=1e-10)
            ee,xx,pp=mode_at(logroot); h=1e-4
            slope=(mode_at(logroot+h)[0].real-mode_at(logroot-h)[0].real)/(2*h)
            candidates.append({'parameter_value':float(getattr(pp,parameter)),
              'factor':float(np.exp(logroot)),'omega':float(abs(ee.imag)),
              'd_real_lambda_d_log_factor':float(slope),'label':'Hopf kandidat; nije potvrđen prvi Ljapunovljev koeficijent'})
    return {'parameter':parameter,'rows':rows,'spectra':np.array(spectra),
            'states':np.array(states),'hopf_candidates':candidates,
            'scope':'Praćenje pronađene fizičke grane; nije globalni dokaz odsustva bifurkacija.'}
