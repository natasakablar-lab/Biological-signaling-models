"""22-state IL-1--NF-kB model derived from supplied reactions, not printed ODEs.
Source: Paper With IL-1.pdf (31 July 2007) and four final-publication tables.
Parameter provenance remains open.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.integrate import solve_ivp

TITLE = 'Mathematical Model of IL-1-NFkB Biological Module'
NAMES = ['IKK_P','IkB','IKK_P:IkB','NFkB','IKK_P:IkB:NFkB','IkB:NFkB',
         'IKK_P:IkB:NFkB (x7)','IkBt','IkBn','NFkBn','IkBn:NFkBn',
         'IL-1','IL-R','C14 (IL-1:IL-R)','MyD88','C16 (MyD88:C14)',
         'IRAK','C18 (IRAK:C16)','TRAF6','C20 (TRAF6:C18)','IKK','C22 (IKK:C20)']
# x5 and x7 have different binding topology despite identical moiety composition.
DEFAULT_SELECTION = [0,1,9,11,20]
@dataclass(frozen=True)
class Parameters:
    a1: float=22.5e-3
    d1: float=1.25e-3
    a2: float=.5
    d2: float=.5e-3
    k2: float=2.04e-2
    a3: float=.5
    d3: float=.5e-3
    k3: float=2.25e-5
    a4: float=1.85e-1
    d4: float=1.25e-3
    k4: float=1.2e-4
    k5: float=4.07e-3
    k6: float=1.54e-6
    k7: float=2.8e-4
    k8: float=4.08e-3
    k9: float=1.13e-4
    a10: float=3e-4
    d10: float=2e-4
    a11: float=.9e-1
    d11: float=.8e-4
    a12: float=.5
    k11: float=1.65e-2
    d12: float=.5e-3
    k13: float=1.38e-2
    a14: float=.5
    d14: float=.5e-3
    a15: float=.5
    d15: float=.5e-3
    a16: float=.5
    d16: float=.5e-3
    a17: float=.5
    d17: float=.5e-3
    a18: float=.5
    d18: float=.5e-3
    k19: float=1.2e-4


def parameter_unit(name):
    if name=='k6': return 'μM/s'
    if name in ['a1','a2','a3','a4','a12','a14','a15','a16','a17','a18','k11']: return 'μM⁻¹ s⁻¹'
    return 's⁻¹'


def initial_conditions(preset='published'):
    y=np.zeros(22)
    y[[0,3,11,12,14,16,18,20]]=[.1,.1,.15,.1,.05,.05,.05,.1]
    if preset=='without_IL1': y[11]=0
    elif preset=='IL1_only': y[0]=0
    elif preset!='published': raise ValueError('Unknown preset')
    return y

# Reaction numbers refer to the updated PDF. R15 is interpreted as basal production.
# R21: 2*x10 -> x8 + 2*x10, hence no NF-kB consumption.
REACTIONS=[
 ('R1f',{11:1,12:1},{13:1},'a14'),('R1r',{13:1},{11:1,12:1},'d14'),
 ('R2f',{14:1,13:1},{15:1},'a15'),('R2r',{15:1},{14:1,13:1},'d15'),
 ('R3f',{16:1,15:1},{17:1},'a16'),('R3r',{17:1},{16:1,15:1},'d16'),
 ('R4f',{18:1,17:1},{19:1},'a17'),('R4r',{19:1},{18:1,17:1},'d17'),
 ('R5f',{20:1,19:1},{21:1},'a18'),('R5r',{21:1},{20:1,19:1},'d18'),
 ('R6',{21:1},{0:1,19:1},'k19'),
 ('R7f',{0:1,1:1},{2:1},'a1'),('R7r',{2:1},{0:1,1:1},'d1'),
 ('R8f',{2:1,3:1},{4:1},'a2'),('R8r',{4:1},{2:1,3:1},'d2'),
 ('R9',{4:1},{0:1,3:1},'k2'),
 ('R10f',{1:1,3:1},{5:1},'a3'),('R10r',{5:1},{1:1,3:1},'d3'),
 ('R11',{5:1},{3:1},'k3'),
 ('R12f',{0:1,5:1},{6:1},'a4'),('R12r',{6:1},{0:1,5:1},'d4'),
 ('R13',{0:1},{},'k4'),('R14',{2:1},{0:1},'k5'),
 ('R15',{}, {7:1},'k6'),('R16',{7:1},{},'k7'),
 ('R17',{7:1},{7:1,1:1},'k8'),('R18',{1:1},{},'k9'),
 ('R19f',{1:1},{8:1},'a10'),('R19r',{8:1},{1:1},'d10'),
 ('R20f',{3:1},{9:1},'a11'),('R20r',{9:1},{3:1},'d11'),
 ('R21',{9:2},{9:2,7:1},'k11'),
 ('R22f',{8:1,9:1},{10:1},'a12'),('R22r',{10:1},{8:1,9:1},'d12'),
 ('R23',{10:1},{5:1},'k13')]
S=np.zeros((22,len(REACTIONS)))
for j,(_,left,right,_) in enumerate(REACTIONS):
    for i,n in left.items(): S[i,j]-=n
    for i,n in right.items(): S[i,j]+=n

# Conserved moieties include every downstream complex containing that component.
CONSERVATION={
 'NFkB':[3,4,5,6,9,10],
 'IL-1':[11,13,15,17,19,21],
 'IL-R':[12,13,15,17,19,21],
 'MyD88':[14,15,17,19,21],
 'IRAK':[16,17,19,21],
 'TRAF6':[18,19,21]}
IKK_INDICES=[0,2,4,6,20,21]


def rates(y,p=Parameters()):
    v=[]
    for _,left,_,key in REACTIONS:
        rate=getattr(p,key)
        for i,n in left.items(): rate*=y[i]**n
        v.append(rate)
    return np.asarray(v)


def rhs(t,y,p=Parameters()): return S@rates(y,p)


def rhs_expanded(t,x,p=Parameters()):
    """Independent expanded equations for auditing the stoichiometric RHS."""
    x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22=x
    return np.array([
      -p.a1*x1*x2+p.d1*x3+p.k2*x5-p.a4*x1*x6+p.d4*x7-p.k4*x1+p.k5*x3+p.k19*x22,
      -p.a1*x1*x2+p.d1*x3-p.a3*x2*x4+p.d3*x6+p.k8*x8-p.k9*x2-p.a10*x2+p.d10*x9,
      p.a1*x1*x2-p.d1*x3-p.a2*x3*x4+p.d2*x5-p.k5*x3,
      -p.a2*x3*x4+p.d2*x5+p.k2*x5-p.a3*x2*x4+p.d3*x6+p.k3*x6-p.a11*x4+p.d11*x10,
      p.a2*x3*x4-p.d2*x5-p.k2*x5,
      p.a3*x2*x4-p.d3*x6-p.k3*x6-p.a4*x1*x6+p.d4*x7+p.k13*x11,
      p.a4*x1*x6-p.d4*x7,
      p.k6-p.k7*x8+p.k11*x10*x10,
      p.a10*x2-p.d10*x9-p.a12*x9*x10+p.d12*x11,
      p.a11*x4-p.d11*x10-p.a12*x9*x10+p.d12*x11,
      p.a12*x9*x10-p.d12*x11-p.k13*x11,
      -p.a14*x12*x13+p.d14*x14,
      -p.a14*x12*x13+p.d14*x14,
      p.a14*x12*x13-p.d14*x14-p.a15*x15*x14+p.d15*x16,
      -p.a15*x15*x14+p.d15*x16,
      p.a15*x15*x14-p.d15*x16-p.a16*x17*x16+p.d16*x18,
      -p.a16*x17*x16+p.d16*x18,
      p.a16*x17*x16-p.d16*x18-p.a17*x19*x18+p.d17*x20,
      -p.a17*x19*x18+p.d17*x20,
      p.a17*x19*x18-p.d17*x20-p.a18*x21*x20+p.d18*x22+p.k19*x22,
      -p.a18*x21*x20+p.d18*x22,
      p.a18*x21*x20-p.d18*x22-p.k19*x22])


def simulate(p=Parameters(),y0=None,t_end=50000.,points=2501,method='Radau',rtol=1e-8,atol=1e-11):
    y0=initial_conditions() if y0 is None else np.array(y0,dtype=float)
    if y0.shape!=(22,) or not np.isfinite(y0).all() or (y0<0).any():
        raise ValueError('Početni uslovi moraju biti 22 konačna nenegativna broja.')
    if any(not np.isfinite(v) or v<0 for v in asdict(p).values()):
        raise ValueError('Parametri moraju biti konačni i nenegativni.')
    if not np.isfinite(t_end) or t_end<=0 or not isinstance(points,int) or not 2<=points<=20001:
        raise ValueError('Neispravan vremenski interval ili broj tačaka (2–20001).')
    sol=solve_ivp(lambda t,y:rhs(t,y,p),(0,t_end),y0,
        t_eval=np.linspace(0,t_end,points),method=method,rtol=rtol,atol=atol)
    if not sol.success: raise RuntimeError(sol.message)
    return sol


def diagnostics(sol,p=Parameters()):
    drift={name:float(np.max(np.abs(sol.y[idx].sum(axis=0)-sol.y[idx,0].sum()))) for name,idx in CONSERVATION.items()}
    return {'minimum_concentration':float(sol.y.min()),'conservation_max_abs_drift':drift,
            'IKK_total_initial':float(sol.y[IKK_INDICES,0].sum()),
            'IKK_total_final':float(sol.y[IKK_INDICES,-1].sum())}


def symbolic_equations():
    """Generate exact symbolic RHS from the same reactions for the equation viewer."""
    import sympy as sp
    x=sp.symbols('x1:23'); ps={key:sp.Symbol(key) for key in asdict(Parameters())}
    v=[]
    for _,left,_,key in REACTIONS:
        term=ps[key]
        for i,n in left.items(): term*=x[i]**n
        v.append(term)
    return x, list(sp.Matrix(S.astype(int))*sp.Matrix(v))
