"""MAPK reconstruction from N. A. Kablar, MTNS 2010, pp. 2363-2368.

Run: python mapk_model.py
Requires: numpy, scipy, matplotlib. Outputs are saved beside this script.
Concentrations: nM; time: seconds. No MATLAB dependency.
This is a phenomenological Michaelis-Menten model, not an explicit
enzyme-complex mass-action model. Assumed parameters are not fitted data.
"""
from dataclasses import dataclass, asdict, replace
from pathlib import Path
import csv
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAMES = ("KKK", "KKKP", "KK", "KKP", "KKPP", "K", "KP", "KPP")
KEEP = np.array([0, 2, 4, 5, 7])
TOTALS = np.array([100., 300., 300.])
Y0 = np.array([100., 0., 300., 0., 0., 300., 0., 0.])


@dataclass(frozen=True)
class Parameters:
    V1: float = 0.25
    V2: float = 0.25
    n: float = 1.
    KI: float = 9.
    KM1: float = 10.
    KM2: float = 8.
    k3: float = 0.025
    k4: float = 0.025
    k5: float = 0.75
    k6: float = 0.75
    k7: float = 0.025
    k8: float = 0.025
    k9: float = 0.5
    k10: float = 0.5
    KM3: float = 15.
    KM4: float = 15.
    KM5: float = 15.
    KM6: float = 15.
    KM7: float = 15.
    KM8: float = 15.
    KM9: float = 15.
    KM10: float = 15.
    # ASSUMED: phosphatase totals, not explicitly tabulated in the paper.
    phosphatase_KK: float = 1.
    phosphatase_K: float = 1.
    # Reference case: additional KKP-mediated reactions switched off.
    k11: float = 0.
    k12: float = 0.
    # ASSUMED by analogy with KM7 and KM8 (inactive when k11=k12=0).
    KM11: float = 15.
    KM12: float = 15.
    # Cross-talk and positive feedback switched off until separately studied.
    # Va,Vd are effective Vmax values (nM/s), not ka1,kd1 (1/s).
    Va: float = 0.
    Vd: float = 0.
    KMa: float = 10.
    KMd: float = 10.
    kp: float = 0.
    KMp: float = 10.  # ASSUMED, inactive while kp=0.

    def validate(self):
        for key, value in asdict(self).items():
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{key} must be finite and nonnegative")
            if (key.startswith("KM") or key in ("KI", "n")) and value <= 0:
                raise ValueError(f"{key} must be positive")


def rates(y, p):
    a, b, c, d, e, f, g, h = y
    v1 = p.V1*a / ((1+(h/p.KI)**p.n)*(p.KM1+a))
    v2 = p.V2*b/(p.KM2+b)
    v3 = p.k3*b*c/(p.KM3+c)
    v4 = p.k4*b*d/(p.KM4+d)
    v5 = p.k5*p.phosphatase_KK*e/(p.KM5+e)
    v6 = p.k6*p.phosphatase_KK*d/(p.KM6+d)
    v7 = p.k7*e*f/(p.KM7+f)
    v8 = p.k8*e*g/(p.KM8+g)
    v9 = p.k9*p.phosphatase_K*h/(p.KM9+h)
    v10 = p.k10*p.phosphatase_K*g/(p.KM10+g)
    v11 = p.k11*d*f/(p.KM11+f)
    v12 = p.k12*d*g/(p.KM12+g)
    va = p.Va*a/(p.KMa+a)
    vd = p.Vd*b/(p.KMd+b)
    vp = p.kp*h*a/(p.KMp+a)
    return np.array([v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,va,vd,vp])


def rhs(t, y, p):
    v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,va,vd,vp = rates(y,p)
    q = v2-v1-va+vd-vp
    return np.array([q, -q, v6-v3, v3+v5-v4-v6, v4-v5,
                     v10-v7-v11, v7+v9-v8-v10+v11-v12, v8-v9+v12])


def expand(z, totals=TOTALS):
    """Five independent states: KKK, KK, KKPP, K, KPP."""
    a,c,e,f,h = z
    return np.array([a, totals[0]-a, c, totals[1]-c-e,
                     e, f, totals[2]-f-h, h])


def reduced_rhs(t, z, p):
    return rhs(t, expand(z), p)[KEEP]


def jacobian(fun, x):
    """Complex-step Jacobian; no clipping/nonanalytic branches in rates."""
    x = np.asarray(x, dtype=float)
    result = np.empty((len(x),len(x)))
    for j in range(len(x)):
        z = x.astype(complex)
        z[j] += 1e-25j
        result[:,j] = np.imag(fun(z))/1e-25
    return result


def simulate(p, end=12000., method="DOP853", rtol=1e-9, atol=1e-11):
    p.validate()
    sol = solve_ivp(lambda t,y: rhs(t,y,p), (0,end), Y0,
                    method=method, t_eval=np.linspace(0,end,6001),
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    sums = np.vstack((sol.y[:2].sum(0), sol.y[2:5].sum(0), sol.y[5:].sum(0)))
    if sol.y.min() < -1e-7 or np.max(np.abs(sums-TOTALS[:,None])) > 1e-6:
        raise RuntimeError("Physical checks failed; concentrations were NOT clipped")
    return sol


def equilibrium(p, guess):
    fun = lambda z: reduced_rhs(0,z,p)
    result = root(fun, guess, jac=lambda z: jacobian(fun,z), tol=1e-10)
    y = expand(result.x)
    residual = float(np.max(np.abs(rhs(0,y,p))))
    if residual > 1e-8 or y.min() < -1e-7:
        raise RuntimeError(f"No admissible equilibrium found: {result.message}")
    J = jacobian(fun, result.x)
    eig = np.linalg.eigvals(J)
    return y, J, eig, residual


def period_estimate(sol):
    """Late-time peak statistics, descriptive only; not proof of a limit cycle."""
    start = int(0.65*len(sol.t))
    signal = sol.y[7,start:]
    span = np.ptp(signal)
    if span < 1e-3:
        return {"late_KPP_span_nM": float(span), "peak_period_s": None}
    idx, _ = find_peaks(signal, prominence=max(1e-3, .05*span))
    intervals = np.diff(sol.t[start:][idx])
    return {"late_KPP_span_nM": float(span),
            "peak_period_s": float(np.mean(intervals)) if len(intervals)>=2 else None,
            "peak_intervals_count": int(len(intervals))}


def main():
    out = Path(__file__).resolve().parent / "results"
    out.mkdir(exist_ok=True)
    base = Parameters()
    # Illustrative weak catalytic activity: 10% of k7=k8.
    # This ratio is a modelling assumption, NOT an experimental estimate.
    weak = replace(base, k11=0.0025, k12=0.0025)
    cases = {"reference_V025": base, "reference_V5": replace(base,V1=5.),
             "estimated_V025": weak, "estimated_V5": replace(weak,V1=5.)}
    report = {}
    plt.rcParams.update({"font.size":10, "axes.spines.top":False,
                         "axes.spines.right":False})
    fig, axs = plt.subplots(2,2,figsize=(12,7),layout="constrained")
    labels = {"reference_V025":"Reference: V1=0.25", "reference_V5":"Reference: V1=5",
              "estimated_V025":"Assumed k11=k12=0.0025: V1=0.25",
              "estimated_V5":"Assumed k11=k12=0.0025: V1=5"}
    for ax,(name,p) in zip(axs.flat,cases.items()):
        sol = simulate(p)
        eq,J,eig,residual = equilibrium(p, sol.y[KEEP,-1])
        # Independent stiff solver and five-state formulation verification.
        alt = simulate(p, method="Radau",rtol=2e-10,atol=2e-12)
        solver_delta = float(np.max(np.abs(sol.y-alt.y)))
        red = solve_ivp(lambda t,z: reduced_rhs(t,z,p),(0,sol.t[-1]),Y0[KEEP],
                        t_eval=sol.t,method="DOP853",rtol=1e-9,atol=1e-11)
        if not red.success:
            raise RuntimeError(red.message)
        reduction_delta = float(np.max(np.abs(sol.y-np.array([expand(z) for z in red.y.T]).T)))
        # Central differences independently check the complex-step derivative.
        x = eq[KEEP]
        finite = np.column_stack([(reduced_rhs(0,x+np.eye(5)[j]*1e-4,p)
                                 -reduced_rhs(0,x-np.eye(5)[j]*1e-4,p))/2e-4
                                 for j in range(5)])
        j_error = float(np.max(np.abs(J-finite)))
        assert solver_delta < 2e-3, (name,solver_delta)
        assert reduction_delta < 2e-3, (name,reduction_delta)
        assert j_error < 1e-6, (name,j_error)
        full_eig = np.linalg.eigvals(jacobian(lambda y:rhs(0,y,p),eq))
        assert np.sum(np.abs(full_eig)<1e-9)==3
        sums = np.vstack((sol.y[:2].sum(0),sol.y[2:5].sum(0),sol.y[5:].sum(0)))
        spectral = float(max(eig.real))
        stability = ("locally asymptotically stable" if spectral < -1e-8
                     else "unstable" if spectral > 1e-8 else "near nonhyperbolic")
        report[name] = {"parameters":asdict(p), "equilibrium_nM":dict(zip(NAMES,eq.tolist())),
                        "rhs_residual_nM_per_s":residual,
                        "eigenvalues_5D_per_s":[[float(v.real),float(v.imag)] for v in eig],
                        "spectral_abscissa_per_s":spectral,"equilibrium_stability":stability,
                        "minimum_concentration_nM":float(sol.y.min()),
                        "max_total_drift_nM":float(np.max(abs(sums-TOTALS[:,None]))),
                        "DOP853_Radau_max_difference_nM":solver_delta,
                        "eight_five_state_max_difference_nM":reduction_delta,
                        "Jacobian_central_difference_error":j_error, **period_estimate(sol)}
        np.savetxt(out/f"{name}_trajectory.csv",np.column_stack((sol.t,sol.y.T)),
                   delimiter=",",header="time_s,"+",".join(NAMES),comments="")
        np.savetxt(out/f"{name}_jacobian_5D.csv",J,delimiter=",")
        for i,label in ((1,"KKKP"),(4,"KKPP"),(7,"KPP")):
            ax.plot(sol.t/60,sol.y[i],label=label,lw=1.15)
        ax.set(title=labels[name],xlabel="Time (min)",ylabel="Concentration (nM)")
        ax.grid(alpha=.2)
        ax.legend(loc="upper right",fontsize=8)
    fig.savefig(out/"time_responses.png",dpi=180)
    plt.close(fig)
    fig,axs = plt.subplots(1,2,figsize=(10,4),layout="constrained")
    for ax,(name,p) in zip(axs,[list(cases.items())[1],list(cases.items())[3]]):
        sol = simulate(p)
        late = sol.t > .65*sol.t[-1]
        ax.plot(sol.y[4,late],sol.y[7,late],lw=1.2)
        ax.set(title=labels[name],xlabel="KKPP (nM)",ylabel="KPP (nM)")
        ax.grid(alpha=.2)
    fig.savefig(out/"late_phase_portraits.png",dpi=180)
    plt.close(fig)
    # Paper's printed equilibrium, same state order as the five-state system.
    printed = np.array([81.938,272.22,7.8503,288.5,2.5738])
    actual = np.array(list(report["reference_V025"]["equilibrium_nM"].values()))[KEEP]
    report["paper_comparison"] = {"state_order":[NAMES[i] for i in KEEP],
        "paper_equilibrium_nM":printed.tolist(),
        "reconstructed_equilibrium_nM":actual.tolist(),
        "difference_nM":(actual-printed).tolist(),
        "note":"k11=k12=0, phosphatase totals=1 nM are reconstruction assumptions."}
    (out/"summary.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    with (out/"equilibria.csv").open("w",newline="",encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["scenario",*NAMES,"spectral_abscissa_per_s","stability"])
        for name in cases:
            r = report[name]
            writer.writerow([name,*r["equilibrium_nM"].values(),r["spectral_abscissa_per_s"],r["equilibrium_stability"]])
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
