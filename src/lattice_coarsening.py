# -*- coding: utf-8 -*-
"""
lattice_coarsening.py  --  Paper 3, Section 6

Direct coarsening simulation in the lattice model, to supply the absolute
length scale that the phase-field calculation cannot.

WHY THIS RUN EXISTS

Section 2.5 explains that the phase-field interface must be artificially widened
for the solver to resolve it, so the domain sizes appearing there are set by the
grid rather than predicted. Here the interface is what it physically is -- one
bond wide -- and no inflation occurs. What this calculation gives up in reach it
gains in meaning: sizes come out in atoms, and hence in nanometres.

WHAT IS MEASURED

Starting from a dispersed solute distribution, the system is held at a fixed
temperature and the domain population is followed as it coarsens. The
observables are the number of clusters, the mean and largest cluster size, and
the radius of gyration of the largest cluster, all as functions of Monte Carlo
sweep. The coarsening exponent is extracted from the growth of the mean cluster
radius.

Because Monte Carlo sweeps are not physical time, the exponent is the
transferable quantity, not the rate. It is the quantity that Section 7
cross-checks against the phase-field result, and the two use interfaces that
differ by more than an order of magnitude in width, so agreement between them
tests something.

TEMPERATURE

Coarsening at 573 K is unobservably slow: the acceptance ratio there is of order
1e-5 and no coarsening occurs within any affordable run. We therefore work at
1400-1800 K, where the chain is mobile, and rely on the exponent being a
property of the transport-limited kinetics rather than of the temperature. That
assumption is itself tested, by running at several temperatures and checking
that the exponent does not move.

    python lattice_coarsening.py                 # default sweep
    python lattice_coarsening.py --quick         # small, fast check
"""

from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

try:
    from sphalerite_mc import SphaleriteLattice, HamiltonianParams, KB_EV
    from samplers import MixedMoveMC
except ModuleNotFoundError as err:
    here = os.path.dirname(os.path.abspath(__file__))
    print(f"Cannot import '{err.name}'. This script needs sphalerite_mc.py and "
          f"samplers.py in the same folder:\n  {here}")
    raise SystemExit(1)

A0 = 5.4093          # ZnS lattice parameter, Angstrom
_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, os.pardir, "data")
_FIG = os.path.join(_HERE, os.pardir, "figures")


# ======================================================================
def cluster_stats(mc, lat, min_size=4):
    """
    Domain population of the current configuration.

    Sizes are converted to physical radii by treating each cluster as compact:
    a cluster of n cation sites occupies n * Omega with Omega = a0^3/4 the
    volume per cation site, so R = (3 n Omega / 4 pi)^(1/3).

    MONOMERS AND SMALL CLUSTERS ARE EXCLUDED from the mean, and min_size sets
    the threshold. This is not cosmetic. At 1800 K some 90 per cent of clusters
    are single solutes evaporated from the domains, and including them makes the
    "mean radius" a measure of the free-solute concentration rather than of the
    domains: it FALLS while the domains grow, because a growing domain sheds
    more monomers. An earlier version did include them and returned coarsening
    exponents of 0.007 to 0.165 across three temperatures, which is not a
    coarsening exponent but an artefact of that mixing.

    Both quantities are reported: mean_radius_nm over clusters of at least
    min_size sites, and monomer_frac, so that the evaporation is visible rather
    than hidden.
    """
    cs = mc.cluster_sizes()
    omega = A0 ** 3 / 4.0                       # Angstrom^3 per cation site
    big = cs[cs >= min_size]
    r_big = (3.0 * big * omega / (4.0 * np.pi)) ** (1.0 / 3.0)   # Angstrom
    r_all = (3.0 * cs * omega / (4.0 * np.pi)) ** (1.0 / 3.0)
    n_sol = int(cs.sum())
    return dict(n_clusters=int(cs.size),
                n_domains=int(big.size),
                n_solute=n_sol,
                max_size=int(cs[0]) if cs.size else 0,
                mean_size=float(big.mean()) if big.size else 0.0,
                mean_radius_nm=float(r_big.mean() / 10.0) if big.size else 0.0,
                mean_radius_all_nm=float(r_all.mean() / 10.0) if cs.size else 0.0,
                max_radius_nm=float(r_all.max() / 10.0) if cs.size else 0.0,
                monomer_frac=float(np.count_nonzero(cs == 1) / cs.size)
                             if cs.size else 0.0,
                frac_in_largest=float(cs[0] / n_sol) if n_sol else 0.0,
                solute_in_domains=float(big.sum() / n_sol) if n_sol else 0.0,
                sizes=[int(v) for v in cs[:20]])


def run_coarsening(L=30, x=0.04, T=1600.0, sweeps=40000, n_record=60,
                   seed=0, verbose=True):
    """
    Hold at fixed T from a dispersed start and follow the domain population.

    x is the total solute fraction, x/2 each of Cu and In.

    THE COMPOSITION MUST BE HIGH ENOUGH TO SUSTAIN COARSENING. A first attempt
    used x = 0.01, which put 1080 solutes in the box; they condensed into a
    single domain within 3000 of the 60,000 sweeps, so more than nine tenths of
    the run measured a structure that had stopped evolving and the fitted
    exponent was dominated by that flat tail. The equilibrium state IS a single
    domain -- that is the result of the preceding paper -- so a coarsening
    measurement must be made before the solute budget is exhausted, which means
    supplying enough solute that exhaustion takes a long time. x = 0.04 gives
    4320 solutes at L = 30 and extends the coarsening window by more than an
    order of magnitude.

    This composition is far above natural indium tenors. It is chosen to make
    the COARSENING EXPONENT measurable, and the exponent is a property of the
    transport-limited kinetics rather than of the concentration; the absolute
    sizes reported are those of the simulated domains, not a prediction of
    natural inclusion sizes.
    """
    lat = SphaleriteLattice(L)
    lat._lut = None
    par = HamiltonianParams()
    mc = MixedMoveMC(lat, par, x_cu=x / 2, x_in=x / 2, seed=seed, p_ss=0.5)

    box_nm = L * A0 / 10.0
    if verbose:
        print(f"  L={L} (N={lat.N:,} sites, {box_nm:.1f} nm cube), "
              f"x={x*100:.2f}%, T={T:.0f} K")
        print(f"    {int(x*lat.N):,} solutes, {sweeps:,} sweeps")

    # logarithmic recording: coarsening is a power law, so uniform spacing
    # wastes most of the points on the late, slow part
    marks = np.unique(np.logspace(0, np.log10(sweeps), n_record).astype(int))
    trace = []
    done = 0
    t0 = time.time()
    st = cluster_stats(mc, lat)
    st.update(sweep=0, acc=np.nan)
    trace.append(st)
    for m in marks:
        acc = mc.run(T, int(m - done), seed_offset=int(m), validate=False)
        done = m
        st = cluster_stats(mc, lat)
        st.update(sweep=int(m), acc=float(acc))
        trace.append(st)
    mc.validate_state()
    if verbose:
        print(f"    {time.time()-t0:.0f} s; final: {trace[-1]['n_clusters']} "
              f"clusters, mean radius {trace[-1]['mean_radius_nm']:.2f} nm, "
              f"largest {trace[-1]['max_radius_nm']:.2f} nm")
    return dict(config=dict(L=L, N=int(lat.N), x=x, T=T, sweeps=sweeps,
                            seed=seed, box_nm=box_nm),
                trace=trace)


def coarsening_exponent(trace, max_frac_in_largest=0.5, min_points=6):
    """
    Fit R ~ t^n over the interval in which coarsening is actually proceeding.

    The fit MUST stop before the solute budget is exhausted. Once one domain
    holds most of the solute there is nothing left to coarsen: the trace goes
    flat, and including that region drags the fitted exponent towards zero. In
    the runs reported here a single domain holds 83 to 99.9 per cent of the
    solute by 3000 sweeps out of 60,000, so more than nine tenths of the trace
    is post-coarsening and must be excluded.

    The window therefore ends when the largest domain first exceeds
    max_frac_in_largest of the solute. It begins after the first few recorded
    points, which cover nucleation rather than coarsening.

    Returns (n, standard error, number of points, sweep at which the window
    closed).
    """
    s = np.array([p["sweep"] for p in trace], float)
    r = np.array([p["mean_radius_nm"] for p in trace], float)
    fl = np.array([p.get("frac_in_largest", 0.0) for p in trace], float)

    stop = np.argmax(fl > max_frac_in_largest)
    if stop == 0:                      # never exceeded, use everything
        stop = len(s)
    ok = (s > 0) & (r > 0)
    idx = np.where(ok)[0]
    idx = idx[idx < stop]
    if len(idx) < min_points:
        return np.nan, np.nan, len(idx), float(s[min(stop, len(s) - 1)])
    idx = idx[2:]                      # drop the earliest, nucleation-dominated
    if len(idx) < min_points:
        return np.nan, np.nan, len(idx), float(s[min(stop, len(s) - 1)])

    ss, rr = s[idx], r[idx]
    A = np.vstack([np.log(ss), np.ones_like(ss)]).T
    coef, *_ = np.linalg.lstsq(A, np.log(rr), rcond=None)
    n = float(coef[0])
    pred = A @ coef
    sigma = np.sqrt(np.sum((np.log(rr) - pred) ** 2) / max(len(ss) - 2, 1))
    se = float(sigma / np.sqrt(np.sum((np.log(ss) - np.log(ss).mean()) ** 2)))
    return n, se, len(ss), float(s[min(stop, len(s) - 1)])


# ======================================================================
def plot_coarsening(results, outfile):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8.5, "pdf.fonttype": 42})

    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.1))
    for key, res in sorted(results.items()):
        tr = res["trace"]
        s = np.array([p["sweep"] for p in tr], float)
        r = np.array([p["mean_radius_nm"] for p in tr])
        nc = np.array([p["n_clusters"] for p in tr])
        mx = np.array([p["max_radius_nm"] for p in tr])
        lab = f"T = {res['config']['T']:.0f} K"
        ok = s > 0
        ax[0].loglog(s[ok], r[ok], "-", lw=1.3, label=lab)
        ax[1].loglog(s[ok], nc[ok], "-", lw=1.3, label=lab)
        ax[2].loglog(s[ok], mx[ok], "-", lw=1.3, label=lab)

    # reference slope
    s_ref = np.array([1e2, 1e4])
    ax[0].loglog(s_ref, 0.25 * (s_ref / s_ref[0]) ** (1 / 3), "k--", lw=0.9,
                 label=r"$t^{1/3}$")
    ax[0].set_xlabel("Monte Carlo sweeps")
    ax[0].set_ylabel("mean domain radius (nm)")
    ax[0].set_title("(a) domain growth", fontsize=8.5, loc="left")
    ax[1].set_xlabel("Monte Carlo sweeps")
    ax[1].set_ylabel("number of domains")
    ax[1].set_title("(b) domain count", fontsize=8.5, loc="left")
    ax[2].set_xlabel("Monte Carlo sweeps")
    ax[2].set_ylabel("largest domain radius (nm)")
    ax[2].set_title("(c) largest domain", fontsize=8.5, loc="left")
    for a in ax:
        a.grid(alpha=0.25, lw=0.5, which="both")
        a.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{outfile}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)


# ======================================================================
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--L", type=int, default=None)
    ap.add_argument("--T", type=float, default=None)
    ap.add_argument("--sweeps", type=int, default=None)
    args = ap.parse_args()

    os.makedirs(_DATA, exist_ok=True)
    os.makedirs(_FIG, exist_ok=True)

    L = args.L if args.L else (16 if args.quick else 30)
    sweeps = args.sweeps if args.sweeps else (4000 if args.quick else 60000)
    temps = [args.T] if args.T else ([1600.0] if args.quick
                                     else [1400.0, 1600.0, 1800.0])

    print("Direct lattice coarsening")
    print(f"  box {L}^3 = {L*A0/10:.1f} nm on a side; interface is one bond "
          f"wide, no inflation\n")

    results, out = {}, os.path.join(_DATA, "lattice_coarsening.json")
    if os.path.exists(out):
        old = json.load(open(out))
        if old.get("meta", {}).get("L") == L and old["meta"]["sweeps"] == sweeps:
            results = {float(k): v for k, v in old["runs"].items()}
            print(f"  resuming; {len(results)} temperature(s) done\n")

    for T in temps:
        if T in results:
            print(f"  T={T:.0f} K already done")
            continue
        res = run_coarsening(L=L, T=T, sweeps=sweeps, seed=11)
        results[T] = res
        json.dump(dict(meta=dict(L=L, sweeps=sweeps, a0=A0),
                       runs={str(k): v for k, v in results.items()}),
                  open(out, "w"), indent=2)

    print(f"\n{'T (K)':>8} {'n':>9} {'+/-':>7} {'pts':>5} {'window ends':>12} "
          f"{'R at end (nm)':>15} {'domains':>9}")
    expos = []
    for T in sorted(results):
        res = results[T]
        n, se, npts, stop = coarsening_exponent(res["trace"])
        if np.isfinite(n):
            expos.append(n)
        # state at the end of the coarsening window, not at the end of the run
        at = min(range(len(res["trace"])),
                 key=lambda i: abs(res["trace"][i]["sweep"] - stop))
        w = res["trace"][at]
        print(f"{T:8.0f} {n:9.3f} {se:7.3f} {npts:5d} {stop:12.0f} "
              f"{w['mean_radius_nm']:15.3f} {w.get('n_domains', -1):9d}")

    if len(expos) > 1:
        print(f"\n  exponent across temperatures: "
              f"{np.mean(expos):.3f} +/- {np.std(expos, ddof=1):.3f}")
        print("  If the exponent is temperature-independent, it is a property")
        print("  of the transport-limited kinetics and can be compared with")
        print("  the phase-field result, which uses a very different interface.")

    plot_coarsening(results, os.path.join(_FIG, "fig_lattice_coarsening"))
    print(f"\nWritten to {out}")
    print(f"Figure written to {_FIG}")
