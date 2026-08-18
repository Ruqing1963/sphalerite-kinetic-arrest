# -*- coding: utf-8 -*-
"""
bulk_free_energy.py  --  Paper 3, Section 2

Measures the bulk free energy density f(c) of the lattice model, which is the
driving force in the Cahn-Hilliard equation.

WHY THIS CANNOT BE ASSUMED

The overwhelmingly common practice in phase-field work is to postulate a
regular-solution or double-well f(c) and fit its coefficients. That is not
available here, and the reason is worth stating because it is specific to this
Hamiltonian.

Under random mixing at x_Cu = x_In = c/2, the pair term of the charge
compensation energy averages to

    <sum_ij dq_i dq_j> = (z/2) [ 2(c/2)^2(-1) + (c/2)^2(+1) + (c/2)^2(+1) ] = 0,

because the heterovalent Cu-In contacts (dq_i dq_j = -1) exactly cancel the
homovalent Cu-Cu and In-In contacts (+1). A mean-field free energy built on
random mixing therefore has NO driving force at all: the entire thermodynamics
of this system lives in the correlations that mean-field theory discards. Any
regular-solution f(c) would describe a different system.

f(c) is therefore measured, by thermodynamic integration over the coupling
strength. Writing H(a) = a * H for a scaling parameter a in [0, 1],

    F(1) - F(0) = int_0^1 da <H>_a,

where F(0) is the ideal solution free energy, known analytically. <H>_a is
sampled at each a. This is exact up to sampling error, and it captures exactly
the correlations that the mean-field treatment misses.
"""

import json
import os
import time

import numpy as np

try:
    from sphalerite_mc import (SphaleriteLattice, HamiltonianParams, KB_EV,
                               _pair_counts)
    from samplers import MixedMoveMC
except ModuleNotFoundError as err:
    here = os.path.dirname(os.path.abspath(__file__))
    print(f"Cannot import '{err.name}'. This script needs sphalerite_mc.py and "
          f"samplers.py in the same folder:\n  {here}")
    raise SystemExit(1)

# ---------------- configuration ----------------
L = 12                       # N = 6912 sites; small enough for a lambda sweep
T = 573.0                    # ore-forming temperature
COMPOSITIONS = [0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20]
N_LAMBDA = 9                 # Gauss-Legendre nodes in the coupling integral
EQUIL, PROD, BLOCKS = 400, 1200, 12
SEED = 31415

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, os.pardir, "data")
os.makedirs(_DATA, exist_ok=True)
OUT = os.path.join(_DATA, "bulk_free_energy.json")


def ideal_entropy_per_site(c):
    """
    Configurational entropy of the ideal solution on the cation sublattice,
    per site, in units of k_B. Species are Zn (1-c), Cu (c/2), In (c/2).
    """
    xs = np.array([1.0 - c, c / 2, c / 2])
    xs = xs[xs > 0]
    return float(-np.sum(xs * np.log(xs)))


def mean_energy(lat, c, scale, seed):
    """
    <H> per site at coupling strength `scale`, by scaling lambda and J together.
    Scaling the whole Hamiltonian is what the integration identity requires;
    scaling only one term would integrate along a different path and give a
    different (wrong) answer.
    """
    par = HamiltonianParams(lam=0.30 * scale)
    par.J = np.array([[0.0, 0.0, 0.0, 0.0],
                      [0.0, +0.05, -0.10, 0.0],
                      [0.0, -0.10, +0.05, 0.0],
                      [0.0, 0.0, 0.0, 0.0]]) * scale
    par.gamma1 = 1.0 * scale
    mc = MixedMoveMC(lat, par, x_cu=c / 2, x_in=c / 2, seed=seed, p_ss=0.5)
    mc.run(T, EQUIL, seed_offset=1, validate=False)
    es = []
    per = max(1, PROD // BLOCKS)
    for b in range(BLOCKS):
        mc.run(T, per, seed_offset=100 + b, validate=False)
        es.append(mc.E / lat.N)
    return float(np.mean(es)), float(np.std(es, ddof=1) / np.sqrt(BLOCKS))


def free_energy(lat, c, seed):
    """
    f(c) per site, by Gauss-Legendre quadrature of <H>_a over a in [0, 1],
    added to the ideal-solution free energy.
    """
    nodes, weights = np.polynomial.legendre.leggauss(N_LAMBDA)
    a = 0.5 * (nodes + 1.0)              # map [-1,1] -> [0,1]
    w = 0.5 * weights
    integral, err2 = 0.0, 0.0
    trace = []
    for k, (ak, wk) in enumerate(zip(a, w)):
        e, se = mean_energy(lat, c, ak, seed + 17 * k)
        # <H>_a is the FULL Hamiltonian evaluated in the ensemble of H(a);
        # the code returns <H(a)> = a<H>, so divide by a to recover <H>_a.
        e_full = e / ak if ak > 0 else 0.0
        se_full = se / ak if ak > 0 else 0.0
        integral += wk * e_full
        err2 += (wk * se_full) ** 2
        trace.append(dict(a=float(ak), E_scaled=e, E_full=e_full, sem=se))
    f_ideal = -KB_EV * T * ideal_entropy_per_site(c)
    return f_ideal + integral, np.sqrt(err2), f_ideal, integral, trace


if __name__ == "__main__":
    results = {}
    if os.path.exists(OUT):
        old = json.load(open(OUT))
        if old.get("config", {}).get("L") == L and old["config"]["T"] == T:
            results = old.get("results", {})
            print(f"resuming; {len(results)} composition(s) done\n")

    lat = SphaleriteLattice(L)
    lat._lut = None
    print(f"Bulk free energy by thermodynamic integration")
    print(f"L={L} (N={lat.N}), T={T:.0f} K, {N_LAMBDA} quadrature nodes\n")
    print(f"{'c':>7} {'f_ideal':>10} {'<H> int':>10} {'f(c)':>11} "
          f"{'+/-':>8} {'min':>6}")

    for c in COMPOSITIONS:
        key = f"{c:.4f}"
        if key in results:
            r = results[key]
            print(f"{c:7.3f} {r['f_ideal']:10.5f} {r['integral']:10.5f} "
                  f"{r['f']:11.5f} {r['f_sem']:8.5f} {'--':>6}")
            continue
        t0 = time.time()
        f, se, f_id, integ, trace = free_energy(lat, c, SEED + int(1e4 * c))
        results[key] = dict(c=c, f=f, f_sem=se, f_ideal=f_id, integral=integ,
                            trace=trace, wall_s=time.time() - t0)
        print(f"{c:7.3f} {f_id:10.5f} {integ:10.5f} {f:11.5f} {se:8.5f} "
              f"{(time.time()-t0)/60:6.1f}")
        json.dump(dict(config=dict(L=L, N=int(lat.N), T=T, n_nodes=N_LAMBDA,
                                   equil=EQUIL, prod=PROD, seed=SEED),
                       results=results), open(OUT, "w"), indent=2)

    # ---- curvature: the spinodal condition ----
    cs = np.array([results[k]["c"] for k in sorted(results, key=float)])
    fs = np.array([results[k]["f"] for k in sorted(results, key=float)])
    if len(cs) >= 4:
        print("\nSecond derivative by cubic-spline fit "
              "(f'' < 0 marks the spinodal region):")
        coef = np.polyfit(cs, fs, min(5, len(cs) - 1))
        d2 = np.polyder(np.poly1d(coef), 2)
        for c in cs:
            flag = "  spinodal" if d2(c) < 0 else ""
            print(f"  c={c:6.3f}   f'' = {d2(c):+10.3f} eV{flag}")
    print(f"\nWritten to {OUT}")
