# -*- coding: utf-8 -*-
"""
cahn_hilliard.py  --  Paper 3, Section 5

Semi-implicit spectral solver for the Cahn-Hilliard equation with a
temperature-dependent mobility, used to follow the morphological evolution of a
cooling sphalerite grain.

    dc/dt = div [ M(T) grad ( df/dc - kappa laplacian c ) ]

WHY SPECTRAL, AND WHY SEMI-IMPLICIT

Cahn-Hilliard is fourth order in space, so an explicit scheme requires
dt < C dx^4 / (M kappa), which at the grid resolutions needed here means
timesteps too small to reach coarsening times. The standard remedy is to treat
the linear fourth-order term implicitly in Fourier space, where it is diagonal,
and the nonlinear bulk term explicitly:

    c_hat^{n+1} = ( c_hat^n - dt M k^2 F[df/dc]^n ) / ( 1 + dt M kappa k^4 )

This is unconditionally stable in the linear term and allows timesteps larger by
orders of magnitude. Periodic boundaries are implicit in the Fourier transform,
which suits a bulk grain interior.

WHAT THIS SOLVER CAN AND CANNOT TELL US

It supplies the PATHWAY and the TOPOLOGY: the sequence of morphologies, whether
coarsening completes or arrests, and how that depends on the ratio of cooling
rate to mobility. Those are dimensionless statements.

It does NOT supply the absolute length scale. The characteristic length is
sqrt(kappa/|f''|), and with kappa fixed from an atomically sharp interface it is
of order one lattice parameter, which no grid can resolve: numerical stability
needs an interface spanning four to six cells, so kappa must be inflated and the
interface artificially widened. Domain sizes here are therefore set by the grid,
not predicted. Absolute sizes come from direct lattice-model coarsening, where
the interface is what it physically is.

Inflating kappa rescales lengths and times together, so dimensionless ratios --
which is what the conclusions rest on -- are unaffected.
"""

from __future__ import annotations

import numpy as np


KB_EV = 8.617333262e-5


# ======================================================================
# Free energy
# ======================================================================
class LatticeFreeEnergy:
    """
    Bulk free energy f(c) interpolated from the thermodynamic-integration
    measurements of bulk_free_energy.py, with its first two derivatives.

    A polynomial is fitted to the measured points, but only to LOW degree and
    only for use as a driving force. Section 2.3 of the manuscript shows that
    the second derivative of such a fit oscillates with the degree chosen and
    must not be read as a measurement of the spinodal; the instability there is
    diagnosed by the common-tangent construction instead. The derivative used
    here is the smooth interpolant, adequate for driving the dynamics but not a
    quantitative claim about f''.
    """

    def __init__(self, c_data, f_data, degree=4):
        self.c = np.asarray(c_data, float)
        self.f = np.asarray(f_data, float)
        self.degree = int(degree)
        self.coef = np.polyfit(self.c, self.f, self.degree)
        self._p = np.poly1d(self.coef)
        self._dp = np.polyder(self._p, 1)
        self._d2p = np.polyder(self._p, 2)

    def f_of_c(self, c):
        return self._p(c)

    def dfdc(self, c):
        return self._dp(c)

    def d2fdc2(self, c):
        return self._d2p(c)

    def spinodal_range(self, lo=None, hi=None, n=400):
        """Composition interval over which the interpolant has f'' < 0."""
        lo = self.c.min() if lo is None else lo
        hi = self.c.max() if hi is None else hi
        cc = np.linspace(lo, hi, n)
        neg = cc[self._d2p(cc) < 0]
        return (float(neg.min()), float(neg.max())) if neg.size else None


class DoubleWell:
    """
    f(c) = W c^2 (1-c)^2, the standard symmetric double well.

    Provided for testing the solver against a case with an analytic interface
    profile and a known interfacial energy, gamma = sqrt(2 kappa W)/6, so that
    numerical error can be separated from the uncertainty in the measured f(c).
    """

    def __init__(self, W=1.0):
        self.W = float(W)

    def f_of_c(self, c):
        return self.W * c ** 2 * (1 - c) ** 2

    def dfdc(self, c):
        return 2 * self.W * c * (1 - c) * (1 - 2 * c)

    def d2fdc2(self, c):
        return 2 * self.W * (1 - 6 * c + 6 * c ** 2)


# ======================================================================
# Mobility
# ======================================================================
class ArrheniusMobility:
    """
    M(T) = M0 exp(-Ea / kB T).

    Ea is bounded below by the model itself: removing a solute from a domain
    breaks heterovalent bonds costing 2*lambda ~ 0.6 eV each, so Ea >= 2*lambda
    for a surface solute. That is a bound, not a determination -- the true
    barrier also contains the migration barrier of whatever diffusion mechanism
    operates, which the configurational model does not represent.

    M0 is left at unity unless measured diffusion data are supplied. With
    M0 = 1 the time axis is in units of the diffusive time at infinite
    temperature, and every result is a RELATIVE calibration: a ratio of domain
    sizes maps onto a ratio of cooling rates without either being absolute.
    """

    def __init__(self, Ea_eV=0.60, M0=1.0):
        self.Ea = float(Ea_eV)
        self.M0 = float(M0)

    def __call__(self, T):
        return self.M0 * np.exp(-self.Ea / (KB_EV * T))


# ======================================================================
# Solver
# ======================================================================
class CahnHilliard:
    """
    Semi-implicit spectral integrator on a periodic square or cubic grid.

    Parameters
    ----------
    shape   : (Nx, Ny) or (Nx, Ny, Nz)
    dx      : grid spacing, in whatever length unit kappa is expressed in
    kappa   : gradient energy coefficient
    free    : object with dfdc(c); see LatticeFreeEnergy or DoubleWell
    mobility: callable M(T), or a float for isothermal runs
    """

    def __init__(self, shape, dx, kappa, free, mobility):
        self.shape = tuple(shape)
        self.dim = len(self.shape)
        if self.dim not in (2, 3):
            raise ValueError("shape must be 2- or 3-dimensional")
        self.dx = float(dx)
        self.kappa = float(kappa)
        self.free = free
        self.mobility = mobility

        # wavevectors; k2 = |k|^2 and k4 = |k|^4 on the real FFT grid
        ks = [2 * np.pi * np.fft.fftfreq(n, d=self.dx) for n in self.shape[:-1]]
        ks.append(2 * np.pi * np.fft.rfftfreq(self.shape[-1], d=self.dx))
        grids = np.meshgrid(*ks, indexing="ij")
        self.k2 = sum(g ** 2 for g in grids)
        self.k4 = self.k2 ** 2

    # ------------------------------------------------------------------
    def initialise(self, c0, noise=0.01, seed=0):
        """Uniform composition c0 plus small random fluctuations."""
        rng = np.random.default_rng(seed)
        self.c = c0 + noise * (rng.random(self.shape) - 0.5)
        return self.c

    def step(self, dt, T=None):
        """
        Advance one timestep.

        The linear fourth-order term is treated implicitly in Fourier space,
        where it is diagonal; the nonlinear bulk term explicitly. Writing
        mu_bulk = df/dc,

            c_hat <- ( c_hat - dt M k^2 F[mu_bulk] ) / ( 1 + dt M kappa k^4 ).
        """
        M = self.mobility(T) if callable(self.mobility) else float(self.mobility)
        c_hat = np.fft.rfftn(self.c)
        mu_hat = np.fft.rfftn(self.free.dfdc(self.c))
        self.c = np.fft.irfftn(
            (c_hat - dt * M * self.k2 * mu_hat) / (1.0 + dt * M * self.kappa * self.k4),
            s=self.shape, axes=tuple(range(self.dim)))
        return self.c

    # ------------------------------------------------------------------
    def free_energy(self):
        """Total free energy, bulk plus gradient, per unit volume."""
        grad2 = np.zeros_like(self.c)
        for ax in range(self.dim):
            g = (np.roll(self.c, -1, ax) - np.roll(self.c, 1, ax)) / (2 * self.dx)
            grad2 += g ** 2
        return float(np.mean(self.free.f_of_c(self.c) + 0.5 * self.kappa * grad2))

    def domain_scale(self):
        """
        Characteristic domain size from the first moment of the structure
        factor, L = 2 pi / <k>, which is the standard measure in coarsening
        studies and is far less noisy than counting domains.
        """
        c = self.c - self.c.mean()
        S = np.abs(np.fft.rfftn(c)) ** 2
        k = np.sqrt(self.k2)
        num = float(np.sum(k * S))
        den = float(np.sum(S))
        if num <= 0 or den <= 0:
            return np.nan
        return float(2 * np.pi * den / num)

    def solute_fraction_in_domains(self, threshold=None):
        """Fraction of the box above the midpoint composition."""
        thr = self.c.mean() if threshold is None else threshold
        return float(np.mean(self.c > thr))


# ======================================================================
# Cooling protocol
# ======================================================================
def cool(sim, T_start, T_end, rate, dt, record_every=100, callback=None):
    """
    Integrate while cooling linearly at `rate` kelvin per unit time.

    Returns a trace of temperature, domain scale and free energy. The cooling
    rate enters only through the dimensionless ratio of the cooling time to the
    diffusive time, so the trace is a relative calibration unless the mobility
    carries absolute units.
    """
    T = float(T_start)
    t = 0.0
    trace = {k: [] for k in ("t", "T", "L", "F", "frac")}
    step = 0
    while T > T_end:
        sim.step(dt, T)
        t += dt
        T -= rate * dt
        step += 1
        if step % record_every == 0:
            trace["t"].append(t)
            trace["T"].append(T)
            trace["L"].append(sim.domain_scale())
            trace["F"].append(sim.free_energy())
            trace["frac"].append(sim.solute_fraction_in_domains())
            if callback is not None:
                callback(step, t, T, sim)
    return {k: np.asarray(v) for k, v in trace.items()}


# ======================================================================
if __name__ == "__main__":
    import json, os, time

    print("Cahn-Hilliard solver self-test\n")

    # --- 1. analytic check against the double well ---
    # For f = W c^2 (1-c)^2 the equilibrium planar interface has
    #   c(x) = 1/2 [1 + tanh(x / (2 xi))],  xi = sqrt(kappa / 2W),
    # and gamma = sqrt(2 kappa W) / 6. Relaxing a step profile should
    # reproduce the tanh width.
    W, kappa = 1.0, 0.5
    xi_exact = np.sqrt(kappa / (2 * W))
    dw = DoubleWell(W)
    sim = CahnHilliard((256, 4), dx=0.25, kappa=kappa, free=dw, mobility=1.0)
    x = np.arange(256) * 0.25
    sim.c = np.where((x[:, None] > 24) & (x[:, None] < 40),
                     np.ones((256, 4)), np.zeros((256, 4)))
    for _ in range(4000):
        sim.step(dt=0.02)
    prof = sim.c[:, 0]
    # width from the maximum gradient of the relaxed profile
    g = np.abs(np.gradient(prof, 0.25))
    xi_num = 1.0 / (4 * g.max()) if g.max() > 0 else np.nan
    print(f"  interface width: numerical {xi_num:.4f}, analytic {xi_exact:.4f}, "
          f"ratio {xi_num/xi_exact:.3f}")

    # --- 2. free energy must decrease monotonically ---
    sim2 = CahnHilliard((128, 128), dx=1.0, kappa=1.0, free=DoubleWell(1.0),
                        mobility=1.0)
    sim2.initialise(0.5, noise=0.02, seed=1)
    Fs, Ls = [], []
    for n in range(2000):
        sim2.step(dt=0.05)
        if n % 200 == 0:
            Fs.append(sim2.free_energy()); Ls.append(sim2.domain_scale())
    Fs = np.asarray(Fs)
    mono = bool(np.all(np.diff(Fs) <= 1e-12))
    print(f"  free energy monotone decreasing: {mono}")
    print(f"    F: {Fs[0]:.5f} -> {Fs[-1]:.5f}")
    print(f"    L: {Ls[0]:.2f} -> {Ls[-1]:.2f}  (should grow)")

    # --- 3. coarsening exponent, L ~ t^(1/3) ---
    sim3 = CahnHilliard((128, 128), dx=1.0, kappa=1.0, free=DoubleWell(1.0),
                        mobility=1.0)
    sim3.initialise(0.5, noise=0.02, seed=2)
    ts, ls = [], []
    t = 0.0
    for n in range(20000):
        sim3.step(dt=0.05); t += 0.05
        if n > 2000 and n % 500 == 0:
            ts.append(t); ls.append(sim3.domain_scale())
    ts, ls = np.asarray(ts), np.asarray(ls)
    good = np.isfinite(ls) & (ls > 0)
    expo = np.polyfit(np.log(ts[good]), np.log(ls[good]), 1)[0]
    print(f"  coarsening exponent: {expo:.3f}  (Lifshitz-Slyozov: 0.333)")

    # --- 4. lattice free energy, if it has been measured ---
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        os.pardir, "data", "bulk_free_energy.json")
    if os.path.exists(path):
        d = json.load(open(path))["results"]
        cs = np.array([d[k]["c"] for k in sorted(d, key=float)])
        fs = np.array([d[k]["f"] for k in sorted(d, key=float)])
        lf = LatticeFreeEnergy(cs, fs, degree=4)
        sp = lf.spinodal_range()
        print(f"\n  lattice f(c) loaded, {len(cs)} points")
        print(f"    interpolant f'' < 0 over c = {sp}"
              if sp else "    interpolant has no unstable region")
        print("    (the manuscript diagnoses the instability by the")
        print("     common-tangent test, not by this interpolant)")
    print("\nself-test complete")
