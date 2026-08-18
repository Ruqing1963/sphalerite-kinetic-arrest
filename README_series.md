# Computable Sphalerite Chemistry

Thermodynamics and kinetics of Cu–In coupled substitution in sphalerite, derived
rather than documented.

Indium is essentially never mined for its own sake. There are no indium ores in
the conventional sense: the overwhelming majority of primary production is
recovered as a by-product of zinc refining, so the global indium endowment is
controlled not by a distinct ore-forming process but by the trace-element
chemistry of sphalerite itself. That chemistry has been documented extensively
and derived hardly at all.

This project derives it, from a configurational Hamiltonian whose central
coefficient is fixed by a measured dielectric constant rather than by fitting,
and follows the consequences to the point where they meet what a microscope
actually shows.

**Authors** — He Yu (Hezhou University, `yuhe@hzxy.edu.cn`), Ruqing Chen (GUT
Geoservice Inc., `ruqing@hotmail.com`), Huanzhang Lu (Université du Québec à
Chicoutimi, `hzlu@uqac.uquebec.ca`)

---

## The three papers

Each answers a question the previous one raised and could not settle.

### Paper 1 — Topology, and an electrostatic ground state

*A statistical thermodynamic lattice model for Cu–In coupled substitution in
sphalerite.*
[code](https://github.com/Ruqing1963/sphalerite-lattice-mc) ·
[10.5281/zenodo.21880502](https://doi.org/10.5281/zenodo.21880502)

Indium in sphalerite is almost invariably accompanied by copper in a near-1:1
molar ratio, across deposit types, and this is universally attributed to the
coupled substitution 2 Zn²⁺ → Cu⁺ + In³⁺. The relation had never been derived.

The device that makes derivation possible is an exact topological property of
zinc blende: **two nearest-neighbour cations share exactly one bridging anion,
and higher neighbours share none.** This reduces the tetrahedral electroneutrality
functional, *without approximation*, to a strictly nearest-neighbour pair
interaction whose coefficient follows from the dielectric constant of ZnS:

&nbsp;&nbsp;&nbsp;&nbsp;λ = e²/(8πε₀ε_r d_NN) ∈ [0.227, 0.368] eV

The Cu–In pair binding energy is then −2λ ≈ −0.6 eV, an order of magnitude above
k_BT at hydrothermal temperatures, with no adjustable parameter anywhere. The
ground state at Cu:In = 1:1 is shown analytically to be exactly the roquesite
cation ordering of CuInS₂ — recovered from a purely electrostatic functional,
with no input describing that structure.

Elastic misfit, retained at full strength in a control calculation, produces no
detectable ordering: the largest strain-mediated pair interaction is some forty
times smaller than k_BT.

### Paper 2 — Equilibrium: indium condenses

*The equilibrium fate of indium in sphalerite: large-scale replica-exchange
simulation and the absence of a dilute solid solution.*
[code](https://github.com/Ruqing1963/sphalerite-equilibrium-mc) ·
[10.5281/zenodo.21911898](https://doi.org/10.5281/zenodo.21911898)

Paper 1 could not settle whether indium sits in solid solution or aggregates: at
4000 sites and 4 at.% solute the solute budget is exhausted by a single domain,
so complete condensation was consistent with, but not diagnostic of, a
thermodynamic preference.

In a box eight times larger, across compositions from 0.10 to 2.00 at.%, the
solute condenses completely at every point. **The monomer fraction is zero
throughout, against a random-solution expectation of 98% at the most dilute
point** — a discrepancy of four orders of magnitude in the largest-cluster
fraction.

The result does not rest on convergence at ore-forming temperatures, which the
sampler cannot guarantee. Each state point is run from a dispersed *and* from a
pre-condensed configuration; these approach equilibrium from opposite sides and
therefore bracket it, whether or not either has converged. At 1800 K, where the
acceptance ratio is near 10%, dispersed solutes spontaneously aggregate while
condensed solutes do not disperse.

Two subsidiary results bear on how trace-element data are read. The pair
enrichment factor rises sixteenfold across the composition scan while the
quantity it describes — the copper content of the first shell around indium —
*falls* slightly; the apparent variation lies almost entirely in a normalisation
whose ceiling is 1/x_Cu. And ideal dilute-mixing entropy, assumed in
conventional partition-coefficient models, is not available if the equilibrium
state is condensed.

### Paper 3 — Kinetics: why equilibrium is not observed

*Kinetic arrest, not equilibrium: why indium nanoinclusions survive in
sphalerite.*
[code](https://github.com/Ruqing1963/sphalerite-kinetic-arrest) · [10.5281/zenodo.21995641](https://doi.org/10.5281/zenodo.21995641)

Paper 2 says the equilibrium state is a **single** domain, and this paper adds
the interfacial energy in closed form,

&nbsp;&nbsp;&nbsp;&nbsp;γ = λ/3 − J_Cu-In/3 − J_like/6 = 0.125 eV per bond,

in which the bulk terms cancel identically, so dividing one domain into two
costs 300 to 1800 k_BT. Natural sphalerite carries multiple domains anyway.

The resolution is that an ore body cools. Coarsening is limited by solute
transport, which is thermally activated and slows *exponentially*; the cooling
timescale falls only as T². **Two such functions cross, and cross once.** Below
the crossing the microstructure stops evolving, and what survives is not
equilibrium but equilibrium interrupted.

The arrest condition is strongly buffered: every quantity entering it, including
the diffusion prefactor this work cannot determine, appears inside a logarithm,
so an order-of-magnitude error shifts the arrest temperature by less than a
quarter. The frozen domain size follows a power law in cooling rate, which
inverts: **inclusion size records thermal history.**

One consequence inverts a standard reading. A smooth LA-ICP-MS depth profile has
often been taken as evidence for solid solution; on this account it is evidence
for *fast cooling*, the indium being aggregated but below the resolution of the
ablation volume.

---

## What is not settled

Stated here rather than buried, because these bound what the work supports.

**The chemical term J is not fixed by first principles.** The campaign that
would fix it is specified and its tooling validated, but not run: a 216-atom
supercell needs 18.7 GB by the code's own estimate, against ~12 GB available on
the workstation used throughout. Everything needed to run it is in the
repositories. Conclusions are classified in each paper by whether they survive
the removal of J — the central ones do, since electrostatics alone gives a
sixteenfold enrichment and removing J *strengthens* the condensation claim.

**The coarsening exponent is not settled.** Four determinations span a factor of
three: 1/3 from Lifshitz–Slyozov theory, 0.259 from phase field, 0.150 from
direct lattice coarsening, ≈0.45 implied by natural end members. They may
describe different regimes — interface-limited coarsening at the atomic scale
crossing over to diffusion-limited as domains grow — but that is untested. The
geospeedometer of Paper 3 therefore supports *ranking* deposits by cooling rate,
not calibrating one.

**Natural tenors lie an order of magnitude below the compositions simulated.**
An entropic estimate places the composition boundary between roughly 3 ppm and
10⁻⁵ ppm — below crustal abundance — implying that no natural sphalerite
crystallises above the solvus, but locating it requires boxes of order 5×10⁵
sites.

**Monte Carlo sweeps are not physical time**, and iron is absent from the model
throughout.

---

## Requirements

```bash
pip install numpy scipy matplotlib numba
```

Numba is optional but gives roughly a fiftyfold speed-up in the lattice code;
without it the Paper 2 composition scan takes days rather than hours. The
phase-field solver uses NumPy's FFT only and does not need Numba.

For the DFT input generator, if the Quantum ESPRESSO route is used:

```bash
pip install ase
```

All calculations reported in the three papers were performed on a single
workstation (8 cores, 16 GB). We regard that as a feature: every result can be
reproduced on ordinary hardware.

---

## Reproducing the results

Scripts resolve paths relative to themselves and run from any working directory.
Long runs checkpoint and resume; deleting the output JSON forces a fresh start.

### Start here — the verification suite

```bash
cd src
python sphalerite_mc.py          # ~5 s
```

Checks lattice topology, the sharing lemma on random site samples, energy
bookkeeping against full recomputation, composition conservation, and every
analytic limit to machine precision — including the binding energy −2λ = −0.6 eV
that the whole series rests on. Everything else depends on this passing.

### Paper 2 — equilibrium condensation

```bash
python composition_scan.py       # ~12 h, checkpointed per composition
python condensation_bound.py     # the bracketing test, ~30 min
python coarsening_test.py        # single vs multiple domains, ~11 h
python scale_profile.py          # verification and profiling to N = 5e5
```

`composition_scan.py` produces Table 1 of Paper 2 and `condensation_bound.py`
the bracketing argument on which its central claim rests.

### Paper 3, Section 6 — lattice coarsening, absolute scale

```bash
python lattice_coarsening.py            # three temperatures, ~1 h
python lattice_coarsening.py --quick    # fast check
python lattice_coarsening.py --L 30 --T 1600 --sweeps 60000
```

This is where absolute sizes come from. The interface is one bond wide and
nothing is inflated, so nanometres are nanometres: at L = 24 (a 13.0 nm cube)
and 1600 K the domain count falls from 173 to 4 while the mean domain radius
grows from 0.44 to 1.66 nm.

Two design points in this script were forced by data and matter to anyone
adapting it. **Monomers must be excluded from the mean radius**: at 1800 K some
90% of clusters are single solutes evaporated from domains, and including them
makes the "mean radius" a measure of free-solute concentration that *falls*
while domains grow. **The fit must stop before the solute budget is exhausted**:
at 1 at.% the solutes condense into one domain within 3000 of 60,000 sweeps, so
most of the run measures a structure that has stopped evolving. Both are
documented in the source.

### Paper 3, Section 5 — phase-field cooling, the four stages

```bash
python cooling_2d.py                 # three cooling rates, ~10 min at 256²
python cooling_2d.py --quick         # single fast run
python cooling_2d.py --steps 32000   # one schedule; larger = slower cooling
```

Produces the morphological sequence — homogeneous → spinodal decomposition →
coarsening → arrest — and the cooling-rate comparison. The composition contrast
rises from 2×10⁻³ to 0.308 as decomposition completes, then the domain count
falls while the surviving domains grow, and finally the curves flatten as the
mobility collapses. That flattening is arrest.

Cooling is specified in **steps** rather than physical time. The timestep
already tracks the mobility, so a fixed step count per kelvin is a fixed amount
of *diffusive progress* per kelvin — the dimensionless ratio that controls the
outcome. A larger step count is slower cooling.

**Read the lengths as grid units, not nanometres.** The Cahn–Hilliard interface
must be widened numerically for any solver to resolve it, so domain sizes here
are set by the grid rather than predicted. Inflating κ rescales lengths and
times together, so the dimensionless conclusions survive; the absolute scale
comes from `lattice_coarsening.py` instead.

### Paper 3, Section 2 — the bulk free energy

```bash
python bulk_free_energy.py       # ~1 h
```

Measures *f(c)* by thermodynamic integration over the coupling strength. This
cannot be replaced by a regular-solution form, and the obstruction is specific
rather than a matter of accuracy: under random mixing the heterovalent and
homovalent contributions **cancel identically**, so a mean-field free energy has
no driving force at all in this system. It would describe a different one.

### First-principles parameterisation — the bridge to J and λ

Not run for these papers, for the memory reason above, but complete and
validated.

```bash
python dft_inputs.py                        # 12 POSCARs + INCAR + KPOINTS
python qe_inputs.py                         # convert to Quantum ESPRESSO
# ... run VASP or pw.x on each structure ...
python collect_energies.py                  # gather into energies.csv
python dft_inputs.py --fit energies.csv     # least squares for λ, J_Cu-In, J_like
```

The design exploits a feature of the problem. The quantities wanted are
*differences* between configurations at the same composition, so the structures
are charge neutral and composition matched within each set: no chemical
potentials enter, no charged-defect corrections are needed, and the band-gap
error that afflicts defect calculations in semiconductors does not propagate.

Three sets of 3×3×3 cells (216 atoms) suffice. **Set A** places one Cu and one
In at increasing separation, determining λ and J_Cu-In jointly while testing the
strict nearest-neighbour range the pair reduction requires — the binding energy
must vanish beyond the first shell. **Set B** uses four topologies of 2Cu + 2In
to separate J_like from J_Cu-In. **Set C** pairs a Zn vacancy with two In,
giving the Cu-versus-vacancy competition at fixed In content, far more robust
than an absolute vacancy formation energy.

Parameters are recovered by least squares over within-set energy differences,
using descriptors computed exactly by the lattice code, so the residual is
itself a test of whether the four-term Hamiltonian is adequate. Two checks are
built in: the generated structures reproduce the analytic Σ(ΔQ)² values of
Paper 1 — 6λ for a first-neighbour Cu–In pair, 8λ for the four-defect ground
state, 16λ and 24λ for the vacancy configurations — before any energy is
computed; and the fitting harness was validated by recovering known parameters
from synthetic energies carrying 5 meV of noise.

`collect_energies.py` **refuses to write the CSV while any job is unconverged**,
since an unconverged energy in a least-squares fit biases every parameter rather
than merely widening the error.

---

## Notes for anyone building on this

Things that cost us time, recorded so they need not cost anyone else's.

**Two agreeing trajectories are not an equilibration test.** In a system with an
acceptance ratio of 10⁻⁵, two chains launched into the same basin agree with
each other whether or not that basin is the equilibrium one. Paper 1 drew a
stronger conclusion from such agreement than it could support; the correction is
Appendix A of Paper 2, and it is why the central claim there rests on bracketing
rather than on convergence.

**Report replica round trips, not swap acceptance.** Reducing the interval
between exchange attempts from five sweeps to one takes the round-trip count
from 0 to 7 at identical cost, while the median swap acceptance moves by less
than 0.01. A study reporting acceptance alone would conclude its ladder was
performing well while no replica traversed it.

**An enrichment factor normalised by bulk concentration will appear to
strengthen on dilution** whether or not anything physical changes, since its
ceiling is set by that concentration. Report the un-normalised quantity
alongside. This applies to LA-ICP-MS datasets spanning orders of magnitude in
tenor as much as to simulation.

**A polynomial fitted to a narrow composition window cannot serve as a global
free energy.** A quartic through our eight measured points has a *negative*
leading coefficient, so f → −∞ and the field diverges within a few hundred
steps.

**In a semi-implicit scheme, set the timestep from the target amplification, not
from an explicit CFL condition.** The latter gives a per-step growth of 3×10⁻⁴,
so the unstable mode needs ~7000 steps to grow tenfold while the run has a few
thousand. The result looks like coarsening and is diffusive smoothing.

**Check f, f′ and f″ against finite differences.** Dropping a term from an
analytic derivative made f′ inconsistent with f″ by a factor of 150 at the
working composition. The solver ran, produced plausible pictures, and described
a free energy that was not the one reported.

**The coarsening exponent does not identify the mechanism.** Two versions of the
phase-field calculation in which decomposition never occurred returned exponents
within 7% of the correct value, because diffusive pattern coarsening obeys a
similar law to domain coarsening. What exposed the error was the composition
contrast, which must grow by orders of magnitude and did not.

---

## Repository layout

The series is split across three repositories, one per paper, each
self-contained. The core lattice module `sphalerite_mc.py` is identical in all
three.

```
src/       simulation and analysis code
data/      numerical results underlying every figure and table (JSON + CSV)
figures/   figures in vector and raster form, with the scripts that make them
paper/     manuscript in plain and Elsevier formats, plus a Chinese
           translation; LaTeX source and compiled PDF for each
dft/       supercell structures and input files for the parameterisation
```

## Citation

Cite the paper whose result you use:

> Yu, H., Chen, R., Lu, H. (2026) *A statistical thermodynamic lattice model for
> Cu–In coupled substitution in sphalerite*, v2.0. Zenodo.
> https://doi.org/10.5281/zenodo.21880502

> Yu, H., Chen, R., Lu, H. (2026) *The equilibrium fate of indium in sphalerite:
> large-scale replica-exchange simulation and the absence of a dilute solid
> solution.* Zenodo. https://doi.org/10.5281/zenodo.21911898

> Yu, H., Chen, R., Lu, H. (2026) *Kinetic arrest, not equilibrium: why indium
> nanoinclusions survive in sphalerite.* Zenodo.
> https://doi.org/10.5281/zenodo.21995641

Each repository carries a `CITATION.cff` with machine-readable metadata.

## Licence

- **Code**: MIT
- **Data and figures**: CC BY 4.0
- **Manuscripts**: © the authors

## Funding

Hezhou Municipal Scientific Research and Development Program, Project Nos.
2024143, 2024141 and 2024104.
