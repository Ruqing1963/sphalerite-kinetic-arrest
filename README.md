# sphalerite-kinetic-arrest

Why indium nanoinclusions survive in sphalerite when equilibrium demands a single
domain. Coarsening is limited by thermally activated transport, cooling is not,
and the two timescales cross exactly once.

Code, data and manuscript for:

> **He Yu**, **Ruqing Chen**, **Huanzhang Lu**, *Kinetic Arrest, Not
> Equilibrium: Why Indium Nanoinclusions Survive in Sphalerite.*

- He Yu — Geological Laboratory, Hezhou University, Hezhou, Guangxi 542899, China (`yuhe@hzxy.edu.cn`)
- Ruqing Chen — GUT Geoservice Inc., Montreal, Quebec, Canada (`ruqing@hotmail.com`)
- Huanzhang Lu — Université du Québec à Chicoutimi, Quebec, Canada (`hzlu@uqac.uquebec.ca`)

Third paper of the *Computable Mineral Deposit Chemistry* series.

| | subject | code | archive |
|---|---|---|---|
| Paper 1 | the Hamiltonian and its analytic ground state | [sphalerite-lattice-mc](https://github.com/Ruqing1963/sphalerite-lattice-mc) | [10.5281/zenodo.21880502](https://doi.org/10.5281/zenodo.21880502) |
| Paper 2 | equilibrium: indium condenses, at every composition | [sphalerite-equilibrium-mc](https://github.com/Ruqing1963/sphalerite-equilibrium-mc) | [10.5281/zenodo.21911898](https://doi.org/10.5281/zenodo.21911898) |
| Paper 3 | kinetics: why the equilibrium state is not observed | this repository | DOI on release |

---

## The problem

Paper 2 established that the equilibrium state of indium in sphalerite is a
*single* condensed domain, at every composition from 0.10 to 2.00 at.%. This
work adds the interfacial energy in closed form,

&nbsp;&nbsp;&nbsp;&nbsp;γ = λ/3 − J<sub>Cu-In</sub>/3 − J<sub>like</sub>/6 = 0.125 eV per bond,

in which the bulk terms cancel identically, so dividing one domain into two costs
300 to 1800 k<sub>B</sub>T. Bulk thermodynamics opposes multiple domains by an
enormous margin.

Natural sphalerite carries them anyway. The resolution is that an ore body cools:
coarsening is limited by solute transport, which is thermally activated and
therefore slows exponentially, while cooling proceeds at a rate set by the
geology. The two timescales cross once, and below the crossing the microstructure
stops evolving. What is observed is not equilibrium but equilibrium interrupted.

## Results

| quantity | value | source |
|---|---|---|
| Interfacial energy γ | 0.125 eV/bond (closed form), 0.128 (relaxed construction) | Section 3 |
| Splitting cost | 296–1769 k<sub>B</sub>T for n = 128–2048 | Section 3.3 |
| Arrest sensitivity | 10× faster cooling raises T<sub>f</sub> by 21% at 573 K | Eq. (13) |
| Coarsening exponent, phase field | 0.259 | Section 5 |
| Coarsening exponent, lattice | 0.150 ± 0.002 | Section 6 |
| Coarsening exponent, natural end members | ≈ 0.45 | Section 8.3 |
| Simulated domain size, absolute | 0.44 → 1.66 nm in a 13 nm box | Section 6.2 |

**The four exponent determinations span a factor of three and we report that
rather than a consensus.** Section 9.5 suggests they may describe different
regimes — interface-limited coarsening at the atomic scale crossing over to
diffusion-limited as domains grow — but that reading is untested. The inversion
of Section 8 therefore supports *ranking* deposits by cooling rate, not
calibrating one.

## Division of labour between the two methods

This governs how every figure should be read, and is stated in Section 2.5
before any result is presented.

**The phase-field calculation gives the pathway and the topology.** Which
morphologies the system passes through, whether coarsening completes or arrests,
and how the outcome depends on the ratio of cooling rate to mobility. It
*cannot* give an absolute length: the Cahn–Hilliard interface must be widened
numerically for any solver to resolve it, so domain sizes there are set by the
grid. Inflating κ rescales lengths and times together, so the dimensionless
conclusions survive.

**The lattice calculation gives the absolute scale.** The interface is one bond
wide and nothing is inflated, so nanometres are nanometres. It reaches far less
time and far fewer domains.

## Requirements

```bash
pip install numpy scipy matplotlib numba
```

Numba is optional but gives roughly a fiftyfold speed-up in the lattice code.
The phase-field solver uses NumPy's FFT only.

The DFT input generator additionally needs ASE if the Quantum ESPRESSO route is
used:

```bash
pip install ase
```

No first-principles calculations were run for these papers; see *Parameters
still to be determined* below.

## Reproducing the results

All scripts resolve paths relative to themselves and can be run from any working
directory. Long runs checkpoint and resume; deleting the output JSON forces a
fresh start.

### Paper 3, Section 2 — bulk free energy

```bash
cd src
python bulk_free_energy.py          # ~1 h
```

Measures *f(c)* by thermodynamic integration over the coupling strength at eight
compositions. This cannot be replaced by a regular-solution form: under random
mixing the heterovalent and homovalent contributions cancel *identically*, so a
mean-field free energy has no driving force at all in this system.

### Paper 3, Section 5 — phase-field cooling

```bash
python cooling_2d.py                # three cooling rates, ~10 min at 256²
python cooling_2d.py --quick        # single fast check
python cooling_2d.py --steps 32000  # one schedule; larger = slower cooling
```

Produces the morphological sequence (homogeneous → spinodal decomposition →
coarsening → arrest) and the cooling-rate comparison. Cooling is specified in
*steps* rather than physical time because the timestep already tracks the
mobility, so a fixed step count per kelvin is a fixed amount of diffusive
progress per kelvin — the dimensionless ratio that controls the outcome.

### Paper 3, Section 6 — direct lattice coarsening

```bash
python lattice_coarsening.py        # three temperatures, ~1 h
python lattice_coarsening.py --quick
```

Gives absolute domain sizes with the interface at its physical width. Two things
in this script were forced by data and matter for anyone adapting it: monomers
must be excluded from the mean radius, and the fit must stop before the solute
budget is exhausted into one domain. Both are documented in the source.

### Paper 2 — equilibrium condensation

In the [Paper 2 repository](https://github.com/Ruqing1963/sphalerite-equilibrium-mc):

```bash
python sphalerite_mc.py             # verification suite, ~5 s
python composition_scan.py          # ~12 h, checkpointed
python condensation_bound.py        # the bracketing test, ~30 min
python coarsening_test.py           # single vs multiple domains, ~11 h
```

### Parameters still to be determined

The chemical term *J*<sub>αβ</sub> is not fixed by first principles in any of
the three papers. The campaign that would fix it is specified and its tooling
validated, but not run: a 216-atom supercell needs 18.7 GB by the code's own
estimate, against ~12 GB available on the workstation used throughout.

```bash
python dft_inputs.py                        # write 12 POSCARs, INCAR, KPOINTS
python qe_inputs.py                         # convert to Quantum ESPRESSO
# ... run VASP or pw.x on each ...
python collect_energies.py                  # gather into energies.csv
python dft_inputs.py --fit energies.csv     # least squares for λ and J
```

The generated structures reproduce the analytic values of Σ(ΔQ)² predicted in
Paper 1 — 6λ for a first-neighbour Cu–In pair, 8λ for the four-defect ground
state, 16λ and 24λ for the vacancy configurations — which checks the generation
before any energy is computed. The fitting harness was validated by recovering
known parameters from synthetic energies carrying 5 meV of noise.

Anyone with sufficient memory can run this and drop the parameters into the
model without repeating any of the sampling reported here.

## Notes for anyone building on this

Four things cost us time and are recorded so they need not cost anyone else's.

**A polynomial fitted to a narrow composition window cannot serve as a global
free energy.** A quartic through our eight measured points has a *negative*
leading coefficient, so f → −∞ and the field diverges within a few hundred
steps. We embed the measurements in a bounded double well instead.

**In a semi-implicit scheme, the timestep should be set from the target
amplification, not from an explicit CFL condition.** Setting it the latter way
gives a per-step growth of 3×10⁻⁴, so the unstable mode needs ~7000 steps to
grow tenfold while the run has a few thousand. The result looks like coarsening
and is diffusive smoothing.

**f, f′ and f″ must be checked against finite differences.** Dropping a term
from an analytic derivative made f′ inconsistent with f″ by a factor of 150 at
the working composition. The solver ran, produced plausible pictures, and
described a free energy that was not the one reported. This is now asserted at
construction.

**The coarsening exponent does not identify the mechanism.** Two versions of the
phase-field calculation in which decomposition never occurred returned exponents
within 7% of the correct value, because diffusive pattern coarsening obeys a
similar law to domain coarsening. What exposed the error was the composition
contrast, which must grow by orders of magnitude and did not. We report it
alongside every exponent.

## Repository layout

```
src/
  sphalerite_mc.py        core lattice module, unchanged from Paper 1
  samplers.py             replica exchange and mixed move set
  calibrate_ladder.py     adaptive ladder construction
  bulk_free_energy.py     thermodynamic integration for f(c)
  cahn_hilliard.py        general Cahn-Hilliard solver, 2D and 3D
  cooling_2d.py           the cooling simulations of Section 5
  lattice_coarsening.py   the direct coarsening of Section 6
data/     numerical results underlying every figure and table
figures/  figures in vector and raster form
paper/    manuscript in plain and Elsevier formats, plus a Chinese
          translation; LaTeX source and compiled PDF for each
```

## Citation

> Yu, H., Chen, R., Lu, H. (2026) *Kinetic Arrest, Not Equilibrium: Why Indium
> Nanoinclusions Survive in Sphalerite.*
> https://github.com/Ruqing1963/sphalerite-kinetic-arrest

## Licence

- **Code** (`src/`): MIT — see `LICENSE`.
- **Data and figures** (`data/`, `figures/`): CC BY 4.0 — see `LICENSE-DATA`.
- **Manuscript** (`paper/`): © the authors.
