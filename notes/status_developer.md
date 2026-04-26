# HBFree — Developer Status

**Date: 2026-04-26 (transient document — update as work progresses)**

---

## What HBFree does today

### Core capability

HBFree is a Harmonic Balance (HB) frequency-domain simulator for nonlinear circuits.
Given a netlist in HBL format (or SPICE via the `s2h` converter), it computes the
steady-state phasor spectrum at all circuit nodes under single-tone or multi-tone
sinusoidal excitation.

**Solver**: Newton-Raphson with LU factorization (`lineq1/lucan/luslv`).
Iteration limit: 50 (configurable via `.hb_options`). Convergence criteria: max step
size (`epsdu`), solution norm (`epsiw`), and residual (`epssol`).

**Output**: `.raw` file with complex phasors. Convention: DC stored as actual value,
AC harmonics stored as V_rms (= V_peak / √2). Verified and documented in
`normalization_chain.md` and `RMS-notes.md`.

### Input pipeline

```
SPICE .ckt  →  s2h  →  .ckt.hbl  →  hbl  →  .raw
```

`s2h` is a two-pass C/flex/bison compiler (`spice2hbl/`). It maps SPICE syntax to HBL
and enforces `Rsmin=1Ω` on all diode series resistance (configured via `s2h.cfg`).

### Implemented element models (`chanes/libelem/`)

| File | Element |
|------|---------|
| `mdsch.for` | Diode (SCHT charge-based model) |
| `lib0.for` | Lossless transmission line (LL0, frequency-domain) |
| `liblin.for` | Norton current source (`jdrive`), series voltage source (`emf`) |
| `lin.for` | Linear two-port |
| `clin.for` | Controlled linear source |
| `biptr.for` | Bipolar transistor |
| `libmod.for` | General model table lookup |
| `junc.for` | Junction capacitance |
| `indsv.for` | Inductor with series resistance |
| `cdiff.for`, `cpoly.for`, `curt.for`, `cusd.for` | Nonlinear capacitor/current variants |
| `icujunc.for`, `icupoly.for` | Nonlinear IC models |
| `mpl.for`, `poly5.for` | Polynomial nonlinear models |
| `ytab.for` | Y-parameter table (S-param import) |
| `stab.for`, `svlutl.for` | Stability, LU utilities |
| others | Supporting utilities |

24 element files total. Elements are registered in `main.for` and dispatched by type code.

---

## What HBFree does NOT do

| Missing capability | Impact | Effort to add |
|--------------------|--------|--------------|
| **Lossless T-line via s2h** | `di2`, `di4`, `di5` cannot be s2h-converted; T element blocks | Medium (s2h parser + mapping) |
| **Continuation method** | `di6` and similar hard-nonlinear circuits diverge | Medium–High |
| **Native SPICE model embedding** | Must re-implement all models in Fortran | High |
| **Verilog-A model embedding** | No VA interpreter interface | High (long-term) |
| **MESFET validation** | No test circuits; element exists but untested vs reference | Low (circuits) + Medium (tuning) |
| **AC small-signal analysis** | No linearization around HB solution | Medium |
| **Noise analysis** | No noise figure / phase-noise output | High |
| **Parametric sweep** | No built-in sweep; manual re-runs only | Medium |
| **Large circuits** | Static array limits (see below) | Medium (refactor) |
| **Dynamic memory** | All arrays are static Fortran COMMON blocks | High (restructure) |

---

## Validated vs ngspice

### Single-tone diode circuits (vs ngspice, `--patch-rs`)

| Circuit | Description | DC error | f1 error | Notes |
|---------|-------------|----------|----------|-------|
| di1 | half-wave, 1V | **0.27%** | 2.96% | reference benchmark |
| di1_lo | half-wave, low drive | — | — | converges, qualitatively matched |
| di1_s2 | half-wave, Rvs=0.1Ω | **0.26%** | 0.98% | used for normalization fix |
| di7 | Graetz bridge (4-diode) | — | — | converges; detailed comparison pending |
| di8 | complex bridge | — | — | converges; detailed comparison pending |

Full comparison tables: `model_validation_diode.md`.

### Two-tone diode circuits (vs ngspice two-tone FFT)

| Circuit | Description | DC error | IF error | T-line |
|---------|-------------|----------|----------|--------|
| di3 | mixer, no T-line, 1V+1V | 10.9% | **6.7%** | no |
| di2 | mixer + T-line, 1V+1V | 3.8% | **1.7%** | yes (LL0 vs SPICE T) |
| di4 | mixer + T-line, 30V LO | 7.4% | 92% (! large LO) | yes |

Full tables + analysis: `model_validation_diode_2tone.md`.

### Key model difference

HBFree uses the SCHT (charge-based) diode. ngspice uses the standard SPICE junction
diode. DC and fundamental agree within 1–11%; higher harmonics diverge predictably.
At large LO (di4), cross-product (intermod) predictions are unreliable. This is a
**model formulation difference**, not a solver bug.

---

## Not yet validated

| Item | Blocker | Next step |
|------|---------|-----------|
| di5 | No HBFree `.raw` (s2h T-element parse fails) | Fix s2h T-element parser |
| di6 | DIVERGED ("CAN NOT MAKE GOOD STEP") | Investigate convergence; implement continuation |
| di1_small | No SPICE `.ckt` (uses HBFree LIB0 T-line internally) | Write `.ckt`, validate |
| MESFET | No test circuits exist | Create `mesfet1.ckt`, run s2h, compare |
| T-line standalone | s2h cannot emit LL0 for T element | Fix s2h |
| di7, di8 detailed | Converge but no column-by-column comparison table | Run `spice_sim.py`, build tables |

---

## Array size limits — need to increase

All sizes are static Fortran COMMON blocks. Current values and risks:

| Constant | Value | Meaning | Risk if exceeded |
|----------|-------|---------|-----------------|
| `ISIZE_MAXNODE` | 200 | nodes before reduction | silent overwrite |
| `ISIZE_MAXVAR` | 1000 | reduced system size | silent overwrite |
| `MAXKN` | 20 | number of harmonics/frequencies | wrong FFT, overwrite |
| `MAXKN1` | 200 = 10×MAXKN | derivative frequency grid | dependent on MAXKN |
| `B1_SIZE` | 256 | time-domain FFT buffer | FFT aliasing |
| `B2_SIZE` | 4096 | sparse 2D FFT buffer | FFT aliasing |
| `MPOINT/NODEL` | 500 | element pointer table | overwrite (not PARAMETER) |
| `BUFFER` | 6000 DCOMPLEX | main sparse matrix | overwrite (magic number) |
| libelem local `Y(15,15)` | 15 | nodes per element | element overwrite |

**Action needed**: before running circuits larger than ~10 nodes or >10 harmonics,
verify no limit is hit. The safe approach is to:
1. Convert all magic-number sizes to named `PARAMETER` constants in `funcsize.i`
2. Verify all dependent arrays use the same constant (not independent copies)
3. Then increase the constants uniformly

See `analysis_array_sizes.md` for full dependency map. The most dangerous issue is
`ISIZE_MAXNODE = 200` colliding numerically with `MAXKN1 = 200` — different semantics,
same value, easy to confuse during refactor.

---

## Continuation method — need to implement

### Problem

`di6` (4-diode bridge, 2-tone) diverges at full drive level. The Newton solver
overshots and cannot recover ("CAN NOT MAKE GOOD STEP" after 18 iterations).
This is expected for circuits with strong nonlinearity where the initial guess
(linearized solution) is far from the true HB solution.

### What continuation does

Ramp a homotopy parameter λ from 0 → 1:
- λ=0: linear circuit (trivially solvable)
- λ=1: full nonlinear circuit (target)

At each λ step, the previous converged solution is the warm-start guess for the next.
Step size can be adaptive (halve on divergence, double on fast convergence).

### Implementation sketch

The cleanest hook is in `main.for` around the Newton loop. A new `.hb_options`
parameter `continuation=1` would enable it. The nonlinear element evaluations
(all in `libelem/`) must accept a λ scale on their current/charge output.

**Risk**: all 24 element files need a λ parameter passed through. The COMMON block
structure makes this non-trivial without a refactor.

**Alternative (simpler)**: source ramping only — scale `distof1`/`distof2` from
near-0 to full value, re-running HBFree at each step. Can be scripted externally
without touching Fortran. Less general but covers many practical cases.

---

## SPICE model embedding — design considerations

### Motivation

Today every device model must be hand-coded in Fortran. Standard SPICE model cards
(`.model nmos nmos level=14 ...`) with hundreds of parameters are impractical to port.

### Options

| Approach | Pros | Cons |
|----------|------|------|
| **Call ngspice as a library** (`libngspice`) | Full SPICE model coverage; no re-implementation | Requires C/Fortran interop; ngspice license (GPL); adds runtime dependency |
| **Tabulate I/V/Q from ngspice** (offline) | No runtime dependency; portable | Large tables; limited accuracy at extremes; workflow overhead |
| **Parse SPICE .model in s2h, evaluate inline** | Self-contained | Huge engineering effort; model coverage always lagging |
| **Verilog-A via OpenVAF/ADMS** | Industry-standard model format | Requires VA compiler toolchain; Fortran/C interop |

**Recommended first step**: tabular approach for a specific model family (e.g.
BSIM4 NMOS), driven by a concrete need. Full libngspice integration is the most
powerful but also the largest commitment — worthwhile only if HBFree is to become
a production tool for modern IC design.

---

## Build and test infrastructure

```bash
cd chanes && make           # gfortran, builds ../hbl
cd hbfree_exmpls && ../hbl di1.ckt.hbl   # run one circuit

# Full pipeline (SPICE → HBL → simulate → compare)
scripts/run_test.sh hbfree_exmpls/di1.ckt --ref hbfree_exmpls/results/di1.raw

# ngspice reference + phasor extraction
python3 scripts/spice_sim.py hbfree_exmpls/di3.ckt   # two-tone auto-detected

# Compare
python3 scripts/cmpraw.py di3.raw di3.spice.raw
```

No automated regression suite yet. All comparisons are run manually.

---

## Known bugs / open issues

| Issue | Severity | Status |
|-------|----------|--------|
| di6 diverges | Medium | Open — continuation method would fix |
| s2h T-element not parsed | Medium | Open — needs parser extension |
| Array size collisions (MAXNODE=MAXKN1=200) | Low (currently) | Open — refactor risk |
| `ISIZE_MAXVAR` not expressed from other constants | Low | Open — magic number |
| libelem local `15` not a PARAMETER | Low | Open — 37 files affected |
| di7/di8 detailed comparison not done | Low | Open |
