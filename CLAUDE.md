# HBFree — Harmonic Balance Simulator

Numerical simulator for nonlinear circuit analysis in the frequency domain.
Author: Gennady Serdyuk, 1996–2004. License: GPL v2.

## Components

| Directory | Language | Role |
|-----------|----------|------|
| `chanes/` | Fortran 77/90 (`.for`) | Main HB simulator |
| `spice2hbl/` | C + flex/bison | SPICE → HBL converter |
| `hbw/` | C++ | Windows wrapper |
| `hbfree_exmpls/` | — | Test circuits + reference results |
| `notes/` | — | Analysis docs |

## Build

```bash
cd chanes
make          # builds ../hbl  (gfortran)
make FC=intel # builds with ifx (Intel oneAPI)
```

Key compiler flags (gfortran):
```
-fdefault-real-8 -fdefault-double-8 -freal-4-real-8
-O2 -malign-double -fdec-char-conversions -fallow-argument-mismatch
```

## Run

```bash
cd hbfree_exmpls
../hbl di1.ckt.hbl
```

Output: `di1.raw` (SPICE-compatible), console log.

## Test

```bash
# Full pipeline: SPICE → HBL → simulate → compare
scripts/run_test.sh hbfree_exmpls/di1.ckt --ref hbfree_exmpls/results/di1.raw

# Compare two .raw files directly
python3 scripts/cmpraw.py ref.raw test.raw [--atol 1e-9] [--rtol 1e-3] [--mag]

# Simulate with ngspice and produce HBFree-compatible .raw + .tran.npz
python3 scripts/spice_sim.py [--patch-rs] [--format hbfree|npz|both] circuit.ckt
```

Converging tests (verified 2026-04-26, after normalization fix):

| Circuit | Status | Iterations | Notes |
|---------|--------|-----------|-------|
| di1 | CONVERGED | 13 | half-wave rectifier |
| di1_lo | CONVERGED | 3 | low-drive variant |
| di1_s2 | CONVERGED | 12 | small Rvs variant |
| di1_small | CONVERGED | 3 | reduced circuit |
| di2 | CONVERGED (step<1e-8) | 34 | fixed by normalization fix |
| di3 | CONVERGED (step<1e-8) | 25 | fixed by normalization fix |
| di4 | CONVERGED (step<1e-8) | 38 | fixed by normalization fix |
| di5 | s2h FAIL | — | lossless T-line element not supported by s2h |
| di6 | DIVERGED | 18 | "CAN NOT MAKE GOOD STEP" — remaining open issue |
| di7 | CONVERGED | 15 | |
| di8 | CONVERGED | 18 | |

di2/di3/di4 use a step-size stop criterion (max step < 1e-8) rather than the
"CONVERGED" banner; final error norms are < 1e-10, which is fully acceptable.

di5 requires a lossless transmission line (`T` element) not yet implemented in s2h.

## Scripts directory (`scripts/`)

| Script | Purpose |
|--------|---------|
| `cmpraw.py` | numdiff-style comparison of two HBFree `.raw` files with atol/rtol |
| `run_test.sh` | Full s2h → hbl → cmpraw pipeline wrapper |
| `spice_sim.py` | Run ngspice `.tran`, FFT results, write HBFree `.raw` + `.npz` |
| `make_conv_summary.py` | Extract convergence summary from `.ckt.run.log` → `.conv` golden file |
| `compare_iv.py` | Compare diode I-V: analytical curve vs ngspice DC sweep vs HBFree DC point |
| `plot_timedomain.py` | Overlay HBFree and ngspice waveforms in the time domain |

### spice_sim.py outputs
- `<stem>.spice.raw` — HBFree-format file with FFT phasors; compare with `cmpraw.py`
- `<stem>.tran.npz` — NumPy archive with time-domain data + phasor arrays per node

### make_conv_summary.py usage
```bash
python3 scripts/make_conv_summary.py di1.ckt.run.log -o results/di1.conv
```
Recognises three outcomes: `CONVERGED`, `STEP_LIMIT`, `DIVERGED`. Write `status: S2H_FAILED` manually for circuits that never ran.

### RS model mismatch (important)
`spice2hbl/s2h.cfg` has `Rsmin=1`, so all HBFree diode models use RS=1Ω minimum.  
SPICE default is RS=0. Use `--patch-rs` with `spice_sim.py` to add RS=1 to ngspice
models for a fairer comparison. Residual differences are due to model formulation
(HBFree SCHT ≠ standard SPICE diode).

## Key files

- `chanes/main.for` — main program + `wrtraw()` subroutine (writes .raw output)
- `chanes/funcsize.i` — array size constants (MAXKN=20, MAXKN1=200, etc.)
- `chanes/circuit.i` — COMMON block declarations
- `chanes/libelem/` — 24 element models (diodes, transmission lines, etc.)

## Array size constants

| Constant | Value | Meaning |
|----------|-------|---------|
| `ISIZE_MAXNODE` | 200 | Max circuit nodes |
| `ISIZE_MAXVAR` | 1000 | Max variables after reduction |
| `MAXKN` | 20 | Max frequencies |
| `MAXKN1` | 200 | Extended frequency grid (= 10×MAXKN) |
| `B1_SIZE` | 256 | FFT buffer |
| `B2_SIZE` | 4096 | Sparse 2D FFT buffer (= 16×B1_SIZE) |

See `notes/analysis_array_sizes.md` for full dependency map.

## Notes

- Source files are in `.for` format (SPAG-restructured Fortran 90 free-form)
- Original `.f` files preserved in git history (commit `931d4d7`)
- `wrtraw()` reads `*.ckt.nodes` for HBL→SPICE node name mapping;
  unmapped nodes are written as `#intNN`
- Known issue: static array limits may cause boundary violations on large circuits

## Claude instructions

- pls read all files in /notes directory - to be aware of project context
- always make plan before doing any programming task - new or on existing code and propose to review it
- always think how to make minimal effort. re-use same library, re-use same tool. make simple, not complex
- think on future - if we may need tool modifiction - foresee that and provide room on that, do not make code rigid, rather flexible
- try to not use magic numbers (you know on that) - use named constants

