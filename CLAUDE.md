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
cd hbfree_exmpls
python cmpraw.py   # compare .raw results with reference (results/)
```

Converging tests: **di1, di7, di8**  
Non-converging (known issue since 2004): di2, di3, di4, di6

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
