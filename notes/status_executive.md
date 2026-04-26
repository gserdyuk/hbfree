# HBFree — Executive Status

**Date: 2026-04-26 (transient document — update as work progresses)**

---

## What HBFree is

HBFree is a specialized circuit simulator for **RF and microwave nonlinear circuit
analysis**. It computes the steady-state frequency-domain response of circuits
containing diodes, transistors, and transmission lines driven by one or two sinusoidal
tones. Typical use cases: mixer conversion gain, harmonic distortion, rectifier
efficiency.

It is an **open-source research tool** (GPL v2), written 1996–2004, now being
actively revived and validated. It is not a commercial product; it requires
engineering effort to use and to extend.

---

## Current capability — what it can do today

**Works reliably:**
- Single-tone harmonic balance: diode rectifiers, bridges (di1, di7, di8 families)
- Two-tone harmonic balance: diode mixers without transmission lines (di3)
- Two-tone with transmission line: small-signal drive (di2), validated at 1–4% DC error
- SPICE netlist input via `s2h` converter (most common elements)
- Comparison infrastructure: `spice_sim.py` + `cmpraw.py` vs ngspice reference

**Validated accuracy (vs ngspice):**
| Regime | DC accuracy | IF accuracy | Verdict |
|--------|-------------|-------------|---------|
| Single-tone, small signal | 0.2–0.3% | 1–3% | Production-quality |
| Two-tone, small signal, no T-line | ~11% | ~7% | Acceptable for feasibility |
| Two-tone, small signal, with T-line | ~4% | ~2% | Good for design guidance |
| Two-tone, large LO (30V), with T-line | ~7% (DC/LO) | ~92% (IF) | **Not reliable for IF** |

The accuracy difference between models is a known theoretical difference (HBFree uses
SCHT charge-based diode; ngspice uses standard SPICE junction diode), not a software
bug. Small-signal results are trustworthy; large-signal intermodulation requires
caution.

---

## What it cannot do today — and what that means for planning

| Limitation | Practical impact | Unblocks when |
|------------|-----------------|---------------|
| **No continuation method** | Hard-driven circuits (e.g. 4-diode bridge, di6) diverge without a result | After continuation is implemented |
| **T-line not in s2h** | Circuits with `T` (lossless transmission line) in SPICE netlist cannot be auto-converted | After s2h T-element fix |
| **Circuit size limits** | Max ~200 nodes, 20 harmonics, 1000 variables. Adequate for small RF sub-circuits; not for chip-level | After array-size refactor |
| **MESFET not validated** | Element exists but no test circuits; accuracy unknown | After test circuits created and validated |
| **No SPICE model embedding** | Models must be hand-coded in Fortran; modern foundry PDK models (BSIM, etc.) not usable | After model-embedding architecture decided |
| **No Verilog-A support** | Cannot use industry-standard behavioral models | Long-term; requires VA toolchain |
| **No noise / AC / parametric analysis** | Only HB steady-state; no noise figure, no S-parameter sweep | Separate feature development |

---

## Activity map — what depends on what

Activities are grouped into independent tracks. Tracks A and B can run in parallel.
Track C depends on Track A being substantially complete.

```
Track A — Validation (can start now, no code changes needed)
  A1. di7, di8 detailed comparison tables       (1–2 days)
  A2. Create MESFET test circuits                (1–2 days circuits + 1 day tuning)
  A3. di1_small: write .ckt, compare            (0.5 day)

Track B — Core simulator improvements (each item mostly independent)
  B1. Fix s2h T-element parser                   (1–2 days)
       → enables: di5 run, T-line standalone validation
  B2. Array size refactor (magic numbers → PARAMETER)   (2–3 days)
       → enables: larger circuits
  B3. Implement continuation method              (3–5 days)
       → enables: di6 and other hard-convergence circuits

Track C — Model embedding (strategic; starts after core is stable)
  C1. Evaluate libngspice vs tabular approach    (1 day — design only)
  C2. Implement chosen approach for one device   (1–4 weeks depending on choice)
  C3. Validate new model vs reference            (1–2 days per device)
```

**No external dependencies** within Tracks A and B — all work is self-contained
in the HBFree repository using gfortran + ngspice + Python.

Track C requires a decision on architecture (see below).

---

## Key decisions needed

### Decision 1: Continuation method — when?

**Context**: di6 (4-diode bridge) diverges. It is the only currently failing circuit
with a clear fix path. Without continuation, any similarly-driven circuit will also fail.

**Options**:
- (a) Implement now, before new test circuits — eliminates known blocker
- (b) Defer until di6 is specifically needed — saves ~3–5 days now
- (c) External source-ramping script — simpler, 1 day, covers ~80% of cases

Recommendation: option (c) first, option (a) later if production use requires it.

### Decision 2: Array size refactor — when?

**Context**: current limits (200 nodes, 20 harmonics) are adequate for the 11 existing
test circuits. A real mixer IC subnetwork may have 50–100 nodes; a PA match network
may need 30–50 harmonics.

**Options**:
- (a) Refactor now to remove magic numbers, then increase limits — 2–3 days
- (b) Defer until a specific circuit hits the limit — reactive

Recommendation: option (a) is low-risk and unblocks future work with no downside.

### Decision 3: SPICE model embedding — strategy

**Context**: foundry PDK models (BSIM4, PSP, etc.) have hundreds of parameters and
are available in SPICE `.model` format. Embedding them into HBFree is the primary
requirement for using HBFree on real IC design tasks.

**Options**:
- (a) **libngspice integration** — call ngspice as a shared library for all I/V/Q
  evaluations. Full model coverage. Adds ngspice runtime dependency. ~2–4 weeks
  of engineering work. Most powerful path if HBFree is to become a production tool.
- (b) **Tabular pre-characterization** — use ngspice to build I/V/Q tables offline;
  HBFree interpolates. 1 week to prototype. Limited to the table's sweep range.
  Suitable if a specific device family is the target.
- (c) **Verilog-A via OpenVAF** — compile VA models to shared libraries. Requires
  OpenVAF toolchain and Fortran/C interop layer. ~1–2 months. Future-proof but
  high upfront cost.

**This decision gates all production IC design work on HBFree.** It should be made
with a concrete target device in mind (what model family, what foundry).

---

## What is NOT on the critical path

The following items are in the notes/validation documents but do not block any
currently planned work:

- di5 (T-line fix needed, but di5 is identical to di4 except harmonic count — di4 is already validated)
- di1_small standalone T-line validation (informational only; LL0 element works in di2)
- Detailed di7/di8 comparison tables (circuits converge; detailed numbers are documentation, not blockers)

---

## Risks

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|-----------|
| di6-type divergence in a real target circuit | Medium | High | Implement continuation (B3) |
| Array overflow silently corrupts results | Low (current circuits) | High | Refactor (B2) before scaling up |
| SPICE model mismatch makes HBFree results unreliable | Medium (large signal) | Medium | Document known regime; decide on embedding (Decision 3) |
| Fortran legacy code fragility | Low | Medium | Incremental changes, run tests after each |

---

## Summary: three sentences for a status meeting

HBFree successfully simulates diode mixer circuits in harmonic balance with 1–11%
accuracy against ngspice for small-signal drive; large-signal intermodulation
predictions are not reliable due to a known model formulation difference.
The simulator is ready for feasibility studies on small diode circuits; two blockers
for broader use are (1) the T-line SPICE converter bug and (2) missing continuation
for hard-driven circuits — both are 1–3 day engineering tasks.
The strategic question that determines HBFree's applicability to IC design work is
whether and how to embed standard SPICE foundry models; that decision should be made
before committing significant development effort to other extensions.
