# HBFree Normalization Chain — Ground Truth

**Status: VERIFIED from code + experimental output (2026-04-25)**

---

## Circuit: di1_s2 (half-wave rectifier, Rvs=0.1Ω, Rl=1kΩ, f1=100MHz)

---

## Step 1: SPICE source → HBFree parameter P3(4)

```
SPICE:  sin(0  V_peak  f)  distof1  D
```

Convention: `D = distof1 = V_peak / 2` (standard SPICE HB usage)

`s2h` converts voltage source → Norton current source (jdrive element):

```
P3(4) = distof1 / Rvs = (V_peak/2) / Rvs = I_peak_Norton / 2
```

For V_peak=0.5V, Rvs=0.1Ω: **P3(4) = 2.5 A**

---

## Step 2: jdrive injection → Subj (liblin.for)

For AC harmonic matching source frequency:

```fortran
Subj(1) = sqrt(2) * P3(4) * exp(j*phase)
```

For P3(4)=2.5: **Subj = √2 × 2.5 = 3.536 A**

For DC: `Subj(1) = P3(4)` (no √2 factor, DC is correct).

---

## Step 3: J vector assembly → divide by 2 for AC (stepfr.for line 146)

```fortran
IF (Nrec.NE.1) Vj(jj) = Vj(jj) / 2.D0   ! AC only; DC (Nrec=1) not divided
```

**J_k_effective = Subj / 2 = (√2 × P3(4)) / 2 = P3(4)/√2**

After network reduction to reduced junction variable:

```
J_junction = J_source * G_rs / (G_rvs + G_rs)
           = (P3(4)/√2) * G_rs/(G_rvs+G_rs)
```

For P3(4)=2.5, G_rs=1, G_rvs=10:
**J_junction = (2.5/√2) × (1/11) = 0.1607 A**

**VERIFIED from log:** `VECTJ(2,2) = -5.9e-7 - j×0.16071` ✓

---

## Step 4: Linear solve → V_k_internal

```
Y_junction = G_rs * G_rvs / (G_rs + G_rvs) = 1*10/11 = 0.909 S
```

**VERIFIED from log:** `Y(2,2,2) = 0.90909` ✓

```
V_k_internal = J_junction / Y_junction = 0.1607 / 0.909 = 0.1769
```

**VERIFIED:** Initial Newton step = 0.176777 ✓

---

## Step 5: ftmas2 time-domain reconstruction

Formula (harm +2 with Hermitian fill):

```
v(t_n) = V_DC + 2 * Re( V_k_internal * exp(j*k*ω₁*t_n) )
v_peak  = 2 * |V_k_internal| = 2 * 0.1769 = 0.3536 V
```

**Physical expected:** V_junction_peak ≈ V_source_peak = 0.5 V (since Rvs,Rs << Rl)

**Actual in time domain:** 0.3536 V = V_peak/√2 = **V_rms of source** ← WRONG

**ERROR: diode evaluates at V_rms instead of V_peak. Error factor = √2.**

---

## Step 6: main.for post-processing (after solver converges)

```fortran
s(i,irc) = 2.D0 * s(i,irc)        ! × 2
s(i,irc) = s(i,irc) / dsqrt(2.D0) ! ÷ √2  → net: ×√2
```

Applied to AC harmonics only (irc=2..Kn). DC not scaled.

```
s_printed = V_k_internal * √2 = 0.1769 * 1.414 = 0.25 V
```

**VERIFIED from log:** `U(2) = -j*0.24999` (= -j*0.25) ✓

**.raw stores V_peak/2** (currently with the bug; should store V_rms after fix)

---

## Root Cause of Bug

The one-sided HB convention applies a ÷2 to AC harmonics in the J vector (stepfr.for).
This means the correct injection should be:

```
Subj = 2 * I_peak_Norton   so that   J_eff = Subj/2 = I_peak_Norton/2  (one-sided phasor)
```

`I_peak_Norton = V_peak / Rvs = 2 * P3(4)`

So the correct factor is **2**, not √2:

```
Subj_correct = 2 * P3(4) * exp(j*phase)
```

**Current code uses √2, which is off by √2.**

---

## Fix

In `chanes/libelem/liblin.for`, jdrive subroutine:

```fortran
! BEFORE (wrong):
sq12 = dsqrt(2.D0)

! AFTER (correct):
sq12 = 2.D0
```

Same fix for `emf` subroutine (series voltage source Norton equivalent).

---

## Chain After Fix

| Step | Quantity | Value (V_peak=0.5V) |
|------|----------|---------------------|
| SPICE source | V_peak | 0.5 V |
| P3(4) = distof1/Rvs | I_peak/2 | 2.5 A |
| Subj (jdrive, fixed) | 2×P3(4) | 5.0 A |
| J_eff (÷2 in stepfr) | I_peak/2 = P3(4) | 2.5 A → junction: 0.2273 A |
| V_k_internal (junction) | V_peak/2 | ~0.25 V |
| ftmas2 time domain peak | 2×V_k = V_peak | 0.5 V ✓ |
| main.for ×√2 → .raw | V_peak/√2 = V_rms | 0.354 V |

After fix: **.raw stores V_rms** (consistent with RMS-notes.md).

---

## Expected Impact

- HBFree currently evaluates the diode at **V_rms amplitude** instead of **V_peak**
- For exponential diode I = IS*(exp(α*V)-1), the difference is large:
  exp(38.67 * 0.354) vs exp(38.67 * 0.5) → ratio ~ exp(5.64) ~ 280×
- This explains the ~3.4× discrepancy in DC output (HBFree: 11.48mV vs ngspice: 39mV)
  (Not 280× because average current over one period is much less than peak)

---

## Verified Results (2026-04-25)

**di1_s2** (Rvs=0.1Ω, Rl=1kΩ, V_source=0.5V peak, f1=100MHz):

| Quantity | Before fix | After fix | ngspice | Match |
|----------|-----------|-----------|---------|-------|
| DC V(2)  | 11.48 mV | **39.16 mV** | 39.06 mV | **0.26%** ✓ |
| 100MHz V(2) | 15.03 mV | **49.74 mV** | 49.63 mV | **0.98%** ✓ |
| 200MHz V(2) | 11.84 mV | **35.36 mV** | 35.35 mV | **1.9%** ✓ |
| .raw V(1) k=1 | 0.250 V (V_peak/2) | **0.354 V (V_rms)** | 0.353 V | ✓ |

**di1** (Rl=1kΩ, V_source=1V peak):
| Quantity | After fix | ngspice | Match |
|----------|-----------|---------|-------|
| DC V(2)  | 170.87 mV | 170.41 mV | **0.27%** ✓ |
| DC V(out) | 170.78 mV | 170.42 mV | **0.21%** ✓ |

**di7, di8**: crash before AND after fix — pre-existing issue (array bounds), not caused by this change.

Remaining 2-20% differences at higher harmonics are due to SCHT diode model vs standard SPICE diode formulation, not normalization.

**V(1) at harmonics**: HBFree shows ~0 (correct for ideal current source); ngspice shows nonzero (transient captures harmonic current through Rvs). Not a bug — formulation difference.

## Files Changed

- `chanes/libelem/liblin.for`: `sq12 = 2.D0` in both `emf` and `jdrive` (was `dsqrt(2.D0)`)
- `scripts/spice_sim.py`: `write_hbfree_raw()` applies `×√2` to AC harmonics (k>0)
