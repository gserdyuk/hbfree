# HBFree Amplitude & Normalization Convention

**Status: CORRECTED after normalization fix (2026-04-25)**

## Overview

HBFree source parameters use **I_peak/2 (one-sided phasor) amplitudes** — not RMS.
The `.raw` output file stores **V_rms** for AC harmonics.  The internal Harmonic
Balance solver works with one-sided complex phasors where `V_peak = 2 × |V_k|`.

**End-to-end chain (concrete example: distof1 = 0.5, Rvs = 1 Ω):**

```
distof1 = 0.5  →  P3(4) = 0.5  →  Subj = 2×0.5 = 1  →  ÷2 in stepfr
  →  J_eff = 0.5  →  V_k_internal = 0.5  →  2·|V_k| = 1 V peak  ✓
```

---

## 1. One-sided spectrum and the factor of 2

For a real signal, the two-sided Fourier spectrum is Hermitian:

```
V_{-k} = V_k*   (complex conjugate)
```

HBFree stores only the **positive** harmonics (k = 0, 1, 2, …, K_max).
The time-domain reconstruction from these one-sided phasors is:

```
v(t) = V_0  +  2·Re( V_1·exp(jω₁t) + V_2·exp(j2ω₁t) + … )
```

The factor of **2** comes from folding the negative-frequency side onto the
positive side.  Consequence:

```
V_peak  =  2 × |V_k|         (one-sided phasor stored in Znn)
|V_k|   =  V_peak / 2
```

**Confirmed in mdsch6 (Newton step limiter, `chanes/libelem/mdsch.for` line 176):**

```fortran
uold = uold + uold + dreal(Val(1,1))
!      ^^^^^^^^^^^ = 2 × Σ|Val(i,1)|  for i=2..Kn
! Total ≈ V_DC + 2×Σ|V_k|  = peak voltage estimate
```

---

## 2. IFFT: phasors → instantaneous time-domain voltages

`ftmas2` (`chanes/ftmas2.for`) converts frequency-domain phasors to time-domain
samples stored in buffer `B1`:

1. Copy sparse one-sided phasors from `Znn` into `B1[k]` for k = 1…K_max
2. Fill the negative-frequency half via Hermitian symmetry:
   `B1[N-k] = conj(B1[k])`
3. Call `harm(+2)` — the "direct" DFT:
   `B1[n] = Σ_k  B1[k] · exp(+j·2π·k·n/N)`

Because both +k and −k are present with conjugate symmetry, the result is:

```
B1[n] = V_k · exp(jωt_n)  +  V_k* · exp(−jωt_n)
       = 2·Re(V_k · exp(jωt_n))
       = V_peak · cos(ωt_n + φ)
```

**B1[n] is the actual instantaneous voltage** — no extra scale factor.

Concrete check: if `V_1 = 0.5` (one-sided phasor), `V_peak = 2×0.5 = 1 V`,
then `B1[n] = cos(2πn/N)` with peak = 1 V. ✓

---

## 3. Nonlinear element evaluation: instantaneous voltages required

`mdsch3` (`chanes/libelem/mdsch.for`) evaluates the diode at each time sample:

```fortran
u = B1(k, 1)                          ! actual instantaneous voltage
IF ( u*alfa .GT. argmax ) u = argmax/alfa
B1(k,1) = tok * (dexp(alfa*u) - 1.0)  ! I = IS·(exp(α·V)−1)
```

The diode I–V characteristic is **instantaneous** (not RMS, not average).
Therefore `u` **must** be the real instantaneous voltage, and the IFFT chain
above guarantees this. ✓

---

## 4. Source models: one-sided phasor (I_peak/2) convention

### jdrive (`chanes/libelem/liblin.for`, jdrive subroutine)

Norton current source (conductance G = P3(1), frequency P3(2), phase P3(3),
amplitude P3(4)):

```fortran
sq12 = 2.D0          ! fixed: was dsqrt(2.D0) — see normalization_chain.md
…
sq1 = sq12 * P3(4) * dcos(P3(3))   ! for Om ≠ 0
sq2 = sq12 * P3(4) * dsin(P3(3))
Subj(1) = dcmplx(sq1, sq2)
```

**P3(4) is I_peak/2** — the one-sided phasor amplitude of the Norton current.
The factor `2` gives the full-amplitude injection needed before `stepfr`'s ÷2:

```
Subj  =  2 × P3(4)  =  I_peak      (injected into J vector)
J_eff =  Subj / 2   =  P3(4)       (after stepfr ÷2 for AC)
```

DC path (`Om = 0`): `Subj = P3(4)` with no factor — DC is not divided by 2 in `stepfr`.

### emf (`chanes/libelem/liblin.for`, emf subroutine)

Series voltage source (resistance P3(1), frequency P3(2), phase P3(3),
amplitude P3(4)) represented as its Norton equivalent:

```fortran
sq1 = sq12 * P3(4) * dcos(P3(3)) / P3(1)
```

**P3(4) is V_peak/2** — one-sided phasor amplitude of the open-circuit voltage.
Norton current = `2 × P3(4) / R`, after stepfr ÷2: `J_eff = P3(4) / R`.

---

## 5. stepfr ÷2 for AC harmonics

`chanes/stepfr.for` divides all AC J-vector entries by 2 before the linear solve:

```fortran
IF (Nrec.NE.1) THEN
   DO jj = 1, k123
      Vj(jj) = Vj(jj) / 2.D0    ! AC only; DC (Nrec=1) is not divided
   ENDDO
ENDIF
```

This implements the one-sided phasor convention: injecting `2×P3(4)` then
dividing by 2 yields `J_eff = P3(4)`, which is the correct one-sided phasor
coefficient for the Norton current.

---

## 6. HB solver: internal phasor representation

The solver works with the one-sided phasor convention: `|V_k_internal| = V_peak / 2`.

With `J_eff = P3(4)` and load conductance G:

```
V_k_internal  =  J_eff / G  =  P3(4) / G
V_peak        =  2 × |V_k_internal|  =  2 × P3(4) / G
```

This is fed back into `ftmas2` for each Newton step, ensuring the nonlinear
evaluation always sees actual instantaneous voltages.

---

## 7. main.for post-processing: converting to .raw output

After the solver converges, `main.for` (lines 337–345) applies:

```fortran
DO irc = 2, Kn          ! only AC harmonics, not DC
  s(i,irc) = 2.D0 * s(i,irc)        ! step A: × 2
  s(i,irc) = s(i,irc) / dsqrt(2.D0) ! step B: ÷ √2   "BRING TO ACTUAL VALUE"
ENDDO
```

Net factor: `× 2 / √2 = × √2`.

Tracing with `V_k_internal = V_peak/2`:

| Step | Value | Meaning |
|------|-------|---------|
| After solver | `V_peak / 2` | one-sided phasor |
| After `× 2` | `V_peak` | actual peak voltage |
| After `÷ √2` | `V_peak / √2 = V_rms` | **RMS voltage** |

**The `.raw` file stores V_rms**, not V_peak.

---

## 8. Complete chain: concrete example with V_peak = 1 V, Rvs = 0.1 Ω

| Stage | Quantity | Value |
|-------|----------|-------|
| Physical signal | V_peak | 1.000 V |
| SPICE `distof1` | V_peak / 2 | 0.500 V |
| HBFree P3(4) = distof1 / Rvs | I_peak / 2 | 5.000 A |
| Source injection Subj = 2 × P3(4) | I_peak | 10.000 A |
| After stepfr ÷2: J_eff = P3(4) | I_peak / 2 | 5.000 A |
| Internal phasor V_k = J_eff / G_rvs | V_peak / 2 | 0.500 V |
| IFFT → B1[n] peak (2 × \|V_k\|) | V_peak | 1.000 V ✓ |
| Diode sees | actual instantaneous voltage | ✓ |
| .raw output (after × √2) | V_rms = V_peak / √2 | 0.707 V |

Reading the `.raw` file: `|S_k| = V_rms`.  To recover V_peak: `V_peak = √2 × |S_k|`.

---

## 9. SPICE `.ckt` → HBFree chain (via s2h)

`spice2hbl/s2h` converts a SPICE voltage source with `distof1` into a Norton
equivalent.  With source resistance `Rvs` (minimum 1 Ω from `s2h.cfg` `Rsmin`):

```
P3(4)  =  distof1 / Rvs        (Norton I_peak/2 parameter)
G      =  1 / Rvs              (Norton shunt conductance)
```

Resulting V_peak at the source node (assuming Rvs ≪ load):

```
V_k_internal  =  P3(4) / G  =  distof1
V_peak        =  2 × V_k_internal  =  2 × distof1
```

**Rule: `V_peak = 2 × distof1`**

In the `.raw` file:

```
|S_k|  =  V_rms  =  V_peak / √2  =  √2 × distof1
```

For SPICE `sin(0  A  f ...)` to match HBFree at the same drive level:

```
sin() amplitude  A  =  V_peak  =  2 × distof1
```

---

## 10. Summary table

| Quantity | Expression | For V_peak=1V, Rvs=0.1Ω |
|----------|------------|--------------------------|
| SPICE `sin()` amplitude | V_peak | 1.000 V |
| SPICE `distof1` | V_peak / 2 | 0.500 V |
| HBFree P3(4) (jdrive) | distof1 / Rvs = I_peak / 2 | 5.000 A |
| jdrive injection Subj | 2 × P3(4) = I_peak | 10.000 A |
| After stepfr ÷2: J_eff | P3(4) = I_peak / 2 | 5.000 A |
| Internal phasor \|V_k\| | V_peak / 2 | 0.500 V |
| Time-domain peak in B1 | V_peak | 1.000 V |
| `.raw` output \|S_k\| | V_rms = V_peak / √2 | 0.707 V |

---

## 11. Implication for spice_sim.py comparison

`spice_sim.py` computes phasors from ngspice transient data as:

```python
S_k = FFT[k] / N       # one-sided, |S_k| = V_peak / 2
```

To match HBFree `.raw` (which stores V_rms):

```python
S_k_hbfree_compatible = S_k * sqrt(2)   # converts V_peak/2 → V_rms
```

Applied in `write_hbfree_raw()` for AC harmonics (k > 0); DC unchanged.

---

## 12. plot_timedomain.py reconstruction

Since `.raw` stores V_rms phasors, time-domain reconstruction uses:

```python
v(t) = V_DC  +  √2 * Re( Σ_{k≥1}  V_k · exp(j · k · ω₁ · t) )
```

`|V_k| = V_rms`, so `√2 × |V_k| = V_peak`, recovering the correct amplitude.

---

*Analysed from: `chanes/libelem/liblin.for`, `chanes/libelem/mdsch.for`,
`chanes/ftmas2.for`, `chanes/harm.for`, `chanes/main.for`, `chanes/stepfr.for`*

*Fix applied 2026-04-25: `sq12 = 2.D0` (was `dsqrt(2.D0)`) in `jdrive` and `emf`.*
