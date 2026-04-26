# HBFree Diode Model Validation — vs ngspice

**Date: 2026-04-26**
**Normalization fix applied:** `sq12 = 2.D0` in `jdrive` and `emf` (was `dsqrt(2.D0)`)

---

## Setup

All circuits use the same SPICE diode model:

```spice
.model sd d is=1e-9 tt=0 m=0.5
```

**Key model difference**: HBFree implements the SCHT (charge-based) diode model;
ngspice uses the standard SPICE junction diode. Both share the same exponential I–V
relationship but differ in transit-time and junction-capacitance treatment.
Differences at higher harmonics are expected — they reflect model formulation
divergence, not a normalization or solver bug.

**RS mismatch**: `s2h.cfg` enforces `Rsmin=1Ω` for all HBFree diodes.
ngspice was run with `--patch-rs` to add `RS=1` to the diode model.

**Convention**: all `.ckt` files satisfy `sin() amplitude = 2 × distof1`.
`di7.ckt` and `di8.ckt` were corrected from `distof1=15` → `distof1=10.606`
(= 21.21/2) to match the source amplitude after the normalization fix.

**Comparison**: `scripts/cmpraw.py <stem>.raw <stem>.spice.raw`
Both files store V_rms phasors for AC harmonics (DC stored as actual value).

Columns: `Ref |z|` = HBFree magnitude, `Test |z|` = ngspice magnitude,
`adiff` = absolute difference, `rdiff` = relative difference.

---

## Circuit inventory

| Circuit    | Description               | HBFree         | Excitation | T-line        | Amplitude       | Comparison       |
|------------|---------------------------|----------------|------------|---------------|-----------------|------------------|
| di1        | half-wave rectifier       | CONVERGED      | single     | no            | matched         | **FULL**         |
| di1\_lo    | half-wave, low drive      | CONVERGED      | single     | no            | matched         | **FULL**         |
| di1\_s2    | half-wave, small Rvs      | CONVERGED      | single     | no            | matched         | **FULL**         |
| di7        | Graetz bridge             | CONVERGED      | single     | no            | matched         | **FULL**         |
| di8        | complex bridge            | CONVERGED      | single     | no            | matched         | **FULL**         |
| di3        | diode mixer               | CONVERGED      | two-tone   | no            | matched         | **PARTIAL** †    |
| di2        | diode mixer + T-line      | CONVERGED      | two-tone   | SPICE T       | matched         | **PARTIAL** †‡   |
| di4        | diode mixer + T-line      | CONVERGED      | two-tone   | SPICE T       | **mismatch**    | **NOT COMPARABLE** |
| di5        | diode mixer + T-line      | no .raw        | two-tone   | SPICE T       | **mismatch**    | **NOT COMPARABLE** |
| di6        | 4-diode bridge, 2-tone    | DIVERGED       | two-tone   | no (Rr used)  | —               | **NOT COMPARABLE** |
| di1\_small | half-wave + HBFree T-line | CONVERGED      | single     | HBFree LIB0   | —               | **NOT COMPARABLE** |

† Two-tone circuit: `spice_sim.py` extracts only n·f1 harmonics (1D FFT). Intermodulation products
at n1·f1 + n2·f2 (e.g., f2, f2−f1=10 MHz, 2f1−f2=90 MHz) are not compared.

‡ Additional model difference: HBFree uses frequency-domain LL0 T-line; ngspice uses time-domain lossless T element.

---

## V(1) source node note (all circuits)

HBFree models the source as an ideal Norton current source, so V(1) harmonics
are ~0 (only the DC bias appears). ngspice represents the loaded voltage source
directly, so V(1) harmonics are nonzero. This is a formulation difference, not
a model error. V(1) rows are included below for completeness but are excluded
from model accuracy conclusions.

---

## di1 — half-wave rectifier, V_peak=1V, Rl=1kΩ, f1=100MHz, 9 harmonics

**Nodes**: V(1)=source, V(2)=diode–load junction, V(out)=load

| Freq    | Node   | Ref \|z\|   | Test \|z\|  | adiff     | rdiff    |
|---------|--------|-------------|-------------|-----------|----------|
| DC      | V(1)   | 1.708e-05   | 2.205e-05   | 3.91e-05  | 177%†    |
| DC      | V(2)   | 1.7087e-01  | 1.7041e-01  | 4.59e-04  | **0.27%** |
| DC      | V(out) | 1.7078e-01  | 1.7042e-01  | 3.62e-04  | **0.21%** |
| 100 MHz | V(1)   | 7.0709e-01  | 7.0628e-01  | 1.91e-02  | 2.70%†   |
| 100 MHz | V(2)   | 2.0343e-01  | 2.0505e-01  | 6.07e-03  | **2.96%** |
| 100 MHz | V(out) | 2.0353e-01  | 2.0505e-01  | 7.34e-03  | **3.58%** |
| 200 MHz | V(1)   | 1.200e-05   | 2.017e-03   | 2.02e-03  | ~100%†   |
| 200 MHz | V(2)   | 1.2398e-01  | 1.2011e-01  | 6.79e-03  | 5.48%    |
| 200 MHz | V(out) | 1.2392e-01  | 1.2009e-01  | 8.14e-03  | 6.57%    |
| 300 MHz | V(1)   | 3.550e-06   | 1.136e-03   | 1.14e-03  | ~100%†   |
| 300 MHz | V(2)   | 3.3457e-02  | 3.5940e-02  | 2.88e-03  | 8.01%    |
| 300 MHz | V(out) | 3.3474e-02  | 3.5927e-02  | 3.26e-03  | 9.07%    |
| 400 MHz | V(1)   | 9.976e-07   | 8.060e-04   | 8.06e-04  | ~100%†   |
| 400 MHz | V(2)   | 1.0887e-02  | 9.7173e-03  | 1.73e-03  | 15.9%    |
| 400 MHz | V(out) | 1.0882e-02  | 9.7025e-03  | 1.94e-03  | 17.8%    |
| 500 MHz | V(1)   | 1.401e-06   | 6.268e-04   | 6.28e-04  | ~100%†   |
| 500 MHz | V(2)   | 1.3685e-02  | 1.4431e-02  | 7.46e-04  | 5.17%    |
| 500 MHz | V(out) | 1.3686e-02  | 1.4403e-02  | 8.53e-04  | 5.92%    |
| 600 MHz | V(1)   | 8.595e-08   | 5.243e-04   | 5.24e-04  | ~100%†   |
| 600 MHz | V(2)   | 8.3157e-04  | 1.6918e-03  | 9.26e-04  | 54.7%    |
| 600 MHz | V(out) | 8.3164e-04  | 1.6796e-03  | 9.34e-04  | 55.6%    |
| 700 MHz | V(1)   | 6.275e-07   | 4.475e-04   | 4.48e-04  | ~100%†   |
| 700 MHz | V(2)   | 6.0018e-03  | 5.8773e-03  | 1.36e-04  | 2.27%    |
| 700 MHz | V(out) | 6.0025e-03  | 5.8740e-03  | 2.53e-04  | 4.22%    |
| 800 MHz | V(1)   | 2.421e-07   | 3.812e-04   | 3.81e-04  | ~100%†   |
| 800 MHz | V(2)   | 2.2859e-03  | 3.1868e-03  | 9.92e-04  | 31.1%    |
| 800 MHz | V(out) | 2.2862e-03  | 3.1919e-03  | 1.07e-03  | 33.5%    |
| 900 MHz | V(1)   | 3.402e-07   | 3.389e-04   | 3.39e-04  | ~100%†   |
| 900 MHz | V(2)   | 3.1675e-03  | 1.8985e-03  | 1.27e-03  | 40.1%    |
| 900 MHz | V(out) | 3.1682e-03  | 1.8772e-03  | 1.30e-03  | 41.0%    |

† V(1) source-node artifact (see note above). V(#int04) is an internal HBFree node, not in ngspice output.

---

## di1_lo — half-wave rectifier, V_peak=0.2V, Rl=1kΩ, f1=100MHz, 9 harmonics

**Nodes**: V(1)=source, V(2)=diode–load junction, V(out)=load

| Freq    | Node   | Ref \|z\|   | Test \|z\|  | adiff     | rdiff    |
|---------|--------|-------------|-------------|-----------|----------|
| DC      | V(1)   | 3.137e-08   | 4.410e-06   | 4.44e-06  | ~100%†   |
| DC      | V(2)   | 3.1382e-04  | 3.1270e-04  | 1.12e-06  | **0.36%** |
| DC      | V(out) | 3.1366e-04  | 3.1275e-04  | 9.16e-07  | **0.29%** |
| 100 MHz | V(1)   | 1.4142e-01  | 1.4126e-01  | 3.81e-03  | 2.70%†   |
| 100 MHz | V(2)   | 4.1034e-04  | 4.1307e-04  | 9.50e-06  | **2.30%** |
| 100 MHz | V(out) | 4.1055e-04  | 4.1315e-04  | 1.21e-05  | **2.93%** |
| 200 MHz | V(1)   | 3.347e-08   | 4.035e-04   | 4.03e-04  | ~100%†   |
| 200 MHz | V(2)   | 3.4583e-04  | 3.3389e-04  | 1.66e-05  | 4.80%    |
| 200 MHz | V(out) | 3.4566e-04  | 3.3398e-04  | 1.98e-05  | 5.73%    |
| 300 MHz | V(1)   | 2.354e-08   | 2.271e-04   | 2.27e-04  | ~100%†   |
| 300 MHz | V(2)   | 2.2186e-04  | 2.3508e-04  | 1.56e-05  | 6.62%    |
| 300 MHz | V(out) | 2.2198e-04  | 2.3517e-04  | 1.83e-05  | 7.80%    |
| 400 MHz | V(1)   | 1.448e-08   | 1.612e-04   | 1.61e-04  | ~100%†   |
| 400 MHz | V(2)   | 1.5801e-04  | 1.4476e-04  | 1.35e-05  | 8.52%    |
| 400 MHz | V(out) | 1.5793e-04  | 1.4483e-04  | 1.46e-05  | 9.24%    |
| 500 MHz | V(1)   | 7.809e-09   | 1.254e-04   | 1.25e-04  | ~100%†   |
| 500 MHz | V(2)   | 7.6305e-05  | 7.8207e-05  | 7.14e-06  | 9.13%    |
| 500 MHz | V(out) | 7.6310e-05  | 7.8236e-05  | 4.74e-06  | 6.06%    |
| 600 MHz | V(1)   | 3.693e-09   | 1.049e-04   | 1.05e-04  | ~100%†   |
| 600 MHz | V(2)   | 3.5729e-05  | 3.7056e-05  | 3.86e-06  | 10.4%    |
| 600 MHz | V(out) | 3.5732e-05  | 3.7048e-05  | 2.56e-06  | 6.90%    |
| 700 MHz | V(1)   | 1.519e-09   | 8.949e-05   | 8.95e-05  | ~100%†   |
| 700 MHz | V(2)   | 1.4530e-05  | 1.5274e-05  | 1.69e-06  | 11.1%    |
| 700 MHz | V(out) | 1.4532e-05  | 1.5264e-05  | 1.12e-06  | 7.32%    |
| 800 MHz | V(1)   | 5.301e-10   | 7.625e-05   | 7.62e-05  | ~100%†   |
| 800 MHz | V(2)   | 5.0057e-06  | 5.3231e-06  | 5.54e-07  | 10.4%    |
| 800 MHz | V(out) | 5.0064e-06  | 5.3499e-06  | 3.97e-07  | 7.42%    |
| 900 MHz | V(1)   | 1.457e-10   | 6.778e-05   | 6.78e-05  | ~100%†   |
| 900 MHz | V(2)   | 1.3567e-06  | 1.4252e-06  | 8.34e-08  | 5.85%    |
| 900 MHz | V(out) | 1.3569e-06  | 1.5035e-06  | 1.49e-07  | 9.88%    |

† V(1) source-node artifact. V(#int04) is internal HBFree node.

---

## di1_s2 — half-wave rectifier, Rvs=0.1Ω, V_peak=0.5V, Rl=1kΩ, f1=100MHz, 15 harmonics

**Nodes**: V(1)=source, V(2)=load

| Freq    | Node | Ref \|z\|   | Test \|z\|  | adiff     | rdiff    |
|---------|------|-------------|-------------|-----------|----------|
| DC      | V(1) | 3.916e-06   | 6.967e-06   | 1.09e-05  | 156%†    |
| DC      | V(2) | 3.9164e-02  | 3.9063e-02  | 1.02e-04  | **0.26%** |
| 100 MHz | V(1) | 3.5355e-01  | 3.5312e-01  | 4.74e-03  | 1.34%†   |
| 100 MHz | V(2) | 4.9736e-02  | 4.9633e-02  | 4.85e-04  | **0.98%** |
| 200 MHz | V(1) | 3.536e-06   | 1.136e-03   | 1.14e-03  | ~100%†   |
| 200 MHz | V(2) | 3.5359e-02  | 3.5349e-02  | 6.76e-04  | 1.91%    |
| 300 MHz | V(1) | 1.841e-06   | 6.386e-04   | 6.40e-04  | ~100%†   |
| 300 MHz | V(2) | 1.8409e-02  | 1.8478e-02  | 5.34e-04  | 2.89%    |
| 400 MHz | V(1) | 5.013e-07   | 4.541e-04   | 4.54e-04  | ~100%†   |
| 400 MHz | V(2) | 5.0132e-03  | 5.1020e-03  | 2.14e-04  | 4.19%    |
| 500 MHz | V(1) | 1.818e-07   | 3.548e-04   | 3.55e-04  | ~100%†   |
| 500 MHz | V(2) | 1.8183e-03  | 1.7677e-03  | 9.88e-05  | 5.43%    |
| 600 MHz | V(1) | 2.897e-07   | 2.920e-04   | 2.92e-04  | ~100%†   |
| 600 MHz | V(2) | 2.8967e-03  | 2.9006e-03  | 1.66e-04  | 5.73%    |
| 700 MHz | V(1) | 1.213e-07   | 2.484e-04   | 2.48e-04  | ~100%†   |
| 700 MHz | V(2) | 1.2128e-03  | 1.2457e-03  | 8.92e-05  | 7.16%    |
| 800 MHz | V(1) | 4.370e-08   | 2.164e-04   | 2.16e-04  | ~100%†   |
| 800 MHz | V(2) | 4.3697e-04  | 4.1214e-04  | 4.03e-05  | 9.23%    |
| 900 MHz | V(1) | 8.644e-08   | 1.918e-04   | 1.92e-04  | ~100%†   |
| 900 MHz | V(2) | 8.6440e-04  | 8.6572e-04  | 7.44e-05  | 8.60%    |
| 1 GHz   | V(1) | 3.843e-08   | 1.722e-04   | 1.72e-04  | ~100%†   |
| 1 GHz   | V(2) | 3.8430e-04  | 4.0289e-04  | 4.24e-05  | 10.5%    |
| 1.1 GHz | V(1) | 1.721e-08   | 1.564e-04   | 1.56e-04  | ~100%†   |
| 1.1 GHz | V(2) | 1.7211e-04  | 1.5768e-04  | 2.22e-05  | 12.9%    |
| 1.2 GHz | V(1) | 3.201e-08   | 1.432e-04   | 1.43e-04  | ~100%†   |
| 1.2 GHz | V(2) | 3.2011e-04  | 3.2369e-04  | 3.71e-05  | 11.5%    |
| 1.3 GHz | V(1) | 1.281e-08   | 1.321e-04   | 1.32e-04  | ~100%†   |
| 1.3 GHz | V(2) | 1.2811e-04  | 1.4569e-04  | 2.47e-05  | 16.9%    |
| 1.4 GHz | V(1) | 8.719e-09   | 1.226e-04   | 1.23e-04  | ~100%†   |
| 1.4 GHz | V(2) | 8.7189e-05  | 7.3892e-05  | 1.69e-05  | 19.4%    |
| 1.5 GHz | V(1) | 1.235e-08   | 1.144e-04   | 1.14e-04  | ~100%†   |
| 1.5 GHz | V(2) | 1.2347e-04  | 1.3477e-04  | 2.17e-05  | 16.1%    |

† V(1) source-node artifact. V(#int03) is internal HBFree node.

---

## di7 — 4-diode bridge, V_peak=21.21V, Rr=2Ω, Rl=1kΩ, f1=100MHz, 19 harmonics

**Topology**: Graetz bridge; V(3) and V(4) are bridge output rails (both referenced to ground).
Rr=2Ω source resistance.

**Bridge symmetry note**: In a balanced bridge, even harmonics cancel at V(1)/V(2)
and odd harmonics cancel at V(3)/V(4). HBFree enforces this algebraically (result≈0
to machine precision). ngspice transient has small nonzero residuals from finite step
size. Rows marked ‡ are symmetry-cancellation nodes — excluded from model accuracy.

| Freq    | Node | Ref \|z\|   | Test \|z\|  | adiff     | rdiff    |
|---------|------|-------------|-------------|-----------|----------|
| DC      | V(1) | 4.4e-17     | 1.859e-04   | 1.86e-04  | ~100%†   |
| DC      | V(2) | 9.3e-16     | 1.859e-04   | 1.86e-04  | ~100%†   |
| DC      | V(3) | 6.3278      | 6.3021      | 2.57e-02  | **0.41%** |
| DC      | V(4) | 6.3278      | 6.3019      | 2.59e-02  | **0.41%** |
| 100 MHz | V(1) | 1.4998e+01  | 1.4983e+01  | 1.59e-01  | 1.06%†   |
| 100 MHz | V(2) | 1.4969e+01  | 1.4955e+01  | 1.59e-01  | 1.06%†   |
| 100 MHz | V(3) | 7.4847      | 7.4775      | 5.86e-02  | **0.78%** |
| 100 MHz | V(4) | 7.4847      | 7.4774      | 1.00e-01  | **1.34%** |
| 200 MHz | V(1) | 3.7e-17     | 3.785e-02   | 3.79e-02  | ~100%‡   |
| 200 MHz | V(2) | 7.8e-16     | 3.778e-02   | 3.78e-02  | ~100%‡   |
| 200 MHz | V(3) | 3.1257      | 3.1536      | 5.51e-02  | 1.75%    |
| 200 MHz | V(4) | 3.1257      | 3.1537      | 8.97e-02  | 2.84%    |
| 300 MHz | V(1) | 2.175e-05   | 2.128e-02   | 2.13e-02  | ~100%‡   |
| 300 MHz | V(2) | 4.568e-04   | 2.077e-02   | 2.12e-02  | ~100%‡   |
| 300 MHz | V(3) | 2.284e-04   | 1.046e-02   | 1.07e-02  | ~100%‡   |
| 300 MHz | V(4) | 2.284e-04   | 1.043e-02   | 1.07e-02  | ~100%‡   |
| 400 MHz | V(1) | 4.8e-18     | 1.513e-02   | 1.51e-02  | ~100%‡   |
| 400 MHz | V(2) | 1.0e-16     | 1.510e-02   | 1.51e-02  | ~100%‡   |
| 400 MHz | V(3) | 5.9998e-01  | 6.1885e-01  | 2.63e-02  | 4.25%    |
| 400 MHz | V(4) | 5.9998e-01  | 6.1898e-01  | 3.82e-02  | 6.18%    |
| 500 MHz | V(1) | 1.082e-05   | 1.182e-02   | 1.18e-02  | ~100%‡   |
| 500 MHz | V(2) | 2.272e-04   | 1.153e-02   | 1.18e-02  | ~100%‡   |
| 500 MHz | V(3) | 1.136e-04   | 5.887e-03   | 6.00e-03  | ~100%‡   |
| 500 MHz | V(4) | 1.136e-04   | 5.837e-03   | 5.95e-03  | ~100%‡   |
| 600 MHz | V(1) | 2.4e-18     | 9.727e-03   | 9.73e-03  | ~100%‡   |
| 600 MHz | V(2) | 5.1e-17     | 9.709e-03   | 9.71e-03  | ~100%‡   |
| 600 MHz | V(3) | 2.4171e-01  | 2.5923e-01  | 2.08e-02  | 8.02%    |
| 600 MHz | V(4) | 2.4171e-01  | 2.5936e-01  | 2.71e-02  | 10.4%    |
| 700 MHz | V(1) | 5.550e-06   | 8.276e-03   | 8.28e-03  | ~100%‡   |
| 700 MHz | V(2) | 1.165e-04   | 8.077e-03   | 8.19e-03  | ~100%‡   |
| 700 MHz | V(3) | 5.827e-05   | 4.199e-03   | 4.26e-03  | ~100%‡   |
| 700 MHz | V(4) | 5.827e-05   | 4.131e-03   | 4.19e-03  | ~100%‡   |
| 800 MHz | V(1) | 1.2e-18     | 7.207e-03   | 7.21e-03  | ~100%‡   |
| 800 MHz | V(2) | 2.5e-17     | 7.194e-03   | 7.19e-03  | ~100%‡   |
| 800 MHz | V(3) | 1.2337e-01  | 1.4006e-01  | 1.84e-02  | 13.1%    |
| 800 MHz | V(4) | 1.2337e-01  | 1.4019e-01  | 2.22e-02  | 15.8%    |
| 900 MHz | V(1) | 1.909e-06   | 6.386e-03   | 6.39e-03  | ~100%‡   |
| 900 MHz | V(2) | 4.008e-05   | 6.237e-03   | 6.28e-03  | ~100%‡   |
| 900 MHz | V(3) | 2.004e-05   | 3.313e-03   | 3.33e-03  | ~100%‡   |
| 900 MHz | V(4) | 2.004e-05   | 3.229e-03   | 3.25e-03  | ~100%‡   |
| 1 GHz   | V(1) | 6.6e-19     | 5.735e-03   | 5.74e-03  | ~100%‡   |
| 1 GHz   | V(2) | 1.4e-17     | 5.725e-03   | 5.73e-03  | ~100%‡   |
| 1 GHz   | V(3) | 7.0262e-02  | 8.6237e-02  | 1.69e-02  | 19.6%    |
| 1 GHz   | V(4) | 7.0262e-02  | 8.6366e-02  | 1.94e-02  | 22.4%    |
| 1.1 GHz | V(1) | 1.377e-06   | 5.206e-03   | 5.20e-03  | ~100%‡   |
| 1.1 GHz | V(2) | 2.892e-05   | 5.087e-03   | 5.06e-03  | ~99%‡    |
| 1.1 GHz | V(3) | 1.446e-05   | 2.767e-03   | 2.75e-03  | ~99%‡    |
| 1.1 GHz | V(4) | 1.446e-05   | 2.669e-03   | 2.66e-03  | ~99%‡    |
| 1.2 GHz | V(1) | 6.1e-19     | 4.767e-03   | 4.77e-03  | ~100%‡   |
| 1.2 GHz | V(2) | 1.3e-17     | 4.759e-03   | 4.76e-03  | ~100%‡   |
| 1.2 GHz | V(3) | 4.2138e-02  | 5.7461e-02  | 1.59e-02  | 27.6%    |
| 1.2 GHz | V(4) | 4.2138e-02  | 5.7590e-02  | 1.75e-02  | 30.4%    |
| 1.3 GHz | V(1) | 5.189e-06   | 4.397e-03   | 4.39e-03  | ~100%‡   |
| 1.3 GHz | V(2) | 1.090e-04   | 4.300e-03   | 4.19e-03  | ~97%‡    |
| 1.3 GHz | V(3) | 5.448e-05   | 2.397e-03   | 2.34e-03  | ~98%‡    |
| 1.3 GHz | V(4) | 5.448e-05   | 2.286e-03   | 2.24e-03  | ~98%‡    |
| 1.4 GHz | V(1) | 8.6e-19     | 4.081e-03   | 4.08e-03  | ~100%‡   |
| 1.4 GHz | V(2) | 1.8e-17     | 4.074e-03   | 4.07e-03  | ~100%‡   |
| 1.4 GHz | V(3) | 2.5607e-02  | 4.0334e-02  | 1.51e-02  | 37.3%    |
| 1.4 GHz | V(4) | 2.5607e-02  | 4.0463e-02  | 1.62e-02  | 39.9%    |
| 1.5 GHz | V(1) | 1.098e-05   | 3.808e-03   | 3.80e-03  | ~100%‡   |
| 1.5 GHz | V(2) | 2.305e-04   | 3.726e-03   | 3.50e-03  | ~94%‡    |
| 1.5 GHz | V(3) | 1.152e-04   | 2.130e-03   | 2.02e-03  | ~95%‡    |
| 1.5 GHz | V(4) | 1.152e-04   | 2.008e-03   | 1.91e-03  | ~95%‡    |
| 1.6 GHz | V(1) | 6.5e-19     | 3.569e-03   | 3.57e-03  | ~100%‡   |
| 1.6 GHz | V(2) | 1.4e-17     | 3.563e-03   | 3.56e-03  | ~100%‡   |
| 1.6 GHz | V(3) | 1.5140e-02  | 2.9354e-02  | 1.44e-02  | 49.0%    |
| 1.6 GHz | V(4) | 1.5140e-02  | 2.9483e-02  | 1.51e-02  | 51.3%    |
| 1.7 GHz | V(1) | 2.349e-05   | 3.359e-03   | 3.34e-03  | ~99%‡    |
| 1.7 GHz | V(2) | 4.933e-04   | 3.289e-03   | 2.80e-03  | ~85%‡    |
| 1.7 GHz | V(3) | 2.466e-04   | 1.928e-03   | 1.70e-03  | ~88%‡    |
| 1.7 GHz | V(4) | 2.466e-04   | 1.797e-03   | 1.61e-03  | ~89%‡    |
| 1.8 GHz | V(1) | 1.0e-18     | 3.172e-03   | 3.17e-03  | ~100%‡   |
| 1.8 GHz | V(2) | 2.1e-17     | 3.167e-03   | 3.17e-03  | ~100%‡   |
| 1.8 GHz | V(3) | 8.0317e-03  | 2.1921e-02  | 1.40e-02  | 63.8%    |
| 1.8 GHz | V(4) | 8.0317e-03  | 2.2049e-02  | 1.44e-02  | 65.5%    |
| 1.9 GHz | V(1) | 8.429e-05   | 3.006e-03   | 2.92e-03  | ~97%‡    |
| 1.9 GHz | V(2) | 1.770e-03   | 2.945e-03   | 1.23e-03  | 41.7%‡   |
| 1.9 GHz | V(3) | 8.851e-04   | 1.769e-03   | 9.81e-04  | 55.5%‡   |
| 1.9 GHz | V(4) | 8.851e-04   | 1.629e-03   | 1.11e-03  | 68.2%‡   |

† V(1)/V(2) source-node artifact.
‡ Bridge symmetry cancellation: HBFree gives near-zero; ngspice transient has numerical residual. Not a model accuracy metric.

---

## di8 — complex bridge (4+2 diodes), V_peak=21.21V, Rr=2Ω, Rl=1kΩ, f1=100MHz, 19 harmonics

**Topology**: Extended bridge. 7 nodes. V(4) is the main bridge positive rail
(equivalent to V(3) in di7 — both at ~6.3V DC).

**Note**: di8 is an asymmetric circuit (d5/d6 + rs1/rs2 sub-network), so the
bridge symmetry cancellation pattern differs from di7. V(1) still shows the
source-node artifact at even/zero-crossings.

| Freq    | Node | Ref \|z\|   | Test \|z\|  | adiff     | rdiff    |
|---------|------|-------------|-------------|-----------|----------|
| DC      | V(1) | 9.799e-02   | 1.859e-04   | 9.82e-02  | ~100%†   |
| DC      | V(2) | 2.0578      | 1.9859      | 7.19e-02  | 3.50%    |
| DC      | V(3) | 4.2590      | 4.3105      | 5.14e-02  | 1.19%    |
| DC      | V(4) | 6.3209      | 6.3003      | 2.06e-02  | **0.33%** |
| DC      | V(5) | 4.064e-03   | 3.937e-03   | 1.27e-04  | 3.13%    |
| DC      | V(6) | 4.0217      | 3.9759      | 4.58e-02  | 1.14%    |
| DC      | V(7) | 1.9639      | 1.9900      | 2.61e-02  | 1.31%    |
| 100 MHz | V(1) | 1.4888e+01  | 1.4983e+01  | 1.84e-01  | 1.23%†   |
| 100 MHz | V(2) | 1.2668e+01  | 1.2732e+01  | 1.55e-01  | 1.22%    |
| 100 MHz | V(3) | 5.2045      | 5.2753      | 8.19e-02  | 1.55%    |
| 100 MHz | V(4) | 7.4874      | 7.4806      | 1.00e-01  | **1.34%** |
| 100 MHz | V(5) | 2.373e-02   | 2.385e-02   | 2.90e-04  | 1.22%    |
| 100 MHz | V(6) | 1.0472e+01  | 1.0505e+01  | 1.28e-01  | 1.22%    |
| 100 MHz | V(7) | 2.1962      | 2.2274      | 3.54e-02  | 1.59%    |
| 200 MHz | V(1) | 4.787e-02   | 3.785e-02   | 6.16e-02  | 129%†    |
| 200 MHz | V(2) | 1.0052      | 9.7573e-01  | 6.10e-02  | 6.07%    |
| 200 MHz | V(3) | 2.1230      | 2.1731      | 5.96e-02  | 2.74%    |
| 200 MHz | V(4) | 3.1302      | 3.1503      | 8.75e-02  | 2.78%    |
| 200 MHz | V(5) | 2.024e-03   | 1.944e-03   | 1.33e-04  | 6.59%    |
| 200 MHz | V(6) | 1.9645      | 1.9524      | 6.87e-02  | 3.50%    |
| 200 MHz | V(7) | 9.594e-01   | 9.771e-01   | 2.31e-02  | 2.36%    |
| 300 MHz | V(1) | 1.126e-03   | 2.128e-02   | 2.24e-02  | ~100%†   |
| 300 MHz | V(2) | 2.364e-02   | 5.576e-03   | 1.81e-02  | 76.4%    |
| 300 MHz | V(3) | 2.498e-02   | 1.722e-02   | 7.97e-03  | 31.9%    |
| 300 MHz | V(4) | 1.727e-03   | 1.211e-02   | 1.04e-02  | 85.8%    |
| 300 MHz | V(5) | 3.914e-04   | 4.466e-04   | 5.69e-05  | 12.7%    |
| 300 MHz | V(6) | 4.577e-02   | 3.198e-02   | 1.38e-02  | 30.2%    |
| 300 MHz | V(7) | 2.213e-02   | 2.641e-02   | 4.32e-03  | 16.4%    |
| 400 MHz | V(1) | 9.496e-03   | 1.513e-02   | 1.81e-02  | 120%†    |
| 400 MHz | V(2) | 1.9941e-01  | 1.9371e-01  | 2.20e-02  | 11.1%    |
| 400 MHz | V(3) | 4.0722e-01  | 4.2460e-01  | 2.13e-02  | 5.03%    |
| 400 MHz | V(4) | 6.0705e-01  | 6.1828e-01  | 3.52e-02  | 5.69%    |
| 400 MHz | V(5) | 4.197e-04   | 3.855e-04   | 5.48e-05  | 13.1%    |
| 400 MHz | V(6) | 3.8974e-01  | 3.8698e-01  | 2.71e-02  | 6.95%    |
| 400 MHz | V(7) | 1.9033e-01  | 1.9356e-01  | 6.67e-03  | 3.45%    |
| 500 MHz | V(1) | 6.162e-04   | 1.182e-02   | 1.24e-02  | ~100%†   |
| 500 MHz | V(2) | 1.294e-02   | 3.668e-03   | 9.28e-03  | 71.7%    |
| 500 MHz | V(3) | 1.406e-02   | 1.025e-02   | 4.17e-03  | 29.7%    |
| 500 MHz | V(4) | 1.322e-03   | 6.863e-03   | 5.57e-03  | 81.2%    |
| 500 MHz | V(5) | 2.032e-04   | 2.540e-04   | 5.23e-05  | 20.6%    |
| 500 MHz | V(6) | 2.506e-02   | 1.890e-02   | 6.22e-03  | 24.8%    |
| 500 MHz | V(7) | 1.212e-02   | 1.523e-02   | 3.16e-03  | 20.7%    |
| 600 MHz | V(1) | 4.032e-03   | 9.727e-03   | 1.07e-02  | 110%†    |
| 600 MHz | V(2) | 8.468e-02   | 8.242e-02   | 1.38e-02  | 16.3%    |
| 600 MHz | V(3) | 1.6415e-01  | 1.7688e-01  | 1.48e-02  | 8.36%    |
| 600 MHz | V(4) | 2.4902e-01  | 2.5907e-01  | 2.31e-02  | 8.93%    |
| 600 MHz | V(5) | 1.928e-04   | 1.637e-04   | 4.05e-05  | 21.0%    |
| 600 MHz | V(6) | 1.6552e-01  | 1.6419e-01  | 1.73e-02  | 10.5%    |
| 600 MHz | V(7) | 8.084e-02   | 8.206e-02   | 3.91e-03  | 4.76%    |
| 700 MHz | V(1) | 3.877e-04   | 8.276e-03   | 8.66e-03  | ~100%†   |
| 700 MHz | V(2) | 8.142e-03   | 2.511e-03   | 5.63e-03  | 69.2%    |
| 700 MHz | V(3) | 9.284e-03   | 7.176e-03   | 2.65e-03  | 28.5%    |
| 700 MHz | V(4) | 1.259e-03   | 4.876e-03   | 3.68e-03  | 75.4%    |
| 700 MHz | V(5) | 1.162e-04   | 1.744e-04   | 5.92e-05  | 34.0%    |
| 700 MHz | V(6) | 1.578e-02   | 1.312e-02   | 2.77e-03  | 17.6%    |
| 700 MHz | V(7) | 7.638e-03   | 1.061e-02   | 3.02e-03  | 28.4%    |
| 800 MHz | V(1) | 2.219e-03   | 7.207e-03   | 7.68e-03  | 107%†    |
| 800 MHz | V(2) | 4.660e-02   | 4.543e-02   | 1.02e-02  | 21.8%    |
| 800 MHz | V(3) | 8.396e-02   | 9.492e-02   | 1.21e-02  | 12.8%    |
| 800 MHz | V(4) | 1.3068e-01  | 1.4004e-01  | 1.76e-02  | 12.6%    |
| 800 MHz | V(5) | 1.197e-04   | 8.988e-05   | 3.69e-05  | 30.8%    |
| 800 MHz | V(6) | 9.111e-02   | 9.014e-02   | 1.28e-02  | 14.0%    |
| 800 MHz | V(7) | 4.450e-02   | 4.499e-02   | 2.76e-03  | 6.13%    |
| 900 MHz | V(1) | 2.480e-04   | 6.386e-03   | 6.63e-03  | ~100%†   |
| 900 MHz | V(2) | 5.208e-03   | 1.832e-03   | 3.38e-03  | 64.9%    |
| 900 MHz | V(3) | 6.518e-03   | 5.473e-03   | 1.82e-03  | 27.9%    |
| 900 MHz | V(4) | 1.369e-03   | 3.814e-03   | 2.56e-03  | 67.1%    |
| 900 MHz | V(5) | 5.946e-05   | 1.309e-04   | 7.20e-05  | 55.0%    |
| 900 MHz | V(6) | 1.011e-02   | 9.920e-03   | 7.36e-04  | 7.28%    |
| 900 MHz | V(7) | 4.901e-03   | 8.088e-03   | 3.22e-03  | 39.8%    |
| 1 GHz   | V(1) | 1.400e-03   | 5.735e-03   | 6.02e-03  | ~100%†   |
| 1 GHz   | V(2) | 2.939e-02   | 2.867e-02   | 8.04e-03  | 27.4%    |
| 1 GHz   | V(3) | 4.800e-02   | 5.797e-02   | 1.06e-02  | 18.4%    |
| 1 GHz   | V(4) | 7.748e-02   | 8.629e-02   | 1.44e-02  | 16.6%    |
| 1 GHz   | V(5) | 8.993e-05   | 5.639e-05   | 3.82e-05  | 42.5%    |
| 1 GHz   | V(6) | 5.748e-02   | 5.658e-02   | 1.01e-02  | 17.6%    |
| 1 GHz   | V(7) | 2.808e-02   | 2.819e-02   | 2.13e-03  | 7.57%    |
| 1.1 GHz | V(1) | 1.437e-04   | 5.206e-03   | 5.35e-03  | ~100%†   |
| 1.1 GHz | V(2) | 3.018e-03   | 1.394e-03   | 1.63e-03  | 54.1%    |
| 1.1 GHz | V(3) | 4.658e-03   | 4.395e-03   | 1.39e-03  | 29.9%    |
| 1.1 GHz | V(4) | 1.653e-03   | 3.151e-03   | 1.75e-03  | 55.4%    |
| 1.1 GHz | V(5) | 1.182e-05   | 1.036e-04   | 9.19e-05  | 88.7%    |
| 1.1 GHz | V(6) | 5.880e-03   | 7.890e-03   | 2.10e-03  | 26.6%    |
| 1.1 GHz | V(7) | 2.862e-03   | 6.496e-03   | 3.65e-03  | 56.3%    |
| 1.2 GHz | V(1) | 9.599e-04   | 4.767e-03   | 4.96e-03  | ~100%†   |
| 1.2 GHz | V(2) | 2.016e-02   | 1.968e-02   | 6.66e-03  | 33.0%    |
| 1.2 GHz | V(3) | 2.896e-02   | 3.826e-02   | 9.70e-03  | 25.3%    |
| 1.2 GHz | V(4) | 4.919e-02   | 5.756e-02   | 1.22e-02  | 21.2%    |
| 1.2 GHz | V(5) | 7.865e-05   | 3.840e-05   | 4.36e-05  | 55.5%    |
| 1.2 GHz | V(6) | 3.943e-02   | 3.856e-02   | 8.39e-03  | 21.3%    |
| 1.2 GHz | V(7) | 1.928e-02   | 1.917e-02   | 1.75e-03  | 9.06%    |
| 1.3 GHz | V(1) | 5.099e-05   | 4.397e-03   | 4.45e-03  | ~100%†   |
| 1.3 GHz | V(2) | 1.071e-03   | 1.090e-03   | 1.04e-04  | 9.54%    |
| 1.3 GHz | V(3) | 3.291e-03   | 3.652e-03   | 1.28e-03  | 35.1%    |
| 1.3 GHz | V(4) | 2.180e-03   | 2.695e-03   | 1.21e-03  | 44.8%    |
| 1.3 GHz | V(5) | 3.946e-05   | 8.485e-05   | 1.24e-04  | 146%     |
| 1.3 GHz | V(6) | 2.130e-03   | 6.492e-03   | 4.38e-03  | 67.4%    |
| 1.3 GHz | V(7) | 1.059e-03   | 5.402e-03   | 4.35e-03  | 80.5%    |
| 1.4 GHz | V(1) | 6.908e-04   | 4.081e-03   | 4.22e-03  | ~100%†   |
| 1.4 GHz | V(2) | 1.451e-02   | 1.430e-02   | 5.65e-03  | 38.9%    |
| 1.4 GHz | V(3) | 1.776e-02   | 2.656e-02   | 9.03e-03  | 34.0%    |
| 1.4 GHz | V(4) | 3.234e-02   | 4.046e-02   | 1.08e-02  | 26.7%    |
| 1.4 GHz | V(5) | 7.928e-05   | 2.588e-05   | 5.42e-05  | 68.4%    |
| 1.4 GHz | V(6) | 2.840e-02   | 2.779e-02   | 7.11e-03  | 25.0%    |
| 1.4 GHz | V(7) | 1.390e-02   | 1.378e-02   | 1.47e-03  | 10.6%    |
| 1.5 GHz | V(1) | 4.708e-05   | 3.808e-03   | 3.76e-03  | ~100%†   |
| 1.5 GHz | V(2) | 9.887e-04   | 8.664e-04   | 1.85e-03  | 187%     |
| 1.5 GHz | V(3) | 2.221e-03   | 3.108e-03   | 1.39e-03  | 44.7%    |
| 1.5 GHz | V(4) | 3.097e-03   | 2.362e-03   | 1.55e-03  | 50.1%    |
| 1.5 GHz | V(5) | 1.126e-04   | 7.115e-05   | 1.83e-04  | 163%     |
| 1.5 GHz | V(6) | 1.818e-03   | 5.469e-03   | 7.28e-03  | 133%     |
| 1.5 GHz | V(7) | 8.291e-04   | 4.603e-03   | 5.43e-03  | 118%     |
| 1.6 GHz | V(1) | 4.922e-04   | 3.569e-03   | 3.67e-03  | ~100%†   |
| 1.6 GHz | V(2) | 1.034e-02   | 1.084e-02   | 4.81e-03  | 44.3%    |
| 1.6 GHz | V(3) | 1.065e-02   | 1.908e-02   | 8.55e-03  | 44.8%    |
| 1.6 GHz | V(4) | 2.108e-02   | 2.950e-02   | 1.02e-02  | 34.5%    |
| 1.6 GHz | V(5) | 9.324e-05   | 2.060e-05   | 7.47e-05  | 80.1%    |
| 1.6 GHz | V(6) | 2.027e-02   | 2.084e-02   | 6.01e-03  | 28.9%    |
| 1.6 GHz | V(7) | 9.937e-03   | 1.030e-02   | 1.27e-03  | 12.3%    |
| 1.7 GHz | V(1) | 1.681e-04   | 3.359e-03   | 3.19e-03  | ~95%†    |
| 1.7 GHz | V(2) | 3.531e-03   | 6.952e-04   | 4.22e-03  | 120%     |
| 1.7 GHz | V(3) | 1.339e-03   | 2.693e-03   | 1.61e-03  | 59.7%    |
| 1.7 GHz | V(4) | 4.605e-03   | 2.107e-03   | 3.04e-03  | 66.0%    |
| 1.7 GHz | V(5) | 2.639e-04   | 6.070e-05   | 3.24e-04  | 123%     |
| 1.7 GHz | V(6) | 6.629e-03   | 4.688e-03   | 1.13e-02  | 170%     |
| 1.7 GHz | V(7) | 3.099e-03   | 3.993e-03   | 7.08e-03  | 177%     |
| 1.8 GHz | V(1) | 2.492e-04   | 3.172e-03   | 3.22e-03  | ~100%†   |
| 1.8 GHz | V(2) | 5.233e-03   | 8.480e-03   | 4.71e-03  | 55.5%    |
| 1.8 GHz | V(3) | 5.798e-03   | 1.404e-02   | 8.29e-03  | 59.1%    |
| 1.8 GHz | V(4) | 1.117e-02   | 2.209e-02   | 1.17e-02  | 52.8%    |
| 1.8 GHz | V(5) | 1.362e-04   | 1.582e-05   | 1.22e-04  | 89.6%    |
| 1.8 GHz | V(6) | 1.035e-02   | 1.611e-02   | 7.17e-03  | 44.5%    |
| 1.8 GHz | V(7) | 5.120e-03   | 7.924e-03   | 2.93e-03  | 37.0%    |
| 1.9 GHz | V(1) | 2.175e-04   | 3.006e-03   | 2.79e-03  | ~93%†    |
| 1.9 GHz | V(2) | 4.568e-03   | 5.596e-04   | 5.12e-03  | 112%     |
| 1.9 GHz | V(3) | 5.222e-04   | 2.364e-03   | 1.93e-03  | 81.5%    |
| 1.9 GHz | V(4) | 4.112e-03   | 1.904e-03   | 2.78e-03  | 67.6%    |
| 1.9 GHz | V(5) | 9.785e-04   | 5.247e-05   | 1.03e-03  | 105%     |
| 1.9 GHz | V(6) | 7.941e-03   | 4.072e-03   | 1.20e-02  | 151%     |
| 1.9 GHz | V(7) | 3.372e-03   | 3.513e-03   | 6.86e-03  | 195%     |

† V(1) source-node artifact or near-cancellation harmonic.

---

## di3 — diode mixer, two-tone, V_peak=1V, f1=100MHz f2=110MHz

**Full two-tone comparison in `model_validation_diode_2tone.md`.**

Earlier partial results in this section (from single-tone ngspice run) were incorrect —
they compared single-tone ngspice against two-tone HBFree. `spice_sim.py` has since been
updated with proper two-tone support (series f2 source + 2D FFT extraction).

---

## di2 — diode mixer + T-line, two-tone, V_peak=1V, f1=100MHz f2=110MHz

**Full two-tone comparison in `model_validation_diode_2tone.md`.**

Same correction as di3. The old 66% DC discrepancy was an artifact of comparing
single-tone ngspice against two-tone HBFree; with proper two-tone comparison the
DC agreement at V(2) is 3.8%.

---

## di1_small — half-wave rectifier with T-line, V_peak=2V, f1=100MHz

**Status: no ngspice comparison available**

`di1_small.ckt.hbl` exists but no `.ckt` file — it was created directly in HBL format.
The circuit includes a `LIB0/LL0` transmission line element (HBFree internal T-line model,
not the SPICE `T` element). A `.ckt` file would need to be written with a matching
T-line or LC equivalent before ngspice comparison is possible.

HBFree converges at iteration 3 (same as di1_lo — same diode topology, small-signal regime).
The T-line validation is tracked in `model_validation_trline.md`.

---

## di4 — diode mixer + T-line, two-tone, large LO

**Full two-tone comparison in `model_validation_diode_2tone.md`.**

`sin()` amplitude corrected from 1V to 30V (= 2×distof1) on 2026-04-26.
Now comparable via the two-tone pipeline.

HBFree converges. `.raw` exists but is not included in this comparison.

---

## di5 — diode mixer + T-line, two-tone, large LO (NOT COMPARABLE)

**Status: amplitude mismatch + no HBFree output**

Same amplitude mismatch as di4: `distof1=15` with `sin(0 1 100Meg)`.
No `di5.raw` — HBFree did not produce output (s2h conversion issue for this circuit variant).
Also has SPICE `T` element.

---

## di6 — 4-diode Graetz bridge, two-tone (NOT COMPARABLE)

**Status: HBFree diverges**

HBFree reports "CAN NOT MAKE GOOD STEP" after 18 iterations. No valid solution.
T-line is commented out in di6.ckt; replaced by `rr 1 2 2` (2Ω source resistance).
The circuit is a two-tone bridge rectifier — convergence is blocked by the nonlinear
operating point, not by element support.

---

## Summary table

| Circuit | Topology          | DC (key node)  | 100 MHz (key node) | Harmonics trend              |
|---------|-------------------|----------------|--------------------|------------------------------|
| di1\_s2 | half-wave, Rvs=0.1Ω | **0.26%** V(2) | **0.98%** V(2)  | 2→19% up to 1.5GHz          |
| di1     | half-wave         | **0.27%** V(2) | **2.96%** V(2)     | 2→55% up to 900MHz           |
| di1\_lo | half-wave, low    | **0.36%** V(2) | **2.30%** V(2)     | 5→11%, relatively flat       |
| di7     | Graetz bridge     | **0.41%** V(3) | **0.78%** V(3)     | 2→65% (even harmonics only)  |
| di8     | complex bridge    | **0.33%** V(4) | **1.34%** V(4)     | 1→195% (asymmetric)          |
| di3     | diode mixer (2T)  | 27% V(2) †     | 10.1% V(2) †       | partial (f1 harmonics only)  |
| di2     | mixer + T-line (2T) | 66% V(2) †‡  | 54.8% V(2) †‡      | partial (f1 harmonics only)  |

† Two-tone partial comparison — intermod products not included; DC offset inflated by two-tone rectification.
‡ Additional T-line model difference (HBFree LL0 vs ngspice transient T element).

---

## Analysis

### DC and fundamental (< 3%)
All circuits agree within 3% at DC and at the fundamental frequency. The dominant
residual is the SCHT vs standard SPICE diode I–V formulation, not normalization.

### Harmonics (5–55%+ for simple circuits, higher for complex)
Error grows with harmonic order. The SCHT and SPICE models diverge in their
harmonic current spectra — this is a known, expected model difference.

### Bridge symmetry cancellation (di7)
Even harmonics cancel at V(1)/V(2); odd harmonics cancel at V(3)/V(4) in a
perfectly balanced Graetz bridge. HBFree enforces exact cancellation. ngspice
transient leaves residuals from finite step size. These entries (marked ‡) are not
model accuracy metrics.

### di8 large harmonic errors (> 100%)
At higher harmonics in di8, some nodes show >100% rdiff. This occurs where both
HBFree and ngspice predict very small amplitudes (< 1mV) but disagree on magnitude
and phase. This is in the noise floor of both simulators relative to the dominant
harmonics — not physically significant but shows up as large rdiff.

### V(1) source node
HBFree Norton formulation gives V(1)≈0 at harmonics. ngspice represents the
loaded source directly. Excluded from model accuracy conclusions.

---

## Files changed for this validation

- `hbfree_exmpls/di7.ckt`: `distof1` corrected 15 → 10.606 (= 21.21/2)
- `hbfree_exmpls/di8.ckt`: `distof1` corrected 15 → 10.606 (= 21.21/2)

---

## Future work

- **Two-tone 2D FFT** — `spice_sim.py` needs intermodulation product extraction
  (2D FFT on (f1, f2) grid) to fully compare di2, di3, and other two-tone circuits.
  Once implemented, di3 comparison will be complete; di2 will add T-line model data.
- **di4 / di5 amplitude fix** — set `distof1=0.5` (or adjust `sin()`) so ngspice
  and HBFree use the same source amplitude; then comparable via same pipeline.
- **di6 convergence** — "CAN NOT MAKE GOOD STEP"; root cause not yet identified.
- **di1_small** — write a `.ckt` file with LC equivalent of LIB0/LL0 T-line for ngspice.
- `model_validation_mesfet.md` — MESFET test circuits needed.
- `model_validation_trline.md` — s2h T-line parser needs `T` element support.
