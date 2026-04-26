# HBFree Diode Model Validation — Two-Tone Circuits

**Date: 2026-04-26**
**Companion to**: `model_validation_diode.md` — setup, models, RS mismatch, convention
and all single-tone results are documented there. This file covers only two-tone circuits.

---

## Two-tone pipeline (new)

`spice_sim.py` was updated on 2026-04-26 with proper two-tone support:

1. **`parse_hb`**: fixed regex `(-?\d+)\s+(-?\d+)` to handle negative harmonic indices
   (e.g. `1 -3`); previously silently produced wrong pairs.
2. **`build_spice`**: for circuits with `distof2`, injects a series voltage source
   `v_f2_{name}` at `f2` with amplitude `2*distof2` (HBFree convention). Previously
   ngspice ran single-tone only, making di2/di3 comparisons in the main doc incorrect.
3. **`fft_phasors_2tone`**: 1D FFT with window `T=1/(f2-f1)` maps each `n1*f1+n2*f2`
   product to a unique FFT bin. Negative-frequency pairs (e.g. IF at `-10 MHz`) use
   conjugate of the positive bin.
4. **`write_hbfree_raw_2tone`**: writes signed `n1*f1+n2*f2` frequencies to match
   HBFree's `.raw` format so `cmpraw.py` matches all pairs by frequency.

---

## Circuit inventory

| Circuit | Description                  | HBFree    | T-line  | V\_peak(LO) | V\_peak(IF) | Comparison |
|---------|------------------------------|-----------|---------|-------------|-------------|------------|
| di3     | diode mixer, no T-line       | CONVERGED | no      | 1V          | 1V          | **FULL**   |
| di2     | diode mixer + T-line         | CONVERGED | SPICE T | 1V          | 1V          | **FULL**   |
| di4     | diode mixer + T-line, lg. LO | CONVERGED | SPICE T | 30V         | 0.4V        | **FULL**   |
| di5     | same as di4, diff harmonics  | no .raw   | SPICE T | 30V         | 0.4V        | NOT COMPARABLE |

di5: s2h conversion failed; no HBFree `.raw` exists. `di5.ckt` amplitude was also
corrected (`sin(0 30 100Meg)`) on 2026-04-26 for future use.

---

## V(1) source node

See main doc note. At f1 and f2, V(1) reflects the loaded source voltage
(ngspice) vs. ideal Norton (HBFree≈0); excluded from model accuracy conclusions.

---

## di3 — diode + Rl=1kΩ, two-tone, V\_peak=1V each, f1=100MHz f2=110MHz

**Nodes**: V(1)=source, V(2)=diode–load, V(out)=load (V(2)≈V(out) here; Rt=1Ω between them)

All 20 harmonic pairs × 3 nodes = 60 comparisons. V(1) rows retained for completeness (†).

| Freq              | Product        | Node   | Ref \|z\|   | Test \|z\|  | rdiff   |
|-------------------|----------------|--------|-------------|-------------|---------|
| DC                | (0,0)          | V(1)   | 2.335e-05   | 2.656e-05   | 12.1%†  |
| DC                | (0,0)          | V(2)   | 2.3368e-01  | 2.6239e-01  | **10.9%** |
| DC                | (0,0)          | V(out) | 2.3345e-01  | 2.6213e-01  | **10.9%** |
| 110 MHz           | f2 (0,1)       | V(1)   | 7.0708e-01  | 7.0123e-01  | 0.8%†   |
| 110 MHz           | f2 (0,1)       | V(2)   | 2.4123e-01  | 2.5285e-01  | **4.6%**  |
| 110 MHz           | f2 (0,1)       | V(out) | 2.4099e-01  | 2.5259e-01  | **4.6%**  |
| 220 MHz           | 2f2 (0,2)      | V(1)   | 8.227e-06   | 1.263e-03   | ~99%†   |
| 220 MHz           | 2f2 (0,2)      | V(2)   | 8.2354e-02  | 6.3509e-02  | 22.9%   |
| 220 MHz           | 2f2 (0,2)      | V(out) | 8.2272e-02  | 6.3445e-02  | 22.9%   |
| 330 MHz           | 3f2 (0,3)      | V(1)   | 9.461e-07   | 7.223e-04   | ~99%†   |
| 330 MHz           | 3f2 (0,3)      | V(2)   | 9.4703e-03  | 4.3951e-03  | 53.6%   |
| 330 MHz           | 3f2 (0,3)      | V(out) | 9.4608e-03  | 4.3907e-03  | 53.6%   |
| 100 MHz           | f1 (1,0)       | V(1)   | 7.0708e-01  | 7.1193e-01  | 0.7%†   |
| 100 MHz           | f1 (1,0)       | V(2)   | 2.2825e-01  | 2.5655e-01  | **11.0%** |
| 100 MHz           | f1 (1,0)       | V(out) | 2.2802e-01  | 2.5630e-01  | **11.0%** |
| 210 MHz           | f1+f2 (1,1)    | V(1)   | 1.518e-05   | 1.362e-03   | ~98%†   |
| 210 MHz           | f1+f2 (1,1)    | V(2)   | 1.5190e-01  | 1.6431e-01  | 7.5%    |
| 210 MHz           | f1+f2 (1,1)    | V(out) | 1.5175e-01  | 1.6415e-01  | 7.5%    |
| 320 MHz           | f1+2f2 (1,2)   | V(1)   | 2.201e-06   | 7.503e-04   | ~99%†   |
| 320 MHz           | f1+2f2 (1,2)   | V(2)   | 2.2029e-02  | 2.5400e-02  | 13.3%   |
| 320 MHz           | f1+2f2 (1,2)   | V(out) | 2.2007e-02  | 2.5375e-02  | 13.3%   |
| -230 MHz          | f1-3f2 (1,-3)  | V(1)   | 1.566e-06   | 1.178e-03   | ~99%†   |
| -230 MHz          | f1-3f2 (1,-3)  | V(2)   | 1.5675e-02  | 1.4937e-02  | 4.7%    |
| -230 MHz          | f1-3f2 (1,-3)  | V(out) | 1.5660e-02  | 1.4922e-02  | 4.7%    |
| -120 MHz          | f1-2f2 (1,-2)  | V(1)   | 2.981e-06   | 8.661e-03   | ~99%†   |
| -120 MHz          | f1-2f2 (1,-2)  | V(2)   | 2.9843e-02  | 2.2519e-02  | 24.5%   |
| -120 MHz          | f1-2f2 (1,-2)  | V(out) | 2.9813e-02  | 2.2496e-02  | 24.5%   |
| -10 MHz           | IF f2-f1 (1,-1)| V(1)   | 1.528e-05   | 1.994e-04   | 92.3%†  |
| -10 MHz           | IF f2-f1 (1,-1)| V(2)   | 1.5296e-01  | 1.6395e-01  | **6.7%**  |
| -10 MHz           | IF f2-f1 (1,-1)| V(out) | 1.5281e-01  | 1.6378e-01  | **6.7%**  |
| 200 MHz           | 2f1 (2,0)      | V(1)   | 6.347e-06   | 1.482e-03   | ~99%†   |
| 200 MHz           | 2f1 (2,0)      | V(2)   | 6.3532e-02  | 6.9734e-02  | 8.9%    |
| 200 MHz           | 2f1 (2,0)      | V(out) | 6.3469e-02  | 6.9665e-02  | 8.9%    |
| 310 MHz           | 2f1+f2 (2,1)   | V(1)   | 1.709e-06   | 7.806e-04   | ~99%†   |
| 310 MHz           | 2f1+f2 (2,1)   | V(2)   | 1.7093e-02  | 2.6673e-02  | 35.9%   |
| 310 MHz           | 2f1+f2 (2,1)   | V(out) | 1.7076e-02  | 2.6647e-02  | 35.9%   |
| -20 MHz           | 2f1-2f2 (2,-2) | V(1)   | — *         | —           | —       |
| -20 MHz           | 2f1-2f2 (2,-2) | V(2)   | — *         | —           | —       |
| 90 MHz            | 2f1-f2 (2,-1)  | V(1)   | 1.562e-06   | 7.299e-03   | ~99%†   |
| 90 MHz            | 2f1-f2 (2,-1)  | V(2)   | 1.5637e-02  | 2.8562e-02  | 45.3%   |
| 90 MHz            | 2f1-f2 (2,-1)  | V(out) | 1.5622e-02  | 2.8533e-02  | 45.3%   |
| 300 MHz           | 3f1 (3,0)      | V(1)   | 4.626e-07   | 8.138e-04   | ~99%†   |
| 300 MHz           | 3f1 (3,0)      | V(2)   | 4.6304e-03  | 2.8082e-03  | 39.4%   |
| 300 MHz           | 3f1 (3,0)      | V(out) | 4.6258e-03  | 2.8054e-03  | 39.4%   |
| 410 MHz           | 3f1+f2 (3,1)   | V(1)   | 1.568e-06   | 6.929e-04   | ~99%†   |
| 410 MHz           | 3f1+f2 (3,1)   | V(2)   | 1.8009e-02  | 1.3273e-02  | 26.3%   |
| 410 MHz           | 3f1+f2 (3,1)   | V(out) | 1.7991e-02  | 1.3260e-02  | 26.3%   |
| 190 MHz           | 3f1-f2 (3,-1)  | V(1)   | 1.564e-06   | 7.077e-04   | ~99%†   |
| 190 MHz           | 3f1-f2 (3,-1)  | V(2)   | 2.0762e-02  | 1.0065e-02  | 51.5%   |
| 190 MHz           | 3f1-f2 (3,-1)  | V(out) | 2.0742e-02  | 1.0055e-02  | 51.5%   |
| 400 MHz           | 4f1 (4,0)      | V(1)   | 8.926e-07   | 5.755e-04   | ~99%†   |
| 400 MHz           | 4f1 (4,0)      | V(2)   | 8.9353e-03  | 2.3966e-04  | 97.3%   |
| 400 MHz           | 4f1 (4,0)      | V(out) | 8.9263e-03  | 2.3942e-04  | 97.3%   |
| 510 MHz           | 4f1+f2 (4,1)   | V(1)   | 6.793e-07   | 4.604e-04   | ~99%†   |
| 510 MHz           | 4f1+f2 (4,1)   | V(2)   | 7.2689e-04  | 4.5462e-04  | 37.5%   |
| 510 MHz           | 4f1+f2 (4,1)   | V(out) | 7.2616e-04  | 4.5417e-04  | 37.5%   |
| 290 MHz           | 4f1-f2 (4,-1)  | V(1)   | 4.960e-09   | 8.461e-04   | ~99%†   |
| 290 MHz           | 4f1-f2 (4,-1)  | V(2)   | 8.4557e-04  | 1.1212e-03  | 24.6%   |
| 290 MHz           | 4f1-f2 (4,-1)  | V(out) | 8.4473e-04  | 1.1201e-03  | 24.6%   |
| 500 MHz           | 5f1 (5,0)      | V(1)   | 9.393e-09   | 4.432e-04   | ~99%†   |
| 500 MHz           | 5f1 (5,0)      | V(2)   | 9.4102e-05  | 3.3747e-04  | 72.1%   |
| 500 MHz           | 5f1 (5,0)      | V(out) | 9.4008e-05  | 3.3713e-04  | 72.1%   |

† V(1) source-node artifact. * (2,-2) = -20 MHz not in di3's harmonic list.

---

## di2 — diode + T-line + Rl=1kΩ, two-tone, V\_peak=1V each, f1=100MHz f2=110MHz

**Nodes**: V(1)=source, V(2)=diode–T-line input, V(out)=T-line output

**T-line**: z0=50Ω, f=100MHz, nl=0.5 (quarter-wave at f1). HBFree uses frequency-domain
LL0; ngspice uses time-domain T element. This adds a second model difference on top
of SCHT vs SPICE diode. `di2` uses the same smaller harmonic set as di4/di5.

Note: di2 harmonic list has (2,-2) = -20 MHz not present in di3.

| Freq              | Product        | Node   | Ref \|z\|   | Test \|z\|  | rdiff   |
|-------------------|----------------|--------|-------------|-------------|---------|
| DC                | (0,0)          | V(2)   | 5.0850e-01  | 5.2870e-01  | **3.8%**  |
| DC                | (0,0)          | V(out) | 5.0825e-01  | 5.2880e-01  | 3.9%    |
| 110 MHz           | f2 (0,1)       | V(2)   | 1.0751e-01  | 1.0602e-01  | **1.4%**  |
| 110 MHz           | f2 (0,1)       | V(out) | 1.1294e-01  | 1.1144e-01  | 1.3%    |
| 100 MHz           | f1 (1,0)       | V(2)   | 4.5433e-01  | 4.5703e-01  | **0.6%**  |
| 100 MHz           | f1 (1,0)       | V(out) | 4.5456e-01  | 4.5709e-01  | 0.6%    |
| -10 MHz           | IF (1,-1)      | V(2)   | 6.7073e-02  | 6.5925e-02  | **1.7%**  |
| -10 MHz           | IF (1,-1)      | V(out) | 7.0511e-02  | 6.9351e-02  | 1.6%    |
| 210 MHz           | f1+f2 (1,1)    | V(2)   | 7.1063e-02  | 6.6237e-02  | 6.8%    |
| 210 MHz           | f1+f2 (1,1)    | V(out) | 7.4600e-02  | 6.9339e-02  | 7.0%    |
| 200 MHz           | 2f1 (2,0)      | V(2)   | 6.6965e-02  | 4.2218e-02  | 37.0%   |
| 200 MHz           | 2f1 (2,0)      | V(out) | 6.6933e-02  | 4.2287e-02  | 36.8%   |
| 90 MHz            | 2f1-f2 (2,-1)  | V(2)   | 1.1686e-02  | 1.4365e-02  | 18.7%   |
| 90 MHz            | 2f1-f2 (2,-1)  | V(out) | 1.2293e-02  | 1.5029e-02  | 18.2%   |
| 310 MHz           | 2f1+f2 (2,1)   | V(2)   | 1.3682e-02  | 1.3185e-02  | 3.6%    |
| 310 MHz           | 2f1+f2 (2,1)   | V(out) | 1.4353e-02  | 1.4036e-02  | 2.2%    |
| -20 MHz           | 2IF (2,-2)     | V(2)   | 7.5524e-03  | 9.3469e-03  | 19.2%   |
| -20 MHz           | 2IF (2,-2)     | V(out) | 9.3262e-03  | 1.1784e-02  | 20.9%   |
| -120 MHz          | f1-2f2 (1,-2)  | — †    | — (not in di2 harmonic list) | | |
| -230 MHz          | f1-3f2 (1,-3)  | — †    | — (not in di2 harmonic list) | | |
| 220 MHz           | 2f2 (0,2)      | V(2)   | 4.5181e-02  | 4.0073e-02  | 11.3%   |
| 220 MHz           | 2f2 (0,2)      | V(out) | 5.5618e-02  | 4.9499e-02  | 11.0%   |
| 330 MHz           | 3f2 (0,3)      | — †    | — (not in di2 harmonic list) | | |
| 300 MHz           | 3f1 (3,0)      | V(2)   | 3.2425e-02  | 3.1252e-02  | 3.6%    |
| 300 MHz           | 3f1 (3,0)      | V(out) | 3.2442e-02  | 3.1151e-02  | 4.0%    |
| 190 MHz           | 3f1-f2 (3,-1)  | V(2)   | 6.0470e-03  | 6.6565e-03  | 9.2%    |
| 190 MHz           | 3f1-f2 (3,-1)  | V(out) | 6.3659e-03  | 6.3236e-03  | 0.7%    |
| 410 MHz           | 3f1+f2 (3,1)   | V(2)   | 5.9461e-03  | 6.4927e-03  | 8.4%    |
| 410 MHz           | 3f1+f2 (3,1)   | V(out) | 6.2335e-03  | 6.8154e-03  | 8.5%    |
| 400 MHz           | 4f1 (4,0)      | V(2)   | 9.5083e-03  | 1.8367e-02  | 48.2%   |
| 400 MHz           | 4f1 (4,0)      | V(out) | 9.5039e-03  | 1.8286e-02  | 48.0%   |
| 290 MHz           | 4f1-f2 (4,-1)  | V(2)   | 1.5083e-03  | 3.4443e-03  | 56.2%   |
| 290 MHz           | 4f1-f2 (4,-1)  | V(out) | 1.5890e-03  | 3.9220e-03  | 59.5%   |
| 510 MHz           | 4f1+f2 (4,1)   | V(2)   | 1.4385e-03  | 3.6426e-03  | 60.5%   |
| 510 MHz           | 4f1+f2 (4,1)   | V(out) | 1.5070e-03  | 3.7134e-03  | 59.4%   |

† These products are in di3's harmonic list but not di2's — di2 uses fewer pairs.
V(1) omitted from di2 table for brevity; follows same source-node artifact pattern.

---

## di4 — diode + T-line + Rl=1kΩ, large LO, V\_peak(LO)=30V, V\_peak(IF)=0.4V, f1=100MHz f2=110MHz

**Nodes**: V(1)=source, V(2)=diode–T-line input, V(out)=T-line output

**Same T-line as di2** (quarter-wave at f1). Same harmonic list as di5.

**Large-signal note**: At 30V LO, the diode conducts extremely strongly over half the
cycle. Small differences in the SCHT vs SPICE nonlinear model produce large relative
errors in the mixing products, even when LO harmonics themselves agree at 7–10%.

| Freq    | Product        | Node   | Ref \|z\|   | Test \|z\|  | rdiff   |
|---------|----------------|--------|-------------|-------------|---------|
| DC      | (0,0)          | V(2)   | 8.6379e+00  | 9.3309e+00  | **7.4%**  |
| DC      | (0,0)          | V(out) | 8.6335e+00  | 9.3308e+00  | 7.5%    |
| 110 MHz | f2 (0,1)       | V(2)   | 1.4933e+00  | 1.2913e-01  | 91.4%   |
| 110 MHz | f2 (0,1)       | V(out) | 1.5687e+00  | 1.3603e-01  | 91.3%   |
| 100 MHz | f1 (1,0)       | V(2)   | 9.4752e+00  | 1.0398e+01  | **8.9%**  |
| 100 MHz | f1 (1,0)       | V(out) | 9.4800e+00  | 1.0398e+01  | 8.8%    |
| -10 MHz | IF (1,-1)      | V(2)   | 1.1046e+00  | 8.8636e-02  | 92.0%   |
| -10 MHz | IF (1,-1)      | V(out) | 1.1612e+00  | 9.3206e-02  | 92.0%   |
| 210 MHz | f1+f2 (1,1)    | V(2)   | 1.1185e+00  | 8.3387e-02  | 92.5%   |
| 210 MHz | f1+f2 (1,1)    | V(out) | 1.1741e+00  | 8.7964e-02  | 92.5%   |
| 200 MHz | 2f1 (2,0)      | V(2)   | 4.0909e+00  | 4.4787e+00  | **8.7%**  |
| 200 MHz | 2f1 (2,0)      | V(out) | 4.0889e+00  | 4.4788e+00  | 8.7%    |
| 90 MHz  | 2f1-f2 (2,-1)  | V(2)   | 3.3913e-01  | 1.1814e-02  | 96.5%   |
| 90 MHz  | 2f1-f2 (2,-1)  | V(out) | 3.5676e-01  | 1.1581e-02  | 96.8%   |
| 310 MHz | 2f1+f2 (2,1)   | V(2)   | 3.4248e-01  | 4.8084e-03  | 98.6%   |
| 310 MHz | 2f1+f2 (2,1)   | V(out) | 3.5927e-01  | 4.9541e-03  | 98.6%   |
| 300 MHz | 3f1 (3,0)      | V(2)   | 1.7318e-01  | 6.0158e-02  | 65.3%   |
| 300 MHz | 3f1 (3,0)      | V(out) | 1.7327e-01  | 6.0158e-02  | 65.3%   |
| 190 MHz | 3f1-f2 (3,-1)  | V(2)   | 2.9161e-01  | 2.2141e-02  | 92.4%   |
| 190 MHz | 3f1-f2 (3,-1)  | V(out) | 3.0699e-01  | 2.3449e-02  | 92.4%   |
| 410 MHz | 3f1+f2 (3,1)   | V(2)   | 3.0463e-01  | 2.6905e-02  | 91.2%   |
| 410 MHz | 3f1+f2 (3,1)   | V(out) | 3.1935e-01  | 2.8417e-02  | 91.1%   |
| 400 MHz | 4f1 (4,0)      | V(2)   | 9.7166e-01  | 8.8585e-01  | **8.8%**  |
| 400 MHz | 4f1 (4,0)      | V(out) | 9.7121e-01  | 8.8589e-01  | 8.8%    |
| 290 MHz | 4f1-f2 (4,-1)  | V(2)   | 2.7776e-01  | 3.3412e-03  | 98.8%   |
| 290 MHz | 4f1-f2 (4,-1)  | V(out) | 2.9262e-01  | 3.6768e-03  | 98.7%   |
| 510 MHz | 4f1+f2 (4,1)   | V(2)   | 2.9672e-01  | 4.3036e-03  | 98.5%   |
| 510 MHz | 4f1+f2 (4,1)   | V(out) | 3.1085e-01  | 4.4662e-03  | 98.6%   |
| 500 MHz | 5f1 (5,0)      | V(2)   | 9.4311e-02  | 3.4068e-02  | 63.9%   |
| 500 MHz | 5f1 (5,0)      | V(out) | 9.4316e-02  | 3.4050e-02  | 63.9%   |
| 390 MHz | 5f1-f2 (5,-1)  | V(2)   | 3.7783e-02  | 1.4613e-02  | 61.3%   |
| 390 MHz | 5f1-f2 (5,-1)  | V(out) | 3.9833e-02  | 1.5387e-02  | 61.4%   |
| 610 MHz | 5f1+f2 (5,1)   | V(2)   | 3.4910e-02  | 1.5680e-02  | 55.1%   |
| 610 MHz | 5f1+f2 (5,1)   | V(out) | 3.6547e-02  | 1.6571e-02  | 54.7%   |
| 600 MHz | 6f1 (6,0)      | V(2)   | 3.3536e-01  | 3.7378e-01  | **10.3%** |
| 600 MHz | 6f1 (6,0)      | V(out) | 3.3539e-01  | 3.7384e-01  | 10.3%   |
| 700 MHz | 7f1 (7,0)      | V(2)   | 1.7224e-01  | 2.3197e-02  | 86.5%   |
| 700 MHz | 7f1 (7,0)      | V(out) | 1.7226e-01  | 2.3190e-02  | 86.5%   |

V(1) omitted; V(1) at f1=100MHz PASS (7.07e-01 ref vs 7.07e-01 test, 0.01%).

---

## Analysis

### Low-order products: SCHT diverges in two-tone mode

For di3 (small signal, no T-line):
- DC: 10.9% — much larger than single-tone di1 (0.27%). Two tones drive more rectification;
  the SCHT model redistributes charge differently from the standard SPICE diode.
- f1, f2 fundamentals: ~5–11%
- IF (10 MHz): **6.7%** — the most practically important result. Acceptable given model mismatch.
- f1+f2 (210 MHz): 7.5%

For di2 (same small signal, with T-line):
- The quarter-wave T-line at f1=100MHz creates an open at f1 and short at 2f1.
- DC: 3.8% — *better* than di3 despite the T-line. Suggests the T-line's impedance
  transformation puts the diode in a regime where SCHT/SPICE models agree more closely.
- f1, f2, IF: 0.6–1.7% — remarkably good agreement for small-signal two-tone.
- 2f1 (200 MHz = half-wave T-line): 37% — the T-line model difference is largest at
  frequencies where LL0 (frequency-domain) and time-domain T disagree most.

### LO harmonic pattern in di4

For di4 (large LO 30V):
- LO harmonics n·f1 for even n (200, 400, 600 MHz) agree at **7–11%**.
- LO harmonics n·f1 for odd n≥3 (300, 500, 700 MHz) show larger errors (55–87%).
- All mixing products involving f2 (n2≠0) show **91–99%** errors.

This two-tier pattern means: HBFree's SCHT+LL0 model reproduces the LO self-mixing
(rectification at even harmonics) reasonably well even at 30V, but the cross-product
generation (diode converting LO harmonics to f2-sideband frequencies) diverges strongly.
The IF and f2 sideband predictions from HBFree are not reliable at large LO for T-line
circuits — use the SPICE results for circuit design in this regime.

### Negative-frequency products

HBFree stores `n1*f1 + n2*f2` as signed. Negative-frequency products (e.g. IF at
−10 MHz, f1−2f2 at −120 MHz) are real intermodulation products stored with their
natural signed frequency. `cmpraw.py` matches them correctly since both .raw files
use the same sign convention after the `write_hbfree_raw_2tone` update.

---

## Files changed for this validation

- `scripts/spice_sim.py`: two-tone pipeline (parse_hb fix, build_spice f2 injection,
  fft_phasors_2tone, write_hbfree_raw_2tone)
- `hbfree_exmpls/di4.ckt`: `sin(0 1 ...)` corrected to `sin(0 30 ...)` (distof1=15 match)
- `hbfree_exmpls/di5.ckt`: same amplitude correction (di5 not yet runnable via HBFree)
