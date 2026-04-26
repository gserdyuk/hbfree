here are test examples.

di1.ckt		1-doide circuit with "t" transmission line under 1-tone excitation
di2.ckt		1-doide circuit with "t" transmission line under 2-tone excitation
di3.ckt		1-doide circuit without transmission line under 2-tone excitation
di4.ckt		as di2.ckt and with extremely large LO signal (15 Volt)
di5.ckt		as di4.ckt, but have a sintax error in input file
di6.ckt   	4-diode bridge rectifier under 2-tone excitation
di7.ckt 	4-diode bridge rectifier under 1-tone excitation
di8.ckt 	too large circuit

test		testing script
./results 	contains reference results


CURRENT TEST STATUS (April, 2004) :
-----------------------------------

TEST        binary release2000          HbFree              
di1.ckt		ok                          OK
di2.ckt		ok                          no convergence (small step)
di3.ckt		ok                          no convergence (small step)
di4.ckt		ok                          no convergence (small step)
di5.ckt		-
di6.ckt   	ok                          no convergence (local failure)
di7.ckt 	ok                          OK
di8.ckt 	-                           converged

NB. That means some problems still persist in open code -
most probably - boundary violations (remember - fortran with static arrays)
Gennady Serdyuk


UPDATE (April, 2026) :
----------------------

Circuit descriptions corrected:
  di1.ckt   1-diode half-wave rectifier, 1-tone, no transmission line
  di2.ckt   1-diode mixer + lossless T-line, 2-tone excitation
  di3.ckt   1-diode mixer, 2-tone excitation, no T-line
  di4.ckt   as di2.ckt, large LO signal (distof1=15, V_peak=30V)
  di5.ckt   as di4.ckt, different harmonic set (di5 has syntax note below)
  di6.ckt   4-diode bridge rectifier, 2-tone excitation
  di7.ckt   Graetz 4-diode bridge rectifier, 1-tone excitation
  di8.ckt   complex 4-diode bridge, 1-tone excitation
  di1_lo    di1 variant with low drive level
  di1_s2    di1 variant with small source resistance (Rvs=0.1 Ohm)
  di1_small di1 with HBFree-internal lossless T-line (no SPICE .ckt)

Normalization fix applied (April 2026):
  sq12 = 2.0 in jdrive and emf (was sqrt(2.0)).
  All .raw golden files in results/ updated after this fix.
  Convention: .raw stores V_rms for AC harmonics, actual value for DC.

CURRENT TEST STATUS (April, 2026) :

TEST        HBFree                      Notes
di1.ckt     CONVERGED (13 iter)         validated vs ngspice, DC err 0.27%
di1_lo      CONVERGED (3 iter)          validated vs ngspice
di1_s2      CONVERGED (12 iter)         validated vs ngspice, DC err 0.26%
di2.ckt     CONVERGED (step<1e-8, 34)   validated vs ngspice 2-tone, DC err 3.8%
di3.ckt     CONVERGED (step<1e-8, 25)   validated vs ngspice 2-tone, DC err 10.9%
di4.ckt     CONVERGED (step<1e-8, 38)   validated vs ngspice 2-tone, DC err 7.4%
di5.ckt     S2H FAIL                    s2h cannot parse lossless T element
di6.ckt     DIVERGED (18 iter)          "CAN NOT MAKE GOOD STEP" - open issue
di7.ckt     CONVERGED (15 iter)         validated vs ngspice
di8.ckt     CONVERGED (18 iter)         validated vs ngspice
di1_small   CONVERGED (3 iter)          uses HBFree LL0 T-line internally

Golden files in results/:
  *.raw          phasor results - primary regression golden
  *.conv         convergence summary (status, iterations, final error/norm)
  *.ckt.run.log  HBFree solver log - reference only, not used for regression
  *.ckt.log      s2h conversion log - reference only

Scripts (../scripts/):
  spice_sim.py        generate ngspice reference .spice.raw for validation
  cmpraw.py           compare two .raw files with atol/rtol tolerances
  make_conv_summary.py extract .conv golden from .ckt.run.log
  run_test.sh         full s2h -> hbl -> cmpraw pipeline wrapper

See ../notes/ for detailed validation results and analysis.

