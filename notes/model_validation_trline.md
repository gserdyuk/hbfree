# HBFree Transmission Line Model Validation — vs ngspice

**Status: PENDING — s2h cannot parse lossless T element; test circuits not validated**

---

## Known blocker

`spice2hbl/s2h` fails on the SPICE lossless transmission line element `T`:

```
#ER1-5543: Error in Lossless Transmission Line 't1':
t1  2 0 out 0 0 z0=50 f=100MEG nl=0.5
```

di2 and di4 use this element and cannot be converted by `s2h` until this is fixed.

---

## Planned approach (after s2h fix)

- Test circuits: di2, di4 (already in `hbfree_exmpls/`, need s2h support)
- di2: diode + T-line, single-tone, two-frequency HB
- di4: diode + T-line, two-tone, high drive (distof1=15, needs distof1 correction to match sin() amplitude)
- Comparison: `scripts/spice_sim.py` + `scripts/cmpraw.py`

---

## To add T-line support to s2h

1. Locate T-line parsing in `spice2hbl/shpars.y` or `shpass1.c`
2. Map to the corresponding HBFree T-line element in `chanes/libelem/`
3. Verify with di2: `s2h -c../spice2hbl/s2h.cfg di2.ckt di2.ckt.hbl`
4. Run and record results following `model_validation_diode.md` format
