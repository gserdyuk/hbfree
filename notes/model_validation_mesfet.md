# HBFree MESFET Model Validation — vs ngspice

**Status: PENDING — no MESFET test circuits exist yet**

---

## Planned setup

- HBFree MESFET element: `chanes/libelem/` (identify relevant element file)
- Reference: ngspice MESFET model (`.model` with `level=1` or `level=3`)
- Test circuits needed: small-signal amplifier, mixer

---

## To add test circuits

1. Create `hbfree_exmpls/mesfet1.ckt` with `.hb` line and matching `.model`
2. Ensure `s2h` can parse the MESFET element (check `shpars.y`)
3. Run `scripts/run_test.sh mesfet1.ckt --keep`
4. Run `scripts/spice_sim.py --patch-rs mesfet1.ckt`
5. Compare with `scripts/cmpraw.py mesfet1.raw mesfet1.spice.raw`
6. Record results in this document following the format in `model_validation_diode.md`
