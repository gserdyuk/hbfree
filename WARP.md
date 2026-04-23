# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

Project overview
- Purpose: Harmonic Balance (HB) simulator for nonlinear RF circuits.
- Flow (big picture):
  1) A SPICE-like netlist (.ckt) is translated to an intermediate “channel” file (.chan) by the C/flex/bison translator s2h found in spice2hbl/.
  2) The Fortran HB solver hbl in chanes/ consumes the .chan and performs the harmonic balance solve, emitting logs/results.
  3) scripts/hbfree is a convenience wrapper that runs s2h then hbl for a given .ckt file and tees hbl output to a .hbl file.
- Key components:
  - chanes/: legacy Fortran 77 HB core (Newton/HB loop, linear solver, element/device library under chanes/libelem/). Builds top-level executable hbl.
  - spice2hbl/: SPICE→HBL translator (C + flex/bison); builds top-level executable s2h and uses s2h.cfg for parsing rules.
  - scripts/hbfree: Bash wrapper for end-to-end run (expects s2h and hbl in PATH; default config path /etc/s2h.cfg).
  - hbfree_exmpls/: Example circuits and simple test scripts.

Prerequisites
- Build tools: make, gcc/g++, gfortran, ar, ranlib, ld
- Generators: flex, bison
- Optional: tee (for scripts), coreutils

Build commands
- Build translator (s2h):
  ```bash
  make -C spice2hbl
  # Output: ./s2h (in repo root)
  ```
- Build solver (hbl):
  ```bash
  make -C chanes
  # Output: ./hbl (in repo root)
  ```
- Clean object files:
  ```bash
  make -C chanes clean && make -C spice2hbl clean
  ```
- Windows wrapper (optional, not required on Linux):
  ```bash
  make -C hbw
  # Output: hbw/hbw.exe
  ```

Configuration (s2h)
- The wrapper script scripts/hbfree uses CONFIG=/etc/s2h.cfg.
- Options to use a local config without system install:
  - Edit scripts/hbfree and set CONFIG to "$REPO_ROOT/spice2hbl/s2h.cfg"; or
  - Manually run s2h with -c path as shown below (bypasses the wrapper).

Run commands
- End-to-end using the wrapper (requires s2h and hbl in PATH):
  ```bash
  # Add repo root (for ./s2h, ./hbl) and scripts/ (for hbfree) to PATH for this shell
  PATH="$PWD:$PWD/scripts:$PATH" \
  hbfree hbfree_exmpls/di1.ckt
  # Produces: di1.chan, di1.elems, di1.nodes, di1.hbl in hbfree_exmpls/
  ```
- End-to-end manually (explicit config path, useful if not modifying the wrapper):
  ```bash
  cd hbfree_exmpls
  ../s2h -c ../spice2hbl/s2h.cfg -m di1.elems -n di1.nodes di1.ckt di1.chan
  ../hbl di1.chan | tee di1.hbl
  ```

Examples and quick tests
- Run a single example:
  ```bash
  PATH="$PWD:$PWD/scripts:$PATH" \
  scripts/hbfree hbfree_exmpls/di1.ckt
  ```
- Run all bundled examples (note: the script references a non-tracked "compare" helper; runs will complete but final compare may be skipped):
  ```bash
  cd hbfree_exmpls
  PATH="$PWD/..:$PWD/../scripts:$PATH" ./test
  # Clean example outputs:
  ./test-clean
  ```

Notes on the Fortran build (chanes/)
- Compiler flags default most REAL to 8 bytes (-fdefault-real-8 etc.) and use legacy Fortran mode; commented options include bounds checking. If debugging array bounds or convergence issues, consider enabling gfortran options like -fbounds-check locally for specific objects.

Repository structure (high-level)
- chanes/: solver main and numerics; libelem/ contains element/device models used during stamping/linearization.
- spice2hbl/: parser/lexer (shpars.y, shlpars.l) and passes; obj/ contains generated/compiled artifacts; s2h.cfg defines translation rules.
- scripts/: operational wrapper (hbfree) orchestrating s2h → hbl and log capture.
- hbfree_exmpls/: sample circuits (.ckt) and reference outputs under results/.
- hbfree_doc/: docs (HTML/PDF) and INSTALL notes (legacy: refers to an old path; there is no top-level install in this repo).

Troubleshooting
- If hbfree_exmpls/test fails at the final "./compare" step, it’s because the helper is not committed; the simulations themselves will have run.
- Convergence behavior varies by example (see hbfree_exmpls/exmpls_readme.txt for historical notes). Small steps or boundary conditions in legacy Fortran arrays can affect results.
