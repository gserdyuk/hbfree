#!/usr/bin/env python3
"""
compare_iv.py  —  Compare diode I-V: analytical model, ngspice DC sweep, HBFree DC point.

Usage:
  python3 scripts/compare_iv.py [--vdc FLOAT] [--out FILE] [--hbl-dir DIR]
  python3 scripts/compare_iv.py --help

Steps performed:
  1. Analytical I-V: solve I = IS*(exp((V - RS*I)/Vt) - 1) by Newton iteration.
     Both HBFree SCHT and SPICE standard diode reduce to this formula with the
     same IS, N, RS, T — so a single curve represents both models.
  2. ngspice DC sweep: write a minimal SPICE deck, run ngspice, overlay result.
  3. HBFree DC verification: write a minimal .hbl (KN=1, DC source), run hbl,
     extract V(1)/V(2) at DC, compute I and overlay on the same plot.
  4. The HBFree and ngspice points should land on the analytical curve if the
     models are equivalent.  Any offset reveals the discrepancy.
"""

import sys
import os
import re
import argparse
import subprocess
import tempfile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ---------------------------------------------------------------------------
# Diode parameters (matching di1_s2 circuit)
# ---------------------------------------------------------------------------
IS   = 1e-9     # A  — saturation current
N    = 1.0      # emission coefficient
T    = 300.0    # K  — temperature
RS   = 1.0      # Ω  — series resistance (Rsmin from s2h.cfg)
RL   = 1000.0   # Ω  — load resistor in DC test circuit
RVS  = 0.001    # Ω  — source series resistance (nearly ideal)

VT   = N * 1.380649e-23 * T / 1.602176634e-19   # N·kT/q  ≈ 0.02585 V at 300 K


# ---------------------------------------------------------------------------
# Analytical I-V  (Newton iteration for implicit I = IS*(exp((V-RS*I)/Vt)-1))
# ---------------------------------------------------------------------------

def diode_current(V_term, max_iter=200, tol=1e-15):
    """Return diode current I for terminal voltage V_term.
    Solves: f(I) = I - IS*(exp((V_term - RS*I)/VT) - 1) = 0
    """
    if V_term <= 0:
        return IS * (np.exp(V_term / VT) - 1.0)   # RS barely matters in reverse

    I = IS * np.exp(V_term / VT / 2)   # starting guess

    for _ in range(max_iter):
        Vj  = V_term - RS * I
        # clamp Vj to avoid overflow
        Vj  = min(Vj, 40.0 * VT)
        exp_val = np.exp(Vj / VT)
        f   = I - IS * (exp_val - 1.0)
        df  = 1.0 + IS * (RS / VT) * exp_val        # df/dI
        step = f / df
        I  -= step
        if abs(step) < tol * (1.0 + abs(I)):
            break

    return I


def analytical_iv(v_range):
    return np.array([diode_current(v) for v in v_range])


# ---------------------------------------------------------------------------
# ngspice DC sweep
# ---------------------------------------------------------------------------

SPICE_DECK = """\
diode IV sweep
* IS={is_} N={n} RS={rs}
Vd  1 0  dc 0
D1  1 2  dmod
Rsense 2 0  1e-3
.model dmod d is={is_} n={n} rs={rs}
.options filetype=ascii
.dc Vd {vmin} {vmax} {vstep}
.save i(Vd)  v(1)
.end
"""

def run_ngspice_dc(vmin=-0.7, vmax=0.65, vstep=0.002, ngspice='ngspice'):
    deck = SPICE_DECK.format(is_=IS, n=N, rs=RS,
                             vmin=vmin, vmax=vmax, vstep=vstep)
    with tempfile.TemporaryDirectory() as td:
        spice_file = os.path.join(td, 'diode_iv.sp')
        raw_file   = os.path.join(td, 'diode_iv.raw')
        with open(spice_file, 'w') as f:
            f.write(deck)
        cmd = [ngspice, '-b', '-r', raw_file, spice_file]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print('ngspice stderr:', r.stderr[:500], file=sys.stderr)
            return None, None
        return _parse_dc_raw(raw_file)


def _parse_dc_raw(rawfile):
    """Parse ngspice DC sweep ASCII raw: return (v_arr, i_arr) numpy arrays."""
    with open(rawfile, encoding='latin-1') as f:
        lines = f.readlines()

    variables = []
    section = 'header'
    data = {}

    for line in lines:
        s = line.strip()
        if section == 'header':
            if s.lower().startswith('variables:'):
                section = 'variables'
            continue
        if section == 'variables':
            if s.lower().startswith('values:'):
                section = 'values'
                var_idx = 0
                continue
            m = re.match(r'\s*\d+\s+(\S+)\s+\S+', line)
            if m:
                variables.append(m.group(1).lower())
        elif section == 'values':
            m_pt = re.match(r'\s*\d+\s+([\d.eE+\-]+)', line)
            if m_pt:
                var_idx = 0
                data.setdefault(variables[0], []).append(float(m_pt.group(1)))
                var_idx = 1
            else:
                m_v = re.match(r'\s+([\d.eE+\-]+)', line)
                if m_v and var_idx < len(variables):
                    data.setdefault(variables[var_idx], []).append(float(m_v.group(1)))
                    var_idx += 1

    v_arr = np.array(data.get('v(1)', data.get('vd', [])))
    # ngspice reports source current as negative (convention: current flows into +)
    i_arr_raw = np.array(data.get('i(vd)', []))
    i_arr = -i_arr_raw   # flip sign: positive when flowing through diode

    return v_arr, i_arr


# ---------------------------------------------------------------------------
# HBFree DC test: write .hbl, run, parse output
# ---------------------------------------------------------------------------

HBL_DC_TEMPLATE = """\
 &SERV EPSIW=1e-06, KITU=0, EPSSOL=1e-12, EPSDU=1e-09, EPSMIN=1e-15,
       MAXDU=100, LIMIT=50, KPRLEN=0, KPRSRT=0, KPRNKR=0,
       KPRLIN=0, KPRSOL=0, MGLOB=1, IAPR=0, KNC=64,
 NAME='diode_dc_test' /

 &CIRCOM /

 &TYP IT='R   ', KOL=2, P=0. /
 &ELEM NE='Rvs ', KNOT=1, 0, PAR={rvs} /
 &ELEM NE='Rl  ', KNOT=2, 0, PAR={rl} /

 &TYP IT='VD  ', 'SCHT', KOL=1, P=1,0.8 /
 &ELEM NE='D1  ', KNOT=1, 2, PAR={rs}, 0, 0, {is_}, 0, {alfa}, 0.5, 0.1 /

 &TYP IT='J   ', KOL=1 /
 &ELEM NE='Vs  ', KNOT=1, 0, PAR=0, 0, 0, {jdc} /

 &TYP IT='END ' /

 &FREQU F1=1, F2=0, MN=0,0, KN=1 /

 &VAR FIN='END ' /

 &QUP IQ='END ' /
"""

ALFA = 1.0 / VT   # = q/(N*k*T)


def run_hbfree_dc(vdc, hbl_dir, hbl_bin='../hbl'):
    """Write a DC-only HBL, run HBFree, return (V1_dc, V2_dc) or None."""
    alfa = ALFA / N
    jdc  = vdc / RVS   # Norton current: V_dc / R_source

    hbl_content = HBL_DC_TEMPLATE.format(
        rvs=RVS, rl=RL, rs=RS,
        is_=IS, alfa=f'{alfa:.6f}',
        jdc=jdc
    )

    hbl_path = os.path.join(hbl_dir, 'diode_dc_test.hbl')
    with open(hbl_path, 'w') as f:
        f.write(hbl_content)

    # run hbl from hbl_dir so it finds diode_dc_test.hbl
    r = subprocess.run(
        [hbl_bin, 'diode_dc_test.hbl'],
        capture_output=True, text=True,
        cwd=hbl_dir
    )
    output = r.stdout + r.stderr

    # parse DC node voltages from output
    # Look for "FREQUENCY  1 COMBINATION(  0,  0)" section
    v1, v2 = _parse_hbfree_dc_output(output)
    return v1, v2, output


def _parse_hbfree_dc_output(text):
    """Extract U(1) and U(2) real parts from the DC frequency block."""
    lines = text.split('\n')
    in_dc = False
    u = {}
    for line in lines:
        if 'COMBINATION(  0,  0)' in line or 'COMBINATION( 0, 0)' in line:
            in_dc = True
            continue
        if in_dc:
            m = re.match(r'\s*U\(\s*(\d+)\)=\s*([\+\-]?[\d\.Ee\+\-]+)\s+([\+\-]?[\d\.Ee\+\-]+)', line)
            if m:
                u[int(m.group(1))] = float(m.group(2))   # real part
            elif line.strip().startswith('FREQUENCY') and u:
                break   # next frequency block
    return u.get(1, None), u.get(2, None)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vdc', type=float, nargs='+', default=[0.40],
                   help='DC verification point voltage(s) (default: 0.40 V); repeat for multiple')
    p.add_argument('--out', default=None,
                   help='Save figure to file (e.g. iv_comparison.png)')
    p.add_argument('--hbl-dir', default=None,
                   help='Directory for HBFree test files (default: auto-detect)')
    p.add_argument('--hbl-bin', default=None,
                   help='Path to hbl executable (default: auto-detect)')
    p.add_argument('--ngspice', default='ngspice',
                   help='Path to ngspice (default: ngspice)')
    p.add_argument('--no-hbfree', action='store_true',
                   help='Skip HBFree DC verification')
    args = p.parse_args()

    # Auto-detect paths relative to this script
    script_dir = os.path.dirname(os.path.realpath(__file__))
    repo_root  = os.path.dirname(script_dir)
    hbl_dir    = args.hbl_dir or os.path.join(repo_root, 'hbfree_exmpls')
    hbl_bin    = args.hbl_bin or os.path.join(repo_root, 'hbl')

    print(f'Parameters: IS={IS:.2e} A  N={N}  RS={RS} Ω  Vt={VT*1e3:.3f} mV  T={T} K')
    print(f'DC verification points: Vdc = {args.vdc}')

    # ---- 1. Analytical I-V ------------------------------------------------
    print('\n[1] Computing analytical I-V...')
    v_analytical = np.linspace(-0.75, 0.65, 800)
    i_analytical = analytical_iv(v_analytical)

    # ---- 2. ngspice DC sweep ----------------------------------------------
    print('[2] Running ngspice DC sweep...')
    v_ng, i_ng = run_ngspice_dc(ngspice=args.ngspice)
    if v_ng is None:
        print('    ngspice failed — skipping')
    else:
        print(f'    Got {len(v_ng)} points, V range [{v_ng.min():.2f}, {v_ng.max():.2f}] V')

    # ---- Analytical DC point solver: Vs -> Rvs -> D1(RS) -> RL -> GND ----
    def full_circuit_current(V_source, tol=1e-15, max_iter=300):
        """Solve I*(Rvs+RS+RL) + Vt*ln(I/IS+1) = V_source via Newton."""
        I = IS * np.exp(min(V_source / (2 * VT), 40.0))   # starting guess
        for _ in range(max_iter):
            Vj    = V_source - (RVS + RS + RL) * I
            Vj    = min(Vj, 40.0 * VT)
            exp_v = np.exp(Vj / VT)
            f     = I - IS * (exp_v - 1.0)
            df    = 1.0 + IS * (RVS + RS + RL) / VT * exp_v
            step  = f / df
            I    -= step
            if abs(step) < tol * (1.0 + abs(I)):
                break
        return I

    # ---- 3. Loop over each requested DC test voltage ----------------------
    dc_points = []   # list of dicts: {vdc, v_term_analyt, i_analyt, hb_v_term, hb_i}
    point_colors = ['tab:blue', 'tab:green', 'tab:purple', 'tab:brown']

    for idx, vdc in enumerate(args.vdc):
        entry = {'vdc': vdc}

        i_analyt_dc   = full_circuit_current(vdc)
        v_term_analyt = vdc - (RVS + RL) * i_analyt_dc
        entry['v_term_analyt'] = v_term_analyt
        entry['i_analyt']      = i_analyt_dc
        print(f'\nAnalytical at Vdc={vdc} V: I={i_analyt_dc*1e6:.4g} µA  V_term={v_term_analyt*1e3:.2f} mV')

        if not args.no_hbfree:
            print(f'[3.{idx+1}] Running HBFree DC test at Vdc = {vdc} V...')
            v1, v2, raw_out = run_hbfree_dc(vdc, hbl_dir, hbl_bin)
            if v1 is not None and v2 is not None:
                entry['hb_i']      = v2 / RL
                entry['hb_v_term'] = v1 - v2
                print(f'    V(1)={v1*1e3:.3f} mV  V(2)={v2*1e3:.3f} mV')
                print(f'    I={entry["hb_i"]*1e6:.4g} µA  V_term={entry["hb_v_term"]*1e3:.3f} mV')
            else:
                entry['hb_i'] = entry['hb_v_term'] = None
                print('    Could not parse HBFree output')
                if raw_out:
                    for line in raw_out.split('\n')[:30]:
                        print('   ', line)
        else:
            entry['hb_i'] = entry['hb_v_term'] = None

        entry['color'] = point_colors[idx % len(point_colors)]
        dc_points.append(entry)

    # ---- 4. Plot ----------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f'Diode I–V comparison  (IS={IS:.0e} A, N={N}, RS={RS} Ω, T={T} K)',
                 fontsize=11)

    for ax, logy in [(ax1, False), (ax2, True)]:
        ax.plot(v_analytical * 1e3, np.abs(i_analytical) * 1e6,
                color='tab:blue', lw=2, label='Analytical (SPICE = HBFree formula)')

        if v_ng is not None:
            ax.plot(v_ng * 1e3, np.abs(i_ng) * 1e6,
                    color='tab:orange', lw=1.5, ls='--', label='ngspice DC sweep')

        for pt in dc_points:
            vdc = pt['vdc']
            col = pt['color']
            ax.scatter([pt['v_term_analyt'] * 1e3], [abs(pt['i_analyt']) * 1e6],
                       color=col, s=80, zorder=5,
                       label=f'Analytical Vdc={vdc:+.2f}V → V_term={pt["v_term_analyt"]*1e3:.0f}mV')
            if pt.get('hb_v_term') is not None:
                ax.scatter([pt['hb_v_term'] * 1e3], [abs(pt['hb_i']) * 1e6],
                           marker='x', s=140, lw=2.5, color=col,
                           zorder=6, label=f'HBFree DC Vdc={vdc:+.2f}V')

        ax.set_xlabel('Terminal voltage  V_anode − V_cathode  (mV)')
        ax.set_ylabel('Diode current  |I|  (µA)')
        ax.grid(True, which='both', alpha=0.35)
        ax.legend(fontsize=9)
        if logy:
            ax.set_yscale('log')
            ax.set_ylim(bottom=1e-6)
            ax.set_title('Log scale')
        else:
            ax.set_title('Linear scale')

    plt.tight_layout()

    if args.out:
        plt.savefig(args.out, dpi=150, bbox_inches='tight')
        print(f'\nSaved: {args.out}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
