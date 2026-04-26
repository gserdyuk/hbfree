#!/usr/bin/env python3
"""
plot_timedomain.py  —  Overlay HBFree and ngspice waveforms in the time domain.

Usage:
  python3 scripts/plot_timedomain.py <base>           # uses <base>.raw + <base>.tran.npz
  python3 scripts/plot_timedomain.py <base>.raw <base>.tran.npz
  python3 scripts/plot_timedomain.py --help

HBFree .raw stores V_rms phasors (|S_k| = V_rms = V_peak/√2 for AC harmonics; DC as actual value).
Time-domain reconstruction: v(t) = V_DC + √2 * Re( Σ_k  V_k · exp(j·k·ω₁·t) )
ngspice transient data comes directly from the .tran.npz window array.
"""

import sys
import re
import argparse
import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Raw file parser (minimal, reused from cmpraw.py logic)
# ---------------------------------------------------------------------------

def parse_raw(filename):
    """Return (variables, harmonics, f1_hz).
    variables: list of {'idx', 'name'}
    harmonics: list of {'freq_hz', 'values': [complex, ...]}
    f1_hz: fundamental frequency (Hz)
    """
    with open(filename, encoding='latin-1') as f:
        lines = f.readlines()

    variables = []
    points = {}
    current_pt = None
    section = 'header'

    for raw_line in lines:
        line = raw_line.rstrip('\n')
        stripped = line.strip()
        if not stripped:
            continue

        if section == 'header':
            if stripped.lower().startswith('variables:'):
                section = 'variables'
            elif stripped.lower().startswith('values:'):
                section = 'values'

        elif section == 'variables':
            if stripped.lower().startswith('values:'):
                section = 'values'
                continue
            m = re.match(r'\s*(\d+)\s+(\S+)\s+(\S+)', line)
            if m:
                variables.append({'idx': int(m.group(1)), 'name': m.group(2)})

        elif section == 'values':
            m_pt = re.match(
                r'\s{0,4}(\d+)\s+([\+\-]?[\d\.Ee\+\-]+),([\+\-]?[\d\.Ee\+\-]+)', line)
            if m_pt:
                pt_idx = int(m_pt.group(1))
                current_pt = {
                    'freq': complex(float(m_pt.group(2)), float(m_pt.group(3))),
                    'values': [],
                }
                points[pt_idx] = current_pt
                continue
            m_val = re.match(
                r'\s+([\+\-]?[\d\.Ee\+\-]+),([\+\-]?[\d\.Ee\+\-]+)', line)
            if m_val and current_pt is not None:
                current_pt['values'].append(
                    complex(float(m_val.group(1)), float(m_val.group(2))))

    # Phasor records at indices 1, 4, 7, ... (3k+1)
    harmonics = []
    k = 0
    while True:
        idx = 3 * k + 1
        if idx not in points:
            break
        pt = points[idx]
        harmonics.append({'freq_hz': pt['freq'].real, 'values': pt['values']})
        k += 1

    # f1 is frequency of harmonic index 1 (0 = DC)
    f1_hz = harmonics[1]['freq_hz'] if len(harmonics) > 1 else 0.0
    return variables, harmonics, f1_hz


# ---------------------------------------------------------------------------
# Time-domain reconstruction from one-sided phasors
# ---------------------------------------------------------------------------

def phasors_to_timedomain(variables, harmonics, node_name, t):
    """Reconstruct v(t) for node_name from HBFree .raw phasors.

    Convention: DC stored as actual value; AC stored as V_rms = V_peak/√2.
    v(t) = V_DC  +  √2 * Re( Σ_{k≥1}  V_k · exp(j · k · ω₁ · t) )
    """
    # Find variable index (0-based into values list = idx-1)
    norm_target = _norm(node_name)
    var_offset = None
    for v in variables:
        if _norm(v['name']) == norm_target:
            var_offset = v['idx'] - 1   # values list is 0-based, idx is 1-based
            break
    if var_offset is None:
        return None

    # f1 from harmonic index 1
    f1_hz = harmonics[1]['freq_hz'] if len(harmonics) > 1 else 0.0
    omega1 = 2 * np.pi * f1_hz

    v_t = np.zeros(len(t), dtype=complex)
    for ki, h in enumerate(harmonics):
        if var_offset >= len(h['values']):
            continue
        V_k = h['values'][var_offset]
        if ki == 0:
            v_t += V_k.real             # DC: actual value (real)
        else:
            v_t += np.sqrt(2.0) * V_k * np.exp(1j * ki * omega1 * t)

    return v_t.real


def _norm(name):
    n = name.strip().lower()
    m = re.match(r'^v\((.+)\)$', n)
    if m:
        n = m.group(1)
    return n


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('raw_or_base',
                   help='HBFree .raw file  OR  base stem (e.g. di1_s2)')
    p.add_argument('npz', nargs='?', default=None,
                   help='.tran.npz file (auto-derived if omitted)')
    p.add_argument('--out', default=None,
                   help='Save figure to file (e.g. --out fig.png); show interactively if omitted')
    p.add_argument('--nodes', default=None,
                   help='Comma-separated node names to plot (default: all shared nodes)')
    p.add_argument('--periods', type=float, default=1.0,
                   help='Number of periods to display (default: 1)')
    args = p.parse_args()

    # Resolve filenames
    raw_file = args.raw_or_base
    npz_file = args.npz
    if not raw_file.endswith('.raw'):
        base = raw_file
        raw_file = base + '.raw'
        if npz_file is None:
            npz_file = base + '.tran.npz'
    elif npz_file is None:
        # derive from raw name: strip .raw, add .tran.npz
        npz_file = re.sub(r'\.raw$', '.tran.npz', raw_file)

    # Load HBFree raw
    try:
        hb_vars, hb_harmonics, f1_hz = parse_raw(raw_file)
    except OSError as e:
        sys.exit(f'Cannot open {raw_file}: {e}')
    hb_node_names = [v['name'] for v in hb_vars if _norm(v['name']) != 'frequency']

    # Load ngspice transient
    try:
        npz = np.load(npz_file, allow_pickle=True)
    except OSError as e:
        sys.exit(f'Cannot open {npz_file}: {e}')

    t_ng = npz['time']
    f1_ng = float(npz['f1'])
    T1 = 1.0 / f1_ng
    # ngspice waveform covers exactly one period in the window
    t_start = t_ng[0]
    t_end   = t_ng[-1]

    # Build list of nodes to plot
    if args.nodes:
        plot_nodes = [n.strip() for n in args.nodes.split(',')]
    else:
        # Intersect HBFree node names with keys in npz
        npz_v_keys = {k for k in npz.keys() if k.startswith('v(')}
        plot_nodes = []
        for name in hb_node_names:
            ng_key = f'v({_norm(name)})'
            if ng_key in npz_v_keys:
                plot_nodes.append(name)
        if not plot_nodes:
            print(f'HBFree nodes : {hb_node_names}')
            print(f'ngspice keys : {sorted(npz_v_keys)}')
            sys.exit('No matching nodes found between .raw and .tran.npz')

    # Time axis for HBFree reconstruction — one period centred on the same window
    N_pts = 500
    t_hb = np.linspace(t_start, t_start + args.periods * T1, N_pts)

    # Plot
    n_rows = len(plot_nodes)
    fig, axes = plt.subplots(n_rows, 1, figsize=(10, 3.5 * n_rows), sharex=False)
    if n_rows == 1:
        axes = [axes]

    for ax, node in zip(axes, plot_nodes):
        # HBFree reconstructed
        v_hb = phasors_to_timedomain(hb_vars, hb_harmonics, node, t_hb)

        # ngspice transient (one period window)
        ng_key = f'v({_norm(node)})'
        v_ng = npz[ng_key]

        t_hb_ns = t_hb * 1e9
        t_ng_ns = t_ng * 1e9

        if v_hb is not None:
            ax.plot(t_hb_ns, v_hb * 1e3, label='HBFree (reconstructed)', color='tab:blue',
                    linewidth=2)
        ax.plot(t_ng_ns, v_ng * 1e3, label='ngspice (transient)', color='tab:orange',
                linewidth=1.5, linestyle='--')

        ax.set_ylabel('Voltage (mV)')
        ax.set_title(f'V({_norm(node)})')
        ax.legend(loc='upper right', fontsize=9)
        ax.grid(True, alpha=0.4)
        ax.set_xlabel('Time (ns)')

    fig.suptitle(f'Time-domain comparison  —  {raw_file}', fontsize=11)
    plt.tight_layout()

    if args.out:
        plt.savefig(args.out, dpi=150, bbox_inches='tight')
        print(f'Saved: {args.out}')
    else:
        plt.show()


if __name__ == '__main__':
    main()
