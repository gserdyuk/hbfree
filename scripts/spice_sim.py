#!/usr/bin/env python3
"""
spice_sim.py - Run ngspice transient simulation on a HBFree circuit file.

Reads the .hb line to determine fundamental frequencies and harmonic count,
computes appropriate .tran parameters, then runs ngspice and saves the
steady-state waveforms in one or both output formats.

Output formats (--format flag):
  hbfree  Write a HBFree-compatible .raw file (<stem>.spice.raw)
          Phasors in V_peak/2 convention (|S_k| = V_peak/2),
          matching HBFree .raw format; usable directly with cmpraw.py
  npz     Write a NumPy .npz archive (<stem>.tran.npz) with time-domain
          waveforms AND phasor arrays for frequency-domain analysis
  both    Write both (default)

Usage:
  python3 spice_sim.py [options] circuit.ckt

Options:
  --out STEM        output file stem (default: <circuit without .ckt>)
  --format FORMAT   hbfree|npz|both (default: both)
  --cycles N        startup cycles before capture window  (default: 100)
  --samples N       samples per period of highest harmonic  (default: 20)
  --patch-rs        add RS=1 to diode models missing RS (matches HBFree default)
  --ngspice PATH    path to ngspice executable  (default: ngspice)
  --keep            keep the generated .spice/.raw transient files after run
  --verbose         print ngspice stdout
"""

import sys
import re
import os
import subprocess
import argparse
import datetime
import numpy as np


# ---------------------------------------------------------------------------
# Parse .hb line from circuit file
# ---------------------------------------------------------------------------

def parse_hb(lines):
    """Return (f1, f2, harmonics_list) from the .hb line(s).

    Handles SPICE-style '+' continuation lines and ignores .hb_options.
    """
    hb_text = []
    in_hb = False
    for line in lines:
        s = line.strip()
        sl = s.lower()
        if sl.startswith('.hb') and not sl.startswith('.hb_'):
            in_hb = True
            hb_text.append(s)
        elif in_hb:
            if s.startswith('+'):
                hb_text.append(s[1:])   # strip leading '+'
            else:
                in_hb = False

    full = ' '.join(hb_text)

    f1 = 0.0
    f2 = 0.0
    m = re.search(r'f1\s*=\s*([\d.eE+\-]+)', full, re.IGNORECASE)
    if m:
        f1 = float(m.group(1))
    m = re.search(r'f2\s*=\s*([\d.eE+\-]+)', full, re.IGNORECASE)
    if m:
        f2 = float(m.group(1))

    # Strip key=value tokens before scanning for harmonic pairs,
    # otherwise "f2=0  0 0 1 0" bleeds into the pair list.
    pairs_text = re.sub(r'\w+\s*=\s*\S+', '', full)
    # Use -?\d+ to handle negative second indices (e.g. "1 -3" in two-tone circuits).
    pairs = re.findall(r'(-?\d+)\s+(-?\d+)', pairs_text)
    harmonics = [(int(a), int(b)) for a, b in pairs]

    return f1, f2, harmonics


_SPICE_KW = re.compile(
    r'^(ac|dc|sin|pulse|exp|pwl|sffm|am|tran|distof\d*)\b', re.I)

_DISTOF2_RE = re.compile(r'\bdistof2\s+([\d.eE+\-]+)', re.I)

# Element types that have a model name as their last node-like token
_HAS_MODEL = set('DQJMKEFGH')


def find_2tone_sources(ckt_lines):
    """Return list of dicts for V-sources that carry a distof2 parameter."""
    result = []
    for line in ckt_lines[1:]:  # skip title line
        s = line.strip()
        if not s or s.startswith('*') or s.startswith('.') or s.startswith('+'):
            continue
        if s[0].upper() != 'V':
            continue
        m = _DISTOF2_RE.search(s)
        if m:
            parts = s.split()
            if len(parts) >= 3:
                result.append(dict(name=parts[0], node1=parts[1],
                                   node2=parts[2], distof2=float(m.group(1))))
    return result

def parse_nodes(ckt_lines):
    """Collect node names from element lines (skip title, comments, directives)."""
    nodes = set()
    for line in ckt_lines[1:]:  # skip line 0 (SPICE title)
        s = line.strip()
        if not s or s.startswith('*') or s.startswith('.') or s.startswith('+'):
            continue
        parts = s.split()
        if len(parts) < 2:
            continue
        has_model = parts[0][0].upper() in _HAS_MODEL
        collected = []
        for p in parts[1:]:
            if _SPICE_KW.match(p) or '=' in p or '(' in p:
                break   # keyword or param — node list ends
            if re.match(r'^\d+$', p):
                collected.append(p)         # pure integer node (e.g. "1", "2")
            elif re.match(r'^\d', p):
                break                       # digit-prefixed value (e.g. "1k") — stop
            elif re.match(r'^[\w#]+$', p):
                collected.append(p.lower()) # alphabetic node name (e.g. "out")
            else:
                break   # float or other value — stop
        if has_model and collected:
            collected = collected[:-1]      # last token was model name, not a node
        for n in collected:
            nodes.add(n)
    nodes.discard('0')
    nodes.discard('gnd')
    return sorted(nodes)


# ---------------------------------------------------------------------------
# Transient parameter computation
# ---------------------------------------------------------------------------

def tran_params(f1, f2, harmonics, cycles_startup, samples_per_period):
    """Return (tstep, tstart, tstop, t_window_len)."""
    if harmonics:
        max_h1 = max(abs(a) for a, b in harmonics)
        max_h2 = max(abs(b) for a, b in harmonics) if f2 > 0 else 0
    else:
        max_h1, max_h2 = 9, 0

    max_freq = max(max_h1 * f1, max_h2 * f2 if f2 > 0 else 0, f1)
    tstep = 1.0 / (samples_per_period * max_freq)

    if f2 > 0 and abs(f2 - f1) > 1:
        t_window = 1.0 / abs(f2 - f1)
    else:
        t_window = 1.0 / f1

    t_startup = cycles_startup / f1
    tstart    = t_startup
    tstop     = t_startup + t_window

    return tstep, tstart, tstop, t_window


# ---------------------------------------------------------------------------
# Patch RS into diode models
# ---------------------------------------------------------------------------

def patch_rs(lines, rsval=1.0):
    """Insert RS=<rsval> into .model D lines that don't already have RS=."""
    out = []
    for line in lines:
        s = line.strip()
        if re.match(r'\.model\b', s, re.IGNORECASE):
            parts = s.split()
            if len(parts) >= 3 and parts[2].upper() == 'D':
                if 'rs' not in s.lower():
                    line = line.rstrip() + f' RS={rsval}\n'
        out.append(line)
    return out


# ---------------------------------------------------------------------------
# Build modified spice deck
# ---------------------------------------------------------------------------

def build_spice(ckt_lines, tstep, tstart, tstop, nodes, patch_rs_val=None, f2=0.0):
    """Return list of lines for the modified spice deck.

    For two-tone circuits (f2 > 0): locates the V-source carrying distof2,
    inserts an intermediate node, and adds a series voltage source at f2 so
    that ngspice sees both tones.  The f2 amplitude is 2*distof2 (HBFree
    convention: V_peak = 2*distof).
    """
    skip_keywords = ('.ac', '.tran', '.dc', '.hb', '.hb_options',
                     '.options', '.end')

    if patch_rs_val is not None:
        ckt_lines = patch_rs(ckt_lines, patch_rs_val)

    # Build map: src_name.lower() -> (mid_node, orig_node2, f2_amp)
    src_mods = {}
    if f2 > 0:
        for s in find_2tone_sources(ckt_lines):
            skey = s['name'].lower()
            src_mods[skey] = (f'_f2int_{skey}', s['node2'], 2.0 * s['distof2'])

    deck = []
    for line in ckt_lines:
        s = line.strip().lower()
        if any(s.startswith(k) for k in skip_keywords):
            continue
        if s == '.end':
            continue
        # Two-tone: redirect node2 of the f1 source and add series f2 source
        if src_mods:
            parts = line.strip().split(None, 3)
            if len(parts) >= 3 and parts[0].lower() in src_mods:
                skey = parts[0].lower()
                mid_node, orig_node2, f2_amp = src_mods[skey]
                rest = parts[3] if len(parts) > 3 else ''
                deck.append(f'{parts[0]} {parts[1]} {mid_node} {rest}\n')
                deck.append(f'v_f2_{skey} {mid_node} {orig_node2}'
                            f' sin(0 {f2_amp:.6g} {f2:.6g})\n')
                continue
        deck.append(line)

    deck.append(f'\n* --- added by spice_sim.py ---\n')
    deck.append(f'.options filetype=ascii\n')
    deck.append(f'.tran {tstep:.6e} {tstop:.6e} {tstart:.6e}\n')
    if nodes:
        deck.append('.save ' + ' '.join(f'v({n})' for n in nodes) + '\n')
    deck.append('.end\n')
    return deck


# ---------------------------------------------------------------------------
# Parse ngspice ASCII raw output
# ---------------------------------------------------------------------------

def parse_ngspice_raw(rawfile):
    """Return dict: {node_name: np.array}, including 'time'."""
    with open(rawfile, encoding='latin-1') as f:
        lines = f.readlines()

    variables = []
    values = {}
    section = 'header'
    var_idx = 0

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
            m = re.match(r'\s*(\d+)\s+(\S+)\s+(\S+)', line)
            if m:
                variables.append(m.group(2).lower())
        elif section == 'values':
            m_pt = re.match(r'\s*(\d+)\s+([\d.eE+\-]+)', line)
            if m_pt:
                var_idx = 0
                val = float(m_pt.group(2))
                name = variables[0] if variables else 'time'
                values.setdefault(name, []).append(val)
                var_idx = 1
            else:
                m_val = re.match(r'\s+([\d.eE+\-]+)', line)
                if m_val and var_idx < len(variables):
                    name = variables[var_idx]
                    values.setdefault(name, []).append(float(m_val.group(1)))
                    var_idx += 1

    return {k: np.array(v) for k, v in values.items()}


# ---------------------------------------------------------------------------
# FFT to phasors
# ---------------------------------------------------------------------------

def fft_phasors(time_arr, volt_arr, f1, max_harmonic, tstart, tstop):
    """
    Extract complex phasors S_k from transient data using FFT.

    Output convention matches HBFree .raw: |S_k| = V_peak / 2.

    Chain:
      rfft gives V[k] = N * V_peak / 2  (one-sided, for a cosine)
      divide by N  -> V_peak / 2

    DC (k=0): S_0 = FFT[0] / N  (real, mean value).

    To recover peak voltage from a phasor: V_peak = 2 * |S_k|.

    Returns dict {k: complex_phasor} for k = 0 .. max_harmonic.
    """
    mask = (time_arr >= tstart - 1e-15) & (time_arr <= tstop + 1e-15)
    t = time_arr[mask]
    v = volt_arr[mask]
    if len(t) < 4:
        return {k: 0j for k in range(max_harmonic + 1)}

    # Interpolate to uniform grid (ngspice uses adaptive timestep)
    N = len(t)
    t_uniform = np.linspace(t[0], t[-1], N)
    v_uniform = np.interp(t_uniform, t, v)

    V = np.fft.rfft(v_uniform)
    T_window = t[-1] - t[0]

    phasors = {}
    for k in range(max_harmonic + 1):
        b = round(k * f1 * T_window)
        if b >= len(V):
            phasors[k] = 0j
        else:
            phasors[k] = complex(V[b]) / N
    return phasors


# ---------------------------------------------------------------------------
# Write HBFree-compatible .raw file
# ---------------------------------------------------------------------------

def write_hbfree_raw(out_path, nodes, phasors_by_node, f1, max_harmonic):
    """
    Write a HBFree-compatible .raw file with phasors from FFT.

    Structure mirrors wrtraw() in main.for: 3 records per frequency group
    (zeros / phasors / zeros), for k = 0 (DC) through max_harmonic.
    Node variables are written as V(node) to match HBFree output.

    Convention (matches HBFree after jdrive fix):
      DC  (k=0): S_0 = mean value (no scaling)
      AC  (k>=1): S_k = (V_peak/2) * sqrt(2) = V_rms  (HBFree stores V_rms)
    """
    Kn = max_harmonic + 1  # DC + harmonics
    npoints = 3 * Kn
    nvar = len(nodes)
    now = datetime.datetime.now().strftime('%d.%m.%Y %H:%M:%S')
    ac_scale = float(np.sqrt(2.0))  # V_peak/2 -> V_rms for AC harmonics

    with open(out_path, 'w') as f:
        f.write('Title: ngspice transient -> harmonic balance\n')
        f.write(f'Date: {now}\n')
        f.write('Plotname:  Harmonic Balance Simulation\n')
        f.write('Flags: complex\n')
        f.write(f'No. Variables:    {nvar + 1}\n')
        f.write(f'No. Points:   {npoints}\n')
        f.write('Command:  version spice_sim\n')
        f.write('Variables:\n')
        f.write('      0    frequency    frequency\n')
        for i, node in enumerate(nodes):
            f.write(f'      {i+1}    V({node})    voltage\n')
        f.write('Values:\n')

        pt = 0
        zero_val = '+0.0000000000000E+00,+0.0000000000000E+00'
        for k in range(Kn):
            freq = k * f1
            freq_hdr = f'{freq:+.13E},+0.0000000000000E+00'
            for rec in range(3):
                f.write(f'{pt:5d}    {freq_hdr}\n')
                for node in nodes:
                    if rec == 1:
                        ph = phasors_by_node.get(node, {}).get(k, 0j)
                        if k > 0:
                            ph = ph * ac_scale  # V_peak/2 -> V_rms
                        re_v = float(np.real(ph))
                        im_v = float(np.imag(ph))
                        f.write(f'         {re_v:+.13E},{im_v:+.13E}\n')
                    else:
                        f.write(f'         {zero_val}\n')
                pt += 1


# ---------------------------------------------------------------------------
# Two-tone FFT and .raw writer
# ---------------------------------------------------------------------------

def fft_phasors_2tone(time_arr, volt_arr, f1, f2, harmonics, tstart, tstop):
    """Extract phasors for all (n1, n2) pairs from a two-tone transient waveform.

    Key insight: with window T = 1/(f2-f1), each n1*f1+n2*f2 product lands on a
    unique FFT bin (bin = round(|freq| / df) where df = 1/T).

    Convention matches fft_phasors: |S| = V_peak/2 for AC, mean for DC.
    For negative-frequency pairs (freq < 0), returns conj of the +|freq| bin.
    """
    mask = (time_arr >= tstart - 1e-15) & (time_arr <= tstop + 1e-15)
    t = time_arr[mask]
    v = volt_arr[mask]
    if len(t) < 4:
        return {pair: 0j for pair in harmonics}

    N = len(t)
    t_uniform = np.linspace(t[0], t[-1], N)
    v_uniform = np.interp(t_uniform, t, v)

    V = np.fft.rfft(v_uniform)
    T_window = t[-1] - t[0]
    df = 1.0 / T_window  # bin spacing ≈ f2 - f1

    phasors = {}
    for (n1, n2) in harmonics:
        if n1 == 0 and n2 == 0:
            phasors[(0, 0)] = complex(V[0]) / N  # DC
            continue
        freq = n1 * f1 + n2 * f2
        bin_idx = int(round(abs(freq) / df))
        if bin_idx >= len(V):
            phasors[(n1, n2)] = 0j
        elif freq < 0:
            phasors[(n1, n2)] = np.conj(complex(V[bin_idx])) / N
        else:
            phasors[(n1, n2)] = complex(V[bin_idx]) / N
    return phasors


def write_hbfree_raw_2tone(out_path, nodes, phasors_by_node, f1, f2, harmonics):
    """Write a HBFree-compatible .raw for a two-tone simulation.

    harmonics: list of (n1, n2) pairs as returned by parse_hb (includes DC (0,0)).
    phasors_by_node: {node: {(n1,n2): complex}} from fft_phasors_2tone.
    Frequencies written as signed n1*f1 + n2*f2, matching HBFree .raw format.
    """
    Kn = len(harmonics)
    npoints = 3 * Kn
    nvar = len(nodes)
    now = datetime.datetime.now().strftime('%d.%m.%Y %H:%M:%S')
    ac_scale = float(np.sqrt(2.0))  # V_peak/2 -> V_rms for AC harmonics

    with open(out_path, 'w') as f:
        f.write('Title: ngspice transient -> harmonic balance (two-tone)\n')
        f.write(f'Date: {now}\n')
        f.write('Plotname:  Harmonic Balance Simulation\n')
        f.write('Flags: complex\n')
        f.write(f'No. Variables:    {nvar + 1}\n')
        f.write(f'No. Points:   {npoints}\n')
        f.write('Command:  version spice_sim\n')
        f.write('Variables:\n')
        f.write('      0    frequency    frequency\n')
        for i, node in enumerate(nodes):
            f.write(f'      {i+1}    V({node})    voltage\n')
        f.write('Values:\n')

        pt = 0
        zero_val = '+0.0000000000000E+00,+0.0000000000000E+00'
        for (n1, n2) in harmonics:
            freq = n1 * f1 + n2 * f2
            is_dc = (n1 == 0 and n2 == 0)
            freq_hdr = f'{freq:+.13E},+0.0000000000000E+00'
            for rec in range(3):
                f.write(f'{pt:5d}    {freq_hdr}\n')
                for node in nodes:
                    if rec == 1:
                        ph = phasors_by_node.get(node, {}).get((n1, n2), 0j)
                        if not is_dc:
                            ph = ph * ac_scale
                        f.write(f'         {float(np.real(ph)):+.13E},'
                                f'{float(np.imag(ph)):+.13E}\n')
                    else:
                        f.write(f'         {zero_val}\n')
                pt += 1


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('ckt', help='Circuit file (.ckt)')
    p.add_argument('--out', default=None,
                   help='Output file stem (default: <circuit without .ckt>)')
    p.add_argument('--format', default='both',
                   choices=['hbfree', 'npz', 'both'],
                   help='Output format: hbfree|npz|both (default: both)')
    p.add_argument('--cycles', type=int, default=100,
                   help='Startup cycles before capture (default: 100)')
    p.add_argument('--samples', type=int, default=20,
                   help='Samples per period of highest harmonic (default: 20)')
    p.add_argument('--patch-rs', action='store_true',
                   help='Add RS=1 to diode models missing RS')
    p.add_argument('--ngspice', default='ngspice',
                   help='Path to ngspice executable')
    p.add_argument('--keep', action='store_true',
                   help='Keep generated transient .spice/.raw files')
    p.add_argument('--verbose', action='store_true',
                   help='Print ngspice stdout')
    args = p.parse_args()

    ckt_path = os.path.realpath(args.ckt)
    if not os.path.isfile(ckt_path):
        sys.exit(f'File not found: {args.ckt}')

    with open(ckt_path) as f:
        ckt_lines = f.readlines()

    f1, f2, harmonics = parse_hb(ckt_lines)
    if f1 == 0:
        sys.exit('Could not find f1 in .hb line')

    print(f'Circuit : {os.path.basename(ckt_path)}')
    print(f'f1      : {f1/1e6:.3f} MHz')
    if f2 > 0:
        print(f'f2      : {f2/1e6:.3f} MHz')
    print(f'Harmonics: {len(harmonics)} pairs')

    tstep, tstart, tstop, t_window = tran_params(
        f1, f2, harmonics, args.cycles, args.samples)

    print(f'tran    : step={tstep:.3e}s  start={tstart:.3e}s  stop={tstop:.3e}s')
    print(f'Window  : {t_window:.3e}s  ({t_window*f1:.0f} periods of f1)')

    nodes = parse_nodes(ckt_lines)
    print(f'Nodes   : {nodes}')

    deck = build_spice(ckt_lines, tstep, tstart, tstop, nodes,
                       patch_rs_val=1.0 if args.patch_rs else None, f2=f2)

    ckt_dir  = os.path.dirname(ckt_path)
    ckt_base = os.path.basename(ckt_path)
    stem = args.out or os.path.join(ckt_dir, ckt_base.replace('.ckt', ''))

    spice_file = stem + '.spice_sim'
    tran_raw   = stem + '.spice_sim.raw'
    hbfree_raw = stem + '.spice.raw'
    npz_file   = stem + '.tran.npz'

    with open(spice_file, 'w') as f:
        f.writelines(deck)

    cmd = [args.ngspice, '-b', '-r', tran_raw, spice_file]
    print(f'\nRunning : {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=not args.verbose, text=True)
    if args.verbose and result.stdout:
        print(result.stdout)
    if result.returncode != 0:
        if not args.verbose and result.stderr:
            print(result.stderr, file=sys.stderr)
        sys.exit(f'ngspice failed (exit {result.returncode})')

    if not os.path.isfile(tran_raw):
        sys.exit(f'ngspice produced no output: {tran_raw}')

    data = parse_ngspice_raw(tran_raw)
    if 'time' not in data:
        sys.exit('Could not parse time from ngspice output')

    time_arr = data['time']
    print(f'Points  : {len(time_arr)}  nodes: {[k for k in data if k != "time"]}')

    if f2 > 0:
        # Two-tone path: extract phasors at every (n1, n2) harmonic pair.
        phasors_by_node = {}
        for node in nodes:
            key = f'v({node})'
            if key in data:
                phasors_by_node[node] = fft_phasors_2tone(
                    time_arr, data[key], f1, f2, harmonics, tstart, tstop)
            else:
                phasors_by_node[node] = {pair: 0j for pair in harmonics}

        if args.format in ('hbfree', 'both'):
            write_hbfree_raw_2tone(hbfree_raw, nodes, phasors_by_node,
                                   f1, f2, harmonics)
            print(f'HBFree  : {hbfree_raw}')

        if args.format in ('npz', 'both'):
            npz_data = dict(f1=f1, f2=f2, t_window=t_window,
                            harmonics_pairs=str(harmonics))
            for k, v in data.items():
                npz_data[k] = v
            np.savez(npz_file, **npz_data)
            print(f'npz     : {npz_file}')
    else:
        # Single-tone path (unchanged).
        max_h = max((abs(a) for a, b in harmonics), default=9)

        phasors_by_node = {}
        for node in nodes:
            key = f'v({node})'
            if key in data:
                phasors_by_node[node] = fft_phasors(
                    time_arr, data[key], f1, max_h, tstart, tstop)
            else:
                phasors_by_node[node] = {k: 0j for k in range(max_h + 1)}

        if args.format in ('hbfree', 'both'):
            write_hbfree_raw(hbfree_raw, nodes, phasors_by_node, f1, max_h)
            print(f'HBFree  : {hbfree_raw}')

        if args.format in ('npz', 'both'):
            npz_data = dict(f1=f1, f2=f2, t_window=t_window,
                            harmonics_max=max_h)
            for k, v in data.items():
                npz_data[k] = v
            for node, ph_dict in phasors_by_node.items():
                ph_arr = np.array([ph_dict.get(k, 0j) for k in range(max_h + 1)])
                safe = re.sub(r'[^\w]', '_', node)
                npz_data[f'ph_re_{safe}'] = np.real(ph_arr)
                npz_data[f'ph_im_{safe}'] = np.imag(ph_arr)
            np.savez(npz_file, **npz_data)
            print(f'npz     : {npz_file}')

    if not args.keep:
        os.unlink(spice_file)
        os.unlink(tran_raw)


if __name__ == '__main__':
    main()
