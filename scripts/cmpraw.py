#!/usr/bin/env python3
"""
cmpraw.py - Numerically compare two HBFree .raw files.

Understands the complex-phasor format written by wrtraw() in main.for:
  each frequency produces 3 records (DC-zeros, +phasor, -zeros),
  so the useful data is at point indices 1, 4, 7, ... (3k+1).

Usage:
  python3 cmpraw.py [options] reference.raw test.raw

Options:
  --atol FLOAT   absolute tolerance          (default: 1e-9)
  --rtol FLOAT   relative tolerance          (default: 1e-3)
  --mag          compare magnitudes only, ignore phase
  --vars A,B,..  restrict comparison to these variable names
  --summary      print summary line only, no per-point table
  --quiet        exit 0 on pass, 1 on fail, minimal output
"""

import sys
import re
import argparse
import math


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def parse_raw(filename):
    """Return dict: {header, variables, harmonics}
    harmonics: list of {freq_hz, values:[complex,...]}  — one entry per frequency.
    variables: list of {idx, name, kind}
    """
    try:
        with open(filename, encoding='latin-1') as f:
            lines = f.readlines()
    except OSError as e:
        sys.exit(f"Cannot open {filename}: {e}")

    header = {}
    variables = []
    points = {}      # int -> {freq: complex, values: [complex]}
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
                continue
            if stripped.lower().startswith('values:'):
                section = 'values'
                continue
            if ':' in stripped:
                key, _, val = stripped.partition(':')
                header[key.strip()] = val.strip()

        elif section == 'variables':
            if stripped.lower().startswith('values:'):
                section = 'values'
                continue
            m = re.match(r'\s*(\d+)\s+(\S+)\s+(\S+)', line)
            if m:
                variables.append({
                    'idx': int(m.group(1)),
                    'name': m.group(2),
                    'kind': m.group(3),
                })

        elif section == 'values':
            # New point header: starts with optional spaces then an integer
            m_pt = re.match(r'\s{0,4}(\d+)\s+([\+\-]?[\d\.Ee\+\-]+),([\+\-]?[\d\.Ee\+\-]+)', line)
            if m_pt:
                pt_idx = int(m_pt.group(1))
                current_pt = {
                    'freq': complex(float(m_pt.group(2)), float(m_pt.group(3))),
                    'values': [],
                }
                points[pt_idx] = current_pt
                continue
            # Continuation line (variable value)
            m_val = re.match(r'\s+([\+\-]?[\d\.Ee\+\-]+),([\+\-]?[\d\.Ee\+\-]+)', line)
            if m_val and current_pt is not None:
                current_pt['values'].append(
                    complex(float(m_val.group(1)), float(m_val.group(2)))
                )

    # Extract harmonics: phasor records are at indices 1, 4, 7, ... (3k+1)
    harmonics = []
    k = 0
    while True:
        idx = 3 * k + 1
        if idx not in points:
            break
        pt = points[idx]
        harmonics.append({'freq_hz': pt['freq'].real, 'values': pt['values']})
        k += 1

    return {'header': header, 'variables': variables, 'harmonics': harmonics}


# ---------------------------------------------------------------------------
# Variable name normalisation
# ---------------------------------------------------------------------------

def _norm(name):
    """Canonical form: lower-case, strip V( ) wrapper."""
    n = name.strip().lower()
    # V(node) -> node
    m = re.match(r'^v\((.+)\)$', n)
    if m:
        n = m.group(1)
    # #internal_variable_N -> #intNN
    m = re.match(r'^#internal_variable_(\d+)$', n)
    if m:
        n = f'#int{int(m.group(1)):02d}'
    return n


def _match_variables(vars1, vars2):
    """Return list of (name_in_1, idx1, name_in_2, idx2) for matching vars."""
    norm1 = {_norm(v['name']): (v['name'], v['idx']) for v in vars1}
    norm2 = {_norm(v['name']): (v['name'], v['idx']) for v in vars2}
    common_norms = set(norm1) & set(norm2)
    matched = []
    for n in sorted(common_norms):
        n1, i1 = norm1[n]
        n2, i2 = norm2[n]
        matched.append((n1, i1, n2, i2))
    only1 = [norm1[n][0] for n in sorted(set(norm1) - common_norms)]
    only2 = [norm2[n][0] for n in sorted(set(norm2) - common_norms)]
    return matched, only1, only2


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def _fmt_freq(f_hz):
    if f_hz == 0:
        return 'DC      '
    for unit, scale in [('GHz', 1e9), ('MHz', 1e6), ('kHz', 1e3), ('Hz', 1)]:
        if f_hz >= scale * 0.999:
            return f'{f_hz/scale:.4g} {unit}'
    return f'{f_hz:.3e} Hz'


def _rel_diff(z1, z2, atol):
    """Relative complex difference: |z1-z2| / max(|z1|, |z2|, atol)."""
    denom = max(abs(z1), abs(z2), atol)
    return abs(z1 - z2) / denom


def compare(ref_file, test_file, atol=1e-9, rtol=1e-3,
            mag_only=False, filter_vars=None, summary_only=False, quiet=False):

    ref = parse_raw(ref_file)
    tst = parse_raw(test_file)

    matched, only_ref, only_tst = _match_variables(ref['variables'], tst['variables'])

    if filter_vars:
        fset = {_norm(v) for v in filter_vars.split(',')}
        matched = [(n1, i1, n2, i2) for n1, i1, n2, i2 in matched if _norm(n1) in fset]

    # Exclude 'frequency' pseudo-variable from value comparison
    cmp_vars = [(n1, i1, n2, i2) for n1, i1, n2, i2 in matched
                if _norm(n1) != 'frequency']

    # Match harmonics by frequency value
    ref_freqs = [h['freq_hz'] for h in ref['harmonics']]
    tst_freqs = [h['freq_hz'] for h in tst['harmonics']]

    freq_pairs = []   # (ref_harmonic_idx, tst_harmonic_idx)
    for ri, rf in enumerate(ref_freqs):
        for ti, tf in enumerate(tst_freqs):
            denom = max(abs(rf), abs(tf), 1.0)
            if abs(rf - tf) / denom < 0.01:   # 1% frequency match
                freq_pairs.append((ri, ti))
                break

    passes = 0
    fails = 0
    fail_details = []

    col_w = [10, 18, 14, 14, 12, 10, 7]

    if not summary_only and not quiet:
        print(f'\n{"Freq":10s}  {"Variable":18s}  {"Ref |z|":>14s}  {"Test |z|":>14s}  {"adiff":>12s}  {"rdiff":>10s}  Status')
        print('-' * 85)

    for ri, ti in freq_pairs:
        rh = ref['harmonics'][ri]
        th = tst['harmonics'][ti]
        freq_str = _fmt_freq(rh['freq_hz'])

        for n1, i1, n2, i2 in cmp_vars:
            if i1 - 1 >= len(rh['values']) or i2 - 1 >= len(th['values']):
                continue

            z1 = rh['values'][i1 - 1]   # idx is 1-based, list is 0-based; idx 0 = frequency
            z2 = th['values'][i2 - 1]

            if mag_only:
                z1_cmp = complex(abs(z1), 0)
                z2_cmp = complex(abs(z2), 0)
            else:
                z1_cmp, z2_cmp = z1, z2

            adiff = abs(z1_cmp - z2_cmp)
            rdiff = _rel_diff(z1_cmp, z2_cmp, atol)

            ok = adiff <= atol or rdiff <= rtol
            if ok:
                passes += 1
                status = 'PASS'
            else:
                fails += 1
                status = 'FAIL'
                fail_details.append(
                    f'  {freq_str:10s}  {n1:18s}  ref={z1:.4e}  test={z2:.4e}  '
                    f'adiff={adiff:.2e}  rdiff={rdiff:.2e}'
                )

            if not summary_only and not quiet:
                flag = '' if ok else ' <--'
                print(f'{freq_str:10s}  {n1:18s}  {abs(z1):14.4e}  {abs(z2):14.4e}  '
                      f'{adiff:12.2e}  {rdiff:10.2e}  {status}{flag}')

    total = passes + fails
    pct = 100 * passes / total if total else 0

    if not quiet:
        print()
        if only_ref:
            print(f'Only in {ref_file}: {", ".join(only_ref)}')
        if only_tst:
            print(f'Only in {test_file}: {", ".join(only_tst)}')
        freq_matched = len(freq_pairs)
        freq_ref = len(ref_freqs)
        freq_tst = len(tst_freqs)
        print(f'Harmonics: {freq_ref} ref, {freq_tst} test, {freq_matched} matched')
        print(f'Result: {passes}/{total} PASS ({pct:.0f}%)  '
              f'[atol={atol:.0e}  rtol={rtol:.0e}{"  mag-only" if mag_only else ""}]')
        if fails:
            print(f'\nFailed values ({fails}):')
            for d in fail_details:
                print(d)

    return fails == 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('ref', help='Reference .raw file')
    p.add_argument('test', help='Test .raw file')
    p.add_argument('--atol', type=float, default=1e-9,
                   help='Absolute tolerance (default: 1e-9)')
    p.add_argument('--rtol', type=float, default=1e-3,
                   help='Relative tolerance (default: 1e-3)')
    p.add_argument('--mag', action='store_true',
                   help='Compare magnitudes only')
    p.add_argument('--vars', default=None,
                   help='Comma-separated variable names to compare')
    p.add_argument('--summary', action='store_true',
                   help='Print summary line only')
    p.add_argument('--quiet', action='store_true',
                   help='Minimal output; exit 0=pass 1=fail')
    args = p.parse_args()

    ok = compare(args.ref, args.test,
                 atol=args.atol, rtol=args.rtol,
                 mag_only=args.mag,
                 filter_vars=args.vars,
                 summary_only=args.summary,
                 quiet=args.quiet)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
