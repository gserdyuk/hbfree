#!/usr/bin/env python3
"""
make_conv_summary.py — extract convergence summary from an HBFree .ckt.run.log.

Usage:
    python3 scripts/make_conv_summary.py <circuit.ckt.run.log> [-o output.conv]

If -o is omitted, prints to stdout.

Three recognised outcomes (matching HBFree log verbatim):
    CONVERGED   — "AT N -TH ITERATION CONVERGED TO"
    STEP_LIMIT  — "AT N -T ITERATION MAX STEP= ..."
    DIVERGED    — "AT N -TH ITERATION CAN NOT MAKE GOOD STEP"
"""

import re
import sys
import os
import argparse

# Terminal iteration line patterns
_RE_CONVERGED  = re.compile(r'AT\s+(\d+)\s+-T[H]?\s+ITERATION\s+CONVERGED\s+TO', re.I)
_RE_STEP_LIMIT = re.compile(
    r'AT\s+(\d+)\s+-T[H]?\s+ITERATION\s+MAX\s+STEP\s*=\s*([\d.E+\-]+)', re.I)
_RE_DIVERGED   = re.compile(r'AT\s+(\d+)\s+-T[H]?\s+ITERATION\s+CAN\s+NOT\s+MAKE\s+GOOD\s+STEP', re.I)

# Follow-on value patterns (matched on lines immediately after the terminal line)
_RE_SOL_ERROR  = re.compile(r'SOLUTION\s+WITH\s+ERROR\s*<=\s*([\d.E+\-]+)', re.I)
_RE_FINAL_ERR  = re.compile(r'^\s*ERROR\s*=\s*([\d.E+\-]+)', re.I)   # final-only (no STEP= after)
_RE_NORM       = re.compile(r'1/2\s+OF\s+SQUARED\s+L-2\s+NORM\s+OF\s+ERROR\s*=\s*([\d.E+\-]+)', re.I)


def _lookahead(lines, start, n=5):
    return lines[start + 1 : start + 1 + n]


def parse_log(log_path):
    with open(log_path, 'r', errors='replace') as f:
        lines = f.readlines()

    result = {}

    for i, line in enumerate(lines):
        m = _RE_CONVERGED.search(line)
        if m:
            result['status'] = 'CONVERGED'
            result['iterations'] = int(m.group(1))
            for ahead in _lookahead(lines, i):
                if 'final_error' not in result:
                    me = _RE_SOL_ERROR.search(ahead)
                    if me:
                        result['final_error'] = me.group(1)
                if 'final_norm' not in result:
                    mn = _RE_NORM.search(ahead)
                    if mn:
                        result['final_norm'] = mn.group(1)
            continue

        m = _RE_STEP_LIMIT.search(line)
        if m:
            result['status'] = 'STEP_LIMIT'
            result['iterations'] = int(m.group(1))
            result['max_step'] = m.group(2)
            for ahead in _lookahead(lines, i):
                if 'final_error' not in result:
                    me = _RE_FINAL_ERR.match(ahead)
                    if me:
                        result['final_error'] = me.group(1)
                if 'final_norm' not in result:
                    mn = _RE_NORM.search(ahead)
                    if mn:
                        result['final_norm'] = mn.group(1)
            continue

        m = _RE_DIVERGED.search(line)
        if m:
            result['status'] = 'DIVERGED'
            result['iterations'] = int(m.group(1))
            for ahead in _lookahead(lines, i):
                if 'final_norm' not in result:
                    mn = _RE_NORM.search(ahead)
                    if mn:
                        result['final_norm'] = mn.group(1)
            continue

    if not result:
        result['status'] = 'UNKNOWN'

    return result


def format_conv(circuit_name, result):
    lines = [f'circuit: {circuit_name}\n']
    lines.append(f'status: {result["status"]}\n')
    if 'iterations' in result:
        lines.append(f'iterations: {result["iterations"]}\n')
    if 'max_step' in result:
        lines.append(f'max_step: {result["max_step"]}\n')
    if 'final_error' in result:
        lines.append(f'final_error: {result["final_error"]}\n')
    if 'final_norm' in result:
        lines.append(f'final_norm: {result["final_norm"]}\n')
    return ''.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='Extract HBFree convergence summary from a .ckt.run.log')
    parser.add_argument('log', help='.ckt.run.log file to parse')
    parser.add_argument('-o', '--out', default='-',
                        help='Output .conv file (default: stdout)')
    args = parser.parse_args()

    log_base = os.path.basename(args.log)
    circuit_name = re.sub(r'\.ckt\.run\.log$|\.run\.log$|\.log$', '', log_base)

    result = parse_log(args.log)
    text = format_conv(circuit_name, result)

    if args.out == '-':
        sys.stdout.write(text)
    else:
        with open(args.out, 'w') as f:
            f.write(text)
        print(f'written: {args.out}')


if __name__ == '__main__':
    main()
