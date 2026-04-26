#!/usr/bin/env bash
# run_test.sh - Convert, simulate, and optionally compare a HBFree test circuit.
#
# Usage:
#   run_test.sh <circuit.ckt> [options]
#
# Options:
#   --ref <file.raw>    compare output against this reference .raw
#   --atol <val>        absolute tolerance for cmpraw (default: 1e-9)
#   --rtol <val>        relative tolerance for cmpraw (default: 1e-3)
#   --keep              keep intermediate .hbl/.nodes/.elems files
#   --log <file>        save simulator stdout to file (default: <circuit>.log)
#
# Exit codes: 0=pass, 1=sim-error, 2=cmpraw-fail

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

S2H="$REPO_DIR/s2h"
HBL="$REPO_DIR/hbl"
CMPRAW="$SCRIPT_DIR/cmpraw.py"
S2H_CFG="$REPO_DIR/spice2hbl/s2h.cfg"

# ---- parse arguments -------------------------------------------------------
CKT=""
REF=""
ATOL="1e-9"
RTOL="1e-3"
KEEP=0
LOGFILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref)   REF="$2";   shift 2 ;;
        --atol)  ATOL="$2";  shift 2 ;;
        --rtol)  RTOL="$2";  shift 2 ;;
        --keep)  KEEP=1;     shift   ;;
        --log)   LOGFILE="$2"; shift 2 ;;
        -*)      echo "Unknown option: $1" >&2; exit 1 ;;
        *)
            if [[ -z "$CKT" ]]; then CKT="$1"
            else echo "Unexpected argument: $1" >&2; exit 1
            fi
            shift ;;
    esac
done

if [[ -z "$CKT" ]]; then
    echo "Usage: run_test.sh <circuit.ckt> [--ref ref.raw] [--atol N] [--rtol N] [--keep] [--log file]"
    exit 1
fi

if [[ ! -f "$CKT" ]]; then
    echo "Error: circuit file not found: $CKT" >&2
    exit 1
fi

CKT_ABS="$(realpath "$CKT")"
CKT_DIR="$(dirname "$CKT_ABS")"
CKT_BASE="$(basename "$CKT_ABS")"

# Derived file names (all in same directory as .ckt)
STEM="${CKT_BASE%.ckt}"          # e.g. di1  (strip .ckt if present)
HBL_FILE="$CKT_DIR/${CKT_BASE}.hbl"
NODES_FILE="$CKT_DIR/${CKT_BASE}.nodes"
ELEMS_FILE="$CKT_DIR/${CKT_BASE}.elems"
RAW_FILE="$CKT_DIR/${STEM}.raw"
[[ -z "$LOGFILE" ]] && LOGFILE="$CKT_DIR/${CKT_BASE}.run.log"

# ---- banner ----------------------------------------------------------------
echo "=== run_test: $CKT_BASE ==="
echo "  hbl    : $HBL"
echo "  s2h    : $S2H"
echo "  cfg    : $S2H_CFG"

# s2h treats '/' as a flag prefix (Windows compat), so all positional
# path arguments must be relative.  Run everything from CKT_DIR.
pushd "$CKT_DIR" > /dev/null

HBL_REL="${CKT_BASE}.hbl"
NODES_REL="${CKT_BASE}.nodes"
ELEMS_REL="${CKT_BASE}.elems"

# ---- step 1: s2h conversion ------------------------------------------------
echo ""
echo "--- [1/3] s2h: converting to HBL ---"
set +e
"$S2H" -c"$S2H_CFG" -n"$NODES_REL" -m"$ELEMS_REL" "$CKT_BASE" "$HBL_REL" 2>&1
S2H_RC=$?
set -e

if [[ $S2H_RC -gt 1 ]]; then
    echo "ERROR: s2h failed (exit $S2H_RC)" >&2
    exit 1
fi
echo "  -> $HBL_FILE"

# ---- step 2: hbl simulation ------------------------------------------------
echo ""
echo "--- [2/3] hbl: running simulation ---"
set +e
"$HBL" "$HBL_REL" 2>&1 | tee "$LOGFILE"
HBL_RC=${PIPESTATUS[0]}
set -e

if [[ $HBL_RC -ne 0 ]]; then
    echo "ERROR: hbl simulation failed (exit $HBL_RC)" >&2
    popd > /dev/null
    exit 1
fi

RAW_REL="${STEM}.raw"
if [[ ! -f "$RAW_REL" ]]; then
    echo "ERROR: expected output not found: $RAW_FILE" >&2
    popd > /dev/null
    exit 1
fi
echo "  -> $RAW_FILE"

# ---- step 3: optional comparison -------------------------------------------
RESULT=0
if [[ -n "$REF" ]]; then
    if [[ ! -f "$REF" ]]; then
        echo "WARNING: reference not found: $REF — skipping comparison" >&2
    else
        echo ""
        echo "--- [3/3] cmpraw: comparing against reference ---"
        echo "  ref  : $REF"
        echo "  test : $RAW_FILE"
        set +e
        python3 "$CMPRAW" --atol "$ATOL" --rtol "$RTOL" "$REF" "$RAW_FILE"
        RESULT=$?
        set -e
        if [[ $RESULT -eq 0 ]]; then
            echo "  PASS"
        else
            echo "  FAIL"
        fi
    fi
fi

# ---- cleanup ----------------------------------------------------------------
if [[ $KEEP -eq 0 ]]; then
    rm -f "$HBL_REL" "$NODES_REL" "$ELEMS_REL"
fi

popd > /dev/null

echo ""
echo "=== done: $CKT_BASE ==="
[[ -n "$REF" ]] && echo "  result : $([ $RESULT -eq 0 ] && echo PASS || echo FAIL)"
echo "  raw    : $RAW_FILE"
echo "  log    : $LOGFILE"

[[ -n "$REF" ]] && exit $RESULT || exit 0
