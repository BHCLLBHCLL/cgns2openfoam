#!/bin/bash
# Run cgns2openfoam on all test cases. Creates out_01_single_zone_tetra, etc. in current directory.
#
# Usage:
#   ./run_all_conversions.sh [path_to_cgns2openfoam]
# If path_to_cgns2openfoam is omitted, uses ./cgns2openfoam (e.g. from build dir).
# Set LD_LIBRARY_PATH if needed: export LD_LIBRARY_PATH="$CGNS_ROOT/lib:$HDF5_ROOT/lib:$LD_LIBRARY_PATH"

CONVERTER="${1:-./cgns2openfoam}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TEST_LIB_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PASS=0
FAIL=0

if ! command -v "$CONVERTER" &>/dev/null && [ ! -x "$CONVERTER" ]; then
  echo "Converter not found or not executable: $CONVERTER" >&2
  echo "Usage: $0 [path_to_cgns2openfoam]" >&2
  echo "Example: $0 /path/to/build/cgns2openfoam" >&2
  exit 1
fi

for dir in "$TEST_LIB_ROOT"/0*_* "$TEST_LIB_ROOT"/1*_*; do
  [ -d "$dir" ] || continue
  case_dir="$(basename "$dir")"
  cgns="$dir/mesh.cgns"
  out_dir="out_$case_dir"
  if [ ! -f "$cgns" ]; then
    echo "Skip $case_dir: no mesh.cgns"
    continue
  fi
  echo "Converting $case_dir -> $out_dir"
  if "$CONVERTER" "$cgns" "$out_dir" 2>&1; then
    echo "[PASS] $case_dir"
    PASS=$((PASS+1))
  else
    echo "[FAIL] $case_dir"
    FAIL=$((FAIL+1))
  fi
done
echo "Done: $PASS passed, $FAIL failed."
[ "$FAIL" -eq 0 ]
