#!/bin/bash
# Generate mesh.cgns for all cases 01-12 using the C++ writer (compatible with cgns2openfoam).
# Run from project root. Requires build/write_cgns_cases and LD_LIBRARY_PATH for CGNS/HDF5.
#
# Usage: ./generate_cgns_with_cpp.sh [path_to_write_cgns_cases]
#   Default: ./write_cgns_cases (e.g. run from build dir: ../cgns_test_lib/scripts/generate_cgns_with_cpp.sh ./write_cgns_cases)

WRITER="${1:-./write_cgns_cases}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TEST_LIB_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ ! -x "$WRITER" ] && ! command -v "$WRITER" &>/dev/null; then
  echo "Writer not found or not executable: $WRITER" >&2
  echo "Build first: cd build && make write_cgns_cases" >&2
  echo "Usage: $0 [path_to_write_cgns_cases]" >&2
  exit 1
fi

for i in $(seq 1 12); do
  case_num=$(printf "%02d" $i)
  dir="$TEST_LIB_ROOT/${case_num}_"
  # Match directory name by number (01_single_zone_tetra, 02_single_zone_hex, ...)
  for d in "$TEST_LIB_ROOT"/${case_num}_*; do
    [ -d "$d" ] || continue
    out="$d/mesh.cgns"
    echo "Generating case $case_num -> $out"
    "$WRITER" "$i" "$out" || exit 1
    break
  done
done
echo "All 12 cases generated."
