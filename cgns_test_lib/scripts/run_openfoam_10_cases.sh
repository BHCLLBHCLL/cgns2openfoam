#!/bin/bash
# 在 10 个 OpenFOAM 案例中依次运行 simpleFoam，并输出验证总结。
# 需已执行: generate_and_verify.sh 且 create_openfoam_cases.py 已生成算例文件。
#
# 用法: ./run_openfoam_10_cases.sh [build_dir]
#   build_dir 默认 ../.. 的 build（即项目 build）
# 需先 source OpenFOAM: export OPENFOAM_BASHRC=/path/to/OpenFOAM-v2412/etc/bashrc

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${1:-$(cd "$SCRIPT_DIR/../.." && pwd)/build}"
CASES=( out_01_single_zone_tetra out_02_single_zone_hex out_03_external_flow_3d out_04_pipe_flow out_05_duct_flow \
        out_06_electronics_solid out_07_electronics_fluid_solid out_08_heat_sink out_09_two_zone_interface out_10_fan_rotor_stator )
PASS=0
FAIL=0

if [ -n "${OPENFOAM_BASHRC}" ]; then
  source "${OPENFOAM_BASHRC}" 2>/dev/null || true
fi

if ! command -v simpleFoam &>/dev/null; then
  echo "simpleFoam not found. Source OpenFOAM first, e.g.: export OPENFOAM_BASHRC=/path/to/OpenFOAM-v2412/etc/bashrc" >&2
  exit 1
fi

echo "=========================================="
echo "  10 个 OpenFOAM 案例验证"
echo "  BUILD_DIR: $BUILD_DIR"
echo "=========================================="

for name in "${CASES[@]}"; do
  dir="$BUILD_DIR/$name"
  if [ ! -d "$dir" ] || [ ! -f "$dir/system/controlDict" ]; then
    echo "[SKIP] $name (case not found or incomplete)"
    continue
  fi
  echo -n "Running $name ... "
  if ( cd "$dir" && simpleFoam > /dev/null 2>&1 ); then
    echo "[PASS]"
    PASS=$((PASS+1))
  else
    echo "[FAIL]"
    FAIL=$((FAIL+1))
  fi
done

echo "=========================================="
echo "  验证总结: 通过 $PASS, 失败 $FAIL / ${#CASES[@]}"
echo "=========================================="
[ "$FAIL" -eq 0 ]
