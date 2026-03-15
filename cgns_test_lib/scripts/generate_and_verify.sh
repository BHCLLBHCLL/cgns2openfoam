#!/bin/bash
# 一键生成 12 个 CGNS 测试网格并运行 cgns2openfoam 转换验证。
#
# 用法（在项目根目录执行）:
#   ./cgns_test_lib/scripts/generate_and_verify.sh
#
# 或指定 build 目录（若不在项目根执行）:
#   ./cgns_test_lib/scripts/generate_and_verify.sh /path/to/build
#
# 依赖: 已构建 cgns2openfoam、write_cgns_cases；CGNS/HDF5 库路径见下方。

set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="${1:-$REPO_ROOT/build}"
CGNS_ROOT="${CGNS_ROOT:-/home/sdcll221/programs/CGNS-4.2.0}"
HDF5_ROOT="${HDF5_ROOT:-/home/sdcll221/programs/hdf5-2.0.0}"
export LD_LIBRARY_PATH="${CGNS_ROOT}/lib:${HDF5_ROOT}/lib:${LD_LIBRARY_PATH}"

echo "=========================================="
echo "  CGNS 测试库 - 一键生成并验证"
echo "=========================================="
echo "  BUILD_DIR: $BUILD_DIR"
echo "  CGNS:      $CGNS_ROOT"
echo "  HDF5:      $HDF5_ROOT"
echo "=========================================="

if [ ! -d "$BUILD_DIR" ]; then
  echo "错误: build 目录不存在: $BUILD_DIR" >&2
  echo "请先执行: mkdir build && cd build && cmake .. -DCGNS_ROOT=$CGNS_ROOT -DHDF5_ROOT=$HDF5_ROOT && make" >&2
  exit 1
fi

CONVERTER="$BUILD_DIR/cgns2openfoam"
WRITER="$BUILD_DIR/write_cgns_cases"

for exe in "$CONVERTER" "$WRITER"; do
  if [ ! -x "$exe" ]; then
    echo "错误: 未找到可执行文件: $exe" >&2
    echo "请在 build 目录执行: make" >&2
    exit 1
  fi
done

echo ""
echo "[1/2] 生成 12 个 CGNS 网格..."
"$SCRIPT_DIR/generate_cgns_with_cpp.sh" "$WRITER" || exit 1

echo ""
echo "[2/2] 运行 cgns2openfoam 转换验证..."
cd "$BUILD_DIR"
"$SCRIPT_DIR/run_all_conversions.sh" "$CONVERTER"
EXIT_CODE=$?

echo ""
if [ "$EXIT_CODE" -eq 0 ]; then
  echo "=========================================="
  echo "  全部完成: 12 个用例已生成并转换通过"
  echo "=========================================="
else
  echo "=========================================="
  echo "  部分失败，请查看上方输出"
  echo "=========================================="
fi
exit $EXIT_CODE
