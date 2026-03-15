# CGNS 测试库 (cgns_test_lib)

本目录包含 10+ 个 CGNS 测试用例，用于验证 **cgns2openfoam** 的转换功能。用例覆盖从最简单的 3D 不可压外流场，到电子产品温度场分析，再到带旋转部件（interface/cyclicAMI）的风扇瞬态仿真。

## 用例列表

| 编号 | 目录名 | 描述 | 类型 |
|------|--------|------|------|
| 01 | `01_single_zone_tetra` | 单区域四面体网格，最小示例 | 不可压流 |
| 02 | `02_single_zone_hex` | 单区域六面体网格 | 不可压流 |
| 03 | `03_external_flow_3d` | 3D 外流场（远场/对称/壁面） | 不可压外流 |
| 04 | `04_pipe_flow` | 圆管/方管流动（入口/出口/壁面） | 不可压内流 |
| 05 | `05_duct_flow` | 矩形通道流动 | 不可压内流 |
| 06 | `06_electronics_solid` | 纯固体区域，电子产品温度场 | 热传导 |
| 07 | `07_electronics_fluid_solid` | 流体+固体两区域（共轭传热） | 共轭传热 |
| 08 | `08_heat_sink` | 简化翅片散热器 | 温度场 |
| 09 | `09_two_zone_interface` | 两区域 1-to-1 对接界面（无旋转） | 多区域 |
| 10 | `10_fan_rotor_stator` | 风扇转子/静子，带 interface | 旋转机械 |
| 11 | `11_fan_transient_interface` | 风扇瞬态，interface/cyclicAMI | 瞬态旋转 |
| 12 | `12_multi_zone_mixed` | 多区域混合单元类型 | 综合 |

## 目录结构

每个用例子目录包含：

- `mesh.cgns` — CGNS 网格文件（可由生成脚本生成或预置）
- `README.md` — 用例说明、边界条件、推荐求解设置

```
cgns_test_lib/
├── README.md                 # 本文件
├── 01_single_zone_tetra/
│   ├── mesh.cgns
│   └── README.md
├── 02_single_zone_hex/
│   └── ...
└── ...
```

## 生成网格

转换器使用 CGNS 4.x C 库读取，要求 CGNS 文件为 C API 写出的 HDF5 布局。**推荐用仓库内 C++ 工具一次性生成全部 12 个用例的兼容网格**。

### 方式一：C++ 工具生成全部 12 例（推荐）

构建后从 `build` 目录执行：

```bash
cd build
export LD_LIBRARY_PATH="/home/sdcll221/programs/CGNS-4.2.0/lib:/home/sdcll221/programs/hdf5-2.0.0/lib:$LD_LIBRARY_PATH"

# 一次性生成 01～12 的 mesh.cgns 到各用例子目录
../cgns_test_lib/scripts/generate_cgns_with_cpp.sh ./write_cgns_cases
```

或按用例单独生成：

```bash
./write_cgns_cases 1 ../cgns_test_lib/01_single_zone_tetra/mesh.cgns
./write_cgns_cases 2 ../cgns_test_lib/02_single_zone_hex/mesh.cgns
# ... 3～12 同理
```

### 方式二：仅 01/02 的简易工具

```bash
./write_cgns_minimal tetra ../cgns_test_lib/01_single_zone_tetra/mesh.cgns
./write_cgns_minimal hex   ../cgns_test_lib/02_single_zone_hex/mesh.cgns
```

### 方式三：Python + meshio（可选）

meshio 写出的 CGNS 与当前 CGNS 4.2 C 库的 HDF5 布局可能不兼容，转换时可能报错，仅作几何参考时可使用。

## 一键生成并验证

在项目根目录已构建（`build` 下存在 `cgns2openfoam` 和 `write_cgns_cases`）时，可直接执行：

```bash
./cgns_test_lib/scripts/generate_and_verify.sh
```

该脚本会依次：生成 12 个 `mesh.cgns` → 对每个用例运行 `cgns2openfoam` 转换 → 输出通过/失败统计。转换结果目录在 `build/out_01_single_zone_tetra`、`build/out_02_single_zone_hex` 等。若 CGNS/HDF5 不在默认路径，可先设置环境变量：

```bash
export CGNS_ROOT=/path/to/CGNS-4.2.0
export HDF5_ROOT=/path/to/hdf5-2.0.0
./cgns_test_lib/scripts/generate_and_verify.sh
```

## 运行转换（手动）

在项目根目录构建并进入 `build` 后，对每个用例执行：

```bash
# 假设在项目根目录，已执行 mkdir build && cd build && cmake .. && make
./cgns2openfoam ../cgns_test_lib/01_single_zone_tetra/mesh.cgns ./out_01
./cgns2openfoam ../cgns_test_lib/02_single_zone_hex/mesh.cgns ./out_02
# ...
```

或使用提供的批量脚本：

```bash
cd build
export LD_LIBRARY_PATH="/home/sdcll221/programs/CGNS-4.2.0/lib:/home/sdcll221/programs/hdf5-2.0.0/lib:$LD_LIBRARY_PATH"
../cgns_test_lib/scripts/run_all_conversions.sh ./cgns2openfoam
```

## 生成可执行的 OpenFOAM 算例（01–10）

转换完成后，可为前 10 个用例生成完整可运行的 OpenFOAM 算例（system/、0/、constant/ 等），便于直接运行 simpleFoam：

```bash
# 1) 先完成转换（见“一键生成并验证”或 run_all_conversions.sh）
# 2) 在项目根目录执行
python3 cgns_test_lib/scripts/create_openfoam_cases.py build
```

生成后，每个 `build/out_01_single_zone_tetra` … `build/out_10_fan_rotor_stator` 均为完整算例，可进入目录运行：

```bash
source /path/to/OpenFOAM-v2412/etc/bashrc
cd build/out_01_single_zone_tetra
./Allrun
# 或直接: simpleFoam
```

批量验证 10 个算例（需已 source OpenFOAM）：

```bash
bash -c 'source /path/to/OpenFOAM-v2412/etc/bashrc; cd build; for d in out_0{1..5}_* out_0{6..9}_* out_10_*; do (cd "$d" && simpleFoam >/dev/null 2>&1) && echo "$d PASS" || echo "$d FAIL"; done'
```

若 CGNS/HDF5 安装在其他路径，可设置：

```bash
export LD_LIBRARY_PATH="$CGNS_ROOT/lib:$HDF5_ROOT/lib:$LD_LIBRARY_PATH"
```

## CGNS 格式要求（与 cgns2openfoam 一致）

- 基（Base）下为非结构化区域（Unstructured）
- 坐标：`CoordinateX`、`CoordinateY`、`CoordinateZ`，类型 `RealDouble`
- 体单元：TETRA_4、PYRA_5、PENTA_6、HEXA_8 等
- 边界面：TRI_3、QUAD_4，Section 名称可作为边界 patch 名
- 边界类型（ZoneBC）：BCWall、BCInflow、BCOutflow、BCSymmetryPlane、BCFarfield 等，会映射为 OpenFOAM 的 wall、inlet、outlet、symmetryPlane、freestream 等
- 多区域界面：1-to-1 或 General 连接，可带周期/旋转信息（用于 cyclicAMI）

## 参考

- 转换器用法见项目根目录 `README.md`
- 设计细节见 `DEVELOPMENT.md`
