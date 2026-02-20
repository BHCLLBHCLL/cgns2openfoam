# cgns2openfoam

将 CGNS 网格文件转换为 OpenFOAM 计算格式。

## 功能

- 读取 CGNS 非结构化网格
- 提取边界条件名称和类型
- 生成 OpenFOAM polyMesh 文件：
  - `points` - 网格点坐标
  - `faces` - 面定义
  - `owner` - 面的所属单元
  - `neighbour` - 内部面的邻居单元
  - `boundary` - 边界定义
  - `cellZones` - 单元区域
  - `faceZones` - 面区域

## 边界类型映射

| CGNS BC Type | OpenFOAM Type |
|--------------|---------------|
| BCWall, BCWallViscous | wall |
| BCInflow, BCInflowSubsonic | inlet |
| BCOutflow, BCOutflowSubsonic | outlet |
| BCSymmetryPlane | symmetryPlane |
| BCFarfield | freestream |
| 其他 | patch |

## 依赖

- CGNS 4.x（需支持 HDF5）
- HDF5 2.x
- CMake 3.14+
- C++17 编译器

## 构建

```bash
mkdir build && cd build
cmake .. -DCGNS_ROOT=/path/to/cgns -DHDF5_ROOT=/path/to/hdf5
make
```

或者设置环境变量后直接构建：

```bash
export CGNS_ROOT=/home/user/Programs/CGNS-4.5.1
export HDF5_ROOT=/home/user/Programs/hdf5-2.0.0
mkdir build && cd build
cmake ..
make
```

## 用法

```bash
./cgns2openfoam <input.cgns> <output_case_dir>
```

示例：

```bash
./cgns2openfoam mesh.cgns ./myCase
```

输出：

```
myCase/
└── constant/
    └── polyMesh/
        ├── points
        ├── faces
        ├── owner
        ├── neighbour
        ├── boundary
        ├── cellZones
        └── faceZones
```

## 项目结构

```
├── CMakeLists.txt
├── include/
│   ├── cgns_reader.h
│   └── openfoam_writer.h
├── src/
│   ├── main.cpp
│   ├── cgns_reader.cpp
│   └── openfoam_writer.cpp
└── README.md
```

## 限制

- 仅支持非结构化网格
- 支持的单元类型：TETRA, PYRA, PENTA, HEXA
- 支持的面类型：TRI, QUAD
