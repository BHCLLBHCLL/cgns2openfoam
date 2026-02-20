# OpenFOAM Sliding Mesh LES Test Case

## 概述
这是一个完整的 OpenFOAM 瞬态 LES 滑移网格算例，由 cgns2openfoam 工具从 CGNS 网格转换生成。

## 配置参数

### 流体属性 (constant/transportProperties)
- 流体：空气 @ 20°C, 1 atm
- 运动粘度：ν = 1.5e-05 m²/s

### 湍流模型 (constant/turbulenceProperties)
- 模拟类型：LES
- LES 模型：Smagorinsky
- 亚格子尺度：cubeRootVol

### 动网格 (constant/dynamicMeshDict)
- 类型：solidBody 旋转运动
- 旋转区域：bladeZone
- 旋转中心：(0 0 0)
- 旋转轴：(0 0 1) - z轴
- 角速度：100 rad/s

### 网格 (constant/polyMesh)
- 总点数：83,450
- 总单元数：415,854
- 内部面：798,664
- 边界面：66,088

### 边界条件

| 边界 | 类型 | p | U | nut |
|------|------|---|---|-----|
| Interface_1_airZone | cyclicAMI | cyclicAMI | cyclicAMI | cyclicAMI |
| Interface_1_bladeZone | cyclicAMI | cyclicAMI | cyclicAMI | cyclicAMI |
| airZone_outlet | patch | fixedValue 0 | pressureInletOutletVelocity | calculated |
| bladeZone_blade | wall | zeroGradient | movingWallVelocity | nutUSpaldingWallFunction |

### 求解器 (system/controlDict)
- 求解器：pimpleFoam
- 时间步长：自适应，maxCo = 0.5
- 初始 deltaT：1e-5 s
- 最大 deltaT：1e-4 s
- 结束时间：0.1 s

## 运行方法

### 串行运行
```bash
source $FOAM_BASHRC
cd testCase
pimpleFoam
```

### 并行运行
```bash
source $FOAM_BASHRC
cd testCase
./Allrun
```

### 清理
```bash
./Allclean
```

## 后处理

### ParaView 查看
```bash
paraview testCase/test.foam
```

### 检查网格质量
```bash
checkMesh
```

## Cell Zones

| Zone | 单元数 | 描述 |
|------|--------|------|
| airZone | 342,533 | 静止区域 |
| bladeZone | 73,321 | 旋转区域 (100 rad/s) |

## 注意事项

1. 确保已正确设置 OpenFOAM 环境 (`source $FOAM_BASHRC`)
2. cyclicAMI 接口需要足够的重叠区域才能正确插值
3. 首次运行建议先检查网格质量 (`checkMesh`)
4. 如需调整转速，修改 `constant/dynamicMeshDict` 中的 `omega` 值
