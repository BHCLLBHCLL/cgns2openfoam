# 03 - 3D 外流场

三维外流场用网格（如绕流/远场）。

- **用途**：不可压外流、远场/对称/壁面边界。
- **边界**：若 CGNS 中设置 ZoneBC，可为 BCFarfield、BCSymmetryPlane、BCWall 等，对应 OpenFOAM 的 freestream、symmetryPlane、wall。
- **建议**：转换后在 OpenFOAM 中设置 0/U、0/p 及 fvSchemes 做简单外流测试。
