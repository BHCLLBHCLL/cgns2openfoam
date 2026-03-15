# 01 - 单区域四面体

最小示例：单区域、纯四面体网格。

- **用途**：验证基本转换（单 zone、TETRA_4）。
- **边界**：由网格生成器得到的外表面，转换后可能为单一 patch（如 defaultFaces）或按 Section 命名。
- **建议**：转换后检查 `constant/polyMesh/points`、`faces`、`boundary` 即可。
