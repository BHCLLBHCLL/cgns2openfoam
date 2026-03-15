# 11 - 风扇瞬态（interface）

与 10 同类几何，用于瞬态仿真（带旋转的 interface）。

- **用途**：瞬态旋转、风扇/泵的瞬态流场。
- **边界**：同 10；瞬态需在 controlDict 中设置 time step 与 endTime，并配置 dynamicMeshDict 等。
- **建议**：转换后可用 icoDyMFoam / pimpleDyMFoam + AMI 或 overset 进行瞬态计算。
