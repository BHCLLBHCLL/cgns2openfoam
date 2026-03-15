# 11 - 风扇瞬态（interface）

与 10 同类几何，用于瞬态仿真（带旋转的 interface）。

- **用途**：瞬态旋转、风扇/泵的瞬态流场。
- **边界**：同 10；瞬态需在 controlDict 中设置 time step 与 endTime。
- **动网格**：已生成 `constant/dynamicMeshDict`，指定 **Rotor** 区绕 z 轴旋转（`omega` 可调）。运行瞬态动网格请使用 **pimpleDyMFoam**（需在 system/controlDict 中把 application 改为 pimpleDyMFoam 并设好 deltaT/endTime）。
- **建议**：转换并生成算例后，用 pimpleDyMFoam 进行瞬态计算。
