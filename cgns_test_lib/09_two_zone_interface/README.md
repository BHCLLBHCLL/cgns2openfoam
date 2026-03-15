# 09 - 两区域界面（无旋转）

两个区域通过共享界面连接，1-to-1 或 General 连接，无周期/旋转。

- **用途**：多区域、interface 对接。
- **边界**：两区域各自外边界 + 中间 interface（在 CGNS 中为 ZoneGridConnectivity）。
- **说明**：转换后 interface 对应 OpenFOAM 的 patch 对或 cyclicAMI（若需旋转则用 10/11）。
