# 10 - 风扇转子/静子（interface）

风扇类几何：转子区 + 静子区，带 interface。

- **用途**：旋转机械、定常/准稳态。
- **边界**：入口、出口、壁面；转子-静子 interface（cyclicAMI 或 sliding mesh 的 interface）。
- **说明**：若 CGNS 中该连接带周期/旋转信息，转换器会写出 cyclicAMI 的 neighbourPatch、rotationAxis 等。
