# 04 - 管流

管道状流动（此处为方管简化）。

- **用途**：不可压内流、入口/出口/壁面。
- **边界**：inlet、outlet、wall（若 CGNS 中已设 ZoneBC）。
- **建议**：转换后可用 simpleFoam 或 pimpleFoam，设置 inlet 速度、outlet 压力。
