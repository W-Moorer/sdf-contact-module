# First-stage gear transmission_02

细网格、无摩擦 RecurDyn 基准工况。该目录中的 `First-stage gear transmission.rmd` 是 NexDyn 细网格建模的唯一来源。

## RMD 摘要

- 模型：`Model2`
- 主动轮：`PART / 3`，`Verbindungsrad_Rechts6`，GGEOM `2`
- 从动轮：`PART / 4`，`Verbindungsrad_Links6`，GGEOM `4`
- 网格：主动轮 `32574` 节点、`65148` 三角面；从动轮 `32569` 节点、`65138` 三角面
- 接触：`GGEOMCONTACT / 1`，无摩擦，`K=100000`，`C=10`，`KORDER=2`，`BPEN=1e-2`
- 驱动：`RevJoint1.RMotion`，`ROTATION VELOCITY FUNCTION = 2*PI`
- 动力学：`END = 1.0`，`STEP = 1000`
- 单位：`kg, mm, N, s`

## NexDyn 对齐

细网格入口复用粗网格实现，通过宏切换到 `_02` 的 RMD 抽取结果：

```powershell
cmake --build build --config Release --target Example_GearFineAlignment -- /m
.\scripts\run_gear_coarse_alignment_backends.ps1 -MeshLevel fine -EndTime 0.005 -OutputStep 0.005
```

最新 0.005 s sweep 结果：7/8 后端通过 smoke，`ManagedKKT` 未通过。`ManagedPenalty` 的从动轮轴向角速度为 `-5.7922317312 rad/s`，相对理论 `-2*pi` 的误差为 `0.4909535760 rad/s`。

## 说明

该工况当前用于网格来源和后端消费链路验证。严格容差 `0.1 rad/s` 尚未达成；后续需要解析 RecurDyn 结果通道并修正 KKT/接触尺度后再作为严格细网格验收。
