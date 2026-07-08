# First-stage gear transmission_01

粗网格、无摩擦 RecurDyn 基准工况。该目录中的 `First-stage gear transmission.rmd` 是 NexDyn 粗网格建模的唯一来源。

## RMD 摘要

- 模型：`Model2`
- 主动轮：`PART / 3`，`Verbindungsrad_Rechts6`，GGEOM `2`
- 从动轮：`PART / 4`，`Verbindungsrad_Links6`，GGEOM `4`
- 网格：每个齿轮 `1606` 节点、`3212` 三角面
- 接触：`GGEOMCONTACT / 1`，无摩擦，`K=100000`，`C=10`，`KORDER=2`，`BPEN=1e-2`
- 驱动：`RevJoint1.RMotion`，`ROTATION VELOCITY FUNCTION = 2*PI`
- 动力学：`END = 1.0`，`STEP = 1000`
- 单位：`kg, mm, N, s`

## NexDyn 对齐

抽取脚本从本目录 RMD 的 GGEOM `2/4` 生成局部 OBJ，并在 NexDyn 中挂到 RMD 的 `RM` 标记位置：

```powershell
python .\scripts\extract_recurdyn_gear_surfaces.py
```

粗网格入口：

```powershell
cmake --build build --config Release --target Example_GearCoarseAlignment -- /m
.\scripts\run_gear_coarse_alignment_backends.ps1 -MeshLevel coarse -EndTime 0.005 -OutputStep 0.005
```

最新 0.005 s sweep 结果：7/8 后端通过 smoke，`ManagedKKT` 未通过。`ManagedPenalty` 的从动轮轴向角速度为 `-5.7921635834 rad/s`，相对理论 `-2*pi` 的误差为 `0.4910217237 rad/s`。

## 说明

RMD `.out` 初始条件中主动轮全局角速度为：

```text
Wx =  3.141592653585
Wy =  0
Wz = -5.441398092706
```

NexDyn 已按该方向完成驱动符号对齐。当前剩余误差主要来自接触传力尺度、短时 penalty 软接触和 KKT 后端尺度建模。
