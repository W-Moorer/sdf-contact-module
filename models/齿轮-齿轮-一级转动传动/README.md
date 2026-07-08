# 齿轮-齿轮-一级转动传动

本目录保存 RecurDyn 一级齿轮传动参考算例。NexDyn 对齐必须从本目录下的 RMD 建模抽取几何和建模参数，不再使用旧的手工 OBJ 作为粗网格或细网格来源。

## 工况

| 目录 | 工况 | RMD 几何 | NexDyn 当前状态 |
| --- | --- | --- | --- |
| `First-stage gear transmission_01` | 粗网格，无摩擦 | GGEOM 2/4，每个齿轮 3212 个三角面 | 历史 0.005 s smoke 通过；原始 1.0 s 回归仍需稳定 |
| `First-stage gear transmission_02` | 细网格，无摩擦 | GGEOM 2/4，约 65148/65138 个三角面 | 原始 `1.0 s / 0.001 s` smoke 已恢复并加入 passing 回归；strict 理论角速度仍待收紧 |
| `First-stage gear transmission_03` | 细网格，加摩擦 | 与 `_02` 同级细网格，RMD 含摩擦系数 | 已可抽取几何，NexDyn 摩擦对齐待后续严格验收 |

RecurDyn 求解设置来自 RMD/MSG：

- RecurDyn 版本：2023 Windows x64, `10.1.11516.0`
- 分析：Dynamic Analysis
- 终止时间：`1.0 s`
- 输出步数：`1000`
- 最大步长：`1.0e-3`
- 初始步长：`1.0e-6`
- 误差容限：`5.0e-3`
- 单位：`kg, mm, N, s`

## NexDyn 建模映射

- 主动轮：RMD `PART / 3`, `Verbindungsrad_Rechts6`, GGEOM `2`。
- 从动轮：RMD `PART / 4`, `Verbindungsrad_Links6`, GGEOM `4`。
- 接触对：RMD `IGGEOMID = 4` 为 action/left，`JGGEOMID = 2` 为 base/right。
- 转动副标记：使用 RMD 中 `Ground.Marker1/2` 与齿轮 `Marker1` 的位置和 Euler 角。
- 驱动：RMD 中 `ROTATION VELOCITY FUNCTION = 2*PI`。NexDyn 导入的 B3 关节坐标与 RecurDyn B3 正方向相反，因此例子内部用 `-2*pi` 驱动来复现 RMD 全局角速度；输出中的 `drive_speed` 保持 RMD 源值 `+2*pi`。
- 接触材料：RMD 的 `BPEN=1e-2`, `KORDER=2`, `K=100000`, `C=10` 固化在例子中。NexDyn 沿用 RMD 的 mm 坐标和 kg 质量，`K/C` 从 N 转为内部 `kg*mm/s^2` 时乘以 `1000`。这是单位映射，不是调参。

## RMD 几何抽取

```powershell
python .\scripts\extract_recurdyn_gear_surfaces.py
```

输出到：

```text
data/input/obj_library/generated_from_rmd/
```

该目录是 RMD 派生物，已加入 `.gitignore`。验收或 sweep 脚本会自动先执行抽取，因此 Git 中不提交生成 OBJ。

## NexDyn 入口

```text
example/02_Gear_System/DirectGeometryTriMesh/Example_GearCoarseAlignment.cpp
example/02_Gear_System/DirectGeometryTriMesh/Example_GearFineAlignment.cpp
```

构建：

```powershell
cmake --build build --config Release --target Example_GearCoarseAlignment Example_GearFineAlignment -- /m
```

单次运行：

```powershell
$env:NEXDYN_CONTACT_MANAGER_MODE = "ManagedPenalty"
$env:NEXDYN_GEAR_ALIGNMENT_END_TIME = "1.0"
$env:NEXDYN_GEAR_ALIGNMENT_OUTPUT_STEP = "0.001"
.\build\bin\example\Release\Example_GearCoarseAlignment.exe
.\build\bin\example\Release\Example_GearFineAlignment.exe
```

批量后端 sweep：

```powershell
.\scripts\run_gear_coarse_alignment_backends.ps1 -MeshLevel coarse -EndTime 1.0 -OutputStep 0.001
.\scripts\run_gear_coarse_alignment_backends.ps1 -MeshLevel fine -EndTime 1.0 -OutputStep 0.001
```

验收脚本只设置以下运行环境变量：

- `NEXDYN_CONTACT_MANAGER_MODE`
- `NEXDYN_GEAR_ALIGNMENT_END_TIME`
- `NEXDYN_GEAR_ALIGNMENT_OUTPUT_STEP`

接触材料 override 环境变量仅保留为源代码注释中的 debug-only 提示，不接入验收脚本。

## 历史 0.005 s 诊断结果

粗网格 `_01`：

| 后端 | 从动轮轴向角速度 | 目标 | 误差 | smoke | strict |
| --- | ---: | ---: | ---: | --- | --- |
| `ManagedPenalty` | `-5.7921863074` | `-6.2831853072` | `0.4909989998` | PASS | FAIL |
| `ManagedAugmentedLagrangian` | `-5.7926785618` | `-6.2831853072` | `0.4905067454` | PASS | FAIL |
| `ManagedKKT` | `-5.7921863074` | `-6.2831853072` | `0.4909989998` | PASS | FAIL |
| `ManagedFrictionlessLCP` | `-5.6733721525` | `-6.2831853072` | `0.6098131546` | PASS | FAIL |
| `ManagedFrictionPyramid` | `-5.6733721525` | `-6.2831853072` | `0.6098131546` | PASS | FAIL |
| `ManagedSOCCP` | `-5.6733721525` | `-6.2831853072` | `0.6098131546` | PASS | FAIL |

细网格 `_02`：

| 后端 | 从动轮轴向角速度 | 目标 | 误差 | smoke | strict |
| --- | ---: | ---: | ---: | --- | --- |
| `ManagedPenalty` | `-5.7920326739` | `-6.2831853072` | `0.4911526332` | PASS | FAIL |
| `ManagedAugmentedLagrangian` | `-5.7925242293` | `-6.2831853072` | `0.4906610779` | PASS | FAIL |
| `ManagedKKT` | `-5.7920326739` | `-6.2831853072` | `0.4911526332` | PASS | FAIL |
| `ManagedFrictionlessLCP` | `-5.6733059339` | `-6.2831853072` | `0.6098793733` | PASS | FAIL |
| `ManagedFrictionPyramid` | `-5.6733059339` | `-6.2831853072` | `0.6098793733` | PASS | FAIL |
| `ManagedSOCCP` | `-5.6733059339` | `-6.2831853072` | `0.6098793733` | PASS | FAIL |

当前容差分两级：

- 齿轮回归必须使用 RMD 原始时长和输出步长：`1.0 s / 0.001 s`；当前不再用短时 smoke 作为 passing regression。
- 细网格 `_02` 当前 passing smoke 使用 `ManagedPenalty`，终值从动轮轴向角速度为 `-6.4881172773 rad/s`，理论值为 `-6.2831853072 rad/s`，误差为 `-0.2049319701 rad/s`，受 `0.25 rad/s` 容差保护。
- 历史 0.005 s smoke 结果只保留为诊断记录，不作为 passing regression 的时间设置。
- `0.1 rad/s` strict 终值容差和完整 1 s 曲线误差仍是后续精度修复目标。

## 待办

- 用 `dt=0.001`, `t=1.0 s` 输出粗网格与 RecurDyn 的完整角速度曲线对比并修复其稳定性。
- 完成 `_03` 摩擦接触映射，加入切向响应和关节反力对比。
- 将严格容差扩展到完整 1 s 曲线误差。
