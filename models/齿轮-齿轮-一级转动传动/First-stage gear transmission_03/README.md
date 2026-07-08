# First-stage gear transmission_03

细网格、加摩擦 RecurDyn 基准工况。该目录中的 `First-stage gear transmission.rmd` 是后续 NexDyn 摩擦齿轮接触建模的唯一来源。

## RMD 摘要

- 模型：`Model2`
- 主动轮：`PART / 3`，`Verbindungsrad_Rechts6`，GGEOM `2`
- 从动轮：`PART / 4`，`Verbindungsrad_Links6`，GGEOM `4`
- 网格：主动轮 `32574` 节点、`65148` 三角面；从动轮 `32569` 节点、`65138` 三角面
- 接触：`GGEOMCONTACT / 1`，`K=100000`，`C=10`，`KORDER=2`，`BPEN=1e-2`
- 摩擦：RMD 中 `D_F_C/S_F_C` 为摩擦相关参数，需在 NexDyn 摩擦后端中单独映射和验证
- 驱动：`RevJoint1.RMotion`，`ROTATION VELOCITY FUNCTION = 2*PI`
- 动力学：`END = 1.0`，`STEP = 1000`
- 单位：`kg, mm, N, s`

## NexDyn 状态

几何抽取已覆盖该工况：

```text
data/input/obj_library/generated_from_rmd/gear_03_fine_friction_right_surface.obj
data/input/obj_library/generated_from_rmd/gear_03_fine_friction_left_surface.obj
```

这些文件由脚本生成，不提交 Git。当前 NexDyn 入口还未启用 `_03` 摩擦工况；下一步需要在现有细网格入口基础上增加摩擦 case 切换，并对比切向响应、关节反力和从动轮角速度波动。
