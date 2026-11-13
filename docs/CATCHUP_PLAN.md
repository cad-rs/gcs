# rcad-gcs 约束求解器调研与补全追赶计划

- 基线：`rcad-gcs 0.1.0`（commit `8094753`，2026-11-13）
- 范围：对照开源与商业 CAD 约束求解器，盘点本项目差距，给出分阶段的补全计划
- 说明：第 1 节以**代码为准**（README / ROADMAP 有几处已落后于实现，见 1.2）
- 进度：**P0 已做完**（2026-11-13），结果、未做项和过程中的新发现见第 8 节。第 1 节描述的是 P0 之前的基线，文中的文件行号对应基线提交

## 0. 结论摘要

**定位。** 诊断（秩式 DOF、冗余 / 冲突分类、依赖组、逐参数可动性）和测试纪律是本项目的强项，与 FreeCAD planegcs 同档。数值内核和交互行为是最弱的两块：固定 λ 的稠密正规方程、没有拖拽、没有选解机制、失败不回滚。词汇上，2D 只有点线圆弧，没有椭圆和样条；装配只有 9 种基础约束。

**先要处理的基线问题**（本机实测，见 1.4）：

- CI 里的 `cargo fmt --check` 和 `cargo clippy -D warnings` 在当前代码上不过；
- 两个清单测试在 Windows 上因 CRLF 失败；
- 装配语料门在 debug 下要跑 411 秒；
- README / ROADMAP 对实现的描述已经落后（解析 Jacobian 早已全覆盖，约束数也对不上）。

**走查出的正确性问题**（见 1.3）：求解失败不回滚；2D 线–圆相切有符号、角度约束分不清 `θ` 与 `θ+π`；自由圆柱多报 2 个 DOF；刚体求解器按全体刚体建列且没有线搜索。

**与业界的差距**（见第 3 节）：对开源，差在内核（信赖域、稀疏、消元）和交互（拖拽、回滚）；对商业，另差更细的分解、面向几何对象的诊断、显式选解、曲线与方程 / 不等式 / 阵列等词汇，以及日志回放和多语言绑定。

**计划**（见第 4、5 节）：九个阶段，1 人全职约 7–9 个月。顺序是 基线止血 → 统一内核 → 数值升级 → 消元与分解 → 诊断产品化 → 拖拽与选解 → 2D 词汇 → 3D 与装配，生态与质量贯穿全程。先内核后词汇，因为三族现在各有一套基础设施，不先统一，每项升级都要做三遍。

**需要你拍板的五件事**在第 6 节：是否引入 `faer`、0.2.0 是否做一次大的破坏性变更、关节 / 运动放在哪一层、2D 词汇与装配谁先、是否投入构造式分解。

---

## 1. 本项目现状（以代码为准）

### 1.1 能力清单

| 维度 | 2D 草图 `sketch` | 3D 几何 `space3d` | 3D 刚体位姿 `rigid3d` |
|---|---|---|---|
| 未知量 | 扁平参数向量；Line 自带 4 个端点坐标，Arc = 圆心 + 半径 + 起止角 | 扁平参数向量；Plane `(n, d)` 4 参数，Cylinder 7 参数 | 每个刚体 6 列（平移 + 四元数流形增量） |
| 几何类型 | Point / Line / Circle / Arc | Point / Line / Plane / Sphere / Cylinder / Circle / Arc（固定标架）/ Bezier 样条 | 刚体（宿主先把几何引用解析成体局部点/方向） |
| 约束变体数 | 24 | 43 | 9（Fixed / Coincident / Distance / Angle / Parallel / Perpendicular / Concentric / PointOnPlane / Lock） |
| Jacobian | 全解析（稀疏行 `(param, ∂)`），CI 里有与有限差分的一致性门 | 全解析，同上 | 点类/方向类解析；Fixed 姿态行与 Lock 旋转行仍为有限差分 |
| 迭代 | 共享内核：正规方程 `(JᵀJ + λI)Δ = −JᵀF`，λ 固定 `1e-10`，步长减半线搜索（≤8 次），收敛后有 polish | 同左，外加平面规范化（gauge fix） | 独立的一套 Gauss–Newton：无线搜索，仅在消元奇异时把 λ 乘 10 重试 |
| 线性代数 | 自写稠密 `Vec<Vec<f64>>` 高斯消元，O(n³) | 同左 | 同左，且每个子系统都按 `6 × 全部刚体数` 建列 |
| 分解 | 参数并查集 → 连通分量；可只解选中的分量 | 同左，另有 `solve_incremental` | 按刚体图连通分量 |
| 诊断 | 逐分量对 `Jᵀ` 做列主元 Householder QR：真实 DOF、冗余 / 冲突 / 部分冗余、依赖组；逐参数可动性（零空间支撑） | 同左，并处理了平面的规范自由度 | 有秩诊断 `diagnose_assembly`，同时保留两套早期启发式（`DOFAnalyzer` 计数式、`SolverDiagnostics` 模式匹配式） |
| 依赖 | 仅 `glam`；`serde` 可选 | 同左 | 同左，但 `rigid3d` 的数据类型没有 serde 派生 |
| 测试 | 生成式语料 320 例（收敛率 / 冲突召回与精度 / 冗余检出 / 无迟滞回环全部卡 100%） | 语料 + 集成测试 | 装配语料 333 例 |

已有的好底子：解析 Jacobian 全覆盖并有 CI 守门；秩诊断比多数开源实现更细（依赖组、逐参数可动性）；语料红线严格；按分量局部求解；跨平台 + wasm32 编译检查。

### 1.2 文档与代码不一致（应先对齐，否则外部读者会低估/误判现状）

| 位置 | 文档说法 | 代码实际 |
|---|---|---|
| `README.md` “How it solves”、`src/lib.rs` 顶部文档 | Jacobian 用中心有限差分 | 2D / 3D 几何两族已全解析（`src/solver.rs:39`、`src/space3d/solver.rs:43`） |
| `ROADMAP.md` §1 第一条 | “只有刚体点类行是解析的” | 同上；刚体方向类行也已解析，只剩 Fixed 姿态与 Lock 旋转行是 FD |
| `README.md` 约束清单 | 2D 19 种 / 3D 27 种 / 刚体 6 种 | 24 / 43 / 9 |
| `README.md` “Known scaling limit” | 归因于“稠密数值 Jacobian” | 真正瓶颈是稠密正规方程 + O(n³) 消元 |
| `ROADMAP.md` §2 “Expose the redundancy analysis” | 诊断还是内部的 | `Sketch::diagnose()` 已公开 |

### 1.3 代码走查发现的具体问题

以下都是读代码确认的，按影响排序。

1. **求解失败不回滚，且线搜索会接受变差的步。** `src/numerics.rs:132-143`：8 次减半都不满足条件时，最后一次试探（α = 1/128）仍留在参数里；未收敛返回时参数停在最后的迭代点。宿主拿到的是被改坏的草图。成熟求解器的约定是“失败则恢复到求解前的状态”。
2. **2D 线–圆相切是有符号的。** `src/constraint.rs:365-381` 的残差是 `有向距离 − r`，圆心只能落在有向直线的固定一侧；初始构形在另一侧时，求解器会把圆翻到对面去。这是解分支（手性）问题的一个具体实例。
3. **2D 角度约束对 `θ` 与 `θ+π` 不可区分。** `src/constraint.rs:310-321` 的残差 `sin(a)·dot − cos(a)·cross` 等价于 `tan θ = tan a`，线段可以整体掉头而残差为零。
4. **刚体 Angle 约束在 0° / 180° 处梯度为零。** `src/rigid3d/solver.rs:868-875` 用 `d1·d2 − cos(target)`，目标角为 0 或 π 时 Jacobian 行退化，且角度无符号。
5. **刚体求解器规模行为差。** `src/rigid3d/solver.rs:220`：每个连通分量的子系统都用 `6 × 全部刚体数` 列，未参与的刚体留下零列，靠 λ 兜底；接地（Fixed）的刚体仍是未知量而不是被消去；`fill_vec3_residual_fd`（`:588-613`）对每一列都整份克隆位姿数组。主循环（`:238-308`）没有线搜索或信赖域。
6. **刚体的平行 / 同轴行结构性亏秩。** Parallel 用 3 行表达秩 2 的条件，Concentric 用 6 行表达秩 4。求解靠正则化能过，但计数式 DOF 全错，诊断必须额外消化这些“自带冗余”。
7. **圆柱的规范自由度未处理。** `src/space3d/diagnosis.rs:326-383` 的 `effective_columns` 只处理了平面；圆柱 7 参数里轴上点沿轴滑动、方向向量缩放是 2 个规范自由度，自由圆柱会多报 2 个 DOF。
8. **DOF 口径不统一。** `Sketch::dof()` / `dof_analysis()`（`src/sketch.rs:180-202`）是“参数数 − 方程数”，与 `diagnose().dof` 的秩口径并存；`rigid3d` 更是三套并存（见 1.1）。
9. **收敛判据不随尺度变化。** `RESIDUAL_TOL = 1e-10` 作用在 RMS 上（`src/numerics.rs:26`）：毫米与米、长度与弧度混在一起没有归一化；RMS 还会随方程数增多而稀释单个大残差。2D / 3D 两族没有 `SolveOptions`，容差和迭代上限都是常量。
10. **没有编辑 API。** `Sketch` 只有 `add_*`；没有删除实体/约束、修改尺寸值、稳定 ID（ID 就是下标）。`tests/persistence.rs:26-42` 的“尺寸往返”实际是每一步从头重建草图再冷启动，并没有覆盖热启动下的尺寸驱动。
11. **三族基础设施是复制粘贴的。** 连通分量 ×3、秩诊断 ×3、求解适配层 ×2 + 刚体独立一套。任何内核升级（稀疏、LM、DM 分解、拖拽）都要做三遍。
12. **样条上点的残差不光滑。** `src/space3d/constraint.rs:2472` 在残差内部做 20 点采样 + 牛顿求最近点；标准做法是把曲线参数 `t` 作为辅助未知量。
13. **零碎性能点。** 每次求解对每个约束调用 2–3 次 `jacobian()` 只为拿“触及的参数集”（`src/graph.rs:60-108`、`src/solver.rs:131-147`），并用 `Vec::contains` 去重（O(n²)）；每次残差求值分配两次 `Vec`。

### 1.4 本机基线（2026-11-13，Windows 11，cargo 1.99.1，debug 构建）

| 检查 | 结果 |
|---|---|
| 单元测试（lib） | 96 通过，2 忽略 |
| `basic` 45、`persistence` 4、`serde_roundtrip` 1、`space3d_corpus` 24、`space3d_integration` 7、文档测试 2 | 全部通过 |
| `tests/corpus`（2D 语料） | 指标门通过（6.4 秒）；`corpus_manifest_matches_generated` **失败** |
| `tests/assembly_corpus` | 指标门通过，但耗时 **411 秒**；`assembly_manifest_matches_generated` **失败** |
| `cargo fmt --all -- --check` | **失败**（大量格式差异） |
| `cargo clippy --all-targets --all-features -- -D warnings` | **失败**（lib 33 条，lib test 40 条） |

- 两个清单测试失败的原因已查明：本机 `core.autocrlf=true`，仓库没有 `.gitattributes`，清单文件检出为 CRLF，而测试与生成器输出逐字节比对。GitHub 的 Windows runner 默认也是 `autocrlf=true`。
- `.github/workflows/ci.yml` 的前两步正是上面这两条 fmt / clippy 命令，按当前代码和这版工具链过不了。推测是从 monorepo 拆出时没有带上工作区级的 rustfmt / lint 配置（仓库根没有 `rustfmt.toml`、`clippy.toml`，`Cargo.toml` 也没有 `[lints]`）。本机没有 `gh`，我没能核对 GitHub 上实际的 CI 状态。
- 装配语料 411 秒。当时我把它归因于 1.3 第 5 条；P0 实测后发现主因另有其事，见第 8 节。

---

## 2. 业界调研

可信度标记：**✔** 表示本次调研抓取到了原始页面或源码确认；**◇** 表示只见于搜索摘要，或来自背景知识、本次未能确认。方括号里的编号对应第 7 节的来源。

### 2.1 开源求解器

**FreeCAD planegcs**（LGPL-2.1+；2026 年仍有求解器提交）[F1][F2][F3]

- 几何与约束：2D 的点、线、圆、弧、椭圆 / 椭圆弧、双曲线弧、抛物线弧、有理 B 样条；36 种求解器级约束，含 AngleViaPoint 系列、PointOnBSpline、SlopeAtBSplineKnot、ArcLength、圆–圆 / 圆–线 / 点–圆距离、Snell。没有不等式，没有偏移曲线图元 ✔
- 数值：DogLeg（默认）、LM、BFGS，另有两子系统的 SQP。回退顺序是 默认算法 → DogLeg → LM → BFGS → 附加“等于初值”临时约束的 SQP；解出的几何无效则拒绝该解。解析梯度；收敛阈值 1e-10，迭代上限 100 ✔
- 分解：连通分量；Equal 约束做参数代换；求解前丢弃冗余约束 ✔
- 诊断：对 `Jᵀ` 做 QR（1000 参数以下用稠密全主元，以上用 SparseQR），DOF = 参数数 − 秩。从 R 读出冲突组，按“出现在最多的组里 → 去掉的 DOF 最少 → 最新添加”挑候选，去掉后重解：能满足就是冗余，否则是冲突。另报部分冗余和依赖参数组 ✔
- 拖拽与选解：被拖点与鼠标哑点重合的临时约束，以较低优先级经 SQP 求解；存在冲突时拒绝拖拽。选解只靠热启动，“DogLeg 不能稳定收敛到同一个解”的问题（#5962）至今未关；有符号点线距离与定向相切是 2025-12 到 2026-04 才加进去的 ✔
- 用户侧的抱怨：中等复杂草图每次操作延迟大（#27319）、过约束漏检（#6201）、对称约束亏秩（#22381）✔

**SolveSpace / libslvs**（GPL-3.0；v3.3 于 2026-10 发布）[SS1]

- 几何与约束：2D / 3D 的点、线段、圆、弧、三次曲线、工作平面、法向（四元数）；38 种约束，含长度比、长度差、等角、弧线等长、SAME_ORIENTATION、WHERE_DRAGGED。没有圆锥曲线和 NURBS，没有不等式 ✔
- 数值：符号表达式树求偏导并做常量折叠；最小范数牛顿步（对 `Aᵀ` 做稀疏 QR）；步长减半最多 8 次；迭代上限 50；要求每个 `|f|` 都低于阈值才算收敛 ✔
- 分解：按 group 求解；`a − b = 0` 做参数合并；单未知量方程先单独解。未知量上限 2048；未见连通分量拆分 ✔
- 诊断：DOF = n − 秩；秩亏时区分“冗余但可解”和不一致；逐条移除约束重测秩来找该删哪条（代码注明只对单条冗余有效）；不收敛时列出残差超限的约束；逐参数测试是否自由 ✔
- 拖拽：被拖参数的列缩放为 1/20，使其变化最小 ✔

**OndselSolver**（FreeCAD 1.0 起的装配求解器；LGPL-2.1；原仓库 2025 年归档，FreeCAD 的分支仍在维护）[F4][F5]

- 模型：刚体，每体 7 个坐标（3 平移 + 4 欧拉参数）。关节库二十余种：Fixed、Revolute、Cylindrical、Translational、Spherical、Universal、Planar、Gear、RackPinion、Screw 等，另有运动驱动和转动 / 平移限位 ✔
- 数值：牛顿法，未知量是坐标加拉格朗日乘子，目标里带“离上一位姿的加权距离”，结果是加权最小移动。解析 Jacobian，自写稀疏高斯消元；迭代上限 25，容差 1e-6 ✔
- 诊断：矩阵奇异时用全主元消元找出冗余方程，移除后重解。不区分冗余与冲突，未见 DOF 计数 ✔
- 拖拽：被拖零件的平移权重 1e6（基线 1e3）；限位由活动集循环事后施加 ✔
- 用户侧的抱怨：多解时结果不确定（#33557）、圆柱关节对初始对齐和选择顺序敏感（#16055、#16057、#17625）、静默失败（#31349）✔

**ezpz**（Zoo / KittyCAD，Rust，MIT，v0.2.29）[E1]

- 几何与约束：2D 的点、线、圆、弧；25 种约束 ✔
- 数值：LM 自适应阻尼；用 `faer` 的稀疏 Cholesky 解 `JᵀJ + λI`，符号分解预先算好，Jacobian 稀疏模式缓存；偏导由 sympy 生成后手写进代码；迭代上限 35 ✔
- 诊断：稠密列主元 QR 求零空间，列出欠约束变量；给出未满足约束列表。没有冗余 / 冲突分类，没有分解，没有拖拽 API ✔
- 特色：约束优先级与权重，返回“能满足的最高优先级集合”；测试设施齐全（proptest、fuzz、criterion、CodSpeed、cargo-mutants）✔

**其他**

- Rust 生态里的同类还有：`solverang` 0.1.0（2026-03，Apache-2.0；牛顿加 LM 回退，`faer` 稀疏，子问题并行，含 2D 草图 / 3D 草图 / 装配三个模块）✔[X9]；`arael-sketch` 0.9.3（2026-10，MIT；LM，编译期符号求导，表达式尺寸，DOF 分析，最小冲突集）✔[X10]。CADmium 的 ISOtope 已于 2025-03 归档 ✔[X1]。
- Dune3D 直接内置 SolveSpace 的求解器并为提速打了补丁（能符号求解的方程先符号求解），加了 Bézier 与等曲率约束 ✔[X2]。Blender 的 CAD Sketcher 用 py-slvs，求解失败时草图保持不变 ✔[X3]。
- JSketcher 把约束转成多项式，线性的先消元 / 代换，再拆成孤立簇逐簇 dog-leg → LM；拖拽时用粗容差（1e-3）✔[X4]。
- GeoSolver（2009）是开源里唯一的构造式实现：簇重写，区分结构性与偶发性的过 / 欠约束，用手性约束或原型草图选解；但只支持点加距离 / 角度 ✔[X5]。
- CadQuery 的装配用 casadi + ipopt 最小化代价和，草图用 nlopt SLSQP 加有限差分，没有 DOF 诊断 ✔[X6]。SALOME SHAPER 的后端是 planegcs ✔[X7]；`@salusoft89/planegcs` 是 planegcs 的 WASM 移植 ✔[X8]。

### 2.2 商业求解器

| 产品 | 归属与用户 | 公开可确认的能力 |
|---|---|---|
| **D-Cubed 2D DCM**（v80.0，2026-09）| Siemens。授权方名单含 Autodesk、PTC、Onshape、Shapr3D、IronCAD、Bentley、中望（2010 年）等 ✔[D25][D28]；AutoCAD、NX、Solid Edge 用 2D DCM ✔[D27] | 几何：点线圆、椭圆、圆锥曲线、插值 / 控制点样条（刚性 / 可缩放 / 柔性）、宿主参数曲线、偏移曲线（经 PGM）、刚性集与可缩放集 ✔[D1][D12]。约束：常规全套，加阵列、定向距离、曲线长 / 弧长、等值组、样条 G1–G3、方程与上下限不等式 ✔[D1][D6][D7][D10]。选解：half-space、alignment、help point ✔[D6]。求解：可选的欠定义策略（含最小移动）、局部化求解、分区多线程、快速“是否已满足”检查、行为版本化 ✔[D1][D11][D10][D5]。诊断：逐几何 / 逐尺寸状态，欠定义几何“怎么动”的描述（绕某点转、沿某向移），约束推断 ✔[D10][D7][D1]。工程：线程安全、WASM 构建、Python / .NET 封装 ✔[D6][D5] |
| **D-Cubed 3D DCM**（公开资料到 v61.0，2023-11）| 同上；用于装配定位、运动学、3D 草图、直接建模 ✔[D2] | 几何含平面、柱、锥、球、环、参数曲线 / 曲面、旋转面、扫掠面 ✔[D2][D23][D24]。限位不等式在拖拽中生效 ✔[D22]。两种欠定义策略：“移动最少的几何”或“相对移动量最小” ✔[D2]。自动分区并行 ✔[D17]。Siemens 明确说优先用特殊情形的代数解法而非通用迭代 ✔[D23]。诊断给出刚性集的滑动方向与转动轴 ✔[D21] |
| **LGS 2D / 3D** | LEDAS 2001 年起开发，2011 年知识产权售予 Bricsys；用于 BricsCAD、Cimatron 等 ✔[L1][L3]。最后一次公开的组件发布是 V9.0（2014）✔[L3] | 自研非线性求解器，基于模式的几何分解加启发式，在 3000+ 工业算例上调优 ✔[L1][L4]。欠定义时取离初始位置最近的解 ✔[L2]。NURBS 在求解器内带解析导数求值，宿主曲线 / 曲面走回调 ✔[L6][L10]。有控制手性的 alignment / orientation 属性 ✔[L8]、阵列、变量与方程、表格 / 枚举约束 ✔[L7][L11][L12]。诊断给出刚性集剩余自由度的方向与轴、过定义回路 ✔[L6]。有 API 日志与版本化二进制转储 ✔[L6] |
| **C3D Solver** | C3D Labs；用于 KOMPAS-3D、nanoCAD、Renga、Altium Designer ✔[C8] | 2D：样条（权重可编辑）、参数曲线、偏移曲线、刚性集、阵列；约束含“曲线某百分比处的点”、角平分、等曲率、固定样条导数；尺寸可为驱动 / 区间 / 变量；一般方程与线性方程 ✔[C2][C5]。3D：齿轮传动、凸轮、区间距离 / 角度、阵列 ✔[C3][C6][C7]。求解：子任务分解、稀疏矩阵、挑选充分约束子集的“Planner” ✔[C1]。诊断：逐约束 良态 / 病态 / 冗余，逐点 DOF 及运动向量 ✔[C2][C5]。有 API 调用日志；2D 有 JavaScript 移植 ✔[C1] |
| **Spatial CDS** | 基于 Dassault Systèmes V5 / V6 的约束技术 ✔[S1] | 12 类几何，含抛物线、双曲线、NURBS 曲线 / 曲面 ✔[S1]。**软约束**（在常规约束之后满足）、带边界的变量、线性方程、help point 与手性控制 ✔[S1][S2]。关节（平面、转动、移动、圆柱、球、万向）加限位与逆运动学 ✔[S1]。三种模式：Update、Move-Under-Constraints、Simulation，后两者在目标不可达时返回最近可达位置而不是失败 ✔[S1]。诊断：逐约束 / 逐实体状态、剩余平移 / 转动 DOF、冗余与冲突定位、自动约束、回放捕获 ✔[S1][S2] |
| **Autodesk VCS** | Autodesk 自研的 2D / 3D 变分约束求解器 ✔[A1] | 公开信息很少。Fusion 的求解器按实体报告是否全约束，并把草图判为过约束 / 不可解 ✔[A2]；AutoConstrain 于 2025-01 发布 ◇[A3] |
| **国内** | — | 华天 CrownCAD 的求解器叫 **DCS**（DGM 是它的建模内核），第三方文章称其对标 DCM 且 API 兼容，未见约束清单或基准 ✔[H1]。中望 2010 年授权了 2D / 3D DCM ✔[D28]，2024-06 宣布自研 ZGS，未见功能细节 ✔[Z1] |

未能确认的部分：D-Cubed 的求解模式名称、状态码名称、回调式（Frustum）API、增量 / 撤销求值等 API 层细节，公开页面上查不到，只能算背景知识 ◇。LGS 与 DCM 的原始论文被付费墙挡住，具体数值方法未确认。

把几家的共同点归拢，商业级求解器与本项目差别最大的是这几类能力：

- **显式选解手段**：half-space、alignment、help point。它们是 API 的一等公民，而不是靠初值碰运气。
- **欠定义行为可选且有名字**：最小移动、移动最少几何、局部化、软约束。
- **诊断面向几何对象**：逐几何状态，加“还能怎么动”的方向与轴。
- **拖拽是独立模式**：跟手，限位在拖拽中生效，目标不可达时给最近可达位置。
- **曲线是完整的一类公民**：样条 / NURBS、宿主参数曲线回调、偏移曲线、G2 / G3。
- **方程、变量、不等式、阵列、刚性集**全部内建。
- **工程设施**：API 日志与回放、行为版本化、线程安全、WASM 与多语言封装、数千个工业模型的回归库。

### 2.3 算法与文献要点

**分解。** 图构造式分解有两条路线：Owen 的自顶向下（在关节点对处拆图直到只剩三角形，再反向闭式求解，D-Cubed 2D DCM 属于这一系 ◇[R1]），与 Fudos–Hoffmann 的自底向上簇合并 ✔[R2]。Hoffmann–Lomonosov–Sitharam 的 DR-planner / Frontier 用网络流找最小稠密子图，2D / 3D 通用 ✔[R3]。方程层面的做法是二部图最大匹配加 Dulmage–Mendelsohn 分解，再做强连通分量得到块三角求解计划，中等工作量就能拿到结构诊断和提速 ✔[R4]。3D 没有刚性的组合刻画（“双香蕉”反例）◇[R5]，实用的 3D 求解器靠见证点上的数值秩加刚性集缩并。

**过 / 欠约束检测。** 结构法快，但漏掉由几何定理引起的依赖；数值秩在任意一点求值可能伪奇异 ✔[R1]。见证构形法（Michelucci–Foufou）在一个满足关联约束的一般位置上求秩，能抓到非结构性依赖；已有秩例程时实现成本低 ✔[R6][R1]。所有已发表的最小过约束集算法都是贪心的 ✔[R1]。

**数值求解。** 正规方程使条件数平方，SolveSpace 的一个分支 PR 指出这会让本来可解的草图停滞 ◇[R7]——这正是本项目当前的做法。成熟实现的共同点：精确导数、消去平凡方程、基于 QR 的最小范数步、信赖域全局化加回退链（细节见 2.1）。

**选根与交互。** 恰定系统的根数随规模指数增长；主流办法是保持草图的手性 / 顺序类型，或加显式选解约束 ✔[R2]。尺寸变化有“临界值”，越过时子问题退化 ◇[R8]；工程上的对策是细分驱动尺寸、逐步热启动、失败或手性翻转时二分。

**曲线。** 曲线上点与相切的标准做法是每个接触点加一个曲线参数未知量 `t`，方程为 `C(t) = P` 和 `C′(t) × d = 0`；弧长要在残差里做数值积分；G2 是曲率表达式相等。

**装配。** 几何上冗余的配合是常态（例如两个面重合再加一个轴重合），求解器必须静默容忍一致性冗余，只报真正的冲突 ◇[R9]。旋转用绕当前位姿的局部指数映射增量（3 个未知量，无规范约束）。FreeCAD 在拖拽时把固定关节连起来的零件缩并成一个求解体 ◇。

**近期动向（2019–2026）。**

- 新的 Rust 求解器：Zoo 的 `ezpz`（MIT，2026-05 公布，README 称 2000×2000 的系统约 2.9 ms）✔[R10]；`arael-sketch`（LM，编译期符号求导，带最小阻塞集分析）✔[R11]。
- 学习式方法都以求解器为后端：SketchGraphs、Vitruvion 生成图元和约束后交给标准求解器 ◇[R12]；Autodesk 的 AutoConstrain 论文（ICCV 2025）用求解器反馈做强化学习，全约束率 93%，纯监督微调只有 34% ◇[R13]。对本项目的含义：快而确定的 DOF / 过约束反馈是这类模型需要的接口。
- 没有找到公开的 GCS 横向基准，需要自建。

**可用的 Rust 构件。**

| crate | 能力 | 备注 |
|---|---|---|
| `faer` 0.24 | 稠密 LLT / LU / QR / **列主元 QR** / SVD；稀疏 LLT / LU / **QR**，带 AMD / COLAMD ✔ | MIT；未到 1.0 但活跃；`rayon` 是默认特性，wasm 下要关 |
| `nalgebra` | ColPivQR、SVD、Cholesky ◇ | 成熟；`nalgebra-sparse` 只有稀疏 Cholesky，没有稀疏 QR ◇ |
| `sprs` / `sprs-ldl` | 稀疏格式；LDL | `sprs-ldl` 是 **LGPL-2.1** ◇，与本库的 MIT / Apache 双许可不合 |
| `levenberg-marquardt` | MINPACK 移植，稠密 ◇ | 只适合小系统 |
| `argmin` | 信赖域（Dogleg、Steihaug）、BFGS ◇ | 通用框架，有额外开销 |
| `num-dual` | 前向模式对偶数求导 ◇ | 可用于表达式约束（P6.5）和导数自检 |
| Enzyme `std::autodiff` | 编译器级自动微分 ◇ | 仅 nightly，目前不能用于稳定版 / wasm 库 |
| `petgraph` | Tarjan SCC、一般图匹配；**没有 DM 分解** ◇ | 二部图 Hopcroft–Karp 加 DM 自己写约 200 行，不值得为此引入依赖 |

---

## 3. 差距分析

### 3.1 差距清单

级别：**急** = 正确性或基线问题，放在 P0；**高** = 不补就谈不上可用的求解器；**中** = 成熟产品的标配；**低** = 锦上添花。

**数值内核**

| 编号 | 差距 | 对标 | 级别 |
|---|---|---|---|
| N1 | 线性步走正规方程（条件数平方），λ 固定，没有信赖域 | planegcs：DogLeg 默认，另有 LM、BFGS；SolveSpace：稀疏 QR 求最小范数步 | 高 |
| N2 | 稠密 O(n³) 消元；迭代内反复分配；刚体子系统按全体刚体建列 | planegcs、SolveSpace 用稀疏 QR；C3D 用稀疏矩阵；ezpz 用 `faer` | 高 |
| N3 | 收敛判据是绝对 RMS，不随模型尺度与量纲变化 | 3D DCM 有模型精度 / 分辨率告警 | 中 |
| N4 | 失败不回滚，线搜索会接受变差的步，没有状态码 | SolveSpace 有区分“冗余但可解 / 未收敛”的状态码；CDS 在不可达时返回最近可达位置 | 急 |
| N5 | 没有回退算法链 | planegcs 在算法之间回退 | 中 |
| N6 | 刚体仍有两处有限差分行 | 各家都是精确导数 | 低 |
| N7 | 没有并行 | DCM 按分区多线程 | 低 |

**结构与分解**

| 编号 | 差距 | 对标 | 级别 |
|---|---|---|---|
| S1 | 重合 / 水平 / 竖直 / 相等 / 固定等平凡方程不消元，全部进牛顿 | SolveSpace 合并 `a = b` 的参数并先解单未知量方程；planegcs 对 Equal 做参数代换；JSketcher 消去线性多项式 | 高 |
| S2 | 只有连通分量，没有 DM 分解与块三角求解 | 文献中的方程级标准做法；商业产品都有更细的分解 | 中 |
| S3 | 没有小块闭式解，更没有构造式分解 | DCM 优先代数特例；LGS 用基于模式的分解 | 中（小块）/ 长期（构造式）|
| S4 | 结构信息不缓存，每次求解重算分量与触及参数集 | — | 中 |
| S5 | 三族基础设施复制粘贴 | — | 高（内部）|

**诊断**

| 编号 | 差距 | 对标 | 级别 |
|---|---|---|---|
| D1 | DOF 口径不统一（计数式与秩式并存，刚体三套并存） | 各家都以秩为准 | 急 |
| D2 | 只有逐参数可动性，没有逐实体状态，也没有“还能怎么动”的几何描述 | DCM、LGS、C3D、CDS 都给出方向与轴 | 高 |
| D3 | 只在当前线性化点诊断，奇异构形会被当成冗余 | 见证构形法 | 中 |
| D4 | 有依赖组，但没有最小冲突集与“建议移除哪条” | planegcs 去掉候选后重解来判定；SolveSpace 逐条移除；CDS 有冲突定位 | 中 |
| D5 | 没有“加这条约束会怎样”的增量判定 | — | 中 |
| D6 | 没有自动约束 / 自动标注建议 | DCM、LGS、CDS、Fusion 都有 | 中 |
| D7 | 没有从动（参考）尺寸 | planegcs 的从动尺寸只求值、不参与求解与秩分析；DCM 的 rigid dimension 不会造成过定义 | 中 |
| D8 | 规范自由度处理不完整：圆柱多报 2 个 DOF；刚体平行 / 同轴行自带冗余 | — | 急 |

**几何与约束词汇**

| 编号 | 差距 | 对标 | 级别 |
|---|---|---|---|
| V1 | 2D 没有椭圆、圆锥曲线、样条 / NURBS、偏移曲线 | 商业产品全部具备 | 高 |
| V2 | 2D 基础约束缺口：点–线 / 线–线 / 圆–线距离、共线、两点水平竖直、端点相切、弧长、等角、长度比等 | 同上 | 高 |
| V3 | 没有不等式与限位 | DCM、C3D、CDS；开源里只有 OndselSolver 的关节限位，planegcs 和 SolveSpace 都没有 | 中 |
| V4 | 没有变量与方程 | DCM、LGS、C3D、CDS | 中 |
| V5 | 没有阵列、镜像、刚性组 | 商业产品全部具备 | 中 |
| V6 | 3D 没有工作平面、圆锥、圆环、距离族 | 3D DCM、LGS、C3D、CDS | 中 |
| V7 | 装配只有 9 种基础约束：没有几何对配合、关节库、限位、齿轮等耦合、运动驱动 | CDS 有关节与逆运动学；C3D 有齿轮、凸轮；Onshape 有成套 mate | 中到高，看产品重心 |

**交互行为**

| 编号 | 差距 | 对标 | 级别 |
|---|---|---|---|
| I1 | 没有拖拽 API，没有软约束或优先级 | planegcs 用两级优先级的 SQP；SolveSpace 对被拖参数做列加权；商业产品都有独立的拖拽模式 | 高 |
| I2 | 没有选解 / 手性保持机制，线–圆相切和角度约束本身就有分支问题 | DCM 的 half-space / alignment / help point；LGS、CDS 的手性控制。开源这边也是短板：planegcs 到 2026 年才补上定向相切，多解不确定是 FreeCAD 的长期抱怨 | 高 |
| I3 | 大的尺寸跳变直接一步求解，没有连续化 | 文献中的驱动尺寸细分 | 中 |
| I4 | 欠约束时的“最小移动”只是最小范数步的副产品，没有权重，也没有保证 | DCM、LGS、CDS 都有具名的欠定义策略 | 中 |

**API 与工程**

| 编号 | 差距 | 对标 | 级别 |
|---|---|---|---|
| A1 | 没有稳定 ID，不能删除实体 / 约束，不能改尺寸值 | 所有可嵌入的求解器都支持增量增删改 | 高 |
| A2 | 2D / 3D 没有 `SolveOptions`；错误靠 `bool` 和 panic | — | 急 |
| A3 | serde 覆盖不全；没有版本化交换格式；没有日志回放 | LGS、C3D、CDS 都有 API 日志或回放 | 中 |
| A4 | 没有绑定 | 2D DCM 有 WASM 构建与 Python / .NET 封装 | 中 |
| A5 | 文档落后于代码；没有示例与指南 | — | 急 |

**质量体系**

| 编号 | 差距 | 对标 | 级别 |
|---|---|---|---|
| Q1 | CI 基线不绿（CRLF 清单、fmt、clippy） | — | 急 |
| Q2 | 没有性能基准与回归门 | — | 急 |
| Q3 | 语料全是生成式的，没有真实草图 | LGS 在 3000+ 工业算例上调优 | 中 |
| Q4 | planegcs 对拍驱动不在本仓库，也不在 CI 里 | — | 中 |
| Q5 | 没有性质测试与模糊测试 | ezpz 有 proptest、fuzz、criterion、CodSpeed、cargo-mutants | 中 |
| Q6 | 持久性测试每步重建草图，不覆盖热启动 | — | 中 |

### 3.2 对比矩阵

● 有；◐ 部分有；○ 没有；`?` 未能确认。“商业”一列指 DCM、LGS、C3D、CDS 中至少两家公开确认。

| 能力 | 本项目 | planegcs | SolveSpace | ezpz | 商业 |
|---|---|---|---|---|---|
| 精确导数 | ●（刚体余 2 处 FD）| ● | ●（符号）| ● | ? |
| 信赖域与回退链 | ○ | ● | ○（牛顿 + 减半）| ◐（仅 LM）| ? |
| 稀疏线性代数 | ○ | ◐（仅诊断）| ● | ● | ● |
| 平凡方程消元 | ○ | ● | ● | ○ | ● |
| 连通分量 | ● | ● | ○（按 group）| ○ | ● |
| 更细的分解（DM / 构造式）| ○ | ○ | ○ | ○ | ● |
| 秩式 DOF | ● | ● | ● | ● | ● |
| 冗余与冲突分类 | ● | ● | ◐ | ○ | ● |
| 依赖组 / 冲突组 | ● | ● | ○ | ○ | ● |
| 逐参数可动性 | ● | ● | ● | ● | ● |
| 逐实体状态与自由度方向 / 轴 | ○ | ○ | ○ | ○ | ● |
| 自动约束建议 | ○ | ○ | ○ | ○ | ● |
| 从动尺寸 | ○ | ● | ● | ○ | ● |
| 拖拽 | ○ | ● | ● | ○ | ● |
| 优先级 / 软约束 | ○ | ◐（两级）| ○ | ● | ● |
| 显式选解 / 手性 | ○ | ◐ | ◐ | ○ | ● |
| 失败回滚与状态码 | ○ | ● | ● | ◐ | ● |
| 椭圆 / 圆锥曲线 | ○ | ● | ○ | ○ | ● |
| 样条 | ◐（仅 3D Bezier）| ●（有理 B 样条）| ◐（三次）| ○ | ● |
| 不等式 / 限位 | ○ | ○ | ○ | ○ | ● |
| 变量与方程 | ○ | ○（在 FreeCAD 表达式层）| ◐（应用层）| ◐ | ● |
| 阵列 / 刚性组 | ○ | ◐（Block）| ○ | ○ | ● |
| 3D 几何约束 | ● | ○ | ● | ○ | ● |
| 装配 | ◐（9 种基础约束）| 另由 OndselSolver 承担 | ◐ | ○ | ● |
| 增删改与稳定句柄 | ○ | ● | ● | ? | ● |
| 绑定（C / WASM / Python）| ○（仅 wasm 编译检查）| ◐（有 WASM 移植）| ● | ◐（WASM）| ● |
| 基准 / 性质 / 模糊测试 | ○ | ◐ | ◐ | ● | ●（工业回归库）|

读这张表的三个结论：

1. **诊断是本项目已经站住的位置。** 秩式 DOF、冗余 / 冲突分类、依赖组、逐参数可动性都有，与 planegcs 同档，比 SolveSpace 和 ezpz 细。
2. **数值内核和交互行为是最弱的两块。** 固定 λ 的稠密正规方程、没有拖拽、没有选解、失败不回滚——这四项在 planegcs 和 SolveSpace 里都是十几年前就有的东西。
3. **与商业产品的结构性差距在三处：** 更细的分解、面向几何对象的诊断（状态 + 方向 / 轴）、显式选解。其中选解在开源阵营里也普遍没做好，是可以做出差异的地方。

定位上要注意：Rust 生态里已经有 ezpz、solverang、arael-sketch 三个同类，其中 ezpz 的内核（稀疏 LM + 符号分解缓存）和测试设施已经领先于本项目。本项目的差异点是三族统一、诊断更细、语料红线更严；要保住这个差异，P1–P5 比 P6 / P7 的词汇扩张更要紧。

---

## 4. 追赶计划

### 4.0 总体思路

三条原则：

1. **先统一内核，再加功能。** 现在三族各有一套分解 / 诊断 / 求解适配，内核升级要做三遍。阶段 1 先把三族收拢到同一条流水线上，后面每项升级只做一次。
2. **每个阶段都有可量化的验收门。** 沿用现有语料红线的做法：新能力先有语料族和门槛，再写实现。文中的性能数字是目标值，阶段 0 测出基线后再校准。
3. **破坏性 API 变更集中在一次（阶段 1，对应 0.2.0）。** 稳定 ID、尺寸表、刚体约束类型化一起改，宿主（rcad-sketcher / rcad-assembly）只适配一次。

阶段依赖：`P0 → P1 → P2 → P3`；`P4`、`P5` 依赖 `P2`/`P3`；`P6`、`P7` 只依赖 `P1`，可与 `P3–P5` 并行；`P8` 贯穿全程。工作项后的 `[N1]` 等标记对应第 3 节的差距编号。

### P0 基线止血与度量（约 1 周）

目标：CI 真实可信，文档与代码一致，后续每一步都有数字可比。

| 编号 | 工作项 |
|---|---|
| P0.1 | 修 CI 基线：加 `.gitattributes`（`* text=auto eol=lf`）解决清单测试在 Windows 上的 CRLF 失败；清掉 fmt / clippy 存量问题（见 1.4）`[Q1]` |
| P0.2 | 对齐 `README.md` / `ROADMAP.md` / `src/lib.rs` 文档（见 1.2）`[A5]` |
| P0.3 | 建 criterion 基准：2D 规模阶梯（50 / 200 / 1k / 5k 参数的矩形链、距离网格）、装配链与闭环、单次尺寸修改后的重解；结果入库作为基线，CI 上超过阈值的退化报错 `[Q2]` |
| P0.4 | 装配语料提速：子系统只为本分量的刚体建列、接地刚体直接消去、FD 不再整份克隆位姿（1.3 第 5 条）。这是 debug 下 411 秒的直接原因，属于改几十行就能见效的问题 `[N2]` |
| P0.5 | 把 planegcs 对拍驱动从父仓库迁入本仓库 `tools/`，再加一个 SolveSpace（py-slvs）驱动；作为不阻塞的 CI 任务产出对比报告（收敛、迭代数、耗时、DOF 是否一致）`[Q4]` |
| P0.6 | 失败回滚 + 线搜索不接受变差的步；`SolveStatus` 枚举（收敛 / 达到迭代上限 / 奇异 / 发散）；2D、3D 两族加 `SolveOptions`（容差、迭代上限）`[N4][A2]` |
| P0.7 | 小修：`Sketch::dof()` 改为秩口径，计数式保留为 `structural_dof_estimate()`；圆柱规范自由度；触及参数集去掉 O(n²) `[D1][D8]` |

验收：三平台 CI 全绿；基准报告入库；装配语料 debug 下 < 60 秒；现有语料红线不变。

### P1 统一内核与数据模型（约 3–4 周，0.2.0，破坏性）

目标：三族共用一条 `预处理 → 分解 → 求解 → 诊断` 流水线；宿主可以增删改而不必重建系统。

| 编号 | 工作项 |
|---|---|
| P1.1 | 抽象出统一的系统接口：参数个数、方程个数、残差、稀疏 Jacobian、**不求值就能拿到的关联结构**（每个方程触及哪些参数）。2D / 3D 几何 / 刚体各实现一次，分解、诊断、求解循环只保留一份 `[S5][S4]` |
| P1.2 | 流形参数化（Ceres 的 Manifold 思路）：把“存储坐标比物理自由度多”的对象统一处理——刚体四元数、平面 `(n, d)`、圆柱轴线、单位方向。迭代在切空间里走，更新后回到流形。这样 `normalize_planes`、`effective_columns`、圆柱多报 DOF、刚体单独一套循环这四处特例都可以删掉 `[D8][S5]` |
| P1.3 | 稳定句柄：代际索引式的 `EntityId` / `ConstraintId`；`remove_entity`、`remove_constraint`、约束启用/抑制 `[A1]` |
| P1.4 | 尺寸表：约束里的数值（距离、角度、半径）移到独立的尺寸表，约束只持有尺寸句柄；`set_dimension(id, value)` 不改变拓扑，结构缓存可以跨求解复用；尺寸带“驱动 / 从动”标志 `[A1][D7][S4]` |
| P1.5 | 刚体约束类型化：用带数据的枚举替换“一堆 `Option` 字段”的结构体；Parallel 用 2 行、Concentric 用 4 行的满秩表达（在方向的切平面基上投影）；补 serde 派生 `[D8][A3]` |
| P1.6 | 下线 `DOFAnalyzer` 与 `SolverDiagnostics` 两套启发式，名字保留为秩诊断之上的薄封装并标记弃用，给 rcad-assembly 一个迁移窗口 `[D1]` |
| P1.7 | 错误模型：`validate()` 对引用了错误类型实体、越界句柄返回类型化错误，替换现在的 `bool` 返回与潜在 panic `[A2]` |

验收：三族经同一个入口求解与诊断；重复的分量 / 诊断代码删除；三套语料红线全部保持；附一份宿主迁移说明。

### P2 数值内核升级（约 3–4 周）

目标：病态和大规模系统都稳定收敛；规模上限从“几十个自由参数”提到几千。

| 编号 | 工作项 |
|---|---|
| P2.1 | Jacobian 改为固定稀疏模式的 CSR 存储；迭代内零分配（工作区复用）`[N2]` |
| P2.2 | 线性步不再走正规方程。小块（稠密）：对 `Jᵀ` 做列主元 QR 并按秩截断，直接得到最小范数步；同一次分解同时给出秩、冗余行、零空间，求解与诊断共用。大块：经 `faer` 做稀疏 QR / 稀疏 Cholesky（可选特性 `sparse`），带 AMD/COLAMD 排序，符号分解跨迭代、跨重解缓存 `[N1][N2]` |
| P2.3 | 全局化：Levenberg–Marquardt（按增益比自适应 μ）与 Powell DogLeg 两种信赖域；失败时按 `DogLeg → LM → BFGS` 回退；线搜索改为对 `½‖F‖²` 的 Armijo 充分下降 `[N1][N5]` |
| P2.4 | 尺度：按约束类型做行均衡（角度类残差乘特征长度，使各行量纲一致）、列缩放；容差改为“线性分辨率 + 角度分辨率 + 模型特征尺寸”，收敛判据用缩放后残差的最大范数加步长范数 `[N3]` |
| P2.5 | 退化几何统一防护（零长线段、零半径、重合点处的梯度），NaN / Inf 检出并转成状态码 `[N4]` |
| P2.6 | 刚体 Fixed 姿态行与 Lock 旋转行改解析（对数映射的 Jacobian），全库不再有 FD 行 `[N6]` |
| P2.7 | 可选特性 `rayon`：按连通分量并行（分量之间天然独立，比并行单个 Jacobian 划算）`[N7]` |

验收：`large_system_makes_progress` 收紧为完全收敛；规模基准 1k 参数 < 5 ms、5k 参数 < 50 ms（release，冷启动，目标值）；语料“大幅漂移”探针的收敛率和分支保持率相对基线只升不降；与 planegcs 对拍，同语料收敛率不低于对方。

### P3 预处理消元与结构分解（约 3–4 周）

目标：把系统变小、变成一串小块；结构信息只算一次。

| 编号 | 工作项 |
|---|---|
| P3.1 | 平凡方程消元：对 `p = q + c` 形式的方程做参数别名（并查集带偏移）。覆盖 Coincident、Horizontal / Vertical、EqualRadius、Concentric、Fixed、Radius / Diameter、DistanceX / Y。别名图里出现的矛盾或重复要带着原约束 ID 回报给诊断，不能因为被消掉就从报告里消失 `[S1]` |
| P3.2 | 结构分解：方程–变量二部图最大匹配（Hopcroft–Karp）→ Dulmage–Mendelsohn 分解得到结构性过约束 / 恰定 / 欠约束三部分 → 恰定部分做强连通分量（Tarjan）和拓扑序，得到块三角求解计划，逐块小牛顿。欠约束部分仍整体求最小范数解 `[S2]` |
| P3.3 | 小块闭式解：单方程、两圆交、线圆交等 1–2 维块直接解析求根，并选离当前构形最近的根。闭式解既快又让选根变成显式决策 `[S3]` |
| P3.4 | 结构缓存：分量、匹配、块顺序按“拓扑版本号”缓存；只改尺寸值时全部复用 `[S4]` |
| P3.5 | 刚性组缩并：2D 的 Block、装配的 Lock 把多个对象并成一个求解体 `[V5]` |

注意点：逐块求解会改变欠约束情形下的选解行为，必须用无迟滞回环语料验证；整体求解保留为回退路径。

验收：真实草图样本上未知量个数下降 40% 以上；500 实体草图改一个尺寸后的重解 p50 < 1 ms（目标值）；结构诊断与数值秩在语料上无矛盾。

### P4 诊断产品化（约 2–3 周）

目标：从“能算出秩”到“能直接驱动 UI”。

| 编号 | 工作项 |
|---|---|
| P4.1 | 三族统一的诊断结果：总 DOF、逐实体状态（定义完全 / 欠定义 / 过定义 / 冲突 / 未求解）、逐约束状态、依赖组 `[D1][D2]` |
| P4.2 | 自由度的几何描述：把零空间运动限制到单个实体上，拟合成“沿某方向平移 / 绕某点转动 / 半径可变”等，而不只是逐参数的布尔值 `[D2]` |
| P4.3 | 见证构形（witness）诊断：除当前点外，在一个仍满足关联类约束的随机扰动构形上再算一次秩，区分“结构性冗余”与“只在当前构形奇异”（现有语料里的 Singular 类）`[D3]` |
| P4.4 | 冲突定位：在依赖组内做删除过滤（去掉候选约束后重解）得到最小冲突集；按“最近添加”排序给出建议移除项 `[D4]` |
| P4.5 | what-if 接口：`would_add(&constraint)` 用缓存的行空间正交基做增量判定（独立 / 冗余 / 冲突），不重做整体 QR `[D5]` |
| P4.6 | 自动约束 / 自动标注：从零空间出发，贪心挑选能提高秩的候选尺寸，给出一组把 DOF 降到 0 的建议 `[D6]` |
| P4.7 | 从动尺寸：测量接口（当前值，可选对驱动尺寸的灵敏度）`[D7]` |

验收：语料新增逐实体状态标注，准确率 100%；种子故障上最小冲突集与预期一致；what-if 与完整诊断结论 100% 一致且快一个数量级以上。

### P5 交互：拖拽与解分支稳定（约 3–4 周）

目标：欠约束草图拖得动、拖得稳；改尺寸不翻面。

| 编号 | 工作项 |
|---|---|
| P5.1 | 拖拽会话 API：`begin_drag → drag_to → end_drag`。硬约束照常满足；被拖对象的目标位置作为软目标；其余参数按权重做最小移动。第一版用加权最小范数步（列缩放）加零空间投影，第二版视效果再上两级优先级的 SQP `[I1][I4]` |
| P5.2 | 权重与优先级：被拖对象、外部参考、新画的几何各有默认参数权重，宿主可覆盖；约束可带优先级，按级逐层加入，返回能满足的最高一级（ezpz 的做法），软约束由此统一表达 `[I4][I1]` |
| P5.3 | 手性守卫：创建约束时记录离散签名（相切在哪一侧、角度在哪个象限、点在线的哪一侧、内切还是外切、弧的转向）；迭代中签名翻转则缩步，无法避免时返回 `BranchFlip` 状态。顺带修掉 1.3 第 2、3、4 条 `[I2]` |
| P5.4 | 尺寸跳变连续化：`set_dimension` 变化量大时自适应细分，逐步热启动；失败或签名翻转时二分 `[I3]` |
| P5.5 | 拖拽不可达时的策略：沿拖拽方向走到最远可行处，保留最后一个可行构形 `[I1]` |
| P5.6 | 帧预算：只重解受影响的块、复用分解；拖拽过程中用粗容差，松手后再精解一次；200 实体草图单帧 < 1 ms（目标值）`[S4]` |
| P5.7 | 持久性测试改为真正的热启动尺寸往返，替换现在的“每步重建”`[Q6]` |

验收：脚本化拖拽语料上零翻面、可行域内零不收敛帧，记录 p99 帧时间；尺寸在 10 倍范围内扫描往返无迟滞。

### P6 2D 几何与约束词汇补全（约 5–7 周）

目标：覆盖主流草图器的常用词汇，能 1:1 导入 FreeCAD 草图的约束集。

| 编号 | 工作项 |
|---|---|
| P6.1 | 曲线抽象：统一的“参数曲线”接口（点、一阶、二阶导数），在其上实现通用约束——曲线上点（曲线参数 `t` 作为辅助未知量）、过点相切、过点夹角、法向、曲线长度。新曲线类型只需实现接口 `[V1][V2]` |
| P6.2 | 新几何：椭圆、椭圆弧、B 样条 / NURBS（控制点 + 权重，节点固定；de Boor 求值带解析导数）、偏移曲线；双曲线 / 抛物线弧低优先级 `[V1]` |
| P6.3 | 补齐基础约束：点–线距离、平行线距离、圆–线 / 圆–圆距离、两点水平 / 竖直、共线、等角、长度比、弧长 / 弧角、带跨度的弧上点、弧中点、线对称、端点相切（G1 接续）、样条 G2 `[V2]` |
| P6.4 | 不等式：尺寸上下限、弧扫角 > 0、半径 > 0、同侧；用活动集或松弛变量实现，手性守卫（P5.3）可复用 `[V3]` |
| P6.5 | 变量与方程：用户变量作为额外参数，尺寸之间的代数关系作为方程约束（表达式树自动求导，或宿主提供带梯度的回调）`[V4]` |
| P6.6 | 阵列（线性 / 圆周）、镜像、刚性组，作为派生几何约束 `[V5]` |
| P6.7 | 每个新变体的固定套路：残差 + 解析 Jacobian + FD 一致性门 + 语料族 + planegcs / SolveSpace 映射 |

验收：给出与 planegcs 约束清单的逐项映射表，除明确排除项外全部覆盖；真实草图样本的约束实例覆盖率 ≥ 99%。

### P7 3D 几何与装配（约 5–7 周）

目标：`space3d` 能承载 3D 草图；`rigid3d` 从“9 种基础约束”长成可用的装配求解器。

| 编号 | 工作项 |
|---|---|
| P7.1 | `space3d`：工作平面实体，2D 约束在任意工作平面上复用（SolveSpace 的模型）；圆锥、圆环；点–线 / 线–线 / 面–面距离族；关于平面对称；与 2D 共用曲线抽象 `[V6]` |
| P7.2 | 几何对配合：面–面（重合 / 距离 / 夹角）、轴–轴、点–线、相切（面–柱、柱–柱、球–面），带同向 / 反向标志 `[V7]` |
| P7.3 | 关节库：转动、移动、圆柱、球、平面、固定、万向、螺旋；耦合：齿轮、齿轮齿条、带传动；限位作为不等式 `[V7][V3]` |
| P7.4 | 求解结构：接地消去、刚性组缩并、子装配（刚性 / 柔性）；初值不再“全部放在原点”，先沿约束树做闭式逐件定位再进牛顿 `[N2]` |
| P7.5 | 运动：驱动某个关节坐标逐步求解（复用 P5.4 的连续化），返回每步位姿；逐刚体的剩余自由度用运动旋量（轴线 + 节距）描述；拖拽零件复用 P5 `[V7][D2]` |
| P7.6 | 标准机构语料：四杆、曲柄滑块、Stewart 平台、行星轮系；DOF 与理论值核对 |

验收：500 刚体 / 1500 配合（链 + 闭环）< 100 ms（目标值）；标准机构语料全过；一致性冗余的配合不报错，只有真正冲突才报。

### P8 生态、互操作与质量体系（贯穿，合计约 4 周）

| 编号 | 工作项 |
|---|---|
| P8.1 | 版本化 JSON 交换格式；求解日志（系统 + 选项 + 初始状态的完整转储）与回放命令行工具，用于复现宿主报来的问题 `[A3]` |
| P8.2 | 绑定：WASM（wasm-bindgen + npm 包 + TS 类型）、C ABI（cbindgen）、Python（PyO3）`[A4]` |
| P8.3 | 真实语料：SketchGraphs 抽样导入、FreeCAD 草图导入；覆盖率与收敛率进 CI 指标 `[Q3]` |
| P8.4 | proptest 性质测试（收敛则残差达标；已解系统再解是不动点；删冗余行不改变解）；cargo-fuzz 打退化几何 `[Q5]` |
| P8.5 | mdbook 指南 + `examples/`；crates.io / docs.rs 按阶段发版 `[A5]` |

### 研究轨（1.0 之后，或有余力时并行）

- **R1 构造式分解。** 2D 刚性簇的图构造分解（Owen、Fudos–Hoffmann 一系）加三角形闭式解。这是商业求解器在速度、选根确定性和诊断精细度上的根本来源，也是投入最大的一项。
- **R2 多解枚举。** 同伦延拓找出草图附近的其他解，支撑“切换到另一个解”的交互。
- **R3 学习式自动约束的对接接口。** 这类模型需要的是快而确定的 DOF / 过约束反馈，P4 做完后本库就是合格的后端。

---

## 5. 时间线与版本

| 版本 | 内容 | 累计周数（1 人全职） |
|---|---|---|
| 0.2.0 | P0 + P1（破坏性；P0 的修正本身已不兼容 0.1.0，所以不再单独发 0.1.1） | 4–5 |
| 0.3.0 | P2 + P3 | 10–13 |
| 0.4.0 | P4 + P5 | 15–20 |
| 0.5.0 | P6 | 20–27 |
| 0.6.0 | P7 | 25–34 |
| 1.0.0 | P8 收尾 + API 冻结 | 29–38 |

合计约 7–9 个人月。P1 之后若有第二个人，P6（2D 词汇）与 P7（装配）可以和 P3–P5 并行，日历时间可压到 5–6 个月。这些是粗估，P0 的基线数据出来后应重新校一遍。

## 6. 需要你拍板的决策点

**2026-11-13 已确认：五项均按下表的建议执行。**

| 编号 | 决策 | 我的建议 |
|---|---|---|
| D1 | 依赖策略：是否引入 `faer` | 引入，但放在可选特性 `sparse` 后面；默认路径仍只依赖 `glam`，保住“零依赖”卖点和 wasm 体积。由基准数据决定以后是否默认开启 |
| D2 | 0.2.0 做一次大的破坏性变更 | 做。稳定 ID、尺寸表、刚体类型化拖得越久，宿主迁移越贵 |
| D3 | 关节 / 齿轮 / 运动是放进本库还是留在 rcad-assembly | 纯位置级的（关节、限位、耦合）进本库；碰撞、动力学留在宿主 |
| D4 | P6（2D 词汇）与 P7（装配）谁先 | 取决于 rcad 近期的产品重心；只有一个人时建议 P6 先，草图是所有建模的入口 |
| D5 | 是否投入 R1 构造式分解 | 1.0 前不做。P3 的消元 + DM 分解 + 小块闭式解能拿到大部分收益 |

主要风险：

- **P1 重构面大。** 缓解：三套语料红线全程不许动，先加新流水线再逐族切换，旧路径保留到切换完成。
- **P3 逐块求解改变选解行为。** 缓解：无迟滞回环与分支保持率作为门槛，整体求解保留为回退。
- **`faer` 仍是 1.0 之前的版本。** 缓解：只在一个薄适配层里用它，稠密路径自有实现兜底。
- **词汇扩张拖垮质量。** 缓解：P6.7 的固定套路不满足不合并。

---

## 7. 参考来源

本次调研抓取并核对过的主要来源。`[R7]` 只见于搜索摘要，没有拿到链接。

**开源**

- F1 planegcs 源码：https://github.com/FreeCAD/FreeCAD/tree/main/src/Mod/Sketcher/App/planegcs
- F2 FreeCAD `Sketch.cpp`：https://github.com/FreeCAD/FreeCAD/blob/main/src/Mod/Sketcher/App/Sketch.cpp
- F3 FreeCAD issue（文中的 #编号）：https://github.com/FreeCAD/FreeCAD/issues
- F4 OndselSolver：https://github.com/FreeCAD/OndselSolver
- F5 FreeCAD Assembly：https://github.com/FreeCAD/FreeCAD/tree/main/src/Mod/Assembly/App
- SS1 SolveSpace（`src/system.cpp`、`include/slvs.h`、`CHANGELOG.md`）：https://github.com/solvespace/solvespace
- E1 ezpz：https://github.com/KittyCAD/ezpz
- X1 ISOtope：https://github.com/CADmium-Co/ISOtope
- X2 Dune3D：https://github.com/dune3d/dune3d
- X3 CAD Sketcher：https://hlorus.github.io/CAD_Sketcher/solver/
- X4 JSketcher：https://github.com/xibyte/jsketcher
- X5 GeoSolver：https://geosolver.sourceforge.net/
- X6 CadQuery：https://github.com/CadQuery/cadquery
- X7 SALOME SHAPER：https://github.com/SalomePlatform/shaper/tree/master/src/SketchSolver
- X8 @salusoft89/planegcs：https://github.com/Salusoft89/planegcs
- X9 solverang：https://docs.rs/crate/solverang/latest
- X10 arael-sketch：https://docs.rs/crate/arael-sketch/latest

**商业**

- D1 2D DCM 产品页：https://www.siemens.com/en-us/products/plm-components/d-cubed/2d-dcm/
- D2 3D DCM 产品页：https://www.siemens.com/en-us/products/plm-components/d-cubed/3d-dcm/
- D5 2D DCM v80.0：https://blogs.sw.siemens.com/plm-components/whats-new-in-d-cubed-2d-dcm-and-pgm-version-80-0/
- D6 v79：https://blogs.sw.siemens.com/plm-components/d-cubed-2d-components-version-79-release/
- D7 v76：https://blogs.sw.siemens.com/plm-components/d-cubed-2d-components-version-76-release/
- D10 2D 发布要点：https://blogs.sw.siemens.com/plm-components/d-cubed-2d-components-release-highlights/
- D11 v72.0：https://blogs.sw.siemens.com/plm-components/d-cubed-2d-dcm-version-72-0/
- D12 v71.0：https://blogs.sw.siemens.com/plm-components/d-cubed-2d-dcm-version-71-0/
- D17 3D v59：https://blogs.sw.siemens.com/plm-components/d-cubed-3d-components-version-59-release/
- D21 3D DCM v55.0：https://blogs.sw.siemens.com/plm-components/d-cubed-3d-dcm-version-55-0/
- D22 3D DCM 有界约束：https://blogs.sw.siemens.com/plm-components/d-cubed-3d-dcm-v50-a-closer-look-at-bounded-constraints/
- D23 “It's all about special cases”：https://blogs.sw.siemens.com/news/geometric-constraint-solving-its-all-about-special-cases/
- D24 旋转面：https://blogs.sw.siemens.com/news/geometric-constraint-solving-and-spin-surfaces/
- D25 约束求解导论（含授权方名单）：https://blogs.sw.siemens.com/plm-components/geometric-constraint-solving-1-introduction/
- D27 2D DCM 集成封装：https://www.digitalengineering247.com/article/d-cubed-2d-dcm-integration-wrapper-available-from-siemens-plm-software
- D28 中望授权 D-Cubed：https://www.digitalengineering247.com/article/zwcad-software-licenses-parasolid-and-additional-d-cubed-components-from-siemens-plm-software/products
- L1 LEDAS 3D 建模：https://ledas.com/expertise/3d-modeling/
- L2 LEDAS 的五个求解器：https://www.ledas.com/post/848-ledas-has-developed-five-constraint-solvers/
- L3 LGS V9：https://ledas.com/news/201-bricsys-released-version-9-of-lgs-technology-constraint-driven-parametric-design-and-optimization
- L4 LGS V8：https://digitalengineering247.com/article/bricsys-announces-version-8-of-bricsys-lgs-2d-and-3d-component-technology
- L6 LGS 3D v4：https://ledas.com/news/145-ledas-releases-lgs-3d-v4-constraint-solver-for-assembly-design-and-direct-modeling/
- L7 LGS 2D v5：https://ledas.com/news/142-ledas-releases-lgs-2d-v5-constraint-solver-with-extended-functions/
- L8 LGS 2D v4.0：https://ledas.com/news/139-ledas-adds-advanced-parametric-drawing-functions-to-version-40-of-its-2d-constraint-solver/
- L10 LGS 2D 更新：https://ledas.com/news/125-ledas-enhances-parametric-drawings-with-updated-2d-geometric-constraint-solver/
- L11 LGS 3D v1.5：https://ledas.com/news/95-ledas-helps-cadcamcae-developers-efficiently-model-large-assemblies-with-lgs-3d-geometric-solver-version-15/
- L12 LGS 2D 2.0：https://ledas.com/news/59-ledas-opens-new-direction-of-geometric-solvers-development-with-release-lgs-2d-20/
- C1 C3D Solver 产品页：http://c3dlabs.com/products/c3d-toolkit/solver/
- C2 C3D 2D API：https://c3d.ascon.net/doc/math/group___constraints2_d___a_p_i.html
- C3 C3D 3D API：https://c3d.ascon.net/doc/math/group___g_c_m__3_d___a_p_i.html
- C5 C3D DOF 分析改进：https://c3dlabs.com/blog/products/dof-analysis-improvements-and-other-new-features-in-c3d-solver/
- C6 C3D 阵列与装配改进：https://c3dlabs.com/blog/products/c3d-solver-parametric-drafting-of-2d-patterns-and-3d-assembly-modeling-improvements/
- C7 C3D 偏移曲线与增量约束：https://c3dlabs.com/blog/tech-tips/c3d-solver-developments-offset-curves-and-increment-constraints/
- C8 C3D Toolkit：https://en.wikipedia.org/wiki/C3D_Toolkit
- S1 Spatial CDS 产品页：https://www.spatial.com/products/cds-constraint-design-solver
- S2 Spatial CDS 方案页：https://www.spatial.com/solutions/3d-modeling/constraint-design-solver
- A1 Autodesk VCS（招聘页中的描述）：https://builtin.com/job/principal-engineer/4068299
- A2 / R13 AutoConstrain 论文：https://arxiv.org/abs/2504.13178
- A3 Fusion AutoConstrain：https://www.autodesk.com/products/fusion-360/blog/autoconstrain-for-fusion-is-here/
- H1 CrownCAD DCS（第三方文章）：https://blog.csdn.net/weixin_43911798/article/details/129500917
- Z1 中望 ZGS：https://www.zwsoft.cn/news/16-13855.html

**算法与文献**

- R1 Zou 等，几何约束求解综述：https://arxiv.org/abs/2202.13795
- R2 Sitharam 等，手册章节（含 Bouma–Fudos–Hoffmann 的选根问题）：https://arxiv.org/abs/1608.05205
- R3 Hoffmann–Lomonosov–Sitharam，DR-planner：https://doi.org/10.1006/jsco.2000.0403
- R4 Ait-Aoudia–Jégou–Michelucci，二部图分解：https://arxiv.org/abs/1405.6131
- R5 3D 刚性的组合刻画问题：https://arxiv.org/abs/1210.0451
- R6 Michelucci–Foufou，见证构形法：https://doi.org/10.1016/j.cad.2006.01.005
- R7 SolveSpace 分支中关于正规方程条件数的 PR 讨论（仅搜索摘要）
- R8 尺寸参数的临界值：https://arxiv.org/abs/2206.08698
- R9 Wojtyra–Frączek，冗余约束的处理：https://doi.org/10.1115/1.4006958
- R10 ezpz（同 E1）
- R11 arael-sketch（同 X10）
- R12 Vitruvion：https://arxiv.org/abs/2109.14124

---

## 8. P0 完成情况（2026-11-13）

### 8.1 逐项状态

| 编号 | 状态 | 做了什么 |
|---|---|---|
| P0.1 | 完成 | 加 `.gitattributes`（LF）；清单测试改为忽略换行符差异后比较；全库 `cargo fmt`（20 个文件，纯格式）；clippy 存量清零（库 42 处、测试 4 处）。`needless_range_loop` 在 `Cargo.toml` 里统一放行，理由写在旁边 |
| P0.2 | 完成 | README、ROADMAP、`src/lib.rs` 与各模块文档按实现重写；CHANGELOG 新增 Unreleased 一节 |
| P0.3 | 部分 | criterion 基准（`benches/solve.rs`）和带基线对照的 `benches/BASELINE.md` 已入库。**CI 回归门没做**：共享 runner 上的墙钟时间噪声大，稳妥的做法是同一个 job 里对比基线与当前提交，而这套 workflow 我在本机无法验证 |
| P0.4 | 部分 | 刚体求解与诊断只为被约束引用的刚体建列；FD 行不再逐列克隆全部位姿；`JᵀJ` 只按非零元累加。**接地刚体消去没做**：它会改变 `ConstraintSolver::components()` 的分组语义，属于 P1 的破坏性变更 |
| P0.5 | 未做 | 被环境挡住：planegcs 对拍驱动在父仓库 `cad-rs/rcad` 里，克隆需要凭据；本机也没有 Python 和 FreeCAD，写了也无法运行验证。已确认先跳过，等另找到机器再做 |
| P0.6 | 完成 | 失败回滚（按连通分量，默认开启，`FailurePolicy::KeepLast` 可关）；`SolveStatus`；2D / 3D 的 `SolveOptions` 与各 `*_with` 入口；没通过线搜索的步不再被采用。线搜索的接受准则本身没动，仍允许残差相对增大 1e-4 |
| P0.7 | 完成 | `dof()` 改为秩口径，计数式保留为 `structural_dof_estimate()`；圆柱规范自由度；触及参数集改为排序去重 |

### 8.2 验收对照

| 验收项 | 结果 |
|---|---|
| 三平台 CI 全绿 | 本机（Windows）上 CI 的全部命令通过：fmt、clippy、两种特性组合的测试、doc、wasm32 检查。Linux / macOS 和 GitHub 上的实际运行**未验证** |
| 基准报告入库 | 完成，见 `benches/BASELINE.md` |
| 装配语料 debug < 60 秒 | **未达到**：411 秒 → 69 秒。release 下 21.0 秒 → 1.7 秒。2D 语料 debug 下 6.2 秒 → 2.4 秒 |
| 现有语料红线不变 | 达到。2D 与装配两套语料的指标与基线逐项相同（含回环次数、最差漂移、大幅漂移翻面数） |

基准摘要（同一台机器，release，中位数）：

| 用例 | 0.1.0 | P0 |
|---|---:|---:|
| 2D 单连通分量冷启动，198 参数 | 68.9 ms | 7.54 ms |
| 2D 单连通分量冷启动，792 参数 | 4.35 s | 357 ms |
| 2D 200 个独立矩形，3600 参数 | 13.7 ms | 8.66 ms |
| 2D 诊断，792 参数 | 737 ms | 767 ms |
| 装配星形 40 个刚体，求解 | 61.1 ms | 4.28 ms |

### 8.3 过程中的新发现

1. **411 秒的主因判断错了。** 我原先归因于“子系统按全体刚体建列”（1.3 第 5 条）。实测只改这一处，release 仅从 21 秒降到 16 秒。真正的主因是 `JᵀJ` 用三重循环稠密累加，而每行其实只有十来个非零元；改成只遍历非零元后降到 1.7 秒，结果逐位不变。2D 求解因此也快了 4–12 倍。
2. **“线搜索不接受变差的步”不能直接照字面实现。** 第一版在减半 8 次都失败时直接判停，结果 2D 语料的大幅漂移探针多了 7 例不收敛（基线 0 例）。旧代码是靠“硬接受一个变差的小步”碰巧逃出去的。最终做法是失败时把 λ 按 `JᵀJ` 的尺度逐级放大重算步长，仍不行才判停；这样没通过线搜索的步不会被采用，探针结果也与基线完全一致。这其实是 P2.3 信赖域的一个简化前身。
3. **单元测试 `large_system_makes_progress` 本身写错了。** 实体布局用了 `param_start = i` 而不是 `2 * i`，相邻点共用坐标，系统不相容，所以它才“只能有进展、不能收敛”。修正布局后完全收敛，测试改名为 `grid_system_converges`。P2 验收里“把这个测试收紧为完全收敛”一条已经提前满足，但原因不是内核变强了。
4. **语料测试的 JSON 报告从来没写出来过。** 代码在运行时读 `CARGO_TARGET_TMPDIR`，而 cargo 只在编译期设置它。已改为 `option_env!`。
5. **没有自由参数的分量会被报告为已收敛。** 全部参数被固定而方程不满足时，0.1.0 返回 `converged = true`。`space3d` 里有一个单元测试（`sphere_tangent_external`）明确断言了这个行为，所以我先保留并提请确认；确认后已改为如实报告（`converged = false`，状态 `Stalled`），该测试的断言也相应更新。
6. **高斯消元遇到 NaN 会 panic。** 选主元时对 `partial_cmp` 的结果 `unwrap`。已改为 `total_cmp`，非有限值走 `Diverged` 状态。
7. **热启动重解几乎不比冷启动快。** 离解 1e-4 处出发与离解 0.05 处出发耗时相同。推测是收敛后的 polish 还要做最多 4 次完整分解，未做剖析确认。这是 P5.6 帧预算的直接输入。
8. **诊断成了较慢的一半。** 792 参数的分量，诊断 767 ms，是求解（357 ms）的两倍。稠密列主元 QR 这次没动，P4 / P2 需要处理。
9. **带冗余约束时诊断会慢到不可用——已修。** 给基准补 3D 用例时发现：诊断为每一条相关行都重建并重新消元同一个 Gram 矩阵。780 参数、255 条冗余行的系统，一次诊断约 113 秒。改为只建一次、只分解一次后是 2.4 秒，结果不变。真实草图里冗余约束很常见，这个问题在 0.1.0 上一直存在，只是语料里没有“大而冗余”的用例所以没暴露。
10. **找到一个现成的病态难例。** 最小刚性的串联桁架（每个点只挂在前三个点上）从 0.05 的扰动出发：20 个点 5 次迭代，120 个点 17 次，160 和 260 个点 100 次内不收敛；每个点加一根冗余斜撑后，各种规模都是 4–5 次。这是误差沿链放大的条件数问题，适合作为 P2 信赖域 / QR 步长的验收用例，建议届时加进语料。

### 8.4 留给后续阶段的事

- **P1：** 接地刚体消去；刚体 `solve()` 改为逐连通分量求解（现在整体求解，星形用例里地面把所有刚体耦合成一个稠密块）。这两项做完，装配语料的 debug 时间才会明显低于 60 秒。
- **P2：** 稀疏分解。单连通分量仍是三次方复杂度：参数 4 倍，耗时约 47 倍。
- **P0 遗留：** CI 基准回归门；planegcs / SolveSpace 对拍驱动（需要父仓库的访问权限和一台装了 FreeCAD 的机器）。

### 8.5 对宿主有影响的行为变化

rcad-sketcher / rcad-assembly 升级时需要留意，细节在 CHANGELOG 的 Unreleased 一节：

- 求解失败后参数回到求解前（原来停在最后一次迭代）；
- `dof()` 不再出现负数，过约束要看 `diagnose()`；
- `SolveResult` 多了 `status` 字段，`DofReport` 多了 `rank` 和 `structural_dof`；
- 参数全部固定且方程不满足的分量，`converged` 由 `true` 变为 `false`；
- 自由圆柱的 DOF 由 7 变为 5。

按 Cargo 对 0.x 的语义化版本约定，这些属于不兼容变更。已确认：版本号定为 0.2.0，P0 与 P1 合在这一版里发布，第 5 节的版本表已相应调整。


---

## 9. P1 进度

分支 `p1-unified-core`，叠在 `p0-baseline` 之上。每一步单独提交，并以语料报告与上一步逐项对比作为验收。

| 步骤 | 状态 | 说明 |
|---|---|---|
| P1.1a 2D 与 3D 几何两族合并到一条流水线 | 完成 | 新增 `src/system.rs`：分量分解、逐分量求解、秩诊断各只有一份实现，按 `Equations` trait 编写。2D 和 3D 各自只剩一个适配层（提供残差与 Jacobian 行；3D 另提供平面、圆柱的规范自由度）。对外 API 不变，只有 `SpaceDiagnosis` 变成了 `Diagnosis` 的类型别名。两套语料的报告与 P0 完全一致，包括全精度的浮点数。删掉约 970 行重复代码，新增约 530 行通用代码 |
| P1.1b 刚体族接入同一条流水线 | 进行中 | **诊断已接入**：`diagnose_assembly` 现在就是通用诊断跑在刚体方程上（每个刚体 7 个存储数、6 个列），不再是第三份 QR 分析；`AssemblyDiagnosis` 多了 `rank` 和 `dependency_groups`。两套语料报告仍与 P0 完全一致。**求解循环还没动**：它仍是独立的一套 Gauss–Newton。换过去之后刚体才有线搜索、λ 升级、回滚和状态码，但迭代路径会变，验收要看装配语料的各项门槛而不是逐位相等 |
| P1.2 流形参数化 | 部分 | 内核已经通过“步空间”取步：普通参数做加法，刚体位姿做回缩（平移相加、四元数旋转后归一化）。平面还没有改成流形（仍是原来的归一化做法）。圆柱不适合改——`Concentric3D` 把轴上点当作圆心用，沿轴滑动不是纯粹的规范自由度，除非同时重新定义这个约束 |
| P1.3 稳定句柄与删除 | 未开始 | |
| P1.4 尺寸表与 `set_dimension` | 未开始 | |
| P1.5 刚体约束类型化、满秩行、serde | 未开始 | |
| P1.6 下线两套启发式分析器 | 未开始 | |
| P1.7 `validate()` 错误模型 | 未开始 | |
