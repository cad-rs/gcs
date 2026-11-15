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
12. **样条上点的残差不光滑。** `src/space3d/constraint.rs:2472` 在残差内部做 20 点采样 + 牛顿求最近点；标准做法是把曲线参数 `t` 作为辅助未知量。（P1 期间补充：这一族约束的导数原先也是错的，而且没有任何测试覆盖，见 §9 的“样条约束”一行。）
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
| 装配语料 debug < 60 秒 | P0 结束时**未达到**（411 秒 → 69 秒；release 下 21.0 秒 → 1.7 秒；2D 语料 debug 下 6.2 秒 → 2.4 秒）。P1.1 做完接地处理后是 38 秒，**已达到** |
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

- **P1：** 接地刚体消去；刚体 `solve()` 改为逐连通分量求解。**这两项已在 P1.1b 做完**（见第 9 节），装配语料的 debug 时间降到 38 秒。
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
| P1.1b 刚体族接入同一条流水线 | 完成 | 分两次提交（诊断一次，求解与接地一次）。**诊断**：`diagnose_assembly` 就是通用诊断跑在刚体方程上（每个刚体 7 个存储数、6 个列），`AssemblyDiagnosis` 多了 `rank` 和 `dependency_groups`；语料报告与 P0 完全一致。**求解**：刚体自己的那套 Gauss–Newton 循环删掉了，改走通用的逐分量求解，因此有了线搜索和 λ 升级；只做这一步时（加接地之前实测）语料报告仍与 P0 逐位相同。**接地**：带完整位姿 `Fixed` 的刚体直接放到目标位姿、不再作为未知量，挂在地面上的各个零件各自独立求解。这一步改变了数值结果，方向是变好：装配语料最差回环漂移 1.2e-7 → 1.0e-13，大幅漂移后翻面的 9 例变成 0 例；语料 debug 耗时 411 秒 → 38 秒（P0 时是 69 秒）。代价是通用流水线的固定开销：5 个刚体的小装配求解慢了 1.2–1.7 倍（109 µs 对 66 µs），40 个刚体的链持平，挂在地面上的 40 件星形快了 6 倍；诊断慢了 1.4–2.3 倍，但现在多报依赖组。刚体求解失败时仍返回最后的位姿（没有改成回滚），状态码也还没有暴露到 `PoseSolveResult`，这两点留给统一选项 / 结果类型时一起做 |
| P1.2 流形参数化 | 部分 | 内核已经通过“步空间”取步：普通参数做加法，刚体位姿做回缩（平移相加、四元数旋转后归一化）。平面还没有改成流形（仍是原来的归一化做法）。圆柱不适合改——`Concentric3D` 把轴上点当作圆心用，沿轴滑动不是纯粹的规范自由度，除非同时重新定义这个约束 |
| P1.3 稳定句柄与删除 | 2D、3D 完成 | 做法比原计划轻：不换 ID 类型，而是保证下标永不移位——删除只留下一个标记为“已删除”的槽位。`Sketch` 和 `SpaceSketch` 都有 `remove_constraint`、`set_constraint_enabled`（抑制）、`replace_constraint`、`remove_entity`（连带删除引用它的约束）。被抑制或删除的约束在流水线里贡献 0 个方程，所以分解、求解、诊断、可动性都自动忽略它们。`add_constraint` 现在返回下标。3D 里删掉的平面、圆柱连同它们的规范自由度一起退出 DOF 计数。`SpaceSketch` 的 `params` / `entities` / `constraints` 仍是公开字段（宿主在直接用），绕过方法直接改它们不会更新 `topology_revision`，文档里写明了。还没有回收槽位的 `compact()`；刚体 `ConstraintSystem` 的同一套接口见下面“刚体部分”一行|
| P1.4 尺寸表与 `set_dimension` | 2D、3D 部分完成 | 3D 与 2D 同一做法：`SpaceSketch::set_constraint_value` / `constraint_value` / `measure`。3D 的角度约束比较的是余弦（线–面是正弦），本身分不清一个角和它的镜像，所以 `measure` 返回主值（线–线、面–面在 `[0, π]`，线–面在 `[0, π/2]`）。**与原计划不同**：没有把数值搬到独立的尺寸表（那要改 `Constraint` 枚举的每个带值变体，宿主的翻译层全部要跟着改），而是原地修改约束里的值：`set_constraint_value` / `constraint_value`。功能目标都达到了——改尺寸不重建、不改拓扑（`topology_revision` 不变，结构缓存可以复用）；从动尺寸 = 抑制 + `measure`。没达到的是“多个约束共享一个值”和“值绑定到变量 / 表达式”，这两样留给 P6.5，届时在上面加一层尺寸表。持久性测试也因此有了真正热启动的尺寸往返用例（`tests/editing.rs`），补上了 Q6 |
| P1.5 刚体约束类型化、满秩行、serde | 完成 | 分两次提交。**类型化**：`rigid3d::Constraint` 从 11 个字段（多数是 `Option`、含义随类型而变，方向向量存在名叫 `point` 的字段里）改成 `{ id, kind }`，`kind` 是带数据的 `ConstraintKind` 枚举；构造函数签名不变，所以只用构造函数和求解器的宿主代码不用改，读字段的地方改用 `constraint_type()` / `component1()` / `component2()` 或直接匹配 `kind`。同时派生 serde，并把 glam 的 serde 支持改成跟随本库的 `serde` 特性（原来是无条件开启，默认构建白带一个 serde）。这一步是纯重构：把装配语料 333 例的全部求解位姿和两次诊断逐位导出，前后字节完全相同。**满秩行**：`Parallel` 从 3 行改成 2 行、`Concentric` 从 6 行改成 4 行，做法是在刚体 1 的方向的垂直平面上取一组随刚体转动的正交基 `e₁, e₂`，方向条件写成 `e_k · d₂ = 0`，同轴的偏移条件写成 `e_k · (p₂ − p₁) = 0`。新旧两种写法的残差平方和相同（都是 `sin²θ` 和垂直偏移的平方），所以误差数值和容差含义不变。**这修掉了一个此前没有记录的诊断缺陷**：改之前 `diagnose_assembly` 把每一个 `Parallel` / `Concentric` 配合都报成冗余，哪怕装配里只有它一个；单独一个 `Concentric` 还会把对方刚体的 `Fixed` 一并列进冗余表。新增 `tests/assembly_rank.rs`：每一种配合单独存在时，秩必须等于它拿掉的自由度数、不得报任何冗余或依赖组；重复的配合仍然要报。装配语料的红线不变（收敛 1.0、584 次回环、最差漂移 1.02e-13、翻面 0）。**没做成的一件事**：本想在装配语料里加一条“冗余精确率”红线，量出来才发现语料的“可解”用例本身就是有意过定的（三对重合点钉住一个刚体，9 个相容方程对 6 个自由度），诊断在那里报冗余是对的，这条指标在现有语料上没有意义，已撤掉；精确率改由上面的逐类型测试保证 |
| 刚体选项 / 结果与另两族对齐（P1.1b 遗留） | 完成 | 刚体求解失败时默认也回滚：没收敛的那个簇的刚体回到求解开始时的位姿（热启动位姿，或 `solve` 的种子位姿），收敛了的簇保留结果，接地刚体始终在固定位姿；`PoseSolveOptions::on_failure = KeepLast` 可以要回原来的“最后迭代点”。`PoseSolveResult` 新增 `status`。选项类型改名 `PoseSolveOptions`（原名 `rigid3d::SolveOptions` 留作弃用别名，它和另两族的 `rcad_gcs::SolveOptions` 同名太容易混）。三族现在共用词汇（`FailurePolicy`、`SolveStatus`、默认回滚），但类型没有合成一个：刚体的容差是 L2 范数而不是 RMS，结果里带的是位姿，硬合并只会多出一堆对某一族无意义的字段 |
| P1.3 / P1.4 / P1.7 的刚体部分 | 完成 | `rigid3d::ConstraintSystem` 有了和两个草图容器同一套编辑接口：`remove_constraint`、`set_constraint_enabled`、`replace_constraint`、`set_constraint_value` / `constraint_value`、`measure`、`topology_revision`、`validate`，以及 `add_component` / `remove_component`。与草图的两点不同：装配系统本身不存位姿，所以 `measure` 要调用方把位姿传进来，改完尺寸后用 `solve_from` 从上一次的结果接着解；`validate` 返回装配自己的 `AssemblyError`（零部件不存在 / 已删除、数值非有限、方向或旋转长度为零、约束两端是同一个零部件）。被删除的零部件不再求解、不计自由度，但在位姿数组里的位置保留，下标不移位。测试里也如实记下了一个已知限制：两个方向起始时平行的 `Angle` 约束没有梯度可循（§1.3 第 4 条），所以用例从“门虚掩着”的位姿热启动。至此三族都能原地编辑；`compact()`（回收已删除的槽位并返回下标映射）三族都还没有 |
| P1.6 下线两套启发式分析器 | 完成 | `DOFAnalyzer` 和 `SolverDiagnostics` 名字保留、标记弃用，内部改成秩诊断之上的封装：先从种子位姿把装配解出来，再在解出的位姿上读诊断（解不出来的簇因为默认回滚而停在种子位姿，冲突正好在那里可见）。于是三套口径不再互相矛盾。`DOFResult` 的总自由度是秩口径；每个刚体的“自由 DOF”改由新增的 `assembly_mobility`（逐刚体、逐坐标的可动性，复用 2D 的零空间支撑分析）给出。`SolverDiagnostics` 按依赖组报冗余、按诊断给出的 id 报冲突、报剩余自由度时列出还能动的零部件；不再报“过约束”（计数口径）和“循环依赖”（闭环配合是正常装配，旧实现把它当缺陷）。顺手导出了 `DiagnosticSeverity` 等三个类型——它们是公开字段的类型，却一直没法在库外写出名字。**两点如实说明**：① 这两个封装现在各要一次求解加一次分解，原来只是扫一遍约束表，宿主若在每帧调用需要注意；② `component_dof` 数的是“能动的坐标个数”，不是“独立运动的个数”（沿对角线滑动的刚体有两个平移坐标可动、但只有一个自由度），后者要对每个刚体各做一次秩计算，留给 P4 的逐实体诊断 |
| P1.7 `validate()` 错误模型 | 2D、3D 完成 | `Sketch::validate()` 和 `SpaceSketch::validate()` 返回类型化的 `SketchError` 列表：引用了不存在或已删除的实体、实体类型不对（给直线加半径、取圆的起点）、尺寸值 / 固定目标 / 方向向量或参数不是有限数；3D 另外检查样条控制点下标是否存在。`SketchError` 对实体类型做了泛型，两族共用一个枚举（3D 是 `SpaceSketchError`）。依据是每个约束新增的 `references()`，并有一条测试保证它和 Jacobian 实际触及的实体逐个变体一致。残差函数里原有的 `bool` 返回值没有动。刚体族见下面“刚体部分”一行 |
| 样条约束（计划外） | 已修 | 给 3D 加 `references()` 的一致性测试时，把此前**没有任何测试覆盖**的四个样条约束放进了“解析 Jacobian 对有限差分”的逐变体测试，查出两个原有缺陷：① `SplineControlPointSymmetric` 对镜像直线完全没有导数，对控制点的导数也不全——直线自由时它既不会动，也不会和样条分到同一个连通分量；② `PointOnSpline` 求导时把曲线上的最近点当成不动，只有点已经在曲线上时才对。debug 构建在每次求解开头会拿解析 Jacobian 和有限差分对比、不一致就 panic，所以点不在曲线上时 debug 构建直接 panic（实测到），release 构建则是带着不准的行迭代。现在 ① 复用 `Symmetric` 的行（方程本来就一样），② 改成精确导数：最近点在曲线内部时跟着参数移动，点落在曲线端点之外时最近点停在端点。另外修了零控制点样条在 debug 下的下标下溢。四个样条约束都进了逐变体测试，并补了求解测试。§1.3 第 12 条说的“残差内部做最近点搜索”这个结构问题没有变，仍留给 P6 |

### 9.1 P1 验收对照（2026-11-13）

P1 的验收条件是四条（见 §4 的 P1 小节），逐条对照：

| 验收条件 | 结果 |
|---|---|
| 三族经同一个入口求解与诊断 | 达到。分解、逐分量求解、秩诊断在 `src/system.rs` 各只有一份，2D、3D 几何、刚体各是一个提供残差和 Jacobian 行的适配层；刚体自己的 Gauss–Newton 循环和它那份 QR 分析已删除 |
| 重复的分量 / 诊断代码删除 | 达到。另外两套启发式的装配分析器改成了秩诊断之上的封装并标记弃用 |
| 三套语料红线全部保持 | 达到。2D：320 例，收敛 1.0，512 次回环，最差漂移 1.59e-4（与 0.1.0 相同）。装配：333 例，收敛 1.0，584 次回环，最差漂移 1.02e-13（P0 时 1.18e-7），大漂移翻面 0（P0 时 9）。3D 几何语料全部通过 |
| 附一份宿主迁移说明 | 达到：[MIGRATION_0.2.md](MIGRATION_0.2.md) |

P1 里**没有做完**的，如实列出：

1. **P1.2 平面没有改成流形。** 平面仍用“求解前后归一化 + 诊断时去掉径向列”的做法；它工作正常，只是比流形写法多一处特例。圆柱按上表的理由不改。
2. **P1.4 没有独立的尺寸表。** 尺寸值原地改；“多个约束共享一个值”和“值绑定到表达式”留给 P6.5。
3. **`compact()` 三族都没有。** 删除只留槽位，长时间编辑后槽位只增不减。
4. **P0 的两项遗留仍在**：CI 里的基准回归门（本机无法验证工作流）、与 planegcs / OndselSolver 的对拍驱动（按你的安排等另一台机器）。
5. **GitHub 上的 CI 从未跑过这些提交。** CI 只在 `main` 的推送和 PR 上触发，而我只推了分支。本机上每个提交都过了与 CI 相同的七个步骤（fmt、clippy、两种特性配置的测试、文档、wasm、基准构建），但 Linux / macOS 上没有验证过。开一个 PR 就能补上这一步。
6. **刚体基准在 P1 末尾重测过，与 P1.1 的数字持平**（各项相差不到 10%，有升有降，结果记在 `benches/BASELINE.md`）。这不说明后面的改动没有性能影响，只说明现有夹具量不到：链和星形夹具只用了重合与角度两种配合，行数变少的平行 / 同轴不在其中，成功的求解也不走回滚路径。要量这两项，得先加含平行 / 同轴配合的夹具。我先前在这里写的“之后的三处改动会影响它们”是没有量过的推测，已改正。

P1 期间在计划之外查出并修掉的原有缺陷有三处，都是靠“把此前没有测试覆盖的东西放进已有的逐变体测试”发现的：3D 样条约束的两处导数错误（其中一处让 debug 构建直接 panic）、装配诊断把每个平行 / 同轴配合都报成冗余、默认构建无条件带上 serde。

## 10. P2 进度

分支仍是 `p1-unified-core`（P2 的提交接在 P1 后面）。

| 步骤 | 状态 | 说明 |
|---|---|---|
| P2.6 刚体 Fixed 姿态行与 Lock 旋转行改解析 | 完成 | 这两类行原来用中心差分（步长 1e-7，每列求两次残差）。现在用旋转对数的导数：残差 `r = log(E)` 在右乘一个小旋转 `ε` 时变化 `A(r)·ε`，`A = I + ½[r]× + c(θ)[r]×²`，`c = 1/θ² − cot(θ/2)/(2θ)`（`θ → 0` 时取级数 `1/12 + θ²/720`）。`Fixed` 的列是 `A(r)·R_tgt`；`Lock` 对刚体 2 是 `A(r)·R_rel`，对刚体 1 是 `−A(r)·Eᵀ`。全库至此没有差分行，`FD_EPS` 只剩测试在用。新增测试把 `A(r)` 直接对差分核对，角度从 0 到 3.1 弧度，跨过级数与闭式的切换点；原有的“解析对差分”总测试里加了一个姿态误差不为零的 `Fixed`（原来只在误差为零处测到它）。基准持平：现有夹具里只有接地那一个约束带姿态行，量不出差别 |
| 姿态残差的两处缺陷（计划外，做 P2.6 时发现） | 已修 | ① **残差依赖四元数的符号。** 四元数 `q` 与 `−q` 是同一个旋转，但旧代码直接取误差四元数的旋转向量：误差角为 `a` 时，一种符号读出 `a`，另一种读出 `2π − a`。实测：离目标 0.001 弧度的 `Lock`，`w > 0` 时报误差 0.001，`w < 0` 时报 6.282。试过的求解仍能到达正确位姿，但对这种位姿报出的 `error` 是错的，`diagnose_assembly` 会把“用相反符号写的同一个姿态要求”报成冲突。宿主把矩阵转成四元数时符号是任意的，所以这不是构造出来的情形。现在先把误差四元数归到 `w ≥ 0` 一侧再取对数，旋转总是走短的一边。② **2e-8 弧度以下的姿态误差读成 0。** glam 的 `to_scaled_axis` 在向量部分小于 1e-8 时直接返回零，于是姿态残差在这个阈值处有一个跳变，再紧的容差也收不到阈值以下。现在自己做这个换算，任何大小都准确。两处都有先在旧代码上确认失败的测试（`tests/assembly_orientation.rs` 和 `solver.rs` 里的单元测试） |
| P2.3 全局化（信赖域） | 部分完成 | 步长控制从“减半线搜索 + 失败后加大 λ”换成 **Powell DogLeg 信赖域**。做法：每次迭代解一次正规方程得到 Gauss–Newton 步，再算出沿负梯度的 Cauchy 步；信赖半径起始为无穷大，所以好解的问题就是原来的整步 Gauss–Newton；一步按“增益比”（残差平方实际下降 ÷ 线性模型预期下降）评判，没有下降的步不采纳，半径减半，下一次试探沿折线向最速下降方向弯——拒绝一步只多算一次残差，不用重新分解。**先量后改**：动手前在串联桁架上量了基线，并做了两个候选（带增益比的 Levenberg–Marquardt、DogLeg）用临时开关对比，数据见下表，据此选了 DogLeg。没做的部分：LM / BFGS 回退链（DogLeg 单独已满足现有指标，暂不需要）、P2.4 的尺度化（半径目前按原始参数向量的范数量，长度和角度同权） |
| P2.4 尺度 | 2D 完成（大坐标、小几何的求解与诊断）；3D 草图与刚体的小几何未做 | **问题**：容差是绝对数（默认 1e-10），而坐标在十亿量级的点只能存到约 1e-7，方程不可能满足到 1e-10，所以远离原点的模型此前一律不收敛；debug 构建还会在大坐标或很小的几何上因“解析对差分”的核对而直接 panic。**做法**：每个方程按它自身数值的舍入精度来要求。取 `size = Σ|∂f/∂q|·|q|`（方程各项的量级），权重 `w = max(1, 4·ε·size/容差)`，残差和 Jacobian 的这一行都除以 `w`；原点附近的模型所有权重都是 1，走的就是原来的路径。权重在第一次算 Jacobian 时顺手取（只看自由列，不多花计算），循环没收敛就结束时再全量重取一次（把被固定的参数也算进去，最多两次）。加权后，只被大权重行碰到的列在 `JᵀJ` 里会小到被正则项 λ 和消元的主元阈值当成奇异，所以把这些列按碰到它的最小权重放大，在放大后的坐标里求步长再换回来。结束时残差按原方程报告；按原方程没到容差、按权重到了的，状态是新的 `ConvergedToPrecision`，`converged` 为真。刚体姿态的列：平移列取位置坐标的大小，旋转列取 1。**走过的弯路**：先试过“逐个方程达到精度就算数”的退出判据，在互相耦合的冗余系统上不成立，才改成加权；加权后先后撞上 λ 淹没被加权的行、主元阈值误判奇异，最后落到列缩放；差分步长一度改成相对步长，结果偏移一大步长就超过了特征尺寸，改回固定小步，核对时算上舍入噪声，并在不一致时用半步长再量一次做外推。**结果**：`tests/scale.rs` 五个红线测试——同一个模型在尺寸 1e-3 到 1e6、偏移到尺寸的 1e9 倍时都收敛；远离原点的冲突仍报冲突；3D 草图和装配同样。2D 语料各项指标与改动前逐项相同；装配语料往返的最差软漂移从 1.02e-13 变为 1.57e-13（允许 4 个舍入单位时；试过 8 个，是 5.08e-13，所以取 4）。**基准**：与改动前的提交交替各跑四轮，中等规模的冷解慢约 3%（`assembly/solve_chain/40` 均值 4.33 → 4.44 ms，`sketch/solve_chain/198` 均值 6.59 → 6.79 ms），而同一个构建各轮之间本身就有约 5% 的波动，所以这个差只在均值上看得出来；小规模夹具持平。两个构建在这两个夹具上的消元次数、正规方程构建次数、残差求值次数逐项相同，把加权关掉差别也还在，所以不是加权本身的开销；原因没有查到，剩下的怀疑是主循环改写后编译出来的代码不同（数据和排查过程记在 `benches/BASELINE.md`），留到 P2.1 / P2.2 换掉分解时一并看。 **没做的（小几何那一半）**：特征尺寸小于约 1e-4 的几何仍然不对。2D 的角度、平行、垂直、等长是按“两个长度的乘积”写的，几何一小，残差先掉到容差以下而关系并没有成立。实测（3-4-5 三角形，`cargo test --test scale -- --ignored report_small_geometry`）：尺寸 1e-5 时报收敛，第三个顶点偏了尺寸的 0.12%；1e-6 时偏 0.19%；把容差收紧到 1e-16，1e-5 要 27 次迭代，1e-6 跑满 100 次不收敛（固定的 λ 把这些行淹没了）。修法是把这些关系改写到单位方向上，让每个方程要么是长度、要么是角度，再让容差随分量的特征尺寸缩放——这会改动 2D 的残差定义，语料指标要重新对一遍，所以单独做。原计划里的“最大范数加步长范数”判据也还没做。 **小几何第一步：残差改写（完成）**。2D 的平行、垂直、角度改成单位方向上的 `sin θ`、`cos θ`、`sin(angle − θ)`，等长改成 `|d1| − |d2|`；3D 的等长线同样改（3D 的平行、垂直、角度原本就是单位化的）。**走过的弯路**：第一版让三个角度关系直接是纯角度（弧度量纲），小尺寸全好了，但尺寸 1e6 的三角形不再收敛——角度行对坐标的导数变成 1/尺寸 = 1e-6，在 `JᵀJ` 里只有 1e-12，被固定的 λ = 1e-10 吞掉，正好是原缺陷的镜像。要让纯角度行和长度行共存，得让家族逐行、逐列报量纲并做双向缩放（圆弧的角度参数也是一列弧度）；改成按计划原文“角度类残差乘特征长度”更简单：乘两条线长度的几何平均 `√(|d1|·|d2|)`，于是 2D 的每个方程都是长度，导数在任何尺寸都是 1 的量级，不需要逐行分类。**结果**：新红线 `the_triangle_solves_exactly_at_any_size`（尺寸 1e-6 到 1e6、容差 1e-10 与 1e-16，第三个顶点相对误差 ≤ 1e-12，实测都在 2e-15 以内、2–3 次迭代）在旧代码上失败（1e-6 时偏 1.89e-3）；`a_3d_equal_length_holds_on_small_lines` 在旧代码上失败（1e-8 时报收敛但线短 10%；1e-6 时旧代码能过，λ 只是让它收敛变慢）。2D 语料（含角度族求解与平行 + 垂直冲突族；等长不在语料里）与装配语料各项指标逐项相同（71 / 1.59e-4 / 9.0；1.57e-13 / 0），3D 语料通过。基准夹具（矩形链、桁架、装配）不含这几种约束，没有做 A/B。**小几何第二步：容差随分量尺寸缩放（2D，完成）**。完全约束的模型靠收敛后的精修（不加正则的整步）到达舍入精度，所以第一步的红线能过；还有自由度的分量没有精修，只满足到绝对的 1e-10——实测（三角形去掉角度约束，可绕固定角转动，`report_small_geometry_free`）：尺寸 1 和 1e-2 时长度相对误差 ≤ 7e-16，1e-4 和 1e-6 时都是 2.7e-8。**做法**：家族通过新的 `Equations::length_scale` 报一个分量的尺寸（2D：约束所引用实体的包围盒长边，圆和圆弧按整圆算；3D 和刚体不报）；尺寸 `L < 1` 时把它作为方程的单位，每个方程的权重从 `max(1, 精度权重)` 改成 `max(L, 精度权重)`，即按 `残差 / L` 对容差。权重在第一次 Jacobian 之前就取 `L`（否则第一次检查会用绝对容差提前放行），第一次 Jacobian 时再换成含精度的权重，残差按新旧权重之比换算。列缩放取“碰到这一列的最小权重”，与均匀的行权重正好抵消，所以步本身不变，只是停的地方变了。`L ≥ 1` 时走原来的路径。顺带：全部参数被固定的分量，原来先按原始残差判“收敛”、只有不满足时才在“权重至少翻倍”的条件下看精度，现在直接按完整权重判（单位和精度都算进去）；`L ≥ 1` 时结论只在残差落在容差的 1–2 倍之间时可能不同。**结果**：新红线 `a_small_model_with_freedom_left_is_held_to_its_own_size`（尺寸 1 到 1e-8，原点与偏离尺寸 1e3 倍两处，长度相对误差 ≤ 1e-12 + 偏离处的舍入）在旧代码上失败（1e-4 时 2.7e-8）；现在 1e-8 到 1e3 都是 3 次迭代、误差 ≤ 1.4e-15。2D 与装配语料各项指标逐项相同。基准：夹具是 4 × 3 的矩形，单位都是 1，新增的只是每次求解一次包围盒；交替三轮在 −4% 到 +6.5% 之间，但这次机器漂移大，同一构建各轮之间差 10–20%，第一轮里成对比较都在 ±3% 以内，没有量出差别（数字在 `benches/BASELINE.md`）。**小几何第三步：2D 诊断（完成）**。`diagnose` 的阈值是 `CONSIST_TOL · max(1, max|f|)`（部分冗余用 `PARTIAL_TOL`），下限 1 是绝对数。先写测试：解好的三角形上加一个同样的距离（应为冗余）或多 0.2% 的距离（应为冲突），尺寸 1e-6 到 1e3，位置在原点、偏离尺寸的 1e3 倍和 1e9 倍。旧代码在 1e-6 时把冲突报成 `partially_redundant`。改法：下限从 1 换成与求解同一个单位 `min(1, 分量尺寸)`。**走过的弯路**：第一版还给阈值加了“4 个舍入单位”的余量，担心远离原点的小模型把舍入误差读成冲突；关掉余量再测，偏离 1e9 倍时也全过——测试里的冗余是逐字相同的副本，两行舍入后的残差完全一样，最小二乘的驱动项正好抵消——所以没有数据支持这个余量，去掉了。不是逐字副本的冗余在远处会不会被误判，没有测到，记在 `ROADMAP.md`。2D 与装配语料各项指标逐项相同。新增开销是诊断时每个分量一次包围盒，和上一步求解里的一样，上一步的 A/B 没有量出差别，这一步没有再测。**还没做**：3D 草图（角度类方程是纯数，和长度混在一起）和刚体（长度与角度混合）仍是绝对容差；圆弧的角度参数（弧度列）在尺度上还没看过 |
| P2.5 退化几何防护 | 部分完成（普查 + 修一处；起点退化时选方向的问题待定） | **先量**：`tests/degenerate.rs` 的 `report_degenerate` 普查 16 种退化起点（零长线段上的长度、点在线上、相切、平行、角度、等长、对称；重合点间距离；点在圆心上的点在圆上；同心圆外切；零半径圆和圆弧；NaN / Inf 坐标、NaN 尺寸、被固定的 NaN 坐标），debug 和 release 各跑一遍，另测 3D 草图（NaN 坐标、零长线段的长度）和刚体（NaN 起始位置）。**结果**：没有一例 panic，没有一例把 NaN 留在参数里，失败的求解都把几何放回原样，NaN / Inf 输入一律报 `Diverged`——计划里“NaN / Inf 检出并转成状态码”在 0.2 先前的改动里已经做到（迁移说明 §2.1）。**唯一的错状态**：零长线段上的 `LineLength`，导数 `dx / |d|` 没有防护，0/0 让求解立刻报 `Diverged`；改成和 `EqualLength` 共用的方向函数，零长时梯度为零，报 `Stalled`，与“重合点间距离”一致（测试先在旧代码上失败）。**看了但没改的**：① 2D 里按长度归一的几处用的是绝对下限 `SOLVER_NORM_TOL = 1e-10`，本想换成只防恰好为零，但先写的测试（1e-12、1e-14 长的线上的点）在旧代码上就能过，没有证据说明它有害，没改。② 起点退化（零长线段给长度、重合点给距离、点在圆心）梯度为零，求解报 `Stalled`；按固定方向取次梯度就能解出来，但这是在替宿主决定几何往哪走，debug 的解析对差分核对也要为它开例外。用户定（2026-11-13）：保持 `Stalled`，不取次梯度。③ 零长线段参与的平行、角度、点在线上、对称是“空满足”（残差为 0），改写前后都如此。④ **两个错答案**：同心的两圆要求外切，求解靠把半径变成 −1 报 `Converged`（没有任何东西让半径保持为正，需要先定半径符号的约定）；零半径圆弧在起点处做诊断，把可以解出的端点约束报成冲突（线性化在奇异点上的局限）。都记在 `ROADMAP.md` |
| P2.2 线性步 | 稠密部分完成（精修复用分解、`LuFactors` 提速、病态步按需走 QR）；稀疏大块未做 | **先量**：在 198 参数的矩形链上数分解次数，冷解（偏 0.05）和热重解（偏 1e-4）都是 2 次迭代、再加收敛后精修的 4 次完整分解——精修占了分解的三分之二，这就是热重解不比冷解便宜的原因（`benches/BASELINE.md` 原来写的猜测是对的）。**`LuFactors` 提速（完成）**：精修改成复用分解时，第一版 A/B 装配链慢了 65–101%；查下去是 `LuFactors::factor` 本身——和 `gaussian_elimination` 同样的算术，却慢 2–3.3 倍（嵌套 Vec 下标访问的边界检查、每列一个乘子 Vec）。改成行切片、乘子存在消元后的下三角里（LAPACK 的存法，后续行交换带着乘子走，求解时先做全部交换），与消元只差几个百分点，结果逐位相同（新增 40 × 40、二十多次行交换的逐位测试）。同样的切片写法用在 `gaussian_elimination` 上没有变快，没改。诊断基准快 4–6%（规模大的两个各轮波动 10–15%，没有量出差别）。**精修复用分解（完成）**：精修的第一步分解正规矩阵并保存，之后的步沿用同一个 Jacobian 和分解（弦迭代：已经收敛，Jacobian 几乎不变），每步只剩回代和一次残差求值；只走一步的精修不再为第二步准备稀疏行。A/B（所有基准，三轮交替，同一构建各轮相差几个百分点）：2D 单分量求解与热重解快 38–49%，20 个独立矩形快 31%，3D 桁架 240 / 780 参数快 36% / 14%；装配夹具（精修都只走一步）在 ±3% 内，例外是最小的装配链 +14%（0.133 → 0.151 ms）和最小的桁架 +8%。**为什么慢了说得清楚的部分**：用临时计数器数过，装配语料 12 592 次求解、这两个小夹具，改前改后的迭代次数和精修步数逐项相同，装配上的精修全是一步；又试了“第一步原样走旧路径（消元、不存分解），第二步才分解”的变体——这两个夹具照样慢（+17%、+9%），而 2D 只快 28–35%（第二步要多一次分解），所以没采用。慢的不是精修做的事，是函数变大后的编译结果，和 P2.4 留下的 3% 同一类，没有再查。语料各项指标逐项相同。**QR 步之前先找证据（完成，等用户定）**：QR 的浮点量约为现在的两倍，又丢掉 `JᵀJ` 的稀疏累加，现有夹具上多半更慢，所以先找正规方程确实失败的例子。找到的是极细的三角形（`tests/hard_cases.rs` 的 `report_needle_triangle`）：A、B 相距 1 且固定，C 到 A 为 2、到 B 为 `1 + δ`，C 离 AB 约 `2√δ`，两行几乎平行，`κ(J) ≈ 1/(2√δ)`，`δ = 0` 时两圆相切。实测分三件事：① 默认容差下 `δ ≤ 1e-10` 全报 `Converged`，C 的高度错到自身大小的量级（`δ = 1e-12` 时 1.7 倍）——相切附近残差小不代表位置准，这是收敛判据的问题（计划里没做的“步长范数”判据），和线性代数无关。② 收紧容差后 `δ = 1e-12`、`1e-14` 跑满 100 次迭代，而同样的起点用精确 2 × 2 求解的牛顿法 12–21 次就收敛；把 λ 从 1e-10 临时改成 1e-20，求解器在 `1e-12` 也 12–18 次收敛——主因是固定的 λ 吞掉了 `JᵀJ` 里特征值约 `4δ` 的方向。③ λ = 1e-20 时 `δ = 1e-14` 仍有一个起点 `Stalled`，精确 2 × 2 却能收敛——这是正规方程把条件数平方的代价。QR 步（不经 `JᵀJ`，正则作用在奇异值而不是它的平方上）同时针对 ② 和 ③；① 要靠步长判据。用户定（2026-11-13）：先做步长判据，再做 QR 步。**步长判据先量（结论与预期相反）**：用临时探针在“宣布收敛”的那一点记下不加正则的 Gauss–Newton 步（精修的第一步）的无穷范数。2D 语料 6807 次收敛里约 860 次这一步大于容差，但几乎都在 1–10 倍之间（矩形族 800 次：循环在第一个低于容差的迭代点停下，下一步本来就是残差的几倍，精修随后把它消掉，无害）；大于 100 倍的 59 次全在 `tangent-*` 族（20 个标为“可解”的相切例子），最大 7.9e-5；其余所有例子最大 1.1e-9。装配 12 018 次里 17 次大于容差，最大 196 倍；3D 语料没有。极细三角形默认容差下那些“错”的答案，这一步是 1e-8 到 1.0e-5——落在相切族的范围之内。也就是说，容差 1e-10 下 `δ ≤ 1e-10` 的极细三角形在数值上就是相切：C 在容差内满足两个距离，位置只能定到约 √容差，和任何相切一样。没有一个阈值能把它和语料里公认“已解”的相切例子分开，加步长判据就要把那 20 个相切例子判成不收敛。① 因此不是缺陷，而是残差容差的本义（相切附近位置精度约为 √(容差·尺寸)）。没有实现。用户定（2026-11-13）：不做步长判据，直接做 QR 步。**QR 步（完成，按需切换）**：写了完全正交分解求最小范数步（`min_norm_step`：按列存的列主元 Householder QR，范数降阶更新，按 `1e-12` 相对阈值截断秩；秩不满时对 `[R11 R12]ᵀ` 再做一次 QR 取最短解，同 LAPACK `xGELSY`），六个单元测试对已知解。先试“每一步都走 QR”：极细三角形 ②③ 都好了，2D 语料最差软漂移 1.59e-4 → 1.21e-5、装配 1.57e-13 → 1.47e-13，其余指标不变；但基准慢了 13%（18 参数）到 4 倍（780 参数）——QR 的算术量是稀疏累加 `JᵀJ` 再消元的 4–6 倍。用户选按需切换：先照旧分解 `JᵀJ + λI`（改用 `LuFactors`，与消元逐位相同），最小主元不到最大主元的 1e-8（秩亏或丢了一半有效位）才改走 QR。临时计数：2D 语料 18 317 步里 75 步走 QR，装配 37 781 步里 16 步，3D 语料 41 步里 17 步。结果：极细三角形 `δ = 1e-12`、`1e-14` 在紧容差下 12–20 次迭代收敛到几何的舍入极限（新红线 `a_needle_triangle_solves_at_a_tight_tolerance`，旧代码 100 次不收敛）；语料各项指标与改动前逐项相同——全 QR 版的漂移改进来自本来良态、按需版不会切换的步，所以没有保留下来；基准持平（±5%，个别离群都是单轮）。宿主可见：`LAMBDA`、`lm_lambda` 只作用在良态的步上（已写进迁移说明）。**没做的**：精修仍用正规方程（λ = 0 的 LU，复用分解），极病态时它可能被拒，但收敛后的精修只影响最后几位；大块的稀疏 QR（`faer`，特性 `sparse`） |
| P2.1 Jacobian 改 CSR、迭代内少分配 | 完成（CSR；正规矩阵与布局不再按行分配。没做到字面上的“零分配”，见末尾） | **先量**：在 `newton_solve_in` 各阶段临时加计时器，另用一个带计数分配器的临时 example 跑小夹具（都没提交）。54 参数矩形链每次求解 954 次分配，稠密 Jacobian 相关的几遍（清零、经 HashMap 填入、三次 `sparse_rows`、加权、`normal_equations` 对每行扫全部列）合计占三分之一，分解占 28%；240 参数桁架分解占 72%，稠密的几遍仍占 21%。本机一次分配加释放约 40 ns。**做法**：每个连通分量在 `solve_one` 里布局一次 CSR 模式（每个方程：其约束列出的自由列，升序去重），同时记下约束每列出一个条目对应的值槽；之后每次 Jacobian 只把导数加到槽里，不再查 HashMap（`system::Assembly`）。加权、线性模型、正规方程、QR 步、精修都直接读同一份行，`sparse_rows` 删掉；正规方程按行的非零元累加，与原来逐项、同序。模式里允许有恰好为 0 的条目，所有乘积跳过它们，和稠密矩阵跳过零一样。差分路径（公开的 `newton_solve`）用全稠密模式。残差闭包直接写进循环的缓冲区，不再新建 Vec 再拷贝。trait 约定约束列出的列与取值无关、顺序固定；条目顺序变了时按行二分查找放置。**验证**：临时在 `newton_solve_in` 末尾把状态、迭代次数、残差和全部参数的位模式写进文件，`--release --all-features --test-threads=1` 跑完整测试集：改前改后 21 437 次求解逐字节相同（含三套语料）。七步门禁全过，语料指标逐项相同。**A/B**（全部基准，三轮交替；改后第三轮所有基准包括不受影响的诊断都慢 10–50%，是机器，取三轮中位数）：单分量求解 18–60 参数快 13–29%（`space/solve_truss/60` −29%、`sketch/solve_chain/54` −18%、`sketch/solve_field/360` −19%），198–240 参数快 5–12%，780–792 参数快 3–6%（剩下的主要是稠密 LU，54 参数链上占 39%、240 参数桁架上占 83%）；装配快 3–6%；诊断 +0–6%，两边各轮重叠，代码没动。分配：54 参数链 954 → 644，240 参数桁架 8 707 → 3 068。**剩下的分配**（计数）：正规矩阵是每行一个 Vec，每次 n + 2 次；布局时每行一个 Vec；按每次 40 ns 估算，在最小的夹具上约值 3–10%。**第二个提交**：正规矩阵和它的 LU 改成一块按行存的缓冲区（`Square`，主元、更新、次序都不变）；布局直接把列写进模式，不再每行一个 Vec，并且和该次求解的所有 Jacobian 构建共用一组草稿行（原来每个约束都把草稿行丢掉重建，刚体族每行要扩容几次）。逐位指纹同样 21 437 次逐字节相同；七步门禁全过，语料指标不变。实测比估算大：A/B 三轮中位数，18–60 参数快 24–34%，装配（除 40 体链）快 23–31%，198–240 参数快 7–12%；最大的几个（40 体链、780–792 参数）中位数 +5% 到 −6%，各轮完全重叠（这一轮机器噪声 15–20%）。为了把大规模说清楚，单独计时了分解本身（同一个带状矩阵，新旧交替各 5 次，取中位数）：n = 18 / 54 / 200 / 790 时扁平存储快 57% / 32% / 25% / 3.3%，没有哪里更慢，所以大规模那几个 +3–5% 是机器。诊断的 Gram 系统也用这个分解，所以诊断基准跟着快 2–25%。只换扁平矩阵、布局照旧的快速对比是小夹具快 3–15%，其余的收益来自布局。分配：54 参数链 644 → 359，240 参数桁架 3 068 → 733，5 体装配链 613 → 311。**没做到“零分配”**：每次迭代仍新建正规矩阵（一次分配，`alloc_zeroed`，大矩阵时比复用后清零还便宜）、几个长度为 n 的向量（步、梯度、残差、保存的参数）。54 参数链每次求解还剩 359 次分配（含迭代外的：克隆、分量划分、布局），全按 40 ns 算上限约 9%，迭代内的只是其中一部分，没有单独量；240 参数桁架 733 次，可以忽略。复用要把工作区穿过 `LinearModel`、精修和 `StepSpace`，目前没有数据说明值得，暂不做（记进 ROADMAP）。|
| P2.7 按连通分量并行（可选特性 `rayon`） | 完成 | **先量**：`solve_one` 在 20 / 200 个独立矩形上占求解时间的 83–87%，在挂在接地件上的星形装配上约 63%，每个分量只有 5–14 µs。**第一版**（每个分量一个任务）实测是亏的：2–10 个小分量的系统慢 20–65 µs，20 个矩形持平，200 个矩形只快 30%（24 线程）；本机唤醒线程池里的一个线程要几十微秒，相当于好几个小分量的工作量。（先用带计数分配器的驱动量过一次更差，那是 24 个线程抢同一个计数器的缓存行，换成不计数的驱动重量。）**做法**：相邻分量按方程数均匀切成任务，每 96 个方程一个、不超过线程数；不到 192 个方程（凑不出两个任务）就照旧顺序求解。每个任务在自己的一份参数副本上依次求解它的分量，结束后把每个分量动过的存储参数（`StepSpace::storage`，刚体是 7 个存储量对 6 列）拷回。分量只读自己约束涉及的参数、只写自己的参数，所以结果与顺序求解逐位相同。门槛先后试过“每个任务至少 128 个方程、凑够两个就并行”（20 体星形被切成 135 + 45 两段，反而慢 13%）、再加“零头并进前一段”（20 个矩形变成 144 + 216，大的一段决定时间，几乎不快），最后改成均匀切分；任务大小试了 64 / 96 / 128，20 个矩形快 21% / 25% / 8%，其余相近，取 96。**线程安全**：`Equations` 只在开启特性时要求 `Sync`（`SharedAcrossThreads`）；刚体族每个约束的求值缓存开启特性时是 `Mutex`（分量之间不共享约束，从不争用），默认构建仍是 `RefCell`。先试过默认也用 `Mutex`，`assembly/solve_chain/5` 三轮都慢约 9%；但改回 `RefCell` 后的一轮里它仍是 +9.8%、另一轮持平，所以这个数和锁分不开，是这个夹具一贯的构建间波动。按特性区分的理由不靠它：默认构建本来就不需要锁。**验证**：强制所有多分量系统走并行跑完整测试集，21 437 次求解排序后与参照逐字节相同；按最终规则再跑一次也相同；默认构建按原顺序逐字节相同。新增三个测试：2D 和刚体各一个比较两条路径的整个参数向量（故意把拷回去掉，两个都失败），一个检查切分规则（24 线程和 2 线程都跑过）。七步门禁全过，语料指标不变。**A/B**（三轮，三种构建交替，中位数）：开启特性后 20 个矩形快 41%，200 个矩形快 64%（2.8 倍），40 体星形快 14%；门槛以下的小装配走顺序路径但仍有代价（5 体 / 20 体星形慢 14% / 7%，5 体链慢 5%，锁和更大的构建）；单分量求解 −1% 到 +6%，在各轮波动内。默认构建大多在 ±4% 内，例外是 5 体装配链一次 +9.8%、另一次持平，5 体装配诊断两次都 +8%，原因没查到（不是缓存布局，两种布局各量了三轮，结果重合）。**没做的**：单个大分量内部不并行（那要等稀疏分解或 P3 的结构分解）；门槛按本机 24 线程量出，线程少的机器上任务数自然更少，没有另外量。 |

**P2.3 的对比数据**（release 构建；“解算”是线性方程求解次数，“翻面”是 2D 语料大幅漂移探针 256 次里回不到原构形的次数）：

| | 原线搜索 | LM（候选） | DogLeg（采用） |
|---|---|---|---|
| 2D 语料红线 | 全过 | 全过 | 全过 |
| 2D 语料：解算次数 | 18 265 | 20 008（+9.5%） | 18 288（+0.1%） |
| 2D 语料：大幅漂移翻面 | 71 | 65 | 71 |
| 装配语料 | 全过，翻面 0 | 与原方案逐次相同 | 与原方案逐次相同 |
| 串联桁架，扰动 0.05，20 / 80 / 120 / 160 / 200 / 260 个点的迭代次数 | 5 / 9 / 17 / **不收敛** / 35 / **不收敛** | 5 / 21 / 25 / 30 / 31 / 53 | 5 / 9 / 22 / 19 / 34 / 84 |
| 串联桁架，扰动 0.2 | 20–200 点“收敛”到别的构形（位移 3–97），260 点不收敛 | 20 点回到原构形，80 点以上 100 次迭代内不收敛 | 20–120 点收敛到别的构形，160 点以上不收敛 |

怎么读这张表：

- **DogLeg 在语料上与原方案等价**（指标逐项相同，工作量差 0.1%），同时修好了扰动 0.05 时 160 点和 260 点两个不收敛的情形。装配语料上两者逐次相同，是因为那里从来没有一步被拒绝。
- **LM 更保守**：翻面少 6 次，扰动 0.2 的 20 点桁架能回到原构形。代价是每拒绝一步就要重新分解一次（解算多 9.5%），中等规模上迭代更多。没有选它，但“保守的步更能保住分支”这个观察对 P5 的拖拽和分支保持有用。
- **扰动 0.2 的长串联桁架是 P3 的题，不是 P2.3 的题。** 三种方案在这里都不好：原方案的“收敛”是跳到了几十个单位之外的另一个解，两个信赖域方案是迭代次数用完。这种系统每个点只依赖前三个点，按块三角顺序逐点求解就是 3 个未知量一组的小问题；全局最小二乘在这里是用错了工具。`tests/hard_cases.rs` 的报告留着，P3 做完后它应当全部收敛。
- 不可解的系统现在报 `Stalled` 并且很快结束（冲突草图不到 20 次迭代），原来要跑满迭代上限后报 `MaxIterations`。

## 11. P3 进度

分支仍是 `p1-unified-core`。

| 步骤 | 状态 | 说明 |
|---|---|---|
| P3.1 平凡方程消元（2D） | 完成（2D 求解；3D 草图、诊断、一般线性系数未做） | **先量**：临时探针（没提交）在 2D 语料 320 例上按约束逐行做带偏移的并查集，数“求解要迭代的列”（未知类里被非消元方程碰到的）：46 223 → 5 371，少 88.4%；方程 46 170 → 5 319。按例子平均少 78.6%、中位数 85.7%，320 例里 23 例不到 40%（角度族 10 例一个也消不掉）。分族：矩形 −88.9%（每个矩形 18 → 2：剩下右边的 x 和上边的 y，对两个长度）、边距 −85.7%、对称 −71.4%、圆链 −66.7%、相切 −50%、欠约束 −50.2%、角度 0%。基准夹具的矩形链 / 矩形阵列与矩形族相同；3D 桁架没有这类方程（留给 P3.2）。探针同时找出了别名图里的矛盾和重复：冲突族 20 例中 7 例的矛盾就在别名图里（半径对直径、重复固定点），冗余族 9 例中 6 例的重复在别名图里（重复的重合、一致的半径与直径）。**做法**（`src/alias.rs`）：`Equations::alias(约束, 行)` 让家族声明一行残差恰好是 `p − q − c` 或 `p − v`（2D：存储点之间的重合、固定、`DistanceX/Y`，水平 / 竖直，等半径，同心，半径，直径；圆弧的派生端点不算）。`solve_one` 对每个分量建并查集：被固定的参数一开始就是已知；一行若把两个不同的类连起来（不能两个都已知）就被吸收，否则（成环、或连两个已知类）留作普通方程——于是一致的重复在舍入内成立，矛盾留下动不了的残差，**没有任何方程被丢掉**。然后把参数放到类说的位置（已知类放到值上；未知类放到各成员所说根值的平均，也就是离当前几何最近的、让被吸收方程成立的点），布局时每个成员的列映射到它的类的列（`Assembly::lay_out`，导数在同一个槽里相加），被吸收的行留空（残差仍参与 RMS，所以容差的含义不变），步空间 `ClassSpace` 让一列移动整个类。约化求解不收敛就把参数放回原样，按全部参数重解（原路径），所以失败的分量结局与以前逐位相同。诊断仍看完整系统，冲突、冗余照常带原约束编号——计划里“别名图里的矛盾不能从报告里消失”由此满足。debug 构建的“解析对差分”核对仍在完整布局上做（在动参数之前）。**三个变体**（临时环境变量开关对比，用完已删）：未知类的起点取根参数现值或成员平均；列范数是否乘 √(成员数)（这样步长、λ、信赖半径与原来按全部参数量的度量完全一致，在一致的线性化上与原来的步在精确算术下相同）。语料上三者的最终几何逐字节相同；√k 度量在大幅漂移下迭代更多（2 977 → 3 940），浮动矩形上移动量小 0.1–0.7%、迭代 2–3 次而不是 1 次。用户定（2026-11-14）：采用、默认开启、取“平均、不乘度量”。**结果**：2D 语料红线全过（收敛、冲突召回 / 精确、冗余、欠约束、无迟滞都是 1.0，singular surfaced 10/10）；**选解变了**：大幅漂移翻面 71 → 5（矩形的镜像翻面全部消失），最差软漂移 1.59e-4 → 6.79e-5，最差硬漂移仍是 9.0；冲突族逐位不变；装配语料逐项相同（584 次回环、1.57e-13、翻面 0）。用一个不计入提交的探针逐例比较最终参数：完全约束例子从画好的状态出发，相切族 20 例差 7.6e-6 到 3.2e-5（相切处位置本来只能定到约 √容差，见 §10 P2.2），其余各族差 ≤ 1e-9，迭代合计 5 020 → 2 825；欠约束的开放距离链在漂移 0.05 后差到 1.6e-3（离起点的移动量相同），漂移 3 后离起点少移动 2–9%。`scale` / `degenerate` 报告只在舍入位上不同。**走过的弯路**：第一版把 NaN 坐标也“放到类说的位置”，一个固定点约束就把 NaN 修好、报 `Converged`，违反了“NaN 输入报 `Diverged`”的约定；改成起点有非有限值时不消元（新测试做了变异检查）。`failure_policy.rs` 里两个测试依赖“一个距离从一个固定点出发要好几次迭代”，消元后固定点先放好，剩下的一行沿步方向是线性的，一次就收敛——换成两段距离定一个顶点的夹具（在旧代码上同样通过）。**新测试** `tests/elimination.rs`：全由别名组成的系统 0 次迭代放好（旧代码 1 次，失败）；漂移 3 的矩形 16 个起点都回到画好的分支（旧代码第 8 个起点翻面，失败）；别名环上的矛盾不被“解掉”、参数放回原样、诊断点名闭环的约束（覆盖性；把“放回原样”去掉后它失败，变异检查）；一致的别名环收敛并报冗余、被消元的坐标严格相等、NaN 坐标仍报 `Diverged`（覆盖性；去掉非有限值检查后最后一个失败）。**A/B**（全部基准，与 `d1bf6e1` 交替三轮，中位数，各轮相差几个百分点）：矩形链 54 / 198 / 792 参数 0.137 → 0.041 ms（−70%）、3.12 → 0.142 ms（22 倍）、198 → 0.77 ms（257 倍），热重解相同；20 / 200 个独立矩形 −42% / −34%（剩下的是每个分量的固定开销）。夹具全是矩形，最占便宜：每个矩形 18 个参数只剩 2 个未知量，792 参数的链只剩 88 列。2D 诊断没变（−0.3% 到 −6.7%）。不消元的两族：3D 桁架 +1.5–2.9%、装配链 +3.7% / +4.7% / +9.5%（5 体），三轮都偏慢；它们多出来的只是每个分量问一遍“有没有别名”、检查一次参数是否有限、布局拆成两个函数，迭代内的代码相同，这部分没有单独计时，原因没查（同 P2.2、P2.7 里这些夹具的构建间差异）。**没做的**：3D 草图族（同样的重合、固定、轴向关系）；系数不是 ±1 的线性方程（中点、`SymmetricPoint`）和圆弧派生端点；诊断的约化；类在拓扑不变时的缓存（P3.4）；约化失败的分量要多解一次（冲突草图付两次求解）。都记在 `ROADMAP.md` |
| P3.2 结构分解与块三角求解 | 完成（第二次试做后用户定采用，2026-11-14）；第一次的实现留在本地分支 `p3.2-blocks-experiment`（`47a8d38`，未推送，已过时） | **第一次试做（暂不采用）**：**先量**：临时探针在每个实际迭代的分量（P3.1 之后）上做最大匹配、Dulmage–Mendelsohn 粗分解和方阵部分的强连通块。串联 3D 桁架（连 1、2、3 个前点）：全是方阵，771 个未知量分成 257 个 3 × 3 块；加了第 5 个前点的“加撑”桁架（基准夹具也是这种）结构上整体过定，分不开。2D 语料（P3.1 之后）：99.1% 的列在方阵部分，最大块只有 2 个未知量。装配语料：结构上整体过定（方程数超过结构秩），分不开。3D 语料小。**做法**（分支上）：`src/structure.rs`（Hopcroft–Karp 匹配，迭代实现；过定部分、方阵部分的 Tarjan 块、欠定部分；6 个单元测试，含 20 万行的链不爆栈）；`solve_one` 在分量分成两块以上时按序逐块求解（过定部分先、欠定部分最后取最小范数），每块只解自己的列、其余固定，`Assembly` 的布局加了“只排某些行 / 只放某些列”（`Shape`）；任何一块不收敛就放回原样走整体求解（P3.1 路径）。布局重构在开关关闭时逐位不变（20 865 次求解排序后逐字节相同）。**结果**：验收题达到——扰动 0.2 的串联桁架 20–260 点全部回到画好的构形（位移约 1e-12），每个最多 2 ms（原来 20–160 点“收敛”到 3–72 个单位外的另一个解，200、260 点 100 次迭代不收敛、要 12–27 s）；扰动 0.05 的 260 点从 2 917 ms 到 2 ms。2D 语料红线全过，翻面 5 → 3、最差硬漂移 9.0 → 4.5，软漂移和欠约束结果不变；但两个画得不对的圆链（k7、k8，圆心间距和半径不符，要大幅移动）落到另一个解：总位移更小（4.2 对 6.0 / 9.6），但有的圆跑到邻圆另一侧、互相重叠。装配语料逐项相同。**代价**：P3.1 之后 2D 的块只有 1–2 个未知量，每块的布局和牛顿开销占了大头，2D 矩形链比 P3.1 慢 1.5–3.6 倍（792 参数 0.79 → 2.9 ms，仍比 P3.1 之前快约 70 倍），矩形阵列 +50%；基准里的加撑桁架分不开，持平。**试过的变体**：每块开始前按起点处的线性化和前面各块的实际位移预测本块起点（相当于整体牛顿步在这些列上的一阶传递，只在降低本块残差时采用）：2D 语料与整体求解逐例相同（差 ≤ 1.1e-17），但桁架的翻面又回来了（位移 6–90），线性传递沿链放大——没有保留。**用户的决定**：暂不采用；将来采用时 `iterations` 报各块之和、`max_iterations` 限每一块。**以后再做时要解决的**：每块的固定开销（布局、分配、牛顿初始化）要先降下来，2D 才不吃亏；块内起点的选择（不预测时圆链翻到另一侧，预测时桁架翻面）需要一个更好的规则，可能就是 P3.3 的“选离当前构形最近的根”；过定部分（装配、加撑桁架）要靠更细的分解或稀疏分解（P2.2 的 `faer`）。**第二次（采用）**：用户选“重做 P3.2 + P3.3”。在 P3.4 的计划上做：分量计划建好时就分好块，每块的布局（只排本块的行、只放本块的列，`Shape`）和步空间（别名类按“类根 → 成员”的有序表直接取，`ClassSpace::of_roots`，不再每块扫一遍分量的全部参数——第一版这一步在 792 参数链上是平方级，冷求解慢 3.3 倍）都存在计划里；不超过 3 个未知量的方块用稠密牛顿（部分选主元消元，步长减半直到残差下降，收敛后像整体求解一样再走最多 `polish_steps` 步严格下降的精修）从当前位置直接解，解不了（奇异、无下降、迭代用完、精度够不着）就交给一般的块牛顿，再不行整个分量放回原样走整体求解。没有做起点预测（它会让桁架翻面）。P3.3 的几何闭式解（两圆交、线圆交）没有做：小块用的是稠密牛顿，“离当前构形最近的根”由从当前位置出发的牛顿给出——圆链 k7、k8 落到邻圆另一侧，正是这条规则的结果（例如 k7 的第 6 个圆，两个根离画的位置 2.5 和 5）。**结果**：验收题达到——扰动 0.2 的串联桁架 20–260 点全部回到画好的构形，每个最多 2 ms（新红线 `serial_trusses_come_back_from_a_large_perturbation`，在 `f5bdd23` 上 20 点就偏 2.98，失败）。2D 语料红线全过，翻面 5 → 3、最差硬漂移 9.0 → 4.5，软漂移与欠约束结果不变，两个圆链如上；装配语料逐项相同。新增 `tests/blocks.rs` 3 个（一块解不出来时 3D / 2D 分量放回原样、链上每点回到自己一侧；都是覆盖性的，3D 那个在去掉“放回原样”后失败），`solve_dense` 两个单元测试。**A/B**（与 `f5bdd23` 交替三轮，中位数）：重解 54 / 198 / 792 参数 −44% / −46% / −67%（792：0.472 → 0.158 ms）；冷求解 54 / 198 参数 +25%、792 参数 −9%，矩形阵列 20 / 200 个 +18% / +22%（第一次求解要分块并给每块排布局，在这些小分量上比整体求解还贵）；3D、装配、诊断 −2.6% 到 +4.3%，最小的两个装配夹具 +7.5% / +8.8%（以前在不相干的改动上也这样动过，没查）。用户拍板时看的是快速测量（冷求解 +10–17%、阵列 +15%），正式三轮更大，已如实写进 CHANGELOG 与 `benches/BASELINE.md`。**P3 的验收目标“500 实体草图改一个尺寸后重解 p50 < 1 ms”按它本身量了**（临时驱动，没提交：100 个矩形链，500 个实体、1 800 个参数，解好后反复改一个宽度再重解 2 000 次）：`f5bdd23`（P3.1 + P3.4）p50 248 µs、p90 318 µs；本步之后 p50 170 µs、p90 213 µs。**没做的**：几何闭式解；更好的块起点规则；冷求解的分块开销；过定部分（装配、加撑桁架）的细分解或稀疏分解 |
| P3.4 结构缓存 | 2D 完成；3D 草图与装配未做 | **先量**：临时计时器（没提交）在 2D 重解（偏离 1e-4，release，含每次克隆草图）上分阶段计时，µs/次：3 个矩形链（54 参数）共 29，其中分量划分 5、`solve_one` 里牛顿之外的准备（列出各方程的列、别名类、放位、布局）15、牛顿 7.7；44 个矩形链（792）共 507：67 / 222 / 204；200 个独立矩形（3600）共 3 228：498 / 1 822 / 684；100 个矩形链（1 800 参数、500 个实体）共 3 581：282 / 840 / 2 394。即这些夹具上结构性的准备占重解的 57–72%（100 个矩形链是 32%），按拓扑版本号缓存最多能省这么多。（当时据此推断“500 实体改一个尺寸后重解 p50 < 1 ms”的验收目标光靠缓存达不到——这个推断不对：上面量的是所有参数都偏 1e-4 的重解，牛顿要在 200 列上迭代整个系统；验收目标说的是改一个尺寸，后来按它本身量过，见 P3.2 行末。）核对过：所有结构性编辑（加删实体与约束、启用 / 抑制、替换约束、`fix_param` / `fix_entity`）都推进拓扑版本号，改尺寸值和移动几何不推进；但实体表是公开字段，宿主可以绕过版本号直接改。**做法**：`system::Plans`（分量列表 + 每个约束属于哪个分量 + 每个分量一个 `OnceLock<Plan>`，分量第一次求解时填入，并行求解时各线程填各自的槽）；`Plan` 保存分量碰到的列、整体布局（没有别名类的分量建计划时就排好，有的只在回退或 debug 核对时现排，所以计划里不再存一份列清单）、别名类的布局和步空间（`ClassPlan`，连同哪些行被吸收）。每次求解照旧重算随数值变的部分：别名类的偏移和已知值（`classes_of` 照跑，结构与计划一致，debug 下断言）、放位、分量尺寸、牛顿迭代。`Sketch` 新增不序列化的字段 `plans`（`Arc`，克隆共享），键是（拓扑版本号，实体表与各列表长度的 FNV 摘要）；`solve_selected_with` 改走 `system::solve_planned`。3D 草图和装配仍每次现做（`solve_components` 不带计划，计划用完即拆，不拷贝）。**验证**：逐位指纹——全测试集 20 865 次求解与 `491d1cc` 排序后逐字节相同（新测试除外）。新增 `tests/plan_cache.rs` 4 个：改尺寸后重解、四种结构性编辑后重解、直接交换实体表里的两个点、克隆后各自编辑，都和“同样状态、从没解过”的草图逐位比较；变异检查：键里去掉版本号时结构性编辑那个失败，去掉实体摘要时实体那个失败（另两个是覆盖性的）。七步门禁全过，语料指标不变。**弯路**：第一版 A/B 里冷求解的 200 个矩形阵列慢了 33%：布局的列映射 `col_of` 按“0 到分量最大列号”分配，200 个 18 列的小分量每个几千项，以前用完即还、分配器复用热内存，现在全保留着——改成只覆盖分量自己的列号区间（`col_base`），冷求解回到持平。同一轮里装配星形慢 4–5%：不带计划的族每次求解时多拷了列清单和模式，改成计划不保留时直接拆开用。修好这两处后的一轮里冷求解矩形链仍慢 +3.5–6.8%，再去掉计划里保存的那份列清单（建计划时要拷贝一次），回到 +0.2–2.1%。**A/B**（与 `491d1cc` 交替三轮，中位数）：重解 54 / 198 / 792 参数 −49% / −50% / −40%（0.0403 → 0.0206 ms、0.148 → 0.073 ms、0.768 → 0.461 ms）；冷求解 −4% 到 +2.4%；装配求解 −3% 到 +2.5%；诊断与 3D 是在去掉列清单之前那一轮量的，−1.6% 到 +3.6%，例外是 `assembly/diagnose_star/5`、`/40` +8.4% / +4.6%（三轮一致，没有重测最终版，原因没查）。**没做的**：3D 草图与装配的缓存；计划随数值变的部分（偏移、尺寸）仍每次重算；验收目标“500 实体改一个尺寸后重解 p50 < 1 ms”在这一步就已达到（见 P3.2 行末的测量） |

## 12. P4 进度

| 步骤 | 状态 | 说明 |
|---|---|---|
| P4.1 逐实体 / 逐约束状态 | 完成（2D、3D 草图、装配） | **用户定（2026-11-14）**：每个实体一个状态、按严重程度取最严重的（冲突 > 未求解 > 过定义 > 欠定义 > 定义完全），新方法 `status()`、`Diagnosis` 不动。**先量**：现有的是 `diagnose()`（dof、秩、冗余 / 冲突 / 部分冗余、依赖组）和逐参数的 `mobility()`，没有逐实体的东西；语料只有整例的标签。**做法**（`src/status.rs`）：约束状态 = 停用 > 冲突 > 未满足 > 部分冗余 > 冗余 > 满足；“未满足”按求解判收敛的同一把尺子（分量尺寸的单位、数值精度的权重），并且只在分量整体不在解上时才算——解好的分量不会因为舍入出现“未满足”。**一个依赖组里的约束取同一种状态**：第一版只用 `diagnose()` 列出的代表，冲突例子里同一个冲突的另一端被报成“未求解”（代表取决于主元顺序）；改成按依赖组后，组里有冲突就全组冲突，否则全组冗余（其中任何一个单独去掉都不改变解）。实体状态取它身上约束的最严重状态，再看它的自由参数里有没有可动的（欠定义）。已删除的实体是 `None`。**顺带发现并单独修了 `mobility()` 的错误**（提交 `7c983ee`）：消元选了主元列却不换行，主元位置上的那一行可能根本不含这一列；语料里 10 个开放距离链全都把被 `Fixed` 固定的起点报成 y 可动。**验证**：语料每个实体按“它所在的例子是怎么搭的”写出期望状态（不是抄求解器的输出）——满足族全部定义完全；欠约束族按缺了哪个约束逐个写（边距例子的孔缺半径 / x / y 时欠定义、钻孔永远欠定义；相切例子的线欠定义；距离链的起点定义完全、其余欠定义）；冲突族全部冲突（按构造每个实体都在冲突里）；冗余族全部过定义；精确相切的奇异族全部过定义（一阶依赖，P4.3 的见证构形要区分的正是这种）。13 075 个实体全部一致（`entity_statuses_match_the_labels`）；去掉按依赖组那一步，13 个实体不一致（变异检查）。另加 `tests/status.rs` 5 个（未求解、删除与停用、冲突各方、冗余只波及自己的几何、状态的次序）。**3D 草图（第二个提交）**：3D 没有 `mobility()`，而且平面存 4 个数只有 3 个自由度、圆柱存 7 个只有 5 个，逐参数的零空间判断会把法向长度当成“可动”。改成逐实体判断（`status::movable_entities`）：在规范的有效列上取零空间基（平面的法向长度不在其中），映回存储参数，限制到实体自己的参数上，再去掉这个实体自己的软规范方向（圆柱轴点沿轴滑动、轴向量缩放），还剩分量才算欠定义；没有规范时它等价于逐参数的可动性，2D 也改用它，13 075 个实体照样全对。零空间基由原来只在测试里用的稠密 RREF 改出（`numerics::null_space_basis`）。验证：被法向和一个点钉住的平面、经过 5 个固定点的圆柱都报定义完全（去掉有效列 / 去掉软规范投影后分别失败，变异检查），少一个点的圆柱、自由平面报欠定义；3D 语料的 8 个冲突例子全部实体冲突、8 个冗余例子解后全部过定义（未解时是“未求解”，按优先级它排在过定义之前）、自由度例子里自由的点 / 线 / 球 / 平面欠定义、固定点定义完全。**装配（第三个提交）**：`rigid3d::assembly_status`，实体是刚体、约束按槽位。**实测推翻了第一版的冗余口径**：按依赖组整组标冗余后，装配语料里几乎所有零件都成了“过定义”——刚体的配合方程之间天然相关（两对重合点共享两点连线方向上的距离，三对重合点 9 个方程秩 6），但这些配合哪一个都去不掉。冗余按文档本来的意思收紧为“单独去掉它，所在分量的秩不降”，逐个候选（一致依赖组的成员、诊断点名的行）去掉后重算秩；2D、3D 的标注不受影响（仍全对）。第一版这样做在装配语料上要 177 s（release）：挂在接地件上的 40 个零件是一个分量，每个候选都要对整个分量重做一次诊断。改法：对自己的列满秩的约束（固定点、带姿态的接地）单独决定了那些列，把它们按住后其余约束的秩差不变（分块：秩 = 被钉住的列数 + 其余部分的秩），于是其余候选在按住后分出的小块上测；锚点自己若是候选，仍在整个分量上测；只在没有规范的族上这样做（有规范时“对存储参数满秩”不是同一个判据）。177 s → 2.7 s。验证：装配语料 4 749 个刚体全部与“按构造写的期望”一致（`body_statuses_match_the_labels`：钉住、链、星形的零件定义完全，轴孔的轴欠定义、地定义完全，两个自由体欠定义，冲突族全冲突，冗余族全过定义）；把“可去掉”判据换成“候选都算”后只有 59 个一致（变异检查）。**没做的**：状态的代价（一次诊断 + 一次可动性分析，各按分量分解一次）没有单独量 |
| P4.2 自由度的几何描述 | 完成（2D、3D 草图、装配） | **先量**（临时探针，没提交）：对每个实体算两种局部运动——“跟随”= 全局零空间限制到这个实体上（其他几何可以跟着动），“独自”= 其他列全按住时这个实体自己的零空间——再按“平移 / 刚体转动 / 形状量”数维数。2D 语料里两者处处相同，出现的形态只有三种：开放距离链的点 65 个各 1 维（垂直于链：链是直的，第一阶上每个点都只能横着动）、边距例子的圆 21 个（圆心 1–2 维、半径）、`under-tangent` 的线 3 条 1 维刚体转动；另有恰好相切的 10 条线（奇异族，状态是过定义）第一阶上绕起点转。**`under-tangent` 的转动中心不可信**：相切点只能定到约 √容差（1e-5），PointOnCircle 与 Tangent 两行在解上几乎相同，秩 3 是近奇异给的，零空间那一维（探针给出“绕 2·起点转”）由残差的方向决定；按构造，有限运动是绕圆心滚动——标注里不查这个中心。装配里两者不同：轴孔例子的轴两者都是“绕孔的轴转”（8/8）；两个自由体用部分配合连起来（`under-partial`）时“跟随”永远是 6 维（整体漂浮），“独自”才是关节：重合 → 绕该点 3 个转动，距离 → 2 平移 3 转动，平行 → 3 平移 1 转动，夹角 → 3 平移 2 转动，同心 → 1 平移 1 转动，重合 + 平行 → 绕 y 的铰链。2D 弯曲的距离链中间点“独自”是 0 维，但它确实欠定义。**用户定（2026-11-14）**：两种都给；新方法 `freedoms()`（`status()` 与它的代价不变）；类型 `EntityFreedom { motions, own }`、`Motion2` / `Motion3`。我把预览里的形状量 StartAngle / EndAngle 改成了 `Sweep`（终止角 − 起始角）：单独的起止角在整体转动时也会变，不是形状量。**做法**：`status::motion_spaces` 是原来 `movable_entities` 的推广（同一个零空间、同样按全局最大分量缩放、同样去掉实体自己的软规范方向、同样的秩容差），每个实体把超过容差的局部向量正交化成一组基（够实体的参数个数就停），`own` 另在“只取完全属于这个实体的有效列”的矩阵上取零空间；`movable_entities` 现在就是“`follow` 非空”，所以 `motions` 为空当且仅当状态不是欠定义（一处细节：软规范方向的投影以前只减在运动里出现过的参数上，现在按整个方向减，是同一个定义的正确实现；语料上状态标注照样全对）。`src/freedom.rs` 描述一组基：每类实体给一个“框架”——沿世界轴平移、绕参考点（线的中点、圆心、平面离原点最近的点、刚体原点）转动时存储参数怎么变，形状量的变化率（线长、半径、2D 圆弧张角），以及要先去掉的存储冗余方向（平面法向与偏移一起缩放；圆柱轴点沿轴滑动、轴向量缩放）。然后按固定次序取基：先求运动空间与平移的交（主角的正弦 < 1e-6，Jacobi 特征分解），平移方向全部可动时就是世界轴；再求与刚体运动的交、去掉平移部分，解出每个转动的旋量（参考点速度、角速度；看不见的分量如线绕自身的转动由最小范数解取零），先试能否找到一个共用的点让每个转动都是纯转动（容许加上可用的平移），能就报“绕过这个点的轴转”（轴取世界轴在转动空间里的投影），不能就逐个报转动或螺旋；剩下的部分逐个形状量看能不能单独改变它，能的报 `Resize`，剩下的报 `Other` 并附单位速度。**验证**：语料按“例子是怎么搭的”写期望（`entity_freedoms_match_the_labels`）：直链的点横着平移、边距例子的孔缺哪个就是哪个（半径 / x / y）、钻孔沿 y 平移并改半径（约束去掉时再加沿 x）、恰好相切的线绕起点 (r, 0) 转、`under-tangent` 的线一个转动（中心不查）、冲突族里半径对直径的圆和冗余族里半径加直径的圆圆心两向可平移、矛盾边距的孔沿 y 平移并改半径，其余实体没有运动；`motions` 与 `own` 都对，且与状态一致——13 075 / 13 075。装配（`body_freedoms_match_the_labels`）：轴绕孔的轴（过孔的轴点、沿 z）转；`under-partial` 七种的 `motions` 都是 3 平移 + 3 个过自身原点的转动，`own` 是上面那些关节（球铰过重合点、距离的转动轴过连线、同心的平移沿轴且转轴就是那条轴线、铰链过重合点沿 y）；冗余族的夹角例子（位置固定的体 0 只能转、体 1 自由）；其余刚体没有运动——4 749 / 4 749，其中冲突族的 40 个按原样读、只核对与状态一致。另加 `tests/freedom.rs` 11 个：自由的点、线（2 平移、绕中点转、改长度）、圆弧（再加半径、张角）；一端固定的线绕它转并改长度，再加水平只剩改长度；圆心在固定线上的圆沿线平移并改半径；两端固定的圆弧只有一个 `Other`（半径与张角反向一起变）；弯曲链中间点有 `motions` 无 `own`、末点 `own` 只剩绕中间点；等半径的两个圆一起改半径而各自 `own` 为空；删除的实体为 `None`；3D 的自由点、线（两个过中点、垂直于线的转动）、球、平面（沿法向平移 + 两个转动，法向长度不算）、圆柱（垂直于轴的 2 平移 2 转动 + 半径）、样条（6 个刚体运动 + 3 个 `Other`），偏移固定的平面仍沿法向平移，固定平面上的点在平面内平移。新功能在旧代码上测不了，做了变异检查（改动后 `cmp` 确认生效、测完恢复）：`own` 换成 `motions` → 弯曲链测试和装配标注（16 个刚体）失败；去掉共用点 → 装配标注 1 个刚体（重合配合的体 1）失败（只有这一处能测出它）；2D 转动生成元反号 → 一端固定的线和 2D 标注（10 个实体）失败。**耗时**（快速测量：临时驱动，release，每项 5 轮取中位数，与 `83e1690` 比；矩形链解好后反复调用）：`status()` 在 3 / 11 个带尺寸的矩形链、11 / 44 个不带尺寸的矩形链上 0.30 → 0.26 ms、8.0 → 7.2 ms、6.3 → 6.0 ms、1.24 → 0.97 s（没变慢；变快的原因没细查，可能是去掉了每个零空间向量建一张按实体分组的哈希表）；`freedoms()` 0.19 ms、3.6 ms、3.1 ms、0.26 s。两者都是整个分量上的稠密零空间，大的欠约束分量很贵（44 个矩形链、792 个参数要近 1 s），只适合求解后调用。装配语料的 debug 门禁测试 71 → 84–87 s（多了新的标注测试）。**没做的**：稀疏的零空间 / 状态；3D 草图没有带标注的自由度语料（只有 `tests/freedom.rs`）；装配冲突族只核对了与状态一致；奇异构形上给的是一阶运动（P4.3 的见证构形要处理的） |
| P4.3 见证构形 | 完成（2D）；3D 草图、装配未做 | **先量**（临时探针，没提交，2D 语料 320 例，扰动幅度 0.1 和 0.5 结果相同）：见证 = 克隆草图、停用尺寸类约束（有 `value()` 的）和 `Fixed`、扰动自由参数、只解剩下的几何关系、再把尺寸类约束和 `Fixed` 的目标改成那里量出的值，然后诊断。奇异族（恰好相切）10/10 dof 1 → 0、冗余消失；冗余族 9/9 仍冗余；可解 256 例、欠约束 25 例全部不变（包括 `under-tangent`：它的近奇异来自几何关系本身——起点在圆上且线与圆相切，起点必是切点——见证上同样存在）；冲突族在见证上变成“冗余”（依赖是结构性的、只是数值矛盾）；平行 + 垂直的 3 例见证退化：`Fixed` 被放开后，“既平行又垂直”只能靠把线缩成零长度满足。**用户定（2026-11-14）**：新增 `ConstraintStatus::Singular`（在 `Met` 与 `Redundant` 之间）和 `EntityStatus::Singular`（在 `UnderDefined` 与 `OverDefined` 之间，压过奇异点上一阶的“可动”）；只做 2D，只在当前点有可去掉的约束时才构造见证。**做法**：`Sketch::witness`（私有）按上面的步骤构造见证，扰动幅度是图的尺寸（尺寸类约束的长度值的平均；没有就用参数的范围）的 0.1，种子固定；几何关系解不出来、或要靠零长度的线 / 零半径的圆才能满足时不用见证。`status_of` 新增钩子：当前点“可去掉”的约束交给见证，在见证上仍可去掉的照旧是冗余，不可去掉的是 `Singular`；见证上的总秩等于当前点的秩时，每个候选在见证上也一定可去掉（一般位置上任何约束子集的秩都不低于当前点），直接返回，不再逐个测。`Diagnosis` 不变（`dof` 仍是当前点的）。**验证**：语料标注中奇异族的实体由 `OverDefined` 改为 `Singular`，13 075 / 13 075；新的 `only_the_exact_tangencies_are_singular`：恰好有这 10 例出现 `Singular` 约束、它们没有 `Redundant`，冗余族仍有 `Redundant`。`tests/status.rs` 加一对对照：圆上的点再加“到圆心 x 距离 = 半径”是 `Singular`（只有在竖直切点处两行才都是径向），换成“到圆心距离 = 半径”是 `Redundant`（处处重复）。`sketch.rs` 两个单元测试：见证离开奇异点 > 1e-3、仍在圆上、尺寸值等于量出的值、所有约束成立（重解 0 次迭代）、见证上 dof 0 而当前点 1；“既平行又垂直”的草图没有见证。变异检查（改动后 `cmp` 确认、测完恢复）：去掉见证 → 语料两个测试失败；见证时不停用尺寸类约束 → 只有单元测试失败（语料与对照测试仍过：见证被解回奇异点，但二重根上牛顿只收敛到约 √容差 = 1e-5 之外，按 1e-9 的秩容差已是独立的——所以必须直接检查见证离得够远）；不改目标值 → 只有单元测试失败（秩不取决于目标值）；不做退化检查 → 只有单元测试失败（退化见证只会让 `Singular` 变少，语料里没有这种组合）。**耗时**（快速测量，临时驱动，release，5 轮中位数，噪声大）：3 / 11 / 44 个矩形链各加一条重复的 `Horizontal`，`status()` 无见证 0.93 / 20 / 4 484 ms，逐个测时 1.2 / 49 / 7 711 ms，加了“秩相等就直接返回”后 1.4 / 32 / 4 617 ms（+50% / +57% / +3%；小夹具上主要是克隆和那次几何关系求解）。没有可去掉的约束时没有额外代价。**没做的**：3D 草图与装配（语料里没有奇异类可以验收；装配几乎每个分量都有相关行，每次都要构造见证）；被固定的几何（`fix_entity` / `fix_param`）不参与扰动，所以“固定几何的位置 + 几何关系”造成的奇异分不出来；见证只取一个固定的伪随机构形；退化检查只看线长和半径 |
| P4.4 冲突定位 | 完成（2D、3D 草图、装配） | **先量**（临时探针，没提交）：2D 语料 20 个冲突例子里，诊断的冲突依赖组已经就是删除过滤得到的最小冲突集（从组开始和从全部约束开始结果相同），并且和按构造推的一致；种子故障里，“重复固定加第三个不同的固定”的依赖组已经分成 {0,1}、{0,2}，“冲突加无关的重复约束”也分得开——**线性冲突上依赖组本身就够**；但**非线性冲突（两个相距 5 的固定点、第三点到两者距离都是 1）诊断完全看不到**：在画好的位置上和求解失败后（失败的求解会放回起点）两行都独立，没有冲突、没有依赖组，`status()` 只报未满足；在整个失败分量上做删除过滤能得到 {两个固定、两段距离}。**用户定（2026-11-14）**：`conflicts()` 返回 `Vec<Conflict>`（`constraints` 按添加顺序最新在前，`suggested()` 是第一个）；三族都做；`status()` 不变。**做法**（`src/conflict.rs`）：对每个连通分量，先在副本上只开这个分量的约束求解一次，解得出就跳过；解不出就做删除过滤——从最新的约束开始逐个停用，停用后仍解不出的不再放回，剩下的就是一个最小冲突集（再去掉任何一个都能解）；诊断里有冲突依赖组且它单独就解不出时从组开始，省求解次数；记下后停用其中最新的一个，再从剩下的约束开始，直到分量能解，所以把每个冲突的建议项都去掉，整个系统就能解。“解不出”= 从当前几何求解不收敛，**或者只靠退化几何才收敛**（2D：零长度的线、零半径的圆 / 圆弧；3D：零长度的线、零长度的平面法向或圆柱轴、零半径；尺度取图的尺寸的 1e-6，2D 与 P4.3 的见证共用 `drawing_size` / `has_degenerate_geometry`）。加这一条的原因是实测：3D 的两个矛盾平面法向被求解器“解”成了零法向量（1 次迭代、残差 0、报 `Converged`），2D 的平行 + 垂直在放开起点的 `Fixed` 后靠把线缩成零长度满足。**验证**：2D 语料 `conflicts_match_the_labels`：20 个冲突例子都报出按构造写的唯一一个最小冲突集（半径对直径、两个固定：这两个；超出跨度的距离：两个固定和距离；矛盾边距：线起点与锚点的重合、线长、两个 x 距离——锚点的固定和线的水平不在内；平行 + 垂直：这两个），去掉建议项后能解，草图参数不动；其余 300 例（从解好的状态）没有冲突。装配语料同名测试：23 个冲突例子，两次固定 {1,0}，重复距离 / 重复角度 / 角度对垂直 / 重合对距离都是两条配合 {2,1}（体 0 的固定不在内：两条配合都是相对的），**两个同心 {1}**：按构造应是 {1,0}，但两个刚体都是单位姿态、第二个同心要对齐的两根轴恰好垂直，这里它的残差没有一阶斜率，求解 0 次迭代就 `Stalled`（第二个会话定的“退化起点保持 Stalled”），按“从这里解不出”的定义它单独就是冲突——标注按这个写并注明原因，ROADMAP 记了“从扰动起点再试”；其余例子从解好的姿态没有冲突。3D 语料 8 个冲突测试各加断言：两个固定 / 两个法向 / 两个球半径 / 两个圆柱半径 {1,0}，超出跨度、线长对两端 {2,1,0}，**平行 + 垂直 {0}**（两条线都被固定且本来垂直，单独的平行就和固定的几何矛盾——我第一版标注写成 {1,0}，错在我），点被固定在平面外 {0}。`tests/conflicts.rs` 7 个种子故障：够不着的两段距离（诊断无冲突无依赖组，`conflicts()` 报 [3,2,1,0]，草图不动）；同一分量再挂一个无关的点（被过滤掉）；重复固定加第三个（[2,0]：保留需要的最旧的那个）；冲突旁边一对一致的重复竖直（[7,2]，重复不在内）；一个分量里两个冲突（{1,0} 和 {4,3}）；两个分量各一个；被抑制的约束不参与。`src/conflict.rs` 4 个单元测试（过滤次序、逐轮去掉建议项、从依赖组开始的求解次数、能解的分量）。变异检查：过滤不删任何约束 → 单元测试、种子故障（无关的点那个）、2D 语料都失败；过滤改成最旧优先 → 只有单元测试失败（语料和种子故障里依赖组提示已把起点缩到最小集）；2D 去掉退化判定 → 2D 语料失败（平行 + 垂直多出那个固定）；3D 去掉退化判定 → 3D 语料失败（两个法向报不出冲突）。七步门禁全过。**没做的**：从扰动起点重试（区分“解不到”与“不可能”）；求解次数（每个冲突约为集合大小次分量求解）在大草图上的代价没量，大分量可以换 QuickXplain；3D 求解本身接受零法向的退化解（`solve()` 报 `Converged`），这是 3D 族的已有问题，修它会改变求解结果，没动，记在 ROADMAP |
| P4.5 what-if `would_add` | 完成（2D、3D 草图、装配） | **用户定（2026-11-14）**：显式探针对象（`Sketch::what_if()` → `WhatIf`、`SpaceSketch::what_if()` → `SpaceWhatIf`、`rigid3d::assembly_what_if(system, positions, rotations)` → `AssemblyWhatIf`），`would_add(&c) -> Addition { verdict, removes }`；口径对照完整诊断（把候选加到副本上诊断：秩增加等于方程数 → `Independent`；否则新增了冲突项 → `Conflicting`，新增了部分冗余项 → `PartiallyRedundant`，否则 `Redundant`；`removes` = `dof` 的减少）；三族都做。**先读**：诊断按行均衡、带列主元的 Householder QR 定秩（剩余列范数 ≤ 1e-9 × 最大初始行范数），相关行用独立行的（λ = 1e-10 正则的）最小二乘步算偏差 e，e 与 `CONSIST_TOL` / `PARTIAL_TOL` 乘尺度比较；分量按 Jacobian 行列出的参数（含固定的）连通。**做法**：先把 `diagnose()` 里逐分量的分析抽成 `system::diagnose_component`（纯重构；2D 语料每例画好时和解好后的诊断与状态、装配语料每例的诊断与状态共 973 行，重构前后逐字节相同），它可以顺带保留 `RowBasis`：有效列（全局参数的组合）、独立行的正交基（MGS 两遍）、最大行范数、最小相对主元（条件指标）、最小二乘步、每个相关行的偏差、最大残差、原始行（算软规范用）。`src/whatif.rs` 的 `RowSpaces`：每个分量第一次被问到时建（`OnceCell`）；问一个候选时，在“追加了候选的方程视图”（2D / 3D 复制一份约束表，装配复制系统）上取候选的行与残差，映到合并后分量的有效列、行均衡，逐行投影到各分量的基与本次新增的方向上：剩余超过秩容差就是独立的（扩充方向，并在正交补里精确修正最小二乘步），否则是相关的，偏差按诊断的规则分类；尺度取合并分量的单位与所有残差的最大值；软规范方向（圆柱轴点滑动、轴缩放）被钉住的个数从 `removes` 里扣掉。**一开始的实现与诊断不一致**：2D 语料 10 096 个候选里，相切族（解在近奇异处）和按原样读的平行 + 垂直冲突例子上有一批不一致——病态分量上，完整诊断的列主元会先选候选行（范数 1）作独立行、换出一条病态行当相关行，用条件好的子集判出冗余，而增量投影用病态的最小二乘步（被放大到约 1e-6）判成部分冗余 / 冲突；本身已不一致的分量上，“新冲突落在哪一行”同样取决于主元顺序，且相关关系的系数会把偏差按 1/|y| 缩放。于是快路径只在结论稳定时下判断：候选有相关行时，要求各分量最小相对主元 ≥ 1e-3、分量本身没有冲突或部分冗余、候选相关行的偏差离两个阈值都在 100 倍以外（部分冗余因此总由精确路径判），候选行的剩余离秩容差在 1000 倍以内也走精确路径；另外固定列变了、有效列对不上（候选碰到平面上分量没碰到的参数）、合并了没有自由列或秩为 0 的分量、候选不连任何分量时，都对涉及的分量做精确诊断（`status::Only` 只开这些约束，前后各诊断一次）。**验证**：2D 语料（矩形族只取 k ≤ 4）每例从前 4 个点、3 条线、2 个圆生成各类候选（固定、重合、距离 / x / y 距离、水平竖直、线长、平行垂直、夹角、半径、点在圆上、相切、等半径；尺寸取量出的值与偏开的值），10 096 / 10 096 与完整诊断一致（冲突 5 039、冗余 4 482、独立 535、部分冗余 40），约九成走快路径（部分冗余全部、冲突分量上的相关候选、相切族的部分走精确路径）；装配语料抽样（小例子隔一个、超过 4 个刚体的每 20 个取一个，77 例）前 3 个刚体之间各类配合（固定位置 / 偏开、带姿态的固定即接地、重合 / 偏开、距离、平行、垂直、夹角、同心）1 443 / 1 443 一致；`tests/whatif.rs`：3D 两组夹具（法向长度 ≠ 1 的平面、平面上的点、线、球、圆柱，以及一组自由实体）129 / 129 一致，连接两个分量 / 只碰固定几何 / 半冗余的固定（`removes` 1）各一例，“同心于固定球”钉住已在分量里的圆柱轴点：秩 +3、`removes` 2。变异检查：不投影 → 2D 语料与 `tests/whatif.rs` 失败；独立行不修正最小二乘步 → 三处都失败；去掉条件 / 一致性护栏 → 2D 语料失败（装配仍过：那里分量条件好且一致）；不扣软规范 → 第一版的软规范测试**没有失败**（自由圆柱不在任何分量里，走的是精确路径），改成圆柱已在分量里后失败。**耗时**（临时驱动，release，与“副本 + 加约束 + 诊断”比）：11 个矩形链（198 参数）带 / 不带尺寸：每问 128 / 107 µs 对 6.8 / 4.5 ms（53× / 42×）；44 个矩形链（792 参数）：3.3 / 3.7 ms 对 976 / 722 ms（299× / 193×），40 个候选全部走快路径；第一次问一个分量要 1.6–2 倍诊断的时间（多建正交基和保留数据）。语料小草图上整体只快约 2.6 倍（分量小，克隆约束表和建基占大头）。新增公开的 `rigid3d::Constraint::equation_count()`（测试的对照口径要用，2D / 3D 早有）。**更正（P4.7 之后的单独提交）**：第一版以为装配里带姿态的 `Fixed` 会让刚体接地、改变未知列，加了“接地集合变了就走精确路径”的判断，并在文档里写“接地 `Fixed` 按秩口径读作冗余”；实际上诊断用的方程视图（`RigidEquations::new`）不做接地，带姿态的 `Fixed` 就是 6 行普通方程（只有求解用的 `grounded()` 才排除接地体），那条判断永远不触发、那句文档是错的——加在自由刚体上是 `Independent`、`removes` 6（新测试 `a_full_fixed_on_a_free_body_takes_its_six_freedoms`）；死代码已删、文档已改。**没做的**：不含 P4.3 的见证与“单独可去掉”的测试（`would_add` 说冗余的，`status()` 里可能是 `Singular`）；投影是稠密的（秩 × 列），大分量上第一问是一次完整分解 |
| P4.6 自动标注 | 完成（2D、3D 草图） | **用户定（2026-11-14）**：返回建议、不改草图（`Sketch::auto_dimensions(datum: Option<PointRef>) -> Vec<Constraint>`、`SpaceSketch::auto_dimensions(datum: Option<SpacePointRef>) -> Vec<SpaceConstraint>`）；只用尺寸；做 2D 和 3D 草图。**做法**（`src/autodim.rs`）：在当前点取物理零空间——按分量在规范的有效列上取零空间基（与 `motion_spaces` 同一套），没被任何约束碰到的自由参数也按有效列给方向（自由平面的法向长度不算），再去掉没被钉住的软规范方向（圆柱轴点滑动、轴缩放），所以维数就是诊断的 `dof`——存成稀疏向量的正交基。候选按固定次序给出、值取量出的值；候选的行（用只含这一个候选的方程视图取，不复制约束表）投影到基上，**每一行**都有之前的行没有的分量（相对行范数 > 1e-6）才采用，采用后只改写与这些行相交的基向量（对它们张成的子空间取去掉这些方向后的正交基，带主元的 Gram–Schmidt）。次序：基准点（给定的；否则第一个动不了的点；否则第一个点），基准点能动时先给它一个 `Fixed`；然后按实体编号：圆 / 圆弧的半径、线长（3D：线长，球 / 圆柱 / 圆 / 圆弧的半径，平面的法向与基准点到平面的距离），再给它的各个点到基准点的 x、y（3D 再加 z，用 `PointDistanceAlongDir`）距离。预览里说的“剩下的转动用线对水平的角度”没有做：2D 的所有参数都能由点的位置和半径定下，用不到。基准点只有一个坐标被约束住时，它的 `Fixed` 只有一行独立、会被拒绝；我一度以为这会漏掉它的另一个坐标，加了补救分支，实测发现主循环会把固定点相对基准点的距离也提议一遍、那一行只含基准点的列，已经覆盖，于是撤掉了补救分支。**验证**：2D 语料 `auto_dimensions_define_every_case`：300 个非冲突例子（解好后）加上全部建议后 dof = 0，新加的约束在 `status()` 里全是 `Met`（没有冗余、在当前几何上成立），建议的方程数之和等于原来的 dof，重解收敛；欠约束族 25 例的建议按构造写出、逐条一致（边距例子：缺半径给孔半径、缺 x / y 给孔心相对锚点的 x / y，钻孔给半径与 y、约束去掉时再给 x；直链：每个点相对锚点的 x；`under-tangent`：固定圆的圆心作基准，线起点的 x），共 117 条建议（含奇异族各 1 条、冗余族里圆心自由的圆各 1 个 `Fixed`）。`tests/autodim.rs` 6 个：自由的点 / 线 / 圆 / 圆弧（先固定第一个点）、给定会动的基准点（先固定它、其余从它量起）、已完全定义的草图（空）、只有一个坐标被约束的基准点（不给 `Fixed`，固定点相对它的 y 距离定下它）、默认基准是第一个动不了的点、3D 自由点 / 线 / 球 / 平面全部定下（无冗余），自由圆柱只定半径、轴的 4 个自由度留下。`src/autodim.rs` 单元测试：零空间维数等于诊断的 dof（自由圆柱 5、自由平面 3、过固定点的圆柱 4）。变异检查：接受只有部分行独立的候选 → `tests/autodim.rs` 3 个与 2D 语料失败；不去掉软规范方向 → 只有单元测试失败（现有词汇里没有只钉软规范的候选，对建议没有影响）；不找“动不了的点”当基准 → 只有新加的默认基准测试失败（语料里第一个点恰好都是定死的）。**没做的**：关系类建议（水平 / 竖直 / 重合）、“人会这样标”的偏好（链式、基线、对称）；3D 圆柱轴与样条控制点；装配的自动配合 |
| P4.7 从动尺寸 | 完成（2D、3D 草图、装配） | **先看现状**：三族都已能量当前值——`Constraint::measure(params, entities)`（不必在草图里）、`Sketch::measure(index)`（任何状态，被抑制的尺寸就是从动尺寸），3D 与装配对应的方法——缺的是“量一个不在草图里的尺寸”的便捷接口和灵敏度。**用户定（2026-11-14）**：`sensitivities(&driven) -> Option<Vec<(usize, f64)>>` 加 `measure_constraint(&c)`；三族都做（装配 `rigid3d::assembly_sensitivities(system, positions, rotations, &driven)`）；冗余的驱动尺寸按最小范数分配。**做法**（`src/sensitivity.rs`）：从动尺寸先取它量出的值，在 `r(x, v) = 0` 处取它的 Jacobian 行与残差对值的导数（中心差分；夹角这类残差对值不是线性的），得 `∇m`；只取 `∇m` 碰到的分量，在规范的有效列上把行均衡；用带主元的 Gram–Schmidt 取行空间的正交基 Q，`∇m` 不在其张成内（投影残差 > 1e-6 相对）就返回 `None`；否则令 L = J Qᵀ，对 L 做薄 QR，解出满足 Jᵀ y = ∇m 的最小范数 y，`∂m/∂v_j = −y_j · ∂r_j/∂v_j`（每个单方程的尺寸类约束，导数同样用中心差分）。**第一版用正规方程**（`solve_normal_equations`）求 y：相切族的分量里有近奇异的行，J Jᵀ 的条件数约 1e10，分解失败时它返回零向量，于是一个显然确定的尺寸（被 x / y 距离钉住的起点到圆心的距离）被判成“不确定”——改成上面不经过正规方程的做法后一致。**装配的一处弯路**：我以为接地体的列不是未知量、要从 `∇m` 里去掉，加了过滤；变异检查（去掉过滤）测不出，查下来是诊断用的 `RigidEquations::new` 不做接地（带姿态的 `Fixed` 在诊断里就是 6 行方程，只有求解用的 `grounded()` 才排除接地体），过滤是死代码，已删掉——同一个误解也在 P4.5 的 what-if 里（见下一个提交）。**验证**：2D 语料 `sensitivities_match_finite_differences`（矩形族取 k ≤ 2）：每例的单值尺寸候选（量出的值）中 3 044 个被确定、213 个不被确定，**恰好与 `what_if().would_add()` 判相关 / 独立一致**；对每个驱动尺寸改 ±h（h = 1e-5 × max(1, |v|)）、重解、量所有从动尺寸，16 609 个导数与中心差分一致（最差相对误差 1.1e-6）；15 个驱动尺寸单独改了解不出来（冗余的），跳过；相切族不做差分（解在近奇异的相切处，改 h 解移动约 √h）。装配语料抽样：179 个确定、15 个不确定，同样与 `would_add` 一致。`tests/sensitivity.rs` 6 个：4 × 3 矩形对角线对宽 0.8、对高 0.6；去掉高后对角线 `None`、x 距离对宽 1；宽给两次各 0.5；两个固定点之间的距离依赖为空；3D 夹具（固定点、平面上的点、沿 y 的距离、固定的球心与半径）与中心差分一致；装配（接地体 + 球铰 + 两个距离 + 一个夹角，取自一个给定姿态的量值）与中心差分一致，去掉夹角后 `None`。变异检查：去掉确定性判定 → 2D 语料与 `tests/sensitivity.rs` 2 个失败；导数反号 → 都失败。**没做的**：保留分解的探针（多个从动尺寸时更快）；奇异构形上给的是线性化的导数 |

**P4 验收对照（2026-11-15）**：计划的三条验收——① 语料逐实体状态 100%：2D 13 075 / 13 075、装配 4 749 / 4 749（P4.1），奇异族改为 `Singular` 后仍全对（P4.3），逐实体运动描述 2D 13 075 / 13 075、装配 4 749 / 4 749（P4.2）；② 种子故障上最小冲突集与预期一致：2D 20 例、装配 23 例、3D 8 例的冲突族全部按构造一致，`tests/conflicts.rs` 7 个种子故障（含线性诊断看不到的非线性冲突）一致（P4.4；装配“两个同心”因退化起点按“从这里解不出”的定义报单个配合，已注明）；③ what-if 与完整诊断 100% 一致且快一个数量级以上：2D 10 096 / 10 096、装配 1 443 / 1 443、3D 129 / 129，792 参数的分量上每问快 193–299 倍（P4.5；第一问要建基，约 1.6–2 次诊断）。另外完成了自动标注（P4.6，2D 语料全部可解例子加上建议后完全定义、无冗余）和从动尺寸的灵敏度（P4.7，16 609 个导数与中心差分一致）。都达到。

## 13. P5 进度

### 13.1 起步基线（2026-11-15，先量后改）

现在**没有拖拽 API**。宿主今天只有两条路，都用现成的公开接口量了一遍（临时探针 `#[ignore] fn report_drag_baseline()`，release 构建，量完删除）：

- **warm-set**：把光标坐标写进被拖点的参数（`set_params`），再 `solve()`。被拖的点是自由的，求解器没有把它拉向光标的任何理由。
- **temp-Fixed**：每帧加一条 `Fixed { point: 手柄, x, y }`，`solve()`，再 `remove_constraint`。

语料是六段脚本化拖拽（都是闭合路径，便于看往返）：自由矩形（原点固定、四边水平 / 垂直、没有尺寸，2 自由度）的角点绕一圈、同一个角点穿过零宽度再回来、四连杆（两个固定铰、三段杆长，1 自由度）的耦合点在可行范围内摆动、同一个耦合点被拉到可达区域之外、五个相切单位圆链（4 自由度）的末圆心绕一圈、40 个自由矩形串成的一个分量（200 实体）的末角点绕一小圈。

| 例子 | 策略 | 帧 | 不收敛 | 翻面 | 回到起点 | 往返漂移 | 跟踪误差 p50 / max | 多余移动 max | 帧时 p50 / p99 (ms) |
|---|---|---|---|---|---|---|---|---|---|
| 矩形角 · 绕圈 | warm-set | 65 | 0 | 0 | 是 | 3.95e-1 | 0.418 / 0.419 | 0.140 | 0.004 / 0.005 |
| 矩形角 · 绕圈 | temp-Fixed | 65 | 0 | 0 | 是 | 4.44e-16 | 0 / 0 | 0.147 | 0.011 / 0.016 |
| 矩形角 · 穿过零宽 | warm-set | 97 | 0 | 2 | 是 | 4.37e-1 | 0.437 / 0.438 | 0.146 | 0.004 / 0.004 |
| 矩形角 · 穿过零宽 | temp-Fixed | 97 | 0 | 2 | 是 | 0 | 0 / 0 | 0.146 | 0.010 / 0.015 |
| 四连杆 · 可行弧 | warm-set | 97 | 0 | 0 | 是 | 1.10e-2 | 1.792 / 3.513 | 0.079 | 0.012 / 0.013 |
| 四连杆 · 可行弧 | temp-Fixed | 97 | **95** | 0 | 是 | 0 | 1.820 / 3.640 | 0 | 0.061 / 0.098 |
| 四连杆 · 不可达 | warm-set | 97 | 0 | **47** | **否** | 2.81e0 | 5.769 / 14.483 | 2.365 | 0.014 / 0.026 |
| 四连杆 · 不可达 | temp-Fixed | 97 | **95** | 0 | 是 | 0 | 7.616 / 15.232 | 0 | 0.161 / 0.216 |
| 相切圆链 · 绕圈 | warm-set | 65 | 0 | 2 | 否 | 7.86e-1 | 0.085 / 0.206 | 0.112 | 0.018 / 0.024 |
| 相切圆链 · 绕圈 | temp-Fixed | 65 | 0 | 4 | 否 | 3.34e-5 | 0 / 0 | 0.233 | 0.022 / 0.042 |
| 200 实体链 · 绕小圈 | warm-set | 33 | 0 | 0 | 是 | 4.01e-1 | 0.487 / 0.498 | 0.166 | 0.098 / 0.164 |
| 200 实体链 · 绕小圈 | temp-Fixed | 33 | 0 | 0 | 是 | 4.44e-16 | 0.000 / 0.000 | 0.195 | 0.252 / 0.351 |

量法说明：**跟踪误差** = 这一帧解完后手柄到光标的距离；**多余移动** = 这一帧里不属于光标的参数中变化最大的那一个（最小移动要压的就是它）；**翻面** = 看着的三角形取向（矩形的两条邻边、连杆的两个肘、圆链的相邻三心）出现 +1 ↔ −1 的真正反号，经过共线（符号 0）本身不算，到了另一侧才算；**回到起点** = 闭合路径走完后取向与起点一致；**往返漂移** = 闭合路径走完后参数与起点的最大差。

结论（P5 要解决的正是这些）：

1. **两条现成的路各坏在一处。** warm-set 永远收敛，但手柄根本不跟着光标走（矩形上固定差 0.42，四连杆上 1.8–3.5），而且闭合拖一圈不回到原构形（矩形 0.40、圆链 0.79、不可达 2.81）——这就是迟滞。temp-Fixed 跟踪是精确的（0）、往返按位回原处（4e-16），但**手柄处自由度少于光标两维时就是过约束**：1 自由度的四连杆上 `Fixed` 多加两条方程，97 帧里 95 帧不收敛；而且每帧改拓扑会废掉结构缓存，帧时是 warm-set 的 2–3 倍。所以 P5.1 的被拖目标必须是**软目标**（最小二乘拉过去），不能是硬约束，并且不能动拓扑。
2. **翻面出现在该出现的地方**：矩形穿过零宽度 2 次、相切圆链绕圈 2–4 次、拉到可达区域之外 47 次（构形被毁、取向没回来、漂移 2.81）。P5.3 的签名守卫与 P5.5 的不可达策略各对应一类。
3. **帧预算目前不是瓶颈**：200 实体一个分量，warm-set p50 0.098 ms / p99 0.164 ms，temp-Fixed 0.252 / 0.351，都在 P5.6 的 1 ms 以内。注意这个夹具没有尺寸约束（否则拖不动），上了尺寸行数会涨，真实拖拽路径定下来后要重量一次。
4. 签名需要对"共线"做迟滞：相切圆链绕圈时取向会经过 0，第一版探针把它当成翻面，多数了一倍（4 → 2）。


### 13.2 P5.1 拖拽会话（完成，2026-11-15）

| 步骤 | 状态 | 说明 |
|---|---|---|
| P5.1 拖拽会话 API | 完成（2D） | 见下 |

**用户定（2026-11-15，四个问题都选了推荐项）**：① API 是 `Sketch` 上的三个方法（`begin_drag` → `drag_to` → `end_drag`），会话状态存在草图里、句柄不借用（拖拽中仍可诊断）；② 一个会话可以带多个被拖点；③ 每帧返回新结构 `DragFrame`（`SolveResult` 的内容 + `tracking_error` + `reached`），`SolveResult` 不动；④ 软目标与硬约束**分层**：硬约束照常用牛顿满足，软目标只在零空间里走最小范数步。我另外加了 `drag_to_with` / `end_drag_with`（和 `solve_with` / `solve_selected_with` 一致的惯例，P5.6 的粗容差要用），汇报里说明了。

**做法**：每帧走"预测 — 校正"对，最多 8 个子步（`src/drag.rs`）。

- **预测**（`system::tangent_step`）：在当前构形处，取让各手柄最贴近其光标的最短运动，运动限制在约束留下的自由方向里。零空间投影器 `P = I − J⁺J` 用 `λ(JᵀJ + λI)⁻¹` 代替——沿奇异值 σ 的方向特征值是 `λ/(σ² + λ)`，自由方向是 1、被约束方向是 `λ/σ²`——这样**不用对雅可比另做一次分解求零空间基**：牛顿步自己的正规方程分解一次，每个被拖参数一次回代就给出投影器的一列（`λ N⁻¹ e_p`）。各手柄的行组成 `A = R P`，用 `min_norm_step_within(A, d, 1e-4)` 取最小范数最小二乘解，它已经落在零空间里，直接就是这一步。方程先按行均衡（最大元为 1），`λ = 1e-8 × JᵀJ 对角最大元`，所以 λ 对几何相对、与方程单位无关。
- **校正**：对会话涉及的分量 `solve_planned` 热启动重解（走现成的别名类与分块计划），把硬约束拉回零。不收敛的子步整步回滚（"没有发生过"）。
- **停下来的条件**：到达光标（`1e-9 × 尺度`）、子步不再把跟踪误差压低 0.1%（此时若变差就回滚，保留较好的那个构形）、或者没有可动方向。
- 没有任何约束碰到的被拖参数直接写成光标值（不用求解）；`fix_param` 按住的不动。
- `tracking_error` = 各手柄到其光标距离的最大值；`residual` 只算涉及分量的方程（拖拽不动其余部分）。过期句柄（再 `begin_drag` 之后）拖不动东西、报 `Stalled`。

**两个必须的护栏，都是实测出来的**：

1. **秩截断要用 `1e-4`，不能用数值秩 `STEP_RCOND = 1e-12`**。投影器给被约束方向留下的是 `λ` 量级而不是零，远在数值秩之上；一旦被当成秩，它就会用"除以这么小的数"去换光标在该方向上的那一份，步长放大上亿倍。实测：四连杆（1 自由度，光标在曲线外）第一帧把耦合点从离光标 15.2 推到 20.3；截断后变成 15.2 → 13.5 → 13.4 沿曲线走，第三个子步变差被回滚。
2. **逐个手柄的"能不能动"判定**（它自己那一项投影值 / 梯度模方 ≤ `1e-4` 就丢掉这一行）：所有手柄都被按住时，`A` 的所有行都是 `λ` 量级，秩截断是相对最大列范数的，照样会取秩。实测（变异检查）：去掉它，被两条距离三角定位的点被拖时会走一个 3.8e6 的步再被校正拉回来——位置上看不出来，但 `frame.iterations` 从 0 变成非 0。

**验证**（`tests/drag.rs`，6 段脚本化拖拽 + 7 个单例，都是新 API，旧代码上无法失败，所以做了变异检查）：

- 每帧硬约束都满足（最差帧残差 0）、可达的光标每帧都到达（最差跟踪误差 3.8e-10）；
- 闭合路径回到原构形：手柄能定住其余构形的四个例子漂移 ≤ 2.8e-10（矩形角 2.83e-10、四连杆 1.3e-11、200 实体链 3.75e-10）；
- 不可达的光标：`reached` 为 `false`、约束仍满足、手柄比开始时更近（走到走不动为止）；
- 翻面：可达的几例 0 次；仍然会翻的两例**断言它会翻**（矩形穿过零宽度——光标自己要求镜像构形；四连杆被拉出可达范围——走到连杆成一条直线再越过去），P5.3 的守卫落地时这两条断言会失败，正好提醒更新；
- 单例：没有约束的点精确跟随、被约束按住的点不动且 `iterations == 0`（用两条距离三角定位，不用与固定点重合——重合会被别名消元吞掉坏步，变异检查就测不出来了）、一帧跨越整张图落在逐帧走到的同一处（≤1e-6）、一帧把相切圆链折过去 4 个半径、两个手柄同时拖、过期句柄拖不动。
- **变异检查**：去掉秩截断 → 2 个测试失败；去掉逐手柄可动判定 → 1 个失败。**两条测不出**（如实记录）：去掉子步的推进上限（`aim` 不设上限）——语料里所有数字一字不差，因为光标每帧挪得小，大跳也被 8 个子步和"不再改善就停"兜住了；去掉"校正不收敛就回滚"——语料里没有一帧的校正失败。"不再改善就停"本身也没有测试失败，但它值 5.5 倍帧时（四连杆可达弧 0.013 → 0.072 ms p50，每帧都把 8 个子步用完）。

**和旧做法的对照**（`report_drag_frames`，数据进了 `benches/BASELINE.md`）：矩形角绕圈，新 API 跟踪 2.8e-10 / 往返 2.8e-10，"把光标写进参数"跟踪 0.42 / 往返 0.40；四连杆可达弧，新 API 每帧收敛、往返 1.3e-11，"临时 `Fixed`" 95/97 帧不收敛；四连杆被拉出范围，新 API 1 次翻面、残差 1.6e-12，"写进参数" 47 次翻面、往返漂移 2.81。相切圆链绕圈三种做法往返都是 0.95–0.98：那是 4 自由度里手柄只定住 2 维带来的路径依赖（物理上的链子拖一圈也不会自己捋直），不是求解器的迟滞，所以测试只对"手柄能定住其余构形"的例子要求往返，对其余要求"漂移不超过光标自己走过的距离"。

**帧时间**：矩形角 0.017 ms、四连杆 0.013 ms、圆链 0.038 ms（p50），都在 P5.6 的 1 ms 以内；**40 个自由矩形（200 实体、一个分量）每帧 19.8 ms，差两个数量级**——预测器对整个分量的正规方程做稠密分解（约 600 列），每个子步一次。P5.6 要做的就是复用分解、只在拖拽涉及的块上做。

**顺带改的**：`numerics::min_norm_step_within(j, f, rcond)`（`min_norm_step` 成了它取 `STEP_RCOND` 的包装）；`SparseJacobian::from_dense` 从 `#[cfg(test)]` 变成正式的 `pub(crate)`。

**没做的**（都在 ROADMAP）：翻面守卫与 `BranchFlip`（P5.3）；真正走到最远可行处（P5.5）；拖拽中的粗容差与帧预算（P5.6）；权重与优先级（P5.2）；3D 草图与装配的入口（`tangent_step` 本身是三族通用的，但没有接 gauge，3D 要先接上有效列）；圆弧端点作手柄（它是算出来的点，不是存储参数；debug 断言，release 按 `param_indices` 的既有回退拖圆心）。

### 13.3 P5.2 第一半：参数权重（完成，2026-11-15）

| 步骤 | 状态 | 说明 |
|---|---|---|
| P5.2 权重 | 完成（2D 拖拽） | 见下 |
| P5.2 约束优先级 | 下一项 | 用户要求跟权重一起做，分开提交 |

**用户定（2026-11-15）**：① 权重**按实体存在草图里**（`set_entity_weight(id, w)` / `entity_weight(id)`，默认 1.0，序列化带着），不是每帧传；② 口径是**连续权重 + 列缩放**（最小化 `Σ wᵢ² Δpᵢ²`，权重翻倍则同样收益下只动四分之一），不是三档角色；③ 约束优先级跟权重一起做（分两个提交）；④ 权重**只影响拖拽**，`solve()` 不读它（两份语料的解一字不改）。

**做法**：`tangent_step` 多一个 `weight_of: &dyn Fn(usize) -> f64`，在缩放坐标 `u = W Δp` 里走原来的那一步——雅可比的列除以权重（复用 `numerics::scale_columns`）、拉动行的系数同样除、最后把步乘回来（`numerics::unscaled`）。非正或非有限的权重算 1.0。`Sketch` 多一个 `weights: Vec<f64>`（缺项 = 1.0，`serde(default)`），每帧由 `param_weights()` 摊成逐参数的权重（实体的权重给它的每个参数）。

**先量（三个夹具，前两个都白量，记下来免得再走一遍）**：

1. **两个矩形串联**（角点的 x 是两个宽度之和）：统一权重下全部运动落在被拖的那个矩形上（第一个矩形动 0.0000），**把它加重到 20 倍也一字不变**。原因：改第一个矩形的宽度会把第二个矩形整体平移，那 8 个参数要按第二个矩形的权重付钱，所以"改第二个"永远便宜。
2. **A—P—B 两段等距拖 P**：统一权重下 A、B 各动 0.2131（对称），**任何权重都是这个数**。原因：P 被光标钉住后，A 只受 `|AP|` 约束、B 只受 `|PB|` 约束，两者互不耦合，各自的最小位移是几何给定的，没有可交换的余量。**权重要起作用，自由运动必须把实体耦合起来**（A 多动能让 B 少动）。
3. **刚性三角形拖顶点**（`|AP|`、`|PB|`、`|AB|` 都给定）：顶点被光标钉住后剩下的自由度正是绕顶点的旋转，它把两个底点的位移对换——这才是权重该做的选择。实测：统一权重 → 三角形平移，两个底点各走 1.0（= 顶点走的量）；A 加重 2 / 4 / 20 倍 → A 走 0.556 / 0.502 / 0.524，B 走 1.553 / 1.809 / 1.909（4 倍以上接近旋转的几何极限）；B 加重 4 倍 → 完全对称。

**没做"被拖实体自动降权"**（当时给用户的预览里提过）：实测它**没有任何效果**——刚性三角形里把顶点的权重改成 0.25 / 1 / 4，两个底点的位移一字不变；"线被一端拖着、另一端在圆上滑"的夹具同样一字不变。被拖点的坐标已经被光标钉住，整体缩放同一个实体的度量换不出东西；真正会被它影响的只有这个实体的**其它**参数（被拖圆的半径、被拖线的远端），让那些先让反而出乎宿主意料（"我拖圆心，半径怎么变了"）。所以默认全是 1.0，要偏向谁由宿主自己设。已在汇报里说明，用户可以推翻。

**验证**（`tests/drag.rs` 新增 3 个）：`a_heavier_entity_gives_less`（统一权重平移、4 倍给 0.50 / 1.81、镜像权重给镜像结果）、`a_weight_is_not_a_hold`（1e6 的权重下 `solve()` 的结果与不加权重逐位相同，并且照样拖得动）、`a_weight_has_to_be_positive_and_finite`（0、负数、NaN、∞、不存在的实体都被拒绝且不改动）。**变异检查**：去掉列缩放（保留行缩放与回缩）→ 2 个测试失败；去掉最后的回缩 → 2 个失败。
