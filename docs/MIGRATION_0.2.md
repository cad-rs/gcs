# 从 0.1.0 迁移到 0.2.0

写给宿主（rcad 的草图层、3D 草图层、装配层）。每一项变化的完整说明在 [CHANGELOG.md](../CHANGELOG.md)；这里只按“你需要做什么”来排，分三类：编译会报错的、不报错但行为变了的、可以顺手换掉的旧做法。

我没有读到宿主的代码，下面写“如果宿主……”的地方都是按库这一侧的改动推出来的，请对照实际用法核对。

## 1. 编译会报错的

| 现象 | 原因 | 改法 |
|---|---|---|
| 读 `constraint.constraint_type`、`component1`、`point1`、`fixed_orientation` 等字段处报“no field” | 刚体 `rigid3d::Constraint` 改成了 `{ id, kind }`，`kind` 是带数据的 `ConstraintKind` 枚举 | 只需要类型和零部件时用 `constraint.constraint_type()`、`component1()`、`component2()`；需要数据时 `match constraint.kind { ConstraintKind::Distance { distance, .. } => … }`。**构造函数签名没有变**，`Constraint::coincident(…)` 等调用不用动 |
| `ConstraintSystem { constraints, num_components }` 结构体字面量报缺字段 / 私有字段 | `ConstraintSystem` 多了记录编辑状态的私有字段 | 改用 `ConstraintSystem::with_components(n)` 再 `add_constraint` |
| `rigid3d::SolveOptions { max_iterations, tolerance, damping, lm_lambda }` 报缺字段 | 多了 `on_failure` 字段；类型改名 `PoseSolveOptions`（旧名是弃用别名） | 写成 `PoseSolveOptions { tolerance: 1e-9, ..Default::default() }` |
| 手写 `PoseSolveResult { … }`、`SolveResult { … }`、`DofReport { … }`、`SpaceDofReport { … }` 字面量报缺字段 | 分别多了 `status`；`status`；`rank` 和 `structural_dof` | 补上字段。只读不建的代码不受影响 |
| 闭包或 `()` 块以 `sk.add_constraint(…)` 结尾时报类型不匹配 | 三个容器的 `add_constraint` 现在返回约束下标 `usize` | 在调用后加分号 |
| 用了 `numerics::newton_solve_with_jacobian` | 它的一个参数类型本来就是库内私有的，现在函数也改成私有 | 用 `newton_solve` / `newton_solve_with` |
| 对 `SolveStatus` 做穷尽 `match` 的地方报缺少分支 | 新增了状态 `ConvergedToPrecision`（见 2.1） | 把它和 `Converged` 归为一类处理；判断成败请用 `converged` / `success` 字段 |
| 用了常量 `numerics::LINE_SEARCH_RELAX`、`LINE_SEARCH_STEPS` | 减半线搜索换成了信赖域，这两个常量随之删除 | 没有对应物；它们只是内部步长控制的参数 |
| 宿主自己的类型对 `glam` 向量派生 serde 时报 `Serialize` 未实现 | 0.1.0 无条件打开了 `glam/serde`，宿主可能在不知情的情况下依赖了它 | 打开本库的 `serde` 特性，或在宿主自己的 `Cargo.toml` 里给 `glam` 开 `serde` |
| `-D warnings` 的构建在 `DOFAnalyzer` / `SolverDiagnostics` / `rigid3d::SolveOptions` 处报弃用 | 三者标了 `#[deprecated]` | 短期加 `#[allow(deprecated)]`，再按第 3 节换掉 |

`solve` / `solve_components` / `solve_space` / `newton_solve` 的参数从 `&mut Vec<f64>` 变成了 `&mut [f64]`。传 `&mut vec` 的调用会自动解引用，一般不用改。

## 2. 不报错，但行为变了

这一类要重点看，因为编译器不会提醒。

### 2.1 求解失败后留下的东西

- **失败不再留下最后的迭代点。** 三族现在都一样：没收敛的连通分量（刚体是簇）回到求解开始时的状态，收敛了的保留结果。如果宿主在失败后自己做了“恢复快照”，现在可以去掉；如果宿主**依赖**失败后的迭代点（比如拿它显示“最接近的位置”），要显式要回来：2D / 3D 用 `SolveOptions { on_failure: FailurePolicy::KeepLast, .. }`，刚体用 `PoseSolveOptions { on_failure: FailurePolicy::KeepLast, .. }`。
- `SolveResult::residual` 和 `PoseSolveResult::error` 是**留下来的那个状态**的残差：默认策略下，失败时它是起始状态的残差，不是迭代到最后的残差。
- 失败原因在 `status` 里：`MaxIterations`、`Stalled`、`Singular`、`Diverged`。以前 NaN 会让消元 panic，现在报 `Diverged`。
- **解不出来的系统现在很快结束，报 `Stalled`。** 冲突的约束组以前会在最小二乘折中点附近转到迭代上限（默认 100 次）才返回；现在一旦任何长度的步都不能再降低残差就停（测试里的冲突草图不到 20 次迭代）。交互中“加了一个矛盾约束”的反馈因此快得多。如果宿主用 `status == MaxIterations` 判断“这是个冲突”，要改成看 `Stalled`（更可靠的做法是看 `diagnose()` 的 `conflicting`）。`MaxIterations` 现在只表示“还在下降，但迭代次数用完了”。
- **远离原点的模型现在能解，状态是 `ConvergedToPrecision`。** 容差是绝对数（默认 1e-10），而十亿量级的坐标只能存到约 1e-7，方程不可能满足到 1e-10；0.1.0 在这种模型上一律失败。现在每个方程按“容差”和“自身数值的 4 个舍入单位”中较大的一个来要求，达到后报 `ConvergedToPrecision`，`converged` 为真，但 `residual` 可能大于容差。宿主如果用 `residual < 某阈值` 自行判断成败，在大坐标模型上会误判，请改看 `converged`。原点附近的模型不受影响。**反方向只解决了一半**，见下一条。
- **2D 的平行、垂直、角度、等长换了残差写法（3D 的等长线也是）。** 解集不变（`Angle` 仍然在 `angle + π` 处也成立），残差的数值变了：平行、垂直、角度原来是两条线向量的叉积、点积和 `sin·点积 − cos·叉积`，是“长度 × 长度”；现在是单位方向上的 `sin θ`、`cos θ`、`sin(angle − θ)` 再乘两条线长度的几何平均，等长是 `|d1| − |d2|`，全都是长度量纲，和其他 2D 方程一致。宿主如果读单个约束的残差值（比如拿它显示“误差”或设阈值），数值会不同：角度类现在大约是“线的端点离满足关系的位置有多远”。效果：完全约束的模型在 1e-6 到 1e6 的尺寸上都解到舍入精度（以前 1e-5 大小的三角形会在顶点偏差尺寸 0.12% 时报收敛）。
- **2D 里小于 1 个单位的分量，容差按它自己的尺寸算。** 分量的尺寸取它的约束所引用实体的包围盒的长边（圆和圆弧按整圆算）；小于 1 时，方程按“残差 ÷ 尺寸”对容差，所以 1e-4 大小、还有自由度的三角形现在和 1 个单位大小时一样，长度误差在尺寸的 1e-15 量级（以前容差是绝对的 1e-10，这种分量拿不到收敛后的精修，误差是尺寸的 2.7e-8）。1 个单位及以上的分量照旧。宿主可见的变化：小分量上报告的 `residual`（仍是模型单位）会比容差小得多；如果宿主自己用“`residual` 接近容差”来判断什么，在小模型上会看到更小的数。2D 的 `diagnose()` 同样按分量尺寸判断：以前阈值有一个绝对下限 1，1e-6 大小的三角形上一个要求多出 0.2% 的距离被报成 `partially_redundant`，现在报 `conflicting`；如果宿主在小模型上靠 `partially_redundant` 做过什么处理，现在这些会落到 `conflicting`。2D 草图不再需要“特征尺寸保持在 1e-4 以上”。**3D 草图和装配还没有做**：3D 的角度类方程是纯数，刚体是长度和角度混合，这两族仍用绝对容差，小于 1e-4 的 3D 模型在还有自由度时只满足到 1e-10 个模型单位，暂时仍请让特征尺寸保持在 1e-4 以上。
- **零长线段上的 `LineLength` 从 `Diverged` 改报 `Stalled`。** 以前它的导数是 0/0，求解立刻以 `Diverged` 结束；现在和“两个重合点之间的距离”一样，梯度为零，报 `Stalled`，几何放回原样。宿主如果用 `Diverged` 识别“画了个零长线段”，改看 `Stalled`（或者在加约束前自己检查零长度）。零长线段、重合点这类起点目前都解不动（没有方向），见 `ROADMAP.md` 的 “Degenerate starts”。
- **病态的步改用 QR，`λ` 只作用在良态的步上。** 牛顿步原来一律解正规方程 `(JᵀJ + λI)·h = −Jᵀf`；现在当分解后的最小主元小于最大主元的 1e-8（秩亏，或者条件数平方后丢了一半有效位）时，这一步改用 J 的列主元 QR 求按秩截断的最小范数解，不加 `λ`。宿主可见的变化：① `numerics::LAMBDA` 和刚体的 `PoseSolveOptions::lm_lambda` 只对良态的步起作用——如果宿主把 `lm_lambda` 调大当作 LM 阻尼来用，这在病态的步上不再生效；② 接近奇异的几何（比如极细的三角形、接近相切）在紧容差下现在能收敛，以前会跑满迭代；③ 其余的步逐位不变，语料指标与基准持平。
- **2D 求解先消掉平凡方程，迭代次数和远处起点的选解会变。** 重合、固定点、存储点之间的 `DistanceX` / `DistanceY`（圆弧的派生端点不算）、水平 / 竖直、等半径、同心、半径、直径，这些方程都只是“一个参数等于另一个加常数”或“一个参数等于某值”。现在每个分量先按它们把参数并成类，把每个参数放到类说的位置（已知的类放到值上，未知的类放到各成员所说位置的平均值，也就是离当前几何最近的、让这些方程成立的位置），再只对剩下的方程迭代。宿主可见的：① `SolveResult::iterations` 变少（语料里从画好的状态出发的完全约束例子合计 5 020 → 2 825），全由这类方程组成的分量 0 次迭代；宿主如果拿迭代次数设过阈值或写过断言，要重新看。② 被这类方程连起来的坐标结束时**严格相等**，而不只是在容差内相等。③ 从远处起点重解时选到的解可能不同：语料的大幅漂移探针里，原来 256 次重解有 71 次回到另一个分支（多数是镜像的矩形），现在 5 次；小幅漂移回环的最差偏差 1.59e-4 → 6.79e-5。④ 欠约束的草图可能停在稍有不同的位置：开放的距离链在漂移 0.05 后与原来最多差 1.6e-3，漂移 3 后离起点的移动量少 2–9%。冲突的系统不受影响：消元后解不出来就把参数放回去、按全部参数重解，结果和以前逐位相同；`diagnose()` 仍然看完整的系统，冲突和冗余照常带原约束编号报告。NaN / Inf 坐标仍报 `Diverged`（不会因为某个固定点约束而被“修好”）。3D 草图和装配没有变。
- **“没有自由参数”不再等于收敛。** 一个分量的参数全被固定而方程不成立时，0.1.0 报 `converged = true`，现在报 `false` / `Stalled`。如果宿主把“全固定”当成快速通过的路径，这里会开始报失败——这些失败是真的。

### 2.2 自由度的数字

- `Sketch::dof()` 和 `SpaceSketch::dof()` 现在是秩口径，和 `diagnose().dof` 是同一个数，**不会为负**。宿主如果用 `dof() < 0` 判断过约束，这个判断从此永远为假：过约束要看 `diagnose()` 的 `redundant` / `conflicting`。原来的“参数数 − 方程数”还在，叫 `structural_dof_estimate()`。
- 自由圆柱报 5 个自由度（原来 7 个）。
- `DOFAnalyzer` 和 `SolverDiagnostics` 的结果变了（见 3.3），而且各自现在要做一次求解加一次分解；原来只是扫一遍约束表。如果宿主每帧调用它们，请改成在装配变化后调用一次。

### 2.3 装配求解

- `ConstraintSolver::solve` / `solve_from` 现在逐簇求解。`success` 的含义是“每个簇自己的残差范数都低于容差”，`iterations` 是各簇之和。
- 带完整位姿的 `Fixed`（`fixed_with_orientation`）是**接地**：刚体被直接放到那个位姿，不再参与迭代；热启动位姿里它在别处也会被纠正过来。只固定位置的 `Fixed` 仍是普通方程。
- `Fixed`（带姿态）和 `Lock` 的姿态误差现在与四元数的符号无关，总是取短的一边（不超过 π）。0.1.0 里同一个位姿用 `q` 和 `−q` 写会报出不同的 `error`（相差可以接近 2π），诊断也可能因此误报冲突。宿主如果为了绕开这一点在传入前把四元数统一成 `w ≥ 0`，现在不必再做；如果宿主拿 `error` 做过阈值判断，原来偶发的“误差约 6.28”不会再出现。
- `Parallel` 现在是 2 个方程、`Concentric` 是 4 个（原来 3 和 6）。解出来的位姿在容差内与原来一致，但不是逐位相同。**诊断结果会明显变干净**：0.1.0 把每一个 `Parallel` / `Concentric` 都报成冗余，哪怕它是装配里唯一的配合。如果宿主为此过滤过诊断结果（例如忽略这两类配合的冗余报告），这层过滤现在应当去掉，否则会漏掉真正的重复配合。

### 2.4 3D 草图的样条约束

`PointOnSpline` 和 `SplineControlPointSymmetric` 的导数原来是错的：前者在点不在曲线上时会让 debug 构建 panic、release 构建带着不准的导数迭代；后者在镜像直线自由时既不会移动直线，也不会把直线和样条分进同一个分量。现在都是精确导数。宿主如果因此避开过某些用法（先把点投到曲线上再加约束、把镜像线固定住），可以不必再绕。

### 2.5 已抑制 / 已删除的约束

三个容器的分解、求解、诊断、可动性都只看生效的约束。这只影响用了新编辑接口的代码。

## 3. 可以顺手换掉的旧做法

### 3.1 改尺寸不必重建

0.1.0 没有编辑接口，改一个尺寸要重建整个草图再冷启动。现在三个容器都可以原地改：

```rust
let width = sk.add_constraint(Constraint::LineLength { line, length: 40.0 });
sk.set_constraint_value(width, 55.0);
sk.solve();                                   // 从当前几何热启动
```

下标是稳定的：删除只把槽位标成“已删除”，后面的下标不移位。因此宿主如果自己维护了“库内下标 ↔ 宿主对象”的映射并在删除后重排，重排可以去掉。还没有回收槽位的 `compact()`，长时间编辑后槽位只增不减。

`SpaceSketch` 和 `ConstraintSystem` 的列表仍是公开字段，直接往里 push 仍然可用，只是 `topology_revision()` 不会知道；直接从中间删除会让下标移位，请改用 `remove_constraint` / `remove_entity` / `remove_component`。

装配系统本身不存位姿，所以：

```rust
system.set_constraint_value(index, 2.5);
let next = solver.solve_from(&system, last.positions.clone(), last.rotations.clone());
let now = system.measure(index, &next.positions, &next.rotations);
```

### 3.2 从动尺寸、量取、结构缓存、合法性检查

- 从动（参考）尺寸 = `set_constraint_enabled(index, false)` + `measure(index)`：约束留着、不参与求解，随时可以读它当前的值。
- `topology_revision()` 只在结构变化时变（增删、启用 / 抑制、固定），改尺寸值和移动几何不会变。宿主按结构缓存的东西（连通分量、UI 上的约束图）可以用它判断是否失效。
- `validate()` 在求解前查出会让求解无意义或 panic 的模型错误。两个草图容器查：引用了不存在或已删除的实体、实体类型不对、数值非有限，3D 另查样条控制点下标越界（`SketchError` / `SpaceSketchError`）；装配查：零部件不存在或已删除、数值非有限、方向或旋转长度为零、约束两端是同一个零部件（`AssemblyError`）。0.1.0 在这些情况下要么 panic，要么静默读错参数。`validate()` 是宿主主动调用的检查，求解本身不会替你调用它。

### 3.3 装配诊断

| 旧 | 新 |
|---|---|
| `DOFAnalyzer::new().analyze(&system).total_dof` | `diagnose_assembly(&system, &positions, &rotations).dof`（在解出的位姿上调用） |
| `DOFResult::free_dofs` / `component_dof` | `assembly_mobility(&system, &positions, &rotations)`：每个刚体六个坐标各自能不能动 |
| `SolverDiagnostics` 的 `Conflict` / `Redundant` 问题列表 | `diagnose_assembly` 的 `conflicting`、`redundant`、`dependency_groups`（都是配合 id）。文案由宿主自己生成，可以本地化 |
| `SolverDiagnostics` 的 `OverConstrained`、`CyclicDependency` | 不再报告。前者是计数口径的误判来源；后者把正常的闭环配合当成了缺陷 |

`assembly_mobility` 报的是“哪些坐标参与了某个剩余运动”，不是“有几个独立运动”：沿对角线滑动的刚体有两个平移坐标可动，但只有一个自由度。总自由度以 `diagnose_assembly().dof` 为准。

诊断与位姿有关：冗余要在解出的位姿上读，冲突在求解失败后留下的位姿上读（默认策略下就是起始位姿，矛盾在那里的线性化里是可见的）。

### 3.4 序列化

刚体的 `Constraint`、`ConstraintKind`、`ConstraintSystem` 现在也在 `serde` 特性下可序列化。0.1.0 存下来的 2D / 3D 草图仍然可以读入（新增的编辑状态字段缺省为“全部生效”）。

一个与本库无关但容易踩到的点：`serde_json` 默认的浮点解析可能有 1 ulp 的误差，JSON 往返后的位姿不保证逐位相同；需要逐位一致时请打开 `serde_json` 的 `float_roundtrip` 特性。

### 3.5 多个连通分量并行求解

新的可选特性 `rayon`：一次求解里的连通分量分到多个线程上（rayon 的全局线程池；如果调用发生在宿主自己的 rayon 池里，就用那个池）。结果与顺序求解逐位相同，宿主不需要改任何调用。只有可分的工作足够多时才会并行：相邻分量按每 96 个方程一个任务打包（不超过线程数，各任务方程数大致相等），总共不到 192 个方程的系统照旧在调用线程上顺序求解；一个大分量不会被拆开。典型的“一张草图一个大分量”不会变快；几十个以上互不相连的小草图、或挂在接地件上的几十个零件才有收益（数字见 CHANGELOG）。不打开这个特性时依赖照旧只有 `glam`。

## 4. 这一版没有动的

以下已知问题仍在，迁移时不要指望它们已经解决（都在 [CATCHUP_PLAN.md](CATCHUP_PLAN.md) 里排了阶段）：

- 2D 线–圆相切有符号、2D 角度分不清 `θ` 与 `θ+π`、刚体 `Angle` 在 0° / 180° 没有梯度（起始时两方向平行的角度约束解不动）。
- 线性代数仍是稠密的，单个分量约 800 个自由参数时一次求解约 0.36 秒。
- 没有拖拽求解，没有解分支选择。
- 刚体求解热启动与冷启动耗时相同。
