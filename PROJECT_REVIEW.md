# 项目审阅与改进建议（Fibroblast Atlas Browser）

基于当前仓库中的 `app.py`、`README.md` 与 `scripts/validate_h5ad.py` 进行的静态审阅。

## 总体评价

这个项目已经具备一个可用的单细胞浏览门户雏形：
- 模块划分清晰（Home / Metadata / Gene Query / Data Dictionary）。
- 大数据量场景下考虑了 downsampling 与 WebGL 渲染。
- 有基础的 H5AD 校验脚本，便于上线前快速检查。

但如果目标是**长期维护 + 稳定生产部署 + 可扩展分析能力**，还存在一些明显短板。

---

## 主要不足（按优先级）

## 1) `app.py` 过于庞大，架构耦合较高（高优先级）

当前大部分业务逻辑、数据读取、绘图、UI 状态都集中在单文件中，后续迭代风险高：
- 回归影响面大；
- 不利于单元测试；
- 新功能（比如差异分析、导出、注释系统）会快速恶化可读性。

**建议：**
- 按职责拆分为：
  - `data_io.py`（加载/校验/表达矩阵抽象）
  - `filters.py`（筛选与采样）
  - `plots.py`（UMAP/violin/dotplot）
  - `pages/`（各模块页面渲染）
  - `config.py`（常量与环境变量）
- 给关键纯函数补单元测试（尤其是 `stratified_sample_indices`、`resolve_gene_index`、`apply_filters`）。

## 2) 错误处理策略偏“静默降级”，生产可观测性不足（高优先级）

`load_adata` 在加载失败时默认回落到 mock 数据，这对 demo 很友好，但在生产场景可能掩盖真实故障。

**建议：**
- 保留 demo fallback，但默认在非开发环境禁用（例如 `APP_ENV=prod` 时强制 `STRICT_DATA=true`）。
- 引入结构化日志（至少 `logging`），将关键异常、筛选规模、渲染模式记录出来。
- 在 UI 顶部展示当前数据源类型（real/demo）、加载时间、数据版本标识。

## 3) 数据契约不完全一致（中高优先级）

README 和校验脚本要求的字段较多，但 app 实际核心依赖主要是 `cell_type` 和 `condition`。这会导致：
- 使用者不知道“必须字段”与“可选字段”的边界；
- 数据接入时反复试错。

**建议：**
- 明确分层：
  - 必须字段（运行阻断）
  - 功能增强字段（影响筛选器/图表）
  - 可选字段（仅展示）
- 在 `scripts/validate_h5ad.py` 支持 `--profile minimal|recommended|strict`。

## 4) 性能策略可再前移到“数据层”（中优先级）

目前主要通过前端点数控制改善交互；对于 70 万级细胞仍可能卡顿。

**建议：**
- 预计算并缓存常用统计（condition/cell_type 层级表达汇总）；
- 对高频查询基因可做轻量缓存（如 session 级 LRU）；
- 针对分类字段提前转 `category`，减少内存与比较成本；
- 在 README 给出推荐硬件和目标响应时间。

## 5) Gene Query 的可用性与生信习惯兼容性一般（中优先级）

当前基因解析仅直接匹配，别名/大小写/ENSEMBL 场景不够友好。

**建议：**
- 增加基因匹配策略：大小写不敏感、symbol ↔ ensembl 双向映射；
- 查询失败时给出 top-k 相似候选；
- 支持一次输入多个基因并输出 heatmap / panel 对比。

## 6) 文档对“部署与运维”覆盖不足（中优先级）

README 已涵盖本地运行与 Docker，但缺少运维关键内容：
- 资源上限与并发建议；
- 反向代理（Nginx）示例；
- 健康检查与监控；
- 数据更新流程（版本化、回滚）。

**建议：**补充 `docs/deployment.md` 和 `docs/data_contract.md`。

## 7) 依赖管理较宽松，复现性一般（中优先级）

`requirements.txt` 仅给下限，长期可能遇到上游 breaking change。

**建议：**
- 使用 `pip-tools` 或 `uv` 生成锁定文件（含 transitive dependencies）；
- CI 中增加“冷启动安装 + smoke test”。

## 8) 测试与 CI 缺失（高优先级）

当前仓库未见自动化测试或 CI 配置。

**建议最小闭环：**
1. `pytest`：核心函数单测；
2. `ruff/black`：静态检查与格式化；
3. GitHub Actions：安装依赖 -> 运行测试 -> 运行 `validate_h5ad.py`（mock 路径场景）。

---

## 可以快速落地的 2 周改进路线图

### Week 1（稳定性）
- 拆出 `data_io.py` / `filters.py` / `plots.py`。
- 给 3~5 个核心函数加单测。
- 引入 `logging` 并增加关键运行指标日志。
- 新增 CI（lint + unit test）。

### Week 2（可用性与性能）
- Gene Query 增加大小写不敏感与候选提示。
- 预计算并缓存 group-level 统计。
- 完善部署文档与数据契约文档。

---

## 结论

这个项目当前定位非常适合“研究组内部演示与探索”。
若要面向更广泛用户稳定运行，优先做三件事：
1. **模块化重构**（降低维护成本）
2. **测试 + CI**（降低回归风险）
3. **数据契约与运维文档补全**（降低接入成本）

以上三项完成后，再继续做高级分析功能，会更稳、更快。
