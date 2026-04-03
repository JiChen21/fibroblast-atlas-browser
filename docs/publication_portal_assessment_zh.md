# Fibroblast Atlas Browser 项目深度解读与论文发表级 Web Portal 改进建议

## 1. 项目定位与当前能力

该项目是一个基于 **Streamlit + AnnData** 的单细胞浏览门户，核心目标是：
- 让研究者在不写代码的情况下查看 UMAP、元数据分布与基因表达；
- 围绕心脏成纤维细胞跨疾病图谱做交互式探索；
- 兼顾一定规模数据（~70 万细胞）的浏览性能。

从代码结构看，应用主入口在 `app.py`，核心算法函数抽离在 `core.py`，并提供最小数据兼容校验脚本 `scripts/validate_h5ad.py`。

---

## 2. 架构拆解（按层）

### 2.1 数据层（AnnData 读取与契约）

- 使用 `H5AD_PATH` 指定数据文件，默认 `./data/FBs_adata.h5ad`。
- 关键依赖包括：
  - `adata.obsm["X_umap"]`（二维坐标）
  - `adata.obs["cell_type"]`
  - `adata.obs["condition"]`
- 若加载失败会退回 mock 数据（可通过 `STRICT_DATA=true` 禁止）。

**评价：**
- 优点：容错好，演示友好；
- 风险：生产环境若误用 demo fallback，可能掩盖数据加载失败。

### 2.2 业务/算法层（core.py）

已拆出的核心函数质量不错，属于“可测、可复用”的纯函数风格：
- 过滤：`apply_filters`
- 分层抽样：`stratified_sample_indices`
- 绘图索引策略：`choose_plot_indices`
- 元数据校验：`validate_core_metadata`
- 基因索引解析：`resolve_gene_index`
- 疾病-亚群组合计数与富集：`build_condition_subtype_counts`、`compute_roe`、`roe_symbol`

**评价：**
- 优点：核心逻辑与 UI 有一定解耦；
- 可提升：统计方法说明（Ro/e、伪计数、阈值）建议写进方法文档并在 UI 提供更强解释。

### 2.3 展示层（Streamlit 页面模块）

当前四个主模块：
1. Atlas Overview
2. Metadata Explore
3. Gene Query
4. Disease–subtype compare

并提供：
- 侧边栏筛选（condition/region/cell_type）
- 渲染模式（Auto/Downsampled/Full）
- 多图联动展示（UMAP + violin + dotplot + heatmap）

**评价：**
- 优点：模块化思路清晰，科研用户上手快；
- 可提升：缺“结果导出”、“会话复现链接”、“图注自动生成”等发表级功能。

---

## 3. 质量与工程化现状评估

### 3.1 测试

当前有 `tests/test_core.py`，覆盖了若干核心函数。

**优点：**
- 对分层采样、过滤、Ro/e 逻辑有基础保护。

**不足：**
- 缺少 `app.py` 页面级行为测试（至少 smoke test）；
- 缺少数据契约回归测试（字段缺失、稀疏矩阵格式变化、极端大数据）。

### 3.2 文档

`README.md` 已包含运行、校验、Docker、性能说明。

**不足：**
- 缺少“论文发表场景”专用规范：版本冻结、图表可追溯、数据引用说明、许可证/伦理声明。

### 3.3 性能

已有 WebGL + 下采样策略。

**不足：**
- 仍主要在“显示层”优化，缺少“统计层预计算与缓存”。

---

## 4. 如果作为论文发表配套 Web Portal，还需要增加什么？

下面按“发表刚需”优先级排序。

### P0（必须补齐）

1. **版本可追溯（Reproducibility）**
   - 在页面显式显示：
     - 数据版本号（例如 DOI / release tag / md5）
     - 应用版本号（git commit SHA）
     - 分析参数（筛选条件、抽样策略、阈值）
   - 提供“复制当前会话参数”按钮（JSON / URL query）。

2. **图表导出与可复现图注**
   - 一键导出 PNG/SVG/PDF；
   - 自动生成图注模板（数据来源、样本量、统计方法、日期、版本号）。

3. **数据与方法透明性页面**
   - 新增 `Methods / Data provenance` 页面，明确：
     - 数据收集标准与纳入排除条件
     - 预处理流程（QC、归一化、批次校正）
     - 注释策略与参考文献

4. **下载与引用规范**
   - 增加 `Download` 页面：
     - 可下载聚合矩阵、marker 表、元数据字典
     - 明确 Citation 文本、许可证、使用限制。

### P1（强烈建议）

5. **更完整的基因检索能力**
   - 大小写不敏感；
   - Symbol / Ensembl ID 双向映射；
   - 同义名（alias）支持与候选提示。

6. **统计严谨性增强**
   - 在 Disease-Subtype 模块中增加显著性检验与多重检验校正（如 FDR）；
   - 提供 effect size + CI + p-value 的表格下载。

7. **多层次导出**
   - 当前图（可视化）导出；
   - 当前筛选后细胞索引导出；
   - 聚合统计导出（condition × cell_type × metric）。

8. **可访问性与国际化**
   - 增加英文/中文切换；
   - 颜色映射考虑色盲友好；
   - 添加键盘导航与替代文本。

### P2（长期价值）

9. **身份与权限分层（若涉受限数据）**
   - 公共摘要可见，受限原始数据需要授权。

10. **使用分析与稳定性监控**
    - 记录匿名访问指标、页面延迟、错误率；
    - 引入健康检查与告警。

11. **可部署性增强**
    - 提供 docker-compose / k8s 清单；
    - 反向代理与 TLS 最佳实践文档。

---

## 5. 建议的“发表级门户”信息架构

建议把导航扩展为：
1. Home
2. Atlas Explorer（现 Metadata Explore）
3. Gene Explorer（现 Gene Query）
4. Disease Comparator（现 compare）
5. Methods & QC
6. Download & Citation
7. Changelog / Version

这样可同时服务：
- 审稿人（看方法和可重复性）
- 研究者（做探索和下载）
- 编辑部（核查发布规范）

---

## 6. 安全与合规（发表时经常被忽略）

若包含潜在敏感临床元数据，需补：
- 去标识策略说明；
- 数据访问分级策略；
- GDPR/HIPAA（视数据来源）合规声明；
- 第三方数据许可证继承说明。

---

## 7. 分阶段落地路线（建议）

### 第一阶段（1–2 周）
- 增加版本展示（数据/代码/参数）；
- 实现图表导出 + 当前参数导出；
- 新增 Methods & Citation 页面（先静态文案）。

### 第二阶段（2–4 周）
- 补全基因 ID 映射与候选提示；
- Disease compare 增加统计检验与 FDR；
- 加入下载中心（聚合结果 + 字段字典）。

### 第三阶段（持续）
- 引入监控、日志、异常告警；
- 完善 CI/CD 与回归测试；
- 对接 DOI 与正式发布流程。

---

## 8. 结论

这个项目已经具备“可用研究浏览器”的核心骨架，尤其在交互模块划分和大规模可视化方面起点不错。若目标是“论文发表配套门户”，重点不是再堆叠图，而是补齐 **可追溯性、导出能力、方法透明性、统计严谨性与发布规范**。这几项到位后，平台的学术说服力会显著提升。
