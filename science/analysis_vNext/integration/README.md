# Versioned integration and annotation safeguards

本目录是2026-10-05本地整合审查的公开代码副本，只包含可执行modules、复现说明和无细胞身份的guard验证。完整审查证据、raw矩阵、参考模型和临床/逐celltables不在这个代码目录内。

`released_rna.py`明确区分full rawRNA与legacy scaledHVG featureaxes，严格验证CSR结构/非负精确整数；完整RNA library归一化10k/log1p后才选genes。旧denseHVG默认拒绝未知列序读取，需要独立验证source scale.data featureaxis；不能将H5内var index再次传入冒充独立source。

`integration_safeguards.R`按unique cell-ID附加embeddings与合并labels，v5/v4兼容RNA counts getter拒绝splitlayers混读，固定UMAPseed参数，并从实际counts/scale.data各自axes导出sourceSHA sidecar。`assemble_B_Plasma_vNext.R`修复旧base变量混用并拒绝覆盖已有outputRDS。

`reannotate_t1d_full_rna_pilot.py`提供fullRNA稀疏、5000cell分块、冻结model候选注释CLI。legacy/source标签只是审阅列，不是模型输入。High/Low相关模型不能当独立验证；其agreement和sigmoidscore不是准确率。候选不自动替换生产labels。

本地实际执行：T1D75746cells pilot；53938cells pairedHVG/fullRNA sensitivity；40cellschunk/model/column invariance；16cells×20690rawcoords与独立source0差异；T1D全部2000scaledaxis与archivedpostprocessedRDS一致，16×2000scaled值差≤5.33e-15。上述T1D证据不自动授权其它release对象。

执行入口见`RUNBOOK.md`；方法边界见`method_protocol.vNext.json`和`schema_sidecar_specification.md`。旧历史脚本弃用风险见`DEPRECATED_ANALYSIS.md`。UMAP只是展示，batch/biology验证应在高维及独立证据下完成。
