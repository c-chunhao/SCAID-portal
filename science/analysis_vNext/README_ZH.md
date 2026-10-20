# SCAID analysis vNext — audited 5 October 2026

这是新的分析入口。历史R脚本、原RDS/base H5、旧文章与图件保留用于溯源；旧细胞marker交集与 pooled CellChat aggregate 不能用于疾病推断。

## 当前可直接执行的入口

- `statistics/recompute_atlas_composition.py`: 明确CLI，从全部release metadata及两条integration exact joins重算；source-positive但无integrated members记NA；真实观察到的亚型0保留。输入只读，输出独立目录。
- `integration/released_rna.py`: 原始counts/完整RNA normalization与scaled HVG轴分开；分块稀疏读取并检查数值及索引。完整library normalization在gene selection之前完成。
- `integration/integration_safeguards.R`: 基于unique cell ID附加embedding、严格metadata join、完整counts normalization与Seurat版本兼容getter；错误轴/missing ID即拒绝。
- `qc/sparse_export_guard.py`: 完整CSR块导出与独立shape检查，不遗漏singleton尾行。
- `release_scientific_gate.py`: typed evidence、hex SHA、技术capture/clinical grouping/annotation/duplicates gate；sample标签不自动成为donor。
- `new_data/review_aps_multimodal.py`: 配对protein/RNA descriptive审查、内部reference flags和source 500-gene sensitivity；不自行覆盖注释/QC。
- `new_data/audit_new_data_sources.py`: 整数计数、CSR shape/indices、完整MatrixMarket记录及feature/barcode复核；整数representation不是源counts来源证明。

Integration详见该目录README/RUNBOOK；执行者须在本地另行保存科学证据、输入路径、版本及SHA；患者数据不随源码发布。其余CellChat/scores/statistical replay脚本是冻结审查输入的证据重放，公开发布时另提供env根目录/resolver，不把本机绝对路径当通用安装方式。

## 执行过的验证

- 48个effective RNA文件的全部2,572,480条count rows与3,646,001,623非零值检查。
- composition六层从原metadata独立重跑：563 T/NK /537 B-plasma defined fractions；148未评估lineage units不填0。
- 48 CellChat对象、144 score H5全数值/结构/身份审查。
- 16 T1D cells×20,690 RNA genes独立source counts精确匹配；16×2000scaled genes与对应历史Seurat对象一致。
- 真正per-capture scDblFinder PSO pilot、前后QC顺序敏感性、源barcode/counts保存；APS配对ADT CLR与Seurat独立比较。
- sparse exporter 0/1/4999/5000/5001/10000/10001边界和人为缺尾检测；raw reader fractional/CSR损坏fail-closed。

## 不可省略的科学条件

1. RNA先从混合特征矩阵中分离；抗体或HTO不进入RNA library size。source filtered counts不等于unfiltered droplets。
2. 仅使用独立验证的technical capture切分doublet；保留所有原barcodes/counts及QC flags。
3. 有source-major lineage但无integration成员应NA。partial retained integration proportions也不等于完整组织绝对丰度。
4. 推断按实际independent source individual进行；须同时满足病例对照、临床key、full-rank和covariate overlap。
5. 2D UMAP mixing/自身cluster标签不是独立biology validation；High/Low同reference预测不是gold standard。
6. cell-type markers不是disease DEGs；评分不是RNA counts/P值；不同CellChat DB snapshot概率不比较疾病强度。
7. 新59候选不自动加入48对象总量，历史缺失数据不补造。全部fine labels尚不能宣称经过盲法准确性验证。
