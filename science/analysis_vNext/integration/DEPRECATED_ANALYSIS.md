# 历史分析脚本的兼容与替代

以下历史脚本保留用于追溯，不应作为新分析的入口：

- `legacy/AID/T_NK/01.QC.Mult/01.QC.Mult.R` 的 benchmark 部分在 UMAP 的二维图上算 LISI，并以每种算法自己生成的 cluster 作为 cell-type label。这个量只能表示二维图和该方法自己的聚类结构，不能证明生物信号保留或选出最佳整合方法。
- `legacy/AID/B_Plasma/02.scAnn/Mult_04/Mult_04.R` 读入的 base object 名为 `scRNA.Mult2`，随后却 `subset(scRNA.Mult1, ...)`。独立运行会找不到 `scRNA.Mult1`；交互式残留对象又可能造成错误来源。
- 历史多算法 reductions/cluster assignment 只调用 `identical(...)` 显示结果而不拒绝错序。新流程必须按 cell ID 匹配，每次 attach 前检查缺失、重复和顺序。
- `T_NK/02.scAnn/Mult_01/Mult_01.R` 的 CellTypist 输入先截取 HVG counts、再重新归一化；不能作为新注释输入规范。整合中屏蔽 TCR/IG 可用于避免克隆信号主导图形，但注释输入应重新读取完整 RNA counts，而不是使用被屏蔽的整合 feature 集。
- `cross_disease_analysis/06_integration_quality.R` 将 raw iLISI 解释成“越高越好”，忽略 condition/tissue/lineage 的研究混杂；其全距离 silhouette 分配过大，且 `geom_text_repel()` 未显式载入 `ggrepel`。新诊断固定条件、组织和谱系，限额抽样，在高维 latent 空间计算；结果注明所用空间与覆盖范围。

`integration_safeguards.R` 替代危险的 positional attach、未命名 base object 合并和 HVG-library normalization；已经用错序、缺失和重复轴以及完整库归一化执行验证。`released_rna.py` 是 H5 schema/feature-axis reader。`audit_metadata_and_latent.py` 是已执行的可复现诊断。`reannotate_t1d_full_rna_pilot.py` 是已经执行的完整 RNA 重注释 pilot，输出候选标签和 review status，不自动把新模型结果覆盖原标签。

未来整合参数必须冻结 seed、软件/参考模型版本、样本与 study keys、HVG 集、维度、cluster resolution 和 cell-order provenance。默认不移除细胞周期；若比较 cell-cycle regression，报告对 Cycling population 的保留影响。study/tissue/disease 无独立重复时，不能用“全部混匀”作为优化目标，也不能把不可识别的混杂重新解释成疾病机制。
