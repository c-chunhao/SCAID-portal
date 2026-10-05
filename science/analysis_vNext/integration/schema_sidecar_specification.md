# RNA 矩阵与 feature axis 规范

`names_obs` 是 cell axis。`names_var` 与 `var/rawvar/_index` 是完整 raw RNA feature axis。`assay/RNA/layers/rawdata` 的 CSR blocks 对应完整 raw RNA counts，顺序必须与前述 axes 一致。

`assay/RNA/layers/data` 的现有 dense matrix 只有 2000 个 HVG features。它的列名必须读取 `var/var/_index`，不能使用 `names_var` 的前 2000 个 gene。实物矩阵 shape 为 cells×HVG，而旧 `shape` attribute 写成 HVG×cells；读取时以真实 dataset shape 及明确 feature index 校验，保留倒序 attribute 作为兼容诊断。名称 `data` 不足以证明它是全基因 normalized RNA；有负值的矩阵更不能当 raw counts 或 CellTypist log1p 输入。

新的 `SCAID-RNA-sidecar-1.0` JSON 明确记录每个 block 的 cell offset、真实 cells×features 形状、feature index、matrix semantics、encoding、倒序 attribute 状态和允许的分析用途。`released_rna.py` 提供可复用 reader，并在 feature 轴不确定时拒绝猜测。

还需区分“有 2000 个轴名”和“每列的 gene 次序已独立验证”。easySCFr 导出代码选择 `rawvar[rownames(rawvar) %in% VariableFeatures(sce), ]` 生成该层 var metadata，未显式按 `rownames(scaleData)` 重排。T1D 的真实 archived postprocessed RDS（20241116/00.data/type_1_diabetes_mellitus/syn53641849/PBMC）已成功读取：RNA scale.data 的全部 2000 rownames 与 H5 的 declared axis 逐项完全一致；16 cells×2000 scaled values 最大差 5.33e-15，符合文本 round-trip 精度；同16 cells×20,690 genes raw counts 与独立 official source 完全一致。因此该 T1D dense axis 已有独立实物支持，未发现错配。历史实际 exporter-call 关联没有存档；其余对象没有由此自动验证。旧 dense/HVG reader 默认拒绝读取，T1D 可显式传入这份独立核验的 `T1D_export_candidate_2_scale_feature_names.txt` 授权读取。Raw counts 正常可读。

新 writer 必须分别写 `raw_feature_names = rownames(raw_counts)` 与 `scaled_HVG_feature_names = rownames(scaleData)`；对应 metadata 必须用 `match()` 按该 axis 重排，并 `identical(rownames(metadata), rownames(matrix))` 校验，不能仅用 membership filtering。每层记录 counts/lognormalized/scaled 的语义，不再把 scaled RNA 仅名为 data。

CellTypist 输入先对完整 raw counts 按全 RNA library 归一化到 10000，再 log1p；若需删选模型 genes，在归一化后删选。不能先截 HVG，再把 HVG library 重新归一化到 10000。单个 CSR block 可独立归一化和分类；不必将百万 cells 转成 dense matrix。

规范只消除错误读取风险，不等于现有网站已发生基因错配。当前网站展示预生成图件，不直接通过这个 dense layer 进行表达矩阵查询。旧文件保持原样，新的 reader 和 JSON sidecar 用于可复现分析和未来导出。

复现：

```bash
python released_rna.py /path/to/released_RNA.h5 --out schema.json
```

返回的 counts CSR 严格检查非负整数值、rows 覆盖和 full gene order。`iter_full_log1p` 保留稀疏矩阵；`read_dense_hvg` 单次最多读 10000 cells，默认关闭，需要 `verified_source_feature_names`（已确认相应 export source 的 scale.data row names）。不把从 legacy H5 读到的 var index 再传进去冒充独立来源。

官方 exporter 源码：[easySCFr save.R](https://github.com/xleizi/easySCF/blob/main/r/R/save.R)。本地历史调用见 `20250706/Fig/SV/GSE198616/PBMC/revise.R:74`；历史 package 版本未保存，当前官方源码仅作语义/实现复核来源。
