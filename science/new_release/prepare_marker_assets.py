"""Validate and copy released within-object marker tables without changing values."""
from pathlib import Path
import csv, hashlib, json, math, shutil

from runtime_paths import WORK, CODE_ROOT, EXPORT_ROOT, RAW_ROOT, ANNOTATION_ROOT, PLAN_PATH, FIGURE_ROOT, INDEX_ROOT, GENE_MANIFEST, BACKEND_ROOT, load_plan
SOURCE = ANNOTATION_ROOT
ROOT = WORK / "markers"
FIELDS = ["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "group1", "group2", "test_group_number", "test_group"]
plan = load_plan()
catalog = json.loads((ROOT / "catalog_snapshot.json").read_text())
if isinstance(catalog, dict):
    catalog = catalog.get("records", catalog.get("objects", catalog.get("data")))
assert isinstance(catalog, list)
manifest = {"schema_version": 1, "comparison": "cell_type_vs_remaining_cells", "columns": FIELDS, "objects": {}}
report = {"objects": [], "total_rows": 0, "source_values_preserved": True}
for obj in plan["objects"]:
    c = obj["catalog"]
    candidates = [x for x in catalog if x["abbreviation"] == c["abbreviation"] and x["dataset_id"] == c["dataset_id"] and x["tissue"] == c["tissue"]]
    assert len(candidates) == 1, (obj["object_id"], len(candidates))
    record = candidates[0]
    record_id = int(record.get("id", record.get("record_id")))
    src = ANNOTATION_ROOT / obj["object_id"] / "DEGs.MajorCellType.txt"
    groups = set()
    count = 0
    negatives = 0
    with src.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == FIELDS, src
        for row in reader:
            assert None not in row and all(row[k] is not None and row[k] != "" for k in FIELDS), (src, count)
            assert all("\x00" not in v for v in row.values())
            for key in ("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"):
                value = float(row[key])
                assert math.isfinite(value), (src, count, key)
                if key != "avg_log2FC":
                    assert 0 <= value <= 1, (src, count, key, value)
            assert row["group2"] == "others", (src, count)
            assert row["test_group_number"].isdigit() and int(row["test_group_number"]) >= 1
            assert row["group1"] in row["test_group"].split(";"), (src, count)
            assert len(row["test_group"].split(";")) == int(row["test_group_number"]), (src, count)
            groups.add(row["group1"])
            negatives += float(row["avg_log2FC"]) < 0
            count += 1
    assert count > 0
    dst = ROOT / "tables" / f"{record_id}.tsv"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    checksum = hashlib.sha256(src.read_bytes()).hexdigest()
    assert hashlib.sha256(dst.read_bytes()).hexdigest() == checksum
    entry = {"dataset_record_id": record_id, "object_id": obj["object_id"], "abbreviation": c["abbreviation"], "dataset_id": c["dataset_id"], "tissue": c["tissue"], "relative_path": f"tables/{record_id}.tsv", "sha256": checksum, "bytes": dst.stat().st_size, "row_count": count, "cell_types": sorted(groups), "comparison": "cell_type_vs_remaining_cells", "filename": f"{c['abbreviation']}_{c['dataset_id']}_{c['tissue']}_cell_type_markers.tsv"}
    assert str(record_id) not in manifest["objects"]
    manifest["objects"][str(record_id)] = entry
    report["objects"].append({**entry, "source_path": str(src), "negative_log2fc_rows": negatives})
    report["total_rows"] += count
assert len(manifest["objects"]) == 31
(ROOT / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
(ROOT / "validation.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"objects": len(manifest["objects"]), "rows": report["total_rows"], "bytes": sum(x["bytes"] for x in report["objects"])}))
