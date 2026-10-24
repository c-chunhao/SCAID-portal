// Cohort descriptions can change while the underlying condition code does not.
// Sample/library labels are summed as catalogue values, not donor identities.
export function summarizeCatalogue(rows) {
  const unique = values => new Set(values.filter(Boolean)).size
  const total = field => rows.reduce((sum, row) => sum + (Number(row[field]) || 0), 0)
  return {
    objects: rows.length,
    sources: unique(rows.map(row => row.dataset_label || row.dataset_id)),
    cells: total('cell_count'),
    samples: total('sample_size'),
    conditions: unique(rows.map(row => String(row.abbreviation || '').trim().toLowerCase())),
    tissues: unique(rows.map(row => String(row.tissue_label || row.tissue || '').trim().toLowerCase())),
  }
}
