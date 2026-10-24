import test from 'node:test'
import assert from 'node:assert/strict'
import { summarizeCatalogue } from '../src/utils/catalogueSummary.js'

test('mixed-cohort wording does not create another condition or independent source', () => {
  const summary = summarizeCatalogue([
    { abbreviation: 'AD', disease_label: 'Atopic dermatitis', dataset_label: 'GSE1', tissue: 'PBMC', cell_count: 12, sample_size: 2 },
    { abbreviation: 'AD', disease_label: 'Atopic dermatitis and healthy references', dataset_label: 'GSE1', tissue: 'Skin', cell_count: 20, sample_size: 3 },
    { abbreviation: 'SLE', disease_label: 'Systemic lupus erythematosus', dataset_label: 'GSE2', tissue: 'pbmc', cell_count: '18' },
  ])
  assert.deepEqual(summary, { objects: 3, sources: 2, cells: 50, samples: 5, conditions: 2, tissues: 2 })
})
