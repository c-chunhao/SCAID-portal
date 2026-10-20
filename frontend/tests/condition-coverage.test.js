import test from 'node:test'
import assert from 'node:assert/strict'
import { buildConditionCoverage } from '../src/utils/conditionCoverage.js'

const conditions = groups => groups.flatMap(group => group.conditions)

test('all current raw condition identities appear exactly once across the ten navigation groups', () => {
  // Independent raw identities from the published 16 October catalogue, not
  // aliases or condition names reconstructed from the implementation groups.
  const abbreviations = ['SLE', 'SSs', 'APS', 'BD', 'pSS', 'DM', 'RA', 'PsA', 'AS', 'JIA', 'Gout', 'PSO', 'AD', 'AP', 'BP', 'DLE', 'LS', 'PN', 'DAHND', 'CD', 'UC', 'CI', 'CIC', 'AIH', 'ATU', 'ITU', 'BD(Uveitis)', 'VKHD(Uveitis)', 'MS', 'MG', 'AIHA', 'ILD', 'SAR', 'T1DM', 'GERD', 'Lyme']
  assert.equal(abbreviations.length, 36)
  assert.equal(new Set(abbreviations).size, 36)
  const rows = abbreviations.flatMap(abbreviation => [
    { abbreviation, disease_label: `Published ${abbreviation}`, dataset_id: 'first' },
    { abbreviation, disease_label: `Published ${abbreviation}`, dataset_id: 'second' },
  ])
  const result = buildConditionCoverage(rows)
  assert.equal(result.length, 10)
  assert.deepEqual(conditions(result).map(condition => condition.abbreviation).sort(), [...abbreviations].sort())
  assert.ok(conditions(result).every(condition => condition.label === `Published ${condition.abbreviation}`))
})

test('distinct source identities keep their literal punctuation even when labels overlap', () => {
  const rows = [
    { abbreviation: 'BD', disease_label: 'Behçet disease' },
    { abbreviation: 'BD(Uveitis)', disease_label: 'Behçet disease' },
    { abbreviation: 'VKHD(Uveitis)', disease_label: 'Vogt–Koyanagi–Harada disease' },
  ]
  const result = conditions(buildConditionCoverage(rows))
  assert.deepEqual(new Set(result.map(condition => condition.abbreviation)), new Set(['BD', 'BD(Uveitis)', 'VKHD(Uveitis)']))
  for (const condition of result) {
    assert.equal(new URLSearchParams({ id: condition.abbreviation }).get('id'), condition.abbreviation)
  }
  assert.equal(result.filter(condition => condition.label === 'Behçet disease').length, 2)
})

test('future conditions remain visible and an incomplete catalogue never creates unpublished links', () => {
  const rows = [{ abbreviation: 'SLE', disease_label: '' }, { abbreviation: 'Future(condition)', disease_label: 'Future published condition' }, { abbreviation: 'SLE', disease_label: 'Systemic lupus erythematosus' }]
  const original = structuredClone(rows)
  const groups = buildConditionCoverage(rows)
  assert.deepEqual(conditions(groups), [{ abbreviation: 'SLE', label: 'Systemic lupus erythematosus' }, { abbreviation: 'Future(condition)', label: 'Future published condition' }])
  assert.equal(groups.at(-1).title, 'Other catalogued conditions')
  assert.deepEqual(rows, original)
  assert.deepEqual(buildConditionCoverage(null), [])
  assert.deepEqual(buildConditionCoverage([null, {}, { abbreviation: '' }]), [])
})
