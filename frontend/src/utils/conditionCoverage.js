// Navigation groups for the 16 October 2026 catalogue snapshot. They are not
// exclusive biological classifications or a dataset's sampled-tissue labels.
export const conditionCoverageGroups = [
  { key: 'systemic', title: 'Systemic & connective', abbreviations: ['SLE', 'SSs', 'APS', 'BD', 'pSS', 'DM'] },
  { key: 'joints', title: 'Joints & spine', abbreviations: ['RA', 'PsA', 'AS', 'JIA', 'Gout'] },
  { key: 'skin', title: 'Skin & soft tissue', abbreviations: ['PSO', 'AD', 'AP', 'BP', 'DLE', 'LS', 'PN', 'DAHND'] },
  { key: 'gut', title: 'Gut & liver', abbreviations: ['CD', 'UC', 'CI', 'CIC', 'AIH'] },
  { key: 'eye', title: 'Eye inflammation', abbreviations: ['ATU', 'ITU', 'BD(Uveitis)', 'VKHD(Uveitis)'] },
  { key: 'neuro', title: 'Neurologic & neuromuscular', abbreviations: ['MS', 'MG'] },
  { key: 'blood', title: 'Blood & marrow', abbreviations: ['AIHA'] },
  { key: 'lung', title: 'Lung & granulomatous', abbreviations: ['ILD', 'SAR'] },
  { key: 'endocrine', title: 'Endocrine', abbreviations: ['T1DM'] },
  { key: 'other', title: 'Other included conditions', abbreviations: ['GERD', 'Lyme'] },
]

// Start from the live published catalogue, rather than create links for every
// condition in this static grouping plan. Preserve raw API identities in links.
export function buildConditionCoverage(records) {
  const catalogue = new Map()
  for (const record of Array.isArray(records) ? records : []) {
    if (typeof record?.abbreviation !== 'string' || !record.abbreviation.trim()) continue
    const abbreviation = record.abbreviation
    const label = typeof record.disease_label === 'string' ? record.disease_label.trim() : ''
    const existing = catalogue.get(abbreviation)
    if (!existing || (label && existing.label === abbreviation)) {
      catalogue.set(abbreviation, { abbreviation, label: label || abbreviation })
    }
  }

  const assigned = new Set()
  const groups = conditionCoverageGroups.map(group => ({
    key: group.key,
    title: group.title,
    conditions: group.abbreviations.flatMap(abbreviation => {
      if (!catalogue.has(abbreviation) || assigned.has(abbreviation)) return []
      assigned.add(abbreviation)
      return [catalogue.get(abbreviation)]
    }),
  })).filter(group => group.conditions.length)

  const remaining = [...catalogue.values()]
    .filter(condition => !assigned.has(condition.abbreviation))
    .sort((first, second) => first.label.localeCompare(second.label, 'en', { sensitivity: 'base' }))
  if (remaining.length) groups.push({ key: 'other-catalogued', title: 'Other catalogued conditions', conditions: remaining })
  return groups
}
