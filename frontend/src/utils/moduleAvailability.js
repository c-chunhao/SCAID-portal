// Only an explicit false hides a module. Existing catalogue entries without
// publication metadata keep their established behaviour.
export const hasPublishedModule = (record, module) => record?.module_availability?.[module] !== false

const analysisModules = {
  Overview: 'overview',
  Main_Marker: 'marker',
  MajorCellType_Go_Ht: 'go',
  MajorCellType_KEGG_Ht: 'kegg',
  CellChatOverview: 'cellchat_overview',
  LRcircle: 'cellchat',
  LRcontribution: 'cellchat',
  SignalingCellRole: 'cellchat',
  SignalingChord: 'cellchat',
  SignalingCircle: 'cellchat',
  SignalingHt: 'cellchat',
  SignalingRole: 'cellchat',
  Single: 'cellchat'
}

export const hasPublishedAnalysis = (record, analysis) => Object.hasOwn(analysisModules, analysis)
  && hasPublishedModule(record, analysisModules[analysis])
  && (analysis !== 'CellChatOverview' || record?.module_availability?.cellchat_overview === true)

export const publicationNotes = record => Array.isArray(record?.availability_notes)
  ? record.availability_notes.filter(note => typeof note === 'string' && note.trim())
  : []
