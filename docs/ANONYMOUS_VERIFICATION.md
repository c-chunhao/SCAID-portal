# Anonymous workflow and source-installation checks

The dated live-site and 52-file production comparison belongs to the
production-source snapshot released as v0.2.0-rc1. v0.2.0-rc2 adds independent
Python 3.12 installations, a local API/proxy check and source-pool code. These
new checks and the older live-site checks have separate scopes.

The live HTTPS site/API were checked on 20 October 2026 without login cookies
or credentials. The compact [machine-readable record](ANONYMOUS_VERIFICATION_20261020.json)
records the tested scope. Catalogue totals were 79 objects, 4,590,932 cell
records, 1,095 entry-specific sample labels and 36 condition labels.

Selected real workflows cover catalogue/IBD retrieval, GSE135779 and GSE174188
overview, CD3D expression, JAK–STAT UCell, existing CellChat Count/ICAM views,
registered downloads, precomputed atlas views/markers and the Information guide.
All 31 additional entries were checked for their advertised module families:
28 expression/score entries and three overview-only entries, with no new
CellChat/H5 availability. Eight registered SLE H5 files were checked using HTTP
byte ranges (first eight bytes) and the HDF5 signature; this is not a full
scientific validation of their contents.

Eleven fresh-browser workflows passed with no observed page error or HTTP5xx.
The rebuilt frontend matched all 52 production files. In isolated local
environments, 56 Django tests, 19 frontend unit tests, lint/build and a separate
SQLite migrate/import synthetic fixture passed. Both current scientific code
packages have explicitly synthetic checks; these check installation and
numerics, without reproducing real biological results in the absence of inputs.

These checks do not claim that every gene/pathway image was inspected, and do
not establish independent annotation accuracy, donor identity, integration
validity for unrestricted cross-disease comparisons or complete historical
raw-data reproduction. No production registrations or data were changed.
