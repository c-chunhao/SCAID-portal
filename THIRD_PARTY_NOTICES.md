# Third-party provenance and notices

The frontend's inspected local Git origin is
[boyane/vue-rabbit](https://gitee.com/boyane/vue-rabbit.git), and the local template
commit is `cb57e6a`. The local template tree has no root LICENSE or NOTICE file.
No upstream licence grant is inferred from the public URL. The SCAID snapshot
contains later project-specific work and excludes the template's Git history.
The package name is now `scaid-web`; that rename does not remove upstream rights.

`frontend/src/styles/bootstrap.scss` retains the Bootstrap **v4.4.1** copyright
header: Bootstrap Authors and Twitter, Inc. Its original MIT licence is included
in `third_party/BOOTSTRAP_LICENSE.txt`, obtained during the earlier source audit
from the official upstream v4.4.1 tag. This applies to Bootstrap, not to the whole
SCAID project.

Source copyright comments and applicable embedded notices are retained. npm
and Python dependency installations are excluded from the source archive;
`frontend/package-lock.json` and backend requirement files record their versions.
Those dependencies retain their own distribution licences and notices.

The current SCAID wordmark, tissue schematic and aggregate atlas figures are
included as frontend assets. Their authorship and reuse rights are not assigned
by this notice. Per-cell matrices, donor/clinical metadata and the full external
research data library are excluded. No project-wide software/data licence or
third-party permission has been invented by preparing this snapshot.

The 24 October update adds project analysis scripts and changes only identified
custom application files. It does not redistribute another template checkout or
assert a new licence for the historical template-derived frontend. The
original-to-public source map documents path-only scientific-source redaction.
External KEGG, CellTypist models and CellChat DB weights are not bundled.

The 24 October interface update bundles the IBM Plex Sans Latin WOFF2, copyright © 2017 IBM Corp., under SIL Open Font License 1.1. The complete licence is in `frontend/public/fonts/OFL-IBM-Plex-Sans.txt`; `frontend/public/fonts/SOURCE.json` records the public source and SHA256. The font is served locally. Its licence applies to the font, without changing the project or historical template licence status.
