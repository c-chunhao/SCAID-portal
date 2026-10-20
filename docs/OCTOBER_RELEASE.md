# October 2026 source release

Version 0.2.0-rc1 records the running portal business logic and the current 79-entry
catalogue. The rebuilt frontend index matches the deployed index exactly. The
backend API source is copied from the working deployment, with environment-only
configuration and relative asset-root adaptation for local installation.

## Data layout

- `SCAID_FIGURE_ROOT`: rendered overview/marker/communication and legacy map tree.
- `SCAID_H5_DOWNLOAD_ROOTS`: reviewed RNA/AUCell/UCell/singscore files.
- `SCAID_PUBLISHED_FIGURES_BASE`: optional prepared new gene/score family tree.
- `SCAID_MEDIA_ROOT` / `SCAID_THUMBNAIL_ROOT`: registered images and preview cache.

Data files and patient metadata are not committed. The aggregate catalogue is
public resource metadata, not a production database dump. A source installation
starts with a separate SQLite demo; use the manifest importer for actual reviewed
assets. Existing production absolute paths must be migrated to the configured
local layout rather than copied blindly. No deployment or patient data changes
were needed to publish this source release.

## Analysis entry points

`science/README.md` indexes the individual-level SLE, source-annotation and
composition workflows. `science/analysis_vNext/README.md` indexes portable RNA
readers, QC, marker, score and integration safeguards. Supply source inputs under
the documented environment variables. The vNext name designates experimental
analysis safeguards; their presence does not mean they replaced the published
RPCA/Harmony atlas or validated a new disease mechanism.

## Licence and persistence

The existing `LICENSE_STATUS.md` and third-party notices remain unchanged. No
project-wide licence has been selected by this release. Data-reuse permissions
remain those of the original sources. The GitHub tag/release provides a versioned
reference and downloadable source archive, without claiming a separate archive DOI.
