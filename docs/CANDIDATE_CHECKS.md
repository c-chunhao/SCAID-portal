# Local source-candidate checks — 21 October 2026

These local checks apply to the v0.3.0-rc1 source candidate. The deployed
website was also accepted on 21 October 2026 within the production scope below.
Historical checks in other documents retain their original release/date scope.

- Backend: 80 selected Django tests passed again after the final source refresh, and one additional dynamic-thumbnail regression test passed (81 total checks) in an isolated test database using the existing web runtime, including markers, strict lazy indices/cache and new module scopes. Django system check reported no issues.
- Frontend: clean locked `npm ci --ignore-scripts`; 19 tests, ESLint and Vite production build passed under Node 24.13.0. The build reported existing Sass legacy-JS-API deprecation warnings. Dependencies/build outputs are excluded from the source package.
- Syntax: 119 Python sources parsed, 43 JSON files decoded, and the R H5AD export script parsed successfully.
- Portability: local path configuration loaded the 31-object aggregate example plan from an unrelated working directory; four scientific CLI help commands succeeded without biological inputs. This does not claim a full scientific re-execution or fresh Python/R dependency installation.
- Secret/path scan: no provider credential, private key, JSON credential payload or private absolute source path was found. One literal-secret pattern was manually classified as the fixed nonproduction Django test key in `backend/tisch_api/test_settings.py:2`; no matching key is reproduced in this report.
- Existing `LICENSE_STATUS.md` and `THIRD_PARTY_NOTICES.md` are preserved byte for byte. No new project licence or template permission is asserted.

The scientific inputs, production database, populated gene registry, generated
PNG cache, raw/processed matrices and individual metadata were excluded.
The source commit and production acceptance reports identify separate scopes;
local source checks alone do not establish deployed scientific validity.


Final source refresh also passed an independent source review: 79 objects equal
48 unchanged legacy objects plus 31 additional objects; every new module flag
is enabled in the final relative asset map. The three matrix-backed Gene
objects and all 59 relative prepared-family roots are present. The public JSON
and TSV download catalogues match across all 125 rows (124 main H5AD files and
one explicitly scoped optional raw-count subset), including names and SHA256.
Five new object records expose study context through the serializer and UI.
Dynamic original and thumbnail routes are covered for `public, no-cache`, ETag,
304 responses and changed-file 200 responses. The thumbnail check also verifies
that regenerated image content changes and dynamic responses bypass static
acceleration headers. The recorded production acceptance scope is summarized below.


## Deployed website acceptance — 21 October 2026

The two acceptance stages passed 518 API checks in total: 320 resolved checks
for catalogue/CellChat/markers/legacy/error/read-only flows, followed by 198
checks for final matrix downloads and matrix-backed Gene views. Stage two
made 374 HTTPS requests. All 125 files passed HEAD and HTTP Range checks;
four small Active Tattoo Uveitis files were downloaded completely and matched
their SHA256 values. These results do not claim complete network downloads or
full network SHA256 verification of the approximately 40 GB matrix library.
Five desktop/mobile scenarios and six actual Gene/KEGG interaction scenarios
also passed. The original two atlases remain unchanged; the PBMC pilot remains
descriptive and the optional raw-count subset keeps its stated scope.
