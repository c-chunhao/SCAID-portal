# Independent scientific environment

Use a separate Python 3.12 environment for the current atlas-diagnostic and
B-cell-state and source-pool packages. Do not reuse the web service environment. A C++17
compiler and the zlib development library are required by the B-cell streaming
example. Python must include the standard `venv`/`ensurepip` support.

From the repository root:

```sh
python3.12 -m venv .venv-science
. .venv-science/bin/activate
python -m pip install -r science/current_revision/requirements-tested-python312.txt
cd science/current_revision/bcell_state
python synthetic_example/run_synthetic.py
cd ../atlas_diagnostics
SCAID_DIAGNOSTIC_DIR=synthetic_example python replay_sampled_diagnostics.py SYNTHETIC
SCAID_DIAGNOSTIC_DIR=synthetic_example python verify_final_diagnostics.py SYNTHETIC
```

A clean environment without inherited system or user site packages was actually
installed from PyPI wheels on 20 October 2026. The B-state and atlas synthetic examples and the
independent atlas calculation passed. Resolved versions are in
`installation_validation.json`; `requirements-tested-python312.txt` records the
complete environment installed for that check. This exact lock is the tested
replay entry point. The package-specific bounded requirements are an alternative
compatibility declaration, not a claim that every combination has been tested.
They allow
NumPy 1.26.4 or 2.x; the fresh installation used NumPy 2.5.4. The code
uses neither NumPy-2-only functionality nor the checked removed NumPy-1 aliases.

The package-specific `tested_versions.json` files record the earlier environment
used to calculate the real scientific results. They are provenance records,
not recommendations to recreate NumPy 1.26 on Python 3.13. No real analysis was
rerun merely to validate this installation. Synthetic inputs are explicitly
synthetic; real replay additionally requires the documented source inputs.

For the existing broader `science/scripts` entry points, install their
additional packages documented in `science/README.md` (including matplotlib,
h5py, PyDESeq2 and the listed R packages). This lightweight current-revision
environment checks the three current-revision packages and is not a lock for every older
scientific script or the web runtime.

The additional source-pool example passed with the same lock (200 invented
people, 20 pools, 56 planned rows and 54 estimable models), including an
independent QR/CR1/BH calculation. Run it from the repository root:

```sh
cd science/current_revision/source_pool_sensitivity
python scripts/synthetic_pool_demo.py --output work/pure_synthetic_pool_example
```
