# Full-stream MatrixMarket audit

`g++ -O3 -std=c++17 audit_matrixmarket_gz.cpp -lz -o audit_matrixmarket_gz`

Use `audit_new_data_sources.py --native-mtx /absolute/path/to/audit_matrixmarket_gz` with the regular inventory, registry, asset and output arguments. Each coordinate/value is scanned; record count, coordinate bounds and integer/finite/nonnegative measurements are checked. Native code avoids materializing the matrix. Python separately validates feature/barcode companions. This changes performance, not audit scope. An integral matrix is still not proof of original count provenance, clinical identity or release eligibility. If the native tool is unavailable the complete bounded Python/numpy stream remains supported.

Validation covers integer/real headers, malformed triples, record count/bounds, fractional counts, nonfinite/negative values and gzip/plain text. MatrixMarket permits duplicate coordinates; they require explicit sparse canonicalization, not silent cell/feature relabelling.
