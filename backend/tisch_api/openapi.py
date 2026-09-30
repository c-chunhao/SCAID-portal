"""Hand-maintained OpenAPI 3.0 description of the public, read-only SCAID API.

Served at /api/schema/. Kept in code (no extra dependency) and covered by
test_openapi.py, which checks that every router path is described here.
"""
from django.http import JsonResponse
from django.views.decorators.http import require_GET

API_VERSION = "1.3.0"

_ENVELOPE = {
    "type": "object",
    "properties": {
        "code": {"type": "integer", "example": 200},
        "message": {"type": "string", "example": "success"},
        "data": {"type": "array", "items": {"type": "object"}},
    },
}

_PAGED_ENVELOPE = {
    "type": "object",
    "properties": {
        "links": {
            "type": "object",
            "properties": {
                "next": {"type": "string", "nullable": True},
                "previous": {"type": "string", "nullable": True},
            },
        },
        "count": {"type": "integer"},
        "results": _ENVELOPE,
    },
}

_CELL_DATA = {
    "type": "object",
    "properties": {
        "id": {"type": "integer"},
        "disease_label": {"type": "string", "description": "Curated condition label from the asset manifest."},
        "disease_full_name": {"type": "string", "nullable": True},
        "disease_name": {"type": "string", "nullable": True},
        "abbreviation": {"type": "string", "nullable": True, "description": "Curated public condition code; historical storage codes do not override PsA / BD classification."},
        "dataset_source": {"type": "string", "example": "GEO"},
        "dataset_id": {"type": "string", "example": "GSE121380"},
        "tissue": {"type": "string", "example": "colon"},
        "sample_size": {"type": "integer", "nullable": True, "description": "Number of sample labels, not verified donors."},
        "cell_count": {"type": "integer", "nullable": True, "description": "Cell records in the released object."},
        "batch": {"type": "string", "nullable": True, "description": "Ingestion batch (YYYYMMDD); not a biological batch."},
        "asset_status": {"type": "string", "enum": ["linked", "unverified"]},
        "asset_message": {"type": "string"},
        "remarks": {"type": "string", "nullable": True},
        "file_name": {"type": "string", "nullable": True},
    },
}

_FIGURE = {
    "type": "object",
    "properties": {
        "id": {"type": "integer"},
        "name": {"type": "string"},
        "relative_path": {"type": "string", "description": "Historical relative asset folder, not a diagnosis. Prefer curated label fields for captions."},
        "page_number": {"type": "integer"},
        "image": {"type": "string", "nullable": True, "example": "/api/pdf-images/123/original/"},
        "image_width": {"type": "string"},
        "image_height": {"type": "string"},
        "is_pdf": {"type": "string"},
        "thumb_url": {"type": "string", "nullable": True},
        "condition_abbreviation": {"type": "string", "nullable": True},
        "disease_label": {"type": "string", "nullable": True},
        "dataset_label": {"type": "string", "nullable": True},
        "tissue_label": {"type": "string", "nullable": True},
    },
}
_FIGURE_ENVELOPE = {**_ENVELOPE, "properties": {**_ENVELOPE["properties"], "data": {"type": "array", "items": _FIGURE}}}
_FIGURE_PAGED_ENVELOPE = {**_PAGED_ENVELOPE, "properties": {**_PAGED_ENVELOPE["properties"], "results": _FIGURE_ENVELOPE}}

_PAGE_PARAMS = [
    {"name": "page", "in": "query", "schema": {"type": "integer", "minimum": 1}},
    {"name": "page_size", "in": "query", "schema": {"type": "integer", "minimum": 1, "maximum": 100},
     "description": "Default 10, maximum 100."},
]


def _q(name, description, schema=None):
    return {"name": name, "in": "query", "description": description, "schema": schema or {"type": "string"}}


SCHEMA = {
    "openapi": "3.0.3",
    "info": {
        "title": "SCAID public API",
        "version": API_VERSION,
        "description": (
            "Read-only endpoints behind the Single-Cell AutoImmune Disease database "
            "(https://scaid01.com/). All views are precomputed; no analysis runs on request. "
            "No registration or token is required. Write methods return 405."
        ),
    },
    "servers": [{"url": "https://scaid01.com/api"}],
    "paths": {
        "/cell-data/": {
            "get": {
                "summary": "Catalogue tree grouped by disease abbreviation",
                "description": "One node per disease abbreviation with nested dataset_id, tissue and batch entries. Used to populate filter menus. Cached server-side for a few minutes.",
                "responses": {"200": {"description": "Tree", "content": {"application/json": {"schema": _ENVELOPE}}}},
            }
        },
        "/cell-data/{id}/": {
            "get": {
                "summary": "One catalogue record rendered as a tree node",
                "parameters": [{"name": "id", "in": "path", "required": True, "schema": {"type": "integer"}}],
                "responses": {"200": {"description": "Node"}, "404": {"description": "Unknown id"}},
            }
        },
        "/cell-data-all/": {
            "get": {
                "summary": "Flat list of disease–tissue data objects",
                "description": "Returns a plain JSON array (not paginated; the catalogue is small). Text filters are case-insensitive substring matches unless stated otherwise.",
                "parameters": [
                    _q("name", "Free-text search across disease names, abbreviation, dataset source and ID, tissue and file name; also resolves manifest aliases such as IBD."),
                    _q("abbreviation", "Exact curated condition code. PsA selects E-MTAB-8207; AS selects true AS. BD selects GSE198616; BD(Uveitis) is a separate condition. SV remains a legacy alias for systemic BD."),
                    _q("disease_name", "Substring match."),
                    _q("disease_full_name", "Substring match."),
                    _q("dataset_source", "Repository name, e.g. GEO, ArrayExpress, Zenodo."),
                    _q("dataset_id", "Repository identifier, e.g. GSE174188."),
                    _q("tissue", "Tissue label, e.g. PBMC."),
                    _q("file_name", "Substring match."),
                    _q("batch", "Ingestion batch (YYYYMMDD)."),
                    _q("sample_size", "Exact integer.", {"type": "integer"}),
                    _q("cell_count", "Exact integer.", {"type": "integer"}),
                ],
                "responses": {"200": {"description": "Array of catalogue records", "content": {"application/json": {"schema": {"type": "array", "items": _CELL_DATA}}}}, "400": {"description": "Invalid integer parameter"}},
            }
        },
        "/cell-data-all/{id}/": {
            "get": {
                "summary": "One catalogue record",
                "parameters": [{"name": "id", "in": "path", "required": True, "schema": {"type": "integer"}}],
                "responses": {"200": {"description": "Record", "content": {"application/json": {"schema": _CELL_DATA}}}, "404": {"description": "Unknown id"}},
            }
        },
        "/pdf-images/": {
            "get": {
                "summary": "Precomputed figure files (rasterized PDF pages and PNG maps)",
                "description": "Scope a query with dataset_record_id (preferred) or with abbreviation / dataset_id / tissue. Select a figure family with pdf_name, or a gene or pathway map with name plus GeneUmap or KEGGScore.",
                "parameters": [
                    _q("dataset_record_id", "id from /cell-data-all/. Restricts results to that object's asset folder.", {"type": "integer", "minimum": 1}),
                    _q("abbreviation", "Disease abbreviation (used when dataset_record_id is absent)."),
                    _q("dataset_id", "Repository identifier."),
                    _q("dataset_source", "Repository name; refines dataset_id."),
                    _q("tissue", "Tissue label."),
                    _q("batch", "Ingestion batch folder."),
                    _q("pdf_name", "Figure family: Overview, Main_Marker, MajorCellType_Go_Ht, MajorCellType_KEGG_Ht, LRcircle, LRcontribution, SignalingCellRole, SignalingChord, SignalingCircle, SignalingHt, SignalingRole, Single, or an exact figure stem."),
                    _q("name", "Gene symbol or KEGG pathway id (e.g. TRIM21, hsa00010); matched against the stored file basename."),
                    _q("GeneUmap", "Present (any value) to restrict to gene expression UMAP maps."),
                    _q("KEGGScore", "Present (any value) to restrict to pathway score maps."),
                    _q("KEGGUmap", "Alias of KEGGScore."),
                    _q("kegg_param", "Scoring method folder: AUCell, UCell or singscore."),
                    *_PAGE_PARAMS,
                ],
                "responses": {"200": {"description": "Paginated figure records; image is a stable same-origin /api/pdf-images/{id}/original/ path, never an internal filesystem path", "content": {"application/json": {"schema": _FIGURE_PAGED_ENVELOPE}}}, "400": {"description": "Invalid page or integer parameter"}, "404": {"description": "Unknown dataset_record_id"}},
            }
        },
        "/pdf-images/{id}/": {
            "get": {
                "summary": "One figure record",
                "parameters": [{"name": "id", "in": "path", "required": True, "schema": {"type": "integer"}}],
                "responses": {"200": {"description": "Record"}, "404": {"description": "Unknown id"}},
            }
        },
        "/pdf-images/{id}/original/": {
            "get": {
                "summary": "Original registered figure image",
                "description": "Streams the PNG or JPEG associated with the figure id. Only files inside configured figure roots are served. GET and HEAD are supported; internal directory paths are not included in figure JSON.",
                "parameters": [{"name": "id", "in": "path", "required": True, "schema": {"type": "integer"}}],
                "responses": {"200": {"description": "image/png or image/jpeg"}, "404": {"description": "Unknown id or unavailable figure"}},
            }
        },
        "/pdf-images/{id}/thumb/": {
            "get": {
                "summary": "JPEG preview of a large figure PNG",
                "description": "Returns a cached JPEG preview (longest side 900 px, roughly 80 KB) of gene, pathway-score and other PNG maps whose originals are about 1 MB. Records whose figure is already a small rasterized page (thumb_url null in the list response) return 404. The original stays available at its stable /api/pdf-images/{id}/original/ path.",
                "parameters": [{"name": "id", "in": "path", "required": True, "schema": {"type": "integer"}}],
                "responses": {"200": {"description": "image/jpeg, cacheable for 30 days"}, "404": {"description": "No preview for this record"}},
            }
        },
        "/by-site/": {
            "get": {
                "summary": "Affected sites with their diseases",
                "description": "Used for the disease browser on the home page. Cached server-side for a few minutes.",
                "responses": {"200": {"description": "Sites", "content": {"application/json": {"schema": _ENVELOPE}}}},
            }
        },
        "/h5-file/": {
            "get": {
                "summary": "Downloadable H5 files for one data object",
                "description": "Requires dataset_record_id or name; without either the list is empty. Each object has an RNA matrix and AUCell, UCell and singscore pathway-score files (field `method`).",
                "parameters": [
                    _q("dataset_record_id", "id from /cell-data-all/.", {"type": "integer", "minimum": 1}),
                    _q("name", "Exact H5 file name (case-insensitive)."),
                    _q("min_size", "Minimum file size in bytes.", {"type": "integer"}),
                    _q("max_size", "Maximum file size in bytes.", {"type": "integer"}),
                    _q("start_date", "Registered on or after (YYYY-MM-DD)."),
                    _q("end_date", "Registered on or before (YYYY-MM-DD)."),
                ],
                "responses": {"200": {"description": "Array; download_url is a same-origin path", "content": {"application/json": {"schema": {"type": "array", "items": {"type": "object", "properties": {"id": {"type": "integer"}, "name": {"type": "string"}, "method": {"type": "string", "enum": ["RNA", "AUCell", "UCell", "singscore"]}, "file_size": {"type": "integer"}, "download_url": {"type": "string", "example": "/api/h5-file/273/download/"}}}}}}}, "400": {"description": "Invalid size or date parameter"}},
            }
        },
        "/h5-file/{id}/": {
            "get": {
                "summary": "One H5 file record",
                "parameters": [{"name": "id", "in": "path", "required": True, "schema": {"type": "integer"}}],
                "responses": {"200": {"description": "Record"}, "404": {"description": "Unknown id"}},
            }
        },
        "/h5-file/{id}/download/": {
            "get": {
                "summary": "Download an H5 file",
                "description": "Streams the registered file as application/x-hdf5 with a Content-Disposition filename. Only files registered inside the configured download roots are served. HEAD is supported for size checks.",
                "parameters": [{"name": "id", "in": "path", "required": True, "schema": {"type": "integer"}}],
                "responses": {"200": {"description": "Binary HDF5 file"}, "404": {"description": "Unknown id or file unavailable"}},
            }
        },
        "/gene-kegg-enum/": {
            "get": {
                "summary": "Indexed gene symbols and KEGG pathway ids for search fields",
                "parameters": [
                    _q("tp", "gene or kegg."),
                    _q("tps", "Comma-separated list of types."),
                    _q("name", "Substring match."),
                    _q("names", "Comma-separated exact names."),
                    _q("name_startswith", "Prefix match (case-sensitive)."),
                    _q("name_endswith", "Suffix match (case-sensitive)."),
                    _q("exclude_tp", "Exclude a type."),
                    _q("ordering", "name, -name, tp or -tp."),
                    *_PAGE_PARAMS,
                ],
                "responses": {"200": {"description": "Paginated names", "content": {"application/json": {"schema": _PAGED_ENVELOPE}}}, "400": {"description": "Invalid page parameter"}},
            }
        },
        "/gene-kegg-enum/{id}/": {
            "get": {
                "summary": "One enumeration record",
                "parameters": [{"name": "id", "in": "path", "required": True, "schema": {"type": "integer"}}],
                "responses": {"200": {"description": "Record"}, "404": {"description": "Unknown id"}},
            }
        },
        "/schema/": {
            "get": {"summary": "This document", "responses": {"200": {"description": "OpenAPI 3.0 JSON"}}}
        },
    },
}


@require_GET
def openapi_schema(request):
    return JsonResponse(SCHEMA, json_dumps_params={"ensure_ascii": False})
