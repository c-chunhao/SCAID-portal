from rest_framework import serializers
from .models import CellData, PDFImage, AffectedSite, Disease, DiseaseSiteRelation, H5File, GeneKeggEnum
from collections import defaultdict
from .asset_paths import effective_abbreviation, manifest_entry, figure_manifest_entry
from .thumbnails import wants_preview
from rest_framework import serializers

class PDFImageSerializer(serializers.ModelSerializer):
    # Customize the image field to return a relative path.
    image = serializers.SerializerMethodField()
    # JPEG preview endpoint for the ~1 MB PNG maps; None for the small rasterized PDF pages.
    thumb_url = serializers.SerializerMethodField()
    condition_abbreviation = serializers.SerializerMethodField()
    disease_label = serializers.SerializerMethodField()
    dataset_label = serializers.SerializerMethodField()
    tissue_label = serializers.SerializerMethodField()

    def get_condition_abbreviation(self, obj):
        entry = figure_manifest_entry(obj)
        return (entry.get('canonical_abbreviation') or entry['abbreviation']) if entry else None

    def get_disease_label(self, obj):
        entry = figure_manifest_entry(obj)
        return (entry.get('display_name') or entry['disease_name']) if entry else None

    def get_dataset_label(self, obj):
        entry = figure_manifest_entry(obj)
        return entry['dataset_id'] if entry else None

    def get_tissue_label(self, obj):
        entry = figure_manifest_entry(obj)
        return entry['tissue'] if entry else None

    class Meta:
        model = PDFImage
        fields = ['id', 'name', 'relative_path', 'page_number', 'image', 'image_width', 'image_height', 'is_pdf', 'thumb_url', 'condition_abbreviation', 'disease_label', 'dataset_label', 'tissue_label']

    def get_thumb_url(self, obj):
        return f'/api/pdf-images/{obj.pk}/thumb/' if wants_preview(obj) else None

    def get_image(self, obj):
        return f'/api/pdf-images/{obj.pk}/original/' if obj.image else None


class DiseaseSerializer(serializers.ModelSerializer):
    class Meta:
        model = Disease
        fields = ['id', 'name_en', 'abbreviation', 'features', 'created_at', 'updated_at']
        extra_kwargs = {
            'created_at': {'format': '%Y-%m-%d %H:%M:%S'},
            'updated_at': {'format': '%Y-%m-%d %H:%M:%S'}
        }

# The catalogue table stores a few tissues in Chinese. Filters keep the stored
# value; Unicode escapes preserve exact keys. English labels are for display
# and are also used by the portal.
TISSUE_LABELS = {
    '\u76ae\u80a4': 'Skin',
    '\u80ba': 'Lung',
    '\u5916\u5468\u8840': 'PBMC',
    '\u6c34\u6ce1': 'Blister fluid',
    '\u773c\u524d\u623f\u6db2\uff08\u623f\u6c34\uff09': 'Aqueous humour',
}


def display_tissue(value):
    return TISSUE_LABELS.get((value or '').strip(), value)


def manifest_disease_label(obj):
    entry = manifest_entry(obj)
    if entry:
        return entry.get('display_name') or entry['disease_name']
    return obj.disease_full_name or obj.disease_name or obj.abbreviation


def manifest_dataset_label(obj):
    entry = manifest_entry(obj)
    return entry['dataset_id'] if entry else obj.dataset_id


def manifest_disease_name(obj):
    entry = manifest_entry(obj)
    return entry.get('canonical_disease_name') if entry and entry.get('canonical_disease_name') else obj.disease_name


def manifest_disease_full_name(obj):
    entry = manifest_entry(obj)
    if entry and entry.get('canonical_disease_name'):
        return entry['canonical_disease_name']
    return obj.disease_full_name or (entry['disease_name'] if entry else None)


class CellDataSerializer(serializers.ModelSerializer):
    disease_label = serializers.SerializerMethodField()
    tissue_label = serializers.SerializerMethodField()
    dataset_label = serializers.SerializerMethodField()
    module_availability = serializers.SerializerMethodField()
    availability_notes = serializers.SerializerMethodField()
    study_context = serializers.SerializerMethodField()

    def get_study_context(self, obj):
        entry = manifest_entry(obj)
        return entry.get('study_context', '') if entry else ''

    def get_module_availability(self, obj):
        entry = manifest_entry(obj)
        # Legacy releases do not yet carry a module-level publication manifest.
        # Null retains their existing behaviour without claiming new validation.
        return entry.get('module_availability') if entry else None

    def get_availability_notes(self, obj):
        entry = manifest_entry(obj)
        return entry.get('availability_notes', []) if entry else []

    def get_disease_label(self, obj):
        return manifest_disease_label(obj)

    def get_tissue_label(self, obj):
        return display_tissue(obj.tissue)

    def get_dataset_label(self, obj):
        return manifest_dataset_label(obj)

    disease_name = serializers.SerializerMethodField()
    disease_full_name = serializers.SerializerMethodField()

    def get_disease_name(self, obj):
        return manifest_disease_name(obj)

    def get_disease_full_name(self, obj):
        return manifest_disease_full_name(obj)

    abbreviation = serializers.SerializerMethodField()
    asset_status = serializers.SerializerMethodField()
    asset_message = serializers.SerializerMethodField()

    def get_abbreviation(self, obj):
        return effective_abbreviation(obj)

    def get_asset_status(self, obj):
        return 'linked' if manifest_entry(obj) else 'unverified'

    def get_asset_message(self, obj):
        return '' if manifest_entry(obj) else 'Asset association is not yet verified for this catalog record.'

    class Meta:
        model = CellData
        fields = '__all__'

class AffectedSiteSerializer(serializers.ModelSerializer):
    diseases = serializers.SerializerMethodField()

    class Meta:
        model = AffectedSite
        fields = ['id', 'site_name', 'name_en', 'description', 'created_at', 'diseases']
        extra_kwargs = {
            'created_at': {'format': '%Y-%m-%d %H:%M:%S'}
        }

    def get_diseases(self, obj):
        # Relations are prefetched by the view (see AffectedSitesWithDiseasesViewSet),
        # so this does not issue one query per site.
        relations = obj.diseasesiterelation_set.all()
        diseases = [relation.disease for relation in relations]
        return DiseaseSerializer(diseases, many=True).data


from collections import defaultdict


class CellDataTreeSerializer(serializers.ModelSerializer):
    abbreviation = serializers.SerializerMethodField()
    disease_label = serializers.SerializerMethodField()
    disease_name = serializers.SerializerMethodField()
    disease_full_name = serializers.SerializerMethodField()

    def get_disease_name(self, obj):
        return manifest_disease_name(obj)

    def get_disease_full_name(self, obj):
        return manifest_disease_full_name(obj)

    def get_abbreviation(self, obj):
        return effective_abbreviation(obj)

    def get_disease_label(self, obj):
        # Group labels by the curated public condition code.
        groups = self.context.get('catalog_groups') or {}
        records = groups.get(effective_abbreviation(obj)) or [obj]
        labels = sorted({manifest_disease_label(record) for record in records if manifest_disease_label(record)})
        return '; '.join(labels) if labels else manifest_disease_label(obj)

    children = serializers.SerializerMethodField()

    class Meta:
        model = CellData
        fields = ['abbreviation', 'disease_label', 'disease_name', 'disease_full_name', 'dataset_source', 'children']

    def get_children(self, obj):
        # Assume obj represents a group of records with the same abbreviation.
        effective = effective_abbreviation(obj)
        catalog_groups = self.context.get('catalog_groups')
        if catalog_groups is None:
            queryset = [record for record in CellData.objects.all() if effective_abbreviation(record) == effective]
        else:
            queryset = catalog_groups.get(effective, [])
        return self.build_tree_structure(queryset)

    @staticmethod
    def build_tree_structure(queryset):
        # Group by dataset_id, tissue, and batch to avoid duplicate base records.
        grouped_data = defaultdict(list)
        for item in queryset:
            key = (item.dataset_id, item.tissue, item.batch)
            grouped_data[key].append(item)

        # Build the nested structure.
        datasets = defaultdict(lambda: defaultdict(list))
        dataset_labels = {}
        for (dataset_id, tissue, batch), items in grouped_data.items():
            # Use the first item because these fields should match within each group.
            item = items[0]
            dataset_labels.setdefault(dataset_id, manifest_dataset_label(item))
            datasets[dataset_id][tissue].append({
                'batch': batch,
                'sample_size': item.sample_size,
                'cell_count': item.cell_count,
                'remarks': item.remarks
            })

        result = []
        for dataset_id, tissues in datasets.items():
            tissue_nodes = []
            for tissue, batch_items in tissues.items():
                tissue_nodes.append({
                    'tissue': tissue,
                    'tissue_label': display_tissue(tissue),
                    'children': batch_items
                })

            result.append({
                'dataset_id': dataset_id,
                'dataset_label': dataset_labels[dataset_id],
                'children': tissue_nodes
            })

        return result


class H5FileSerializer(serializers.ModelSerializer):
    download_url = serializers.SerializerMethodField()
    method = serializers.SerializerMethodField()
    file_format = serializers.SerializerMethodField()
    has_raw_counts = serializers.SerializerMethodField()
    content_description = serializers.SerializerMethodField()

    def _file_metadata(self, obj):
        if not hasattr(obj, '_scaid_download_metadata'):
            from .asset_paths import load_manifest
            value = {}
            for entry in load_manifest():
                for method, path in entry.get('h5_paths', {}).items():
                    if path == obj.path:
                        value = entry.get('h5_file_metadata', {}).get(method, {})
                        break
                if value:
                    break
            obj._scaid_download_metadata = value
        return obj._scaid_download_metadata

    def get_has_raw_counts(self, obj):
        return self._file_metadata(obj).get('has_raw_counts')

    def get_content_description(self, obj):
        return self._file_metadata(obj).get('content_description', '')

    def get_file_format(self, obj):
        from pathlib import Path
        return 'H5AD' if Path(obj.path).suffix.lower() == '.h5ad' else 'H5'

    def get_method(self, obj):
        return obj.relative_path.strip("/\\").replace("\\", "/").split("/")[-1]

    def get_download_url(self, obj):
        return f"/api/h5-file/{obj.pk}/download/"

    create_time = serializers.DateTimeField(format='%Y-%m-%d %H:%M:%S')
    last_modified = serializers.DateTimeField(format='%Y-%m-%d %H:%M:%S')

    class Meta:
        model = H5File
        fields = ['id', 'name', 'relative_path', 'file_size', 'create_time', 'last_modified', 'download_url', 'method', 'file_format', 'has_raw_counts', 'content_description']


class GeneKeggEnumSerializer(serializers.ModelSerializer):
    tp_display = serializers.CharField(source='get_tp_display', read_only=True)

    class Meta:
        model = GeneKeggEnum
        fields = ['id', 'tp', 'tp_display', 'name']
