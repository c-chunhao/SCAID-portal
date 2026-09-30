from datetime import datetime, time, timedelta

from django.db.models import Q
from django.conf import settings
from django.http import FileResponse, Http404, HttpResponse
from django.utils.http import content_disposition_header
from pathlib import Path
from urllib.parse import quote
from django.shortcuts import get_object_or_404
from django.utils import timezone
from django.utils.decorators import method_decorator
from django.views.decorators.cache import cache_page
from rest_framework import viewsets
from rest_framework.decorators import action
from rest_framework.exceptions import ValidationError
from rest_framework.pagination import PageNumberPagination
from rest_framework.response import Response

from .asset_paths import (
    DISEASE_FOLDERS, dataset_names, path_sequence_condition, effective_abbreviation, matching_catalog_records, manifest_search_ids, canonical_query_ids,
    manifest_image_condition, manifest_h5_condition, manifest_folders, folder_condition, resolved_h5_path, tissue_names,
)
from .models import PDFImage, CellData, AffectedSite, GeneKeggEnum, H5File
from .thumbnails import ensure_thumbnail, thumbnail_relative, wants_preview, resolved_figure_path
from .serializers import (
    PDFImageSerializer, CellDataTreeSerializer, AffectedSiteSerializer,
    CellDataSerializer, GeneKeggEnumSerializer, H5FileSerializer,
)


def integer_param(params, name, minimum=0, maximum=9223372036854775807):
    raw = params.get(name)
    if raw is None:
        return None
    try:
        if not str(raw).strip().isdigit():
            raise ValueError
        value = int(raw)
        if not minimum <= value <= maximum:
            raise ValueError
        return value
    except (ValueError, TypeError):
        raise ValidationError({name: 'Enter a valid integer between %s and %s.' % (minimum, maximum)}) from None


def requested_dataset(params):
    record_id = integer_param(params, 'dataset_record_id', minimum=1)
    return get_object_or_404(CellData, pk=record_id) if record_id is not None else None


class PublicReadOnlyViewSet(viewsets.ReadOnlyModelViewSet):
    http_method_names = ['get', 'head', 'options']


class ValidatedPageNumberPagination(PageNumberPagination):
    page_size = 10
    page_size_query_param = 'page_size'
    max_page_size = 100

    def paginate_queryset(self, queryset, request, view=None):
        integer_param(request.query_params, 'page_size', minimum=1)
        if request.query_params.get('page') not in self.last_page_strings:
            integer_param(request.query_params, 'page', minimum=1)
        return super().paginate_queryset(queryset, request, view)


class CustomPageNumberPagination(ValidatedPageNumberPagination):
    def get_paginated_response(self, data):
        return Response({
            'links': {'next': self.get_next_link(), 'previous': self.get_previous_link()},
            'count': self.page.paginator.count,
            'results': data,
        })


class PDFImageViewSet(PublicReadOnlyViewSet):
    queryset = PDFImage.objects.all().order_by('id')
    serializer_class = PDFImageSerializer
    pagination_class = CustomPageNumberPagination

    def get_queryset(self):
        queryset = super().get_queryset()
        params = self.request.query_params
        record = requested_dataset(params)
        # Manifest-linked records are scoped through the indexed relative_path
        # column; regex matching on the image path stays as the fallback only.
        folders = manifest_folders(record) if record is not None else []
        conditions = Q()
        if record is not None:
            # Catalog batch describes ingestion history; asset paths can be from
            # the current consolidated release and must not use that old batch.
            conditions &= folder_condition(folders, prefix=True) if folders else manifest_image_condition(record)
        else:
            scopes = {key: params.get(key) for key in ('abbreviation', 'dataset_id', 'tissue')}
            candidates = matching_catalog_records(CellData.objects.all(), **scopes) if any(scopes.values()) else []
            if candidates:
                scoped = Q(pk__in=[])
                for candidate in candidates:
                    candidate_folders = manifest_folders(candidate)
                    scoped |= folder_condition(candidate_folders, prefix=True) if candidate_folders else manifest_image_condition(candidate)
                conditions &= scoped
            else:
                # A known condition plus a mismatched accession/tissue is an
                # empty scientific scope, not permission to fall back to a
                # legacy directory (AS/E-MTAB-8207 contains reclassified PsA).
                if params.get('abbreviation') and canonical_query_ids(CellData.objects.all(), params['abbreviation']) is not None:
                    return queryset.none()
                if params.get('abbreviation'):
                    abbreviation = params['abbreviation'].strip()
                    conditions &= path_sequence_condition('image', (abbreviation, *DISEASE_FOLDERS.get(abbreviation, ())))
                if params.get('dataset_id'):
                    conditions &= path_sequence_condition('image', dataset_names(params['dataset_id'], params.get('dataset_source')))
                if params.get('tissue'):
                    conditions &= path_sequence_condition('image', tissue_names(params['tissue']))
            if params.get('batch'):
                conditions &= path_sequence_condition('image', (params['batch'],))

        if params.get('kegg_param'):
            method = params['kegg_param'].strip()
            conditions &= folder_condition(folders, f'KEGGScore/{method}') if folders else path_sequence_condition('image', (method,))
        if params.get('name'):
            clean_name = params['name'].strip()
            if clean_name.lower().endswith(('.png', '.jpg')):
                clean_name = clean_name[:-4]
            # Ingested PNGs retain their basename in name. Exact equality also
            # allows a filename index instead of scanning every image path.
            conditions &= Q(name__in=[f'{clean_name}.png', f'Umap_{clean_name}.png', f'{clean_name}.jpg', f'Umap_{clean_name}.jpg'])

        pdf_name = params.get('pdf_name')
        if pdf_name and any(key in params for key in ('GeneUmap', 'KEGGScore', 'KEGGUmap')):
            clean_name = pdf_name.strip()
            if clean_name.lower().endswith(('.png', '.jpg')):
                clean_name = clean_name[:-4]
            conditions &= Q(name__in=[f'{clean_name}.png', f'Umap_{clean_name}.png', f'{clean_name}.jpg', f'Umap_{clean_name}.jpg'])
        elif pdf_name == 'Overview':
            stems = ('Umap_Sample', 'Umap_RNA_snn_res.1', 'Umap_MajorCellType', 'MajorCellType_CellPro')
            names = [stem + suffix for stem in stems for suffix in ('.pdf', '.png', '.jpg', '_page_1.jpg')]
            conditions &= Q(name__in=names)
            conditions &= folder_condition(folders, 'Dataset/Overview') if folders else path_sequence_condition('image', ('Overview',))
        elif pdf_name in {
            'LRcircle', 'LRcontribution', 'SignalingCellRole', 'SignalingChord',
            'SignalingCircle', 'SignalingHt', 'SignalingRole', 'Single',
        }:
            # CellChat folders have no filename filter, so the regex fallback scans
            # the whole table (3.4 s live); the indexed folder match takes milliseconds.
            conditions &= folder_condition(folders, f'CellChat/{pdf_name}') if folders else path_sequence_condition('image', (pdf_name,))
        elif pdf_name:
            # Support original image names and rendered PDF pages, with a true
            # basename boundary (e.g. Count must not match OtherCount).
            import re
            clean_name = pdf_name.strip()
            names = [prefix + clean_name + suffix for prefix in ('', 'Umap_') for suffix in ('.pdf', '.png', '.jpg', '_page_1.jpg')]
            conditions &= Q(name__in=names)
            pattern = r'(?:^|[/\\])(?:Umap_)?' + re.escape(clean_name) + r'(?:_page_[0-9]+)?\.(?:jpg|png)$'
            conditions &= Q(image__iregex=pattern)

        if 'GeneUmap' in params:
            conditions &= folder_condition(folders, 'GeneUmap') if folders else path_sequence_condition('image', ('GeneUmap',))
        if 'KEGGScore' in params or 'KEGGUmap' in params:
            conditions &= folder_condition(folders, 'KEGGScore', prefix=True) if folders else path_sequence_condition('image', ('KEGGScore',))
        return queryset.filter(conditions)

    def list(self, request, *args, **kwargs):
        queryset = self.filter_queryset(self.get_queryset())
        page = self.paginate_queryset(queryset)
        data = self.get_serializer(page if page is not None else queryset, many=True).data
        payload = {'code': 200, 'message': 'success', 'data': data}
        return self.get_paginated_response(payload) if page is not None else Response(payload)

    @action(detail=True, methods=['get', 'head'], url_path='original')
    def original(self, request, pk=None):
        """Original registered figure, without disclosing a server directory."""
        record = get_object_or_404(PDFImage, pk=pk)
        source = resolved_figure_path(record)
        if source is None or source.suffix.lower() not in {'.png', '.jpg', '.jpeg'}:
            raise Http404('Figure file unavailable.')
        import mimetypes
        response = FileResponse(source.open('rb'), content_type=mimetypes.guess_type(source.name)[0] or 'application/octet-stream')
        response['Content-Disposition'] = content_disposition_header(False, source.name)
        response['Cache-Control'] = 'public, max-age=86400'
        response['X-Content-Type-Options'] = 'nosniff'
        return response

    @action(detail=True, methods=['get', 'head'], url_path='thumb')
    def thumb(self, request, pk=None):
        """JPEG preview (max 900 px) of a large figure PNG; see thumbnails.py."""
        record = PDFImage.objects.filter(pk=pk).first()
        if record is None or not wants_preview(record):
            raise Http404('No preview for this figure.')
        target = ensure_thumbnail(record)
        if target is None:
            raise Http404('Figure file unavailable.')
        prefix = getattr(settings, 'THUMBNAIL_ACCEL_PREFIX', '')
        if prefix:
            response = HttpResponse(content_type='image/jpeg')
            response['X-Accel-Redirect'] = prefix.rstrip('/') + '/' + quote(thumbnail_relative(record.pk).as_posix(), safe='/')
        else:
            response = FileResponse(target.open('rb'), content_type='image/jpeg')
        response['Cache-Control'] = 'public, max-age=2592000'
        response['X-Content-Type-Options'] = 'nosniff'
        return response


CATALOG_CACHE_SECONDS = getattr(settings, 'CATALOG_CACHE_SECONDS', 0)


class CellDataViewSet(PublicReadOnlyViewSet):
    queryset = CellData.objects.all().order_by('id')
    serializer_class = CellDataTreeSerializer

    @method_decorator(cache_page(CATALOG_CACHE_SECONDS))
    def list(self, request, *args, **kwargs):
        representatives = {}
        catalog_groups = {}
        for record in self.get_queryset():
            # Null abbreviations must not merge unrelated diseases into one tree.
            key = effective_abbreviation(record)
            representatives.setdefault(key, record)
            catalog_groups.setdefault(key, []).append(record)
        context = {**self.get_serializer_context(), 'catalog_groups': catalog_groups}
        return Response({'code': 200, 'message': 'success', 'data': self.get_serializer(representatives.values(), many=True, context=context).data})


class CellDataViewSetAll(PublicReadOnlyViewSet):
    queryset = CellData.objects.all().order_by('id')
    serializer_class = CellDataSerializer
    pagination_class = None

    def get_queryset(self):
        queryset = super().get_queryset()
        params = self.request.query_params
        conditions = Q()
        searchable = ('disease_name', 'disease_full_name', 'abbreviation', 'dataset_source', 'dataset_id', 'tissue', 'file_name')
        if params.get('name'):
            search = Q()
            for field in searchable:
                search |= Q(**{field + '__icontains': params['name'].strip()})
            exact_condition_ids = canonical_query_ids(CellData.objects.all(), params['name'])
            if exact_condition_ids is not None:
                search = Q(pk__in=exact_condition_ids)
            else:
                search |= Q(pk__in=manifest_search_ids(CellData.objects.all(), params['name']))
            conditions &= search
        for field in (*searchable, 'batch'):
            if params.get(field):
                if field == 'abbreviation':
                    records = matching_catalog_records(CellData.objects.all(), abbreviation=params[field])
                    conditions &= Q(pk__in=[record.pk for record in records])
                else:
                    conditions &= Q(**{field + '__icontains': params[field]})
        for field in ('sample_size', 'cell_count'):
            value = integer_param(params, field, maximum=2147483647)
            if value is not None:
                conditions &= Q(**{field: value})
        return queryset.filter(conditions)


class AffectedSitesWithDiseasesViewSet(viewsets.ViewSet):
    http_method_names = ['get', 'head', 'options']

    @method_decorator(cache_page(CATALOG_CACHE_SECONDS))
    def list(self, request):
        sites = AffectedSite.objects.order_by('id').prefetch_related('diseasesiterelation_set__disease')
        data = AffectedSiteSerializer(sites, many=True).data
        return Response({'code': 200, 'message': 'success', 'data': data})


class H5FileViewSet(PublicReadOnlyViewSet):
    queryset = H5File.objects.all().order_by('id')
    serializer_class = H5FileSerializer
    pagination_class = None

    def get_queryset(self):
        queryset = super().get_queryset()
        if self.action in {'retrieve', 'download'}:
            return queryset
        params = self.request.query_params
        record = requested_dataset(params)
        name = (params.get('name') or '').strip()
        if record is not None:
            # Match a complete catalog-derived filename, including disease and
            # tissue. Never use nullable file_name as a whole-library fallback.
            queryset = queryset.filter(manifest_h5_condition(record))
        elif name and name.lower() not in {'null', 'undefined'}:
            queryset = queryset.filter(name__iexact=name)
        else:
            queryset = queryset.none()

        minimum = integer_param(params, 'min_size')
        maximum = integer_param(params, 'max_size')
        if minimum is not None and maximum is not None and minimum > maximum:
            raise ValidationError({'max_size': 'Must be greater than or equal to min_size.'})
        if minimum is not None:
            queryset = queryset.filter(file_size__gte=minimum)
        if maximum is not None:
            queryset = queryset.filter(file_size__lte=maximum)
        dates = {}
        for key in ('start_date', 'end_date'):
            if key in params:
                try:
                    dates[key] = datetime.strptime(params[key], '%Y-%m-%d').date()
                except (ValueError, TypeError):
                    raise ValidationError({key: 'Use YYYY-MM-DD.'}) from None
        if len(dates) == 2 and dates['start_date'] > dates['end_date']:
            raise ValidationError({'end_date': 'Must be on or after start_date.'})
        if 'start_date' in dates:
            queryset = queryset.filter(create_time__gte=timezone.make_aware(datetime.combine(dates['start_date'], time.min)))
        if 'end_date' in dates:
            try:
                end = dates['end_date'] + timedelta(days=1)
            except OverflowError:
                raise ValidationError({'end_date': 'Date is out of range.'}) from None
            queryset = queryset.filter(create_time__lt=timezone.make_aware(datetime.combine(end, time.min)))
        return queryset

    @action(detail=True, methods=['get', 'head'], url_path='download')
    def download(self, request, pk=None):
        record = self.get_object()
        path = resolved_h5_path(record)
        prefix = getattr(settings, 'H5_DOWNLOAD_ACCEL_PREFIX', '')
        if prefix:
            roots = [Path(root).resolve() for root in settings.H5_DOWNLOAD_ROOTS]
            root = next(root for root in roots if path.is_relative_to(root))
            response = HttpResponse(content_type='application/x-hdf5')
            response['X-Accel-Redirect'] = prefix.rstrip('/') + '/' + quote(path.relative_to(root).as_posix(), safe='/')
            response['Content-Disposition'] = content_disposition_header(True, path.name)
            response['Content-Length'] = path.stat().st_size
            response['X-Content-Type-Options'] = 'nosniff'
            return response
        try:
            stream = path.open('rb')
        except OSError:
            raise Http404('File unavailable.') from None
        response = FileResponse(stream, as_attachment=True, filename=path.name, content_type='application/x-hdf5')
        response['X-Content-Type-Options'] = 'nosniff'
        return response


class GeneKeggEnumPagination(ValidatedPageNumberPagination):
    pass


class GeneKeggEnumViewSet(PublicReadOnlyViewSet):
    queryset = GeneKeggEnum.objects.all().order_by('name', 'id')
    serializer_class = GeneKeggEnumSerializer
    pagination_class = GeneKeggEnumPagination

    def get_queryset(self):
        queryset = super().get_queryset()
        params = self.request.query_params
        conditions = Q()
        if params.get('tp'):
            conditions &= Q(tp__iexact=params['tp'].strip().lower())
        if params.get('name'):
            conditions &= Q(name__icontains=params['name'])
        if params.get('tps'):
            types = Q(pk__in=[])
            for tp in params['tps'].split(','):
                types |= Q(tp__iexact=tp.strip().lower())
            conditions &= types
        if params.get('names'):
            conditions &= Q(name__in=[name.strip() for name in params['names'].split(',')])
        if params.get('exclude_tp'):
            conditions &= ~Q(tp__iexact=params['exclude_tp'].strip().lower())
        if params.get('name_startswith'):
            conditions &= Q(name__startswith=params['name_startswith'])
        if params.get('name_endswith'):
            conditions &= Q(name__endswith=params['name_endswith'])
        queryset = queryset.filter(conditions)
        ordering = params.get('ordering')
        if ordering in ['name', '-name', 'tp', '-tp']:
            queryset = queryset.order_by(ordering, 'id')
        return queryset

    def list(self, request, *args, **kwargs):
        queryset = self.filter_queryset(self.get_queryset())
        page = self.paginate_queryset(queryset)
        data = self.get_serializer(page if page is not None else queryset, many=True).data
        payload = {'code': 200, 'message': 'success', 'data': data}
        return self.get_paginated_response(payload) if page is not None else Response(payload)
