# api/urls.py
from django.urls import path
from rest_framework.routers import DefaultRouter
from .openapi import openapi_schema
from .views import CellDataViewSet, PDFImageViewSet, AffectedSitesWithDiseasesViewSet, CellDataViewSetAll, H5FileViewSet, GeneKeggEnumViewSet

router = DefaultRouter()
router.register(r'cell-data', CellDataViewSet, basename='cell-data')
router.register(r'cell-data-all', CellDataViewSetAll, basename='cell-data-all')
router.register(r'pdf-images', PDFImageViewSet)
router.register(r'by-site', AffectedSitesWithDiseasesViewSet, basename='affected-sites')
router.register(r'h5-file', H5FileViewSet)
router.register(r'gene-kegg-enum', GeneKeggEnumViewSet)

urlpatterns = [path('schema/', openapi_schema, name='api-schema'), *router.urls]