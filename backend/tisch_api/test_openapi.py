"""The hand-written OpenAPI document must describe every routed endpoint."""
import re

from django.test import TestCase
from django.urls import URLPattern, URLResolver, get_resolver
from rest_framework.test import APIClient

from .models import AffectedSite, Disease, DiseaseSiteRelation
from .openapi import SCHEMA


def routed_paths():
    found = set()

    def walk(patterns, prefix=''):
        for entry in patterns:
            if isinstance(entry, URLResolver):
                walk(entry.url_patterns, prefix + str(entry.pattern))
            elif isinstance(entry, URLPattern):
                route = prefix + str(entry.pattern)
                # Skip DRF's optional ".json" format-suffix variants and the API root.
                if route.startswith('api/') and 'format' not in route:
                    found.add('/' + route[len('api/'):])
    walk(get_resolver().url_patterns)
    return {p for p in found if p not in {'/', ''}}


class OpenApiSchemaTests(TestCase):
    def test_schema_served_as_json(self):
        response = APIClient().get('/api/schema/')
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response['Content-Type'].split(';')[0], 'application/json')
        body = response.json()
        self.assertEqual(body['openapi'], '3.0.3')
        self.assertIn('/cell-data-all/', body['paths'])

    def test_every_routed_path_is_documented(self):
        documented = set(SCHEMA['paths'])
        for route in routed_paths():
            normalized = re.sub(r'\(\?P<pk>[^)]*\)', '{id}', route).replace('<pk>', '{id}')
            normalized = '/' + normalized.lstrip('/').lstrip('^').rstrip('$')
            self.assertIn(normalized, documented, route)

    def test_schema_rejects_writes(self):
        self.assertEqual(APIClient().post('/api/schema/', {}).status_code, 405)


class AffectedSiteQueryTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        for i in range(5):
            site = AffectedSite.objects.create(site_name=f'部位{i}', name_en=f'site{i}')
            for j in range(3):
                disease = Disease.objects.create(name_en=f'disease{i}{j}', abbreviation=f'D{i}{j}')
                DiseaseSiteRelation.objects.create(site=site, disease=disease)

    def test_by_site_uses_bounded_queries(self):
        # One query for sites, one for relations, one for diseases; not one per site.
        with self.assertNumQueries(3):
            response = APIClient().get('/api/by-site/', {'format': 'json'})
        self.assertEqual(response.status_code, 200)
        data = response.json()['data']
        self.assertEqual(len(data), 5)
        self.assertEqual(sorted(d['abbreviation'] for d in data[0]['diseases']), ['D00', 'D01', 'D02'])
