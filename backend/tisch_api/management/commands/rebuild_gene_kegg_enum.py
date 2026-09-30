from collections import Counter

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from tisch_api.catalog import collect_catalog
from tisch_api.models import GeneKeggEnum


class Command(BaseCommand):
    help = 'Add missing gene/KEGG search terms from published figures; existing terms are preserved.'

    def add_arguments(self, parser):
        parser.add_argument('--figure-root', default=str(settings.SCAID_FIGURE_ROOT))
        parser.add_argument('--apply', action='store_true', help='Persist new terms (default is a read-only preview).')

    def handle(self, *args, **options):
        try:
            terms = collect_catalog(options['figure_root'])
        except ValueError as exc:
            raise CommandError(str(exc)) from exc
        if not terms:
            raise CommandError('No published gene/KEGG figures found; refusing to populate an empty catalog.')
        if any(len(name) > 50 for _, name in terms):
            raise CommandError('At least one figure name exceeds the database field length; no changes made.')
        counts = Counter(tp for tp, _ in terms)
        self.stdout.write(f"Found {counts['gene']} genes and {counts['kegg']} KEGG pathways.")
        if not options['apply']:
            self.stdout.write('Preview only. Use --apply after the migration has been applied.')
            return
        with transaction.atomic():
            existing = set(GeneKeggEnum.objects.values_list('tp', 'name'))
            missing = sorted(terms - existing)
            GeneKeggEnum.objects.bulk_create(
                [GeneKeggEnum(tp=tp, name=name) for tp, name in missing],
                batch_size=1000,
            )
        self.stdout.write(self.style.SUCCESS(f'Added {len(missing)} terms; existing terms preserved.'))
