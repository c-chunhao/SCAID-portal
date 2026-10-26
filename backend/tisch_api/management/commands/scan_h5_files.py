import os
from django.core.management.base import BaseCommand
from django.utils import timezone
from tisch_api.models import H5File  # Model used to record H5 file metadata.


class Command(BaseCommand):
    help = 'Scan a directory for .h5 files and record their paths in the database.'

    def add_arguments(self, parser):
        parser.add_argument(
            'root_directory',
            type=str,
            help='Root directory containing .h5 files.'
        )

    def handle(self, *args, **options):
        root_dir = options['root_directory']

        if not os.path.isdir(root_dir):
            self.stdout.write(self.style.ERROR(f'Directory does not exist: {root_dir}'))
            return

        for root, dirs, files in os.walk(root_dir):
            for file in files:
                if file.lower().endswith('.h5'):
                    file_path = os.path.join(root, file)
                    self.process_h5_file(file_path, root, root_dir)

    def process_h5_file(self, file_path, file_dir, root_dir):
        try:
            relative_path = os.path.relpath(file_dir, root_dir)
            file_name = os.path.basename(file_path)

            # Create and save the database record.
            H5File.objects.create(
                name=file_name,
                path=file_path,
                relative_path=relative_path,
                create_time=timezone.now(),
                file_size=os.path.getsize(file_path)
            )

            self.stdout.write(self.style.SUCCESS(
                f'Recorded .h5 file: {file_path}'
            ))

        except Exception as e:
            self.stdout.write(self.style.ERROR(
                f'Error processing .h5 file {file_path}: {str(e)}'
            ))