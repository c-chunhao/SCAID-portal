import os
import fitz
from django.core.management.base import BaseCommand
from django.core.files import File
from PIL import Image
import io
from django.utils import timezone
from tisch_api.models import PDFImage


class Command(BaseCommand):
    help = 'Process PDF and image files while preserving the directory structure.'

    def add_arguments(self, parser):
        parser.add_argument(
            'root_directory',
            type=str,
            help='Root directory containing PDF and image files.'
        )

    def handle(self, *args, **options):
        root_dir = options['root_directory']

        if not os.path.isdir(root_dir):
            self.stdout.write(self.style.ERROR(f'Directory does not exist: {root_dir}'))
            return

        # Supported image extensions.
        image_extensions = ['.jpg', '.jpeg', '.png', '.gif', '.bmp']

        for root, dirs, files in os.walk(root_dir):
            for file in files:
                file_lower = file.lower()
                file_path = os.path.join(root, file)

                if file_lower.endswith('.pdf'):
                    self.process_pdf(file_path, root, root_dir)
                elif any(file_lower.endswith(ext) for ext in image_extensions):
                    self.process_image(file_path, root, root_dir)

    def process_pdf(self, file_path, pdf_dir, root_dir):
        try:
            relative_path = os.path.relpath(pdf_dir, root_dir)
            doc = fitz.open(file_path)
            name = os.path.basename(file_path)

            for page_num in range(len(doc)):
                page = doc.load_page(page_num)
                pix = page.get_pixmap()

                img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
                img_io = io.BytesIO()
                img.save(img_io, format='JPEG', quality=85)
                img_io.seek(0)

                save_path = os.path.join('', relative_path)
                image_name = f"{os.path.splitext(name)[0]}_page_{page_num + 1}.jpg"

                pdf_image = PDFImage(
                    name=name,
                    path=file_path,
                    relative_path=relative_path,
                    page_number=page_num + 1,
                    is_pdf=True,
                    create_time=timezone.now(),
                    image_width=pix.width,
                    image_height=pix.height
                )

                pdf_image.image.save(
                    os.path.join(save_path, image_name),
                    File(img_io)
                )
                pdf_image.save()

                self.stdout.write(self.style.SUCCESS(
                    f'Extracted page {page_num + 1} from {file_path}'
                ))

            doc.close()

        except Exception as e:
            self.stdout.write(self.style.ERROR(
                f'Error processing PDF {file_path}: {str(e)}'
            ))

    def process_image(self, file_path, img_dir, root_dir):
        try:
            relative_path = os.path.relpath(img_dir, root_dir)
            img_name = os.path.basename(file_path)

            # Read the image dimensions.
            with Image.open(file_path) as img:
                width, height = img.size

            # Create a record with the original image path in both path and image.
            pdf_image = PDFImage(
                name=img_name,
                path=file_path,  # Original path.
                relative_path=relative_path,
                page_number=1,
                is_pdf=False,
                create_time=timezone.now(),
                image_width=width,
                image_height=height,
                image=file_path  # Also store the original path in the image field.
            )
            pdf_image.save()

            self.stdout.write(self.style.SUCCESS(
                f'Recorded metadata for image {file_path} (file was not copied)'
            ))

        except Exception as e:
            self.stdout.write(self.style.ERROR(
                f'Error processing image {file_path}: {str(e)}'
            ))