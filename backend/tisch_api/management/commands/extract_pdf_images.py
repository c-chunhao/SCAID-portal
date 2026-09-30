import os
import fitz
from django.core.management.base import BaseCommand
from django.core.files import File
from PIL import Image
import io
from django.utils import timezone
from tisch_api.models import PDFImage


class Command(BaseCommand):
    help = '处理PDF和图片文件，保留目录结构'

    def add_arguments(self, parser):
        parser.add_argument(
            'root_directory',
            type=str,
            help='包含PDF和图片文件的根目录'
        )

    def handle(self, *args, **options):
        root_dir = options['root_directory']

        if not os.path.isdir(root_dir):
            self.stdout.write(self.style.ERROR(f'目录不存在: {root_dir}'))
            return

        # 支持的图片扩展名
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
                    f'从 {file_path} 提取了第 {page_num + 1} 页'
                ))

            doc.close()

        except Exception as e:
            self.stdout.write(self.style.ERROR(
                f'处理PDF {file_path} 时出错: {str(e)}'
            ))

    def process_image(self, file_path, img_dir, root_dir):
        try:
            relative_path = os.path.relpath(img_dir, root_dir)
            img_name = os.path.basename(file_path)

            # 获取图片尺寸
            with Image.open(file_path) as img:
                width, height = img.size

            # 创建记录，将原始图片路径同时保存到path和image字段
            pdf_image = PDFImage(
                name=img_name,
                path=file_path,  # 原始路径
                relative_path=relative_path,
                page_number=1,
                is_pdf=False,
                create_time=timezone.now(),
                image_width=width,
                image_height=height,
                image=file_path  # 也将原始路径保存到image字段
            )
            pdf_image.save()

            self.stdout.write(self.style.SUCCESS(
                f'记录了图片文件 {file_path} 的信息(未复制文件)'
            ))

        except Exception as e:
            self.stdout.write(self.style.ERROR(
                f'处理图片 {file_path} 时出错: {str(e)}'
            ))