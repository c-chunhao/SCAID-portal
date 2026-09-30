import os
from django.core.management.base import BaseCommand
from django.utils import timezone
from tisch_api.models import H5File  # 假设您已经或将要创建这个模型


class Command(BaseCommand):
    help = '扫描目录中的.h5文件并将其路径记录到数据库中'

    def add_arguments(self, parser):
        parser.add_argument(
            'root_directory',
            type=str,
            help='包含.h5文件的根目录路径'
        )

    def handle(self, *args, **options):
        root_dir = options['root_directory']

        if not os.path.isdir(root_dir):
            self.stdout.write(self.style.ERROR(f'目录不存在: {root_dir}'))
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

            # 创建记录并保存到数据库
            H5File.objects.create(
                name=file_name,
                path=file_path,
                relative_path=relative_path,
                create_time=timezone.now(),
                file_size=os.path.getsize(file_path)
            )

            self.stdout.write(self.style.SUCCESS(
                f'已记录.h5文件: {file_path}'
            ))

        except Exception as e:
            self.stdout.write(self.style.ERROR(
                f'处理.h5文件{file_path}时出错: {str(e)}'
            ))