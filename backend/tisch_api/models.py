from django.db import models

from django.db import models
# forms.py
from django import forms


class PDFImage(models.Model):
    name = models.CharField(max_length=255)
    path = models.CharField(max_length=512)
    relative_path = models.CharField(max_length=512)  # 新增字段保存相对路径
    page_number = models.IntegerField()
    image = models.FileField(upload_to='pdf_images/', null=True, blank=True)
    extracted_at = models.DateTimeField(auto_now_add=True)
    image_width = models.CharField(max_length=512)
    image_height = models.CharField(max_length=512)
    is_pdf = models.CharField(max_length=512)
    create_time = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'tisch_api_pdfimage'
        indexes = [
            models.Index(fields=['name'], name='scaid_pdf_name_idx'),
            # Full image paths are not indexed: catalogue queries use the
            # basename and relative folder indexes, including legacy aliases.
            # Folder-scoped figure queries match this column exactly or by prefix.
            models.Index(fields=['relative_path'], name='scaid_pdf_relpath_idx'),
            # A PostgreSQL trigram GinIndex used to be declared here; the
            # production database is MySQL and the migrations never created it.
        ]

    def __str__(self):
        return f"{self.name} - Page {self.page_number}"

class CellData(models.Model):
    disease_name = models.CharField(max_length=100, blank=True, null=True, verbose_name='疾病标准名称')
    disease_full_name = models.CharField(max_length=255, blank=True, null=True, verbose_name='疾病完整名称/子类型')
    abbreviation = models.CharField(max_length=50, blank=True, null=True, verbose_name='疾病缩写')
    dataset_source = models.CharField(max_length=50, blank=True, null=True, verbose_name='数据集来源')
    dataset_id = models.CharField(max_length=100, blank=True, null=True, verbose_name='数据集唯一标识')
    tissue = models.CharField(max_length=50, blank=True, null=True, verbose_name='组织类型')
    sample_size = models.IntegerField(blank=True, null=True, verbose_name='样本量')
    cell_count = models.IntegerField(blank=True, null=True, verbose_name='细胞总数')
    remarks = models.TextField(blank=True, null=True, verbose_name='备注信息')
    batch = models.CharField(max_length=8, blank=True, null=True, verbose_name='批次日期')
    file_name = models.CharField(max_length=8, blank=True, null=True, verbose_name='文件名称')

    class Meta:
        db_table = 'cell_data'
        verbose_name = '细胞研究数据'
        verbose_name_plural = '细胞研究数据'

    def __str__(self):
        return f"{self.disease_name} ({self.dataset_id})"

class H5File(models.Model):
    name = models.CharField(max_length=255, verbose_name="文件名")
    path = models.CharField(max_length=1024, verbose_name="完整路径")
    relative_path = models.CharField(max_length=1024, verbose_name="相对路径")
    create_time = models.DateTimeField(auto_now_add=True)
    file_size = models.BigIntegerField(verbose_name="文件大小(字节)")
    last_modified = models.DateTimeField(auto_now=True, verbose_name="最后修改时间")

    class Meta:
        verbose_name = "H5文件"
        verbose_name_plural = "H5文件"
        db_table = "h5_files"
        indexes = [
            models.Index(fields=['name'], name='h5file_name_idx'),
            models.Index(fields=['path'], name='h5file_path_idx'),
        ]

    def __str__(self):
        return f"{self.name} ({self.relative_path})"

class CellDataSearchForm(forms.Form):
    disease_name = forms.CharField(label='疾病名称', required=False)
    dataset_source = forms.ChoiceField(
        label='数据集来源',
        choices=[
            ('', '所有来源'),
            ('GEO', 'GEO'),
            ('Zenodo', 'Zenodo'),
            ('SCP', 'SCP'),
        ],
        required=False
    )
    tissue = forms.CharField(label='组织类型', required=False)
    dataset_id = forms.CharField(label='数据集ID', required=False)
    min_sample_size = forms.IntegerField(label='最小样本量', required=False)
    max_sample_size = forms.IntegerField(label='最大样本量', required=False)


class AffectedSite(models.Model):
    site_name = models.CharField(max_length=100, unique=True, verbose_name='发病部位名称')
    name_en = models.CharField(max_length=100, unique=True, verbose_name='发病部位英文名称')
    description = models.TextField(blank=True, null=True, verbose_name='部位描述')
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'h_affected_site'
        verbose_name = '发病部位'
        verbose_name_plural = '发病部位'

    def __str__(self):
        return self.site_name


class Disease(models.Model):
    name_en = models.CharField(max_length=100, verbose_name='疾病英文名称')
    abbreviation = models.CharField(max_length=20, verbose_name='疾病缩写')
    features = models.TextField(blank=True, null=True, verbose_name='病变特征描述')
    created_at = models.DateTimeField(auto_now_add=True, verbose_name='创建时间')
    updated_at = models.DateTimeField(auto_now=True, verbose_name='更新时间')

    class Meta:
        db_table = 'h_disease'
        verbose_name = '疾病'
        verbose_name_plural = '疾病'
        indexes = [
            models.Index(fields=['abbreviation'], name='idx_abbreviation'),
            models.Index(fields=['name_en'], name='idx_name_en'),
        ]

    def __str__(self):
        return f"{self.name_en} ({self.abbreviation})"


class DiseaseSiteRelation(models.Model):
    disease = models.ForeignKey(Disease, on_delete=models.CASCADE, verbose_name='疾病ID')
    site = models.ForeignKey(AffectedSite, on_delete=models.CASCADE, verbose_name='部位ID')
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'h_disease_site_relation'
        verbose_name = '疾病与发病部位关联'
        verbose_name_plural = '疾病与发病部位关联'
        unique_together = (('disease', 'site'),)

    def __str__(self):
        return f"{self.disease} - {self.site}"


class GeneKeggEnum(models.Model):
    TP_CHOICES = [
        ('gene', 'Gene'),
        ('kegg', 'KEGG'),
    ]

    tp = models.CharField(
        max_length=20,
        choices=TP_CHOICES,
        null=True,
        blank=True,
        verbose_name='类型（gene、kegg）'
    )
    name = models.CharField(
        max_length=50,
        null=True,
        blank=True,
        verbose_name='名称'
    )

    class Meta:
        verbose_name = 'Gene/KEGG枚举'
        verbose_name_plural = 'Gene/KEGG枚举'
        db_table = 'gene_kegg_enum'
        indexes = [
            models.Index(fields=['tp'], name='gene_kegg_tp_idx'),
            models.Index(fields=['name'], name='gene_kegg_name_idx'),
        ]

    def __str__(self):
        return f"{self.get_tp_display()}: {self.name}"
