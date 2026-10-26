from django.db import models

from django.db import models
# forms.py
from django import forms


class PDFImage(models.Model):
    name = models.CharField(max_length=255)
    path = models.CharField(max_length=512)
    relative_path = models.CharField(max_length=512)  # Store the relative path.
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
    disease_name = models.CharField(max_length=100, blank=True, null=True, verbose_name='Standard disease name')
    disease_full_name = models.CharField(max_length=255, blank=True, null=True, verbose_name='Full disease name/subtype')
    abbreviation = models.CharField(max_length=50, blank=True, null=True, verbose_name='Disease abbreviation')
    dataset_source = models.CharField(max_length=50, blank=True, null=True, verbose_name='Dataset source')
    dataset_id = models.CharField(max_length=100, blank=True, null=True, verbose_name='Unique dataset identifier')
    tissue = models.CharField(max_length=50, blank=True, null=True, verbose_name='Tissue type')
    sample_size = models.IntegerField(blank=True, null=True, verbose_name='Sample size')
    cell_count = models.IntegerField(blank=True, null=True, verbose_name='Total cell count')
    remarks = models.TextField(blank=True, null=True, verbose_name='Notes')
    batch = models.CharField(max_length=8, blank=True, null=True, verbose_name='Batch date')
    file_name = models.CharField(max_length=8, blank=True, null=True, verbose_name='File name')

    class Meta:
        db_table = 'cell_data'
        verbose_name = 'Cell study data'
        verbose_name_plural = 'Cell study data'

    def __str__(self):
        return f"{self.disease_name} ({self.dataset_id})"

class H5File(models.Model):
    name = models.CharField(max_length=255, verbose_name="File name")
    path = models.CharField(max_length=1024, verbose_name="Full path")
    relative_path = models.CharField(max_length=1024, verbose_name="Relative path")
    create_time = models.DateTimeField(auto_now_add=True)
    file_size = models.BigIntegerField(verbose_name="File size (bytes)")
    last_modified = models.DateTimeField(auto_now=True, verbose_name="Last modified time")

    class Meta:
        verbose_name = "H5 file"
        verbose_name_plural = "H5 files"
        db_table = "h5_files"
        indexes = [
            models.Index(fields=['name'], name='h5file_name_idx'),
            models.Index(fields=['path'], name='h5file_path_idx'),
        ]

    def __str__(self):
        return f"{self.name} ({self.relative_path})"

class CellDataSearchForm(forms.Form):
    disease_name = forms.CharField(label='Disease name', required=False)
    dataset_source = forms.ChoiceField(
        label='Dataset source',
        choices=[
            ('', 'All sources'),
            ('GEO', 'GEO'),
            ('Zenodo', 'Zenodo'),
            ('SCP', 'SCP'),
        ],
        required=False
    )
    tissue = forms.CharField(label='Tissue type', required=False)
    dataset_id = forms.CharField(label='Dataset ID', required=False)
    min_sample_size = forms.IntegerField(label='Minimum sample size', required=False)
    max_sample_size = forms.IntegerField(label='Maximum sample size', required=False)


class AffectedSite(models.Model):
    site_name = models.CharField(max_length=100, unique=True, verbose_name='Affected site name')
    name_en = models.CharField(max_length=100, unique=True, verbose_name='Affected site name in English')
    description = models.TextField(blank=True, null=True, verbose_name='Site description')
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'h_affected_site'
        verbose_name = 'Affected site'
        verbose_name_plural = 'Affected sites'

    def __str__(self):
        return self.site_name


class Disease(models.Model):
    name_en = models.CharField(max_length=100, verbose_name='Disease name in English')
    abbreviation = models.CharField(max_length=20, verbose_name='Disease abbreviation')
    features = models.TextField(blank=True, null=True, verbose_name='Pathological feature description')
    created_at = models.DateTimeField(auto_now_add=True, verbose_name='Created at')
    updated_at = models.DateTimeField(auto_now=True, verbose_name='Updated at')

    class Meta:
        db_table = 'h_disease'
        verbose_name = 'Disease'
        verbose_name_plural = 'Diseases'
        indexes = [
            models.Index(fields=['abbreviation'], name='idx_abbreviation'),
            models.Index(fields=['name_en'], name='idx_name_en'),
        ]

    def __str__(self):
        return f"{self.name_en} ({self.abbreviation})"


class DiseaseSiteRelation(models.Model):
    disease = models.ForeignKey(Disease, on_delete=models.CASCADE, verbose_name='Disease ID')
    site = models.ForeignKey(AffectedSite, on_delete=models.CASCADE, verbose_name='Site ID')
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'h_disease_site_relation'
        verbose_name = 'Disease and site association'
        verbose_name_plural = 'Disease and site associations'
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
        verbose_name='Type (gene, kegg)'
    )
    name = models.CharField(
        max_length=50,
        null=True,
        blank=True,
        verbose_name='Name'
    )

    class Meta:
        verbose_name = 'Gene/KEGG entry'
        verbose_name_plural = 'Gene/KEGG entries'
        db_table = 'gene_kegg_enum'
        indexes = [
            models.Index(fields=['tp'], name='gene_kegg_tp_idx'),
            models.Index(fields=['name'], name='gene_kegg_name_idx'),
        ]

    def __str__(self):
        return f"{self.get_tp_display()}: {self.name}"
