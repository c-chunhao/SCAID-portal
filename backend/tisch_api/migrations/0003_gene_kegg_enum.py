from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [('tisch_api', '0002_populate_relative_path')]

    operations = [
        migrations.CreateModel(
            name='GeneKeggEnum',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('tp', models.CharField(blank=True, choices=[('gene', 'Gene'), ('kegg', 'KEGG')], max_length=20, null=True, verbose_name='Type (gene, kegg)')),
                ('name', models.CharField(blank=True, max_length=50, null=True, verbose_name='Name')),
            ],
            options={
                'db_table': 'gene_kegg_enum',
                'verbose_name': 'Gene/KEGG entry',
                'verbose_name_plural': 'Gene/KEGG entries',
                'indexes': [
                    models.Index(fields=['tp'], name='gene_kegg_tp_idx'),
                    models.Index(fields=['name'], name='gene_kegg_name_idx'),
                ],
            },
        ),
    ]
