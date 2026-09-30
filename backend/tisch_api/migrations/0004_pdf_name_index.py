from django.db import migrations, models


INDEX_NAME = 'scaid_pdf_name_idx'


def add_index(apps, schema_editor):
    model = apps.get_model('tisch_api', 'PDFImage')
    with schema_editor.connection.cursor() as cursor:
        constraints = schema_editor.connection.introspection.get_constraints(cursor, model._meta.db_table)
        if INDEX_NAME in constraints:
            return
        if schema_editor.connection.vendor == 'mysql':
            # Fail quickly if an unrelated transaction holds a metadata lock;
            # index construction itself permits concurrent reads and writes.
            cursor.execute('SET SESSION lock_wait_timeout = 10')
            cursor.execute(
                'ALTER TABLE tisch_api_pdfimage ADD INDEX scaid_pdf_name_idx (name), '
                'ALGORITHM=INPLACE, LOCK=NONE'
            )
        else:
            schema_editor.add_index(model, models.Index(fields=['name'], name=INDEX_NAME))


def remove_index(apps, schema_editor):
    model = apps.get_model('tisch_api', 'PDFImage')
    schema_editor.remove_index(model, models.Index(fields=['name'], name=INDEX_NAME))


class Migration(migrations.Migration):
    dependencies = [('tisch_api', '0003_gene_kegg_enum')]
    atomic = False
    operations = [
        migrations.SeparateDatabaseAndState(
            database_operations=[migrations.RunPython(add_index, remove_index)],
            state_operations=[migrations.AddIndex('pdfimage', models.Index(fields=['name'], name=INDEX_NAME))],
        ),
    ]
