"""Index relative_path so folder-scoped figure queries stop scanning the table.

The production schema was partly created by hand, so the database operation
checks SHOW INDEX first and only creates the index when it is missing.
"""
from django.db import migrations, models

INDEX = 'scaid_pdf_relpath_idx'
TABLE = 'tisch_api_pdfimage'


def create_index(apps, schema_editor):
    # Early migrations omitted this manually provisioned production column.
    # A clean install must create it before building the index; existing tables
    # keep their current values and column unchanged.
    model = apps.get_model('tisch_api', 'PDFImage')
    with schema_editor.connection.cursor() as cursor:
        columns = {column.name for column in schema_editor.connection.introspection.get_table_description(cursor, TABLE)}
    if 'relative_path' not in columns:
        field = models.CharField(max_length=512, default='')
        field.set_attributes_from_name('relative_path')
        schema_editor.add_field(model, field)
    with schema_editor.connection.cursor() as cursor:
        if schema_editor.connection.vendor == 'mysql':
            cursor.execute(f"SHOW INDEX FROM `{TABLE}` WHERE Key_name = %s", [INDEX])
            if cursor.fetchone():
                return
            cursor.execute(f"CREATE INDEX `{INDEX}` ON `{TABLE}` (`relative_path`)")
        else:
            cursor.execute(f'CREATE INDEX IF NOT EXISTS "{INDEX}" ON "{TABLE}" ("relative_path")')


def drop_index(apps, schema_editor):
    with schema_editor.connection.cursor() as cursor:
        if schema_editor.connection.vendor == 'mysql':
            cursor.execute(f"SHOW INDEX FROM `{TABLE}` WHERE Key_name = %s", [INDEX])
            if cursor.fetchone():
                cursor.execute(f"DROP INDEX `{INDEX}` ON `{TABLE}`")
        else:
            cursor.execute(f'DROP INDEX IF EXISTS "{INDEX}"')


class Migration(migrations.Migration):
    dependencies = [('tisch_api', '0004_pdf_name_index')]

    operations = [
        migrations.SeparateDatabaseAndState(
            state_operations=[
                migrations.AddField('pdfimage', 'relative_path', models.CharField(max_length=512)),
                migrations.AddIndex(
                    model_name='pdfimage',
                    index=models.Index(fields=['relative_path'], name=INDEX),
                ),
            ],
            database_operations=[migrations.RunPython(create_index, drop_index)],
        ),
    ]
