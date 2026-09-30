"""Portable SCAID configuration: environment variables only, no site credentials.

Local defaults use an empty SQLite catalogue. Production mode requires explicit
hosts and enables HTTPS safeguards; use a read-only DB identity for web serving.
"""
import os
from pathlib import Path
from django.core.exceptions import ImproperlyConfigured

BASE_DIR = Path(__file__).resolve().parent.parent
PRODUCTION = os.environ.get('SCAID_DEPLOYMENT', 'local').lower() == 'production'
SECRET_KEY = os.environ.get('SCAID_SECRET_KEY', '')
if len(SECRET_KEY) < 50 or SECRET_KEY.startswith(('django-insecure-', '<GENERATE')):
    raise ImproperlyConfigured('Set SCAID_SECRET_KEY to a new generated secret of at least 50 characters.')
DEBUG = os.environ.get('SCAID_DEBUG', 'false').lower() == 'true'
if PRODUCTION and DEBUG:
    raise ImproperlyConfigured('DEBUG must be false in production.')
ALLOWED_HOSTS = [host.strip() for host in os.environ.get('SCAID_ALLOWED_HOSTS', '' if PRODUCTION else '127.0.0.1,localhost').split(',') if host.strip()]
if not ALLOWED_HOSTS or '*' in ALLOWED_HOSTS:
    raise ImproperlyConfigured('Configure explicit SCAID_ALLOWED_HOSTS.')
INSTALLED_APPS = ['django.contrib.admin', 'django.contrib.auth', 'django.contrib.contenttypes', 'django.contrib.sessions', 'django.contrib.messages', 'django.contrib.staticfiles', 'tisch_api', 'rest_framework', 'csp']
MIDDLEWARE = ['django.middleware.security.SecurityMiddleware', 'django.contrib.sessions.middleware.SessionMiddleware', 'django.middleware.common.CommonMiddleware', 'django.middleware.csrf.CsrfViewMiddleware', 'django.contrib.auth.middleware.AuthenticationMiddleware', 'django.contrib.messages.middleware.MessageMiddleware', 'django.middleware.clickjacking.XFrameOptionsMiddleware', 'csp.middleware.CSPMiddleware']
ROOT_URLCONF = 'tisch_backend.urls'
WSGI_APPLICATION = 'tisch_backend.wsgi.application'
TEMPLATES = [{'BACKEND': 'django.template.backends.django.DjangoTemplates', 'DIRS': [], 'APP_DIRS': True, 'OPTIONS': {'context_processors': ['django.template.context_processors.debug', 'django.template.context_processors.request', 'django.contrib.auth.context_processors.auth', 'django.contrib.messages.context_processors.messages']}}]
REST_FRAMEWORK = {'DEFAULT_PERMISSION_CLASSES': ['rest_framework.permissions.IsAuthenticatedOrReadOnly'], 'DEFAULT_RENDERER_CLASSES': ['rest_framework.renderers.JSONRenderer'], 'DEFAULT_PAGINATION_CLASS': 'rest_framework.pagination.PageNumberPagination'}
CONTENT_SECURITY_POLICY = {'DIRECTIVES': {'default-src': ["'self'"], 'img-src': ["'self'", 'data:', 'blob:']}}
if os.environ.get('SCAID_DB_ENGINE', 'sqlite').lower() == 'mysql':
    for variable in ('SCAID_DB_NAME', 'SCAID_DB_USER', 'SCAID_DB_PASSWORD'):
        if not os.environ.get(variable):
            raise ImproperlyConfigured(f'Configure {variable} outside the source repository.')
    DATABASES = {'default': {'ENGINE': 'django.db.backends.mysql', 'NAME': os.environ['SCAID_DB_NAME'], 'USER': os.environ['SCAID_DB_USER'], 'PASSWORD': os.environ['SCAID_DB_PASSWORD'], 'HOST': os.environ.get('SCAID_DB_HOST', '127.0.0.1'), 'PORT': os.environ.get('SCAID_DB_PORT', '3306'), 'OPTIONS': {'charset': 'utf8mb4'}, 'CONN_MAX_AGE': 60}}
else:
    DATABASES = {'default': {'ENGINE': 'django.db.backends.sqlite3', 'NAME': os.environ.get('SCAID_SQLITE_PATH', str(BASE_DIR / 'db.sqlite3'))}}
AUTH_PASSWORD_VALIDATORS = [{'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator'}, {'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator'}, {'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator'}, {'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator'}]
LANGUAGE_CODE = 'en-us'
TIME_ZONE = 'UTC'
USE_I18N = True
USE_TZ = True
DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'
STATIC_URL = '/static/'
STATIC_ROOT = BASE_DIR / 'staticfiles'
MEDIA_URL = '/media/'
MEDIA_ROOT = Path(os.environ.get('SCAID_MEDIA_ROOT', str(BASE_DIR / 'media')))
FILE_UPLOAD_PERMISSIONS = 0o644
H5_DOWNLOAD_ROOTS = [Path(path).resolve() for path in os.environ.get('SCAID_H5_DOWNLOAD_ROOTS', '/srv/scaid/data/Data.H5').split(os.pathsep) if path]
SCAID_FIGURE_ROOT = Path(os.environ.get('SCAID_FIGURE_ROOT', '/srv/scaid/data/Fig'))
SCAID_EXTRA_FIGURE_ROOTS = [Path(path) for path in os.environ.get('SCAID_EXTRA_FIGURE_ROOTS', '').split(os.pathsep) if path]
H5_DOWNLOAD_ACCEL_PREFIX = os.environ.get('SCAID_H5_ACCEL_PREFIX', '')
THUMBNAIL_ROOT = Path(os.environ.get('SCAID_THUMBNAIL_ROOT', str(MEDIA_ROOT / 'thumbs')))
THUMBNAIL_ACCEL_PREFIX = os.environ.get('SCAID_THUMB_ACCEL_PREFIX', '')
THUMBNAIL_MAX_SIDE = 900
CATALOG_CACHE_SECONDS = 300 if PRODUCTION else 0
SESSION_COOKIE_SECURE = PRODUCTION
CSRF_COOKIE_SECURE = PRODUCTION
SECURE_SSL_REDIRECT = PRODUCTION
SECURE_HSTS_SECONDS = 31536000 if PRODUCTION else 0
SECURE_HSTS_INCLUDE_SUBDOMAINS = PRODUCTION
SECURE_HSTS_PRELOAD = PRODUCTION
SECURE_CONTENT_TYPE_NOSNIFF = True
X_FRAME_OPTIONS = 'DENY'
# Trust this header only when the app is reachable exclusively through a proxy
# that strips untrusted client forwarding headers and sets this value itself.
if os.environ.get('SCAID_TRUST_PROXY_SSL', 'false').lower() == 'true':
    SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')
