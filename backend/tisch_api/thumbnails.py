"""On-demand JPEG previews for the large (3300x3000, ~1 MB) figure PNGs.

The public link runs at a few hundred KB/s, so a results page of twelve
full-size PNGs takes close to a minute. Previews are generated once per figure
(~80 KB, 900 px), cached under THUMBNAIL_ROOT and then served by Nginx through
an internal X-Accel location. Originals stay available at their /system/media/
URLs. Only files inside the configured figure roots are ever opened.
"""
from __future__ import annotations

import os
import tempfile
import fcntl
from contextlib import contextmanager
from pathlib import Path

from django.conf import settings
from PIL import Image

PREVIEW_SUFFIXES = {'.png'}


def figure_roots():
    roots = [Path(settings.SCAID_FIGURE_ROOT), Path(settings.MEDIA_ROOT) / 'pdf_images']
    roots += [Path(root) for root in getattr(settings, 'SCAID_EXTRA_FIGURE_ROOTS', []) if root]
    return [root.resolve() for root in roots]


def resolved_figure_path(record):
    """Absolute source path for a PDFImage row, or None if it is outside the roots."""
    name = (record.image.name if record.image else '') or ''
    if not name:
        return None
    candidates = [Path(name)] if os.path.isabs(name) else [Path(settings.MEDIA_ROOT) / name, Path(name)]
    roots = figure_roots()
    if not os.path.isabs(name):
        # Older rendered records store pdf_images/<folder>/<page>.jpg under
        # an explicitly configured legacy figure root. Never infer private roots.
        relative = Path(name)
        if relative.parts and relative.parts[0] == 'pdf_images':
            relative = Path(*relative.parts[1:])
        candidates.extend(root / relative for root in roots)
    for candidate in candidates:
        try:
            resolved = candidate.resolve(strict=True)
        except (OSError, RuntimeError):
            continue
        if resolved.is_file() and any(resolved.is_relative_to(root) for root in roots):
            return resolved
    return None


def wants_preview(record):
    return Path(record.image.name or '').suffix.lower() in PREVIEW_SUFFIXES if record.image else False


def thumbnail_relative(pk):
    # Two-level fan-out keeps directories small (1.16M figures in the catalogue).
    return Path(f'{pk % 1000:03d}') / f'{pk}.jpg'


def make_thumbnail(source: Path, target: Path, max_side=None):
    max_side = max_side or getattr(settings, 'THUMBNAIL_MAX_SIDE', 900)
    target.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(source) as image:
        image.draft('RGB', (max_side * 2, max_side * 2))
        if image.mode in ('RGBA', 'LA', 'P'):
            rgba = image.convert('RGBA')
            flat = Image.new('RGB', rgba.size, (255, 255, 255))
            flat.paste(rgba, mask=rgba.getchannel('A'))
            image = flat
        else:
            image = image.convert('RGB')
        image.thumbnail((max_side, max_side), Image.LANCZOS)
        # Each writer owns a distinct temporary file on the target filesystem.
        # Readers see a complete old or new JPEG, never a partial write.
        descriptor, temporary = tempfile.mkstemp(prefix=f'.{target.name}.', suffix='.tmp', dir=target.parent)
        tmp = Path(temporary)
        try:
            with os.fdopen(descriptor, 'wb') as output:
                image.save(output, 'JPEG', quality=82, optimize=True, progressive=True)
            os.replace(tmp, target)
        finally:
            tmp.unlink(missing_ok=True)
    return target


@contextmanager
def thumbnail_lock(target):
    """Serialize a cache miss across Gunicorn processes on this Linux host."""
    target.parent.mkdir(parents=True, exist_ok=True)
    # Keep the lock inode in place: unlinking it would let another process lock
    # a different inode while an existing writer still owns the original one.
    with target.with_suffix('.jpg.lock').open('a+b') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def is_fresh(source, target):
    try:
        return target.stat().st_mtime_ns >= source.stat().st_mtime_ns
    except FileNotFoundError:
        return False


def ensure_thumbnail(record):
    """Return the cached preview path; recheck freshness after acquiring the lock."""
    source = resolved_figure_path(record)
    if source is None:
        return None
    target = Path(settings.THUMBNAIL_ROOT) / thumbnail_relative(record.pk)
    if not is_fresh(source, target):
        with thumbnail_lock(target):
            if not is_fresh(source, target):
                make_thumbnail(source, target)
    return target
