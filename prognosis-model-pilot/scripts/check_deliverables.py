#!/usr/bin/env python3
"""Check a project's declared deliverables, not scientific validity or page layout.

Requires Python 3.10+ and only its standard library. No analysis, downloads, or file mutation.
Exit 0: requested files satisfy these checks. Exit 1: missing/incomplete/invalid.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import struct
import sys
import xml.etree.ElementTree as ET
import zipfile


FIELDS = {
    'artifact_id', 'scope', 'module', 'kind', 'figure_id', 'required', 'status',
    'path', 'reference_project', 'reference_path', 'source_result',
    'report_anchor', 'reason',
}
SCOPES = {'proposal', 'analysis', 'report', 'shared'}
STATUSES = {'planned', 'completed', 'blocked', 'failed', 'not_applicable', 'not_requested'}
PACKAGE_PARTS = {
    '.docx': 'word/document.xml', '.pptx': 'ppt/presentation.xml',
    '.xlsx': 'xl/workbook.xml',
}


def file_problem(path: Path) -> str | None:
    """Lightweight integrity check; rendering and model loading remain separate."""
    try:
        if not path.is_file():
            return 'file_missing'
        if path.stat().st_size == 0:
            return 'file_empty'
        suffix = path.suffix.lower()
        with path.open('rb') as stream:
            head = stream.read(32)
        if suffix == '.png':
            if head[:8] != b'\x89PNG\r\n\x1a\n' or head[12:16] != b'IHDR':
                return 'invalid_png_header'
            if len(head) < 24 or min(struct.unpack('>II', head[16:24])) <= 0:
                return 'invalid_png_dimensions'
        elif suffix == '.pdf' and not head.startswith(b'%PDF-'):
            return 'invalid_pdf_header'
        elif suffix in {'.jpg', '.jpeg'} and not head.startswith(b'\xff\xd8\xff'):
            return 'invalid_jpeg_header'
        elif suffix == '.svg':
            if ET.parse(path).getroot().tag.rsplit('}', 1)[-1] != 'svg':
                return 'invalid_svg_root'
        elif suffix in PACKAGE_PARTS:
            with zipfile.ZipFile(path) as package:
                if PACKAGE_PARTS[suffix] not in package.namelist():
                    return 'missing_office_main_part'
                if package.testzip() is not None:
                    return 'corrupt_office_package'
                ET.fromstring(package.read(PACKAGE_PARTS[suffix]))
    except (OSError, ValueError, struct.error, ET.ParseError, zipfile.BadZipFile) as exc:
        return f'unreadable_file:{type(exc).__name__}'
    return None


def extensions(value: str) -> set[str]:
    return {part.strip().lower().lstrip('.') for part in value.split(',') if part.strip()}


def audit(root: Path, manifest: Path, scope: str, figure_formats: set[str],
          document_formats: set[str], workflow_formats: set[str]) -> dict:
    issues = []

    def issue(code, row=None, **details):
        issues.append({'code': code, **({'artifact_id': row.get('artifact_id', '')} if row else {}), **details})

    try:
        with manifest.open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            missing = FIELDS - set(reader.fieldnames or [])
            if missing:
                return {'status': 'incomplete', 'scope': scope,
                        'issues': [{'code': 'missing_manifest_columns', 'columns': sorted(missing)}]}
            rows = [{key: (value or '').strip() for key, value in row.items() if key is not None}
                    for row in reader if any((value or '').strip() for key, value in row.items() if key is not None)]
    except (OSError, UnicodeError, csv.Error) as exc:
        return {'status': 'incomplete', 'scope': scope,
                'issues': [{'code': 'manifest_unreadable', 'error_type': type(exc).__name__}]}
    if not rows:
        issue('manifest_has_no_items')

    selected = []
    seen = set()
    for row in rows:
        artifact = row['artifact_id']
        if not artifact or artifact in seen:
            issue('missing_or_duplicate_artifact_id', row)
        seen.add(artifact)
        if row['scope'] not in SCOPES:
            issue('invalid_scope', row, value=row['scope'])
        if row['required'].lower() not in {'true', 'false'}:
            issue('invalid_required_value', row)
        if row['status'] not in STATUSES:
            issue('invalid_status', row, value=row['status'])
        if row['status'] in {'not_applicable', 'not_requested', 'blocked', 'failed'} and not row['reason']:
            issue('missing_status_reason', row)
        if scope == 'full' or row['scope'] in {scope, 'shared'}:
            selected.append(row)
    if not selected:
        issue('no_items_in_requested_scope')

    required_rows = [row for row in selected if row['required'].lower() == 'true']
    if not required_rows:
        issue('no_required_items_in_requested_scope')
    completed = 0
    groups = {}
    for row in selected:
        required = row['required'].lower() == 'true'
        if required and row['status'] != 'completed':
            issue('required_item_incomplete', row, state=row['status'], reason=row['reason'])
        if required and row['kind'] in {'figure', 'workflow'}:
            if not row['figure_id']:
                issue('missing_figure_id', row)
            else:
                key = (row['kind'], row['figure_id'])
                groups.setdefault(key, set()).add(Path(row['path']).suffix.lower().lstrip('.'))
        if row['status'] != 'completed':
            continue
        completed += 1
        relative = Path(row['path'])
        target = (root / relative).resolve()
        if not row['path'] or relative.is_absolute() or not target.is_relative_to(root):
            issue('output_path_must_be_project_relative', row, path=row['path'])
            continue
        problem = file_problem(target)
        if problem:
            issue(problem, row, path=row['path'])
        sources = [part.strip() for part in row['source_result'].split('|') if part.strip()]
        if row['kind'] == 'figure' and not sources:
            issue('figure_missing_source_result', row)
        distinct_sources = 0
        for source in sources:
            source_path = Path(source.split('#', 1)[0])
            source_path = (root / source_path).resolve()
            if source_path == target:
                issue('source_result_points_to_output', row)
                continue
            distinct_sources += 1
            try:
                available = source_path.is_file() and source_path.stat().st_size > 0
            except OSError:
                available = False
            if not available:
                issue('source_result_missing_or_empty', row, source=source)
        if row['kind'] == 'figure' and sources and not distinct_sources:
            issue('figure_has_no_distinct_source', row)

    for (kind, identifier), actual in groups.items():
        needed = figure_formats if kind == 'figure' else workflow_formats
        if needed - actual:
            issue('missing_format_in_manifest', kind=kind, figure_id=identifier,
                  missing_formats=sorted(needed - actual))

    doc_scopes = {'proposal', 'report'} if scope == 'full' else ({scope} if scope in {'proposal', 'report'} else set())
    for doc_scope in sorted(doc_scopes):
        actual = {Path(row['path']).suffix.lower().lstrip('.') for row in required_rows
                  if row['scope'] == doc_scope and row['kind'] == 'document'}
        if document_formats - actual:
            issue('missing_formal_document_format', document_scope=doc_scope,
                  missing_formats=sorted(document_formats - actual))
    if scope in {'full', 'proposal'} and not any(row['kind'] == 'workflow' for row in required_rows):
        issue('missing_required_workflow')
    if scope == 'full':
        for kind in ('code', 'table', 'figure', 'model'):
            if not any(row['kind'] == kind and row['scope'] == 'analysis' for row in required_rows):
                issue('missing_analysis_deliverable_category', kind=kind)

    return {
        'status': 'pass' if not issues else 'incomplete', 'scope': scope,
        'manifest_rows': len(rows), 'selected_rows': len(selected),
        'required_rows': len(required_rows), 'completed_rows': completed,
        'issues': issues,
        'limits': 'File/state/format/source checks only. Scientific validity, R object loading, image content and Word/PDF layout require separate review.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--manifest', required=True, type=Path)
    parser.add_argument('--scope', choices=['full', 'proposal', 'analysis', 'report'], default='full')
    parser.add_argument('--figure-formats', default='png,pdf')
    parser.add_argument('--document-formats', default='docx,pdf')
    parser.add_argument('--workflow-formats', default='png,pptx')
    args = parser.parse_args()
    result = audit(args.root.resolve(), args.manifest.resolve(), args.scope,
                   extensions(args.figure_formats), extensions(args.document_formats),
                   extensions(args.workflow_formats))
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    sys.exit(main())
