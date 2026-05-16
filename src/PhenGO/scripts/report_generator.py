#!/usr/bin/env python3
"""HTML summariser for PhenGO result directories.

Given a PhenGO output directory (or a PhenGO-Predict output directory),
parse every known artefact and generate a self-contained HTML report.

Recognised artefacts
--------------------
PhenGO core
  PhenGO_params.txt               – run parameters
  *.arff                           – gene / GO-term matrix  → dataset stats
  GO_term_validation_report.txt   – GO term validity
  Unique_GO_Nodes.txt             – distinct GO nodes used
  go_enrichment/                  – enrichment input files  → row counts
  GO_Children&Parents.json        – GO tree (depth summary)
  PhenGO.log                      – last lines

PhenGO-Predict (Predict/ sub-dir or standalone)
  PhenGO_Predict.log              – model metrics
  */final_report.txt              – per-model text report
  */roc_curve.png, */training_history.png, */feature_importance*.{png,csv}
"""
import os as _os, sys as _sys
if __name__ == "__main__" and not __package__:
    import importlib.util as _ilu
    _here   = _os.path.dirname(_os.path.abspath(__file__))
    _src    = _os.path.normpath(_os.path.join(_here, '..', '..'))
    _phengo = _os.path.dirname(_here)
    if _src not in _sys.path:
        _sys.path.insert(0, _src)
    for _n, _f, _d in [
        ('PhenGO',         _os.path.join(_phengo, '__init__.py'), _phengo),
        ('PhenGO.scripts', _os.path.join(_here,   '__init__.py'), _here),
    ]:
        if _n not in _sys.modules:
            _sp = _ilu.spec_from_file_location(_n, _f, submodule_search_locations=[_d])
            _mo = _ilu.module_from_spec(_sp)
            _mo.__path__ = [_d]
            _sys.modules[_n] = _mo
            _sp.loader.exec_module(_mo)
    __package__ = 'PhenGO.scripts'
    del _ilu, _here, _src, _phengo, _n, _f, _d, _sp, _mo

import os
import re
import json
import base64
import argparse
import logging
import statistics
from collections import Counter

from ..constants import configure_logger

logger = logging.getLogger('PhenGO.report_generator')

# ── helpers ──────────────────────────────────────────────────────────────────

CSS = """
body{font-family:Arial,sans-serif;max-width:1200px;margin:auto;padding:20px;background:#f9f9f9}
h1{background:#2c3e50;color:#fff;padding:15px;border-radius:6px}
h2{background:#34495e;color:#fff;padding:10px;border-radius:4px;margin-top:30px}
h3{color:#2c3e50;border-bottom:2px solid #bdc3c7;padding-bottom:5px}
table{border-collapse:collapse;width:100%;margin:10px 0}
th{background:#2c3e50;color:#fff;padding:8px 12px;text-align:left}
td{padding:6px 12px;border-bottom:1px solid #ddd}
tr:nth-child(even){background:#ecf0f1}
.metric-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(200px,1fr));gap:10px;margin:10px 0}
.metric-card{background:#fff;border:1px solid #bdc3c7;border-radius:6px;padding:12px;text-align:center}
.metric-card .val{font-size:1.8em;font-weight:bold;color:#2980b9}
.metric-card .lbl{font-size:.85em;color:#7f8c8d;margin-top:4px}
details>summary{cursor:pointer;padding:6px;background:#ecf0f1;border-radius:4px;margin:4px 0}
pre{background:#2c3e50;color:#ecf0f1;padding:12px;border-radius:4px;overflow-x:auto;font-size:.85em}
img{max-width:100%;border:1px solid #bdc3c7;border-radius:4px;margin:6px 0}
.warn{color:#c0392b;font-weight:bold}
.ok{color:#27ae60;font-weight:bold}
"""


def embed_image(path: str) -> str | None:
    try:
        with open(path, 'rb') as fh:
            data = fh.read()
        b64 = base64.b64encode(data).decode('ascii')
        ext = os.path.splitext(path)[1].lstrip('.').lower()
        mime = {'jpg': 'jpeg', 'jpeg': 'jpeg', 'png': 'png', 'svg': 'svg+xml'}.get(ext, ext)
        return f'data:image/{mime};base64,{b64}'
    except Exception:
        return None


def _read(path: str, errors='replace') -> str:
    try:
        with open(path, encoding='utf-8', errors=errors) as fh:
            return fh.read()
    except Exception:
        return ''


def metric_card(value, label: str) -> str:
    return (f'<div class="metric-card">'
            f'<div class="val">{value}</div>'
            f'<div class="lbl">{label}</div>'
            f'</div>')


def section(title: str, body: str) -> str:
    return f'<h2>{title}</h2>\n{body}\n'


def kv_table(rows: list[tuple]) -> str:
    """Render a two-column key/value table."""
    html = '<table><tr><th>Parameter</th><th>Value</th></tr>'
    for k, v in rows:
        html += f'<tr><td>{k}</td><td>{v}</td></tr>'
    return html + '</table>'


# ── PhenGO core parsers ───────────────────────────────────────────────────────

def parse_params(path: str) -> dict:
    params = {}
    for line in _read(path).splitlines():
        if ':' in line:
            k, _, v = line.partition(':')
            params[k.strip()] = v.strip()
    return params


def parse_arff_stats(path: str) -> dict:
    """Return gene count, GO term count, class distribution and avg GO per gene."""
    attrs = []
    data_lines = []
    in_data = False
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('%'):
                continue
            if line.lower().startswith('@attribute'):
                parts = line.split()
                attrs.append(parts[1].strip("'\""))
            elif line.lower() == '@data':
                in_data = True
            elif in_data:
                data_lines.append(line)

    go_terms = attrs[1:-1]
    class_counts: Counter = Counter()
    go_per_gene: list[int] = []

    for line in data_lines:
        parts = [p.strip() for p in line.split(',')]
        if len(parts) < 3:
            continue
        label = parts[-1].lower()
        class_counts[label] += 1
        n_go = sum(1 for v in parts[1:-1] if v not in ('0', 'NA', '', '0.0'))
        go_per_gene.append(n_go)

    return {
        'n_genes': len(data_lines),
        'n_go_terms': len(go_terms),
        'class_counts': dict(class_counts),
        'avg_go_per_gene': round(statistics.mean(go_per_gene), 2) if go_per_gene else 0,
        'median_go_per_gene': round(statistics.median(go_per_gene), 1) if go_per_gene else 0,
        'max_go_per_gene': max(go_per_gene) if go_per_gene else 0,
        'pct_zero_go': round(go_per_gene.count(0) / len(go_per_gene) * 100, 1) if go_per_gene else 0,
    }


def parse_go_validation(path: str) -> dict:
    text = _read(path)
    result = {}
    for line in text.splitlines():
        for key in ('Total GO terms checked', 'Valid terms', 'Missing terms',
                    'Obsolete terms found', 'Genes affected'):
            if line.startswith(key):
                _, _, val = line.partition(':')
                result[key] = val.strip()
    return result


def parse_unique_go_nodes(path: str) -> int:
    return sum(1 for line in _read(path).splitlines() if line.strip())


def enrichment_file_stats(enrich_dir: str) -> list[tuple]:
    rows = []
    for fname in sorted(os.listdir(enrich_dir)):
        fpath = os.path.join(enrich_dir, fname)
        if os.path.isfile(fpath):
            try:
                with open(fpath, 'rb') as fh:
                    n = sum(1 for _ in fh)
                rows.append((fname, f'{n:,} lines', f'{os.path.getsize(fpath) / 1024:.1f} KB'))
            except Exception:
                rows.append((fname, '?', '?'))
    return rows


def go_depth_stats(json_path: str) -> dict | None:
    """Compute depth statistics from GO_Children&Parents.json.

    Supports two formats:
      1. {GO_ID: {"depth": N, ...}} – depth pre-computed
      2. {GO_ID: {"p": [...parents...], "c": [...children...]}} – tree structure,
         depth computed via BFS from roots (nodes with no parents).
    """
    try:
        with open(json_path, 'r', errors='replace') as fh:
            data = json.load(fh)
    except Exception:
        return None

    if not isinstance(data, dict):
        return None

    depths: list[int] = []

    # Heuristic 1: pre-computed depth field
    sample = next(iter(data.values()), {}) if data else {}
    if isinstance(sample, dict) and any(k in sample for k in ('depth', 'level', 'min_depth', 'max_depth')):
        for val in data.values():
            if isinstance(val, dict):
                for key in ('depth', 'level', 'min_depth', 'max_depth'):
                    if key in val and isinstance(val[key], (int, float)):
                        depths.append(int(val[key]))
                        break

    # Heuristic 2: p/c tree structure → BFS to find min-depth from root
    if not depths and isinstance(sample, dict) and ('p' in sample or 'c' in sample):
        from collections import deque
        # Roots = nodes whose parent list is empty
        depth_map: dict[str, int] = {}
        roots = [go_id for go_id, val in data.items()
                 if isinstance(val, dict) and not val.get('p')]
        queue: deque = deque()
        for r in roots:
            depth_map[r] = 0
            queue.append(r)
        while queue:
            node = queue.popleft()
            children = data.get(node, {}).get('c', [])
            for child in children:
                if child not in depth_map and child in data:
                    depth_map[child] = depth_map[node] + 1
                    queue.append(child)
        # Any term not reached (disconnected) gets depth from its minimum parent depth
        for go_id, val in data.items():
            if go_id not in depth_map and isinstance(val, dict):
                parent_depths = [depth_map[p] for p in val.get('p', []) if p in depth_map]
                if parent_depths:
                    depth_map[go_id] = min(parent_depths) + 1
        depths = list(depth_map.values())

    if not depths:
        return None

    return {
        'n_terms': len(depths),
        'mean_depth': round(statistics.mean(depths), 2),
        'median_depth': round(statistics.median(depths), 1),
        'min_depth': min(depths),
        'max_depth': max(depths),
        'depth_distribution': dict(sorted(Counter(depths).items())),
    }


# ── PhenGO-Predict parsers ────────────────────────────────────────────────────

def parse_predict_log(path: str) -> dict:
    """Extract headline metrics from PhenGO_Predict.log."""
    text = _read(path)
    result: dict = {}

    # Look for lines like: Loss=0.312  Accuracy=0.831  AUC=0.916
    summary_match = re.search(
        r'Loss=([\d.]+)\s+Accuracy=([\d.]+)\s+AUC=([\d.]+)', text)
    if summary_match:
        result['Loss'] = summary_match.group(1)
        result['Accuracy'] = summary_match.group(2)
        result['AUC'] = summary_match.group(3)

    # Also try individual lines
    for metric in ('Loss', 'Accuracy', 'Precision', 'Recall', 'AUC', 'F1'):
        m = re.search(rf'{metric}[=:\s]+([\d.]+)', text, re.IGNORECASE)
        if m and metric not in result:
            result[metric] = m.group(1)

    # Threshold line
    m = re.search(r'Using threshold ([\d.]+)', text)
    if m:
        result['Threshold'] = m.group(1)

    # Training set / test set sizes
    for tag in ('Training set', 'Test set'):
        m = re.search(rf'{tag}:\s*([\d,]+)\s*samples', text)
        if m:
            result[tag] = m.group(1)

    return result


def parse_final_report(path: str) -> str:
    return _read(path)


def find_images(directory: str, recurse: bool = True) -> list[str]:
    imgs = []
    ext = {'.png', '.jpg', '.jpeg', '.svg'}
    if recurse:
        for root, _, files in os.walk(directory):
            for f in sorted(files):
                if os.path.splitext(f)[1].lower() in ext:
                    imgs.append(os.path.join(root, f))
    else:
        for f in sorted(os.listdir(directory)):
            if os.path.splitext(f)[1].lower() in ext:
                imgs.append(os.path.join(directory, f))
    return imgs


# ── HTML builders ─────────────────────────────────────────────────────────────

def build_params_section(params: dict) -> str:
    rows = [(k, v) for k, v in params.items()]
    return section('Run Parameters', kv_table(rows))


def build_arff_section(stats: dict) -> str:
    cards = [
        metric_card(f"{stats['n_genes']:,}", 'Total genes'),
        metric_card(f"{stats['n_go_terms']:,}", 'GO feature columns'),
        metric_card(f"{stats['avg_go_per_gene']}", 'Mean GO terms / gene'),
        metric_card(f"{stats['median_go_per_gene']}", 'Median GO terms / gene'),
        metric_card(f"{stats['max_go_per_gene']}", 'Max GO terms / gene'),
        metric_card(f"{stats['pct_zero_go']}%", 'Genes with 0 GO terms'),
    ]
    for cls, cnt in sorted(stats['class_counts'].items()):
        pct = round(cnt / stats['n_genes'] * 100, 1) if stats['n_genes'] else 0
        cards.append(metric_card(f"{cnt:,} ({pct}%)", cls.capitalize()))

    html = '<div class="metric-grid">' + ''.join(cards) + '</div>'
    return section('ARFF Dataset Statistics', html)


def build_go_validation_section(val: dict) -> str:
    rows = list(val.items())
    missing = int(val.get('Missing terms', 0) or 0)
    obsolete = int(val.get('Obsolete terms found', 0) or 0)
    status = '<span class="ok">✔ All terms valid</span>' if (missing + obsolete) == 0 \
        else f'<span class="warn">⚠ {missing} missing, {obsolete} obsolete</span>'
    return section('GO Term Validation', status + kv_table(rows))


def build_go_depth_section(depth: dict) -> str:
    cards = [
        metric_card(f"{depth['n_terms']:,}", 'GO terms with depth info'),
        metric_card(f"{depth['mean_depth']}", 'Mean depth'),
        metric_card(f"{depth['median_depth']}", 'Median depth'),
        metric_card(depth['min_depth'], 'Min depth'),
        metric_card(depth['max_depth'], 'Max depth'),
    ]
    dist = depth['depth_distribution']
    table_rows = [(str(d), f"{cnt:,}") for d, cnt in sorted(dist.items())]
    html = '<div class="metric-grid">' + ''.join(cards) + '</div>'
    html += '<details><summary>Depth distribution</summary>'
    html += '<table><tr><th>Depth</th><th>GO terms</th></tr>'
    for d, cnt in table_rows:
        html += f'<tr><td>{d}</td><td>{cnt}</td></tr>'
    html += '</table></details>'
    return section('GO Term Tree Depth / Granularity', html)


def build_enrichment_section(rows: list) -> str:
    html = '<table><tr><th>File</th><th>Lines</th><th>Size</th></tr>'
    for fname, lines, size in rows:
        html += f'<tr><td>{fname}</td><td>{lines}</td><td>{size}</td></tr>'
    html += '</table>'
    return section('GO Enrichment Input Files', html)


def build_predict_metrics_section(metrics: dict) -> str:
    if not metrics:
        return section('Prediction Performance', '<p>No metrics found in predict log.</p>')

    priority = ['Accuracy', 'AUC', 'Precision', 'Recall', 'Loss', 'F1', 'Threshold',
                'Training set', 'Test set']
    cards = []
    seen = set()
    for key in priority:
        if key in metrics:
            cards.append(metric_card(metrics[key], key))
            seen.add(key)
    for key, val in metrics.items():
        if key not in seen:
            cards.append(metric_card(val, key))

    html = '<div class="metric-grid">' + ''.join(cards) + '</div>'
    return section('Prediction Performance', html)


def build_images_section(image_paths: list[str], title: str = 'Plots') -> str:
    if not image_paths:
        return ''
    parts = [f'<h2>{title}</h2>']
    for img_path in image_paths:
        uri = embed_image(img_path)
        label = os.path.relpath(img_path)
        if uri:
            parts.append(f'<div><h3>{label}</h3><img src="{uri}" alt="{label}"></div>')
        else:
            parts.append(f'<div><h3>{label}</h3><p class="warn">Could not embed image.</p></div>')
    return '\n'.join(parts)


def build_log_section(log_path: str, title: str = 'Log (last 60 lines)') -> str:
    lines = _read(log_path).splitlines()
    tail = '\n'.join(lines[-60:])
    return section(title, f'<details><summary>Show log tail</summary><pre>{tail}</pre></details>')


# ── main report builder ───────────────────────────────────────────────────────

def generate_report(input_dir: str, out_path: str) -> None:
    input_dir = os.path.abspath(input_dir)
    logger.info(f"Generating HTML report for: {input_dir}")

    parts: list[str] = []
    title = os.path.basename(input_dir)

    parts.append(f'<html><head><meta charset="utf-8">'
                 f'<title>PhenGO Report – {title}</title>'
                 f'<style>{CSS}</style></head><body>')
    parts.append(f'<h1>PhenGO Report – {title}</h1>')

    entries = os.listdir(input_dir)

    # ── PhenGO_params.txt ─────────────────────────────────────────────────────
    params_path = os.path.join(input_dir, 'PhenGO_params.txt')
    if os.path.exists(params_path):
        params = parse_params(params_path)
        parts.append(build_params_section(params))
    else:
        logger.debug("PhenGO_params.txt not found")

    # ── ARFF file ─────────────────────────────────────────────────────────────
    arff_files = [f for f in entries if f.endswith('.arff')]
    if arff_files:
        arff_path = os.path.join(input_dir, arff_files[0])
        logger.info(f"  Parsing ARFF: {arff_path}")
        try:
            stats = parse_arff_stats(arff_path)
            parts.append(build_arff_section(stats))
        except Exception as e:
            logger.warning(f"  Could not parse ARFF: {e}")
            parts.append(section('ARFF Dataset Statistics',
                                  f'<p class="warn">Parse error: {e}</p>'))

    # ── GO term validation ───────────────────────────────────────────────────
    val_path = os.path.join(input_dir, 'GO_term_validation_report.txt')
    if os.path.exists(val_path):
        val = parse_go_validation(val_path)
        parts.append(build_go_validation_section(val))

    # ── Unique GO nodes ──────────────────────────────────────────────────────
    uniq_path = os.path.join(input_dir, 'Unique_GO_Nodes.txt')
    if os.path.exists(uniq_path):
        n_nodes = parse_unique_go_nodes(uniq_path)
        parts.append(section('Unique GO Nodes',
                              '<div class="metric-grid">' +
                              metric_card(f"{n_nodes:,}", 'Unique GO nodes used') +
                              '</div>'))

    # ── GO tree depth (GO_Children&Parents.json) ─────────────────────────────
    go_json_path = os.path.join(input_dir, 'GO_Children&Parents.json')
    if os.path.exists(go_json_path):
        depth_stats = go_depth_stats(go_json_path)
        if depth_stats:
            parts.append(build_go_depth_section(depth_stats))
        else:
            parts.append(section('GO Term Tree Depth',
                                  '<p>Could not extract depth information from '
                                  'GO_Children&amp;Parents.json '
                                  '(depth field not found in expected format).</p>'))

    # ── go_enrichment/ ───────────────────────────────────────────────────────
    enrich_dir = os.path.join(input_dir, 'go_enrichment')
    if os.path.isdir(enrich_dir):
        rows = enrichment_file_stats(enrich_dir)
        parts.append(build_enrichment_section(rows))

    # ── PhenGO main log ──────────────────────────────────────────────────────
    main_log = os.path.join(input_dir, 'PhenGO.log')
    if os.path.exists(main_log):
        parts.append(build_log_section(main_log, 'PhenGO Log (last 60 lines)'))

    # ── Predict subdirectory ─────────────────────────────────────────────────
    predict_dir = os.path.join(input_dir, 'Predict')
    if not os.path.isdir(predict_dir):
        # Maybe we're *inside* a Predict directory already
        predict_dir = input_dir

    predict_log = os.path.join(predict_dir, 'PhenGO_Predict.log')
    if os.path.exists(predict_log):
        parts.append('<h2>PhenGO-Predict Results</h2>')
        metrics = parse_predict_log(predict_log)
        parts.append(build_predict_metrics_section(metrics))

        # Per-model final reports
        for root, dirs, files in os.walk(predict_dir):
            for fname in sorted(files):
                if fname == 'final_report.txt':
                    rp = os.path.join(root, fname)
                    model_name = os.path.basename(root)
                    content = parse_final_report(rp)
                    parts.append(
                        section(f'Final Report – {model_name}',
                                f'<details><summary>Show full report</summary>'
                                f'<pre>{content}</pre></details>'))

        # Images (plots)
        imgs = find_images(predict_dir, recurse=True)
        if imgs:
            parts.append(build_images_section(imgs, 'Prediction Plots'))

        # Predict log tail
        parts.append(build_log_section(predict_log, 'Predict Log (last 60 lines)'))

    # ── Any remaining CSVs in the root ──────────────────────────────────────
    csvs = [f for f in entries if f.endswith('.csv')]
    if csvs:
        parts.append('<h2>Output Tables (preview)</h2>')
        try:
            import csv as _csv
            for csv_file in sorted(csvs):
                full = os.path.join(input_dir, csv_file)
                rows_html = []
                try:
                    with open(full, newline='', encoding='utf-8', errors='replace') as fh:
                        reader = _csv.reader(fh)
                        all_rows = list(reader)
                    if all_rows:
                        header = all_rows[0]
                        rows_html.append('<table>')
                        rows_html.append('<tr>' + ''.join(f'<th>{h}</th>' for h in header) + '</tr>')
                        for row in all_rows[1:21]:
                            rows_html.append('<tr>' + ''.join(f'<td>{c}</td>' for c in row) + '</tr>')
                        rows_html.append('</table>')
                        if len(all_rows) > 21:
                            rows_html.append(f'<p><em>… {len(all_rows)-1:,} rows total (showing first 20)</em></p>')
                except Exception as e:
                    rows_html.append(f'<p class="warn">Could not read {csv_file}: {e}</p>')
                parts.append(f'<h3>{csv_file}</h3>' + '\n'.join(rows_html))
        except ImportError:
            parts.append('<p>csv module unavailable.</p>')

    parts.append('</body></html>')

    html = '\n'.join(parts)
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(html)
    logger.info(f"HTML report written to: {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Generate a rich HTML summary report from a PhenGO results directory')
    parser.add_argument('-i', '--input-dir', required=True,
                        help='PhenGO output directory to summarise')
    parser.add_argument('-o', '--output',
                        help='Output HTML file path (default: <input-dir>/PhenGO_report.html)')
    args = parser.parse_args()

    configure_logger('PhenGO.report_generator', enable_file=False)

    input_dir = os.path.abspath(args.input_dir)
    out_path = args.output if args.output else os.path.join(input_dir, 'PhenGO_report.html')

    generate_report(input_dir, out_path)
    print(f"Report written to: {out_path}")


if __name__ == '__main__':
    main()

