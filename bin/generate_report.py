#!/usr/bin/env python3
"""
HTML Report Generator for BaitBench

Generates a simple HTML report summarizing capture results.
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

# HTML template
HTML_TEMPLATE = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BaitBench Report - {run_name}</title>
    <style>
        * {{
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            max-width: 1000px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .metric-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            text-align: center;
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #2c3e50;
        }}
        .metric-label {{
            color: #7f8c8d;
            font-size: 0.9em;
            text-transform: uppercase;
        }}
        .metric-good {{
            color: #27ae60;
        }}
        .metric-warn {{
            color: #f39c12;
        }}
        .metric-bad {{
            color: #e74c3c;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #ecf0f1;
        }}
        th {{
            background: #3498db;
            color: white;
            font-weight: 600;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .classification-TP {{
            color: #27ae60;
            font-weight: bold;
        }}
        .classification-FP {{
            color: #e74c3c;
            font-weight: bold;
        }}
        .classification-FN {{
            color: #f39c12;
            font-weight: bold;
        }}
        .footer {{
            text-align: center;
            color: #95a5a6;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
        }}
        .info-box {{
            background: #e8f4f8;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 20px 0;
            border-radius: 0 8px 8px 0;
        }}
    </style>
</head>
<body>
    <h1>BaitBench Report</h1>

    <div class="info-box">
        <strong>Run:</strong> {run_name}<br>
        <strong>Generated:</strong> {timestamp}
    </div>

    <h2>Summary Metrics</h2>
    <div class="summary-grid">
        <div class="metric-card">
            <div class="metric-value">{reads_generated:,}</div>
            <div class="metric-label">Reads Generated</div>
        </div>
        <div class="metric-card">
            <div class="metric-value">{reads_captured:,}</div>
            <div class="metric-label">Reads Captured</div>
        </div>
        <div class="metric-card">
            <div class="metric-value">{capture_rate_pct:.1f}%</div>
            <div class="metric-label">Capture Rate</div>
        </div>
    </div>

    <h2>Captured Reads by Source</h2>
    <div class="summary-grid">
        <div class="metric-card">
            <div class="metric-value metric-good">{target_captured:,}</div>
            <div class="metric-label">Target Reads Captured</div>
        </div>
        <div class="metric-card">
            <div class="metric-value metric-bad">{distractor_captured:,}</div>
            <div class="metric-label">Distractor Reads Captured</div>
        </div>
    </div>

    <h2>Detection Performance</h2>
    <div class="summary-grid">
        <div class="metric-card">
            <div class="metric-value {sensitivity_class}">{sensitivity_pct:.1f}%</div>
            <div class="metric-label">Sensitivity</div>
        </div>
        <div class="metric-card">
            <div class="metric-value {specificity_class}">{specificity_pct:.1f}%</div>
            <div class="metric-label">Specificity</div>
        </div>
        <div class="metric-card">
            <div class="metric-value {precision_class}">{precision_pct:.1f}%</div>
            <div class="metric-label">Precision</div>
        </div>
        <div class="metric-card">
            <div class="metric-value {f1_class}">{f1_pct:.1f}%</div>
            <div class="metric-label">F1 Score</div>
        </div>
    </div>

    <h2>Confusion Matrix</h2>
    <div class="summary-grid">
        <div class="metric-card">
            <div class="metric-value metric-good">{tp_count}</div>
            <div class="metric-label">True Positives</div>
        </div>
        <div class="metric-card">
            <div class="metric-value metric-bad">{fp_count}</div>
            <div class="metric-label">False Positives</div>
        </div>
        <div class="metric-card">
            <div class="metric-value metric-warn">{fn_count}</div>
            <div class="metric-label">False Negatives</div>
        </div>
        <div class="metric-card">
            <div class="metric-value metric-good">{tn_count}</div>
            <div class="metric-label">True Negatives</div>
        </div>
    </div>

    <h2>Detection Details</h2>
    {detail_table}

    <div class="footer">
        Generated by <a href="https://github.com/yourusername/BaitBench">BaitBench</a>
    </div>
</body>
</html>
'''


def parse_summary_tsv(summary_file):
    """Parse summary TSV file."""
    with open(summary_file, 'r') as f:
        headers = f.readline().strip().split('\t')
        values = f.readline().strip().split('\t')

    return dict(zip(headers, values))


def parse_detail_tsv(detail_file):
    """Parse detail TSV file."""
    rows = []
    with open(detail_file, 'r') as f:
        headers = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            rows.append(dict(zip(headers, values)))
    return rows


def get_metric_class(value, thresholds=(0.9, 0.7)):
    """Return CSS class based on metric value."""
    if value >= thresholds[0]:
        return 'metric-good'
    elif value >= thresholds[1]:
        return 'metric-warn'
    else:
        return 'metric-bad'


def generate_detail_table(details):
    """Generate HTML table for detection details."""
    if not details:
        return '<p>No detection details available.</p>'

    rows_html = []
    for row in details:
        classification = row.get('classification', '')
        css_class = f'classification-{classification}'

        rows_html.append(f'''
        <tr>
            <td>{row.get('reference_id', '')}</td>
            <td>{row.get('category', '')}</td>
            <td>{row.get('read_count', '')}</td>
            <td class="{css_class}">{classification}</td>
        </tr>
        ''')

    return f'''
    <table>
        <thead>
            <tr>
                <th>Reference ID</th>
                <th>Category</th>
                <th>Read Count</th>
                <th>Classification</th>
            </tr>
        </thead>
        <tbody>
            {''.join(rows_html)}
        </tbody>
    </table>
    '''


def main():
    parser = argparse.ArgumentParser(
        description='Generate HTML report for BaitBench results'
    )

    parser.add_argument('--summary', required=True,
                        help='Summary TSV file (results.tsv)')
    parser.add_argument('--detail', required=True,
                        help='Detail TSV file (detected_detail.tsv)')
    parser.add_argument('--output', required=True,
                        help='Output HTML file')
    parser.add_argument('--run-name', default='BaitBench Run',
                        help='Run name for report title')

    args = parser.parse_args()

    # Parse input files
    print(f"Reading summary from {args.summary}...", file=sys.stderr)
    summary = parse_summary_tsv(args.summary)

    print(f"Reading details from {args.detail}...", file=sys.stderr)
    details = parse_detail_tsv(args.detail)

    # Extract and convert values
    reads_generated = int(summary.get('reads_generated', 0))
    reads_captured = int(summary.get('reads_captured', 0))
    capture_rate = float(summary.get('capture_rate', 0))
    target_captured = int(summary.get('target_captured', 0))
    distractor_captured = int(summary.get('distractor_captured', 0))

    sensitivity = float(summary.get('sensitivity', 0))
    specificity = float(summary.get('specificity', 0))
    precision = float(summary.get('precision', 0))
    f1_score = float(summary.get('f1_score', 0))

    tp_count = int(summary.get('tp_count', 0))
    fp_count = int(summary.get('fp_count', 0))
    fn_count = int(summary.get('fn_count', 0))
    tn_count = int(summary.get('tn_count', 0))

    # Generate HTML
    html = HTML_TEMPLATE.format(
        run_name=args.run_name,
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        reads_generated=reads_generated,
        reads_captured=reads_captured,
        capture_rate_pct=capture_rate * 100,
        target_captured=target_captured,
        distractor_captured=distractor_captured,
        sensitivity_pct=sensitivity * 100,
        sensitivity_class=get_metric_class(sensitivity),
        specificity_pct=specificity * 100,
        specificity_class=get_metric_class(specificity),
        precision_pct=precision * 100,
        precision_class=get_metric_class(precision),
        f1_pct=f1_score * 100,
        f1_class=get_metric_class(f1_score),
        tp_count=tp_count,
        fp_count=fp_count,
        fn_count=fn_count,
        tn_count=tn_count,
        detail_table=generate_detail_table(details)
    )

    # Write output
    print(f"Writing report to {args.output}...", file=sys.stderr)
    with open(args.output, 'w') as f:
        f.write(html)

    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
