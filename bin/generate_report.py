import base64
import pathlib
import argparse
import pandas as pd
from jinja2 import Template

HTML_TEMPLATE = """<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>RNA‑pair report</title>
  <!-- DataTables -->
  <link  href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css" rel="stylesheet">
  <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
  <script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
  <style>
      body{ font-family: sans-serif; margin:1rem;}
      table.dataTable tbody td{ vertical-align:top; }
      .svg-cell{ min-width:260px; }
      /* Monospace for sequence columns */
      {% for col in sequence_cols %}td.{{col|replace(' ', '_')}}{ font-family: monospace; }
      {% endfor %}
  </style>
</head>
<body>
  <h1>Top‑N structure pairs</h1>
  <table id="report" class="display">
    <thead>
      <tr>
        {% for col in cols %}<th>{{ col }}</th>{% endfor %}
        <th>Structure 1</th>
        <th>Structure 2</th>
      </tr>
    </thead>
    <tbody>
    {% for row in rows %}
      <tr>
        {% for col in cols %}
          <td class="{{col|replace(' ', '_')}}">
            {% if col in sequence_cols %}
              {{ row[col]|safe }}  {# Safe render for <br> tags #}
            {% else %}
              {{ row[col] }}
            {% endif %}
          </td>
        {% endfor %}
        <td class="svg-cell">{{ row["__svg1"]|safe }}</td>
        <td class="svg-cell">{{ row["__svg2"]|safe }}</td>
      </tr>
    {% endfor %}
    </tbody>
  </table>
<script>
$(document).ready(function(){
    $('#report').DataTable({
        pageLength: 20
        {% if score_col_index is not none %},
        order: [[ {{ score_col_index }}, 'desc' ]]
        {% endif %}
    });
});
</script>
</body>
</html>
"""

def inline_svg(path: pathlib.Path) -> str:
    try:
        text = path.read_text()
        if text.lstrip().startswith('<?xml'):
            text = text.split('?>',1)[1]
        b64 = base64.b64encode(text.encode('utf-8')).decode('ascii')
        return f'<img src="data:image/svg+xml;base64,{b64}" width="250"/>'
    except Exception:
        return f'<span style="color:red">Missing {path.name}</span>'


def make_report(pairs_tsv, svg_dir, output_html, id_column):
    df = pd.read_csv(pairs_tsv, sep='\t')
    svg_dir = pathlib.Path(svg_dir)
    
    # Dynamic ID columns
    id1_col = f"{id_column}_1"
    id2_col = f"{id_column}_2"

    # Auto-detect sequence columns (case-insensitive)
    sequence_cols = [col for col in df.columns if 'sequence' in col.lower()]
    
    # Determine score column index if exists
    score_col_index = None
    if 'score' in df.columns:
        score_col_index = df.columns.get_loc('score')
    
    rows = []
    for idx, rec in df.iterrows():
        i = idx + 1
        row_data = rec.to_dict()
        
        # Handle sequence columns
        for col in sequence_cols:
            val = row_data.get(col)
            if pd.isna(val):
                processed = ''
            else:
                seq = str(val).strip()
                chunks = [seq[j:j+30] for j in range(0, len(seq), 30)]
                processed = '<br>'.join(chunks)
            row_data[col] = processed
        
        # Extract dynamic IDs
        id1 = row_data.get(id1_col)
        id2 = row_data.get(id2_col)

        # Inlined SVGs
        safe1 = ''.join(c if c.isalnum() else '_' for c in f"{id1}_{i}")
        safe2 = ''.join(c if c.isalnum() else '_' for c in f"{id2}_{i}")
        row_data["__svg1"] = inline_svg(svg_dir / f"{safe1}.svg")
        row_data["__svg2"] = inline_svg(svg_dir / f"{safe2}.svg")
        
        rows.append(row_data)
    
    html = Template(HTML_TEMPLATE).render(
        cols=list(df.columns),
        rows=rows,
        score_col_index=score_col_index,
        sequence_cols=sequence_cols
    )
    pathlib.Path(output_html).write_text(html, encoding='utf-8')
    print(f"[ok] Wrote report with formatted sequences to {output_html}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate an HTML report of top‑N structure pairs."
    )
    parser.add_argument('--pairs',     required=True,
                        help='Top‑N TSV file containing pair data')
    parser.add_argument('--svg-dir',   required=True,
                        help='Directory with individual SVG files')
    parser.add_argument('--output',    required=True,
                        help='Output HTML path')
    parser.add_argument('--id-column', default='exon_id',
                        help='Base name of the original ID column (without _1/_2)')
    args = parser.parse_args()
    make_report(args.pairs, args.svg_dir, args.output, args.id_column)
