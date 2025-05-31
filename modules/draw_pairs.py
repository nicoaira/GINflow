#!/usr/bin/env python3
"""
draw_pairs.py – generate RNArtistCore plots for RNA pairs in parallel,
highlight a window, strip the Reactivity legend, stitch SVGs side-by-side,
and export PNGs.

Usage:
  python3 draw_pairs.py --tsv top_N.tsv --outdir outdir \
        --id-column <id_column> [--width W] [--height H] \
        [--highlight-colour COLOUR] [--num-workers N] \
        [--keep-temp] [--debug]
"""

import os
import csv
import subprocess
import tempfile
import shutil
import xml.etree.ElementTree as ET
import argparse
import sys
import re
import threading
import concurrent.futures

# ensure SVG namespace
ET.register_namespace('', 'http://www.w3.org/2000/svg')

# load CairoSVG
try:
    import cairosvg
except ImportError:
    sys.exit("[fatal] Please install 'cairosvg' (pip install cairosvg)")

# This will be set based on --id-column
ID_COLUMN = None

def parse_args():
    p = argparse.ArgumentParser(
        description='Draw & highlight RNA pairs – PNG output (parallel)')
    p.add_argument('--tsv',            required=True,  help='Input TSV')
    p.add_argument('--outdir',         required=True,  help='Output directory')
    p.add_argument('--id-column',      default='exon_id',
                       help='Name of the original ID column (without _1/_2).')
    p.add_argument('--width',    type=float, default=500,  help='SVG width')
    p.add_argument('--height',   type=float, default=500,  help='SVG height')
    p.add_argument('--highlight-colour', default="#FFD700",
                   help='HTML colour for window highlight')
    p.add_argument('--num-workers', type=int, default=1,
                   help='Parallel workers for rendering')
    p.add_argument('--keep-temp', action='store_true', help='Keep temp dir')
    p.add_argument('--debug',     action='store_true', help='Verbose')
    return p.parse_args()

def make_highlight(start, end, colour):
    if not start or not end:
        return ""
    return f'''        color {{
            value = "{colour}"
            location {{ {start} to {end} }}
        }}
'''

KTS_TEMPLATE = '''import io.github.fjossinet.rnartist.core.*

rnartist {{
    ss {{
        bn {{
            seq   = "{seq}"
            value = "{dot}"
            name  = "{name}"
        }}
    }}

    theme {{
        details {{ value = 4 }}
{highlight_block}    }}

    svg {{
        path   = "{path}"
        width  = {w}
        height = {h}
    }}
}}'''

def combine(svg1, svg2, dst_svg):
    try:
        t1, t2 = ET.parse(svg1), ET.parse(svg2)
        r1, r2 = t1.getroot(), t2.getroot()
        w1 = float(r1.get('width','0').rstrip('px'))
        h1 = float(r1.get('height','0').rstrip('px'))
        w2 = float(r2.get('width','0').rstrip('px'))
        h2 = float(r2.get('height','0').rstrip('px'))
        W, H = w1 + w2, max(h1, h2)

        svg = ET.Element('{http://www.w3.org/2000/svg}svg',
                         {'width':str(W),'height':str(H)})
        g1 = ET.SubElement(svg,'g')
        for el in list(r1): g1.append(el)
        g2 = ET.SubElement(svg,'g',{'transform':f'translate({w1},0)'})
        for el in list(r2): g2.append(el)
        ET.ElementTree(svg).write(dst_svg, xml_declaration=True, encoding='utf-8')
        return True,(W,H)
    except Exception as e:
        print(f"[combine] {e}", file=sys.stderr)
        return False,(0,0)

def svg_to_png(src, dst):
    try:
        cairosvg.svg2png(url=src, write_to=dst)
        return True
    except Exception as e:
        print(f"[png] {e}", file=sys.stderr)
        return False

SVG_REACTIVITY_PATTERN = re.compile(
    r'<defs>.*?id="reactivities_scale".*?</text>\s*</svg>',
    flags=re.DOTALL
)
def strip_reactivity_legend(svg_path):
    try:
        txt = open(svg_path,'r',encoding='utf-8').read()
        cleaned, n = SVG_REACTIVITY_PATTERN.subn('</svg>', txt, count=1)
        if n:
            open(svg_path,'w',encoding='utf-8').write(cleaned)
    except Exception as e:
        print(f"[strip] {svg_path}: {e}", file=sys.stderr)

def validate(seq, dot):
    return len(seq)==len(dot)

def process_pair(i, row):
    # Dynamically pick up ID_COLUMN_1 and ID_COLUMN_2
    id1 = row.get(f'{ID_COLUMN}_1') or row.get('id1') or f"RNA_{i}_1"
    id2 = row.get(f'{ID_COLUMN}_2') or row.get('id2') or f"RNA_{i}_2"
    seq1 = row.get('sequence_1') or row.get('seq1') or ""
    seq2 = row.get('sequence_2') or row.get('seq2') or ""
    dot1 = row.get('secondary_structure_1') or row.get('structure1') or "."*len(seq1)
    dot2 = row.get('secondary_structure_2') or row.get('structure2') or "."*len(seq2)

    seq1 = seq1.replace('T','U').replace('t','u')
    seq2 = seq2.replace('T','U').replace('t','u')

    if not (validate(seq1,dot1) and validate(seq2,dot2)):
        return

    w1s = int(row.get('window_start_1',0) or 0)
    w1e = int(row.get('window_end_1',  0) or 0)
    w2s = int(row.get('window_start_2',0) or 0)
    w2e = int(row.get('window_end_2',  0) or 0)
    hl1 = make_highlight(w1s,w1e,ARGS.highlight_colour)
    hl2 = make_highlight(w2s,w2e,ARGS.highlight_colour)

    n1 = re.sub(r'[^A-Za-z0-9_]', '_', f"{id1}_{i}")
    n2 = re.sub(r'[^A-Za-z0-9_]', '_', f"{id2}_{i}")

    w_str = f"{ARGS.width:.1f}"
    h_str = f"{ARGS.height:.1f}"

    s1 = KTS_TEMPLATE.format(
        path=TMPDIR, w=w_str, h=h_str,
        seq=seq1, dot=dot1, name=n1,
        highlight_block=hl1
    )
    s2 = KTS_TEMPLATE.format(
        path=TMPDIR, w=w_str, h=h_str,
        seq=seq2, dot=dot2, name=n2,
        highlight_block=hl2
    )

    for name,script in ((n1,s1),(n2,s2)):
        with open(os.path.join(KTS_DIR,f"{name}.kts"),'w') as fh: fh.write(script)
        with open(os.path.join(TMPDIR,f"{name}.kts"),'w') as fh: fh.write(script)

    if ARGS.debug: print(f"[rnartist] {n1}")
    subprocess.run(['rnartistcore', os.path.join(TMPDIR,f"{n1}.kts")])
    if ARGS.debug: print(f"[rnartist] {n2}")
    subprocess.run(['rnartistcore', os.path.join(TMPDIR,f"{n2}.kts")])

    svg1 = os.path.join(TMPDIR, f"{n1}.svg")
    svg2 = os.path.join(TMPDIR, f"{n2}.svg")
    if not (os.path.exists(svg1) and os.path.exists(svg2)):
        with LOG_LOCK: LOG_FH.write(f"{n1}\n{n2}\n")
        return

    strip_reactivity_legend(svg1)
    strip_reactivity_legend(svg2)

    indiv = os.path.join(ARGS.outdir,"individual_svgs"); os.makedirs(indiv, exist_ok=True)
    shutil.copy(svg1, os.path.join(indiv,f"{n1}.svg"))
    shutil.copy(svg2, os.path.join(indiv,f"{n2}.svg"))

    base    = f"pair_{i}_{id1}_{id2}"
    pair_svg= os.path.join(ARGS.outdir, base + ".svg")
    ok,_    = combine(svg1,svg2,pair_svg)
    if ok: strip_reactivity_legend(pair_svg)

    pair_png = os.path.join(ARGS.outdir, base + ".png")
    if svg_to_png(pair_svg, pair_png):
        print(f"[ok] {pair_png}")
    else:
        with LOG_LOCK: LOG_FH.write(f"{n1}\n{n2}\n")

if __name__ == '__main__':
    ARGS   = parse_args()
    global ID_COLUMN
    ID_COLUMN = ARGS.id_column

    # prepare output
    os.makedirs(ARGS.outdir, exist_ok=True)
    KTS_DIR  = os.path.join(ARGS.outdir,"kts_scripts"); os.makedirs(KTS_DIR, exist_ok=True)
    LOG_FH   = open(os.path.join(ARGS.outdir,'failures.log'),'w')
    LOG_LOCK = threading.Lock()
    TMPDIR   = tempfile.mkdtemp(prefix='rn_')
    if not ARGS.keep_temp:
        def cleanup():
            shutil.rmtree(TMPDIR, ignore_errors=True)
        import atexit
        atexit.register(cleanup)

    # load all rows
    with open(ARGS.tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows   = list(reader)

    # dispatch
    if ARGS.num_workers>1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=ARGS.num_workers) as exe:
            futures = [
                exe.submit(process_pair, i, row)
                for i,row in enumerate(rows,1)
            ]
            for _ in concurrent.futures.as_completed(futures):
                pass
    else:
        for i,row in enumerate(rows,1):
            process_pair(i,row)

    LOG_FH.close()
    if ARGS.keep_temp and ARGS.debug:
        print(f"[info] kept temp dir {TMPDIR}")
