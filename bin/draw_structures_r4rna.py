#!/usr/bin/env python3
"""
draw_structures_r4rna.py – generate R4RNA arc plots for RNA pairs in parallel,
highlighting specific regions.

Usage:
  python3 draw_structures_r4rna.py --tsv top_N.tsv --outdir outdir \
        --id-column <id_column> [--highlight-colour COLOUR] [--num-workers N] [--debug]
"""

import os
import csv
import subprocess
import tempfile
import shutil
import argparse
import sys
import re
import threading
import concurrent.futures

# This will be set based on --id-column
ID_COLUMN = None

def parse_args():
    p = argparse.ArgumentParser(
        description='Draw & highlight RNA pairs using R4RNA – SVG output (parallel)')
    p.add_argument('--tsv',            required=True,  help='Input TSV')
    p.add_argument('--outdir',         required=True,  help='Output directory')
    p.add_argument('--id-column',      default='exon_id',
                       help='Name of the original ID column (without _1/_2).')
    p.add_argument('--pair-type',      choices=['window', 'contig', 'alignment'], default='window',
                       help='Type of pairs in TSV: "window", "contig", or "alignment"')
    p.add_argument('--highlight-colour', default="#FF0000",
                   help='HTML colour for region highlight (default: red)')
    p.add_argument('--num-workers', type=int, default=1,
                   help='Parallel workers for rendering')
    p.add_argument('--keep-temp', action='store_true', help='Keep temp dir')
    p.add_argument('--debug',     action='store_true', help='Verbose')
    # R4RNA specific parameters
    p.add_argument('--width-per-100nt', type=float, default=2.0,
                   help='Plot width per 100 nucleotides')
    p.add_argument('--height', type=float, default=4.0,
                   help='Plot height in inches')
    p.add_argument('--arrowhead-size', type=float, default=0.08,
                   help='Arrowhead size in inches')
    return p.parse_args()

def _parse_int_list(v):
    """Parse an integer or a comma/semicolon-separated list of integers.
    Returns a list of ints. Empty / None -> [].
    """
    if v is None:
        return []
    if isinstance(v, (int, float)):
        try:
            iv = int(v)
            return [iv]
        except Exception:
            return []
    s = str(v).strip()
    if s == "":
        return []
    s = s.replace(";", ",").replace("|", ",")
    parts = [p.strip() for p in s.split(",") if p.strip() != ""]
    out = []
    for p in parts:
        try:
            out.append(int(p))
        except Exception:
            continue
    return out

def validate(seq, dot):
    return len(seq) == len(dot)

R_TEMPLATE = '''# Load the R4RNA library
library(R4RNA)

# Data for structure
seq <- "{seq}"
ss <- "{ss}"

# Scaling parameters
len <- nchar(seq)
width_per_100_nt <- {width_per_100nt}
plot_width <- (len / 100) * width_per_100_nt
plot_height <- {plot_height}
arrowhead_size <- {arrowhead_size}
line_width <- 1.5

# Export SVG
svg("{output_svg}", width = plot_width, height = plot_height)

helix_df <- viennaToHelix(ss)
helix_df$col <- "grey"

# Highlight regions
{highlight_code}

# Plot arcs without default line and arrow
plotHelix(helix_df, line = FALSE, arrow = FALSE, lwd = line_width)

# Draw baseline with highlighted segments
{baseline_code}

dev.off()
'''

def make_r_highlight(starts_list, ends_list):
    """Generate R code to highlight multiple regions"""
    if not starts_list or not ends_list:
        return ""
    
    highlight_lines = []
    n = min(len(starts_list), len(ends_list))
    for i in range(n):
        s, e = starts_list[i], ends_list[i]
        if s is not None and e is not None:
            # R uses 1-based indexing, so add 1 to 0-based coordinates
            highlight_lines.append(
                f'highlight_indices <- which(helix_df$i >= {s+1} & helix_df$j <= {e+1})\n'
                f'helix_df$col[highlight_indices] <- "{ARGS.highlight_colour}"'
            )
    
    return "\n".join(highlight_lines)

def make_baseline_code(starts_list, ends_list):
    """Generate R code to draw the baseline with highlighted segments"""
    if not starts_list or not ends_list:
        # No highlights - draw simple arrow from start to end
        return '''arrows(
  x0 = 1, y0 = 0, x1 = len, y1 = 0,
  length = arrowhead_size,
  code = 2,
  col = "black",
  lwd = line_width,
  xpd = TRUE
)'''
    
    # Build segments for the baseline
    n = min(len(starts_list), len(ends_list))
    segments = []
    
    # Get all breakpoints (start and end of each highlight range)
    breakpoints = [(1, starts_list[0] + 1)]  # +1 for R's 1-based indexing
    
    for i in range(n):
        s, e = starts_list[i], ends_list[i]
        if s is not None and e is not None:
            s_r = s + 1  # Convert to 1-based
            e_r = e + 1  # Convert to 1-based
            
            # Add highlighted segment
            breakpoints.append((s_r, e_r))
            
            # If there's another highlight after this one, add the gap
            if i + 1 < n and ends_list[i] is not None and starts_list[i + 1] is not None:
                breakpoints.append((e_r, starts_list[i + 1] + 1))
    
    # Add final segment from last highlight to end
    if ends_list[-1] is not None:
        breakpoints.append((ends_list[-1] + 1, None))  # None means "len"
    
    # Generate R code for each segment
    code_lines = []
    
    # First segment (before first highlight) - black
    first_highlight_start = starts_list[0] + 1 if starts_list else None
    if first_highlight_start and first_highlight_start > 1:
        code_lines.append(f'segments(x0 = 1, y0 = 0, x1 = {first_highlight_start}, y1 = 0, col = "black", lwd = line_width)')
    
    # Highlighted segments and gaps between them
    for i in range(n):
        s, e = starts_list[i], ends_list[i]
        if s is not None and e is not None:
            s_r = s + 1
            e_r = e + 1
            
            # Highlighted segment - red
            code_lines.append(f'segments(x0 = {s_r}, y0 = 0, x1 = {e_r}, y1 = 0, col = "{ARGS.highlight_colour}", lwd = line_width)')
            
            # Gap to next highlight (if exists) - black
            if i + 1 < n and starts_list[i + 1] is not None:
                next_start = starts_list[i + 1] + 1
                code_lines.append(f'segments(x0 = {e_r}, y0 = 0, x1 = {next_start}, y1 = 0, col = "black", lwd = line_width)')
    
    # Final segment (after last highlight) with arrow - black
    last_highlight_end = ends_list[-1] + 1 if ends_list else 1
    code_lines.append(f'''arrows(
  x0 = {last_highlight_end}, y0 = 0, x1 = len, y1 = 0,
  length = arrowhead_size,
  code = 2,
  col = "black",
  lwd = line_width,
  xpd = TRUE
)''')
    
    return '\n'.join(code_lines)

def process_pair(i, row):
    global ID_COLUMN
    
    # Get IDs
    id1 = row.get('query_id') or row.get(f'{ID_COLUMN}_1') or row.get('id1') or f"RNA_{i}_1"
    id2 = row.get('target_id') or row.get('subject_id') or row.get(f'{ID_COLUMN}_2') or row.get('id2') or f"RNA_{i}_2"

    # Get sequences and structures
    seq1 = row.get('query_sequence') or row.get('sequence_1') or row.get('seq1') or ""
    seq2 = row.get('target_sequence') or row.get('subject_sequence') or row.get('sequence_2') or row.get('seq2') or ""
    
    dot1 = (
        row.get('query_secondary_structure') or
        row.get('query_structure') or
        row.get('query_rnafold_dotbracket') or
        row.get('secondary_structure_1') or
        row.get('structure_1') or
        row.get('structure1') or "."*len(seq1)
    )
    dot2 = (
        row.get('target_structure') or
        row.get('subject_secondary_structure') or
        row.get('subject_structure') or
        row.get('subject_rnafold_dotbracket') or
        row.get('secondary_structure_2') or
        row.get('structure_2') or
        row.get('structure2') or "."*len(seq2)
    )

    # Convert T to U
    seq1 = seq1.replace('T','U').replace('t','u')
    seq2 = seq2.replace('T','U').replace('t','u')

    # Validate
    if not (validate(seq1, dot1) and validate(seq2, dot2)):
        if ARGS.debug:
            print(f"[warning] Validation failed for row {i}: id1={id1}, id2={id2}", file=sys.stderr)
        return

    # Get highlight coordinates based on pair type
    if ARGS.pair_type == 'contig':
        w1s_list = _parse_int_list(row.get('query_contig_start') or row.get('contig_start_1') or None)
        w1e_list = _parse_int_list(row.get('query_contig_end') or row.get('contig_end_1') or None)
        w2s_list = _parse_int_list(row.get('subject_contig_start') or row.get('contig_start_2') or None)
        w2e_list = _parse_int_list(row.get('subject_contig_end') or row.get('contig_end_2') or None)
    elif ARGS.pair_type == 'alignment':
        w1s_list = _parse_int_list(row.get('query_start') or None)
        w1e_list = _parse_int_list(row.get('query_end') or None)
        w2s_list = _parse_int_list(row.get('target_start') or None)
        w2e_list = _parse_int_list(row.get('target_end') or None)
    else:  # window
        w1s_list = _parse_int_list(row.get('query_window_start') or row.get('window_start_1') or None)
        w1e_list = _parse_int_list(row.get('query_window_end') or row.get('window_end_1') or None)
        w2s_list = _parse_int_list(row.get('subject_window_start') or row.get('window_start_2') or None)
        w2e_list = _parse_int_list(row.get('subject_window_end') or row.get('window_end_2') or None)
    
    if ARGS.debug:
        print(f"[debug] Row {i}: {ARGS.pair_type}_1 starts={w1s_list} ends={w1e_list}; {ARGS.pair_type}_2 starts={w2s_list} ends={w2e_list}")
    
    # Create sanitized names
    n1 = re.sub(r'[^A-Za-z0-9_]', '_', f"{id1}_{i}")
    n2 = re.sub(r'[^A-Za-z0-9_]', '_', f"{id2}_{i}")

    # Generate R scripts
    hl1 = make_r_highlight(w1s_list, w1e_list)
    hl2 = make_r_highlight(w2s_list, w2e_list)
    
    # Generate baseline code with highlighted segments
    baseline1 = make_baseline_code(w1s_list, w1e_list)
    baseline2 = make_baseline_code(w2s_list, w2e_list)

    svg1 = os.path.join(TMPDIR, f"{n1}.svg")
    svg2 = os.path.join(TMPDIR, f"{n2}.svg")

    r_script1 = R_TEMPLATE.format(
        seq=seq1,
        ss=dot1,
        output_svg=svg1,
        width_per_100nt=ARGS.width_per_100nt,
        plot_height=ARGS.height,
        arrowhead_size=ARGS.arrowhead_size,
        highlight_code=hl1 if hl1 else "# No highlights",
        baseline_code=baseline1
    )
    
    r_script2 = R_TEMPLATE.format(
        seq=seq2,
        ss=dot2,
        output_svg=svg2,
        width_per_100nt=ARGS.width_per_100nt,
        plot_height=ARGS.height,
        arrowhead_size=ARGS.arrowhead_size,
        highlight_code=hl2 if hl2 else "# No highlights",
        baseline_code=baseline2
    )

    # Write R scripts
    r_file1 = os.path.join(R_SCRIPTS_DIR, f"{n1}.R")
    r_file2 = os.path.join(R_SCRIPTS_DIR, f"{n2}.R")
    
    with open(r_file1, 'w') as fh:
        fh.write(r_script1)
    with open(r_file2, 'w') as fh:
        fh.write(r_script2)

    # Execute R scripts
    if ARGS.debug:
        print(f"[R4RNA] {n1}")
    
    result1 = subprocess.run(['Rscript', r_file1], capture_output=True, text=True)
    if result1.returncode != 0:
        print(f"[error] R4RNA failed for {n1}: {result1.stderr}", file=sys.stderr)
        with LOG_LOCK:
            LOG_FH.write(f"{n1}\tR4RNA error: {result1.stderr}\n")
        return
    
    if ARGS.debug:
        print(f"[R4RNA] {n2}")
    
    result2 = subprocess.run(['Rscript', r_file2], capture_output=True, text=True)
    if result2.returncode != 0:
        print(f"[error] R4RNA failed for {n2}: {result2.stderr}", file=sys.stderr)
        with LOG_LOCK:
            LOG_FH.write(f"{n2}\tR4RNA error: {result2.stderr}\n")
        return

    # Copy SVGs to output directory
    if os.path.exists(svg1) and os.path.exists(svg2):
        os.makedirs(ARGS.outdir, exist_ok=True)
        shutil.copy(svg1, os.path.join(ARGS.outdir, f"{n1}.svg"))
        shutil.copy(svg2, os.path.join(ARGS.outdir, f"{n2}.svg"))
    else:
        if ARGS.debug:
            print(f"[warning] SVG files not found: {svg1}={os.path.exists(svg1)}, {svg2}={os.path.exists(svg2)}", file=sys.stderr)
        with LOG_LOCK:
            LOG_FH.write(f"{n1}\n{n2}\n")

if __name__ == '__main__':
    ARGS = parse_args()
    ID_COLUMN = ARGS.id_column

    # Prepare output directories
    os.makedirs(ARGS.outdir, exist_ok=True)
    R_SCRIPTS_DIR = os.path.join(ARGS.outdir, "r_scripts")
    os.makedirs(R_SCRIPTS_DIR, exist_ok=True)
    
    LOG_FH = open(os.path.join(ARGS.outdir, 'failures.log'), 'w')
    LOG_LOCK = threading.Lock()
    
    TMPDIR = tempfile.mkdtemp(prefix='r4rna_')
    if not ARGS.keep_temp:
        def cleanup():
            shutil.rmtree(TMPDIR, ignore_errors=True)
        import atexit
        atexit.register(cleanup)

    # Load all rows
    with open(ARGS.tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        rows = list(reader)

    # Process pairs in parallel
    if ARGS.num_workers > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=ARGS.num_workers) as exe:
            futures = [exe.submit(process_pair, i, row) for i, row in enumerate(rows, 1)]
            for _ in concurrent.futures.as_completed(futures):
                pass
    else:
        for i, row in enumerate(rows, 1):
            process_pair(i, row)
    
    LOG_FH.close()
    print(f"[R4RNA] Processing complete. Output in {ARGS.outdir}")
