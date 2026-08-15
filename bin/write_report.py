#!/usr/bin/env python3
"""Build a self-contained HTML search report from ginflow outputs."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from datetime import datetime, timezone
from html import escape
from pathlib import Path


GRAY = "#B0B0B0"
SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")
CLUSTER_HEAD = re.compile(
    r"^# cluster\s+(\S+)\s+(\S+)\s+vs\s+(\S+)",
)
IDENTITY = re.compile(
    r"base identity\s*=\s*([0-9.]+)%.*structure identity\s*=\s*([0-9.]+)%",
    re.I,
)


def safe_name(value: str) -> str:
    cleaned = SAFE_NAME.sub("_", value).strip("_")
    return cleaned or "rna"


def parse_float(value: str | None, default: float = 0.0) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def parse_int(value: str | None, default: int = 0) -> int:
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


def fmt_evalue(value: float) -> str:
    if math.isinf(value):
        return "inf"
    if value == 0.0:
        return "0"
    if value >= 0.001 and value < 1000:
        text = f"{value:.4g}"
        return text
    return f"{value:.3e}".replace("e-0", "e-").replace("e+0", "e+")


def log10_e(value: float) -> float:
    if value <= 0.0:
        return 300.0
    if math.isinf(value):
        return 0.0
    return max(0.0, -math.log10(value))


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_alignment_text(path: Path | None) -> dict[tuple[str, str, str], dict]:
    if path is None or not path.is_file():
        return {}
    text = path.read_text()
    blocks: dict[tuple[str, str, str], dict] = {}
    current: list[str] = []
    key: tuple[str, str, str] | None = None
    meta: dict = {}
    for line in text.splitlines():
        match = CLUSTER_HEAD.match(line)
        if match:
            if key is not None:
                blocks[key] = {**meta, "text": "\n".join(current).rstrip()}
            key = (match.group(1), match.group(2), match.group(3))
            current = [line]
            meta = {}
            continue
        if key is None:
            continue
        current.append(line)
        ident = IDENTITY.search(line)
        if ident:
            meta["base_identity"] = float(ident.group(1))
            meta["structure_identity"] = float(ident.group(2))
    if key is not None:
        blocks[key] = {**meta, "text": "\n".join(current).rstrip()}
    return blocks


def load_svg_file(path: Path) -> str:
    raw = path.read_text(errors="replace")
    raw = re.sub(r"<script[\s\S]*?</script>", "", raw, flags=re.I)
    raw = re.sub(r"<\?xml[^>]*\?>", "", raw)
    return raw.strip()


def load_svgs(sources: Path | list[Path] | None) -> dict[str, str]:
    if sources is None:
        return {}
    paths = sources if isinstance(sources, list) else [sources]
    found: dict[str, str] = {}
    for source in paths:
        if source is None:
            continue
        if source.is_file() and source.suffix.lower() == ".svg":
            found[source.name] = load_svg_file(source)
            continue
        if not source.is_dir():
            continue
        for svg in sorted(source.rglob("*.svg")):
            found[svg.name] = load_svg_file(svg)
    return found


def match_plot(
    svgs: dict[str, str],
    index: int,
    query_id: str,
    target_id: str,
    kind: str,
    cluster_id: str = "",
) -> str | None:
    safe_q = safe_name(query_id)
    safe_t = safe_name(target_id)
    names = []
    if cluster_id != "":
        names.append(f"{safe_name(cluster_id)}_{safe_q}__{safe_t}_{kind}.svg")
    names.append(f"{index:04d}_{safe_q}__{safe_t}_{kind}.svg")
    for exact in names:
        if exact in svgs:
            return svgs[exact]
    suffix = f"_{kind}.svg"
    for name, body in svgs.items():
        if name.endswith(suffix) and safe_q in name and safe_t in name:
            return body
    return None


def fmt_pct(value: float | None) -> str:
    if value is None:
        return "—"
    return f"{value:.1f}%"


def e_class(evalue: float) -> str:
    if evalue < 1e-10:
        return "strong"
    if evalue < 1e-3:
        return "good"
    if evalue < 1e-2:
        return "weak"
    return "ns"


def build_hits(
    rows: list[dict[str, str]],
    texts: dict[tuple[str, str, str], dict],
    rn_svgs: dict[str, str],
    r4_svgs: dict[str, str],
    sw_svgs: dict[str, str] | None = None,
) -> list[dict]:
    hits = []
    for index, row in enumerate(rows):
        query_id = row.get("query_id", f"query_{index}")
        target_id = row.get("target_id", f"target_{index}")
        cluster_id = row.get("cluster_id", str(index))
        evalue = parse_float(row.get("evalue"), math.inf)
        q_len = parse_int(row.get("query_length"))
        t_len = parse_int(row.get("target_length"))
        q_start = parse_int(row.get("query_start"))
        q_end = parse_int(row.get("query_end"))
        t_start = parse_int(row.get("target_start"))
        t_end = parse_int(row.get("target_end"))
        aligned = parse_int(row.get("aligned_columns"))
        matches = parse_int(row.get("match_count"))
        block = texts.get((cluster_id, query_id, target_id), {})
        cover_q = (q_end - q_start) / q_len if q_len else 0.0
        cover_t = (t_end - t_start) / t_len if t_len else 0.0
        pair_id = matches / aligned if aligned else 0.0
        hits.append({
            "index": index,
            "rank": index + 1,
            "cluster_id": cluster_id,
            "query_id": query_id,
            "target_id": target_id,
            "self": query_id == target_id,
            "score": parse_float(row.get("score")),
            "bit_score": parse_float(row.get("bit_score")),
            "evalue": evalue,
            "evalue_pair": parse_float(row.get("evalue_pair"), math.inf),
            "evalue_label": fmt_evalue(evalue),
            "eclass": e_class(evalue),
            "neglog10e": log10_e(evalue),
            "query_start": q_start,
            "query_end": q_end,
            "target_start": t_start,
            "target_end": t_end,
            "query_length": q_len,
            "target_length": t_len,
            "match_count": matches,
            "aligned_columns": aligned,
            "seed_count": parse_int(row.get("seed_count")),
            "coverage_query": cover_q,
            "coverage_target": cover_t,
            "pair_identity": pair_id,
            "base_identity": block.get("base_identity"),
            "structure_identity": block.get("structure_identity"),
            "query_sequence": row.get("query_sequence", ""),
            "query_structure": row.get("query_structure", ""),
            "target_sequence": row.get("target_sequence", ""),
            "target_structure": row.get("target_structure", ""),
            "alignment_text": block.get("text", ""),
            "plot_rn_query": match_plot(rn_svgs, index, query_id, target_id, "query", cluster_id),
            "plot_rn_target": match_plot(rn_svgs, index, query_id, target_id, "target", cluster_id),
            "plot_r4_query": match_plot(r4_svgs, index, query_id, target_id, "query", cluster_id),
            "plot_r4_target": match_plot(r4_svgs, index, query_id, target_id, "target", cluster_id),
            "plot_sw_similarity": match_plot(sw_svgs or {}, index, query_id, target_id, "similarity", cluster_id),
            "plot_sw_scores": match_plot(sw_svgs or {}, index, query_id, target_id, "scores", cluster_id),
        })
    return hits


def query_summaries(hits: list[dict]) -> list[dict]:
    grouped: dict[str, list[dict]] = {}
    for hit in hits:
        grouped.setdefault(hit["query_id"], []).append(hit)
    summaries = []
    for query_id, group in grouped.items():
        best = min(group, key=lambda item: (item["evalue"], -item["score"]))
        summaries.append({
            "query_id": query_id,
            "n_hits": len(group),
            "n_nonself": sum(1 for item in group if not item["self"]),
            "best_evalue": best["evalue"],
            "best_evalue_label": best["evalue_label"],
            "best_target": best["target_id"],
            "length": group[0]["query_length"],
        })
    summaries.sort(key=lambda item: (item["best_evalue"], item["query_id"]))
    return summaries


def gel_svg(hits: list[dict], height: int = 168) -> str:
    if not hits:
        return ""
    width = 46
    values = [min(item["neglog10e"], 80.0) for item in hits]
    vmax = max(values) if values else 1.0
    vmax = max(vmax, 2.0)
    ticks = []
    for item, value in zip(hits, values):
        y = 12 + (1.0 - value / vmax) * (height - 24)
        ticks.append(
            f'<line class="gel-tick {escape(item["eclass"])}" data-hit="{item["index"]}" '
            f'x1="18" x2="34" y1="{y:.1f}" y2="{y:.1f}"/>'
        )
    labels = [
        f'<text x="14" y="14" class="gel-lab">{-0:.0f}</text>',
        f'<text x="14" y="{height/2 + 4:.0f}" class="gel-lab">{vmax/2:.0f}</text>',
        f'<text x="14" y="{height - 4}" class="gel-lab">{vmax:.0f}</text>',
    ]
    return (
        f'<svg class="gel" viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
        f'aria-hidden="true"><rect class="gel-lane" x="20" y="8" width="12" height="{height-16}" rx="2"/>'
        + "".join(ticks) + "".join(labels) + "</svg>"
    )


def rail_html(start: int, end: int, length: int, colour: str) -> str:
    if length <= 0:
        return '<div class="rail empty"></div>'
    left = 100.0 * start / length
    width = 100.0 * max(end - start, 0) / length
    return (
        f'<div class="rail" title="{start}–{end} / {length}">'
        f'<div class="rail-track"></div>'
        f'<div class="rail-hit" style="left:{left:.3f}%;width:{width:.3f}%;background:{escape(colour)}"></div>'
        f'<div class="rail-meta"><span>0</span><span>{length}</span></div>'
        f"</div>"
    )


def marked_seq(seq: str, start: int, end: int) -> str:
    if not seq:
        return '<span class="empty-seq">no sequence in this table</span>'
    start = max(0, min(start, len(seq)))
    end = max(start, min(end, len(seq)))
    return (
        f'<span class="dim">{escape(seq[:start])}</span>'
        f'<span class="lit">{escape(seq[start:end])}</span>'
        f'<span class="dim">{escape(seq[end:])}</span>'
    )


def plot_cell(svg: str | None, caption: str) -> str:
    if svg:
        return (
            f'<figure class="plot"><figcaption>{escape(caption)}</figcaption>{svg}</figure>'
        )
    return f'<div class="plot plot-empty"><span class="muted">No {escape(caption)}</span></div>'


def plot_panel(hit: dict) -> str:
    rows = (
        ("RNArtistCore", "plot_rn_query", "plot_rn_target", "query", "target"),
        ("R4RNA", "plot_r4_query", "plot_r4_target", "query", "target"),
        ("Alignment", "plot_sw_similarity", "plot_sw_scores", "cosine", "SW scores"),
    )
    body = []
    for label, q_key, t_key, q_cap, t_cap in rows:
        query_svg = hit.get(q_key)
        target_svg = hit.get(t_key)
        if not query_svg and not target_svg:
            continue
        body.append(
            "<tr>"
            f'<th scope="row">{escape(label)}</th>'
            f"<td>{plot_cell(query_svg, q_cap)}</td>"
            f"<td>{plot_cell(target_svg, t_cap)}</td>"
            "</tr>"
        )
    if not body:
        return (
            '<p class="muted">No plots for this hit. Re-run with --plot_backend '
            "rnartistcore/r4rna/both and/or --plot_sw.</p>"
        )
    return (
        '<section class="plot-panel" aria-label="Structure and alignment plots">'
        '<table class="plot-grid">'
        "<thead><tr><th></th><th>Query</th><th>Target</th></tr></thead>"
        f"<tbody>{''.join(body)}</tbody>"
        "</table></section>"
    )


def hit_article(hit: dict, colour: str) -> str:
    plot_block = plot_panel(hit)
    base_id = hit["base_identity"]
    if base_id is None and hit["aligned_columns"]:
        base_id = 100.0 * hit["pair_identity"]
    aln = hit["alignment_text"]
    aln_html = f"<pre class='aln'>{escape(aln)}</pre>" if aln else ""
    self_pill = '<span class="pill self">self</span>' if hit["self"] else ""
    return f"""
<article class="hit" id="hit-{hit['index']}" data-query="{escape(hit['query_id'])}" data-self="{str(hit['self']).lower()}" data-e="{hit['evalue']}" data-eclass="{hit['eclass']}" hidden>
  <header class="hit-head">
    <div class="who">
      <p class="eyebrow">hit {hit['rank']:03d} · cluster {escape(str(hit['cluster_id']))} {self_pill}</p>
      <h2><span class="qid">{escape(hit['query_id'])}</span><span class="vs">→</span><span class="tid">{escape(hit['target_id'])}</span></h2>
    </div>
    <dl class="metrics">
      <div><dt>E</dt><dd class="e {hit['eclass']}">{escape(hit['evalue_label'])}</dd></div>
      <div><dt>bits</dt><dd>{hit['bit_score']:.1f}</dd></div>
      <div><dt>score</dt><dd>{hit['score']:.2f}</dd></div>
      <div><dt>base id</dt><dd>{escape(fmt_pct(base_id))}</dd></div>
      <div><dt>structure id</dt><dd>{escape(fmt_pct(hit["structure_identity"]))}</dd></div>
    </dl>
  </header>
  <section class="spans" aria-label="Aligned spans">
    <div class="span-row">
      <span class="span-lab">query {hit['query_start']}–{hit['query_end']}</span>
      {rail_html(hit['query_start'], hit['query_end'], hit['query_length'], colour)}
    </div>
    <div class="span-row">
      <span class="span-lab">target {hit['target_start']}–{hit['target_end']}</span>
      {rail_html(hit['target_start'], hit['target_end'], hit['target_length'], colour)}
    </div>
  </section>
  <section class="letters" aria-label="Sequences">
    <p class="seq-lab">query sequence</p>
    <p class="seq">{marked_seq(hit['query_sequence'], hit['query_start'], hit['query_end'])}</p>
    <p class="seq-lab">query structure</p>
    <p class="seq">{marked_seq(hit['query_structure'], hit['query_start'], hit['query_end'])}</p>
    <p class="seq-lab">target sequence</p>
    <p class="seq">{marked_seq(hit['target_sequence'], hit['target_start'], hit['target_end'])}</p>
    <p class="seq-lab">target structure</p>
    <p class="seq">{marked_seq(hit['target_structure'], hit['target_start'], hit['target_end'])}</p>
  </section>
  {aln_html}
  {plot_block}
</article>
"""


CSS = """
:root {
  --bench: #1e2a2e;
  --bench-2: #263338;
  --glass: #f3f7f6;
  --card: #ffffff;
  --ink: #1a2426;
  --mute: #5e7074;
  --line: #d5e0de;
  --pair: #0e8f78;
  --pair-soft: #d7f3ec;
  --warn: #c45c26;
  --ns: #8a9698;
  --good: #1f6f8b;
  --strong: #0e8f78;
  --weak: #c45c26;
  --self: #7a5c2e;
  --mono: ui-monospace, "Cascadia Mono", "Source Code Pro", "Noto Sans Mono", Menlo, Consolas, monospace;
  --sans: "Segoe UI", "Source Sans 3", "IBM Plex Sans", system-ui, sans-serif;
  --serif: Palatino, "Palatino Linotype", "Book Antiqua", "Iowan Old Style", Georgia, serif;
}
* { box-sizing: border-box; }
html, body { margin: 0; padding: 0; background: var(--glass); color: var(--ink); font-family: var(--sans); }
body { min-height: 100vh; }
.mast {
  background:
    radial-gradient(1200px 280px at 80% -40%, rgba(14,143,120,.28), transparent 60%),
    linear-gradient(180deg, var(--bench) 0%, var(--bench-2) 100%);
  color: #e7eeec;
  padding: 2.2rem 6vw 1.6rem;
  border-bottom: 4px solid var(--pair);
}
.mast-kicker {
  font-family: var(--mono);
  letter-spacing: .28em;
  text-transform: uppercase;
  font-size: .68rem;
  color: #9db5b0;
  margin: 0 0 .4rem;
}
.mast h1 {
  font-family: var(--serif);
  font-weight: 500;
  font-size: clamp(1.8rem, 3vw, 2.6rem);
  margin: 0 0 .35rem;
  letter-spacing: -.02em;
}
.mast p.lede { margin: 0; color: #c5d4d0; max-width: 46rem; line-height: 1.45; }
.stats {
  display: flex; flex-wrap: wrap; gap: 1.6rem 2.4rem;
  margin: 1.4rem 0 0; padding: 0; list-style: none;
}
.stats li { min-width: 6.5rem; }
.stats b {
  display: block; font-family: var(--mono); font-size: 1.15rem; font-weight: 600; color: #fff;
}
.stats span { display: block; font-size: .72rem; letter-spacing: .08em; text-transform: uppercase; color: #9db5b0; margin-top: .15rem; }
.shell { display: grid; grid-template-columns: 220px minmax(0, 1fr) 64px; gap: 0; }
.queries {
  background: #e7eeec;
  border-right: 1px solid var(--line);
  padding: 1.2rem 1rem 2rem;
  position: sticky; top: 0; align-self: start; max-height: 100vh; overflow: auto;
}
.queries h2, .controls h2, .gel-wrap h2 {
  font-size: .68rem; letter-spacing: .16em; text-transform: uppercase; color: var(--mute);
  margin: 0 0 .7rem; font-weight: 600;
}
.query-btn {
  display: block; width: 100%; text-align: left; border: 0; background: transparent;
  padding: .45rem .5rem; margin: 0 0 .15rem; border-radius: 4px; cursor: pointer;
  font: inherit; color: inherit;
}
.query-btn:hover, .query-btn:focus-visible { background: #fff; outline: 2px solid var(--pair); }
.query-btn.active { background: var(--card); box-shadow: inset 3px 0 0 var(--pair); }
.query-btn .name { display: block; font-family: var(--mono); font-size: .78rem; word-break: break-all; }
.query-btn .meta { display: block; color: var(--mute); font-size: .72rem; margin-top: .1rem; }
.main { padding: 1.2rem 4vw 4rem; min-width: 0; }
.controls {
  display: flex; flex-wrap: wrap; gap: .7rem 1rem; align-items: end;
  margin-bottom: 1rem; padding-bottom: .9rem; border-bottom: 1px solid var(--line);
}
label.ctl { display: flex; flex-direction: column; gap: .25rem; font-size: .75rem; color: var(--mute); }
input[type="search"], select {
  font: inherit; border: 1px solid var(--line); background: var(--card); color: var(--ink);
  padding: .35rem .5rem; min-width: 12rem; border-radius: 3px;
}
.check { flex-direction: row; align-items: center; gap: .4rem; min-height: 2rem; }
.table-wrap { overflow: auto; border: 1px solid var(--line); background: var(--card); }
table.hits { width: 100%; border-collapse: collapse; font-size: .86rem; }
table.hits th {
  text-align: left; font-size: .68rem; letter-spacing: .1em; text-transform: uppercase;
  color: var(--mute); font-weight: 600; padding: .55rem .65rem; border-bottom: 1px solid var(--line);
  background: #eef3f2; position: sticky; top: 0;
}
table.hits td { padding: .5rem .65rem; border-bottom: 1px solid #edf2f1; vertical-align: top; }
table.hits tr { cursor: pointer; }
table.hits tr:hover td { background: var(--pair-soft); }
table.hits tr.selected td { background: #e4f6f1; }
table.hits tr[hidden] { display: none; }
.idc { font-family: var(--mono); font-size: .78rem; word-break: break-all; }
.e { font-family: var(--mono); font-variant-numeric: tabular-nums; }
.e.strong { color: var(--strong); font-weight: 650; }
.e.good { color: var(--good); }
.e.weak { color: var(--warn); }
.e.ns { color: var(--ns); }
.pill {
  display: inline-block; font-size: .65rem; letter-spacing: .08em; text-transform: uppercase;
  padding: .1rem .35rem; border: 1px solid currentColor; border-radius: 99px; margin-left: .35rem;
}
.pill.self { color: var(--self); }
.gel-wrap {
  padding: 1.2rem .4rem 2rem 0;
  border-left: 1px solid var(--line);
  position: sticky; top: 0; align-self: start;
  background: var(--glass);
}
.gel-wrap h2 { writing-mode: vertical-rl; transform: rotate(180deg); margin: 0 auto .6rem; }
.gel { display: block; margin: 0 auto; }
.gel-lane { fill: #d5e0de; }
.gel-tick { stroke: var(--ns); stroke-width: 1.4; }
.gel-tick.strong { stroke: var(--strong); stroke-width: 2; }
.gel-tick.good { stroke: var(--good); }
.gel-tick.weak { stroke: var(--warn); }
.gel-tick.on { stroke: var(--ink); stroke-width: 2.4; }
.gel-lab { fill: var(--mute); font-size: 7px; font-family: var(--mono); text-anchor: end; }
.hit {
  margin-top: 1.4rem; background: var(--card); border: 1px solid var(--line);
  padding: 1.2rem 1.3rem 1.4rem;
}
.hit[hidden] { display: none; }
.hit-head { display: flex; flex-wrap: wrap; justify-content: space-between; gap: 1rem; }
.eyebrow { margin: 0 0 .25rem; font-size: .7rem; letter-spacing: .14em; text-transform: uppercase; color: var(--mute); }
.hit h2 { margin: 0; font-family: var(--serif); font-weight: 500; font-size: 1.25rem; line-height: 1.25; }
.vs { color: var(--mute); padding: 0 .35rem; font-family: var(--sans); }
.metrics { display: flex; flex-wrap: wrap; gap: .8rem 1.2rem; margin: 0; }
.metrics div { margin: 0; }
.metrics dt { font-size: .68rem; letter-spacing: .12em; text-transform: uppercase; color: var(--mute); }
.metrics dd { margin: .1rem 0 0; font-family: var(--mono); font-size: 1rem; }
.spans { margin: 1.1rem 0 .4rem; display: grid; gap: .55rem; }
.span-row { display: grid; grid-template-columns: 11rem minmax(0,1fr); gap: .7rem; align-items: center; }
.span-lab { font-family: var(--mono); font-size: .75rem; color: var(--mute); }
.rail { position: relative; height: 1.35rem; }
.rail-track { position: absolute; left: 0; right: 0; top: .45rem; height: .28rem; background: #d5e0de; border-radius: 99px; }
.rail-hit { position: absolute; top: .28rem; height: .62rem; border-radius: 2px; min-width: 2px; opacity: .95; }
.rail-meta { position: absolute; inset: auto 0 0 0; display: flex; justify-content: space-between; font-size: .65rem; color: var(--mute); font-family: var(--mono); }
.letters { margin-top: .9rem; }
.seq-lab { margin: .65rem 0 .15rem; font-size: .68rem; letter-spacing: .12em; text-transform: uppercase; color: var(--mute); }
.seq {
  margin: 0; font-family: var(--mono); font-size: .72rem; line-height: 1.55;
  word-break: break-all; background: #f7faf9; padding: .4rem .5rem; border-left: 3px solid var(--line);
}
.seq .dim { color: #9aa8aa; }
.seq .lit { color: var(--ink); background: var(--pair-soft); box-shadow: inset 0 -2px 0 var(--pair); }
.aln {
  margin: 1rem 0 0; overflow: auto; background: var(--bench); color: #d5e4e0;
  padding: .8rem 1rem; font-size: .72rem; line-height: 1.35;
}
.plot-panel { margin-top: 1.1rem; overflow-x: auto; }
.plot-grid { width: 100%; border-collapse: collapse; table-layout: fixed; min-width: 36rem; }
.plot-grid th, .plot-grid td { border: 1px solid var(--line); vertical-align: top; padding: .45rem; }
.plot-grid thead th {
  text-align: center; font-size: .68rem; letter-spacing: .1em; text-transform: uppercase;
  color: var(--mute); font-weight: 600; background: #eef3f2;
}
.plot-grid tbody th {
  width: 7.4rem; font-size: .68rem; letter-spacing: .08em; text-transform: uppercase;
  color: var(--mute); font-weight: 600; background: #f7faf9; vertical-align: middle;
}
.plot-grid td { width: calc((100% - 7.4rem) / 2); background: #fff; }
.plot { margin: 0; background: #fff; padding: .2rem; }
.plot figcaption { font-size: .7rem; letter-spacing: .08em; text-transform: uppercase; color: var(--mute); margin: .1rem .2rem .35rem; }
.plot svg { width: 100%; height: auto; display: block; background: #fff; }
.plot-empty { min-height: 4.5rem; display: flex; align-items: center; justify-content: center; background: #f7faf9; }
.pager { display: flex; align-items: center; gap: .45rem; margin-left: auto; }
.pager button {
  font: inherit; border: 1px solid var(--line); background: var(--card); color: var(--ink);
  padding: .3rem .55rem; border-radius: 3px; cursor: pointer;
}
.pager button:disabled { opacity: .45; cursor: default; }
.muted, .empty { color: var(--mute); }
.methods {
  margin: 2rem 6vw 3rem; padding-top: 1rem; border-top: 1px solid var(--line);
  color: var(--mute); font-size: .82rem; max-width: 50rem; line-height: 1.5;
}
.methods code { font-family: var(--mono); font-size: .78rem; color: var(--ink); }
@media (max-width: 900px) {
  .shell { grid-template-columns: 1fr; }
  .queries, .gel-wrap { position: static; max-height: none; border: 0; }
  .gel-wrap { display: none; }
  .span-row { grid-template-columns: 1fr; }
}
@media print {
  .queries, .controls, .gel-wrap, .mast { break-inside: avoid; }
  .hit { break-inside: avoid; display: block !important; }
  .aln { background: #fff; color: #000; border: 1px solid #ccc; }
  body { background: #fff; }
}
@media (prefers-reduced-motion: reduce) {
  * { transition: none !important; }
}
"""

JS = """
(function () {
  const state = { query: "all", emax: Infinity, hideSelf: false, q: "", selected: 0, page: 0, pageSize: 10 };
  const rows = Array.from(document.querySelectorAll("tr[data-hit]"));
  const cards = Array.from(document.querySelectorAll("article.hit"));
  const ticks = Array.from(document.querySelectorAll(".gel-tick"));
  const qbtns = Array.from(document.querySelectorAll(".query-btn"));
  function num(tr) { return Number(tr.getAttribute("data-hit")); }
  function visible(tr) {
    const e = Number(tr.getAttribute("data-e"));
    const self = tr.getAttribute("data-self") === "true";
    const query = tr.getAttribute("data-query");
    const hay = (tr.getAttribute("data-hay") || "").toLowerCase();
    if (state.query !== "all" && query !== state.query) return false;
    if (e > state.emax) return false;
    if (state.hideSelf && self) return false;
    if (state.q && !hay.includes(state.q)) return false;
    return true;
  }
  function filtered() {
    return rows.filter(visible).map(num);
  }
  function pageCount(shown) {
    return Math.max(1, Math.ceil(shown.length / state.pageSize) || 1);
  }
  function render() {
    const shown = filtered();
    const pages = pageCount(shown);
    if (state.page >= pages) state.page = pages - 1;
    if (state.page < 0) state.page = 0;
    const start = state.page * state.pageSize;
    const pageHits = shown.slice(start, start + state.pageSize);
    if (!pageHits.includes(state.selected)) state.selected = pageHits[0] ?? -1;
    rows.forEach((tr) => {
      const i = num(tr);
      tr.hidden = !pageHits.includes(i);
      tr.classList.toggle("selected", i === state.selected);
    });
    cards.forEach((card) => {
      const i = Number(card.id.replace("hit-", ""));
      card.hidden = !pageHits.includes(i);
    });
    ticks.forEach((tick) => tick.classList.toggle("on", Number(tick.getAttribute("data-hit")) === state.selected));
    qbtns.forEach((btn) => btn.classList.toggle("active", btn.getAttribute("data-query") === state.query));
    const n = document.getElementById("shown-count");
    if (n) n.textContent = String(shown.length);
    const status = document.getElementById("page-status");
    if (status) status.textContent = (shown.length ? state.page + 1 : 0) + " / " + (shown.length ? pages : 0);
    const prev = document.getElementById("page-prev");
    const next = document.getElementById("page-next");
    if (prev) prev.disabled = state.page <= 0;
    if (next) next.disabled = state.page >= pages - 1 || shown.length === 0;
  }
  function resetPage() { state.page = 0; }
  document.getElementById("flt-query")?.addEventListener("click", (ev) => {
    const btn = ev.target.closest(".query-btn");
    if (!btn) return;
    state.query = btn.getAttribute("data-query");
    resetPage();
    render();
  });
  document.getElementById("flt-e")?.addEventListener("change", (ev) => {
    const v = ev.target.value;
    state.emax = v === "all" ? Infinity : Number(v);
    resetPage();
    render();
  });
  document.getElementById("flt-self")?.addEventListener("change", (ev) => {
    state.hideSelf = ev.target.checked;
    resetPage();
    render();
  });
  document.getElementById("flt-q")?.addEventListener("input", (ev) => {
    state.q = ev.target.value.trim().toLowerCase();
    resetPage();
    render();
  });
  document.getElementById("flt-page-size")?.addEventListener("change", (ev) => {
    state.pageSize = Number(ev.target.value) || 10;
    resetPage();
    render();
  });
  document.getElementById("page-prev")?.addEventListener("click", () => {
    state.page -= 1;
    render();
  });
  document.getElementById("page-next")?.addEventListener("click", () => {
    state.page += 1;
    render();
  });
  document.getElementById("hit-table")?.addEventListener("click", (ev) => {
    const tr = ev.target.closest("tr[data-hit]");
    if (!tr || tr.hidden) return;
    state.selected = num(tr);
    render();
    document.getElementById("hit-" + state.selected)?.scrollIntoView({ block: "nearest" });
  });
  document.addEventListener("keydown", (ev) => {
    if (ev.target.matches("input, select, textarea")) return;
    const shown = filtered();
    const at = shown.indexOf(state.selected);
    if (ev.key === "j" && at < shown.length - 1) {
      state.selected = shown[at + 1];
      state.page = Math.floor((at + 1) / state.pageSize);
      render();
      document.getElementById("hit-" + state.selected)?.scrollIntoView({ block: "nearest" });
    }
    if (ev.key === "k" && at > 0) {
      state.selected = shown[at - 1];
      state.page = Math.floor((at - 1) / state.pageSize);
      render();
      document.getElementById("hit-" + state.selected)?.scrollIntoView({ block: "nearest" });
    }
  });
  render();
})();
"""


def render_html(hits: list[dict], queries: list[dict], evd: dict, meta: dict, colour: str) -> str:
    n_align = len(hits)
    n_query = len(queries)
    n_self = sum(1 for hit in hits if hit["self"])
    n_sig = sum(1 for hit in hits if hit["evalue"] < 0.01)
    best = hits[0]["evalue_label"] if hits else "—"
    generated = meta.get("generated", "")
    query_btns = [
        '<button type="button" class="query-btn active" data-query="all"><span class="name">all queries</span>'
        f'<span class="meta">{n_align} hits</span></button>'
    ]
    for item in queries:
        query_btns.append(
            f'<button type="button" class="query-btn" data-query="{escape(item["query_id"])}">'
            f'<span class="name">{escape(item["query_id"])}</span>'
            f'<span class="meta">{item["n_hits"]} hits · best {escape(item["best_evalue_label"])}</span></button>'
        )
    rows_html = []
    for hit in hits:
        base_id = hit["base_identity"]
        if base_id is None and hit["aligned_columns"]:
            base_id = 100.0 * hit["pair_identity"]
        self_mark = " · self" if hit["self"] else ""
        rows_html.append(
            f'<tr data-hit="{hit["index"]}" data-query="{escape(hit["query_id"])}" '
            f'data-self="{str(hit["self"]).lower()}" data-e="{hit["evalue"]}" '
            f'data-hay="{escape(hit["query_id"] + " " + hit["target_id"])}">'
            f'<td>{hit["rank"]}</td>'
            f'<td class="idc">{escape(hit["query_id"])}</td>'
            f'<td class="idc">{escape(hit["target_id"])}{self_mark}</td>'
            f'<td class="e {hit["eclass"]}">{escape(hit["evalue_label"])}</td>'
            f'<td>{hit["bit_score"]:.1f}</td>'
            f'<td>{escape(fmt_pct(base_id))}</td>'
            f'<td>{escape(fmt_pct(hit["structure_identity"]))}</td>'
            f'<td>{hit["query_start"]}–{hit["query_end"]}</td>'
            f'<td>{hit["target_start"]}–{hit["target_end"]}</td>'
            f"</tr>"
        )
    articles = "\n".join(hit_article(hit, colour) for hit in hits)
    empty = "" if hits else '<p class="empty">No alignments were produced. Check seed thresholds or cluster filters.</p>'
    lam = evd.get("lambda")
    k_value = evd.get("K")
    residues = evd.get("database_residues")
    n_db = evd.get("database_sequences")
    methods = (
        f"E = K m N exp(−λS), with m = query length and N = {residues} database residues"
        if residues else
        "E = K m N exp(−λS), with m = query length and N = searchable residues"
    )
    lam_txt = f"{lam:.4g}" if isinstance(lam, (int, float)) else "—"
    k_txt = f"{k_value:.4g}" if isinstance(k_value, (int, float)) else "—"
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>ginflow search report</title>
<style>{CSS}</style>
</head>
<body>
<header class="mast">
  <p class="mast-kicker">ginflow</p>
  <h1>Search report</h1>
  <p class="lede">Local GINFINITY-SW alignments ranked like BLAST, lowest database E-value first. The teal mark is the aligned span; gray is the rest of the molecule.</p>
  <ul class="stats">
    <li><b>{n_query}</b><span>queries</span></li>
    <li><b>{n_align}</b><span>alignments</span></li>
    <li><b>{n_sig}</b><span>E &lt; 0.01</span></li>
    <li><b>{n_self}</b><span>self matches</span></li>
    <li><b>{escape(str(best))}</b><span>best E</span></li>
    <li><b>{escape(str(n_db or "—"))}</b><span>database RNAs</span></li>
    <li><b>{escape(str(evd.get("n_clusters", "—")))}</b><span>clusters</span></li>
    <li><b>{escape(str(evd.get("n_seeds", "—")))}</b><span>seeds</span></li>
  </ul>
</header>
<div class="shell">
  <nav class="queries" id="flt-query" aria-label="Queries">{''.join(query_btns)}</nav>
  <main class="main">
    <div class="controls">
      <label class="ctl">E-value at most
        <select id="flt-e">
          <option value="all">no cutoff</option>
          <option value="1">1</option>
          <option value="0.01">0.01</option>
          <option value="0.001">1e-3</option>
          <option value="1e-5">1e-5</option>
          <option value="1e-10">1e-10</option>
        </select>
      </label>
      <label class="ctl">Find
        <input id="flt-q" type="search" placeholder="query or target id" autocomplete="off"/>
      </label>
      <label class="ctl check"><input id="flt-self" type="checkbox"/> Hide self matches</label>
      <label class="ctl">Results per page
        <select id="flt-page-size">
          <option value="10" selected>10</option>
          <option value="25">25</option>
          <option value="50">50</option>
          <option value="100">100</option>
          <option value="150">150</option>
        </select>
      </label>
      <nav class="pager" aria-label="Results pages">
        <button type="button" id="page-prev">Prev</button>
        <span id="page-status" class="muted">1 / 1</span>
        <button type="button" id="page-next">Next</button>
      </nav>
      <p class="muted" style="margin:0 0 .2rem">showing <strong id="shown-count">{n_align}</strong> · j / k moves</p>
    </div>
    <div class="table-wrap">
      <table class="hits" id="hit-table">
        <thead><tr><th>#</th><th>query</th><th>target</th><th>E</th><th>bits</th><th>base id</th><th>structure id</th><th>q span</th><th>t span</th></tr></thead>
        <tbody>{''.join(rows_html)}</tbody>
      </table>
    </div>
    {empty}
    {articles}
  </main>
  <aside class="gel-wrap" aria-label="E-value lane">
    <h2>−log10 E</h2>
    {gel_svg(hits)}
  </aside>
</div>
<footer class="methods">
  <p>{escape(methods)}. Fitted λ = <code>{escape(lam_txt)}</code>, K = <code>{escape(k_txt)}</code>
  from reverse-sequence null alignments ({escape(str(evd.get("fit_method", "length-aware Gumbel MLE")))}).
  Generated {escape(generated)}.</p>
</footer>
<script>{JS}</script>
</body>
</html>
"""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--alignment-text", type=Path)
    parser.add_argument("--evd", type=Path)
    parser.add_argument("--clusters", type=Path)
    parser.add_argument("--seeds", type=Path)
    parser.add_argument("--plots-rnartist", type=Path, nargs="*")
    parser.add_argument("--plots-r4rna", type=Path, nargs="*")
    parser.add_argument("--plots-sw", type=Path, nargs="*")
    parser.add_argument("--highlight-colour", default="#0E8F78")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = read_tsv(args.alignments)
    texts = parse_alignment_text(args.alignment_text)
    evd: dict = {}
    if args.evd and args.evd.is_file():
        payload = json.loads(args.evd.read_text())
        if isinstance(payload, dict):
            evd = payload
    rn_svgs = load_svgs(args.plots_rnartist)
    r4_svgs = load_svgs(args.plots_r4rna)
    sw_svgs = load_svgs(args.plots_sw)
    hits = build_hits(rows, texts, rn_svgs, r4_svgs, sw_svgs)
    queries = query_summaries(hits)
    if args.clusters and args.clusters.is_file():
        evd.setdefault("n_clusters", len(read_tsv(args.clusters)))
    if args.seeds and args.seeds.is_file():
        evd.setdefault("n_seeds", len(read_tsv(args.seeds)))
    html = render_html(
        hits,
        queries,
        evd,
        {"generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")},
        args.highlight_colour,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(html)
    print(f"wrote {args.output} ({len(hits)} hits, {len(queries)} queries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
