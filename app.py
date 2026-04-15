"""
NylB Variant Enrichment Pipeline — Interactive Dashboard
Dark scientific theme · Dash + Plotly
"""

import logging
import re
import threading
import subprocess
from pathlib import Path
from datetime import datetime

import yaml
from Bio import SeqIO

import pandas as pd
import pysam
import dash
from dash import dcc, html, Input, Output, State, no_update
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

# ─── Paths ────────────────────────────────────────────────────────────────────

ROOT        = Path(__file__).parent
RESULTS     = ROOT / "results"
QC_DIR      = RESULTS / "qc"
ALIGNED_DIR = RESULTS / "aligned"
REPORT_DIR  = RESULTS / "report"

ENRICHMENT_TSV  = REPORT_DIR / "enrichment_report.tsv"
VALIDATION_TXT  = REPORT_DIR / "validation_report.txt"
WT_FASTA        = ROOT / "reference" / "nylB_wildtype.fasta"
VARIANTS_FASTA  = ROOT / "reference" / "nylB_variants.fasta"

CONFIG_YAML       = ROOT / "config.yaml"
SAMPLES           = ["pre_selection", "post_selection"]
WINNER_THRESHOLD  = 1.2

# NylB gene length (bp) and estimated mean Nanopore read length — used for
# coverage-to-read-count approximation in the config hints.
NYLB_GENE_LEN    = 1179
MEAN_READ_LEN_EST = 1500

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s — %(name)s — %(levelname)s — %(message)s",
)
logger = logging.getLogger(__name__)

# ─── Config I/O ───────────────────────────────────────────────────────────────

def load_config() -> dict:
    with open(CONFIG_YAML) as f:
        return yaml.safe_load(f)


def save_config(cfg: dict) -> None:
    with open(CONFIG_YAML, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, sort_keys=False)

# ─── NylB Active Site ─────────────────────────────────────────────────────────
# Catalytic residues from Negoro et al. crystal structure studies of NylB
# (6-aminohexanoate cyclic dimer hydrolase, Arthrobacter sp. KI72).

ACTIVE_SITE_CODONS = {
    110: "Asn110 — oxyanion hole",
    112: "Ser112 — catalytic nucleophile",
    115: "Lys115 — activates Ser112",
    187: "Arg187 — substrate binding",
    215: "Tyr215 — general base",
}

# Map every 0-indexed nucleotide position → residue label
ACTIVE_SITE_NT: dict[int, str] = {}
for _codon, _label in ACTIVE_SITE_CODONS.items():
    _start = (_codon - 1) * 3
    for _p in range(_start, _start + 3):
        ACTIVE_SITE_NT[_p] = _label

# ─── Theme ────────────────────────────────────────────────────────────────────

C = {
    "bg":      "#0d1117",
    "surface": "#161b22",
    "border":  "#21262d",
    "blue":    "#58a6ff",
    "green":   "#3fb950",
    "yellow":  "#e3b341",
    "red":     "#f85149",
    "purple":  "#bc8cff",
    "text":    "#e6edf3",
    "muted":   "#8b949e",
}

MONO = "'IBM Plex Mono', 'JetBrains Mono', 'Fira Code', monospace"

PLOT_BASE = dict(
    paper_bgcolor=C["surface"],
    plot_bgcolor=C["bg"],
    font=dict(color=C["text"], family=MONO, size=11),
)

def plot_layout(**kwargs) -> dict:
    """Merge PLOT_BASE with per-chart overrides (margin, height, etc.)."""
    return {**PLOT_BASE, **kwargs}

# ─── Data Loaders ─────────────────────────────────────────────────────────────

def load_enrichment() -> pd.DataFrame | None:
    if not ENRICHMENT_TSV.exists():
        return None
    df = pd.read_csv(ENRICHMENT_TSV, sep="\t")
    df["is_winner"] = df["enrichment_ratio"] >= WINNER_THRESHOLD
    return df.sort_values("enrichment_ratio", ascending=False).reset_index(drop=True)


def parse_nanostat(path: Path) -> dict:
    if not path.exists():
        return {}
    text = path.read_text()
    patterns = {
        "mean_length":    r"Mean read length:\s+([\d,\.]+)",
        "mean_quality":   r"Mean read quality:\s+([\d,\.]+)",
        "median_length":  r"Median read length:\s+([\d,\.]+)",
        "median_quality": r"Median read quality:\s+([\d,\.]+)",
        "n_reads":        r"Number of reads:\s+([\d,\.]+)",
        "n50":            r"Read length N50:\s+([\d,\.]+)",
        "total_bases":    r"Total bases:\s+([\d,\.]+)",
        "q10_pct":        r">Q10:\s+\d+ \(([\d\.]+)%\)",
        "q15_pct":        r">Q15:\s+\d+ \(([\d\.]+)%\)",
        "q20_pct":        r">Q20:\s+\d+ \(([\d\.]+)%\)",
    }
    out = {}
    for key, pat in patterns.items():
        m = re.search(pat, text)
        if m:
            out[key] = float(m.group(1).replace(",", ""))
    return out


def load_qc() -> dict[str, dict]:
    return {s: parse_nanostat(QC_DIR / f"{s}_stats.txt") for s in SAMPLES}


def get_flagstat(bam: Path) -> dict:
    result = {"total": 0, "mapped": 0, "pct": 0.0}
    if not bam.exists():
        return result
    try:
        r = subprocess.run(["samtools", "flagstat", str(bam)],
                           capture_output=True, text=True, timeout=30)
        m_total  = re.search(r"(\d+) \+ \d+ in total", r.stdout)
        m_mapped = re.search(r"(\d+) \+ \d+ mapped",   r.stdout)
        if m_total:  result["total"]  = int(m_total.group(1))
        if m_mapped: result["mapped"] = int(m_mapped.group(1))
        if result["total"]:
            result["pct"] = 100.0 * result["mapped"] / result["total"]
    except Exception as e:
        logger.error("flagstat failed for %s: %s", bam, e)
    return result


def get_idxstats(bam: Path) -> dict[str, int]:
    if not bam.exists():
        return {}
    try:
        raw = pysam.idxstats(str(bam))
        counts = {}
        for line in raw.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 3 and parts[0] != "*":
                counts[parts[0]] = int(parts[2])
        return counts
    except Exception as e:
        logger.error("idxstats failed for %s: %s", bam, e)
        return {}


def load_alignment() -> dict:
    return {
        s: {
            "flagstat":  get_flagstat(ALIGNED_DIR / f"{s}.bam"),
            "idxstats":  get_idxstats(ALIGNED_DIR / f"{s}.bam"),
        }
        for s in SAMPLES
    }


def read_fasta_sequences(path: Path) -> dict[str, str]:
    """Return {name: sequence} from a multi-FASTA file."""
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}


def get_mutations() -> dict[str, list[dict]]:
    if not WT_FASTA.exists() or not VARIANTS_FASTA.exists():
        return {}
    wt_seqs  = read_fasta_sequences(WT_FASTA)
    var_seqs = read_fasta_sequences(VARIANTS_FASTA)
    if not wt_seqs:
        return {}
    # Wildtype is the first (and only) entry in the WT file
    wt = next(iter(wt_seqs.values()))

    result = {}
    for vname, vseq in var_seqs.items():
        muts = []
        for i, (w, v) in enumerate(zip(wt, vseq)):
            if w != v:
                muts.append({
                    "nt_pos":           i + 1,          # 1-indexed
                    "codon":            i // 3 + 1,
                    "wt_base":          w,
                    "var_base":         v,
                    "is_active_site":   i in ACTIVE_SITE_NT,
                    "active_site_label": ACTIVE_SITE_NT.get(i, ""),
                })
        result[vname] = muts
    return result


# ─── Pipeline State ───────────────────────────────────────────────────────────

RULES = ["simulate_reads", "merge_fastqs", "qc", "trim", "align", "count_variants", "validate"]

_lock    = threading.Lock()
_log:    list[str] = []
_status: dict = {"running": False, "rules": {r: "pending" for r in RULES}}


def _stream(proc: subprocess.Popen, rule_context: str | None = None) -> int:
    """Stream a subprocess's stdout into _log, optionally tracking rule status."""
    for line in proc.stdout:
        with _lock:
            _log.append(line)
            if len(_log) > 10_000:
                del _log[:5_000]
            if rule_context is None:
                continue
            m_start = re.search(r"rule (\w+):", line)
            if m_start and m_start.group(1) in RULES:
                _status["rules"][m_start.group(1)] = "running"
            if re.search(r"Finished job \d+", line):
                for r in RULES:
                    if _status["rules"][r] == "running":
                        _status["rules"][r] = "complete"
            if re.search(r"\bError\b|\bException\b", line):
                for r in RULES:
                    if _status["rules"][r] == "running":
                        _status["rules"][r] = "failed"
    proc.wait()
    return proc.returncode


def _log_step(msg: str) -> None:
    with _lock:
        _log.append(f"[{datetime.now():%H:%M:%S}] {msg}\n")


def _run_pipeline():
    with _lock:
        _status["running"] = True
        _status["rules"]   = {r: "pending" for r in RULES}
        _log.clear()

    try:
        # ── Step 1: generate variant sequences ────────────────────────────
        _log_step("Step 1/3 — Generating variant sequences…")
        proc = subprocess.Popen(
            ["python3", "scripts/generate_variants.py"],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, cwd=str(ROOT),
        )
        rc = _stream(proc)
        if rc != 0:
            _log_step(f"generate_variants.py failed (exit {rc}). Aborting.")
            with _lock:
                _status["running"] = False
            return

        # ── Step 2: write per-variant FASTAs for both libraries ───────────
        _log_step("Step 2/3 — Preparing individual library FASTAs…")
        proc = subprocess.Popen(
            ["python3", "scripts/prepare_libraries.py"],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, cwd=str(ROOT),
        )
        rc = _stream(proc)
        if rc != 0:
            _log_step(f"prepare_libraries.py failed (exit {rc}). Aborting.")
            with _lock:
                _status["running"] = False
            return

        # ── Step 3: run Snakemake ─────────────────────────────────────────
        # --forceall: rerun every rule regardless of existing output files,
        # so that config changes always propagate fully.
        _log_step("Step 3/3 — Running Snakemake (--forceall)…")
        proc = subprocess.Popen(
            ["snakemake", "--cores", "all", "--forceall",
             "--rerun-incomplete", "--printshellcmds"],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, cwd=str(ROOT),
        )
        rc = _stream(proc, rule_context="snakemake")

        with _lock:
            if rc == 0:
                _log_step("Pipeline finished successfully.")
                for r in RULES:
                    if _status["rules"][r] != "failed":
                        _status["rules"][r] = "complete"
            else:
                _log_step(f"Snakemake exited with code {rc}.")
            _status["running"] = False

    except Exception as e:
        _log_step(f"ERROR: {e}")
        with _lock:
            _status["running"] = False


# ─── Re-usable Components ─────────────────────────────────────────────────────

def stat_block(value: str, label: str, color: str = None) -> html.Div:
    return html.Div([
        html.Div(value, className="stat-val", style={"color": color or C["blue"]}),
        html.Div(label, className="stat-lbl"),
    ], className="stat-block")


def rule_badge(rule: str, status: str) -> html.Span:
    return html.Span(rule, className=f"rbadge rb-{status}", id=f"rbadge-{rule}")


def sec_header(text: str) -> html.Div:
    return html.Div(text, className="sec-hdr")


# ─── Tab: Pipeline Control ────────────────────────────────────────────────────

def _qty_hint(val: str) -> str:
    """Return a human-readable note like '≈ 153 reads' for a coverage string."""
    try:
        cov = int(str(val).rstrip("xX"))
        reads = round(cov * NYLB_GENE_LEN / MEAN_READ_LEN_EST)
        return f"≈ {reads:,} reads per variant"
    except (TypeError, ValueError):
        return ""


def tab_pipeline() -> html.Div:
    cfg  = load_config()
    sim  = cfg.get("simulation", {})
    n_v  = sim.get("n_variants",    10)
    n_m  = sim.get("n_mutations",   15)
    seed = sim.get("seed",          42)
    base = sim.get("base_quantity", "200x")
    win  = sim.get("winner_quantity","500x")
    winners = sim.get("winners", [])

    input_style = {
        "background": "#0d1117", "color": C["text"],
        "border": f"1px solid {C['border']}", "borderRadius": "4px",
        "fontFamily": MONO, "fontSize": "0.82rem",
        "padding": "0.3rem 0.5rem", "width": "100%",
    }
    label_style = {
        "color": C["muted"], "fontFamily": MONO,
        "fontSize": "0.68rem", "textTransform": "uppercase",
        "letterSpacing": "0.08em", "marginBottom": "3px",
        "display": "block",
    }

    def field(label, input_el):
        return html.Div([
            html.Span(label, style=label_style),
            input_el,
        ], style={"marginBottom": "0.9rem"})

    winner_options = [{"label": str(i), "value": i} for i in range(1, n_v + 1)]

    config_form = html.Div([
        sec_header("Configuration"),
        dbc.Row([
            # ── Simulation column ──────────────────────────────────────────
            dbc.Col([
                html.Div("Simulation", style={**label_style, "color": C["blue"],
                                              "marginBottom": "0.75rem", "fontSize": "0.72rem"}),
                field("Number of variants",
                    dbc.Input(id="cfg-n-variants", type="number", value=n_v,
                              min=1, max=50, step=1, style=input_style)),
                field("Mutations per variant",
                    dbc.Input(id="cfg-n-mutations", type="number", value=n_m,
                              min=1, max=100, step=1, style=input_style)),
                field("Random seed",
                    dbc.Input(id="cfg-seed", type="number", value=seed,
                              min=0, step=1, style=input_style)),
            ], md=3),
            # ── Coverage column ────────────────────────────────────────────
            dbc.Col([
                html.Div("Coverage", style={**label_style, "color": C["blue"],
                                            "marginBottom": "0.75rem", "fontSize": "0.72rem"}),
                field("Base quantity (non-winners)",
                    dbc.Input(id="cfg-base-qty", type="text", value=base,
                              placeholder="e.g. 200x", style=input_style)),
                html.Div(id="cfg-base-hint",
                         children=_qty_hint(base),
                         style={"color": C["muted"], "fontFamily": MONO,
                                "fontSize": "0.68rem", "marginTop": "-0.6rem",
                                "marginBottom": "0.9rem"}),
                field("Winner quantity",
                    dbc.Input(id="cfg-winner-qty", type="text", value=win,
                              placeholder="e.g. 500x", style=input_style)),
                html.Div(id="cfg-winner-hint",
                         children=_qty_hint(win),
                         style={"color": C["muted"], "fontFamily": MONO,
                                "fontSize": "0.68rem", "marginTop": "-0.6rem",
                                "marginBottom": "0.9rem"}),
            ], md=3),
            # ── Winners column ─────────────────────────────────────────────
            dbc.Col([
                html.Div("Winner variants", style={**label_style, "color": C["blue"],
                                                   "marginBottom": "0.75rem", "fontSize": "0.72rem"}),
                html.Div(
                    "Variants simulated at winner quantity — i.e. ground-truth positives.",
                    style={"color": C["muted"], "fontSize": "0.75rem",
                           "marginBottom": "0.65rem", "lineHeight": "1.5"},
                ),
                dbc.Checklist(
                    id="cfg-winners",
                    options=winner_options,
                    value=winners,
                    inline=True,
                    style={"fontFamily": MONO, "fontSize": "0.8rem", "color": C["text"]},
                    inputStyle={"marginRight": "3px"},
                    labelStyle={"marginRight": "10px", "color": C["text"]},
                ),
            ], md=4),
            # ── Save button ────────────────────────────────────────────────
            dbc.Col([
                html.Div(style={"height": "1.8rem"}),   # align with inputs
                dbc.Button(
                    "💾  Save Config", id="save-cfg-btn",
                    color="primary", size="sm", outline=True,
                    style={"fontFamily": MONO, "fontSize": "0.78rem",
                           "letterSpacing": "0.04em", "whiteSpace": "nowrap"},
                ),
                html.Div(id="cfg-msg", style={"marginTop": "0.6rem"}),
            ], md=2, style={"display": "flex", "flexDirection": "column",
                            "justifyContent": "flex-start"}),
        ]),
    ], className="sci-card", style={"marginBottom": "1.2rem"})

    run_panel = html.Div([
        dbc.Row([
            dbc.Col([
                sec_header("Run Pipeline"),
                dbc.Button(
                    "▶  Run Snakemake", id="run-btn", color="success", size="sm",
                    className="mb-2",
                    style={"fontFamily": MONO, "letterSpacing": "0.05em", "fontSize": "0.82rem"},
                ),
                html.Div(id="run-msg"),
            ], md=3),
            dbc.Col([
                sec_header("Rule Status"),
                html.Div(id="rule-badges",
                         children=[rule_badge(r, "pending") for r in RULES]),
            ], md=9),
        ], className="mb-3"),
        dbc.Row([
            dbc.Col([
                sec_header("Live Output"),
                html.Div(id="log-out", className="log-win",
                         children="No pipeline run yet."),
            ]),
        ]),
    ])

    return html.Div([
        config_form,
        run_panel,
        dcc.Interval(id="poll", interval=800, disabled=True),
    ], className="tab-pane-inner")


# ─── Tab: QC ─────────────────────────────────────────────────────────────────

def tab_qc() -> html.Div:
    qc = load_qc()

    sample_cards = []
    for sample in SAMPLES:
        s     = qc.get(sample, {})
        label = sample.replace("_", " ").title()

        warnings = []
        if s.get("mean_quality", 99) < 10:
            warnings.append("Mean quality below Q10 — low-confidence basecalls.")
        if s.get("q15_pct", 100) < 50:
            warnings.append("Fewer than 50% of reads pass Q15.")
        if s.get("n_reads", 9999) < 500:
            warnings.append("Low read count — consider re-sequencing.")

        n_reads   = f"{s['n_reads']:,.0f}"     if "n_reads"   in s else "—"
        mq        = f"{s['mean_quality']:.1f}" if "mean_quality" in s else "—"
        ml        = f"{s['mean_length']:,.0f}" if "mean_length" in s else "—"
        medq      = f"{s['median_quality']:.1f}" if "median_quality" in s else "—"
        n50       = f"{s['n50']:,.0f}"         if "n50"       in s else "—"
        tb        = f"{s['total_bases']/1e6:.2f} Mb" if "total_bases" in s else "—"

        content = [
            html.H6(label, style={"fontFamily": MONO, "color": C["text"], "marginBottom": "0.9rem"}),
            dbc.Row([
                dbc.Col(stat_block(n_reads, "Reads",         C["blue"]),   md=4, className="mb-2"),
                dbc.Col(stat_block(mq,      "Mean Quality",  C["green"]),  md=4, className="mb-2"),
                dbc.Col(stat_block(ml,      "Mean Len (bp)", C["purple"]), md=4, className="mb-2"),
            ]),
            dbc.Row([
                dbc.Col(stat_block(medq, "Median Quality"), md=4, className="mb-2"),
                dbc.Col(stat_block(n50,  "N50 (bp)"),       md=4, className="mb-2"),
                dbc.Col(stat_block(tb,   "Total Bases"),    md=4, className="mb-2"),
            ]),
        ]
        if warnings:
            content.append(html.Div([html.Div(f"⚠  {w}") for w in warnings], className="warn-box"))

        sample_cards.append(
            dbc.Col(html.Div(content, className="sci-card"), md=6, className="mb-3")
        )

    # Quality cutoff bar chart
    thresholds = ["Q10", "Q15", "Q20"]
    keys       = ["q10_pct", "q15_pct", "q20_pct"]
    fig = go.Figure()
    for i, sample in enumerate(SAMPLES):
        s = qc.get(sample, {})
        vals = [s.get(k, 0) for k in keys]
        fig.add_trace(go.Bar(
            name=sample.replace("_", " ").title(),
            x=thresholds, y=vals,
            text=[f"{v:.1f}%" for v in vals],
            textposition="auto",
            marker_color=[C["blue"], C["green"]][i],
        ))
    fig.update_layout(
        **PLOT_BASE,
        barmode="group",
        yaxis=dict(title="% Reads", range=[0, 108],
                   gridcolor=C["border"], zerolinecolor=C["border"]),
        xaxis=dict(gridcolor=C["border"]),
        legend=dict(bgcolor="rgba(0,0,0,0)"),
        height=270,
    )

    return html.Div([
        sec_header("Sequencing Quality Control"),
        dbc.Row(sample_cards),
        sec_header("Reads Passing Quality Cutoffs"),
        dcc.Graph(figure=fig, config={"displayModeBar": False}),
    ], className="tab-pane-inner")


# ─── Tab: Alignment ───────────────────────────────────────────────────────────

def tab_alignment() -> html.Div:
    data = load_alignment()

    rate_cards = []
    for sample in SAMPLES:
        fs    = data[sample]["flagstat"]
        pct   = fs["pct"]
        color = C["green"] if pct >= 90 else C["yellow"] if pct >= 70 else C["red"]
        label = sample.replace("_", " ").title()
        rate_cards.append(dbc.Col([
            sec_header(label),
            html.Div([
                dbc.Row([
                    dbc.Col(stat_block(f"{fs['total']:,}",  "Total Reads",   C["blue"]),  md=4),
                    dbc.Col(stat_block(f"{fs['mapped']:,}", "Mapped Reads",  C["green"]), md=4),
                    dbc.Col(stat_block(f"{pct:.1f}%",       "Mapping Rate",  color),      md=4),
                ]),
            ], className="sci-card"),
        ], md=6, className="mb-4"))

    # Per-variant read count chart
    all_variants = sorted({
        v for s in SAMPLES for v in data[s]["idxstats"]
    })
    fig = go.Figure()
    for i, sample in enumerate(SAMPLES):
        counts = [data[sample]["idxstats"].get(v, 0) for v in all_variants]
        fig.add_trace(go.Bar(
            name=sample.replace("_", " ").title(),
            x=[v.replace("nylB_", "") for v in all_variants],
            y=counts,
            marker_color=[C["blue"], C["green"]][i],
        ))
    fig.update_layout(
        **PLOT_BASE,
        barmode="group",
        yaxis=dict(title="Read Count",
                   gridcolor=C["border"], zerolinecolor=C["border"]),
        xaxis=dict(gridcolor=C["border"], tickangle=-40),
        legend=dict(bgcolor="rgba(0,0,0,0)"),
        height=320, margin=dict(t=24, b=70, l=56, r=16),
    )

    return html.Div([
        sec_header("Alignment Summary"),
        dbc.Row(rate_cards),
        sec_header("Reads per Variant"),
        dcc.Graph(figure=fig, config={"displayModeBar": False}),
    ], className="tab-pane-inner")


# ─── Tab: Enrichment Results ──────────────────────────────────────────────────

def tab_enrichment() -> html.Div:
    df = load_enrichment()
    if df is None or df.empty:
        return html.Div("Enrichment report not found. Run the pipeline first.",
                        style={"color": C["muted"], "padding": "2rem", "fontFamily": MONO})

    n_winners   = int(df["is_winner"].sum())
    top_ratio   = f'{df.iloc[0]["enrichment_ratio"]:.3f}×'
    top_variant = df.iloc[0]["variant"].replace("nylB_", "")

    # ── Bar chart ──────────────────────────────────────────────────────────
    bar_colors = [C["green"] if w else C["blue"] for w in df["is_winner"]]
    bar = go.Figure(go.Bar(
        x=[v.replace("nylB_", "") for v in df["variant"]],
        y=df["enrichment_ratio"],
        marker_color=bar_colors,
        customdata=df[["variant", "pre_count", "post_count", "enrichment_ratio"]].values,
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "Enrichment: %{y:.3f}×<br>"
            "Pre-count: %{customdata[1]:,d}<br>"
            "Post-count: %{customdata[2]:,d}<br>"
            "<extra></extra>"
        ),
    ))
    bar.add_hline(
        y=WINNER_THRESHOLD,
        line=dict(color=C["yellow"], dash="dash", width=1.5),
        annotation_text=f"≥{WINNER_THRESHOLD}× threshold",
        annotation_font_color=C["yellow"],
        annotation_position="top right",
    )
    bar.update_layout(
        **PLOT_BASE,
        yaxis=dict(title="Enrichment Ratio",
                   gridcolor=C["border"], zerolinecolor=C["border"]),
        xaxis=dict(gridcolor=C["border"], tickangle=-40),
        showlegend=False,
        height=340, margin=dict(t=30, b=70, l=60, r=16),
    )

    # ── Scatter ────────────────────────────────────────────────────────────
    scatter = go.Figure()
    for is_w, name, marker in [
        (False, "Not enriched", dict(color=C["blue"],  size=9,  opacity=0.75)),
        (True,  "Enriched ★",  dict(color=C["green"], size=13, symbol="star", opacity=0.9)),
    ]:
        mask = df["is_winner"] == is_w
        sub  = df[mask]
        scatter.add_trace(go.Scatter(
            x=sub["pre_freq"], y=sub["post_freq"],
            mode="markers+text" if is_w else "markers",
            name=name,
            marker=marker,
            text=[v.replace("nylB_", "") for v in sub["variant"]] if is_w else None,
            textposition="top center",
            textfont=dict(color=C["green"], size=10),
            customdata=sub["variant"].values,
            hovertemplate="<b>%{customdata}</b><br>Pre: %{x:.4f}<br>Post: %{y:.4f}<extra></extra>",
        ))
    max_freq = max(df["pre_freq"].max(), df["post_freq"].max()) * 1.12
    scatter.add_trace(go.Scatter(
        x=[0, max_freq], y=[0, max_freq],
        mode="lines",
        line=dict(color=C["border"], dash="dot", width=1),
        showlegend=False, hoverinfo="skip",
    ))
    scatter.update_layout(
        **PLOT_BASE,
        xaxis=dict(title="Pre-selection Frequency",
                   gridcolor=C["border"], zerolinecolor=C["border"]),
        yaxis=dict(title="Post-selection Frequency",
                   gridcolor=C["border"], zerolinecolor=C["border"]),
        legend=dict(bgcolor="rgba(0,0,0,0)", bordercolor=C["border"],
                    font=dict(size=10)),
        height=340,
    )

    return html.Div([
        dbc.Row([
            dbc.Col(stat_block(str(len(df)),    "Variants Tested",             C["blue"]),   md=3, className="mb-3"),
            dbc.Col(stat_block(str(n_winners),  f"Enriched  ≥{WINNER_THRESHOLD}×", C["green"]), md=3, className="mb-3"),
            dbc.Col(stat_block(top_ratio,       "Top Enrichment",              C["yellow"]), md=3, className="mb-3"),
            dbc.Col(stat_block(top_variant,     "Top Variant",                 C["purple"]), md=3, className="mb-3"),
        ]),
        dbc.Row([
            dbc.Col([
                sec_header("Enrichment Ratios"),
                dcc.Graph(figure=bar, config={"displayModeBar": False}),
            ], md=7),
            dbc.Col([
                sec_header("Pre vs Post Selection Frequency"),
                dcc.Graph(figure=scatter, config={"displayModeBar": False}),
            ], md=5),
        ], className="mb-3"),
        dbc.Row([
            dbc.Col([
                sec_header("Export"),
                dbc.Button("↓  Download TSV", id="dl-btn", color="secondary", size="sm",
                           style={"fontFamily": MONO, "fontSize": "0.78rem"}),
                dcc.Download(id="dl-tsv"),
            ]),
        ]),
    ], className="tab-pane-inner")


# ─── Tab: Mutation Viewer ─────────────────────────────────────────────────────

def tab_mutations() -> html.Div:
    mutations = get_mutations()
    df        = load_enrichment()
    winners   = set(df[df["is_winner"]]["variant"].tolist()) if df is not None else set()

    if not mutations:
        return html.Div("Variant sequences not found.",
                        style={"color": C["muted"], "padding": "2rem", "fontFamily": MONO})

    wt_seqs = read_fasta_sequences(WT_FASTA)
    if not wt_seqs:
        return html.Div("Wildtype sequence not found.",
                        style={"color": C["muted"], "padding": "2rem", "fontFamily": MONO})
    wt_len  = len(next(iter(wt_seqs.values())))

    variants   = sorted(mutations.keys())
    y_labels   = [
        f"★ {v.replace('nylB_', '')}" if v in winners else v.replace("nylB_", "")
        for v in variants
    ]

    # Build z and hover text matrices (variants × nucleotide positions)
    z, hover = [], []
    for var in variants:
        mut_map = {m["nt_pos"] - 1: m for m in mutations[var]}   # 0-indexed
        row_z, row_h = [], []
        for pos in range(wt_len):
            if pos in mut_map:
                m = mut_map[pos]
                row_z.append(2 if m["is_active_site"] else 1)
                tip = f"nt{m['nt_pos']}: {m['wt_base']}→{m['var_base']}  (codon {m['codon']})"
                if m["is_active_site"]:
                    tip += f"<br>⚠ {m['active_site_label']}"
                row_h.append(tip)
            else:
                row_z.append(0)
                row_h.append("")
        z.append(row_z)
        hover.append(row_h)

    # Discrete colorscale: 0→bg, 1→blue, 2→red  (zmin=0, zmax=2 → norm=0,0.5,1)
    colorscale = [
        [0.0,  C["bg"]],
        [0.45, C["bg"]],
        [0.5,  C["blue"]],
        [0.95, C["blue"]],
        [1.0,  C["red"]],
    ]

    heatmap_fig = go.Figure(go.Heatmap(
        z=z, x=list(range(1, wt_len + 1)), y=y_labels,
        text=hover,
        hovertemplate="<b>%{y}</b><br>%{text}<extra></extra>",
        colorscale=colorscale,
        showscale=False, zmin=0, zmax=2,
    ))
    # Dashed vertical lines at active site nt positions
    for pos0 in ACTIVE_SITE_NT:
        heatmap_fig.add_vline(
            x=pos0 + 1,
            line=dict(color=C["yellow"], width=1, dash="dot"),
        )
    heatmap_fig.update_layout(
        **PLOT_BASE,
        xaxis=dict(title="Nucleotide Position", gridcolor=C["border"]),
        yaxis=dict(gridcolor="rgba(0,0,0,0)"),
        height=len(variants) * 34 + 80,
        margin=dict(t=24, b=60, l=110, r=16),
    )

    # Active site reference table
    tbl_rows = [
        html.Tr([
            html.Td(f"Codon {codon}",
                    style={"fontFamily": MONO, "color": C["yellow"], "paddingRight": "1.5rem"}),
            html.Td(f"nt {(codon-1)*3+1}–{(codon-1)*3+3}",
                    style={"fontFamily": MONO, "color": C["muted"], "paddingRight": "1.5rem"}),
            html.Td(label, style={"color": C["text"]}),
        ])
        for codon, label in ACTIVE_SITE_CODONS.items()
    ]

    legend = html.Div([
        html.Span("■", style={"color": C["blue"],  "marginRight": "4px"}),
        html.Span("Mutation",         style={"color": C["muted"], "fontSize": "0.8rem", "marginRight": "1rem"}),
        html.Span("■", style={"color": C["red"],   "marginRight": "4px"}),
        html.Span("Active-site mutation", style={"color": C["muted"], "fontSize": "0.8rem", "marginRight": "1rem"}),
        html.Span("★", style={"color": C["green"], "marginRight": "4px"}),
        html.Span("Enriched variant", style={"color": C["muted"], "fontSize": "0.8rem"}),
    ], style={"marginBottom": "0.6rem"})

    return html.Div([
        sec_header("Mutation Map — variants vs wildtype"),
        legend,
        dcc.Graph(figure=heatmap_fig, config={"displayModeBar": False}),
        html.Hr(style={"borderColor": C["border"], "margin": "1.5rem 0"}),
        sec_header("NylB Active Site Residues"),
        html.P(
            "Residues from Negoro et al. crystal structure of NylB "
            "(6-aminohexanoate cyclic dimer hydrolase, Arthrobacter sp. KI72). "
            "Mutations at these positions are most likely to affect catalytic activity.",
            style={"color": C["muted"], "fontSize": "0.82rem", "maxWidth": "640px",
                   "lineHeight": "1.65", "marginBottom": "0.75rem"},
        ),
        dbc.Table(
            [html.Thead(html.Tr([
                html.Th("Residue",    style={"color": C["muted"], "fontWeight": "400", "fontSize": "0.75rem"}),
                html.Th("Nucleotides",style={"color": C["muted"], "fontWeight": "400", "fontSize": "0.75rem"}),
                html.Th("Role",       style={"color": C["muted"], "fontWeight": "400", "fontSize": "0.75rem"}),
            ])),
             html.Tbody(tbl_rows)],
            size="sm", bordered=False,
            style={"background": "transparent", "fontSize": "0.83rem", "maxWidth": "580px"},
        ),
    ], className="tab-pane-inner")


# ─── Tab: About ───────────────────────────────────────────────────────────────

def tab_about() -> html.Div:
    def item(title, body):
        return dbc.ListGroupItem([
            html.Strong(title, style={"color": C["text"]}),
            html.Span(f" — {body}", style={"color": C["muted"]}),
        ], style={"background": C["surface"], "border": f"1px solid {C['border']}",
                  "fontSize": "0.84rem"})

    stack = [
        ("Snakemake",  "workflow manager, runs all rules in order"),
        ("Badread",    "simulates realistic Nanopore sequencing reads"),
        ("NanoStat",   "reports read quality and length statistics"),
        ("Porechop",   "trims sequencing adapters from reads"),
        ("Minimap2",   "aligns long reads to the variant reference"),
        ("samtools",   "sorts, indexes, and queries BAM files"),
        ("pysam",      "Python interface to alignment data"),
        ("Dash + Plotly", "this interactive dashboard"),
    ]

    return html.Div([
        dbc.Row([
            dbc.Col([
                html.H5("What is this pipeline?",
                        style={"color": C["blue"], "fontFamily": MONO, "marginBottom": "0.9rem"}),
                html.P(
                    "This dashboard visualises the NylB directed evolution pipeline — a tool for simulating "
                    "and analysing selection experiments on NylB, the enzyme from Arthrobacter sp. KI72 "
                    "that degrades nylon-6 oligomers (6-aminohexanoate cyclic dimer hydrolase).",
                    style={"color": C["text"], "lineHeight": "1.72", "marginBottom": "1.2rem"},
                ),
                html.H6("How directed evolution works",
                        style={"color": C["text"], "marginBottom": "0.6rem"}),
                html.P(
                    "You create a library of enzyme variants — each carrying a handful of random mutations — "
                    "then apply a selection pressure (e.g., substrate availability, temperature). "
                    "Variants that survive become more abundant. "
                    "The enrichment ratio (post-selection frequency ÷ pre-selection frequency) "
                    f"quantifies how much a variant was favoured. Variants with ratio ≥ {WINNER_THRESHOLD}× "
                    "are flagged as enriched.",
                    style={"color": C["muted"], "lineHeight": "1.72", "marginBottom": "1.2rem"},
                ),
                html.H6("Pipeline steps",
                        style={"color": C["text"], "marginBottom": "0.6rem"}),
                dbc.ListGroup([
                    item("Simulate reads",  "Badread generates synthetic Nanopore reads for each variant at defined coverage"),
                    item("Merge",           "per-variant FASTQs are concatenated into library files"),
                    item("QC",              "NanoStat reports read count, quality, and length distributions"),
                    item("Trim",            "Porechop removes sequencing adapters"),
                    item("Align",           "Minimap2 maps reads to the 10-variant reference; samtools sorts and indexes the BAM"),
                    item("Count & Enrich",  "read counts per variant are tallied and normalised to frequencies; enrichment ratio = post / pre"),
                    item("Validate",        "checks that the known winner variants (2, 5, 8) appear at the top of the ranked list"),
                ], flush=True, style={"marginBottom": "1.4rem"}),
                html.H6("Reading the Mutation Viewer",
                        style={"color": C["text"], "marginBottom": "0.6rem"}),
                html.P(
                    "Each row is a variant; each column is a nucleotide position in the gene. "
                    "Blue cells are mutations relative to the wildtype. "
                    "Red cells fall within a known NylB active-site residue (Ser112, Lys115, Tyr215, Arg187, Asn110) "
                    "— mutations there are most likely to change enzymatic activity. "
                    "Dashed yellow lines mark active-site codon boundaries.",
                    style={"color": C["muted"], "lineHeight": "1.72"},
                ),
            ], md=8),
            dbc.Col([
                html.Div([
                    sec_header("Pipeline Stack"),
                    dbc.ListGroup([item(t, b) for t, b in stack], flush=True),
                ], className="sci-card"),
            ], md=4),
        ]),
    ], className="tab-pane-inner")


# ─── App Layout ───────────────────────────────────────────────────────────────

app = dash.Dash(
    __name__,
    external_stylesheets=[
        dbc.themes.DARKLY,
        "https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600"
        "&family=IBM+Plex+Sans:wght@300;400;600&display=swap",
    ],
    suppress_callback_exceptions=True,
    title="NylB Enrichment Dashboard",
)

HEADER = html.Div([
    dbc.Container([
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.Span("NylB", className="brand-primary"),
                    html.Span(" Enrichment", className="brand-sub"),
                ]),
                html.Div(
                    "Directed evolution dashboard  ·  Arthrobacter sp. KI72  ·  6-Ahx cyclic dimer hydrolase",
                    className="brand-tagline",
                ),
            ], md=9),
            dbc.Col([
                html.Div(id="last-updated",
                         style={"textAlign": "right", "color": C["muted"],
                                "fontFamily": MONO, "fontSize": "0.68rem",
                                "marginTop": "0.4rem"}),
            ], md=3),
        ], align="center"),
    ], fluid=True),
], className="dash-header")

app.layout = html.Div([
    HEADER,
    dbc.Container([
        dbc.Tabs(
            id="tabs",
            active_tab="tab-pipeline",
            children=[
                dbc.Tab(label="Pipeline Control", tab_id="tab-pipeline"),
                dbc.Tab(label="QC",               tab_id="tab-qc"),
                dbc.Tab(label="Alignment",        tab_id="tab-alignment"),
                dbc.Tab(label="Enrichment",       tab_id="tab-enrichment"),
                dbc.Tab(label="Mutation Viewer",  tab_id="tab-mutations"),
                dbc.Tab(label="About",            tab_id="tab-about"),
            ],
        ),
        html.Div(id="tab-content"),
    ], fluid=True),
    dcc.Interval(id="results-watcher", interval=15_000, n_intervals=0),
], style={"background": C["bg"], "minHeight": "100vh"})


# ─── Callbacks ────────────────────────────────────────────────────────────────

@app.callback(
    Output("tab-content",   "children"),
    Output("last-updated",  "children"),
    Input("tabs",           "active_tab"),
    Input("results-watcher","n_intervals"),
)
def render_tab(active_tab, _):
    tab_map = {
        "tab-pipeline":   tab_pipeline,
        "tab-qc":         tab_qc,
        "tab-alignment":  tab_alignment,
        "tab-enrichment": tab_enrichment,
        "tab-mutations":  tab_mutations,
        "tab-about":      tab_about,
    }
    content = tab_map.get(active_tab, tab_about)()

    ts = "No results yet"
    if ENRICHMENT_TSV.exists():
        mtime = datetime.fromtimestamp(ENRICHMENT_TSV.stat().st_mtime)
        ts = f"Results  {mtime:%Y-%m-%d %H:%M}"

    return content, ts


@app.callback(
    Output("run-btn",  "disabled"),
    Output("poll",     "disabled"),
    Output("run-msg",  "children"),
    Input("run-btn",   "n_clicks"),
    prevent_initial_call=True,
)
def start_pipeline(n_clicks):
    with _lock:
        already_running = _status["running"]
    if n_clicks and not already_running:
        threading.Thread(target=_run_pipeline, daemon=True).start()
        msg = html.Div("Pipeline started…",
                       style={"color": C["blue"], "fontFamily": MONO,
                              "fontSize": "0.8rem", "marginTop": "0.4rem"})
        return True, False, msg
    return no_update, no_update, no_update


@app.callback(
    Output("log-out",     "children"),
    Output("rule-badges", "children"),
    Output("run-btn",     "disabled",  allow_duplicate=True),
    Output("poll",        "disabled",  allow_duplicate=True),
    Input("poll",         "n_intervals"),
    prevent_initial_call=True,
)
def poll_pipeline(_):
    with _lock:
        log_text = "".join(_log[-300:])
        rules    = dict(_status["rules"])
        running  = _status["running"]

    badges          = [rule_badge(r, rules.get(r, "pending")) for r in RULES]
    btn_disabled    = running
    poll_disabled   = not running

    return log_text or "Waiting for output…", badges, btn_disabled, poll_disabled


@app.callback(
    Output("dl-tsv", "data"),
    Input("dl-btn",  "n_clicks"),
    prevent_initial_call=True,
)
def download_tsv(n_clicks):
    if n_clicks and ENRICHMENT_TSV.exists():
        return dcc.send_file(str(ENRICHMENT_TSV))
    return no_update


# ── Config: live coverage hints ───────────────────────────────────────────────

@app.callback(
    Output("cfg-base-hint",   "children"),
    Input("cfg-base-qty",     "value"),
    prevent_initial_call=True,
)
def update_base_hint(val):
    return _qty_hint(val or "")


@app.callback(
    Output("cfg-winner-hint", "children"),
    Input("cfg-winner-qty",   "value"),
    prevent_initial_call=True,
)
def update_winner_hint(val):
    return _qty_hint(val or "")


# ── Config: save ──────────────────────────────────────────────────────────────

@app.callback(
    Output("cfg-msg",       "children"),
    Input("save-cfg-btn",   "n_clicks"),
    State("cfg-n-variants", "value"),
    State("cfg-n-mutations","value"),
    State("cfg-seed",       "value"),
    State("cfg-base-qty",   "value"),
    State("cfg-winner-qty", "value"),
    State("cfg-winners",    "value"),
    prevent_initial_call=True,
)
def save_config_callback(n_clicks, n_variants, n_mutations, seed,
                         base_qty, winner_qty, winners):
    errors = []

    # ── Validate ─────────────────────────────────────────────────────────────
    try:
        n_variants = int(n_variants)
        assert 1 <= n_variants <= 50
    except (TypeError, ValueError, AssertionError):
        errors.append("Variants must be an integer 1–50.")

    try:
        n_mutations = int(n_mutations)
        assert 1 <= n_mutations <= 100
    except (TypeError, ValueError, AssertionError):
        errors.append("Mutations must be an integer 1–100.")

    try:
        seed = int(seed)
        assert seed >= 0
    except (TypeError, ValueError, AssertionError):
        errors.append("Seed must be a non-negative integer.")

    for qty_name, qty_val in [("Base quantity", base_qty), ("Winner quantity", winner_qty)]:
        if not re.fullmatch(r"\d+[xX]", str(qty_val or "").strip()):
            errors.append(f"{qty_name} must be a number followed by 'x' (e.g. 200x).")

    winners = sorted(set(int(w) for w in (winners or [])))
    if not errors and any(w < 1 or w > n_variants for w in winners):
        errors.append(f"All winner indices must be between 1 and {n_variants}.")

    if errors:
        return html.Div(
            [html.Div(f"✗  {e}") for e in errors],
            style={"color": C["red"], "fontFamily": MONO,
                   "fontSize": "0.75rem", "lineHeight": "1.7"},
        )

    # ── Write ─────────────────────────────────────────────────────────────────
    try:
        cfg = load_config()
        cfg["simulation"]["n_variants"]     = n_variants
        cfg["simulation"]["n_mutations"]    = n_mutations
        cfg["simulation"]["seed"]           = seed
        cfg["simulation"]["base_quantity"]  = str(base_qty).strip()
        cfg["simulation"]["winner_quantity"]= str(winner_qty).strip()
        cfg["simulation"]["winners"]        = winners
        save_config(cfg)
    except Exception as e:
        return html.Div(
            f"✗  Failed to write config.yaml: {e}",
            style={"color": C["red"], "fontFamily": MONO, "fontSize": "0.75rem"},
        )

    ts = datetime.now().strftime("%H:%M:%S")
    return html.Div(
        f"✓  Saved at {ts}",
        style={"color": C["green"], "fontFamily": MONO, "fontSize": "0.75rem"},
    )


# ─── Entry Point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)
