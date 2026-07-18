#!/usr/bin/env python3
"""
YithCOMPASM.py — Compare two genome assemblies: metrics, whole-genome dot
plots colored by percent identity, sequence correspondence, unaligned
regions, and lightweight rearrangement flagging.

Subcommands
-----------
    compare_assemblies    Compare two FASTA assemblies (metrics, alignment,
                          dot plot, correspondence, unaligned regions,
                          rearrangement candidates)

Usage
-----
    YithCOMPASM.py compare_assemblies \\
        --assembly_a genomeA.fasta --assembly_b genomeB.fasta \\
        --output results_run/ --threads 8
"""

import argparse
import csv
import getpass
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

VERSION = "v0.1.0"

MAX_LABELED_SEQS = 30   # only the largest N sequences per axis get gridlines/labels

_LOG_FH = None


# ── Logging ──────────────────────────────────────────────────────────────────

def _log(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, file=sys.stderr)
    if _LOG_FH is not None:
        print(line, file=_LOG_FH, flush=True)


def _banner(title: str) -> None:
    bar = "─" * (len(title) + 4)
    _log(f"┌{bar}┐")
    _log(f"│  {title}  │")
    _log(f"└{bar}┘")


_QUOTE_LINES = [
    "\"In its vast libraries were volumes of texts and",
    " pictures holding the whole of earth's annals—",
    " histories and descriptions of every species that",
    " had ever been or that ever would be.\"",
    "          — H.P. Lovecraft, The Shadow Out of Time (1936)",
]


def _print_quote() -> None:
    width = max(len(l) for l in _QUOTE_LINES) + 4
    border = "─" * width
    _log(f"┌{border}┐")
    for line in _QUOTE_LINES:
        padding = width - len(line) - 1
        _log(f"│ {line}{' ' * padding}│")
    _log(f"└{border}┘")


# ── External tool helpers ─────────────────────────────────────────────────────

def _require_tool(name: str) -> str:
    tool = shutil.which(name)
    if tool is None:
        print(f"ERROR: '{name}' not found in PATH.\n"
              f"       Install with:  conda install -c bioconda {name}",
              file=sys.stderr)
        sys.exit(1)
    return tool


def _run(cmd: list, stdout_path: Path = None, cwd: Path = None) -> subprocess.CompletedProcess:
    _log(f"  $ {' '.join(str(c) for c in cmd)}" + (f" > {stdout_path}" if stdout_path else ""))
    stdout_fh = open(stdout_path, "w") if stdout_path else subprocess.PIPE
    try:
        result = subprocess.run(cmd, stdout=stdout_fh, stderr=subprocess.PIPE, text=True, cwd=cwd)
    finally:
        if stdout_path:
            stdout_fh.close()
    if result.returncode != 0:
        print(f"ERROR: command failed (exit {result.returncode}):\n"
              f"{result.stderr[-3000:]}", file=sys.stderr)
        sys.exit(1)
    return result


# ── Checkpoint / resume ───────────────────────────────────────────────────────

def _checkpoint(path: Path, label: str, force: bool) -> bool:
    if not force and path.exists() and path.stat().st_size > 0:
        _log(f"  [checkpoint] {label} — {path.name} already exists, skipping")
        return True
    return False


# ── Input validation ──────────────────────────────────────────────────────────

def _validate_inputs(pairs: list) -> None:
    ok = True
    for flag, path in pairs:
        if not path.exists():
            print(f"ERROR: {flag} not found: {path}", file=sys.stderr)
            ok = False
    if not ok:
        sys.exit(1)


# ── FASTA length scan (streaming, no full sequence held in memory) ──────────

def scan_fasta_lengths(fasta_path: Path) -> dict:
    """Return {seq_id: length}, insertion order preserved."""
    lengths = {}
    seq_id = None
    length = 0
    n_count = 0
    gc_count = 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_id is not None:
                    lengths[seq_id] = length
                seq_id = line[1:].split()[0]
                length = 0
            else:
                s = line.strip().upper()
                length += len(s)
        if seq_id is not None:
            lengths[seq_id] = length
    return lengths


def scan_fasta_composition(fasta_path: Path) -> dict:
    """Return dict of assembly-wide metrics: n_seqs, total_length, min/max/mean/median,
    n50, l50, gc_content, n_count, n_percent."""
    lengths = []
    seq_id = None
    length = 0
    gc_count = 0
    n_count = 0
    total_bases = 0

    def _flush():
        nonlocal length
        if seq_id is not None:
            lengths.append(length)

    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_id is not None:
                    lengths.append(length)
                seq_id = line[1:].split()[0]
                length = 0
            else:
                s = line.strip().upper()
                length += len(s)
                total_bases += len(s)
                gc_count += s.count("G") + s.count("C")
                n_count += s.count("N")
        if seq_id is not None:
            lengths.append(length)

    if not lengths:
        return {"n_seqs": 0, "total_length": 0, "min_length": 0, "max_length": 0,
                "mean_length": 0, "median_length": 0, "n50": 0, "l50": 0,
                "gc_content": 0.0, "n_percent": 0.0}

    lengths_desc = sorted(lengths, reverse=True)
    total = sum(lengths)
    cumsum = 0
    n50 = l50 = 0
    for i, length_i in enumerate(lengths_desc, 1):
        cumsum += length_i
        if cumsum >= total * 0.5:
            n50, l50 = length_i, i
            break
    mid = len(lengths) // 2
    sorted_lens = sorted(lengths)
    median = (sorted_lens[mid] if len(lengths) % 2 else
              (sorted_lens[mid - 1] + sorted_lens[mid]) / 2)

    return {
        "n_seqs": len(lengths),
        "total_length": total,
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": round(total / len(lengths), 2),
        "median_length": median,
        "n50": n50,
        "l50": l50,
        "gc_content": round(gc_count / total_bases * 100, 4) if total_bases else 0.0,
        "n_percent": round(n_count / total_bases * 100, 4) if total_bases else 0.0,
    }


# ── Module 1: assembly metrics comparison ─────────────────────────────────────

def run_metrics_comparison(assembly_a: Path, assembly_b: Path, out_path: Path, force: bool) -> dict:
    if _checkpoint(out_path, "metrics comparison", force):
        with open(out_path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return {r["metric"]: (r["assembly_a"], r["assembly_b"]) for r in rows}

    metrics_a = scan_fasta_composition(assembly_a)
    metrics_b = scan_fasta_composition(assembly_b)
    order = ["n_seqs", "total_length", "min_length", "max_length", "mean_length",
             "median_length", "n50", "l50", "gc_content", "n_percent"]

    with open(out_path, "w") as fh:
        fh.write("metric\tassembly_a\tassembly_b\tdelta\n")
        for m in order:
            va, vb = metrics_a[m], metrics_b[m]
            delta = round(va - vb, 4) if isinstance(va, (int, float)) else "NA"
            fh.write(f"{m}\t{va}\t{vb}\t{delta}\n")

    _log(f"  Written: {out_path.name}")
    return {m: (metrics_a[m], metrics_b[m]) for m in order}


# ── Alignment ──────────────────────────────────────────────────────────────────

def run_minimap2(assembly_a: Path, assembly_b: Path, paf_path: Path,
                 preset: str, threads: int, force: bool) -> None:
    """assembly_a is used as query, assembly_b as target (PAF query-first convention)."""
    if _checkpoint(paf_path, "minimap2 alignment", force):
        return
    minimap2 = _require_tool("minimap2")
    # -c forces base-level CIGAR generation so PAF columns 10/11 (nmatch/alen)
    # reflect exact alignment identity, not the default minimizer-chain heuristic.
    cmd = [minimap2, "-x", preset, "-c", "-t", str(threads), str(assembly_b), str(assembly_a)]
    _run(cmd, stdout_path=paf_path)
    _log(f"  Written: {paf_path.name}")


def parse_paf(paf_path: Path, min_align_len: int, min_identity: float,
             include_secondary: bool) -> list:
    """
    Parse a PAF file into alignment records, sorted by (qname, qstart).
    Returns list of dicts: qname, qlen, qstart, qend, strand,
                            tname, tlen, tstart, tend, identity, alen
    """
    records = []
    with open(paf_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname, qlen, qstart, qend, strand = fields[0], int(fields[1]), int(fields[2]), int(fields[3]), fields[4]
            tname, tlen, tstart, tend = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            nmatch, alen = int(fields[9]), int(fields[10])

            if not include_secondary:
                tp = next((f[5:] for f in fields[12:] if f.startswith("tp:A:")), None)
                if tp == "S":
                    continue

            if alen < min_align_len:
                continue
            identity = (nmatch / alen * 100) if alen > 0 else 0.0
            if identity < min_identity:
                continue

            records.append({
                "qname": qname, "qlen": qlen, "qstart": qstart, "qend": qend,
                "strand": strand, "tname": tname, "tlen": tlen,
                "tstart": tstart, "tend": tend, "identity": identity, "alen": alen,
            })
    records.sort(key=lambda r: (r["qname"], r["qstart"]))
    return records


# ── Module 2: alignment summary ───────────────────────────────────────────────

def compute_alignment_summary(records: list, lens_a: dict, lens_b: dict) -> dict:
    total_len_a = sum(lens_a.values())
    total_len_b = sum(lens_b.values())

    if not records:
        return {
            "n_alignments": 0, "total_len_a": total_len_a, "total_len_b": total_len_b,
            "n_seqs_a": len(lens_a), "n_seqs_b": len(lens_b),
            "aligned_bp_query": 0, "coverage_a_pct": 0.0, "coverage_b_pct": 0.0,
            "mean_identity": 0.0, "weighted_identity": 0.0,
            "min_identity": 0.0, "max_identity": 0.0,
        }

    aligned_bp_query = sum(r["qend"] - r["qstart"] for r in records)
    aligned_bp_target = sum(r["tend"] - r["tstart"] for r in records)
    identities = [r["identity"] for r in records]
    weights = [r["alen"] for r in records]
    weighted_identity = sum(i * w for i, w in zip(identities, weights)) / sum(weights)

    return {
        "n_alignments": len(records),
        "total_len_a": total_len_a, "total_len_b": total_len_b,
        "n_seqs_a": len(lens_a), "n_seqs_b": len(lens_b),
        "aligned_bp_query": aligned_bp_query,
        "coverage_a_pct": round(aligned_bp_query / total_len_a * 100, 4) if total_len_a else 0.0,
        "coverage_b_pct": round(aligned_bp_target / total_len_b * 100, 4) if total_len_b else 0.0,
        "mean_identity": round(sum(identities) / len(identities), 4),
        "weighted_identity": round(weighted_identity, 4),
        "min_identity": round(min(identities), 4),
        "max_identity": round(max(identities), 4),
    }


def write_alignment_summary(path: Path, stats: dict, assembly_a: Path, assembly_b: Path,
                            preset: str, command: str) -> None:
    with open(path, "w") as fh:
        fh.write(f"YithCOMPASM.py {VERSION} — compare_assemblies\n")
        fh.write(f"Date       : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        fh.write(f"Command    : {command}\n\n")
        fh.write(f"Assembly A : {assembly_a}\n")
        fh.write(f"  sequences: {stats['n_seqs_a']:,}\n")
        fh.write(f"  length   : {stats['total_len_a']:,} bp\n")
        fh.write(f"Assembly B : {assembly_b}\n")
        fh.write(f"  sequences: {stats['n_seqs_b']:,}\n")
        fh.write(f"  length   : {stats['total_len_b']:,} bp\n\n")
        fh.write(f"minimap2 preset      : {preset}\n")
        fh.write(f"Alignment blocks kept: {stats['n_alignments']:,}\n")
        fh.write(f"Aligned bp (query)   : {stats['aligned_bp_query']:,}\n")
        fh.write(f"Coverage of A        : {stats['coverage_a_pct']:.2f}%\n")
        fh.write(f"Coverage of B        : {stats['coverage_b_pct']:.2f}%\n")
        fh.write(f"Mean identity        : {stats['mean_identity']:.2f}%\n")
        fh.write(f"Length-weighted identity: {stats['weighted_identity']:.2f}%\n")
        fh.write(f"Min / Max identity   : {stats['min_identity']:.2f}% / {stats['max_identity']:.2f}%\n")


# ── Module 3: dot plot ─────────────────────────────────────────────────────────

def _ordered_offsets(lengths: dict) -> tuple:
    """Sort sequences by length descending, return (offsets dict, ordered list)."""
    ordered = sorted(lengths.items(), key=lambda kv: kv[1], reverse=True)
    offsets = {}
    cum = 0
    for name, length in ordered:
        offsets[name] = cum
        cum += length
    return offsets, ordered


def build_dot_plot(records: list, lens_a: dict, lens_b: dict,
                   out_path: Path, color_map: str, plot_formats: list) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize

    matplotlib.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "figure.facecolor": "white",
    })

    offsets_a, ordered_a = _ordered_offsets(lens_a)
    offsets_b, ordered_b = _ordered_offsets(lens_b)

    segments, colors = [], []
    for r in records:
        x0 = offsets_a[r["qname"]] + r["qstart"]
        x1 = offsets_a[r["qname"]] + r["qend"]
        if r["strand"] == "+":
            y0 = offsets_b[r["tname"]] + r["tstart"]
            y1 = offsets_b[r["tname"]] + r["tend"]
        else:
            y0 = offsets_b[r["tname"]] + r["tend"]
            y1 = offsets_b[r["tname"]] + r["tstart"]
        segments.append([(x0, y0), (x1, y1)])
        colors.append(r["identity"])

    fig, ax = plt.subplots(figsize=(10, 10))

    if segments:
        norm = Normalize(vmin=min(colors), vmax=100.0)
        lc = LineCollection(segments, cmap=color_map, norm=norm, linewidths=1.3)
        lc.set_array(colors)
        ax.add_collection(lc)
        cbar = fig.colorbar(lc, ax=ax, fraction=0.04, pad=0.02)
        cbar.set_label("% identity")

    total_a = sum(lens_a.values()) or 1
    total_b = sum(lens_b.values()) or 1
    ax.set_xlim(0, total_a)
    ax.set_ylim(0, total_b)

    for name, length in ordered_a[:MAX_LABELED_SEQS]:
        ax.axvline(offsets_a[name], color="#cccccc", linewidth=0.4, zorder=0)
    for name, length in ordered_b[:MAX_LABELED_SEQS]:
        ax.axhline(offsets_b[name], color="#cccccc", linewidth=0.4, zorder=0)

    ax.set_xlabel("Assembly A (query) position (bp)")
    ax.set_ylabel("Assembly B (target) position (bp)")
    ax.set_title(f"Dot plot: {len(records):,} alignment blocks")

    fig.tight_layout()
    for fmt in plot_formats:
        fmt_path = out_path.with_suffix(f".{fmt}")
        fig.savefig(fmt_path)
        _log(f"  Written: {fmt_path.name}")
    plt.close(fig)


# ── Module 4: sequence correspondence table ───────────────────────────────────

def build_correspondence_table(records: list, lens_a: dict, out_path: Path) -> None:
    """For each query sequence, report its best-matching target sequence
    (by total aligned bp) and the fraction of the query it covers."""
    by_query = {}
    for r in records:
        key = (r["qname"], r["tname"])
        by_query.setdefault(r["qname"], {}).setdefault(r["tname"], 0)
        by_query[r["qname"]][r["tname"]] += r["qend"] - r["qstart"]

    rows = []
    for qname in sorted(by_query):
        targets = by_query[qname]
        best_target, best_bp = max(targets.items(), key=lambda kv: kv[1])
        qlen = lens_a.get(qname, 0)
        coverage_pct = round(best_bp / qlen * 100, 2) if qlen else 0.0
        n_targets = len(targets)
        rows.append((qname, best_target, best_bp, coverage_pct, n_targets))

    with open(out_path, "w") as fh:
        fh.write("query_seq\tbest_match_target_seq\taligned_bp\tquery_coverage_pct\tn_target_seqs_hit\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row) + "\n")
    _log(f"  Written: {out_path.name} ({len(rows)} query sequences)")


# ── Module 5: unaligned / assembly-specific regions ───────────────────────────

def _merge_intervals(intervals: list) -> list:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def build_unaligned_regions(records: list, lens_a: dict, lens_b: dict,
                            out_a: Path, out_b: Path, min_gap: int = 500) -> None:
    """Write BED-style TSVs of regions in A and B with no alignment coverage."""
    def _write(intervals_by_seq: dict, lengths: dict, out_path: Path):
        n_regions = 0
        total_bp = 0
        with open(out_path, "w") as fh:
            fh.write("seq_id\tstart\tend\tlength\n")
            for seq_id in sorted(lengths):
                seq_len = lengths[seq_id]
                covered = _merge_intervals(intervals_by_seq.get(seq_id, []))
                cursor = 0
                for start, end in covered:
                    if start - cursor >= min_gap:
                        fh.write(f"{seq_id}\t{cursor}\t{start}\t{start - cursor}\n")
                        n_regions += 1
                        total_bp += start - cursor
                    cursor = max(cursor, end)
                if seq_len - cursor >= min_gap:
                    fh.write(f"{seq_id}\t{cursor}\t{seq_len}\t{seq_len - cursor}\n")
                    n_regions += 1
                    total_bp += seq_len - cursor
        return n_regions, total_bp

    query_intervals, target_intervals = {}, {}
    for r in records:
        query_intervals.setdefault(r["qname"], []).append((r["qstart"], r["qend"]))
        target_intervals.setdefault(r["tname"], []).append((r["tstart"], r["tend"]))

    n_a, bp_a = _write(query_intervals, lens_a, out_a)
    n_b, bp_b = _write(target_intervals, lens_b, out_b)
    _log(f"  Written: {out_a.name} ({n_a} regions, {bp_a:,} bp specific to A)")
    _log(f"  Written: {out_b.name} ({n_b} regions, {bp_b:,} bp specific to B)")


# ── Module 6: lightweight rearrangement flagging ──────────────────────────────

REARRANGEMENT_SLACK_BP = 200   # tolerance for target-coordinate order breaks


def flag_rearrangements(records: list, out_path: Path, min_block_len: int = 1000) -> int:
    """
    Heuristic, alignment-derived (not a full SV caller): for each (qname, tname)
    pair with >= 2 sizeable alignment blocks, flag likely inversions (a run of
    blocks on the minority strand flanked by majority-strand blocks) and likely
    rearrangements (blocks whose target coordinate breaks the expected monotonic
    order along the query, beyond REARRANGEMENT_SLACK_BP of tolerance).
    """
    by_pair = {}
    for r in records:
        if r["alen"] < min_block_len:
            continue
        by_pair.setdefault((r["qname"], r["tname"]), []).append(r)

    flags = []
    for (qname, tname), blocks in by_pair.items():
        if len(blocks) < 2:
            continue
        blocks = sorted(blocks, key=lambda r: r["qstart"])
        plus_bp = sum(r["alen"] for r in blocks if r["strand"] == "+")
        minus_bp = sum(r["alen"] for r in blocks if r["strand"] == "-")
        majority = "+" if plus_bp >= minus_bp else "-"

        for r in blocks:
            if r["strand"] != majority:
                flags.append(("inversion", qname, r["qstart"], r["qend"],
                             tname, r["tstart"], r["tend"]))

        majority_blocks = [r for r in blocks if r["strand"] == majority]
        for prev, curr in zip(majority_blocks, majority_blocks[1:]):
            if majority == "+":
                if curr["tstart"] < prev["tend"] - REARRANGEMENT_SLACK_BP:
                    flags.append(("rearrangement", qname, curr["qstart"], curr["qend"],
                                 tname, curr["tstart"], curr["tend"]))
            else:
                if curr["tend"] > prev["tstart"] + REARRANGEMENT_SLACK_BP:
                    flags.append(("rearrangement", qname, curr["qstart"], curr["qend"],
                                 tname, curr["tstart"], curr["tend"]))

    flags.sort(key=lambda f: (f[1], f[2]))
    with open(out_path, "w") as fh:
        fh.write("type\tquery_seq\tquery_start\tquery_end\ttarget_seq\ttarget_start\ttarget_end\n")
        for f in flags:
            fh.write("\t".join(str(v) for v in f) + "\n")
    _log(f"  Written: {out_path.name} ({len(flags)} candidate rearrangements/inversions)")
    return len(flags)


# ── Carbon footprint / resource usage ─────────────────────────────────────────

def _start_tracker(logs_dir: Path, out_file: str, project_name: str, disabled: bool):
    if disabled:
        _log("  Carbon footprint tracking disabled (--disable_co2_tracking)")
        return None
    try:
        from codecarbon import EmissionsTracker
        tracker = EmissionsTracker(
            output_dir=str(logs_dir), output_file=out_file,
            project_name=project_name, log_level="warning",
        )
        tracker.start()
        _log("  codecarbon tracker started")
        return tracker
    except ImportError:
        _log("  codecarbon not installed — carbon tracking skipped "
             "(conda install -c conda-forge codecarbon)")
        return None


# ── compare_assemblies subcommand ─────────────────────────────────────────────

def run_compare_assemblies(args) -> None:
    global _LOG_FH

    args.assembly_a = args.assembly_a.resolve()
    args.assembly_b = args.assembly_b.resolve()

    run_dir = Path(args.output)
    results = run_dir / "results"
    workdir = run_dir / "workdir"
    logs_dir = run_dir / "logs"
    for d in (results, workdir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)
    prefix = run_dir.name

    log_path = logs_dir / "Run_YithCOMPASM.log"
    _LOG_FH = open(log_path, "w")
    sep = "=" * 62
    _LOG_FH.write(f"{sep}\n  YithCOMPASM {VERSION}  —  Run Log\n{sep}\n")
    _LOG_FH.write(f"Date      : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    _LOG_FH.write(f"User      : {getpass.getuser()}\n")
    _LOG_FH.write(f"Server    : {platform.node()}\n")
    _LOG_FH.write(f"OS        : {platform.system()} {platform.release()} ({platform.machine()})\n")
    _LOG_FH.write(f"Directory : {os.getcwd()}\n")
    _LOG_FH.write(f"Command   : {' '.join(sys.argv)}\n")
    _LOG_FH.write(f"{sep}\n\n")
    _LOG_FH.flush()

    _print_quote()

    _validate_inputs([("--assembly_a", args.assembly_a), ("--assembly_b", args.assembly_b)])

    if args.force:
        _log("--force set: all steps will rerun regardless of existing outputs")
    elif workdir.exists() and any(workdir.iterdir()):
        _log("Existing workdir found — resuming from checkpoints (use --force to rerun all steps)")

    plot_formats = [f.strip().lstrip(".") for f in args.format.split(",")]

    if args.dry_run:
        _banner("Dry run — no steps will be executed")
        _log(f"  Assembly A  : {args.assembly_a}")
        _log(f"  Assembly B  : {args.assembly_b}")
        _log(f"  Output      : {run_dir}/")
        _log("  Steps that would run:")
        if not args.skip_metrics:
            _log("    [1] Metrics comparison    → results/mod01_metrics_*.tsv")
        _log("    [2] minimap2 alignment    → workdir/*.paf")
        if not args.skip_dotplot:
            _log("    [3] Dot plot              → results/mod03_dotplot_*")
        _log("    [4] Alignment summary     → results/mod04_alignment_summary_*.txt")
        if not args.skip_correspondence:
            _log("    [5] Sequence correspondence → results/mod05_correspondence_*.tsv")
        if not args.skip_unaligned:
            _log("    [6] Unaligned regions     → results/mod06_unaligned_*.tsv")
        if not args.skip_rearrangements:
            _log("    [7] Rearrangement flags   → results/mod07_rearrangements_*.tsv")
        _log("  Exiting (--dry_run).")
        if _LOG_FH:
            _LOG_FH.close()
        sys.exit(0)

    t_start = time.monotonic()
    tracker = _start_tracker(logs_dir, f"{prefix}.emissions.csv", "YithCOMPASM",
                             args.disable_co2_tracking)

    command = " ".join(sys.argv)
    summary = {
        "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "version": VERSION,
        "input_assembly_a": str(args.assembly_a),
        "input_assembly_b": str(args.assembly_b),
        "parameters": {
            "preset": args.preset, "threads": args.threads,
            "min_align_len": args.min_align_len, "min_identity": args.min_identity,
            "include_secondary": args.include_secondary, "color_map": args.color_map,
        },
    }

    if not args.skip_metrics:
        _banner("Module 1 — Assembly metrics comparison")
        metrics_path = results / f"mod01_metrics_comparison_{prefix}.tsv"
        run_metrics_comparison(args.assembly_a, args.assembly_b, metrics_path, args.force)

    _banner("Module 2 — minimap2 alignment")
    paf_path = workdir / f"{prefix}.paf"
    run_minimap2(args.assembly_a, args.assembly_b, paf_path, args.preset, args.threads, args.force)

    _log("Scanning assembly sequence lengths...")
    lens_a = scan_fasta_lengths(args.assembly_a)
    lens_b = scan_fasta_lengths(args.assembly_b)
    _log(f"  Assembly A: {len(lens_a):,} sequences, {sum(lens_a.values()):,} bp")
    _log(f"  Assembly B: {len(lens_b):,} sequences, {sum(lens_b.values()):,} bp")

    _log("Parsing PAF alignment...")
    records = parse_paf(paf_path, args.min_align_len, args.min_identity, args.include_secondary)
    _log(f"  {len(records):,} alignment blocks kept "
         f"(min_align_len={args.min_align_len}, min_identity={args.min_identity})")

    if not args.skip_dotplot:
        _banner("Module 3 — Dot plot")
        dotplot_path = results / f"mod03_dotplot_{prefix}"
        build_dot_plot(records, lens_a, lens_b, dotplot_path, args.color_map, plot_formats)

    _banner("Module 4 — Alignment summary")
    alignment_stats = compute_alignment_summary(records, lens_a, lens_b)
    summary_path = results / f"mod04_alignment_summary_{prefix}.txt"
    write_alignment_summary(summary_path, alignment_stats, args.assembly_a, args.assembly_b,
                            args.preset, command)
    _log(f"  Written: {summary_path.name}")
    _log(f"  Coverage A: {alignment_stats['coverage_a_pct']:.2f}%  "
         f"Coverage B: {alignment_stats['coverage_b_pct']:.2f}%  "
         f"Weighted identity: {alignment_stats['weighted_identity']:.2f}%")
    summary["results"] = alignment_stats

    if not args.skip_correspondence:
        _banner("Module 5 — Sequence correspondence")
        corr_path = results / f"mod05_correspondence_{prefix}.tsv"
        build_correspondence_table(records, lens_a, corr_path)

    if not args.skip_unaligned:
        _banner("Module 6 — Unaligned / assembly-specific regions")
        unaligned_a = results / f"mod06_unaligned_A_{prefix}.tsv"
        unaligned_b = results / f"mod06_unaligned_B_{prefix}.tsv"
        build_unaligned_regions(records, lens_a, lens_b, unaligned_a, unaligned_b)

    if not args.skip_rearrangements:
        _banner("Module 7 — Rearrangement / inversion flags")
        rearr_path = results / f"mod07_rearrangements_{prefix}.tsv"
        n_flags = flag_rearrangements(records, rearr_path, args.min_rearrangement_len)
        summary["n_rearrangement_flags"] = n_flags

    emissions_kg = None
    if tracker is not None:
        try:
            emissions_kg = tracker.stop()
        except Exception:
            pass

    elapsed_s = time.monotonic() - t_start
    ru = resource.getrusage(resource.RUSAGE_SELF)
    peak_mem_mb = (ru.ru_maxrss / (1024 * 1024) if platform.system() == "Darwin"
                  else ru.ru_maxrss / 1024)
    summary["resource_usage"] = {
        "wall_clock_s": round(elapsed_s, 1),
        "peak_mem_mb": round(peak_mem_mb, 1),
        "emissions_kg_CO2eq": emissions_kg,
    }

    summary_json_path = results / f"{prefix}.run_summary.json"
    with open(summary_json_path, "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    _log(f"Written: {summary_json_path.name}")

    _banner("Done")
    if _LOG_FH is not None:
        _LOG_FH.close()


# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="YithCOMPASM",
        description="Compare genome assemblies: metrics, alignment-based dot "
                    "plots, sequence correspondence, and structural comparison.",
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    sub = ap.add_subparsers(dest="command", required=True)

    cmp_ap = sub.add_parser("compare_assemblies",
                            help="Compare two FASTA assemblies")
    cmp_ap.add_argument("--assembly_a", required=True, type=Path, help="First assembly FASTA (used as query)")
    cmp_ap.add_argument("--assembly_b", required=True, type=Path, help="Second assembly FASTA (used as target)")
    cmp_ap.add_argument("--output", required=True, help="Output directory")
    cmp_ap.add_argument("--threads", type=int, default=4, help="Threads passed to minimap2 (default: 4)")
    cmp_ap.add_argument("--preset", default="asm5", choices=["asm5", "asm10", "asm20"],
                        help="minimap2 preset: asm5/asm10/asm20 ~ up to 5%%/10%%/20%% divergence (default: asm5)")
    cmp_ap.add_argument("--min_align_len", type=int, default=1000,
                        help="Minimum alignment block length (bp) to include (default: 1000)")
    cmp_ap.add_argument("--min_identity", type=float, default=0.0,
                        help="Minimum %% identity to include (default: 0.0, no filtering)")
    cmp_ap.add_argument("--include_secondary", action="store_true",
                        help="Include minimap2 secondary alignments (default: primary only)")
    cmp_ap.add_argument("--color_map", default="RdYlGn",
                        help="Matplotlib colormap name for the dot plot identity scale (default: RdYlGn)")
    cmp_ap.add_argument("--format", default="jpeg",
                        help="Comma-separated dot plot formats: jpeg, png, pdf, svg (default: jpeg)")
    cmp_ap.add_argument("--skip_metrics", action="store_true", help="Skip Module 1 — assembly metrics comparison")
    cmp_ap.add_argument("--skip_dotplot", action="store_true", help="Skip Module 3 — dot plot")
    cmp_ap.add_argument("--skip_correspondence", action="store_true", help="Skip Module 5 — sequence correspondence table")
    cmp_ap.add_argument("--skip_unaligned", action="store_true", help="Skip Module 6 — unaligned region report")
    cmp_ap.add_argument("--skip_rearrangements", action="store_true", help="Skip Module 7 — rearrangement flagging")
    cmp_ap.add_argument("--min_rearrangement_len", type=int, default=1000,
                        help="Minimum alignment block length (bp) considered for Module 7 "
                             "inversion/rearrangement flagging (default: 1000)")
    cmp_ap.add_argument("--force", action="store_true",
                        help="Rerun all steps from scratch even if intermediate outputs exist in workdir/")
    cmp_ap.add_argument("--dry_run", action="store_true",
                        help="Validate inputs and print the steps that would run, then exit without executing anything")
    cmp_ap.add_argument("--disable_co2_tracking", action="store_true",
                        help="Disable carbon footprint tracking even if codecarbon is installed")
    cmp_ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")

    args = ap.parse_args(argv)

    if args.command == "compare_assemblies":
        run_compare_assemblies(args)


if __name__ == "__main__":
    main()
