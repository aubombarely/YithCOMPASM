# YithCOMPASM

<p align="center">
  <img src="assets/yithcompasm_logo.svg" width="260" alt="YithCOMPASM logo"/>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/version-v0.4.1-teal" alt="Version v0.4.1"/>
  <img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+"/>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey" alt="Platform"/>
</p>

<p align="center">
  <a href="CHANGELOG.md">Changelog</a>
</p>

Genome assembly comparison toolkit. Named for the **Great Race of Yith**
(H.P. Lovecraft, *The Shadow Out of Time*) — beings whose defining trait is
documenting and comparing civilizations and eras across time — combined
with **COMPASM** (Compare Assemblies), echoing a compass for navigating
between two genomes.

## Overview

| Module | Output | Description |
|---|---|---|
| 1 | `mod01_metrics_comparison_*.tsv` | Side-by-side assembly metrics (N50, GC%, sequence counts...) with deltas |
| 2 | `workdir/*.paf` | minimap2 whole-genome alignment (checkpointed) |
| 3 | `mod03_dotplot_*.jpeg`, `mod03_dotplot_*.html` | Whole-genome dot plot (static image + zoomable/interactive HTML), alignment blocks colored by % identity |
| 4 | `mod04_alignment_summary_*.txt` | Coverage and identity statistics |
| 5 | `mod05_correspondence_*.tsv` | Best-matching target sequence per query sequence |
| 6 | `mod06_unaligned_{A,B}_*.tsv` | Regions with no alignment coverage (assembly-specific sequence) |
| 7 | `mod07_redundancy_{A,B}_*.tsv` | Regions covered ≥2x (e.g. a collapsed haplotype hit by multiple contigs of the other assembly) |
| 8 | `mod08_rearrangements_*.tsv` | Alignment-derived inversion/rearrangement candidates (heuristic, not a full SV caller) |

## Requirements

### Conda installation

```bash
conda env create -f envs/YithCOMPASM.yaml
conda activate YithCOMPASM
```

### Dependencies

| Tool | Purpose |
|---|---|
| Python ≥ 3.10 | pipeline orchestration |
| minimap2 | whole-genome pairwise alignment |
| matplotlib | dot plot rendering |
| codecarbon | carbon footprint tracking *(optional)* |

## Usage

```bash
python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a genomeA.fasta \
    --assembly_b genomeB.fasta \
    --output results_run/ \
    --threads 8
```

## Options

| Argument | Required | Description |
|---|---|---|
| `--assembly_a` | Yes | First assembly FASTA (used as query) |
| `--assembly_b` | Yes | Second assembly FASTA (used as target) |
| `--output` | Yes | Output directory |
| `--threads` | No | Threads passed to minimap2 (default: 4) |
| `--preset` | No | minimap2 preset: `asm5`/`asm10`/`asm20` (up to 5%/10%/20% divergence, large near-collinear genome-vs-genome) or `map-hifi`/`map-ont`/`map-pb` (use when `--assembly_b` is a small reference, e.g. a gene/organelle panel — the `asm*` presets are not sensitive enough for that; default: `asm5`) |
| `--min_align_len` | No | Minimum alignment block length (bp) to include (default: 1000) |
| `--min_identity` | No | Minimum % identity to include (default: 0.0) |
| `--include_secondary` | No | Include minimap2 secondary alignments (default: primary only) |
| `--color_map` | No | Matplotlib colormap for the dot plot identity scale (default: `RdYlGn`) |
| `--format` | No | Comma-separated dot plot formats: jpeg, png, pdf, svg (default: `jpeg`) |
| `--min_rearrangement_len` | No | Minimum block length considered for Module 8 flagging (default: 1000) |
| `--skip_metrics` | No | Skip Module 1 |
| `--skip_dotplot` | No | Skip Module 3 (both the static image and the interactive HTML) |
| `--skip_dotplot_html` | No | Keep the static dot plot image but skip generating the interactive HTML version |
| `--skip_correspondence` | No | Skip Module 5 |
| `--skip_unaligned` | No | Skip Module 6 |
| `--skip_redundancy` | No | Skip Module 7 |
| `--skip_rearrangements` | No | Skip Module 8 |
| `--force` | No | Rerun all steps even if checkpointed outputs exist |
| `--dry_run` | No | Validate inputs and print steps without executing |
| `--disable_co2_tracking` | No | Disable carbon footprint tracking |
| `--version` | No | Show version and exit |

## Output directory layout

```
{output}/
├── results/
│   ├── mod01_metrics_comparison_{prefix}.tsv
│   ├── mod03_dotplot_{prefix}.jpeg
│   ├── mod03_dotplot_{prefix}.html   (interactive: zoom/pan/hover/filter/basket)
│   ├── mod04_alignment_summary_{prefix}.txt
│   ├── mod05_correspondence_{prefix}.tsv
│   ├── mod06_unaligned_A_{prefix}.tsv
│   ├── mod06_unaligned_B_{prefix}.tsv
│   ├── mod07_redundancy_A_{prefix}.tsv
│   ├── mod07_redundancy_B_{prefix}.tsv
│   ├── mod08_rearrangements_{prefix}.tsv
│   └── {prefix}.run_summary.json
├── workdir/
│   └── {prefix}.paf              (safe to delete after a successful run)
└── logs/
    ├── Run_YithCOMPASM.log
    └── {prefix}.emissions.csv    (if codecarbon is installed)
```

## Percent identity and the `-c` flag

Alignment is run as `minimap2 -x {preset} -c ...`. The `-c` flag forces
base-level CIGAR generation — without it, minimap2's default PAF `nmatch`/
`alen` columns are approximate minimizer-chain estimates, not exact
alignment identity, and will visibly understate real similarity (verified
during development: omitting `-c` on a 1.5%-divergence test pair reported
~55-78% identity instead of the correct ~98.4%).

## Interactive dot plot

`mod03_dotplot_{prefix}.html` is a self-contained interactive companion to
the static image — vanilla JS + SVG, no CDN dependencies, opens directly in
any browser, works fully offline. It exists because static dot plots can't
label every sequence ID for fragmented assemblies (dozens to hundreds of
contigs) without becoming unreadable. Features:

- **Zoom / pan**: scroll to zoom, drag to pan.
- **Hover** an alignment line for its exact query/target coordinates,
  identity, and strand.
- **Filter** by minimum sequence length and minimum % identity, and by an
  adjustable SeqID label size (8-60px) so labels are legible at whatever
  zoom level and screen you're using.
- **Search** by sequence ID to jump straight to it.
- **Basket**: click an alignment line to add both its query and target
  sequence IDs to a basket shown in the side panel; each entry is
  individually removable.
- **Export**: download the basket as a plain ID list (`.txt`), or load the
  original assembly FASTA files locally in the browser (nothing is
  uploaded anywhere) and download just the basketed sequences as a ready-to-
  use FASTA — e.g. for a follow-up BLAST of a region flagged as redundant or
  rearranged.

See the [UserCase01 interactive dot plot](examples/SPSC01_HiFiasm_vs_Flye/dotplot_interactive.html)
for a real example (open it locally after cloning the repo).

## Coverage, redundancy, and collapsed haplotypes

`coverage_a_pct`/`coverage_b_pct` are computed from **deduplicated**
(interval-merged) alignment coverage — each base of an assembly is counted
at most once, so these percentages can never exceed 100%. This matters
specifically when comparing two assemblies of the same sample built with
different collapsing behavior (e.g. a haplotype-collapsing assembler like
early-mode HiFiasm vs. a haplotype-preserving one like Flye): if the naive
sum of aligned base-pairs were used instead, a target region hit by two
different query contigs (one per haplotype) would be counted twice,
inflating coverage past 100% — which is exactly what an earlier version of
this tool did before this was caught by testing against real HiFiasm/Flye
assemblies of a hybrid yeast strain.

That redundancy signal is genuinely useful, so it isn't just discarded —
it's reported explicitly:
- `redundant_bp_a`/`redundant_bp_b` and `multiplicity_a`/`multiplicity_b`
  in `mod04_alignment_summary_*.txt` and the run summary JSON give a
  whole-assembly-level number (e.g. `multiplicity_b: 1.74` means assembly B
  is, on average, hit 1.74x over its aligned length — consistent with a
  substantial fraction of it being present as two separate, only partially
  merged haplotype copies in assembly A).
- **Module 7** (`mod07_redundancy_{A,B}_*.tsv`) gives the region-level
  breakdown: a sweep-line coverage-depth report listing every region
  covered ≥2x, with its exact depth. This is what actually lets you find
  *which* contigs/regions are collapsed, not just that some redundancy
  exists somewhere.

## Rearrangement/inversion flagging — what it is and isn't

Module 8 is a lightweight heuristic derived directly from the alignment
already computed for the dot plot: for each query/target sequence pair with
multiple alignment blocks, it flags blocks on the minority strand
(candidate inversions) and blocks that break the expected monotonic target-
coordinate order along the query (candidate rearrangements). It is **not**
a structural variant caller (no SyRI/AsmSV-style breakpoint refinement) and
it can only flag what minimap2 actually reports as separate alignment
records — a small inversion sandwiched between two colinear anchors can get
bridged by minimap2 into one degenerate alignment instead of being split,
in which case Module 8 will not see it. See `test/README.md` for a worked
example of this behavior.

## Use cases

- [UserCase01 — Is the HiFiasm assembly "really duplicated"?](UserCase01_SPSC01_hybrid_yeast_HiFiasm_vs_Flye.md)
  Real HiFiasm vs. Flye assemblies of a hybrid yeast strain (SPSC01,
  *S. cerevisiae* × *S. pombe* protoplast fusant). Uses the redundancy and
  correspondence modules to explain a 70%-duplicated-BUSCO discrepancy
  between the two assemblers with base-pair-resolved evidence.

## Run log / run summary / carbon footprint

Every run writes `logs/Run_YithCOMPASM.log`, `results/{prefix}.run_summary.json`
(parameters, alignment statistics, wall-clock time, peak memory, and
estimated CO2 emissions if `codecarbon` is installed), and an optional
`logs/{prefix}.emissions.csv`.

## Recommended workflow

1. Run `compare_assemblies` with default settings first.
2. Inspect `mod04_alignment_summary_*.txt` for overall coverage/identity.
3. Open `mod03_dotplot_*.jpeg` to visually scan for large-scale structural
   differences.
4. Cross-check anything visually striking against `mod08_rearrangements_*.tsv`,
   `mod07_redundancy_*.tsv`, and `mod06_unaligned_*.tsv` for a tabular breakdown.
5. Use `mod05_correspondence_*.tsv` to confirm scaffold-to-chromosome
   assignment, e.g. as input to `GenoToolBoxPlus/FastaRename.py` if you need
   to harmonize sequence naming between the two assemblies.

## Future functionality

Additional subcommands are planned under the same `YithCOMPASM.py` entry
point (e.g. `identify_syntenic_blocks`), and a `compare_simplerepeats_landscapes`
subcommand leveraging UbboTELORNA's telomere/repeat output is under
consideration.

## Quicktest

```bash
conda activate YithCOMPASM

python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a test/test_assembly_A.fasta \
    --assembly_b test/test_assembly_B.fasta \
    --output test_run/ \
    --threads 2
```

See `test/README.md` for the exact expected output of this run.

## Third-party tools and citations

- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
  *Bioinformatics*, 34:3094-3100.
- Hunter, J. D. (2007). Matplotlib: A 2D graphics environment.
  *Computing in Science & Engineering*, 9:90-95.
