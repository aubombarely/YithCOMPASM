<img src="https://img.shields.io/badge/version-v0.1.0-teal"/>
<img src="https://img.shields.io/badge/python-3.10%2B-blue"/>
<img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey"/>

See [CHANGELOG.md](CHANGELOG.md) for release history.

# YithCOMPASM

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
| 3 | `mod03_dotplot_*.jpeg` | Whole-genome dot plot, alignment blocks colored by % identity |
| 4 | `mod04_alignment_summary_*.txt` | Coverage and identity statistics |
| 5 | `mod05_correspondence_*.tsv` | Best-matching target sequence per query sequence |
| 6 | `mod06_unaligned_{A,B}_*.tsv` | Regions with no alignment coverage (assembly-specific sequence) |
| 7 | `mod07_rearrangements_*.tsv` | Alignment-derived inversion/rearrangement candidates (heuristic, not a full SV caller) |

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
| `--preset` | No | minimap2 preset: `asm5`/`asm10`/`asm20` ~ up to 5%/10%/20% divergence (default: `asm5`) |
| `--min_align_len` | No | Minimum alignment block length (bp) to include (default: 1000) |
| `--min_identity` | No | Minimum % identity to include (default: 0.0) |
| `--include_secondary` | No | Include minimap2 secondary alignments (default: primary only) |
| `--color_map` | No | Matplotlib colormap for the dot plot identity scale (default: `RdYlGn`) |
| `--format` | No | Comma-separated dot plot formats: jpeg, png, pdf, svg (default: `jpeg`) |
| `--min_rearrangement_len` | No | Minimum block length considered for Module 7 flagging (default: 1000) |
| `--skip_metrics` | No | Skip Module 1 |
| `--skip_dotplot` | No | Skip Module 3 |
| `--skip_correspondence` | No | Skip Module 5 |
| `--skip_unaligned` | No | Skip Module 6 |
| `--skip_rearrangements` | No | Skip Module 7 |
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
│   ├── mod04_alignment_summary_{prefix}.txt
│   ├── mod05_correspondence_{prefix}.tsv
│   ├── mod06_unaligned_A_{prefix}.tsv
│   ├── mod06_unaligned_B_{prefix}.tsv
│   ├── mod07_rearrangements_{prefix}.tsv
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

## Rearrangement/inversion flagging — what it is and isn't

Module 7 is a lightweight heuristic derived directly from the alignment
already computed for the dot plot: for each query/target sequence pair with
multiple alignment blocks, it flags blocks on the minority strand
(candidate inversions) and blocks that break the expected monotonic target-
coordinate order along the query (candidate rearrangements). It is **not**
a structural variant caller (no SyRI/AsmSV-style breakpoint refinement) and
it can only flag what minimap2 actually reports as separate alignment
records — a small inversion sandwiched between two colinear anchors can get
bridged by minimap2 into one degenerate alignment instead of being split,
in which case Module 7 will not see it. See `test/README.md` for a worked
example of this behavior.

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
4. Cross-check anything visually striking against `mod07_rearrangements_*.tsv`
   and `mod06_unaligned_*.tsv` for a tabular breakdown.
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
