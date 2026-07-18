# Test data

`make_test_data.py` generates a small synthetic pair of assemblies
(`random.seed(42)`) with known, deliberate differences so every module can
be verified against a concrete expected outcome. Regenerate with:

```bash
python3 make_test_data.py
```

## Assemblies

**test_assembly_A.fasta** — 3 sequences, 28,000 bp total
- `ctg1` 12,000 bp
- `ctg2` 8,000 bp
- `ctg4` 8,000 bp — a second, independently-mutated (~2%) copy of `ctg2`,
  simulating an uncollapsed second haplotype that both align back to the
  same target region in B (see Module 7 below)

**test_assembly_B.fasta** — 3 sequences, 21,300 bp total, derived from A
- `ctg1` 12,000 bp — ~1.5% SNPs; single breakpoint at 8,000 bp with the
  second half (A:ctg1:8000-12000) reverse-complemented
- `ctg2` 7,300 bp — ~1.5% SNPs; the two halves swapped
  (A:ctg2:4000-8000 moved before A:ctg2:0-4000); the last 700 bp of the
  first half (A:ctg2:3300-4000) dropped entirely
- `ctg3` 2,000 bp — unrelated random sequence, present only in B

## Quicktest

```bash
conda activate YithCOMPASM

python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a test/test_assembly_A.fasta \
    --assembly_b test/test_assembly_B.fasta \
    --output test_run/ \
    --threads 2
```

## Expected outputs (verified 2026-07-18, minimap2 2.31)

| Module | File | Expected content |
|---|---|---|
| 1 | `mod01_metrics_comparison_*.tsv` | A: 3 seqs/28,000bp; B: 3 seqs/21,300bp |
| 2 | `workdir/*.paf` | 6 alignment records (2 each for ctg1, ctg2, ctg4) |
| 3 | `mod03_dotplot_*.jpeg` | Diagonal dot series for ctg1/ctg2/ctg4 against B, one showing a clean anti-diagonal segment past the ctg1 breakpoint (the inversion), and ctg2+ctg4 both landing on the same B:ctg2 target band (the redundancy) |
| 4 | `mod04_alignment_summary_*.txt` | Coverage A ≈94.5%, Coverage B ≈90.3%, weighted identity ≈97.9%, **redundant bp in B ≈7,246 (multiplicity ≈1.38x)** |
| 5 | `mod05_correspondence_*.tsv` | ctg1→ctg1 (~99.7%), ctg2→ctg2 (~90.8%), **ctg4→ctg2 (~90.6%)** — two different query sequences both landing on B:ctg2 is the correspondence-table signature of a collapsed target |
| 6 | `mod06_unaligned_A_*.tsv` | One region: `ctg2  3300  4000  700` |
| 6 | `mod06_unaligned_B_*.tsv` | One region: `ctg3  0  2000  2000` |
| 7 | `mod07_redundancy_A_*.tsv` | Empty (0 regions) |
| 7 | `mod07_redundancy_B_*.tsv` | Two regions on `ctg2`, depth 2, covering ≈7,246 bp total — the region hit by both `ctg2` and `ctg4` |
| 8 | `mod08_rearrangements_*.tsv` | One `inversion` on `(ctg1, ctg1)`; one `rearrangement` on `(ctg2, ctg2)`; one `rearrangement` on `(ctg4, ctg2)` |
| 9 | `mod09_identity_histogram_*.jpeg`, `.tsv` | 4 bins, bp-weighted mean identity ≈97.9% (matches Module 4) |

**Note on Module 8 and minimap2 chaining:** the test's inversion uses a
single clean breakpoint (not an internal block sandwiched between two
colinear anchors). minimap2's `asm5`/`asm10`/`asm20` presets use a large
default chaining gap tuned for assembly-level structural tolerance — a
*sandwiched* few-kb inversion can get bridged into one degenerate
match/deletion/insertion alignment instead of being split into a separate
reverse-strand PAF record, which Module 8 would then miss (it only sees
the split records minimap2 actually reports). A single-breakpoint
inversion reliably splits into two clean records and is what this test
data uses; keep this in mind when interpreting Module 8 results on real
assemblies with small, internally-nested rearrangements.

## Real-world validation

Beyond this synthetic dataset, Module 4/5/7 were validated against real
genome assemblies from a hybrid yeast genome-assembly class exercise:
comparing a Flye assembly (99 contigs, 21.5 Mb, kept both haplotype copies
largely separate) against a HiFiasm assembly of the same sample (17
contigs, 12.6 Mb, collapsed both copies into one representation) showed
`multiplicity_b ≈ 1.74x` and every one of the 16 HiFiasm contigs receiving
2-14 best-matching Flye contigs in the Module 5 correspondence table — a
clean, independent confirmation that the redundancy metrics correctly
detect real haplotype-collapse events, not just the synthetic construction
above.
