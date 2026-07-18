# Test data

`make_test_data.py` generates a small synthetic pair of assemblies
(`random.seed(42)`) with known, deliberate differences so every module can
be verified against a concrete expected outcome. Regenerate with:

```bash
python3 make_test_data.py
```

## Assemblies

**test_assembly_A.fasta** — 2 sequences, 20,000 bp total
- `ctg1` 12,000 bp
- `ctg2` 8,000 bp

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
| 1 | `mod01_metrics_comparison_*.tsv` | A: 2 seqs/20,000bp; B: 3 seqs/21,300bp |
| 2 | `workdir/*.paf` | 4 alignment records (2 for ctg1, 2 for ctg2) |
| 3 | `mod03_dotplot_*.jpeg` | Two mostly-diagonal series of dots (ctg1, ctg2), one series showing a clean anti-diagonal segment past the ctg1 breakpoint (the inversion) |
| 4 | `mod04_alignment_summary_*.txt` | Coverage A ≈96.5%, Coverage B ≈90.6%, weighted identity ≈98.4% |
| 5 | `mod05_correspondence_*.tsv` | ctg1→ctg1 (~100% coverage), ctg2→ctg2 (~91% coverage) |
| 6 | `mod06_unaligned_A_*.tsv` | One region: `ctg2  3300  4000  700` |
| 6 | `mod06_unaligned_B_*.tsv` | One region: `ctg3  0  2000  2000` |
| 7 | `mod07_rearrangements_*.tsv` | One `inversion` at `ctg1 8000-12000`; one `rearrangement` at `ctg2 4000-8000` |

**Note on Module 7 and minimap2 chaining:** the test's inversion uses a
single clean breakpoint (not an internal block sandwiched between two
colinear anchors). minimap2's `asm5`/`asm10`/`asm20` presets use a large
default chaining gap tuned for assembly-level structural tolerance — a
*sandwiched* few-kb inversion can get bridged into one degenerate
match/deletion/insertion alignment instead of being split into a separate
reverse-strand PAF record, which Module 7 would then miss (it only sees
the split records minimap2 actually reports). A single-breakpoint
inversion reliably splits into two clean records and is what this test
data uses; keep this in mind when interpreting Module 7 results on real
assemblies with small, internally-nested rearrangements.
