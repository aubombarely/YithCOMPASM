# Changelog

## [v0.3.0] — 2026-07-18

### Added
- `--preset` now also accepts `map-hifi`/`map-ont`/`map-pb`, alongside the
  existing `asm5`/`asm10`/`asm20`. Discovered while searching real
  assemblies against a small (19-86 kb) mitochondrial genome reference
  panel: the `asm*` presets are tuned for large, near-collinear
  genome-vs-genome alignment and silently miss real hits when
  `--assembly_b` is much smaller than `--assembly_a` (e.g. a single
  gene/organelle panel) — `minimap2 -x asm20` found zero alignments where
  `-x map-ont` correctly found real (if short) matches.

## [v0.2.0] — 2026-07-18

### Added
- Module 7 — redundancy / coverage-depth report (`mod07_redundancy_{A,B}_*.tsv`):
  a sweep-line coverage-depth analysis flagging every region covered ≥2x by
  alignment blocks, for both assemblies. Rearrangement flagging renumbered
  from Module 7 to Module 8 (`mod08_rearrangements_*.tsv`) to make room.
- `redundant_bp_a`/`redundant_bp_b` and `multiplicity_a`/`multiplicity_b`
  fields in the alignment summary and run summary JSON.
- `--skip_redundancy` flag.
- Test dataset extended with a fourth query sequence (`ctg4`, a divergent
  duplicate of `ctg2`) to exercise the new redundancy detection.

### Fixed
- `coverage_a_pct`/`coverage_b_pct` could previously exceed 100% when a
  target region was hit by multiple query alignments (raw base-pair sums,
  not deduplicated) — discovered while validating against real assembly
  data (a haplotype-collapsing HiFiasm assembly vs. a haplotype-preserving
  Flye assembly of the same hybrid yeast sample, which legitimately
  produced `coverage_b_pct = 169.52%`). Coverage is now computed from
  merged, deduplicated alignment intervals and is correctly bounded to
  [0, 100]; the redundancy this previously (accidentally) signaled is now
  reported explicitly instead, in more detail, via the fields and module
  above.

## [v0.1.0] — 2026-07-18

### Added
- Initial release: `compare_assemblies` subcommand.
- Module 1 — assembly metrics comparison (N50/L50, GC%, sequence counts, deltas).
- Module 2 — minimap2 whole-genome alignment (`-c` for exact identity, checkpointed).
- Module 3 — whole-genome dot plot, alignment blocks colored by % identity.
- Module 4 — coverage/identity summary statistics.
- Module 5 — per-sequence best-match correspondence table.
- Module 6 — unaligned/assembly-specific region report for both assemblies.
- Module 7 — lightweight, alignment-derived inversion/rearrangement flagging.
- Synthetic test dataset (`test/make_test_data.py`, `random.seed(42)`)
  exercising every module with verified expected output.
