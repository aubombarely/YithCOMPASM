# Changelog

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
