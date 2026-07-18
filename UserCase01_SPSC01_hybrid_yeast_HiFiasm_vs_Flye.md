# Use Case 01 — Is the HiFiasm assembly "really duplicated"? Comparing HiFiasm and Flye assemblies of a hybrid yeast strain

**Goal:** Use `YithCOMPASM compare_assemblies` to explain a striking BUSCO
discrepancy observed between two PacBio HiFi assemblies of the same sample —
built with HiFiasm and Flye — and determine whether it reflects a technical
assembly artifact or genuine biology.

This example originates from the *Escuela de Verano BCB 2026* genome
assembly practical (Exercise 5–6, `P01_AssemblyPractice.pdf`). The class QC
table (BUSCO, Merqury) flags a large duplication discrepancy between the two
assemblers and explicitly asks: *"Is the HIFIASM really duplicated? Is a
technical problem or do it have a biological reason?"* — this use case
answers that question with base-pair-resolved, whole-genome evidence rather
than gene-sampling metrics alone.

---

## Overview

| Item | Detail |
|---|---|
| Sample | SPSC01 — an industrial hybrid yeast strain created by protoplast fusion of *Saccharomyces cerevisiae* and *Schizosaccharomyces pombe* |
| BioProject | [PRJNA1075684](https://www.ncbi.nlm.nih.gov/bioproject/1075684) — "*Saccharomyces cerevisiae* SPSC01 Genome sequencing and assembly" |
| SRA run | SRR27947616 (PacBio HiFi) |
| Read filtering | Chopper, Q20, minimum length 30 kb |
| Assemblers compared | HiFiasm (primary contigs) vs. Flye |
| SPSC01 published assembly | GCA_047651925.1 (17 chromosomes) |
| Parental references used (Step 6) | *S. cerevisiae* S288C (GCF_000146045.2), *S. pombe* 972h⁻ (GCF_000002945.2) |

---

## Background — why this sample is a good stress test

SPSC01 is not a simple diploid strain: it is a laboratory-made **interspecies
fusant** combining an industrial *S. cerevisiae* strain with *S. pombe*,
selected for improved flocculation and stress tolerance. This matters for
interpreting assembly output — *S. cerevisiae* and *S. pombe* are only
distantly related (different genera), so any duplication signal here is
*not* explained by the usual "two haplotypes of one diploid genome" story
in the same simple way it would be for an intraspecific cross. Whatever the
assemblies disagree about needs to be checked directly rather than assumed.

---

## Step 1 — Data source

PacBio HiFi reads for SPSC01 were downloaded from SRA run **SRR27947616**,
associated with BioProject
[PRJNA1075684](https://www.ncbi.nlm.nih.gov/bioproject/1075684).

---

## Step 2 — Read filtering with Chopper

Raw HiFi reads were filtered with Chopper (quality ≥ Q20, minimum length
30 kb) before assembly:

```bash
chopper -t 3 --quality 20 --minlength 30000 \
    --input SRR27947616_PBHiFi_SPSC01.fastq \
    > SRR27947616_PBHiFi_SPSC01.fQ20l30k.fastq
```

---

## Step 3 — Assembly

Both assemblers were run on the same filtered HiFi read set.

### HiFiasm

```bash
hifiasm -t 3 -o SacerSPSC01HIFIASM01 \
    SRR27947616_PBHiFi_SPSC01.fQ20l30k.fastq

# primary-contig GFA converted to FASTA with GenoToolBoxPlus
python3 scripts/GFA2FASTA.py \
    --input   SacerSPSC01HIFIASM01.bp.p_ctg.gfa \
    --output  SacerSPSC01HIFIASM01.bp.p_ctg.fasta \
    --summary SacerSPSC01HIFIASM01.bp.p_ctg.summary.txt
```

### Flye

```bash
flye --genome-size 12m --threads 3 \
    --out-dir SacerFLYE01_HIFI \
    --pacbio-hifi SRR27947616_PBHiFi_SPSC01.fQ20l30k.fastq
```

Note the `--genome-size 12m` estimate — sized for a single *S. cerevisiae*-like
haploid genome, not the ~24.7 Mb a full *S. cerevisiae* + *S. pombe* fusant
would need if both parental genomes were present intact. This mainly affects
Flye's internal coverage heuristics rather than acting as a hard cap, but is
worth keeping in mind given the results below.

### Basic assembly metrics (QUAST, from the class practical)

| Metric | Flye | HiFiasm |
|---|---|---|
| Total assembly size (Mb) | 21.49 | 12.62 |
| Total contigs (N) | 99 | 17 |
| Longest contig (Mb) | 1.21 | 1.53 |
| N90 (Mb) | 0.13 | 0.94 |
| L90 (N) | 65 | 13 |

---

## Step 4 — The BUSCO discrepancy that motivates this use case

BUSCO (`saccharomycetes_odb10`, 2137 orthologs) was run on both assemblies
as part of the class QC exercise. The exact output (`short_summary.*.txt`)
shows:

| Metric | Flye | HiFiasm |
|---|---|---|
| Complete BUSCOs (C) | 99.2% | 98.7% |
| **Single-copy (S)** | **28.9%** (618 genes) | **95.6%** (2044 genes) |
| **Duplicated (D)** | **70.3%** (1502 genes) | **3.0%** (65 genes) |
| Fragmented (F) | 0.2% | 0.0% |
| Missing (M) | 0.6% | 1.3% |

This is the crux of the class's Question 3. 70% of near-universal
single-copy orthologs appear in duplicate in the Flye assembly — a huge
signal — but BUSCO only samples ~2,137 loci and can't say *where* in the
genome this is happening, *how much* sequence is affected in total, or
whether it's isolated to specific regions/contigs. A same-assembly
duplication scan (e.g. Mash-based, as used in `YuggASMoth`) also doesn't
resolve this, because it only flags redundant sequence *within* a single
assembly — and neither assembly is internally self-redundant here (each
one is an internally consistent representation; running YuggASMoth's
duplication check on the HiFiasm assembly alone reports none, which is
expected and not actually informative about this question). Answering the
class's question needs a direct, base-pair-resolved comparison *between*
the two assemblies — which is exactly what `YithCOMPASM compare_assemblies`
does.

---

## Step 5 — Compare the assemblies with YithCOMPASM

```bash
python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a Sacer_HIFI_FLYE.fasta \
    --assembly_b Sacer_HIFI_HIFIASM.fasta \
    --output     results/SPSC01_flye_vs_hifiasm \
    --preset     asm20 \
    --threads    4
```

`asm20` was used instead of the default `asm5` because the divergence
between the two parental subgenomes of an interspecies fusant is expected
to exceed the ~5% `asm5` is tuned for.

---

## Results

| Metric | Value |
|---|---|
| Alignment blocks kept | 136 |
| Coverage of Flye (query) | 99.28% |
| Coverage of HiFiasm (target) | 97.23% |
| Length-weighted identity | 97.98% (range 72.5–100%) |
| **Redundant bp in HiFiasm** | **9,124,531 bp** |
| **Multiplicity in HiFiasm** | **1.74x** |
| Redundant regions (depth ≥2) in HiFiasm | 132 (48 of which are depth ≥3) |
| Rearrangement/inversion candidates | 3 |

The `multiplicity_b = 1.74x` figure lines up almost exactly with the raw
assembly-size ratio (21.49 Mb / 12.62 Mb ≈ 1.70) — independent, whole-genome
confirmation of the same signal BUSCO's gene sampling picked up, but now
quantified in base pairs and localized to 132 explicit regions
(`examples/SPSC01_HiFiasm_vs_Flye/redundancy_HIFIASM.tsv`) instead of an
aggregate percentage over ~2,100 genes.

The correspondence table (`mod05_correspondence_*.tsv`) shows the same
pattern from a different angle: every one of HiFiasm's 17 contigs is the
best match for somewhere between 2 and 14 different Flye contigs — the
signature of many query sequences collapsing onto one target, consistent
with Flye resolving two divergent copies of most loci that HiFiasm's
primary-contig output merges into one.

### Dot plot

![SPSC01 Flye vs HiFiasm dot plot](examples/SPSC01_HiFiasm_vs_Flye/dotplot.png)

Flye (query, x-axis) vs. HiFiasm (target, y-axis), alignment blocks colored
by percent identity. The repeated near-parallel diagonal bands — multiple
Flye contigs tracking the same y-range on a single HiFiasm contig — are the
visual signature of the redundancy quantified above.

An [interactive version of this plot](examples/SPSC01_HiFiasm_vs_Flye/dotplot_interactive.html)
is also generated alongside the static image (`mod03_dotplot_*.html`, open it
directly in any browser — no server needed). It's the more useful view for
this dataset specifically because Flye's 99 contigs are far too many to
label on a static image: zoom in to read individual sequence IDs, filter by
minimum sequence length or %identity to declutter the fragmented tail of
Flye's assembly, click an alignment line to add both its Flye and HiFiasm
sequence to a "basket", and export the basketed sequences straight to FASTA
(loading the two assembly files locally in the browser) for follow-up, e.g.
a targeted BLAST of one of the redundant regions.

---

## Step 6 — Resolving the origin of the redundancy: comparison against both parental references

The comparison above establishes *that* HiFiasm's assembly is redundant
relative to Flye's, but not *which* parental subgenome (*S. cerevisiae* or
*S. pombe*) the retained sequence belongs to, or whether the "hybrid"
nature of SPSC01 is actually reflected in the sequenced genome content.
Both assemblies were compared against both parental reference genomes:

```bash
# S. cerevisiae reference: GCF_000146045.2 (S288C)
# S. pombe reference:      GCF_000002945.2 (972h-)

python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a Sacer_HIFI_FLYE.fasta \
    --assembly_b GCF_000146045.2_R64_genomic.fna \
    --output results/flye_vs_scer --preset asm20 --threads 4

python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a Sacer_HIFI_FLYE.fasta \
    --assembly_b GCF_000002945.2_ASM294v3_genomic.fna \
    --output results/flye_vs_spombe --preset asm20 --threads 4

python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a Sacer_HIFI_HIFIASM.fasta \
    --assembly_b GCF_000146045.2_R64_genomic.fna \
    --output results/hifiasm_vs_scer --preset asm20 --threads 4

python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a Sacer_HIFI_HIFIASM.fasta \
    --assembly_b GCF_000002945.2_ASM294v3_genomic.fna \
    --output results/hifiasm_vs_spombe --preset asm20 --threads 4
```

### Results

| Comparison | Coverage of assembly | Coverage of reference | Multiplicity |
|---|---|---|---|
| Flye vs. *S. cerevisiae* (S288C) | **98.38%** | 97.98% | **1.79x** |
| Flye vs. *S. pombe* (972h⁻) | 0.18% | 0.10% | — |
| HiFiasm vs. *S. cerevisiae* (S288C) | **96.23%** | 96.73% | **1.04x** |
| HiFiasm vs. *S. pombe* (972h⁻) | 1.45% | 0.08% | — |

This is decisive. Essentially the entire Flye assembly (98%) and HiFiasm
assembly (96%) is *S. cerevisiae* sequence; almost nothing aligns to
*S. pombe* (the handful of hits are a few kb at most, most plausibly
conserved rRNA/tRNA loci rather than real synteny). The internal
consistency check confirms it numerically, not just qualitatively: the
*S. cerevisiae* reference is 12.16 Mb; 12.16 × 1.79 ≈ 21.77 Mb (Flye's real
total is 21.49 Mb), and 12.16 × 1.04 ≈ 12.65 Mb (HiFiasm's real total is
12.62 Mb) — both match almost exactly.

![Flye vs S. cerevisiae reference dot plot](examples/SPSC01_HiFiasm_vs_Flye/dotplot_flye_vs_Scerevisiae_ref.png)
![HiFiasm vs S. cerevisiae reference dot plot](examples/SPSC01_HiFiasm_vs_Flye/dotplot_hifiasm_vs_Scerevisiae_ref.png)

---

## Discussion — answering the class's Question 3

**"Is the HIFIASM really duplicated? Is a technical problem or do it have a
biological reason?"**

The `YithCOMPASM` result supports a biological explanation over a technical
error, with one important caveat worth investigating further:

1. **The redundancy is real, large, and structured, not noise.** 9.12 Mb —
   nearly three-quarters of the HiFiasm assembly's own 12.62 Mb length — is
   covered ≥2x by Flye sequence, concentrated in 132 discrete regions rather
   than being smeared thinly and randomly across the genome. A technical
   assembly bug would be far less likely to produce this consistent,
   genome-scale, roughly-2x pattern.

2. **This is consistent with expected assembler behavior, not a bug in
   either tool.** HiFiasm's primary-contig output (`bp.p_ctg`, what
   `GFA2FASTA.py` converted here) is designed to collapse heterozygous
   variation into one representative haploid mosaic per locus by default —
   it is not attempting haplotype-resolved output unless run with
   additional phasing information (trio binning, Hi-C). Flye, run here
   without any haplotype-collapsing step, is retaining both divergent
   copies as separate contigs. Seen this way, "HiFiasm is duplicated" is
   the wrong framing — rather, **Flye is retaining real biological
   variation that HiFiasm's primary assembly mode is designed to discard.**
   This directly matches why a same-assembly duplication scan (YuggASMoth)
   reports nothing wrong with the HiFiasm assembly on its own: there
   genuinely isn't redundant sequence *within* it — the "duplication" only
   becomes visible when compared against an assembly that made a different
   representational choice.

3. **The redundancy is *S. cerevisiae* haplotype duplication, not retained
   *S. pombe* content — and the "hybrid" genome isn't really in these
   assemblies at all.** Step 6 resolves what was initially left as an open
   question: both assemblies are >96% *S. cerevisiae* sequence, and
   *S. pombe* is essentially absent (<2% coverage, driven by a handful of
   short conserved loci rather than real synteny). So the redundancy isn't
   "Flye kept both parental subgenomes separate, HiFiasm merged them" — it
   is specifically **two divergent haplotype copies of the *S. cerevisiae*
   subgenome**, retained separately by Flye (1.79x vs. the pure
   *S. cerevisiae* reference) and collapsed to essentially one copy by
   HiFiasm (1.04x). The *S. pombe* side of this engineered fusant strain is
   not represented in the HiFi sequencing data used for either assembly —
   consistent with the well-documented genomic instability of protoplast
   fusants, where one parental genome is frequently lost over generations
   of industrial strain propagation and selection. Confirming *why* it was
   lost (loss during fusant construction vs. loss during subsequent culture,
   vs. simply not present in the DNA extraction used for this sequencing
   run) would need information outside the scope of what these two
   assemblies alone can answer.

**Practical takeaway for the class exercise:** don't discard or "fix" the
HiFiasm duplication as an error — it's the expected behavior of a
primary/collapsed assembly mode applied to a genuinely heterozygous
diploid *S. cerevisiae* subgenome. If a haplotype-resolved assembly is the
actual goal, that requires rerunning HiFiasm with phasing data (trio,
Hi-C) rather than treating this primary-contig output as final, or working
from the Flye assembly directly with the understanding that it is *not*
single-copy per locus. Separately, if recovering the *S. pombe* half of
this strain's genome is actually the goal, that is a sample/sequencing
question (is *S. pombe* DNA present in the culture at all?), not an
assembly-parameter question — no amount of reassembly of this particular
read set will recover subgenome content that isn't in the reads.

---

## Addendum — hunting for a hidden organism via candidate contigs and BLAST

Step 6 shows *S. pombe* is essentially absent by whole-genome alignment,
but alignment against a *specific* reference only tests "does this match
*this* organism" — it can't identify an unexpected one. A cheap way to
probe for that: mitochondria exist at very high copy number per cell, so
even if a second organism's *nuclear* genome were too low-abundance to
assemble, its mitochondrial genome could still assemble cleanly as a
small, distinct contig. Two candidates stood out:

- **`contig_147`** (Flye, 45 kb) — only 50% explained by the
  *S. cerevisiae* reference in the Step 6 correspondence table, the
  worst-explained contig in either assembly.
- **`ptg000017l`** (HiFiasm, 92.7 kb) — the smallest HiFiasm contig, close
  to *S. cerevisiae* mitochondrial size (~85.8 kb); notably, **no contig in
  either assembly has the real *S. cerevisiae* mitochondrion
  (`NC_001224.1`) as its best match at all**, which is itself a small open
  question.

Both were submitted to NCBI's remote BLAST (`blastn`, `nt`, megablast) via
the CGI URL API, since no local BLAST+ installation was available:

```bash
# submit
curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi" \
    --data-urlencode "CMD=Put" \
    --data-urlencode "PROGRAM=blastn" \
    --data-urlencode "MEGABLAST=on" \
    --data-urlencode "DATABASE=nt" \
    --data-urlencode "QUERY=$(grep -v '^>' contig_147.fasta | tr -d '\n')"
# poll (no more than once/minute) using the returned RID, then fetch results:
curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=<RID>&FORMAT_TYPE=Text&DESCRIPTIONS=15"
```

**Result: both are ordinary *S. cerevisiae*, not a hidden organism.**

- `contig_147` top hits: the *AQY1* aquaporin gene region and several
  *S. cerevisiae* strain chromosomes (IV, XII, XVI, VII), 98–100% identity.
  Aquaporin loci and subtelomeres are known to be paralogous and
  strain-variable — this explains why a single-reference (S288C) alignment
  couldn't cleanly assign it to one chromosome, without needing a foreign
  organism to explain it.
- `ptg000017l` top hit: *S. cerevisiae* chromosome XII at 99% identity
  (very high bitscore), plus several 35S rRNA / rDNA-adjacent hits — this
  region carries the large tandem rDNA repeat array (~100–200 copies),
  a classic repeat-collapse assembly failure mode. Most likely a
  mis-assembled rDNA-array fragment, not foreign DNA. (A secondary
  *S. paradoxus* hit at 94% identity is expected background conservation
  in that region across the genus, not a distinct signal — the
  *S. cerevisiae* hit is both higher-identity and far higher-scoring.)

**This reinforces, rather than overturns, the Step 6 conclusion.** The two
strongest candidates for "unexplained sequence that might belong to a
different organism" both turned out to have mundane, well-known
explanations specific to *S. cerevisiae* genome assembly difficulty
(paralogous subtelomeres, rDNA repeat collapse). There is no BLAST
evidence anywhere in either assembly of *S. pombe* or a third organism.
This still left one loose end open: why does no contig in either assembly
have the real *S. cerevisiae* mitochondrion as its best match? The next
section answers that directly.

### Follow-up — mapping a mitogenome reference panel directly

Rather than guessing which contig might be a mitochondrial genome,
a more direct test: download known mitochondrial genomes for several
relevant species and align them against *both whole assemblies* at once.
A hit would show up regardless of which contig it landed in or how large
that contig is.

| Species | Accession | Length |
|---|---|---|
| *S. cerevisiae* | NC_001224.1 | 85,779 bp |
| *S. pombe* | NC_088682.1 | 19,433 bp |
| *S. paradoxus* | CP125617.1 | 73,893 bp |
| *S. uvarum* | OP499835.1 | 64,015 bp |
| *S. japonicus* | NC_004332.1 | 80,059 bp |
| *S. octosporus* | NC_004312.1 | 44,227 bp |

```bash
python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a Sacer_HIFI_FLYE.fasta \
    --assembly_b mitogenome_panel.fasta \
    --output results/flye_vs_mitopanel \
    --preset map-ont --min_align_len 30 --threads 4
```

**Important methodological catch:** the first attempt at this used
`--preset asm20` (matching the rest of this use case) and found **zero**
alignments — not filtered out, minimap2 itself reported nothing, even
against *S. cerevisiae*'s own mitochondrion. This turned out to be the
wrong tool setting, not a biological result: the `asm*` presets are tuned
for large, near-collinear genome-vs-genome alignment and are not sensitive
enough when the reference is tiny (here, 19-86 kb) relative to a
multi-megabase assembly. Switching to `--preset map-ont` (added to
YithCOMPASM specifically because of this) recovered real alignments.

**With the right preset, the result is decisive and different from what
this document previously speculated:**

| | Flye | HiFiasm |
|---|---|---|
| Alignment blocks | 12 | 6 |
| Aligned bp (deduplicated) | 1,089 | 478 |
| Coverage of mitopanel (367,406 bp total) | 0.24% | 0.13% |
| Longest single hit | 169 bp | 93 bp |

Every hit is a short (46-169 bp) fragment, and the same fragment on a
given nuclear contig frequently matches *multiple different species'*
mitogenomes at similar length and identity — the signature of conserved
tRNA-like elements or ancient nuclear-embedded mitochondrial fragments
(NUMTs), not a real assembled mitochondrial genome. There is no long,
contiguous, high-identity block anywhere, for any of the six species
tested.

**Corrected conclusion:** mitochondrial DNA is not meaningfully assembled
in either genome — not merged into a larger contig as earlier speculated,
essentially absent. Given the nuclear genome is unambiguously ~96-98%
*S. cerevisiae* (Step 6), this is most plausibly a library-prep/size-selection
artifact — PacBio HiFi SMRTbell library preparation is typically optimized
around the bulk nuclear fragment size distribution, and small circular
mitochondrial molecules often do not survive that process the same way —
rather than anything to do with *S. pombe* or a hidden organism.

### Follow-up — a divergent haplotype pair confirmed by BLAST (contig_127 / contig_128)

While exploring the redundancy pattern in `YithCOMPASM`'s interactive dot
plot (`examples/SPSC01_HiFiasm_vs_Flye/dotplot_interactive.html`), two Flye
contigs stood out: `contig_127` and `contig_128` both map onto the same
HiFiasm contig, `ptg000015l` — a candidate for "maybe one of this pair is
actually the missing *S. pombe* subgenome, collapsed by HiFiasm onto a
*S. cerevisiae* contig by mistake." The alignment identity already computed
by `YithCOMPASM` argued against that before any BLAST was run:

| Query | Aligned to `ptg000015l` | Identity |
|---|---|---|
| `contig_127` | 328,455 bp | 100.00% |
| `contig_128` | 246,346 bp + 89,904 bp | 93.29% / 91.57% |

*S. cerevisiae* and *S. pombe* are ~300-450 million years diverged and share
only ~50-70% nucleotide identity even in conserved coding regions — nowhere
near what minimap2 (even at the permissive `asm20`, ≤20% divergence) would
align across hundreds of kb. 91-93% identity is squarely within
same-species/divergent-haplotype range, not cross-species.

BLAST confirmed it directly. A first attempt submitting the full ~330 kb
contigs to NCBI's remote `blastn`/`nt`/megablast (same CGI method as the
addendum above) failed server-side (`Error: Process size limit exceeded,
SIGXFSZ`) — a near-100%-identity, 330 kb query matches so many near-identical
*S. cerevisiae* strain assemblies in the redundant `nt` database that the
result itself is too large for NCBI to build. Since only species-level
identification was needed (not a full structural mapping — `YithCOMPASM`
already has that), a representative 5 kb slice from the middle of each
contig was submitted instead:

| Contig slice | Top hit | Identity | Bitscore |
|---|---|---|---|
| `contig_127` | `LN907791.1` *S. cerevisiae* N85, chromosome VIII | 100% | 9201 |
| `contig_128` | `LN907791.1` *S. cerevisiae* N85, chromosome VIII | 100% | 9217 |

Both slices land on the *same* *S. cerevisiae* chromosome VIII locus, with
essentially tied scores. Every hit past the top ~7 near-100%-identity
*S. cerevisiae* strain assemblies falls into other *Saccharomyces sensu
stricto* species (*S. paradoxus* 91-92%, *S. mikatae* 88%, *S. kudriavzevii*
86%, *S. eubayanus* 85%) — ordinary genus-level background conservation at a
moderately conserved locus, not a hybrid-origin signal.

**Confirms, rather than overturns, Step 6 and the addendum above:**
`contig_127`/`contig_128` are two divergent *S. cerevisiae* haplotypes of the
same chromosome VIII locus — exactly the kind of pair Flye keeps separate
and HiFiasm's primary-contig mode collapses onto one (`ptg000015l`) — not a
second species hiding inside a "duplicated" HiFiasm contig.

### Follow-up — cross-checking against a GenomeScope2 k-mer estimate

A GenomeScope2 run on the same HiFi reads gave an independent, k-mer-based
estimate: genome size 11.74 Mb, heterozygosity 0.37%. Comparing that against
what `YithCOMPASM` measured directly from the alignment:

**Genome size — good agreement.**

| | Size | vs. GenomeScope2 (11.74 Mb) |
|---|---|---|
| GenomeScope2 estimate | 11.74 Mb | — |
| Real *S. cerevisiae* S288C reference | 12.16 Mb | +3.6% |
| HiFiasm collapsed assembly | 12.62 Mb | +7.5% |
| Flye assembly | 21.49 Mb | 1.83× (not the full 2× a completely haplotype-separated diploid would give) |

HiFiasm's primary-contig output approximates one representative copy per
locus, and 12.62 Mb landing within single-digit percent of both the k-mer
estimate and the real reference is a solid match — the small excess is the
same residual under-collapsing already quantified above (HiFiasm's own
multiplicity, 1.04x). Flye's 1.83× (rather than a clean 2×) foreshadows the
next point: haplotype separation isn't uniform across the genome.

**Heterozygosity — consistent for most of the genome, but a single number
hides real structure.** This is what Module 9 (`mod09_identity_histogram_*`)
is for: a bp-weighted histogram of alignment identity across all 136
Flye-vs-HiFiasm blocks, generated with:

```bash
python3 scripts/YithCOMPASM.py compare_assemblies \
    --assembly_a Sacer_HIFI_FLYE.fasta \
    --assembly_b Sacer_HIFI_HIFIASM.fasta \
    --output     results/SPSC01_flye_vs_hifiasm \
    --preset     asm20 --threads 4
```

![Identity distribution: Flye vs HiFiasm](examples/SPSC01_HiFiasm_vs_Flye/identity_histogram.jpeg)

| Identity | % of aligned bp |
|---|---|
| ≥99.5% (≤0.5% divergence) | **58.4%** |
| 99.0-99.5% | 5.6% |
| 98.0-99.0% | 8.5% |
| 95.0-98.0% | 16.8% |
| 90.0-95.0% | 6.6% |
| <90.0% | 4.2% |

(Full per-1%-bin breakdown: `examples/SPSC01_HiFiasm_vs_Flye/identity_histogram.tsv`.)

The dominant mass — 58.4% of aligned bp — sits at ≥99.5% identity, the same
regime as GenomeScope2's 0.37% heterozygosity call. That agreement is also
*why* Flye doesn't reach a clean 2× assembly size: where haplotype
divergence is genuinely this low, there's essentially no k-mer/read signal
for either assembler to phase on, low-heterozygosity or not.

But **41.6% of the aligned genome is meaningfully more divergent than a
uniform 0.37% would predict** — including 4.2% below 90% identity, which is
exactly the range `contig_128` fell into in the follow-up above, and matches
what the earlier BLAST addendum already flagged as rDNA-array and
paralogous-subtelomere regions. GenomeScope2's %het is a single genome-wide
average from a k-mer coverage mixture model tuned to a roughly uniform SNP
signal — it's known to underrepresent exactly this kind of localized,
structurally complex divergence, since those k-mers tend to get absorbed
into the model's "repetitive/error" component rather than cleanly
contributing to the heterozygosity peak.

**Conclusion: the two estimates agree where they're measuring the same
thing** (genome size, and the bulk near-identical fraction of the genome)
**but answer different questions.** GenomeScope2 gives one number for the
whole genome; `YithCOMPASM`'s alignment-based identity distribution shows
that number is really an average masking a genome that's mostly
near-identical between haplotypes with a substantial, much more divergent
minority concentrated in repetitive/complex regions — useful context for
deciding how much to trust a single heterozygosity percentage when planning
downstream analysis (e.g. variant calling parameters, or how aggressively to
expect haplotype-aware assembly to actually separate the two copies).

---

## Overall conclusion — hybrid strain, or just a heterozygous *S. cerevisiae*?

Every independent test run against this specific HiFi dataset — whole-genome
alignment to both parental references (Step 6), a targeted BLAST hunt for
the two least-explained contigs in either assembly (Addendum), a
six-species mitogenome reference panel mapped directly against both
assemblies (Follow-up), and now BLAST confirmation of the specific
haplotype pair that first looked like it could be a hidden second genome
(this section) — converges on the same answer: **there is no detectable
*S. pombe* sequence, nuclear or mitochondrial, anywhere in either assembly.**
What *is* present and abundant is heterozygosity: two divergent haplotypes
of the *S. cerevisiae* subgenome that Flye keeps as separate contigs and
HiFiasm's primary-contig mode collapses into one, which is the entire
redundancy pattern this use case set out to explain.

That is a statement about *this sequencing run and these two assemblies*,
not necessarily about the strain's actual history. Two explanations are
consistent with everything observed here, and sequence data alone can't
distinguish between them:

- **The strain genuinely was built by *S. cerevisiae* × *S. pombe*
  protoplast fusion, as described, but has since lost the *S. pombe*
  subgenome.** Genome instability and preferential loss of one parental
  genome is a well-documented outcome of interspecies protoplast fusants,
  especially across a divergence as large as *Saccharomyces*/
  *Schizosaccharomyces* — industrial strain propagation and selection over
  many generations commonly erodes or eliminates the less-adapted parental
  genome entirely. Under this explanation, "hybrid" correctly describes the
  strain's *construction history*, but the genome actually present in the
  cells that were sequenced for this dataset is no longer hybrid.
- **The isolate sequenced for this dataset is simply a heterozygous/
  divergent-haplotype *S. cerevisiae***, and the "hybrid" label reflects
  strain lineage/documentation that doesn't precisely describe the specific
  culture or DNA extraction used here (e.g. a related but non-hybrid
  sub-isolate, a labeling carryover from an earlier stage of strain
  development, or fusion that never produced stable *S. pombe* retention in
  the first place).

Both scenarios produce an identical genome to the one observed. Resolving
which actually happened would need information this dataset can't supply
on its own — strain provenance records, karyotyping, or ideally a fresh DNA
extraction sequenced and checked the same way. **What can be said with
confidence from the data in hand: treat this specific assembly as a
heterozygous diploid *S. cerevisiae* genome, not a two-subgenome hybrid
genome** — any downstream analysis (variant calling, haplotype-resolved
assembly, annotation) should be planned accordingly, regardless of which
historical explanation is correct.

---

## Reproducing this example

All commands above use the exact assembly files and parameters run for
this example. Reference genomes: *S. cerevisiae* S288C
([GCF_000146045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/))
and *S. pombe* 972h⁻
([GCF_000002945.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002945.2/)),
both fetched via the NCBI Datasets API. The six-species mitogenome panel
(fetched via NCBI E-utilities `efetch`) is committed as
[`mitogenome_panel.fasta`](examples/SPSC01_HiFiasm_vs_Flye/mitogenome_panel.fasta).
Output files are committed under
[`examples/SPSC01_HiFiasm_vs_Flye/`](examples/SPSC01_HiFiasm_vs_Flye/):
dot plots and the SPSC01-pair summary/redundancy/rearrangement files at
the top level, and all six reference-comparison alignment summaries
(*S. cerevisiae*, *S. pombe*, mitogenome panel × both assemblers) under
`reference_comparisons/` — enough to check every number in this document
without rerunning anything.
