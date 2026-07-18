#!/usr/bin/env python3
"""
make_test_data.py — generate a small synthetic pair of assemblies for
YithCOMPASM testing (random.seed(42) for reproducibility).

Assembly A (test_assembly_A.fasta):
    ctg1  12,000 bp
    ctg2   8,000 bp

Assembly B (test_assembly_B.fasta), derived from A with known, deliberate
differences so the tool's modules have something concrete to detect:
    ctg1  ~1.5% SNPs throughout; single breakpoint at 8,000 bp, with the
          second half (originally A:ctg1:8000-12000) reverse-complemented
          -> expect one Module 7 "inversion" flag on (ctg1, ctg1). A single,
             non-sandwiched breakpoint is used deliberately: minimap2's
             asm* presets have a large default chaining gap and will bridge
             a *sandwiched* internal inversion into one degenerate
             match/deletion/insertion block instead of splitting it into a
             separate reverse-strand PAF record. A single clean breakpoint
             is what reliably produces two separate, cleanly split records.
    ctg2  ~1.5% SNPs throughout; the two halves swapped (A:ctg2:4000-8000
          moved before A:ctg2:0-4000); the last 700 bp of the first half
          (A:ctg2:3300-4000) dropped entirely
          -> expect exactly one Module 7 "rearrangement" flag, and one
             Module 6 "unaligned in A" region around ctg2:3300-4000
    ctg3  2,000 bp of unrelated random sequence, present only in B
          -> expect a Module 6 "unaligned in B" region covering ~all of ctg3
"""

import random
from pathlib import Path

random.seed(42)

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def random_seq(length: int) -> str:
    return "".join(random.choice(BASES) for _ in range(length))


def mutate(seq: str, rate: float) -> str:
    seq = list(seq)
    for i in range(len(seq)):
        if random.random() < rate:
            seq[i] = random.choice([b for b in BASES if b != seq[i]])
    return "".join(seq)


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def write_fasta(path: Path, records: dict, width: int = 70) -> None:
    with open(path, "w") as fh:
        for seq_id, seq in records.items():
            fh.write(f">{seq_id}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def main():
    out_dir = Path(__file__).parent

    # Assembly A
    ctg1_a = random_seq(12000)
    ctg2_a = random_seq(8000)
    write_fasta(out_dir / "test_assembly_A.fasta", {"ctg1": ctg1_a, "ctg2": ctg2_a})

    # Assembly B — ctg1: SNPs + one inverted block after a single breakpoint at 8000bp
    s1, s2 = ctg1_a[0:8000], ctg1_a[8000:12000]
    ctg1_b = s1 + revcomp(s2)
    ctg1_b = mutate(ctg1_b, 0.015)

    # Assembly B — ctg2: two halves swapped, last 700bp of the first half dropped
    t1, t2 = ctg2_a[0:3300], ctg2_a[4000:8000]   # t1 is A:ctg2:0-3300 (3300-4000 dropped)
    ctg2_b = t2 + t1
    ctg2_b = mutate(ctg2_b, 0.015)

    # Assembly B — ctg3: unrelated sequence, private to B
    ctg3_b = random_seq(2000)

    write_fasta(out_dir / "test_assembly_B.fasta",
               {"ctg1": ctg1_b, "ctg2": ctg2_b, "ctg3": ctg3_b})

    print("Written test_assembly_A.fasta (ctg1=12000bp, ctg2=8000bp)")
    print("Written test_assembly_B.fasta (ctg1=12000bp w/ inversion, "
         "ctg2=7300bp w/ swap+drop, ctg3=2000bp private)")


if __name__ == "__main__":
    main()
