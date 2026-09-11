#!/bin/bash

# Stage 1: haplotype frequency spectra (and H12/H2H1) in sliding windows.
# --calc-spec is what writes the .spectra.gz file the second command reads;
# without it this stage writes .stats.gz instead. --lassi belongs to stage 2
# and has no effect here (lassip now warns if you pass it).
lassip --vcf YRI.chr22.vcf.gz --pop YRI.ids.pop.txt --calc-spec --hapstats \
       --k 10 --winsize 117 --winstep 12 --out YRI.chr22

# Stage 2: the LASSI likelihood ratio, reading the spectra written above.
# Pass every contig of a population in one command; --salti runs saltiLASSI
# instead, and --threads parallelises over windows.
lassip --spectra YRI.chr22.lassip.hap.spectra.gz --lassi --out YRI.chr22
