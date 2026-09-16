#!/bin/bash

# Stage 1: haplotype frequency spectra (and H12/H2H1) in sliding windows.
# --calc-spec is what writes the .spectra.gz file the second command reads;
# without it this stage writes .stats.gz instead. --lassi belongs to stage 2
# and has no effect here (lassip now warns if you pass it).
lassip --vcf YRI.chr22.vcf.gz --pop YRI.ids.pop.txt --calc-spec --hapstats \
       --k 10 --winsize 117 --winstep 12 --out YRI.chr22

# Several contigs go in one stage-1 run, either as one file per contig or as
# files holding several, and the run writes a single spectra file per
# population covering all of them. The example ships only chr22, so this form
# is shown rather than run:
#
#   lassip --vcf chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz --pop YRI.ids.pop.txt \
#          --calc-spec --hapstats --k 10 --winsize 117 --winstep 12 --out YRI
#   lassip --spectra YRI.YRI.lassip.hap.spectra.gz --salti --out YRI
#
# Passing the per-contig spectra files separately to stage 2 gives the same
# answer; what you must not do is concatenate them by hand, since the rows of
# one contig then sit inside another's block.

# Stage 2: the LASSI likelihood ratio, reading the spectra written above.
# At the default --filter-level 2 each population is filtered separately and
# gets its own file, named <out>.<pop>.lassip.hap.spectra.gz -- hence the YRI
# in the middle. Pass every contig of a population in one command; --salti
# runs saltiLASSI instead, and --threads parallelises over windows.
lassip --spectra YRI.chr22.YRI.lassip.hap.spectra.gz --lassi --out YRI.chr22
