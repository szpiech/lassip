#!/bin/bash

lassip --vcf YRI.chr22.vcf.gz --hapstats --lassi --winsize 117 --winstep 12 --out YRI.chr22 --pop YRI.ids.pop.txt
lassip --spectra YRI.chr22.lassip.hap.spectra.gz --lassi --out YRI.chr22
