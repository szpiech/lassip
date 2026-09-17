/* lassip -- a program to calculate haplotype frequency spectrum statistics
   Copyright (C) 2020  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/

#ifndef __LASSIP_CLI_H__
#define __LASSIP_CLI_H__

const string VERSION = "1.3.1";

const string USAGE = "\
Usage: lassip --vcf <file> --pop <file> --calc-spec [--hapstats] --winsize <int> --winstep <int> --out <prefix>\n\
       lassip --spectra <file> [<file> ...] (--lassi | --salti | --avg-spec) --out <prefix>";

const string PREAMBLE = "\n\
lassip computes haplotype frequency spectrum statistics in two stages.\n\
\n\
Stage 1 reads genotypes and writes a haplotype frequency spectrum per window:\n\
\n\
  lassip --vcf YRI.chr22.vcf.gz --pop YRI.ids.pop.txt --calc-spec --hapstats \\\n\
         --winsize 117 --winstep 12 --k 10 --out YRI.chr22\n\
\n\
Stage 2 reads those spectra (all contigs of a population at once) and writes\n\
the likelihood ratio statistics:\n\
\n\
  lassip --spectra YRI.chr*.lassip.hap.spectra.gz --salti --out YRI\n\
\n\
Flags belonging to the other stage are ignored, so check which stage you are in\n\
if an option seems to have no effect.\n\
\n\
Methods: saltiLASSI (DeGiorgio and Szpiech 2022, PLoS Genetics 18:e1010134),\n\
LASSI (Harris and DeGiorgio 2020, MBE doi.org/10.1093/molbev/msaa115),\n\
H12 and H2/H1 (Garud et al. 2015, PLoS Genetics 11:e1005004),\n\
G123 and G2/G1, reported in place of H12 and H2/H1 under --unphased\n\
(Harris et al. 2018, Genetics 210:1429-1452).";

const string ARG_THREADS = "--threads";
const int DEFAULT_THREADS = 1;
const string HELP_THREADS = "The number of threads to spawn during computations.";

// I/O flags

const string ARG_FILENAME_MAP = "--map";
const string DEFAULT_FILENAME_MAP = "";
const string HELP_FILENAME_MAP = "A map file formatted <chr#> <locusID> <genetic pos> <physical pos>.\n\
\tSites in VCF not in map file will be interpolated.";

const string ARG_FILENAME_POP1_VCF = "--vcf";
const string DEFAULT_FILENAME_POP1_VCF = "";
const string HELP_FILENAME_POP1_VCF = "One or more VCF files containing haplotype data.\n\
\tVariants should be coded 0/1 and every file must contain all of the samples\n\
\tnamed by --pop. Several contigs may be analysed in one run, given either as\n\
\tone file per contig or as files holding several; the run then writes a single\n\
\tspectra file per population covering all of them. A contig may not appear\n\
\ttwice, and a file's records must be grouped by contig.";

const string ARG_FILENAME_POPFILE = "--pop";
const string DEFAULT_FILENAME_POPFILE = "";
const string HELP_FILENAME_POPFILE = "A file containing <ind ID> <pop ID>.";

const string ARG_FILENAME_SPECFILES = "--spectra";
const string DEFAULT_FILENAME_SPECFILES = "";
const string HELP_FILENAME_SPECFILES = "A list of spectra files for finalization.\n\
\tContigs are delimited by the chr column rather than by the file, so a file may\n\
\thold any number of them and the results are the same either way. A file's rows\n\
\tmust be grouped by contig and ascending within one, and no contig may appear\n\
\ttwice across the files given here.";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The basename for all output files.";

// Window control flags

const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 0;
const string HELP_WINSIZE = "The window size within which to calculate statistics.";

const string ARG_WINSTEP = "--winstep";
const int DEFAULT_WINSTEP = 0;
const string HELP_WINSTEP = "The sliding window step size.";

// Statistics flags

const string ARG_CALC_SPEC = "--calc-spec";
const bool DEFAULT_CALC_SPEC = false;
const string HELP_CALC_SPEC = "Set this flag to compute K-truncated haplotype\n\
frequency spectra.";

const string ARG_AVG_SPEC = "--avg-spec";
const bool DEFAULT_AVG_SPEC = false;
const string HELP_AVG_SPEC = "Set this flag to compute and output the average\n\
K-truncated haplotype frequency spectrum from a set of .spectra files.";

const string ARG_NULL_SPEC = "--null-spec";
const string DEFAULT_NULL_SPEC = "";
const string HELP_NULL_SPEC = "A file containing a null K-truncated\n\
haplotype spectrum for use computing LASSI or saltiLASSI.";

const string ARG_LASSI = "--lassi";
const bool DEFAULT_LASSI = false;
const string HELP_LASSI = "Set this flag to use the LASSI method.";

const string ARG_LASSI_CHOICE = "--lassi-choice";
const int DEFAULT_LASSI_CHOICE = 4;
const string HELP_LASSI_CHOICE = "Set this flag to change the way LASSI\n\
\tdistributes mass across sweeping haplotype classes. Takes an integer in {1..5}.";

const string ARG_HAPSTATS = "--hapstats";
const bool DEFAULT_HAPSTATS = false;
const string HELP_HAPSTATS = "Set this flag to calculate haplotype statistics.";

const string ARG_SALTI = "--salti";
const bool DEFAULT_SALTI = false;
const string HELP_SALTI = "Set this flag to use the saltiLASSI method.";


// Other flags

const string ARG_K = "--k";
const int DEFAULT_K = 10;
const string HELP_K = "Top K haplotypes for LASSI computations.";

const string ARG_FILTER_LEVEL = "--filter-level";
const int DEFAULT_FILTER_LEVEL = 2;
const string HELP_FILTER_LEVEL = "Filter monomorphic sites and sites\nwith missing data: 0-no filtering, 1-compute freq for all samples, 2-compute freq per pop.";

const string ARG_KEEP_MONO = "--keep-monomorphic";
const bool DEFAULT_KEEP_MONO = false;
const string HELP_KEEP_MONO = "Set this flag to retain any monomorphic sites in the data.";

const string ARG_FILTER_LMISS = "--max-lmiss";
const double DEFAULT_FILTER_LMISS = 0.1;
const string HELP_FILTER_LMISS = "Filter loci with > this proportion of missing data.";

const string ARG_FILTER_HMISS = "--max-hmiss";
const double DEFAULT_FILTER_HMISS = 0.2;
const string HELP_FILTER_HMISS = "Drop haplotypes with > this proportion of missing data when computing the HFS.";

const string HAP_CLUSTER_GARUD = "garud-shuffle";
const string HAP_CLUSTER_BESTCOMP = "best-comp";
const string HAP_CLUSTER_SOFTEM = "soft-em";

const string ARG_HAP_CLUSTER = "--hap-cluster";
const string DEFAULT_HAP_CLUSTER = HAP_CLUSTER_BESTCOMP;
const string HELP_HAP_CLUSTER = "How haplotypes carrying missing genotypes are grouped into\n\
classes. Missing sites act as wildcards, so such a haplotype can be compatible with\n\
several distinct haplotypes and something has to choose. One of:\n\
  best-comp      (default) each haplotype joins the most frequent class it is\n\
                 compatible with, taking the best-observed haplotypes first.\n\
                 Deterministic; --seed has no effect on it.\n\
  garud-shuffle  the pre-1.3 behaviour: shuffle the haplotypes, then let each in\n\
                 turn absorb every compatible one not yet claimed. The result\n\
                 depends on the shuffle, so runs differ unless --seed is fixed.\n\
  soft-em        EXPERIMENTAL. Divide an ambiguous haplotype\'s count across the\n\
                 classes it could belong to in proportion to their frequencies,\n\
                 iterated to convergence. Class sizes become fractional.\n\
Without missing genotypes and at --match-tol 0 all three give the same classes.";

const string ARG_MATCH_TOL = "--match-tol";
const int DEFAULT_MATCH_TOL = 0;
const string HELP_MATCH_TOL = "Group haplotypes with missing data into\nthe same class as a haplotype with no missing data if they have <= this many pairwise differences.";

const string ARG_DIST_TYPE = "--dist-type";
const string DEFAULT_DIST_TYPE = "bp";
const string HELP_DIST_TYPE = "Distance measure for saltiLASSI: bp, cm, nw.";

const string ARG_MAX_GAP = "--max-gap";
const double DEFAULT_MAX_GAP = 3000000;
const string HELP_MAX_GAP = "Widest gap in basepairs between two genetic map positions that\n\
--dist-type cm will interpolate across. A window falling in a wider gap cannot be\n\
placed and lassip stops rather than guess its genetic position. Real maps have\n\
gaps at centromeres and other low-recombination regions, so raise this if your\n\
map is sparse where your data are dense. Set 0 to interpolate across any gap.";

const string ARG_MAX_EXTEND_BP = "--max-extend-bp";
const double DEFAULT_MAX_EXTEND_BP = 100000;
const string HELP_MAX_EXTEND_BP = "Maximum distance in basepairs from core window to consider for saltiLASSI.";

const string ARG_MAX_EXTEND_CM = "--max-extend-cm";
const double DEFAULT_MAX_EXTEND_CM = 0.05;
const string HELP_MAX_EXTEND_CM = "Maximum distance in centimorgans from core window to consider for saltiLASSI.";

const string ARG_MAX_EXTEND_NW = "--max-extend-nw";
const double DEFAULT_MAX_EXTEND_NW = 5;
const string HELP_MAX_EXTEND_NW = "Maximum distance in number of windows from core window to consider for saltiLASSI.";

const string ARG_SEED = "--seed";
const int DEFAULT_SEED = 1;
const string HELP_SEED = "Seed for the shuffle used when clustering haplotypes that\n\
carry missing data. Runs with the same seed, inputs and flags are reproducible\n\
regardless of --threads. Set 0 to seed from the system clock instead, which\n\
reproduces the non-reproducible behaviour of lassip <= 1.2.2.";

const string ARG_UNPHASED = "--unphased";
const bool DEFAULT_UNPHASED = false;
const string HELP_UNPHASED = "Set this flag to indicate data are unphased.";


//#define NOPTS 7

//const string STATS[NOPTS] = {ARG_PI, ARG_PIK, ARG_SEGSITES, ARG_EHH, ARG_EHHK, ARG_TAJ_D, ARG_FAY_WU_H};

#endif