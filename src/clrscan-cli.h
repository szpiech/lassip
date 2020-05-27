#ifndef __CLRSCAN_CLI_H__
#define __CLRSCAN_CLI_H__

const string VERSION = "0.0.1";

const string PREAMBLE = "";

const string ARG_THREADS = "--threads";
const int DEFAULT_THREADS = 1;
const string HELP_THREADS = "The number of threads to spawn during computations.";

// I/O flags
//const string ARG_FILENAME_TPED = "--tped";
//const string DEFAULT_FILENAME_TPED = "__hapfile1";
//const string HELP_FILENAME_TPED = "A TPED file containing haplotype and map data.\n\
//\tVariants should be coded 0/1";

const string ARG_FILENAME_POP1_VCF = "--vcf";
const string DEFAULT_FILENAME_POP1_VCF = "__vcffile1";
const string HELP_FILENAME_POP1_VCF = "A VCF file containing haplotype data.\n\
\tVariants should be coded 0/1";

const string ARG_FILENAME_POPFILE = "--pop";
const string DEFAULT_FILENAME_POPFILE = "__popfile1";
const string HELP_FILENAME_POPFILE = "A file containing <ind ID> <pop ID>.";

const string ARG_FILENAME_SPECFILES = "--spectra";
const string DEFAULT_FILENAME_SPECFILES = "__specfile1";
const string HELP_FILENAME_SPECFILES = "A list of spectra files for finalization.";

/*
const string ARG_FILENAME_MAP = "--map";
const string DEFAULT_FILENAME_MAP = "__mapfile";
const string HELP_FILENAME_MAP = "A mapfile with one row per variant site.\n\
\tFormatted <chr#> <locusID> <genetic pos> <physical pos>.";
*/
const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The basename for all output files.";

// Window control flags
/*
const string ARG_BP = "--bp";
const bool DEFAULT_BP = false;
const string HELP_BP = "Use bps for window sizes.";

const string ARG_SITES = "--sites";
const bool DEFAULT_SITES = false;
const string HELP_SITES = "Use sites for window sizes.";
*/

const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 0;
const string HELP_WINSIZE = "The window size within which to calculate statistics.";

const string ARG_WINSTEP = "--winstep";
const int DEFAULT_WINSTEP = 0;
const string HELP_WINSTEP = "The sliding window step size.";

/*
const string ARG_PARTITION = "--partition";
const int DEFAULT_PARTITION = 0;
const string HELP_PARTITION = "Partition the sliding window into non-overlapping sub windows\n\
\tof varying sizes within which all statistics (except EHH-based ones) are calculated separately.\n\
\te.g. For a sliding window of 100kb, --partition 25000 50000 25000 would instruct\n\
\tdivstats to calculate statistics separately within a central 50kb and within the\n\
\ttwo flanking 25kb regions for each 100kb window. Partitions must add up to --winsize.\n\
\tSet to 0 to simply calculate within the entire window.";
const int MAX_PARTITION = 20;
*/
// Statistics flags
const string ARG_LASSI = "--lassi";
const bool DEFAULT_LASSI = false;
const string HELP_LASSI = "Set this flag to use the LASSI method.";

const string ARG_LASSI_CHOICE = "--lassi-choice";
const int DEFAULT_LASSI_CHOICE = 4;
const string HELP_LASSI_CHOICE = "Set this flag to change the way LASSI\n\
\tdistributes mass across sweeping haplotype classes. Take an integer in {1..5}.";

const string ARG_SFINDER = "--sweepfinder";
const bool DEFAULT_SFINDER = false;
const string HELP_SFINDER = "Set this flag to use the SweepFinder method.";

const string ARG_SFINDER2 = "--sweepfinder2";
const bool DEFAULT_SFINDER2 = false;
const string HELP_SFINDER2 = "Set this flag to use the SweepFinder2 method.";
/*
const string ARG_PI = "--pi";
const bool DEFAULT_PI = false;
const string HELP_PI = "Set this flag to calculate mean pairwise sequence difference.";

const string ARG_PIK = "--pik";
const int DEFAULT_PIK = 0;
const string HELP_PIK = "Set this flag to calculate mean pairwise sequence difference amongst \n\
\tthe k most frequent haplotypes. You can choose more than one, e.g. --pik 2 3 4 will\n\
\tcalculate pi amongst the top 2, 3, and 4 most frequent haplotypes.  If set to 0, does\n\
\tnot calculate.";

const string ARG_SEGSITES = "--s";
const bool DEFAULT_SEGSITES = false;
const string HELP_SEGSITES = "Set this flag to calculate the number of segregating sites.";

const string ARG_EHH = "--ehh";
const int DEFAULT_EHH = 0;
const string HELP_EHH = "A list of window sizes within which to calculate EHH.\n\
\tThese subwindows will be centered on the middle of the current\n\
\twindow and may not be larger than --winsize.\n";

const string ARG_EHHK = "--ehhk";
const int DEFAULT_EHHK = 0;
const string HELP_EHHK = "Calculates EHH, after collapsing\n\
\tthe k most frequent haplotypes into a single identity class using the windows\n\
\tdefined by --ehh. This flag requires --ehh to be set.\n\
\tIf set to 0 does not calculate.";

const string ARG_TAJ_D = "--d";
const bool DEFAULT_TAJ_D = false;
const string HELP_TAJ_D = "Set this flag to calculate Tajima's D.";

const string ARG_FAY_WU_H = "--h";
const bool DEFAULT_FAY_WU_H = false;
const string HELP_FAY_WU_H = "Set this flag to calculate Fay and Wu's H.";
*/
// Other flags
const string ARG_INIT = "--initial";
const bool DEFAULT_INIT = false;
const string HELP_INIT = "Initial computations.";

const string ARG_K = "--k";
const int DEFAULT_K = 10;
const string HELP_K = "Top K haplotypes for LASSI computations.";

const string ARG_FINALIZE = "--finalize";
const bool DEFAULT_FINALIZE = false;
const string HELP_FINALIZE = "Finalize computations.";
/*
const string ARG_EHH_PART = "--ehh-part";
const bool DEFAULT_EHH_PART = false;
const string HELP_EHH_PART = "Calculates EHH/EHHK in any partitions of the main window.\n\
\tTo be distinguished from --ehh, which calculates EHH in sub-windows\n\
\tcentered on the main window.  Requires --ehh/--ehhk to be set.";

const string ARG_NO_SFS_SUB = "--no-sfs-sub";
const bool DEFAULT_NO_SFS_SUB = false;
const string HELP_NO_SFS_SUB = "Do not subsample the SFS to handle missing data.\n\
\tEffectively treats missing data as 0/0 in the SFS. Can speed up computation of\n\
sfs-based statistics substantially.";

const string ARG_PMAP = "--pmap";
const bool DEFAULT_PMAP = false;
const string HELP_PMAP = "Use physical map instead of a genetic map.";
*/

//#define NOPTS 7

//const string STATS[NOPTS] = {ARG_PI, ARG_PIK, ARG_SEGSITES, ARG_EHH, ARG_EHHK, ARG_TAJ_D, ARG_FAY_WU_H};

#endif