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

const string VERSION = "1.0.1";

const string PREAMBLE = "";

const string ARG_THREADS = "--threads";
const int DEFAULT_THREADS = 1;
const string HELP_THREADS = "The number of threads to spawn during computations.";

// I/O flags

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

// Other flags
/*
const string ARG_INIT = "--initial";
const bool DEFAULT_INIT = false;
const string HELP_INIT = "Initial computations.";
*/
const string ARG_K = "--k";
const int DEFAULT_K = 10;
const string HELP_K = "Top K haplotypes for LASSI computations.";

const string ARG_FILTER_LEVEL = "--filter-level";
const int DEFAULT_FILTER_LEVEL = 0;
const string HELP_FILTER_LEVEL = "Filter level: 0-none, 1-poly in analysis, 2-poly in each pop.";

/*
const string ARG_FINALIZE = "--finalize";
const bool DEFAULT_FINALIZE = false;
const string HELP_FINALIZE = "Finalize computations.";
*/

const string ARG_UNPHASED = "--unphased";
const bool DEFAULT_UNPHASED = false;
const string HELP_UNPHASED = "Set this flag to indicate data are unphased.";


//#define NOPTS 7

//const string STATS[NOPTS] = {ARG_PI, ARG_PIK, ARG_SEGSITES, ARG_EHH, ARG_EHHK, ARG_TAJ_D, ARG_FAY_WU_H};

#endif