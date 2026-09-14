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
#include <iostream>
#include <fstream>
#include <string>
#include "param_t.h"
#include "lassip-wintools.h"
#include "lassip-winstats.h"
#include "lassip-data.h"
#include "lassip-cli.h"

using namespace std;

namespace {

//Everything the command line said, unpacked once. INIT selects the stage:
//--vcf computes spectra from genotypes, --spectra computes statistics from
//spectra, and most flags are meaningful in only one of them.
struct Config
{
    int numThreads;
    string mapFilename;
    bool MAP;
    string vcfFilename;
    bool VCF;
    string outfileBase;
    string popFilename;
    bool POP;
    vector<string> spectraFiles;
    int WINSIZE;
    int WINSTEP;
    bool LASSI;
    int LASSI_CHOICE;
    bool HAPSTATS;
    bool SALTI;
    bool CALC_SPEC;
    bool AVG_SPEC;
    string nullSpecFile;
    int K;
    bool PHASED;
    int FILTER_LEVEL;
    bool KEEP_MONO;
    string DIST_TYPE;
    double FILTER_LMISS;
    double FILTER_HMISS;
    int MATCH_TOL;
    double MAX_GAP;
    double MAX_EXTEND_BP;
    double MAX_EXTEND_NW;
    double MAX_EXTEND_CM;
    bool INIT;
    param_t *params;
};

//Flags are registered in the order they should appear under their help
//section; the section name is the label passed to addFlag.
void registerFlags(param_t &params)
{
  params.addFlag(ARG_THREADS, DEFAULT_THREADS, "General", HELP_THREADS);

  // I/O flags
  params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "Input and output", HELP_OUTFILE);
  params.addFlag(ARG_FILENAME_MAP, DEFAULT_FILENAME_MAP, "Input and output", HELP_FILENAME_MAP);
  params.addFlag(ARG_FILENAME_POP1_VCF, DEFAULT_FILENAME_POP1_VCF, "Input and output", HELP_FILENAME_POP1_VCF);
  params.addFlag(ARG_FILENAME_POPFILE, DEFAULT_FILENAME_POPFILE, "Input and output", HELP_FILENAME_POPFILE);
  params.addListFlag(ARG_FILENAME_SPECFILES, DEFAULT_FILENAME_SPECFILES, "Input and output", HELP_FILENAME_SPECFILES);
  
  // Window control flags
  params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "Windows", HELP_WINSIZE);
  params.addFlag(ARG_WINSTEP, DEFAULT_WINSTEP, "Windows", HELP_WINSTEP);
  
  // Statistics flags
  params.addFlag(ARG_CALC_SPEC, DEFAULT_CALC_SPEC, "Statistics", HELP_CALC_SPEC);
  params.addFlag(ARG_AVG_SPEC, DEFAULT_AVG_SPEC, "Statistics", HELP_AVG_SPEC);
  params.addFlag(ARG_NULL_SPEC, DEFAULT_NULL_SPEC, "Statistics", HELP_NULL_SPEC);
  params.addFlag(ARG_LASSI, DEFAULT_LASSI, "Statistics", HELP_LASSI);
  params.addFlag(ARG_LASSI_CHOICE, DEFAULT_LASSI_CHOICE, "Statistics", HELP_LASSI_CHOICE);
  params.addFlag(ARG_HAPSTATS, DEFAULT_HAPSTATS, "Statistics", HELP_HAPSTATS);
  params.addFlag(ARG_SALTI, DEFAULT_SALTI, "Statistics", HELP_SALTI);
  
  // Other flags
  params.addFlag(ARG_K, DEFAULT_K, "Statistics", HELP_K);
  params.addFlag(ARG_UNPHASED, DEFAULT_UNPHASED, "Input and output", HELP_UNPHASED);
  params.addFlag(ARG_FILTER_LEVEL, DEFAULT_FILTER_LEVEL, "Filtering", HELP_FILTER_LEVEL);
  params.addFlag(ARG_FILTER_LMISS, DEFAULT_FILTER_LMISS, "Filtering", HELP_FILTER_LMISS);
  params.addFlag(ARG_FILTER_HMISS, DEFAULT_FILTER_HMISS, "Filtering", HELP_FILTER_HMISS);
  params.addFlag(ARG_MATCH_TOL, DEFAULT_MATCH_TOL, "Filtering", HELP_MATCH_TOL);
  params.addFlag(ARG_SEED, DEFAULT_SEED, "General", HELP_SEED);
  params.addFlag(ARG_DIST_TYPE, DEFAULT_DIST_TYPE, "saltiLASSI", HELP_DIST_TYPE);
  params.addFlag(ARG_MAX_GAP, DEFAULT_MAX_GAP, "saltiLASSI", HELP_MAX_GAP);
  params.addFlag(ARG_MAX_EXTEND_BP, DEFAULT_MAX_EXTEND_BP, "saltiLASSI", HELP_MAX_EXTEND_BP);
  params.addFlag(ARG_MAX_EXTEND_CM, DEFAULT_MAX_EXTEND_CM, "saltiLASSI", HELP_MAX_EXTEND_CM);
  params.addFlag(ARG_MAX_EXTEND_NW, DEFAULT_MAX_EXTEND_NW, "saltiLASSI", HELP_MAX_EXTEND_NW);
  params.addFlag(ARG_KEEP_MONO, DEFAULT_KEEP_MONO, "Filtering", HELP_KEEP_MONO);
}

Config readConfig(param_t &params)
{
  int numThreads = params.getIntFlag(ARG_THREADS);

  // I/O
  string mapFilename = params.getStringFlag(ARG_FILENAME_MAP);
  bool MAP = params.isFlagSet(ARG_FILENAME_MAP);
  string vcfFilename = params.getStringFlag(ARG_FILENAME_POP1_VCF);
  bool VCF = params.isFlagSet(ARG_FILENAME_POP1_VCF);
  string outfileBase = params.getStringFlag(ARG_OUTFILE);
  string popFilename = params.getStringFlag(ARG_FILENAME_POPFILE);
  bool POP = params.isFlagSet(ARG_FILENAME_POPFILE);
  vector<string> spectraFiles = params.getStringListFlag(ARG_FILENAME_SPECFILES);
  
  // Window control
  int WINSIZE = params.getIntFlag(ARG_WINSIZE);
  int WINSTEP = params.getIntFlag(ARG_WINSTEP);

  // Statistics
  bool LASSI = params.getBoolFlag(ARG_LASSI);
  int LASSI_CHOICE = params.getIntFlag(ARG_LASSI_CHOICE);
  bool HAPSTATS = params.getBoolFlag(ARG_HAPSTATS);
  bool SALTI = params.getBoolFlag(ARG_SALTI);
  bool CALC_SPEC = params.getBoolFlag(ARG_CALC_SPEC);
  bool AVG_SPEC = params.getBoolFlag(ARG_AVG_SPEC);
  string nullSpecFile = params.getStringFlag(ARG_NULL_SPEC);

  // Other flags
  int K = params.getIntFlag(ARG_K);
  bool PHASED = !(params.getBoolFlag(ARG_UNPHASED));
  int FILTER_LEVEL = params.getIntFlag(ARG_FILTER_LEVEL);
  bool KEEP_MONO = params.getBoolFlag(ARG_KEEP_MONO);
  string DIST_TYPE = params.getStringFlag(ARG_DIST_TYPE);
  double FILTER_LMISS = params.getDoubleFlag(ARG_FILTER_LMISS);
  double FILTER_HMISS = params.getDoubleFlag(ARG_FILTER_HMISS);
  int MATCH_TOL = params.getIntFlag(ARG_MATCH_TOL);
  //string DIST_TYPE = "bp";
  double MAX_GAP = params.getDoubleFlag(ARG_MAX_GAP);
  double MAX_EXTEND_BP = params.getDoubleFlag(ARG_MAX_EXTEND_BP);
  double MAX_EXTEND_NW = params.getDoubleFlag(ARG_MAX_EXTEND_NW);
  double MAX_EXTEND_CM = params.getDoubleFlag(ARG_MAX_EXTEND_CM);

  Config cfg;
  cfg.numThreads = numThreads;
  cfg.mapFilename = mapFilename;
  cfg.MAP = MAP;
  cfg.vcfFilename = vcfFilename;
  cfg.VCF = VCF;
  cfg.outfileBase = outfileBase;
  cfg.popFilename = popFilename;
  cfg.POP = POP;
  cfg.spectraFiles = spectraFiles;
  cfg.WINSIZE = WINSIZE;
  cfg.WINSTEP = WINSTEP;
  cfg.LASSI = LASSI;
  cfg.LASSI_CHOICE = LASSI_CHOICE;
  cfg.HAPSTATS = HAPSTATS;
  cfg.SALTI = SALTI;
  cfg.CALC_SPEC = CALC_SPEC;
  cfg.AVG_SPEC = AVG_SPEC;
  cfg.nullSpecFile = nullSpecFile;
  cfg.K = K;
  cfg.PHASED = PHASED;
  cfg.FILTER_LEVEL = FILTER_LEVEL;
  cfg.KEEP_MONO = KEEP_MONO;
  cfg.DIST_TYPE = DIST_TYPE;
  cfg.FILTER_LMISS = FILTER_LMISS;
  cfg.FILTER_HMISS = FILTER_HMISS;
  cfg.MATCH_TOL = MATCH_TOL;
  cfg.MAX_GAP = MAX_GAP;
  cfg.MAX_EXTEND_BP = MAX_EXTEND_BP;
  cfg.MAX_EXTEND_NW = MAX_EXTEND_NW;
  cfg.MAX_EXTEND_CM = MAX_EXTEND_CM;
  cfg.INIT = VCF;
  cfg.params = &params;
  return cfg;
}

//Flags that belong to the other stage are warned about rather than ignored.
void warnCrossStageFlags(const Config &cfg)
{
  param_t &params = *cfg.params;
  const bool INIT = cfg.INIT;
  //Flags are shared by both stages but most only act in one of them. Silence
  //here is how example/do_lassip_YRI.bash came to pass --lassi to stage 1,
  //where it does nothing, and stopped reproducing the file committed beside it.
  if(INIT){
    const string stage2Only[] = {ARG_LASSI, ARG_SALTI, ARG_AVG_SPEC, ARG_NULL_SPEC,
                                 ARG_LASSI_CHOICE, ARG_MAX_EXTEND_BP, ARG_MAX_EXTEND_CM,
                                 ARG_MAX_EXTEND_NW, ARG_MAX_GAP};  //--map already draws a hard error here
    for (unsigned int i = 0; i < sizeof(stage2Only)/sizeof(stage2Only[0]); i++){
      if(params.isFlagSet(stage2Only[i])){
        cerr << "WARNING: " << stage2Only[i] << " has no effect with --vcf; it applies to the --spectra stage.\n";
      }
    }
  }
  else{
    const string stage1Only[] = {ARG_CALC_SPEC, ARG_HAPSTATS, ARG_WINSIZE, ARG_WINSTEP,
                                 ARG_K, ARG_UNPHASED, ARG_FILENAME_POPFILE, ARG_FILTER_LEVEL,
                                 ARG_FILTER_LMISS, ARG_FILTER_HMISS, ARG_KEEP_MONO,
                                 ARG_MATCH_TOL, ARG_SEED};
    for (unsigned int i = 0; i < sizeof(stage1Only)/sizeof(stage1Only[0]); i++){
      if(params.isFlagSet(stage1Only[i])){
        cerr << "WARNING: " << stage1Only[i] << " has no effect with --spectra; it applies to the --vcf stage.\n";
      }
    }
  }
}

//True if the command line is usable. Every problem is reported, not just the
//first, so one run tells the user everything that is wrong.
bool validate(const Config &cfg)
{
  const int numThreads = cfg.numThreads;
  const bool MAP = cfg.MAP;
  const bool POP = cfg.POP;
  const int WINSIZE = cfg.WINSIZE;
  const int WINSTEP = cfg.WINSTEP;
  const bool LASSI = cfg.LASSI;
  const int LASSI_CHOICE = cfg.LASSI_CHOICE;
  const bool HAPSTATS = cfg.HAPSTATS;
  const bool SALTI = cfg.SALTI;
  const bool CALC_SPEC = cfg.CALC_SPEC;
  const bool AVG_SPEC = cfg.AVG_SPEC;
  const int K = cfg.K;
  const int FILTER_LEVEL = cfg.FILTER_LEVEL;
  const string &DIST_TYPE = cfg.DIST_TYPE;
  const double FILTER_LMISS = cfg.FILTER_LMISS;
  const double FILTER_HMISS = cfg.FILTER_HMISS;
  const int MATCH_TOL = cfg.MATCH_TOL;
  const double MAX_GAP = cfg.MAX_GAP;
  const double MAX_EXTEND_BP = cfg.MAX_EXTEND_BP;
  const double MAX_EXTEND_NW = cfg.MAX_EXTEND_NW;
  const double MAX_EXTEND_CM = cfg.MAX_EXTEND_CM;
  param_t &params = *cfg.params;
  const bool INIT = cfg.INIT;
  const bool FINALIZE = !cfg.INIT;
  bool ERROR = false;

  //--dist-type selects the coordinate the statistics are reported and extended
  //along; it is meaningful at both stages, and used to be checked only at the
  //second, so `--calc-spec --dist-type cm` silently wrote a nameless column.
  if(DIST_TYPE.compare("bp") != 0 &&
     DIST_TYPE.compare("cm") != 0 &&
     DIST_TYPE.compare("nw") != 0){
    cerr << "ERROR: --dist-type must be one of bp, cm or nw.\n";
    ERROR = true;
  }

  if(!INIT && !FINALIZE){
    cerr << "ERROR: Must specify either --vcf or --spectra.\n";
    ERROR = true;
  }

  if(INIT && FINALIZE){
    cerr << "ERROR: Must specify either --vcf or --spectra.\n";
    ERROR = true;
  }

  if (numThreads <= 0) {
    cerr << "ERROR: Must specify a positive number of threads.\n";
    ERROR = true;
  }

  if(INIT){
    if (!POP){
      cerr << "ERROR: Must provide a map from ind to pop with --pop.\n";
      ERROR = true;
    }

    if (WINSIZE < 1) {
      cerr << "ERROR: Window size needs to be greater than 0.\n";
      ERROR = true;
    }

    if (WINSTEP < 1) {
      cerr << "ERROR: Window step size needs to be greater than 0.\n";
      ERROR = true;
    }

    if(FILTER_LEVEL != 0 && FILTER_LEVEL != 1 && FILTER_LEVEL != 2){
      cerr << "ERROR: Filter level must be 0, 1, or 2.\n";
      ERROR = true;
    }

    if(FILTER_LMISS < 0 || FILTER_LMISS > 1){
      cerr << "ERROR: Missing data locus filter must be in [0,1].\n";
      ERROR = true;
    }

    if(FILTER_HMISS < 0 || FILTER_HMISS > 1){
      cerr << "ERROR: Missing data halotype filter must be in [0,1].\n";
      ERROR = true;
    }

    if(MATCH_TOL < 0){
      cerr << "ERROR: Haplotype match tolerance must be an integer >= 0.\n";
      ERROR = true;
    }

    if (!CALC_SPEC && !HAPSTATS){
      cerr << "ERROR: Must use --calc-spec or --hapstats.\n";
      ERROR = true;
    }    
    if (CALC_SPEC && K < 1){
      cerr << "ERROR: K must be >= 1.\n";
      ERROR = true;
    }

    if(MAP){
      cerr << "ERROR: Map file not required at this stage.\n";
      ERROR = true;
    }
  }

  warnCrossStageFlags(cfg);

  if(FINALIZE){
    if(!params.isFlagSet(ARG_FILENAME_SPECFILES)){
      cerr << "ERROR: Must provide spectra files to calculate LASSI statistic.\n";
      ERROR = true;
    }
    if(LASSI_CHOICE < 1 || LASSI_CHOICE > 5){
      cerr << "ERROR: --lassi-choice must be an integer in {1..5}.\n";
      ERROR = true;
    }
  
    if (!LASSI && !SALTI && !AVG_SPEC){
      cerr << "ERROR: Must use --lassi or --salti for analyzing haplotype spectra.\n";
      cerr << "\tOr use --avg-spec to compute average spectra from *.spectra files.\n";
      ERROR = true;
    }
    if(LASSI && SALTI){
      cerr << "ERROR: Must choose only one of --lassi or --salti for analyzing haplotype spectra.\n";
      ERROR = true;
    }

    if(DIST_TYPE.compare("bp") != 0 &&
      DIST_TYPE.compare("cm") != 0 &&
      DIST_TYPE.compare("nw") != 0){
      cerr << "ERROR: Must choose bp, cm, or nw for distance measure.\n";
      ERROR = true;
    }

    if(LASSI && MAP){
      cerr << "ERROR: Map file only used for saltiLASSI computations.\n";
      ERROR = true;
    }

    if(MAP && DIST_TYPE.compare("cm") != 0){
      cerr << "ERROR: Must choose --dist-type cm when providing a map file.\n";
      ERROR = true;
    }

    if(MAX_GAP < 0){
      cerr << "ERROR: --max-gap must be >= 0 (0 interpolates across any gap).\n";
      ERROR = true;
    }

    if(params.isFlagSet(ARG_MAX_GAP) && DIST_TYPE.compare("cm") != 0){
      cerr << "WARNING: --max-gap has no effect without --dist-type cm.\n";
    }

    if(!MAP && DIST_TYPE.compare("cm") == 0){
      cerr << "ERROR: Must provide a map file when choosing --dist-type cm.\n";
      ERROR = true;
    }

    if(SALTI && DIST_TYPE.compare("bp") == 0 && MAX_EXTEND_BP < 1){
      cerr << "ERROR: MAX_EXTEND (bp) must be >= 1.\n";
      ERROR = true;
    }

    if(SALTI && DIST_TYPE.compare("cm") == 0 && MAX_EXTEND_CM <= 0){
      cerr << "ERROR: MAX_EXTEND (cm) must be > 0.\n";
      ERROR = true;
    }

    if(SALTI && DIST_TYPE.compare("nw") == 0 && MAX_EXTEND_NW <= 0){
      cerr << "ERROR: MAX_EXTEND (nw) must be > 0.\n";
      ERROR = true;
    }
  }

  return !ERROR;
}

//Stage 1: genotypes in, one haplotype frequency spectrum per window out.
int runSpectra(const Config &cfg)
{
  const int numThreads = cfg.numThreads;
  const string &vcfFilename = cfg.vcfFilename;
  const string &outfileBase = cfg.outfileBase;
  const string &popFilename = cfg.popFilename;
  const int WINSIZE = cfg.WINSIZE;
  const int WINSTEP = cfg.WINSTEP;
  const bool HAPSTATS = cfg.HAPSTATS;
  const bool CALC_SPEC = cfg.CALC_SPEC;
  const int K = cfg.K;
  const bool PHASED = cfg.PHASED;
  const int FILTER_LEVEL = cfg.FILTER_LEVEL;
  const bool KEEP_MONO = cfg.KEEP_MONO;
  const string &DIST_TYPE = cfg.DIST_TYPE;
  const double FILTER_LMISS = cfg.FILTER_LMISS;
  param_t &params = *cfg.params;

    PopData *popData = readPopData(popFilename);

    if(PHASED) checkK(popData,double(K)/2.0);
    else if(!PHASED) checkK(popData,double(K));

    map< string, HaplotypeData* > *hapDataByPop = readHaplotypeDataVCF(vcfFilename, popData, PHASED, (FILTER_LEVEL < 2));

    if(FILTER_LEVEL > 0){
      hapDataByPop = filterHaplotypeData(hapDataByPop, popData, FILTER_LEVEL, FILTER_LMISS, KEEP_MONO, PHASED);
    } 

    LASSIInitialResults *results = initResults(hapDataByPop, popData, WINSIZE, WINSTEP, K, HAPSTATS, DIST_TYPE);
    //One work cursor per population; threads claim chunks of windows from it.
    WorkCursor cursor;
    cursor.nunits = popData->npops;
    cursor.next = new std::atomic<unsigned int>[cursor.nunits];
    for (unsigned int u = 0; u < cursor.nunits; u++) cursor.next[u] = 0;

    vector<LASSI_work_order_t> orders(numThreads);
    vector<std::thread> peer;
    for (int i = 0; i < numThreads; i++){
      orders[i].id = i;
      orders[i].cursor = &cursor;
      orders[i].nullWins.assign(popData->npops, 0);
      orders[i].hapDataByPop = hapDataByPop;
      orders[i].popData = popData;
      orders[i].params = &params;
      orders[i].results = results;
      peer.push_back(std::thread(calc_LASSI_stats, &orders[i]));
    }
    for (int i = 0; i < numThreads; i++) peer[i].join();
    for (int i = 0; i < numThreads; i++)
      for (int pop = 0; pop < popData->npops; pop++)
        results->pops[pop].nullWins += orders[i].nullWins[pop];
    delete [] cursor.next;
    cerr << "Done.\n";
    writeLASSIInitialResults(outfileBase, results, hapDataByPop, popData, K, CALC_SPEC, HAPSTATS, PHASED, FILTER_LEVEL, DIST_TYPE);

  return 0;
}

//Stage 2: spectra in, LASSI or saltiLASSI statistics out.
int runStatistics(const Config &cfg)
{
  const int numThreads = cfg.numThreads;
  const string &mapFilename = cfg.mapFilename;
  const double MAX_GAP = cfg.MAX_GAP;
  const string &outfileBase = cfg.outfileBase;
  const vector<string> &spectraFiles = cfg.spectraFiles;
  const bool LASSI = cfg.LASSI;
  const int LASSI_CHOICE = cfg.LASSI_CHOICE;
  const bool SALTI = cfg.SALTI;
  const bool AVG_SPEC = cfg.AVG_SPEC;
  const string &nullSpecFile = cfg.nullSpecFile;
  int K = cfg.K;
  const string &DIST_TYPE = cfg.DIST_TYPE;
  param_t &params = *cfg.params;
  bool ERROR = false;
  string ending;

    map<string, vector<SpectrumData *>* > *specDataByPopByChr = readSpecData(spectraFiles);
    map<string, SpectrumData* > *avgSpecByPop;
    
    if(params.isFlagSet(ARG_NULL_SPEC)){
      avgSpecByPop = averageSpec(nullSpecFile);
      if(!checkNull(avgSpecByPop,specDataByPopByChr)) return EXIT_DATAERR;
    }
    else avgSpecByPop = averageSpec(specDataByPopByChr);

    if(AVG_SPEC){
      writeAverageSpec(outfileBase,avgSpecByPop);
      return 0;
    }

    //Check null spectra to make sure p_K >= 1/(100*K) otherwise grid search won't work.
    //This condition can happen when there is not enough data to estimate the truncated
    //HFS well, and the Kth highest frequency haplotype is < 1/(100*K)
    for(map<string, SpectrumData* >::iterator it = avgSpecByPop->begin(); it != avgSpecByPop->end(); it++){
      SpectrumData *spec = it->second;
      //K comes from the spectra file at this stage; --k is not read here, and
      //using it made this threshold wrong unless the user happened to repeat
      //the same --k they used to build the spectra.
      if(spec->freq[0][spec->K-1] < 1.0/(100.0*double(spec->K))){
        cerr << "ERROR: Null spectrum " << it->first << " likely not well-estimated. Kth frequency class is < " << 1.0/(100.0*double(spec->K)) << ", preventing grid search.\n";
        ERROR = true;
      }
    }

    if(ERROR) return EXIT_DATAERR;

    map<string, vector<LASSIResults *>* > *resultsByPopByChr = initResults(specDataByPopByChr, SALTI);

    if(LASSI){
      //One cursor per population-contig pair, visited in the same order by the
      //workers as it is counted here.
      WorkCursor cursor;
      cursor.nunits = 0;
      for (map<string, vector<SpectrumData *>* >::iterator u = specDataByPopByChr->begin();
           u != specDataByPopByChr->end(); u++) cursor.nunits += u->second->size();
      cursor.next = new std::atomic<unsigned int>[cursor.nunits];
      for (unsigned int u = 0; u < cursor.nunits; u++) cursor.next[u] = 0;

      vector<LASSI_work_order2_t> orders(numThreads);
      vector<std::thread> peer;
      for (int i = 0; i < numThreads; i++){
        orders[i].id = i;
        orders[i].cursor = &cursor;
        orders[i].specDataByPopByChr = specDataByPopByChr;
        orders[i].avgSpecByPop = avgSpecByPop;
        orders[i].resultsByPopByChr = resultsByPopByChr;
        orders[i].params = &params;
        peer.push_back(std::thread(calc_LASSI_stats2, &orders[i]));
      }
      for (int i = 0; i < numThreads; i++) peer[i].join();
      delete [] cursor.next;
      cerr << "Done.\n";
    }
    else if (SALTI){
      cerr << "saltiLASSI\n";

      if(DIST_TYPE.compare("nw") == 0){
        //populate specDataByPopByChr->chr->pop->dist[] with integers
        fillNWDistance(specDataByPopByChr);
      }
      else if (DIST_TYPE.compare("cm") == 0){
        //load genetic map from file
        GMapData geneticMap(mapFilename,MAX_GAP);
        //populate specDataByPopByChr->chr->pop->dist[] with genetic distances 
        fillCMDistance(specDataByPopByChr,geneticMap);
      }

      map<string, vector<SpectrumData *>* >::iterator it;
      //pop
      for (it = specDataByPopByChr->begin(); it != specDataByPopByChr->end(); it++){
        string popName = it->first;
        cerr << popName << endl;
        SpectrumData *avgSpec = avgSpecByPop->at(popName);
        vector<SpectrumData *> *specDataByChr = specDataByPopByChr->at(popName);
        vector<LASSIResults *> *resultsByChr = resultsByPopByChr->at(popName);

        double dmin = getDMin(specDataByChr);

        //chr
        for(unsigned int c = 0; c < specDataByChr->size(); c++){
          K = specDataByChr->at(c)->K;
          double U = avgSpec->freq[0][K-1];

          //The sweep spectra depend on the null spectrum, the scaling choice, m
          //and epsilon, but not on the window, so one table serves every window.
          double **f = calcF(LASSI_CHOICE, K);
          double ***q = initQ(K, U);
          calcQ(q, avgSpec, f);
          for(int i = 0; i < K; i++) delete [] f[i];
          delete [] f;

          WorkCursor cursor;
          cursor.nunits = 1;
          cursor.next = new std::atomic<unsigned int>[1];
          cursor.next[0] = 0;

          vector<SALTI_work_order_t> orders(numThreads);
          vector<std::thread> peer;
          for (int i = 0; i < numThreads; i++){
            orders[i].id = i;
            orders[i].cursor = &cursor;
            orders[i].specData = specDataByChr->at(c);
            orders[i].avgSpec = avgSpec;
            orders[i].results = resultsByChr->at(c);
            orders[i].params = &params;
            orders[i].q = q;
            orders[i].dmin = dmin;
            peer.push_back(std::thread(calc_SALTI_stats, &orders[i]));
          }
          for (int i = 0; i < numThreads; i++) peer[i].join();
          delete [] cursor.next;

          releaseQ(q,K,U);

          cerr << "Done with contig " << specDataByChr->at(c)->info[0][0] << ".\n";
        }
      }
    }


    if(specDataByPopByChr->begin()->second->at(0)->PHASED){
      ending = ".lassip.hap.";
    }
    else{
      ending = ".lassip.mlg.";
    }
    string outfile = outfileBase + ending + "out.gz";
    writeLASSIFinalResults(outfile, resultsByPopByChr, specDataByPopByChr, SALTI);

  return 0;
}

} // namespace

//The body of the program. main() at the bottom of this file is a thin wrapper
//that turns the exceptions thrown from here and from the data layer into
//distinct exit codes.
int lassipMain(int argc, char *argv[])
{
  param_t params;
  params.setPreamble(PREAMBLE);
  params.setUsage(USAGE);
  params.setVersion("lassip v" + VERSION);
  registerFlags(params);

  if (argc == 1){
    cerr << USAGE << "\n";
    cerr << "Run lassip --help for the full list of options.\n";
    return EXIT_USAGE;
  }

  params.parseCommandLine(argc, argv);

  cerr << "lassip v" + VERSION + "\n";

  Config cfg = readConfig(params);
  if (!validate(cfg)) return EXIT_USAGE;

  return cfg.INIT ? runSpectra(cfg) : runStatistics(cfg);
}

int main(int argc, char *argv[])
{
  //Distinct exit codes so that a caller can tell a bad command line from bad
  //input data. The data layer signals failure by throwing an int; --help and
  //--version unwind through ParamExit with code 0.
  try {
    return lassipMain(argc, argv);
  }
  catch (const ParamExit &e){
    return e.code;
  }
  catch (int){
    cerr << "lassip: exiting after the error above.\n";
    return EXIT_DATAERR;
  }
  catch (const exception &e){
    cerr << "lassip: " << e.what() << "\n";
    return EXIT_INTERNAL;
  }
  catch (...){
    cerr << "lassip: unknown error.\n";
    return EXIT_INTERNAL;
  }
}
