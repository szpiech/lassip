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

int main(int argc, char *argv[])
{
  cerr << "lassip v" + VERSION + "\n";
  param_t params;
  params.setPreamble(PREAMBLE);

  params.addFlag(ARG_THREADS, DEFAULT_THREADS, "", HELP_THREADS);

  // I/O flags
  params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
  params.addFlag(ARG_FILENAME_POP1_VCF, DEFAULT_FILENAME_POP1_VCF, "", HELP_FILENAME_POP1_VCF);
  params.addFlag(ARG_FILENAME_POPFILE, DEFAULT_FILENAME_POPFILE, "", HELP_FILENAME_POPFILE);
  params.addListFlag(ARG_FILENAME_SPECFILES, DEFAULT_FILENAME_SPECFILES, "", HELP_FILENAME_SPECFILES);
  
  // Window control flags
  params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
  params.addFlag(ARG_WINSTEP, DEFAULT_WINSTEP, "", HELP_WINSTEP);
  
  // Statistics flags
  params.addFlag(ARG_CALC_SPEC, DEFAULT_CALC_SPEC, "", HELP_CALC_SPEC);
  params.addFlag(ARG_AVG_SPEC, DEFAULT_AVG_SPEC, "", HELP_AVG_SPEC);
  params.addFlag(ARG_NULL_SPEC, DEFAULT_NULL_SPEC, "", HELP_NULL_SPEC);
  params.addFlag(ARG_LASSI, DEFAULT_LASSI, "", HELP_LASSI);
  params.addFlag(ARG_LASSI_CHOICE, DEFAULT_LASSI_CHOICE, "", HELP_LASSI_CHOICE);
  params.addFlag(ARG_HAPSTATS, DEFAULT_HAPSTATS, "", HELP_HAPSTATS);
  params.addFlag(ARG_SALTI, DEFAULT_SALTI, "", HELP_SALTI);
  
  // Other flags
  params.addFlag(ARG_K, DEFAULT_K, "", HELP_K);
  params.addFlag(ARG_UNPHASED, DEFAULT_UNPHASED, "", HELP_UNPHASED);
  params.addFlag(ARG_FILTER_LEVEL, DEFAULT_FILTER_LEVEL, "", HELP_FILTER_LEVEL);
  //params.addFlag(ARG_DIST_TYPE, DEFAULT_DIST_TYPE, "", HELP_DIST_TYPE);
  params.addFlag(ARG_MAX_EXTEND_BP, DEFAULT_MAX_EXTEND_BP, "", HELP_MAX_EXTEND_BP);
  //params.addFlag(ARG_MAX_EXTEND_CM, DEFAULT_MAX_EXTEND_CM, "", HELP_MAX_EXTEND_CM);
  //params.addFlag(ARG_MAX_EXTEND_NW, DEFAULT_MAX_EXTEND_NW, "", HELP_MAX_EXTEND_NW);

  try {
    params.parseCommandLine(argc, argv);
  }
  catch (...) {
    return 1;
  }

  int numThreads = params.getIntFlag(ARG_THREADS);

  // I/O
  string vcfFilename = params.getStringFlag(ARG_FILENAME_POP1_VCF);
  bool VCF = (vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) == 0) ? false : true;
  string outfileBase = params.getStringFlag(ARG_OUTFILE);
  string popFilename = params.getStringFlag(ARG_FILENAME_POPFILE);
  bool POP = (popFilename.compare(DEFAULT_FILENAME_POPFILE) == 0) ? false : true;
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
  //string DIST_TYPE = params.getStringFlag(ARG_DIST_TYPE);
  string DIST_TYPE = "bp";
  double MAX_EXTEND_BP = params.getDoubleFlag(ARG_MAX_EXTEND_BP);
  //double MAX_EXTEND_NW = params.getDoubleFlag(ARG_MAX_EXTEND_NW);

  // Check for consistency errors within flags
  bool ERROR = false;

  bool INIT = VCF;
  bool FINALIZE = !VCF;

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

    if(DIST_TYPE.compare("bp") != 0 &&
      DIST_TYPE.compare("cm") != 0 &&
      DIST_TYPE.compare("ns") != 0 &&
      DIST_TYPE.compare("nw") != 0){
      //cerr << "ERROR: Must choose bp, cm, ns, or nw for distance measure.\n";
      cerr << "ERROR: Must choose bp for distance measure.\n";
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
  }

  if(FINALIZE){
    if(spectraFiles.size() == 1 && spectraFiles[0].compare(DEFAULT_FILENAME_SPECFILES) == 0){
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
    if(SALTI && MAX_EXTEND_BP < 1){
      cerr << "ERROR: MAX_EXTEND (bp) must be >= 1.\n";
      ERROR = true;
    }
  }

  if (ERROR) {
    return 1;
  }

  string ending;
  

 
  if(INIT){
    PopData *popData = readPopData(popFilename);

    if(PHASED) checkK(popData,double(K)/2.0);
    else if(!PHASED) checkK(popData,double(K));

    map< string, HaplotypeData* > *hapDataByPop = readHaplotypeDataVCF(vcfFilename, popData, PHASED, (FILTER_LEVEL < 2));

    if(FILTER_LEVEL > 0){
      hapDataByPop = filterHaplotypeData(hapDataByPop, popData, FILTER_LEVEL);
    } 

    LASSIInitialResults *results = initResults(hapDataByPop, popData, WINSIZE, WINSTEP, K, HAPSTATS, DIST_TYPE);
    LASSI_work_order_t *order;
    pthread_t *peer = new pthread_t[numThreads];
    
    for (int i = 0; i < numThreads; i++){
      order = new LASSI_work_order_t;
      order->id = i;
      order->hapDataByPop = hapDataByPop;
      order->popData = popData;
      order->params = &params;
      order->results = results;
      pthread_create(&(peer[i]),
                     NULL,
                     (void *(*)(void *))calc_LASSI_stats,
                     (void *)order);
    }
    for (int i = 0; i < numThreads; i++) pthread_join(peer[i], NULL);
    delete [] peer;
    cerr << "Done.\n";
    writeLASSIInitialResults(outfileBase, results, hapDataByPop, popData, K, CALC_SPEC, HAPSTATS, PHASED, FILTER_LEVEL, DIST_TYPE);
  }
  else{//finalize
    

    map<string, vector<SpectrumData *>* > *specDataByPopByChr = readSpecData(spectraFiles);
    map<string, SpectrumData* > *avgSpecByPop;
    
    if(nullSpecFile.compare(DEFAULT_NULL_SPEC) != 0){
      avgSpecByPop = averageSpec(nullSpecFile);
      if(!checkNull(avgSpecByPop,specDataByPopByChr)) return 1;
    }
    else avgSpecByPop = averageSpec(specDataByPopByChr);

    if(AVG_SPEC){
      writeAverageSpec(outfileBase,avgSpecByPop);
      return 0;
    }

    map<string, vector<LASSIResults *>* > *resultsByPopByChr = initResults(specDataByPopByChr, SALTI);

    if(LASSI){
      LASSI_work_order2_t *order;
      pthread_t *peer = new pthread_t[numThreads];
    
      for (int i = 0; i < numThreads; i++){
        order = new LASSI_work_order2_t;
        order->id = i;
        order->specDataByPopByChr = specDataByPopByChr;
        order->avgSpecByPop = avgSpecByPop;
        order->resultsByPopByChr = resultsByPopByChr;
        order->params = &params;

        pthread_create(&(peer[i]),
                        NULL,
                       (void *(*)(void *))calc_LASSI_stats2,
                       (void *)order);      
      }

      for (int i = 0; i < numThreads; i++) pthread_join(peer[i], NULL);
      delete [] peer;
      cerr << "Done.\n";
    }
    else if (SALTI){
      cerr << "saltiLASSI\n";
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
          SALTI_work_order_t *order;
          pthread_t *peer = new pthread_t[numThreads];
          int nwins = specDataByChr->at(c)->nwins;
          K = specDataByChr->at(c)->K;
          double U = avgSpec->freq[0][K-1];
          double ****q = initQ(nwins,K,U);
          //PRECOMPUTE all possible Qs
          for (int i = 0; i < numThreads; i++){
            order = new SALTI_work_order_t;
            order->id = i;
            order->specData = specDataByChr->at(c);
            order->avgSpec = avgSpec;
            order->results = resultsByChr->at(c);
            order->params = &params;
            order->q = q;

            pthread_create(&(peer[i]),
                            NULL,
                           (void *(*)(void *))calc_SALTI_stats1,
                           (void *)order);      
          }

          for (int i = 0; i < numThreads; i++) pthread_join(peer[i], NULL);
          delete [] peer;

          peer = new pthread_t[numThreads];
          for (int i = 0; i < numThreads; i++){
            order = new SALTI_work_order_t;
            order->id = i;
            order->specData = specDataByChr->at(c);
            order->avgSpec = avgSpec;
            order->results = resultsByChr->at(c);
            order->params = &params;
            order->q = q;
            order->dmin = dmin;

            pthread_create(&(peer[i]),
                            NULL,
                           (void *(*)(void *))calc_SALTI_stats2,
                           (void *)order);      
          }

          for (int i = 0; i < numThreads; i++) pthread_join(peer[i], NULL);
          delete [] peer;

          releaseQ(q,nwins,K,U);

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
  
  }

  return 0;
}

