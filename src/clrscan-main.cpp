/* clrscan -- a program to calculate window-based diversity statistics
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
#include "clrscan-wintools.h"
#include "clrscan-winstats.h"
#include "clrscan-data.h"
#include "clrscan-cli.h"

using namespace std;

int main(int argc, char *argv[])
{
  cerr << "clrscan v" + VERSION + "\n";
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
  params.addFlag(ARG_LASSI, DEFAULT_LASSI, "", HELP_LASSI);
  params.addFlag(ARG_LASSI_CHOICE, DEFAULT_LASSI_CHOICE, "", HELP_LASSI_CHOICE);
  params.addFlag(ARG_HAPSTATS, DEFAULT_HAPSTATS, "", HELP_HAPSTATS);
  
  // Other flags
  //params.addFlag(ARG_INIT, DEFAULT_INIT, "", HELP_INIT);
  params.addFlag(ARG_K, DEFAULT_K, "", HELP_K);
  //params.addFlag(ARG_FINALIZE, DEFAULT_FINALIZE, "", HELP_FINALIZE);  
  params.addFlag(ARG_UNPHASED, DEFAULT_UNPHASED, "", HELP_UNPHASED);

  
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
  bool POP = (vcfFilename.compare(DEFAULT_FILENAME_POPFILE) == 0) ? false : true;
  vector<string> spectraFiles = params.getStringListFlag(ARG_FILENAME_SPECFILES);
  
  // Window control
  int WINSIZE = params.getIntFlag(ARG_WINSIZE);
  int WINSTEP = params.getIntFlag(ARG_WINSTEP);

  // Statistics
  bool LASSI = params.getBoolFlag(ARG_LASSI);
  int LASSI_CHOICE = params.getIntFlag(ARG_LASSI_CHOICE);
  bool HAPSTATS = params.getBoolFlag(ARG_HAPSTATS);

  // Other flags
  
  int K = params.getIntFlag(ARG_K);
  
  bool PHASED = !(params.getBoolFlag(ARG_UNPHASED));

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
  }

  if (!LASSI && !HAPSTATS){
    cerr << "ERROR: Must use --lassi, --hapstats, or both.\n";
    ERROR = true;
  }

  if (LASSI){
    if(K < 1){
      cerr << "ERROR: K must be >= 1.\n";
      ERROR = true;
    }
  }
  
/*
  if (TPED && VCF) {
    cerr << "ERROR: Must provide a TPED or VCF not both.\n";
    ERROR = true;
  }
*/
  if (ERROR) {
    return 1;
  }

  string ending;
  

 
  if(INIT){
    PopData *popData = readPopData(popFilename);

    if(PHASED) checkK(popData,double(K)/2.0);
    else if(!PHASED) checkK(popData,double(K));

    if(PHASED) ending = ".hap.";
    if(!PHASED) ending = ".mlg.";

    string outfile;
    if(LASSI) outfile = outfileBase + ending + "clrscan.spectra.gz";
    if(!LASSI && HAPSTATS) outfile = outfileBase + ending + "clrscan.hapstats.gz";

    map< string, HaplotypeData* > *hapDataByPop = readHaplotypeDataVCF(vcfFilename, popData, PHASED); 

    LASSIInitialResults *results = initResults(hapDataByPop, popData, WINSIZE, WINSTEP, K, HAPSTATS);
    LASSI_work_order_t *order;
    pthread_t *peer = new pthread_t[numThreads];
    int prev_index = 0;
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
    writeLASSIInitialResults(outfile, results, hapDataByPop->at(popData->popOrder[0])->map, popData, K, LASSI, HAPSTATS, PHASED);
  }
  else{//finalize
    

    map<string, vector<SpectrumData *>* > *specDataByPopByChr = readSpecData(spectraFiles);
    map<string, SpectrumData* > *avgSpecByPop = averageSpec(specDataByPopByChr);
    map<string, vector<LASSIResults *>* > *resultsByPopByChr = initResults(specDataByPopByChr);

    LASSI_work_order2_t *order;
    pthread_t *peer = new pthread_t[numThreads];
    int prev_index = 0;
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
    if(specDataByPopByChr->begin()->second->at(0)->PHASED){
      ending = ".hap.";
    }
    else{
      ending = ".mlg.";
    }
    string outfile = outfileBase + ending + "clrscan.out.gz";
    writeLASSIFinalResults(outfile, resultsByPopByChr, specDataByPopByChr);
  
  }

  return 0;
}

