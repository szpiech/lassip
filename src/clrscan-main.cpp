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
  //params.addFlag(ARG_FILENAME_TPED, DEFAULT_FILENAME_TPED, "", HELP_FILENAME_TPED);
  params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
  params.addFlag(ARG_FILENAME_POP1_VCF, DEFAULT_FILENAME_POP1_VCF, "", HELP_FILENAME_POP1_VCF);
  params.addFlag(ARG_FILENAME_POPFILE, DEFAULT_FILENAME_POPFILE, "", HELP_FILENAME_POPFILE);
  //params.addFlag(ARG_FILENAME_MAP, DEFAULT_FILENAME_MAP, "", HELP_FILENAME_MAP);
  params.addListFlag(ARG_FILENAME_SPECFILES, DEFAULT_FILENAME_SPECFILES, "", HELP_FILENAME_SPECFILES);
  
  // Window control flags
  //params.addFlag(ARG_BP, DEFAULT_BP, "", HELP_BP);
  //params.addFlag(ARG_SITES, DEFAULT_SITES, "", HELP_SITES);
  params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
  params.addFlag(ARG_WINSTEP, DEFAULT_WINSTEP, "", HELP_WINSTEP);
  //params.addListFlag(ARG_PARTITION, DEFAULT_PARTITION, "", HELP_PARTITION);

  // Statistics flags
  params.addFlag(ARG_LASSI, DEFAULT_LASSI, "", HELP_LASSI);
  params.addFlag(ARG_LASSI_CHOICE, DEFAULT_LASSI_CHOICE, "", HELP_LASSI_CHOICE);
  /*
  params.addFlag(ARG_PI, DEFAULT_PI, "", HELP_PI);
  params.addListFlag(ARG_PIK, DEFAULT_PIK, "", HELP_PIK);
  params.addFlag(ARG_SEGSITES, DEFAULT_SEGSITES, "", HELP_SEGSITES);
  params.addListFlag(ARG_EHH, DEFAULT_EHH, "", HELP_EHH);
  params.addListFlag(ARG_EHHK, DEFAULT_EHHK, "", HELP_EHHK);
  params.addFlag(ARG_TAJ_D, DEFAULT_TAJ_D, "", HELP_TAJ_D);
  params.addFlag(ARG_FAY_WU_H, DEFAULT_FAY_WU_H, "", HELP_FAY_WU_H);
  */
  // Other flags
  params.addFlag(ARG_INIT, DEFAULT_INIT, "", HELP_INIT);
  params.addFlag(ARG_K, DEFAULT_K, "", HELP_K);
  params.addFlag(ARG_FINALIZE, DEFAULT_FINALIZE, "", HELP_FINALIZE);
  //params.addFlag(ARG_EHH_PART, DEFAULT_EHH_PART, "", HELP_EHH_PART);
  //params.addFlag(ARG_NO_SFS_SUB, DEFAULT_NO_SFS_SUB, "", HELP_NO_SFS_SUB);
  //params.addFlag(ARG_PMAP, DEFAULT_PMAP, "", HELP_PMAP);
  
  try {
    params.parseCommandLine(argc, argv);
  }
  catch (...) {
    return 1;
  }

  int numThreads = params.getIntFlag(ARG_THREADS);

  // I/O
  //string tpedFilename = params.getStringFlag(ARG_FILENAME_TPED);
  //bool TPED = (tpedFilename.compare(DEFAULT_FILENAME_TPED) == 0) ? false : true;
  string vcfFilename = params.getStringFlag(ARG_FILENAME_POP1_VCF);
  bool VCF = (vcfFilename.compare(DEFAULT_FILENAME_POP1_VCF) == 0) ? false : true;
  string outfileBase = params.getStringFlag(ARG_OUTFILE);
  string popFilename = params.getStringFlag(ARG_FILENAME_POPFILE);
  bool POP = (vcfFilename.compare(DEFAULT_FILENAME_POPFILE) == 0) ? false : true;
  vector<string> spectraFiles = params.getStringListFlag(ARG_FILENAME_SPECFILES);
  //string mapFilename = params.getStringFlag(ARG_FILENAME_MAP);
  //bool MAP = (mapFilename.compare(DEFAULT_FILENAME_MAP) == 0) ? false : true;
  
  // Window control
  //bool USE_BP = params.getBoolFlag(ARG_BP);
  //bool USE_SITES = params.getBoolFlag(ARG_SITES);
  int WINSIZE = params.getIntFlag(ARG_WINSIZE);
  int WINSTEP = params.getIntFlag(ARG_WINSTEP);
  //vector<int> PARTITIONS = params.getIntListFlag(ARG_PARTITION);
  //bool DO_PARTITION = false;

  // Statistics
  bool LASSI = params.getBoolFlag(ARG_LASSI);
  int LASSI_CHOICE = params.getIntFlag(ARG_LASSI_CHOICE);
  //bool CALC_PI = params.getBoolFlag(ARG_PI);
  //vector<int> PIK_CHOICE = params.getIntListFlag(ARG_PIK);
  //bool CALC_PIK = false;
  //bool CALC_S = params.getBoolFlag(ARG_SEGSITES);
  //vector<int> EHH_WINS = params.getIntListFlag(ARG_EHH);
  //bool CALC_EHH = false;
  //vector<int> EHHK_CHOICES = params.getIntListFlag(ARG_EHHK);
  //bool CALC_EHHK = false;
  //bool CALC_TAJ_D = params.getBoolFlag(ARG_TAJ_D);
  //bool CALC_FAY_WU_H = params.getBoolFlag(ARG_FAY_WU_H);

  // Other flags
  bool INIT = params.getBoolFlag(ARG_INIT);
  int K = params.getIntFlag(ARG_K);
  bool FINALIZE = params.getBoolFlag(ARG_FINALIZE);
  //bool EHH_PART = params.getBoolFlag(ARG_EHH_PART);
  //bool PMAP = params.getBoolFlag(ARG_PMAP);
  //bool SFS_SUB = !(params.getBoolFlag(ARG_NO_SFS_SUB));

  // Check for consistency errors within flags
  bool ERROR = false;

  
/*
  if (!USE_SITES && !USE_BP){
    cerr << "ERROR: Must choose to measure windows in either sites or bps.\n";
    ERROR = true;
  }

  if (USE_SITES && USE_BP){
    cerr << "ERROR: Must choose to measure windows in either sites or bps not both.\n";
    ERROR = true;
  }
*/
  if(!INIT && !FINALIZE){
    cerr << "ERROR: Must specify either --initial or --finalize.\n";
    ERROR = true;
  }

  if(INIT && FINALIZE){
    cerr << "ERROR: Must specify either --initial or --finalize.\n";
    ERROR = true;
  }

  if (numThreads <= 0) {
    cerr << "ERROR: Must specify a positive number of threads.\n";
    ERROR = true;
  }

  if(INIT){
    if (!POP){
      cerr << "ERROR: Must provide a map from ind to pop with --pop for --initial runs.\n";
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
    if (/*!TPED &&*/ !VCF) {
      cerr << "ERROR: Must provide a file with genetic data.\n";
      ERROR = true;
    }
  }

  if(FINALIZE){
    if(spectraFiles.size() == 1 && spectraFiles[0].compare(DEFAULT_FILENAME_SPECFILES) == 0){
      cerr << "ERROR: Must provide spectra files when using --finalize.\n";
      ERROR = true;
    }
    if(LASSI_CHOICE < 1 || LASSI_CHOICE > 5){
      cerr << "ERROR: --lassi-choice must be an integer in {1..5}.\n";
      ERROR = true;
    }
  }

  if (!LASSI){
    cerr << "ERROR: Must specify which statistic to use.\n";
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
  if(LASSI) ending = ".lassi.";

  vector<string> outfiles; 
  if(INIT){

    PopData *popData;
    popData = readPopData(popFilename);
    for (int i = 0; i < popData->popOrder.size(); i++){
      outfiles.push_back(outfileBase + ending + popData->popOrder[i] + ".clrscan.spectra.gz");
    }
    
    // Load VCF data
    map< string, HaplotypeData* > *hapDataByPop;
    hapDataByPop = readHaplotypeDataVCF(vcfFilename, popData); 
    MapData *mapData;
    mapData = hapDataByPop->at(popData->popOrder[0])->map;  

    vector< pair_t* > *windows = findAllWindows(mapData, WINSIZE, WINSTEP);

    cerr << "Calculating haplotype frequency spectra in " << windows->size() << " windows.\n";

    string names = "";

    map<string,double ** > results;
    for(int j = 0;j < popData->popOrder.size(); j++){
      double ** x = new double*[windows->size()];
      for (int i = 0; i < windows->size(); i++) x[i] = new double[K];
      results[popData->popOrder[j]] = x;
    }

    work_order_t *order;
    pthread_t *peer = new pthread_t[numThreads];
    int prev_index = 0;
    for (int i = 0; i < numThreads; i++){
      order = new work_order_t;
      order->id = i;
      //order->numStats = numStats;
      order->hapDataByPop = hapDataByPop;
      //order->mapData = mapData;
      //order->freqData = freqData;
      order->popData = popData;
      //order->flog = &flog;
      //order->bar = &pbar;
      order->params = &params;
      order->results = &results;
      order->windows = windows;
      order->names = &names;
      //order->USE_BP = USE_BP;
      pthread_create(&(peer[i]),
                     NULL,
                     (void *(*)(void *))calc_stats,
                     (void *)order);
    }

    for (int i = 0; i < numThreads; i++) pthread_join(peer[i], NULL);

    delete [] peer;
    cerr << "Done.\n";
    
    ogzstream fout;

    for (int i = 0; i < outfiles.size(); i++){
      string popName = popData->popOrder[i];
      fout.open(outfiles[i].c_str());
      if (fout.fail()) {
        cerr << "ERROR: Failed to open " << outfiles[i] << " for writing.\n";
        return 1;
      }
      fout << "#wins " << windows->size() << " K " << K << endl;
      fout << "chr\tstart\tend\tnSNPs\tnhaps\t" << names << endl;
      for (int w = 0; w < windows->size(); w++) {
        int st = windows->at(w)->start;
        int en = windows->at(w)->end;
        fout << mapData->chr << "\t" 
          << mapData->physicalPos[st] << "\t" 
          << mapData->physicalPos[en] << "\t" 
          << windows->at(w)->end - windows->at(w)->start + 1 << "\t"
          << popData->pop2inds[popName].size()*2;
        for (int s = 0; s < K; s++) {
          double **x = results.at(popName);
          fout << "\t" << x[w][s];
        }
        fout << endl;
      }
      fout.close();
    }
  }
  else{//finalize
    string outfile = outfileBase + ending + "clrscan.out.gz";
    vector<SpectrumData *> *specDataByChr = readSpecData(spectraFiles);
    
    SpectrumData *avgSpec = averageSpec(specDataByChr);
    for (int i = 0; i < avgSpec->K; i++) cerr << avgSpec->freq[0][i] << " ";
    cerr << endl;

    vector<LASSIResults *> *resultsByChr = initResults(specDataByChr);

    work_order2_t *order;
    pthread_t *peer = new pthread_t[numThreads];
    int prev_index = 0;
    for (int i = 0; i < numThreads; i++){
      order = new work_order2_t;
      order->id = i;
      order->specDataByChr = specDataByChr;
      order->avgSpec = avgSpec;
      order->resultsByChr = resultsByChr;
      order->params = &params;

      pthread_create(&(peer[i]),
                     NULL,
                     (void *(*)(void *))calc_stats2,
                     (void *)order);
    }

    for (int i = 0; i < numThreads; i++) pthread_join(peer[i], NULL);

    delete [] peer;

    cerr << "Done.\n";

    ogzstream fout;
    fout.open(outfile.c_str());
    if (fout.fail()) {
      cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
      return 1;
    }
    fout << "chr\tstart\tend\tnSNPs\tnhaps\tm\tT\n";
    for(int c = 0; c < resultsByChr->size(); c++){
      SpectrumData *specData = specDataByChr->at(c);
      LASSIResults *results = resultsByChr->at(c);
      for(int w = 0; w < specData->nwins; w++){
        for(int i = 0; i < 5; i++) fout << specData->info[w][i] << "\t";
        fout << results->m[w] << "\t" << results->T[w] << "\n";
      }
    }

  }

  return 0;
}

