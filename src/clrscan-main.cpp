/* grail -- a program to calculate window-based diversity statistics
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
#include "grail-wintools.h"
#include "grail-winstats.h"
#include "grail-data.h"
#include "grail-cli.h"

using namespace std;

int main(int argc, char *argv[])
{
  cerr << "grail v" + VERSION + "\n";
  param_t params;
  params.setPreamble(PREAMBLE);

  params.addFlag(ARG_THREADS, DEFAULT_THREADS, "", HELP_THREADS);

  // I/O flags
  //params.addFlag(ARG_FILENAME_TPED, DEFAULT_FILENAME_TPED, "", HELP_FILENAME_TPED);
  params.addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
  params.addFlag(ARG_FILENAME_POP1_VCF, DEFAULT_FILENAME_POP1_VCF, "", HELP_FILENAME_POP1_VCF);
  params.addFlag(ARG_FILENAME_POPFILE, DEFAULT_FILENAME_POPFILE, "", HELP_FILENAME_POPFILE);
  //params.addFlag(ARG_FILENAME_MAP, DEFAULT_FILENAME_MAP, "", HELP_FILENAME_MAP);
    
  // Window control flags
  params.addFlag(ARG_BP, DEFAULT_BP, "", HELP_BP);
  params.addFlag(ARG_SITES, DEFAULT_SITES, "", HELP_SITES);
  params.addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
  params.addFlag(ARG_WINSTEP, DEFAULT_WINSTEP, "", HELP_WINSTEP);
  //params.addListFlag(ARG_PARTITION, DEFAULT_PARTITION, "", HELP_PARTITION);

  // Statistics flags
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
  //string mapFilename = params.getStringFlag(ARG_FILENAME_MAP);
  //bool MAP = (mapFilename.compare(DEFAULT_FILENAME_MAP) == 0) ? false : true;
  
  // Window control
  bool USE_BP = params.getBoolFlag(ARG_BP);
  bool USE_SITES = params.getBoolFlag(ARG_SITES);
  int WINSIZE = params.getIntFlag(ARG_WINSIZE);
  int WINSTEP = params.getIntFlag(ARG_WINSTEP);
  //vector<int> PARTITIONS = params.getIntListFlag(ARG_PARTITION);
  //bool DO_PARTITION = false;

  // Statistics
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
  //bool EHH_PART = params.getBoolFlag(ARG_EHH_PART);
  //bool PMAP = params.getBoolFlag(ARG_PMAP);
  //bool SFS_SUB = !(params.getBoolFlag(ARG_NO_SFS_SUB));

  // Check for consistency errors within flags
  bool ERROR = false;

  if (!POP){
    cerr << "ERROR: Must provide a map from ind to pop with --pop.\n";
    ERROR = true;
  }

  if (!USE_SITES && !USE_BP){
    cerr << "ERROR: Must choose to measure windows in either sites or bps.\n";
    ERROR = true;
  }

  if (USE_SITES && USE_BP){
    cerr << "ERROR: Must choose to measure windows in either sites or bps not both.\n";
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

  if (numThreads <= 0) {
    cerr << "ERROR: Must specify a positive number of threads.\n";
    ERROR = true;
  }

  if (/*!TPED &&*/ !VCF) {
    cerr << "ERROR: Must provide a file with genetic data.\n";
    ERROR = true;
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

  string outfile = outfileBase + ".grail.out";
  ofstream fout;
  fout.open(outfile.c_str());
  if (fout.fail()) {
    cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
    return 1;
  }

  PopData *popData;
  popData = readPopData(popFilename);

/*
  cout << "npops = " << popData->npops << endl;
  cout << "ninds = " << popData->nind << endl;
  for (int i = 0; i < popData->popOrder.size(); i++){
    cout << popData->popOrder[i] << " at " << popData->pop2index[popData->popOrder[i]] << endl;
    for (int j = 0; j < popData->pop2inds[popData->popOrder[i]].size(); j++){
      cout << popData->pop2inds[popData->popOrder[i]][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
  for (int i = 0; i < popData->indOrder.size(); i++){
    cout << popData->indOrder[i] << endl;
  }


  return 0;
  */
  map< string, HaplotypeData* > *hapDataByPop;
  //HaplotypeData *hapData;
  
  //FreqData *freqData;
  /*if (TPED){
    hapData = readHaplotypeDataTPED(tpedFilename);
  }
  else */if (VCF){
    hapDataByPop = readHaplotypeDataVCF(vcfFilename, popData); 
  }
  //load physical positions
  findAllAlleles(hapDataByPop, popData);
  /*
cout << "g = " << mapData->g << endl;

for (int i = 0; i < popData->popOrder.size(); i++){
    cout << popData->popOrder[i] << " with " << popData->pop2inds[popData->popOrder[i]].size() << endl;
    hapData = hapDataByPop->at(popData->popOrder[i]);
    for (int j = 0; j < hapData->nloci; j++){
      for (int k = 0; k < hapData->nhaps; k++){
        cout << hapData->data[k][j] << " ";
      }
      map<char,unsigned int>::iterator it;
      for(it = hapData->freq->count[j].begin(); it != hapData->freq->count[j].end(); it++){
        cout << "| " << it->first << " (n=" << it->second << ") ";
      }
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;

  
  for (int i = 0; i < mapData->nloci; i++){
    for (int j = 0; j < mapData->alleles[i].size(); j++){
      cout << mapData->alleles[i][j] << " ";
    }
    cout << endl;
  }

  return 0;
*/
  //if(TPED) mapData = readMapDataTPED(tpedFilename, hapData->nloci, hapData->nhaps);
  //else if (VCF) mapData = readMapDataVCF(vcfFilename, hapData->nloci);

  //freqData = initFreqData(hapData);

  MapData *mapData;
  mapData = hapDataByPop->at(popData->popOrder[0])->map;
  vector< pair_t* > *windows = findAllWindows(mapData, WINSIZE, WINSTEP, USE_BP);

  //
  int numStats = popData->npops * 2; //richness and private richness per pop

  /*(CALC_PI +
                  CALC_PIK * PIK_CHOICE.size() +
                  CALC_S +
                  CALC_TAJ_D +
                  CALC_FAY_WU_H) *
                 (DO_PARTITION * PARTITIONS.size() + 1) +
                 (EHH_PART + EHH_PART * CALC_EHHK * EHHK_CHOICES.size()) * (DO_PARTITION * PARTITIONS.size()) +
                 (CALC_EHH * EHH_WINS.size() +
                  CALC_EHHK * EHHK_CHOICES.size() * EHH_WINS.size());
*/
  cerr << "Calculating " << numStats << " statistics in " << windows->size() << " windows.\n";

  string names = "";
  double **results = new double*[windows->size()];
  for (int i = 0; i < windows->size(); i++) results[i] = new double[numStats];

  work_order_t *order;
  pthread_t *peer = new pthread_t[numThreads];
  int prev_index = 0;
  for (int i = 0; i < numThreads; i++)
  {
    order = new work_order_t;
    order->id = i;
    order->numStats = numStats;
    order->hapDataByPop = hapDataByPop;
    //order->mapData = mapData;
    //order->freqData = freqData;
    order->popData = popData;
    //order->flog = &flog;
    //order->bar = &pbar;
    order->params = &params;
    order->results = results;
    order->windows = windows;
    order->names = &names;
    order->USE_BP = USE_BP;
    pthread_create(&(peer[i]),
                   NULL,
                   (void *(*)(void *))calc_stats,
                   (void *)order);
  }

  for (int i = 0; i < numThreads; i++)
  {
    pthread_join(peer[i], NULL);
  }

  delete [] peer;
  cerr << "Done.\n";
  
  fout << "chr\tstart\tend\tnSNPs\t" << names << endl;
  for (int w = 0; w < windows->size(); w++) {
    fout << mapData->chr << "\t" 
      << windows->at(w)->winStart << "\t" 
      << windows->at(w)->winStart + WINSIZE << "\t" 
      << windows->at(w)->end - windows->at(w)->start + 1;
    for (int s = 0; s < numStats; s++) {
      fout << "\t" << results[w][s];
    }
    fout << endl;
  }

  fout.close();

  //releaseHapData(hapData);
  //releaseMapData(mapData);
  //releaseFreqData(freqData);

  return 0;
}

