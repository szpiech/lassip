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

#ifndef __LASSIP_DATA_H__
#define __LASSIP_DATA_H__
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "gzstream.h"
#include <map>

using namespace std;

const double MISSING = -999;
const char MISSING_CHAR = '?';
const char MISSING_ALLELE = '-';
const string TPED_MISSING = "-9";
const char VCF_MISSING = '.';

struct MapData
{
  unsigned int *physicalPos;
  //double *geneticPos;
  string *locusName;
  int nloci;
  string chr;
  //vector<char> *alleles;
  //int g;
};

struct FreqData
{
  map<char,unsigned int> *count;
  //vector<char> alleles;
  //int *nmissing;
  int nloci;
  int nhaps;
};

struct HaplotypeData
{
  //map<char, double> *Q;
  char **data;
  int nhaps;
  int nloci;
  MapData *map;
  //FreqData *freq;
};

struct array_t
{
  double *data;
  int size;
};

struct PopData
{
  map<string,string> ind2pop;
  map<string, vector<string> > pop2inds;
  vector<string> popOrder;
  vector<string> indOrder;
  map<string,int> pop2index;
  int npops;
  int nind;
};

struct HaplotypeFrequencySpectrum {
  map<string,int> hap2count;
  //multimap<int,string> count2hap;
  int *sortedCount;
  int size;
  int numClasses;
};

struct pair_t //guess it's a triplet...
{
  int start;
  int end;
  int winStart;
};

struct SpectrumData {
  double **freq;
  int nwins;
  int K;
  unsigned int **info;
  unsigned int *nhaps;
  unsigned int *uhaps;
  double *dist;
  double *h12;
  double *h2h1;
  bool HAPSTATS;
  bool PHASED;
};

struct LASSIResults {
  int *m;
  double *A;
  double *T;
  double *h12;
  double *h2h1;
  int nwins;
  bool HAPSTATS;
};

struct LASSIInitialResults{
  map<string,vector< pair_t* > *> *windows;
  map<string,double ** > *data;
  map<string,string> *names;
  map<string,double *> *h12;
  map<string,double *> *h2h1;
  map<string,double *> *dist;
};


double ****initQ(int nwins,int K, double U);
void releaseQ(double ****q, int nwins,int K, double U);

vector< pair_t* > *findAllWindows(MapData *mapData, int WINSIZE, int WINSTEP, bool USE_BP = false);
void releaseAllWindows(vector< pair_t* > *windows);

LASSIInitialResults *initResults(map< string, HaplotypeData* > *hapDataByPop, PopData *popData, 
                                int WINSIZE, int WINSTEP, int K, bool HAPSTATS, string DIST_TYPE);
void writeLASSIInitialResults(string outfile, LASSIInitialResults *results, map< string, HaplotypeData* > *hapDataByPop,
                              PopData *popData, int K, bool SPECFILE, bool HAPSTATS, bool PHASED, int FILTER_LEVEL, string DIST_TYPE);

void writeLASSIFinalResults(string outfile, map<string, vector<LASSIResults *>* > *resultsByPopByChr,
                            map<string, vector<SpectrumData *>* > *specDataByPopByChr);

LASSIResults *initResults(int nwins, bool HAPSTATS, bool SALTI);
vector<LASSIResults *> *initResults(vector<SpectrumData *> *specDataByChr, bool SALTI);
map<string, vector<LASSIResults *>* > *initResults(map<string, vector<SpectrumData *>* > *specDataByPopByChr, bool SALTI);

void releaseResults(LASSIResults *data);

SpectrumData *initSpecData(int nwins, int K, bool doinfo = true, bool HAPSTATS = true);
void releaseSpecData(SpectrumData *data);

map<string, SpectrumData *> *readSpecData(string filename);
map<string, vector<SpectrumData *>* > *readSpecData(vector<string> filenames);

SpectrumData *averageSpec(vector<SpectrumData *> *specDataByChr);
map<string, SpectrumData* > *averageSpec(map<string, vector<SpectrumData *>* > *specDataByPopByChr);

//map<string,char> storeMap();

//void extractAlleleStrs(string gt, string &string1, string &string2);
//void codeAlleles(string string1, string string2, char &allele1, char &allele2);

HaplotypeFrequencySpectrum *initHaplotypeFrequencySpectrum();
void releaseHaplotypeFrequencySpectrum(HaplotypeFrequencySpectrum *data);

array_t *initArray(int size, double fill = 0);
void releaseArray(array_t* data);

PopData *initPopData();
void releasePopData(PopData *data);
PopData *readPopData(string filename);
void checkK(PopData *data, double K);

//allocates the arrays and populates them with -9 or "--" depending on type
MapData *initMapData(int nloci);
void releaseMapData(MapData *data);

//allocates the arrays and populates them with MISSING
FreqData *initFreqData(int nloci);
//FreqData *initFreqData(HaplotypeData* data);
void releaseFreqData(FreqData *data);

//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
//MapData *readMapData(string filename, int expected_loci);
//MapData *readMapDataTPED(string filename, int expected_loci, int expected_haps);
//MapData *readMapDataVCF(string filename, int expected_loci); //Physical positions only

//allocates the 2-d array and populated it with -9
HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci, bool domap = true);
void releaseHapData(HaplotypeData *data);

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
//HaplotypeData *readHaplotypeData(string filename);
//HaplotypeData *readHaplotypeDataTPED(string filename);
HaplotypeData *readHaplotypeDataVCF(string filename);

//vector< HaplotypeData* > *readHaplotypeDataTPED(string filename, PopData *data);
map< string, HaplotypeData* > *readHaplotypeDataVCF(string filename, PopData *data, bool PHASED, bool SHARED_MAP);
//void findAllAlleles(map< string, HaplotypeData* > *hapDataByPop, PopData *popData);
map< string, HaplotypeData* > *filterHaplotypeData(map< string, HaplotypeData* > *hapDataByPop, PopData *popData, int FILTER_LEVEL);


//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);

#endif
