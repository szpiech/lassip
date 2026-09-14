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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include "gzstream.h"
#include <map>
#include <cstdio>

using namespace std;

const double MISSING = -999;
const unsigned int MISSING_UINT = -999;
const char MISSING_CHAR = '?';
const char MISSING_ALLELE = '-';
const string TPED_MISSING = "-9";
const char VCF_MISSING = '.';

class GMapData
{
public:
    GMapData(string filename, double mGap);
    ~GMapData();
    bool getMapInfo(double queryPos, double &gPos, string &locName, string c, int &current_locus);

private:
    map<string,int * > physicalPos;
    map<string,double * > geneticPos;
    map<string,string * > locusName;
    map<string,map<int, int> > ppos2index;
    map<string,int> nloci;
    vector<string> chr;

    double MAXGAP;

    int countFields(const string &str);
    //void initGMapData(int n);
    void initGMapData(vector<string> chrstr, map<string,int> nloci, bool ADD = false);
    double interpolate(double x0, double y0, double x1, double y1, double query);
};

inline double GMapData::interpolate(double x0, double y0, double x1, double y1, double query)
{
    return ( ( (y1 - y0) / (x1 - x0) ) * query + ( y0 - ((y1 - y0) / (x1 - x0)) * x0 ) );
}


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

//Genotypes are drawn from a four-symbol alphabet, so they are stored two bits
//per locus rather than one char: 60 MB becomes 15 MB for the YRI chr22 example,
//and a 2,000-haplotype x 5M-locus dataset needs 2.5 GB instead of 10 GB.
//Rows are haplotype-major and padded by one byte so that a window extraction
//straddling the end of a row can read row[b+1] unconditionally.
const unsigned char GT_0    = 0;
const unsigned char GT_1    = 1;
const unsigned char GT_2    = 2;
const unsigned char GT_MISS = 3;

inline unsigned char gtCode(char allele){
  if (allele == '0') return GT_0;
  if (allele == '1') return GT_1;
  if (allele == '2') return GT_2;
  return GT_MISS;
}

inline char gtChar(unsigned char code){
  if (code == GT_0) return '0';
  if (code == GT_1) return '1';
  if (code == GT_2) return '2';
  return MISSING_ALLELE;
}

inline int gtStride(int nloci){ return (nloci + 3) >> 2; }

inline unsigned char getGT(const unsigned char *row, int locus){
  return (row[locus >> 2] >> ((locus & 3) << 1)) & 3;
}

inline void setGT(unsigned char *row, int locus, unsigned char code){
  int shift = (locus & 3) << 1;
  row[locus >> 2] = (unsigned char)((row[locus >> 2] & ~(3u << shift)) | ((unsigned int)code << shift));
}

//Set a genotype given a pointer to the byte that holds it (used while reading,
//where rows are block lists rather than one contiguous array).
inline void setGTInByte(unsigned char *byte, int locus, unsigned char code){
  int shift = (locus & 3) << 1;
  *byte = (unsigned char)((*byte & ~(3u << shift)) | ((unsigned int)code << shift));
}

//Copy haplen genotypes starting at locus `start` into `out`, itself packed two
//bits per locus, with the unused codes of the final byte cleared so that two
//equal windows always give equal bytes.
inline void extractWindow(const unsigned char *row, int start, int haplen, string &out){
  int nb = (haplen + 3) >> 2;
  out.assign(nb, '\0');
  const unsigned char *p = row + (start >> 2);
  int shift = (start & 3) << 1;
  if (shift == 0){
    for (int b = 0; b < nb; b++) out[b] = (char)p[b];
  }
  else{
    for (int b = 0; b < nb; b++)
      out[b] = (char)((p[b] >> shift) | (unsigned int)(p[b + 1] << (8 - shift)));
  }
  int used = haplen & 3;
  if (used) out[nb - 1] = (char)((unsigned char)out[nb - 1] & (unsigned char)((1u << (used * 2)) - 1));
}

inline int countMissingWindow(const string &packed, int haplen){
  int n = 0;
  for (int i = 0; i < haplen; i++)
    if ((((unsigned char)packed[i >> 2]) >> ((i & 3) << 1) & 3) == GT_MISS) n++;
  return n;
}

//Expand a packed window back to one char per genotype, for the clustering
//path, which compares and rewrites haplotypes site by site.
inline void unpackWindow(const string &packed, int haplen, string &out){
  out.resize(haplen);
  for (int i = 0; i < haplen; i++)
    out[i] = gtChar((((unsigned char)packed[i >> 2]) >> ((i & 3) << 1)) & 3);
}

struct HaplotypeData
{
  unsigned char **data;   //packed two bits per locus; index with getGT/setGT
  int nhaps;
  int nloci;
  int stride;             //bytes per haplotype, excluding the padding byte
  MapData *map;
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
  map<string,double> hap2count;
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
  //unsigned int **info;
  string **info;
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
  map<string,int> *nullWins;
  map<string,double ** > *data;
  map<string,string> *names;
  map<string,double *> *h12;
  map<string,double *> *h2h1;
  map<string,double *> *dist;
};

void writeAverageSpec(string outfileBase, map<string, SpectrumData* > *avgSpecByPop);
bool checkNull(map<string, SpectrumData* > *avgSpecByPop,map<string, vector<SpectrumData *>* > *specDataByPopByChr);
map<string, SpectrumData* > *averageSpec(string nullSpecFile);

vector< pair_t* > *findAllWindows(MapData *mapData, int WINSIZE, int WINSTEP, bool USE_BP = false);
void releaseAllWindows(vector< pair_t* > *windows);

LASSIInitialResults *initResults(map< string, HaplotypeData* > *hapDataByPop, PopData *popData, 
                                int WINSIZE, int WINSTEP, int K, bool HAPSTATS, string DIST_TYPE);
void writeLASSIInitialResults(string outfile, LASSIInitialResults *results, map< string, HaplotypeData* > *hapDataByPop,
                              PopData *popData, int K, bool SPECFILE, bool HAPSTATS, bool PHASED, int FILTER_LEVEL, string DIST_TYPE);

void writeLASSIFinalResults(string outfile, map<string, vector<LASSIResults *>* > *resultsByPopByChr,
                            map<string, vector<SpectrumData *>* > *specDataByPopByChr, bool SALTI);

LASSIResults *initResults(int nwins, bool HAPSTATS, bool SALTI);
vector<LASSIResults *> *initResults(vector<SpectrumData *> *specDataByChr, bool SALTI);
map<string, vector<LASSIResults *>* > *initResults(map<string, vector<SpectrumData *>* > *specDataByPopByChr, bool SALTI);

void fillNWDistance(map<string, vector<SpectrumData *>* > *specDataByPopByChr);
void fillCMDistance(map<string, vector<SpectrumData *>* > *specDataByPopByChr, GMapData &geneticMap);

void releaseResults(LASSIResults *data);

SpectrumData *initSpecData(int nwins, int K, bool doinfo = true, bool HAPSTATS = true);
void releaseSpecData(SpectrumData *data);

map<string, SpectrumData *> *readSpecData(string filename);
map<string, vector<SpectrumData *>* > *readSpecData(vector<string> filenames);

SpectrumData *averageSpec(vector<SpectrumData *> *specDataByChr);
map<string, SpectrumData* > *averageSpec(map<string, vector<SpectrumData *>* > *specDataByPopByChr);
map<string, SpectrumData* > *averageSpec(string nullSpecFile);
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

//allocates the 2-d array and populated it with -9
HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci, bool domap = true);
HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci, bool domap, bool allocRows);
void releaseHapData(HaplotypeData *data);

//reads haplotype data from a VCF, splitting samples into the populations named
//in the pop file; throws on malformed input
map< string, HaplotypeData* > *readHaplotypeDataVCF(string filename, PopData *data, bool PHASED, bool SHARED_MAP);
void accumulateLocusCounts(HaplotypeData *hapData, int *count, int *count2, int *nmissing);
map< string, HaplotypeData* > *compactLoci(map< string, HaplotypeData* > *hapDataByPop, PopData *popData, const char *keep, int keepLoci, bool perPopMap);
map< string, HaplotypeData* > *filterHaplotypeData(map< string, HaplotypeData* > *hapDataByPop, PopData *popData, int FILTER_LEVEL, double FILTER_LMISS, bool KEEP_MONOMORPHIC, bool PHASED);

//counts the number of "fields" in a string
//where a field is defined as a contiguous set of non whitespace
//characters and fields are delimited by whitespace
int countFields(const string &str);

#endif
