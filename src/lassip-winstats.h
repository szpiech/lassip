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

#ifndef __LASSIP_WINSTATS_H__
#define __LASSIP_WINSTATS_H__

#include "lassip-data.h"
#include <map>
#include <vector>
#include <cstdlib>
#include <cmath>

using namespace std;

double getDMin(vector<SpectrumData *> *specDataByChr);

double ***initQ(int K, double U);
void releaseQ(double ***q, int K, double U);
void calcMTA(LASSIResults *results, double ***q, SpectrumData *specData, SpectrumData *avgSpec, int w, double dmin, double MAX_EXTEND);

int compare (const void *a, const void *b);
double calcH12(HaplotypeFrequencySpectrum *hfs, bool PHASED);
double calcH2H1(HaplotypeFrequencySpectrum *hfs);
double **calcF(int type, int K);

int nEpsGrid(int K, double U);
double epsAt(int ei, int K);
void calcQ(double ***q, SpectrumData *avgSpec, double **f);
void calcQ(double *q, SpectrumData *avgSpec, double **f, double U, int m, double e);

void calcMandT(LASSIResults *results, SpectrumData *specData, SpectrumData *avgSpec, double **f, int w);
double calcLASSINullLikelihood(SpectrumData *specData,SpectrumData *avgSpec,int w);
double calcLASSIAltLikelihood(SpectrumData *specData, SpectrumData *avgSpec, double **f, double U, int m, double e, int w);

unsigned int windowSeed(int seed, int start, int end);
HaplotypeFrequencySpectrum *hfs_window(HaplotypeData *hapData, pair_t* snpIndex, double FILTER_HMISS, int MATCH_TOL, int SEED);
void garud_match_haps_w_missing_shuffle(map<string,double> &hap2count,map<string,double> &miss_hap2count, int len, int MATCH_TOL, unsigned int seed);
int garud_ndiff_str(const string &str1, const string &str2, string &str3, int MATCH_TOL);
int numSitesInDataWin(pair_t* win);


#endif