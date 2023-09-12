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

void calcMTA(LASSIResults *results, double ****q, SpectrumData *specData, SpectrumData *avgSpec, int w, double dmin, double MAX_EXTEND);
double calcSALTINullLikelihood(SpectrumData *specData,SpectrumData *avgSpec,int w,int rightLim, int leftLim);
double calcSALTIAltLikelihood(SpectrumData *specData,SpectrumData *avgSpec,double ****q,int e, int m, double A, int w,int rightLim,int leftLim);

int compare (const void *a, const void *b);
int *uniqInt(int *array, int size, int &newSize);
double calcH12(HaplotypeFrequencySpectrum *hfs, bool PHASED);
double calcH2H1(HaplotypeFrequencySpectrum *hfs);
double **calcF(int type, int K);

void calcQ(double ***q, SpectrumData *avgSpec, double **f, int w);
void calcQ(double *q, SpectrumData *avgSpec, double **f, double U, int m, double e, int w);

void calcMandT(LASSIResults *results, SpectrumData *specData, SpectrumData *avgSpec, double **f, int w);
double calcLASSINullLikelihood(SpectrumData *specData,SpectrumData *avgSpec,int w);
double calcLASSIAltLikelihood(SpectrumData *specData, SpectrumData *avgSpec, double **f, double U, int m, double e, int w);

HaplotypeFrequencySpectrum *hfs_window(HaplotypeData *hapData, pair_t* snpIndex);

int numSitesInDataWin(pair_t* win);


#endif