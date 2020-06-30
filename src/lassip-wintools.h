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

#ifndef __LASSIP_WINTOOLS_H__
#define __LASSIP_WINTOOLS_H__

#include <vector>
#include <map>
#include <sstream>
#include "lassip-data.h"
#include "lassip-winstats.h"
#include "param_t.h"
#include "lassip-cli.h"

using namespace std;

struct LASSI_work_order_t
{
    int id;
    map< string, HaplotypeData* > *hapDataByPop;
    PopData *popData;
    LASSIInitialResults *results;    
    param_t *params;
};

struct LASSI_work_order2_t
{
    int id;
    map<string, vector<SpectrumData *>* > *specDataByPopByChr;
    map<string, SpectrumData* > *avgSpecByPop;
    map<string, vector<LASSIResults *>* > *resultsByPopByChr;
    param_t *params;
};

//Use this for first phase precomputing the Qs
struct SALTI_work_order_t
{
    int id;
    SpectrumData *specData;
    SpectrumData *avgSpec;
    double ****q; //win->e->m->i
    param_t *params;
    LASSIResults *results;
    double dmin;
};



pair_t* findInclusiveSNPIndicies(unsigned int startSnpIndex, unsigned int currWinStart, int WINSIZE, MapData* mapData);

//vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData, bool USE_BP);

//vector< pair_t* > *getEHHWindows(int snpStart, int winStart, int WINSIZE, vector<int> &EHH_WINS, MapData *mapData, bool USE_BP);

void calc_LASSI_stats2(void *work_order);
void calc_LASSI_stats(void *work_order);

void calc_SALTI_stats1(void *order);
void calc_SALTI_stats2(void *order);

string int2str(int i);

#endif