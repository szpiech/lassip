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
#include <atomic>
#include <thread>
#include "lassip-data.h"
#include "lassip-winstats.h"
#include "param_t.h"
#include "lassip-cli.h"

using namespace std;

//Work is claimed dynamically: each thread takes the next chunk of windows with
//one atomic fetch_add rather than walking a fixed stride (i = id; i += nthreads).
//Window cost varies with the number of SNPs, the number of unique haplotypes and
//the size of the flanking set, so a fixed stride leaves threads idle at the tail.
//One cursor per unit of work: a population in stage 1, a population-contig pair
//in stage 2. Workers visit units in the same order as the caller allocated them.
struct WorkCursor
{
    std::atomic<unsigned int> *next;
    unsigned int nunits;
};

//Enough chunks per thread to balance, few enough that the atomic is amortised.
inline unsigned int chunkFor(unsigned int total, int numThreads)
{
    if (numThreads <= 1) return (total > 0) ? total : 1;
    unsigned int c = total / (unsigned int)(numThreads * 16);
    if (c < 1) c = 1;
    if (c > 64) c = 64;
    return c;
}

inline bool claimChunk(std::atomic<unsigned int> &cursor, unsigned int total,
                       unsigned int chunk, unsigned int &begin, unsigned int &end)
{
    begin = cursor.fetch_add(chunk);
    if (begin >= total) return false;
    end = begin + chunk;
    if (end > total) end = total;
    return true;
}

struct LASSI_work_order_t
{
    int id;
    WorkCursor *cursor;
    //counted per thread and merged after the join; incrementing the shared
    //results->nullWins map from every thread was a data race
    vector<int> nullWins;
    map< string, HaplotypeData* > *hapDataByPop;
    PopData *popData;
    LASSIInitialResults *results;    
    param_t *params;
};

struct LASSI_work_order2_t
{
    int id;
    WorkCursor *cursor;
    map<string, vector<SpectrumData *>* > *specDataByPopByChr;
    map<string, SpectrumData* > *avgSpecByPop;
    map<string, vector<LASSIResults *>* > *resultsByPopByChr;
    param_t *params;
};

//Use this for first phase precomputing the Qs
struct SALTI_work_order_t
{
    int id;
    WorkCursor *cursor;
    SpectrumData *specData;
    SpectrumData *avgSpec;
    double ***q; //e->m->i
    param_t *params;
    LASSIResults *results;
    double dmin;
};



pair_t* findInclusiveSNPIndicies(unsigned int startSnpIndex, unsigned int currWinStart, int WINSIZE, MapData* mapData);

//vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData, bool USE_BP);

//vector< pair_t* > *getEHHWindows(int snpStart, int winStart, int WINSIZE, vector<int> &EHH_WINS, MapData *mapData, bool USE_BP);

void calc_LASSI_stats2(LASSI_work_order2_t *p);
void calc_LASSI_stats(LASSI_work_order_t *p);

void calc_SALTI_stats(SALTI_work_order_t *p);


#endif