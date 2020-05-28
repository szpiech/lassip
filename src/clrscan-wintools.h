#ifndef __CLRSCAN_WINTOOLS_H__
#define __CLRSCAN_WINTOOLS_H__

#include <vector>
#include <map>
#include <sstream>
#include "clrscan-data.h"
#include "clrscan-winstats.h"
#include "param_t.h"
#include "clrscan-cli.h"

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


pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData);

//vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData, bool USE_BP);

//vector< pair_t* > *getEHHWindows(int snpStart, int winStart, int WINSIZE, vector<int> &EHH_WINS, MapData *mapData, bool USE_BP);

void calc_LASSI_stats2(void *work_order);
void calc_LASSI_stats(void *work_order);

string int2str(int i);

#endif