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

struct work_order_t
{
    int id;
    //int numStats;
    //HaplotypeData *hapData;
    map< string, HaplotypeData* > *hapDataByPop;
    //MapData *mapData;
    //FreqData *freqData;
    PopData *popData;

	vector< pair_t* > *windows;

    map<string,double ** > *results;
    string *names;
    //ofstream *flog;
    //Bar *bar;

    param_t *params;

    //bool USE_BP;
};

struct work_order2_t
{
    int id;
    vector<SpectrumData *> *specDataByChr;
    SpectrumData *avgSpec;
    vector<LASSIResults *> *resultsByChr;
    param_t *params;
};


pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData);

vector< pair_t* > *findAllWindows(MapData *mapData, int WINSIZE, int WINSTEP, bool USE_BP = false);
void releaseAllWindows(vector< pair_t* > *windows);

//vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData, bool USE_BP);

//vector< pair_t* > *getEHHWindows(int snpStart, int winStart, int WINSIZE, vector<int> &EHH_WINS, MapData *mapData, bool USE_BP);

void calc_stats2(void *work_order);
void calc_stats(void *work_order);

string int2str(int i);

#endif