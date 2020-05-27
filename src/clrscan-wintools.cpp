#include "clrscan-wintools.h"


vector< pair_t* > *findAllWindows(MapData *mapData, int WINSIZE, int WINSTEP, bool USE_BP) {
	vector< pair_t* > *windows = new vector< pair_t* >;
	int numSnps = mapData->nloci;
	if (USE_BP){
		int endOfData = mapData->physicalPos[numSnps - 1];	
		int snpIndexStart = 0;

		for (int currWinStart = 0; currWinStart < endOfData; currWinStart += WINSTEP/*, currWinEnd += WINSTEP*/) {	
			//Find SNP index boundaries for the whole window
			pair_t *snps = findInclusiveSNPIndicies(snpIndexStart, currWinStart, WINSIZE, mapData);
			windows->push_back(snps);
			snpIndexStart = snps->start;
		}
	}
	else{//USE_SITES
		for (int i = 0; i < numSnps; i += WINSTEP){
			pair_t* snps = new pair_t;
			snps->start = i;
			if (i+WINSIZE-1 >= numSnps){
				delete snps;
				return windows;
			}
			else{
				snps->end = i+WINSIZE-1;
			}
			snps->winStart = mapData->physicalPos[i];
			windows->push_back(snps);
		}
	}
	return windows;
}

void releaseAllWindows(vector< pair_t* > *windows) {
	for (int i = 0; i < windows->size(); i++) delete windows->at(i);
	delete windows;
	return;
}

void calc_stats(void *order) {
	work_order_t *p = (work_order_t *)order;
	map< string, HaplotypeData* > *hapDataByPop = p->hapDataByPop;
	//MapData *mapData = p->mapData;
	//FreqData *freqData = p->freqData;
	PopData *popData = p->popData;
	param_t *params = p->params;
	vector< pair_t* > *windows = p->windows;
	int WINSIZE = params->getIntFlag(ARG_WINSIZE);
	map<string,double **> *results = p->results;
	int id = p->id;
	//int numStats = p->numStats;
	int K = p->params->getIntFlag(ARG_K);
	string *names = p->names;
	//vector<int> PARTITIONS = params->getIntListFlag(ARG_PARTITION);
	//bool USE_BP = p->USE_BP;


	int numThreads = params->getIntFlag(ARG_THREADS);
	//array_t *sfs, *partition_sfs;
	HaplotypeFrequencySpectrum *hfs;//, *partition_hfs, *pik_hfs;
	pair_t *snps;//, *partition_snps;


	/*
	bool NEED_SFS = false;
	//Do we need to calculate the SFS for every window?
	for (int j = 0; j < NOPTS; j++){
		if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) NEED_SFS = true;
		else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) NEED_SFS = true;
		else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) NEED_SFS = true;
		else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) NEED_SFS = true;
	}
	*/
	//Cycle over all windows and calculate stats

	for (int i = id; i < windows->size(); i += numThreads) {
		snps = windows->at(i);
		string popName;

		for (int p = 0; p < popData->npops; p++){
			popName = popData->popOrder[p];
			//cerr << popName << endl;
			hfs = hfs_window(hapDataByPop->at(popName), snps);
			double tot = 0;
			for (int s = 0; s < K; s++){
				if(s < hfs->numUniq) tot+=double(hfs->sortedCount[s]);
			}
			//cerr << tot << endl;
			for (int s = 0; s < K; s++){
				if (i == 0 && p == 0){
					stringstream ss;
					ss << s+1;
					(*names) += "hfs_" + ss.str(); 
					if (s != K-1) (*names) += "\t";
				}
				double **x = results->at(popName);
				if(s < hfs->numUniq) x[i][s] = double(hfs->sortedCount[s])/tot;
				else x[i][s] = 0;
			}
			releaseHaplotypeFrequencySpectrum(hfs);
		}
		
	}
	return;
}


void calc_stats2(void *order) {
	work_order2_t *p = (work_order2_t *)order;
	vector<SpectrumData *> *specDataByChr = p->specDataByChr;
	param_t *params = p->params;
	SpectrumData *avgSpec = p->avgSpec;
	SpectrumData *specData;
    vector<LASSIResults *> *resultsByChr = p->resultsByChr;
    LASSIResults *results;
	int id = p->id;
	int LASSI_CHOICE = params->getIntFlag(ARG_LASSI_CHOICE); 
	int numThreads = params->getIntFlag(ARG_THREADS);
	
	double **f = calcF(LASSI_CHOICE,avgSpec->K);
	//Cycle over all windows and calculate stats
	for (int c = 0; c < specDataByChr->size(); c++){
		specData = specDataByChr->at(c);
		results = resultsByChr->at(c);
		for (int i = id; i < specDataByChr->at(c)->nwins; i += numThreads) {
			calcMandT(results, specData, avgSpec, f, i);
		}
	}
	for(int i = 0; i < avgSpec->K; i++) delete [] f[i];
	delete [] f;

/*
	for (int i = id; i < windows->size(); i += numThreads) {
		snps = windows->at(i);
		string popName;

		for (int p = 0; p < popData->npops; p++){
			popName = popData->popOrder[p];
			//cerr << popName << endl;
			hfs = hfs_window(hapDataByPop->at(popName), snps);
			double tot = 0;
			for (int s = 0; s < K; s++){
				if(s < hfs->numUniq) tot+=double(hfs->sortedCount[s]);
			}
			//cerr << tot << endl;
			for (int s = 0; s < K; s++){
				if (i == 0 && p == 0){
					stringstream ss;
					ss << s+1;
					(*names) += "hfs_" + ss.str(); 
					if (s != K-1) (*names) += "\t";
				}
				double **x = results->at(popName);
				if(s < hfs->numUniq) x[i][s] = double(hfs->sortedCount[s])/tot;
				else x[i][s] = 0;
			}
			releaseHaplotypeFrequencySpectrum(hfs);
		}
		
	}
	*/
	return;
}

string int2str(int i) {
	char buffer[10];
	sprintf(buffer, "%d", i);
	return string(buffer);
}

/*
vector< pair_t* > *getPartitionWindows(int snpStart, int winStart, vector<int> &PARTITIONS, MapData *mapData, bool USE_BP) {
	vector< pair_t* > *partition_windows = new vector< pair_t* >;
	if (USE_BP){
		int partitionSnpIndexStart = snpStart;
		int partitionCurrWinStart = winStart;
		for (int i = 0; i < PARTITIONS.size(); i++) {
			pair_t *partition_snps = findInclusiveSNPIndicies(partitionSnpIndexStart, partitionCurrWinStart, PARTITIONS[i], mapData);
			partition_windows->push_back(partition_snps);
			partitionSnpIndexStart = partition_snps->end;
			partitionCurrWinStart += PARTITIONS[i];
		}
	}
	else{//USE_SITES
		int numSnps = mapData->nloci;
		int currStart = snpStart;
		int currEnd = -1;
		for (int i = 0; i < PARTITIONS.size(); i++){
			pair_t *partition_snps = new pair_t;
			currEnd = currStart + PARTITIONS[i] - 1;
			partition_snps->start = currStart;
			partition_snps->end = (currEnd >= numSnps) ? numSnps -1 : currEnd;
			partition_windows->push_back(partition_snps);
			currStart = currEnd + 1;
		}
	}
	return partition_windows;
}
*/
/*
vector< pair_t* > *getEHHWindows(int snpStart, int winStart, int WINSIZE, vector<int> &EHH_WINS, MapData *mapData, bool USE_BP) {
	vector< pair_t* > *ehh_windows = new vector< pair_t* >;
	if(USE_BP){
		int currWinStart = winStart;
		for (int i = 0; i < EHH_WINS.size(); i++) {
			pair_t *snps = findInclusiveSNPIndicies(snpStart, ( winStart + (WINSIZE * 0.5) - (EHH_WINS[i] * 0.5) ) , EHH_WINS[i], mapData);
			ehh_windows->push_back(snps);
		}
	}
	else{//USE_SITES
		double mid = (WINSIZE - 1) * 0.5 + snpStart;
		for (int i = 0; i < EHH_WINS.size(); i++) {
			pair_t *snps = new pair_t;
			if (EHH_WINS[i] % 2 == 0){
				snps->start = int(mid - (EHH_WINS[i] * 0.5) + 0.5);
				snps->end = int(mid + (EHH_WINS[i] * 0.5));
			}
			else{
				if (WINSIZE % 2 == 0){
					snps->start = int(mid - (EHH_WINS[i] * 0.5) + 1);
					snps->end = int(mid + (EHH_WINS[i] * 0.5) - 1);
				}
				else{
					snps->start = int(mid - (EHH_WINS[i] * 0.5) + 0.5);
					snps->end = int(mid + (EHH_WINS[i] * 0.5));
				}
			}
			ehh_windows->push_back(snps);
		}
	}
	return ehh_windows;
}
*/

pair_t* findInclusiveSNPIndicies(int startSnpIndex, int currWinStart, int WINSIZE, MapData* mapData) {

	int currWinEnd = currWinStart + WINSIZE - 1;
	int endSnpIndex = startSnpIndex;
	int numSnps = mapData->nloci;

	pair_t* snps = new pair_t;
	snps->winStart = currWinStart;
	if (mapData->physicalPos[numSnps - 1] < currWinStart) {
		snps->start = numSnps;
		snps->end = numSnps - 1;
		return snps;
	}

	while (mapData->physicalPos[startSnpIndex] < currWinStart) {
		startSnpIndex++;
	}
	while (mapData->physicalPos[endSnpIndex] < currWinEnd) {
		endSnpIndex++;
	}
	endSnpIndex--;
	endSnpIndex = (endSnpIndex >= numSnps) ? numSnps - 1 : endSnpIndex;

	snps->start = startSnpIndex;
	snps->end = endSnpIndex;
	return snps;
}