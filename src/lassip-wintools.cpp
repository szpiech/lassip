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

#include "lassip-wintools.h"


void calc_LASSI_stats(void *order) {
	LASSI_work_order_t *p = (LASSI_work_order_t *)order;
	map< string, HaplotypeData* > *hapDataByPop = p->hapDataByPop;
	PopData *popData = p->popData;
	param_t *params = p->params;
	vector< pair_t* > *windows = p->results->windows;
	int WINSIZE = params->getIntFlag(ARG_WINSIZE);
	map<string,double **> *results = p->results->data;
	map<string,double *> *h12ByPop = p->results->h12;
  	map<string,double *> *h2h1ByPop = p->results->h2h1;
	int id = p->id;
	int K = p->params->getIntFlag(ARG_K);
	bool HAPSTATS = p->params->getBoolFlag(ARG_HAPSTATS);
	bool PHASED = !(p->params->getBoolFlag(ARG_UNPHASED));
	map<string,string> *names = p->results->names;

	int numThreads = params->getIntFlag(ARG_THREADS);
	HaplotypeFrequencySpectrum *hfs;
	pair_t *snps;

	string popName;
	for (int p = 0; p < popData->npops; p++){
		popName = popData->popOrder[p];
		for (int i = id; i < windows->size(); i += numThreads) {
			snps = windows->at(i);		
			hfs = hfs_window(hapDataByPop->at(popName), snps);
			double *h12; 
			double *h2h1;
			if(HAPSTATS){
				h12 = h12ByPop->at(popName);
				h2h1 = h2h1ByPop->at(popName);
			}
			double **x = results->at(popName);
			double tot = 0;
			for (int s = 0; s < K; s++){
				if(s < hfs->numClasses) tot+=double(hfs->sortedCount[s]);
			}
			//cerr << tot << endl;
			for (int s = 0; s < K; s++){
				if (i == 0){
					stringstream ss;
					ss << s+1;
					names->at(popName) += popName + "_hfs_" + ss.str(); 
					if (s != K-1) names->at(popName) += "\t";
				}
				if(s < hfs->numClasses) x[i][s] = double(hfs->sortedCount[s])/tot;
				else x[i][s] = 0;
			}
			x[i][K] = hfs->size;
			x[i][K+1] = hfs->hap2count.size();
			if(HAPSTATS){
				h12[i] = calcH12(hfs, PHASED);
				//h2h1[i] = x[i][1]/x[i][0];
				h2h1[i] = calcH2H1(hfs);
			}
			releaseHaplotypeFrequencySpectrum(hfs);
		}
		
	}
	return;
}


void calc_LASSI_stats2(void *order) {
	LASSI_work_order2_t *p = (LASSI_work_order2_t *)order;

	map<string, vector<SpectrumData *>* > *specDataByPopByChr = p->specDataByPopByChr;
    map<string, SpectrumData* > *avgSpecByPop = p->avgSpecByPop;
    map<string, vector<LASSIResults *>* > *resultsByPopByChr = p->resultsByPopByChr;
	param_t *params = p->params;    
	int id = p->id;
	int LASSI_CHOICE = params->getIntFlag(ARG_LASSI_CHOICE); 
	int numThreads = params->getIntFlag(ARG_THREADS);
	int K = avgSpecByPop->begin()->second->K;
	double **f = calcF(LASSI_CHOICE,K);
	
	vector<LASSIResults *> *resultsByChr;
	LASSIResults *results;
	vector<SpectrumData *> *specDataByChr;
	SpectrumData *specData;
	SpectrumData *avgSpec;
	string popName;

	map<string, vector<SpectrumData *>* >::iterator it;
	for(it = specDataByPopByChr->begin(); it != specDataByPopByChr->end(); it++){
		popName = it->first;
		specDataByChr = it->second;
		resultsByChr = resultsByPopByChr->at(popName);
		avgSpec = avgSpecByPop->at(popName);
		for (int c = 0; c < specDataByChr->size(); c++){
			specData = specDataByChr->at(c);
			results = resultsByChr->at(c);
			for (int i = id; i < specDataByChr->at(c)->nwins; i += numThreads) {
				calcMandT(results, specData, avgSpec, f, i);
			}
		}
	}

	for(int i = 0; i < K; i++) delete [] f[i];
	delete [] f;

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