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


void calc_LASSI_stats(LASSI_work_order_t *p) {
	map< string, HaplotypeData* > *hapDataByPop = p->hapDataByPop;
	PopData *popData = p->popData;
	param_t *params = p->params;
	
	//int WINSIZE = params->getIntFlag(ARG_WINSIZE);
	int K = p->params->getIntFlag(ARG_K);
	bool HAPSTATS = p->params->getBoolFlag(ARG_HAPSTATS);
	bool PHASED = !(p->params->getBoolFlag(ARG_UNPHASED));

	double FILTER_HMISS = p->params->getDoubleFlag(ARG_FILTER_HMISS);
	int MATCH_TOL = p->params->getIntFlag(ARG_MATCH_TOL);
	int SEED = p->params->getIntFlag(ARG_SEED);
	//resolved once per worker, not per window
	int CLUSTER = clusterMethodCode(p->params->getStringFlag(ARG_HAP_CLUSTER));
	
	int numThreads = params->getIntFlag(ARG_THREADS);
	HaplotypeFrequencySpectrum *hfs;
	pair_t *snps;

	string popName;
	for (int pop = 0; pop < popData->npops; pop++){
		popName = popData->popOrder[pop];
		PopResults &pr = p->results->pops[pop];
		vector< pair_t* > *windows = pr.windows;
		unsigned int nwin = windows->size();
		unsigned int chunk = chunkFor(nwin, numThreads);
		unsigned int begin, end;

		while (claimChunk(p->cursor->next[pop], nwin, chunk, begin, end))
		for (unsigned int i = begin; i < end; i++) {
			snps = windows->at(i);		
			hfs = hfs_window(hapDataByPop->at(popName), snps, FILTER_HMISS, MATCH_TOL, SEED, CLUSTER);
			if(hfs == NULL) p->nullWins[pop]++;
			double **x = pr.data;
			double tot = 0;
			for (int s = 0; s < K; s++){
				if(hfs == NULL) break;
				if(s < hfs->numClasses) tot+=double(hfs->sortedCount[s]);
			}
			//cerr << tot << endl;
			for (int s = 0; s < K; s++){
				if (i == 0){
					stringstream ss;
					ss << s+1;
					pr.header += popName + "_hfs_" + ss.str(); 
					if (s != K-1) pr.header += "\t";
				}
				if(hfs == NULL) x[i][s] = 0;
				else if(s < hfs->numClasses) x[i][s] = double(hfs->sortedCount[s])/tot;
				else x[i][s] = 0;
			}
			if(hfs == NULL){
				x[i][K] = 0;
				x[i][K+1] = 0;
			}
			else{
				x[i][K] = hfs->size;
				//numClasses is the number of distinct haplotypes after any clustering;
				//hap2count is only populated when the clustering path runs
				x[i][K+1] = hfs->numClasses;
			}
			if(HAPSTATS){
				if(hfs == NULL){
					pr.h12[i] = 0;
					pr.h2h1[i] = 0;
				}
				else{
					pr.h12[i] = calcH12(hfs, PHASED);
					pr.h2h1[i] = calcH2H1(hfs);
				}
			}
			releaseHaplotypeFrequencySpectrum(hfs);
		}
		
	}
	return;
}


void calc_LASSI_stats2(LASSI_work_order2_t *p) {
	map<string, vector<SpectrumData *>* > *specDataByPopByChr = p->specDataByPopByChr;
    map<string, SpectrumData* > *avgSpecByPop = p->avgSpecByPop;
    map<string, vector<LASSIResults *>* > *resultsByPopByChr = p->resultsByPopByChr;
	param_t *params = p->params;
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

	unsigned int unit = 0;
	map<string, vector<SpectrumData *>* >::iterator it;
	for(it = specDataByPopByChr->begin(); it != specDataByPopByChr->end(); it++){
		popName = it->first;
		specDataByChr = it->second;
		resultsByChr = resultsByPopByChr->at(popName);
		avgSpec = avgSpecByPop->at(popName);
		for (unsigned int c = 0; c < specDataByChr->size(); c++, unit++){
			specData = specDataByChr->at(c);
			results = resultsByChr->at(c);
			unsigned int nwin = specData->nwins;
			unsigned int chunk = chunkFor(nwin, numThreads);
			unsigned int begin, end;
			while (claimChunk(p->cursor->next[unit], nwin, chunk, begin, end))
				for (unsigned int i = begin; i < end; i++) calcMandT(results, specData, avgSpec, f, i);
		}
	}

	for(int i = 0; i < K; i++) delete [] f[i];
	delete [] f;

	return;
}

void calc_SALTI_stats(SALTI_work_order_t *p) {
	SpectrumData *specData = p->specData;
    SpectrumData *avgSpec = p->avgSpec;
    LASSIResults *results = p->results;
	param_t *params = p->params;
	//int LASSI_CHOICE = params->getIntFlag(ARG_LASSI_CHOICE); 
	int numThreads = params->getIntFlag(ARG_THREADS);
	//validate() rejects any other --dist-type before a thread is started, so the
	//else below is unreachable -- but left uninitialised the variable is
	//indeterminate if that ever stops being true, and calcMTA would silently
	//extend over a garbage distance.
	double MAX_EXTEND = 0;
	string DIST_TYPE = params->getStringFlag(ARG_DIST_TYPE);
	if(DIST_TYPE.compare("bp") == 0){
		MAX_EXTEND = params->getDoubleFlag(ARG_MAX_EXTEND_BP);
	}
	else if(DIST_TYPE.compare("nw") == 0){
		MAX_EXTEND = params->getDoubleFlag(ARG_MAX_EXTEND_NW);
	}
	else if(DIST_TYPE.compare("cm") == 0){
		MAX_EXTEND = params->getDoubleFlag(ARG_MAX_EXTEND_CM);
	}
	else{
		cerr << "ERROR: internal: unhandled " << ARG_DIST_TYPE << " " << DIST_TYPE << ".\n";
		exit(EXIT_INTERNAL);
	}
	
	//int K = avgSpec->K;
	//double **f = calcF(LASSI_CHOICE,K);
	double ***q = p->q;
	//int width = 100;

	unsigned int nwin = specData->nwins;
	unsigned int chunk = chunkFor(nwin, numThreads);
	unsigned int begin, end;
	while (claimChunk(p->cursor->next[0], nwin, chunk, begin, end))
		for (unsigned int i = begin; i < end; i++)
			calcMTA(results, q, specData, avgSpec, i, p->dmin, MAX_EXTEND);

	//for(int i = 0; i < K; i++) delete [] f[i];
	//delete [] f;

	return;
}


pair_t* findInclusiveSNPIndicies(unsigned int startSnpIndex, unsigned int currWinStart, int WINSIZE, MapData* mapData) {

	unsigned int currWinEnd = currWinStart + WINSIZE - 1;
	unsigned int endSnpIndex = startSnpIndex;
	unsigned int numSnps = mapData->nloci;

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