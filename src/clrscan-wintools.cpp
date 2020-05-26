#include "grail-wintools.h"


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
			snps->end = (i+WINSIZE-1 >= numSnps) ? numSnps - 1 : i+WINSIZE-1;
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
	double **results = p->results;
	int id = p->id;
	int numStats = p->numStats;
	string *names = p->names;
	//vector<int> PARTITIONS = params->getIntListFlag(ARG_PARTITION);
	bool USE_BP = p->USE_BP;


	int numThreads = params->getIntFlag(ARG_THREADS);
	//array_t *sfs, *partition_sfs;
	HaplotypeFrequencySpectrum *hfs, *partition_hfs, *pik_hfs;
	pair_t *snps, *partition_snps;


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
		int s = 0;
		for (int p = 0; p < popData->npops; p++){
			popName = popData->popOrder[p];
			calc_Q(hapDataByPop, popName, snps);
		}
		for (int p = 0; p < popData->npops; p++){
			popName = popData->popOrder[p];
			if (i == 0) (*names) += "alpha_" + popName + "\t";
			//calc alphas
			results[i][s] = calc_alpha(hapDataByPop, popName, snps);
			s++;
			if (i == 0) (*names) += "pi_" + popName;
			if (p != popData->npops-1) (*names) += "\t";
			//cals pis
			results[i][s] = calc_pi(hapDataByPop, popName, snps);
			s++;
			
		}
		//if (NEED_SFS) sfs = sfs_window(freqData, snps, SFS_SUB);

		//int s = 0;
		//int s_pi = MISSING; //note storage location of pi if it exists
		//useful for calculating Taj's D or F&W's H
		//int s_S = MISSING;
		/*
		for (int j = 0; j < NOPTS; j++) {
			if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) {
				if (i == 0) (*names) += "pi ";
				results[i][s] = pi_from_sfs(sfs);
				s_pi = s;
				s++;
			}
			else if (STATS[j].compare(ARG_PIK) == 0 && PIK_CHOICE[0] != 0) {
				pik_hfs = hfs_window(hapData, snps);
				for (int k = 0; k < PIK_CHOICE.size(); k++) {
					if (i == 0) (*names) += "pi" + int2str(PIK_CHOICE[k]) + " ";
					results[i][s] = pi_k2(pik_hfs, PIK_CHOICE[k]);
					s++;
				}
			}
			else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) {
				if (i == 0) (*names) += "S ";
				results[i][s] = segsites(sfs);
				s_S = s;
				s++;
			}
			else if (STATS[j].compare(ARG_EHH) == 0 && EHH_WINS[0] != 0) {
				vector< pair_t* > *ehh_windows = getEHHWindows(snps->start, snps->winStart, WINSIZE, EHH_WINS, mapData, USE_BP);
				for (int w = 0; w < ehh_windows->size(); w++) {
					if (i == 0) (*names) += "ehh_" + int2str(EHH_WINS[w]) + " ";
					hfs = hfs_window(hapData, ehh_windows->at(w));
					results[i][s] = ehh_from_hfs(hfs);
					s++;
					releaseHaplotypeFrequencySpectrum(hfs);
				}
				releaseAllWindows(ehh_windows);
			}
			else if (STATS[j].compare(ARG_EHHK) == 0 && EHHK_CHOICES[0] != 0) {
				vector< pair_t* > *ehh_windows = getEHHWindows(snps->start, snps->winStart, WINSIZE, EHH_WINS, mapData, USE_BP);
				for (int w = 0; w < ehh_windows->size(); w++) {
					hfs = hfs_window(hapData, ehh_windows->at(w));
					for (int k = 0; k < EHHK_CHOICES.size(); k++) {
						if (i == 0) (*names) += "ehh" + int2str(EHHK_CHOICES[k]) + "_" + int2str(EHH_WINS[w]) + " ";
						results[i][s] = ehhk_from_hfs(hfs, EHHK_CHOICES[k]);
						s++;
					}
					releaseHaplotypeFrequencySpectrum(hfs);
				}
				releaseAllWindows(ehh_windows);
			}
			else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) {
				if (i == 0) (*names) += "D ";
				if (s_pi >= 0 && s_S >= 0) results[i][s] = tajimaD_from_sfs(sfs, results[i][s_pi], results[i][s_S]);
				else if (s_pi < 0 && s_S >= 0) results[i][s] = tajimaD_from_sfs(sfs, s_pi, results[i][s_S]);
				else if (s_pi >= 0 && s_S < 0) results[i][s] = tajimaD_from_sfs(sfs, results[i][s_pi], s_S);
				else results[i][s] = tajimaD_from_sfs(sfs);
				s++;
			}
			else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) {
				if (i == 0) (*names) += "H ";
				if (s_pi >= 0) results[i][s] = fayWuH_from_sfs(sfs, results[i][s_pi]);
				else results[i][s] = fayWuH_from_sfs(sfs);
				s++;
			}
		}
		*/
		//releaseArray(sfs);
		/*
		if (DO_PARTITION) {
			char part[2];
			part[0] = 'A';
			part[1] = '\0';
			vector< pair_t* > *partition_windows = getPartitionWindows(snps->start, snps->winStart, PARTITIONS, mapData, USE_BP);
			for (int p = 0; p < partition_windows->size(); p++) {
				int s_pi0 = MISSING;
				int s_S0 = MISSING;
				string partStr(part);
				partition_snps = partition_windows->at(p);
				if (NEED_SFS) partition_sfs = sfs_window(freqData, partition_snps, SFS_SUB);

				for (int j = 0; j < NOPTS; j++) {
					if (STATS[j].compare(ARG_PI) == 0 && params->getBoolFlag(ARG_PI)) {
						if (i == 0) (*names) += "pi_" + partStr + " ";
						results[i][s] = pi_from_sfs(partition_sfs);
						s_pi0 = s;
						s++;
					}
					else if (STATS[j].compare(ARG_PIK) == 0 && PIK_CHOICE[0] != 0) {
						pair_t *shifted_snps = new pair_t;
						shifted_snps->start = partition_snps->start - snps->start;
						shifted_snps->end = partition_snps->end - snps->start;
						for (int k = 0; k < PIK_CHOICE.size(); k++) {
							if (i == 0) (*names) += "pi" +  int2str(PIK_CHOICE[k]) + "_" + partStr + " ";
							results[i][s] = pi_k2(pik_hfs, PIK_CHOICE[k], shifted_snps);
							s++;
						}
						delete shifted_snps;
					}
					else if (STATS[j].compare(ARG_SEGSITES) == 0 && params->getBoolFlag(ARG_SEGSITES)) {
						if (i == 0) (*names) += "S_" + partStr + " ";
						results[i][s] = segsites(partition_sfs);
						s_S0 = s;
						s++;
					}
					else if (STATS[j].compare(ARG_TAJ_D) == 0 && params->getBoolFlag(ARG_TAJ_D)) {
						if (i == 0) (*names) += "D_" + partStr + " ";
						if (s_pi0 >= 0 && s_S0 >= 0) results[i][s] = tajimaD_from_sfs(partition_sfs, results[i][s_pi0], results[i][s_S0]);
						else if (s_pi0 < 0 && s_S0 >= 0) results[i][s] = tajimaD_from_sfs(partition_sfs, s_pi0, results[i][s_S0]);
						else if (s_pi0 >= 0 && s_S0 < 0) results[i][s] = tajimaD_from_sfs(partition_sfs, results[i][s_pi0], s_S0);
						else results[i][s] = tajimaD_from_sfs(partition_sfs);
						s++;
					}
					else if (STATS[j].compare(ARG_FAY_WU_H) == 0 && params->getBoolFlag(ARG_FAY_WU_H)) {
						if (i == 0) (*names) += "H_" + partStr + " ";
						if (s_pi0 >= 0) results[i][s] = fayWuH_from_sfs(partition_sfs, results[i][s_pi0]);
						else results[i][s] = fayWuH_from_sfs(partition_sfs);
						s++;
					}
					else if (STATS[j].compare(ARG_EHH) == 0 && EHH_WINS[0] != 0 && params->getBoolFlag(ARG_EHH_PART)) {
						if (i == 0) (*names) += "ehh_" + partStr + " ";
						hfs = hfs_window(hapData, partition_snps);
						results[i][s] = ehh_from_hfs(hfs);
						s++;
						releaseHaplotypeFrequencySpectrum(hfs);
					}
					else if (STATS[j].compare(ARG_EHHK) == 0 && EHHK_CHOICES[0] != 0 && params->getBoolFlag(ARG_EHH_PART)) {

						hfs = hfs_window(hapData, partition_snps);
						for (int k = 0; k < EHHK_CHOICES.size(); k++) {
							if (i == 0) (*names) += "ehh" + int2str(EHHK_CHOICES[k]) + "_" + partStr + " ";
							results[i][s] = ehhk_from_hfs(hfs, EHHK_CHOICES[k]);
							s++;
						}
						releaseHaplotypeFrequencySpectrum(hfs);

					}
				}
				part[0]++;
				releaseArray(partition_sfs);
			}
		}
		releaseHaplotypeFrequencySpectrum(pik_hfs);
		*/
	}
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