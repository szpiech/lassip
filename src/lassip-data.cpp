/* lassip -- a program to calculate haploytpe frequency spetrum statistics
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
#include "lassip-data.h"
#include "lassip-wintools.h"


void writeAverageSpec(string outfileBase, map<string, SpectrumData* > *avgSpecByPop){
    ogzstream fout;
    string outfile = outfileBase + ".lassip.null.spectra.gz";
    fout.open(outfile.c_str());
    if (fout.fail()) {
      cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
      throw 1;
    }
    fout << "#K " << avgSpecByPop->begin()->second->K << " npop " << avgSpecByPop->size() << endl;
    map<string, SpectrumData* >::iterator it;
    for(it = avgSpecByPop->begin(); it != avgSpecByPop->end(); it++){
        fout << it->first;
        for(int k = 0; k < it->second->K; k++) fout << "\t" << it->second->freq[0][k];
        fout << endl;
    }
    return;
}
bool checkNull(map<string, SpectrumData* > *avgSpecByPop,map<string, vector<SpectrumData *>* > *specDataByPopByChr){

    if(avgSpecByPop->begin()->second->K != specDataByPopByChr->begin()->second->at(0)->K){
        cerr << "ERROR: Mismatching K between provided null spectrum and provided empirical spectrum data.\n";
        return false;
    }

    if(avgSpecByPop->size() != specDataByPopByChr->size()){
        cerr << "ERROR: Mismatching populations between provided null spectrum and provided empirical spectrum data.\n";
        return false;
    }
    map<string, SpectrumData* >::iterator it;
    for(it = avgSpecByPop->begin(); it != avgSpecByPop->end(); it++){
        if(specDataByPopByChr->count(it->first) == 0){
            cerr << "ERROR: Mismatching populations between provided null spectrum and provided empirical spectrum data.\n";
            return false;
        }
    }

    map<string, vector<SpectrumData *>* >::iterator it2;
    for(it2 = specDataByPopByChr->begin(); it2 != specDataByPopByChr->end(); it2++){
        if(avgSpecByPop->count(it2->first) == 0){
            cerr << "ERROR: Mismatching populations between provided null spectrum and provided empirical spectrum data.\n";
            return false;
        }
    }

    return true;
}

map<string, SpectrumData* > *averageSpec(string nullSpecFile){
    igzstream fin;
    fin.open(nullSpecFile.c_str());
    if (fin.fail()) {
      cerr << "ERROR: Failed to open " << nullSpecFile << " for writing.\n";
      throw 1;
    }

    //unsigned int nwins = 0;
    int K;
    int npops;
    string junk, popName;
    //SpectrumData *avgSpec = initSpecData(1,K);
    
    stringstream ss;
    getline(fin,junk);
    ss.str(junk);
    //K 10 npop 1
    ss >> junk >> K >> junk >> npops;

    cerr << "Loading null spectrum from " << nullSpecFile << " for npops = " << npops << " K = " << K << endl;

    map<string, SpectrumData* > *avgSpecByPop = new map<string, SpectrumData* >;
    for(int p = 0; p < npops; p++){
        getline(fin,junk);
        ss.clear();
        ss.str(junk);
        ss >> popName;
        avgSpecByPop->operator[](popName) = initSpecData(1,K,false,false);
        for(int k = 0; k < K; k++) ss >> avgSpecByPop->at(popName)->freq[0][k];
    }
    return avgSpecByPop;
}

double ****initQ(int nwins,int K, double U){
    double ****q = new double***[nwins];
    int nEps = int(U*100.0*K);
    for(int w = 0; w < nwins; w++){
        q[w] = new double **[nEps];
        for(int e = 0; e < nEps; e++){
            q[w][e] = new double *[K-1];
            for(int m = 0; m < K-1; m++) q[w][e][m] = new double[K];
        }
    }
    return q;
}

void releaseQ(double ****q, int nwins,int K, double U){
    int nEps = int(U*100.0*K);
    for(int w = 0; w < nwins; w++){
        for(int e = 0; e < nEps; e++){
            for(int m = 0; m < K-1; m++) delete [] q[w][e][m];
            delete [] q[w][e];
        }
        delete [] q[w];
    }
    delete [] q;
    return;
}

map< string, HaplotypeData* > *filterHaplotypeData(map< string, HaplotypeData* > *hapDataByPop, PopData *popData, int FILTER_LEVEL){
    if(FILTER_LEVEL == 1){//filter sites monomorphic across all pops
        int nOriginalLoci = 0;
        int totHaps = 0;
        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            totHaps += hapDataByPop->at(popData->popOrder[p])->nhaps;
            if(p == 0) nOriginalLoci = hapDataByPop->at(popData->popOrder[p])->nloci;
        }
        int *count = new int[nOriginalLoci];
        int *nmissing = new int[nOriginalLoci];
        for(int i = 0; i < nOriginalLoci; i++){
            count[i] = 0;
            nmissing[i] = 0;
        }

        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            string popName = popData->popOrder[p];
            HaplotypeData *hapData = hapDataByPop->at(popName);
            for(int l = 0; l < hapData->nloci; l++){
                for(int h = 0; h < hapData->nhaps; h++){
                    count[l] += (hapData->data[h][l] == '1') ? 1 : 0;
                    nmissing[l] += (hapData->data[h][l] == MISSING_ALLELE) ? 1 : 0;
                }
            }
        }
        int keepLoci = 0;
        for(int i = 0; i < nOriginalLoci; i++){
            if(count[i] > 0 && count[i] < totHaps - nmissing[i]) keepLoci++;
        }

        cerr << "Filtering " << nOriginalLoci - keepLoci << " loci.\n";

        map< string, HaplotypeData* > *newHapDataByPop = new map< string, HaplotypeData* >;
        MapData *oldMapData = hapDataByPop->begin()->second->map;
        MapData *newMapData = initMapData(keepLoci);
        newMapData->chr = oldMapData->chr;
        int l0 = 0;
        for(int l = 0; l < oldMapData->nloci; l++){
            if(count[l] > 0 && count[l] < totHaps - nmissing[l]){
                newMapData->physicalPos[l0] = oldMapData->physicalPos[l];
                newMapData->locusName[l0] = oldMapData->locusName[l];
                l0++;
            }
        }

        releaseMapData(oldMapData);

        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            string popName = popData->popOrder[p];
            HaplotypeData *hapData = hapDataByPop->at(popName);
            newHapDataByPop->operator[](popName) = initHaplotypeData(hapData->nhaps,keepLoci,false);
            l0 = 0;
            for(int l = 0; l < hapData->nloci; l++){
                if(count[l] > 0 && count[l] < totHaps - nmissing[l]){
                    for(int h = 0; h < hapData->nhaps; h++){
                        newHapDataByPop->at(popName)->data[h][l0] = hapData->data[h][l];
                    }
                    l0++;
                }
            }
            newHapDataByPop->at(popName)->map = newMapData;
            hapData->map = NULL;
            releaseHapData(hapData);
        }

        delete [] count;
        delete [] nmissing;
        delete hapDataByPop;
        return newHapDataByPop;

    }
    else{//FILTER_LEVEL == 2, filter sites monomorphic within pops
        map< string, HaplotypeData* > *newHapDataByPop = new map< string, HaplotypeData* >;

        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            string popName = popData->popOrder[p];
            HaplotypeData *hapData = hapDataByPop->at(popName);
            
            int nOriginalLoci = hapDataByPop->at(popName)->nloci;
            int totHaps = hapDataByPop->at(popName)->nhaps;

            int *count = new int[nOriginalLoci];
            int *nmissing = new int[nOriginalLoci];
            for(int i = 0; i < nOriginalLoci; i++){
                count[i] = 0;
                nmissing[i] = 0;
            }

            for(int l = 0; l < hapData->nloci; l++){
                for(int h = 0; h < hapData->nhaps; h++){
                    count[l] += (hapData->data[h][l] == '1') ? 1 : 0;
                    nmissing[l] += (hapData->data[h][l] == MISSING_ALLELE) ? 1 : 0;
                }
            }
            
            int keepLoci = 0;
            for(int i = 0; i < nOriginalLoci; i++){
                if(count[i] > 0 && count[i] < totHaps - nmissing[i]) keepLoci++;
            }

            cerr << "Filtering " << nOriginalLoci - keepLoci << " loci in " << popName << ".\n";
            newHapDataByPop->operator[](popName) = initHaplotypeData(hapData->nhaps,keepLoci,true);
            newHapDataByPop->at(popName)->map->chr = hapData->map->chr;
            int l0 = 0;
            for(int l = 0; l < hapData->nloci; l++){
                if(count[l] > 0 && count[l] < totHaps - nmissing[l]){
                    newHapDataByPop->at(popName)->map->physicalPos[l0] = hapData->map->physicalPos[l];
                    newHapDataByPop->at(popName)->map->locusName[l0] = hapData->map->locusName[l];
                    for(int h = 0; h < hapData->nhaps; h++){
                        newHapDataByPop->at(popName)->data[h][l0] = hapData->data[h][l];
                    }
                    l0++;
                }
            }
            releaseHapData(hapData);
            delete [] count;
            delete [] nmissing;
        }
        delete hapDataByPop;
        return newHapDataByPop;

    }
}


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
    for (unsigned int i = 0; i < windows->size(); i++) delete windows->at(i);
    delete windows;
    return;
}

void writeLASSIFinalResults(string outfile, map<string, vector<LASSIResults *>* > *resultsByPopByChr, map<string, vector<SpectrumData *>* > *specDataByPopByChr, bool SALTI){
    ogzstream fout;
    fout.open(outfile.c_str());
    if (fout.fail()) {
      cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
      throw 1;
    }
    
    bool HAPSTATS = resultsByPopByChr->begin()->second->at(0)->HAPSTATS;
    bool PHASED = specDataByPopByChr->begin()->second->at(0)->PHASED;

    string h12 = "h12";
    string h2h1 = "h2h1";
    if(!PHASED){
        h12 = "g123";
        h2h1 = "g2g1";
    } 
    
    int nchr = 0;
    vector<int> nwins;
    map<string, vector<LASSIResults *>* >::iterator it;
    fout << "chr\tstart\tend\tnSNPs\tpos";
    for(it = resultsByPopByChr->begin(); it != resultsByPopByChr->end(); it++){
        nchr = it->second->size();
        fout << "\t" << it->first << "_nhaps";
        fout << "\t" << it->first << "_uhaps\t";
        if(HAPSTATS){
            fout << it->first << "_" << h12 << "\t" 
                << it->first << "_" << h2h1 << "\t";
        }  
        fout << it->first << "_m\t";
        if(SALTI) fout << it->first << "_A\t";
        fout << it->first << "_T";
    }
    fout << endl;

    for(int c = 0; c < nchr; c++) nwins.push_back(resultsByPopByChr->begin()->second->at(c)->nwins);

    SpectrumData *specData;
    LASSIResults *results;
    unsigned int **info;
    unsigned int *nhaps;
    unsigned int *uhaps;
    double *dist;
    
    for(int c = 0; c < nchr; c++){
        for(int w = 0; w < nwins[c]; w++){
            info = specDataByPopByChr->begin()->second->at(c)->info;
            dist = specDataByPopByChr->begin()->second->at(c)->dist;
            for(int i = 0; i < 4; i++) fout << info[w][i] << "\t";
            fout << setprecision(10) << dist[w] << setprecision(6);
            for(it = resultsByPopByChr->begin(); it != resultsByPopByChr->end(); it++){
                nhaps = specDataByPopByChr->at(it->first)->at(c)->nhaps;
                uhaps = specDataByPopByChr->at(it->first)->at(c)->uhaps;
                results = it->second->at(c);
                fout << "\t" << nhaps[w] << "\t";
                fout << uhaps[w] << "\t";
                if(HAPSTATS){
                    fout << results->h12[w] << "\t"
                        << results->h2h1[w] << "\t";
                }
                fout << results->m[w] << "\t";
                if(SALTI) fout << results->A[w] << "\t";
                fout << results->T[w];
            }
            fout << endl;
        }
    }
    return;
}

void writeLASSIInitialResults(string outfileBase, LASSIInitialResults *results, map< string, HaplotypeData* > *hapDataByPop, PopData *popData, int K, bool SPECFILE, bool HAPSTATS, bool PHASED, int FILTER_LEVEL, string DIST_TYPE){
    string ending, outfile;
    ogzstream fout;
    if(PHASED) ending = ".lassip.hap.";
    if(!PHASED) ending = ".lassip.mlg.";
    string h12 = "h12";
    string h2h1 = "h2h1";
    if(!PHASED){
        h12 = "g123";
        h2h1 = "g2g1";
    }    

    string distStr;
    if(DIST_TYPE.compare("bp") == 0) distStr = "ppos";
    else if (DIST_TYPE.compare("gm") == 0) distStr = "gpos";
    else if (DIST_TYPE.compare("ns") == 0) distStr = "siteNum";
    else if (DIST_TYPE.compare("nw") == 0) distStr = "winNum";

    //bool SPECFILE = LASSI || SALTI;

    if(FILTER_LEVEL < 2){
        if(SPECFILE) outfile = outfileBase + ending + "spectra.gz";
        if(!SPECFILE && HAPSTATS) outfile = outfileBase + ending + "stats.gz";

        fout.open(outfile.c_str());
        if (fout.fail()) {
            cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
            throw 1;
        }
        
        vector< pair_t* > *windows = results->windows->begin()->second;
        MapData *mapData = hapDataByPop->begin()->second->map;

        if (SPECFILE){
            fout << "#phased " << PHASED << " hapstats " << HAPSTATS << " wins " << windows->size() << " K " << K << " npop " << popData->popOrder.size();
            for(unsigned int p = 0; p < popData->popOrder.size(); p++) fout << " " << popData->popOrder[p];
            fout << endl;
        }
        fout << "chr\tstart\tend\tnSNPs\t" << distStr;
        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            fout << "\t" << popData->popOrder[p] << "_nhaps\t" << popData->popOrder[p] << "_uhaps\t";
            if(HAPSTATS) fout << popData->popOrder[p] << "_" << h12 << "\t" << popData->popOrder[p] << "_" << h2h1 << "\t";
            if(SPECFILE) fout << results->names->at(popData->popOrder[p]);
        }
        fout << endl;
    
        double *dist = results->dist->begin()->second;
        
        for (unsigned int w = 0; w < windows->size(); w++) {
            int st = windows->at(w)->start;
            int en = windows->at(w)->end;
            fout << mapData->chr << "\t" 
                << mapData->physicalPos[st] << "\t" 
                << mapData->physicalPos[en] << "\t"
                << windows->at(w)->end - windows->at(w)->start + 1 << "\t"
                << setprecision(10) 
                << dist[w]
                << setprecision(6);
            for(unsigned int p = 0; p < popData->popOrder.size(); p++){
                string popName = popData->popOrder[p];
                double **x = results->data->at(popName);
                double *h12;
                double *h2h1;
                if(HAPSTATS){
                    h12 = results->h12->at(popName);
                    h2h1 = results->h2h1->at(popName);
                }
                fout << "\t" << x[w][K];
                fout << "\t" << x[w][K+1];
                if(HAPSTATS){
                    fout << "\t" << h12[w];
                    fout << "\t" << h2h1[w];
                }
                if(SPECFILE){
                    for (int s = 0; s < K; s++) fout << "\t" << x[w][s];
                }       
            }
            fout << endl;
        }
        fout.close();
    }
    else{
        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            string popName = popData->popOrder[p];
            MapData *mapData = hapDataByPop->at(popName)->map;
            vector< pair_t* > *windows = results->windows->at(popName);

            if(SPECFILE) outfile = outfileBase + "." + popName + ending + "spectra.gz";
            if(!SPECFILE && HAPSTATS) outfile = outfileBase + "." + popName  + ending + "stats.gz";

            fout.open(outfile.c_str());
            if (fout.fail()) {
                cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
                throw 1;
            }

            if (SPECFILE){
                fout << "#phased " << PHASED << " hapstats " << HAPSTATS << " wins " << windows->size() << " K " << K << " npop " << 1;
                fout << " " << popName;
                fout << endl;
            }
        
            fout << "chr\tstart\tend\tnSNPs\t" << distStr;
            fout << "\t" << popName << "_nhaps\t" << popName << "_uhaps\t";
            if(HAPSTATS) fout << popName << "_" << h12 << "\t" << popName << "_" << h2h1 << "\t";
            if(SPECFILE) fout << results->names->at(popName);
            fout << endl;
            
            double *dist = results->dist->at(popName);
            double *h12;
            double *h2h1;
            if(HAPSTATS){
                h12 = results->h12->at(popName);
                h2h1 = results->h2h1->at(popName);
            }
            double **x = results->data->at(popName);

            for (unsigned int w = 0; w < windows->size(); w++) {
                int st = windows->at(w)->start;
                int en = windows->at(w)->end; 
                fout << mapData->chr << "\t" 
                    << mapData->physicalPos[st] << "\t" 
                    << mapData->physicalPos[en] << "\t"
                    << windows->at(w)->end - windows->at(w)->start + 1 << "\t"
                    << setprecision(10) 
                    << dist[w]
                    << setprecision(6)
                    << "\t" << x[w][K]
                    << "\t" << x[w][K+1];
                if(HAPSTATS){
                    fout << "\t" << h12[w];
                    fout << "\t" << h2h1[w];
                }
                if(SPECFILE){
                    for (int s = 0; s < K; s++) fout << "\t" << x[w][s];
                }       
                fout << endl;
            }
            fout.close();
            fout.clear();
        }
    }
    return;
}


LASSIInitialResults *initResults(map< string, HaplotypeData* > *hapDataByPop, PopData *popData, int WINSIZE, int WINSTEP, int K, bool HAPSTATS, string DIST_TYPE){
    

    LASSIInitialResults *results = new LASSIInitialResults;
    results->windows = new map<string,vector< pair_t* > *>;
    results->names = new map<string,string>;
    results->data = new map<string,double ** >;
    if(HAPSTATS){
        results->h12 = new map<string,double *>;
        results->h2h1 = new map<string,double *>;
    }
    results->dist = new map<string,double *>;

    for(unsigned int j = 0;j < popData->popOrder.size(); j++){
        string popName = popData->popOrder[j];
        MapData *mapData = hapDataByPop->at(popName)->map;

        results->windows->operator[](popName) = findAllWindows(mapData, WINSIZE, WINSTEP);
        cerr << "Calculating haplotype frequency spectra in " << results->windows->at(popName)->size() << " windows ";
        cerr << "in pop " << popName << ".\n";

        double ** x = new double*[results->windows->at(popName)->size()];
        double *h12;
        double *h2h1;
        double *dist = new double[results->windows->at(popName)->size()];;

        if(HAPSTATS){
            h12 = new double[results->windows->at(popName)->size()];
            h2h1 = new double[results->windows->at(popName)->size()];
        }
        
        vector< pair_t* > *wins = results->windows->at(popName);
        for (unsigned int i = 0; i < results->windows->at(popName)->size(); i++){
            x[i] = new double[K+2];
            int st = wins->at(i)->start;
            int en = wins->at(i)->end;
            if(DIST_TYPE.compare("bp") == 0) dist[i] = (mapData->physicalPos[en]-mapData->physicalPos[st]+1)*0.5+mapData->physicalPos[st];
            else if(DIST_TYPE.compare("cm") == 0) dist[i] = i; //NEED TO UPDATE TO USING GENETIC MAP
            else if(DIST_TYPE.compare("ns") == 0) dist[i] = (en-st+1)*0.5+st;
            else if(DIST_TYPE.compare("nw") == 0) dist[i] = i;
        }

        results->data->operator[](popName) = x;
        results->names->operator[](popName) = "";

        if(HAPSTATS){
            results->h12->operator[](popName) = h12;
            results->h2h1->operator[](popName) = h2h1;
        }
        results->dist->operator[](popName) = dist;
    }
    return results;
}

LASSIResults *initResults(int nwins, bool HAPSTATS, bool SALTI){
    LASSIResults *data = new LASSIResults;
    data->m = new int[nwins];
    data->nwins = nwins;
    data->T = new double[nwins];
    if(SALTI) data->A = new double[nwins];
    else data->A = NULL;
    if(HAPSTATS){
        data->h12 = new double[nwins];
        data->h2h1 = new double[nwins];
    }
    else{
        data->h12 = NULL;
        data->h2h1 = NULL;
    }
    for(int i = 0; i < nwins; i++){
        data->m[i] = 0;
        data->T[i] = 0;
        if(SALTI) data->A[i] = 0;
        if(HAPSTATS){
            data->h12[i] = 0;
            data->h2h1[i] = 0;
        }
    }
    return data;
}

map<string, vector<LASSIResults *>* > *initResults(map<string, vector<SpectrumData *>* > *specDataByPopByChr, bool SALTI){
    map<string, vector<LASSIResults *>* > *resultsByPopByChr = new map<string, vector<LASSIResults *>* >;
    map<string, vector<SpectrumData *>* >::iterator it;
    for(it = specDataByPopByChr->begin(); it != specDataByPopByChr->end(); it++){
        resultsByPopByChr->operator[](it->first) = initResults(it->second, SALTI);
    }
    return resultsByPopByChr;
}

vector<LASSIResults *> *initResults(vector<SpectrumData *> *specDataByChr, bool SALTI){
    vector<LASSIResults *> *dataByChr = new vector<LASSIResults *>;
    LASSIResults *data;
    for(unsigned int i = 0; i < specDataByChr->size(); i++){
        data = initResults(specDataByChr->at(i)->nwins, false, SALTI);
        //data->dist = specDataByChr->at(i)->dist;
        //specDataByChr->at(i)->dist = NULL;
        if(specDataByChr->at(i)->HAPSTATS){
            data->h12 = specDataByChr->at(i)->h12;
            specDataByChr->at(i)->h12 = NULL;
            data->h2h1 = specDataByChr->at(i)->h2h1;
            specDataByChr->at(i)->h2h1 = NULL;
            data->HAPSTATS = specDataByChr->at(i)->HAPSTATS;
        }
        dataByChr->push_back(data);
    }
    return dataByChr;
}

void releaseResults(LASSIResults *data){
    if (data == NULL) return;
    if (data->m != NULL) delete [] data->m;
    if (data->T != NULL) delete [] data->T;
    if (data->A != NULL) delete [] data->A;
    if (data->h12 != NULL) delete [] data->h12;
    if (data->h2h1 != NULL) delete [] data->h2h1;

    return;
}

//need for each pop
map<string, SpectrumData* > *averageSpec(map<string, vector<SpectrumData *>* > *specDataByPopByChr){
    map<string, SpectrumData* > *avgSpecByPop = new map<string, SpectrumData* >;

    map<string, vector<SpectrumData *>* >::iterator it;
    for(it = specDataByPopByChr->begin(); it != specDataByPopByChr->end(); it++){
        avgSpecByPop->operator[](it->first) = averageSpec(it->second);
    }
    return avgSpecByPop;
}

SpectrumData *averageSpec(vector<SpectrumData *> *specDataByChr){
    int K = specDataByChr->at(0)->K;
    unsigned int nwins = 0;
    SpectrumData *avgSpec = initSpecData(1,K,false,false);
    for (unsigned int i = 0; i < specDataByChr->size(); i++) nwins += specDataByChr->at(i)->nwins;
    for (unsigned int i = 0; i < specDataByChr->size(); i++){
        for (int w = 0; w < specDataByChr->at(i)->nwins; w++){
            for (int j = 0; j < K; j++){
                avgSpec->freq[0][j] += specDataByChr->at(i)->freq[w][j]/nwins;
            }
        }
    }
    return avgSpec;
}

map<string, vector<SpectrumData *>* > *readSpecData(vector<string> filenames){

    map<string, vector<SpectrumData *>* > *specDataByPopByChr = new map<string, vector<SpectrumData *>* >;
    map<string, SpectrumData *> *specDataByPop = readSpecData(filenames[0]);
    
    int K = specDataByPop->begin()->second->K;
    bool PHASED = specDataByPop->begin()->second->PHASED;
    bool HAPSTATS = specDataByPop->begin()->second->HAPSTATS;

    unsigned int npops = specDataByPop->size();
    cerr << "Loading " << filenames[0] << " with " << npops << " pops and K = " << K << endl;

    map<string, SpectrumData *>::iterator it;
    for(it = specDataByPop->begin(); it != specDataByPop->end(); it++){
        specDataByPopByChr->operator[](it->first) = new vector<SpectrumData *>;
        specDataByPopByChr->at(it->first)->push_back(it->second);
    }
    
    for (unsigned int i = 1; i < filenames.size(); i++){
        specDataByPop = readSpecData(filenames[i]);

        cerr << "Loading " << filenames[i] << " with " << specDataByPop->size() 
            << " pops and K = " << specDataByPop->begin()->second->K << endl;
      
        if (K != specDataByPop->begin()->second->K || 
            npops != specDataByPop->size() || 
            PHASED != specDataByPop->begin()->second->PHASED || 
            HAPSTATS != specDataByPop->begin()->second->HAPSTATS){
            
            cerr << "ERROR: Spectra files don't match.\n";
            throw 0;
        }
      
        K = specDataByPop->begin()->second->K;

        for(it = specDataByPop->begin(); it != specDataByPop->end(); it++){
            if(specDataByPopByChr->count(it->first) == 0){
                cerr << "ERROR: Not all files have the same set of populations.\n";
                throw 0;
            }
            specDataByPopByChr->at(it->first)->push_back(it->second);
        }
    }
    return specDataByPopByChr;
}

map<string, SpectrumData *> *readSpecData(string filename){
    igzstream fin;
    stringstream ss;
    string junk, junk0;
    unsigned int nwins;
    int K, npop;
    bool HAPSTATS;
    bool PHASED;
    vector<string> popNames;

    fin.open(filename.c_str());
    if (fin.fail()) {
        cerr << "ERROR: Failed to open " << filename << " for writing.\n";
        throw 0;
    }
    getline(fin,junk);
    ss.str(junk);
    ss >> junk0 >> PHASED >> junk >> HAPSTATS >> junk >> nwins >> junk >> K >> junk >> npop;
    if(junk0.compare("#phased") != 0){
        cerr << "ERROR: Must provide valid spectra files.\n";
        throw 0;
    }
    for(int i = 0; i < npop; i++){
        ss >> junk;
        popNames.push_back(junk);
    }
    ss.clear();

    map<string, SpectrumData *> *specDataByPop = new map<string, SpectrumData *>;
    SpectrumData *data;
    unsigned int **info = new unsigned int*[nwins];
    double *dist = new double[nwins];
    for(int p = 0; p < npop; p++){
        data = initSpecData(nwins,K,false, HAPSTATS);
        data->info = info;
        data->dist = dist;
        data->HAPSTATS = HAPSTATS;
        data->PHASED = PHASED;
        specDataByPop->operator[](popNames[p]) = data;
    }
    
    getline(fin,junk);
    for(unsigned int w = 0; w < nwins; w++){
        info[w] = new unsigned int[4];
        getline(fin,junk);
        ss.str(junk);
        for(int i = 0; i < 4; i++) ss >> info[w][i];
        ss >> dist[w];
        for(int p = 0; p < npop; p++){
            data = specDataByPop->at(popNames[p]);
            ss >> data->nhaps[w];
            ss >> data->uhaps[w];
            if(HAPSTATS){
                ss >> data->h12[w];
                ss >> data->h2h1[w];
            }
            for(int i = 0; i < K; i++) ss >> data->freq[w][i];
        }
        ss.clear();
    }

    fin.close();
    return specDataByPop;
}

SpectrumData *initSpecData(int nwins, int K, bool doinfo, bool HAPSTATS){
    SpectrumData *data = new SpectrumData;
    data->K = K;
    data->nwins = nwins;
    data->freq = new double*[nwins];
    if(doinfo) data->info = new unsigned int*[nwins];
    else data->info = NULL;
    data->nhaps = new unsigned int[nwins];
    data->uhaps = new unsigned int[nwins];
    if(HAPSTATS){
        data->h12 = new double[nwins];
        data->h2h1 = new double[nwins];
    }
    for (int j = 0; j < nwins; j++){
        data->freq[j] = new double[K];
        for (int i = 0; i < K; i++) data->freq[j][i] = 0;
        if(doinfo){
            data->info[j] = new unsigned int[4];
            for (int i = 0; i < 4; i++) data->info[j][i] = 0;
        }
    }
    return data;
}

void releaseSpecData(SpectrumData *data){
    if (data == NULL) return;
    for (int j = 0; j < data->nwins; j++){
        if(data->freq[j] != NULL){
            delete [] data->freq[j];
        }
        if(data->info[j] != NULL){
            delete [] data->info[j];
        }
    }
    if (data->freq != NULL) delete [] data->freq;
    if (data->info != NULL) delete [] data->info;
    if (data->nhaps != NULL) delete [] data->nhaps;
    if (data->uhaps != NULL) delete [] data->uhaps;
    return;
}


HaplotypeFrequencySpectrum *initHaplotypeFrequencySpectrum(){
    HaplotypeFrequencySpectrum *hfs = new HaplotypeFrequencySpectrum;
    hfs->sortedCount = NULL;
    hfs->size = 0;
    hfs->numClasses = 0;
    return hfs;
}

void releaseHaplotypeFrequencySpectrum(HaplotypeFrequencySpectrum *hfs){
    if(hfs == NULL){
        return;
    }

    if(hfs->sortedCount != NULL){
        delete [] hfs->sortedCount;
    }

    delete hfs;
    return;
}

array_t *initArray(int size, double fill){
    array_t *data = new array_t;
    data->size = size;
    data->data = new double[size];
    for(int i = 0; i < size; i++){
        data->data[i] = fill;
    }
    return data;
}
void releaseArray(array_t* data){
    if(data == NULL){
        return;
    }

    if(data->data != NULL){
        delete [] data->data;
    }
    delete data;
    return;
}


PopData *initPopData(){
    PopData *data = new PopData;
    data->npops = 0;
    data->nind = 0;
    return data;
}
void releasePopData(PopData *data){
    delete data;
}


PopData *readPopData(string filename){
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    PopData *popData = initPopData();
    stringstream ss;
    string line, ind, pop;
    while (getline(fin, line)){
        if(countFields(line) != 2){
            cerr << "ERROR: Population file format is <ind ID> <pop ID>.\n";
            throw 0;
        }
        ss.str(line);
        ss >> ind >> pop;
        if(popData->ind2pop.count(ind) == 0){
            popData->ind2pop[ind] = pop;
            popData->nind++;
            popData->indOrder.push_back(ind);
        }
        else{
            cerr << "ERROR: Duplicate individual ID found " << ind << endl;
            throw 0;
        }
        if(popData->pop2inds.count(pop) == 0){
            popData->npops++;
            popData->popOrder.push_back(pop);
            popData->pop2index[pop] = popData->npops-1;
        }
        popData->pop2inds[pop].push_back(ind);
        ss.clear();
    }
    fin.close();

    return popData;
}

void checkK(PopData *popData, double K){
    bool ERROR = false;
    for(unsigned int p = 0; p < popData->popOrder.size(); p++){
      if(popData->pop2inds[popData->popOrder[p]].size() < K){
        cerr << "ERROR: K is greater than total number of haplotypes in " << popData->popOrder[p] << ".\n";
        ERROR = true;
      }
    }
    if(ERROR) throw 1;
    return;
}

FreqData *initFreqData(int nhaps, int nloci) {
    FreqData *freqData = new FreqData;
    freqData->count = new map<char,unsigned int>[nloci];
    for (int i = 0; i < nloci; i++) freqData->count[i][MISSING_ALLELE] = 0;
    /*
    freqData->nmissing = new int[nloci];
    for (int i = 0; i < nloci; i++) {
        freqData->nmissing[i] = 0;
    }
    */
    freqData->nloci = nloci;
    freqData->nhaps = nhaps;
    //freqData->storeAs = storeMap();
    return freqData;
}

/*
FreqData *initFreqData(HaplotypeData* data) {
    if (data == NULL) {
        cerr << "ERROR: Can not compute frequencies on NULL data.\n";
        throw 0;
    }
    FreqData *freqData = initFreqData(data->nhaps,data->nloci);
    //freqData->storeAs = storeMap();

    for (int locus = 0; locus < data->nloci; locus++)
    {
        freqData->count[locus] = 0;
        for (int hap = 0; hap < data->nhaps; hap++)
        {
            freqData->count[locus] += ( data->data[hap][locus] == '1' ? 1 : 0 );
            freqData->nmissing[locus] += ( data->data[hap][locus] == MISSING_ALLELE ? 1 : 0 );
            if (data->data[hap][locus] != '0' && data->data[hap][locus] != '1' && data->data[hap][locus] != MISSING_ALLELE)
            {
                cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                throw 0;
            }
        }
    }

    return freqData;
}
*/
void releaseFreqData(FreqData *data){
    if(data == NULL){
        return;
    }
    if(data->count != NULL){
        delete [] data->count;
    }
    /*
    if(data->nmissing != NULL){
        delete [] data->nmissing;
    }
    */
    delete data;
    return;
}


//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
/*
MapData *readMapData(string filename, int expected_loci)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols << ".\n";
            throw 0;
        }
    }

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci << ".\n";
        throw 0;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading map data for " << nloci << " loci\n";

    MapData *data = initMapData(nloci);

    string chr, junk;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> data->chr;
        fin >> data->locusName[locus];
        fin >> data->geneticPos[locus];
        fin >> data->physicalPos[locus];
        //if (USE_PMAP) data->geneticPos[locus] = data->physicalPos[locus];
    }

    fin.close();
    return data;
}

MapData *readMapDataTPED(string filename, int expected_loci, int expected_haps)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int num_cols = 4;
    int current_cols = 0;
    while (getline(fin, line))
    {
        nloci++;
        current_cols = countFields(line);
        if (current_cols != num_cols + expected_haps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols + expected_haps << ".\n";
            throw 0;
        }
    }

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in map file but found " << nloci << ".\n";
        throw 0;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading map data for " << nloci << " loci\n";

    MapData *data = initMapData(nloci);

    string chr;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> data->chr;
        fin >> data->locusName[locus];
        fin >> data->geneticPos[locus];
        fin >> data->physicalPos[locus];
        getline(fin, line);
    }

    fin.close();
    return data;
}

MapData *readMapDataVCF(string filename, int expected_loci) {
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    string line;
    int nloci = 0;
    int numCommentedLines = 0;
    while (getline(fin, line))
    {
        if (line[0] == '#') {
            numCommentedLines++;
        }
        else {
            nloci++;
        }
    }

    if (nloci != expected_loci)
    {
        cerr << "ERROR: Expected " << expected_loci << " loci in file but found " << nloci << ".\n";
        throw 0;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading map data for " << nloci << " loci\n";

    for (int i = 0; i < numCommentedLines; i++) {
        getline(fin, line);
    }

    MapData *data = initMapData(nloci);

    string chr;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> data->chr;
        fin >> data->physicalPos[locus];
        fin >> data->locusName[locus];
        getline(fin, line);
        data->geneticPos[locus] = data->physicalPos[locus];
    }

    fin.close();
    return data;
}
*/
//allocates the arrays and populates them with MISSING or "--" depending on type

MapData *initMapData(int nloci)
{
    if (nloci < 1)
    {
        cerr << "ERROR: number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    MapData *data = new MapData;
    data->nloci = nloci;
    data->locusName = new string[nloci];
    data->physicalPos = new unsigned int[nloci];
    //data->alleles = new vector<char>[nloci];
    //data->geneticPos = new double[nloci];

    for (int locus = 0; locus < nloci; locus++)
    {
        data->locusName[locus] = "--";
        data->physicalPos[locus] = MISSING_UINT;
        //data->geneticPos[locus] = MISSING;
    }

    return data;
}

void releaseMapData(MapData *data)
{
    if (data == NULL) return;
    data->nloci = -9;
    delete [] data->locusName;
    delete [] data->physicalPos;
    //delete [] data->alleles;
    //delete [] data->geneticPos;
    delete data;
    data = NULL;
    return;
}

//reads in haplotype data and also does basic checks on integrity of format
//returns a populated HaplotypeData structure if successful
//throws an exception otherwise
/*
HaplotypeData *readHaplotypeData(string filename)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    //int fileStart = fin.tellg();
    string line;
    int nhaps = 0;
    int previous_nloci = -1;
    int current_nloci = 0;
    //Counts number of haps (rows) and number of loci (cols)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nhaps++;
        current_nloci = countFields(line);
        //cout << "nloci: " << current_nloci << endl;
        if (previous_nloci < 0)
        {
            previous_nloci = current_nloci;
            continue;
        }
        else if (previous_nloci != current_nloci)
        {
            cerr << "ERROR: line " << nhaps << " of " << filename << " has " << current_nloci
                 << ", but the previous line has " << previous_nloci << ".\n";
            throw 0;
        }
        previous_nloci = current_nloci;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading " << nhaps << " haplotypes and " << current_nloci << " loci...\n";

    HaplotypeData *data = initHaplotypeData(nhaps, current_nloci);


    for (int hap = 0; hap < data->nhaps; hap++)
    {
        for (int locus = 0; locus < data->nloci; locus++)
        {
            fin >> data->data[hap][locus];
            if (data->data[hap][locus] != '0' && data->data[hap][locus] != '1')
            {
                cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                throw 0;
            }
        }
    }

    fin.close();

    return data;
}
*/
/*
HaplotypeData *readHaplotypeDataTPED(string filename)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 4;
    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;
        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nhaps << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading " << current_nhaps - numMapCols << " haplotypes and " << nloci << " loci...\n";

    HaplotypeData *data = initHaplotypeData(current_nhaps - numMapCols, nloci);

    string junk, chr, name, gmap;
    unsigned int pos;
    string allele;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> chr >> name >> gmap >> pos;
        if(locus == 0) data->map->chr = chr;
        data->map->locusName[locus] = name;
        data->map->physicalPos[locus] = pos;
        for (int hap = 0; hap < data->nhaps; hap++)
        {
            fin >> allele;
            if (allele.compare("0") != 0 && allele.compare("1") != 0 && allele.compare(TPED_MISSING) != 0)
            {
                cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                throw 0;
            }
            if (allele.compare(TPED_MISSING) == 0) data->data[hap][locus] = MISSING_ALLELE;
            else data->data[hap][locus] = allele[0];
        }
    }

    fin.close();

    return data;
}
*/
//not done
/*
vector< HaplotypeData* > *readHaplotypeDataTPED(string filename, PopData *popData){
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 4;
    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;
        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nhaps << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    cerr << "Loading " << current_nhaps - numMapCols << " haplotypes and " << nloci << " loci...\n";

    HaplotypeData *data = initHaplotypeData(current_nhaps - numMapCols, nloci);

    string junk, chr, name, gmap;
    unsigned int pos;
    string allele;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> chr >> name >> gmap >> pos;
        if(locus == 0) data->map->chr = chr;
        data->map->locusName[locus] = name;
        data->map->physicalPos[locus] = pos;
        for (int hap = 0; hap < data->nhaps; hap++)
        {
            fin >> allele;
            if (allele.compare("0") != 0 && allele.compare("1") != 0 && allele.compare(TPED_MISSING) != 0)
            {
                cerr << "ERROR:  Alleles must be coded 0/1 only.\n";
                throw 0;
            }
            if (allele.compare(TPED_MISSING) == 0) data->data[hap][locus] = MISSING_ALLELE;
            else data->data[hap][locus] = allele[0];
        }
    }

    fin.close();

    return data;
}
*/
HaplotypeData *readHaplotypeDataVCF(string filename)
{
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 9;
    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        if (line[0] == '#') {
            continue;
        }
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;
        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nhaps << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int nhaps = (current_nhaps - numMapCols) * 2;
    int nfields = (current_nhaps - numMapCols);
    cerr << "Loading " << nhaps << " haplotypes and " << nloci << " loci...\n";

    HaplotypeData *data = initHaplotypeData(nhaps, nloci);

    string junk, chr, name, ref, alt, qual, filter, info, format;
    unsigned int pos;
    char allele1, allele2;
    //bool skipLine = false;
    for (int locus = 0; locus < data->nloci; locus++)
    {
        fin >> chr; 
        if (chr[0] == '#'){
            getline(fin, junk);
            //skipLine = false;
            locus--;
            continue;
        }
        fin >> pos >> name >> ref >> alt >> qual >> filter >> info >> format;
        if (locus == 0) data->map->chr = chr;
        data->map->locusName[locus] = name;
        data->map->physicalPos[locus] = pos;
        for (int field = 0; field < nfields; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            //separator = junk[1];
            allele2 = junk[2];
            if ( (allele1 != '0' && allele1 != '1' && allele1 != VCF_MISSING) || (allele2 != '0' && allele2 != '1' && allele2 != VCF_MISSING) )
            {
                cerr << "ERROR: Alleles must be coded 0/1 only.\n";
                cerr << allele1 << " " << allele2 << endl;
                throw 0;
            }

            if (allele1 == VCF_MISSING) data->data[2 * field][locus] = MISSING_ALLELE;
            else data->data[2 * field][locus] = allele1;
            if (allele2 == VCF_MISSING) data->data[2 * field + 1][locus] = MISSING_ALLELE;
            else data->data[2 * field + 1][locus] = allele2;
        }
    }

    fin.close();

    return data;
}


map< string, HaplotypeData* > *readHaplotypeDataVCF(string filename, PopData *popData, bool PHASED, bool SHARED_MAP){
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail()){
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 9;
    string line;
    int nloci = 0;
    //int previous_nhaps = -1;
    //int current_nhaps = 0;
    while (getline(fin, line))
    {
        if (line[0] == '#') continue;
        else nloci++;
    }
    fin.clear();
    fin.close();

    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    string junk;
    map<string,bool> checkInd;
    
    while(getline(fin, junk)) if(junk[0] == '#' && junk[1] == 'C') break;
    
    int nfields = (countFields(junk) - numMapCols);
    string *inds = new string[nfields];
    stringstream ss;
    ss.str(junk);
    for (int i = 0; i < numMapCols; i++) ss >> junk;
    for (int i = 0; i < nfields; i++){
        ss >> inds[i];
        checkInd[inds[i]] = true;
    }

    for (unsigned int i = 0; i < popData->indOrder.size(); i++){
        if(checkInd.count(popData->indOrder[i]) == 0){
            cerr << "ERROR: " << popData->indOrder[i] << " does not exist in " << filename << endl;
            throw 0;
        }
    }



    int nhaps = nfields;
    if(PHASED) nhaps *= 2;
    
    int nload = 0;
    for (int i = 0; i < popData->npops; i++) nload += (popData->pop2inds[popData->popOrder[i]].size());
    if(PHASED) nload *= 2;
    cerr << "Loading " << nload << "/" << nhaps << " ";
    if(PHASED) cerr << "phased";
    else if (!PHASED) cerr << "unphased"; 
    cerr << " haplotypes with " << nloci << " loci across " << popData->npops << " pops.\n";

    
    MapData *mapData;
    if(SHARED_MAP) mapData = initMapData(nloci);

    map<string,int> pop2indIndex;
    map<string, HaplotypeData* > * dataByPop = new map<string, HaplotypeData* >;
    for (int i = 0; i < popData->npops; i++){
        string popName = popData->popOrder[i];
        if(PHASED) nhaps = (popData->pop2inds[popName].size()) * 2;
        if(!PHASED) nhaps = (popData->pop2inds[popName].size());
        if(SHARED_MAP){
            dataByPop->operator[](popName) = initHaplotypeData(nhaps,nloci,false);
            dataByPop->at(popName)->map = mapData;
        }
        else{
            dataByPop->operator[](popName) = initHaplotypeData(nhaps,nloci,true);
        }
        pop2indIndex[popName] = 0;
    }

    //HaplotypeData *data = initHaplotypeData(nhaps, nloci);

    string chr, name, ref, alt, qual, filter, info, format, alleleStr1, alleleStr2;
    unsigned int pos;
    char allele1, allele2;

    //map<string,char> storeAs = storeMap();

    for (int locus = 0; locus < nloci; locus++)
    {
        //for (int i = 0; i < popData->npops; i++) pop2indIndex[popData->popOrder[i]] = 0;

        fin >> chr >> pos >> name >> ref >> alt >> qual >> filter >> info >> format;
        if(SHARED_MAP){
            for (int i = 0; i < popData->npops; i++) pop2indIndex[popData->popOrder[i]] = 0;
            if (locus == 0) mapData->chr = chr;
            mapData->locusName[locus] = name;
            mapData->physicalPos[locus] = pos;
        }
        else{
            for (int i = 0; i < popData->npops; i++){
                string popName = popData->popOrder[i];
                pop2indIndex[popName] = 0;
                if (locus == 0) dataByPop->at(popName)->map->chr = chr;
                dataByPop->at(popName)->map->locusName[locus] = name;
                dataByPop->at(popName)->map->physicalPos[locus] = pos;
            }
        }
        //cerr << mapData->physicalPos[locus] << " ";
        for (int field = 0; field < nfields; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            allele2 = junk[2];
            if (popData->ind2pop.count(inds[field]) != 0){
                //extractAlleleStrs(junk,alleleStr1,alleleStr2);
                //if(storeAs.count(alleleStr1) == 0 || storeAs.count(alleleStr2) == 0){
                if((allele1 != '1' && allele1 != '0' && allele1 != MISSING_ALLELE) || 
                    (allele2 != '1' && allele2 != '0' && allele2 != MISSING_ALLELE)){
                    cerr << "ERROR: Alleles must be coded 0/1/. only.\n";
                    throw 0;
                }
                //allele1 = storeAs[alleleStr1];
                //allele2 = storeAs[alleleStr2];

                string p = popData->ind2pop[inds[field]];
                int f = pop2indIndex[p];
                //if(dataByPop->at(p)->freq->count[locus].count(allele1) == 0) dataByPop->at(p)->freq->count[locus][allele1] = 1;
                //else dataByPop->at(p)->freq->count[locus][allele1]++;
                //if(dataByPop->at(p)->freq->count[locus].count(allele2) == 0) dataByPop->at(p)->freq->count[locus][allele2] = 1;
                //else dataByPop->at(p)->freq->count[locus][allele2]++;
                //if(allele1 == MISSING_ALLELE) dataByPop->at(p)->freq->nmissing[locus]++;
                //if(allele2 == MISSING_ALLELE) dataByPop->at(p)->freq->nmissing[locus]++;
                if(PHASED){
                    if (allele1 == VCF_MISSING) dataByPop->at(p)->data[2 * f][locus] = MISSING_ALLELE;
                    else dataByPop->at(p)->data[2 * f][locus] = allele1;
                    if (allele2 == VCF_MISSING) dataByPop->at(p)->data[2 * f + 1][locus] = MISSING_ALLELE;
                    else dataByPop->at(p)->data[2 * f + 1][locus] = allele2;
                }
                else if (!PHASED){
                    if(allele1 == '1'){
                        if(allele2 == '1') dataByPop->at(p)->data[f][locus] = '2';
                        else if(allele2 == '0') dataByPop->at(p)->data[f][locus] = '1';
                        else if(allele2 == VCF_MISSING) dataByPop->at(p)->data[f][locus] = MISSING_ALLELE;
                    }
                    else if (allele1 == '0'){
                        if(allele2 == '1') dataByPop->at(p)->data[f][locus] = '1';
                        else if(allele2 == '0') dataByPop->at(p)->data[f][locus] = '0';
                        else if(allele2 == VCF_MISSING) dataByPop->at(p)->data[f][locus] = MISSING_ALLELE;
                    }
                    else if (allele1 == VCF_MISSING) dataByPop->at(p)->data[f][locus] = MISSING_ALLELE;
                }
                pop2indIndex[p]++;
            }
        }
    }

    fin.close();

    return dataByPop;
}

/*
map< string, HaplotypeData* > *readHaplotypeDataVCF(string filename, PopData *popData, bool PHASED){
    igzstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }

    int numMapCols = 9;
    //int fileStart = fin.tellg();
    string line;
    int nloci = 0;
    int previous_nhaps = -1;
    int current_nhaps = 0;
    //Counts number of haps (cols) and number of loci (rows)
    //if any lines differ, send an error message and throw an exception
    while (getline(fin, line))
    {
        if (line[0] == '#') {
            continue;
        }
        //getline(fin,line);
        //if(fin.eof()) break;
        nloci++;
        current_nhaps = countFields(line);
        //cout << "nloci: " << current_nhaps << endl;
        if (previous_nhaps < 0)
        {
            previous_nhaps = current_nhaps;
            continue;
        }
        else if (previous_nhaps != current_nhaps)
        {
            cerr << "ERROR: line " << nloci << " of " << filename << " has " << current_nhaps
                 << " fields, but the previous line has " << previous_nhaps << " fields.\n";
            throw 0;
        }
        previous_nhaps = current_nhaps;
    }

    fin.clear(); // clear error flags
    //fin.seekg(fileStart);
    fin.close();
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw 0;
    }





    int nhaps;
    if(PHASED) nhaps = (current_nhaps - numMapCols) * 2;
    if(!PHASED) nhaps = (current_nhaps - numMapCols);
    int nfields = (current_nhaps - numMapCols);
    int nload = 0;
    for (int i = 0; i < popData->npops; i++) nload += (popData->pop2inds[popData->popOrder[i]].size());
    if(PHASED) nload *= 2;
    cerr << "Loading " << nload << "/" << nhaps << " ";
    if(PHASED) cerr << "phased";
    else if (!PHASED) cerr << "unphased"; 
    cerr << " haplotypes with " << nloci << " loci across " << popData->npops << " pops.\n";

    MapData *mapData = initMapData(nloci);
    //mapData->g = nhaps;
    map<string,int> pop2indIndex;
    map<string, HaplotypeData* > * dataByPop = new map<string, HaplotypeData* >;
    for (int i = 0; i < popData->npops; i++){
        string popName = popData->popOrder[i];
        if(PHASED) nhaps = (popData->pop2inds[popName].size()) * 2;
        if(!PHASED) nhaps = (popData->pop2inds[popName].size());
        dataByPop->operator[](popName) = initHaplotypeData(nhaps,nloci,false);
        dataByPop->at(popName)->map = mapData;
        pop2indIndex[popName] = 0;
    }

    //HaplotypeData *data = initHaplotypeData(nhaps, nloci);

    string junk, chr, name, ref, alt, qual, filter, info, format, alleleStr1, alleleStr2;
    string *inds = new string[nfields];
    unsigned int pos;
    char allele1, allele2, separator;
    map<string,bool> checkInd;
    //bool skipLine = false;
    while(getline(fin, junk)) if(junk[0] == '#' && junk[1] == 'C') break;
    stringstream ss;
    ss.str(junk);
    for (int i = 0; i < 9; i++) ss >> name;
    for (int i = 0; i < nfields; i++){
        ss >> inds[i];
        checkInd[inds[i]] = true;
    }

    for (int i = 0; i < popData->indOrder.size(); i++){
        if(checkInd.count(popData->indOrder[i]) == 0){
            cerr << "ERROR: " << popData->indOrder[i] << " does not exist in " << filename << endl;
            throw 0;
        }
    }

    //map<string,char> storeAs = storeMap();

    for (int locus = 0; locus < nloci; locus++)
    {
        for (int i = 0; i < popData->npops; i++) pop2indIndex[popData->popOrder[i]] = 0;

        fin >> chr >> pos >> name >> ref >> alt >> qual >> filter >> info >> format;
        if (locus == 0) mapData->chr = chr;
        mapData->locusName[locus] = name;
        mapData->physicalPos[locus] = pos;
        //cerr << mapData->physicalPos[locus] << " ";
        for (int field = 0; field < nfields; field++)
        {
            fin >> junk;
            allele1 = junk[0];
            allele2 = junk[2];
            if (popData->ind2pop.count(inds[field]) != 0){
                //extractAlleleStrs(junk,alleleStr1,alleleStr2);
                //if(storeAs.count(alleleStr1) == 0 || storeAs.count(alleleStr2) == 0){
                if((allele1 != '1' && allele1 != '0' && allele1 != MISSING_ALLELE) || 
                    (allele2 != '1' && allele2 != '0' && allele2 != MISSING_ALLELE)){
                    cerr << "ERROR: Alleles must be coded 0/1/. only.\n";
                    throw 0;
                }
                //allele1 = storeAs[alleleStr1];
                //allele2 = storeAs[alleleStr2];

                string p = popData->ind2pop[inds[field]];
                int f = pop2indIndex[p];
                //if(dataByPop->at(p)->freq->count[locus].count(allele1) == 0) dataByPop->at(p)->freq->count[locus][allele1] = 1;
                //else dataByPop->at(p)->freq->count[locus][allele1]++;
                //if(dataByPop->at(p)->freq->count[locus].count(allele2) == 0) dataByPop->at(p)->freq->count[locus][allele2] = 1;
                //else dataByPop->at(p)->freq->count[locus][allele2]++;
                //if(allele1 == MISSING_ALLELE) dataByPop->at(p)->freq->nmissing[locus]++;
                //if(allele2 == MISSING_ALLELE) dataByPop->at(p)->freq->nmissing[locus]++;
                if(PHASED){
                    if (allele1 == VCF_MISSING) dataByPop->at(p)->data[2 * f][locus] = MISSING_ALLELE;
                    else dataByPop->at(p)->data[2 * f][locus] = allele1;
                    if (allele2 == VCF_MISSING) dataByPop->at(p)->data[2 * f + 1][locus] = MISSING_ALLELE;
                    else dataByPop->at(p)->data[2 * f + 1][locus] = allele2;
                }
                else if (!PHASED){
                    if(allele1 == '1'){
                        if(allele2 == '1') dataByPop->at(p)->data[f][locus] = '2';
                        else if(allele2 == '0') dataByPop->at(p)->data[f][locus] = '1';
                        else if(allele2 == VCF_MISSING) dataByPop->at(p)->data[f][locus] = MISSING_ALLELE;
                    }
                    else if (allele1 == '0'){
                        if(allele2 == '1') dataByPop->at(p)->data[f][locus] = '1';
                        else if(allele2 == '0') dataByPop->at(p)->data[f][locus] = '0';
                        else if(allele2 == VCF_MISSING) dataByPop->at(p)->data[f][locus] = MISSING_ALLELE;
                    }
                    else if (allele1 == VCF_MISSING) dataByPop->at(p)->data[f][locus] = MISSING_ALLELE;
                }
                pop2indIndex[p]++;
            }
        }
    }

    fin.close();

    return dataByPop;
}
*/
/*
void extractAlleleStrs(string gt, string &string1, string &string2){
    size_t loc = gt.find('/');
    if (loc == string::npos) loc = gt.find('|');
    if (loc == string::npos){
        cerr << "ERROR: Genotype (" << gt << ") not coded with / or |.\n";
        throw 0;
    }
    string1 = gt.substr(0,loc);
    string2 = gt.substr(loc+1,string::npos);
    return;
}
*/
/*
void findAllAlleles(map< string, HaplotypeData* > *hapDataByPop, PopData *popData){
    map< string, HaplotypeData* >::iterator it;
    string popName;
    FreqData *freqData;
    MapData *mapData = hapDataByPop->at(popData->popOrder[0])->map;
    for (it = hapDataByPop->begin(); it != hapDataByPop->end(); it++) it->second->Q = new map<char, double>[mapData->nloci];
    
    for (int locus = 0; locus < mapData->nloci; locus++){
        
        map<char, bool> a;
        for (it = hapDataByPop->begin(); it != hapDataByPop->end(); it++){
            popName = it->first;
            freqData = it->second->freq;
            map<char,unsigned int>::iterator it2;
            for (it2 = freqData->count[locus].begin(); it2 != freqData->count[locus].end(); it2++){
                if (it2->first != MISSING_ALLELE) a[it2->first] = true;
                else if (mapData->g > freqData->nhaps - it2->second) mapData->g = freqData->nhaps - it2->second;
            }
        }

        map<char, bool>::iterator it3;
        for (it3 = a.begin(); it3 != a.end(); it3++){
            mapData->alleles[locus].push_back(it3->first);
            for (it = hapDataByPop->begin(); it != hapDataByPop->end(); it++){
                it->second->Q[locus][it3->first] = 0;
            }
        }
    }
}
*/
HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci, bool domap)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    HaplotypeData *data = new HaplotypeData;
    data->nhaps = nhaps;
    data->nloci = nloci;

    data->data = new char *[nhaps];
    for (unsigned int i = 0; i < nhaps; i++)
    {
        data->data[i] = new char[nloci];
        for (unsigned int j = 0; j < nloci; j++)
        {
            data->data[i][j] = MISSING_CHAR;
        }
    }

    if (domap) data->map = initMapData(nloci);
    //data->freq = initFreqData(nhaps, nloci);
    //data->Q = NULL;

    return data;
}

void releaseHapData(HaplotypeData *data)
{
    if (data == NULL) return;
    for (int i = 0; i < data->nhaps; i++)
    {
        delete [] data->data[i];
    }

    delete [] data->data;

    if(data->map != NULL) releaseMapData(data->map);
    //if(data->freq != NULL) releaseFreqData(data->freq);
    data->map = NULL;
    data->data = NULL;
    data->nhaps = -9;
    data->nloci = -9;
    delete data;
    data = NULL;
    return;
}


int countFields(const string &str)
{
    string::const_iterator it;
    int result;
    int numFields = 0;
    int seenChar = 0;
    for (it = str.begin() ; it < str.end(); it++)
    {
        result = isspace(*it);
        if (result == 0 && seenChar == 0)
        {
            numFields++;
            seenChar = 1;
        }
        else if (result != 0)
        {
            seenChar = 0;
        }
    }
    return numFields;
}
/*
map<string,char> storeMap(){
    map<string,char> x;
    x["0"] = '0';
    x["1"] = '1';
    x["2"] = '2';
    x["3"] = '3';
    x["4"] = '4';
    x["5"] = '5';
    x["6"] = '6';
    x["7"] = '7';
    x["8"] = '8';
    x["9"] = '9';
    x["10"] = 'a';
    x["11"] = 'b';
    x["12"] = 'c';
    x["13"] = 'd';
    x["14"] = 'e';
    x["15"] = 'f';
    x["16"] = 'g';
    x["17"] = 'h';
    x["18"] = 'i';
    x["19"] = 'j';
    x["20"] = 'k';
    x["21"] = 'l';
    x["22"] = 'm';
    x["23"] = 'n';
    x["24"] = 'o';
    x["25"] = 'p';
    x["26"] = 'q';
    x["27"] = 'r';
    x["28"] = 's';
    x["29"] = 't';
    x["30"] = 'u';
    x["31"] = 'v';
    x["32"] = 'w';
    x["33"] = 'x';
    x["34"] = 'y';
    x["35"] = 'z';
    x["36"] = 'A';
    x["37"] = 'B';
    x["38"] = 'C';
    x["39"] = 'D';
    x["40"] = 'E';
    x["41"] = 'F';
    x["42"] = 'G';
    x["43"] = 'H';
    x["44"] = 'I';
    x["45"] = 'J';
    x["46"] = 'K';
    x["47"] = 'L';
    x["48"] = 'M';
    x["49"] = 'N';
    x["50"] = 'O';
    x["51"] = 'P';
    x["52"] = 'Q';
    x["53"] = 'R';
    x["54"] = 'S';
    x["55"] = 'T';
    x["56"] = 'U';
    x["57"] = 'V';
    x["58"] = 'W';
    x["59"] = 'X';
    x["60"] = 'Y';
    x["61"] = 'Z';
    x["."] = MISSING_ALLELE;
    return x;
}
*/
