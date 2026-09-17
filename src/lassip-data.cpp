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
#include <cstring>   //memset, memcpy -- libc++ pulls these in transitively, libstdc++ does not
#include <set>      //contig uniqueness while grouping spectra rows

void writeAverageSpec(string outfileBase, map<string, SpectrumData* > *avgSpecByPop){
    ogzstream fout;
    string outfile = outfileBase + ".lassip.null.spectra.gz";
    fout.open(outfile.c_str());
    if (fout.fail()) {
      cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
      throw IOError(outfile, true);
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
      cerr << "ERROR: Failed to open " << nullSpecFile << " for reading.\n";
      throw IOError(nullSpecFile, false);
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

//Per-locus allele and missing-genotype counts, accumulated over one population.
//The genotype matrix is haplotype-major, so the haplotype index has to be the
//outer loop: with l outer, the inner loop strides by nloci bytes and touches a
//new cache line on every iteration.
void accumulateLocusCounts(HaplotypeData *hapData, int *count, int *count2, int *nmissing){
    for(int h = 0; h < hapData->nhaps; h++){
        const unsigned char *row = hapData->data[h];
        int l = 0;
        //one byte holds four genotypes, so the row is read a quarter as often
        for(int b = 0; b < hapData->stride && l < hapData->nloci; b++){
            unsigned int packed = row[b];
            for(int k = 0; k < 4 && l < hapData->nloci; k++, l++, packed >>= 2){
                unsigned int a = packed & 3;
                if(a == GT_1) count[l]++;
                else if(a == GT_2) count2[l]++;
                else if(a == GT_MISS) nmissing[l]++;
            }
        }
    }
    return;
}

//Copy the loci flagged in keep into a freshly sized matrix. perPopMap selects
//between one MapData shared by every population (FILTER_LEVEL 1, where all
//populations keep the same loci) and one MapData per population (FILTER_LEVEL 2).
map< string, HaplotypeData* > *compactLoci(map< string, HaplotypeData* > *hapDataByPop, PopData *popData,
                                           const char *keep, int keepLoci, bool perPopMap){
    map< string, HaplotypeData* > *newHapDataByPop = new map< string, HaplotypeData* >;
    MapData *newMapData = NULL;

    if(!perPopMap){
        MapData *oldMapData = hapDataByPop->begin()->second->map;
        newMapData = initMapData(keepLoci);
        newMapData->chr = oldMapData->chr;
        int l0 = 0;
        for(int l = 0; l < oldMapData->nloci; l++){
            if(keep[l]){
                newMapData->physicalPos[l0] = oldMapData->physicalPos[l];
                newMapData->locusName[l0] = oldMapData->locusName[l];
                l0++;
            }
        }
        releaseMapData(oldMapData);
    }

    for(unsigned int p = 0; p < popData->popOrder.size(); p++){
        string popName = popData->popOrder[p];
        HaplotypeData *hapData = hapDataByPop->at(popName);
        HaplotypeData *newHapData = initHaplotypeData(hapData->nhaps, keepLoci, perPopMap);
        newHapDataByPop->operator[](popName) = newHapData;

        if(perPopMap){
            newHapData->map->chr = hapData->map->chr;
            int l0 = 0;
            for(int l = 0; l < hapData->nloci; l++){
                if(keep[l]){
                    newHapData->map->physicalPos[l0] = hapData->map->physicalPos[l];
                    newHapData->map->locusName[l0] = hapData->map->locusName[l];
                    l0++;
                }
            }
        }
        else{
            newHapData->map = newMapData;
        }

        for(int h = 0; h < hapData->nhaps; h++){
            const unsigned char *src = hapData->data[h];
            unsigned char *dst = newHapData->data[h];
            int l0 = 0;
            for(int l = 0; l < hapData->nloci; l++){
                if(keep[l]) setGT(dst, l0++, getGT(src, l));
            }
        }

        if(!perPopMap) hapData->map = NULL;
        releaseHapData(hapData);
    }

    delete hapDataByPop;
    return newHapDataByPop;
}

//Drop loci whose missing-genotype fraction exceeds FILTER_LMISS and, unless
//KEEP_MONOMORPHIC, loci that carry no variation. Both criteria are evaluated on
//the input matrix -- removing one locus cannot change whether another is
//monomorphic -- so they are applied in a single compaction rather than the two
//full copies of the genotype matrix this used to take.
//
//FILTER_LEVEL 1 pools every population and keeps one locus set; FILTER_LEVEL 2
//evaluates each population separately and lets them keep different loci.
map< string, HaplotypeData* > *filterHaplotypeData(map< string, HaplotypeData* > *hapDataByPop, PopData *popData,
                                                   int FILTER_LEVEL, double FILTER_LMISS, bool KEEP_MONOMORPHIC, bool PHASED){
    if(FILTER_LEVEL == 1){
        int nOriginalLoci = hapDataByPop->begin()->second->nloci;
        int totHaps = 0;
        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            totHaps += hapDataByPop->at(popData->popOrder[p])->nhaps;
        }

        int *count = new int[nOriginalLoci]();
        int *count2 = new int[nOriginalLoci]();
        int *nmissing = new int[nOriginalLoci]();
        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            accumulateLocusCounts(hapDataByPop->at(popData->popOrder[p]), count, count2, nmissing);
        }

        char *keep = new char[nOriginalLoci];
        int nDroppedMissing = 0, nDroppedMono = 0, keepLoci = 0;
        for(int l = 0; l < nOriginalLoci; l++){
            bool passMiss = (double(nmissing[l])/double(totHaps) <= FILTER_LMISS);
            bool polymorphic = (count[l]+count2[l] > 0 && count[l]+count2[l] < totHaps - nmissing[l]);
            if(!passMiss) nDroppedMissing++;
            else if(!KEEP_MONOMORPHIC && !polymorphic) nDroppedMono++;
            keep[l] = (passMiss && (KEEP_MONOMORPHIC || polymorphic)) ? 1 : 0;
            keepLoci += keep[l];
        }

        cerr << "Filtering " << nDroppedMissing << " loci with missing data fraction >" << FILTER_LMISS << ".\n";
        if(!KEEP_MONOMORPHIC) cerr << "Filtering " << nDroppedMono << " monomorphic loci.\n";

        map< string, HaplotypeData* > *out = compactLoci(hapDataByPop, popData, keep, keepLoci, false);
        delete [] count;
        delete [] count2;
        delete [] nmissing;
        delete [] keep;
        //No release here: compactLoci already frees the input it compacts,
        //including the MapData every population shares at this filter level.
        return out;
    }
    else{//FILTER_LEVEL == 2: evaluate, and report, each population separately
        map< string, HaplotypeData* > *newHapDataByPop = new map< string, HaplotypeData* >;

        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            string popName = popData->popOrder[p];
            HaplotypeData *hapData = hapDataByPop->at(popName);
            int nOriginalLoci = hapData->nloci;
            int totHaps = hapData->nhaps;

            int *count = new int[nOriginalLoci]();
            int *count2 = new int[nOriginalLoci]();
            int *nmissing = new int[nOriginalLoci]();
            accumulateLocusCounts(hapData, count, count2, nmissing);

            char *keep = new char[nOriginalLoci];
            int keepLoci = 0;
            for(int l = 0; l < nOriginalLoci; l++){
                bool passMiss = (double(nmissing[l])/double(totHaps) <= FILTER_LMISS);
                bool polymorphic = (count[l]+count2[l] > 0 && count[l]+count2[l] < totHaps - nmissing[l]);
                keep[l] = (passMiss && (KEEP_MONOMORPHIC || polymorphic)) ? 1 : 0;
                keepLoci += keep[l];
            }

            cerr << "Filtering " << nOriginalLoci - keepLoci << " loci in " << popName << ".\n";

            HaplotypeData *newHapData = initHaplotypeData(hapData->nhaps, keepLoci, true);
            newHapDataByPop->operator[](popName) = newHapData;
            newHapData->map->chr = hapData->map->chr;
            int l0 = 0;
            for(int l = 0; l < nOriginalLoci; l++){
                if(keep[l]){
                    newHapData->map->physicalPos[l0] = hapData->map->physicalPos[l];
                    newHapData->map->locusName[l0] = hapData->map->locusName[l];
                    l0++;
                }
            }
            for(int h = 0; h < hapData->nhaps; h++){
                const unsigned char *src = hapData->data[h];
                unsigned char *dst = newHapData->data[h];
                l0 = 0;
                for(int l = 0; l < nOriginalLoci; l++){
                    if(keep[l]) setGT(dst, l0++, getGT(src, l));
                }
            }

            releaseHapData(hapData);
            delete [] count;
            delete [] count2;
            delete [] nmissing;
            delete [] keep;
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
      throw IOError(outfile, true);
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
        if(!SALTI) fout << it->first << "_T";
        if(SALTI) fout << it->first << "_L";
    }
    fout << endl;

    for(int c = 0; c < nchr; c++) nwins.push_back(resultsByPopByChr->begin()->second->at(c)->nwins);

    //SpectrumData *specData;
    LASSIResults *results;
    string **info;
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

namespace {

//Rows this contig will actually contribute: the windows that survive the same
//test the row loop applies. The header count has to be exact, because a spectra
//file whose declared window count disagrees with its rows is refused on read.
unsigned int writtenRows(const LASSIInitialResults *r, int K, int onlyPop){
    //At --filter-level 2 each population has its own window list, so the count
    //must be taken over that population's windows, not population 0's.
    const vector< pair_t* > *windows = r->pops[onlyPop >= 0 ? onlyPop : 0].windows;
    unsigned int n = 0;
    for (unsigned int w = 0; w < windows->size(); w++){
        bool skip = false;
        if(onlyPop >= 0){
            if(r->pops[onlyPop].data[w][K] == 0) skip = true;
        }
        else{
            for (unsigned int p = 0; p < r->pops.size(); p++)
                if(r->pops[p].data[w][K] == 0) skip = true;
        }
        if(!skip) n++;
    }
    return n;
}

//'contigs <n> <name> <rows> ...', appended after the population names so that
//the positional header parse in earlier versions simply stops before it. Only
//written when there is more than one contig, which keeps single-contig output
//byte-identical to what lassip has always produced.
string contigField(const vector<LASSIInitialResults *> &byContig, int K, int onlyPop){
    if(byContig.size() < 2) return "";
    stringstream ss;
    ss << " contigs " << byContig.size();
    for (unsigned int c = 0; c < byContig.size(); c++)
        ss << " " << byContig[c]->pops[0].chr << " " << writtenRows(byContig[c], K, onlyPop);
    return ss.str();
}

} //namespace

void writeLASSIInitialResults(string outfileBase, const vector<LASSIInitialResults *> &resultsByContig, PopData *popData, int K, bool SPECFILE, bool HAPSTATS, bool PHASED, int FILTER_LEVEL, string DIST_TYPE){
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
    //The stage-1 column is always the window's physical midpoint; --dist-type
    //only changes what stage 2 measures distance in. Label it accordingly
    //rather than leaving the header blank, which is what --dist-type cm did:
    //the old table answered to gm and ns, which are not accepted values, and
    //had no entry for cm.
    if(DIST_TYPE.compare("nw") == 0) distStr = "winNum";
    else distStr = "ppos";

    if(FILTER_LEVEL < 2){
        if(SPECFILE) outfile = outfileBase + ending + "spectra.gz";
        if(!SPECFILE && HAPSTATS) outfile = outfileBase + ending + "stats.gz";

        fout.open(outfile.c_str());
        if (fout.fail()) {
            cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
            throw IOError(outfile, true);
        }

        if (SPECFILE){
            unsigned int total = 0;
            for (unsigned int c = 0; c < resultsByContig.size(); c++)
                total += writtenRows(resultsByContig[c], K, -1);
            fout << "#phased " << PHASED << " hapstats " << HAPSTATS << " wins " << total << " K " << K << " npop " << popData->popOrder.size();
            for(unsigned int p = 0; p < popData->popOrder.size(); p++) fout << " " << popData->popOrder[p];
            fout << contigField(resultsByContig, K, -1);
            fout << endl;
        }
        fout << "chr\tstart\tend\tnSNPs\t" << distStr;
        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            fout << "\t" << popData->popOrder[p] << "_nhaps\t" << popData->popOrder[p] << "_uhaps\t";
            if(HAPSTATS) fout << popData->popOrder[p] << "_" << h12 << "\t" << popData->popOrder[p] << "_" << h2h1 << "\t";
            if(SPECFILE) fout << resultsByContig[0]->pops[p].header;
        }
        fout << endl;

        for (unsigned int c = 0; c < resultsByContig.size(); c++){
            const LASSIInitialResults *results = resultsByContig[c];
            const PopResults &ref = results->pops[0];
            for (unsigned int w = 0; w < ref.windows->size(); w++) {
                bool skip = false;
                for(unsigned int p = 0; p < results->pops.size(); p++){
                    if(results->pops[p].data[w][K] == 0) skip = true;
                }

                if(skip) continue;

                fout << ref.chr << "\t" 
                    << ref.startPos[w] << "\t" 
                    << ref.endPos[w] << "\t"
                    << ref.nsnps[w] << "\t"
                    << setprecision(10) 
                    << ref.dist[w]
                    << setprecision(6);
                for(unsigned int p = 0; p < results->pops.size(); p++){
                    double **x = results->pops[p].data;
                    double *h12v = results->pops[p].h12;
                    double *h2h1v = results->pops[p].h2h1;
                    fout << "\t" << x[w][K];
                    fout << "\t" << x[w][K+1];
                    if(HAPSTATS){
                        fout << "\t" << h12v[w];
                        fout << "\t" << h2h1v[w];
                    }
                    if(SPECFILE){
                        for (int s = 0; s < K; s++) fout << "\t" << x[w][s];
                    }       
                }
                fout << endl;
            }
        }
        fout.close();
    }
    else{
        for(unsigned int p = 0; p < popData->popOrder.size(); p++){
            string popName = resultsByContig[0]->pops[p].name;

            if(SPECFILE) outfile = outfileBase + "." + popName + ending + "spectra.gz";
            if(!SPECFILE && HAPSTATS) outfile = outfileBase + "." + popName  + ending + "stats.gz";

            fout.open(outfile.c_str());
            if (fout.fail()) {
                cerr << "ERROR: Failed to open " << outfile << " for writing.\n";
                throw IOError(outfile, true);
            }

            if (SPECFILE){
                unsigned int total = 0;
                for (unsigned int c = 0; c < resultsByContig.size(); c++)
                    total += writtenRows(resultsByContig[c], K, p);
                fout << "#phased " << PHASED << " hapstats " << HAPSTATS << " wins " << total << " K " << K << " npop " << 1;
                fout << " " << popName;
                fout << contigField(resultsByContig, K, p);
                fout << endl;
            }
        
            fout << "chr\tstart\tend\tnSNPs\t" << distStr;
            fout << "\t" << popName << "_nhaps\t" << popName << "_uhaps\t";
            if(HAPSTATS) fout << popName << "_" << h12 << "\t" << popName << "_" << h2h1 << "\t";
            if(SPECFILE) fout << resultsByContig[0]->pops[p].header;
            fout << endl;

            for (unsigned int c = 0; c < resultsByContig.size(); c++){
                const PopResults &pr = resultsByContig[c]->pops[p];
                double *dist = pr.dist;
                double *h12v = pr.h12;
                double *h2h1v = pr.h2h1;
                double **x = pr.data;

                for (unsigned int w = 0; w < pr.windows->size(); w++) {
                    if(x[w][K] == 0) continue;

                    fout << pr.chr << "\t" 
                        << pr.startPos[w] << "\t" 
                        << pr.endPos[w] << "\t"
                        << pr.nsnps[w] << "\t"
                        << setprecision(10) 
                        << dist[w]
                        << setprecision(6)
                        << "\t" << x[w][K]
                        << "\t" << x[w][K+1];
                    if(HAPSTATS){
                        fout << "\t" << h12v[w];
                        fout << "\t" << h2h1v[w];
                    }
                    if(SPECFILE){
                        for (int s = 0; s < K; s++) fout << "\t" << x[w][s];
                    }       
                    fout << endl;
                }
            }
            fout.close();
            fout.clear();
        }
    }
    return;
}



LASSIInitialResults *initResults(map< string, HaplotypeData* > *hapDataByPop, PopData *popData, int WINSIZE, int WINSTEP, int K, bool HAPSTATS, string DIST_TYPE){

    LASSIInitialResults *results = new LASSIInitialResults;
    results->pops.resize(popData->popOrder.size());

    for(unsigned int j = 0; j < popData->popOrder.size(); j++){
        PopResults &pr = results->pops[j];
        pr.name = popData->popOrder[j];
        MapData *mapData = hapDataByPop->at(pr.name)->map;

        pr.windows = findAllWindows(mapData, WINSIZE, WINSTEP);
        unsigned int nwin = pr.windows->size();
        cerr << "Calculating haplotype frequency spectra in " << nwin << " windows ";
        cerr << "in pop " << pr.name << ".\n";

        pr.data = new double*[nwin];
        pr.dist = new double[nwin];
        pr.h12 = HAPSTATS ? new double[nwin] : NULL;
        pr.h2h1 = HAPSTATS ? new double[nwin] : NULL;
        pr.header = "";
        pr.nullWins = 0;

        pr.chr = mapData->chr;
        pr.startPos = new unsigned int[nwin];
        pr.endPos = new unsigned int[nwin];
        pr.nsnps = new int[nwin];

        for (unsigned int i = 0; i < nwin; i++){
            pr.data[i] = new double[K+2];
            int st = pr.windows->at(i)->start;
            int en = pr.windows->at(i)->end;
            pr.dist[i] = (mapData->physicalPos[en]-mapData->physicalPos[st]+1)*0.5+mapData->physicalPos[st];
            pr.startPos[i] = mapData->physicalPos[st];
            pr.endPos[i] = mapData->physicalPos[en];
            pr.nsnps[i] = en - st + 1;
        }
    }

    return results;
}

LASSIResults *initResults(int nwins, bool HAPSTATS, bool SALTI){
    LASSIResults *data = new LASSIResults;
    data->m = new int[nwins];
    data->nwins = nwins;
    data->T = new double[nwins];
    data->HAPSTATS = HAPSTATS;
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

namespace {

//Everything on a spectra file's '#phased ...' line.
struct SpecHeader {
    bool PHASED;
    bool HAPSTATS;
    unsigned int nwins;
    int K;
    int npop;
    vector<string> popNames;
    //Optional trailing 'contigs <n> <name> <nwins> ...', written by lassip 1.3+
    //for a file holding more than one contig. Absent from older files and from
    //single-contig files, where the layout is taken from the chr column alone.
    //The field sits after the population names because the parse is positional:
    //appending to the line leaves it readable by earlier versions.
    bool haveContigs;
    vector<string> contigNames;
    vector<unsigned int> contigWins;
};

SpecHeader parseSpecHeader(const string &line, const string &filename){
    SpecHeader h;
    stringstream ss(line);
    string junk, tag;
    ss >> tag >> h.PHASED >> junk >> h.HAPSTATS >> junk >> h.nwins >> junk >> h.K >> junk >> h.npop;
    if(tag.compare("#phased") != 0 || ss.fail()){
        cerr << "ERROR: " << filename << " is not a valid spectra file.\n";
        throw 0;
    }
    for(int i = 0; i < h.npop; i++){
        string name;
        if(!(ss >> name)){
            cerr << "ERROR: " << filename << " header names fewer than " << h.npop << " populations.\n";
            throw 0;
        }
        h.popNames.push_back(name);
    }
    h.haveContigs = false;
    string next;
    if(ss >> next && next.compare("contigs") == 0){
        int ncontig = 0;
        ss >> ncontig;
        for(int i = 0; i < ncontig; i++){
            string name;
            unsigned int n = 0;
            if(!(ss >> name >> n)){
                cerr << "ERROR: " << filename << " header declares " << ncontig
                     << " contigs but lists " << i << ".\n";
                throw 0;
            }
            h.contigNames.push_back(name);
            h.contigWins.push_back(n);
        }
        h.haveContigs = true;
    }
    return h;
}

//First pass over the rows, reading the chr and start columns only, to learn the
//contig layout before anything is allocated. Also where a file that was
//concatenated wrongly is caught: rows must be grouped by contig and ascending
//within one, because the flanking scan walks by index and assumes dist is
//monotone.
vector< pair<string, unsigned int> > scanContigRuns(const string &filename, unsigned int nwins){
    igzstream fin;
    fin.open(filename.c_str());
    if(fin.fail()){
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw IOError(filename, false);
    }
    string line, chr, prev;
    getline(fin, line);   //#phased ...
    getline(fin, line);   //column names
    vector< pair<string, unsigned int> > runs;
    set<string> seen;
    long lastStart = -1;
    unsigned int nrow = 0;
    while(getline(fin, line)){
        if(line.length() == 0) continue;
        size_t t1 = line.find('\t');
        if(t1 == string::npos){
            cerr << "ERROR: " << filename << " row " << nrow + 1 << " is not tab separated.\n";
            throw 0;
        }
        chr = line.substr(0, t1);
        size_t t2 = line.find('\t', t1 + 1);
        long start = atol(line.substr(t1 + 1, t2 - t1 - 1).c_str());
        if(nrow == 0 || chr.compare(prev) != 0){
            if(seen.count(chr) > 0){
                cerr << "ERROR: " << filename << " returns to contig " << chr
                     << " after leaving it. Rows must be grouped by contig.\n";
                throw 0;
            }
            seen.insert(chr);
            runs.push_back(make_pair(chr, 0u));
            prev = chr;
            lastStart = -1;
        }
        if(start <= lastStart){
            cerr << "ERROR: " << filename << " contig " << chr << " has a window starting at "
                 << start << " after one starting at " << lastStart
                 << ". Rows must ascend within a contig.\n";
            throw 0;
        }
        lastStart = start;
        runs.back().second++;
        nrow++;
    }
    fin.close();
    if(nrow != nwins){
        cerr << "ERROR: " << filename << " header says " << nwins << " windows but the file has "
             << nrow << " rows.\n";
        throw 0;
    }
    return runs;
}

} //namespace

//Read the spectra files into one block per CONTIG per population.
//
//Contigs are delimited by the chr column, not by the file, so one file may hold
//any number of contigs and N files may hold one each. Grouping this way is what
//makes the contig boundary structural: calcMTA bounds its flanking scan by
//SpectrumData::nwins, and for --dist-type nw the distance IS the window's index
//within the block, so a block spanning two contigs would draw neighbouring
//contig windows into the likelihood with no coordinate discontinuity to stop it.
map<string, vector<SpectrumData *>* > *readSpecData(vector<string> filenames){

    map<string, vector<SpectrumData *>* > *specDataByPopByChr = new map<string, vector<SpectrumData *>* >;
    set<string> contigsSeen;
    int K = 0, npop = 0;
    bool PHASED = false, HAPSTATS = false;
    vector<string> popNames;

    for (unsigned int f = 0; f < filenames.size(); f++){
        igzstream fin;
        fin.open(filenames[f].c_str());
        if(fin.fail()){
            cerr << "ERROR: Failed to open " << filenames[f] << " for reading.\n";
            throw IOError(filenames[f], false);
        }
        string line;
        getline(fin, line);
        SpecHeader h = parseSpecHeader(line, filenames[f]);
        getline(fin, line);   //column names

        vector< pair<string, unsigned int> > runs = scanContigRuns(filenames[f], h.nwins);

        if(h.haveContigs){
            bool ok = (h.contigNames.size() == runs.size());
            for (size_t i = 0; ok && i < runs.size(); i++)
                ok = (h.contigNames[i].compare(runs[i].first) == 0 && h.contigWins[i] == runs[i].second);
            if(!ok){
                cerr << "ERROR: " << filenames[f] << " header lists a contig layout its rows do not match.\n";
                throw 0;
            }
        }

        if(f == 0){
            K = h.K; npop = h.npop; PHASED = h.PHASED; HAPSTATS = h.HAPSTATS; popNames = h.popNames;
            for (int p = 0; p < npop; p++)
                specDataByPopByChr->operator[](popNames[p]) = new vector<SpectrumData *>;
        }
        else{
            if(K != h.K || npop != h.npop || PHASED != h.PHASED || HAPSTATS != h.HAPSTATS){
                cerr << "ERROR: Spectra files don't match.\n";
                throw 0;
            }
            for (int p = 0; p < h.npop; p++){
                if(specDataByPopByChr->count(h.popNames[p]) == 0){
                    cerr << "ERROR: Not all files have the same set of populations.\n";
                    throw 0;
                }
            }
        }

        cerr << "Loading " << filenames[f] << " with " << npop << " pops, K = " << K << " and "
             << runs.size() << (runs.size() == 1 ? " contig" : " contigs") << endl;

        for (size_t c = 0; c < runs.size(); c++){
            if(contigsSeen.count(runs[c].first) > 0){
                cerr << "ERROR: contig " << runs[c].first << " appears more than once across the spectra files.\n";
                throw 0;
            }
            contigsSeen.insert(runs[c].first);

            unsigned int n = runs[c].second;
            //info and dist describe the windows, not any one population, so one
            //copy is shared by every population's block for this contig. Any
            //cleanup must free them once per contig, NOT once per population.
            string **info = new string*[n];
            double *dist = new double[n];
            vector<SpectrumData *> blocks;
            for (int p = 0; p < npop; p++){
                SpectrumData *data = initSpecData(n, K, false, HAPSTATS);
                data->info = info;
                data->dist = dist;
                data->HAPSTATS = HAPSTATS;
                data->PHASED = PHASED;
                data->chr = runs[c].first;
                blocks.push_back(data);
                specDataByPopByChr->at(popNames[p])->push_back(data);
            }

            for (unsigned int w = 0; w < n; w++){
                info[w] = new string[4];
                if(!getline(fin, line)){
                    cerr << "ERROR: " << filenames[f] << " ended before its rows were read.\n";
                    throw 0;
                }
                stringstream ss(line);
                for (int i = 0; i < 4; i++) ss >> info[w][i];
                ss >> dist[w];
                for (int p = 0; p < npop; p++){
                    SpectrumData *data = blocks[p];
                    ss >> data->nhaps[w];
                    ss >> data->uhaps[w];
                    if(HAPSTATS){
                        ss >> data->h12[w];
                        ss >> data->h2h1[w];
                    }
                    for (int i = 0; i < K; i++) ss >> data->freq[w][i];
                }
                if(ss.fail()){
                    cerr << "ERROR: " << filenames[f] << " has a malformed row at "
                         << info[w][0] << ":" << info[w][1] << ".\n";
                    throw 0;
                }
            }
        }
        fin.close();
    }
    return specDataByPopByChr;
}

SpectrumData *initSpecData(int nwins, int K, bool doinfo, bool HAPSTATS){
    SpectrumData *data = new SpectrumData;
    data->K = K;
    data->nwins = nwins;
    data->freq = new double*[nwins];
    if(doinfo) data->info = new string*[nwins];
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
            data->info[j] = new string[4];
            for (int i = 0; i < 4; i++) data->info[j][i] = '0';
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
        throw IOError(filename, false);
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


/*
*/


//reads in map data and also does basic checks on integrity of format
//returns a populated MapData structure if successful
//throws an exception otherwise
/*


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

void releaseHapDataByPop(map< string, HaplotypeData* > *hapDataByPop)
{
    if (hapDataByPop == NULL) return;
    set<MapData *> maps;
    map< string, HaplotypeData* >::iterator it;
    for(it = hapDataByPop->begin(); it != hapDataByPop->end(); it++){
        if(it->second == NULL) continue;
        //Every population may point at the same MapData; free it once.
        if(it->second->map != NULL && maps.count(it->second->map) > 0) it->second->map = NULL;
        else if(it->second->map != NULL) maps.insert(it->second->map);
        releaseHapData(it->second);
    }
    delete hapDataByPop;
    return;
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
*/
/*
*/
//not done
/*
*/


//Locate the byte holding a locus inside a block-list row.
static inline unsigned char *cellPtr(vector<unsigned char*> &blocks, int locus){
    int idx = locus >> 2;
    return &blocks[idx / 16384][idx % 16384];
}

//A VCF is read in a single pass and handed back one contig at a time, so a
//file holding several contigs costs no more memory than one holding a single
//contig. Everything derived from the header -- the column to row mapping and
//the per-population haplotype counts -- lives in the handle and is computed
//once; only the genotype blocks and the locus names are per contig.
VCFReader *openVCF(string filename, PopData *popData, bool PHASED){
    VCFReader *r = new VCFReader;
    r->filename = filename;
    r->exhausted = false;
    r->ncontigs = 0;

    cerr << "Opening " << filename << "...\n";
    r->fin.open(filename.c_str());
    if (r->fin.fail()){
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw IOError(filename, false);
    }

    const int numMapCols = 9;
    string junk;
    map<string,bool> checkInd;

    while(getline(r->fin, junk)) if(junk[0] == '#' && junk[1] == 'C') break;

    r->nfields = (countFields(junk) - numMapCols);
    r->inds = new string[r->nfields];
    stringstream ss;
    ss.str(junk);
    for (int i = 0; i < numMapCols; i++) ss >> junk;
    for (int i = 0; i < r->nfields; i++){
        ss >> r->inds[i];
        checkInd[r->inds[i]] = true;
    }

    for (unsigned int i = 0; i < popData->indOrder.size(); i++){
        if(checkInd.count(popData->indOrder[i]) == 0){
            cerr << "ERROR: " << popData->indOrder[i] << " does not exist in " << filename << endl;
            throw 0;
        }
    }

    int nhaps = r->nfields;
    if(PHASED) nhaps *= 2;

    int nload = 0;
    for (int i = 0; i < popData->npops; i++) nload += (popData->pop2inds[popData->popOrder[i]].size());
    if(PHASED) nload *= 2;
    cerr << "Loading " << nload << "/" << nhaps << " ";
    if(PHASED) cerr << "phased";
    else if (!PHASED) cerr << "unphased";
    cerr << " haplotypes across " << popData->npops << " pops.\n";

    map<string,int> pop2indIndex;
    for (int i = 0; i < popData->npops; i++){
        string popName = popData->popOrder[i];
        r->pop2nhaps[popName] = (popData->pop2inds[popName].size()) * (PHASED ? 2 : 1);
        pop2indIndex[popName] = 0;
    }

    //Resolve each VCF sample column to the rows it writes, once. The genotype
    //loop used to look up ind2pop (twice), pop2indIndex and dataByPop for every
    //genotype of every locus -- five red-black-tree lookups keyed on a string,
    //all of them determined by the column index alone.
    r->row1 = new int[r->nfields];
    r->row2 = new int[r->nfields];
    r->nrows = 0;
    for (int field = 0; field < r->nfields; field++){
        r->row1[field] = -1;
        r->row2[field] = -1;
        if (popData->ind2pop.count(r->inds[field]) == 0) continue;
        string p = popData->ind2pop[r->inds[field]];
        int f = pop2indIndex[p]++;
        r->row1[field] = r->nrows++;
        r->rowPop.push_back(p);
        r->rowIndexInPop.push_back(PHASED ? 2*f : f);
        if (PHASED){
            r->row2[field] = r->nrows++;
            r->rowPop.push_back(p);
            r->rowIndexInPop.push_back(2*f + 1);
        }
    }

    return r;
}

void closeVCF(VCFReader *r){
    if (r == NULL) return;
    r->fin.close();
    delete [] r->inds;
    delete [] r->row1;
    delete [] r->row2;
    delete r;
    return;
}

//Read the next contig's records, or return NULL when the file is spent. The
//record that ends a contig is the first record of the next one, so it is held
//in the handle rather than re-read.
map< string, HaplotypeData* > *readContigVCF(VCFReader *r, PopData *popData, bool PHASED, bool SHARED_MAP){
    if (r->exhausted) return NULL;

    //Bytes per storage block while reading; a power of two so cellPtr divides
    //with shifts. 16 KB holds 65,536 loci, so slack is at most 16 KB per row.
    const int GT_BLOCK = 16384;

    //Rows grow as fixed blocks rather than as std::vector: a vector doubling to
    //hold n bytes peaks at 2n, and with one row per haplotype that overshoot is
    //the whole genotype matrix again. Blocks are never copied or reallocated.
    vector< vector<unsigned char*> > rows(r->nrows);

    //deque, not vector: a vector doubling to hold one entry per locus copies
    //and briefly holds two arrays of every locus name
    deque<string> locusNames;
    deque<unsigned int> physicalPos;
    string contig;

    string line, chr, name;
    unsigned int pos;
    int nloci = 0;

    while (true)
    {
        if (!r->pending.empty()){
            line.swap(r->pending);
            r->pending.clear();
        }
        else if (!getline(r->fin, line)){
            r->exhausted = true;
            break;
        }

        if (line.size() == 0 || line[0] == '#') continue;

        const char *c = line.c_str();
        const char *lineEnd = c + line.size();

        //CHROM. Leading whitespace is skipped before the first token, not only
        //between tokens. Reading the token straight from line[0] ended it
        //immediately on a record beginning with a space or tab, giving an empty
        //contig name and shifting every later field by one: the FORMAT column
        //was then read as the first sample's genotype and the run died on
        //"Alleles must be coded 0/1/. only". lassip <= 1.2.2 took the fixed
        //columns by stream extraction (fin >> chr >> pos >> ...), which skips
        //whitespace, so the fields were not shifted there; that is about field
        //extraction alone and not a claim that such a file completed a run.
        while (c < lineEnd && (*c == '\t' || *c == ' ')) c++;
        const char *tok = c;
        while (c < lineEnd && *c != '\t' && *c != ' ') c++;
        chr.assign(tok, c - tok);
        if (chr.empty()){
            cerr << "ERROR: " << r->filename << " has a record with no CHROM field.\n";
            throw 0;
        }
        while (c < lineEnd && (*c == '\t' || *c == ' ')) c++;
        //POS
        pos = 0;
        while (c < lineEnd && *c >= '0' && *c <= '9'){ pos = pos*10 + (*c - '0'); c++; }
        while (c < lineEnd && (*c == '\t' || *c == ' ')) c++;
        //ID
        tok = c;
        while (c < lineEnd && *c != '\t' && *c != ' ') c++;
        name.assign(tok, c - tok);
        //REF ALT QUAL FILTER INFO FORMAT
        for (int skip = 0; skip < 6; skip++){
            while (c < lineEnd && (*c == '\t' || *c == ' ')) c++;
            while (c < lineEnd && *c != '\t' && *c != ' ') c++;
        }

        if (nloci == 0){
            //Contigs are returned one at a time in the order they appear, so a
            //contig whose records are split into separate runs would have to be
            //buffered to be assembled. Refuse instead: sorting the file is the
            //user's one-line fix, and silently treating the runs as separate
            //contigs would put two blocks under one name in the spectra.
            if (r->seenContigs.count(chr) > 0){
                cerr << "ERROR: the records for contig " << chr << " in " << r->filename
                     << " are not all together. lassip reads a VCF in one pass, so its\n"
                     << "\trecords must be grouped by contig (sorting by position does this).\n";
                throw 0;
            }
            r->seenContigs.insert(chr);
            contig = chr;
        }
        else if (chr != contig){
            //first record of the next contig: keep it for the next call
            r->pending = line;
            break;
        }

        int locus = nloci;
        //every fourth locus opens a new byte, and every GT_BLOCK bytes a new block
        if ((locus & 3) == 0 && (((locus >> 2) & (GT_BLOCK - 1)) == 0)){
            for (int rw = 0; rw < r->nrows; rw++){
                rows[rw].push_back(new unsigned char[GT_BLOCK]);
                memset(rows[rw].back(), 0xFF, GT_BLOCK);
            }
        }

        locusNames.push_back(name);
        physicalPos.push_back(pos);

        for (int field = 0; field < r->nfields; field++)
        {
            while (c < lineEnd && (*c == '\t' || *c == ' ')) c++;
            const char *gt = c;
            while (c < lineEnd && *c != '\t' && *c != ' ') c++;
            size_t gtlen = c - gt;
            if (gtlen == 0){
                cerr << "ERROR: " << r->filename << " has fewer genotype fields than samples at "
                     << chr << ":" << pos << ".\n";
                throw 0;
            }
            if (r->row1[field] < 0) continue;

            char allele1, allele2;
            if (gtlen == 1 && gt[0] == VCF_MISSING){
                allele1 = VCF_MISSING;
                allele2 = VCF_MISSING;
            }
            else{
                allele1 = gt[0];
                //a haploid or truncated GT has no second allele; treat it as missing
                //rather than reading past the end of the field
                allele2 = (gtlen > 2) ? gt[2] : VCF_MISSING;
            }

            if((allele1 != '1' && allele1 != '0' && allele1 != VCF_MISSING) ||
               (allele2 != '1' && allele2 != '0' && allele2 != VCF_MISSING)){
                //name the offending field: when a malformed record shifts the
                //columns, this error is where it surfaces, and the bare message
                //gave nothing to look at
                cerr << "ERROR: Alleles must be coded 0/1/. only. Read \""
                     << string(gt, gtlen) << "\" for sample " << r->inds[field]
                     << " at " << chr << ":" << pos << " in " << r->filename << ".\n";
                throw 0;
            }

            if(PHASED){
                setGTInByte(cellPtr(rows[r->row1[field]], locus), locus, (allele1 == VCF_MISSING) ? GT_MISS : gtCode(allele1));
                setGTInByte(cellPtr(rows[r->row2[field]], locus), locus, (allele2 == VCF_MISSING) ? GT_MISS : gtCode(allele2));
            }
            else{
                unsigned char code;
                if (allele1 == VCF_MISSING || allele2 == VCF_MISSING) code = GT_MISS;
                else if (allele1 == '1' && allele2 == '1') code = GT_2;
                else if (allele1 == '0' && allele2 == '0') code = GT_0;
                else code = GT_1;
                setGTInByte(cellPtr(rows[r->row1[field]], locus), locus, code);
            }
        }

        nloci++;
    }

    if (nloci < 1){
        if (r->ncontigs == 0){
            cerr << "ERROR: " << r->filename << " contains no variant records.\n";
            throw 0;
        }
        return NULL;
    }
    r->ncontigs++;
    cerr << "Read " << nloci << " loci from " << contig << ".\n";

    MapData *mapData = NULL;
    if(SHARED_MAP) mapData = initMapData(nloci);

    map<string, HaplotypeData* > *dataByPop = new map<string, HaplotypeData* >;
    for (int i = 0; i < popData->npops; i++){
        string popName = popData->popOrder[i];
        HaplotypeData *hd = initHaplotypeData(r->pop2nhaps[popName], nloci, !SHARED_MAP, false);
        if(SHARED_MAP) hd->map = mapData;
        dataByPop->operator[](popName) = hd;
    }

    //Hand each grown row to its population and free it immediately, so the
    //growable copy and the final matrix never both hold the whole dataset.
    int stride = gtStride(nloci);
    for (int rw = 0; rw < r->nrows; rw++){
        HaplotypeData *hd = dataByPop->at(r->rowPop[rw]);
        unsigned char *dst = new unsigned char[stride + 1];
        memset(dst, 0xFF, stride + 1);
        hd->data[r->rowIndexInPop[rw]] = dst;
        for (unsigned int b = 0; b < rows[rw].size(); b++){
            int off = b * GT_BLOCK;
            int n = (stride - off < GT_BLOCK) ? stride - off : GT_BLOCK;
            if (n > 0) memcpy(dst + off, rows[rw][b], n);
            delete [] rows[rw][b];
        }
        dst[stride] = 0xFF;
        vector<unsigned char*>().swap(rows[rw]);
    }

    for (int i = 0; i < popData->npops; i++){
        MapData *md = dataByPop->at(popData->popOrder[i])->map;
        if (md->nloci != nloci) continue;
        md->chr = contig;
        for (int l = 0; l < nloci; l++){
            //the last map to be filled can take the names rather than copy them
            if (SHARED_MAP || i + 1 == popData->npops) md->locusName[l].swap(locusNames[l]);
            else md->locusName[l] = locusNames[l];
            md->physicalPos[l] = physicalPos[l];
        }
        if (SHARED_MAP) break;
    }

    return dataByPop;
}


/*
*/
/*
*/
HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci, bool domap)
{
    return initHaplotypeData(nhaps, nloci, domap, true);
}

HaplotypeData *initHaplotypeData(unsigned int nhaps, unsigned int nloci, bool domap, bool allocRows)
{
    if (nhaps < 1 || nloci < 1)
    {
        cerr << "ERROR: number of haplotypes (" << nhaps << ") and number of loci (" << nloci << ") must be positive.\n";
        throw 0;
    }

    HaplotypeData *data = new HaplotypeData;
    data->nhaps = nhaps;
    data->nloci = nloci;

    data->stride = gtStride(nloci);
    data->data = new unsigned char *[nhaps];
    for (unsigned int i = 0; i < nhaps; i++) data->data[i] = NULL;
    //allocRows false lets the caller fill rows one at a time and free its own
    //storage as it goes, so the two copies never coexist
    if (allocRows){
        for (unsigned int i = 0; i < nhaps; i++)
        {
            //one padding byte so extractWindow can read one byte past the last
            data->data[i] = new unsigned char[data->stride + 1];
            for (int j = 0; j <= data->stride; j++) data->data[i][j] = 0xFF;  //all missing
        }
    }

    if (domap) data->map = initMapData(nloci);
    //data->Q = NULL;

    return data;
}

void releaseHapData(HaplotypeData *data)
{
    if (data == NULL) return;
    for (int i = 0; i < data->nhaps; i++)
    {
        if (data->data[i] != NULL) delete [] data->data[i];
    }

    delete [] data->data;

    if(data->map != NULL) releaseMapData(data->map);
    //if(data->freq != NULL) releaseFreqData(data->freq);
    data->map = NULL;
    data->data = NULL;
    data->nhaps = -9;
    data->nloci = -9;
    data->stride = -9;
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
*/

bool GMapData::getMapInfo(double queryPos, double &gPos, string &locName, string c, int &current_index)
{
    bool success = true;

    //map<>::operator[] would insert a null entry for a contig the map does not
    //cover, and the branches below would then dereference it
    if (nloci.count(c) == 0 || nloci[c] < 2) return false;

    if (queryPos < physicalPos[c][0])
    {
        gPos = this->interpolate(physicalPos[c][0], geneticPos[c][0],
                                 physicalPos[c][1], geneticPos[c][1],
                                 queryPos);
        success = true;
    }
    else if (queryPos > physicalPos[c][nloci[c] - 1])
    {
        gPos = this->interpolate(physicalPos[c][nloci[c] - 2], geneticPos[c][nloci[c] - 2],
                                 physicalPos[c][nloci[c] - 1], geneticPos[c][nloci[c] - 1],
                                 queryPos);
        char buffer[50];
        snprintf(buffer, sizeof(buffer), "chr%s:%f", c.c_str(), queryPos);
        locName = buffer;
    }
    else if (ppos2index[c].count(queryPos) > 0)
    {
        int index = ppos2index[c][queryPos];
        gPos = geneticPos[c][index];
        locName = locusName[c][index];
    }
    else
    {
        //Map positions are sorted, so find the bracketing pair by binary search.
        //This used to scan forward from current_index, and if the scan found
        //nothing -- which happens whenever a query arrives below where a
        //previous one left the cursor -- it fell through with startIndex and
        //endIndex uninitialised and interpolated between two arbitrary indices
        //of physicalPos.
        int lo = 0;
        int hi = nloci[c] - 1;
        while (hi - lo > 1)
        {
            int mid = lo + (hi - lo) / 2;
            if (physicalPos[c][mid] <= queryPos) lo = mid;
            else hi = mid;
        }

        if (!(queryPos > physicalPos[c][lo] && queryPos < physicalPos[c][hi])) return false;

        int startIndex = lo;
        int endIndex = hi;
        current_index = startIndex;

        //MAXGAP <= 0 means the user asked for no limit (--max-gap 0)
        if (MAXGAP > 0 && physicalPos[c][endIndex] - physicalPos[c][startIndex] > MAXGAP) return false;

        gPos = this->interpolate(physicalPos[c][startIndex], geneticPos[c][startIndex],
                                 physicalPos[c][endIndex], geneticPos[c][endIndex],
                                 queryPos);
        char buffer[50];
        snprintf(buffer, sizeof(buffer), "chr%s:%f", c.c_str(), queryPos);
        locName = buffer;
    }
    return success;
}

GMapData::GMapData(string filename, double mGap)
{
    ifstream fin;
    cerr << "Opening " << filename << "...\n";
    fin.open(filename.c_str());

    if (fin.fail())
    {
        cerr << "ERROR: Failed to open " << filename << " for reading.\n";
        throw IOError(filename, false);
    }


    map<string,int> nloci;
    vector<string> chrstr;
    int fileStart = fin.tellg();
    string line;
    int n = 0;
    int num_cols = 4;
    int current_cols = 0;
    string currChr = "-1";
    stringstream ss;

    while (getline(fin, line))
    {
        ss.str(line);
        ss >> currChr;

        if(nloci.count(currChr) == 0){
            nloci[currChr] = 0;
            chrstr.push_back(currChr);
        }
        nloci[currChr]++;

        n++;
        current_cols = countFields(line);
        if (current_cols != num_cols)
        {
            cerr << "ERROR: line " << n << " of " << filename << " has " << current_cols
                 << ", but expected " << num_cols << ".\n";
            throw 0;
        }
    }

    fin.clear();
    fin.seekg(fileStart);

    cerr << "Loading map data for " << n << " loci.\n";

    //this->initGMapData(n+1);
    this->initGMapData(chrstr,nloci,false);

    string c, chrCheck;

    for(long unsigned int i = 0; i != chrstr.size(); i++){
        c = chrstr[i];
        for (int locus = 0; locus < nloci[c]; locus++)
        {
            /*
            if(locus == 0){
                locusName[c][0] = "0";
                geneticPos[c][0] = 0;
                physicalPos[c][0] = 0;
                ppos2index[c][physicalPos[c][0]] = 0;
                cout << c << " " << physicalPos[c][locus] << " " << geneticPos[c][locus] << endl;
                continue;
            }
            */
            fin >> chrCheck;
            if(c != chrCheck){
                cerr << "ERROR: " << c << " " << chrCheck << ": Mismatch among chromosomes when reading genetic map.\n";
                throw 0;
            }
            fin >> locusName[c][locus];
            fin >> geneticPos[c][locus];
            fin >> physicalPos[c][locus];
            //cout << c << " " << physicalPos[c][locus] << " " << geneticPos[c][locus] << endl;
            ppos2index[c][physicalPos[c][locus]] = locus;
        }
    }

    MAXGAP = mGap;
    fin.close();
    return;
}

void GMapData::initGMapData(vector<string> chrstr, map<string,int> nl, bool ADD){

    for(long unsigned int i = 0; i != chrstr.size(); i++){
        if(ADD) nl[chrstr[i]]++;
    }

    chr = chrstr;
    nloci = nl;

    int n;
    string c;

    for(long unsigned int i = 0; i != chrstr.size(); i++){
        c = chrstr[i];
        n = nloci[c];
        physicalPos[c] = new int[n];
        geneticPos[c] = new double[n];
        locusName[c] = new string[n];

        for(int l = 0; l < n; l++){
            physicalPos[c][l] = MISSING;
            geneticPos[c][l] = MISSING;
            locusName[c][l] = "--";
        }
    }
    return;
}
GMapData::~GMapData()
{
    for(long unsigned int i = 0; i != chr.size(); i++){
        string c = chr[i];
        delete [] locusName[c];
        delete [] physicalPos[c];
        delete [] geneticPos[c];
    }
}

int GMapData::countFields(const string &str)
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

void fillNWDistance(map<string, vector<SpectrumData *>* > *specDataByPopByChr){
    map<string, vector<SpectrumData *>* >::iterator it;
    for(it = specDataByPopByChr->begin(); it != specDataByPopByChr->end(); it++){
        vector<SpectrumData *> *specDataByChr = it->second;
        for(unsigned int i = 0; i < specDataByChr->size(); i++){
            for (int w = 0; w < specDataByChr->at(i)->nwins; w++){
                specDataByChr->at(i)->dist[w] = w;
            }
        }
    }
    return;
}

void fillCMDistance(map<string, vector<SpectrumData *>* > *specDataByPopByChr, GMapData &geneticMap){
    string c;
    string locName;
    double gPos;
    int current_locus;

    //getMapInfo returns false when it cannot place a window: the contig is not
    //in the map, or the window falls in a gap wider than MAXGAP. Its gPos is
    //then untouched, and this loop used to write it anyway -- the previous
    //window's genetic position, or an uninitialised double for the first.
    long unplaced = 0;
    string firstBadChr;
    double firstBadPos = 0;

    map<string, vector<SpectrumData *>* >::iterator it;
    for(it = specDataByPopByChr->begin(); it != specDataByPopByChr->end(); it++){
        vector<SpectrumData *> *specDataByChr = it->second;
        for(unsigned int i = 0; i < specDataByChr->size(); i++){
            current_locus = 0;
            for (int w = 0; w < specDataByChr->at(i)->nwins; w++){
                c = specDataByChr->at(i)->info[w][0];
                gPos = 0;
                if(!geneticMap.getMapInfo(specDataByChr->at(i)->dist[w],gPos,locName,c,current_locus)){
                    if(unplaced == 0){
                        firstBadChr = c;
                        firstBadPos = specDataByChr->at(i)->dist[w];
                    }
                    unplaced++;
                    continue;
                }
                specDataByChr->at(i)->dist[w] = gPos;
            }
        }
    }

    if(unplaced > 0){
        cerr << "ERROR: the genetic map does not place " << unplaced << " window(s), the first at "
             << firstBadChr << ":" << firstBadPos << ".\n";
        cerr << "\tThe map must cover every contig in the spectra, and gaps wider than "
             << geneticMap.maxGap() << " bp are not interpolated across (--max-gap).\n";
        throw 0;
    }

    return;
}
