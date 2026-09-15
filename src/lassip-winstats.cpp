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
#include "lassip-winstats.h"
#include "lassip-cli.h"  //--hap-cluster values
#include <algorithm> //shuffle
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

double getDMin(vector<SpectrumData *> *specDataByChr){
   double dmin = 9999999999;

   for(unsigned int c = 0; c < specDataByChr->size(); c++){
      for(int w = 1; w < specDataByChr->at(c)->nwins; w++){
         double diff = fabs(specDataByChr->at(c)->dist[w]-specDataByChr->at(c)->dist[w-1]);
         if(dmin > diff && diff > 0) dmin = diff;
      }
   }
   return dmin;
}

//Number of points on the epsilon grid. initQ, releaseQ, calcQ and calcMTA must
//all agree on this: the array is allocated to this size and indexed by the loops
//below, and the previous code derived it two different ways -- int(U*100*K) for
//the allocation, and by accumulating e += epsStep for the loops -- so rounding
//could in principle have given the loops one iteration more than the allocation.
int nEpsGrid(int K, double U){
   int n = int(U * 100.0 * double(K));
   return (n > 0) ? n : 0;
}

//Allocate the sweep-spectrum table. Its values depend only on the null spectrum,
//the scaling choice, m and epsilon -- not on the window -- so one table serves
//every window of a contig. It used to be allocated per window, which for the
//23,208-window YRI chr22 example meant 852 MB of identical copies (of a 1.0 GB
//peak RSS) and a cache miss on every window; the distinct data is 37 KB.
double ***initQ(int K, double U){
   int nEps = nEpsGrid(K, U);
   double ***q = new double **[nEps];
   for(int e = 0; e < nEps; e++){
      q[e] = new double *[K-1];
      for(int m = 0; m < K-1; m++) q[e][m] = new double[K];
   }
   return q;
}

void releaseQ(double ***q, int K, double U){
   int nEps = nEpsGrid(K, U);
   for(int e = 0; e < nEps; e++){
      for(int m = 0; m < K-1; m++) delete [] q[e][m];
      delete [] q[e];
   }
   delete [] q;
   return;
}

//The epsilon grid point at index ei (0-based), i.e. the ei+1'th step of 1/(100K).
double epsAt(int ei, int K){
   return double(ei + 1) / (100.0 * double(K));
}

void calcQ(double ***q, SpectrumData *avgSpec, double **f){
   int K = avgSpec->K;
   double U = avgSpec->freq[0][K-1];
   int nEps = nEpsGrid(K, U);
   for (int ei = 0; ei < nEps; ei++){
      for (int m = 1; m < K; m++){
         calcQ(q[ei][m-1], avgSpec, f, U, m, epsAt(ei, K));
      }
   }
   return;
}
void calcQ(double *q, SpectrumData *avgSpec, double **f, double U, int m, double e){
   int K = avgSpec->K;

   for(int i = m+1; i <= K; i++){
      if(m == K-1) q[i-1] = e;
      else q[i-1] = U - double(i-m-1.0)/double(K-m-1.0) * (U - e);
   }
   for(int i = 1; i <= m; i++){
      q[i-1] = 0;
      for(int j = m+1; j <= K; j++) q[i-1] += avgSpec->freq[0][j-1] - q[j-1];
      q[i-1] *= f[m-1][i-1];
      q[i-1] += avgSpec->freq[0][i-1];
   }

   return;
}
double calcH12(HaplotypeFrequencySpectrum *hfs, bool PHASED){
   double tot = hfs->size;
   const double *c = hfs->sortedCount;
   if(PHASED){
      if(hfs->numClasses == 1){
         return (double(c[0])/tot)*(double(c[0])/tot);
      }
      else if(hfs->numClasses == 2){
         return ((double(c[0])/tot)+(double(c[1])/tot))*((double(c[0])/tot)+(double(c[1])/tot));
      }
      else if(hfs->numClasses > 2){
         double res = ((double(c[0])/tot)+(double(c[1])/tot))*((double(c[0])/tot)+(double(c[1])/tot));
         for(int i = 2; i < hfs->numClasses; i++){
            res += (double(c[i])/tot)*(double(c[i])/tot);
         }
         return res;
      }
   }
   else{
      if(hfs->numClasses == 1){
         return (double(c[0])/tot)*(double(c[0])/tot);
      }
      else if(hfs->numClasses == 2){
         return ((double(c[0])/tot)+(double(c[1])/tot))*((double(c[0])/tot)+(double(c[1])/tot));
      }
      else if(hfs->numClasses == 3){
         return ((double(c[0])/tot)+(double(c[1])/tot)+(double(c[2])/tot))*((double(c[0])/tot)+(double(c[1])/tot)+(double(c[2])/tot));
      }
      else if(hfs->numClasses > 3){
         double res = ((double(c[0])/tot)+(double(c[1])/tot)+(double(c[2])/tot))*((double(c[0])/tot)+(double(c[1])/tot)+(double(c[2])/tot));
         for(int i = 3; i < hfs->numClasses; i++){
            res += (c[i]/tot)*(c[i]/tot);
         }
         return res;
      }
   }
   
   return -1;
}

double calcH2H1(HaplotypeFrequencySpectrum *hfs){
   const double *c = hfs->sortedCount;
   double tot = hfs->size;
   double first = (double(c[0])/tot)*(double(c[0])/tot);
   double res = first;
   for(int i = 1; i < hfs->numClasses; i++){
      res += (double(c[i])/tot)*(double(c[i])/tot);
   }
   return (res-first)/res;
}


void calcMTA(LASSIResults *results, double ***q, SpectrumData *specData, SpectrumData *avgSpec, int w, double dmin,double MAX_EXTEND){
   int rightLim, leftLim;
   double *dist = specData->dist;
   int d = w;
   while(fabs(dist[w] - dist[d]) <= MAX_EXTEND){
      d++;
      if(d >= specData->nwins){
         d--;
         break;
      }
   }
   leftLim = d;
   d = w;
   while(fabs(dist[w] - dist[d]) <= MAX_EXTEND){
      d--;
      if(d < 0){
         d++;
         break;
      }
   }
   rightLim = d;

   int K = avgSpec->K;
   double U = avgSpec->freq[0][K-1];
   int nEps = nEpsGrid(K, U);
   int nloc = leftLim - rightLim + 1;

   //The likelihood being maximised is
   //   L(A,e,m) = sum_win [ Pr(A,win)*sum_i n_win f_win,i log q_e,m,i
   //                        + (1-Pr(A,win))*sum_i n_win f_win,i log p_i ]
   //Neither logarithm depends on A, and q does not depend on the window, so the
   //grid search below only needs the per-window sums, not the logs themselves.
   //Computing them inside the (A,m,e) loops, as this function used to, evaluates
   //log() ~10^8 times per window where a few thousand distinct values exist.

   //log of the null spectrum
   double *logp = new double[K];
   for (int i = 0; i < K; i++) logp[i] = log(avgSpec->freq[0][i]);

   //null contribution of each window in the flanking region
   double *P = new double[nloc];
   double nullLikelihood = 0;
   for (int j = 0; j < nloc; j++){
      int win = rightLim + j;
      double n = double(specData->nhaps[win]);
      double s = 0;
      for (int i = 0; i < K; i++) s += n*specData->freq[win][i]*logp[i];
      P[j] = s;
      nullLikelihood += s;
   }

   //sweep contribution of each window, for every (eps, m), minus the null one
   double *Q = new double[(size_t)nEps*(K-1)*nloc];
   double *lq = new double[K];
   for (int e = 0; e < nEps; e++){
      for (int m = 0; m < K-1; m++){
         for (int i = 0; i < K; i++) lq[i] = log(q[e][m][i]);
         double *Qem = Q + ((size_t)e*(K-1) + m)*nloc;
         for (int j = 0; j < nloc; j++){
            int win = rightLim + j;
            double n = double(specData->nhaps[win]);
            double s = 0;
            for (int i = 0; i < K; i++) s += n*specData->freq[win][i]*lq[i];
            Qem[j] = s - P[j];
         }
      }
   }

   double Amin = -log(0.99999)/dmin;
   double Amax = -log(0.00001)/dmin;
   double lAmin = log(Amin);
   double lAmax = log(Amax);
   double lstep = (lAmax-lAmin)/100;

   int maxM = -1;
   double maxA = -1;
   double maxAltLikelihood = -99999999;
   double *Pr = new double[nloc];

   //NOTE: where the likelihood ridge in A is flat -- which it is for a few percent
   //of windows -- which of two adjacent grid points attains the maximum is decided
   //by rounding, so the reported A moves under any change to summation order,
   //compiler or optimisation level. The maximised likelihood itself, and with it
   //m and T, are unaffected. Reporting the interval of A within some delta of the
   //optimum, or an explicit tie-breaking rule, would remove the arbitrariness.
   for (double A = lAmin; A <= lAmax; A += lstep){
      double expA = exp(A);
      for (int j = 0; j < nloc; j++) Pr[j] = exp(-expA*fabs(dist[w]-dist[rightLim+j]));
      for (int m = 1; m <= K; m++){
         if(m == K){
            //a sweep involving all K classes is the neutral background
            if(nullLikelihood > maxAltLikelihood){
               maxAltLikelihood = nullLikelihood;
               maxM = m;
               maxA = 1.0/expA;
            }
            continue;
         }
         for (int e = 0; e < nEps; e++){
            double *Qem = Q + ((size_t)e*(K-1) + (m-1))*nloc;
            double alt = nullLikelihood;
            for (int j = 0; j < nloc; j++) alt += Pr[j]*Qem[j];
            if(alt > maxAltLikelihood){
               maxAltLikelihood = alt;
               maxM = m;
               maxA = 1.0/expA;
            }
         }
      }
   }

   //m == K is identical to the neutral background
   //so set the number of sweeping haplotypes to 0
   if(maxM == K){
      maxM = 0;
      maxA = 0;
   }
   results->A[w] = maxA;
   results->m[w] = maxM;
   results->T[w] = 2.0 * (maxAltLikelihood - nullLikelihood);

   delete [] logp;
   delete [] P;
   delete [] Q;
   delete [] lq;
   delete [] Pr;
   return;
}

void calcMandT(LASSIResults *results, SpectrumData *specData, SpectrumData *avgSpec, double **f, int w){
   double nullLikelihood = calcLASSINullLikelihood(specData,avgSpec,w);
   //cerr << "null: " << nullLikelihood << endl;
   int K = avgSpec->K;
   double U = avgSpec->freq[0][K-1];

   int maxM = -1;
   //double maxE = -1;
   double maxAltLikelihood = -99999999;
   double altLikelihood = -99999999;
   double epsStep = 1.0/(100.0*double(K));
   for (int m = 1; m <= K; m++){
      for (double e = epsStep; e <= U; e += epsStep){
         altLikelihood = calcLASSIAltLikelihood(specData, avgSpec, f, U, m, e, w);
         if(altLikelihood > maxAltLikelihood){
            maxAltLikelihood = altLikelihood;
            maxM = m;
            //maxE = e;
         }
      }
   }

   //m == K is identical to the neutral background
   //so set the number of sweeping haplotypes to 0
   if(maxM == K) maxM = 0;

   results->m[w] = maxM;
   results->T[w] = 2.0 * (maxAltLikelihood - nullLikelihood);
   return;
}

double calcLASSINullLikelihood(SpectrumData *specData,SpectrumData *avgSpec,int w){
   double res = 0;
   for (int i = 0; i < avgSpec->K; i++){
      res += double(specData->nhaps[w])*specData->freq[w][i]*log(avgSpec->freq[0][i]);
   }
   return res;
}

double calcLASSIAltLikelihood(SpectrumData *specData, SpectrumData *avgSpec, double **f, double U, int m, double e, int w){
   int K = avgSpec->K;
   if(m == K) return calcLASSINullLikelihood(specData,avgSpec,w);

   double *q = new double[K];

   for(int i = m+1; i <= K; i++){
      if(m == K-1) q[i-1] = e;
      else q[i-1] = U - double(i-m-1.0)/double(K-m-1.0) * (U - e);
   }
   for(int i = 1; i <= m; i++){
      q[i-1] = 0;
      for(int j = m+1; j <= K; j++) q[i-1] += avgSpec->freq[0][j-1] - q[j-1];
      q[i-1] *= f[m-1][i-1];
      q[i-1] += avgSpec->freq[0][i-1];
   }

   double res = 0;
   for (int i = 0; i < avgSpec->K; i++){
      res += double(specData->nhaps[w])*specData->freq[w][i]*log(q[i]);
   }
   delete [] q;
   return res;
}


double **calcF(int type, int K){
   double **f = new double*[K];
   
   for(int i = 0; i < K; i++) f[i] = new double[i+1];
   
   double sum = 0;
   if(type == 1){
      for(int i = 0; i < K; i++){
         sum = 0;
         for(int j = 0; j < i+1; j++){
            f[i][j] = 1.0 / double(i+1);
            sum += f[i][j];
         }
         for(int j = 0; j < i+1; j++) f[i][j] /= sum;
      } 
   }
   else if (type == 2){
      for(int i = 0; i < K; i++){
         sum = 0;
         for(int j = 0; j < i+1; j++){
            f[i][j] = 1.0 / double(j+1);
            sum += f[i][j];
         }
         for(int j = 0; j < i+1; j++) f[i][j] /= sum;
      } 
   }
   else if (type == 3){
      for(int i = 0; i < K; i++){
         sum = 0;
         for(int j = 0; j < i+1; j++){
            f[i][j] = 1.0 / ( double(j+1) * double(j+1) );
            sum += f[i][j];
         }
         for(int j = 0; j < i+1; j++) f[i][j] /= sum;
      } 
   }
   else if (type == 4){
      for(int i = 0; i < K; i++){
         sum = 0;
         for(int j = 0; j < i+1; j++){
            f[i][j] = 1.0 / exp(double(j+1));
            sum += f[i][j];
         }
         for(int j = 0; j < i+1; j++) f[i][j] /= sum;
      } 
   }
   else if (type == 5){
      for(int i = 0; i < K; i++){
         sum = 0;
         for(int j = 0; j < i+1; j++){
            f[i][j] = 1.0 / exp(double(j+1) * double(j+1));
            sum += f[i][j];
         }
         for(int j = 0; j < i+1; j++) f[i][j] /= sum;
      } 
   }
   else{
      cerr << "ERROR: Invalid scaling choice.\n";
      throw 0;
   }
   return f;
}

//This function immediately stops counting differences once MATCH_TOL is exceeded
//It also combines loci where str1 is missing but str2 is not, stored in str3
//Intended usage is to use str3 to replace str1 iff ndiff == 0.
int garud_ndiff_str(const string &str1, const string &str2, string &str3, int MATCH_TOL){
   int ndiff = 0;

   str3 = str1;

   if(str1.length() != str2.length()){
      cerr << "WARNING: haplotypes not of same length!\n";
      return -1;
   }

   for (size_t i = 0; i < str1.length(); i++){
      if(str1[i] != str2[i]){
         if (str2[i] != MISSING_ALLELE){
            if (str1[i] != MISSING_ALLELE){
               ndiff++;
               if(ndiff > MATCH_TOL){
                  return ndiff;
               }
            }
            else{
               str3.replace(i,1,1,str2[i]);
            }
         }
      }
   }
   return ndiff;
}



//Derive a per-window seed. Mixing the window's SNP boundaries into the user seed
//makes the shuffle below depend only on the seed and the window itself, so results
//do not change with --threads or with the order in which windows are processed.
unsigned int windowSeed(int seed, int start, int end){
   if (seed == 0) return (unsigned int) chrono::system_clock::now().time_since_epoch().count();
   unsigned long long z = (unsigned long long)(unsigned int)seed  * 0x9E3779B97F4A7C15ULL
                        + (unsigned long long)(unsigned int)start * 0xBF58476D1CE4E5B9ULL
                        + (unsigned long long)(unsigned int)end   * 0x94D049BB133111EBULL;
   z ^= z >> 31;
   z *= 0xBF58476D1CE4E5B9ULL;
   z ^= z >> 29;
   return (unsigned int)(z ^ (z >> 32));
}

//The haplotype order fed to the clustering below changes which incomplete
//haplotype gets merged into which complete one, so the shuffle is a real
//degree of freedom in the result, not just a test harness. The seed is
//therefore a user-visible parameter (--seed); see windowSeed above.
void garud_match_haps_w_missing_shuffle(map<string,double> &hap2count,map<string,double> &miss_hap2count, int len, int MATCH_TOL, unsigned int seed){

   map<string, double>::iterator it1;
   map<string, double>::iterator it2;
   //map<string, int>::iterator it3;

   vector<string> hapIDs;
   //Combining them, this is a little hacky, as I originally planned to handle them differently
   map<string, double> hap2countCombined;
   for (it1 = hap2count.begin(); it1 != hap2count.end(); it1++){
      hap2countCombined[it1->first] = it1->second;
      hapIDs.push_back(it1->first);
   }
   for (it1 = miss_hap2count.begin(); it1 != miss_hap2count.end(); it1++){
      hap2countCombined[it1->first] = it1->second;
      hapIDs.push_back(it1->first);
   }

   shuffle(hapIDs.begin(),hapIDs.end(),default_random_engine(seed));

   map<string, int> compared;
   string hap1, hap2, mergedhap;
   double count1, count2;
   hap2count.clear();

   for (size_t i = 0; i < hapIDs.size(); i++){
      hap1 = hapIDs[i];
      count1 = hap2countCombined[hap1];
      if(compared.count(hap1) == 0){
         compared[hap1] = 1;
         hap2count[hap1] = count1;
      }
      for (size_t j = 0; j < hapIDs.size(); j++){
         hap2 = hapIDs[j];
         count2 = hap2countCombined[hap2];
         if(compared.count(hap2) == 0){
            int d = garud_ndiff_str(hap1,hap2,mergedhap,MATCH_TOL);
            if(d == 0 && mergedhap.compare(hap1) != 0){
               hap2count[mergedhap] = hap2count[hap1];
               hap2count.erase(hap1);
               hap1 = mergedhap;
            }
            //<= , not < : --match-tol is documented as "<= this many pairwise
            //differences". With < , --match-tol 0 pooled nothing and every
            //setting behaved as the one below it.
            if(d <= MATCH_TOL){
               hap2count[hap1] += count2;
               compared[hap2] = 1;
            }
         }
      }
   }

/*
   vector<string> to_delete;
   for (it1 = hap2count.begin(); it1 != hap2count.end(); it1++){
      if(it1->second < 2){
         to_delete.push_back(it1->first);
      }
   }

   for (int i = 0; i < to_delete.size(); i++){
      hap2count.erase(to_delete[i]);
   }
*/
   return;
}




int clusterMethodCode(const string &name){
   if (name.compare(HAP_CLUSTER_GARUD) == 0) return CLUSTER_GARUD_SHUFFLE;
   if (name.compare(HAP_CLUSTER_BESTCOMP) == 0) return CLUSTER_BEST_COMP;
   if (name.compare(HAP_CLUSTER_SOFTEM) == 0) return CLUSTER_SOFT_EM;
   return -1;
}

//Number of sites actually called in a window's haplotype.
static int nObservedSites(const string &hap){
   int n = 0;
   for (size_t i = 0; i < hap.length(); i++) if (hap[i] != MISSING_ALLELE) n++;
   return n;
}

//Deterministic processing order for the rules below: most sites observed first,
//then most frequent, then the haplotype itself. Information before
//arbitrariness -- the best-observed haplotypes establish the classes before the
//ambiguous ones are placed, and the order is a function of the data alone, so
//no seed is involved.
static void clusterOrder(const map<string,double> &combined, vector<string> &order){
   vector< pair< pair<int,double>, string > > keyed;
   keyed.reserve(combined.size());
   for (map<string,double>::const_iterator it = combined.begin(); it != combined.end(); it++)
      keyed.push_back(make_pair(make_pair(-nObservedSites(it->first), -it->second), it->first));
   sort(keyed.begin(), keyed.end());
   order.clear();
   order.reserve(keyed.size());
   for (size_t i = 0; i < keyed.size(); i++) order.push_back(keyed[i].second);
   return;
}

static void combineCounts(const map<string,double> &hap2count, const map<string,double> &miss_hap2count,
                          map<string,double> &combined){
   combined = hap2count;
   for (map<string,double>::const_iterator it = miss_hap2count.begin(); it != miss_hap2count.end(); it++)
      combined[it->first] = it->second;
   return;
}

//--hap-cluster best-comp. Each haplotype, in the deterministic order above,
//joins the MOST FREQUENT class it is compatible with rather than whichever one
//reaches it first, and only ever joins an existing class -- two class
//representatives are never merged with each other, so nothing chains
//transitively through a partial observation. Needs no fully observed haplotype
//to anchor on, which matters because in a long window at even a few percent
//missing there may not be one.
void match_haps_best_compatible(map<string,double> &hap2count, map<string,double> &miss_hap2count,
                                int len, int MATCH_TOL){
   map<string,double> combined;
   combineCounts(hap2count, miss_hap2count, combined);

   vector<string> order;
   clusterOrder(combined, order);

   vector<string> rep;      //class representative, filled in as the class absorbs
   vector<double> cnt;
   string merged;

   for (size_t i = 0; i < order.size(); i++){
      const string &hap = order[i];
      int best = -1;
      for (size_t j = 0; j < rep.size(); j++){
         if (garud_ndiff_str(rep[j], hap, merged, MATCH_TOL) <= MATCH_TOL){
            //ties broken on the representative, so the choice never depends on
            //the order the classes happen to sit in
            if (best < 0 || cnt[j] > cnt[best] ||
                (cnt[j] == cnt[best] && rep[j].compare(rep[best]) < 0)) best = (int)j;
         }
      }
      if (best < 0){
         rep.push_back(hap);
         cnt.push_back(combined[hap]);
         continue;
      }
      //when the two agree at every jointly observed site, the class label takes
      //on the sites the newcomer fills in
      if (garud_ndiff_str(rep[best], hap, merged, MATCH_TOL) == 0) rep[best] = merged;
      cnt[best] += combined[hap];
   }

   hap2count.clear();
   for (size_t j = 0; j < rep.size(); j++) hap2count[rep[j]] += cnt[j];
   return;
}

//--hap-cluster soft-em (EXPERIMENTAL). Rather than commit an ambiguous
//haplotype to one class, divide its count across every class it could have come
//from, in proportion to how common those classes are, and iterate until the
//frequencies stop moving. This is EM for the spectrum when genotypes are
//missing at random: for observation i with count c_i and compatible classes
//C(i), the E step is w_it = pi_t / sum_{t' in C(i)} pi_t' and the M step is
//pi_t = sum_i c_i w_it.
//
//The classes are the filled-in representatives from the deterministic pass, not
//the raw observed patterns. Letting a partly observed pattern be its own class
//collapses the estimate: such a pattern is compatible with the most
//observations, so the likelihood is maximised by putting all the mass on it.
//
//Class sizes come out fractional, which is why HaplotypeFrequencySpectrum
//stores doubles.
void match_haps_soft_em(map<string,double> &hap2count, map<string,double> &miss_hap2count,
                        int len, int MATCH_TOL){
   map<string,double> combined;
   combineCounts(hap2count, miss_hap2count, combined);

   map<string,double> seedHap = hap2count, seedMiss = miss_hap2count;
   match_haps_best_compatible(seedHap, seedMiss, len, MATCH_TOL);

   vector<string> types;
   vector<double> weight;
   for (map<string,double>::iterator it = seedHap.begin(); it != seedHap.end(); it++){
      types.push_back(it->first);
      weight.push_back(it->second);
   }

   vector<string> obs;
   vector<double> obsCount;
   for (map<string,double>::iterator it = combined.begin(); it != combined.end(); it++){
      obs.push_back(it->first);
      obsCount.push_back(it->second);
   }

   vector< vector<int> > compat(obs.size());
   string merged;
   double total = 0;
   for (size_t i = 0; i < obs.size(); i++){
      total += obsCount[i];
      for (size_t j = 0; j < types.size(); j++){
         if (garud_ndiff_str(types[j], obs[i], merged, MATCH_TOL) <= MATCH_TOL)
            compat[i].push_back((int)j);
      }
   }

   const int MAXIT = 200;
   const double EPS = 1e-12;
   vector<double> next(types.size(), 0.0);
   for (int iter = 0; iter < MAXIT; iter++){
      for (size_t j = 0; j < next.size(); j++) next[j] = 0;
      for (size_t i = 0; i < obs.size(); i++){
         if (compat[i].size() == 0) continue;
         double z = 0;
         for (size_t k = 0; k < compat[i].size(); k++) z += weight[compat[i][k]];
         if (z <= 0){
            double share = obsCount[i]/double(compat[i].size());
            for (size_t k = 0; k < compat[i].size(); k++) next[compat[i][k]] += share;
         }
         else{
            for (size_t k = 0; k < compat[i].size(); k++)
               next[compat[i][k]] += obsCount[i]*weight[compat[i][k]]/z;
         }
      }
      double moved = 0;
      for (size_t j = 0; j < next.size(); j++) moved += fabs(next[j] - weight[j]);
      weight = next;
      if (moved < EPS*total) break;
   }

   hap2count.clear();
   for (size_t j = 0; j < types.size(); j++)
      if (weight[j] > 0) hap2count[types[j]] += weight[j];
   //A haplotype can end up compatible with no class: the class it was placed in
   //may have had sites filled in afterwards that it disagrees with. It keeps its
   //own class, so no count is lost.
   for (size_t i = 0; i < obs.size(); i++)
      if (compat[i].size() == 0) hap2count[obs[i]] += obsCount[i];
   return;
}

HaplotypeFrequencySpectrum *hfs_window(HaplotypeData * hapData, pair_t* snpIndex, double FILTER_HMISS, int MATCH_TOL, int SEED, int CLUSTER) {
   if (numSitesInDataWin(snpIndex) <= 0) return NULL;

   HaplotypeFrequencySpectrum *hfs = initHaplotypeFrequencySpectrum();
   int haplen = snpIndex->end - snpIndex->start + 1;

   //Extract each haplotype's window in packed form: a quarter of the bytes of
   //the char form, and it is the key itself on the common path below. The
   //buffers are per-thread and kept across windows, so the string and map
   //allocations this used to make for every haplotype of every window -- all
   //of them hitting one allocator arena from every thread -- are gone.
   static thread_local vector<string> packed;
   static thread_local vector<int> nmiss;
   if ((int)packed.size() < hapData->nhaps){
      packed.resize(hapData->nhaps);
      nmiss.resize(hapData->nhaps);
   }

   bool anyMissing = false;
   for (int hap = 0; hap < hapData->nhaps; hap++) {
      extractWindow(hapData->data[hap], snpIndex->start, haplen, packed[hap]);
      nmiss[hap] = countMissingWindow(packed[hap], haplen);
      if (nmiss[hap] > 0 && double(nmiss[hap])/double(haplen) <= FILTER_HMISS) anyMissing = true;
   }

   //Common path: no haplotype in this window carries missing data and no match
   //tolerance was asked for, so nothing can be merged and only the counts are
   //needed. Tally the packed windows directly -- no char strings, no ordered
   //map, and the keys are exact, so this is not a hashing approximation.
   if (!anyMissing && MATCH_TOL == 0){
      static thread_local unordered_map<string,double> counts;
      counts.clear();
      for (int hap = 0; hap < hapData->nhaps; hap++) {
         if (double(nmiss[hap])/double(haplen) > FILTER_HMISS) continue;
         counts[packed[hap]]++;
      }
      if (counts.size() == 0) return NULL;

      double *sortedCount = new double[counts.size()];
      hfs->numClasses = counts.size();
      int i = 0;
      for (unordered_map<string,double>::iterator it = counts.begin(); it != counts.end(); it++, i++) {
         sortedCount[i] = it->second;
         hfs->size += it->second;
      }
      qsort(sortedCount, hfs->numClasses, sizeof(double), compare);
      hfs->sortedCount = sortedCount;
      return hfs;
   }

   //Clustering path: haplotypes have to be compared and rewritten site by site,
   //so unpack the windows and proceed exactly as before.
   map<string,double> miss_hap2count;
   string haplotype;

   for (int hap = 0; hap < hapData->nhaps; hap++) {
      unpackWindow(packed[hap], haplen, haplotype);

      if (double(nmiss[hap])/double(haplen) > FILTER_HMISS) continue;

      if (nmiss[hap] == 0) hfs->hap2count[haplotype]++;
      else miss_hap2count[haplotype]++;
   }

   if (CLUSTER == CLUSTER_GARUD_SHUFFLE)
      garud_match_haps_w_missing_shuffle(hfs->hap2count, miss_hap2count, haplen, MATCH_TOL,
                                         windowSeed(SEED, snpIndex->start, snpIndex->end));
   else if (CLUSTER == CLUSTER_SOFT_EM)
      match_haps_soft_em(hfs->hap2count, miss_hap2count, haplen, MATCH_TOL);
   else
      match_haps_best_compatible(hfs->hap2count, miss_hap2count, haplen, MATCH_TOL);

   if(hfs->hap2count.size() == 0) return NULL;

   double *sortedCount = new double[hfs->hap2count.size()];
   hfs->numClasses = hfs->hap2count.size();
   map<string, double>::iterator it;
   int i = 0;
   for (it = hfs->hap2count.begin(); it != hfs->hap2count.end(); it++, i++) {
      sortedCount[i] = it->second;
      hfs->size += it->second;
   }

   qsort(sortedCount, hfs->hap2count.size(), sizeof(double), compare);
   hfs->sortedCount = sortedCount;

   return hfs;
}


//Descending, for qsort over the class sizes.
int compare (const void *a, const void *b)
{
   double x = *(const double *)a, y = *(const double *)b;
   if (y > x) return 1;
   if (y < x) return -1;
   return 0;
}

int numSitesInDataWin(pair_t* win) {
   return (win->end - win->start + 1);
}
