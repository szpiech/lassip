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

double getDMin(vector<SpectrumData *> *specDataByChr){
   double dmin = 9999999999;

   for(unsigned int c = 0; c < specDataByChr->size(); c++){
      for(int w = 1; w < specDataByChr->at(c)->nwins; w++){
         double diff = abs(specDataByChr->at(c)->dist[w]-specDataByChr->at(c)->dist[w-1]);
         if(dmin > diff && diff > 0) dmin = diff;
      }
   }
   return dmin;
}

void calcQ(double ***q, SpectrumData *avgSpec, double **f, int w){
   int K = avgSpec->K;
   double U = avgSpec->freq[0][K-1];
   double epsStep = 1.0/(100.0*double(K));
   int ei = 0;
   for (double e = epsStep; e <= U; e += epsStep){
      for (int m = 1; m < K; m++){
         calcQ(q[ei][m-1], avgSpec, f, U, m, e, w);
      }
      ei++;
   }
   return;
}
void calcQ(double *q, SpectrumData *avgSpec, double **f, double U, int m, double e, int w){
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
   int *c = hfs->sortedCount;
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
   int *c = hfs->sortedCount;
   double tot = hfs->size;
   double first = (double(c[0])/tot)*(double(c[0])/tot);
   double res = first;
   for(int i = 1; i < hfs->numClasses; i++){
      res += (double(c[i])/tot)*(double(c[i])/tot);
   }
   return (res-first)/res;
}


void calcMTA(LASSIResults *results, double ****q, SpectrumData *specData, SpectrumData *avgSpec, int w, double dmin,double MAX_EXTEND){
   //int MAX_EXTEND = 2500000;
   int rightLim, leftLim;
   double *dist = specData->dist;
   int d = w;
   while(abs(dist[w] - dist[d]) <= MAX_EXTEND){
      d++;
      if(d >= specData->nwins){
         d--;
         break;
      }
   }
   leftLim = d;
   d = w;
   while(abs(dist[w] - dist[d]) <= MAX_EXTEND){
      d--;
      if(d < 0){
         d++;
         break;
      }
   }
   rightLim = d;
   
   double nullLikelihood = calcSALTINullLikelihood(specData,avgSpec,w,rightLim,leftLim);
   //cerr << "null: " << nullLikelihood << endl;
   int K = avgSpec->K;
   double U = avgSpec->freq[0][K-1];

   //for (int i = 0; i < K; i++) cerr << avgSpec->freq[0][i] << " ";
   //cerr << endl;

   int maxM = -1;
   //double maxE = -1;
   double maxA = -1;
   double maxAltLikelihood = -99999999;
   double altLikelihood = -99999999;
   double epsStep = 1.0/(100.0*double(K));

   double Amin = -log(0.99999)/dmin;
   double Amax = -log(0.00001)/dmin;
   double lAmin = log(Amin);
   double lAmax = log(Amax);
   double lstep = (lAmax-lAmin)/100;

   //cerr << "lAmax " << lAmax << " lAmin " << lAmin << " lstep " << lstep << endl;

   //double p = 0-log(1e-8);
   //double d;

   for (double A = lAmin; A <= lAmax; A += lstep){
      //cerr << "p " << p << " A " << A << " d " << d << endl;
      //double currNullLikelihood = calcSALTINullLikelihood(specData,avgSpec,w,d);
      for (int m = 1; m <= K; m++){
         int ei = 0;
         for (double e = epsStep; e <= U; e += epsStep){
            altLikelihood = calcSALTIAltLikelihood(specData, avgSpec, q, ei, m-1, exp(A), w, rightLim, leftLim);
            //cerr << "A " << A << " m " << m << " e " << " alt " << altLikelihood << endl;
            if(altLikelihood > maxAltLikelihood){
               maxAltLikelihood = altLikelihood;
               //nullLikelihood = currNullLikelihood;
               maxM = m;
               //maxE = e;
               maxA = 1.0/exp(A);
            }
            ei++;
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
   return;
}


double calcSALTINullLikelihood(SpectrumData *specData,SpectrumData *avgSpec,int w,int rightLim, int leftLim){
   double res = 0;
   //int start = (w-width >= 0) ? w-width : 0;
   //int end = (w+width < specData->nwins) ? w+width : specData->nwins-1;
   //double center = specData->dist[w];
   /*
   //to the left
   int win = w-1;
   while(center - specData->dist[win] <= d && win >= 0){
      for (int i = 0; i < avgSpec->K; i++){
         res += double(specData->nhaps[win])*specData->freq[win][i]*log(avgSpec->freq[0][i]);
      }
      win--;
   }
   //to the right
   win = w;
   while(specData->dist[win] - center  <= d && win < specData->nwins){
      for (int i = 0; i < avgSpec->K; i++){
         res += double(specData->nhaps[win])*specData->freq[win][i]*log(avgSpec->freq[0][i]);
      }
      win++;
   }
   */

   for (int win = rightLim; win <= leftLim; win++){
      for (int i = 0; i < avgSpec->K; i++){
         res += double(specData->nhaps[win])*specData->freq[win][i]*log(avgSpec->freq[0][i]);
      }
   }
   return res;
}

double calcSALTIAltLikelihood(SpectrumData *specData,SpectrumData *avgSpec,double ****q,int e, int m, double A, int w,int rightLim, int leftLim){
   if(m+1 == avgSpec->K) return calcSALTINullLikelihood(specData,avgSpec,w,rightLim,leftLim);
   double res = 0;
   //int start = (w-width >= 0) ? w-width : 0;
   //int end = (w+width < specData->nwins) ? w+width : specData->nwins-1;
   //double center = specData->dist[w];

   /*
   //to the left
   int win = w-1;
   while(center - specData->dist[win] <= d && win >= 0){
      double Pr = exp(-A*abs(center-specData->dist[win]));
      for (int i = 0; i < avgSpec->K; i++){
         double pi = double(specData->nhaps[win])*specData->freq[win][i]*log(avgSpec->freq[0][i]);
         double qi = double(specData->nhaps[win])*specData->freq[win][i]*log(q[win][e][m][i]);
         res += Pr*qi+(1-Pr)*pi;
      }
      win--;
   }
   //to the right
   win = w;
   while(specData->dist[win] - center <= d && win < specData->nwins){
      double Pr = exp(-A*abs(center-specData->dist[win]));
      for (int i = 0; i < avgSpec->K; i++){
         double pi = double(specData->nhaps[win])*specData->freq[win][i]*log(avgSpec->freq[0][i]);
         double qi = double(specData->nhaps[win])*specData->freq[win][i]*log(q[win][e][m][i]);
         res += Pr*qi+(1-Pr)*pi;
      }
      win++;
   }
   */

   for (int win = rightLim; win <= leftLim; win++){
      double Pr = exp(-A*abs(specData->dist[w]-specData->dist[win]));
      for (int i = 0; i < avgSpec->K; i++){
         double pi = double(specData->nhaps[win])*specData->freq[win][i]*log(avgSpec->freq[0][i]);
         double qi = double(specData->nhaps[win])*specData->freq[win][i]*log(q[win][e][m][i]);
         res += Pr*qi+(1-Pr)*pi;
      }
   }
   return res;
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

int ndiff_str(string str1, string str2){
   int ndiff = 0;
   bool all_missing = true;
   string::iterator i;
   string::iterator j;
   for (i = str1.begin(); i != str1.end(); i++){
      for (j = str2.begin(); j != str2.end(); j++){
         if (*i == MISSING_ALLELE || *j == MISSING_ALLELE) continue;
         if (*i != *j){
            ndiff++;
            all_missing = false;
         }
      }
   }
   return all_missing ? -1 : ndiff;
}


void match_haps_w_missing(map<string,double> &hap2count,map<string,double> &miss_hap2count, int len, int MATCH_TOL){
   map<string, double>::iterator it1;
   map<string, double>::iterator it2;
   
   vector<string> best_matches;
   vector<string> to_delete;
   int mindiff = len+1;

   for (it2 = miss_hap2count.begin(); it2 != miss_hap2count.end(); it2++){
      for (it1 = hap2count.begin(); it1 != hap2count.end(); it1++){
         
         int d = ndiff_str(it2->first,it1->first);

         if(d <= MATCH_TOL && d >= 0){
            if(d == mindiff){
               best_matches.push_back(it1->first);
            }
            else if (d < mindiff){
               mindiff = d;
               best_matches.clear();
               best_matches.push_back(it1->first);
            }
         }
         else{
            continue;
         }
      
      }
   
      if(best_matches.size() > 0){
         for (int i = 0; i < best_matches.size(); i++){
            hap2count[best_matches[i]] += (it2->second/double(best_matches.size()));
         }
         to_delete.push_back(it2->first);
      }
      best_matches.clear();
      mindiff = len+1;
   }

   for (int i = 0; i < to_delete.size(); i++){
      miss_hap2count.erase(to_delete[i]);
   }

   if (miss_hap2count.empty()) return;

   //naive clustering of incomplete haplotypes with each other
   //In principle results could be dependent on order processed and clusters may not strictly conform to MATCH_TOL
   //Can be improved at the expense of computational burden

   map<string,double> miss_hap2count2 = miss_hap2count;
   map<string,double> miss_clustered2count;
   best_matches.clear();
   //to_delete.clear();
   mindiff = len+1;
   
   int i0 = 0;
   int j0 = 0;
   for (it2 = miss_hap2count2.begin(); it2 != miss_hap2count2.end(); it2++){
      for (it1 = miss_hap2count.begin(); it1 != miss_hap2count.end(); it1++){
         
         if (i0 == j0) continue;

         int d = ndiff_str(it2->first,it1->first);

         if(d <= MATCH_TOL && d >= 0){
            if(d == mindiff){
               best_matches.push_back(it1->first);
            }
            else if (d < mindiff){
               mindiff = d;
               best_matches.clear();
               best_matches.push_back(it1->first);
            }
         }
         else{
            continue;
         }
         j0++;
      }
      
      if(best_matches.size() > 0){
         for (int i = 0; i < best_matches.size(); i++){
            miss_clustered2count[best_matches[i]] = hap2count[best_matches[i]];
            miss_clustered2count[best_matches[i]] += (it2->second/double(best_matches.size()));
         }
         miss_hap2count.erase(it2->first);
      }
      best_matches.clear();
      mindiff = len+1;
      i0++;
   }

   for (it2 = miss_clustered2count.begin(); it2 != miss_clustered2count.end(); it2++){
      hap2count[it2->first] = it2->second;
   }

   return;

}

void printHFS(map<string,double> hap2count){
   map<string,double>::iterator it;
   for (it = hap2count.begin(); it != hap2count.end(); it++){
      cout << it->first << "\t" << it->second << endl;
   }
}

HaplotypeFrequencySpectrum *hfs_window(HaplotypeData * hapData, pair_t* snpIndex, double FILTER_HMISS, int MATCH_TOL) {
   if (numSitesInDataWin(snpIndex) <= 0) return NULL;

   HaplotypeFrequencySpectrum *hfs = initHaplotypeFrequencySpectrum();

   bool skip = false;
   //Generate haplotypes and populate hap2count

   //map<string,double> nomiss_hap2count;
   map<string,double> miss_hap2count;
   int nmissing = 0;
   int haplen = snpIndex->end - snpIndex->start + 1;
   
   for (int hap = 0; hap < hapData->nhaps; hap++) {
      string haplotype;

      nmissing = 0;
      
      for (int site = snpIndex->start; site <= snpIndex->end; site++) {
         if (hapData->data[hap][site] == MISSING_ALLELE){
            //skip = true;
            //break;
            nmissing++;
         }
         if (site == snpIndex->start) {
            haplotype = hapData->data[hap][site];
         }
         else {
            haplotype += hapData->data[hap][site];
         }
      }

      skip = (double(nmissing)/double(haplen) > FILTER_HMISS);

      if (!skip){
         if(nmissing == 0){
            if (hfs->hap2count.count(haplotype) == 0) {
               hfs->hap2count[haplotype] = 1;
            }
            else {
               hfs->hap2count[haplotype]++;
            }
         }
         else{
            if(miss_hap2count.count(haplotype) == 0){
               miss_hap2count[haplotype] = 1;
            }
            else{
               miss_hap2count[haplotype]++;
            }
         }
      }
      else{
         skip = false;
      }
   }

   match_haps_w_missing(hfs->hap2count, miss_hap2count, haplen, MATCH_TOL);

   if(hfs->hap2count.size() == 0) return NULL;

   //Populate count2hap and sortedCounts
   int *sortedCount = new int[hfs->hap2count.size()]; //could contain duplicates
   hfs->numClasses = hfs->hap2count.size();
   map<string, double>::iterator it;
   int i = 0;
   for (it = hfs->hap2count.begin(); it != hfs->hap2count.end(); it++, i++) {
      sortedCount[i] = it->second;//unsorted
      hfs->size += it->second;
      //hfs->count2hap.insert(pair<int, string>(it->second, it->first));
   }

   qsort(sortedCount, hfs->hap2count.size(), sizeof(int), compare);//sorted but with possible duplicates
   //hfs->sortedCount = uniqInt(sortedCount, hfs->numClasses, hfs->size);//remove duplicates
   hfs->sortedCount = sortedCount;
   //delete [] sortedCount;

   //printHFS(hfs->hap2count);

   return hfs;
}

int *uniqInt(int *array, int size, int &newSize) {
   map<int, int> uniq;
   for (int i = 0; i < size; i++) {
      uniq[array[i]] = 1;
   }
   newSize = uniq.size();
   int *newArray = new int[newSize];
   int prev = array[0];
   int j = 0;
   newArray[j] = prev;
   j++;
   for (int i = 1; i < size; i++) {
      if (array[i] != prev) {
         prev = array[i];
         newArray[j] = prev;
         j++;
      }
   }
   return newArray;
}

int compare (const void *a, const void *b)
{
   return ( *(int *)b - * (int *)a );
}

int numSitesInDataWin(pair_t* win) {
   return (win->end - win->start + 1);
}
