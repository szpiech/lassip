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

void calcMandT(LASSIResults *results, SpectrumData *specData, SpectrumData *avgSpec, double **f, int w){
   double nullLikelihood = calcLASSINullLikelihood(specData,avgSpec,w);
   //cerr << "null: " << nullLikelihood << endl;
   int K = avgSpec->K;
   double U = avgSpec->freq[0][K-1];

   int maxM = -1;
   double maxE = -1;
   double maxAltLikelihood = -99999999;
   double altLikelihood = -99999999;
   double epsStep = 1.0/(100.0*double(K));
   for (int m = 1; m <= K; m++){
      for (double e = epsStep; e <= U; e += epsStep){
         altLikelihood = calcLASSIAltLikelihood(specData, avgSpec, f, U, m, e, w);
         if(altLikelihood > maxAltLikelihood){
            maxAltLikelihood = altLikelihood;
            maxM = m;
            maxE = e;
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

/*
int hamming_dist_ptr(short *one, short *two, int length)
{
   if(length == 0) return 0;

   int diff = 0;

   for(int i = 0; i < length; i++)
   {
      if(one[i] != two[i]) diff++;
   }

   return diff;
}

int hamming_dist_ptr(char *one, char *two, int length)
{
   if(length == 0) return 0;

   int diff = 0;

   for(int i = 0; i < length; i++)
   {
      if(one[i] != two[i]) diff++;
   }

   return diff;
}

int hamming_dist_str(string one, string two, pair_t *subset_snps) {
   int diff = 0, start, end;
   if(subset_snps == NULL){
      start = 0;
      end = one.length()-1;
   }
   else{
      start = subset_snps->start;
      end = subset_snps->end;
   }
   for (int i = start; i <= end; i++)
   {
      if (one[i] != two[i]) diff++;
   }
   return diff;
}

double ehh_from_hfs(HaplotypeFrequencySpectrum * hfs) {
   if (hfs == NULL) return MISSING;
   map<string, int>::iterator it;
   double tot = 0;
   double homozygosity = 0;
   for (it = hfs->hap2count.begin(); it != hfs->hap2count.end(); it++) {
      tot += it->second;
      homozygosity += (it->second > 1) ? nCk(it->second, 2) : 0;
   }
   homozygosity /= nCk(tot, 2);
   return homozygosity;
}

double ehhk_from_hfs(HaplotypeFrequencySpectrum * hfs, int k) {
   if (hfs == NULL) return MISSING;
   double res = 0;
   double homozygosity = 0;
   double tot = 0;
   map<string, int>::iterator it;
   int *sortedCounts = new int[hfs->hap2count.size()];
   int i = 0;
   for (it = hfs->hap2count.begin(); it != hfs->hap2count.end(); it++) {
      tot += it->second;
      homozygosity += (it->second > 1) ? nCk(it->second, 2) : 0;
      sortedCounts[i] = it->second;
      i++;
   }

   qsort(sortedCounts, hfs->hap2count.size(), sizeof(int), compare);

   res = homozygosity;
   int combined = 0;
   int maxK = k = (hfs->numUniq < k) ? hfs->numUniq : k;
   for (int i = 0; i < maxK; i++) {
      combined += sortedCounts[i];
      res -= (sortedCounts[i] > 1) ? nCk(sortedCounts[i], 2) : 0;
   }

   res += (combined > 1) ? nCk(combined, 2) : 0;

   delete [] sortedCounts;

   return (res / nCk(tot, 2));
}

double pi_numerator_btw_pools(string * haps1, int length1, string * haps2, int length2, map<string, int> &hap2count, pair_t *subset_snps) {
   double num = 0;

   for (int i = 0; i < length1; i++) {
      for (int j = 0; j < length2; j++) {
         num += hamming_dist_str(haps1[i], haps2[j], subset_snps) * hap2count[haps1[i]] * hap2count[haps2[j]];
      }
   }

   return num;
}

double pi_numerator(string * haps, int length, map<string, int> &hap2count, pair_t *subset_snps) {
   double num = 0;

   for (int i = 0; i < length; i++) {
      for (int j = i + 1; j < length; j++) {
         num += hamming_dist_str(haps[i], haps[j], subset_snps) * hap2count[haps[i]] * hap2count[haps[j]];
      }
   }

   return num;
}
*/
/*
   Here we define pi_k as pi restricted to haplotypes belonging to the first
   k haplotype classes.  For example consider the following haplotype
   spectrum:

   count   hap
   10      011
    5      001
    3      111
    1      000

   Then for k = 2 we calculate pi amongst the following subset of haplotypes:

   count   hap
   10      011
    5      001

   And for k = 3 we calculate pi amonst the following subset of haplotypes:

   count   hap
   10      011
    5      001
    3      111

   What do we do if there are ties in a frequency class?  For example:

   count   hap
   10      011
    5      001
    5      100
    3      111
    1      000

   Then for k = 2 we calculate the average pi amongst the following pairs subset of haplotypes:

   count   hap
   10      011
    5      001

    and

    count   hap
   10      011
    5      100

   What if there are > 2 ties? For example:

   count   hap
   10      011
    5      001
    5      100
    5      111
    1      000

   Here there is a three-way tie for the 2nd most frequent haplotype.  For k = 2, we
   calculate pi amongst the top two most frequent classes for each of
   the three choices of haplotypes and report the mean.  For k = 3, we must choose 2
   of the tied haplotype classes to represent the 2nd and 3rd most frequent haplotypes.
   There are 3 choose 2 ways to make this choice, and we take the mean across them.

*/
/*
double pi_k2(HaplotypeFrequencySpectrum * hfs, int k, pair_t *subset_snps) {//this is assumed to be offset
   if (hfs == NULL || (hfs->numUniq < k)) return MISSING;

   pair <multimap<int, string>::iterator, multimap<int, string>::iterator> ret;
   multimap<int, string>::iterator it;

   //k = (hfs->numUniq < k) ? hfs->numUniq : k;

   string *haps = new string[k];

   //counts the number of unique hap classes upto k
   int howmanyUniqHaps = 0;
   //If the kth most frequent haplotype class has > 1 haplotype associated with it, this counts how many
   int numNextClass = 0;
   //If the kth most frequent haplotype class has > 1 haplotype associated with it, this stores the ties
   string *equalFreqHaps;//length == numNextClass
   int h = 0;
   for (int i = 0; i < k; i++) {
      int key = hfs->sortedCount[i];
      numNextClass = hfs->count2hap.count(key);

      if (howmanyUniqHaps + numNextClass > k) {
         equalFreqHaps = new string[numNextClass];
         int j = 0;
         ret = hfs->count2hap.equal_range(key);
         for (it = ret.first; it != ret.second; it++) {
            equalFreqHaps[j] = it->second;
            j++;
         }
         break;
      }

      howmanyUniqHaps += numNextClass;
      ret = hfs->count2hap.equal_range(key);

      for (it = ret.first; it != ret.second; it++) {
         haps[h] = it->second;
         h++;
      }
      numNextClass = 0;
   }

   int numHapsMissing = k - howmanyUniqHaps;

      //cout << "-----\n";
      //for (int i = 0; i < howmanyUniqHaps; i++) {
      //   cout << "  " << haps[i] << " " << hfs->hap2count[haps[i]] << endl;
      //}
   

   double pi = 0;
   int nhaps = 0;
   double denominator;

   if (numHapsMissing == 0) {
      for (int i = 0; i < k; i++) nhaps += hfs->hap2count[haps[i]];
      denominator = (nhaps) * (nhaps - 1) * 0.5;
      //cout << pi / denominator << endl;
      return pi_numerator(haps, k, hfs->hap2count, subset_snps) / denominator;
   }
   else {
      for (int i = 0; i < k; i++) {
         if (i < howmanyUniqHaps) {
            nhaps += hfs->hap2count[haps[i]];
         }
         else {
            nhaps += hfs->hap2count[equalFreqHaps[i - howmanyUniqHaps]];
         }
      }
      denominator = (nhaps) * (nhaps - 1) * 0.5;
      double pi_partial = pi_numerator(haps, howmanyUniqHaps, hfs->hap2count, subset_snps);
      gsl_combination * c;
      c = gsl_combination_calloc (numNextClass, numHapsMissing);

      //cout << "--next class--\n";
      string *chosenHaps = new string[numHapsMissing];
      //double pi_combo = 0;
      do
      {
         //cout << "[ ";
         for (int i = 0; i < numHapsMissing; i++) {
            chosenHaps[i] = equalFreqHaps[gsl_combination_get(c, i)];
            //cout << chosenHaps[i] << "\n  ";
         }
         //cout << "] " << hfs->hap2count[chosenHaps[0]] << " : ";
         pi += pi_partial +
               pi_numerator_btw_pools(haps, howmanyUniqHaps, chosenHaps, numHapsMissing, hfs->hap2count, subset_snps) +
               pi_numerator(chosenHaps, numHapsMissing, hfs->hap2count, subset_snps);
         
         //pi_combo = pi_partial +
         //           pi_numerator_btw_pools(haps, howmanyUniqHaps, chosenHaps, numHapsMissing, hfs->hap2count) +
         //           pi_numerator(chosenHaps, numHapsMissing, hfs->hap2count);
         //cout << pi_combo << " / " << denominator << " -> " << pi_combo / denominator << endl;
         

      } while (gsl_combination_next (c) == GSL_SUCCESS);
      pi /= nCk(numNextClass, numHapsMissing);
      gsl_combination_free (c);
      delete [] chosenHaps;
   }

   //cout << pi / denominator << endl;


   if (equalFreqHaps != NULL) {
      delete [] equalFreqHaps;
   }

   delete [] haps;

   return pi / denominator;
}
*/
/*
double pi_k(HaplotypeFrequencySpectrum * hfs, int k) {
   if (hfs == NULL) return MISSING;
   k = (hfs->numUniq < k) ? hfs->numUniq : k;
   string *haps = new string[k];

   //Grab the first k most frequent haplotypes
   //Right now if there are ties in the final one
   //we only take the first few upto k total unique haplotypes
   //In the future, will take the mean
   int i = 0;
   for (int j = 0; j < hfs->size; j++) {
      int key = hfs->sortedCount[j];
      pair <multimap<int, string>::iterator, multimap<int, string>::iterator> ret;
      ret = hfs->count2hap.equal_range(key);
      multimap<int, string>::iterator it;
      for (it = ret.first; it != ret.second; it++) {
         haps[i] = it->second;
         i++;
         if (i >= k) break;
      }
      if (it != ret.second) break;
   }

   double pi_k = 0;
   int nhaps = 0;
   for (i = 0; i < k; i++) nhaps += hfs->hap2count[haps[i]];

   double denominator = (nhaps) * (nhaps - 1) * 0.5;

   for (i = 0; i < k; i++) {
      for (int j = i + 1; j < k; j++) {
         pi_k += hamming_dist_str(haps[i], haps[j]) * hfs->hap2count[haps[i]] * hfs->hap2count[haps[j]];
      }
   }

   delete [] haps;

   return pi_k / denominator;
}
*/
HaplotypeFrequencySpectrum *hfs_window(HaplotypeData * hapData, pair_t* snpIndex) {
   if (numSitesInDataWin(snpIndex) <= 0) return NULL;

   HaplotypeFrequencySpectrum *hfs = initHaplotypeFrequencySpectrum();

   bool skip = false;
   //Generate haplotypes and populate hap2count
   for (int hap = 0; hap < hapData->nhaps; hap++) {
      string haplotype;

      for (int site = snpIndex->start; site <= snpIndex->end; site++) {
         if (hapData->data[hap][site] == MISSING_ALLELE){
            skip = true;
            break;
         }
         if (site == snpIndex->start) {
            //haplotypeList[hap] = data[hap][site];
            haplotype = hapData->data[hap][site];
         }
         else {
            //haplotypeList[hap] += data[hap][site];
            haplotype += hapData->data[hap][site];
         }
      }
      //cerr << haplotype << endl;
      if (!skip){
         if (hfs->hap2count.count(haplotype) == 0) {
            hfs->hap2count[haplotype] = 1;
         }
         else {
            hfs->hap2count[haplotype]++;
         }
      }
      else{
         skip = false;
      }
   }

   if(hfs->hap2count.size() == 0) return NULL;

   //Populate count2hap and sortedCounts
   int *sortedCount = new int[hfs->hap2count.size()]; //could contain duplicates
   hfs->numClasses = hfs->hap2count.size();
   map<string, int>::iterator it;
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
/*
double pi_window(HaplotypeData * hapData, pair_t* snpIndex) {
   if (numSitesInDataWin(snpIndex) <= 0) return MISSING;
   //int startSnpIndex; int endSnpIndex;
   double pi = 0;
   double denominator = (hapData->nhaps) * (hapData->nhaps - 1) * 0.5;
   int length = snpIndex->end - snpIndex->start + 1;
   if (length == 0) {
      pi = 0;
   }
   else {
      for (int i = 0; i < hapData->nhaps; i++) {
         for (int j = i + 1; j < hapData->nhaps; j++) {
            pi += hamming_dist_ptr(hapData->data[i] + snpIndex->start, hapData->data[j] + snpIndex->start, length);
         }
      }
   }
   return (pi / denominator);
}

double subsample_sfs(array_t *sfs, int H, int j){
   double res = 0;
   int n = sfs->size;
   for (int i = j; i < n-1; i++){
      res += sfs->data[i] * nCk(i,j)*nCk(n-i,H-j)/nCk(n,H);
   }
   return res;
}

array_t *sfs_window(FreqData * freqData, pair_t* snpIndex, bool SFS_SUB) {
   if (numSitesInDataWin(snpIndex) <= 0) return NULL;
   int nTargetHaps = freqData->nhaps;
   vector<int> nhaps;
   int n;
   if(SFS_SUB){
      for (int i = snpIndex->start; i <= snpIndex->end; i++){
         n = freqData->nhaps - freqData->nmissing[i];
         nhaps.push_back(n);
         if (nTargetHaps > n){
            nTargetHaps = n;
         }
      }
   }

   array_t *sfs = initArray(nTargetHaps + 1);
   
   if (nTargetHaps != freqData->nhaps){
      array_t *s;
      map<int,array_t*> multiSFS;
      int j;
      for (int i = 0; i < nhaps.size(); i++){
         n = nhaps[i];
         j = i + snpIndex->start;
         if (multiSFS.count(n) == 0) multiSFS[n] = initArray(n + 1);
         multiSFS[n]->data[freqData->count[j]]++;
      }
      map<int,array_t*>::iterator it;
      for (it = multiSFS.begin(); it != multiSFS.end(); it++){
         n = it->first;
         s = it->second;
         if(n == nTargetHaps){
            for (int i = 0; i < sfs->size; i++){
               sfs->data[i] += s->data[i];
            }
         }
         else{
            for (int i = 0; i < sfs->size; i++){
               sfs->data[i] += subsample_sfs(s,nTargetHaps,i);
            }
         }
         releaseArray(s);
      }
   }
   else{
      for (int i = snpIndex->start; i <= snpIndex->end; i++) {
         sfs->data[freqData->count[i]]++;
      }
   }

   return sfs;
}
*/

/*
array_t *sfs_window(FreqData * freqData, pair_t* snpIndex) {
   if (numSitesInDataWin(snpIndex) <= 0) return NULL;
   array_t *sfs = initArray(freqData->nhaps + 1);

   for (int i = snpIndex->start; i <= snpIndex->end; i++) {
      sfs->data[freqData->count[i]]++;
   }

   return sfs;
}
*/
/*
double pi_from_sfs(array_t *sfs) {
   if (sfs == NULL) return MISSING;
   double pi = 0;
   int n = sfs->size - 1;
   double denominator = n * (n - 1) * 0.5;

   for (int i = 1; i < n; i++) {
      pi += i * (n - i) * sfs->data[i];
   }
   return pi / denominator;
}

double fayWuH_from_sfs(array_t *sfs, double pi) {
   if (sfs == NULL) return MISSING;
   if (pi < 0) {
      pi = pi_from_sfs(sfs);
   }
   return (pi - thetaH_from_sfs(sfs));
}

double thetaH_from_sfs(array_t *sfs) {
   if (sfs == NULL) return MISSING;
   double thetaH = 0;
   int n = sfs->size - 1;
   double denominator = n * (n - 1) * 0.5;

   for (int i = 1; i < n; i++) {
      thetaH += i * i * sfs->data[i];
   }
   return thetaH / denominator;
}

double tajimaD_from_sfs(array_t *sfs, double pi, double S) {
   if (sfs == NULL) return MISSING;
   if (pi < 0) {
      pi = pi_from_sfs(sfs);
   }

   if (S < 0) {
      S = segsites(sfs);
   }
   int n = sfs->size - 1;
   double e1, e2, a1, a2, denominator;

   a1 = calc_a1(n);
   a2 = calc_a2(n);
   e1 = calc_e1(n, a1);
   e2 = calc_e2(n, a1, a2);

   denominator = e1 * S + e2 * S * (S - 1);

   return (pi - S / a1) / denominator;
}

int segsites(array_t *sfs) {
   if (sfs == NULL) return MISSING;
   double s = 0;
   int n = sfs->size - 1;
   for (int i = 1; i < n; i++) {
      s += sfs->data[i];
   }
   return s;
}

double s_from_sfs(array_t *sfs) {
   if (sfs == NULL) return MISSING;
   double s = 0;
   int n = sfs->size - 1;

   for (int i = 1; i < n; i++) {
      s += sfs->data[i];
   }
   return s;
}

double calc_a1(int n) {
   double a = 0;
   for (double i = 1; i < n; i++)
      a += 1.0 / i;
   return a;
}

double calc_a2(int n) {
   double a = 0;
   for (double i = 1; i < n; i++)
      a += 1.0 / (i * i);
   return a;
}

double calc_e1(int n, double a1) {
   return ( (n + 1.0) / (3.0 * n - 3.0) - (1.0 / a1) ) / a1;
}

double calc_e2(int n, double a1, double a2) {
   return ( (2.0 * n * n + 2.0 * n + 6.0) / (9.0 * n * (n - 1)) - (n + 2) / (n * a1) + a2 / (a1 * a1)) / (a1 * a1 + a2);
}
*/
int numSitesInDataWin(pair_t* win) {
   return (win->end - win->start + 1);
}
