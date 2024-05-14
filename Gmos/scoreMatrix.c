/***** matrix.c *************************************************************
* Description: Source file for computing DNA scoring matrix (GTR model).
*
* Author: Mirjana Domazet-Loso
*
* This file is part of gmos.
*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

#include "commonSC.h"
#include "eprintf.h"
#include "stringUtil.h"
#include "intervalStack.h"
#include "interval.h"

#include "interface.h"
#include "nw.h"
#include "queryBTNode.h"
#include "lNode.h"
#include "align.h"
#include "mergeList.h"
#include "scoreMatrix.h"

/****  Simple example for computing the rate matrix:
* Seq1 A C G T
* Seq2 A A A T
* Seq3 A A G T
*
* For each pair of nucleotides, freq. of occurence is computed (also, AC is treated as CA); 
* e.g. freq(AA) = 3 x 2 / 2 + 2 x 1 / 2 = 4; freq(AC) = 2 x 1 = 2;
* Observed_frequency(N1, N2) = freq(N1, N2) / sum_of_all_observed_freq
* 
* Expected_frequency(AA) = 6 * 6 / (12 * 12)
* Expected_frequency(AC) = 2 * (6 / 12) * (1 / 12); a factor of 2 because A can be aligned to C and vice versa
* Matrix(N1, N2) = rounded off value of (2 x log2(Observed_frequency(N1, N2) / Expected_frequency(N1, N2)))
*/

/*** General idea for computation (dna scoring matrix, gtr model, is computed from the existing alignments):
* (1) array queryMatch[queryLength] --> counts the number of nucleotides in subject regions matching each query nucleotide ( + 1 for a query nucleotide)
* (2) listMismatches --> list of mismathces for query positions (only query positions for which exist mismatches are included in this list)
* (3) freqNucl[alphabetSize] --> count each nucleotide over all msa-s along the query genome
* (4) sumFreqNucl --> sum of all elements in freqNucl
* (5) observedFreq --> matrix of observed freq-s
* (6) expectedFreq --> matrix of expected freq-s
* (7) scoreMatrix --> dna scoring matrix; gtr model

* (1) for each ar (aligned_region in mNode) in fnalList
*    for each match within ar
*       for each nucleotide N within match (N matches a query nucleotide querySeq[i])
*         ++ queryMatch[i]
*    for each aligned region between matches 
*      for each nucleotide N
*         if N matches a query nucleotide querySeq[i] then
*            ++ queryMatch[i]
*         else
*            ++ listMismatches[i][N]
*
* (2) compute total number of nucleotides in msa: sum of all elements in queryMatch + sum of all elements in listMismatches
* (3) compute matrix of observed freq-s from queryMatch and listMismatches;
*    for each i in queryMatch 
*      N = querySeq[i]
*      observedFreq(N, N) += queryMatch[i] * (queryMatch[i] - 1) / 2
*    for each i in listMismatches[i]
*       M = querySeq[i]
*       for each N in listMismatches[i]
*         observedFreq(M, N) += listMismatches[i][N] * queryMatch[i];
*
* (4) compute sum_of_all_observed_freq: sum of all elements in upper triangular matrix (including diagonal) of observed frequencies
* (assuming that observedFreq(M, N) = observedFreq(N, M))
*   
* (5) compute expected freq-s
*   expectedFreq(N, N) = freqNucl[N] * freqNucl[N] / (sumFreqNucl * sumFreqNucl)
*   expectedFreq(M, N) = 2 * freqNucl[M] * freqNucl[N] / (sumFreqNucl * sumFreqNucl) 
*
* (6) compute gtr matrix: scoreMatrix(N1, N2) = rounded off value of (2 x log2(observedFreq(N1, N2) / expectedFreq(N1, N2)))
*/


/* compute |Q| scoring matrices and nucleotide freq-s
* number of columns in scoreMatrix is |alphabetSize|^2 */
void getScoringMatrices(mNode **finalList, Int64 numOfQueries, char *seq, Int64 *leftBorders, Int64 *strandBorders,
						Int64 *strandLen, float **scoreMatrixAll, 
						int	alphabetSize, double ***nuclFreq, float **expectedFreq) {
  Int64 i, j, k, l;
	mNode *m = NULL;
	int maxListMismatches = 1000; // resize when necessary
	Mismatch *listMismatches = NULL;
  long long totalNumNucleotides = 0;
  double observedFreq[ALPHABET_DNA_SIZE][ALPHABET_DNA_SIZE] = { { 0.0 } }; // observed frequencies
  //double expectedFreq[ALPHABET_DNA_SIZE][ALPHABET_DNA_SIZE] = { { 0.0 } }; // expected frequencies
  
  long long sumObservedFreq = 0; // sum of all observed freq-s (sum of all elements in upper triangular matrix (including diagonal) of observed frequencies)
	Int64 occNucl[ALPHABET_DNA_SIZE] = { 0 }; // number of occurences of a nucleotide (A-T --> indices: 0-4)

  listMismatches = getListMismatches(maxListMismatches, alphabetSize);
	*nuclFreq = (double **)/*e*/malloc(sizeof(double *) * numOfQueries);

	for (i = 0; i < numOfQueries; i++) { // for each query
		m = finalList[i];
		(*nuclFreq)[i] = (double *)/*e*/malloc(sizeof(double) * alphabetSize);
		// queryMatch: count nucleotides in the forward query strand
		short *queryMatch = (short *)/*e*/malloc(sizeof(short) * strandLen[i]);
		int numMismatches = 0;		
		for (j = 0; j < strandLen[i]; j++) {
			queryMatch[j] = 0; 
		}
		// initializiaton
		initializeListMismatches(maxListMismatches, alphabetSize, listMismatches);
		for (j = 0; j < alphabetSize; j++) {
			occNucl[j] = 0;
			for (k = 0; k < alphabetSize; k++) {
				observedFreq[j][k] = 0;
				//expectedFreq[j][k] = 0;
				expectedFreq[i][j * alphabetSize + k] = 0;
			}
		}
    // (1): for each query: compute listMismatches and queryMatch
    while (m && m->elem->endPosQ <= strandLen[i]) { // for each mNode / aRegion in finalList[i]
			aNode *a = m->elem->start;
			int iAlignedSeg = 0; // index of aligned segment between exact matches
			while (a) { // for each exact match in m and each aligned region between two exact matches (if exists)
				// if this is the first match for a query column; then +2 (+1 to count a query nucleotide and +1 for subject), otherwise only +1 for subject nucleotide
				for (k = a->lbQ; k < a->lbQ + a->sl - 1; k++) {
					if (seq[leftBorders[i] + k] != toupper(AMBIGUOUS_CHAR_QUERY)) { // ignore 'X'
						queryMatch[k] += (queryMatch[k] == 0) ? 2 : 1;
					}
				}
				if (a->next) { // if there is an aligned region between two matches
          //	if a subject nucleotide (N) matches a query nucleotide at the position i, then ++queryMatch[i]
          //	else ++listMismatches[i][N]
          int k2 = k; // k equals the first position after the last match (in a)
					int lenAlignedSegment = strlen(m->elem->alignedSegment[iAlignedSeg]); // query segment: iAlignedSeg; subject segment: iAlignedSeg + 1
					for (l = 0; l < lenAlignedSegment; l++) {
            if (m->elem->alignedSegment[iAlignedSeg][l] != '-' && m->elem->alignedSegment[iAlignedSeg + 1][l] != '-' 
							&& m->elem->alignedSegment[iAlignedSeg][l] != toupper(AMBIGUOUS_CHAR_QUERY) 
							&& m->elem->alignedSegment[iAlignedSeg + 1][l] != toupper(AMBIGUOUS_CHAR_SUBJECT)) { // if it is not a gap or X/Y character on either side
              if (queryMatch[k2 + l] == 0) {
                ++queryMatch[k2 + l];
              }
							if (m->elem->alignedSegment[iAlignedSeg][l] == m->elem->alignedSegment[iAlignedSeg + 1][l]) { // if subject nucleotide matches a query nucleotide
								++ queryMatch[k2 + l];
							}
							else { // if subject nucleotide does not match a query nucleotide, store a mismatch in listMismatches for a query position
								++numMismatches;
								if (numMismatches >= maxListMismatches) {
                  /* resize listMismatches to new max size */
                  resizeListMismatches(&listMismatches, &maxListMismatches, numMismatches, alphabetSize);
								}
								// add mismatch at the end of the list
								listMismatches[numMismatches - 1].queryPos = k2 + l;
                int iNucl = strchr(strNucl, m->elem->alignedSegment[iAlignedSeg + 1][l]) - strNucl; // strchr returns pointer to char
								++ listMismatches[numMismatches - 1].mismatches[iNucl]; 
							  // iNucl: A --> 0, C --> 1, G --> 2, T --> 3
							}
						}
						else if (m->elem->alignedSegment[iAlignedSeg][l] == '-') { // if it is a gap, decrease index
							-- k2;
						}
					} // end for
          iAlignedSeg += 2; // next aligned segment between two exact matches
				} // end if a->next
				a = a->next;
			} // end while(a)
			m = m->next;
		} // while (m)
		
    // (2): compute total number of nucleotides in all msa = sum of all elements in queryMatch + sum of all elements in listMismatches
    // (3) compute matrix of observed freq-s from queryMatch and listMismatches;
    // (4) compute sum_of_all_observed_freq : sum of all elements in upper triangular matrix(including diagonal) of observed frequencies
    // (assuming that observedFreq(M, N) = observedFreq(N, M))
    totalNumNucleotides = getTotalNumNucleotidesObservedFreq(strandLen[i], alphabetSize, queryMatch,
      listMismatches, numMismatches, &observedFreq[0][0], seq, leftBorders[i], &sumObservedFreq, occNucl);

    // (5) compute expected freq-s
    //     expectedFreq(M, N) = 2 * occNucl[M] * occNucl[N] / (totalNumNucleotides * totalNumNucleotides)
    for (j = 0; j < alphabetSize; j++) {
			//expectedFreq[j][j] = (double)occNucl[j] / totalNumNucleotides * occNucl[j] / totalNumNucleotides;
			expectedFreq[i][j * alphabetSize + j] = (float)occNucl[j] / totalNumNucleotides * occNucl[j] / totalNumNucleotides;
      k = j + 1;
      for (; k < alphabetSize; k++) {
				//expectedFreq[j][k] += (double)2.0 * occNucl[j] / totalNumNucleotides * occNucl[k] / totalNumNucleotides;
				//    expectedFreq[k][j] = expectedFreq[j][k];
				expectedFreq[i][j * alphabetSize + k] += (float)2.0 * occNucl[j] / totalNumNucleotides * occNucl[k] / totalNumNucleotides;
				expectedFreq[i][k * alphabetSize + j] = expectedFreq[i][j * alphabetSize + k];
			}
    }
    // (6) compute gtr matrix : scoreMatrix(N1, N2) = rounded off value of(2 x log2(observedFreq(N1, N2) / expectedFreq(N1, N2)))
    for (j = 0; j < alphabetSize; j++) {
      double minScore = 0;
      for (k = j; k < alphabetSize; k++) {
        // if observedFreq[j][k] == 0, then set big penalty --> square of the maximal found penalty among all nucleotides
				//double score = ( observedFreq[j][k] > 0 ) ? (2 * log(observedFreq[j][k] / expectedFreq[j][k]) / log(2.0)) : 0; 
				double score = ( observedFreq[j][k] > 0 ) ? (2 * log(observedFreq[j][k] / expectedFreq[i][j * alphabetSize + k]) / log(2.0)) : 0; 
				if (score < 0) {
					score = (fabs(score) - fabs((int)score)) >= 0.5 ? (int)score - 1 : (int)score;
          minScore = score < minScore ? score : minScore;
        }
				else {
					score = (score - (int)score) >= 0.5 ? (int)score + 1 : (int)score;
				}
				(scoreMatrixAll)[i][k * alphabetSize + j] = (scoreMatrixAll)[i][j * alphabetSize + k] = (float)score;
      }
      for (k = j; k < alphabetSize; k++) { // update scores where observedFreq[j][k] == 0
        if ((scoreMatrixAll)[i][k * alphabetSize + j] == 0) {
          (scoreMatrixAll)[i][k * alphabetSize + j] = (scoreMatrixAll)[i][j * alphabetSize + k] = (float)(minScore + minScore); // set big penalty --> double the maximal found penalty 
        }
      }
			(*nuclFreq)[i][j] = (double)occNucl[j] / totalNumNucleotides;
    }
		
    free(queryMatch);
	} // end for i
	free(listMismatches);
}

/* allocate memory for scoring matrices and initialize their elem-s to 0 */
float **allocateScoreMatrixAll(Int64 numOfQueries, int alphabetSize) {
  Int64 i, j;
  float **scoreMatrixAll = NULL;
  // initialize matrices
  scoreMatrixAll = (float **)/*e*/malloc(sizeof(float *) * numOfQueries);
  for (i = 0; i < numOfQueries; i++) {
    scoreMatrixAll[i] = (float *)/*e*/malloc(sizeof(float) * alphabetSize * alphabetSize);
    for (j = 0; j < alphabetSize * alphabetSize; j++) {
      scoreMatrixAll[i][j] = 0;
    }
  }
  return scoreMatrixAll;
}

/* initialize listMismatches */
void initializeListMismatches(int maxListMismatches, int alphabetSize, Mismatch *listMismatches) {
	int i, j;
	for (i = 0; i < maxListMismatches; i++) {
		listMismatches[i].queryPos = -1;
		for (j = 0; j < alphabetSize; j++) {
			listMismatches[i].mismatches[j] = 0;
		}
	}
}

/* allocate and initialize listMismatches */
Mismatch *getListMismatches(int maxListMismatches, int alphabetSize) {
  Mismatch *listMismatches = (Mismatch *)/*e*/malloc(sizeof(Mismatch) * maxListMismatches);
	initializeListMismatches(maxListMismatches, alphabetSize, listMismatches);
  return listMismatches;
}

/* resize listMismatches to new max size */
void resizeListMismatches(Mismatch **listMismatches, int *maxListMismatches, int numMismatches, int alphabetSize) {
  int i, j;
  
  *maxListMismatches *= 2;
  *listMismatches = (Mismatch *)/*e*/realloc(*listMismatches, *maxListMismatches * sizeof(Mismatch));
  for (i = numMismatches; i < *maxListMismatches; i ++) {
    for (j = 0; j < alphabetSize; j ++) {
      (*listMismatches)[i].mismatches[j] = 0;
    }
  }
}

// compute 
// (a) total number of nucleotides in all msa = sum of all elements in queryMatch + sum of all elements in listMismatches (step 2)
// (b) observed frequencies (step 3)
// (c) sum_of_all_observed_freq : sum of all elements in upper triangular matrix(including diagonal) of observed frequencies;  observedFreq(M, N) = observedFreq(N, M)
long long getTotalNumNucleotidesObservedFreq(Int64 queryLen, int alphabetSize, short *queryMatch, 
	Mismatch *listMismatches, int numMismatches, 
  double *observedFreq, char *seq, Int64 leftBorderQ, long long *sumObservedFreq, Int64 *occNucl) {
  long long totalNumNucleotides = 0;
  Int64 j, k;
	*sumObservedFreq = 0;
  for (j = 0; j < queryLen; j++) { // column matches 
		if (queryMatch[j] && seq[leftBorderQ + j] != toupper(AMBIGUOUS_CHAR_QUERY)) {
      totalNumNucleotides += queryMatch[j];
      int qNucl = strchr(strNucl, seq[leftBorderQ + j]) - strNucl;
      observedFreq[qNucl * alphabetSize + qNucl] += (double)queryMatch[j] / 2.0 * (queryMatch[j] - 1);
      occNucl[qNucl] += queryMatch[j];
    }
  }

  for (j = 0; j < numMismatches; j++) { // column mismatches 
    Int64 qPos = listMismatches[j].queryPos; // query position
    char qNucl = strchr(strNucl, seq[leftBorderQ + qPos]) - strNucl; // index of query nucleotide (A-T --> 0-3)
    
    for (k = 0; k < alphabetSize; k++) {
      short numMismatchesK = listMismatches[j].mismatches[k];
      if (k != qNucl && numMismatchesK) {
        totalNumNucleotides += numMismatchesK;
        observedFreq[qNucl * alphabetSize + k] += (double)numMismatchesK * queryMatch[qPos];
        observedFreq[k * alphabetSize + qNucl] += (double)numMismatchesK * queryMatch[qPos];
        occNucl[k] += numMismatchesK;
      }
    } // for k
  } // for js
  
  // sum_of_all_observed_freq: sum of all elements in upper triangular matrix(including diagonal) of observed frequencies;  observedFreq(M, N) = observedFreq(N, M)
  for (j = 0; j < alphabetSize; j++) {
    for (k = j; k < alphabetSize; k++) {
      *sumObservedFreq += (long long)observedFreq[j * alphabetSize + k]; // != totalNumNucleotides
    }
  }
	for (j = 0; j < alphabetSize; j++) {
		for (k = 0; k < alphabetSize; k++) {
			observedFreq[j * alphabetSize + k] /= *sumObservedFreq; 
		}
	}
	return totalNumNucleotides;
}

/* query specific: calculate lambda --> used for computing E-value; complexity adjused scoring 
 * Reference: BLAST Scoring Parameters, E. Michael Gertz, 16 March 2005. 
 * http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/algo/blast/core/blast_stat.c
 */
double getLambda(float *scoreMatrixSingleQuery, int alphabetSize, double *nuclFreqs) {
	int i, j;
	double sum, check;
	double lambdaSingleQuery = 0.5, lambdaUpper = 0, lambdaLower = 0;
	do {
		sum = 0;
		check = 0;
		for (i = 0; i < alphabetSize; i++) {
			for (j = 0; j < alphabetSize; j++) {
				if (nuclFreqs[i] && nuclFreqs[j]) {
          sum += nuclFreqs[i] * nuclFreqs[j] * exp(lambdaSingleQuery * scoreMatrixSingleQuery[i * alphabetSize + j]);
					check += nuclFreqs[i] * nuclFreqs[j];
				}
			}
		}
		assert((check < (double)1.001) && (check >(double)0.999));
		if (sum < 1.0) {
      lambdaLower = lambdaSingleQuery;
      lambdaSingleQuery *= 2.0;
		}
	} while (sum < 1.0);

  lambdaUpper = lambdaSingleQuery;
	while (lambdaUpper - lambdaLower > (double)0.00001) {
    lambdaSingleQuery = (lambdaLower + lambdaUpper) / 2.0;
		sum = 0;
		check = 0;
		for (i = 0; i < alphabetSize; i++) {
			for (j = 0; j < alphabetSize; j++) {
				if (nuclFreqs[i] && nuclFreqs[j]) {
          sum += nuclFreqs[i] * nuclFreqs[j] * exp(lambdaSingleQuery * scoreMatrixSingleQuery[i * alphabetSize + j]);
					check += nuclFreqs[i] * nuclFreqs[j];
				}
			}
		}
		assert((check < (double)1.001) && (check >(double).999));
		if (sum >= 1.0) {
      lambdaUpper = lambdaSingleQuery;
		}
		else {
      lambdaLower = lambdaSingleQuery;
		}
	}
  return lambdaSingleQuery;
}

/* free memory - nuclFreq */
void freeNuclFreq(double **nuclFreq, int numOfQueries) {
	int i;
	for (i = 0; i < numOfQueries; i++) {
		free(nuclFreq[i]);
	}
	free(nuclFreq);
}

/* free memory - scoreMatrix */
void freeScoringMatrices(float **scoreMatrixAll, int numOfQueries) {
	int i;
	for (i = 0; i < numOfQueries; i++) {
    free(scoreMatrixAll[i]);
	}
  free(scoreMatrixAll);
}

