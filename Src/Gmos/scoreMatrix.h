/***** matrix.h *************************************************************
* Description: Header file for dna scoring matrix (GTR model).
*
* Author: Mirjana Domazet-Loso
*
* This file is part of gmos.
*****************************************************************************/

#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

typedef struct mismatch {
	Int64 queryPos;                       /* query position */
	short mismatches[ALPHABET_DNA_SIZE];  /* number of mismatches for each nucleotide A,C,G,T over all subjects aligned to queryPos column*/
} Mismatch;

/* calculate scoreMatrix and nuclFreq for each query */
void getScoringMatrices(mNode **finalList, Int64 numOfQueries, char *seq, Int64 *leftBorders, Int64 *strandBorders,
	Int64 *strandLen, float **scoreMatrix, int alphabetSize, double ***nuclFreq, float **expectedFreq);

/* allocate memory for scoring matrices and initialize their elem-s to 0 */
float **allocateScoreMatrixAll(Int64 numOfQueries, int alphabetSize);

void initializeListMismatches(int maxListMismatches, int alphabetSize, Mismatch *listMismatches);

/* allocate and initialize listMismatches */
Mismatch *getListMismatches(int maxListMismatches, int alphabetSize);

/* resize listMismatches to new max size */
void resizeListMismatches(Mismatch **listMismatches, int *maxListMismatches, int numMismatches, int alphabetSize);

// compute (a) total number of nucleotides in all msa and (b) observed frequencies
long long getTotalNumNucleotidesObservedFreq(Int64 queryLen, int alphabetSize, short *queryMatch, 
	Mismatch *listMismatches, int numMismatches,
	double *observedFreq, char *seq, Int64 leftBorderQ, long long *sumObservedFreq, Int64 *occNucl);

double getLambda(float *scoreMatrixSingleQuery, int alphabetSize, double *nuclFreqs);

void freeNuclFreq(double **nuclFreq, int numOfQueries);
void freeScoringMatrices(float **scoreMatrixAll, int numOfQueries);

#endif // SCOREMATRIX_H
