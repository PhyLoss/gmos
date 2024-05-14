/***** resolveOverlaps.c *************************************************************
* Description: Source file for resolving overlaps between aligned segments.
*
* Author: Mirjana Domazet-Loso
*
* This file is part of gmos.
*****************************************************************************/

#ifndef RESOLVEOVERLAPS_H
#define RESOLVEOVERLAPS_H

#define RATIO_MINSEGLENGTH 0.8

//extern double lambda; // query-specific 
//extern float ambCharPenalty;
//extern float scoreMatrix[ALPHABET_DNA_SIZE][ALPHABET_DNA_SIZE];
double gapPenalty;
double excellentEValue; // e.g. EValueLow * EValueLow
double EValueLow;

/* compute e-value */
double getEValueScoreBit(double score, double queryLen, long long allSeqLength);

/* compute score of a segment within m starting at startSeg and ending at endSeg */
double getScoreSeg(mNode *m, char *seq, Int64 iQuery, Int64 *leftBorders, Int64 startSeg, Int64 endSeg, double gapPenaltyLocal);

//void freeARegionFromMToN(aRegion *r, aNode *m, aNode *n, long long allSeqLength, aRegion *r2/*, double lambda*/);
void freeARegionFromMToN(aRegion *r, aNode *m, aNode *n, long long allSeqLength, aRegion *r2/*, double lambda*/, char *seq, Int64 *leftBorders, Int64 iQuery);

/* returns 1 if e-values are equal or equally good, and otherwise 0 */
int compareEValues(double eValueM, double eValueN);

/* find new leftMost and its predecessor, such that leftMost includes startPosQ */
void findLeftMost(mNode **leftMost, mNode **prevLM, Int64 startPosQ);

/* find new rightMost node up to n */
void findRightMost(mNode *leftMost, mNode *n, mNode **rightMost);

void connectPrev(int leaveM, mNode *m, mNode *n, mNode **prev, mNode ***finalList, int i);

int shortenM(mNode **m, mNode *n, Int64 bpStart, Int64 bpEnd, mNode ***finalList, Int64 iQuery, Args *a,
  long long allSeqLength, mNode **m2, int leaveLeftPart/*, double lambda*/, char *seq, Int64 *leftBorders);

void shortenN(mNode **nn, Int64 start, Int64 end, long long allSeqLength, Args *a, Int64 *leftBorders, Int64 iQuery /*, double lambda*/, char *seq);

void resolveOverlapsFinalList(FILE *fpout, mNode ***finalList, Int64 numOfQueries, Args *args, char *seq, Int64 *leftBorders, Int64 *strandBorders,
  char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat,
  long long allSeqLength, Int64 numOfSubjects, double *lambdaAll, float **scoreMatrixAll, int alphabetSize, 
	float ambCharPenaltyAll, char fastScore, double *arrayK);

void connectM2(mNode *tprev, mNode **m2);

void connectNNext(mNode **m, mNode *n);

int roundNeg(double value);

#endif // RESOLVEOVERLAPS_H