/***** mergeList.h *********************************************************************************
 * Description: Header file for functions for merging elements of linked lists of aligned regions.
 * Author: Mirjana Domazet-Loso
 * 
 * This file is part of gmos.
 *
 ***************************************************************************************************/

#ifndef MERGELIST_H
#define MERGELIST_H

#define MAXLEN 10000
#define MAXGAP 2000

#define MAXSEQNAME 1000
#define DIST 100L // block size
#define DIST1 0L
#define THRESHOLD_McNEMAR 3.84 // p < 0.05 (If the test statistic is > 3.84, the p-value will be < 0.05)
#define RATIO_DIFF 0.05f
#define GAP_PENALTY 1

typedef struct mnode {
	//struct aregion *elem;			/* pointer to an element of a aregion list */  
	aRegion *elem;				      /* pointer to an element of a aregion list */  
  Int64 subjectId;            /* subject id */
	struct mnode *next;         /* pointer to the next element of a merged list of pointers */
} mNode;

//extern double EValueLow;
//extern double lambda; // query-specific 
//extern float ambCharPenalty;
//extern float scoreMatrix[ALPHABET_DNA_SIZE][ALPHABET_DNA_SIZE];

mNode *getMNode(aRegion *p, Int64 subjectId);
void freeMNode(mNode *m);

/* construct |Q| lists of pointers to first elements (heads) of each list; list must be sorted in asc order of lb of first elem-s */
void getCandidateList(aRegion **head, Int64 numOfQueries, Int64 numOfSubjects, mNode ***candidateList, short argB);

void printRow(FILE *fpout, aRegion *p, Int64 queryStart, Int64 queryEnd, char queryStrand, char **seqHeaders, Int64 *strandLen, Int64 subjectId, Int64 numOfQueries);

void printCandidateList(mNode **mList, FILE *fpout, Int64 numOfQueries, char *header, char **seqHeaders, short argB, Int64 *strandLen);

/* moving elem-s from candidateList to finalList */
void addCandidateElem(mNode **candidateListI, mNode *p);

void constructFinalList(mNode ***candidateList, aRegion **head, Int64 numOfQueries, Int64 numOfSubjects, mNode ***finalList, short argB,  
      char **seqHeaders, Int64 *strandBorders, Int64 *leftBorders, FILE *fqout, long long allSeqLength, Int64 *strandLen, 
			double *lambdaAll);

void freeMList(mNode **mList, Int64 numOfQu, Int64 numOfSubjects);

/* extend segment up to n (n's begin part - before the first match stays as it is) */
void extendMToN(mNode *m, mNode *n, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned,
	Int64 numOfQueries, Int64 i, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength,
	Int64 nStartPosQ, char fastScore);
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty);

/* extend segment n towards m, which is its predecessor (m's end part, after the last match, stays as it is) */
void extendNToM(mNode *n, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned,
	Int64 numOfQueries, Int64 Qi, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength,
	Int64 mEndPosQ, char fastScore);
	//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty);

void deleteAlignmentBeforeStart(mNode *n);
void deleteAlignmentAfterEnd(mNode *n);

/* delete alignedSegments starting from the z-th element */
void deleteAlignedSegments(mNode *n, int z);
/* delete aligned segments before the position z */
void deleteAlignedSegmentsBeforeZ(mNode *n, int z);

/* checkDeleteLast4 returns 1 (keep new, delete last), -1 (keep last and delete new) or 0 (keep both last and new) */
int checkDeleteLast4(mNode *n, mNode *last, long long minTotalLen);

/* extend segment m over n */
void extendMToNOverlap(mNode **m, mNode **n, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned,
	Int64 numOfQueries, Int64 i, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength,
	mNode **mprev, char fastScore);
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty);

/* extend segment n over m */
void extendNToMOverlap(mNode **m, mNode **n, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned,
	Int64 numOfQueries, Int64 i, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength,
	mNode **mprev, mNode ***finalList, char fastScore);
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty);

void splitMNMiddle(mNode **m, mNode **n, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned,
	Int64 numOfQueries, Args *args, int maxCols, char *seq, element *mat,
	long long allSeqLength, mNode **mprev, Int64 iQuery, mNode ***finalList, char fastScore);
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty);

void connectFinalList(mNode ***finalList, Int64 numOfQueries, Args *args, char *seq, Int64 *leftBorders, Int64 *strandBorders,
                      char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat, 
											long long allSeqLength, double *lambdaAll, float **scoreMatrixAll, 
											int alphabetSize, float ambCharPenaltyAll, char fastScore, double *arrayK);

/* function returns 1 if m is "better" than n, -1 if if n is "better" than m, and 0 if bot are equal when ratios and e-scores are compared */
int getBetterOverallScore(mNode *m, mNode *n);

void concatenateFinalList(mNode **finalList, Int64 numOfQueries, Args *args, char *seq, Int64 *leftBorders, Int64 *strandBorders,
                      char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat, 
											long long allSeqLength, Int64 numOfSubjects, double *lambdaAll, float **scoreMatrixAll, 
											int alphabetSize, float ambCharPenaltyAll, char fastScore, double *arrayK);

// count number of nucleotides in seq; gaps not included
int getNumNucleotides(char *seq);

//void alignRecSmSn(FILE *fpout, char *seqRec1, char *seqS1, char *seqRec2, char *seqS2, Int64 *cntAll1, Int64 *cntAll2, Args *a, Int64 *bpStart, Int64 *bpEnd);
void alignRecSmSn1(char *seqRec1, char *seqS1, char *seqRec2, char *seqS2, Int64 *cntAll1, Int64 *cntAll2, Args *a, Int64 *bpStart, Int64 *bpEnd);

void findBP(mNode **finalList, Int64 numOfQueries, Args *args, char *seq, Int64 *leftBorders, Int64 *strandBorders,
	char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat,
	long long allSeqLength, Int64 numOfSubjects, int maxGap,
	double *lambdaAll, float **scoreMatrixAll, int alphabetSize, float ambCharPenaltyAll, char fastScore, 
	Int64 *strandLen, double *arrayK);

/* return 1 if both start and end for query and subject are on the same strand */
int checkBorders(Int64 startQ, Int64 endQ, Int64 startS, Int64 endS, Int64 *strandLen, Int64 i, Int64 j);

// copy n char-s from src to dest; count non-gap characters
int strncpyNoGapCnt(char *dest, char *src, Int64 n);

/* compute alignment score of two aligned sequences */
double getAlignmentScore(char *seq1, char *seq2, Args *a);

double getMcNemar(double b, double c);
double getMcNemarNormalized(double b, double c, double r);

Int64 strncpyStart(char *dest, char *src, Int64 start);

/* print a single row in alignment */
void printAlignmentRow(aRegion *r, int *k, int *lenTemp, char **tempSeq, Int64 iSeq, Int64 a_sl, Int64 lb, Int64 *leftBorders, char *seq, FILE *fout);
/* print alignment in multi fasta format */
void printAlignment(mNode **mList, FILE *fout, Int64 numOfQueries, char **seqHeaders, Int64 *strandLen, char *seq, Int64 *leftBorders);

void refineAlignment(mNode **finalList, Int64 numOfQueries, char *seq, Int64 *leftBorders, Int64 *strandBorders,
	Int64 *strandLen, float **scoreMatrixAll, int alphabetSize, float ambCharPenaltyAll, double *lambdaAll, int maxCols,
	char *first, char *second, char *firstAligned, char *secondAligned, element *mat, Args *arguments, 
	long long allSeqLength, char fastScore, double *arrayK);

// coount gaps in string
int cntGaps(char *str);

void updateBP(mNode **finalList, Int64 numOfQueries, Args *args, Int64 *leftBorders, char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat,
	long long allSeqLength, Int64 numOfSubjects, char *seq, char fastScore, Int64 *strandBorders, Int64 *seqBorders, double *arrayK);

#endif // MERGELIST_H
