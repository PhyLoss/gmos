/***** align.h *************************************************************
 * Description: Header file for alignment functions.
 * Author: Mirjana Domazet-Loso
 * 
 * This file is part of gmos.
 *****************************************************************************/

#ifndef ALIGN_H
#define ALIGN_H
  
typedef struct anode {
  Int64 sl;                         /* shulen attached to the left-most position of this interval */
  Int64 lbQ;                        /* left border of a query segment */
  Int64 lbS;                        /* left border of a subject segment */

	struct anode *next;							  /* next element in the list of an aligned region */
} aNode;

typedef struct aregion {
  aNode *start;                     /* pointer to the first anode of an aligned region */
  aNode *end;                       /* pointer to the last anode of an aligned region */
  double score;                      /* total score over all aligned segments in a region */

  char **alignedSegment;            /* list of strings representing segments between anodes which are aligned using NW algorithm; there are 2 * numAlignedSegment
                                     * strings; *2 since a pair of segments is aligned: a segment of Q and a segment of S */
  double *scoreAlignedSegment;      /* array of scores of aligned segments */
  int numAlignedSegment;            /* number of aligned segments = number of anodes in aregion - 1*/
  int maxAlignedSegment;            /* currently maximal elements in scoreAlignedSegment */

////////////// NEW --> aligned segments just before/after aRegion
  char **alignedStart; /* aligned segment just before aRegion --> Q: alignedStart[0], S: alignedStart[1] */
  char **alignedEnd;   /* aligned segment just after aRegion --> Q: alignedEnd[0], S: alignedEnd[1] */
  double scoreStart, scoreEnd;
////////////// end - NEW

  ////////////// NEW --> final start and end position of a region - when aligned segments just before/after start-anode and end-anode are included
  Int64 startPosQ;
  Int64 endPosQ;
  Int64 startPosS;
  Int64 endPosS;

  double ratio;                     /* ratio of total score over total length of the region */
  long long totalLen;                   /* total length of the region (including gaps) */
  double score_blast;               /* blast score (bits) */
  double e_value;                   /* e-value blast */
  ////////////// end - NEW
  char active;                      /* region is still open, i.e. new aNodes can be added to its end (1), or closed (0) */
  struct aregion *next;             /* next region */  
} aRegion;

#define LEN_OVERLAP 20 // ?????????????????????????

//extern double K; // query-specific 
//extern double lambda; // query-specific 
//extern float ambCharPenalty;
//extern float scoreMatrix[ALPHABET_DNA_SIZE][ALPHABET_DNA_SIZE];

Int64 absInt64(Int64 a);

aNode *getANode(lNode *p);
aRegion *getARegion(aNode *a, float match);

// testing purpose
//aNode *getANodeTest(lNode *p, int iSubject);

/* pointers to first elements of |Q| * |S| aregions' lists */
aRegion **getARegionLists(Int64 numOfQueries, Int64 numOfSubjects);

void checkMaxAlignedSegments(aRegion *r);

/* pointers to first elements of |Q| * |S| lists of aligned regions */
aRegion **getARHead(Int64 numOfQueries, Int64 numOfSubjects);

/* free aRegion*/
void freeARegion(aRegion *r);
void _freeARegionFromMToN(aRegion *r, aNode *m, aNode *n, long long allSeqLength, aRegion *r2/*, double lambda*/);

/* free all aregion lists */
void freeARList(aRegion **arList, Int64 numOfQueries/*, Int64 numOfSubjects*/);

/* close list of segments in region r and compute all alignments between anode-s in a region using Needleman-Wunsch alg. */
void computeARegion(aRegion *r, Args *arguments, char *seq, element *mat, char *firstAligned, char *secondAligned, int maxCols,
	int iQuery, int iSubject, Int64 *leftBorders, char *first, char *second, Int64 *strandBorders,
	Int64 numOfQueries, aRegion *tail, long long allSeqLength, char fastScore/*, double lambda, float *scoreMatrix,
	int alphabetSize, float ambCharPenalty*/);


/* check whether aRegion is shorter than set by -f option, or has score less than set by -F option; in that case return 0, else return 1 */
int checkARegion (aRegion *r, Args *arguments);

/* delete last region from aRegion list, and redirect pointers */
void deleteLastARegion(aRegion **ahead, aRegion **tail, aRegion **prev, aRegion *r);

/* find minimal distance between any subject position of p and any subject position of q */
void findMinDistance(lNode *p, lNode *q, Args *arguments, Int64 *minDistance, int iSubject, Int64 *leftBorders, Int64 *strandBorders, aRegion *r);

/* construct all lists of aligned regions for each combination (Qi, Sj) */
void constructAlignment(char *seq, lNode **head, Args *arguments, aRegion ***ahead, Int64 numOfQueries, Int64 numOfSubjects,
	Int64 *leftBorders, Int64 *strandBorders, Int64 *seqBorders, long long allSeqLength,
	double *lambdaAll, float **scoreMatrixAll, int alphabetSize, float ambCharPenalty, char fastScore, double *arrayK);

/* print complete aregion list */
void printARegionList(aRegion **arList, Int64 numOfQueries, Int64 numOfSubjects, FILE *fpout, short argB, char **seqHeaders, Int64 *strandLen);

/* get beginning position of an aligned segment just before r */
Int64 getBeginPos(Int64 beginRPos, Int64 tailEndPos, aRegion *tail, Int64 afterTailLen, Int64 argB, Int64 *leftBorders, Int64 *strandBorders, int i, aRegion *r);

/* get end position of an aligned segment just after r */
Int64 getEndPos(Int64 endRPos, Int64 argB, Int64 *leftBorders, Int64 *strandBorders, int i);

/* align segments before the beginning of r and after the end of r */
void addAlignmentBeginEnd(aRegion *r, aRegion *tail, Args *arguments, char *seq, element *mat, char *first, char *second, int maxCols,
	int iQuery, int iSubject, Int64 numOfQueries, Int64 *leftBorders, Int64 *strandBorders,
	char *firstAligned, char *secondAligned, char fastScore);
	//float *scoreMatrix, int alphabetSize, float ambCharPenalty);

// connect r2 to tail->end, re-compute alignment between tail->end and r2->start, and then delete r2
void connectARegion(aRegion *tail, aRegion *r2, Args *arguments, char *seq, Int64 *leftBorders, int iSubject, int iQuery, Int64 numOfQueries,
	char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat,
	long long allSeqLength/*, double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty*/, char fastScore);

void addAlignedSegments(aRegion *tail, aRegion *r2);

/* free aRegion partly - when part of its items is re-connected to another aRegion */
void freeARegionPartly(aRegion *r);

/* allocate space for the variables needed for the computation of nw alignment */
void allocateElemsNW(element **mat, char **firstAligned, char **secondAligned, int *maxCols, char **first, char **second, Args *arguments, Int64 extraLen);

/* free elements needed for the computation of nw alignments */
void freeElemsNW(element *mat, char *firstAligned, char *secondAligned, char *first, char *second);

/* align a segment between region from the start position in the i-th sequence (length m) to the start position of the j-th sequence (length n) */
void alignSegment(Int64 startPosISeq, Int64 startPosJSeq, Int64 *leftBorders, Int64 m, Int64 n, char *first, char *second, char *firstAligned, char *secondAligned, char *seq
	, Int64 iSeq, Int64 jSeq, element *mat, Args *arguments, int maxCols, char fastScore);

/* copy aligned segments to aRegion's elements alignedSeg1 and alignedSeg2 and update appropriate scores */
void storeAlignedSegments(double *newScore, double score, double *totalScore, char *firstAligned, char *secondAligned, char **alignedSeg1, char **alignedSeg2);

//void computeEValueScoreBit(aRegion *r, long long allSeqLength);
void computeEValueScoreBit(aRegion *r, long long allSeqLength/*, double lambda*/);

int connectARegion2(aRegion *tail, aRegion *r2, Args *arguments, char *seq, Int64 *leftBorders, int iSubject, int iQuery, Int64 numOfQueries,
	char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat,
	long long allSeqLength, char fastScore);

void storeAlignedSegments2(double *newScore, double score, double *totalScore, char *firstAligned, char *secondAligned, char **alignedSeg1, char **alignedSeg2);

void aNodeToM(aRegion *p, aRegion *p2, aNode *nn);

// new functions
aNode *getANodeNextPosition(lNode *p, int i);
void deleteCurrARegion(aRegion **ahead, aRegion **tail, aRegion **prev, aRegion *r/*, aRegion **activeLast, aRegion **activeFirst*/);

/* construct a single list of aligned regions for (Qi, Sj) */
void constructARegionList2(lNode *head, aRegion **ahead, Args *arguments, char *seq, element *mat,
  char *firstAligned, char *secondAligned, int maxCols, int iQuery, int iSubject, Int64 *leftBorders,
  char *first, char *second, Int64 *strandBorders, Int64 *seqBorders, Int64 numOfQueries, long long allSeqLength, char fastScore);

void initializeHeadTail(lNode *p, aRegion **ahead, aRegion **tail, aRegion **activeFirst, aRegion **prev, aRegion **prevActiveFirst, Args *arguments);
void addANode(aRegion *rCurr, aNode *a, Args *arguments);
void deleteANode(aRegion *rCurr, Args *arguments);

/* compute start and end coordinate on forward strand */
void getSubjectStartEndFR(Int64 *subjectStart, Int64 *subjectEnd, Int64 *strandLen, Int64 subjectId, Int64 numOfQueries, aRegion *p, char strainFR);

/* set query specific lambda and scoreMatrix */
void setQueryLambdaMatrix(double *lambdaAll, float **scoreMatrixAll, int alphabetSize, float ambCharPenaltyAll, Int64 i, double *arrayK);

/* calculate gap penalty */
float getGapPenalty();

/* calculate score between two exact matches using scoreMatrix*/
double scoreMatch(Int64 start, Int64 end, Int64 iQuery, char *seq, Int64 *leftBorders);

/* calculate score between two aligned regions using scoreMatrix*/
double scoreAlignedRegion(Int64 len, char *qSeq, char *sSeq, double gapPenaltyLocal);

#endif // ALIGN_H
