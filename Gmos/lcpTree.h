/***** lcpTree.h *************************************************************
 * Description: Functions for lcp-interval tree processing.
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
 *
 *****************************************************************************/

#ifndef LCPTREE_H
#define LCPTREE_H

#define MAXNUMLEAVES 9 /* size of the alphabet: A,C,T,G,Z,$,x,y --> x/y are added as masks to all non-A,C,G,T characters in query/subject input; Z is the border character */
#define MAXNUMLEAVESALL 1000

//#define BSEARCH 1 /* binary search is set */
#define MTREES 1 /* m binary trees --> for each of m subjects, construct a binary tree of intervals */
#define STEP 10

double lambda; // query-specific 
double blastK; // query specific
float ambCharPenalty;
float scoreMatrix[ALPHABET_DNA_SIZE][ALPHABET_DNA_SIZE];
double EValueLow;

/* main function in this file (the entry point) */
void getLcpTreeShulens(FILE *fpout, Args *a, SequenceUnion *seqUnion, FILE *fqout, FILE *fQout, FILE *faout, FILE *fLen);

// traverse lcp tree and compute all intervals for all queries 
void traverseLcpTree(Int64 *lcpTab, Int64 *sa, Sequence *seq, Int64 numOfSubjects, Int64 numOfQueries, Int64 *seqBorders
											, Int64 *leftBorders, Int64 *strandBorders
											, Int64 maxDepth, Int64 *maxShulens
											, qNode ***root, Args *a);

//void updateShulen(Int64 **shulens, Int64 intervalLcp, Int64 suffixArrayValue, Int64 j, Sequence *query, short iQuery);
void checkLeaves(Interval *interval, Int64 *seqBorders, Int64 numOfSubjects);

void addQueryLeaf(Int64 *queryLeaves, int maxCols, Int64 *numQueryLeaves, Int64 queryIndex, Int64 pos);

void checkLeaves2 (Interval *interval, Int64 *seqBorders, Int64 lb, Int64 rb, short *QS, Int64 *queryLeaves, int maxCols,
										Int64 *numQueryLeaves, Int64 *leftBorders, Int64 *sa,
                    Int64 **subjectLeaves, Int64 *numSLeaves, Int64 *maxSLeaves, Int64 numOfQueries);

void determineQS(Interval *interval, Int64 numOfSubjects, Int64 *seqBorders, Int64 *leftBorders, short *subjects, Int64 *queryLeaves, 
								 //int maxCols, Int64 *numQueryLeaves, Int64 *sa, Int64 *subjectList, Int64 *numOfSubjInterval);
								 int maxCols, Int64 *numQueryLeaves, Int64 *sa,
                 Int64 **subjectLeaves, Int64 *numSLeaves, Int64 *maxSLeaves, Int64 numOfQueries);

void reallocUnresLeaves(Int64 **unresolvedQLeaves, Int64 *maxUnresQLeaves, Int64 numUnresolvedQLeaves, int multiplyBy);

void process2(Interval *interval, Int64 *seqBorders, 
							Int64 *leftBorders, Int64 *strandBorders, Int64 numOfSubjects, short *QS, Int64 *queryLeaves, int maxCols, Int64 *sa, 
							Int64 **unresolvedQLeaves, Int64 *maxUnresQLeaves, Stack *reserveQIStack, Int64 *maxShulens, 
							Int64 **subjectLeaves, Int64 *maxSLeaves, Int64 *numSLeaves, Int64 numOfQueries, qNode ***root,
							Int64 **subjectLeavesAll, Int64 *maxSLeavesAll, Int64 *numSLeavesAll);

void resolveQueryIntervals(Int64 numQueryLeaves, Int64 *queryLeaves, int maxCols, Int64 *seqBorders, Int64 *leftBorders, Int64 *strandBorders,
													 Int64 numOfSubjects, Interval *interval, Stack *reserveQIStack, Int64 *maxShulens, qNode ***root,
													 Int64 **subjectLeaves, Int64 *numSLeaves/*, Int64 numOfQueries*/);


/* auxiliary function to print out sa and lcp */
void printSA_LCP(FILE *fpout, Int64 *sa, Int64 *lcp, Int64 seqUnionLen);

/* find a subject/query for each position in the suffix array */
short *getQS(Int64 *seqBorders, Int64 *leftBorders, Int64 numOfSubjects, Int64 numOfQueries, int step);
void initSubjectLeaves(Int64 ***subjectLeaves, Int64 **maxSLeaves, Int64 **numSLeaves, Int64 maxNumLeaves, Int64 numOfSubjects);

void freeSubjectLeaves(Int64 **subjectLeaves, Int64 numOfSubjects);

void printSeqLen(FILE *fLen, Int64 *strandLen, Int64 numOfQueries, Int64 numOfSubjects);

void addQueryInterval(Int64 i, Int64 query_lb, Int64 query_rb, Int64 query_sl, Int64 k, Int64 numOfSubjects,
	Stack *reserveQIStack, Int64 *maxShulens, qNode ***root, Int64 **subjectLeaves, Int64 *numSLeaves, 
	Int64 *strandBorders, Int64 *leftBorders, Int64 numOfQueries);

void checkQueryInterval(Int64 i, Int64 *query_lb, Int64 *query_rb, Int64 *query_sl, Int64 *leftBorders,
	Int64 *strandBorders, Int64 *seqBorders);

int existNegValues(float *scoreMat, int alphabetSize);

#endif // LCPTREE_H

