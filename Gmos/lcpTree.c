/***** lcpTree.c *************************************************************
 * Description: Functions for lcp-interval tree processing.
 * Reference: Abouelhoda, M. I., Kurtz, S., and Ohlebusch, E. (2002).
 *   The enhanced suffix array and its applications to genome analysis.
 *   Proceedings of the Second Workshop on Algorithms in Bioinformatics,
 *   Springer-Verlag, Lecture Notes in Computer Science.
 *
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "commonSC.h"
#include "eprintf.h"
#include "stringUtil.h"
#include "expectedShulen.h"
#include "sequenceData.h"
#include "sequenceUnion.h"
#include "intervalStack.h"
#include "interval.h"

#include "shulen.h"
#if VER32
#include "divsufsort.h"
#endif
#include "interface.h"
#include "queryBTNode.h"
#include "lNode.h"
#include "nw.h"
#include "align.h"
#include "mergeList.h"
#include "scoreMatrix.h"
#include "resolveOverlaps.h"
#include "blast_stat.h"
#include "lcpTree.h"


#if defined(_DEBUG) && defined(WIN) 
	#include "leakWatcher.h"
#endif

#if defined(_DEBUG) && defined(WIN)  
  #define new DEBUG_NEW
  #undef THIS_FILE
  static char THIS_FILE[] = __FILE__;
#endif

double lambda = DEFAULT_LAMBDA; // query-specific 
float ambCharPenalty = DEFAULT_G;
float scoreMatrix[ALPHABET_DNA_SIZE][ALPHABET_DNA_SIZE] = { { 0.0 } };
double EValueLow;

/* getLcpTreeShulens: compute intervals for each query; 
 * This is the only entry point to the functions in this file.
 */
void getLcpTreeShulens(FILE *fpout, Args *a, SequenceUnion *seqUnion, FILE *fqout, FILE *fQout, FILE *faout, FILE *fLen) {

	Int64 *sa = NULL, *lcpTab = NULL;
	Int64 i, j, k, ns;
	Int64 maxDepth;
	Int64 *leftBorders = NULL, lb;
	Int64 *strandBorders = NULL;
  Int64 *strandLen = NULL; // new
	Int64 *maxShulens = NULL, maxs = 0, lS0 = 0;
	Int64 *minSumWin = NULL; // minimal sum (threshold) for which winner-sequences are considered to have strong signal

  qNode **root = NULL; // there are numOfQueries * numOfSubjects interval trees
  lNode **head = NULL; // there are numOfQueries * numOfSubjects linked lists construced from binary interval trees
  mNode **candidateList = NULL; // list of candidates to be added to a merged list; initially: first elements of all subject lists for a single query
  mNode **finalList = NULL; // final list of - merged list of "winning" nodes among subject linked lists

  aRegion **aRegionList = NULL; // lists of aligned regions; for each combination (Qi, Sj)

	float **scoreMatrixAll = NULL; // numOfQueries x (ALPHABET_DNA_SIZE x ALPHABET_DNA_SIZE) 
	float **expectedFreq = NULL; // numOfQueries x (ALPHABET_DNA_SIZE x ALPHABET_DNA_SIZE) 
	int alphabetSize = ALPHABET_DNA_SIZE;  
	double **nuclFreq = NULL; // numOfQueries x ALPHABET_DNA_SIZE
	double *arrayLambda = NULL; // lambda for each query
	double *arrayK = NULL; // K for each query
	float ambCharPenaltyAll;
	char fastScore = 1;
#if DEBUG
	clock_t end, start, end2, end3;
	double elapsed_time1, elapsed_time2, elapsed_time3;
#endif	
  element *mat = NULL;
  char *firstAligned = NULL;
  char *secondAligned = NULL;
  char *first = NULL;
  char *second = NULL;
  int maxCols;
  int maxGap = 0;
  long long allSeqLength = 0;

  maxDepth = a->D;
  EValueLow = a->E;
	// array of left borders of each sequence
	leftBorders = (Int64 *)/*e*/malloc(sizeof(Int64) * (seqUnion->numOfSubjects + seqUnion->numOfQueries)); 
	// array of fwd strand borders of each sequence
	strandBorders = (Int64 *)/*e*/malloc(sizeof(Int64) * (seqUnion->numOfSubjects + seqUnion->numOfQueries));
	strandLen = (Int64 *)/*e*/malloc(sizeof(Int64) * (seqUnion->numOfSubjects + seqUnion->numOfQueries));
	lb = 0;
	for (i = 0; i < seqUnion->numOfSubjects + seqUnion->numOfQueries; i++) {
		leftBorders[i] = lb;
		strandBorders[i] = leftBorders[i] + (seqUnion->seqBorders[i] - leftBorders[i]) / 2;
		strandLen[i] = strandBorders[i] - leftBorders[i] + 1;
		lb = seqUnion->seqBorders[i] + 1;
	}

	/* for each query form an array of pointers; each pointer points to the query interval 
	* whose right border is closest to the upper bound in terms of args->q, 
	* e.g. when qi.rb = 978, then an element [qi][0] points to qi, that is 978 is closest to 999=upper bound for [qi][0]
	*/

	// compute suffix array
#if DEBUG  	
	start = clock();
#endif
#ifdef VER32
  sa = getSuffixArray2(seqUnion->seqUnion);
#else
  sa = getSuffixArray(seqUnion->seqUnion);
#endif
#if DEBUG  	
	end = clock();
	elapsed_time1 = (double)(end - start) / CLOCKS_PER_SEC;
#endif	
	if (!sa) {
		eprintf("[ERROR]Suffix array: out of memory!\n");
	}	
	
	// compute lcp array
	lcpTab = getLcp(seqUnion->seqUnion, sa);
	if (!lcpTab) {
		eprintf("[ERROR]LCP: out of memory!\n");
	}
	lcpTab[0] = 0; /* or -1 */
  lcpTab[1] = 0; /* since data in sa and lcpTab start with position 1, and not 0 */
#if DEBUG  	
	end2 = clock();
	elapsed_time2 = (double)(end2 - end) / CLOCKS_PER_SEC;
#endif	
  // print sa, lcp 
  #if DEBUG
	  printSA_LCP(fpout, sa, lcpTab, seqUnion->len);
  #endif

	// print run-time
#if DEBUG
	if (a->t) {
		printf( "\nSA calculation: %.2f seconds.\n", elapsed_time1);
		printf( "\nLCP calculation: %.2f seconds.\n", elapsed_time2);
	}
#endif
  /* calculate max shulens expected only by chance for each query */
	/* using both subject's and query's gc-content */
	maxShulens = (Int64 *)/*e*/malloc(seqUnion->numOfQueries * sizeof(Int64));
	minSumWin = (Int64 *)/*e*/malloc(seqUnion->numOfQueries * sizeof(Int64));
	lS0 = seqUnion->seqBorders[seqUnion->numOfQueries] - leftBorders[seqUnion->numOfQueries] + 1; // length of subject = S0
	for (i = 0; i < seqUnion->numOfQueries; i++) {
		//arguments: args->P, lS, gcQ, gcS for query=Qi and subject=S0 --> important: depends on args->P e [0, 1]
		maxShulens[i] = maxShulenNew(a->P, lS0, seqUnion->gc[i], seqUnion->gc[seqUnion->numOfQueries]);
		for (j = 1; j < seqUnion->numOfSubjects; j++) {
			maxs = maxShulenNew(a->P, seqUnion->seqBorders[j + seqUnion->numOfQueries] - leftBorders[j + seqUnion->numOfQueries] + 1
															, seqUnion->gc[i], seqUnion->gc[seqUnion->numOfQueries + j]);
			if (maxs > maxShulens[i]) {
				maxShulens[i] = maxs; 
			}
		}		
		maxShulens[i] = (Int64)(a->m * maxShulens[i] + 1); 
    if (a->L == -1) {
      a->L = maxShulens[i];
    }

  #if DEBUG
    //fprintf(fpout, "\nQ%d maxshulen = %d", i, maxShulens[i]);
  #endif
    if (a->T > 0) { /* in this case, threshold for a "strong" interval is set to a->T, otherwise it is either maxShulen[i] or 0 (0, if all intervals are included) */
      a->L = maxShulens[i] = a->T;
    }
    else if (a->s == 0) { /* all intervals are included, i.e. no threshold is set */
      a->L = maxShulens[i] = 0;
    }

  #if DEBUG
    fprintf(fpout, " final maxshulen/threshold = %d\n", i, maxShulens[i]);
  #endif
	}		

	// compute lists of query intervals
	traverseLcpTree(lcpTab, sa, seqUnion->seqUnion, seqUnion->numOfSubjects, seqUnion->numOfQueries, seqUnion->seqBorders, leftBorders, strandBorders
		, maxDepth, maxShulens, &root, a);
	
#if DEBUG  
	end3 = clock();
	elapsed_time3 = (double)(end3 - end2) / CLOCKS_PER_SEC;
	if (a->t) {
		printf( "\nLCP-tree traversal calculation: %.2f seconds.\n", elapsed_time3);
	}
#endif
	free(sa);
	free(lcpTab);
	free(maxShulens);

	// print lists of intervals for each query
	ns = seqUnion->numOfSubjects;
	for (i = 0; i < seqUnion->numOfQueries; i++) {			
	    if (fpout) { // suppress printing of interval analysis on stdout as default action
		    fprintf(fpout, "\n-------------- %s --------------\n", seqUnion->seqUnion->headers[i] + 1);	
        fflush(fpout);
      }
      for (j = 0; j < ns; j++) {
        if (root[i * ns + j]) { // not an empty tree
          correctBT(root[i * ns + j], -1, strandBorders[i] - leftBorders[i]);
        }
      }
	}
	
  /* construct a linked list for each interval tree, i.e. a single list for each pair (Qi, Sj), i=1,..,m; j=1,..,n */
  head = getHead(seqUnion->numOfQueries, seqUnion->numOfSubjects);
  constructLists(root, &head, seqUnion->numOfQueries, seqUnion->numOfSubjects);
  //freeBTQueryIntervals(root, seqUnion->numOfQueries, seqUnion->numOfSubjects);
  free(root); //free root of all binary trees (binary trees were already deallocated in constructLists!!)
#if DEBUG  
  printList(head, seqUnion->numOfQueries, seqUnion->numOfSubjects, fpout, strandBorders, leftBorders);
#endif   
	arrayLambda = (double *)/*e*/malloc(sizeof(double) * seqUnion->numOfQueries);
	arrayK = (double *)/*e*/malloc(sizeof(double) * seqUnion->numOfQueries);
	for (i = 0; i < seqUnion->numOfQueries; i++) {
		arrayLambda[i] = DEFAULT_LAMBDA;
		arrayK[i] = DEFAULT_K;
	}
	scoreMatrixAll = allocateScoreMatrixAll(seqUnion->numOfQueries, alphabetSize);
	expectedFreq = allocateScoreMatrixAll(seqUnion->numOfQueries, alphabetSize);

	for (i = 0; i < seqUnion->numOfQueries; i++) {
		for (j = 0; j < alphabetSize; j++) {
			for (k = j + 1; k < alphabetSize; k++) {
				scoreMatrixAll[i][j * alphabetSize + k] = scoreMatrixAll[i][k * alphabetSize + j] = DEFAULT_S; // mismatch
			}
			scoreMatrixAll[i][j * alphabetSize + j] = DEFAULT_M; // match
		}
	}
	ambCharPenaltyAll = DEFAULT_S;

  /* align space between exact matches which are listed in linked lists pointed by head-s*/
	// compute initial alignments with simple match/mismatch scoring to speed up the process????????????????
  aRegionList = getARegionLists(seqUnion->numOfQueries, seqUnion->numOfSubjects);
	fastScore = 1; // fast alignment using fixed match/mismatch reward/penalty
  constructAlignment(seqUnion->seqUnion->seq, head, a, &aRegionList, seqUnion->numOfQueries, seqUnion->numOfSubjects, leftBorders
                    , strandBorders, seqUnion->seqBorders, (long long)seqUnion->seqBorders[seqUnion->numOfQueries + seqUnion->numOfSubjects - 1], 
										arrayLambda, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, fastScore, arrayK);
#if !DEBUG  
  if (fpout) {
    printARegionList(aRegionList, seqUnion->numOfQueries, seqUnion->numOfSubjects, fpout, a->B, seqUnion->seqUnion->headers, strandLen);
  }
#endif
  candidateList = (mNode **)/*e*/malloc(sizeof(mNode *) * seqUnion->numOfQueries);
  // first elements
  getCandidateList(aRegionList, seqUnion->numOfQueries, seqUnion->numOfSubjects, &candidateList, a->B);
  
  finalList = (mNode **)/*e*/malloc(sizeof(mNode *) * seqUnion->numOfQueries);
  constructFinalList(&candidateList, aRegionList, seqUnion->numOfQueries, seqUnion->numOfSubjects, &finalList, a->B, seqUnion->seqUnion->headers, 
          strandBorders, leftBorders, fqout, (long long)seqUnion->seqBorders[seqUnion->numOfQueries + seqUnion->numOfSubjects - 1], strandLen, arrayLambda); // prints final results in fqout
  if (fpout) {
    printCandidateList(finalList, fpout, seqUnion->numOfQueries, " final list of best intervals ", seqUnion->seqUnion->headers, a->B, strandLen);
  }
	maxCols = 0;
	allocateElemsNW(&mat, &firstAligned, &secondAligned, &maxCols, &first, &second, a, 2 * a->f + 2 * a->B);
	/////////////// NEW - 30/01/2015
	getScoringMatrices(finalList, seqUnion->numOfQueries, seqUnion->seqUnion->seq, leftBorders, strandBorders,
		strandLen, scoreMatrixAll, alphabetSize, &nuclFreq, expectedFreq);
	//arrayLambda = (double *)/*e*/malloc(sizeof(double) * seqUnion->numOfQueries);
	for (i = 0; i < seqUnion->numOfQueries; i++) {
    // if nucl. freq-s are 0, then there was no match for a query so set arraylambda[i] to 0 and skip alignment refinement
    if (finalList[i]) {
      arrayLambda[i] = getLambda(&scoreMatrixAll[i][0], alphabetSize, nuclFreq[i]);
			/* calculate K */
			arrayK[i] = getK(&scoreMatrixAll[i][0], expectedFreq[i], arrayLambda[i], alphabetSize);
    }
	}
	// calcualte: (1) refined alignment scores using scoreMatrixAll and (2) E-value based on scoreMatrixAll and lambda
	fastScore = 0; // refined alignment scoring system 
	refineAlignment(finalList, seqUnion->numOfQueries, seqUnion->seqUnion->seq, leftBorders, strandBorders,
		strandLen, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, arrayLambda, maxCols,
		first, second, firstAligned, secondAligned, mat, a, 
		seqUnion->seqBorders[seqUnion->numOfQueries + seqUnion->numOfSubjects - 1], fastScore, arrayK);
	if (fpout) {
		printCandidateList(finalList, fpout, seqUnion->numOfQueries, " final best list - new scoring ", seqUnion->seqUnion->headers, a->B, strandLen);
	}
	/////////////
	/* concatenate elements of the same subject type if less than maxLenApart apart and with positive alignment score
	 * is final merged list with overlapping) */
	fastScore = 0;
	concatenateFinalList(finalList, seqUnion->numOfQueries, a, seqUnion->seqUnion->seq, leftBorders, strandBorders,
                      firstAligned, secondAligned, maxCols, first, second, mat, 
											seqUnion->seqBorders[seqUnion->numOfQueries + seqUnion->numOfSubjects - 1], 
											seqUnion->numOfSubjects, arrayLambda, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, fastScore, arrayK);
  if (fpout) {
    printCandidateList(finalList, fpout, seqUnion->numOfQueries, " final merged list - with overlapping ", seqUnion->seqUnion->headers, a->B, strandLen);
  }

  /************************************************ critical: resolveOverlapsFinalList ***********************/
  allSeqLength = seqUnion->seqBorders[seqUnion->numOfQueries + seqUnion->numOfSubjects - 1] - seqUnion->seqBorders[seqUnion->numOfQueries - 1];
  resolveOverlapsFinalList(fpout, &finalList, seqUnion->numOfQueries, a, seqUnion->seqUnion->seq, leftBorders, strandBorders,
                      firstAligned, secondAligned, maxCols, first, second, mat, allSeqLength/*seqUnion->seqBorders[seqUnion->numOfQueries + seqUnion->numOfSubjects - 1]*/, 
											seqUnion->numOfSubjects, arrayLambda, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, fastScore, arrayK);
  if (fpout) {
    printCandidateList(finalList, fpout, seqUnion->numOfQueries, " final merged list - overlaps excluded ", seqUnion->seqUnion->headers, a->B, strandLen);
  }
  freeElemsNW(mat, firstAligned, secondAligned, first, second);

  /* maxGap should never exceed some percentage of length of the original sequence, since matrix size depends on maxGap^2. Alternative: Hirschberg's algorithm */
  maxGap = 2 * a->f;
  allocateElemsNW(&mat, &firstAligned, &secondAligned, &maxCols, &first, &second, a, maxGap);
  maxCols = maxGap;
  findBP(finalList, seqUnion->numOfQueries, a, seqUnion->seqUnion->seq, leftBorders, strandBorders,
                      firstAligned, secondAligned, maxCols, first, second, mat, 
											seqUnion->seqBorders[seqUnion->numOfQueries + seqUnion->numOfSubjects - 1], seqUnion->numOfSubjects, 
											maxGap, arrayLambda, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, 0, strandLen, arrayK);
  if (fpout) {
    printCandidateList(finalList, fpout, seqUnion->numOfQueries, " final merged list - found break points", seqUnion->seqUnion->headers, a->B, strandLen);
  }
  
  //updateBP(finalList, seqUnion->numOfQueries, a, leftBorders, firstAligned, secondAligned, maxCols, first, second, mat, allSeqLength, seqUnion->numOfSubjects,
		//seqUnion->seqUnion->seq, fastScore, strandBorders, seqUnion->seqBorders, arrayK);
  //if (fpout) {
  //  printCandidateList(finalList, fpout, seqUnion->numOfQueries, " final merged list - added break points", seqUnion->seqUnion->headers, a->B, strandLen);
  //}
  
  ///* print query as a mosaic structure */
  if (fQout) {
    printCandidateList(finalList, fQout, seqUnion->numOfQueries, " mosaic structure ", seqUnion->seqUnion->headers, a->B, strandLen);
    //print seq-s lengths
    printSeqLen(fLen, strandLen, seqUnion->numOfQueries, seqUnion->numOfSubjects);
  }
  if (faout) {
    printAlignment(finalList, faout, seqUnion->numOfQueries, seqUnion->seqUnion->headers, strandLen, seqUnion->seqUnion->seq, leftBorders);
  }

  freeElemsNW(mat, firstAligned, secondAligned, first, second);
  freeMList(finalList, seqUnion->numOfQueries, seqUnion->numOfSubjects); 
  free(candidateList);

  freeARList(aRegionList, seqUnion->numOfQueries); 
	freeNuclFreq(nuclFreq, seqUnion->numOfQueries);
	freeScoringMatrices(scoreMatrixAll, seqUnion->numOfQueries);
	freeScoringMatrices(expectedFreq, seqUnion->numOfQueries);
	free(arrayLambda);
	free(arrayK);
  free(leftBorders);
  free(strandBorders);
  free(strandLen);
  free(minSumWin);
}

/* traverseLcpTree: bottom-up traversal lcp-interval tree */
#if defined(WIN)
static 
#endif
void traverseLcpTree(Int64 *lcpTab, Int64 *sa, Sequence *seq, Int64 numOfSubjects, Int64 numOfQueries, Int64 *seqBorders
											, Int64 *leftBorders, Int64 *strandBorders
											, Int64 maxDepth, Int64 *maxShulens
											, qNode ***root, Args *a) {

  Interval *lastInterval, *interval;
  Int64 i, j, lb, rightEnd, lastIsNull;
  Stack *treeStack = NULL;
	Stack *reserveStack = NULL; // reserve stack for intervals
	Stack *reserveQIStack = NULL; // reserve stack for query intervals
	short *QS; /* array of subject and query indexes corresponding to each position in the SA; pos. values --> subjects, neg.values --> queries */
	int step = STEP;

	//queryLeaves - maximal number is the maximal number of leaves at each interval; queryLeaves[i][0]-query index, queryLeaves[i][1]-position in a query
	Int64 queryLeaves[MAXNUMLEAVES][2]; //queryLeaves[i][0]-query index, queryLeaves[i][1]-position in a query
	//Int64 subjectLeaves[MAXNUMLEAVES][2]; //subjectLeaves[i][0]-subject index, subjectLeaves[i][1]-position in a subject
	int maxCols = 2;
	qNode *qn = NULL;

  // list of subjects; list contains subjects of a single interval; number of subjects in an interval is denoted as numOfSubjInterval <= numOfSubjects
  //Int64 *subjectList = NULL;

	// unresolvedQLeaves - unresolved leaves from the interval's subtree; unresolvedQLeaves[i][0]-query index, unresolvedQLeaves[i][1]-position in a query
	Int64 *unresolvedQLeaves = NULL;
	Int64 maxUnresQLeaves = 100;

  Int64 **subjectLeaves = NULL; /* subject leaves - n arrays; a single 1D array contains subject positions for a single subject in an interval */
	Int64 *maxSLeaves = NULL; /* array of n elements; maximal elements of each of subjectLeaves arrays */
  Int64 *numSLeaves = NULL; /*array of n elements; current number of elements of each of subjectLeaves arrays */
  
  Int64 **subjectLeavesAll = NULL; /* all subject leaves - in a subtree of an interval (not just immediate children) */
	Int64 *maxSLeavesAll = NULL; /* array of n elements; maximal elements of each of subjectLeaves arrays */
  Int64 *numSLeavesAll = NULL; /*array of n elements; current number of elements of each of subjectLeaves arrays */  
  
  // allocate and initialize numSLeaves, maxSLeaves, subjectLeave
  initSubjectLeaves(&subjectLeaves, &maxSLeaves, &numSLeaves, MAXNUMLEAVES, numOfSubjects);
  // allocate and initialize numSLeavesAll, maxSLeavesAll, subjectLeaveAll
  initSubjectLeaves(&subjectLeavesAll, &maxSLeavesAll, &numSLeavesAll, MAXNUMLEAVESALL, numOfSubjects);

  // unresolved query leaves
  unresolvedQLeaves = (Int64 *) /*e*/malloc(maxUnresQLeaves * sizeof(Int64) * 2); // * 2 since each leaf is represented by its query index and a position within a query
	
	*root = getBTQueryIntervals(numOfQueries, numOfSubjects);

	/* initialize auxiliary arrays: subjects and queryLeaves */
  QS = getQS(seqBorders, leftBorders, numOfSubjects, numOfQueries, step);
	
  sa++; /* since data in sa and lcpTab start with position 1, and not 0 */
  lcpTab++;	
  rightEnd = seq->len-1; 
  //lcpTab[0] = 0; /* or -1 */ // moved to getLcpTreeShulens

  treeStack = createStack();    /* true stack */
  lastIsNull = 1; 
  lastInterval = NULL;
  interval = NULL;
  push(treeStack, getInterval(0, 0, rightEnd, NULL, numOfSubjects, NULL, maxDepth)); /* push root node to the stack */

  /* auxiliary stack for (Interval *); used for saving used allocated locations; since stack operations pop/push are faster than malloc */
	reserveStack = createStack(); 
  
	/* auxiliary stack for (queryInterval *) */
	reserveQIStack = createStack(); 

	for (i = 1; i < seq->len; i++){
    lb = i - 1;
    while(lcpTab[i] < ((Interval *)(treeStack->top))->lcp) { /* end of the previous interval, so the interval should be popped from the stack*/						
			/* if the current top is the child of the new node, then pop the top */			
      ((Interval *)(treeStack->top))->rb = i - 1;
      lastInterval = (Interval *)pop(treeStack);			
      if (lastInterval->lcp >= a->L) {
			  process2(lastInterval, seqBorders, 
								leftBorders, strandBorders, numOfSubjects, QS, &queryLeaves[0][0], maxCols, sa, 
                &unresolvedQLeaves, &maxUnresQLeaves, reserveQIStack, maxShulens, 
								subjectLeaves, maxSLeaves, numSLeaves, numOfQueries, root, 
                subjectLeavesAll, maxSLeavesAll, numSLeavesAll);
                //*root); 			
      }
			lb = lastInterval->lb;

			/* save child intervals of popped interval */
      for(j = 0; j < lastInterval->numChildren; j ++) {
				lastInterval->children[j]->numChildren = 0;
				lastInterval->children[j]->parent = NULL;
				push(reserveStack, (void *)lastInterval->children[j]);
      }
      lastIsNull = 0;
			if (lcpTab[i] <= ((Interval *)treeStack->top)->lcp) { /* the new top is the parent of the ex top*/
				lastInterval->parent = (Interval *)treeStack->top;
				addChild((Interval *)treeStack->top, lastInterval);
				lastIsNull = 1;
			}			
    } /* end while */

 		if (lcpTab[i] > ((Interval *)treeStack->top)->lcp) { /* add interval to the stack */
			if (isEmpty(reserveStack)) { 
				//interval = getInterval(lcpTab[i], lb, rightEnd, NULL, numOfSubjects, treeStack->top, maxDepth);// treeStack->top or null
				interval = getInterval(lcpTab[i], lb, rightEnd, NULL, numOfSubjects, NULL, maxDepth);
			}
			else { /* use locations from the reserveStack */ 
				interval = (Interval *)pop(reserveStack); // pop the last child of the lastInterval from the reservestack
				interval->lcp = lcpTab[i];
				interval->lb = lb;
				interval->rb = rightEnd;
				interval->numChildren = 0;
				interval->parent = NULL;
      }
      if (!lastIsNull){ 
				lastInterval->parent = interval;
				addChild(interval, lastInterval);
				lastIsNull = 1;
      }
      push(treeStack, interval);
    }
  }

#if DEBUG  
	//printf("Empting stack...\n");
	//printf("Number of intervals allocated:%lld\n\n", (long long)numInterval);
#endif

	while(!isEmpty(treeStack)) { 
    interval = (Interval *) pop(treeStack);    
		/* when the stack is not empty and the current interval has no parent, then his parent is on the top of the stack*/
		if (!isEmpty(treeStack) && !interval->parent && interval->lcp > ((Interval *)(treeStack->top))->lcp) { 
			interval->parent = (Interval *) treeStack->top;
			addChild((Interval *)treeStack->top, interval);
		}
		/* process: 
		 * 1) it labels an interval according to whether it has ST leaves from query (isQuery) or from subject (isSubject) or both (isQuery && isSubject).
		 * 2) if(isQuery && isSubject) it determines the corresponding query shustring lengths (if any) */
    //if (lastInterval->lcp > a->L) { // bug!!
    if (interval->lcp > a->L) {
		  process2(interval, seqBorders, leftBorders, strandBorders, numOfSubjects, QS, &queryLeaves[0][0], maxCols, sa, &unresolvedQLeaves, 
						  &maxUnresQLeaves, reserveQIStack, maxShulens, subjectLeaves, maxSLeaves, numSLeaves, numOfQueries, root,
              subjectLeavesAll, maxSLeavesAll, numSLeavesAll);

    }
		for (i = 0; i < interval->numChildren; i++) {
      freeInterval((void *)interval->children[i]);
		}
		if (interval->parent == NULL) {
			freeInterval((void *)interval); // root
		}
    else {
		  interval->numChildren = 0;
    }
	}  
	//freeInterval((void *)interval); /* free the root interval; it has to be done here, what is different from kr 2!! */
  freeStack(treeStack, freeInterval);
  
	/* free memory allocated from the reserveStack */
	while (!isEmpty(reserveStack)) {
    interval = (Interval *) pop(reserveStack);
		for (i = 0; i < interval->numChildren; i++) {
    	freeInterval((void *)interval->children[i]);
		}
    //freeInterval((void *)interval); /* free the root interval */
		if (!interval->parent) {
			freeInterval((void *)interval); 
		}
  }
	/* free reserve stack */
	freeStack(reserveStack, freeInterval);
  //printf("\n6");

	/* free reserve queryInterval stack */
	while (!isEmpty(reserveQIStack)) {
		qn = (qNode *) pop(reserveQIStack);
		qn->left = qn->right = NULL;
		freeQNode(qn); 
	}
	freeStack(reserveQIStack, freeQNode);
	

#if DEBUG  
	//printf("Finished.. Number of intervals allocated:%lld\n\n", (long long)numInterval);
#endif
	
	free(QS);
	free(unresolvedQLeaves);
	freeSubjectLeaves(subjectLeaves, numOfSubjects);
	freeSubjectLeaves(subjectLeavesAll, numOfSubjects);
	free(maxSLeaves);
  free(numSLeaves);
	free(maxSLeavesAll);
  free(numSLeavesAll);
}


/* adding a new query (or subject) leaf to the interval*/
#if defined(WIN)
static 
#endif
void addQueryLeaf(Int64 *queryLeaves, int maxCols, Int64 *numQueryLeaves, Int64 queryIndex, Int64 pos) {

	queryLeaves[*numQueryLeaves * maxCols /*+ 0*/] = queryIndex;
	queryLeaves[*numQueryLeaves * maxCols + 1] = pos; // relative within that query sequence, that is, positions always start with 0
	++ (*numQueryLeaves);
}



/* add query leaves or denote subject leaves */
#if defined(WIN)
static 
#endif
void checkLeaves2 (Interval *interval, Int64 *seqBorders, Int64 lb, Int64 rb, short *QS, Int64 *queryLeaves, int maxCols,
										//Int64 *numQueryLeaves, Int64 *leftBorders, Int64 *sa, Int64 *subjectList, Int64 *numOfSubjInterval) {
										Int64 *numQueryLeaves, Int64 *leftBorders, Int64 *sa, 
                    Int64 **subjectLeaves, Int64 *numSLeaves, Int64 *maxSLeaves, Int64 numOfQueries) {
 
	Int64 i, k;
	for (i = lb; i <= rb; i ++) { /* leaf of an lcp tree in the terms of ST is the inner node which contains only leaves */
		k = QS[sa[i]]; //k = findSubject(seqBorders, suffixArray[i], subjects);
		if (queryLeaves != NULL && k < 0) { // query --> (suffixArray[i] >= queryStart && suffixArray[i] <= queryEnd) {
			k = -k; // k starts with 1, 2, ...
			addQueryLeaf(queryLeaves, maxCols, numQueryLeaves, k - 1, sa[i] - leftBorders[k - 1]); // position within a query - relative to the left border
		}
		else if (k > 0) { // subject
			-- k; // k starts with 1, 2..
			/* subject is now a bit in a bit-vector, and not a value on its own */
			interval->subjectIndex[k / WORDSIZE] |= (MASK_ONE << k % WORDSIZE); 
      if (numSLeaves[k] >= maxSLeaves[k]) {
        maxSLeaves[k] *= 2;
        subjectLeaves[k] = (Int64 *)/*e*/realloc(subjectLeaves[k], maxSLeaves[k] * sizeof(Int64));
      }
      subjectLeaves[k][numSLeaves[k]] = sa[i] - leftBorders[k + numOfQueries];
      ++ numSLeaves[k];
			//addQueryLeaf(subjectLeaves, maxCols, numSubjectLeaves, k, sa[i] - leftBorders[k + numOfQueries]); // position within a subject - relative to the left border
		}
	}
}

/* determine whether interval has Q and/or S leaves; for subject leaves just set a flag for Sj subject, and
 * for query leaves add thm to the list of query leaves (queryLeaves) */
void determineQS(Interval *interval, Int64 numOfSubjects, Int64 *seqBorders, Int64 *leftBorders, short *subjects, Int64 *queryLeaves, 
								 //int maxCols, Int64 *numQueryLeaves, Int64 *sa, Int64 *subjectList, Int64 *numOfSubjInterval) {
								 int maxCols, Int64 *numQueryLeaves, Int64 *sa,
                 Int64 **subjectLeaves, Int64 *numSLeaves, Int64 *maxSLeaves, Int64 numOfQueries) {

  Int64 i;
	
    // /* a leaf of lcp-interval tree */
	if (interval->numChildren == 0) {
    checkLeaves2(interval, seqBorders, interval->lb, interval->rb, subjects, queryLeaves, maxCols, numQueryLeaves, leftBorders, sa, subjectLeaves, numSLeaves, maxSLeaves, numOfQueries);
	}
  /* an internal node of lcp-interval tree */
	else {		
		// find query leaves left to the left-most interval
		if (interval->lb < interval->children[0]->lb) { // is there a leaf before the left-most child
			checkLeaves2(interval, seqBorders, interval->lb, interval->children[0]->lb - 1, subjects, queryLeaves, maxCols, numQueryLeaves, leftBorders, sa,
        subjectLeaves, numSLeaves, maxSLeaves, numOfQueries);

		}    
		// is there a leaf between intervals
		for (i = 0; i < interval->numChildren - 1; i++) {
			if (interval->children[i]->rb + 1 < interval->children[i + 1]->lb) { 
				checkLeaves2(interval, seqBorders, interval->children[i]->rb + 1, interval->children[i + 1]->lb - 1, subjects, queryLeaves, maxCols, numQueryLeaves, leftBorders, sa, 
			//checkLeaves2(interval, seqBorders, interval->lb, interval->children[0]->lb - 1, subjects, queryLeaves, maxCols, numQueryLeaves, leftBorders, sa, subjectList, numOfSubjInterval);
          subjectLeaves, numSLeaves, maxSLeaves, numOfQueries);
			}
		}

    // right-most child node
		if (interval->children[interval->numChildren - 1]->rb < interval->rb) {			
				checkLeaves2(interval, seqBorders, interval->children[interval->numChildren - 1]->rb + 1, interval->rb, subjects, queryLeaves, maxCols, numQueryLeaves, leftBorders, sa,
          subjectLeaves, numSLeaves, maxSLeaves, numOfQueries);
		}
	}
}

// reallocate unresolved leaves
void reallocUnresLeaves(Int64 **unresolvedQLeaves, Int64 *maxUnresQLeaves, Int64 numUnresolvedQLeaves, int multiplyBy) {

  //memory is reallocated only when the previous amount is exceeded (to avoid reallocation in each iteration)
	while (numUnresolvedQLeaves > *maxUnresQLeaves) {
		*maxUnresQLeaves *= 2;
  }
	*unresolvedQLeaves = (Int64 *)/*e*/realloc(*unresolvedQLeaves, (*maxUnresQLeaves) * multiplyBy * sizeof(Int64)); // * 2 since both queryIndex and the position are stored for each leaf
	if (*unresolvedQLeaves == NULL) {
		eprintf("Memory allocation failed when reallocating unresolvedQLeaves!\n");
	}
}


// process each interval
void process2(Interval *interval, Int64 *seqBorders, 
							Int64 *leftBorders, Int64 *strandBorders, Int64 numOfSubjects, short *QS, 
							Int64 *queryLeaves, int maxCols, Int64 *sa, Int64 **unresolvedQLeaves, Int64 *maxUnresQLeaves, Stack *reserveQIStack, Int64 *maxShulens, 
              Int64 **subjectLeaves, Int64 *maxSLeaves, Int64 *numSLeaves, Int64 numOfQueries, qNode ***root,
              Int64 **subjectLeavesAll, Int64 *maxSLeavesAll, Int64 *numSLeavesAll) {

	Int64 i, j/*, query_lb, query_rb, query_sl*/;
	Int64 numOfSubjWords = numOfSubjects / WORDSIZE + 1;
  Int64 numQueryLeaves = 0;
	Int64 numUnresolvedQLeaves = 0; // unresolved query leaves from child subtrees
 // Int64 numSubjectLeaves = 0;
	Int64 numUnresolvedSLeaves = 0; // unresolved subject leaves from child subtrees
  //Int64 numOfSubjInterval = 0;

  // subject list is initally empty

  /* initialize - begin */
	interval->isQuery = 0;
  interval->isSubject = 0;
	for (j = 0; j < numOfSubjWords; j++) {
		interval->subjectIndex[j] = (Word)0;
	}
  for (i = 0; i < numOfSubjects; i++) {
    numSLeaves[i] = 0;
    numSLeavesAll[i] = 0;
  }
  /* initialize - end */

	/* determine query and subject leaves, and all the subjects that belong to this interval */
	determineQS (interval, numOfSubjects, seqBorders, leftBorders, QS, queryLeaves, maxCols, &numQueryLeaves, sa, subjectLeaves, numSLeaves, maxSLeaves, numOfQueries);
  
  /* NEW: determine all subject leaves in the subtree of the current interval including its own immediate leaves */
  checkLeaves2(interval, seqBorders, interval->lb, interval->rb, QS, NULL, maxCols, &numQueryLeaves, leftBorders, sa, subjectLeavesAll, numSLeavesAll, maxSLeavesAll, numOfQueries);

	// does interval have query leaves anywhere in the subtree 
	if (numQueryLeaves > 0) {
		interval->isQuery = 1;
	}
	else {
		for (i = 0; i < interval->numChildren; i++) {
			if (interval->children[i]->isQuery == 1) {
				interval->isQuery = 1;
				break;
			}
		}
	}
		
	// does interval have subject leaves in the subtree (in determineQS, subjects coming from immediate leaves were already set)
	for (j = 0; j < numOfSubjWords; j++) {
		for (i = 0; i < interval->numChildren; i++) {
			interval->subjectIndex[j] |= interval->children[i]->subjectIndex[j];
		}
	}
	for (j = 0; j < numOfSubjWords; j++) {
		if (interval->subjectIndex[j]) {
			interval->isSubject = 1;
			break;
		}
	}

	// if an interval doesn't have both query and subject leaves
	if (interval->isQuery == 0 || interval->isSubject == 0) {
    return;
	}
  
	// If exists a child interval, that has only query or only subject intervals in its subtree, then these leaves need to be resolved
	numUnresolvedSLeaves = 0; // redundant??????'
  numUnresolvedQLeaves = 0;
	for (i = 0; i < interval->numChildren; i++) {
		if (interval->children[i]->isSubject == 0) { //only query leaves
			numUnresolvedQLeaves += interval->children[i]->rb - interval->children[i]->lb + 1; 
		}
		if (interval->children[i]->isQuery == 0) { // only subject leaves --????? what if only Q1 leaves exist and not Q2, then subject leaves cannot be resolved vs Q2 BUG!!!
			numUnresolvedSLeaves += interval->children[i]->rb - interval->children[i]->lb + 1; 
		}
	}
	
  /* unresolved subject positions - start */
  for (i = 0; i < numOfSubjects; i++) {
	  if (numUnresolvedSLeaves + numSLeaves[i] > maxSLeaves[i]) { // although numUnresolvedSLeaves for a particular Si can be < total numUnresolvedSLeaves
	  //if (numSLeaves[i] > maxSLeaves[i]) { // although numUnresolvedSLeaves for a particular Si can be < total numUnresolvedSLeaves
      // (this is done here to avoid later reallocation for each subject in particular)
		  
      //memory is reallocated only when the previous amount is exceeded (to avoid reallocation in each iteration)
      reallocUnresLeaves(&subjectLeaves[i], &maxSLeaves[i], numUnresolvedSLeaves + numSLeaves[i], 1); // ?????
    }
  }

  numUnresolvedSLeaves = 0; // redundant?
  for (i = 0; i < interval->numChildren; i++) {
	  if (interval->children[i]->isQuery == 0) {
      checkLeaves2(interval, seqBorders, interval->children[i]->lb, interval->children[i]->rb, QS, NULL, maxCols, NULL,
								  leftBorders, sa, subjectLeaves, numSLeaves, maxSLeaves, numOfQueries);

      //for (j = interval->children[i]->lb; j <= interval->children[i]->rb; j ++) {

      //}
    }
  }		
  /* unresolved subject positions - end */


  /* when unresolved query leaves exist */
	if (numUnresolvedQLeaves > 0) {
		//memory is reallocated only when the previous amount is exceeded (to avoid reallocation in each iteration)
    reallocUnresLeaves(unresolvedQLeaves, maxUnresQLeaves, numUnresolvedQLeaves, 2);

		numUnresolvedQLeaves = 0;
		for (i = 0; i < interval->numChildren; i++) {
			if (interval->children[i]->isSubject == 0) {
				checkLeaves2 (interval, seqBorders, interval->children[i]->lb, interval->children[i]->rb, QS, *unresolvedQLeaves, maxCols,
										&numUnresolvedQLeaves, leftBorders, sa, NULL, NULL, NULL, numOfQueries);
			}
		}		
		//resolveQueryIntervals(numUnresolvedQLeaves, *unresolvedQLeaves, maxCols, seqBorders, leftBorders, strandBorders, numOfSubjects, 
		//											interval, reserveQIStack, maxShulens, root, subjectLeaves, numSLeaves);
		resolveQueryIntervals(numUnresolvedQLeaves, *unresolvedQLeaves, maxCols, seqBorders, leftBorders, strandBorders, numOfSubjects, 
													interval, reserveQIStack, maxShulens, root, subjectLeavesAll, numSLeavesAll);
	}
	
	// immediate query leaves
	if (numQueryLeaves > 0) {
		resolveQueryIntervals(numQueryLeaves, queryLeaves, maxCols, seqBorders, leftBorders, strandBorders,
														//numOfSubjects, interval, reserveQIStack, maxShulens, root, subjectLeaves, numSLeaves);	
														numOfSubjects, interval, reserveQIStack, maxShulens, root, subjectLeavesAll, numSLeavesAll);	
	}
	
#if DEBUG
	//if(interval->isQuery && interval->isSubject) {
 //   printf("%lld-[%lld..%lld]sq\n",(long long)interval->lcp,(long long)interval->lb,(long long)interval->rb);
	//	if (interval->parent) printf("\tparent:%lld-[%lld..%lld]sq\n",(long long)interval->parent->lcp,(long long)interval->parent->lb,(long long)interval->parent->rb);
	//	printf("subjects: "); 
	//	for (i = 0; i <= numOfSubjects / WORDSIZE; i++) printf("%d ",interval->subjectIndex[i]);
	//	printf("\n"); 
	//}
	//else if(interval->isQuery) {
 //   printf("%lld-[%lld..%lld]q\n",(long long)interval->lcp,(long long)interval->lb,(long long)interval->rb);
	//	if (interval->parent) printf("\tparent:%lld-[%lld..%lld]sq\n",(long long)interval->parent->lcp,(long long)interval->parent->lb,(long long)interval->parent->rb);
	//}
	//else if(interval->isSubject) {
 //   printf("%lld-[%lld..%lld]s\n",(long long)interval->lcp,(long long)interval->lb,(long long)interval->rb);
	//	if (interval->parent) printf("\tparent:%lld-[%lld..%lld]sq\n",(long long)interval->parent->lcp,(long long)interval->parent->lb,(long long)interval->parent->rb);
	//	printf("subjects: "); 
	//	for (i = 0; i <= numOfSubjects / WORDSIZE; i++) printf("%d ",interval->subjectIndex[i]);
	//	printf("\n"); 
	//}
#endif
}

/* add new query intervals to list of query intervals --> calling addNode */
void resolveQueryIntervals(Int64 numQueryLeaves, Int64 *queryLeaves, int maxCols, Int64 *seqBorders, Int64 *leftBorders, Int64 *strandBorders,
													Int64 numOfSubjects, Interval *interval, Stack *reserveQIStack, Int64 *maxShulens, qNode ***root,
                          Int64 **subjectLeaves, Int64 *numSLeaves) {
	
	Int64 j, i; // i = queryIndex
	Int64 query_lb, query_rb, query_sl;
	//Int64 strandBorder;
	qNode *newNode = NULL;
  Int64 k, indIT;

	query_sl = interval->lcp + 1; // sl =shulen

	for (j = 0; j < numQueryLeaves; j++) {
		i = queryLeaves[j * maxCols];
		if (query_sl < maxShulens[i]) {
		  continue; // ignore intervals that have max shulens shorter than expected by chance alone or whatever the threshold is set to
		}  

		// shulen = interval->lcp + 1
		query_lb = queryLeaves[j * maxCols + 1]; // lb = position in a query, rb = where the value of lcp is 1 
		query_rb = query_lb + interval->lcp; // + 1 - 1;

		/* check for the border of the strand or the whole query !!!! */
		//if (query_rb + leftBorders[i] >= seqBorders[i]) {
		//	query_sl = seqBorders[i] - (leftBorders[i] + query_lb) + 1;
		//	query_rb = seqBorders[i] - leftBorders[i]; // relative end of the sequence
		//	if (query_sl == 0) {
		//		printf("query_sl == 0");
		//	}
		//}
		if (query_lb + leftBorders[i] <= strandBorders[i] && query_rb + leftBorders[i] > strandBorders[i]) { // it also covers query_rb + leftBorders[i] > seqBorders[i]
			query_sl = strandBorders[i] - (leftBorders[i] + query_lb) + 1;
			query_rb = strandBorders[i] - leftBorders[i]; // relative end of the strand		
			if (query_sl == 0) {
				printf("query_sl == 0");
			}
		}
		else if (query_lb + leftBorders[i] > strandBorders[i] && query_lb + leftBorders[i] <= seqBorders[i] && query_rb + leftBorders[i] > seqBorders[i]) {
			query_sl = seqBorders[i] - (leftBorders[i] + query_lb) + 1;
			query_rb = seqBorders[i] - leftBorders[i]; // relative end of the strand		
			if (query_sl == 0) {
				printf("query_sl == 0");
			}
		}
		//query_slAvg = (double)(2 * query_sl - (query_rb - query_lb)) / 2; // (query_sl + query_sl - (query_rb - query_lb)) / 2
    if (query_lb > query_rb) {
	    eprintf("[ERROR]: Query interval borders are incorrect!");
    }

    for (k = 0; k < numOfSubjects; k++) {
      if (numSLeaves[k]) { // if there are some Sk positions attached to this interval
        /* get interval from a reserve stack and not allocate memory each time */
        if (isEmpty(reserveQIStack)) { 
          //newNode	= getQNode(query_sl, query_lb, query_rb, subjectLeaves, numSLeaves, maxCols, k);
          newNode	= getQNode(query_sl, query_lb, query_rb, subjectLeaves[k], numSLeaves[k]);
        }
        else { /* use locations from the reserveStack */ 
          newNode	= (qNode *)pop(reserveQIStack); // pop the last child of the lastInterval from the reservestack
          setQNode(newNode, query_sl, query_lb, query_rb, NULL, NULL, subjectLeaves[k], numSLeaves[k]);
        }
        /*****************************/    
        indIT = i * numOfSubjects + k; // index of an interval tree for Qi and Sk
        (*root)[indIT] = addNode((*root)[indIT], newNode, reserveQIStack/*, numOfSubjects*/);
      }
    }
	}
}

/* auxiliary function to print out sa and lcp */
#if defined(WIN)
static 
#endif
void printSA_LCP(FILE *fpout, Int64 *sa, Int64 *lcp, Int64 seqUnionLen) {
	Int64 i;

	fprintf(fpout, "\n---------- SA -----------\n");
	for (i = 1; i <= seqUnionLen; i++) {
		fprintf(fpout, "sa[%lld]=%lld\n", (long long)(i - 1), (long long)sa[i]);
	}
	fprintf(fpout, "\n---------- LCP -----------\n");
	for (i = 1; i <= seqUnionLen; i ++) {
		fprintf(fpout, "lcp[%lld]=%lld\n", (long long)(i - 1), (long long)lcp[i]);
	}	
  fflush(fpout);
}

/* find a subject/query for each position in the suffix array */
#if defined(WIN)
static 
#endif
short *getQS(Int64 *seqBorders, Int64 *leftBorders, Int64 numOfSubjects, Int64 numOfQueries, int step) {
	
	Int64 i, j;
	short *QS = NULL;

	QS = (short *)/*e*/malloc(sizeof(short) * (seqBorders[numOfSubjects + numOfQueries - 1] + 1));
	
	// query - query strands have negative values to distinct them from subject indexes
	for (i = 0; i < numOfQueries; i ++) {
		for (j = leftBorders[i]; j <= seqBorders[i]; j ++) {
			QS[j] = -(short)(i + 1);
		}
	}

	// subject - positive values
	for (i = 0; i < numOfSubjects; i ++) {
		//for (; j <= seqBorders[i] / step; j ++) {
		for (j = leftBorders[numOfQueries + i]; j <= seqBorders[numOfQueries + i]; j ++) {
			QS[j] = (short)i + 1;
		}
	}
//#if DEBUG
#if 0
	printf("QS:\n");
	for (j = 0; j <= seqBorders[numOfSubjects + numOfQueries - 1]; j ++) {
		printf("[%lld] = %hd\n", j, QS[j]);
	}
#endif
	return QS;
}

// initialize subjectLeaves and accompanying arrays
void initSubjectLeaves(Int64 ***subjectLeaves, Int64 **maxSLeaves, Int64 **numSLeaves, Int64 maxNumLeaves, Int64 numOfSubjects) {

  Int64 i;
  
  *maxSLeaves = (Int64 *)/*e*/malloc(numOfSubjects * sizeof(Int64));
  *numSLeaves = (Int64 *)/*e*/malloc(numOfSubjects * sizeof(Int64));
  *subjectLeaves = (Int64 **)/*e*/malloc(numOfSubjects * sizeof(Int64 *));

  for (i = 0; i < numOfSubjects; i++) {
    (*maxSLeaves)[i] = maxNumLeaves * 2; // initially, but it might be changed in the program
    (*numSLeaves)[i] = 0;
    (*subjectLeaves)[i] = (Int64 *)/*e*/malloc((*maxSLeaves)[i] * sizeof(Int64));
  }
}

// deallocate subjectLeaves
void freeSubjectLeaves(Int64 **subjectLeaves, Int64 numOfSubjects) {

  Int64 i;

  for (i = 0; i < numOfSubjects; i++) {
    free(subjectLeaves[i]);
  }
  free(subjectLeaves);
}

// print sequences lengths
void printSeqLen(FILE *fLen, Int64 *strandLen, Int64 numOfQueries, Int64 numOfSubjects) {
  Int64 i;
  
  // print number of queries and queries' lengths
  fprintf(fLen, "#%lld\n", (long long) numOfQueries);
  for (i = 0; i < numOfQueries; i ++) {
    fprintf(fLen, "%lld\n", (long long) strandLen[i]);    
  }
  // print number of subjects and subjects' lengths
  fprintf(fLen, "#%lld\n", (long long) numOfSubjects);
  for (i = 0; i < numOfSubjects; i ++) {
    fprintf(fLen, "%lld\n", (long long) strandLen[i + numOfQueries]);    
  }
}


