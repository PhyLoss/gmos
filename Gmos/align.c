/***** align.c *************************************************************
 * Description: Functions for alignment.
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
 *
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
#include "expectedShulen.h"
#include "intervalStack.h"
#include "interval.h"

#include "interface.h"
#include "nw.h"
#include "queryBTNode.h"
#include "lNode.h"
#include "align.h"

#if defined(_DEBUG) && defined(WIN) 
#include "leakWatcher.h"
#endif

#if defined(_DEBUG) && defined(WIN) 
  #define new DEBUG_NEW
  #undef THIS_FILE
  static char THIS_FILE[] = __FILE__;
#endif


/* aNode describes an exact match between [Q.lb, Q.lb + Q.sl - 2] and [S.lb, S.lb + Q.sl - 2] */
aNode *getANode(lNode *p) {
  aNode *a;

  // new node in list
  a = (aNode *)/*e*/malloc(sizeof(aNode));
  /* find starting positions (lb) of matching segments in Q and S */
  a->lbQ = p->lb;
  a->lbS = p->listPositions[0]; // just take first element in subject list; this could be changed later
  a->sl = p->sl;
  a->next = NULL;
  return a;
}

/* aRegion contains list of 1 or more anode-s. If more than one anode, then the anodes are connected with aligned positions (see alignedSegment). 
 * Positions between anode-s are aligned Needleman-Wunsch algorithm and have to be shorter than arg->A */
aRegion *getARegion(aNode *a, float match) {
  
  aRegion *r;
  int i;

  // new node in list
  r = (aRegion *)/*e*/malloc(sizeof(aRegion));
  r->start = a; /* points to the first anode in region */
  r->end = a;   /* points to the last anode in region */
  r->score = (a->sl - 1) * match; /* initial score */
  
  r->numAlignedSegment = 0;
  r->maxAlignedSegment = 1;
  r->scoreAlignedSegment = (double *) /*e*/malloc(sizeof(double) * r->maxAlignedSegment);
  r->alignedSegment = (char **) /*e*/malloc(sizeof(char *) * 2 * r->maxAlignedSegment); /* *2 since there are 2 parts that are aligned: a query part and a subject part */
  for (i = 0; i < 2 * r->maxAlignedSegment; i++) {
    r->alignedSegment[i] = NULL;
  }

  ////////// NEW ////////////
  r->alignedStart = (char **)/*e*/malloc(sizeof(char *) * 2); /* aligned segment just before aRegion --> Q: alignedStart[0], S: alignedStart[1] */
  r->alignedEnd = (char **)/*e*/malloc(sizeof(char *) * 2); /* aligned segment just after aRegion --> Q: alignedEnd[0], S: alignedEnd[1] */
  for (i = 0; i < 2; i++) {
    r->alignedStart[i] = NULL;
    r->alignedEnd[i] = NULL;
  }
  r->scoreStart = r->scoreEnd = 0;
  ////////// ?? why is this not initialized here?
  //r->startPosQ = r->endPosQ = -1;
  //r->startPosS = r->endPosS = -1;
  r->startPosQ = r->start->lbQ;
  r->endPosQ = r->end->lbQ + r->end->sl - 2;
  r->startPosS = r->start->lbS;
  r->endPosS = r->end->lbS + r->end->sl - 2;

  r->totalLen = a->sl - 1;
  r->ratio = r->score / r->totalLen;
  r->next = NULL;
  r->active = 1;
  return r;
}

/* pointers to first elements of |Q| * |S| aregions' lists */
aRegion **getARegionLists(Int64 numOfQueries, Int64 numOfSubjects) {
	
  aRegion **ar; // head of |Q| * |S| lists
	Int64 i, j;

	ar = (aRegion **)/*e*/malloc(numOfQueries * numOfSubjects * sizeof(aRegion *));
	/* each element in the list is a pointer to the first element of a linked list */
	for (i = 0; i < numOfQueries; i++) {
    for (j = 0; j < numOfSubjects; j++) {
      ar[i * numOfSubjects + j] = NULL; // initially: all lists are empty	
    }
  }
	return ar;
}


void freeARegion(aRegion *r) {
  
  int i;
  aNode *a, *b;
  
  a = r->start;
  while (a) {
    b = a->next;
    free(a);
    a = b;
  }

  for (i = 0; i < 2 * r->maxAlignedSegment; i++) {
    free(r->alignedSegment[i]);
  }
  free(r->alignedSegment);
  free(r->scoreAlignedSegment);
  
  //// NEW
  for (i = 0; i < 2; i++) {
    free(r->alignedStart[i]);
    free(r->alignedEnd[i]);
  }
  free(r->alignedStart);
  free(r->alignedEnd);  
  free(r);
}

/* free aNode-s up from m to n; if something remains after n --> construct new aRegion r2 */
void _freeARegionFromMToN(aRegion *r, aNode *m, aNode *n, long long allSeqLength, aRegion *r2/*, double lambda*/) {
  
  int i, j;
  aNode *a = NULL, *b = NULL, *prev;
  int cntAlignedSegments = 0;
  int cntAlignedSegmentsLeft = 0;
  int cntAlignedSegmentsRight = 0;
  int oldNumAlignedSegment = r->numAlignedSegment;

  a = r->start;
  prev = NULL;
  r->score = 0;
  while (a != m) { // find m
    prev = a;
    r->score += a->sl - 1;
    if (cntAlignedSegmentsLeft < r->numAlignedSegment) {
      r->score += r->scoreAlignedSegment[cntAlignedSegmentsLeft];
      ++ cntAlignedSegmentsLeft;
    }
    a = a->next;
  }
  if (cntAlignedSegmentsLeft >= 1) {
    r->score -= r->scoreAlignedSegment[cntAlignedSegmentsLeft - 1]; // last alignment is not included
  }

  while (a) { // delete from m up to and including n
    if (a == n) {
      if (!b || m != n) { // b could be possibly equal to a (when n == m); in that case, skip the following while-loop
        b = a->next;
        free(a);
        a = NULL;
      }
      break;
    }
    // update totalLen and score for match and following aligned segment
    b = a->next;
    free(a);
    a = b;
    b = b->next;
    ++ cntAlignedSegments;
  } // end while
  
  // transfer aligned segments and matches starting with b to m2    
  if (r2 != NULL) {
    //if (b) { // b could be possibly equal to a (when n == m); in that case, skip the following while-loop
    b = b->next; // since b was already added to r2 (before calling this function)
    while (b != NULL) {
      r2->end->next = b; 
      r2->end = b;
      r2->score += b->sl - 1;
      ++ cntAlignedSegmentsRight;
      b = (b == r->end) ? NULL : b->next;
    } // end while
    //}

    r2->maxAlignedSegment = r->maxAlignedSegment;
    checkMaxAlignedSegments(r2);
    for (i = 0; i < r2->maxAlignedSegment * 2; i++) {
      r2->alignedSegment[i] = NULL;
    }
    // copy rest of aligned segments from r
    r2->numAlignedSegment = cntAlignedSegmentsRight;
    for (i = oldNumAlignedSegment - cntAlignedSegmentsRight, j = 0; i < oldNumAlignedSegment; i++, j++) {
      r2->alignedSegment[2 * j] = (char *)/*e*/realloc(r2->alignedSegment[2 * j], strlen(r->alignedSegment[2 * i]) + 1);
      r2->alignedSegment[2 * j + 1] = (char *)/*e*/realloc(r2->alignedSegment[2 * j + 1], strlen(r->alignedSegment[2 * i + 1]) + 1);
      strcpy(r2->alignedSegment[2 * j], r->alignedSegment[2 * i]);
      strcpy(r2->alignedSegment[2 * j + 1], r->alignedSegment[2 * i + 1]);
      r2->scoreAlignedSegment[j] = r->scoreAlignedSegment[i];
      r2->score += r2->scoreAlignedSegment[j];
    }

    // update r2 values -- start -- 
    r2->startPosQ = r2->start->lbQ;
    r2->startPosS = r2->start->lbS;
    r2->endPosQ = r2->end->lbQ + r2->end->sl - 2;
    r2->endPosS = r2->end->lbS + r2->end->sl - 2;
    r2->totalLen = r2->endPosQ - r2->startPosQ + 1;
    r2->ratio = r2->score / r2->totalLen;
    computeEValueScoreBit(r2, allSeqLength/*, lambda*/);
    // update r2 values -- end -- 
  }

  // update r (if r remained) values -- start -- 
  if (prev) { // connect prev (last from the left to be kept)
    r->end = prev;
    r->end->next = prev->next = NULL;
    r->startPosQ = r->start->lbQ;
    r->startPosS = r->start->lbS;
    r->endPosQ = r->end->lbQ + r->end->sl - 2; 
    r->endPosS = r->end->lbS + r->end->sl - 2; 
    r->totalLen = r->endPosQ - r->startPosQ + 1;
    r->scoreStart = r->scoreEnd = 0;
    r->numAlignedSegment = cntAlignedSegmentsLeft - 1;
    for (i = 0; i < 2; i++) {
      free(r->alignedStart[i]);
      free(r->alignedEnd[i]);
      r->alignedStart[i] = r->alignedEnd[i] = NULL;
    }
    r->ratio = r->score / r->totalLen;
    computeEValueScoreBit(r, allSeqLength/*, lambda*/);
  }
  else { // if m was starting node, then m will be deleted completely (in shortenM)
    r->end = r->start = NULL;
  }
  // update r values -- end --

}



/* free list of aRegion-s -- numOfQueries lists for Sj */
void freeARList(aRegion **arList, Int64 numOfQueries/*, Int64 numOfSubjects*/) {
  
  Int64 i;
  aRegion *m, *n;

	for (i = 0; i < numOfQueries; i++) {
    m = arList[i];
    while (m) {
      n = m->next;
      freeARegion(m);
      m = n;
    } // end while
  }
  free(arList);
}

/* check whether aRegion is shorter than set by -f option, or has score less than set by -F option; in that case return 0, else return 1 */
int checkARegion (aRegion *r, Args *arguments) {

  int retValue = 1;
  Int64 rStart, rEnd, rLen;
  rStart = r->startPosQ;
  rEnd = r->endPosQ;
  rLen = rEnd - rStart + 1;

  if (rLen < arguments->f || r->score < arguments->F || r->ratio < arguments->R) {
    retValue = 0;
  }
  return retValue;
}

/* check whether number of aligned segments in aRegion are > maximum; in that case, increase maximum and the size of scoreAlignedSegment */
void checkMaxAlignedSegments(aRegion *r) {
  
  int i;
  int oldMaxAlignedSegment = r->maxAlignedSegment;

  while (r->numAlignedSegment > r->maxAlignedSegment) {
    r->maxAlignedSegment *= 2;
  }
  r->scoreAlignedSegment = (double *)/*e*/realloc(r->scoreAlignedSegment, r->maxAlignedSegment * sizeof(double));
  /* *2 since there are 2 parts that need to be aligned: a Q part and a S part */
  r->alignedSegment = (char **)/*e*/realloc(r->alignedSegment, r->maxAlignedSegment * 2 * sizeof(char *));
  /* set to NULL uninitialized segments */
  for (i = oldMaxAlignedSegment * 2; i < r->maxAlignedSegment * 2; i++) {
    r->alignedSegment[i] = NULL;
  }
}

/* close list of segments in aregion r and compute all alignments between anode-s using Needleman-Wunsch alg. */
void computeARegion(aRegion *r, Args *arguments, char *seq, element *mat, char *firstAligned, char *secondAligned, int maxCols,
                    int iQuery, int iSubject, Int64 *leftBorders, char *first, char *second, Int64 *strandBorders,
										Int64 numOfQueries, aRegion *tail, long long allSeqLength, char fastScore/*, double lambda, float *scoreMatrix,
										int alphabetSize, float ambCharPenalty*/) {

  aNode *a = r->start;
  aNode *b = NULL;
  int m, n, a_rbQ, a_rbS;
	float gapPenalty = fastScore ? arguments->G : getGapPenalty();
	ambCharPenalty = gapPenalty;
  r->totalLen = 0;

  while (a != r->end) {
    b = a->next;
    // align characters between anodes a and b (segments in anodes are exact matches): 
    // first - query positions, second - subject positions
    /* Note: right borders (e.g. a_rbQ) should be within a strand --> this should have been taken care of in constructARegionList! */
    a_rbQ = a->lbQ + a->sl - 2;
    r->totalLen += a->sl - 1;

    m = b->lbQ - (a_rbQ + 1); /* length of space between query parts of a and b; m could be 0 */
    if (m < 0) {
      m = 0;
      a->sl = a_rbQ - a->lbQ + 2;
    }
    if (m > arguments->A) {
      eprintf("[ERROR] computeARegion - Length of query and subject segments to be aligned have to be < %d!", arguments->A);
    }
    strncpy(first, seq + leftBorders[iQuery] + a_rbQ + 1, m);
    first[m] = '\0';

    /* extract subject segment into second */
    a_rbS = a->lbS + a->sl - 2;
    n = b->lbS - (a_rbS + 1); /* n could be 0 */
    if (n < 0) {
      n = 0;
      a->sl = a_rbS - a->lbS + 2;
    }
    if (n > arguments->A) {
      eprintf("[ERROR] computeARegion - Length of query and subject segments to be aligned have to be < %d!", arguments->A);
    }
    strncpy(second, seq + leftBorders[iSubject + numOfQueries] + a_rbS + 1, n); ////// ???????????????????????????
    second[n] = '\0';
    firstAligned[0] = secondAligned[0] = '\0';

		//computeMatrixNW(first, second, mat, m, n, arguments->G, arguments->M, arguments->S, maxCols,
		computeMatrixNW(first, second, mat, m, n, gapPenalty, arguments->M, arguments->S, maxCols,
			fastScore/*, alphabetSize, ambCharPenalty*/); // gap penalty, match, mismatch penalty
    tracebackBestScoreNW(first, second, firstAligned, secondAligned, mat, m, n, maxCols);    
    
    /* store new aligned segments in the list r->alignedSegment*/    
    ++ r->numAlignedSegment;
    checkMaxAlignedSegments(r); // allocate space for aligned segments
    r->scoreAlignedSegment[r->numAlignedSegment - 1] = mat[m * maxCols + n].score;
    r->score += r->scoreAlignedSegment[r->numAlignedSegment - 1];
    r->alignedSegment[r->numAlignedSegment * 2 - 2] = (char *)/*e*/realloc(r->alignedSegment[r->numAlignedSegment * 2 - 2], strlen(firstAligned) + 1);
    r->alignedSegment[r->numAlignedSegment * 2 - 1] = (char *)/*e*/realloc(r->alignedSegment[r->numAlignedSegment * 2 - 1], strlen(secondAligned) + 1);
    strcpy(r->alignedSegment[r->numAlignedSegment * 2 - 2], firstAligned);
    strcpy(r->alignedSegment[r->numAlignedSegment * 2 - 1], secondAligned);
    a = b;
    r->totalLen += strlen(firstAligned);
  } // end while
  r->totalLen += r->end->sl - 1; // add length of the last match

  r->startPosQ = r->start->lbQ;
  r->endPosQ = r->end->lbQ + r->end->sl - 2;
  r->startPosS = r->start->lbS;
  r->endPosS = r->end->lbS + r->end->sl - 2;

  // NEW: add alignment before first aNode and after aNode
  //addAlignmentBeginEnd(r, tail, arguments, seq, mat, first, second, maxCols, iQuery, iSubject, numOfQueries, leftBorders, strandBorders, firstAligned, secondAligned);      
  if (r->alignedStart[0] != NULL) {
    r->totalLen += strlen(r->alignedStart[0]);
  }
  if (r->alignedEnd[0] != NULL) {
    r->totalLen += strlen(r->alignedEnd[0]);
  }
  
  // final ratio
  r->ratio = r->score / r->totalLen; // length also includes gaps within aligned regions
  computeEValueScoreBit(r, allSeqLength/*, lambda*/);

  r->active = 0; // closed region
}

/* find minimal distance between any end subject position of p and any start subject position of q; iS = iSubject + numofQueries */
void findMinDistance(lNode *p, lNode *q, Args *arguments, Int64 *minDistance, int iS, Int64 *leftBorders, Int64 *strandBorders, aRegion *r) {

  int i, j;
  Int64 distance, p_rbS, q_rbS;
  Int64 p_lbS = -1, q_lbS = -1;

  *minDistance = -1;
  for (i = 0; i < p->sizeListPositions; i++) {
    for (j = 0; j < q->sizeListPositions; j++) {
      
      //distance = absInt64(q->listPositions[i] - (p->listPositions[j] + p->sl - 2)) - 1; /* lb ad rb are not included in distance length */
      if (q->listPositions[j] < p->listPositions[i]) { /* distance should be positive if p and q are expected to be within the same region */
        continue;
      }
      distance = q->listPositions[j] - (p->listPositions[i] + p->sl - 2) - 1; /* lb ad rb are not included in distance length; */
      if (distance < 0) { // overlapping on the subject side; e.g. when 2 copies of query fragments are matching the same or almost the same subject fragment
        distance = 0;
        //continue;
      }
      if (distance <= arguments->A && (*minDistance == -1 || distance < *minDistance)) {
        *minDistance = distance;
        p_lbS = p->listPositions[i];
        q_lbS = q->listPositions[j];
        /* check subject segments are not within the same strand */
        p_rbS = p_lbS + p->sl - 2;
        q_rbS = q_lbS + q->sl - 2;
				//if (! ((leftBorders[iS] + q_rbS <= strandBorders[iS] && leftBorders[iS] + p_rbS <= strandBorders[iS])
				if (! ((leftBorders[iS] + q_rbS < strandBorders[iS] && leftBorders[iS] + p_rbS < strandBorders[iS]) // new 04032015
              || (leftBorders[iS] + p_lbS > strandBorders[iS] && leftBorders[iS] + q_lbS > strandBorders[iS])) ) {
          *minDistance = -1; 
        }
      }
    }
  }

  // when the solution is found, i.e., *minDistance is set to some real distance, remove unused candidate positions from the list of subject positions to avoid later confusion
  if (*minDistance != -1 && (p->sizeListPositions > 1 || q->sizeListPositions > 1)) {
    p->sizeListPositions = q->sizeListPositions = 1;
    p->listPositions[0] = p_lbS;
    q->listPositions[0] = q_lbS;
    r->end->lbS = p_lbS;
  } // end if
}

/* delete last region from aRegion list, and redirect pointers */
void deleteLastARegion(aRegion **ahead, aRegion **tail, aRegion **prev, aRegion *r) {
  
  if (*tail == *ahead) { // delete r from the beginning of the list
    *tail = *ahead = NULL;
  }        
  else { // delete from the end of the list
    (*prev)->next = NULL;
    *tail = *prev;
  }
  freeARegion(r);
}

/* pointers to first elements of |Q| * |S| lists of aligned regions */
aRegion **getARHead(Int64 numOfQueries, Int64 numOfSubjects) {
	
  aRegion **arHead; // head of |Q| * |S| lists of aligned regions 
	Int64 i, j;

	arHead = (aRegion **)/*e*/malloc(numOfQueries * numOfSubjects * sizeof(aRegion *));
	/* each element in the list is a pointer to the first element of a linked list of aligned regions */
	for (i = 0; i < numOfQueries; i++) {
    for (j = 0; j < numOfSubjects; j++) {
      arHead[i * numOfSubjects + j] = NULL; // initially: all lists are empty	
    }
  }
	return arHead;
}

/* allocate space for the variables needed for the computation of nw alignment; 
 * extraLen is used when final merged segments are connected (to allocate space for additional alignment between segments); otherwise is 0
 */
void allocateElemsNW(element **mat, char **firstAligned, char **secondAligned, int *maxCols, char **first, char **second, Args *arguments, Int64 extraLen) {

  // compute maxCols
  *maxCols = 0;
  *maxCols = arguments->A > (2 * arguments->B) ? (arguments->A + 1) : (2 * arguments->B + 1);
  *maxCols = extraLen > *maxCols ? (extraLen + 1) : *maxCols;
  *mat = (element *)/*e*/malloc(*maxCols * (*maxCols) * sizeof(element));
  if (mat == NULL) {
    eprintf("[gmos]No more memory available. Please report unexpected behavior.");
  }

  /* since the maximal distance between matches should be < maxcols, then the maximal size of an aligned regions should be at most 
  2 * maxcols */
  *firstAligned = (char *)/*e*/malloc(2 * (*maxCols) + 1);
  *secondAligned = (char *)/*e*/malloc(2 * (*maxCols) + 1);
  *first = (char *)/*e*/malloc(*maxCols + 1);
  *second = (char *)/*e*/malloc(*maxCols + 1);
}

/* free elements needed for the computation of nw alignments */
void freeElemsNW(element *mat, char *firstAligned, char *secondAligned, char *first, char *second) {
  
  freeMatrix(mat);
  free(first);
  free(second);
  free(firstAligned);
  free(secondAligned);

}

/* construct all lists of aligned regions for each combination (Qi, Sj) */
void constructAlignment(char *seq, lNode **head, Args *arguments, aRegion ***ahead, Int64 numOfQueries, Int64 numOfSubjects, 
                        Int64 *leftBorders, Int64 *strandBorders, Int64 *seqBorders, long long allSeqLength, 
												double *lambdaAll, float **scoreMatrixAll, int alphabetSize, float ambCharPenaltyAll, 
												char fastScore, double *arrayK) {

	Int64 i, j;
  element *mat = NULL;
  char *firstAligned = NULL;
  char *secondAligned = NULL;
  char *first = NULL;
  char *second = NULL;
  int maxCols;

  allocateElemsNW(&mat, &firstAligned, &secondAligned, &maxCols, &first, &second, arguments, 0);

  for (i = 0; i < numOfQueries; i++) {
    // set query specific lambda and scoreMatrix
		setQueryLambdaMatrix(lambdaAll, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, i, arrayK);

    for (j = 0; j < numOfSubjects; j++) {
      /* construct a single list of aligned regions for (Qi, Sj) */
      // wrong: mat is of type elem!! memset(mat, 0, maxCols * maxCols * sizeof(double));
      setMatElemsToZero (mat, maxCols);
      if (head[i * numOfSubjects + j]) {
        // constructARegionList(head[i * numOfSubjects + j], &((*ahead)[i * numOfSubjects + j]), arguments, seq, mat,
                          // firstAligned, secondAligned, maxCols, (int)i, (int)j, leftBorders, first, second, strandBorders, seqBorders, numOfQueries, allSeqLength);
				constructARegionList2(head[i * numOfSubjects + j], &((*ahead)[i * numOfSubjects + j]), arguments, seq, mat,
					firstAligned, secondAligned, maxCols, (int)i, (int)j, leftBorders, first, second,
					strandBorders, seqBorders, numOfQueries, allSeqLength, fastScore);
													/*, lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty*/
      }
      /* possibly: aRegionList may be empty for Sj, then set pointer to the list to null */
      if ((*ahead)[i * numOfSubjects + j] == NULL) {
        // do nothing
      }
      freeList(head[i * numOfSubjects + j]);
    }
  }
  
  //free(mat);
  free(head);
  freeElemsNW(mat, firstAligned, secondAligned, first, second);
}


Int64 absInt64(Int64 a) {
  if (a >= 0) {
    return a;
  }
  else {
    return -a;
  }
}

/* print complete aregion list */
void printARegionList(aRegion **arList, Int64 numOfQueries, Int64 numOfSubjects, FILE *fpout, short argB, char **seqHeaders, Int64 *strandLen) {

  Int64 i, j;
  int k;
  aRegion *r;
  aNode *a;
  char strainFR = '+'; // forward straind: '+'; reverse strand: '-'
  Int64 subjectStart = 0, subjectEnd = 0;	
  Int64 a_lbS = 0, a_rbS = 0;

  for (i = 0; i < numOfQueries; i++) {
    for (j = 0; j < numOfSubjects; j++) {
      r = arList[i * numOfSubjects + j];
      fprintf(fpout, "\n------------- Regions' list - query: %lld (%s), subject: %lld (%s) ---------------------------\n", 
              (long long)i + 1, seqHeaders[i] + 1, (long long)j + 1, seqHeaders[numOfQueries + j] + 1);
      fflush(fpout);
      while (r && r->endPosQ <= strandLen[i]) { // forward strand
        a = r->start;
		    strainFR = (r->startPosS >= strandLen[numOfQueries + j]) ? '-' : '+';
        getSubjectStartEndFR(&subjectStart, &subjectEnd, strandLen, j, numOfQueries, r, strainFR);
		    //if (strainFR == '-') { // reverse strand --> compute corresponding positions on fwd strand
		    //  subjectStart = 2 * strandLen[numOfQueries + j] - r->startPosS;
		    //  subjectEnd = 2 * strandLen[numOfQueries + j] - r->endPosS;
	     // }
	     // else {
		    //  subjectStart = r->startPosS;
		    //  subjectEnd = r->endPosS;	  
	     // }
        fprintf(fpout, "Query: %lld - %lld Subject: %lld - %lld (%c) Score: %.1lf\n", (long long)r->startPosQ + 1, (long long)r->endPosQ + 1, 
          (long long)subjectStart + 1, (long long)subjectEnd + 1, strainFR, r->score);			
        fflush(fpout);
        if (argB && r->alignedStart[0]) {
          // fprintf(fpout, "Q*: %s S*: %s\n", r->alignedStart[0], r->alignedStart[1]);	
          // fflush(fpout);
        }
        k = 0;
        while (a) {
          fprintf(fpout, "\t%lld %lld %lld\t", (long long)a->lbQ + 1, (long long)a->lbQ + a->sl + 1 - 2, (long long)a->sl - 1);		
          fflush(fpout);
          fprintf(fpout, "\tSubject: ");
          fflush(fpout);
		      a_lbS = (strainFR == '-') ? (2 * strandLen[numOfQueries + j] - a->lbS) : a->lbS;
		      a_rbS = (strainFR == '-') ? (2 * strandLen[numOfQueries + j] - (a->lbS + a->sl - 2)) : a->lbS + a->sl - 2;
		      //fprintf(fpout, "\t%lld - %lld\n", (long long)a->lbS + 1, (long long)a->lbS + a->sl + 1 - 2);
		      fprintf(fpout, "\t%lld - %lld\n", (long long)a_lbS + 1, (long long)a_rbS + 1);		  
          fflush(fpout);
          if (k < r->numAlignedSegment * 2) {
            fprintf(fpout, "Q: %s S: %s\n", r->alignedSegment[k], r->alignedSegment[k + 1]);			
            fflush(fpout);
            k += 2;
          }
          a = a->next;
        } // end while
        if (argB && r->alignedStart[0]) {
          // fprintf(fpout, "Q*: %s S*: %s\n", r->alignedEnd[0], r->alignedEnd[1]);			
          // fflush(fpout);
        }
        fprintf(fpout, "\n");
        fflush(fpout);
        r = r->next;
      } // end while      
    }
  }
}


/* get beginning position of an aligned segment just before r */
Int64 getBeginPos(Int64 beginRPos, Int64 tailEndPos, aRegion *tail, Int64 afterTailLen, Int64 argB, Int64 *leftBorders, Int64 *strandBorders, int i, aRegion *r) {
  
  Int64 beginPos = -1, lastTailPos = -1;
  
  beginRPos += leftBorders[i]; // absolute position in concatenation of all sequences
  beginPos = beginRPos - argB;

  if (tail && tail != r) { /* if there is a predecessor of r */
    lastTailPos = tailEndPos + afterTailLen + leftBorders[i];
  }
  if (lastTailPos >= beginPos) {
    beginPos = lastTailPos + 1;
  }

  if (beginPos < leftBorders[i]) { /* check sequence start */
    beginPos = leftBorders[i];
  }
	else if (beginPos <= strandBorders[i] && beginRPos > strandBorders[i]) { /* check start of the reverse strand */
	//else if (beginPos < strandBorders[i] && beginRPos > strandBorders[i]) { /* check start of the reverse strand */ // new 04032015
    beginPos = strandBorders[i] + 1;
  }

  return beginPos - leftBorders[i]; // relative position within i-th sequence
}


/* get end position of an aligned segment just after r */
Int64 getEndPos(Int64 endRPos, Int64 argB, Int64 *leftBorders, Int64 *strandBorders, int i) {

  Int64 endPos = -1;

  // compute end position of an alignment just after r; 
  // check strand or sequence end (the begininning of the next aregion is not checked, since the next aregion does not exist at this point)
  endPos = endRPos + argB + leftBorders[i]; // absolute position in concatenation of all sequences
  endRPos += leftBorders[i];
  if (endPos > strandBorders[i] && endRPos <= strandBorders[i]) { // strand end
    endPos = strandBorders[i];
  }
  else if (endPos > leftBorders[i] + 2 * (strandBorders[i] - leftBorders[i]) - 1) { // sequence end
    endPos = leftBorders[i] + 2 * (strandBorders[i] - leftBorders[i]) - 1;
  }
  return endPos - leftBorders[i]; // relative position within i-th sequence
}


/* align segments before the beginning of r and after the end of r */
void addAlignmentBeginEnd(aRegion *r, aRegion *tail, Args *arguments, char *seq, element *mat, char *first, char *second, int maxCols, 
                          int iQuery, int iSubject, Int64 numOfQueries, Int64 *leftBorders, Int64 *strandBorders, 
													char *firstAligned, char *secondAligned, char fastScore
													/*, float *scoreMatrix, int alphabetSize, float ambCharPenalty*/) {

  Int64 beginPosQ = -1, endPosQ = -1;
  Int64 beginPosS = -1, endPosS = -1;
  Int64 tailAlignedEndLenQ = 0, tailAlignedEndLenS = 0;
  int m = 0, n = 0;
	float gapPenalty = fastScore ? arguments->G : getGapPenalty();
	ambCharPenalty = gapPenalty;

  /* align segments before the beginning of r, and after the end of r;
  * at most a->B positions should be aligned (take care of the beginning and the end of both subject and query sequence!)
  */
  
  /* (1) retrieve subsequences of both subject and query which are just before r */
  if (tail && tail->alignedEnd[0]) {
    tailAlignedEndLenQ = strlen(tail->alignedEnd[0]);
    tailAlignedEndLenS = strlen(tail->alignedEnd[1]);
  }
  beginPosQ = getBeginPos(r->start->lbQ, tail->end->lbQ + tail->end->sl - 2, tail, tailAlignedEndLenQ, arguments->B, leftBorders, strandBorders, iQuery, r);
  beginPosS = getBeginPos(r->start->lbS, tail->end->lbS + tail->end->sl - 2, tail, tailAlignedEndLenS, arguments->B, leftBorders, strandBorders, numOfQueries + iSubject, r);
  
  m = r->start->lbQ - beginPosQ;
  if (m > arguments->B) {
    m = arguments->B;
  }
  n = r->start->lbS - beginPosS;
  if (n > arguments->B) {
    n = arguments->B;
  }
  strncpy(first, seq + leftBorders[iQuery] + beginPosQ, m); // query
  strncpy(second, seq + leftBorders[numOfQueries + iSubject] + beginPosS, n); // subject
  first[m] = '\0';
  second[n] = '\0'; 
  firstAligned[0] = secondAligned[0] = '\0';
	//computeMatrixNW(first, second, mat, m, n, arguments->G, arguments->M, arguments->S, maxCols, fastScore); // gap penalty, match, mismatch penalty
	computeMatrixNW(first, second, mat, m, n, gapPenalty, arguments->M, arguments->S, maxCols, fastScore); 
  tracebackBestScoreNW(first, second, firstAligned, secondAligned, mat, m, n, maxCols);    
    
  /* store new aligned segments in r->alignedStart */    
  r->scoreStart = mat[m * maxCols + n].score;
  if (r->scoreStart < 0) {
    r->scoreStart = 0;
    strcpy(firstAligned, "");
    strcpy(secondAligned, "");
  }
  r->score += r->scoreStart;  
  r->alignedStart[0] = (char *)/*e*/malloc(strlen(firstAligned) + 1);
  r->alignedStart[1] = (char *)/*e*/malloc(strlen(secondAligned) + 1);
  strcpy(r->alignedStart[0], firstAligned);
  strcpy(r->alignedStart[1], secondAligned);

  /* (2) retrieve subsequences of both subject and query which are just after r */
  endPosQ = getEndPos(r->end->lbQ + r->end->sl - 2, arguments->B, leftBorders, strandBorders, iQuery);
  endPosS = getEndPos(r->end->lbS + r->end->sl - 2, arguments->B, leftBorders, strandBorders, iSubject + numOfQueries);
  
  m = endPosQ - (r->end->lbQ + r->end->sl - 2);
  if (m > arguments->B) {
    m = arguments->B;
  }
  else if (m < 0) {
    m = 0;
    r->end->sl = endPosQ - r->end->lbQ + 2;
  }

  n = endPosS - (r->end->lbS + r->end->sl - 2);
  if (n > arguments->B) {
    n = arguments->B;
  }
  else if (n < 0) {
    n = 0;
    r->end->sl = endPosS - r->end->lbS + 2;
  }
  strncpy(first, seq + leftBorders[iQuery] + r->end->lbQ + r->end->sl - 1, m); // query
  strncpy(second, seq + leftBorders[iSubject + numOfQueries] + r->end->lbS + r->end->sl - 1, n); // subject
  first[m] = '\0';
  second[n] = '\0'; 
  firstAligned[0] = secondAligned[0] = '\0';
	//computeMatrixNW(first, second, mat, m, n, arguments->G, arguments->M, arguments->S, maxCols, fastScore);
		//scoreMatrix, alphabetSize, ambCharPenalty); // gap penalty, match, mismatch penalty
	computeMatrixNW(first, second, mat, m, n, gapPenalty, arguments->M, arguments->S, maxCols, fastScore);
	tracebackBestScoreNW(first, second, firstAligned, secondAligned, mat, m, n, maxCols);    
    
  /* store new aligned segments in r->alignedEnd */    
  r->scoreEnd = mat[m * maxCols + n].score;
  if (r->scoreEnd < 0) {
    r->scoreEnd = 0;
    strcpy(firstAligned, "");
    strcpy(secondAligned, "");
  }
  r->score += r->scoreEnd;  
  r->alignedEnd[0] = (char *)/*e*/malloc(strlen(firstAligned) + 1);
  r->alignedEnd[1] = (char *)/*e*/malloc(strlen(secondAligned) + 1);
  strcpy(r->alignedEnd[0], firstAligned);
  strcpy(r->alignedEnd[1], secondAligned);

  /* final positions of an aligned query and subject region, when segments just before and after start-anode and end-anode are added */
  //r->startPosQ = r->start->lbQ - strlen(r->alignedStart[0]);
  //r->endPosQ = r->end->lbQ + r->end->sl - 2 + strlen(r->alignedEnd[0]);
  //r->startPosS = r->start->lbS - strlen(r->alignedStart[1]);
  //r->endPosS = r->end->lbS + r->end->sl - 2 + strlen(r->alignedEnd[1]);
  // The upper code was wrong --> since strlen(r->alignedEnd[1]) could have been longer than args->B, although only args->B nucleotides were involved in computation
  r->startPosQ = (strlen(r->alignedStart[0]) > (unsigned int)arguments->B) ? (r->start->lbQ - arguments->B) : (r->start->lbQ - strlen(r->alignedStart[0]));
	r->endPosQ = (strlen(r->alignedEnd[0]) > (unsigned int)arguments->B) ? (r->end->lbQ + r->end->sl - 2 + arguments->B) : (r->end->lbQ + r->end->sl - 2 + strlen(r->alignedEnd[0]));
	r->startPosS = (strlen(r->alignedStart[1]) > (unsigned int)arguments->B) ? (r->start->lbS - arguments->B) : (r->start->lbS - strlen(r->alignedStart[1]));
	r->endPosS = (strlen(r->alignedEnd[1]) > (unsigned int)arguments->B) ? (r->end->lbS + r->end->sl - 2 + arguments->B) : (r->end->lbS + r->end->sl - 2 + strlen(r->alignedEnd[1]));
}

/* align a segment from the start position in the i-th sequence (length m) to the start position of the j-th sequence (length n) */
void alignSegment(Int64 startPosISeq, Int64 startPosJSeq, Int64 *leftBorders, Int64 m, Int64 n,  char *first, char *second, char *firstAligned, char *secondAligned, char *seq
	, Int64 iSeq, Int64 jSeq, element *mat, Args *arguments, int maxCols, char fastScore) {
	//, float *scoreMatrix, int alphabetSize, float ambCharPenalty) {
	int i;
	float gapPenalty = fastScore ? arguments->G : getGapPenalty();
	ambCharPenalty = gapPenalty;

  assert(m >= 0);
  assert(n >= 0);
  if (m >= 0 && n >= 0) {
		if (m) {
			strncpy(first, seq + leftBorders[iSeq] + startPosISeq, m); 
		}
		if (n) {
			strncpy(second, seq + leftBorders[jSeq] + startPosJSeq, n); 
		}
		first[m] = second[n] = '\0';
		firstAligned[0] = secondAligned[0] = '\0';
		if (m && n) {
			//computeMatrixNW(first, second, mat, m, n, arguments->G, arguments->M, arguments->S, maxCols, fastScore);
				//scoreMatrix, alphabetSize, ambCharPenalty); // gap penalty, match, mismatch penalty
			computeMatrixNW(first, second, mat, m, n, gapPenalty, arguments->M, arguments->S, maxCols, fastScore);
			tracebackBestScoreNW(first, second, firstAligned, secondAligned, mat, m, n, maxCols);    
		}
		else if (m) {
			strcpy(firstAligned, first);
			// new 31/01/2015
			for (i = 0; i < m; i++) {
				secondAligned[i] = '-';
			}
			secondAligned[m] = '\0';
			///
			mat->score = strlen(first) * gapPenalty;
		}
		else if (n) {
			strcpy(secondAligned, second);
			// new 31/01/2015
			for (i = 0; i < n; i++) {
				firstAligned[i] = '-';
			}
			firstAligned[n] = '\0';
			///
			mat->score = strlen(second) * gapPenalty;
		}
  }
	else {
		fprintf(stderr, "Error alignSegment m or n < 0: m = %lld n = %lld, i = %lld j = %lld\n", (long long)m, (long long)n, (long long)leftBorders[iSeq] + startPosISeq, 
				(long long)leftBorders[jSeq] + startPosJSeq);
    fflush(stderr);
	}  

}

/* copy aligned segments to aRegion's elements alignedSeg1 and alignedSeg2 and update appropriate s
cores */
void storeAlignedSegments(double *newScore, double score, double *totalScore, char *firstAligned, char *secondAligned, char **alignedSeg1, char **alignedSeg2) {

  *newScore = score;
  *totalScore += score;  
  if (*newScore <= 0) {
    *newScore = 0;
    //(*alignedSeg1)[0] = (*alignedSeg2)[0] = '\0';
    if (*alignedSeg1 && *alignedSeg2) {
      strcpy(*alignedSeg1, "");
      strcpy(*alignedSeg2, "");
    }
  }
  else {
    *alignedSeg1 = (char *)/*e*/realloc(*alignedSeg1, strlen(firstAligned) + 1);
    *alignedSeg2 = (char *)/*e*/realloc(*alignedSeg2, strlen(secondAligned) + 1);
    strcpy(*alignedSeg1, firstAligned);
    strcpy(*alignedSeg2, secondAligned);
  }
}

// connect r2 to tail->end, re-compute alignment between tail->end and r2->start, and then delete r2
void connectARegion(aRegion *tail, aRegion *r2, Args *arguments, char *seq, Int64 *leftBorders, int iSubject, int iQuery, Int64 numOfQueries,
                    char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat, 
										long long allSeqLength, char fastScore) {
										//, double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty, char fastScore) {

  int m, n;  
  double mat_score = 0.0;

  /* align region between tail and r2 */
  m = r2->start->lbQ - (tail->end->lbQ + tail->end->sl - 1);
  n = r2->start->lbS - (tail->end->lbS + tail->end->sl - 1);
  alignSegment(tail->end->lbQ + tail->end->sl - 1, tail->end->lbS + tail->end->sl - 1, leftBorders, m, n,  first, second, firstAligned, secondAligned, seq
                  , iQuery, iSubject + numOfQueries, mat, arguments, maxCols/*, scoreMatrix, alphabetSize, ambCharPenalty*/, fastScore);
  if (m && n) {
    mat_score = mat[m * maxCols + n].score;
  }
  else {
    mat_score = mat->score;
  }

  /* new totalLen and score */
  tail->totalLen = tail->totalLen + r2->totalLen + strlen(firstAligned); // subtract tail->after-end and r2->before-start
  if (r2->alignedStart[0] != NULL) {
    tail->totalLen -= strlen(r2->alignedStart[0]);
  }
  if (tail->alignedEnd[0] != NULL) {
    tail->totalLen -= strlen(tail->alignedEnd[0]);
  }

  /* add aligned new region to tail */
  ++ tail->numAlignedSegment;
  checkMaxAlignedSegments(tail); // allocate space for aligned segment beetween tail and r2
  tail->score = tail->score - tail->scoreEnd + r2->score - r2->scoreStart;  
  //storeAlignedSegments(&tail->scoreAlignedSegment[tail->numAlignedSegment - 1], mat[m * maxCols + n].score, &tail->score, firstAligned, secondAligned, 
  //  &tail->alignedSegment[tail->numAlignedSegment * 2 - 2], &tail->alignedSegment[tail->numAlignedSegment * 2 - 1]);
  storeAlignedSegments(&tail->scoreAlignedSegment[tail->numAlignedSegment - 1], mat_score, &tail->score, firstAligned, secondAligned, 
    &tail->alignedSegment[tail->numAlignedSegment * 2 - 2], &tail->alignedSegment[tail->numAlignedSegment * 2 - 1]);

  tail->ratio = tail->score / tail->totalLen;
 
  /* add r2 aligned segments and their scores to tail and re-compute numAlignedSegment and maxAlignedSegment*/
  addAlignedSegments(tail, r2);
 
  /* connect tail and r2, and tail to r2's successor */
  tail->end->next = r2->start;
  tail->end = r2->end;
  tail->next = r2->next;
  
  // new end after the last node
  free(tail->alignedEnd[0]);
  free(tail->alignedEnd[1]);
  tail->alignedEnd[0] = r2->alignedEnd[0];
  tail->alignedEnd[1] = r2->alignedEnd[1];
  tail->scoreEnd = r2->scoreEnd;

  tail->endPosQ = r2->endPosQ;
  tail->endPosS = r2->endPosS;
  tail->ratio = tail->score / tail->totalLen;
  computeEValueScoreBit(tail, allSeqLength/*, lambda*/);
  freeARegionPartly(r2);

}

/* add aligned segments of r2 to tail */
void addAlignedSegments(aRegion *tail, aRegion *r2) {
  
  int i;
  int oldNumAlignedSegment = tail->numAlignedSegment;
  
  tail->numAlignedSegment += r2->numAlignedSegment;
  checkMaxAlignedSegments(tail); // allocate space for r2 aligned segments and their scores

  for (i = 1; i <= r2->numAlignedSegment; i++) {
    /* to avoid double allocating, i.e. allocating nodes for r2 and then deleting them, all aligned segments are just added to tail->alignedSegment list */
    tail->alignedSegment[(oldNumAlignedSegment + i) * 2 - 2] = r2->alignedSegment[i * 2 - 2];
    tail->alignedSegment[(oldNumAlignedSegment + i) * 2 - 1] = r2->alignedSegment[i * 2 - 1];
    tail->scoreAlignedSegment[oldNumAlignedSegment + i - 1] = r2->scoreAlignedSegment[i - 1];
  }
}

void freeARegionPartly(aRegion *r) {
  
  int i;  
  free(r->alignedSegment);
  free(r->scoreAlignedSegment);
  
  //// NEW
  for (i = 0; i < 2; i++) {
    free(r->alignedStart[i]);
  }
  free(r->alignedStart);  
  free(r->alignedEnd);  
  free(r);
  r = NULL;
}

/* compute S bit score and E-value for a segment n; 
 * References: http://www.ncbi.nlm.nih.gov/BLAST/tutorial/ */
void computeEValueScoreBit(aRegion *r, long long allSeqLength/*, double lambda*/) {
  //sdouble p;
	//r->score_blast = (DEFAULT_LAMBDA * (double)r->score - log(DEFAULT_K)) / log(2.0);s
	//r->score_blast = (lambda * (double)r->score - log(DEFAULT_K)) / log(2.0);
	r->score_blast = (lambda * (double)r->score - log(blastK)) / log(2.0);
  r->e_value = (double)r->totalLen * (double)allSeqLength * pow(2.0, -r->score_blast);
}

/* copy aligned segments to aRegion's elements alignedSeg1 and alignedSeg2 and update appropriate scores */
void storeAlignedSegments2(double *newScore, double score, double *totalScore, char *firstAligned, char *secondAligned, char **alignedSeg1, char **alignedSeg2) {

  *newScore = score;
  *totalScore += score;  
  *alignedSeg1 = (char *)/*e*/realloc(*alignedSeg1, strlen(firstAligned) + 1);
  *alignedSeg2 = (char *)/*e*/realloc(*alignedSeg2, strlen(secondAligned) + 1);
  strcpy(*alignedSeg1, firstAligned);
  strcpy(*alignedSeg2, secondAligned);
}


  
// connect r2 to tail->end, re-compute alignment between tail->end and r2->start, and then delete r2 (tail stays)
int connectARegion2(aRegion *tail, aRegion *r2, Args *arguments, char *seq, Int64 *leftBorders, int iSubject, int iQuery, Int64 numOfQueries,
                    char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat, 
										long long allSeqLength, char fastScore) {
										//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty) {

  int m, n;  
  double mat_score = 0.0;

  /* align region between tail and r2; tail and r2 do not have start and end */
  m = r2->start->lbQ - (tail->end->lbQ + tail->end->sl - 1);
  n = r2->start->lbS - (tail->end->lbS + tail->end->sl - 1);
	alignSegment(tail->end->lbQ + tail->end->sl - 1, tail->end->lbS + tail->end->sl - 1, leftBorders, m, n, first, second, firstAligned, secondAligned, seq
		, iQuery, iSubject + numOfQueries, mat, arguments, maxCols, fastScore);
									//scoreMatrix, alphabetSize, ambCharPenalty);
  if (m && n) {
    mat_score = mat[m * maxCols + n].score;
  }
  else {
    mat_score = mat->score;
  }
	if (mat_score < 0) {
		return 0;
	}
  /* new totalLen and score */
  //tail->totalLen = tail->totalLen + r2->totalLen + strlen(firstAligned); // subtract tail->after-end and r2->before-start
  tail->totalLen = r2->end->lbQ + r2->end->sl - 2 - tail->start->lbQ + 1;
  tail->score = tail->score + r2->score;

  /* add aligned new region to tail and update score */
  ++ tail->numAlignedSegment;
  checkMaxAlignedSegments(tail); // allocate space for aligned segment beetween tail and r2
  //storeAlignedSegments2(&tail->scoreAlignedSegment[tail->numAlignedSegment - 1], mat[m * maxCols + n].score, &tail->score, firstAligned, secondAligned, 
  //  &tail->alignedSegment[tail->numAlignedSegment * 2 - 2], &tail->alignedSegment[tail->numAlignedSegment * 2 - 1]);  
  storeAlignedSegments2(&tail->scoreAlignedSegment[tail->numAlignedSegment - 1], mat_score, &tail->score, firstAligned, secondAligned, 
    &tail->alignedSegment[tail->numAlignedSegment * 2 - 2], &tail->alignedSegment[tail->numAlignedSegment * 2 - 1]);  
 
  /* add r2 aligned segments and their scores to tail and re-compute numAlignedSegment and maxAlignedSegment*/
  addAlignedSegments(tail, r2);
 
  // free old aligned ends
  free(tail->alignedEnd[0]);
  free(tail->alignedEnd[1]);
  tail->alignedEnd[0] = r2->alignedEnd[0];
  tail->alignedEnd[1] = r2->alignedEnd[1];
  tail->scoreEnd = r2->scoreEnd;

  // free old aligned r2->start, if any --> freeARegionPartly

  /* connect tail and r2, and tail to r2's successor */
  tail->end->next = r2->start;
  tail->end = r2->end;
  tail->next = r2->next;
  
  // new end after the last node
  tail->endPosQ = r2->endPosQ;
  tail->endPosS = r2->endPosS;
  tail->ratio = tail->score / tail->totalLen;
  if (tail->score < 0) {
    tail->score = tail->score;
  }
  computeEValueScoreBit(tail, allSeqLength/*, lambda*/);
  freeARegionPartly(r2);
	return 1;
}

void initializeHeadTail(lNode *p, aRegion **ahead, aRegion **tail, aRegion **activeFirst, aRegion **prev, aRegion **prevActiveFirst, Args *arguments) {
  
  aNode *a = NULL;
  aRegion *r = NULL;
  int i;
  a = getANode(p); // always takes listPositions[0]
  r = getARegion(a, arguments->M);
  if (*ahead == NULL) { 
    *ahead = *tail = r;
    *prevActiveFirst = *prev = NULL;
  }
  else if (*tail == NULL) { // kada se *tail obriše, jer je bio zadnji *activeFirst, koji se onda obrisao
    if (*prev) { // ??
      (*prev)->next = r;
    }
    if (*prevActiveFirst) { // ??
      (*prevActiveFirst)->next = r;
    }
    *tail = r;
  }
  *activeFirst = r;
	for (i = 1; i < p->sizeListPositions; i++) { // typically sizeListPositions is 1, rarely 2 and almost never > 2
	  a = getANodeNextPosition(p, i);
	  r = getARegion(a, arguments->M); // each new aNode is added to a new aRegion
    (*tail)->next = r; /* connect to existing aregion-s in the list */
	  *tail = r;
	}
}


// check that at least one q->listPositions[i] is less than A positions apart from the current end on the subject side of rCurr
int chkSubjectList(aRegion *rCurr, lNode *q, short A) {
	int i;
	for (i = 0; i < q->sizeListPositions; i++) {
		if (q->listPositions[i] - rCurr->endPosS - 1 <= A) {
			return 1;
		}
	}
	return 0;
}

/* construct a single list of aligned regions for (Qi, Sj); function is called only if head != NULL */
void constructARegionList2(lNode *head, aRegion **ahead, Args *arguments, char *seq, element *mat,
                          char *firstAligned, char *secondAligned, int maxCols, int iQuery, int iSubject, Int64 *leftBorders,
                          char *first, char *second, Int64 *strandBorders, Int64 *seqBorders, Int64 numOfQueries, 
													long long allSeqLength, char fastScore
													/*, double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty*/) {
  
  lNode *q = NULL, *p = NULL;
  int j;
  aNode *a = NULL;
  aRegion *tail = NULL, *r2 = NULL, *prev = NULL;
  aRegion *activeFirst = NULL; //, *activeLast = NULL; // first and last aRegion-s that should be checked
  aRegion *pActiveFirst = NULL; //, *activeLast = NULL; // first and last aRegion-s that should be checked
  aRegion *ppActiveFirst = NULL; //, *activeLast = NULL; // first and last aRegion-s that should be checked
  aRegion *prevActiveFirst = NULL; 
  aRegion *rCurr = NULL;
  int cnt = 0, leaveRCurr = 0;

  p = head;
  // first element in this list of aligned regions
  if (*ahead == NULL) { 
    initializeHeadTail(p, ahead, &tail, &activeFirst, &prev, &prevActiveFirst, arguments);
  }

  q = p->next; // next segment
  // Note: all aNode-s (p, q, ...) are represented using relative coordinates for fwd and rev strand, i.e. even on rev, strand lbStart < lbEnd
  //prevActiveFirst = NULL; // last closed region preceding activeFirst
  while (q && q->rb + leftBorders[iQuery] < strandBorders[iQuery]) { // at this moment: only for forward strand    
		// check whether a will be added to the end of some existing active aRegion or to the beginning of a new one
    rCurr = activeFirst; 
    prev = prevActiveFirst; // prev is previous to rCurr
    cnt = 0; // counting all subject regions in listPositions
    while (rCurr) { // if activeFirst is null, this loop is skipped
      if (rCurr->active) { // not closed aRegion
        if ((q->lb - rCurr->endPosQ - 1) > arguments->A /*||
//////////////// critical line - if uncommented, all lower signals are ignored /////////////
						!chkSubjectList(rCurr, q, arguments->A)*/) { // close rCurr (a is too far from rCurr, and so are all subsequent q-s)          
          // close rCurr --> computeARegion sets active to 0
          computeARegion(rCurr, arguments, seq, mat, firstAligned, secondAligned, maxCols, iQuery, iSubject, leftBorders, 
						first, second, strandBorders, numOfQueries, tail, allSeqLength, fastScore
						/*, lambda, scoreMatrix, alphabetSize, ambCharPenalty*/);
          leaveRCurr = checkARegion(rCurr, arguments); // delete rCurr if shorter than f, has ratio less than R, etc; otherwise, leave it

          // find new prevActiveFirst i activeFirst; note: prevActiveFirst is just a pointer to a predecessor of activeFirst; 
          // doesn't make any connections using ->next (it's done using prev) ??
          if (rCurr == activeFirst) { // if rCurr is activeFirst, then find new activeFirst            
            pActiveFirst = activeFirst->next; // find next activeFirst: pActiveFirst
            while (pActiveFirst && pActiveFirst->active == 0) {
              ppActiveFirst = pActiveFirst;
              pActiveFirst = pActiveFirst->next;
            }
            // in the following two cases, it doesn't matter whether activeFirst will be deleted or not
            if (pActiveFirst && pActiveFirst != activeFirst->next) { 
              prevActiveFirst = ppActiveFirst;
            }
            else if (!pActiveFirst) {
              if (leaveRCurr == 0 && prevActiveFirst) {
                prevActiveFirst->next = activeFirst->next; // ????
              }
              else {
                //prevActiveFirst = NULL; do nothing
              }
            }
						// when rCurr is not deleted
            else if (leaveRCurr && pActiveFirst && pActiveFirst == activeFirst->next) { 
							prevActiveFirst = activeFirst;  
            }
            activeFirst = pActiveFirst;
          } // end if rCurr == activeFirst    	    

          if (leaveRCurr == 0) { // delete rCurr if shorter than f, has ratio less than R, etc.
            r2 = rCurr->next;
            deleteCurrARegion(ahead, &tail, &prev, rCurr);
            rCurr = r2;
            if (rCurr && activeFirst && rCurr->startPosQ < activeFirst->startPosQ) { // if there are aRegions between rCurr and activeFirst which are all closed, skip them
              prev = prevActiveFirst;
              rCurr = activeFirst;
            }
            continue;
          }
        }
        else { // if (a->lbQ - rCurr->endPosQ - 1 <= arguments->A) {
					for (j = 0; j < q->sizeListPositions; j++) { // for each subject part from listPositions
            if (q->listPositions[j] > 0) {
							// new 04032015
							if ((q->listPositions[j] + leftBorders[iSubject + numOfQueries] < strandBorders[iSubject + numOfQueries] &&
								q->listPositions[j] + q->sl + leftBorders[iSubject + numOfQueries] > strandBorders[iSubject + numOfQueries]) ||
								(q->listPositions[j] + q->sl > 2 * (strandBorders[iSubject + numOfQueries] - leftBorders[iSubject + numOfQueries] + 1) - 1)) {
								continue;
							}
							else if (q->listPositions[j] + leftBorders[iSubject + numOfQueries] > strandBorders[iSubject + numOfQueries] &&
								rCurr->endPosS + leftBorders[iSubject + numOfQueries] < strandBorders[iSubject + numOfQueries]) {
								continue;
							}
							else if (q->listPositions[j] + leftBorders[iSubject + numOfQueries] < strandBorders[iSubject + numOfQueries] &&
								rCurr->endPosS + leftBorders[iSubject + numOfQueries] > strandBorders[iSubject + numOfQueries]) {
								continue;
							}
							// new 04032015 - end

							// new - changed 21/01/2015
              if (((q->listPositions[j] - rCurr->endPosS - 1) >= 0 && (q->listPositions[j] - rCurr->endPosS - 1) <= arguments->A) // add aNode to the end of rCurr, if they are close enough
                || ((q->listPositions[j] - rCurr->endPosS - 1) < 0 && fabs(q->listPositions[j] - rCurr->endPosS - 1) <= LEN_OVERLAP)) { // add aNode to the end of rCurr, if they are close enough
                // overlapping of two exact matches over a very short region, typically <= 10 bp; in this case, shorten new node from left
                if ((q->listPositions[j] - rCurr->endPosS - 1) < 0 && fabs(q->listPositions[j] - rCurr->endPosS - 1) <= LEN_OVERLAP) { // add aNode to the end of rCurr, if they are close enough
                  // this happens e.g. when there is a very short gap (1 bp or so) between long matches
                  int len = (rCurr->end->lbQ + rCurr->end->sl - 2) - q->lb + 1; // end - q->start + 1
									///////////// new 14022015 ////////////
									// if the distance between subject and query would now increase > A, then do not add this subject position
									if (q->listPositions[j] + len - rCurr->endPosS - 1 > arguments->A) {
										continue;
									}
									//////////////////// end //////////////
									q->lb = rCurr->end->lbQ + rCurr->end->sl - 1;
                  q->sl -= len; //= q->rb - q->lb + 2;
                  q->listPositions[j] += len;
								}
                //// if a gap should be added on the subject or query side
                if (q->listPositions[j] == rCurr->endPosS || q->lb == rCurr->endPosQ) {
									///////////// new 14022015 ////////////
									// if the distance between subject and query would now increase > A, then do not add this subject position
									if ((q->listPositions[j] + 1 - rCurr->endPosS - 1 > arguments->A) || 
										(q->lb + 1 - rCurr->endPosQ - 1 > arguments->A)) {
										continue;
									}
									//////////////////// end //////////////
									q->lb += 1;
                  q->sl -= 1; //= q->rb - q->lb + 2;
                  q->listPositions[j] += 1;
                }
                /////////// end - 22/01/2015
                a = getANodeNextPosition(q, j);
                addANode(rCurr, a, arguments);
                q->listPositions[j] = -(q->listPositions[j]); // setting to negative value means that q->listPositions[j] was added; otherwise it needs to be added in a new aRegion 
                ++cnt;
								break; // ???????? added 03/02/2015
              }
            }
          } // end for
          if (cnt == q->sizeListPositions) {
            break;
          }
        } // end else-if
      } // end if
      prev = rCurr;
      rCurr = rCurr->next;
    } // end while - over rCurr from activeFirst to activeLast
    
    // // add remaining matches stored in q->sizeListPositions to new aRegions 
    if (cnt < q->sizeListPositions) {
      if (tail) { // if there is something in the list
	      for (j = 0; j < q->sizeListPositions; j++) { // for each subject part from listPositions check whether it was already added to rCurr, if not, add it to a new aRegion
          if (q->listPositions[j] > 0) { // add aNode in a new aRegion
	          a = getANodeNextPosition(q, j); // this is not deallocated!!!!!!!!!!!!!!
            
            // Adding a new aregion r2 after tail 
            // Alternative: deleting old end of tail, and adding a new element if it is closer than the old end (did not work well, so I left an easier solution)
            r2 = getARegion(a, arguments->M); // sets active to 1
            if (prev) {
              prev = tail;
            }
            tail->next = r2; /* connect to existing aregion-s in the list */
            tail = r2;
            if (!activeFirst) { // ??
              activeFirst = tail;
              prevActiveFirst = prev;
            }
          }
	      } // end for
      }
      else { // list was empty, so new head and tail should be added
        initializeHeadTail(q, ahead, &tail, &activeFirst, &prev, &prevActiveFirst, arguments);
      }
    }
    q = q->next;
  } // end while over lNode-s
  
  // close all previously unclosed aRegions  
  prev = prevActiveFirst;
  rCurr = activeFirst; 
  while (rCurr) {
    if (rCurr->active) {
      computeARegion(rCurr, arguments, seq, mat, firstAligned, secondAligned, maxCols, iQuery, iSubject, leftBorders, 
				first, second, strandBorders, numOfQueries, tail, allSeqLength, fastScore
				/*, lambda, scoreMatrix, alphabetSize, ambCharPenalty*/);
      // computeARegion sets active to 0
			if (checkARegion(rCurr, arguments) == 0) { // delete rCurr if shorter than f, has ratio less than R, etc.
        r2 = rCurr->next; 
        deleteCurrARegion(ahead, &tail, &prev, rCurr); // connects prev and rCurr->next
        rCurr = r2;
        continue;
      }
    }
    prev = rCurr; // skip inactive
    if (rCurr) { // rCurr could have been set to NULL in above deleteCurrARegion
	    rCurr = rCurr->next;
    }    
  } // end while
}

/* delete current region from aRegion list, and redirect pointers */
void deleteCurrARegion(aRegion **ahead, aRegion **tail, aRegion **prev, aRegion *r/*, aRegion **activeLast, aRegion **activeFirst*/) {
  
  if (*tail == *ahead) { // delete the only element in the list
    *tail = *ahead = NULL;
    *prev = NULL; // ?? redundant
  }        
  else if (*ahead == r) { // delete r from the beginning of the list 
    *ahead = r->next;    
    *prev = NULL; // ?? redundant
  }
  else { // delete from the middle of the list
    if (*prev) {
      (*prev)->next = r->next;
    }
    if (*tail == r) {
      *tail = *prev;
    }
  }
  freeARegion(r);
  r = NULL;
}

aNode *getANodeNextPosition(lNode *p, int i) {
  aNode *a;

  // new node in list
  a = (aNode *)/*e*/malloc(sizeof(aNode));
  /* find starting positions (lb) of matching segments in Q and S */
  a->lbQ = p->lb;
  a->lbS = p->listPositions[i]; // take i-th element in subject list; this could be changed later
  a->sl = p->sl;
  a->next = NULL;
  return a;
}

// add aNode to the end of aRegion
void addANode(aRegion *rCurr, aNode *a, Args *arguments) {
  
  if (rCurr->start == NULL) {
    rCurr->start = a;
    rCurr->startPosQ = rCurr->start->lbQ;
    rCurr->startPosS = rCurr->start->lbS;
  }
  else {
    rCurr->end->next = a;
  }
  rCurr->end = a;
  rCurr->score += (a->sl - 1) * arguments->M;	// alignment part is added later in computeARegion 
  rCurr->endPosQ = rCurr->end->lbQ + rCurr->end->sl - 2;
  rCurr->endPosS = rCurr->end->lbS + rCurr->end->sl - 2;
}

// delete aNode from the end of aRegion
void deleteANode(aRegion *rCurr, Args *arguments) { 
  aNode *oldEnd = NULL, *prev = NULL;  
  
  oldEnd = rCurr->end;
  prev = rCurr->start;
  if (prev != rCurr->end) {
    while (prev->next != oldEnd) {
      prev = prev->next;
    }
    rCurr->end = prev;
    rCurr->end->next = NULL;
  }
  else {
    rCurr->end = rCurr->start = NULL;
  }
  free(oldEnd);
  rCurr->score -= (oldEnd->sl - 1) * arguments->M;	// alignment part is added later in computeARegion 
  if (rCurr->end) {
    rCurr->endPosQ = rCurr->end->lbQ + rCurr->end->sl - 2;
    rCurr->endPosS = rCurr->end->lbS + rCurr->end->sl - 2;
  }
}

/* compute start and end coordinate on forward strand */
void getSubjectStartEndFR(Int64 *subjectStart, Int64 *subjectEnd, Int64 *strandLen, Int64 subjectId, Int64 numOfQueries, aRegion *p, char strainFR) {

  if (strainFR == '-') { // reverse strand --> compute corresponding positions on fwd strand
	  *subjectStart = 2 * strandLen[numOfQueries + subjectId] - p->startPosS;
	  *subjectEnd = 2 * strandLen[numOfQueries + subjectId] - p->endPosS;
  }
  else {
	  *subjectStart = p->startPosS;
	  *subjectEnd = p->endPosS;	  
  }
}

/* set query specific lambda and scoreMatrix */
void setQueryLambdaMatrix(double *lambdaAll, float **scoreMatrixAll, int alphabetSize, float ambCharPenaltyAll, 
	Int64 i, double *arrayK) {
	int k, l;
	lambda = lambdaAll[i];
	blastK = arrayK[i];

	for (k = 0; k < alphabetSize; k++) {
		for (l = 0; l < alphabetSize; l++) {
			scoreMatrix[k][l] = scoreMatrixAll[i][k * alphabetSize + l];
		}
	}
	ambCharPenalty = ambCharPenaltyAll;
}

/* calculate gap penalty */
float getGapPenalty () {
	int i, j;
	float gapPenalty = 0;
	for (i = 0; i < ALPHABET_DNA_SIZE; i++) {
		for (j = i + 1; j < ALPHABET_DNA_SIZE; j++) {
			if (scoreMatrix[i][j] < gapPenalty) {
				gapPenalty = scoreMatrix[i][j];
			}
		}
	}
	return gapPenalty;
}

/* calculate score between two exact matches using scoreMatrix*/
double scoreMatch(Int64 start, Int64 end, Int64 iQuery, char *seq, Int64 *leftBorders) {
  Int64 j;
  char nucl;
  int pos;
  double score = 0.0;

  for (j = start; j <= end; j++) {
    nucl = seq[leftBorders[iQuery] + j];
    pos = strchr(strNucl, nucl) - strNucl;
    score += scoreMatrix[pos][pos];
  }
  return score;
}

/* calculate score between two aligned regions using scoreMatrix*/
double scoreAlignedRegion(Int64 len, char *qSeq, char *sSeq, double gapPenaltyLocal) {
  Int64 j;
  int posQ, posS;
  double score = 0.0;
  char qNucl, sNucl;
  for (j = 0; j < len; j++) {
    qNucl = qSeq[j];
    sNucl = sSeq[j];
    if (qNucl == toupper(AMBIGUOUS_CHAR_QUERY) || qNucl == toupper(AMBIGUOUS_CHAR_SUBJECT)
      || sNucl == toupper(AMBIGUOUS_CHAR_SUBJECT) || sNucl == toupper(AMBIGUOUS_CHAR_QUERY)
      || qNucl == '-' || sNucl == '-') { // treat ambiguous nucleotides and gaps as gaps
      score += gapPenaltyLocal;
      //if (qNucl == '-') { // do not count gap on query side in length
      //  --j;
      //}
    }
    else {
      posQ = strchr(strNucl, qNucl) - strNucl;
      posS = strchr(strNucl, sNucl) - strNucl;
      score += scoreMatrix[posQ][posS];
    }
  }
  return score;
}
