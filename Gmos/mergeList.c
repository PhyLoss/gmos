/***** mergeList.c *************************************************************
 * Description: Functions for merging elements of linked lists.
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "commonSC.h"
#include "eprintf.h"
#include "stringUtil.h"
#include "expectedShulen.h"
#include "intervalStack.h"
#include "interval.h"
#include <assert.h>

#include "interface.h"
#include "nw.h"
#include "queryBTNode.h"
#include "lNode.h"
#include "align.h"
#include "mergeList.h"

#if defined(_DEBUG) && defined(WIN) 
#include "leakWatcher.h"
#endif

#if defined(_DEBUG) && defined(WIN) 
  #define new DEBUG_NEW
  #undef THIS_FILE
  static char THIS_FILE[] = __FILE__;
#endif

void freeMNode(mNode *m) {
  freeARegion(m->elem);
  free(m);
}
  
mNode *getMNode(aRegion *p, Int64 subjectId) {
  mNode *m;
  // new node in list
  m = (mNode *)/*e*/malloc(sizeof(mNode));
  m->elem = p;
  m->subjectId = subjectId;
  m->next = NULL;
  return m;
}

/* construct |Q| lists of pointers to first elements (heads) of each list; list must be sorted in asc order of lb of first elem-s */
void getCandidateList(aRegion **head, Int64 numOfQueries, Int64 numOfSubjects, mNode ***candidateList, short argB) {

  Int64 i, j;
  //lNode *p;
  aRegion *p = NULL;
  mNode *first = NULL, *m = NULL, *prev = NULL, *n = NULL;

	for (i = 0; i < numOfQueries; i++) {
    first = NULL;
    p = head[i * numOfSubjects + 0];
    if (p) {
      first = getMNode(p, 0);
    }
    for (j = 1; j < numOfSubjects; j++) {
      // insert every new node in asc order of lb
      p = head[i * numOfSubjects + j];
      
      // find a position in a merged list for a new candidate node
      m = first;
      if (p) {
        //while (m && p->start->lbQ > m->elem->start->lbQ) {
        while (m && p->startPosQ > m->elem->startPosQ) {
          prev = m;
          m = m->next;
        } // end while
        n = getMNode(p, j);
        if (m == first) { // adding a new node at the beginning of a merged list        
          n->next = first;
          first = n;
        }
        else { // adding a new node in the middle or at the end of a merged list
          prev->next = n;
          n->next = m;
        }
      } // end if
    } // end for
    
    (*candidateList)[i] = first;
  } // end for 
}

/* adding a new element to a candidateListI (candidateList od Qi): 
 * when an element of subject Sj is removed from candidateList, the next element of Sj (p) is added 
 * to candidateList at the appropriate position (in asc order of lb) */
void addCandidateElem(mNode **candidateListI, mNode *p) {
//void addCandidateElem(mNode ***candidateListI, mNode *p) {

  mNode *m = NULL, *prev = NULL;;

  // find a position in a merged list for a new candidate node
  m = *candidateListI;
  //m = **candidateListI;
  //while (m && p->elem->start->lbQ > m->elem->start->lbQ) {
  while (m && p->elem->startPosQ > m->elem->startPosQ) {
    prev = m;
    m = m->next;
  }
  if (m == *candidateListI) { // adding a new node at the beginning of candidateList        
  //if (m == **candidateListI) { // adding a new node at the beginning of candidateList        
    p->next = *candidateListI;
    *candidateListI = p;
    //p->next = **candidateListI;
    //**candidateListI = p;
  }
  else { // adding a new node in the middle or at the end of candidateList
    prev->next = p;
    p->next = m;
  }
}


/* checkDeleteLast4 returns 1 (keep new, delete last), -1 (keep last and delete new) or 0 (keep both last and new) */
int checkDeleteLast4(mNode *n, mNode *last, long long minTotalLen) {
  
  int retValue = 0;
  Int64 nStart, nEnd, lastStart, lastEnd;
  long long nLen, lastLen;
  double nScore, lastScore;

  return 0;
  nStart = n->elem->startPosQ;
  nEnd = n->elem->endPosQ;
  lastStart = last->elem->startPosQ;
  lastEnd = last->elem->endPosQ;  
  nLen = n->elem->totalLen;
  lastLen = last->elem->totalLen;
  
  /* alternatively, compare e-values */
	nScore = n->elem->ratio;
	lastScore = last->elem->ratio;

  // (1) if scores are equal, then keep both if they are of the same length, or keep longer if it contains the other one
  if (lastScore == nScore) {    
    if (nLen == lastLen) {
      retValue = 0;
    }
    else if (nEnd <= lastEnd) {
      retValue = -1;
    }
    else if (nLen > lastLen && nStart == lastStart) {
      retValue = 1;
    }
  }
  else if (lastScore >= nScore && nEnd <= lastEnd) { // n is shorter and has lower score (n contained within last)
    retValue = -1;
  }
  else if (lastScore < nScore && lastStart == nStart && nEnd >= lastEnd) { // last is shorter and has lower score(last contained within n)
    retValue = 1;
  }
  else { // keep both in all other options
    retValue = 0;
  }
  //retValue = 0;
  return retValue;
}


/* construct finalList: starting from the initial candidateLIst and subsequently adding other elements from linked lists */
void constructFinalList(mNode ***candidateList, aRegion **head, Int64 numOfQueries, Int64 numOfSubjects, mNode ***finalList, short argB, 
      char **seqHeaders, Int64 *strandBorders, Int64 *leftBorders, FILE *fqout, long long allSeqLength, 
			Int64 *strandLen, double *lambdaAll) {
	
  Int64 i, j, subjectId;
  mNode *last = NULL, *n = NULL;
  aRegion *p = NULL;
  char delimiter = '\t';
  long long minTotalLen = 0;
  int retCheckDeleteLast = 0;
  char strainFR = '+'; // forward/reverse strand (+/-)
  Int64 subjectStart = 0, subjectEnd = 0;

  if (fqout) {
    fprintf(fqout, "------------- Candidate List -------------- \n");
    fflush(fqout);
  }

  for (i = 0; i < numOfQueries; i++) {
    (*finalList)[i] = NULL;
    // print iQuery and all subject lengths in file, so it can be later used for display
    if (fqout) {
      fprintf(fqout, "\n------------- Query %lld -------------- \n", (long long)i + 1);
      fflush(fqout);
      fprintf(fqout, "%lld\n", (long long int) strandLen[i]); // query length - fwd strand
      fflush(fqout);
      for (j = numOfQueries; j < numOfQueries + numOfSubjects; j++) {
        fprintf(fqout, "%lld\n", (long long int) strandLen[j]); // query length - fwd strand    
        fflush(fqout);
      }
    }
		lambda = lambdaAll[i];
    while ((*candidateList)[i]) { // while there are elements in candidateList, add them to finalList; at the beginning: only first elem-s are in candidateList
      n = (*candidateList)[i]; // n will be moved from candidateList to finalList
      // NEW - print all candidates - query fwd strand
      //if (leftBorders[i] + n->elem->endPosQ <= strandBorders[i]) {
      if (n->elem->endPosQ <= strandLen[i]) {
        // add ratio = score / fragment_length
        computeEValueScoreBit(n->elem, allSeqLength/*, lambda[i]*/);
        if (fqout) {
			    fprintf(fqout, "Q%lld(+)%12lld%12lld%12lld%9.1lf%6.3lf%c%10.6e%c", (long long)(i + 1), (long long)n->elem->startPosQ + 1, (long long)n->elem->endPosQ + 1,
				        (long long)n->elem->endPosQ - n->elem->startPosQ + 1, n->elem->score, n->elem->ratio, delimiter,
				          n->elem->e_value, delimiter); // 24/01/2015
          fflush(fqout);          
          strainFR = (n->elem->startPosS >= strandLen[numOfQueries + n->subjectId]) ? '-' : '+';
          getSubjectStartEndFR(&subjectStart, &subjectEnd, strandLen, n->subjectId, numOfQueries, n->elem, strainFR);
          fprintf(fqout, "S%-6lld%-65.65s%12lld%12lld%c%c\n", (long long)n->subjectId + 1, seqHeaders[numOfQueries + n->subjectId] + 1, 
                  (long long)subjectStart + 1, (long long)subjectEnd + 1, delimiter, strainFR); // last position        
          fflush(fqout);
        }
      }
      // END NEW

      (*candidateList)[i] = (*candidateList)[i]->next; // move to the next elem in candidateList
      n->next = NULL;
      subjectId = n->subjectId;

      //head[i * numOfSubjects + subjectId] = head[i * numOfSubjects + subjectId]->next; // ?????? if == null
      if (head[i * numOfSubjects + subjectId]) { // this has to be done be here because of possible deletion of n in freeMNode(n) (*)
        head[i * numOfSubjects + subjectId] = head[i * numOfSubjects + subjectId]->next;
      }

      if ((*finalList)[i] == NULL) { // add the first element in finalList
        (*finalList)[i] = n;        
        last = n;
      }
      else { // add at the end of finalList
        // if n overlapps with one or more elements in finalList, then leave element(s) with the longest sl, and others should be removed from the list.
        // Note: Only the last added elem in finalList (prev) could be overlapping with n; previous to prev should have been already taken care of when prev was added ???
        // NEW: overlapping allowed for at most -O option percentage between two regions to be both displayed in final list
        minTotalLen = (last->elem->totalLen < n->elem->totalLen) ? last->elem->totalLen : n->elem->totalLen; // totalLen is >= (elem->endPosQ - elem->startPosQ + 1)
    
        if (n->elem->startPosQ > last->elem->endPosQ) { // new: 1) if new and last do not overlap, keep them both
          last->next = n;
          last = n;
        }
        else { // overlapping, i.e. n->elem->lb <= last->elem->rb
          //if (checkDeleteLast(n, last) == 1) { 
          retCheckDeleteLast = checkDeleteLast4(n, last, minTotalLen); // delete shorter segment within longer and with lower score
          // new: checkDeleteLast returns 1 (keep new, delete last), -1 (keep last and delete new) or 0 (keep both last and new)
          if (retCheckDeleteLast == 1) { // delete last, leave new
            freeARegion(last->elem);
            last->elem = n->elem;
            last->subjectId = n->subjectId;
            free(n); // not freeMNode!
            n = NULL;
          }
          else if (retCheckDeleteLast == -1) { // leave last, delete n            
            freeMNode(n); // (*)
            n = NULL;
          }
          else { // retCheckDeleteLast == 0 --> leave both
            last->next = n;
            last = n;          
          }
        } // end else
      } // end else
      
      // get next elem from linked list of subjectId and add it to candidateList
      if (head[i * numOfSubjects + subjectId]) { // if there are still elements in subjectId linked list
        p = head[i * numOfSubjects + subjectId]; // move to the next elem in subjectId linked list (p)
        n = getMNode(p, subjectId);    
        addCandidateElem(&(*candidateList)[i], n); // add p in asc order of lb in candidate list ?? candidateList by ref??
      } // end if

    } // end while

  }
}

// print a single row
//void printRow(FILE *fpout, aRegion *p, Int64 queryStart, Int64 queryEnd, char queryStrand, char **seqHeaders, Int64 *strandBorders, Int64 *leftBorders, Int64 subjectId, Int64 numOfQueries) {
void printRow(FILE *fpout, aRegion *p, Int64 queryStart, Int64 queryEnd, char queryStrand, char **seqHeaders, Int64 *strandLen, Int64 subjectId, Int64 numOfQueries) {

  char strainFR = '+'; // forward straind: '+'; reverse strand: '-'
  char delimiter = '\t';
  Int64 subjectStart = 0, subjectEnd = 0;	
  
  //fprintf(fpout, "%12lld%12lld\t%c%12lld%12.1lf%c%4.2lf%c", (long long)queryStart + 1, (long long)queryEnd + 1, queryStrand, (long long)queryEnd - queryStart + 1, 
  fprintf(fpout, "%12lld%12lld\t%c%12lld%12.1lf%c", (long long)queryStart + 1, (long long)queryEnd + 1, queryStrand, (long long)queryEnd - queryStart + 1, 
	        p->score, delimiter/*, p->ratio, delimiter*/);			
	fprintf(fpout, "%10.6e\t", (double)p->e_value); // new - print E-value
  fflush(fpout);
  // check forward/reverse strand
  strainFR = (p->startPosS >= strandLen[numOfQueries + subjectId]) ? '-' : '+';
  getSubjectStartEndFR(&subjectStart, &subjectEnd, strandLen, subjectId, numOfQueries, p, strainFR);
  fprintf(fpout, "S%-5lld%c%80.80s%c%12lld%12lld%c%c\n", (long long)subjectId + 1, delimiter, seqHeaders[numOfQueries + subjectId] + 1, delimiter,
				  (long long)subjectStart + 1, (long long)subjectEnd + 1, delimiter, strainFR); // last position
  fflush(fpout);
}

//void printCandidateList(mNode **mList, FILE *fpout, Int64 numOfQueries, char *header, char **seqHeaders, Int64 *strandBorders, Int64 *leftBorders, short argB) {
void printCandidateList(mNode **mList, FILE *fpout, Int64 numOfQueries, char *header, char **seqHeaders, short argB, Int64 *strandLen) {
  
  Int64 i;
  mNode *m;
  aRegion *p;
  
  for (i = 0; i < numOfQueries; i++) {    
    m = mList[i];
    fprintf(fpout, "\n------------- %s: (%lld) %s -------------\n", header, (long long)i + 1, seqHeaders[i] + 1);
    fflush(fpout);
    //while (m && m->elem->endPosQ + leftBorders[i] <= strandBorders[i]) {
    while (m && m->elem->endPosQ <= strandLen[i]) {
      p = m->elem;
	    printRow(fpout, p, p->startPosQ, p->endPosQ, '+', seqHeaders, strandLen, m->subjectId, numOfQueries);
      m = m->next;
    } // end while
  }
  fflush(fpout);
}

/* print a single row in alignment */
void printAlignmentRow(aRegion *r, int *k, int *lenTemp, char **tempSeq, Int64 iSeq, Int64 a_sl, Int64 lb, Int64 *leftBorders, char *seq, FILE *fout) {

  if (a_sl - 1 > *lenTemp) { // resize tempSeq
    *lenTemp = a_sl - 1;
    *tempSeq = (char *)realloc(*tempSeq, sizeof(char) * (*lenTemp + 1));
  }
  strncpy(*tempSeq, seq + leftBorders[iSeq] + lb, (size_t)a_sl - 1); // query match
  (*tempSeq)[a_sl - 1] = '\0';
  fprintf(fout, "%s", *tempSeq);
  fflush(fout);
  if (*k < r->numAlignedSegment * 2) {
    fprintf(fout, "%s", r->alignedSegment[*k]); // query aligned fragment			
    fflush(fout);
    *k += 2;
  }
}

/* print alignment in multi fasta format, e.g.
 * >Q[1..1000]
 * ACGGGT...   
 * >S1[2057..3057]
 * CCGGGT...   
*/
void printAlignment(mNode **mList, FILE *fout, Int64 numOfQueries, char **seqHeaders, Int64 *strandLen, char *seq, Int64 *leftBorders) {
  
  Int64 i;
  mNode *m;
  aRegion *p;
  char queryName[MAXSEQNAME + 1];
  char subjectName[MAXSEQNAME + 1];
  char *tempSeq = NULL;
  int lenTemp = 1000; // initial value
  char strainFR = '+';
  aNode *a;
  int k = 0;
  Int64 subjectStart = 0, subjectEnd = 0;

  tempSeq = (char *)malloc(sizeof(char) * (lenTemp + 1));
  tempSeq[0] = '\0';
  for (i = 0; i < numOfQueries; i++) {    
    m = mList[i];
    while (m && m->elem->endPosQ <= strandLen[i]) { // print fwd strand
      strcpy(queryName, seqHeaders[i] + 1);
      strcpy(subjectName, seqHeaders[numOfQueries + m->subjectId] + 1);
      p = m->elem;
      
      // (1a) print query name and aligned part before first match
      fprintf(fout, ">%s[%lld .. %lld]+\n", queryName, (long long)p->startPosQ + 1, (long long)p->endPosQ + 1);
      fflush(fout);
      if (p->alignedStart[0]) {
        fprintf(fout, "%s", p->alignedStart[0]);                        
        fflush(fout);
      }
      // (1b) print query sequence (exact matches + aligned fragments)
      k = 0;
      a = p->start;
      while (a) { 
        if (a->sl > 0) { // a->sl <= 0 possibly at the end or the beginning of aRegion (due to split while looking for breakpoints)
          printAlignmentRow(p, &k, &lenTemp, &tempSeq, i, a->sl, a->lbQ, leftBorders, seq, fout);
        }
        a = a->next;
      } // end while
      // (1c) print aligned part after last query match
      if (p->alignedEnd[0]) {
        fprintf(fout, "%s", p->alignedEnd[0]);			
      }
      fprintf(fout, "\n");			
      fflush(fout);
      
      // (2a) print subject name and aligned part before first match
      // check forward/reverse strand
      strainFR = (p->startPosS >= strandLen[numOfQueries + m->subjectId]) ? '-' : '+';
      getSubjectStartEndFR(&subjectStart, &subjectEnd, strandLen, m->subjectId, numOfQueries, p, strainFR);
      fprintf(fout, ">%s[%lld .. %lld]%c\n", subjectName, (long long)subjectStart + 1, (long long)subjectEnd + 1, strainFR);
      fflush(fout);
      if (p->alignedStart[1]) {
        fprintf(fout, "%s", p->alignedStart[1]);                        
        fflush(fout);
      }
      // (2b) print subject sequence (exact matches + aligned fragments)
      k = 1;
      a = p->start;
      while (a) {
        if (a->sl > 0) { // a->sl <= 0 possibly at the end or the beginning of aRegion (due to split while looking for breakpoints)
          printAlignmentRow(p, &k, &lenTemp, &tempSeq, m->subjectId + numOfQueries, a->sl, a->lbS, leftBorders, seq, fout);
        }
        a = a->next;
      } // end while
      // (2c) print aligned part after last subject match
      if (p->alignedEnd[1]) {
        fprintf(fout, "%s", p->alignedEnd[1]);			
      }
      fprintf(fout, "\n");			
      fflush(fout);
      m = m->next;
    } // end while
  }
  free(tempSeq);
}


/* free candidateList*/
void freeMList(mNode **mList, Int64 numOfQueries, Int64 numOfSubjects) {
  
  Int64 i;
  mNode *m, *n;

	for (i = 0; i < numOfQueries; i++) {
    m = mList[i];
    while (m) {
      n = m->next;
      freeMNode(m);
      m = n;
    } // end while
  }
  free(mList);
}

/* extend segment n towards m - its predecessor (m's end part, i.e. the segment after the last match, stays as it is) */
void extendNToM(mNode *n, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned, 
                Int64 numOfQueries, Int64 Qi, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength, 
								Int64 mEndPosQ, char fastScore) {
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty) {
  Int64 lenQ = 0, lenS = 0;
  Int64 Sj = n->subjectId + numOfQueries; // relative subject index corrected for numOfQueries

  lenQ = n->elem->start->lbQ - mEndPosQ;
  lenS = lenQ;
  if (n->elem->start->lbS - lenS < 0) {
    lenS = n->elem->start->lbS;
  }
	alignSegment(mEndPosQ + 1, n->elem->start->lbS - lenS + 1, leftBorders, lenQ, lenS, first, second, firstAligned,
		secondAligned, seq, Qi, Sj, mat, args, maxCols, fastScore);
		//, scoreMatrix, alphabetSize, ambCharPenalty);

  // correct score, totalLen and new right border 
  n->elem->score -= n->elem->scoreStart;
  if (n->elem->alignedStart[0] != NULL) {
    n->elem->totalLen -= strlen(n->elem->alignedStart[0]);  
  }
  storeAlignedSegments(&n->elem->scoreStart, mat[lenQ * maxCols + lenS].score, &n->elem->score, firstAligned, secondAligned, &n->elem->alignedStart[0], &n->elem->alignedStart[1]);
  
  n->elem->totalLen += strlen(firstAligned);
  n->elem->startPosQ = mEndPosQ + 1;
  n->elem->startPosS = n->elem->end->lbS - lenS + 1;
  n->elem->ratio = n->elem->score / n->elem->totalLen;
  computeEValueScoreBit(n->elem, allSeqLength/*, lambda*/);
}

/* extend segment m up to n (n's begin part, i.e.the part before the first match, stays as it is) */
void extendMToN(mNode *m, mNode *n, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned, 
                Int64 numOfQueries, Int64 Qi, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength, 
								Int64 nStartPosQ, char fastScore) { 
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty) {

  Int64 lenQ = 0, lenS = 0;
  Int64 Sj = m->subjectId + numOfQueries; // relative subject index corrected for numOfQueries

  lenQ = nStartPosQ - (m->elem->end->lbQ + m->elem->end->sl - 1);
  lenS = lenQ;
  if (leftBorders[Sj] + m->elem->end->lbS + m->elem->end->sl - 1 + lenS > strandBorders[Sj]) { // Sj -- subject of m 
    lenS = strandBorders[Sj] - (m->elem->end->lbS + m->elem->end->sl - 1 + leftBorders[Sj]);
  }
	alignSegment(m->elem->end->lbQ + m->elem->end->sl - 1, m->elem->end->lbS + m->elem->end->sl - 1, leftBorders, lenQ, lenS, first, second, firstAligned, secondAligned,
		seq, Qi, Sj, mat, args, maxCols, fastScore);
							//scoreMatrix, alphabetSize, ambCharPenalty);

  // correct score, totalLen and new right border 
  m->elem->score -= m->elem->scoreEnd;
  if (m->elem->alignedEnd[0] != NULL) {
    m->elem->totalLen -= strlen(m->elem->alignedEnd[0]); 
  }
  storeAlignedSegments(&m->elem->scoreEnd, mat[lenQ * maxCols + lenS].score, &m->elem->score, firstAligned, secondAligned, &m->elem->alignedEnd[0], &m->elem->alignedEnd[1]);

  m->elem->totalLen += strlen(firstAligned);
  m->elem->endPosQ = nStartPosQ - 1;
  m->elem->endPosS = m->elem->end->lbS + m->elem->end->sl + lenS - 2; //??? lenS or strlen(m->elem->alignedEnd[1])??
  m->elem->ratio = m->elem->score / m->elem->totalLen;
  computeEValueScoreBit(m->elem, allSeqLength/*, lambda*/);
}

/* delete alignment before the first anode (start) of n and adjust score and totalLen */
void deleteAlignmentBeforeStart(mNode *n) {

  if (n->elem->alignedStart[0] != NULL) {
    n->elem->totalLen -= strlen(n->elem->alignedStart[0]);
    n->elem->score -= n->elem->scoreStart;
  }
  n->elem->scoreStart = 0;
  free(n->elem->alignedStart[0]);
  free(n->elem->alignedStart[1]);
  n->elem->alignedStart[0] = n->elem->alignedStart[1] = NULL;
}

/* delete alignment after the last anode (end ) of n and adjust score and totalLen */
void deleteAlignmentAfterEnd(mNode *n) {

  if (n->elem->alignedEnd[0] != NULL) {
    n->elem->totalLen -= strlen(n->elem->alignedEnd[0]);
  }
  n->elem->score -= n->elem->scoreEnd;
  n->elem->scoreEnd = 0;
  free(n->elem->alignedEnd[0]);
  free(n->elem->alignedEnd[1]);
  n->elem->alignedEnd[0] = n->elem->alignedEnd[1] = NULL;
}

/* delete alignedSegments starting from the z-th element onwards */
void deleteAlignedSegments(mNode *n, int z) {
  
  int k;
  for (k = z; k < n->elem->numAlignedSegment; k++) {
    n->elem->score -= n->elem->scoreAlignedSegment[k];
    n->elem->totalLen -= strlen(n->elem->alignedSegment[2 * k]);
    free(n->elem->alignedSegment[2 * k]);
    free(n->elem->alignedSegment[2 * k + 1]);
    n->elem->alignedSegment[2 * k] = NULL;
    n->elem->alignedSegment[2 * k + 1] = NULL;
    -- n->elem->numAlignedSegment;
  }
}

/* delete aligned segments before the position z */
void deleteAlignedSegmentsBeforeZ(mNode *n, int z) {
  
  int k;  
  for (k = z; k < n->elem->numAlignedSegment; k++) {
      n->elem->alignedSegment[(k - z) * 2] = (char *)realloc(n->elem->alignedSegment[(k - z) * 2], strlen(n->elem->alignedSegment[k * 2]) + 1);
      n->elem->alignedSegment[(k - z) * 2 + 1] = (char *)realloc(n->elem->alignedSegment[(k - z) * 2 + 1], strlen(n->elem->alignedSegment[k * 2 + 1]) + 1);
      strcpy(n->elem->alignedSegment[(k - z) * 2], n->elem->alignedSegment[k * 2]);
      strcpy(n->elem->alignedSegment[(k - z) * 2 + 1], n->elem->alignedSegment[k * 2 + 1]);
      n->elem->scoreAlignedSegment[k - z] = n->elem->scoreAlignedSegment[k];
    }
    for (k = n->elem->numAlignedSegment - z; k < n->elem->numAlignedSegment; k++) {
      n->elem->scoreAlignedSegment[k] = 0;
      free(n->elem->alignedSegment[k * 2]);
      free(n->elem->alignedSegment[k * 2 + 1]);
      n->elem->alignedSegment[k * 2] = n->elem->alignedSegment[k * 2 + 1] = NULL;
    }
    n->elem->numAlignedSegment -= z;
}


/* extend segment m over n */
void extendMToNOverlap(mNode **mm, mNode **nn, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned, 
                Int64 numOfQueries, Int64 i, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength, 
								mNode **mprev, char fastScore) {
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty) {

  aNode *p = NULL;
  aNode *oldStart = NULL;
  int cnt = 0;
  mNode *m, *n;
  
  m = *mm;
  n = *nn;
  p = n->elem->start;
  cnt = -1;
  while (p != NULL) {
    if (p->lbQ > m->elem->endPosQ) {
      break;
    }
    else {
      // delete previous aNode and aligned part
      oldStart = n->elem->start;
      n->elem->totalLen -= oldStart->sl - 1;
      n->elem->score -= (oldStart->sl - 1) * args->M; // exact match
      n->elem->start = oldStart->next; // new start
      free(oldStart);
      oldStart = NULL;
      
      if (cnt == -1) { // delete alignment part before the first match (alignedStart)
        cnt = 0;
        deleteAlignmentBeforeStart(n);
      }
      // later delete aligned part(s) between matches; now just count how many of them need to be deleted;
      // also resolve scores and lengths of deleted alignments
      ++ cnt; // number of deleted aligned segments; deleted a.s. is always after p
      if (cnt <= n->elem->numAlignedSegment) {
        n->elem->score -= n->elem->scoreAlignedSegment[cnt - 1];
        n->elem->totalLen -= strlen(n->elem->alignedSegment[cnt * 2 - 2]);
      }
    } // end else
    p = n->elem->start;
  } // end while
          
  // m and n overlap only in alignment segment before anode start of n, so only delete that part
  if (cnt == -1) { 
    deleteAlignmentBeforeStart(n);
  }
  // delete cnt aligned parts between segments and close n from the left side;
  // move aligned segments from the positions cnt to numAlignedSegment to the left and update numAlignedSegment
  else {
    if (cnt >= 1) {
      deleteAlignedSegmentsBeforeZ(n, cnt - 1);
    }
  }
  // close the left border of n at p->lbQ
  if (p != NULL) {
    n->elem->start = p;
    n->elem->startPosQ = n->elem->start->lbQ;
    n->elem->startPosS = n->elem->start->lbS;
    n->elem->ratio = n->elem->score / n->elem->totalLen;
    computeEValueScoreBit(n->elem, allSeqLength/*, lambda*/);
  }
  // if n is now shorter <= args->f or has lower ratio then args->R, discard it
  if (p == NULL || (p != NULL && (n->elem->ratio < args->R || n->elem->totalLen < args->f))) {
    m->next = n->next;
    freeMNode(n);
    n = m->next;
  }
  else {
    // extend m to the new start of n (p->lbQ - 1); this could be a bit longer then before
		extendMToN(m, n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, i, args,
			maxCols, seq, mat, allSeqLength, n->elem->startPosQ, fastScore);
			//lambda, scoreMatrix, alphabetSize, ambCharPenalty);  
    *mprev = m;
    m = n;
    n = n->next;
  }
  *mm = m;
  *nn = n;
}


/* extend segment n over m */
void extendNToMOverlap(mNode **mm, mNode **nn, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned, 
                Int64 numOfQueries, Int64 i, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength, 
								mNode **mprev, mNode ***finalList, char fastScore) {
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty) {

  aNode *p = NULL, *pprev = NULL, *old;
  int cnt = 0;
  Int64 mPosEnd = 0, newNPosStartQ = 0;
  mNode *m = NULL, *n = NULL;
  m = *mm;
  n = *nn;

  p = m->elem->start;
  while (p != NULL) {
    if (p->lbQ + p->sl - 2 >= n->elem->startPosQ) {
      mPosEnd = (pprev != NULL) ? (pprev->lbQ + pprev->sl - 2) : 0;
      m->elem->endPosQ = mPosEnd;
      m->elem->endPosS = (pprev != NULL) ? (pprev->lbS + pprev->sl - 2) : 0;
      m->elem->end = pprev;
      if (pprev) {
        pprev->next = NULL;
      }
      break;
    }
    else {
      ++ cnt; // counting anodes that are before overlap (if cnt is 0, then the overlap is before the end of the start node)
      pprev = p;
      p = p->next;
    }
  } // end while

  // p can be the last anode in m, i.e. the end node, or can be any node between start and end (including start)
  
  // (1) if p is start, then delete m and connect n to m's predecessor (n is already covering overlapping region)
  if (p == m->elem->start) {
    if (*mprev) {
      (*mprev)->next = n;
    }
    else {
      *finalList[i] = n;
    }
    freeMNode(m);
    m = NULL;
  }
  // (2) if p is between start and end (including end), then delete all anodes from p to end (including p and end) and aligned regions between anodes
  else { 
    while (p != NULL) {
      old = p;
      m->elem->totalLen -= p->sl - 1;
      m->elem->score -= (p->sl - 1) * args->M; // exact match
      p = p->next; // new start
      free(old);
      old = NULL;
    } // end while

    if (cnt >= 1) { // if cnt > 1, then delete all alignedSegments after p
      deleteAlignedSegments(m, cnt - 1);
    }
    deleteAlignmentAfterEnd(m);
    
    // define new end of m
    m->elem->ratio = m->elem->score / m->elem->totalLen;
    computeEValueScoreBit(m->elem, allSeqLength/*, lambda*/);
  } // end else

  // if m is now shorter <= args->f or has lower ratio then args->R, discard it
  newNPosStartQ = n->elem->startPosQ - 1;
  if (m != NULL) {
    if (m->elem->ratio < args->R || m->elem->totalLen < args->f) {
      if (*mprev) {
        (*mprev)->next = n;
      }
      else {
        *finalList[i] = n;
      }
      freeMNode(m);
      m = NULL;
    }
    else {
      newNPosStartQ = m->elem->endPosQ;
    }
  }

  // extend n to the new end of m (pprev->lbS + pprev->sl - 2); this could be a bit longer then before
	extendNToM(n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, i, args, maxCols,
		seq, mat, allSeqLength, newNPosStartQ, fastScore);
		//lambda, scoreMatrix, alphabetSize, ambCharPenalty);
  *mm = n;
  *nn = n->next;

}


/* split m and n somewhere around middle point (currently: n is "preferred" a bit more, i.e. aNode containing m.P. stays in n, and m is shortened up to the m.P.) */
void splitMNMiddle(mNode **mm, mNode **nn, Int64 *strandBorders, Int64 *leftBorders, char *first, char *second, char *firstAligned, char *secondAligned, 
                Int64 numOfQueries, Args *args, int maxCols, char *seq, element *mat, long long allSeqLength, 
								mNode **mprev, Int64 iQuery, mNode ***finalList, char fastScore) {
								//double lambda, float *scoreMatrix, int alphabetSize, float ambCharPenalty) {

  // Start from m->elem->start and n->elem->start. 
  // The new end of m, denoted p, is the last aNode that is completely before the middle point .
  // The new start element of n, denoted r, is an aNode that includes 
  // a) the middle point or (if the middle point goes through r)
  // b) the first region just after the middle point (if the middle point goes between two aNode-s)
  // If p should not be overlapping with r
  // In all cases, recompute them new start/end elements of m/n, their scores, totalLen-s, etc. 
  // If m and/or is shorter than args->f after, for now on --> leave them both????
  aNode *p = NULL, *r = NULL, *prev = NULL, *pnext = NULL;
  Int64 middlePoint = 0;
  Int64 pLen = 0, rLen = 0;
  int k = 1, z = 0;
  //Int64 newNPosStartQ = 0;
  mNode *m, *n;
  
  m = *mm;
  n = *nn;
  middlePoint = n->elem->startPosQ + (m->elem->endPosQ - n->elem->startPosQ) / 2;

  // find p
  m->elem->score = m->elem->scoreStart; // compute new score of m
  p = m->elem->start->next;
  while (p && p->lbQ + p->sl - 2 < middlePoint) {
    // add alignedSegment just before p
    m->elem->score += args->M * (p->sl - 1) + m->elem->scoreAlignedSegment[k - 1]; 
    prev = p;
    p = p->next;
    ++ k;
  }
  deleteAlignedSegments(m, k - 1);

  // delete all nodes from p to including m->elem->end
  while (p) {
    pnext = p->next;
    free(p);
    p = pnext;
  }
  p = prev;
  m->elem->end = p;
  p->next = NULL;

  // find r
  prev = NULL;
  //n->elem->score -= n->elem->scoreStart;
  // delete aligned part just before n->elem->start
  deleteAlignmentBeforeStart(n);
  r = n->elem->start;
  while (r && r->lbQ + r->sl - 2 < middlePoint) {
    n->elem->score -= args->M * (r->sl - 1);
    //if (prev != NULL) { // in deleteAlignedSegmentsBeforeZ
    //  n->elem->score -= n->elem->scoreAlignedSegment[z];
    //}
    ++ z;
    prev = r;
    r = r->next;
    free(prev);
    prev = NULL;
  }
  if (z > 0) {    
    deleteAlignedSegmentsBeforeZ(n, z);
  }
  // (1) if r-> lb is after the middle point, then shorten m to the middle point (add aligned region from the end of p to the middle point-1) and let n begin from the middle point 
  // (extend start of n from the middle pint to the start of r)
  n->elem->start = r;
  // align segment after p up to the middle point --> this will be new end of m (the same for both cases (1) and (2))
  
  if (r->lbQ > middlePoint) {
    pLen = middlePoint - (p->lbQ + p->sl - 1);
    m->elem->endPosQ = middlePoint - 1;
    m->elem->endPosS = p->lbS + p->sl - 2 + pLen;
		alignSegment(p->lbQ + p->sl - 1, p->lbS + p->sl - 1, leftBorders, pLen, pLen, first, second,
			firstAligned, secondAligned, seq, iQuery, numOfQueries + m->subjectId, mat, args, maxCols, fastScore);
			//scoreMatrix, alphabetSize, ambCharPenalty);
    storeAlignedSegments(&m->elem->scoreEnd, mat[pLen * maxCols + pLen].score, &m->elem->score, firstAligned, secondAligned, &m->elem->alignedEnd[0], &m->elem->alignedEnd[1]);
    // align segment just before r
    rLen = r->lbQ - middlePoint;
		alignSegment(middlePoint, r->lbS - rLen + 1, leftBorders, rLen, rLen, first, second, firstAligned, secondAligned, seq, iQuery,
			numOfQueries + n->subjectId, mat, args, maxCols, fastScore);
			//scoreMatrix, alphabetSize, ambCharPenalty);
    storeAlignedSegments(&n->elem->scoreStart, mat[rLen * maxCols + rLen].score, &n->elem->score, firstAligned, secondAligned, &n->elem->alignedStart[0], &n->elem->alignedStart[1]);
    n->elem->startPosQ = middlePoint;
    n->elem->startPosS = r->lbS - rLen + 1;
  }
  // 2) if the middle point is within r, then shorten m up to the left border of r (p may also be shortened if it overlaps with r? --> NO! this could not happen
  // since p has already been computed to be BEFORE the middle point)
  else {
    m->elem->endPosQ = r->lbQ - 1;
    pLen = r->lbQ - (p->lbQ + p->sl - 1);
		alignSegment(p->lbQ + p->sl - 1, p->lbS + p->sl - 1, leftBorders, pLen, pLen, first, second, firstAligned,
			secondAligned, seq, iQuery, numOfQueries + m->subjectId, mat, args, maxCols, fastScore);
			//scoreMatrix, alphabetSize, ambCharPenalty);
    storeAlignedSegments(&m->elem->scoreEnd, mat[pLen * maxCols + pLen].score, &m->elem->score, firstAligned, secondAligned, &m->elem->alignedEnd[0], &m->elem->alignedEnd[1]);
    m->elem->endPosS = p->lbS + p->sl - 2 + pLen;      
    n->elem->startPosQ = r->lbQ;
    n->elem->startPosS = r->lbS - pLen + 1;
    storeAlignedSegments(&n->elem->scoreStart, 0, &n->elem->score, "", "", &n->elem->alignedStart[0], &n->elem->alignedStart[1]); // there is no aligned start, since n starts witm r          }
  }
  
  m->elem->totalLen = m->elem->endPosQ - m->elem->startPosQ + 1;
  n->elem->totalLen = n->elem->endPosQ - n->elem->startPosQ + 1;
  m->elem->ratio = m->elem->score / m->elem->totalLen;
  n->elem->ratio = n->elem->score / n->elem->totalLen;
  computeEValueScoreBit(n->elem, allSeqLength/*, lambda*/);
  computeEValueScoreBit(m->elem, allSeqLength/*, lambda*/);

  // if m is now shorter <= args->f or has lower ratio then args->R, discard it
  //newNPosStartQ = n->elem->startPosQ - 1;
  if (m->elem->ratio < args->R || m->elem->totalLen < args->f) {
    if (*mprev) {
      (*mprev)->next = n;
    }
    else {
      *finalList[iQuery] = n;
    }
    freeMNode(m);
  }

  // if n is now shorter <= args->f or has lower ratio then args->R, discard it, but only if m is not null (so that at least one of intervals stays)
  *mm = n;
  *nn = n->next;
  if (m && (n->elem->ratio < args->R || n->elem->totalLen < args->f)) {
    m->next = n->next;
    freeMNode(n);
    *mm = m;
  }
  
}

/* function returns 1 if m is "better" than n, -1 if if n is "better" than m, and 0 if bot are equal when ratios and e-scores are compared */
int getBetterOverallScore(mNode *m, mNode *n) {
  int retValue = 0;

  if (m->elem->e_value < n->elem->e_value) { 
    retValue = 1;
  }
  else if (m->elem->e_value > n->elem->e_value) {
    retValue = -1;
  }
  else if (m->elem->ratio > n->elem->ratio) { // e-scores are equal
    retValue = 1;
  }
  else if (m->elem->ratio < n->elem->ratio) { // e-scores are equal
    retValue = -1;
  }
  return retValue;
}

/* connect mnode-s of final lists so there is no gap or overlap between them --> used when only the best answer (recombinant) needs to be dispayed */
void connectFinalList(mNode ***finalList, Int64 numOfQueries, Args *args, char *seq, Int64 *leftBorders, Int64 *strandBorders,
                      char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat, 
											long long allSeqLength, double *lambdaAll, float **scoreMatrixAll, int alphabetSize, 
											float ambCharPenaltyAll, char fastScore, double *arrayK) {

  mNode *m = NULL, *n = NULL, *prev = NULL;
  Int64 i = 0;
  Int64 lenApartQ = 0; // how far apart are two neighboring intervals on query part 
  Int64 lenApartS = 0; // how far apart are two neighboring intervals on subject part
  Int64 middlePoint = 0;
  Int64 maxLenApart = args->C;
  int betterOverallScore = 0;

  for (i = 0; i < numOfQueries; i++) {    
		setQueryLambdaMatrix(lambdaAll, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, i, arrayK);

    m = (*finalList)[i];
    prev = NULL;
    n = m->next;
    while (n != NULL && n->elem->endPosQ + leftBorders[i] <= strandBorders[i]) {
      // (1) if there is a gap between m and n, connect them if they are both query and subject at most args->C apart; otherwise let them stay apart
      // a) if m and are of the same subject type and both Q and S part at most args->C apart, connect them into a single interval
      // b) if they are of different subject types or of the same type but not in close asc order, 
      // connect them so the "better" interval (higher e-value or ratio) is extended towards the "worse" one
      // c) if e-values are equal, split in the middle

      lenApartQ = n->elem->startPosQ - m->elem->endPosQ - 1;
      lenApartS = n->elem->startPosS - m->elem->endPosS - 1;
      //minTotalLen = (m->elem->totalLen < n->elem->totalLen) ? m->elem->totalLen : n->elem->totalLen;
      betterOverallScore = getBetterOverallScore(m, n);

      if (lenApartQ <= maxLenApart && lenApartQ > 0) {
        // a) segments are in asc order and close enough both on query and subject part
        if (n->subjectId == m->subjectId && lenApartS > 0 && lenApartS <= maxLenApart) { 
          // delete begin of n and end of m, connect m->end and n->start, re-compute alignment between them
          // --> connect n and m in m and delete n->elem (in connectARegion) and n (here)
					connectARegion(m->elem, n->elem, args, seq, leftBorders, n->subjectId, i, numOfQueries, firstAligned,
						secondAligned, maxCols, first, second, mat, allSeqLength, fastScore);
						//lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty);
          computeEValueScoreBit(m->elem, allSeqLength/*, lambda[i]*/);
          m->next = n->next; // m and prev stay the same, but n moves to the next node in finalList          
          free(n); // not freeMNode(n), since n->next is deallocated in connectARegion->freeARegion
          n = m->next;
        }
        else {
          // b) if m and n are of different subject types or of the same subject type, but not too close, 
          // then connect them so the "better" interval is extended towards the "worse" one (lower e-value or higher ratio is better)
          //if (m->elem->e_value < n->elem->e_value) { // extend m towards n 
          //if (m->elem->ratio > n->elem->ratio//void computeEValueScoreBit(mNode *n, long long allSeqLength);
          if (betterOverallScore == 1) {
						extendMToN(m, n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, i,
							args, maxCols, seq, mat, allSeqLength, n->elem->startPosQ, fastScore);
							//lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty);
          }
          //else if (m->elem->e_value > n->elem->e_value) { // extend n towards m
          //else if (m->elem->ratio < n->elem->ratio) { // extend n towards m
          else if (betterOverallScore == -1) {
            extendNToM(n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, i, args, 
							maxCols, seq, mat, allSeqLength, m->elem->endPosQ, fastScore);
          }
          // c) if m and n are of different subject types or of the same subject type, but not too close, and they have equal e-values,
          // then split in the middle
          else {
            middlePoint = m->elem->endPosQ + (n->elem->startPosQ - m->elem->endPosQ) / 2;
            // m ends at middlePoint-1 and n starts at middlePoint
						extendMToN(m, n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, i, args,
							maxCols, seq, mat, allSeqLength, middlePoint, fastScore);
							// lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty);
						extendNToM(n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, i, args,
							maxCols, seq, mat, allSeqLength, middlePoint - 1, fastScore);
							//lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty);
          }
          //if (prev == NULL) { / this could not happen, since in part (1) there is no discarding of m
          //  *finalList[i] = m;
          //}
          prev = m;
          m = n;
          n = n->next;
        } // end else
      } // end if (1)

      // (2) if m and n are overlapping 
      // a) stretch the "better" interval (lower e-value) over an overlapping segment
      // and shorten the "worse" one; if the worse one is now shorter <= arg->f, discard it (done in function extendMToNOverlap/extendNToMOverlap)
      // b) if they have equal e-values, try to split the overlapping segment somewhere around the middle, i.e. so that both segments include their closest exact matches around middle
      // (Note: m and n are of different subject types, since otherwise they would have already been connected)
      else if (lenApartQ < 0) {
        //if (m->elem->e_value < n->elem->e_value) { // extend m towards n (take care of m, prev and n in function)
        //if (m->elem->ratio > n->elem->ratio) { // extend m towards n (take care of m, prev and n in function)
        if (betterOverallScore == 1) {
					extendMToNOverlap(&m, &n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, i,
						args, maxCols, seq, mat, allSeqLength, &prev, fastScore);
						//lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty);
        }
        //else if (m->elem->e_value > n->elem->e_value) { // extend n towards m (take care of m, prev and n in function)
        //else if (m->elem->ratio < n->elem->ratio) { // extend n towards m (take care of m, prev and n in function)
        else if (betterOverallScore == -1) {
					extendNToMOverlap(&m, &n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, i,
						args, maxCols, seq, mat, allSeqLength, &prev, finalList, fastScore);
						//lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty);
        }
        else { // if e-values are equal, split around the middle
					splitMNMiddle(&m, &n, strandBorders, leftBorders, first, second, firstAligned, secondAligned, numOfQueries, args,
						maxCols, seq, mat, allSeqLength, &prev, i, finalList, fastScore);
						//lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty);
          // what about checking whether m and/or n are shorter than args->f or have lower ratio then threshold? --> and then deleting them from the list  - done in splitMNMiddle
        }
      }
      else { // they are exactly on border of each other or a gap is too long, do nothing; move to next pair of neighboring intervals
        //if (prev == NULL) { // this could not happen here, since there is no discarding of m in this part (only in part 2)
        //  *finalList[i] = m;
        //}
        prev = m;
        m = n;
        n = n->next;
      }
      // end while
    }
  }
}


/* concatenate elements of the same subject type if less than maxLenApart (args->C/args->f ?) apart; 
 * --> final merged list with overlapping
*/
void concatenateFinalList(mNode **finalList, Int64 numOfQueries, Args *args, char *seq, Int64 *leftBorders, Int64 *strandBorders,
                      char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, 
											element *mat, long long allSeqLength, Int64 numOfSubjects, 
											double *lambdaAll, float **scoreMatrixAll, int alphabetSize, float ambCharPenaltyAll, 
											char fastScore, double *arrayK) {

  mNode *m = NULL, *n = NULL, *prev = NULL;
  Int64 i = 0, j;
  Int64 lenApartQ = 0; // how far apart are two neighboring intervals on query part 
  Int64 lenApartS = 0; // how far apart are two neighboring intervals on subject part
	Int64 maxLenApart = 0; // args->f * 2;
  
  mNode *prevAll = NULL; // previous to current element (m), considering all elements in finalList
  // auxiliary list of pointers to mNode; j-th pointer points to the previous element of Sj in finalList
  mNode **prevSubject = (mNode **)/*e*/malloc(sizeof(mNode *) * numOfSubjects);
	char sameStrand = 0;

  for (i = 0; i < numOfQueries; i++) {    
		setQueryLambdaMatrix(lambdaAll, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, i, arrayK);

    m = finalList[i];
    for (j = 0; j < numOfSubjects; j++) {
      prevSubject[j] = NULL;
    }

    while (m != NULL && m->elem->endPosQ + leftBorders[i] <= strandBorders[i]) {
      n = m->next; // successor of m
      prev = prevSubject[m->subjectId];
      if (prev != NULL) {        
        // if m and prev (both of the same subject type) are both Q and S part at most maxLenApart apart, then connect them into a single interval
        lenApartQ = m->elem->startPosQ - prev->elem->endPosQ - 1;
        lenApartS = m->elem->startPosS - prev->elem->endPosS - 1;
				// both the same strand
				if ((m->elem->startPosS + leftBorders[m->subjectId] < strandBorders[m->subjectId]
					&& prev->elem->endPosS + leftBorders[m->subjectId] < strandBorders[m->subjectId]) ||
					(m->elem->startPosS + leftBorders[m->subjectId] >= strandBorders[m->subjectId]
					&& prev->elem->endPosS + leftBorders[m->subjectId] >= strandBorders[m->subjectId])) {
					sameStrand = 1;
				}
				else {
					sameStrand = 0;
				}
						// both subject regions should be on the same strand, either reverse or forward
				if (sameStrand && lenApartQ <= maxLenApart && lenApartQ >= 0 && lenApartS <= maxLenApart && lenApartS >= 0) {
          // connect m->start (--> r2) and prev->end (--> tail) and compute alignment between them
					if (connectARegion2(prev->elem, m->elem, args, seq, leftBorders, m->subjectId, i, numOfQueries, firstAligned,
						secondAligned, maxCols, first, second, mat, allSeqLength, fastScore)) {
						//lambda[i], &scoreMatrix[i][0], alphabetSize, ambCharPenalty)) { // connect only if alignment score is positive
						//computeEValueScoreBit(prev->elem, allSeqLength, lambda[i]); // done in connectARegion2             
						// prev points to m->next --> done in connectARegion2
						// m->elem is resolved/deallocated in connectARegion2->freeARegionPartly
						prevAll->next = m->next;
						if (m == finalList[i]) {
							m = m;
						}
						free(m); // not freeMNode(m) since freeARegion(m->elem) is done in freeARegionPartly
						m = NULL;
					}
        }
        else {
          prevSubject[m->subjectId] = m;
        }
      }
      else {
        prevSubject[m->subjectId] = m;
      }
      if (m) {
        prevAll = m;
      }
      m = n;
    } // end while
    if (prevAll) {
      prevAll->next = NULL;
    }
    
    // free rest of the list -- reverse strand
    if (m == finalList[i]) {
      m = m;
    }
    while (m) { 
      n = m->next;
      if (m == finalList[i]) {
        finalList[i] = NULL;
      }
      freeMNode(m);
      m = n;
    } // end while
    
  } // end for
  free(prevSubject);
}


// count number of nucleotides in seq; gaps not included
int getNumNucleotides(char *seq) {  
  int retValue = 0;

  while (*seq) {
    if (*seq != '-') {
      ++ retValue;
    }
    ++ seq;
  }
  return retValue;
}

/* find breakpoints between 2 nodes of different subject types */
void findBP(mNode **finalList, Int64 numOfQueries, Args *args, char *seq, Int64 *leftBorders, Int64 *strandBorders,
                      char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat, 
											long long allSeqLength, Int64 numOfSubjects, int maxGap, 
											double *lambdaAll, float **scoreMatrixAll, int alphabetSize, float ambCharPenaltyAll, 
											char fastScore, Int64 *strandLen, double *arrayK) {

  Int64 i;
  mNode *m = NULL, *n = NULL;
  //mNode *prev = NULL;
  Int64 mStartQ = 0, mStartS = 0, len = 0;
  Int64 nStartQ = 0, nStartS = 0;
  char seqRec1[MAXGAP * 2 + 1], seqS1[MAXGAP * 2 + 1]; // !!!!!!!!!!!!!!!!!!!
  char seqRec2[MAXGAP * 2 + 1], seqS2[MAXGAP * 2 + 1];
  Int64 bpStart, bpEnd;
  Int64 cntAll1, cntAll2;
  Int64 posStart = 0, cntNonGap = 0, cntNonGapR1 = 0;
  double scoreSm = 0.0, scoreSn = 0.0;
  Int64 lenQ = 0;

  for (i = 0; i < numOfQueries; i++) {    
		setQueryLambdaMatrix(lambdaAll, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, i, arrayK);

    m = finalList[i];
    if (m == NULL) {
      continue;
    }
    //prev = NULL;
    while (m->next != NULL /*&& m->elem->endPosQ + leftBorders[i] <= strandBorders[i]*/) { // strand borders were resolved in concatenateFinalList
			n = m->next;
			if (m->subjectId != n->subjectId) { // when neigboring segments belong to different subjects and are closer than 2 * args->f
				// find BP between m and m->next; look into segments of size DIST starting from m's end and check for each segment whether m or m->next wins
				mStartQ = m->elem->end->lbQ + m->elem->end->sl - 1; // starting position of a query part that needs to be aligned
				mStartS = m->elem->end->lbS + m->elem->end->sl - 1;
				len = n->elem->start->lbQ - mStartQ; // gap length
				nStartQ = n->elem->start->lbQ - len; // Q: always fwd strand
				nStartS = n->elem->start->lbS - len;
				/* return 1 if both start and end for query and subject are on the same strand */
				if ((nStartS < 0 || len <= 1) ||
					(!checkBorders(mStartQ, n->elem->start->lbQ - 1, mStartS, mStartS + len - 1, strandLen, i, 
					m->subjectId + numOfQueries)) ||
					(!checkBorders(mStartQ, n->elem->start->lbQ - 1, nStartS, nStartS + len - 1, strandLen, i, 
					n->subjectId + numOfQueries))) {
					m = n; // improve: take only one side of the border
					continue;
				}

        if (len > 0 && len < maxGap) {
          // potential problem - if gap on subject side goes over strand border
					alignSegment(mStartQ, mStartS, leftBorders, len, len, first, second, firstAligned, secondAligned, seq,
						i, m->subjectId + numOfQueries, mat, args, maxCols, fastScore);

          scoreSm = mat[len * maxCols + len].score;
          if (scoreSm > 0) { // copy only if alignment scores above 0 (since this gap can be very short, no ratio is looked at)
            strcpy(seqRec1, firstAligned);
            strcpy(seqS1, secondAligned);
            m->elem->alignedEnd[0] = (char *)/*e*/realloc(m->elem->alignedEnd[0], strlen(seqRec1) + 1); // allocate max. length
            m->elem->alignedEnd[1] = (char *)/*e*/realloc(m->elem->alignedEnd[1], strlen(seqRec1) + 1);        
          }

					alignSegment(nStartQ, nStartS, leftBorders, len, len, first, second, firstAligned, secondAligned, seq, i,
						n->subjectId + numOfQueries, mat, args, maxCols, fastScore);
						//&scoreMatrix[i][0], alphabetSize, ambCharPenalty);
          scoreSn = mat[len * maxCols + len].score;
          if (scoreSn > 0) { // copy only if alignment scores above 0 (since this gap can be very short, no ratio is looked at)
            strcpy(seqRec2, firstAligned);
            strcpy(seqS2, secondAligned);
            n->elem->alignedStart[0] = (char *)/*e*/realloc(n->elem->alignedStart[0], strlen(seqRec2) + 1); // allocate max. length
            n->elem->alignedStart[1] = (char *)/*e*/realloc(n->elem->alignedStart[1], strlen(seqRec2) + 1);
          }

          // align s1 and s2 to recombinant and count differences --> sort of distance, but not converted to distance since the real distance is not important, just
          // which seq (s1 or s2) is closer to recombinant
          if (scoreSm > 0 && scoreSn > 0) {
            alignRecSmSn1(seqRec2, seqS2, seqRec1, seqS1, &cntAll2, &cntAll1, args, &bpStart, &bpEnd); // find breakpoint bpStart --> first bp within fragment between m and n        
            lenQ = bpStart;
          }
          else if (scoreSm > 0) { // only add fragment at the end of m
            lenQ = len;
          }
          else { // (scoreSn > 0) -  only add fragment at the beginning of n
          }

          //cntNonGapR1 = strncpyNoGapCnt(m->elem->alignedEnd[0], seqRec1, bpStart); // copy bpStart nucleotides; gaps are not included into counting
          if (scoreSm > 0) {
            cntNonGapR1 = strncpyNoGapCnt(m->elem->alignedEnd[0], seqRec1, lenQ); // copy bpStart nucleotides; gaps are not included into counting
            m->elem->endPosQ = m->elem->end->lbQ + m->elem->end->sl - 2 + cntNonGapR1;
            //cntNonGap = strncpyNoGapCnt(m->elem->alignedEnd[1], seqS1, bpStart);
            cntNonGap = strncpyNoGapCnt(m->elem->alignedEnd[1], seqS1, lenQ);
            m->elem->endPosS = m->elem->end->lbS + m->elem->end->sl - 2 + cntNonGap;
            m->elem->scoreEnd = getAlignmentScore(m->elem->alignedEnd[0], m->elem->alignedEnd[1], args);
            // add scoreStart scoreEnd and other statistics /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            m->elem->totalLen += strlen(m->elem->alignedEnd[0]);
            m->elem->score += m->elem->scoreEnd;
          }
          else {
            cntNonGapR1 = 0;
          }
          if (scoreSn > 0) {
            posStart = strncpyStart(n->elem->alignedStart[0], seqRec2, cntNonGapR1); // copy char-s from seqRec2 skipping first cntNonGapR1 nucleotides
            cntNonGap = strncpyNoGapCnt(n->elem->alignedStart[1], seqS2 + posStart, strlen(seqRec2) - posStart);        
            n->elem->startPosQ = m->elem->endPosQ + 1;
            n->elem->startPosS = n->elem->start->lbS - cntNonGap;
            n->elem->scoreStart = getAlignmentScore(n->elem->alignedStart[0], n->elem->alignedStart[1], args);
            n->elem->totalLen += strlen(n->elem->alignedStart[0]);
            n->elem->score += n->elem->scoreStart;          
          }

          //// add scoreStart scoreEnd and other statistics /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          //m->elem->totalLen += strlen(m->elem->alignedEnd[0]);
          //m->elem->score += m->elem->scoreEnd;
          //n->elem->totalLen += strlen(n->elem->alignedStart[0]);
          //n->elem->score += n->elem->scoreStart;          
        }
        else if (len < 0) { // since S1 and S2 are overlapping over exact matches - leave as it is
        }
        m->elem->ratio = m->elem->score / m->elem->totalLen;
        computeEValueScoreBit(m->elem, allSeqLength/*, lambda[i]*/);
        n->elem->ratio = n->elem->score / n->elem->totalLen;
        computeEValueScoreBit(n->elem, allSeqLength/*, lambda[i]*/);
      }
      m = n;
    } // end while
  }

}

/* return 1 if both start and end for query and subject are on the same strand */
int checkBorders(Int64 startQ, Int64 endQ, Int64 startS, Int64 endS, Int64 *strandLen, Int64 i, Int64 j) {
	int retValue = 1;
	if (endQ >= strandLen[i]) {
		retValue = 0;
	}
	else if (startS == strandLen[j] || (endS >= strandLen[j] && startS < strandLen[j])) {
		retValue = 0;
	}
	else if (endS >= 2 * strandLen[j]) {
		retValue = 0;
	}
	return retValue;
}


// copy n char-s from src to dest; count non-gap characters
int strncpyNoGapCnt(char *dest, char *src, Int64 n) {

  int retValue = 0;
  while (*src && n) {
    if (*src != '-') {
      ++ retValue;
    }
    -- n;
    *dest = *src;
    ++ dest;
    ++ src;
  }
  *dest = '\0';
  return retValue;
}

// copy n-chars from src to dest starting at nucleotide start; 
// returns position in src from which copying started
Int64 strncpyStart(char *dest, char *src, Int64 start) {

  int retValue = 0;
  while (*src && start) {
    if (*src != '-') {
      -- start;
    }
    ++ retValue;
    ++ src;
    //++ dest;
  }

  while (*src) {
    *dest = *src;
    ++ dest;
    ++ src;
  }

  *dest = '\0';
  return retValue;
}


// align s1 and s2 to recombinant and count differences --> sort of distance, but not converted to distance since the real distance is not important, just
// which seq (s1 or s2) is closer to recombinant
void alignRecSmSn1(char *seqRec1, char *seqS1, char *seqRec2, char *seqS2, Int64 *cntAll1, Int64 *cntAll2, Args *a, Int64 *bpStart, Int64 *bpEnd) {

  Int64 r1 = 0, r2 = 0;
  Int64 j = 0;
  Int64 cnt1 = 0, cnt2 = 0, cntPlusGap1 = 0, cntPlusGap2 = 0;
  Int64 gap1 = 0, gap2 = 0;
  //Int64 prevEnd = 0;
  int SnWins = 0, SnWinsOld = 0; // when Sn wins, SnWins is 1, otherwise 0
  Int64 newBPStart = -1, newBPEnd = -1;
  Int64 firstBPStart = -1 /*, lastBPEnd = -1*/;
  double p_McNemar = 0.;
  Int64 cntSegment1 = 0, cntSegment2 = 0; // counts over current segment (e.g. several consequent blocks where Sn wins)

  *cntAll1 = *cntAll2 = 0;  
  *bpStart = *bpEnd = -1; // *bpStart when n (1) first wins over m (2), and bpEnd when n last timme wins over m (doesn't matter what happens in between)
  // but this only makes sense when *cntAll1 > *cntAll2, otherwise m is considered winner over the whole segment

  while (seqRec1[r1] && seqRec2[r2]) { // seqRec1 && seqRec2 contain same nucleotides, except maybe for gaps
    // count differences/informative sites: +1 for match, 0 for mismatch, -1 for gap
    while (seqRec1[r1] == '-') { // while r1 has gaps, r2 waits; if S1 is gap and r1 isn't, then there is no skipping
      ++ r1;
      ++ gap1;
      //cnt2 += GAP_PENALTY;
    }
    if (seqRec1[r1] == '\0') {
      break;
    }
    while (seqRec2[r2] == '-') {
      ++ r2;
      ++ gap2;
      //cnt1 += GAP_PENALTY;
    }
    if (seqRec1[r2] == '\0') {
      break;
    }
    // compare informative sites (sites where all are equal are non-informative)
    if (seqRec1[r1] == seqS1[r1] && seqRec2[r2] != seqS2[r2]) {
      if (seqS2[r2] == '-') {
        ++ gap2;
				cnt1 += GAP_PENALTY;
      }
      else {
        ++ cnt1;
      }
    }
    if (seqRec1[r1] != seqS1[r1] && seqRec2[r2] == seqS2[r2]) {
      if (seqS1[r1] == '-') {
        ++ gap1;
				cnt2 += GAP_PENALTY;
      }
      else {
        ++ cnt2;
      }
    }
    ++j;
    if (j % DIST == 0) {
      cntPlusGap1 = cnt1 /*+ gap2*/; // gaps are treated the same as mismatches; 
      // TWO OPTIONS TO BE TESTED: (i) all gaps are treated the same (not just those on the S1 or S2 side, but also on R side), and (ii) gaps are treated with extra penalty
      cntPlusGap2 = cnt2 /*+ gap1*/;
      *cntAll1 += cntPlusGap1; // S1 wins on its matches (which are S2's mismatches) and S2's gaps
      *cntAll2 += cntPlusGap2;
      SnWinsOld = SnWins;

      if (cntPlusGap1 > cntPlusGap2 || (cntPlusGap1 == cntPlusGap2 && SnWinsOld == 1)) { // Sn wins (worse results: if cnt1 >= cnt2)
        SnWins = 1;
      }
      else {
        SnWins = 0;
      } 

      if (SnWins) {
        if (SnWinsOld == 0) { // previous segment: Sm won
          newBPStart = j - DIST; // break point - 1st approximation: at the border of regions the transition took place and up to the end of the interval
          // Note: this bp will make sense only if *cntAll1 > *cntAll2 after all DIST regions have been checked
        }
        newBPEnd = j - 1;
        firstBPStart = (firstBPStart == -1) ? newBPStart : firstBPStart;
        //lastBPEnd = newBPEnd;
      }
      else { // Sm wins
        if (SnWinsOld == 1) { // previous segment was Sn
          if (newBPEnd - newBPStart + 1 >= a->f) { 
            p_McNemar = getMcNemar(cntSegment1, cntSegment2); // old values of cntSegment1 and cntSegment2 (up to the last block)
            //p_McNemar = getMcNemarNormalized(cntSegment1, cntSegment2, a->R * (newBPEnd - newBPStart + 1));
            if (p_McNemar >= THRESHOLD_McNEMAR) { // p <= 0.05
              if (*bpStart == -1) { // first Sn segment over complete overlapping region
                *bpStart = newBPStart;
              }
              *bpEnd = newBPEnd; // stretch Sn segment (connect to previous Sn segment)
            }
          }
          //newBPStart = newBPEnd = -1;
        } // end if
      } // end if
      if (SnWins != SnWinsOld) {
        cntSegment1 = cntSegment2 = 0;
      }
      cntSegment1 += cntPlusGap1;
      cntSegment2 += cntPlusGap2;
      cnt1 = cnt2 = 0;
      gap1 = gap2 = 0;
      //prevEnd = j - 1;
    }
    ++ r1;
    ++ r2;
  } // end while
  
  // check last segment (which can be shorter than DIST)
  cntPlusGap1 = cnt1 /*+ gap2*/;
  cntPlusGap2 = cnt2 /*+ gap1*/;
  if (SnWins) { // only if last block was won by Sn - then it makes sense to add this small part
    cntSegment1 += cntPlusGap1;
    cntSegment2 += cntPlusGap2;  
    p_McNemar = getMcNemar(cntSegment1, cntSegment2);
    //p_McNemar = getMcNemarNormalized(cntSegment1, cntSegment2, a->R * (newBPEnd - newBPStart + 1));
    newBPEnd = j - 1; // or j????
    if (newBPEnd - newBPStart + 1 >= a->f || newBPStart == 0) { // long enough region or Sn wins over whole overlapping region
      if (p_McNemar >= THRESHOLD_McNEMAR) { // p < 0.05
        if (*bpStart == -1) { // first Sn segment over long enough overlapping region
          *bpStart = newBPStart;
        } // else - stretch Sn segment (connect to previous Sn segment)
        *bpEnd = newBPEnd; // stretch Sn segment (connect to previous Sn segment)
      }
    }
  }
  else if (SnWins == 0 && *cntAll1 > *cntAll2) { // Sm wins over the last whole block 
  // --> check if Sn won over some previous segment; if not, check whether it wins over whole overlapping region
    if (*bpStart == -1) {
      *cntAll1 += cntPlusGap1;
      *cntAll2 += cntPlusGap2;  
      p_McNemar = getMcNemar(*cntAll2, *cntAll2);
      //p_McNemar = getMcNemarNormalized(*cntAll2, *cntAll2, a->R * (newBPEnd + 1));
      if (p_McNemar >= THRESHOLD_McNEMAR) { // p < 0.05
        *bpStart = 0; // ????????????????? or firstbpstart?
        *bpEnd = newBPEnd; // alternative: from the beginning to the end of the overlapping region
      }
    }
  }

#if DEBUG
  if (fpout) {
    fprintf(fpout, "Breakpoints ==> start: %lld end: %lld\n", (long long) *bpStart, (long long) *bpEnd);
    //fprintf(fpout, "Number of changes: %d\n", numChanges);
    fprintf(fpout, "Final score - S1 aligned to R: %lld\n", (long long) *cntAll1);
    fprintf(fpout, "Final score - S2 aligned to R: %lld\n\n", (long long) *cntAll2);
  }
#endif
}


double getMcNemar(double b, double c) {
  return (b - c) * (b - c) / (b + c);
}

double getMcNemarNormalized(double b, double c, double r) {
  return (b - c) * (b - c) / (b + c) / r;
}

/* compute alignment score of two aligned sequences */
double getAlignmentScore(char *q, char *s, Args *args) {

  int j;
  double score = 0;
  int len = strlen(q);

  for (j = 0; j < len; j++) {
    if (q[j] == '-' || s[j] == '-') {
      score += args->G;
    }
    else if (q[j] == s[j]) {
      score += args->M;
    }
    else {
      score += args->S;
    }
  }
  return score;
}


/* refine alignment score and E-value using new score matrix and lambda */
void refineAlignment(mNode **finalList, Int64 numOfQueries, char *seq, Int64 *leftBorders, Int64 *strandBorders,
	Int64 *strandLen, float **scoreMatrixAll, int alphabetSize, float ambCharPenaltyAll, double *lambdaAll, int maxCols, 
	char *first, char *second, char *firstAligned, char *secondAligned, element *mat, Args *arguments, long long allSeqLength,
	char fastScore, double *arrayK) {
	// calcualte: (1) refined alignment scores using scoreMatrix and (2) E-value based on scoreMatrix and lambda
	Int64 i, j;
  mNode *m = NULL;
	aRegion *ar = NULL;
	aNode *an = NULL;
	double score = 0.0;
	short nIndex = 0;
	for (i = 0; i < numOfQueries; i++) {
		m = finalList[i];
		if (m == NULL) {
			continue;
		}
    setQueryLambdaMatrix(lambdaAll, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, i, arrayK);
    while (m && m->elem->endPosQ + leftBorders[i] <= strandBorders[i]) {
			ar = m->elem;
			an = ar->start;
			score = 0;
			int iAlignedSegment = 0;
			
			while (an) { // scan all aNode-s in m->elem
				// calculate the score of the exact match; it's enough to determine query nucleotide, since the subject nucleotide is its match
				for (j = 0; j < an->sl - 1; j++) {
					nIndex = strchr(strNucl, seq[leftBorders[i] + an->lbQ + j]) - strNucl;
					//score += scoreMatrixAll[i][nIndex * alphabetSize + nIndex];
					score += scoreMatrix[nIndex][nIndex];
				}	
				// calculate the score of aligned segment between two exact matches
				if (iAlignedSegment < ar->numAlignedSegment) {
					Int64 mLen = an->next->lbQ - (an->lbQ + an->sl - 2) - 1;
					Int64 nLen = an->next->lbS - (an->lbS + an->sl - 2) - 1;
					Int64 min = 1;
					if (mLen <= 0 || nLen <= 0) {
						if (mLen <= 0 && nLen <= 0) {
							min = mLen < nLen ? mLen : nLen;
						}
						else {
							min = mLen <= 0 ? mLen : nLen;
						}
						if (min == mLen) {
							an->next->lbQ = an->lbQ + an->sl - 1;
							an->next->sl -= (-mLen /*+ 1*/);
							an->next->lbS += (-mLen /*+ 1*/);
						}
						else if (min == nLen) {
							an->next->lbS = an->lbS + an->sl - 1;
							an->next->sl -= (-nLen /*+ 1*/);
							an->next->lbQ += (-nLen /*+ 1*/);
						}
					}
					mLen = an->next->lbQ - (an->lbQ + an->sl - 2) - 1;
					nLen = an->next->lbS - (an->lbS + an->sl - 2) - 1;
					alignSegment(an->lbQ + an->sl - 1, an->lbS + an->sl - 1, leftBorders,
						mLen, nLen, first, second, firstAligned, secondAligned, seq, i, m->subjectId + numOfQueries, mat,
						arguments, maxCols, fastScore);
						//&scoreMatrix[i][0], alphabetSize, ambCharPenalty);
					ar->alignedSegment[iAlignedSegment * 2] = 
						(char *)/*e*/realloc(ar->alignedSegment[iAlignedSegment * 2], strlen(firstAligned) + 1);
					ar->alignedSegment[iAlignedSegment * 2 + 1] = 
						(char *)/*e*/realloc(ar->alignedSegment[iAlignedSegment * 2 + 1], strlen(secondAligned) + 1);
					strcpy(ar->alignedSegment[iAlignedSegment * 2], firstAligned);
					strcpy(ar->alignedSegment[iAlignedSegment * 2 + 1], secondAligned);
					if (mLen && nLen) {
						ar->scoreAlignedSegment[iAlignedSegment] = mat[mLen * maxCols + nLen].score;
					}
					else {
						ar->scoreAlignedSegment[iAlignedSegment] = mat->score;
					}
					score += ar->scoreAlignedSegment[iAlignedSegment];
					++iAlignedSegment;
				}
				an = an->next;
			} // end while (an)
			ar->score = score;
      ar->totalLen = ar->endPosQ - ar->startPosQ + 1;
			computeEValueScoreBit(ar, allSeqLength/*, lambda[i]*/);
			
			/* Future work:
			 * If ar's E-value has risen and it's not longer significant, then try to shorten ar from the left and/or right
			 * to increase the score, while at the same time keeping the ar's length > args->f 
			 * The problem of finding optimal subregion within ar! */ 
			m = m->next;
		} // end while (m)
	} // end for
}

// extend some regions to left or right towards the breakpoint
void updateBP(mNode **finalList, Int64 numOfQueries, Args *args, Int64 *leftBorders, char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat,
	long long allSeqLength, Int64 numOfSubjects, char *seq, char fastScore, Int64 *strandBorders, Int64 *seqBorders, double *arrayK) {

  mNode *m = NULL, *mStart = NULL, *mEnd = NULL;
  mNode *n = NULL, *p = NULL;
  Int64 i;

  for (i = 0; i < numOfQueries; i++) {
    m = finalList[i];
    if (m == NULL) {
      continue;
    }
    // extend some regions to the rigth, if necessary
    mStart = mEnd = m;
    while (m->next != NULL) {
      n = m->next;
      if (n->elem->end->lbQ + n->elem->end->sl - 2 == m->elem->end->lbQ + m->elem->end->sl - 2 &&
        n->elem->endPosQ != m->elem->endPosQ) { // n was extended to the rigth, so try to extend every region from mStart to mEnd (predecessors of n)
        p = mStart;
        while (p != n) { // check all predecessors of n that are candidates for the extension towards the breakpoint
          int len = strlen(n->elem->alignedEnd[0]) - cntGaps(n->elem->alignedEnd[0]); // strlen without gaps
          double score = 0;
          if (len < maxCols) {
						// check that subject alignment does not cross strand or sequence border
						if ((p->elem->endPosS + len + leftBorders[p->subjectId + numOfQueries] >= strandBorders[p->subjectId + numOfQueries] &&
							p->elem->endPosS + leftBorders[p->subjectId + numOfQueries] < strandBorders[p->subjectId + numOfQueries]) ||
							p->elem->endPosS + len + leftBorders[p->subjectId + numOfQueries] >= seqBorders[p->subjectId + numOfQueries]) {
							// do nothing
						}
						else {
							alignSegment(p->elem->endPosQ, p->elem->endPosS, leftBorders, len, len, first, second, firstAligned, secondAligned, seq,
								i, p->subjectId + numOfQueries, mat, args, maxCols, fastScore);
							score = mat[len * maxCols + len].score;
							if (score > 0) { // copy only if alignment scores above 0 (since this gap can be very short, no ratio is looked at)
								p->elem->alignedEnd[0] = (char *)/*e*/realloc(p->elem->alignedEnd[0], strlen(firstAligned) + 1); // allocate max. length
								p->elem->alignedEnd[1] = (char *)/*e*/realloc(p->elem->alignedEnd[1], strlen(secondAligned) + 1);
								strcpy(p->elem->alignedEnd[0], firstAligned);
								strcpy(p->elem->alignedEnd[1], secondAligned);
								p->elem->scoreEnd = score; // getAlignmentScore(m->elem->alignedEnd[0], m->elem->alignedEnd[1], args);
								// add scoreStart scoreEnd and other statistics /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
								p->elem->totalLen += strlen(p->elem->alignedEnd[0]);
								p->elem->score += p->elem->scoreEnd;
								p->elem->endPosQ = n->elem->endPosQ;
								p->elem->endPosS = n->elem->endPosS;
								computeEValueScoreBit(p->elem, allSeqLength);
							}
						}
          }
          p = p->next;
        } // end while
      }
      else if (n->elem->endPosQ == m->elem->endPosQ) { // if both m and n end at the same position
        mEnd = n;
      }
      else { // new block
        mStart = mEnd = n;
      } // end if
      m = m->next;
    } // end while

    m = finalList[i];
    while (m && m->next) { // extend to the left
      n = m->next;
      if (n->elem->start->lbQ == m->elem->start->lbQ && n->elem->startPosQ != m->elem->startPosQ) {
        // m was extended to the left, so try to extend every region n and its successors towards the breakpoint
        while (n && n->elem->start->lbQ == m->elem->start->lbQ) {
          int len = strlen(m->elem->alignedStart[0]) - cntGaps(m->elem->alignedStart[0]); // number of no gap positions in m->elem->alignedStart[0]
          double score = 0;
          if (len < maxCols) {
						// check that on the subject side the strand or seq border is
						if (n->elem->start->lbS - len < 0 ||
							(n->elem->start->lbS + leftBorders[n->subjectId + numOfQueries] >= strandBorders[n->subjectId + numOfQueries] &&
							n->elem->start->lbS - len + leftBorders[n->subjectId + numOfQueries] < strandBorders[n->subjectId + numOfQueries])) {
							// do nothing
						}
						else {
							alignSegment(n->elem->start->lbQ - len, n->elem->start->lbS - len, leftBorders, len, len, first, second, firstAligned, secondAligned, seq,
								i, n->subjectId + numOfQueries, mat, args, maxCols, fastScore);
							score = mat[len * maxCols + len].score;
							if (score > 0) { // copy only if alignment scores above 0 (since this gap can be very short, no ratio is looked at)
								n->elem->alignedStart[0] = (char *)/*e*/realloc(n->elem->alignedStart[0], strlen(firstAligned) + 1); // allocate max. length
								n->elem->alignedStart[1] = (char *)/*e*/realloc(n->elem->alignedStart[1], strlen(secondAligned) + 1);
								strcpy(n->elem->alignedStart[0], firstAligned);
								strcpy(n->elem->alignedStart[1], secondAligned);
								n->elem->scoreStart = score; // getAlignmentScore(m->elem->alignedEnd[0], m->elem->alignedEnd[1], args);
								// add scoreStart scoreEnd and other statistics /// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
								n->elem->totalLen += strlen(n->elem->alignedStart[0]);
								n->elem->score += n->elem->scoreStart;
								n->elem->startPosQ = m->elem->startPosQ;
								n->elem->startPosS = m->elem->startPosS;
								computeEValueScoreBit(n->elem, allSeqLength);
							} // end if
						}
          }
          n = n->next;
        } // end while
        m = n;
      }
      else {
        m = m->next;
      }
    } // end while

  } // end for
}

// coount gaps in string
int cntGaps(char *str) {
  int cnt = 0;
  while (*str) {
    if (*str == ' ') {
      ++cnt;
    }
    ++str;
  }
  return cnt;
}

