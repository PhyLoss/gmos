/***** resolveOverlaps.c *************************************************************
* Description: Source file for resolving overlaps between aligned segments.
*
* Author: Mirjana Domazet-Loso
*
* This file is part of gmos.
*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
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
#include "resolveOverlaps.h"

/* ************ check all shorten-ings, and at the end, delete all remaining mNodes whose e-value is > EValueLow */
/* resolve overlaps remained after concatenation */
void resolveOverlapsFinalList(FILE *fpout, mNode ***finalList, Int64 numOfQueries, Args *args, char *seq, Int64 *leftBorders, Int64 *strandBorders,
	char *firstAligned, char *secondAligned, int maxCols, char *first, char *second, element *mat,
	long long allSeqLength, Int64 numOfSubjects, double *lambdaAll, float **scoreMatrixAll, int alphabetSize, 
	float ambCharPenaltyAll, char fastScore, double *arrayK) {

  Int64 i = 0;
  mNode *m = NULL, *n = NULL, *prev = NULL, *p = NULL;
	Int64 bpStart = -1, bpEnd = -1;
	int leaveM = 0; // 1 - if m stays, 0 otherwise
	int leaveLeftPart = 0;
	int minSegmentLen = args->f;
  
  mNode *leftMost = NULL; /* the leftmost mNode in the finalList that overlaps with n; the leftmost meaning with the smallest value of startPosQ */
  //mNode *oldLeftMost = NULL; 
  mNode *rightMost = NULL; /* the rightmost mNode in the finalList that overlaps with n; the rightmost meaning with the highest value of endPosQ  */
  mNode *curr = NULL, *curr2 = NULL, *currNext = NULL;
  mNode *prevLM = NULL, *prevCurr = NULL;
  int leaveCurr = 0; // 1 - if curr stays, 0 otherwise

  /*********************** initialize ************************/
  gapPenalty = getGapPenalty();
  EValueLow = args->E;
	//excellentEValue = EValueLow * EValueLow;
	excellentEValue = args->X;

	for (i = 0; i < numOfQueries; i++) {
		m = (*finalList)[i];
		if (m == NULL) {
			continue;
		}

    /* set query specific lambda and scoreMatrix */
		setQueryLambdaMatrix(lambdaAll, scoreMatrixAll, alphabetSize, ambCharPenaltyAll, i, arrayK);

		prev = NULL;
    prevLM = NULL;
    leftMost = rightMost = m;
		while (m->next != NULL && m->elem->endPosQ + leftBorders[i] <= strandBorders[i]) { // strand borders were resolved in concatenateFinalList
			n = m->next;
			// if n covers exactly the same query region in the ocmpletely same alignment as m, do nothing; i.e. go to the next 
			if (n->elem->startPosQ == m->elem->startPosQ && n->elem->startPosQ == m->elem->startPosQ &&
				compareEValues(m->elem->e_value, n->elem->e_value)) {
				prev = m;
				m = n;
				continue;
			}

			//oldLeftMost = leftMost;
      findLeftMost(&leftMost, &prevLM, n->elem->startPosQ); /* find new leftmost that includes startPosQ */
      if (leftMost) {
        findRightMost(leftMost, n, &rightMost); /* find new rightMost node up to n */
      }
      else {
        rightMost = NULL;
      }

      if ((!rightMost /*n->elem->startPosQ > rightMost->elem->endPosQ*/) || // CASE 0: n and rightMost do not overlap, keep them both - or -
        (rightMost->subjectId == n->subjectId)) { // the same subject type, so keep them both, i.e. do nothing
        prev = m;
      	m = n;
      } // end CASE 0
      /********************************************************************************************************************************/

      else if (rightMost->elem->startPosQ <= n->elem->startPosQ && n->elem->endPosQ <= rightMost->elem->endPosQ) { // CASE 1: n within rightMost
				// check rightMost's (before: m's) score over the overlapping region; if both score better than E-valueLow, keep them both;
				// otherwise, keep the better e-value
        double scoreOverlap = getScoreSeg(rightMost, seq, i, leftBorders, n->elem->startPosQ, n->elem->endPosQ, gapPenalty);
				double eValueRM = getEValueScoreBit(scoreOverlap, (double)n->elem->endPosQ - n->elem->startPosQ + 1, allSeqLength);
        if (compareEValues(eValueRM, n->elem->e_value)) { // CASE 1a) keep them both (if scores are equally good or both scores are excellent)
					prev = m;
					m = n;
				}
        else if (eValueRM < n->elem->e_value) { // CASE 1b) if rightMost (previously: m) scores better (lower E-value) than n over the overlapping region, free n
					m->next = n->next;
					freeMNode(n);
					// m stays the same 			
				}
				else { // CASE 1c) if n scores better over the overlapping region, shorten all mNodes from leftmost to m from left and right ************************
					bpStart = n->elem->startPosQ;
					bpEnd = n->elem->endPosQ;
					leaveLeftPart = 0;
          curr = leftMost;
          prevCurr = prevLM;
					char currIsLeftMost = 1;
          while (curr) {
						currNext = curr->next;
						if (curr != leftMost) { // leftmost and rightmost --- redundant??
							currIsLeftMost = 0;
						}
            if (curr->elem->endPosQ >= n->elem->startPosQ) {
              leaveCurr = shortenM(&curr, n, bpStart, bpEnd, finalList, i, args, allSeqLength, &curr2, leaveLeftPart, seq, leftBorders);
              /* new scores and evalues for curr and curr2 (previously: m and m2) are computed in shortenM --> freeARegionFromMToN */
              // Note: curr2 (right part of curr) is connected to its predecessor/successor (whom ever they may be) in shortenM ???????????????
              if (leaveCurr) { // left part of curr stays; it will still point to its old successor
                prevCurr = curr;
                if (prevCurr == leftMost) { // leftMost may not be the true leftMost any more, but it doesn't matter, it wll be updated later
								}
							}
              else { // if leaveCurr == 0; its predecessor should be connected to its successor
								if (currIsLeftMost) { // if curr was also leftMost (rightMost is not important here, since it will be udpated wth the next m)
                  leftMost = currNext;
									if (prevLM) {
										prevLM->next = currNext;
									}
                }
                connectPrev(leaveCurr, NULL, currNext, &prevCurr, finalList, i); // connect prevCurr to curr->next
              } // end else
            }
            else {
              prevCurr = curr;
            }
            curr = currNext;
            if (curr == n) {
              m = n;
              break;
            }
          } // end while
					////////////////
				}
			} // end CASE 1
      /********************************************************************************************************************************/

      else if (rightMost->elem->startPosQ <= n->elem->startPosQ && rightMost->elem->endPosQ < n->elem->endPosQ) { // CASE 2: m and n overlap
				// if they score equally well over the overlapping region, keep them both;
				// otherwise, shorten the mNode which has lower score and keep it if its shorten part is longer than minSegmentLen
        double scoreOverlap1 = getScoreSeg(rightMost, seq, i, leftBorders, n->elem->startPosQ, rightMost->elem->endPosQ, gapPenalty); // before: m instead of rightMost
				double scoreOverlap2 = getScoreSeg(n, seq, i, leftBorders, n->elem->startPosQ, rightMost->elem->endPosQ, gapPenalty);
        Int64 segLen = rightMost->elem->endPosQ - n->elem->startPosQ + 1; // before: m instead of rightMost
				double eValue1 = getEValueScoreBit(scoreOverlap1, (double)segLen, allSeqLength);
				double eValue2 = getEValueScoreBit(scoreOverlap2, (double)segLen, allSeqLength);
				//************************
				if ((scoreOverlap1 == scoreOverlap2 && segLen < minSegmentLen)
          || (compareEValues(eValue1, eValue2) && segLen >= minSegmentLen)) { // CASE 2a)
          // keep both m and n (for short overlaps: if scores are equal; for long overlaps: e-values are low)
					prev = m;
					m = n;
        }
				else if ((scoreOverlap1 < scoreOverlap2 && segLen < minSegmentLen)
							|| (eValue1 > eValue2 && segLen >= minSegmentLen)) { // // CASE 2b) n wins; shorten all nodes from leftmost to m up to n ?????????????????
					bpStart = n->elem->startPosQ;
					bpEnd = rightMost->elem->endPosQ + 1; // right part stays; before: m instead of rightMost
					leaveLeftPart = 0;
          curr = leftMost;
          prevCurr = prevLM;
					char currIsLeftMost = 1;
					while (curr && curr != n) { /*********************/
						currNext = curr->next;
						if (curr != leftMost) { // leftmost and rightmost --- redundant??
							currIsLeftMost = 0;
						}
						if (curr->elem->endPosQ >= bpStart) {
              leaveM = shortenM(&curr, n, bpStart, bpEnd, finalList, i, args, allSeqLength, &curr2, leaveLeftPart, seq, leftBorders);
              if (leaveM) { // curr->elem->e_value already checked in shortenM
                if (curr->next == n) {
                  connectPrev(1, m, NULL, &prev, finalList, i); // here - m is curr? redundant?
                }
                prevCurr = curr;
              }
							else { //(leaveM == 0) {
								connectPrev(leaveM, curr, currNext, &prevCurr, finalList, i);
								if (currIsLeftMost) { // when leftMost was deleted
									leftMost = currNext;
									if (prevLM) {
										prevLM->next = currNext;
									}
								}
							}
							curr = currNext;
              if (curr == n) {
                m = n;
                break;
              }
            }
            else {
              prevCurr = curr;
              curr = currNext;
            }
          } // end while
				}
				else { // CASE 2c
					// shorten n from the left side (n should start after rightMost's right border)
					// shorten n so it covers the region from start to end; start is the first position within n that stays
					p = n->next;
          shortenN(&n, rightMost->elem->endPosQ, n->elem->endPosQ, allSeqLength, args, leftBorders, i, seq); // before: m instead of rightMost
					if (n) { // n should be kept if e-value is still low, i.e. != null 
						connectNNext(&m, n); // connect m and n->next or n, since n doesn't have to be directly connected to m anymore
						if (m->next == n) {
							connectPrev(1, m, NULL, &prev, finalList, i);
							m = n;
						}
					}
					if (!n) { // when n is deleted (too short), m stays as it is
						m->next = p;
					}
				} // end CASE 2c
			} // end CASE 2
      /********************************************************************************************************************************/
		} // end while
	} // end for

}


/* when n is within m; compare their e-values over n; if both score < EValueLow, then keep them both;
 * else keep the better score (lower e-value) */
double getScoreSeg(mNode *m, char *seq, Int64 iQuery, Int64 *leftBorders, Int64 startSeg, Int64 endSeg, double gapPenaltyLocal) {
	aRegion *arM = m->elem;
	aNode *am = arM->start;
	int i = 0;
	double score = 0;
	char found = 0;
	Int64 amEnd, end;
	while (am) {
		amEnd = am->lbQ + am->sl - 2;
		if (am->lbQ <= startSeg && startSeg <= amEnd) { // (case 1) found startSeg
			end = amEnd < endSeg ? amEnd : endSeg; // end of exact match to be included in the score
			score += scoreMatch(startSeg, end, iQuery, seq, leftBorders);
			// check next aligned region (between two matches)
			//if (amEnd < endSeg) {
			//	Int64 len = am->next->lbQ < endSeg ? strlen(arM->alignedSegment[2 * i]) : (endSeg - amEnd);
			//	/* calculate score between two aligned regions using scoreMatrix*/
			//	score += scoreAlignedRegion(len, arM->alignedSegment[2 * i], arM->alignedSegment[2 * i + 1], gapPenaltyLocal);
			//}
			found = 1; // found the beginning of m 
		}
		if ((found && amEnd < endSeg) // check next aligned region (between two matches)
			|| (!found && i < arM->numAlignedSegment && amEnd < startSeg && startSeg < am->next->lbQ)) {
			Int64 len = endSeg < am->next->lbQ ? (endSeg - amEnd) : strlen(arM->alignedSegment[2 * i]);
			end = endSeg < am->next->lbQ ? endSeg : am->next->lbQ - 1;
			/* calculate score between two aligned regions using scoreMatrix*/
			score += scoreAlignedRegion(len, arM->alignedSegment[2 * i], arM->alignedSegment[2 * i + 1], gapPenaltyLocal);
			found = 1; // found the beginning of m 
		}
		++i; // next pair of segments
		am = am->next; // next match
		if (found) { 
			break; // continue score computation in the next while-loop
		}
	} // end while

	// add to score every match and aligned region up to segEnd;
	// in both cases (1 and 2): continue calculation at am
	while (end < endSeg && am) {
		amEnd = am->lbQ + am->sl - 2;
		end = amEnd < endSeg ? amEnd : endSeg; // exact match to be included in the score
		score += scoreMatch(am->lbQ, end, iQuery, seq, leftBorders);
		if (amEnd >= endSeg) { // the end of the overlap
			break;
		}
		else if (i < arM->numAlignedSegment) { // aligned region (between two matches)
			Int64 len = am->next->lbQ < endSeg ? strlen(arM->alignedSegment[2 * i]) : (endSeg - amEnd);
			/* calculate score between two exact matches using scoreMatrix*/
      score += scoreAlignedRegion(len, arM->alignedSegment[2 * i], arM->alignedSegment[2 * i + 1], gapPenaltyLocal);
		}

		if (am->next && am->next->lbQ > endSeg) { // the end of the overlap
			break;
		}
		++i; // next pair of segments
		am = am->next;
	} // end while
	return score;
}

/* compute e-value */
double getEValueScoreBit(double score, double queryLen, long long allSeqLength) {
	double score_blast = (lambda * score - log(DEFAULT_K)) / log(2.0); // S' = (lambda x S - ln K) / ln 2
	double eValue = queryLen * (double)allSeqLength * pow(2.0, -score_blast); // e-value = mn x 2^(-S')
	return eValue;
}


// shorten m up to bpStart and shorten m from bpEnd
// if any of these segments too short, delete it, otherwise - keep shorten part;
// shortenM returns 1 if m stays, and 0 otherwise
int shortenM(mNode **m, mNode *n, Int64 bpStart, Int64 bpEnd, mNode ***finalList, Int64 iQuery, Args *a,
  long long allSeqLength, mNode **m2, int leaveLeftPart/*, double lambda*/, char *seq, Int64 *leftBorders) {

  aNode *mm = (*m)->elem->start, *mprev = NULL;
  aNode *nn = NULL, *nprev = NULL;
  aRegion *p2 = NULL;
  int mmAlignedSegments = 0;
  int nnAlignedSegments = 0;
  Int64 lenLeft = 0, lenRight = 0;
  mNode *tprev = NULL;
  int retValue = 1; // 1 if m stays, and 0 otherwise
  int minSegmentLen = (int)(a->r * a->f);
  double m_score = 0.0;
  //double m_eValue = 0;

  *m2 = NULL;
  //m_score += scoreMatch(m->elem->start->lbQ, m->elem->start->lbQ + m->elem->start->sl - 2, iQuery, seq, leftBorders);
  while (mm != NULL) { // mprev - last match from the left *before* recStart; mm - first match *in* overlap region
    //if (mm->lbQ + mm->sl - 2 < bpStart) {
    if (mm->lbQ < bpStart) {
      mprev = mm;
      mm = (mm == (*m)->elem->end) ? NULL : mm->next;
      if (mm && mm->lbQ < bpStart) {
        m_score += scoreMatch(mprev->lbQ, mprev->lbQ + mprev->sl - 2, iQuery, seq, leftBorders);
        m_score += (*m)->elem->scoreAlignedSegment[mmAlignedSegments];
      }
      ++mmAlignedSegments;
    }
    else {
      break;
    }
  } // end while
  --mmAlignedSegments; // number of aligned segments within left part of m that needs to be kept

	/// new 16022015
	if (mm && mm->lbQ + mm->sl - 2 > bpEnd) {
		mm = NULL;
		return 1;
	}
	/// new

  nprev = nn = mm;
  while (nn != NULL) { // find nn - first match from the right *after* recEnd; nprev - *last* match in overlap region
    //if (nn->lbQ + nn->sl - 2 > bpEnd) {
    if (nn != mm && nn->lbQ + nn->sl - 2 > bpEnd) {
      break;
    }
    else {
      nprev = nn;
      nn = (nn == (*m)->elem->end) ? NULL : nn->next;
      ++nnAlignedSegments;
    }
  } // end while
  --nnAlignedSegments; // number of aligned segments within right part of m that needs to be kept

  // check mm - left part of m (up to bpStart)
  lenLeft = (mprev != NULL) ? (mprev->lbQ + mprev->sl - 1 - (*m)->elem->start->lbQ) : 0; // length of m's left part (before overlapping)
  lenRight = (nn != NULL) ? ((*m)->elem->end->lbQ + (*m)->elem->end->sl - 1 - nn->lbQ) : 0; // length of m's right part (after overlapping)

  //m_eValue = getEValueScoreBit(m_score, (double)lenLeft, allSeqLength);

  // keep m2 (m's right part) 
  if (lenRight >= minSegmentLen) { // construct m2; nn - starting aNode in m2
    p2 = getARegion(nn, a->M);
    // add all the nodes from nn to end are added later in freeARegionFromMToN and new e-value, startPosQ, endPosQ, etc. are computed there!!!!!
    *m2 = getMNode(p2, (*m)->subjectId);
    // Note: m2 might not be connected directly to n->next, because n->next may come before m2; 
    // so find m2's place in linked list
  }

  if (leaveLeftPart || // m stays; thus: prev -> m -> n -> ... [m2] -> ... **************** this should be better written so m2's score is computed in advance
    (lenLeft >= minSegmentLen /*&& m_eValue < EValueLow*/)) { 
		if (lenRight >= minSegmentLen /*&& (*m2)->elem->e_value < EValueLow*/) { // possibly keep m2 (m's right part) ==> CASE 2
      freeARegionFromMToN((*m)->elem, mm, nprev, allSeqLength, (*m2)->elem, seq, leftBorders, iQuery); // delete middle part of m
			tprev = n;
			connectM2(tprev, m2); // connects m2 to both its predecessor and its successor
    }
    else { // no m2, i.e. m2 is null
      freeARegionFromMToN((*m)->elem, mm, (*m)->elem->end, allSeqLength, NULL, seq, leftBorders, iQuery); // delete m's right part
    } // end if
  }
  else { // (lenLeft < minSegmentLen or new m's e-value is too high, i.e. m is deleted
    retValue = 0;
		if (*m2 && lenRight >= minSegmentLen /*&& (*m2)->elem->e_value < EValueLow*/) { // possibly m2 stays, i.e. prev -> n -> ... -> [m2] -> ..  ==> CASE 3 (m2 is in fact n2_main, and m is n_main)
      freeARegionFromMToN((*m)->elem, (*m)->elem->start, nprev, allSeqLength, (*m2)->elem/*, lambda*/, seq, leftBorders, iQuery); // delete left part of m, up to and including nprev 
			tprev = n;
			connectM2(tprev, m2);
    }
    else { // delete complete m (m2 included)
    }
    freeMNode(*m); // delete complete m (with or without m2)
    *m = NULL;
  }
  return retValue;
}


/* free aNode-s up from m to n; if something remains after n --> construct new aRegion r2 */
void freeARegionFromMToN(aRegion *r, aNode *m, aNode *n, long long allSeqLength, aRegion *r2/*, double lambda*/, char *seq, Int64 *leftBorders, Int64 iQuery) {

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
    //r->score += a->sl - 1;
    r->score += scoreMatch(a->lbQ, a->lbQ + a->sl - 2, iQuery, seq, leftBorders); // new
    if (cntAlignedSegmentsLeft < r->numAlignedSegment) {
      r->score += r->scoreAlignedSegment[cntAlignedSegmentsLeft];
      ++cntAlignedSegmentsLeft;
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
    ++cntAlignedSegments;
  } // end while

  // transfer aligned segments and matches starting with b to m2    
  if (r2 != NULL) {
    //if (b) { // b could be possibly equal to a (when n == m); in that case, skip the following while-loop
    b = b->next; // since b was already added to r2 (before calling this function)
    while (b != NULL) {
      r2->end->next = b;
      r2->end = b;
      //r2->score += b->sl - 1;
      r2->score += scoreMatch(b->lbQ, b->lbQ + b->sl - 2, iQuery, seq, leftBorders); // new

      ++cntAlignedSegmentsRight;
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

int roundNeg(double value) {
	return (int)(value - ceil(value) > -0.5 ? ceil(value) : floor(value));
}

/* returns 1 if e-values are equal or equally good, and otherwise 0 */
int compareEValues(double eValueM, double eValueN) {
	int logM = roundNeg(log10(eValueM));
	int logN = roundNeg(log10(eValueN));
  if (fabs(logM - logN) < 1 ||
    (eValueM <= excellentEValue && eValueN <= excellentEValue)) { // if scores are equally good or both scores are excellent
    return 1;
  }
  else {
    return 0;
  }
}

/* find new leftMost and its predecessor, such that leftMost includes startPosQ */
void findLeftMost(mNode **leftMost, mNode **prevLM, Int64 startPosQ) {
  mNode *curr = *leftMost;
  mNode *prevCurr = *prevLM;
  while (curr) { // find new leftmost
    if (startPosQ > curr->elem->endPosQ) {
      prevCurr = curr;
      curr = curr->next;
    }
    else {
      break;
    }
  } // end while
  if (curr) {
    *leftMost = curr;
    *prevLM = prevCurr;
  }
  else {
    *leftMost = NULL;
    *prevLM = NULL;
  }
}

/* find new rightMost node up to n */
void findRightMost(mNode *leftMost, mNode *n, mNode **rightMost) {
  mNode *curr = leftMost;
  Int64 maxRB = leftMost->elem->endPosQ;
  *rightMost = leftMost;

  while (curr) { // find new rightmost before n
    if (curr == n) {
      break;
    }
    if (n->elem->startPosQ < curr->elem->endPosQ) {
      if (curr->elem->endPosQ > maxRB) {
        maxRB = curr->elem->endPosQ;
        *rightMost = curr;
      }
    }
    curr = curr->next;
  } // end while
}

// shorten n so it covers the region from start to end; start is the first position within n that stays
// (it is either added to alignedStart or aNode is shortened)
void shortenN(mNode **nn, Int64 start, Int64 end, long long allSeqLength, Args *a, Int64 *leftBorders, Int64 iQuery /*, double lambda*/, char *seq) {

  aNode *p = NULL, *p2 = NULL;
  int i, iOld;
  int cntAlignedSegments = 0, keepAlignedSegments = 0;
  mNode *n = *nn;
  int flagNoANode = 0;
  int minSegmentLen = (int)(a->r * a->f);

  //Int64 endPrev = -1;  

  n->elem->score = 0;
  // delete aNode-s and aligned parts before start
  p = n->elem->start;
  while (p->lbQ + p->sl - 2 < start) { // keep aNode which contains start (i.e. final n comprises aNode which covers start)
    //endPrev = p->lbQ + p->sl - 2;
    p2 = p->next;
    free(p);
    p = p2;
    ++cntAlignedSegments;
  } // end while

  if (p->lbQ <= end) { // at least one aNode between start and end, or overlapping start and end
    n->elem->start = p;
    n->elem->startPosQ = n->elem->start->lbQ;
    n->elem->startPosS = n->elem->start->lbS;
    //while (p && p->lbQ + p->sl - 2 <= end) { // keep aNodes from start till end
    p2 = NULL;
    while (p && p->lbQ <= end) { // keep aNodes from start till end
      ///////////////////////////////////////////////////////////////////////n->elem->score += p->sl - 1; 
      n->elem->score += scoreMatch(p->lbQ, p->lbQ + p->sl - 2, iQuery, seq, leftBorders);
      p2 = p;
      p = p->next;
      ++keepAlignedSegments;
    }

    --keepAlignedSegments;
    assert(p2 > 0);
    p2->next = NULL;
    n->elem->end = p2;
    n->elem->endPosQ = n->elem->end->lbQ + n->elem->end->sl - 2;
    n->elem->endPosS = n->elem->end->lbS + n->elem->end->sl - 2;

    // delete aNode-s until end node
    while (p != NULL) {
      p2 = p->next;
      free(p);
      p = p2;
    } // end while
  }
  else {
    n->elem->start = p; // ?
    flagNoANode = 1;
  }

  // if, after deleting of the middle part, n becomes too short
  //if (!n || /*(n->elem->start == n->elem->end)*/ (n->elem->endPosQ - n->elem->startPosQ + 1 < a->f)) { 
  //if (!n || (n->elem->endPosQ - n->elem->startPosQ + 1 < a->f)) { // n should always be ! = null
  if (flagNoANode == 1 || (n->elem && n->elem->endPosQ - n->elem->startPosQ + 1 < minSegmentLen)) {
    freeMNode(n);
    n = NULL;
    *nn = NULL;
  }
  else {
    assert(cntAlignedSegments > 0);
    // delete first cntAlignedSegments segments from n
    if (cntAlignedSegments > 0) { // ?????????????????????
      for (i = 0; i < keepAlignedSegments; i++) { // ??? Are they some problems with overlapping of segments - reallocation ??
        iOld = i + cntAlignedSegments;
        n->elem->alignedSegment[2 * i] = (char *)/*e*/realloc(n->elem->alignedSegment[2 * i], strlen(n->elem->alignedSegment[2 * iOld]) + 1);
        n->elem->alignedSegment[2 * i + 1] = (char *)/*e*/realloc(n->elem->alignedSegment[2 * i + 1], strlen(n->elem->alignedSegment[2 * iOld + 1]) + 1);
        strcpy(n->elem->alignedSegment[2 * i], n->elem->alignedSegment[2 * iOld]);
        strcpy(n->elem->alignedSegment[2 * i + 1], n->elem->alignedSegment[2 * iOld + 1]);
        n->elem->scoreAlignedSegment[i] = n->elem->scoreAlignedSegment[iOld]; // copy old score
        n->elem->score += n->elem->scoreAlignedSegment[i];
      }
    }
    // NOT LIKE THIS! this will be deleted when freeMNode --> in this way: this allocations are lost and never deleted
    //for (; i < n->elem->numAlignedSegment; i++) {
    //  n->elem->alignedSegment[2 * i] = n->elem->alignedSegment[2 * i + 1] = NULL; 
    //}
    n->elem->numAlignedSegment = keepAlignedSegments;

    // update other values
    n->elem->totalLen = n->elem->endPosQ - n->elem->startPosQ + 1;
    n->elem->ratio = n->elem->score / n->elem->totalLen;
    computeEValueScoreBit(n->elem, allSeqLength/*, lambda*/);
    //if (n->elem->e_value >= EValueLow) { // *********** leave out?
    //  freeMNode(n);
    //  n = NULL;
    //  *nn = NULL;
    //}
  }
}

/* connect prev to m or n, depending on wheter m is deleted or not (leaveM == 0 or 1 respectively) */
void connectPrev(int leaveM, mNode *m, mNode *n, mNode **prev, mNode ***finalList, int i) {
  if (leaveM == 1) {
    if (*prev == NULL) {
      (*finalList)[i] = m;
    }
    *prev = m;
  }
  else { // m was deleted
    if (*prev == NULL) {
      (*finalList)[i] = n;
    }
    else {
      (*prev)->next = n;
    }
  }
}

/* connect m2 to tprev (predecessor) and tprev->next (successor) */
void connectM2(mNode *tprev, mNode **m2) {

  mNode *t = NULL;

  for (t = tprev->next; t && t->elem->start->lbQ < (*m2)->elem->start->lbQ;) {
    tprev = t;
    t = t->next;
  }
  (*m2)->next = t;
  tprev->next = *m2;
}

void connectNNext(mNode **m, mNode *n) {

  mNode *p = NULL, *pp = NULL;
  p = (*m)->next = n->next;
  while (p && (p->elem->start->lbQ < n->elem->start->lbQ)) {
    pp = p;
    p = p->next;
  } // end while

  if (p == n->next) {
    (*m)->next = n;
  }
  else {
    pp->next = n;
    n->next = p;
  }
}
