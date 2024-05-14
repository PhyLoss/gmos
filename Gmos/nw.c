/***** nw.c *************************************************************************************************
 * Description: Functions for sequence alignment using Needleman-Wunsch algorithm.
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
 *
 ************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "commonSC.h"
#include "eprintf.h"
#include "stringUtil.h"

#include "interface.h"
#include "nw.h"

#if defined(_DEBUG) && defined(WIN) 
#include "leakWatcher.h"
#endif

#if defined(_DEBUG) && defined(WIN) 
  #define new DEBUG_NEW
  #undef THIS_FILE
  static char THIS_FILE[] = __FILE__;
#endif

/* allocate matrix and set all elements to 0 */
element *getMatrix(int rows, int cols) {  
  element *mat = NULL;
  mat = (element *)/*e*/malloc(sizeof(element) * rows * cols);
  return mat;
}

/* free matrix */
void freeMatrix(element *mat) {
  free(mat);
}

/* calculate match/mismatch score between two aligned nucleotides */
//float getScore(float *scoreMatrix, int alphabetSize, char qNucl, char sNucl, float ambCharPenalty) {
float getScore(char qNucl, char sNucl) {
	float retValue = 0;
	if (qNucl == toupper(AMBIGUOUS_CHAR_QUERY) || qNucl == toupper(AMBIGUOUS_CHAR_SUBJECT)
		|| sNucl == toupper(AMBIGUOUS_CHAR_SUBJECT) || sNucl  == toupper(AMBIGUOUS_CHAR_QUERY)) { // treat ambiguous nucleotides as mismatches
		retValue = ambCharPenalty;
	}
	else {
		int i = strchr(strNucl, qNucl) - strNucl;
		int j = strchr(strNucl, sNucl) - strNucl;
		//retValue = scoreMatrix[i * alphabetSize + j];
		retValue = scoreMatrix[i][j];
	}
	return retValue;
}

/* compute elements of matrix SW alg. */
void computeMatrixNW(char *first, char *second, element *mat, int m, int n, float gap, float match, float mismatch, int cols,
	char fastScore) {
	//											float *scoreMatrix, int alphabetSize, float ambCharPenalty, 
//                   element *bestScore, int *bestScore_i, int *bestScore_j, int cols) {
  int i, j, k;
  float match_mismatch;
  float scores[3];
  int prev;

  /* set upper-left corner element */
  mat[0].score = 0;
  mat[0].prev = START; // no prev. elem;

  /* set scores of 0th row and 0th column to -gap*j and -gap*i and their previous elements to left (L) or up (U), respectively */  
  for (j = 1; j <= n; j++) { // 0-th row
    mat[j].score = j * gap;
    mat[j].prev = LEFT;
  } 
  for (i = 1; i <= m; i++) { // 0-th column
    mat[i * cols].score = i * gap;
    mat[i * cols].prev = UPPER;
  } 
  
  /* set all other elements of matrix:
   * mat(i, j) = max { 0, mat(i-1, j-1) + w(a_i, b_j), mat(i-1, j) + w(a_i, -), mat(i, j-1) + w(-, b_j)}
   */
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++) {
			if (fastScore) {
				if (first[i - 1] == AMBIGUOUS_CHAR_QUERY || first[i - 1] == AMBIGUOUS_CHAR_SUBJECT
					|| second[j - 1] == AMBIGUOUS_CHAR_SUBJECT || second[j - 1] == AMBIGUOUS_CHAR_QUERY) { // treat ambiguous nucleotides as mismatches
					match_mismatch = mismatch;
				}
				else {
					match_mismatch = (first[i - 1] == second[j - 1]) ? match : mismatch; /* i-th position in matrix corresponds to (i-1)st char in first */
				}
			} // end fastScore
			else { // compute match/mismatch score based on score matrix 
				match_mismatch = getScore(first[i - 1], second[j - 1]);
			}
      scores[DIAGONAL] = mat[(i - 1) * cols + (j - 1)].score + match_mismatch;

      scores[UPPER] = mat[(i - 1) * cols + j].score + gap; // when backtracing: from mat(i, j) goes up to mat(i-1, j)
      scores[LEFT] = mat[i * cols + (j - 1)].score + gap; // when backtracing: from mat(i, j) goes left to mat(i, j-1)
      
      // find direction (diagonal, upper, left) that should be taken when backtracing from mat(i, j) (in order to get max score)
      prev = DIAGONAL; 
      for (k = 1; k < 3; k ++) { // find direction leading to max score
        if (scores[k] > scores[prev]) {
          prev = k;
        }
      } // end for
      mat[i * cols + j].score = scores[prev];
      mat[i * cols + j].prev = prev;
    }    
  } // end for  
}

/* tracing back the path from (m, n) until (0, 0) is encountered;
 * function computes aligned sequences firstNW and secondNW
*/
void tracebackBestScoreNW(char *first, char *second, char *firstNW, char *secondNW, element *mat, int m, int n, int cols) {
  int i = 0, j = 0;
  //int prev_i, prev_j;
  int curr_i, curr_j;
  
  /* starting from the bootom right element of the matrix */
  curr_i = m;
  curr_j = n;

  while (curr_i > 0 && curr_j > 0) { 
    //printf("%d %d\n", curr_i, curr_j);
    // find previous
    if (mat[curr_i * cols + curr_j].prev == DIAGONAL) {
      firstNW[i] = first[curr_i - 1];    
      secondNW[j] = second[curr_j - 1];      
      -- curr_j;
      -- curr_i;
    }
    else if (mat[curr_i * cols + curr_j].prev == UPPER) {
      firstNW[i] = first[curr_i - 1];    
      secondNW[j] = '-';      
      -- curr_i;
    }
    else if (mat[curr_i * cols + curr_j].prev == LEFT) {
      firstNW[i] = '-';    
      secondNW[j] = second[curr_j - 1];     
      -- curr_j;
    }
    ++ i;
    ++ j;
  } // end while

  while (curr_j > 0) { 
    //printf("%d %d\n", curr_i, curr_j);
    firstNW[i] = '-'; /* deletion */
    secondNW[j] = second[curr_j - 1];      
    -- curr_j;
    ++ i;
    ++ j;
  }

  while (curr_i > 0) { 
    //printf("%d %d\n", curr_i, curr_j);
    firstNW[i] = first[curr_i - 1]; /* insertion */
    secondNW[j] = '-';      
    -- curr_i;
    ++ i;
    ++ j;
  }

  firstNW[i] = '\0';
  secondNW[j] = '\0';
  reverse(firstNW);
  reverse(secondNW);
}

///* reverse a string, return str */ - already in stringUtil
//char *reverse(char *str) {
//
//  char *left  = str;
//  char *right = left + strlen(str) - 1;
//  char tmp;
//
//  while (left < right) {
//    tmp = *left;
//    *(left++) = *right;
//    *(right--) = tmp;
//  }
//  return str;
//}

/* compute elements of matrix NW algorithm */
//void computeMatrixNW(char *first, char *second, element *mat, int m, int n, double gap, double match, double mismatch, int cols) {
//  
//  int i, j;
//  double match_mismatch;
//  double diagonalScore, leftScore, upScore;
//  element e;
//  
//  /* set scores of 0th row and 0th column and their previous elements */
//  mat[0].score = 0.;
//  mat[0].prev = -1;
//
//  for (j = 1; j <= n; j++) {
//    mat[j].score = gap * j;
//    mat[j].i = 0;
//    mat[j].j = j - 1;
//  } 
//  for (i = 1; i <= m; i++) {
//    mat[i * cols].score = gap * i;
//    mat[i * cols].i = i - 1;
//    mat[i * cols].j = 0;
//  } 
//  
//  /* set all other elements of matrix:
//   * mat(i, j) = max { mat(i-1, j-1) + w(a_i, b_j), mat(i-1, j) + w(a_i, -), mat(i, j-1) + w(-, b_j)}
//   */
//  for (i = 1; i <= m; i++) {
//    for (j = 1; j <= n; j++) {
//      match_mismatch = (first[i - 1] == second[j - 1]) ? match : mismatch; /* i-th position in matrix corresponds to (i-1)st char in first */
//      diagonalScore = mat[(i - 1) * cols + (j - 1)].score + match_mismatch;
//      upScore = mat[(i - 1) * cols + j].score + gap; 
//      leftScore = mat[i * cols + (j - 1)].score + gap;
//      e = findMax(diagonalScore, i - 1, j - 1, upScore, i - 1 , j);
//      e = findMax(e.score, e.i, e.j, leftScore, i, j - 1);
//      mat[i * cols + j] = e;
//    }    
//  }
//}

/* tracing back the path from (m, n) until (0, 0) is encountered;
 * function computes aligned sequences firstNW and secondNW
*/
//void tracebackBestScoreNW(char *first, char *second, char *firstNW, char *secondNW, element *mat, int m, int n, int cols) {
//  
//  int i = 0, j = 0;
//  int prev_i, prev_j;
//  int curr_i, curr_j;
//  
//  /* starting from the bootom right element of the matrix */
//  curr_i = m;
//  curr_j = n;
//
//  while (curr_i > 0 && curr_j > 0) { 
//    //printf("%d %d\n", curr_i, curr_j);
//    prev_i = mat[curr_i * cols + curr_j].i;
//    prev_j = mat[curr_i * cols + curr_j].j;
//
//  if (curr_i == prev_i) {
//      firstNW[i] = '-'; /* deletion */
//      secondNW[j] = second[curr_j - 1];      
//    }
//    else if (curr_j == prev_j ) {
//      firstNW[i] = first[curr_i - 1]; /* insertion */
//      secondNW[j] = '-';      
//    }
//    else {
//      firstNW[i] = first[curr_i - 1]; /* i-th position in matrix corresponds to (i-1)st char in string first */
//      secondNW[j] = second[curr_j - 1];
//    }
//    ++ i;
//    ++ j;
//    curr_i = prev_i;
//    curr_j = prev_j;
//  } // end while
//
//  while (curr_j > 0) { 
//    //printf("%d %d\n", curr_i, curr_j);
//    firstNW[i] = '-'; /* deletion */
//    secondNW[j] = second[curr_j - 1];      
//    -- curr_j;
//    ++ i;
//    ++ j;
//  }
//
//  while (curr_i > 0) { 
//    //printf("%d %d\n", curr_i, curr_j);
//    firstNW[i] = first[curr_i - 1]; /* insertion */
//    secondNW[j] = '-';      
//    -- curr_i;
//    ++ i;
//    ++ j;
//  }
//
//  firstNW[i] = '\0';
//  secondNW[j] = '\0';
//  reverse(firstNW);
//  reverse(secondNW);
//}

/* set all score-s in matrix to 0 */
void setMatElemsToZero (element *mat, int maxCols) {

  int i, j;
  for (i = 0; i < maxCols; i++) {
    for (j = 0; j < maxCols; j++) {
      mat[i * maxCols + j].score = 0.;
    }
  }
}
