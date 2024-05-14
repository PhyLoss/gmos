/***** sw.c *************************************************************************************************
 * Description: Functions for sequence alignment using Smith-Waterman algorithm.
 * Author: Mirjana Domazet-Loso, March 19, 2012
 *
 * This file is part of qa.
 *
 ************************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "commonSC.h"
#include "eprintf.h"
#include "stringUtil.h"

#include "sw.h"

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
  
  element *mat;

  mat = (element *)/*e*/malloc(sizeof(element) * rows * cols);
  return mat;

}

void freeMatrix(element *mat) {
  free(mat);
}

element findMax(double score1, int i1, int j1, double score2, int i2, int j2) {
  
  element e;
  
  if (score1 >= score2) {
    e.score = score1;
    e.i = i1;
    e.j = j1;
  }
  else {
    e.score = score2;
    e.i = i2;
    e.j = j2;  
  }
  return e;
}

/* compute elements of matrix SW alg. */
void computeMatrix(char *first, char *second, element *mat, int m, int n, double gap, double match, double mismatch, 
                   element *bestScore, int *bestScore_i, int *bestScore_j, int cols) {
  
  int i, j;
  double match_mismatch;
  double diagonalScore, leftScore, upScore;
  element e;
  
  /* set scores of 0th row and 0th column to 0 and their previous elements to (-1, -1) */
  for (j = 0; j <= n; j++) {
    mat[j].score = 0.;
    mat[j].i = mat[j].j = -1;
  } 
  for (i = 1; i <= m; i++) {
    mat[i * cols].score = 0.;
    mat[i * cols].i = mat[i * cols].j = -1;
  } 
  
  /* initialize best score element to (0, 0) */
  (*bestScore).score = (*bestScore).i = (*bestScore).j = 0;

  /* set all other elements of matrix:
   * mat(i, j) = max { 0, mat(i-1, j-1) + w(a_i, b_j), mat(i-1, j) + w(a_i, -), mat(i, j-1) + w(-, b_j)}
   */
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++) {
      match_mismatch = (first[i - 1] == second[j - 1]) ? match : mismatch; /* i-th position in matrix corresponds to (i-1)st char in first */
      diagonalScore = mat[(i - 1) * cols + (j - 1)].score + match_mismatch;
      upScore = mat[(i - 1) * cols + j].score + gap; 
      leftScore = mat[i * cols + (j - 1)].score + gap;
      e = findMax(diagonalScore, i - 1, j - 1, upScore, i - 1 , j);
      e = findMax(e.score, e.i, e.j, leftScore, i, j - 1);
      if (e.score < 0) { /* every score in mat has to be non-negative */      
        e.score = 0;
        e.i = e.j = -1; /* pointing nowhere, backward search stops here */
      }
      mat[i * cols + j] = e;
      if ((*bestScore).score < e.score) { /* new element with best score */
        *bestScore_i = i;
        *bestScore_j = j;
        *bestScore = e;
      }
    }    
  }
}

/* tracing back the path od local alignment starting from (bestScore.i, bestScore.j) until first 0 is encountered;
 * function computes aligned sequences firstNW and secondNW. BUG!!!!!! CHECH THIS LATER!!!
*/
void tracebackBestScore(char *first, char *second, element bestScore, int bestScore_i, int bestScore_j, 
                        char *firstNW, char *secondNW, element *mat, int cols) {
  
  int i = 0, j = 0;
  element e;
  //int prev_i, prev_j;
  int prev_i = bestScore_i;
  int prev_j = bestScore_j;
  
  /* starting from element (bestScore_i, bestScore_j) */
  firstNW[i++] = first[bestScore_i - 1]; /* i-th position in matrix corresponds to (i-1)st char in first */
  secondNW[j++] = second[bestScore_j - 1];
  e = bestScore;

  while (mat[e.i * cols + e.j].score) { /* if predecessor's score is > 0*/
  //while (e.score) {
    e = mat[e.i * cols + e.j];
    if (e.i == prev_i || e.i == 0) {
      firstNW[i] = '-'; /* deletion */
      secondNW[j] = second[e.j - 1];      
    }
    //else if (e.j == prev_i || e.j == 0) {
    else if (e.j == prev_j || e.j == 0) {
      firstNW[i] = first[e.i - 1]; /* insertion */
      secondNW[j] = '-';      
    }
    else {
      firstNW[i] = first[e.i - 1]; /* i-th position in matrix corresponds to (i-1)st char in string first */
      secondNW[j] = second[e.j - 1];
    }
    ++ i;
    ++ j;
    prev_i = e.i;
    prev_j = e.j;
    //e = mat[e.i * cols + e.j]; /* trace back to predecessor of e */  
  } // end while

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
void computeMatrixNW(char *first, char *second, element *mat, int m, int n, double gap, double match, double mismatch, int cols) {
  
  int i, j;
  double match_mismatch;
  double diagonalScore, leftScore, upScore;
  element e;
  
  /* set scores of 0th row and 0th column and their previous elements */
  mat[0].score = 0.;
  mat[0].i = mat[0].j = -1;

  for (j = 1; j <= n; j++) {

    mat[j].score = gap * j;
    mat[j].i = 0;
    mat[j].j = j - 1;
  } 
  for (i = 1; i <= m; i++) {
    mat[i * cols].score = gap * i;
    mat[i * cols].i = i - 1;
    mat[i * cols].j = 0;
  } 
  
  /* set all other elements of matrix:
   * mat(i, j) = max { mat(i-1, j-1) + w(a_i, b_j), mat(i-1, j) + w(a_i, -), mat(i, j-1) + w(-, b_j)}
   */
  for (i = 1; i <= m; i++) {
    for (j = 1; j <= n; j++) {
      match_mismatch = (first[i - 1] == second[j - 1]) ? match : mismatch; /* i-th position in matrix corresponds to (i-1)st char in first */
      diagonalScore = mat[(i - 1) * cols + (j - 1)].score + match_mismatch;
      upScore = mat[(i - 1) * cols + j].score + gap; 
      leftScore = mat[i * cols + (j - 1)].score + gap;
      e = findMax(diagonalScore, i - 1, j - 1, upScore, i - 1 , j);
      e = findMax(e.score, e.i, e.j, leftScore, i, j - 1);
      mat[i * cols + j] = e;
    }    
  }
}

/* tracing back the path from (m, n) until (0, 0) is encountered;
 * function computes aligned sequences firstNW and secondNW
*/
void tracebackBestScoreNW(char *first, char *second, char *firstNW, char *secondNW, element *mat, int m, int n, int cols) {
  
  int i = 0, j = 0;
  int prev_i, prev_j;
  int curr_i, curr_j;
  
  /* starting from the bootom right element of the matrix */
  curr_i = m;
  curr_j = n;

  while (curr_i > 0 && curr_j > 0) { 
    //printf("%d %d\n", curr_i, curr_j);
    prev_i = mat[curr_i * cols + curr_j].i;
    prev_j = mat[curr_i * cols + curr_j].j;

  if (curr_i == prev_i) {
      firstNW[i] = '-'; /* deletion */
      secondNW[j] = second[curr_j - 1];      
    }
    else if (curr_j == prev_j ) {
      firstNW[i] = first[curr_i - 1]; /* insertion */
      secondNW[j] = '-';      
    }
    else {
      firstNW[i] = first[curr_i - 1]; /* i-th position in matrix corresponds to (i-1)st char in string first */
      secondNW[j] = second[curr_j - 1];
    }
    ++ i;
    ++ j;
    curr_i = prev_i;
    curr_j = prev_j;
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

/* set all score-s in matrix to 0 */
void setMatElemsToZero (element *mat, int maxCols) {

  int i, j;
  for (i = 0; i < maxCols; i++) {
    for (j = 0; j < maxCols; j++) {
      mat[i * maxCols + j].score = 0.;
    }
  }
}
