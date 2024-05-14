/***** sw.h *************************************************************************************************
 * Description: Header file for local sequence alignment functions - using Smith–Waterman algorithm.
 * Author: Mirjana Domazet-Loso, March 19, 2012
 *
 ************************************************************************************************************/

#ifndef SW_H
#define SW_H

typedef struct elem {
  int i;                         /* previous element (when backtracking) - row index */
  int j;                         /* previous element (when backtracking) - column index */
  double score;                  /* current score at this position in matrix */
} element;

element *getMatrix(int m, int n);
void freeMatrix(element *mat);

element findMax(double score1, int i1, int j1, double score2, int i2, int j2);

/* reverse a string, return str; already defined in stringUtil.h */
//char *reverse(char *str);

/* set all score-s in matrix to 0 */
void setMatElemsToZero (struct elem *mat, int maxCols);

void computeMatrixNW(char *first, char *second, element *mat, int m, int n, double gap, double match, double mismatch, int cols);
void tracebackBestScoreNW(char *first, char *second, char *firstSW, char *secondSW, element *mat, int m, int n, int cols);

#endif // SW_H
