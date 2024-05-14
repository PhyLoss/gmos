/***** nw.h *************************************************************************************************
 * Description: Header file for local sequence alignment functions - using Needleman–WUnsch algorithm.
 * Author: Mirjana Domazet-Loso
 * 
 * This file is part of gmos.
 *
 ************************************************************************************************************/

#ifndef NW_H
#define NW_H

#define START -1
#define DIAGONAL 0
#define LEFT 1
#define UPPER 2

typedef struct elem {
  //int i;                         /* previous element (when backtracking) - row index */
  //int j;                         /* previous element (when backtracking) - column index */
  char prev;                    /* previous element: D(iagonal), U(p), L(eft) or 0 (starting position, i.e. upper-left corner) */
  float score;                  /* current score at this position in matrix */
} element;

extern double lambda; // query-specific 
extern double blastK; // query-specific 
extern float ambCharPenalty;
extern float scoreMatrix[ALPHABET_DNA_SIZE][ALPHABET_DNA_SIZE];

element *getMatrix(int m, int n);
void freeMatrix(element *mat);
/* set all score-s in matrix to 0 */
void setMatElemsToZero (struct elem *mat, int maxCols);

/* reverse a string, return str; already defined in stringUtil.h */
//char *reverse(char *str);

//void computeMatrixNW(char *first, char *second, element *mat, int m, int n, float gap, float match, float mismatch, int cols);
//void computeMatrixNW(char *first, char *second, element *mat, int m, int n, float gap, float match, float mismatch, int cols,
//	float *scoreMatrix, int alphabetSize, float ambCharPenalty);
void computeMatrixNW(char *first, char *second, element *mat, int m, int n, float gap, float match, float mismatch, 
	int cols, char fastScore);

void tracebackBestScoreNW(char *first, char *second, char *firstNW, char *secondNW, element *mat, int m, int n, int cols);

//float getScore(float *scoreMatrix, int alphabetSize, char qNucl, char sNucl, float ambCharPenalty);
float getScore(char qNucl, char sNucl);

#endif // NW_H
