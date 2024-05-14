/***** blast_stat.h *************************************************************
* Description: Header file for blast functions.
* Author: Mirjana Domazet-Loso
* Adapted from http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/
*****************************************************************************/
#ifndef BLAST_STAT_H
#define BLAST_STAT_H

typedef short Int2;
typedef long  Int4;
#ifndef INT2_MAX
/** largest number represented by signed (two byte) short */
#define INT2_MAX 32767
#endif

#ifndef INT2_MIN
/** smallest (most negative) number represented by signed (two byte) short */
#define INT2_MIN (-32768)
#endif

#define BLAST_SCORE_MIN INT2_MIN   /**< minimum allowed score (for one letter comparison). */
#define BLAST_SCORE_MAX INT2_MAX   /**< maximum allowed score (for one letter comparison). */
#define BLAST_SCORE_RANGE_MAX   (BLAST_SCORE_MAX - BLAST_SCORE_MIN) /**< maximum allowed range of BLAST scores. */

#define BLAST_KARLIN_K_SUMLIMIT_DEFAULT 0.0001 /**< K_SUMLIMIT_DEFAULT == sumlimit used in BlastKarlinLHtoK() */
#define BLAST_KARLIN_K_ITER_MAX 100 /**< upper limit on iterations for BlastKarlinLHtoK */
//#define ADD_GEOMETRIC_TERMS_TO_K
/** Holds score frequencies used in calculation of Karlin-Altschul parameters for an ungapped search. */
typedef struct Blast_ScoreFreq {
	Int4         score_min; /**< lowest allowed scores */
	Int4         score_max; /**< highest allowed scores */
	Int4         obs_min;   /**< lowest observed (actual) scores */
	Int4         obs_max;   /**< highest observed (actual) scores */
	double       score_avg; /**< average score, must be negative for local alignment. */
	double*      sprob0;    /**< arrays for frequency of given score */
	double*      sprob;     /**< arrays for frequency of given score, shifted down by score_min. */
} Blast_ScoreFreq;

#ifndef ABS
/** returns absolute value of a (|a|) */
#define ABS(a)  ((a)>=0?(a):-(a))
#endif

/********* functions *******************/
Int2 BlastScoreChk(Int4 lo, Int4 hi);

Blast_ScoreFreq* Blast_ScoreFreqNew(Int4 score_min, Int4 score_max);
Int2 BlastScoreFreqCalc(float *scoreMatrix, int	alphabetSize, Blast_ScoreFreq* sfp, float *expectedFreq);

double BLAST_Powi(double x, Int4 n);
double BlastKarlinLtoH(Blast_ScoreFreq* sfp, double lambda);

double BlastKarlinLHtoK(Blast_ScoreFreq* sfp, double lambda, double H);

Int4 BLAST_Gcd(Int4 a, Int4 b);
Blast_ScoreFreq* Blast_ScoreFreqFree(Blast_ScoreFreq* sfp);

double BLAST_Expm1(double x);
void sfree(void *x);

/* calculate K */
double getK(float *scoreMatrix, float *expectedFreq, double lambda, int alphabetSize);

#endif // BLAST_STAT_H
