/* Adopted from http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/algo/blast/core/blast_stat.c */

/* $Id: blast_stat.c 66158 2015-02-06 15:18:12Z boratyng $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
* Author: Tom Madden
*
*/

/** @file blast_stat.c
* Functions to calculate BLAST probabilities etc.
* Detailed Contents:
*
* - allocate and deallocate structures used by BLAST to calculate
* probabilities etc.
*
* - calculate residue frequencies for query and "average" database.
*
* - read in matrix or load it from memory.
*
*  - calculate sum-p from a collection of HSP's, for both the case
*   of a "small" gap and a "large" gap, when give a total score and the
*   number of HSP's.
*
* - calculate expect values for p-values.
*
* - calculate pseuod-scores from p-values.
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "blast_stat.h"

/** Check that the lo and hi score are within the allowed ranges
* @param lo the lowest permitted value [in]
* @param hi the highest permitted value [in]
* @return zero on success, 1 otherwise
*/

 Int2 BlastScoreChk(Int4 lo, Int4 hi) {
	if (lo >= 0 || hi <= 0 ||
		lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
		return 1;

	if (hi - lo > BLAST_SCORE_RANGE_MAX)
		return 1;

	return 0;
}

/* initialize Blast_ScoreFreq */
Blast_ScoreFreq* Blast_ScoreFreqNew(Int4 score_min, Int4 score_max) {
	Blast_ScoreFreq*  sfp;
	Int4  range;

	if (BlastScoreChk(score_min, score_max) != 0)
		return NULL;

	sfp = (Blast_ScoreFreq*)calloc(1, sizeof(Blast_ScoreFreq));
	if (sfp == NULL)
		return NULL;

	range = score_max - score_min + 1;
	sfp->sprob = (double*)calloc(range, sizeof(double));
	if (sfp->sprob == NULL)
	{
		Blast_ScoreFreqFree(sfp);
		return NULL;
	}

	sfp->sprob0 = sfp->sprob;
	sfp->sprob -= score_min;        /* center around 0 */
	sfp->score_min = score_min;
	sfp->score_max = score_max;
	sfp->obs_min = sfp->obs_max = 0;
	sfp->score_avg = 0.0;
	return sfp;
}

/** Calculates the score frequencies.
*
* @param sbp object with scoring information [in]
* @param sfp object to hold frequency information [in|out]
* @param rfp1 letter frequencies for first sequence (query) [in]
* @param rfp2 letter frequencies for second sequence (database) [in]
* @return zero on success
*/
//static Int2 BlastScoreFreqCalc(const BlastScoreBlk* sbp, Blast_ScoreFreq* sfp, Blast_ResFreq* rfp1, Blast_ResFreq* rfp2) {
Int2 BlastScoreFreqCalc(float *scoreMatrix, int	alphabetSize, Blast_ScoreFreq* sfp, float *expectedFreq) {
//	Int4 **matrix;
	Int4 score, obs_min, obs_max;
	double score_sum, score_avg;
	Int2 index1, index2;
	double lowScore = scoreMatrix[0], highScore = scoreMatrix[0];
	for (index1 = 0; index1 < alphabetSize; index1++) {
		for (index2 = 0; index2 < alphabetSize; index2++) {
			if (scoreMatrix[index1 * alphabetSize + index2] < lowScore) {
				lowScore = scoreMatrix[index1 * alphabetSize + index2];
			}
			if (scoreMatrix[index1 * alphabetSize + index2] > highScore) {
				highScore = scoreMatrix[index1 * alphabetSize + index2];
			}
		}
	}
	
	if (scoreMatrix == NULL || sfp == NULL)
		return 1;

	if (lowScore < sfp->score_min || highScore > sfp->score_max)
		return 1;

	for (score = sfp->score_min; score <= sfp->score_max; score++)
		sfp->sprob[score] = 0.0;

	//matrix = sbp->matrix->data;
	for (index1 = 0; index1 < alphabetSize; index1 ++) {
		for (index2 = 0; index2 < alphabetSize; index2 ++) {
			score = (Int4)scoreMatrix[index1 * alphabetSize + index2];
			if (score >= lowScore) {
				if (index1 == index2)
					sfp->sprob[score] += expectedFreq[index1 * alphabetSize + index2];
				else 
					sfp->sprob[score] += expectedFreq[index1 * alphabetSize + index2] / 2;
			}
		}
	}

	score_sum = 0.;
	obs_min = obs_max = BLAST_SCORE_MIN;
	for (score = sfp->score_min; score <= sfp->score_max; score++) {
		if (sfp->sprob[score] > 0.) {
			score_sum += sfp->sprob[score];
			obs_max = score;
			if (obs_min == BLAST_SCORE_MIN)
				obs_min = score;
		}
	}
	sfp->obs_min = obs_min;
	sfp->obs_max = obs_max;

	score_avg = 0.0;
	if (score_sum > 0.0001 || score_sum < -0.0001) {
		for (score = obs_min; score <= obs_max; score++) {
			sfp->sprob[score] /= score_sum;
			score_avg += score * sfp->sprob[score];
		}
	}
	sfp->score_avg = score_avg;
	return 0;
}

/** Calculate H, the relative entropy of the p's and q's
*
* @param sfp object containing scoring frequency information [in]
* @param lambda a Karlin-Altschul parameter [in]
* @return H, a Karlin-Altschul parameter
*/
double BlastKarlinLtoH(Blast_ScoreFreq* sfp, double lambda) {
	Int4 score;
	double H, etonlam, sum, scale;

	double *probs = sfp->sprob;
	Int4 low = sfp->obs_min, high = sfp->obs_max;

	if (lambda < 0.) {
		return -1.;
	}
	if (BlastScoreChk(low, high) != 0) return -1.;

	etonlam = exp(-lambda);
	sum = low * probs[low];
	for (score = low + 1; score <= high; score++) {
		sum = score * probs[score] + etonlam * sum;
	}

	scale = BLAST_Powi(etonlam, high);
	if (scale > 0.0) {
		H = lambda * sum / scale;
	}
	else { /* Underflow of exp( -lambda * high ) */
		H = lambda * exp(lambda * high + log(sum));
	}
	return H;
}

double BLAST_Powi(double x, Int4 n) {
	double y;

	if (n == 0)
		return 1.;

	if (x == 0.) {
		if (n < 0) {
			return DBL_MAX;
		}
		return 0.;
	}

	if (n < 0) {
		x = 1. / x;
		n = -n;
	}

	y = 1.;
	while (n > 0) {
		if (n & 1)
			y *= x;
		n /= 2;
		x *= x;
	}
	return y;
}

/**************** Statistical Significance Parameter Subroutine ****************

Version 1.0     February 2, 1990
Version 1.2     July 6,     1990

Program by:     Stephen Altschul

Address:        National Center for Biotechnology Information
National Library of Medicine
National Institutes of Health
Bethesda, MD  20894

Internet:       altschul@ncbi.nlm.nih.gov

See:  Karlin, S. & Altschul, S.F. "Methods for Assessing the Statistical
Significance of Molecular Sequence Features by Using General Scoring
Schemes,"  Proc. Natl. Acad. Sci. USA 87 (1990), 2264-2268.

Computes the parameters lambda and K for use in calculating the
statistical significance of high-scoring segments or subalignments.

The scoring scheme must be integer valued.  A positive score must be
possible, but the expected (mean) score must be negative.

A program that calls this routine must provide the value of the lowest
possible score, the value of the greatest possible score, and a pointer
to an array of probabilities for the occurrence of all scores between
these two extreme scores.  For example, if score -2 occurs with
probability 0.7, score 0 occurs with probability 0.1, and score 3
occurs with probability 0.2, then the subroutine must be called with
low = -2, high = 3, and pr pointing to the array of values
{ 0.7, 0.0, 0.1, 0.0, 0.0, 0.2 }.  The calling program must also provide
pointers to lambda and K; the subroutine will then calculate the values
of these two parameters.  In this example, lambda=0.330 and K=0.154.

The parameters lambda and K can be used as follows.  Suppose we are
given a length N random sequence of independent letters.  Associated
with each letter is a score, and the probabilities of the letters
determine the probability for each score.  Let S be the aggregate score
of the highest scoring contiguous segment of this sequence.  Then if N
is sufficiently large (greater than 100), the following bound on the
probability that S is greater than or equal to x applies:

P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ].

In other words, the p-value for this segment can be written as
1-exp[-KN*exp(-lambda*S)].

This formula can be applied to pairwise sequence comparison by assigning
scores to pairs of letters (e.g. amino acids), and by replacing N in the
formula with N*M, where N and M are the lengths of the two sequences
being compared.

In addition, letting y = KN*exp(-lambda*S), the p-value for finding m
distinct segments all with score >= S is given by:

2             m-1           -y
1 - [ 1 + y + y /2! + ... + y   /(m-1)! ] e

Notice that for m=1 this formula reduces to 1-exp(-y), which is the same
as the previous formula.

*******************************************************************************/

/** The following procedure computes K. The input includes Lambda, H,
*  and an array of probabilities for each score.
*  There are distinct closed form for three cases:
*  1. high score is 1 low score is -1
*  2. high score is 1 low score is not -1
*  3. low score is -1, high score is not 1
*
* Otherwise, in most cases the value is computed as:
* -exp(-2.0*outerSum) / ((H/lambda)*(exp(-lambda) - 1)
* The last term (exp(-lambda) - 1) can be computed in two different
* ways depending on whether lambda is small or not.
* outerSum is a sum of the terms
* innerSum/j, where j is denoted by iterCounter in the code.
* The sum is truncated when the new term innersum/j i sufficiently small.
* innerSum is a weighted sum of the probabilities of
* of achieving a total score i in a gapless alignment,
* which we denote by P(i,j).
* of exactly j characters. innerSum(j) has two parts
* Sum over i < 0  P(i,j)exp(-i * lambda) +
* Sum over i >=0  P(i,j)
* The terms P(i,j) are computed by dynamic programming.
* An earlier version was flawed in that ignored the special case 1
* and tried to replace the tail of the computation of outerSum
* by a geometric series, but the base of the geometric series
* was not accurately estimated in some cases.
*
* @param sfp object holding scoring frequency information [in]
* @param lambda a Karlin-Altschul parameter [in]
* @param H a Karlin-Altschul parameter [in]
* @return K, another Karlin-Altschul parameter
*/

double BlastKarlinLHtoK(Blast_ScoreFreq* sfp, double lambda, double H)
{
	/*The next array stores the probabilities of getting each possible
	score in an alignment of fixed length; the array is shifted
	during part of the computation, so that
	entry 0 is for score 0.  */
	double         *alignmentScoreProbabilities = NULL;
	Int4            low;    /* Lowest score (must be negative) */
	Int4            high;   /* Highest score (must be positive) */
	Int4            range;  /* range of scores, computed as high - low*/
	double          K;      /* local copy of K  to return*/
	int             i;   /*loop index*/
	int             iterCounter; /*counter on iterations*/
	Int4            divisor; /*candidate divisor of all scores with
													 non-zero probabilities*/
	/*highest and lowest possible alignment scores for current length*/
	Int4            lowAlignmentScore, highAlignmentScore;
	Int4            first, last; /*loop indices for dynamic program*/
	register double innerSum;
	double          oldsum, oldsum2;  /* values of innerSum on previous
																		iterations*/
	double          outerSum;        /* holds sum over j of (innerSum
																	 for iteration j/j)*/

	double          score_avg; /*average score*/
	/*first term to use in the closed form for the case where
	high == 1 or low == -1, but not both*/
	double          firstTermClosedForm;  /*usually store H/lambda*/
	int             iterlimit; /*upper limit on iterations*/
	double          sumlimit; /*lower limit on contributions
														to sum over scores*/

	/*array of score probabilities reindexed so that low is at index 0*/
	double         *probArrayStartLow;

	/*pointers used in dynamic program*/
	double         *ptrP, *ptr1, *ptr2, *ptr1e;
	double          expMinusLambda; /*e^^(-Lambda) */

	if (lambda <= 0. || H <= 0.) {
		/* Theory dictates that H and lambda must be positive, so
		* return -1 to indicate an error */
		return -1.;
	}

	/*Karlin-Altschul theory works only if the expected score
	is negative*/
	if (sfp->score_avg >= 0.0) {
		return -1.;
	}

	low = sfp->obs_min;
	high = sfp->obs_max;
	range = high - low;

	probArrayStartLow = &sfp->sprob[low];
	/* Look for the greatest common divisor ("delta" in Appendix of PNAS 87 of
	Karlin&Altschul (1990) */
	for (i = 1, divisor = -low; i <= range && divisor > 1; ++i) {
		if (probArrayStartLow[i] != 0.0)
			divisor = BLAST_Gcd(divisor, i);
	}

	high /= divisor;
	low /= divisor;
	lambda *= divisor;

	range = high - low;

	firstTermClosedForm = H / lambda;
	expMinusLambda = exp((double)-lambda);

	if (low == -1 && high == 1) {
		K = (sfp->sprob[low*divisor] - sfp->sprob[high*divisor]) *
			(sfp->sprob[low*divisor] - sfp->sprob[high*divisor]) / sfp->sprob[low*divisor];
		return(K);
	}

	if (low == -1 || high == 1) {
		if (high != 1) {
			score_avg = sfp->score_avg / divisor;
			firstTermClosedForm
				= (score_avg * score_avg) / firstTermClosedForm;
		}
		return firstTermClosedForm * (1.0 - expMinusLambda);
	}

	sumlimit = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
	iterlimit = BLAST_KARLIN_K_ITER_MAX;

	alignmentScoreProbabilities =
		(double *)calloc((iterlimit*range + 1), sizeof(*alignmentScoreProbabilities));
	if (alignmentScoreProbabilities == NULL)
		return -1.;

	outerSum = 0.;
	lowAlignmentScore = highAlignmentScore = 0;
	alignmentScoreProbabilities[0] = innerSum = oldsum = oldsum2 = 1.;

	for (iterCounter = 0;
		((iterCounter < iterlimit) && (innerSum > sumlimit));
		outerSum += innerSum /= ++iterCounter) {
		first = last = range;
		lowAlignmentScore += low;
		highAlignmentScore += high;
		/*dynamic program to compute P(i,j)*/
		for (ptrP = alignmentScoreProbabilities +
			(highAlignmentScore - lowAlignmentScore);
			ptrP >= alignmentScoreProbabilities;
			*ptrP-- = innerSum) {
			ptr1 = ptrP - first;
			ptr1e = ptrP - last;
			ptr2 = probArrayStartLow + first;
			for (innerSum = 0.; ptr1 >= ptr1e;) {
				innerSum += *ptr1  *  *ptr2;
				ptr1--;
				ptr2++;
			}
			if (first)
				--first;
			if (ptrP - alignmentScoreProbabilities <= range)
				--last;
		}
		/* Horner's rule */
		innerSum = *++ptrP;
		for (i = lowAlignmentScore + 1; i < 0; i++) {
			innerSum = *++ptrP + innerSum * expMinusLambda;
		}
		innerSum *= expMinusLambda;

		for (; i <= highAlignmentScore; ++i)
			innerSum += *++ptrP;
		oldsum2 = oldsum;
		oldsum = innerSum;
	}

#ifdef ADD_GEOMETRIC_TERMS_TO_K
	/*old code assumed that the later terms in sum were
	asymptotically comparable to those of a geometric
	progression, and tried to speed up convergence by
	guessing the estimated ratio between sucessive terms
	and using the explicit terms of a geometric progression
	to speed up convergence. However, the assumption does not
	always hold, and convergenece of the above code is fast
	enough in practice*/
	/* Terms of geometric progression added for correction */
		{
			double     ratio;  /* fraction used to generate the
												 geometric progression */

			ratio = oldsum / oldsum2;
			if (ratio >= (1.0 - sumlimit*0.001)) {
				K = -1.;
				if (alignmentScoreProbabilities != NULL)
					sfree(alignmentScoreProbabilities);
				return K;
			}
			sumlimit *= 0.01;
			while (innerSum > sumlimit) {
				oldsum   *= ratio;
				outerSum += innerSum = oldsum / ++iterCounter;
	}
}
#endif

	K = -exp((double)-2.0*outerSum) /
		(firstTermClosedForm*BLAST_Expm1(-(double)lambda));

	if (alignmentScoreProbabilities != NULL)
		sfree(alignmentScoreProbabilities);

	return K;
}



Int4 BLAST_Gcd(Int4 a, Int4 b) {
	Int4   c;

	b = ABS(b);
	if (b > a)
		c = a, a = b, b = c;

	while (b != 0) {
		c = a%b;
		a = b;
		b = c;
	}
	return a;
}


Blast_ScoreFreq* Blast_ScoreFreqFree(Blast_ScoreFreq* sfp) {
	if (sfp == NULL)
		return NULL;

	if (sfp->sprob0 != NULL)
		sfree(sfp->sprob0);
	sfree(sfp);
	return sfp;
}

double BLAST_Expm1(double x) {
	double absx = ABS(x);

	if (absx > .33)
		return exp(x) - 1.;

	if (absx < 1.e-16)
		return x;

	return x * (1. + x *
		(1. / 2. + x *
		(1. / 6. + x *
		(1. / 24. + x *
		(1. / 120. + x *
		(1. / 720. + x *
		(1. / 5040. + x *
		(1. / 40320. + x *
		(1. / 362880. + x *
		(1. / 3628800. + x *
		(1. / 39916800. + x *
		(1. / 479001600. +
		x / 6227020800.))))))))))));
}

void sfree(void *x) {
	free(x);
	x = NULL;
	return;
}

/* calculate K */
double getK(float *scoreMatrix, float *expectedFreq, double lambda, int alphabetSize) {
	Int4 score_min = BLAST_SCORE_MIN;
	Int4 score_max = BLAST_SCORE_MAX;
	double H, K;

	Blast_ScoreFreq* sfp = Blast_ScoreFreqNew(score_min, score_max);

	if (BlastScoreFreqCalc(scoreMatrix, alphabetSize, sfp, expectedFreq)) {
		printf("[ERROR] Karlin-Altschul parameter K cannot be calculated.\n");
	}
	H = BlastKarlinLtoH(sfp, lambda);
	K = BlastKarlinLHtoK(sfp, lambda, H);
	sfp = Blast_ScoreFreqFree(sfp);
	return K;
}


