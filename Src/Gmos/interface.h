/***** interface.h **********************************************************
 * Description: Header file for user interface.
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
 *
 *****************************************************************************/ 
#ifndef INTERFACE_H
#define INTERFACE_H

#define DEFAULT_D 1000000
#define DEFAULT_m 1.f
#define DEFAULT_P 0.4
#define DEFAULT_f 200
#define DEFAULT_RATIO 0.3 // F = ratio * f
#define DEFAULT_SHORTER_LENGTH 0.25 // alllowed decrease in region's length compensated by 0.01 increase in region's ratio
#define DEFAULT_F 0	
#define DEFAULT_O 0.5
#define DEFAULT_T 0

#define DEFAULT_A 70
#define DEFAULT_B 0
#define DEFAULT_G -2.0f 
#define DEFAULT_M 1.0f
#define DEFAULT_S -1.0f
#define DEFAULT_K 0.1
#define DEFAULT_LAMBDA 1.58
#define ALPHABET_DNA_SIZE 4

#define EVALUE_LOW 1.e-100
#define RATIO_MINFRAGMENT 0.80

#define MAX_STRNUCL 4 // "ACGT"
//#define NAMELENGTH 50 // already defined in stringUtil.h

/* define argument container */
typedef struct args {
	char **i;               /* name of query files */
	int queryFileNumber;    /* number of input query files */

	int subjectFileNumber;	/* number of input subject files */
	char **j;				  /* list of names of subject files (1 or more); -j would include subject files from the working directory; 
							      * Option -j and option -d are not allowed both at the same time!*/	
  
	//char *d;				  /* name of the directory where the subject files are located; 
	//						      * Option -j and option -d are not allowed both at the same time!*/	

	//char *o;          /* name of output file - all results */
	char *O;          /* name of output file - all results */
	char *q;          /* name of output file - only final results in column mode, values separated by a delimiter */

	int h;          /* help message? */
	int v;          /* version? */
	int p;          /* print program information */	  
	int e;          /* error message? */
//	int r;          /* print reverse strand (default: only forward strand) */
	int t;          /* print run-time */	  
  //char *Q;        /* compute and print final list (query as a mosaic structure) */
  char *o;        /* compute and print final list (query as a mosaic structure) */
  char *a;        /* print alignments (multi-fasta format) */

	//int n;          /* number of characters printed from header line */
	int D;          /* maximum depth of suffix tree */

	double P;       /* fraction of random shustrings excluded from analysis; P e [0, 1] */
	double m;       /* multiplier for maxshulen - defines the threshold of a significant interval: threshold = m * maxshulen_by_chance-alone */
  	
	int f;			    /* minimal length of fragment */
  double F;			  /* minimal score of fragment */
  double R;			  /* minimal ratio = score / length of fragment*/

  //double O;       /* precentage of overlapping area between two regions - to be displayed in final list; [0, 1] */
  Int64 T;      /* threshold for a strong signal; only those intervals starting with sl > T are considered */
  Int64 L;      /* lcp threshold - function process2 is called only when lcp of an interval > L; default = max shulen by chance alone */    

	int s;			/* include only strong signal (default: not/0) */

  short A;      /* maximal distance between two segments to be aligned, e.g. 20 */
  short B;      /* maximal distance of the aligned segment just before and just after a region, e.g. 10 */
  int C;        /* maximal allowed distance between two neighboring intervals in the connected final merged list */

  float G;      /* initial gap penalty, e.g. -10 */
  float M;      /* initial match, e.g. +1 */
  float S;      /* initial mismatch penalty, e.g. -1 */
  
  double E;     /* E-value of an aligned segment to be considered relevant */
  double r;     /* proportion of the minimal fragment length to be kept in case of overlapping */
	double X;     /* "excellent" E-value of an aligned segment; if two or more overlapping segments all have excellent
								*  e-values, then they are all kept */
} Args;

const char strNucl[MAX_STRNUCL + 1]; 
Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
void printSplash(char *version);

void freeArgs(Args *args);

// check multiple definition of an argument
void checkMultipleArgDef(char list[], char c);

/* check whether an argument is a null pointer */
void checkNullPointer(char c, char *arg);

/* get all file names */
void getFileNames(int argc, char *argv[], char list[], char c, char ***fileNames, int *cntFiles, int *arg);

#endif // INTERFACE_H
