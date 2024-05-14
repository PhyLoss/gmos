/***** interface.c ************************************************
 * Description: Functions for gathering arguments from the command line.
 * Author: Mirjana Domazet-Loso
*
* This file is part of gmos.
*****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include "commonSC.h"
#include "interface.h"
#include "eprintf.h"
#include "stringUtil.h"

#if defined(_DEBUG) && defined(WIN) 
#include "leakWatcher.h"
#endif

#if defined(_DEBUG) && defined(WIN) 
  #define new DEBUG_NEW
  #undef THIS_FILE
  static char THIS_FILE[] = __FILE__;
#endif

const char strNucl[] = "ACGT";

// check multiple definition of an argument
// if the argument is already defined, then produce warning??
#if defined(WIN)
static 
#endif
void checkMultipleArgDef(char list[], char c) {

  if(!list[(int)c]){
    list[(int)c] = 1;
  }
  else {
    printf("Warning: -%c multiply defined\n", c); // or error ???????
  }
}

/* check whether an argument is a null pointer */
#if defined(WIN)
static 
#endif
void checkNullPointer(char c, char *arg) {	

  if (!arg) {
    eprintf("[ERROR] Option -%c: argument is NULL\n", c, arg);
  }
}

/* read one or more query or subject files specified with an option -i or -j (-d) */
#if defined(WIN)
static 
#endif
void getFileNames(int argc, char *argv[], char list[], char c, char ***fileNames, int *cntFiles, int *arg) {
  char **r;
  int numberOfFiles = 1;

  // get all file names (supposition: file name cannot start with -)
  while ((*arg + numberOfFiles) < argc && argv[*arg + numberOfFiles][0] != '-') {
    *fileNames = (char **)erealloc(*fileNames, numberOfFiles * sizeof (char *));
    r = *fileNames;
    r[numberOfFiles - 1] = argv[*arg + numberOfFiles];
    checkNullPointer(c, r[numberOfFiles - 1]);
    ++ numberOfFiles;
  }
  *cntFiles = numberOfFiles - 1; // set the argument
	
  // position of the new argument (-1 because arg is increased by 1 as the last step of while loop)
  *arg = *arg + numberOfFiles - 1;
}

Args *getArgs(int argc, char *argv[]){
  Args *args;
  char c;
	char flagExcellentValue = 0;
  int arg;
  char list[256] = {0}; // options are not yet set

  args = (Args *)/*e*/malloc(sizeof(Args));

  args->i = NULL;
  args->j = NULL;
//  args->d = NULL;
  //args->o = NULL;
  args->O = NULL;
  args->q = NULL;

  args->h = 0;
  args->v = 0;
  args->p = 0;  
  args->e = 0;
//  args->r = 0;
  //args->n = NAMELENGTH;
  args->D = DEFAULT_D;
  args->P = 1-DEFAULT_P;
  args->m = DEFAULT_m;
  args->t = 0;
  //args->Q = NULL;
  args->o = NULL;
  args->a = NULL;
  args->s = 1;
  args->T = DEFAULT_T; /* change to maxshulen_by_chance_alone + 1 ??? */
  args->L = -1; /* lcp threshold for calling "process2" */
  arg = 1;

  args->A = DEFAULT_A; 
  args->B = DEFAULT_B; 
  args->C = 0; 
  args->G = DEFAULT_G; 
  args->M = DEFAULT_M; 
  args->S = DEFAULT_S; 
  args->f = DEFAULT_f;
  args->F = DEFAULT_F;
  args->R = DEFAULT_RATIO;
	args->E = EVALUE_LOW;
  args->X = EVALUE_LOW;
  args->r = RATIO_MINFRAGMENT;

  while (arg < argc){
    c = argv[arg][1];
    switch (c){    
    checkMultipleArgDef(list, c);

    case 'i':                           /* query file */ 
      checkNullPointer(c, argv[arg + 1]);      
      // read one or more query files specified with an option -i
      getFileNames(argc, argv, list, c, &args->i, &args->queryFileNumber, &arg);
      break;
				      
    case 'j':                           /* subject files */
      /* check that -d and -j option are used mutually exclusive! */
      //if (args->d) {
		    //eprintf("Options -d and -j are mutually exclusive");
      //}
      //else {
        checkNullPointer(c, argv[arg + 1]);
		    getFileNames(argc, argv, list, c, &args->j, &args->subjectFileNumber, &arg);
		    break;
      //}

    //case 'd':                           /* directory where the subject files are saved and possibly the query file */
    //  if (args->j) {
		  //  eprintf("Options -d and -j are mutually exclusive");
    //  }
    //  else {
		  //  args->d = argv[++arg];
		  //  checkNullPointer(c, args->d);
		  //  break;
    //  }

    //case 'o':     /* output file */
    case 'O':       /* output file */
      checkNullPointer(c, argv[++arg]);
      args->O = argv[arg];
      break;

    case 'q':                           /* output file - final results separated by delimiter */
      checkNullPointer(c, argv[++arg]);
      args->q = argv[arg];
      break;

    //case 'n':                            /* number of characters from header line printed */
    //  checkNullPointer(c, argv[++arg]);
    //  args->n = atoi(argv[arg]);
    //  if (args->n < 0) {
		  //  args->n = INT_MAX;
    //  }
    //  break;

    //case 'D':                            /* maximum depth of suffix tree */
    //  checkNullPointer(c, argv[arg]);
    //  args->D = atoi(argv[++arg]);
    //  if (args->D < 0) {
		  //  args->D = DEFAULT_D;
    //  }
    //  break;
				
    //case 'T':                            /* shulen threshold for a "strong" interval */
    //  checkNullPointer(c, argv[arg]);
    //  args->T = atoi(argv[++arg]);
    //  if (args->T < 0) {
		  //  args->T = DEFAULT_T;
    //  }
    //  break;

    //case 'L':                            /* lcp threshold for calling "process2" */
    //  checkNullPointer(c, argv[arg]);
    //  args->L = atoi(argv[++arg]);
    //  if (args->L < 0) {
		  //  args->L = -1;
    //  }
    //  break;

    case 'v':                           /* print program version test */
      args->v = 1;
      break;
      
    case 'h':                           /* print help */
      args->h = 1;
      break;

    // case 'r':                           /* print reverse strand (default: only forward strand) */
      // args->r = 1;
      // break;

	  case 'P':                            /* fraction of random shustrings excluded from analysis */
      args->P = atof(argv[++arg]);
      checkNullPointer(c, argv[arg]);
      if (args->P < 0) {
		    args->P = 0;
      }
      else if(args->P > 1) {
		    args->P = 1;
      }
      args->P = 1 - args->P;
      break;

    case 'm':                            /* multiplier for maxshulen by chance alone -- used for defining significant intervals */
      args->m = atol(argv[++arg]);
      checkNullPointer(c, argv[arg]);
      if (args->m < 0) {
		    args->m = DEFAULT_m;
      }
      break;

    case 't':                           /* print run-time information */
      args->t = 1;
      break;

    //case 'Q':                           /* compute and print query as a mosaic structure */
    case 'o':                           /* compute and print query as a mosaic structure */
      checkNullPointer(c, argv[++arg]);
      //args->Q = argv[arg];
      args->o = argv[arg];
      break;

    //case 'a':                           /* print alignments (multi-fasta) */
    //  checkNullPointer(c, argv[++arg]);
    //  args->a = argv[arg];
    //  break;

    case 'f':                           /* minimal length of recombination fragment */
      args->f = atoi(argv[++arg]);
      checkNullPointer(c, argv[arg]);
      if (args->f < 0) {
		    args->f = DEFAULT_f;
      }
      break;

    case 'F':                           /* minimal score of recombination fragment */
      args->F = atof(argv[++arg]);
      checkNullPointer(c, argv[arg]);
      if (args->F < 0) {
		    args->F = DEFAULT_F;
      }
      break;

    case 'R':                           /* minimal ratio of score and length of recombination fragment */
      args->R = atof(argv[++arg]);
      checkNullPointer(c, argv[arg]);
      if (args->R < 0 || args->R > 1.) {
		    args->R = DEFAULT_RATIO;
      }
      break;

    case 's':                           /* only "strong" signal (interval) included */
      args->s = 1;
      break;

    case 'A':                            /* maximal distance between two segments to be aligned */
      checkNullPointer(c, argv[arg]);
      args->A = atoi(argv[++arg]);
      if (args->A < 0) {
		    args->A = DEFAULT_A;
      }
      break;

    case 'B':                            /* maximal distance of the segment jsut before or just after a region to be aligned */
      checkNullPointer(c, argv[arg]);
      args->B = atoi(argv[++arg]);
      if (args->B < 0) {
		    args->B = DEFAULT_B;
      }
      break;

    case 'C':                            /* maximal allowed distance between two neighboring intervals in the connected final merged list */
      checkNullPointer(c, argv[arg]);
      args->C = atoi(argv[++arg]);
      if (args->C < 0) {
		    args->C = 0;
      }
      break;

    case 'G':                            /* gap penalty */
      checkNullPointer(c, argv[arg]);
      args->G = (float)atof(argv[++arg]);
      if (args->G < 0) {
		    args->G = DEFAULT_G;
      }
      break;

    case 'M':                            /* match reward */
      checkNullPointer(c, argv[arg]);
      args->M = (float)atof(argv[++arg]);
      if (args->M < 0) {
		    args->M = DEFAULT_M;
      }
      break;

    case 'S':                            /* mismatch penalty */
      checkNullPointer(c, argv[arg]);
      args->S = (float)atof(argv[++arg]);
      if (args->S < 0) {
		    args->S = DEFAULT_S;
      }
      break;

    case 'E':                           /* maximal e-value of an aligned segment to be considered relevant */
      args->E = atof(argv[++arg]);
      checkNullPointer(c, argv[arg]);
      if (args->E < 0) {
        args->E = EVALUE_LOW;
      }
      break;

		case 'X':                           /* maximal e-value of an aligned segment to be considered relevant */
			args->X = atof(argv[++arg]);
			flagExcellentValue = 1;
			checkNullPointer(c, argv[arg]);
			if (args->X < 0) {
        args->X = EVALUE_LOW;
			}
			break;
		
		case 'r':                           /* proportion of the minimal fragment length to be kept in case of overlapping = args->r * maximal fragment */
      args->r = atof(argv[++arg]);
      checkNullPointer(c, argv[arg]);
      if (args->r < 0) {
        args->r = RATIO_MINFRAGMENT;
      }
      break;

    default:
      printf("# unknown argument: %c\n", c);
      args->e = 1;
      return args;
    } // end switch
		
    arg++;  
  } // end while

  //if (args->f != DEFAULT_f && args->F == DEFAULT_F) {
  //  args->F = DEFAULT_RATIO * args->f;
  //}
  
  // allowed stdin?
  if (!args->h && !args->v) {
    if (!args->i) {
       printf("[ERROR]: Query file must be specified using -i option.\n");
       args->e = 1;
    }
    if (!args->j) {
      printf("[ERROR]: Subject file(s) must be specified using -j option.\n");
      args->e = 1;
    }
    //if (!args->Q) {
    if (!args->o) {
      printf("[ERROR]: Output file must be specified using -o option.\n");
      args->e = 1;
    }
    if (args->e) {
        //printf("Usage: %s -i <FILE> -j <FILE(S)> -Q <FILE> [options]\n", progname());
        printf("Usage: %s -i <FILE> -j <FILE(S)> -o <FILE> [options]\n", progname());
        exit(EXIT_FAILURE);
    }
  }
  // args->B should be <= args->A/2
  if (args->A && args->B && args->B > args->A / 2) {

    args->B = args->A / 2;
    printf("Warning[gmos]: Parameter -B has been changed to its upper limit, i.e. A/2: %d.\n", args->B);
  }
  // if args->C is not set, set it to 2 * args->f
  if (args->C == 0) {
    args->C = args->f;
  }

	if (!flagExcellentValue) {
    args->X = EVALUE_LOW;
	}

  return args;
}

void printUsage(char *version){
  printf("Purpose: Compute the mosaic structure of a query when compared to a set of subjects.\n");
  printf("Usage: %s -i <FILE> -j <FILE(S)> -o <FILE> [options]\n", progname());
  printf("Options:\n");
  printf("\t[-f <NUM> minimal length of local alignment; default: %d]\n", DEFAULT_f);	  
  printf("\t[-F <NUM> minimal score of aligned region; default: %.1f]\n", (float)DEFAULT_F);	  
  printf("\t[-R <NUM> minimal ratio of the alignment score and length; default: %.1f]\n", DEFAULT_RATIO);	  
  printf("\t[-A <NUM> maximal distance between two neighboring exact matches in a local alignment; default: %d]\n", DEFAULT_A);	  
  printf("\t[-G <NUM> initial gap penalty; default: %.1f]\n", DEFAULT_G);	  
  printf("\t[-M <NUM> initial match reward; default: %.1f]\n", DEFAULT_M);	  
  printf("\t[-S <NUM> initial mismatch penalty; default: %.1f]\n", DEFAULT_S);	  
	//printf("\t[-E <NUM> e-value threshold; default: %.6e]\n", EVALUE_LOW);
	printf("\t[-X <NUM> e-value threshold (all local alignments, whose E-value is less than the threshold, will be reported in the output); default: %.6e]\n", EVALUE_LOW);
  printf("\t[-r <NUM> scaling factor for the minimal local alignment length; default: %.2lf]\n", RATIO_MINFRAGMENT);
  //  printf("\t[-r print reverse strand; default: only forward strand]\n");			     
  printf("\t[-t print run-time]\n");			     
  printf("\t[-h print this help message]\n");
  printf("\t[-v print information about program]\n");		
}

void printSplash(char *version){
  printf("gmos %s\n",version);
  printf("Written by Mirjana Domazet-Loso\n");
  printf("Distributed under the GNU General Public License\n");
  printf("Please send bug reports to Mirjana.Domazet-Loso@fer.hr\n");
}

/* clean memory allocated by Args structure*/
void freeArgs(Args *args) {
  /* deallocate memory allocated for array of pointers (to subject filenames) */
  if (args->j != NULL) {
    /* when the subject files are read from the directory specified with -d, then the memory allocated for the names 
     * of subject files should be manually allocated;
     * when the subject files are specified with -j, then args->j is just an array of pointers to
     * memory handled by argv[]!!
     */
    //if (args->d) { 
    //  for (i = 0; i < args->subjectFileNumber; i++) {
	   //   free((args->j)[i]);
    //  }
    //}
    free(args->j);
  }

  if (args->i != NULL) { // deallocate the memory for the pointer
    free(args->i);
  }
  free(args);
}




