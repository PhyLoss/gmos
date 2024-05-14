/***** mainGmos.c *************************************************************************************************************
 * Description: Main program - aligining query to subject sequences -- looking only for strong local alignment.
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
 ****************************************************************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "commonSC.h"
#include "eprintf.h"
#include "interface.h"
#include "sequenceData.h"
#include "sequenceUnion.h"
#include "interval.h"
#include "intervalStack.h"
//#include "shulen.h"

#include "interface.h"
#include "queryBTNode.h"
#include "lNode.h"
#include "nw.h"
#include "align.h"
#include "mergeList.h"
#include "lcpTree.h"
#include "mainGmos.h"

#if defined(WIN)
	#include <io.h>
  #include <conio.h>
	#include "fileListWin.h"
#elif defined(UNIX)
	#include <unistd.h>
	#include <fcntl.h>
	#include "fileListUnix.h"
#endif

#if defined(_DEBUG) && defined(WIN) 
	#include "leakWatcher.h"

	#define new DEBUG_NEW
  #undef THIS_FILE
  static char THIS_FILE[] = __FILE__;
#endif

#define MAXLEN_SUBJECTNAME 1024
#define MAXLEN_VERSION 256
#define MAX_FILENAME 1024

  
int main(int argc, char *argv[]){

	FILE *fpout = NULL;
	FILE *fqout = NULL;
	FILE *fQout = NULL;
	FILE *faout = NULL;
	FILE *fLen = NULL;
  char fLenName[MAX_FILENAME + 1];
  char fAName[MAX_FILENAME + 1];
  char fQName[MAX_FILENAME + 1];
	int subjectDscr, queryDscr;
	Int64 subjectNumSeq, queryNumSeq;
	char version[MAXLEN_VERSION + 1];
	Sequence **query = NULL, **subject = NULL; 
	Sequence **sArrayAll = NULL; 
	Args *args;
	SequenceUnion *seqUnion = NULL;
  int i = 0;
	clock_t end, start;

	#if defined(_DEBUG) && defined(WIN) 
		_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
	#endif  	

	#if VER32
		setprogname2("gmos");
	#else // VER64
		setprogname2("gmos64");
	#endif
#if defined(WIN)
  strcpy(version, "1.0");	
#else
  strcpy(version, "1.0");	
#endif
  start = clock();

	// get arguments (including query and subject files)
  args = getArgs(argc, argv);
  if (args->h == 1) { 
		printUsage(version);
    #if defined(WIN) 
      getchar();
    #endif
		return 0; 
	}
	else if (args->v){
		printSplash(version);
    #if defined(WIN) 
      getchar();
    #endif
		return 0;
  }
	else if (args->e){
    #if defined(WIN) 
      getchar();
    #endif
		exit(EXIT_FAILURE);
  }

	/*
	1) open fasta input files: 1 or more query files (Q1 .. Qm) and 1 or more subject files (S1.. Sn)
	2) read from input files and generate sequences
	--> concatenate all the seq-s in one big sequenceUnion
	3) form suffix array (sa) from all sequences (from sequenceUnion)
	4) form longest common prefix array (lcp) from sa
	5) simulate sufix tree traversal using algorithm by Abouelhoda et al.
	--> every node in the tree has to have flag that is set to Q and/or to one or more Si; 
	--> if the Q/Si flag is set, it means that the sequence Q/Si is passing through that node
	--> construct n interval trees for each Qj: Ii for each Si; 
	    each Ii can contain only "strong" intervals - list of strong local alignment intervals between Qj and Si 
		(store positions of both Qj and Sj corresponding intervals)

	6) the basic output(looking at local shustrings): for each Qj (j = 1,..,m) print list of "strong" intervals
		(strong local alignment) between Qj and every Si (i = 1,..,n)
	7) print for each Sj the best suggested solution as a combinaiton of local alignments between Qj and all Si
		(alternative? : enable a user to create her own best solution)
	*/		
	
	
	// STEP (1) open input files; the number of files expected is >= 2 : 1 query file and 1 (or more) subject files

	/* open query file and read data */
	queryDscr = -1; /* initialize */
	if (args->queryFileNumber == 0) { /* this means that the query comes from the stdin; however, this is not yet allowed ! */
    queryDscr = 0;
		args->i = 0;
		args->queryFileNumber = 1; /* one input file, but in fact it is stdin */
	}
	/* read data from query file(s); array of queries consists of Sequence objects where each object represents one file */
	query = (Sequence **)/*e*/malloc(sizeof(Sequence *) * args->queryFileNumber);
  if (readFiles(query, &queryNumSeq, queryDscr, args->queryFileNumber, args->i, AMBIGUOUS_CHAR_QUERY) == 0) { // error reading query file(s)
    memoryDeallocation2(args->queryFileNumber, args, query); // not queryNumSeq, but args->queryFileNumber, since at this point multiple query seq-s from a single file haven't been divided in multiple Sequence-s
    exit(EXIT_FAILURE);
  }

	/* check whether subject file(s) exist; if it doesn't, and the subject subdirectory is not specified, then consider the stdin as the subject file */
	subjectDscr = -1; /* initialize */
	if (args->subjectFileNumber == 0 /*&& !args->d*/) { /* this means that the subject comes from the stdin; however, this is not yet allowed ! */
    subjectDscr = 0;
		args->j = 0;
		args->subjectFileNumber = 1; /* one input file, but in fact it is stdin */
	}
	//else if (args->d) { /* subject files are in the subdirectory (of the current directory) specified with -d option */
	//	#if WIN
	//		args->j = listFilesWin(args->d, &i); //(char **)/*e*/realloc(args->j, numberOfFiles * sizeof (char *));
	//	#elif UNIX
	//		args->j = listFilesUnix(args->d, &i); //(char **)/*e*/realloc(args->j, numberOfFiles * sizeof (char *));
	//	#endif
	//	args->subjectFileNumber = i;		
	//}
	else { /* subject files are listed with -j option, so args->subjectFileNumber and arg->j are already set */
	}

#if DEBUG
  if (fqout) {
    int j = 0;
		for (i = 0; i < queryNumSeq; i++) { // print all headers
			//fprintf(fqout,"%-80.80s\n%-80.80s\n", sArrayAll[i]->headers[0], sArrayAll[i]->seq); 
			fprintf(fqout,"%-120.120s\n", query[0]->headers[i]); 
      for (; query[0]->seq[j] != BORDER; j++) {
        fputc(query[0]->seq[j], fqout);
      }
      fflush(fqout);        
      if (query[0]->seq[j] == BORDER) {
        fputc('\n', fqout);
        ++ j;
      }
		}
    fputc('\n', fqout);
  }
#endif


	/* read data from subject file(s); array of subjects consists of Sequence objects where each object represents one file */
	subject = (Sequence **)/*e*/malloc(sizeof(Sequence *) * args->subjectFileNumber);
  for (i = 0; i < args->subjectFileNumber; i++) {
    subject[i] = NULL;
  }
  if (readFiles (subject, &subjectNumSeq, subjectDscr, args->subjectFileNumber, args->j, AMBIGUOUS_CHAR_SUBJECT) == 0) { // error reading subject file(s)
    freeSequenceArray(args->subjectFileNumber, subject); // not subjectNumSeq, but args->subjectFileNumber, since at this point multiple subject seq-s from a single file haven't been divided in multiple Sequence-s
    memoryDeallocation2(args->queryFileNumber, args, query); // not queryNumSeq, but args->queryFileNumber, since at this point multiple query seq-s from a single file haven't been divided in multiple Sequence-s
    exit(EXIT_FAILURE);
  }

//#if 1
//  if (fqout) {
//    int j = 0, k;
//    i = 0;
//    for (k = 0; k < args->subjectFileNumber; k++) {
//		  for (; i < subjectNumSeq; i++) { // print all headers
//			  //fprintf(fqout,"%-80.80s\n%-80.80s\n", sArrayAll[i]->headers[0], sArrayAll[i]->seq); 
//			  fprintf(fqout,"%-120.120s\n", subject[k]->headers[i]); 
//        for (; subject[k]->seq[j] != BORDER; j++) {
//          fputc(subject[k]->seq[j], fqout);
//        }
//        fflush(fqout);        
//        if (query[0]->seq[j] == BORDER) {
//          fputc('\n', fqout);
//          ++ j;
//        }
//		  }
//      fputc('\n', fqout);
//    }
//  }
//#endif
//

	/* STEP (2) form a sequence union:
	 * 1) concatenate all the subject seq-s to the query
	 * 2) define borders in seqBorders and bordersWithinSeq
	*/
	prepareSeqUnion(&seqUnion, query, args, subject, subjectNumSeq, queryNumSeq, &sArrayAll);

	// open output files
	//if (args->o) { 
	if (args->O) { 
	  fpout = efopen(args->O, "w");
	}
	if (args->q) { 
	  fqout = efopen(args->q, "w");
	}
	//if (args->Q) { 
	if (args->o) { // fasta file
    // filename of the file containing seq. lengths should be the same as fQout except for the extension ".len"
    i = strlen(args->o) - strlen(strchr(args->o, '.'));
    // len filename and alignment filename are constructed from args->o and these files are automatically added
    strncpy(fQName, args->o, i);
    fQName[i] = '\0';
    strcat(fQName, ".txt");
    fQout = efopen(fQName, "w"); // txt file!

    strncpy(fLenName, args->o, i);
    fLenName[i] = '\0';
    strcat(fLenName, ".len");
    fLen = efopen(fLenName, "w");

    strncpy(fAName, args->o, i);
    fAName[i] = '\0';
    strcat(fAName, ".fasta");
    faout = efopen(fAName, "w");
	}

	/* print to the output stream */
#if DEBUG
    if (fpout) {
		  fprintf(fpout,"\nQUERY(QUERIES):\n");
      fflush(fpout);
		  for (i = 0; i < queryNumSeq; i++) { // print all headers
			  fprintf(fpout,"%-80.80s\n%-80.80s\n", sArrayAll[i]->headers[0], sArrayAll[i]->seq); 
        fflush(fpout);
		  }

		  fprintf(fpout,"SUBJECT(S):\n");
		  for (i = 0; i < subjectNumSeq; i++) {
			  fprintf(fpout,"%-80.80s\n%-80.80s\n", sArrayAll[i + queryNumSeq]->headers[0], sArrayAll[i + queryNumSeq]->seq); 
        fflush(fpout);
		  }
    }
	#endif

	// fpout - intervals output
	getLcpTreeShulens(fpout, args, seqUnion, fqout, fQout, faout, fLen); 
	
  /* closing files */
	//if (args->o) { 
	if (args->O) { 
		fclose(fpout);
	}
	//if (args->Q) { 
	if (args->o) { 
		fclose(fQout);
    fclose(fLen);
		fclose(faout);
	}
	if (args->q) { 
		fclose(fqout);
	}

	// print run-time
	end = clock();
	if (args->t) {
		printf( "\nRunning time: %.2f seconds.\n", (double)(end - start) / CLOCKS_PER_SEC);
	}

	/* memory deallocation */ 
	memoryDeallocation(args, &seqUnion, &sArrayAll);

  #if defined(_DEBUG) && defined(WIN) 
		_CrtDumpMemoryLeaks();
	#endif  	
  
  #if defined(WIN) 
    printf("Finished.\n");
    getchar();
  #endif
	return 0;
}

/****************************************************************************************************************************************/

#if defined(WIN)
  static 
#endif
/* auxiliary memory deallocation - when query or subject input fails */
void memoryDeallocation2(int queryNumSeq, Args *args, Sequence **query) {

    freeSequenceArray(queryNumSeq, query);
    freeArgs(args);
    free(progname());
		#if WIN
      fgetc(stdin);
    #endif
}

/* dealllocate memory for all the objects called in main */ 
#if defined(WIN)
static 
#endif
void memoryDeallocation(Args *args, SequenceUnion **seqUnion, Sequence ***sArrayAll) {
	
	freeSequenceArray((*seqUnion)->numOfSubjects + (*seqUnion)->numOfQueries, *sArrayAll);
	//freeSequenceArray(args->subjectFileNumber, subject); /* moved to prepareSeqUnion */
	//freeSequenceArray(args->queryFileNumber, query);
	
	/* the list of subject file names is deallocated with the call of freeArgs(args)!! */
	freeArgs(args);
	freeSequenceUnion(*seqUnion, (*seqUnion)->numOfSubjects + (*seqUnion)->numOfQueries);	
	free(progname());
}

/* read data from query (or subject) file(s); array of queries ("query") (or subjects) consists 
 * of Sequence objects where each object represents one file */
#if defined(WIN)
static 
#endif
//int readFiles (Sequence **seqArray, Int64 *numSeq, int fileDscr, int numOfFiles, char **fileNames) {
int readFiles (Sequence **seqArray, Int64 *numSeq, int fileDscr, int numOfFiles, char **fileNames, char ambiguous_char) {
	int i;
  //int j;
  int retValue = 1;

	*numSeq = 0;
	for (i = 0; i < numOfFiles; i++) {		
		if (fileDscr == -1) {
			fileDscr = open((const char *)fileNames[i], 0);
		}
		//else if (queryDscr == 0) { //stdin
		//	// do nothing
		//}
		
    /* if a file cannot be open, deallocate memory and terminate the program */
		if (fileDscr < 0) { // memory deallocation is done later
			//for (j = 0; j < i; j++) {
			//	freeSequence(seqArray[j]);					
			//}
			printf("[ERROR] Could not open file %s\n", fileNames[i]);	
      retValue = 0;
      break;
		} 
		
    /* if it's an empty file or an error file, then skip the rest of the loop */   
		if (lseek(fileDscr, 0L, SEEK_END) <= 0L) { /* -1 for error, 0 for the file beginning */
			if (fileDscr > 0) {
				printf("[ERROR] File: %s is empty or file error!\n", (const char *)fileNames[i]);
			}
			else { // no input from stdin
				printf("[ERROR] No input file specified.\nUsage: %s -i <FILE> -j <FILE(S)> -Q <FILE> [options]\n", progname());
			}
      retValue = 0;
      break;
		}
		else { /*	if a file is ok, then read its sequence(s) */		
			lseek(fileDscr, 0L, SEEK_SET); /* return to the file beginning */
			seqArray[i] = readFasta(fileDscr, ambiguous_char); 
			close(fileDscr);
			*numSeq += seqArray[i]->numSeq;
		}
		fileDscr = -1;
	} /* end for */	
  return retValue;
}


/* form a sequence union as a concatenation of all subject seq-s to the query 
 * and define its borders (seqBorders and bordersWithinSeq)
 */
#if defined(WIN)
static 
#endif
void prepareSeqUnion(SequenceUnion **seqUnion, Sequence **query, Args *args, Sequence **subject, Int64 subjectNumSeq, Int64 queryNumSeq, Sequence ***sArrayAll) {

	Sequence **sArray, *cat;
	Int64 i, j, k;

	/* form a sequence union:
	 * 1) concatenate all the subject seq-s to the query sequences
	 * 2) define borders in seqBorders and bordersWithinSeq
	*/
	*seqUnion = getSequenceUnion(subjectNumSeq + queryNumSeq, subjectNumSeq, queryNumSeq); /* subjectNumSeq = actual number of subjects in subject files */ 

	/* if subject[i] or query[i] is consisted of multiple strains, then divide those strains in separate Sequence objects - each strain in one Sequence object */
	sArray = NULL; /* array of strains from one subject */
	*sArrayAll = (Sequence **)/*e*/malloc((size_t)(subjectNumSeq + queryNumSeq) * sizeof(Sequence *)); /* array of strains from all subjects and queries */
	
	// queries
	for (i = 0, j = 0; i < args->queryFileNumber; i++) {		/* j: position in the sArrayAll which contains strains from all queries */
		if (query[i]) { /* not an empty file */			
			sArray = getArrayOfSeq(query[i]); /* sArray is the list of strains (Sequence objects) from i-th query */
			for (k = 0; k < query[i]->numSeq; k++) {
				(*sArrayAll)[j++] = sArray[k];
			}
			free(sArray); /* sArray[k] should not be deallocated, since they are just moved to *sArrayAll */
		}
	}
	freeSequenceArray(args->queryFileNumber, query);
	query = NULL;

	// subjects
	for (i = 0; i < args->subjectFileNumber; i++) {		/* j: position in the sArrayAll which contains strains from all subjects */
		if (subject[i]) { /* not an empty file */			
			sArray = getArrayOfSeq(subject[i]); /* sArray is the list of strains (Sequence objects) from i-th subject */
			for (k = 0; k < subject[i]->numSeq; k++) {
				(*sArrayAll)[j++] = sArray[k];
			}
			free(sArray); /* sArray[k] should not be deallocated, since they are just moved to *sArrayAll */
		}
	}
	freeSequenceArray(args->subjectFileNumber, subject);
	subject = NULL;

	/* cat = concatenation of all subject seq-s; also set borders within seq */
	cat = getSeqUnionBorders(sArrayAll, subjectNumSeq, queryNumSeq, seqUnion);
	//(*seqUnion)->seqUnion = catSeq(query, cat, 'Q');	
	(*seqUnion)->seqUnion = cat;	

	(*seqUnion)->numOfSubjects = subjectNumSeq;
	(*seqUnion)->numOfQueries = queryNumSeq;
	if (queryNumSeq + subjectNumSeq > 1) {
		//freeSequence(cat);
	}
	
	(*seqUnion)->len = (*seqUnion)->seqBorders[subjectNumSeq + queryNumSeq - 1] + 1;
	(*seqUnion)->seqUnion->queryStart = 0;
	(*seqUnion)->seqUnion->queryEnd = (*seqUnion)->seqBorders[queryNumSeq - 1];
}

/* form seq. union borders */
#if defined(WIN)
static 
#endif
Sequence *getSeqUnionBorders(Sequence ***sArrayAll, Int64 subjectNumSeq, Int64 queryNumSeq, SequenceUnion **seqUnion) {

	Int64 j, k, numChar, numStrain; 
	Sequence *cat = NULL, *catOld;
	char type;
	Int64 lb = 0; /* left border */

	for (j = 0; j < subjectNumSeq + queryNumSeq; j++) {
		prepareSeq((*sArrayAll)[j]);
    /* concatenate the strain[j] to the concatenated strains of all the subjects so far */
		if (j == 0) {
			cat = (*sArrayAll)[0];
		}
		else {
			catOld = cat;
			// first add queries, and then subjects 
			type = (j < queryNumSeq) ? 'Q' : 'S';
			cat = catSeq(catOld, (*sArrayAll)[j], type);
			if (j > 1) { /* when i == 1, catOld points to sArrayAll[0] and this should not be deallocated */
				freeSequence(catOld);
			}
		}
		/* set borders of the sequence union */
		(*seqUnion)->seqBorders[j] = lb + (*sArrayAll)[j]->len - 1; // len icludes border

		numStrain = (*sArrayAll)[j]->numSeq * 2; // include reverse strand
		(*seqUnion)->bordersWithinSeq[j] = (Int64 *)/*e*/malloc(sizeof(Int64) * (size_t)numStrain); /* fwd [ + rev] strand */
		for (k = 0; k < numStrain; k++) {
			//(*seqUnion)->bordersWithinSeq[j + 1][k] = (*seqUnion)->seqBorders[j] + (*sArrayAll)[j]->borders[k] + 1; // len icludes border
			(*seqUnion)->bordersWithinSeq[j][k] = lb + (*sArrayAll)[j]->borders[k]; // len icludes border
		}		
		lb = (*seqUnion)->seqBorders[j] + 1; /* new left border */
		
		/* gc content */
		(*seqUnion)->gc[j] = (double)(*sArrayAll)[j]->freqTab[0]['G'] + (*sArrayAll)[j]->freqTab[0]['C'];
		numChar = (Int64)(*seqUnion)->gc[j] + (*sArrayAll)[j]->freqTab[0]['A'] + (*sArrayAll)[j]->freqTab[0]['T'];
		(*seqUnion)->gc[j] /= numChar;
	}
	return cat;
}

