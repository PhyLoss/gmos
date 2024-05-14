/***** mainGmos.h *******************************************************
 * Description: Header file for functions in main.
 * Author: Mirjana Domazet-Loso
 * 
 * This file is part of gmos.
 *
 *****************************************************************************/ 
#ifndef MAINGMOS_H
#define MAINGMOS_H

//#define MAXLEN 500

/* auxiliary memory deallocation - when query or subject input fails */
void memoryDeallocation2(int queryNumSeq, Args *args, Sequence **query);

/* dealllocate memory for all the objects called in main */ 
void memoryDeallocation(Args *args, SequenceUnion **seqUnion, Sequence ***sArrayAll);

/* read a query fasta file(s) and form a query Sequence object; returns 1 if ok, 0 otherwise */
int readFiles (Sequence **seqArray, Int64 *numSeq, int fileDscr, int numOfFiles, char **fileNames, char ambiguous_char);

/* form a sequence union as a concatenation of all subject seq-s and query 
 * and define its borders (seqBorders and bordersWithinSeq) */
void prepareSeqUnion(SequenceUnion **seqUnion, Sequence **query, Args *args, Sequence **subject, Int64 subjectNumSeq, Int64 queryNumSeq, Sequence ***sArrayAll);

/* form seq. union borders */
Sequence *getSeqUnionBorders(Sequence ***sArrayAll, Int64 subjectNumSeq, Int64 queryNumSeq, SequenceUnion **seqUnion);

#endif // MAINGMOS_H
