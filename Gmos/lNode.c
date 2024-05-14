/***** lNode.c *************************************************************
 * Description: Functions for linked lists of intervals.
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "commonSC.h"
#include "eprintf.h"
#include "expectedShulen.h"
#include "intervalStack.h"
#include "interval.h"

#include "queryBTNode.h"
#include "lNode.h"

#if defined(_DEBUG) && defined(WIN) 
#include "leakWatcher.h"
#endif

#if defined(_DEBUG) && defined(WIN) 
  #define new DEBUG_NEW
  #undef THIS_FILE
  static char THIS_FILE[] = __FILE__;
#endif

lNode *getLNode(qNode *p) {
  Int64 i;
  lNode *l;

  // new node in list
  l = (lNode *)/*e*/malloc(sizeof(lNode));
  l->lb = p->lb;
  l->rb = p->rb;
  l->sl = p->sl;

  l->listPositions = (Int64 *)/*e*/malloc(p->sizeListPositions * sizeof(Int64));
  l->sizeListPositions = p->sizeListPositions;
  for (i = 0; i < l->sizeListPositions; i++) {
    l->listPositions[i] = p->listPositions[i];
  } 
  l->next = NULL;
  return l;
}

/* elements of binary interval tree copy to a linked list */
void treeToList(qNode *p, lNode **head, lNode **tail) {
  
  lNode *l;

  // inorder tree traversal
  if (p) {
    // left child
    if (p->left) {
      treeToList(p->left, head, tail); 
    }
    // new node in list
    l = getLNode(p);
    /* correct rb to fit the exact match and not shulen! */
    -- l->rb;
    
    // first element in the list
    if (*head == NULL) {
      *head = *tail = l;
    }
    // every other element of the list
    else {
      (*tail)->next = l;
      *tail = l;
    }
    
    // right child
    if (p->right) {
      treeToList(p->right, head, tail); 
    }
  } // end if
}

/* construct a linked list for each interval tree, i.e. a single list for each pair (Qi, Sj), i=1,..,m; j=1,..,n */
void constructLists(qNode **root, lNode ***head, Int64 numOfQueries, Int64 numOfSubjects) {

  Int64 i, j;
  lNode *tail = NULL;

  for (i = 0; i < numOfQueries; i++) {
    for (j = 0; j < numOfSubjects; j++) {
      treeToList(root[i * numOfSubjects + j], &((*head)[i * numOfSubjects + j]), &tail);
      freeBinTree(root[i * numOfSubjects + j]);
    }
  }
}


/* print each new list, i.e. a single list for each pair (Qi, Sj), i=1,..,m; j=1,..,n */
void printList(lNode **list,  Int64 numOfQueries, Int64 numOfSubjects, FILE *fpout, Int64 *strandBorders, Int64 *leftBorders) {

  Int64 i, j, k;
  lNode *p = NULL;

  for (i = 0; i < numOfQueries; i++) {
    for (j = 0; j < numOfSubjects; j++) {
      p = list[i * numOfSubjects + j];
      fprintf(fpout, "\n------------- Linked list - query %lld, subject %lld -------------\n", (long long)i + 1, (long long)j + 1);
      fflush(fpout);
      //while (p) {
      while (p && p->rb + leftBorders[i] < strandBorders[i]) {
			  fprintf(fpout, "%lld %lld %lld\t", (long long)p->lb, (long long)p->rb, (long long)(p->sl - 1));			
        fflush(fpout);
        // print list of subject "winner" positions
        fprintf(fpout, "Subject: ");
        fflush(fpout);
        for (k = 0; k < p->sizeListPositions - 1; k++) {
          fprintf(fpout, "%lld - %lld, ", (long long)p->listPositions[k], (long long)(p->listPositions[k] + p->rb - p->lb));
          fflush(fpout);
        }
        fprintf(fpout, "%lld - %lld\n", (long long)p->listPositions[k], (long long)(p->listPositions[k] + p->rb - p->lb)); // last position
        fflush(fpout);
        p = p->next;
      } // end while      
    }
  }
}

/* pointers to first elements of |Q| * |S| lists */
lNode **getHead(Int64 numOfQueries, Int64 numOfSubjects) {
	
  lNode **head; // head of |Q| * |S| lists
	Int64 i, j;

	head = /*e*/malloc(numOfQueries * numOfSubjects * sizeof(lNode *));
	/* each element in the list is a pointer to the first element of a linked list */
	for (i = 0; i < numOfQueries; i++) {
    for (j = 0; j < numOfSubjects; j++) {
      head[i * numOfSubjects + j] = NULL; // initially: all lists are empty	
    }
  }
	return head;
}

void freeLNode(lNode *p) {
  free(p->listPositions);
  free(p);
} 

/* free list */
void freeList(lNode *p) {
	
  lNode *q;

  while (p) {
		q = p;
    p = p->next;
		// free node
    freeLNode(q);
	}
	q = NULL;
}


/* free all linked lists */
void freeAllLists(lNode **head, Int64 numOfQueries, Int64 numOfSubjects) {
	
	Int64 i, j;

	for (i = 0; i < numOfQueries; i++) {
    for (j = 0; j < numOfSubjects; j++) {
		  freeList(head[i * numOfSubjects + j]);
    }	
  }
	free(head);
}


