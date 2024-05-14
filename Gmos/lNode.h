/***** lNode.h *************************************************************
 * Description: Header file for lists of intervals.
 * Author: Mirjana Domazet-Loso, 2012
 *
 * This file is part of gmos.
 *
 *****************************************************************************/

#ifndef LNODE_H
#define LNODE_H

typedef struct lnode {
  Int64 sl;                         /* shulen attached to the left-most position of this interval */
  Int64 lb;                         /* left border */
  Int64 rb;                         /* right border */

  struct lnode *next;               /* next element in the list */
  Int64 *listPositions;             /* array of subject start positions; it could be > 1 for a query segment starting at lb, but mostly it would be = 1*/
  int sizeListPositions;            /* exact number of elements in listPositions */
} lNode;

lNode *getLNode(qNode *p);
void freeLNode(lNode *p);

/* elements of binary interval tree copy to a linked list */
void treeToList(qNode *p, lNode **head, lNode **tail);

/* construct a linked list for each interval tree, i.e. a single list for each pair (Qi, Sj), i=1,..,m; j=1,..,n */
void constructLists(qNode **root, lNode ***head, Int64 numOfQueries, Int64 numOfSubjects);

/* print each new list, i.e. a single list for each pair (Qi, Sj), i=1,..,m; j=1,..,n */
void printList(lNode **list,  Int64 numOfQueries, Int64 numOfSubjects, FILE *fpout, Int64 *strandBorders, Int64 *leftBorders);

/* het pointers to first elements of |Q| * |S| lists */
lNode **getHead(Int64 numOfQueries, Int64 numOfSubjects);

/* free list */
void freeList(lNode *p);

/* free all linked lists */
void freeAllLists(lNode **head, Int64 numOfQueries, Int64 numOfSubjects);

#endif //LNODE_H
