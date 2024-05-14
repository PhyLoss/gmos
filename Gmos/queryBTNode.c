/***** queryBTNode.c *************************************************************
 * Description: Functions for query interval processing - each query is 
 * a node in a binary tree
 * Author: Mirjana Domazet-Loso
 *
 * This file is part of gmos.
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "commonSC.h"
#include "eprintf.h"
#include "intervalStack.h"
#include "expectedShulen.h"

#define VAR2_DECLS /* this is defined only once in this file; all the others just include qBTNode.h without this definition */
#include "queryBTNode.h"

#if defined(_DEBUG) && defined(WIN) 
	#include "leakWatcher.h"
#endif

#if defined(_DEBUG) && defined(WIN)  
  #define new DEBUG_NEW
  #undef THIS_FILE
  static char THIS_FILE[] = __FILE__;
#endif

#define JOIN_SAMESUBJECTS 0 // seems to have no effect

/* get minimal value */
__inline Int64 getMin(Int64 a, Int64 b) {
	return (a <= b ? a : b);
}
	
/* get maximal value */
__inline Int64 getMax(Int64 a, Int64 b) {
	return (a >= b ? a : b);
}

/* allocate new interval */
//qNode *getQNode(Int64 sl, Int64 lb, Int64 rb, Int64 numOfSubjects, Word *subjectIndex) {
qNode *getQNode(Int64 sl, Int64 lb, Int64 rb, Int64 *subjectIdLeaves, Int64 numSIdLeaves) {
  qNode *interval;
  Int64 i;
	//Int64 ns;

  interval = (qNode *)/*e*/malloc(sizeof(qNode));
	++ numQNode;
  interval->sl = sl;  /* longest common prefix */
  interval->lb = lb;    /* left border */
  interval->rb = rb;    /* right border */
  interval->id = ++ nextQNodeId;
  //interval->slAvg = slAvg;  /* longest common prefix */

	/* allocating bit-vector for subject(s) */
	//ns = numOfSubjects / WORDSIZE + 1; // number of words to be allocated for subjects, each word contains WORDSIZE bits
	//interval->subjectIndex = (Word *)/*e*/malloc(sizeof(Word) * ns); // bit-vector of subject-winners
	//for (i = 0; i < ns; i++) {
	//	interval->subjectIndex[i] = subjectIndex[i];
	//}
	interval->left = NULL;
	interval->right = NULL;

  interval->maxListPositions = numSIdLeaves;
  interval->listPositions = (Int64 *)/*e*/malloc(interval->maxListPositions * sizeof(Int64));
  interval->sizeListPositions = numSIdLeaves;
  for (i = 0; i < numSIdLeaves; i ++) {
    interval->listPositions[i] = subjectIdLeaves[i];     
  }
	return interval;
}

/* free interval */
//void freeQNode(qNode *interval) {
void freeQNode(void *_interval) {

	qNode *interval = (qNode *)_interval;
  if (interval) {
//		free(interval->subjectIndex);
    free(interval->listPositions);
		free(interval);
		interval = NULL;
		-- numQNode;
	}
}

/* set values of a query interval */
//void setQNode(qNode *qi, Int64 sl, Int64 lb, Int64 rb, Int64 numOfSubjects, Word *subjectIndex, qNode *left, qNode *right) {
void setQNode(qNode *qi, Int64 sl, Int64 lb, Int64 rb, qNode *left, qNode *right, Int64 *subjectIdLeaves, Int64 numSIdLeaves) {

	Int64 i;

	qi->lb = lb;
	qi->rb = rb;
	qi->sl = sl;
	//qi->slAvg = 0;//slAvg;

	//ns = numOfSubjects / WORDSIZE + 1; // number of words to be allocated for subjects, each word contains WORDSIZE bits
	//for (i = 0; i < ns; i++) {
	//	qi->subjectIndex[i] = subjectIndex[i];
	//}
	qi->left = left;
	qi->right = right;
  
  if (qi->maxListPositions < numSIdLeaves) {
    qi->maxListPositions = numSIdLeaves;
    qi->listPositions = (Int64 *)erealloc(qi->listPositions, qi->maxListPositions * sizeof(Int64));
  }
  qi->sizeListPositions = numSIdLeaves;
  for (i = 0; i < qi->sizeListPositions; i ++) {
    qi->listPositions[i] = subjectIdLeaves[i];     
  }
}

/* return numOfQueries * numOfSubjects binary tree roots */
qNode **getBTQueryIntervals(Int64 numOfQueries, Int64 numOfSubjects) {
	qNode **root; // lists of query intervals, there are |Q| * |S| lists
	Int64 i, j;

	//root = /*e*/malloc(numOfQueries * sizeof(qNode *));
	root = (qNode **)/*e*/malloc(numOfQueries * numOfSubjects * sizeof(qNode *));
	/* each element in the list is a pointer to the first interval in a query-interval list */
	for (i = 0; i < numOfQueries; i++) {
    for (j = 0; j < numOfSubjects; j++) {
      root[i * numOfSubjects + j] = NULL; // initially: all lists are empty	
    }
  }
	return root;
}

/* deallocate tree */
void freeBinTree(qNode *p) {
	if (p) {
		freeBinTree(p->left);
		freeBinTree(p->right);
		freeQNode(p);
		p = NULL;
	}
}

/* free all binary trees */
void freeBTQueryIntervals(qNode **root, Int64 numOfQueries, Int64 numOfSubjects) {
	
	Int64 i, j;

	for (i = 0; i < numOfQueries; i++) {
	    for (j = 0; j < numOfSubjects; j++) {
		    freeBinTree(root[i * numOfSubjects + j]);
      }	
  }
	free(root);
}

/* traverse a binary tree to list all the intervals */
void binTreeTraverse(qNode *p, char **headers, Int64 leftBorder, Int64 ns, Int64 iQuery, Int64 numOfQueries, FILE *fpout) {

	Int64 i;

	if (p) {
		binTreeTraverse(p->left, headers, leftBorder, ns, iQuery, numOfQueries, fpout);

		if (p->lb <= leftBorder) { // left border -- comment this line if both strands need to be displayed
			fprintf(fpout, "%lld %lld %lld\t", (long long)p->lb, (long long)p->rb, (long long)p->sl);
			fflush(fpout);
      // print list of subject "winner" positions
      fprintf(fpout, "Subject: ");
      fflush(fpout);
      for (i = 0; i < p->sizeListPositions - 1; i++) {
        fprintf(fpout, "%lld - %lld, ", (long long)p->listPositions[i], (long long)(p->listPositions[i] + p->rb - p->lb));
        fflush(fpout);
      }
      fprintf(fpout, "%lld - %lld\n", (long long)p->listPositions[i], (long long)(p->listPositions[i] + p->rb - p->lb)); // last position
      fflush(fpout);
			//fprintf(fpout, "\n");
		} // left border -- comment this line if both strands need to be displayed

		binTreeTraverse(p->right, headers, leftBorder, ns, iQuery, numOfQueries, fpout);
    fflush(fpout);
	}	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* function returns 1 if the subject lists of p and p2 are the same */
//int checkSubjectList(qNode *p, qNode *p2, Int64 numOfSubjects) {
//	
//	int retValue = 1;
//	Int64 ns = numOfSubjects / WORDSIZE + 1;
//	Int64 j;
//
//	for (j = 0; j < ns; j++) { 
//		if (p->subjectIndex[j] != p2->subjectIndex[j]) {
//			retValue = 0;
//			break;
//		}
//	}
//	return retValue;
//}

/* join p to its left child l, and delete l */
//void joinLeftChild(qNode *p, qNode *l, Int64 numOfSubjects, Stack *reserveQIStack) {
void joinLeftChild(qNode *p, qNode *l, Stack *reserveQIStack) {
	
	//Int64 i;
	//Int64 numOfSubjWords = numOfSubjects / WORDSIZE + 1;

	//i = numOfSubjWords;
	/* adjust avg shulen across joined intervals */
	//p->slAvg = 0;
	p->lb = l->lb;	
	p->sl = l->sl;	
	//for (i = 0; i < numOfSubjWords; i++) {
	//	p->subjectIndex[i] = l->subjectIndex[i];
	//}		
	push(reserveQIStack, (void *)l);		
}

/* join p to its right child r, and delete r */
void joinRightChild(qNode *p, qNode *r, Stack *reserveQIStack) {
	
	p->rb = r->rb;	
	push(reserveQIStack, (void *)r);		
}

/* join node p to l (l is p's left child) */
//void joinLeftSubTree(qNode *l, qNode *p, Int64 numOfSubjects, Stack *reserveQIStack) {
void joinLeftSubTree(qNode *l, qNode *p, Stack *reserveQIStack) {
	
	l->rb = getMin(l->rb, p->lb - 1); // possibly adjust right border

	if (p->lb - 1 == l->rb && p->sl == (l->sl - (p->lb - l->lb)) ) { // p and n are intervals/nodes that should be joined
		printf("joinLeftSubTree\n");
		if (JOIN_SAMESUBJECTS == 0) {
			//joinLeftChild(p, l, numOfSubjects, reserveQIStack);
			joinLeftChild(p, l, reserveQIStack);
		}
		else {
			// check whether subject lists are the same, and only then join them
			//if (checkSubjectList(p, l /*, numOfSubjects */) == 1) { 
//      ??? different concept now - in a single IT all intervals are of the same subject
				//joinLeftChild(p, l, numOfSubjects, reserveQIStack);
				joinLeftChild(p, l, reserveQIStack);
			//}
		}
	}
	else {
		p->left = l; // p and n are NOT intervals/nodes that should be joined
	}
}

/* join node p to r (r is p's right child) */
//void joinRightSubTree(qNode *r, qNode *p, Stack *reserveQIStack, Int64 numOfSubjects) {
void joinRightSubTree(qNode *r, qNode *p, Stack *reserveQIStack) {
	
	p->rb = getMin(p->rb, r->lb - 1); // possibly adjust right border
	if (r->lb - 1 == p->rb && r->sl == (p->sl - (r->lb - p->lb)) ) { // p and n are intervals/nodes that should be joined
		printf("joinRightSubTree\n");
	
		if (JOIN_SAMESUBJECTS == 0) {
			joinRightChild(p, r, reserveQIStack); // leave p, delete r
		}
		else {
			// check whether subject lists are the same, and only then join them
			//if (checkSubjectList(p, r /*, numOfSubjects*/) == 1) { 
//        ??? different concept now - in a single IT all intervals are of the same subject
				joinRightChild(p, r, reserveQIStack);
			//}
		}
	}
	else {
		p->right = r; // p and n are NOT intervals/nodes that should be joined
	}
}


/* add a new interval/node (n) to the binary tree */
qNode *addNode(qNode *p, qNode *n, Stack *reserveQIStack/*, Int64 numOfSubjects*/) {
	
	if (p == NULL) { // add a new node, also the first node as the root of the tree
		p = n;
	}
	// (1) left subtree
	else if (n->lb < p->lb) { 
		//if (n->rb >= p->rb && (n->sl == p->sl + (p->lb - n->lb))) { 
		if ((n->sl == p->sl + (p->lb - n->lb))  /*&& checkSubjectList(p, n, numOfSubjects)*/) { 
			// n is a superinterval of p ending at the same rb, but only in theory,
			// since the originalp->rb could have been changed before n came along): copy information of n to p, but leave p's children nodes as they are
			// and leave p->rb!!
			p->rb = getMin(p->rb, n->rb);
			//setQNode(p, n->sl, n->lb, p->rb, numOfSubjects, n->subjectIndex, p->left, p->right);
      //setQNode(p, n->sl, n->lb, p->rb, p->left, p->right, n->listPositions, n->sizeListPositions, 0, 0);
      setQNode(p, n->sl, n->lb, p->rb, p->left, p->right, n->listPositions, n->sizeListPositions);
			if (p->left) {
				p->left->rb = getMin(p->left->rb, p->lb - 1);
			}
			/////////////////// check subtree --> !!!!!!!!!!!! rb of any node in the subtree should be < than p->lb
			// THIS IS DONE ONCE THE WHOLE TREE HAS BEEN BUILT --> see correctBT!!
			push(reserveQIStack, (void *)n); // delete new interval
		}
		else {
			// rb in the left subtree cannot be greater than the lb of the root of the subtree
			n->rb = getMin(p->lb - 1, n->rb); //p->rb = min lb in p's subtree - 1 OR its original value

			if (p->left == NULL) { // add a new node at the empty place
				//joinLeftSubTree(n, p, numOfSubjects, reserveQIStack);
				p->left = n; // p and n are NOT intervals/nodes that should be joined
			}
			else { // add a new node between p and its left child
				//p->left = addNode(p->left, n, reserveQIStack, numOfSubjects);
				p->left = addNode(p->left, n, reserveQIStack);
			}
		}
	}
	
	// (2) right subtree
	else if (n->lb > p->lb) { 
		if ((p->sl == n->sl + (n->lb - p->lb)) /*&& checkSubjectList(p, n, numOfSubjects)*/) {
			// p is a superinterval of n (ending at the same rb only if no other node came along that should overlap with p): delete n
			// This should be done here, and not later, since then the connection between these two intervals can be lost in the tree
			// (Otherwise: in some other step p and n would be joined)
			push(reserveQIStack, (void *)n); // delete new interval
		}
		else {
			p->rb = getMin(p->rb, n->lb - 1); //p->rb cannot be greater than any lb in its right subtree
			if (p->right == NULL) { // add a new node at the empty place
				//joinRightSubTree(n, p, reserveQIStack, numOfSubjects);
				p->right = n;
			}

			// add a new node in the right subtree of a parent node
			else { //if (p->right != NULL) {
				//p->right = addNode(p->right, n, reserveQIStack, numOfSubjects);
				p->right = addNode(p->right, n, reserveQIStack);
			}
		}
	}
	// error
	else if (p->lb == n->lb) {
		printf("ERROR when allocating new node!\n"); // THIS IS BECAUSE THE ALLOCATION OF newInterval IS ALSO SOME OLD INTERVAL p from the tree
	}
	return p;
}


/* traverse a binary tree to correct rb of all intervals */
void correctBT(qNode *p, Int64 maxrb, Int64 leftBorder) {

	if (p) { // all nodes, including rev. strand!!
		correctBT(p->left, p->lb - 1, leftBorder);
		if (maxrb != -1) {
			p->rb = getMin(p->rb, maxrb);
		}
		correctBT(p->right, maxrb, leftBorder);
	}	
}

////////////////////// testing purposes functions /////////////////////////////////

/* add a new interval/node (n) to the binary tree (p is root)*/
qNode *insertNode(qNode *p, qNode *n) {

	if (p == NULL) { // add a new node, also the first node as the root of the tree
		p = n;
	}
	else if (n->lb < p->lb) { 
		p->left = insertNode(p->left, n);
	}
	else {
		p->right = insertNode(p->right, n);			
	}

	return p;
}

/* inorder tree traversal */
void inorder(qNode *p, char **subjectNames, int numOfSubjects) {
	
	int i;
	Int64 temp = 1LL;
	
	if (p->left) {
		inorder(p->left, subjectNames, numOfSubjects);
	}
	
	printf("%lld %lld %lld ", (long long)p->lb, (long long)p->rb, (long long)p->sl);
	for (i = 0; i < numOfSubjects; i++) {
		if (/*p->subjectIndex[0] & */  temp) {
			printf("%s ", subjectNames[i]);
		}
		temp = temp << 1;
	}
	printf("\n");
	if (p->right) {
		inorder(p->right, subjectNames, numOfSubjects);
	} 
}

