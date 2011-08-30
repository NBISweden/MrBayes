/*
 *  MrBayes 3
 *
 *  (c) 2002-2010
 *
 *  John P. Huelsenbeck
 *  Dept. Integrative Biology
 *  University of California, Berkeley
 *  Berkeley, CA 94720-3140
 *  johnh@berkeley.edu
 *
 *  Fredrik Ronquist
 *  Swedish Museum of Natural History
 *  Box 50007
 *  SE-10405 Stockholm, SWEDEN
 *  fredrik.ronquist@nrm.se
 *
 *  With important contributions by
 *
 *  Paul van der Mark (paulvdm@sc.fsu.edu)
 *  Maxim Teslenko (maxim.teslenko@nrm.se)
 *
 *  and by many users (run 'acknowledgements' to see more info)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details (www.gnu.org).
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include "bayes.h"
#include "best.h"
#include "mb.h"
#include "mcmc.h"
#include "mbmath.h"
#include "model.h"
#include "globals.h"
#include "tree.h"
#include "utils.h"
#include "command.h"

const char* const svnRevisionTreeC="$Rev$";   /* Revision keyword which is expended/updated by svn on each commit/update*/

/* local prototypes */
void    DatedNodeDepths (TreeNode *p, MrBFlt *nodeDepths, int *index);
void    DatedNodes (TreeNode *p, TreeNode **datedTips, int *index);
void    MarkDatedSubtreeNodes (TreeNode *p);
int     NConstrainedTips (TreeNode *p);
int     NDatedTips (TreeNode *p);
void	PrintNode (char **s, int *len, TreeNode *p, int isRooted);
void    ResetPolyNode (PolyNode *p);
void	ResetTreeNode (TreeNode *p);
void    SetNodeDepths (Tree *t);

/* local global variable */
char    noLabel[] = "";


/* AddToTreeList: Add tree at end of tree list */
int AddToTreeList (TreeList *treeList, Tree *tree)

{

	TreeListElement		*listElement = (TreeListElement *) SafeCalloc (1, sizeof(TreeListElement));
	if (!listElement)
		return (ERROR);

	listElement->order = (int *) SafeCalloc (tree->nIntNodes-1, sizeof(int));
	if (!listElement->order)
		return (ERROR);
	listElement->next = NULL;

	if (treeList->last == NULL)
		treeList->last = treeList->first = listElement;
	else
		{
		treeList->last->next = listElement;
		treeList->last = listElement;
		}

	if (tree->isRooted)
		StoreRTopology (tree, listElement->order);
	else
		StoreUTopology (tree, listElement->order);

	return (NO_ERROR);
}





/* AllocatePolyTree: Allocate memory space for a polytomous tree */
PolyTree *AllocatePolyTree (int numTaxa)
{
	int			i;
	PolyTree	*pt;

	pt = (PolyTree *) SafeCalloc (1, sizeof (PolyTree));
	if (!pt)
		return (NULL);

	pt->memNodes = 2*numTaxa;  
	pt->nodes = (PolyNode *) SafeCalloc (2*numTaxa, sizeof(PolyNode));
	pt->allDownPass = (PolyNode **) SafeCalloc (3*numTaxa, sizeof (PolyNode *));
	pt->intDownPass = pt->allDownPass + 2*numTaxa;
	if (pt->nodes == NULL || pt->allDownPass == NULL)
		{
		free (pt->nodes);
		free (pt->allDownPass);
		free (pt);
		return (NULL);
		}

	/* initialize nodes and set index and memoryIndex */
	for (i=0; i<2*numTaxa; i++)
		{
		ResetPolyNode(&pt->nodes[i]);
        pt->nodes[i].memoryIndex = i;
        pt->nodes[i].index = i;
		}

    /* initialize tree properties */
    pt->nNodes = pt->nIntNodes = 0;
    pt->root = NULL;
    pt->brlensDef = NO;
    pt->isRooted = NO;
    pt->isClock = NO;
    pt->isCalibrated = NO;
    pt->isRelaxed = NO;
    pt->clockRate = 0.0;
    strcpy(pt->name,"");

    /* initialize bitsets */
    pt->bitsets = NULL;
    
    /* initialize relaxed clock parameters */
    pt->nESets = 0;
    pt->nEvents = NULL;
    pt->position = NULL;
    pt->rateMult = NULL;
    pt->eSetName = NULL;

    pt->nBSets = 0;
    pt->effectiveBrLen = NULL;
    pt->bSetName = NULL;

    /* initialize population size set parameters */
    pt->popSizeSet = NO;
    pt->popSize = NULL;
    pt->popSizeSetName = NULL;

	return (pt);
}





/* AllocatePolyTreeRelClockParams: Allocate space for relaxed clock parameters */
int AllocatePolyTreeRelClockParams (PolyTree *pt, int nBSets, int nESets)
{
    int     i;

    /* free previous clock params if any */
    FreePolyTreeRelClockParams (pt);

    /* set number of sets */
    pt->nBSets = nBSets;
    pt->nESets = nESets;

    /* we do not allocate space for the actual names here; these will be NULL pointers */

    /* take care of branch params */
    if (pt->nBSets > 0)
		{
		pt->bSetName = (char **) SafeCalloc (pt->nBSets, sizeof (char *));
        pt->effectiveBrLen = (MrBFlt **) SafeCalloc (pt->nBSets, sizeof (MrBFlt *));
		for (i=0; i<pt->nBSets; i++)
            pt->effectiveBrLen[i] = (MrBFlt *) SafeCalloc (pt->memNodes, sizeof(MrBFlt));
		}
	
    /* take care of breakpoint params */
    if (pt->nESets > 0)
		{
		pt->eSetName = (char **) SafeCalloc (pt->nESets, sizeof(char *));
		pt->nEvents = (int **) SafeCalloc (pt->nESets, sizeof(int *));
        pt->position = (MrBFlt ***) SafeCalloc (pt->nESets, sizeof(MrBFlt **));
        pt->rateMult = (MrBFlt ***) SafeCalloc (pt->nESets, sizeof(MrBFlt **));
		for (i=0; i<pt->nESets; i++)
			{
    		pt->nEvents[i] = (int *) SafeCalloc (pt->memNodes, sizeof(int));
            pt->position[i] = (MrBFlt **) SafeCalloc (pt->memNodes, sizeof(MrBFlt *));
            pt->rateMult[i] = (MrBFlt **) SafeCalloc (pt->memNodes, sizeof(MrBFlt *));
		    }
        }

    return (NO_ERROR);
}





/* AllocatePolyTreePartitions: Allocate space for and set partitions for polytomous tree */
int AllocatePolyTreePartitions (PolyTree *pt)
{
	int			i, nLongsNeeded, numTaxa;

    /* get some handy numbers */
    numTaxa = pt->memNodes/2;
	nLongsNeeded = (numTaxa -1) / nBitsInALong + 1;

    /* allocate space */
    pt->bitsets = (SafeLong *) SafeRealloc ((void *)pt->bitsets, pt->memNodes*nLongsNeeded*sizeof(SafeLong));
	if (pt->bitsets == NULL)
		return (ERROR);
    for (i=0; i<pt->memNodes*nLongsNeeded; i++)
        pt->bitsets[i] = 0;
	
    /* set node partition pointers */
    for (i=0; i<pt->memNodes; i++)
		pt->nodes[i].partition = pt->bitsets + i*nLongsNeeded;

	/* clear and set partitions; if the tree is empty, nothing is set */
	ResetPolyTreePartitions(pt);
    
	return (NO_ERROR);
}





/* AllocateTree: Allocate memory space for a tree (unrooted or rooted) */
Tree *AllocateTree (int numTaxa)
{
	int		i;
	Tree	*t;
	
	t = (Tree *) SafeCalloc (1, sizeof (Tree));
	if (t == NULL)
		return NULL;

	/* initialize basic tree properties */
    t->memNodes = 2*numTaxa;
	strcpy (t->name, "");
    
    t->isRooted = NO;
    t->isClock = NO;

    t->checkConstraints = NO;
    t->nConstraints = 0;
    t->nLocks = 0;
    t->isCalibrated = NO;
    t->nNodes = t->nIntNodes = 0;
    t->nRelParts = 0;
    t->relParts = NULL;

    /* initialize pointers */
	t->bitsets = NULL;
	t->flags = NULL;
    t->constraints = NULL;

    /* allocate and initialize nodes and node arrays (enough for both rooted and unrooted trees) */
    t->nNodes = 0;
    t->nIntNodes = 0;
    if ((t->nodes = (TreeNode *) SafeCalloc (2*numTaxa, sizeof (TreeNode))) == NULL)
		{
		free (t);
		return NULL;
		}
	if ((t->allDownPass = (TreeNode **) SafeCalloc (3*numTaxa, sizeof (TreeNode *))) == NULL)
		{
		free (t->nodes);
		free (t);
		return NULL;
		}
	t->intDownPass = t->allDownPass + t->memNodes;
	
	/* initialize nodes and set index and memoryIndex */
	for (i=0; i<t->memNodes; i++)
		{
		ResetTreeNode(&t->nodes[i]);
		t->nodes[i].memoryIndex = i;
		t->nodes[i].index = i;
        }

	return t;
}





/* AllocateFixedTree: Allocate memory space for a fixed unrooted or rooted tree */
Tree *AllocateFixedTree (int numTaxa, int isRooted)
{
	int		i;
	Tree	*t;
	
	t = (Tree *) SafeCalloc (1, sizeof (Tree));
	if (t == NULL)
		return NULL;

	/* initialize basic tree properties */
    if (isRooted == YES)
        t->memNodes = 2*numTaxa;
    else
        t->memNodes = 2*numTaxa - 2;
	strcpy (t->name, "");
    
    t->isRooted = isRooted;
    t->isClock = NO;

    t->checkConstraints = NO;
    t->nConstraints = 0;
    t->nLocks = 0;
    t->isCalibrated = NO;
    t->nNodes = t->nIntNodes = 0;
    t->nRelParts = 0;
    t->relParts = NULL;

    /* initialize pointers */
	t->bitsets = NULL;
	t->flags = NULL;
    t->constraints = NULL;

    /* allocate and initialize nodes and node arrays (enough for both rooted and unrooted trees) */
    if (t->isRooted)
        {
        t->nNodes = 2*numTaxa;
        t->nIntNodes = numTaxa - 1;
        }
    else
        {
        t->nNodes = 2*numTaxa - 2;
        t->nIntNodes = numTaxa - 2;
        }
    if ((t->nodes = (TreeNode *) SafeCalloc (t->nNodes, sizeof (TreeNode))) == NULL)
		{
		free (t);
		return NULL;
		}
	if ((t->allDownPass = (TreeNode **) SafeCalloc (t->nNodes + t->nIntNodes, sizeof (TreeNode *))) == NULL)
		{
		free (t->nodes);
		free (t);
		return NULL;
		}
	t->intDownPass = t->allDownPass + t->nNodes;
	
	/* initialize nodes and set index and memoryIndex */
	for (i=0; i<t->memNodes; i++)
		{
		ResetTreeNode(&t->nodes[i]);
		t->nodes[i].memoryIndex = i;
		t->nodes[i].index = i;
        }

	return t;
}





/* AllocateTreePartitions: Allocate space for and set partitions for tree */
int AllocateTreePartitions (Tree *t)
{
	int			i, nLongsNeeded, numTaxa;
	TreeNode	*p;
	
    /* get some handy numbers */
	if (t->isRooted == YES)
		numTaxa = t->nNodes - t->nIntNodes - 1;
	else
		numTaxa = t->nNodes - t->nIntNodes;
	nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;

	/* reallocate space */
    t->bitsets = (SafeLong *) SafeRealloc ((void *) t->bitsets, (size_t)(t->nNodes*nLongsNeeded*sizeof(SafeLong)));
    if (!t->bitsets)
		return (ERROR);
	
    /* clear bit fields */
    for (i=0; i<t->nNodes*nLongsNeeded; i++)
		t->bitsets[0] = 0;
        
    /* set node pointers to bit fields */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
		p->partition = t->bitsets + i*nLongsNeeded;
		}

	/* set partition specifiers for terminals */
	ResetTreePartitions(t);
    
	return (NO_ERROR);
}





int AreTopologiesSame (Tree *t1, Tree *t2)

{
	int			i, j, k, nLongsNeeded, nTaxa, revertBits;
	SafeLong	*bitsets, *mask;
	TreeNode	*p, *q;

    if (t1->nNodes != t2->nNodes)
		return (NO);
	if (t1->nIntNodes != t2->nIntNodes)
		return (NO);
	
	if (t1->isRooted == NO && (t1->root->index != t2->root->index))
		revertBits = YES;
	else
		revertBits = NO;
	
	if (t1->isRooted == YES)
		nTaxa = t1->nNodes - t1->nIntNodes - 1;
	else
		nTaxa = t1->nNodes - t1->nIntNodes;
	
    /* allocate space */
    nLongsNeeded = (nTaxa - 1) / nBitsInALong + 1;
    bitsets = (SafeLong *) SafeCalloc (4*nLongsNeeded*nTaxa+nLongsNeeded, sizeof(SafeLong));
    mask = bitsets + 4*nLongsNeeded*nTaxa;
	
    /* set mask */
    for (i=0; i<nTaxa; i++)
        SetBit(i, mask);
	
    /* set partition pointers */
	for (i=0; i<t1->nNodes; i++) 
		{
		p = t1->allDownPass[i];
		q = t2->allDownPass[i];
		p->partition = bitsets + i*2*nLongsNeeded;
		q->partition = p->partition + nLongsNeeded;
		}
    
    /* set terminal partitions */
    for (i=0; i<t1->nNodes; i++)
        {
        p = t1->allDownPass[i];
        q = t2->allDownPass[i];
        if (p->right == NULL)
            SetBit (p->index, p->partition);
        if (q->right == NULL)
            SetBit (q->index, q->partition);
        }
	
    /* set internal partitions */
    for (i=0; i<t1->nIntNodes; i++)
        {
        p = t1->intDownPass[i];
        q = t2->intDownPass[i];
        for (j=0; j<nLongsNeeded; j++)
            {
            p->partition[j] = p->left->partition[j] | p->right->partition[j];
            q->partition[j] = q->left->partition[j] | q->right->partition[j];
            }
        }

	/* check for congruence */
	for (i=0; i<t1->nIntNodes; i++)
		{
		p = t1->intDownPass[i];
		if (t1->isRooted == NO && IsBitSet(t2->root->index,p->partition))
                    FlipBits(p->partition,nLongsNeeded, mask);
		for (j=0; j<t2->nIntNodes; j++)
			{
			q = t2->intDownPass[j];
			for (k=0; k<nLongsNeeded; k++)
				{
				if (p->partition[k] != q->partition[k])
					break;
				}
			if (k == nLongsNeeded)
				break;
			}
		if (j == t2->nIntNodes)
			{
			free (bitsets);
			return (NO);			
			}
		}

    free (bitsets);
    return (YES);
}





int AreTreesSame (Tree *t1, Tree *t2)

{
	int			i, j, k, nLongsNeeded, nTaxa, revertBits;
	SafeLong	*bitsets, *mask;
	TreeNode	*p, *q;

	extern void ShowNodes(TreeNode*, int, int);
	
	if (t1->nNodes != t2->nNodes)
		return (NO);
	if (t1->nIntNodes != t2->nIntNodes)
		return (NO);
	
	if (t1->isRooted == NO && (t1->root->index != t2->root->index))
		revertBits = YES;
	else
		revertBits = NO;
	
	if (t1->isRooted == YES)
		nTaxa = t1->nNodes - t1->nIntNodes - 1;
	else
		nTaxa = t1->nNodes - t1->nIntNodes;
	
	/* allocate space */
    nLongsNeeded = (nTaxa - 1) / nBitsInALong + 1;
    bitsets = (SafeLong *) SafeCalloc (4*nLongsNeeded*nTaxa+nLongsNeeded, sizeof(SafeLong));
    mask = bitsets + 4*nLongsNeeded*nTaxa;
	
    /* set mask */
    for (i=0; i<nTaxa; i++)
        SetBit(i, mask);

	/* set partition pointers */
	for (i=0; i<t1->nNodes; i++) 
		{
		p = t1->allDownPass[i];
		q = t2->allDownPass[i];
		p->partition = bitsets + i*2*nLongsNeeded;
		q->partition = p->partition + nLongsNeeded;
		}
	
	/* set terminal partitions */
	for (i=0; i<t1->nNodes; i++)
		{
		p = t1->allDownPass[i];
		q = t2->allDownPass[i];
		if (p->right == NULL)
			SetBit (p->index, p->partition);
		if (q->right == NULL)
			SetBit (q->index, q->partition);
		}

    /* set internal partitions */
    for (i=0; i<t1->nIntNodes; i++)
        {
        p = t1->intDownPass[i];
        q = t2->intDownPass[i];
        for (j=0; j<nLongsNeeded; j++)
            {
            p->partition[j] = p->left->partition[j] | p->right->partition[j];
            q->partition[j] = q->left->partition[j] | q->right->partition[j];
            }
        }

	/* check for congruence */
	for (i=0; i<t1->nNodes; i++)
		{
		p = t1->allDownPass[i];
		if (p->anc == NULL && t1->isRooted == YES)
			continue;
		if (t1->isRooted == NO && IsBitSet(t2->root->index,p->partition))
            FlipBits(p->partition,nLongsNeeded, mask);
		for (j=0; j<t2->nNodes; j++)
			{
			q = t2->allDownPass[j];
			for (k=0; k<nLongsNeeded; k++)
				{
				if (p->partition[k] != q->partition[k])
					break;
				}
			if (k == nLongsNeeded && AreDoublesEqual (p->length, q->length, 0.000001) == YES)
				break;
            else if (k == nLongsNeeded)
                {
                free (bitsets);
                return (NO);
                }
			}
		if (j == t2->nNodes)
			{
			free (bitsets);
			return (NO);			
			}
		}
	free (bitsets);
	return (YES);
}





/*----------------------------------------------------------------
|
|	BuildConstraintTree: Build constraint tree. The tree t is
|      needed only to hold information about constraints and
|      included taxa.
|
----------------------------------------------------------------*/
int BuildConstraintTree (Tree *t, PolyTree *pt, char **localTaxonNames)

{

	int				i, j, k, constraintId, nLongsNeeded, nextNode;
	SafeLong		*constraintPartition, *mask;
	PolyNode		*pp, *qq, *rr, *ss, *tt;
	
	pt->isRooted = t->isRooted;

    nLongsNeeded = (numLocalTaxa - 1) / nBitsInALong + 1;
	constraintPartition = (SafeLong *) SafeCalloc (2*nLongsNeeded, sizeof(SafeLong));
	if (!constraintPartition)
		{
		MrBayesPrint ("%s   Problems allocating constraintPartition in BuildConstraintTree", spacer);
		return (ERROR);
		}
	mask = constraintPartition + nLongsNeeded;

	/* calculate mask (needed to take care of unused bits when flipping partitions) */
	for (i=0; i<numLocalTaxa; i++)
		SetBit (i, mask); 

	/* reset all nodes */
	for (i=0; i<2*numLocalTaxa; i++)
		{
		pp = &pt->nodes[i];
		pp->isDated = NO;
		pp->calibration = NULL;
		pp->age = -1.0;
		pp->isLocked = NO;
		pp->lockID = -1;
		pp->index = i;
		}

	/* build a bush */
	pt->root = &pt->nodes[numLocalTaxa];
	for (i=0; i<numLocalTaxa; i++)
		{
		pp = &pt->nodes[i];
		pp->index = i;
		pp->left = NULL;
		if (i == numLocalTaxa - 1)
			pp->sib = NULL;
		else
			pp->sib = &pt->nodes[i+1];
		pp->anc = pt->root;
		}
	pp = pt->root;
	pp->left = &pt->nodes[0];
	pp->anc = pp->sib = NULL;
	pt->nNodes = numLocalTaxa + 1;
	pt->nIntNodes = 1;

	/* make sure the outgroup is the right-most node */
	pt->nodes[localOutGroup].index = numLocalTaxa - 1;
	pt->nodes[numLocalTaxa - 1].index = localOutGroup;

	/* allocate and set partition specifiers in bush */
	GetPolyDownPass(pt);
	AllocatePolyTreePartitions(pt);

	/* set terminal taxon labels */
	for (i=0; i<pt->nNodes; i++)
		{
		pp = pt->allDownPass[i];
		if (pp->index < numLocalTaxa)
			strcpy (pp->label, localTaxonNames[pp->index]);
		}

	/* resolve the bush according to constraints */
	/* for now, satisfy all constraints */
	/* for now, bail out if constraints are not compatible */
	/* Eventually, we might want to be build a parsimony (WAB) or compatibility (WIB) matrix and
	   draw a starting tree from the universe according to the score of the tree. A simple way of accomplishing
	   approximately this is to use sequential addition, with probabilities in each step determined
	   by the parsimony or compatibility score of the different possibilities. */ 
	nextNode = numLocalTaxa + 1;
	for (constraintId=0; constraintId<numDefinedConstraints; constraintId++)
		{
		if (t->constraints[constraintId] == NO || definedConstraintsType[constraintId] != HARD )
			continue;

        /* initialize bits in partition to add; get rid of deleted taxa in the process */
        ClearBits(constraintPartition, nLongsNeeded);
        for (i=j=0; i<numTaxa; i++)
            {
            if (taxaInfo[i].isDeleted == YES)
                continue;
            if (IsBitSet(i, definedConstraint[constraintId]) == YES)
                SetBit(j, constraintPartition);
            j++;
            }
        assert (j == numLocalTaxa);
				
		/* make sure outgroup is outside constrained partition if the tree is unrooted */
		if (t->isRooted == NO && IsBitSet(localOutGroup, constraintPartition))
			FlipBits(constraintPartition, nLongsNeeded, mask);

		/* check that partition should be included */
        k = NumBits(constraintPartition, nLongsNeeded);
        if (k == 0)
			{
			MrBayesPrint ("%s   WARNING: Constraint '%s' refers only to deleted taxa\n", spacer, constraintNames[constraintId]);
			MrBayesPrint ("%s            and will be disregarded\n", spacer);
			t->constraints[constraintId] = NO;
			continue;
            }
        if (k == 1)
			{
			MrBayesPrint ("%s   WARNING: Constraint '%s' refers to a single tip and\n", spacer, constraintNames[constraintId]);
			MrBayesPrint ("%s            will be disregarded\n", spacer);
			t->constraints[constraintId] = NO;
			continue;
            }

		/* check if root in rooted tree (we allow this to enable inference of ancestral states) */
        if (k == numLocalTaxa && t->isRooted == YES)
            {
            pt->root->isLocked = YES;
            pt->root->lockID = constraintId;
            continue;
            }

		/* check if interior root in unrooted tree (we allow this to enable inference of ancestral states) */
        if ((k == numLocalTaxa - 1 || k == numLocalTaxa) && t->isRooted == NO)
            {
            pt->root->isLocked = YES;
            pt->root->lockID = constraintId;
            continue;
            }

		/* find first included terminal */
		k = FirstTaxonInPartition (constraintPartition, nLongsNeeded);
		for (i=0; pt->nodes[i].index != k; i++)
			;
		pp = &pt->nodes[i];

		/* go down until node is not included in constraint */
		do {
			qq = pp;
			pp = pp->anc;		
		} while (IsPartNested(pp->partition, constraintPartition, nLongsNeeded));	

		/* check that the node has not yet been included */
		for (i=0; i<nLongsNeeded; i++)
			{
			if (qq->partition[i] != constraintPartition[i])
				break;
			}
		if (i==nLongsNeeded)
			{
			MrBayesPrint ("%s   WARNING: Constraint '%s' is a duplicate of another constraint\n", spacer, constraintNames[constraintId]);
			MrBayesPrint ("%s            and will be ignored\n", spacer);
			t->constraints[i] = NO;
			continue;
			}

        /* create a new node */
		tt = &pt->nodes[nextNode++];
		tt->anc = pp;
		tt->isLocked = YES;
		tt->lockID = constraintId;
		for (i=0; i<nLongsNeeded; i++)
			tt->partition[i] = constraintPartition[i];
		pt->nIntNodes++;
		pt->nNodes++;

		/* sort descendant nodes in two connected groups: included and excluded */
		/* if there is a descendant that overlaps (incompatible) then return error */
		rr = ss = NULL;
		qq = pp->left;
		do {
			if (IsPartNested(qq->partition, constraintPartition, nLongsNeeded))
				{
				if (ss != NULL)
					ss->sib = qq;
				else
					tt->left = qq;
				ss = qq;
				qq->anc = tt;
				}
			else if (IsPartCompatible(qq->partition, constraintPartition, nLongsNeeded))
				{
				if (rr != NULL)
					rr->sib = qq;
				else
					tt->sib = qq;
				rr = qq;
				}
			else
				{
				free (constraintPartition);
				return (ERROR);
				}
			qq = qq->sib;
			} while (qq != NULL);
		pp->left = tt;
		rr->sib = ss->sib = NULL;
		}
	
	/* relabel interior nodes */
    GetPolyDownPass(pt);
	for (i=0; i<pt->nIntNodes; i++)
		pt->intDownPass[i]->index = i + numLocalTaxa;

	/* exit */
	free (constraintPartition);
	FreePolyTreePartitions(pt);
	return (NO_ERROR);
}





/*----------------------------------------------
|
|   BuildRandomRTopology: Builds a random rooted
|      topology. Will set indices in t->nodes
|      such that they are from 0 to n-1 for the n tips
|      and from n to 2n-2 for the n-1 interior 
|      nodes. Last is root. Does not touch labels
|      of tips.
|
----------------------------------------------*/
int BuildRandomRTopology (Tree *t, SafeLong *seed)
{
	int			i, j, nTips;
	TreeNode	*p, *q, *r;

	nTips = t->nNodes - t->nIntNodes - 1;
	
	for (i=0; i<t->nNodes; i++)
		{
		p = &t->nodes[i];
		p->index = i;
		p->left = p->right = p->anc = NULL;
		}

	/* connect the first two tip nodes */
	q = &t->nodes[0];
	r = &t->nodes[1];
	p = &t->nodes[nTips];
	q->anc = r->anc = p;
	p->left = q;
	p->right = r;
	q = &t->nodes[2*nTips-1];
	p->anc = q;
	q->left = p;

	for (i=2; i<nTips; i++)
		{
		q = &t->nodes[i];
		r = &t->nodes[i-2+nTips+1];
		q->anc = r;
		r->left = q;
		j = (int) (RandomNumber(seed) * (2 * i - 1));
		if (j < i)
			p = &t->nodes[j];
		else
			p = &t->nodes[j-i + nTips];
		r->right = p;
		r->anc = p->anc;
		if (p->anc != NULL)
			{
			if (p->anc->left == p)
				p->anc->left = r;
			else
				p->anc->right = r;
			}
		p->anc = r;
		}

	/* set root and get downpass */
	t->root = &t->nodes[2*nTips-1];
	GetDownPass(t);

	/* relabel interior nodes */
	for (i=0; i<t->nIntNodes; i++)
		t->intDownPass[i]->index = i+nTips;

	return (NO_ERROR);
}





/*----------------------------------------------
|
|   BuildRandomUTopology: Builds a random unrooted
|      topology. Assumes that indices are set
|      in t->nodes from 0 to n-1 for the n tips
|      and from n to 2n-3 for the n-2 interior 
|      nodes. Move the calculation root after
|      this routine to get the right root.
|
----------------------------------------------*/
int BuildRandomUTopology (Tree *t, SafeLong *seed)
{
	int			i, j, nTips;
	TreeNode	*p, *q, *r;

	nTips = t->nNodes - t->nIntNodes;
	
	for (i=0; i<t->nNodes; i++)
		{
		p = &t->nodes[i];
		p->index = i;
		p->left = p->right = p->anc = NULL;
		}
	
	/* connect the first three nodes, assuming 0 is calc root */
	q = &t->nodes[1];
	r = &t->nodes[2];
	p = &t->nodes[nTips];
	q->anc = r->anc = p;
	p->left = q;
	p->right = r;
	q = &t->nodes[0];
	p->anc = q;
	q->left = p;

	for (i=3; i<nTips; i++)
		{
		q = &t->nodes[i];
		r = &t->nodes[i - 3 + nTips + 1];
		q->anc = r;
		r->left = q;
		j = (int) (RandomNumber(seed) * (2 * i - 3));
		if (j < i - 1)
			p = &t->nodes[j+1];
		else
			p = &t->nodes[j+1-i + nTips];
		r->right = p;
		r->anc = p->anc;
		if (p->anc->left == p)
			p->anc->left = r;
		else
			p->anc->right = r;
		p->anc = r;
		}

	t->root = &t->nodes[0];
	
	/* get downpass */
	GetDownPass(t);

	/* relabel interior nodes */
	for (i=0; i<t->nIntNodes; i++)
		t->intDownPass[i]->index = i+nTips;

	return (NO_ERROR);
}





/*----------------------------------------------------------------
|
|	CheckConstraints: Check that tree complies with constraints
|
----------------------------------------------------------------*/
int CheckConstraints (Tree *t)

{

	int				a, i, j, k, nLongsNeeded;
	SafeLong		*constraintPartition, *mask;
	TreeNode		*p=NULL;
	   	
	if (t->checkConstraints == NO)
		return (NO_ERROR);

	/* allocate space */
	nLongsNeeded = (numLocalTaxa - 1) / nBitsInALong + 1;
	constraintPartition = (SafeLong *) SafeCalloc (2*nLongsNeeded, sizeof(SafeLong));
	if (!constraintPartition)
		{
		MrBayesPrint ("%s   Problems allocating constraintPartition in CheckConstraints", spacer);
		return (ERROR);
		}
	mask = constraintPartition + nLongsNeeded;

    /* set mask (needed to reset unused bits when flipping partitions) */
	for (i=0; i<numLocalTaxa; i++) 
	  SetBit (i, mask); 
	
	if (AllocateTreePartitions(t) == ERROR)
		{
		MrBayesPrint ("%s   Problems allocating tree partitions in CheckConstraints", spacer);
		return (ERROR);
		}

	for (a=0; a<numDefinedConstraints; a++)
		{
        if (t->constraints[a] == NO  || definedConstraintsType[a] != HARD)
            continue;

		/* set bits in partition to check */
        ClearBits(constraintPartition, nLongsNeeded);
		for (j=k=0; j<numTaxa; j++)
            {
            if (taxaInfo[j].isDeleted == YES)
                continue;
            if (IsBitSet(j, definedConstraint[a]) == YES)
                SetBit(k, constraintPartition);
            k++;
            }

		/* make sure outgroup is outside constrained partition if unrooted tree */
		if (t->isRooted == NO && IsBitSet(localOutGroup, constraintPartition))
			FlipBits(constraintPartition, nLongsNeeded, mask);

		/* find the locked node */
		for (i=j=0; i<t->nNodes; i++)
			{
			if (t->allDownPass[i]->isLocked == YES && t->allDownPass[i]->lockID == a)
				{
				p = t->allDownPass[i];
				j++;
				}
			}
	
		if (j != 1)
			{
            MrBayesPrint ("%s   Tree has %d locks with id %d identifying constraint '%s'\n", spacer, j, a, constraintNames[a]);
			free (constraintPartition);
			FreeTreePartitions(t);
            return (ERROR);
			}

		/* check that locked node is correct */
		for (i=0; i<nLongsNeeded; i++)
			{
			if (p->partition[i] != constraintPartition[i]) 
				{
				MrBayesPrint ("%s   Lock %d is set for the wrong node [this is a bug]\n", spacer, a);
				free (constraintPartition);
				FreeTreePartitions(t);
                return (ERROR);
				}
			}
		}
	
	FreeTreePartitions (t);
	free (constraintPartition);
	return (NO_ERROR);
}




/*----------------------------------------------------------------
|
|	CheckSetConstraints: Check and set tree constraints
|
----------------------------------------------------------------*/
int CheckSetConstraints (Tree *t)

{

	int				a, i, j, k, nLongsNeeded, foundIt;
	SafeLong		*constraintPartition, *mask;
	TreeNode		*p;
	   	
    if (t->checkConstraints == NO)
	    return (NO_ERROR);

    /* reset all existing locks, if any */
    t->nLocks=0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->isLocked = NO;
        p->lockID = -1;
        if (p->left != NULL)
            {
            p->calibration = NULL;
            p->isDated = NO;
            p->age = -1.0;
            }
        }

    /* allocate space */
	if (AllocateTreePartitions (t) == ERROR)
		{
		MrBayesPrint ("%s   Problems allocating tree bitsets", spacer);
		return ERROR;
		}

	nLongsNeeded = ((numLocalTaxa - 1) / nBitsInALong) + 1;
	constraintPartition = (SafeLong *) SafeCalloc (2*nLongsNeeded, sizeof(SafeLong));
	if (!constraintPartition)
		{
		MrBayesPrint ("%s   Problems allocating constraintPartition", spacer);
		FreeTreePartitions(t);
		return ERROR;
		}
	mask = constraintPartition + nLongsNeeded;

    /* set mask (needed to take care of unused bits when flipping partitions) */
	for (i=0; i<numLocalTaxa; i++)
		SetBit (i, mask);
	
	for (a=0; a<numDefinedConstraints; a++)
		{
		if (modelParams[t->relParts[0]].activeConstraints[a] == NO || definedConstraintsType[a] != HARD )
			continue;

		/* set bits in partition to add */
        ClearBits(constraintPartition, nLongsNeeded);
		for (i=j=0; i<numTaxa; i++)
            {
            if (taxaInfo[i].isDeleted == YES)
                continue;
            if (IsBitSet(i, definedConstraint[a]) == YES)
                SetBit(j, constraintPartition);
            j++;
            }

        k = NumBits(constraintPartition, nLongsNeeded);
        if (k == 0 || k == 1)
            continue;

		/* make sure outgroup is outside constrained partition (marked 0) */
		if (t->isRooted == NO && IsBitSet(localOutGroup, constraintPartition) == YES)
			FlipBits(constraintPartition, nLongsNeeded, mask);

		/* find the node that should be locked */
		foundIt = NO;
		for (i=0; i<t->nIntNodes; i++)
			{
			p = t->intDownPass[i];
			for (j=0; j<nLongsNeeded; j++)
				{
				if (p->partition[j] != constraintPartition[j])
					break;
				}

			if (j == nLongsNeeded)
				{
				foundIt = YES;
                p->isLocked = YES;
				p->lockID = a;
				if (nodeCalibration[a].prior != unconstrained)
					{
					p->isDated = YES;
					p->calibration = &nodeCalibration[a];
					}
                t->nLocks++;
				break;
				}				
			}
	
		if (foundIt == NO)
			{
            MrBayesPrint ("%s   Tree breaks constraint '%s'\n", spacer, constraintNames[a]);
			FreeTreePartitions (t);
			free (constraintPartition);
			return (ERROR);
			}
		}
	
	/* exit */
	FreeTreePartitions(t);
	free (constraintPartition);
	return (NO_ERROR);
}





/*-----------------------------------------------------------------------
|
|   ColorClusters: Recursive function to color the clusters in a tree by
|      assigning numbers to them in their variable x
|
------------------------------------------------------------------------*/
void ColorClusters (TreeNode *p, int *index)
{
    if (p!=NULL)
        {
        if (p->isLocked == YES || p->anc == NULL || p->anc->anc == NULL)
            p->x = (++(*index));
        else
            p->x = p->anc->x;
        ColorClusters(p->left, index);
        ColorClusters(p->right, index);
        }
}





/* CopyPolyNodes: Copies everything except pointers and memoryIndex */
void CopyPolyNodes (PolyNode *p, PolyNode *q, int nLongsNeeded)
{
	p->index                  = q->index; 
	p->mark                   = q->mark;
	p->length                 = q->length;
	p->x                      = q->x;
	p->y                      = q->y;
	p->isDated				  = q->isDated;
	p->calibration			  = q->calibration;
	p->age					  = q->age;
	p->isLocked				  = q->isLocked;
	p->lockID				  = q->lockID;
	strcpy (p->label, q->label);
    if( nLongsNeeded!=0 )
        {
        assert(p->partition);
        assert(q->partition);
        memcpy(p->partition,q->partition, nLongsNeeded*sizeof(SafeLong));
        }
	p->support                = q->support;
	p->f                      = q->f;
}





void CopySubtreeToTree (Tree *subtree, Tree *t)

{
	
        int			i, /*j,*/ k;
	TreeNode	*p, *q=NULL, *r;

	for (i=/*j=*/0; i<subtree->nNodes - 1; i++)
		{
		p = subtree->allDownPass[i];

		for (k=0; k<t->nNodes; k++)
			{
			q = t->allDownPass[k];
			if (q->index == p->index)
				break;
			}
		q->length = p->length;
		q->marked = YES;
		if (p->left != NULL && p->right != NULL)
			{
			for (k=0; k<t->nNodes; k++)
				{
				r = t->allDownPass[k];
				if (r->index == p->left->index)
					{
					q->left = r;
					r->anc = q;
					}
				else if (r->index == p->right->index)
					{
					q->right = r;
					r->anc = q;
					}
				}
			}
		}

	p = subtree->root;

	for (k=0; k<t->nNodes; k++)
		{
		q = t->allDownPass[k];
		if (q->index == p->index)
			break;
		}

	if (q->left->marked == YES)
		{
		for (k=0; k<t->nIntNodes; k++)
			{
			r = t->intDownPass[k];
			if (r->index == p->left->index)
				{
				q->left = r;
				r->anc = q;
				}
			}
		}
	else if (q->right->marked == YES)
		{
		for (k=0; k<t->nIntNodes; k++)
			{
			r = t->intDownPass[k];
			if (r->index == p->left->index)
				{
				q->right = r;
				r->anc = q;
				}
			}
		}
}





/*-----------------------------------------------------------------
|
|	CopyToPolyTreeFromPolyTree: copies second tree to first tree
|
-----------------------------------------------------------------*/
int CopyToPolyTreeFromPolyTree (PolyTree *to, PolyTree *from)
{
	int			i, j, k, nLongsNeeded;
	PolyNode	*p, *q;

    /* check we have enough memory */
    assert (to->memNodes >= from->nNodes);
    if( from->bitsets==NULL || to->bitsets==NULL )
        {
        nLongsNeeded=0;
        }
    else
        {
        assert (to->memNodes >= from->memNodes);/*Otherwise partition length woould not be long enough for nodes in "to" */
        nLongsNeeded = (from->memNodes/2 - 1) / nBitsInALong + 1;
        }


	/* copy nodes */
	for (i=0; i<from->nNodes; i++)
		{
		/* copy pointers */
		p  = from->nodes + i;
		q  = to->nodes + i;

		if (p->anc != NULL)
			q->anc = to->nodes + p->anc->memoryIndex;
		else
            {
			q->anc = NULL;
            to->root = q;
            }

		if (p->left != NULL)
			q->left = to->nodes + p->left->memoryIndex;
		else
			q->left = NULL;

		if (p->sib != NULL)
 			q->sib = to->nodes + p->sib->memoryIndex;
		else
			q->sib = NULL;

		/* Copy everything else except memoryIndex */
        CopyPolyNodes (q, p, nLongsNeeded);
		}

    /* fill node arrays */
	/* copy tree properties */
    to->nNodes = from->nNodes;
	to->nIntNodes = from->nIntNodes;
    to->isRooted = from->isRooted;
    to->isClock = from->isClock;
    to->isRelaxed = from->isRelaxed;
    to->clockRate = from->clockRate;
	strcpy (to->name, from->name);
  
    GetPolyDownPass (to);

	/* copy partitions */
    if (from->bitsets)
		{
		if (!to->bitsets)
            AllocatePolyTreePartitions(to);
        else
            ResetPolyTreePartitions(to);
		}

	/* copy relaxed clock parameters */
    FreePolyTreeRelClockParams (to);
	
    if (from->nBSets + from->nESets > 0)
        AllocatePolyTreeRelClockParams (to, from->nBSets, from->nESets);

	for (i=0; i<to->nBSets; i++)
		{
		to->bSetName[i] = (char *) SafeCalloc (strlen(from->bSetName[i])+2, sizeof(char));
		strcpy (to->bSetName[i], from->bSetName[i]);
        for (j=0; j<from->nNodes; j++)
            to->effectiveBrLen[i][j] = from->effectiveBrLen[i][j];
		}
	
	for (i=0; i<to->nESets; i++)
		{
		to->eSetName[i] = (char *) SafeCalloc (strlen(from->eSetName[i])+2, sizeof(char));
		strcpy (to->eSetName[i], from->eSetName[i]);
        for (j=0; j<from->nNodes; j++)
            {
            to->nEvents[i][j] = from->nEvents[i][j];
            if (to->nEvents[i][j] > 0)
                {
                to->position[i][j] = (MrBFlt *) SafeCalloc (to->nEvents[i][j], sizeof (MrBFlt));
                to->rateMult[i][j] = (MrBFlt *) SafeCalloc (to->nEvents[i][j], sizeof (MrBFlt));
                for (k=0; k<to->nEvents[i][j]; k++)
                    {
                    to->position[i][j][k] = from->position[i][j][k];
                    to->rateMult[i][j][k] = from->rateMult[i][j][k];
                    }
                }
            }
        }
	
    /* copy population size parameters */
    FreePolyTreePopSizeParams(to);
    to->popSizeSet = from->popSizeSet;
    if (to->popSizeSet == YES)
        {
        to->popSize = (MrBFlt *) SafeCalloc (to->nNodes-1, sizeof(MrBFlt));
        for (i=0; i<to->nNodes-1; i++)
            to->popSize[i] = from->popSize[i];
        to->popSizeSetName = (char *) SafeCalloc (strlen(from->popSizeSetName) + 1, sizeof(char));
        strcpy (to->popSizeSetName, from->popSizeSetName);
        }

    return (NO_ERROR);

}





/*-----------------------------------------------------------------
|
|	CopyToSpeciesTreeFromPolyTree: copies second tree (polytomous) to
|		first tree, which is a species tree. The species tree needs to
|       be allocated enough space first to hold the resulting tree.
|
-----------------------------------------------------------------*/
int CopyToSpeciesTreeFromPolyTree (Tree *to, PolyTree *from)

{

	int			i;
	PolyNode	*p;
	TreeNode	*q, *q1;
#if !defined (NDEBUG)
    int         j;
#endif

    /* make sure assumptions are correct */
    assert (from->isRooted == YES);
    assert (from->isClock == YES);
    assert (from->nNodes - from->nIntNodes == numSpecies);
    assert (to->memNodes == 2*numSpecies);
    assert (to->nIntNodes == from->nIntNodes);
    assert (to->nNodes == from->nNodes + 1);

    /* make sure indices are set correctly for from nodes */
#if !defined (NDEBUG)
    for (i=0; i<from->nNodes; i++)
        {
        for (j=0; j<from->nNodes; j++)
            {
            p = from->allDownPass[j];
            if (p->index == i)
                break;
            }
        assert (j != from->nNodes);
        assert (!(p->left == NULL && p->index >= numSpecies));
        }
#endif

    /* copy nodes */
	for (i=0; i<from->nNodes; i++)
		{
		/* copy pointers */
		p  = from->allDownPass[i];
		q  = to->nodes + p->index;

		if (p->anc != NULL)
			q->anc = to->nodes + p->anc->index;
		else
			q->anc = NULL;

		if (p->left != NULL)	
			q->left = to->nodes + p->left->index;
		else
			q->left = NULL;

		if (p->left != NULL)
			q->right = to->nodes + p->left->sib->index;
		else
			q->right = NULL;

        q->nodeDepth              = p->depth;
        q->age					  = p->age;
		q->length                 = p->length;
        q->index                  = p->index;
        if (q->index < numSpecies)
            q->label = speciesNameSets[speciespartitionNum].names[q->index];
        else
            q->label = noLabel;
		}

	/* fix root */
	p = from->root;
	q = to->nodes + p->index;
	q1 = to->nodes + from->nNodes;      /* get the 'extra' root node that polytomous trees do not use */
	q->anc = q1;
    q1->index = from->nNodes;
	q1->left = q;
	q1->right = q1->anc = NULL;
	q1->isLocked = NO;
	q1->lockID = -1;
	q1->isDated = NO;
	q1->calibration = NULL;
	q1->age = -1.0;
	to->root = q1;

	/* get downpass */
	GetDownPass (to);
	
    /* a user tree might not come with node depths set */
    if (to->root->left->nodeDepth == 0.0)
        SetNodeDepths(to);

    /* set partitions */
    if (to->bitsets)
        ResetTreePartitions(to);

    return (NO_ERROR);		
}





/*-----------------------------------------------------------------
|
|	CopyToTreeFromPolyTree: copies second tree (polytomous) to first
|		tree (used to initialize constrained starting trees, e.g.).
|		An unrooted source tree will be rooted on outgroup
|		An unrooted source tree that needs to be copied to
|       a rooted target tree will be randomly rooted on a node below
|       all defined constraints. The to tree needs to be allocated
|       enough space first to hold the resulting tree.
|
-----------------------------------------------------------------*/
int CopyToTreeFromPolyTree (Tree *to, PolyTree *from)

{

	int			i, j, nNodesNeeded;
	PolyNode	*p;
	TreeNode	*q, *q1;

    /* refuse to arbitrarily root an input tree */
    assert (!(from->isRooted == NO && to->isRooted == YES));
    if ( (from->isRooted == NO) && (to->isRooted == YES) )
        {
        MrBayesPrint ("%s   Fail to copy trees due to difference in rootedness of source and destination. \n", spacer);
        return (ERROR);
        }
    /* calculate space needed */
    if (from->isRooted == YES && to->isRooted == YES)
		nNodesNeeded    = from->nNodes + 1;
	else if (from->isRooted == YES && to->isRooted == NO)
        nNodesNeeded    = from->nNodes;
    else /* if (from->isRooted = NO && to->isRooted == NO) */
		nNodesNeeded    = from->nNodes;

    /* make sure assumptions are in order */
    assert (to->memNodes >= nNodesNeeded);
    assert (localOutGroup >= 0 && localOutGroup < numLocalTaxa);
    assert (numLocalTaxa == from->nNodes - from->nIntNodes);
    assert (!(from->isRooted == YES && from->nNodes != 2*from->nIntNodes + 1));
    assert (!(from->isRooted == NO  && from->nNodes != 2*from->nIntNodes + 2));
    
    /* make sure indices are set correctly for from nodes */
    for (i=0; i<from->nNodes; i++)
        {
        for (j=0; j<from->nNodes; j++)
            {
            p = from->allDownPass[j];
            if (p->index == i)
                break;
            }
        assert (j != from->nNodes);
        assert (!(p->left == NULL && p->index >= numLocalTaxa));
        }
                
	/* deal with root */
    if (to->isRooted == NO && from->isRooted == YES)
        Deroot(from);

    /* make sure calculation root is set correctly */
    if (to->isRooted == NO && MovePolyCalculationRoot (from, localOutGroup) == ERROR)
        return ERROR;

    /* copy nodes */
	for (i=0; i<from->nNodes; i++)
		{
		/* copy pointers */
		p  = from->allDownPass[i];
		q  = to->nodes + p->index;

		if (p->anc != NULL)
			q->anc = to->nodes + p->anc->index;
		else
			q->anc = NULL;

		if (p->left != NULL)	
			q->left = to->nodes + p->left->index;
		else
			q->left = NULL;

		if (p->left != NULL)
			q->right = to->nodes + p->left->sib->index;
		else
			q->right = NULL;

		q->isLocked				  = p->isLocked;
		q->lockID				  = p->lockID;
		q->isDated				  = p->isDated;
		q->calibration			  = p->calibration;
		q->age					  = p->age;
        q->nodeDepth              = p->depth;
		q->length                 = p->length;
        q->index                  = p->index;
        if (q->index < numLocalTaxa)
            q->label = localTaxonNames[q->index];
        else
            q->label = noLabel;
		}

	/* fix root */
	if (to->isRooted == NO)
		{
		p = from->root;
		q = to->nodes + p->index;
		q->anc = to->root = to->nodes + p->left->sib->sib->index;
		q->length = to->root->length;
		to->root->length = 0.0;
		to->root->left = q;
		to->root->right = to->root->anc = NULL;
		}
	else
		{
		p = from->root;
		q = to->nodes + p->index;
		q1 = to->nodes + from->nNodes;      /* get the 'extra' root node that polytomous trees do not use */
		q->anc = q1;
        q1->index = from->nNodes;
		q1->left = q;
		q1->right = q1->anc = NULL;
		q1->isLocked = NO;
		q1->lockID = -1;
		q1->isDated = NO;
		q1->calibration = NULL;
		q1->age = -1.0;
		to->root = q1;
		}

	/* get downpass */
	GetDownPass (to);
	
	/* set node depths */
    if (to->isRooted == YES && to->root->left->nodeDepth == 0.0)
	    SetNodeDepths(to);

    /* set partitions */
    if (to->bitsets)
        ResetTreePartitions(to);

    /* relaxed clock parameters are not stored in binary trees but in separate parameters */

    return (NO_ERROR);		
}





/*-----------------------------------------------------------------
|
|	CopyToTreeFromTree: copies second tree to first tree
|		(used to initialize brlen sets for same topology)
|       Note: partition information of nodes are not copied if
|       either "from" or "to" tree does not have bitsets allocated
|
-----------------------------------------------------------------*/
int CopyToTreeFromTree (Tree *to, Tree *from)

{

	int			i, numTaxa, nLongsNeeded;
	TreeNode	*p, *q;

    numTaxa = from->nNodes - from->nIntNodes - (from->isRooted == YES ? 1 : 0);
    nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;
    if( from->bitsets==NULL || to->bitsets==NULL )
        nLongsNeeded=0;

	/* check that there is enough memory */
    assert (to->memNodes >= from->nNodes);
    
    /* copy nodes (use index of p as memoryIndex for q) */
	for (i=0; i<from->nNodes; i++)
		{
		/* copy pointers */
		p  = from->nodes + i;
		q  = to->nodes + p->index;

		if (p->anc != NULL)
			q->anc = to->nodes + p->anc->index;
		else
			{
			q->anc = NULL;
			to->root = q;
			}

		if (p->left != NULL)
			q->left = to->nodes + p->left->index;
		else
			q->left = NULL;

		if (p->right != NULL)
 			q->right = to->nodes + p->right->index;
		else
			q->right = NULL;

		CopyTreeNodes (q, p, nLongsNeeded);
		}

	/* create new node arrays */
    to->nNodes = from->nNodes;
	to->nIntNodes = from->nIntNodes;
    GetDownPass(to);

	/* copy tree properties (these should be constant most of them) */
	strcpy (to->name, from->name);
    to->isRooted = from->isRooted;
    to->isClock = from->isClock;
    to->isCalibrated = from->isCalibrated;
    to->checkConstraints = from->checkConstraints;
    to->nConstraints = from->nConstraints;
    to->constraints = from->constraints;
    to->nLocks = from->nLocks;
    to->nRelParts = from->nRelParts;
    to->relParts = from->relParts;

	/* copy partitions */
    if (from->bitsets)
		{
		if (!to->bitsets)
            AllocateTreePartitions(to);
        else
            ResetTreePartitions(to);
		}

	return (NO_ERROR);

}





/* Copeing node q to node p. nLongsNeeded is the length of partition that the node represents. */
void CopyTreeNodes (TreeNode *p, TreeNode *q, int nLongsNeeded)
{

    /* copies everything except pointers and memoryIndex */
	p->index                  = q->index;
	p->scalerNode			  = q->scalerNode;			
	p->upDateCl               = q->upDateCl;
	p->upDateTi				  = q->upDateTi;
	p->marked                 = q->marked;
	p->length                 = q->length;
	p->nodeDepth              = q->nodeDepth;
	p->x                      = q->x;
	p->y                      = q->y;
	p->isDated				  = q->isDated;
	p->calibration			  = q->calibration;
	p->age					  = q->age;
	p->isLocked				  = q->isLocked;
	p->lockID				  = q->lockID;
	p->d					  = q->d;
	//p->partition			  = q->partition;//the content should be copied not the pointers. Otherwise we are risking to have segmentation faults
    if( nLongsNeeded!=0 )
        {
        assert(p->partition);
        assert(q->partition);
        memcpy(p->partition,q->partition, nLongsNeeded*sizeof(SafeLong));
        }
    p->label                  = q->label;
}





void CopyTreeToSubtree (Tree *t, Tree *subtree)

{
	
	int			i, j, k;
	TreeNode	*p, *q, *r;

	for (i=j=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->marked == NO)
			continue;

		q = &subtree->nodes[j++];
		q->index = p->index;
		q->length = p->length;
		if (p->left == NULL || p->left->marked == NO)
			q->left = q->right = NULL;
		else
			{
			for (k=0; k<j-1; k++)
				{
				r = &subtree->nodes[k];
				if (r->index == p->left->index)
					{
					q->left = r;
					r->anc = q;
					}
				else if (r->index == p->right->index)
					{
					q->right = r;
					r->anc = q;
					}
				}
			}
		
		if (p->anc->marked == NO)
			{
			r = &subtree->nodes[j++];
			subtree->root = r;
			r->anc = r->right = NULL;
			r->left = q;
			q->anc = r;
			r->length = 0.0;
			r->index = p->anc->index;
			}

		}

	GetDownPass (subtree);

	subtree->isRooted = t->isRooted;
	subtree->nRelParts = t->nRelParts;
	subtree->relParts = t->relParts;
}





/* DatedNodeDepths: Recursive function to get node depths */
void DatedNodeDepths (TreeNode *p, MrBFlt *nodeDepths, int *index)
{
    if (p != NULL)
        {
        if (p->left == NULL || p->isDated == YES)
            nodeDepths[(*index)++] = p->nodeDepth;
        else
            {
            DatedNodeDepths (p->left,  nodeDepths, index);
            DatedNodeDepths (p->right, nodeDepths, index);
            }
        }
}





/* DatedNodes: Recursive function to get dated tips or interior nodes */
void DatedNodes (TreeNode *p, TreeNode **datedNodes, int *index)
{
    if (p != NULL)
        {
        if (p->left != NULL && p->isDated == NO)
            {
            DatedNodes (p->left,  datedNodes, index);
            DatedNodes (p->right, datedNodes, index);
            }
        datedNodes[(*index)++] = p;
        }
}





/* Deroot: Deroot a rooted polytomous tree with branch lengths */
int Deroot (PolyTree *pt)
{
    PolyNode    *p, *q, *r, tempNode;
    int         i;

    p = pt->root;

    if (p->left->sib->sib != NULL)
        return (ERROR);      /* tree is not rooted or it is polytomous */

    if (p != &pt->nodes[pt->nNodes-1])
        {
        q = &pt->nodes[pt->nNodes-1];
        /* now swap content of p and q including pointers */
        tempNode = *q;
        *q = *p;
        *p = tempNode;
        /* swap back memoryindex */
        p->memoryIndex = q->memoryIndex;
        q->memoryIndex = tempNode.memoryIndex;
        /* all pointers to q should be pointers to p */
        for (i=0; i<pt->nNodes; i++)
            {
            r = &pt->nodes[i];
            if (r->left == q)
                r->left = p;
            if (r->anc == q)
                r->anc = p;
            if (r->sib == q)
                r->sib = p;
            }
        /* all pointers to p should be pointers to q; all these are anc pointers from the descendants of the root */
        pt->root = q;
        for (r=q->left; r!=NULL; r=r->sib)
            r->anc = q;
        /* finally set p to the new root */
        p = pt->root;
        }

    /* make sure the left of the old root is interior and can be used as new root */
    if (p->left->left == NULL)
        {
        q = p->left;
        r = q->sib;
        p->left = r;
        r->sib = q;
        q->sib = NULL;
        }
    
    pt->root = p->left;
    pt->root->left->sib->sib = p->left->sib;
    p->left->sib->length += pt->root->length;
    pt->root->length = 0.0;
    pt->root->sib = NULL;
    pt->root->anc = NULL;

    pt->nNodes--;
    pt->nIntNodes--;

    GetPolyDownPass(pt);

    return (NO_ERROR);
}





/* EraseTreeList: Erase all trees in treeList */
void EraseTreeList (TreeList *treeList)
{
	TreeListElement *listElement;
	TreeListElement *previous;

	listElement = treeList->first;
	if (listElement != NULL)
		do 
			{
			free (listElement->order);
			previous = listElement;
			listElement = listElement->next;
			free (previous);
			} 
		while (listElement != NULL);

	treeList->first = treeList->last = NULL;
}





/* FreePolyTree: Free memory space for a polytomous tree (unrooted or rooted) */
void FreePolyTree (PolyTree *pt)
{	
	if (pt != NULL)
		{
        FreePolyTreePartitions(pt);
        FreePolyTreeRelClockParams(pt);
        FreePolyTreePopSizeParams(pt);
		free (pt->allDownPass);
		free (pt->nodes);
        free (pt);
		}
}





/* FreePolyTreePartitions: Free memory space for polytomous tree partitions */
void FreePolyTreePartitions (PolyTree *pt)
{
    int i;
	if (pt != NULL && pt->bitsets != NULL)
		{
        for (i=0; i<pt->memNodes; i++)
            pt->nodes[i].partition = NULL;
		free (pt->bitsets);
		pt->bitsets = NULL;
		}
}





/* FreePolyTreePopSizeParams: Free population size set parameters of polytree */
void FreePolyTreePopSizeParams (PolyTree *pt)
{
    if (pt->popSizeSet == YES)
        {
        free (pt->popSize);
        free (pt->popSizeSetName);
        }
    pt->popSizeSet = NO;
    pt->popSize = NULL;
    pt->popSizeSetName = NULL;
}





/* FreePolyTreeRelClockParams: Free relaxed clock parameters of polytree */
void FreePolyTreeRelClockParams (PolyTree *pt)
{
    int i, j;

    /* free breakpoint clock parameters */
    for (i=0; i<pt->nESets; i++)
        {
        for (j=0; j<pt->memNodes; j++)
            {
            if (pt->nEvents[i][j] > 0)
                {
                free (pt->position[i][j]);
                free (pt->rateMult[i][j]);
                }
            }
		free (pt->eSetName[i]);
        free (pt->nEvents[i]);
        free (pt->position[i]);
        free (pt->rateMult[i]);
        }
    free (pt->nEvents);
    free (pt->position);
    free (pt->rateMult);
    free (pt->eSetName);
	pt->nESets = 0;
    pt->nEvents = NULL;
    pt->position = NULL;
    pt->rateMult = NULL;
    pt->eSetName = NULL;

    /* free branch clock parameters */
	for (i=0; i<pt->nBSets; i++)
        {
		free (pt->bSetName[i]);
        free (pt->effectiveBrLen[i]);
        }
    free (pt->effectiveBrLen);
    free (pt->bSetName);
    pt->nBSets = 0;
    pt->effectiveBrLen = NULL;
    pt->bSetName = NULL;
}





/* FreeTree: Free memory space for a tree (unrooted or rooted) */
void FreeTree (Tree *t)
{
	if (t != NULL)
		{
		free (t->bitsets);
		free (t->flags);
		free (t->allDownPass);
		free (t->nodes);
		free (t);
		}
}





/* FreeTreePartitions: Free memory space for tree partitions */
void FreeTreePartitions (Tree *t)
{
    int     i;

	if (t != NULL && t->bitsets != NULL)
		{
		free (t->bitsets);
		t->bitsets = NULL;
        for (i=0; i<t->memNodes; i++)
            t->nodes[i].partition = NULL;
		}
}





/*-------------------------------------------------------------------------------------------
|
|   GetDatedNodeDepths: Get an array containing the node depths of the dated tips,
|       internal or external, plus dated root
|
---------------------------------------------------------------------------------------------*/
void GetDatedNodeDepths (TreeNode *p, MrBFlt *nodeDepths)
{
    int index = 0;
    
    assert (p != NULL);

    nodeDepths[index++] = p->nodeDepth;     /* include root node depth */
    if (p->left != NULL)
        {
        DatedNodeDepths (p->left, nodeDepths, &index);
        DatedNodeDepths (p->right, nodeDepths, &index);
        }
}





/*-------------------------------------------------------------------------------------------
|
|   GetDatedNodes: Get an array containing the dated tips,
|       internal or external, and all interior nodes in the same subtree
|
---------------------------------------------------------------------------------------------*/
void GetDatedNodes (TreeNode *p, TreeNode **datedNodes)
{
    int     index = 0;
    
    assert (p != NULL);

    if (p->left!= NULL)
        {
        DatedNodes (p->left,  datedNodes, &index);
        DatedNodes (p->right, datedNodes, &index);
        }
}





/* get down pass for tree t (wrapper function) */
void GetDownPass (Tree *t)

{

	int i, j;

    i = j = 0;
	GetNodeDownPass (t, t->root, &i, &j);
		
}





/* get the actual down pass sequences */
void GetNodeDownPass (Tree *t, TreeNode *p, int *i, int *j)

{
	
	if (p != NULL )
		{
		GetNodeDownPass (t, p->left,  i, j);
		GetNodeDownPass (t, p->right, i, j);
		if (p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			t->intDownPass[(*i)++] = p;
			t->allDownPass[(*j)++] = p;
			}
		else if (p->left == NULL && p->right == NULL && p->anc != NULL)
			{
			t->allDownPass[(*j)++] = p;
			}
		else if (p->left != NULL && p->right == NULL && p->anc == NULL)
			{
			t->allDownPass[(*j)++] = p;
			}
		}
		
}





/* GetPolyAges: Get PolyTree node ages */
void GetPolyAges (PolyTree *t)
{
    int         i;
    PolyNode    *p;

    GetPolyDepths (t); /* just to make sure... */
    
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->age = p->depth / t->clockRate;
        }
}





/* GetPolyDepths: Get PolyTree node depths */
void GetPolyDepths (PolyTree *t)
{
    int         i;
    MrBFlt      maxDepth;
    PolyNode    *p;

    maxDepth = t->root->depth = 0.0;

    for (i=t->nNodes-2; i>=0; i--)
        {
        p = t->allDownPass[i];
        p->depth = p->anc->depth + p->length;
        if (p->depth > maxDepth)
            maxDepth = p->depth;
        }

    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->depth = maxDepth - p->depth;
        }
}





/* get down pass for polytomous tree t (wrapper function) */
void GetPolyDownPass (PolyTree *t)

{

	int i, j;

	i = j = 0;
	GetPolyNodeDownPass (t, t->root, &i, &j);
    assert( t->nIntNodes==j );
		
}





/* get the actual down pass sequences for a polytomous tree */
void GetPolyNodeDownPass (PolyTree *t, PolyNode *p, int *i, int *j)

{
	
	PolyNode	*q;
	
	if (p->left != NULL)
		{
		for (q=p->left; q!=NULL; q=q->sib)
			GetPolyNodeDownPass(t, q, i, j);
		}

	t->allDownPass[(*i)++] = p;
	if (p->left != NULL )
		t->intDownPass[(*j)++] = p;

}





/* GetFromTreeList: Get first tree from a tree list */
int GetFromTreeList (TreeList *treeList, Tree *tree)

{
	TreeListElement	*listElement;

	if (treeList->first == NULL)
		{
		MrBayesPrint ("%s   Tree list empty\n", spacer);
		return (ERROR);
		}
	if (tree->isRooted == YES)
		RetrieveRTopology (tree, treeList->first->order);
	else
        {
        RetrieveUTopology (tree, treeList->first->order);
        if (localOutGroup != 0)
            MoveCalculationRoot (tree, localOutGroup);
        }

	listElement = treeList->first;
	treeList->first = listElement->next;

	free (listElement->order);
	free (listElement);

	return (NO_ERROR);
}





/*------------------------------------------------------------------
|
|   InitBrlens: This routine will set all branch lengths of a
|      nonclock tree to the value given by 'v'.
|
------------------------------------------------------------------*/
int InitBrlens (Tree *t, MrBFlt v)

{
	int			i;
	TreeNode	*p;

	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->anc != NULL && !(t->isRooted == YES && p->anc->anc == NULL))
			p->length = v;
		else
			p->length = 0.0;
		}

	return (NO_ERROR);
}





/* 
@param levUp		is the number of edges between the "node" and the most resent calibrated predecessor +1,
					so for the calibrated nodes it should be 1
@param calibrUp		is the age of the most resent calibrated predecessor
@return				age of the node	
*/
MrBFlt SetNodeCalibratedAge(TreeNode *node, unsigned levUp, MrBFlt calibrUp )
{
	MrBFlt r,l;

	if( node->age != -1.0 )
		{
		if(node->right != NULL)
			SetNodeCalibratedAge( node->right, 2, node->age );
		if(node->left != NULL)
			SetNodeCalibratedAge( node->left, 2, node->age );
		return node->age;
		}

	r = SetNodeCalibratedAge( node->right, levUp+1, calibrUp );
	l = SetNodeCalibratedAge( node->left, levUp+1, calibrUp );

	if( r>l )
		{
		assert( calibrUp - r > 0.0 );
		return node->age = r + ( calibrUp - r )/levUp;
		}
	else
		{
		assert( calibrUp - l> 0.0 );
		return node->age = l + ( calibrUp - l )/levUp;
		}

}





/*-------------------------------------------------------------------
|
|	InitCalibratedBrlens: This routine will build a clock tree
|		consistent with calibration constraints on terminal
|		taxa and/or constrained interior nodes. At least one
|		node should be calibrated.
|    	If not possible to build such a tree, ERROR
|       is returned.
|
--------------------------------------------------------------------*/
int InitCalibratedBrlens (Tree *t, MrBFlt clockRate, SafeLong *seed)

{

	int				i, treeAgeFixed;
	TreeNode		*p;
    Model           *mp;
    MrBFlt          treeAge;

#if 0
    printf ("Before initializing calibrated brlens\n");
    ShowNodes(t->root, 0, YES);
#endif
    
    if (t->isRooted == NO)
		{
		MrBayesPrint ("%s   Tree is unrooted\n", spacer);
		return (ERROR);
		}

    /* Check whether root is fixed indirectly */
    mp = &modelParams[t->relParts[0]];
    if (!strcmp(mp->clockPr, "Uniform") && !strcmp(mp->treeAgePr,"Fixed"))
        {
        treeAgeFixed = YES;
        treeAge = mp->treeAgeFix;
        }
    else
        {
        treeAgeFixed = NO;
        treeAge = -1.0;
        }
	
	/* date all nodes from top to bottom with min. age as nodeDepth*/
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->anc != NULL)
			{
			if (p->left == NULL && p->right == NULL)
				{
				if (p->isDated == NO)
					{
					p->nodeDepth = 0.0;
					p->age = 0.0;
					}
				else
					{
					if (p->calibration->prior == fixed)
						p->nodeDepth = p->age = p->calibration->age;
					else if (p->calibration->prior == uniform)
						p->nodeDepth = p->age = p->calibration->min;
					else /* if (p->calibration->prior == offsetExponential) */
						p->nodeDepth = p->age = p->calibration->offset;
					}
				}
			else
				{
				if (p->left->nodeDepth > p->right->nodeDepth)
					p->nodeDepth = p->left->nodeDepth;
				else
					p->nodeDepth = p->right->nodeDepth;
				if (p->isDated == YES || (p->anc->anc == NULL && treeAgeFixed))
					{
                    if (p->isDated == NO)
                        {
                        if (treeAge <= p->nodeDepth)
                            {
						    MrBayesPrint ("%s   Calibration inconsistency for root node\n", spacer);
                            return (ERROR);
                            }
                        else
                            p->age = p->nodeDepth = treeAge;
                        }
					else if (p->calibration->prior == fixed)
						{
						if (p->calibration->age <= p->nodeDepth)
							{
                            if (p->anc->anc == NULL)
							    MrBayesPrint ("%s   Calibration inconsistency for root node\n", spacer);
                            else
                                MrBayesPrint ("%s   Calibration inconsistency for node '%s'\n", spacer, constraintNames[p->lockID]);
							return (ERROR);
							}
						else
							p->age = p->nodeDepth = p->calibration->age;
						}
					else if (p->calibration->prior == uniform)
						{
						if (p->calibration->max <= p->nodeDepth)
							{
							MrBayesPrint ("%s   Calibration inconsistency for node '%s'\n", spacer, constraintNames[p->lockID]);
							return (ERROR);
							}
						else
							{
							if (p->calibration->min < p->nodeDepth)
								p->age = p->nodeDepth;
							else
								p->age = p->nodeDepth = p->calibration->min;
							}
						}
					else /* if (p->calibration.prior == offsetExponential) */
						{
						if (p->calibration->offset < p->nodeDepth)
							p->age = p->nodeDepth;
						else
							p->age = p->nodeDepth = p->calibration->offset;
						}
					}
				else
					p->age = -1.0;
				}
			}
		}
	

	/* try to make root node deeper than minimum age */
    p = t->root->left;
    if (p->calibration == NULL && treeAgeFixed == NO)
        {
		if( p->nodeDepth==0.0 ) p->nodeDepth = 1.0;
        p->age = p->nodeDepth *= 1.5;
        }
    else if (treeAgeFixed == YES || p->calibration->prior == fixed)
        {
        /* can't do much ... */
        }
    else if (p->calibration->prior == uniform)
        {
        if (p->nodeDepth * 1.5 < p->calibration->max)
            p->nodeDepth = p->age = 1.5 * p->nodeDepth;
        else
            p->nodeDepth = p->age = p->calibration->max;
        }
    else /* if (t->root->calibration->prior == offsetExponential */
        {
		assert( p->calibration->prior == offsetExponential );
        /* Make random draw */
        p->age = p->calibration->offset - log (RandomNumber(seed)) / p->calibration->lambda;
        if (p->age > 1.5 * p->nodeDepth)
            p->nodeDepth = p->age;
        else
            p->age = p->nodeDepth = 1.5 * p->nodeDepth;
        }

	SetNodeCalibratedAge( p, 1, p->age );


    /* Setup node depths */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		p->nodeDepth = p->age * clockRate;
        assert (!(p->left == NULL && p->calibration == NULL && p->nodeDepth != 0.0));
		}

	/* calculate branch lengths */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->anc != NULL)
			{
			if (p->anc->anc != NULL)
				{
				p->length = p->anc->nodeDepth - p->nodeDepth;
				assert( p->length > 0.0 );
				}
			else
				p->length = 0.0;
			}
		}

#if 0	
	printf ("after\n");
    ShowNodes (t->root, 0, YES);
	getchar();
#endif

	return (NO_ERROR);
	
}





/*-------------------------------------------------------
|
|   InitClockBrlens: This routine will initialize
|      a clock tree by setting the root to depth 1.0
|      and then assigning node depths according to
|      the relative node depth measured in terms of the
|      maximum number of branches to the tip from each
|      node.
|
--------------------------------------------------------*/
int InitClockBrlens (Tree *t)

{

	int				i, maxBrSegments=0;
	TreeNode		*p;

	if (t->isRooted == NO)
		{
		MrBayesPrint ("%s   Tree is unrooted\n", spacer);
		return (ERROR);
		}
	
	/* calculate maximum number of branch segments above root */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->anc != NULL)
			{
			if (p->left == NULL && p->right == NULL)
				{
				p->x = 0;
				}
			else
				{
				if (p->left->x > p->right->x)
					p->x = p->left->x + 1;
				else
					p->x = p->right->x + 1;
				}
			if (p->anc->anc == NULL)
				maxBrSegments = p->x;
			}
		}

	/* assign node depths */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->anc != NULL)
			p->nodeDepth = (MrBFlt) (p->x) / (MrBFlt) maxBrSegments;
		else
			p->nodeDepth = 0.0;
		}
		
	/* calculate branch lengths */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->anc != NULL)
			{
			if (p->anc->anc != NULL)
				p->length = p->anc->nodeDepth - p->nodeDepth;
			else
				p->length = 0.0;
			}
		}

	return (NO_ERROR);
	
}





int GetRandomEmbeddedSubtree (Tree *t, int nTerminals, SafeLong *seed, int *nEmbeddedTrees)

{
	
	int			i, j, k, n, ran, *pP, *pL, *pR, nLeaves, *nSubTrees;
	TreeNode	*p=NULL, **leaf;

	/* Calculate number of leaves in subtree (number of terminals minus the root) */
	nLeaves = nTerminals - 1;
	
	/* Initialize all flags */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		p->marked = NO;
		p->x = 0;
		p->y = 0;
		}
	
	/* Allocate memory */
	nSubTrees = (int *) SafeCalloc (nTerminals * t->nNodes, sizeof(int));
	if (!nSubTrees)
		return (ERROR);
	leaf = (TreeNode **) SafeMalloc (nLeaves * sizeof (TreeNode *));
	if (!leaf)
		{
		free (nSubTrees);
		return (ERROR);
		}

	/* Calculate how many embedded trees rooted at each node */
	(*nEmbeddedTrees) = 0;
	for (i=0; i<t->nNodes-1; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			{
			p->x = 0;
			nSubTrees[p->index*nTerminals + 1] = 1;
			}
		else
			{
			pL = nSubTrees + p->left->index*nTerminals;
			pR = nSubTrees + p->right->index*nTerminals;
			pP = nSubTrees + p->index*nTerminals;
			pP[1] = 1;
			for (j=2; j<=nLeaves; j++)
				{
				for (k=1; k<j; k++)
					{
					pP[j] += pL[k] * pR[j-k];
					}
				}
			p->x = pP[nLeaves];
			(*nEmbeddedTrees) += p->x;
			}
		}

	/* Randomly select one embedded tree of the right size */
	ran = (int) (RandomNumber(seed) * (*nEmbeddedTrees));

	/* Find the interior root corresponding to this tree */
	for (i=j=0; i<t->nIntNodes; i++)
		{
		p = t->intDownPass[i];
		j += p->x;
		if (j>ran)
			break;
		}

	/* Find one random embedded tree with this root */
	p->y = nLeaves;
	p->marked = YES;
	leaf[0] = p;
	n = 1;
	while (n < nLeaves)
		{
		/* select a node with more than one descendant */
		for (i=0; i<n; i++)
			{
			p = leaf[i];
			if (p->y > 1)
				break;
			}

		/* break it into descendants */
		pL = nSubTrees + p->left->index*nTerminals;
		pR = nSubTrees + p->right->index*nTerminals;
		pP = nSubTrees + p->index*nTerminals;
		ran = (int) (RandomNumber (seed) * pP[p->y]);
		k = 0;
		for (j=1; j<p->y; j++)
			{
			k += pL[j] * pR[p->y-j];
			if (k > ran)
				break;
			}

			p->left->y = j;
		p->right->y = p->y - j;
		p->left->marked = YES;
		p->right->marked = YES;
		leaf[i] = p->left;
		leaf[n++] = p->right;
		}

	free (nSubTrees);
	free (leaf);

	return (NO_ERROR);
}


		
		
		
/*-----------------------------------------------------------------------------
|
| IsCalibratedClockSatisfied: This routine sets calibrated clock tree params
|     and checks that user defined brlens satisfy the specified calibration(s)
|     up to tolerance tol
|
|------------------------------------------------------------------------------*/
int IsCalibratedClockSatisfied (Tree *t, MrBFlt tol)

{

	int				i, j, maxRateConstrained, minRateConstrained, isViolated;
	MrBFlt			f, maxHeight, minRate=0, maxRate=0, ageToAdd, *x, *y, clockRate;
	TreeNode		*p, *q, *r, *s;

	if (t->isRooted == NO)
		return (NO);
		
	x = (MrBFlt *) SafeCalloc (2*t->nNodes, sizeof (MrBFlt));
	if (x == NULL)
		{
		MrBayesPrint ("%s   Out of memory in IsCalibratedClockSatisfied\n", spacer);
		free (x);
		return (NO);
		}
	y = x + t->nNodes;

	/* reset node depth and age, and set minimum (x) and maximum (y) age of each node */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		p->age = -1.0;
		p->nodeDepth = -1.0;
		if (p->isDated == YES)
			{
			if (p->calibration->prior == fixed)
				x[p->index] = y[p->index] = p->calibration->age;
			else if (p->calibration->prior == uniform)
				{
				x[p->index] = p->calibration->min;
				y[p->index] = p->calibration->max;
				}
			else /* if (p->calibration->prior == offsetExponential) */
				{	
				x[p->index] = p->calibration->offset;
				y[p->index] = -1.0;
				}
			}
		else if (p->left == NULL && p->right == NULL)
			x[p->index] = y[p->index] = 0.0;
		else
			{
			x[p->index] = y[p->index] = -1.0;
			}
		}

	/* calculate node heights in branch length units */
	p = t->root->left;
	p->nodeDepth = 0.0;
	for (i=t->nNodes-3; i>=0; i--)
		{
		p = t->allDownPass[i];
		p->nodeDepth = p->anc->nodeDepth + p->length;
		}

	/* find maximum height of tree */	
	maxHeight = -1.0;
	for (i=0; i<t->nNodes-1; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL && p->right == NULL)
			{
			if (p->nodeDepth > maxHeight)
				{
				maxHeight = p->nodeDepth;
				}
			}
		}
	
	/* calculate node depth from tip of tree */
	for (i=0; i<t->nNodes-1; i++)
		{
		p = t->allDownPass[i];
		p->nodeDepth = maxHeight - p->nodeDepth;
		}
		

	/* check potentially constraining calibrations */
	/* and find minum and maximum possible rate */
	maxRateConstrained = NO;
	minRateConstrained = NO;
	isViolated = NO;
	for (i=0; i<t->nNodes-1; i++)
		{
		p = t->allDownPass[i];
		if (x[p->index] < 0.0 && y[p->index] < 0.0)
			continue;
		for (j=i+1; j<t->nNodes-1; j++)
			{
			q = t->allDownPass[j];
			if (x[q->index] < 0.0 && y[q->index] < 0.0)
				continue;
			if (AreDoublesEqual (p->nodeDepth, q->nodeDepth, tol) == YES)
				{
				/* same depth so they must share a possible age */
				if ((AreDoublesEqual (x[p->index], y[q->index], tol) == NO && x[p->index] > y[q->index])
					|| (AreDoublesEqual (y[p->index], x[q->index], tol) == NO && y[p->index] < x[q->index]))
					{
					isViolated = YES;
					break;
					}
				}
			else
				{
				if (p->nodeDepth > q->nodeDepth)
					{
					r = p;
					s = q;
					}
				else
					{
					r = q;
					s = p;
					}
				if (x[r->index] >= 0.0 && y[s->index] >= 0.0)
					{
					f = (r->nodeDepth - s->nodeDepth) / (x[r->index] - y[s->index]);
					if (f <= 0.0)
						{
						isViolated = YES;
						break;
						}
					if (maxRateConstrained == NO)
						{
						maxRateConstrained = YES;
						maxRate = f;
						}
					else if (f < maxRate)
						maxRate = f;
					}
				if (y[r->index] >= 0.0 && x[s->index] >= 0.0)
					{
					f = (r->nodeDepth - s->nodeDepth) / (y[r->index] - x[s->index]);
					if (f <= 0.0)
						{
						isViolated = YES;
						break;
						}
					if (minRateConstrained == NO)
						{
						minRateConstrained = YES;
						minRate = f;
						}
					else if (f > minRate)
						minRate = f;
					}
				}
			}
		if (isViolated == YES)
			break;
		}

	/* check if outright violation */
	if (isViolated == YES)
		{
		MrBayesPrint ("%s   Branch lengths do not satisfy the calibration(s)\n", spacer);
		free (x);
		return (NO);
		}
	
	/* check that minimum and maximum rates are consistent */
	if (minRateConstrained == YES && maxRateConstrained == YES
		&& AreDoublesEqual (minRate, maxRate, tol) == NO && minRate > maxRate)
		{
		MrBayesPrint ("%s   Branch lengths do not satisfy the calibration(s)\n", spacer);
		free (x);
		return (NO);
		}

	/* date all nodes based on a suitable rate */
	if (minRateConstrained == YES)
		clockRate = minRate;
	else if (maxRateConstrained == YES)
		clockRate = 0.5 * maxRate;
	else
		clockRate = 1.0;
	for (i=0; i<t->nNodes-1; i++)
		{
		p = t->allDownPass[i];
		p->age = p->nodeDepth / clockRate;
		}

	/* check if there is an age to add */
	ageToAdd = 0.0;
	for (i=0; i<t->nNodes-1; i++)
		{
		p = t->allDownPass[i];
		if (x[p->index] > 0.0 && x[p->index] > p->age)
			{
			f = x[p->index] - p->age;
			if (f > ageToAdd)
				ageToAdd = f;
			}
		}
	
	/* add extra length if any */
	if (AreDoublesEqual (ageToAdd, 0.0, 0.00000001) == NO)
		{
		for (i=0; i<t->nNodes-1; i++)
			{
			p = t->allDownPass[i];
			p->age += ageToAdd;
			}
		}

	free (x);

    /* reset node depths to ensure that non-dated tips have node depth 0.0 */
    SetNodeDepths(t);

	return (YES);
	
}





int IsClockSatisfied (Tree *t, MrBFlt tol)

{

	int				i, foundFirstLength, isClockLike;
	MrBFlt			firstLength=0.0, length;
	TreeNode		*p, *q;

	if (t->isRooted == NO)
		return (NO);
		
	foundFirstLength = NO;
	isClockLike = YES;
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL && p->right == NULL)
			{
			length = 0.0;
			q = p;
			while (q->anc != NULL)
				{
				if (q->anc->anc != NULL)
					length += q->length;
				q = q->anc;
				}
			if (foundFirstLength == NO)
				{
				firstLength = length;
				foundFirstLength = YES;
				}
			else
				{
				if (AreDoublesEqual (firstLength, length, tol) == NO)
					isClockLike = NO;
				}
			}
		}
	if (firstLength < BRLENS_MIN)
		isClockLike = NO;

	return (isClockLike);
	
}





/* Check that tree obeys topology constraints and that node depths and ages are consistent */
int IsTreeConsistent (Param *param, int chain, int state)
{
    Tree        *tree;
    TreeNode    *p;
    int         i, j;
    MrBFlt      b, r, rAnc, clockRate;

    if (param->paramType != P_TOPOLOGY && param->paramType != P_BRLENS && param->paramType != P_SPECIESTREE)
        return YES;

    tree      = GetTree(param, chain, state);
    if (modelSettings[param->relParts[0]].clockRate != NULL)
        clockRate = *GetParamVals(modelSettings[param->relParts[0]].clockRate, chain, state);
    else
        clockRate = 0.0;

    if (CheckConstraints(tree)==ERROR) {
        printf ("Tree does not obey constraints\n");
        return NO;
    }

    /* check that the last few indices are not taken in a rooted tree */
    if (tree->isRooted == YES && tree->root->index != tree->nNodes - 1)
		{
		printf("Problem with root index\n"); 
        return NO;
		}
    if (tree->isRooted == YES && tree->root->left->index != tree->nNodes - 2)
        {
		printf("Problem with interior root index\n");
		return NO;
		}

    if (tree->isClock == NO)
        {
        for (i=0; i<tree->nNodes-1; i++)
            {
            p = tree->allDownPass[i];
            if (p->length <= 0.0)
                {
                if(p->length == 0.0)
                    printf ("Node %d has zero branch length %f\n", p->index, p->length);
                else
                    printf ("Node %d has negative branch length %f\n", p->index, p->length);
                return NO;
                }
            }
        return YES;
        }

    /* Clock trees */

    /* Check that lengths and depths are consistent */
    for (i=0; i<tree->nNodes-2; i++) {
        p = tree->allDownPass[i];
        if (p->length <= 0.0) {
             if(p->length == 0.0)
                printf ("Node %d has zero branch length %f\n", p->index, p->length);
             else
                printf ("Node %d has negative branch length %f\n", p->index, p->length);
            return NO;
        }
        if (fabs(p->anc->nodeDepth - p->nodeDepth - p->length) > 0.000001) {
            printf ("Node %d has length %f but nodeDepth %f and ancNodeDepth %f\n",
                p->index, p->length, p->nodeDepth, p->anc->nodeDepth);
            return NO;
        }
        if (p->left == NULL && p->isDated == NO && p->nodeDepth != 0.0) {
                printf ("Node %d is a nondated tip but has node depth %f\n",
                    p->index, p->nodeDepth);
                return NO;
        }
    }

    /* Check that ages and calibrations are consistent */
    if (tree->isCalibrated == YES)
        {
        for (i=0; i<tree->nNodes-1; i++)
            {
            p = tree->allDownPass[i];
            if (p->isDated == YES) {
                if (fabs(p->age - p->nodeDepth/clockRate) > 0.000001)
                    {
                    printf ("Node %d has age %f but nodeDepth %f when clock rate is %f\n",
                        p->index, p->age, p->nodeDepth, clockRate);
                    return NO;
                    }
                if (p->calibration->prior == fixed && fabs(p->age - p->calibration->age) > 0.000001)
                    {
                    printf ("Node %d has age %f but should be fixed to age %f\n",
                        p->index, p->age, p->calibration->age);
                    return NO;
                    }
                else if (p->calibration->prior == offsetExponential && p->age < p->calibration->offset)
                    {
                    printf ("Node %d has age %f but should be minimally of age %f\n",
                        p->index, p->age, p->calibration->offset);
                    return NO;
                    }
                else if (p->calibration->prior == uniform && (p->age < p->calibration->min || p->age > p->calibration->max))
                    {
                    printf ("Node %d has age %f but should be in the interval (%f,%f)\n",
                        p->index, p->age, p->calibration->min, p->calibration->max);
                    return NO;
                    }
                }
            }
        }


    for (i=0; i<param->nSubParams; i++)
        {
        if (param->subParams[i]->paramId == TK02BRANCHRATES)
            {
            rAnc = GetParamVals(param->subParams[i], chain, state)[tree->root->left->index];
            if (fabs(rAnc - 1.0) > 1E-6)
                {
                MrBayesPrint("%s   TK02 relaxed clock mismatch in root rate, which is %e\n", spacer, rAnc);
                return NO;
                }
            for (j=0; j<tree->nNodes-2; j++)
                {
                p = tree->allDownPass[j];
                b = GetParamSubVals(param->subParams[i], chain, state)[p->index];
                r = GetParamVals(param->subParams[i], chain, state)[p->index];
                rAnc = GetParamVals(param->subParams[i], chain, state)[p->anc->index];
                if (fabs(p->length * (r + rAnc) / 2.0 - b) > 0.000001)
                    {
                    MrBayesPrint("%s   TK02 relaxed clock mismatch in branch %d\n", spacer, p->index);
                    return NO;
                    }
                }
            }
        else if (param->subParams[i]->paramId == IGRBRANCHLENS)
            {
            for (j=0; j<tree->nNodes-2; j++)
                {
                p = tree->allDownPass[j];
                b = GetParamSubVals(param->subParams[i], chain, state)[p->index];
                r = GetParamVals(param->subParams[i], chain, state)[p->index];
                if (fabs(p->length * r - b) > 0.000001)
                    {
                    MrBayesPrint("%s   Igr relaxed clock mismatch in branch %d\n", spacer, p->index);
                    return NO;
                    }
                }
            }
        }

    if (param->paramType == P_SPECIESTREE)
        return (IsSpeciesTreeConsistent(GetTree(param, chain, state), chain));

    return YES;
}





/* LabelTree: Label tree; remove previous labels if any */
int LabelTree (Tree *t, char **taxonNames)
{
	int			i, nTaxa;
	TreeNode	*p = NULL;

    nTaxa = t->nNodes - t->nIntNodes;
    if (t->isRooted == YES)
        nTaxa--;
    
    /* erase previous labels, if any */
	for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;
		t->nodes[i].label = noLabel;
        }

	/* add labels */
    for (i=0; i<t->nNodes; i++)
        {
		p = &t->nodes[i];
		if (p->left == NULL || (t->isRooted == NO && p->anc == NULL))
            {
            if (p->marked == YES || p->index < 0 || p->index >= nTaxa)
                {
                MrBayesPrint ("%s   Taxon node index repeated or out of range\n", spacer);
                return (ERROR);
                }
            else
                p->label = taxonNames[p->index];
            p->marked = YES;
            }
        else if (p->index > 0 && p->index < nTaxa)
            {
            MrBayesPrint ("%s   Terminal taxon index set for interior node\n", spacer);
            return (ERROR);
            }
        }

	return (NO_ERROR);
}





/*-------------------------------------------------------------------------------------------
|
|   Mark: This routine will mark up a subtree rooted at p
|
---------------------------------------------------------------------------------------------*/
void Mark (TreeNode *p)
{
    if (p != NULL)
        {
        p->marked = YES;
        Mark (p->left);
        Mark (p->right);
        }
}





/*-------------------------------------------------------------------------------------------
|
|   MarkDatedSubtree: This routine will mark up a subtree rooted at p with dated tips, whether
|       internal or external
|
---------------------------------------------------------------------------------------------*/
void MarkDatedSubtree (TreeNode *p)
{
    if (p != NULL)
        {
        p->marked = YES;
        MarkDatedSubtreeNodes (p->left);
        MarkDatedSubtreeNodes (p->right);
        }
}





/* MarkDatedSubtreeNodes: Recursive function to mark parts of a dated subtree */
void MarkDatedSubtreeNodes (TreeNode *p)
{
    if (p != NULL)
        {
        p->marked = YES;
        if (p->isDated == NO && p->left != NULL)
            {
            MarkDatedSubtreeNodes (p->left);
            MarkDatedSubtreeNodes (p->right);
            }
        }
}






/*-------------------------------------------------------------------------------------------
|
|   MoveCalculationRoot: This routine will move the calculation root to the terminal with 
|      index outgroup
|
---------------------------------------------------------------------------------------------*/
int MoveCalculationRoot (Tree *t, int outgroup)

{

	int 			i;
	TreeNode		*p, *q, *r;
    
    if (t->isRooted == YES || outgroup < 0 || outgroup > t->nNodes - t->nIntNodes - (t->isRooted == YES ? 1 : 0))
		{
		MrBayesPrint ("%s   Problem moving calculation root\n", spacer);
		return (ERROR);
		}

    if (t->root->index == outgroup)
		return (NO_ERROR);    /* nothing to do */

	/* mark the path to the new calculation root */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL && p->right == NULL)
			{
			if (p->index == outgroup)
				p->marked = YES;
			else
				p->marked = NO;
			}
		else
			{
			if (p->left->marked == YES || p->right->marked == YES)
				p->marked = YES;
			else
				p->marked = NO;
			}
		}

	/* rotate the tree to use the specified calculation root */
	p = t->root->left;
	q = t->root;
	q->anc = p;
	q->left = q->right = NULL;
	q->length = p->length;
	while (p->left != NULL && p->right != NULL)
		{
		if (p->left->marked == YES)
			{
			r = p->left;
			p->anc = r;
			p->left = q;
			p->length = r->length;
			q = p;
			p = r;
			}
		else /* if (p->right->marked == YES) */
			{
			r = p->right;
			p->anc = r;
			p->right = q;
			p->length = r->length;
			q = p;
			p = r;
			}
		}
	p->left = p->anc;
	p->right = p->anc = NULL;
	t->root = p;
	p->length = 0.0;

	GetDownPass (t);

	return (NO_ERROR);
}





/*-------------------------------------------------------------------------------------------
|
|   MovePolyCalculationRoot: This routine will move the calculation root to the terminal with 
|      index outgroup and place it as the right-most descendant of the root node
|
---------------------------------------------------------------------------------------------*/
int MovePolyCalculationRoot (PolyTree *t, int outgroup)

{

	int 			i;
	PolyNode		*p = NULL, *q, *r;

	/* check if tree is rooted, in which case calculation root is irrelevant */
	if (t->root->left->sib->sib == NULL)
		return (NO_ERROR);

	if (outgroup < 0 || outgroup > t->nNodes - t->nIntNodes)
		{
		MrBayesPrint ("%s   Outgroup index is out of range\n", spacer);
		return (ERROR);
		}

	if (t->root->left->sib->sib->sib != NULL)
		{
		MrBayesPrint ("%s   Root has more than three descendants\n", spacer);
		return (ERROR);
		}

    /* check if rerooting actually necessary */
	if (t->root->left->sib->sib->index == outgroup)
        return (NO_ERROR);
    
    /* mark the path to the new calculation root */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->index == outgroup)
			break;
		}
	if (p->left != NULL)
		{
		MrBayesPrint ("%s   Outgroup index set for internal node\n", spacer);
        for (i=0; i<t->nNodes; i++)
            printf ("%d -- %d\n", i, t->allDownPass[i]->index);
		getchar();
        return (ERROR);
		}

	/* mark path to current root */
	for (i=0; i<t->nNodes; i++)
		t->allDownPass[i]->mark = NO;
	q = p;
	while (q != NULL)
		{
		q->mark = YES;
		q = q->anc;
		}

	/* rotate the tree to use the specified calculation root */
	p = t->root;
	for (;;)
		{
		/* find marked descendant */
		for (q=p->left; q->mark == NO; q=q->sib)
			;
		if (q->index == outgroup)
			break;
		/* add old root to descendants of that node */
		for (r=q->left; r->sib!=NULL; r=r->sib)
			;
		r->sib = p;
		p->sib = NULL;   /* should not be needed */
		p->anc = q;
		p->length = q->length;
		/* remove that node from descendants of old root node */
		if (p->left == q)
			p->left = q->sib;
		else
			{
			for (r=p->left; r->sib!=q; r=r->sib)
				;
			r->sib = r->sib->sib;
			}
		/* make new node root */
		q->sib = NULL;
		q->anc = NULL;
		q->length = 0.0;
		p = q;
		}
	
	/* p is now the new root */
	t->root = p;

	/* finally make sure calculation root is last node among root's descendants */
	for (q=p->left; q->sib!=NULL; q=q->sib)
		;
	if (q->index != outgroup)
		{
		if (p->left->index == outgroup)
			{
			q->sib = p->left;
			p->left = p->left->sib;
			q->sib->sib = NULL;
			}
		else
			{
			for (r=p->left; r->sib->index!=outgroup; r=r->sib)
				;
			q->sib = r->sib;
			r->sib = r->sib->sib;
			q->sib->sib = NULL;
			}
		}

	GetPolyDownPass (t);

	return (NO_ERROR);
}





/* 
@return the number of levels for the tree rooted at the "node" 
*/
int NrSubTreeLevels(TreeNode *node)
{
	int r,l;

	if( node == NULL )
		{
		return -1;
		}

	r = NrSubTreeLevels( node->right );
	l = NrSubTreeLevels( node->left );

	return ((r>l)?(r):(l))+1;
}





/*-------------------------------------------------------------------------------------------
|
|   NumConstrainedTips: This routine will return the number of constrained tips, internal or external
|
---------------------------------------------------------------------------------------------*/
int NumConstrainedTips (TreeNode *p)
{
    int     i = 0;

    if (p == NULL)
        return i;
    if (p->left == NULL)
        return 1;

    i += NConstrainedTips (p->left);
    i += NConstrainedTips (p->right);

    return i;
}





/* NConstrainedTips: Recursive function to get the number of constrained tips */
int NConstrainedTips (TreeNode *p)
{
    int     i=0;
    
    if (p!=NULL)
        {
        if (p->left == NULL || p->isLocked == YES)
            return 1;
        else
            {
            i += NConstrainedTips (p->left);
            i += NConstrainedTips (p->right);
            }
        }
    return i;
}





/*-------------------------------------------------------------------------------------------
|
|   NumDatedTips: This routine will return the number of dated tips, internal or external
|
---------------------------------------------------------------------------------------------*/
int NumDatedTips (TreeNode *p)
{
    int     i = 0;

    assert (p != NULL && p->left != NULL);

    i += NDatedTips (p->left);
    i += NDatedTips (p->right);

    return i;
}





/* NDatedTips: recursive function to get the number of dated tips */
int NDatedTips (TreeNode *p)
{
    int     i=0;
    
    assert(p!=NULL);

    if (p->left == NULL || p->isDated == YES)
        return 1;
    else
        {
        i += NDatedTips (p->left);
        i += NDatedTips (p->right);
        return i;
        }
}





/* OrderTips: Order tips in a polytomous tree */
void OrderTips (PolyTree *t)
{
    int         i, j;
    PolyNode    *p, *q, *r, *pl, *ql, *rl;

	/* label by minimum index */
	for (i=0; i<t->nNodes; i++)
		{
	    p = t->allDownPass[i];
		if (p->left == NULL)
			{
			if (t->isRooted == NO && p->index == localOutGroup)
				p->x = -1;
			else
				p->x = p->index;
			}
		else
			{
			j = t->nNodes;
			for (q=p->left; q!=NULL; q=q->sib)
				{
				if (q->x < j)
					j = q->x;
				}
			p->x = j;
			}
		}

    /* and rearrange */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL || p->anc == NULL)
			continue;
		for (ql=NULL, q=p->left; q->sib!=NULL; ql=q, q=q->sib)
			{
			for (rl=q, r=q->sib; r!=NULL; rl=r, r=r->sib)
				{
				if (r->x < q->x)
					{
					if (ql == NULL)
						p->left = r;
					if (r == q->sib) /* swap adjacent q and r */
						{
						if (ql != NULL)
							ql->sib = r;
						pl = r->sib;
						r->sib = q;
						q->sib = pl;
						}
					else	/* swap separated q and r */
						{
						if (ql != NULL)
							ql->sib = r;
						pl = r->sib;
						r->sib = q->sib;
						rl->sib = q;
						q->sib = pl;
						}
					pl = q;
					q = r;
					r = pl;
					}
				}
			}
		}
    GetPolyDownPass(t);
}





/* PrintNodes: Print a list of tree nodes, pointers and length */
void PrintNodes (Tree *t)
{
	int			i;
	TreeNode	*p;

	printf ("Node\tleft\tright\tanc\tlength\n");
	for (i=0; i<t->nNodes; i++)
		{
		p = &t->nodes[i];
		MrBayesPrint ("%d\t%d\t%d\t%d\t%f\t%f\n",
			p->index,
			p->left == NULL ? -1 : p->left->index,
			p->right == NULL ? -1 : p->right->index,
			p->anc == NULL ? -1 : p->anc->index,
			p->length,
			p->nodeDepth);
		}

    if (t->root == NULL)
        MrBayesPrint ("root: NULL\n");
    else
        MrBayesPrint ("root: %d\n", t->root->index);

    MrBayesPrint ("allDownPass:");
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p!=NULL)
            MrBayesPrint ("  %d", p->index);
        else
            MrBayesPrint ("  NULL");
        }
    MrBayesPrint ("\nintDownPass:  ");
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        if (p!=NULL)
            MrBayesPrint ("  %d\t", p->index);
        else
            MrBayesPrint ("  NULL\t");
        }
    MrBayesPrint ("\n");
}





/* PrintPolyNodes: Print a list of polytomous tree nodes, pointers and length */
void PrintPolyNodes (PolyTree *pt)
{
	int			i, j, k;
	PolyNode	*p;

	printf ("Node\tleft\tsib\tanc\tlength\tlabel\n");
	for (i=0; i<pt->memNodes; i++)
		{
		p = &pt->nodes[i];
		MrBayesPrint ("%d\t%d\t%d\t%d\t%f\t%s\n",
			p->index,
			p->left == NULL ? -1 : p->left->index,
			p->sib == NULL ? -1 : p->sib->index,
			p->anc == NULL ? -1 : p->anc->index,
			p->length,
            p->label);
		}
	MrBayesPrint ("root: %d\n", pt->root->index); 
	fflush(stdout);

    if (pt->nBSets > 0)
        {
        for (i=0; i<pt->nBSets; i++)
            {
            printf ("Effective branch length set '%s'\n", pt->bSetName[i]);
            for (j=0; j<pt->nNodes; j++)
                {
                printf ("%d:%f", j, pt->effectiveBrLen[pt->nBSets][j]);
                if (j != pt->nNodes-1)
                    printf (", ");
                }
            printf ("\n");
            }
        }
    
    if (pt->nESets > 0)
        {
        for (i=0; i<pt->nESets; i++)
            {
            printf ("Cpp event set '%s'\n", pt->eSetName[i]);
            for (j=0; j<pt->nNodes; j++)
                {
                if (pt->nEvents[i*pt->nNodes+j] > 0)
                    {
                    printf ("\tNode %d -- %d:(", j, pt->nEvents[i][j]);
                    for (k=0; k<pt->nEvents[i][j]; k++)
                        {
                        printf ("%f %f", pt->position[i][j][k], pt->rateMult[i][j][k]);
                        if (k != pt->nEvents[i][j]-1)
                            printf (", ");
                        }
                    printf (")\n");
                    }
                }
            printf ("\n");
            }
        }

    fflush(stdout);
}





/**
Update relaxed clock paramiter of the branch of a node with index "b" after node with index "a" is removed. 
i.e. make branch of node with index "b" be a concatenation of its original branch and the branch of node with index "a"
Relaxed clock paramiter of node with index "a" become invalid in the process.
Note: For Non-clock models the routine has no effect. 

|       |
|       |
a       |
|   ->  |
|       |
|       b
b

*/
void AppendRelaxedBranch (int a,int b,PolyTree *t){

    int i,len;

    for (i=0; i<t->nBSets; i++)
        {
        t->effectiveBrLen[i][b] += t->effectiveBrLen[i][a];
        }

    for (i=0; i<t->nESets; i++)
        {
        len=t->nEvents[i][a]+t->nEvents[i][b];
        t->position[i][a] = (MrBFlt *) SafeRealloc ((void *)t->position[i][a], len*sizeof(MrBFlt));
	    t->rateMult[i][a] = (MrBFlt *) SafeRealloc ((void *)t->rateMult[i][a], len*sizeof(MrBFlt));
        memcpy( t->position[i][a]+t->nEvents[i][a], t->position[i][b], t->nEvents[i][b]*sizeof(MrBFlt));
        memcpy( t->rateMult[i][a]+t->nEvents[i][a], t->rateMult[i][b], t->nEvents[i][b]*sizeof(MrBFlt));
        free(t->position[i][b]);
        free(t->rateMult[i][b]);
        t->position[i][b] = t->position[i][a];
        t->rateMult[i][b] = t->rateMult[i][a];
        t->position[i][a] = NULL;
        t->rateMult[i][a] = NULL;
        t->nEvents[i][a] = 0;
        t->nEvents[i][b] = len;
        }

}





/**
Swap relaxed clock paramiters of the branch of nodes with index "a" and "b".
*/
void SwapRelaxedBranchInfo (int a,int b,PolyTree *t){

    int i,j;
    MrBFlt tmp, *tmpp;

    for (i=0; i<t->nBSets; i++)
        {
        tmp = t->effectiveBrLen[i][a];
        t->effectiveBrLen[i][a] = t->effectiveBrLen[i][b];
        t->effectiveBrLen[i][b] = tmp;
        }
    if (t->popSizeSet == YES)
        {
        tmp = t->popSize[a];
        t->popSize[a]=t->popSize[b];
        t->popSize[b] = tmp;
        }

    for (i=0; i<t->nESets; i++)
        {
        tmpp = t->position[i][a];
        t->position[i][a] = t->position[i][b];
        t->position[i][b] = tmpp;
        tmpp = t->rateMult[i][a];
	    t->rateMult[i][a] = t->rateMult[i][b];
        t->rateMult[i][b] = tmpp;
        j = t->nEvents[i][a];
        t->nEvents[i][a] = t->nEvents[i][b];
        t->nEvents[i][b] = j;
        }
}





/*-------------------------------------------------------------------------------------------
|
|   PrunePolyTree: This routine will prune a polytomous tree according to the currently
|      included taxa.  NB! All tree nodes cannot be accessed by cycling over the
|      pt->nodes array after the deletion, because some spaces will be occupied by deleted
|      nodes and pt->nNodes is no longer the length of this array
|      (if proper re-arangment of pt->nodes needed then remove "#if 0" and make changes to p->x, see below).
|
---------------------------------------------------------------------------------------------*/
int PrunePolyTree (PolyTree *pt)

{

	int 			i, j, numDeleted, numTermPruned, numIntPruned, index;
	PolyNode		*p = NULL, *q=NULL, *r=NULL, *qa;

	numDeleted = 0;
	for (i=0; i<pt->nNodes; i++)
		{
		p = pt->allDownPass[i];
        CheckString (taxaNames, numTaxa, p->label, &index);
        if (p->left == NULL && taxaInfo[index].isDeleted == YES)
			numDeleted++;
		}
		
    if (numDeleted == 0)
		{
		/* nothing to do */
		return (NO_ERROR);
		}
	if (pt->nNodes - pt->nIntNodes - numDeleted < 3)
		{
		MrBayesPrint ("%s   Pruned tree has less than three taxa in it\n", spacer);
		return (ERROR);
		}
	if (pt->nNodes - pt->nIntNodes < numLocalTaxa)
		{
		MrBayesPrint ("%s   Tree to be pruned does not include all taxa\n", spacer);
		return (ERROR);
		}

	/* prune away one node at a time */
	numIntPruned = 0;
    numTermPruned = 0;
	for (i=0; i<pt->nNodes; i++)
		{
		p = pt->allDownPass[i];
		if (p->left != NULL)
			continue;
        CheckString (taxaNames, numTaxa, p->label, &index);
		if (taxaInfo[index].isDeleted == YES)
			{
            numTermPruned++;
			for (q=p->anc->left; q!=NULL; q=q->sib)
				{
				if (q->sib == p)
					break;
				}
			if (q == NULL)
				{
				/* p is the left of its ancestor */
                assert( p->anc->left == p);
				p->anc->left = p->sib;
				}
			else
				{
				/* p is q->sib; this also works if p->sib is NULL */
				q->sib = p->sib;
				}
			/* if only one child left, delete ancestral node */
			j = 0;
			for (q=p->anc->left; q!=NULL; q=q->sib)
				j++;
			if (j == 1)
				{
				/* p->anc->left is only child left; make p->anc be p->anc->left and accommodate its length */
				numIntPruned++;
                qa= p->anc;
                q = qa->left;
				if (q->left == NULL)
                    {
                    AppendRelaxedBranch ( qa->index, q->index, pt );
                    qa->index = q->index;
                    qa->length += q->length;
                    strcpy(qa->label, q->label);
                    qa->left = NULL;
                    /* To make sure that q is not treated as the representer of the tip it represented before. i.e. make  condition if (p->left != NULL) true */
                    q->left = (struct pNode*)1; 
                    }
                else
                    {
                    if (qa->anc != NULL)
                        {
                        AppendRelaxedBranch ( qa->index, q->index, pt );
                        qa->length += q->length;
                        }
                    qa->index   = q->index;
                    qa->left = q->left;
                    for(r=q->left; r!= NULL; r=r->sib)
                        r->anc = qa;
                    }
				}
            /* if unrooted, then root node has to have more then 2 children, thus the following check */
            if (j == 2 && pt->isRooted == NO && p->anc->anc == NULL)
                {
                numIntPruned++;
                r=p->anc; /*r is the root with only 2 children*/
                if ( r->left->left != NULL )
                    {/* Make r->left new root by attaching "right" child of r to children of r->left */
                    for (q=r->left->left; q->sib!=NULL; q=q->sib)
                        ;
                    q->sib = r->left->sib;
                    r->left->sib->anc = q->anc;
                    r->left->sib->length += q->anc->length;
                    r->left->sib = NULL;
                    r->left->anc = NULL;
                    pt->root = r->left;
                    }
                else
                    {/* Make "right" child of r (r->left->sib) the new root by attaching r->left to children of r->"right" */
                    for ( q=r->left->sib->left; q->sib!=NULL; q=q->sib )
                        ;
                    q->sib = r->left;
                    r->left->anc = q->anc;
                    r->left->length += q->anc->length;
                    r->left->sib = NULL;
                    q->anc->anc = NULL;
                    pt->root = q->anc;
                    }
                }
			}
		}

#if 0 
	/* place unused space at the end of pt->nodes array. If activated this code p->x has to be set to non 0 value for all p that are deleted. */
	for (i=0; i<pt->nNodes; i++)
		{
		p = &pt->nodes[i];
		if (p->x != 0)
			{
			for (j=i+1; j<pt->nNodes; j++)
				{
				q = &pt->nodes[j];
				if (q->x == 0)
					break;
				}
			if (j != pt->nNodes)
				{
				/* swap nodes; quite difficult! */
				CopyPolyNodes (p, q, nLongsNeeded);
				p->left = q->left;
				p->sib = q->sib;
				p->anc = q->anc;
				for (k=0; k<pt->nNodes; k++)
					{
					r = &pt->nodes[k];
					if (r->left == q)
						r->left = p;
					if (r->sib == q)
						r->sib = p;
					if (r->anc == q)
						r->anc = p;
					}
				}
			}
		}
#endif

	/* correct number of nodes */
	pt->nNodes -= (numTermPruned + numIntPruned);
	pt->nIntNodes -= numIntPruned;
	
	/* get downpass; note that the deletion procedure does not change the root in rooted case */
    i=j=0;
    GetPolyNodeDownPass ( pt, pt->root, &i, &j );
    assert(i==pt->nNodes);
    assert(j==pt->nIntNodes);

	return (NO_ERROR);
	
}





/*--------------------------------------------------------------------
|
|		RandPerturb: Randomly perturb a tree by nPert NNIs
|
---------------------------------------------------------------------*/
int RandPerturb (Tree *t, int nPert, SafeLong *seed)
{
	
	int			i, whichNode;
	TreeNode	*p, *q, *a, *b, *c;
	
	if (t->nConstraints >= t->nIntNodes)
		{
		MrBayesPrint ("%s   User tree cannot be perturbed because all nodes are locked\n", spacer);
		return (ERROR);
		}

	for (i=0; i<nPert; i++)
		{
		do
			{
			whichNode = (int)(RandomNumber(seed) * (t->nIntNodes - 1));
			p = t->intDownPass[whichNode];
			} while (p->isLocked == YES);
		
		q = p->anc;
		a  = p->left;
		b  = p->right;
		if (q->left == p)
			c  = q->right;
		else	
			c  = q->left;
		
		if (RandomNumber(seed) < 0.5)
			{
			/* swap b and c */
			p->right = c;
			c->anc  = p;

			if (q->left == c)
				q->left = b;
			else
				q->right = b;
			b->anc = q;
			}
		else
			{
			/* swap a and c */
			p->left = c;
			c->anc  = p;

			if (q->left == c)
				q->left = a;
			else
				q->right = a;
			a->anc = q;
			}

		if (t->isCalibrated == YES)
			InitCalibratedBrlens (t, 0.0001, seed);
		else if (t->isClock == YES)
			InitClockBrlens (t);
		}
	
	GetDownPass (t);

	if (t->checkConstraints == YES && CheckConstraints (t) == NO_ERROR)
		{
		MrBayesPrint ("%s   Broke constraints when perturbing tree\n", spacer);
		return (ERROR);
		}

	return (NO_ERROR);
}





/*
|       Reorder array of nodes "nodeArray" such that first nodes in it could be paired with "w" to create imediat common ancestor and this ancesor node would not vialate any constraint.  
|       
| @param w                      Reference node as described 
| @param nodeArray              A set of node to order as described
| @param activeConstraints      Array containing indeces of active constraints in the set of defined constraints
| @param nLongsNeeded           Length of partition information (in SafeLong) in a node and constraint deffinition.
| @param isRooted               Do constraints apply to rootet tree YES or NO
|
| @return                       Number of nodes in "nodeArray" that could be paired  with "w" to create imediat common ancestor and this ancesor node would not vialate any constraint
*/
int ConstraintAllowedSet(PolyNode *w, PolyNode **nodeArray, int nodeArraySize, int *activeConstraints, int activeConstraintsSize, int nLongsNeeded, int isRooted )
{

    int             i, j,  k, FirstEmpty;
    SafeLong        **constraintPartition;
    PolyNode        *tmp;

    for (j=0; j<activeConstraintsSize; j++)
        {
        k=activeConstraints[j];

        if( definedConstraintsType[k] == PARTIAL )
            {
            if( ( IsPartNested(definedConstraintPruned[k], w->partition, nLongsNeeded) == YES ) ||
                ( isRooted == NO && IsPartNested(definedConstraintTwoPruned[k], w->partition, nLongsNeeded) == YES )
              )
                continue;/* all nodes are compartable because condition of the constraint has to be sutsfied in the subtree rooted at w*/

            FirstEmpty = IsSectionEmpty(definedConstraintPruned[k], w->partition, nLongsNeeded);
            if( FirstEmpty == YES &&  IsSectionEmpty(definedConstraintTwoPruned[k], w->partition, nLongsNeeded) == YES )
                continue; /* all nodes are compartable becouse w does not contain any constraint taxa*/

            assert(FirstEmpty^IsSectionEmpty(definedConstraintTwoPruned[k], w->partition, nLongsNeeded));

            if( FirstEmpty == YES )
                {/*w->partition has intersection with definedConstraintTwoPruned[k], thus remove all nodes from nodeArray that intersect with definedConstraintPruned[k]*/
                constraintPartition=definedConstraintPruned;
                }
            else
                {/*w->partition has intersection with definedConstraintPruned[k], thus remove all nodes from nodeArray that intersect with definedConstraintTwoPruned[k]*/
                constraintPartition=definedConstraintTwoPruned;
                }

            for(i=0;i<nodeArraySize;i++)
                {
                if( IsSectionEmpty(constraintPartition[k], nodeArray[i]->partition, nLongsNeeded) == NO &&
                    ( ( FirstEmpty == NO && isRooted== YES ) ||  IsPartNested(constraintPartition[k], nodeArray[i]->partition, nLongsNeeded) == NO  )
                  )
                  /*second part of if statment is to bail out "nodeArray[i]" when "w" contains nodes for example from definedConstraintPruned and "nodeArray[i]" have definedConstraintTwoPruned fully nested in it
                  This bail out not applicable if t->isRooted== YES Since we should create a rooting node for the first set of taxa in the constraint.
                  Note that such case possible because we may have hard constraint node that fully nest definedConstraintTwoPruned but also having taxa from definedConstraintPruned keeping constraint active.*/
                    {
                    tmp = nodeArray[i];
                    nodeArray[i]=nodeArray[--nodeArraySize];
                    nodeArray[nodeArraySize]=tmp;
                    i--;
                    }
                }
            }/*end if PARTIAL*/
        else 
            {
            assert( definedConstraintsType[k] == NEGATIVE );
            if( isRooted == YES || IsBitSet(localOutGroup, definedConstraintPruned[k])==NO )
                constraintPartition=definedConstraintPruned;
            else
                constraintPartition=definedConstraintTwoPruned;
            
            if( IsSectionEmpty(constraintPartition[k], w->partition, nLongsNeeded)==YES )
                continue;

            for(i=0;i<nodeArraySize;i++)
                {
                if( IsUnionEqThird (w->partition, nodeArray[i]->partition, constraintPartition[k], nLongsNeeded) == YES )
                    {
                    tmp = nodeArray[i];
                    nodeArray[i]=nodeArray[--nodeArraySize];
                    nodeArray[nodeArraySize]=tmp;
                    i--;
                    }
                }

            }/*end if NEGATIVE*/
        }

   return nodeArraySize;
}




/*
|               Check if "partition" violate any constraint.  
|       
| @param partiton               Partition to check 
| @param activeConstraints      Array containing indeces of active constraints in the set of defined constraints  
| @param nLongsNeeded           Length of partition information (in SafeLong) in a node and constraint deffinition
| @param isRooted               Do constraints apply to rootet tree YES or NO
|
| @return                       Index of first violated constraint in activeConstraints array, -1 if no constraint is violated.
*/
int ViolatedConstraint(SafeLong *partition, int *activeConstraints, int activeConstraintsSize, int nLongsNeeded, int isRooted )
{

    int             j, k;
    SafeLong        **constraintPartition;


    for (j=0; j<activeConstraintsSize; j++)
        {
        k=activeConstraints[j];
        assert(definedConstraintsType[k] != HARD);

        if( definedConstraintsType[k] == PARTIAL )
            {
            if( ( IsSectionEmpty(definedConstraintPruned[k], partition, nLongsNeeded) == NO ) &&
                ( IsSectionEmpty(definedConstraintTwoPruned[k], partition, nLongsNeeded) == NO ) &&
                ( IsPartNested(definedConstraintPruned[k], partition, nLongsNeeded) == NO) &&
                !( isRooted == NO && IsPartNested(definedConstraintTwoPruned[k], partition, nLongsNeeded) == YES)
              )
                return j;

            }/*end if PARTIAL*/
        else 
            {
            assert( definedConstraintsType[k] == NEGATIVE );
            if( isRooted == YES || IsBitSet(localOutGroup, definedConstraintPruned[k])==NO )
                constraintPartition=definedConstraintPruned;
            else
                constraintPartition=definedConstraintTwoPruned;

            if( IsUnionEqThird (partition, partition, constraintPartition[k], nLongsNeeded) == YES )
                return j;

            }/*end if NEGATIVE*/
        }

   return -1;
}






/*
|         Remove from activeConstraints references to constraints that become satisfied if PolyNode "w" exist, i.e. they do not need to be checked furter thus become not active  
|
| @param activeConstraints      Array containing indeces of active constraints in the set of defined constraints
| @param nLongsNeeded           Length of partition information (in SafeLong) in a node and constraint deffinition.
| @param isRooted               Do constraints apply to rootet tree YES or NO
|
| @return                       Size of pruned "activeConstraints" array
*/
int PruneActiveConstraints (PolyNode *w, int *activeConstraints, int activeConstraintsSize, int nLongsNeeded, int isRooted )
{

    int             j,  k;
    SafeLong        **constraintPartition;
    //PolyNode        *tmp;

    for (j=0; j<activeConstraintsSize; j++)
        {
        k=activeConstraints[j];

        if( definedConstraintsType[k] == PARTIAL )
            {
            if( (IsPartNested(definedConstraintPruned[k], w->partition, nLongsNeeded) == YES && IsSectionEmpty(definedConstraintTwoPruned[k], w->partition, nLongsNeeded)) ||
                ( isRooted == NO && IsPartNested(definedConstraintTwoPruned[k], w->partition, nLongsNeeded) == YES && IsSectionEmpty(definedConstraintPruned[k], w->partition, nLongsNeeded))
              )
                {
                //tmp = activeConstraints[j];
                activeConstraints[j]=activeConstraints[--activeConstraintsSize];
                //activeConstraints[activeConstraintsSize]=tmp;
                j--;
                }
            }/*end if PARTIAL*/
        else 
            {
            assert( definedConstraintsType[k] == NEGATIVE );
            if( isRooted == YES || IsBitSet(localOutGroup, definedConstraintPruned[k])==NO )
                constraintPartition=definedConstraintPruned;
            else
                constraintPartition=definedConstraintTwoPruned;
            
            if( IsPartNested(constraintPartition[k], w->partition, nLongsNeeded)==NO && IsSectionEmpty(constraintPartition[k], w->partition, nLongsNeeded)==NO )
                {
                //tmp = activeConstraints[j];
                activeConstraints[j]=activeConstraints[--activeConstraintsSize];
                //activeConstraints[activeConstraintsSize]=tmp;
                j--;
                }
            }/*end if NEGATIVE*/
        }

   return activeConstraintsSize;
}





/*--------------------------------------------------------------------
|
|		    RandResolve: Randomly resolve a polytomous tree
|
| @param    tt is a tree which contains information about applicable constraints. If it is set to NULL then no constraints will be used. 
            If t!=NULL then partitions of nodes of polytree should be allocated for example by AllocatePolyTreePartitions (t);
| @return   NO_ERROR on succes, ABORT if could not resolve a tree without vialating some consraint, ERROR if any other error occur 
---------------------------------------------------------------------*/
int RandResolve (Tree *tt, PolyTree *t, SafeLong *seed, int destinationIsRooted)

{

	int			i, j, k, nextNode, stopNode, rand1, rand2, nLongsNeeded, tmp;
	PolyNode	*p=NULL, *q, *r, *u, *w1, *w2;
    int         nodeArrayAllowedSize, nodeArraySize, activeConstraintsSize;
    PolyNode	**nodeArray;
    int         *activeConstraints;

    assert(tt==NULL || t->bitsets!=NULL); /* partition fields of t nodes need to be allocated if constraints are used*/
    assert( numLocalTaxa <= t->memNodes/2); /* allocated tree has to be big enough*/
    nLongsNeeded = (numLocalTaxa - 1) / nBitsInALong + 1; /* allocated lenght of partitions is t->memNodes/2 bits but only first numLocalTaxa bits are used */

    nodeArray = t->allDownPass; /*temporary use t->allDownPass for different purpose. It get properly reset at the end. */
    activeConstraints = tempActiveConstraints;
    activeConstraintsSize = 0;

    /* collect constraints to consider if applicable*/
    if( tt!=NULL && tt->constraints!=NULL )
        {
        for (k=0; k<numDefinedConstraints; k++)
            {
            if( tt->constraints[k] == YES && definedConstraintsType[k] != HARD )
                activeConstraints[activeConstraintsSize++]=k;
            }
        }

	/* count immediate descendants */
	GetPolyDownPass(t);
	for (i=0; i<t->nIntNodes; i++)
		{
		p = t->intDownPass[i];
        tmp=ViolatedConstraint(p->partition, activeConstraints, activeConstraintsSize, nLongsNeeded, t->isRooted);
        if ( tmp != -1)
            {
            assert( p->isLocked == YES );
            MrBayesPrint ("%s   Could not build a constraint tree since hard constraint \"%s\" and constraint \"%s\" are incompatible\n", spacer, constraintNames[p->lockID], constraintNames[activeConstraints[tmp]]);
            return (ERROR);
            }
        activeConstraintsSize = PruneActiveConstraints (p, activeConstraints, activeConstraintsSize, nLongsNeeded, t->isRooted);
		j = 0;
		for (q=p->left; q!=NULL; q=q->sib)
			j++;
		p->x = j;
		}

	/* add one node at a time */
	if (destinationIsRooted == NO)
		stopNode = 2*numLocalTaxa - 2;
	else
		stopNode = 2*numLocalTaxa - 1;
	for (nextNode=t->nNodes; nextNode < stopNode; nextNode++)
		{
		/* find a polytomy to break */
		for (i=0; i<t->nIntNodes; i++)
			{
			p = t->intDownPass[i];
			if (destinationIsRooted == YES && p->x > 2)
				break;
			if (destinationIsRooted == NO && ((p->anc != NULL && p->x > 2) || (p->anc == NULL && p->x > 3)))
				break;
			}

		/* if we can't find one, there's an error */
		if (i == t->nIntNodes)
            {
			return  ERROR;
            }

        nodeArraySize=0;
        /*Collect initial list of condidate nodes to join*/
        for (q = p->left; q!= NULL; q = q->sib)
            {
            nodeArray[nodeArraySize++]=q;
            }
        assert(nodeArraySize==p->x);

        /* identify two descendants randomly */
		/* make sure we do not select outgroup if it is an unrooted tree */
		if (p->anc == NULL && destinationIsRooted == NO)
			nodeArraySize--;

        do
            {
            /*Pick first node*/
            rand1 = (int) (RandomNumber(seed) * nodeArraySize);
            w1 = nodeArray[rand1];
            nodeArray[rand1] = nodeArray[--nodeArraySize];

            if( nodeArraySize==0 )
                return ABORT; /*Potentaily here we could instead revert by removing last added node and try again. */

            /*Move all nodes in nodeArray which can be paired with w to the begining of array*/
            nodeArrayAllowedSize=ConstraintAllowedSet(w1, nodeArray, nodeArraySize, activeConstraints, activeConstraintsSize, nLongsNeeded, t->isRooted );
            /*TODO optimization for Maxim(if not Maxim remove it if you still see it): if nodeArrayAllowedSize==0 then set w1->y */
            }while( nodeArrayAllowedSize==0 );

		rand2 = (int) (RandomNumber(seed) *nodeArrayAllowedSize);
        w2 = nodeArray[rand2];


		/* create a new node */
		u = &t->nodes[nextNode];
		u->anc = p;
		u->x = 2;
		p->x--;

        if( tt!=NULL ){
            for (j=0; j<nLongsNeeded; j++)
			    u->partition[j] = w1->partition[j] | w2->partition[j] ;
            activeConstraintsSize = PruneActiveConstraints (u, activeConstraints, activeConstraintsSize, nLongsNeeded, t->isRooted );
        }

        u->left = w1;
		t->nNodes++;
		t->nIntNodes++;
		

		/* connect tree together */
        r = u;
        for (q = p->left; q!= NULL; q = q->sib)
			{
			if ( q != w1 && q!=w2)
				{
                r->sib=q;
                r = q;
				}
			}
		r->sib = NULL;
        w1->sib = w2;
        w2->sib = NULL;
        w1->anc = u;
        w2->anc = u;
        p->left = u;

		/* update tree */
		GetPolyDownPass (t);
		}

    /* relabel interior nodes (important that last indices are at the bottom!) */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        p->index = numLocalTaxa + i;
        }
	return NO_ERROR;
}





/* ResetTreeNode: Reset tree node except for memory index */
void ResetTreeNode (TreeNode *p)
{
	/* do not change memoryIndex; that is set once and for all when tree is allocated */
	p->index                  = 0; 
	p->scalerNode			  = NO;			
	p->upDateCl               = NO;
	p->upDateTi				  = NO;
	p->marked                 = NO;
	p->length                 = 0.0;
	p->nodeDepth              = 0.0;
	p->x                      = 0;
	p->y                      = 0;
	p->index				  = 0;
	p->isDated				  = NO;
	p->calibration			  = NULL;
	p->age					  = -1.0;
	p->isLocked				  = NO;
	p->lockID				  = -1;
	p->label                  = noLabel;
	p->d					  = 0.0;
	p->partition			  = NULL;
}





/* ResetPolyNode: Reset all values of one node in a polytree */
void ResetPolyNode (PolyNode *p)
{
    /* we reset everything here except memoryIndex, which should be immutable */
    p->length = 0.0;
    p->depth = 0.0;
    p->age = 0.0;
    p->anc = p->left = p->sib = NULL;
    p->calibration = NULL;
    p->f = 0.0;
    p->index = 0;
    p->isDated = NO;
    p->isLocked = NO;
    strcpy (p->label,"");
    p->lockID = 0;
    p->partition = NULL;
    p->partitionIndex = 0;
    p->support = 0.0;
    p->x = p->y = 0;
}





/* ResetPolyTree: Reset polytomous tree to pristine state but keep relevant memory. */
void ResetPolyTree (PolyTree *pt)
{
	int	    i, maxTaxa, nLongsNeeded;

	/* clear nodes */
	for (i=0; i<pt->memNodes; i++)
		ResetPolyNode (&pt->nodes[i]);

    /* empty node arrays and tree properties but keep space */
    for (i=0; i<pt->nNodes; i++)
        pt->allDownPass[i] = NULL;
    for (i=0; i<pt->nIntNodes; i++)
        pt->intDownPass[i] = NULL;
    pt->nNodes = 0;
    pt->nIntNodes = 0;
    pt->root = NULL;
    pt->brlensDef = NO;
    pt->isRooted = NO;
    pt->isClock = NO;
    pt->isRelaxed = NO;
    pt->clockRate = 0.0;

    /* empty bitsets but keep space and pointers */
    if (pt->bitsets)
        {
        maxTaxa = pt->memNodes / 2;
        nLongsNeeded = (maxTaxa - 1) / nBitsInALong + 1;
        for (i=0; i<pt->memNodes*nLongsNeeded; i++)
            pt->bitsets[i] = 0;
        for (i=0; i<pt->memNodes; i++)
            pt->nodes[i].partition = pt->bitsets + i*nLongsNeeded;
        }

    /* empty relaxed clock parameters */
    FreePolyTreeRelClockParams (pt);

    /* empty population size set parameters */
    FreePolyTreePopSizeParams (pt);
}





/* ResetPolyTreePartitions: Reset and set bit patterns describing partitions */
void ResetPolyTreePartitions (PolyTree *pt)
{
    int         i, j, numTaxa, nLongsNeeded;
    PolyNode    *pp;

    /* get some handy numbers */
    numTaxa = pt->memNodes/2;
    nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;
    
    /* reset bits describing partitions */
    for (i=0; i<pt->memNodes*nLongsNeeded; i++)
		{
        pt->bitsets[i] = 0;
        }

    /* set bits describing partitions */
    for (i=0; i<pt->nNodes; i++)
		{
		assert (pt->allDownPass != NULL && pt->allDownPass[i] != NULL);
        assert (pt->allDownPass[i]->partition != NULL);
        
        pp = pt->allDownPass[i];
		if (pp->left == NULL)
            {
			SetBit (pp->index, pp->partition);
            }
		if (pp->anc != NULL)
			{
			for (j=0; j<nLongsNeeded; j++)
				pp->anc->partition[j] |= pp->partition[j];
			}
		}
}





/*----------------------------------------------
|
|   ResetRootHeight: Reset node heights in a clock
|      tree to fit a new root height. Assumes
|      node depths and lengths set correctly.
|
-----------------------------------------------*/
int ResetRootHeight (Tree *t, MrBFlt rootHeight)
{
	int         i;
    TreeNode	*p;
    MrBFlt      factor, x, y;

    if (t->isClock == NO)
        return ERROR;
    
    /* make sure node depths are set */
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->nodeDepth = 0.0;
        else
            {
            x = p->left->nodeDepth + p->left->length;
            y = p->right->nodeDepth + p->right->length;
            if (x > y)
                p->nodeDepth = x;
            else
                p->nodeDepth = y;
            }
        }
    for (i=t->nNodes-3; i>=0; i--)
        {
        p = t->allDownPass[i];
        p->nodeDepth = p->anc->nodeDepth - p->length;
        }

    /* now reset node depths and branch lengths */
    factor = rootHeight / t->root->left->nodeDepth;
    t->root->left->nodeDepth = rootHeight;
    for (i=t->nNodes-2; i>=0; i--)
        {
        p = t->allDownPass[i];
        p->nodeDepth *= factor;
        p->length *= factor;
        }

    return NO_ERROR;
}





/*----------------------------------------------
|
|   ResetTipIndices: reset tip indices to be from 
|      0 to number of included taxa, in same order
|      as in the original taxon set.
|
-----------------------------------------------*/
void ResetTipIndices(PolyTree *pt)
{
    int         i, j, k, m;
    PolyNode    *p;


    for (i=j=0; i<numTaxa; i++)
		{
		for (k=0; k<pt->nNodes; k++)
			{
			p = pt->allDownPass[k];
			if (StrCmpCaseInsensitive(p->label,taxaNames[i]) == 0)
				break;
			}
        if (k < pt->nNodes)
            {
            assert (p->left == NULL);
            if(p->index!=j){
                SwapRelaxedBranchInfo( p->index, j, pt);
                for(m=0; m<pt->nNodes; m++)
                    {
                    if(pt->allDownPass[m]->index==j)
                        {
                        pt->allDownPass[m]->index=p->index;
                        break;
                        }
                    }
		        p->index = j;
                }
            j++;
            }
		}
}





/*----------------------------------------------
|
|   ResetTopology: rebuild the tree t to fit the 
|      Newick string s. Everyting except topology
|      is left in the same state in t.
|
-----------------------------------------------*/
int ResetTopology (Tree *t, char *s)
{
	TreeNode	*p, *q;
	int			i, j, k, inLength;
	char		temp[30];
	
    /* set all pointers to NULL */
	for (i=0; i<t->memNodes; i++)
		{
		p = &t->nodes[i];
		p->anc = p->right = p->left = NULL;
		p->index = -1;
		}
	p = &t->nodes[0];

    /* start out assuming that the tree is rooted; we will detect below if it is not */
    t->isRooted = YES;
	inLength = NO;
	for (i=0, j=1; *s!='\0'; s++)
		{
		if (*s == ',' || *s == ')' || *s == ':')
			{
			if (p->right == NULL && inLength == NO)
				{
				temp[i] = '\0';
				k = atoi (temp);
				p->index = k-1;
				i = 0;
				}
			else
				inLength = NO;
			}
		if (*s == '(')
			{
			q = p;
			p = &t->nodes[j++];
			q->left = p;
			p->anc = q;
			}
		else if (*s == ',')
			{
			if (p->anc->right == NULL)
				{
				q = p->anc;
				p = &t->nodes[j++];
				p->anc = q;
				q->right = p;
				}
			else /* if p->anc->right == p (near 'root' of unrooted trees) */
				{
				q = p->anc;
				p = &t->nodes[j++];
				q->anc = p;
				p->left = q;
                t->isRooted = NO;
				}
			}
		else if (*s == ')')
			{
			p = p->anc;
			}
		else if (*s == ':')
			{
			inLength = YES;
			}
		else if (inLength == NO)
			{
			temp[i++] = *s;
			}
		}

	/* attach root to rooted tree */
	if (t->isRooted == YES)
		{
		p = &t->nodes[0];
		q = &t->nodes[j++];
		q->left = p;
		p->anc = q;
		}

    /* relabel interior nodes, find number of nodes and root */
    t->nNodes = j;
    t->nIntNodes = t->nNodes/2 - 1;

    if (t->isRooted == NO)
        j = t->nNodes - t->nIntNodes;
    else
        j = t->nNodes - t->nIntNodes - 1;

	for (i=0; i<t->nNodes; i++)
		{
		p = &t->nodes[i];
		if (p->index == -1)
			p->index = j++;
		if (p->anc == NULL)
			t->root = p;
		}

	GetDownPass (t);

	return NO_ERROR;
}





/*-----------------------------------------------------------------
|
|	ResetBrlensFromTree: copies brlens from second tree (vTree) to
|       first tree (used to initialize brlen sets for same topology)
|
-----------------------------------------------------------------*/
int ResetBrlensFromTree (Tree *tree, Tree *vTree)

{

	int			i, j, k, nLongsNeeded, numTips;
    MrBFlt      d1, d2;
	TreeNode	*p, *q;

	if (tree->isRooted != vTree->isRooted)
		return (ERROR);
	
	if (AreTopologiesSame (tree, vTree) == NO)
		return (ERROR);

	/* allocate and set up partitions */
	AllocateTreePartitions (vTree);
	AllocateTreePartitions (tree);
	numTips = tree->nNodes - tree->nIntNodes - (tree->isRooted == YES ? 1 : 0);
	nLongsNeeded = (int) ((numTips - 1) / nBitsInALong) + 1;

	for (i=0; i<vTree->nNodes; i++)
		{
		p  = vTree->allDownPass[i];
		for (j=0; j<tree->nNodes; j++)
			{
			q  = tree->allDownPass[j];
			for (k=0; k<nLongsNeeded; k++)
				if (p->partition[k] != q->partition[k])
					break;
			if (k==nLongsNeeded)
				{
				q->length = p->length;
                if (tree->isRooted == YES)
                    q->nodeDepth = p->nodeDepth;
				}
			}
		}
    if (tree->isRooted == YES)
        {
    	for (i=0; i<tree->nNodes-1; i++)
	    	{
	    	p  = tree->allDownPass[i];
            if (p->left == NULL)
                p->nodeDepth = 0.0;
            else
                {
                d1 = p->left->nodeDepth + p->left->length;
                d2 = p->right->nodeDepth + p->right->length;
			    if (d1 > d2)
                    p->nodeDepth = d1;
                else
                    p->nodeDepth = d2;
                }
            }
        for (i=tree->nNodes-3; i>=0; i--)
            {
            p = tree->allDownPass[i];
            if (p->left==NULL && p->calibration==NULL)
                continue;    /* leave at 0.0 */
            p->nodeDepth = p->anc->nodeDepth - p->length;
            }
		}

	FreeTreePartitions(tree);
	FreeTreePartitions(vTree);
	
	return (NO_ERROR);
		
}





/* ResetIntNodeIndices: Set int node indices in downpass order from numTaxa to 2*numTaxa-2 */
void ResetIntNodeIndices (PolyTree *t)
{
    int i, m, index;

    index = t->nNodes - t->nIntNodes;

    for (i=0; i<t->nIntNodes; i++)
        {
        if(t->intDownPass[i]->index != index)
            {
            SwapRelaxedBranchInfo( t->intDownPass[i]->index, index, t);
            for(m=0; m<t->nIntNodes; m++)
                {
                if(t->intDownPass[m]->index==index)
                    {
                    t->intDownPass[m]->index=t->intDownPass[i]->index;
                    break;
                    }
                }
	        t->intDownPass[i]->index = index;
            }
        index++;
        }
}





/* ResetTopologyFromTree: use top to set topology in tree */
int ResetTopologyFromTree (Tree *tree, Tree *top)
{
	int			i, j, k;
	TreeNode	*p, *q, *r, *p1;

    /* adopt rooting */
    tree->isRooted = top->isRooted;
    tree->nNodes = top->nNodes;
    tree->nIntNodes = top->nIntNodes;
	
    /* set all pointers to NULL */
	for (i=0; i<tree->nNodes; i++)
		{
		p = &tree->nodes[i];
		p->anc = p->right = p->left = NULL;
        }

    /* now copy topology */
    for (i=0; i<top->nIntNodes; i++)
		{
		p1 = top->intDownPass[i];
		
		k = p1->index;
		for (j=0; j<tree->nNodes; j++)
			if (tree->nodes[j].index == k)
				break;
		p = &tree->nodes[j];

		k = p1->left->index;
		for (j=0; j<tree->nNodes; j++)
			if (tree->nodes[j].index == k)
				break;
		q = &tree->nodes[j];

		k = p1->right->index;
		for (j=0; j<tree->nNodes; j++)
			if (tree->nodes[j].index == k)
				break;
		r = &tree->nodes[j];

		p->left = q;
		p->right= r;
		q->anc = r->anc = p;
		}

	/* arrange the root */
	k = top->root->index;
	for (j=0; j<tree->nNodes; j++)
		if (tree->nodes[j].index == k)
			break;
	p = &tree->nodes[j];

	k = top->root->left->index;
	for (j=0; j<tree->nNodes; j++)
		if (tree->nodes[j].index == k)
			break;
	q = &tree->nodes[j];
	p->left = q;
	q->anc = p;
	p->right = p->anc = NULL;
	tree->root = p;

	GetDownPass(tree);

	return (NO_ERROR);
}





/* ResetTopologyFromPolyTree: use polytree top to set topology in tree */
int ResetTopologyFromPolyTree (Tree *tree, PolyTree *top)
{
	int			i, j, k;
	TreeNode	*p, *q, *r;
    PolyNode    *p1;

    if (tree->isRooted != top->isRooted)
		return (ERROR);
	
    /* set all pointers to NULL */
	for (i=0; i<tree->nNodes; i++)
		{
		p = &tree->nodes[i];
		p->anc = p->right = p->left = NULL;
        }

    /* now copy topology */
    for (i=0; i<top->nIntNodes; i++)
		{
		p1 = top->intDownPass[i];
		
		k = p1->index;
		for (j=0; j<tree->nNodes; j++)
			if (tree->nodes[j].index == k)
				break;
		p = &tree->nodes[j];

		k = p1->left->index;
		for (j=0; j<tree->nNodes; j++)
			if (tree->nodes[j].index == k)
				break;
		q = &tree->nodes[j];

		k = p1->left->sib->index;
		for (j=0; j<tree->nNodes; j++)
			if (tree->nodes[j].index == k)
				break;
		r = &tree->nodes[j];

		p->left = q;
		p->right= r;
		q->anc = r->anc = p;
		}

	/* arrange the root */
	if (top->isRooted == YES)
        {
        k = top->root->index;
	    for (j=0; j<tree->nNodes; j++)
		    if (tree->nodes[j].index == k)
			    break;
	    p = &tree->nodes[j];

	    k = top->nNodes;
	    for (j=0; j<tree->nNodes; j++)
		    if (tree->nodes[j].index == k)
			    break;
	    q = &tree->nodes[j];

        q->left = p;
	    q->anc = NULL;
        q->right = NULL;
	    tree->root = q;
        }
    else /* if (top->isRooted == NO) */
    {
        k = top->root->index;
	    for (j=0; j<tree->nNodes; j++)
		    if (tree->nodes[j].index == k)
			    break;
	    p = &tree->nodes[j];

        k = localOutGroup;
        for (p1=top->root->left; p1!=NULL; p1=p1->sib)
            if (p1->index == k)
                break;

        assert (p1 != NULL);
        if (p1 == NULL)
            return (ERROR);

        q = &tree->nodes[p1->index];
        k = p1->anc->left->sib->sib->index;     /* index of missing child */
        if (p->left == q)
            p->left = &tree->nodes[k];
        else if (p->right == q)
            p->right = &tree->nodes[k];

        q->anc = q->right = NULL;
        p->anc = q;
        q->left = p;
    }

	GetDownPass(tree);

	return (NO_ERROR);
}





/* ResetTreePartitions: Reset bitsets describing tree partitions */
void ResetTreePartitions (Tree *t)
{
    int         i, j, numTaxa, nLongsNeeded;
    TreeNode    *p;

    /* get some handy numbers */
    numTaxa = t->nNodes - t->nIntNodes - (t->isRooted == YES ? 1 : 0);
    nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;
    
    /* reset bits describing partitions */
    for (i=0; i<t->nNodes; i++)
		{
		assert (t->allDownPass != NULL && t->allDownPass[i] != NULL);
        assert (t->allDownPass[i]->partition != NULL);
        
        p = t->allDownPass[i];
        for (j=0; j<nLongsNeeded; j++)
            p->partition[j] = 0;
        }

    /* set bits describing partitions */
    for (i=0; i<t->nNodes; i++)
		{
        p = t->allDownPass[i];
		if (p->left == NULL || (p->anc == NULL && t->isRooted == NO))
			SetBit (p->index, p->partition);
        else if (p->anc != NULL)
			{
			for (j=0; j<nLongsNeeded; j++)
				p->partition[j] = p->left->partition[j] | p->right->partition[j];
			}
		}
}





/*-------------------------------------------------------
|
|   RetrieveRTopology: This routine will rebuild a rooted
|      tree from the order array created by StoreRTopology.
|      All tree information except the structure will
|      remain unaltered.
|
--------------------------------------------------------*/
int RetrieveRTopology (Tree *t, int *order)

{
	int			i, numTaxa;
	TreeNode	*p, *q, *r;
	
	numTaxa = t->nNodes - t->nIntNodes - 1;
	
	/* sort the tips in the t->allDownPass array */
	p = t->nodes;
	for (i=0; i<t->nNodes; i++, p++)
		t->allDownPass[p->index] = p;

	/* make sure the root has index 2*numTaxa-1 */
	q = t->allDownPass[t->nNodes-1];
	q->anc = q->right = NULL;
	t->root = q;

	/* connect the first two tips */
	p = t->allDownPass[numTaxa];
	p->anc = q;
	q->left = p;
	p->length = 0.0;
	q = t->allDownPass[0];
	r = t->allDownPass[1];
	p->left = q;
	p->right = r;
	q->anc = r->anc = p;

	/* add one tip at a time */
	for (i=2; i<numTaxa; i++)
		{
		p = t->allDownPass[i];
		q = t->allDownPass[numTaxa-1+i];
		r = t->allDownPass[*(order++)];
		p->anc = q;
		q->left = p;
		q->right = r;
		q->anc = r->anc;
		if (r->anc->left == r)
			r->anc->left = q;
		else
			r->anc->right = q;
		r->anc = q;
		}

	/* get downpass */
	GetDownPass (t);

	/* relabel interior nodes (root is correctly labeled already) */
	for (i=0; i<t->nIntNodes; i++)
		t->intDownPass[i]->index = i+numTaxa;

	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   RetrieveRTree: This routine will rebuild a rooted
|      tree from the arrays created by StoreRTree.
|      All tree information except the structure and
|      branch lengths will remain unaltered.
|
--------------------------------------------------------*/
int RetrieveRTree (Tree *t, int *order, MrBFlt *brlens)

{
	int			i, numTaxa;
	TreeNode	*p, *q, *r;

	numTaxa = t->nNodes - t->nIntNodes - 1;
	
	/* sort the tips in the t->allDownPass array */
	p = t->nodes;
	for (i=0; i<t->nNodes; i++, p++)
		t->allDownPass[p->index] = p;

	/* make sure that root has index 2*numTaxa-1 */
	q = t->allDownPass[t->nNodes-1];
	q->anc = q->right = NULL;
	q->length = 0.0;
	t->root = q;

	/* connect the first three tips */
	p = t->allDownPass[numTaxa];
	p->anc = q;
	q->left = p;
	p->length = 0.0;
	q = t->allDownPass[0];
	r = t->allDownPass[1];
	p->left = q;
	p->right = r;
	q->anc = r->anc = p;
	q->length = *(brlens++);
	r->length = *(brlens++);

	/* add one tip at a time */
	for (i=2; i<numTaxa; i++)
		{
		p = t->allDownPass[i];
		q = t->allDownPass[numTaxa-1+i];
		r = t->allDownPass[*(order++)];
		p->anc = q;
		q->left = p;
		q->right = r;
		q->anc = r->anc;
		if (r->anc->left == r)
			r->anc->left = q;
		else
			r->anc->right = q;
		r->anc = q;
		if (q->anc->anc != NULL)
			q->length = *(brlens++);
		else
			{
			r->length = *(brlens++);
			q->length = 0.0;
			}
		p->length = *(brlens++);
		}

	/* get downpass */
	GetDownPass (t);

	/* relabel interior nodes (root is correctly labeled already) */
	for (i=0; i<t->nIntNodes; i++)
		t->intDownPass[i]->index = i+numTaxa;

	/* set the node depths */
	SetNodeDepths (t);
	
	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   RetrieveRTreeWithIndices: This routine will rebuild a rooted
|      tree from the arrays created by StoreRTreeWithIndices.
|      All tree information except the structure, branch lengths
|      and node indices will remain unaltered.
|
--------------------------------------------------------*/
int RetrieveRTreeWithIndices (Tree *t, int *order, MrBFlt *brlens)

{
	int			i, numTaxa;
	TreeNode	*p, *q, *r;

    extern void ShowNodes (TreeNode *, int, int);

	numTaxa = t->nNodes - t->nIntNodes - 1;
	
	/* sort the tips in the t->allDownPass array */
	p = t->nodes;
	for (i=0; i<t->nNodes; i++, p++)
		t->allDownPass[p->index] = p;

	/* make sure that root has index 2*numTaxa-1 */
	q = t->allDownPass[t->nNodes-1];
	q->anc = q->right = NULL;
	q->length = 0.0;
	t->root = q;

	/* connect the first three 'tips' with interior node, index from order array */
	p = t->allDownPass[numTaxa];
	p->x = *(order++);
    p->anc = q;
	q->left = p;
	p->length = 0.0;
	q = t->allDownPass[0];
	r = t->allDownPass[1];
	p->left = q;
	p->right = r;
	q->anc = r->anc = p;
	q->length = *(brlens++);
	r->length = *(brlens++);

	/* add one tip at a time */
	for (i=2; i<numTaxa; i++)
		{
		p = t->allDownPass[i];
        assert (*order >= numTaxa && *order < 2*numTaxa - 1);
		q = t->allDownPass[numTaxa-1+i];
		q->x = *(order++);
		r = t->allDownPass[*(order++)];
		p->anc = q;
		q->left = p;
		q->right = r;
		q->anc = r->anc;
		if (r->anc->left == r)
			r->anc->left = q;
		else
			r->anc->right = q;
		r->anc = q;
		if (q->anc->anc != NULL)
			q->length = *(brlens++);
		else
			{
			r->length = *(brlens++);
			q->length = 0.0;
			}
		p->length = *(brlens++);
		}

	/* get downpass */
	GetDownPass (t);

    /* relabel interior nodes using labels in scratch variable x */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        p->index = p->x;
        }

	/* set the node depths */
	SetNodeDepths (t);
	
	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   RetrieveUTopology: This routine will rebuild an unrooted
|      tree from the order array created by StoreUTopology.
|      All tree information except the structure
|      will remain unaltered.
|
--------------------------------------------------------*/
int RetrieveUTopology (Tree *t, int *order)

{
	int			i, numTips;
	TreeNode	*p, *q, *r;
	
    /* preliminaries */
    numTips = t->nNodes - t->nIntNodes;
	for (i=0; i<t->nNodes; i++)
        t->nodes[i].left = t->nodes[i].right = t->nodes[i].anc = NULL;

	/* sort the tips in the t->allDownPass array */
	p = t->nodes;
	for (i=0; i<t->nNodes; i++, p++)
		t->allDownPass[p->index] = p;

	/* make sure root has index 0 */
	q = t->allDownPass[0];
	q->anc = q->right = NULL;
	t->root = q;

	/* connect the first three tips */
	p = t->allDownPass[numTips];
	p->anc = q;
	q->left = p;
	q = t->allDownPass[1];
	r = t->allDownPass[2];
	p->left = q;
	p->right = r;
	q->anc = r->anc = p;

	/* add one tip at a time */
	for (i=3; i<numTips; i++)
		{
		p = t->allDownPass[i];
		q = t->allDownPass[numTips-2+i];
		r = t->allDownPass[order[i-3]];
		p->anc = q;
		q->left = p;
		q->right = r;
		q->anc = r->anc;
		if (r->anc->left == r)
			r->anc->left = q;
		else
			r->anc->right = q;
		r->anc = q;
		}

	/* get downpass */
	GetDownPass (t);
	
	/* relabel interior nodes (root is correctly labeled already) */
	for (i=0; i<t->nIntNodes; i++)
		t->intDownPass[i]->index = i+numTips;

	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   RetrieveUTree: This routine will rebuild an unrooted
|      tree from the arrays created by StoreUTree.
|      All tree information except the structure and
|      branch lengths will remain unaltered.
|
--------------------------------------------------------*/
int RetrieveUTree (Tree *t, int *order, MrBFlt *brlens)

{
	int			i, numTips;
	TreeNode	*p, *q, *r;
	
    /* preliminaries */
    numTips = t->nNodes - t->nIntNodes;
	for (i=0; i<t->nNodes; i++)
        t->nodes[i].left = t->nodes[i].right = t->nodes[i].anc = NULL;
	
	/* sort the tips in the t->allDownPass array */
	p = t->nodes;
	for (i=0; i<t->nNodes; i++, p++)
		t->allDownPass[p->index] = p;

	/* make sure that root has index 0 */
	q = t->allDownPass[0];
	q->anc = q->right = NULL;
	t->root = q;

	/* connect the first three tips */
	p = t->allDownPass[numTips];
	p->anc = q;
	q->left = p;
	p->length = *(brlens++);
	q = t->allDownPass[1];
	r = t->allDownPass[2];
	p->left = q;
	p->right = r;
	q->anc = r->anc = p;
	q->length = *(brlens++);
	r->length = *(brlens++);

	/* add one tip at a time */
	for (i=3; i<numTips; i++)
		{
		p = t->allDownPass[i];
		q = t->allDownPass[numTips-2+i];
		r = t->allDownPass[order[i-3]];
		p->anc = q;
		q->left = p;
		q->right = r;
		q->anc = r->anc;
		if (r->anc->left == r)
			r->anc->left = q;
		else
			r->anc->right = q;
		r->anc = q;
		q->length = *(brlens++);
		p->length = *(brlens++);
		}

	/* get downpass */
	GetDownPass (t);

	/* relabel interior nodes (root is correctly labeled already) */
	for (i=0; i<t->nIntNodes; i++)
		t->intDownPass[i]->index = i+numTips;

	return (NO_ERROR);
}





void SetDatedNodeAges (Param *param, int chain, int state)

{

	int		    i;
    MrBFlt      clockRate;
    ModelInfo   *m;
	TreeNode	*p;
    Tree        *t;

	extern void ShowNodes(TreeNode *,int,int);

    t = GetTree (param, chain, state);
    m = &modelSettings[t->relParts[0]];

    if (m->clockRate == NULL)
        clockRate = 1.0;
    else
        clockRate = *GetParamVals(m->clockRate, chain, state);

    for (i=0; i<t->nNodes-1; i++)
		{
		p = t->allDownPass[i];
		if (p->isDated == YES)
            p->age = p->nodeDepth / clockRate;
        else
            p->age = -1.0;
        }
}





void SetNodeDepths (Tree *t)

{

	int		i;
	MrBFlt		d1, d2;
	TreeNode	*p;

	extern void ShowNodes(TreeNode *,int,int);

	for (i=0; i<t->nNodes-1; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			p->nodeDepth = 0.0;
		else
			{
			d1 = p->left->nodeDepth  + p->left->length;
			d2 = p->right->nodeDepth + p->right->length;
            //assert (!(t->isCalibrated == NO && AreDoublesEqual(d1,d2,0.00001)==NO)); // may not work if we set startval topology of strict clock tree by non clock tree. 
			if (d1 > d2)
				p->nodeDepth = d1;
			else
				p->nodeDepth = d2;
			}
		}

	for (i=t->nNodes-3; i>=0; i--)
		{
		p = t->allDownPass[i];
		if (p->left == NULL && p->calibration == NULL)
			p->nodeDepth = 0.0;
		else
			p->nodeDepth = p->anc->nodeDepth - p->length;
		}
}





int ShowPolyNodes (PolyTree *pt)

{

	int 			i;
	PolyNode		*p;

	/* this is the tree, on a node-by-node basis */
    printf ("   memnodes = %d  nNodes = %d  nIntNodes = %d  root = %d\n", pt->memNodes, pt->nNodes, pt->nIntNodes, pt->root->index);
    printf ("   isRooted = %d\n", pt->isRooted);
    printf ("   no. index (left sib anc) -- locked/free -- label (p->x)\n");
	for (i=0; i<pt->memNodes; i++)
		{
		p = &pt->nodes[i];
		if (!(p->left == NULL && p->sib == NULL && p->anc == NULL))
			{
			printf ("%4d -- %4d ", i, p->index);
			if (p->left != NULL)
				printf ("(%4d ", p->left->index);
			else
				printf ("(null ");

			if (p->sib != NULL)
				printf ("%4d ", p->sib->index);
			else
				printf ("null ");
				
			if (p->anc != NULL)
				printf ("%4d)", p->anc->index);
			else
				printf ("null)");
			
			if (p->isLocked == YES)
				printf ("-- locked -- ");
            else
                printf ("-- free --");

			if (p->left == NULL && p->anc != NULL)
				printf ("  \"%s\" (%d)\n", p->label, p->x);
			else
				printf (" \"\" (%d)\n", p->x);
			}
		}

	return NO_ERROR;
}





/* ShowTree: Show tree on screen */
int ShowTree (Tree *t)

{

	int 			i, j, k, x, nLines, nLevels, levelDepth, from, to;
	char			treeLine[SCREENWIDTH2], labelLine[100];
	TreeNode		*p;
	
	/* get coordinates */
	x = 0;
	nLines = 0;
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL && p->right == NULL)
			{
			p->x = x;
			x += 2;
			p->y = 0;
			nLines += 2;
			}
		else if (p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			p->x = p->left->x + (p->right->x - p->left->x) / 2;
			if (p->left->y > p->right->y)
				p->y = p->left->y + 1;
			else
				p->y = p->right->y + 1;
			}
		else
			{
			p->x = x;
			x += 2;
			p->y = 0;
			}
		} 

    /* print tree out, line-by-line */
	levelDepth = SCREENWIDTH / t->root->left->y;
	nLevels = t->root->left->y;
	for (j=0; j<=nLines-2; j++)
		{
        for (i=0; i<SCREENWIDTH2-2; i++)
            treeLine[i] = ' ';
        treeLine[SCREENWIDTH-1] = '\n';
		if (j % 2 == 0)
			{
			for (i=0; i<t->nNodes; i++)
				{
				p = t->allDownPass[i];
				if (p->left == NULL && p->x == j)
					{
					strcpy (labelLine, p->label);
					}
				}
			}
		for (i=0; i<t->nNodes; i++)
			{
			p = t->allDownPass[i];
			if (p->anc != NULL)
				{
				if (p->anc->anc != NULL)
					{
					if (p->x == j)
						{
						from = (nLevels - p->anc->y) * levelDepth;
						to   = (nLevels - p->y) * levelDepth;
						if (p->y == 0)
							to = SCREENWIDTH-1;
						if (to >= SCREENWIDTH)
							to = SCREENWIDTH-1;
							
						for (k=from; k<to; k++)
							treeLine[k] = '-';
						if (p->anc->left == p)
							treeLine[from] = '/';
						else
							treeLine[from] = '\\';
						if (p->left != NULL)
							{
							treeLine[to] = '+';
							}
						if (p->anc->anc == t->root && p->anc->right == p)
							{
							if (t->isRooted == NO)
								treeLine[to] = '+';
							else
								treeLine[from] = '\\';
							}
						}
					else
						{
						if (p->left != NULL && p->right != NULL)
							{
							if (j < p->x && j > p->left->x)
								{
								from = (nLevels - p->y) * levelDepth;
								treeLine[from] = '|';
								}
							else if (j > p->x && j < p->right->x && p->left != NULL)
								{
								from = (nLevels - p->y) * levelDepth;
								treeLine[from] = '|';
								}
							}
						}
					}
				else
					{
					if (p->x == j)
						{
						treeLine[0] = '|'; /* temp */
						}
					else if (j < p->x && j > p->left->x)
						{
						treeLine[0] = '|';
						}
					else if (j > p->x && j < p->right->x)
						{
						treeLine[0] = '|';
						}
					if (t->isRooted == NO)
						{
						if (j > p->x && j <= nLines-2)
							treeLine[0] = '|';
						if (j == p->right->x)
							treeLine[0] = '+';
						}
					else
						{
						if (j == p->x)
							treeLine[0] = '+';
						}
					}
				}
			}
		treeLine[SCREENWIDTH-1] = '\0';
		if (j % 2 == 0)
			MrBayesPrint ("   %s %s\n", treeLine, labelLine);
		else
			MrBayesPrint ("   %s \n", treeLine);
		}

	if (t->isRooted == NO)
		{
		for (i=0; i<SCREENWIDTH; i++)
			treeLine[i] = ' ';
		treeLine[SCREENWIDTH-1] = '\0';
		MrBayesPrint ("   |\n");
		for (k=0; k<SCREENWIDTH; k++)
			treeLine[k] = '-';
		treeLine[SCREENWIDTH-1] = '\0';
		treeLine[0] = '\\';
		strcpy (labelLine, t->root->label);
		labelLine[19] = '\0';
		MrBayesPrint ("   %s %s\n", treeLine, labelLine);
		}
	
#if defined (DEBUG_CONSTRAINTS)
	for (i=0; i<t->nNodes; i++)
		printf ("%d -- %s\n", t->allDownPass[i]->index + 1, t->allDownPass[i]->isLocked == YES ? "locked" : "free");
#endif

	return (NO_ERROR);
	   
}





/*-------------------------------------------------------
|
|   StoreRPolyTopology: Same as StoreRTopology but for
|   binary polytree source trees.
|
--------------------------------------------------------*/
int StoreRPolyTopology (PolyTree *t, int *order)

{
	int			i, numTaxa;
	PolyNode	*p, *q;
	
	/* find number of taxa */
	numTaxa = t->nNodes - t->nIntNodes;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first two taxa */
		if (p->index > 1 && p->index < numTaxa)
			order[p->index-2] = i;
		}

	/* label the interior nodes with the correct index */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			p->x = p->y = p->index;
		else
			{
			if (p->left->y < p->left->sib->y)
				{
				p->y = p->left->y;
				p->x = p->left->sib->y + numTaxa - 1;
				}
			else
				{
				p->y = p->left->sib->y;
				p->x = p->left->y + numTaxa - 1;
				}
			}
		}

	/* break the tree into pieces */
	for (i=0; i<numTaxa-2; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTaxa-3-i]];
		q = p->anc;
		if (q->left == p)
			{
			order[numTaxa-3-i] = q->left->sib->x;
			p->sib->anc = q->anc;
            if (q->anc == NULL)
                {
                p->sib->left->sib->sib = p->sib->sib;
                p->sib->sib = NULL;
                }
			else if (q->anc->left == q)
                {
				q->anc->left = q->left->sib;
                p->sib->sib = q->sib;
                }
			else
				q->anc->left->sib = q->left->sib;
			}
		else
			{
			order[numTaxa-3-i] = q->left->x;
			q->left->anc = q->anc;
			if (q->anc == NULL)
                {
                q->left->left->sib->sib = p->sib;
                q->left->sib = NULL;
                }
			else if (q->anc->left == q)
                {
				q->anc->left = q->left;
                q->anc->left->sib = q->sib;
                }
			else
                {
				q->anc->left->sib = q->left;
                q->left->sib = NULL;
                }
			}
		}

	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   StoreRPolyTree: Same as StoreRTree but for
|      binary rooted polytree source trees.
|
--------------------------------------------------------*/
int StoreRPolyTree (PolyTree *t, int *order, MrBFlt *brlens)

{
	int			i, j, numTaxa;
	PolyNode	*p, *q;
	
	/* find number of taxa */
	numTaxa = t->nNodes - t->nIntNodes;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first two taxa */
		if (p->index > 1 && p->index < numTaxa)
			order[p->index-2] = i;
		}

	/* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			p->x = p->y = p->index;
		else
			{
			if (p->left->y < p->left->sib->y)
				{
				p->y = p->left->y;
				p->x = p->left->sib->y + numTaxa - 1;
				}
			else
				{
				p->y = p->left->sib->y;
				p->x = p->left->y + numTaxa - 1;
				}
			}
		}

	/* break the tree into pieces */
    j = t->nNodes - 2;     /* index of first branch length */
	for (i=0; i<numTaxa-2; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTaxa-3-i]];
		q = p->anc;
        brlens[j--] = p->length;
        brlens[j--] = q->length;
		if (q->left == p)
			{
			order[numTaxa-3-i] = q->left->sib->x;
			p->sib->anc = q->anc;
            if (q->anc == NULL)
                {
                p->sib->left->sib->sib = p->sib->sib;
                p->sib->sib = NULL;
                }
			else if (q->anc->left == q)
                {
				q->anc->left = q->left->sib;
                p->sib->sib = q->sib;
                }
			else
				q->anc->left->sib = q->left->sib;
			}
		else
			{
			order[numTaxa-3-i] = q->left->x;
			q->left->anc = q->anc;
			if (q->anc == NULL)
                {
                q->left->left->sib->sib = p->sib;
                q->left->sib = NULL;
                }
			else if (q->anc->left == q)
                {
				q->anc->left = q->left;
                q->anc->left->sib = q->sib;
                }
			else
                {
				q->anc->left->sib = q->left;
                q->left->sib = NULL;
                }
			}
		}

    /* store the last two lengths; index 0 and 1 */
    p = t->root;
    brlens[p->left->index] = p->left->length;
    brlens[p->left->sib->index] = p->left->sib->length;

	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   StoreRTopology: This routine will break a rooted tree
|      into an array of ints describing the structure
|      of the tree. The tree will be destroyed
|      in the process (the node pointers, that is).
|      However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreRTopology (Tree *t, int *order)

{
	int			i, numTaxa;
	TreeNode	*p, *q;
	
	/* find number of taxa */
	numTaxa = t->nNodes - t->nIntNodes - 1;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first two taxa */
		if (p->index > 1 && p->index < numTaxa)
			order[p->index-2] = i;
		}

	/* label the interior nodes with the correct index */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			p->x = p->y = p->index;
		else if (p->right != NULL)
			{
			if (p->left->y < p->right->y)
				{
				p->y = p->left->y;
				p->x = p->right->y + numTaxa - 1;
				}
			else
				{
				p->y = p->right->y;
				p->x = p->left->y + numTaxa - 1;
				}
			}
		}

	/* break the tree into pieces */
	for (i=0; i<numTaxa-2; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTaxa-3-i]];
		q = p->anc;
		if (q->left == p)
			{
			order[numTaxa-3-i] = q->right->x;
			q->right->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->right;
			else
				q->anc->right = q->right;
			}
		else
			{
			order[numTaxa-3-i] = q->left->x;
			q->left->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->left;
			else
				q->anc->right = q->left;
			}
		}

	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   StoreRTree: This routine will break a rooted tree
|      into an array of ints describing the structure
|      of the tree and an array of doubles storing
|      the branch lengths. The tree will be
|      destroyed in the process (the node pointers,
|      that is). However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreRTree (Tree *t, int *order, MrBFlt *brlens)

{
	int			i, j, numTaxa;
	TreeNode	*p, *q;

	extern void ShowNodes (TreeNode *p, int indent, int isRooted);

	/* find number of taxa */
	numTaxa = t->nNodes - t->nIntNodes - 1;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first two taxa */
		if (p->index > 1 && p->index < numTaxa)
			order[p->index-2] = i;
		}

	/* label the interior nodes with the correct index */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			p->x = p->y = p->index;
		else if (p->right != NULL)
			{
			if (p->left->y < p->right->y)
				{
				p->y = p->left->y;
				p->x = p->right->y + numTaxa - 1;
				}
			else
				{
				p->y = p->right->y;
				p->x = p->left->y + numTaxa - 1;
				}
			}
		}

	/* break the tree into pieces */
	j = 2 * numTaxa - 3;
	for (i=0; i<numTaxa-2; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTaxa-3-i]];
		q = p->anc;
		brlens[j--] = p->length;
		if (q->left == p)
			{
			if (q->anc->anc != NULL)
				brlens[j--] = q->length;
			else
				brlens[j--] = q->right->length;
			order[numTaxa-3-i] = q->right->x;
			q->right->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->right;
			else
				q->anc->right = q->right;
			}
		else
			{
			if (q->anc->anc != NULL)
				brlens[j--] = q->length;
			else
				brlens[j--] = q->left->length;
			order[numTaxa-3-i] = q->left->x;
			q->left->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->left;
			else
				q->anc->right = q->left;
			}
		}

    /* store the final two branch lengths in the right order; they have indices 0 and 1 */
    p = t->root->left;
    brlens[p->left->index] = p->left->length;
    brlens[p->right->index] = p->right->length;

    return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   StoreRTreeWithIndices: This routine will break a rooted
|      tree into an array of ints describing the structure
|      of the tree and the interior node indices, and an array
|      of doubles storing the branch lengths. The tree will be
|      destroyed in the process (the node pointers,
|      that is). However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreRTreeWithIndices (Tree *t, int *order, MrBFlt *brlens)

{
	int			i, j, k, numTaxa;
	TreeNode	*p, *q;

	extern void ShowNodes (TreeNode *p, int indent, int isRooted);

	/* find number of taxa */
	numTaxa = t->nNodes - t->nIntNodes - 1;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first two taxa */
		if (p->index > 1 && p->index < numTaxa)
			order[p->index-2] = i;
		}

	/* label the interior nodes with the correct index */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			p->x = p->y = p->index;
		else if (p->right != NULL)
			{
			if (p->left->y < p->right->y)
				{
				p->y = p->left->y;
				p->x = p->right->y + numTaxa - 1;
				}
			else
				{
				p->y = p->right->y;
				p->x = p->left->y + numTaxa - 1;
				}
			}
		}

	/* break the tree into pieces */
	j = 2 * numTaxa - 3;
    k = 2*(numTaxa - 2);
	for (i=0; i<numTaxa-2; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTaxa-3-i]];
		q = p->anc;
		brlens[j--] = p->length;
		if (q->left == p)
			{
			if (q->anc->anc != NULL)
				brlens[j--] = q->length;
			else
				brlens[j--] = q->right->length;
			order[k--] = q->right->x;
            order[k--] = q->index;
			q->right->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->right;
			else
				q->anc->right = q->right;
			}
		else
			{
			if (q->anc->anc != NULL)
				brlens[j--] = q->length;
			else
				brlens[j--] = q->left->length;
			order[k--] = q->left->x;
            order[k--] = q->index;
			q->left->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->left;
			else
				q->anc->right = q->left;
			}
		}

    /* store the final two branch lengths in the right order; they have indices 0 and 1 */
    p = t->root->left;
    order[k] = p->index;
    brlens[p->left->index] = p->left->length;
    brlens[p->right->index] = p->right->length;

    return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   StoreUPolyTopology: Same as StoreUTopology but for
|      binary polytree source.
|
--------------------------------------------------------*/
int StoreUPolyTopology (PolyTree *t, int *order)

{
	int			i, numTips;
	PolyNode	*p, *q;

    /* check if the tree is rooted on taxon 0 */
    if (t->root->left->sib->sib->index != 0)
		MovePolyCalculationRoot (t, 0);

    /* rearrange the root */
    t->root->anc = t->root->left->sib->sib;
    t->root->left->sib->sib = NULL;
    t->root->anc->left = t->root;
    t->root->anc->sib = NULL;
    t->root->anc->anc = NULL;
    t->root = t->root->anc;

	/* find number of tips */
	numTips = t->nNodes - t->nIntNodes;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first three taxa */
		if (p->index > 2 && p->index < numTips)
			order[p->index-3] = i;
		}

	/* label the interior nodes with the correct index */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL || p->anc == NULL)
			p->x = p->y = p->index;
		else
			{
			if (p->left->y < p->left->sib->y)
				{
				p->y = p->left->y;
				p->x = p->left->sib->y + numTips - 2;
				}
			else
				{
				p->y = p->left->sib->y;
				p->x = p->left->y + numTips - 2;
				}
			}
		}

    /* break the tree into pieces */
	for (i=0; i<numTips-3; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTips-4-i]];
		q = p->anc;
		if (q->left == p)
			{
			order[numTips-4-i] = q->left->sib->x;
			p->sib->anc = q->anc;
            if (q->anc->left == q)
                {
                q->anc->left = p->sib;
                p->sib->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = p->sib;
                p->sib->sib = q->sib;
                }
			}
		else
			{
			order[numTips-4-i] = q->left->x;
			q->left->anc = q->anc;
            if (q->anc->left == q)
                {
                q->anc->left = q->left;
                q->left->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = q->left;
                q->left->sib = q->sib;
                }
			}
		}

	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   StoreUPolyTree: Same as StoreUTopology but for
|      binary polytree source.
|
--------------------------------------------------------*/
int StoreUPolyTree (PolyTree *t, int *order, MrBFlt *brlens)

{
	int			i, j, numTips;
	PolyNode	*p, *q;

    /* check if the tree is rooted on taxon 0 */
    if (t->root->left->sib->sib->index != 0)
		MovePolyCalculationRoot (t, 0);

    /* rearrange the root */
    t->root->anc = t->root->left->sib->sib;
    t->root->left->sib->sib = NULL;
    t->root->anc->left = t->root;
    t->root->anc->sib = NULL;
    t->root->anc->anc = NULL;
    t->root = t->root->anc;

    /* find number of tips */
	numTips = t->nNodes - t->nIntNodes;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first three taxa */
		if (p->index > 2 && p->index < numTips)
			order[p->index-3] = i;
		}

	/* label the interior nodes with the correct index */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL || p->anc == NULL)
			p->x = p->y = p->index;
		else
			{
			if (p->left->y < p->left->sib->y)
				{
				p->y = p->left->y;
				p->x = p->left->sib->y + numTips - 2;
				}
			else
				{
				p->y = p->left->sib->y;
				p->x = p->left->y + numTips - 2;
				}
			}
		}

	/* break the tree into pieces */
    j = 2*numTips - 4;
	for (i=0; i<numTips-3; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTips-4-i]];
        assert (p->index > 2 && p->index < numTips);
        assert (p->anc->anc != NULL);
		q = p->anc;
        brlens[j--] = p->length;
        brlens[j--] = q->length;
		if (q->left == p)
			{
			order[numTips-4-i] = q->left->sib->x;
			p->sib->anc = q->anc;
            if (q->anc->left == q)
                {
                q->anc->left = p->sib;
                p->sib->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = p->sib;
                p->sib->sib = q->sib;
                }
			}
		else
			{
			order[numTips-4-i] = q->left->x;
			q->left->anc = q->anc;
            if (q->anc->left == q)
                {
                q->anc->left = q->left;
                q->left->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = q->left;
                q->left->sib = q->sib;
                }
			}
		}

    /* store last three branch lengths, index 0, 1, 2 */
    q = t->root;
    assert (q->index == 0);
    brlens[q->index] = q->length;
    q = q->left->left;
    assert (q->index == 1 || q->index == 2);
    brlens[q->index] = q->length;
    q = q->sib;
    assert (q->index == 1 || q->index == 2);
    brlens[q->index] = q->length;

    return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   StoreUTopology: This routine will break an unrooted tree
|      into an array of ints describing the structure
|      of the tree. The tree will be destroyed
|      in the process (the node pointers, that is).
|      However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreUTopology (Tree *t, int *order)

{
	int			i, numTips;
	TreeNode	*p, *q;

    /* check if the tree is rooted on taxon 0 */
    if (t->root->index != 0)
		MoveCalculationRoot (t, 0);

	/* find number of tips */
	numTips = t->nNodes - t->nIntNodes;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first three taxa */
		if (p->index > 2 && p->index < numTips)
			order[p->index-3] = i;
		}

	/* label the interior nodes with the correct index */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			p->x = p->y = p->index;
		else if (p->right != NULL)
			{
			if (p->left->y < p->right->y)
				{
				p->y = p->left->y;
				p->x = p->right->y + numTips - 2;
				}
			else
				{
				p->y = p->right->y;
				p->x = p->left->y + numTips - 2;
				}
			}
		}

	/* break the tree into pieces */
	for (i=0; i<numTips-3; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTips-4-i]];
		q = p->anc;
		if (q->left == p)
			{
			order[numTips-4-i] = q->right->x;
			q->right->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->right;
			else
				q->anc->right = q->right;
			}
		else
			{
			order[numTips-4-i] = q->left->x;
			q->left->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->left;
			else
				q->anc->right = q->left;
			}
		}

	return (NO_ERROR);
}





/*-------------------------------------------------------
|
|   StoreUTree: This routine will break an unrooted tree
|      into an array of ints describing the structure
|      of the tree and an array of doubles storing
|      the branch lengths. The tree will be
|      destroyed in the process (the node pointers,
|      that is). However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreUTree (Tree *t, int *order, MrBFlt *brlens)

{
	int			i, j, numTips;
	TreeNode	*p, *q;

	/* check if the tree is rooted on taxon 0 */
	if (t->root->index != 0)
		MoveCalculationRoot(t, 0);

	/* find number of tips */
	numTips = t->nNodes - t->nIntNodes;

	/* first get the terminal taxon positions and store
	   them in the order array. */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		/* we do not need to worry about the first three taxa */
		if (p->index > 2 && p->index < numTips)
			order[p->index-3] = i;
		}

	/* label the interior nodes with the correct index */
	for (i=0; i<t->nNodes; i++)
		{
		p = t->allDownPass[i];
		if (p->left == NULL)
			p->x = p->y = p->index;
		else if (p->right != NULL)
			{
			if (p->left->y < p->right->y)
				{
				p->y = p->left->y;
				p->x = p->right->y + numTips - 2;
				}
			else
				{
				p->y = p->right->y;
				p->x = p->left->y + numTips - 2;
				}
			}
		}

	/* break the tree into pieces */
	j = 2 * numTips - 4;
	for (i=0; i<numTips-3; i++)
		{
		/* find the next node to remove */
		p = t->allDownPass[order[numTips-4-i]];
		q = p->anc;
		brlens[j--] = p->length;
		brlens[j--] = q->length;
		if (q->left == p)
			{
			order[numTips-4-i] = q->right->x;
			q->right->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->right;
			else
				q->anc->right = q->right;
			}
		else
			{
			order[numTips-4-i] = q->left->x;
			q->left->anc = q->anc;
			if (q->anc->left == q)
				q->anc->left = q->left;
			else
				q->anc->right = q->left;
			}
		}

	/* store the final three branch lengths */
	/* we need to check the rotation of the tree to 
           store the brlens in the right order (after node index) */
	p = t->root->left;
	if (p->right->index == 2)
            {
            brlens[j--] = p->right->length;
	    brlens[j--] = p->left->length;
            }
	else
            {
            brlens[j--] = p->left->length;
            brlens[j--] = p->right->length;
            }
	brlens[j--] = p->length;

	return (NO_ERROR);
}





/* TreeLength: Calculate tree length */
MrBFlt TreeLen (Tree *t)
{
    int     i, numLenNodes;
    MrBFlt  len = 0.0;

    if (t->isRooted == NO)
        numLenNodes = t->nNodes - 1;
    else
        numLenNodes = t->nNodes - 2;

    for (i=0; i<numLenNodes; i++)
        len += t->allDownPass[i]->length;

    return len;
}





/*-------------------------------------------------------------------------------------------
|
|   Unmark: This routine will unmark a subtree rooted at p
|
---------------------------------------------------------------------------------------------*/
void Unmark (TreeNode *p)
{
    if (p != NULL)
        {
        p->marked = NO;
        Unmark (p->left);
        Unmark (p->right);
        }
}





void WriteEventTree (TreeNode *p, int chain, Param *param)

{
	int             j, nEvents;
	MrBFlt			brlen, *position, *rateMult;

	if (p != NULL)
		{
		if (p->left == NULL && p->right == NULL)
			{
			printf ("%d:%s", p->index + 1, MbPrintNum(p->length));
            if (param->paramType == P_CPPEVENTS)
				{
				nEvents = param->nEvents[2*chain+state[chain]][p->index];
				if (nEvents > 0)
					{
                    printf ("[&E %s %d: (", param->name, nEvents);
					position = param->position[2*chain+state[chain]][p->index];
					rateMult = param->rateMult[2*chain+state[chain]][p->index];
					for (j=0; j<nEvents; j++)
						{
						printf ("%s", MbPrintNum(position[j]) );
                        printf (" %s", MbPrintNum(rateMult[j]) );
						if (j != nEvents-1)
							printf (", ");
						}
					printf (")]");
					}
                else
                    printf ("[&E %s 0]", param->name);
				}
	        brlen = (GetParamSubVals (param, chain, state[chain])[p->index] + GetParamVals (param, chain, state[chain])[p->anc->index]) / 2.0;
	        printf ("[&B %s %s]", param->name, MbPrintNum(brlen));
			}
		else
			{
			if (p->anc != NULL)
				printf ("(");
			WriteEventTree(p->left, chain, param);
			printf (",");
			WriteEventTree(p->right, chain, param);
			if (p->anc != NULL)
				{				
				if (p->anc->anc != NULL)
					{
					printf ("):%s", MbPrintNum(p->length));
                    if (param->paramType == P_CPPEVENTS)
				        {
				        nEvents = param->nEvents[2*chain+state[chain]][p->index];
				        if (nEvents > 0)
					        {
                            printf ("[&E %s %d: (", param->name, nEvents);
					        position = param->position[2*chain+state[chain]][p->index];
					        rateMult = param->rateMult[2*chain+state[chain]][p->index];
					        for (j=0; j<nEvents; j++)
						        {
                                printf ("%s", MbPrintNum(position[j]) );
                                printf (" %s", MbPrintNum(rateMult[j]) );
						        if (j != nEvents-1)
							        printf (", ");
						        }
					        printf (")]");
					        }
                        else
                            printf ("[&E %s 0]", param->name);
				        }
			        brlen = (GetParamSubVals (param, chain, state[chain])[p->index] + GetParamVals (param, chain, state[chain])[p->anc->index]) / 2.0;
			        printf ("[&B %s %s]", param->name, MbPrintNum(brlen));
					}
				else
					printf(")");
				}
			}
		}
}





void WriteEventTreeToPrintString (TreeNode *p, int chain, Param *param, int printAll)

{
	char			*tempStr;
	int             i, j, nEvents, tempStrSize = TEMPSTRSIZE;
	MrBFlt			brlen, *position, *rateMult;

	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
		MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));

	if (p != NULL)
		{
		if (p->left == NULL && p->right == NULL)
			{
			SafeSprintf (&tempStr, &tempStrSize, "%d:%s", p->index + 1, MbPrintNum(p->length));
			AddToPrintString (tempStr);
			for (i=0; i<param->nSubParams; i++)
				{
                if (param->subParams[i]->paramType == P_CPPEVENTS)
					{
					nEvents = param->subParams[i]->nEvents[2*chain+state[chain]][p->index];
					if (nEvents > 0)
						{
                        SafeSprintf (&tempStr, &tempStrSize, "[&E %s %d", param->subParams[i]->name, nEvents);
						AddToPrintString (tempStr);
						position = param->subParams[i]->position[2*chain+state[chain]][p->index];
						rateMult = param->subParams[i]->rateMult[2*chain+state[chain]][p->index];
						if (printAll == YES)
                            {
                            SafeSprintf (&tempStr, &tempStrSize, ": (");
						    AddToPrintString (tempStr);
                            for (j=0; j<nEvents; j++)
							    {
						        SafeSprintf (&tempStr, &tempStrSize, "%s", MbPrintNum(position[j]));
							    AddToPrintString (tempStr);
                                SafeSprintf (&tempStr, &tempStrSize, " %s",  MbPrintNum(rateMult[j]));
							    AddToPrintString (tempStr);
							    if (j != nEvents-1)
								    AddToPrintString (",");
                                else
                                    AddToPrintString (")");
							    }
                            }
						AddToPrintString ("]");
						}
                    else
                        {
                        SafeSprintf (&tempStr, &tempStrSize, "[&E %s 0]", param->subParams[i]->name);
						AddToPrintString (tempStr);
                        }
					}
                else if (param->subParams[i]->paramType != P_CPPEVENTS)
                    {
                    /* other relaxed clock models */
                    brlen = GetParamSubVals (param->subParams[i], chain, state[chain])[p->index];
				    SafeSprintf (&tempStr, &tempStrSize, "[&B %s %s]", param->subParams[i]->name, MbPrintNum(brlen));
				    AddToPrintString (tempStr);
                    }
				}
			}
		else
			{
			if (p->anc != NULL)
				AddToPrintString ("(");
			WriteEventTreeToPrintString (p->left, chain, param, printAll);
			AddToPrintString (",");
			WriteEventTreeToPrintString (p->right, chain, param, printAll);	
			if (p->anc != NULL)
				{				
				if (p->anc->anc != NULL)
					{
                    SafeSprintf (&tempStr, &tempStrSize, "):%s", MbPrintNum(p->length));
					AddToPrintString (tempStr);
					for (i=0; i<param->nSubParams; i++)
						{
                        if (param->subParams[i]->paramType == P_CPPEVENTS)
							{
					        nEvents = param->subParams[i]->nEvents[2*chain+state[chain]][p->index];
					        if (nEvents > 0)
						        {
                                SafeSprintf (&tempStr, &tempStrSize, "[&E %s %d", param->subParams[i]->name, nEvents);
						        AddToPrintString (tempStr);
						        position = param->subParams[i]->position[2*chain+state[chain]][p->index];
						        rateMult = param->subParams[i]->rateMult[2*chain+state[chain]][p->index];
						        if (printAll == YES)
                                    {
                                    SafeSprintf (&tempStr, &tempStrSize, ": (");
						            AddToPrintString (tempStr);
                                    for (j=0; j<nEvents; j++)
							            {
						                SafeSprintf (&tempStr, &tempStrSize, "%s", MbPrintNum(position[j]));
							            AddToPrintString (tempStr);
                                        SafeSprintf (&tempStr, &tempStrSize, " %s",  MbPrintNum(rateMult[j]));
							            AddToPrintString (tempStr);
							            if (j != nEvents-1)
								            AddToPrintString (",");
                                        else
                                            AddToPrintString (")");
							            }
                                    }
						        AddToPrintString ("]");
						        }
                            else
                                {
                                SafeSprintf (&tempStr, &tempStrSize, "[&E %s 0]", param->subParams[i]->name);
						        AddToPrintString (tempStr);
                                }
					        }
                        else if (param->subParams[i]->paramType != P_CPPEVENTS)
                            {
                            /* other relaxed clock models */
                            brlen = GetParamSubVals (param->subParams[i], chain, state[chain])[p->index];
				            SafeSprintf (&tempStr, &tempStrSize, "[&B %s %s]", param->subParams[i]->name, MbPrintNum(brlen));
				            AddToPrintString (tempStr);
                            }
						}
					}
				else
					AddToPrintString(")");
				}
			}
		}
	free (tempStr);
}





void WriteEvolTree (TreeNode *p, int chain, Param *param)

{
	MrBFlt			*length;

	if (p != NULL)
		{
		length = GetParamSubVals(param, chain, state[chain]);
        if (p->left == NULL && p->right == NULL)
			{
			printf ("%d:%s", p->index + 1, MbPrintNum(length[p->index]));
			}
		else
			{
			if (p->anc != NULL)
				printf ("(");
			WriteEvolTree(p->left, chain, param);
			printf (",");
			WriteEvolTree(p->right, chain, param);
			if (p->anc != NULL)
				{				
				if (p->anc->anc != NULL)
					printf ("):%s", MbPrintNum(length[p->index]));
				else
					printf(")");
				}
			}
		}
}





void WriteTreeToPrintString (Param *param, int chain, TreeNode *p, int showBrlens, int isRooted)

{
	char			*tempStr;
	int             i, tempStrSize = TEMPSTRSIZE, nEvents;
    MrBFlt          brlen, N;

	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
		MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));

	if (p != NULL)
		{
		if (p->left == NULL && p->right == NULL)
			{
			if (showBrlens == YES)
                {
				SafeSprintf (&tempStr, &tempStrSize, "%d:%s", p->index + 1, MbPrintNum(p->length));
                }
			else
				SafeSprintf (&tempStr, &tempStrSize, "%d", p->index + 1);
			AddToPrintString (tempStr);
			if (param->paramType == P_BRLENS)
                {
                for (i=0; i<param->nSubParams; i++)
				    {
                    if (param->subParams[i]->paramType == P_CPPEVENTS)
					    {
					    nEvents = param->subParams[i]->nEvents[2*chain+state[chain]][p->index];
                        SafeSprintf (&tempStr, &tempStrSize, "[&E %s %d]", param->subParams[i]->name, nEvents);
				        AddToPrintString (tempStr);
                        }
                    brlen = GetParamSubVals (param->subParams[i], chain, state[chain])[p->index];
				    SafeSprintf (&tempStr, &tempStrSize, "[&B %s %s]", param->subParams[i]->name, MbPrintNum(brlen));
				    AddToPrintString (tempStr);
				    }
                }
            else if (param->paramType == P_SPECIESTREE && modelSettings[param->relParts[0]].popSize->nValues > 1)
                {
                N = GetParamVals (modelSettings[param->relParts[0]].popSize, chain, state[chain])[p->index];
				SafeSprintf (&tempStr, &tempStrSize, "[&N %s %s]", modelSettings[param->relParts[0]].popSize->name, MbPrintNum(N));
				AddToPrintString (tempStr);
                }
			}
		else
			{
			if (p->anc != NULL)
				AddToPrintString ("(");
			WriteTreeToPrintString (param, chain, p->left,  showBrlens, isRooted);
			if (p->anc != NULL)
				AddToPrintString (",");
			WriteTreeToPrintString (param, chain, p->right, showBrlens, isRooted);	
			if (p->anc != NULL)
				{
				if (p->anc->anc == NULL && isRooted == NO)
					{
					if (showBrlens == YES)
				        SafeSprintf (&tempStr, &tempStrSize, ",%d:%s)", p->anc->index + 1, MbPrintNum(p->length));
					else
						SafeSprintf (&tempStr, &tempStrSize, ",%d)", p->anc->index + 1);
					AddToPrintString (tempStr);
					}
                else if (p->anc->anc != NULL)
                    {
				    if (showBrlens == YES)
					    SafeSprintf (&tempStr, &tempStrSize, "):%s", MbPrintNum(p->length));
				    else
					    SafeSprintf (&tempStr, &tempStrSize, ")");
			        AddToPrintString (tempStr);
			        if (param->paramType == P_BRLENS)
                        {
                        for (i=0; i<param->nSubParams; i++)
				            {
                            if (param->subParams[i]->paramType == P_CPPEVENTS)
					            {
					            nEvents = param->subParams[i]->nEvents[2*chain+state[chain]][p->index];
                                SafeSprintf (&tempStr, &tempStrSize, "[&E %s %d]", param->subParams[i]->name, nEvents);
				                AddToPrintString (tempStr);
                                }
                            brlen = GetParamSubVals (param->subParams[i], chain, state[chain])[p->index];
				            SafeSprintf (&tempStr, &tempStrSize, "[&B %s %s]", param->subParams[i]->name, MbPrintNum(brlen));
				            AddToPrintString (tempStr);
				            }
                        }
                    else if (param->paramType == P_SPECIESTREE && modelSettings[param->relParts[0]].popSize->nValues > 1)
                        {
                        N = GetParamVals (modelSettings[param->relParts[0]].popSize, chain, state[chain])[p->index];
				        SafeSprintf (&tempStr, &tempStrSize, "[&N %s %s]", modelSettings[param->relParts[0]].popSize->name, MbPrintNum(N));
				        AddToPrintString (tempStr);
                        }
                    }
				else if (param->paramType == P_SPECIESTREE && modelSettings[param->relParts[0]].popSize->nValues > 1)
                    {
                    N = GetParamVals (modelSettings[param->relParts[0]].popSize, chain, state[chain])[p->index];
                    SafeSprintf (&tempStr, &tempStrSize, ")[&N %s %s]", modelSettings[param->relParts[0]].popSize->name, MbPrintNum(N));
			        AddToPrintString (tempStr);
                    }
                else
                    AddToPrintString(")");
				}
			}
		}
	free (tempStr);
}





/* WriteTopologyToFile: Simply write topology to file */
void WriteTopologyToFile (FILE *fp, TreeNode *p, int isRooted)

{
	if (p != NULL)
		{
		if (p->left == NULL && p->right == NULL)
			fprintf (fp, "%d", p->index + 1);
		else
			{
			if (p->anc != NULL)
				fprintf (fp, "(");
			WriteTopologyToFile (fp, p->left, isRooted);
			if (p->anc != NULL)
				fprintf (fp, ",");
			WriteTopologyToFile (fp, p->right, isRooted);	
			if (p->anc != NULL)
				{
				if (p->anc->anc == NULL && isRooted == NO)
                    fprintf (fp, ",%d", p->anc->index + 1);
			    fprintf (fp, ")");
				}
			}
		}
}
