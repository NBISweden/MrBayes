/*
 *  Best 2.2
 *
 *  This file contains the functions 
 *  for calculating the probability of 
 *  gene trees given the species tree 
 *  and the prior probability of the 
 *  species tree
 *
 *  Liang Liu
 *  Department of Statistics
 *  The Ohio State University
 *  Columbus, Ohio
 *  
 *  liuliang@stat.ohio-state.edu
 */

#include    "bayes.h"
#include    "best.h"
#include    "command.h"
#include    "mcmc.h"
#include    "model.h"
#include    "proposal.h"
#include    "utils.h"

/****************************** Local functions converted by Fredrik from BEST code *****************************/
int         CompareDepths (const void *x, const void *y);
int         CompareDoubles (const void *x, const void *y);
int         CompareNodes (const void *x, const void *y);
int         CompareNodesByX (const void *x, const void *y);
int         GetSpeciesTreeFromMinDepths (Tree* speciesTree, double *depthMatrix);
int         GetDepthMatrix(Tree * speciesTree, double *depthMatrix);
int         GetMeanDist (Tree *speciesTree, double *depthMatrix, double *mean);
int         GetMinDepthMatrix (Tree **geneTrees, int numGeneTrees, double *depthMatrix);
void        LineagesIn (TreeNode* geneTreeNode, TreeNode* speciesTreeNode);
double      LnPriorProbGeneTree (Tree *geneTree, double mu, Tree *speciesTree, double *popSizePtr);
double      LnProposalProbSpeciesTree (Tree *speciesTree, double *depthMatrix, double expRate);
void        MapGeneTreeToSpeciesTree (Tree *geneTree, Tree *speciesTree);
int         ModifyDepthMatrix (double expRate, double *depthMatrix, RandLong *seed);

/* Global BEST variables */
BitsLong    **speciesPairSets;
double      *depthMatrix;

/* Allocate variables used by best code during mcmc */
void AllocateBestChainVariables (void)
{
    int     i, j, index, numUpperTriang, nLongsNeeded;

    // Free if by mistake variables are already allocated
    if (memAllocs[ALLOC_BEST] == YES)
        FreeBestChainVariables ();

    // Allocate space for upper triangular pair sets
    numUpperTriang     = (numSpecies * (numSpecies-1)) / 2;
    nLongsNeeded       = ((numSpecies - 1) / nBitsInALong) + 1;
    speciesPairSets    = (BitsLong **) SafeCalloc (numUpperTriang, sizeof(BitsLong *));
    speciesPairSets[0] = (BitsLong *)  SafeCalloc (numUpperTriang*nLongsNeeded, sizeof(BitsLong));
    for (i=1; i<numUpperTriang; i++)
        speciesPairSets[i] = speciesPairSets[0] + i*nLongsNeeded;

    // Set upper triangular pair partitions once and for all
    index = 0;
    for (i=0; i<numSpecies; i++) {
        for (j=i+1; j<numSpecies; j++) {
            SetBit(i, speciesPairSets[index]);
            SetBit(j, speciesPairSets[index]);
            index++;
            }
        }

    /* allocate species for depthMatrix */
    depthMatrix = SafeCalloc (numUpperTriang, sizeof(double));

    memAllocs[ALLOC_BEST] = YES;
}


/** Compare function (Depth struct) for qsort */
int CompareDepths (const void *x, const void *y) {

    if ((*((Depth *)(x))).depth < (*((Depth *)(y))).depth)
        return -1;
    else if ((*((Depth *)(x))).depth > (*((Depth *)(y))).depth)
        return 1;
    else
        return 0;
}


/** Compare function (doubles) for qsort */
int CompareDoubles (const void *x, const void *y) {

    if (*((double *)(x)) < *((double *)(y)))
        return -1;
    else if (*((double *)(x)) > *((double *)(y)))
        return 1;
    else
        return 0;
}


/** Compare function (TreeNode struct) for qsort */
int CompareNodes (const void *x, const void *y) {

    if ((*((TreeNode **)(x)))->nodeDepth < (*((TreeNode**)(y)))->nodeDepth)
        return -1;
    else if ((*((TreeNode **)(x)))->nodeDepth > (*((TreeNode**)(y)))->nodeDepth)
        return 1;
    else
        return 0;
}


/** Compare function (TreeNode struct; sort by x, then by nodeDepth) for qsort */
int CompareNodesByX (const void *x, const void *y) {

    if ((*((TreeNode **)(x)))->x < (*((TreeNode**)(y)))->x)
        return -1;
    else if ((*((TreeNode **)(x)))->x > (*((TreeNode**)(y)))->x)
        return 1;
    else {
        if ((*((TreeNode **)(x)))->nodeDepth < (*((TreeNode**)(y)))->nodeDepth)
            return -1;
        else if ((*((TreeNode **)(x)))->nodeDepth > (*((TreeNode**)(y)))->nodeDepth)
            return 1;
        else
            return 0;
        }
}


/**-----------------------------------------------------------------
|
|   FillSpeciesTreeParams: Fill in species trees (start value)
|
------------------------------------------------------------------*/
int FillSpeciesTreeParams (RandLong *seed, int fromChain, int toChain)
{
    int         i, k, chn, numGeneTrees, freeBestChainVars;
    Param       *p;
    Tree        *speciesTree, **geneTrees;

    // Allocate space for global best model variables used in this function, in case they are not allocated
    if (memAllocs[ALLOC_BEST] == NO)
        {
        freeBestChainVars = YES;
        AllocateBestChainVariables();
        }
    else
        freeBestChainVars = NO;

    // Use global variable numTopologies to calculate number of gene trees
    // There is one topology for the species tree, the other ones are gene trees
    // The number of current divisions is not safe because one gene tree can have
    // several partitions, for instance if we assign the different genes on the
    // mitochondrion different substitution models, or if we assign different rates
    // to the codon site positions in a sequence
    numGeneTrees = numTopologies - 1;
    geneTrees   = (Tree **) SafeCalloc (numGeneTrees, sizeof(Tree*));

    // Build species trees for state 0
    for (chn=fromChain; chn<toChain; chn++)
        {
        for (k=0; k<numParams; k++)
            {
            p = &params[k];
            if (p->paramType == P_SPECIESTREE)
                {
                // Find species tree and gene trees
                speciesTree = GetTree(p, chn, 0);
                for (i=0; i<p->nSubParams; i++)
                    geneTrees[i] = GetTree(p->subParams[i], chn, 0);

                // Get minimum depth matrix for species tree
                GetMinDepthMatrix (geneTrees, numGeneTrees, depthMatrix);

                // Get a species tree from min depth matrix
                GetSpeciesTreeFromMinDepths(speciesTree, depthMatrix);

                assert (IsSpeciesTreeConsistent(speciesTree, chn) == YES);
                
                // Label the tips
                if (LabelTree (speciesTree, speciesNameSets[speciespartitionNum].names) == ERROR)
                    {
                    FreeBestChainVariables();
                    return (ERROR);
                    }
                }
            }
        }

    // Free gene trees
    free (geneTrees);

    // Free best model variables if appropriate
    if (freeBestChainVars == YES)
        FreeBestChainVariables();

    return (NO_ERROR);
    MrBayesPrint ("%lf", *seed); /* just because I am tired of seeing the unused parameter error msg */
}


/**-----------------------------------------------------------------
|
|   FreeBestChainVariables: Free best variables used during an mcmc
|   run.
|
------------------------------------------------------------------*/
void FreeBestChainVariables(void)
{
    if (memAllocs[ALLOC_BEST] == YES) {
        free (speciesPairSets[0]);
        free (speciesPairSets);
        speciesPairSets = NULL;
        }

    free (depthMatrix);
    depthMatrix = NULL;

    memAllocs[ALLOC_BEST] = NO;
}


/**---------------------------------------------------------------------
|
|   GetDepthMatrix:
|
|   This algorithm calculates the upper triangular depth matrix for the
|   species tree. Time complexity O(n^2).
|
|   @param      speciesTree     The species tree (in)
|   @param      depthMatrix     The minimum depth matrix, upper triangular array (out)
|   @returns    Returns ERROR or NO_ERROR
----------------------------------------------------------------------*/
int GetDepthMatrix (Tree *speciesTree, double *depthMatrix) {

    int         i, left, right, numUpperTriang, index, nLongsNeeded, freeBitsets;
    double      maxDepth;
    TreeNode    *p;

    // Make sure we have bitfields allocated and set
    if (speciesTree->bitsets == NULL)
        {
        AllocateTreePartitions(speciesTree);
        freeBitsets = YES;
        }
    else
        {
        ResetTreePartitions(speciesTree);   // just in case
        freeBitsets = NO;
        }

    // Calculate number of values in the upper triangular matrix
    numUpperTriang = numSpecies * (numSpecies - 1) / 2;

    // Number of longs needed in a bitfield representing a species set
    nLongsNeeded   = ((numSpecies -1) / nBitsInALong) + 1;

    // Set all cells to max
    maxDepth = speciesTree->root->left->nodeDepth;  // root depth
    for (i=0; i<numUpperTriang; i++)
        depthMatrix[i] = maxDepth;

    // Loop over interior nodes
    for (i=0; i<speciesTree->nIntNodes; i++)
        {
        p = speciesTree->intDownPass[i];
        for (left = FirstTaxonInPartition(p->left->partition, nLongsNeeded); left < numSpecies; left = NextTaxonInPartition(left, p->left->partition, nLongsNeeded))
            {
            for (right = FirstTaxonInPartition(p->right->partition, nLongsNeeded); right < numSpecies; right = NextTaxonInPartition(right, p->right->partition, nLongsNeeded))
                {
                index = UpperTriangIndex(left, right, numSpecies);
                depthMatrix[index] = p->nodeDepth;
                }
            }
        }

    // Free partitions if appropriate
    if (freeBitsets == YES)
        FreeTreePartitions(speciesTree);

    return (NO_ERROR);
}


/**---------------------------------------------------------------------
|
|   GetMeanDist
|
|   This algorithm calculates the mean distance between a distance matrix
|   and the minimum depths that define a species tree.
|
|   @param      speciesTree     The species tree (in)
|   @param      minDepthMatrix  The minimum depth matrix, upper triangular array (in)
|   @param      mean            The mean distance (out)
|   @returns    Returns ERROR or NO_ERROR
----------------------------------------------------------------------*/
int GetMeanDist (Tree *speciesTree, double *minDepthMatrix, double *mean) {

    int         i, left, right, index, nLongsNeeded, freeBitsets;
    double      dist, minDist=0.0, distSum;
    TreeNode    *p;

    // Make sure we have bitfields allocated and set
    if (speciesTree->bitsets == NULL)
        {
        AllocateTreePartitions(speciesTree);
        freeBitsets = YES;
        }
    else
        {
        ResetTreePartitions(speciesTree);   // just in case
        freeBitsets = NO;
        }

    // Number of longs needed in a bitfield representing a species set
    nLongsNeeded   = ((numSpecies -1) / nBitsInALong) + 1;

    // Loop over interior nodes
    distSum = 0.0;
    for (i=0; i<speciesTree->nIntNodes; i++)
        {
        p = speciesTree->intDownPass[i];
        p->x = 0;
        while ((left=FirstTaxonInPartition(p->left->partition, nLongsNeeded)) < numSpecies)
            {
            while ((right=FirstTaxonInPartition(p->right->partition, nLongsNeeded)) < numSpecies)
                {
                p->x++;
                index = UpperTriangIndex(left, right, numSpecies);
                dist = depthMatrix[index] - p->nodeDepth;
                if (p->x == 1)
                    minDist = dist;
                else if (dist < minDist)
                    minDist = dist;
                ClearBit(right, p->right->partition);
                }
            ClearBit(left, p->left->partition);
            }
        distSum += minDist;
        }

    (*mean) = distSum / speciesTree->nIntNodes;

    // Reset partitions that were destroyed above or free partitions, as appropriate
    if (freeBitsets == YES)
        FreeTreePartitions(speciesTree);
    else
        ResetTreePartitions(speciesTree);

    return (NO_ERROR);
    MrBayesPrint ("%lf", *minDepthMatrix); /* just because I am tired of seeing the unused parameter error msg */
}


/**---------------------------------------------------------------------
|
|   GetMinDepthMatrix: converted from GetMinDists.
|
|   This algorithm scans the gene trees and calculates the minimum depth
|   (height) separating species across gene trees. The complexity of the
|   original algorithm was O(mn^3), where m is the number of gene trees and
|   n is the number of taxa in each gene tree. I think this algorithm has
|   complexity that is better on average, but the difference is small.
|
|   I have rewritten the algorithm also to show alternative techniques that
|   could be used in this and other BEST algorithms.
|
|   @param      geneTrees       The gene trees (in)
|   @param      depthMatrix     The minimum depth matrix, upper triangular array (out)
|   @returns    Returns ERROR or NO_ERROR
----------------------------------------------------------------------*/
int GetMinDepthMatrix (Tree **geneTrees, int numGeneTrees, double *depthMatrix) {

    int         i, j, w, nLongsNeeded, numUpperTriang, index, trace=0;
    double      maxDepth;
    TreeNode    *p;
    BitsLong    **speciesSets;

    // Allocate space for species partitions
    nLongsNeeded   = ((numSpecies -1) / nBitsInALong) + 1;   // number of longs needed in a bitfield representing a species set
    speciesSets    = (BitsLong **) SafeCalloc ((2*(size_t)numLocalTaxa-1), sizeof(BitsLong *));
    speciesSets[0] = (BitsLong *)  SafeCalloc ((2*(size_t)numLocalTaxa-1)*nLongsNeeded, sizeof(BitsLong));
    for (i=1; i<2*numLocalTaxa-1; i++)
        speciesSets[i] = speciesSets[0] + i*nLongsNeeded;

    // Set tip species partitions once and for all
    for (i=0; i<numLocalTaxa; i++)
        SetBit(speciespartitionId[i][speciespartitionNum]-1, speciesSets[i]);

    // Set initial max depth for upper triangular matrix
    numUpperTriang = (numSpecies * (numSpecies - 1)) / 2;
    maxDepth       = geneTrees[0]->root->left->nodeDepth;
    for (i=0; i<numUpperTriang; i++)
        depthMatrix[i] = maxDepth;

    // Now we are ready to cycle over gene trees
    for (w=0; w<numGeneTrees; w++)
        {
        if (trace) {
            printf ("\nGene %d\n",w);
            ShowTree(geneTrees[w]);
            }

        // Set species sets for interior nodes. O(n)
        for (i=0; i<geneTrees[w]->nIntNodes; i++) {
            p = geneTrees[w]->intDownPass[i];
            for (j=0; j<nLongsNeeded; j++)
                speciesSets[p->index][j] = speciesSets[p->left->index][j] | speciesSets[p->right->index][j];       
            }

        // Now order the interior nodes in terms of node depth. We rely on the fact that the
        // ordered sequence is a valid downpass sequence. O(log n).
        qsort((void *)(geneTrees[w]->intDownPass), (size_t)(geneTrees[w]->nIntNodes), sizeof(TreeNode *), CompareNodes);

        // Finally find the minimum for each cell in the upper triangular matrix
        // This is the time critical step with complexity O(n^3) in the simplest
        // algorithm version. This algorithm should do a little better in most cases.
        for (i=0; i<numUpperTriang; i++) {
            
            // Find shallowest node that has the pair
            for (j=0; j<geneTrees[w]->nIntNodes; j++) {
                p = geneTrees[w]->intDownPass[j];
                
                // Because nodes are ordered in time, if this test is true then we cannot beat the minimum
                if (p->nodeDepth > depthMatrix[i])
                    break;

                // Check whether the node is a candidate minimum for the species pair
                // If the test is true, we know from the test above that p->nodeDepth is
                // either a tie or the new minimum
                if (IsPartNested(speciesPairSets[i], speciesSets[p->index], nLongsNeeded) == YES) {
                    depthMatrix[i] = p->nodeDepth;
                    break;
                    }
                }
            }
        }   // Next gene tree

    if (trace)
        {
        index = 0;
        printf ("Mindepth matrix\n");
        for (i=0;i<numSpecies;i++) {
            for (j=0; j<i; j++)
                printf ("         ");
            for (j=i+1;j<numSpecies;j++) {
                printf ("%.6f ",depthMatrix[index]);
                index++;
                }
            printf ("\n");
            }
        printf ("\n");
        }

    free (speciesSets[0]);
    free (speciesSets);

    return (NO_ERROR);
}


/**---------------------------------------------------------------------
|
|   GetSpeciesTreeFromMinDepths: converted from GetConstraints, Startsptree,
|   and MaximumTree.
|
|   This is a clustering algorithm based on minimum depths for species pairs.
|   It reduces an n choose 2 upper triangular min depth matrix to an array
|   of n-1 node depths, which fit onto a tree.
|
|   @param      speciesTree     The species tree to be filled  (out)
|   @param      depthMatrix     The min depth matrix, upper triangular array (in)
|   @returns    Returns NO_ERROR if success, ERROR if negative brlens occur
----------------------------------------------------------------------*/
int GetSpeciesTreeFromMinDepths (Tree* speciesTree, double *depthMatrix) {

    int         i, j, numUpperTriang, nLongsNeeded, index, nextNodeIndex;
    Depth       *minDepth;
    PolyTree    *polyTree;
    PolyNode    *p, *q, *r, *u, *qPrev, *rPrev;

    nLongsNeeded    = ((numSpecies - 1) / nBitsInALong) + 1;
    numUpperTriang  = numSpecies*(numSpecies - 1) / 2;
    minDepth        = (Depth *) SafeCalloc (numUpperTriang, sizeof(Depth));

    // Convert depthMatrix to an array of Depth structs
    index = 0;
    for (i=0; i<numSpecies; i++) {
        for (j=i+1; j<numSpecies; j++) {
            minDepth[index].depth   = depthMatrix[index];
            minDepth[index].pairSet = speciesPairSets[index];
            index++;
            }
        }

    // Sort the array of distance structs (O(log n^2))
    qsort((void *)(minDepth), (size_t)(numUpperTriang), sizeof(Depth), CompareDepths);

    // The algorithm below reduces the upper triangular matrix (n choose 2) to an n-1
    // array in O(n^2log(n)) time. We build the tree at the same time, since we can
    // find included pairs in the tree in log(n) time. We use a polytomous tree for this.
    
    // Allocate space for polytomous tree and set up partitions
    polyTree = AllocatePolyTree(numSpecies);
    AllocatePolyTreePartitions(polyTree);

    // Build initial tree (a bush)
    polyTree->isRooted = YES;
    polyTree->isClock = YES;
    polyTree->root = &polyTree->nodes[2*numSpecies-2];
    for (i=0; i<numSpecies; i++) {
        p = &polyTree->nodes[i];
        p->index = i;
        p->depth = 0.0;
        p->left = NULL;
        if (i<numSpecies-1)
            p->sib = &polyTree->nodes[i+1];
        else
            p->sib = NULL;
        p->anc = polyTree->root;
        }
    p = polyTree->root;
    p->index = 2*numSpecies - 2;
    p->left = &polyTree->nodes[0];
    p->sib = NULL;
    p->anc = NULL;
    p->depth = -1.0;
    polyTree->nNodes = numSpecies + 1;
    polyTree->nIntNodes = 1;
    GetPolyDownPass(polyTree);
    ResetPolyTreePartitions(polyTree);      /* set bitsets (partitions) for initial tree */

    // Resolve bush using sorted depth structs
    nextNodeIndex = numSpecies;
    for (i=0; i<numUpperTriang; i++) {
            
        // Find tip corresponding to first taxon in pair
        p = &polyTree->nodes[FirstTaxonInPartition(minDepth[i].pairSet, nLongsNeeded)];
        
        // Descend tree until we find a node within which the pair set is nested
        do  {
            p = p->anc;
            }
        while (!IsPartNested(minDepth[i].pairSet, p->partition, nLongsNeeded));

        if (p->left->sib->sib != NULL) {
            // This node is still a polytomy
            
            // Find left and right descendants of new node
            qPrev = NULL;
            for (q=p->left; IsSectionEmpty(q->partition, minDepth[i].pairSet, nLongsNeeded); q=q->sib)
                qPrev = q;
            rPrev = q;
            for (r=q->sib;  IsSectionEmpty(r->partition, minDepth[i].pairSet, nLongsNeeded); r=r->sib)
                rPrev = r;
            
            // Introduce the new node
            u = &polyTree->nodes[nextNodeIndex];
            u->index = nextNodeIndex;
            nextNodeIndex++;
            polyTree->nIntNodes++;
            polyTree->nNodes++;
            u->left = q;
            u->anc = p;
            if (p->left == q)
                p->left = u;
            else
                qPrev->sib = u;
            // former upstream sibling to r should point to r->sib
            if (rPrev == q)
                u->sib = r->sib;
            else
                rPrev->sib = r->sib;
            if (q->sib == r)
                u->sib = r->sib;
            else
                u->sib = q->sib;
            u->depth = minDepth[i].depth;   // because minDepth structs are sorted, we know this is the min depth
            assert (u->depth > 0.0);

            // Create new taxon set with bitfield operations
            for (j=0; j<nLongsNeeded; j++)
                u->partition[j] = q->partition[j] | r->partition[j];

            // Patch the tree together with the new node added
            q->sib  = r;
            r->sib = NULL;
            q->anc = u;
            r->anc = u;
            }
        else if (p == polyTree->root && p->depth < 0.0) {
            // This is the first time we hit the root of the tree && it is resolved
            p->depth = minDepth[i].depth;
            assert (p->depth > 0.0);
            }
        // other cases should not be added to tree
        }

    // Make sure we have a complete species tree
    assert (polyTree->nIntNodes == numSpecies - 1);

    // Set traversal sequences
    GetPolyDownPass(polyTree);

    // Set branch lengths from node depths (not done automatically for us)
    // Make sure all branch lengths are nonnegative (we can have 0.0 brlens, they
    // should not be problematic in a species tree; they occur when there are
    // ties in the min depth matrix that have not been modified by the move)
    for (i=0; i<polyTree->nNodes; i++) {
        p = polyTree->allDownPass[i];
        if (p->anc == NULL)
            p->length = 0.0;
        else
            p->length = p->anc->depth - p->depth;
        if (p->length < 0.0) {
            FreePolyTree(polyTree);
            free (minDepth);
            return (ERROR); 
            }           
        }

    // Copy to species tree from polytomous tree
    CopyToSpeciesTreeFromPolyTree (speciesTree, polyTree);

    // Free locally allocated variables
    FreePolyTree(polyTree);
    free (minDepth);

    return(NO_ERROR);
}


/**---------------------------------------------------------------------------------------
|
|   IsSpeciesTreeConsistent: Called when user tries to set a species tree or when
|      attempting to use a species tree from a check point as starting value.
|
-----------------------------------------------------------------------------------------*/
int IsSpeciesTreeConsistent (Tree *speciesTree, int chain)
{
    int     i, answer, numGeneTrees, numUpperTriang, freeBestVars;
    double  *speciesTreeDepthMatrix;
    Tree    **geneTrees;

    freeBestVars = NO;
    if (memAllocs[ALLOC_BEST] == NO)
        {
        AllocateBestChainVariables();
        freeBestVars = YES;
        }

    numGeneTrees = numTrees - 1;
    geneTrees = (Tree **) SafeCalloc (numGeneTrees, sizeof(Tree *));
    for (i=0; i<numTrees-1; i++)
        geneTrees[i] = GetTreeFromIndex(i, chain, state[chain]);

    numUpperTriang = numSpecies * (numSpecies - 1) / 2;
    speciesTreeDepthMatrix = (double *) SafeCalloc (numUpperTriang, sizeof(double));

    GetMinDepthMatrix(geneTrees, numGeneTrees, depthMatrix);
    GetDepthMatrix(speciesTree, speciesTreeDepthMatrix);

    for (i=0; i<numUpperTriang; i++)
        {
        if (depthMatrix[i] < speciesTreeDepthMatrix[i])
            break;
        }

    if (i == numUpperTriang)
        answer = YES;
    else
        answer = NO;

    if (answer == NO)
        ShowNodes(speciesTree->root, 0, YES);

    if (freeBestVars == YES)
        FreeBestChainVariables();

    free (speciesTreeDepthMatrix);
    free (geneTrees);

    return answer;
}


/**---------------------------------------------------------------------------------------
|
|   LineagesIn: Recursive function to get number of gene tree lineages coming into each
|      branch of the species tree (in ->x of speciestree nodes). We also assign each gene
|      tree lineage to the corresponding species tree lineage (in ->x of genetree nodes).
|      Finally, number of coalescent events is recorded (in ->y of speciestree nodes).
|      Time complexity is O(n).
|
-----------------------------------------------------------------------------------------*/
void LineagesIn (TreeNode *geneTreeNode, TreeNode *speciesTreeNode)
{
    int nLongsNeeded;
    
    if (geneTreeNode->nodeDepth < speciesTreeNode->nodeDepth) {
        // climb up species tree
        if (speciesTreeNode->left == NULL) {
            assert (geneTreeNode->left == NULL);
            speciesTreeNode->x++;
            }
        else {
            nLongsNeeded = (numSpecies - 1) / nBitsInALong + 1;
            speciesTreeNode->x++;
            if (IsPartNested(geneTreeNode->partition, speciesTreeNode->left->partition, nLongsNeeded) == YES)
                LineagesIn (geneTreeNode, speciesTreeNode->left);
            else if (IsPartNested(geneTreeNode->partition, speciesTreeNode->right->partition, nLongsNeeded) == YES)
                LineagesIn (geneTreeNode, speciesTreeNode->right);
            }
        }
    else {
        // climb up gene tree
        if (geneTreeNode->left != NULL)
            LineagesIn(geneTreeNode->left, speciesTreeNode);
        if (geneTreeNode->right != NULL)
            LineagesIn(geneTreeNode->right, speciesTreeNode);
        if (geneTreeNode->left == NULL) {
            speciesTreeNode->x++;
            assert (speciesTreeNode->left == NULL);
            }
        else {
            speciesTreeNode->y++;
            }
        geneTreeNode->x = speciesTreeNode->index;
        }
}


/**-----------------------------------------------------------------
|
|   LnSpeciesTreeProb: Wrapper for LnJointGeneTreeSpeciesTreePr to
|       free calling functions from retrieving gene and species trees.
|
------------------------------------------------------------------*/
double LnSpeciesTreeProb(int chain)
{
    int         i, numGeneTrees;
    double      lnProb;
    Tree        **geneTrees, *speciesTree;
    ModelInfo   *m;

    m = &modelSettings[0];

    speciesTree = GetTree(m->speciesTree, chain, state[chain]);

    numGeneTrees = m->speciesTree->nSubParams;
    geneTrees = (Tree **) SafeCalloc (numGeneTrees, sizeof(Tree *));

    for (i=0; i<m->speciesTree->nSubParams; i++)
        geneTrees[i] = GetTree(m->speciesTree->subParams[i], chain, state[chain]);

    lnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, speciesTree, chain);

    free (geneTrees);

    return lnProb;
}


/**-----------------------------------------------------------------
|
|   LnJointGeneTreeSpeciesTreePr: Converted from LnJointGenetreePr,
|   SPLogLike, SPLogPrior.
|
|   In this function we calculate the entire probability of the species
|   tree, including its probability given its priors, and the probability
|   of the gene trees given the species tree.
|
------------------------------------------------------------------*/
double LnJointGeneTreeSpeciesTreePr(Tree **geneTrees, int numGeneTrees, Tree *speciesTree, int chain)
{
    double      lnPrior, lnLike, clockRate, mu, *popSizePtr, sR, eR, sF;
    int         i;
    ModelInfo   *m;
    ModelParams *mp;

    // Get model info for species tree
    m = &modelSettings[speciesTree->relParts[0]];

    // Get model params for species tree
    mp = &modelParams[speciesTree->relParts[0]];

    // Get popSize ptr
    popSizePtr = GetParamVals(m->popSize, chain, state[chain]);

    // Get clock rate
    if (speciesTree->isCalibrated == YES)
        clockRate = *GetParamVals(m->clockRate, chain, state[chain]);
    else
        clockRate = 1.0;

    // Calculate probability of gene trees given species tree
    // ShowNodes(speciesTree->root, 0, YES);
    lnLike = 0.0;
    mu = clockRate;
    for (i=0; i<numGeneTrees; i++) {
        lnLike += LnPriorProbGeneTree(geneTrees[i], mu, speciesTree, popSizePtr);
        }

    // Calculate probability of species tree given its priors
    if (strcmp(mp->speciesTreeBrlensPr, "Birthdeath") == 0) {
        sR = *GetParamVals(m->speciationRates, chain, state[chain]);
        eR = *GetParamVals(m->extinctionRates, chain, state[chain]);
        sF = mp->sampleProb;
        lnPrior = 0.0;
        LnBirthDeathPriorPr(speciesTree, clockRate, &lnPrior, sR, eR, mp->sampleStrat, sF);
        }
    else
        lnPrior = 0.0;

    // The population size is taken care of elsewhere

    return lnLike + lnPrior;
}


/**-----------------------------------------------------------------
|
|   LnPriorProbGeneTree: Calculate the prior probability of a gene
|   tree.
|
------------------------------------------------------------------*/
double LnPriorProbGeneTree (Tree *geneTree, double mu, Tree *speciesTree, double *popSizePtr)
{ 
    int         i, k, index, nEvents, trace=0;
    double      N, lnProb, ploidyFactor, theta, timeInterval;
    TreeNode    *p, *q=NULL, *r;
    ModelParams *mp;

    // Get model params
    mp = &modelParams[speciesTree->relParts[0]];

    // Find ploidy setting
    if (strcmp(mp->ploidy, "Diploid") == 0)
        ploidyFactor = 4.0;
    else if (strcmp(mp->ploidy, "Haploid") == 0)
        ploidyFactor = 2.0;
    else /* if (strcmp(mp->ploidy, "Zlinked") == 0) */
        ploidyFactor = 3.0;

    // Initialize species tree with theta in d
    for (i=0; i<speciesTree->nNodes-1; i++) {
        p = speciesTree->allDownPass[i];
        if (strcmp(mp->popVarPr, "Equal") != 0)
            N = popSizePtr[p->index];
        else
            N = popSizePtr[0];
        p->d = ploidyFactor * N * mu;
        }
    
    // Map gene tree to species tree
    MapGeneTreeToSpeciesTree(geneTree, speciesTree);

    // Sort gene tree interior nodes first by speciestree branch on which they coalesce, then in time order
    qsort((void *)(geneTree->intDownPass), (size_t)(geneTree->nIntNodes), sizeof(TreeNode *), CompareNodesByX);

    // Debug output of qsort result
    if (trace) {
        printf ("index -- x -- nodeDepth for gene tree\n");
        for (i=0; i<geneTree->nIntNodes; i++)
            printf ("%d -- %d -- %e\n", geneTree->intDownPass[i]->index, geneTree->intDownPass[i]->x, geneTree->intDownPass[i]->nodeDepth);
        }

    // Now calculate probability after making sure species tree nodes appear in index order
    // (the order does not have to be a correct downpass sequence)
    for (i=0; i<speciesTree->memNodes; i++)
        {
        p = &(speciesTree->nodes[i]);
        speciesTree->allDownPass[p->index] = p;
        }
    index = 0;
    lnProb = 0.0;
    for (i=0; i<speciesTree->nNodes-1; i++) {

        p = speciesTree->allDownPass[i];

        // Get theta
        theta = p->d;

        // Get number of events
        nEvents = p->y;

        // Calculate probability
        lnProb += nEvents * log (2.0 / theta);

        for (k=p->x; k > p->x - p->y; k--) {

            q = geneTree->intDownPass[index];
            assert (q->x == p->index);

            if (k == p->x)
                timeInterval = q->nodeDepth - p->nodeDepth;
            else {
                r = geneTree->intDownPass[index-1];
                timeInterval = q->nodeDepth - r->nodeDepth;
            }

            lnProb -= (k * (k - 1) * timeInterval) / theta;
            index++;
            }

        if (p->x - p->y > 1) {

            if (nEvents == 0)
                timeInterval = p->anc->nodeDepth - p->nodeDepth;
            else
                timeInterval = p->anc->nodeDepth - q->nodeDepth;

            assert (p->anc->anc != NULL);
            assert (timeInterval >= 0.0);

            k = p->x - p->y;
            lnProb -= (k * (k - 1) * timeInterval) / theta;
            }
        }

    // Restore downpass sequences (probably not necessary for gene tree, but may be if some
    // code relies on intDownPass and allDownPass to be in same order)
    GetDownPass(speciesTree);
    GetDownPass(geneTree);

    // Free space
    FreeTreePartitions(speciesTree);
    FreeTreePartitions(geneTree);

    return lnProb;
}


/**---------------------------------------------------------------------
|
|   LnProposalProbSpeciesTree:
|
|   This algorithm calculates the probability of proposing a particular
|   species tree given a distance matrix modified using a scheme based on
|   truncated exponential distributions with rate expRate.
|
|   @param      speciesTree     The species tree (in)
|   @param      depthMatrix     The minimum depth matrix, upper triangular array (in)
|   @param      expRate         Rate of truncated exponential distribution
|   @returns    Returns probability of proposing the species tree
----------------------------------------------------------------------*/
double LnProposalProbSpeciesTree (Tree *speciesTree, double *depthMatrix, double expRate)
{
    int         i, left, right, index, nLongsNeeded, freeBitsets;
    double      dist, normConst=1.0, negLambdaX=0.0, eNegLambdaX, density, prob,
                sumDensRatio, prodProb, lnProb;
    TreeNode    *p;

    // Make sure we have bitfields allocated and set
    if (speciesTree->bitsets == NULL)
        freeBitsets = YES;
    else
        freeBitsets = NO;
    AllocateTreePartitions(speciesTree);

    // Number of longs needed in a bitfield representing a species set
    nLongsNeeded   = ((numSpecies -1) / nBitsInALong) + 1;

    // Loop over interior nodes
    lnProb = 0.0;
    for (i=0; i<speciesTree->nIntNodes; i++)
        {
        p = speciesTree->intDownPass[i];
        p->x = 0;
        sumDensRatio = 0.0;
        prodProb = 1.0;
        for (left = FirstTaxonInPartition(p->left->partition, nLongsNeeded); left < numSpecies; left = NextTaxonInPartition(left, p->left->partition, nLongsNeeded))
            {
            for (right = FirstTaxonInPartition(p->right->partition, nLongsNeeded); right < numSpecies; right = NextTaxonInPartition(right, p->right->partition, nLongsNeeded))
                {
                p->x++;
                index         = UpperTriangIndex(left, right, numSpecies);  assert (index < numSpecies*(numSpecies - 1) / 2);
                dist          = depthMatrix[index] - p->nodeDepth;          // distance between depth matrix entry and actual species-tree node
                normConst     = 1.0 - exp(-expRate * depthMatrix[index]);   // normalization constant because of truncation of exp distribution
                negLambdaX    = - expRate * dist;
                eNegLambdaX   = exp(negLambdaX);
                density       = expRate * eNegLambdaX / normConst;      // density for x == dist, f(dist)
                prob          = (1.0 - eNegLambdaX) / normConst;        // cumulative prob for x <= dist, F(dist)
                sumDensRatio += density / prob;                         // warning: if dist==0, prob is ZERO!
                prodProb     *= prob;
                }
            }
        if (p->x == 1)
            lnProb += log(expRate) + negLambdaX - log(normConst);
        else
            lnProb += log(sumDensRatio * prodProb);
        }

    // to avoid lnProposalProb is NaN at initial steps
    if (lnProb != lnProb)  lnProb = 0.0;
    
    // Free partitions if appropriate
    if (freeBitsets == YES)
        FreeTreePartitions(speciesTree);

    return (lnProb);
}


/**-----------------------------------------------------------------
|
|   MapGeneTreeToSpeciesTree: Fold gene tree into species tree. We
|      are going to use ->x of gene tree to give index of the
|      corresponding node in the species tree. ->x in the species
|      tree will give the number of lineages into the corresponding
|      branch, and ->y will give the number of coalescent events on
|      that branch.
|
------------------------------------------------------------------*/
void MapGeneTreeToSpeciesTree (Tree *geneTree, Tree *speciesTree)
{ 
    int         i, j, nLongsNeeded, trace=0;
    TreeNode    *p;

    // Initialize species partitions for both gene tree and species tree
    // This will set the partitions to reflect the partitions in the tree itself,
    // which is OK for the species tree, but we want the gene tree partitions to
    // reflect the species partitions and not the gene partitions, so we need to
    // set them here
    AllocateTreePartitions(geneTree);
    AllocateTreePartitions(speciesTree);
    nLongsNeeded = (numSpecies - 1) / nBitsInALong + 1;
    for (i=0; i<geneTree->nNodes-1; i++) {
        p = geneTree->allDownPass[i];
        ClearBits(p->partition, nLongsNeeded);
        if (p->left == NULL)
            SetBit(speciespartitionId[p->index][speciespartitionNum]-1, p->partition);
        else {
            for (j=0; j<nLongsNeeded; j++)
                p->partition[j] = p->left->partition[j] | p->right->partition[j];
            }
        }
    // Species tree partitions already set by call to AllocateTreePartitions

    // Reset ->x and ->y of species tree (->x of gene tree does not need to be initialized)
    for (i=0; i<speciesTree->nNodes; i++)
        {
        p = speciesTree->allDownPass[i];
        p->x = 0;
        p->y = 0;
        }

    // Call recursive function to match gene tree and species tree
    LineagesIn(geneTree->root->left, speciesTree->root->left);

    if (trace) {
        printf ("index -- x -- y   for species tree\n");
        for (i=0; i<speciesTree->nNodes-1; i++)
            printf ("%-2d -- %d -- %d\n", speciesTree->allDownPass[i]->index, speciesTree->allDownPass[i]->x, speciesTree->allDownPass[i]->y);
        }

    if (trace) {
        printf ("index -- x -- nodeDepth for gene tree\n");
        for (i=0; i<geneTree->nIntNodes; i++)
            printf ("%-2d -- %d -- %e\n", geneTree->intDownPass[i]->index, geneTree->intDownPass[i]->x, geneTree->intDownPass[i]->nodeDepth);
        }

    // Free space
    FreeTreePartitions(speciesTree);
    FreeTreePartitions(geneTree);
}


/**---------------------------------------------------------------------
|
|   ModifyDepthMatrix:
|
|   This algorithm uses a truncated exponential distribution to modify
|   a depth matrix.
|
|   @param      expRate         The rate of the exponential distribution (in)
|   @param      depthMatrix     The minimum depth matrix to be modified, upper triangular array (in/out)
|   @param      seed            Pointer to seed for random number generator (in/ut)
|   @returns    Returns ERROR or NO_ERROR
----------------------------------------------------------------------*/
int ModifyDepthMatrix (double expRate, double *depthMatrix, RandLong *seed)
{
    int     i, numUpperTriang;
    double  u, interval, delta;

    numUpperTriang = numSpecies * (numSpecies - 1) / 2;
    for (i=0; i<numUpperTriang; i++)
        {
        interval = depthMatrix[i];
        u = RandomNumber(seed);
        delta = log (1.0 - u*(1.0 - exp(-expRate*interval))) / (-expRate);
        assert (delta <= interval);
        depthMatrix[i] -= delta;
        }

    return (NO_ERROR);
}


/**-----------------------------------------------------------------
|
|   Move_GeneTree1: Propose a new gene tree using ExtSPRClock
|
|   @param param            The parameter (gene tree) to change
|   @param chain            The chain number
|   @param seed             Pointer to the seed of the random number gen.
|   @param lnPriorRatio     Pointer to the log prior ratio (out)
|   @param lnProposalRatio  Pointer to the log proposal (Hastings) ratio (out)
|   @param mvp              Pointer to tuning parameter(s)
------------------------------------------------------------------*/
int Move_GeneTree1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             i, numGeneTrees, numUpperTriang;
    double          newLnProb, oldLnProb, backwardLnProposalProb, forwardLnProposalProb,
                    *oldMinDepths, *modMinDepths, forwardLambda, backwardLambda, mean;
    Tree            *newSpeciesTree, *oldSpeciesTree, **geneTrees;
    ModelInfo       *m;

    // Calculate number of gene trees
    numGeneTrees = numTopologies - 1;

    // Get model settings
    m = &modelSettings[param->relParts[0]];

    // Get species tree (this trick is possible because we always copy tree params)
    newSpeciesTree = GetTree (m->speciesTree, chain, state[chain]);
    oldSpeciesTree = GetTree (m->speciesTree, chain, state[chain] ^ 1);

    // Get gene trees
    geneTrees = (Tree **) SafeCalloc (2*numGeneTrees, sizeof(Tree *));
    for (i=0; i<m->speciesTree->nSubParams; i++) {
        geneTrees[i] = GetTree(m->speciesTree->subParams[i], chain, state[chain]);
        }

    // Allocate space for depth matrix copy
    numUpperTriang = numSpecies * (numSpecies - 1) / 2;
    oldMinDepths   = (double *) SafeCalloc (2*numUpperTriang, sizeof(double));
    modMinDepths   = oldMinDepths + numUpperTriang;

    // Get min depth matrix for old gene trees
    GetMinDepthMatrix(geneTrees, numTopologies-1, depthMatrix);

    // Save a copy
    for (i=0; i<numUpperTriang; i++)
        oldMinDepths[i] = depthMatrix[i];

    // Get forward lambda
    GetMeanDist(oldSpeciesTree, depthMatrix, &mean);
    // if (mean < 1E-6) mean = 1E-6;
    forwardLambda = 1.0 / mean;

    // Calculate joint probability of old gene trees and old species tree
    oldLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, oldSpeciesTree, chain);

    // Modify the picked gene tree using code from a regular MrBayes move
    Move_ExtSPRClock(param, chain, seed, lnPriorRatio, lnProposalRatio, mvp);

    // Update the min depth matrix
    GetMinDepthMatrix(geneTrees, numTopologies-1, depthMatrix);

    // Copy the min depth matrix
    for (i=0; i<numUpperTriang; i++)
        modMinDepths[i] = depthMatrix[i];

    // Modify the min depth matrix
    ModifyDepthMatrix (forwardLambda, modMinDepths, seed);

    // Get a new species tree
    if (GetSpeciesTreeFromMinDepths (newSpeciesTree, modMinDepths) == ERROR) {
        abortMove = YES;
        free (geneTrees);
        free (oldMinDepths);
        return (NO_ERROR);
        }
    
    // Calculate joint probability of new gene trees and new species tree
    newLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, newSpeciesTree, chain);

    // Get backward lambda
    GetMeanDist(newSpeciesTree, depthMatrix, &mean);
    // if (mean < 1E-6) mean = 1E-6;
    backwardLambda = 1.0 / mean;

    // Get proposal probability of old species tree
    backwardLnProposalProb = LnProposalProbSpeciesTree (oldSpeciesTree, oldMinDepths, backwardLambda);

    // Get proposal probability of new species tree
    forwardLnProposalProb = LnProposalProbSpeciesTree (newSpeciesTree, depthMatrix, forwardLambda);

    // Update prior ratio taking species tree into account
    (*lnPriorRatio) += (newLnProb - oldLnProb);
        
    // Update proposal ratio based on this move
    (*lnProposalRatio) += (backwardLnProposalProb - forwardLnProposalProb);

    // Free allocated memory
    free (geneTrees);
    free (oldMinDepths);

    return (NO_ERROR);
}


/**-----------------------------------------------------------------
|
|   Move_GeneTree2: Propose a new gene tree using NNIClock
|
|   @param param            The parameter to change
|   @param chain            The chain number
|   @param seed             Pointer to the seed of the random number gen.
|   @param lnPriorRatio     Pointer to the log prior ratio (out)
|   @param lnProposalRatio  Pointer to the log proposal (Hastings) ratio (out)
|   @param mvp              Pointer to tuning parameter(s)
------------------------------------------------------------------*/
int Move_GeneTree2 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             i, numGeneTrees, numUpperTriang;
    double          newLnProb, oldLnProb, backwardLnProposalProb, forwardLnProposalProb,
                    *oldMinDepths, *modMinDepths, forwardLambda, backwardLambda, mean;
    Tree            *newSpeciesTree, *oldSpeciesTree, **geneTrees;
    ModelInfo       *m;

    // Calculate number of gene trees
    numGeneTrees = numTopologies - 1;

    // Get model settings
    m = &modelSettings[param->relParts[0]];

    // Get species tree (this trick is possible because we always copy tree params)
    newSpeciesTree = GetTree (m->speciesTree, chain, state[chain]);
    oldSpeciesTree = GetTree (m->speciesTree, chain, state[chain] ^ 1);

    // Get gene trees
    geneTrees = (Tree **) SafeCalloc (2*numGeneTrees, sizeof(Tree *));
    for (i=0; i<m->speciesTree->nSubParams; i++) {
        geneTrees[i] = GetTree(m->speciesTree->subParams[i], chain, state[chain]);
        }

    // Allocate space for depth matrix copy
    numUpperTriang = numSpecies * (numSpecies - 1) / 2;
    oldMinDepths   = (double *) SafeCalloc (2*numUpperTriang, sizeof(double));
    modMinDepths   = oldMinDepths + numUpperTriang;

    // Get min depth matrix for old gene trees
    GetMinDepthMatrix(geneTrees, numTopologies-1, depthMatrix);

    // Save a copy
    for (i=0; i<numUpperTriang; i++)
        oldMinDepths[i] = depthMatrix[i];

    // Get forward lambda
    GetMeanDist(oldSpeciesTree, depthMatrix, &mean);
    // if (mean < 1E-6) mean = 1E-6;
    forwardLambda = 1.0 / mean;

    // Calculate joint probability of old gene trees and old species tree
    oldLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, oldSpeciesTree, chain);

    // Modify the picked gene tree using code from a regular MrBayes move (no tuning parameter, so passing on mvp is OK)
    Move_NNIClock(param, chain, seed, lnPriorRatio, lnProposalRatio, mvp);

    // Update the min depth matrix
    GetMinDepthMatrix(geneTrees, numTopologies-1, depthMatrix);

    // Copy the min depth matrix
    for (i=0; i<numUpperTriang; i++)
        modMinDepths[i] = depthMatrix[i];

    // Modify the min depth matrix
    ModifyDepthMatrix (forwardLambda, modMinDepths, seed);

    // Get a new species tree
    if (GetSpeciesTreeFromMinDepths (newSpeciesTree, modMinDepths) == ERROR) {
        abortMove = YES;
        free (geneTrees);
        free (oldMinDepths);
        return (NO_ERROR);
        }
    
    // Calculate joint probability of new gene trees and new species tree
    newLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, newSpeciesTree, chain);

    // Get backward lambda
    GetMeanDist(newSpeciesTree, depthMatrix, &mean);
    // if (mean < 1E-6) mean = 1E-6;
    backwardLambda = 1.0 / mean;

    // Get proposal probability of old species tree
    backwardLnProposalProb = LnProposalProbSpeciesTree (oldSpeciesTree, oldMinDepths, backwardLambda);

    // Get proposal probability of new species tree
    forwardLnProposalProb = LnProposalProbSpeciesTree (newSpeciesTree, depthMatrix, forwardLambda);

    // Update prior ratio taking species tree into account
    (*lnPriorRatio) += (newLnProb - oldLnProb);
        
    // Update proposal ratio based on this move
    (*lnProposalRatio) += (backwardLnProposalProb - forwardLnProposalProb);

    // Free allocated memory
    free (geneTrees);
    free (oldMinDepths);

    return (NO_ERROR);
}


/**-----------------------------------------------------------------
|
|   Move_GeneTree3: Propose a new gene tree using ParsSPRClock
|
|   @param param            The parameter to change
|   @param chain            The chain number
|   @param seed             Pointer to the seed of the random number gen.
|   @param lnPriorRatio     Pointer to the log prior ratio (out)
|   @param lnProposalRatio  Pointer to the log proposal (Hastings) ratio (out)
|   @param mvp              Pointer to tuning parameter(s)
------------------------------------------------------------------*/
int Move_GeneTree3 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             i, numGeneTrees, numUpperTriang;
    double          newLnProb, oldLnProb, backwardLnProposalProb, forwardLnProposalProb,
                    *oldMinDepths, *modMinDepths, forwardLambda, backwardLambda, mean;
    Tree            *newSpeciesTree, *oldSpeciesTree, **geneTrees;
    ModelInfo       *m;

    // Calculate number of gene trees
    numGeneTrees = numTopologies - 1;

    // Get model settings
    m = &modelSettings[param->relParts[0]];

    // Get species tree (this trick is possible because we always copy tree params)
    newSpeciesTree = GetTree (m->speciesTree, chain, state[chain]);
    oldSpeciesTree = GetTree (m->speciesTree, chain, state[chain] ^ 1);

    // Get gene trees
    geneTrees = (Tree **) SafeCalloc (2*numGeneTrees, sizeof(Tree *));
    for (i=0; i<m->speciesTree->nSubParams; i++) {
        geneTrees[i] = GetTree(m->speciesTree->subParams[i], chain, state[chain]);
        }

    // Allocate space for depth matrix copy
    numUpperTriang = numSpecies * (numSpecies - 1) / 2;
    oldMinDepths   = (double *) SafeCalloc (2*numUpperTriang, sizeof(double));
    modMinDepths   = oldMinDepths + numUpperTriang;

    // Get min depth matrix for old gene trees
    GetMinDepthMatrix(geneTrees, numTopologies-1, depthMatrix);

    // Save a copy
    for (i=0; i<numUpperTriang; i++)
        oldMinDepths[i] = depthMatrix[i];

    // Get forward lambda
    GetMeanDist(oldSpeciesTree, depthMatrix, &mean);
    // if (mean < 1E-6) mean = 1E-6;
    forwardLambda = 1.0 / mean;

    // Calculate joint probability of old gene trees and old species tree
    oldLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, oldSpeciesTree, chain);

    // Modify the picked gene tree using code from a regular MrBayes move
    Move_ParsSPRClock(param, chain, seed, lnPriorRatio, lnProposalRatio, mvp);

    // Update the min depth matrix
    GetMinDepthMatrix(geneTrees, numTopologies-1, depthMatrix);

    // Copy the min depth matrix
    for (i=0; i<numUpperTriang; i++)
        modMinDepths[i] = depthMatrix[i];

    // Modify the min depth matrix
    ModifyDepthMatrix (forwardLambda, modMinDepths, seed);

    // Get a new species tree
    if (GetSpeciesTreeFromMinDepths (newSpeciesTree, modMinDepths) == ERROR) {
        abortMove = YES;
        free (geneTrees);
        free (oldMinDepths);
        return (NO_ERROR);
        }
   
    // Calculate joint probability of new gene trees and new species tree
    newLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, newSpeciesTree, chain);

    // Get backward lambda
    GetMeanDist(newSpeciesTree, depthMatrix, &mean);
    // if (mean < 1E-6) mean = 1E-6;
    backwardLambda = 1.0 / mean;

    // Get proposal probability of old species tree
    backwardLnProposalProb = LnProposalProbSpeciesTree (oldSpeciesTree, oldMinDepths, backwardLambda);

    // Get proposal probability of new species tree
    forwardLnProposalProb = LnProposalProbSpeciesTree (newSpeciesTree, depthMatrix, forwardLambda);

    // Update prior ratio taking species tree into account
    (*lnPriorRatio) += (newLnProb - oldLnProb);
        
    // Update proposal ratio based on this move
    (*lnProposalRatio) += (backwardLnProposalProb - forwardLnProposalProb);

    // Free allocated memory
    free (geneTrees);
    free (oldMinDepths);

    return (NO_ERROR);
}


/*-----------------------------------------------------------------------------------
|
|   Move_NodeSliderGeneTree: Move the position of one (root or nonroot) node in a
|      gene tree inside a species tree.
|
-------------------------------------------------------------------------------------*/
int Move_NodeSliderGeneTree (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, *nEvents;
    MrBFlt      window, minDepth, maxDepth, oldDepth, newDepth,
                oldLeftLength=0.0, oldRightLength=0.0, clockRate,
                oldPLength=0.0, lambda=0.0, nu=0.0, igrvar=0.0,
                *brlens=NULL, *tk02Rate=NULL, *igrRate=NULL, *popSizePtr;
    TreeNode    *p, *q;
    ModelInfo   *m;
    Tree        *geneTree, *speciesTree;
    Param       *subParm;

    window = mvp[0]; /* window size */
 
    m = &modelSettings[param->relParts[0]];

    /* get gene tree and species tree */
    geneTree    = GetTree (param, chain, state[chain]);
    speciesTree = GetTree (m->speciesTree, chain, state[chain]);

    /* get population size(s) */
    popSizePtr = GetParamVals(m->popSize, chain, state[chain]);

    /* get clock rate */
    if (m->clockRate == NULL)
        clockRate = 1.0;
    else
        clockRate = *GetParamVals(m->clockRate, chain, state[chain]);

    /* pick a node to be changed */
    p = geneTree->intDownPass[(int)(RandomNumber(seed)*geneTree->nIntNodes)];

#   if defined (DEBUG_CSLIDER)
    printf ("Before node slider (gene tree):\n");
    printf ("Picked branch with index %d and depth %f\n", p->index, p->nodeDepth);
    if (p->anc->anc == NULL)
        printf ("Old clock rate: %f\n", clockRate);
    ShowNodes (t->root, 0, t->isRooted);
    getchar();
#   endif

    /* get gene tree prior prob before move */
    (*lnPriorRatio) -= LnPriorProbGeneTree(geneTree, clockRate, speciesTree, popSizePtr);

    /* store values needed later for prior calculation (relaxed clocks) */
    oldPLength = p->length;
    if (p->left != NULL)
        {
        oldLeftLength = p->left->length;
        oldRightLength = p->right->length;
        }
    else
        oldLeftLength = oldRightLength = 0.0;

    /* find species tree branch to which the gene tree node belongs */
    MapGeneTreeToSpeciesTree(geneTree, speciesTree);
    q = NULL;
    for (i=0; i<speciesTree->nNodes-1; i++)
        {
        q = speciesTree->allDownPass[i];
        if (p->x == q->index)
            break;
        }
    assert (q != NULL && p->x == q->index);

    /* determine lower and upper bound */
    minDepth = p->left->nodeDepth + POS_MIN;
    if (p->right->nodeDepth + POS_MIN > minDepth)
        minDepth = p->right->nodeDepth + POS_MIN;
    if (q->nodeDepth + POS_MIN > minDepth)
        minDepth = q->nodeDepth + POS_MIN;
    if (p->anc->anc == NULL)
        maxDepth = TREEHEIGHT_MAX;
    else
        maxDepth = p->anc->nodeDepth - POS_MIN;
    
    /* abort if impossible */
    if (minDepth >= maxDepth)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    if (maxDepth-minDepth < window)
        {
        window = maxDepth-minDepth;
        }

    /* pick the new node depth */
    oldDepth = p->nodeDepth;
    newDepth = oldDepth + (RandomNumber(seed) - 0.5) * window;
    
    /* reflect the new node depth */
    while (newDepth < minDepth || newDepth > maxDepth)
        {
        if (newDepth < minDepth)
            newDepth = 2.0 * minDepth - newDepth;
        if (newDepth > maxDepth)
            newDepth = 2.0 * maxDepth - newDepth;
        }
    p->nodeDepth = newDepth;

    /* determine new branch lengths around p and set update of transition probabilities */
    if (p->left != NULL)
        {
        p->left->length = p->nodeDepth - p->left->nodeDepth;
        assert (p->left->length >= POS_MIN);
        p->left->upDateTi = YES;
        p->right->length = p->nodeDepth - p->right->nodeDepth;
        assert (p->right->length >= POS_MIN);
        p->right->upDateTi = YES;
        }
    if (p->anc->anc != NULL)
        {
        p->length = p->anc->nodeDepth - p->nodeDepth;
        assert (p->length >= POS_MIN);
        p->upDateTi = YES;
        }

    /* set flags for update of cond likes from p and down to root */
    q = p;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }

    /* calculate proposal ratio */
    (*lnProposalRatio) = 0.0;

    /* calculate prior ratio */
    (*lnPriorRatio) += LnPriorProbGeneTree (geneTree, clockRate, speciesTree, popSizePtr);

    /* adjust proposal and prior ratio for relaxed clock models */
    for (i=0; i<param->nSubParams; i++)
        {
        subParm = param->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);

            /* proposal ratio */
            if (p->left != NULL)
                {
                (*lnProposalRatio) += nEvents[p->left->index ] * log (p->left->length  / oldLeftLength);
                (*lnProposalRatio) += nEvents[p->right->index] * log (p->right->length / oldRightLength);
                }
            if (p->anc->anc != NULL)
                (*lnProposalRatio) += nEvents[p->index] * log (p->length / oldPLength);

            /* prior ratio */
            if (p->anc->anc == NULL) // two branches changed in same direction
                (*lnPriorRatio) += lambda * (2.0 * (oldDepth - newDepth));
            else if (p->left != NULL) // two branches changed in one direction, one branch in the other direction
                (*lnPriorRatio) += lambda * (oldDepth - newDepth);
            else /* if (p->left == NULL) */ // one branch changed
                (*lnPriorRatio) += lambda * (newDepth - oldDepth);

            /* update effective evolutionary lengths */
            if (UpdateCppEvolLengths (subParm, p, chain) == ERROR)
                {
                abortMove = YES;
                return (NO_ERROR);
                }
            }
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);

            /* no proposal ratio effect */

            /* prior ratio */
            if (p->left != NULL)
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[p->index], nu*oldLeftLength, tk02Rate[p->left->index]);
                (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[p->index], nu*oldRightLength, tk02Rate[p->right->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->index], nu*p->left->length, tk02Rate[p->left->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->index], nu*p->right->length, tk02Rate[p->right->index]);
                }
            if (p->anc->anc != NULL)
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[p->anc->index], nu*oldPLength, tk02Rate[p->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->anc->index], nu*p->length, tk02Rate[p->index]);
                }

            /* update effective evolutionary lengths */
            if (p->left != NULL)
                {
                brlens[p->left->index] = p->left->length * (tk02Rate[p->left->index]+tk02Rate[p->index])/2.0;
                brlens[p->right->index] = p->right->length * (tk02Rate[p->right->index]+tk02Rate[p->index])/2.0;
                }
            if (p->anc->anc != NULL)
                {
                brlens[p->index] = p->length * (tk02Rate[p->index]+tk02Rate[p->anc->index])/2.0;
                }
            }
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            
            if (p->left != NULL)
                {
                (*lnPriorRatio) -= LnProbGamma (oldLeftLength/igrvar, oldLeftLength/igrvar, igrRate[p->left->index ]);
                (*lnPriorRatio) -= LnProbGamma (oldRightLength/igrvar, oldRightLength/igrvar, igrRate[p->right->index]);
                (*lnPriorRatio) += LnProbGamma (p->left->length/igrvar, p->left->length/igrvar, igrRate[p->left->index ]);
                (*lnPriorRatio) += LnProbGamma (p->right->length/igrvar, p->right->length/igrvar, igrRate[p->right->index]);
                }
            if (p->anc->anc != NULL)
                {
                (*lnPriorRatio) -= LnProbGamma (oldPLength/igrvar, oldPLength/igrvar, igrRate[p->index]);
                (*lnPriorRatio) += LnProbGamma (p->length /igrvar, p->length /igrvar, igrRate[p->index]);
                }

            if (p->left != NULL)
                {
                brlens[p->left->index ] = igrRate[p->left->index ] * p->left->length;
                brlens[p->right->index] = igrRate[p->right->index] * p->right->length;
                }
            if (p->anc->anc != NULL)
                {
                brlens[p->index] = igrRate[p->index] * p->length;
                }
            }
        }
    
#   if defined (DEBUG_CSLIDER)
    printf ("After node slider (gene tree):\n");
    printf ("Old depth: %f -- New depth: %f -- LnPriorRatio %f -- LnProposalRatio %f\n",
        oldDepth, newDepth, (*lnPriorRatio), (*lnProposalRatio));
    ShowNodes (t->root, 0, t->isRooted);
    getchar();
#   endif

    return (NO_ERROR);
    
}


/*------------------------------------------------------------------
|
|   Move_SpeciesTree: Propose a new species tree
|
------------------------------------------------------------------*/
int Move_SpeciesTree (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             i, numGeneTrees, numUpperTriang;
    double          newLnProb, oldLnProb, backwardLnProposalProb, forwardLnProposalProb, *modMinDepths,
                    forwardLambda, backwardLambda, lambdaDiv, mean;
    Tree            *newSpeciesTree, *oldSpeciesTree, **geneTrees;
    ModelInfo       *m;

    /* get tuning parameter (lambda divider) */
    lambdaDiv = mvp[0];

    /* calculate number of gene trees */
    numGeneTrees = param->nSubParams;

    /* get model settings */
    m = &modelSettings[param->relParts[0]];

    /* get new and old species trees */
    newSpeciesTree = GetTree (m->speciesTree, chain, state[chain]);
    oldSpeciesTree = GetTree (m->speciesTree, chain, state[chain] ^ 1);

    /* get gene trees */
    geneTrees = (Tree **) SafeCalloc (numGeneTrees, sizeof(Tree*));
    for (i=0; i<param->nSubParams; i++)
        geneTrees[i] = GetTree(param->subParams[i], chain, state[chain]);

    /* get minimum depth matrix */
    GetMinDepthMatrix(geneTrees, numGeneTrees, depthMatrix);

    /* get forward lambda */
    GetMeanDist(oldSpeciesTree, depthMatrix, &mean);
    forwardLambda = 1.0 / (mean * lambdaDiv);

    /* make a copy for modification */
    numUpperTriang = numSpecies * (numSpecies - 1) / 2;
    modMinDepths = (double *) SafeCalloc (numUpperTriang, sizeof(double));
    for (i=0; i<numUpperTriang; i++)
        modMinDepths[i] = depthMatrix[i];

    /* modify minimum depth matrix */
    ModifyDepthMatrix (forwardLambda, modMinDepths, seed);

    /* construct a new species tree from the modified constraints */
    if (GetSpeciesTreeFromMinDepths(newSpeciesTree, modMinDepths) == ERROR) {
        abortMove = YES;
        free (modMinDepths);
        free (geneTrees);
        return (NO_ERROR);
        }

    /* get lambda for back move */
    GetMeanDist(newSpeciesTree, depthMatrix, &mean);
    backwardLambda = 1.0 / (mean * lambdaDiv);

    /* calculate proposal ratio */
    backwardLnProposalProb = LnProposalProbSpeciesTree (oldSpeciesTree, depthMatrix, backwardLambda);
    forwardLnProposalProb  = LnProposalProbSpeciesTree (newSpeciesTree, depthMatrix, forwardLambda);
    (*lnProposalRatio) = backwardLnProposalProb - forwardLnProposalProb;

#   if defined (BEST_MPI_ENABLED)
    // Broadcast the proposed species tree to all processors if MPI version
#   endif

#   if defined (BEST_MPI_ENABLED)
    // Let each processor calculate the ln probability ratio of its current gene tree(s)
    //    given the new and old species tree in the MPI version

    // Assemble the ln probability ratios across the processors and to lnPriorRatio
#   else
    /* calculate the ln probability ratio of the current gene trees
       given the new and old species trees */
    newLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, newSpeciesTree, chain);
    oldLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, oldSpeciesTree, chain);
#   endif

    /* set (*lnPriorRatio) to ln probability ratio */
    (*lnPriorRatio) = (newLnProb - oldLnProb);
    
    /* free allocated space */
    free (modMinDepths);
    free (geneTrees);

    return (NO_ERROR);
}


/** Show upper triangular matrix */
void ShowUpperTriangMatrix (double *values, int squareSize)
{
    int     i, j, index;

    index = 0;
    printf ("Upper triang matrix:\n");
    for (i=0; i<squareSize; i++) {
        for (j=0; j<i; j++)
            printf ("         ");
        for (j=i+1; j<squareSize; j++) {
            printf ("%.6f ", values[index]);
            index++;
            }
        printf ("\n");
        }
    printf ("\n");
}

