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

#include    <assert.h>

#include	"best.h"
#include	"command.h"
#include    "globals.h"
#include    "mb.h"
#include    "mbmath.h"
#include    "mcmc.h"
#include    "model.h"
#include    "sumt.h"
#include    "tree.h"
#include    "utils.h"

#if 0
int         Constraint(SPTree *genetree, int numgenetree, SPTree *speciestree, Distance *constraint);
double      CalNodeAge(int node, SPTree *tree);
double 	    ChangeBrlen(SPTree *speciestree, int spnode, Tree *genetree, TreeNode *p);
int         ChangeConstraint(Distance *dist, int nconstraints);
int 	    CheckConstraint(SPTree *genetrees, int ngene, Distance *constraint, int nconstraints);
long int    ClockTreeNodeDist(SPTree *clocktree, int ngene, Distance *dist);
void        CopyParam(int chn);
void        DelAllmissingspecies(int nummissing,int *missnode, SPTree *speciestree);
void        DeleteaSpecies(int inode,SPTree *speciestree);
void        FindDescendantTaxa(SPTree *tree, int inode, int *taxa, int *index);
int         FindSpnodeDownGenenode(SPTree *speciestree, int spnode, TreeNode *p);
int         FindspNodeBelow(SPTree *speciestree, int spnode, double dis_gene);
int         FindaPosition(int nodenumber,int root,double currentdistance,double *position,SPTree *speciestree);  
int         GetCoaltime(SPTree *genetree,SPTree *speciestree, CoalTime *coal);
int         GetNcoal(int inode, CoalTime *coal, SPTree *genetree, SPTree *speciestree, CoalTime *coal_t,int *index);
int         GetConstraints(double **distanceMatrix, Distance *constr);
long int    GetMinDists(SPTree *geneTrees, int numGeneTrees, double **distanceMatrix);
int         LnLikehood1Tree(SPTree *genetree, SPTree *speciestree, double *lnp);
int         LnLikehood1Tree_invgamma(SPTree *genetree, SPTree *speciestree, double *a, double *b);
int         Move_SPSpeciation (SParam *param, int chain, long *seed, double *lnLikeRatio, double *lnPriorRatio, double *lnProposalRatio, double *mvp);
int         Move_SPExtinction (SParam *param, int chain, long *seed, double *lnLikeRatio, double *lnPriorRatio, double *lnProposalRatio, double *mvp);
int         MaximumTree(Distance *dist, SPTree *speciestree, int nconstraints);
double      NodeDistance(SPTree *tree, int inode, int jnode);
long int    OneClockTreeNodeDist(SPTree *clocktree, Distance *dist);
double      PopNodeDistance(int inode, int jnode, SPTree *tree);
int 		populationMutation (Tree *genetree, SPTree *speciestree, MrBFlt genemu);
MrBFlt 	    Prob_sptree(Distance *tau,int ntime);
int         poisson(double x);
void        PrintInf(FILE *fp);
void        quick_struct(Distance *item,int count);
void        qs_struct(Distance *item,int left,int right);
int	        ReadaTree (FILE *fTree,SPTree *tree);
void        ReadGeneTree(FILE *fTree);
int         StartSptree(SPTree *speciestree, int numchange);
double 	    SPLnP1 (double t, double l, double m, double r);
double 	    SPLnVt (double t, double l, double m, double r);
int         SPLogLike_invgamma(SPTree *genetree, SPTree *speciestree, double *lnl);
int         SPPickProposal (void);
void        SPPreparePrintFiles (void);
int			SPSaveSprintf(char **target, int *targetLen, char *fmt, ...);
void        SPMrBayesPrint (char *format, ...);
int         SPLogLike(SPTree *genetree, SPTree *speciestree, double *lnl);
int         SPTreeConstraint(Distance *minimumdistance, Distance *distance, long int nconstraints, int nspecies);
void        SPWriteTreeToFile (SPTree *tree, int inode, int showBrlens, int showTheta, int showMu, int isRooted);
void        SPPrintToScreen (int curGen, time_t endingT, time_t startingT);
int         SPLogPrior(SPTree *speciestree, double *lnprior);
double 	    Toclocktree(SPTree *t, int node);
void        ToSingleGenetree(Tree *file,int i,MrBFlt GeneMu);
void        ToGenetree(Tree *file[], int *updatedtreeid, int nfile, double *GeneMu);
double 	    TreeL(Tree *t);
#endif

/****************************** Local functions converted by Fredrik from best code *****************************/
int         CompareDepths (const void *x, const void *y);
int         CompareDoubles (const void *x, const void *y);
int         CompareNodes (const void *x, const void *y);
int         ConvertToClockTrees(Tree **geneTrees, int numGeneTrees);
int         GetSpeciesTreeFromMinDepths (Tree* speciesTree, double *depthMatrix);
int         GetMinDepthMatrix (Tree **geneTrees, int numGeneTrees, double *depthMatrix);
double      LnJointGeneTreeSpeciesTreePr(Tree **geneTrees, int numGeneTrees, Tree *speciesTree, int chain, int state);
double      LnCoalescenceProb (TreeNode *geneTreeNode, TreeNode *speciesTreeNode, int *nEvents, double *coalTimes, int ploidy);
double      LnPriorProbGeneTree (Tree *geneTree, double mu, Tree *speciesTree, double *popSizePtr);


#if 0
McmcPara 	mcmc;
SPTree 	    sptree;
SPTree	    *gtree;
ModelParam	modelParam;
int		    nGene;
char		spacer[10]="  ";
FILE		*fptree;     		/* Output tree file */
FILE		*fpparm;     		/* Output parameter file */          
int		    *spnode;       /*Vector of taxaIDs for each tip */
long int 	seed;
int		    curGeneration=0;  /*should be defined in mcmc.c and used to print species tree*/
double	    speciationR, extinctionR, sampleF; 

/*global variables*/
int		jointGenePr=0 ; 		/*best or mrbayes*/
int		spMupr = 0; 		/*non-clock species tree or clock species tree*/
int		spMupralpha = 5; 		/*for non-clock species tree model, gamma(spMupralpha,r/spMupralpha) for mutation rate ratio r*/
int		speciestreePr = 0; 	/*birth-death prior or uniform prior*/
double 	poissonMean = 1.0; 	/*the number of nodes changed for proposing a new species tree*/
double  ***speciesTreeDistanceMatrix;      // distance matrixces for species trees
double  ***geneTreeDistanceMatrix[];         // distance matrices for gene trees
#endif



/**************************** Functions that interact directly with MrBayes ****************************/

// TODO: Set the global variables below in main mrbayes code

/* Global best variables */
SafeLong    **speciesPairSets;
double      *depthMatrix;

// P_GENERATE, P_PARTRATE

// Check consistency of best model in CheckModel


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
    speciesPairSets    = (SafeLong **) calloc (numUpperTriang, sizeof(SafeLong *));
    speciesPairSets[0] = (SafeLong *)  calloc (numUpperTriang*nLongsNeeded, sizeof(int));
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





/** Convert nonclock gene trees to clock gene trees */
int ConvertToClockTrees(Tree **geneTrees, int numGeneTrees)
{
    return (NO_ERROR);
}





/**-----------------------------------------------------------------
|
|	FillSpeciesTreeParams: Fill in species trees (start value)
|
------------------------------------------------------------------*/
int FillSpeciesTreeParams (SafeLong *seed, int fromChain, int toChain)

{
    int			i, k, chn, numGeneTrees, geneTreesAreClock, numUpperTriang;
	Param		*p;
	Tree		*speciesTree, **geneTrees;
    double      *depthMatrix;

    // Allocate space for the temporary depthMatrix variable used in this function.
    // We cannot use the global depthMatrix variable because it is not allocated until
    // the chain is run.
    numUpperTriang     = (numSpecies * (numSpecies-1)) / 2;
    depthMatrix = (double *) calloc (numUpperTriang, sizeof(double));

    // Use global variable numTopologies to calculate number of gene trees
    // There is one topology for the species tree, the other ones are gene trees
    // The number of current divisions is not safe because one gene tree can have
    // several partitions, for instance if we model the different genes on the
    // mitochondrion differently, or if we assign different rates to the codon
    // site positions in a sequence
    numGeneTrees = numTopologies - 1;
    geneTrees   = (Tree **) calloc (numGeneTrees, sizeof(Tree*));

    // Build species trees for state 0
	for (chn=fromChain; chn<toChain; chn++)
		{
		for (k=0; k<numParams; k++)
			{
			p = &params[k];
			if (p->paramType == P_SPECIESTREE && p->fill == YES)
	            {
                // Find species tree and gene trees
                speciesTree = GetTree(p, chn, 0);
                for (i=0; i<p->nSubParams; i++)
                    geneTrees[i] = GetTree(p->subParams[i], chn, 0);

                // If gene trees are non-clock, convert them to clock trees
                if (geneTrees[0]->isClock == NO)
                    geneTreesAreClock = NO;
                else
                    geneTreesAreClock = YES;
                if (geneTreesAreClock == NO)
                    ConvertToClockTrees(geneTrees, numGeneTrees);

                // Get minimum depth matrix for species tree
                GetMinDepthMatrix (geneTrees, numGeneTrees, depthMatrix);

                // Get a species tree from min depth matrix
                GetSpeciesTreeFromMinDepths(speciesTree, depthMatrix);

                // Label the tips
                if (LabelTree (speciesTree, speciesNameSets[speciespartitionNum].names) == ERROR)
                    {
                    free (depthMatrix);
					return (ERROR);
                    }
    
                // Flip fill flag to NO if last chain to fill
                if (chn == toChain-1)	/* last chain to fill */
					p->fill = NO;

                // Free converted trees
                if (geneTreesAreClock == NO)
                    {
                    for (i=0; i<numGeneTrees; i++)
                        FreeTree(geneTrees[i]);
                    }
                }
            }
        }

    // Free gene trees
    free (geneTrees);

    // Free depthMatrix
    free (depthMatrix);

    return (NO_ERROR);
}





/**-----------------------------------------------------------------
|
|	FreeBestChainVariables: Free best variables used during an mcmc
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

    SafeFree ((void**)(&depthMatrix));    // sets depthMatrix to NULL

    memAllocs[ALLOC_BEST] = NO;
}





/**---------------------------------------------------------------------------------------
|
|   LnCoalescenceProb: recursive function to calculate coalescence ln probability. The
|   time complexity is O(n), where n is the number of nodes in the gene tree.
|
-----------------------------------------------------------------------------------------*/
double LnCoalescenceProb (TreeNode *geneTreeNode, TreeNode *speciesTreeNode, int *nEvents, double *coalTimes, int ploidy)
{
    int     i, k, nLongsNeeded, nRemaining, nIn, nOut;
    double  lnProb, theta, timeInterval;

    lnProb = 0.0;

    if (geneTreeNode->nodeDepth <= speciesTreeNode->nodeDepth) {
        // climb up species tree (or stop)
        if (speciesTreeNode->left == NULL) {
            assert (geneTreeNode->left == NULL);
            speciesTreeNode->x++;       // increase number of lineages in == starting lineages
            (*nEvents) = 0;
            return 0.0;
        }
        else /* if (speciesTreeNode->left != NULL) */ {
            nLongsNeeded = (numSpecies - 1) / nBitsInALong + 1;
            speciesTreeNode->x++;   // increase number of lineages in
            if (IsPartNested(geneTreeNode->partition, speciesTreeNode->left->partition, nLongsNeeded))
                return LnCoalescenceProb (geneTreeNode, speciesTreeNode->left, nEvents, coalTimes, ploidy);
            else
                return LnCoalescenceProb (geneTreeNode, speciesTreeNode->right, nEvents, coalTimes, ploidy);
        }
    }
    else {
        // climb up gene tree
        lnProb += LnCoalescenceProb (geneTreeNode->left , speciesTreeNode, nEvents, coalTimes, ploidy);
        lnProb += LnCoalescenceProb (geneTreeNode->right, speciesTreeNode, nEvents, coalTimes, ploidy);

        // Add current coalescent event
        coalTimes[(*nEvents)] = geneTreeNode->nodeDepth;
        (*nEvents)++;
        
        // If the last coalescent event, we add to probability and reset nEvents
        nIn = speciesTreeNode->x;
        nRemaining = nIn  - (*nEvents);
        if (speciesTreeNode->anc->anc == NULL)
            nOut = 1;
        else
            nOut = speciesTreeNode->anc->x;

        if (nRemaining == nOut) {

            // Get theta and adjust if not diploid gene
            theta = speciesTreeNode->y;
            if (ploidy == HAPLOID)
                theta *= 0.25;      // Shouldn't this be 0.50 ?? -- Fredrik
            else if (ploidy == ZLINKED)
                theta *= 0.75;

            // Sort coalescent times
            qsort((void *)(coalTimes), (size_t)(*nEvents), sizeof(double), CompareDoubles);

            // Adjust probability
            lnProb += (*nEvents) * log (2.0 / theta);
            for (k=nIn, i=0; k>nOut; k--, i++) {
                if (i == 0)
                    timeInterval = coalTimes[0] - speciesTreeNode->nodeDepth;
                else
                    timeInterval = coalTimes[i] - coalTimes[i-1];
                lnProb -= (k * (k - 1)) / (theta * timeInterval);
            }

            if (nOut > 1) {
                timeInterval = speciesTreeNode->anc->nodeDepth - coalTimes[(*nEvents)-1];
                lnProb -= (nOut * (nOut - 1)) / (theta * timeInterval);
            }

            // Reset number of events
            (*nEvents) = 0;
        }
    }

    return lnProb;
}





/**-----------------------------------------------------------------
|
|	LnJointGeneTreeSpeciesTreePr: Converted from LnJointGenetreePr,
|   SPLogLike, SPLogPrior.
|
|   In this function we calculate the entire probability of the species
|   tree, including its probability given its priors, and the probability
|   of the gene trees given the species tree.
|
------------------------------------------------------------------*/
double LnJointGeneTreeSpeciesTreePr(Tree **geneTrees, int numGeneTrees, Tree *speciesTree, int chain, int state)
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
    popSizePtr = GetParamVals(m->popSize, chain, state);

    // Get clock rate
    if (speciesTree->isCalibrated == YES)
        clockRate = *GetParamVals(m->clockRate, chain, state);
    else
        clockRate = 1.0;

    // Calculate probability of gene trees given species tree
    lnLike = 0.0;
    for (i=0; i<numGeneTrees; i++) {
        // Get mu for the gene tree
        m = &modelSettings[geneTrees[i]->relParts[0]];
        mu = *GetParamVals(m->geneTreeRateMult, chain, state);
        mu *= clockRate;
        lnLike += LnPriorProbGeneTree(geneTrees[i], mu, speciesTree, popSizePtr);
    }

    // Calculate probability of species tree given its priors
    if (strcmp(mp->speciesTreeBrlensPr, "Birthdeath") == 0) {
        sR = *GetParamVals(m->speciationRates, chain, state);
        eR = *GetParamVals(m->extinctionRates, chain, state);
        sF = mp->sampleProb;
        lnPrior = 0.0;
        LnBirthDeathPriorPr(speciesTree, clockRate, &lnPrior, sR, eR, sF);
    }
    else
        lnPrior = 0.0;

    // The population size is taken care of elsewhere

    return lnLike + lnPrior;
}





/**-----------------------------------------------------------------
|
|	LnPriorProbGeneTree: Calculate the prior probability of a gene
|   tree.
|
------------------------------------------------------------------*/
double LnPriorProbGeneTree (Tree *geneTree, double mu, Tree *speciesTree, double *popSizePtr)
{ 
   	int         i, j, nLongsNeeded, nEvents, ploidy;
   	double      N, *coalTimes, lnProb;
    TreeNode    *p;
    ModelInfo   *m;
    ModelParams *mp;

    // Get model info
    m = &modelSettings[speciesTree->relParts[0]];

    // Get model params
    mp = &modelParams[speciesTree->relParts[0]];

    // Initialize species tree with # lineages in in x (initially 0) and diploid theta in d
    for (i=0; i<speciesTree->nNodes-1; i++) {
        p = speciesTree->allDownPass[i];
        p->x = 0;
        if (strcmp(mp->popVarPr, "equal") == 0)
            N = popSizePtr[p->index];
        else
            N = popSizePtr[0];
        p->d = 4.0 * N * mu;
    }
    
    // Initialize species partitions for both gene tree and species tree */
    AllocateTreePartitions(geneTree);
    AllocateTreePartitions(speciesTree);
    nLongsNeeded = (numSpecies - 1) / nBitsInALong + 1;
    for (i=0; i<geneTree->nNodes; i++) {
        p = geneTree->allDownPass[i];
        if (p->left == NULL)
            SetBit(speciespartitionId[p->index][speciespartitionNum]-1, p->partition);
        else {
            for (j=0; j<nLongsNeeded; j++)
                p->partition[j] = p->left->partition[j] | p->right->partition[j];
        }
    }
    for (i=0; i<speciesTree->nNodes; i++) {
        p = geneTree->allDownPass[i];
        if (p->left == NULL)
            SetBit(p->index, p->partition);
        else {
            for (j=0; j<nLongsNeeded; j++)
                p->partition[j] = p->left->partition[j] | p->right->partition[j];
        }
        
    }

    // Allocate space for coalescence times, number of int nodes in gene tree in the worst case
    coalTimes = (double *) calloc (geneTree->nIntNodes, sizeof(double));

    // Find ploidy setting
    if (strcmp(mp->ploidy, "Diploid"))
        ploidy = DIPLOID;
    else if (strcmp(mp->ploidy, "Haploid"))
        ploidy = HAPLOID;
    else
        ploidy = ZLINKED;

    // Finally call recursive function to calculate ln probability
    nEvents = 0;
    lnProb = LnCoalescenceProb (geneTree->root->left, speciesTree->root->left, &nEvents, coalTimes, ploidy);

    // Free space
    FreeTreePartitions(speciesTree);
    FreeTreePartitions(geneTree);
    free (coalTimes);

    return lnProb;
}





/**-----------------------------------------------------------------
|
|	Move_GeneTree: Propose a new gene tree
|
|   @param param            The parameter to change
|   @param chain            The chain number
|   @param seed             Pointer to the seed of the random number gen.
|   @param lnPriorRatio     Pointer to the log prior ratio (out)
|   @param lnProposalRatio  Pointer to the log proposal (Hastings) ratio (out)
|   @param mvp              Pointer to tuning parameter(s)
------------------------------------------------------------------*/
int Move_GeneTree (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)

{
    int             i, numGeneTrees;
    double          extSprClockTuningParam, newLnProb, oldLnProb;
    Tree			*newSpeciesTree, *oldSpeciesTree, **newGeneTrees, **oldGeneTrees;
    ModelInfo       *m;
    ModelParams     *mp;

    /* calculate number of gene trees */
    numGeneTrees = numTopologies - 1;

    /* set tuning param for extending SPR clock move */
    extSprClockTuningParam = 0.5;

    /* get model params */
	mp = &modelParams[param->relParts[0]];
	
	/* get model settings */
    m = &modelSettings[param->relParts[0]];

    /* get new and old species trees */
    newSpeciesTree = GetTree (m->speciesTree, chain, state[chain]);
    oldSpeciesTree = GetTree (m->speciesTree, chain, state[chain] ^ 1);

    // Get new and old gene trees
    newGeneTrees = (Tree **) calloc (2*numGeneTrees, sizeof(Tree *));
    oldGeneTrees = newGeneTrees + numGeneTrees;
    for (i=0; i<param->nSubParams; i++) {
        newGeneTrees[i] = GetTree(param->subParams[i], chain, state[chain]);
        oldGeneTrees[i] = GetTree(param->subParams[i], chain, state[chain]^1);
    }

    // Modify the picked gene tree using code from a regular MrBayes move
    Move_ExtSPRClock(param, chain, seed, lnPriorRatio, lnProposalRatio, &extSprClockTuningParam);

    // Update the min depth matrix
    GetMinDepthMatrix(newGeneTrees, numTopologies-1, depthMatrix);

    // Modify the min depth matrix

    // Get a new species tree
    GetSpeciesTreeFromMinDepths(newSpeciesTree, depthMatrix);
    
    // Calculate joint probability of new gene trees and new species tree
    newLnProb = LnJointGeneTreeSpeciesTreePr(newGeneTrees, numGeneTrees, newSpeciesTree, chain, state[chain]);

    // Calculate joint probability of old gene trees and old species tree
    oldLnProb = LnJointGeneTreeSpeciesTreePr(oldGeneTrees, numGeneTrees, oldSpeciesTree, chain, state[chain]^1);

    // Update prior ratio taking species tree into account
    (*lnPriorRatio) += (newLnProb - oldLnProb);
        
    // Update proposal ratio based on this move (it should be 0.0)
    (*lnProposalRatio) += 0.0;

    // Free allocated memory
    free (newGeneTrees);

    return (NO_ERROR);
}





/*------------------------------------------------------------------
|
|	Move_SpeciesTree: Propose a new species tree
|
------------------------------------------------------------------*/
int Move_SpeciesTree (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             i, numGeneTrees;
    double          newLnProb, oldLnProb;
    Tree			*newSpeciesTree, *oldSpeciesTree, **geneTrees;
    ModelInfo       *m;
    ModelParams     *mp;

    /* calculate number of gene trees */
    numGeneTrees = numTopologies - 1;

    /* get model params */
	mp = &modelParams[param->relParts[0]];
	
	/* get model settings */
    m = &modelSettings[param->relParts[0]];

    /* get new and old species trees */
    newSpeciesTree = GetTree (m->speciesTree, chain, state[chain]);
    oldSpeciesTree = GetTree (m->speciesTree, chain, state[chain] ^ 1);

    // get gene trees
    geneTrees = (Tree **) calloc (numGeneTrees, sizeof(Tree*));
    for (i=0; i<param->nSubParams; i++)
        geneTrees[i] = GetTree(param->subParams[i], chain, state[chain]);

    // Modify the species tree, given info on the gene trees
    // Based on the current minimum distance matrix, find a Poisson
    // distributed number of distances and modify them using
    GetMinDepthMatrix(geneTrees, numGeneTrees, depthMatrix);

    // ChangeConstraint();

    // Construct a new species tree from the new constraints
    GetSpeciesTreeFromMinDepths(newSpeciesTree, depthMatrix);
    
    // Calculate proposal ratio, should be 0.0
    (*lnProposalRatio) = 0.0;

#if defined (BEST_MPI_ENABLED)
    // Broadcast the proposed species tree to all processors if MPI version
#endif

    // Calculate the ln prior probability ratio of the new to old species trees from hyperpriors
    // Should be 0.0 for the uniform
    (*lnPriorRatio) = 0.0;

#if defined (BEST_MPI_ENABLED)
    // Let each processor calculate the ln probability ratio of its current gene tree(s)
    //    given the new and old species tree in the MPI version

    // Assemble the ln probability ratios across the processors and to lnPriorRatio
#else
    // Calculate the ln probability ratio of the current gene trees
    //    given the new and old species trees
    newLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, newSpeciesTree, chain, state[chain]);
    oldLnProb = LnJointGeneTreeSpeciesTreePr(geneTrees, numGeneTrees, oldSpeciesTree, chain, state[chain]);
#endif

    // Add ln probability ratio to (*lnPriorRatio)
    (*lnPriorRatio) += (newLnProb - oldLnProb);
    
    // Free allocated space
    free (geneTrees);

    return (NO_ERROR);
}





/******************** best functions converted by Fredrik *********************/


/* Compare function (Depth struct) for qsort */
int CompareDepths (const void *x, const void *y) {

    if (((Depth *)(x))->depth < ((Depth *)(y))->depth)
        return YES;
    else
        return NO;
}


/* Compare function (doubles) for qsort */
int CompareDoubles (const void *x, const void *y) {

    if (*((double *)(x)) < *((double *)(y)))
        return YES;
    else
        return NO;
}


/* Compare function (TreeNode struct) for qsort */
int CompareNodes (const void *x, const void *y) {

    if (((TreeNode *)(x))->nodeDepth < ((TreeNode*)(y))->nodeDepth)
        return YES;
    else
        return NO;
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
|   @returns    Returns NO_ERROR
----------------------------------------------------------------------*/
int GetSpeciesTreeFromMinDepths (Tree* speciesTree, double *depthMatrix) {

    int         i, j, numUpperTriang, nLongsNeeded, index, nextNodeIndex;
    Depth       *minDepth;
    PolyTree    *polyTree;
    PolyNode    *p, *q, *r, *u;

    nLongsNeeded    = ((numSpecies - 1) / nBitsInALong) + 1;
    numUpperTriang  = numSpecies*(numSpecies - 1) / 2;
    minDepth        = (Depth *) calloc (numUpperTriang, sizeof(Depth));

	// Convert depthMatrix to an array of Depth structs
    index = 0;
    for(i=0; i<numSpecies; i++) {
        for(j=i+1; j<numSpecies; j++) {
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
    polyTree->root = &polyTree->nodes[numSpecies];
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
    p->left = &polyTree->nodes[0];
    p->sib = NULL;
    p->anc = NULL;
    p->depth = -1.0;

    // Resolve bush using sorted depth structs
    nextNodeIndex = numSpecies+1;
    for(i=0; i<numUpperTriang; i++) {
            
        // Find tip corresponding to first taxon in pair
        p = &polyTree->nodes[FirstTaxonInPartition(minDepth[i].pairSet, nLongsNeeded)];
        
        // Descend tree until we find a node within which the pair set is nested
        do {
            p = p->anc;
        } while (!IsPartNested(minDepth[i].pairSet, p->partition, nLongsNeeded));

        if (p == polyTree->root && p->depth < 0.0) {

            // This is the first time we hit the root of the tree
            p->depth = minDepth[i].depth;

        }
        else if (p->left->sib->sib != NULL) {

            // This node is still a polytomy
            
            // Find left and right descendants of new node
            for (q=p->left; IsSectionEmpty(q->partition, minDepth[i].pairSet, nLongsNeeded); q=q->sib)
                ;
            for (r=q->sib;  IsSectionEmpty(r->partition, minDepth[i].pairSet, nLongsNeeded); r=r->sib)
                ;
            
            // Introduce the new node
            u = &polyTree->nodes[nextNodeIndex++];
            u->left = q;
            u->anc = p;
            if (q->sib == r)
                u->sib = r->sib;
            else
                u->sib = q->sib;

            // Create new taxon set with bitfield operations
            for (j=0; j<nLongsNeeded; j++)
                u->partition[j] = q->partition[j] | r->partition[j];

            // Patch the tree together with the new node added
            q->sib  = r;
            r->sib = NULL;
            q->anc = u;
            r->anc = u;
            if (p->left == q)
                p->left = u;
            else {
                // we need to find the former upstream sibling to q
                for (r=p->left; r!=q; r=r->sib)
                    ;
                r->sib = u;
            }
        }
        // other cases should not be added to tree
    }

    // Set traversal sequences
    GetPolyDownPass(polyTree);

    // Set branch lengths from node depths (not done automatically for us)
    for (i=0; i<polyTree->nNodes; i++) {
        p = polyTree->allDownPass[i];
        if (p->anc == NULL)
            p->length = 0.0;
        else
            p->length = p->anc->depth - p->depth;
    }

    // Copy to species tree from polytomous tree
    CopyToTreeFromPolyTree (speciesTree, polyTree);

    // Free locally allocated variables
    FreePolyTree(polyTree);
    free (minDepth);

    return(NO_ERROR);
}





/**---------------------------------------------------------------------
|
|   GetMinDepthMatrix: converted from GetMinDists (the copy of GetMinDists
|   below is mangled, sorry about that).
|
|   This algorithm scans the gene trees and calculates the minimum depth
|   (height) separating species across gene trees. The complexity of the
|   original algorithm was O(mn^3), where m is the number of gene trees and
|   n is the number of taxa in each gene tree. I think this algorithm has
|   complexity that is better on average, but the difference is small.
|
|   I have rewritten the algorithm also to show alternative techniques that
|   could be used in this and other best algorithms.
|
|   @param      geneTrees       The gene trees (in)
|   @param      depthMatrix     The minimum depth matrix, upper triangular array (out)
|   @returns    Returns ERROR or NO_ERROR
----------------------------------------------------------------------*/
int GetMinDepthMatrix (Tree **geneTrees, int numGeneTrees, double *depthMatrix) {

    int         trace=0;
	int         i, j, w, nLongsNeeded, numUpperTriang, index;
    double      maxDepth;
    TreeNode    *p;
    SafeLong    **speciesSets;

    // Allocate space for species partitions
    nLongsNeeded   = ((numSpecies -1) / nBitsInALong) + 1;   // number of longs needed in a bitfield representing a species set
    speciesSets    = (SafeLong **) calloc (2*numLocalTaxa-1, sizeof(SafeLong *));
    speciesSets[0] = (SafeLong *)  calloc ((2*numLocalTaxa-1)*nLongsNeeded, sizeof(int));
    for (i=1; i<2*numLocalTaxa-1; i++)
        speciesSets[i] = speciesSets[0] + i*nLongsNeeded;

    // Set tip species partitions once and for all
    for (i=0; i<numLocalTaxa; i++)
        SetBit(i, speciesSets[i]);

    // Set initial max depth for upper triangular matrix
    numUpperTriang = (numSpecies * (numSpecies - 1)) / 2;
    maxDepth       = geneTrees[0]->root->left->nodeDepth;
    for (i=0; i<numUpperTriang; i++)
        depthMatrix[i] = maxDepth;

    // Now we are ready to cycle over gene trees
	for (w=0; w<numGeneTrees; w++) {
		if (trace) {
            printf("\nGene %d\n",w);
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
        qsort((void *)(geneTrees[w]->intDownPass), (size_t) geneTrees[w]->nIntNodes, sizeof(TreeNode *), CompareNodes);

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

    if(trace) {
        index = 0;
        for(i=0;i<numSpecies;i++) {
            printf("\n");
	        for(j=i+1;j<numSpecies;j++) {
	            printf("%05f ",depthMatrix[index]);
	            if(depthMatrix[index]<.000001) system("PAUSE");
            }
        }
    }

    free (speciesSets[0]);
    free (speciesSets);

    return (NO_ERROR);
}




/**************************original  best functions **************************/
#if 0
//recursive function that loads the IDs of the (index) tips descending from inode into taxa
void FindDescendantTaxa(SPTree *tree, int inode, int *taxa, int *index) {
	if(inode < tree->nTaxa) {
		taxa[(*index)++] = inode;
	} else {
		FindDescendantTaxa(tree, tree->nodes[inode].sons[0], taxa, index);
		FindDescendantTaxa(tree, tree->nodes[inode].sons[1], taxa, index);
	}
}

int poisson(double x)
{
	int    poi_value;             
  	double t_sum;                 

  	poi_value = 0;
  	t_sum = 0.0;

  	while(1)
  	{
    		t_sum = t_sum - x * log(RandomNumber(&globalSeed));
    		if (t_sum >= 1.0) break;
    			poi_value++;
  	}
 	return(poi_value);
}

void ToSingleGenetree(Tree *file,int i,MrBFlt GeneMu)
{
  	TreeNode *p;
  	int j,nspecies;
  	MrBFlt tlg,tl;
  
  	nspecies = gtree[i].nTaxa;
      
   if(file->isRooted == 1)
   {
	gtree[i].root = file->root->index;
	gtree[i].nodes[gtree[i].root].sons[0] = file->root->left->index;
	gtree[i].nodes[gtree[i].root].sons[1] = file->root->right->index;
   }
   else
   {
   	gtree[i].root = nspecies*2-2;
   	gtree[i].nodes[gtree[i].root].sons[0] = file->root->left->index;
   	gtree[i].nodes[gtree[i].root].sons[1] = file->root->index;
   }
   gtree[i].nodes[gtree[i].root].nson = 2;
   gtree[i].nodes[gtree[i].root].father = -1;
   gtree[i].nodes[gtree[i].root].brlens = 0.0;
   for(j=0;j<file->nNodes;j++)
     { 
       p=file->allDownPass[j];
       if(p != file->root)
          {
              
              if(p->anc != NULL)
                     {
                      
                      if(p != file->root->left) 
                           {
                              gtree[i].nodes[p->index].brlens = p->length;
                              gtree[i].nodes[p->index].father = p->anc->index;
                           }
                      else {
                            gtree[i].nodes[p->index].brlens = p->length/2;
                            gtree[i].nodes[p->index].father = gtree[i].root;
                           }
                     }
              else gtree[i].nodes[p->index].father = -1;
              if(p->left != NULL)
                    { 
                     gtree[i].nodes[p->index].nson = 2;
                     gtree[i].nodes[p->index].sons[0] = p->left->index;
                     gtree[i].nodes[p->index].sons[1] = p->right->index;
                    }
              else
                    { 
                       gtree[i].nodes[p->index].nson = 0;
                       gtree[i].nodes[p->index].sons[0] = -2;
                       gtree[i].nodes[p->index].sons[1] = -2;
                    }
          }
       else
          {
                       gtree[i].nodes[p->index].nson = 0;
                       gtree[i].nodes[p->index].sons[0] = -2;
                       gtree[i].nodes[p->index].sons[1] = -2;
                       gtree[i].nodes[p->index].brlens = p->left->length/2;
                       gtree[i].nodes[p->index].father = gtree[i].root;
          }
             
     }
  

 
  Toclocktree(&gtree[i],gtree[i].root);
 
    tl = TreeL(file);
    tlg = 0.0;
    for(j=0;j<nspecies*2-1;j++) tlg += gtree[i].nodes[j].brlens;
 
  for(j=0;j<nspecies*2-1;j++)
	{
		gtree[i].nodes[j].brlens /= (GeneMu*tlg/tl);
            gtree[i].nodes[j].age /= (GeneMu*tlg/tl);
	} 
  
}

void        ToGenetree(Tree *file[],int *updatedtreeid, int nfile,double *GeneMu)
{
  	TreeNode *p;
  	int i,j;
  	double tlg,tl;

  	for(i=0;i<nfile;i++)
  	{
		if(file[i]->isRooted != 1)
		{
			gtree[i].root = gtree[i].nTaxa*2-2;
   			gtree[i].nodes[gtree[i].root].sons[0] = file[i]->root->left->index;
   			gtree[i].nodes[gtree[i].root].sons[1] = file[i]->root->index;
   			gtree[i].nodes[gtree[i].root].nson = 2;
   			gtree[i].nodes[gtree[i].root].father = -1;
   			gtree[i].nodes[gtree[i].root].brlens = 0.0;

   			for(j=0;j<file[i]->nNodes;j++)
     			{ 
       				p=file[i]->allDownPass[j];
       				if(p != file[i]->root)
          			{
              				if(p->anc != NULL)
                     			{
                       				if(p != file[i]->root->left) 
                           			{
                              			gtree[i].nodes[p->index].brlens = p->length;
                              			gtree[i].nodes[p->index].father = p->anc->index;
                           			}
                      				else 
						{
                            			gtree[i].nodes[p->index].brlens = p->length/2;
                            			gtree[i].nodes[p->index].father = gtree[i].root;
                           			}
                     			}	
              				else 
						gtree[i].nodes[p->index].father = -1;
              				if(p->left != NULL)
                    			{ 
                      			gtree[i].nodes[p->index].nson = 2;
                     			gtree[i].nodes[p->index].sons[0] = p->left->index;
                     			gtree[i].nodes[p->index].sons[1] = p->right->index;
                    			}
               				else
                    			{ 
                       			gtree[i].nodes[p->index].nson = 0;
                       			gtree[i].nodes[p->index].sons[0] = -2;
                       			gtree[i].nodes[p->index].sons[1] = -2;
                    			}
          			}
       				else
          			{
                       		gtree[i].nodes[p->index].nson = 0;
                       		gtree[i].nodes[p->index].sons[0] = -2;
                       		gtree[i].nodes[p->index].sons[1] = -2;
                       		gtree[i].nodes[p->index].brlens = p->left->length/2;
                       		gtree[i].nodes[p->index].father = gtree[i].root;
          			}		
         		}
			/*for(j=0;j<file[i]->nNodes;j++)
                                printf("node %d %d %d %d %lf\n",j,gtree[i].nodes[j].sons[0],gtree[i].nodes[j].sons[1],gtree[i].nodes[j].father,gtree[i].nodes[j].brlens);*/
		}
		else
		{
			for(j=0;j<file[i]->nNodes-1;j++)
                	{
 				p=file[i]->allDownPass[j];
                        	gtree[i].nodes[p->index].brlens = p->length;
				gtree[i].nodes[p->index].father = p->anc->index;

                        	if(p->left != NULL)
				{
                                gtree[i].nodes[p->index].sons[0] = p->left->index;
				gtree[i].nodes[p->index].sons[1] = p->right->index;
				gtree[i].nodes[p->index].nson = 2;
				}
				else
                                {
                                gtree[i].nodes[p->index].nson = 0;
                                gtree[i].nodes[p->index].sons[0] = -2;
                                gtree[i].nodes[p->index].sons[1] = -2;
                                }
			}
			p=file[i]->allDownPass[file[i]->nNodes-1];
			gtree[i].root = p->left->index;
			gtree[i].nodes[gtree[i].root].father = -1;

			/*for(j=0;j<file[i]->nNodes-1;j++)
				printf("node %d %d %d %d %lf\n",j,gtree[i].nodes[j].sons[0],gtree[i].nodes[j].sons[1],gtree[i].nodes[j].father,gtree[i].nodes[j].brlens);*/
		}  	
	}
  	
	for(i=0;i<nfile;i++) Toclocktree(&gtree[i],gtree[i].root);
  
	/*normalized by total branch length*/
 	if (!strcmp(modelParams[0].brlensPr,"Clock"))
      	{
		for(i=0;i<nfile;i++)
                {
                                for(j=0;j<gtree[i].nTaxa*2-1;j++)
                                {
                                gtree[i].nodes[j].brlens /= GeneMu[i];
                                gtree[i].nodes[j].age    /= GeneMu[i];
                                }
		}
      	}
	else
	{
  		for(i=0;i<nfile;i++)
   		{
    			tl = TreeL(file[i]);
    			tlg = 0.0;
    			for(j=0;j<gtree[i].nTaxa*2-1;j++) tlg += gtree[i].nodes[j].brlens;
    			for(j=0;j<gtree[i].nTaxa*2-1;j++)
			{
			gtree[i].nodes[j].brlens /= (GeneMu[i]*tlg/tl);
			gtree[i].nodes[j].age   /= (GeneMu[i]*tlg/tl);
			}
   		}
	}

}

long int ClockTreeNodeDist(SPTree *clocktree, int ngene, Distance *dist) {
	int i, j, k, w, indexa, indexb, son0, son1;
	long int indexc=0;
	int *taxa0, *taxa1;

    taxa0 = (int *) calloc (clocktree[0].nTaxa, sizeof(int));
    taxa1 = (int *) calloc (clocktree[0].nTaxa, sizeof(int));

	for(w=0; w<ngene; w++) {
		for(j=clocktree[w].nTaxa; j<2*clocktree[w].nTaxa-1;j++) {
			son0 = clocktree[w].nodes[j].sons[0];
			son1 = clocktree[w].nodes[j].sons[1];
			indexa = 0; 
			FindDescendantTaxa(&clocktree[w], son0, taxa0, &indexa);
			indexb = 0;
			FindDescendantTaxa(&clocktree[w], son1, taxa1, &indexb);
		
			for(i=0; i<indexa; i++) for(k=0; k<indexb; k++)
				if(spnode[taxa0[i]] != spnode[taxa1[k]]) { //if from different species
					dist[indexc].nodes[0] = spnode[taxa0[i]];
					dist[indexc].nodes[1] = spnode[taxa1[k]];
					dist[indexc].dist = 2*(clocktree[w].nodes[j].age);
					indexc++;
				}
		}
	}

    free (taxa0);
    free (taxa1);

    return (indexc);
}

long int OneClockTreeNodeDist(SPTree *clocktree, Distance *dist)
{
	int i, j, k, indexa, indexb, son0, son1;
	long int indexc=0;
	int taxa0[NSPECIES], taxa1[NSPECIES]; 

	for(j=clocktree->nTaxa; j<2*clocktree->nTaxa-1;j++)
	{	
		son0 = clocktree->nodes[j].sons[0];
		son1 = clocktree->nodes[j].sons[1];
		indexa = 0; 
		FindDescendantTaxa(clocktree, son0, taxa0, &indexa);
		indexb = 0;
		FindDescendantTaxa(clocktree, son1, taxa1, &indexb);
	
		for(i=0; i<indexa; i++)
			for(k=0; k<indexb; k++)
			{
				if(spnode[taxa0[i]] != spnode[taxa1[k]])
				{
					dist[indexc].nodes[0] = spnode[taxa0[i]];
					dist[indexc].nodes[1] = spnode[taxa1[k]];
					dist[indexc].dist = 2*(clocktree->nodes[j].age);
					indexc++; 
				}
			}
    }
	return (indexc);
}

#if 0
void OutNodes(treenode n, int id, int lvl) {
	int i;
    for(i=0; i<lvl; i++)
        printf("  ");
	printf("NODE %d: Age=%f, %d sons, species %d\n", id, n.age, n.nson, spnode[id]);
}

void PrintSPTree(SPTree *s,int nid,int lvl) {
	OutNodes(s->nodes[nid],nid,lvl);
	if(s->nodes[nid].nson>0) PrintSPTree(s,s->nodes[nid].sons[0],lvl+1);
	if(s->nodes[nid].nson>1) PrintSPTree(s,s->nodes[nid].sons[1],lvl+1);
}
#endif

void PrintMinMat(double **md, int nsp) {
	int i, j;

    for(i=0; i<nsp; i++)
    {
        printf("\n");
        for(j=0; j<nsp; j++)
            printf("%05f ",md[0][i*nsp+j]);
	}
	printf("\n");
}

/*Scans across the ngene clocktrees to find mindist between all species*/
long int GetMinDists(SPTree *geneTrees, int numGeneTrees, double **distanceMatrix) {

    int         trace=0;
	int         i, j, w;
    TreeNode    *father, *son0, *son1;
	long        indexc=0;

    // TODO
    int         numSpecies=0;
    int         *speciesPartition = NULL;

    for (i=0; i<numSpecies; i++)
        for (j=0;j<numSpecies; j++)
            distanceMatrix[0][i*numSpecies+j]=0.0;

	for (w=0; w<numGeneTrees; w++) {
		if (trace) {
            printf("\nGene %d\n",w);
            //ShowTree(&geneTrees[w]);
        }

        // AllocateTreePartitions(&geneTrees[w]);
		for (j=0; j<0 /*geneTrees[w].nIntNodes */; j++) {

            //for each internode
            //father = geneTrees[w].intDownPass[j];
            father = NULL;
			son0   = father->left;
			son1   = father->right;

            if(trace) {
                ShowParts(stdout, son0->partition, numLocalTaxa);
                ShowParts(stdout, son1->partition, numLocalTaxa);
            }
			
			/*
            for(i=0; i<nR; i++) for(k=0; k<nL; k++)
				if(spnode[taxa0[i]] != spnode[taxa1[k]]) { //if from different species
				   if(spnode[taxa0[i]]<spnode[taxa1[k]]) indexc=spnode[taxa0[i]]*nsp+spnode[taxa1[k]];
				   else indexc=spnode[taxa1[k]]*nsp+spnode[taxa0[i]]; //figure out where to put this min
				   if((w==0 && distanceMatrix[indexc]==0.0) || 2*(clocktree[w].nodes[j].age)<distanceMatrix[0][indexc]) {
						distanceMatrix[indexc] = 2*(clocktree[w].nodes[j].age);
						if(trace) PrintMinMat(distanceMatrix,nsp);
						if(distanceMatrix[indexc]<.000001) printf("\n0 BL ERROR: Gene %d %ld.%ld",w,indexc/nsp,indexc%nsp);
					}
				}
            */
		}
	}

    /*
	if(trace) for(i=0;i<nsp;i++) { printf("\n");
		for(j=i+1;j<nsp;j++) {
			printf("%05f ",md[0][i*nsp+j]);
			if(md[0][i*nsp+j]<.000001) system("PAUSE");
		}
	}

    				    if ((w==0 && distanceMatrix[0][indexc] == 0.0) || 2.0*(father->nodeDepth) < distanceMatrix[0][indexc])
    						distanceMatrix[0][indexc] = 2*(father->nodeDepth);
						if (trace) PrintMinMat(distanceMatrix, numSpecies);
						assert (distanceMatrix[0][indexc]>.000001);
                    }
                }
            }
            // FreeTreePartitions(&geneTrees[w]);
        }
    }

    if(trace) {
        for(i=0;i<numSpecies;i++) {
            printf("\n");
	        for(j=i+1;j<numSpecies;j++) {
	            printf("%05f ",distanceMatrix[0][i*numSpecies+j]);
	            if(distanceMatrix[0][i*numSpecies+j]<.000001) system("PAUSE");
            }
        }
    }
    */

    return (NO_ERROR);
}




double TreeL(Tree *t)
{
	int			i;
	double		tl;
	TreeNode		*p;

	tl = 0.0;
	for (i=0; i<t->nNodes; i++)
	{
		p = t->allDownPass[i];
		if (p->anc != NULL)
		{
			if (p->anc->anc == NULL)
			{
				if (t->isRooted == NO)
					tl += p->length;
			}
			else
			{
				tl += p->length;
			}
		}
	}				
	return (tl);	
}

double Toclocktree(SPTree *t, int node)
{
  	double a,b;
  
  
  	if(t->nodes[node].nson==0)
     	{
         t->nodes[node].age = 0.0;
         return(0.0);
     	}
  	else
     	{   
         a = Toclocktree(t,t->nodes[node].sons[0]) + t->nodes[t->nodes[node].sons[0]].brlens;
         b = Toclocktree(t,t->nodes[node].sons[1]) + t->nodes[t->nodes[node].sons[1]].brlens;
         
         if(a>b)
            { 
               t->nodes[t->nodes[node].sons[1]].brlens += (a-b);
               t->nodes[node].age = a;
            }
         else if(a<b)
             { 
               t->nodes[t->nodes[node].sons[0]].brlens += (b-a);
               t->nodes[node].age = b;
             }
         else
             {
               t->nodes[node].age = a;
             }
         return(t->nodes[node].age);
     	}
}

int populationMutation (Tree *genetree,SPTree *speciestree, MrBFlt genemu)
{  
	int inode, i, k, inode_gene, stop, index=0;
	int *genetreenodes;
	TreeNode *p=NULL;

	genetreenodes = (int*)malloc(2*sptree.nTaxa*sizeof(int));
	if(!genetreenodes)
		{
		printf("allocating problem for genetreenodes\n");
		return(ERROR);
		}
	FOR(i,sptree.nSpecies)
      { 
        	FOR(inode,sptree.speciesIndex[i][0])
        	{        
        		inode_gene = sptree.speciesIndex[i][inode+1];
			FOR(k,(genetree->nNodes-1))
				{
				if(genetree->allDownPass[k]->index == inode_gene)
					{
					p = genetree->allDownPass[k];
					break;
					}
				}
   
        		stop=0;
        		do{
				/*check if the node is already taken care of*/
				FOR(k,index)
              			if(p->index == genetreenodes[k]) {stop=1;break;}
				if(stop == 1)
					break;
				/*change the branch length of node p*/				
				p->length = ChangeBrlen(speciestree, i, genetree, p);

				/*copy p to genetreenode*/
				genetreenodes[index] = p->index;
                 		index++;

				/*reset p*/
				p = p->anc;
			}while(p->index != genetree->root->left->index);
		}
	} 

	/*genes may have different mutation rates affecting all branches*/
	FOR(k,(genetree->nNodes-1))
	{
		p = genetree->allDownPass[k];
		p->length *= genemu;
	}
   	
    	if(index != 2*sptree.nTaxa-2)
    	{
      	printf("Error for the coalescence time calculation\n");
       	return(ERROR);
    	}
      free(genetreenodes);
	return(NO_ERROR);
}
				
double ChangeBrlen(SPTree *speciestree, int spnode, Tree *genetree, TreeNode *p)
{
	double length;
	int 	inode_sp, jnode_sp, father;

	/*find the node in the species tree that is immidiately down the gene node p*/
    	inode_sp = FindSpnodeDownGenenode(speciestree,spnode,p);
	jnode_sp = FindSpnodeDownGenenode(speciestree,spnode,p->anc);

	if(inode_sp == jnode_sp)
	{
		length = (p->length) * (speciestree->nodes[inode_sp].mu);
	}
	else
	{
		father = speciestree->nodes[inode_sp].father;
		length = (speciestree->nodes[father].age - p->nodeDepth)*(speciestree->nodes[inode_sp].mu);
		while(father != jnode_sp)
		{
			inode_sp = father;
			father = speciestree->nodes[father].father;
			length += (speciestree->nodes[father].age - speciestree->nodes[inode_sp].age)*(speciestree->nodes[inode_sp].mu);
		}
		length += (p->anc->nodeDepth - speciestree->nodes[father].age)*(speciestree->nodes[father].mu);
  	}
	return(length);		
}

int FindSpnodeDownGenenode(SPTree *speciestree, int spnode, TreeNode *p)
{
	int father, findnode=spnode;
	double depth;

	depth = p->nodeDepth;
	father = speciestree->nodes[spnode].father;

	if(p->index < sptree.nTaxa)
		return findnode;

	while(speciestree->nodes[father].age <= depth)
	{
		if(father == speciestree->root)
		{
			findnode = father;
			break;
		}
		else
		{
			findnode = father;
			father = speciestree->nodes[father].father;
		}
	}
	
	return (findnode);

}
			
int LnJointGenetreePr(Tree *t[], int *updatedtreeid, int num_tree, double *lncoalPrior, double *GeneMu, SPTree *speciestree)
{ 	
    	double	lnLike, lnPrior; 
 	int 		k, numchange;

	/*convert non-clock gene trees to clock gene trees*/
    	ToGenetree(t,updatedtreeid, num_tree,GeneMu);

	/*generate a random species tree*/
	numchange = poisson(1/poissonMean);
	if(StartSptree(&sptree, numchange)==ERROR)
	   {
		   printf(" Problem with startSptree \n");
		   return ERROR;

	   }  

   /*Calculate likelihood*/
	if(modelParam.thetainvgamma == 1) {
		if(SPLogLike_invgamma(gtree,&sptree,&lnLike) == ERROR) {
     		printf("Error for Lnlike\n");
     		return ERROR;
      }
	} else {
		if(SPLogLike(gtree,&sptree,&lnLike) == ERROR)
          	{
           		printf("Error for Lnlike\n");
           		return ERROR;
      	}
	}

	/*Calculate prior probability*/
    	if(SPLogPrior(&sptree,&lnPrior) == ERROR)
    	{
     		printf("Error for LnPrior\n");
           	return ERROR;
    	}                                
 
	/*Copy the species tree*/
	FOR(k, (2*sptree.nSpecies-1))
   	{
         	speciestree->nodes[k].nson = sptree.nodes[k].nson;
          	speciestree->nodes[k].sons[0] = sptree.nodes[k].sons[0];
          	speciestree->nodes[k].sons[1] = sptree.nodes[k].sons[1];
		speciestree->nodes[k].father = sptree.nodes[k].father;
          	speciestree->nodes[k].brlens = sptree.nodes[k].brlens;
          	speciestree->nodes[k].age = sptree.nodes[k].age;
		speciestree->nodes[k].theta = sptree.nodes[k].theta;
		speciestree->root = sptree.root;    
		speciestree->nSpecies = sptree.nSpecies;
		/*if(k == speciestree->root) speciestree->nodes[k].theta = 0.0;*/
   	}

	/*lnproposal = (-poissonMean + numchange * log(poissonMean));
	for(k = 1; k <= numchange; k++)
		lnproposal -= log(k); */

	*lncoalPrior = (lnLike + lnPrior);
    	return NO_ERROR;		       
}

int SPPrintTreeTitle (int curGen, FILE *fout)
{

	int				i;
	char			*name, spfilename[100];
	FILE			*sumt;

	/* print the translate block information and the top of the file */
	if (curGen == 1)
		{
		/* print #NEXUS and translation block information */
		fprintf (fout, "#NEXUS\n");
		fprintf (fout, "[ID: %ld]\n", seed);
		fprintf (fout, "begin trees;\n");
		fprintf (fout, "   translate\n");
		for (i=0; i<sptree.nSpecies-1; i++)
			{
			name = taxaSetNames[i];
			fprintf (fout, "      %2d %s,\n", i+1, name);
			}
		name = taxaSetNames[i];
		fprintf (fout, "      %2d %s;\n", i+1, name);
		}

	/* write the tree in Newick format */
	fprintf (fout, "   tree rep.%d = ", curGen);
	
	/* create a file to summarize species trees in the .sptree files */
	if (curGen == 1)
		{
		sprintf(spfilename, "%s.sumt", chainParams.chainFileName);
		sumt = fopen(spfilename, "w");
		sprintf(spfilename, "%s.sptree", chainParams.chainFileName);
	
		/* print #NEXUS and translation block information */
		fprintf (sumt, "#NEXUS\n");
		fprintf (sumt, "Begin data;\n");
		fprintf (sumt, "Dimensions ntax = %d nchar = 1;\n", sptree.nSpecies);
		fprintf (sumt, "FORMAT DATATYPE = dna gap=-  MISSING=?  interleave=yes;\nMatrix\n");
		for (i=0; i<sptree.nSpecies; i++)
			{
			name = taxaSetNames[i];
			fprintf (sumt, " %s A\n", name);
			}
		fprintf (sumt, "\n;\nEND;\n");
		fprintf (sumt, "Begin mrbayes;\n  prset best = 1;\n");
		fprintf (sumt, "  sumt burnin = %d nruns = %d filename = %s contype = allcompat;\nEND;\n", chainParams.numGen/chainParams.sampleFreq/2, chainParams.numRuns, spfilename);
		
		fclose (sumt);
		}
 		   		
   	return (NO_ERROR);		
}


#if 0
int SPAddToPrintString (char *tempStr)
{
	size_t			len1, len2;
	
	len1 = (int) strlen(printString);
	len2 = (int) strlen(tempStr);
	if (len1 + len2 + 5 > printStringSize)
		{
		printStringSize += len1 + len2 - printStringSize + 200;
		printString = (char*)realloc((void *)printString, printStringSize * sizeof(char));
		if (!printString)
			{
			SPMrBayesPrint ("%s   Problem reallocating printString (%d)\n", spacer, printStringSize * sizeof(char));
			goto errorExit;
			}
		}
	strcat(printString, tempStr);	
#	if 0
	printf ("printString(%d) -> \"%s\"\n", printStringSize, printString);
#	endif	
	return (NO_ERROR);
	
	errorExit:
		return (ERROR);
}
#endif

void SPPreparePrintFiles (void)
{
	char		localFileName[100], fileName[100];

	/* Get root of local file name */
	strcpy (localFileName, mcmc.chainFileName);

	/* Prepare the .p, .t */
	sprintf (fileName, "%s.p", localFileName);
	fpparm =(FILE*)OpenTextFileW(fileName);
       sprintf (fileName, "%s.t", localFileName);
	fptree =(FILE*)OpenTextFileW(fileName);
}

void InitiateParam(void)
{ 
 	modelParam.sRprior[0]=0;
  	modelParam.sRprior[1]=10;
  	modelParam.eRprior[0]=0;
  	modelParam.eRprior[1]=10;
  	modelParam.sF=1.0;
  	modelParam.treeHeightExp = 1.0;
}

int SPLogPrior(SPTree *speciestree, double *lnprior)
{
   	/* prior for the species tree*/
   	if(speciestreePr == 0)
		*lnprior = 0.0;
   	else
 	{
  		if (SPLnBirthDeathPriorPr (lnprior, speciationR, extinctionR, sampleF,speciestree) == ERROR)
   		{	
      		printf ("   Problem calculating prior for birth-death process\n");
      		return(ERROR);
   		}
		*lnprior += log(modelParam.treeHeightExp) - modelParam.treeHeightExp*speciestree->nodes[speciestree->root].age;
	}
   	return(NO_ERROR);
} 

void SPMrBayesPrint (char *format, ...)
{

	va_list                 ptr;

	va_start (ptr, format);

	vprintf (format, ptr);
	fflush (stdout);

	va_end(ptr);
}


int ReadaTree (FILE *fTree,SPTree *tree)
{
/* Read a tree from fTree, using the parenthesis node representation of trees.
   Branch lengths are read in nodes[].branch, and branch (node) labels 
   (integers) are preceeded by # and read in nodes[].label.  If the clade label
   $ is used, the label is read into nodes[].divtime first and then moved into
   nodes[].label in the routine DownTreeCladeLabel().

   Both names and numbers for species are accepted.  
   Species names are considered case-sensitive, with trailing blanks ignored.

   copyname = 0: species numbers and names are both accepted, but names have to match 
                 the names in taxaName[], which are from the sequence data file.
              1: species names are copied into taxaName[], but species numbers
                 are accepted.  
              2: the tree can have species names only, which are copied into 
                 taxaName[].
*/
   int cnode, cfather=-1;  /* current node and father */
   int inodeb=0;  /* node number that will have the next branch length */
   int i, level=0, ch=' ';
   char  skips[]="\"\'";
   int nnode;   

   nnode=tree->nTaxa;  
   FOR(i,2*tree->nTaxa-1) {
      tree->nodes[i].father=-1;
      tree->nodes[i].nson=0; 
   }


   while(isspace(ch)) 
     {
      ch=fgetc(fTree);  
     }
   ungetc(ch,fTree);

   for (;;) {
      ch = fgetc (fTree);
      if (ch==EOF) return(-1);
      else if (!isgraph(ch) || ch==skips[0] || ch==skips[1]) continue;
      else if (ch=='(') {
         level++;
         cnode=nnode++;
   
         if(nnode>2*tree->nTaxa-1)
              {
                  printf("check tree: perhaps too many '('s");
                  exit(-1);
              }
         if (cfather>=0) {
            tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = cnode;
            tree->nodes[cnode].father=cfather;
         }
         else
            tree->root=cnode;
         cfather=cnode;
      }
      else if (ch==')') { level--;  inodeb=cfather; cfather=tree->nodes[cfather].father; }
	  else if (ch==':') if( fscanf(fTree,"%lf",&tree->nodes[inodeb].brlens) < 1 ) 
							{ printf("Error in treefile: unexpected charactor after : in place where floating point number representing branch length should be"); exit(-1); } 
      else if (ch==',') ;
      else if (ch==';' && level!=0) 
         {
            printf("; in treefile");
            exit(-1);
         }
      else if (isdigit(ch))
      { 
         ungetc(ch, fTree); 
         if( fscanf(fTree,"%d",&inodeb) < 1 )  //practicly should not happend but put it here to avoid compiler warning
		 { printf("Error in treefile."); exit(-1); }
         inodeb--;
         tree->nodes[inodeb].father=cfather;
         tree->nodes[cfather].sons[tree->nodes[cfather].nson++]=inodeb;
      }
      if (level<=0) break;
   }
   
   for ( ; ; ) {
      while(isspace(ch=fgetc(fTree)) && ch!=';' ) 
         ;
      if (ch==':')       if( fscanf(fTree, "%lf", &tree->nodes[tree->root].brlens) < 1 ) 
							{ printf("Error in treefile: unexpected charactor after : in place where floating point number representing branch length should be"); exit(-1); }
      else if (ch==';')  break;
      else  { ungetc(ch,fTree);  break; }
   }

   if(nnode!=2*tree->nTaxa-1) { printf(" # of nodes != %d\n",2*tree->nTaxa-1); exit(-1);}

   return (0);
}


int ReadControlfile(FILE *fdata)
{ 

  	int i,j,k,m,*species,max;

	 
  	/* parameters for MCMC */
    /*
    seed=swapSeed;
   	SetSeed(seed);
    */

   	/* parameters for prior */
    /*
	if(!strcmp(modelParams[0].popSizePr,"Gamma"))
  		modelParam.thetainvgamma = 0;
	else
		modelParam.thetainvgamma = 1;
	modelParam.thetaprior[0] = modelParams[0].thetaGamma[0];
	modelParam.thetaprior[1] = modelParams[0].thetaGamma[1];	
    */

	/* get information for gene trees  */
  	nGene = numCurrentDivisions;  
  	gtree = (SPTree*)malloc(nGene*sizeof(SPTree)); 
  	FOR(i,nGene) 
  	{
	  	gtree[i].nTaxa = numTaxa;     
	  	gtree[i].nodes =(Treenode*)malloc((2*gtree[i].nTaxa-1)*sizeof(Treenode));
      		if(!gtree[i].nodes)
      		{ 
		   	printf(" allocating problem for gtree.nodes\n");
           		return(ERROR);
      		}
	  	gtree[i].taxaName = (char **)calloc(gtree[i].nTaxa , sizeof(char *));
	  	gtree[i].taxaName[0]=(char *)calloc(gtree[i].nTaxa*LSPNAME,sizeof(char));

	  	for(j = 0; j < gtree[i].nTaxa; j++) 
			gtree[i].taxaName[j] = gtree[i].taxaName[0]+j*LSPNAME;

	  	if(!gtree[i].taxaName)
	  	{
	   		printf(" allocating problem for spname\n");
	   		return(ERROR);
	  	} 
     
  	}

  	FOR(i,nGene)
      		FOR(j,gtree[i].nTaxa) 
			strcpy(gtree[i].taxaName[j], taxaNames[j]);

  	FOR(i,nGene) 
	{
		if (!strcmp(modelParams[i].ploidy,"Haploid"))
			gtree[i].haploid = 1;
		else if (!strcmp(modelParams[i].ploidy,"Zlink"))
			gtree[i].haploid = -1;
		else
			gtree[i].haploid = 0;
	}

  	/* allocate space for the species tree */
   	sptree.nSpecies = numTaxaSets;
   	sptree.nodes = (Treenode*)malloc((2*sptree.nSpecies-1)*sizeof(Treenode));
   	if(!sptree.nodes)
   	{
       		printf("allocating problem for sptree.nodes\n");
       		return(ERROR);
   	}

   	species = (int*)malloc(sptree.nSpecies*sizeof(int));
	sptree.nTaxa = 0;
 	FOR(i,sptree.nSpecies)  
        {
	   species[i] = 0;
	   FOR(k,numTaxa)
		{
		if (IsBitSet(k, taxaSet[i]) == YES)
			species[i]++;
		}
	  
	   sptree.nTaxa += species[i];
	}

	if(sptree.nTaxa != gtree[0].nTaxa)
	{
		printf("Wrong Taxaset!\n");
                return(ERROR);
        }

   	spnode = (int*)malloc(sptree.nTaxa*sizeof(int));

   	sptree.taxaName = (char **)calloc(sptree.nTaxa,sizeof(char *));
   	sptree.taxaName[0]=(char *)calloc(sptree.nTaxa*LSPNAME,sizeof(char));


	for(j = 0; j < sptree.nTaxa; j++) 
		sptree.taxaName[j] = sptree.taxaName[0]+j*LSPNAME;
	if(!sptree.taxaName)
	  {
	   	printf(" allocating problem for spname\n");
	   	return(ERROR);
	  } 

   	max = 0;
   	FOR(i,sptree.nSpecies) 
         	if(max<species[i]) 
			max = species[i];
   	max++;

   	sptree.speciesIndex = (int**)calloc(sptree.nSpecies,sizeof(int*));
   	if(!sptree.speciesIndex)
	   {
	   printf(" allocating problem for sptree.speciesIndex\n");
	   return(ERROR);
	   } 
   	sptree.speciesIndex[0] = (int*)calloc(sptree.nSpecies*max,sizeof(int));
   	FOR(i,sptree.nSpecies)
        	sptree.speciesIndex[i] = sptree.speciesIndex[0] + i*max;

   	FOR(i,sptree.nSpecies)
   	{
      		sptree.speciesIndex[i][0] = species[i];
		m = 1;
		FOR(k,numTaxa)
		{
			if(IsBitSet(k, taxaSet[i]) == YES)
			{
				spnode[k] = i;
				sptree.speciesIndex[i][m] = k;
				m++;
			}
		}

		/*FOR(k,sptree.speciesIndex[i][0])
			printf("index %d ",sptree.speciesIndex[i][k+1]);
		printf("\n");*/
   	}  
   	FOR(i,sptree.nTaxa)   
		strcpy(sptree.taxaName[i], taxaNames[i]);		

	sptree.nPop = 0;
   	FOR(i,sptree.nSpecies)
       		if(species[i]>1) 
			sptree.nPop++;
   	sptree.nPop += (sptree.nSpecies-1);

   	sptree.popIndex = (int*)malloc(sptree.nPop*sizeof(int));
   	j = 0;
   	FOR(i,sptree.nSpecies)
       		if(species[i]>1) 
          	{
               		sptree.popIndex[j] = i;
               		j++;
          	}

   	for(i=1;i<=sptree.nSpecies-1;i++)
       		sptree.popIndex[j++] = sptree.nSpecies-1+i;
   
   	sptree.nconstraint = sptree.nSpecies*(sptree.nSpecies-1)/2;

   	sptree.constraint = (double*)malloc(sptree.nconstraint*sizeof(double));
   	
	if(!sptree.constraint)
       	{
          	printf("allocating problem for sptree.constraint\n");
          	return(ERROR);
       	}

   	sptree.treeConstraint = (Distance*)malloc((sptree.nSpecies-1)*sizeof(Distance));
   	
	if(!sptree.treeConstraint)
	{
		printf("allocating problem for treecontraint\n");
		return(ERROR);
	}

  	/* allocating space for gene tree */

       	FOR(i,nGene)
        {


           gtree[i].taxaIndex = (int*)malloc(sptree.nTaxa*sizeof(int));
           if(!gtree[i].taxaIndex)
            {
                printf("allocating problem for gtree.taxaIndex\n");
                return(ERROR);
             }


           gtree[i].speciesIndex = (int**)calloc(sptree.nSpecies,sizeof(int*));
           if(!gtree[i].speciesIndex)
	   {
	   printf(" allocating problem for gtree[i].speciesIndex\n");
	   return(ERROR);
	   } 
           gtree[i].speciesIndex[0] = (int*)calloc(sptree.nSpecies*max,sizeof(int));
           FOR(j,sptree.nSpecies)
                gtree[i].speciesIndex[j] = gtree[i].speciesIndex[0] + j*max;


           gtree[i].popIndex = (int*)malloc(sptree.nPop*sizeof(int));
           if(!gtree[i].popIndex)
	   {
	   printf(" allocating problem for gtree[i].popIndex\n");
	   return(ERROR);
	   } 
        }

	FOR(i,nGene)
        	FOR(j,sptree.nTaxa)
         		gtree[i].taxaIndex[j] = j;

	FOR(i,nGene)
        {
 		gtree[i].nSpecies = sptree.nSpecies;
           	gtree[i].nPop = sptree.nPop;

		/*copy species index*/
           	FOR(j,sptree.nSpecies)
           	{
             		gtree[i].speciesIndex[j][0] = sptree.speciesIndex[j][0];
             		for(m=1;m<sptree.speciesIndex[j][0]+1;m++)           
             			gtree[i].speciesIndex[j][m] = sptree.speciesIndex[j][m];

             	}
		/*copy population index*/  
        	for(j=0;j<sptree.nPop;j++)
              		gtree[i].popIndex[j] = sptree.popIndex[j];      
     
        }
  	free(species); 
  	return(NO_ERROR);  
}

void PrintInf(FILE *fp)
{
  int i,j,k;
    
  /* write information into log.txt */
  fprintf(fp,"%s      Setting MCMC parameters\n",spacer);
  fprintf(fp,"%s      seed = %ld;\n",spacer, seed); 
  fprintf(fp,"%s      Theta = gamma(%lf, %lf); \n",spacer, modelParam.thetaprior[0],modelParam.thetaprior[1]);
  fprintf(fp,"\n\n\n");

  /* print Taxa information */
  fprintf(fp,"%s      Taxa Matrix\n",spacer);
  fprintf(fp,"%s          \t",spacer);
  FOR(i,sptree.nTaxa) fprintf(fp,"%s\t",sptree.taxaName[i]);
  fprintf(fp,"\n");

  FOR(j,nGene)
  { 
      fprintf(fp,"%s      Gene%d\t",spacer, j+1);
      FOR(i,sptree.nTaxa)
        fprintf(fp,"%d\t",gtree[j].taxaIndex[i]+1);
      fprintf(fp,"\n");
  }
  fprintf(fp,"\n\n");



  /*print species information */
  fprintf(fp,"%s      Species Matrix in Species SPTree\n",spacer);
  fprintf(fp,"%s      Species\tnTaxa\tList\n",spacer);
  FOR(i,sptree.nSpecies) 
      {
         fprintf(fp,"%s      %d.\t%d\t", spacer, i+1, sptree.speciesIndex[i][0]);
         FOR(j,sptree.speciesIndex[i][0])
          fprintf(fp,"%d\t",sptree.speciesIndex[i][j+1]+1);
         fprintf(fp,"\n");
      }
  fprintf(fp,"\n");
 
  FOR(j,nGene)
  {
      fprintf(fp,"%s      Species Matrix in Gene Tree%d\n",spacer, j+1);
      fprintf(fp,"%s      Species\tnTaxa\tList\n",spacer);
      FOR(i,sptree.nSpecies) 
       {
            fprintf(fp,"%s      %d.\t%d\t",spacer, i+1,gtree[j].speciesIndex[i][0]);
         FOR(k,gtree[j].speciesIndex[i][0])
          fprintf(fp,"%d\t",gtree[j].speciesIndex[i][k+1]+1);
         fprintf(fp,"\n");

        }
     fprintf(fp,"\n");
  }
  fprintf(fp,"\n\n");

  /*print population information*/
  fprintf(fp,"%s      Population Matrix in Species Tree\n",spacer);
  fprintf(fp,"%s      Species tree Populations:\t",spacer);
  FOR(i,sptree.nPop) 
      {
         fprintf(fp,"%d\t", sptree.popIndex[i]);
      }
  fprintf(fp,"\n");
   
  fprintf(fp,"%s      Gene tree Populations\n",spacer);
  FOR(j,nGene)
  {FOR(i,gtree[j].nPop) 
      {
         fprintf(fp,"%d\t", gtree[j].popIndex[i]);
      }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n\n");



}
  

/*
void FindMissingname(Tree *genetree)
{ 
  
	int i,j,k,index,stop,missingindex,totalname;
	char **name;

  
	totalname = 0;
	FOR(i, nGene) totalname += gtree[i].nTaxa;

	name = (char**)calloc(totalname*sizeof(char*));
  
	if(!name)
	{
	  printf("allocating problem for name\n");
	  return(ERROR);
	}

 	For(i,gtree[0].nTaxa) name[i] = gtree[0].taxaName[i];
 
	totalname = gtree[0].nTaxa;

	for(i=1;i<nGene;i++)
	{  
	 
		index = 0;
		FOR(j,gtree[i].nTaxa)
		{ 
		  find = 0;
          	for(k=0;k<totalname;k++)
		  {
                if(!strcmp(gtree[i].taxaName[j],name[k])) 
				{
					gtree[i].taxaIndex[index] = k;
					find = 1;
					index++;
					break;
				}
		  }               
		  if(find == 0) 
		  {			  
			  name[totalname] = gtree[i].taxaName[j]);
			  gtree[i].taxaIndex[index] = totalname;
			  index++;
			  totalname++;	             
		  }
		}
	}

	sptree.nTaxa = totalname;
	sptree.taxaName = (char **)malloc(sptree.nTaxa * sizeof(char *));
    	if(!sptree.taxaName)
	{
	  printf(" allocating problem for sptree.taxaName\n");
	  return(ERROR);
	} 
	for(i=0;i<sptree.nTaxa;i++)
		sptree.taxaName[i] = name[i];

	for(i=0;i<nGene;i++)
	{
		index = 0;
		gtree[i].nmissTaxa = sptree.nTaxa - gtree[i].nTaxa;
		gtree[i].mtaxaIndex = (int*)malloc(gtree[i].nmissTaxa*(int));
		for(j=0;j<sptree.nTaxa-1;j++)
		{
			k=0;
			if(j == gtree[i].taxaIndex[k])
			{
				k++;
				continue;
			}
			else
			{
				gtree[i].mtaxaIndex[index];
				index++;
			}
	}

        free(name);
	return(NO_ERROR);  
}
*/

void ReadGeneTree(FILE *fTree)
{
  int i;  

  /* Read in Gene Trees */
  FOR(i,nGene) ReadaTree(fTree,&gtree[i]);

 
  /* Calculate the age of nodes */
  FOR(i,nGene) CalNodeAge(gtree[i].root, &gtree[i]);


}
  
int SPPrintTree(int curGen, SPTree *tree, int showBrlens, int showTheta, int showMu, int isRooted)
{
#if 0
    char			*tempStr;
	int             tempStrSize, i;
	FILE			*sumt;
	char			*name, spfilename[100];

	/* allocate the print string */
	printStringSize = 200;
	printString = (char *)SafeMalloc((size_t) (printStringSize * sizeof(char)));
	if (!printString)
		{
		SPMrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
		return (ERROR);
		}
	*printString = '\0';

	tempStrSize = 200;
	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
		{
		SPMrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
		return (ERROR);
		}

	/*print tree file header*/
	if (curGen == 1)
	{
		/* print #NEXUS and translation block information */
		sprintf (tempStr, "#NEXUS\n");
		if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
		sprintf (tempStr, "[ID: %ld]\n", seed);
		if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
		sprintf (tempStr, "begin trees;\n");
		if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
		sprintf (tempStr, "   translate\n");
		if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
		for (i=0; i<sptree.nSpecies-1; i++)
			{
			name = taxaSetNames[i];
			sprintf (tempStr, "      %2d %s,\n", i+1, name);
			if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
			}
		name = taxaSetNames[i];
		sprintf (tempStr, "      %2d %s;\n", i+1, name);
        if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);

		/*create .sumt file*/
# if defined (MPI_ENABLED)
		if(proc_id == 0)
		{
# endif
		sprintf(spfilename, "%s.sumt", chainParams.chainFileName);
		sumt = fopen(spfilename, "w");
		sprintf(spfilename, "%s.sptree", chainParams.chainFileName);
	
		/* print #NEXUS and translation block information */
		sprintf(spfilename, "%s.sumt", chainParams.chainFileName);
		sumt = fopen(spfilename, "w");
		sprintf(spfilename, "%s.sptree", chainParams.chainFileName);
		fprintf (sumt, "#NEXUS\n");
		fprintf (sumt, "Begin data;\n");
		fprintf (sumt, "Dimensions ntax = %d nchar = 1;\n", sptree.nSpecies);
		fprintf (sumt, "FORMAT DATATYPE = dna gap=-  MISSING=?  interleave=yes;\nMatrix\n");
		for (i=0; i<sptree.nSpecies; i++)
			{
			name = taxaSetNames[i];
			fprintf (sumt, " %s A\n", name);
			}
		fprintf (sumt, "\n;\nEND;\n");
		fprintf (sumt, "Begin mrbayes;\n  prset best = 1;\n");
		fprintf (sumt, "  sumt burnin = %d nruns = %d filename = %s contype = allcompat;\nEND;\n", chainParams.numGen/chainParams.sampleFreq/2, chainParams.numRuns, spfilename);
		fclose (sumt);
#if defined (MPI_ENABLED)
		}
# endif
	}

	/* write the tree in Newick format */
	sprintf (tempStr, "   tree rep.%d = ", curGen);
	if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);

    SPSaveSprintf (&tempStr, &tempStrSize,"(");
	SPAddToPrintString (tempStr);
					
    SPWriteTreeToFile (tree, tree->root, showBrlens, showTheta, showMu, isRooted);

    if(showTheta == YES) 
		SPSaveSprintf (&tempStr, &tempStrSize,")[#%lf];\n",tree->nodes[tree->root].theta);
    else 
		SPSaveSprintf (&tempStr, &tempStrSize,");\n");
	SPAddToPrintString (tempStr);
	free (tempStr); 

#endif

	return (NO_ERROR);					
}

void SPWriteTreeToFile (SPTree *tree, int inode, int showBrlens, int showTheta, int showMu, int isRooted)

{
#if 0
        char			*tempStr;
		int                      tempStrSize = 200;

		tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
		if (!tempStr)
			SPMrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));


			
		if (tree->nodes[inode].nson == 0)
			{
				if (showBrlens == YES)
				{
    				SPSaveSprintf (&tempStr, &tempStrSize, "%d:%lf", inode+1, tree->nodes[inode].brlens);
					SPAddToPrintString (tempStr);
					if((tree->nodes[inode].theta>0) && showTheta == YES) 
					{
						SPSaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[inode].theta);
						SPAddToPrintString (tempStr);
					}
				}
				else
				{
					SPSaveSprintf (&tempStr, &tempStrSize, "%d", inode+1);
					SPAddToPrintString (tempStr);
				}
			}
		else
			{
				if (inode != tree->root)
				{
					SPSaveSprintf (&tempStr, &tempStrSize, "(");
					SPAddToPrintString (tempStr);
				}
				SPWriteTreeToFile (tree,tree->nodes[inode].sons[0],  showBrlens, showTheta, showMu, isRooted);
				SPSaveSprintf (&tempStr, &tempStrSize, ",");
				SPAddToPrintString (tempStr);
				SPWriteTreeToFile (tree,tree->nodes[inode].sons[1], showBrlens, showTheta, showMu, isRooted);	
				if (inode != tree->root)
				{
					if (tree->nodes[inode].father == tree->root && isRooted == NO)
					{
						if (showBrlens == YES)
						{
							SPSaveSprintf (&tempStr, &tempStrSize, ",%d:%lf", tree->nodes[inode].father + 1, tree->nodes[tree->nodes[inode].father].brlens);
							SPAddToPrintString (tempStr);
							if((tree->nodes[tree->nodes[inode].father].theta>0) && showTheta == YES) 
							{
								SPSaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[tree->nodes[inode].father].theta);
								SPAddToPrintString (tempStr);
							}
						}
						else
						{
							SPSaveSprintf (&tempStr, &tempStrSize, ",%d", tree->nodes[inode].father + 1);
							SPAddToPrintString (tempStr);
						}
					}
				
					if (showBrlens == YES && isRooted == YES) /*tree->nodes[inode].father != tree->root)*/
					{
						SPSaveSprintf (&tempStr, &tempStrSize,"):%lf", tree->nodes[inode].brlens);
						SPAddToPrintString (tempStr);
						if((tree->nodes[inode].theta > 0) && showTheta == YES)
						{
							SPSaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[inode].theta);
							SPAddToPrintString (tempStr);
						}
					}
					else
					{
						SPSaveSprintf (&tempStr, &tempStrSize, ")");
						SPAddToPrintString (tempStr);
					}					
				}
			}
	free (tempStr);
#endif
}


double SpeciesDistance(int inode, int jnode, SPTree *genetree)
{
  int i,j;
  double distance;
  double con = 0.0;

  if(genetree->speciesIndex[inode][0] == 0 || genetree->speciesIndex[jnode][0] == 0)
     return(NA);

  FOR(i,genetree->speciesIndex[inode][0])
     FOR(j,genetree->speciesIndex[jnode][0])
          {      
                 distance = NodeDistance(genetree,genetree->speciesIndex[inode][i+1],genetree->speciesIndex[jnode][j+1]);
                 if(j==0&&i==0) 
                 {
                     con = distance;
                     continue;
                 }
                 if(con>distance) con = distance;
           }
  return(con);
}

double NodeDistance(SPTree *tree, int inode, int jnode)
{
  int i, *ancestor,father,nancestor,stop=0;

  ancestor=(int*)malloc(tree->nTaxa*sizeof(int));

  i=0;
  father = inode;
  do{
      father = tree->nodes[father].father;
      ancestor[i] = father;
      i++;
    }while( father != tree->root);

  nancestor = i;

  father = jnode;
  do{
      father = tree->nodes[father].father;
      FOR(i,nancestor) 
         if(ancestor[i] == father) 
             {
                stop = 1;
                break;
              }
     }while(stop == 0);
  free(ancestor);
  return(2*tree->nodes[father].age);
}
 
int Constraint(SPTree *genetree, int numgenetree, SPTree *speciestree, Distance *constraint)
{
  	int i, j, w, k;
  	double distance;

  	w = 0;
  	FOR(j,(speciestree->nSpecies-1))
    		for(k=j+1;k<speciestree->nSpecies;k++)
     		{
			constraint[w].nodes[0] = j;
             	constraint[w].nodes[1] = k;
		 	constraint[w].dist = 100000.0;
      		FOR(i,numgenetree)
          		{  
             		distance = SpeciesDistance(j,k,&genetree[i]);
             		if(constraint[w].dist > distance) constraint[w].dist = distance;
          		}
      		if(constraint[w].dist <= 0.0 || constraint[w].dist >= 100000.0)
          		{
             		printf("constraint%d is wrong\n",w);
             		return(ERROR);
           		}
       		w++;
     		}
  	return(NO_ERROR);
}

double CalNodeAge(int node, SPTree *tree)
{ 
  if(tree->nodes[node].nson == 0) 
      {
         tree->nodes[node].age = 0.0;
         return(0.0);
      }
  else
      {
         tree->nodes[node].age = 0.5*(CalNodeAge(tree->nodes[node].sons[0],tree) + tree->nodes[tree->nodes[node].sons[0]].brlens + CalNodeAge(tree->nodes[node].sons[1],tree) + tree->nodes[tree->nodes[node].sons[1]].brlens);
         return( tree->nodes[node].age );
      }
}
       
int StartSptree(SPTree *speciestree, int numchange) { //*speciestree is the address of the global sptree
   int i, j;
   Distance* onetreeConstraint;
   
   onetreeConstraint = (Distance *) calloc (speciestree->nSpecies, sizeof(Distance));

	//CNKA 10/10 we don't want all cophenetic distances for all genes, just the min across genes
	speciestree->mindist = (double*)malloc(speciestree->nSpecies*speciestree->nSpecies*sizeof(double));
   GetMinDists(gtree, nGene, &speciestree->mindist);
   //GetConstraints(speciestree, onetreeConstraint);

	/*jiggle polytomys*/
  	for(i=0;i<sptree.nSpecies-2;i++)
		for(j=i+1;j<sptree.nSpecies-1;j++) {
			if(onetreeConstraint[i].dist == onetreeConstraint[j].dist) {
				onetreeConstraint[j].dist -= 1e-10;
				//if(onetreeConstraint[j].dist < 0.0) onetreeConstraint[j].dist = 0.0;
			}
		}

	//change constraints
	for(i=0; i<numchange; i++)
		ChangeConstraint(onetreeConstraint, speciestree->nSpecies-1);

  	//propose a species tree
  	MaximumTree(onetreeConstraint, speciestree, speciestree->nSpecies-1);

  	// initial theta
  	FOR(i, 2*speciestree->nSpecies-1) speciestree->nodes[i].theta = -1.0;
  	FOR(i, speciestree->nPop) speciestree->nodes[speciestree->popIndex[i]].theta = 0.1;

    free (onetreeConstraint);

    return(NO_ERROR); 
}

int CheckConstraint(SPTree *genetrees, int ngene, Distance *constraint, int nconstraints) {
	int i, j, k, l, node0, node1;
	double dist;

	FOR(i, nconstraints)
	{
		node0 = constraint[i].nodes[0];
		node1 = constraint[i].nodes[1];
		FOR(j, ngene)
		{
			FOR(k, sptree.speciesIndex[node0][0])
				FOR(l,sptree.speciesIndex[node1][0])
				{
					dist = NodeDistance(&genetrees[j],sptree.speciesIndex[node0][k+1],sptree.speciesIndex[node1][l+1]);
					if(dist <  constraint[i].dist)
					{
						printf("constrain %d %d %f genedist %d %d %d %f\n",node0,node1,constraint[i].dist,j,sptree.speciesIndex[node0][k+1],sptree.speciesIndex[node1][l+1],dist);
						return ERROR;
					}
				}
		}
	}
	return (NO_ERROR);
}

int ChangeConstraint(Distance *dist, int nconstraints)
{
	int i;
	double window=0.0001;

	i = (int)(RandomNumber(&globalSeed)*nconstraints);
	dist[i].dist -= RandomNumber(&globalSeed)*window;
	
  	if(dist[i].dist < 1e-10) 
     	{ 
      	dist[i].dist = 1e-10;
     	}
   	return(NO_ERROR);
}


int MaximumTree(Distance *dist, SPTree *speciestree, int nconstraints)
{
	int i, j, k, father, son;

	/*sort the distances*/	
	quick_struct(dist, nconstraints);

	for(i=0; i<2*(nconstraints+1)-1; i++)
	{
		speciestree->nodes[i].father = -2;
	}

	k = nconstraints + 1;
	for(i=0; i<nconstraints; i++)
	{
		speciestree->nodes[k].nson = 2;
		speciestree->nodes[k].age = dist[i].dist/2;
		for(j=0; j<2; j++)
		{
			speciestree->nodes[dist[i].nodes[j]].age = 0.0;
			speciestree->nodes[dist[i].nodes[j]].nson = 0;
			speciestree->nodes[dist[i].nodes[j]].sons[0] = -2;
			speciestree->nodes[dist[i].nodes[j]].sons[1] = -2;

			father = speciestree->nodes[dist[i].nodes[j]].father;
			son = dist[i].nodes[j];
			while(father != -2)
			{
				son = father;
				father = speciestree->nodes[father].father;
			} 
			speciestree->nodes[k].sons[j] = son;
			speciestree->nodes[son].father = k;
			speciestree->nodes[son].brlens = speciestree->nodes[k].age - speciestree->nodes[son].age; 
		}
		k++;
	}
	speciestree->root = k-1;
	speciestree->nodes[k-1].father = -1;
	return(NO_ERROR);
}

int SPTreeConstraint(Distance *minimumdistance, Distance *distance, long int nconstraints, int nspecies) {
	int i, k;
	int node[NSPECIES], index[2];
   long int j;
	
	quick_struct(minimumdistance, nconstraints);

  	distance[0].dist = minimumdistance[0].dist;
  	distance[0].nodes[0]=node[0]=minimumdistance[0].nodes[0];
  	distance[0].nodes[1]=node[1]=minimumdistance[0].nodes[1];
  
  	for(i=2; i<nspecies; i++)
  	{ 
	  	FOR(j, nconstraints)  
	  	{ 
		  	index[0]=0;index[1]=0;
		  	FOR(k,i)
		  	{
             		if (minimumdistance[j].nodes[0] == node[k]) {index[0]++;}
             		if (minimumdistance[j].nodes[1] == node[k]) {index[1]++;}
		  	}

		  	if((index[0]==0)&&(index[1]==1)) 
		  	{
			   	distance[i-1].nodes[0]=minimumdistance[j].nodes[1];
               distance[i-1].nodes[1]=minimumdistance[j].nodes[0];
               distance[i-1].dist=minimumdistance[j].dist;
               node[i]=minimumdistance[j].nodes[0];
               break;
		  	}        
		  	if((index[0]==1)&&(index[1]==0))
		  	{
			   	distance[i-1].nodes[0]=minimumdistance[j].nodes[0];
               distance[i-1].nodes[1]=minimumdistance[j].nodes[1];
               distance[i-1].dist = minimumdistance[j].dist;
               node[i] = minimumdistance[j].nodes[1];
               break;
		  	}         
	  	}  
  	}
	return(NO_ERROR);
}

//loads constraints on the species tree given that s->mindist is already set
int GetConstraints (double **distanceMatrix, Distance *constr) {

	int trace=0;
	int i,j,k=0,nsp=numSpecies;
	int *node, index[2];
    Distance *minimumdistance;

    node = (int *) calloc (nsp, sizeof(int));
    minimumdistance = (Distance *) calloc (nsp*(nsp-1)/2, sizeof(Distance));

	//convert s->mindist to a distance structure
	for(i=0; i<nsp; i++) for(j=i+1;j<nsp;j++) {
		minimumdistance[k].dist=distanceMatrix[0][i*nsp+j];
		minimumdistance[k].nodes[0]=i;
		minimumdistance[k++].nodes[1]=j;
	}
	/*for(i=0; i<nsp*(nsp-1)/2; i++) {
		printf("min %lf ",minimumdistance[i].dist);
	}*/
	quick_struct(minimumdistance, nsp*(nsp-1)/2);

  	constr[0].dist = minimumdistance[0].dist;
  	constr[0].nodes[0]=node[0]=minimumdistance[0].nodes[0];
  	constr[0].nodes[1]=node[1]=minimumdistance[0].nodes[1];

  	for(i=2; i<nsp; i++)
  	{
	  	FOR(j, nsp*(nsp-1)/2)
	  	{
		  	index[0]=index[1]=0;
		  	FOR(k,i)
		  	{
       		if (minimumdistance[j].nodes[0] == node[k]) {index[0]++;}
       		if (minimumdistance[j].nodes[1] == node[k]) {index[1]++;}
		  	}

		  	if((index[0]==0)&&(index[1]==1)) {
		   	constr[i-1].nodes[0]=minimumdistance[j].nodes[1];
	         constr[i-1].nodes[1]=minimumdistance[j].nodes[0];
	         constr[i-1].dist=minimumdistance[j].dist;
	         node[i]=minimumdistance[j].nodes[0];
	         break;
		  	}
		  	if((index[0]==1)&&(index[1]==0)) {
		   	constr[i-1].nodes[0]=minimumdistance[j].nodes[0];
            constr[i-1].nodes[1]=minimumdistance[j].nodes[1];
            constr[i-1].dist = minimumdistance[j].dist;
            node[i] = minimumdistance[j].nodes[1];
            break;
		  	}
	  	}
  	}
  	if(trace) for(i=0;i<nsp-1;i++)
	  printf("\n%d. (%d.%d:%f)",i+1,constr[i].nodes[0],constr[i].nodes[1],constr[i].dist);

    free (minimumdistance);
    free (node);
    return(NO_ERROR);
}


int FindaPosition(int nodenumber,int root,double currentdistance,double *position,SPTree *speciestree)  /*ancestor contains nodes and distance */
{
  

  position[0]=0.0;
  position[1]=0.0;

  while(nodenumber!=root)
  { 
	  position[1] += speciestree->nodes[nodenumber].brlens;
	  if (currentdistance<=position[1])
	  {
		  position[0]=nodenumber;
		  position[1]=currentdistance-(position[1]-speciestree->nodes[nodenumber].brlens);
		  break;
	  } 
	  if(currentdistance == position[1])
	  {
		  printf("Equal distance for position %f\n",currentdistance);
		  return(ERROR);
	  }
	  nodenumber=speciestree->nodes[nodenumber].father;
  }
  if(nodenumber==root)
  {
	  position[0]=root;
	  position[1]=currentdistance-position[1];
  }
  
  return(NO_ERROR);
}

void quick_struct(Distance item[],int count)
         {qs_struct(item,0,count-1);}

void qs_struct(Distance item[],int left,int right)
{

   register int i, j;
   double x;
   Distance temp;

  i=left; j=right;
  x=item[(left+right)/2].dist;


  do {
       while(item[i].dist<x && i<right) i++;
       while(item[j].dist>x && j>left) j--;

       if(i<=j)  {
          temp=item[i];
          item[i]=item[j];
          item[j]=temp;
          i++; j--;
          }
     } while(i<=j);


    if(left<j)  qs_struct(item,left, j);
    if(i<right)  qs_struct(item,i, right);
}

int SPLogLike(SPTree *genetree, SPTree *speciestree, double *lnl)
{
  	int k;
  	double lnp;

  	*lnl = 0.0;

  	FOR(k,nGene) 
    	{
       	if( LnLikehood1Tree(&genetree[k],speciestree, &lnp) == ERROR)
          	{
            	printf("Error for calculating the log likelhood of gene trees\n");
               	return(ERROR);
          	}
       	*lnl += lnp;
      }
  	return(NO_ERROR);
}

int LnLikehood1Tree(SPTree *genetree, SPTree *speciestree, double *lnp)
{ 
   	int i,j, ipop,k=0, nt;
   	double *t,y,theta;
   	CoalTime coal,coal_t;

  	if( GetCoaltime(genetree,speciestree,&coal) == ERROR)
     	{
         printf(" bugs in GetCoaltime \n");
         return(ERROR);
     	}
 
  	i=0;
  	GetNcoal(speciestree->root, &coal, genetree, speciestree, &coal_t, &i);
   
  	for(ipop=0; ipop<genetree->nPop; ipop++)
  	{      
	  nt=coal_t.ncoal[ipop]; 
	  t=coal_t.tj[ipop];   
	  for(i=0; i<nt-1; i++)  for (j=i+1; j<nt; j++)
	  {
            if (t[j]<t[i])  
			{ 
				y=t[i]; 
				t[i]=t[j];
				t[j]=y;
			}
	  }
  	}

  	*lnp = 0.0;
 	for(i=0; i<genetree->nPop; i++) 
	{
      	ipop=i; t=coal_t.tj[ipop];
         	/* handle the haploid gene*/
         	if(genetree->haploid == 1) /*haploid*/
              	theta = 0.25*speciestree->nodes[coal_t.nodes[ipop]].theta;
         	else if (genetree->haploid == -1) /*zlink*/
				theta = 0.75*speciestree->nodes[coal_t.nodes[ipop]].theta;
			else  /*diploid*/
              	theta = speciestree->nodes[coal_t.nodes[ipop]].theta;

#	if defined (DEBUG)
         	if(theta<1.e-6)
            {
               printf("theta is wrong for node %d  theta %lf\n",coal_t.nodes[ipop],theta);
               return(ERROR);
            }
#	endif
 
         	if(coal_t.nin[ipop]>coal_t.nout[ipop]) 
		{
         		*lnp+=(coal_t.nin[ipop]-coal_t.nout[ipop])*log(2/theta);
         		for(j=coal_t.nin[ipop],k=0; j>coal_t.nout[ipop]; j--,k++) 
			{
            	y = t[k]- (k==0? 0:t[k-1]);       
            	*lnp -= j*(j-1)/theta*y;
         		}
      	}
      	if(coal_t.nout[ipop]>1) 
		{
         	y = (coal_t.nout[ipop]==coal_t.nin[ipop] ? speciestree->nodes[coal_t.nodes[ipop]].brlens : speciestree->nodes[coal_t.nodes[ipop]].brlens-t[k-1]);
         	*lnp -= coal_t.nout[ipop]*(coal_t.nout[ipop]-1)/theta*y;
         
      	}
   	}  
   	return(NO_ERROR);
}

int SPLogLike_invgamma(SPTree *genetree, SPTree *speciestree, double *lnl)
{
  	int i,j,k;
  	double a[NSPECIES],b[NSPECIES], totalcoal=0.0, alpha=modelParam.thetaprior[0],beta=modelParam.thetaprior[1];
  

  	/*initial a and b*/
  	for(i=0; i<speciestree->nPop;i++)
	{
		a[i] = 0.0;
 		b[i] = 0.0;
	}

  	*lnl = 0.0;

  	FOR(k,nGene) {
		if( LnLikehood1Tree_invgamma(&genetree[k],speciestree, a,b) == ERROR) {
      	printf("Error for calculating the log likelhood of gene trees\n");
      	return(ERROR);
    	}
   }

   for(i=0; i<speciestree->nPop;i++){
     	*lnl += alpha*log(beta)-(a[i]+alpha)*log(b[i]+beta);
		totalcoal += a[i];
		for(j=0;j<a[i];j++) *lnl += log(alpha+j);
		speciestree->nodes[sptree.popIndex[i]].theta = (b[i]+beta)/(alpha+a[i]-1); 
	}
  	if(totalcoal != (speciestree->nTaxa-1)*nGene)
  	{
		printf("the number of coalescence is %f\n",totalcoal);
		return(ERROR);
  	}
  	return(NO_ERROR);
}

int LnLikehood1Tree_invgamma(SPTree *genetree, SPTree *speciestree, double *a, double *b)
{ 
   	int i,j, ipop, k=0,nt, popnode=0;
   	double *t,y;
   	CoalTime coal,coal_t;

  	if( GetCoaltime(genetree,speciestree,&coal) == ERROR)
     	{
         printf(" bugs in GetCoaltime \n");
         return(ERROR);
     	}
  
  	i=0;
  	GetNcoal(speciestree->root, &coal, genetree, speciestree, &coal_t, &i);
  
  	for(ipop=0; ipop<genetree->nPop; ipop++)
 	{      
	  nt=coal_t.ncoal[ipop]; 
	  t=coal_t.tj[ipop];   
	  for(i=0; i<nt-1; i++)  for (j=i+1; j<nt; j++)
	  {
            if (t[j]<t[i])  
			{ 
				y=t[i]; 
				t[i]=t[j];
				t[j]=y;
			}
	  }
  	}

 	for(i=0; i<genetree->nPop; i++) 
	{
      	ipop=i; t=coal_t.tj[ipop];
    
	   	/* find the population node*/
	   	for(j=0;j<sptree.nPop;j++)
		{
			if(sptree.popIndex[j] == coal_t.nodes[i])
			{
				popnode = j;
				break;
			}
		}

         if(coal_t.nin[ipop]>coal_t.nout[ipop]) {
			 a[popnode] += (coal_t.nin[ipop]-coal_t.nout[ipop]);

         for(j=coal_t.nin[ipop],k=0; j>coal_t.nout[ipop]; j--,k++) {
            y = t[k]- (k==0? 0:t[k-1]);         
		if(genetree->haploid == 1) /*haploid*/
		  	 b[popnode] += 4*j*(j-1)*y;
		else if (genetree->haploid == -1) /*zlink*/
			b[popnode] += 1.3333333*j*(j-1)*y;
		else /*diploid*/
			 b[popnode] += j*(j-1)*y;

         }
      }
      if(coal_t.nout[ipop]>1) {
         y = (coal_t.nout[ipop]==coal_t.nin[ipop] ? speciestree->nodes[coal_t.nodes[ipop]].brlens : speciestree->nodes[coal_t.nodes[ipop]].brlens-t[k-1]);
         if(genetree->haploid == 1) /*haploid*/
		  	 b[popnode] += 4*coal_t.nout[ipop]*(coal_t.nout[ipop]-1)*y;
		 else if(genetree->haploid == -1) /*zlink*/
		  	 b[popnode] += 1.3333333*coal_t.nout[ipop]*(coal_t.nout[ipop]-1)*y;
		else /*diploid*/
			 b[popnode] += coal_t.nout[ipop]*(coal_t.nout[ipop]-1)*y;
         
      }
   }
   return(NO_ERROR);
}

int GetCoaltime(SPTree *genetree,SPTree *speciestree, CoalTime *coal)
{  
	int inode,i,k,inode_gene,inode_sp,stop,index=0;
	int *genetreenodes;
	double distance_gene;

	genetreenodes = (int*)malloc(2*genetree->nTaxa*sizeof(int));
	if(!genetreenodes)
	{
		printf("allocating problem for genetreenodes\n");
		return(ERROR);
	}
	FOR(i,genetree->nSpecies)
      { 
        	FOR(inode,genetree->speciesIndex[i][0])
        	{        
        		inode_gene=genetree->speciesIndex[i][inode+1];      
        		stop=0;        
        		do{
				inode_gene = genetree->nodes[inode_gene].father; 
				FOR(k,index)
              			if(inode_gene==genetreenodes[k]) 
					{
						stop=1;
						break;
					}
           			if(stop==1) break;

				/*find the node in the species tree that is right below the inode_gene*/
				distance_gene = genetree->nodes[inode_gene].age;
            		inode_sp = FindspNodeBelow(speciestree,i,distance_gene);
				/*record the coalescence*/
				coal->ncoal[index]=1;
                 		coal->nodes[index]=inode_sp;
                 		coal->tj[index][0]= distance_gene - speciestree->nodes[inode_sp].age;
                 		genetreenodes[index]=inode_gene;
                 		index++;					
           		} while(inode_gene!=genetree->root);
		}
    	}
	
    	/* check if all the Taxa coalesce */
    	if(index != genetree->nTaxa-1)
    	{
       	printf("Error for the coalescence time calculation\n");
       	return(ERROR);
    	}
     	free(genetreenodes);
	return(NO_ERROR);
}

int FindspNodeBelow(SPTree *speciestree, int spnode, double dis_gene)
{
	int inode=spnode,father;
	double dis_sp;

	father = speciestree->nodes[inode].father;
	dis_sp = speciestree->nodes[father].age;
	while(dis_gene >= dis_sp)
	{
		if(father == speciestree->root)
		{
			inode = speciestree->root;
			break;
		}
		inode = speciestree->nodes[inode].father;
		father = speciestree->nodes[inode].father;
		dis_sp = speciestree->nodes[father].age;
	}
	return(inode);
}


int GetNcoal(int inode, CoalTime *coal, SPTree *genetree, SPTree *speciestree, CoalTime *coal_t,int *index)
{
  	int k,i,kson,w;

  	if(inode<speciestree->nSpecies)   
  	{
     		if(genetree->speciesIndex[inode][0] == 1 )    
        		return(1);
    		else
       	{
         		coal_t->nodes[*index] = inode;
         		coal_t->nin[*index] = genetree->speciesIndex[inode][0];
         		w=0;
         		FOR(i,genetree->nTaxa-1)
         		{
           			if(coal->nodes[i]==inode) 
				{  
					coal_t->tj[*index][w]=coal->tj[i][0];
					w++;
				}
	  		}
  
	  		coal_t->ncoal[*index]=w;
         		w=*index;
	  		(*index)++;
         		coal_t->nout[w] = coal_t->nin[w] - coal_t->ncoal[w];
         		return(coal_t->nout[w]);
      	}
  	}
  	else
  	{
	 	w=0;
	 	FOR(i,genetree->nTaxa-1) 
	 	{
        		if(coal->nodes[i]==inode) 
			{  
				coal_t->tj[*index][w]=coal->tj[i][0];
				w++;
			}
	 	}
  
	 	coal_t->ncoal[*index]=w;
	 	coal_t->nodes[*index]=inode;
	 	w=*index;
	 	(*index)++;
	 	coal_t->nin[w]=0;
	 	FOR(k,speciestree->nodes[inode].nson)
	 	{
		 	kson=speciestree->nodes[inode].sons[k];
		 	coal_t->nin[w]+=GetNcoal(kson,coal,genetree,speciestree,coal_t,index);
	 	}
	 	coal_t->nout[w]=coal_t->nin[w]-coal_t->ncoal[w];
 
	 	return(coal_t->nout[w]);
  	}
}

double SPLnP1 (double t, double l, double m, double r) {
	double		p0t;
	p0t = r*(l-m) / (r*l + (l*(1.0-r)-m)*exp((m-l)*t) );
	return (log(1.0/r) + 2.0*log(p0t) + (m-l)*t);
}

double SPLnVt (double t, double l, double m, double r) {
	double		p0t;
	p0t = r*(l-m) / (r*l + (l*(1.0-r)-m)*exp((m-l)*t) );
	return (log(1.0 - (1.0/r) * p0t * exp((m-l)*t)));
}

int SPLnBirthDeathPriorPr (double *prob, double sR, double eR, double sF,SPTree *speciestree) {
	int		i, j;
	double	*nt,rootTime;
   nt = (double*)malloc((speciestree->nSpecies-1)*sizeof(double));
   if(!nt) {
   	printf("allocating problem for nt\n");
      return(ERROR);
   }
	// rescale all of the node times on the tree
	rootTime = speciestree->nodes[speciestree->root].age;
   j = 0;
   for (i=speciestree->nSpecies; i<2*speciestree->nSpecies-1; i++)
   	if(i != speciestree->root) nt[j++] = speciestree->nodes[i].age;
   nt[speciestree->nSpecies-2] = rootTime;
	for (i=0; i<speciestree->nSpecies-1; i++) {
		nt[i] /= rootTime;
      if(nt[i]>1.0) {
			printf("nt time error\n");
			return(ERROR);
		}
	}
	/* I think this is correct. It looks as if Yang and Rannala (1997)
	   have the root time constrained to be 1.0. */
	rootTime = 1.0;
							
	/* calculate probabilities of tree */
	if ((sR - eR)> 1e-8) {
		(*prob) = (speciestree->nSpecies - 2.0) * log(sR);
		for (i=0; i<speciestree->nSpecies-2; i++)
			(*prob) += SPLnP1 (nt[i], sR, eR, sF) - SPLnVt (rootTime, sR, eR, sF);
	} else {
		(*prob) = 0.0;
		for (i=0; i<speciestree->nSpecies-2; i++)
			(*prob) += log (1.0 + sF * eR) - (2.0 * log(1.0 + sF * eR * nt[i]));
	}

	/*plus the probability for the root*/
	*prob -= modelParam.treeHeightExp*speciestree->nodes[speciestree->root].age;
	free(nt);
	return (NO_ERROR);
}

int FindConstraint(int node, int *anode)
{
	int i;

	if(node == sptree.treeConstraint[0].nodes[0])
	{
		*anode = sptree.treeConstraint[0].nodes[1];
		return(0);
	}
	
	for(i=0;i<sptree.nSpecies;i++)
		if(node == sptree.treeConstraint[i].nodes[1])
			break;
	*anode = sptree.treeConstraint[i].nodes[0];
	return(i);
}

MrBFlt Prob_sptree(Distance *tau,int ntime)
{
  int i,j,round=100;
  MrBFlt time,*lnp,lamda,roottime,randomroottime,max,lnprior2=0.0;
 
  lnp = (MrBFlt *)malloc((size_t) ((round) * sizeof(MrBFlt)));
		
  if (!lnp)
			{
			printf ("%s   Problem allocating lnp \n", spacer);
			exit(-1);
			}

 
  /* find the root depth and save it into tau[ntime-1]*/
  for(i=(ntime-2);i>=0;i--)
  {
	  if(tau[ntime-1].dist < tau[i].dist)
	  {
		  time = tau[ntime-1].dist;
		  tau[ntime-1].dist = tau[i].dist;
		  tau[i].dist = time;
	  }
  }

 
  for(i=0;i<round;i++)
   {
     lnp[i] = 0.0;
     lamda = (double)(i+1)/10;
    
     /*generate a random roottime */
     roottime = tau[ntime-1].dist*0.5;
     randomroottime = (1-exp(-roottime))*RandomNumber(&globalSeed);
     
     /* given the randomroottime, calculate the density function */
     for(j=0;j<ntime-1;j++)
       {
         time = tau[j].dist/2;
         if(time < randomroottime) lnp[i] += log(1-exp(-lamda*time))-log((1-exp(-lamda*randomroottime)));
        
       }  
   }
   
   /*find the maximum*/
   max = lnp[0];
 
 
   for(i=1;i<round;i++) if(max<lnp[i]) max = lnp[i];
   for(i=0;i<round;i++) if(max-lnp[i]<=10) 
   {
	   lnprior2 += exp(lnp[i]-max);
   }   

   /*log prior is equal to the log of the average*/
   lnprior2 = -log(round) + max + log(lnprior2);

   /* plus the probability of the roottime*/    
   lnprior2 += log(1-exp(-roottime)); 
   
   free(lnp);
  
   return(lnprior2);
}

#define TARGETLENDELTA (100)


int SPSaveSprintf(char **target, int *targetLen, char *fmt, ...) {
  va_list    argp;
  int        len,retval;

  va_start(argp, fmt);
#ifdef VISUAL
  len = _vsnprintf(NULL, 0, fmt, argp);
#else
  len = vsnprintf(NULL, 0, fmt, argp);
#endif

  va_end(argp);

  if(len>*targetLen)
        {
/*        fprintf(stderr, "* increasing buffer from %d to %d bytes\n", *targetLen, len+TARGETLENDELTA); */
        *targetLen = len+TARGETLENDELTA; /* make it a little bigger */
            *target = (char*)realloc(*target, *targetLen);
        }

  va_start(argp,fmt);
  retval=vsprintf(*target, fmt, argp);
  va_end(argp);

/*   fprintf(stderr, "* savesprintf /%s/\n",*target); */
  return retval;
}

#endif
