/*
 *  MrBayes 3
 *
 *  (c) 2002-2013
 *
 *  This file originally contributed by:
 *
 *  Marc A. Suchard
 *  Department of Biomathematics
 *  University of California, Los Angeles
 *  Los Angeles, CA 90095
 *  msuchard@ucla.edu
 *
 *  With important contributions by
 *
 *  Maxim Teslenko (maxim.teslenko@nrm.se)
 *
 *  and by many users (run 'acknowledgments' to see more info)
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

#include "bayes.h"
#include "mbbeagle.h"
#include "mcmc.h"
#include "model.h"
#include "utils.h"

const char* const svnRevisionMbbeagleC = "$Rev$";   /* Revision keyword which is expended/updated by svn on each commit/update */

/* Functions and variables defined in mcmc.c that are not exported in mcmc.h */
void    LaunchLogLikeForDivision(int chain, int d, MrBFlt* lnL);

void    FlipCondLikeSpace (ModelInfo *m, int chain, int nodeIndex);
void    FlipNodeScalerSpace (ModelInfo *m, int chain, int nodeIndex);
void    FlipSiteScalerSpace (ModelInfo *m, int chain);
void    FlipTiProbsSpace (ModelInfo *m, int chain, int nodeIndex);
void    ResetSiteScalers (ModelInfo *m, int chain);
void    CopySiteScalers (ModelInfo *m, int chain);

int     TreeCondLikes_Beagle_No_Rescale (Tree *t, int division, int chain);
int     TreeCondLikes_Beagle_Rescale_All (Tree *t, int division, int chain);

extern int *chainId;
extern int numLocalChains;


#if defined (BEAGLE_ENABLED)
/*------------------------------------------------------------------------
|
|   InitBeagleInstance: create and initialize a beagle instance
|
-------------------------------------------------------------------------*/
int InitBeagleInstance (ModelInfo *m, int division)
{
    int                     i, j, k, c, s, *inStates, numPartAmbigTips;
    double                  *inPartials;
    BitsLong                *charBits;
    BeagleInstanceDetails   details;
    long preferedFlags, requiredFlags;
    int resource;
    
    if (m->useBeagle == NO)
        return ERROR;
    
    /* at least one eigen buffer needed */
    if (m->nCijkParts == 0)
        m->nCijkParts = 1;

    /* allocate memory used by beagle */
    m->logLikelihoods          = (MrBFlt *) SafeCalloc ((numLocalChains)*m->numChars, sizeof(MrBFlt));
    m->inRates                 = (MrBFlt *) SafeCalloc (m->numGammaCats, sizeof(MrBFlt));
    m->branchLengths           = (MrBFlt *) SafeCalloc (2*numLocalTaxa, sizeof(MrBFlt));
    m->tiProbIndices           = (int *) SafeCalloc (2*numLocalTaxa, sizeof(int));
    m->inWeights               = (MrBFlt *) SafeCalloc (m->numGammaCats*m->nCijkParts, sizeof(MrBFlt));
    m->bufferIndices           = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->eigenIndices            = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->childBufferIndices      = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->childTiProbIndices      = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->cumulativeScaleIndices  = (int *) SafeCalloc (m->nCijkParts, sizeof(int));

    numPartAmbigTips = 0;
    if (m->numStates != m->numModelStates)
        numPartAmbigTips = numLocalTaxa;
    else
        {
        for (i=0; i<numLocalTaxa; i++)
            {
            if (m->isPartAmbig[i] == YES)
                numPartAmbigTips++;
            }
        }

    if (beagleResourceNumber >= 0 && beagleResourceNumber != 99)
        {
        resource = beagleResourceNumber;
        beagleResourceCount = 1;
        }
    else if (beagleResourceCount != 0) 
        {
        resource = beagleResource[beagleInstanceCount % beagleResourceCount];
        }
    preferedFlags = beagleFlags;
    
    requiredFlags = 0L;
    
    if (beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS)
        requiredFlags |= BEAGLE_FLAG_SCALERS_LOG; //BEAGLE_FLAG_SCALERS_RAW;

    /* TODO: allocate fewer buffers when nCijkParts > 1 */
    /* create beagle instance */
    m->beagleInstance = beagleCreateInstance(numLocalTaxa,
                                             m->numCondLikes * m->nCijkParts,
                                             numLocalTaxa - numPartAmbigTips,
                                             m->numModelStates,
                                             m->numChars,
                                            (numLocalChains + 1) * m->nCijkParts,
                                             m->numTiProbs*m->nCijkParts,
                                             m->numGammaCats,
                                             m->numScalers * m->nCijkParts,
                                             (beagleResourceCount == 0 ? NULL : &resource),
                                             (beagleResourceCount == 0 ? 0 : 1),                                             
                                             preferedFlags,
                                             requiredFlags,
                                             &details);

    if (m->beagleInstance < 0)
        {
        MrBayesPrint ("%s   Failed to start BEAGLE instance\n", spacer);
        return (ERROR);
        }
    else
        {
        MrBayesPrint ("\n%s   Using BEAGLE resource %i for division %d:", spacer, details.resourceNumber, division+1);
#   if defined (THREADS_ENABLED)
        MrBayesPrint (" (%s)\n", (tryToUseThreads ? "threaded" : "non-threaded"));
#   else
        MrBayesPrint (" (non-threaded)\n");
#   endif
        MrBayesPrint ("%s      Rsrc Name : %s\n", spacer, details.resourceName);
        MrBayesPrint ("%s      Impl Name : %s\n", spacer, details.implName);
        MrBayesPrint ("%s      Flags:", spacer);
        BeaglePrintFlags(details.flags);
        MrBayesPrint ("\n");
        beagleInstanceCount++;          
        }

    /* initialize tip data */
    inStates = (int *) SafeMalloc (m->numChars * sizeof(int));
    if (!inStates)
        return ERROR;
    inPartials = (double *) SafeMalloc (m->numChars * m->numModelStates * sizeof(double));
    if (!inPartials)
        return ERROR;
    for (i=0; i<numLocalTaxa; i++)
        {
        if (m->isPartAmbig[i] == NO)
            {
            charBits = m->parsSets[i];
            for (c=0; c<m->numChars; c++)
                {
                for (s=j=0; s<m->numModelStates; s++)
                    {
                    if (IsBitSet(s, charBits))
                        {
                        inStates[c] = s;
                        j++;
                        }
                    }
                if (j == m->numModelStates)
                    inStates[c] = j;
                else
                    assert (j==1);
                charBits += m->nParsIntsPerSite;
                }
            beagleSetTipStates(m->beagleInstance, i, inStates);
            }
        else /* if (m->isPartAmbig == YES) */
            {
            k = 0;
            charBits = m->parsSets[i];
            for (c=0; c<m->numChars; c++)
                {
                for (s=0; s<m->numModelStates; s++)
                    {
                    if (IsBitSet(s%m->numStates, charBits))
                        inPartials[k++] = 1.0;
                    else
                        inPartials[k++] = 0.0;
                    }
                charBits += m->nParsIntsPerSite;
                }
            beagleSetTipPartials(m->beagleInstance, i, inPartials);
            }
        }
    free (inStates);
    free (inPartials);

    return NO_ERROR;
}


/*-----------------------------------------------------------------
|
|   LaunchBEAGLELogLikeForDivision: calculate the log likelihood  
|       of the new state of the chain for a single division
|
-----------------------------------------------------------------*/
void LaunchBEAGLELogLikeForDivision(int chain, int d, ModelInfo* m, Tree* tree, MrBFlt* lnL)
{
    int i, rescaleFreqNew;
    int *isScalerNode;
    TreeNode *p;
    
    if (beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS) 
        {
    
#   if defined (DEBUG_MB_BEAGLE_FLOW)
        MrBayesPrint ("ALWAYS RESCALING\n");
#   endif
        /* Flip and copy or reset site scalers */
        FlipSiteScalerSpace(m, chain);
        if (m->upDateAll == YES) {
            for (i=0; i<m->nCijkParts; i++) {           
                beagleResetScaleFactors(m->beagleInstance, m->siteScalerIndex[chain] + i);
                }
            }
        else
            CopySiteScalers(m, chain);

        TreeTiProbs_Beagle(tree, d, chain);
        TreeCondLikes_Beagle(tree, d, chain);
        TreeLikelihood_Beagle(tree, d, chain, lnL, (chainId[chain] % chainParams.numChains));
        } 
    else 
        { /* MB_BEAGLE_SCALE_DYNAMIC */
    
        /* This flag is only valid within this block */
        m->rescaleBeagleAll = NO;        
        TreeTiProbs_Beagle(tree, d, chain);
        if (m->succesCount[chain] > 1000)
            {
            m->succesCount[chain] = 10;
            m->rescaleFreq[chain]++; /* increase rescaleFreq independent of whether we accept or reject new state */
            m->rescaleFreqOld = rescaleFreqNew = m->rescaleFreq[chain];
            for (i=0; i<tree->nIntNodes; i++)
                {
                p = tree->intDownPass[i];
                if (p->upDateCl == YES) {
                     /* flip to the new workspace since TreeCondLikes_Beagle_Rescale_All() does not do it for
                        (p->upDateCl == YES) since it assumes that TreeCondLikes_Beagle_No_Rescale() did it */
                    FlipCondLikeSpace (m, chain, p->index);
                   }
                }
            goto rescale_all;
            }

        if (beagleScalingFrequency != 0 && 
            m->beagleComputeCount[chain] % (beagleScalingFrequency) == 0)
            {
            m->rescaleFreqOld = rescaleFreqNew = m->rescaleFreq[chain];
            for (i=0; i<tree->nIntNodes; i++)
                {
                p = tree->intDownPass[i];
                if (p->upDateCl == YES) {
                     /* flip to the new workspace since TreeCondLikes_Beagle_Rescale_All() does not do it for
                        (p->upDateCl == YES) since it assumes that TreeCondLikes_Beagle_No_Rescale() did it */
                    FlipCondLikeSpace (m, chain, p->index);
                   }
                }
            goto rescale_all;
            }

        TreeCondLikes_Beagle_No_Rescale(tree, d, chain);

        /* Check if likelihood is valid */      
        if (TreeLikelihood_Beagle(tree, d, chain, lnL, (chainId[chain] % chainParams.numChains)) == BEAGLE_ERROR_FLOATING_POINT) 
            {
            m->rescaleFreqOld = rescaleFreqNew = m->rescaleFreq[chain];
            if (rescaleFreqNew > 1 && m->succesCount[chain] < 40)
                {
                if (m->succesCount[chain] < 10)
                    {
                    if (m->succesCount[chain] < 4)
                        {
                        rescaleFreqNew -= rescaleFreqNew >> 3; /* <== we cut up to 12,5% of rescaleFreq */
                        if (m->succesCount[chain] < 2)
                            {
                            rescaleFreqNew -= rescaleFreqNew >> 3;
                            /* to avoid situation when we may stack at high rescaleFreq when new states do not get accepted because of low liklihood but there proposed frequency is high we reduce rescaleFreq even if we reject the last move*/
                            /* basically the higher probability of proposing of low liklihood state which needs smaller rescaleFreq would lead to higher probability of hitting this code which should reduce rescaleFreqOld thus reduce further probability of hitting this code */
                            /* at some point this negative feedback mechanism should get in balance with the mechanism of periodically increasing rescaleFreq when long sequence of successes is achieved*/
                            m->rescaleFreqOld -= m->rescaleFreqOld >> 3;
                            }
                        m->rescaleFreqOld -= m->rescaleFreqOld >> 3;
                        m->rescaleFreqOld--;
                        m->rescaleFreqOld = (m->rescaleFreqOld ? m->rescaleFreqOld:1);
                        m->recalculateScalers = YES; 
                        recalcScalers = YES;
                        }
                    }
                rescaleFreqNew--;
                rescaleFreqNew = (rescaleFreqNew ? rescaleFreqNew : 1);
                }
            m->succesCount[chain] = 0;
    rescale_all:
#   if defined (DEBUG_MB_BEAGLE_FLOW)
            MrBayesPrint ("NUMERICAL RESCALING\n");
#   endif

            m->rescaleBeagleAll = YES;
            FlipSiteScalerSpace(m, chain);
            isScalerNode = m->isScalerNode[chain];
    while_loop:
            ResetScalersPartition (isScalerNode, tree, rescaleFreqNew);
            for (i=0; i<m->nCijkParts; i++) 
                {           
                beagleResetScaleFactors(m->beagleInstance, m->siteScalerIndex[chain] + i);
                }
            
            TreeCondLikes_Beagle_Rescale_All (tree, d, chain);
            if (TreeLikelihood_Beagle(tree, d, chain, lnL, (chainId[chain] % chainParams.numChains)) == BEAGLE_ERROR_FLOATING_POINT)
                {
                if (rescaleFreqNew > 1)
                    {
                    /* Swap back scalers which were swapped in TreeCondLikes_Beagle_Rescale_All() */
                    for (i=0; i<tree->nIntNodes; i++)
                        {
                        p = tree->intDownPass[i];
                        if (isScalerNode[p->index] == YES)
                            FlipNodeScalerSpace (m, chain, p->index);
                        }
                    rescaleFreqNew -= rescaleFreqNew >> 3; /* <== we cut up to 12,5% of rescaleFreq */
                    rescaleFreqNew--;                      /* we cut extra 1 of rescaleFreq */
                    goto while_loop;
                    }
                }
            m->rescaleFreq[chain] = rescaleFreqNew;
            }
        }
    
    /* Count number of evaluations */
    m->beagleComputeCount[chain]++;
}


void recalculateScalers(int chain)
{
    int         i, d, rescaleFreqNew;
    int         *isScalerNode;
    ModelInfo*  m;
    Tree        *tree;

    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        if (m->recalculateScalers == YES)
            {
            m->recalculateScalers = NO;
            tree = GetTree(m->brlens, chain, state[chain]);

            rescaleFreqNew = m->rescaleFreq[chain];
            isScalerNode = m->isScalerNode[chain];

            ResetScalersPartition (isScalerNode, tree, rescaleFreqNew);
            for (i=0; i<m->nCijkParts; i++) {           
                beagleResetScaleFactors(m->beagleInstance, m->siteScalerIndex[chain] + i);
            }
            /* here it does not matter if we flip CL space or not */
            TreeCondLikes_Beagle_Rescale_All (tree, d, chain);
            }
        }
}


void BeagleAddGPUDevicesToList(int **newResourceList, int *beagleResourceCount)
{
    BeagleResourceList* beagleResources;
    int i, gpuCount;
    
    beagleResources = beagleGetResourceList();
    if (*newResourceList == NULL) {
        *newResourceList = (int*) SafeCalloc(sizeof(int), beagleResources->length);
    }
    gpuCount = 0;
    for (i = 0; i < beagleResources->length; i++) {
        if (beagleResources->list[i].supportFlags & BEAGLE_FLAG_PROCESSOR_GPU) {
            (*newResourceList)[gpuCount] = i;
            gpuCount++;
        }
    }
    *beagleResourceCount = gpuCount;            
}


void BeagleRemoveGPUDevicesFromList(int **beagleResource, int *beagleResourceCount)
{
    *beagleResourceCount = 0;
}


/*-----
|
| BeaglePrintResources: outputs the available BEAGLE resources
|
----------*/
void BeaglePrintResources()
{
    int i;
    BeagleResourceList* beagleResources;
    
    beagleResources = beagleGetResourceList();
    MrBayesPrint ("%s   Available resources reported by beagle library:\n", spacer);
    for (i=0; i<beagleResources->length; i++) 
        {
        MrBayesPrint ("\tResource %i:\n", i);       
        MrBayesPrint ("\tName: %s\n", beagleResources->list[i].name);
        if (i > 0) 
            {
            MrBayesPrint ("\tDesc: %s\n", beagleResources->list[i].description);
            }
        MrBayesPrint ("\tFlags:");
        BeaglePrintFlags(beagleResources->list[i].supportFlags);
        MrBayesPrint ("\n\n");
        }
    MrBayesPrint ("%s   BEAGLE version: %s\n", spacer, beagleGetVersion());
}


int BeagleCheckFlagCompatability(long inFlags)
{
    if (inFlags & BEAGLE_FLAG_PROCESSOR_GPU) {
        if (inFlags & BEAGLE_FLAG_VECTOR_SSE) {
            MrBayesPrint ("%s   Simultaneous use of GPU and SSE not available.\n", spacer);
            return NO;
        }
        if (inFlags & BEAGLE_FLAG_THREADING_OPENMP) {
            MrBayesPrint ("%s   Simultaneous use of GPU and OpenMP not available.\n", spacer);
            return NO;
        }
    }

    return YES;
}


/*-------------------
|
|  BeaglePrintFlags: outputs beagle instance details
|
______________________*/
void BeaglePrintFlags(long inFlags) 
{
    int     i, k;
    char *names[] = { "PROCESSOR_CPU",
                      "PROCESSOR_GPU",
                      "PROCESSOR_FPGA",
                      "PROCESSOR_CELL",
                      "PRECISION_DOUBLE",
                      "PRECISION_SINGLE",
                      "COMPUTATION_ASYNCH",
                      "COMPUTATION_SYNCH",
                      "EIGEN_REAL",
                      "EIGEN_COMPLEX",
                      "SCALING_MANUAL",
                      "SCALING_AUTO",
                      "SCALING_ALWAYS",
                      "SCALING_DYNAMIC",
                      "SCALERS_RAW",
                      "SCALERS_LOG",
                      "VECTOR_NONE",
                      "VECTOR_SSE",
                      "THREADING_NONE",
                      "THREADING_OPENMP"
                    };
    long flags[] = { BEAGLE_FLAG_PROCESSOR_CPU,
                     BEAGLE_FLAG_PROCESSOR_GPU,
                     BEAGLE_FLAG_PROCESSOR_FPGA,
                     BEAGLE_FLAG_PROCESSOR_CELL,
                     BEAGLE_FLAG_PRECISION_DOUBLE,
                     BEAGLE_FLAG_PRECISION_SINGLE,
                     BEAGLE_FLAG_COMPUTATION_ASYNCH,
                     BEAGLE_FLAG_COMPUTATION_SYNCH,
                     BEAGLE_FLAG_EIGEN_REAL,
                     BEAGLE_FLAG_EIGEN_COMPLEX,
                     BEAGLE_FLAG_SCALING_MANUAL,
                     BEAGLE_FLAG_SCALING_AUTO,
                     BEAGLE_FLAG_SCALING_ALWAYS,
                     BEAGLE_FLAG_SCALING_DYNAMIC,
                     BEAGLE_FLAG_SCALERS_RAW,
                     BEAGLE_FLAG_SCALERS_LOG,
                     BEAGLE_FLAG_VECTOR_NONE,
                     BEAGLE_FLAG_VECTOR_SSE,
                     BEAGLE_FLAG_THREADING_NONE,
                     BEAGLE_FLAG_THREADING_OPENMP
                    };

    for (i=k=0; i<20; i++)
        {
        if (inFlags & flags[i])
            {
            if (k%4 == 0 && k > 0)
                MrBayesPrint ("\n%s            ", spacer);
            MrBayesPrint (" %s", names[i]);
            k++;
            }
        }
}


#if defined(THREADS_ENABLED)
void *LaunchThreadLogLikeForDivision(void *arguments)
{
    int d, chain;
    MrBFlt *lnL;
    LaunchStruct* launchStruct;
    
    launchStruct = (LaunchStruct*) arguments;
    chain = launchStruct->chain;
    d = launchStruct->division;
    lnL = launchStruct->lnL;
    LaunchLogLikeForDivision(chain, d, lnL);    
    return 0;
}


MrBFlt LaunchLogLikeForAllDivisionsInParallel(int chain)
{
    int d;
    int threadError;
    pthread_t* threads;
    LaunchStruct* launchValues;
    int* wait;
    ModelInfo* m;
    MrBFlt chainLnLike;
    
    chainLnLike = 0.0;

    /* TODO Initialize only once */
    threads = (pthread_t*) SafeMalloc(sizeof(pthread_t) * numCurrentDivisions);
    launchValues = (LaunchStruct*) SafeMalloc(sizeof(LaunchStruct) * numCurrentDivisions);
    wait = (int*) SafeMalloc(sizeof(int) * numCurrentDivisions);
    
    /* Cycle through divisions and recalculate tis and cond likes as necessary. */
    /* Code below does not try to avoid recalculating ti probs for divisions    */
    /* that could share ti probs with other divisions.                          */
    for (d=0; d<numCurrentDivisions; d++)
        {
        
#   if defined (BEST_MPI_ENABLED)
        if (isDivisionActive[d] == NO)
            continue;
#   endif
        m = &modelSettings[d];
        
        if (m->upDateCl == YES) 
            {                   
            launchValues[d].chain = chain;
            launchValues[d].division = d;
            launchValues[d].lnL = &(m->lnLike[2*chain + state[chain]]);
            /* Fork */
            threadError = pthread_create(&threads[d], NULL, 
                                         LaunchThreadLogLikeForDivision, 
                                         (void*) &launchValues[d]);
            assert (0 == threadError);
            wait[d] = 1;                    
            }           
        else 
            {
            wait[d] = 0;
            }
        }
    
    for (d = 0; d < numCurrentDivisions; d++)
        {
        /* Join */
        if (wait[d]) 
            {
            threadError = pthread_join(threads[d], NULL);
            assert (0 == threadError);
            }               
        m = &modelSettings[d];
        chainLnLike += m->lnLike[2*chain + state[chain]];
        }
            
    /* TODO Free these once */
    free(threads);
    free(launchValues);
    free(wait);
    
    return chainLnLike;
}
#endif


int ScheduleLogLikeForAllDivisions()
{
    int d;
    int divisionsToLaunch = 0;
    ModelInfo       *m;
    
    if (numCurrentDivisions < 2) {
        return 0;
        }
    
    for (d=0; d<numCurrentDivisions; d++) {
        m = &modelSettings[d];
        if (m->upDateCl == YES) {
            divisionsToLaunch++;
            }
        }
    
    return (divisionsToLaunch > 1);
}


/*----------------------------------------------------------------
 |
 |  TreeCondLikes_Beagle: This routine updates all conditional
 |       (partial) likelihoods of a beagle instance while doing no rescaling.
 |      That potentialy can make final liklihood bad then calculation with rescaling needs to be done.
 |
 -----------------------------------------------------------------*/
int TreeCondLikes_Beagle_No_Rescale (Tree *t, int division, int chain)
{
    int                 i, j, cumulativeScaleIndex;
    BeagleOperation     operations;
    TreeNode            *p;
    ModelInfo           *m;
    unsigned            chil1Step, chil2Step;
    int                 *isScalerNode;
    
    m = &modelSettings[division];
    isScalerNode = m->isScalerNode[chain];
    
    for (i=0; i<t->nIntNodes; i++)
    {
        p = t->intDownPass[i];
        
        /* check if conditional likelihoods need updating */
        if (p->upDateCl == YES)
        {
            /* flip to the new workspace */
            FlipCondLikeSpace (m, chain, p->index);
            
            /* update conditional likelihoods */
            operations.destinationPartials    = m->condLikeIndex[chain][p->index       ];
            operations.child1Partials         = m->condLikeIndex[chain][p->left->index ];
            operations.child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ];
            operations.child2Partials         = m->condLikeIndex[chain][p->right->index];
            operations.child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index];
            
            /* All partials for tips are the same across omega categories, thus we are doing the following two if statments.*/
            if (p->left->left== NULL)
                chil1Step=0;
            else
                chil1Step=1;
            
            if (p->right->left== NULL)
                chil2Step=0;
            else
                chil2Step=1;
            
            operations.destinationScaleWrite = BEAGLE_OP_NONE;
            cumulativeScaleIndex  = BEAGLE_OP_NONE;
            if (isScalerNode[p->index] == YES)
                {
                operations.destinationScaleRead  = m->nodeScalerIndex[chain][p->index];
                }
            else
                {
                operations.destinationScaleRead  = BEAGLE_OP_NONE;
                }
            
            for (j=0; j<m->nCijkParts; j++)
            {
                beagleUpdatePartials(m->beagleInstance,
                                     &operations,
                                     1,
                                     cumulativeScaleIndex);
                
                operations.destinationPartials++;
                operations.child1Partials+=chil1Step;
                operations.child1TransitionMatrix++;
                operations.child2Partials+=chil2Step;
                operations.child2TransitionMatrix++;
                
                if (isScalerNode[p->index] == YES)
                    operations.destinationScaleRead++;
            }
        }
    }
    
    return NO_ERROR;
}


/*----------------------------------------------------------------
 |
 |  TreeCondLikes_Beagle: This routine updates all conditional
 |       (partial) likelihoods of a beagle instance while rescaling at every node.
 |        Note: all nodes get recalculated, not only tached by move.
 |
 -----------------------------------------------------------------*/
int TreeCondLikes_Beagle_Rescale_All (Tree *t, int division, int chain)
{
    int                 i, j, cumulativeScaleIndex;
    BeagleOperation     operations;
    TreeNode            *p;
    ModelInfo           *m;
    unsigned            chil1Step, chil2Step;
    int                 *isScalerNode;
    
    m = &modelSettings[division];
    isScalerNode = m->isScalerNode[chain];
    
    for (i=0; i<t->nIntNodes; i++)
    {
        p = t->intDownPass[i];
        
        if (p->upDateCl == NO) {
            //p->upDateCl = YES;
            /* flip to the new workspace */
            FlipCondLikeSpace (m, chain, p->index);
        }
        
        /* update conditional likelihoods */
        operations.destinationPartials    = m->condLikeIndex[chain][p->index       ];
        operations.child1Partials         = m->condLikeIndex[chain][p->left->index ];
        operations.child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ];
        operations.child2Partials         = m->condLikeIndex[chain][p->right->index];
        operations.child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index];
        
        /* All partials for tips are the same across omega catigoris, thus we are doing the following two if statments.*/
        if (p->left->left== NULL)
            chil1Step=0;
        else
            chil1Step=1;
        
        if (p->right->left== NULL)
            chil2Step=0;
        else
            chil2Step=1;

        operations.destinationScaleRead = BEAGLE_OP_NONE;
        if (isScalerNode[p->index] == YES)
            {
            FlipNodeScalerSpace (m, chain, p->index);
            operations.destinationScaleWrite = m->nodeScalerIndex[chain][p->index];
            cumulativeScaleIndex  = m->siteScalerIndex[chain];
            }
        else
            {
            operations.destinationScaleWrite = BEAGLE_OP_NONE;
            cumulativeScaleIndex  = BEAGLE_OP_NONE;
            }
        
        
        
        for (j=0; j<m->nCijkParts; j++)
        {
            beagleUpdatePartials(m->beagleInstance,
                                 &operations,
                                 1,
                                 cumulativeScaleIndex);
            
            operations.destinationPartials++;
            operations.child1Partials+=chil1Step;
            operations.child1TransitionMatrix++;
            operations.child2Partials+=chil2Step;
            operations.child2TransitionMatrix++;
            
            if (isScalerNode[p->index] == YES) {
                operations.destinationScaleWrite++;
                cumulativeScaleIndex++;
                }

            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   TreeCondLikes_Beagle: This routine updates all conditional
|       (partial) likelihoods of a beagle instance.
|
-----------------------------------------------------------------*/
int TreeCondLikes_Beagle (Tree *t, int division, int chain)
{
    int                 i, j, destinationScaleRead, cumulativeScaleIndex;
    BeagleOperation     operations;
    TreeNode            *p;
    ModelInfo           *m;
    unsigned            chil1Step, chil2Step;
    int                 *isScalerNode;
    
    m = &modelSettings[division];
    isScalerNode = m->isScalerNode[chain];
    
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        
        /* check if conditional likelihoods need updating */
        if (p->upDateCl == YES)
            {
            /* remove old scalers */
            if (m->unscaledNodes[chain][p->index] == 0 && m->upDateAll == NO)
                {
                destinationScaleRead = m->nodeScalerIndex[chain][p->index];
                cumulativeScaleIndex = m->siteScalerIndex[chain];
                for (j=0; j<m->nCijkParts; j++)
                    {
                    beagleRemoveScaleFactors(m->beagleInstance,
                                             &destinationScaleRead,
                                             1,
                                             cumulativeScaleIndex);
                    destinationScaleRead++;
                    cumulativeScaleIndex++;
                    }
                }

            /* flip to the new workspace */
            FlipCondLikeSpace (m, chain, p->index);
            FlipNodeScalerSpace (m, chain, p->index);
            
            /* update conditional likelihoods */
            operations.destinationPartials    = m->condLikeIndex[chain][p->index       ];
            operations.child1Partials         = m->condLikeIndex[chain][p->left->index ];
            operations.child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ];
            operations.child2Partials         = m->condLikeIndex[chain][p->right->index];
            operations.child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index];

            /* All partials for tips are the same across omega categories, thus we are doing the following two if statments.*/
            if (p->left->left== NULL && p->left->right== NULL)
                chil1Step=0;
            else
                chil1Step=1;

            if (p->right->left== NULL && p->right->right== NULL)
                chil2Step=0;
            else
                chil2Step=1;

            /* Now deal with scalers */
            m->unscaledNodes[chain][p->index] = 1 + m->unscaledNodes[chain][p->left->index] + m->unscaledNodes[chain][p->right->index];

            if (m->unscaledNodes[chain][p->index] >= m->rescaleFreq[chain])
                {
                m->unscaledNodes[chain][p->index] = 0;
                operations.destinationScaleWrite = m->nodeScalerIndex[chain][p->index];
                cumulativeScaleIndex  = m->siteScalerIndex[chain];
                }
            else
                {
                operations.destinationScaleWrite = BEAGLE_OP_NONE;
                cumulativeScaleIndex  = BEAGLE_OP_NONE;
                }
            operations.destinationScaleRead = BEAGLE_OP_NONE;

            for (j=0; j<m->nCijkParts; j++)
                {
                beagleUpdatePartials(m->beagleInstance,
                                     &operations,
                                     1,
                                     cumulativeScaleIndex);

                operations.destinationPartials++;
                operations.child1Partials+=chil1Step;
                operations.child1TransitionMatrix++;
                operations.child2Partials+=chil2Step;
                operations.child2TransitionMatrix++;
                if (isScalerNode[p->index] == YES)
                    {
                    operations.destinationScaleWrite++;
                    cumulativeScaleIndex++;
                    }
                }

            }
        } /* end of for */

    return NO_ERROR;
}


/**---------------------------------------------------------------------------
|
|   TreeLikelihood_Beagle: Accumulate the log likelihoods calculated by Beagle
|      at the root.
|
---------------------------------------- -------------------------------------*/
int TreeLikelihood_Beagle (Tree *t, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int         i, j, c = 0, nStates, hasPInvar, beagleReturn;
    MrBFlt      *swr, s01, s10, probOn, probOff, covBF[40], pInvar=0.0, *bs, freq, likeI, lnLikeI, diff, *omegaCatFreq;
    CLFlt       *clInvar=NULL, *nSitesOfPat;
    double      *nSitesOfPat_Beagle;
    TreeNode    *p;
    ModelInfo   *m;
    double      pUnobserved;

#   if defined (MB_PRINT_DYNAMIC_RESCALE_FAIL_STAT)
    static unsigned countBeagleDynamicFail=0;
    static unsigned countALL=0;
#   endif

    /* find root node */
    p = t->root->left;

    /* find model settings and nStates, pInvar, invar cond likes */
    m = &modelSettings[division];
    
    nStates = m->numModelStates;
    if (m->pInvar == NULL)
        {
        hasPInvar = NO;
        }
    else
        {
        hasPInvar = YES;
        pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
        clInvar = m->invCondLikes;
        }
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* if covarion model, adjust base frequencies */
    if (m->switchRates != NULL)
        {
        /* find the stationary frequencies */
        swr = GetParamVals(m->switchRates, chain, state[chain]);
        s01 = swr[0];
        s10 = swr[1];
        probOn = s01 / (s01 + s10);
        probOff =  1.0 - probOn;

        /* now adjust the base frequencies; on-state stored first in cond likes */
        for (j=0; j<nStates/2; j++)
            {
            covBF[j] = bs[j] * probOn;
            covBF[j+nStates/2] = bs[j] * probOff;
            }

        /* finally set bs pointer to adjusted values */
        bs = covBF;
        }

    /* TODO Really only need to check if state frequencies have changed */
    if (m->upDateCijk == YES)
        {
        /* set base frequencies in BEAGLE instance */
        for (i=0; i<m->nCijkParts; i++)
            beagleSetStateFrequencies(m->beagleInstance,
                                      m->cijkIndex[chain] + i,
                                      bs);
        }

    /* find category frequencies */
    if (hasPInvar == NO)
        freq = 1.0 / m->numGammaCats;
    else
        freq = (1.0 - pInvar) / m->numGammaCats;

    /* TODO: cat weights only need to be set when they change */
    /* set category frequencies in beagle instance */
    if (m->numOmegaCats > 1)
        {
        omegaCatFreq = GetParamSubVals(m->omega, chain, state[chain]);
        for (i=0; i<m->nCijkParts; i++)
            {
            for (j=0; j<m->numGammaCats; j++)
                m->inWeights[j] = freq * omegaCatFreq[i];
            beagleSetCategoryWeights(m->beagleInstance,
                                     m->cijkIndex[chain] + i,
                                     m->inWeights);
            }
        }
    else if (hasPInvar == YES)
        {
        for (i=0; i<m->numGammaCats; i++)
            m->inWeights[i] = freq;
        beagleSetCategoryWeights(m->beagleInstance,
                                 m->cijkIndex[chain],
                                 m->inWeights);
        }

    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* TODO: pattern weights only need to be set when they change */
    /* set pattern weights in beagle instance if using dynamic reweighting */
    if (chainParams.weightScheme[0] + chainParams.weightScheme[1] > ETA)
        {
        nSitesOfPat_Beagle = (double *) SafeMalloc (m->numChars * sizeof(double));
        for (c=0; c<m->numChars; c++)
            nSitesOfPat_Beagle[c] = numSitesOfPat[m->compCharStart + c];
        beagleSetPatternWeights(m->beagleInstance,
                                nSitesOfPat_Beagle);
        SafeFree ((void **)(&nSitesOfPat_Beagle));
        }

    /* find root log likelihoods and scalers */
    for (i=0; i<m->nCijkParts; i++)
        {
        m->bufferIndices[i] = m->condLikeIndex[chain][p->index] + i;
        m->eigenIndices[i]  = m->cijkIndex[chain] + i;
        m->cumulativeScaleIndices[i] = m->siteScalerIndex[chain] + i;
        if (t->isRooted == NO)
            {
            m->childBufferIndices[i]     = m->condLikeIndex  [chain][p->anc->index];
            m->childTiProbIndices[i]     = m->tiProbsIndex   [chain][p->index] + i;
            }
        }

    /* reset lnL */
    *lnL = 0.0;

    /* get root log likelihoods */
    if (t->isRooted == YES)
        {        
        beagleReturn = beagleCalculateRootLogLikelihoods(m->beagleInstance,
                                          m->bufferIndices,
                                          m->eigenIndices,
                                          m->eigenIndices,
                                          m->cumulativeScaleIndices,
                                          m->nCijkParts,
                                          lnL);

        }
    else
        {
        beagleReturn = beagleCalculateEdgeLogLikelihoods(m->beagleInstance,
                                          m->bufferIndices,
                                          m->childBufferIndices,
                                          m->childTiProbIndices,
                                          NULL,
                                          NULL,
                                          m->eigenIndices,
                                          m->eigenIndices,
                                          m->cumulativeScaleIndices,
                                          m->nCijkParts,
                                          lnL,
                                          NULL,
                                          NULL);       
        }
#   if defined (MB_PRINT_DYNAMIC_RESCALE_FAIL_STAT)
    countALL++;
#   endif
    if (beagleReturn == BEAGLE_ERROR_FLOATING_POINT)
        {
#   if defined (MB_PRINT_DYNAMIC_RESCALE_FAIL_STAT)
        countBeagleDynamicFail++;
        MrBayesPrint ("DEBUG INFO (not an error) countBeagleDynamicFail:%d countALL:%d\n", countBeagleDynamicFail, countALL);
#   endif
        return beagleReturn;
        }
    assert (beagleReturn == BEAGLE_SUCCESS);
    m->succesCount[chain]++;
    
    /* accumulate logs across sites */
    if (hasPInvar == NO)
        {
        if (m->dataType == RESTRICTION)
            {
            beagleGetSiteLogLikelihoods(m->beagleInstance, m->logLikelihoods);
            (*lnL) = 0.0;
            pUnobserved = 0.0;
            for (c=0; c<m->numDummyChars; c++)
                {
                pUnobserved +=  exp((double)m->logLikelihoods[c]);
                }
            /* correct for absent characters */
            (*lnL) -= log (1-pUnobserved) * (m->numUncompressedChars);
            for (; c<m->numChars; c++)
                {
                (*lnL) += m->logLikelihoods[c] * nSitesOfPat[c];
                }
            }
        /* already done, just check for numerical errors */
        assert ((*lnL) == (*lnL));
        }
    else
        {
        /* has invariable category */
        beagleGetSiteLogLikelihoods(m->beagleInstance,
                                    m->logLikelihoods);
        (*lnL) = 0.0;
        for (c=0; c<m->numChars; c++)
            {
            likeI = 0.0;
            for (j=0; j<nStates; j++)
                likeI += (*(clInvar++)) * bs[j];
            if (likeI != 0.0)
                {
                lnLikeI = log(likeI * pInvar);
                diff = lnLikeI - m->logLikelihoods[c];
                }
            else
                diff = -1000.0;
            if (diff < -200.0)
                (*lnL) += m->logLikelihoods[c] * nSitesOfPat[c];
            else if (diff > 200.0)
                (*lnL) += lnLikeI * nSitesOfPat[c];
            else
                {
                (*lnL) += (m->logLikelihoods[c] + log(1.0 + exp(diff))) * nSitesOfPat[c];
                }

            /* check for numerical errors */
            assert ((*lnL) == (*lnL));
            }       
        }
        
    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   TreeTiProbs_Beagle: This routine updates all transition
|       probability matrices of a beagle instance.
|
-----------------------------------------------------------------*/
int TreeTiProbs_Beagle (Tree *t, int division, int chain)
{
    int         i, j, k, count;
    MrBFlt      correctionFactor, theRate, baseRate, *catRate, length;
    TreeNode    *p;
    ModelInfo   *m;
    
    /* get model settings */
    m = &modelSettings[division];
    
    /* find the correction factor to make branch lengths
       in terms of expected number of substitutions per character */
    correctionFactor = 1.0;
    if (m->dataType == DNA || m->dataType == RNA)
        {
        if (m->nucModelId == NUCMODEL_DOUBLET)
            correctionFactor = 2.0;
        else if (m->nucModelId == NUCMODEL_CODON)
            correctionFactor = 3.0;
        }

    /* get rate multipliers (for gamma & partition specific rates) */
    theRate = 1.0;
    baseRate = GetRate (division, chain);
    
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
        
    /* get category rates for gamma */
    if (m->shape == NULL)
        catRate = &theRate;
    else
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    
    /* get effective category rates */
    for (k=0; k<m->numGammaCats; k++)
        m->inRates[k] = baseRate * catRate[k] * correctionFactor;

    /* TODO: only need to set category rates when they change */
    /* set category rates */
    beagleSetCategoryRates(m->beagleInstance, m->inRates);
    
    /* get ti prob indices and branch lengths to update */
    for (i=count=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        
        /* check if transition probs need updating */
        if (p->upDateTi == YES)
            {
            /* flip transition probability */
            FlipTiProbsSpace (m, chain, p->index);
            
            /* find length */
            if (m->cppEvents != NULL)
                {
                length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
                }
            else if (m->tk02BranchRates != NULL)
                {
                length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
                }
            else if (m->igrBranchRates != NULL)
                {
                length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
                }
            else if (m->mixedBrchRates != NULL)
                {
                length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
                }
            else
                length = p->length;

            /* numerical errors might ensue if we allow very large or very small branch lengths, which might
               occur in relaxed clock models; an elegant solution would be to substitute the stationary
               probs and initial probs but for now we truncate lengths at small or large values */
            if (length > BRLENS_MAX)
                length = BRLENS_MAX;
            else if (length < BRLENS_MIN)
                length = BRLENS_MIN;

            m->branchLengths[count] = length;
            
            /* find index */
            m->tiProbIndices[count] = m->tiProbsIndex[chain][p->index];
            count++;
            }
        }

    /* TODO: only need to update branches that have changed */
    /* calculate transition probabilities */
    if (count > 0) {
        for (i=0; i<m->nCijkParts; i++)
            {
            beagleUpdateTransitionMatrices(m->beagleInstance,
                                           m->cijkIndex[chain] + i,
                                           m->tiProbIndices,
                                           NULL,
                                           NULL,
                                           m->branchLengths,
                                           count);
            for (j=0; j<count; j++)
                m->tiProbIndices[j]++;
            }
        }

    /* return success */
    return NO_ERROR;
}

#endif // BEAGLE_ENABLED


void BeagleNotLinked()
{
    MrBayesPrint ("%s   BEAGLE library is not linked to this executable.\n", spacer);
}


void BeagleThreadsNotLinked()
{
    MrBayesPrint ("%s   Pthreads library is not linked to this executable.\n", spacer);
}

