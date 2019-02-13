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

// #define DEBUG_MB_BEAGLE_MULTIPART
// #define DEBUG_MB_BEAGLE_MULTIPART_SITELNL

const char* const svnRevisionMbbeagleC = "$Rev$";   /* Revision keyword which is expended/updated by svn on each commit/update */

/* Functions and variables defined in mcmc.c that are not exported in mcmc.h */
void    LaunchLogLikeForDivision(int chain, int d, MrBFlt* lnL);

void    FlipCondLikeSpace (ModelInfo *m, int chain, int nodeIndex);
void    FlipNodeScalerSpace (ModelInfo *m, int chain, int nodeIndex);
void    FlipSiteScalerSpace (ModelInfo *m, int chain);
void    FlipTiProbsSpace (ModelInfo *m, int chain, int nodeIndex);
void    ResetSiteScalers (ModelInfo *m, int chain);
void    CopySiteScalers (ModelInfo *m, int chain);


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
    
    if (m->beagleInstance == -99) {
        MrBayesPrint ("\n%s   Error: Consecutive MCMC runs not currently supported with BEAGLE. Please restart MrBayes.\n\n", spacer);
        exit(1);
    }


    if (m->useBeagle == NO)
        return ERROR;
    
    /* at least one eigen buffer needed */
    if (m->nCijkParts == 0)
        m->nCijkParts = 1;

    /* allocate memory used by beagle */
    m->logLikelihoods          = (MrBFlt *) SafeCalloc ((numLocalChains)*m->numChars, sizeof(MrBFlt));
    m->inRates                 = (MrBFlt *) SafeCalloc (m->numRateCats, sizeof(MrBFlt));
    m->branchLengths           = (MrBFlt *) SafeCalloc (2*numLocalTaxa, sizeof(MrBFlt));
    m->tiProbIndices           = (int *) SafeCalloc (2*numLocalTaxa, sizeof(int));
    m->inWeights               = (MrBFlt *) SafeCalloc (m->numRateCats*m->nCijkParts, sizeof(MrBFlt));
    m->bufferIndices           = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->eigenIndices            = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->childBufferIndices      = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->childTiProbIndices      = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->cumulativeScaleIndices  = (int *) SafeCalloc (m->nCijkParts, sizeof(int));
    m->operations              = (BeagleOperation *) SafeCalloc(m->numCondLikes*m->nCijkParts, sizeof(BeagleOperation));
    m->scaleFactorsOps         = (int *) SafeCalloc (2*numLocalTaxa*m->nCijkParts, sizeof(int));

    numPartAmbigTips = 0;
    if (m->numStates != m->numModelStates)
        numPartAmbigTips = numLocalTaxa;
    else
        {
        for (i=0; i<numLocalTaxa; i++)
            {
#if defined(BEAGLE_V3_ENABLED)
            if (m->isPartAmbig[i] == YES || beagleAllFloatTips)
#else 
            if (m->isPartAmbig[i] == YES)
#endif
                {
                numPartAmbigTips++;
                }
            }
        }

    m->beagleInstance = createBeagleInstance(m, m->nCijkParts, m->numRateCats, m->numModelStates, m->numCondLikes, m->numScalers, m->numChars, m->numTiProbs, numPartAmbigTips, division);

    if (m->beagleInstance < 0)
        return ERROR;

    /* initialize tip data */
    inStates = (int *) SafeMalloc (m->numChars * sizeof(int));
    if (!inStates)
        return ERROR;
    inPartials = (double *) SafeMalloc (m->numChars * m->numModelStates * sizeof(double));
    if (!inPartials)
        return ERROR;
    for (i=0; i<numLocalTaxa; i++)
        {
#if defined(BEAGLE_V3_ENABLED)
        if (m->isPartAmbig[i] == NO && !(beagleAllFloatTips))
#else
        if (m->isPartAmbig[i] == NO)
#endif
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


int createBeagleInstance(ModelInfo *m, int nCijkParts, int numRateCats, int numModelStates, int numCondLikes, int numScalers, int numChars, int numTiProbs, int numPartAmbigTips, int division)
{
    int                     resource, resourceCount, beagleInstance;
    BeagleInstanceDetails   details;
    long                    preferredFlags, requiredFlags;

#if defined (BEAGLE_V3_ENABLED)
    int                     i;
    Tree*                   t;
#endif

    resourceCount = beagleResourceCount;

    preferredFlags = beagleFlags;
    requiredFlags = 0L;

    if (beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS)
        {
        requiredFlags |= BEAGLE_FLAG_SCALERS_LOG;
        }

    if (beagleResourceNumber >= 0 && beagleResourceNumber != 99)
        {
        resource = beagleResourceNumber;
        resourceCount = 1;
        }
    else if (beagleResourceCount != 0) 
        {
#   if defined (MPI_ENABLED)
      resource = beagleResource[(beagleInstanceCount + proc_id) % beagleResourceCount];
#   else
      resource = beagleResource[beagleInstanceCount % beagleResourceCount];
#   endif
        }
#if defined (BEAGLE_V3_ENABLED)
    else if (beagleResourceNumber == 99)
        {
        int numInstancePartitions = 1;
        if (division == -1)
            {
            numInstancePartitions = numCurrentDivisions;
            }

        MrBayesPrint ("\n%s   Running benchmarks to automatically select fastest BEAGLE resource... ", spacer);

        long benchmarkFlags = BEAGLE_BENCHFLAG_SCALING_NONE;
        if (beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS) 
            {
            benchmarkFlags = BEAGLE_BENCHFLAG_SCALING_ALWAYS;
            }
        else if (beagleScalingScheme == MB_BEAGLE_SCALE_DYNAMIC)
            {
            benchmarkFlags = BEAGLE_BENCHFLAG_SCALING_DYNAMIC;
            }

        // select fastest resource
        BeagleBenchmarkedResourceList* rBList;
        rBList = beagleGetBenchmarkedResourceList(
                    numLocalTaxa,
                    numLocalTaxa - numPartAmbigTips,
                    numModelStates,
                    numChars,
                    numRateCats,
                    NULL,
                    0,
                    preferredFlags,
                    requiredFlags,
                    nCijkParts,
                    numInstancePartitions,
                    0,
                    benchmarkFlags);

#       if defined (DEBUG_MB_BEAGLE_FLOW)
            fprintf(stdout, "Resource benchmarks:\n");
            fprintf(stdout, "\ttipCount : %d\n", numLocalTaxa);
            fprintf(stdout, "\tcompactBufferCount : %d\n", numLocalTaxa - numPartAmbigTips);
            fprintf(stdout, "\tstateCount : %d\n", numModelStates);
            fprintf(stdout, "\tpatternCount : %d\n", numChars);
            fprintf(stdout, "\tcategoryCount : %d\n", numRateCats);
            fprintf(stdout, "\teigenModelCount : %d\n", nCijkParts);
            fprintf(stdout, "\tpartitionCount : %d\n", numInstancePartitions);
            fprintf(stdout, "\tPreferred Flags:");
            BeaglePrintFlags(preferredFlags);
            fprintf(stdout, "\n");
            fprintf(stdout, "\tRequired Flags:");
            BeaglePrintFlags(requiredFlags);
            fprintf(stdout, "\n");
            fprintf(stdout, "\tBenchmark Flags: ");
            (benchmarkFlags & BEAGLE_BENCHFLAG_SCALING_ALWAYS ? fprintf(stdout, "BEAGLE_BENCHFLAG_SCALING_ALWAYS") : fprintf(stdout, "BEAGLE_BENCHFLAG_SCALING_DYNAMIC"));
            fprintf(stdout, "\n");
            fprintf(stdout, "\n");



            if (rBList != NULL)
                {
                for (int i = 0; i < rBList->length; i++)
                    {
                    fprintf(stdout, "\tResource %i:\n\t\tName : %s\n", i, rBList->list[i].name);
                    fprintf(stdout, "\t\tDesc : %s\n", rBList->list[i].description);
                    fprintf(stdout, "\t\tSupport Flags:");
                    BeaglePrintFlags(rBList->list[i].supportFlags);
                    fprintf(stdout, "\n");
                    fprintf(stdout, "\t\tRequired Flags:");
                    BeaglePrintFlags(rBList->list[i].requiredFlags);
                    fprintf(stdout, "\n");
                    fprintf(stdout, "\t\tBenchmark Results:\n");
                    fprintf(stdout, "\t\t\tNmbr : %d\n", rBList->list[i].number);
                    fprintf(stdout, "\t\t\tImpl : %s\n", rBList->list[i].implName);
                    fprintf(stdout, "\t\t\tFlags:");
                    BeaglePrintFlags(rBList->list[i].benchedFlags);
                    fprintf(stdout, "\n");
                    fprintf(stdout, "\t\t\tPerf : %.4f ms (%.2fx CPU)\n", rBList->list[i].benchmarkResult, rBList->list[i].performanceRatio);
                   }
                }
            fprintf(stdout, "\n");
#       endif

        if (rBList != NULL)
            {
            double fastestTime = rBList->list[0].benchmarkResult;
            resource = rBList->list[0].number;
            resourceCount = 1;

            for (int i = 1; i < rBList->length; i++)
                {
                if (rBList->list[i].benchmarkResult < fastestTime)
                    {
                    fastestTime = rBList->list[i].benchmarkResult;
                    resource = rBList->list[i].number;
                    }
                }
            }
        }
#endif /* BEAGLE_V3_ENABLED */

    /* TODO: allocate fewer buffers when nCijkParts > 1 */
    /* create beagle instance */
    beagleInstance = beagleCreateInstance(numLocalTaxa,
                                          numCondLikes * nCijkParts,
                                          numLocalTaxa - numPartAmbigTips,
                                          numModelStates,
                                          numChars,
                                          ((numLocalChains + 1) * nCijkParts) * numCurrentDivisions,
                                          numTiProbs*nCijkParts,
                                          numRateCats,
                                          numScalers * nCijkParts,
                                          (resourceCount == 0 ? NULL : &resource),
                                          (resourceCount == 0 ? 0 : 1),
                                          preferredFlags,
                                          requiredFlags,
                                          &details);

    if (beagleInstance < 0)
        {
        MrBayesPrint ("%s   Failed to start BEAGLE instance\n", spacer);
        }
    else 
        {
        if (division < 0)
            {
            /* do not use multipartition mode (division < 0) for CPU resource */
            if (details.flags & BEAGLE_FLAG_FRAMEWORK_CPU)
                {
                MrBayesPrint ("\n%s   Selected resource is the CPU, changing to multi-instance BEAGLE mode\n", spacer);
                beagleFinalizeInstance(beagleInstance);
                return -1;
                }
            else
                {
                MrBayesPrint ("\n%s   Using BEAGLE v%s resource %i for %d division%s:\n", spacer, beagleGetVersion(), details.resourceNumber, numCurrentDivisions, (numCurrentDivisions > 1 ? "s" : ""));
                }
            }
        else
            {
            MrBayesPrint ("\n%s   Using BEAGLE v%s resource %i for division %d:\n", spacer, beagleGetVersion(), details.resourceNumber, division+1);
            }
        MrBayesPrint ("%s      Rsrc Name : %s\n", spacer, details.resourceName);
        MrBayesPrint ("%s      Impl Name : %s\n", spacer, details.implName);
        MrBayesPrint ("%s      Flags:", spacer);
        BeaglePrintFlags(details.flags);
MrBayesPrint ("%s      MODEL STATES: %d", spacer, numModelStates);
        MrBayesPrint ("\n");
        beagleInstanceCount++;
        }

#if defined (BEAGLE_V3_ENABLED)
    /* use level-order traversal with CUDA implementation or OpenCL with multi-partition */
    if(((details.flags & BEAGLE_FLAG_FRAMEWORK_CUDA) && division < 1 ) ||
        ((details.flags & BEAGLE_FLAG_FRAMEWORK_OPENCL) && division < 0))
        {
        for (i=0; i<(numTrees * 2 * numGlobalChains); i++)
            {
            t = mcmcTree[i];
            if (t != NULL)
                {
                if ((t->intDownPassLevel = (TreeNode **) SafeCalloc (numTaxa, sizeof (TreeNode *))) == NULL)
                    {
                    free (t->nodes);
                    free (t->allDownPass);
                    free (t);
                    return -1;
                    }
                t->levelPassEnabled = 1;
                GetDownPass(t);
                }
            }
        }
    /* set max number of CPU threads to be used by BEAGLE with CPU implementation */
    if ((details.flags & BEAGLE_FLAG_FRAMEWORK_CPU) && beagleThreadCount > 1 && beagleThreadCount != 99)
        {
        beagleSetCPUThreadCount(beagleInstance, beagleThreadCount);
        }
#endif /* BEAGLE_V3_ENABLED */

    return beagleInstance;
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
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
        printf("m->upDateAll %d\n", m->upDateAll);
#endif
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
        TreeCondLikes_Beagle_Always_Rescale(tree, d, chain);
        TreeLikelihood_Beagle(tree, d, chain, lnL, (chainId[chain] % chainParams.numChains));
        } 
    else 
        { /* MB_BEAGLE_SCALE_DYNAMIC */
    
        /* This flag is only valid within this block */
        m->rescaleBeagleAll = NO;        
        TreeTiProbs_Beagle(tree, d, chain);
        if (m->successCount[chain] > 1000)
            {
            m->successCount[chain] = 10;
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
            if (rescaleFreqNew > 1 && m->successCount[chain] < 40)
                {
                if (m->successCount[chain] < 10)
                    {
                    if (m->successCount[chain] < 4)
                        {
                        rescaleFreqNew -= rescaleFreqNew >> 3; /* <== we cut up to 12,5% of rescaleFreq */
                        if (m->successCount[chain] < 2)
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
            m->successCount[chain] = 0;
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
#if defined (BEAGLE_V3_ENABLED)
    int         divisionCount;
    int         *divisions;
#endif

    if (modelSettings[0].useBeagleMultiPartitions == NO)
        {
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
#if defined (BEAGLE_V3_ENABLED)
    else 
        {
        divisions = (int *) SafeCalloc (numCurrentDivisions, sizeof(int));
        divisionCount = 0;

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
                    beagleResetScaleFactorsByPartition(m->beagleInstance, m->siteScalerIndex[chain] + i, m->divisionIndex);
                }

                divisions[divisionCount++] = d;
                }
            }

        /* here it does not matter if we flip CL space or not */
        TreeCondLikes_BeagleMultiPartition_Rescale_All(divisions, divisionCount, chain);

        free(divisions);
        }
#endif
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
#if defined (BEAGLE_V3_ENABLED)
                      "THREADING_CPP",
#endif
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
#if defined (BEAGLE_V3_ENABLED)
                     BEAGLE_FLAG_THREADING_CPP,
#endif
                     BEAGLE_FLAG_THREADING_OPENMP
                    };

    int flagCount = 20;

#if defined (BEAGLE_V3_ENABLED)
    flagCount = 21;
#endif

    for (i=k=0; i<flagCount; i++)
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
    int                 i, j, cumulativeScaleIndex, op;
    TreeNode            *p;
    ModelInfo           *m;
    unsigned            chil1Step, chil2Step;
    int                 *isScalerNode;
    
    m = &modelSettings[division];
    isScalerNode = m->isScalerNode[chain];
    
    op = 0;

    for (i=0; i<t->nIntNodes; i++)
        {

#if defined (BEAGLE_V3_ENABLED)
    if (t->levelPassEnabled)
        {
        p = t->intDownPassLevel[i];
        }
    else
#endif
        {
        p = t->intDownPass[i];
        }
        
        /* check if conditional likelihoods need updating */
        if (p->upDateCl == YES)
            {
            /* flip to the new workspace */
            FlipCondLikeSpace (m, chain, p->index);
            
            /* update conditional likelihoods */
            m->operations[op].destinationPartials    = m->condLikeIndex[chain][p->index       ];
            m->operations[op].child1Partials         = m->condLikeIndex[chain][p->left->index ];
            m->operations[op].child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ];
            m->operations[op].child2Partials         = m->condLikeIndex[chain][p->right->index];
            m->operations[op].child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index];
            
            /* All partials for tips are the same across omega categories, thus we are doing the following two if statments.*/
            if (p->left->left== NULL)
                chil1Step=0;
            else
                chil1Step=1;
            
            if (p->right->left== NULL)
                chil2Step=0;
            else
                chil2Step=1;
            
            m->operations[op].destinationScaleWrite = BEAGLE_OP_NONE;
            cumulativeScaleIndex  = BEAGLE_OP_NONE;
            if (isScalerNode[p->index] == YES)
                {
                m->operations[op].destinationScaleRead  = m->nodeScalerIndex[chain][p->index];
                }
            else
                {
                m->operations[op].destinationScaleRead  = BEAGLE_OP_NONE;
                }
            
            for (j=1; j<m->nCijkParts; j++)
                {   
                m->operations[op+1].destinationPartials    = m->operations[op].destinationPartials + 1;
                m->operations[op+1].child1Partials         = m->operations[op].child1Partials + chil1Step;
                m->operations[op+1].child1TransitionMatrix = m->operations[op].child1TransitionMatrix + 1;
                m->operations[op+1].child2Partials         = m->operations[op].child2Partials + chil2Step;
                m->operations[op+1].child2TransitionMatrix = m->operations[op].child2TransitionMatrix + 1;
                
                m->operations[op+1].destinationScaleWrite = BEAGLE_OP_NONE;
                if ( isScalerNode[p->index] == YES )
                    {
                    m->operations[op+1].destinationScaleRead = m->operations[op].destinationScaleRead + 1;
                    }
                else 
                    {
                    m->operations[op+1].destinationScaleRead = BEAGLE_OP_NONE;
                    }
                op++;
                }
            op++;
            }
        }

    beagleUpdatePartials(m->beagleInstance,
                         m->operations,
                         op,
                         cumulativeScaleIndex);
    
    return NO_ERROR;
}


/*----------------------------------------------------------------
 |
 |  TreeCondLikes_Beagle: This routine updates all conditional
 |       (partial) likelihoods of a beagle instance while rescaling at every node.
 |        Note: all nodes get recalculated.
 |
 -----------------------------------------------------------------*/
int TreeCondLikes_Beagle_Rescale_All (Tree *t, int division, int chain)
{
    int                 i, j, cumulativeScaleIndex, op, opJ, opPrev;
    TreeNode            *p;
    ModelInfo           *m;
    unsigned            chil1Step, chil2Step;
    int                 *isScalerNode;
    
    m = &modelSettings[division];
    isScalerNode = m->isScalerNode[chain];
    
    op = 0;

    for (i=0; i<t->nIntNodes; i++)
        {

#if defined (BEAGLE_V3_ENABLED)
    if (t->levelPassEnabled)
        {
        p = t->intDownPassLevel[i];
        }
    else
#endif
        {
        p = t->intDownPass[i];
        }

        if (p->upDateCl == NO)
            {
            //p->upDateCl = YES;
            /* flip to the new workspace */
            FlipCondLikeSpace (m, chain, p->index);
            }
        
        /* update conditional likelihoods */
        m->operations[op].destinationPartials    = m->condLikeIndex[chain][p->index       ];
        m->operations[op].child1Partials         = m->condLikeIndex[chain][p->left->index ];
        m->operations[op].child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ];
        m->operations[op].child2Partials         = m->condLikeIndex[chain][p->right->index];
        m->operations[op].child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index];
        
        /* All partials for tips are the same across omega catigoris, thus we are doing the following two if statments.*/
        if (p->left->left== NULL)
            chil1Step=0;
        else
            chil1Step=1;
        
        if (p->right->left== NULL)
            chil2Step=0;
        else
            chil2Step=1;

        m->operations[op].destinationScaleRead = BEAGLE_OP_NONE;
        if (isScalerNode[p->index] == YES)
            {
            FlipNodeScalerSpace (m, chain, p->index);
            m->operations[op].destinationScaleWrite = m->nodeScalerIndex[chain][p->index];
            cumulativeScaleIndex  = m->siteScalerIndex[chain];
            }
        else
            {
            m->operations[op].destinationScaleWrite = BEAGLE_OP_NONE;
            cumulativeScaleIndex  = BEAGLE_OP_NONE;
            }
        
        
        
        for (j=1; j<m->nCijkParts; j++)
            {
            opJ = op + j * t->nIntNodes;
            opPrev = op + (j-1) * t->nIntNodes;

            m->operations[opJ].destinationPartials    = m->operations[opPrev].destinationPartials + 1;
            m->operations[opJ].child1Partials         = m->operations[opPrev].child1Partials + chil1Step;
            m->operations[opJ].child1TransitionMatrix = m->operations[opPrev].child1TransitionMatrix + 1;
            m->operations[opJ].child2Partials         = m->operations[opPrev].child2Partials + chil2Step;
            m->operations[opJ].child2TransitionMatrix = m->operations[opPrev].child2TransitionMatrix + 1;

            m->operations[opJ].destinationScaleRead = BEAGLE_OP_NONE;
            if ( isScalerNode[p->index] == YES )
                {
                m->operations[opJ].destinationScaleWrite = m->operations[opPrev].destinationScaleWrite + 1;
                }
            else 
                {
                m->operations[opJ].destinationScaleWrite = BEAGLE_OP_NONE;
                }
            }
        op++;
        }

    cumulativeScaleIndex  = m->siteScalerIndex[chain];
    for (j=0; j<m->nCijkParts; j++)
        {
        beagleUpdatePartials(m->beagleInstance,
                             &m->operations[j*t->nIntNodes],
                             t->nIntNodes,
                             cumulativeScaleIndex);
        cumulativeScaleIndex++;
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   TreeCondLikes_Beagle: This routine updates and rescales
|        all conditional (partial) likelihoods of a beagle instance.
|
-----------------------------------------------------------------*/
int TreeCondLikes_Beagle_Always_Rescale (Tree *t, int division, int chain)
{
    int                 i, j, op, opJ, opPrev, scaleOp;
    int                 destinationScaleRead, cumulativeScaleIndex;
    TreeNode            *p;
    ModelInfo           *m;
    unsigned            chil1Step, chil2Step;
    
    m = &modelSettings[division];
    
    op = 0;
    scaleOp = 0;

    for (i=0; i<t->nIntNodes; i++)
        {

#if defined (BEAGLE_V3_ENABLED)
    if (t->levelPassEnabled)
        {
        p = t->intDownPassLevel[i];
        }
    else
#endif
        {
        p = t->intDownPass[i];
        }
        
        /* check if conditional likelihoods need updating */
        if (p->upDateCl == YES)
            {
            /* remove old scalers */
            if (m->upDateAll == NO)
                {
                destinationScaleRead = m->nodeScalerIndex[chain][p->index];
                for (j=0; j<m->nCijkParts; j++)
                    {
                    m->scaleFactorsOps[scaleOp+j*t->nIntNodes] = destinationScaleRead;
                    destinationScaleRead++;
                    }
                scaleOp++;
                }

            /* flip to the new workspace */
            FlipCondLikeSpace (m, chain, p->index);
            FlipNodeScalerSpace (m, chain, p->index);
            
            /* update conditional likelihoods */
            m->operations[op].destinationPartials    = m->condLikeIndex[chain][p->index       ];
            m->operations[op].child1Partials         = m->condLikeIndex[chain][p->left->index ];
            m->operations[op].child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ];
            m->operations[op].child2Partials         = m->condLikeIndex[chain][p->right->index];
            m->operations[op].child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index];

            /* All partials for tips are the same across omega catigoris, thus we are doing the following two if statments.*/
            if (p->left->left== NULL && p->left->right== NULL)
                chil1Step=0;
            else
                chil1Step=1;

            if (p->right->left== NULL && p->right->right== NULL)
                chil2Step=0;
            else
                chil2Step=1;

            m->operations[op].destinationScaleWrite = m->nodeScalerIndex[chain][p->index];
            m->operations[op].destinationScaleRead = BEAGLE_OP_NONE;

            for (j=1; j<m->nCijkParts; j++)
                {
                opJ = op + j * t->nIntNodes;
                opPrev = op + (j-1) * t->nIntNodes;

                m->operations[opJ].destinationPartials    = m->operations[opPrev].destinationPartials + 1;
                m->operations[opJ].child1Partials         = m->operations[opPrev].child1Partials + chil1Step;
                m->operations[opJ].child1TransitionMatrix = m->operations[opPrev].child1TransitionMatrix + 1;
                m->operations[opJ].child2Partials         = m->operations[opPrev].child2Partials + chil2Step;
                m->operations[opJ].child2TransitionMatrix = m->operations[opPrev].child2TransitionMatrix + 1;

                m->operations[opJ].destinationScaleWrite = m->operations[opPrev].destinationScaleWrite + 1;
                m->operations[opJ].destinationScaleRead = BEAGLE_OP_NONE;
                }
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
            printf("%d %d ", op, m->operations[op].destinationPartials);
            printf("%d ", m->operations[op].destinationScaleWrite);
            printf("%d ", m->operations[op].destinationScaleRead);
            printf("%d ", m->operations[op].child1Partials);
            printf("%d ", m->operations[op].child1TransitionMatrix);
            printf("%d ", m->operations[op].child2Partials);
            printf("%d\n", m->operations[op].child2TransitionMatrix);
#endif
            op++;
            }
        } /* end of for */

    cumulativeScaleIndex  = m->siteScalerIndex[chain];
    for (j=0; j<m->nCijkParts; j++)
        {
        if (m->upDateAll == NO)
            {
            beagleRemoveScaleFactors(m->beagleInstance,
                                     &m->scaleFactorsOps[j*t->nIntNodes],
                                     scaleOp,
                                     cumulativeScaleIndex);
            }

        beagleUpdatePartials(m->beagleInstance,
                             &m->operations[j*t->nIntNodes],
                             op,
                             cumulativeScaleIndex);
        cumulativeScaleIndex++;
        }

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
        freq = 1.0 / m->numRateCats;
    else
        freq = (1.0 - pInvar) / m->numRateCats;
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
    printf("freq = %f\n", freq);
#endif
    /* TODO: cat weights only need to be set when they change */
    /* set category frequencies in beagle instance */
    if (m->numOmegaCats > 1)
        {
        omegaCatFreq = GetParamSubVals(m->omega, chain, state[chain]);
        for (i=0; i<m->nCijkParts; i++)
            {
            for (j=0; j<m->numRateCats; j++)
                m->inWeights[j] = freq * omegaCatFreq[i];
            beagleSetCategoryWeights(m->beagleInstance,
                                     m->cijkIndex[chain] + i,
                                     m->inWeights);
            }
        }
    else if (hasPInvar == YES)
        {
        for (i=0; i<m->numRateCats; i++)
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
        SAFEFREE (nSitesOfPat_Beagle);
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
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
    printf("bufferIndices = %d, eigenIndices = %d, cumulativeScaleIndices = %d\n", m->bufferIndices[0], m->eigenIndices[0], m->cumulativeScaleIndices[0]);
#endif
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
    if (*lnL > DBL_MAX || *lnL < -DBL_MAX) {
        beagleReturn = BEAGLE_ERROR_FLOATING_POINT;
    }
    if (beagleReturn == BEAGLE_ERROR_FLOATING_POINT)
        {
#   if defined (MB_PRINT_DYNAMIC_RESCALE_FAIL_STAT)
        countBeagleDynamicFail++;
        MrBayesPrint ("DEBUG INFO (not an error) countBeagleDynamicFail:%d countALL:%d\n", countBeagleDynamicFail, countALL);
#   endif
        return beagleReturn;
        }
    assert (beagleReturn == BEAGLE_SUCCESS);
    m->successCount[chain]++;
    
#if defined (DEBUG_MB_BEAGLE_MULTIPART_SITELNL)
    beagleGetSiteLogLikelihoods(m->beagleInstance, m->logLikelihoods);
    printf("lnL = ");
    for (c=0; c<m->numChars; c++)
        {
        printf("[%d] %f ", c, m->logLikelihoods[c]);
        }
    printf("\n");
#endif

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
    for (k=0; k<m->numRateCats; k++)
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
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
            printf("%d %f ", count, m->branchLengths[count]);
            printf(" %d", m->tiProbIndices[count]);
            printf(" %d", 0);
            printf(" %d\n", m->cijkIndex[chain]);
#endif
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

#if defined (BEAGLE_V3_ENABLED)

/*------------------------------------------------------------------------
|
|   InitBeagleMultiPartitionInstance: create and initialize a beagle instance for multiple partitions
|
-------------------------------------------------------------------------*/
int InitBeagleMultiPartitionInstance ()
{
    int                     i, j, k, c, s, d, *inStates, numPartAmbigTips, sizePD;
    int                     nCijkParts, numRateCats, numModelStates, numCondLikes, numScalers;
    int                     numChars, numTiProbs;
    int                     *anyDivPartAmbigTip;
    double                  *inPartials;
    ModelInfo               *m, *mPrev;
    BitsLong                *charBits;
    int                     beagleInstance;

    m = &modelSettings[0];

    if (m->beagleInstance == -99) {
        MrBayesPrint ("\n%s   Error: Consecutive MCMC runs not currently supported with BEAGLE. Please restart MrBayes.\n\n", spacer);
        exit(1);
    }

    if (m->useBeagle == NO)
        return ERROR;

    /* check if divisions can be merged into single instance */
    for (d=1; d<numCurrentDivisions; d++)
        {
        mPrev = m;
        m     = &modelSettings[d];
        if (m->useBeagle      == NO                    ||
            m->nCijkParts     != mPrev->nCijkParts     ||
            m->numRateCats   != mPrev->numRateCats   ||
            m->numModelStates != mPrev->numModelStates ||
            m->numCondLikes   != mPrev->numCondLikes   ||
            m->numScalers     != mPrev->numScalers     ||
            m->numTiProbs     != mPrev->numTiProbs)
            {
#ifdef DEBUG_MB_BEAGLE_MULTIPART
            printf("m->useBeagle      = %d\n m->nCijkParts     = %d | mPrev->nCijkParts     = %d\n m->numRateCats   = %d | mPrev->numRateCats   = %d\n m->numModelStates = %d | mPrev->numModelStates = %d\n m->numCondLikes   = %d | mPrev->numCondLikes   = %d\n m->numScalers     = %d | mPrev->numScalers     = %d\n m->numTiProbs     = %d | mPrev->numTiProbs     = %d\n",
                    m->useBeagle,
                    m->nCijkParts    ,  mPrev->nCijkParts    ,
                    m->numRateCats  ,  mPrev->numRateCats  ,
                    m->numModelStates,  mPrev->numModelStates,
                    m->numCondLikes  ,  mPrev->numCondLikes  ,
                    m->numScalers    ,  mPrev->numScalers    ,
                    m->numTiProbs    ,  mPrev->numTiProbs);
#endif
            return ERROR;
            }
        }
    
    anyDivPartAmbigTip = (int *) SafeCalloc (numLocalTaxa * numCurrentDivisions, sizeof(int));

    m = &modelSettings[0];

    nCijkParts     = (m->nCijkParts == 0 ? 1 : m->nCijkParts);
    numRateCats   = m->numRateCats;
    numModelStates = m->numModelStates;
    numCondLikes   = m->numCondLikes;
    numScalers     = m->numScalers;
    numChars   = 0;
    numTiProbs = 0;

    sizePD = nCijkParts*numCurrentDivisions;
    m->cijkIndicesAll         = (int *)    SafeCalloc (sizePD*2*numLocalTaxa, sizeof(int));
    m->categoryRateIndicesAll = (int *)    SafeCalloc (sizePD*2*numLocalTaxa, sizeof(int));
    m->operationsAll          = (BeagleOperationByPartition *) SafeCalloc(sizePD*m->numCondLikes, sizeof(BeagleOperationByPartition));
    m->tiProbIndices          = (int *)    SafeCalloc (sizePD*2*numLocalTaxa, sizeof(int));
    m->branchLengths          = (MrBFlt *) SafeCalloc (sizePD*2*numLocalTaxa, sizeof(MrBFlt));
    m->bufferIndices          = (int *)    SafeCalloc (sizePD, sizeof(int));
    m->childBufferIndices     = (int *)    SafeCalloc (sizePD, sizeof(int));
    m->childTiProbIndices     = (int *)    SafeCalloc (sizePD, sizeof(int));
    m->eigenIndices           = (int *)    SafeCalloc (sizePD, sizeof(int));
    m->cumulativeScaleIndices = (int *)    SafeCalloc (sizePD, sizeof(int));

    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];

        /* at least one eigen buffer needed */
        if (m->nCijkParts == 0)
            m->nCijkParts = 1;

        /* allocate memory used by beagle */
        m->inRates                 = (MrBFlt *) SafeCalloc (m->numRateCats, sizeof(MrBFlt));
        m->inWeights               = (MrBFlt *) SafeCalloc (m->numRateCats*m->nCijkParts, sizeof(MrBFlt));
        m->operationsByPartition   = (BeagleOperationByPartition *) SafeCalloc(m->numCondLikes*m->nCijkParts, sizeof(BeagleOperationByPartition));
        m->scaleFactorsOps         = (int *) SafeCalloc (2*numLocalTaxa*m->nCijkParts, sizeof(int));

        numChars     += m->numChars;
        numTiProbs   += m->numTiProbs;
        }

    m = &modelSettings[0];
    m->numCharsAll = numChars;
    m->logLikelihoodsAll = (MrBFlt *) SafeCalloc ((numLocalChains)*m->numCharsAll, sizeof(MrBFlt));

    numPartAmbigTips = 0;
    if (m->numStates != m->numModelStates)
        numPartAmbigTips = numLocalTaxa;
    else
        {
        for (i=0; i<numLocalTaxa; i++)
            {
#if defined(BEAGLE_V3_ENABLED)
            if (!beagleAllFloatTips)
#else
            if (1)
#endif
                {
                for (d=0; d<numCurrentDivisions; d++)
                    {
                    m = &modelSettings[d];
                    if (m->isPartAmbig[i] == YES)
                        {
                        anyDivPartAmbigTip[i] = YES;
                        numPartAmbigTips++;
                        break;
                        }
                    }
                }
            else
                {
                anyDivPartAmbigTip[i] = YES;
                numPartAmbigTips++;
                }
            }
        }

    beagleInstance = createBeagleInstance(m, nCijkParts, numRateCats, numModelStates, numCondLikes, numScalers, numChars, numTiProbs, numPartAmbigTips, -1);
    
    if (beagleInstance < 0)
        return ERROR;

    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        m->beagleInstance = beagleInstance;
        m->useBeagleMultiPartitions = YES;
        }

    /* initialize tip data */
    inStates = (int *) SafeMalloc (numChars * sizeof(int));
    if (!inStates)
        return ERROR;
    inPartials = (double *) SafeMalloc (numChars * numModelStates * sizeof(double));
    if (!inPartials)
        return ERROR;
    for (i=0; i<numLocalTaxa; i++)
        {
        if (anyDivPartAmbigTip[i] == NO)
            {
            k=0;
            for (d=0; d<numCurrentDivisions; d++)
                {
                m = &modelSettings[d];
                charBits = m->parsSets[i];
                for (c=0; c<m->numChars; c++)
                    {
                    for (s=j=0; s<m->numModelStates; s++)
                        {
                        if (IsBitSet(s, charBits))
                            {
                            inStates[k] = s;
                            j++;
                            }
                        }
                    if (j == m->numModelStates)
                        inStates[k] = j;
                    else
                        assert (j==1);
                    charBits += m->nParsIntsPerSite;
                    k++;
                    }
                }
            beagleSetTipStates(beagleInstance, i, inStates);
            }
        else /* if (m->isPartAmbig == YES) */
            {
            k = 0;
            for (d=0; d<numCurrentDivisions; d++)
                {
                m = &modelSettings[d];
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
                }
            beagleSetTipPartials(beagleInstance, i, inPartials);
            }
        }

    free (anyDivPartAmbigTip);
    free (inStates);
    free (inPartials);

    return NO_ERROR;
}


void LaunchBEAGLELogLikeMultiPartition(int* divisions, int divisionCount, int chain, MrBFlt* lnL)
{
    int       i, d, dIndex, rescaleDivisionCount, norescaleDivisionCount, beagleReturn;
    int       *isScalerNode, *rescaleDivisions, *norescaleDivisions;
    MrBFlt    rescaleLnL;
    Tree      *tree;
    TreeNode  *p;
    ModelInfo *m;
    
    if (beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS) 
        {
    
#   if defined (DEBUG_MB_BEAGLE_FLOW)
        MrBayesPrint ("ALWAYS RESCALING\n");
#   endif

        for (d=0; d<divisionCount; d++)
            {
            dIndex = divisions[d];
            m = &modelSettings[dIndex];
        
            /* Flip and copy or reset site scalers */
            FlipSiteScalerSpace(m, chain);
            
            /* rescale if there is at least one partition requiring it */
            if (m->upDateAll == YES)
                {
                for (i=0; i<m->nCijkParts; i++)
                    {           
                    beagleResetScaleFactorsByPartition(m->beagleInstance, m->siteScalerIndex[chain] + i, m->divisionIndex);
                    }
                }
            else
                {
                CopySiteScalers(m, chain);
                }
            }

        TreeTiProbs_BeagleMultiPartition(divisions, divisionCount, chain);
        TreeCondLikes_BeagleMultiPartition_Always_Rescale(divisions, divisionCount, chain);
        /* TODO integrate only divisions being updated */
        TreeLikelihood_BeagleMultiPartition(divisions, divisionCount, chain, lnL, (chainId[chain] % chainParams.numChains));
        } 
    else 
        { /* MB_BEAGLE_SCALE_DYNAMIC */

        (*lnL) = 0.0;
        rescaleLnL = 0.0;

        rescaleDivisions = (int *) SafeCalloc (numCurrentDivisions, sizeof(int));
        norescaleDivisions = (int *) SafeCalloc (numCurrentDivisions, sizeof(int));
        rescaleDivisionCount = 0;
        norescaleDivisionCount = 0;

        TreeTiProbs_BeagleMultiPartition(divisions, divisionCount, chain);

        for (d=0; d<divisionCount; d++)
            {
            dIndex = divisions[d];
            norescaleDivisions[d] = dIndex;
            m = &modelSettings[dIndex];
            tree = GetTree(m->brlens, chain, state[chain]);

            /* This flag is only valid within this block */
            m->rescaleBeagleAll = NO;        

            if (m->successCount[chain] > 1000)
                {
                m->successCount[chain] = 10;
                m->rescaleFreq[chain]++; /* increase rescaleFreq independent of whether we accept or reject new state */
                m->rescaleFreqOld = m->rescaleFreqNew = m->rescaleFreq[chain];
                for (i=0; i<tree->nIntNodes; i++)
                    {
                    p = tree->intDownPass[i];
                    if (p->upDateCl == YES) {
                         /* flip to the new workspace since TreeCondLikes_Beagle_Rescale_All() does not do it for
                            (p->upDateCl == YES) since it assumes that TreeCondLikes_Beagle_No_Rescale() did it */
                        FlipCondLikeSpace (m, chain, p->index);
                       }
                    }
                rescaleDivisions[rescaleDivisionCount++] = dIndex;
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
                printf ("NUMERICAL RESCALING FOR DIVISION %d due to successes\n", dIndex);
#endif

                norescaleDivisions[d] = -1;
                }
            else if (beagleScalingFrequency != 0 && m->beagleComputeCount[chain] % (beagleScalingFrequency) == 0)
                {
                m->rescaleFreqOld = m->rescaleFreqNew = m->rescaleFreq[chain];
                for (i=0; i<tree->nIntNodes; i++)
                    {
                    p = tree->intDownPass[i];
                    if (p->upDateCl == YES) {
                         /* flip to the new workspace since TreeCondLikes_Beagle_Rescale_All() does not do it for
                            (p->upDateCl == YES) since it assumes that TreeCondLikes_Beagle_No_Rescale() did it */
                        FlipCondLikeSpace (m, chain, p->index);
                       }
                    }
                rescaleDivisions[rescaleDivisionCount++] = dIndex;
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
                printf ("NUMERICAL RESCALING FOR DIVISION %d due to beagleScalingFrequency\n", dIndex);
#endif
                norescaleDivisions[d] = -1;
                }
            }

        if (rescaleDivisionCount < divisionCount)
            {
            i = 0;
            for (d=0; d<divisionCount; d++)
                {
                dIndex = norescaleDivisions[d];
                if (dIndex != -1)
                    {
                    norescaleDivisions[i++] = dIndex;
                    }
                }
            norescaleDivisionCount = i;

            TreeCondLikes_BeagleMultiPartition_No_Rescale(norescaleDivisions, norescaleDivisionCount, chain);

            beagleReturn = TreeLikelihood_BeagleMultiPartition(norescaleDivisions, norescaleDivisionCount, chain, lnL, (chainId[chain] % chainParams.numChains));

            /* Check if likelihood is valid */      
            if (beagleReturn == BEAGLE_ERROR_FLOATING_POINT) 
                {
                (*lnL)=0;
                for (d=0; d<norescaleDivisionCount; d++)
                    {
                    dIndex = norescaleDivisions[d];
                    m = &modelSettings[dIndex];
                    if (m->lnLike[2*chain + state[chain]] > DBL_MAX ||
                        m->lnLike[2*chain + state[chain]] < -DBL_MAX ||
                        m->lnLike[2*chain + state[chain]] != m->lnLike[2*chain + state[chain]])
                        {

                        m->rescaleFreqOld = m->rescaleFreqNew = m->rescaleFreq[chain];
                        if (m->rescaleFreqNew > 1 && m->successCount[chain] < 40)
                            {
                            if (m->successCount[chain] < 10)
                                {
                                if (m->successCount[chain] < 4)
                                    {
                                    m->rescaleFreqNew -= m->rescaleFreqNew >> 3; /* <== we cut up to 12,5% of rescaleFreq */
                                    if (m->successCount[chain] < 2)
                                        {
                                        m->rescaleFreqNew -= m->rescaleFreqNew >> 3;
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
                            m->rescaleFreqNew--;
                            m->rescaleFreqNew = (m->rescaleFreqNew ? m->rescaleFreqNew : 1);
                            }
                        m->successCount[chain] = 0;
                        rescaleDivisions[rescaleDivisionCount++] = dIndex;
                        norescaleDivisions[d] = -1;
                        } 
                    else
                        {
                            (*lnL)+= m->lnLike[2*chain + state[chain]];
                        }
                    }
                // /* rescaling all divisions */
                // for (d=0; d<divisionCount; d++)
                //     {
                //     rescaleDivisions[d] = divisions[d];
                //     }
                // rescaleDivisionCount = divisionCount;
                }
            }

        if (rescaleDivisionCount > 0) 
            {
            for (d=0; d<rescaleDivisionCount; d++)
                {
                dIndex = rescaleDivisions[d];
                m = &modelSettings[dIndex];

#   if defined (DEBUG_MB_BEAGLE_FLOW)
            MrBayesPrint ("NUMERICAL RESCALING FOR DIVISION %d\n", dIndex);
#   endif

                m->rescaleBeagleAll = YES;
                FlipSiteScalerSpace(m, chain);
                }
        
        while_loop:
            for (d=0; d<rescaleDivisionCount; d++)
                {
                dIndex = rescaleDivisions[d];
                m = &modelSettings[dIndex];

                isScalerNode = m->isScalerNode[chain];
                ResetScalersPartition (isScalerNode, tree, m->rescaleFreqNew);
                for (i=0; i<m->nCijkParts; i++) 
                    {           
                    beagleResetScaleFactorsByPartition(m->beagleInstance, m->siteScalerIndex[chain] + i, m->divisionIndex);
                    }
                }
                
                TreeCondLikes_BeagleMultiPartition_Rescale_All (rescaleDivisions, rescaleDivisionCount, chain);
                if (TreeLikelihood_BeagleMultiPartition(rescaleDivisions, rescaleDivisionCount, chain, &rescaleLnL, (chainId[chain] % chainParams.numChains)) == BEAGLE_ERROR_FLOATING_POINT)
                    {

                    for (d=0; d<rescaleDivisionCount; d++)
                        {
                        dIndex = rescaleDivisions[d];
                        m = &modelSettings[dIndex];
                        tree = GetTree(m->brlens, chain, state[chain]);
                        isScalerNode = m->isScalerNode[chain];

                        if (m->rescaleFreqNew > 1)
                            {
                            /* Swap back scalers which were swapped in TreeCondLikes_Beagle_Rescale_All() */
                            for (i=0; i<tree->nIntNodes; i++)
                                {
                                p = tree->intDownPass[i];
                                if (isScalerNode[p->index] == YES)
                                    FlipNodeScalerSpace (m, chain, p->index);
                                }
                            m->rescaleFreqNew -= m->rescaleFreqNew >> 3; /* <== we cut up to 12,5% of rescaleFreq */
                            m->rescaleFreqNew--;                      /* we cut extra 1 of rescaleFreq */
                            }
                        else
                            {
                            m->rescaleFreq[chain] = m->rescaleFreqNew;
                            rescaleDivisions[d] = -1;
                            }
                        }

                        i = 0;
                        for (d=0; d<rescaleDivisionCount; d++)
                            {
                            dIndex = rescaleDivisions[d];
                            if (dIndex != -1)
                                {
                                rescaleDivisions[i++] = dIndex;
                                }
                            }
                        rescaleDivisionCount = i;

                        if (rescaleDivisionCount > 0)
                            {
                            goto while_loop;
                            }
                    }

            for (d=0; d<rescaleDivisionCount; d++)
                {
                dIndex = rescaleDivisions[d];
                m = &modelSettings[dIndex];
                m->rescaleFreq[chain] = m->rescaleFreqNew;
                }

            (*lnL) += rescaleLnL;
            }

        for (d=0; d<divisionCount; d++)
            {
            dIndex = divisions[d];
            m = &modelSettings[dIndex];
            /* Count number of evaluations */
            m->beagleComputeCount[chain]++;
            }

        free(rescaleDivisions);
        free(norescaleDivisions);

        }

    return;
}


/*----------------------------------------------------------------
|
|   TreeTiProbs_BeagleMultiPartition: This routine updates all transition
|       probability matrices of a beagle instance across all divisions.
|
-----------------------------------------------------------------*/
int TreeTiProbs_BeagleMultiPartition (int* divisions, int divisionCount, int chain)
{
    int         i, j, k, d, dIndex, count, divisionOffset, divisionCijkIndex;
    MrBFlt      correctionFactor, theRate, baseRate, *catRate, length, *branchLengthsAll;
    TreeNode    *p;
    Tree        *t;
    ModelInfo   *m;
    int         *cijkIndicesAll, *tiProbIndicesAll, *categoryRateIndicesAll;

    cijkIndicesAll         = modelSettings[0].cijkIndicesAll;
    categoryRateIndicesAll = modelSettings[0].categoryRateIndicesAll;
    tiProbIndicesAll       = modelSettings[0].tiProbIndices;
    branchLengthsAll       = modelSettings[0].branchLengths;

    count = 0;
    divisionOffset = 0;
    for (d=0; d<divisionCount; d++)
        {
        dIndex = divisions[d];

        /* get model settings */
        m = &modelSettings[dIndex];
        t = GetTree(m->brlens, chain, state[chain]);
        divisionOffset = m->numTiProbs * m->nCijkParts * m->divisionIndex;
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
        printf("m->upDateAll %d\n", m->upDateAll);
#endif
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
        baseRate = GetRate (dIndex, chain);
        
        /* compensate for invariable sites if appropriate */
        if (m->pInvar != NULL)
            baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
            
        /* get category rates for gamma */
        if (m->shape == NULL)
            catRate = &theRate;
        else
            catRate = GetParamSubVals (m->shape, chain, state[chain]);
        
        /* get effective category rates */
        for (k=0; k<m->numRateCats; k++)
            m->inRates[k] = baseRate * catRate[k] * correctionFactor;

        /* TODO: only need to set category rates when they change */
        /* set category rates */
        beagleSetCategoryRatesWithIndex(m->beagleInstance, m->divisionIndex, m->inRates);
        
        divisionCijkIndex = m->cijkIndex[chain] + (numLocalChains + 1) * m->nCijkParts * m->divisionIndex;

        /* get ti prob indices and branch lengths to update */
        for (i=0; i<t->nNodes; i++)
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

                branchLengthsAll[count]       = length;
                tiProbIndicesAll[count]       = m->tiProbsIndex[chain][p->index] + divisionOffset;
                categoryRateIndicesAll[count] = m->divisionIndex;
                cijkIndicesAll[count]         = divisionCijkIndex;
#if defined (DEBUG_MB_BEAGLE_MULTIPART)
                {
                printf("%2d %f ", count, branchLengthsAll[count]);
                printf(" %3d", tiProbIndicesAll[count]);
                printf(" %d", categoryRateIndicesAll[count]);
                printf(" %f %f %f %f", m->inRates[0], m->inRates[1], m->inRates[2], m->inRates[3]);
                printf(" %d\n", cijkIndicesAll[count]);
                }
#endif
                count++;
                }
            }
        }

    /* TODO: only need to update branches that have changed */
    /* calculate transition probabilities */
    if (count > 0)
        {
        for (i=1; i<m->nCijkParts; i++)
            {
            for (j=0; j<count; j++)
                {
                k = j + i*count;
                branchLengthsAll[k]       = branchLengthsAll[j];
                tiProbIndicesAll[k]       = tiProbIndicesAll[j] + i;
                categoryRateIndicesAll[k] = categoryRateIndicesAll[j];
                cijkIndicesAll[k]         = cijkIndicesAll[j] + i;
                }
            }

        beagleUpdateTransitionMatricesWithMultipleModels(m->beagleInstance,
                                                         cijkIndicesAll,
                                                         categoryRateIndicesAll,
                                                         tiProbIndicesAll,
                                                         NULL,
                                                         NULL,
                                                         branchLengthsAll,
                                                         count * m->nCijkParts);

        }

    /* return success */
    return NO_ERROR;
}


/*----------------------------------------------------------------
 |
 |  TreeCondLikes_Beagle: This routine updates all conditional
 |       (partial) likelihoods of a beagle instance across all divisions while doing no rescaling.
 |      That potentialy can make final liklihood bad then calculation with rescaling needs to be done.
 |
 -----------------------------------------------------------------*/
int TreeCondLikes_BeagleMultiPartition_No_Rescale (int* divisions, int divisionCount, int chain)
{
    int                        i, j, d, dIndex, divisionOffset;
    int                        opJ, opPrev, opCountMax, opCountTotal;
    Tree                       *t;
    TreeNode                   *p;
    ModelInfo                  *m;
    unsigned                   chil1Step, chil2Step;
    int                        *isScalerNode;
    BeagleOperationByPartition *operationsAll;

    m = &modelSettings[0];
    operationsAll = m->operationsAll;
    opCountMax = 0;

    for (d=0; d<divisionCount; d++)
        {
        dIndex = divisions[d];
        m = &modelSettings[dIndex];
        t = GetTree(m->brlens, chain, state[chain]);
        divisionOffset = m->numTiProbs * m->nCijkParts * m->divisionIndex;

        isScalerNode = m->isScalerNode[chain];
        
        m->opCount = 0;

        for (i=0; i<t->nIntNodes; i++)
            {

            if (t->levelPassEnabled)
                {
                p = t->intDownPassLevel[i];
                }
            else
                {
                p = t->intDownPass[i];
                }
            
            /* check if conditional likelihoods need updating */
            if (p->upDateCl == YES)
                {
                /* flip to the new workspace */
                FlipCondLikeSpace (m, chain, p->index);
                
                /* update conditional likelihoods */
                m->operationsByPartition[m->opCount].destinationPartials    = m->condLikeIndex[chain][p->index       ];
                m->operationsByPartition[m->opCount].child1Partials         = m->condLikeIndex[chain][p->left->index ];
                m->operationsByPartition[m->opCount].child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ] + divisionOffset;
                m->operationsByPartition[m->opCount].child2Partials         = m->condLikeIndex[chain][p->right->index];
                m->operationsByPartition[m->opCount].child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index] + divisionOffset;
                
                /* All partials for tips are the same across omega categories, thus we are doing the following two if statments.*/
                if (p->left->left== NULL)
                    chil1Step=0;
                else
                    chil1Step=1;
                
                if (p->right->left== NULL)
                    chil2Step=0;
                else
                    chil2Step=1;
                
                m->operationsByPartition[m->opCount].destinationScaleWrite = BEAGLE_OP_NONE;
                m->operationsByPartition[m->opCount].cumulativeScaleIndex  = BEAGLE_OP_NONE;
                if (isScalerNode[p->index] == YES)
                    {
                    m->operationsByPartition[m->opCount].destinationScaleRead  = m->nodeScalerIndex[chain][p->index];
                    }
                else
                    {
                    m->operationsByPartition[m->opCount].destinationScaleRead  = BEAGLE_OP_NONE;
                    }
                
                for (j=1; j<m->nCijkParts; j++)
                    {
                    opJ = m->opCount + j * t->nIntNodes;
                    opPrev = m->opCount + (j-1) * t->nIntNodes;

                    m->operationsByPartition[opJ].destinationPartials    = m->operationsByPartition[opPrev].destinationPartials + 1;
                    m->operationsByPartition[opJ].child1Partials         = m->operationsByPartition[opPrev].child1Partials + chil1Step;
                    m->operationsByPartition[opJ].child1TransitionMatrix = m->operationsByPartition[opPrev].child1TransitionMatrix + 1;
                    m->operationsByPartition[opJ].child2Partials         = m->operationsByPartition[opPrev].child2Partials + chil2Step;
                    m->operationsByPartition[opJ].child2TransitionMatrix = m->operationsByPartition[opPrev].child2TransitionMatrix + 1;
                    
                    m->operationsByPartition[opJ].destinationScaleWrite = BEAGLE_OP_NONE;
                    m->operationsByPartition[opJ].cumulativeScaleIndex  = BEAGLE_OP_NONE;
                    if ( isScalerNode[p->index] == YES )
                        {
                        m->operationsByPartition[opJ].destinationScaleRead = m->operationsByPartition[opPrev].destinationScaleRead + 1;
                        }
                    else 
                        {
                        m->operationsByPartition[opJ].destinationScaleRead = BEAGLE_OP_NONE;
                        }
                    }
                m->opCount++;
                }
            }
            if (m->opCount > opCountMax)
                opCountMax = m->opCount;
        }

    opCountTotal = 0;
    for (i=0; i<opCountMax; i++)
        {
        for (j=0; j<m->nCijkParts; j++)
            {
            opJ = i + j*t->nIntNodes;
            for (d=0; d<divisionCount; d++)
                {
                dIndex = divisions[d];
                m = &modelSettings[dIndex];
                if (i<m->opCount)
                    {
                    operationsAll[opCountTotal].destinationPartials    = m->operationsByPartition[opJ].destinationPartials;
                    operationsAll[opCountTotal].destinationScaleWrite  = m->operationsByPartition[opJ].destinationScaleWrite;
                    operationsAll[opCountTotal].destinationScaleRead   = m->operationsByPartition[opJ].destinationScaleRead;
                    operationsAll[opCountTotal].child1Partials         = m->operationsByPartition[opJ].child1Partials;
                    operationsAll[opCountTotal].child1TransitionMatrix = m->operationsByPartition[opJ].child1TransitionMatrix;
                    operationsAll[opCountTotal].child2Partials         = m->operationsByPartition[opJ].child2Partials;
                    operationsAll[opCountTotal].child2TransitionMatrix = m->operationsByPartition[opJ].child2TransitionMatrix;
                    operationsAll[opCountTotal].partition              = m->divisionIndex;
                    operationsAll[opCountTotal].cumulativeScaleIndex   = m->operationsByPartition[opJ].cumulativeScaleIndex;
                    opCountTotal++;
                    }
                }
            }
        }

    beagleUpdatePartialsByPartition(m->beagleInstance,
                                    operationsAll,
                                    opCountTotal);

    return NO_ERROR;
}


/*----------------------------------------------------------------
 |
 |  TreeCondLikes_Beagle: This routine updates all conditional
 |       (partial) likelihoods of a beagle instance across all divisions while rescaling at every node.
 |        Note: all nodes get recalculated, not only tached by move.
 |
 -----------------------------------------------------------------*/
int TreeCondLikes_BeagleMultiPartition_Rescale_All (int* divisions, int divisionCount, int chain)
{
    int                        i, j, d, dIndex, divisionOffset;
    int                        opJ, opPrev, opCountMax, opCountTotal;
    Tree                       *t;
    TreeNode                   *p;
    ModelInfo                  *m;
    unsigned                   chil1Step, chil2Step;
    int                        *isScalerNode;
    BeagleOperationByPartition *operationsAll;

    m = &modelSettings[0];
    operationsAll = m->operationsAll;
    opCountMax = 0;

    t = GetTree(m->brlens, chain, state[chain]);

    for (d=0; d<divisionCount; d++)
        {
        dIndex = divisions[d];
        m = &modelSettings[dIndex];
        t = GetTree(m->brlens, chain, state[chain]);
        divisionOffset = m->numTiProbs * m->nCijkParts * m->divisionIndex;

        isScalerNode = m->isScalerNode[chain];
        
        m->opCount = 0;
        
        for (i=0; i<t->nIntNodes; i++)
            {

            if (t->levelPassEnabled)
                {
                p = t->intDownPassLevel[i];
                }
            else
                {
                p = t->intDownPass[i];
                }
            
            if (p->upDateCl == NO)
                {
                //p->upDateCl = YES;
                /* flip to the new workspace */
                FlipCondLikeSpace (m, chain, p->index);
                }
            
            /* update conditional likelihoods */
            m->operationsByPartition[m->opCount].destinationPartials    = m->condLikeIndex[chain][p->index       ];
            m->operationsByPartition[m->opCount].child1Partials         = m->condLikeIndex[chain][p->left->index ];
            m->operationsByPartition[m->opCount].child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ] + divisionOffset;
            m->operationsByPartition[m->opCount].child2Partials         = m->condLikeIndex[chain][p->right->index];
            m->operationsByPartition[m->opCount].child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index] + divisionOffset;
            
            /* All partials for tips are the same across omega catigoris, thus we are doing the following two if statments.*/
            if (p->left->left== NULL)
                chil1Step=0;
            else
                chil1Step=1;
            
            if (p->right->left== NULL)
                chil2Step=0;
            else
                chil2Step=1;

            m->operationsByPartition[m->opCount].destinationScaleRead = BEAGLE_OP_NONE;
            if (isScalerNode[p->index] == YES)
                {
                FlipNodeScalerSpace (m, chain, p->index);
                m->operationsByPartition[m->opCount].destinationScaleWrite = m->nodeScalerIndex[chain][p->index];
                m->operationsByPartition[m->opCount].cumulativeScaleIndex = m->siteScalerIndex[chain];
                }
            else
                {
                m->operationsByPartition[m->opCount].destinationScaleWrite = BEAGLE_OP_NONE;
                m->operationsByPartition[m->opCount].cumulativeScaleIndex = BEAGLE_OP_NONE;
                }
                        
            for (j=1; j<m->nCijkParts; j++)
                {
                opJ = m->opCount + j * t->nIntNodes;
                opPrev = m->opCount + (j-1) * t->nIntNodes;

                m->operationsByPartition[opJ].destinationPartials    = m->operationsByPartition[opPrev].destinationPartials + 1;
                m->operationsByPartition[opJ].child1Partials         = m->operationsByPartition[opPrev].child1Partials + chil1Step;
                m->operationsByPartition[opJ].child1TransitionMatrix = m->operationsByPartition[opPrev].child1TransitionMatrix + 1;
                m->operationsByPartition[opJ].child2Partials         = m->operationsByPartition[opPrev].child2Partials + chil2Step;
                m->operationsByPartition[opJ].child2TransitionMatrix = m->operationsByPartition[opPrev].child2TransitionMatrix + 1;

                m->operationsByPartition[opJ].destinationScaleRead = BEAGLE_OP_NONE;
                if ( isScalerNode[p->index] == YES )
                    {
                    m->operationsByPartition[opJ].destinationScaleWrite = m->operationsByPartition[opPrev].destinationScaleWrite + 1;
                    m->operationsByPartition[opJ].cumulativeScaleIndex = m->operationsByPartition[opPrev].cumulativeScaleIndex + 1;
                    }
                else 
                    {
                    m->operationsByPartition[opJ].destinationScaleWrite = BEAGLE_OP_NONE;
                    m->operationsByPartition[opJ].cumulativeScaleIndex = BEAGLE_OP_NONE;
                    }
                }
            m->opCount++;
            }
            if (m->opCount > opCountMax)
                opCountMax = m->opCount;
        }


    for (j=0; j<m->nCijkParts; j++)
        {
        opCountTotal = 0;
        for (i=0; i<opCountMax; i++)
            {
            opJ = i + j*t->nIntNodes;
            for (d=0; d<divisionCount; d++)
                {
                dIndex = divisions[d];
                m = &modelSettings[dIndex];
                if (i<m->opCount)
                    {
                    operationsAll[opCountTotal].destinationPartials    = m->operationsByPartition[opJ].destinationPartials;
                    operationsAll[opCountTotal].destinationScaleWrite  = m->operationsByPartition[opJ].destinationScaleWrite;
                    operationsAll[opCountTotal].destinationScaleRead   = m->operationsByPartition[opJ].destinationScaleRead;
                    operationsAll[opCountTotal].child1Partials         = m->operationsByPartition[opJ].child1Partials;
                    operationsAll[opCountTotal].child1TransitionMatrix = m->operationsByPartition[opJ].child1TransitionMatrix;
                    operationsAll[opCountTotal].child2Partials         = m->operationsByPartition[opJ].child2Partials;
                    operationsAll[opCountTotal].child2TransitionMatrix = m->operationsByPartition[opJ].child2TransitionMatrix;
                    operationsAll[opCountTotal].partition              = m->divisionIndex;
                    operationsAll[opCountTotal].cumulativeScaleIndex   = m->operationsByPartition[opJ].cumulativeScaleIndex;
                    opCountTotal++;
                    }
                }
            }
        beagleUpdatePartialsByPartition(m->beagleInstance,
                                        operationsAll,
                                        opCountTotal);
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   TreeCondLikes_BeagleMultiPartition: This routine updates and rescales all conditional
|       (partial) likelihoods of a beagle instance across all divisions.
|
-----------------------------------------------------------------*/
int TreeCondLikes_BeagleMultiPartition_Always_Rescale (int* divisions, int divisionCount, int chain)
{
    int                        i, j, d, dIndex, destinationScaleRead, cumulativeScaleIndex;
    int                        opJ, opPrev, opCountMax, opCountTotal, divisionOffset, scaleOp;
    Tree                       *t;
    TreeNode                   *p;
    ModelInfo                  *m;
    unsigned                   chil1Step, chil2Step;
    BeagleOperationByPartition *operationsAll;
    
    m = &modelSettings[0];
    operationsAll = m->operationsAll;
    opCountMax = 0;

    for (d=0; d<divisionCount; d++)
        {
        dIndex = divisions[d];

        /* get model settings */
        m = &modelSettings[dIndex];
        t = GetTree(m->brlens, chain, state[chain]);
        divisionOffset = m->numTiProbs * m->nCijkParts * m->divisionIndex;

        m->opCount = 0;
        scaleOp = 0;

        for (i=0; i<t->nIntNodes; i++)
            {

            if (t->levelPassEnabled)
                {
                p = t->intDownPassLevel[i];
                }
            else
                {
                p = t->intDownPass[i];
                }

            
            /* check if conditional likelihoods need updating */
            if (p->upDateCl == YES)
                {

                /* remove old scalers */
                if (m->upDateAll == NO)
                    {
                    destinationScaleRead = m->nodeScalerIndex[chain][p->index];
                    for (j=0; j<m->nCijkParts; j++)
                        {
                        m->scaleFactorsOps[scaleOp+j*t->nIntNodes] = destinationScaleRead;
                        destinationScaleRead++;
                        }
                    scaleOp++;
                    }


                /* flip to the new workspace */
                FlipCondLikeSpace (m, chain, p->index);
                FlipNodeScalerSpace (m, chain, p->index);
                
                /* update conditional likelihoods */
                m->operationsByPartition[m->opCount].destinationPartials    = m->condLikeIndex[chain][p->index       ];
                m->operationsByPartition[m->opCount].child1Partials         = m->condLikeIndex[chain][p->left->index ];
                m->operationsByPartition[m->opCount].child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ] + divisionOffset;
                m->operationsByPartition[m->opCount].child2Partials         = m->condLikeIndex[chain][p->right->index];
                m->operationsByPartition[m->opCount].child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index] + divisionOffset;

                /* All partials for tips are the same across omega catigoris, thus we are doing the following two if statments.*/
                if (p->left->left== NULL && p->left->right== NULL)
                    chil1Step=0;
                else
                    chil1Step=1;

                if (p->right->left== NULL && p->right->right== NULL)
                    chil2Step=0;
                else
                    chil2Step=1;

                m->operationsByPartition[m->opCount].destinationScaleWrite = m->nodeScalerIndex[chain][p->index];
                m->operationsByPartition[m->opCount].cumulativeScaleIndex  = m->siteScalerIndex[chain];
                m->operationsByPartition[m->opCount].destinationScaleRead = BEAGLE_OP_NONE;

                for (j=1; j<m->nCijkParts; j++)
                    {
                    opJ = m->opCount + j * t->nIntNodes;
                    opPrev = m->opCount + (j-1) * t->nIntNodes;

                    m->operationsByPartition[opJ].destinationPartials    = m->operationsByPartition[opPrev].destinationPartials + 1;
                    m->operationsByPartition[opJ].child1Partials         = m->operationsByPartition[opPrev].child1Partials + chil1Step;
                    m->operationsByPartition[opJ].child1TransitionMatrix = m->operationsByPartition[opPrev].child1TransitionMatrix + 1;
                    m->operationsByPartition[opJ].child2Partials         = m->operationsByPartition[opPrev].child2Partials + chil2Step;
                    m->operationsByPartition[opJ].child2TransitionMatrix = m->operationsByPartition[opPrev].child2TransitionMatrix + 1;

                    m->operationsByPartition[opJ].destinationScaleWrite = m->operationsByPartition[opPrev].destinationScaleWrite + 1;
                    m->operationsByPartition[opJ].cumulativeScaleIndex  = m->operationsByPartition[opPrev].cumulativeScaleIndex + 1;
                    m->operationsByPartition[opJ].destinationScaleRead = BEAGLE_OP_NONE;
                    }
                m->opCount++;
                }
            } /* end of nIntNodes for */

            if (m->opCount > opCountMax)
                {
                opCountMax = m->opCount;
                }

            if (m->upDateAll == NO)
                {
                cumulativeScaleIndex  = m->siteScalerIndex[chain];
                for (j=0; j<m->nCijkParts; j++)
                    {
                    beagleRemoveScaleFactorsByPartition(m->beagleInstance,
                                             &m->scaleFactorsOps[j*t->nIntNodes],
                                             scaleOp,
                                             cumulativeScaleIndex,
                                             m->divisionIndex);
                    cumulativeScaleIndex++;
                    }
                }
        }

    for (j=0; j<m->nCijkParts; j++)
        {
        opCountTotal = 0;
        for (i=0; i<opCountMax; i++)
            {
            opJ = i + j*t->nIntNodes;
            for (d=0; d<divisionCount; d++)
                {
                dIndex = divisions[d];
                m = &modelSettings[dIndex];
                if (i<m->opCount)
                    {
                    operationsAll[opCountTotal].destinationPartials    = m->operationsByPartition[opJ].destinationPartials;
                    operationsAll[opCountTotal].destinationScaleWrite  = m->operationsByPartition[opJ].destinationScaleWrite;
                    operationsAll[opCountTotal].destinationScaleRead   = m->operationsByPartition[opJ].destinationScaleRead;
                    operationsAll[opCountTotal].child1Partials         = m->operationsByPartition[opJ].child1Partials;
                    operationsAll[opCountTotal].child1TransitionMatrix = m->operationsByPartition[opJ].child1TransitionMatrix;
                    operationsAll[opCountTotal].child2Partials         = m->operationsByPartition[opJ].child2Partials;
                    operationsAll[opCountTotal].child2TransitionMatrix = m->operationsByPartition[opJ].child2TransitionMatrix;
                    operationsAll[opCountTotal].partition              = m->divisionIndex;
                    operationsAll[opCountTotal].cumulativeScaleIndex   = m->operationsByPartition[opJ].cumulativeScaleIndex;

#if defined (DEBUG_MB_BEAGLE_MULTIPART)
                    {
                    printf("%2d %3d ", opCountTotal, operationsAll[opCountTotal].destinationPartials);
                    printf("%3d ", operationsAll[opCountTotal].destinationScaleWrite);
                    printf("%3d ", operationsAll[opCountTotal].destinationScaleRead);
                    printf("%3d ", operationsAll[opCountTotal].child1Partials);
                    printf("%3d ", operationsAll[opCountTotal].child1TransitionMatrix);
                    printf("%3d ", operationsAll[opCountTotal].child2Partials);
                    printf("%3d ", operationsAll[opCountTotal].child2TransitionMatrix);
                    printf("%d ", operationsAll[opCountTotal].partition);
                    printf("%3d\n", operationsAll[opCountTotal].cumulativeScaleIndex);
                    }
#endif

                    opCountTotal++;
                    }
                }
            }
        beagleUpdatePartialsByPartition(m->beagleInstance,
                                        operationsAll,
                                        opCountTotal);
        }

    return NO_ERROR;
}

/**---------------------------------------------------------------------------
|
|   TreeLikelihood_BeagleMultiPartition: Accumulate the log likelihoods calculated by Beagle
|      at the root across all divisions.
|
---------------------------------------- -------------------------------------*/
int TreeLikelihood_BeagleMultiPartition (int* divisions, int divisionCount, int chain, MrBFlt *lnL, int whichSitePats)
{
    int         i, j, d, c = 0, nStates, beagleReturn, site, dIndex, divisionOffset;
    int         hasPInvar, hasAnyPInvar, hasAnyDataRestriction;
    MrBFlt      *swr, s01, s10, probOn, probOff, covBF[40], pInvar=0.0, *bs, freq, likeI, lnLikeI, diff, *omegaCatFreq, *lnLDiv;
    CLFlt       *clInvar=NULL, *nSitesOfPat;
    double      *nSitesOfPat_Beagle;
    TreeNode    *p;
    Tree        *t, *tZero;
    ModelInfo   *m;
    double      pUnobserved;

#   if defined (MB_PRINT_DYNAMIC_RESCALE_FAIL_STAT)
    static unsigned countBeagleDynamicFail=0;
    static unsigned countALL=0;
#   endif

    hasAnyPInvar = NO;
    hasAnyDataRestriction = NO;

    m = &modelSettings[0];

    for (d=0; d<divisionCount; d++)
        {
        dIndex = divisions[d];
        
        /* find model settings and nStates, pInvar, invar cond likes */
        m = &modelSettings[dIndex];
        
        divisionOffset = (numLocalChains + 1) * m->nCijkParts * m->divisionIndex;

        t = GetTree(m->brlens, chain, state[chain]);

        /* find root node */
        p = t->root->left;
        
        nStates = m->numModelStates;
        if (m->pInvar == NULL)
            {
            hasPInvar = NO;
            }
        else
            {
            hasPInvar = YES;
            hasAnyPInvar = YES;
            pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));

#if defined (DEBUG_MB_BEAGLE_MULTIPART)
            printf("pInvar[%d] = %f\n", dIndex, *(GetParamVals (m->pInvar, chain, state[chain])));
#endif
            }
    
        if (m->dataType == RESTRICTION)
            hasAnyDataRestriction = YES;

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
                                          m->cijkIndex[chain] + i + divisionOffset,
                                          bs);
            }

        /* find category frequencies */
        if (hasPInvar == NO)
            freq = 1.0 / m->numRateCats;
        else
            freq = (1.0 - pInvar) / m->numRateCats;

#if defined (DEBUG_MB_BEAGLE_MULTIPART)
        printf("freq = %f\n", freq);
#endif

        /* TODO: cat weights only need to be set when they change */
        /* set category frequencies in beagle instance */
        if (m->numOmegaCats > 1)
            {
            omegaCatFreq = GetParamSubVals(m->omega, chain, state[chain]);
            for (i=0; i<m->nCijkParts; i++)
                {
                for (j=0; j<m->numRateCats; j++)
                    m->inWeights[j] = freq * omegaCatFreq[i];
                beagleSetCategoryWeights(m->beagleInstance,
                                         m->cijkIndex[chain] + i + divisionOffset,
                                         m->inWeights);
                }
            }
        else if (hasPInvar == YES)
            {
            for (i=0; i<m->numRateCats; i++)
                m->inWeights[i] = freq;
            beagleSetCategoryWeights(m->beagleInstance,
                                     m->cijkIndex[chain] + divisionOffset,
                                     m->inWeights);
            }
        }
    
    /* TODO: pattern weights only need to be set when they change */
    /* set pattern weights in beagle instance if using dynamic reweighting */
    if (chainParams.weightScheme[0] + chainParams.weightScheme[1] > ETA)
        {
        nSitesOfPat_Beagle = (double *) SafeMalloc (m->numCharsAll * sizeof(double));
        site = 0;
        for (d=0; d<numCurrentDivisions; d++)
            {
            for (c=0; c < modelSettings[d].numChars; c++)
                {
                nSitesOfPat_Beagle[site++] = numSitesOfPat[modelSettings[d].compCharStart + c];
                }
            }
        beagleSetPatternWeights(m->beagleInstance,
                                nSitesOfPat_Beagle);
        SAFEFREE (nSitesOfPat_Beagle);
        }

    tZero = GetTree(modelSettings[0].brlens, chain, state[chain]);

    /* find root log likelihoods and scalers */
    j=0;
    for (i=0; i<m->nCijkParts; i++)
        {
        for (d=0; d<divisionCount; d++)
            {
            dIndex = divisions[d];
            m = &modelSettings[dIndex];
            t = GetTree(m->brlens, chain, state[chain]);
            p = t->root->left;


            divisionOffset = (numLocalChains + 1) * m->nCijkParts * m->divisionIndex;

            modelSettings[0].bufferIndices[j] = m->condLikeIndex[chain][p->index] + i;
            modelSettings[0].eigenIndices[j]  = m->cijkIndex[chain] + i + divisionOffset;
            modelSettings[0].cumulativeScaleIndices[j] = m->siteScalerIndex[chain] + i;

#if defined (DEBUG_MB_BEAGLE_MULTIPART)
            printf("dIndex = %d, bufferIndices = %d, eigenIndices = %d, cumulativeScaleIndices = %d\n", dIndex, modelSettings[0].bufferIndices[j], modelSettings[0].eigenIndices[j], modelSettings[0].cumulativeScaleIndices[j]);
#endif      
            if (tZero->isRooted == NO)
                {
                divisionOffset = m->numTiProbs * m->nCijkParts * m->divisionIndex;

                modelSettings[0].childBufferIndices[j] = m->condLikeIndex[chain][p->anc->index];
                modelSettings[0].childTiProbIndices[j] = m->tiProbsIndex [chain][p->index] + i + divisionOffset;

#if defined (DEBUG_MB_BEAGLE_MULTIPART)
                printf("childBufferIndices = %d, childTiProbIndices = %d\n", modelSettings[0].childBufferIndices[j], modelSettings[0].childTiProbIndices[j]);
#endif      
                }
            j++;
            }

        }

    /* reset lnL */
    *lnL = 0.0;

    /* get root log likelihoods */
    if (tZero->isRooted == YES)
        {        
        beagleReturn = beagleCalculateRootLogLikelihoodsByPartition(modelSettings[0].beagleInstance,
                                                                    modelSettings[0].bufferIndices,
                                                                    modelSettings[0].eigenIndices,
                                                                    modelSettings[0].eigenIndices,
                                                                    modelSettings[0].cumulativeScaleIndices,
                                                                    divisions,
                                                                    divisionCount,
                                                                    modelSettings[0].nCijkParts,
                                                                    modelSettings[0].logLikelihoodsAll, /* assuming numChars >= numDivisions */
                                                                    lnL);

        }
    else
        {
        beagleReturn = beagleCalculateEdgeLogLikelihoodsByPartition(
                                                    modelSettings[0].beagleInstance,
                                                    modelSettings[0].bufferIndices,
                                                    modelSettings[0].childBufferIndices,
                                                    modelSettings[0].childTiProbIndices,
                                                    NULL,
                                                    NULL,
                                                    modelSettings[0].eigenIndices,
                                                    modelSettings[0].eigenIndices,
                                                    modelSettings[0].cumulativeScaleIndices,
                                                    divisions,
                                                    divisionCount,
                                                    modelSettings[0].nCijkParts,
                                                    modelSettings[0].logLikelihoodsAll,
                                                    lnL,
                                                    NULL,
                                                    NULL,
                                                    NULL,
                                                    NULL);

        }

#   if defined (MB_PRINT_DYNAMIC_RESCALE_FAIL_STAT)
    countALL++;
#   endif

    if (*lnL > DBL_MAX || *lnL < -DBL_MAX) {
        beagleReturn = BEAGLE_ERROR_FLOATING_POINT;
    }

    for (d=0; d<divisionCount; d++)
        {
        dIndex = divisions[d];
        m = &modelSettings[dIndex];
        m->lnLike[2*chain + state[chain]] = modelSettings[0].logLikelihoodsAll[d];
        if (m->lnLike[2*chain + state[chain]] > DBL_MAX || 
            m->lnLike[2*chain + state[chain]] < -DBL_MAX ||
            m->lnLike[2*chain + state[chain]] != m->lnLike[2*chain + state[chain]])
            {
            beagleReturn = BEAGLE_ERROR_FLOATING_POINT;
            }
        else
            {
            m->successCount[chain]++;
            }
        }

    if (beagleReturn == BEAGLE_ERROR_FLOATING_POINT)
        {
#   if defined (MB_PRINT_DYNAMIC_RESCALE_FAIL_STAT)
        countBeagleDynamicFail++;
        MrBayesPrint ("DEBUG INFO (not an error) countBeagleDynamicFail:%d countALL:%d\n", countBeagleDynamicFail, countALL);
#   endif
        }

#if defined (DEBUG_MB_BEAGLE_MULTIPART_SITELNL)
        beagleGetSiteLogLikelihoods(modelSettings[0].beagleInstance, modelSettings[0].logLikelihoodsAll);
        printf("lnL = ");
        site = 0;
        for (d=0; d<numCurrentDivisions; d++)
            {
            for (c=0; c<modelSettings[d].numChars; c++)
                {
                printf("[%d] %f ", c, modelSettings[0].logLikelihoodsAll[site]);
                site++;
                }
            }
        printf("\n");
#endif

    /* accumulate logs across sites */
    if (hasAnyPInvar || hasAnyDataRestriction)
        {
        beagleGetSiteLogLikelihoods(modelSettings[0].beagleInstance, modelSettings[0].logLikelihoodsAll);

        (*lnL) = 0.0;
        for (d=0; d<divisionCount; d++)
            {
            dIndex = divisions[d];    
            m = &modelSettings[dIndex];
            lnLDiv = &m->lnLike[2*chain + state[chain]];

            if (!(*lnLDiv > DBL_MAX || *lnLDiv < -DBL_MAX || *lnLDiv != *lnLDiv))
                {
                /* find nSitesOfPat */
                nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
                
                (*lnLDiv) = 0.0;

                site = 0;
                for (i=0; i<dIndex; i++) {
                    site += modelSettings[i].numChars;
                }
                if (m->pInvar == NULL)
                    {
                    if (m->dataType == RESTRICTION)
                        {
                        pUnobserved = 0.0;
                        for (c=0; c<m->numDummyChars; c++)
                            {
                            pUnobserved +=  exp((double)modelSettings[0].logLikelihoodsAll[site]);
                            site++;
                            }
                        /* correct for absent characters */
                        (*lnLDiv) -= log (1-pUnobserved) * (m->numUncompressedChars);
                        for (; c<m->numChars; c++)
                            {
                            (*lnLDiv) += modelSettings[0].logLikelihoodsAll[site] * nSitesOfPat[c];
                            site++;
                            }
                        }
                    }
                else
                    {
                    /* has invariable category */
                    pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
                    clInvar = m->invCondLikes;

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

                    for (c=0; c<m->numChars; c++)
                        {
                        likeI = 0.0;
                        for (j=0; j<nStates; j++)
                            likeI += (*(clInvar++)) * bs[j];
                        if (likeI != 0.0)
                            {
                            lnLikeI = log(likeI * pInvar);
                            diff = lnLikeI - modelSettings[0].logLikelihoodsAll[site];
                            }
                        else
                            diff = -1000.0;
                        if (diff < -200.0)
                            (*lnLDiv) += modelSettings[0].logLikelihoodsAll[site] * nSitesOfPat[c];
                        else if (diff > 200.0)
                            (*lnLDiv) += lnLikeI * nSitesOfPat[c];
                        else
                            {
                            (*lnLDiv) += (modelSettings[0].logLikelihoodsAll[site] + log(1.0 + exp(diff))) * nSitesOfPat[c];
                            }
                        site++;
                        }       
                    }
                    (*lnL) += m->lnLike[2*chain + state[chain]];
                }
            }
            /* check for numerical errors */
            assert ((*lnL) == (*lnL));
        }
        
    return beagleReturn;
}

#endif /* BEAGLE_V3_ENABLED */

#endif /* BEAGLE_ENABLED */


void BeagleNotLinked()
{
    MrBayesPrint ("%s   BEAGLE library is not linked to this executable.\n", spacer);
}


void BeagleThreadsNotAvailable()
{
    MrBayesPrint ("%s   BEAGLE CPU threading requires v3.1 and higher of the library.\n", spacer);
}
