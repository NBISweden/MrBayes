/*
 *  MrBayes 3
 *
 *  (c) 2002-2010
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
#include <ctype.h>
#include <assert.h>
#include "mb.h"
#include "globals.h"
#include "mbbeagle.h"
#include "utils.h"
#include "mcmc.h"

void    FlipCondLikeSpace (ModelInfo *m, int chain, int nodeIndex);
void    FlipNodeScalerSpace (ModelInfo *m, int chain, int nodeIndex);
void    FlipSiteScalerSpace (ModelInfo *m, int chain);
void    ResetSiteScalers (ModelInfo *m, int chain);
void    CopySiteScalers (ModelInfo *m, int chain);

int     TreeCondLikes_Beagle (Tree *t, int division, int chain);
int     TreeLikelihood_Beagle (Tree *t, int division, int chain, MrBFlt *lnL, int whichSitePats);
int     TreeTiProbs_Beagle (Tree *t, int division, int chain);

int		TreeCondLikes_Beagle_No_Rescale (Tree *t, int division, int chain);
int		TreeCondLikes_Beagle_Rescale_All (Tree *t, int division, int chain);

extern int *chainId;
extern int numLocalChains;

//#define DEBUG_MB_BEAGLE_FLOW

#if defined (BEAGLE_ENABLED)
/*------------------------------------------------------------------------
|
|	InitBeagleInstance: create and initialize a beagle instance
|
-------------------------------------------------------------------------*/
int InitBeagleInstance (ModelInfo *m)
{
    int                     i, j, k, c, s, *inStates, numPartAmbigTips;
    double                  *inPartials;
    SafeLong                *charBits;
    BeagleInstanceDetails   details;
    long preferedFlags, requiredFlags;
	int resource;
    
    if (m->useBeagle == NO)
        return ERROR;
    
    /* at least one eigen buffer needed */
    if (m->nCijkParts == 0)
        m->nCijkParts = 1;

    /* allocate memory used by beagle */
    m->logLikelihoods          = (MrBFlt *) calloc ((numLocalChains)*m->numChars, sizeof(MrBFlt));
    m->inRates                 = (MrBFlt *) calloc (m->numGammaCats, sizeof(MrBFlt));
    m->branchLengths           = (MrBFlt *) calloc (2*numLocalTaxa, sizeof(MrBFlt));
    m->tiProbIndices           = (int *) calloc (2*numLocalTaxa, sizeof(int));
    m->inWeights               = (MrBFlt *) calloc (m->numGammaCats*m->nCijkParts, sizeof(MrBFlt));
    m->bufferIndices           = (int *) calloc (m->nCijkParts, sizeof(int));
    m->eigenIndices            = (int *) calloc (m->nCijkParts, sizeof(int));
    m->childBufferIndices      = (int *) calloc (m->nCijkParts, sizeof(int));
    m->childTiProbIndices      = (int *) calloc (m->nCijkParts, sizeof(int));
    m->cumulativeScaleIndices  = (int *) calloc (m->nCijkParts, sizeof(int));

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

	if (beagleResourceCount == 0) 
		{
		preferedFlags = beagleFlags;
		} 
		else 
		{
		resource = beagleResource[beagleInstanceCount % beagleResourceCount];
		preferedFlags = beagleFlags;		
		}
	
    requiredFlags = 0L;
    
    if (beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS)
        requiredFlags |= BEAGLE_FLAG_SCALERS_LOG;

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
        MrBayesPrint ("%s   Failed to start beagle instance\n", spacer);
        return (ERROR);
        }
    else
        {
		MrBayesPrint( "\n   Using BEAGLE resource %i:", details.resourceNumber);
#if defined (THREADS_ENABLED)
        MrBayesPrint( " (%s)\n", (tryToUseThreads ? "threaded" : "non-threaded"));
#else
		MrBayesPrint( " (non-threaded)\n");
#endif
    	MrBayesPrint( "\tRsrc Name : %s\n", details.resourceName);
    	MrBayesPrint( "\tImpl Name : %s\n", details.implName);    
    	MrBayesPrint( "\tFlags:");
    	BeaglePrintFlags(details.flags);
    	MrBayesPrint( "\n");
		beagleInstanceCount++;
		m->beagleComputeCount = (long*) calloc(sizeof(long), numLocalChains); // TODO Need to free.				
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
|	LaunchBEAGLELogLikeForDivision: calculate the log likelihood  
|		of the new state of the chain for a single division
|
-----------------------------------------------------------------*/
void LaunchBEAGLELogLikeForDivision(int chain, int d, ModelInfo* m, Tree* tree, MrBFlt* lnL)  {	
    
    if (beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS) {
	
#if defined (DEBUG_MB_BEAGLE_FLOW)
		printf("ALWAYS RESCALING\n");
#endif
        /* Flip and copy or reset site scalers */
        FlipSiteScalerSpace(m, chain);
        if (m->upDateAll == YES)
            ResetSiteScalers(m, chain);
        else
            CopySiteScalers(m, chain);

        TreeTiProbs_Beagle(tree, d, chain);
        TreeCondLikes_Beagle(tree, d, chain);
        TreeLikelihood_Beagle(tree, d, chain, lnL, (chainId[chain] % chainParams.numChains));
    } else { /* MB_BEAGLE_SCALE_DYNAMIC */
	
		/* This flag is only valid within this block */
        m->rescaleBeagleAll = NO;
        
        TreeTiProbs_Beagle(tree, d, chain);
        if (beagleScalingFrequency != 0 && 
			m->beagleComputeCount[chain] % beagleScalingFrequency == 0) { // Force recompute
#if defined (DEBUG_MB_BEAGLE_FLOW)
			printf("FORCED RESCALING\n");
#endif		
			m->rescaleBeagleAll = YES;
            FlipSiteScalerSpace(m, chain);
            ResetSiteScalers(m, chain);
            TreeCondLikes_Beagle_Rescale_All(tree, d, chain);           
        } else {
            TreeCondLikes_Beagle_No_Rescale(tree, d, chain);
        }

		/* Check if likelihood is valid */
        if (TreeLikelihood_Beagle(tree, d, chain, lnL, (chainId[chain] % chainParams.numChains)) == BEAGLE_ERROR_FLOATING_POINT) {
#if defined (DEBUG_MB_BEAGLE_FLOW)
			printf("NUMERICAL RESCALING\n");
#endif
			
            m->rescaleBeagleAll = YES;
            FlipSiteScalerSpace(m, chain);
            ResetSiteScalers(m, chain);
            TreeCondLikes_Beagle_Rescale_All (tree, d, chain);
            TreeLikelihood_Beagle(tree, d, chain, lnL, (chainId[chain] % chainParams.numChains));
        }		        
    }
	
	/* Count number of evaluations */
	m->beagleComputeCount[chain]++;
}


void BeagleAddGPUDevicesToList(int **newResourceList, int *beagleResourceCount) {		
	BeagleResourceList* beagleResources;
	int i, gpuCount;
	
	beagleResources = beagleGetResourceList();
	if (*newResourceList == NULL) {
		*newResourceList = (int*) calloc(sizeof(int), beagleResources->length);
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

void BeagleRemoveGPUDevicesFromList(int **beagleResource, int *beagleResourceCount) {
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
        MrBayesPrint("\tFlags:");
        BeaglePrintFlags(beagleResources->list[i].supportFlags);
		MrBayesPrint("\n\n");
		}
}

int BeagleCheckFlagCompatability(long inFlags) {
    if (inFlags & BEAGLE_FLAG_PROCESSOR_GPU) {
        if (inFlags & BEAGLE_FLAG_VECTOR_SSE) {
            MrBayesPrint("%s   Simultaneous use of GPU and SSE not available.\n", spacer);
            return NO;
        }
        if (inFlags & BEAGLE_FLAG_THREADING_OPENMP) {
            MrBayesPrint("%s   Simultaneous use of GPU and OpenMP not available.\n", spacer);
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

    if (inFlags & BEAGLE_FLAG_PROCESSOR_CPU)      
    	MrBayesPrint( " PROCESSOR_CPU", spacer);
    if (inFlags & BEAGLE_FLAG_PROCESSOR_GPU)      
    	MrBayesPrint( " PROCESSOR_GPU", spacer);
    if (inFlags & BEAGLE_FLAG_PROCESSOR_FPGA)     
    	MrBayesPrint( " PROCESSOR_FPGA", spacer);
    if (inFlags & BEAGLE_FLAG_PROCESSOR_CELL)     
    	MrBayesPrint( " PROCESSOR_CELL", spacer);
    if (inFlags & BEAGLE_FLAG_PRECISION_DOUBLE)   
    	MrBayesPrint( " PRECISION_DOUBLE", spacer);
    if (inFlags & BEAGLE_FLAG_PRECISION_SINGLE)   
    	MrBayesPrint( " PRECISION_SINGLE", spacer);
    if (inFlags & BEAGLE_FLAG_COMPUTATION_ASYNCH) 
    	MrBayesPrint( " COMPUTATION_ASYNCH", spacer);
    if (inFlags & BEAGLE_FLAG_COMPUTATION_SYNCH)  
    	MrBayesPrint( " COMPUTATION_SYNCH", spacer);
    if (inFlags & BEAGLE_FLAG_EIGEN_REAL)         
    	MrBayesPrint( " EIGEN_REAL", spacer);
    if (inFlags & BEAGLE_FLAG_EIGEN_COMPLEX)      
    	MrBayesPrint( " EIGEN_COMPLEX", spacer);
    if (inFlags & BEAGLE_FLAG_SCALING_MANUAL)     
    	MrBayesPrint( " SCALING_MANUAL", spacer);
    if (inFlags & BEAGLE_FLAG_SCALING_AUTO)       
    	MrBayesPrint( " SCALING_AUTO", spacer);
    if (inFlags & BEAGLE_FLAG_SCALING_ALWAYS)     
    	MrBayesPrint( " SCALING_ALWAYS", spacer);
    if (inFlags & BEAGLE_FLAG_SCALING_DYNAMIC)    
    	MrBayesPrint( " SCALING_DYNAMIC", spacer);
    if (inFlags & BEAGLE_FLAG_SCALERS_RAW)        
    	MrBayesPrint( " SCALERS_RAW", spacer);
    if (inFlags & BEAGLE_FLAG_SCALERS_LOG)        
    	MrBayesPrint( " SCALERS_LOG", spacer);
    if (inFlags & BEAGLE_FLAG_VECTOR_NONE)        
    	MrBayesPrint( " VECTOR_NONE", spacer);
    if (inFlags & BEAGLE_FLAG_VECTOR_SSE)         
    	MrBayesPrint( " VECTOR_SSE", spacer);
    if (inFlags & BEAGLE_FLAG_THREADING_NONE)     
    	MrBayesPrint( " THREADING_NONE", spacer);
    if (inFlags & BEAGLE_FLAG_THREADING_OPENMP)   
    	MrBayesPrint( " THREADING_OPENMP", spacer);

}
    
int ScheduleLogLikeForAllDivisions() {
	int d;
	int divisionsToLaunch = 0;
	ModelInfo		*m;
		
	if (numCurrentDivisions < 2) {
		return 0;
	}

	for (d=0; d<numCurrentDivisions; d++) {		
		m = &modelSettings[d];		
		if (m->upDateCl == YES)	{
			divisionsToLaunch++;
		}
	}
	return (divisionsToLaunch > 1);
}

#if defined(THREADS_ENABLED)
void *LaunchThreadLogLikeForDivision(void *arguments) {
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

MrBFlt LaunchLogLikeForAllDivisionsInParallel(int chain) {	
	int d;
	int threadError;
	pthread_t* threads;
	LaunchStruct* launchValues;
	int* wait;
	ModelInfo* m;
	MrBFlt chainLnLike;
	
	chainLnLike = 0.0;

	/* TODO Initialize only once */
	threads = (pthread_t*) malloc(sizeof(pthread_t) * numCurrentDivisions);
	launchValues = (LaunchStruct*) malloc(sizeof(LaunchStruct) * numCurrentDivisions);
	wait = (int*) malloc(sizeof(int) * numCurrentDivisions);
	
	/* Cycle through divisions and recalculate tis and cond likes as necessary. */
	/* Code below does not try to avoid recalculating ti probs for divisions    */
	/* that could share ti probs with other divisions.                          */
	for (d=0; d<numCurrentDivisions; d++)
		{
		
#if defined (BEST_MPI_ENABLED)
        if (isDivisionActive[d] == NO)
            continue;
#endif
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
			assert(0 == threadError);
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
			assert(0 == threadError);
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

/*----------------------------------------------------------------
 |
 |	TreeCondLikes_Beagle: This routine updates all conditional
 |       (partial) likelihoods of a beagle instance while doing no rescaling.
 |		That potentialy can make final liklihood bad then calculation with rescaling needs to be done.
 |
 -----------------------------------------------------------------*/
int TreeCondLikes_Beagle_No_Rescale (Tree *t, int division, int chain)
{
    int                 i, j, cumulativeScaleIndex;
    BeagleOperation     operations;
    TreeNode            *p;
    ModelInfo           *m;
	unsigned			chil1Step, chil2Step;
    
    m = &modelSettings[division];
    
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
            
			/* All partials for tips are the same across omega catigoris, thus we are doing the following two if statments.*/
			if(p->left->left== NULL && p->left->right== NULL)
				chil1Step=0;
			else
				chil1Step=1;
            
			if(p->right->left== NULL && p->right->right== NULL)
				chil2Step=0;
			else
				chil2Step=1;
            
            operations.destinationScaleWrite = BEAGLE_OP_NONE;
            cumulativeScaleIndex  = BEAGLE_OP_NONE;
            operations.destinationScaleRead  = m->nodeScalerIndex[chain][p->index];
            
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
                
                operations.destinationScaleRead++;
            }
        }
    }
    
    return NO_ERROR;
}




/*----------------------------------------------------------------
 |
 |	TreeCondLikes_Beagle: This routine updates all conditional
 |       (partial) likelihoods of a beagle instance while rescaling at every node.
 Note: all nodes get recalculated, not only tached by move.
 |
 -----------------------------------------------------------------*/
int TreeCondLikes_Beagle_Rescale_All (Tree *t, int division, int chain)
{
    int                 i, j, cumulativeScaleIndex;
    BeagleOperation     operations;
    TreeNode            *p;
    ModelInfo           *m;
	unsigned			chil1Step, chil2Step;
    
    m = &modelSettings[division];
    
    for (i=0; i<t->nIntNodes; i++)
    {
        p = t->intDownPass[i];
        
        FlipNodeScalerSpace (m, chain, p->index);
        if (p->upDateCl == NO ) {
			p->upDateCl = YES;
            /* flip to the new workspace */
            FlipCondLikeSpace (m, chain, p->index);
        }
        
        {
            /* update conditional likelihoods */
            operations.destinationPartials    = m->condLikeIndex[chain][p->index       ];
            operations.child1Partials         = m->condLikeIndex[chain][p->left->index ];
            operations.child1TransitionMatrix = m->tiProbsIndex [chain][p->left->index ];
            operations.child2Partials         = m->condLikeIndex[chain][p->right->index];
            operations.child2TransitionMatrix = m->tiProbsIndex [chain][p->right->index];
            
			/* All partials for tips are the same across omega catigoris, thus we are doing the following two if statments.*/
			if(p->left->left== NULL && p->left->right== NULL)
				chil1Step=0;
			else
				chil1Step=1;
            
			if(p->right->left== NULL && p->right->right== NULL)
				chil2Step=0;
			else
				chil2Step=1;
            
            
			//if ( 1 || p->scalerNode == YES)
            //             {
            //             //m->scalersSet[chain][p->index] = YES;
            //             operations.destinationScaleWrite = m->nodeScalerIndex[chain][p->index];
            //             cumulativeScaleIndex  = m->siteScalerIndex[chain];
            //             }
            //         else
            //             {
            //             operations.destinationScaleWrite = BEAGLE_OP_NONE;
            //             cumulativeScaleIndex  = BEAGLE_OP_NONE;
            //             }
			operations.destinationScaleRead = BEAGLE_OP_NONE;
			operations.destinationScaleWrite = m->nodeScalerIndex[chain][p->index];
            cumulativeScaleIndex  = m->siteScalerIndex[chain];
            
            
            
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
                
                operations.destinationScaleWrite++;
                cumulativeScaleIndex++;

                }
            }
        }

    return NO_ERROR;
}

#endif // BEAGLE_ENABLED

void BeagleNotLinked()
{
    MrBayesPrint("%s   BEAGLE library is not linked to this executable.\n", spacer);
}

void BeagleThreadsNotLinked()
{
	MrBayesPrint("%s   Pthreads library is not linked to this executable.\n", spacer);
} 
