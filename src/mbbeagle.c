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

void    FlipSiteScalerSpace (ModelInfo *m, int chain);
void    ResetSiteScalers (ModelInfo *m, int chain);
void    CopySiteScalers (ModelInfo *m, int chain);

int     TreeCondLikes_Beagle (Tree *t, int division, int chain);
int     TreeLikelihood_Beagle (Tree *t, int division, int chain, MrBFlt *lnL, int whichSitePats);
int     TreeTiProbs_Beagle (Tree *t, int division, int chain);
extern int *chainId;


/*-----------------------------------------------------------------
|
|	LaunchBEAGLELogLikeForDivision: calculate the log likelihood  
|		of the new state of the chain for a single division
|
-----------------------------------------------------------------*/
void LaunchBEAGLELogLikeForDivision(int chain, int d, ModelInfo* m, Tree* tree, MrBFlt* lnL) 
{	

	/* Flip and copy or reset site scalers */
	FlipSiteScalerSpace(m, chain);
	if (m->upDateAll == YES)
		ResetSiteScalers(m, chain);
	else
		CopySiteScalers(m, chain);

	TreeTiProbs_Beagle(tree, d, chain);
	TreeCondLikes_Beagle(tree, d, chain);
	TreeLikelihood_Beagle(tree, d, chain, lnL, (chainId[chain] % chainParams.numChains));
}


void BeagleAddGPUDevicesToList(int **newResourceList, int *beagleResourceCount) {
#if defined (BEAGLE_ENABLED)		
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
#else		
	BeagleNotLinked();
#endif		
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
#if defined (BEAGLE_ENABLED)
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
#else
	BeagleNotLinked();
#endif
}

void BeagleNotLinked()
{
    MrBayesPrint("%s   BEAGLE library is not linked to this executable.\n", spacer);
}

void BeagleThreadsNotLinked()
{
	MrBayesPrint("%s   Pthreads library is not linked to this executable.\n", spacer);
}
    

int BeagleCheckFlagCompatability(long inFlags) {
#if defined (BEAGLE_ENABLED)    
#if defined (BEAGLE_ENABLED)
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
#endif

    return YES;
#else
    BeagleNotLinked();
    return 0;
#endif
}


/*-------------------
|
|  BeaglePrintFlags: outputs beagle instance details
|
______________________*/
void BeaglePrintFlags(long inFlags) 
{
#if defined (BEAGLE_ENABLED)
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
#else
	// Do nothing
#endif
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