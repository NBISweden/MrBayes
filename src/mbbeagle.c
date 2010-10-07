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
#include "mb.h"
#include "globals.h"
#include "mbbeagle.h"
#include "utils.h"

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
    

int BeagleCheckFlagCompatability(long inFlags) {
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
