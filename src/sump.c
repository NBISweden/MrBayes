/*
 *  MrBayes 3
 *
 *  (c) 2002-2013
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
#include <ctype.h>
#include "mb.h"
#include "globals.h"
#include "command.h"
#include "bayes.h"
#include "sump.h"
#include "mbmath.h"
#include "mcmc.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

const char* const svnRevisionSumpC="$Rev$";   /* Revision keyword which is expended/updated by svn on each commit/update*/

/* local prototypes */
int      CompareModelProbs (const void *x, const void *y);
int		 PrintModelStats (char *fileName, char **headerNames, int nHeaders, ParameterSample *parameterSamples, int nRuns, int nSamples);
int		 PrintOverlayPlot (MrBFlt **xVals, MrBFlt **yVals, int nRows, int startingFrom, int nSamples);
int		 PrintParamStats (char *fileName, char **headerNames, int nHeaders, ParameterSample *parameterSamples, int nRuns, int nSamples);
void	 PrintPlotHeader (void);



/* AllocateParameterSamples: Allocate space for parameter samples */
int AllocateParameterSamples (ParameterSample **parameterSamples, int numRuns, int numRows, int numColumns)
{
    int     i, j;
    
    (*parameterSamples) = (ParameterSample *) SafeCalloc (numColumns, sizeof(ParameterSample));
    if (!(*parameterSamples))
        return ERROR;
    (*parameterSamples)[0].values = (MrBFlt **) SafeCalloc (numColumns * numRuns, sizeof (MrBFlt *));
    if (!((*parameterSamples)[0].values))
        {
        FreeParameterSamples(*parameterSamples);
        return ERROR;
        }
    (*parameterSamples)[0].values[0] = (MrBFlt *) SafeCalloc (numColumns * numRuns * numRows, sizeof (MrBFlt));
    for (i=1; i<numColumns; i++)
        (*parameterSamples)[i].values = (*parameterSamples)[0].values + i*numRuns;
    for (i=1; i<numRuns; i++)
        (*parameterSamples)[0].values[i] = (*parameterSamples)[0].values[0] + i*numRows;
    for (i=1; i<numColumns; i++)
        {
        for (j=0; j<numRuns; j++)
            {
            (*parameterSamples)[i].values[j] = (*parameterSamples)[0].values[0] + i*numRuns*numRows + j*numRows;
            }
        }

    return NO_ERROR;
}





/** Compare function (ModelProb) for qsort. Note reverse sort order (from larger to smaller probs) */
int CompareModelProbs (const void *x, const void *y) {

    if ((*((ModelProb *)(x))).prob > (*((ModelProb *)(y))).prob)
        return -1;
    else if ((*((ModelProb *)(x))).prob < (*((ModelProb *)(y))).prob)
        return 1;
    else
        return 0;
}





int DoSump (void)

{

	int			    i, n, nHeaders=0, numRows, numColumns, numRuns, whichIsX, whichIsY,
				    unreliable, oneUnreliable, burnin, longestHeader, len;
	MrBFlt		    mean, harm_mean;
	char		    **headerNames=NULL, temp[120];
    SumpFileInfo    fileInfo, firstFileInfo;
    ParameterSample *parameterSamples=NULL;
    FILE            *fpLstat=NULL;


#	if defined (MPI_ENABLED)
    if (proc_id != 0)
		return NO_ERROR;
#	endif


	/* tell user we are ready to go */
	if (sumpParams.numRuns == 1)
		MrBayesPrint ("%s   Summarizing parameters in file %s.p\n", spacer, sumpParams.sumpFileName);
	else if (sumpParams.numRuns == 2)
		MrBayesPrint ("%s   Summarizing parameters in files %s.run1.p and %s.run2.p\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
	else /* if (sumpParams.numRuns > 2) */
		{
		MrBayesPrint ("%s   Summarizing parameters in %d files (%s.run1.p,\n", spacer, sumpParams.numRuns, sumpParams.sumpFileName);
		MrBayesPrint ("%s      %s.run2.p, etc)\n", spacer, sumpParams.sumpFileName);
		}
	MrBayesPrint ("%s   Writing summary statistics to file %s.pstat\n", spacer, sumpParams.sumpFileName);

    if (chainParams.relativeBurnin == YES)
        MrBayesPrint ("%s   Using relative burnin ('relburnin=yes'), discarding the first %.0f %% of samples\n",
            spacer, chainParams.burninFraction*100.0, chainParams.burninFraction);
    else
        MrBayesPrint ("%s   Using absolute burnin ('relburnin=no'), discarding the first %d samples\n",
            spacer, chainParams.chainBurnIn, chainParams.chainBurnIn);

    /* Initialize to silence warning. */
	firstFileInfo.numRows = 0;
	firstFileInfo.numColumns = 0;

    /* examine input file(s) */
    for (i=0; i<sumpParams.numRuns; i++)
        {
        if (sumpParams.numRuns == 1)
            sprintf (temp, "%s.p", sumpParams.sumpFileName);
        else
            sprintf (temp, "%s.run%d.p", sumpParams.sumpFileName, i+1);

        if (ExamineSumpFile (temp, &fileInfo, &headerNames, &nHeaders) == ERROR)
	        goto errorExit;

        if (i==0)
            {
        	if (fileInfo.numRows == 0 || fileInfo.numColumns == 0)
				{
				MrBayesPrint ("%s   The number of rows or columns in file %d is equal to zero\n", spacer, temp);
				goto errorExit;
				}
            firstFileInfo = fileInfo;
            }
        else
            {
            if (firstFileInfo.numRows != fileInfo.numRows || firstFileInfo.numColumns != fileInfo.numColumns)
                {
                MrBayesPrint ("%s   First file had %d rows and %d columns while file %s had %d rows and %d columns\n",
                    spacer, firstFileInfo.numRows, firstFileInfo.numColumns, temp, fileInfo.numRows, fileInfo.numColumns);
                MrBayesPrint ("%s   MrBayes expects the same number of rows and columns in all files\n", spacer);
                goto errorExit;
                }
            }
        }

    numRows = fileInfo.numRows;
    numColumns = fileInfo.numColumns;
    numRuns = sumpParams.numRuns;

    /* allocate space to hold parameter information */
    if (AllocateParameterSamples (&parameterSamples, numRuns, numRows, numColumns) == ERROR)
        return ERROR;

    /* read samples */
    for (i=0; i<sumpParams.numRuns; i++)
        {
        /* derive file name */
        if (sumpParams.numRuns == 1)
            sprintf (temp, "%s.p", sumpParams.sumpFileName);
        else
            sprintf (temp, "%s.run%d.p", sumpParams.sumpFileName, i+1);
        
        /* read samples */    
        if (ReadParamSamples (temp, &fileInfo, parameterSamples, i) == ERROR)
            goto errorExit;
        }

	/* get length of longest header */
	longestHeader = 9; /* 9 is the length of the word "parameter" (for printing table) */
	for (i=0; i<nHeaders; i++)
		{
		len = (int) strlen(headerNames[i]);
		if (len > longestHeader)
			longestHeader = len;
		}


	/* Print header */
	PrintPlotHeader ();

    /* Print trace plots */
    if (FindHeader("Gen", headerNames, nHeaders, &whichIsX) == ERROR)
        {
		MrBayesPrint ("%s   Could not find the 'Gen' column\n", spacer);
        return ERROR;
        }
    if (FindHeader("LnL", headerNames, nHeaders, &whichIsY) == ERROR)
        {
		MrBayesPrint ("%s   Could not find the 'LnL' column\n", spacer);
        return ERROR;
        }
                    

    if (sumpParams.numRuns > 1)
		{
		if (sumpParams.allRuns == YES)
			{
			for (i=0; i<sumpParams.numRuns; i++)
				{
				MrBayesPrint ("\n%s   Samples from run %d:\n", spacer, i+1);
				if (PrintPlot (parameterSamples[whichIsX].values[i], parameterSamples[whichIsY].values[i], numRows) == ERROR)
					goto errorExit;
				}
			}
		else
            {
            if (PrintOverlayPlot (parameterSamples[whichIsX].values, parameterSamples[whichIsY].values, numRuns, 0, numRows) == ERROR)
			    goto errorExit;
            }
		}
	else
		{
		if (PrintPlot (parameterSamples[whichIsX].values[0], parameterSamples[whichIsY].values[0], numRows) == ERROR)
			goto errorExit;
        }
			
	/* calculate arithmetic and harmonic means of likelihoods */

    /* open output file */
    strncpy (temp, sumpParams.sumpOutfile, 90);
    strcat (temp, ".lstat");
    fpLstat = OpenNewMBPrintFile (temp);
    if (!fpLstat)
        goto errorExit;

    /* print unique identifier to the output file */
	if (strlen(stamp) > 1)
		fprintf (fpLstat, "[ID: %s]\n", stamp);

    /* print header */
    if (sumpParams.numRuns == 1)
        MrBayesPrintf(fpLstat, "arithmetic_mean\tharmonic_mean\tvalues_discarded\n");
    else
        MrBayesPrintf(fpLstat, "run\tarithmetic_mean\tharmonic_mean\tvalues_discarded\n");

    oneUnreliable = NO;
	for (n=0; n<sumpParams.numRuns; n++)
		{
		unreliable = NO;
		if (HarmonicArithmeticMeanOnLogs (parameterSamples[whichIsY].values[n], numRows, &mean, &harm_mean) == ERROR)
			{
			unreliable = YES;
			oneUnreliable = YES;
			}
		if (sumpParams.numRuns == 1)
			{
			MrBayesPrint ("\n");
			MrBayesPrint ("%s   Estimated marginal likelihoods for run sampled in file \"%s.p\":\n", spacer, sumpParams.sumpFileName);
			MrBayesPrint ("%s      (Use the harmonic mean for Bayes factor comparisons of models)\n", spacer, sumpParams.sumpFileName);
            MrBayesPrint ("%s      (Values are saved to the file %s.lstat)\n\n", spacer, sumpParams.sumpOutfile);
			MrBayesPrint ("%s   Arithmetic mean   Harmonic mean\n", spacer);
			MrBayesPrint ("%s   --------------------------------\n", spacer);
			if (unreliable == NO)
				MrBayesPrint ("%s     %9.2lf        %9.2lf\n", spacer, mean, harm_mean);
			else
				MrBayesPrint ("%s     %9.2lf *      %9.2lf *\n", spacer, mean, harm_mean);

            /* print to file */
            MrBayesPrintf(fpLstat, "%s\t", MbPrintNum(mean));
            MrBayesPrintf(fpLstat, "%s\t", MbPrintNum(harm_mean));
            if (unreliable == YES)
                MrBayesPrintf(fpLstat, "yes\n");
            else
                MrBayesPrintf(fpLstat, "no\n");
			}
		else
			{
			if (n == 0)
				{
				MrBayesPrint ("\n");
				MrBayesPrint ("%s   Estimated marginal likelihoods for runs sampled in files\n", spacer);
				if (sumpParams.numRuns > 2)
					MrBayesPrint ("%s      \"%s.run1.p\", \"%s.run2.p\", etc:\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
				else /* if (sumpParams.numRuns == 2) */
					MrBayesPrint ("%s      \"%s.run1.p\" and \"%s.run2.p\":\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
			    MrBayesPrint ("%s      (Use the harmonic mean for Bayes factor comparisons of models)\n\n", spacer, sumpParams.sumpFileName);
                MrBayesPrint ("%s      (Values are saved to the file %s.lstat)\n\n", spacer, sumpParams.sumpOutfile);
				MrBayesPrint ("%s   Run   Arithmetic mean   Harmonic mean\n", spacer);
				MrBayesPrint ("%s   --------------------------------------\n", spacer);
				}
			if (unreliable == NO)
				MrBayesPrint ("%s   %3d     %9.2lf        %9.2lf\n", spacer, n+1, mean, harm_mean);
			else
				MrBayesPrint ("%s   %3d     %9.2lf *      %9.2lf *\n", spacer, n+1, mean, harm_mean);

            /* print to file */
            MrBayesPrintf(fpLstat, "%d\t", n+1);
            MrBayesPrintf(fpLstat, "%s\t", MbPrintNum(mean));
            MrBayesPrintf(fpLstat, "%s\t", MbPrintNum(harm_mean));
            if (unreliable == YES)
                MrBayesPrintf(fpLstat, "yes\n");
            else
                MrBayesPrintf(fpLstat, "no\n");
			}					
		}	/* next run */
	if (sumpParams.numRuns == 1)
		{
		MrBayesPrint ("%s   --------------------------------\n", spacer);
		}
	else
		{
		if (HarmonicArithmeticMeanOnLogs (parameterSamples[whichIsY].values[0], sumpParams.numRuns*numRows, &mean, &harm_mean) == ERROR)
			{
			unreliable = YES;
			oneUnreliable = YES;
			}
		else
			unreliable = NO;
		MrBayesPrint ("%s   --------------------------------------\n", spacer);
		if (unreliable == YES)
			MrBayesPrint ("%s   TOTAL   %9.2lf *      %9.2lf *\n", spacer, mean, harm_mean);
		else
			MrBayesPrint ("%s   TOTAL   %9.2lf        %9.2lf\n", spacer, mean, harm_mean);
		MrBayesPrint ("%s   --------------------------------------\n", spacer);

        /* print total to file */
        MrBayesPrintf(fpLstat, "all\t");
        MrBayesPrintf(fpLstat, "%s\t", MbPrintNum(mean));
        MrBayesPrintf(fpLstat, "%s\t", MbPrintNum(harm_mean));
        if (unreliable == YES)
            MrBayesPrintf(fpLstat, "yes\n");
        else
            MrBayesPrintf(fpLstat, "no\n");
		}
	if (oneUnreliable == YES)
		{
		MrBayesPrint ("%s   * These estimates may be unreliable because \n", spacer);
		MrBayesPrint ("%s     some extreme values were excluded\n\n", spacer);
		}
	else
		{
		MrBayesPrint ("\n");
		}
    SafeFclose(&fpLstat);

    /* Calculate burnin */
    burnin = fileInfo.firstParamLine - fileInfo.headerLine;

    /* Print parameter information to screen and to file. */
	if (sumpParams.numRuns > 1 && sumpParams.allRuns == YES)
		{
		for (i=0; i<sumpParams.numRuns; i++)
			{
			/* print table header */
			MrBayesPrint ("\n");
			MrBayesPrint ("%s   Model parameter summaries for run sampled in file \"%s.run%d.p\":\n", spacer, sumpParams.sumpFileName, i+1);
			MrBayesPrint ("%s      (Based on %d samples out of a total of %d samples from this analysis)\n\n", spacer, numRows, numRows + burnin);
            if (PrintParamStats (sumpParams.sumpOutfile, headerNames, nHeaders, parameterSamples, numRuns, numRows) == ERROR)
				goto errorExit;
			if (PrintModelStats (sumpParams.sumpOutfile, headerNames, nHeaders, parameterSamples, numRuns, numRows) == ERROR)
				goto errorExit;
			}
        }

	MrBayesPrint ("\n");
	if (sumpParams.numRuns == 1)
		MrBayesPrint ("%s   Model parameter summaries for run sampled in file \"%s\":\n", spacer, sumpParams.sumpFileName);
	else if (sumpParams.numRuns == 2)
		{
		MrBayesPrint ("%s   Model parameter summaries over the runs sampled in files\n", spacer);
		MrBayesPrint ("%s      \"%s.run1.p\" and \"%s.run2.p\":\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
		}
	else
		{
		MrBayesPrint ("%s   Model parameter summaries over all %d runs sampled in files\n", spacer, sumpParams.numRuns);
		MrBayesPrint ("%s      \"%s.run1.p\", \"%s.run2.p\" etc:\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
		}

	if (sumpParams.numRuns == 1)
        {
		MrBayesPrint ("%s      Based on a total of %d samples out of a total of %d samples\n", spacer, numRows, numRows + burnin);
		MrBayesPrint ("%s         from this analysis.\n", spacer);
        }
	else
		{
		MrBayesPrint ("%s      Summaries are based on a total of %d samples from %d runs.\n", spacer, sumpParams.numRuns*numRows, sumpParams.numRuns);
		MrBayesPrint ("%s      Each run produced %d samples of which %d samples were included.\n", spacer, numRows + burnin, numRows);
		}
	MrBayesPrint ("%s      Parameter summaries saved to file \"%s.pstat\".\n", spacer, sumpParams.sumpOutfile);

    if (PrintParamStats (sumpParams.sumpOutfile, headerNames, nHeaders, parameterSamples, numRuns, numRows) == ERROR)
		goto errorExit;
    if (PrintModelStats (sumpParams.sumpOutfile, headerNames, nHeaders, parameterSamples, numRuns, numRows) == ERROR)
		goto errorExit;

    /* free memory */
    FreeParameterSamples(parameterSamples);
    for (i=0; i<nHeaders; i++)
        free (headerNames[i]);
    free (headerNames);

    expecting = Expecting(COMMAND);
	strcpy (spacer, "");
	
	return (NO_ERROR);
	
errorExit:

    /* free memory */
    FreeParameterSamples (parameterSamples);
    for (i=0; i<nHeaders; i++)
        free (headerNames[i]);
    free (headerNames);

    if (fpLstat)
        SafeFclose (&fpLstat);

    expecting = Expecting(COMMAND);    
	strcpy (spacer, "");

	return (ERROR);
}




int DoSumSs (void)
{

	int			    i, nHeaders=0, numRows, numColumns, numRuns, whichIsX, whichIsY,
				    longestHeader, len;
	char		    **headerNames=NULL, temp[120];
    SumpFileInfo    fileInfo, firstFileInfo;
    ParameterSample *parameterSamples=NULL;
    int             stepIndexSS,numSamplesInStepSS, stepBeginSS, stepBurnin;
    MrBFlt          *lnlp, *nextSteplnlp, *firstlnlp;
    MrBFlt          *marginalLnLSS=NULL,stepScalerSS,stepAcumulatorSS, stepLengthSS, tmpMfl; 
    int             beginPrint, countPrint;
    float           tmpf;
    MrBFlt          **plotArrayY=NULL,**plotArrayX=NULL;
    int             j, k, count;
    MrBFlt          sum;
    int             firstPass = YES;

#	if defined (MPI_ENABLED)
    if (proc_id != 0)
		return NO_ERROR;
#	endif

    chainParams.isSS=YES;

	/* tell user we are ready to go */
	if (sumssParams.numRuns == 1)
		MrBayesPrint ("%s   Summarizing parameters in file %s.p\n", spacer, sumpParams.sumpFileName);
	else if (sumssParams.numRuns == 2)
		MrBayesPrint ("%s   Summarizing parameters in files %s.run1.p and %s.run2.p\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
	else /* if (sumssParams.numRuns > 2) */
		{
		MrBayesPrint ("%s   Summarizing parameters in %d files (%s.run1.p,\n", spacer, sumssParams.numRuns, sumpParams.sumpFileName);
		MrBayesPrint ("%s      %s.run2.p, etc)\n", spacer, sumpParams.sumpFileName);
		}
	//MrBayesPrint ("%s   Writing summary statistics to file %s.pstat\n", spacer, sumpParams.sumpFileName);

    if (chainParams.relativeBurnin == YES)
        MrBayesPrint ("%s   Using relative burnin ('relburnin=yes'), discarding the first %.0f %% of samples within each step.\n",
            spacer, chainParams.burninFraction*100.0, chainParams.burninFraction);
    else
        MrBayesPrint ("%s   Using absolute burnin ('relburnin=no'), discarding the first %d samples within each step.\n",
            spacer, chainParams.chainBurnIn, chainParams.chainBurnIn);

    /* Initialize to silence warning. */
	firstFileInfo.numRows = 0;
	firstFileInfo.numColumns = 0;

    /* examine input file(s) */
    for (i=0; i<sumssParams.numRuns; i++)
        {
        if (sumssParams.numRuns == 1)
            sprintf (temp, "%s.p", sumpParams.sumpFileName);
        else
            sprintf (temp, "%s.run%d.p", sumpParams.sumpFileName, i+1);

        if (ExamineSumpFile (temp, &fileInfo, &headerNames, &nHeaders) == ERROR)
	        goto errorExit;

        if (i==0)
            {
        	if (fileInfo.numRows == 0 || fileInfo.numColumns == 0)
				{
				MrBayesPrint ("%s   The number of rows or columns in file %d is equal to zero\n", spacer, temp);
				goto errorExit;
				}
            firstFileInfo = fileInfo;
            }
        else
            {
            if (firstFileInfo.numRows != fileInfo.numRows || firstFileInfo.numColumns != fileInfo.numColumns)
                {
                MrBayesPrint ("%s   First file had %d rows and %d columns while file %s had %d rows and %d columns\n",
                    spacer, firstFileInfo.numRows, firstFileInfo.numColumns, temp, fileInfo.numRows, fileInfo.numColumns);
                MrBayesPrint ("%s   MrBayes expects the same number of rows and columns in all files\n", spacer);
                goto errorExit;
                }
            }
        }

    numRows = fileInfo.numRows;
    numColumns = fileInfo.numColumns;
    numRuns = sumssParams.numRuns;

    /* allocate space to hold parameter information */
    if (AllocateParameterSamples (&parameterSamples, numRuns, numRows, numColumns) == ERROR)
        goto errorExit;

    /* read samples */
    for (i=0; i<sumssParams.numRuns; i++)
        {
        /* derive file name */
        if (sumssParams.numRuns == 1)
            sprintf (temp, "%s.p", sumpParams.sumpFileName);
        else
            sprintf (temp, "%s.run%d.p", sumpParams.sumpFileName, i+1);
        
        /* read samples */    
        if (ReadParamSamples (temp, &fileInfo, parameterSamples, i) == ERROR)
            goto errorExit;
        }

	/* get length of longest header */
	longestHeader = 9; /* 9 is the length of the word "parameter" (for printing table) */
	for (i=0; i<nHeaders; i++)
		{
		len = (int) strlen(headerNames[i]);
		if (len > longestHeader)
			longestHeader = len;
		}


    /* Print trace plots */
    if (FindHeader("Gen", headerNames, nHeaders, &whichIsX) == ERROR)
        {
		MrBayesPrint ("%s   Could not find the 'Gen' column\n", spacer);
        goto errorExit;
        }
    if (FindHeader("LnL", headerNames, nHeaders, &whichIsY) == ERROR)
        {
		MrBayesPrint ("%s   Could not find the 'LnL' column\n", spacer);
        goto errorExit;
        }
                    


    if(chainParams.burninSS > 0)
        {
        stepBeginSS = chainParams.burninSS + 1;
        }
    else
        {
        numSamplesInStepSS = (numRows-1)/(chainParams.numStepsSS-chainParams.burninSS);
        stepBeginSS = numSamplesInStepSS + 1;
        }

    numSamplesInStepSS = (numRows - stepBeginSS)/chainParams.numStepsSS;
    if( (numRows - stepBeginSS)%chainParams.numStepsSS!=0 )
        {
        MrBayesPrint ("%s   Error:  Number of samples could not be evenly devided among steps (%d samples among %d steps). \n", spacer,(numRows - stepBeginSS),chainParams.numStepsSS);
        goto errorExit;
        }


    if( chainParams.relativeBurnin == YES )
        {
        stepBurnin = (int)(numSamplesInStepSS*chainParams.burninFraction);
        }
    else
        {
        stepBurnin = chainParams.chainBurnIn;
        if(stepBurnin >= numSamplesInStepSS )
            {
            MrBayesPrint ("%s   Error:  Burnin in each step(%d) is longer then the step itself(%d). \n", spacer,stepBurnin, numSamplesInStepSS );
            goto errorExit;               
            }
        }

    marginalLnLSS = (MrBFlt *) SafeCalloc (sumssParams.numRuns, sizeof(MrBFlt));
        /*Preparing and printing joined plot.*/
    plotArrayY = (MrBFlt **) SafeCalloc (sumssParams.numRuns+1, sizeof(MrBFlt*));
    for(i=0; i<sumssParams.numRuns+1; i++)
        plotArrayY[i] = (MrBFlt *) SafeCalloc (numSamplesInStepSS, sizeof(MrBFlt));

    plotArrayX = (MrBFlt **) SafeCalloc (sumssParams.numRuns, sizeof(MrBFlt*));
    for(i=0; i<sumssParams.numRuns; i++)
        {
        plotArrayX[i] = (MrBFlt *) SafeCalloc (numSamplesInStepSS, sizeof(MrBFlt));
        for(j=0; j<numSamplesInStepSS; j++)
            plotArrayX[i][j]=j+1;
        }

    MrBayesPrint ("%s   In total %d sampls are red from .p files.\n", spacer, numRows );
    MrBayesPrint ("\n");
    MrBayesPrint ("%s   Marginal likelihood (in natural log units) is estimated using stepping-stone sampling\n", spacer );
    MrBayesPrint ("%s   based on %d steps with %d samples within each step. \n", spacer, chainParams.numStepsSS, numSamplesInStepSS );
    MrBayesPrint ("%s   First %d samples (including generation 0) are discarded as initial burn-in.\n", spacer, stepBeginSS);
        if(chainParams.startFromPriorSS==YES)
            MrBayesPrint ("%s   Sampling is assumed have being done from prior to posterior.\n", spacer);
        else
            {
            MrBayesPrint ("%s   Sampling is assumed have being done from posterior to prior.\n", spacer);
            }

sumssTable:

    MrBayesPrint ("\n\n%s   Step contribution table.\n\n",spacer);
    MrBayesPrint ("   Columns in the table: \n");
    MrBayesPrint ("   Step -- Index of the step \n");
    MrBayesPrint ("   runX -- Contribution to the marginal log likelihood of run X, i.e. marginal \n"); 
    MrBayesPrint ("           log likelihood for run X is the sum across all steps in column runX.\n\n");

    if( firstPass == YES && chainParams.relativeBurnin == YES )
        MrBayesPrint ("%s   The table entrances are based on samples excluding burn-in %d samples  (%d%%)    \n", spacer, stepBurnin,(int)(100*chainParams.burninFraction) );
    else
        MrBayesPrint ("%s   The table entrances are based on samples excluding burn-in %d samples      \n", spacer, stepBurnin);
    MrBayesPrint ("%s   discarded at the begining of each step.  \n\n", spacer);

    //MrBayesPrint ("%s       Run   Marginal likelihood (ln)\n",spacer);
    //MrBayesPrint ("%s       ------------------------------\n",spacer);
    MrBayesPrint ("   Step");
    for (j=0; j<sumssParams.numRuns ; j++)
        {
        if(j<9)
            MrBayesPrint (" ");
        MrBayesPrint ("      run%d", j+1);
        }
    MrBayesPrint ("\n");
    for(i=0; i<sumssParams.numRuns; i++)
        {
        marginalLnLSS[i] = 0.0;  
        }
    for(stepIndexSS = chainParams.numStepsSS-1; stepIndexSS>=0; stepIndexSS--)   
        {
        if(chainParams.startFromPriorSS==YES)
            {
            stepLengthSS = BetaQuantile( chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-stepIndexSS)/(MrBFlt)chainParams.numStepsSS)-BetaQuantile( chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-1-stepIndexSS)/(MrBFlt)chainParams.numStepsSS);
            }
        else
            {
            stepLengthSS = BetaQuantile ( chainParams.alphaSS, 1.0, (MrBFlt)(stepIndexSS+1)/(MrBFlt)chainParams.numStepsSS) - BetaQuantile ( chainParams.alphaSS, 1.0, (MrBFlt)stepIndexSS/(MrBFlt)chainParams.numStepsSS);
            }
        MrBayesPrint ("   %3d   ", chainParams.numStepsSS-stepIndexSS);
        for(i=0; i<sumssParams.numRuns; i++)
            {
            lnlp = parameterSamples[whichIsY].values[i] + stepBeginSS + (chainParams.numStepsSS-stepIndexSS-1)*numSamplesInStepSS;
            nextSteplnlp = lnlp+numSamplesInStepSS;
            lnlp+= stepBurnin;
            stepAcumulatorSS = 0.0;
            stepScalerSS = *lnlp*stepLengthSS;
            while( lnlp<nextSteplnlp )
                {
               if( *lnlp*stepLengthSS > stepScalerSS + 200.0 )
                    {
                    // adjust scaler;
                    stepAcumulatorSS /= exp( *lnlp*stepLengthSS - 100.0 - stepScalerSS ); 
                    stepScalerSS= *lnlp*stepLengthSS - 100.0;
                    }
                stepAcumulatorSS += exp( *lnlp*stepLengthSS - stepScalerSS );
                lnlp++;
                }
            tmpMfl = (log( stepAcumulatorSS/(numSamplesInStepSS-stepBurnin) ) + stepScalerSS);
            MrBayesPrint (" %10.3lf", tmpMfl);
            marginalLnLSS[i] += tmpMfl;
            }
        MrBayesPrint ("\n");
        //MrBayesPrint ("%s       %3d    %9.2f   \n", spacer, i+1, marginalLnLSS );
        }
    MrBayesPrint ("         ");
    for (j=0; j<sumssParams.numRuns ; j++)
        {
        if(j<9)
            MrBayesPrint ("-");
        MrBayesPrint ("----------");
        }
    MrBayesPrint ("\n");
    MrBayesPrint ("   Sum:  ");
    for (j=0; j<sumssParams.numRuns ; j++)
        MrBayesPrint (" %10.3lf", marginalLnLSS[j]);
        
	MrBayesPrint ("\n");
/*
            if (sumssParams.numRuns > 1)
		    {
		    MrBayesPrint ("%s   Below are rough plots of the generations (x-axis) during burn in  \n", spacer);
		    MrBayesPrint ("%s   phase versus the log probability of observing the data (y-axis).  \n", spacer);
		    MrBayesPrint ("%s   You can use these graphs to determine if the burn in for your SS  \n", spacer);
		    MrBayesPrint ("%s   analysis was sufficiant. The log probability suppose to plateau   \n", spacer);
		    MrBayesPrint ("%s   indicating that you may be at stationarity by the time you finish \n", spacer);
		    MrBayesPrint ("%s   burn in phase. This burn in, unlike burn in within each step, is  \n", spacer);
		    MrBayesPrint ("%s   fixed and can not be changed.                                     \n", spacer);
		    }
	    else
		    {
		    MrBayesPrint ("%s   Below is a rough plot of the generations (x-axis) during burn in  \n", spacer);
		    MrBayesPrint ("%s   phase versus the log probability of observing the data (y-axis).  \n", spacer);
		    MrBayesPrint ("%s   You can use these graph to determine if the burn in for your SS   \n", spacer);
		    MrBayesPrint ("%s   analysis was sufficiant. The log probability suppose to plateau   \n", spacer);
		    MrBayesPrint ("%s   indicating that you may be at stationarity by the time you finish \n", spacer);
		    MrBayesPrint ("%s   burn in phase. This burn in, unlike burn in within each step, is  \n", spacer);
		    MrBayesPrint ("%s   fixed and can not be changed.                                     \n", spacer);
		    }
            */

    if( firstPass == NO )
        goto sumssExitOptions;

    sumssStepPlot:

    MrBayesPrint ("\n\n%s   Step plot(s).\n",spacer);
    while(1)
        {
         MrBayesPrint ("\n");
        if( sumssParams.stepToPlot == 0 )
            {
            beginPrint=(int)(sumssParams.discardFraction*stepBeginSS);
            countPrint=stepBeginSS-beginPrint;
            MrBayesPrint ("%s   Ploting step 0, i.e initial burn-in phase consisting of %d samples.\n", spacer,stepBeginSS);
            MrBayesPrint ("%s   According to 'Discardfrac=%.2f', first %d samples are not ploted.\n", spacer,sumssParams.discardFraction,beginPrint);
            }
        else
            {
            if( sumssParams.stepToPlot > chainParams.numStepsSS )
                {
                MrBayesPrint ("%s   Chosen index of step to print %d is out of range of step indices[0,...,%d].\n", spacer,sumssParams.stepToPlot,chainParams.numStepsSS);
                goto errorExit;
                }
            beginPrint=stepBeginSS+(sumssParams.stepToPlot-1)*numSamplesInStepSS + (int)(sumssParams.discardFraction*numSamplesInStepSS);
            countPrint=numSamplesInStepSS-(int)(sumssParams.discardFraction*numSamplesInStepSS);
            MrBayesPrint ("%s   Ploting step %d consisting of %d samples.\n", spacer,sumssParams.stepToPlot,numSamplesInStepSS);
            MrBayesPrint ("%s   According to 'Discardfrac=%.2f', first %d samples are not ploted.\n", spacer,sumssParams.discardFraction,(int)(sumssParams.discardFraction*numSamplesInStepSS));
            }


        if (sumssParams.numRuns > 1)
		    {
		    if (sumpParams.allRuns == YES)
			    {
			    for (i=0; i<sumssParams.numRuns; i++)
				    {
				    MrBayesPrint ("\n%s   Samples from run %d:\n", spacer, i+1);
				    if (PrintPlot (parameterSamples[whichIsX].values[i]+beginPrint, parameterSamples[whichIsY].values[i]+beginPrint, countPrint) == ERROR)
					    goto errorExit;
				    }
			    }
		    else
                {
                if (PrintOverlayPlot (parameterSamples[whichIsX].values, parameterSamples[whichIsY].values, numRuns, beginPrint, countPrint) == ERROR)
			        goto errorExit;
                }
		    }
	    else
		    {
		    if (PrintPlot (parameterSamples[whichIsX].values[0]+beginPrint, parameterSamples[whichIsY].values[0]+beginPrint, countPrint) == ERROR)
			    goto errorExit;
            }

        if( sumssParams.askForMorePlots == NO || firstPass == YES )
            break;

        MrBayesPrint (" You can choose to print new step plots for different steps or discard fractions.\n");
        MrBayesPrint (" Allowed range of 'Steptoplot' are from 0 to %d.\n", chainParams.numStepsSS);
        MrBayesPrint (" If the next entered value is negative, 'sumss' will stop printing step plots.\n");
        MrBayesPrint (" If the next entered value is positive, but out of range, you will be offered\n");
        MrBayesPrint (" to change paramiter 'Discardfrac' of 'sumss'.\n");
        MrBayesPrint (" Enter new step number 'Steptoplot':");
        k = scanf("%d",&j);
        if(j < 0 )
            break;
        if(j > chainParams.numStepsSS)
            {
            do
                {
                MrBayesPrint (" Enter new value for 'Discardfrac', should be in range 0.0 to 1.0:");
                k = scanf("%f",&tmpf);
                sumssParams.discardFraction =  (MrBFlt)tmpf;
                }
            while(sumssParams.discardFraction < 0.0 || sumssParams.discardFraction > 1.0);
            }
        else
            sumssParams.stepToPlot=j;
    }

    if( firstPass == NO )
        goto sumssExitOptions;

	sumssJoinedPlot:		

    MrBayesPrint ("\n\n%s   Joined plot(s).\n",spacer);
    while(1)
        {
        MrBayesPrint ("\n");
        MrBayesPrint ("%s   Joined plot of %d samples of all steps together. 'smoothing' is set to:%d\n", spacer,numSamplesInStepSS,sumssParams.smoothing);
        MrBayesPrint ("%s   According to step burn-in, first %d samples are not ploted.\n", spacer,stepBurnin);

        for(i=0; i<sumssParams.numRuns; i++)
            {
            for(j=stepBurnin;j<numSamplesInStepSS;j++)
                plotArrayY[sumssParams.numRuns][j]=0.0;
            lnlp= parameterSamples[whichIsY].values[i] + stepBeginSS;
            nextSteplnlp=lnlp;
            for(stepIndexSS = chainParams.numStepsSS-1; stepIndexSS>0; stepIndexSS--)
                {
                firstlnlp=plotArrayY[sumssParams.numRuns] + stepBurnin;
                lnlp+=stepBurnin;
                nextSteplnlp += numSamplesInStepSS;
                while( lnlp<nextSteplnlp )
                    {
                    *firstlnlp+=*lnlp;
                    firstlnlp++;
                    lnlp++;
                    }
                }
            for(j=stepBurnin;j<numSamplesInStepSS;j++)
                {
                sum=0.0;
                count=0;
                for(k=j-sumssParams.smoothing;k<=j+sumssParams.smoothing;k++)
                    {
                    if(k>=stepBurnin && k<numSamplesInStepSS)
                        {
                        sum += plotArrayY[sumssParams.numRuns][k];
                        count++;
                        }
                    }
                plotArrayY[i][j] = sum/count;
                /*
                if( max < plotArrayY[i][j])
                    max=plotArrayY[i][j];
                    */
                }
        /*  for(j=stepBurnin;j<numSamplesInStepSS;j++)
                {
                plotArrayY[i][j] /= max;
                }*/
            }

        beginPrint=stepBurnin;
        countPrint=numSamplesInStepSS-stepBurnin;

        if (sumssParams.numRuns > 1)
		    {
		    if (sumpParams.allRuns == YES)
			    {
			    for (i=0; i<sumssParams.numRuns; i++)
				    {
				    MrBayesPrint ("\n%s   Samples from run %d:\n", spacer, i+1);
				    if (PrintPlot (plotArrayX[i]+beginPrint, plotArrayY[i]+beginPrint, countPrint) == ERROR)
					    goto errorExit;
				    }
			    }
		    else
                {
                if (PrintOverlayPlot (plotArrayX, plotArrayY, numRuns, beginPrint, countPrint) == ERROR)
			        goto errorExit;
                }
		    }
	    else
		    {
		    if (PrintPlot (plotArrayX[0]+beginPrint, plotArrayY[0]+beginPrint, countPrint) == ERROR)
			    goto errorExit;
            }

        if( sumssParams.askForMorePlots == NO || firstPass == YES )
            break;

        MrBayesPrint (" You can choose to print new joined plots with different step burn-in or smoothing.\n");
        MrBayesPrint (" Allowed range of step burn-in values are from 0 to %d.\n", numSamplesInStepSS-1);
        MrBayesPrint (" If the next entered value is negative, 'sumss' will stop printing joined plots.\n");
        MrBayesPrint (" If the next entered value is positive, but out of range, you will be offered\n");
        MrBayesPrint (" to change 'Smoothing'.\n");
        MrBayesPrint (" Enter new step burn-in:");
        k = scanf("%d",&j);
        if(j < 0 )
            break;
        if(j >= numSamplesInStepSS)
            {
            MrBayesPrint (" Enter new value for 'Smoothing':");
            k = scanf("%d",&j);
            sumssParams.smoothing = abs(j);
            }
        else
            stepBurnin=j;
    }

    firstPass = NO;
sumssExitOptions:
    if(sumssParams.askForMorePlots == YES )
        {
        MrBayesPrint ("\n");
        MrBayesPrint (" Sumss is interactive, because of parameter 'Askmore=YES' setting. \n");
        MrBayesPrint (" What would you like to do next?\n");
        MrBayesPrint ("   1) Print updated table according to new step burn-in.\n");
        MrBayesPrint ("   2) Print Step plot(s).\n");
        MrBayesPrint ("   3) Print Joined plot(s).\n");
        MrBayesPrint ("   4) Exit 'sumss'.\n");
        MrBayesPrint (" Enter a number that corresponds to one of the options:");
        do
            {
            k = scanf("%d",&j);
            }while(j<1 || j>4);

        if(j == 1)
            {
            MrBayesPrint (" Allowed range of step burn-in values are from 0 to %d\n", numSamplesInStepSS-1);
            MrBayesPrint (" Current step burn-in value is:%d\n", stepBurnin);
            MrBayesPrint (" Enter new step burn-in:");
            do
                {
                k = scanf("%d",&stepBurnin);
                }
            while(stepBurnin < 0 || stepBurnin > numSamplesInStepSS-1);
            MrBayesPrint ("\n"); 
            goto sumssTable;
            }
        else if(j == 2)
            {
            goto sumssStepPlot;
            }
        else if(j == 3)
            goto sumssJoinedPlot; 

        }
 
    /* free memory */
    FreeParameterSamples(parameterSamples);
    for (i=0; i<nHeaders; i++)
        free (headerNames[i]);
    free (headerNames);

    expecting = Expecting(COMMAND);
	strcpy (spacer, "");
    chainParams.isSS=NO;
    for(i=0; i<sumssParams.numRuns+1; i++)
        free(plotArrayY[i]);
    free(plotArrayY);
    for(i=0; i<sumssParams.numRuns; i++)
        free(plotArrayX[i]);
    free(plotArrayX);
    free(marginalLnLSS);
	
	return (NO_ERROR);
	
errorExit:

    /* free memory */
    FreeParameterSamples (parameterSamples);
    if( headerNames!=NULL )
        for (i=0; i<nHeaders; i++)
            free (headerNames[i]);
    free (headerNames);

    expecting = Expecting(COMMAND);    
	strcpy (spacer, "");
    chainParams.isSS=NO;
    if( plotArrayY!=NULL )
        for(i=0; i<sumssParams.numRuns+1; i++)
            free(plotArrayY[i]);
    free(plotArrayY);
    if( plotArrayX!=NULL )
        for(i=0; i<sumssParams.numRuns; i++)
            free(plotArrayX[i]);
    free(plotArrayX);
    free(marginalLnLSS);

	return (ERROR);
}





int DoSumpParm (char *parmName, char *tkn)

{

	int			tempI;
    MrBFlt      tempD;
	char		tempStr[100];

	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else
		{
		if (!strcmp(parmName, "Xxxxxxxxxx"))
			{
			expecting  = Expecting(PARAMETER);
			expecting |= Expecting(SEMICOLON);
			}
		/* set Filename (sumpParams.sumpFileName) ***************************************************/
		else if (!strcmp(parmName, "Filename"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
                if(strlen(tkn)>99 && (strchr(tkn,' ')-tkn) > 99 )
                    {
                    MrBayesPrint ("%s   Maximum allowed length of file name is 99 characters. The given name:\n", spacer);
                    MrBayesPrint ("%s      '%s'\n", spacer,tkn);
                    return (ERROR);
                    } 
				sscanf (tkn, "%s", tempStr);
				strcpy (sumpParams.sumpFileName, tempStr);
				strcpy (sumpParams.sumpOutfile, tempStr);
				MrBayesPrint ("%s   Setting sump filename and output file name to %s\n", spacer, sumpParams.sumpFileName);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Outputname (sumpParams.sumpOutfile) *******************************************************/
		else if (!strcmp(parmName, "Outputname"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
                if(strlen(tkn)>99 && (strchr(tkn,' ')-tkn) > 99 )
                    {
                    MrBayesPrint ("%s   Maximum allowed length of file name is 99 characters. The given name:\n", spacer);
                    MrBayesPrint ("%s      '%s'\n", spacer,tkn);
                    return (ERROR);
                    }
				sscanf (tkn, "%s", tempStr);
				strcpy (sumpParams.sumpOutfile, tempStr);
				MrBayesPrint ("%s   Setting sump output file name to \"%s\"\n", spacer, sumpParams.sumpOutfile);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Relburnin (chainParams.relativeBurnin) ********************************************************/
		else if (!strcmp(parmName, "Relburnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						chainParams.relativeBurnin = YES;
					else
						chainParams.relativeBurnin = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for Relburnin\n", spacer);
					return (ERROR);
					}
				if (chainParams.relativeBurnin == YES)
					MrBayesPrint ("%s   Using relative burnin (a fraction of samples discarded).\n", spacer);
				else
					MrBayesPrint ("%s   Using absolute burnin (a fixed number of samples discarded).\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				{
				return (ERROR);
				}
			}
		/* set Burnin (chainParams.chainBurnIn) ***********************************************************/
		else if (!strcmp(parmName, "Burnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
                chainParams.chainBurnIn = tempI;
				MrBayesPrint ("%s   Setting burn-in to %d\n", spacer, chainParams.chainBurnIn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				{
				return (ERROR);
				}
			}
		/* set Burninfrac (chainParams.burninFraction) ************************************************************/
		else if (!strcmp(parmName, "Burninfrac"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);
				if (tempD < 0.01)
					{
					MrBayesPrint ("%s   Burnin fraction too low (< 0.01)\n", spacer);
					return (ERROR);
					}
				if (tempD > 0.50)
					{
					MrBayesPrint ("%s   Burnin fraction too high (> 0.50)\n", spacer);
					return (ERROR);
					}
                chainParams.burninFraction = tempD;
				MrBayesPrint ("%s   Setting burnin fraction to %.2f\n", spacer, chainParams.burninFraction);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else 
				{
				return (ERROR);
				}
			}
		/* set Minprob (sumpParams.minProb) ************************************************************/
		else if (!strcmp(parmName, "Minprob"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);
				if (tempD > 0.50)
					{
					MrBayesPrint ("%s   Minprob too high (it should be smaller than 0.50)\n", spacer);
					return (ERROR);
					}
                sumpParams.minProb = tempD;
				MrBayesPrint ("%s   Setting minprob to %1.3f\n", spacer, sumpParams.minProb);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else 
				{
				return (ERROR);
				}
			}
		/* set Nruns (sumpParams.numRuns) *******************************************************/
		else if (!strcmp(parmName, "Nruns"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
				if (tempI < 1)
					{
					MrBayesPrint ("%s   Nruns must be at least 1\n", spacer);
					return (ERROR);
					}
				else
					{
					sumpParams.numRuns = tempI;
					MrBayesPrint ("%s   Setting sump nruns to %d\n", spacer, sumpParams.numRuns);
					expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
					}
				}
			else
				return (ERROR);
			}
		/* set Hpd (sumpParams.HPD) ********************************************************/
		else if (!strcmp(parmName, "Hpd"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumpParams.HPD = YES;
					else
						sumpParams.HPD = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for Hpd\n", spacer);
					return (ERROR);
					}
				if (sumpParams.HPD == YES)
					MrBayesPrint ("%s   Reporting 95 %% region of Highest Posterior Density (HPD).\n", spacer);
				else
					MrBayesPrint ("%s   Reporting median interval containing 95 %% of values.\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				{
				return (ERROR);
				}
			}
		/* set Allruns (sumpParams.allRuns) ********************************************************/
		else if (!strcmp(parmName, "Allruns"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumpParams.allRuns = YES;
					else
						sumpParams.allRuns = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for allruns (valid arguments are 'yes' and 'no')\n", spacer);
					return (ERROR);
					}
				if (sumpParams.allRuns == YES)
					MrBayesPrint ("%s   Setting sump to print information for each run\n", spacer);
				else
					MrBayesPrint ("%s   Setting sump to print only summary information for all runs\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		else
			return (ERROR);
		}

	return (NO_ERROR);

}





int DoSumSsParm (char *parmName, char *tkn)

{

	int			tempI;
    MrBFlt      tempD;
	char		tempStr[100];

	if (expecting == Expecting(PARAMETER))
		{
		expecting = Expecting(EQUALSIGN);
		}
	else
		{
		if (!strcmp(parmName, "Xxxxxxxxxx"))
			{
			expecting  = Expecting(PARAMETER);
			expecting |= Expecting(SEMICOLON);
			}
		/* set Filename (sumpParams.sumpFileName) ***************************************************/
		else if (!strcmp(parmName, "Filename"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
                if(strlen(tkn)>99 && (strchr(tkn,' ')-tkn) > 99 )
                    {
                    MrBayesPrint ("%s   Maximum allowed length of file name is 99 characters. The given name:\n", spacer);
                    MrBayesPrint ("%s      '%s'\n", spacer,tkn);
                    return (ERROR);
                    } 
				sscanf (tkn, "%s", tempStr);
				strcpy (sumpParams.sumpFileName, tempStr);
				strcpy (sumpParams.sumpOutfile, tempStr);
				MrBayesPrint ("%s   Setting sump filename and output file name to %s\n", spacer, sumpParams.sumpFileName);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Outputname (sumpParams.sumpOutfile) *******************************************************/
		/*else if (!strcmp(parmName, "Outputname"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
                if(strlen(tkn)>99 && (strchr(tkn,' ')-tkn) > 99 )
                    {
                    MrBayesPrint ("%s   Maximum allowed length of file name is 99 characters. The given name:\n", spacer);
                    MrBayesPrint ("%s      '%s'\n", spacer,tkn);
                    return (ERROR);
                    }
				sscanf (tkn, "%s", tempStr);
				strcpy (sumpParams.sumpOutfile, tempStr);
				MrBayesPrint ("%s   Setting sump output file name to \"%s\"\n", spacer, sumpParams.sumpOutfile);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}*/
		/* set Relburnin (chainParams.relativeBurnin) ********************************************************/
		else if (!strcmp(parmName, "Relburnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						chainParams.relativeBurnin = YES;
					else
						chainParams.relativeBurnin = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for Relburnin\n", spacer);
					return (ERROR);
					}
				if (chainParams.relativeBurnin == YES)
					MrBayesPrint ("%s   Using relative burnin (a fraction of samples discarded).\n", spacer);
				else
					MrBayesPrint ("%s   Using absolute burnin (a fixed number of samples discarded).\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				{
				return (ERROR);
				}
			}
		/* set Burnin (chainParams.chainBurnIn) ***********************************************************/
		else if (!strcmp(parmName, "Burnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
                chainParams.chainBurnIn = tempI;
				MrBayesPrint ("%s   Setting burn-in to %d\n", spacer, chainParams.chainBurnIn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				{
				return (ERROR);
				}
			}
		/* set Burninfrac (chainParams.burninFraction) ************************************************************/
		else if (!strcmp(parmName, "Burninfrac"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);
				if (tempD < 0.01)
					{
					MrBayesPrint ("%s   Burnin fraction too low (< 0.01)\n", spacer);
					return (ERROR);
					}
				if (tempD > 0.50)
					{
					MrBayesPrint ("%s   Burnin fraction too high (> 0.50)\n", spacer);
					return (ERROR);
					}
                chainParams.burninFraction = tempD;
				MrBayesPrint ("%s   Setting burnin fraction to %.2f\n", spacer, chainParams.burninFraction);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else 
				{
				return (ERROR);
				}
			}
		/* set Nruns (sumssParams.numRuns) *******************************************************/
		else if (!strcmp(parmName, "Nruns"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
				if (tempI < 1)
					{
					MrBayesPrint ("%s   Nruns must be at least 1\n", spacer);
					return (ERROR);
					}
				else
					{
					sumssParams.numRuns = tempI;
					MrBayesPrint ("%s   Setting sumss nruns to %d\n", spacer, sumssParams.numRuns);
					expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
					}
				}
			else
				return (ERROR);
			}
		/* set Allruns (sumssParams.allRuns) ********************************************************/
		else if (!strcmp(parmName, "Allruns"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumssParams.allRuns = YES;
					else
						sumssParams.allRuns = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for allruns (valid arguments are 'yes' and 'no')\n", spacer);
					return (ERROR);
					}
				if (sumssParams.allRuns == YES)
					MrBayesPrint ("%s   Setting sump to print information for each run\n", spacer);
				else
					MrBayesPrint ("%s   Setting sump to print only summary information for all runs\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
        /* set Steptoplot (sumssParams.stepToPlot) *******************************************************/
		else if (!strcmp(parmName, "Steptoplot"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
				if (tempI < 0)
					{
					MrBayesPrint ("%s   Steptoplot must be at least 0\n", spacer);
					return (ERROR);
					}
				else
					{
					sumssParams.stepToPlot = tempI;
					MrBayesPrint ("%s   Setting sumss steptoplot to %d\n", spacer, sumssParams.stepToPlot);
					expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
					}
				}
			else
				return (ERROR);
			}
        /* set Smoothing (sumssParams.smoothing ) *******************************************************/
		else if (!strcmp(parmName, "Smoothing"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
				if (tempI < 0)
					{
					MrBayesPrint ("%s   Smoothing must be at least 0\n", spacer);
					return (ERROR);
					}
				else
					{
					sumssParams.smoothing  = tempI;
					MrBayesPrint ("%s   Setting sumss smoothing to %d\n", spacer, sumssParams.smoothing );
					expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
					}
				}
			else
				return (ERROR);
			}
		/* set Allruns (sumssParams.askForMorePlots) ********************************************************/
		else if (!strcmp(parmName, "Askmore"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumssParams.askForMorePlots = YES;
					else
						sumssParams.askForMorePlots = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for askmore (valid arguments are 'yes' and 'no')\n", spacer);
					return (ERROR);
					}
				if (sumssParams.askForMorePlots == YES)
					MrBayesPrint ("%s   Setting sumss to be interactiva by asking for more plots of burn-in or individual steps.\n", spacer);
				else
					MrBayesPrint ("%s   Setting sumss not to be interactive. It will not ask to print more plots.\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Discardfrac (sumssParams.discardFraction) ************************************************************/
		else if (!strcmp(parmName, "Discardfrac"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);
				if (tempD < 0.00)
					{
					MrBayesPrint ("%s   Discard fraction too low (< 0.00)\n", spacer);
					return (ERROR);
					}
				if (tempD > 1.00)
					{
					MrBayesPrint ("%s   Discard fraction too high (> 1.00)\n", spacer);
					return (ERROR);
					}
                sumssParams.discardFraction = tempD;
				MrBayesPrint ("%s   Setting discard fraction to %.2f\n", spacer, sumssParams.discardFraction);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else 
				{
				return (ERROR);
				}
			}
		else
			return (ERROR);
		}

	return (NO_ERROR);

}





/* ExamineSumpFile: Collect info on the parameter samples in the file */
int ExamineSumpFile (char *fileName, SumpFileInfo *fileInfo, char ***headerNames, int *nHeaders)
{
    char    *sumpTokenP, sumpToken[CMD_STRING_LENGTH], *s=NULL, *headerLine, *t;
    int     i, lineTerm, inSumpComment, lineNum, lastNonDigitLine, numParamLines, allDigitLine,
            lastTokenWasDash, nNumbersOnThisLine, tokenType, burnin, nLines, firstNumCols;
    MrBFlt  tempD;
    FILE    *fp = NULL;


	/* open binary file */
	if ((fp = OpenBinaryFileR(fileName)) == NULL)
        {
        /* test for some simple errors */
        if (strlen(sumpParams.sumpFileName) > 2)
            {
            s = sumpParams.sumpFileName + (int) strlen(sumpParams.sumpFileName) - 2;
		    if (strcmp(s, ".p") == 0)
		        {
		        MrBayesPrint ("%s   It appears that you need to remove '.p' from the 'Filename' parameter\n", spacer);
		        MrBayesPrint ("%s   Also make sure that 'Nruns' is set correctly\n", spacer);
		        return ERROR;
                }
            }
        MrBayesPrint ("%s   Make sure that 'Nruns' is set correctly\n", spacer);
		return ERROR;
        }
    
	/* find out what type of line termination is used */
	lineTerm = LineTermType (fp);
	if (lineTerm != LINETERM_MAC && lineTerm != LINETERM_DOS && lineTerm != LINETERM_UNIX)
		{
		MrBayesPrint ("%s   Unknown line termination\n", spacer);
		goto errorExit;
		}
		
	/* find length of longest line */
	fileInfo->longestLineLength = LongestLine (fp);
	fileInfo->longestLineLength += 10;      /* better safe than sorry; if you fgets with raw longestLineLength, you run into problems */

	/* allocate string long enough to hold a line */
	s = (char *)SafeMalloc((size_t) (2*(fileInfo->longestLineLength + 10) * sizeof(char)));
	if (!s)
		{
		MrBayesPrint ("%s   Problem allocating string for reading sump file\n", spacer);
		goto errorExit;
		}
    headerLine = s + fileInfo->longestLineLength + 10;

	/* close binary file */
	SafeFclose (&fp);
	
	/* open text file */
	if ((fp = OpenTextFileR(fileName)) == NULL)
		goto errorExit;
	
	/* Check file for appropriate blocks. We want to find the last block
	   in the file and start from there. */
	inSumpComment = NO;
	lineNum = lastNonDigitLine = numParamLines = 0;
	while (fgets (s, fileInfo->longestLineLength + 2, fp) != NULL)
        {
		sumpTokenP = &s[0];
		allDigitLine = YES;
		lastTokenWasDash = NO;
		nNumbersOnThisLine = 0;
		do {
			if(GetToken (sumpToken, &tokenType, &sumpTokenP))
                goto errorExit;
            /*printf ("%s (%d)\n", sumpToken, tokenType);*/
			if (IsSame("[", sumpToken) == SAME)
				inSumpComment = YES;
			if (IsSame("]", sumpToken) == SAME)
				inSumpComment = NO;
					
			if (inSumpComment == NO)
				{
				if (tokenType == NUMBER)
					{
					sscanf (sumpToken, "%lf", &tempD);
					if (lastTokenWasDash == YES)
						tempD *= -1.0;
					nNumbersOnThisLine++;
					lastTokenWasDash = NO;
					}
				else if (tokenType == DASH)
					{
					lastTokenWasDash = YES;
					}
				else if (tokenType != UNKNOWN_TOKEN_TYPE)
					{
					allDigitLine = NO;
					lastTokenWasDash = NO;
					}
				}
			} while (*sumpToken);
		lineNum++;
		
		if (allDigitLine == NO)
			{
			lastNonDigitLine = lineNum;
			numParamLines = 0;
			}
		else
			{
			if (nNumbersOnThisLine > 0)
				numParamLines++;
			}
		}
		
	/* Now, check some aspects of the .p file. */
	if (inSumpComment == YES)
		{
		MrBayesPrint ("%s   Unterminated comment in file \"%s\"\n", spacer, fileName);
			goto errorExit;
		}
	if (numParamLines <= 0)
		{
		MrBayesPrint ("%s   No parameters were found in file or there characters not representing a number in last string of the file.\"%s\"\n", spacer, fileName);
			goto errorExit;
		}

    /* calculate burnin */
    if ( chainParams.isSS == YES )
        {
        burnin = 0;
        }
    else
        {
        if (chainParams.relativeBurnin == YES)
            burnin = (int) (chainParams.burninFraction * numParamLines);
        else
            burnin = chainParams.chainBurnIn;
        }
    
	/* check against burnin */
	if (burnin > numParamLines)
		{
		MrBayesPrint ("%s   No parameters can be sampled from file %s as the burnin (%d) exceeds the number of lines in last block (%d)\n",
            spacer, fileName, burnin, numParamLines);
		MrBayesPrint ("%s   Try setting burnin to a number less than %d\n", spacer, numParamLines);
		goto errorExit;
		}

	/* Set some info in fileInfo */
    fileInfo->firstParamLine = lastNonDigitLine + burnin;
    fileInfo->headerLine = lastNonDigitLine;

    /* Calculate and check the number of columns and rows for the file; get header line at the same time */
	(void)fseek(fp, 0L, 0);
	for (lineNum=0; lineNum<lastNonDigitLine; lineNum++)
	    if(fgets (s, fileInfo->longestLineLength + 2, fp)==NULL)
            goto errorExit;
    strcpy(headerLine, s);
    for (; lineNum < lastNonDigitLine+burnin; lineNum++)
	    if(fgets (s, fileInfo->longestLineLength + 2, fp)==NULL)
            goto errorExit;

	inSumpComment = NO;
	nLines = 0;
	firstNumCols = 0;
	while (fgets (s, fileInfo->longestLineLength + 2, fp) != NULL)
        {
		sumpTokenP = &s[0];
		allDigitLine = YES;
		lastTokenWasDash = NO;
		nNumbersOnThisLine = 0;
		do {
			if(GetToken (sumpToken, &tokenType, &sumpTokenP))
                goto errorExit;
			if (IsSame("[", sumpToken) == SAME)
				inSumpComment = YES;
			if (IsSame("]", sumpToken) == SAME)
				inSumpComment = NO;
			if (inSumpComment == NO)
				{
				if (tokenType == NUMBER)
					{
					nNumbersOnThisLine++;
					lastTokenWasDash = NO;
					}
				else if (tokenType == DASH)
					{
					lastTokenWasDash = YES;
					}
				else if (tokenType != UNKNOWN_TOKEN_TYPE)
					{
					allDigitLine = NO;
					lastTokenWasDash = NO;
					}
				}
			} while (*sumpToken);
		lineNum++;
		if (allDigitLine == NO)
			{
			MrBayesPrint ("%s   Found a line with non-digit characters (line %d) in file %s\n", spacer, lineNum, fileName);
			goto errorExit;
			}
		else
			{
			if (nNumbersOnThisLine > 0)
				{
				nLines++;
				if (nLines == 1)
					firstNumCols = nNumbersOnThisLine;
				else
					{
					if (nNumbersOnThisLine != firstNumCols)
						{
						MrBayesPrint ("%s   Number of columns is not even (%d in first line and %d in %d line of file %s)\n", spacer, firstNumCols, nNumbersOnThisLine, lineNum, fileName);
						goto errorExit;
						}
					}
				}
            }
		}
	fileInfo->numRows = nLines;
	fileInfo->numColumns = firstNumCols;

    /* set or check headers */
    if ((*headerNames) == NULL)
        {
        GetHeaders (headerNames, headerLine, nHeaders);
        if (*nHeaders != fileInfo->numColumns)
            {
			MrBayesPrint ("%s   Expected %d headers but found %d headers\n", spacer, fileInfo->numColumns, *nHeaders);
            for (i=0; i<*nHeaders; i++)
                SafeFree ((void **) &((*headerNames)[i]));
            SafeFree ((void **) &(*headerNames));
			*nHeaders=0;
            goto errorExit;
            }
        }
    else
        {
        if (*nHeaders != fileInfo->numColumns)
            {
			MrBayesPrint ("%s   Expected %d columns but found %d columns\n", spacer, *nHeaders, fileInfo->numColumns);
            goto errorExit;
            }
        for (i=0, t=strtok(headerLine,"\t\n\r"); t!=NULL; t=strtok(NULL,"\t\n\r"), i++)
            {
            if (i == *nHeaders)
                {
                MrBayesPrint ("%s   Expected %d headers but found more headers.\n",
                spacer, fileInfo->numColumns);
                goto errorExit;
                }             
            if (strcmp(t,(*headerNames)[i])!=0)
                {
                MrBayesPrint ("%s   Expected header '%s' for column %d but the header for this column was '%s' in file '%s'\n", spacer, (*headerNames)[i], i+1, t, fileName);
                MrBayesPrint ("%s   It could be that some paramiter values are not numbers and the whole string containing \n",spacer); 
                MrBayesPrint ("%s   this wrongly formated paramiter is treated as a header.\n",spacer);
                goto errorExit;
                }
            }
        if (t != NULL)
            {
            MrBayesPrint ("%s   Expected %d headers but found more headers.\n",spacer, fileInfo->numColumns);
            goto errorExit;
            }
        if (i < *nHeaders)
            {
            MrBayesPrint ("%s   Expected header '%s' for column %d but the header for this column was '%s' in file '%s'\n",
                spacer, (*headerNames)[i], i+1, t, fileName);
            goto errorExit;
            }
        }

    free (s);
    fclose(fp);
    return (NO_ERROR);

errorExit:

    free(s);
    fclose(fp);
    return (ERROR);
}





/***************************************************
|
|   FindHeader: Find token in list
|
----------------------------------------------------*/
int FindHeader (char *token, char **headerNames, int nHeaders, int *index)
{
    int         i, match=0, nMatches;

    *index = -1;
    nMatches = 0;
    for (i=0; i<nHeaders; i++)
        {
        if (!strcmp(token,headerNames[i]))
            {
            nMatches++;
            match = i;
            }
        }

    if (nMatches != 1)
        return (ERROR);

    *index = match;
    return (NO_ERROR);
}





/* FreeParameterSamples: Free parameter samples space */
void FreeParameterSamples (ParameterSample *parameterSamples)
{
    if (parameterSamples != NULL)
        {
        free (parameterSamples[0].values[0]);
        free (parameterSamples[0].values);
        free (parameterSamples);
        }
}





/***************************************************
|
|   GetHeaders: Get headers from headerLine and put
|      them in list while updating nHeaders to reflect
|      the number of headers
|
----------------------------------------------------*/
int GetHeaders (char ***headerNames, char *headerLine, int *nHeaders)
{
	char		*s;

    (*nHeaders) = 0;
	for (s=strtok(headerLine," \t\n\r"); s!=NULL; s=strtok(NULL," \t\n\r"))
		{
		if (AddString (headerNames, *nHeaders, s) == ERROR)
			{
			MrBayesPrint ("%s   Error adding header to list of headers \n", spacer, s);
			return ERROR;
			}
		(*nHeaders)++;
        }
				
	return (NO_ERROR);	
}





/* PrintMargLikes: Print marginal likelihoods to screen and to .lstat file */
int PrintMargLikes (char *fileName, char **headerNames, int nHeaders, ParameterSample *parameterSamples, int nRuns, int nSamples)
{
	int		i, j, len, longestHeader, *sampleCounts=NULL;
	char	temp[100];
    Stat    theStats;
    FILE    *fp;
	
	/* calculate longest header */
	longestHeader = 9;	/* length of 'parameter' */
	for (i=0; i<nHeaders; i++)
		{
        strcpy (temp, headerNames[i]);
		len = (int) strlen(temp);
        for (j=0; modelIndicatorParams[j][0]!='\0'; j++)
            if (IsSame (temp,modelIndicatorParams[j]) != DIFFERENT)
                break;
        if (modelIndicatorParams[j][0]!='\0')
            continue;
		if (!strcmp (temp, "Gen"))
			continue;
		if (!strcmp (temp, "LnL") == SAME)
			continue;
        if (!strcmp (temp, "LnPr") == SAME)
            continue;
		if (len > longestHeader)
			longestHeader = len;
		}
	
    /* open output file */
    strncpy (temp, fileName, 90);
    strcat (temp, ".pstat");
    fp = OpenNewMBPrintFile (temp);
    if (!fp)
        return ERROR;

    /* print unique identifier to the output file */
	if (strlen(stamp) > 1)
		fprintf (fp, "[ID: %s]\n", stamp);

    /* allocate and set nSamples */
    sampleCounts = (int *) SafeCalloc (nRuns, sizeof(int));
    if (!sampleCounts)
        {
        fclose(fp);
        return ERROR;
        }
    for (i=0; i<nRuns; i++)
        sampleCounts[i] = nSamples;

    /* print the header rows */
    MrBayesPrint("\n");
	if (sumpParams.HPD == YES)
        MrBayesPrint ("%s   %*c                             95%% HPD Interval\n", spacer, longestHeader, ' ');
    else
        MrBayesPrint ("%s   %*c                            95%% Cred. Interval\n", spacer, longestHeader, ' ');
	MrBayesPrint ("%s   %*c                           --------------------\n", spacer, longestHeader, ' ');

	MrBayesPrint ("%s   Parameter%*c      Mean     Variance     Lower       Upper       Median", spacer, longestHeader-9, ' ');
	if (nRuns > 1)
		MrBayesPrint ("     PSRF+ ");
	MrBayesPrint ("\n");

	MrBayesPrint ("%s   ", spacer);
	for (j=0; j<longestHeader+1; j++)
		MrBayesPrint ("-");
	MrBayesPrint ("-----------------------------------------------------------");
	if (nRuns > 1)
		MrBayesPrint ("----------");
	MrBayesPrint ("\n");
    if (nRuns > 1)
        MrBayesPrintf (fp, "Parameter\tMean\tVariance\tLower\tUpper\tMedian\tPSRF\n");
    else
        MrBayesPrintf (fp, "Parameter\tMean\tVariance\tLower\tUpper\tMedian\n");

	/* print table values */
	for (i=0; i<nHeaders; i++)
		{
        strcpy (temp, headerNames[i]);
		len = (int) strlen(temp);
        for (j=0; modelIndicatorParams[j][0]!='\0'; j++)
            if (IsSame (temp,modelIndicatorParams[j]) != DIFFERENT)
                break;
		if (IsSame (temp, "Gen") == SAME)
			continue;
		if (IsSame (temp, "LnL") == SAME)
			continue;
        if (!strcmp (temp, "LnPr") == SAME)
            continue;

		GetSummary (parameterSamples[i].values, nRuns, sampleCounts, &theStats, sumpParams.HPD);
		
		MrBayesPrint ("%s   %-*s ", spacer, longestHeader, temp);
		MrBayesPrint ("%10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf", theStats.mean, theStats.var, theStats.lower, theStats.upper, theStats.median);
		MrBayesPrintf (fp, "%s", temp);
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.mean));
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.var));
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.lower));
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.upper));
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.median));
		if (nRuns > 1)
			{
			if (theStats.PSRF < 0.0)
                {
				MrBayesPrint ("     NA   ");
                MrBayesPrintf (fp, "NA");
                }
			else
                {
				MrBayesPrint ("  %7.3lf", theStats.PSRF);
                MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.PSRF));
                }
			}
		MrBayesPrint ("\n");
		MrBayesPrintf (fp, "\n");
		}
	MrBayesPrint ("%s   ", spacer);
	for (j=0; j<longestHeader+1; j++)
		MrBayesPrint ("-");
	MrBayesPrint ("-----------------------------------------------------------");
	if (nRuns > 1)
		MrBayesPrint ("----------");
	MrBayesPrint ("\n");
	if (nRuns > 1)
		{
		MrBayesPrint ("%s   + Convergence diagnostic (PSRF = Potential Scale Reduction Factor; Gelman\n", spacer);
		MrBayesPrint ("%s     and Rubin, 1992) should approach 1.0 as runs converge.\n", spacer);
		}

    fclose (fp);
    free (sampleCounts);

    return (NO_ERROR);
}





/* PrintModelStats: Print model stats to screen and to .mstat file */
int PrintModelStats (char *fileName, char **headerNames, int nHeaders, ParameterSample *parameterSamples, int nRuns, int nSamples)
{
	int         i, j, j1, j2, k, longestName, nElements, *modelCounts=NULL;
	MrBFlt	    f, *prob=NULL, *sum=NULL, *ssq=NULL, *min=NULL, *max=NULL, *stddev=NULL;
	char	    temp[100];
    FILE        *fp;
    ModelProb   *elem = NULL;

    /* nHeaders - is a convenient synonym for number of column headers */

    /* check if we have any model indicator variables and also check for longest header */
    k = 0;
    longestName = 0;
    for (i=0; i<nHeaders; i++)
        {
        for (j=0; strcmp(modelIndicatorParams[j],"")!=0; j++)
            {
            if (IsSame (headerNames[i], modelIndicatorParams[j]) != DIFFERENT)
                {
                k++;
                for (j1=0; strcmp(modelElementNames[j][j1],"")!=0; j1++)
                    {
                    j2 = (int)(strlen(headerNames[i]) + 2 + strlen(modelElementNames[j][j1]));
                    if (j2 > longestName)
                        longestName = j2;
                    }
                break;
                }
            }
        }

    /* return if nothing to do */
    if (k==0)
        return NO_ERROR;

    /* open output file */
	MrBayesPrint ("%s   Model probabilities above %1.3lf\n", spacer, sumpParams.minProb);
    MrBayesPrint ("%s   Estimates saved to file \"%s.mstat\".\n", spacer, sumpParams.sumpOutfile);
    strncpy (temp,fileName,90);
    strcat (temp, ".mstat");
    fp = OpenNewMBPrintFile(temp);
    if (!fp)
        return ERROR;
    MrBayesPrint ("\n");

    /* print unique identifier to the output file */
	if (strlen(stamp) > 1)
		fprintf (fp, "[ID: %s]\n", stamp);

    /* print header */
	MrBayesPrintf (fp, "\n\n");
    if (nRuns == 1)
        {
        MrBayesPrint ("%s        %*c        Posterior\n", spacer, longestName-5, ' ');
        MrBayesPrint ("%s   Model%*c       Probability\n", spacer, longestName-5, ' ');
        MrBayesPrint ("%s   -----", spacer);
        for (i=0; i<longestName-5; i++)
            MrBayesPrint ("-");
        MrBayesPrint ("------------------\n");
        MrBayesPrintf (fp, "Model\tProbability\n");
        }
    else
        {
		MrBayesPrint ("%s        %*c        Posterior      Standard         Min.           Max.   \n", spacer, longestName-5, ' ');
		MrBayesPrint ("%s   Model%*c       Probability     Deviation     Probability    Probability\n", spacer, longestName-5, ' ');
        MrBayesPrint ("%s   -----", spacer);
        for (i=0; i<longestName-5; i++)
            MrBayesPrint ("-");
        MrBayesPrint ("---------------------------------------------------------------\n");
        MrBayesPrintf (fp, "Model\tProbability\tStd_dev\tMin_prob\tMax_prob\n");
        }

    /* calculate and print values */
    for (i=0; i<nHeaders; i++)
		{
        for (j=0; modelIndicatorParams[j][0]!='\0'; j++)
            if (IsSame (headerNames[i], modelIndicatorParams[j]) != DIFFERENT)
                break;
        if (modelIndicatorParams[j][0] == '\0')
            continue;

        for (nElements=0; modelElementNames[j][nElements][0]!='\0'; nElements++)
            ;
        
        modelCounts = (int *) SafeCalloc (nElements, sizeof(int));
        if (!modelCounts)
            {
            fclose(fp);
            return ERROR;
            }
        prob = (MrBFlt *) SafeCalloc (6*nElements, sizeof(MrBFlt));
        if (!prob)
            {
            free (modelCounts);
            fclose (fp);
            return ERROR;
            }
        sum    = prob + nElements;
        ssq    = prob + 2*nElements;
        stddev = prob + 3*nElements;
        min    = prob + 4*nElements;
        max    = prob + 5*nElements;

        for (j1=0; j1<nElements; j1++)
            min[j1] = 1.0;

        for (j1=0; j1<nRuns; j1++)
            {
            for (j2=0; j2<nElements; j2++)
                modelCounts[j2] = 0;
            for (j2=0; j2<nSamples; j2++)
                modelCounts[(int)(parameterSamples[i].values[j1][j2] + 0.1)]++;
            for (j2=0; j2<nElements; j2++)
                {
				f = (MrBFlt) modelCounts[j2] / (MrBFlt) nSamples;
                sum[j2] += f;
                ssq[j2] += f*f;
                if (f<min[j2])
                    min[j2] = f;
                if (f > max[j2])
                    max[j2] = f;
                }
            }

		for (j1=0; j1<nElements; j1++)
			{
            prob[j1] = sum[j1] / (MrBFlt) nRuns;
			f = ssq[j1] - (sum[j1] * sum[j1] / (MrBFlt) nRuns);
			f /= (nRuns - 1);
			if (f <= 0.0)
				stddev[j1] = 0.0;
			else
				stddev[j1] = sqrt (f);
			}

        elem = (ModelProb *) SafeCalloc (nElements, sizeof(ModelProb));
        for (j1=0; j1<nElements; j1++)
            {
            elem[j1].index = j1;
            elem[j1].prob = prob[j1];
            }

        /* sort in terms of decreasing probabilities */
        qsort((void *) elem, (size_t) nElements, (size_t) sizeof(ModelProb), CompareModelProbs);

        for (j1=0; j1<nElements; j1++)
			{
        	if (elem[j1].prob <= sumpParams.minProb)
                break;

            if (nRuns == 1)
        	    {
        	    sprintf (temp, "%s[%s]", headerNames[i], modelElementNames[j][elem[j1].index]);
			    MrBayesPrint ("%s   %-*s          %1.3lf\n", spacer, longestName, temp, prob[elem[j1].index]);
        	    MrBayesPrintf (fp, "%s\t%s\n", temp, MbPrintNum(prob[elem[j1].index])); 
        	    }
    	    else /* if (nRuns > 1) */
        	    {
        	    sprintf (temp, "%s[%s]", headerNames[i], modelElementNames[j][elem[j1].index]);
			    MrBayesPrint ("%s   %-*s          %1.3lf          %1.3lf          %1.3lf          %1.3lf\n", 
            	    spacer, longestName, temp, prob[elem[j1].index], stddev[elem[j1].index], min[elem[j1].index], max[elem[j1].index]);
        	    MrBayesPrintf (fp, "%s", temp);
        	    MrBayesPrintf (fp, "\t%s", MbPrintNum(prob[elem[j1].index]));
        	    MrBayesPrintf (fp, "\t%s", MbPrintNum(stddev[elem[j1].index]));
        	    MrBayesPrintf (fp, "\t%s", MbPrintNum(min[elem[j1].index]));
        	    MrBayesPrintf (fp, "\t%s", MbPrintNum(max[elem[j1].index]));
        	    MrBayesPrintf (fp, "\n");
        	    }
			}
        free(elem);
        elem = NULL;
        free(modelCounts);
        modelCounts = NULL;
        free (prob);
        prob = NULL;
		}

	/* print footer */
    if (nRuns == 1)
        {
        MrBayesPrint ("%s   -----", spacer);
        for (i=0; i<longestName-5; i++)
            MrBayesPrint ("-");
        MrBayesPrint ("------------------\n\n");
        }
    else
        {
        MrBayesPrint ("%s   -----", spacer);
        for (i=0; i<longestName-5; i++)
            MrBayesPrint ("-");
        MrBayesPrint ("---------------------------------------------------------------\n\n");
        }

    /* close output file */
    fclose (fp);

	return (NO_ERROR);
}





/* PrintOverlayPlot: Print overlay x-y plot of log likelihood vs. generation for several runs */
int PrintOverlayPlot (MrBFlt **xVals, MrBFlt **yVals, int nRuns,  int startingFrom, int nSamples)
{
	int		i, j, k, k2, n, screenHeight, screenWidth, numY[60], width;
	char	plotSymbol[15][60];
	MrBFlt	x, y, minX, maxX, minY, maxY, meanY[60];
	
	if (nRuns == 2)
		MrBayesPrint ("\n%s   Overlay plot for both runs:\n", spacer);
	else
		MrBayesPrint ("\n%s   Overlay plot for all %d runs:\n", spacer, sumpParams.numRuns);
	if (nRuns > 9)
		MrBayesPrint ("%s   (1 = Run number 1; 2 = Run number 2 etc.; x = Run number 10 or above; * = Several runs)\n", spacer);
	else if (nRuns > 2)
		MrBayesPrint ("%s   (1 = Run number 1; 2 = Run number 2 etc.; * = Several runs)\n", spacer);
	else
		MrBayesPrint ("%s   (1 = Run number 1; 2 = Run number 2; * = Both runs)\n", spacer);

	/* print overlay x-y plot of log likelihood vs. generation for all runs */
	screenWidth = 60; /* don't change this without changing numY, meanY, and plotSymbol declared above */
	screenHeight = 15;

	/* find minX, minY, maxX, and maxY over all runs */
	minX = minY = 1000000000.0;
	maxX = maxY = -1000000000.0;
	for (n=0; n<nRuns; n++)
		{
		for (i=startingFrom; i<startingFrom+nSamples; i++)
			{
			x = xVals[n][i];
			if (x < minX)
				minX = x;
			if (x > maxX)
				maxX = x;
			}
		}
	for (n=0; n<nRuns; n++)
		{
		y = 0.0;
		j = 0;
		k2 = 0;
		for (i=startingFrom; i<startingFrom+nSamples; i++)
			{
			x = xVals[n][i];
			k = (int) (((x - minX) / (maxX - minX)) * screenWidth);
			if (k <= 0)
				k = 0;
			if (k >= screenWidth)
				k = screenWidth - 1;
			if (k == j)
				{
				y += yVals[n][i];
				k2 ++;
				}
			else
				{
				y /= k2;
				if (y < minY)
					minY = y;
				if (y > maxY)
					maxY = y;
				k2 = 1;
				y = yVals[n][i];
				j++;
				}
			}
		if (k2 > 0)
			{
			y /= k2;
			if (y < minY)
				minY = y;
			if (y > maxY)
				maxY = y;
			}
		}
	
	/* initialize the plot symbols */
	for (i=0; i<screenHeight; i++)
		for (j=0; j<screenWidth; j++)
			plotSymbol[i][j] = ' ';

	/* assemble the plot symbols */
	for (n=0; n<nRuns; n++)
		{
		/* find the plot points for this run */
		for (i=0; i<screenWidth; i++)
			{
			numY[i] = 0;
			meanY[i] = 0.0;
			}
		for (i=startingFrom; i<startingFrom+nSamples; i++)
			{
			x = xVals[n][i];
			y = yVals[n][i];
			k = (int)(((x - minX) / (maxX - minX)) * screenWidth);
			if (k >= screenWidth)
				k = screenWidth - 1;
			if (k <= 0)
				k = 0;
			meanY[k] += y;
			numY[k]++;
			}

		/* transfer these points to the overlay */
		for (i=0; i<screenWidth; i++)
			{
			if (numY[i] > 0)
				{
				k = (int) ((((meanY[i] / numY[i]) - minY)/ (maxY - minY)) * screenHeight);
				if (k < 0)
					k = 0;
				else if (k >= screenHeight)
					k = screenHeight - 1;
				if (plotSymbol[k][i] == ' ')
					{
					if (n <= 8)
						plotSymbol[k][i] = '1' + n;
					else
						plotSymbol[k][i] = 'x';
					}
				else
					plotSymbol[k][i] = '*';
				}
			}
		} /* next run */

	/* now print the overlay plot */
	MrBayesPrint ("\n%s   +", spacer);
	for (i=0; i<screenWidth; i++)
		MrBayesPrint ("-");
	MrBayesPrint ("+ %1.2lf\n", maxY);
	for (i=screenHeight-1; i>=0; i--)
		{
		MrBayesPrint ("%s   |", spacer);
		for (j=0; j<screenWidth; j++)
			{
			MrBayesPrint ("%c", plotSymbol[i][j]);
			}
		MrBayesPrint ("|\n");
		}
	MrBayesPrint ("%s   +", spacer);
	for (i=0; i<screenWidth; i++)
		{
		if (i % (screenWidth/10) == 0 && i != 0)
			MrBayesPrint ("+");
		else
			MrBayesPrint ("-");
		}
	MrBayesPrint ("+ %1.2lf\n", minY);
	MrBayesPrint ("%s   ^", spacer);
	for (i=0; i<screenWidth; i++)
		MrBayesPrint (" ");
	MrBayesPrint ("^\n");
	MrBayesPrint ("%s   %1.0lf", spacer, minX);
    if((int)minX>0)
        width=(int)(log10(minX));
    else if((int)minX==0)
        width=1;
    else
        width=(int)(log10(-minX))+1;
	for (i=0; i<screenWidth-width; i++)
		MrBayesPrint (" ");
	MrBayesPrint ("%1.0lf\n\n", maxX);

	return (NO_ERROR);
}





/* PrintParamStats: Print parameter table (not model indicator params) to screen and .pstat file */
int PrintParamStats (char *fileName, char **headerNames, int nHeaders, ParameterSample *parameterSamples, int nRuns, int nSamples)
{
	int		i, j, len, longestHeader, *sampleCounts=NULL;
    static char	*temp=NULL;
    char    tempf[100];
    Stat    theStats;
    FILE    *fp;
	
	/* calculate longest header */
	longestHeader = 9;	/* length of 'parameter' */
	for (i=0; i<nHeaders; i++)
		{
        SafeStrcpy (&temp, headerNames[i]);
		len = (int) strlen(temp);
        for (j=0; modelIndicatorParams[j][0]!='\0'; j++)
            if (IsSame (temp,modelIndicatorParams[j]) != DIFFERENT)
                break;
        if (modelIndicatorParams[j][0]!='\0')
            continue;
		if (!strcmp (temp, "Gen"))
			continue;
		if (!strcmp (temp, "LnL") == SAME)
			continue;
        if (!strcmp (temp, "LnPr") == SAME)
            continue;
        if (len > longestHeader)
			longestHeader = len;
		}
	
    /* open output file */
    strncpy (tempf, fileName, 90);
    strcat (tempf, ".pstat");
    fp = OpenNewMBPrintFile (tempf);
    if (!fp)
        return ERROR;

    /* print unique identifier to the output file */
	if (strlen(stamp) > 1)
		fprintf (fp, "[ID: %s]\n", stamp);

    /* allocate and set nSamples */
    sampleCounts = (int *) SafeCalloc (nRuns, sizeof(int));
    if (!sampleCounts)
        {
        fclose(fp);
        return ERROR;
        }
    for (i=0; i<nRuns; i++)
        sampleCounts[i] = nSamples;

    /* print the header rows */
    MrBayesPrint("\n");
	if (sumpParams.HPD == YES)
        MrBayesPrint ("%s   %*c                             95%% HPD Interval\n", spacer, longestHeader, ' ');
    else
        MrBayesPrint ("%s   %*c                            95%% Cred. Interval\n", spacer, longestHeader, ' ');
	MrBayesPrint ("%s   %*c                           --------------------\n", spacer, longestHeader, ' ');

	if (nRuns > 1)
        MrBayesPrint ("%s   Parameter%*c     Mean      Variance     Lower       Upper       Median    min ESS*  avg ESS    PSRF+ ", spacer, longestHeader-9, ' ');
    else
        MrBayesPrint ("%s   Parameter%*c     Mean      Variance     Lower       Upper       Median     ESS*", spacer, longestHeader-9, ' ');
	MrBayesPrint ("\n");

	MrBayesPrint ("%s   ", spacer);
	for (j=0; j<longestHeader+1; j++)
		MrBayesPrint ("-");
	MrBayesPrint ("---------------------------------------------------------------------");
	if (nRuns > 1)
		MrBayesPrint ("-------------------");
	MrBayesPrint ("\n");
    if (nRuns > 1)
        MrBayesPrintf (fp, "Parameter\tMean\tVariance\tLower\tUpper\tMedian\tminESS\tavgESS\tPSRF\n");
    else
        MrBayesPrintf (fp, "Parameter\tMean\tVariance\tLower\tUpper\tMedian\tESS\n");

	/* print table values */
	for (i=0; i<nHeaders; i++)
		{
        SafeStrcpy(&temp, headerNames[i]);
		len = (int) strlen(temp);
        for (j=0; modelIndicatorParams[j][0]!='\0'; j++)
            if (IsSame (temp,modelIndicatorParams[j]) != DIFFERENT)
                break;
        if (modelIndicatorParams[j][0]!='\0')
            continue;
		if (IsSame (temp, "Gen") == SAME)
			continue;
		if (IsSame (temp, "LnL") == SAME)
			continue;
        if (!strcmp (temp, "LnPr") == SAME)
            continue;

		GetSummary (parameterSamples[i].values, nRuns, sampleCounts, &theStats, sumpParams.HPD);
		
		MrBayesPrint ("%s   %-*s ", spacer, longestHeader, temp);
		MrBayesPrint ("%10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf", theStats.mean, theStats.var, theStats.lower, theStats.upper, theStats.median);
		MrBayesPrintf (fp, "%s", temp);
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.mean));
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.var));
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.lower));
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.upper));
		MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.median));

        if(theStats.minESS == theStats.minESS)
            {
            MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.minESS));
            MrBayesPrint ("  %8.2lf", theStats.minESS);
            }
        else
            {
            MrBayesPrint ("       NA ");
            MrBayesPrintf (fp, "NA");
            }
		if (nRuns > 1)
			{
            if(theStats.minESS == theStats.minESS)
                {
                MrBayesPrint ("  %8.2lf", theStats.avrESS);
                MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.avrESS));
                }
            else
                {
                MrBayesPrint ("       NA ");
                MrBayesPrintf (fp, "NA");
                }
			if (theStats.PSRF < 0.0)
                {
				MrBayesPrint ("     NA   ");
                MrBayesPrintf (fp, "NA");
                }
			else
                {
				MrBayesPrint ("  %7.3lf", theStats.PSRF);
                MrBayesPrintf (fp, "\t%s", MbPrintNum(theStats.PSRF));
                }
			}
		MrBayesPrint ("\n");
		MrBayesPrintf (fp, "\n");
		}
	MrBayesPrint ("%s   ", spacer);
	for (j=0; j<longestHeader+1; j++)
		MrBayesPrint ("-");
	MrBayesPrint ("---------------------------------------------------------------------");
	if (nRuns > 1)
		MrBayesPrint ("-------------------");
	MrBayesPrint ("\n");
	if (nRuns > 1)
		{
		MrBayesPrint ("%s   * Convergence diagnostic (ESS = Estimated Sample Size); min and avg values\n", spacer);
        MrBayesPrint ("%s     correspond to minimal and average ESS among runs. \n", spacer); 
        MrBayesPrint ("%s     ESS value below 100 may indicate that the parameter is undersampled. \n", spacer);
        MrBayesPrint ("%s   + Convergence diagnostic (PSRF = Potential Scale Reduction Factor; Gelman\n", spacer);
		MrBayesPrint ("%s     and Rubin, 1992) should approach 1.0 as runs converge.\n", spacer);
		}
    else
        {
        MrBayesPrint ("%s   * Convergence diagnostic (ESS = Estimated Sample Size); ESS value \n", spacer);
        MrBayesPrint ("%s     below 100 may indicate that the parameter is undersampled. \n", spacer);
        }
	MrBayesPrint ("\n\n");

    fclose (fp);
    free (sampleCounts);
    SafeFree ((void **)&temp);

    return (NO_ERROR);
}





/* PrintPlot: Print x-y plot of log likelihood vs. generation */
int PrintPlot (MrBFlt *xVals, MrBFlt *yVals, int numVals)
{
	int		i, j, k, numY[60], screenWidth, screenHeight;
	MrBFlt	x, y, minX, maxX, minY, maxY, meanY[60], diff;
					
	/* print x-y plot of log likelihood vs. generation */
	screenWidth = 60; /* don't change this without changing numY and meanY, declared above */
	screenHeight = 15;

	/* find minX and maxX */
	minX = xVals[0];
	maxX = xVals[0];
	for (i=0; i<numVals; i++)
		{
		x = xVals[i];
		if (x < minX)
			minX = x;
		if (x > maxX)
			maxX = x;
		}

	/* collect Y data */
	for (i=0; i<screenWidth; i++)
		{
		numY[i] = 0;
		meanY[i] = 0.0;
		}
	for (i=0; i<numVals; i++)
		{
		x = xVals[i];
		y = yVals[i];
		k = (int)(((x - minX) / (maxX - minX)) * screenWidth);
		if (k >= screenWidth)
			k = screenWidth - 1;
		if (k < 0)
			k = 0;
		meanY[k] += y;
		numY[k]++;
		}

	/* find minY and maxY */
	minY = maxY = meanY[0] / numY[0];
	for (i=0; i<screenWidth; i++)
		{
        if( meanY[i] == 0) /* with some compilers if( NaN < 1 ) is equal true !!! so we realy need this check*/
            continue;
		meanY[i] /= numY[i];
		if (meanY[i] < minY)
			minY = meanY[i];
		if (meanY[i] > maxY)
			maxY = meanY[i];
		}

    /* find difference */
    diff = maxY - minY;

    /* print plot */
	MrBayesPrint ("\n   +");
	for (i=0; i<screenWidth; i++)
		MrBayesPrint ("-");
	MrBayesPrint ("+ %1.3lf\n", maxY);
	for (j=screenHeight-1; j>=0; j--)
		{
		MrBayesPrint ("   |");
		for (i=0; i<screenWidth; i++)
			{
			if (numY[i] > 0)
				{
				if ((meanY[i] > ((diff/screenHeight)*j)+minY && meanY[i] <= ((diff/screenHeight)*(j+1))+minY) ||
                    (j == 0 && meanY[i] <= minY))
					MrBayesPrint ("*");
				else
					MrBayesPrint (" ");
				}
			else
				{
				MrBayesPrint (" ");
				}
			}
		MrBayesPrint ("|\n");
		}
	MrBayesPrint ("   +");
	for (i=0; i<screenWidth; i++)
		{
		if (i % (screenWidth/10) == 0 && i != 0)
			MrBayesPrint ("+");
		else
			MrBayesPrint ("-");
		}
	MrBayesPrint ("+ %1.3lf\n", minY);
	MrBayesPrint ("   ^");
	for (i=0; i<screenWidth; i++)
		MrBayesPrint (" ");
	MrBayesPrint ("^\n");
	MrBayesPrint ("   %1.0lf", minX);
    if (minX == 0)
        j = 1;
    else
        j = (int)(log10(minX)) + 1;
	for (i=0; i<screenWidth-j; i++)
		MrBayesPrint (" ");
	MrBayesPrint ("%1.0lf\n\n", maxX);

	return (NO_ERROR);
}





void PrintPlotHeader (void)
{
	MrBayesPrint ("\n");
	if (sumpParams.numRuns > 1)
		{
		MrBayesPrint ("%s   Below are rough plots of the generation (x-axis) versus the log   \n", spacer);
		MrBayesPrint ("%s   probability of observing the data (y-axis). You can use these     \n", spacer);
		MrBayesPrint ("%s   graphs to determine what the burn in for your analysis should be. \n", spacer);
		MrBayesPrint ("%s   When the log probability starts to plateau you may be at station- \n", spacer);
		MrBayesPrint ("%s   arity. Sample trees and parameters after the log probability      \n", spacer);
		MrBayesPrint ("%s   plateaus. Of course, this is not a guarantee that you are at sta- \n", spacer);
		MrBayesPrint ("%s   tionarity. Also examine the convergence diagnostics provided by   \n", spacer);
		MrBayesPrint ("%s   the 'sump' and 'sumt' commands for all the parameters in your     \n", spacer);
		MrBayesPrint ("%s   model. Remember that the burn in is the number of samples to dis- \n", spacer);
		MrBayesPrint ("%s   card. There are a total of ngen / samplefreq samples taken during \n", spacer);
		MrBayesPrint ("%s   a MCMC analysis.                                                  \n", spacer);
		}
	else
		{
		MrBayesPrint ("%s   Below is a rough plot of the generation (x-axis) versus the log   \n", spacer);
		MrBayesPrint ("%s   probability of observing the data (y-axis). You can use this      \n", spacer);
		MrBayesPrint ("%s   graph to determine what the burn in for your analysis should be.  \n", spacer);
		MrBayesPrint ("%s   When the log probability starts to plateau you may be at station- \n", spacer);
		MrBayesPrint ("%s   arity. Sample trees and parameters after the log probability      \n", spacer);
		MrBayesPrint ("%s   plateaus. Of course, this is not a guarantee that you are at sta- \n", spacer);
		MrBayesPrint ("%s   analysis should be. When the log probability starts to plateau    \n", spacer);
		MrBayesPrint ("%s   tionarity. When possible, run multiple analyses starting from dif-\n", spacer);
		MrBayesPrint ("%s   ferent random trees; if the inferences you make for independent   \n", spacer);
		MrBayesPrint ("%s   analyses are the same, this is reasonable evidence that the chains\n", spacer);
		MrBayesPrint ("%s   have converged. You can use MrBayes to run several independent    \n", spacer);
		MrBayesPrint ("%s   analyses simultaneously. During such a run, MrBayes will monitor  \n", spacer);
		MrBayesPrint ("%s   the convergence of topologies. After the run has been completed,  \n", spacer);
		MrBayesPrint ("%s   the 'sumt' and 'sump' functions will provide additional conver-   \n", spacer);
		MrBayesPrint ("%s   gence diagnostics for all the parameters in your model. Remember  \n", spacer);
		MrBayesPrint ("%s   that the burn in is the number of samples to discard. There are   \n", spacer);
		MrBayesPrint ("%s   a total of ngen / samplefreq samples taken during a MCMC analysis.\n", spacer);
		}
}





/* ReadParamSamples: Read parameter samples from .p file */
int ReadParamSamples (char *fileName, SumpFileInfo *fileInfo, ParameterSample *parameterSamples, int runNo)
{
    char    sumpToken[CMD_STRING_LENGTH], *s=NULL, *p;
    int     inSumpComment, lineNum, numLinesRead, numLinesToRead, column, lastTokenWasDash,
            tokenType;
    MrBFlt  tempD;
    FILE    *fp;

    /* open file */
    if ((fp = OpenTextFileR (fileName)) == NULL)
        return ERROR;

    /* allocate space for reading lines */
    s = (char *) SafeCalloc (fileInfo->longestLineLength + 10, sizeof(char));

	/* fast forward to beginning of last unburned parameter line. */
	for (lineNum=0; lineNum<fileInfo->firstParamLine; lineNum++)
	  if (fgets (s, fileInfo->longestLineLength + 5, fp) == 0)
          goto errorExit;

    /* parse file, line-by-line. We are only parsing lines that have digits that should be read. */
	inSumpComment = NO;
	numLinesToRead = fileInfo->numRows;
	numLinesRead = 0;
	while (fgets (s, fileInfo->longestLineLength + 1, fp) != NULL)
		{
		lastTokenWasDash = NO;
        column = 0;
        p = s;
		do {
			if(GetToken (sumpToken, &tokenType, &p))
                goto errorExit;
            if (IsSame("[", sumpToken) == SAME)
				inSumpComment = YES;
			if (IsSame("]", sumpToken) == SAME)
				inSumpComment = NO;
			if (inSumpComment == NO)
				{
				if (tokenType == NUMBER)
					{
					/* read the number */
					if (column >= fileInfo->numColumns)
						{
						MrBayesPrint ("%s   Too many values read on line %d of file %s\n", spacer, lineNum, fileName);
						goto errorExit;
						}
					sscanf (sumpToken, "%lf", &tempD);
					if (lastTokenWasDash == YES)
						tempD *= -1.0;
					parameterSamples[column].values[runNo][numLinesRead] = tempD;
					column++;
					lastTokenWasDash = NO;
					}
				else if (tokenType == DASH)
					{
					lastTokenWasDash = YES;
					}
				else if (tokenType != UNKNOWN_TOKEN_TYPE)
					{
					/* we have a problem */
					MrBayesPrint ("%s   Line %d of file %s has non-digit characters\n", spacer, lineNum, fileName);
					goto errorExit;
					}
				}
			} while (*sumpToken);

        lineNum++;
        if (column == fileInfo->numColumns)
			numLinesRead++;
        else if (column != 0)
            {
            MrBayesPrint ("%s   Too few values on line %d of file %s\n", spacer, lineNum, fileName);
			goto errorExit;
            }
        }

		
	/* Check how many parameter line was read in. */
	if (numLinesRead != numLinesToRead)
		{
		MrBayesPrint ("%s   Unable to read all lines that should contain parameter samples\n", spacer);
		goto errorExit;
		}

    fclose (fp);
    free (s);

    return (NO_ERROR);

errorExit:

    fclose (fp);
    free (s);

    return ERROR;
}

