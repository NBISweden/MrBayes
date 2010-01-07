/*
 *  MrBayes 3.1.2
 *
 *  copyright 2002-2005
 *
 *  John P. Huelsenbeck
 *  Section of Ecology, Behavior and Evolution
 *  Division of Biological Sciences
 *  University of California, San Diego
 *  La Jolla, CA 92093-0116
 *
 *  johnh@biomail.ucsd.edu
 *
 *	Fredrik Ronquist
 *  Paul van der Mark
 *  School of Computational Science
 *  Florida State University
 *  Tallahassee, FL 32306-4120
 *
 *  ronquist@scs.fsu.edu
 *  paulvdm@scs.fsu.edu
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
 * Paul 20070501:
 * Applied a pach from Julian Catchen <catchen@cs.uoregon.edu> for a bug
 * concerning non-matching headers.
 * "I tracked down the problem and fixed it in the attached patch. It was a
 * simple problem where a string parsing variable (headerNames) would be
 * reallocated when reading the second (or any subsequent file) file that
 * had lines longer than those in the first file. Problem was, it would
 * reallocate the variable and then clear it, causing any further header
 * comparisons to fail with the above message."
 */
/* id-string for ident, do not edit: cvs will update this string */
const char sumpID[]="$Id: sump.c,v 3.37 2009/01/06 21:40:10 ronquist Exp $";

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
#include "mcmc.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif



/* local prototypes */
int      HarmonicArithmeticMean (MrBFlt *paramVals, int nRows, int nCols, int col, MrBFlt *mean, MrBFlt *harm_mean);
void     MeanVariance(MrBFlt *paramVals, int nRows, int nCols, int col, MrBFlt *m, MrBFlt *v);
MrBFlt   PotentialScaleReduction (MrBFlt *paramVals, int nRows, int nCols, int nRuns, int col);
int		 PrintAATable (FILE *fp, MrBFlt *parmVals, int nRows, int nCols, int nRuns);
int		 PrintOverlayPlot (FILE *fp);
int		 PrintParamTable (FILE *fp, MrBFlt *parmVals, int nRows, int nCols, int nRuns);
int      PrintPlot (FILE *fp, MrBFlt *paramVals);
void	 PrintPlotHeader (FILE *fp);
void     Sort2 (MrBFlt *item, int left, int right);
int 	 SortColumn (MrBFlt *paramValues, int nRows, int nCols, int column);
int      SortParameters (MrBFlt *paramValues);

/* global parameters declared in this file */
char		*sumpTokenP, sumpToken[CMD_STRING_LENGTH];

/* local (but global to this file) parameters */
int			numColumns, numRows;
MrBFlt		*parameterValues;
char		*headerNames;



int DoSump (void)

{

	int			i, j, n, lineTerm=LINETERM_UNIX, longestLineLength=0, tokenType, lineNum, lastTokenWasDash,
				inSumpComment, allDigitLine, nNumbersOnThisLine, lastNonDigitLine,
				numParamLines, numLinesToRead, numLinesRead, firstNumCols=0, nLines,
				nHeaders=0, len, longestHeader, firstFileNumParamLines=0, firstFileNumRows=0, firstFileNumColumns=0,
				unreliable, oneUnreliable, burnin=0;
	MrBFlt		tempD, *parameterPointer, mean, harm_mean;
	char		*s=NULL, *s1, *headerLine=NULL, temp[100];
	FILE		*fp=NULL, *fpOut=NULL;
	
#	if defined (MPI_ENABLED)
	int			ierror;
#	endif
	
#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
#	endif

	/* set file pointers to NULL */
	fp = fpOut = NULL;

	/* check if there is anything to do */
	if (sumpParams.margLike == NO && sumpParams.plot == NO && sumpParams.table == NO)
		{
		MrBayesPrint ("%s   Nothing to do for sump. Try setting Plot, Marglike or Table to yes.\n", spacer);
		return NO_ERROR;
		}

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

	if (sumpParams.printToFile == NO)
		MrBayesPrint ("%s   Writing output to screen but not to file ('Printtofile = No')\n", spacer);
	else
		MrBayesPrint ("%s   Writing output to file %s\n", spacer, sumpParams.sumpOutfile);

	/* open output text file */
	if (sumpParams.printToFile == YES)
		{
		if ((fpOut = OpenNewMBPrintFile (sumpParams.sumpOutfile)) == NULL)
			goto errorExit;
		/* print unique identifier to the output file */
		len = (int) strlen (stamp);
		if (len <= 1)
			fprintf (fpOut, "[ID: None Available]\n");
		else
			fprintf (fpOut, "[ID: %s]\n", stamp);
		}

	/* open input file(s) */
	for (n=0; n<sumpParams.numRuns; n++)
		{
		  DEBUG("reading file %d\n",n);
		/* open binary file */
		if (sumpParams.numRuns == 1)
			sprintf (temp, "%s.p", sumpParams.sumpFileName);
		else
			sprintf (temp, "%s.run%d.p", sumpParams.sumpFileName, n+1);
		DEBUG("strlen temp = %u\n",strlen(temp));

		if ((fp = OpenBinaryFileR(temp)) == NULL)
			{
			s = sumpParams.sumpFileName + (int) strlen(sumpParams.sumpFileName) - 2;
			if (strcmp(s, ".p") == 0)
				{
				MrBayesPrint ("%s   It appears that you need to remove '.p' from the 'Filename' parameter\n", spacer);
				MrBayesPrint ("%s   Also make sure that 'Nruns' is set correctly\n", spacer);
				}
			else
				MrBayesPrint ("%s   Make sure that 'Nruns' is set correctly\n", spacer);
			goto errorExit;
			}
		
		/* find out what type of line termination is used */
		if (n == 0)
			{
			lineTerm = LineTermType (fp);
			if (lineTerm == LINETERM_MAC)
				MrBayesPrint ("%s   Macintosh line termination\n", spacer);
			else if (lineTerm == LINETERM_DOS)
				MrBayesPrint ("%s   DOS line termination\n", spacer);
			else if (lineTerm == LINETERM_UNIX)
				MrBayesPrint ("%s   UNIX line termination\n", spacer);
			else
				{
				MrBayesPrint ("%s   Unknown line termination\n", spacer);
				goto errorExit;
				}
			}
		else
			{
			if (LineTermType (fp) != lineTerm)
				{
				MrBayesPrint ("%s   Nonconforming line termination\n", spacer);
				MrBayesPrint ("%s   Make sure all input files have the same line termination (Mac, DOS, or UNIX)\n", spacer);
				goto errorExit;
				}
			}
		
		if (n == 0)
			{
			/* find length of longest line */
			longestLineLength = LongestLine (fp);
			MrBayesPrint ("%s   Longest line length = %d\n", spacer, longestLineLength);
			longestLineLength += 10;

			/* allocate a string long enough to hold a line */
			if (memAllocs[ALLOC_SUMPSTRING] == YES)
				{
				MrBayesPrint ("%s   Sump string is already allocated\n", spacer);
				goto errorExit;
				}
			DEBUG("allocating s size: %u\n",longestLineLength*sizeof(char));
			s = (char *)SafeMalloc((size_t) (longestLineLength * sizeof(char)));
			if (!s)
				{
				MrBayesPrint ("%s   Problem allocating string for reading sump file\n", spacer);
				goto errorExit;
				}
			DEBUG("allocating headerLine size: %u\n",longestLineLength*sizeof(char));
			headerLine = (char *)SafeMalloc((size_t) ((longestLineLength)* sizeof(char)));
			if (!headerLine)
				{
				MrBayesPrint ("%s   Problem allocating headerLine for reading sump file\n", spacer);
				goto errorExit;
				}
			DEBUG("allocating headerNames size: %u\n",(longestLineLength+40) * sizeof(char));
			headerNames = (char *)SafeMalloc((size_t) ((longestLineLength+40) * sizeof(char)));
			if (!headerNames)
				{
				MrBayesPrint ("%s   Problem allocating headerNames for reading sump file\n", spacer);
				goto errorExit;
				}
			for (i=0; i<longestLineLength+40; i++)
				headerNames[i] = ' ';
			headerNames[longestLineLength+40-1] = '\0';
			memAllocs[ALLOC_SUMPSTRING] = YES;
			}
		else
			{
			/* find length of longest line */
			i = LongestLine (fp) + 10;
			/* reallocate memory if necessary */
			DEBUG("longestLineLength=%d\n",i-10);

			if (i > longestLineLength)
				{
				longestLineLength = i;

				DEBUG("reallocating s new size: %u\n",longestLineLength*sizeof(char));
				s = (char *) realloc((void *) s, (size_t) (longestLineLength * sizeof(char)));
				if (!s)
					{
					MrBayesPrint ("%s   Problem reallocating string for reading sump file\n", spacer);
					goto errorExit;
					}
				DEBUG("reallocating headerLine newsize: %u\n",longestLineLength*sizeof(char));
				headerLine = (char *)realloc((void *) headerLine, (size_t) ((longestLineLength)* sizeof(char)));
				if (!headerLine)
					{
					MrBayesPrint ("%s   Problem reallocating headerLine for reading sump file\n", spacer);
					goto errorExit;
					}

				/* We do not want to reallocate the headerNames character array -- this string holds the header
				 * tokens from the first file and its contents should not change when we read in a new file.
				 * When this code does execute, after reallocation, the contents of the array are cleared,
				 * this causes the header comparison at line 612 to fail, since it compares the headers from
				 * the new file to the now blank headerNames array.
				 * -- catchen@cs.uoregon.edu, may 1, 2007
				 */

				/*
				DEBUG("reallocating headerNames newsize: %u\n",(longestLineLength+40) * sizeof(char));
				headerNames = (char *)realloc((void *)headerNames, (size_t) ((longestLineLength+40) * sizeof(char)));
				if (!headerNames)
				  {
				    MrBayesPrint ("%s   Problem allocating headerNames for reading sump file\n", spacer);
				    goto errorExit;
				  }
				for (i=0; i<longestLineLength+40; i++)
				  headerNames[i] = ' ';
				headerNames[longestLineLength+40-1] = '\0';
				*/
				}

			}

		/* close binary file */
		SafeFclose (&fp);
	
		/* open text file */
		if (sumpParams.numRuns == 1)
			sprintf (temp, "%s.p", sumpParams.sumpFileName);
		else
			sprintf (temp, "%s.run%d.p", sumpParams.sumpFileName, n+1);
		
		if ((fp = OpenTextFileR(temp)) == NULL)
			goto errorExit;
	
		/* Check file for appropriate blocks. We want to find the last block
		   in the file and start from there. */
		inSumpComment = NO;
		lineNum = lastNonDigitLine = numParamLines = 0;
		while (fgets (s, longestLineLength, fp) != NULL)
			{
			sumpTokenP = &s[0];
			allDigitLine = YES;
			lastTokenWasDash = NO;
			nNumbersOnThisLine = 0;
			do
				{
				DEBUG("sumpTokenP offset from s: %d\n",(int) sumpTokenP - (int) s)
				GetSumpToken (&tokenType, &sumpTokenP, sumpToken);
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
				strcpy (headerLine, s);
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
			if (sumpParams.numRuns > 1)
				MrBayesPrint ("%s   Unterminated comment in file \"%s.run%d.p\"\n", spacer, sumpParams.sumpFileName, n+1);
			else
				MrBayesPrint ("%s   Unterminated comment in file \"%s\"\n", spacer, sumpParams.sumpFileName);
			goto errorExit;
			}
		if (numParamLines <= 0)
			{
			if (sumpParams.numRuns > 1)
				MrBayesPrint ("%s   No parameters were found in file \"%s.run%d.p\"\n", spacer, sumpParams.sumpFileName, n+1);
			else
				MrBayesPrint ("%s   No parameters were found in file \"%s\"\n", spacer, sumpParams.sumpFileName);
			goto errorExit;
			}

		/* calculate burnin */
        if (sumpParams.relativeBurnin == NO)
            burnin = sumpParams.sumpBurnIn;
        else
            burnin = (int) (sumpParams.sumpBurnInFraction * numParamLines);

        if (n == 0)
			{
			if (burnin > numParamLines)
				{
				MrBayesPrint ("%s   No parameters are sampled as the burnin exceeds the number of lines in last block\n", spacer);
				MrBayesPrint ("%s   Try setting burnin to a number less than %d\n", spacer, numParamLines);
				goto errorExit;
				}
			
			/* tell the user that everything is fine */
			if (sumpParams.numRuns > 1)
				MrBayesPrint ("%s   Found %d parameter lines in file \"%s.run%d.p\"\n", spacer, numParamLines, sumpParams.sumpFileName, n+1);
			else
				MrBayesPrint ("%s   Found %d parameter lines in file \"%s\"\n", spacer, numParamLines, sumpParams.sumpFileName);
			if (burnin > 0)
				MrBayesPrint ("%s   Of the %d lines, %d of them will be summarized (starting at line %d)\n", spacer, numParamLines, numParamLines - burnin, lastNonDigitLine + burnin + 1);
			else
				MrBayesPrint ("%s   All %d lines will be summarized (starting at line %d)\n", spacer, numParamLines, lastNonDigitLine+1);
			MrBayesPrint ("%s   (Only the last set of lines will be read, in case multiple\n", spacer);
			MrBayesPrint ("%s   parameter blocks are present in the same file.)\n", spacer);

			firstFileNumParamLines = numParamLines;
			}

		else
			{
			/* check that we find the same number of parameter lines in all files */
			if (numParamLines != firstFileNumParamLines)
				{
				MrBayesPrint ("%s   Incompatible files: First file had %d while file %d had %d post-burnin samples\n", spacer, firstFileNumParamLines, n+1, numParamLines); 
				}
			}

		/* Calculate and check the number of columns and rows for the file */
		(void)fseek(fp, 0L, 0);	
		for (lineNum=0; lineNum<lastNonDigitLine+burnin; lineNum++)
		  if(fgets (s, longestLineLength, fp)==NULL) 
#ifdef PAUL
		    printf("%s:%d: fgets failed\n",__FILE__,__LINE__);
#else
		    {} 
#endif
		inSumpComment = NO;
		nLines = 0;
		numRows = numColumns = 0;
		while (fgets (s, longestLineLength, fp) != NULL)
			{
			sumpTokenP = &s[0];
			allDigitLine = YES;
			lastTokenWasDash = NO;
			nNumbersOnThisLine = 0;
			do
				{
				GetSumpToken (&tokenType, &sumpTokenP, sumpToken);
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
				if (sumpParams.numRuns > 1)
					MrBayesPrint ("%s   Found a line with non-digit characters (line %d) in file %d\n", spacer, lineNum, n+1);
				else
					MrBayesPrint ("%s   Found a line with non-digit characters (line %d)\n", spacer, lineNum);
				goto errorExit;
				}
			else
				{
				if (nNumbersOnThisLine > 0)
					{
					nLines++;
					if (nLines == 1 && n == 0)
						firstNumCols = nNumbersOnThisLine;
					else
						{
						if (nNumbersOnThisLine != firstNumCols)
							{
							if (sumpParams.numRuns == 1)
								MrBayesPrint ("%s   Number of columns is not even (%d in first line and %d in %d line)\n", spacer, firstNumCols, nNumbersOnThisLine, lineNum);
							else
								MrBayesPrint ("%s   Number of columns is not even (%d in first line of first file and %d in %d line of file %d)\n", spacer, firstNumCols, nNumbersOnThisLine, lineNum, i+1);
							goto errorExit;
							}
						}
					}
				}
			}
		numRows = nLines;
		numColumns = firstNumCols;
		if (n==0)
			{
			MrBayesPrint ("%s   %d rows and %d columns in each row\n", spacer, numRows, numColumns);
			if (sumpParams.numRuns > 1)
				MrBayesPrint ("%s   Expecting the same layout in all subsequent files\n", spacer);

			/* allocate space to hold parameter information */
			if (numRows == 0 || numColumns == 0)
				{
				MrBayesPrint ("%s   The number of rows or columns is equal to zero\n", spacer);
				goto errorExit;
				}
			if (memAllocs[ALLOC_SUMPINFO] == YES)
				{
				MrBayesPrint ("%s   Sump string is already allocated\n", spacer);
				goto errorExit;
				}
			parameterValues = (MrBFlt *)SafeMalloc((size_t) (sumpParams.numRuns * numRows * numColumns * sizeof(MrBFlt)));
			if (!parameterValues)
				{
				MrBayesPrint ("%s   Problem allocating parameterValues\n", spacer);
				goto errorExit;
				}
			memAllocs[ALLOC_SUMPINFO] = YES;
			for (i=0; i<sumpParams.numRuns*numRows*numColumns; i++)
				parameterValues[i] = 0.0;
			firstFileNumRows = numRows;
			firstFileNumColumns = numColumns;
			}
		else
			{
			if (numRows != firstFileNumRows || numColumns != firstFileNumColumns)
				{
				MrBayesPrint ("%s   Incompatible files: Different number of columns or rows in file %d\n", spacer, n+1);
				}
			}

		/* Now we read the file for real. First, rewind file pointer to beginning of file... */
		(void)fseek(fp, 0L, 0);	
	
		/* ...and fast forward to beginning of last unburned parameter line. */
		for (lineNum=0; lineNum<lastNonDigitLine+burnin; lineNum++)
		  if(fgets (s, longestLineLength, fp)==0) 
#ifdef PAUL
		    printf("%s:%d error in fgets\n",__FILE__,__LINE__);
#else
		    {} 
#endif
		/* ...and parse file, line-by-line. We are only parsing lines that have digits that should be read. */
		parameterPointer = parameterValues + n*numRows*numColumns;
		inSumpComment = NO;
		numLinesToRead = numParamLines - burnin;
		numLinesRead = j = 0;
		while (fgets (s, longestLineLength, fp) != NULL)
			{
			sumpTokenP = &s[0];
			allDigitLine = YES;
			lastTokenWasDash = NO;
			nNumbersOnThisLine = 0;
			do
				{
				GetSumpToken (&tokenType, &sumpTokenP, sumpToken);
				if (IsSame("[", sumpToken) == SAME)
					inSumpComment = YES;
				if (IsSame("]", sumpToken) == SAME)
					inSumpComment = NO;
				if (inSumpComment == NO)
					{
					if (tokenType == NUMBER)
						{
						/* read the information from this line */
						if (j >= numRows * numColumns)
							{
							if (sumpParams.numRuns > 1)
								MrBayesPrint ("%s   Too many parameter values read in (%d) in file %d\n", spacer, j, n+1);
							else
								MrBayesPrint ("%s   Too many parameter values read in (%d)\n", spacer, j);
							goto errorExit;
							}
						sscanf (sumpToken, "%lf", &tempD);
						if (lastTokenWasDash == YES)
							tempD *= -1.0;
						parameterPointer[numLinesRead * numColumns + nNumbersOnThisLine] = tempD;
						j++;
						nNumbersOnThisLine++;
						lastTokenWasDash = NO;
						}
					else if (tokenType == DASH)
						{
						lastTokenWasDash = YES;
						}
					else if (tokenType != UNKNOWN_TOKEN_TYPE)
						{
						/* we have a problem */
						if (sumpParams.numRuns > 1)
							MrBayesPrint ("%s   Found a line with non-digit characters (line %d) in file %d\n", spacer, lineNum, n+1);
						else
							MrBayesPrint ("%s   Found a line with non-digit characters (line %d)\n", spacer, lineNum);
						goto errorExit;
						}
					}
				} while (*sumpToken);
			lineNum++;
			if (nNumbersOnThisLine > 0)
				numLinesRead++;
			}
	
		/* tell user how many lines were successfully read */
		if (sumpParams.numRuns == 1)
			MrBayesPrint ("%s   Successfully read %d lines from last parameter block\n", spacer, numLinesRead);
		else
			MrBayesPrint ("%s   Successfully read %d lines from last parameter block of file %d\n", spacer, numLinesRead, n+1);
	
		/* Check how many parameter line was read in. */
		if (numLinesRead != numParamLines - burnin)
			{
			MrBayesPrint ("%s   Unable to read all lines that should contain parameter samples\n", spacer);
			goto errorExit;
			}
				
		if (n == 0)
			{
		    /* separate header line into titles for each column. GetHeaders function will
			   place the parsed headers into the 'headerNames' string, separated by the '|' character. */ 
			if (GetHeaders (headerLine, &nHeaders) == ERROR)
			    goto errorExit;

			/* get length of longest header */
			longestHeader = 9; /* 9 is the length of the word "parameter" (for printing table) */
			for (i=0; i<nHeaders; i++)
				{
				if (GetNameFromString (headerNames, temp, i+1) == ERROR)
					{
					MrBayesPrint ("%s   Error getting header names \n", spacer);
					goto errorExit;
					}
				len = (int) strlen(temp);
				if (len > longestHeader)
					longestHeader = len;
				}
			}
		else
			{
			/* Check the the headers in file N match the headers in the first file */
			for (i=0; i<nHeaders; i++)
				{
				GetNameFromString (headerNames, temp, i+1);

				if (i==0)
					s1 = strtok (headerLine, " \t");
				else
					s1 = strtok (NULL, " \t\n");

				if (s1 == NULL || strcmp (temp, s1) != 0)
					{
					MrBayesPrint ("%s   Non-matching headers (%s / %s) in file %d\n", spacer, temp, s1, n+1);
					goto errorExit;
					}
				}
			}
		SafeFclose (&fp);
		}	/* next file */
			
		
	if (sumpParams.plot == YES)
		{
		/* Print header */
		PrintPlotHeader (stdout);
		if (sumpParams.printToFile == YES)
			PrintPlotHeader (fpOut);

		if (sumpParams.numRuns > 1)
			{
			if (sumpParams.allRuns == YES)
				{
				for (i=0; i<sumpParams.numRuns; i++)
					{
					MrBayesPrint ("\n%s   Samples from run %d:\n", spacer, i+1);
					if (PrintPlot (stdout, parameterValues + i*numRows*numColumns) == ERROR)
						goto errorExit;
					if (sumpParams.printToFile == YES)
						{
						MrBayesPrintf (fpOut, "\n%s   Samples from run %d:\n", spacer, i+1);
						if (PrintPlot (fpOut, parameterValues + i*numRows*numColumns) == ERROR)
							goto errorExit;
						}
					}
				}
			
			if (PrintOverlayPlot (stdout) == ERROR)
				goto errorExit;
			if (sumpParams.printToFile == YES)
				{
				if (PrintOverlayPlot (fpOut) == ERROR)
					goto errorExit;
				}
			}
		else
			{
			if (PrintPlot (stdout, parameterValues) == ERROR)
				goto errorExit;
			if (sumpParams.printToFile == YES)
				{
				if (PrintPlot (fpOut, parameterValues) == ERROR)
					goto errorExit;
				}
			}
		}
			
	if (sumpParams.margLike == YES)
		{
		/* deal with lnL */
		for (i=0; i<nHeaders; i++)
			{
			if (GetNameFromString (headerNames, temp, i+1) == ERROR)
				{
				MrBayesPrint ("%s   Error getting header names \n", spacer);
				goto errorExit;
				}
			len = (int) strlen(temp);
			if (IsSame (temp, "lnL") != SAME)
				continue;

			oneUnreliable = NO;
			for (n=0; n<sumpParams.numRuns; n++)
				{
				if (SortColumn (parameterValues + n*numRows*numColumns, numRows, numColumns, i) == ERROR)
					goto errorExit;
				unreliable = NO;
				if (HarmonicArithmeticMean(parameterValues+n*numRows*numColumns, numRows, numColumns, i, &mean, &harm_mean) == ERROR)
					{
					unreliable = YES;
					oneUnreliable = YES;
					}
				if (sumpParams.numRuns == 1)
					{
					MrBayesPrint ("\n");
					MrBayesPrint ("%s   Estimated marginal likelihoods for run sampled in file \"%s.p\":\n", spacer, sumpParams.sumpFileName);
					MrBayesPrint ("%s      (Use the harmonic mean for Bayes factor comparisons of models)\n\n", spacer, sumpParams.sumpFileName);
					MrBayesPrint ("%s   Arithmetic mean   Harmonic mean\n", spacer);
					MrBayesPrint ("%s   --------------------------------\n", spacer);
					if (unreliable == NO)
						MrBayesPrint ("%s     %10.2lf        %10.2lf\n", spacer, mean, harm_mean);
					else
						MrBayesPrint ("%s     %10.2lf *      %10.2lf *\n", spacer, mean, harm_mean);
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
						MrBayesPrint ("%s   Run   Arithmetic mean   Harmonic mean\n", spacer);
						MrBayesPrint ("%s   --------------------------------------\n", spacer);
						}
					if (unreliable == NO)
						MrBayesPrint ("%s   %3d     %10.2lf        %10.2lf\n", spacer, n+1, mean, harm_mean);
					else
						MrBayesPrint ("%s   %3d     %10.2lf *      %10.2lf *\n", spacer, n+1, mean, harm_mean);
					}					

				if (sumpParams.printToFile == YES)
					{
					if (sumpParams.numRuns == 1)
						{
						MrBayesPrintf (fpOut, "\n%s   Estimated marginal likelihoods for run sampled in file \"%s\"\n\n", spacer, sumpParams.sumpFileName);
						MrBayesPrintf (fpOut, "%s   Arithmetic mean   Harmonic mean\n", spacer);
						MrBayesPrintf (fpOut, "%s   --------------------------------\n", spacer);
						if (unreliable == NO)
							MrBayesPrintf (fpOut, "%s     %10.2lf        %10.2lf\n", spacer, mean, harm_mean);
						else
							MrBayesPrintf (fpOut, "%s     %10.2lf *      %10.2lf *\n", spacer, mean, harm_mean);
						}
					else
						{
						if (n == 0)
							{
							MrBayesPrintf (fpOut, "%s   Estimated marginal likelihoods for runs sampled in files\n", spacer);
							if (sumpParams.numRuns > 2)
								MrBayesPrintf (fpOut, "%s      \"%s.run1.p\", \"%s.run2.p\", etc:\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
							else /* if (sumpParams.numRuns == 2) */
								MrBayesPrintf (fpOut, "%s      \"%s.run1.p\" and \"%s.run2.p\":\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
							MrBayesPrintf (fpOut, "%s   Run   Arithmetic mean   Harmonic mean\n", spacer);
							MrBayesPrintf (fpOut, "%s   --------------------------------------\n", spacer);
							}
						if (unreliable == NO)
							MrBayesPrintf (fpOut, "%s   %3d     %10.2lf        %10.2lf\n", spacer, n+1, mean, harm_mean);
						else
							MrBayesPrintf (fpOut, "%s   %3d     %10.2lf *      %10.2lf *\n", spacer, n+1, mean, harm_mean);
						}					
					}
				}	/* next run */
			if (sumpParams.numRuns == 1)
				{
				MrBayesPrint ("%s   --------------------------------\n", spacer);
				if (sumpParams.printToFile == YES)
					MrBayesPrintf (fpOut, "      --------------------------------\n");
				}
			else
				{
				if (SortColumn (parameterValues, sumpParams.numRuns*numRows, numColumns, i) == ERROR)
					goto errorExit;
				if (HarmonicArithmeticMean (parameterValues, sumpParams.numRuns*numRows, numColumns, i, &mean, &harm_mean) == ERROR)
					{
					unreliable = YES;
					oneUnreliable = YES;
					}
				else
					unreliable = NO;
				MrBayesPrint ("%s   --------------------------------------\n", spacer);
				if (unreliable == YES)
					MrBayesPrint ("%s   TOTAL   %10.2lf *      %10.2lf *\n", spacer, mean, harm_mean);
				else
					MrBayesPrint ("%s   TOTAL   %10.2lf        %10.2lf\n", spacer, mean, harm_mean);
				MrBayesPrint ("%s   --------------------------------------\n", spacer);
				if (sumpParams.printToFile == YES)
					{
					MrBayesPrintf (fpOut, "%s   --------------------------------------\n", spacer);
					MrBayesPrintf (fpOut, "%s   --------------------------------------\n", spacer);
					if (unreliable == YES)
						MrBayesPrintf (fpOut, "%s   TOTAL   %10.2lf *      %10.2lf *\n", spacer, mean, harm_mean);
					else
						MrBayesPrintf (fpOut, "%s   TOTAL   %10.2lf        %10.2lf\n", spacer, mean, harm_mean);
					MrBayesPrintf (fpOut, "%s   --------------------------------------\n", spacer);
					}
				}
			if (oneUnreliable == YES)
				{
				MrBayesPrint ("%s   * These estimates may be unreliable because \n", spacer);
				MrBayesPrint ("%s     some extreme values were excluded\n\n", spacer);
				if (sumpParams.printToFile == YES)
					{
					MrBayesPrintf (fpOut, "%s   * These estimates may be unreliable because \n", spacer);
					MrBayesPrintf (fpOut, "%s     some extreme values were excluded\n\n", spacer);
					}
				}
			else
				{
				MrBayesPrint ("\n");
				if (sumpParams.printToFile == YES)
					{
					MrBayesPrintf (fpOut, "\n");
					}
				}

			}	/* check next header */
		}

	if (sumpParams.table == YES)
		{
		/* Print parameter information to screen and possibly to file. */
		if (sumpParams.numRuns > 1 && sumpParams.allRuns == YES)
			{
			for (i=0; i<sumpParams.numRuns; i++)
				{
				/* print table header */
				MrBayesPrint ("\n");
				MrBayesPrint ("%s   Model parameter summaries for run sampled in file \"%s.run%d.p\":\n", spacer, sumpParams.sumpFileName, i+1);
				MrBayesPrint ("%s      (Based on %d samples out of a total of %d samples from this analysis)\n\n", spacer, numRows, numRows + burnin);
				if (PrintParamTable (stdout, parameterValues + i*numRows*numColumns, numRows, numColumns, 1) == ERROR)
					goto errorExit;
				if (PrintAATable (stdout, parameterValues + i*numRows*numColumns, numRows, numColumns, 1) == ERROR)
					goto errorExit;
				if (sumpParams.printToFile == YES)
					{
					MrBayesPrintf (fpOut, "\n");
					MrBayesPrintf (fpOut, "%s   Model parameter summaries for run sampled in file \"%s.run%d.p\":\n", spacer, sumpParams.sumpFileName, i+1);
					MrBayesPrintf (fpOut, "%s   (Based on %d samples out of a total of %d samples from this analysis)\n\n", spacer, numRows, numRows + burnin);
					if (PrintParamTable (fpOut, parameterValues + i*numRows*numColumns, numRows, numColumns, 1) == ERROR)
						goto errorExit;
					if (PrintAATable (fpOut, parameterValues + i*numRows*numColumns, numRows, numColumns, 1) == ERROR)
						goto errorExit;
					}
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
			MrBayesPrint ("%s      (Based on a total of %d samples out of a total of %d samples from this analysis)\n\n", spacer, numRows, numRows + burnin);
		else
			{
			MrBayesPrint ("%s      (Summaries are based on a total of %d samples from %d runs)\n", spacer, sumpParams.numRuns*numRows, sumpParams.numRuns);
			MrBayesPrint ("%s      (Each run produced %d samples of which %d samples were included)\n\n", spacer, numRows + burnin, numRows);
			}

		if (PrintParamTable (stdout, parameterValues, numRows, numColumns, sumpParams.numRuns) == ERROR)
			goto errorExit;
		if (PrintAATable (stdout, parameterValues, numRows, numColumns, sumpParams.numRuns) == ERROR)
			goto errorExit;

		if (sumpParams.printToFile == YES)
			{
			MrBayesPrintf (fpOut, "\n");
			if (sumpParams.numRuns == 1)
				MrBayesPrintf (fpOut, "%s   Model parameter summaries for run sampled in file \"%s\":\n\n", spacer, sumpParams.sumpFileName);
			else if (sumpParams.numRuns == 2)
				{
				MrBayesPrint ("%s   Model parameter summaries over the runs sampled in files\n", spacer);
				MrBayesPrint ("%s      \"%s.run1.p\" and \"%s.run2.p\":\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
				}
			else
				{
				MrBayesPrintf (fpOut, "%s   Model parameter summaries over all %d runs sampled in files\n", spacer, sumpParams.numRuns);
				MrBayesPrintf (fpOut, "%s      \"%s.run1.p\", \"%s.run2.p\" etc:\n", spacer, sumpParams.sumpFileName, sumpParams.sumpFileName);
				}
			if (sumpParams.numRuns == 1)
				MrBayesPrintf (fpOut, "%s   (Based on a total of %d samples out of a total of %d samples from this analysis)\n\n", spacer, numRows, numRows + burnin);
			else
				{
				MrBayesPrintf (fpOut, "%s      (Summaries are based on a total of %d samples from %d runs)\n", spacer, sumpParams.numRuns*numRows, sumpParams.numRuns);
				MrBayesPrintf (fpOut, "%s      (Each run produced %d samples of which %d samples were included)\n\n", spacer, numRows + burnin, numRows);
				}
		
			if (PrintParamTable (fpOut, parameterValues, numRows, numColumns, sumpParams.numRuns) == ERROR)
				goto errorExit;
			if (PrintAATable (fpOut, parameterValues, numRows, numColumns, sumpParams.numRuns) == ERROR)
				goto errorExit;
			}
		}

	/* free memory and file pointers */
	if (memAllocs[ALLOC_SUMPSTRING] == YES)
		{
		free (s);
		free (headerLine);
		free (headerNames);
		memAllocs[ALLOC_SUMPSTRING] = NO;
		}
	if (memAllocs[ALLOC_SUMPINFO] == YES)
		{
		free (parameterValues);
		memAllocs[ALLOC_SUMPINFO] = NO;
		}
	SafeFclose (&fp);
	SafeFclose (&fpOut);
	expecting = Expecting(COMMAND);
	
#	if defined (MPI_ENABLED)
		}
	ierror = MPI_Barrier (MPI_COMM_WORLD);
	if (ierror != MPI_SUCCESS)
		{
		MrBayesPrint ("%s   Problem at chain barrier\n", spacer);
		goto errorExit;
		}
#	endif

	return (NO_ERROR);
	
	errorExit:
		expecting = Expecting(COMMAND);
		if (memAllocs[ALLOC_SUMPSTRING] == YES)
			{
			if (s)
				free (s);
			if (headerLine)
				free (headerLine);
			if (headerNames)
				free (headerNames);
			memAllocs[ALLOC_SUMPSTRING] = NO;
			}
		if (memAllocs[ALLOC_SUMPINFO] == YES)
			{
			free (parameterValues);
			memAllocs[ALLOC_SUMPINFO] = NO;
			}
		SafeFclose (&fp);
		SafeFclose (&fpOut);
		strcpy (spacer, "");
		strcpy (sumpToken, "Sump");
		i = 0;
		if (FindValidCommand (sumpToken, &i) == ERROR)
			MrBayesPrint ("%s   Could not find sump\n", spacer);
	
#	if defined (MPI_ENABLED)
		ierror = MPI_Barrier (MPI_COMM_WORLD);
		if (ierror != MPI_SUCCESS)
			{
			MrBayesPrint ("%s   Problem at chain barrier\n", spacer);
			/* We are already returning an error anyway */
			}
#	endif

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
				strcpy (sumpParams.sumpFileName, tkn);
				MrBayesPrint ("%s   Setting sump filename to %s\n", spacer, sumpParams.sumpFileName);
				strcpy (sumpParams.sumpOutfile, tkn);
				strcat (sumpParams.sumpOutfile, ".stat");
				MrBayesPrint ("%s   Setting sump output file to %s\n", spacer, sumpParams.sumpOutfile);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Relburnin (sumpParams.relativeBurnin) ********************************************************/
		else if (!strcmp(parmName, "Relburnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumpParams.relativeBurnin = YES;
					else
						sumpParams.relativeBurnin = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for Relburnin\n", spacer);
					free(tempStr);
					return (ERROR);
					}
				if (sumpParams.relativeBurnin == YES)
					MrBayesPrint ("%s   Using relative burnin (a fraction of samples discarded).\n", spacer);
				else
					MrBayesPrint ("%s   Using absolute burnin (a fixed number of samples discarded).\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				{
				free (tempStr);
				return (ERROR);
				}
			}
		/* set Burnin (sumpParams.sumpBurnIn) ***********************************************************/
		else if (!strcmp(parmName, "Burnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
                sumpParams.sumpBurnIn = tempI;
				MrBayesPrint ("%s   Setting sump burn-in to %d\n", spacer, sumpParams.sumpBurnIn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				{
				free(tempStr);
				return (ERROR);
				}
			}
		/* set Burninfrac (sumpParams.sumpBurnInFraction) ************************************************************/
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
					free(tempStr);
					return (ERROR);
					}
				if (tempD > 0.50)
					{
					MrBayesPrint ("%s   Burnin fraction too high (> 0.50)\n", spacer);
					free(tempStr);
					return (ERROR);
					}
                sumpParams.sumpBurnInFraction = tempD;
				MrBayesPrint ("%s   Setting burnin fraction to %.2f\n", spacer, sumpParams.sumpBurnInFraction);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else 
				{
				free(tempStr);
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
				sscanf (tkn, "%s", tempStr);
				strcpy (sumpParams.sumpOutfile, tempStr);
				MrBayesPrint ("%s   Setting sump output file name to \"%s\"\n", spacer, sumpParams.sumpOutfile);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Plot (sumpParams.plot) ********************************************************/
		else if (!strcmp(parmName, "Plot"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumpParams.plot = YES;
					else
						sumpParams.plot = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for plot (valid arguments are 'yes' and 'no')\n", spacer);
					return (ERROR);
					}
				if (sumpParams.plot == YES)
					MrBayesPrint ("%s   Setting sump to print likelihood generation plot\n", spacer);
				else
					MrBayesPrint ("%s   Setting sump not to print likelihood generation plot\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Table (sumpParams.table) ********************************************************/
		else if (!strcmp(parmName, "Table"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumpParams.table = YES;
					else
						sumpParams.table = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for 'table' (valid arguments are 'yes' and 'no')\n", spacer);
					return (ERROR);
					}
				if (sumpParams.plot == YES)
					MrBayesPrint ("%s   Setting sump to print parameter table\n", spacer);
				else
					MrBayesPrint ("%s   Setting sump not to print parameter table\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Marglike (sumpParams.marglike) ********************************************************/
		else if (!strcmp(parmName, "Marglike"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumpParams.margLike = YES;
					else
						sumpParams.margLike = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for plot (valid arguments are 'yes' and 'no')\n", spacer);
					return (ERROR);
					}
				if (sumpParams.plot == YES)
					MrBayesPrint ("%s   Setting sump to print marginal likelihood estimates\n", spacer);
				else
					MrBayesPrint ("%s   Setting sump not to print marginal likelihood estimates\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Printtofile (sumpParams.printToFile) ********************************************************/
		else if (!strcmp(parmName, "Printtofile"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumpParams.printToFile = YES;
					else
						sumpParams.printToFile = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for printtofile (valid arguments are 'yes' and 'no')\n", spacer);
					return (ERROR);
					}
				if (sumpParams.printToFile == YES)
					MrBayesPrint ("%s   Setting sump to print output to file\n", spacer);
				else
					MrBayesPrint ("%s   Setting sump not to print output to file\n", spacer);
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





int GetHeaders (char *h, int *nHeaders)

{

	int			i, j, x, wasLastCharSpace;
	char		temp[100];

	wasLastCharSpace = YES;
	i = j = (*nHeaders) = 0;
	while (h[i] != '\0')
		{
		if (h[i] == ' ' || h[i] == '\t' || h[i] == '\n' || h[i] == '\r')
			{
			if (wasLastCharSpace == NO)
				{
				temp[j++] = '\0';
				if (AddToString (temp, headerNames, &x) == ERROR)
					{
					MrBayesPrint ("%s   Error adding header to list of headers \n", spacer);
					return ERROR;
					}
				(*nHeaders)++;
				}
			wasLastCharSpace = YES;
			j = 0;
			}
		else
			{
			temp[j++] = h[i];
			wasLastCharSpace = NO;
			}
		i++;
		}

	if (wasLastCharSpace == NO)
		{
		temp[j++] = '\0';
		if (AddToString (temp, headerNames, &x) == ERROR)
			{
			MrBayesPrint ("%s   Error adding header to list of headers \n", spacer);
			return ERROR;
			}
		(*nHeaders)++;
		}
		
#	if 0
	for (i=0; i<(*nHeaders); i++)
		{
		if (GetNameFromString (headerNames, temp, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Error getting header names \n", spacer);
			return (ERROR);
			}
		printf ("%4d -> '%s'\n", i, temp);
		}
#	endif
		
	return (NO_ERROR);
	
}





void GetSumpToken (int *tokenType, char **tknSource, char *token)

{
		
	int		allNumbers, foundExp, foundExpSign;
	
	(*tokenType) = 0;
	
	while (IsWhite(**tknSource) == 1 || IsWhite(**tknSource) == 2)
		{
		if (IsWhite(**tknSource) == 2)
			{
			*tokenType = RETURNSYMBOL;
			foundNewLine = YES;
			/* MrBayesPrint ("RETURN\n"); */
			}
		(*tknSource)++;
		}
	
	*tokenType = UNKNOWN_TOKEN_TYPE;
	if (IsIn(**tknSource,"="))
		{
		*token++ = *(*tknSource)++;
		*tokenType = EQUALSIGN;
		}
	else if (IsIn(**tknSource,";"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = SEMICOLON;
		}
	else if (IsIn(**tknSource,":"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = COLON;
		}
	else if (IsIn(**tknSource,","))
		{
		*token++ = *(*tknSource)++;
		*tokenType = COMMA;
		}
	else if (IsIn(**tknSource,"#"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = POUNDSIGN;
		}
	else if (IsIn(**tknSource,"("))
		{
		*token++ = *(*tknSource)++;
		*tokenType = LEFTPAR;
		}
	else if (IsIn(**tknSource,")"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = RIGHTPAR;
		}
	else if (IsIn(**tknSource,"{"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = LEFTCURL;
		}
	else if (IsIn(**tknSource,"}"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = RIGHTCURL;
		}
	else if (IsIn(**tknSource,"["))
		{
		*token++ = *(*tknSource)++;
		*tokenType = LEFTCOMMENT;
		}
	else if (IsIn(**tknSource,"]"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = RIGHTCOMMENT;
		}
	else if (IsIn(**tknSource,"?"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = QUESTIONMARK;
		}
	else if (IsIn(**tknSource,"-"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = DASH;
		}
	else if (IsIn(**tknSource,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789."))
		{
		if (!IsIn(**tknSource,"0123456789."))
			allNumbers = FALSE;
        else
            allNumbers = TRUE;
        foundExp = foundExpSign = FALSE;
		*token++ = *(*tknSource)++;
		while(IsIn(**tknSource,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.-+"))
			{
			if(allNumbers == TRUE && !IsIn((*tknSource)[-1],"Ee") && **tknSource=='-')
                break;
            else if (allNumbers == TRUE && IsIn(**tknSource,"Ee") && foundExp == NO)
                foundExp = TRUE;
            else if (allNumbers == TRUE && IsIn(**tknSource,"+-") && IsIn((*tknSource)[-1],"Ee"))
                foundExpSign = TRUE;
            else if (!IsIn(**tknSource,"0123456789."))
				allNumbers = FALSE;
		    *token++ = *(*tknSource)++;
			}
		if (allNumbers == TRUE)
			*tokenType = NUMBER;
		else
			*tokenType = ALPHA;
		}
	else if (IsIn(**tknSource,"*"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = ASTERISK;
		}
	else if (IsIn(**tknSource,"/"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = FORWARDSLASH;
		}
	else if (IsIn(**tknSource,"'\\'"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = BACKSLASH;
		}
	else if (IsIn(**tknSource,"!"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = EXCLAMATIONMARK;
		}
	else if (IsIn(**tknSource,"%"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = PERCENT;
		}
	else if (IsIn(**tknSource,"\""))
		{
		*token++ = *(*tknSource)++;
		*tokenType = QUOTATIONMARK;
		}
	else if (IsIn(**tknSource,"&~+^$@|{}`><"))
		{
		*token++ = *(*tknSource)++;
		*tokenType = WEIRD;
		}

	*token = '\0';
	
}





int HarmonicArithmeticMean (MrBFlt *paramVals, int nRows, int nCols, int col, MrBFlt *mean, MrBFlt *harm_mean)

{

	int				i, reliable;
	MrBFlt			a, aOld, x, y, scaler, n;

	reliable = YES;
	
	scaler = paramVals[(nRows-1) * nCols + col] - 100.0;
	a = aOld = n = 0.0;
	for (i=0; i<nRows; i++)
		{
		x = paramVals[i * nCols + col];
		y = x;
		y -= scaler;
		if (y < -100.0)
			{
			reliable = NO;
			continue;
			}
		else
			x = (MrBFlt) exp(y);
			
		if (n < 0.5)
			a = x;
		else
			{
			aOld = a;
			a = aOld + (x - aOld) / (n+(MrBFlt)1.0);
			}
		n += 1.0;
		}

	/* arithmetic mean */
	(*mean) = (MrBFlt) log(a) + scaler;
	
	scaler = (MrBFlt) (0.0 - paramVals[col]) - 100.0;
	a = aOld = n = 0.0;
	for (i=0; i<nRows; i++)
		{
		x = (MrBFlt) (0.0 - paramVals[i * nCols + col]);
		y = x;
		y -= scaler;
		if (y < -100.0)
			{
			reliable = NO;
			continue;
			}
		else
			x = (MrBFlt) exp(y);
			
		if (n < 0.5)
			a = x;
		else
			{
			aOld = a;
			a = aOld + (x - aOld) / (n+(MrBFlt)1.0);
			}
		n += (MrBFlt) 1.0;
		}

	/* harmonic mean */
	(*harm_mean) = - (MrBFlt) log(a) - scaler;

	if (reliable == YES)
		return (NO_ERROR);
	else
		return (ERROR);
	
}





void MeanVariance (MrBFlt *paramVals, int nRows, int nCols, int col, MrBFlt *m, MrBFlt *v)

{

	int				i;
	MrBFlt			a, aOld, s, x;

	a = s = 0.0;
	for (i=0; i<nRows; i++)
		{
		x = paramVals[i * nCols + col];
		aOld = a;
		a += (x - a) / (MrBFlt) (i + 1);
		s += (x - a) * (x - aOld);
		}

	/* mean */
	(*m) = a;
	
	/* variance */
	if (nRows <= 1)
		(*v) = 0.0;
	else
		(*v) = s / (nRows - 1);
			
}





MrBFlt PotentialScaleReduction (MrBFlt *paramVals, int nRows, int nCols, int nRuns, int col)

{

	int				i, k;
	MrBFlt			aW, aOldW, sW, aB, aOldB, sB, x, R;

	aB = sB = sW = 0.0;
	for (k=0; k<nRuns; k++)
		{
		aW = 0.0;
		for (i=0; i<nRows; i++)
			{
			x = paramVals[(k*nRows + i) * nCols + col];
			aOldW = aW;
			aW += (x - aW) / (MrBFlt) (i + 1);
			if (i != 0)
				sW += (x - aW) * (x - aOldW);
			}
		x = aW;
		aOldB = aB;
		aB += (x - aB) / (MrBFlt) (k + 1);
		if (k!=0)
			sB += (x - aB) * (x - aOldB);
		}

	sB = sB / (MrBFlt) (nRuns - 1);
	sW = sW / (MrBFlt) (nRuns * nRows - nRuns);

	if (sW > 0.0)
		{
		R = ((MrBFlt)(nRows - 1) / (MrBFlt) (nRows)) + ((MrBFlt)(nRuns + 1) / (MrBFlt) (nRuns)) * (sB / sW);
		return sqrt(R);
		}
	else
		return -1.0;
}




/* PrintAATable: Print amino acid table to file fp */
int PrintAATable (FILE *fp, MrBFlt *parmVals, int nRows, int nCols, int nRuns)
{
	int		i, j, k, len, startReading;
	int		aaModelCounts[20];
	MrBFlt	f, sum[20], ssq[20], stddev[20], prob[20];
	char	parts[100], temp[100];
	
	for (i=0; GetNameFromString (headerNames, temp, i+1) != ERROR; i++)
		{
		len = (int) strlen(temp);
		if (IsSame (temp, "Aamodel") == SAME || IsSame (temp, "Aamodel") == CONSISTENT_WITH)
			{
			MrBayesPrintf (fp, "\n\n");
			if (IsSame (temp, "Aamodel") == SAME)
				MrBayesPrintf (fp, "%s   Amino acid model probabilities:\n\n", spacer);
			else
				{
				startReading = NO;
				j = k = 0;
				while (temp[j] != '\0')
					{
					if (startReading == YES && temp[j] != '}')
						parts[k++] = temp[j];
					if (temp[j] == '{')
						startReading = YES;
					else if (temp[j] == '}')
						startReading = NO;
					j++;
					}
				parts[k] = '\0';
				MrBayesPrintf (fp, "%s   Amino acid model probabilities for partition(s) %s:\n\n", spacer, parts);
				}
			if (nRuns == 1)
				{
				for (j=0; j<20; j++)
					aaModelCounts[j] = 0;
				for (j=0; j<numRows; j++)
					aaModelCounts[(int)(parmVals[j * nCols + i])]++;
				MrBayesPrintf (fp, "%s              Posterior\n", spacer);
				MrBayesPrintf (fp, "%s   Model     Probability\n", spacer);
				MrBayesPrintf (fp, "%s   -------------------\n", spacer);
				MrBayesPrintf (fp, "%s   Poisson      %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_POISSON] / numRows);
				MrBayesPrintf (fp, "%s   Jones        %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_JONES]   / numRows);
				MrBayesPrintf (fp, "%s   Dayhoff      %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_DAY]     / numRows);
				MrBayesPrintf (fp, "%s   Mtrev        %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_MTREV]   / numRows);
				MrBayesPrintf (fp, "%s   Mtmam        %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_MTMAM]   / numRows);
				MrBayesPrintf (fp, "%s   Wag          %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_WAG]     / numRows);
				MrBayesPrintf (fp, "%s   Rtrev        %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_RTREV]   / numRows);
				MrBayesPrintf (fp, "%s   Cprev        %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_CPREV]   / numRows);
				MrBayesPrintf (fp, "%s   Vt           %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_VT]      / numRows);
				MrBayesPrintf (fp, "%s   Blosum       %1.3lf\n", spacer, (MrBFlt)aaModelCounts[AAMODEL_BLOSUM]  / numRows);
				MrBayesPrintf (fp, "%s   -------------------\n", spacer);
				}
			else /* if (nRuns > 1) */
				{
				for (j=0; j<20; j++)
					prob[j] = sum[j] = ssq[j] = 0.0;
				for (j=0; j<nRuns; j++)
					{
					for (k=0; k<20; k++)
						aaModelCounts[k] = 0;
					for (k=0; k<nRows; k++)
						{
						aaModelCounts[(int)(parmVals[j*nRows*nCols + k*nCols + i])]++;
						}
					for (k=0; k<20; k++)
						{
						f = (aaModelCounts[k] / (MrBFlt) nRows);
						prob[k] += (f - prob[k])/ (MrBFlt) (j + 1);
						sum[k] += f;
						ssq[k] += f * f;
						}
					}
				for (j=0; j<20; j++)
					{
					f = ssq[j] - sum[j] * sum[j] / (MrBFlt) nRuns;
					f /= (nRuns - 1);
					if (f <= 0.0)
						stddev[j] = 0.0;
					else
						stddev[j] = sqrt (f);
					}
				MrBayesPrintf (fp, "%s              Posterior     Standard\n", spacer);
				MrBayesPrintf (fp, "%s   Model     Probability    Deviation\n", spacer);
				MrBayesPrintf (fp, "%s   ----------------------------------\n", spacer);
				MrBayesPrintf (fp, "%s   Poisson      %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_POISSON], stddev[AAMODEL_POISSON]);
				MrBayesPrintf (fp, "%s   Jones        %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_JONES], stddev[AAMODEL_JONES]);
				MrBayesPrintf (fp, "%s   Dayhoff      %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_DAY], stddev[AAMODEL_DAY]);
				MrBayesPrintf (fp, "%s   Mtrev        %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_MTREV], stddev[AAMODEL_MTREV]);
				MrBayesPrintf (fp, "%s   Mtmam        %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_MTMAM], stddev[AAMODEL_MTMAM]);
				MrBayesPrintf (fp, "%s   Wag          %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_WAG], stddev[AAMODEL_WAG]);
				MrBayesPrintf (fp, "%s   Rtrev        %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_RTREV], stddev[AAMODEL_RTREV]);
				MrBayesPrintf (fp, "%s   Cprev        %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_CPREV], stddev[AAMODEL_CPREV]);
				MrBayesPrintf (fp, "%s   Vt           %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_VT], stddev[AAMODEL_VT]);
				MrBayesPrintf (fp, "%s   Blosum       %1.3lf       %1.6lf\n", spacer, prob[AAMODEL_BLOSUM], stddev[AAMODEL_BLOSUM]);
				MrBayesPrintf (fp, "%s   ----------------------------------\n", spacer);
				}
			}
		}

	return (NO_ERROR);
}





/* PrintOverlayPlot: Print overlay x-y plot of log likelihood vs. generation for several runs */
int PrintOverlayPlot (FILE *fp)
{
	int		i, j, k, k2, n, len, screenHeight, screenWidth, numY[60], whichIsX, whichIsY;
	char	temp[100], plotSymbol[15][60];
	MrBFlt	x, y, minX, maxX, minY, maxY, meanY[60];
	
	if (sumpParams.numRuns == 2)
		MrBayesPrintf (fp, "\n%s   Overlay plot for both runs:\n", spacer);
	else
		MrBayesPrintf (fp, "\n%s   Overlay plot for all %d runs:\n", spacer, sumpParams.numRuns);
	if (sumpParams.numRuns > 9)
		MrBayesPrintf (fp, "%s   (1 = Run number 1; 2 = Run number 2 etc.; x = Run number 10 or above; * = Several runs)\n", spacer);
	else if (sumpParams.numRuns > 2)
		MrBayesPrintf (fp, "%s   (1 = Run number 1; 2 = Run number 2 etc.; * = Several runs)\n", spacer);
	else
		MrBayesPrintf (fp, "%s   (1 = Run number 1; 2 = Run number 2; * = Both runs)\n", spacer);

	/* print overlay x-y plot of log likelihood vs. generation for all runs */
	screenWidth = 60; /* don't change this without changing numY, meanY, and plotSymbol declared above */
	screenHeight = 15;

	/* find the x and y columns of the first run */
	whichIsX = whichIsY = -1;
	for (i=0; GetNameFromString (headerNames, temp, i+1) != ERROR; i++)
		{
		len = (int) strlen(temp);
		if (IsSame (temp, "Gen") == SAME)
			whichIsX = i;
		if (IsSame (temp, "lnL") == SAME)
			whichIsY = i;
		}		
	if (whichIsX < 0 || whichIsY < 0)
		{
		MrBayesPrint ("%s   Could not find the 'Gen' and 'lnL' columns of the first file\n", spacer);
		return (ERROR);
		}

	/* find minX, minY, maxX, and maxY over all runs */
	minX = minY = 1000000000.0;
	maxX = maxY = -1000000000.0;
	for (n=0; n<sumpParams.numRuns; n++)
		{
		for (i=0; i<numRows; i++)
			{
			x = parameterValues[(n*numRows + i) * numColumns + whichIsX];
			if (x < minX)
				minX = x;
			if (x > maxX)
				maxX = x;
			}
		}
	for (n=0; n<sumpParams.numRuns; n++)
		{
		y = 0.0;
		j = 0;
		k2 = 0;
		for (i=0; i<numRows; i++)
			{
			x = parameterValues[(n*numRows + i) * numColumns + whichIsX];
			k = (int) (((x - minX) / (maxX - minX)) * screenWidth);
			if (k <= 0)
				k = 0;
			if (k >= screenWidth)
				k = screenWidth - 1;
			if (k == j)
				{
				y += parameterValues[(n*numRows + i) * numColumns + whichIsY];
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
				y = parameterValues[(n*numRows + i) * numColumns + whichIsY];
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
	for (n=0; n<sumpParams.numRuns; n++)
		{
		whichIsX = whichIsY = -1;
		for (i=0; GetNameFromString (headerNames, temp, i+1) != ERROR; i++)
			{
			len = (int) strlen(temp);
			if (IsSame (temp, "Gen") == SAME)
				whichIsX = i;
			else if (IsSame (temp, "lnL") == SAME)
				whichIsY = i;
			}		
		if (whichIsX < 0 || whichIsY < 0)
			{
			MrBayesPrint ("%s   Could not find the 'Gen' and 'lnL' columns of file %d\n", spacer, n+1);
			return (ERROR);
			}

		/* find the plot points for this run */
		for (i=0; i<screenWidth; i++)
			{
			numY[i] = 0;
			meanY[i] = 0.0;
			}
		for (i=0; i<numRows; i++)
			{
			x = parameterValues[(n*numRows + i) * numColumns + whichIsX];
			y = parameterValues[(n*numRows + i) * numColumns + whichIsY];
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
	MrBayesPrintf (fp, "\n   +");
	for (i=0; i<screenWidth; i++)
		MrBayesPrintf (fp, "-");
	MrBayesPrintf (fp, "+ %1.2lf\n", maxY);
	for (i=screenHeight-1; i>=0; i--)
		{
		MrBayesPrintf (fp, "   |");
		for (j=0; j<screenWidth; j++)
			{
			MrBayesPrintf (fp, "%c", plotSymbol[i][j]);
			}
		MrBayesPrintf (fp, "|\n");
		}
	MrBayesPrintf (fp, "   +");
	for (i=0; i<screenWidth; i++)
		{
		if (i % (screenWidth/10) == 0 && i != 0)
			MrBayesPrintf (fp, "+");
		else
			MrBayesPrintf (fp, "-");
		}
	MrBayesPrintf (fp, "+ %1.2lf\n", minY);
	MrBayesPrintf (fp, "   ^");
	for (i=0; i<screenWidth; i++)
		MrBayesPrintf (fp, " ");
	MrBayesPrintf (fp, "^\n");
	MrBayesPrintf (fp, "   %1.0lf", minX);
	for (i=0; i<screenWidth-(int)(log10(minX)); i++)
		MrBayesPrintf (fp, " ");
	MrBayesPrintf (fp, "%1.0lf\n\n", maxX);

	return (NO_ERROR);
}




/* PrintParamTable: Print parameter table (not amino acids) to file fp */
int PrintParamTable (FILE *fp, MrBFlt *parmVals, int nRows, int nCols, int nRuns)
{
	int		i, j, len, longestHeader;
	char	temp[100];
	MrBFlt	lower, upper, median, mean, var, sqrt_R=0.0;
	/* PAUL 20051206 We need to keep the original values from parmVals otherwise a second call
	 * to this function will give an incorrect PSRF. This is the case with 
	 * printtofile=yes. This might be a bug in the PotentialScaleReduction function? */
	MrBFlt  *parmVals2;	
	
	/* calculate longest header */
	longestHeader = 9;	/* length of 'parameter' */
	for (i=0; GetNameFromString (headerNames, temp, i+1) != ERROR; i++)
		{
		len = (int) strlen(temp);
		if (IsSame (temp, "Gen") == SAME)
			continue;
		if (IsSame (temp, "Aamodel") == SAME || IsSame (temp, "Aamodel") == CONSISTENT_WITH)
			continue;
		if (IsSame (temp, "lnL") == SAME)
			continue;
		if (len > longestHeader)
			longestHeader = len;
		}
	
	/* print the header rows */
	MrBayesPrintf (fp, "%s   %*c                                  95%% Cred. Interval\n", spacer, longestHeader, ' ');
	MrBayesPrintf (fp, "%s   %*c                                 ----------------------\n", spacer, longestHeader, ' ');

	MrBayesPrintf (fp, "%s   Parameter%*c      Mean        Variance       Lower         Upper         Median", spacer, longestHeader-9, ' ');
	if (nRuns > 1)
		MrBayesPrintf (fp, "       PSRF *");
	MrBayesPrintf (fp, "\n");

	MrBayesPrintf (fp, "%s   ", spacer);
	for (j=0; j<longestHeader+1; j++)
		MrBayesPrintf (fp, "-");
	MrBayesPrintf (fp, "--------------------------------------------------------------------");
	if (nRuns > 1)
		MrBayesPrintf (fp, "-----------");
	MrBayesPrintf (fp, "\n");

	/* print table values */
	for (i=0; GetNameFromString (headerNames, temp, i+1) != ERROR; i++)
		{
		len = (int) strlen(temp);
		if (IsSame (temp, "Gen") == SAME)
			continue;
		if (IsSame (temp, "Aamodel") == SAME || IsSame (temp, "Aamodel") == CONSISTENT_WITH)
			continue;
		if (IsSame (temp, "lnL") == SAME)
			continue;

		if (nRuns > 1)
			sqrt_R = PotentialScaleReduction (parmVals, nRows, nCols, nRuns, i);

		parmVals2 = (MrBFlt *)SafeMalloc((size_t) (nRuns * numRows * numColumns * sizeof(MrBFlt)));
        if (!parmVals2)
            {
            MrBayesPrint ("%s   Problem allocating parameterValues\n", spacer);
            return ERROR;
            }
        memcpy(parmVals2, parmVals, (nRuns * numRows * numColumns * sizeof(MrBFlt)));
		if (SortColumn (parmVals2, nRuns*nRows, nCols, i) == ERROR) 
		    {
		    free(parmVals2);
			return ERROR;
		    }
		    
		MeanVariance (parmVals2, nRuns*nRows, nCols, i, &mean, &var);
		lower = parmVals2[((int)(nRuns * nRows * 0.025) * nCols + i)];
		upper = parmVals2[((int)(nRuns * nRows * 0.975) * nCols + i)];
		median = parmVals2[((int)(nRuns * nRows * 0.500) * nCols + i)];
		free(parmVals2);
		
		MrBayesPrintf (fp, "%s   %-*s ", spacer, longestHeader, temp);
		MrBayesPrintf (fp, "%12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf", mean, var, lower, upper, median);
		if (nRuns > 1)
			{
			if (sqrt_R < 0.0)
				MrBayesPrintf (fp, "       N/A  ");
			else
				MrBayesPrintf (fp, "  %9.3lf", sqrt_R);
			}
		MrBayesPrintf (fp, "\n");
		}
	MrBayesPrintf (fp, "%s   ", spacer);
	for (j=0; j<longestHeader+1; j++)
		MrBayesPrintf (fp, "-");
	MrBayesPrintf (fp, "--------------------------------------------------------------------");
	if (nRuns > 1)
		MrBayesPrintf (fp, "-----------");
	MrBayesPrintf (fp, "\n");
	if (nRuns > 1)
		{
		MrBayesPrintf (fp, "%s   * Convergence diagnostic (PSRF = Potential scale reduction factor [Gelman\n", spacer);
		MrBayesPrintf (fp, "%s     and Rubin, 1992], uncorrected) should approach 1 as runs converge. The\n", spacer);
		MrBayesPrintf (fp, "%s     values may be unreliable if you have a small number of samples. PSRF should\n", spacer);
		MrBayesPrintf (fp, "%s     only be used as a rough guide to convergence since all the assumptions\n", spacer);
		MrBayesPrintf (fp, "%s     that allow one to interpret it as a scale reduction factor are not met in\n", spacer);
		MrBayesPrintf (fp, "%s     the phylogenetic context.\n", spacer);
		}
	return (NO_ERROR);
}




/* PrintPlot: Print x-y plot of log likelihood vs. generation */
int PrintPlot (FILE *fp, MrBFlt *paramVals)
{
	int		i, j, k, len, numY[60], screenWidth, screenHeight, whichIsX, whichIsY;
	char	temp[100];
	MrBFlt	x, y, minX, maxX, minY, maxY, meanY[60];
	
	/* print x-y plot of log likelihood vs. generation */
	screenWidth = 60; /* don't change this without changing numY and meanY, declared above */
	screenHeight = 15;
	whichIsX = whichIsY = -1;
	for (i=0; GetNameFromString (headerNames, temp, i+1) != ERROR; i++)
		{
		len = (int) strlen(temp);
		if (IsSame (temp, "Gen") == SAME)
			whichIsX = i;
		else if (IsSame (temp, "lnL") == SAME)
			whichIsY = i;
		}		
	if (whichIsX < 0 || whichIsY < 0)
		{
		MrBayesPrint ("%s   Could not find the 'Gen' and 'lnL' columns\n", spacer);
		return (ERROR);
		}

	/* find minX and maxX */
	minX = 1000000000.0;
	maxX = -1000000000.0;
	for (i=0; i<numRows; i++)
		{
		x = paramVals[i * numColumns + whichIsX];
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
	for (i=0; i<numRows; i++)
		{
		x = paramVals[i * numColumns + whichIsX];
		y = paramVals[i * numColumns + whichIsY];
		k = (int)(((x - minX) / (maxX - minX)) * screenWidth);
		if (k >= screenWidth)
			k = screenWidth - 1;
		if (k < 0)
			k = 0;
		meanY[k] += y;
		numY[k]++;
		}

	/* find minY and maxY */
	minY = 1000000000.0;
	maxY = -1000000000.0;
	for (i=0; i<screenWidth; i++)
		{
		meanY[i] /= numY[i];
		if (meanY[i] < minY)
			minY = meanY[i];
		if (meanY[i] > maxY)
			maxY = meanY[i];
		}

	/* print plot */
	MrBayesPrintf (fp, "\n   +");
	for (i=0; i<screenWidth; i++)
		MrBayesPrintf (fp, "-");
	MrBayesPrintf (fp, "+ %1.2lf\n", maxY);
	for (j=screenHeight-1; j>=0; j--)
		{
		MrBayesPrintf (fp, "   |");
		for (i=0; i<screenWidth; i++)
			{
			if (numY[i] > 0)
				{
				k = (int) (((meanY[i] - minY) / (maxY - minY)) * screenHeight);
				if (k < 0)
					k = 0;
				if (k >= screenHeight)
					k = screenHeight - 1;
				if (k == j)
					MrBayesPrintf (fp, "*");
				else
					MrBayesPrintf (fp, " ");
				}
			else
				{
				MrBayesPrintf (fp, " ");
				}
			}
		MrBayesPrintf (fp, "|\n");
		}
	MrBayesPrintf (fp, "   +");
	for (i=0; i<screenWidth; i++)
		{
		if (i % (screenWidth/10) == 0 && i != 0)
			MrBayesPrintf (fp, "+");
		else
			MrBayesPrintf (fp, "-");
		}
	MrBayesPrintf (fp, "+ %1.2lf\n", minY);
	MrBayesPrintf (fp, "   ^");
	for (i=0; i<screenWidth; i++)
		MrBayesPrintf (fp, " ");
	MrBayesPrintf (fp, "^\n");
	MrBayesPrintf (fp, "   %1.0lf", minX);
	for (i=0; i<screenWidth-(int)(log10(minX)); i++)
		MrBayesPrintf (fp, " ");
	MrBayesPrintf (fp, "%1.0lf\n\n", maxX);

	return (NO_ERROR);
}





void PrintPlotHeader (FILE *fp)
{
	MrBayesPrintf (fp, "\n");
	if (sumpParams.numRuns > 1)
		{
		MrBayesPrintf (fp, "%s   Below are rough plots of the generation (x-axis) versus the log   \n", spacer);
		MrBayesPrintf (fp, "%s   probability of observing the data (y-axis). You can use these     \n", spacer);
		MrBayesPrintf (fp, "%s   graphs to determine what the burn in for your analysis should be. \n", spacer);
		MrBayesPrintf (fp, "%s   When the log probability starts to plateau you may be at station- \n", spacer);
		MrBayesPrintf (fp, "%s   arity. Sample trees and parameters after the log probability      \n", spacer);
		MrBayesPrintf (fp, "%s   plateaus. Of course, this is not a guarantee that you are at sta- \n", spacer);
		MrBayesPrintf (fp, "%s   tionarity. Also examine the convergence diagnostics provided by   \n", spacer);
		MrBayesPrintf (fp, "%s   the 'sump' and 'sumt' commands for all the parameters in your     \n", spacer);
		MrBayesPrintf (fp, "%s   model. Remember that the burn in is the number of samples to dis- \n", spacer);
		MrBayesPrintf (fp, "%s   card. There are a total of ngen / samplefreq samples taken during \n", spacer);
		MrBayesPrintf (fp, "%s   a MCMC analysis.                                                  \n", spacer);
		}
	else
		{
		MrBayesPrintf (fp, "%s   Below is a rough plot of the generation (x-axis) versus the log   \n", spacer);
		MrBayesPrintf (fp, "%s   probability of observing the data (y-axis). You can use this      \n", spacer);
		MrBayesPrintf (fp, "%s   graph to determine what the burn in for your analysis should be.  \n", spacer);
		MrBayesPrintf (fp, "%s   When the log probability starts to plateau you may be at station- \n", spacer);
		MrBayesPrintf (fp, "%s   arity. Sample trees and parameters after the log probability      \n", spacer);
		MrBayesPrintf (fp, "%s   plateaus. Of course, this is not a guarantee that you are at sta- \n", spacer);
		MrBayesPrintf (fp, "%s   analysis should be. When the log probability starts to plateau    \n", spacer);
		MrBayesPrintf (fp, "%s   tionarity. When possible, run multiple analyses starting from dif-\n", spacer);
		MrBayesPrintf (fp, "%s   ferent random trees; if the inferences you make for independent   \n", spacer);
		MrBayesPrintf (fp, "%s   analyses are the same, this is reasonable evidence that the chains\n", spacer);
		MrBayesPrintf (fp, "%s   have converged. You can use MrBayes to run several independent    \n", spacer);
		MrBayesPrintf (fp, "%s   analyses simultaneously. During such a run, MrBayes will monitor  \n", spacer);
		MrBayesPrintf (fp, "%s   the convergence of topologies. After the run has been completed,  \n", spacer);
		MrBayesPrintf (fp, "%s   the 'sumt' and 'sump' functions will provide additional conver-   \n", spacer);
		MrBayesPrintf (fp, "%s   gence diagnostics for all the parameters in your model. Remember  \n", spacer);
		MrBayesPrintf (fp, "%s   that the burn in is the number of samples to discard. There are   \n", spacer);
		MrBayesPrintf (fp, "%s   a total of ngen / samplefreq samples taken during a MCMC analysis.\n", spacer);
		}
}

void Sort (MrBFlt *item, int count)

{

	Sort2 (item, 0, count-1);

}





void Sort2 (MrBFlt *item, int left, int right)

{

	register int	i, j;
	MrBFlt			x, y;

	i = left;
	j = right;
	x = item[(left+right)/2];
	do 
		{
		while (item[i] < x && i < right)
			i++;
		while (x < item[j] && j > left)
			j--;
		if (i <= j)
			{
			y = item[i];
			item[i] = item[j];
			item[j] = y;
			
			i++;
			j--;
			}
		} while (i <= j);
	if (left < j)
		Sort2 (item, left, j);
	if (i < right)
		Sort2 (item, i, right);

}





int SortColumn (MrBFlt *paramValues, int nRows, int nCols, int column)

{

	int			i;
	MrBFlt		*tempParamValues;
	
	/* allocate information for sorting */
	tempParamValues = (MrBFlt *)SafeMalloc((size_t) (nRows * sizeof(MrBFlt)));
	if (!tempParamValues)
		{
		MrBayesPrint ("%s   Problem allocating tempParamValues\n", spacer);
		return ERROR;
		}
		
	for (i=0; i<nRows; i++)
		tempParamValues[i] = paramValues[i * nCols + column];
	
	Sort (tempParamValues, nRows);
	
	for (i=0; i<nRows; i++)
		paramValues[i * nCols + column] = tempParamValues[i];

#	if 0
	for (i=0; i<nRows; i++)
		printf ("%2.3lf ", paramValues[i * nCols + column]);
	printf ("\n");
#	endif

	/* free memory */
	free (tempParamValues);
	return (NO_ERROR);
	
}




int SortParameters (MrBFlt *paramValues)

{

	int			i, j;
	MrBFlt		*tempParamValues;
	
	/* allocate information for sorting */
	tempParamValues = (MrBFlt *)SafeMalloc((size_t) (numRows * sizeof(MrBFlt)));
	if (!tempParamValues)
		{
		MrBayesPrint ("%s   Problem allocating tempParamValues\n", spacer);
		goto errorExit;
		}
		
	/* sort each column */
	for (j=0; j<numColumns; j++)
		{
		for (i=0; i<numRows; i++)
			tempParamValues[i] = 0.0;
		for (i=0; i<numRows; i++)
			tempParamValues[i] = paramValues[i * numColumns + j];
		Sort (tempParamValues, numRows);
		for (i=0; i<numRows; i++)
			paramValues[i * numColumns + j] = tempParamValues[i];
		}

#	if 0
	for (i=0; i<numRows; i++)
		{
		for (j=0; j<numColumns; j++)
			printf ("%2.3lf ", parameterValues[i * numColumns + j]);
		printf ("\n");
		}
#	endif

	/* free memory */
	free (tempParamValues);
	return (NO_ERROR);
	
	errorExit:
		if (tempParamValues)
			free (tempParamValues);
		return (ERROR);

}








