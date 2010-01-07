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
 */
/* id-string for ident, do not edit: cvs will update this string */
const char plotID[]="$Id: plot.c,v 3.20 2009/01/06 21:40:10 ronquist Exp $";

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
#include "plot.h"
#include "sump.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif


/* local (to this file) */
int		foundCurly;
char	*plotTokenP, plotToken[CMD_STRING_LENGTH];






int DoPlot (void)

{


	int			i, j, k, n, lineTerm, longestLineLength, tokenType, lineNum, lastTokenWasDash,
				inPlotComment, allDigitLine, nNumbersOnThisLine, lastNonDigitLine,
				numParamLines, numLinesToRead, numLinesRead, firstNumCols=0, nLines,
				nHeaders, len, longestHeader, whichIsX, whichIsY, screenWidth, screenHeigth, numY[60],
				numPlotted;
	MrBFlt		tempD, minX, minY, maxX, maxY,
				meanY[60], x, y, diff;
	char		*s=NULL, *headerLine=NULL, temp[100];
	FILE		*fp;
	
#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
#	endif

	/* set file pointer to NULL */
	fp = NULL;

	/* tell user we are ready to go */
	MrBayesPrint ("%s   Plotting parameters in file %s\n", spacer, plotParams.plotFileName);
	
	/* open binary file */
	if ((fp = OpenBinaryFileR(plotParams.plotFileName)) == NULL)
		goto errorExit;
		
	/* find out what type of line termination is used */
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
		
	/* find length of longest line */
	longestLineLength = LongestLine (fp);
	MrBayesPrint ("%s   Longest line length = %d\n", spacer, longestLineLength);
	longestLineLength += 10;
	
	/* allocate a string long enough to hold a line */
	if (memAllocs[ALLOC_SUMPSTRING] == YES)
		{
		MrBayesPrint ("%s   Plot string is already allocated\n", spacer);
		goto errorExit;
		}
	s = (char *)SafeMalloc((size_t) (longestLineLength * sizeof(char)));
	if (!s)
		{
		MrBayesPrint ("%s   Problem allocating string for reading plot file\n", spacer);
		goto errorExit;
		}
	headerLine = (char *)SafeMalloc((size_t) (longestLineLength * sizeof(char)));
	if (!headerLine)
		{
		MrBayesPrint ("%s   Problem allocating headerLine for reading plot file\n", spacer);
		goto errorExit;
		}
	headerNames = (char *)SafeMalloc((size_t) ((longestLineLength+40) * sizeof(char)));
	if (!headerNames)
		{
		MrBayesPrint ("%s   Problem allocating headerNames for reading plot file\n", spacer);
		goto errorExit;
		}
	for (i=0; i<longestLineLength+40; i++)
		headerNames[i] = ' ';
	headerNames[longestLineLength+40-1] = '\0';
	memAllocs[ALLOC_SUMPSTRING] = YES;
		
	/* close binary file */
	SafeFclose (&fp);
	
	/* open text file */
	if ((fp = OpenTextFileR(plotParams.plotFileName)) == NULL)
		goto errorExit;
	
	/* Check file for appropriate blocks. We want to find the last block
	   in the file and start from there. */
	inPlotComment = NO;
	lineNum = lastNonDigitLine = numParamLines = 0;
	while (fgets (s, longestLineLength, fp) != NULL)
		{
		plotTokenP = &s[0];
		allDigitLine = YES;
		lastTokenWasDash = NO;
		nNumbersOnThisLine = 0;
		do
			{
			GetSumpToken (&tokenType, &plotTokenP, plotToken);
			/*printf ("%s (%d)\n", plotToken, tokenType);*/
			if (IsSame("[", plotToken) == SAME)
				inPlotComment = YES;
			if (IsSame("]", plotToken) == SAME)
				inPlotComment = NO;
				
			if (inPlotComment == NO)
				{
				if (tokenType == NUMBER)
					{
					sscanf (plotToken, "%lf", &tempD);
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
				
			} while (*plotToken);
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
	if (inPlotComment == YES)
		{
		MrBayesPrint ("%s   Unterminated comment in file %s\n", spacer, plotParams.plotFileName);
		goto errorExit;
		}
	if (numParamLines <= 0)
		{
		MrBayesPrint ("%s   No parameters were found in file %s\n", spacer, plotParams.plotFileName);
		goto errorExit;
		}
	if (plotParams.plotBurnIn > numParamLines)
		{
		MrBayesPrint ("%s   No parameters are sampled as the burnin exceeds the number of lines in last block\n", spacer);
		MrBayesPrint ("%s   Try setting burnin to a number less than %d\n", spacer, numParamLines);
		goto errorExit;
		}
		
	/* tell the user that everything is fine */
	MrBayesPrint ("%s   Found %d parameter lines in file \"%s\"\n", spacer, numParamLines, plotParams.plotFileName);
	if (plotParams.plotBurnIn > 0)
		MrBayesPrint ("%s   Of the %d lines, %d of them will be summarized (starting at line %d)\n", spacer, numParamLines, numParamLines - plotParams.plotBurnIn, lastNonDigitLine + plotParams.plotBurnIn + 1);
	else
		MrBayesPrint ("%s   All %d lines will be summarized (starting at line %d)\n", spacer, numParamLines, lastNonDigitLine+1);
	MrBayesPrint ("%s   (Only the last set of lines will be read, in case multiple\n", spacer);
	MrBayesPrint ("%s   parameter blocks are present in the same file.)\n", spacer);
	
	/* Calculate and check the number of columns and rows for the file */
	(void)fseek(fp, 0L, 0);	
	for (lineNum=0; lineNum<lastNonDigitLine+plotParams.plotBurnIn; lineNum++)
		fgets (s, longestLineLength, fp);
	inPlotComment = NO;
	nLines = 0;
	numRows = numColumns = 0;
	while (fgets (s, longestLineLength, fp) != NULL)
		{
		plotTokenP = &s[0];
		allDigitLine = YES;
		lastTokenWasDash = NO;
		nNumbersOnThisLine = 0;
		do
			{
			GetSumpToken (&tokenType, &plotTokenP, plotToken);
			if (IsSame("[", plotToken) == SAME)
				inPlotComment = YES;
			if (IsSame("]", plotToken) == SAME)
				inPlotComment = NO;
			if (inPlotComment == NO)
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
			} while (*plotToken);
		lineNum++;
		if (allDigitLine == NO)
			{
			MrBayesPrint ("%s   Found a line with non-digit characters (line %d)\n", spacer, lineNum);
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
						MrBayesPrint ("%s   Number of lines is not even (%d in first line and %d in %d line)\n", spacer, firstNumCols, nNumbersOnThisLine, lineNum);
						goto errorExit;
						}
					}
				}
			}
		}
	numRows = nLines;
	numColumns = firstNumCols;
	MrBayesPrint ("%s   %d rows and %d columns in each row\n", spacer, numRows, numColumns);
	
	/* allocate space to hold parameter information */
	if (numRows == 0 || numColumns == 0)
		{
		MrBayesPrint ("%s   The number of rows or columns is equal to zero\n", spacer);
		goto errorExit;
		}
	if (memAllocs[ALLOC_SUMPINFO] == YES)
		{
		MrBayesPrint ("%s   Plot string is already allocated\n", spacer);
		goto errorExit;
		}
	parameterValues = (MrBFlt *)SafeMalloc((size_t) (numRows * numColumns * sizeof(MrBFlt)));
	if (!parameterValues)
		{
		MrBayesPrint ("%s   Problem allocating parameterValues\n", spacer);
		goto errorExit;
		}
	memAllocs[ALLOC_SUMPINFO] = YES;
	for (i=0; i<numRows*numColumns; i++)
		parameterValues[i] = 0.0;

	/* Now we read the file for real. First, rewind file pointer to beginning of file... */
	(void)fseek(fp, 0L, 0);	
	
	/* ...and fast forward to beginning of last unburned parameter line. */
	for (lineNum=0; lineNum<lastNonDigitLine+plotParams.plotBurnIn; lineNum++)
		fgets (s, longestLineLength, fp);
		
	/* ...and parse file, line-by-line. We are only parsing lines that have digits that should be read. */
	inPlotComment = NO;
	numLinesToRead = numParamLines - plotParams.plotBurnIn;
	numLinesRead = j = 0;
	while (fgets (s, longestLineLength, fp) != NULL)
		{
		plotTokenP = &s[0];
		allDigitLine = YES;
		lastTokenWasDash = NO;
		nNumbersOnThisLine = 0;
		do
			{
			GetSumpToken (&tokenType, &plotTokenP, plotToken);
			if (IsSame("[", plotToken) == SAME)
				inPlotComment = YES;
			if (IsSame("]", plotToken) == SAME)
				inPlotComment = NO;
			if (inPlotComment == NO)
				{
				if (tokenType == NUMBER)
					{
					/* read the information from this line */
					if (j >= numRows * numColumns)
						{
						MrBayesPrint ("%s   Too many parameter values read in (%d)\n", spacer, j);
						goto errorExit;
						}
					sscanf (plotToken, "%lf", &tempD);
					if (lastTokenWasDash == YES)
						tempD *= -1.0;
					parameterValues[numLinesRead * numColumns + nNumbersOnThisLine] = tempD;
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
					MrBayesPrint ("%s   Found a line with non-digit characters (line %d)\n", spacer, lineNum);
					goto errorExit;
					}
				}
			} while (*plotToken);
		lineNum++;
		if (nNumbersOnThisLine > 0)
			numLinesRead++;
		}
	
	/* tell user how many lines were successfully read */
	MrBayesPrint ("%s   Successfully read %d lines from last parameter block\n", spacer, numLinesRead);
	
	/* Check that at least one parameter line was read in. */
	if (numLinesRead <= 0)
		{
		MrBayesPrint ("%s   No parameters read in\n", spacer);
		goto errorExit;
		}
				
	/* separate header line into titles for each column */
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
		
	/* print x-y plot of parameter vs. generation */
	screenWidth = 60; /* don't change this without changing numY and meanY, declared above */
	screenHeigth = 15;
	whichIsX = -1;
	for (i=0; i<nHeaders; i++)
		{
		if (GetNameFromString (headerNames, temp, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Error getting header names \n", spacer);
			goto errorExit;
			}
		len = (int) strlen(temp);
		if (IsSame (temp, "Gen") == SAME)
			whichIsX = i;
		}		
		
	if (whichIsX < 0)
		{
		MrBayesPrint ("%s   Could not find a column labelled \"Gen\" \n", spacer);
		goto errorExit;
		}
		
	numPlotted = 0;
	for (n=0; n<nHeaders; n++)
		{
		if (GetNameFromString (headerNames, temp, n+1) == ERROR)
			{
			MrBayesPrint ("%s   Error getting header names \n", spacer);
			goto errorExit;
			}
		whichIsY = -1;
		if (!strcmp(plotParams.match, "Perfect"))
			{
			if (IsSame (temp, plotParams.parameter) == SAME)
				whichIsY = n;
			}
		else if (!strcmp(plotParams.match, "All"))
			{
			whichIsY = n;
			}
		else
			{
			if (IsSame (temp, plotParams.parameter) == CONSISTENT_WITH)
				whichIsY = n;
			}
			
		if (whichIsY >= 0 && whichIsX != whichIsY)
			{			
			minX = minY = 1000000000.0;
			maxX = maxY = -1000000000.0;
			for (i=0; i<numRows; i++)
				{
				x = parameterValues[i * numColumns + whichIsX];
				y = parameterValues[i * numColumns + whichIsY];
				if (x < minX)
					minX = x;
				if (y < minY)
					minY = y;
				if (x > maxX)
					maxX = x;
				if (y > maxY)
					maxY = y;
				}
			for (i=0; i<screenWidth; i++)
				{
				numY[i] = 0;
				meanY[i] = 0.0;
				}
			for (i=0; i<numRows; i++)
				{
				x = parameterValues[i * numColumns + whichIsX];
				y = parameterValues[i * numColumns + whichIsY];
				k = (int)(((x - minX) / (maxX - minX)) * screenWidth);
				if (k >= screenWidth)
					k = screenWidth - 1;
				meanY[k] += y;
				numY[k]++;
				}
			if (maxY - minY < 0.000001)
				{
				maxY = meanY[0]/numY[0] + (MrBFlt) 0.1;
				minY = meanY[0]/numY[0] - (MrBFlt) 0.1;
				} 
			else
				{
				diff = maxY - minY;
				maxY += diff * (MrBFlt) 0.025;
				minY -= diff * (MrBFlt) 0.025;
				}
			MrBayesPrint ("\n");
			MrBayesPrint ("%s   Rough plot of parameter %s \n", spacer, temp);
				
			MrBayesPrint ("\n   +");
			for (i=0; i<screenWidth; i++)
				MrBayesPrint ("-");
			MrBayesPrint ("+ %1.2lf\n", maxY);
			for (j=screenHeigth-1; j>=0; j--)
				{
				MrBayesPrint ("   |");
				for (i=0; i<screenWidth; i++)
					{
					if (numY[i] > 0)
						{
						if (meanY[i] / numY[i] > (((maxY - minY)/screenHeigth)*j)+minY && meanY[i] / numY[i] <= (((maxY - minY)/screenHeigth)*(j+1))+minY)
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
			MrBayesPrint ("+ %1.2lf\n", minY);
			MrBayesPrint ("   ^");
			for (i=0; i<screenWidth; i++)
				MrBayesPrint (" ");
			MrBayesPrint ("^\n");
			MrBayesPrint ("   %1.0lf", minX);
			for (i=0; i<screenWidth; i++)
				MrBayesPrint (" ");
			MrBayesPrint ("%1.0lf\n\n", maxX);
			numPlotted++;
			}
				
		}
		
	if (numPlotted == 0)
		{
		MrBayesPrint ("%s   Did not find any parameters matching \"%s\" to plot\n", spacer, plotParams.parameter);
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
	expecting = Expecting(COMMAND);
	
#	if defined (MPI_ENABLED)
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
		strcpy (spacer, "");
		strcpy (plotToken, "Plot");
		i = 0;
		if (FindValidCommand (plotToken, &i) == ERROR)
			MrBayesPrint ("%s   Could not find plot\n", spacer);
		return (ERROR);	
	
}





int DoPlotParm (char *parmName, char *tkn)

{

	int			tempI;
	char		tempStr[100];

	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before sumt can be used\n", spacer);
		return (ERROR);
		}

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
		/* set Filename (plotParams.plotFileName) ***************************************************/
		else if (!strcmp(parmName, "Filename"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
				strcpy (plotParams.plotFileName, tkn);
				MrBayesPrint ("%s   Setting plot filename to %s\n", spacer, plotParams.plotFileName);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Burnin (plotParams.plotBurnIn) *******************************************************/
		else if (!strcmp(parmName, "Burnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
				plotParams.plotBurnIn = tempI;
				MrBayesPrint ("%s   Setting plot burnin to %d\n", spacer, plotParams.plotBurnIn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Parameter (plotParams.parameter) *******************************************************/
		else if (!strcmp(parmName, "Parameter"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
				strcpy (plotParams.parameter, tkn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Parameter (plotParams.match) *******************************************************/
		else if (!strcmp(parmName, "Match"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					strcpy (plotParams.match, tempStr);
				else
					return (ERROR);
				
				MrBayesPrint ("%s   Setting plot matching to %s\n", spacer, plotParams.match);
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

