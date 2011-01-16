/*
 *  MrBayes 3
 *
 *  (c) 2002-2010
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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "mb.h"
#include "globals.h"
#include "utils.h"

const char* const svnRevisionUtilsC="$Rev$";   /* Revision keyword which is expended/updated by svn on each commit/update*/

/* AddBitfield: Add bitfield to list of bitfields */
int AddBitfield (SafeLong ***list, int listLen, int *set, int setLen)
{
    int     i, nLongsNeeded;

    nLongsNeeded = (setLen - 1) / nBitsInALong + 1;

    (*list) = (SafeLong **) SafeRealloc ((void *)(*list), (size_t)((listLen+1)*sizeof(SafeLong *)));
    if (!(*list))
        return ERROR;
    
    (*list)[listLen] = (SafeLong *) SafeMalloc ((size_t)(nLongsNeeded*sizeof(SafeLong)));
    if (!(*list)[listLen])
        return ERROR;

    ClearBits ((*list)[listLen], nLongsNeeded);
    for (i=0; i<setLen; i++)
        if (set[i] == YES)
            SetBit(i, (*list)[listLen]);

    return NO_ERROR;
}





#if defined (SSE_ENABLED)
/* Aligned safe free */
void AlignedSafeFree (void **ptr)
{

    ALIGNED_FREE (*ptr);

    (*ptr) = NULL;
}
#endif





/* ClearBit: Clear one bit in a bitfield */
void ClearBit (int i, SafeLong *bits)
{
	SafeLong		x;

	bits += i / nBitsInALong;

	x = 1 << (i % nBitsInALong);
    x ^= safeLongWithAllBitsSet;

	(*bits) &= x;
}





/* ClearBits: Clear all bits in a bitfield */
void ClearBits (SafeLong *bits, int nLongs)
{
	int     i;
    
    for (i=0; i<nLongs; i++)
        bits[i] = 0;
}





/* CopyResults: copy results from one file to another*/
int CopyResults (FILE *toFile, char *fromFileName, int lastGen)
{
    int     longestLine;
    char    *strBuf, *strCpy, *word;
    FILE    *fromFile;

    if ((fromFile = OpenBinaryFileR(fromFileName)) == NULL)
        return ERROR;

    longestLine = LongestLine(fromFile)+10;
    SafeFclose(&fromFile);
    strBuf = (char *) calloc (2*(longestLine+2),sizeof(char));
    strCpy = strBuf + longestLine + 2;

    if ((fromFile = OpenTextFileR(fromFileName)) == NULL)
        return ERROR;
    
    while (fgets(strBuf,longestLine,fromFile)!=NULL)
        {
        strncpy(strCpy,strBuf,longestLine);
        word = strtok(strCpy," ");
        /* atoi returns 0 when word is not integer number */
        if (atoi(word)>lastGen)
            break;
        fprintf (toFile,"%s",strBuf);
        fflush (toFile);
        }
    
    SafeFclose(&fromFile);
    free(strBuf);
    return (NO_ERROR);
}




/* CopyTreeResults: copy tree results from one file to another*/
int CopyTreeResults (FILE *toFile, char *fromFileName, int lastGen, int *numTrees)
{
    int     longestLine;
    char    *strBuf, *strCpy, *word;
    FILE    *fromFile;
    
    (*numTrees) = 0;

    if ((fromFile = OpenBinaryFileR(fromFileName)) == NULL)
        return ERROR;

    longestLine = LongestLine(fromFile)+10;
    SafeFclose(&fromFile);
    strBuf = (char *) calloc (2*(longestLine+2),sizeof(char));
    strCpy = strBuf + longestLine + 2;

    if ((fromFile = OpenTextFileR(fromFileName)) == NULL)
        return ERROR;
    
    while (fgets(strBuf,longestLine,fromFile)!=NULL)
        {
        strncpy(strCpy,strBuf,longestLine);
        word = strtok(strCpy," ");
        if (strcmp(word,"tree")==0)
            {
            word = strtok(NULL," ");
            /* atoi returns 0 when word is not integer number,
               4 is offset to get rid of "rep." in tree name */
            if (atoi(word+4)>lastGen)
                break;
            (*numTrees)++;
            fprintf (toFile,"%s",strBuf);
            }
        else if (*numTrees == 0)   /* do not print the end statement */
            fprintf (toFile,"%s",strBuf);
        fflush (toFile);
        }
        
    SafeFclose(&fromFile);
    free(strBuf);
    return (NO_ERROR);
}





/* FirstTaxonInPartition: Find index of first taxon in partition */
int FirstTaxonInPartition (SafeLong *partition, int length)

{

    int         i, j, nBits, taxon;
    SafeLong    x;

    nBits = sizeof(SafeLong) * 8;

    taxon = 0;
    for (i=0; i<length; i++)
		{
		x = 1;
		for (j=0; j<nBits; j++)
			{
			if (partition[i] & x)
				return taxon;
			taxon++;
			x <<= 1;
			}
		}

	return taxon;

}





/* FirstTree: Return file position of first tree after current position */
SafeLong FirstTree (FILE *fp, char *lineBuf, int longestLine)
{
	SafeLong	firstTree;
	char		*word;
	
	do {
		firstTree = ftell(fp);
		if ((fgets (lineBuf, longestLine, fp)) == NULL)
			return 0;
		word = strtok (lineBuf, " ");
		} while (strcmp(word,"tree")!=0);

	return (firstTree);
}





int Flip01 (int x)

{

	if (x == 0)
		return (1);
	else
		return (0);
		
}





void FlipBits (SafeLong *partition, int length, SafeLong *mask)

{

	int			i;
	
	for (i=0; i<length; i++)
		{
		partition[i] ^= mask[i];
		}
}





/*-----------------------------------------------------------------
|
|	FlipOneBit: flip bit n in SafeLong *p
|
------------------------------------------------------------------*/
void FlipOneBit (int n, SafeLong *p)

{

	SafeLong		x;

	p += n/nBitsInALong;
	x = 1 << (n % nBitsInALong);
	(*p) ^= x;

}





/* Convert from 0-based growth function over six states to model index */
int      FromGrowthFxnToIndex(int *growthFxn)
{
    /* TODO: Bret - silly return for now */
    return 0;
}





/* Convert from model index to 0-based growth function over six states */
void     FromIndexToGrowthFxn(int index, int *growthFxn)
{

    /* TODO: Bret - silly stuff for now */
    growthFxn[0] = 0;
    growthFxn[1] = 0;
    growthFxn[2] = 0;
    growthFxn[3] = 0;
    growthFxn[4] = 0;
    growthFxn[5] = 0;

}





/* GetIntSummary: Get summary statistics for a number of runs (int version) */
void GetIntSummary (int **vals, int nRows, int *rowCount, Stat *theStats, int HPD)

{

    int     i, j, nVals;
    MrBFlt  *theValues, *p;

    nVals = 0;
    for (i=0; i<nRows; i++)
        nVals += rowCount[i];

    theValues = (MrBFlt *) calloc (nVals, sizeof(MrBFlt));

    /* extract values */
    p = theValues;
    for (i=0; i<nRows; i++)
        {
        for (j=0; j<rowCount[i]; j++)
            {
            (*p++) = (MrBFlt) (vals[i][j]);
            }
        }
    
    /* get statistics */
    MeanVariance (theValues, nVals, &(theStats->mean), &(theStats->var));
    if (HPD == YES)
        LowerUpperMedian (theValues, nVals, &(theStats->lower), &(theStats->upper), &(theStats->median));
    else
        LowerUpperMedian (theValues, nVals, &(theStats->lower), &(theStats->upper), &(theStats->median));

    free (theValues);
}





/* Get k from 0-based growth function */
int GetKFromGrowthFxn(int *growthFxn)
{
    int i, k=0;

    for (i=0; i<6; i++)
        if (growthFxn[i] > k)
            k = growthFxn[i];
    
    return k+1;
}





/* GetSummary: Get summary statistics for a number of runs */
void GetSummary (MrBFlt **vals, int nRows, int *rowCount, Stat *theStats, int HPD)

{

    int     i, nVals;
    MrBFlt  *theValues, *p;

    nVals = 0;
    for (i=0; i<nRows; i++)
        nVals += rowCount[i];

    theValues = (MrBFlt *) calloc (nVals, sizeof(MrBFlt));

    /* extract values */
    p = theValues;
    for (i=0; i<nRows; i++)
        {
        memcpy((void *)(p),(void *)(vals[i]),(size_t)(rowCount[i]*sizeof(MrBFlt)));
        p += rowCount[i];
        }
    
    /* get statistics */
    MeanVariance (theValues, nVals, &(theStats->mean), &(theStats->var));
    if (HPD == YES)
        LowerUpperMedianHPD (theValues, nVals, &(theStats->lower), &(theStats->upper), &(theStats->median));
    else
        LowerUpperMedian (theValues, nVals, &(theStats->lower), &(theStats->upper), &(theStats->median));
    if (nRows > 1)
        theStats->PSRF = PotentialScaleReduction (vals, nRows, rowCount);

    free (theValues);
}





/* HarmonicArithmeticMean: Calculate harmonic and arithmetic mean from log values */
int HarmonicArithmeticMeanOnLogs (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *harm_mean)
{
	int				i, reliable;
	MrBFlt			a, aOld, x, y, scaler, n;

	reliable = YES;
	
	scaler = vals[nVals-1];
	a = aOld = n = 0.0;
	for (i=0; i<nVals; i++)
		{
		y = vals[i];
		y -= scaler;
		if (y < -100.0 || y > 100.0)
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
	
	scaler = (MrBFlt) (0.0 - vals[nVals-1]);
	a = aOld = n = 0.0;
	for (i=0; i<nVals; i++)
		{
		y = (MrBFlt) (0.0 - vals[i]);
		y -= scaler;
		if (y < -100.0 || y > 100.0)
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





/* IsBitSet: Is bit i set in SafeLong *bits ? */
int IsBitSet (int i, SafeLong *bits)

{

	SafeLong		x;

	bits += i / nBitsInALong;

	x = 1 << (i % nBitsInALong);

	if ((*bits) & x)
		return (YES);
	else
		return (NO);
		
}





/* IsConsistentWith: Is token consistent with expected word, case insensitive ? */
int IsConsistentWith (const char *token, const char *expected)
{
    int     i, len;

    if (strlen(token) > strlen(expected))
        return NO;

    len = (int) strlen (token);

    for (i=0; i<len; i++)
        {
        if (tolower(token[i]) != tolower(expected[i]))
            return NO;
        }

    return YES;
}





/* IsPartCompatible: Determine whether two partitions are nonoverlapping or nested (compatible) or
        incompatible (partially overlapping) */
int IsPartCompatible (SafeLong *smaller, SafeLong *larger, int length)
{
	int i;

    /* test first if they overlap */
	for (i=0; i<length; i++)
		if ((smaller[i]&larger[i]) != 0)
			break;

	/* if they overlap, they must be nested */
    if (i != length)	/* potentially incompatible */
		{
		for (i=0; i<length; i++)
			if ((smaller[i]|larger[i]) != larger[i])
				break;
		}
		
	if (i == length)	/* passed either one of the tests */
		return YES;
	else
		return NO;
}




/* IsPartNested: Test whether smaller partition is nested in larger partition */
int IsPartNested (SafeLong *smaller, SafeLong *larger, int length)

{

	int i;

	for (i=0; i<length; i++)
		if ((smaller[i] | larger[i]) != larger[i])
			break;
		
	if (i == length)
		return YES;
	else
		return NO;

}





/* IsSectionEmpty: Test whether section of two bitfields is empty */
int IsSectionEmpty (SafeLong *bitField1, SafeLong *bitField2, int length)
{
	int i;

	for (i=0; i<length; i++)
		if ((bitField1[i] & bitField2[i]) != 0)
			return NO;
		
	return YES;
}





/* LastBlock: Return file position of last block in file */
SafeLong LastBlock (FILE *fp, char *lineBuf, int longestLine)
{
	SafeLong	lastBlock;
	char	*word;
	
	lastBlock = 0L;
	rewind (fp);

	while ((fgets (lineBuf, longestLine, fp)) != NULL)
		{
		word = strtok (lineBuf, " ");
		if (strcmp (word, "begin") == 0)
			lastBlock = ftell (fp);
		}

	return lastBlock;
}





int LineTermType (FILE *fp)

{

	int			ch, nextCh, term;

	term = LINETERM_UNIX;	/* default if no line endings are found */
	while ((ch = getc(fp)) != EOF)
		{
		if ((ch == '\n') || (ch == '\r'))
			{
			if (ch == '\n')
				term = LINETERM_UNIX;
			else /* ch = '\r' */
				{
				/* First test below handles one-line MAC file */
				if (((nextCh = getc(fp)) == EOF) || (nextCh != '\n'))
					term = LINETERM_MAC;
				else
					term = LINETERM_DOS;
				}
			break;
			}
		}
	(void)fseek(fp, 0L, 0);		/* rewind */
	
	return (term);

}




/*The longest line in a file including line terminating characters present in binary mode.*/
int LongestLine (FILE *fp)

{

	int			ch, lineLength, longest;
	
	longest = 0;
	lineLength = 0;
	ch = fgetc(fp);
	while ( ch != EOF)
		{
		if((ch != '\n') && (ch != '\r'))
			{
			ch = fgetc(fp);
			lineLength++;
			continue;
			}
		if (ch == '\r')
			{
			if( (ch = fgetc(fp)) == '\n' )
				{
				/* windows \r\n */
				lineLength++;
				ch = fgetc(fp);
				}
			else
				{
				/* old mac \r */
				}
			}
		else  /*unix, linux,new mac or text mode read \n*/
			{
				ch = fgetc(fp);
			}

		if (lineLength > longest)
				longest = lineLength;
			lineLength = 0;
		/*
		if ((ch == '\n') || (ch == '\r'))
			{
			if (lineLength > longest)
				longest = lineLength;
			lineLength = 0;
			}
		else
			lineLength++;
			*/
		}
	rewind (fp);		/* rewind */
	
	return (longest+1); /*+1 to accommodate last character*/

}





/* LowerUpperMedian: Determine median and 95 % credible interval */
void LowerUpperMedian (MrBFlt *vals, int nVals, MrBFlt *lower, MrBFlt *upper, MrBFlt *median)

{    
    SortMrBFlt (vals, 0, nVals-1);
    
    *lower  = vals[(int)(0.025*nVals)];
    *upper  = vals[(int)(0.975*nVals)];
    *median = vals[nVals/2];

}




/* LowerUpperMedianHPD: Use a simple way to determine HPD */
void LowerUpperMedianHPD (MrBFlt *vals, int nVals, MrBFlt *lower, MrBFlt *upper, MrBFlt *median)
{
    int     i, width, theStart;
    MrBFlt  f, g, interval;

	SortMrBFlt (vals, 0, nVals-1);
    
    width = (int)(nVals * 0.95 + 0.5);
    theStart = 0;
    interval = vals[width-1] - vals[0];
    for (i=1; i<nVals-width; i++)
    {
        f = vals[i];
        g = vals[i+width];
        if (g - f < interval)
        {
            interval = g - f;
            theStart = i;
        }
    }

    *lower  = vals[theStart];
    *upper  = vals[theStart+width-1];
    *median = vals[nVals/2];
}





/*NOTE!!!! The result of this function should be used before consequtive call to it again. It means NVER use it like this:  printf( "%s %s", MbPrintNum (a),MbPrintNum (b) ) */
char *MbPrintNum (MrBFlt num)
{
    static char s[40];

    if (scientific == YES)
        sprintf(s,"%.*le", precision, num);
    else
        sprintf(s,"%.*lf", precision, num);

    return s;
}





void MeanVariance (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *var)

{

	int				i;
	MrBFlt			a, aOld, s, x;

	a = s = 0.0;
	for (i=0; i<nVals; i++)
		{
		x = vals[i];
		aOld = a;
		a += (x - a) / (MrBFlt) (i + 1);
		s += (x - a) * (x - aOld);
		}

	/* mean */
	(*mean) = a;
	
	/* variance */
	if (nVals <= 1)
		(*var) = 0.0;
	else
		(*var) = s / (nVals - 1);
			
}





void MrBayesPrint (char *format, ...)

{
	va_list ptr;

#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
		if (echoMB == YES)
			{
			va_start (ptr, format);
			vprintf (format, ptr);
			va_end(ptr);
			fflush (stdout);
			}
		if (logToFile == YES)
			{
			if (logFileFp == NULL)
				printf ("%s   Could not print log output to file\n", spacer);
			else
				{
				va_start (ptr, format);
				vfprintf (logFileFp, format, ptr);
				va_end(ptr);
				fflush (logFileFp);
				}
			}
		}
#	else
	if (chainParams.redirect == NO)
		{
		if (echoMB == YES)
			{
			va_start (ptr, format);
			vprintf (format, ptr);
			va_end(ptr);
			fflush (stdout);
			}
		if (logToFile == YES)
			{
			if (logFileFp == NULL)
				{
				printf ("%s   Could not print log output to file\n", spacer);
				logToFile = NO;
				}
			else
				{
				va_start (ptr, format);
				vfprintf (logFileFp, format, ptr);
				va_end(ptr);
				fflush (logFileFp);
				}
			}
		}
#	endif
}




void MrBayesPrintf (FILE *f, char *format, ...)

{
	va_list                 ptr;

#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
		va_start (ptr, format);
		vfprintf (f, format, ptr);
		va_end(ptr);
		fflush(f);
		}
#	else
	va_start (ptr, format);
	vfprintf (f, format, ptr);
	va_end(ptr);
	fflush(f);
#	endif
}





/** Next taxon in partition, for cycling over set bits in bit fields */
int NextTaxonInPartition(int currentTaxon, SafeLong *partition, int length)
{
    int         i, j, taxon;
    SafeLong    x;

    taxon = currentTaxon + 1;
    i = taxon / nBitsInALong;
    x = (1 << taxon % nBitsInALong);
    for (j=taxon%nBitsInALong; j<nBitsInALong; j++)
        {
        if (partition[i] & x)
            return taxon;
        taxon++;
        x <<= 1;
        }

    for (i++; i<length; i++)
		{
		x = 1;
		for (j=0; j<nBitsInALong; j++)
			{
			if (partition[i] & x)
				return taxon;
			taxon++;
			x <<= 1;
			}
		}    

    return taxon;
}





/* NumBits: Count bits in a bitfield */
int NumBits (SafeLong *x, int len)
{
	int         i, n=0;
    SafeLong    y;

	for (i=0; i<len; i++)
        {
        y = x[i];
        while (y != 0)
            {
		    y &= (y-1);
            n++;
            }
        }
	return n;
}





FILE *OpenBinaryFileR (char *name)

{

	FILE		*fp;
    char        fileName[100];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 99 - strlen(fileName));

    if ((fp = fopen (fileName, "rb")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileR (char *name)

{

	FILE		*fp;
    char        fileName[100];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 99 - strlen(fileName));

    if ((fp = fopen (fileName, "r")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileA (char *name)

{

	FILE		*fp;
    char        fileName[100];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 99 - strlen(fileName));

    if ((fp = fopen (fileName, "a+")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileW (char *name)

{

	FILE		*fp;
    char        fileName[100];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 99 - strlen(fileName));

	if ((fp = fopen (fileName, "w+")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





MrBFlt PotentialScaleReduction (MrBFlt **vals, int nRuns, int *count)

{

	int				i, j, nVals;
	MrBFlt			aW, aOldW, sW, sWj, aB, aOldB, sB, x, R2, weight;

	aB = sB = sW = sWj = 0.0;
    nVals = 0;
	for (j=0; j<nRuns; j++)
		{
		if(count[j]==0)
			{
			return -1.0;
			}
		nVals += count[j];
		aW = vals[j][0];
		for (i=1; i<count[j]; i++)
			{
			x = vals[j][i];
			aOldW = aW;
			aW += (x - aW) / (MrBFlt) (i + 1);
			sWj += (x - aW) * (x - aOldW);
			}
        sW += sWj / (MrBFlt)(count[j] - 1);
		x = aW;
		aOldB = aB;
		aB += (x - aB) / (MrBFlt) (j + 1);
		if (j!=0)
			sB += (x - aB) * (x - aOldB);
		}

	sB = sB / (MrBFlt) (nRuns - 1);
	sW = sW / (MrBFlt) (nRuns);

	weight = (MrBFlt) nVals / (MrBFlt) nRuns;
    if (sW > 0.0)
		{
		R2 = ((weight - 1.0) / weight) + ((MrBFlt)(nRuns + 1) / (MrBFlt) (nRuns)) * (sB / sW);
		return sqrt(R2);
		}
	else
		return -1.0;
}





/* SafeCalloc: Print error if out of memory */
void *SafeCalloc(size_t n, size_t s) {

    void *ptr = calloc(n, s);

    if(ptr==NULL)
        {
        MrBayesPrint ("%s   Out of memory\n", spacer);
        return NULL;
        }

    return ptr;
}





int SafeFclose(FILE **fp) {
	int retval=-1;
#if defined MPI_ENABLED
	if (proc_id == 0) {
#endif
	if( fp!=NULL && (*fp)!=NULL ) 
		retval=fclose(*fp);
	*fp = NULL;
#if defined MPI_ENABLED
	}
#endif
	return retval;	
}





/* SafeFree: Set pointer to freed space to NULL */
void SafeFree (void **ptr)
{
    free (*ptr);

    (*ptr) = NULL;
}





/* SafeMalloc: Print error if out of memory; clear memory */
void *SafeMalloc(size_t s) {

    void *ptr = malloc(s);

    if(ptr==NULL)
        {
        MrBayesPrint ("%s   Out of memory\n", spacer);
        return NULL;
        }

    return memset(ptr,0,s);
}





/* SafeRealloc: Print error if out of memory */
void *SafeRealloc(void *ptr, size_t s) {

    if (ptr == NULL)
        {
        ptr = malloc (s);
        memset(ptr, 0, s);
        }
    else
        ptr = realloc (ptr, s);

    if(ptr==NULL)
        {
        MrBayesPrint ("%s   Out of memory\n", spacer);
        return NULL;
        }

    return ptr;
}





/* SafeStrcat: Allocate or reallocate target to fit result; assumes ptr is NULL if not allocated */
char *SafeStrcat (char **target, const char *source)
{
    if (*target == NULL)
        *target = (char *) calloc (strlen(source)+1, sizeof(char));
    else
        *target = (char *) realloc ((void *)*target, (size_t)(strlen(source)+strlen(*target)+1)*sizeof(char));

    if (*target)
        strcat(*target,source);

    return (*target);
}





/* SetBit: Set a particular bit in a series of longs */
void SetBit (int i, SafeLong *bits)
{
	SafeLong		x;

	bits += i / nBitsInALong;

	x = 1 << (i % nBitsInALong);

	(*bits) |= x;
}





void SortInts (int *item, int *assoc, int count, int descendingOrder)

{

	SortInts2 (item, assoc, 0, count-1, descendingOrder);

}





void SortInts2 (int *item, int *assoc, int left, int right, int descendingOrder)

{

	register int	i, j, x, y;

	if (descendingOrder == YES)
		{
		i = left;
		j = right;
		x = item[(left+right)/2];
		do 
			{
			while (item[i] > x && i < right)
				i++;
			while (x > item[j] && j > left)
				j--;
			if (i <= j)
				{
				y = item[i];
				item[i] = item[j];
				item[j] = y;
				
				if (assoc)
					{
					y = assoc[i];
					assoc[i] = assoc[j];
					assoc[j] = y;
					}				
				i++;
				j--;
				}
			} while (i <= j);
		if (left < j)
			SortInts2 (item, assoc, left, j, descendingOrder);
		if (i < right)
			SortInts2 (item, assoc, i, right, descendingOrder);
		}
	else
		{
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
				
				if (assoc)
					{
					y = assoc[i];
					assoc[i] = assoc[j];
					assoc[j] = y;
					}				
				i++;
				j--;
				}
			} while (i <= j);
		if (left < j)
			SortInts2 (item, assoc, left, j, descendingOrder);
		if (i < right)
			SortInts2 (item, assoc, i, right, descendingOrder);
		}

}





/* SortMrBFlt: Sort in increasing order */
void SortMrBFlt (MrBFlt *item, int left, int right)

{

	register int	i, j;
	MrBFlt			x, temp;

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
			temp = item[i];
			item[i] = item[j];
			item[j] = temp;
				
			i++;
			j--;
			}
		} while (i <= j);
	if (left < j)
		SortMrBFlt (item, left, j);
	if (i < right)
		SortMrBFlt (item, i, right);
}





/* StrCmpCaseInsensitive: Case insensitive string comparison */
int StrCmpCaseInsensitive (char *s, char *t)
{
    int i, minLen;

    if (strlen(s) < strlen(t))
        minLen = (int) strlen(s);
    else
        minLen = (int) strlen(t);

    for (i=0; i<minLen; i++)
        if (tolower(s[i])!= tolower(t[i]))
            break;

    if (s[i] == '\0' && t[i] == '\0')
        return 0;
    else if (tolower(s[i]) > tolower(t[i]))
        return 1;
    else
        return -1;
}





/* StripComments: Strip possibly nested comments from the string s.
	Example: s="text1[text2[text3]]"-> s="text1" */
void StripComments (char *s)
{
	char	*t;
	int		inComment;

	inComment = 0;
	for (t=s; *s != '\0'; s++)
		{
		if (inComment == 0)
			{
			if (*s == '[')
				inComment++;
			else
				*t++ = *s;
			}
		else
			{
			if (*s == ']')
				inComment--;
            else if (*s == '[')
                inComment++;
			}
 		}
    *t = '\0';
}





FILE *TestOpenTextFileR (char *name)

{
    char        fileName[100];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 99 - strlen(fileName));

    return fopen (fileName, "r");	
}





int UpperTriangIndex(int i, int j, int size)
{
    if (i < j)
        return (2*size - i - 3) * i / 2 + j - 1;
    else
        return (2*size - j - 3) * j / 2 + i - 1;
}





int WantTo (const char *msg)
{
    char    s[100];
    int     i;

    MrBayesPrint ("%s   %s? (yes/no): ", spacer, msg);

	for (i=0; i<10; i++)
	    {
	    if (fgets (s, 98, stdin) == NULL)
		    {
		    MrBayesPrint ("%s   Failed to retrieve answer; will take that as a no\n", spacer);
		    return NO;
		    }

        /* Strip away the newline */
        s[strlen(s)-1] = '\0';

        /* Check answer */
        if (IsConsistentWith (s, "yes") == YES)
            return YES;
        else if (IsConsistentWith (s, "no") == YES)
            return NO;

        MrBayesPrint ("%s   Enter yes or no: ", spacer);
	    }

    MrBayesPrint ("%s   MrBayes does not understand; will take that as a no\n", spacer);

    return NO;
}



