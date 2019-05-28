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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "mb.h"
#include "globals.h"
#include "utils.h"

const char* const svnRevisionUtilsC="$Rev: 551 $";   /* Revision keyword which is expended/updated by svn on each commit/update*/

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





/* CopyResults: copy results from one file to another up to lastGen*/
int CopyResults (FILE *toFile, char *fromFileName, int lastGen)
{
    int     longestLine;
    char    *strBuf, *strCpy, *word;
    FILE    *fromFile;

    if ((fromFile = OpenBinaryFileR(fromFileName)) == NULL)
        return ERROR;

    longestLine = LongestLine(fromFile)+10;
    SafeFclose(&fromFile);
    strBuf = (char *) SafeCalloc (2*(longestLine+2),sizeof(char));
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


/* CopyProcessSsFile: copy results from one file to another up to lastStep. Also marginalLnLSS is collected for processed steps*/
int CopyProcessSsFile (FILE *toFile, char *fromFileName, int lastStep, MrBFlt *marginalLnLSS, MrBFlt * splitfreqSS)
{
    int     longestLine, run, curStep, i;
    double  tmp;
    char    *strBuf, *strCpy, *word, *tmpcp;
    FILE    *fromFile;

    if ((fromFile = OpenBinaryFileR(fromFileName)) == NULL)
        return ERROR;

    longestLine = LongestLine(fromFile)+10;
    SafeFclose(&fromFile);
    strBuf = (char *) SafeCalloc (2*(longestLine+2),sizeof(char));
    strCpy = strBuf + longestLine + 2;

    if ((fromFile = OpenTextFileR(fromFileName)) == NULL)
        return ERROR;
    
    while (fgets(strBuf,longestLine,fromFile)!=NULL)
        {
        strncpy(strCpy,strBuf,longestLine);
        word = strtok(strCpy," \t\n");
        /* atoi returns 0 when word is not integer number */
        if (atoi(word)>lastStep)
            break;
        fprintf (toFile,"%s",strBuf);
        fflush (toFile);
        curStep = atoi(word);
        if ( curStep > 0 )
            {
            strtok(NULL,"\t\n"); /*skip power*/
            for (run=0; run<chainParams.numRuns; run++)
                {
                tmpcp = strtok(NULL,"\t\n");
                if(tmpcp == NULL )
                    {
                    MrBayesPrint ("%s   Error: In .ss file not enough ellements on the string :%s        \n", spacer, strBuf);
                    return ERROR;
                    }
                tmp = atof(tmpcp);
                if(tmp == 0.0 )
                    {
                    MrBayesPrint ("%s   Error: Value of some step contribution is 0.0 or not a number in .ss file. Sting:%s        \n", spacer, strBuf);
                    return ERROR;
                    }
                marginalLnLSS[run]+=tmp;
                }
			for (i=0; i<numTopologies; i++)
				{
                tmpcp = strtok(NULL,"\t\n");
                if(tmpcp == NULL )
                    {
                    MrBayesPrint ("%s   Error: In .ss file not enough ellements on the string :%s        \n", spacer, strBuf);
                    return ERROR;
                    }
                tmp = atof(tmpcp);
                splitfreqSS[i*chainParams.numStepsSS + curStep-1] = tmp;                
    			}
            }
        }
    
    SafeFclose(&fromFile);
    free(strBuf);
    return (NO_ERROR);
}



/* CopyTreeResults: copy tree results upto lastGen from one file to another. numTrees is return containing number of trees that were copied. */
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
    strBuf = (char *) SafeCalloc (2*(longestLine+2),sizeof(char));
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
		firstTree = (int) ftell(fp);
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
int FromGrowthFxnToIndex(int *growthFxn)
{
    int     i, j, k, max, fxn[6];

    /* set local growth fxn to lexicographical max */
    for (i=0; i<6; i++)
        fxn[i] = i;

    /* decrease until we reach growthFxn */
    for (k=202; k>=0; k--)
        {
        for (i=0; i<6; i++)
            {
            if (fxn[i] != growthFxn[i])
                break;
            }
        if (i == 6)
            break;

        /* get next growth fxn */
        for (i=5; i>=0; i--)
            {
            fxn[i]--;
            if (fxn[i] >= 0)
                break;
            }

        if (i < 0)
            return -1;  /* error */
        else if (i < 5)
            {
            max = 0;
            for (j=0; j<=i; j++)
                {
                if (fxn[j] > max)
                    max = fxn[j];
                }
            fxn[++i] = max + 1;
            for (++i; i<6; i++)
                fxn[i] = fxn[i-1] + 1;
            }
        }

    return k;
}





/* Convert from model index to 0-based growth function over six states */
void FromIndexToGrowthFxn(int index, int *growthFxn)
{
    int     i, j, max, k;

    /* set growth fxn to lexicographical max */
    for (i=0; i<6; i++)
        growthFxn[i] = i;

    /* decrease until we reach index */
    for (k=202; k>index; k--)
        {
        for (i=5; i>=0; i--)
            {
            growthFxn[i]--;
            if (growthFxn[i] >= 0)
                break;
            }

        if (i < 0)
            return; /* ERROR */
        else if (i < 5)
            {
            max = 0;
            for (j=0; j<=i; j++)
                {
                if (growthFxn[j] > max)
                    max = growthFxn[j];
                }
            growthFxn[++i] = max + 1;
            for (++i; i<6; i++)
                growthFxn[i] = growthFxn[i-1] + 1;
            }
        }
}





/* GetIntSummary: Get summary statistics for a number of runs (int version) */
void GetIntSummary (int **vals, int nRows, int *rowCount, Stat *theStats, int HPD)

{

    int     i, j, nVals;
    MrBFlt  *theValues, *p;

    nVals = 0;
    for (i=0; i<nRows; i++)
        nVals += rowCount[i];

    theValues = (MrBFlt *) SafeCalloc (nVals, sizeof(MrBFlt));

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
    MrBFlt  *theValues, *p, *ESS;

    nVals = 0;
    for (i=0; i<nRows; i++)
        nVals += rowCount[i];

    theValues = (MrBFlt *) SafeMalloc ( (size_t) nVals* sizeof(MrBFlt));

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

    ESS =   (MrBFlt *) SafeMalloc ((size_t) nRows * sizeof(MrBFlt));

    EstimatedSampleSize (vals, nRows, rowCount, ESS);
    theStats->avrESS = theStats->minESS = ESS[0];
    for(i=1; i<nRows; i++)
        {
        theStats->avrESS += ESS[i];
        if(theStats->minESS > ESS[i] )
            {
            theStats->minESS = ESS[i];
            }
        }
    theStats->avrESS /=nRows;

    free (ESS);
    free (theValues);
}



/* HarmonicArithmeticMean: Calculate harmonic and arithmetic mean from log values */
int HarmonicArithmeticMeanOnLogs (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *harm_mean)
{
	int				i, reliable;
	MrBFlt			a, x, y, scaler, n;

	reliable = YES;
	
	scaler = vals[nVals-1];
	a  = n = 0.0;
	for (i=0; i<nVals; i++)
		{
		y = vals[i];
		y -= scaler;
		if (y > 400.0)
			{
            if (y > 5000.0)
                {
			    reliable = NO;
			    continue;
                }
            a /= exp( y - 100.0 ); 
            scaler += y - 100.0;
            y = 100.0;
			}
		
	    x = (MrBFlt) exp(y);
			
		if (n < 0.5)
			a = x;
		else
			{
            a += x;
			}
		n += 1.0;
		}

	/* arithmetic mean */
	(*mean) = (MrBFlt) log(a/n) + scaler;
	
	scaler = (MrBFlt) (0.0 - vals[nVals-1]);
	a  = n = 0.0;
	for (i=0; i<nVals; i++)
		{
		y = (MrBFlt) (0.0 - vals[i]);
		y -= scaler;
		if (y > 400.0)
			{
            if (y > 5000.0)
                {
			    reliable = NO;
			    continue;
                }
            a /= exp( y - 100.0 ); 
            scaler += y - 100.0;
            y = 100.0;
			}
		
	    x = (MrBFlt) exp(y);
			
		if (n < 0.5)
			a = x;
		else
			{
            a += x;
			}
		n += (MrBFlt) 1.0;
		}

	/* harmonic mean */
	(*harm_mean) = - (MrBFlt) log(a/n) - scaler;

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





/* IsSectionEmpty: Test whether union of bitField1 and bitField2 equal to bitField3*/
int IsUnionEqThird (SafeLong *bitField1, SafeLong *bitField2, SafeLong *bitField3, int length)
{
	int i;

	for (i=0; i<length; i++)
		if ((bitField1[i] | bitField2[i]) != bitField3[i] )
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
			lastBlock = (int) ftell (fp);
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





/*NOTE!!!! The result of this function should be used before consequtive call to it again. It means NEVER use it like this:  printf( "%s %s", MbPrintNum (a),MbPrintNum (b) ) */
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





/*  Compute mean and variance of log scaled values.

@param vals    pointer to values in log scale
@param nVals   number of "vals", minimum 1
@param mean    adress of variable where computed mean is returned by the function
@param var     adress of variable where computed variance is returned by the function. Could be set to NULL if this value need not to be returened. 
@param varEst  adress of variable where computed estimate of the population variance is returned, could be set to NULL if this value need not to be returened. 
               Could be set to NULL if this value need not to be returened.

Note: We devide by nVals or by (nVals-1) when var and varEst is calculated from the sum of square differences.
    */
void MeanVarianceLog (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *var, MrBFlt *varEst )

{

	int				i;
	MrBFlt			a, aOld, s, x, y, scaler;

	a = s = 0.0;
    scaler = vals[nVals-1];
	for (i=0; i<nVals; i++)
		{
		y = vals[i];
		y -= scaler;
		if (y > 200.0)
			{
            a /= exp( y - 100.0 );
            s /= exp( 2*(y - 100));
            scaler += y - 100.0;
            y = 100.0;
			}

        x=(MrBFlt)exp(y);

		aOld = a;
		a += (x - a) / (MrBFlt) (i + 1);
		s += (x - a) * (x - aOld);
		}

	/* mean */
	(*mean) = log(a) + scaler;
	
	/* variance */
    if( var!=NULL )
        {
	    if (nVals <= 1)
		    (*var) = 0.0;
	    else
		    (*var) = log( s / (nVals)) + 2*scaler;
        }

	/* variance */
    if( varEst!=NULL )
        {
	    if (nVals <= 1)
		    (*varEst) = 0.0;
	    else
		    (*varEst) = log( s / (nVals+1)) + 2*scaler;
        }
			
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
    char        fileName[200];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 199 - strlen(fileName));

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
    char        fileName[200];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 199 - strlen(fileName));

    if ((fp = fopen (fileName, "r")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, fileName);
		return (NULL);
		}
	else
		return (fp);
	
}




FILE *OpenTextFileRQuait (char *name)
{

	FILE		*fp;
    char        fileName[200];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 199 - strlen(fileName));

    if ((fp = fopen (fileName, "r")) == NULL)  
		{   
		return (NULL);
		}
	else
		return (fp);
	
}





FILE *OpenTextFileA (char *name)

{

	FILE		*fp;
    char        fileName[200];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 199 - strlen(fileName));

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
    char        fileName[200];

    strcpy(fileName, workingDir);
    strncat(fileName, name, 199 - strlen(fileName));

	if ((fp = fopen (fileName, "w+")) == NULL)  
		{   
		MrBayesPrint ("%s   Could not open file \"%s\"\n", spacer, name);
		return (NULL);
		}
	else
		return (fp);
	
}





/*!
\param vals[0..nRuns][count[]]   All records for all runs 
\param nRuns                     Number of runs
\param count[0..nRuns]           Number of records in each run
\return PSRF
*/
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





/*!
\param vals[0..nRuns][count[]]   All records for all runs 
\param nRuns                     Number of runs
\param count[0..nRuns]           Number of records in each run
\param returnESS[0..nRuns]       Is an arry in which the routine returns ESS values for each run.
*/
void EstimatedSampleSize (MrBFlt **vals, int nRuns, int *count, MrBFlt *returnESS)
{

	int		    i, j, lag, maxLag, samples;
    MrBFlt      *values, mean, del1, del2, varStat=0.0;
    MrBFlt      gammaStat[2000];
        
    for( i=0; i<nRuns; i++)
        {
        samples=count[i];
        values=vals[i];
        mean=0.0;
        for(j=0; j<samples; j++ )
            {
            mean+=values[j];
            }
        mean /=samples;

        maxLag = ((samples - 1) > 2000)?2000:(samples - 1);

        for (lag = 0; lag < maxLag; lag++)
            {
            gammaStat[lag]=0;
            for (j = 0; j < samples - lag; j++) 
                {
                del1 = values[j] - mean;
                del2 = values[j + lag] - mean;
                gammaStat[lag] += (del1 * del2);
                }

            gammaStat[lag] /= ((MrBFlt) (samples - lag));

            if (lag == 0) 
                {
                varStat = gammaStat[0];
                } 
            else if (lag % 2 == 0) 
                {
                if (gammaStat[lag - 1] + gammaStat[lag] > 0) 
                    {
                    varStat += 2.0 * (gammaStat[lag - 1] + gammaStat[lag]);
                    }
                else
                    maxLag = lag;
                }
            }
        returnESS[i] = (gammaStat[0] * samples) / varStat;
        }

}





/* SafeCalloc: Print error if out of memory */
void *SafeCalloc(size_t n, size_t s) {

    void *ptr;
    
    if( s*n == 0 )
        {
        //return NULL;
        }

    ptr= calloc(n, s);

    if(ptr==NULL)
        {
        MrBayesPrint ("%s   Out of memory. Most probable course for the problem is that MrBayes reached\n", spacer);
        MrBayesPrint ("%s   the limit of allowed memory for a process in your Operating System. Consult\n", spacer);
        MrBayesPrint ("%s   documentation of your OS how to extend the limit, or use 64 bit version OS \n", spacer);
        MrBayesPrint ("%s   and compile 64 bit version of MrBayes.                                     \n", spacer);
        MrBayesPrint ("%s   Segmentation fault may follow.                                             \n", spacer);
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

    void *ptr;

    if( s==0 )
        {
        return NULL;
        }

    ptr= malloc(s);

    if(ptr==NULL)
        {
        MrBayesPrint ("%s   Out of memory. Most probable course for the problem is that MrBayes reached\n", spacer);
        MrBayesPrint ("%s   the limit of allowed memory for a process in your Operating System. Consult\n", spacer);
        MrBayesPrint ("%s   documentation of your OS how to extend the limit, or use 64 bit version OS \n", spacer);
        MrBayesPrint ("%s   and compile 64 bit version of MrBayes.                                     \n", spacer);
        MrBayesPrint ("%s   Segmentation fault may follow.                                             \n", spacer);
        return NULL;
        }

    return memset(ptr,0,s);
}





/* SafeRealloc: Print error if out of memory */
void *SafeRealloc(void *ptr, size_t s) {

    if( s==0 )
        {
        free(ptr);
        return NULL;
        }

    if (ptr == NULL)
        {
        ptr = malloc (s);
        memset(ptr, 0, s);
        }
    else
        ptr = realloc (ptr, s);

    if(ptr==NULL)
        {
        MrBayesPrint ("%s   Out of memory. Most probable course for the problem is that MrBayes reached\n", spacer);
        MrBayesPrint ("%s   the limit of allowed memory for a process in your Operating System. Consult\n", spacer);
        MrBayesPrint ("%s   documentation of your OS how to extend the limit, or use 64 bit version OS \n", spacer);
        MrBayesPrint ("%s   and compile 64 bit version of MrBayes.                                     \n", spacer);
        MrBayesPrint ("%s   Segmentation fault may follow.                                             \n", spacer);
        return NULL;
        }

    return ptr;
}





/* SafeStrcat: Allocate or reallocate target to fit result; assumes ptr is NULL if not allocated */
char *SafeStrcat (char **target, const char *source)
{
    if (*target == NULL)
        *target = (char *) SafeCalloc (strlen(source)+1, sizeof(char));
    else
        *target = (char *) SafeRealloc ((void *)*target, (size_t)(strlen(source)+strlen(*target)+1)*sizeof(char));

    if (*target)
        strcat(*target,source);

    return (*target);
}




/* SafeStrcpy: Allocate or reallocate target to fit result; assumes ptr is NULL if not allocated */
char *SafeStrcpy (char **target, const char *source)
{

    *target = (char *) SafeRealloc ((void *)*target, (size_t)(strlen(source)+1)*sizeof(char));

    if (*target)
        strcpy(*target,source);

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





/*---------
|
|   UpdateGrowthFxn: We expect a set of unique indexes from 0 to 5
|      indicating a partition of 6 rates into sets. We make sure
|      the indices correspond to a restricted growth function here.
|
-----------------------*/
void UpdateGrowthFxn(int *growthFxn)
{
    int     i, j, max, fxn[6];

    for (i=0; i<6; i++)
        fxn[i] = -1;

    max = 0;
    for (i=0; i<6; i++)
        {
        if (fxn[i] != -1)
            continue;
        for (j=i; j<6; j++)
            {
            if (growthFxn[j] == growthFxn[i])
                fxn[j] = max;
            }
        max++;
        }

    for (i=0; i<6; i++)
        growthFxn[i] = fxn[i];   
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



