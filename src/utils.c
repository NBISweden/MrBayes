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
 *  Maxim Teslenko (maxkth@gmail.com)
 *  Chi Zhang (zhangchicool@gmail.com)
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
#include "best.h"
#include "command.h"
#include "mcmc.h"
#include "model.h"
#include "utils.h"

#define MAX_GAMMA_CATS                      20
#define PI                                  3.14159265358979324
#define PIOVER2                             1.57079632679489662
#define POINTGAMMA(prob,alpha,beta)         PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define PAI2                                6.283185307
#define TINY                                1.0e-20
#define EVALUATE_COMPLEX_NUMBERS            2
#if !defined(MAX)
#define MAX(a,b)                            (((a) > (b)) ? (a) : (b))
#endif
#if !defined(MIN)
#define MIN(a,b)                            (((a) < (b)) ? (a) : (b))
#endif
#define SQUARE(a)                           ((a)*(a))

/* local global variable */
char    noLabel[] = "";

/* local prototypes */
void    DatedNodeDepths (TreeNode *p, MrBFlt *nodeDepths, int *index);
void    DatedNodes (TreeNode *p, TreeNode **datedTips, int *index);
int     NConstrainedTips (TreeNode *p);
int     NDatedTips (TreeNode *p);
void    PrintNode (char **s, int *len, TreeNode *p, int isRooted);
void    ResetPolyNode (PolyNode *p);
void    ResetTreeNode (TreeNode *p);
void    SetNodeDepths (Tree *t);

void    AddTwoMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result);
void    BackSubstitutionRow (int dim, MrBFlt **u, MrBFlt *b);
void    Balanc (int dim, MrBFlt **a, int *low, int *high, MrBFlt *scale);
void    BalBak (int dim, int low, int high, MrBFlt *scale, int m, MrBFlt **z);
MrBFlt  BetaCf (MrBFlt a, MrBFlt b, MrBFlt x);
MrBFlt  BetaQuantile (MrBFlt alpha, MrBFlt beta, MrBFlt x);
MrBFlt  CdfBinormal (MrBFlt h1, MrBFlt h2, MrBFlt r);
MrBFlt  CdfNormal (MrBFlt x);
complex Complex (MrBFlt a, MrBFlt b);
MrBFlt  ComplexAbsoluteValue (complex a);
complex ComplexAddition (complex a, complex b);
complex ComplexConjugate (complex a);
complex ComplexDivision (complex a, complex b);
void    ComplexDivision2 (MrBFlt ar, MrBFlt ai, MrBFlt br, MrBFlt bi, MrBFlt *cr, MrBFlt *ci);
complex ComplexExponentiation (complex a);
int     ComplexInvertMatrix (int dim, complex **a, MrBFlt *dwork, int *indx, complex **aInverse, complex *col);
complex ComplexLog (complex a);
void    ComplexLUBackSubstitution (int dim, complex **a, int *indx, complex *b);
int     ComplexLUDecompose (int dim, complex **a, MrBFlt *vv, int *indx, MrBFlt *pd);
complex ComplexMultiplication (complex a, complex b);
complex ComplexSquareRoot (complex a);
complex ComplexSubtraction (complex a, complex b);
int     ComputeEigenSystem (int dim, MrBFlt **a, MrBFlt *v, MrBFlt *vi, MrBFlt **u, int *iwork, MrBFlt *dwork);
void    ComputeLandU (int dim, MrBFlt **aMat, MrBFlt **lMat, MrBFlt **uMat);
void    ComputeMatrixExponential (int dim, MrBFlt **a, int qValue, MrBFlt **f);
void    DivideByTwos (int dim, MrBFlt **a, int power);
MrBFlt  D_sign (MrBFlt a, MrBFlt b);
int     EigensForRealMatrix (int dim, MrBFlt **a, MrBFlt *wr, MrBFlt *wi, MrBFlt **z, int *iv1, MrBFlt *fv1);
void    ElmHes (int dim, int low, int high, MrBFlt **a, int *interchanged);
void    ElTran (int dim, int low, int high, MrBFlt **a, int *interchanged, MrBFlt **z);
void    Exchange (int j, int k, int l, int m, int n, MrBFlt **a, MrBFlt *scale);
MrBFlt  Factorial (int x);
void    ForwardSubstitutionRow (int dim, MrBFlt **L, MrBFlt *b);
MrBFlt  GammaRandomVariable (MrBFlt a, MrBFlt b, RandLong *seed);
void    GaussianElimination (int dim, MrBFlt **a, MrBFlt **bMat, MrBFlt **xMat);
int     Hqr2 (int dim, int low, int high, MrBFlt **h, MrBFlt *wr, MrBFlt *wi, MrBFlt **z);
MrBFlt  IncompleteBetaFunction (MrBFlt alpha, MrBFlt beta, MrBFlt x);
MrBFlt  IncompleteGamma (MrBFlt x, MrBFlt alpha, MrBFlt LnGamma_alpha);
int     InvertMatrix (int dim, MrBFlt **a, MrBFlt *col, int *indx, MrBFlt **aInv);
MrBFlt  LBinormal (MrBFlt h1, MrBFlt h2, MrBFlt r);
int     LogBase2Plus1 (MrBFlt x);
void    LUBackSubstitution (int dim, MrBFlt **a, int *indx, MrBFlt *b);
int     LUDecompose (int dim, MrBFlt **a, MrBFlt *vv, int *indx, MrBFlt *pd);
void    MultiplyMatrixByScalar (int dim, MrBFlt **a, MrBFlt scalar, MrBFlt **result);
MrBFlt  PointChi2 (MrBFlt prob, MrBFlt v);
void    PrintComplexVector (int dim, complex *vec);
void    PrintSquareComplexMatrix (int dim, complex **m);
void    PrintSquareDoubleMatrix (int dim, MrBFlt **matrix);
void    PrintSquareIntegerMatrix (int dim, int **matrix);
complex ProductOfRealAndComplex (MrBFlt a, complex b);
MrBFlt  RndGamma (MrBFlt s, RandLong *seed);
MrBFlt  RndGamma1 (MrBFlt s, RandLong *seed);
MrBFlt  RndGamma2 (MrBFlt s, RandLong *seed);
int     SetQvalue (MrBFlt tol);
void    SetToIdentity (int dim, MrBFlt **matrix);
MrBFlt  Tha (MrBFlt h1, MrBFlt h2, MrBFlt a1, MrBFlt a2);
void    TiProbsUsingEigens (int dim, MrBFlt *cijk, MrBFlt *eigenVals, MrBFlt v, MrBFlt r, MrBFlt **tMat, MrBFlt **fMat, MrBFlt **sMat);
void    TiProbsUsingPadeApprox (int dim, MrBFlt **qMat, MrBFlt v, MrBFlt r, MrBFlt **tMat, MrBFlt **fMat, MrBFlt **sMat);

MrBFlt  QuantileLogNormal (MrBFlt prob, MrBFlt mu, MrBFlt sigma);
int     DiscreteLogNormal (MrBFlt *rK, MrBFlt sigma, int K, int median);
MrBFlt  LogNormalPoint (MrBFlt x, MrBFlt mu, MrBFlt sigma);

/* AddBitfield: Add bitfield to list of bitfields. The function uses global variable nLongsNeeded. */
int AddBitfield (BitsLong ***list, int listLen, int *set, int setLen)
{
    int     i, nLongsNeeded;

    nLongsNeeded = (setLen - 1) / nBitsInALong + 1;

    (*list) = (BitsLong **) SafeRealloc ((void *)(*list), ((size_t)listLen+1)*sizeof(BitsLong *));
    if (!(*list))
        return ERROR;
    
    (*list)[listLen] = (BitsLong *) SafeMalloc ((size_t)nLongsNeeded*sizeof(BitsLong));
    if (!(*list)[listLen])
        return ERROR;

    ClearBits ((*list)[listLen], nLongsNeeded);
    for (i=0; i<setLen; i++)
        if (set[i] == YES)
            SetBit(i, (*list)[listLen]);

    return NO_ERROR;
}


#if defined (SSE_ENABLED)   /* SSE or more advanced SIMD */
void * AlignedMalloc (size_t size, size_t alignment)
{
    void *mem;

    #if defined GCC_SIMD    /* gcc compiler */
    if (posix_memalign (&mem, alignment, size))
        return 0;
    #elif defined ICC_SIMD   /* icc compiler */
    mem = _mm_malloc (size, alignment);
    #elif defined MS_VCPP_SIMD  /* ms visual */
    mem = _aligned_malloc (size, alignment);
    #else
    mem = malloc (size);
    #endif

    return mem;
}


void AlignedSafeFree (void **ptr)
{

    #if defined ICC_VEC     /* icc compiler */
    _mm_free (*ptr);
    #elif defined MS_VCPP_VEC  /* ms visual */
    _aligned_free (*ptr);
    #else
    free (*ptr);
    #endif
    
    (*ptr) = NULL;
}
#endif


int AreBitfieldsEqual (BitsLong *p, BitsLong *q, int length)
{
    int i;
    
    for (i=0; i<length; i++)
        {
        if (p[i] != q[i])
            return NO;
        }
    
    return YES;
}


/*----------------------------------------------------------------
|
|   Bit: return 1 if bit n is set in BitsLong *p
|       else return 0
|
-----------------------------------------------------------------*/
int Bit (int n, BitsLong *p)
{
    BitsLong        x, bitsLongOne;

    bitsLongOne = 1;

    p += n / nBitsInALong;
    x = bitsLongOne << (n % nBitsInALong);

    if ((x & (*p)) == 0)
        return 0;
    else
        return 1;

}


/* ClearBit: Clear one bit in a bitfield */
void ClearBit (int i, BitsLong *bits)
{
    BitsLong        x, bitsLongOne=1;

    bits += i / nBitsInALong;

    x = bitsLongOne << (i % nBitsInALong);
    x ^= bitsLongWithAllBitsSet;

    (*bits) &= x;
}


/* ClearBits: Clear all bits in a bitfield */
void ClearBits (BitsLong *bits, int nLongs)
{
    int     i;
    
    for (i=0; i<nLongs; i++)
        bits[i] = 0;
}


/* Copy bitfields */
void CopyBits (BitsLong *dest, BitsLong *source, int length)
{
    int     i;

    for (i=0; i<length; i++)
        dest[i] = source[i];
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
        strncpy (strCpy,strBuf,longestLine);
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
        strncpy (strCpy,strBuf,longestLine);
        word = strtok(strCpy," \t\n");
        /* atoi returns 0 when word is not integer number */
        if (atoi(word)>lastStep)
            break;
        fprintf (toFile,"%s",strBuf);
        fflush (toFile);
        curStep = atoi(word);
        if (curStep > 0)
            {
            strtok(NULL,"\t\n"); /*skip power*/
            for (run=0; run<chainParams.numRuns; run++)
                {
                tmpcp = strtok(NULL,"\t\n");
                if (tmpcp == NULL)
                    {
                    MrBayesPrint ("%s   Error: In .ss file not enough ellements on the string :%s        \n", spacer, strBuf);
                    return ERROR;
                    }
                tmp = atof(tmpcp);
                if (tmp == 0.0)
                    {
                    MrBayesPrint ("%s   Error: Value of some step contribution is 0.0 or not a number in .ss file. Sting:%s        \n", spacer, strBuf);
                    return ERROR;
                    }
                marginalLnLSS[run]+=tmp;
                }
            for (i=0; i<numTopologies; i++)
                {
                tmpcp = strtok(NULL,"\t\n");
                if (tmpcp == NULL)
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
        strncpy (strCpy,strBuf,longestLine);
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
int FirstTaxonInPartition (BitsLong *partition, int length)
{
    int         i, j, nBits, taxon;
    BitsLong    x, bitsLongOne=1;

    nBits = sizeof(BitsLong) * 8;

    taxon = 0;
    for (i=0; i<length; i++)
        {
        x = bitsLongOne;
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
long FirstTree (FILE *fp, char *lineBuf, int longestLine)
{
    long    firstTree;
    char    *word;
    
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


void FlipBits (BitsLong *partition, int length, BitsLong *mask)
{
    int         i;
    
    for (i=0; i<length; i++)
        {
        partition[i] ^= mask[i];
        }
}


/*-----------------------------------------------------------------
|
|   FlipOneBit: flip bit n in BitsLong *p
|
------------------------------------------------------------------*/
void FlipOneBit (int n, BitsLong *p)
{
    BitsLong        x, bitsLongOne=1;

    p += n/nBitsInALong;
    x = bitsLongOne << (n % nBitsInALong);
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

    theValues = (MrBFlt *) SafeMalloc ((size_t)nVals * sizeof(MrBFlt));

    /* extract values */
    p = theValues;
    for (i=0; i<nRows; i++)
        {
        memcpy ((void *)(p), (void *)(vals[i]), (size_t)rowCount[i] * sizeof(MrBFlt));
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

    ESS = (MrBFlt *) SafeMalloc ((size_t)nRows * sizeof(MrBFlt));

    EstimatedSampleSize (vals, nRows, rowCount, ESS);
    theStats->avrESS = theStats->minESS = ESS[0];
    for (i=1; i<nRows; i++)
        {
        theStats->avrESS += ESS[i];
        if (theStats->minESS > ESS[i])
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
    int             i, reliable;
    MrBFlt          a, x, y, scaler, n;

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
            a /= exp(y - 100.0);
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
            a /= exp(y - 100.0);
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


/* IsBitSet: Is bit i set in BitsLong *bits ? */
int IsBitSet (int i, BitsLong *bits)
{
    BitsLong        x, bitsLongOne=1;

    bits += i / nBitsInALong;

    x = bitsLongOne << (i % nBitsInALong);

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
int IsPartCompatible (BitsLong *smaller, BitsLong *larger, int length)
{
    int i;

    /* test first if they overlap */
    for (i=0; i<length; i++)
        if ((smaller[i]&larger[i]) != 0)
            break;

    /* if they overlap, they must be nested */
    if (i != length)    /* potentially incompatible */
        {
        for (i=0; i<length; i++)
            if ((smaller[i]|larger[i]) != larger[i])
                break;
        }
        
    if (i == length)    /* passed either one of the tests */
        return YES;
    else
        return NO;
}


/* IsPartNested: Test whether smaller partition is nested in larger partition */
int IsPartNested (BitsLong *smaller, BitsLong *larger, int length)
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
int IsSectionEmpty (BitsLong *bitField1, BitsLong *bitField2, int length)
{
    int i;

    for (i=0; i<length; i++)
        if ((bitField1[i] & bitField2[i]) != 0)
            return NO;
        
    return YES;
}


/* IsSectionEmpty: Test whether union of bitField1 and bitField2 equal to bitField3*/
int IsUnionEqThird (BitsLong *bitField1, BitsLong *bitField2, BitsLong *bitField3, int length)
{
    int i;

    for (i=0; i<length; i++)
        if ((bitField1[i] | bitField2[i]) != bitField3[i])
            return NO;
        
    return YES;
}


/* LastBlock: Return file position of last block in file */
long LastBlock (FILE *fp, char *lineBuf, int longestLine)
{
    long    lastBlock;
    char    *word;
    
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
    int         ch, nextCh, term;

    term = LINETERM_UNIX;   /* default if no line endings are found */
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
    (void)fseek(fp, 0L, 0);     /* rewind */
    
    return (term);
}


/*The longest line in a file including line terminating characters present in binary mode.*/
int LongestLine (FILE *fp)
{
    int         ch, lineLength, longest;
    
    longest = 0;
    lineLength = 0;
    ch = fgetc(fp);
    while (ch != EOF)
        {
        if ((ch != '\n') && (ch != '\r'))
            {
            ch = fgetc(fp);
            lineLength++;
            continue;
            }
        if (ch == '\r')
            {
            if ((ch = fgetc(fp)) == '\n')
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
    rewind (fp);        /* rewind */
    
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


MrBFlt MaximumValue (MrBFlt x, MrBFlt y)
{
    if (x > y)
        return (x);
    else
        return (y);
}


MrBFlt MinimumValue (MrBFlt x, MrBFlt y)
{
    if (x < y)
        return (x);
    else
        return (y);
}


/* NOTE!!!! The result of this function should be used before consequtive call to it again.
   It means NEVER use it like this:  printf ("%s %s", MbPrintNum (a),MbPrintNum (b)) */
char *MbPrintNum (MrBFlt num)
{
    static char s[40];

    if (scientific == YES)
        sprintf (s,"%.*le", precision, num);
    else
        sprintf (s,"%.*lf", precision, num);

    return s;
}


void MeanVariance (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *var)
{
    int             i;
    MrBFlt          a, aOld, s, x;

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
Note: We devide by nVals or by (nVals-1) when var and varEst is calculated from the sum of square differences. */
void MeanVarianceLog (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *var, MrBFlt *varEst)
{
    int             i;
    MrBFlt          a, aOld, s, x, y, scaler;

    a = s = 0.0;
    scaler = vals[nVals-1];
    for (i=0; i<nVals; i++)
        {
        y = vals[i];
        y -= scaler;
        if (y > 200.0)
            {
            a /= exp(y - 100.0);
            s /= exp(2*(y - 100));
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
    if (var!=NULL)
        {
        if (nVals <= 1)
            (*var) = 0.0;
        else
            (*var) = log(s / nVals) + 2*scaler;
        }

    /* variance */
    if (varEst!=NULL)
        {
        if (nVals <= 1)
            (*varEst) = 0.0;
        else
            (*varEst) = log(s / (nVals+1)) + 2*scaler;
        }
}


void MrBayesPrint (char *format, ...)
{
    va_list ptr;

#   if defined (MPI_ENABLED)
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
#   else
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
#   endif
}


void MrBayesPrintf (FILE *f, char *format, ...)
{
    va_list                 ptr;

#   if defined (MPI_ENABLED)
    if (proc_id == 0)
        {
        va_start (ptr, format);
        vfprintf (f, format, ptr);
        va_end(ptr);
        fflush(f);
        }
#   else
    va_start (ptr, format);
    vfprintf (f, format, ptr);
    va_end(ptr);
    fflush(f);
#   endif
}


/** Next taxon in partition, for cycling over set bits in bit fields */
int NextTaxonInPartition(int currentTaxon, BitsLong *partition, int length)
{
    int         i, j, taxon;
    BitsLong    x, bitsLongOne=1;

    taxon = currentTaxon + 1;
    i = taxon / nBitsInALong;
    x = (bitsLongOne << taxon % nBitsInALong);
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


/* NBits: count bits in an int */
int NBits (int x)
{
    int n=0;
    
    for (n=0; x != 0; n++)
        x &= (x-1);
    
    return n;
}


/* NumBits: Count bits in a bitfield */
int NumBits (BitsLong *x, int len)
{
    int         i, n=0;
    BitsLong    y;

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
    FILE        *fp;
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
    FILE        *fp;
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
    FILE        *fp;
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
    FILE        *fp;
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
    FILE        *fp;
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
    int             i, j, nVals;
    MrBFlt          aW, aOldW, sW, sWj, aB, aOldB, sB, x, R2, weight;

    aB = sB = sW = sWj = 0.0;
    nVals = 0;
    for (j=0; j<nRuns; j++)
        {
        if (count[j]==0)
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
    int         i, j, lag, maxLag, samples;
    MrBFlt      *values, mean, del1, del2, varStat=0.0;
    MrBFlt      gammaStat[2000];
        
    for (i=0; i<nRuns; i++)
        {
        samples=count[i];
        values=vals[i];
        mean=0.0;
        for (j=0; j<samples; j++)
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

    ptr = calloc(n, s);

    if (ptr==NULL && n*s > 0)
        {
        MrBayesPrint ("%s   Out of memory. Most probable course for the problem is that MrBayes reached\n", spacer);
        MrBayesPrint ("%s   the limit of allowed memory for a process in your Operating System. Consult\n", spacer);
        MrBayesPrint ("%s   documentation of your OS how to extend the limit, or use 64 bit version OS \n", spacer);
        MrBayesPrint ("%s   and compile 64 bit version of MrBayes.                                     \n", spacer);
        MrBayesPrint ("%s   Segmentation fault may follow.                                             \n", spacer);
        }

    return ptr;
}


int SafeFclose(FILE **fp) {
    int retval=-1;
#   if defined MPI_ENABLED
    if (proc_id == 0) {
#   endif
    if (fp!=NULL && (*fp)!=NULL) 
        retval=fclose(*fp);
    *fp = NULL;
#   if defined MPI_ENABLED
    }
#   endif
    return retval;  
}


/* SafeFree: Set pointer to freed space to NULL */
void SafeFree (void **ptr)
{
    free (*ptr);

    (*ptr) = NULL;
}


/* SafeMalloc: Print error if out of memory; clear memory */
void *SafeMalloc (size_t s)
{
    void *ptr;

    ptr = calloc(1, s);

    if (ptr==NULL && s > 0)
        {
        MrBayesPrint ("%s   Out of memory. Most probable course for the problem is that MrBayes reached\n", spacer);
        MrBayesPrint ("%s   the limit of allowed memory for a process in your Operating System. Consult\n", spacer);
        MrBayesPrint ("%s   documentation of your OS how to extend the limit, or use 64 bit version OS \n", spacer);
        MrBayesPrint ("%s   and compile 64 bit version of MrBayes.                                     \n", spacer);
        MrBayesPrint ("%s   Segmentation fault may follow.                                             \n", spacer);
        }

    return ptr;
}


/* SafeRealloc: Print error if out of memory */
void *SafeRealloc (void *ptr, size_t s)
{
    void *tmp;

    tmp = realloc(ptr, s);

    if (tmp != NULL) {
        if (ptr == NULL)
            memset(tmp, 0, s);
        ptr = tmp;
    }

    if (tmp == NULL && s > 0) {
        MrBayesPrint ("%s   Out of memory. Most probable course for the problem is that MrBayes reached\n", spacer);
        MrBayesPrint ("%s   the limit of allowed memory for a process in your Operating System. Consult\n", spacer);
        MrBayesPrint ("%s   documentation of your OS how to extend the limit, or use 64 bit version OS \n", spacer);
        MrBayesPrint ("%s   and compile 64 bit version of MrBayes.                                     \n", spacer);
        MrBayesPrint ("%s   Segmentation fault may follow.                                             \n", spacer);
        }

    return ptr;
}


/* SafeStrcat: Allocate or reallocate target to fit result; assumes ptr is NULL if not allocated */
char *SafeStrcat (char **target, const char *source)
{
    if (*target == NULL)
        *target = (char *) SafeCalloc (strlen(source)+1, sizeof(char));
    else
        *target = (char *) SafeRealloc ((void *)*target, (strlen(source)+strlen(*target)+1)*sizeof(char));

    if (*target)
        strcat(*target, source);

    return (*target);
}


/* SafeStrcpy: Allocate or reallocate target to fit result; assumes ptr is NULL if not allocated */
char *SafeStrcpy (char **target, const char *source)
{
    *target = (char *) SafeRealloc ((void *)*target, (strlen(source)+1)*sizeof(char));

    if (*target)
        strcpy(*target,source);

    return (*target);
}


/* SetBit: Set a particular bit in a series of longs */
void SetBit (int i, BitsLong *bits)
{
    BitsLong        x, bitsLongOne=1;

    bits += i / nBitsInALong;

    x = bitsLongOne << (i % nBitsInALong);

    (*bits) |= x;
}


/* This routine is not called from anywhere */
#if 0
void SortInts (int *item, int *assoc, int count, int descendingOrder)
{
    SortInts2 (item, assoc, 0, count-1, descendingOrder);
}


void SortInts2 (int *item, int *assoc, int left, int right, int descendingOrder)
{
    register int    i, j, x, y;

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
#endif


int MrBFlt_cmp(const void *a, const void *b)
{
    MrBFlt x = *(MrBFlt *)a;
    MrBFlt y = *(MrBFlt *)b;

    if (x < y) {
        return -1;
    }
    else if (x > y) {
        return 1;
    }
    return 0;
}

/* SortMrBFlt: Sort in increasing order */
void SortMrBFlt (MrBFlt *item, int left, int right)
{
    qsort((void *)item, right - left + 1, sizeof(MrBFlt), &MrBFlt_cmp);

#if 0
    register int    i, j;
    MrBFlt          x, temp;

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
#endif
}


/* StrCmpCaseInsensitive: Case insensitive string comparison */
int StrCmpCaseInsensitive (char *s, char *t)
{
    return strcasecmp(s, t);
#if 0
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
#endif
}


/* StripComments: Strip possibly nested comments from the string s.
    Example: s="text1[text2[text3]]"-> s="text1" */
void StripComments (char *s)
{
    char    *t;
    int     inComment;

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


/* the following are moved from tree.c */
/* AddToTreeList: Add tree at end of tree list */
int AddToTreeList (TreeList *treeList, Tree *tree)
{
    TreeListElement     *listElement = (TreeListElement *) SafeCalloc (1, sizeof(TreeListElement));
    if (!listElement)
        return (ERROR);

    listElement->order = (int *) SafeCalloc (tree->nIntNodes-1, sizeof(int));
    if (!listElement->order)
        return (ERROR);
    listElement->next = NULL;

    if (treeList->last == NULL)
        treeList->last = treeList->first = listElement;
    else
        {
        treeList->last->next = listElement;
        treeList->last = listElement;
        }

    if (tree->isRooted)
        StoreRTopology (tree, listElement->order);
    else
        StoreUTopology (tree, listElement->order);

    return (NO_ERROR);
}


/* AllocatePolyTree: Allocate memory space for a polytomous tree */
PolyTree *AllocatePolyTree (int numTaxa)
{
    int         i;
    PolyTree    *pt;

    pt = (PolyTree *) SafeCalloc (1, sizeof (PolyTree));
    if (!pt)
        return (NULL);

    pt->memNodes = 2*numTaxa;  
    pt->nodes = (PolyNode *) SafeCalloc (2*numTaxa, sizeof(PolyNode));
    pt->allDownPass = (PolyNode **) SafeCalloc (3*numTaxa, sizeof (PolyNode *));
    pt->intDownPass = pt->allDownPass + 2*numTaxa;
    if (pt->nodes == NULL || pt->allDownPass == NULL)
        {
        free (pt->nodes);
        free (pt->allDownPass);
        free (pt);
        return (NULL);
        }

    /* initialize nodes and set index and memoryIndex */
    for (i=0; i<2*numTaxa; i++)
        {
        ResetPolyNode(&pt->nodes[i]);
        pt->nodes[i].memoryIndex = i;
        pt->nodes[i].index = i;
        }

    /* initialize tree properties */
    pt->nNodes = pt->nIntNodes = 0;
    pt->root = NULL;
    pt->brlensDef = NO;
    pt->isRooted = NO;
    pt->isClock = NO;
    pt->isCalibrated = NO;
    pt->isRelaxed = NO;
    pt->clockRate = 0.0;
    strcpy(pt->name,"");

    /* initialize bitsets */
    pt->bitsets = NULL;
    
    /* initialize relaxed clock parameters */
    pt->nESets = 0;
    pt->nEvents = NULL;
    pt->position = NULL;
    pt->rateMult = NULL;
    pt->eSetName = NULL;

    pt->nBSets = 0;
    pt->effectiveBrLen = NULL;
    pt->bSetName = NULL;

    /* initialize population size set parameters */
    pt->popSizeSet = NO;
    pt->popSize = NULL;
    pt->popSizeSetName = NULL;

    return (pt);
}


/* AllocatePolyTreeRelClockParams: Allocate space for relaxed clock parameters */
int AllocatePolyTreeRelClockParams (PolyTree *pt, int nBSets, int nESets)
{
    int     i;

    /* free previous clock params if any */
    FreePolyTreeRelClockParams (pt);

    /* set number of sets */
    pt->nBSets = nBSets;
    pt->nESets = nESets;

    /* we do not allocate space for the actual names here; these will be NULL pointers */

    /* take care of branch params */
    if (pt->nBSets > 0)
        {
        pt->bSetName = (char **) SafeCalloc (pt->nBSets, sizeof (char *));
        pt->effectiveBrLen = (MrBFlt **) SafeCalloc (pt->nBSets, sizeof (MrBFlt *));
        for (i=0; i<pt->nBSets; i++)
            pt->effectiveBrLen[i] = (MrBFlt *) SafeCalloc (pt->memNodes, sizeof(MrBFlt));
        }
    
    /* take care of breakpoint params */
    if (pt->nESets > 0)
        {
        pt->eSetName = (char **) SafeCalloc (pt->nESets, sizeof(char *));
        pt->nEvents = (int **) SafeCalloc (pt->nESets, sizeof(int *));
        pt->position = (MrBFlt ***) SafeCalloc (pt->nESets, sizeof(MrBFlt **));
        pt->rateMult = (MrBFlt ***) SafeCalloc (pt->nESets, sizeof(MrBFlt **));
        for (i=0; i<pt->nESets; i++)
            {
            pt->nEvents[i] = (int *) SafeCalloc (pt->memNodes, sizeof(int));
            pt->position[i] = (MrBFlt **) SafeCalloc (pt->memNodes, sizeof(MrBFlt *));
            pt->rateMult[i] = (MrBFlt **) SafeCalloc (pt->memNodes, sizeof(MrBFlt *));
            }
        }

    return (NO_ERROR);
}


/* AllocatePolyTreePartitions: Allocate space for and set partitions for polytomous tree */
int AllocatePolyTreePartitions (PolyTree *pt)
{
    int         i, nLongsNeeded, numTaxa;

    /* get some handy numbers */
    numTaxa = pt->memNodes/2;
    nLongsNeeded = (numTaxa -1) / nBitsInALong + 1;

    /* allocate space */
    pt->bitsets = (BitsLong *) SafeRealloc ((void *)pt->bitsets, pt->memNodes*nLongsNeeded*sizeof(BitsLong));
    if (pt->bitsets == NULL)
        return (ERROR);
    for (i=0; i<pt->memNodes*nLongsNeeded; i++)
        pt->bitsets[i] = 0;
    
    /* set node partition pointers */
    for (i=0; i<pt->memNodes; i++)
        pt->nodes[i].partition = pt->bitsets + i*nLongsNeeded;

    /* clear and set partitions; if the tree is empty, nothing is set */
    ResetPolyTreePartitions(pt);
    
    return (NO_ERROR);
}


/* AllocateTree: Allocate memory space for a tree (unrooted or rooted) */
Tree *AllocateTree (int numTaxa)
{
    int     i;
    Tree    *t;
    
    t = (Tree *) SafeCalloc (1, sizeof (Tree));
    if (t == NULL)
        return NULL;

    /* initialize basic tree properties */
    t->memNodes = 2*numTaxa;
    strcpy (t->name, "");
    
    t->isRooted = NO;
    t->isClock = NO;

    t->checkConstraints = NO;
    t->nConstraints = 0;
    t->nLocks = 0;
    t->isCalibrated = NO;
    t->nNodes = t->nIntNodes = 0;
    t->nRelParts = 0;
    t->relParts = NULL;

    /* initialize pointers */
    t->bitsets = NULL;
    t->flags = NULL;
    t->constraints = NULL;

    /* allocate and initialize nodes and node arrays (enough for both rooted and unrooted trees) */
    t->nNodes = 0;
    t->nIntNodes = 0;
    if ((t->nodes = (TreeNode *) SafeCalloc (2*numTaxa, sizeof (TreeNode))) == NULL)
        {
        free (t);
        return NULL;
        }
    if ((t->allDownPass = (TreeNode **) SafeCalloc (3*numTaxa, sizeof (TreeNode *))) == NULL)
        {
        free (t->nodes);
        free (t);
        return NULL;
        }
    t->intDownPass = t->allDownPass + t->memNodes;
    
    /* initialize nodes and set index and memoryIndex */
    for (i=0; i<t->memNodes; i++)
        {
        ResetTreeNode(&t->nodes[i]);
        t->nodes[i].memoryIndex = i;
        t->nodes[i].index = i;
        }

    return t;
}


/* AllocateFixedTree: Allocate memory space for a fixed unrooted or rooted tree */
Tree *AllocateFixedTree (int numTaxa, int isRooted)
{
    int     i;
    Tree    *t;
    
    t = (Tree *) SafeCalloc (1, sizeof (Tree));
    if (t == NULL)
        return NULL;

    /* initialize basic tree properties */
    if (isRooted == YES)
        t->memNodes = 2*numTaxa;
    else
        t->memNodes = 2*numTaxa - 2;
    strcpy (t->name, "");
    
    t->isRooted = isRooted;
    t->isClock = NO;

    t->checkConstraints = NO;
    t->nConstraints = 0;
    t->nLocks = 0;
    t->isCalibrated = NO;
    t->nNodes = t->nIntNodes = 0;
    t->nRelParts = 0;
    t->relParts = NULL;

    /* initialize pointers */
    t->bitsets = NULL;
    t->flags = NULL;
    t->constraints = NULL;

    /* allocate and initialize nodes and node arrays (enough for both rooted and unrooted trees) */
    if (t->isRooted)
        {
        t->nNodes = 2*numTaxa;
        t->nIntNodes = numTaxa - 1;
        }
    else
        {
        t->nNodes = 2*numTaxa - 2;
        t->nIntNodes = numTaxa - 2;
        }
    if ((t->nodes = (TreeNode *) SafeCalloc (t->nNodes, sizeof (TreeNode))) == NULL)
        {
        free (t);
        return NULL;
        }
    if ((t->allDownPass = (TreeNode **) SafeCalloc (t->nNodes + t->nIntNodes, sizeof (TreeNode *))) == NULL)
        {
        free (t->nodes);
        free (t);
        return NULL;
        }
    t->intDownPass = t->allDownPass + t->nNodes;
    
    /* initialize nodes and set index and memoryIndex */
    for (i=0; i<t->memNodes; i++)
        {
        ResetTreeNode(&t->nodes[i]);
        t->nodes[i].memoryIndex = i;
        t->nodes[i].index = i;
        }

    return t;
}


/* AllocateTreePartitions: Allocate space for and set partitions for tree */
int AllocateTreePartitions (Tree *t)
{
    int         i, nLongsNeeded, numTaxa;
    TreeNode    *p;
    
    /* get some handy numbers */
    if (t->isRooted == YES)
        numTaxa = t->nNodes - t->nIntNodes - 1;
    else
        numTaxa = t->nNodes - t->nIntNodes;
    nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;

    /* reallocate space */
    t->bitsets = (BitsLong *) SafeRealloc ((void *)(t->bitsets), (size_t)(t->nNodes) * (size_t)nLongsNeeded * sizeof(BitsLong));
    if (!t->bitsets)
        return (ERROR);
    
    /* clear bit fields */
    for (i=0; i<t->nNodes*nLongsNeeded; i++)
        t->bitsets[i] = 0;
        
    /* set node pointers to bit fields */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->partition = t->bitsets + i*nLongsNeeded;
        }

    /* set partition specifiers for terminals */
    ResetTreePartitions(t);
    
    return (NO_ERROR);
}


int AreTopologiesSame (Tree *t1, Tree *t2)
{
    int         i, j, k, nLongsNeeded, nTaxa;
    BitsLong    *mask;
    TreeNode    *p, *q;

    if (t1->nNodes != t2->nNodes)
        return (NO);
    if (t1->nIntNodes != t2->nIntNodes)
        return (NO);
    
    if (t1->isRooted == YES)
        nTaxa = t1->nNodes - t1->nIntNodes - 1;
    else
        nTaxa = t1->nNodes - t1->nIntNodes;
    
    /* allocate space for mask */
    nLongsNeeded = (nTaxa - 1) / nBitsInALong + 1;
    mask = (BitsLong *) SafeCalloc (nLongsNeeded, sizeof(BitsLong));
    
    /* set mask */
    for (i=0; i<nTaxa; i++)
        SetBit(i, mask);
    
    /* allocate and set partition pointers */
    AllocateTreePartitions (t1);
    AllocateTreePartitions (t2);

    /* check for congruence */
    for (i=0; i<t1->nIntNodes; i++)
        {
        p = t1->intDownPass[i];
        if (t1->isRooted == NO && IsBitSet(t2->root->index,p->partition))
                    FlipBits(p->partition,nLongsNeeded, mask);
        for (j=0; j<t2->nIntNodes; j++)
            {
            q = t2->intDownPass[j];
            for (k=0; k<nLongsNeeded; k++)
                {
                if (p->partition[k] != q->partition[k])
                    break;
                }
            if (k == nLongsNeeded)
                break;
            }
        if (j == t2->nIntNodes)
            {
            FreeTreePartitions (t1);
            FreeTreePartitions (t2);
            free (mask);
            return (NO);            
            }
        }

    FreeTreePartitions (t1);
    FreeTreePartitions (t2);
    free (mask);
    return (YES);
}


int AreTreesSame (Tree *t1, Tree *t2)
{
    int         i, j, k, nLongsNeeded, nTaxa;
    BitsLong    *mask;
    TreeNode    *p, *q;

    extern void ShowNodes(TreeNode*, int, int);
    
    if (t1->nNodes != t2->nNodes)
        return (NO);
    if (t1->nIntNodes != t2->nIntNodes)
        return (NO);
    
    if (t1->isRooted == YES)
        nTaxa = t1->nNodes - t1->nIntNodes - 1;
    else
        nTaxa = t1->nNodes - t1->nIntNodes;
    
    /* allocate space for mask */
    nLongsNeeded = (nTaxa - 1) / nBitsInALong + 1;
    mask = (BitsLong *) SafeCalloc (nLongsNeeded, sizeof(BitsLong));
    
    /* set mask */
    for (i=0; i<nTaxa; i++)
        SetBit(i, mask);

    /* allocate and set partition pointers */
    AllocateTreePartitions (t1);
    AllocateTreePartitions (t2);

    /* check for congruence */
    for (i=0; i<t1->nNodes; i++)
        {
        p = t1->allDownPass[i];
        if (p->anc == NULL && t1->isRooted == YES)
            continue;
        if (t1->isRooted == NO && IsBitSet(t2->root->index,p->partition))
            FlipBits(p->partition,nLongsNeeded, mask);
        for (j=0; j<t2->nNodes; j++)
            {
            q = t2->allDownPass[j];
            for (k=0; k<nLongsNeeded; k++)
                {
                if (p->partition[k] != q->partition[k])
                    break;
                }
            if (k == nLongsNeeded && AreDoublesEqual (p->length, q->length, 0.000001) == YES)
                break;
            else if (k == nLongsNeeded)
                {
                FreeTreePartitions (t1);
                FreeTreePartitions (t2);
                free (mask);
                return (NO);
                }
            }
        if (j == t2->nNodes)
            {
            FreeTreePartitions (t1);
            FreeTreePartitions (t2);
            free (mask);
            return (NO);            
            }
        }

    FreeTreePartitions (t1);
    FreeTreePartitions (t2);
    free (mask);
    return (YES);
}


/*----------------------------------------------------------------
|
|   BuildConstraintTree: Build constraint tree. The tree t is
|      needed only to hold information about constraints and
|      included taxa.
|
----------------------------------------------------------------*/
int BuildConstraintTree (Tree *t, PolyTree *pt, char **localTaxonNames)
{
    int             i, j, k, constraintId, nLongsNeeded, nextNode;
    BitsLong        *constraintPartition, *mask;
    PolyNode        *pp, *qq, *rr, *ss, *tt;
    
    pt->isRooted = t->isRooted;

    nLongsNeeded = (numLocalTaxa - 1) / nBitsInALong + 1;
    constraintPartition = (BitsLong *) SafeCalloc (2*nLongsNeeded, sizeof(BitsLong));
    if (!constraintPartition)
        {
        MrBayesPrint ("%s   Problems allocating constraintPartition in BuildConstraintTree", spacer);
        return (ERROR);
        }
    mask = constraintPartition + nLongsNeeded;

    /* calculate mask (needed to take care of unused bits when flipping partitions) */
    for (i=0; i<numLocalTaxa; i++)
        SetBit (i, mask);

    /* reset all nodes */
    for (i=0; i<2*numLocalTaxa; i++)
        {
        pp = &pt->nodes[i];
        pp->isDated = NO;
        pp->calibration = NULL;
        pp->age = -1.0;
        pp->isLocked = NO;
        pp->lockID = -1;
        pp->index = i;
        }

    /* build a bush */
    pt->root = &pt->nodes[numLocalTaxa];
    for (i=0; i<numLocalTaxa; i++)
        {
        pp = &pt->nodes[i];
        pp->index = i;
        pp->left = NULL;
        if (i == numLocalTaxa - 1)
            pp->sib = NULL;
        else
            pp->sib = &pt->nodes[i+1];
        pp->anc = pt->root;
        }
    pp = pt->root;
    pp->left = &pt->nodes[0];
    pp->anc = pp->sib = NULL;
    pt->nNodes = numLocalTaxa + 1;
    pt->nIntNodes = 1;

    /* make sure the outgroup is the right-most node */
    pt->nodes[localOutGroup].index = numLocalTaxa - 1;
    pt->nodes[numLocalTaxa - 1].index = localOutGroup;

    /* allocate and set partition specifiers in bush */
    GetPolyDownPass(pt);
    AllocatePolyTreePartitions(pt);

    /* set terminal taxon labels */
    for (i=0; i<pt->nNodes; i++)
        {
        pp = pt->allDownPass[i];
        if (pp->index < numLocalTaxa)
            strcpy (pp->label, localTaxonNames[pp->index]);
        }

    /* resolve the bush according to constraints */
    /* for now, satisfy all constraints */
    /* for now, bail out if constraints are not compatible */
    /* Eventually, we might want to be build a parsimony (WAB) or compatibility (WIB) matrix and
       draw a starting tree from the universe according to the score of the tree. A simple way of accomplishing
       approximately this is to use sequential addition, with probabilities in each step determined
       by the parsimony or compatibility score of the different possibilities. */ 
    nextNode = numLocalTaxa + 1;
    t->nLocks=0;
    for (constraintId=0; constraintId<numDefinedConstraints; constraintId++)
        {
        if (t->constraints[constraintId] == NO || definedConstraintsType[constraintId] != HARD)
            continue;

        /* initialize bits in partition to add; get rid of deleted taxa in the process */
        ClearBits(constraintPartition, nLongsNeeded);
        for (i=j=0; i<numTaxa; i++)
            {
            if (taxaInfo[i].isDeleted == YES)
                continue;
            if (IsBitSet(i, definedConstraint[constraintId]) == YES)
                SetBit(j, constraintPartition);
            j++;
            }
        assert (j == numLocalTaxa);
                
        /* make sure outgroup is outside constrained partition if the tree is unrooted */
        if (t->isRooted == NO && IsBitSet(localOutGroup, constraintPartition))
            FlipBits(constraintPartition, nLongsNeeded, mask);

        /* check that partition should be included */
        k = NumBits(constraintPartition, nLongsNeeded);
        if (k == 0)
            {
            MrBayesPrint ("%s   WARNING: Constraint '%s' refers only to deleted taxa\n", spacer, constraintNames[constraintId]);
            MrBayesPrint ("%s            and will be disregarded\n", spacer);
            t->constraints[constraintId] = NO;
            continue;
            }
        if (k == 1)
            {
            MrBayesPrint ("%s   WARNING: Constraint '%s' refers to a single tip and\n", spacer, constraintNames[constraintId]);
            MrBayesPrint ("%s            will be disregarded\n", spacer);
            t->constraints[constraintId] = NO;
            continue;
            }

        /* check if root in rooted tree (we allow this to enable inference of ancestral states) */
        if (k == numLocalTaxa && t->isRooted == YES)
            {
            if (pt->root->isLocked == YES) {
                MrBayesPrint ("%s   WARNING: Constraint '%s' is a duplicate of another constraint\n", spacer, constraintNames[constraintId]);
                MrBayesPrint ("%s            and will be ignored\n", spacer);
                t->constraints[constraintId] = NO;
                continue;
                }
            pt->root->isLocked = YES;
            pt->root->lockID = constraintId;
            t->nLocks++;
            continue;
            }

        /* check if interior root in unrooted tree (we allow this to enable inference of ancestral states) */
        if ((k == numLocalTaxa - 1 || k == numLocalTaxa) && t->isRooted == NO)
            {
            if (pt->root->isLocked == YES) {
                MrBayesPrint ("%s   WARNING: Constraint '%s' is a duplicate of another constraint\n", spacer, constraintNames[constraintId]);
                MrBayesPrint ("%s            and will be ignored\n", spacer);
                t->constraints[constraintId] = NO;
                continue;
                }
            pt->root->isLocked = YES;
            pt->root->lockID = constraintId;
            t->nLocks++;
            continue;
            }

        /* find first included terminal */
        k = FirstTaxonInPartition (constraintPartition, nLongsNeeded);
        for (i=0; pt->nodes[i].index != k; i++)
            ;
        pp = &pt->nodes[i];

        /* go down until node is not included in constraint */
        do {
            qq = pp;
            pp = pp->anc;       
        } while (IsPartNested(pp->partition, constraintPartition, nLongsNeeded));

        /* check that the node has not yet been included */
        for (i=0; i<nLongsNeeded; i++)
            {
            if (qq->partition[i] != constraintPartition[i])
                break;
            }
        if (i==nLongsNeeded)
            {
            MrBayesPrint ("%s   WARNING: Constraint '%s' is a duplicate of another constraint\n", spacer, constraintNames[constraintId]);
            MrBayesPrint ("%s            and will be ignored\n", spacer);
            t->constraints[constraintId] = NO;
            continue;
            }

        /* create a new node */
        tt = &pt->nodes[nextNode++];
        tt->anc = pp;
        tt->isLocked = YES;
        tt->lockID = constraintId;
        t->nLocks++;
        for (i=0; i<nLongsNeeded; i++)
            tt->partition[i] = constraintPartition[i];
        pt->nIntNodes++;
        pt->nNodes++;

        /* sort descendant nodes in two connected groups: included and excluded */
        /* if there is a descendant that overlaps (incompatible) then return error */
        rr = ss = NULL;
        qq = pp->left;
        do {
            if (IsPartNested(qq->partition, constraintPartition, nLongsNeeded))
                {
                if (ss != NULL)
                    ss->sib = qq;
                else
                    tt->left = qq;
                ss = qq;
                qq->anc = tt;
                }
            else if (IsPartCompatible(qq->partition, constraintPartition, nLongsNeeded))
                {
                if (rr != NULL)
                    rr->sib = qq;
                else
                    tt->sib = qq;
                rr = qq;
                }
            else
                {
                free (constraintPartition);
                return (ERROR);
                }
            qq = qq->sib;
            } while (qq != NULL);
        pp->left = tt;
        rr->sib = ss->sib = NULL;
        }
    
    /* relabel interior nodes */
    GetPolyDownPass(pt);
    for (i=0; i<pt->nIntNodes; i++)
        pt->intDownPass[i]->index = i + numLocalTaxa;

    /* exit */
    free (constraintPartition);
    FreePolyTreePartitions(pt);
    return (NO_ERROR);
}


/*----------------------------------------------
|
|   BuildRandomRTopology: Builds a random rooted
|      topology. Will set indices in t->nodes
|      such that they are from 0 to n-1 for the n tips
|      and from n to 2n-2 for the n-1 interior 
|      nodes. Last is root. Does not touch labels
|      of tips.
|
----------------------------------------------*/
int BuildRandomRTopology (Tree *t, RandLong *seed)
{
    int         i, j, nTips;
    TreeNode    *p, *q, *r;

    nTips = t->nNodes - t->nIntNodes - 1;
    
    for (i=0; i<t->nNodes; i++)
        {
        p = &t->nodes[i];
        p->index = i;
        p->left = p->right = p->anc = NULL;
        }

    /* connect the first two tip nodes */
    q = &t->nodes[0];
    r = &t->nodes[1];
    p = &t->nodes[nTips];
    q->anc = r->anc = p;
    p->left = q;
    p->right = r;
    q = &t->nodes[2*nTips-1];
    p->anc = q;
    q->left = p;

    for (i=2; i<nTips; i++)
        {
        q = &t->nodes[i];
        r = &t->nodes[i-2+nTips+1];
        q->anc = r;
        r->left = q;
        j = (int) (RandomNumber(seed) * (2 * i - 1));
        if (j < i)
            p = &t->nodes[j];
        else
            p = &t->nodes[j-i + nTips];
        r->right = p;
        r->anc = p->anc;
        if (p->anc != NULL)
            {
            if (p->anc->left == p)
                p->anc->left = r;
            else
                p->anc->right = r;
            }
        p->anc = r;
        }

    /* set root and get downpass */
    t->root = &t->nodes[2*nTips-1];
    GetDownPass (t);

    /* relabel interior nodes */
    for (i=0; i<t->nIntNodes; i++)
        t->intDownPass[i]->index = i+nTips;

    return (NO_ERROR);
}


/*----------------------------------------------
|
|   BuildRandomUTopology: Builds a random unrooted
|      topology. Assumes that indices are set
|      in t->nodes from 0 to n-1 for the n tips
|      and from n to 2n-3 for the n-2 interior 
|      nodes. Move the calculation root after
|      this routine to get the right root.
|
----------------------------------------------*/
int BuildRandomUTopology (Tree *t, RandLong *seed)
{
    int         i, j, nTips;
    TreeNode    *p, *q, *r;

    nTips = t->nNodes - t->nIntNodes;
    
    for (i=0; i<t->nNodes; i++)
        {
        p = &t->nodes[i];
        p->index = i;
        p->left = p->right = p->anc = NULL;
        }
    
    /* connect the first three nodes, assuming 0 is calc root */
    q = &t->nodes[1];
    r = &t->nodes[2];
    p = &t->nodes[nTips];
    q->anc = r->anc = p;
    p->left = q;
    p->right = r;
    q = &t->nodes[0];
    p->anc = q;
    q->left = p;

    for (i=3; i<nTips; i++)
        {
        q = &t->nodes[i];
        r = &t->nodes[i - 3 + nTips + 1];
        q->anc = r;
        r->left = q;
        j = (int) (RandomNumber(seed) * (2 * i - 3));
        if (j < i - 1)
            p = &t->nodes[j+1];
        else
            p = &t->nodes[j+1-i + nTips];
        r->right = p;
        r->anc = p->anc;
        if (p->anc->left == p)
            p->anc->left = r;
        else
            p->anc->right = r;
        p->anc = r;
        }

    t->root = &t->nodes[0];
    
    /* get downpass */
    GetDownPass (t);

    /* relabel interior nodes */
    for (i=0; i<t->nIntNodes; i++)
        t->intDownPass[i]->index = i+nTips;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   CheckConstraints: Check that tree complies with constraints
|
----------------------------------------------------------------*/
int CheckConstraints (Tree *t)
{
    int             a, i, j, k, nLongsNeeded;
    BitsLong        *constraintPartition, *mask;
    TreeNode        *p=NULL;
        
    if (t->checkConstraints == NO)
        return (NO_ERROR);

    /* allocate space */
    nLongsNeeded = (numLocalTaxa - 1) / nBitsInALong + 1;
    constraintPartition = (BitsLong *) SafeCalloc (2*nLongsNeeded, sizeof(BitsLong));
    if (!constraintPartition)
        {
        MrBayesPrint ("%s   Problems allocating constraintPartition in CheckConstraints", spacer);
        return (ERROR);
        }
    mask = constraintPartition + nLongsNeeded;

    /* set mask (needed to reset unused bits when flipping partitions) */
    for (i=0; i<numLocalTaxa; i++) 
      SetBit (i, mask); 
    
    if (AllocateTreePartitions(t) == ERROR)
        {
        MrBayesPrint ("%s   Problems allocating tree partitions in CheckConstraints", spacer);
        return (ERROR);
        }

    for (a=0; a<numDefinedConstraints; a++)
        {
        if (t->constraints[a] == NO  || definedConstraintsType[a] != HARD)
            continue;

        /* set bits in partition to check */
        ClearBits(constraintPartition, nLongsNeeded);
        for (j=k=0; j<numTaxa; j++)
            {
            if (taxaInfo[j].isDeleted == YES)
                continue;
            if (IsBitSet(j, definedConstraint[a]) == YES)
                SetBit(k, constraintPartition);
            k++;
            }

        /* make sure outgroup is outside constrained partition if unrooted tree */
        if (t->isRooted == NO && IsBitSet(localOutGroup, constraintPartition))
            FlipBits(constraintPartition, nLongsNeeded, mask);

        /* find the locked node */
        for (i=j=0; i<t->nNodes; i++)
            {
            if (t->allDownPass[i]->isLocked == YES && t->allDownPass[i]->lockID == a)
                {
                p = t->allDownPass[i];
                j++;
                }
            }
    
        if (j != 1)
            {
            MrBayesPrint ("%s   Tree has %d locks with id %d identifying constraint '%s'\n", spacer, j, a, constraintNames[a]);
            free (constraintPartition);
            FreeTreePartitions(t);
            return (ERROR);
            }

        /* check that locked node is correct */
        for (i=0; i<nLongsNeeded; i++)
            {
            if (p->partition[i] != constraintPartition[i]) 
                {
                MrBayesPrint ("%s   Lock %d is set for the wrong node [this is a bug]\n", spacer, a);
                free (constraintPartition);
                FreeTreePartitions(t);
                return (ERROR);
                }
            }
        }
    
    FreeTreePartitions (t);
    free (constraintPartition);
    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   CheckSetConstraints: Check and set tree constraints
|
----------------------------------------------------------------*/
int CheckSetConstraints (Tree *t)
{
    int             a, i, j, k, nLongsNeeded, foundIt, numLocks;
    BitsLong        *constraintPartition, *mask;
    TreeNode        *p;
        
    if (t->checkConstraints == NO)
        return (NO_ERROR);

    /* reset all existing locks, if any */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->isLocked = NO;
        p->lockID = -1;
        if (p->left != NULL)
            {
            p->calibration = NULL;
            p->isDated = NO;
            p->age = -1.0;
            }
        }

    /* allocate space */
    if (AllocateTreePartitions (t) == ERROR)
        {
        MrBayesPrint ("%s   Problems allocating tree bitsets", spacer);
        return ERROR;
        }

    nLongsNeeded = ((numLocalTaxa - 1) / nBitsInALong) + 1;
    constraintPartition = (BitsLong *) SafeCalloc (2*nLongsNeeded, sizeof(BitsLong));
    if (!constraintPartition)
        {
        MrBayesPrint ("%s   Problems allocating constraintPartition", spacer);
        FreeTreePartitions(t);
        return ERROR;
        }
    mask = constraintPartition + nLongsNeeded;

    /* set mask (needed to take care of unused bits when flipping partitions) */
    for (i=0; i<numLocalTaxa; i++)
        SetBit (i, mask);

    numLocks = 0;
    for (a=0; a<numDefinedConstraints; a++)
        {
        if (modelParams[t->relParts[0]].activeConstraints[a] == NO || definedConstraintsType[a] != HARD)
            continue;

        /* set bits in partition to add */
        ClearBits(constraintPartition, nLongsNeeded);
        for (i=j=0; i<numTaxa; i++)
            {
            if (taxaInfo[i].isDeleted == YES)
                continue;
            if (IsBitSet(i, definedConstraint[a]) == YES)
                SetBit(j, constraintPartition);
            j++;
            }

        /* make sure outgroup is outside constrained partition (marked 0) */
        if (t->isRooted == NO && IsBitSet(localOutGroup, constraintPartition) == YES)
            FlipBits(constraintPartition, nLongsNeeded, mask);

        /* skip partition if uninformative */
        k = NumBits(constraintPartition, nLongsNeeded);
        if (k == 0 || k == 1)
            continue;
            
        /* find the node that should be locked */
        foundIt = NO;
        for (i=0; i<t->nIntNodes; i++)
            {
            p = t->intDownPass[i];
            for (j=0; j<nLongsNeeded; j++)
                {
                if (p->partition[j] != constraintPartition[j])
                    break;
                }

            if (j == nLongsNeeded)
                {
                foundIt = YES;
                p->isLocked = YES;
                p->lockID = a;
                if (nodeCalibration[a].prior != unconstrained)
                    {
                    p->isDated = YES;
                    p->calibration = &nodeCalibration[a];
                    }
                numLocks++;
                break;
                }
            }
    
        if (foundIt == NO)
            {
            MrBayesPrint ("%s   Tree breaks constraint '%s'\n", spacer, constraintNames[a]);
            FreeTreePartitions (t);
            free (constraintPartition);
            return (ERROR);
            }
        }

    if (numLocks != t->nLocks)
        {
        MrBayesPrint ("%s   Inconsistent lock settings. This is a bug, please report it.\n", spacer);
        FreeTreePartitions (t);
        free (constraintPartition);
        return (ERROR);
        }
    
    /* exit */
    FreeTreePartitions(t);
    free (constraintPartition);
    return (NO_ERROR);
}


/*-----------------------------------------------------------------------
|
|   ColorClusters: Recursive function to color the clusters in a tree by
|      assigning numbers to them in their variable x
|
------------------------------------------------------------------------*/
void ColorClusters (TreeNode *p, int *index)
{
    if (p!=NULL)
        {
        if (p->isLocked == YES || p->anc == NULL || p->anc->anc == NULL)
            p->x = (++(*index));
        else
            p->x = p->anc->x;
        ColorClusters(p->left, index);
        ColorClusters(p->right, index);
        }
}


/* CopyPolyNodes: Copies everything except pointers and memoryIndex */
void CopyPolyNodes (PolyNode *p, PolyNode *q, int nLongsNeeded)
{
    p->index                  = q->index; 
    p->mark                   = q->mark;
    p->length                 = q->length;
    p->x                      = q->x;
    p->y                      = q->y;
    p->isDated                = q->isDated;
    p->calibration            = q->calibration;
    p->age                    = q->age;
    p->isLocked               = q->isLocked;
    p->lockID                 = q->lockID;
    strcpy (p->label, q->label);
    if (nLongsNeeded!=0)
        {
        assert (p->partition);
        assert (q->partition);
        memcpy (p->partition,q->partition, nLongsNeeded*sizeof(BitsLong));
        }
    p->support                = q->support;
    p->f                      = q->f;
}


void CopySubtreeToTree (Tree *subtree, Tree *t)
{
    int         i, /*j,*/ k;
    TreeNode    *p, *q=NULL, *r;

    for (i=/*j=*/0; i<subtree->nNodes - 1; i++)
        {
        p = subtree->allDownPass[i];

        for (k=0; k<t->nNodes; k++)
            {
            q = t->allDownPass[k];
            if (q->index == p->index)
                break;
            }
        q->length = p->length;
        q->marked = YES;
        if (p->left != NULL && p->right != NULL)
            {
            for (k=0; k<t->nNodes; k++)
                {
                r = t->allDownPass[k];
                if (r->index == p->left->index)
                    {
                    q->left = r;
                    r->anc = q;
                    }
                else if (r->index == p->right->index)
                    {
                    q->right = r;
                    r->anc = q;
                    }
                }
            }
        }

    p = subtree->root;

    for (k=0; k<t->nNodes; k++)
        {
        q = t->allDownPass[k];
        if (q->index == p->index)
            break;
        }

    if (q->left->marked == YES)
        {
        for (k=0; k<t->nIntNodes; k++)
            {
            r = t->intDownPass[k];
            if (r->index == p->left->index)
                {
                q->left = r;
                r->anc = q;
                }
            }
        }
    else if (q->right->marked == YES)
        {
        for (k=0; k<t->nIntNodes; k++)
            {
            r = t->intDownPass[k];
            if (r->index == p->left->index)
                {
                q->right = r;
                r->anc = q;
                }
            }
        }
}


/*-----------------------------------------------------------------
|
|   CopyToPolyTreeFromPolyTree: copies second tree to first tree
|
-----------------------------------------------------------------*/
int CopyToPolyTreeFromPolyTree (PolyTree *to, PolyTree *from)
{
    int         i, j, k, nLongsNeeded;
    PolyNode    *p, *q;

    /* check we have enough memory */
    assert (to->memNodes >= from->nNodes);
    if (from->bitsets==NULL || to->bitsets==NULL)
        {
        nLongsNeeded=0;
        }
    else
        {
        assert (to->memNodes >= from->memNodes);/*Otherwise partition length woould not be long enough for nodes in "to" */
        nLongsNeeded = (from->memNodes/2 - 1) / nBitsInALong + 1;
        }

    /* copy nodes */
    for (i=0; i<from->nNodes; i++)
        {
        /* copy pointers */
        p  = from->nodes + i;
        q  = to->nodes + i;

        if (p->anc != NULL)
            q->anc = to->nodes + p->anc->memoryIndex;
        else
            {
            q->anc = NULL;
            to->root = q;
            }

        if (p->left != NULL)
            q->left = to->nodes + p->left->memoryIndex;
        else
            q->left = NULL;

        if (p->sib != NULL)
            q->sib = to->nodes + p->sib->memoryIndex;
        else
            q->sib = NULL;

        /* Copy everything else except memoryIndex */
        CopyPolyNodes (q, p, nLongsNeeded);
        }

    /* fill node arrays */
    /* copy tree properties */
    to->nNodes = from->nNodes;
    to->nIntNodes = from->nIntNodes;
    to->isRooted = from->isRooted;
    to->isClock = from->isClock;
    to->isRelaxed = from->isRelaxed;
    to->clockRate = from->clockRate;
    strcpy (to->name, from->name);
  
    GetPolyDownPass (to);

    /* copy partitions */
    if (from->bitsets)
        {
        if (!to->bitsets)
            AllocatePolyTreePartitions(to);
        else
            ResetPolyTreePartitions(to);
        }

    /* copy relaxed clock parameters */
    FreePolyTreeRelClockParams (to);
    
    if (from->nBSets + from->nESets > 0)
        AllocatePolyTreeRelClockParams (to, from->nBSets, from->nESets);

    for (i=0; i<to->nBSets; i++)
        {
        to->bSetName[i] = (char *) SafeCalloc (strlen(from->bSetName[i])+2, sizeof(char));
        strcpy (to->bSetName[i], from->bSetName[i]);
        for (j=0; j<from->nNodes; j++)
            to->effectiveBrLen[i][j] = from->effectiveBrLen[i][j];
        }
    
    for (i=0; i<to->nESets; i++)
        {
        to->eSetName[i] = (char *) SafeCalloc (strlen(from->eSetName[i])+2, sizeof(char));
        strcpy (to->eSetName[i], from->eSetName[i]);
        for (j=0; j<from->nNodes; j++)
            {
            to->nEvents[i][j] = from->nEvents[i][j];
            if (to->nEvents[i][j] > 0)
                {
                to->position[i][j] = (MrBFlt *) SafeCalloc (to->nEvents[i][j], sizeof (MrBFlt));
                to->rateMult[i][j] = (MrBFlt *) SafeCalloc (to->nEvents[i][j], sizeof (MrBFlt));
                for (k=0; k<to->nEvents[i][j]; k++)
                    {
                    to->position[i][j][k] = from->position[i][j][k];
                    to->rateMult[i][j][k] = from->rateMult[i][j][k];
                    }
                }
            }
        }
    
    /* copy population size parameters */
    FreePolyTreePopSizeParams(to);
    to->popSizeSet = from->popSizeSet;
    if (to->popSizeSet == YES)
        {
        to->popSize = (MrBFlt *) SafeCalloc (to->nNodes, sizeof(MrBFlt));
        for (i=0; i<to->nNodes; i++)
            to->popSize[i] = from->popSize[i];
        to->popSizeSetName = (char *) SafeCalloc (strlen(from->popSizeSetName) + 1, sizeof(char));
        strcpy (to->popSizeSetName, from->popSizeSetName);
        }

    return (NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   CopyToSpeciesTreeFromPolyTree: copies second tree (polytomous) to
|       first tree, which is a species tree. The species tree needs to
|       be allocated enough space first to hold the resulting tree.
|
-----------------------------------------------------------------*/
int CopyToSpeciesTreeFromPolyTree (Tree *to, PolyTree *from)
{
    int         i;
    PolyNode    *p;
    TreeNode    *q, *q1;
#   if defined (DEBUG_SPECIESTREE)
    int         j;
#   endif

    /* make sure assumptions are correct */
    assert (from->isRooted == YES);
    assert (from->isClock == YES);
    assert (from->nNodes - from->nIntNodes == numSpecies);
    assert (to->memNodes == 2*numSpecies);
    assert (to->nIntNodes == from->nIntNodes);
    assert (to->nNodes == from->nNodes + 1);

    /* make sure indices are set correctly for from nodes */
#   if defined (DEBUG_SPECIESTREE)
    for (i=0; i<from->nNodes; i++)
        {
        for (j=0; j<from->nNodes; j++)
            {
            p = from->allDownPass[j];
            if (p->index == i)
                break;
            }
        assert (j != from->nNodes);
        assert (!(p->left == NULL && p->index >= numSpecies));
        }
#   endif

    /* copy nodes */
    for (i=0; i<from->nNodes; i++)
        {
        /* copy pointers */
        p  = from->allDownPass[i];
        q  = to->nodes + p->index;

        if (p->anc != NULL)
            q->anc = to->nodes + p->anc->index;
        else
            q->anc = NULL;

        if (p->left != NULL)    
            q->left = to->nodes + p->left->index;
        else
            q->left = NULL;

        if (p->left != NULL)
            q->right = to->nodes + p->left->sib->index;
        else
            q->right = NULL;

        q->nodeDepth              = p->depth;
        q->age                    = p->age;
        q->length                 = p->length;
        q->index                  = p->index;
        if (q->index < numSpecies)
            q->label = speciesNameSets[speciespartitionNum].names[q->index];
        else
            q->label = noLabel;
        }

    /* fix root */
    p = from->root;
    q = to->nodes + p->index;
    q1 = to->nodes + from->nNodes;      /* get the 'extra' root node that polytomous trees do not use */
    q->anc = q1;
    q1->index = from->nNodes;
    q1->left = q;
    q1->right = q1->anc = NULL;
    q1->isLocked = NO;
    q1->lockID = -1;
    q1->isDated = NO;
    q1->calibration = NULL;
    q1->age = -1.0;
    to->root = q1;

    /* get downpass */
    GetDownPass (to);
    
    /* a user tree might not come with node depths set */
    if (to->root->left->nodeDepth == 0.0)
        SetNodeDepths(to);

    /* set partitions */
    if (to->bitsets)
        ResetTreePartitions(to);

    return (NO_ERROR);      
}


/*-----------------------------------------------------------------
|
|   CopyToTreeFromPolyTree: copies second tree (polytomous) to first
|       tree (used to initialize constrained starting trees, e.g.).
|       An unrooted source tree will be rooted on outgroup
|       An unrooted source tree that needs to be copied to
|       a rooted target tree will be randomly rooted on a node below
|       all defined constraints. The to tree needs to be allocated
|       enough space first to hold the resulting tree.
|
-----------------------------------------------------------------*/
int CopyToTreeFromPolyTree (Tree *to, PolyTree *from)
{
    int         i, j;
    PolyNode    *p=NULL;
    TreeNode    *q, *q1;

    /* refuse to arbitrarily root an input tree */
    assert (!(from->isRooted == NO && to->isRooted == YES));
    if ((from->isRooted == NO) && (to->isRooted == YES))
        {
        MrBayesPrint ("%s   Failed to copy trees due to difference in rootedness of source and destination. \n", spacer);
        return (ERROR);
        }

    /* make sure assumptions are in order */
    assert (to->memNodes >= from->nNodes + (to->isRooted == NO ? 0 : 1));
    assert (localOutGroup >= 0 && localOutGroup < numLocalTaxa);
    assert (numLocalTaxa == from->nNodes - from->nIntNodes);
    assert (!(from->isRooted == YES && from->nNodes != 2*from->nIntNodes + 1));
    assert (!(from->isRooted == NO  && from->nNodes != 2*from->nIntNodes + 2));
    
    /* make sure indices are set correctly for from nodes */
    for (i=0; i<from->nNodes; i++)
        {
        for (j=0; j<from->nNodes; j++)
            {
            p = from->allDownPass[j];
            if (p->index == i)
                break;
            }
        assert (j != from->nNodes);
        assert (!(p->left == NULL && p->index >= numLocalTaxa));
        }
                
    /* deal with root */
    if (to->isRooted == NO && from->isRooted == YES)
        Deroot(from);

    /* make sure calculation root is set correctly */
    if (to->isRooted == NO && MovePolyCalculationRoot (from, localOutGroup) == ERROR)
        return ERROR;

    /* copy nodes */
    for (i=0; i<from->nNodes; i++)
        {
        /* copy pointers */
        p  = from->allDownPass[i];
        q  = to->nodes + p->index;

        if (p->anc != NULL)
            q->anc = to->nodes + p->anc->index;
        else
            q->anc = NULL;

        if (p->left != NULL)    
            q->left = to->nodes + p->left->index;
        else
            q->left = NULL;

        if (p->left != NULL)
            q->right = to->nodes + p->left->sib->index;
        else
            q->right = NULL;

        q->isLocked               = p->isLocked;
        q->lockID                 = p->lockID;
        q->isDated                = p->isDated;
        q->calibration            = p->calibration;
        q->age                    = p->age;
        q->nodeDepth              = p->depth;
        q->length                 = p->length;
        q->index                  = p->index;
        if (q->index < numLocalTaxa)
            q->label = localTaxonNames[q->index];
        else
            q->label = noLabel;
        }

    /* fix root */
    if (to->isRooted == NO)
        {
        p = from->root;
        q = to->nodes + p->index;
        q->anc = to->root = to->nodes + p->left->sib->sib->index;
        q->length = to->root->length;
        to->root->length = 0.0;
        to->root->left = q;
        to->root->right = to->root->anc = NULL;
        }
    else
        {
        p = from->root;
        q = to->nodes + p->index;
        q1 = to->nodes + from->nNodes;      /* get the 'extra' root node that polytomous trees do not use */
        q->anc = q1;
        q1->index = from->nNodes;
        q1->left = q;
        q1->right = q1->anc = NULL;
        q1->isLocked = NO;
        q1->lockID = -1;
        q1->isDated = NO;
        q1->calibration = NULL;
        q1->age = -1.0;
        to->root = q1;
        }

    /* get downpass */
    GetDownPass (to);
    
    /* set node depths */
    if (to->isRooted == YES && to->root->left->nodeDepth == 0.0)
        SetNodeDepths(to);

    /* set partitions */
    if (to->bitsets)
        ResetTreePartitions(to);

    /* relaxed clock parameters are not stored in binary trees but in separate parameters */

    return (NO_ERROR);      
}


/*-----------------------------------------------------------------
|
|   CopyToTreeFromTree: copies second tree to first tree
|       (used to initialize brlen sets for same topology)
|       Note: partition information of nodes are not copied if
|       either "from" or "to" tree does not have bitsets allocated
|
-----------------------------------------------------------------*/
int CopyToTreeFromTree (Tree *to, Tree *from)
{
    int         i, numTaxa, nLongsNeeded;
    TreeNode    *p, *q;

    numTaxa = from->nNodes - from->nIntNodes - (from->isRooted == YES ? 1 : 0);
    nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;
    if (from->bitsets==NULL || to->bitsets==NULL)
        nLongsNeeded=0;

    /* check that there is enough memory */
    assert (to->memNodes >= from->nNodes);
    
    /* copy nodes (use index of p as memoryIndex for q) */
    for (i=0; i<from->nNodes; i++)
        {
        /* copy pointers */
        p  = from->nodes + i;
        q  = to->nodes + p->index;

        if (p->anc != NULL)
            q->anc = to->nodes + p->anc->index;
        else
            {
            q->anc = NULL;
            to->root = q;
            }

        if (p->left != NULL)
            q->left = to->nodes + p->left->index;
        else
            q->left = NULL;

        if (p->right != NULL)
            q->right = to->nodes + p->right->index;
        else
            q->right = NULL;

        CopyTreeNodes (q, p, nLongsNeeded);
        }

    /* create new node arrays */
    to->nNodes = from->nNodes;
    to->nIntNodes = from->nIntNodes;
    GetDownPass (to);

    /* copy tree properties (these should be constant most of them) */
    strcpy (to->name, from->name);
    to->isRooted = from->isRooted;
    to->isClock = from->isClock;
    to->isCalibrated = from->isCalibrated;
    to->checkConstraints = from->checkConstraints;
    to->nConstraints = from->nConstraints;
    to->constraints = from->constraints;
    to->nLocks = from->nLocks;
    to->nRelParts = from->nRelParts;
    to->relParts = from->relParts;

    /* copy partitions */
    if (from->bitsets)
        {
        if (!to->bitsets)
            AllocateTreePartitions(to);
        else
            ResetTreePartitions(to);
        }

    return (NO_ERROR);
}


/* Copy node q to node p */
void CopyTreeNodes (TreeNode *p, TreeNode *q, int nLongsNeeded)
{
    /* copies everything except pointers and memoryIndex */
    p->label                  = q->label;
    p->index                  = q->index;
    p->upDateCl               = q->upDateCl;
    p->upDateTi               = q->upDateTi;
    p->isLocked               = q->isLocked;
    p->lockID                 = q->lockID;
    p->isDated                = q->isDated;
    p->marked                 = q->marked;
    p->x                      = q->x;
    p->y                      = q->y;
    p->d                      = q->d;
    p->length                 = q->length;
    p->nodeDepth              = q->nodeDepth;
    p->calibration            = q->calibration;
    p->age                    = q->age;
    if (nLongsNeeded != 0)
        {
        assert (p->partition);
        assert (q->partition);
        memcpy (p->partition, q->partition, nLongsNeeded*sizeof(BitsLong));
        }
}


void CopyTreeToSubtree (Tree *t, Tree *subtree)
{
    int         i, j, k;
    TreeNode    *p, *q, *r;

    for (i=j=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO)
            continue;

        q = &subtree->nodes[j++];
        q->index = p->index;
        q->length = p->length;
        if (p->left == NULL || p->left->marked == NO)
            q->left = q->right = NULL;
        else
            {
            for (k=0; k<j-1; k++)
                {
                r = &subtree->nodes[k];
                if (r->index == p->left->index)
                    {
                    q->left = r;
                    r->anc = q;
                    }
                else if (r->index == p->right->index)
                    {
                    q->right = r;
                    r->anc = q;
                    }
                }
            }
        
        if (p->anc->marked == NO)
            {
            r = &subtree->nodes[j++];
            subtree->root = r;
            r->anc = r->right = NULL;
            r->left = q;
            q->anc = r;
            r->length = 0.0;
            r->index = p->anc->index;
            }

        }

    GetDownPass (subtree);

    subtree->isRooted = t->isRooted;
    subtree->nRelParts = t->nRelParts;
    subtree->relParts = t->relParts;
}


/* DatedNodeDepths: Recursive function to get node depths */
void DatedNodeDepths (TreeNode *p, MrBFlt *nodeDepths, int *index)
{
    if (p != NULL)
        {
        if (p->left == NULL || p->isDated == YES)
            nodeDepths[(*index)++] = p->nodeDepth;
        else
            {
            DatedNodeDepths (p->left,  nodeDepths, index);
            DatedNodeDepths (p->right, nodeDepths, index);
            }
        }
}


/* DatedNodes: Recursive function to get dated tips or interior nodes */
void DatedNodes (TreeNode *p, TreeNode **datedNodes, int *index)
{
    if (p != NULL)
        {
        if (p->left != NULL && p->isDated == NO)
            {
            DatedNodes (p->left,  datedNodes, index);
            DatedNodes (p->right, datedNodes, index);
            }
        datedNodes[(*index)++] = p;
        }
}


/* Deroot: Deroot a rooted polytomous tree with branch lengths */
int Deroot (PolyTree *pt)
{
    PolyNode    *p, *q, *r, tempNode;
    int         i;

    p = pt->root;

    if (p->left->sib->sib != NULL)
        return (ERROR);      /* tree is not rooted or it is polytomous */

    if (p != &pt->nodes[pt->nNodes-1])
        {
        q = &pt->nodes[pt->nNodes-1];
        /* now swap content of p and q including pointers */
        tempNode = *q;
        *q = *p;
        *p = tempNode;
        /* swap back memoryindex */
        p->memoryIndex = q->memoryIndex;
        q->memoryIndex = tempNode.memoryIndex;
        /* all pointers to q should be pointers to p */
        for (i=0; i<pt->nNodes; i++)
            {
            r = &pt->nodes[i];
            if (r->left == q)
                r->left = p;
            if (r->anc == q)
                r->anc = p;
            if (r->sib == q)
                r->sib = p;
            }
        /* all pointers to p should be pointers to q; all these are anc pointers from the descendants of the root */
        pt->root = q;
        for (r=q->left; r!=NULL; r=r->sib)
            r->anc = q;
        /* finally set p to the new root */
        p = pt->root;
        }

    /* make sure the left of the old root is interior and can be used as new root */
    if (p->left->left == NULL)
        {
        q = p->left;
        r = q->sib;
        p->left = r;
        r->sib = q;
        q->sib = NULL;
        }
    
    pt->root = p->left;
    pt->root->left->sib->sib = p->left->sib;
    p->left->sib->length += pt->root->length;
    pt->root->length = 0.0;
    pt->root->sib = NULL;
    pt->root->anc = NULL;

    pt->nNodes--;
    pt->nIntNodes--;

    GetPolyDownPass(pt);

    return (NO_ERROR);
}


/* EraseTreeList: Erase all trees in treeList */
void EraseTreeList (TreeList *treeList)
{
    TreeListElement *listElement;
    TreeListElement *previous;

    listElement = treeList->first;
    if (listElement != NULL)
        do 
            {
            free (listElement->order);
            previous = listElement;
            listElement = listElement->next;
            free (previous);
            } 
        while (listElement != NULL);

    treeList->first = treeList->last = NULL;
}


void UpdateTreeWithClockrate (Tree *t, MrBFlt clockRate)
{
    int i;
    TreeNode *p;

    if (t->fromUserTree == NO)
        {
        /*Set nodeDepth*/
        for (i=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            p->nodeDepth = p->age * clockRate;
            }

        /* calculate branch lengths */
        for (i=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            if (p->anc != NULL)
                {
                if (p->anc->anc != NULL)
                    {
                    p->length = p->anc->nodeDepth - p->nodeDepth;
                    }
                else
                    p->length = 0.0; //not a problem for root node. 
                }
            }
        }
    else
        {
        for (i=0; i<t->nNodes-1; i++)
            {
            p = t->allDownPass[i];
            p->age = p->nodeDepth / clockRate;
            }
        }
}


/*----------------------------------------------------------------
|
|   findAllowedClockrate: Finds the range of clock rates allowed for the tree.
|
|   @param t        - tree to check (IN)  
|   @minClockRate   - adress where minimum allowed clock rate is stored (OUT)
|   @maxClockRate   - adress where maximum allowed clock rate is stored (OUT)
|
----------------------------------------------------------------*/
void findAllowedClockrate (Tree *t, MrBFlt *minClockRate, MrBFlt *maxClockRate)
{
    int i;
    TreeNode *p;
    MrBFlt min, max, tmp;

    min=0.0;
    max=MRBFLT_MAX;

    *minClockRate = 2.0;
    *maxClockRate = 1.0;

    if (t->fromUserTree == NO)
        {
        for (i=0; i<t->nNodes-1; i++)
            {
            p = t->allDownPass[i];
            if (p->anc->anc != NULL)
                {
                tmp = BRLENS_MIN/(p->anc->age - p->age);
                assert (tmp > 0);
                if (tmp > min)
                    min = tmp;

                tmp = BRLENS_MAX/(p->anc->age - p->age);
                assert (tmp > 0);
                if (tmp > max)
                    max = tmp;
                }
            }
        *minClockRate= min;
        *maxClockRate= max;
        }
    else
        {
        IsCalibratedClockSatisfied (t,minClockRate,maxClockRate, 0.001);
        }
}


/* FreePolyTree: Free memory space for a polytomous tree (unrooted or rooted) */
void FreePolyTree (PolyTree *pt)
{   
    if (pt != NULL)
        {
        FreePolyTreePartitions(pt);
        FreePolyTreeRelClockParams(pt);
        FreePolyTreePopSizeParams(pt);
        free (pt->allDownPass);
        free (pt->nodes);
        free (pt);
        }
}


/* FreePolyTreePartitions: Free memory space for polytomous tree partitions */
void FreePolyTreePartitions (PolyTree *pt)
{
    int i;
    if (pt != NULL && pt->bitsets != NULL)
        {
        for (i=0; i<pt->memNodes; i++)
            pt->nodes[i].partition = NULL;
        free (pt->bitsets);
        pt->bitsets = NULL;
        }
}


/* FreePolyTreePopSizeParams: Free population size set parameters of polytree */
void FreePolyTreePopSizeParams (PolyTree *pt)
{
    if (pt->popSizeSet == YES)
        {
        free (pt->popSize);
        free (pt->popSizeSetName);
        }
    pt->popSizeSet = NO;
    pt->popSize = NULL;
    pt->popSizeSetName = NULL;
}


/* FreePolyTreeRelClockParams: Free relaxed clock parameters of polytree */
void FreePolyTreeRelClockParams (PolyTree *pt)
{
    int i, j;

    /* free breakpoint clock parameters */
    for (i=0; i<pt->nESets; i++)
        {
        for (j=0; j<pt->memNodes; j++)
            {
            if (pt->nEvents[i][j] > 0)
                {
                free (pt->position[i][j]);
                free (pt->rateMult[i][j]);
                }
            }
        free (pt->eSetName[i]);
        free (pt->nEvents[i]);
        free (pt->position[i]);
        free (pt->rateMult[i]);
        }
    free (pt->nEvents);
    free (pt->position);
    free (pt->rateMult);
    free (pt->eSetName);
    pt->nESets = 0;
    pt->nEvents = NULL;
    pt->position = NULL;
    pt->rateMult = NULL;
    pt->eSetName = NULL;

    /* free branch clock parameters */
    for (i=0; i<pt->nBSets; i++)
        {
        free (pt->bSetName[i]);
        free (pt->effectiveBrLen[i]);
        }
    free (pt->effectiveBrLen);
    free (pt->bSetName);
    pt->nBSets = 0;
    pt->effectiveBrLen = NULL;
    pt->bSetName = NULL;
}


/* FreeTree: Free memory space for a tree (unrooted or rooted) */
void FreeTree (Tree *t)
{
    if (t != NULL)
        {
        free (t->bitsets);
        free (t->flags);
        free (t->allDownPass);
        free (t->nodes);
        free (t);
        }
}


/* FreeTreePartitions: Free memory space for tree partitions */
void FreeTreePartitions (Tree *t)
{
    int     i;

    if (t != NULL && t->bitsets != NULL)
        {
        free (t->bitsets);
        t->bitsets = NULL;
        for (i=0; i<t->memNodes; i++)
            t->nodes[i].partition = NULL;
        }
}


/*-------------------------------------------------------------------------------------------
|
|   GetDatedNodeDepths: Get an array containing the node depths of the dated tips,
|       internal or external, plus dated root
|
---------------------------------------------------------------------------------------------*/
void GetDatedNodeDepths (TreeNode *p, MrBFlt *nodeDepths)
{
    int index = 0;
    
    assert (p != NULL);

    nodeDepths[index++] = p->nodeDepth;     /* include root node depth */
    if (p->left != NULL)
        {
        DatedNodeDepths (p->left, nodeDepths, &index);
        DatedNodeDepths (p->right, nodeDepths, &index);
        }
}


/*-------------------------------------------------------------------------------------------
|
|   GetDatedNodes: Get an array containing the dated tips,
|       internal or external, and all interior nodes in the same subtree
|
---------------------------------------------------------------------------------------------*/
void GetDatedNodes (TreeNode *p, TreeNode **datedNodes)
{
    int     index = 0;
    
    assert (p != NULL);

    if (p->left!= NULL)
        {
        DatedNodes (p->left,  datedNodes, &index);
        DatedNodes (p->right, datedNodes, &index);
        }
}


/* get down pass for tree t (wrapper function) */
void GetDownPass (Tree *t)
{
    int i, j;

    i = j = 0;
    GetNodeDownPass (t, t->root, &i, &j);
}


/* get the actual down pass sequences */
void GetNodeDownPass (Tree *t, TreeNode *p, int *i, int *j)
{
    if (p != NULL)
        {
        GetNodeDownPass (t, p->left,  i, j);
        GetNodeDownPass (t, p->right, i, j);
        if (p->left != NULL && p->right != NULL && p->anc != NULL)
            {
            t->intDownPass[(*i)++] = p;
            t->allDownPass[(*j)++] = p;
            }
        else if (p->left == NULL && p->right == NULL && p->anc != NULL)
            {
            t->allDownPass[(*j)++] = p;
            }
        else if (p->left != NULL && p->right == NULL && p->anc == NULL)
            {
            t->allDownPass[(*j)++] = p;
            }
        }
}


/* GetPolyAges: Get PolyTree node ages */
void GetPolyAges (PolyTree *t)
{
    int         i;
    PolyNode    *p;

    GetPolyDepths (t); /* just to make sure... */
    
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->age = p->depth / t->clockRate;
        }
}


/* GetPolyDepths: Get PolyTree node depths */
void GetPolyDepths (PolyTree *t)
{
    int         i;
    MrBFlt      maxDepth;
    PolyNode    *p;

    maxDepth = t->root->depth = 0.0;

    for (i=t->nNodes-2; i>=0; i--)
        {
        p = t->allDownPass[i];
        p->depth = p->anc->depth + p->length;
        if (p->depth > maxDepth)
            maxDepth = p->depth;
        }

    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->depth = maxDepth - p->depth;
        }
}


/* get down pass for polytomous tree t (wrapper function) */
void GetPolyDownPass (PolyTree *t)
{
    int i, j;

    i = j = 0;
    GetPolyNodeDownPass (t, t->root, &i, &j);
    assert (t->nIntNodes==j);
}


/* get the actual down pass sequences for a polytomous tree */
void GetPolyNodeDownPass (PolyTree *t, PolyNode *p, int *i, int *j)
{
    PolyNode    *q;
    
    if (p->left != NULL)
        {
        for (q=p->left; q!=NULL; q=q->sib)
            GetPolyNodeDownPass(t, q, i, j);
        }

    t->allDownPass[(*i)++] = p;
    if (p->left != NULL)
        t->intDownPass[(*j)++] = p;
}


/* GetFromTreeList: Get first tree from a tree list and remove it from the list*/
int GetFromTreeList (TreeList *treeList, Tree *tree)
{
    TreeListElement *listElement;

    if (treeList->first == NULL)
        {
        MrBayesPrint ("%s   Tree list empty\n", spacer);
        return (ERROR);
        }
    if (tree->isRooted == YES)
        RetrieveRTopology (tree, treeList->first->order);
    else
        {
        RetrieveUTopology (tree, treeList->first->order);
        if (localOutGroup != 0)
            MoveCalculationRoot (tree, localOutGroup);
        }

    listElement = treeList->first;
    treeList->first = listElement->next;

    free (listElement->order);
    free (listElement);

    return (NO_ERROR);
}


/*------------------------------------------------------------------
|
|   InitBrlens: This routine will set all branch lengths of a
|      nonclock tree to the value given by 'v'.
|
------------------------------------------------------------------*/
int InitBrlens (Tree *t, MrBFlt v)
{
    int         i;
    TreeNode    *p;

    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->anc != NULL && !(t->isRooted == YES && p->anc->anc == NULL))
            p->length = v;
        else
            p->length = 0.0;
        }

    return (NO_ERROR);
}


/* 
@param levUp        is the number of edges between the "node" and the most resent calibrated predecessor +1,
                    for the calibrated nodes it should be 1
@param calibrUp     is the age of the most resent calibrated predecessor
@return             age of the node 
*/
MrBFlt SetNodeCalibratedAge(TreeNode *node, unsigned levUp, MrBFlt calibrUp)
{
    MrBFlt r,l;

    if (node->age != -1.0)
        {
        if (node->right != NULL)
            SetNodeCalibratedAge (node->right, 2, node->age);
        if (node->left != NULL)
            SetNodeCalibratedAge (node->left,  2, node->age);
        return node->age;
        }

    r = SetNodeCalibratedAge (node->right, levUp+1, calibrUp);
    l = SetNodeCalibratedAge (node->left,  levUp+1, calibrUp);

    if (r > l)
        {
        assert (calibrUp - r > 0.0);
        return node->age = r + (calibrUp - r)/levUp;
        }
    else
        {
        assert (calibrUp - l > 0.0);
        return node->age = l + (calibrUp - l)/levUp;
        }
}


/*-------------------------------------------------------------------
|
|   InitCalibratedBrlens: This routine will build a clock tree
|       consistent with calibration constraints on terminal
|       taxa and/or constrained interior nodes. At least one
|       node should be calibrated.
|       If not possible to build such a tree, ERROR
|       is returned.
|
--------------------------------------------------------------------*/
int InitCalibratedBrlens (Tree *t, MrBFlt clockRate, RandLong *seed)
{
    int             i;
    TreeNode        *p;
    Model           *mp;
    MrBFlt          treeAgeMin, treeAgeMax;
    Calibration     *calibrationPtr;

#   ifdef DEBUG_CALIBRATION
    printf ("Before initializing calibrated brlens\n");
    ShowNodes(t->root, 0, YES);
#   endif
    
    if (t->isRooted == NO)
        {
        MrBayesPrint ("%s   Tree is unrooted\n", spacer);
        return (ERROR);
        }

    /* Check whether root has age constraints */
    mp = &modelParams[t->relParts[0]];
    treeAgeMin = 0.0;
    treeAgeMax = POS_INFINITY;
    if (t->root->left->isDated == YES)
        {
        treeAgeMin = t->root->left->calibration->min;
        treeAgeMax = t->root->left->calibration->max;
        }
    else if (!strcmp(mp->clockPr, "Uniform") || !strcmp(mp->clockPr, "Fossilization"))
        {
        if (mp->treeAgePr.min > treeAgeMin)
            treeAgeMin = mp->treeAgePr.min;
        if (mp->treeAgePr.max < treeAgeMax)
            treeAgeMax = mp->treeAgePr.max;
        }

    /* date all nodes from top to bottom with min. age as nodeDepth*/
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->anc != NULL)
            {
            if (p->left == NULL && p->right == NULL)
                {
                if (p->isDated == NO)
                    {
                    p->nodeDepth = 0.0;
                    p->age = 0.0;
                    }
                else
                    {
                    if (p->calibration->prior == fixed)
                        p->nodeDepth = p->age = p->calibration->priorParams[0];
                    else
                        p->nodeDepth = p->age = p->calibration->min;
                    }
                }
            else
                {
                if (p->left->nodeDepth > p->right->nodeDepth)
                    p->nodeDepth = p->left->nodeDepth;
                else
                    p->nodeDepth = p->right->nodeDepth;
                if (p->isDated == YES || (p->anc->anc == NULL && (!strcmp(mp->clockPr,"Uniform") || !strcmp(mp->clockPr,"Fossilization"))))
                    {
                    if (p->isDated == NO)
                        calibrationPtr = &mp->treeAgePr;
                    else
                        calibrationPtr = p->calibration;

                    if (calibrationPtr->max <= p->nodeDepth)
                        {
                        if (p->isDated == NO)
                            MrBayesPrint ("%s   Calibration inconsistency for root node\n", spacer);
                        else
                            MrBayesPrint ("%s   Calibration inconsistency for node '%s'\n", spacer, constraintNames[p->lockID]);
                        MrBayesPrint ("%s   Cannot make a tree where the node is %s\n", spacer, calibrationPtr->name);
                        return (ERROR);
                        }
                    else
                        {
                        if (calibrationPtr->min < p->nodeDepth)
                            p->age = p->nodeDepth;
                        else
                            p->age = p->nodeDepth = calibrationPtr->min;
                        }
                    }
                else
                    p->age = -1.0;
                }
            }
        }

    /* try to make root node deeper than minimum age */
    p = t->root->left;
    if (p->nodeDepth==0.0)  p->nodeDepth = 1.0;
    if (p->nodeDepth * 1.5 < treeAgeMax)
        p->nodeDepth = p->age = 1.5 * p->nodeDepth;
    else
        p->nodeDepth = p->age = treeAgeMax;

    SetNodeCalibratedAge (p, 1, p->age);

    /* Setup node depths */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->nodeDepth = p->age * clockRate;
        assert (!(p->left == NULL && p->calibration == NULL && p->nodeDepth != 0.0));
        }

    /* calculate branch lengths */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->anc != NULL)
            {
            if (p->anc->anc != NULL)
                {
                p->length = p->anc->nodeDepth - p->nodeDepth;
                if (p->length < BRLENS_MIN)
                    {
                    //MrBayesPrint ("%s   Restrictions of node calibration and clockrate makes some branch lenghts too small.\n", spacer);
                    //return (ERROR);
                    }
                if (p->length > BRLENS_MAX)
                    {
                    //MrBayesPrint ("%s   Restrictions of node calibration and clockrate makes some branch lenghts too long.\n", spacer);
                    //return (ERROR);
                    }
                }
            else
                p->length = 0.0; //not a problem for root node. 
            }
        }

#   ifdef DEBUG_CALIBRATION
    printf ("after\n");
    ShowNodes (t->root, 0, YES);
    getchar();
#   endif

    return (NO_ERROR);
    MrBayesPrint ("%lf", *seed); /* just because I am tired of seeing the unused parameter error msg */
}


/*-------------------------------------------------------
|
|   InitClockBrlens: This routine will initialize
|      a clock tree by setting the root to depth 1.0
|      and then assigning node depths according to
|      the relative node depth measured in terms of the
|      maximum number of branches to the tip from each
|      node.
|
--------------------------------------------------------*/
int InitClockBrlens (Tree *t)
{
    int             i, maxBrSegments=0;
    TreeNode        *p;

    if (t->isRooted == NO)
        {
        MrBayesPrint ("%s   Tree is unrooted\n", spacer);
        return (ERROR);
        }
    
    /* calculate maximum number of branch segments above root */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->anc != NULL)
            {
            if (p->left == NULL && p->right == NULL)
                {
                p->x = 0;
                }
            else
                {
                if (p->left->x > p->right->x)
                    p->x = p->left->x + 1;
                else
                    p->x = p->right->x + 1;
                }
            if (p->anc->anc == NULL)
                maxBrSegments = p->x;
            }
        }

    /* assign node depths */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->anc != NULL)
            p->nodeDepth = (MrBFlt) (p->x) / (MrBFlt) maxBrSegments;
        else
            p->nodeDepth = 0.0;
        }
        
    /* calculate branch lengths */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->anc != NULL)
            {
            if (p->anc->anc != NULL)
                p->length = p->anc->nodeDepth - p->nodeDepth;
            else
                p->length = 0.0;
            }
        }

    return (NO_ERROR);
}


int GetRandomEmbeddedSubtree (Tree *t, int nTerminals, RandLong *seed, int *nEmbeddedTrees)
{
    int         i, j, k, n, ran, *pP, *pL, *pR, nLeaves, *nSubTrees;
    TreeNode    *p=NULL, **leaf;

    /* Calculate number of leaves in subtree (number of terminals minus the root) */
    nLeaves = nTerminals - 1;
    
    /* Initialize all flags */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;
        p->x = 0;
        p->y = 0;
        }
    
    /* Allocate memory */
    nSubTrees = (int *) SafeCalloc (nTerminals * t->nNodes, sizeof(int));
    if (!nSubTrees)
        return (ERROR);
    leaf = (TreeNode **) SafeMalloc (nLeaves * sizeof (TreeNode *));
    if (!leaf)
        {
        free (nSubTrees);
        return (ERROR);
        }

    /* Calculate how many embedded trees rooted at each node */
    (*nEmbeddedTrees) = 0;
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            {
            p->x = 0;
            nSubTrees[p->index*nTerminals + 1] = 1;
            }
        else
            {
            pL = nSubTrees + p->left->index*nTerminals;
            pR = nSubTrees + p->right->index*nTerminals;
            pP = nSubTrees + p->index*nTerminals;
            pP[1] = 1;
            for (j=2; j<=nLeaves; j++)
                {
                for (k=1; k<j; k++)
                    {
                    pP[j] += pL[k] * pR[j-k];
                    }
                }
            p->x = pP[nLeaves];
            (*nEmbeddedTrees) += p->x;
            }
        }

    /* Randomly select one embedded tree of the right size */
    ran = (int) (RandomNumber(seed) * (*nEmbeddedTrees));

    /* Find the interior root corresponding to this tree */
    for (i=j=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        j += p->x;
        if (j>ran)
            break;
        }

    /* Find one random embedded tree with this root */
    p->y = nLeaves;
    p->marked = YES;
    leaf[0] = p;
    n = 1;
    while (n < nLeaves)
        {
        /* select a node with more than one descendant */
        for (i=0; i<n; i++)
            {
            p = leaf[i];
            if (p->y > 1)
                break;
            }

        /* break it into descendants */
        pL = nSubTrees + p->left->index*nTerminals;
        pR = nSubTrees + p->right->index*nTerminals;
        pP = nSubTrees + p->index*nTerminals;
        ran = (int) (RandomNumber(seed) * pP[p->y]);
        k = 0;
        for (j=1; j<p->y; j++)
            {
            k += pL[j] * pR[p->y-j];
            if (k > ran)
                break;
            }

            p->left->y = j;
        p->right->y = p->y - j;
        p->left->marked = YES;
        p->right->marked = YES;
        leaf[i] = p->left;
        leaf[n++] = p->right;
        }

    free (nSubTrees);
    free (leaf);

    return (NO_ERROR);
}

        
/*-----------------------------------------------------------------------------
|
| IsCalibratedClockSatisfied: This routine SETS (not just checks as name suggested) calibrated clock tree nodes age, depth. based on branch lengthes
|     and checks that user defined brlens satisfy the specified calibration(s) up to tolerance tol
| TODO: clock rate is devived here and used to set ages but clockrate parameter is not updated here (make sure that it does not produce inconsistancy)
|
|------------------------------------------------------------------------------*/
int IsCalibratedClockSatisfied (Tree *t,MrBFlt *minClockRate,MrBFlt *maxClockRate , MrBFlt tol)
{
    int             i, j, maxRateConstrained, minRateConstrained, isViolated;
    MrBFlt          f, maxHeight, minRate=0, maxRate=0, ageToAdd, *x, *y, clockRate;
    TreeNode        *p, *q, *r, *s;

    /* By defauult assume the tree does not have allowed range of clockrate */
    *minClockRate = 2.0;
    *maxClockRate = 1.0;

    if (t->isRooted == NO)
        return (NO);
        
    x = (MrBFlt *) SafeCalloc (2*t->nNodes, sizeof (MrBFlt));
    if (x == NULL)
        {
        MrBayesPrint ("%s   Out of memory in IsCalibratedClockSatisfied\n", spacer);
        free (x);
        return (NO);
        }
    y = x + t->nNodes;

    /* reset node depth and age, and set minimum (x) and maximum (y) age of each node */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->age = -1.0;
        p->nodeDepth = -1.0;
        if (p->isDated == YES)
            {
            assert (p->calibration->prior != unconstrained);
            x[p->index] = p->calibration->min;
            y[p->index] = p->calibration->max;
            }
        else if (p->left == NULL && p->right == NULL)
            x[p->index] = y[p->index] = 0.0;
        else
            {
            x[p->index] = y[p->index] = -1.0;
            }
        }

    /* calculate node heights in branch length units */
    /* node depth will be set from the root for now  */
    p = t->root->left;
    p->nodeDepth = 0.0;
    for (i=t->nNodes-3; i>=0; i--)
        {
        p = t->allDownPass[i];
        p->nodeDepth = p->anc->nodeDepth + p->length;
        }

    /* find maximum height of tree */   
    maxHeight = -1.0;
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL && p->right == NULL)
            {
            if (p->nodeDepth > maxHeight)
                {
                maxHeight = p->nodeDepth;
                }
            }
        }
    
    /* calculate node depth from tip of tree */
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        p->nodeDepth = maxHeight - p->nodeDepth;
        }

    /* check potentially constraining calibrations */
    /* and find minimum and maximum possible rate */
    maxRateConstrained = NO;
    minRateConstrained = NO;
    isViolated = NO;
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (x[p->index] < 0.0 && y[p->index] < 0.0)
            continue;
        for (j=i+1; j<t->nNodes-1; j++)
            {
            q = t->allDownPass[j];
            if (x[q->index] < 0.0 && y[q->index] < 0.0)
                continue;
            if (p->nodeDepth == q->nodeDepth) // becouse clock rate could be as low as possible we can not take approximate equality. 
                {
                /* same depth so they must share a possible age */
                if ((x[p->index] != -1.0 && y[q->index] !=-1.0 && AreDoublesEqual (x[p->index], y[q->index], tol) == NO && x[p->index] > y[q->index]) ||
                    (y[p->index] != -1.0 && x[q->index] !=-1.0 && AreDoublesEqual (y[p->index], x[q->index], tol) == NO && y[p->index] < x[q->index]))
                    {
                    isViolated = YES;
                    break;
                    }
                }
            else
                {
                if (p->nodeDepth > q->nodeDepth)
                    {
                    r = p;
                    s = q;
                    }
                else
                    {
                    r = q;
                    s = p;
                    }
                if (x[r->index] >= 0.0 && y[s->index] >= 0.0)
                    {
                    f = (r->nodeDepth - s->nodeDepth) / (x[r->index] - y[s->index]);
                    if (f <= 0.0 || x[r->index] == y[s->index])
                        {
                        if (AreDoublesEqual (r->nodeDepth, s->nodeDepth, tol*0.1) == YES)
                            continue;
                        if ((r->calibration != NULL && r->calibration->prior != fixed) || (s->calibration != NULL && s->calibration->prior != fixed))
                            continue;
                        isViolated = YES;
                        break;
                        }
                    if (maxRateConstrained == NO)
                        {
                        maxRateConstrained = YES;
                        maxRate = f;
                        }
                    else if (f < maxRate)
                        maxRate = f;
                    }
                if (y[r->index] >= 0.0 && x[s->index] >= 0.0)
                    {
                    f = (r->nodeDepth - s->nodeDepth) / (y[r->index] - x[s->index]);
                    if (f <= 0.0 || y[r->index] == x[s->index])
                        {
                        if (AreDoublesEqual (r->nodeDepth, s->nodeDepth, tol*0.1) == YES)
                            continue;
                        isViolated = YES;
                        break;
                        }
                    if (minRateConstrained == NO)
                        {
                        minRateConstrained = YES;
                        minRate = f;
                        }
                    else if (f > minRate)
                        minRate = f;
                    }
                }
            }
        if (isViolated == YES)
            break;
        }

    /* check if outright violation */
    if (isViolated == YES)
        {
        MrBayesPrint ("%s   Branch lengths do not satisfy the calibration(s)\n", spacer);
        free (x);
        return (NO);
        }
    
    /* Allow tollerance */
    if (minRateConstrained == YES && maxRateConstrained == YES && AreDoublesEqual (minRate, maxRate, tol) == YES && minRate > maxRate) 
        {
        maxRate = minRate;
        }

    if (minRateConstrained == YES)
        *minClockRate = minRate;
    else
        *minClockRate = 0.0;

    if (maxRateConstrained == YES)
        *maxClockRate = maxRate;
    else
        *maxClockRate = MRBFLT_MAX;

    /* check that minimum and maximum rates are consistent */
    if (minRateConstrained == YES && maxRateConstrained == YES && minRate > maxRate)
        {
        MrBayesPrint ("%s   Branch lengths do not satisfy the calibration(s)\n", spacer);
        free (x);
        return (NO);
        }

    /* date all nodes based on a suitable rate */
    if (minRateConstrained == YES)
        clockRate = minRate;
    else if (maxRateConstrained == YES)
        clockRate = 0.5 * maxRate;
    else
        clockRate = 1.0;
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        p->age = p->nodeDepth / clockRate;
        }

    /* check if there is an age to add (I guess this is here because when max rate is close to minrate and we have numerical precision inacuracy) */
    ageToAdd = 0.0;
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (x[p->index] > 0.0 && x[p->index] > p->age)
            {
            f = x[p->index] - p->age;
            if (f > ageToAdd)
                ageToAdd = f;
            }
        }
    
    /* add extra length if any */
    if (AreDoublesEqual (ageToAdd, 0.0, 0.00000001) == NO)
        {
        for (i=0; i<t->nNodes-1; i++)
            {
            p = t->allDownPass[i];
            p->age += ageToAdd;
            }
        }

    free (x);

    /* reset node depths to ensure that non-dated tips have node depth 0.0 */
    SetNodeDepths(t);

    return (YES);
}


int IsClockSatisfied (Tree *t, MrBFlt tol)
{
    int             i, foundFirstLength, isClockLike;
    MrBFlt          firstLength=0.0, length;
    TreeNode        *p, *q;

    if (t->isRooted == NO)
        return (NO);
        
    foundFirstLength = NO;
    isClockLike = YES;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL && p->right == NULL)
            {
            if (p->isDated == YES)
                {
                //continue;
                length = p->nodeDepth;
                }
            else
                length = 0.0;
            q = p;
            while (q->anc != NULL)
                {
                if (q->anc->anc != NULL)
                    length += q->length;
                q = q->anc;
                }
            if (foundFirstLength == NO)
                {
                firstLength = length;
                foundFirstLength = YES;
                }
            else
                {
                if (AreDoublesEqual (firstLength, length, tol) == NO)
                    {
                    MrBayesPrint ("%s   Node (%s) is not at the same depth as some other tip taking colibration into account. \n", spacer, p->label);
                    isClockLike = NO;
                    }
                }
            }
        }
    if (firstLength < BRLENS_MIN)
        isClockLike = NO;

    return (isClockLike);
}


/* Check that tree obeys topology constraints and that node depths and ages are consistent */
int IsTreeConsistent (Param *param, int chain, int state)
{
    Tree        *tree;
    TreeNode    *p;
    int         i, j;
    MrBFlt      b, r, rAnc, clockRate;
    Param       *subParm;

    if (param->paramType != P_TOPOLOGY && param->paramType != P_BRLENS && param->paramType != P_SPECIESTREE)
        return YES;

    tree      = GetTree(param, chain, state);
    if (modelSettings[param->relParts[0]].clockRate != NULL)
        clockRate = *GetParamVals(modelSettings[param->relParts[0]].clockRate, chain, state);
    else
        clockRate = 1.0;

    if (CheckConstraints(tree)==ERROR) {
        printf ("Tree does not obey constraints\n");
        return NO;
    }

    /* check that the last few indices are not taken in a rooted tree */
    if (tree->isRooted == YES && tree->root->index != tree->nNodes - 1)
        {
        printf ("Problem with root index\n");
        return NO;
        }
    if (tree->isRooted == YES && tree->root->left->index != tree->nNodes - 2)
        {
        printf ("Problem with interior root index\n");
        return NO;
        }

    if (tree->isClock == NO)
        {
        for (i=0; i<tree->nNodes-1; i++)
            {
            p = tree->allDownPass[i];
            if (p->length <= 0.0)
                {
                if (p->length == 0.0)
                    printf ("Node %d has zero branch length %f\n", p->index, p->length);
                else
                    printf ("Node %d has negative branch length %f\n", p->index, p->length);
                return NO;
                }
            }
        return YES;
        }

    /* Clock trees */

    /* Check that lengths and depths are consistent */
    for (i=0; i<tree->nNodes-2; i++) {
        p = tree->allDownPass[i];
        if (p->length < 0.0) {
            printf ("Node %d has negative branch length %f\n", p->index, p->length);
            return NO;
        }
        if (fabs(p->anc->nodeDepth - p->nodeDepth - p->length) > 0.000001) {
            printf ("Node %d has length %f but nodeDepth %f and ancNodeDepth %f\n",
                p->index, p->length, p->nodeDepth, p->anc->nodeDepth);
            return NO;
        }
        if (p->left == NULL && p->isDated == NO && p->nodeDepth != 0.0) {
                printf ("Node %d is an autodated tip (0.0) but has node depth %lf\n",
                    p->index, p->nodeDepth);
                return NO;
        }
    }

    /* Check that ages and calibrations are consistent */
    if (tree->isCalibrated == YES)
        {
        for (i=0; i<tree->nNodes-1; i++)
            {
            p = tree->allDownPass[i];
            if (p->isDated == YES) {
                if (fabs((p->age - p->nodeDepth/clockRate)/p->age) > 0.000001)
                    {
                    printf ("Node %d has age %f but nodeDepth %f when clock rate is %f\n",
                        p->index, p->age, p->nodeDepth, clockRate);
                    return NO;
                    }
                if (p->calibration->prior == fixed && fabs((p->age - p->calibration->priorParams[0])/p->age) > 0.000001)
                    {
                    printf ("Node %d has age %f but should be fixed to age %f\n",
                        p->index, p->age, p->calibration->priorParams[0]);
                    return NO;
                    }
                else if (p->calibration->prior == uniform &&
                        ((p->age - p->calibration->min)/p->age < -0.000001 || (p->age - p->calibration->max)/p->age > 0.000001))
                    {
                    printf ("Node %d has age %f but should be in the interval [%f,%f]\n",
                        p->index, p->age, p->calibration->min, p->calibration->max);
                    return NO;
                    }
                else if ((p->age - p->calibration->min)/p->age < -0.000001)
                    {
                    printf ("Node %d has age %f but should be at least of age %f\n",
                        p->index, p->age, p->calibration->min);
                    return NO;
                    }
                else if ((p->age - p->calibration->max)/p->age > 0.000001)
                    {
                    printf ("Node %d has age %f but should be no older than %f\n",
                        p->index, p->age, p->calibration->max);
                    return NO;
                    }
                }
            }
        }

    for (i=0; i<param->nSubParams; i++)
        {
        subParm = param->subParams[i];
        if (subParm->paramId == TK02BRANCHRATES || (subParm->paramId == MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state) == RCL_TK02))
            {
            rAnc = GetParamVals(subParm, chain, state)[tree->root->left->index];
            if (fabs(rAnc - 1.0) > 1E-6)
                {
                printf ("%s   TK02 relaxed clock mismatch in root rate, which is %e\n", spacer, rAnc);
                return NO;
                }
            for (j=0; j<tree->nNodes-2; j++)
                {
                p = tree->allDownPass[j];
                b = GetParamSubVals(subParm, chain, state)[p->index];
                r = GetParamVals(subParm, chain, state)[p->index];
                rAnc = GetParamVals(subParm, chain, state)[p->anc->index];
                if (fabs(p->length * (r + rAnc) / 2.0 - b) > 0.000001)
                    {
                    printf ("%s   TK02 relaxed clock mismatch in branch %d\n", spacer, p->index);
                    return NO;
                    }
                }
            }
        else if (subParm->paramId == IGRBRANCHRATES || (subParm->paramId == MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state) == RCL_IGR))
            {
            for (j=0; j<tree->nNodes-2; j++)
                {
                p = tree->allDownPass[j];
                b = GetParamSubVals(subParm, chain, state)[p->index];
                r = GetParamVals(subParm, chain, state)[p->index];
                if (fabs(p->length * r - b) > 0.000001)
                    {
                    printf ("%s   Igr relaxed clock mismatch in branch %d\n", spacer, p->index);
                    return NO;
                    }
                }
            }
        }

    if (param->paramType == P_SPECIESTREE)
        return (IsSpeciesTreeConsistent(GetTree(param, chain, state), chain));

    return YES;
}


/* LabelTree: Label tree; remove previous labels if any */
int LabelTree (Tree *t, char **taxonNames)
{
    int         i, nTaxa;
    TreeNode    *p = NULL;

    nTaxa = t->nNodes - t->nIntNodes;
    if (t->isRooted == YES)
        nTaxa--;
    
    /* erase previous labels, if any */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;
        t->nodes[i].label = noLabel;
        }

    /* add labels */
    for (i=0; i<t->nNodes; i++)
        {
        p = &t->nodes[i];
        if (p->left == NULL || (t->isRooted == NO && p->anc == NULL))
            {
            if (p->marked == YES || p->index < 0 || p->index >= nTaxa)
                {
                MrBayesPrint ("%s   Taxon node index repeated or out of range\n", spacer);
                return (ERROR);
                }
            else
                p->label = taxonNames[p->index];
            p->marked = YES;
            }
        else if (p->index > 0 && p->index < nTaxa)
            {
            MrBayesPrint ("%s   Terminal taxon index set for interior node\n", spacer);
            return (ERROR);
            }
        }

    return (NO_ERROR);
}


/*-------------------------------------------------------------------------------------------
|
|   Mark: This routine will mark up a subtree rooted at p
|
---------------------------------------------------------------------------------------------*/
void Mark (TreeNode *p)
{
    if (p != NULL)
        {
        p->marked = YES;
        Mark (p->left);
        Mark (p->right);
        }
}


/*-------------------------------------------------------------------------------------------
 |
 |   MarkDistance: This routine will mark up an unconstrained subtree rooted at p within dist
 |      The distance will be positive in the crown part and negative in the root part.
 |
 ---------------------------------------------------------------------------------------------*/
void MarkDistance (TreeNode *p, int YESorNO, int dist, int *n)
{
    if (p == NULL || p->anc == NULL)
        return;
    
    p->marked = YES;
    if (YESorNO == YES) // in root part
        p->x = p->anc->x -1;
    else               // in crown part
        p->x = p->anc->x +1;
    (*n)++;
        
    if (p->isLocked == NO && abs(p->x) < dist)
        {
        MarkDistance (p->left, YESorNO, dist, n);
        MarkDistance (p->right,YESorNO, dist, n);
        }
}


/*-------------------------------------------------------------------------------------------
 |
 |   MarkUnconstrained: This routine will mark up an unconstrained subtree rooted at p
 |
 ---------------------------------------------------------------------------------------------*/
void MarkUnconstrained (TreeNode *p)
{
    if (p != NULL)
        {
        p->marked = YES;
        if (p->isLocked == NO)
            {
            MarkUnconstrained (p->left);
            MarkUnconstrained (p->right);
            }
        }
}


/*-------------------------------------------------------------------------------------------
|
|   MoveCalculationRoot: This routine will move the calculation root to the terminal with 
|      index outgroup
|
---------------------------------------------------------------------------------------------*/
int MoveCalculationRoot (Tree *t, int outgroup)
{
    int             i;
    TreeNode        *p, *q, *r;
    
    if (t->isRooted == YES || outgroup < 0 || outgroup > t->nNodes - t->nIntNodes - (t->isRooted == YES ? 1 : 0))
        {
        MrBayesPrint ("%s   Problem moving calculation root\n", spacer);
        return (ERROR);
        }

    if (t->root->index == outgroup)
        return (NO_ERROR);    /* nothing to do */

    /* mark the path to the new calculation root */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL && p->right == NULL)
            {
            if (p->index == outgroup)
                p->marked = YES;
            else
                p->marked = NO;
            }
        else
            {
            if (p->left->marked == YES || p->right->marked == YES)
                p->marked = YES;
            else
                p->marked = NO;
            }
        }

    /* rotate the tree to use the specified calculation root */
    p = t->root->left;
    q = t->root;
    q->anc = p;
    q->left = q->right = NULL;
    q->length = p->length;
    while (p->left != NULL && p->right != NULL)
        {
        if (p->left->marked == YES)
            {
            r = p->left;
            p->anc = r;
            p->left = q;
            p->length = r->length;
            q = p;
            p = r;
            }
        else /* if (p->right->marked == YES) */
            {
            r = p->right;
            p->anc = r;
            p->right = q;
            p->length = r->length;
            q = p;
            p = r;
            }
        }
    p->left = p->anc;
    p->right = p->anc = NULL;
    t->root = p;
    p->length = 0.0;

    GetDownPass (t);

    return (NO_ERROR);
}


/*-------------------------------------------------------------------------------------------
|
|   MovePolyCalculationRoot: This routine will move the calculation root to the terminal with 
|      index outgroup and place it as the right-most descendant of the root node
|
---------------------------------------------------------------------------------------------*/
int MovePolyCalculationRoot (PolyTree *t, int outgroup)
{
    int             i;
    PolyNode        *p = NULL, *q, *r;

    /* check if tree is rooted, in which case calculation root is irrelevant */
    if (t->root->left->sib->sib == NULL)
        return (NO_ERROR);

    if (outgroup < 0 || outgroup > t->nNodes - t->nIntNodes)
        {
        MrBayesPrint ("%s   Outgroup index is out of range\n", spacer);
        return (ERROR);
        }

    if (t->root->left->sib->sib->sib != NULL)
        {
        MrBayesPrint ("%s   Root has more than three descendants\n", spacer);
        return (ERROR);
        }

    /* check if rerooting actually necessary */
    if (t->root->left->sib->sib->index == outgroup)
        return (NO_ERROR);
    
    /* mark the path to the new calculation root */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->index == outgroup)
            break;
        }
    if (p->left != NULL)
        {
        MrBayesPrint ("%s   Outgroup index set for internal node\n", spacer);
        for (i=0; i<t->nNodes; i++)
            printf ("%d -- %d\n", i, t->allDownPass[i]->index);
        getchar();
        return (ERROR);
        }

    /* mark path to current root */
    for (i=0; i<t->nNodes; i++)
        t->allDownPass[i]->mark = NO;
    q = p;
    while (q != NULL)
        {
        q->mark = YES;
        q = q->anc;
        }

    /* rotate the tree to use the specified calculation root */
    p = t->root;
    for (;;)
        {
        /* find marked descendant */
        for (q=p->left; q->mark == NO; q=q->sib)
            ;
        if (q->index == outgroup)
            break;
        /* add old root to descendants of that node */
        for (r=q->left; r->sib!=NULL; r=r->sib)
            ;
        r->sib = p;
        p->sib = NULL;   /* should not be needed */
        p->anc = q;
        p->length = q->length;
        /* remove that node from descendants of old root node */
        if (p->left == q)
            p->left = q->sib;
        else
            {
            for (r=p->left; r->sib!=q; r=r->sib)
                ;
            r->sib = r->sib->sib;
            }
        /* make new node root */
        q->sib = NULL;
        q->anc = NULL;
        q->length = 0.0;
        p = q;
        }
    
    /* p is now the new root */
    t->root = p;

    /* finally make sure calculation root is last node among root's descendants */
    for (q=p->left; q->sib!=NULL; q=q->sib)
        ;
    if (q->index != outgroup)
        {
        if (p->left->index == outgroup)
            {
            q->sib = p->left;
            p->left = p->left->sib;
            q->sib->sib = NULL;
            }
        else
            {
            for (r=p->left; r->sib->index!=outgroup; r=r->sib)
                ;
            q->sib = r->sib;
            r->sib = r->sib->sib;
            q->sib->sib = NULL;
            }
        }

    GetPolyDownPass (t);

    return (NO_ERROR);
}


/* 
@return the number of levels for the tree rooted at the "node" 
*/
int NrSubTreeLevels(TreeNode *node)
{
    int r,l;

    if (node == NULL)
        {
        return -1;
        }

    r = NrSubTreeLevels (node->right);
    l = NrSubTreeLevels (node->left);

    return ((r>l)?(r):(l))+1;
}


/*-------------------------------------------------------------------------------------------
|
|   NumConstrainedTips: This routine will return the number of constrained tips, internal or external
|
---------------------------------------------------------------------------------------------*/
int NumConstrainedTips (TreeNode *p)
{
    int     i = 0;

    if (p == NULL)
        return i;
    if (p->left == NULL)
        return 1;

    i += NConstrainedTips (p->left);
    i += NConstrainedTips (p->right);

    return i;
}


/* NConstrainedTips: Recursive function to get the number of constrained tips */
int NConstrainedTips (TreeNode *p)
{
    int     i=0;
    
    if (p!=NULL)
        {
        if (p->left == NULL || p->isLocked == YES)
            return 1;
        else
            {
            i += NConstrainedTips (p->left);
            i += NConstrainedTips (p->right);
            }
        }
    return i;
}


/*-------------------------------------------------------------------------------------------
|
|   NumDatedTips: This routine will return the number of dated tips, internal or external
|
---------------------------------------------------------------------------------------------*/
int NumDatedTips (TreeNode *p)
{
    int     i = 0;

    assert (p != NULL && p->left != NULL);

    i += NDatedTips (p->left);
    i += NDatedTips (p->right);

    return i;
}


/* NDatedTips: recursive function to get the number of dated tips */
int NDatedTips (TreeNode *p)
{
    int     i=0;
    
    assert (p!=NULL);

    if (p->left == NULL || p->isDated == YES)
        return 1;
    else
        {
        i += NDatedTips (p->left);
        i += NDatedTips (p->right);
        return i;
        }
}


/* OrderTips: Order tips in a polytomous tree */
void OrderTips (PolyTree *t)
{
    int         i, j;
    PolyNode    *p, *q, *r, *pl, *ql, *rl;

    /* label by minimum index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            {
            if (t->isRooted == NO && p->index == localOutGroup)
                p->x = -1;
            else
                p->x = p->index;
            }
        else
            {
            j = t->nNodes;
            for (q=p->left; q!=NULL; q=q->sib)
                {
                if (q->x < j)
                    j = q->x;
                }
            p->x = j;
            }
        }

    /* and rearrange */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL || p->anc == NULL)
            continue;
        for (ql=NULL, q=p->left; q->sib!=NULL; ql=q, q=q->sib)
            {
            for (rl=q, r=q->sib; r!=NULL; rl=r, r=r->sib)
                {
                if (r->x < q->x)
                    {
                    if (ql == NULL)
                        p->left = r;
                    if (r == q->sib) /* swap adjacent q and r */
                        {
                        if (ql != NULL)
                            ql->sib = r;
                        pl = r->sib;
                        r->sib = q;
                        q->sib = pl;
                        }
                    else    /* swap separated q and r */
                        {
                        if (ql != NULL)
                            ql->sib = r;
                        pl = r->sib;
                        r->sib = q->sib;
                        rl->sib = q;
                        q->sib = pl;
                        }
                    pl = q;
                    q = r;
                    r = pl;
                    }
                }
            }
        }
    GetPolyDownPass(t);
}


/* PrintNodes: Print a list of tree nodes, pointers and length */
void PrintNodes (Tree *t)
{
    int         i;
    TreeNode    *p;

    printf ("Node\tleft\tright\tanc\tlength\n");
    for (i=0; i<t->nNodes; i++)
        {
        p = &t->nodes[i];
        printf ("%d\t%d\t%d\t%d\t%f\t%f\n",
            p->index,
            p->left == NULL ? -1 : p->left->index,
            p->right == NULL ? -1 : p->right->index,
            p->anc == NULL ? -1 : p->anc->index,
            p->length,
            p->nodeDepth);
        }

    if (t->root == NULL)
        printf ("root: NULL\n");
    else
        printf ("root: %d\n", t->root->index);

    printf ("allDownPass:");
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p!=NULL)
            printf ("  %d", p->index);
        else
            printf ("  NULL");
        }
    printf ("\nintDownPass:  ");
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        if (p!=NULL)
            printf ("  %d\t", p->index);
        else
            printf ("  NULL\t");
        }
    printf ("\n");
}


/* PrintPolyNodes: Print a list of polytomous tree nodes, pointers and length */
void PrintPolyNodes (PolyTree *pt)
{
    int         i, j, k;
    PolyNode    *p;

    printf ("Node\tleft\tsib\tanc\tlength\tlabel\n");
    for (i=0; i<pt->memNodes; i++)
        {
        p = &pt->nodes[i];
        printf ("%d\t%d\t%d\t%d\t%f\t%s\n",
            p->index,
            p->left == NULL ? -1 : p->left->index,
            p->sib == NULL ? -1 : p->sib->index,
            p->anc == NULL ? -1 : p->anc->index,
            p->length,
            p->label);
        }
    printf ("root: %d\n", pt->root->index);
    fflush(stdout);

    if (pt->nBSets > 0)
        {
        for (i=0; i<pt->nBSets; i++)
            {
            printf ("Effective branch length set '%s'\n", pt->bSetName[i]);
            for (j=0; j<pt->nNodes; j++)
                {
                printf ("%d:%f", j, pt->effectiveBrLen[pt->nBSets][j]);
                if (j != pt->nNodes-1)
                    printf (", ");
                }
            printf ("\n");
            }
        }
    
    if (pt->nESets > 0)
        {
        for (i=0; i<pt->nESets; i++)
            {
            printf ("Cpp event set '%s'\n", pt->eSetName[i]);
            for (j=0; j<pt->nNodes; j++)
                {
                if (pt->nEvents[i*pt->nNodes+j] > 0)
                    {
                    printf ("\tNode %d -- %d:(", j, pt->nEvents[i][j]);
                    for (k=0; k<pt->nEvents[i][j]; k++)
                        {
                        printf ("%f %f", pt->position[i][j][k], pt->rateMult[i][j][k]);
                        if (k != pt->nEvents[i][j]-1)
                            printf (", ");
                        }
                    printf (")\n");
                    }
                }
            printf ("\n");
            }
        }

    fflush(stdout);
}


/* PrintTranslateBlock: Print a translate block to file fp for tree t */
void PrintTranslateBlock (FILE *fp, Tree *t)
{
    int     i, j, nTaxa;

    if (t->isRooted == NO)
        nTaxa = t->nNodes - t->nIntNodes;
    else
        nTaxa = t->nNodes - t->nIntNodes - 1;

    fprintf (fp, "\ttranslate\n");

    for (i=0; i<nTaxa; i++)
        {
        for (j=0; j<t->nNodes; j++)
            if (t->allDownPass[j]->index == i)
                break;
        if (i == nTaxa-1)
            fprintf (fp, "\t\t%d\t%s;\n", i+1, t->allDownPass[j]->label);
        else
            fprintf (fp, "\t\t%d\t%s,\n", i+1, t->allDownPass[j]->label);
        }
}


/**
Update relaxed clock parameter of the branch of a node with index "b" after node with index "a" is removed. 
i.e. make branch of node with index "b" be a concatenation of its original branch and the branch of node with index "a"
Relaxed clock parameter of node with index "a" become invalid in the process.
Note: For Non-clock models the routine has no effect.

|       |
|       |
a       |
|   ->  |
|       |
|       b
b                */
void AppendRelaxedBranch (int a,int b,PolyTree *t)
{
    int i,len;

    for (i=0; i<t->nBSets; i++)
        {
        t->effectiveBrLen[i][b] += t->effectiveBrLen[i][a];
        }

    for (i=0; i<t->nESets; i++)
        {
        len=t->nEvents[i][a]+t->nEvents[i][b];
        t->position[i][a] = (MrBFlt *) SafeRealloc ((void *)t->position[i][a], len*sizeof(MrBFlt));
        t->rateMult[i][a] = (MrBFlt *) SafeRealloc ((void *)t->rateMult[i][a], len*sizeof(MrBFlt));
        memcpy (t->position[i][a]+t->nEvents[i][a], t->position[i][b], t->nEvents[i][b]*sizeof(MrBFlt));
        memcpy (t->rateMult[i][a]+t->nEvents[i][a], t->rateMult[i][b], t->nEvents[i][b]*sizeof(MrBFlt));
        free(t->position[i][b]);
        free(t->rateMult[i][b]);
        t->position[i][b] = t->position[i][a];
        t->rateMult[i][b] = t->rateMult[i][a];
        t->position[i][a] = NULL;
        t->rateMult[i][a] = NULL;
        t->nEvents[i][a] = 0;
        t->nEvents[i][b] = len;
        }

}


/**
Swap relaxed clock paramiters of the branch of nodes with index "a" and "b".
*/
void SwapRelaxedBranchInfo (int a,int b,PolyTree *t)
{
    int i,j;
    MrBFlt tmp, *tmpp;

    for (i=0; i<t->nBSets; i++)
        {
        tmp = t->effectiveBrLen[i][a];
        t->effectiveBrLen[i][a] = t->effectiveBrLen[i][b];
        t->effectiveBrLen[i][b] = tmp;
        }
    if (t->popSizeSet == YES)
        {
        tmp = t->popSize[a];
        t->popSize[a]=t->popSize[b];
        t->popSize[b] = tmp;
        }

    for (i=0; i<t->nESets; i++)
        {
        tmpp = t->position[i][a];
        t->position[i][a] = t->position[i][b];
        t->position[i][b] = tmpp;
        tmpp = t->rateMult[i][a];
        t->rateMult[i][a] = t->rateMult[i][b];
        t->rateMult[i][b] = tmpp;
        j = t->nEvents[i][a];
        t->nEvents[i][a] = t->nEvents[i][b];
        t->nEvents[i][b] = j;
        }
}


/*-------------------------------------------------------------------------------------------
|
|   PrunePolyTree: This routine will prune a polytomous tree according to the currently
|      included taxa.  NB! All tree nodes cannot be accessed by cycling over the
|      pt->nodes array after the deletion, because some spaces will be occupied by deleted
|      nodes and pt->nNodes is no longer the length of this array
|      (if proper re-arangment of pt->nodes needed then remove "#if 0" and make changes to p->x, see below).
|
---------------------------------------------------------------------------------------------*/
int PrunePolyTree (PolyTree *pt)
{
    int             i, j, numDeleted, numTermPruned, numIntPruned, index;
    PolyNode        *p = NULL, *q=NULL, *r=NULL, *qa;

    numDeleted = 0;
    for (i=0; i<pt->nNodes; i++)
        {
        p = pt->allDownPass[i];
        CheckString (taxaNames, numTaxa, p->label, &index);
        if (p->left == NULL && taxaInfo[index].isDeleted == YES)
            numDeleted++;
        }
        
    if (numDeleted == 0)
        {
        /* nothing to do */
        return (NO_ERROR);
        }
    if (pt->nNodes - pt->nIntNodes - numDeleted < 3)
        {
        MrBayesPrint ("%s   Pruned tree has less than three taxa in it\n", spacer);
        return (ERROR);
        }
    if (pt->nNodes - pt->nIntNodes < numLocalTaxa)
        {
        MrBayesPrint ("%s   Tree to be pruned does not include all taxa\n", spacer);
        return (ERROR);
        }

    /* prune away one node at a time */
    numIntPruned = 0;
    numTermPruned = 0;
    for (i=0; i<pt->nNodes; i++)
        {
        p = pt->allDownPass[i];
        if (p->left != NULL)
            continue;
        CheckString (taxaNames, numTaxa, p->label, &index);
        if (taxaInfo[index].isDeleted == YES)
            {
            numTermPruned++;
            for (q=p->anc->left; q!=NULL; q=q->sib)
                {
                if (q->sib == p)
                    break;
                }
            if (q == NULL)
                {
                /* p is the left of its ancestor */
                assert (p->anc->left == p);
                p->anc->left = p->sib;
                }
            else
                {
                /* p is q->sib; this also works if p->sib is NULL */
                q->sib = p->sib;
                }
            /* if only one child left, delete ancestral node */
            j = 0;
            for (q=p->anc->left; q!=NULL; q=q->sib)
                j++;
            if (j == 1)
                {
                /* p->anc->left is only child left; make p->anc be p->anc->left and accommodate its length */
                numIntPruned++;
                qa= p->anc;
                q = qa->left;
                if (q->left == NULL)
                    {
                    AppendRelaxedBranch (qa->index, q->index, pt);
                    qa->index = q->index;
                    qa->length += q->length;
                    strcpy(qa->label, q->label);
                    qa->left = NULL;
                    /* To make sure that q is not treated as the representer of the tip it represented before. i.e. make  condition if (p->left != NULL) true */
                    q->left = (struct pNode*)1; 
                    }
                else
                    {
                    if (qa->anc != NULL)
                        {
                        AppendRelaxedBranch (qa->index, q->index, pt);
                        qa->length += q->length;
                        }
                    qa->index   = q->index;
                    qa->left = q->left;
                    for (r=q->left; r!= NULL; r=r->sib)
                        r->anc = qa;
                    }
                }
            /* if unrooted, then root node has to have more then 2 children, thus the following check */
            if (j == 2 && pt->isRooted == NO && p->anc->anc == NULL)
                {
                numIntPruned++;
                r=p->anc; /*r is the root with only 2 children*/
                if (r->left->left != NULL)
                    {/* Make r->left new root by attaching "right" child of r to children of r->left */
                    for (q=r->left->left; q->sib!=NULL; q=q->sib)
                        ;
                    q->sib = r->left->sib;
                    r->left->sib->anc = q->anc;
                    r->left->sib->length += q->anc->length;
                    r->left->sib = NULL;
                    r->left->anc = NULL;
                    pt->root = r->left;
                    }
                else
                    {/* Make "right" child of r (r->left->sib) the new root by attaching r->left to children of r->"right" */
                    for (q=r->left->sib->left; q->sib!=NULL; q=q->sib)
                        ;
                    q->sib = r->left;
                    r->left->anc = q->anc;
                    r->left->length += q->anc->length;
                    r->left->sib = NULL;
                    q->anc->anc = NULL;
                    pt->root = q->anc;
                    }
                }
            }
        }

#if 0 
    /* place unused space at the end of pt->nodes array. If activated this code p->x has to be set to non 0 value for all p that are deleted. */
    for (i=0; i<pt->nNodes; i++)
        {
        p = &pt->nodes[i];
        if (p->x != 0)
            {
            for (j=i+1; j<pt->nNodes; j++)
                {
                q = &pt->nodes[j];
                if (q->x == 0)
                    break;
                }
            if (j != pt->nNodes)
                {
                /* swap nodes; quite difficult! */
                CopyPolyNodes (p, q, nLongsNeeded);
                p->left = q->left;
                p->sib = q->sib;
                p->anc = q->anc;
                for (k=0; k<pt->nNodes; k++)
                    {
                    r = &pt->nodes[k];
                    if (r->left == q)
                        r->left = p;
                    if (r->sib == q)
                        r->sib = p;
                    if (r->anc == q)
                        r->anc = p;
                    }
                }
            }
        }
#endif

    /* correct number of nodes */
    pt->nNodes -= (numTermPruned + numIntPruned);
    pt->nIntNodes -= numIntPruned;
    
    /* get downpass; note that the deletion procedure does not change the root in rooted case */
    i=j=0;
    GetPolyNodeDownPass (pt, pt->root, &i, &j);
    assert (i==pt->nNodes);
    assert (j==pt->nIntNodes);

    return (NO_ERROR);
}


/*--------------------------------------------------------------------
|
|       RandPerturb: Randomly perturb a tree by nPert NNIs
|
---------------------------------------------------------------------*/
int RandPerturb (Tree *t, int nPert, RandLong *seed)
{
    int         i, whichNode;
    TreeNode    *p, *q, *a, *b, *c;
    
    if (t->nConstraints >= t->nIntNodes)
        {
        MrBayesPrint ("%s   User tree cannot be perturbed because all nodes are locked\n", spacer);
        return (ERROR);
        }

    for (i=0; i<nPert; i++)
        {
        do
            {
            whichNode = (int)(RandomNumber(seed) * (t->nIntNodes - 1));
            p = t->intDownPass[whichNode];
            } while (p->isLocked == YES);
        
        q = p->anc;
        a  = p->left;
        b  = p->right;
        if (q->left == p)
            c  = q->right;
        else    
            c  = q->left;
        
        if (RandomNumber(seed) < 0.5)
            {
            /* swap b and c */
            p->right = c;
            c->anc  = p;

            if (q->left == c)
                q->left = b;
            else
                q->right = b;
            b->anc = q;
            }
        else
            {
            /* swap a and c */
            p->left = c;
            c->anc  = p;

            if (q->left == c)
                q->left = a;
            else
                q->right = a;
            a->anc = q;
            }

        if (t->isCalibrated == YES)
            InitCalibratedBrlens (t, 0.0001, seed);
        else if (t->isClock == YES)
            InitClockBrlens (t);
        }
    
    GetDownPass (t);

    if (t->checkConstraints == YES && CheckConstraints (t) == NO_ERROR)
        {
        MrBayesPrint ("%s   Broke constraints when perturbing tree\n", spacer);
        return (ERROR);
        }

    return (NO_ERROR);
}


/*
|       Reorder array of nodes "nodeArray" such that first nodes in it could be paired with "w" to create imediat common ancestor and this ancesor node would not vialate any constraint.  
|       
| @param w                      Reference node as described 
| @param nodeArray              A set of node to order as described
| @param activeConstraints      Array containing indeces of active constraints in the set of defined constraints
| @param nLongsNeeded           Length of partition information (in BitsLong) in a node and constraint deffinition.
| @param isRooted               Do constraints apply to rootet tree YES or NO
|
| @return                       Number of nodes in "nodeArray" that could be paired  with "w" to create imediat common ancestor and this ancesor node would not vialate any constraint
*/
int ConstraintAllowedSet(PolyNode *w, PolyNode **nodeArray, int nodeArraySize, int *activeConstraints, int activeConstraintsSize, int nLongsNeeded, int isRooted)
{
    int             i, j,  k, FirstEmpty;
    BitsLong        **constraintPartition;
    PolyNode        *tmp;

    for (j=0; j<activeConstraintsSize; j++)
        {
        k=activeConstraints[j];

        if (definedConstraintsType[k] == PARTIAL)
            {
            if ((IsPartNested(definedConstraintPruned[k], w->partition, nLongsNeeded) == YES) ||
                (isRooted == NO && IsPartNested(definedConstraintTwoPruned[k], w->partition, nLongsNeeded) == YES))
                continue;/* all nodes are compartable because condition of the constraint has to be sutsfied in the subtree rooted at w*/

            FirstEmpty = IsSectionEmpty(definedConstraintPruned[k], w->partition, nLongsNeeded);
            if (FirstEmpty == YES &&  IsSectionEmpty(definedConstraintTwoPruned[k], w->partition, nLongsNeeded) == YES)
                continue; /* all nodes are compartable becouse w does not contain any constraint taxa*/

            assert (FirstEmpty^IsSectionEmpty(definedConstraintTwoPruned[k], w->partition, nLongsNeeded));

            if (FirstEmpty == YES)
                {/*w->partition has intersection with definedConstraintTwoPruned[k], thus remove all nodes from nodeArray that intersect with definedConstraintPruned[k]*/
                constraintPartition=definedConstraintPruned;
                }
            else
                {/*w->partition has intersection with definedConstraintPruned[k], thus remove all nodes from nodeArray that intersect with definedConstraintTwoPruned[k]*/
                constraintPartition=definedConstraintTwoPruned;
                }

            for (i=0;i<nodeArraySize;i++)
                {
                if (IsSectionEmpty(constraintPartition[k], nodeArray[i]->partition, nLongsNeeded) == NO &&
                    ((FirstEmpty == NO && isRooted== YES) ||  IsPartNested(constraintPartition[k], nodeArray[i]->partition, nLongsNeeded) == NO))
                  /*second part of if statment is to bail out "nodeArray[i]" when "w" contains nodes for example from definedConstraintPruned and "nodeArray[i]" have definedConstraintTwoPruned fully nested in it
                  This bail out not applicable if t->isRooted== YES Since we should create a rooting node for the first set of taxa in the constraint.
                  Note that such case possible because we may have hard constraint node that fully nest definedConstraintTwoPruned but also having taxa from definedConstraintPruned keeping constraint active.*/
                    {
                    tmp = nodeArray[i];
                    nodeArray[i]=nodeArray[--nodeArraySize];
                    nodeArray[nodeArraySize]=tmp;
                    i--;
                    }
                }
            }/*end if PARTIAL*/
        else 
            {
            assert (definedConstraintsType[k] == NEGATIVE);
            if (isRooted == YES || IsBitSet(localOutGroup, definedConstraintPruned[k])==NO)
                constraintPartition=definedConstraintPruned;
            else
                constraintPartition=definedConstraintTwoPruned;
            
            if (IsSectionEmpty(constraintPartition[k], w->partition, nLongsNeeded)==YES)
                continue;

            for (i=0;i<nodeArraySize;i++)
                {
                if (IsUnionEqThird (w->partition, nodeArray[i]->partition, constraintPartition[k], nLongsNeeded) == YES)
                    {
                    tmp = nodeArray[i];
                    nodeArray[i]=nodeArray[--nodeArraySize];
                    nodeArray[nodeArraySize]=tmp;
                    i--;
                    }
                }

            }/*end if NEGATIVE*/
        }

   return nodeArraySize;
}


/*
|               Check if "partition" violate any constraint.  
|       
| @param partiton               Partition to check 
| @param activeConstraints      Array containing indeces of active constraints in the set of defined constraints  
| @param nLongsNeeded           Length of partition information (in BitsLong) in a node and constraint deffinition
| @param isRooted               Do constraints apply to rootet tree YES or NO
|
| @return                       Index of first violated constraint in activeConstraints array, -1 if no constraint is violated.
*/
int ViolatedConstraint(BitsLong *partition, int *activeConstraints, int activeConstraintsSize, int nLongsNeeded, int isRooted)
{
    int             j, k;
    BitsLong        **constraintPartition;

    for (j=0; j<activeConstraintsSize; j++)
        {
        k=activeConstraints[j];
        assert (definedConstraintsType[k] != HARD);

        if (definedConstraintsType[k] == PARTIAL)
            {
            if ((IsSectionEmpty(definedConstraintPruned[k], partition, nLongsNeeded) == NO) &&
                (IsSectionEmpty(definedConstraintTwoPruned[k], partition, nLongsNeeded) == NO) &&
                (IsPartNested(definedConstraintPruned[k], partition, nLongsNeeded) == NO) &&
                !(isRooted == NO && IsPartNested(definedConstraintTwoPruned[k], partition, nLongsNeeded) == YES))
                return j;
            }/*end if PARTIAL*/
        else 
            {
            assert (definedConstraintsType[k] == NEGATIVE);
            if (isRooted == YES || IsBitSet(localOutGroup, definedConstraintPruned[k])==NO)
                constraintPartition=definedConstraintPruned;
            else
                constraintPartition=definedConstraintTwoPruned;

            if (IsUnionEqThird (partition, partition, constraintPartition[k], nLongsNeeded) == YES)
                return j;
            }/*end if NEGATIVE*/
        }

   return -1;
}


/*
|         Remove from activeConstraints references to constraints that become satisfied if PolyNode "w" exist, i.e. they do not need to be checked furter thus become not active  
|
| @param activeConstraints      Array containing indeces of active constraints in the set of defined constraints
| @param nLongsNeeded           Length of partition information (in BitsLong) in a node and constraint deffinition.
| @param isRooted               Do constraints apply to rootet tree YES or NO
|
| @return                       Size of pruned "activeConstraints" array
*/
int PruneActiveConstraints (PolyNode *w, int *activeConstraints, int activeConstraintsSize, int nLongsNeeded, int isRooted)
{
    int             j,  k;
    BitsLong        **constraintPartition;
    //PolyNode        *tmp;

    for (j=0; j<activeConstraintsSize; j++)
        {
        k=activeConstraints[j];

        if (definedConstraintsType[k] == PARTIAL)
            {
            if ((IsPartNested(definedConstraintPruned[k], w->partition, nLongsNeeded) == YES && IsSectionEmpty(definedConstraintTwoPruned[k], w->partition, nLongsNeeded)) ||
               (isRooted == NO && IsPartNested(definedConstraintTwoPruned[k], w->partition, nLongsNeeded) == YES && IsSectionEmpty(definedConstraintPruned[k], w->partition, nLongsNeeded)))
                {
                //tmp = activeConstraints[j];
                activeConstraints[j]=activeConstraints[--activeConstraintsSize];
                //activeConstraints[activeConstraintsSize]=tmp;
                j--;
                }
            }/*end if PARTIAL*/
        else 
            {
            assert (definedConstraintsType[k] == NEGATIVE);
            if (isRooted == YES || IsBitSet(localOutGroup, definedConstraintPruned[k])==NO)
                constraintPartition=definedConstraintPruned;
            else
                constraintPartition=definedConstraintTwoPruned;
            
            if (IsPartNested(constraintPartition[k], w->partition, nLongsNeeded)==NO && IsSectionEmpty(constraintPartition[k], w->partition, nLongsNeeded)==NO)
                {
                //tmp = activeConstraints[j];
                activeConstraints[j]=activeConstraints[--activeConstraintsSize];
                //activeConstraints[activeConstraintsSize]=tmp;
                j--;
                }
            }/*end if NEGATIVE*/
        }

   return activeConstraintsSize;
}


/*--------------------------------------------------------------------
|
|           RandResolve: Randomly resolve a polytomous tree
|
| @param    tt is a tree which contains information about applicable constraints. If it is set to NULL then no constraints will be used. 
|           If t!=NULL then partitions of nodes of polytree should be allocated for example by AllocatePolyTreePartitions (t);
| @return   NO_ERROR on succes, ABORT if could not resolve a tree without vialating some consraint, ERROR if any other error occur 
---------------------------------------------------------------------*/
int RandResolve (Tree *tt, PolyTree *t, RandLong *seed, int destinationIsRooted)
{
    int         i, j, k, nextNode, stopNode, rand1, rand2, nTaxa, nLongsNeeded, tmp;
    PolyNode    *p=NULL, *q, *r, *u, *w1, *w2;
    int         nodeArrayAllowedSize, nodeArraySize, activeConstraintsSize;
    PolyNode    **nodeArray;
    int         *activeConstraints;

    assert (tt==NULL || t->bitsets!=NULL); /* partition fields of t nodes need to be allocated if constraints are used*/
    nTaxa = t->nNodes - t->nIntNodes;     /* different from numLocalTaxa potentially if a species tree */
    assert (nTaxa <= t->memNodes/2); /* allocated tree has to be big enough*/
    nLongsNeeded = (nTaxa - 1) / nBitsInALong + 1; /* allocated lenght of partitions is t->memNodes/2 bits but only first nTaxa bits are used */

    nodeArray = t->allDownPass; /*temporary use t->allDownPass for different purpose. It get properly reset at the end. */
    activeConstraints = tempActiveConstraints;
    activeConstraintsSize = 0;

    /* collect constraints to consider if applicable*/
    if (tt!=NULL && tt->constraints!=NULL)
        {
        for (k=0; k<numDefinedConstraints; k++)
            {
            if (tt->constraints[k] == YES && definedConstraintsType[k] != HARD)
                activeConstraints[activeConstraintsSize++]=k;
            }
        }

    /* count immediate descendants */
    GetPolyDownPass(t);
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        tmp=ViolatedConstraint(p->partition, activeConstraints, activeConstraintsSize, nLongsNeeded, t->isRooted);
        if (tmp != -1)
            {
            assert (p->isLocked == YES);
            MrBayesPrint ("%s   Could not build a constraint tree since hard constraint \"%s\" and constraint \"%s\" are incompatible\n", spacer, constraintNames[p->lockID], constraintNames[activeConstraints[tmp]]);
            return (ERROR);
            }
        activeConstraintsSize = PruneActiveConstraints (p, activeConstraints, activeConstraintsSize, nLongsNeeded, t->isRooted);
        j = 0;
        for (q=p->left; q!=NULL; q=q->sib)
            j++;
        p->x = j;
        }

    /* add one node at a time */
    if (destinationIsRooted == NO)
        stopNode = 2*nTaxa - 2;
    else
        stopNode = 2*nTaxa - 1;
    for (nextNode=t->nNodes; nextNode < stopNode; nextNode++)
        {
        /* find a polytomy to break */
        for (i=0; i<t->nIntNodes; i++)
            {
            p = t->intDownPass[i];
            if (destinationIsRooted == YES && p->x > 2)
                break;
            if (destinationIsRooted == NO && ((p->anc != NULL && p->x > 2) || (p->anc == NULL && p->x > 3)))
                break;
            }

        /* if we can't find one, there's an error */
        if (i == t->nIntNodes)
            {
            return  ERROR;
            }

        nodeArraySize=0;
        /*Collect initial list of candidate nodes to join*/
        for (q = p->left; q!= NULL; q = q->sib)
            {
            nodeArray[nodeArraySize++]=q;
            }
        assert (nodeArraySize==p->x);

        /* identify two descendants randomly */
        /* make sure we do not select outgroup if it is an unrooted tree */
        if (p->anc == NULL && destinationIsRooted == NO)
            nodeArraySize--;

        do
            {
            /* Pick first node */
            rand1 = (int) (RandomNumber(seed) * nodeArraySize);
            w1 = nodeArray[rand1];
            nodeArray[rand1] = nodeArray[--nodeArraySize];

            if (nodeArraySize==0)
                return ABORT; /* Potentaily here we could instead revert by removing last added node and try again. */

            /* Move all nodes in nodeArray which can be paired with w to the begining of array */
            nodeArrayAllowedSize=ConstraintAllowedSet(w1, nodeArray, nodeArraySize, activeConstraints, activeConstraintsSize, nLongsNeeded, t->isRooted);
            /* TODO optimization for Maxim (if not Maxim remove it if you still see it): if nodeArrayAllowedSize==0 then set w1->y */
            } while (nodeArrayAllowedSize == 0);

        rand2 = (int) (RandomNumber(seed) *nodeArrayAllowedSize);
        w2 = nodeArray[rand2];

        /* create a new node */
        u = &t->nodes[nextNode];
        u->anc = p;
        u->x = 2;
        p->x--;

        if (tt != NULL) {
            for (j=0; j<nLongsNeeded; j++)
                u->partition[j] = w1->partition[j] | w2->partition[j] ;
            activeConstraintsSize = PruneActiveConstraints (u, activeConstraints, activeConstraintsSize, nLongsNeeded, t->isRooted);
        }

        u->left = w1;
        t->nNodes++;
        t->nIntNodes++;

        /* connect tree together */
        r = u;
        for (q = p->left; q!= NULL; q = q->sib)
            {
            if (q != w1 && q != w2)
                {
                r->sib=q;
                r = q;
                }
            }
        r->sib = NULL;
        w1->sib = w2;
        w2->sib = NULL;
        w1->anc = u;
        w2->anc = u;
        p->left = u;

        /* update tree */
        GetPolyDownPass (t);
        }

    /* relabel interior nodes (important that last indices are at the bottom!) */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        p->index = nTaxa + i;
        }
    return NO_ERROR;
}


/* ResetTreeNode: Reset tree node except for memory index */
void ResetTreeNode (TreeNode *p)
{
    /* do not change memoryIndex; that is set once and for all when tree is allocated */
    p->index                  = 0; 
    p->upDateCl               = NO;
    p->upDateTi               = NO;
    p->marked                 = NO;
    p->length                 = 0.0;
    p->nodeDepth              = 0.0;
    p->x                      = 0;
    p->y                      = 0;
    p->index                  = 0;
    p->isDated                = NO;
    p->calibration            = NULL;
    p->age                    = -1.0;
    p->isLocked               = NO;
    p->lockID                 = -1;
    p->label                  = noLabel;
    p->d                      = 0.0;
    p->partition              = NULL;
}


/* ResetPolyNode: Reset all values of one node in a polytree */
void ResetPolyNode (PolyNode *p)
{
    /* we reset everything here except memoryIndex, which should be immutable */
    p->length = 0.0;
    p->depth = 0.0;
    p->age = 0.0;
    p->anc = p->left = p->sib = NULL;
    p->calibration = NULL;
    p->f = 0.0;
    p->index = 0;
    p->isDated = NO;
    p->isLocked = NO;
    strcpy (p->label,"");
    p->lockID = 0;
    p->partition = NULL;
    p->partitionIndex = 0;
    p->support = 0.0;
    p->x = p->y = 0;
}


/* ResetPolyTree: Reset polytomous tree to pristine state but keep relevant memory. */
void ResetPolyTree (PolyTree *pt)
{
    int     i, maxTaxa, nLongsNeeded;

    /* clear nodes */
    for (i=0; i<pt->memNodes; i++)
        ResetPolyNode (&pt->nodes[i]);

    /* empty node arrays and tree properties but keep space */
    for (i=0; i<pt->nNodes; i++)
        pt->allDownPass[i] = NULL;
    for (i=0; i<pt->nIntNodes; i++)
        pt->intDownPass[i] = NULL;
    pt->nNodes = 0;
    pt->nIntNodes = 0;
    pt->root = NULL;
    pt->brlensDef = NO;
    pt->isRooted = NO;
    pt->isClock = NO;
    pt->isRelaxed = NO;
    pt->clockRate = 0.0;

    /* empty bitsets but keep space and pointers */
    if (pt->bitsets)
        {
        maxTaxa = pt->memNodes / 2;
        nLongsNeeded = (maxTaxa - 1) / nBitsInALong + 1;
        for (i=0; i<pt->memNodes*nLongsNeeded; i++)
            pt->bitsets[i] = 0;
        for (i=0; i<pt->memNodes; i++)
            pt->nodes[i].partition = pt->bitsets + i*nLongsNeeded;
        }

    /* empty relaxed clock parameters */
    FreePolyTreeRelClockParams (pt);

    /* empty population size set parameters */
    FreePolyTreePopSizeParams (pt);
}


/* ResetPolyTreePartitions: Reset and set bit patterns describing partitions */
void ResetPolyTreePartitions (PolyTree *pt)
{
    int         i, j, numTaxa, nLongsNeeded;
    PolyNode    *pp;

    /* get some handy numbers */
    numTaxa = pt->memNodes/2;
    nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;
    
    /* reset bits describing partitions */
    for (i=0; i<pt->memNodes*nLongsNeeded; i++)
        {
        pt->bitsets[i] = 0;
        }

    /* set bits describing partitions */
    for (i=0; i<pt->nNodes; i++)
        {
        assert (pt->allDownPass != NULL && pt->allDownPass[i] != NULL);
        assert (pt->allDownPass[i]->partition != NULL);
        
        pp = pt->allDownPass[i];
        if (pp->left == NULL)
            {
            SetBit (pp->index, pp->partition);
            }
        if (pp->anc != NULL)
            {
            for (j=0; j<nLongsNeeded; j++)
                pp->anc->partition[j] |= pp->partition[j];
            }
        }
}


/*----------------------------------------------
|
|   ResetRootHeight: Reset node heights in a clock
|      tree to fit a new root height. Assumes
|      node depths and lengths set correctly.
|
-----------------------------------------------*/
int ResetRootHeight (Tree *t, MrBFlt rootHeight)
{
    int         i;
    TreeNode    *p;
    MrBFlt      factor, x, y;

    if (t->isClock == NO)
        return ERROR;
    
    /* make sure node depths are set */
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->nodeDepth = 0.0;
        else
            {
            x = p->left->nodeDepth + p->left->length;
            y = p->right->nodeDepth + p->right->length;
            if (x > y)
                p->nodeDepth = x;
            else
                p->nodeDepth = y;
            }
        }
    for (i=t->nNodes-3; i>=0; i--)
        {
        p = t->allDownPass[i];
        p->nodeDepth = p->anc->nodeDepth - p->length;
        }

    /* now reset node depths and branch lengths */
    factor = rootHeight / t->root->left->nodeDepth;
    t->root->left->nodeDepth = rootHeight;
    for (i=t->nNodes-2; i>=0; i--)
        {
        p = t->allDownPass[i];
        p->nodeDepth *= factor;
        p->length *= factor;
        }

    return NO_ERROR;
}


/*----------------------------------------------
|
|   ResetTipIndices: reset tip indices to be from 
|      0 to number of included taxa, in same order
|      as in the original taxon set.
|
-----------------------------------------------*/
void ResetTipIndices(PolyTree *pt)
{
    int         i, j, k, m;
    PolyNode    *p = NULL;

    for (i=j=0; i<numTaxa; i++)
        {
        for (k=0; k<pt->nNodes; k++)
            {
            p = pt->allDownPass[k];
            if (StrCmpCaseInsensitive(p->label,taxaNames[i]) == 0)
                break;
            }
        if (k < pt->nNodes)
            {
            assert (p->left == NULL);
            if (p->index!=j) {
                SwapRelaxedBranchInfo (p->index, j, pt);
                for (m=0; m<pt->nNodes; m++)
                    {
                    if (pt->allDownPass[m]->index==j)
                        {
                        pt->allDownPass[m]->index=p->index;
                        break;
                        }
                    }
                p->index = j;
                }
            j++;
            }
        }
}


/*----------------------------------------------
|
|   ResetTopology: rebuild the tree t to fit the 
|      Newick string s. Everyting except topology
|      is left in the same state in t.
|
-----------------------------------------------*/
int ResetTopology (Tree *t, char *s)
{
    TreeNode    *p, *q;
    int         i, j, k, inLength;
    char        temp[30];
    
    /* set all pointers to NULL */
    for (i=0; i<t->memNodes; i++)
        {
        p = &t->nodes[i];
        p->anc = p->right = p->left = NULL;
        p->index = -1;
        }
    p = &t->nodes[0];

    /* start out assuming that the tree is rooted; we will detect below if it is not */
    t->isRooted = YES;
    inLength = NO;
    for (i=0, j=1; *s!='\0'; s++)
        {
        if (*s == ',' || *s == ')' || *s == ':')
            {
            if (p->right == NULL && inLength == NO)
                {
                temp[i] = '\0';
                k = atoi (temp);
                p->index = k-1;
                i = 0;
                }
            else
                inLength = NO;
            }
        if (*s == '(')
            {
            q = p;
            p = &t->nodes[j++];
            q->left = p;
            p->anc = q;
            }
        else if (*s == ',')
            {
            if (p->anc->right == NULL)
                {
                q = p->anc;
                p = &t->nodes[j++];
                p->anc = q;
                q->right = p;
                }
            else /* if p->anc->right == p (near 'root' of unrooted trees) */
                {
                q = p->anc;
                p = &t->nodes[j++];
                q->anc = p;
                p->left = q;
                t->isRooted = NO;
                }
            }
        else if (*s == ')')
            {
            p = p->anc;
            }
        else if (*s == ':')
            {
            inLength = YES;
            }
        else if (inLength == NO)
            {
            temp[i++] = *s;
            }
        }

    /* attach root to rooted tree */
    if (t->isRooted == YES)
        {
        p = &t->nodes[0];
        q = &t->nodes[j++];
        q->left = p;
        p->anc = q;
        }

    /* relabel interior nodes, find number of nodes and root */
    t->nNodes = j;
    t->nIntNodes = t->nNodes/2 - 1;

    if (t->isRooted == NO)
        j = t->nNodes - t->nIntNodes;
    else
        j = t->nNodes - t->nIntNodes - 1;

    for (i=0; i<t->nNodes; i++)
        {
        p = &t->nodes[i];
        if (p->index == -1)
            p->index = j++;
        if (p->anc == NULL)
            t->root = p;
        }

    GetDownPass (t);

    return NO_ERROR;
}


/*-----------------------------------------------------------------
|
|   ResetBrlensFromTree: copies brlens and depths from second tree (vTree) to
|       first tree (used to initialize brlen sets for same topology)
|
-----------------------------------------------------------------*/
int ResetBrlensFromTree (Tree *tree, Tree *vTree)
{
    int         i, j, k, nLongsNeeded, numTips;
    MrBFlt      d1, d2;
    TreeNode    *p, *q;

    if (tree->isRooted != vTree->isRooted)
        return (ERROR);
    
    if (AreTopologiesSame (tree, vTree) == NO)
        return (ERROR);

    /* allocate and set up partitions */
    AllocateTreePartitions (vTree);
    AllocateTreePartitions (tree);
    numTips = tree->nNodes - tree->nIntNodes - (tree->isRooted == YES ? 1 : 0);
    nLongsNeeded = (int) ((numTips - 1) / nBitsInALong) + 1;

    /*copy lengths and nodeDepthes*/
    for (i=0; i<vTree->nNodes; i++)
        {
        p  = vTree->allDownPass[i];
        for (j=0; j<tree->nNodes; j++)
            {
            q  = tree->allDownPass[j];
            for (k=0; k<nLongsNeeded; k++)
                if (p->partition[k] != q->partition[k])
                    break;
            if (k==nLongsNeeded)
                {
                q->length = p->length;
                if (tree->isRooted == YES)
                    q->nodeDepth = p->nodeDepth;
                }
            }
        }

    if (tree->isRooted == YES)
        {
        /*Next compute height for the root. */
        for (i=0; i<tree->nNodes-1; i++)
            {
            p  = tree->allDownPass[i];
            if (p->left == NULL)
                p->nodeDepth = 0.0;
            else
                {
                d1 = p->left->nodeDepth + p->left->length;
                d2 = p->right->nodeDepth + p->right->length;
                if (d1 > d2)
                    p->nodeDepth = d1;
                else
                    p->nodeDepth = d2;
                }
            }
        for (i=tree->nNodes-3; i>=0; i--)
            {
            p = tree->allDownPass[i];
            if (p->left==NULL && p->calibration==NULL) 
                continue;    /* leave at 0.0 */
            p->nodeDepth = p->anc->nodeDepth - p->length;
            }
        }

    FreeTreePartitions(tree);
    FreeTreePartitions(vTree);
    
    return (NO_ERROR);
}


/* ResetIntNodeIndices: Set int node indices in downpass order from numTaxa to 2*numTaxa-2 */
void ResetIntNodeIndices (PolyTree *t)
{
    int i, m, index;

    index = t->nNodes - t->nIntNodes;

    for (i=0; i<t->nIntNodes; i++)
        {
        if (t->intDownPass[i]->index != index)
            {
            SwapRelaxedBranchInfo (t->intDownPass[i]->index, index, t);
            for (m=0; m<t->nIntNodes; m++)
                {
                if (t->intDownPass[m]->index==index)
                    {
                    t->intDownPass[m]->index=t->intDownPass[i]->index;
                    break;
                    }
                }
            t->intDownPass[i]->index = index;
            }
        index++;
        }
}


/* ResetTopologyFromTree: use top to set topology in tree */
int ResetTopologyFromTree (Tree *tree, Tree *top)
{
    int         i, j, k;
    TreeNode    *p, *q, *r, *p1;

    /* adopt rooting */
    tree->isRooted = top->isRooted;
    tree->nNodes = top->nNodes;
    tree->nIntNodes = top->nIntNodes;
    
    /* set all pointers to NULL */
    for (i=0; i<tree->nNodes; i++)
        {
        p = &tree->nodes[i];
        p->anc = p->right = p->left = NULL;
        }

    /* now copy topology */
    for (i=0; i<top->nIntNodes; i++)
        {
        p1 = top->intDownPass[i];
        
        k = p1->index;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        p = &tree->nodes[j];

        k = p1->left->index;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        q = &tree->nodes[j];

        k = p1->right->index;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        r = &tree->nodes[j];

        p->left = q;
        p->right= r;
        q->anc = r->anc = p;
        }

    /* arrange the root */
    k = top->root->index;
    for (j=0; j<tree->nNodes; j++)
        if (tree->nodes[j].index == k)
            break;
    p = &tree->nodes[j];

    k = top->root->left->index;
    for (j=0; j<tree->nNodes; j++)
        if (tree->nodes[j].index == k)
            break;
    q = &tree->nodes[j];
    p->left = q;
    q->anc = p;
    p->right = p->anc = NULL;
    tree->root = p;

    GetDownPass (tree);

    return (NO_ERROR);
}


/* ResetTopologyFromPolyTree: use polytree top to set topology in tree */
int ResetTopologyFromPolyTree (Tree *tree, PolyTree *top)
{
    int         i, j, k;
    TreeNode    *p, *q, *r;
    PolyNode    *p1;

    if (tree->isRooted != top->isRooted)
        return (ERROR);
    
    /* set all pointers to NULL */
    for (i=0; i<tree->nNodes; i++)
        {
        p = &tree->nodes[i];
        p->anc = p->right = p->left = NULL;
        }

    /* now copy topology */
    for (i=0; i<top->nIntNodes; i++)
        {
        p1 = top->intDownPass[i];
        
        k = p1->index;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        p = &tree->nodes[j];

        k = p1->left->index;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        q = &tree->nodes[j];

        k = p1->left->sib->index;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        r = &tree->nodes[j];

        p->left = q;
        p->right= r;
        q->anc = r->anc = p;
        }

    /* arrange the root */
    if (top->isRooted == YES)
        {
        k = top->root->index;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        p = &tree->nodes[j];

        k = top->nNodes;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        q = &tree->nodes[j];

        q->left = p;
        q->anc = NULL;
        q->right = NULL;
        tree->root = q;
        }
    else /* if (top->isRooted == NO) */
    {
        k = top->root->index;
        for (j=0; j<tree->nNodes; j++)
            if (tree->nodes[j].index == k)
                break;
        p = &tree->nodes[j];

        k = localOutGroup;
        for (p1=top->root->left; p1!=NULL; p1=p1->sib)
            if (p1->index == k)
                break;

        assert (p1 != NULL);
        if (p1 == NULL)
            return (ERROR);

        q = &tree->nodes[p1->index];
        k = p1->anc->left->sib->sib->index;     /* index of missing child */
        if (p->left == q)
            p->left = &tree->nodes[k];
        else if (p->right == q)
            p->right = &tree->nodes[k];

        q->anc = q->right = NULL;
        p->anc = q;
        q->left = p;
    }

    GetDownPass (tree);

    return (NO_ERROR);
}


/* ResetTreePartitions: Reset bitsets describing tree partitions */
void ResetTreePartitions (Tree *t)
{
    int         i, j, numTaxa, nLongsNeeded;
    TreeNode    *p;

    /* get some handy numbers */
    numTaxa = t->nNodes - t->nIntNodes - (t->isRooted == YES ? 1 : 0);
    nLongsNeeded = (numTaxa - 1) / nBitsInALong + 1;
    
    /* reset bits describing partitions */
    for (i=0; i<t->nNodes; i++)
        {
        assert (t->allDownPass != NULL && t->allDownPass[i] != NULL);
        assert (t->allDownPass[i]->partition != NULL);
        
        p = t->allDownPass[i];
        for (j=0; j<nLongsNeeded; j++)
            p->partition[j] = 0;
        }

    /* set bits describing partitions */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL || (p->anc == NULL && t->isRooted == NO))
            SetBit (p->index, p->partition);
        else if (p->anc != NULL)
            {
            for (j=0; j<nLongsNeeded; j++)
                p->partition[j] = p->left->partition[j] | p->right->partition[j];
            }
        }
}


/*-------------------------------------------------------
|
|   RetrieveRTopology: This routine will rebuild a rooted
|      tree from the order array created by StoreRTopology.
|      All tree information except the structure will
|      remain unaltered.
|
--------------------------------------------------------*/
int RetrieveRTopology (Tree *t, int *order)
{
    int         i, numTaxa;
    TreeNode    *p, *q, *r;
    
    numTaxa = t->nNodes - t->nIntNodes - 1;
    
    /* sort the tips in the t->allDownPass array */
    p = t->nodes;
    for (i=0; i<t->nNodes; i++, p++)
        t->allDownPass[p->index] = p;

    /* make sure the root has index 2*numTaxa-1 */
    q = t->allDownPass[t->nNodes-1];
    q->anc = q->right = NULL;
    t->root = q;

    /* connect the first two tips */
    p = t->allDownPass[numTaxa];
    p->anc = q;
    q->left = p;
    p->length = 0.0;
    q = t->allDownPass[0];
    r = t->allDownPass[1];
    p->left = q;
    p->right = r;
    q->anc = r->anc = p;

    /* add one tip at a time */
    for (i=2; i<numTaxa; i++)
        {
        p = t->allDownPass[i];
        q = t->allDownPass[numTaxa-1+i];
        r = t->allDownPass[*(order++)];
        p->anc = q;
        q->left = p;
        q->right = r;
        q->anc = r->anc;
        if (r->anc->left == r)
            r->anc->left = q;
        else
            r->anc->right = q;
        r->anc = q;
        }

    /* get downpass */
    GetDownPass (t);

    /* relabel interior nodes (root is correctly labeled already) */
    for (i=0; i<t->nIntNodes; i++)
        t->intDownPass[i]->index = i+numTaxa;

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   RetrieveRTree: This routine will rebuild a rooted
|      tree from the arrays created by StoreRTree.
|      All tree information except the structure and
|      branch lengths will remain unaltered.
|
--------------------------------------------------------*/
int RetrieveRTree (Tree *t, int *order, MrBFlt *brlens)
{
    int         i, numTaxa;
    TreeNode    *p, *q, *r;

    numTaxa = t->nNodes - t->nIntNodes - 1;
    
    /* sort the tips in the t->allDownPass array */
    p = t->nodes;
    for (i=0; i<t->nNodes; i++, p++)
        t->allDownPass[p->index] = p;

    /* make sure that root has index 2*numTaxa-1 */
    q = t->allDownPass[t->nNodes-1];
    q->anc = q->right = NULL;
    q->length = 0.0;
    t->root = q;

    /* connect the first three tips */
    p = t->allDownPass[numTaxa];
    p->anc = q;
    q->left = p;
    p->length = 0.0;
    q = t->allDownPass[0];
    r = t->allDownPass[1];
    p->left = q;
    p->right = r;
    q->anc = r->anc = p;
    q->length = *(brlens++);
    r->length = *(brlens++);

    /* add one tip at a time */
    for (i=2; i<numTaxa; i++)
        {
        p = t->allDownPass[i];
        q = t->allDownPass[numTaxa-1+i];
        r = t->allDownPass[*(order++)];
        p->anc = q;
        q->left = p;
        q->right = r;
        q->anc = r->anc;
        if (r->anc->left == r)
            r->anc->left = q;
        else
            r->anc->right = q;
        r->anc = q;
        if (q->anc->anc != NULL)
            q->length = *(brlens++);
        else
            {
            r->length = *(brlens++);
            q->length = 0.0;
            }
        p->length = *(brlens++);
        }

    /* get downpass */
    GetDownPass (t);

    /* relabel interior nodes (root is correctly labeled already) */
    for (i=0; i<t->nIntNodes; i++)
        t->intDownPass[i]->index = i+numTaxa;

    /* set the node depths */
    SetNodeDepths (t);
    
    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   RetrieveRTreeWithIndices: This routine will rebuild a rooted
|      tree from the arrays created by StoreRTreeWithIndices.
|      All tree information except the structure, branch lengths
|      and node indices will remain unaltered.
|
--------------------------------------------------------*/
int RetrieveRTreeWithIndices (Tree *t, int *order, MrBFlt *brlens)
{
    int         i, numTaxa;
    TreeNode    *p, *q, *r;

    extern void ShowNodes (TreeNode *, int, int);

    numTaxa = t->nNodes - t->nIntNodes - 1;
    
    /* sort the tips in the t->allDownPass array */
    p = t->nodes;
    for (i=0; i<t->nNodes; i++, p++)
        t->allDownPass[p->index] = p;

    /* make sure that root has index 2*numTaxa-1 */
    q = t->allDownPass[t->nNodes-1];
    q->anc = q->right = NULL;
    q->length = 0.0;
    t->root = q;

    /* connect the first three 'tips' with interior node, index from order array */
    p = t->allDownPass[numTaxa];
    p->x = *(order++);
    p->anc = q;
    q->left = p;
    p->length = 0.0;
    q = t->allDownPass[0];
    r = t->allDownPass[1];
    p->left = q;
    p->right = r;
    q->anc = r->anc = p;
    q->length = *(brlens++);
    r->length = *(brlens++);

    /* add one tip at a time */
    for (i=2; i<numTaxa; i++)
        {
        p = t->allDownPass[i];
        assert (*order >= numTaxa && *order < 2*numTaxa - 1);
        q = t->allDownPass[numTaxa-1+i];
        q->x = *(order++);
        r = t->allDownPass[*(order++)];
        p->anc = q;
        q->left = p;
        q->right = r;
        q->anc = r->anc;
        if (r->anc->left == r)
            r->anc->left = q;
        else
            r->anc->right = q;
        r->anc = q;
        if (q->anc->anc != NULL)
            q->length = *(brlens++);
        else
            {
            r->length = *(brlens++);
            q->length = 0.0;
            }
        p->length = *(brlens++);
        }

    /* get downpass */
    GetDownPass (t);

    /* relabel interior nodes using labels in scratch variable x */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        p->index = p->x;
        }

    /* set the node depths */
    SetNodeDepths (t);
    
    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   RetrieveUTopology: This routine will rebuild an unrooted
|      tree from the order array created by StoreUTopology.
|      All tree information except the structure
|      will remain unaltered.
|
--------------------------------------------------------*/
int RetrieveUTopology (Tree *t, int *order)
{
    int         i, numTips;
    TreeNode    *p, *q, *r;
    
    /* preliminaries */
    numTips = t->nNodes - t->nIntNodes;
    for (i=0; i<t->nNodes; i++)
        t->nodes[i].left = t->nodes[i].right = t->nodes[i].anc = NULL;

    /* sort the tips in the t->allDownPass array */
    p = t->nodes;
    for (i=0; i<t->nNodes; i++, p++)
        t->allDownPass[p->index] = p;

    /* make sure root has index 0 */
    q = t->allDownPass[0];
    q->anc = q->right = NULL;
    t->root = q;

    /* connect the first three tips */
    p = t->allDownPass[numTips];
    p->anc = q;
    q->left = p;
    q = t->allDownPass[1];
    r = t->allDownPass[2];
    p->left = q;
    p->right = r;
    q->anc = r->anc = p;

    /* add one tip at a time */
    for (i=3; i<numTips; i++)
        {
        p = t->allDownPass[i];
        q = t->allDownPass[numTips-2+i];
        r = t->allDownPass[order[i-3]];
        p->anc = q;
        q->left = p;
        q->right = r;
        q->anc = r->anc;
        if (r->anc->left == r)
            r->anc->left = q;
        else
            r->anc->right = q;
        r->anc = q;
        }

    /* get downpass */
    GetDownPass (t);
    
    /* relabel interior nodes (root is correctly labeled already) */
    for (i=0; i<t->nIntNodes; i++)
        t->intDownPass[i]->index = i+numTips;

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   RetrieveUTree: This routine will rebuild an unrooted
|      tree from the arrays created by StoreUTree.
|      All tree information except the structure and
|      branch lengths will remain unaltered.
|
--------------------------------------------------------*/
int RetrieveUTree (Tree *t, int *order, MrBFlt *brlens)
{
    int         i, numTips;
    TreeNode    *p, *q, *r;
    
    /* preliminaries */
    numTips = t->nNodes - t->nIntNodes;
    for (i=0; i<t->nNodes; i++)
        t->nodes[i].left = t->nodes[i].right = t->nodes[i].anc = NULL;
    
    /* sort the tips in the t->allDownPass array */
    p = t->nodes;
    for (i=0; i<t->nNodes; i++, p++)
        t->allDownPass[p->index] = p;

    /* make sure that root has index 0 */
    q = t->allDownPass[0];
    q->anc = q->right = NULL;
    t->root = q;

    /* connect the first three tips */
    p = t->allDownPass[numTips];
    p->anc = q;
    q->left = p;
    p->length = *(brlens++);
    q = t->allDownPass[1];
    r = t->allDownPass[2];
    p->left = q;
    p->right = r;
    q->anc = r->anc = p;
    q->length = *(brlens++);
    r->length = *(brlens++);

    /* add one tip at a time */
    for (i=3; i<numTips; i++)
        {
        p = t->allDownPass[i];
        q = t->allDownPass[numTips-2+i];
        r = t->allDownPass[order[i-3]];
        p->anc = q;
        q->left = p;
        q->right = r;
        q->anc = r->anc;
        if (r->anc->left == r)
            r->anc->left = q;
        else
            r->anc->right = q;
        r->anc = q;
        q->length = *(brlens++);
        p->length = *(brlens++);
        }

    /* get downpass */
    GetDownPass (t);

    /* relabel interior nodes (root is correctly labeled already) */
    for (i=0; i<t->nIntNodes; i++)
        t->intDownPass[i]->index = i+numTips;

    return (NO_ERROR);
}


void SetDatedNodeAges (Param *param, int chain, int state)
{
    int         i;
    MrBFlt      clockRate;
    ModelInfo   *m;
    TreeNode    *p;
    Tree        *t;

    extern void ShowNodes(TreeNode *,int,int);

    t = GetTree (param, chain, state);
    m = &modelSettings[t->relParts[0]];

    if (m->clockRate == NULL)
        clockRate = 1.0;
    else
        clockRate = *GetParamVals(m->clockRate, chain, state);

    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (p->isDated == YES)
            p->age = p->nodeDepth / clockRate;
        else
            p->age = -1.0;
        }
}


void SetNodeDepths (Tree *t)
{
    int     i;
    MrBFlt      d1, d2;
    TreeNode    *p;

    extern void ShowNodes(TreeNode *,int,int);

    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->nodeDepth = 0.0;
        else
            {
            d1 = p->left->nodeDepth  + p->left->length;
            d2 = p->right->nodeDepth + p->right->length;
            //assert (!(t->isCalibrated == NO && AreDoublesEqual(d1,d2,0.00001)==NO)); // may not work if we set startval topology of strict clock tree by non clock tree. 
            if (d1 > d2)
                p->nodeDepth = d1;
            else
                p->nodeDepth = d2;
            }
        }

    for (i=t->nNodes-3; i>=0; i--)
        {
        p = t->allDownPass[i];
        if (p->left == NULL && p->calibration == NULL)
            p->nodeDepth = 0.0;
        else
            p->nodeDepth = p->anc->nodeDepth - p->length;
        }
}


/* Set ages of a clock tree according to depth and clockrate. Check that resulting ages are consistant with calibration.
|  return YES if tree is age consistent, No otherwise.
*/
int SetTreeNodeAges (Param *param, int chain, int state)
{
    Tree        *tree;
    TreeNode    *p;
    int         i;
    MrBFlt      clockRate;

    if (param->paramType != P_TOPOLOGY && param->paramType != P_BRLENS && param->paramType != P_SPECIESTREE)
        return YES;

    tree = GetTree(param, chain, state);
    if (modelSettings[param->relParts[0]].clockRate != NULL)
        clockRate = *GetParamVals(modelSettings[param->relParts[0]].clockRate, chain, state);
    else
        return YES;

    /* Clock trees */

    /* Check that lengths and depths are consistant. That would work for the case when we set up branch length from starting tree  */
    for (i=0; i<tree->nNodes-1; i++) {
        p = tree->allDownPass[i];
        p->age =  p->nodeDepth / clockRate;
    }

    /* Check that ages and calibrations are consistent */
    if (tree->isCalibrated == YES)
        {
        for (i=0; i<tree->nNodes-1; i++)
            {
            p = tree->allDownPass[i];
            if (p->isDated == YES) {
                if (p->calibration->prior == fixed && fabs((p->age - p->calibration->priorParams[0])/p->age) > 0.000001)
                    {
                    printf ("Node %d has age %f but should be fixed to age %f\n",
                        p->index, p->age, p->calibration->priorParams[0]);
                    return NO;
                    }
                else if (p->calibration->prior == uniform && (p->age < p->calibration->min || p->age > p->calibration->max))
                    {
                    printf ("Node %d has age %f but should be in the interval [%f,%f]\n",
                        p->index, p->age, p->calibration->min, p->calibration->max);
                    return NO;
                    }
                else if (p->age < p->calibration->min)
                    {
                    printf ("Node %d has age %f but should be minimally of age %f\n",
                        p->index, p->age, p->calibration->min);
                    return NO;
                    }
                else if (p->age > p->calibration->max)
                    {
                    printf ("Node %d has age %f but should be maximally of age %f\n",
                        p->index, p->age, p->calibration->max);
                    return NO;
                    }
                }
            }
        }

    return YES;
}


int ShowPolyNodes (PolyTree *pt)
{
    int             i;
    PolyNode        *p;

    /* this is the tree, on a node-by-node basis */
    printf ("   memnodes = %d  nNodes = %d  nIntNodes = %d  root = %d\n", pt->memNodes, pt->nNodes, pt->nIntNodes, pt->root->index);
    printf ("   isRooted = %d\n", pt->isRooted);
    printf ("   no. index (left sib anc) -- locked/free -- label (p->x)\n");
    for (i=0; i<pt->memNodes; i++)
        {
        p = &pt->nodes[i];
        if (!(p->left == NULL && p->sib == NULL && p->anc == NULL))
            {
            printf ("%4d -- %4d ", i, p->index);
            if (p->left != NULL)
                printf ("(%4d ", p->left->index);
            else
                printf ("(null ");

            if (p->sib != NULL)
                printf ("%4d ", p->sib->index);
            else
                printf ("null ");
                
            if (p->anc != NULL)
                printf ("%4d)", p->anc->index);
            else
                printf ("null)");
            
            if (p->isLocked == YES)
                printf ("-- locked -- ");
            else
                printf ("-- free --");

            if (p->left == NULL && p->anc != NULL)
                printf ("  \"%s\" (%d)\n", p->label, p->x);
            else
                printf (" \"\" (%d)\n", p->x);
            }
        }

    return NO_ERROR;
}


/* ShowTree: Show tree on screen */
int ShowTree (Tree *t)
{
    int             i, j, k, x, nLines, nLevels, levelDepth, from, to;
    char            treeLine[SCREENWIDTH2], labelLine[100];
    TreeNode        *p;
    
    /* get coordinates */
    x = 0;
    nLines = 0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL && p->right == NULL)
            {
            p->x = x;
            x += 2;
            p->y = 0;
            nLines += 2;
            }
        else if (p->left != NULL && p->right != NULL && p->anc != NULL)
            {
            p->x = p->left->x + (p->right->x - p->left->x) / 2;
            if (p->left->y > p->right->y)
                p->y = p->left->y + 1;
            else
                p->y = p->right->y + 1;
            }
        else
            {
            p->x = x;
            x += 2;
            p->y = 0;
            }
        }

    /* print tree out, line-by-line */
    levelDepth = SCREENWIDTH / t->root->left->y;
    nLevels = t->root->left->y;
    for (j=0; j<=nLines-2; j++)
        {
        for (i=0; i<SCREENWIDTH2-2; i++)
            treeLine[i] = ' ';
        treeLine[SCREENWIDTH-1] = '\n';
        if (j % 2 == 0)
            {
            for (i=0; i<t->nNodes; i++)
                {
                p = t->allDownPass[i];
                if (p->left == NULL && p->x == j)
                    {
                    strcpy (labelLine, p->label);
                    }
                }
            }
        for (i=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            if (p->anc != NULL)
                {
                if (p->anc->anc != NULL)
                    {
                    if (p->x == j)
                        {
                        from = (nLevels - p->anc->y) * levelDepth;
                        to   = (nLevels - p->y) * levelDepth;
                        if (p->y == 0)
                            to = SCREENWIDTH-1;
                        if (to >= SCREENWIDTH)
                            to = SCREENWIDTH-1;
                            
                        for (k=from; k<to; k++)
                            treeLine[k] = '-';
                        if (p->anc->left == p)
                            treeLine[from] = '/';
                        else
                            treeLine[from] = '\\';
                        if (p->left != NULL)
                            {
                            treeLine[to] = '+';
                            }
                        if (p->anc->anc == t->root && p->anc->right == p)
                            {
                            if (t->isRooted == NO)
                                treeLine[to] = '+';
                            else
                                treeLine[from] = '\\';
                            }
                        }
                    else
                        {
                        if (p->left != NULL && p->right != NULL)
                            {
                            if (j < p->x && j > p->left->x)
                                {
                                from = (nLevels - p->y) * levelDepth;
                                treeLine[from] = '|';
                                }
                            else if (j > p->x && j < p->right->x && p->left != NULL)
                                {
                                from = (nLevels - p->y) * levelDepth;
                                treeLine[from] = '|';
                                }
                            }
                        }
                    }
                else
                    {
                    if (p->x == j)
                        {
                        treeLine[0] = '|'; /* temp */
                        }
                    else if (j < p->x && j > p->left->x)
                        {
                        treeLine[0] = '|';
                        }
                    else if (j > p->x && j < p->right->x)
                        {
                        treeLine[0] = '|';
                        }
                    if (t->isRooted == NO)
                        {
                        if (j > p->x && j <= nLines-2)
                            treeLine[0] = '|';
                        if (j == p->right->x)
                            treeLine[0] = '+';
                        }
                    else
                        {
                        if (j == p->x)
                            treeLine[0] = '+';
                        }
                    }
                }
            }
        treeLine[SCREENWIDTH-1] = '\0';
        if (j % 2 == 0)
            MrBayesPrint ("   %s %s\n", treeLine, labelLine);
        else
            MrBayesPrint ("   %s \n", treeLine);
        }

    if (t->isRooted == NO)
        {
        for (i=0; i<SCREENWIDTH; i++)
            treeLine[i] = ' ';
        treeLine[SCREENWIDTH-1] = '\0';
        MrBayesPrint ("   |\n");
        for (k=0; k<SCREENWIDTH; k++)
            treeLine[k] = '-';
        treeLine[SCREENWIDTH-1] = '\0';
        treeLine[0] = '\\';
        strcpy (labelLine, t->root->label);
        labelLine[19] = '\0';
        MrBayesPrint ("   %s %s\n", treeLine, labelLine);
        }
    
#   if defined (DEBUG_CONSTRAINTS)
    for (i=0; i<t->nNodes; i++)
        printf ("%d -- %s\n", t->allDownPass[i]->index + 1, t->allDownPass[i]->isLocked == YES ? "locked" : "free");
#   endif

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreRPolyTopology: Same as StoreRTopology but for
|   binary polytree source trees.
|
--------------------------------------------------------*/
int StoreRPolyTopology (PolyTree *t, int *order)
{
    int         i, numTaxa;
    PolyNode    *p, *q;
    
    /* find number of taxa */
    numTaxa = t->nNodes - t->nIntNodes;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first two taxa */
        if (p->index > 1 && p->index < numTaxa)
            order[p->index-2] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->x = p->y = p->index;
        else
            {
            if (p->left->y < p->left->sib->y)
                {
                p->y = p->left->y;
                p->x = p->left->sib->y + numTaxa - 1;
                }
            else
                {
                p->y = p->left->sib->y;
                p->x = p->left->y + numTaxa - 1;
                }
            }
        }

    /* break the tree into pieces */
    for (i=0; i<numTaxa-2; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTaxa-3-i]];
        q = p->anc;
        if (q->left == p)
            {
            order[numTaxa-3-i] = q->left->sib->x;
            p->sib->anc = q->anc;
            if (q->anc == NULL)
                {
                p->sib->left->sib->sib = p->sib->sib;
                p->sib->sib = NULL;
                }
            else if (q->anc->left == q)
                {
                q->anc->left = q->left->sib;
                p->sib->sib = q->sib;
                }
            else
                q->anc->left->sib = q->left->sib;
            }
        else
            {
            order[numTaxa-3-i] = q->left->x;
            q->left->anc = q->anc;
            if (q->anc == NULL)
                {
                q->left->left->sib->sib = p->sib;
                q->left->sib = NULL;
                }
            else if (q->anc->left == q)
                {
                q->anc->left = q->left;
                q->anc->left->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = q->left;
                q->left->sib = NULL;
                }
            }
        }

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreRPolyTree: Same as StoreRTree but for
|      binary rooted polytree source trees.
|
--------------------------------------------------------*/
int StoreRPolyTree (PolyTree *t, int *order, MrBFlt *brlens)
{
    int         i, j, numTaxa;
    PolyNode    *p, *q;
    
    /* find number of taxa */
    numTaxa = t->nNodes - t->nIntNodes;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first two taxa */
        if (p->index > 1 && p->index < numTaxa)
            order[p->index-2] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->x = p->y = p->index;
        else
            {
            if (p->left->y < p->left->sib->y)
                {
                p->y = p->left->y;
                p->x = p->left->sib->y + numTaxa - 1;
                }
            else
                {
                p->y = p->left->sib->y;
                p->x = p->left->y + numTaxa - 1;
                }
            }
        }

    /* break the tree into pieces */
    j = t->nNodes - 2;     /* index of first branch length */
    for (i=0; i<numTaxa-2; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTaxa-3-i]];
        q = p->anc;
        brlens[j--] = p->length;
        brlens[j--] = q->length;
        if (q->left == p)
            {
            order[numTaxa-3-i] = q->left->sib->x;
            p->sib->anc = q->anc;
            if (q->anc == NULL)
                {
                p->sib->left->sib->sib = p->sib->sib;
                p->sib->sib = NULL;
                }
            else if (q->anc->left == q)
                {
                q->anc->left = q->left->sib;
                p->sib->sib = q->sib;
                }
            else
                q->anc->left->sib = q->left->sib;
            }
        else
            {
            order[numTaxa-3-i] = q->left->x;
            q->left->anc = q->anc;
            if (q->anc == NULL)
                {
                q->left->left->sib->sib = p->sib;
                q->left->sib = NULL;
                }
            else if (q->anc->left == q)
                {
                q->anc->left = q->left;
                q->anc->left->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = q->left;
                q->left->sib = NULL;
                }
            }
        }

    /* store the last two lengths; index 0 and 1 */
    p = t->root;
    brlens[p->left->index] = p->left->length;
    brlens[p->left->sib->index] = p->left->sib->length;

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreRTopology: This routine will break a rooted tree
|      into an array of ints describing the structure
|      of the tree. The tree will be destroyed
|      in the process (the node pointers, that is).
|      However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreRTopology (Tree *t, int *order)
{
    int         i, numTaxa;
    TreeNode    *p, *q;
    
    /* find number of taxa */
    numTaxa = t->nNodes - t->nIntNodes - 1;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first two taxa */
        if (p->index > 1 && p->index < numTaxa)
            order[p->index-2] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->x = p->y = p->index;
        else if (p->right != NULL)
            {
            if (p->left->y < p->right->y)
                {
                p->y = p->left->y;
                p->x = p->right->y + numTaxa - 1;
                }
            else
                {
                p->y = p->right->y;
                p->x = p->left->y + numTaxa - 1;
                }
            }
        }

    /* break the tree into pieces */
    for (i=0; i<numTaxa-2; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTaxa-3-i]];
        q = p->anc;
        if (q->left == p)
            {
            order[numTaxa-3-i] = q->right->x;
            q->right->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->right;
            else
                q->anc->right = q->right;
            }
        else
            {
            order[numTaxa-3-i] = q->left->x;
            q->left->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->left;
            else
                q->anc->right = q->left;
            }
        }

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreRTree: This routine will break a rooted tree
|      into an array of ints describing the structure
|      of the tree and an array of doubles storing
|      the branch lengths. The tree will be
|      destroyed in the process (the node pointers,
|      that is). However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreRTree (Tree *t, int *order, MrBFlt *brlens)
{
    int         i, j, numTaxa;
    TreeNode    *p, *q;

    extern void ShowNodes (TreeNode *p, int indent, int isRooted);

    /* find number of taxa */
    numTaxa = t->nNodes - t->nIntNodes - 1;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first two taxa */
        if (p->index > 1 && p->index < numTaxa)
            order[p->index-2] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->x = p->y = p->index;
        else if (p->right != NULL)
            {
            if (p->left->y < p->right->y)
                {
                p->y = p->left->y;
                p->x = p->right->y + numTaxa - 1;
                }
            else
                {
                p->y = p->right->y;
                p->x = p->left->y + numTaxa - 1;
                }
            }
        }

    /* break the tree into pieces */
    j = 2 * numTaxa - 3;
    for (i=0; i<numTaxa-2; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTaxa-3-i]];
        q = p->anc;
        brlens[j--] = p->length;
        if (q->left == p)
            {
            if (q->anc->anc != NULL)
                brlens[j--] = q->length;
            else
                brlens[j--] = q->right->length;
            order[numTaxa-3-i] = q->right->x;
            q->right->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->right;
            else
                q->anc->right = q->right;
            }
        else
            {
            if (q->anc->anc != NULL)
                brlens[j--] = q->length;
            else
                brlens[j--] = q->left->length;
            order[numTaxa-3-i] = q->left->x;
            q->left->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->left;
            else
                q->anc->right = q->left;
            }
        }

    /* store the final two branch lengths in the right order; they have indices 0 and 1 */
    p = t->root->left;
    brlens[p->left->index] = p->left->length;
    brlens[p->right->index] = p->right->length;

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreRTreeWithIndices: This routine will break a rooted
|      tree into an array of ints describing the structure
|      of the tree and the interior node indices, and an array
|      of doubles storing the branch lengths. The tree will be
|      destroyed in the process (the node pointers,
|      that is). However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreRTreeWithIndices (Tree *t, int *order, MrBFlt *brlens)
{
    int         i, j, k, numTaxa;
    TreeNode    *p, *q;

    extern void ShowNodes (TreeNode *p, int indent, int isRooted);

    /* find number of taxa */
    numTaxa = t->nNodes - t->nIntNodes - 1;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first two taxa */
        if (p->index > 1 && p->index < numTaxa)
            order[p->index-2] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->x = p->y = p->index;
        else if (p->right != NULL)
            {
            if (p->left->y < p->right->y)
                {
                p->y = p->left->y;
                p->x = p->right->y + numTaxa - 1;
                }
            else
                {
                p->y = p->right->y;
                p->x = p->left->y + numTaxa - 1;
                }
            }
        }

    /* break the tree into pieces */
    j = 2 * numTaxa - 3;
    k = 2*(numTaxa - 2);
    for (i=0; i<numTaxa-2; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTaxa-3-i]];
        q = p->anc;
        brlens[j--] = p->length;
        if (q->left == p)
            {
            if (q->anc->anc != NULL)
                brlens[j--] = q->length;
            else
                brlens[j--] = q->right->length;
            order[k--] = q->right->x;
            order[k--] = q->index;
            q->right->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->right;
            else
                q->anc->right = q->right;
            }
        else
            {
            if (q->anc->anc != NULL)
                brlens[j--] = q->length;
            else
                brlens[j--] = q->left->length;
            order[k--] = q->left->x;
            order[k--] = q->index;
            q->left->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->left;
            else
                q->anc->right = q->left;
            }
        }

    /* store the final two branch lengths in the right order; they have indices 0 and 1 */
    p = t->root->left;
    order[k] = p->index;
    brlens[p->left->index] = p->left->length;
    brlens[p->right->index] = p->right->length;

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreUPolyTopology: Same as StoreUTopology but for
|      binary polytree source.
|
--------------------------------------------------------*/
int StoreUPolyTopology (PolyTree *t, int *order)
{
    int         i, numTips;
    PolyNode    *p, *q;

    /* check if the tree is rooted on taxon 0 */
    if (t->root->left->sib->sib->index != 0)
        MovePolyCalculationRoot (t, 0);

    /* rearrange the root */
    t->root->anc = t->root->left->sib->sib;
    t->root->left->sib->sib = NULL;
    t->root->anc->left = t->root;
    t->root->anc->sib = NULL;
    t->root->anc->anc = NULL;
    t->root = t->root->anc;

    /* find number of tips */
    numTips = t->nNodes - t->nIntNodes;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first three taxa */
        if (p->index > 2 && p->index < numTips)
            order[p->index-3] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL || p->anc == NULL)
            p->x = p->y = p->index;
        else
            {
            if (p->left->y < p->left->sib->y)
                {
                p->y = p->left->y;
                p->x = p->left->sib->y + numTips - 2;
                }
            else
                {
                p->y = p->left->sib->y;
                p->x = p->left->y + numTips - 2;
                }
            }
        }

    /* break the tree into pieces */
    for (i=0; i<numTips-3; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTips-4-i]];
        q = p->anc;
        if (q->left == p)
            {
            order[numTips-4-i] = q->left->sib->x;
            p->sib->anc = q->anc;
            if (q->anc->left == q)
                {
                q->anc->left = p->sib;
                p->sib->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = p->sib;
                p->sib->sib = q->sib;
                }
            }
        else
            {
            order[numTips-4-i] = q->left->x;
            q->left->anc = q->anc;
            if (q->anc->left == q)
                {
                q->anc->left = q->left;
                q->left->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = q->left;
                q->left->sib = q->sib;
                }
            }
        }

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreUPolyTree: Same as StoreUTopology but for
|      binary polytree source.
|
--------------------------------------------------------*/
int StoreUPolyTree (PolyTree *t, int *order, MrBFlt *brlens)
{
    int         i, j, numTips;
    PolyNode    *p, *q;

    /* check if the tree is rooted on taxon 0 */
    if (t->root->left->sib->sib->index != 0)
        MovePolyCalculationRoot (t, 0);

    /* rearrange the root */
    t->root->anc = t->root->left->sib->sib;
    t->root->left->sib->sib = NULL;
    t->root->anc->left = t->root;
    t->root->anc->sib = NULL;
    t->root->anc->anc = NULL;
    t->root = t->root->anc;

    /* find number of tips */
    numTips = t->nNodes - t->nIntNodes;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first three taxa */
        if (p->index > 2 && p->index < numTips)
            order[p->index-3] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL || p->anc == NULL)
            p->x = p->y = p->index;
        else
            {
            if (p->left->y < p->left->sib->y)
                {
                p->y = p->left->y;
                p->x = p->left->sib->y + numTips - 2;
                }
            else
                {
                p->y = p->left->sib->y;
                p->x = p->left->y + numTips - 2;
                }
            }
        }

    /* break the tree into pieces */
    j = 2*numTips - 4;
    for (i=0; i<numTips-3; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTips-4-i]];
        assert (p->index > 2 && p->index < numTips);
        assert (p->anc->anc != NULL);
        q = p->anc;
        brlens[j--] = p->length;
        brlens[j--] = q->length;
        if (q->left == p)
            {
            order[numTips-4-i] = q->left->sib->x;
            p->sib->anc = q->anc;
            if (q->anc->left == q)
                {
                q->anc->left = p->sib;
                p->sib->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = p->sib;
                p->sib->sib = q->sib;
                }
            }
        else
            {
            order[numTips-4-i] = q->left->x;
            q->left->anc = q->anc;
            if (q->anc->left == q)
                {
                q->anc->left = q->left;
                q->left->sib = q->sib;
                }
            else
                {
                q->anc->left->sib = q->left;
                q->left->sib = q->sib;
                }
            }
        }

    /* store last three branch lengths, index 0, 1, 2 */
    q = t->root;
    assert (q->index == 0);
    brlens[q->index] = q->length;
    q = q->left->left;
    assert (q->index == 1 || q->index == 2);
    brlens[q->index] = q->length;
    q = q->sib;
    assert (q->index == 1 || q->index == 2);
    brlens[q->index] = q->length;

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreUTopology: This routine will break an unrooted tree
|      into an array of ints describing the structure
|      of the tree. The tree will be destroyed
|      in the process (the node pointers, that is).
|      However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreUTopology (Tree *t, int *order)
{
    int         i, numTips;
    TreeNode    *p, *q;

    /* check if the tree is rooted on taxon 0 */
    if (t->root->index != 0)
        MoveCalculationRoot (t, 0);

    /* find number of tips */
    numTips = t->nNodes - t->nIntNodes;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first three taxa */
        if (p->index > 2 && p->index < numTips)
            order[p->index-3] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->x = p->y = p->index;
        else if (p->right != NULL)
            {
            if (p->left->y < p->right->y)
                {
                p->y = p->left->y;
                p->x = p->right->y + numTips - 2;
                }
            else
                {
                p->y = p->right->y;
                p->x = p->left->y + numTips - 2;
                }
            }
        }

    /* break the tree into pieces */
    for (i=0; i<numTips-3; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTips-4-i]];
        q = p->anc;
        if (q->left == p)
            {
            order[numTips-4-i] = q->right->x;
            q->right->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->right;
            else
                q->anc->right = q->right;
            }
        else
            {
            order[numTips-4-i] = q->left->x;
            q->left->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->left;
            else
                q->anc->right = q->left;
            }
        }

    return (NO_ERROR);
}


/*-------------------------------------------------------
|
|   StoreUTree: This routine will break an unrooted tree
|      into an array of ints describing the structure
|      of the tree and an array of doubles storing
|      the branch lengths. The tree will be
|      destroyed in the process (the node pointers,
|      that is). However, the tree is not deleted.
|
--------------------------------------------------------*/
int StoreUTree (Tree *t, int *order, MrBFlt *brlens)
{
    int         i, j, numTips;
    TreeNode    *p, *q;

    /* check if the tree is rooted on taxon 0 */
    if (t->root->index != 0)
        MoveCalculationRoot(t, 0);

    /* find number of tips */
    numTips = t->nNodes - t->nIntNodes;

    /* first get the terminal taxon positions and store
       them in the order array. */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        /* we do not need to worry about the first three taxa */
        if (p->index > 2 && p->index < numTips)
            order[p->index-3] = i;
        }

    /* label the interior nodes with the correct index */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            p->x = p->y = p->index;
        else if (p->right != NULL)
            {
            if (p->left->y < p->right->y)
                {
                p->y = p->left->y;
                p->x = p->right->y + numTips - 2;
                }
            else
                {
                p->y = p->right->y;
                p->x = p->left->y + numTips - 2;
                }
            }
        }

    /* break the tree into pieces */
    j = 2 * numTips - 4;
    for (i=0; i<numTips-3; i++)
        {
        /* find the next node to remove */
        p = t->allDownPass[order[numTips-4-i]];
        q = p->anc;
        brlens[j--] = p->length;
        brlens[j--] = q->length;
        if (q->left == p)
            {
            order[numTips-4-i] = q->right->x;
            q->right->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->right;
            else
                q->anc->right = q->right;
            }
        else
            {
            order[numTips-4-i] = q->left->x;
            q->left->anc = q->anc;
            if (q->anc->left == q)
                q->anc->left = q->left;
            else
                q->anc->right = q->left;
            }
        }

    /* store the final three branch lengths */
    /* we need to check the rotation of the tree to 
           store the brlens in the right order (after node index) */
    p = t->root->left;
    if (p->right->index == 2)
            {
            brlens[j--] = p->right->length;
        brlens[j--] = p->left->length;
            }
    else
            {
            brlens[j--] = p->left->length;
            brlens[j--] = p->right->length;
            }
    brlens[j--] = p->length;

    return (NO_ERROR);
}


/* TreeLength: Calculate tree length */
MrBFlt TreeLen (Tree *t)
{
    int     i, numLenNodes;
    MrBFlt  len = 0.0;

    if (t->isRooted == NO)
        numLenNodes = t->nNodes - 1;
    else
        numLenNodes = t->nNodes - 2;

    for (i=0; i<numLenNodes; i++)
        len += t->allDownPass[i]->length;

    return len;
}


/*-------------------------------------------------------------------------------------------
|
|   Unmark: This routine will unmark a subtree rooted at p
|
---------------------------------------------------------------------------------------------*/
void Unmark (TreeNode *p)
{
    if (p != NULL)
        {
        p->marked = NO;
        Unmark (p->left);
        Unmark (p->right);
        }
}


void WriteEventTree (TreeNode *p, int chain, Param *param)
{
    int             j, nEvents;
    MrBFlt          brlen, *position, *rateMult;

    if (p != NULL)
        {
        if (p->left == NULL && p->right == NULL)
            {
            printf ("%d:%s", p->index + 1, MbPrintNum(p->length));
            if (param->paramType == P_CPPEVENTS)
                {
                nEvents = param->nEvents[2*chain+state[chain]][p->index];
                if (nEvents > 0)
                    {
                    printf ("[&E %s %d: (", param->name, nEvents);
                    position = param->position[2*chain+state[chain]][p->index];
                    rateMult = param->rateMult[2*chain+state[chain]][p->index];
                    for (j=0; j<nEvents; j++)
                        {
                        printf ("%s %s", MbPrintNum(position[j]), MbPrintNum(rateMult[j]));
                        if (j != nEvents-1)
                            printf (", ");
                        }
                    printf (")]");
                    }
                else
                    printf ("[&E %s 0]", param->name);
                }
            brlen = GetParamSubVals (param, chain, state[chain])[p->index];
            // brlen = (GetParamSubVals (param, chain, state[chain])[p->index] + GetParamVals (param, chain, state[chain])[p->anc->index]) / 2.0;
            printf ("[&B %s %s]", param->name, MbPrintNum(brlen));
            }
        else
            {
            if (p->anc != NULL)
                printf ("(");
            WriteEventTree(p->left, chain, param);
            printf (",");
            WriteEventTree(p->right, chain, param);
            if (p->anc != NULL)
                {               
                if (p->anc->anc != NULL)
                    {
                    printf ("):%s", MbPrintNum(p->length));
                    if (param->paramType == P_CPPEVENTS)
                        {
                        nEvents = param->nEvents[2*chain+state[chain]][p->index];
                        if (nEvents > 0)
                            {
                            printf ("[&E %s %d: (", param->name, nEvents);
                            position = param->position[2*chain+state[chain]][p->index];
                            rateMult = param->rateMult[2*chain+state[chain]][p->index];
                            for (j=0; j<nEvents; j++)
                                {
                                printf ("%s %s", MbPrintNum(position[j]), MbPrintNum(rateMult[j]));
                                if (j != nEvents-1)
                                    printf (", ");
                                }
                            printf (")]");
                            }
                        else
                            printf ("[&E %s 0]", param->name);
                        }
                    brlen = GetParamSubVals (param, chain, state[chain])[p->index];
                    // brlen = (GetParamSubVals (param, chain, state[chain])[p->index] + GetParamVals (param, chain, state[chain])[p->anc->index]) / 2.0;
                    printf ("[&B %s %s]", param->name, MbPrintNum(brlen));
                    }
                else
                    printf (")");
                }
            }
        }
}


void WriteEventTreeToPrintString (TreeNode *p, int chain, Param *param, int printAll)
{
    char            *tempStr;
    int             i, j, nEvents, tempStrSize = TEMPSTRSIZE;
    MrBFlt          brlen, *position, *rateMult;

    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));

    if (p != NULL)
        {
        if (p->left == NULL && p->right == NULL)
            {
            SafeSprintf (&tempStr, &tempStrSize, "%d:%s", p->index + 1, MbPrintNum(p->length));
            AddToPrintString (tempStr);
            for (i=0; i<param->nSubParams; i++)
                {
                if (param->subParams[i]->paramType == P_CPPEVENTS)
                    {
                    nEvents = param->subParams[i]->nEvents[2*chain+state[chain]][p->index];
                    if (nEvents > 0)
                        {
                        SafeSprintf (&tempStr, &tempStrSize, "[&E %s %d", param->subParams[i]->name, nEvents);
                        AddToPrintString (tempStr);
                        position = param->subParams[i]->position[2*chain+state[chain]][p->index];
                        rateMult = param->subParams[i]->rateMult[2*chain+state[chain]][p->index];
                        if (printAll == YES)
                            {
                            SafeSprintf (&tempStr, &tempStrSize, ": (");
                            AddToPrintString (tempStr);
                            for (j=0; j<nEvents; j++)
                                {
                                SafeSprintf (&tempStr, &tempStrSize, "%s", MbPrintNum(position[j]));
                                AddToPrintString (tempStr);
                                SafeSprintf (&tempStr, &tempStrSize, " %s",  MbPrintNum(rateMult[j]));
                                AddToPrintString (tempStr);
                                if (j != nEvents-1)
                                    AddToPrintString (",");
                                else
                                    AddToPrintString (")");
                                }
                            }
                        AddToPrintString ("]");
                        }
                    else
                        {
                        SafeSprintf (&tempStr, &tempStrSize, "[&E %s 0]", param->subParams[i]->name);
                        AddToPrintString (tempStr);
                        }
                    }
                else if (param->subParams[i]->paramType != P_CPPEVENTS)
                    {
                    /* other relaxed clock models */
                    brlen = GetParamSubVals (param->subParams[i], chain, state[chain])[p->index];
                    SafeSprintf (&tempStr, &tempStrSize, "[&B %s %s]", param->subParams[i]->name, MbPrintNum(brlen));
                    AddToPrintString (tempStr);
                    }
                }
            }
        else
            {
            if (p->anc != NULL)
                AddToPrintString ("(");
            WriteEventTreeToPrintString (p->left, chain, param, printAll);
            AddToPrintString (",");
            WriteEventTreeToPrintString (p->right, chain, param, printAll); 
            if (p->anc != NULL)
                {               
                if (p->anc->anc != NULL)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "):%s", MbPrintNum(p->length));
                    AddToPrintString (tempStr);
                    for (i=0; i<param->nSubParams; i++)
                        {
                        if (param->subParams[i]->paramType == P_CPPEVENTS)
                            {
                            nEvents = param->subParams[i]->nEvents[2*chain+state[chain]][p->index];
                            if (nEvents > 0)
                                {
                                SafeSprintf (&tempStr, &tempStrSize, "[&E %s %d", param->subParams[i]->name, nEvents);
                                AddToPrintString (tempStr);
                                position = param->subParams[i]->position[2*chain+state[chain]][p->index];
                                rateMult = param->subParams[i]->rateMult[2*chain+state[chain]][p->index];
                                if (printAll == YES)
                                    {
                                    SafeSprintf (&tempStr, &tempStrSize, ": (");
                                    AddToPrintString (tempStr);
                                    for (j=0; j<nEvents; j++)
                                        {
                                        SafeSprintf (&tempStr, &tempStrSize, "%s", MbPrintNum(position[j]));
                                        AddToPrintString (tempStr);
                                        SafeSprintf (&tempStr, &tempStrSize, " %s",  MbPrintNum(rateMult[j]));
                                        AddToPrintString (tempStr);
                                        if (j != nEvents-1)
                                            AddToPrintString (",");
                                        else
                                            AddToPrintString (")");
                                        }
                                    }
                                AddToPrintString ("]");
                                }
                            else
                                {
                                SafeSprintf (&tempStr, &tempStrSize, "[&E %s 0]", param->subParams[i]->name);
                                AddToPrintString (tempStr);
                                }
                            }
                        else if (param->subParams[i]->paramType != P_CPPEVENTS)
                            {
                            /* other relaxed clock models */
                            brlen = GetParamSubVals (param->subParams[i], chain, state[chain])[p->index];
                            SafeSprintf (&tempStr, &tempStrSize, "[&B %s %s]", param->subParams[i]->name, MbPrintNum(brlen));
                            AddToPrintString (tempStr);
                            }
                        }
                    }
                else
                    AddToPrintString(")");
                }
            }
        }
    free (tempStr);
}


void WriteEvolTree (TreeNode *p, int chain, Param *param)
{
    MrBFlt          *length;

    if (p != NULL)
        {
        length = GetParamSubVals(param, chain, state[chain]);
        if (p->left == NULL && p->right == NULL)
            {
            printf ("%d:%s", p->index + 1, MbPrintNum(length[p->index]));
            }
        else
            {
            if (p->anc != NULL)
                printf ("(");
            WriteEvolTree(p->left, chain, param);
            printf (",");
            WriteEvolTree(p->right, chain, param);
            if (p->anc != NULL)
                {               
                if (p->anc->anc != NULL)
                    printf ("):%s", MbPrintNum(length[p->index]));
                else
                    printf (")");
                }
            }
        }
}


void WriteNoEvtTreeToPrintString (TreeNode *p, int chain, Param *param, int showBrlens, int isRooted)
{
    char            *tempStr;
    int             i, tempStrSize = TEMPSTRSIZE, nEvents;
    MrBFlt          brlen, N;

    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));

    if (p != NULL)
        {
        if (p->left == NULL && p->right == NULL)
            {
            if (showBrlens == YES)
                {
                SafeSprintf (&tempStr, &tempStrSize, "%d:%s", p->index + 1, MbPrintNum(p->length));
                }
            else
                SafeSprintf (&tempStr, &tempStrSize, "%d", p->index + 1);
            AddToPrintString (tempStr);
            if (param->paramType == P_BRLENS)
                {
                for (i=0; i<param->nSubParams; i++)
                    {
                    if (param->subParams[i]->paramType == P_CPPEVENTS)
                        {
                        nEvents = param->subParams[i]->nEvents[2*chain+state[chain]][p->index];
                        SafeSprintf (&tempStr, &tempStrSize, "[&E %s %d]", param->subParams[i]->name, nEvents);
                        AddToPrintString (tempStr);
                        }
                    brlen = GetParamSubVals (param->subParams[i], chain, state[chain])[p->index];
                    SafeSprintf (&tempStr, &tempStrSize, "[&B %s %s]", param->subParams[i]->name, MbPrintNum(brlen));
                    AddToPrintString (tempStr);
                    }
                }
            else if (param->paramType == P_SPECIESTREE && modelSettings[param->relParts[0]].popSize->nValues > 1)
                {
                N = GetParamVals (modelSettings[param->relParts[0]].popSize, chain, state[chain])[p->index];
                SafeSprintf (&tempStr, &tempStrSize, "[&N %s %s]", modelSettings[param->relParts[0]].popSize->name, MbPrintNum(N));
                AddToPrintString (tempStr);
                }
            }
        else
            {
            if (p->anc != NULL)
                AddToPrintString ("(");
            WriteNoEvtTreeToPrintString (p->left,  chain, param, showBrlens, isRooted);
            if (p->anc != NULL)
                AddToPrintString (",");
            WriteNoEvtTreeToPrintString (p->right, chain, param, showBrlens, isRooted);
            if (p->anc != NULL)
                {
                if (p->anc->anc == NULL && isRooted == NO)
                    {
                    if (showBrlens == YES)
                        SafeSprintf (&tempStr, &tempStrSize, ",%d:%s)", p->anc->index + 1, MbPrintNum(p->length));
                    else
                        SafeSprintf (&tempStr, &tempStrSize, ",%d)", p->anc->index + 1);
                    AddToPrintString (tempStr);
                    }
                else if (p->anc->anc != NULL)
                    {
                    if (showBrlens == YES)
                        SafeSprintf (&tempStr, &tempStrSize, "):%s", MbPrintNum(p->length));
                    else
                        SafeSprintf (&tempStr, &tempStrSize, ")");
                    AddToPrintString (tempStr);
                    if (param->paramType == P_BRLENS)
                        {
                        for (i=0; i<param->nSubParams; i++)
                            {
                            if (param->subParams[i]->paramType == P_CPPEVENTS)
                                {
                                nEvents = param->subParams[i]->nEvents[2*chain+state[chain]][p->index];
                                SafeSprintf (&tempStr, &tempStrSize, "[&E %s %d]", param->subParams[i]->name, nEvents);
                                AddToPrintString (tempStr);
                                }
                            brlen = GetParamSubVals (param->subParams[i], chain, state[chain])[p->index];
                            SafeSprintf (&tempStr, &tempStrSize, "[&B %s %s]", param->subParams[i]->name, MbPrintNum(brlen));
                            AddToPrintString (tempStr);
                            }
                        }
                    else if (param->paramType == P_SPECIESTREE && modelSettings[param->relParts[0]].popSize->nValues > 1)
                        {
                        N = GetParamVals (modelSettings[param->relParts[0]].popSize, chain, state[chain])[p->index];
                        SafeSprintf (&tempStr, &tempStrSize, "[&N %s %s]", modelSettings[param->relParts[0]].popSize->name, MbPrintNum(N));
                        AddToPrintString (tempStr);
                        }
                    }
                else if (param->paramType == P_SPECIESTREE && modelSettings[param->relParts[0]].popSize->nValues > 1)
                    {
                    N = GetParamVals (modelSettings[param->relParts[0]].popSize, chain, state[chain])[p->index];
                    SafeSprintf (&tempStr, &tempStrSize, ")[&N %s %s]", modelSettings[param->relParts[0]].popSize->name, MbPrintNum(N));
                    AddToPrintString (tempStr);
                    }
                else
                    AddToPrintString(")");
                }
            }
        }
    free (tempStr);
}


/* WriteTopologyToFile: Simply write topology to file */
void WriteTopologyToFile (FILE *fp, TreeNode *p, int isRooted)
{
    if (p != NULL)
        {
        if (p->left == NULL && p->right == NULL)
            fprintf (fp, "%d", p->index + 1);
        else
            {
            if (p->anc != NULL)
                fprintf (fp, "(");
            WriteTopologyToFile (fp, p->left, isRooted);
            if (p->anc != NULL)
                fprintf (fp, ",");
            WriteTopologyToFile (fp, p->right, isRooted);   
            if (p->anc != NULL)
                {
                if (p->anc->anc == NULL && isRooted == NO)
                    fprintf (fp, ",%d", p->anc->index + 1);
                fprintf (fp, ")");
                }
            }
        }
}


/* the following are moved from mbmath.c */
/*---------------------------------------------------------------------------------
|
|   AddTwoMatrices
|
|   Takes the sum of two matrices, "a" and "b", and puts the results in a matrix
|   called "result".
|
---------------------------------------------------------------------------------*/
void AddTwoMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result)
{
    int         row, col;

    for (row=0; row<dim; row++)
        {
        for (col=0; col<dim; col++) 
            {
            result[row][col] = a[row][col] + b[row][col];
            }
        }
}


/*---------------------------------------------------------------------------------
|
|   AllocateSquareComplexMatrix
|
|   Allocate memory for a square (dim X dim) complex matrix.
|
---------------------------------------------------------------------------------*/
complex **AllocateSquareComplexMatrix (int dim)
{
    int         i;
    complex     **m;

    m = (complex **) SafeMalloc ((size_t)dim * sizeof(complex*));
    if (!m) 
        {
        MrBayesPrint ("%s   Error: Problem allocating a square complex matrix.\n", spacer);
        exit (0);
        }
    m[0]=(complex *) SafeMalloc ((size_t)dim * (size_t)dim *sizeof(complex));
    if (!m[0]) 
        {
        MrBayesPrint ("%s   Error: Problem allocating a square complex matrix.\n", spacer);
        exit (0);
        }
    for (i=1;i<dim;i++) 
        {
        m[i] = m[i-1] + dim;
        }
        
    return (m);
}


/*---------------------------------------------------------------------------------
|
|   AllocateSquareDoubleMatrix
|
|   Allocate memory for a square (dim X dim) matrix of doubles.
|
---------------------------------------------------------------------------------*/
MrBFlt **AllocateSquareDoubleMatrix (int dim)
{
    int         i;
    MrBFlt      **m;
    
    m = (MrBFlt **) SafeMalloc ((size_t)dim * sizeof(MrBFlt*));
    if (!m)
        {
        MrBayesPrint ("%s   Error: Problem allocating a square matrix of doubles.\n", spacer);
        exit(1);
        }
    m[0] = (MrBFlt *) SafeMalloc ((size_t)dim * (size_t)dim * sizeof(MrBFlt));
    if (!m[0])
        {
        MrBayesPrint ("%s   Error: Problem allocating a square matrix of doubles.\n", spacer);
        exit(1);
        }
    for (i=1; i<dim; i++)
        {
        m[i] = m[i-1] + dim;
        }

    return (m);
}


/*---------------------------------------------------------------------------------
|
|   AllocateSquareIntegerMatrix
|
|   Allocate memory for a square (dim X dim) matrix of integers.
|
---------------------------------------------------------------------------------*/
int **AllocateSquareIntegerMatrix (int dim)
{
    int     i, **m;
    
    m = (int **) SafeMalloc ((size_t)dim * sizeof(int*));
    if (!m)
        {
        MrBayesPrint ("%s   Error: Problem allocating a square matrix of integers.\n", spacer);
        exit(1);
        }
    m[0] = (int *) SafeMalloc ((size_t)dim * (size_t)dim * sizeof(int));
    if (!m[0])
        {
        MrBayesPrint ("%s   Error: Problem allocating a square matrix of integers.\n", spacer);
        exit(1);
        }
    for (i=1; i<dim; i++)
        {
        m[i] = m[i-1] + dim;
        }

    return (m);
}


/*---------------------------------------------------------------------------------
|
|   AutodGamma
|
|   Auto-discrete-gamma distribution of rates over sites, K equal-probable
|   categories, with the mean for each category used.         
|   This routine calculates M[], using rho and K (numGammaCats)     
|
---------------------------------------------------------------------------------*/
int AutodGamma (MrBFlt *M, MrBFlt rho, int K)
{
    int         i, j, i1, i2;
    MrBFlt      point[MAX_GAMMA_CATS], x, y, large = 20.0, sum;

    for (i=0; i<K-1; i++) 
        point[i] = PointNormal ((i + 1.0) / K);
    for (i=0; i<K; i++) 
        {
        for (j=0; j<K; j++) 
            {
            x = (i < K-1 ? point[i]:large);
            y = (j < K-1 ? point[j]:large);
            M[i * K + j] = CdfBinormal (x, y, rho);
            }
        }
    for (i1=0; i1<2*K-1; i1++) 
        {
        for (i2=0; i2<K*K; i2++) 
            {
            i = i2 / K; 
            j = i2 % K;
            if (AreDoublesEqual(i+j, 2*(K-1.0)-i1, ETA)==NO)
                continue;
            y = 0;
            if (i > 0) 
                y -= M[(i-1)*K+j];
            if (j > 0) 
                y -= M[i*K+(j-1)];
            if (i > 0 && j > 0) 
                y += M[(i-1)*K+(j-1)];
            M[i*K+j] = (M[i*K+j] + y) * K;
            }
        }
    for (i=0; i<K; i++)
        {
        sum = 0.0;
        for (j=0; j<K; j++)
            {
            if (M[i*K+j] < 0.0)
                M[i*K+j] = 0.0;
            sum += M[i*K+j];
            }
        for (j=0; j<K; j++)
            M[i*K+j] /= sum;
        }
    
//    MrBayesPrint ("rho = %lf\n", rho);
//    for (i=0; i<K; i++)
//        {
//        for (j=0; j<K; j++)
//            MrBayesPrint ("%lf ", M[i*K + j]);
//        MrBayesPrint ("\n");
//        }
    
    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
|
|   BackSubstitutionRow
|
---------------------------------------------------------------------------------*/
void BackSubstitutionRow (int dim, MrBFlt **u, MrBFlt *b)
{
    int             i, j;
    MrBFlt          dotProduct;

    b[dim-1] /= u[dim-1][dim-1];
    for (i=dim-2; i>=0; i--) 
        {
        dotProduct = 0.0;
        for (j=i+1; j<dim; j++)
            dotProduct += u[i][j] * b[j];
        b[i] = (b[i] - dotProduct) / u[i][i];
        }
}


/*---------------------------------------------------------------------------------
|
|   Balanc
|
|   This subroutine balances a real matrix and isolates
|   eigenvalues whenever possible.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * a contains the input matrix to be balanced
|
|   On output:
|
|    * a contains the balanced matrix.
|
|    * low and high are two integers such that a(i,j)
|      is equal to zero if
|         (1) i is greater than j and
|         (2) j=1,...,low-1 or i=igh+1,...,n.
|
|    * scale contains information determining the
|      permutations and scaling factors used.
|
|   Suppose that the principal submatrix in rows pLow through pHigh
|   has been balanced, that p(j) denotes the index interchanged
|   with j during the permutation step, and that the elements
|   of the diagonal matrix used are denoted by d(i,j). Then
|      scale(j) = p(j),    for j = 1,...,pLow-1
|               = d(j,j),      j = pLow,...,pHigh
|               = p(j)         j = pHigh+1,...,dim.
|   The order in which the interchanges are made is dim to pHigh+1,
|   then 1 to pLow-1.
|
|   Note that 1 is returned for pHigh if pHigh is zero formally.
|
|   The algol procedure exc contained in balance appears in
|   balanc in line.  (Note that the algol roles of identifiers
|   k,l have been reversed.)
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|
|   This function was converted from FORTRAN by D. L. Swofford.
|   
---------------------------------------------------------------------------------*/
void Balanc (int dim, MrBFlt **a, int *low, int *high, MrBFlt *scale)
{
    int         i, j, k, l, m, noconv;
    MrBFlt      c, f, g, r, s, b2;

    b2 = FLT_RADIX * FLT_RADIX;
    k = 0;
    l = dim - 1;
    
    for (j=l; j>=0; j--)
        {
        for (i=0; i<=l; i++)
            {
            if (i != j)
                {
                  if (AreDoublesEqual(a[j][i],0.0, ETA)==NO)
                    goto next_j1;
                }
            }
            
        /* bug that DLS caught */
        /*m = l;
        Exchange(j, k, l, m, dim, a, scale);
        if (l < 0)
            goto leave;
        else
            j = --l;*/
        m = l;
        Exchange(j, k, l, m, dim, a, scale);
        if (--l < 0)
            goto leave;
        next_j1:
            ;
        }

    for (j=k; j<=l; j++)
        {
        for (i=k; i<=l; i++)
            {
            if (i != j)
                {
                if (AreDoublesEqual(a[i][j], 0.0, ETA)==NO)
                    goto next_j;
                }
            }
        m = k;
        Exchange(j, k, l, m, dim, a, scale);
        k++;
        next_j:
            ;
        }

    for (i=k; i<=l; i++)
        scale[i] = 1.0;
    
    do  {
        noconv = FALSE;
        for (i=k; i<=l; i++)
            {
            c = 0.0;
            r = 0.0;
            for (j=k; j<=l; j++)
                {
                if (j != i)
                    {
                    c += fabs(a[j][i]);
                    r += fabs(a[i][j]);
                    }
                }
            if (AreDoublesEqual(c,0.0,ETA)==NO && AreDoublesEqual(r,0.0,ETA)==NO)
                {
                g = r / FLT_RADIX;
                f = 1.0;
                s = c + r;
                while (c < g)
                    {
                    f *= FLT_RADIX;
                    c *= b2;
                    }
                g = r * FLT_RADIX;
                while (c >= g)
                    {
                    f /= FLT_RADIX;
                    c /= b2;
                    }
                if ((c + r) / f < s * .95)
                    {
                    g = 1.0 / f;
                    scale[i] *= f;
                    noconv = TRUE;              
                    for (j=k; j<dim; j++)
                        a[i][j] *= g;
                    for (j=0; j<=l; j++)
                        a[j][i] *= f;
                    }
                }
            }   
        }
        while (noconv);
    leave:
        *low = k;
        *high = l;
    
#   if 0 
/* begin f2c version of code:
   balanc.f -- translated by f2c (version 19971204) */
int balanc (int *nm, int *n, MrBFlt *a, int *low, int *igh, MrBFlt *scale)

{

    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;
    MrBFlt d__1;

    /* Local variables */
    static MrBFlt iexc;
    static MrBFlt c__, f, g;
    static MrBFlt i__, j, k, l, m;
    static MrBFlt r__, s, radix, b2;
    static MrBFlt jj;
    static logical noconv;

    /* parameter adjustments */
    --scale;
    a_dim1 = *nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* function Body */
    radix = 16.0;

    b2 = radix * radix;
    k = 1;
    l = *n;
    goto L100;
    
    /* .......... in-line procedure for row and column exchange .......... */
    L20:
    scale[m] = (MrBFlt) j;
    if (j == m) 
        goto L50;

    i__1 = l;
    for (i__ = 1; i__ <= i__1; ++i__) 
        {
        f = a[i__ + j * a_dim1];
        a[i__ + j * a_dim1] = a[i__ + m * a_dim1];
        a[i__ + m * a_dim1] = f;
        /* L30: */
        }

    i__1 = *n;
    for (i__ = k; i__ <= i__1; ++i__) 
        {
        f = a[j + i__ * a_dim1];
        a[j + i__ * a_dim1] = a[m + i__ * a_dim1];
        a[m + i__ * a_dim1] = f;
        /* L40: */
        }

    L50:
    switch (iexc) 
        {
        case 1:  
            goto L80;
        case 2:  
            goto L130;
        }
        
    /* .......... search for rows isolating an eigenvalue and push them down .......... */
    L80:
    if (l == 1) 
        goto L280;
    --l;
    
    /* .......... for j=l step -1 until 1 do -- .......... */
    L100:
    i__1 = l;
    for (jj = 1; jj <= i__1; ++jj) 
        {
        j = l + 1 - jj;
        i__2 = l;
        for (i__ = 1; i__ <= i__2; ++i__) 
            {
            if (i__ == j) 
                goto L110;
            if (a[j + i__ * a_dim1] != 0.) 
                goto L120;
            L110:
            ;
            }
        m = l;
        iexc = 1;
        goto L20;
        L120:
        ;
        }

    goto L140;
    /* .......... search for columns isolating an eigenvalue and push them left .......... */
    L130:
    ++k;

    L140:
    i__1 = l;
    for (j = k; j <= i__1; ++j) 
        {
        i__2 = l;
        for (i__ = k; i__ <= i__2; ++i__) 
            {
            if (i__ == j) 
                goto L150;
            if (a[i__ + j * a_dim1] != 0.) 
                goto L170;
            L150:
            ;
            }
        m = k;
        iexc = 2;
        goto L20;
        L170:
        ;
        }
        
    /* .......... now balance the submatrix in rows k to l .......... */
    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) 
        {
        /* L180: */
        scale[i__] = 1.0;
        }
    /* .......... iterative loop for norm reduction .......... */
    L190:
    noconv = FALSE;

    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) 
        {
        c__ = 0.0;
        r__ = 0.0;
        i__2 = l;
        for (j = k; j <= i__2; ++j) 
            {
            if (j == i__) 
                goto L200;
            c__ += (d__1 = a[j + i__ * a_dim1], abs(d__1));
            r__ += (d__1 = a[i__ + j * a_dim1], abs(d__1));
            L200:
            ;
            }
        
        /* .......... guard against zero c or r due to underflow .......... */
        if (c__ == 0. || r__ == 0.) 
            goto L270;
        g = r__ / radix;
        f = 1.0;
        s = c__ + r__;
        L210:
        if (c__ >= g) 
            goto L220;
        f *= radix;
        c__ *= b2;
        goto L210;
        L220:
        g = r__ * radix;
        L230:
        if (c__ < g) 
            goto L240;
        f /= radix;
        c__ /= b2;
        goto L230;
        
        /*     .......... now balance .......... */
        L240:
        if ((c__ + r__) / f >= s * .95) 
            goto L270;
        g = 1.0 / f;
        scale[i__] *= f;
        noconv = TRUE;
        
        i__2 = *n;
        for (j = k; j <= i__2; ++j) 
            {
            /* L250: */
            a[i__ + j * a_dim1] *= g;
            }

        i__2 = l;
        for (j = 1; j <= i__2; ++j) 
            {
            /* L260: */
            a[j + i__ * a_dim1] *= f;
            }

        L270:
        ;
        }

    if (noconv) 
        goto L190;

    L280:
    *low = k;
    *igh = l;
    return 0;
    
} 
/* end f2c version of code */
#   endif
    
}


/*---------------------------------------------------------------------------------
|
|   BalBak
|
|   This subroutine forms the eigenvectors of a real general 
|   matrix by back transforming those of the corresponding 
|   balanced matrix determined by  balance. 
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * low and high are integers determined by  balance
|
|    * scale contains information determining the permutations 
|      and scaling factors used by balance
|
|    * m is the number of columns of z to be back transformed
|
|    * z contains the real and imaginary parts of the eigen-
|      vectors to be back transformed in its first m columns
|
|   On output:
|
|    * z contains the real and imaginary parts of the
|      transformed eigenvectors in its first m columns
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|   
---------------------------------------------------------------------------------*/
void BalBak (int dim, int low, int high, MrBFlt *scale, int m, MrBFlt **z)
{
    int         i, j, k, ii;
    MrBFlt      s;

    if (m != 0) /* change "==" to "!=" to eliminate a goto statement */
        {
        if (high != low) /* change "==" to "!=" to eliminate a goto statement */
            {
            for (i=low; i<=high; i++)
                {
                s = scale[i];
                for (j=0; j<m; j++)
                    z[i][j] *= s;
                }
            }
        for (ii=0; ii<dim; ii++)
            {
            i = ii;
            if ((i < low) || (i > high)) /* was (i >= lo) && (i<= hi) but this */
                {                        /* eliminates a goto statement        */
                if (i < low)
                    i = low - ii;
                k = (int)scale[i];
                if (k != i) /* change "==" to "!=" to eliminate a goto statement */
                    {
                    for (j = 0; j < m; j++)
                        {
                        s = z[i][j];
                        z[i][j] = z[k][j];
                        z[k][j] = s;
                        }
                    }
                }
            }
        }

#if 0
/* begin f2c version of code:
   balbak.f -- translated by f2c (version 19971204) */
int balbak (int *nm, int *n, int *low, int *igh, MrBFlt *scale, int *m, MrBFlt *z__)

{

    /* system generated locals */
    int z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static int i__, j, k;
    static MrBFlt s;
    static int ii;

    /* parameter adjustments */
    --scale;
    z_dim1 = *nm;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;

    /* function Body */
    if (*m == 0) 
        goto L200;
    if (*igh == *low) 
        goto L120;

    i__1 = *igh;
    for (i__ = *low; i__ <= i__1; ++i__) 
        {
        s = scale[i__];
        /* .......... left hand eigenvectors are back transformed */
        /*            if the foregoing statement is replaced by */
        /*            s=1.0d0/scale(i) ........... */
        i__2 = *m;
        for (j = 1; j <= i__2; ++j) 
            {
            /* L100: */
            z__[i__ + j * z_dim1] *= s;
            }

        /* L110: */
        }
        
    /* .........for i=low-1 step -1 until 1, igh+1 step 1 until n do -- .......... */
    L120:
    i__1 = *n;
    for (ii = 1; ii <= i__1; ++ii) 
        {
        i__ = ii;
        if (i__ >= *low && i__ <= *igh) 
            goto L140;
    if (i__ < *low) 
        i__ = *low - ii;
    k = (integer) scale[i__];
    if (k == i__) 
        goto L140;

    i__2 = *m;
    for (j = 1; j <= i__2; ++j) 
        {
        s = z__[i__ + j * z_dim1];
        z__[i__ + j * z_dim1] = z__[k + j * z_dim1];
        z__[k + j * z_dim1] = s;
        /* L130: */
        }
    L140:
    ;
    }

    L200:
    return 0;
    
}
/* end f2c version of code */
#endif
        
}


void BetaBreaks (MrBFlt alpha, MrBFlt beta, MrBFlt *values, int K)
{
    int             i;
    MrBFlt          r, quantile, lower, upper;
            
    r = (1.0 / K) * 0.5;
    lower = 0.0;
    upper = (1.0 / K);
    r = (upper - lower) * 0.5 + lower;
    for (i=0; i<K; i++)
        {
        quantile = BetaQuantile (alpha, beta, r);
        values[i] = quantile;
        lower += (1.0/K);
        upper += (1.0/K);
        r += (1.0/K);
        }
        
#   if 0
    for (i=0; i<K; i++)
        {
        MrBayesPrint ("%4d %lf %lf\n", i, values[i]);
        }
#   endif
}


MrBFlt BetaCf (MrBFlt a, MrBFlt b, MrBFlt x)
{
    int         m, m2;
    MrBFlt      aa, c, d, del, h, qab, qam, qap;
    
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;
    d = 1.0 - qab * x / qap;
    if (fabs(d) < (1.0e-30))
        d = (1.0e-30);
    d = 1.0 / d;
    h = d;
    for (m=1; m<=100; m++)
        {
        m2 = 2 * m;
        aa = m * (b-m) * x / ((qam+m2) * (a+m2));
        d = 1.0 + aa * d;
        if (fabs(d) < (1.0e-30))
            d = (1.0e-30);
        c = 1.0 + aa / c;
        if (fabs(c) < (1.0e-30))
            c = (1.0e-30);
        d = 1.0 / d;
        h *= d * c;
        aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2));
        d = 1.0 + aa * d;
        if (fabs(d) < (1.0e-30))
            d = (1.0e-30);
        c = 1.0 + aa / c;
        if (fabs(c) < (1.0e-30))
            c = (1.0e-30);
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < (3.0e-7))
            break;
        }
    if (m > 100)
        {
        MrBayesPrint ("%s   Error in BetaCf.\n", spacer);
        exit(0);
        }
    return (h);
}


MrBFlt BetaQuantile (MrBFlt alpha, MrBFlt beta, MrBFlt x)
{
    int     i, stopIter, direction, nswitches;
    MrBFlt  curPos, curFraction, increment;
    
    i = nswitches = 0;
    curPos = 0.5;
    stopIter = NO;
    increment = 0.25;
    curFraction = IncompleteBetaFunction (alpha, beta, curPos);
    if (curFraction > x)
        direction = DOWN;
    else
        direction = UP;

    while (stopIter == NO)
        {
        curFraction = IncompleteBetaFunction (alpha, beta, curPos);
        if (curFraction > x && direction == DOWN)
            {
            /* continue going down */
            while (curPos - increment <= 0.0)
                {
                increment /= 2.0;
                }
            curPos -= increment;
            }
        else if (curFraction > x && direction == UP)
            {
            /* switch directions, and go down */
            nswitches++;
            direction = DOWN;
            while (curPos - increment <= 0.0)
                {
                increment /= 2.0;
                }
            increment /= 2.0;
            curPos -= increment;
            }
        else if (curFraction < x && direction == UP)
            {
            /* continue going up */
            while (curPos + increment >= 1.0)
                {
                increment /= 2.0;
                }
            curPos += increment;
            }
        else if (curFraction < x && direction == DOWN)
            {
            /* switch directions, and go up */
            nswitches++;
            direction = UP;
            while (curPos + increment >= 1.0)
                {
                increment /= 2.0;
                }
            increment /= 2.0;
            curPos += increment;
            }
        else
            {
            stopIter = YES;
            }
        if (i > 1000 || nswitches > 20)
            stopIter = YES;
        i++;
        }
        
    return (curPos);
}


/*---------------------------------------------------------------------------------
|
|   CalcCijk
|
|   This function precalculates the product of the eigenvectors and their
|   inverse for faster calculation of transition probabilities. The output
|   is a vector of precalculated values. The input is the eigenvectors (u) and
|   the inverse of the eigenvector matrix (v).
|
---------------------------------------------------------------------------------*/
void CalcCijk (int dim, MrBFlt *c_ijk, MrBFlt **u, MrBFlt **v)
{
    register int    i, j, k;
    MrBFlt          *pc;

    pc = c_ijk;
    for (i=0; i<dim; i++)
        for (j=0; j<dim; j++)
            for (k=0; k<dim; k++)
                *pc++ = u[i][k] * v[k][j];
}


/*---------------------------------------------------------------------------------
|
|   CdfBinormal
|
|   F(h1,h2,r) = prob(x<h1, y<h2), where x and y are standard binormal.
|
---------------------------------------------------------------------------------*/
MrBFlt CdfBinormal (MrBFlt h1, MrBFlt h2, MrBFlt r)
{
    return (LBinormal(h1, h2, r) + CdfNormal(h1) + CdfNormal(h2) - 1.0);
}


/*---------------------------------------------------------------------------------
|
|   CdfNormal
|
|   Calculates the cumulative density distribution (CDF) for the normal using:
|
|   Hill, I. D.  1973.  The normal integral.  Applied Statistics, 22:424-427.
|      (AS66)                                                  
|
---------------------------------------------------------------------------------*/
MrBFlt CdfNormal (MrBFlt x)
{
    int             invers = 0;
    MrBFlt          p, limit = 10.0, t = 1.28, y = x*x/2.0;

    if (x < 0.0) 
        {  
        invers = 1;  
        x  *= -1.0; 
        }
    if (x > limit)  
        return (invers ? 0 : 1);
    if (x < t)  
        p = 0.5 - x * (0.398942280444 - 0.399903438504 * y /
            (y + 5.75885480458 - 29.8213557808 /
            (y + 2.62433121679 + 48.6959930692 /
            (y + 5.92885724438))));
    else 
        p = 0.398942280385 * exp(-y) /
            (x - 3.8052e-8 + 1.00000615302 /
            (x + 3.98064794e-4 + 1.98615381364 /
            (x - 0.151679116635 + 5.29330324926 /
            (x + 4.8385912808 - 15.1508972451 /
            (x + 0.742380924027 + 30.789933034 /
            (x + 3.99019417011))))));
            
    return (invers ? p : 1-p);
}


/*---------------------------------------------------------------------------------
|
|   Complex
|
|   Returns a complex number with specified real and imaginary parts.
|
---------------------------------------------------------------------------------*/
complex Complex (MrBFlt a, MrBFlt b)
{
    complex c;
    
    c.re = a;
    c.im = b;
    
    return (c);
}


/*---------------------------------------------------------------------------------
|
|   ComplexAbsoluteValue
|
|   Returns the complex absolute value (modulus) of a complex number.
|
---------------------------------------------------------------------------------*/
MrBFlt ComplexAbsoluteValue (complex a)
{
    MrBFlt      x, y, answer, temp;
    
    x = fabs(a.re);
    y = fabs(a.im);
    if (AreDoublesEqual(x, 0.0, ETA)==YES)  /* x == 0.0 */
        answer = y;
    else if (AreDoublesEqual(y, 0.0, ETA)==YES) /* y == 0.0 */
        answer = x;
    else if (x > y) 
        {
        temp = y / x;
        answer = x * sqrt(1.0 + temp * temp);
        }
    else
        {
        temp = x / y;
        answer = y * sqrt(1.0 + temp * temp);
        }

    return (answer);
}


/*---------------------------------------------------------------------------------
|
|   ComplexAddition
|
|   Returns the complex sum of two complex numbers.
|
---------------------------------------------------------------------------------*/
complex ComplexAddition (complex a, complex b)
{
    complex     c;
    
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    
    return (c);
}


/*---------------------------------------------------------------------------------
|
|   ComplexConjugate
|
|   Returns the complex conjugate of a complex number.
|
---------------------------------------------------------------------------------*/
complex ComplexConjugate (complex a)
{
    complex     c;
    
    c.re = a.re;
    c.im = -a.im;
    
    return (c);
}


/*---------------------------------------------------------------------------------
|
|   ComplexDivision
|
|   Returns the complex quotient of two complex numbers.
|
---------------------------------------------------------------------------------*/
complex ComplexDivision (complex a, complex b)
{
    complex     c;
    MrBFlt      r, den;
    
    if (fabs(b.re) >= fabs(b.im)) 
        {
        r = b.im / b.re;
        den = b.re + r * b.im;
        c.re = (a.re + r * a.im) / den;
        c.im = (a.im - r * a.re) / den;
        } 
    else
        {
        r = b.re / b.im;
        den = b.im + r * b.re;
        c.re = (a.re * r + a.im) / den;
        c.im = (a.im * r - a.re) / den;
        }
    
    return (c);
}


/*---------------------------------------------------------------------------------
|
|   ComplexDivision2
|
|   Returns the complex quotient of two complex numbers. It does not require that
|   the numbers be in a complex structure.
|
---------------------------------------------------------------------------------*/
void ComplexDivision2 (MrBFlt ar, MrBFlt ai, MrBFlt br, MrBFlt bi, MrBFlt *cr, MrBFlt *ci)
{
    MrBFlt      s, ais, bis, ars, brs;

    s = fabs(br) + fabs(bi);
    ars = ar / s;
    ais = ai / s;
    brs = br / s;
    bis = bi / s;
    s = brs*brs + bis*bis;
    *cr = (ars*brs + ais*bis) / s;
    *ci = (ais*brs - ars*bis) / s;
}


/*---------------------------------------------------------------------------------
|
|   ComplexExponentiation
|
|   Returns the complex exponential of a complex number.
|
---------------------------------------------------------------------------------*/
complex ComplexExponentiation (complex a)
{
    complex     c;

    c.re = exp(a.re);
    if (AreDoublesEqual(a.im,0.0, ETA)==YES) /* == 0 */
        c.im = 0; 
    else
        { 
        c.im = c.re*sin(a.im); 
        c.re *= cos(a.im); 
        }

    return (c);
}


/*---------------------------------------------------------------------------------
|
|   ComplexInvertMatrix
|
|   Inverts a matrix of complex numbers using the LU-decomposition method. 
|   The program has the following variables:
|
|      a        -- the matrix to be inverted
|      aInverse -- the results of the matrix inversion
|      dim      -- the dimension of the square matrix a and its inverse
|      dwork    -- a work vector of doubles
|      indx     -- a work vector of integers
|      col      -- carries the results of the back substitution
|      
|   The function returns YES (1) or NO (0) if the results are singular.
|
---------------------------------------------------------------------------------*/
int ComplexInvertMatrix (int dim, complex **a, MrBFlt *dwork, int *indx, complex **aInverse, complex *col)
{
    int             isSingular, i, j;

    isSingular = ComplexLUDecompose (dim, a, dwork, indx, (MrBFlt *)NULL);

    if (isSingular == 0) 
        {
        for (j=0; j<dim; j++) 
            {
            for (i=0; i<dim; i++)
                col[i] = Complex (0.0, 0.0);
            col[j] = Complex (1.0, 0.0);
            ComplexLUBackSubstitution (dim, a, indx, col);
            for (i=0; i<dim; i++)
                aInverse[i][j] = col[i];
            }
        }

    return (isSingular);
}


/*---------------------------------------------------------------------------------
|
|   ComplexExponentiation
|
|   Returns the complex exponential of a complex number.
|
---------------------------------------------------------------------------------*/
complex ComplexLog (complex a)
{
    complex     c;
    
    c.re = log(ComplexAbsoluteValue(a));
    if (AreDoublesEqual(a.re,0.0,ETA)==YES) /* == 0.0 */ 
        {
        c.im = PIOVER2;
        } 
    else 
        {
        c.im = atan2(a.im, a.re);
        }
        
    return (c);
}


/*---------------------------------------------------------------------------------
|
|   ComplexLUBackSubstitution
|
|   Perform back-substitution into a LU-decomposed matrix to obtain
|   the inverse.
|      
---------------------------------------------------------------------------------*/
void ComplexLUBackSubstitution (int dim, complex **a, int *indx, complex *b)
{
    int             i, ip, j, ii = -1;
    complex         sum;

    for (i = 0; i < dim; i++) 
        {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii >= 0) 
            {
            for (j = ii; j <= i - 1; j++)
                sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][j], b[j]));
            } 
        else if (AreDoublesEqual(sum.re,0.0,ETA)==NO || AreDoublesEqual(sum.im, 0.0, ETA)==NO) /* 2x != 0.0 */
            ii = i;
        b[i] = sum;
        }
    for (i = dim - 1; i >= 0; i--) 
        {
        sum = b[i];
        for (j = i + 1; j < dim; j++)
            sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][j], b[j]));
        b[i] = ComplexDivision (sum, a[i][i]);
        }
}


/*---------------------------------------------------------------------------------
|
|   ComplexLUDecompose
|
|   Replaces the matrix a with its LU-decomposition. 
|   The program has the following variables:
|
|      a        -- the matrix
|      dim      -- the dimension of the square matrix a and its inverse
|      vv       -- a work vector of doubles
|      indx     -- row permutation according to partitial pivoting sequence
|      pd       -- 1 if number of row interchanges was even, -1 if number of
|                  row interchanges was odd. Can be NULL.
|      
|   The function returns YES (1) or NO (0) if the results are singular.
|
---------------------------------------------------------------------------------*/
int ComplexLUDecompose (int dim, complex **a, MrBFlt *vv, int *indx, MrBFlt *pd)
{
    int             i, imax, j, k;
    MrBFlt          big, dum, temp, d;
    complex         sum, cdum;

    d = 1.0;
    imax = 0;

    for (i = 0; i < dim; i++) 
        {
        big = 0.0;
        for (j = 0; j < dim; j++) 
            {
            if ((temp = ComplexAbsoluteValue (a[i][j])) > big)
                big = temp;
            }
        if (AreDoublesEqual(big, 0.0, ETA)==YES) /* == 0.0 */
            {
            MrBayesPrint ("%s   Error: Problem in ComplexLUDecompose\n", spacer);
            return (1);
            }
        vv[i] = 1.0 / big;
        }

    for (j = 0; j < dim; j++) 
        {
        for (i = 0; i < j; i++) 
            {
            sum = a[i][j];
            for (k = 0; k < i; k++) 
                sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][k], a[k][j]));
            a[i][j] = sum;
            }
        big = 0.0;
        for (i = j; i < dim; i++) 
            {
            sum = a[i][j];
            for (k = 0; k < j; k++)
            sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][k], a[k][j]));
            a[i][j] = sum;
            dum = vv[i] * ComplexAbsoluteValue (sum);
            if (dum >= big) 
                {
                big = dum;
                imax = i;
                }
            }
        if (j != imax) 
            {
            for (k = 0; k < dim; k++) 
                {
                cdum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = cdum;
                }       
            d = -d;
            vv[imax] = vv[j];
            }
        indx[j] = imax;
        if (AreDoublesEqual(a[j][j].re, 0.0, ETA)==YES && AreDoublesEqual(a[j][j].im, 0.0, ETA)==YES) /* 2x == 0.0 */
            a[j][j] = Complex (1.0e-20, 1.0e-20);
        if (j != dim - 1)
            {
            cdum = ComplexDivision (Complex(1.0, 0.0), a[j][j]);
            for (i = j + 1; i < dim; i++)
            a[i][j] = ComplexMultiplication (a[i][j], cdum);
            }
        }

    if (pd != NULL)
        *pd = d;
        
    return (0);
}


/*---------------------------------------------------------------------------------
|
|   ComplexMultiplication
|
|   Returns the complex product of two complex numbers.
|
---------------------------------------------------------------------------------*/
complex ComplexMultiplication (complex a, complex b)
{
    complex     c;
    
    c.re = a.re * b.re - a.im * b.im;
    c.im = a.im * b.re + a.re * b.im;
    
    return (c);
}


/*---------------------------------------------------------------------------------
|
|   ComplexSquareRoot
|
|   Returns the complex square root of a complex number.
|
---------------------------------------------------------------------------------*/
complex ComplexSquareRoot (complex a)
{
    complex         c;
    MrBFlt          x, y, w, r;
    
    if (AreDoublesEqual(a.re, 0.0, ETA)==YES && AreDoublesEqual(a.im, 0.0, ETA)==YES) /* 2x == 0.0 */
        {
        c.re = 0.0;
        c.im = 0.0;
        return (c);
        }
    else
        {
        x = fabs(a.re);
        y = fabs(a.im);
        if (x >= y)
            {
            r = y / x;
            w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
            }
        else
            {
            r = x / y;
            w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
            }
        if (a.re >= 0.0)
            {
            c.re = w;
            c.im = a.im / (2.0 * w);
            }
        else
            {
            c.im = (a.im >= 0.0) ? w : -w;
            c.re = a.im / (2.0 * c.im);
            }
        return (c);
        }
}


/*---------------------------------------------------------------------------------
|
|   ComplexSubtraction
|
|   Returns the complex difference of two complex numbers.
|
---------------------------------------------------------------------------------*/
complex ComplexSubtraction (complex a, complex b)
{
    complex     c;
    
    c.re = a.re - b.re;
    c.im = a.im - b.im;
    
    return (c);
}


/*---------------------------------------------------------------------------------
|
|   ComputeEigenSystem
|
|   Calculates the eigenvalues, eigenvectors, and the inverse of the eigenvectors
|   for a matrix of real numbers.
|
---------------------------------------------------------------------------------*/
int ComputeEigenSystem (int dim, MrBFlt **a, MrBFlt *v, MrBFlt *vi, MrBFlt **u, int *iwork, MrBFlt *dwork)
{
    int         i, rc;

    rc = EigensForRealMatrix (dim, a, v, vi, u, iwork, dwork);
    if (rc != NO_ERROR)
        {
        MrBayesPrint ("%s   Error in ComputeEigenSystem.\n", spacer);
        return (ERROR);
        }
    for (i=0; i<dim; i++)
        {
        if (AreDoublesEqual(vi[i], 0.0, ETA)==NO) /* != 0.0 */
            return (EVALUATE_COMPLEX_NUMBERS);
        }

    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
|
|   ComputeLandU
|
|   This function computes the L and U decomposition of a matrix. Basically,
|   we find matrices lMat and uMat such that
|
|      lMat * uMat = aMat
|
---------------------------------------------------------------------------------*/
void ComputeLandU (int dim, MrBFlt **aMat, MrBFlt **lMat, MrBFlt **uMat)
{
    int         i, j, k, m, row, col;

    for (j=0; j<dim; j++) 
        {
        for (k=0; k<j; k++)
            for (i=k+1; i<j; i++)
                aMat[i][j] = aMat[i][j] - aMat[i][k] * aMat[k][j];

        for (k=0; k<j; k++)
            for (i=j; i<dim; i++)
                aMat[i][j] = aMat[i][j] - aMat[i][k]*aMat[k][j];

        for (m=j+1; m<dim; m++)
            aMat[m][j] /= aMat[j][j]; 
        }

    for (row=0; row<dim; row++)
        {
        for (col=0; col<dim; col++) 
            {
            if (row <= col) 
                {
                uMat[row][col] = aMat[row][col];
                lMat[row][col] = (row == col ? 1.0 : 0.0);
                }
            else 
                {
                lMat[row][col] = aMat[row][col];
                uMat[row][col] = 0.0;
                }
            }
        }
}


/*---------------------------------------------------------------------------------
|
|   ComputeMatrixExponential
|
|   The method approximates the matrix exponential, f = e^a, using
|   the algorithm 11.3.1, described in:
|  
|   Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
|      The Johns Hopkins University Press, Baltimore, Maryland.
|
|   The method has the advantage of error control. The error is controlled by
|   setting qValue appropriately (using the function SetQValue).
|
---------------------------------------------------------------------------------*/
void ComputeMatrixExponential (int dim, MrBFlt **a, int qValue, MrBFlt **f)
{
    int         i, j, k, negativeFactor;
    MrBFlt      maxAValue, c, **d, **n, **x, **cX;

    d  = AllocateSquareDoubleMatrix (dim);
    n  = AllocateSquareDoubleMatrix (dim);
    x  = AllocateSquareDoubleMatrix (dim);
    cX = AllocateSquareDoubleMatrix (dim);

    SetToIdentity (dim, d);
    SetToIdentity (dim, n);
    SetToIdentity (dim, x);

    maxAValue = 0;
    for (i=0; i<dim; i++)
        maxAValue = MAX (maxAValue, a[i][i]);

    j = MAX (0, LogBase2Plus1 (maxAValue));

    DivideByTwos (dim, a, j);
    
    c = 1;
    for (k=1; k<=qValue; k++) 
        {
        c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

        /* X = AX */
        MultiplyMatrices (dim, a, x, x);

        /* N = N + cX */
        MultiplyMatrixByScalar (dim, x, c, cX);
        AddTwoMatrices (dim, n, cX, n);

        /* D = D + (-1)^k*cX */
        negativeFactor = (k % 2 == 0 ? 1 : -1);
        if (negativeFactor == -1)
            MultiplyMatrixByScalar (dim, cX, negativeFactor, cX);
        AddTwoMatrices (dim, d, cX, d);      
        }

    GaussianElimination (dim, d, n, f);

    for (k = 0; k < j; k++)
        MultiplyMatrices (dim, f, f, f);
    
    for (i=0; i<dim; i++)
        {
        for (j=0; j<dim; j++)
            {
            if (f[i][j] < 0.0)
                f[i][j] = fabs(f[i][j]);
            }
        }
        
    FreeSquareDoubleMatrix (d);
    FreeSquareDoubleMatrix (n);
    FreeSquareDoubleMatrix (x);
    FreeSquareDoubleMatrix (cX);
}


/*---------------------------------------------------------------------------------
|
|   CopyComplexMatrices
|
|   Copies the contents of one matrix of complex numbers to another matrix.
|
---------------------------------------------------------------------------------*/
void CopyComplexMatrices (int dim, complex **from, complex **to)
{
    int         i, j;
    
    for (i=0; i<dim; i++)
        {
        for (j=0; j<dim; j++) 
            {
            to[i][j].re = from[i][j].re;
            to[i][j].im = from[i][j].im;
            }
        }
}


/*---------------------------------------------------------------------------------
|
|   CopyDoubleMatrices
|
|   Copies the contents of one matrix of doubles to another matrix.
|
---------------------------------------------------------------------------------*/
void CopyDoubleMatrices (int dim, MrBFlt **from, MrBFlt **to)
{
    int         i, j;
    
    for (i=0; i<dim; i++)
        {
        for (j=0; j<dim; j++) 
            {
            to[i][j] = from[i][j];
            }
        }
}


/*---------------------------------------------------------------------------------
|
|   DirichletRandomVariable
|
|   Generate a Dirichlet-distributed random variable. The parameter of the
|   Dirichlet is contained in the vector alp. The random variable is contained
|   in the vector z.
|      
---------------------------------------------------------------------------------*/
void DirichletRandomVariable (MrBFlt *alp, MrBFlt *z, int n, RandLong *seed)
{
    int     i;
    MrBFlt  sum;

    sum = 0.0;
    for (i=0; i<n; i++)
        {
        z[i] = RndGamma (alp[i], seed) / 1.0;
        sum += z[i];
        }
    for (i=0; i<n; i++)
        z[i] /= sum;
}


/*---------------------------------------------------------------------------------
|
|   DiscreteGamma
|
|   Discretization of gamma distribution with equal proportions in each
|   category.
|
---------------------------------------------------------------------------------*/
int DiscreteGamma (MrBFlt *rK, MrBFlt alfa, MrBFlt beta, int K, int median)
{
    int             i;
    MrBFlt          gap05 = 1.0/(2.0*K), t, factor = alfa/beta*K, lnga1;

    if (median) 
        {
        for (i=0; i<K; i++) 
            rK[i] = POINTGAMMA((i*2.0+1.0)*gap05, alfa, beta);
        for (i=0,t=0; i<K; i++) 
            t += rK[i];
        for (i=0; i<K; i++)     
            rK[i] *= factor / t;
        }
    else 
        {
        lnga1 = LnGamma(alfa+1);
        /* calculate the points in the gamma distribution */
        for (i=0; i<K-1; i++) 
            rK[i] = POINTGAMMA((i+1.0)/K, alfa, beta);
        /* calculate the cumulative values */
        for (i=0; i<K-1; i++) 
            rK[i] = IncompleteGamma(rK[i] * beta, alfa + 1.0, lnga1);
        rK[K-1] = 1.0;
        /* calculate the relative values and rescale */
        for (i=K-1; i>0; i--)
            {
            rK[i] -= rK[i-1];
            rK[i] *= factor;
            }
        rK[0] *= factor;
        }

    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
 |
 |   DiscreteLogNormal
 |
 |   Discretization of lognormal distribution with equal proportions in each
 |   category.
 |
 |   LBH Notes:     K = # of rate classes
 |                *rK = pointer to output rate class matrix
 |               alfa = alpha param
 |               beta = beta param
 |             median = flag to use media or not (1 = use median, 0 = mean?)
 |
 ---------------------------------------------------------------------------------*/
int DiscreteLogNormal (MrBFlt *rK, MrBFlt sigma, int K, int median)
{
    int i;
    MrBFlt t, factor;
    MrBFlt sigmaL = sqrt(sigma);
    MrBFlt mu = -0.5*sigmaL*sigmaL;
    if (median)
        {
        for (i=0; i<K; i++) {
            rK[i] = QuantileLogNormal( ((2.0*i + 1) / (2.0 * K)), mu, sigmaL);
            }
        for (i=0,t=0.0; i<K; i++) {
            t = t+rK[i];
            }
        t /= K;
        for (i=0; i<K; i++)
            rK[i] /= t;
        }
    else
        {
        mu = -0.5*sigmaL*sigmaL;
        /* Mean set to 1.0 so factor = K */
        factor = 1.0*K;
        for (i=0; i<K-1; i++) {
            rK[i] = QuantileLogNormal(((i + 1.0) / (K)), mu, sigmaL);
            }
        for (i=0; i<K-1; i++) {
            //rK[i] = LogNormalPoint(rK[i], mu, sigma);
            //rK[i] = QuantileLogNormal(rK[i], mu, sigma);
            //rK[i] = CdfNormal((log(rK[i])-mu)/sigma);
            rK[i] = 1 - (1.0 * CdfNormal((mu + sigmaL*sigmaL - log(rK[i]))/sigmaL));
            }
        rK[K-1] = 1.0;
        for (i=K-1; i>0; i--) {
            rK[i] -= rK[i-1];
            rK[i] *= factor;
            }
        rK[0] *= factor;
        }

    return (NO_ERROR);
}


/* LogNormal Quantile Function */
MrBFlt QuantileLogNormal (MrBFlt prob, MrBFlt mu, MrBFlt sigma)
{
    MrBFlt a = 0.0, b = 0.0;
    a = PointNormal((0.5*(2.0*prob-1.0))+0.5) / sqrt(2.0);
    b = mu+(sqrt(2.0)* sigma * a);
    return exp(b);
}


/* LogNormal Point Function */
MrBFlt LogNormalPoint (MrBFlt x, MrBFlt mu, MrBFlt sigma)
{
    if(x <= 0.0) return(0.0);
    MrBFlt a = LnProbLogNormal(mu, sigma, x);
    return exp(a);
}


/*---------------------------------------------------------------------------------
|
|   DivideByTwos
|
|   Divides all of the elements of the matrix a by 2^power.
|      
---------------------------------------------------------------------------------*/
void DivideByTwos (int dim, MrBFlt **a, int power)
{
    int         divisor = 1, i, row, col;

    for (i=0; i<power; i++)
        divisor = divisor * 2;

    for (row=0; row<dim; row++)
        for (col=0; col<dim; col++)
            a[row][col] /= divisor;
}


/*---------------------------------------------------------------------------------
|
|   D_sign
|
|   This function is called from "Hqr2".
|
---------------------------------------------------------------------------------*/
MrBFlt D_sign (MrBFlt a, MrBFlt b)
{
    MrBFlt      x;

    x = (a >= 0 ? a : -a);
    
    return (b >= 0 ? x : -x);
}


/*---------------------------------------------------------------------------------
|
|   Eigens
|
|   The matrix of interest is a. The ouptut is the real and imaginary parts of the 
|   eigenvalues (wr and wi). z contains the real and imaginary parts of the 
|   eigenvectors. iv2 and fv1 are working vectors.
|      
---------------------------------------------------------------------------------*/
int EigensForRealMatrix (int dim, MrBFlt **a, MrBFlt *wr, MrBFlt *wi, MrBFlt **z, int *iv1, MrBFlt *fv1)
{
    static int  is1, is2;
    int         ierr;

    Balanc (dim, a, &is1, &is2, fv1);
    ElmHes (dim, is1, is2, a, iv1);
    ElTran (dim, is1, is2, a, iv1, z);
    ierr = Hqr2 (dim, is1, is2, a, wr, wi, z);
    if (ierr == 0)
        BalBak (dim, is1, is2, fv1, dim, z);

    return (ierr);
}


/*---------------------------------------------------------------------------------
|
|   ElmHes
|
|   Given a real general matrix, this subroutine
|   reduces a submatrix situated in rows and columns
|   low through high to upper Hessenberg form by
|   stabilized elementary similarity transformations.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc.  if  balanc  has not been used,
|      set low=1, high=dim.
|
|    * a contains the input matrix.
|
|   On output:
|
|    * a contains the hessenberg matrix.  The multipliers
|      which were used in the reduction are stored in the
|      remaining triangle under the hessenberg matrix.
|
|    * interchanged contains information on the rows and columns
|      interchanged in the reduction.
|
|   Only elements low through high are used.
|
---------------------------------------------------------------------------------*/
void ElmHes (int dim, int low, int high, MrBFlt **a, int *interchanged)
{
    int         i, j, m, la, mm1, kp1, mp1;
    MrBFlt      x, y;
    
    la = high - 1;
    kp1 = low + 1;
    if (la < kp1)
        return; /* remove goto statement, which exits at bottom of function */

    for (m=kp1; m<=la; m++)
        {
        mm1 = m - 1;
        x = 0.0;
        i = m;
    
        for (j=m; j<=high; j++)
            {
            if (fabs(a[j][mm1]) > fabs(x)) /* change direction of inequality */
                {                          /* remove goto statement          */
                x = a[j][mm1];
                i = j;
                }
            }
    
        interchanged[m] = i;
        if (i != m) /* change "==" to "!=", eliminating goto statement */
            {
            /* interchange rows and columns of a */
            for (j=mm1; j<dim; j++)
                {
                y = a[i][j];
                a[i][j] = a[m][j];
                a[m][j] = y;
                }
            for (j=0; j<=high; j++)
                {
                y = a[j][i];
                a[j][i] = a[j][m];
                a[j][m] = y;
                }
            }

        if (AreDoublesEqual(x, 0.0, ETA)==NO) /* change "==" to "!=", eliminating goto statement */
            {
            mp1 = m + 1;
        
            for (i=mp1; i<=high; i++)
                {
                y = a[i][mm1];
                if (AreDoublesEqual(y, 0.0, ETA)==NO) /* != 0.0 */
                    {
                    y /= x;
                    a[i][mm1] = y;
                    for (j = m; j < dim; j++)
                        a[i][j] -= y * a[m][j];
                    for (j = 0; j <= high; j++)
                        a[j][m] += y * a[j][i];
                    }
                }
            }
        }

#if 0
/* begin f2c version of code:
   elmhes.f -- translated by f2c (version 19971204) */
int elmhes (int *nm, int *n, int *low, int *igh, MrBFlt *a, int *int__)

{

    /*system generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    MrBFlt d__1;

    /* local variables */
    static int i__, j, m;
    static MrBFlt x, y;
    static int la, mm1, kp1, mp1;

    /* parameter adjustments */
    a_dim1 = *nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --int__;

    /* function body */
    la = *igh - 1;
    kp1 = *low + 1;
    if (la < kp1) 
        goto L200;

    i__1 = la;
    for (m = kp1; m <= i__1; ++m) 
        {
        mm1 = m - 1;
        x = 0.;
        i__ = m;
        i__2 = *igh;
        for (j = m; j <= i__2; ++j) 
            {
            if ((d__1 = a[j + mm1 * a_dim1], abs(d__1)) <= abs(x)) 
                goto L100;
            x = a[j + mm1 * a_dim1];
            i__ = j;
            L100:
            ;
        }

    int__[m] = i__;
    if (i__ == m) 
        goto L130;

    /* .......... interchange rows and columns of a.......... */
    i__2 = *n;
    for (j = mm1; j <= i__2; ++j) 
        {
        y = a[i__ + j * a_dim1];
        a[i__ + j * a_dim1] = a[m + j * a_dim1];
        a[m + j * a_dim1] = y;
        /* L110: */
        }

    i__2 = *igh;
    for (j = 1; j <= i__2; ++j) 
        {
        y = a[j + i__ * a_dim1];
        a[j + i__ * a_dim1] = a[j + m * a_dim1];
        a[j + m * a_dim1] = y;
        /* L120: */
        }
        
    /* .......... end interchange .......... */
    L130:
    if (x == 0.) 
        goto L180;
    mp1 = m + 1;

    i__2 = *igh;
    for (i__ = mp1; i__ <= i__2; ++i__) 
        {
        y = a[i__ + mm1 * a_dim1];
        if (y == 0.) 
            goto L160;
        y /= x;
        a[i__ + mm1 * a_dim1] = y;

        i__3 = *n;
        for (j = m; j <= i__3; ++j) 
            {
            /* L140: */
            a[i__ + j * a_dim1] -= y * a[m + j * a_dim1];
            }

        i__3 = *igh;
        for (j = 1; j <= i__3; ++j) 
            {
            /* L150: */
            a[j + m * a_dim1] += y * a[j + i__ * a_dim1];
            }

        L160:
            ;
        }

    L180:
        ;
    }

    L200:
    return 0;
    
}
/* end f2c version of code */
#endif
        
}


/*---------------------------------------------------------------------------------
|
|   ElTran
|
|   This subroutine accumulates the stabilized elementary
|   similarity transformations used in the reduction of a
|   real general matrix to upper Hessenberg form by ElmHes.
|
|   On input:
|
|    * dim is the order of the matrix.
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc. If Balanc has not been used,
|      set low=0, high=dim-1.
|
|    * a contains the multipliers which were used in the
|      reduction by  ElmHes in its lower triangle
|      below the subdiagonal.
|
|    * interchanged contains information on the rows and columns
|      interchanged in the reduction by ElmHes.
|      only elements low through high are used.
|
|   On output:
|
|    * z contains the transformation matrix produced in the
|      reduction by ElmHes.
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|   
---------------------------------------------------------------------------------*/
void ElTran (int dim, int low, int high, MrBFlt **a, int *interchanged, MrBFlt **z)
{
    int         i, j, mp;

    /* initialize z to identity matrix */
    for (j=0; j<dim; j++)
        {
        for (i=0; i<dim; i++)
            z[i][j] = 0.0;
        z[j][j] = 1.0;
        }
    for (mp=high-1; mp>=low+1; mp--) /* there were a number of additional    */
        {                            /* variables (kl, la, m, mm, mp1) that  */
        for (i=mp+1; i<=high; i++)   /* have been eliminated here simply by  */
            z[i][mp] = a[i][mp-1];   /* initializing variables appropriately */
        i = interchanged[mp];        /* in the loops                         */
        if (i != mp) /* change "==" to "!=" to eliminate a goto statement */
            {
            for (j=mp; j<=high; j++)
                {
                z[mp][j] = z[i][j];
                z[i][j] = 0.0;
                }
            z[i][mp] = 1.0;
            }
        }
    
#if 0
/* begin f2c version of code:
   eltran.f -- translated by f2c (version 19971204) */
int eltran (int *nm, int *n, int *low, int *igh, MrBFlt *a, int *int__, MrBFlt *z__)

{

    /* system generated locals */
    int a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;

    /* local variables */
    static int i__, j, kl, mm, mp, mp1;

    /*     .......... initialize z to identity matrix .......... */
    
    /* parameter adjustments */
    z_dim1 = *nm;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --int__;
    a_dim1 = *nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) 
        {
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) 
            {
            /* L60: */
            z__[i__ + j * z_dim1] = 0.0;
            }
        z__[j + j * z_dim1] = 1.0;
        /* L80: */
        }

    kl = *igh - *low - 1;
    if (kl < 1) 
        goto L200;

    /* .......... for mp=igh-1 step -1 until low+1 do -- .......... */
    i__1 = kl;
    for (mm = 1; mm <= i__1; ++mm) 
        {
        mp = *igh - mm;
        mp1 = mp + 1;
        i__2 = *igh;
        for (i__ = mp1; i__ <= i__2; ++i__) 
            {
            /* L100: */
            z__[i__ + mp * z_dim1] = a[i__ + (mp - 1) * a_dim1];
            }
        i__ = int__[mp];
        if (i__ == mp) 
            goto L140;
        i__2 = *igh;
        for (j = mp; j <= i__2; ++j) 
            {
            z__[mp + j * z_dim1] = z__[i__ + j * z_dim1];
            z__[i__ + j * z_dim1] = 0.;
            /* L130: */
            }
        z__[i__ + mp * z_dim1] = 1.;
        L140:
            ;
        }

    L200:
    return 0;

}
/* end f2c version of code */
#endif
    
}


/*---------------------------------------------------------------------------------
|
|   Exchange
|
---------------------------------------------------------------------------------*/
void Exchange (int j, int k, int l, int m, int n, MrBFlt **a, MrBFlt *scale)
{
    int         i;
    MrBFlt      f;

    scale[m] = (MrBFlt)j;
    if (j != m)
        {
        for (i = 0; i <= l; i++)
            {
            f = a[i][j];
            a[i][j] = a[i][m];
            a[i][m] = f;
            }   
        for (i = k; i < n; i++)
            {
            f = a[j][i];
            a[j][i] = a[m][i];
            a[m][i] = f;
            }
        }
}


/*---------------------------------------------------------------------------------
|
|   Factorial
|
|   Returns x!
|      
---------------------------------------------------------------------------------*/
MrBFlt Factorial (int x)
{
    int         i;
    MrBFlt      fac;
    
    fac = 1.0;
    for (i=0; i<x; i++)
        {
        fac *= (i+1);
        }
        
    return (fac);
}


/*---------------------------------------------------------------------------------
|
|   ForwardSubstitutionRow
|
---------------------------------------------------------------------------------*/
void ForwardSubstitutionRow (int dim, MrBFlt **L, MrBFlt *b)
{
    int         i, j;
    MrBFlt      dotProduct;

    b[0] = b[0] / L[0][0];
    for (i=1; i<dim; i++) 
        {
        dotProduct = 0.0;
        for (j=0; j<i; j++)
            dotProduct += L[i][j] * b[j];
        b[i] = (b[i] - dotProduct) / L[i][i];
        }
}


/*---------------------------------------------------------------------------------
|
|   FreeSquareComplexMatrix
|
|   Frees a matrix of complex numbers.
|      
---------------------------------------------------------------------------------*/
void FreeSquareComplexMatrix (complex **m)
{
    free((char *) (m[0]));
    free((char *) (m));
}


/*---------------------------------------------------------------------------------
|
|   FreeSquareDoubleMatrix
|
|   Frees a matrix of doubles.
|      
---------------------------------------------------------------------------------*/
void FreeSquareDoubleMatrix (MrBFlt **m)
{
    free((char *) (m[0]));
    free((char *) (m));
}


/*---------------------------------------------------------------------------------
|
|   FreeSquareIntegerMatrix
|
|   Frees a matrix of integers.
|      
---------------------------------------------------------------------------------*/
void FreeSquareIntegerMatrix (int **m)
{
    free((char *) (m[0]));
    free((char *) (m));
}


/*---------------------------------------------------------------------------------
|
|   GammaRandomVariable
|
|   This function generates a gamma-distributed random variable with parameters
|   a and b. The mean is E(X) = a / b and the variance is Var(X) = a / b^2.
|      
---------------------------------------------------------------------------------*/
MrBFlt GammaRandomVariable (MrBFlt a, MrBFlt b, RandLong *seed)
{
    return (RndGamma (a, seed) / b);
}


/*---------------------------------------------------------------------------------
|
|   GaussianElimination
|      
---------------------------------------------------------------------------------*/
void GaussianElimination (int dim, MrBFlt **a, MrBFlt **bMat, MrBFlt **xMat)
{
    int         i, k;
    MrBFlt      *bVec, **lMat, **uMat;

    lMat = AllocateSquareDoubleMatrix (dim);
    uMat = AllocateSquareDoubleMatrix (dim);
    bVec = (MrBFlt *) SafeMalloc ((size_t)dim * sizeof(MrBFlt));
    if (!bVec)
        {
        MrBayesPrint ("%s   Error: Problem allocating bVec\n", spacer);
        exit (0);
        }

    ComputeLandU (dim, a, lMat, uMat);

    for (k=0; k<dim; k++) 
        {
        
        for (i=0; i<dim; i++)
            bVec[i] = bMat[i][k];

        /* Answer of Ly = b (which is solving for y) is copied into b. */
        ForwardSubstitutionRow (dim, lMat, bVec);

        /* Answer of Ux = y (solving for x and the y was copied into b above) 
           is also copied into b. */
        BackSubstitutionRow (dim, uMat, bVec);

        for (i=0; i<dim; i++)
            xMat[i][k] = bVec[i];

        }
    
    FreeSquareDoubleMatrix (lMat);
    FreeSquareDoubleMatrix (uMat);
    free (bVec);
}


/*---------------------------------------------------------------------------------
|
|   GetEigens
|
|   returns NO if non complex eigendecomposition, YES if complex eigendecomposition,  ABORT if an error has occured
|
---------------------------------------------------------------------------------*/
int GetEigens (int dim, MrBFlt **q, MrBFlt *eigenValues, MrBFlt *eigvalsImag, MrBFlt **eigvecs, MrBFlt **inverseEigvecs, complex **Ceigvecs, complex **CinverseEigvecs)
{
    int         i, j, rc, *iWork, isComplex;
    MrBFlt      **tempWork, *dWork;
    complex     **cWork, *Ccol;

    /* allocate memory */
    dWork = (MrBFlt *) SafeMalloc ((size_t)dim * sizeof(MrBFlt));
    iWork = (int *) SafeMalloc ((size_t)dim * sizeof(int));
    if (!dWork || !iWork)
        {
        MrBayesPrint ("%s   Error: Problem in GetEigens\n", spacer);
        exit (0);
        }

    /* calculate eigenvalues and eigenvectors */
    isComplex = NO;
    rc = ComputeEigenSystem (dim, q, eigenValues, eigvalsImag, eigvecs, iWork, dWork);
    if (rc != NO_ERROR)
        {
        if (rc == EVALUATE_COMPLEX_NUMBERS)
            isComplex = YES;
        else
            isComplex = ABORT;
        }

    /* invert eigenvectors */
    if (isComplex == NO)
        {
        tempWork = AllocateSquareDoubleMatrix (dim);
        CopyDoubleMatrices (dim, eigvecs, tempWork);
        InvertMatrix (dim, tempWork, dWork, iWork, inverseEigvecs);
        FreeSquareDoubleMatrix (tempWork);
        }
    else if (isComplex == YES)
        {
        for (i=0; i<dim; i++)
            {
              if (fabs(eigvalsImag[i])<1E-20) /* == 0.0 */
                { 
                for (j=0; j<dim; j++)
                    {
                    Ceigvecs[j][i].re = eigvecs[j][i];
                    Ceigvecs[j][i].im = 0.0;
                    }
                }
            else if (eigvalsImag[i] > 0)
                { 
                for (j=0; j<dim; j++)
                    {
                    Ceigvecs[j][i].re = eigvecs[j][i];
                    Ceigvecs[j][i].im = eigvecs[j][i + 1];
                    }
                }
            else if (eigvalsImag[i] < 0)
                { 
                for (j=0; j<dim; j++)
                    {
                    Ceigvecs[j][i].re =  eigvecs[j][i-1];
                    Ceigvecs[j][i].im = -eigvecs[j][i];
                    }
                }
            }
        Ccol = (complex *) SafeMalloc ((size_t)dim * sizeof(complex));
        if (!Ccol)
            {
            MrBayesPrint ("%s   Error: Problem in GetEigens\n", spacer);
            exit (0);
            }
        cWork = AllocateSquareComplexMatrix (dim);
        CopyComplexMatrices (dim, Ceigvecs, cWork);
        ComplexInvertMatrix (dim, cWork, dWork, iWork, CinverseEigvecs, Ccol);
        free (Ccol);
        FreeSquareComplexMatrix (cWork);
        }

    free (dWork);
    free (iWork);

    return (isComplex);
}


/*---------------------------------------------------------------------------------
|
|   Hqr2
|
|   This subroutine finds the eigenvalues and eigenvectors
|   of a real upper Hessenberg matrix by the QR method. The
|   eigenvectors of a real general matrix can also be found
|   if ElmHes  and ElTran or OrtHes and OrTran have
|   been used to reduce this general matrix to Hessenberg form
|   and to accumulate the similarity transformations.
|
|   On input:
|
|    * dim is the order of the matrix.
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc. If  balanc has not been used,
|      set low=0, high=dim-1.
|
|    * h contains the upper hessenberg matrix. Information about
|      the transformations used in the reduction to Hessenberg
|      form by  ElmHes  or OrtHes, if performed, is stored
|      in the remaining triangle under the Hessenberg matrix.
|
|   On output:
|
|    * h has been destroyed.
|
|    * wr and wi contain the real and imaginary parts,
|      respectively, of the eigenvalues. The eigenvalues
|      are unordered except that complex conjugate pairs
|      of values appear consecutively with the eigenvalue
|      having the positive imaginary part first. If an
|      error exit is made, the eigenvalues should be correct
|      for indices j,...,dim-1.
|
|    * z contains the transformation matrix produced by ElTran
|      after the reduction by ElmHes, or by OrTran after the
|      reduction by OrtHes, if performed. If the eigenvectors
|      of the Hessenberg matrix are desired, z must contain the
|      identity matrix.
|
|   Calls ComplexDivision2 for complex division.
|
|   This function returns:
|      zero       for normal return,
|      j          if the limit of 30*n iterations is exhausted
|                 while the j-th eigenvalue is being sought.
|
|   This subroutine is a translation of the ALGOL procedure HQR2,
|   Num. Math. 14, 219,231(1970) by Martin, Peters, and Wilkinson.
|   Handbook for Automatic Computation, vol. II - Linear Algebra,
|   pp. 357-391 (1971).
|   
---------------------------------------------------------------------------------*/
int Hqr2 (int dim, int low, int high, MrBFlt **h, MrBFlt *wr, MrBFlt *wi, MrBFlt **z)
{
    int         i, j, k, l, m, na, en, notlas, mp2, itn, its, enm2, twoRoots;
    MrBFlt      norm, p=0.0, q=0.0, r=0.0, s=0.0, t, w=0.0, x, y=0.0, ra, sa, vi, vr, zz=0.0, tst1, tst2;

    norm = 0.0;
    k = 0;  /* used for array indexing. FORTRAN version: k = 1 */
    
    /* store roots isolated by balance, and compute matrix norm */
    for (i=0; i<dim; i++)
        {
        for (j=k; j<dim; j++)
            norm += fabs(h[i][j]);

        k = i;
        if ((i < low) || (i > high))
            {
            wr[i] = h[i][i];
            wi[i] = 0.0;
            }
        }
    en = high;
    t = 0.0;
    itn = dim * 30;

    /* search for next eigenvalues */
    while (en >= low) /* changed from an "if (en < lo)" to eliminate a goto statement */
        {
        its = 0;
        na = en - 1;
        enm2 = na - 1;
        twoRoots = FALSE;

        for (;;)
            {
            for (l=en; l>low; l--) /* changed indexing, got rid of lo, ll */
                {
                s = fabs(h[l-1][l-1]) + fabs(h[l][l]);
                if (AreDoublesEqual(s, 0.0, ETA)==YES) /* == 0.0 */
                    s = norm;
                tst1 = s;
                tst2 = tst1 + fabs(h[l][l-1]);
                if (fabs(tst2 - tst1) < ETA) /* tst2 == tst1 */
                    break; /* changed to break to remove a goto statement */
                }
    
            /* form shift */
            x = h[en][en];
            if (l == en) /* changed to break to remove a goto statement */
                break;
            y = h[na][na];
            w = h[en][na] * h[na][en];
            if (l == na)         /* used to return to other parts of the code */
                {
                twoRoots = TRUE;
                break;
                }
            if (itn == 0)
                return (en);
                
            /* form exceptional shift */
            if ((its == 10) || (its == 20)) /* changed to remove a goto statement */
                {
                t += x;
                for (i = low; i <= en; i++)
                    h[i][i] -= x;
                s = fabs(h[en][na]) + fabs(h[na][enm2]);
                x = 0.75 * s;
                y = x;
                w = -0.4375 * s * s;
                }
            its++;
            itn--;
            
            /* look for two consecutive small sub-diagonal elements */
            for (m=enm2; m>=l; m--)
                {
                /* removed m = enm2 + l - mm and above loop to remove variables */
                zz = h[m][m];
                r = x - zz;
                s = y - zz;
                p = (r * s - w) / h[m+1][m] + h[m][m+1];
                q = h[m+1][m+1] - zz - r - s;
                r = h[m+2][m+1];
                s = fabs(p) + fabs(q) + fabs(r);
                p /= s;
                q /= s;
                r /= s;
                if (m == l)
                    break; /* changed to break to remove a goto statement */
                tst1 = fabs(p) * (fabs(h[m-1][m-1]) + fabs(zz) + fabs(h[m+1][m+1]));
                tst2 = tst1 + fabs(h[m][m-1]) * (fabs(q) + fabs(r));
                if (fabs(tst2 - tst1) < ETA) /* tst2 == tst1 */
                    break; /* changed to break to remove a goto statement */
                }
        
            mp2 = m + 2;
            for (i = mp2; i <= en; i++)
                {
                h[i][i-2] = 0.0;
                if (i != mp2) /* changed "==" to "!=" to remove a goto statement */
                    h[i][i-3] = 0.0;
                }
    
            /* MrBFlt QR step involving rows l to en and columns m to en */
            for (k=m; k<=na; k++)
                {
                notlas = (k != na);
                if (k != m) /* changed "==" to "!=" to remove a goto statement */
                    {
                    p = h[k][k-1];
                    q = h[k+1][k-1];
                    r = 0.0;
                    if (notlas)
                        r = h[k+2][k-1];
                    x = fabs(p) + fabs(q) + fabs(r);
                    if (x < ETA) /* == 0.0 */
                        continue; /* changed to continue remove a goto statement */
                    p /= x;
                    q /= x;
                    r /= x;
                    }
    
                /*s = sqrt(p*p+q*q+r*r);
                sgn = (p<0)?-1:(p>0);
                s = sgn*sqrt(p*p+q*q+r*r);*/
                s = D_sign(sqrt(p*p + q*q + r*r), p);
                if (k != m) /* changed "==" to "!=" to remove a goto statement */
                    h[k][k-1] = -s * x;
                else if (l != m) /* else if gets rid of another goto statement */
                    h[k][k-1] = -h[k][k-1];
                p += s;
                x = p / s;
                y = q / s;
                zz = r / s;
                q /= p;
                r /= p;
                if (!notlas) /* changed to !notlas to remove goto statement (see **) */
                    {
                    /* row modification */
                    for (j=k; j<dim; j++)
                        {
                        p = h[k][j] + q * h[k+1][j];
                        h[k][j] -= p * x;
                        h[k+1][j] -= p * y;
                        } 
                    j = MIN(en, k + 3);
                    
                    /* column modification */
                    for (i=0; i<=j; i++)
                        {
                        p = x * h[i][k] + y * h[i][k+1];
                        h[i][k] -= p;
                        h[i][k+1] -= p * q;
                        }
                        
                    /* accumulate transformations */
                    for (i=low; i<=high; i++)
                        {
                        p = x * z[i][k] + y * z[i][k+1];
                        z[i][k] -= p;
                        z[i][k+1] -= p * q;
                        }
                    }
                else /* (**) also put in else */
                    {
                    /* row modification */
                    for (j=k; j<dim; j++)
                        {
                        p = h[k][j] + q * h[k+1][j] + r * h[k+2][j];
                        h[k][j] -= p * x;
                        h[k+1][j] -= p * y;
                        h[k+2][j] -= p * zz;
                        }
                    j = MIN(en, k + 3);
                    
                    /* column modification */
                    for (i = 0; i <= j; i++)
                        {
                        p = x * h[i][k] + y * h[i][k+1] + zz * h[i][k+2];
                        h[i][k] -= p;
                        h[i][k+1] -= p * q;
                        h[i][k+2] -= p * r;
                        }
                        
                    /* accumulate transformations */
                    for (i = low; i <= high; i++)
                        {
                        p = x * z[i][k] + y * z[i][k+1] + zz * z[i][k+2];
                        z[i][k] -= p;
                        z[i][k+1] -= p * q;
                        z[i][k+2] -= p * r;
                        }
                    }
                }
            }

        if (twoRoots)
            {
            /* two roots found */
            p = (y - x) / 2.0;
            q = p * p + w;
            zz = sqrt(fabs(q));
            h[en][en] = x + t;
            x = h[en][en];
            h[na][na] = y + t;
            if (q >= -1e-12) /* change "<" to ">=", and also change "0.0" to */
                {            /* a small number (Swofford's change)           */
                /* real pair */
                zz = p + D_sign(zz, p);
                wr[na] = x + zz;
                wr[en] = wr[na];
                if (fabs(zz) > ETA) /* != 0.0 */
                    wr[en] = x - w/zz;
                wi[na] = 0.0;
                wi[en] = 0.0;
                x = h[en][na];
                s = fabs(x) + fabs(zz);
                p = x / s;
                q = zz / s;
                r = sqrt(p*p + q*q);
                p /= r;
                q /= r;
                
                /* row modification */
                for (j=na; j<dim; j++)
                    {
                    zz = h[na][j];
                    h[na][j] = q * zz + p * h[en][j];
                    h[en][j] = q * h[en][j] - p * zz;
                    }
                    
                /* column modification */
                for (i = 0; i <= en; i++)
                    {
                    zz = h[i][na];
                    h[i][na] = q * zz + p * h[i][en];
                    h[i][en] = q * h[i][en] - p * zz;
                    }
                    
                /* accumulate transformations */
                for (i = low; i <= high; i++)
                    {
                    zz = z[i][na];
                    z[i][na] = q * zz + p * z[i][en];
                    z[i][en] = q * z[i][en] - p * zz;
                    }
                }
            else
                {
                /* complex pair */
                wr[na] = x + p;
                wr[en] = x + p;
                wi[na] = zz;
                wi[en] = -zz;
                }
            en = enm2;
            }
        else
            {
            /* one root found */
            h[en][en] = x + t;
            wr[en] = h[en][en];
            wi[en] = 0.0;
            en = na;
            }
        }
    
    if (fabs(norm) < ETA) /* == 0.0 */
        return (0); /* was a goto end of function */

    for (en=dim-1; en>=0; en--)
        {
        /*en = n - nn - 1; and change for loop */
        p = wr[en];
        q = wi[en];
        na = en - 1;

        if (q < -1e-12)
            {
            /* last vector component chosen imaginary so that eigenvector
               matrix is triangular */
            m = na;
            if (fabs(h[en][na]) > fabs(h[na][en]))
                {
                h[na][na] = q / h[en][na];
                h[na][en] = -(h[en][en] - p) / h[en][na];
                }
            else
                ComplexDivision2 (0.0, -h[na][en], h[na][na] - p, q, &h[na][na], &h[na][en]);

            h[en][na] = 0.0;
            h[en][en] = 1.0;
            enm2 = na - 1;
            if (enm2 >= 0) /* changed direction to remove goto statement */
                {
                for (i=enm2; i>=0; i--)
                    {
                    w = h[i][i] - p;
                    ra = 0.0;
                    sa = 0.0;
            
                    for (j=m; j<=en; j++)
                        {
                        ra += h[i][j] * h[j][na];
                        sa += h[i][j] * h[j][en];
                        }
            
                    if (wi[i] < 0.0) /* changed direction to remove goto statement */
                        {
                        zz = w;
                        r = ra;
                        s = sa;
                        }
                    else
                        {
                        m = i;
                        if (fabs(wi[i])<ETA) /* == 0.0 */ /* changed direction to remove goto statement */
                            ComplexDivision2 (-ra, -sa, w, q, &h[i][na], &h[i][en]);
                        else
                            {
                            /* solve complex equations */
                            x = h[i][i+1];
                            y = h[i+1][i];
                            vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
                            vi = (wr[i] - p) * 2.0 * q;
                            if ((fabs(vr)<ETA) && (fabs(vi)<ETA))
                                {
                                tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
                                vr = tst1;
                                do  {
                                    vr *= .01;
                                    tst2 = tst1 + vr;
                                    }
                                    while (tst2 > tst1); /* made into a do/while loop */
                                }
                            ComplexDivision2 (x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, &h[i][na], &h[i][en]);
                            if (fabs(x) > fabs(zz) + fabs(q)) /* changed direction to remove goto statement */
                                {
                                h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;
                                h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;
                                }
                            else
                                ComplexDivision2 (-r - y * h[i][na], -s - y * h[i][en], zz, q, &h[i+1][na], &h[i+1][en]);
                            }
                            
                        /* overflow control */
                        tst1 = fabs(h[i][na]);
                        tst2 = fabs(h[i][en]);
                        t = MAX(tst1, tst2);
                        if (t > ETA) /* t != 0.0 */
                            {
                            tst1 = t;
                            tst2 = tst1 + 1.0 / tst1;
                            if (tst2 <= tst1)
                                {
                                for (j = i; j <= en; j++)
                                    {
                                    h[j][na] /= t;
                                    h[j][en] /= t;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        else if (fabs(q)<ETA)
            {
            /* real vector */
            m = en;
            h[en][en] = 1.0;
            if (na >= 0)
                {
                for (i=na; i>=0; i--)
                    {
                    w = h[i][i] - p;
                    r = 0.0;
                    for (j = m; j <= en; j++)
                        r += h[i][j] * h[j][en];
                    if (wi[i] < 0.0) /* changed direction to remove goto statement */
                        {
                        zz = w;
                        s = r;
                        continue;  /* changed to continue to remove goto statement */
                        }
                    else
                        {
                        m = i;
                        if (fabs(wi[i])<ETA) /* changed to remove goto statement */
                            {
                            t = w;
                            if (fabs(t)<ETA)  /* changed to remove goto statement */
                                {
                                tst1 = norm;
                                t = tst1;
                                do  {
                                    t *= .01;
                                    tst2 = norm + t;
                                    }
                                    while (tst2 > tst1);
                                }           
                            h[i][en] = -r / t;
                            }
                        else
                            {
                            /* solve real equations */
                            x = h[i][i+1];
                            y = h[i+1][i];
                            q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
                            t = (x * s - zz * r) / q;
                            h[i][en] = t;
                            if (fabs(x) > fabs(zz))  /* changed direction to remove goto statement */
                                h[i+1][en] = (-r - w * t) / x;
                            else
                                h[i+1][en] = (-s - y * t) / zz;
                            }
                
                        /* overflow control */
                        t = fabs(h[i][en]);
                        if (t > ETA)
                            {
                            tst1 = t;
                            tst2 = tst1 + 1. / tst1;
                            if (tst2 <= tst1)
                                {
                                for (j = i; j <= en; j++)
                                    h[j][en] /= t;
                                }
                            }
                        }
                    }
                }
            }
        }
    
    for (i=0; i<dim; i++)
        {
        if ((i < low) || (i > high)) /* changed to rid goto statement */
            {
            for (j=i; j<dim; j++)
                z[i][j] = h[i][j];
            }
        }

    /* multiply by transformation matrix to give vectors of original
       full matrix */
    for (j=dim-1; j>=low; j--)
        {
        m = MIN(j, high);
        for (i=low; i<=high; i++)
            {
            zz = 0.0;
            for (k = low; k <= m; k++)
                zz += z[i][k] * h[k][j];
            z[i][j] = zz;
            }
        }

    return (0);
    
#if 0
int hqr2 (int *nm, int *n, int *low, int *igh, MrBFlt *h__, MrBFlt *wr, MrBFlt *wi, MrBFlt *z__, int *ierr)
    
{

    /* system generated locals */
    int h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3;
    MrBFlt d__1, d__2, d__3, d__4;

    /* builtin functions */
    MrBFlt sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static MrBFlt norm;
    static int i__, j, k, l, m;
    static MrBFlt p, q, r__, s, t, w, x, y;
    static int na, ii, en, jj;
    static MrBFlt ra, sa;
    static int ll, mm, nn;
    static MrBFlt vi, vr, zz;
    static logical notlas;
    static int mp2, itn, its, enm2;
    static MrBFlt tst1, tst2;

    /* parameter adjustments */
    z_dim1 = *nm;
    z_offset = z_dim1 + 1;
    z__ -= z_offset;
    --wi;
    --wr;
    h_dim1 = *nm;
    h_offset = h_dim1 + 1;
    h__ -= h_offset;

    /* function Body */
    *ierr = 0;
    norm = 0.;
    k = 1;
    
    /* .......... store roots isolated by balanc and compute matrix norm .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) 
        {
        i__2 = *n;
        for (j = k; j <= i__2; ++j) 
            {
            /* L40: */
            norm += (d__1 = h__[i__ + j * h_dim1], abs(d__1));
            }
        k = i__;
        if (i__ >= *low && i__ <= *igh) 
            goto L50;
        wr[i__] = h__[i__ + i__ * h_dim1];
        wi[i__] = 0.;
        L50:
            ;
        }

    en = *igh;
    t = 0.;
    itn = *n * 30;
    
    /* ..........search for next eigenvalues.......... */
    L60:
    if (en < *low) 
        goto L340;
    its = 0;
    na = en - 1;
    enm2 = na - 1;
    
    /* ..........look for single small sub-diagonal element for l=en step -1 until low do -- .......... */
    L70:
    i__1 = en;
    for (ll = *low; ll <= i__1; ++ll) 
        {
        l = en + *low - ll;
        if (l == *low) 
            goto L100;
        s = (d__1 = h__[l - 1 + (l - 1) * h_dim1], abs(d__1)) + (d__2 = h__[l + l * h_dim1], abs(d__2));
        if (s == 0.0) 
            s = norm;
        tst1 = s;
        tst2 = tst1 + (d__1 = h__[l + (l - 1) * h_dim1], abs(d__1));
        if (tst2 == tst1) 
            goto L100;
        /* L80: */
        }
        
    /* .......... form shift .......... */
    L100:
    x = h__[en + en * h_dim1];
    if (l == en) 
        goto L270;
    y = h__[na + na * h_dim1];
    w = h__[en + na * h_dim1] * h__[na + en * h_dim1];
    if (l == na) 
        goto L280;
    if (itn == 0) 
        goto L1000;
    if (its != 10 && its != 20) 
        goto L130;

    /* .......... form exceptional shift .......... */
    t += x;

    i__1 = en;
    for (i__ = *low; i__ <= i__1; ++i__) 
        {
        /* L120: */
        h__[i__ + i__ * h_dim1] -= x;
        }

    s = (d__1 = h__[en + na * h_dim1], abs(d__1)) + (d__2 = h__[na + enm2 * h_dim1], abs(d__2));
    x = s * 0.75;
    y = x;
    w = s * -0.4375 * s;
    L130:
    ++its;
    --itn;
    
    /* .......... look for two consecutive small sub-diagonal elements for m=en-2 step -1 until l do -- .......... */
    i__1 = enm2;
    for (mm = l; mm <= i__1; ++mm) 
        {
        m = enm2 + l - mm;
        zz = h__[m + m * h_dim1];
        r__ = x - zz;
        s = y - zz;
        p = (r__ * s - w) / h__[m + 1 + m * h_dim1] + h__[m + (m + 1) * h_dim1];
        q = h__[m + 1 + (m + 1) * h_dim1] - zz - r__ - s;
        r__ = h__[m + 2 + (m + 1) * h_dim1];
        s = abs(p) + abs(q) + abs(r__);
        p /= s;
        q /= s;
        r__ /= s;
        if (m == l) 
            goto L150;
        tst1 = abs(p) * ((d__1 = h__[m - 1 + (m - 1) * h_dim1], abs(d__1)) + 
        abs(zz) + (d__2 = h__[m + 1 + (m + 1) * h_dim1], abs(d__2)));
        tst2 = tst1 + (d__1 = h__[m + (m - 1) * h_dim1], abs(d__1)) * (abs(q) + abs(r__));
        if (tst2 == tst1) 
            goto L150;
        /* L140: */
        }
    L150:
    mp2 = m + 2;

    i__1 = en;
    for (i__ = mp2; i__ <= i__1; ++i__) 
        {
        h__[i__ + (i__ - 2) * h_dim1] = 0.0;
        if (i__ == mp2)
            goto L160;
        h__[i__ + (i__ - 3) * h_dim1] = 0.;
        L160:
            ;
        }
        
    /*     .......... MrBFlt qr step involving rows l to en and columns m to en .......... */
    i__1 = na;
    for (k = m; k <= i__1; ++k) 
        {
        notlas = k != na;
        if (k == m) 
            goto L170;
        p = h__[k + (k - 1) * h_dim1];
        q = h__[k + 1 + (k - 1) * h_dim1];
        r__ = 0.;
        if (notlas) 
            r__ = h__[k + 2 + (k - 1) * h_dim1];
        x = abs(p) + abs(q) + abs(r__);
        if (x == 0.) 
            goto L260;
        p /= x;
        q /= x;
        r__ /= x;
        L170:
        d__1 = sqrt(p * p + q * q + r__ * r__);
        s = d_sign(&d__1, &p);
        if (k == m) 
            goto L180;
        h__[k + (k - 1) * h_dim1] = -s * x;
        goto L190;
        L180:
        if (l != m) 
            {
            h__[k + (k - 1) * h_dim1] = -h__[k + (k - 1) * h_dim1];
            }
        L190:
        p += s;
        x = p / s;
        y = q / s;
        zz = r__ / s;
        q /= p;
        r__ /= p;
        if (notlas) 
            goto L225;
        
        /* .......... row modification .......... */
        i__2 = *n;
        for (j = k; j <= i__2; ++j) 
            {
            p = h__[k + j * h_dim1] + q * h__[k + 1 + j * h_dim1];
            h__[k + j * h_dim1] -= p * x;
            h__[k + 1 + j * h_dim1] -= p * y;
            /* L200: */
            }

        /* computing MIN */
        i__2 = en, i__3 = k + 3;
        j = min(i__2,i__3);
        
        /* .......... column modification .......... */
        i__2 = j;
        for (i__ = 1; i__ <= i__2; ++i__) 
            {
            p = x * h__[i__ + k * h_dim1] + y * h__[i__ + (k + 1) * h_dim1];
            h__[i__ + k * h_dim1] -= p;
            h__[i__ + (k + 1) * h_dim1] -= p * q;
            /* L210: */
            }
            
        /* .......... accumulate transformations .......... */
        i__2 = *igh;
        for (i__ = *low; i__ <= i__2; ++i__) 
            {
            p = x * z__[i__ + k * z_dim1] + y * z__[i__ + (k + 1) * z_dim1];
            z__[i__ + k * z_dim1] -= p;
            z__[i__ + (k + 1) * z_dim1] -= p * q;
            /* L220: */
            }
        goto L255;
        L225:
        
        /* .......... row modification .......... */
        i__2 = *n;
        for (j = k; j <= i__2; ++j) 
            {
            p = h__[k + j * h_dim1] + q * h__[k + 1 + j * h_dim1] + r__ * h__[k + 2 + j * h_dim1];
            h__[k + j * h_dim1] -= p * x;
            h__[k + 1 + j * h_dim1] -= p * y;
            h__[k + 2 + j * h_dim1] -= p * zz;
            /* L230: */
            }

        /* computing MIN */
        i__2 = en, i__3 = k + 3;
        j = min(i__2,i__3);
        
        /* .......... column modification .......... */
        i__2 = j;
        for (i__ = 1; i__ <= i__2; ++i__) 
            {
            p = x * h__[i__ + k * h_dim1] + y * h__[i__ + (k + 1) * h_dim1] + 
            zz * h__[i__ + (k + 2) * h_dim1];
            h__[i__ + k * h_dim1] -= p;
            h__[i__ + (k + 1) * h_dim1] -= p * q;
            h__[i__ + (k + 2) * h_dim1] -= p * r__;
            /* L240: */
            }
        
        /* .......... accumulate transformations .......... */
        i__2 = *igh;
        for (i__ = *low; i__ <= i__2; ++i__) 
            {
            p = x * z__[i__ + k * z_dim1] + y * z__[i__ + (k + 1) * z_dim1] + zz * z__[i__ + (k + 2) * z_dim1];
            z__[i__ + k * z_dim1] -= p;
            z__[i__ + (k + 1) * z_dim1] -= p * q;
            z__[i__ + (k + 2) * z_dim1] -= p * r__;
            /* L250: */
            }
        L255:
        L260:
            ;
        }
    goto L70;
    
    /* .......... one root found .......... */
    L270:
    h__[en + en * h_dim1] = x + t;
    wr[en] = h__[en + en * h_dim1];
    wi[en] = 0.;
    en = na;
    goto L60;
    
    /* .......... two roots found .......... */
    L280:
    p = (y - x) / 2.;
    q = p * p + w;
    zz = sqrt((abs(q)));
    h__[en + en * h_dim1] = x + t;
    x = h__[en + en * h_dim1];
    h__[na + na * h_dim1] = y + t;
    if (q < 0.) 
        goto L320;
    
    /* .......... real pair .......... */
    zz = p + d_sign(&zz, &p);
    wr[na] = x + zz;
    wr[en] = wr[na];
    if (zz != 0.) 
        {
        wr[en] = x - w / zz;
        }
    wi[na] = 0.0;
    wi[en] = 0.0;
    x = h__[en + na * h_dim1];
    s = abs(x) + abs(zz);
    p = x / s;
    q = zz / s;
    r__ = sqrt(p * p + q * q);
    p /= r__;
    q /= r__;
    
    /* .......... row modification .......... */
    i__1 = *n;
    for (j = na; j <= i__1; ++j) 
        {
        zz = h__[na + j * h_dim1];
        h__[na + j * h_dim1] = q * zz + p * h__[en + j * h_dim1];
        h__[en + j * h_dim1] = q * h__[en + j * h_dim1] - p * zz;
        /* L290: */
        }
    
    /* .......... column modification .......... */
    i__1 = en;
    for (i__ = 1; i__ <= i__1; ++i__) 
        {
        zz = h__[i__ + na * h_dim1];
        h__[i__ + na * h_dim1] = q * zz + p * h__[i__ + en * h_dim1];
        h__[i__ + en * h_dim1] = q * h__[i__ + en * h_dim1] - p * zz;
        /* L300: */
        }
        
    /* .......... accumulate transformations .......... */
    i__1 = *igh;
    for (i__ = *low; i__ <= i__1; ++i__) 
        {
        zz = z__[i__ + na * z_dim1];
        z__[i__ + na * z_dim1] = q * zz + p * z__[i__ + en * z_dim1];
        z__[i__ + en * z_dim1] = q * z__[i__ + en * z_dim1] - p * zz;
        /* L310: */
        }
    goto L330;
    
    /* .......... complex pair .......... */
    L320:
    wr[na] = x + p;
    wr[en] = x + p;
    wi[na] = zz;
    wi[en] = -zz;
    L330:
    en = enm2;
    goto L60;
    
    /* .......... all roots found.  backsubstitute to find vectors of upper triangular form .......... */
    L340:
    if (norm == 0.0) 
        goto L1001;

    /* .......... for en=n step -1 until 1 do -- .......... */
    i__1 = *n;
    for (nn = 1; nn <= i__1; ++nn) 
        {
        en = *n + 1 - nn;
        p = wr[en];
        q = wi[en];
        na = en - 1;
        if (q < 0.) 
            goto L710;
        else if (q == 0) 
            goto L600;
        else 
            goto L800;
            
        /* .......... real vector .......... */
        L600:
        m = en;
        h__[en + en * h_dim1] = 1.0;
        if (na == 0) 
            goto L800;
        
        /*     .......... for i=en-1 step -1 until 1 do -- .......... */
        i__2 = na;
        for (ii = 1; ii <= i__2; ++ii) 
            {
            i__ = en - ii;
            w = h__[i__ + i__ * h_dim1] - p;
            r__ = 0.0;

            i__3 = en;
            for (j = m; j <= i__3; ++j) 
                {
                /* L610: */
                r__ += h__[i__ + j * h_dim1] * h__[j + en * h_dim1];
                }

            if (wi[i__] >= 0.0) 
                goto L630;
            zz = w;
            s = r__;
            goto L700;
            L630:
            m = i__;
            if (wi[i__] != 0.0) 
                goto L640;
            t = w;
            if (t != 0.0)
                goto L635;
            tst1 = norm;
            t = tst1;
            L632:
            t *= 0.01;
            tst2 = norm + t;
            if (tst2 > tst1) 
                goto L632;
            L635:
            h__[i__ + en * h_dim1] = -r__ / t;
            goto L680;
            
            /* .......... solve real equations .......... */
            L640:
            x = h__[i__ + (i__ + 1) * h_dim1];
            y = h__[i__ + 1 + i__ * h_dim1];
            q = (wr[i__] - p) * (wr[i__] - p) + wi[i__] * wi[i__];
            t = (x * s - zz * r__) / q;
            h__[i__ + en * h_dim1] = t;
            if (abs(x) <= abs(zz)) 
                goto L650;
            h__[i__ + 1 + en * h_dim1] = (-r__ - w * t) / x;
            goto L680;
            L650:
            h__[i__ + 1 + en * h_dim1] = (-s - y * t) / zz;

            /*     .......... overflow control .......... */
            L680:
            t = (d__1 = h__[i__ + en * h_dim1], abs(d__1));
            if (t == 0.0) 
                goto L700;
            tst1 = t;
            tst2 = tst1 + 1.0 / tst1;
            if (tst2 > tst1) 
                goto L700;
            i__3 = en;
            for (j = i__; j <= i__3; ++j) 
                {
                h__[j + en * h_dim1] /= t;
                /* L690: */
                }

            L700:
                ;
            }
            
        /* .......... end real vector .......... */
        goto L800;
        
        /* .......... complex vector .......... */
        L710:
        m = na;
        
        /* .......... last vector component chosen imaginary so that eigenvector matrix is triangular .......... */
        if ((d__1 = h__[en + na * h_dim1], abs(d__1)) <= (d__2 = h__[na + en *
        h_dim1], abs(d__2))) 
            goto L720;
        h__[na + na * h_dim1] = q / h__[en + na * h_dim1];
        h__[na + en * h_dim1] = -(h__[en + en * h_dim1] - p) / h__[en + na * h_dim1];
        goto L730;
        L720:
        d__1 = -h__[na + en * h_dim1];
        d__2 = h__[na + na * h_dim1] - p;
        cdiv_(&c_b49, &d__1, &d__2, &q, &h__[na + na * h_dim1], &h__[na + en *
        h_dim1]);
        L730:
        h__[en + na * h_dim1] = 0.0;
        h__[en + en * h_dim1] = 1.0;
        enm2 = na - 1;
        if (enm2 == 0) 
            goto L800;

        /*     .......... for i=en-2 step -1 until 1 do -- .......... */
        i__2 = enm2;
        for (ii = 1; ii <= i__2; ++ii) 
            {
            i__ = na - ii;
            w = h__[i__ + i__ * h_dim1] - p;
            ra = 0.0;
            sa = 0.0;

            i__3 = en;
            for (j = m; j <= i__3; ++j) 
                {
                ra += h__[i__ + j * h_dim1] * h__[j + na * h_dim1];
                sa += h__[i__ + j * h_dim1] * h__[j + en * h_dim1];
                /* L760: */
                }

            if (wi[i__] >= 0.0) 
                goto L770;
            zz = w;
            r__ = ra;
            s = sa;
            goto L795;
            L770:
            m = i__;
            if (wi[i__] != 0.0) 
                goto L780;
            d__1 = -ra;
            d__2 = -sa;
            cdiv_(&d__1, &d__2, &w, &q, &h__[i__ + na * h_dim1], &h__[i__ + en * h_dim1]);
            goto L790;
            
            /*     .......... solve complex equations .......... */
            L780:
            x = h__[i__ + (i__ + 1) * h_dim1];
            y = h__[i__ + 1 + i__ * h_dim1];
            vr = (wr[i__] - p) * (wr[i__] - p) + wi[i__] * wi[i__] - q * q;
            vi = (wr[i__] - p) * 2.0 * q;
            if (vr != 0.0 || vi != 0.0) 
                goto L784;
            tst1 = norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz));
            vr = tst1;
            L783:
            vr *= 0.01;
            tst2 = tst1 + vr;
            if (tst2 > tst1) 
                goto L783;
            L784:
            d__1 = x * r__ - zz * ra + q * sa;
            d__2 = x * s - zz * sa - q * ra;
            cdiv_(&d__1, &d__2, &vr, &vi, &h__[i__ + na * h_dim1], &h__[i__ + en * h_dim1]);
            if (abs(x) <= abs(zz) + abs(q)) 
                goto L785;
            h__[i__ + 1 + na * h_dim1] = (-ra - w * h__[i__ + na * h_dim1] + q * h__[i__ + en * h_dim1]) / x;
            h__[i__ + 1 + en * h_dim1] = (-sa - w * h__[i__ + en * h_dim1] - q * h__[i__ + na * h_dim1]) / x;
            goto L790;
            L785:
            d__1 = -r__ - y * h__[i__ + na * h_dim1];
            d__2 = -s - y * h__[i__ + en * h_dim1];
            cdiv_(&d__1, &d__2, &zz, &q, &h__[i__ + 1 + na * h_dim1], &h__[i__ + 1 + en * h_dim1]);

            /*     .......... overflow control .......... */
            L790:
            /* Computing MAX */
            d__3 = (d__1 = h__[i__ + na * h_dim1], abs(d__1)), d__4 = (d__2 = h__[i__ + en * h_dim1], abs(d__2));
            t = max(d__3,d__4);
            if (t == 0.0) 
                goto L795;
            tst1 = t;
            tst2 = tst1 + 1.0 / tst1;
            if (tst2 > tst1) 
                goto L795;
            i__3 = en;
            for (j = i__; j <= i__3; ++j) 
                {
                h__[j + na * h_dim1] /= t;
                h__[j + en * h_dim1] /= t;
                /* L792: */
                }
            L795:
                ;
            }
        /*     .......... end complex vector .......... */
        L800:
            ;
        }
    /*     .......... end back substitution vectors of isolated roots .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) 
        {
        if (i__ >= *low && i__ <= *igh) 
            goto L840;
        i__2 = *n;
        for (j = i__; j <= i__2; ++j) 
            {
            /* L820: */
            z__[i__ + j * z_dim1] = h__[i__ + j * h_dim1];
            }
        L840:
        ;
        }
        
    /* .......... multiply by transformation matrix to give vectors of original full matrix. */
    /*            for j=n step -1 until low do -- .......... */
    i__1 = *n;
    for (jj = *low; jj <= i__1; ++jj) 
        {
        j = *n + *low - jj;
        m = min(j,*igh);

        i__2 = *igh;
        for (i__ = *low; i__ <= i__2; ++i__) 
            {
            zz = 0.0;
            i__3 = m;
            for (k = *low; k <= i__3; ++k) 
                {
                /* L860: */
                zz += z__[i__ + k * z_dim1] * h__[k + j * h_dim1];
                }

            z__[i__ + j * z_dim1] = zz;
            /* L880: */
            }
        }

    goto L1001;
    /* .......... set error -- all eigenvalues have not converged after 30*n iterations .......... */
    L1000:
    *ierr = en;
    L1001:
    return 0;
    
}
/* end f2c version of code */
#endif

}


MrBFlt IncompleteBetaFunction (MrBFlt alpha, MrBFlt beta, MrBFlt x)
{
    MrBFlt      bt, gm1, gm2, gm3, temp;
    
    if (x < 0.0 || x > 1.0) 
        {
        MrBayesPrint ("%s   Error: Problem in IncompleteBetaFunction.\n", spacer);
        exit (0);
        }
    if (fabs(x) < ETA || fabs(x-1.0)<ETA) /* x == 0.0 || x == 1.0 */
        {
        bt = 0.0;
        }
    else
        {
        gm1 = LnGamma (alpha + beta);
        gm2 = LnGamma (alpha);
        gm3 = LnGamma (beta);
        temp = gm1 - gm2 - gm3 + (alpha) * log(x) + (beta) * log(1.0 - x);
        bt = exp(temp);
        }
    if (x < (alpha + 1.0)/(alpha + beta + 2.0))
        return (bt * BetaCf(alpha, beta, x) / alpha);
    else
        return (1.0 - bt * BetaCf(beta, alpha, 1.0-x) / beta);
}


/*---------------------------------------------------------------------------------
|
|   IncompleteGamma
|
|   Returns the incomplete gamma ratio I(x,alpha) where x is the upper
|   limit of the integration and alpha is the shape parameter.  Returns (-1)
|   if in error.   
|
|   Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied
|      Statistics, 19:285-287 (AS32)
|
---------------------------------------------------------------------------------*/
MrBFlt IncompleteGamma (MrBFlt x, MrBFlt alpha, MrBFlt LnGamma_alpha)
{
    int             i;
    MrBFlt      p = alpha, g = LnGamma_alpha,
                    accurate = 1e-8, overflow = 1e30,
                    factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0, 
                    dif = 0.0, term = 0.0, pn[6];

    if (fabs(x) < ETA) 
        return (0.0);
    if (x < 0 || p <= 0) 
        return (-1.0);

    factor = exp(p*log(x)-x-g);   
    if (x>1 && x>=p) 
        goto l30;
    gin = 1.0;  
    term = 1.0;  
    rn = p;
    l20:
        rn++;
        term *= x/rn;   
        gin += term;
        if (term > accurate) 
            goto l20;
        gin *= factor/p;
        goto l50;
    l30:
        a = 1.0-p;   
        b = a+x+1.0;  
        term = 0.0;
        pn[0] = 1.0;  
        pn[1] = x;  
        pn[2] = x+1;  
        pn[3] = x*b;
        gin = pn[2]/pn[3];
    l32:
        a++;  
        b += 2.0;  
        term++;   
        an = a*term;
        for (i=0; i<2; i++) 
            pn[i+4] = b*pn[i+2]-an*pn[i];
        if (fabs(pn[5]) < ETA) 
            goto l35;
        rn = pn[4]/pn[5];   
        dif = fabs(gin-rn);
        if (dif>accurate) 
            goto l34;
        if (dif<=accurate*rn) 
            goto l42;
    l34:
        gin = rn;
    l35:
        for (i=0; i<4; i++) 
            pn[i] = pn[i+2];
        if (fabs(pn[4]) < overflow) 
            goto l32;
        for (i=0; i<4; i++) 
            pn[i] /= overflow;
        goto l32;
    l42:
        gin = 1.0-factor*gin;
    l50:
        return (gin);
}


/*---------------------------------------------------------------------------------
|
|   InvertMatrix
|
|   Calculates aInv = a^{-1} using LU-decomposition. The input matrix a is
|   destroyed in the process. The program returns an error if the matrix is
|   singular. col and indx are work vectors.
|      
---------------------------------------------------------------------------------*/
int InvertMatrix (int dim, MrBFlt **a, MrBFlt *col, int *indx, MrBFlt **aInv)
{
    int         rc, i, j;
    
    rc = LUDecompose (dim, a, col, indx, (MrBFlt *)NULL);
    if (rc == FALSE)
        {
        for (j = 0; j < dim; j++)
            {
            for (i = 0; i < dim; i++)
                col[i] = 0.0;
            col[j] = 1.0;
            LUBackSubstitution (dim, a, indx, col);
            for (i = 0; i < dim; i++)
                aInv[i][j] = col[i];
            }
        }
        
    return (rc);
}


/*---------------------------------------------------------------------------------
|
|   LBinormal
|
|   L(h1,h2,r) = prob(x>h1, y>h2), where x and y are standard binormal, 
|   with r=corr(x,y),  error < 2e-7.
|
|   Drezner Z., and G.O. Wesolowsky (1990) On the computation of the
|      bivariate normal integral.  J. Statist. Comput. Simul. 35:101-107.
|
---------------------------------------------------------------------------------*/
MrBFlt LBinormal (MrBFlt h1, MrBFlt h2, MrBFlt r)
{
    int i;
    MrBFlt      x[]={0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992};
    MrBFlt      w[]={0.018854042, 0.038088059, 0.0452707394,0.038088059,0.018854042};
    MrBFlt      Lh=0.0, r1, r2, r3, rr, aa, ab, h3, h5, h6, h7, h12, temp1, temp2, exp1, exp2;

    h12 = (h1 * h1 + h2 * h2) / 2.0;
    if (fabs(r) >= 0.7) 
        {
        r2 = 1.0 - r * r;   
        r3 = sqrt(r2);
        if (r < 0) 
            h2 *= -1;
        h3 = h1 * h2;   
        h7 = exp(-h3 / 2.0);
        if (fabs(r-1.0)>ETA)  /* fabs(r) != 1.0 */
            {
            h6 = fabs(h1-h2);  
            h5 = h6 * h6 / 2.0; 
            h6 /= r3; 
            aa = 0.5 - h3 / 8;  
            ab = 3.0 - 2.0 * aa * h5;
            temp1 = -h5 / r2;
            if (temp1 < -100.0)
                exp1 = 0.0;
            else
                exp1 = exp(temp1);
            Lh = 0.13298076 * h6 * ab * (1.0 - CdfNormal(h6)) - exp1 * (ab + aa * r2) * 0.053051647;
            for (i=0; i<5; i++) 
                {
                r1 = r3 * x[i];
                rr = r1 * r1;   
                r2 = sqrt(1.0 - rr);
                temp1 = -h5 / rr;
                if (temp1 < -100.0)
                    exp1 = 0.0;
                else
                    exp1 = exp(temp1);
                temp2 = -h3 / (1.0 + r2);
                if (temp2 < -100.0)
                    exp2 = 0.0;
                else
                    exp2 = exp(temp2);
                Lh -= w[i] * exp1 * (exp2 / r2 / h7 - 1.0 - aa * rr);
                }
            }
        if (r > 0) 
            Lh = Lh * r3 * h7 + (1.0 - CdfNormal(MAX(h1, h2)));
        else if (r<0) 
            Lh = (h1 < h2 ? CdfNormal(h2) - CdfNormal(h1) : 0) - Lh * r3 * h7;
        }
    else 
        {
        h3 = h1 * h2;
        if (fabs(r)>ETA) 
            {
            for (i=0; i<5; i++) 
                {
                r1 = r * x[i]; 
                r2 = 1.0 - r1 * r1;
                temp1 = (r1 * h3 - h12) / r2;
                if (temp1 < -100.0)
                    exp1 = 0.0;
                else
                    exp1 = exp(temp1);
                Lh += w[i] * exp1 / sqrt(r2);
                }
            }
        Lh = (1.0 - CdfNormal(h1)) * (1.0 - CdfNormal(h2)) + r * Lh;
        }
    return (Lh);
}


/*---------------------------------------------------------------------------------
|
|   LnFactorial: Calculates the log of the factorial for an integer
|
---------------------------------------------------------------------------------*/
MrBFlt  LnFactorial (int value)
{
    int     i;
    MrBFlt  result;

    result = 0.0;

    for (i = 2; i<=value; i++)
        result += log(i);

    return result;
}


/*---------------------------------------------------------------------------------
|
|   LnGamma
|
|   Calculates the log of the gamma function. The Gamma function is equal
|   to:
|
|      Gamma(alp) = {integral from 0 to infinity} t^{alp-1} e^-t dt
|
|   The result is accurate to 10 decimal places. Stirling's formula is used
|   for the central polynomial part of the procedure.
|
|   Pike, M. C. and I. D. Hill.  1966.  Algorithm 291: Logarithm of the gamma
|      function.  Communications of the Association for Computing
|      Machinery, 9:684.
|      
---------------------------------------------------------------------------------*/
MrBFlt LnGamma (MrBFlt alp)
{
    MrBFlt      x = alp, f = 0.0, z;
    
    if (x < 7) 
        {
        f = 1.0;
        z = x-1.0;
        while (++z < 7.0)  
            f *= z;
        x = z;   
        f = -log(f);
        }
    z = 1.0 / (x*x);
    return  (f + (x-0.5)*log(x) - x + 0.918938533204673 + 
            (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +0.083333333333333)/x);
}


/* Calculate probability of a realization for exponential random variable */
MrBFlt LnPriorProbExponential (MrBFlt val, MrBFlt *params)
{
    return log(params[0]) - params[0] * val;
}


/* Calculate probability of a realization for exponential random variable; parameter mean and not rate */
MrBFlt LnPriorProbExponential_Param_Mean (MrBFlt val, MrBFlt *params)
{
    return - log(params[0]) - val / params[0];
}


/* Calculate probability of a realization for a fixed variable */
MrBFlt LnPriorProbFix (MrBFlt val, MrBFlt *params)
{
    if (fabs((val - params[0])/val) < 1E-5)
        return 0.0;
    else
        return NEG_INFINITY;
}


/* Calculate probability of a realization for gamma random variable */
MrBFlt LnPriorProbGamma (MrBFlt val, MrBFlt *params)
{
    return (params[0] - 1) * log(val) + params[0] * log(params[1]) - params[1] * val - LnGamma(params[0]);
}


/* Calculate probability of a realization for gamma random variable; parameters mean and sd */
MrBFlt LnPriorProbGamma_Param_Mean_Sd (MrBFlt val, MrBFlt *params)
{
    MrBFlt  alpha, beta;

    beta  = params[0] / (params[1]*params[1]);
    alpha = params[0] * beta;

    return (alpha - 1) * log(val) + alpha * log(beta) - beta * val - LnGamma(alpha);
}


/* Calculate probability of a realization for lognormal random variable */
MrBFlt LnPriorProbLognormal (MrBFlt val, MrBFlt *params)
{
    MrBFlt z;

    z = (log(val) - params[0]) / params[1];

    return - log(params[1] * val * sqrt(2.0 * PI)) - z * z / 2.0;
}


/* Calculate probability of a realization for lognormal random variable; parameters mean and sd on linear scale */
MrBFlt LnPriorProbLognormal_Param_Mean_Sd (MrBFlt val, MrBFlt *params)
{
    MrBFlt z, mean_log, sd_log;

    sd_log      = sqrt (log((params[1]*params[1])/(params[0]*params[0]) + 1));
    mean_log    = log(params[0]) - sd_log * sd_log / 2.0;

    z= (log(val) - mean_log) / sd_log;

    return - log(sd_log * val * sqrt(2.0 * PI)) - z * z / 2.0;
}


/* Calculate probability of a realization for normal random variable */
MrBFlt LnPriorProbNormal (MrBFlt val, MrBFlt *params)
{
    MrBFlt z;

    z = (val - params[0]) / params[1];

    return - log(params[1] * sqrt(2.0 * PI)) - z * z / 2.0;
}


/* Calculate probability of a realization for an offset exponential random variable */
MrBFlt LnPriorProbOffsetExponential (MrBFlt val, MrBFlt *params)
{
    return log(params[1]) - params[1] * (val - params[0]);
}


/* Calculate probability of a realization for an offset exponential random variable; parameters offset and mean */
MrBFlt LnPriorProbOffsetExponential_Param_Offset_Mean (MrBFlt val, MrBFlt *params)
{
    MrBFlt  x, rate;

    x    = val - params[0];
    rate = 1.0 / (params[1] - params[0]);

    return log(rate) - rate * x;
}


/* Calculate probability of a realization for an offset gamma random variable */
MrBFlt LnPriorProbOffsetGamma (MrBFlt val, MrBFlt *params)
{
    MrBFlt x, alpha, beta;
    
    x     = val - params[0];
    alpha = params[1];
    beta  = params[2];

    return (alpha - 1) * log(x) + alpha * log(beta) - beta * x - LnGamma(alpha);
}


/* Calculate probability of a realization for an offset gamma random variable; parameters offset, mean and sd */
MrBFlt LnPriorProbOffsetGamma_Param_Offset_Mean_Sd (MrBFlt val, MrBFlt *params)
{
    MrBFlt  x, mean, sd, alpha, beta;

    x     = val - params[0];
    mean  = params[1] - params[0];
    sd    = params[2];

    beta  = mean / (sd*sd);
    alpha = mean * beta;

    return (alpha - 1) * log(x) + alpha * log(beta) - beta * x - LnGamma(alpha);
}


/* Calculate probability of a realization for an offset lognormal random variable */
MrBFlt LnPriorProbOffsetLognormal (MrBFlt val, MrBFlt *params)
{
    MrBFlt x, mean_log, sd_log, z;

    x        = val - params[0];
    mean_log = params[1] - params[0];
    sd_log   = params[2];

    z = (log(x) - mean_log) / sd_log;

    return - log(sd_log * x * sqrt(2.0 * PI)) - z * z / 2.0;
}


/* Calculate probability of a realization for an offset lognormal random variable; parameters offset, mean and sd */
MrBFlt LnPriorProbOffsetLognormal_Param_Offset_Mean_Sd (MrBFlt val, MrBFlt *params)
{
    MrBFlt x, mean, sd, mean_log, sd_log, z;

    x        = val - params[0];
    mean     = params[1] - params[0];
    sd       = params[2];
    sd_log   = sqrt (log((sd*sd)/(mean*mean) + 1));
    mean_log = log(mean) - sd_log * sd_log / 2.0;

    z = (log(x) - mean_log) / sd_log;

    return - log(sd_log * x * sqrt(2.0 * PI)) - z * z / 2.0;
}


/* Calculate probability of a realization for truncated (only positive values) normal random variable */
MrBFlt LnPriorProbTruncatedNormal (MrBFlt val, MrBFlt *params)
{
    MrBFlt z, z_0, normConst;

    z = (val - params[0]) / params[1];
    z_0 = (0.0 - params[0]) / params[1];
    normConst = CdfNormal(z_0);
    
    return - log(params[1] * sqrt(2.0 * PI)) - z * z / 2.0 - log(1.0 - normConst);
}


/* Calculate probability of a realization for arbitrarily truncated normal random variable; parameters truncation point, mean and sd */
MrBFlt LnPriorProbTruncatedNormal_Param_Trunc_Mean_Sd (MrBFlt val, MrBFlt *params)
{
    MrBFlt z, z_trunc, normConst;

    z = (val - params[1]) / params[2];
    z_trunc = (params[0] - params[1]) / params[2];
    normConst = CdfNormal(z_trunc);

    return - log(params[2] * sqrt(2.0 * PI)) - z * z / 2.0 - log(1.0 - normConst);
}


/* Calculate probability of a realization for uniform random variable */
MrBFlt LnPriorProbUniform (MrBFlt val, MrBFlt *params)
{
    return - log(params[1] - params[0]);
    MrBayesPrint ("%lf", val); /* just because I am tired of seeing the unused parameter error msg */
}


/* Calculate probability ratio of realizations for exponential random variable */
MrBFlt LnProbRatioExponential (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    return params[0] * (oldX - newX);
}


/* Calculate probability ratio of realizations for exponential random variable; parameter mean and not rate */
MrBFlt LnProbRatioExponential_Param_Mean (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    return (oldX - newX) / params[0];
}


/* Calculate probability of a realization for a fixed variable */
MrBFlt LnProbRatioFix (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    if (fabs((newX - params[0])/newX) < 1E-5 && fabs((oldX - params[0])/oldX) < 1E-5)
        return 0.0;
    else
        return NEG_INFINITY;
}


/* Calculate probability ratio of realizations for gamma random variable */
MrBFlt LnProbRatioGamma (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt  alpha, beta;

    alpha   = params[0];
    beta    = params[1];

    return (alpha - 1.0) * (log(newX) - log(oldX)) - beta * (newX - oldX);
}


/* Calculate probability ratio of realizations for gamma random variable; parameters mean and sd */
MrBFlt LnProbRatioGamma_Param_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt  alpha, beta;

    beta  = params[0] / (params[1]*params[1]);
    alpha = params[0] * beta;

    return (alpha - 1.0) * (log(newX) - log(oldX)) - beta * (newX - oldX);
}


/* Calculate probability ratio of realizations for log normal random variable */
MrBFlt LnProbRatioLognormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt  newZ, oldZ;

    newZ = (log(newX) - params[0]) / params[1];
    oldZ = (log(oldX) - params[0]) / params[1];

    return (oldZ * oldZ - newZ * newZ) / 2.0 + log(oldX) - log(newX);
}


/* Calculate probability ratio of realizations for log normal random variable; parameters mean and sd */
MrBFlt LnProbRatioLognormal_Param_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{    
    MrBFlt newZ, oldZ, mean_log, sd_log;

    sd_log      = sqrt (log((params[1]*params[1])/(params[0]*params[0]) + 1));
    mean_log    = log(params[0]) - sd_log * sd_log / 2.0;

    newZ = (log(newX) - mean_log) / sd_log;
    oldZ = (log(oldX) - mean_log) / sd_log;

    return (oldZ * oldZ - newZ * newZ) / 2.0 + log(oldX) - log(newX);
}


/* Calculate probability ratio of realizations for normal random variable */
MrBFlt LnProbRatioNormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt  newZ, oldZ;

    newZ = (newX - params[0]) / params[1];
    oldZ = (oldX - params[0]) / params[1];

    return (oldZ * oldZ - newZ * newZ) / 2.0;
}


/* Calculate probability ratio of realizations for offset exponential random variable */
MrBFlt LnProbRatioOffsetExponential (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    return params[1] * (oldX - newX);
}


/* Calculate probability ratio of realizations for offset exponential random variable; parameters offset and mean */
MrBFlt LnProbRatioOffsetExponential_Param_Offset_Mean (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    return (oldX - newX) / (params[1] - params[0]);
}


/* Calculate probability ratio of realizations for offset gamma random variable */
MrBFlt LnProbRatioOffsetGamma (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt  alpha, beta, newZ, oldZ;

    alpha = params[1];
    beta  = params[2];
    newZ  = newX - params[0];
    oldZ  = oldX - params[0];

    return (alpha - 1.0) * (log(newZ) - log(oldZ)) - beta * (newZ - oldZ);
}


/* Calculate probability ratio of realizations for offset gamma random variable; parameters offset, mean and sd */
MrBFlt LnProbRatioOffsetGamma_Param_Offset_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt  mean, sd, alpha, beta;

    mean  = params[1] - params[0];
    sd    = params[2];

    beta  = mean / (sd*sd);
    alpha = mean * beta;

    newX  -= params[0];
    oldX  -= params[0];

    return (alpha - 1.0) * (log(newX) - log(oldX)) - beta * (newX - oldX);
}


/* Calculate probability ratio of realizations for offset lognormal random variable */
MrBFlt LnProbRatioOffsetLognormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt newZ, oldZ, mean_log, sd_log;

    sd_log      = params[2];
    mean_log    = params[1];

    newZ = (log(newX-params[0]) - mean_log) / sd_log;
    oldZ = (log(oldX-params[0]) - mean_log) / sd_log;

    return (oldZ * oldZ - newZ * newZ) / 2.0 + log(oldX-params[0]) - log(newX-params[0]);
}


/* Calculate probability ratio of realizations for offset lognormal random variable; parameters offset, mean and sd */
MrBFlt LnProbRatioOffsetLognormal_Param_Offset_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt newZ, oldZ, mean, sd, mean_log, sd_log;

    mean        = params[1] - params[0];
    sd          = params[2];
    sd_log      = sqrt (log((sd*sd)/(mean*mean) + 1));
    mean_log    = log(mean) - sd_log * sd_log / 2.0;

    newX -= params[0];
    oldX -= params[0];
    newZ = (log(newX) - mean_log) / sd_log;
    oldZ = (log(oldX) - mean_log) / sd_log;

    return (oldZ * oldZ - newZ * newZ) / 2.0 - log(newX / oldX);
}


/* Calculate probability ratio of realizations for truncated normal random variable */
MrBFlt LnProbRatioTruncatedNormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt  newZ, oldZ;

    newZ = (newX - params[0]) / params[1];
    oldZ = (oldX - params[0]) / params[1];

    return (oldZ * oldZ - newZ * newZ) / 2.0;
}


/* Calculate probability ratio of realizations for arbitrarily truncated normal random variable; parameters truncation point, mean and sd */
MrBFlt LnProbRatioTruncatedNormal_Param_Trunc_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    MrBFlt  newZ, oldZ;

    newZ = (newX - params[1]) / params[2];
    oldZ = (oldX - params[1]) / params[2];

    return (oldZ * oldZ - newZ * newZ) / 2.0;
}


/* Calculate probability ratio of realizations for uniform random variable */
MrBFlt LnProbRatioUniform (MrBFlt newX, MrBFlt oldX, MrBFlt *params)
{
    return 0.0;
    MrBayesPrint ("%lf %lf", newX, oldX); /* just because I am tired of seeing the unused parameter error msg */
    MrBayesPrint ("%lf", *params);
}


/* Log probability for a value drawn from a gamma distribution */
MrBFlt LnProbGamma (MrBFlt alpha, MrBFlt beta, MrBFlt x)
{
    MrBFlt lnProb;

    lnProb = (alpha-1.0)*log(x) + alpha*log(beta) - x*beta - LnGamma(alpha);

    return lnProb;
}


/* Log probability for a value drawn from a truncated gamma distribution */
MrBFlt LnProbTruncGamma (MrBFlt alpha, MrBFlt beta, MrBFlt x, MrBFlt min, MrBFlt max)
{
    MrBFlt lnProb;

    lnProb = (alpha-1.0)*log(x) + alpha*log(beta) - x*beta - LnGamma(alpha);

    lnProb -= log (IncompleteGamma (max*beta, alpha, LnGamma(alpha)) - IncompleteGamma (min*beta, alpha, LnGamma(alpha)));

    return lnProb;
}


/* Log probability for a value drawn from a lognormal distribution */
MrBFlt LnProbLogNormal (MrBFlt exp, MrBFlt sd, MrBFlt x)
{
    MrBFlt lnProb, z;
    
    z = (log(x) - exp) / sd;
    
    lnProb = - log (x * sd * sqrt (2.0 * PI)) - (z * z / 2.0);
    
    return lnProb;
}


/* Log ratio for two values drawn from a lognormal distribution */
MrBFlt LnRatioLogNormal (MrBFlt exp, MrBFlt sd, MrBFlt xNew, MrBFlt xOld)
{
    MrBFlt  newZ, oldZ;
    
    newZ = (log(xNew) - exp) / sd;
    oldZ = (log(xOld) - exp) / sd;
    
    return (oldZ * oldZ - newZ * newZ) / 2.0 + log(xOld) - log(xNew);
}


/* Log probability for a value drawn from a lognormal distribution;
   parameters are mean and variance of value (not of log value) */
MrBFlt LnProbTK02LogNormal (MrBFlt mean, MrBFlt var, MrBFlt x)
{
    MrBFlt  z, lnProb, mu, sigma;
    
    sigma = sqrt(log(1.0 + (var / (mean*mean))));
    mu    = log(mean) - sigma * sigma / 2.0;
    
    z = (log(x) - mu) / sigma;
    
    lnProb = - log (x * sigma * sqrt (2.0 * PI)) - (z * z / 2.0);
    
    return lnProb;
}


/* Log ratio for two values drawn from a lognormal distribution */
MrBFlt LnRatioTK02LogNormal (MrBFlt mean, MrBFlt var, MrBFlt xNew, MrBFlt xOld)
{
    MrBFlt  newZ, oldZ, mu, sigma;

    sigma = sqrt(log(1.0 + (var / (mean*mean))));
    mu    = log(mean) - sigma * sigma / 2.0;

    newZ = (log(xNew) - mu) / sigma;
    oldZ = (log(xOld) - mu) / sigma;

    return (oldZ * oldZ - newZ * newZ) / 2.0 + log(xOld) - log(xNew);
}


/*---------------------------------------------------------------------------------
|
|   LogBase2Plus1
|
|   This function is called from ComputeMatrixExponential.
|      
---------------------------------------------------------------------------------*/
int LogBase2Plus1 (MrBFlt x)
{
    int     j = 0;

    while (x > 1.0 - 1.0e-07) 
        {
        x /= 2.0;
        j++;
        }
        
    return (j);
}


/*---------------------------------------------------------------------------------
|
|   LogNormalRandomVariable
|
|   Draw a random variable from a lognormal distribution.
|      
---------------------------------------------------------------------------------*/
MrBFlt LogNormalRandomVariable (MrBFlt mean, MrBFlt sd, RandLong *seed)
{
    MrBFlt      x;
    
    x = PointNormal(RandomNumber(seed));

    x*= sd;
    x += mean;
    
    return exp(x);
}


/*---------------------------------------------------------------------------------
|
|   LUBackSubstitution
|
|   Back substitute into an LU-decomposed matrix.
|      
---------------------------------------------------------------------------------*/
void LUBackSubstitution (int dim, MrBFlt **a, int *indx, MrBFlt *b)
{
    int         i, ip, j, ii = -1;
    MrBFlt      sum;

    for (i=0; i<dim; i++)
        {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii >= 0)
            {
            for (j=ii; j<=i-1; j++)
                sum -= a[i][j] * b[j];
            }
        else if (fabs(sum)>ETA)
            ii = i;
        b[i] = sum;
        }
    for (i=dim-1; i>=0; i--)
        {
        sum = b[i];
        for (j=i+1; j<dim; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
        }
}


/*---------------------------------------------------------------------------------
|
|   LUDecompose
|
|   Calculate the LU-decomposition of the matrix a. The matrix a is replaced.
|      
---------------------------------------------------------------------------------*/
int LUDecompose (int dim, MrBFlt **a, MrBFlt *vv, int *indx, MrBFlt *pd)
{
    int         i, imax=0, j, k;
    MrBFlt      big, dum, sum, temp, d;

    d = 1.0;
    for (i=0; i<dim; i++)
        {
        big = 0.0;
        for (j = 0; j < dim; j++)
            {
            if ((temp = fabs(a[i][j])) > big)
                big = temp;
            }
        if (fabs(big)<ETA)
            {
            MrBayesPrint ("%s   Error: Problem in LUDecompose\n", spacer);
            return (ERROR);
            }
        vv[i] = 1.0 / big;
        }
    for (j=0; j<dim; j++)
        {
        for (i = 0; i < j; i++)
            {
            sum = a[i][j];
            for (k = 0; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            }
        big = 0.0;
        for (i=j; i<dim; i++)
            {
            sum = a[i][j];
            for (k = 0; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            dum = vv[i] * fabs(sum);
            if (dum >= big)
                {
                big = dum;
                imax = i;
                }
            }
        if (j != imax)
            {
            for (k=0; k<dim; k++)
                {
                dum = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = dum;
                }   
            d = -d;
            vv[imax] = vv[j];
            }
        indx[j] = imax;
        if (fabs(a[j][j])<ETA)
            a[j][j] = TINY;
        if (j != dim - 1)
            {
            dum = 1.0 / (a[j][j]);
            for (i=j+1; i<dim; i++)
                a[i][j] *= dum;
            }
        }
    if (pd != NULL)
        *pd = d;
        
    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
|
|   MultiplyMatrices
|
|   Multiply matrix a by matrix b and put the results in matrix result.
|
---------------------------------------------------------------------------------*/
void MultiplyMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result)
{
    register int    i, j, k;
    MrBFlt          **temp;

    temp = AllocateSquareDoubleMatrix (dim);

    for (i=0; i<dim; i++)
        {
        for (j=0; j<dim; j++) 
            {
            temp[i][j] = 0.0;
            for (k=0; k<dim; k++) 
                {
                temp[i][j] += a[i][k] * b[k][j];
                }
            }
        }
    for (i=0; i<dim; i++)
        {
        for (j=0; j<dim; j++) 
            {
            result[i][j] = temp[i][j];
            }
        }
        
    FreeSquareDoubleMatrix (temp);
}


/*---------------------------------------------------------------------------------
|
|   MultiplyMatrixByScalar
|
|   Multiply the elements of matrix a by a scalar.
|
---------------------------------------------------------------------------------*/
void MultiplyMatrixByScalar (int dim, MrBFlt **a, MrBFlt scalar, MrBFlt **result)
{
    int         row, col;

    for (row=0; row<dim; row++)
        for (col=0; col<dim; col++)
             result[row][col] = a[row][col] * scalar;
}


/*---------------------------------------------------------------------------------
|
|   MultiplyMatrixNTimes
|
---------------------------------------------------------------------------------*/
int MultiplyMatrixNTimes (int dim, MrBFlt **Mat, int power, MrBFlt **Result)
{
    register int    i, j;
    int             k, numSquares, numRemaining;
    MrBFlt          **TempIn, **TempOut;

    if (power < 0)
        {
        MrBayesPrint ("%s   Error: Power cannot be a negative number.\n", spacer);
        return (ERROR);
        }
    else if (power == 0)
        {
        for (i=0; i<dim; i++)
            for (j=0; j<dim; j++)
                Result[i][j] = 1.0;
        }
    else
        {
        TempIn  = AllocateSquareDoubleMatrix (dim);
        TempOut = AllocateSquareDoubleMatrix (dim);

        /* how many times can I multiply the matrices together */
        numSquares = 0;
        while ((1 << numSquares) < power)
            numSquares++;
        numRemaining = power - (1 << numSquares);
        
        /* now, multiply matrix by power of 2's */
        CopyDoubleMatrices (dim, Mat, TempIn);
        for (k=0; k<numSquares; k++)
            {
            MultiplyMatrices (dim, TempIn, TempIn, TempOut);
            CopyDoubleMatrices (dim, TempOut, TempIn);
            }
            
        /* TempIn is Mat^numSquares. Now, multiply it by Mat numRemaining times */
        for (k=0; k<numRemaining; k++)
            {
            MultiplyMatrices (dim, TempIn, Mat, TempOut);
            CopyDoubleMatrices (dim, TempOut, TempIn);
            }
            
        /* copy result */
        CopyDoubleMatrices (dim, TempIn, Result);
        
        FreeSquareDoubleMatrix (TempIn);
        FreeSquareDoubleMatrix (TempOut);
        }

    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
|
|   PointChi2
|
|   Returns z so that Prob(x < z) = prob where x is Chi2 distributed with df=v. 
|   Returns -1 if in error.   0.000002 < prob < 0.999998.
|
---------------------------------------------------------------------------------*/
MrBFlt PointChi2 (MrBFlt prob, MrBFlt v)
{
    MrBFlt      e = 0.5e-6, aa = 0.6931471805, p = prob, g,
                    xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
                    x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6,
                    tmp;

    if (p < 0.000002 || p > 0.999998 || v <= 0.0) 
        return (-1.0);
    g = LnGamma (v/2.0);
    xx = v/2.0;   
    c = xx - 1.0;
    if (v >= -1.24*log(p)) 
        goto l1;
    ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
    if (ch-e<0) 
        return (ch);
    goto l4;
    l1:
        if (v > 0.32) 
            goto l3;
        ch = 0.4;   
        a = log(1.0-p);
    l2:
        q = ch;  
        p1 = 1.0+ch*(4.67+ch);  
        p2 = ch*(6.73+ch*(6.66+ch));
        t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
        ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
        if (fabs(q/ch-1.0)-0.01 <= 0.0) 
            goto l4;
        else                       
            goto l2;
    l3: 
        x = PointNormal (p);
        p1 = 0.222222/v;   
        tmp = (x*sqrt(p1)+1.0-p1);
        ch = v*tmp*tmp*tmp;
        if (ch > 2.2*v+6.0)  
            ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
    l4:
        q = ch;   
        p1 = 0.5*ch;
        if ((t = IncompleteGamma (p1, xx, g)) < 0.0) 
            {
            MrBayesPrint ("%s   Error: Problem in PointChi2", spacer);
            return (-1.0);
            }
        p2 = p-t;
        t = p2*exp(xx*aa+g+p1-c*log(ch));   
        b = t/ch;  
        a = 0.5*t-b*c;
        s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
        s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
        s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
        s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
        s5 = (84.0+264.0*a+c*(175.0+606.0*a)) / 2520.0;
        s6 = (120.0+c*(346.0+127.0*c)) / 5040.0;
        ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
        if (fabs(q/ch-1.0) > e) 
            goto l4;
        return (ch);
}


/*---------------------------------------------------------------------------------
|
|   PointNormal
|
|   Returns z so That Prob{x<z} = prob where x ~ N(0,1) and
|   (1e-12) < prob < 1-(1e-12).  Returns (-9999) if in error.
|
|   Odeh, R. E. and J. O. Evans.  1974.  The percentage points of the normal
|     distribution.  Applied Statistics, 22:96-97 (AS70).
|
|   Newer methods:
|
|   Wichura, M. J.  1988.  Algorithm AS 241: The percentage points of the
|      normal distribution.  37:477-484.
|   Beasley, JD & S. G. Springer.  1977.  Algorithm AS 111: The percentage
|      points of the normal distribution.  26:118-121.
|
---------------------------------------------------------------------------------*/
MrBFlt PointNormal (MrBFlt prob)
{
    MrBFlt      a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
                    a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
                    b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
                    y, z = 0, p = prob, p1;

    p1 = (p<0.5 ? p : 1-p);
    if (p1<1e-20) 
       return (-9999);
    y = sqrt (log(1/(p1*p1)));   
    z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
    
    return (p<0.5 ? -z : z);
}


/*---------------------------------------------------------------------------------
|
|   PrintComplexVector
|
|   Prints a vector of dim complex numbers.
|
---------------------------------------------------------------------------------*/
void PrintComplexVector (int dim, complex *vec)
{
    int     i;

    MrBayesPrint ("{");
    for (i = 0; i < (dim - 1); i++) 
        {
        MrBayesPrint ("%lf + %lfi, ", vec[i].re, vec[i].im);
        if (i == 1) 
            MrBayesPrint("\n    ");
        }
    MrBayesPrint ("%lf + %lfi}\n", vec[dim - 1].re, vec[dim - 1].im);
}


/*---------------------------------------------------------------------------------
|
|   PrintSquareComplexMatrix
|
|   Prints a square matrix of complex numbers.
|
---------------------------------------------------------------------------------*/
void PrintSquareComplexMatrix (int dim, complex **m)
{
    int     row, col;

    MrBayesPrint ("{");
    for (row = 0; row < (dim - 1); row++) 
        {
        MrBayesPrint ("{");
        for (col = 0; col < (dim - 1); col++) 
            {
            MrBayesPrint ("%lf + %lfi, ", m[row][col].re, m[row][col].im);
            if (col == 1) 
                MrBayesPrint ("\n    ");
            }
        MrBayesPrint ("%lf + %lfi},\n", 
        m[row][dim - 1].re, m[row][dim - 1].im);
        }
    MrBayesPrint ("{");
    for (col = 0; col < (dim - 1); col++) 
        {
        MrBayesPrint ("%lf + %lfi, ", m[dim - 1][col].re, m[dim - 1][col].im);
        if (col == 1) 
            MrBayesPrint ("\n    ");
        }
    MrBayesPrint ("%lf + %lfi}}", m[dim - 1][dim - 1].re, m[dim - 1][dim - 1].im);
    MrBayesPrint ("\n");
}


/*---------------------------------------------------------------------------------
|
|   PrintSquareDoubleMatrix
|
|   Prints a square matrix of doubles.
|
---------------------------------------------------------------------------------*/
void PrintSquareDoubleMatrix (int dim, MrBFlt **matrix)
{
    int         i, j;
    
    for (i=0; i<dim; i++) 
        {
        for (j=0; j<dim; j++)
            MrBayesPrint ("%1.6lf ", matrix[i][j]);
        MrBayesPrint ("\n");
        }
}


/*---------------------------------------------------------------------------------
|
|   PrintSquareIntegerMatrix
|
|   Prints a square matrix of integers.
|
---------------------------------------------------------------------------------*/
void PrintSquareIntegerMatrix (int dim, int **matrix)
{
    int         i, j;
    
    for (i=0; i<dim; i++) 
        {
        for (j=0; j<dim; j++)
            MrBayesPrint ("%d ", matrix[i][j]);
        MrBayesPrint ("\n");
        }
}


/*---------------------------------------------------------------------------------
|
|   ProductOfRealAndComplex
|
|   Returns the complex product of a real and complex number.
|
---------------------------------------------------------------------------------*/
complex ProductOfRealAndComplex (MrBFlt a, complex b)
{
    complex     c;
    
    c.re = a * b.re;
    c.im = a * b.im;
    
    return (c);
}


/*---------------------------------------------------------------------------------
|
|   PsiExp: Returns psi (also called digamma) exponentiated
|       Algorithm from http://lib.stat.cmu.edu/apstat/103
|
---------------------------------------------------------------------------------*/
MrBFlt  PsiExp (MrBFlt alpha)
{
    MrBFlt      digamma, y, r, s, c, s3, s4, s5, d1;
    
    s = 1.0e-05;
    c = 8.5;
    s3 = 8.333333333333333333333333e-02;
    s4 = 8.333333333333333333333333e-03;
    s5 = 3.968253968e-03;
    d1 = -0.577215664901532860606512;   /* negative of Euler's constant */

    digamma = 0.0;
    y = alpha;
    if (y <= 0.0)
        return (0.0);
    
    if (y <= s)
        {
        digamma = d1 - 1.0 / y;
        return (exp (digamma));
        }
    
    while (y < c)
        {
        digamma -= 1.0 / y;
        y += 1.0;
        }

    r = 1.0 / y;
    digamma += (log (y) - 0.5 * r);
    r *= r;
    digamma -= r * (s3 - r * (s4 - r * s5));
    
    return (exp (digamma));
}


/*---------------------------------------------------------------------------------
|
|   PsiGammaLnProb: Calculates the log probability of a PsiGamma distributed
|      variable
|
---------------------------------------------------------------------------------*/
MrBFlt  PsiGammaLnProb (MrBFlt alpha, MrBFlt value)
{
    MrBFlt  beta, lnProb;

    beta = PsiExp (alpha);

    lnProb = alpha * log (beta) - LnGamma (alpha) + (alpha - 1.0) * log (value) - beta * value;

    return lnProb;
}


/*---------------------------------------------------------------------------------
|
|   PsiGammaLnRatio: Calculates the log prob ratio of two PsiGamma distributed
|      variables
|
---------------------------------------------------------------------------------*/
MrBFlt  PsiGammaLnRatio (MrBFlt alpha, MrBFlt numerator, MrBFlt denominator)
{
    MrBFlt beta, lnRatio;

    beta = PsiExp (alpha);

    lnRatio = (alpha - 1.0) * (log (numerator) - log (denominator)) - beta * (numerator - denominator);
    
    return (lnRatio);
}


/*---------------------------------------------------------------------------------
|
|   PsiGammaRandomVariable: Returns a random draw from the PsiGamma
|
---------------------------------------------------------------------------------*/
MrBFlt  PsiGammaRandomVariable (MrBFlt alpha, RandLong *seed)
{
    return GammaRandomVariable (alpha, PsiExp(alpha), seed);
}


/*---------------------------------------------------------------------------------
|
|   QuantileGamma
|
---------------------------------------------------------------------------------*/
MrBFlt QuantileGamma (MrBFlt x, MrBFlt alfa, MrBFlt beta)
{
    MrBFlt      quantile;

    quantile = POINTGAMMA(x, alfa, beta);
    
    return (quantile);
}


/*---------------------------------------------------------------------------------
|
|   RandomNumber
|
|   This pseudorandom number generator is described in:
|   Park, S. K. and K. W. Miller.  1988.  Random number generators: good
|      ones are hard to find.  Communications of the ACM, 31(10):1192-1201.
|
---------------------------------------------------------------------------------*/
MrBFlt RandomNumber (RandLong *seed)
{
    RandLong    lo, hi, test;

    hi = (*seed) / 127773;
    lo = (*seed) % 127773;
    test = 16807 * lo - 2836 * hi;
    if (test > 0)
        *seed = test;
    else
        *seed = test + 2147483647;
    return ((MrBFlt)(*seed) / (MrBFlt)2147483647);
}


/*---------------------------------------------------------------------------------
|
|   RndGamma
|
---------------------------------------------------------------------------------*/
MrBFlt RndGamma (MrBFlt s, RandLong *seed)
{
    MrBFlt  r=0.0;
    
    if (s <= 0.0)    
        puts ("Gamma parameter less than zero\n");

    else if (s < 1.0)  
        r = RndGamma1 (s, seed);
    else if (s > 1.0)  
        r = RndGamma2 (s, seed);
    else    /* 0-log() == -1 * log(), but =- looks confusing */
        r -= log(RandomNumber(seed));
        
    return (r);
}


/*---------------------------------------------------------------------------------
|
|   RndGamma1
|
---------------------------------------------------------------------------------*/
MrBFlt RndGamma1 (MrBFlt s, RandLong *seed)
{
    MrBFlt          r, x=0.0, tiny=1e-37, w;
    static MrBFlt   a, p, uf, ss=10.0, d;
    
    if (fabs(s-ss)>ETA) /* s != ss */ 
        {
        a  = 1.0 - s;
        p  = a / (a + s * exp(-a));
        uf = p * pow(tiny / a, s);
        d  = a * log(a);
        ss = s;
        }
    for (;;) 
        {
        r = RandomNumber(seed);
        if (r > p)        
            x = a - log((1.0 - r) / (1.0 - p)), w = a * log(x) - d;
        else if (r>uf)  
            x = a * pow(r / p, 1.0 / s), w = x;
        else            
            return (0.0);
        r = RandomNumber(seed);
        if (1.0 - r <= w && r > 0.0)
        if (r*(w + 1.0) >= 1.0 || -log(r) <= w)  
            continue;
        break;
        }
        
    return (x);
}


/*---------------------------------------------------------------------------------
|
|   RndGamma2
|
---------------------------------------------------------------------------------*/
MrBFlt RndGamma2 (MrBFlt s, RandLong *seed)
{
    MrBFlt          r , d, f, g, x;
    static MrBFlt   b, h, ss=0.0;
    
    if (fabs(s-ss)>ETA) /* s != ss */
        {
        b  = s - 1.0;
        h  = sqrt(3.0 * s - 0.75);
        ss = s;
        }
    for (;;) 
        {
        r = RandomNumber(seed);
        g = r - r * r;
        f = (r - 0.5) * h / sqrt(g);
        x = b + f;
        if (x <= 0.0) 
            continue;
        r = RandomNumber(seed);
        d = 64 * r * r * g * g * g;
        if (d * x < x - 2.0 * f * f || log(d) < 2.0 * (b * log(x / b) - f))  
            break;
        }
        
    return (x);
}


/*---------------------------------------------------------------------------------
|
|   SetQvalue
|
|   The Pade method for calculating the matrix exponential, tMat = e^{qMat * v}, 
|   has an error, e(p,q), that can be controlled by setting p and q to appropriate
|   values. The error is:
|
|      e(p,q) = 2^(3-(p+q)) * ((p!*q!) / (p+q)! * (p+q+1)!)
|
|   Setting p = q will minimize the error for a given amount of work. This function 
|   assumes that p = q. The function takes in as a parameter the desired tolerance
|   for the accuracy of the matrix exponentiation, and returns qV = p = q, that
|   will achieve the tolerance. The Pade approximation method is described in:
|  
|   Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
|      The Johns Hopkins University Press, Baltimore, Maryland.
|
|   The function is called from TiProbsUsingPadeApprox.
|
---------------------------------------------------------------------------------*/
int SetQvalue (MrBFlt tol)
{
    int         qV;
    MrBFlt      x;
    
    /*
    x = pow(2.0, 3.0 - (0 + 0)) * Factorial(0) * Factorial (0) / (Factorial(0+0) * Factorial (0+0+1));
    */
    x = 4.0;
    qV = 0;
    while (x > tol)
        {
        qV++;
        x = pow(2.0, 3.0 - (qV + qV)) * Factorial(qV) * Factorial (qV) / (Factorial(qV+qV) * Factorial (qV+qV+1));
        }
        
    return (qV);
}


/*---------------------------------------------------------------------------------
|
|   SetToIdentity
|
|   Make a dim X dim identity matrix.
|
---------------------------------------------------------------------------------*/
void SetToIdentity (int dim, MrBFlt **matrix)
{
    int         row, col;

    for (row=0; row<dim; row++)
        for (col=0; col<dim; col++)
            matrix[row][col] = (row == col ? 1.0 : 0.0);
}


/*---------------------------------------------------------------------------------
|
|   Tha
|
|   Calculate Owen's (1956) T(h,a) function, -inf <= h, a <= inf, 
|   where h = h1/h2, a = a1/a2, from the program of: 
|
|   Young, J. C. and C. E. Minder.  1974.  Algorithm AS 76.  An integral  
|      useful in calculating non-central t and bivariate normal  
|      probabilities.  Appl. Statist., 23:455-457.  [Correction: Appl.  
|      Statist., 28:113 (1979).  Remarks: Appl. Statist. 27:379 (1978),
|      28: 113 (1979), 34:100-101 (1985), 38:580-582 (1988)]  
|
|   See also: 
|
|   Johnson, N. L.  and S. Kotz.  1972.  Distributions in statistics:
|      multivariate distributions.  Wiley and Sons.  New York.  pp. 93-100.
|
---------------------------------------------------------------------------------*/
MrBFlt Tha (MrBFlt h1, MrBFlt h2, MrBFlt a1, MrBFlt a2)
{
    int             ng = 5, i;
    MrBFlt          U[] = {0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533},
                    R[] = {0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357},
                    pai2 = 6.283185307, tv1 = 1e-35, tv2 = 15.0, tv3 = 15.0, tv4 = 1e-5,
                    a, h, rt, t, x1, x2, r1, r2, s, k, sign = 1.0;

    if (fabs(h2) < tv1) 
        return (0.0);
    h = h1 / h2;
    if (fabs(a2) < tv1) 
        {
        t = CdfNormal(h);
        if (h >= 0.0) 
            t = (1.0 - t) / 2.0;
        else      
            t /= 2.0;
        return (t*(a1 >= 0.0 ? 1.0 : -1.0));
        }
    a = a1 / a2;
    if (a < 0.0) 
        sign = -1.0;  
    a = fabs(a);  
    h = fabs(h);   
    k = h*a;
    if (h > tv2 || a < tv1) 
        return (0.0);
    if (h < tv1) 
        return (atan(a)/pai2*sign);
    if (h < 0.3 && a > 7.0) /* (Boys RJ, 1989) */
        {             
        x1 = exp(-k*k/2.0)/k;
        x2 = (CdfNormal(k)-0.5)*sqrt(pai2);
        t = 0.25 - (x1+x2)/pai2*h + ((1.0+2.0/(k*k))*x1+x2)/(6.0*pai2)*h*h*h;
        return (MAX(t,0)*sign);
        }
    t = -h*h / 2.0;  
    x2 = a;  
    s = a*a;
    if (log(1.0+s)-t*s >= tv3) 
        {
        x1 = a/2;  
        s /= 4.0;
    for (;;) /* truncation point by Newton iteration */
        {        
        x2 = x1 + (t*s+tv3-log(s+1.0)) / (2.0*x1*(1.0/(s+1.0)-t));
        s = x2*x2;
        if (fabs(x2-x1) < tv4) 
            break;
        x1 = x2;
        }
    }
    for (i=0,rt=0; i<ng; i++) /* Gauss quadrature */
        {          
        r1 = 1.0+s*SQUARE(0.5+U[i]);
        r2 = 1.0+s*SQUARE(0.5-U[i]);
        rt+= R[i]*(exp(t*r1)/r1 + exp(t*r2)/r2);
        }
        
    return (MAX(rt*x2/pai2,0)*sign);
}


/*---------------------------------------------------------------------------------
|
|   TiProbsUsingEigens
|
---------------------------------------------------------------------------------*/
void TiProbsUsingEigens (int dim, MrBFlt *cijk, MrBFlt *eigenVals, MrBFlt v, MrBFlt r, MrBFlt **tMat, MrBFlt **fMat, MrBFlt **sMat)
{
    int             i, j, s;
    MrBFlt          sum, sumF, sumS, *ptr, EigValexp[192];

    for (s=0; s<dim; s++)
        EigValexp[s] = exp(eigenVals[s] * v * r);

    ptr = cijk;
    for (i=0; i<dim; i++)
        {
        for (j=0; j<dim; j++)
            {
            sum = 0.0;
            for (s=0; s<dim; s++)
                sum += (*ptr++) * EigValexp[s];
            tMat[i][j] = (sum < 0.0) ? 0.0 : sum;
            }
        }
        
#   if 0
    for (i=0; i<dim; i++)
        {
        sum = 0.0;
        for (j=0; j<dim; j++)
            {
            sum += tMat[i][j];
            }
        if (sum > 1.0001 || sum < 0.9999)
            {
            MrBayesPrint ("%s   Warning: Transition probabilities do not sum to 1.0 (%lf)\n", spacer, sum);
            }
        }
#   endif
    
    if (fMat != NULL && sMat != NULL)
        {
        ptr = cijk;
        for (i=0; i<dim; i++)
            {
            for (j=0; j<dim; j++)
                {
                sumF = sumS = 0.0;
                for (s=0; s<dim; s++)
                    {
                    sumF += (*ptr) * eigenVals[s] * r * EigValexp[s];
                    sumS += (*ptr++) * eigenVals[s] * eigenVals[s] * r * r * EigValexp[s];
                    }
                fMat[i][j] = sumF;
                sMat[i][j] = sumS;
                }
            }
        }
}


/*---------------------------------------------------------------------------------
|
|   TiProbsUsingPadeApprox
|
|   The method approximates the matrix exponential, tMat = e^{qMat * v}, using
|   the Pade approximation method, described in:
|  
|   Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
|      The Johns Hopkins University Press, Baltimore, Maryland.
|
|   The method approximates the matrix exponential with accuracy tol.
|
---------------------------------------------------------------------------------*/
void TiProbsUsingPadeApprox (int dim, MrBFlt **qMat, MrBFlt v, MrBFlt r, MrBFlt **tMat, MrBFlt **fMat, MrBFlt **sMat)
{
    int         qValue;
    MrBFlt      **a, tol;
    
    tol = 0.0000001;
    
    a = AllocateSquareDoubleMatrix (dim);
    
    MultiplyMatrixByScalar (dim, qMat, v * r, a);

    qValue = SetQvalue (tol);

    ComputeMatrixExponential (dim, a, qValue, tMat);
    
    FreeSquareDoubleMatrix (a);

    if (fMat != NULL && sMat != NULL)
        {
        MultiplyMatrices (dim, qMat, tMat, fMat);
        MultiplyMatrices (dim, qMat, fMat, sMat);
        }
}

