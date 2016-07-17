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
#include "command.h"
#include "mcmc.h"
#include "sumpt.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

typedef struct partctr
    {
    struct partctr  *left, *right;
    BitsLong        *partition;
    int             totCount;
    int             *count;
    MrBFlt          **length;
    MrBFlt          **height;
    MrBFlt          **age;
    int             ***nEvents; /* nEvents[0,nESets][0,numRuns][0,count[RunID]] */
    MrBFlt          ***bRate;   /* bRate  [0,nBSets][0,numRuns][0,count[RunID]] */
    MrBFlt          ***bLen;    /* bLen   [0,nBSets][0,numRuns][0,count[RunID]] */
    MrBFlt          **popSize;  /* popSize[0,numRuns][0,count[RunID]]           */
    }
    PartCtr;

typedef struct treectr
    {
    struct treectr  *left, *right;
    int             count;
    int             *order;
    }
    TreeCtr;

typedef struct
    {
    int     longestLineLength;
    int     numTreeBlocks;
    int     lastTreeBlockBegin;
    int     lastTreeBlockEnd;
    int     numTreesInLastBlock;
    }
    SumtFileInfo;

#define ALLOC_LEN               100     /* number of values to allocate each time in partition counter nodes */

#if defined (PRINT_RATEMUL_CPP)
FILE     *rateMultfp=NULL;
#endif

#undef  DEBUG_CONTREE

extern int inSumtCommand;
extern int inComparetreeCommand;
extern int DoUserTree (void);
extern int DoUserTreeParm (char *parmName, char *tkn);
extern int SafeFclose(FILE **);

extern int  SetUpPartitionCounters (void);
extern int  AddTreeToPartitionCounters (Tree *tree, int treeId, int runId);
extern void CalcTopoConvDiagn2 (int *nTrees);
extern void FreeChainMemory (void);

/* local prototypes */
int      CompareModelProbs (const void *x, const void *y);
int      PrintModelStats (char *fileName, char **headerNames, int nHeaders, ParameterSample *parameterSamples, int nRuns, int nSamples);
int      PrintOverlayPlot (MrBFlt **xVals, MrBFlt **yVals, int nRows, int startingFrom, int nSamples);
int      PrintParamStats (char *fileName, char **headerNames, int nHeaders, ParameterSample *parameterSamples, int nRuns, int nSamples);
void     PrintPlotHeader (void);

PartCtr *AddSumtPartition (PartCtr *r, PolyTree *t, PolyNode *p, int runId);
TreeCtr *AddSumtTree (TreeCtr *r, int *order);
PartCtr *AllocPartCtr (void);
TreeCtr *AllocTreeCtr (void);
void     CalculateTreeToTreeDistance (Tree *tree1, Tree *tree2, MrBFlt *d1, MrBFlt *d2, MrBFlt *d3);
int      ConTree (PartCtr **treeParts, int numTreeParts);
MrBFlt   CppEvolRate (PolyTree *t, PolyNode *p, int eSet);
int      ExamineSumtFile (char *fileName, SumtFileInfo *sumtFileInfo, char *treeName, int *brlensDef);
void     FreePartCtr (PartCtr *r);
void     FreeSumtParams (void);
void     FreeTreeCtr (TreeCtr *r);
int      Label (PolyNode *p, int addIndex, char *label, int maxLength);
int      OpenBrlensFile (int treeNo);
int      OpenComptFiles (void);
int      OpenSumtFiles (int treeNo);
void     PartCtrUppass (PartCtr *r, PartCtr **uppass, int *index);
int      PrintBrParamsToFile (PartCtr **treeParts, int numTreeParts, int treeNo, char *divString);
void     PrintConTree (FILE *fp, PolyTree *t);
void     PrintFigTreeConTree (FILE *fp, PolyTree *t, PartCtr **treeParts);
void     PrintFigTreeNodeInfo (FILE *fp, PartCtr *x, MrBFlt length);
void     PrintSumtTableLine (int numRuns, int *rowCount, Stat *theStats, MrBFlt *numPSRFSamples, MrBFlt *maxPSRF, MrBFlt *sumPSRF);
void     PrintSumtTaxaInfo (void);
void     Range (MrBFlt *vals, int nVals, MrBFlt *min, MrBFlt *max);
void     ResetTaxonSet (void);
int      ShowConPhylogram (FILE *fp, PolyTree *t, int screenWidth);
void     ShowSomeParts (FILE *fp, BitsLong *p, int offset, int nTaxaToShow);
void     SortPartCtr (PartCtr **item, int left, int right);
void     SortTerminalPartCtr (PartCtr **item, int len);
void     SortTreeCtr (TreeCtr **item, int left, int right);
int      StoreSumtTree (PackedTree *treeList, int index, PolyTree *t);
void     TreeCtrUppass (TreeCtr *r, TreeCtr **uppass, int *index);
int      TreeProb (void);
void     WriteConTree (PolyNode *p, FILE *fp, int showSupport);
void     WriteFigTreeConTree (PolyNode *p, FILE *fp, PartCtr **treeParts);

/* local (to this file) */
static int          numUniqueSplitsFound, numUniqueTreesFound, numPackedTrees[2], numAsterices;  /* length of local to this file variables */
static FILE        *fpParts=NULL, *fpTstat=NULL, *fpVstat, *fpCon=NULL, *fpTrees=NULL, *fpDists=NULL;  /* file pointers */
static PartCtr     *partCtrRoot = NULL;        /* binary tree for holding splits info      */
static TreeCtr     *treeCtrRoot = NULL;        /* binary tree for holding unique tree info */
static PackedTree  *packedTreeList[2];         /* list of trees in packed format           */


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
    int             i, n, nHeaders=0, numRows, numColumns, numRuns, whichIsX, whichIsY,
                    unreliable, oneUnreliable, burnin, longestHeader, len;
    MrBFlt          mean, harm_mean;
    char            **headerNames=NULL, temp[120];
    SumpFileInfo    fileInfo, firstFileInfo;
    ParameterSample *parameterSamples=NULL;
    FILE            *fpLstat=NULL;

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return NO_ERROR;
#   endif

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
        MrBayesPrintf (fpLstat, "[ID: %s]\n", stamp);

    /* print header */
    if (sumpParams.numRuns == 1)
        MrBayesPrintf (fpLstat, "arithmetic_mean\tharmonic_mean\tvalues_discarded\n");
    else
        MrBayesPrintf (fpLstat, "run\tarithmetic_mean\tharmonic_mean\tvalues_discarded\n");

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
            MrBayesPrintf (fpLstat, "%s\t", MbPrintNum(mean));
            MrBayesPrintf (fpLstat, "%s\t", MbPrintNum(harm_mean));
            if (unreliable == YES)
                MrBayesPrintf (fpLstat, "yes\n");
            else
                MrBayesPrintf (fpLstat, "no\n");
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
            MrBayesPrintf (fpLstat, "%d\t", n+1);
            MrBayesPrintf (fpLstat, "%s\t", MbPrintNum(mean));
            MrBayesPrintf (fpLstat, "%s\t", MbPrintNum(harm_mean));
            if (unreliable == YES)
                MrBayesPrintf (fpLstat, "yes\n");
            else
                MrBayesPrintf (fpLstat, "no\n");
            }                   
        }   /* next run */
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
        MrBayesPrintf (fpLstat, "all\t");
        MrBayesPrintf (fpLstat, "%s\t", MbPrintNum(mean));
        MrBayesPrintf (fpLstat, "%s\t", MbPrintNum(harm_mean));
        if (unreliable == YES)
            MrBayesPrintf (fpLstat, "yes\n");
        else
            MrBayesPrintf (fpLstat, "no\n");
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
    int             i, nHeaders=0, numRows, numColumns, numRuns, whichIsX, whichIsY,
                    longestHeader, len;
    char            **headerNames=NULL, temp[120];
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

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return NO_ERROR;
#   endif

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

    if (chainParams.burninSS > 0)
        {
        stepBeginSS = chainParams.burninSS + 1;
        }
    else
        {
        numSamplesInStepSS = (numRows-1)/(chainParams.numStepsSS-chainParams.burninSS);
        stepBeginSS = numSamplesInStepSS + 1;
        }

    numSamplesInStepSS = (numRows - stepBeginSS)/chainParams.numStepsSS;
    if ((numRows - stepBeginSS)%chainParams.numStepsSS!=0)
        {
        MrBayesPrint ("%s   Error:  Number of samples could not be evenly devided among steps (%d samples among %d steps). \n", spacer,(numRows - stepBeginSS),chainParams.numStepsSS);
        goto errorExit;
        }

    if (chainParams.relativeBurnin == YES)
        {
        stepBurnin = (int)(numSamplesInStepSS*chainParams.burninFraction);
        }
    else
        {
        stepBurnin = chainParams.chainBurnIn;
        if (stepBurnin >= numSamplesInStepSS)
            {
            MrBayesPrint ("%s   Error:  Burnin in each step(%d) is longer then the step itself(%d). \n", spacer,stepBurnin, numSamplesInStepSS);
            goto errorExit;               
            }
        }

    marginalLnLSS = (MrBFlt *) SafeCalloc (sumssParams.numRuns, sizeof(MrBFlt));
        /*Preparing and printing joined plot.*/
    plotArrayY = (MrBFlt **) SafeCalloc (sumssParams.numRuns+1, sizeof(MrBFlt*));
    for (i=0; i<sumssParams.numRuns+1; i++)
        plotArrayY[i] = (MrBFlt *) SafeCalloc (numSamplesInStepSS, sizeof(MrBFlt));

    plotArrayX = (MrBFlt **) SafeCalloc (sumssParams.numRuns, sizeof(MrBFlt*));
    for (i=0; i<sumssParams.numRuns; i++)
        {
        plotArrayX[i] = (MrBFlt *) SafeCalloc (numSamplesInStepSS, sizeof(MrBFlt));
        for (j=0; j<numSamplesInStepSS; j++)
            plotArrayX[i][j]=j+1;
        }

    MrBayesPrint ("%s   In total %d sampls are red from .p files.\n", spacer, numRows);
    MrBayesPrint ("\n");
    MrBayesPrint ("%s   Marginal likelihood (in natural log units) is estimated using stepping-stone sampling\n", spacer);
    MrBayesPrint ("%s   based on %d steps with %d samples within each step. \n", spacer, chainParams.numStepsSS, numSamplesInStepSS);
    MrBayesPrint ("%s   First %d samples (including generation 0) are discarded as initial burn-in.\n", spacer, stepBeginSS);
        if (chainParams.startFromPriorSS==YES)
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

    if (firstPass == YES && chainParams.relativeBurnin == YES)
        MrBayesPrint ("%s   The table entrances are based on samples excluding burn-in %d samples  (%d%%)    \n", spacer, stepBurnin,(int)(100*chainParams.burninFraction));
    else
        MrBayesPrint ("%s   The table entrances are based on samples excluding burn-in %d samples      \n", spacer, stepBurnin);
    MrBayesPrint ("%s   discarded at the begining of each step.  \n\n", spacer);

    //MrBayesPrint ("%s       Run   Marginal likelihood (ln)\n",spacer);
    //MrBayesPrint ("%s       ------------------------------\n",spacer);
    MrBayesPrint ("   Step");
    for (j=0; j<sumssParams.numRuns ; j++)
        {
        if (j<9)
            MrBayesPrint (" ");
        MrBayesPrint ("      run%d", j+1);
        }
    MrBayesPrint ("\n");
    for (i=0; i<sumssParams.numRuns; i++)
        {
        marginalLnLSS[i] = 0.0;  
        }
    for (stepIndexSS = chainParams.numStepsSS-1; stepIndexSS>=0; stepIndexSS--)   
        {
        if (chainParams.startFromPriorSS==YES)
            {
            stepLengthSS = BetaQuantile(chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-stepIndexSS)/(MrBFlt)chainParams.numStepsSS) - BetaQuantile(chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-1-stepIndexSS)/(MrBFlt)chainParams.numStepsSS);
            }
        else
            {
            stepLengthSS = BetaQuantile(chainParams.alphaSS, 1.0, (MrBFlt)(stepIndexSS+1)/(MrBFlt)chainParams.numStepsSS) - BetaQuantile(chainParams.alphaSS, 1.0, (MrBFlt)stepIndexSS/(MrBFlt)chainParams.numStepsSS);
            }
        MrBayesPrint ("   %3d   ", chainParams.numStepsSS-stepIndexSS);
        for (i=0; i<sumssParams.numRuns; i++)
            {
            lnlp = parameterSamples[whichIsY].values[i] + stepBeginSS + (chainParams.numStepsSS-stepIndexSS-1)*numSamplesInStepSS;
            nextSteplnlp = lnlp+numSamplesInStepSS;
            lnlp+= stepBurnin;
            stepAcumulatorSS = 0.0;
            stepScalerSS = *lnlp*stepLengthSS;
            while (lnlp < nextSteplnlp)
                {
               if (*lnlp*stepLengthSS > stepScalerSS + 200.0)
                    {
                    // adjust scaler;
                    stepAcumulatorSS /= exp(*lnlp*stepLengthSS - 100.0 - stepScalerSS); 
                    stepScalerSS= *lnlp*stepLengthSS - 100.0;
                    }
                stepAcumulatorSS += exp(*lnlp*stepLengthSS - stepScalerSS);
                lnlp++;
                }
            tmpMfl = (log(stepAcumulatorSS/(numSamplesInStepSS-stepBurnin)) + stepScalerSS);
            MrBayesPrint (" %10.3lf", tmpMfl);
            marginalLnLSS[i] += tmpMfl;
            }
        MrBayesPrint ("\n");
        //MrBayesPrint ("%s       %3d    %9.2f   \n", spacer, i+1, marginalLnLSS);
        }
    MrBayesPrint ("         ");
    for (j=0; j<sumssParams.numRuns ; j++)
        {
        if (j<9)
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

    if (firstPass == NO)
        goto sumssExitOptions;

    sumssStepPlot:

    MrBayesPrint ("\n\n%s   Step plot(s).\n",spacer);
    while (1)
        {
         MrBayesPrint ("\n");
        if (sumssParams.stepToPlot == 0)
            {
            beginPrint=(int)(sumssParams.discardFraction*stepBeginSS);
            countPrint=stepBeginSS-beginPrint;
            MrBayesPrint ("%s   Ploting step 0, i.e initial burn-in phase consisting of %d samples.\n", spacer,stepBeginSS);
            MrBayesPrint ("%s   According to 'Discardfrac=%.2f', first %d samples are not ploted.\n", spacer,sumssParams.discardFraction,beginPrint);
            }
        else
            {
            if (sumssParams.stepToPlot > chainParams.numStepsSS)
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

        if (sumssParams.askForMorePlots == NO || firstPass == YES)
            break;

        MrBayesPrint (" You can choose to print new step plots for different steps or discard fractions.\n");
        MrBayesPrint (" Allowed range of 'Steptoplot' are from 0 to %d.\n", chainParams.numStepsSS);
        MrBayesPrint (" If the next entered value is negative, 'sumss' will stop printing step plots.\n");
        MrBayesPrint (" If the next entered value is positive, but out of range, you will be offered\n");
        MrBayesPrint (" to change parameter 'Discardfrac' of 'sumss'.\n");
        MrBayesPrint (" Enter new step number 'Steptoplot':");
        k = scanf("%d",&j);
        if (j < 0)
            break;
        if (j > chainParams.numStepsSS)
            {
            do
                {
                MrBayesPrint (" Enter new value for 'Discardfrac', should be in range 0.0 to 1.0:");
                k = scanf("%f",&tmpf);
                sumssParams.discardFraction =  (MrBFlt)tmpf;
                }
            while (sumssParams.discardFraction < 0.0 || sumssParams.discardFraction > 1.0);
            }
        else
            sumssParams.stepToPlot=j;
    }

    if (firstPass == NO)
        goto sumssExitOptions;

    sumssJoinedPlot:

    MrBayesPrint ("\n\n%s   Joined plot(s).\n",spacer);
    while (1)
        {
        MrBayesPrint ("\n");
        MrBayesPrint ("%s   Joined plot of %d samples of all steps together. 'smoothing' is set to:%d\n", spacer,numSamplesInStepSS,sumssParams.smoothing);
        MrBayesPrint ("%s   According to step burn-in, first %d samples are not ploted.\n", spacer,stepBurnin);

        for (i=0; i<sumssParams.numRuns; i++)
            {
            for (j=stepBurnin;j<numSamplesInStepSS;j++)
                plotArrayY[sumssParams.numRuns][j]=0.0;
            lnlp= parameterSamples[whichIsY].values[i] + stepBeginSS;
            nextSteplnlp=lnlp;
            for (stepIndexSS = chainParams.numStepsSS-1; stepIndexSS>0; stepIndexSS--)
                {
                firstlnlp=plotArrayY[sumssParams.numRuns] + stepBurnin;
                lnlp+=stepBurnin;
                nextSteplnlp += numSamplesInStepSS;
                while (lnlp < nextSteplnlp)
                    {
                    *firstlnlp+=*lnlp;
                    firstlnlp++;
                    lnlp++;
                    }
                }
            for (j=stepBurnin;j<numSamplesInStepSS;j++)
                {
                sum=0.0;
                count=0;
                for (k=j-sumssParams.smoothing;k<=j+sumssParams.smoothing;k++)
                    {
                    if (k>=stepBurnin && k<numSamplesInStepSS)
                        {
                        sum += plotArrayY[sumssParams.numRuns][k];
                        count++;
                        }
                    }
                plotArrayY[i][j] = sum/count;
                /*
                if (max < plotArrayY[i][j])
                    max=plotArrayY[i][j];
                    */
                }
        /*  for (j=stepBurnin;j<numSamplesInStepSS;j++)
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

        if (sumssParams.askForMorePlots == NO || firstPass == YES)
            break;

        MrBayesPrint (" You can choose to print new joined plots with different step burn-in or smoothing.\n");
        MrBayesPrint (" Allowed range of step burn-in values are from 0 to %d.\n", numSamplesInStepSS-1);
        MrBayesPrint (" If the next entered value is negative, 'sumss' will stop printing joined plots.\n");
        MrBayesPrint (" If the next entered value is positive, but out of range, you will be offered\n");
        MrBayesPrint (" to change 'Smoothing'.\n");
        MrBayesPrint (" Enter new step burn-in:");
        k = scanf("%d",&j);
        if (j < 0)
            break;
        if (j >= numSamplesInStepSS)
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

    if (sumssParams.askForMorePlots == YES)
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
            }while (j<1 || j>4);

        if (j == 1)
            {
            MrBayesPrint (" Allowed range of step burn-in values are from 0 to %d\n", numSamplesInStepSS-1);
            MrBayesPrint (" Current step burn-in value is:%d\n", stepBurnin);
            MrBayesPrint (" Enter new step burn-in:");
            do
                {
                k = scanf("%d",&stepBurnin);
                }
            while (stepBurnin < 0 || stepBurnin > numSamplesInStepSS-1);
            MrBayesPrint ("\n"); 
            goto sumssTable;
            }
        else if (j == 2)
            {
            goto sumssStepPlot;
            }
        else if (j == 3)
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
    for (i=0; i<sumssParams.numRuns+1; i++)
        free(plotArrayY[i]);
    free(plotArrayY);
    for (i=0; i<sumssParams.numRuns; i++)
        free(plotArrayX[i]);
    free(plotArrayX);
    free(marginalLnLSS);
    
    return (NO_ERROR);
    
errorExit:

    /* free memory */
    FreeParameterSamples (parameterSamples);
    if (headerNames!=NULL)
        for (i=0; i<nHeaders; i++)
            free (headerNames[i]);
    free (headerNames);

    expecting = Expecting(COMMAND);    
    strcpy (spacer, "");
    chainParams.isSS=NO;
    if (plotArrayY!=NULL)
        for (i=0; i<sumssParams.numRuns+1; i++)
            free(plotArrayY[i]);
    free(plotArrayY);
    if (plotArrayX!=NULL)
        for (i=0; i<sumssParams.numRuns; i++)
            free(plotArrayX[i]);
    free(plotArrayX);
    free(marginalLnLSS);

    return (ERROR);
}


int DoSumpParm (char *parmName, char *tkn)
{
    int         tempI;
    MrBFlt      tempD;
    char        tempStr[100];

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
                if (strlen(tkn)>99 && (strchr(tkn,' ')-tkn) > 99)
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
                if (strlen(tkn)>99 && (strchr(tkn,' ')-tkn) > 99)
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
    int         tempI;
    MrBFlt      tempD;
    char        tempStr[100];

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
                if (strlen(tkn)>99 && (strchr(tkn,' ')-tkn) > 99)
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
                if (strlen(tkn)>99 && (strchr(tkn,' ')-tkn) > 99)
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
        /* set Smoothing (sumssParams.smoothing) *******************************************************/
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
                    MrBayesPrint ("%s   Setting sumss smoothing to %d\n", spacer, sumssParams.smoothing);
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
    s = (char *)SafeMalloc(2 * ((size_t)(fileInfo->longestLineLength) +10) * sizeof(char));
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
            if (GetToken (sumpToken, &tokenType, &sumpTokenP))
                goto errorExit;
            /* printf ("%s (%d)\n", sumpToken, tokenType); */
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
    if (chainParams.isSS == YES)
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
        if (fgets (s, fileInfo->longestLineLength + 2, fp)==NULL)
            goto errorExit;
    strcpy(headerLine, s);
    for (; lineNum < lastNonDigitLine+burnin; lineNum++)
        if (fgets (s, fileInfo->longestLineLength + 2, fp)==NULL)
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
            if (GetToken (sumpToken, &tokenType, &sumpTokenP))
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
                        MrBayesPrint ("%s   Number of columns is not even (%d in first line and %d in line %d of file %s)\n", spacer, firstNumCols, nNumbersOnThisLine, lineNum, fileName);
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
                MrBayesPrint ("%s   It could be that some parameter values are not numbers and the whole string containing \n",spacer); 
                MrBayesPrint ("%s   this wrongly formated parameter is treated as a header.\n",spacer);
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
    char        *s;

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
    int     i, j, len, longestHeader, *sampleCounts=NULL;
    char    temp[100];
    Stat    theStats;
    FILE    *fp;
    
    /* calculate longest header */
    longestHeader = 9;  /* length of 'parameter' */
    for (i=0; i<nHeaders; i++)
        {
        strcpy (temp, headerNames[i]);
        len = (int) strlen(temp);
        for (j=0; modelIndicatorParams[j][0]!='\0'; j++)
            if (IsSame (temp,modelIndicatorParams[j]) != DIFFERENT)
                break;
        if (modelIndicatorParams[j][0]!='\0')
            continue;
        if (!strcmp (temp, "Gen") || !strcmp (temp, "LnL") || !strcmp (temp, "LnPr"))
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
        MrBayesPrintf (fp, "[ID: %s]\n", stamp);

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
    MrBayesPrint ("\n");
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
        if (!strcmp (temp, "Gen") || !strcmp (temp, "LnL") || !strcmp (temp, "LnPr"))
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
    MrBFlt      f, *prob=NULL, *sum=NULL, *ssq=NULL, *min=NULL, *max=NULL, *stddev=NULL;
    char        temp[100];
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
        MrBayesPrintf (fp, "[ID: %s]\n", stamp);

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
        qsort((void *)elem, (size_t)nElements, sizeof(ModelProb), CompareModelProbs);

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
    int     i, j, k, k2, n, screenHeight, screenWidth, numY[60], width;
    char    plotSymbol[15][60];
    MrBFlt  x, y, minX, maxX, minY, maxY, meanY[60];
    
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
    if ((int)minX>0)
        width=(int)(log10(minX));
    else if ((int)minX==0)
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
    int     i, j, len, longestHeader, *sampleCounts=NULL;
    static char *temp=NULL;
    char    tempf[100];
    Stat    theStats;
    FILE    *fp;
    
    /* calculate longest header */
    longestHeader = 9;  /* length of 'parameter' */
    for (i=0; i<nHeaders; i++)
        {
        SafeStrcpy (&temp, headerNames[i]);
        len = (int) strlen(temp);
        for (j=0; modelIndicatorParams[j][0]!='\0'; j++)
            if (IsSame (temp,modelIndicatorParams[j]) != DIFFERENT)
                break;
        if (modelIndicatorParams[j][0]!='\0')
            continue;
        if (!strcmp (temp, "Gen") || !strcmp (temp, "LnL") || !strcmp (temp, "LnPr"))
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
        MrBayesPrintf (fp, "[ID: %s]\n", stamp);

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
    MrBayesPrint ("\n");
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
        if (!strcmp (temp, "Gen") || !strcmp (temp, "LnL") || !strcmp (temp, "LnPr"))
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

        if (theStats.minESS == theStats.minESS)
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
            if (theStats.minESS == theStats.minESS)
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
    int     i, j, k, numY[60], screenWidth, screenHeight;
    MrBFlt  x, y, minX, maxX, minY, maxY, meanY[60], diff;
                    
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
        if (meanY[i] == 0) /* with some compilers if (NaN < 1) is equal true !!! so we realy need this check*/
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
            if (GetToken (sumpToken, &tokenType, &p))
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


/* the following are moved from sumt.c to combine with sump.c */
PartCtr *AddSumtPartition (PartCtr *r, PolyTree *t, PolyNode *p, int runId)
{
    int     i, n, comp, nLongsNeeded = sumtParams.BitsLongsNeeded;
    
    if (r == NULL)
        {
        /* new partition */
        /* create a new node */
        r = AllocPartCtr ();
        if (r == NULL)
            return NULL;
        numUniqueSplitsFound++;
        for (i=0; i<nLongsNeeded; i++)
            r->partition[i] = p->partition[i];
        for (i=0; i<sumtParams.numRuns; i++)
            r->count[i] = 0;
        r->left = r->right = NULL;
        /* record values */
        if (sumtParams.brlensDef == YES)
            r->length[runId][0]= p->length;
        if (sumtParams.isClock == YES)
            r->height[runId][0]= p->depth;
        if (sumtParams.isCalibrated == YES)
            r->age[runId][0]= p->age;
        for (i=0; i<sumtParams.nESets; i++)
            r->nEvents[i][runId][0] = t->nEvents[i][p->index];
        for (i=0; i<sumtParams.nBSets; i++)
            {
            r->bLen [i][runId][0] = t->effectiveBrLen[i][p->index];
            r->bRate[i][runId][0] = t->effectiveBrLen[i][p->index] / p->length;
            }
        if (t->popSizeSet == YES)
            r->popSize[runId][0] = t->popSize[p->index];
        r->count[runId] ++;
        r->totCount++;
        }
    else
        {
        for (i=0; i<nLongsNeeded; i++)
            {
            if (r->partition[i] != p->partition[i])
                break;
            }
        
        if (i == nLongsNeeded)
            comp = 0;
        else if (r->partition[i] < p->partition[i])
            comp = -1;
        else
            comp = 1;
        
        if (comp == 0)          /* repeated partition */
            {
            n = r->count[runId];
            /* check if we need to allocate more space */
            if (n % ALLOC_LEN == 0)
                {
                /* allocate more space */
                if (sumtParams.brlensDef == YES)
                    r->length[runId] = (MrBFlt *) SafeRealloc ((void *)r->length[runId], ((size_t)n+ALLOC_LEN)*sizeof(MrBFlt));
                if (sumtParams.isClock == YES)
                    r->height[runId] = (MrBFlt *) SafeRealloc ((void *)r->height[runId], ((size_t)n+ALLOC_LEN)*sizeof(MrBFlt));
                if (sumtParams.isCalibrated == YES)
                    r->age[runId] = (MrBFlt *) SafeRealloc ((void *)r->age[runId], ((size_t)n+ALLOC_LEN)*sizeof(MrBFlt));
                if (sumtParams.nESets > 0)
                    {
                    for (i=0; i<sumtParams.nESets; i++)
                        r->nEvents[i][runId] = (int *) SafeRealloc ((void *)r->nEvents[i][runId], ((size_t)n+ALLOC_LEN)*sizeof(int));
                    }
                if (sumtParams.nBSets > 0)
                    {
                    for (i=0; i<sumtParams.nBSets; i++)
                        {
                        r->bRate[i][runId]   = (MrBFlt *) SafeRealloc ((void *)r->bRate[i][runId], ((size_t)n+ALLOC_LEN)*sizeof(MrBFlt));
                        r->bLen [i][runId]   = (MrBFlt *) SafeRealloc ((void *)r->bLen [i][runId], ((size_t)n+ALLOC_LEN)*sizeof(MrBFlt));
                        }
                    }
                if (sumtParams.popSizeSet == YES)
                    r->popSize[runId] = (MrBFlt *) SafeRealloc ((void *)r->popSize[runId], ((size_t)n+ALLOC_LEN)*sizeof(MrBFlt));
                }
            /* record values */
            r->count[runId]++;
            r->totCount++;
            if (sumtParams.brlensDef == YES)
                r->length[runId][n]= p->length;
            if (sumtParams.isClock == YES)
                r->height[runId][n]= p->depth;
            if (sumtParams.isCalibrated == YES)
                r->age[runId][n]= p->age;
            if (sumtParams.nESets > 0)
                {
                for (i=0; i<sumtParams.nESets; i++)
                    r->nEvents[i][runId][n] = t->nEvents[i][p->index];
                }
            if (sumtParams.nBSets > 0)
                {
                for (i=0; i<sumtParams.nBSets; i++)
                    {
                    r->bLen [i][runId][n]   = t->effectiveBrLen[i][p->index];
                    r->bRate[i][runId][n]   = t->effectiveBrLen[i][p->index] / p->length;
                    }
                }
            if (sumtParams.popSizeSet == YES)
                r->popSize[runId][n] = t->popSize[p->index];
            }
        else if (comp < 0)      /* greater than -> into left subtree */
            {
            if ((r->left = AddSumtPartition (r->left, t, p, runId)) == NULL)
                {
                FreePartCtr (r);
                return NULL;
                }
            }
        else
            {
            /* smaller than -> into right subtree */
            if ((r->right = AddSumtPartition (r->right, t, p, runId)) == NULL)
                {
                FreePartCtr (r);
                return NULL;
                }
            }
        }

    return r;
}


TreeCtr *AddSumtTree (TreeCtr *r, int *order)
{
    int     i, comp;

    if (r == NULL)
        {
        /* new tree */
        /* create a new node */
        r = AllocTreeCtr();
        if (!r)
            return NULL;
        numUniqueTreesFound++;
        for (i=0; i<sumtParams.orderLen; i++)
            r->order[i] = order[i];
        r->count = 1;
        }
    else
        {
        for (i=0; i<sumtParams.orderLen; i++)
            if (r->order[i] != order[i])
                break;
        
        if (i==sumtParams.orderLen)
            comp = 0;
        else if (order[i] < r->order[i])
            comp = 1;
        else
            comp = -1;
        
        if (comp == 0)          /* repeated partition */
            r->count++;
        else if (comp < 0)      /* greater than -> into left subtree */
            {
            if ((r->left = AddSumtTree (r->left, order)) == NULL)
                {
                FreeTreeCtr (r);
                return NULL;
                }
            }
        else
            {
            /* smaller than -> into right subtree */
            if ((r->right = AddSumtTree (r->right, order)) == NULL)
                {
                FreeTreeCtr (r);
                return NULL;
                }
            }
        }

    return r;
}


/* AllocPartCtr: Allocate space for one partition counter node using info in sumtParams */
PartCtr *AllocPartCtr ()
{
    int             i, j;
    PartCtr         *r;
    
    /* allocate basic stuff */
    r = (PartCtr *) SafeCalloc (1, sizeof(PartCtr));
    r->left = r->right = NULL;
    r->partition = (BitsLong *) SafeCalloc ((size_t)(sumtParams.BitsLongsNeeded), sizeof(BitsLong));
    r->count = (int *) SafeCalloc ((size_t)(sumtParams.numRuns), sizeof (int));
    if (sumtParams.brlensDef)
        {
        r->length = (MrBFlt **) SafeCalloc ((size_t)(sumtParams.numRuns), sizeof (MrBFlt *));
        for (i=0; i<sumtParams.numRuns; i++)
            r->length[i] = (MrBFlt *) SafeCalloc ((size_t)ALLOC_LEN, sizeof(MrBFlt));
        }
    if (sumtParams.isClock)
        {
        r->height = (MrBFlt **) SafeCalloc ((size_t)(sumtParams.numRuns), sizeof (MrBFlt *));
        for (i=0; i<sumtParams.numRuns; i++)
            r->height[i] = (MrBFlt *) SafeCalloc ((size_t)ALLOC_LEN, sizeof(MrBFlt));
        r->age = (MrBFlt **) SafeCalloc ((size_t)(sumtParams.numRuns), sizeof (MrBFlt *));
        for (i=0; i<sumtParams.numRuns; i++)
            r->age[i] = (MrBFlt *) SafeCalloc ((size_t)ALLOC_LEN, sizeof(MrBFlt));
        }

    /* allocate relaxed clock parameters: eRate, nEvents, bRate */
    if (sumtParams.nESets > 0)
        r->nEvents = (int    ***) SafeCalloc ((size_t)(sumtParams.nESets), sizeof(int **));
    for (i=0; i<sumtParams.nESets; i++)
        {
        r->nEvents[i] = (int    **) SafeCalloc ((size_t)(sumtParams.numRuns), sizeof(int *));
        for (j=0; j<sumtParams.numRuns; j++)
            r->nEvents[i][j] = (int    *) SafeCalloc ((size_t) ALLOC_LEN, sizeof(int));
        }
    if (sumtParams.nBSets > 0)
        {
        r->bLen  = (MrBFlt ***) SafeCalloc ((size_t)(sumtParams.nBSets), sizeof(MrBFlt **));
        r->bRate = (MrBFlt ***) SafeCalloc ((size_t)(sumtParams.nBSets), sizeof(MrBFlt **));
        }
    for (i=0; i<sumtParams.nBSets; i++)
        {
        r->bLen[i]    = (MrBFlt **) SafeCalloc ((size_t)(sumtParams.numRuns), sizeof(MrBFlt *));
        r->bRate[i]   = (MrBFlt **) SafeCalloc ((size_t)(sumtParams.numRuns), sizeof(MrBFlt *));
        for (j=0; j<sumtParams.numRuns; j++)
            {
            r->bLen[i][j]    = (MrBFlt *) SafeCalloc ((size_t)ALLOC_LEN, sizeof(MrBFlt));
            r->bRate[i][j]   = (MrBFlt *) SafeCalloc ((size_t)ALLOC_LEN, sizeof(MrBFlt));
            }
        }
    if (sumtParams.popSizeSet == YES)
        {
        r->popSize = (MrBFlt **) SafeCalloc ((size_t)(sumtParams.numRuns), sizeof (MrBFlt *));
        for (i=0; i<sumtParams.numRuns; i++)
            r->popSize[i] = (MrBFlt *) SafeCalloc ((size_t)ALLOC_LEN, sizeof(MrBFlt));
        }

    return r;
}


/* AllocTreeCtr: Allocate space for a tree counter node using info in sumtParams struct*/
TreeCtr *AllocTreeCtr ()
{
    TreeCtr     *r;

    r = (TreeCtr *) SafeCalloc (1, sizeof(TreeCtr));
    
    r->left = r->right = NULL;
    
    r->order = (int *) SafeCalloc ((size_t)(sumtParams.orderLen), sizeof(int));

    return r;
}


void CalculateTreeToTreeDistance (Tree *tree1, Tree *tree2, MrBFlt *d1, MrBFlt *d2, MrBFlt *d3)
{
    int         i, j, k;
    MrBFlt      treeLen1=0.0, treeLen2=0.0;
    TreeNode    *p, *q=NULL;
    
    (*d1) = (*d2) = (*d3) = 0.0;

    /* set distance-based measures to max value */
    if (sumtParams.brlensDef == YES)
        {
        treeLen1 = TreeLen(tree1);
        treeLen2 = TreeLen(tree2);
        (*d2) = treeLen1 + treeLen2;
        (*d3) = 2.0;
        }

    /* now we can get distances in a single pass */
    for (i=0; i<tree1->nNodes; i++)
        {
        p = tree1->allDownPass[i];
        for (j=0; j<tree2->nNodes; j++)
            {
            q = tree2->allDownPass[j];
            for (k=0; k<sumtParams.BitsLongsNeeded; k++)
                if (p->partition[k] != q->partition[k])
                    break;
            if (k == sumtParams.BitsLongsNeeded)
                break;
            }
        if (j < tree2->nNodes)
            {
            /* match */
            if (sumtParams.brlensDef == YES)
                {
                (*d2) -= (p->length + q->length - fabs(p->length - q->length));
                (*d3) -= (p->length/treeLen1 + q->length/treeLen2 - fabs(p->length/treeLen1 - q->length/treeLen2));
                }
            }
        else /* if (k < sumtParams.BitsLongsNeeded) */
            {
            /* no match */
            (*d1) += 2.0;
            }
        }

#   if 0        
    printf ("DISTANCES: %lf %lf %lf (%lf %lf)\n", *d1, *d2, *d3, tl1, tl2);
    for (i=0; i<nnds; i++)
        {
        printf ("%4d -- %4d (%lf) %4d (%lf)\n", i, list1[i], lengths1[i], list2[i], lengths2[i]);
        }
#   endif
}


/* ConTree: Construct consensus tree FIXME: numTreeParts is not used*/
int ConTree (PartCtr **treeParts, int numTreeParts)
{
    int         i, j, targetNode, nBits, isCompat, numTerminalsEncountered;
    BitsLong    x, *partition = NULL, bitsLongOne;
    MrBFlt      freq, freqInterapted=0;
    PolyTree    *t, *t2=NULL;
    PolyNode    *p, *q, *r, *ql, *pl;
    PartCtr     *part;
    Stat        theStats;
    int         isFirstLoop=1, isInterapted=0;
    
    /* we use a BitsLong variable set to one to escape dependence on interpretation of constant (signed?) int/long '1L' */
    bitsLongOne = 1;
    
    /* check that we have at least three species */
    if (sumtParams.numTaxa < 3)
        {
        MrBayesPrint ("%s   Too few taxa included to show consensus trees\n", spacer);
        return ERROR;
        }
    
treeConstruction:
    /* now, make a consensus tree */
    /* first allocate and initialize consensus tree */
    t = AllocatePolyTree(sumtParams.numTaxa);
    if (!t)
        {
        MrBayesPrint ("%s   Could not allocate consensus tree\n", spacer);
        return(ERROR);
        }
    t->isRooted = sumtParams.isRooted;
    t->isClock = sumtParams.isClock;
    t->isRelaxed = sumtParams.isRelaxed;

    /* initialize consensus tree nodes */
    for (i=0; i<sumtParams.numTaxa; i++)
        {
        t->nodes[i].left = NULL;
        t->nodes[i].sib = NULL;
        t->nodes[i].index = i;
        t->nodes[i].partitionIndex = -1;     /* partition ID */
        t->nodes[i].age = 0.0;              /* temporally set to minimum value to allow any insertion in front of the terminal before actual
                                               values of age and depth are available */
        strcpy(t->nodes[i].label, sumtParams.taxaNames[i]);
        t->nodes[i].depth = 0.0;
        }
    for (; i<t->memNodes; i++)
        {
        t->nodes[i].left = NULL;
        t->nodes[i].sib = NULL;
        t->nodes[i].index = i;
        t->nodes[i].partitionIndex = -1;     /* partition ID */
        strcpy (t->nodes[i].label, "");
        }

    /* create bush 
       ->x counts number of subtended terminals 
       make sure t->root->left is in outgroup */
    p = t->root = &t->nodes[sumtParams.numTaxa];
    p->anc = p->sib = NULL;
    p->x = sumtParams.numTaxa;
    p->age = MRBFLT_MAX; /* temporarily set to maximum value to allow any insertion in front of the root before actual values of age and depth are available */
    p->depth = MRBFLT_MAX;
    j = localOutGroup;
    q = &t->nodes[j];
    p->left = q;
    q->anc = p;
    q->x = 1;
    for (i=0; i<sumtParams.numTaxa; i++)
        {
        if (i != j)
            {
            q->sib = &t->nodes[i];
            q = q->sib;
            q->anc = p;
            q->x = 1;
            }
        }
    q->sib = NULL;

    /* Resolve bush according to partitions.
       Partitions may include incompatible ones.
       Partitions must be sorted from most frequent to least frequent 
       for quit test to work when a 50% majority rule tree is requested
       and in general for consensus tree to be correct. */
    t->nNodes = sumtParams.numTaxa + 1;
    t->nIntNodes = 1;
    if (sumtParams.isRooted == YES)
        targetNode = 2 * sumtParams.numTaxa - 2;
    else
        targetNode = 2 * sumtParams.numTaxa - 3;

    numTerminalsEncountered = 0;
    for (i=0; i<numUniqueSplitsFound; i++)
        {
        /* get partition */
        part = treeParts[i];

        /* calculate frequency and test if time to quit */
        if (t->nNodes > targetNode && numTerminalsEncountered == sumtParams.numTaxa)
            break;
        freq = (MrBFlt)(part->totCount) / (MrBFlt)(sumtParams.numTreesSampled);
        if (freq < 0.50 && !strcmp(sumtParams.sumtConType, "Halfcompat"))
            break;
        
        /* get partition */
        partition = part->partition;

        /* count bits in this partition */
        for (j=nBits=0; j<sumtParams.BitsLongsNeeded; j++)
            {
            x = partition[j];
            for (x = partition[j]; x != 0; x &= (x - 1))
                nBits++;
            }

        /* find out if this is an informative partition */
        if (nBits == sumtParams.numTaxa || nBits == 0)
            {
            /* this is the root (for setting age of root node when tree is dated) */
            q = t->root;
            q->partitionIndex = i;
            if (sumtParams.isClock == YES)
                {
                GetSummary(part->height, sumtParams.numRuns, part->count, &theStats, sumtParams.HPD);
                q->depth = theStats.median;
                for (p = q->left; p!=NULL; p = p->sib)
                    {
                    if (q->depth <= p->depth)
                        break;
                    }
                assert (p==NULL); /* Root always has 100% freq and it should be older than any other node that has 100% freq. */
                }
            if (sumtParams.isCalibrated == YES)
                {
                GetSummary(part->age, sumtParams.numRuns, part->count, &theStats, sumtParams.HPD);
                q->age = theStats.median;
                for (p = q->left; p!=NULL; p = p->sib)
                    {
                    if (q->age <= p->age)
                        break;
                    }
                assert (p==NULL); /* Root always has 100% freq and it should be older than any other node that has 100% freq. */
                }
            }
        else if (nBits > 1 && !(nBits == sumtParams.numTaxa - 1 && sumtParams.isRooted == NO))
            {
            /* this is an informative partition */
            /* find anc of partition */
            j = FirstTaxonInPartition (partition, sumtParams.BitsLongsNeeded);
            for (p = &t->nodes[j]; p!=NULL; p = p->anc)
                if (p->x > nBits)
                    break;
                    
            /* do not include if incompatible with ancestor or any of descendants
               do not check terminals or root because it is
               redundant and partitions have not necessarily been set for those */
            isCompat = YES;
            if (p->anc != NULL && IsPartNested(partition, p->partition, sumtParams.BitsLongsNeeded)==NO)
                isCompat = NO;
            else 
                {
                for (q=p->left; q!=NULL; q=q->sib)
                    {
                    if (q->x > 1 && IsPartCompatible(q->partition, partition, sumtParams.BitsLongsNeeded)==NO)
                        break;
                    }
                if (q!=NULL)
                    isCompat = NO;
                }

            if (isCompat == NO)
                continue;

            /* set new node */
            q = &t->nodes[t->nNodes];
            q->partitionIndex = i;
            q->x = nBits;
            q->partition = partition;
            GetSummary(part->length, sumtParams.numRuns, part->count, &theStats, sumtParams.HPD);
            q->support = freq;
            q->length = theStats.median;
            r=NULL;
            if (sumtParams.isClock == YES)
                {
                GetSummary(part->height, sumtParams.numRuns, part->count, &theStats, sumtParams.HPD);
                q->depth = theStats.median;
                if (freq < 1.00)
                    {
                    for (r = p->left; r!=NULL; r = r->sib)
                        {
                        if (IsPartNested(r->partition, partition, sumtParams.BitsLongsNeeded) &&  r->depth >= q->depth)
                            break; /* child is older then the node we try to add. Not good.*/
                        }
                    if (p->depth <= q->depth)
                        {      /* New node older than the parent. Not good.*/
                        r = p; /* Just to make r!=NULL*/
                        }
                    }
                }
            if (sumtParams.isCalibrated == YES)
                {
                GetSummary(part->age, sumtParams.numRuns, part->count, &theStats, sumtParams.HPD);
                q->age = theStats.median;
                if (freq < 1.00)
                    {
                    for (r = p->left; r!=NULL; r = r->sib)
                        {
                        if (freq < 1.00 && IsPartNested(r->partition, partition, sumtParams.BitsLongsNeeded) && r->age >= q->age)
                            break; /* child is older then the node we try to add. Not good.*/
                        }
                    if (p->age <= q->age)
                        {      /* New node older than the parent. Not good.*/
                        r = p; /* Just to make r!=NULL*/
                        }
                    }
                }

            if (r!=NULL && isFirstLoop)
                {
                 /* cancel the addition of the new node*/
                isInterapted =1;
                freqInterapted=freq;
                break; /* Finish creating the polytree */
                }
            t->nNodes++;
            t->nIntNodes++;

            /* go through descendants of anc */
            ql = pl = NULL;
            for (r=p->left; r!=NULL; r=r ->sib)
                {
                /* test if r is in the new partition or not */
                if ((r->x > 1 && IsPartNested(r->partition, partition, sumtParams.BitsLongsNeeded)) || (r->x == 1 && (partition[r->index / nBitsInALong] & (bitsLongOne << (r->index % nBitsInALong))) != 0))
                    {
                    /* r is in the partition */
                    if (ql == NULL)
                        q->left = r;
                    else
                        ql->sib = r;
                    ql = r;
                    r->anc = q;
                    }
                else
                    {
                    /* r is not in the partition */
                    if (pl == NULL)
                        p->left = r;
                    else
                        pl->sib = r;
                    pl = r;
                    }
                }
            /* terminate new sib-node chain */
            ql->sib = NULL;
            /* new node is last in old sib-node chain */
            pl->sib = q;
            q->sib = NULL;
            q->anc = p;
            }
        else
            /* singleton partition */
            {
            if (nBits == sumtParams.numTaxa - 1)
                j = localOutGroup;
            else
                j = FirstTaxonInPartition(partition, sumtParams.BitsLongsNeeded); /* nbits == 1 */
            q = &t->nodes[j];
            q->partitionIndex = i;
            q->partition = partition;
            numTerminalsEncountered++;
            GetSummary(part->length, sumtParams.numRuns, part->count, &theStats, sumtParams.HPD);
            q->length = theStats.median;
            if (sumtParams.isClock == YES)
                {
                GetSummary(part->height, sumtParams.numRuns, part->count, &theStats, sumtParams.HPD);
                q->depth = theStats.median;
                if (q->anc->depth < q->depth)
                    {
                    /* We never should get here because terminals always have 100% freq and they are younger than any other node that has 100% freq. */
                    /* We should be careful with the trees with 0-brl generated under the fossilized birth-death prior! (<= is changed to <) */
                    assert (0);
                    }
                }
            if (sumtParams.isCalibrated == YES)
                {
                GetSummary(part->age, sumtParams.numRuns, part->count, &theStats, sumtParams.HPD);
                q->age = theStats.median;
                if (q->anc->age < q->age)
                    {
                    /* We never should get here because terminals always have 100% freq and they are younger than any other node that has 100% freq. */
                    /* We should be careful with the trees with 0-brl generated under the fossilized birth-death prior! (<= is changed to <) */
                    assert (0);
                    }
                }
            }
        }

    if (isFirstLoop)
        {
        t2 = t;
        if (isInterapted)
            {
            isFirstLoop = 0;
            goto treeConstruction;
            }
        }

    /* get downpass arrays */
    GetPolyDownPass(t);

    /* order tips */
    if (sumtParams.orderTaxa == YES)
        OrderTips (t);

    if (t!=t2)
        {
        /* get downpass arrays */
        GetPolyDownPass(t2);

        /* order tips */
        if (sumtParams.orderTaxa == YES)
            OrderTips (t2);
        }
        
    /* draw tree to stdout and fp */
    MrBayesPrint ("\n%s   Clade credibility values:\n\n", spacer);
    ShowConTree (stdout, t, 80, YES);
    if (logToFile == YES)
        ShowConTree (logFileFp, t, 80, YES);
    if (sumtParams.brlensDef == YES)
        {
        MrBayesPrint ("\n");
        if (sumtParams.isClock == YES)
            MrBayesPrint ("%s   Phylogram (based on median node depths):\n", spacer);
        else
            MrBayesPrint ("%s   Phylogram (based on average branch lengths):\n", spacer);
        if (isInterapted)
            {
            MrBayesPrint ("%s   Warning. Phylogram containing all nodes with credibility values exceeding\n",spacer);
            MrBayesPrint ("%s   the level set by Contype could not be constructed.\n",spacer);
            MrBayesPrint ("%s   Only nodes with credibility values exceeding %.2f%% (percentage of trees\n", spacer, freqInterapted*100);
            MrBayesPrint ("%s   where the node is present) were included in the phylogram.\n", spacer);
            }
        MrBayesPrint ("\n");
        ShowConPhylogram (stdout, t2, 80);
        if (logToFile == YES)
            ShowConPhylogram (logFileFp, t2, 80);
        }

    /* print taxa block */
    MrBayesPrintf (fpCon, "begin taxa;\n");
    MrBayesPrintf (fpCon, "\tdimensions ntax=%d;\n", sumtParams.numTaxa);
    MrBayesPrintf (fpCon, "\ttaxlabels\n", sumtParams.numTaxa);
    for (i=0; i<sumtParams.numTaxa; i++)
        {
        for (j=0; j<t2->nNodes; j++)
            if (t2->nodes[j].index == i)
                break;
        MrBayesPrintf (fpCon, "\t\t%s\n", t2->nodes[j].label);
        }
    MrBayesPrintf (fpCon, "\t\t;\nend;\n");
    
    MrBayesPrintf (fpCon, "begin trees;\n");
    MrBayesPrintf (fpCon, "\ttranslate\n");
    for (i=0; i<sumtParams.numTaxa; i++)
        {
        for (j=0; j<t2->nNodes; j++)
            if (t2->nodes[j].index == i)
                break;
        if (i == sumtParams.numTaxa-1)
            MrBayesPrintf (fpCon, "\t\t%d\t%s\n", t2->nodes[i].index+1, t2->nodes[i].label);
        else
            MrBayesPrintf (fpCon, "\t\t%d\t%s,\n", t2->nodes[i].index+1, t2->nodes[i].label);
        }
    MrBayesPrintf (fpCon, "\t\t;\n");
    if (sumtParams.consensusFormat == SIMPLE)
        PrintConTree(fpCon, t2);
    else if (sumtParams.consensusFormat == FIGTREE)
        PrintFigTreeConTree(fpCon, t2, treeParts);
    MrBayesPrintf (fpCon, "end;\n");

    if (t!=t2)
        {
        FreePolyTree (t2);
        }
    /* free memory */
    FreePolyTree (t);

    return (NO_ERROR);
    MrBayesPrint ("%d", numTreeParts); /* just because I am tired of seeing the unused parameter error msg */
}


MrBFlt CppEvolRate (PolyTree *t, PolyNode *p, int eSet)
{
    int         i, nEvents;
    MrBFlt      ancRate, branchRate, *rate, *pos;
    PolyNode    *q;

    nEvents = t->nEvents[eSet][p->index];
    pos = t->position[eSet][p->index];
    rate = t->rateMult[eSet][p->index];

    /* note that event positions are from top of branch (more recent, descendant tip) */
    ancRate = 1.0;
//  if (t->eType[eSet] == CPPm)
//      {
    for (q=p; q->anc != NULL; q=q->anc)
        {
        for (i=0; i<t->nEvents[eSet][p->index]; i++)
            ancRate *= t->rateMult[eSet][p->index][i];
        }
    if (nEvents > 0)
        {
        branchRate = rate[0] * pos[0];
        for (i=1; i<nEvents; i++)
            {
            branchRate += (pos[i] - pos[i-1]);
            branchRate *= rate[i];
            }
        branchRate += 1.0 - pos[nEvents-1];
        branchRate *= ancRate;
        }
    else
        branchRate = ancRate;
//      }
/*
    else if (t->eType[eSet] == CPPi)
        {
        for (q=p; q->anc != NULL; q=q->anc)
            {
            if (t->nEvents[eSet][p->index]>0)
                {
                ancRate = t->rateMult[eSet][p->index][0];
                break;
                }
            }
        if (nEvents > 0)
            {
            branchRate = ancRate * (1.0 - pos[nEvents-1]);
            for (i=nEvents-2; i>=0; i--)
                {
                branchRate += (rate[i+1] * (pos[i+1] - pos[i]));
                }
            branchRate += (rate[0] * pos[0]);
            }
        else
            branchRate = ancRate;
        }
*/

    return branchRate;
}


int DoCompareTree (void)
{
    int             i, j, k, n, longestLineLength, brlensDef[2], numTreesInLastBlock[2],
                    lastTreeBlockBegin[2], lastTreeBlockEnd[2], xaxis, yaxis, starHolder[80],
                    minNumTrees, screenWidth, screenHeigth, numY[60], nSamples;
    RandLong        temporarySeed;
    BitsLong        *mask;
    PartCtr         *x;
    MrBFlt          xProb, yProb, xInc, yInc, xUpper, xLower, yUpper, yLower, *dT1=NULL, *dT2=NULL, *dT3=NULL, d1, d2, d3, 
                    meanY[60], xVal, yVal, minX, minY, maxX, maxY, sums[3];
    char            *s=NULL, prCh, treeName[2][100];
    FILE            *fp;
    time_t          curTime;
    PartCtr         **treeParts=NULL;
    Tree            *tree1=NULL, *tree2=NULL;
    SumtFileInfo    sumtFileInfo;
    
#   if defined (MPI_ENABLED)
    if (proc_id == 0)
        {
#   endif

    /* Make sure we read trees using DoSumtTree() code instead of with the user tree code */
    inComparetreeCommand = YES;

    /* set file pointer to NULL */
    fp = NULL;

    strcpy(treeName[0],"tree"); //in case if parameter is not specified in a .t file
    strcpy(treeName[1],"tree");

    /* Check that a data set has been read in. We check taxon names against
       those read in. */
    if (isTaxsetDef == NO)
        {
        MrBayesPrint ("%s   A matrix or set of taxon labels must be specified before comparetree can be used\n", spacer);
        goto errorExit;
        }

    /* open output files for summary information (two files); check if we want to overwrite previous results */
    if (OpenComptFiles () == ERROR)
        goto errorExit;

    MrBayesPrint ("%s   Examining files ...\n", spacer);

    /* Examine first file */
    if (ExamineSumtFile(comptreeParams.comptFileName1, &sumtFileInfo, treeName[0], &(brlensDef[0])) == ERROR)
        return ERROR;

    /* Capture info */
    longestLineLength      = sumtFileInfo.longestLineLength;
    numTreesInLastBlock[0] = sumtFileInfo.numTreesInLastBlock;
    lastTreeBlockBegin[0]  = sumtFileInfo.lastTreeBlockBegin;
    lastTreeBlockEnd[0]    = sumtFileInfo.lastTreeBlockEnd;

    /* Examine second file */
    if (ExamineSumtFile(comptreeParams.comptFileName2, &sumtFileInfo, treeName[1], &brlensDef[1]) == ERROR)
        return ERROR;

    /* Capture info */
    if (longestLineLength < sumtFileInfo.longestLineLength)
        longestLineLength = sumtFileInfo.longestLineLength;
    numTreesInLastBlock[1] = sumtFileInfo.numTreesInLastBlock;
    lastTreeBlockBegin[1]  = sumtFileInfo.lastTreeBlockBegin;
    lastTreeBlockEnd[1]    = sumtFileInfo.lastTreeBlockEnd;

    /* Check whether we should work with brlens */
    if (brlensDef[0] == YES && brlensDef[1] == YES)
        sumtParams.brlensDef = YES;
    else
        sumtParams.brlensDef = NO;
    
    /* Allocate space for command string */
    longestLineLength += 10;
    s = (char *)SafeMalloc((size_t)longestLineLength * sizeof(char));
    if (!s)
        {
        MrBayesPrint ("%s   Problem allocating string for reading tree file\n", spacer);
        goto errorExit;
        }

    /* Allocate space for packed trees */
    if (chainParams.relativeBurnin == YES)
        {
        numPackedTrees[0] = numTreesInLastBlock[0] - (int)(chainParams.burninFraction * numTreesInLastBlock[0]);
        numPackedTrees[1] = numTreesInLastBlock[1] - (int)(chainParams.burninFraction * numTreesInLastBlock[1]);
        }
    else
        {
        numPackedTrees[0] = numTreesInLastBlock[0] - chainParams.chainBurnIn;
        numPackedTrees[1] = numTreesInLastBlock[1] - chainParams.chainBurnIn;
        }
    if (memAllocs[ALLOC_PACKEDTREES] == YES)
        {
        MrBayesPrint ("%s   packedTreeList is already allocated\n", spacer);
        goto errorExit;
        }
    packedTreeList[0] = (PackedTree *) SafeCalloc(numPackedTrees[0]+numPackedTrees[1], sizeof(PackedTree));
    packedTreeList[1] = packedTreeList[0] + numPackedTrees[0];
    if (!packedTreeList[0])
        {
        MrBayesPrint ("%s   Problem allocating packed tree list\n", spacer);
        goto errorExit;
        }
    memAllocs[ALLOC_PACKEDTREES] = YES;

    /* Tell user we are ready to go */
    MrBayesPrint ("%s   Summarizing trees in files \"%s\" and \"%s\"\n", spacer,
        comptreeParams.comptFileName1,
        comptreeParams.comptFileName2);
    
    if (chainParams.relativeBurnin == YES)
        MrBayesPrint ("%s   Using relative burnin ('relburnin=yes'), discarding the first %.0f %% ('burninfrac=%1.2f') of sampled trees\n",
            spacer, chainParams.burninFraction*100.0, chainParams.burninFraction);
    else
        MrBayesPrint ("%s   Using absolute burnin ('relburnin=no'), discarding the first %d ('burnin=%d') sampled trees\n",
            spacer, chainParams.chainBurnIn, chainParams.chainBurnIn);

    MrBayesPrint ("%s   Writing statistics to file %s.<dists|pairs>\n", spacer, comptreeParams.comptOutfile);

    /* Set up cheap status bar. */
    MrBayesPrint ("\n%s   Tree reading status:\n\n", spacer);
    MrBayesPrint ("%s   0      10      20      30      40      50      60      70      80      90     100\n", spacer);
    MrBayesPrint ("%s   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v\n", spacer);
    MrBayesPrint ("%s   *", spacer);
    numAsterices = 0;
        
    /* Read file 1 for real */
    if ((fp = OpenTextFileR(comptreeParams.comptFileName1)) == NULL)
        goto errorExit;
        
    /* ...and fast forward to beginning of last tree block (skipping begin trees). */
    for (i=0; i<lastTreeBlockBegin[0] + 1; i++)
        {
        if (fgets (s, longestLineLength, fp) == NULL)
            {
            printf ("Error in function: %s at line: %d in file: %s", __func__, __LINE__, __FILE__);
            }
        }
        
    /* Calculate burnin */
    if (chainParams.relativeBurnin == YES)
        comptreeParams.burnin = (int)(chainParams.burninFraction * numTreesInLastBlock[0]);
    else
        comptreeParams.burnin = chainParams.chainBurnIn;

    /* Initialize sumtParams struct */
    numUniqueSplitsFound = numUniqueTreesFound = 0;
    sumtParams.runId = 0;
    strcpy(sumtParams.curFileName, comptreeParams.comptFileName1);
    sumtParams.tree = AllocatePolyTree (numTaxa);
    AllocatePolyTreePartitions (sumtParams.tree);
    sumtParams.numTreesEncountered = sumtParams.numTreesSampled = 0;
    sumtParams.numFileTrees = (int *) SafeCalloc (2*2+2*numTaxa, sizeof(int));
    sumtParams.numFileTreesSampled = sumtParams.numFileTrees + sumtParams.numRuns;
    sumtParams.order = sumtParams.numFileTrees + 2*sumtParams.numRuns;
    sumtParams.absentTaxa = sumtParams.numFileTrees + 2*sumtParams.numRuns + numTaxa;
    sumtParams.numTreesInLastBlock = numTreesInLastBlock[0];
    if (!sumtParams.numFileTrees)
        {
        MrBayesPrint ("%s   Problems allocating sumtParams.numFileTrees in DoSumt()\n", spacer);
        goto errorExit;
        }
    else
        memAllocs[ALLOC_SUMTPARAMS] = YES;

    /* ... and parse the file */
    expecting = Expecting(COMMAND);
    inTreesBlock = YES;
    ResetTranslateTable();
    for (i=0; i<lastTreeBlockEnd[0] - lastTreeBlockBegin[0] - 1; i++)
        {
        if (fgets (s, longestLineLength, fp) == NULL)
            {
            printf ("Error in function: %s at line: %d in file: %s", __func__, __LINE__, __FILE__);
            }
        /*MrBayesPrint ("%s", s);*/
        if (ParseCommand (s) == ERROR)
            goto errorExit;
        }
    inTreesBlock = NO;
    ResetTranslateTable();

    /* Check that at least one tree was read in. */
    if (sumtParams.numFileTreesSampled[0] <= 0)
        {
        MrBayesPrint ("%s   No trees read in\n", spacer);
        goto errorExit;
        }
        
    /* ... and close file */
    SafeFclose (&fp);

    /* Read file 2 for real */
    if ((fp = OpenTextFileR(comptreeParams.comptFileName2)) == NULL)
        goto errorExit;
        
    /* ...and fast forward to beginning of last tree block. */
    for (i=0; i<lastTreeBlockBegin[1] + 1; i++)
        {
        if (fgets (s, longestLineLength, fp) == NULL)
            {
            printf ("Error in function: %s at line: %d in file: %s", __func__, __LINE__, __FILE__);
            }
        }
        
    /* Renitialize sumtParams struct */
    sumtParams.runId = 1;
    strcpy (sumtParams.curFileName, comptreeParams.comptFileName2);

    /* Calculate burnin */
    if (chainParams.relativeBurnin == YES)
        comptreeParams.burnin = (int)(chainParams.burninFraction * numTreesInLastBlock[1]);
    else
        comptreeParams.burnin = chainParams.chainBurnIn;

    /* ... and parse the file */
    expecting = Expecting(COMMAND);
    inTreesBlock = YES;
    ResetTranslateTable();
    for (i=0; i<lastTreeBlockEnd[1] - lastTreeBlockBegin[1] - 1; i++)
        {
        if (fgets (s, longestLineLength, fp) == NULL)
            {
            printf ("Error in function: %s at line: %d in file: %s", __func__, __LINE__, __FILE__);
            }
        /*MrBayesPrint ("%s", s);*/
        if (ParseCommand (s) == ERROR)
            goto errorExit;
        }
    inTreesBlock = NO;
    ResetTranslateTable();

    /* Check that at least one tree was read in. */
    if (sumtParams.numFileTreesSampled[1] <= 0)
        {
        MrBayesPrint ("%s   No trees read in\n", spacer);
        goto errorExit;
        }
        
    /* ... and close file */
    SafeFclose (&fp);

    /* Now finish cheap status bar. */
    if (numAsterices < 80)
        for (i=0; i<80 - numAsterices; i++)
            MrBayesPrint ("*");
    MrBayesPrint ("\n\n");
    
    /* tell user how many trees were successfully read */
    MrBayesPrint ("%s   Read %d trees from last tree block of file \"%s\" (sampling %d of them)\n", spacer,
        sumtParams.numFileTrees[0],
        comptreeParams.comptFileName1,
        sumtParams.numFileTreesSampled[0]);
    MrBayesPrint ("%s   Read %d trees from last tree block of file \"%s\" (sampling %d of them)\n", spacer,
        sumtParams.numFileTrees[1],
        comptreeParams.comptFileName2,
        sumtParams.numFileTreesSampled[1]);
    
    /* Extract partition counter pointers */
    treeParts = (PartCtr **) SafeCalloc ((size_t)(numUniqueSplitsFound), sizeof(PartCtr *));
    i = 0;
    PartCtrUppass(partCtrRoot, treeParts, &i);

    /* Sort taxon partitions (clades, splits) ... */
    SortPartCtr (treeParts, 0, numUniqueSplitsFound-1);
        
    /* print to screen */
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   General explanation:                                                          \n", spacer);
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   In an unrooted tree, a taxon bipartition (split) is specified by removing a   \n", spacer);
    MrBayesPrint ("%s   branch, thereby dividing the species into those to the left and those to the  \n", spacer);
    MrBayesPrint ("%s   right of the branch. Here, taxa to one side of the removed branch are denoted \n", spacer);
    MrBayesPrint ("%s   '.' and those to the other side are denoted '*'. Specifically, the '.' symbol \n", spacer);
    MrBayesPrint ("%s   is used for the taxa on the same side as the outgroup.                        \n", spacer);
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   In a rooted or clock tree, the '*' symbol is simply used to denote the taxa   \n", spacer);
    MrBayesPrint ("%s   that are included in a particular group (clade), that is, all the descendants \n", spacer);
    MrBayesPrint ("%s   of a particular branch in the tree.  Taxa that are not included are denoted   \n", spacer);
    MrBayesPrint ("%s   using the '.' symbol.                                                         \n", spacer);
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   The output includes the ID of the encountered clades or splits (sorted from   \n", spacer);
    MrBayesPrint ("%s   highest to lowest probability), the bipartition or clade in '.*' format, \n", spacer);
    MrBayesPrint ("%s   number of times the bipartition or clade was observed in the first tree file  \n", spacer);
    MrBayesPrint ("%s   served in the first tree file, the number of times the bipartition was,       \n", spacer);
    MrBayesPrint ("%s   observed in the second tree file, the proportion of the time the bipartition  \n", spacer);
    MrBayesPrint ("%s   was found in the first tree file, and the proportion of the time the bi-      \n", spacer);
    MrBayesPrint ("%s   partition was found in the second tree file.                                  \n", spacer);
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   List of taxa in bipartitions:                                                 \n", spacer);
    MrBayesPrint ("                                                                                   \n");
    j = 1;
    for (k=0; k<numTaxa; k++)
        {
        if (sumtParams.absentTaxa[k] == NO && taxaInfo[k].isDeleted == NO)
            {
            MrBayesPrint ("%s   %4d -- %s\n", spacer, j++, taxaNames[k]);
            }
        }
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   List of taxon bipartitions found in tree files:                               \n\n", spacer);

    i = (int)(log10(sumtParams.numTreesSampled)) - 1;
    if (i<1)
        i = 1;
    j = sumtParams.numTaxa - 8;
    if (j < 1)
        j = 1;        
    MrBayesPrint ("%s     ID -- Partition%*c  No1%*c  No2%*c  Freq1   Freq2\n",
        spacer, j, ' ', i, ' ', i, ' ');

    mask = SafeCalloc (sumtParams.BitsLongsNeeded, sizeof(BitsLong));
    for (i=0; i<sumtParams.numTaxa; i++)
        SetBit (i, mask);
    for (i=0; i<numUniqueSplitsFound; i++)
        {
        x = treeParts[i];
        if (IsBitSet(localOutGroup, x->partition) == YES && sumtParams.isRooted == NO)
            FlipBits(x->partition, sumtParams.BitsLongsNeeded, mask);
        if ((MrBFlt)x->totCount/(MrBFlt)sumtParams.numTreesSampled >= comptreeParams.minPartFreq)
            {
            MrBayesPrint ("%s   %4d -- ", spacer, i+1);
            ShowParts (stdout, x->partition, sumtParams.numTaxa);

            j = (int)(log10(sumtParams.numTreesSampled)) + 1;
            if (j < 3)
                j = 3;
            MrBayesPrint ("   %*d   %*d   %1.3lf   %1.3lf\n", 
            j, x->count[0], j, x->count[1],
            (MrBFlt)x->count[0]/(MrBFlt)sumtParams.numFileTreesSampled[0], 
            (MrBFlt)x->count[1]/(MrBFlt)sumtParams.numFileTreesSampled[1]);
            
            MrBayesPrintf (fpParts, "%d\t%d\t%d\t%1.3lf\t%1.3lf\n", 
            i+1, x->count[0], x->count[1], 
            (MrBFlt)x->count[0]/(MrBFlt)sumtParams.numFileTreesSampled[0], 
            (MrBFlt)x->count[1]/(MrBFlt)sumtParams.numFileTreesSampled[1]);
            }
        }
    free (mask);

    /* make a nifty graph plotting frequencies of clades found in the two tree files */
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   Bivariate plot of clade probabilities:                                        \n", spacer);
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   This graph plots the probabilities of clades found in file 1 (the x-axis)     \n", spacer);
    MrBayesPrint ("%s   against the probabilities of the same clades found in file 2 (the y-axis).    \n", spacer);
    MrBayesPrint ("                                                                                   \n");
    xInc = (MrBFlt) (1.0 / 80.0);
    yInc = (MrBFlt) (1.0 / 40.0);
    yUpper = 1.0;
    yLower = yUpper - yInc;
    for (yaxis=39; yaxis>=0; yaxis--)
        {
        xLower = 0.0;
        xUpper = xLower + xInc;
        for (xaxis=0; xaxis<80; xaxis++)
            {
            starHolder[xaxis] = 0;
            for (i=0; i<numUniqueSplitsFound; i++)
                {
                x = treeParts[i];
                xProb = (MrBFlt)x->count[0]/(MrBFlt)sumtParams.numFileTreesSampled[0];
                yProb = (MrBFlt)x->count[1]/(MrBFlt)sumtParams.numFileTreesSampled[1];
                if ((xProb > xLower || (xProb == 0.0 && xaxis == 0)) && (xProb <= xUpper || (xProb == 1.0 && xaxis == 79))
                 && (yProb > yLower || (yProb == 0.0 && yaxis == 0)) && (yProb <= yUpper || (yProb == 1.0 && yaxis == 39)))
                    starHolder[xaxis] = 1;
                }
            xLower += xInc;
            xUpper = xLower + xInc;
            }
        
        MrBayesPrint ("%s   ", spacer);
        for (xaxis=0; xaxis<80; xaxis++)
            {
            prCh = ' ';
            if ((xaxis == 0 && yaxis == 0) || (xaxis == 79 && yaxis == 39))
                prCh = '+';
            else if ((xaxis == 0 && yaxis == 39) || (xaxis == 79 && yaxis == 0))
                prCh = '+';
            else if ((yaxis == 0 || yaxis == 39) && xaxis > 0 && xaxis < 79)
                prCh = '-';
            else if ((xaxis == 0 || xaxis == 79) && yaxis > 0 && yaxis < 39)
                prCh = '|';
            if (starHolder[xaxis] == 1)
                prCh = '*';
            MrBayesPrint ("%c", prCh);
            }
            if (yaxis == 39)
                MrBayesPrint (" 1.00\n");
            else if (yaxis == 0)
                MrBayesPrint (" 0.00\n");
            else
                MrBayesPrint ("\n");

        yUpper -= yInc;
        yLower = yUpper - yInc;
        }

    MrBayesPrint ("%s   ^                                                                              ^\n", spacer);
    MrBayesPrint ("%s  0.00                                                                          1.00\n", spacer);
        
    /* get tree-to-tree distances: first allocate some space */
    minNumTrees = sumtParams.numFileTreesSampled[0];
    if (sumtParams.numFileTreesSampled[1] < minNumTrees)
        minNumTrees = sumtParams.numFileTreesSampled[1];
    dT1 = (MrBFlt *) SafeMalloc (3 * (size_t)minNumTrees * sizeof(MrBFlt));
    tree1 = AllocateFixedTree (sumtParams.numTaxa, sumtParams.isRooted);
    tree2 = AllocateFixedTree (sumtParams.numTaxa, sumtParams.isRooted);
    if (!dT1 || !tree1 || !tree2)
        {
        MrBayesPrint ("%s   Problem allocating topological distances\n", spacer);
        goto errorExit;
        }
    dT2 = dT1 + minNumTrees;
    dT3 = dT2 + minNumTrees;
    
    for (i=0; i<minNumTrees; i++)
        {
        if (sumtParams.isRooted == NO)
            {
            RetrieveUTree (tree1, packedTreeList[0][i].order, packedTreeList[0][i].brlens);
            RetrieveUTree (tree2, packedTreeList[1][i].order, packedTreeList[1][i].brlens);
            }
        else
            {
            RetrieveRTree (tree1, packedTreeList[0][i].order, packedTreeList[0][i].brlens);
            RetrieveRTree (tree2, packedTreeList[1][i].order, packedTreeList[1][i].brlens);
            }
        /* Allocate and set partitions now that we have a tree */
        if (i == 0)
            {
            AllocateTreePartitions(tree1);
            AllocateTreePartitions(tree2);
            }
        else
            {
            ResetTreePartitions(tree1);
            ResetTreePartitions(tree2);
            }
        CalculateTreeToTreeDistance (tree1, tree2, &d1, &d2, &d3);
        dT1[i] = d1;
        dT2[i] = d2;
        dT3[i] = d3;
        }
        
    for (i=0; i<minNumTrees; i++)
        {
        /*MrBayesPrint ("%s   %4d -- %lf %lf %lf\n", spacer, i+1, dT1[i], dT2[i], dT3[i]);*/    
        if (sumtParams.brlensDef == YES)
            MrBayesPrintf (fpDists, "%d\t%lf\t%lf\t%lf\n", i+1, dT1[i], dT2[i], dT3[i]);    
        else
            MrBayesPrintf (fpDists, "%d\t%lf\n", i+1, dT1[i]);  
        }
        
    /* print x-y plot of log likelihood vs. generation */
    MrBayesPrint ("\n");
    MrBayesPrint ("%s   Rough plots of generation (x-axis) versus the measures of tree-   \n", spacer);
    MrBayesPrint ("%s   to-tree distances (y-axis).                                       \n", spacer);
    MrBayesPrint ("\n");
    MrBayesPrint ("%s   Distance(Robinson-Foulds):\n", spacer);
    MrBayesPrint ("\n");
    screenWidth = 60; /* don't change this without changing numY and meanY, declared above */
    screenHeigth = 15;
    minX = minY = 1000000000.0;
    maxX = maxY = -1000000000.0;
    for (i=0; i<minNumTrees; i++)
        {
        xVal = (MrBFlt) (i + comptreeParams.burnin);
        yVal = dT1[i];
        if (xVal < minX)
            minX = xVal;
        if (yVal < minY)
            minY = yVal;
        if (xVal > maxX)
            maxX = xVal;
        if (yVal > maxY)
            maxY = yVal;
        }
    for (i=0; i<screenWidth; i++)
        {
        numY[i] = 0;
        meanY[i] = 0.0;
        }
    for (i=0; i<minNumTrees; i++)
        {
        xVal = (MrBFlt) (i + comptreeParams.burnin);
        yVal = dT1[i];
        k = (int)(((xVal - minX) / (maxX - minX)) * screenWidth);
        if (k >= screenWidth)
            k = screenWidth - 1;
        meanY[k] += yVal;
        numY[k]++;
        }
    MrBayesPrint ("\n%s   +", spacer);
    for (i=0; i<screenWidth; i++)
        MrBayesPrint ("-");
    MrBayesPrint ("+ %1.2lf\n", maxY);
    for (j=screenHeigth-1; j>=0; j--)
        {
        MrBayesPrint ("%s   |", spacer);
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
    MrBayesPrint ("%s   +", spacer);
    for (i=0; i<screenWidth; i++)
        {
        if (numY[i] > 0 && meanY[i] / numY[i] <= minY)
            MrBayesPrint ("*");
        else if (i % (screenWidth/10) == 0 && i != 0)
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
    for (i=0; i<screenWidth; i++)
        MrBayesPrint (" ");
    MrBayesPrint ("%1.0lf\n\n", maxX);

    if (sumtParams.brlensDef == YES)
        {
        for (n=0; n<2; n++)
            {
            MrBayesPrint ("\n");
            if (n == 0)
                MrBayesPrint ("%s   Distance(Robinson-Foulds with branch lengths):\n", spacer);
            else
                MrBayesPrint ("%s   Distance(Robinson-Foulds with scaled branch lengths):\n", spacer);
            MrBayesPrint ("\n");
            screenWidth = 60; /* don't change this without changing numY and meanY, declared above */
            screenHeigth = 15;
            minX = minY = 1000000000.0;
            maxX = maxY = -1000000000.0;
            for (i=0; i<minNumTrees; i++)
                {
                xVal = (MrBFlt) (i + comptreeParams.burnin);
                if (n == 0)
                    yVal = dT2[i];
                else
                    yVal = dT3[i];
                if (xVal < minX)
                    minX = xVal;
                if (yVal < minY)
                    minY = yVal;
                if (xVal > maxX)
                    maxX = xVal;
                if (yVal > maxY)
                    maxY = yVal;
                }
            for (i=0; i<screenWidth; i++)
                {
                numY[i] = 0;
                meanY[i] = 0.0;
                }
            for (i=0; i<minNumTrees; i++)
                {
                xVal = (MrBFlt) (i + comptreeParams.burnin);
                if (n == 0)
                    yVal = dT2[i];
                else
                    yVal = dT3[i];
                k = (int)(((xVal - minX) / (maxX - minX)) * screenWidth);
                if (k >= screenWidth)
                    k = screenWidth - 1;
                meanY[k] += yVal;
                numY[k]++;
                }
            MrBayesPrint ("\n%s   +", spacer);
            for (i=0; i<screenWidth; i++)
                MrBayesPrint ("-");
            MrBayesPrint ("+ %1.2lf\n", maxY);
            for (j=screenHeigth-1; j>=0; j--)
                {
                MrBayesPrint ("%s   |", spacer);
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
            MrBayesPrint ("%s   +", spacer);
            for (i=0; i<screenWidth; i++)
                {
                if (numY[i] > 0 && meanY[i] / numY[i] <= minY)
                    MrBayesPrint ("*");
                else if (i % (screenWidth/10) == 0 && i != 0)
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
            for (i=0; i<screenWidth; i++)
                MrBayesPrint (" ");
            MrBayesPrint ("%1.0lf\n\n", maxX);
            }
        }

    /* calculate average tree-to-tree distances */
    curTime = time(NULL);
    temporarySeed  = (RandLong)curTime;
    if (temporarySeed < 0)
        temporarySeed = -temporarySeed;
    sums[0] = sums[1] = sums[2] = 0.0;
    nSamples = 1000;
    for (n=0; n<nSamples; n++)
        {
        i = (int) RandomNumber(&temporarySeed) * minNumTrees;
        j = (int) RandomNumber(&temporarySeed) * minNumTrees;
        if (sumtParams.isRooted == NO)
            {
            RetrieveUTree (tree1, packedTreeList[0][i].order, packedTreeList[0][i].brlens);
            RetrieveUTree (tree2, packedTreeList[1][j].order, packedTreeList[1][j].brlens);
            }
        else /* if (sumtParams.isRooted == YES) */
            {
            RetrieveRTree (tree1, packedTreeList[0][i].order, packedTreeList[0][i].brlens);
            RetrieveRTree (tree2, packedTreeList[1][j].order, packedTreeList[1][j].brlens);
            }
        CalculateTreeToTreeDistance (tree1, tree2, &d1, &d2, &d3);
        sums[0] += d1;
        sums[1] += d2;
        sums[2] += d3;
        }
    MrBayesPrint ("%s   Mean tree-to-tree distances, based on %d trees randomly sampled from both files:\n\n", spacer, nSamples);
    MrBayesPrint ("%s                                 Mean(Robinson-Foulds) = %1.3lf\n", spacer, sums[0]/nSamples);
    if (sumtParams.brlensDef == YES)
        {
        MrBayesPrint ("%s             Mean(Robinson-Foulds with branch lengths) = %1.3lf\n", spacer, sums[1]/nSamples);
        MrBayesPrint ("%s      Mean(Robinson-Foulds with scaled branch lengths) = %1.3lf\n", spacer, sums[2]/nSamples);
        }
    MrBayesPrint ("\n");

    /* free memory and file pointers */
    free(s);    
    free (dT1);
    FreeTree (tree1);
    FreeTree (tree2);
    if (memAllocs[ALLOC_PACKEDTREES] == YES)
        {
        for (i=0; i<numPackedTrees[0]+numPackedTrees[1]; i++)
            {
            free(packedTreeList[0][i].order);
            free(packedTreeList[0][i].brlens);
            }
        free(packedTreeList[0]);
        memAllocs[ALLOC_PACKEDTREES] = NO;
        }

    /* free sumtParams */
    FreeSumtParams();

    /* close files */
    SafeFclose (&fp);
    SafeFclose (&fpParts);
    SafeFclose (&fpDists);
    
    /* free pointer array to partitions, part and tree counters */
    free (treeParts);
    FreePartCtr (partCtrRoot);
    FreeTreeCtr (treeCtrRoot);
    partCtrRoot = NULL;
    treeCtrRoot = NULL;

    /* reset taxon set */
    ResetTaxonSet();

#   if defined (MPI_ENABLED)
        }
#   endif

    expecting = Expecting(COMMAND);
    inComparetreeCommand = NO;
    return (NO_ERROR);
    
    /* error exit */            
    errorExit:
        if (s) free(s);

        /* free sumtParams */
        FreeSumtParams();

        free (dT1);
        FreeTree (tree1);
        FreeTree (tree2);
        if (memAllocs[ALLOC_PACKEDTREES] == YES)
            {
            for (i=0; i<numPackedTrees[0]+numPackedTrees[1]; i++)
                {
                free(packedTreeList[0][i].order);
                free(packedTreeList[0][i].brlens);
                }
            free(packedTreeList[0]);
            memAllocs[ALLOC_PACKEDTREES] = NO;
            }

        /* reset taxon set */
        ResetTaxonSet();

        /* free pointer array to partitions, part and tree counters */
        free (treeParts);
        FreePartCtr (partCtrRoot);
        FreeTreeCtr (treeCtrRoot);
        partCtrRoot = NULL;
        treeCtrRoot = NULL;

        SafeFclose (&fp);
        SafeFclose (&fpParts);
        SafeFclose (&fpDists);
        strcpy (spacer, "");

        expecting = Expecting(COMMAND);
        inComparetreeCommand = NO;

        return (ERROR);
}


int DoCompareTreeParm (char *parmName, char *tkn)
{
    int         tempI;
    MrBFlt      tempD;
    char        tempStr[100];

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
        /* set Filename (comptreeParams.comptFileName1) *************************************************/
        else if (!strcmp(parmName, "Filename1"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (comptreeParams.comptFileName1, tkn);
                MrBayesPrint ("%s   Setting comparetree filename 1 to %s\n", spacer, comptreeParams.comptFileName1);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Filename (comptreeParams.comptFileName2) *************************************************/
        else if (!strcmp(parmName, "Filename2"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (comptreeParams.comptFileName2, tkn);
                MrBayesPrint ("%s   Setting comparetree filename 2 to %s\n", spacer, comptreeParams.comptFileName2);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Filename (comptreeParams.comptOutfile) ***************************************************/
        else if (!strcmp(parmName, "Outputname"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (comptreeParams.comptOutfile, tkn);
                MrBayesPrint ("%s   Setting comparetree output file to %s\n", spacer, comptreeParams.comptOutfile);
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
                //free (tempStr);
                return (ERROR);
                }
            }
        /* set Burnin (chainParams.chainBurnIn) *******************************************************/
        else if (!strcmp(parmName, "Burnin"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                chainParams.chainBurnIn = tempI;
                MrBayesPrint ("%s   Setting burnin to %d\n", spacer, chainParams.chainBurnIn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
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
        /* set Minpartfreq (comptreeParams.minPartFreq) *******************************************************/
        else if (!strcmp(parmName, "Minpartfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                comptreeParams.minPartFreq = tempD;
                MrBayesPrint ("%s   Including partitions with probability greater than or equal to %lf in summary statistics\n", spacer, comptreeParams.minPartFreq);
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


int DoCompRefTree (void)
{
    /* Compare a tree file with the reference tree files to generate the SDSFs.
       Use parameters in CompareTree and MCMCP (lazy option) */
    
    char         outName[130], inName[130], inRefName[130], treeName[100], *lineBuf=NULL, *s;
    FILE         *fpTre=NULL, *fpOut=NULL;
    int          i, n, longestL=0, burnin, gen, nRefRun, nTre[100]={0};
    SumtFileInfo tFileInfo;
    Tree         *t;

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return NO_ERROR;
#   endif

    /* Check that a data set has been read in. We check taxon names against those read in. */
    if (isTaxsetDef == NO)
        {
        MrBayesPrint ("%s   A matrix or set of taxon labels must be specified before compareref can be used\n", spacer);
        return (ERROR);
        }

    /* this is a hack to account for the additional comparing tree sample 
       chainParams.numRuns is used in AddTreeToPartitionCounters to alloc mem correctly */
    nRefRun = chainParams.numRuns;
    chainParams.numRuns += 1;

    /* initialize */
    if ((t = AllocateTree (numLocalTaxa)) == NULL)
        {
        MrBayesPrint ("%s   Problem allocating diagn tree\n", spacer);
        goto errorExit;
        }
    if (SetUpPartitionCounters () == ERROR)
        goto errorExit;
    if ((chainParams.stat = (STATS *) SafeCalloc (numTopologies, sizeof (STATS))) == NULL)
        goto errorExit;
    else
        memAllocs[ALLOC_STATS] = YES;

    /* deal with the reference tree files */
    // for (k=0; k<numTopologies; k++)
    for (n=0; n<nRefRun; n++)
        {
        if (nRefRun == 1)
            sprintf (inRefName, "%s.t", comptreeParams.comptFileName2);
        else
            sprintf (inRefName, "%s.run%d.t", comptreeParams.comptFileName2, n+1);
        
        /* Examine each ref tree file, save info to tFileInfo */
        if (ExamineSumtFile(inRefName, &tFileInfo, treeName, &sumtParams.brlensDef) == ERROR)
            goto errorExit;
        if (longestL < tFileInfo.longestLineLength)
            {
            longestL = tFileInfo.longestLineLength;
            lineBuf = (char *) SafeRealloc (lineBuf, (size_t)longestL * sizeof(char));
            if (!lineBuf)
                {
                MrBayesPrint ("%s   Problem allocating string for reading tree file\n", spacer);
                goto errorExit;
                }
            }
        
        /* calculate burnin */
        if (chainParams.relativeBurnin == YES)
            burnin = (int)(chainParams.burninFraction * tFileInfo.numTreesInLastBlock);
        else
            burnin = chainParams.chainBurnIn;
        if (burnin >= tFileInfo.numTreesInLastBlock)
            {
            MrBayesPrint ("%s   Burnin should be smaller than the total number of trees\n", spacer);
            goto errorExit;
            }
       
        /* open the ref tree file */
        if ((fpTre = OpenTextFileR (inRefName)) == NULL)
            goto errorExit;
        /* ...and fast forward to beginning in last tree block */
        for (i=0; i <= tFileInfo.lastTreeBlockBegin; i++)
            {
            if (fgets(lineBuf, longestL-2, fpTre) == NULL)
                goto errorExit;
            }
        
        /* process each ref tree */
        for (i=1; i <= tFileInfo.numTreesInLastBlock; i++)
            {
            do {
                if (fgets (lineBuf, longestL-2, fpTre) == NULL)
                    goto errorExit;
                s = strtok (lineBuf, " ");
                }
            while (strcmp (s, "tree") != 0);
            
            /* add reference trees to partition counters, discarding burnin */
            if (i > burnin)
                {
                s = strtok (NULL, ";");
                while (*s != '(')
                    s++;
                StripComments(s);

                if (ResetTopology (t, s) == ERROR)
                    goto errorExit;
                if (AddTreeToPartitionCounters (t, 0, n) == ERROR)
                    goto errorExit;
                nTre[n]++;
                }
            }

        /* close the tree file */
        SafeFclose (&fpTre);
        }
    /* end reference tree files */
    
    /* open output file */
    strcpy (outName, comptreeParams.comptOutfile);
    strcat (outName, ".sdsf");
    if ((fpOut = OpenNewMBPrintFile (outName)) == NULL)
        goto errorExit;
    /* print stamp and header */
    if ((int)strlen(stamp) > 1)
        MrBayesPrintf (fpOut, "[ID: %s]\n", stamp);
    if (chainParams.diagnStat == AVGSTDDEV)
        MrBayesPrintf (fpOut, "Gen\tASDSF\n");
    else  // MAXSTDDEV
        MrBayesPrintf (fpOut, "Gen\tMSDSF\n");

    /* Examine the tree file to be compared, save info to tFileInfo */
    strcpy(inName, comptreeParams.comptFileName1);
    if (ExamineSumtFile(inName, &tFileInfo, treeName, &sumtParams.brlensDef) == ERROR)
        goto errorExit;
    if (longestL < tFileInfo.longestLineLength)
        {
        longestL = tFileInfo.longestLineLength;
        lineBuf = (char *) SafeRealloc (lineBuf, (size_t)longestL * sizeof(char));
        if (!lineBuf)
            {
            MrBayesPrint ("%s   Problem allocating string for reading tree file\n", spacer);
            goto errorExit;
            }
        }

    /* open the tree file to be compared */
    if ((fpTre = OpenTextFileR (inName)) == NULL)
        goto errorExit;
    /* ...and fast forward to beginning in last tree block */
    for (i=0; i <= tFileInfo.lastTreeBlockBegin; i++)
        {
        if (fgets(lineBuf, longestL-2, fpTre) == NULL)
            goto errorExit;
        }

    /* process each tree to be compared and print SDSF to file */
    for (i=1; i <= tFileInfo.numTreesInLastBlock; i++)
        {
        do {
            if (fgets (lineBuf, longestL-2, fpTre) == NULL)
                goto errorExit;
            s = strtok (lineBuf, " ");
            }
        while (strcmp (s, "tree") != 0);
        
        s = strtok (NULL, ";");
        gen = atoi(s+4);  // 4 is offset to get rid of "rep." in tree name
        while (*s != '(')
            s++;
        StripComments(s);

        /* add the tree to partition counters */
        if (ResetTopology (t, s) == ERROR)
            goto errorExit;
        if (AddTreeToPartitionCounters (t, 0, nRefRun) == ERROR)
            goto errorExit;
        nTre[nRefRun]++;
            
        /* calculate and write stdev of split freq */
        CalcTopoConvDiagn2 (nTre);
        if (chainParams.stat[0].numPartitions == 0)
            {
            MrBayesPrintf (fpOut, "%d\tNA\n", gen);
            }
        else if (chainParams.diagnStat == AVGSTDDEV)
            {
            MrBayesPrintf (fpOut, "%d\t%lf\n", gen, chainParams.stat[0].avgStdDev);
            }
        else  // MAXSTDDEV
            {
            MrBayesPrintf (fpOut, "%d\t%lf\n", gen, chainParams.stat[0].max);
            }
        }
    
    /* change back to the actual numRuns, end of hack */
    chainParams.numRuns -= 1;
    
    /* close tree file */
    SafeFclose (&fpTre);
    /* close output file */
    SafeFclose (&fpOut);
    /* free memory */
    free(lineBuf);
    FreeTree (t);
    FreeChainMemory();
    
    return (NO_ERROR);
    
    /* error exit */
errorExit:
    MrBayesPrint ("%s   Error in DoCompRefTree\n", spacer);
    
    chainParams.numRuns -= 1;
    SafeFclose (&fpTre);
    SafeFclose (&fpOut);
    
    FreeTree (t);
    FreeChainMemory();
    free(lineBuf);
    
    return (ERROR);
}


#if defined (PRINT_RATEMUL_CPP)
int DELETE_ME_count_taxa(PolyNode *p)
{
    int sum=0;

    if (p->left==NULL) {
        if (p->depth > 0.1)
            return 1;
        else
            return 0;
    }

    p=p->left;
    while (p != NULL) {
        sum+=DELETE_ME_count_taxa(p);
        p=p->sib;
    }
    return sum;
}


void DELETE_ME_dump_depth(PolyNode *p)
{
    /*print depth of two taxa clade*/
    /*
    if (p->left != NULL && p->left->left == NULL && p->left->sib != NULL && p->left->sib->left == NULL){
        fprintf (rateMultfp,"%f\n",p->depth);
        }
    if (p->left != NULL && p->left->left == NULL && p->left->sib != NULL && p->left->sib->left == NULL){
        if (p->left->depth > 0.1 && p->left->sib->depth > 0.1)
            fprintf (rateMultfp,"%f\n",p->depth);
        }
    */
    /* print depth of three taxa clade */
    if (((p->left != NULL && p->left->left == NULL) && p->left->sib != NULL && p->left->sib->left != NULL &&  p->left->sib->left->left == NULL && p->left->sib->left->sib->left == NULL) ||
       (p->left != NULL && p->left->left != NULL && p->left->left->left == NULL && p->left->left->sib->left == NULL && (p->left->sib->left == NULL))) {
        if (DELETE_ME_count_taxa(p)==2)
            fprintf (rateMultfp,"%f\n",p->depth);
        }
 
    p=p->left;
    while (p != NULL) {
        DELETE_ME_dump_depth(p);
        p=p->sib;
    }
}
#endif


int DoSumt (void)
{
    int             i, j=0, k, n, len, longestName, treeNo, numTreePartsToPrint,
                    maxWidthID, maxWidthNumberPartitions, maxNumTaxa, tableWidth=0, unreliable, oneUnreliable,
                    longestHeader;
    MrBFlt          f, var_s, sum_s, stddev_s=0.0, sumsq_s, sumStdDev=0.0, maxStdDev=0.0, sumPSRF=0.0,
                    maxPSRF=0.0, avgStdDev=0.0, avgPSRF=0.0, min_s=0.0, max_s=0.0, numPSRFSamples=0, min;
    PartCtr         *x;
    char            *s=NULL, tempName[120], fileName[120], treeName[100], divString[100];
    char            *tempStr=NULL; /*not static because error ext is handeled*/
    int             tempStrLength;
    FILE            *fp=NULL;
    PartCtr         **treeParts=NULL,*tmp;
    SumtFileInfo    sumtFileInfo;
    Stat            theStats;
    BitsLong        *mask;

#define SCREEN_WIDTH 80
    
#   if defined (MPI_ENABLED)
    if (proc_id == 0)
        {
#   endif

    /* Ensure that we read trees with sumt code and not user tree code */
    inSumtCommand = YES;

    /* set file pointers to NULL */
    fp = fpParts = fpTstat = fpVstat = fpCon = fpTrees = NULL;

    strcpy(treeName,"tree"); //in case if parameter is not specified in a .t file

    /* Check if there is anything to do */
    if (sumtParams.table == NO && sumtParams.summary == NO && sumtParams.showConsensus == NO)
        {
        MrBayesPrint ("%s   Nothing to do, all output parameters (Table, Summary, Consensus) set to 'NO'\n", spacer);
        goto errorExit;
        }

    SafeStrcpy(&tempStr, " ");

    /* Initialize sumtParams struct */
    sumtParams.numTaxa = 0;
    sumtParams.taxaNames = NULL;
    sumtParams.BitsLongsNeeded = 0;
    sumtParams.tree = AllocatePolyTree (numTaxa);
    sumtParams.nESets = 0;
    sumtParams.nBSets = 0;
    sumtParams.eSetName = NULL;
    sumtParams.bSetName = NULL;
    sumtParams.popSizeSet = NO;
    sumtParams.popSizeSetName = NULL;
    AllocatePolyTreePartitions (sumtParams.tree);
    sumtParams.numFileTrees = (int *) SafeCalloc (2*sumtParams.numRuns+2*numTaxa, sizeof(int));
    sumtParams.numFileTreesSampled = sumtParams.numFileTrees + sumtParams.numRuns;
    sumtParams.order = sumtParams.numFileTrees + 2*sumtParams.numRuns;
    sumtParams.absentTaxa = sumtParams.numFileTrees + 2*sumtParams.numRuns + numTaxa;
    if (!sumtParams.numFileTrees)
        {
        MrBayesPrint ("%s   Problems allocating sumtParams.numFileTrees in DoSumt()\n", spacer);
        goto errorExit;
        }
    else
        memAllocs[ALLOC_SUMTPARAMS] = YES;

    /* Make sure outgroup is set correctly */
    
    for (treeNo = 0; treeNo < sumtParams.numTrees; treeNo++)
        {
        /* initialize across-file tree and partition counters */
        sumtParams.numTreesSampled = sumtParams.numTreesEncountered = 0;
        numUniqueSplitsFound = numUniqueTreesFound = 0;

        /* initialize oneUnreliable && unreliable */
        oneUnreliable = unreliable = NO;

        /* initialize summary statistics */
        sumStdDev = 0.0;
        sumPSRF = 0.0;
        numPSRFSamples = 0;
        maxStdDev = 0.0;
        maxPSRF = 0.0;

        /* tell user we are ready to go */
        if (sumtParams.numTrees > 1)
            sprintf (fileName,"%s.tree%d", sumtParams.sumtFileName, treeNo+1);
        else
            strcpy (fileName, sumtParams.sumtFileName);

        if (sumtParams.numRuns == 1)
            MrBayesPrint ("%s   Summarizing trees in file \"%s.t\"\n", spacer, fileName);
        else if (sumtParams.numRuns == 2)
            MrBayesPrint ("%s   Summarizing trees in files \"%s.run1.t\" and \"%s.run2.t\"\n", spacer, fileName, fileName);
        else if (sumtParams.numRuns > 2)
            MrBayesPrint ("%s   Summarizing trees in files \"%s.run1.t\", \"%s.run2.t\",...,\"%s.run%d.t\"\n", spacer, fileName, fileName,fileName,sumtParams.numRuns);

        if (chainParams.relativeBurnin == YES)
            MrBayesPrint ("%s   Using relative burnin ('relburnin=yes'), discarding the first %.0f %% of sampled trees\n",
                spacer, chainParams.burninFraction*100.0, chainParams.burninFraction);
        else
            MrBayesPrint ("%s   Using absolute burnin ('relburnin=no'), discarding the first %d sampled trees\n",
                spacer, chainParams.chainBurnIn, chainParams.chainBurnIn);
        
        MrBayesPrint ("%s   Writing statistics to files %s.<parts|tstat|vstat|trprobs|con>\n", spacer, sumtParams.sumtOutfile);

        for (sumtParams.runId=0; sumtParams.runId < sumtParams.numRuns; sumtParams.runId++)
            {
            /* initialize tree counters */
            sumtParams.numFileTrees[sumtParams.runId] = sumtParams.numFileTreesSampled[sumtParams.runId] = 0;

            /* open binary file */
            if (sumtParams.numRuns == 1)
                sprintf (tempName, "%s.t", fileName);
            else
                sprintf (tempName, "%s.run%d.t", fileName, sumtParams.runId+1);
            strcpy(sumtParams.curFileName, tempName);

            /* tell user we are examining files if for the first run */
            if (sumtParams.runId == 0)
                {
                if (sumtParams.numRuns > 1 && sumtParams.numTrees > 1)
                    MrBayesPrint ("%s   Examining first file for tree %d ...\n", spacer, treeNo);
                else if (sumtParams.numRuns > 1 && sumtParams.numTrees == 1)
                    MrBayesPrint ("%s   Examining first file ...\n", spacer);
                else if (sumtParams.numRuns == 1 && sumtParams.numTrees > 1)
                    MrBayesPrint ("%s   Examining file for tree %d ...\n", spacer, treeNo);
                else
                    MrBayesPrint ("%s   Examining file ...\n", spacer);
                }

            /* examine file */
            if (ExamineSumtFile(tempName, &sumtFileInfo, treeName, &sumtParams.brlensDef) == ERROR)
                goto errorExit;

            /* capture values */
            if (sumtParams.runId == 0)

            /* catch lack of sampled trees */
            if (chainParams.relativeBurnin == NO && chainParams.chainBurnIn > sumtFileInfo.numTreesInLastBlock)
                {
                MrBayesPrint ("%s   No trees are sampled as the burnin exceeds the number of trees in last block\n", spacer);
                MrBayesPrint ("%s   Try setting burnin to a number less than %d\n", spacer, sumtFileInfo.numTreesInLastBlock);
                goto errorExit;
                }
                
            /* tell the user that everything is fine */
            if (sumtParams.runId == 0)
                {
                if (sumtFileInfo.numTreeBlocks == 1)
                    MrBayesPrint ("%s   Found one tree block in file \"%s\" with %d trees in last block\n",
                        spacer, tempName, sumtFileInfo.numTreesInLastBlock);
                else
                    {
                    MrBayesPrint ("%s   Found %d tree blocks in file \"%s\" with %d trees in last block\n",
                        spacer, sumtFileInfo.numTreeBlocks, tempName, sumtFileInfo.numTreesInLastBlock);
                    MrBayesPrint ("%s   Only the %d trees in last tree block will be summarized\n", spacer, sumtFileInfo.numTreesInLastBlock);
                    }
                sumtParams.numTreesInLastBlock = sumtFileInfo.numTreesInLastBlock;
                if (sumtParams.numRuns > 1)
                    MrBayesPrint ("%s   Expecting the same number of trees in the last tree block of all files\n", spacer);
                if (chainParams.relativeBurnin == NO)
                    sumtParams.burnin = chainParams.chainBurnIn;
                else
                    sumtParams.burnin = (int) (sumtFileInfo.numTreesInLastBlock * chainParams.burninFraction);
                }
            else
                {
                if (sumtFileInfo.numTreesInLastBlock != sumtParams.numFileTrees[0])
                    {
                    MrBayesPrint ("%s   Found %d trees in first file but %d trees in file \"%s\"\n", spacer,
                        sumtParams.numFileTrees[0],
                        sumtFileInfo.numTreesInLastBlock,
                        tempName);
                    goto errorExit;
                    }
                }
        
            /* Now we read the file for real. First, allocate a string for reading the file... */
            if (sumtParams.runId == 0 && treeNo == 0)
                {
                s = (char *) SafeMalloc ((size_t)(sumtFileInfo.longestLineLength) * sizeof(char));
                if (!s)
                    {
                    MrBayesPrint ("%s   Problem allocating string for reading sumt file\n", spacer);
                    goto errorExit;
                    }
                }
            else
                {
                free (s);
                s = (char *) SafeMalloc ((size_t)(sumtFileInfo.longestLineLength) * sizeof (char));
                if (!s)
                    {
                    MrBayesPrint ("%s   Problem reallocating string for reading sumt file\n", spacer);
                    goto errorExit;
                    }
                }
        
            /* ... open the file ... */
            if ((fp = OpenTextFileR(tempName)) == NULL)
                goto errorExit;
    
            /* ...and fast forward to beginning of last tree block. */
            for (i=0; i<sumtFileInfo.lastTreeBlockBegin+1; i++)
                {
                if (fgets (s, sumtFileInfo.longestLineLength-2, fp) == NULL)
                    {
                    printf ("Error in function: %s at line: %d in file: %s", __func__, __LINE__, __FILE__);
                    }
                }

#   if defined (PRINT_RATEMUL_CPP)
            sprintf (tempName, "%s.ratemult", chainParams.chainFileName);
            if ((rateMultfp=OpenNewMBPrintFile (tempName)) == NULL)
                {
                printf ("Error oppening file: %s to write", tempName);
                goto errorExit;
                }
            fprintf (rateMultfp,"rateMult_CPP\n");
#   endif

            /* Set up cheap status bar. */
            if (sumtParams.runId == 0)
                {
                if (sumtParams.numTrees > 1)
                    MrBayesPrint ("\n%s   Tree reading status for tree %d:\n\n", spacer, treeNo+1);
                else
                    MrBayesPrint ("\n%s   Tree reading status:\n\n", spacer);
                MrBayesPrint ("%s   0      10      20      30      40      50      60      70      80      90     100\n", spacer);
                MrBayesPrint ("%s   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v\n", spacer);
                MrBayesPrint ("%s   *", spacer);
                numAsterices = 0;
                }
        
            /* ...and parse file, tree-by-tree. We are only parsing lines between the "begin trees" and "end" statements.
            We don't actually get those lines, however, but rather the lines between those statements. */
            expecting = Expecting(COMMAND);
            /* We skip the begin trees statement so we need to set up some variables here */
            inTreesBlock = YES;
            ResetTranslateTable();
            for (i=0; i<sumtFileInfo.lastTreeBlockEnd - sumtFileInfo.lastTreeBlockBegin - 1; i++)
                {
                if (fgets (s, sumtFileInfo.longestLineLength-2, fp) == NULL)
                    {
                    printf ("Error in function: %s at line: %d in file: %s", __func__, __LINE__, __FILE__);
                    }
                /*MrBayesPrint ("%s", s);*/
                if (ParseCommand (s) == ERROR)
                    goto errorExit;
                }
            inTreesBlock = NO;
            ResetTranslateTable();
    
            /* Finish cheap status bar. */
            if (sumtParams.runId == sumtParams.numRuns - 1)
                {
                if (numAsterices < 80)
                    {
                    for (i=0; i<80 - numAsterices; i++)
                        MrBayesPrint ("*");
                    }
                MrBayesPrint ("\n\n");
                }
    
            /* print out information on absent and pruned taxa */
            if (sumtParams.runId == sumtParams.numRuns - 1 && treeNo == 0)
                PrintSumtTaxaInfo ();

            /* tell user how many trees were successfully read */
            if (sumtParams.numRuns == 1)
                MrBayesPrint ("%s   Read %d trees from last tree block (sampling %d of them)\n", spacer,
                    sumtParams.numTreesEncountered, sumtParams.numTreesSampled);
            else if (sumtParams.numRuns > 1)
                {
                if (sumtParams.runId != 0 && sumtParams.numFileTreesSampled[sumtParams.runId]!=sumtParams.numFileTreesSampled[0])
                    {
                    if (sumtParams.runId == sumtParams.numRuns - 1)
                        MrBayesPrint ("\n\n");
                    MrBayesPrint ("%s   Found %d post-burnin trees in the first file but %d post-burnin trees in file %d\n",
                            spacer,
                            sumtParams.numFileTreesSampled[0],
                            sumtParams.numFileTreesSampled[sumtParams.runId],
                            sumtParams.runId+1);
                    goto errorExit;
                    }
                if (sumtParams.runId == sumtParams.numRuns - 1)
                    {
                    MrBayesPrint ("%s   Read a total of %d trees in %d files (sampling %d of them)\n", spacer,
                        sumtParams.numTreesEncountered,
                        sumtParams.numRuns, sumtParams.numTreesSampled);
                    MrBayesPrint ("%s      (Each file contained %d trees of which %d were sampled)\n", spacer,
                        sumtParams.numFileTrees[0],
                        sumtParams.numFileTreesSampled[0]);
                    }
                }

            /* Check that at least one tree was read in. */
            if (sumtParams.numTreesSampled <= 0)
                {
                MrBayesPrint ("%s   No trees read in\n", spacer);
                goto errorExit;
                }

            SafeFclose (&fp);
            }   /* next run for this tree */
                
        /* Extract partition counter pointers */
        treeParts = (PartCtr **) SafeCalloc ((size_t)numUniqueSplitsFound, sizeof(PartCtr *));
        i = 0;
        PartCtrUppass(partCtrRoot, treeParts, &i);

        min = (sumtParams.minPartFreq * (sumtParams.numTreesSampled/sumtParams.numRuns));
        numTreePartsToPrint=numUniqueSplitsFound;
        for (i=0; i<numTreePartsToPrint;)
            {
            for (j=0; j<sumtParams.numRuns;j++)
                {
                if (treeParts[i]->count[j]>=min)
                    break;
                }
            if (j==sumtParams.numRuns)
                {
                numTreePartsToPrint--;
                tmp=treeParts[numTreePartsToPrint];
                treeParts[numTreePartsToPrint]=treeParts[i];
                treeParts[i]=tmp;
                }
            else
                {
                i++;
                }
            }

        /* Sort taxon partitions (clades, splits) ... */
        SortPartCtr (treeParts, 0, numTreePartsToPrint-1);

        /* Sort root and tips among those splits always present */
        SortTerminalPartCtr (treeParts, numUniqueSplitsFound);

        /* open output files for summary information (three files) */
        if (OpenSumtFiles (treeNo) == ERROR)
            goto errorExit;

        /* Print partitions to screen. */
        if (treeNo == 0)
            {
            longestName = 0;
            for (k=0; k<numTaxa; k++)
                {
                if (taxaInfo[k].isDeleted == NO && sumtParams.absentTaxa[k] == NO)
                    continue;
                len = (int) strlen (taxaNames[k]);
                if (len > longestName)
                    longestName = len;
                }
            if (sumtParams.table == YES)
                {
                MrBayesPrint ("                                                                                   \n");
                MrBayesPrint ("%s   General explanation:                                                          \n", spacer);
                MrBayesPrint ("                                                                                   \n");
                MrBayesPrint ("%s   In an unrooted tree, a taxon bipartition (split) is specified by removing a   \n", spacer);
                MrBayesPrint ("%s   branch, thereby dividing the species into those to the left and those to the  \n", spacer);
                MrBayesPrint ("%s   right of the branch. Here, taxa to one side of the removed branch are denoted \n", spacer);
                MrBayesPrint ("%s   '.' and those to the other side are denoted '*'. Specifically, the '.' symbol \n", spacer);
                MrBayesPrint ("%s   is used for the taxa on the same side as the outgroup.                        \n", spacer);
                MrBayesPrint ("                                                                                   \n");
                MrBayesPrint ("%s   In a rooted or clock tree, the tree is rooted using the model and not by      \n", spacer);
                MrBayesPrint ("%s   reference to an outgroup. Each bipartition therefore corresponds to a clade,  \n", spacer);
                MrBayesPrint ("%s   that is, a group that includes all the descendants of a particular branch in  \n", spacer);
                MrBayesPrint ("%s   the tree.  Taxa that are included in each clade are denoted using '*', and    \n", spacer);
                MrBayesPrint ("%s   taxa that are not included are denoted using the '.' symbol.                  \n", spacer);
                MrBayesPrint ("                                                                                   \n");
                MrBayesPrint ("%s   The output first includes a key to all the bipartitions with frequency larger \n", spacer);
                MrBayesPrint ("%s   or equual to (Minpartfreq) in at least one run. Minpartfreq is a parameter to \n", spacer);
                MrBayesPrint ("%s   sumt command and currently it is set to %1.2lf.  This is followed by a table  \n", spacer, sumtParams.minPartFreq);
                MrBayesPrint ("%s   with statistics for the informative bipartitions (those including at least    \n", spacer);
                MrBayesPrint ("%s   two taxa), sorted from highest to lowest probability. For each bipartition,   \n", spacer);
                MrBayesPrint ("%s   the table gives the number of times the partition or split was observed in all\n", spacer);
                MrBayesPrint ("%s   runs (#obs) and the posterior probability of the bipartition (Probab.), which \n", spacer);
                MrBayesPrint ("%s   is the same as the split frequency. If several runs are summarized, this is   \n", spacer);
                MrBayesPrint ("%s   followed by the minimum split frequency (Min(s)), the maximum frequency       \n", spacer);
                MrBayesPrint ("%s   (Max(s)), and the standard deviation of frequencies (Stddev(s)) across runs.  \n", spacer);
                MrBayesPrint ("%s   The latter value should approach 0 for all bipartitions as MCMC runs converge.\n", spacer);
                MrBayesPrint ("                                                                                   \n");
                MrBayesPrint ("%s   This is followed by a table summarizing branch lengths, node heights (if a    \n", spacer);
                MrBayesPrint ("%s   clock model was used) and relaxed clock parameters (if a relaxed clock model  \n", spacer);
                MrBayesPrint ("%s   was used). The mean, variance, and 95 %% credible interval are given for each \n", spacer);
                MrBayesPrint ("%s   of these parameters. If several runs are summarized, the potential scale      \n", spacer);
                MrBayesPrint ("%s   reduction factor (PSRF) is also given; it should approach 1 as runs converge. \n", spacer);
                MrBayesPrint ("%s   Node heights will take calibration points into account, if such points were   \n", spacer);
                MrBayesPrint ("%s   used in the analysis.                                                         \n", spacer);
                MrBayesPrint ("%s                                                                                 \n", spacer);
                MrBayesPrint ("%s   Note that Stddev may be unreliable if the partition is not present in all     \n", spacer);
                MrBayesPrint ("%s   runs (the last column indicates the number of runs that sampled the partition \n", spacer);
                MrBayesPrint ("%s   if more than one run is summarized). The PSRF is not calculated at all if     \n", spacer); 
                MrBayesPrint ("%s   the partition is not present in all runs.The PSRF is also sensitive to small  \n", spacer);
                MrBayesPrint ("%s   sample sizes and it should only be considered a rough guide to convergence    \n", spacer);
                MrBayesPrint ("%s   since some of the assumptions allowing one to interpret it as a true potential\n", spacer);
                MrBayesPrint ("%s   scale reduction factor are violated in MrBayes.                               \n", spacer);
                MrBayesPrint ("%s                                                                                 \n", spacer);
                MrBayesPrint ("%s   List of taxa in bipartitions:                                                 \n", spacer);
                MrBayesPrint ("                                                                                   \n");
                j = 1;
                for (k=0; k<numTaxa; k++)
                    {
                    if (taxaInfo[k].isDeleted == NO && sumtParams.absentTaxa[k] == NO)
                        {
                        MrBayesPrint ("%s   %4d -- %s\n", spacer, j++, taxaNames[k]);
                        }
                    }
                }
            }
    
        if (sumtParams.numTrees > 1 && (sumtParams.table == YES || sumtParams.summary == YES))
            {
            MrBayesPrint ("\n\n");
            MrBayesPrint ("%s   Results for tree number %d\n", spacer, treeNo+1);
            MrBayesPrint ("%s   ==========================\n\n", spacer);
            }

        /* First print key to taxon bipartitions */
        if (sumtParams.table == YES)
            {
            MrBayesPrint ("\n");
            if (sumtParams.numTrees == 1)
                MrBayesPrint ("%s   Key to taxon bipartitions (saved to file \"%s.parts\"):\n\n", spacer, sumtParams.sumtOutfile);
            else
                MrBayesPrint ("%s   Key to taxon bipartitions (saved to file \"%s.tree%d.parts\"):\n\n", spacer,  sumtParams.sumtOutfile, treeNo+1);
            }

        /* calculate a couple of numbers that are handy to have */
        /*numTreePartsToPrint = 0;
        for (i=0; i<numUniqueSplitsFound; i++)
            {
            if ((MrBFlt)treeParts[i]->totCount/(MrBFlt)sumtParams.numTreesSampled >= sumtParams.minPartFreq)
                numTreePartsToPrint++;
            }
            */
        maxWidthID = (int) (log10 (numTreePartsToPrint)) + 1;
        if (maxWidthID < 2)
            maxWidthID = 2;
        maxNumTaxa = SCREEN_WIDTH - 9;

        for (j=0; j<sumtParams.numTaxa; j+=maxNumTaxa)
            {
            /* print header to screen and to parts file simultaneously */

            /* first print header to screen */
            MrBayesPrint ("%s   ", spacer);
            if (j == 0)
                MrBayesPrint ("%*s -- Partition\n", maxWidthID, "ID");
            else
                MrBayesPrint ("%*s -- Partition (continued)\n", maxWidthID, "ID");
            tableWidth = maxWidthID + 4 + sumtParams.numTaxa;
            if (tableWidth > SCREEN_WIDTH)
                tableWidth = SCREEN_WIDTH;
            MrBayesPrint ("%s   ", spacer);
            for (i=0; i<tableWidth; i++)
                MrBayesPrint ("-");
            MrBayesPrint ("\n");

            /* now print header to file */
            MrBayesPrintf (fpParts, "ID\tPartition\n");

            /* now, show partitions that were found on screen; print to .parts file simultaneously */
            mask = SafeCalloc (sumtParams.BitsLongsNeeded, sizeof(BitsLong));
            for (i=0; i<sumtParams.numTaxa; i++)
                SetBit (i, mask);
            for (i=0; i<numTreePartsToPrint; i++)
                {
                x = treeParts[i];
                if (IsBitSet(localOutGroup, x->partition) == YES && sumtParams.isRooted == NO)
                    FlipBits(x->partition, sumtParams.BitsLongsNeeded, mask);

                if ((NumBits(x->partition, sumtParams.BitsLongsNeeded) == numLocalTaxa || NumBits(x->partition, sumtParams.BitsLongsNeeded) == 0) && sumtParams.isClock == NO)
                    continue;

                MrBayesPrint ("%s   %*d -- ", spacer, maxWidthID, i);
                if (sumtParams.numTaxa <= maxNumTaxa)
                    ShowParts (stdout, x->partition, sumtParams.numTaxa);
                else
                    {
                    if (sumtParams.numTaxa - j > maxNumTaxa)
                        ShowSomeParts (stdout, x->partition, j, maxNumTaxa);
                    else
                        ShowSomeParts (stdout, x->partition, j, sumtParams.numTaxa - j);
                    }
                fflush(stdout);
                MrBayesPrint ("\n");

                MrBayesPrintf (fpParts, "%d\t", i);
                ShowParts (fpParts, x->partition, sumtParams.numTaxa);
                MrBayesPrintf (fpParts, "\n");
                }
            free (mask);

            /* finish screen table */
            if (sumtParams.table == YES)
                {
                MrBayesPrint ("%s   ", spacer);
                for (i=0; i<tableWidth; i++)
                    {
                    MrBayesPrint ("-");
                    }
                MrBayesPrint ("\n");
                if (oneUnreliable == YES)
                    {
                    MrBayesPrint ("%s   * The partition was not found in all runs so the values are unreliable\n\n", spacer);
                    }
                else
                    {
                    MrBayesPrint ("\n");
                    }
                }
            }

        /* Second, print statitistics for taxon bipartitions */
        if (sumtParams.table == YES)
            {
            if (sumtParams.isRooted == NO)
                MrBayesPrint ("%s   Summary statistics for informative taxon bipartitions\n", spacer);
            else
                MrBayesPrint ("%s   Summary statistics for informative taxon bipartitions (clades)\n", spacer);
            MrBayesPrint ("%s      (saved to file \"%s.tstat\"):\n\n", spacer, sumtParams.sumtOutfile);
            }

        /* calculate a couple of numbers that are handy to have */
        /*numTreePartsToPrint = 0;
        for (i=0; i<numUniqueSplitsFound; i++)
            {
            if ((MrBFlt)treeParts[i]->totCount/(MrBFlt)sumtParams.numTreesSampled >= sumtParams.minPartFreq)
                numTreePartsToPrint++;
            }
            */
        maxWidthID = (int) (log10 (numTreePartsToPrint)) + 1;
        if (maxWidthID < 2)
            maxWidthID = 2;
        maxWidthNumberPartitions = (int) (log10 (treeParts[0]->totCount)) + 1;
        if (maxWidthNumberPartitions < 4)
            maxWidthNumberPartitions = 4;

        /* print header to screen and to parts file simultaneously */
        if (sumtParams.table == YES)
            {
            /* first print header to screen */
            MrBayesPrint ("%s   ", spacer);
            MrBayesPrint ("%*s   ", maxWidthID, "ID");
            tableWidth = maxWidthID + 3;
            MrBayesPrint ("#obs");
            tableWidth += 4;
            for (i=4; i<maxWidthNumberPartitions; i++)
                {
                MrBayesPrint (" ");
                tableWidth++;
                }
            MrBayesPrint ("    Probab.");
            tableWidth += 11;
            if (sumtParams.numRuns > 1)
                {
                MrBayesPrint ("     Sd(s)+ ");
                MrBayesPrint ("     Min(s)      Max(s) ");
                tableWidth += 36;
                MrBayesPrint ("  Nruns ");
                tableWidth += 8;
                }
            MrBayesPrint ("\n%s   ", spacer);
            for (i=0; i<tableWidth; i++)
                {
                MrBayesPrint ("-");
                }
            MrBayesPrint ("\n");

            /* now print header to file */
            if (sumtParams.numRuns > 1)
                MrBayesPrintf (fpTstat, "ID\t#obs\tProbability(=s)\tStddev(s)\tMin(s)\tMax(s)\tNruns\n");
            else
                MrBayesPrintf (fpTstat, "ID\t#obs\tProbability(=s)\n");
            }

        /* now, show informative partitions that were found on screen; print to .tstat file simultaneously */
        for (i=0; i<numTreePartsToPrint; i++)
            {
            x = treeParts[i];

            /* skip uninformative partitions */
            if (NumBits(x->partition, sumtParams.BitsLongsNeeded) <= 1 || NumBits(x->partition, sumtParams.BitsLongsNeeded) == sumtParams.numTaxa)
                continue;
            if (NumBits(x->partition, sumtParams.BitsLongsNeeded) == sumtParams.numTaxa - 1 && sumtParams.isRooted == NO)
                continue;

            if (sumtParams.table == YES)
                {
                MrBayesPrint ("%s   %*d", spacer, maxWidthID, i);
                fflush(stdout);
                MrBayesPrintf (fpTstat, "%d\t", i);
                }
            if (sumtParams.numRuns > 1)
                {
                sum_s = 0.0;
                sumsq_s = 0.0;
                min_s = 1.0;
                max_s = 0.0;
                for (n=j=0; n<sumtParams.numRuns; n++)
                    {
                    if (x->count[n] > 0)
                        j++;
                    f = (MrBFlt) x->count[n] / (MrBFlt) sumtParams.numFileTreesSampled[n];
                    sum_s += f;
                    sumsq_s += f * f;
                    if (f < min_s)
                        min_s = f;
                    if (f  > max_s)
                        max_s = f;
                    }
                var_s = sumsq_s - sum_s * sum_s / (MrBFlt) sumtParams.numRuns;
                var_s /= (sumtParams.numRuns - 1);
                if (var_s > 0.0)
                    stddev_s = sqrt (var_s);
                else
                    stddev_s = 0.0;
                if (j == sumtParams.numRuns)
                    unreliable = NO;
                else
                    {
                    unreliable = YES;
                    oneUnreliable = YES;
                    }
                }
            if (sumtParams.table == YES)
                {
                f = (MrBFlt) x->totCount / (MrBFlt) sumtParams.numTreesSampled;
                MrBayesPrint ("  %*d    %1.6lf", maxWidthNumberPartitions, x->totCount, f);
                MrBayesPrintf (fpTstat, "\t%d\t%s", x->totCount, MbPrintNum(f));
                if (sumtParams.numRuns > 1)
                    {
                    MrBayesPrint ("    %1.6lf    %1.6lf    %1.6lf", stddev_s, min_s, max_s);
                    MrBayesPrint ("  %3d", j);
                    MrBayesPrintf (fpTstat, "\t%s", MbPrintNum(stddev_s));
                    MrBayesPrintf (fpTstat, "\t%s", MbPrintNum(min_s));
                    MrBayesPrintf (fpTstat, "\t%s", MbPrintNum(max_s));
                    MrBayesPrintf (fpTstat, "\t%d", j);
                    }
                MrBayesPrintf (fpTstat, "\n");
                if (unreliable == YES)
                    MrBayesPrint (" *\n");
                else
                    MrBayesPrint ("\n");
                sumStdDev += stddev_s;
                if (stddev_s > maxStdDev)
                    maxStdDev = stddev_s;
                }
            }

        /* finish screen table */
        if (sumtParams.table == YES)
            {
            MrBayesPrint ("%s   ", spacer);
            for (i=0; i<tableWidth; i++)
                {
                MrBayesPrint ("-");
                }
            MrBayesPrint ("\n");
            if (sumtParams.numRuns > 1)
                {
                MrBayesPrint ("%s   + Convergence diagnostic (standard deviation of split frequencies)\n", spacer);
                MrBayesPrint ("%s     should approach 0.0 as runs converge.\n\n", spacer);
                }
            if (oneUnreliable == YES)
                MrBayesPrint ("%s   * The partition was not found in all runs so the values are unreliable\n", spacer);
            }

        /* Third, print statitistics for branch and node parameters */
        if (sumtParams.table == YES)
            {
            MrBayesPrint ("\n");
            MrBayesPrint ("%s   Summary statistics for branch and node parameters\n", spacer);
            MrBayesPrint ("%s      (saved to file \"%s.vstat\"):\n", spacer, sumtParams.sumtOutfile);
            }
        
        if (sumtParams.table == YES)
            {
            /* calculate longest header */
            longestHeader = 9;  /* length of 'parameter' */
            i = (int)(log10(numTreePartsToPrint)) + 3;   /* length of partition specifier including [] */
            len = i + (int)(strlen(treeName)) + 2;   /* length of length{m}[n] or height{m}[n] */
            if (len > longestHeader)
                longestHeader = len;
            for (j=0; j<sumtParams.nBSets; j++)
                {
                len = (int) strlen(sumtParams.tree->bSetName[j]) + 7 + i;
                if (len > longestHeader)
                    longestHeader = len;
                }
            for (j=0; j<sumtParams.nESets; j++)
                {
                len = (int) strlen(sumtParams.tree->eSetName[j]) + 8 + i;
                if (len > longestHeader)
                    longestHeader = len;
                }
            if (sumtParams.popSizeSet == YES)
                {
                len = (int) strlen(sumtParams.tree->popSizeSetName) + i;
                if (len > longestHeader)
                    longestHeader = len;
                }
    
            /* print the header rows */
            MrBayesPrint ("\n");
            if (sumtParams.HPD == NO)
                MrBayesPrint ("%s   %*c                             95%% Cred. Interval\n", spacer, longestHeader, ' ');
            else
                MrBayesPrint ("%s   %*c                              95%% HPD Interval\n", spacer, longestHeader, ' ');
            MrBayesPrint ("%s   %*c                            --------------------\n", spacer, longestHeader, ' ');

            MrBayesPrint ("%s   Parameter%*c     Mean       Variance     Lower       Upper       Median", spacer, longestHeader-9, ' ');
            tableWidth = 68 + longestHeader - 9;
            if (sumtParams.HPD == YES)
                MrBayesPrintf (fpVstat, "Parameter\tMean\tVariance\tCredInt_Lower\tCredInt_Upper\tMedian", spacer, longestHeader-9, ' ');
            else
                MrBayesPrintf (fpVstat, "Parameter\tMean\tVariance\tHPD_Lower\tHPD_Upper\tMedian", spacer, longestHeader-9, ' ');
            if (sumtParams.numRuns > 1)
                {
                    MrBayesPrint ("     PSRF+  Nruns");
                tableWidth += 17;
                MrBayesPrintf (fpVstat, "\tPSRF\tNruns");
                }
            MrBayesPrint ("\n");
            MrBayesPrintf (fpVstat, "\n");

            MrBayesPrint ("%s   ", spacer);
            for (j=0; j<tableWidth; j++)
                {
                MrBayesPrint ("-");
                }
            MrBayesPrint ("\n");

            /* print lengths */
            strcpy (divString, treeName+4);
            for (i=1; i<numTreePartsToPrint; i++)
                {
                x = treeParts[i];
                tempStrLength=(int)strlen(tempStr);
                SafeSprintf (&tempStr,&tempStrLength, "length%s[%d]", divString, i);
                len = (int) strlen(tempStr);

                GetSummary (x->length, sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);

                MrBayesPrint ("%s   %-*s  ", spacer, longestHeader, tempStr);
                MrBayesPrintf (fpVstat, "%s", tempStr);

                PrintSumtTableLine(sumtParams.numRuns, x->count, &theStats, &numPSRFSamples, &maxPSRF, &sumPSRF);

                }

            /* print heights */
           if (sumtParams.isClock == YES)
                {
                strcpy (divString, treeName+4);
                for (i=0; i<numTreePartsToPrint; i++)
                    {
                    x = treeParts[i];
                    tempStrLength=(int)strlen(tempStr);
                    SafeSprintf (&tempStr,&tempStrLength, "height%s[%d]", divString, i);
                    len = (int) strlen(tempStr);

                    GetSummary (x->height, sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);

                    MrBayesPrint ("%s   %-*s  ", spacer, longestHeader, tempStr);
                    MrBayesPrintf (fpVstat, "%s", tempStr);

                    PrintSumtTableLine(sumtParams.numRuns, x->count, &theStats, &numPSRFSamples, &maxPSRF, &sumPSRF);
                    }
                }

            /* print ages */
            if (sumtParams.isCalibrated == YES)
                {
                strcpy (divString, treeName+4);
                for (i=0; i<numTreePartsToPrint; i++)
                    {
                    x = treeParts[i];
                    tempStrLength=(int)strlen(tempStr);
                    SafeSprintf (&tempStr,&tempStrLength, "age%s[%d]", divString, i);
                    len = (int) strlen(tempStr);

                    GetSummary (x->age, sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);

                    MrBayesPrint ("%s   %-*s  ", spacer, longestHeader, tempStr);
                    MrBayesPrintf (fpVstat, "%s", tempStr);

                    PrintSumtTableLine(sumtParams.numRuns, x->count, &theStats, &numPSRFSamples, &maxPSRF, &sumPSRF);
                    }
                }

            /* print effective branch lengths */
            if (sumtParams.isRelaxed == YES)
                {
                for (i=0; i<sumtParams.nBSets; i++)
                    {
                    for (j=1; j<numTreePartsToPrint; j++)
                        {
                        x = treeParts[j];
                        tempStrLength=(int)strlen(tempStr);
                        SafeSprintf (&tempStr,&tempStrLength, "%s_length[%d]", sumtParams.bSetName[i], j);
                        len = (int) strlen(tempStr);

                        GetSummary (x->bLen[i], sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);

                        MrBayesPrint ("%s   %-*s  ", spacer, longestHeader, tempStr);
                        MrBayesPrintf (fpVstat, "%s", tempStr);

                        PrintSumtTableLine(sumtParams.numRuns, x->count, &theStats, &numPSRFSamples, &maxPSRF, &sumPSRF);
                        }
                    for (j=1; j<numTreePartsToPrint; j++)
                        {
                        x = treeParts[j];
                        tempStrLength=(int)strlen(tempStr);
                        SafeSprintf (&tempStr,&tempStrLength, "%s_rate[%d]", sumtParams.bSetName[i], j);
                        len = (int) strlen(tempStr);

                        GetSummary (x->bRate[i], sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);

                        MrBayesPrint ("%s   %-*s  ", spacer, longestHeader, tempStr);
                        MrBayesPrintf (fpVstat, "%s", tempStr);

                        PrintSumtTableLine(sumtParams.numRuns, x->count, &theStats, &numPSRFSamples, &maxPSRF, &sumPSRF);
                        }
                    }
                for (i=0; i<sumtParams.nESets; i++)
                    {
                    for (j=1; j<numTreePartsToPrint; j++)
                        {
                        x = treeParts[j];
                        tempStrLength=(int)strlen(tempStr);
                        SafeSprintf (&tempStr,&tempStrLength, "%s_nEvents[%d]", sumtParams.eSetName[i], j);
                        len = (int) strlen(tempStr);

                        GetIntSummary (x->nEvents[i], sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);

                        MrBayesPrint ("%s   %-*s  ", spacer, longestHeader, tempStr);
                        MrBayesPrintf (fpVstat, "%s", tempStr);

                        PrintSumtTableLine(sumtParams.numRuns, x->count, &theStats, &numPSRFSamples, &maxPSRF, &sumPSRF);
                        }
                    }
                }

            /* print population size sets */
            if (sumtParams.popSizeSet == YES)
                {
                for (j=1; j<numTreePartsToPrint; j++)
                    {
                    x = treeParts[j];
                    tempStrLength=(int)strlen(tempStr);
                    SafeSprintf (&tempStr,&tempStrLength, "%s[%d]", sumtParams.popSizeSetName, j);
                    len = (int) strlen(tempStr);

                    GetSummary (x->popSize, sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);

                    MrBayesPrint ("%s   %-*s  ", spacer, longestHeader, tempStr);
                    MrBayesPrintf (fpVstat, "%s", tempStr);

                    PrintSumtTableLine(sumtParams.numRuns, x->count, &theStats, &numPSRFSamples, &maxPSRF, &sumPSRF);
                    }
                }

            /* finish table */
            MrBayesPrint ("%s   ", spacer);
            for (j=0; j<tableWidth; j++)
                MrBayesPrint ("-");
            MrBayesPrint ("\n");

            if (sumtParams.numRuns > 1)
                {
                MrBayesPrint ("%s   + Convergence diagnostic (PSRF = Potential Scale Reduction Factor; Gelman\n", spacer);
                MrBayesPrint ("%s     and Rubin, 1992) should approach 1.0 as runs converge. NA is reported when\n", spacer);
                MrBayesPrint ("%s     deviation of parameter values within all runs is 0 or when a parameter\n", spacer);
                MrBayesPrint ("%s     value (a branch length, for instance) is not sampled in all runs.\n", spacer);
                }

            if (oneUnreliable == YES)
                {
                MrBayesPrint ("%s   * The partition was not found in all runs so the values are unreliable.\n", spacer);
                }
            MrBayesPrint ("\n\n");
            }
            
            /* Print branch parameters to file if requested */
            strcpy( divString, treeName+4);
            if (sumtParams.saveBrParams == YES && PrintBrParamsToFile (treeParts, numUniqueSplitsFound, treeNo, divString) == ERROR)
                goto errorExit;
            
            /* Exclude trivial splits when calculating average standard deviation of split frequencies. */
            avgStdDev = sumStdDev / (numTreePartsToPrint-sumtParams.numTaxa-1);
            avgPSRF   = sumPSRF / numPSRFSamples;

            if (sumtParams.numRuns > 1 && sumtParams.summary == YES)
                {
                MrBayesPrint ("%s   Summary statistics for partitions with frequency >= %1.2lf in at least one run:\n", spacer, sumtParams.minPartFreq);
                MrBayesPrint ("%s       Average standard deviation of split frequencies = %1.6lf\n", spacer, avgStdDev);
                MrBayesPrint ("%s       Maximum standard deviation of split frequencies = %1.6lf\n", spacer, maxStdDev);
                }
            if (sumtParams.brlensDef == YES && sumtParams.numRuns > 1 && sumtParams.summary == YES)
                {
                MrBayesPrint ("%s       Average PSRF for parameter values (excluding NA and >10.0) = %1.3lf\n", spacer, avgPSRF);
                if (maxPSRF == 10)
                    MrBayesPrint ("%s       Maximum PSRF for parameter values = NA\n", spacer);
                else
                    MrBayesPrint ("%s       Maximum PSRF for parameter values = %1.3lf\n", spacer, maxPSRF);
                }
            MrBayesPrint ("\n");

        SortPartCtr (treeParts, 0, numUniqueSplitsFound-1); /* We sort again but this time we sort all partitions instead of just first numTreePartsToPrintNow */
        /* make the majority rule consensus tree */
        if (sumtParams.showConsensus == YES && ConTree (treeParts, numUniqueSplitsFound) == ERROR) 
            goto errorExit;
            
        /* get probabilities of individual trees */
        if (TreeProb () == ERROR)
            goto errorExit;
        
        /* close files */
        SafeFclose (&fpParts);
        SafeFclose (&fpTstat);
        SafeFclose (&fpVstat);
        SafeFclose (&fpCon);
        SafeFclose (&fpTrees);

#   if defined (PRINT_RATEMUL_CPP)
        SafeFclose (&rateMultfp);
#   endif

        /* free pointer array to partitions */
        free (treeParts);
        treeParts = NULL;
        FreePartCtr (partCtrRoot);
        partCtrRoot = NULL;
        FreeTreeCtr (treeCtrRoot);
        treeCtrRoot = NULL;
        } /* next tree */

    /* free memory and file pointers */
    if (s) free(s);
    FreeSumtParams();

    /* reset numLocalTaxa and localOutGroup */
    ResetTaxonSet();

#   if defined (MPI_ENABLED)
        }
#   endif

    expecting = Expecting(COMMAND);
    inSumtCommand = NO;
    SafeFree ((void **)&tempStr);

    return (NO_ERROR);
    
    /* error exit */
    errorExit:
        /* free sumtParams */
        if (s) free(s);
        FreeSumtParams();
        
        /* close files in case they are open*/
        SafeFclose (&fp);
        SafeFclose (&fpParts);
        SafeFclose (&fpTstat);
        SafeFclose (&fpVstat);
        SafeFclose (&fpCon);
        SafeFclose (&fpTrees);

#   if defined (PRINT_RATEMUL_CPP)
        SafeFclose (&rateMultfp);
#   endif

        /* free pointer array to partitions, part and tree counters */
        free (treeParts);
        FreePartCtr (partCtrRoot);
        FreeTreeCtr (treeCtrRoot);
        partCtrRoot = NULL;
        treeCtrRoot = NULL;

        /* reset taxon set */
        ResetTaxonSet();

        expecting = Expecting(COMMAND);
        inSumtCommand = NO;
        SafeFree ((void **)&tempStr);

        return (ERROR);
}


int DoSumtParm (char *parmName, char *tkn)
{
    int         tempI;
    MrBFlt      tempD;
    char        tempStr[100];
    
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
        /* set Filename (sumtParams.sumtFileName) ***************************************************/
        else if (!strcmp(parmName, "Filename"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                if (strlen(tkn)>99)
                    {
                    MrBayesPrint ("%s   Maximum allowed length of file name is 99 characters. The given name:\n", spacer);
                    MrBayesPrint ("%s      '%s'\n", spacer,tkn);
                    MrBayesPrint ("%s   has %d characters.\n", spacer,strlen(tkn));
                    return (ERROR);
                    }
                strcpy (sumtParams.sumtFileName, tkn);
                strcpy(sumtParams.sumtOutfile, tkn);
                MrBayesPrint ("%s   Setting sumt filename and outputname to %s\n", spacer, sumtParams.sumtFileName);
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
                MrBayesPrint ("%s   Setting urn-in to %d\n", spacer, chainParams.chainBurnIn);
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
        /* set Nruns (sumtParams.numRuns) *******************************************************/
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
                    sumtParams.numRuns = tempI;
                    MrBayesPrint ("%s   Setting sumt nruns to %d\n", spacer, sumtParams.numRuns);
                    expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                    }
                }
            else
                return (ERROR);
            }
        /* set Ntrees (sumtParams.numTrees) *******************************************************/
        else if (!strcmp(parmName, "Ntrees"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   Ntrees must be at least 1\n", spacer);
                    return (ERROR);
                    }
                else
                    {
                    sumtParams.numTrees = tempI;
                    MrBayesPrint ("%s   Setting sumt ntrees to %d\n", spacer, sumtParams.numTrees);
                    expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                    }
                }
            else
                return (ERROR);
            }
        /* set Contype (sumtParams.sumtConType) *****************************************************/
        else if (!strcmp(parmName, "Contype"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    strcpy (sumtParams.sumtConType, tempStr);
                    MrBayesPrint ("%s   Setting sumt contype to %s\n", spacer, sumtParams.sumtConType);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Conformat (sumtParams.consensusFormat) *****************************************************/
        else if (!strcmp(parmName, "Conformat"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr,"Figtree"))
                        {
                        MrBayesPrint ("%s   Setting sumt conformat to Figtree\n", spacer);
                        sumtParams.consensusFormat = FIGTREE;
                        }
                    else
                        {
                        MrBayesPrint ("%s   Setting sumt conformat to Simple\n", spacer);
                        sumtParams.consensusFormat = SIMPLE;
                        }
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for calctreeprobs\n", spacer);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Calctreeprobs (sumtParams.calcTreeprobs) *********************************************/
        else if (!strcmp(parmName, "Calctreeprobs"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        sumtParams.calcTreeprobs = YES;
                    else
                        sumtParams.calcTreeprobs = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for calctreeprobs\n", spacer);
                    return (ERROR);
                    }
                if (sumtParams.calcTreeprobs == YES)
                    MrBayesPrint ("%s   Setting calctreeprobs to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting calctreeprobs to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Showtreeprobs (sumtParams.showSumtTrees) *********************************************/
        else if (!strcmp(parmName, "Showtreeprobs"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        sumtParams.showSumtTrees = YES;
                    else
                        sumtParams.showSumtTrees = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for showtreeprobs\n", spacer);
                    return (ERROR);
                    }
                if (sumtParams.showSumtTrees == YES)
                    MrBayesPrint ("%s   Setting showtreeprobs to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting showtreeprobs to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
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
                        sumtParams.HPD = YES;
                    else
                        sumtParams.HPD = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Hpd\n", spacer);
                    return (ERROR);
                    }
                if (sumtParams.HPD == YES)
                    MrBayesPrint ("%s   Reporting 95 %% region of Highest Posterior Density (HPD).\n", spacer);
                else
                    MrBayesPrint ("%s   Reporting median interval containing 95 %% of values.\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Savebrparams (sumtParams.saveBrParams) *********************************************/
        else if (!strcmp(parmName, "Savebrparams"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        sumtParams.saveBrParams = YES;
                    else
                        sumtParams.saveBrParams = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for savebrparams\n", spacer);
                    return (ERROR);
                    }
                if (sumtParams.saveBrParams == YES)
                    MrBayesPrint ("%s   Setting savebrparams to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting savebrparams to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Minbrparamfreq (sumtParams.minBrParamFreq) *******************************************************/
        else if (!strcmp(parmName, "Minbrparamfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                sumtParams.minBrParamFreq = tempD;
                MrBayesPrint ("%s   Printing branch parameters to file for partitions with probability >= %lf\n", spacer, sumtParams.minBrParamFreq);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Ordertaxa (sumtParams.orderTaxa) *********************************************/
        else if (!strcmp(parmName, "Ordertaxa"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        sumtParams.orderTaxa = YES;
                    else
                        sumtParams.orderTaxa = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for ordertaxa\n", spacer);
                    return (ERROR);
                    }
                if (sumtParams.orderTaxa == YES)
                    MrBayesPrint ("%s   Setting ordertaxa to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting ordertaxa to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Outputname (sumtParams.sumtOutfile) *******************************************************/
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
                strcpy (sumtParams.sumtOutfile, tempStr);
                MrBayesPrint ("%s   Setting sumt output file name to \"%s\"\n", spacer, sumtParams.sumtOutfile);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Table (sumtParams.table) ********************************************************/
        else if (!strcmp(parmName, "Table"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        sumtParams.table = YES;
                    else
                        sumtParams.table = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Table (valid arguments are 'yes' and 'no')\n", spacer);
                    return (ERROR);
                    }
                if (sumtParams.table == YES)
                    MrBayesPrint ("%s   Setting sumt to compute table of partition frequencies\n", spacer);
                else
                    MrBayesPrint ("%s   Setting sumt not to compute table of partition frequencies\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Summary (sumtParams.summary) ********************************************************/
        else if (!strcmp(parmName, "Summary"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        sumtParams.summary = YES;
                    else
                        sumtParams.summary = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for 'Summary' (valid arguments are 'yes' and 'no')\n", spacer);
                    return (ERROR);
                    }
                if (sumtParams.summary == YES)
                    MrBayesPrint ("%s   Setting sumt to summary statistics\n", spacer);
                else
                    MrBayesPrint ("%s   Setting sumt not to compute summary statistics\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Consensus (sumtParams.showConsensus) ********************************************************/
        else if (!strcmp(parmName, "Consensus"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        sumtParams.showConsensus = YES;
                    else
                        sumtParams.showConsensus = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Consensus (valid arguments are 'yes' and 'no')\n", spacer);
                    return (ERROR);
                    }
                if (sumtParams.showConsensus == YES)
                    MrBayesPrint ("%s   Setting sumt to show consensus trees\n", spacer);
                else
                    MrBayesPrint ("%s   Setting sumt not to show consensus trees\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Minpartfreq (sumtParams.minPartFreq) *******************************************************/
        else if (!strcmp(parmName, "Minpartfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                sumtParams.minPartFreq = tempD;
                MrBayesPrint ("%s   Including partitions with probability greater than or equal to %lf in summary statistics\n", spacer, sumtParams.minPartFreq);
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


int DoSumtTree (void)
{
    int             i, j, z, printEvery, nAstPerPrint, burnin;
    MrBFlt          x, y;
    PolyTree        *t;
    PolyNode        *p;

#   if defined (PRINT_RATEMUL_CPP)
    /* get depths if relevant */
    if (sumtParams.tree->isClock)
        GetPolyDepths (sumtParams.tree);
    if (rateMultfp !=NULL && sumtParams.tree->root !=NULL)
        DELETE_ME_dump_depth(sumtParams.tree->root);
#   endif

    /* increment number of trees read in */
    sumtParams.numFileTrees[sumtParams.runId]++;
    sumtParams.numTreesEncountered++;

    /*  update status bar */
    if (sumtParams.numTreesInLastBlock * sumtParams.numRuns < 80)
        {
        printEvery = 1;
        nAstPerPrint = 80 / (sumtParams.numTreesInLastBlock * sumtParams.numRuns);
        if (sumtParams.numTreesEncountered % printEvery == 0)
            {
            for (i=0; i<nAstPerPrint; i++)
                {
                MrBayesPrint ("*");
                numAsterices++;
                }
            }
        }
    else
        {
        x = (MrBFlt)(sumtParams.numTreesInLastBlock * sumtParams.numRuns) / (MrBFlt) (80);
        y = (MrBFlt)(sumtParams.numFileTrees[sumtParams.runId] + sumtParams.numTreesInLastBlock * sumtParams.runId) / x;
        z = (int)y;
        if (numAsterices < z)
            {
            MrBayesPrint ("*");
            numAsterices++;
            }
        }
    
    /* get burnin */
    if (inComparetreeCommand == YES)
        burnin = comptreeParams.burnin;
    else
        burnin = sumtParams.burnin;

    if (sumtParams.numFileTrees[sumtParams.runId] > burnin)
        {
        /* increment the number of trees sampled */
        sumtParams.numFileTreesSampled[sumtParams.runId]++;
        sumtParams.numTreesSampled++;

        /* get the tree we just read in */
        t = sumtParams.tree;
        
        /* move calculation root for nonrooted trees if necessary */
        MovePolyCalculationRoot (t, localOutGroup);
        
        /* check taxon set and outgroup */
        if (sumtParams.runId == 0 && sumtParams.numFileTreesSampled[0] == 1)
            {
            if (isTranslateDef == YES && isTranslateDiff == YES)
                {
                /* we are using a translate block with different taxa set (not all taxa) */
                if (t->nNodes - t->nIntNodes != numTranslates)
                    {
                    MrBayesPrint ("%s   ERROR: Expected %d taxa; found %d taxa\n", spacer, numTranslates, sumtParams.numTaxa);
                    return (ERROR);
                    }
                /* figure out which taxa are absent */
                for (i=0; i<numTaxa; i++)
                    sumtParams.absentTaxa[i] = YES;
                for (i=0; i<numTranslates; i++)
                    {
                    for (j=0; j<numTaxa; j++)
                        {
                        if (strcmp(transFrom[i], taxaNames[j]) == 0)
                            break;
                        }
                    if (j < numTaxa)
                        sumtParams.absentTaxa[j] = NO;
                    }
                /* find local outgroup */
                for (i=0; i<numTranslates; i++)
                    if (strcmp(transFrom[i], taxaNames[outGroupNum]) == 0)
                        break;
                if (i == numTranslates)
                    localOutGroup = 0;      /* no previous outgroup assignment is valid */
                else
                    localOutGroup = i;
                }
            else
                {
                /* we are using the current taxa set */
                for (i=0; i<numTaxa; i++)
                    sumtParams.absentTaxa[i] = YES;
                for (i=0; i<t->nNodes; i++)
                    {
                    p = t->allDownPass[i];
                    if (p->left == NULL)
                        sumtParams.absentTaxa[p->index] = NO;
                    }
                localOutGroup = 0;
                for (i=j=0; i<numTaxa; i++)
                    {
                    if (sumtParams.absentTaxa[i] == NO && taxaInfo[i].isDeleted == NO)
                        {
                        if (i == outGroupNum)
                            localOutGroup = j;
                        j++;
                        }
                    }
                }

            /* now we can safely prune the tree based on taxaInfo[].isDeleted */
            /* note that PrunePolyTree relies on labels to recognize tips */
            PrunePolyTree(t);

            /* reset tip and int node indices in case some taxa deleted */
            ResetTipIndices (t);
            ResetIntNodeIndices(t);

            /* set basic parameters */
            sumtParams.numTaxa = t->nNodes - t->nIntNodes;
            numLocalTaxa = sumtParams.numTaxa;
            sumtParams.BitsLongsNeeded = ((numLocalTaxa-1) / nBitsInALong) + 1;
            if (t->isRooted == YES)
                sumtParams.orderLen = numLocalTaxa - 2;
            else
                sumtParams.orderLen = numLocalTaxa - 3;
            }
        else
            {
            /* check whether we have some taxa that should not be in the tree */
            for (i=0; i<t->nNodes; i++)
                {
                p = t->allDownPass[i];
                if (p->left == NULL && taxaInfo[p->index].isDeleted == NO && sumtParams.absentTaxa[p->index] == YES)
                    {
                    MrBayesPrint ("%s   Taxon %d should not be in sampled tree\n", spacer, p->index + 1);
                    return (ERROR);
                    }
                }

            /* now we can safely prune the tree based on taxaInfo[].isDeleted */
            PrunePolyTree (t);

            /* reset tip and int node indices in case some taxa deleted */
            ResetTipIndices (t);
            ResetIntNodeIndices(t);

            /* check that all taxa are included */
            if (t->nNodes - t->nIntNodes != sumtParams.numTaxa)
                {
                MrBayesPrint ("%s   Expecting %d taxa but tree '%s' in file '%s' has %d taxa\n",
                    spacer, sumtParams.numTaxa, t->name, sumtParams.curFileName, t->nNodes-t->nIntNodes);
                return ERROR;
                }
            }

        if (sumtParams.runId == 0 && sumtParams.numFileTreesSampled[0] == 1)
            {
            /* harvest labels (can only be done safely after pruning) */
            for (i=0; i<sumtParams.numTaxa; i++)
                {
                for (j=0; j<t->nNodes; j++)
                    {
                    p = t->allDownPass[j];
                    if (p->index == i) {
                         if (strlen(p->label)>99)
                            {
                            MrBayesPrint ("%s   Taxon name %s is too long. Maximun 99 characters is allowed.\n", spacer, p->label);
                            return (ERROR);
                            }
                        AddString(&sumtParams.taxaNames, i, p->label);
                        }
                    }
                }
            }

        /* check that tree agrees with template */
        if (sumtParams.numTreesSampled == 1)
            {
            sumtParams.brlensDef = t->brlensDef;
            sumtParams.isRooted = t->isRooted;
            sumtParams.isClock = t->isClock;
            sumtParams.isCalibrated = t->isCalibrated;
            sumtParams.isRelaxed = t->isRelaxed;
            sumtParams.nBSets = 0;
            sumtParams.nESets = 0;
            for (i=0; i<t->nBSets; i++)
                AddString(&sumtParams.bSetName,sumtParams.nBSets++,t->bSetName[i]);
            for (i=0; i<t->nESets; i++)
                AddString(&sumtParams.eSetName,sumtParams.nESets++,t->eSetName[i]);
            if (t->popSizeSet == YES)
                {
                sumtParams.popSizeSet = YES;
                sumtParams.popSizeSetName = (char *) SafeCalloc (strlen(t->popSizeSetName)+1, sizeof(char));
                strcpy(sumtParams.popSizeSetName, t->popSizeSetName);
                }
            else
                sumtParams.popSizeSet = NO;
            }
        else /* if (sumtParams.numTreesSampled > 1) */
            {
            if (sumtParams.brlensDef != t->brlensDef)
                {
                MrBayesPrint ("%s   Trees with and without branch lengths mixed\n", spacer);
                return ERROR;
                }
            if (sumtParams.isRooted != t->isRooted)
                {
                if (sumtParams.isRooted == YES)
                    MrBayesPrint ("%s   Expected rooted tree but tree '%s' in file '%s' is not rooted\n",
                        spacer, t->name, sumtParams.curFileName);
                else if (sumtParams.isRooted == NO)
                    MrBayesPrint ("%s   Expected unrooted tree but tree '%s' in file '%s' is rooted\n",
                        spacer, t->name, sumtParams.curFileName);
                return ERROR;
                }
            if (sumtParams.isClock != t->isClock)
                {
                if (sumtParams.isClock == YES)
                    MrBayesPrint ("%s   Expected clock tree but tree '%s' in file '%s' is not clock\n",
                        spacer, t->name, sumtParams.curFileName);
                else if (sumtParams.isClock == NO)
                    MrBayesPrint ("%s   Expected nonclock tree but tree '%s' in file '%s' is clock\n",
                        spacer, t->name, sumtParams.curFileName);
                return ERROR;
                }
            if (sumtParams.isCalibrated != t->isCalibrated)
                {
                if (sumtParams.isCalibrated == YES)
                    MrBayesPrint ("%s   Expected calibrated tree but tree '%s' in file '%s' is not calibrated\n",
                        spacer, t->name, sumtParams.curFileName);
                else if (sumtParams.isCalibrated == NO)
                    MrBayesPrint ("%s   Expected noncalibrated tree but tree '%s' in file '%s' is calibrated\n",
                        spacer, t->name, sumtParams.curFileName);
                return ERROR;
                }
            if (inComparetreeCommand == NO && sumtParams.isRelaxed != t->isRelaxed)
                {
                if (sumtParams.isRelaxed == YES)
                    MrBayesPrint ("%s   Expected relaxed clock tree but tree '%s' in file '%s' is not relaxed\n",
                        spacer, t->name, sumtParams.curFileName);
                else if (sumtParams.isRelaxed == NO)
                    MrBayesPrint ("%s   Expected unrooted tree but tree '%s' in file '%s' is rooted\n",
                        spacer, t->name, sumtParams.curFileName);
                return ERROR;
                }
            if (inComparetreeCommand == NO && (sumtParams.nESets != t->nESets || sumtParams.nBSets != t->nBSets))
                {
                MrBayesPrint ("%s   Tree '%s' in file '%s' does not have the expected relaxed clock parameters\n",
                        spacer, t->name, sumtParams.curFileName);
                return ERROR;
                }
            }

        /* set partitions for tree */
        ResetPolyTreePartitions(t);

        /* get depths if relevant */
        if (t->isClock)
            GetPolyDepths (t);

        /* get ages if relevant */
        if (t->isCalibrated)
            GetPolyAges (t);
        
        /* add partitions to counters */
        for (i=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            partCtrRoot = AddSumtPartition (partCtrRoot, t, p, sumtParams.runId);
            }
            
        /* add the tree to relevant tree list */
        if (inSumtCommand == YES)
            {
            if (t->isRooted == YES)
                StoreRPolyTopology (t, sumtParams.order);
            else /* if (sumtParams.isRooted == NO) */
                StoreUPolyTopology (t, sumtParams.order);
            treeCtrRoot = AddSumtTree (treeCtrRoot, sumtParams.order);
            }
        else
            {
            i = sumtParams.numFileTreesSampled[sumtParams.runId] - 1;
            if (StoreSumtTree (packedTreeList[sumtParams.runId], i, t) == ERROR)
                return (ERROR);
            }

        /* Display the tree nodes. */
#       if 0
        ShowPolyNodes(t);
#       endif
        }
    
    return (NO_ERROR);  
}


int ExamineSumtFile (char *fileName, SumtFileInfo *sumtFileInfo, char *treeName, int *brlensDef)
{
    int     i, foundBegin, lineTerm, inTreeBlock, blockErrors, inSumtComment, lineNum, numTreesInBlock,
            tokenType;
    char    sumtToken[100], *s, *sumtTokenP;
    FILE    *fp;

    /* open binary file */
    if ((fp = OpenBinaryFileR(fileName)) == NULL)
        return ERROR;
        
    /* find out what type of line termination is used for file 1 */
    lineTerm = LineTermType (fp);
    if (lineTerm != LINETERM_MAC && lineTerm != LINETERM_DOS && lineTerm != LINETERM_UNIX)
        {
        MrBayesPrint ("%s   Unknown line termination for file  \"%s\"\n", spacer, fileName);
        return ERROR;
        }

    /* find length of longest line in either file */
    sumtFileInfo->longestLineLength = LongestLine (fp);
    sumtFileInfo->longestLineLength += 10;
    
    /* allocate a string long enough to hold a line */
    s = (char *)SafeMalloc((size_t)(sumtFileInfo->longestLineLength) * sizeof(char));
    if (!s)
        {
        MrBayesPrint ("%s   Problem allocating string for examining file \"%s\"\n", spacer, fileName);
        return (ERROR);
        }
        
    /* close binary file */
    SafeFclose (&fp);
    
    foundBegin = inTreeBlock = blockErrors = inSumtComment = NO;
    lineNum = numTreesInBlock = 0;
    sumtFileInfo->numTreeBlocks = 0;
    sumtFileInfo->lastTreeBlockBegin = 0;
    sumtFileInfo->lastTreeBlockEnd = 0;
    sumtFileInfo->numTreesInLastBlock = 0;

    /* open text file */
    if ((fp = OpenTextFileR(fileName))==NULL)
        {
        MrBayesPrint ("%s   Could not read file \"%s\" in text mode \n", spacer, fileName);
        return (ERROR);
        }

    /* read file */
    while (fgets (s, sumtFileInfo->longestLineLength-2, fp) != NULL)
        {
        sumtTokenP = &s[0];
        do
            {
            if (GetToken (sumtToken, &tokenType, &sumtTokenP))
                goto errorExit;
            if (IsSame("[", sumtToken) == SAME)
                inSumtComment = YES;
            if (IsSame("]", sumtToken) == SAME)
                inSumtComment = NO;
            
            if (inSumtComment == YES)
                {
                if (IsSame ("Param", sumtToken) == SAME)
                    {
                    /* extract the tree name */
                    if (GetToken (sumtToken, &tokenType, &sumtTokenP))   /* get the colon */
                        goto errorExit;
                    if (GetToken (sumtToken, &tokenType, &sumtTokenP))   /* get the tree name */
                        goto errorExit;
                    strcpy (treeName, sumtToken);
                    if (GetToken (sumtToken, &tokenType, &sumtTokenP))
                        goto errorExit;
                    while (IsSame("]", sumtToken) != SAME)
                        {
                        strcat (treeName, sumtToken);
                        if (GetToken (sumtToken, &tokenType, &sumtTokenP))
                            goto errorExit;
                        }
                    inSumtComment = NO;
                    }
                }
            else /* if (inSumtComment == NO) */
                {
                if (foundBegin == YES)
                    {
                    if (IsSame("Trees", sumtToken) == SAME)
                        {
                        numTreesInBlock = 0;
                        inTreeBlock = YES;
                        foundBegin = NO;
                        sumtFileInfo->lastTreeBlockBegin = lineNum;
                        }
                    }
                else
                    {
                    if (IsSame("Begin", sumtToken) == SAME)
                        {
                        if (foundBegin == YES)
                            {
                            MrBayesPrint ("%s   Found inappropriate \"Begin\" statement in file\n", spacer);
                            blockErrors = YES;
                            }
                        foundBegin = YES;
                        }
                    else if (IsSame("End", sumtToken) == SAME)
                        {
                        if (inTreeBlock == YES)
                            {
                            sumtFileInfo->numTreeBlocks++;
                            inTreeBlock = NO;
                            sumtFileInfo->lastTreeBlockEnd = lineNum;
                            }
                        else
                            {
                            MrBayesPrint ("%s   Found inappropriate \"End\" statement in file\n", spacer);
                            blockErrors = YES;
                            }
                        sumtFileInfo->numTreesInLastBlock = numTreesInBlock;
                        }
                    else if (IsSame("Tree", sumtToken) == SAME)
                        {
                        if (inTreeBlock == YES)
                            {
                            numTreesInBlock++;
                            if (numTreesInBlock == 1)
                                {
                                *brlensDef = NO;
                                for (i=0; s[i]!='\0'; i++)
                                    {
                                    if (s[i] == ':')
                                        {
                                        *brlensDef = YES;
                                        break;
                                        }
                                    }
                                }
                            }
                        else
                            {
                            MrBayesPrint ("%s   Found a \"Tree\" statement that is not in a tree block\n", spacer);
                            blockErrors = YES;
                            }
                        }
                    }
                }
                
            } while (*sumtToken);
        lineNum++;
        }

    /* Now, check some aspects of the tree file, such as the number of tree blocks and whether they are properly terminated. */
    if (inTreeBlock == YES)
        {
        MrBayesPrint ("%s   Unterminated tree block in file %s. You probably need to\n", spacer, fileName);
        MrBayesPrint ("%s   add a new line to the end of the file with \"End;\" on it.\n", spacer);
        goto errorExit;
        }
    if (inSumtComment == YES)
        {
        MrBayesPrint ("%s   Unterminated comment in file %s\n", spacer, fileName);
        goto errorExit;
        }
    if (blockErrors == YES)
        {
        MrBayesPrint ("%s   Found formatting errors in file %s\n", spacer, fileName);
        goto errorExit;
        }
    if (sumtFileInfo->lastTreeBlockEnd < sumtFileInfo->lastTreeBlockBegin)
        {
        MrBayesPrint ("%s   Problem reading tree file %s\n", spacer, fileName);
        goto errorExit;
        }
    if (sumtFileInfo->numTreesInLastBlock <= 0)
        {
        MrBayesPrint ("%s   No trees were found in last tree block of file %s\n", spacer, fileName);
        goto errorExit;
        }
    free (s);
    return (NO_ERROR);

errorExit:
    free (s);
    return (ERROR);
}


/* FreePartCtr: Recursively free partition counter nodes */
void FreePartCtr (PartCtr *r)
{
    int     i, j;

    if (r==NULL)
        return;
    
    FreePartCtr (r->left);
    FreePartCtr (r->right);

    /* free relaxed clock parameters: eRate, nEvents, bRate */
    if (sumtParams.nESets > 0)
        {
        for (i=0; i<sumtParams.nESets; i++)
            {
            for (j=0; j<sumtParams.numRuns; j++)
                free (r->nEvents[i][j]);
            free (r->nEvents[i]);
            }
        free (r->nEvents);
        }
    if (sumtParams.nBSets > 0)
        {
        for (i=0; i<sumtParams.nBSets; i++)
            {
            for (j=0; j<sumtParams.numRuns; j++)
                {
                free (r->bLen [i][j]);
                free (r->bRate[i][j]);
                }
            free (r->bLen [i]);
            free (r->bRate[i]);
            }
        free (r->bLen);
        free (r->bRate);
        }

    /* free basic parameters */
    for (i=0; i<sumtParams.numRuns; i++)
        free (r->length[i]);

    free (r->length);
    free (r->count);
    free (r->partition);
    free (r);
    numUniqueSplitsFound--;
    r = NULL;
}


/* FreeSumtParams: Free parameters allocated in sumtParams struct */
void FreeSumtParams(void)
{
    int     i;

    if (memAllocs[ALLOC_SUMTPARAMS] == YES)
        {
        for (i=0; i<sumtParams.numTaxa; i++)
            free(sumtParams.taxaNames[i]);
        free (sumtParams.taxaNames);
        sumtParams.taxaNames = NULL;
        if (sumtParams.numFileTrees) free (sumtParams.numFileTrees);
        sumtParams.numFileTrees = NULL;
        FreePolyTree (sumtParams.tree);
        sumtParams.tree = NULL;
        if (sumtParams.nBSets > 0)
            {
            for (i=0; i<sumtParams.nBSets; i++)
                free(sumtParams.bSetName[i]);
            free (sumtParams.bSetName);
            sumtParams.bSetName = NULL;
            sumtParams.nBSets = 0;
            }
        if (sumtParams.nESets > 0)
            {
            for (i=0; i<sumtParams.nESets; i++)
                free(sumtParams.eSetName[i]);
            free (sumtParams.eSetName);
            sumtParams.eSetName = NULL;
            sumtParams.nESets = 0;
            }
        if (sumtParams.popSizeSet == YES)
            {
            free (sumtParams.popSizeSetName);
            sumtParams.popSizeSetName = NULL;
            sumtParams.popSizeSet = NO;
            }
        memAllocs[ALLOC_SUMTPARAMS] = NO;
        }
}


/* FreeTreeCtr: Recursively free tree counter nodes */
void FreeTreeCtr (TreeCtr *r)
{
    if (r==NULL)
        return;
    
    FreeTreeCtr (r->left);
    FreeTreeCtr (r->right);

    free (r->order);
    free (r);
    numUniqueTreesFound--;
    r = NULL;
}


/* Label: Calculate length of label and fill in char *label if not NULL */
int Label (PolyNode *p, int addIndex, char *label, int maxLength)
{
    int     i, j0, j1, k, n, length, nameLength, index;

    if (p == NULL)
        return 0;

    /* first calculate length */
    if (inSumtCommand == YES && isTranslateDiff == NO)
        {
        for (index=i=0; index<numTaxa; index++)
            {
            if (sumtParams.absentTaxa[index] == YES || taxaInfo[index].isDeleted == YES)
                continue;
            if (p->index == i)
                break;
            else
                i++;
            }
        }
    else
        index = p->index;

    if (addIndex != NO)
        length = (int)(strlen(p->label)) + 4 + (int)(log10(index+1));
    else
        length = (int)(strlen(p->label));
    length = (length > maxLength ? maxLength : length);

    /* fill in label if label != NULL */
    if (label != NULL)
        {
        if (addIndex != NO)
            nameLength = length - 4 - (int)(log10(index+1));
        else
            nameLength = length;

        for (i=0; i<nameLength-1; i++)
            label[i] = p->label[i];
        if ((int)strlen(p->label) > nameLength)
            label[i] = '~';
        else
            label[i] = p->label[i];

        if (addIndex != NO)
            {
            label[++i] = ' ';
            label[++i] = '(';
            n = index + 1;
            k = (int)(log10(n)) + 1;
            while (n != 0)
                {
                j0 = (int)(log10(n));
                j1 = (int)(pow(10,j0));
                label[++i] = '0' + n/j1;
                n = n % j1;
                k--;
                }
            while (k!=0)
                {
                label[++i] = '0';
                k--;
                }
            label[++i] = ')';
            }
        label[++i] = '\0';
        }

    return length;
}


int OpenComptFiles (void)
{
    int         len, previousFiles, oldNoWarn, oldAutoOverwrite;
    char        pFilename[120], dFilename[120];
    FILE        *fpTemp;

    oldNoWarn = noWarn;
    oldAutoOverwrite = autoOverwrite;

    /* set file names */
    strcpy (pFilename, comptreeParams.comptOutfile);
    strcpy (dFilename, comptreeParams.comptOutfile);
    strcat (pFilename, ".pairs");
    strcat (dFilename, ".dists");

    /* one overwrite check for both files */
    previousFiles = NO;
    if (noWarn == NO)
        {
        if ((fpTemp = OpenTextFileR(pFilename)) != NULL)
            {
            previousFiles = YES;
            fclose(fpTemp);
            }
        if ((fpTemp = OpenTextFileR(dFilename)) != NULL)
            {
            previousFiles = YES;
            fclose(fpTemp);
            }
        if (previousFiles == YES)
            {
            MrBayesPrint("%s   There are previous compare results saved using the same filenames.\n", spacer);
            if (WantTo("Do you want to overwrite these results") == YES)
                {
                MrBayesPrint("\n");
                noWarn = YES;
                autoOverwrite = YES;
                }
            else
                {
                MrBayesPrint("\n");
                MrBayesPrint("%s   Please specify a different output file name before running the comparetree command.\n", spacer);
                MrBayesPrint("%s      You can do that using 'comparetree outputfile=<name>'. You can also move or\n", spacer);
                MrBayesPrint("%s      rename the old result files.\n", spacer);
                return ERROR;
                }
            }
        }

    if ((fpParts = OpenNewMBPrintFile (pFilename)) == NULL)
        {
        noWarn = oldNoWarn;
        autoOverwrite = oldAutoOverwrite;
        return ERROR;
        }
    if ((fpDists = OpenNewMBPrintFile (dFilename)) == NULL)
        {
        noWarn = oldNoWarn;
        autoOverwrite = oldAutoOverwrite;
        return ERROR;
        }
        
    /* Reset file flags */
    noWarn = oldNoWarn;
    autoOverwrite = oldAutoOverwrite;

    /* print unique identifiers to each file */
    len = (int) strlen (stamp);
    if (len > 1)
        {
        MrBayesPrintf (fpParts, "[ID: %s]\n", stamp);
        MrBayesPrintf (fpDists, "[ID: %s]\n", stamp);
        }

    return (NO_ERROR);
}


int OpenSumtFiles (int treeNo)
{
    int         i, len,  oldNoWarn, oldAutoOverwrite, previousFiles;
    char        pFilename[120], sFilename[120], vFilename[120], cFilename[120], tFilename[120];
    FILE        *fpTemp;

    oldNoWarn = noWarn;
    oldAutoOverwrite = autoOverwrite;

    /* one overwrite check for all files */
    if (noWarn == NO && treeNo == 0)
        {
        previousFiles = NO;
        for (i=0; i<sumtParams.numTrees; i++)
            {
            if (sumtParams.numTrees > 1)
                {
                sprintf (pFilename, "%s.tree%d.parts", sumtParams.sumtOutfile, i+1);
                sprintf (sFilename, "%s.tree%d.tstat", sumtParams.sumtOutfile, i+1);
                sprintf (vFilename, "%s.tree%d.vstat", sumtParams.sumtOutfile, i+1);
                sprintf (cFilename, "%s.tree%d.con.tre", sumtParams.sumtOutfile, i+1);
                sprintf (tFilename, "%s.tree%d.trprobs", sumtParams.sumtOutfile, i+1);
                }
            else
                {
                sprintf (pFilename, "%s.parts", sumtParams.sumtOutfile);
                sprintf (sFilename, "%s.tstat", sumtParams.sumtOutfile);
                sprintf (vFilename, "%s.vstat", sumtParams.sumtOutfile);
                sprintf (cFilename, "%s.con.tre", sumtParams.sumtOutfile);
                sprintf (tFilename, "%s.trprobs", sumtParams.sumtOutfile);
                }
            if ((fpTemp = TestOpenTextFileR(pFilename)) != NULL)
                {
                previousFiles = YES;
                fclose(fpTemp);
                }
            if ((fpTemp = TestOpenTextFileR(sFilename)) != NULL)
                {
                previousFiles = YES;
                fclose(fpTemp);
                }
            if ((fpTemp = TestOpenTextFileR(vFilename)) != NULL)
                {
                previousFiles = YES;
                fclose(fpTemp);
                }
            if ((fpTemp = TestOpenTextFileR(cFilename)) != NULL)
                {
                previousFiles = YES;
                fclose(fpTemp);
                }
            if ((fpTemp = TestOpenTextFileR(tFilename)) != NULL)
                {
                previousFiles = YES;
                fclose(fpTemp);
                }
            if (previousFiles == YES)
                {
                MrBayesPrint("\n");
                MrBayesPrint("%s   There are previous tree sample summaries saved using the same filenames.\n", spacer);
                if (WantTo("Do you want to overwrite these results") == YES)
                    {
                    MrBayesPrint("\n");
                    noWarn = YES;
                    autoOverwrite = YES;
                    }
                else
                    {
                    MrBayesPrint("\n");
                    MrBayesPrint("%s   Please specify a different output file name before running the sumt command.\n", spacer);
                    MrBayesPrint("%s      You can do that using 'sumt outputfile=<name>'. You can also move or\n", spacer);
                    MrBayesPrint("%s      rename the old result files.\n", spacer);
                    return ERROR;
                    }
                }
            }
        }

    /* set file names */
    if (sumtParams.numTrees > 1)
        {
        sprintf (pFilename, "%s.tree%d.parts", sumtParams.sumtOutfile, treeNo+1);
        sprintf (sFilename, "%s.tree%d.tstat", sumtParams.sumtOutfile, treeNo+1);
        sprintf (vFilename, "%s.tree%d.vstat", sumtParams.sumtOutfile, treeNo+1);
        sprintf (cFilename, "%s.tree%d.con.tre", sumtParams.sumtOutfile, treeNo+1);
        sprintf (tFilename, "%s.tree%d.trprobs", sumtParams.sumtOutfile, treeNo+1);
        }
    else
        {
        sprintf (pFilename, "%s.parts", sumtParams.sumtOutfile);
        sprintf (sFilename, "%s.tstat", sumtParams.sumtOutfile);
        sprintf (vFilename, "%s.vstat", sumtParams.sumtOutfile);
        sprintf (cFilename, "%s.con.tre", sumtParams.sumtOutfile);
        sprintf (tFilename, "%s.trprobs", sumtParams.sumtOutfile);
        }

    /* open files checking for over-write as appropriate */
    if ((fpParts = OpenNewMBPrintFile(pFilename)) == NULL)
        return ERROR;
    if ((fpTstat = OpenNewMBPrintFile(sFilename)) == NULL)
        {
        SafeFclose (&fpParts);
        return ERROR;
        }
    if ((fpVstat = OpenNewMBPrintFile(vFilename)) == NULL)
        {
        SafeFclose (&fpParts);
        SafeFclose (&fpTstat);
        return ERROR;
        }
    if ((fpCon = OpenNewMBPrintFile(cFilename)) == NULL)
        {
        SafeFclose (&fpParts);
        SafeFclose (&fpTstat);
        SafeFclose (&fpVstat);
        return ERROR;
        }
    if (sumtParams.calcTreeprobs == YES)
        {
        if ((fpTrees = OpenNewMBPrintFile(tFilename)) == NULL)
            {
            SafeFclose (&fpParts);
            SafeFclose (&fpTstat);
            SafeFclose (&fpVstat);
            SafeFclose (&fpCon);
            return ERROR;
            }
        }

    /* print #NEXUS if appropriate */
    MrBayesPrintf (fpCon,   "#NEXUS\n");
    if (sumtParams.calcTreeprobs == YES)
        MrBayesPrintf (fpTrees, "#NEXUS\n");

    /* print unique identifiers to each file */
    len = (int) strlen (stamp);
    if (len > 1)
        {
        MrBayesPrintf (fpParts, "[ID: %s]\n", stamp);
        MrBayesPrintf (fpTstat, "[ID: %s]\n", stamp);
        MrBayesPrintf (fpVstat, "[ID: %s]\n", stamp);
        MrBayesPrintf (fpCon,   "[ID: %s]\n", stamp);
        if (sumtParams.calcTreeprobs == YES)
            MrBayesPrintf (fpTrees, "[ID: %s]\n", stamp);
        }

    /* Reset noWarn and autoOverwrite */
    if (treeNo == sumtParams.numTrees - 1)
        {
        noWarn = oldNoWarn;
        autoOverwrite = oldAutoOverwrite;
        }

    return (NO_ERROR);      
}


void PartCtrUppass (PartCtr *r, PartCtr **uppass, int *index)
{
    if (r != NULL)
        {
        uppass[(*index)++] = r;

        PartCtrUppass (r->left, uppass, index);
        PartCtrUppass (r->right, uppass, index);
        }
}


/* PrintBrParamsToFile: Print branach parameters to file */
int PrintBrParamsToFile (PartCtr **treeParts, int numTreeParts, int treeNo, char *divString)
{
    int     i, j, numPartitions, min, treeSample, runNo;
    char    filename[100];
    PartCtr *x;
    FILE    *fp;
    
    /* set file name */
    if (sumtParams.numTrees > 1)
        sprintf (filename, "%s.tree%d.brparams", sumtParams.sumtOutfile, treeNo+1);
    else
        sprintf (filename, "%s.brparams", sumtParams.sumtOutfile);
    MrBayesPrint( "%s   Saving branch parameters to file \"%s\"\n", spacer, filename);
    
    /* Open file checking for over-write as appropriate */
    if ((fp = OpenNewMBPrintFile(filename)) == NULL)
    {
        MrBayesPrint ("\n");
        return ERROR;
    }

    /* count number of branch params to print */
    min = (sumtParams.minBrParamFreq * (sumtParams.numTreesSampled/sumtParams.numRuns));
    for (i=0; i<numTreeParts; i++)
    {
        if (treeParts[i]->totCount < min)
            break;
    }
    numPartitions = i;
    
    /* print header */

    /* -- header for branch lengths; NB! We skip the root branch length (i starts at 1) */
    for (i=1; i<numPartitions; ++i)
    {
        MrBayesPrintf (fp, "length%s[%d]", divString, i);
        if (i!=numPartitions-1)
            MrBayesPrintf (fp, "\t");
    }
    
    /* -- header for node heights */
    if (sumtParams.isClock == YES )
        {
        MrBayesPrintf (fp, "\t");
        for (i=0; i<numPartitions; i++)
            {
            MrBayesPrintf (fp, "height%s[%d]", divString, i);
            if (i!=numPartitions-1)
                MrBayesPrintf (fp, "\t");
            }
        }
    
    /* -- header for node ages */
    if (sumtParams.isCalibrated == YES )
    {
        MrBayesPrintf (fp, "\t");
        for (i=0; i<numPartitions; i++)
        {
            MrBayesPrintf (fp, "age%s[%d]", divString, i);
            if (i!=numPartitions-1)
                MrBayesPrintf (fp, "\t");
        }
    }

    /* -- header for effective branch lengths and branch rates; NB! We skip the root branch (j starts at 1) */
    if (sumtParams.isRelaxed == YES)
    {
        for (i=0; i<sumtParams.nBSets; i++)
        {
            MrBayesPrintf (fp, "\t");
            for (j=1; j<numPartitions; j++)
            {
                MrBayesPrintf (fp, "%s_length[%d]", sumtParams.bSetName[i], j);
                MrBayesPrintf (fp, "\t");
            }

            for (j=1; j<numPartitions; j++)
            {
                MrBayesPrintf (fp, "%s_rate[%d]", sumtParams.bSetName[i], j);
                if (j!=numPartitions-1)
                    MrBayesPrintf (fp, "\t");
            }
        }
        for (i=0; i<sumtParams.nESets; i++)
        {
            MrBayesPrintf (fp, "\t");
            for (j=1; j<numPartitions; j++)
            {
                MrBayesPrintf (fp, "%s_nEvents[%d]", sumtParams.eSetName[i], j);
                if (j!=numPartitions-1)
                    MrBayesPrintf (fp, "\t");
            }
        }
    }

    /* finalize header line */
    MrBayesPrintf (fp, "\n");
    
    /* print values */
    for (runNo=0; runNo<sumtParams.numRuns; runNo++)
    {
        for (treeSample=0; treeSample<sumtParams.numTreesSampled/sumtParams.numRuns; treeSample++)
        {
            /* print branch lengths */
            for (i=1; i<numPartitions; i++)
            {
                x = treeParts[i];
                if ( x->count[runNo] > treeSample )
                    MrBayesPrintf (fp, "%s", MbPrintNum (x->length[runNo][treeSample]));
                else
                    MrBayesPrintf (fp, "NA");
                if (i!=numPartitions-1)
                    MrBayesPrintf (fp, "\t");
            }
            /* print node heights */
            if (sumtParams.isClock == YES)
            {
                MrBayesPrintf (fp, "\t");
                for (i=0; i<numPartitions; i++)
                {
                    x = treeParts[i];
                    if ( x->count[runNo] > treeSample )
                        MrBayesPrintf (fp, "%s", MbPrintNum (x->height[runNo][treeSample]));
                    else
                        MrBayesPrintf (fp, "NA");
                    if (i!=numPartitions-1)
                        MrBayesPrintf (fp, "\t");
                }
            }
            /* print node ages */
            if (sumtParams.isCalibrated == YES)
            {
                MrBayesPrintf (fp, "\t");
                for (i=0; i<numPartitions; i++)
                {
                    x = treeParts[i];
                    if ( x->count[runNo] > treeSample )
                        MrBayesPrintf (fp, "%s", MbPrintNum (x->age[runNo][treeSample]));
                    else
                        MrBayesPrintf (fp, "NA");
                    if (i!=numPartitions-1)
                        MrBayesPrintf (fp, "\t");
                }
            }
            /* print effective branch lengths */
            if (sumtParams.isRelaxed)
            {
                for (i=0; i<sumtParams.nBSets; i++)
                {
                    MrBayesPrintf (fp, "\t");
                    for (j=1; j<numPartitions; j++)
                    {
                        x = treeParts[j];
                        if ( x->count[runNo] > treeSample )
                            MrBayesPrintf (fp, "%s", MbPrintNum (x->bLen[i][runNo][treeSample]));
                        else
                            MrBayesPrintf (fp, "NA");
                        MrBayesPrintf (fp, "\t");
                    }
                    for (j=1; j<numPartitions; j++)
                    {
                        x = treeParts[j];
                        if ( x->count[runNo] > treeSample )
                            MrBayesPrintf (fp, "%s", MbPrintNum (x->bRate[i][runNo][treeSample]));
                        else
                            MrBayesPrintf (fp, "NA");
                        if (j!=numPartitions-1)
                            MrBayesPrintf (fp, "\t");
                    }
                }
                for (i=0; i<sumtParams.nESets; i++)
                {
                    MrBayesPrintf (fp, "\t");
                    for (j=1; j<numPartitions; j++)
                    {
                        x = treeParts[j];
                        if ( x->count[runNo] > treeSample )
                            MrBayesPrintf (fp, "%d", x->nEvents[i][runNo][treeSample]);
                        else
                            MrBayesPrintf (fp, "NA");
                        if (j!=numPartitions-1)
                            MrBayesPrintf (fp, "\t");
                    }
                }
            }
            MrBayesPrintf (fp, "\n");
        }   // next tree sample
    }   // next run

    MrBayesPrint("\n");
    return NO_ERROR;
}


/* PrintConTree: Print consensus tree in standard format readable by TreeView, Paup etc */
void PrintConTree (FILE *fp, PolyTree *t)
{
    MrBayesPrintf (fp, "   [Note: This tree contains information on the topology, \n");
    MrBayesPrintf (fp, "          branch lengths (if present), and the probability\n");
    MrBayesPrintf (fp, "          of the partition indicated by the branch.]\n");
    if (!strcmp(sumtParams.sumtConType, "Halfcompat"))
        MrBayesPrintf (fp, "   tree con_50_majrule = ");
    else
        MrBayesPrintf (fp, "   tree con_all_compat = ");
    WriteConTree (t->root, fp, YES);
    MrBayesPrintf (fp, ";\n");
    if (sumtParams.brlensDef == YES)
        {
        MrBayesPrintf (fp, "\n");
        MrBayesPrintf (fp, "   [Note: This tree contains information only on the topology\n");
        MrBayesPrintf (fp, "          and branch lengths (median of the posterior probability density).]\n");
        if (!strcmp(sumtParams.sumtConType, "Halfcompat"))
            MrBayesPrintf (fp, "   tree con_50_majrule = ");
        else
            MrBayesPrintf (fp, "   tree con_all_compat = ");
        WriteConTree (t->root, fp, NO);
        MrBayesPrintf (fp, ";\n");
        }
}


/* PrintFigTreeConTree: Print consensus tree in rich format for FigTree */
void PrintFigTreeConTree (FILE *fp, PolyTree *t, PartCtr **treeParts)
{
    if (!strcmp(sumtParams.sumtConType, "Halfcompat"))
        MrBayesPrintf (fp, "   tree con_50_majrule = ");
    else
        MrBayesPrintf (fp, "   tree con_all_compat = ");
    if (t->isRooted == YES)
        MrBayesPrintf (fp, "[&R] ");
    else
        MrBayesPrintf (fp, "[&U] ");

    WriteFigTreeConTree (t->root, fp, treeParts);
    MrBayesPrintf (fp, ";\n");
}


void PrintFigTreeNodeInfo (FILE *fp, PartCtr *x, MrBFlt length)
{
    int     i, postProbPercent, postProbSdPercent;
    MrBFlt  *support, mean, var, min, max;
    Stat    theStats;

    support = SafeCalloc (sumtParams.numRuns, sizeof(MrBFlt));
    for (i=0; i<sumtParams.numRuns; i++)
        {
        support[i] = (MrBFlt) x->count[i] / (MrBFlt) sumtParams.numFileTreesSampled[i];
        }
    if (sumtParams.numRuns > 1)
        {
        MeanVariance (support, sumtParams.numRuns, &mean, &var);
        Range (support, sumtParams.numRuns, &min, &max);
        postProbPercent = (int) (100.0*mean + 0.5);
        postProbSdPercent = (int) (100.0 * sqrt(var) + 0.5);
        fprintf (fp, "[&prob=%.8le,prob_stddev=%.8le,prob_range={%.8le,%.8le},prob(percent)=\"%d\",prob+-sd=\"%d+-%d\"",
            mean, sqrt(var), min, max, postProbPercent, postProbPercent, postProbSdPercent);
        }
    else
        {
        postProbPercent = (int) (100.0*support[0] + 0.5);
        fprintf (fp, "[&prob=%.8le,prob(percent)=\"%d\"", support[0], postProbPercent);
        }
    if (sumtParams.isClock == YES)
        {
        GetSummary (x->height, sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);
        if (sumtParams.HPD == YES)
            fprintf (fp, ",height_mean=%.8le,height_median=%.8le,height_95%%HPD={%.8le,%.8le}", theStats.mean, theStats.median, theStats.lower, theStats.upper);
        else
            fprintf (fp, ",height_mean=%.8le,height_median=%.8le,height_95%%CredInt={%.8le,%.8le}", theStats.mean, theStats.median, theStats.lower, theStats.upper);
        }
    if (sumtParams.isCalibrated == YES)
        {
        GetSummary (x->age, sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);
        if (sumtParams.HPD == YES)
            fprintf (fp, ",age_mean=%.8le,age_median=%.8le,age_95%%HPD={%.8le,%.8le}", theStats.mean, theStats.median, theStats.lower, theStats.upper);
        else
            fprintf (fp, ",age_mean=%.8le,age_median=%.8le,age_95%%CredInt={%.8le,%.8le}", theStats.mean, theStats.median, theStats.lower, theStats.upper);
        }
    fprintf (fp, "]");
    if (length >= 0.0)
        fprintf (fp, ":%s", MbPrintNum(length));
    if (sumtParams.brlensDef == YES)
        {
        GetSummary (x->length, sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);
        if (sumtParams.HPD == YES)
            fprintf (fp, "[&length_mean=%.8le,length_median=%.8le,length_95%%HPD={%.8le,%.8le}", theStats.mean, theStats.median, theStats.lower, theStats.upper);
        else
            fprintf (fp, "[&length_mean=%.8le,length_median=%.8le,length_95%%CredInt={%.8le,%.8le}", theStats.mean, theStats.median, theStats.lower, theStats.upper);
        }
    if (sumtParams.isClock == YES && sumtParams.isRelaxed == YES)
        {
        for (i=0; i<sumtParams.nBSets; i++)
            {
            GetSummary (x->bLen[i], sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);
            if (sumtParams.HPD == YES)
                fprintf (fp, ",effectivebrlen%s_mean=%lf,effectivebrlen%s_median=%lf,effectivebrlen%s_95%%HPD={%lf,%lf}",
                    sumtParams.tree->bSetName[i], theStats.mean,
                    sumtParams.tree->bSetName[i], theStats.median,
                    sumtParams.tree->bSetName[i], theStats.lower,
                    theStats.upper);
            else
                fprintf (fp, ",effectivebrlen%s_mean=%lf,effectivebrlen%s_median=%lf,effectivebrlen%s_95%%CredInt={%lf,%lf}",
                    sumtParams.tree->bSetName[i], theStats.mean,
                    sumtParams.tree->bSetName[i], theStats.median,
                    sumtParams.tree->bSetName[i], theStats.lower,
                    theStats.upper);
            GetSummary (x->bRate[i], sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);
            if (sumtParams.HPD == YES)
                fprintf (fp, ",rate%s_mean=%lf,rate%s_median=%lf,rate%s_95%%HPD={%lf,%lf}",
                    sumtParams.tree->bSetName[i], theStats.mean,
                    sumtParams.tree->bSetName[i], theStats.median,
                    sumtParams.tree->bSetName[i], theStats.lower,
                    theStats.upper);
            else
                fprintf (fp, ",rate%s_mean=%lf,rate%s_median=%lf,rate%s_95%%CredInt={%lf,%lf}",
                    sumtParams.tree->bSetName[i], theStats.mean,
                    sumtParams.tree->bSetName[i], theStats.median,
                    sumtParams.tree->bSetName[i], theStats.lower,
                    theStats.upper);
            }
        for (i=0; i<sumtParams.nESets; i++)
            {
            GetIntSummary (x->nEvents[i], sumtParams.numRuns, x->count, &theStats, sumtParams.HPD);
            if (sumtParams.HPD == YES)
                fprintf (fp, ",nEvents%s_mean=%lf,nEvents%s_median=%lf,nEvents%s_95%%HPD={%lf,%lf}",
                    sumtParams.tree->eSetName[i], theStats.mean,
                    sumtParams.tree->eSetName[i], theStats.median,
                    sumtParams.tree->eSetName[i], theStats.lower,
                    theStats.upper);
            else
                fprintf (fp, ",nEvents%s_mean=%lf,nEvents%s_median=%lf,nEvents%s_95%%CredInt={%lf,%lf}",
                    sumtParams.tree->eSetName[i], theStats.mean,
                    sumtParams.tree->eSetName[i], theStats.median,
                    sumtParams.tree->eSetName[i], theStats.lower,
                    theStats.upper);
            }
        }
    if (sumtParams.brlensDef == YES)
        fprintf (fp, "]");

    free (support);
}


void PrintSumtTableLine(int numRuns, int *rowCount, Stat *theStats, MrBFlt *numPSRFSamples, MrBFlt *maxPSRF, MrBFlt *sumPSRF)
{
    int j,k;

    MrBayesPrint ("%10.6lf  %10.6lf  %10.6lf  %10.6lf  %10.6lf", theStats->mean, theStats->var, theStats->lower, theStats->upper, theStats->median);

    MrBayesPrintf (fpVstat, "\t%s", MbPrintNum(theStats->mean));
    MrBayesPrintf (fpVstat, "\t%s", MbPrintNum(theStats->var));
    MrBayesPrintf (fpVstat, "\t%s", MbPrintNum(theStats->lower));
    MrBayesPrintf (fpVstat, "\t%s", MbPrintNum(theStats->upper));
    MrBayesPrintf (fpVstat, "\t%s", MbPrintNum(theStats->median));

    if (numRuns > 1)
        {
        for (j=k=0; j<numRuns; j++)
            if (rowCount[j] > 0)
                k++;
            if (theStats->PSRF < 0.0)
                {
                MrBayesPrint ("     NA    %3d", k);
                MrBayesPrintf (fpVstat, "\tNA\t%d", k);
                }
            else
                {
                if (theStats->PSRF > 10.0)
                    {
                    MrBayesPrint ("    >10.0  %3d", k);
                    MrBayesPrintf (fpVstat, "\tNA\t%d", k);
                    (*maxPSRF) = 10.0;
                    }
                else
                    {
                    MrBayesPrint ("  %7.3lf  %3d", theStats->PSRF, k);
                    MrBayesPrintf (fpVstat, "\t%s\t%d", MbPrintNum(theStats->PSRF), k);
                    (*sumPSRF) += theStats->PSRF;
                    (*numPSRFSamples)++;
                    if (theStats->PSRF > *maxPSRF)
                         (*maxPSRF) = theStats->PSRF;
                    }
                }

                if (k != numRuns)
                    MrBayesPrint (" *");
        }

        MrBayesPrintf (fpVstat, "\n");
        MrBayesPrint ("\n");
}


/* PrintSumtTaxaInfo: Print information on pruned and absent taxa */
void PrintSumtTaxaInfo (void)
{
    int     i, j, lineWidth, numExcludedTaxa, len;
    char    tempStr[100];

    /* print out information on absent taxa */
    numExcludedTaxa = 0;
    for (i=0; i<numTaxa; i++)
        if (sumtParams.absentTaxa[i] == YES)
            numExcludedTaxa++;

    if (numExcludedTaxa > 0)
        {
        if (numExcludedTaxa == 1)
            MrBayesPrint ("%s   The following taxon was absent from trees:\n", spacer);
        else
            MrBayesPrint ("%s   The following %d taxa were absent from trees:\n", spacer, numExcludedTaxa);
        MrBayesPrint ("%s      ", spacer);
        j = lineWidth = 0;
        for (i=0; i<numTaxa; i++)
            {
            if (sumtParams.absentTaxa[i] == YES)
                {
                j++;
                strcpy (tempStr, taxaNames[i]);
                len = (int) strlen(tempStr);
                lineWidth += len+2;
                if (lineWidth > 60)
                    {
                    MrBayesPrint ("\n%s      ", spacer);
                    lineWidth = 0;
                    }
                if (numExcludedTaxa == 1)
                    MrBayesPrint ("%s\n", tempStr);
                else if (numExcludedTaxa == 2 && j == 1)
                    MrBayesPrint ("%s ", tempStr);
                else if (j == numExcludedTaxa)
                    MrBayesPrint ("and %s\n", tempStr);
                else
                    MrBayesPrint ("%s, ", tempStr);
                }
            }
        MrBayesPrint ("\n");
        }

    /* print out information on pruned taxa */
    numExcludedTaxa = 0;
    for (i=0; i<numTaxa; i++)
        if (taxaInfo[i].isDeleted == YES && sumtParams.absentTaxa[i] == NO)
            numExcludedTaxa++;
    
    if (numExcludedTaxa > 0)
        {
        if (numExcludedTaxa == 1)
            MrBayesPrint ("%s   The following taxon was pruned from trees:\n", spacer);
        else
            MrBayesPrint ("%s   The following %d taxa were pruned from trees:\n", spacer, numExcludedTaxa);
        MrBayesPrint ("%s      ", spacer);
        j = lineWidth = 0;
        for (i=0; i<numTaxa; i++)
            {
            if (taxaInfo[i].isDeleted == YES && sumtParams.absentTaxa[i] == NO)
                {
                j++;
                strcpy (tempStr, taxaNames[i]);
                len = (int) strlen(tempStr);
                lineWidth += len+2;
                if (lineWidth > 60)
                    {
                    MrBayesPrint ("\n%s      ", spacer);
                    lineWidth = 0;
                    }
                if (numExcludedTaxa == 1)
                    MrBayesPrint ("%s\n", tempStr);
                else if (numExcludedTaxa == 2 && j == 1)
                    MrBayesPrint ("%s ", tempStr);
                else if (j == numExcludedTaxa)
                    MrBayesPrint ("and %s\n", tempStr);
                else
                    MrBayesPrint ("%s, ", tempStr);
                }
            }
        MrBayesPrint ("\n");
        }
}


/* Range: Determine range for a vector of MrBFlt values */
void Range (MrBFlt *vals, int nVals, MrBFlt *min, MrBFlt *max)
{    
    SortMrBFlt (vals, 0, nVals-1);
    
    *min  = vals[0];
    *max  = vals[nVals-1];
}


/* ResetTaxonSet: Reset included taxa and local outgroup number */
void ResetTaxonSet (void)
{
    int     i, j;

    /* reset numLocalTaxa and localOutGroup */
    localOutGroup = 0;
    numLocalTaxa = 0;
    for (i=j=0; i<numTaxa; i++)
        {
        if (taxaInfo[i].isDeleted == NO)
            {
            if (i == outGroupNum)
                localOutGroup = numLocalTaxa;
            numLocalTaxa++;
            }
        }
}


void ResetTranslateTable (void)
{
    int i;

    for (i=0; i<numTranslates; i++)
        {
        free (transFrom[i]);
        free (transTo[i]);
        }
    free (transFrom);
    free (transTo);
    transFrom = NULL;
    transTo = NULL;
    numTranslates = 0;
    isTranslateDef = NO;
    isTranslateDiff = NO;
}


int ShowConPhylogram (FILE *fp, PolyTree *t, int screenWidth)
{
    int             i, j, k, nLines, from, to, treeWidth=0, barLength, printExponential,
                    precision, width, newPos, curPos, nTimes, numSpaces, maxLabelLength;
    char            *printLine, *markLine, temp[30], *label;
    MrBFlt          scale, f, scaleBar;
    PolyNode        *p, *q;

    /* set max label length */
    maxLabelLength = 20;

    /* allocate space for label, printLine and markLine */
    printLine = (char *) SafeCalloc ((2*screenWidth+2),sizeof(char)); 
    label = (char *) SafeCalloc (maxLabelLength+1, sizeof(char));
    if (!printLine || !label)
        return ERROR;
    markLine = printLine + screenWidth + 1;

    /* calculate scale */
    scale = 0.0;
    t->root->f = 0.0;
    for (i=t->nNodes-2; i>=0; i--)
        {
        p = t->allDownPass[i];
        /* find distance to root in relevant units */
        if (sumtParams.isClock == YES && sumtParams.isCalibrated == NO)
            p->f = t->root->depth - p->depth;
        else if (sumtParams.isClock == YES && sumtParams.isCalibrated == YES)
            p->f = t->root->age - p->age;
        else
            p->f = p->anc->f + p->length;
        if (p->left == NULL)
            {
            f = p->f / (screenWidth - Label(p,YES,NULL,maxLabelLength) - 2);
            if (f > scale)
                {
                scale = f;
                treeWidth = screenWidth - Label(p,YES,NULL,maxLabelLength) - 2;
                }
            }
        }
    
    /* calculate x coordinates */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->x = (int) (0.5 + (p->f / scale));
        }

    /* calculate y coordinates and lines to print */
    for (i=nLines=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left != NULL)
            {
            /* internal node */
            for (q=p->left->sib; q->sib!=NULL; q=q->sib)
                ;
            p->y = (int) (0.5 + ((p->left->y + q->y) / 2.0));
            }
        else 
            {
            /* terminal node */
            p->y = nLines;
            nLines += 2;
            }
        }

    /* print tree line by line */

    for (i=0; i<nLines; i++)
        {
        MrBayesPrint ("%s   ", spacer);
        /* empty printLine */
        for (j=0; j<screenWidth; j++)
            {
            printLine[j] = ' ';
            }   
        printLine[j]='\0';

        for (j=0; j<t->nNodes; j++)
            {
            p = t->allDownPass[j];
            if (p->y != i)
                continue;

            /* this branch should be printed */
            /* add branch */
            if (p->anc == NULL)
                {
                /* this is the root of the whole tree */
                printLine[p->x] = '+';
                }
            else
                {
                /* this is an ordinary branch */
                to = p->x;
                from = p->anc->x;
                for (k=from+1; k<=to; k++)
                    printLine[k] = '-';
                if (p == p->anc->left)
                    {
                    if (markLine[from] == 0)
                        printLine[from] = '/';
                    else
                        printLine[from] = '|';
                    markLine[from] ++;
                    }
                else if (p->sib == NULL)
                    {
                    if (markLine[from] == 1)
                        printLine[from] = '\\';
                    else
                        printLine[from] = '|';
                    markLine[from] --;
                    }
                if (p->left!=NULL)
                    {
                    if (from != to)
                        printLine[to] = '+';
                    else
                        printLine[to] = '|';
                    }
                else
                    {
                    /* add label if the branch is terminal */
                    Label(p,YES,label,maxLabelLength);
                    sprintf (printLine+to+2,"%s", label);
                    }
                }
            }

        /* check for cross branches */
        for (j=0; j<screenWidth; j++)
            {
            if (markLine[j] >= 1 && printLine[j] == ' ')
                printLine[j] = '|';
            }
        MrBayesPrintf (fp, "%s\n",printLine);
        }

    /* print scale */
    k = (int) (floor (log10 (scale * 80)));
    scaleBar = pow (10, k);
    barLength = (int) (scaleBar / scale);
    if (barLength > 80)
        {
        barLength /= 10;
        scaleBar /= 10.0;
        }
    else if (barLength > 40)
        {
        barLength /= 5;
        scaleBar /= 5.0;
        }
    else if (barLength > 16)
        {
        barLength /= 2;
        scaleBar /= 2.0;
        }

    if (t->isClock == YES)
        {
        MrBayesPrint ("%s   ", spacer);
        for (i=0; i<treeWidth; i++)
            printLine[i] = '-';
        nTimes = (int) (treeWidth / (scaleBar / scale));
        for (i=0; i<=nTimes; i++)
            printLine[treeWidth - (int)(i*scaleBar/scale)] = '|';
        MrBayesPrint ("%s\n", printLine);
        
        MrBayesPrint ("%s   ", spacer);
        f = treeWidth * scale;
        if (f >= 1000.0 || f < 0.10)
            {
            printExponential = YES;
            precision = 0;
            }
        else
            {
            printExponential = NO;
            precision = 2 - (int) (log10 (f));
            }
        
        curPos = 0;
        f = nTimes * scaleBar;
        while (curPos < treeWidth)
            {
            /* print the number */
            if (printExponential == YES)
                sprintf (temp, "%.2e", f);
            else
                sprintf (temp, "%.*lf", precision, f);
            
            /* room to print ? if so, print */
            width = (int) strlen (temp);
            newPos = treeWidth - (int) (nTimes * (scaleBar / scale));
            numSpaces = newPos - width / 2 - curPos;
            if (numSpaces >= 0 || (numSpaces >= -2 && curPos == 0))
                {
                while (numSpaces > 0)
                    {
                    printLine[curPos++] = ' ';
                    numSpaces--;
                    }
                for (i=0; temp[i]!='\0'; i++)
                    printLine[curPos++] = temp[i];
                }

            /* get new number */
            f -= scaleBar;
            nTimes--;
            }

        MrBayesPrint ("%s\n", printLine);
            
        if (sumtParams.isCalibrated == YES)
            MrBayesPrint ("\n%s   [User-defined time units]\n\n", spacer);
        else
            MrBayesPrint ("\n%s   [Expected changes per site]\n\n", spacer);
        }
    else
        {
        MrBayesPrintf (fp, "%s   |", spacer);
        for (i=0; i<barLength-1; i++)
            MrBayesPrintf (fp, "-");
        MrBayesPrintf (fp, "| %1.3lf expected changes per site\n\n", scaleBar);
        }

    free (printLine);

    return NO_ERROR;
}


int ShowConTree (FILE *fp, PolyTree *t, int screenWidth, int showSupport)
{
    int             i, j, k, treeWidth, minBranchLength, maxWidth, isTreeDivided,
                    printWidth, nLines, nodesToBePrinted, from, to, maxLabelLength,
                    maxLength;
    char            *printLine, *markLine, temp[20], *label;
    PolyNode        *p=NULL, *q;

    maxLength = 20;         /* max length of label */
    minBranchLength = 5;    /* min length of branch in tree */
    isTreeDivided = NO;
    
    /* allocate space for printLine, markLine and label */
    printLine = (char *) SafeCalloc (maxLength+1+(2*screenWidth+2),sizeof(char));
    if (!printLine)
        return ERROR;
    markLine = printLine + screenWidth + 1;
    label = markLine + screenWidth + 1;

    /* get fresh internal node indices */
    k = t->nNodes - t->nIntNodes;
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        p->index = k++;
        }
    
    /* calculate max length of labels including taxon index number */
    maxLabelLength = 0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL)
            {
            j = Label(p,YES,NULL,maxLength);
            if (j > maxLabelLength)
                maxLabelLength = j;
            }
        }

    /* make sure label can hold an interior node index number */
    j = (int) (3.0 + log10((MrBFlt)t->nNodes)); 
    maxLabelLength = (maxLabelLength > j ? maxLabelLength : j);

    /* calculate remaining screen width for tree
       and maxWidth in terms of branches */
    treeWidth = screenWidth - maxLabelLength - 1;
    maxWidth = treeWidth / minBranchLength;
    
    /* unmark whole tree */
    for (i=0; i<t->nNodes; i++)
        t->allDownPass[i]->mark = 0;
    nodesToBePrinted = t->nNodes;

    while (nodesToBePrinted > 0)
        {
        /* count depth of nodes in unprinted tree */
        for (i=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            if (p->mark == 0)   /* the node has not been printed yet */
                {
                p->x = 0;
                /* if it is an interior node in the tree that will be printed
                   compute the depth of the node */
                if (p->left != NULL && p->left->mark == 0)
                    {
                    for (q = p->left; q!=NULL; q=q->sib)
                        {
                        if (q->x > p->x)
                            p->x = q->x;
                        }
                    p->x++;
                    /* break when root of print subtree has been found */
                    if (p->x >= maxWidth)
                        break;
                    }
                }
            }

        /* if internal node then find largest nonprinted subtree among descendant nodes */
        if (p->anc != NULL)
            {
            for (q=p->left; q!=NULL; q=q->sib)
                {
                if (q->x == p->x - 1 && q->mark == 0)
                    p = q;
                }
            MrBayesPrintf (fp, "%s   Subtree rooted at node %d:\n\n", spacer, p->index);
            isTreeDivided = YES;
            }
        else if (isTreeDivided == YES)
            MrBayesPrintf (fp, "%s   Root part of tree:\n\n", spacer);

        /* mark subtree for printing and
           translate x coordinates from depth to position */
        if (p->anc == NULL)
            printWidth = p->x;
        else
            printWidth = p->x + 1;
        p->mark = 1;
        p->x = (int) (treeWidth - 0.5 - ((treeWidth - 1) * (p->x / (MrBFlt) printWidth)));
        for (i=t->nNodes-2; i>=0; i--)
            {
            p = t->allDownPass[i];
            if (p->mark == 0 && p->anc->mark == 1)
                {   
                p->mark = 1;
                p->x = (int) (treeWidth - 0.5 - ((treeWidth - 1) * (p->x / (MrBFlt) printWidth)));
                }
            }

        /* calculate y coordinates of nodes to be printed and lines to print */
        for (i=nLines=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            if (p->mark == 1)
                {
                if (p->left != NULL && p->left->mark == 1)
                    {
                    /* internal node */
                    for (q=p->left->sib; q->sib!=NULL; q=q->sib)
                        ;
                    p->y = (int) (0.5 + ((p->left->y + q->y) / 2.0));
                    }
                else 
                    {
                    /* terminal node */
                    p->y = nLines;
                    nLines += 2;
                    }
                }
            }

        /* print subtree line by line */
        for (i=0; i<nLines; i++)
            {
            MrBayesPrintf (fp, "%s   ", spacer);
            /* empty printLine */
            for (j=0; j<screenWidth; j++)
                {
                printLine[j] = ' ';
                }   
            printLine[j]='\0';

            for (j=0; j<t->nNodes; j++)
                {
                p = t->allDownPass[j];
                if (p->mark != 1 || p->y != i)
                    continue;

                /* this branch should be printed
                   add label if the branch is terminal in tree to be printed */
                if (p->left == NULL)
                    {
                    Label (p,YES,label,maxLength);
                    sprintf (printLine+treeWidth+1,"%s", label);
                    }
                else if (p->left->mark == 2)
                    sprintf (printLine+treeWidth+1,"(%d)", p->index);
                
                /* add branch */
                if (p->anc == NULL)
                    {
                    /* this is the root of the whole tree */
                    printLine[p->x] = '+';
                    nodesToBePrinted--;
                    }
                else if (p->anc->mark == 0)
                    {
                    /* this is a root of a subtree
                       this branch will have to be printed again so do
                       not decrease nodesToBePrinted */
                    to = p->x;
                    from = 0;
                    for (k=from; k<to; k++)
                        printLine[k] = '-';
                    printLine[to] = '+';
                    if (showSupport == YES)
                        sprintf (temp, "%d", (int) (p->support*100.0 + 0.5));
                    else
                        *temp='\0';
                    from = (int)(from + 1.5 + ((to - from - 1 - strlen(temp)) / 2.0));
                    for (k=0; temp[k]!='\0'; k++)
                        printLine[from++] = temp[k];
                    }
                else
                    {
                    /* this is an ordinary branch */
                    to = p->x;
                    from = p->anc->x;
                    for (k=from+1; k<=to; k++)
                        printLine[k] = '-';
                    if (p == p->anc->left)
                        {
                        printLine[from] = '/';
                        markLine[from] = 1;
                        }
                    else if (p->sib == NULL)
                        {
                        printLine[from] = '\\';
                        markLine[from] = 0;
                        }
                    if (p->left!=NULL && p->left->mark!=2)
                        {
                        printLine[to] = '+';
                        if (showSupport == YES)
                            sprintf (temp, "%d", (int) (p->support*100.0 + 0.5));
                        else
                            *temp='\0';
                        from = (int)(from + 1.5 + ((to - from - 1 - strlen(temp)) / 2.0));
                        for (k=0; temp[k]!='\0'; k++)
                            printLine[from++] = temp[k];
                        }
                    nodesToBePrinted--;
                    }
                }

            /* check for cross branches */
            for (j=0; j<treeWidth; j++)
                {
                if (markLine[j] == 1 && printLine[j] == ' ')
                    printLine[j] = '|';
                }
        
            MrBayesPrintf (fp, "%s\n",printLine);
            }

        /* mark printed branches */
        for (i=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            if (p->mark == 1)
                {
                if (p->anc == NULL)
                    p->mark = 2;
                else if (p->anc->mark == 0)
                    p->mark = 0;    /* this branch will have to be printed again */
                else
                    p->mark = 2;
                }
            }

        }   /* next subtree */
    
    free (printLine);

    return NO_ERROR;
}


void ShowParts (FILE *fp, BitsLong *p, int nTaxaToShow)
{
    int         i;
    BitsLong    x, y, bitsLongOne;

    bitsLongOne = 1;
    
    for (i=0; i<nTaxaToShow; i++)
        {
        y = p[i / nBitsInALong];
        x = bitsLongOne << (i % nBitsInALong);
        if ((x & y) == 0)
            MrBayesPrintf (fp, ".");
        else
            MrBayesPrintf (fp, "*");
        }
}


void ShowSomeParts (FILE *fp, BitsLong *p, int offset, int nTaxaToShow)
{
    int         i;
    BitsLong    x, y, bitsLongOne;

    bitsLongOne = 1;
    
    for (i=offset; i<offset+nTaxaToShow; i++)
        {
        y = p[i / nBitsInALong];
        x = bitsLongOne << (i % nBitsInALong);
        if ((x & y) == 0)
            MrBayesPrintf (fp, ".");
        else
            MrBayesPrintf (fp, "*");
        }
}


void SortPartCtr (PartCtr **item, int left, int right)
{
    register int    i, j;
    PartCtr         *tempPartCtr;
    int             x;

    assert (left >= 0);
    assert (right >= 0);

    i = left;
    j = right;
    x = item[(left+right)/2]->totCount;
    do 
        {
        while (item[i]->totCount > x && i < right)
            i++;
        while (x > item[j]->totCount && j > left)
            j--;
        if (i <= j)
            {
            tempPartCtr = item[i];
            item[i] = item[j];
            item[j] = tempPartCtr;
                
            i++;
            j--;
            }
        } while (i <= j);
    if (left < j)
        SortPartCtr (item, left, j);
    if (i < right)
        SortPartCtr (item, i, right);
}


void SortTerminalPartCtr (PartCtr **item, int len)
{
    register int    i, j, maxCount;
    PartCtr         *temp;

    maxCount = item[0]->totCount;
    
    /* put root first */
    for (i=0; item[i]->totCount == maxCount; i++)
        if (NumBits(item[i]->partition, sumtParams.BitsLongsNeeded) == sumtParams.numTaxa)
            break;

    if (i!=0)
        {
        temp = item[0];
        item[0] = item[i];
        item[i] = temp;
        }

    /* then find terminals in index order */
    for (i=1; i<=sumtParams.numTaxa; i++)
        {
        for (j=i; item[j]->totCount == maxCount && j<len; j++)
            if (NumBits(item[j]->partition, sumtParams.BitsLongsNeeded) == 1 && 
                FirstTaxonInPartition(item[j]->partition, sumtParams.BitsLongsNeeded) == i-1)
                break;

        if (j!=i)
            {
            temp = item[i];
            item[i] = item[j];
            item[j] = temp;
            }
        }
}


void SortTreeCtr (TreeCtr **item, int left, int right)
{
    register int    i, j;
    TreeCtr         *tempTreeCtr;
    int             x;

    i = left;
    j = right;
    x = item[(left+right)/2]->count;
    do 
        {
        while (item[i]->count > x && i < right)
            i++;
        while (x > item[j]->count && j > left)
            j--;
        if (i <= j)
            {
            tempTreeCtr = item[i];
            item[i] = item[j];
            item[j] = tempTreeCtr;
                
            i++;
            j--;
            }
        } while (i <= j);
    if (left < j)
        SortTreeCtr (item, left, j);
    if (i < right)
        SortTreeCtr (item, i, right);
}


/* StoreSumtTree: Store tree in treeList in packed format */
int StoreSumtTree (PackedTree *treeList, int index, PolyTree *t)
{
    int orderLen, numBrlens;

    assert (treeList[index].brlens == NULL);
    assert (treeList[index].order == NULL);

    /* get tree dimensions */
    numBrlens = t->nNodes - 1;
    orderLen = t->nIntNodes - 1;

    /* allocate space */
    treeList[index].brlens = (MrBFlt *) SafeCalloc (numBrlens, sizeof(MrBFlt));
    treeList[index].order  = (int *) SafeCalloc (orderLen, sizeof(MrBFlt));
    if (!treeList[index].order || !treeList[index].brlens)
        {
        MrBayesPrint ("%s   Could not store packed representation of tree '%s'\n", spacer, t->name);
        return (ERROR);
        }

    /* store tree */
    if (t->isRooted == YES)
        StoreRPolyTree (t, treeList[index].order, treeList[index].brlens);
    else
        StoreUPolyTree (t, treeList[index].order, treeList[index].brlens);

    return (NO_ERROR);
}


/* TreeCtrUppass: extract TreeCtr nodes in uppass sequence */
void TreeCtrUppass (TreeCtr *r, TreeCtr **uppass, int *index)
{
    if (r != NULL)
        {
        uppass[(*index)++] = r;

        TreeCtrUppass (r->left, uppass, index);
        TreeCtrUppass (r->right, uppass, index);
        }
}


int TreeProb (void)
{
    int         i, num, nInSets[5];
    MrBFlt      treeProb, cumTreeProb;
    TreeCtr     **trees;
    Tree        *theTree;
    
    /* check if we need to do this */
    if (sumtParams.calcTreeprobs == NO)
        return (NO_ERROR);

    MrBayesPrint ("%s   Calculating tree probabilities...\n\n", spacer);

    /* allocate space for tree counters and trees */
    trees = (TreeCtr **) SafeCalloc ((size_t)numUniqueTreesFound, sizeof(TreeCtr *));
    theTree = AllocateTree (sumtParams.numTaxa);
    if (!trees || !theTree)
        {
        MrBayesPrint ("%s   Problem allocating trees or theTree in TreeProb\n", spacer);
        return (ERROR);
        }
    
    /* extract trees */
    i = 0;
    TreeCtrUppass (treeCtrRoot, trees, &i);

    /* sort trees */
    SortTreeCtr (trees, 0, numUniqueTreesFound-1);
    
    /* set basic params in receiving tree */
    theTree->isRooted = sumtParams.isRooted;
    if (theTree->isRooted)
        {
        theTree->nNodes = 2 * sumtParams.numTaxa;
        theTree->nIntNodes = sumtParams.numTaxa - 1;
        }
    else
        {
        theTree->nNodes = 2 * sumtParams.numTaxa - 2;
        theTree->nIntNodes = sumtParams.numTaxa - 2;
        }

    /* show tree data */
    cumTreeProb = 0.0;
    nInSets[0] = nInSets[1] = nInSets[2] = nInSets[3] = nInSets[4] = 0;
    for (num=0; num<numUniqueTreesFound; num++)   /* loop over all of the trees that were found */
        {
        /* get probability of tree */
        treeProb = (MrBFlt)trees[num]->count / (MrBFlt)sumtParams.numTreesSampled;
        cumTreeProb += treeProb;
        if (cumTreeProb >= 0.0 && cumTreeProb < 0.5)
            nInSets[0]++;
        else if (cumTreeProb >= 0.5 && cumTreeProb < 0.9)
            nInSets[1]++;
        else if (cumTreeProb >= 0.9 && cumTreeProb < 0.95)
            nInSets[2]++;
        else if (cumTreeProb >= 0.95 && cumTreeProb < 0.99)
            nInSets[3]++;
        else
            nInSets[4]++;
        
        /* draw tree to stdout */
        if (theTree->isRooted == YES)
            RetrieveRTopology (theTree, trees[num]->order);
        else
            RetrieveUTopology (theTree, trees[num]->order);
        if (sumtParams.showSumtTrees == YES)
            {
            MrBayesPrint ("\n%s   Tree %d (p = %1.3lf, P = %1.3lf):\n\n", spacer, num+1, treeProb, cumTreeProb);
            ShowTree (theTree);
            }

        /* draw tree to file */
        if (num == 0)
            {
            MrBayesPrintf (fpTrees, "[This file contains the trees that were found during the MCMC\n");
            MrBayesPrintf (fpTrees, "search, sorted by posterior probability. \"p\" indicates the\n");
            MrBayesPrintf (fpTrees, "posterior probability of the tree whereas \"P\" indicates the\n");
            MrBayesPrintf (fpTrees, "cumulative posterior probability.]\n\n");
            MrBayesPrintf (fpTrees, "begin trees;\n");
            MrBayesPrintf (fpTrees, "   translate\n");
            for (i=0; i<sumtParams.numTaxa; i++)
                {
                if (i == sumtParams.numTaxa - 1)
                    MrBayesPrintf (fpTrees, "   %2d %s;\n", i+1, sumtParams.taxaNames[i]);
                else
                    MrBayesPrintf (fpTrees, "   %2d %s,\n", i+1, sumtParams.taxaNames[i]);
                }
            }
        MrBayesPrintf (fpTrees, "   tree tree_%d [p = %1.3lf, P = %1.3lf] = [&W %1.6lf] ", num+1, treeProb, cumTreeProb, treeProb);
        WriteTopologyToFile (fpTrees, theTree->root->left, theTree->isRooted);
        MrBayesPrintf (fpTrees, ";\n");
        if (num == numUniqueTreesFound - 1)
            MrBayesPrintf (fpTrees, "end;\n");  
        }
        
    /* print out general information on credible sets of trees */
    i = nInSets[0] + nInSets[1] + nInSets[2] + nInSets[3] + nInSets[4];
    MrBayesPrint ("%s   Credible sets of trees (%d tree%s sampled):\n", spacer, i, i > 1 ? "s" : "");
    i = nInSets[0] + 1;
    if (i > 1)
        MrBayesPrint ("%s      50 %% credible set contains %d trees\n", spacer, i);
    i += nInSets[1];
    if (i > 1)
        MrBayesPrint ("%s      90 %% credible set contains %d trees\n", spacer, i);
    i += nInSets[2];
    if (i > 1)
        MrBayesPrint ("%s      95 %% credible set contains %d trees\n", spacer, i);
    i += nInSets[3];
    MrBayesPrint ("%s      99 %% credible set contains %d tree%s\n\n", spacer, i, i > 1 ? "s" : "");
    
    /* free memory */
    free (trees);
        
    return (NO_ERROR);  
}


void WriteConTree (PolyNode *p, FILE *fp, int showSupport)
{
    PolyNode        *q;

    if (p->anc != NULL)
        if (p->anc->left == p)
            fprintf (fp, "(");

    for (q = p->left; q != NULL; q = q->sib)
        {
        if (q->anc->left != q) /* Note that q->anc always exists (it is p) */
            fprintf (fp, ",");
        WriteConTree (q, fp, showSupport);
        }
    if (p->left == NULL)
        {
        if (sumtParams.brlensDef == YES)
            {
            if (sumtParams.isClock == NO)
                fprintf (fp, "%d:%s", p->index+1, MbPrintNum(p->length));
            else
                fprintf (fp, "%d:%s", p->index+1, MbPrintNum(p->anc->depth - p->depth));
            }
        else
            fprintf (fp, "%d", p->index+1);
        }
        
    if (p->sib == NULL && p->anc != NULL)
        {
        if (p->anc->anc != NULL)
            {
            if (sumtParams.brlensDef == YES && showSupport == NO)
                {
                if (sumtParams.isClock == NO)
                    fprintf (fp, "):%s", MbPrintNum(p->anc->length));
                else
                    fprintf (fp, "):%s", MbPrintNum(p->anc->anc->depth - p->anc->depth));
                }
            else if (sumtParams.brlensDef == NO && showSupport == YES)
                fprintf (fp, ")%1.3lf", p->anc->support); 
            else if (sumtParams.brlensDef == YES && showSupport == YES)
                {
                if (sumtParams.isClock == NO)
                    fprintf (fp, ")%1.3lf:%s", p->anc->support, MbPrintNum(p->anc->length));
                else
                    fprintf (fp, ")%1.3lf:%s", p->anc->support, MbPrintNum(p->anc->anc->depth - p->anc->depth));
                }
            else
                fprintf (fp, ")");
            }
        else
            fprintf (fp, ")");
        }
}


/* WriteFigTreeConTree: Include rich information for each node in a consensus tree */
void WriteFigTreeConTree (PolyNode *p, FILE *fp, PartCtr **treeParts)
{
    PolyNode        *q;

    if (p->left == NULL)
        {
        fprintf (fp, "%d", p->index+1);
        if (sumtParams.isClock == NO)
            PrintFigTreeNodeInfo(fp, treeParts[p->partitionIndex], p->length);
        else if (sumtParams.isCalibrated == YES)
            PrintFigTreeNodeInfo(fp, treeParts[p->partitionIndex], p->anc->age - p->age);
        else
            PrintFigTreeNodeInfo(fp, treeParts[p->partitionIndex], p->anc->depth - p->depth);
        }
    else
        {
        fprintf  (fp, "(");
        for (q = p->left; q != NULL; q = q->sib)
            {
            WriteFigTreeConTree (q, fp, treeParts);
            if (q->sib != NULL)
                fprintf (fp, ",");
            }
        fprintf (fp, ")");
        if (p->partitionIndex >= 0 && p->partitionIndex < numUniqueSplitsFound)
            {
            if (p->anc == NULL)
                {
                if (sumtParams.isClock == YES)
                    PrintFigTreeNodeInfo(fp,treeParts[p->partitionIndex], -1.0);
                }
            else if (sumtParams.isClock == NO)
                PrintFigTreeNodeInfo(fp, treeParts[p->partitionIndex], p->length);
            else if (sumtParams.isCalibrated == YES)
                PrintFigTreeNodeInfo(fp, treeParts[p->partitionIndex], p->anc->age - p->age);
            else
                PrintFigTreeNodeInfo(fp, treeParts[p->partitionIndex], p->anc->depth - p->depth);
            }
        }
}

