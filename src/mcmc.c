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
#include "likelihood.h"
#include "mbbeagle.h"
#include "mcmc.h"
#include "model.h"
#include "proposal.h"
#include "sumpt.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif
#include <signal.h>

#if defined (WIN_VERSION) && !defined (__GNUC__)
#define VISUAL
#else
typedef void (*sighandler_t) (int);
#endif

#define GIBBS_SAMPLE_FREQ           100         /* generations between gibbs sampling of gamma cats */
#define MAX_SMALL_JUMP              10          /* threshold for precalculating trans probs of adgamma model */
#define BIG_JUMP                    100         /* threshold for using stationary approximation */
#define MAX_RUNS                    120         /* maximum number of independent runs */
#define PFILE                       0
#define TFILE                       1
#define CALFILE                     2
#define MCMCFILE                    3
#define MAXTUNINGPARAM              10000       /* limit to ensure convergence for autotuning */
#define SAMPLE_ALL_SS                           /* if defined makes ss sample every generation instead of every sample frequency */
#define BEAGLE_RESCALE_FREQ         160
#define BEAGLE_RESCALE_FREQ_DOUBLE  10          /* The factor by which BEAGLE_RESCALE_FREQ get multiplied if double presicion is used */
#define TARGETLENDELTA              100

/* debugging compiler statements */
#undef  DEBUG_SETUPTERMSTATE
#undef  DEBUG_RUNCHAIN
#undef  DEBUG_NOSHORTCUTS
#undef  DEBUG_NOSCALING
#undef  DEBUG_TIPROBS_STD
#undef  DEBUG_RUN_WITHOUT_DATA
#undef  DEBUG_CONSTRAINTS
#undef  DEBUG_LNLIKELIHOOD  /* slow if defined!! */
#undef  DEBUG_LIKELIHOOD
#undef  DEBUG_FBDPR       // #undef  FBDPR_CondOnN
#undef  SHOW_MOVE

#if defined (MPI_ENABLED)
#define ERROR_TEST2(failString,X1,X2) \
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);\
    if (sumErrors > 0)\
        {\
        MrBayesPrint ("%s   "failString"\n", spacer);\
        X1;X2;\
        }
#else
#define ERROR_TEST2(failString,X1,X2) \
    if (nErrors > 0)\
        {\
        MrBayesPrint ("%s   "failString"\n", spacer);\
        X1;X2;\
        }
#endif

/* local (to this file) data types */
typedef struct pfnode
    {
    struct pfnode   *left;
    struct pfnode   *right;
    int             *count;
    BitsLong        *partition;
    } PFNODE;

/* local prototypes */
int       AddTreeSamples (int from, int to, int saveToList);
PFNODE   *AddPartition (PFNODE *r, BitsLong *p, int runId);
int       AddTreeToPartitionCounters (Tree *tree, int treeId, int runId);
int       AttemptSwap (int swapA, int swapB, RandLong *seed);
void      BuildExhaustiveSearchTree (Tree *t, int chain, int nTaxInTree, TreeInfo *tInfo);
int       BuildStepwiseTree (Tree *t, int chain, RandLong *seed);
int       CalcLikeAdgamma (int d, Param *param, int chain, MrBFlt *lnL);
void      CalcPartFreqStats (PFNODE *p, STATS *stat);
void      CalcTopoConvDiagn (int numSamples);
#ifdef    VISUAL
BOOL      WINAPI CatchInterrupt (DWORD signum);
#else  
void      CatchInterrupt (int signum);
#endif  
int       CheckTemperature (void);
void      CloseMBPrintFiles (void);
PFNODE   *CompactTree (PFNODE *p);
int       ConfirmAbortRun(void);
void      CopyParams (int chain);
void      CopyPFNodeDown (PFNODE *p);
void      CopySiteScalers (ModelInfo *m, int chain);
void      CopyTrees (int chain);
int       ExtendChainQuery (void);
int       FillNumSitesOfPat (void);
TreeNode *FindBestNode (Tree *t, TreeNode *p, TreeNode *addNode, CLFlt *minLength, int chain);
void      FlipCijkSpace (ModelInfo *m, int chain);
void      FlipCondLikeSpace (ModelInfo *m, int chain, int nodeIndex);
void      FlipNodeScalerSpace (ModelInfo *m, int chain, int nodeIndex);
void      FlipSiteScalerSpace (ModelInfo *m, int chain);
void      FlipTiProbsSpace (ModelInfo *m, int chain, int nodeIndex);
void      FreeChainMemory (void);
MrBFlt    GetFitchPartials (ModelInfo *m, int chain, int source1, int source2, int destination);
void      GetStamp (void);
void      GetSwappers (int *swapA, int *swapB, int curGen);
void      GetTempDownPassSeq (TreeNode *p, int *i, TreeNode **dp);
MrBFlt    GibbsSampleGamma (int chain, int division, RandLong *seed);
int       InitAdGamma(void);
int       InitChainCondLikes (void);
int       InitClockBrlens (Tree *t);
int       InitEigenSystemInfo (ModelInfo *m);
int       InitInvCondLikes (void);
int       InitParsSets (void);
int       InitPrintParams (void);
int       IsPFNodeEmpty (PFNODE *p);
PFNODE   *LargestNonemptyPFNode (PFNODE *p, int *i, int j);
MrBFlt    LogLike (int chain);
MrBFlt    LogPrior (int chain);
int       LnBirthDeathPriorPrRandom    (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt sR, MrBFlt eR, MrBFlt sF);
int       LnBirthDeathPriorPrDiversity (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt sR, MrBFlt eR, MrBFlt sF);
int       LnBirthDeathPriorPrCluster   (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt sR, MrBFlt eR, MrBFlt sF);
int       LnFossilizedBDPriorFossilTip (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt *sR, MrBFlt *eR, MrBFlt sF, MrBFlt *fR);
int       LnFossilizedBDPriorRandom    (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt *sR, MrBFlt *eR, MrBFlt sF, MrBFlt *fR);
int       LnFossilizedBDPriorDiversity (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt *sR, MrBFlt *eR, MrBFlt sF, MrBFlt *fR);
MrBFlt    LnP0 (MrBFlt t, MrBFlt l, MrBFlt m);
MrBFlt    LnP0Subsample (MrBFlt t, MrBFlt l, MrBFlt m, MrBFlt f);
MrBFlt    LnP1 (MrBFlt t, MrBFlt l, MrBFlt m);
MrBFlt    LnP1Subsample (MrBFlt t, MrBFlt l, MrBFlt m, MrBFlt f);
MrBFlt    LnP0_fossil (MrBFlt t, MrBFlt lambda, MrBFlt mu, MrBFlt psi, MrBFlt c1, MrBFlt c2);
MrBFlt    LnP1_fossil (MrBFlt t, MrBFlt rho, MrBFlt c1, MrBFlt c2);
MrBFlt    LnQi_fossil (MrBFlt t, MrBFlt *t_f, int sl, MrBFlt *c1, MrBFlt *c2);
MrBFlt    LnPi_fossil (MrBFlt t, MrBFlt *t_f, int sl, MrBFlt *c1, MrBFlt *c2, MrBFlt *lambda, MrBFlt *mu, MrBFlt *psi);
int       NewtonRaphsonBrlen (Tree *t, TreeNode *p, int chain);
void      NodeToNodeDistances (Tree *t, TreeNode *fromNode);
int       PickProposal (RandLong *seed, int chainIndex);
int       NumCppEvents (Param *p, int chain);
int       PosSelProbs (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       PosSelProbs_SSE (TreeNode *p, int division, int chain);
#endif
int       PreparePrintFiles (void);
int       PrintAncStates_Bin (TreeNode *p, int division, int chain);
int       PrintAncStates_Gen (TreeNode *p, int division, int chain);
int       PrintAncStates_NUC4 (TreeNode *p, int division, int chain);
int       PrintAncStates_Std (TreeNode *p, int division, int chain);
int       PrintCalTree (int curGen, Tree *tree);
int       PrintCheckPoint (int gen);
int       PrintMCMCDiagnosticsToFile (int curGen);
#if defined (MPI_ENABLED)
int       PrintMPISlaves (FILE *fp);
#endif
void      PrintParamValues (Param *p, int chain, char *s);
int       PrintParsMatrix (void);
int       PrintSiteRates_Gen (TreeNode *p, int division, int chain);
int       PrintSiteRates_Std (TreeNode *p, int division, int chain);
int       PrintStates (int curGen, int coldId);
int       PrintStatesToFiles (int n);
int       PrintSwapInfo (void);
int       PrintTermState (void);
void      PrintTiProbs (CLFlt *tP, MrBFlt *bs, int nStates);
int       PrintTopConvInfo (void);
void      PrintToScreen (int curGen, int startGen, time_t endingT, time_t startingT);
int       PrintTree (int curGen, Param *treeParam, int chain, int showBrlens, MrBFlt clockRate);
MrBFlt    PropAncFossil (Param *param, int chain);
#if defined (MPI_ENABLED)
int       ReassembleMoveInfo (void);
int       ReassembleParamVals (int *curId);
int       ReassembleSwapInfo (void);
int       ReassembleTuningParams (void);
void      RedistributeMoveInfo (void);
int       RedistributeParamVals (void);
int       RedistributeTuningParams (void);
#endif
int       RemovePartition (PFNODE *r, BitsLong *p, int runId);
int       RemoveTreeFromPartitionCounters (Tree *tree, int treeId, int runId);
int       RemoveTreeSamples (int from, int to);
int       ReopenMBPrintFiles (void);
void      ResetChainIds (void);
void      ResetFlips(int chain);
void      ResetSiteScalers (ModelInfo *m, int chain);
int       ReusePreviousResults(int *numSamples, int);
int       RunChain (RandLong *seed);
int       SafeSprintf (char **target, int *targetLen, char *fmt, ...);
void      SetChainIds (void);
void      SetFileNames (void);
int       SetLikeFunctions (void);
int       SetLocalChainsAndDataSplits (void);
int       SetModelInfo (void);
int       SetMoves (void);
int       SetBinaryQMatrix (MrBFlt **a, int whichChain, int division);
int       SetNucQMatrix (MrBFlt **a, int n, int whichChain, int division, MrBFlt rateMult, MrBFlt *rA, MrBFlt *rS);
int       SetProteinQMatrix (MrBFlt **a, int n, int whichChain, int division, MrBFlt rateMult);
int       SetStdQMatrix (MrBFlt **a, int nStates, MrBFlt *bs, int cType);
int       SetUpPartitionCounters (void);
int       SetUpTermState (void);
int       SetUsedMoves (void);
int       ShowMoveSummary (void);
void      ShowValuesForChain (int chn);
int       SiteOmegas (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       SiteOmegas_SSE (TreeNode *p, int division, int chain);
#endif
PFNODE   *SmallestNonemptyPFNode (PFNODE *p, int *i, int j);
PFNODE   *Talloc (void);
void      Tfree (PFNODE *r);
MrBFlt    Temperature (int x);
void      TouchAllCijks (int chain);
void      TouchAllPartitions (void);
void      TouchAllTrees (int chain);
void      TouchEverything (int chain);

/* globals declared here and used elsewhere */
int             *bsIndex;                    /* compressed std stat freq index               */
Chain           chainParams;                 /* parameters of Markov chain                   */
int             *compCharPos;                /* char position in compressed matrix           */
int             *compColPos;                 /* column position in compressed matrix         */
BitsLong        *compMatrix;                 /* compressed character matrix                  */
int             compMatrixRowSize;           /* row size of compressed matrix                */
char            inputFileName[100];          /* input (NEXUS) file name                      */
MoveType        moveTypes[NUM_MOVE_TYPES];   /* holds information on the move types          */
int             numCompressedChars;          /* number of compressed characters              */
int             numMoveTypes;                /* the number of move types                     */
CLFlt           *numSitesOfPat;              /* no. sites of each pattern                    */
int             *origChar;                   /* index from compressed char to original char  */
char            stamp[11];                   /* holds a unique identifier for each analysis  */
MrBFlt          *stdStateFreqs;              /* std char state frequencies                   */
int             *stdType;                    /* compressed std char type: ord, unord, irrev  */
int             *tiIndex;                    /* compressed std char ti index                 */

#if defined (BEAGLE_ENABLED)
int             recalcScalers;               /* shoud we recalculate scalers for current state YES/NO */
#endif

/* globals used here but declared elsewhere (in likelihood.c) */
extern CLFlt     *preLikeL;                  /* precalculated cond likes for left descendant */
extern CLFlt     *preLikeR;                  /* precalculated cond likes for right descendant*/
extern CLFlt     *preLikeA;                  /* precalculated cond likes for ancestor        */

/* local (to this file) variables */
int             numLocalChains;              /* number of Markov chains                      */
int             *chainId = NULL;             /* information on the id (0 ...) of the chain   */
MrBFlt          *curLnL = NULL;              /* stores log likelihood                        */
MrBFlt          *curLnPr = NULL;             /* stores log prior probability                 */
int             stepRelativeBurninSS;        /* Should we use relative burn in within each step or not    */
MrBFlt          powerSS;                     /* power (betta) in power posterior destribution used in SS  */
MrBFlt          *marginalLnLSS = NULL;       /* marginal liklihood obtained using stepppingstone sampling */
MrBFlt          *stepAcumulatorSS = NULL;    /* accumulates liklihoods for current step in SS             */
MrBFlt          *stepScalerSS = NULL;        /* scaler of stepAcumulatorSS in log scale in SS             */
MrBFlt          *splitfreqSS = NULL;         /* array holding split frequencis for each step in SS        */
int             *sympiIndex;                 /* sympi state freq index for multistate chars  */
int             stdStateFreqsRowSize;        /* row size for std state frequencies           */
int             *weight;                     /* weight of each compressed char               */
int             *chainTempId;                /* info ton temp, change to float holding temp? */
int             state[MAX_CHAINS];           /* state of chain                               */
int             augmentData;                 /* are data being augmented for any division?   */
int             *nAccepted;                  /* counter of accepted moves                    */
int             *termState = NULL;           /* stores character states of tips              */
int             *isPartAmbig = NULL;         /* records whether tips are partially ambiguous */
BitsLong        **parsPtrSpace = NULL;       /* space holding pointers to parsimony sets     */
BitsLong        ***parsPtr = NULL;           /* pointers to pars state sets for chain & node */
CLFlt           *parsNodeLengthSpace = NULL; /* space for parsimony node lengths             */
CLFlt           **parsNodeLen = NULL;        /* pointers to pars node lengths for chains     */
char            *printString;                /* string for printing to a file                */
size_t          printStringSize;             /* length of printString                        */
MCMCMove        **usedMoves;                 /* vector of pointers to used moves             */
int             numUsedMoves;                /* the number of moves used by chain            */
Param           **printParam;                /* vector of pointers to normal params to print */
int             numPrintParams;              /* the number of normal params to print         */
Param           **printTreeParam;            /* vector of pointers to tree params to print   */
Param           **topologyParam;             /* vector of pointers to topology params        */
int             numPrintTreeParams;          /* the number of tree params to print           */
int             codon[6][64];                /* holds info on amino acids coded in code      */
int             chainHasAdgamma;             /* indicates if chain has adgamma HMMs          */
int             inferPosSel;                 /* indicates if positive selection is inferred  */
MrBFlt          *posSelProbs;                /* probs. for positive selection                */
int             hasMarkovTi[MAX_SMALL_JUMP]; /* vector marking size of observed HMM jumps    */
int             *siteJump;                   /* vector of sitejumps for adgamma model        */
MrBFlt          **rateProbs;                 /* pointers to rate probs used by adgamma model */
MrBFlt          *rateProbSpace;              /* space for rate probs used by adgamma model   */
int             rateProbRowSize;             /* size of rate probs for one chain one state   */
MrBFlt          **markovTi[MAX_SMALL_JUMP];  /* trans prob matrices used in calc of adgamma  */
MrBFlt          **markovTiN;                 /* trans prob matrices used in calc of adgamma  */
int             whichReweightNum;            /* used for setting reweighting of char pats    */
int             ***swapInfo;                 /* keeps track of attempts & successes of swaps */
int             tempIndex;                   /* keeps track of which user temp is specified  */
int             abortMove;                   /* flag determining whether to abort move       */
PFNODE          **partFreqTreeRoot;          /* root of tree(s) holding partition freqs      */
int             nLongsNeeded;                /* number of longs needed for partitions        */
BitsLong        **partition;                 /* matrix holding partitions                    */
MrBFlt          *maxLnL0 = NULL;             /* maximum likelihood                           */
FILE            *fpMcmc = NULL;              /* pointer to .mcmc file                        */
FILE            **fpParm = NULL;             /* pointer to .p file(s)                        */
FILE            ***fpTree = NULL;            /* pointer to .t file(s)                        */
FILE            *fpSS = NULL;                /* pointer to .ss file                          */
static int      requestAbortRun;             /* flag for aborting mcmc analysis              */
int             *topologyPrintIndex;         /* print file index of each topology            */
int             *printTreeTopologyIndex;     /* topology index of each tree print file       */
int             numPreviousGen;              /* number of generations in run to append to    */

#if defined (MPI_ENABLED)
int             lowestLocalRunId;            /* lowest local run Id                          */
int             highestLocalRunId;           /* highest local run Id                         */
#endif

#if defined (PRINT_DUMP)
FILE            **fpDump = NULL;             /* pointer to .dump file(s)                     */
#endif

/* AddPartition: Add a partition to the tree keeping track of partition frequencies */
PFNODE *AddPartition (PFNODE *r, BitsLong *p, int runId)
{
    int     i, comp;
    
    if (r == NULL)
        {
        /* new partition */
        r = Talloc ();                  /* create a new node */
        if (r == NULL)
            return NULL;
        for (i=0; i<nLongsNeeded; i++)
            r->partition[i] = p[i];
        for (i=0; i<chainParams.numRuns; i++)
            r->count[i] = 0;
        r->count[runId] = 1;
        r->left = r->right = NULL;
        }
    else
        {
        for (i=0; i<nLongsNeeded; i++)
            {
            if (r->partition[i] != p[i])
                break;
            }
        
        if (i == nLongsNeeded)
            comp = 0;
        else if (r->partition[i] < p[i])
            comp = -1;
        else
            comp = 1;
        
        if (comp == 0)          /* repeated partition */
            r->count[runId]++;
        else if (comp < 0)      /* greater than -> into left subtree */
            {
            if ((r->left = AddPartition (r->left, p, runId)) == NULL)
                {
                Tfree (r);
                return NULL;
                }
            }
        else
            {
            /* smaller than -> into right subtree */
            if ((r->right = AddPartition (r->right, p, runId)) == NULL)
                {
                Tfree (r);
                return NULL;
                }
            }
        }

    return r;
}


int AddToPrintString (char *tempStr)
{
    size_t  len1, len2;
    
    len1 = strlen(printString);
    len2 = strlen(tempStr);
    if (len1 + len2 + 5 > printStringSize)
        {
        printStringSize += len1 + len2 - printStringSize + 200;
        printString = (char*)SafeRealloc((void*)printString, printStringSize * sizeof(char));
        if (!printString)
            {
            MrBayesPrint ("%s   Problem reallocating printString (%d)\n", spacer, printStringSize * sizeof(char));
            goto errorExit;
            }
        }
    strcat(printString, tempStr);   
    // printf ("printString(%d) -> \"%s\"\n", printStringSize, printString);
    return (NO_ERROR);
    
    errorExit:
        return (ERROR);
}


/* AddTreeSamples: Add tree samples from .t files to partition counters. if saveToList == YES then also save trees in tree list */
int AddTreeSamples (int from, int to, int saveToList)
{
    int     i, j, k, longestLine;
    BitsLong    lastBlock;
    char    *word, *s, *lineBuf;
    FILE    *fp;
    Tree    *t;
    char    *tempStr;
    int     tempStrSize = TEMPSTRSIZE;

    if (from > to)
        return (NO_ERROR);

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return (NO_ERROR);
#   endif

    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    for (i=0; i<numTopologies; i++)
        {
        t = chainParams.dtree;

        for (j=0; j<chainParams.numRuns; j++)
            {
            if (numPrintTreeParams == 1)
                {
                if (chainParams.numRuns == 1)
                    SafeSprintf (&tempStr, &tempStrSize, "%s.t", chainParams.chainFileName);
                else
                    SafeSprintf (&tempStr, &tempStrSize, "%s.run%d.t", chainParams.chainFileName, j+1);
                }
            else
                {
                if (chainParams.numRuns == 1)
                    SafeSprintf (&tempStr, &tempStrSize, "%s.tree%d.t", chainParams.chainFileName, i+1);
                else
                    SafeSprintf (&tempStr, &tempStrSize, "%s.tree%d.run%d.t", chainParams.chainFileName, i+1, j+1);
                }

            if ((fp = OpenBinaryFileR (tempStr)) == NULL) 
                {
                MrBayesPrint ("%s   Problem openning file %s.\n", spacer, tempStr);
                free (tempStr);
                return (ERROR);
                }

            longestLine = LongestLine (fp);
            SafeFclose (&fp);

            if ((fp = OpenTextFileR (tempStr)) == NULL) 
                {
                free (tempStr);
                return (ERROR);
                }

            lineBuf = (char *) SafeCalloc (longestLine + 10, sizeof (char));
            if (!lineBuf)
                {
                SafeFclose (&fp);
                free (tempStr);
                return (ERROR);
                }

            lastBlock = LastBlock (fp, lineBuf, longestLine);
            fseek (fp, lastBlock, SEEK_SET);

            for (k=1; k<=to; k++)
                {
                do {
                    if (fgets (lineBuf, longestLine, fp) == NULL) 
                        {
                        SafeFclose (&fp);
                        free (lineBuf);
                        free (tempStr);
                        return ERROR;
                        }
                    word = strtok (lineBuf, " ");
                    } while (strcmp (word, "tree") != 0);
                if (k>=from)
                    {
                    s = strtok (NULL, ";");
                    while (*s != '(')
                        s++;
                    StripComments(s);
                    if (ResetTopology (t, s) == ERROR)
                        {
                        SafeFclose (&fp);
                        free (lineBuf);
                        free (tempStr);
                        return ERROR;
                        }
                    if (AddTreeToPartitionCounters (t, i, j) == ERROR)
                        {
                        SafeFclose (&fp);
                        free (lineBuf);
                        free (tempStr);
                        return ERROR;
                        }
                    if (saveToList == YES)
                        if (AddToTreeList(&chainParams.treeList[numTopologies*j+i],t) == ERROR)
                            return (ERROR);
                    }
                }
            SafeFclose (&fp);
            free (lineBuf);
            } /* next run */
        } /* next tree */
    free (tempStr);
    return (NO_ERROR);
}


/* AddTreeToPartitionCounters: Break a tree into partitions and add those to counters */
int AddTreeToPartitionCounters (Tree *tree, int treeId, int runId)
{
    int         i, j, nTaxa;
    TreeNode    *p;

    if (tree->isRooted == YES)
        nTaxa = tree->nNodes - tree->nIntNodes - 1;
    else
        nTaxa = tree->nNodes - tree->nIntNodes;

    for (i=0; i<nTaxa; i++)
        {
        ClearBits(partition[i], nLongsNeeded);
        SetBit(i, partition[i]);
        }

    for (i=0; i<tree->nIntNodes-1; i++)
        {
        p = tree->intDownPass[i];
        assert (p->index >= tree->nNodes - tree->nIntNodes - (tree->isRooted == YES ? 1 : 0));
        for (j=0; j<nLongsNeeded; j++)
            {
            partition[p->index][j] = partition[p->left->index][j] | partition[p->right->index][j];
            }

        if ((partFreqTreeRoot[treeId] = AddPartition (partFreqTreeRoot[treeId], partition[p->index], runId)) == NULL)
            {
            MrBayesPrint ("%s   Could not allocate space for new partition in AddTreeToPartitionCounters\n", spacer);
            return ERROR;
            }
        }

    return NO_ERROR;
}


int AttemptSwap (int swapA, int swapB, RandLong *seed)
{
    int             d, tempX, reweightingChars, isSwapSuccessful, chI, chJ, runId;
    MrBFlt          tempA, tempB, lnLikeA, lnLikeB, lnPriorA, lnPriorB, lnR, r,
                    lnLikeStateAonDataB=0.0, lnLikeStateBonDataA=0.0, lnL;
    ModelInfo       *m;
    Tree            *tree;
#   if defined (MPI_ENABLED)
    int             numChainsForProc, tempIdA=0, tempIdB=0, proc, procIdForA=0, procIdForB=0, tempIdMy=0, procIdPartner=0,
                    whichElementA=0, whichElementB=0, lower, upper, areWeA, doISwap, ierror,
                    myId, partnerId, i, run;
    MrBFlt          swapRan;
    MPI_Status      status[2];
    MPI_Request     request[2];
#   endif

#   if defined (MPI_ENABLED)

    /* get the number of chains handled by this proc */
    /* the number will be corrected further down for unbalanced scenarios */
    numChainsForProc = (int) (chainParams.numChains * chainParams.numRuns / num_procs);

#   endif

    /* are we using character reweighting? */
    reweightingChars = NO;
    if ((chainParams.weightScheme[0] + chainParams.weightScheme[1]) > 0.00001)
        reweightingChars = YES;
            
#   if defined (MPI_ENABLED)

    /* figure out processors involved in swap */
    lower = upper = 0;
    for (proc=0; proc<num_procs; proc++)
        {
        /* assign or increment chain id */
        if (proc < (chainParams.numChains * chainParams.numRuns) % num_procs)
            upper += numChainsForProc+1;
        else
            upper += numChainsForProc;

        /* if swapA lies between lower and upper
            * chain id's we know that this is the proc
            * swapA is in */
        if (swapA >= lower && swapA < upper)
            {
            procIdForA = proc;
            whichElementA = swapA - lower;
            }
        if (swapB >= lower && swapB < upper)
            {
            procIdForB = proc;
            whichElementB = swapB - lower;
            }
        lower = upper;
        }

    /* NOTE: at this point, procIdForA and procIdForB *
        * store the proc id's of swapping procs. Also,   *
        * whichElementA and whichElementB store the      *
        * chainId[] index of swapping procs              */

    /* figure out if I am involved in the swap */
    doISwap = areWeA = NO;
    if (proc_id == procIdForA)
        {
        doISwap = YES;
        areWeA = YES;
        }
    else if (proc_id == procIdForB)
        {
        doISwap = YES;
        }

    /* chain's that do not swap, continue to the next iteration */  
    if (doISwap == YES)
        {
        
        /* no need to communicate accross processors if swapping chains are in the same proc */
        if (procIdForA == procIdForB)
            {
            if (reweightingChars == YES)
                {
                /* use character reweighting */
                lnLikeStateAonDataB = 0.0;
                for (d=0; d<numCurrentDivisions; d++)
                    {
                    m = &modelSettings[d];
                    tree = GetTree(m->brlens, whichElementA, state[whichElementA]);
                    lnL = 0.0;
                    m->Likelihood (tree->root->left, d, whichElementA, &lnL, chainId[whichElementB] % chainParams.numChains);
                    lnLikeStateAonDataB += lnL;
                    }
                lnLikeStateBonDataA = 0.0;
                for (d=0; d<numCurrentDivisions; d++)
                    {
                    m = &modelSettings[d];
                    tree = GetTree(m->brlens, whichElementB, state[whichElementB]);
                    lnL = 0.0;
                    m->Likelihood (tree->root->left, d, whichElementB, &lnL, chainId[whichElementA] % chainParams.numChains);
                    lnLikeStateBonDataA += lnL;
                    }
                }

            /*curLnPr[whichElementA] = LogPrior(whichElementA);
            curLnPr[whichElementB] = LogPrior(whichElementB);*/

            /* then do the serial thing - simply swap chain id's */
            tempA = Temperature (chainId[whichElementA]);
            tempB = Temperature (chainId[whichElementB]);
            lnLikeA = curLnL[whichElementA];
            lnLikeB = curLnL[whichElementB];
            if (chainParams.isSS == YES)
                {
                lnLikeA *= powerSS;
                lnLikeB *= powerSS;
                }
            lnPriorA = curLnPr[whichElementA];
            lnPriorB = curLnPr[whichElementB];
            if (reweightingChars == YES)
                {
                if (chainParams.isSS == YES)
                    lnR = (tempB * (lnLikeStateAonDataB*powerSS + lnPriorA) + tempA * (lnLikeStateBonDataA*powerSS + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
                else
                    lnR = (tempB * (lnLikeStateAonDataB + lnPriorA) + tempA * (lnLikeStateBonDataA + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
                }
            else
                lnR = (tempB * (lnLikeA + lnPriorA) + tempA * (lnLikeB + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
            if (lnR <  -100.0)
                r =  0.0;
            else if (lnR > 0.0)
                r =  1.0;
            else
                r =  exp(lnR);

            isSwapSuccessful = NO;
            if (RandomNumber(seed) < r)
                {
                /* swap chain id's (heats) */
                tempX = chainId[whichElementA];
                chainId[whichElementA] = chainId[whichElementB];
                chainId[whichElementB] = tempX;
                if (reweightingChars == YES)
                    {
                    curLnL[whichElementA] = lnLikeStateAonDataB;
                    curLnL[whichElementB] = lnLikeStateBonDataA;
                    }
                isSwapSuccessful = YES;
                }
                
            chI = chainId[whichElementA];
            chJ = chainId[whichElementB];
            if (chainId[whichElementB] < chainId[whichElementA])
                {
                chI = chainId[whichElementB];
                chJ = chainId[whichElementA];
                }
            runId = chI / chainParams.numChains;
            chI = chI % chainParams.numChains;
            chJ = chJ % chainParams.numChains;
            swapInfo[runId][chJ][chI]++;
            if (isSwapSuccessful == YES)
                swapInfo[runId][chI][chJ]++;
            }
        /* we need to communicate across processors */
        else
            {
            if (reweightingChars == YES)
                {
                /* If we are reweighting characters, then we need to do an additional communication to
                    figure out the chainId's of the partner. We need to have this information so we can
                    properly calculate likelihoods with switched observations. */
                if (areWeA == YES)
                    {
                    lnLikeStateAonDataB = 0.0;
                    myId = chainId[whichElementA];
                    ierror = MPI_Isend (&myId, 1, MPI_INT, procIdForB, 0, MPI_COMM_WORLD, &request[0]);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Irecv (&partnerId, 1, MPI_INT, procIdForB, 0, MPI_COMM_WORLD, &request[1]);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (2, request, status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    for (d=0; d<numCurrentDivisions; d++)
                        {
                        m = &modelSettings[d];
                        tree = GetTree(m->brlens, whichElementA, state[whichElementA]);
                        lnL = 0.0;
                        m->Likelihood (tree->root->left, d, whichElementA, &lnL, partnerId);
                        lnLikeStateAonDataB = lnL;
                        }
                    }
                else
                    {
                    lnLikeStateBonDataA = 0.0;
                    myId = chainId[whichElementB];
                    ierror = MPI_Isend (&myId, 1, MPI_INT, procIdForA, 0, MPI_COMM_WORLD, &request[0]);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Irecv (&partnerId, 1, MPI_INT, procIdForA, 0, MPI_COMM_WORLD, &request[1]);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (2, request, status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    for (d=0; d<numCurrentDivisions; d++)
                        {
                        m = &modelSettings[d];
                        tree = GetTree(m->brlens, whichElementB, state[whichElementB]);
                        lnL = 0.0;
                        m->Likelihood (tree->root->left, d, whichElementB, &lnL, partnerId);
                        lnLikeStateBonDataA = lnL;
                        }
                    }
                }

            if (areWeA == YES)
                {
                /*curLnPr[whichElementA] = LogPrior(whichElementA);*/

                /* we are processor A */
                tempIdA = chainId[whichElementA];
                lnLikeA = curLnL[whichElementA];
                lnPriorA = curLnPr[whichElementA];
                swapRan = RandomNumber(seed);

                myStateInfo[0] = lnLikeA;
                myStateInfo[1] = lnPriorA;
                myStateInfo[2] = tempIdA;
                myStateInfo[3] = swapRan;
                myStateInfo[4] = 0.0;
                if (reweightingChars == YES)
                    {
                    myStateInfo[2] = lnLikeStateAonDataB;
                    tempIdB = partnerId;
                    }
                    
                ierror = MPI_Isend (&myStateInfo, 5, MPI_DOUBLE, procIdForB, 0, MPI_COMM_WORLD, &request[0]);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Irecv (&partnerStateInfo, 5, MPI_DOUBLE, procIdForB, 0, MPI_COMM_WORLD, &request[1]);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Waitall (2, request, status);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }

                lnLikeA = curLnL[whichElementA];
                lnLikeB = partnerStateInfo[0];
                if (chainParams.isSS == YES)
                    {
                    lnLikeA *= powerSS;
                    lnLikeB *= powerSS;
                    }
                lnPriorA = curLnPr[whichElementA];
                lnPriorB = partnerStateInfo[1];
                if (reweightingChars == YES)
                    lnLikeStateBonDataA = partnerStateInfo[2];
                else
                    tempIdB = partnerStateInfo[2];
                
                tempA = Temperature (tempIdA);
                tempB = Temperature (tempIdB);

                if (reweightingChars == YES)
                    {
                    if (chainParams.isSS == YES)
                        lnR = (tempB * (lnLikeStateAonDataB*powerSS + lnPriorA) + tempA * (lnLikeStateBonDataA*powerSS + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
                    else
                        lnR = (tempB * (lnLikeStateAonDataB + lnPriorA) + tempA * (lnLikeStateBonDataA + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
                    }
                else
                    lnR = (tempB * (lnLikeA + lnPriorA) + tempA * (lnLikeB + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
                if (lnR < -100.0)
                    r = 0.0;
                else if (lnR > 0.0)
                    r = 1.0;
                else
                    r = exp(lnR);

                /* process A's random number is used to make the swap decision */
                isSwapSuccessful = NO;
                if (swapRan < r)
                    {
                    /* swap chain id's (heats) */
                    isSwapSuccessful = YES;
                    tempIdMy = chainId[whichElementA];
                    procIdPartner = procIdForB;

                    if (reweightingChars == YES)
                        chainId[whichElementA] = tempIdB;
                    else
                        chainId[whichElementA] = (int)(partnerStateInfo[2]);
                    if (reweightingChars == YES)
                        {
                        curLnL[whichElementA] = lnLikeStateAonDataB;
                        }
                    }
                    
                /* only processor A keeps track of the swap success/failure */
                chI = tempIdA;
                chJ = tempIdB;
                if (tempIdB < tempIdA)
                    {
                    chI = tempIdB;
                    chJ = tempIdA;
                    }
                runId = chI / chainParams.numChains;
                chI = chI % chainParams.numChains;
                chJ = chJ % chainParams.numChains;
                swapInfo[runId][chJ][chI]++;
                if (isSwapSuccessful == YES)
                    {
                    swapInfo[runId][chI][chJ]++;
                    /* exchange the move info */
                    for (i=0; i<numUsedMoves; i++)
                        {
                        myStateInfo[0] = usedMoves[i]->nAccepted[tempIdA];
                        myStateInfo[1] = usedMoves[i]->nTried[tempIdA];
                        myStateInfo[2] = usedMoves[i]->nBatches[tempIdA];
                        myStateInfo[3] = usedMoves[i]->nTotAccepted[tempIdA];
                        myStateInfo[4] = usedMoves[i]->nTotTried[tempIdA];
                        myStateInfo[5] = usedMoves[i]->lastAcceptanceRate[tempIdA];
                        if (usedMoves[i]->moveType->numTuningParams > 0)
                            myStateInfo[6] = usedMoves[i]->tuningParam[tempIdA][0];
                        else
                            myStateInfo[6] = 0.0;

                        ierror = MPI_Isend (&myStateInfo, 7, MPI_DOUBLE, procIdForB, 0, MPI_COMM_WORLD, &request[0]);
                        if (ierror != MPI_SUCCESS)
                            {
                            return (ERROR);
                            }
                        ierror = MPI_Irecv (&partnerStateInfo, 7, MPI_DOUBLE, procIdForB, 0, MPI_COMM_WORLD, &request[1]);
                        if (ierror != MPI_SUCCESS)
                            {
                            return (ERROR);
                            }
                        ierror = MPI_Waitall (2, request, status);
                        if (ierror != MPI_SUCCESS)
                            {
                            return (ERROR);
                            }

                        usedMoves[i]->nAccepted[tempIdB]          = (int)partnerStateInfo[0];
                        usedMoves[i]->nTried[tempIdB]             = (int)partnerStateInfo[1];
                        usedMoves[i]->nBatches[tempIdB]           = (int)partnerStateInfo[2];
                        usedMoves[i]->nTotAccepted[tempIdB]       = (int)partnerStateInfo[3];
                        usedMoves[i]->nTotTried[tempIdB]          = (int)partnerStateInfo[4];
                        usedMoves[i]->lastAcceptanceRate[tempIdB] = partnerStateInfo[5];
                        if (usedMoves[i]->moveType->numTuningParams > 0)
                            usedMoves[i]->tuningParam[tempIdB][0]     = partnerStateInfo[6];

                        usedMoves[i]->nAccepted[tempIdA]          = 0;
                        usedMoves[i]->nTried[tempIdA]             = 0;
                        usedMoves[i]->nBatches[tempIdA]           = 0;
                        usedMoves[i]->lastAcceptanceRate[tempIdA] = 0.0;
                        usedMoves[i]->nTotAccepted[tempIdA]       = 0;
                        usedMoves[i]->nTotTried[tempIdA]          = 0;
                        if (usedMoves[i]->moveType->numTuningParams > 0)
                            usedMoves[i]->tuningParam[tempIdA][0]     = 0.0;
                
                        }
                    }
                }
            else
                {
                /*curLnPr[whichElementB] = LogPrior(whichElementB);*/

                /* we are processor B */
                tempIdB  = chainId[whichElementB];
                lnLikeB  = curLnL[whichElementB];
                lnPriorB = curLnPr[whichElementB];
                swapRan  = -1.0;

                myStateInfo[0] = lnLikeB;
                myStateInfo[1] = lnPriorB;
                myStateInfo[2] = tempIdB;
                myStateInfo[3] = swapRan;
                myStateInfo[4] = 0.0;
                if (reweightingChars == YES)
                    {
                    myStateInfo[2] = lnLikeStateBonDataA;
                    tempIdA = partnerId;
                    }

                ierror = MPI_Isend (&myStateInfo, 5, MPI_DOUBLE, procIdForA, 0, MPI_COMM_WORLD, &request[0]);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Irecv (&partnerStateInfo, 5, MPI_DOUBLE, procIdForA, 0, MPI_COMM_WORLD, &request[1]);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Waitall (2, request, status);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }

                lnLikeB = curLnL[whichElementB];
                lnLikeA = partnerStateInfo[0];
                lnPriorB = curLnPr[whichElementB];
                lnPriorA = partnerStateInfo[1];
                if (reweightingChars == YES)
                    lnLikeStateAonDataB = partnerStateInfo[2];
                else
                    tempIdA = partnerStateInfo[2];

                tempB = Temperature (tempIdB);
                tempA = Temperature (tempIdA);

                if (chainParams.isSS == YES)
                    {
                    lnLikeA *= powerSS;
                    lnLikeB *= powerSS;
                    }

                if (reweightingChars == YES)
                    {
                    if (chainParams.isSS == YES)
                        lnR = (tempB * (lnLikeStateAonDataB*powerSS + lnPriorA) + tempA * (lnLikeStateBonDataA*powerSS + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
                    else
                        lnR = (tempB * (lnLikeStateAonDataB + lnPriorA) + tempA * (lnLikeStateBonDataA + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
                    }
                else
                    lnR = (tempB * (lnLikeA + lnPriorA) + tempA * (lnLikeB + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
                if (lnR < -100.0)
                    r = 0.0;
                else if (lnR > 0.0)
                    r = 1.0;
                else
                    r = exp(lnR);

                /* we use process A's random number to make the swap decision */
                isSwapSuccessful = NO;
                if (partnerStateInfo[3] < r)
                    {
                    isSwapSuccessful = YES;
                    tempIdMy = chainId[whichElementB];
                    procIdPartner = procIdForA;

                    if (reweightingChars == YES)
                        chainId[whichElementB] = tempIdA;
                    else
                        chainId[whichElementB] = (int)(partnerStateInfo[2]);
                    if (reweightingChars == YES)
                        {
                        curLnL[whichElementB] = lnLikeStateBonDataA;
                        }
                    /* swap the move info */
                    for (i=0; i<numUsedMoves; i++)
                        {
                        myStateInfo[0] = usedMoves[i]->nAccepted[tempIdB];
                        myStateInfo[1] = usedMoves[i]->nTried[tempIdB];
                        myStateInfo[2] = usedMoves[i]->nBatches[tempIdB];
                        myStateInfo[3] = usedMoves[i]->nTotAccepted[tempIdB];
                        myStateInfo[4] = usedMoves[i]->nTotTried[tempIdB];
                        myStateInfo[5] = usedMoves[i]->lastAcceptanceRate[tempIdB];
                        if (usedMoves[i]->moveType->numTuningParams > 0)
                            myStateInfo[6] = usedMoves[i]->tuningParam[tempIdB][0];
                        else
                            myStateInfo[6] = 0.0;

                        ierror = MPI_Isend (&myStateInfo, 7, MPI_DOUBLE, procIdForA, 0, MPI_COMM_WORLD, &request[0]);
                        if (ierror != MPI_SUCCESS)
                            {
                            return (ERROR);
                            }
                        ierror = MPI_Irecv (&partnerStateInfo, 7, MPI_DOUBLE, procIdForA, 0, MPI_COMM_WORLD, &request[1]);
                        if (ierror != MPI_SUCCESS)
                            {
                            return (ERROR);
                            }
                        ierror = MPI_Waitall (2, request, status);
                        if (ierror != MPI_SUCCESS)
                            {
                            return (ERROR);
                            }

                        usedMoves[i]->nAccepted[tempIdA]          = (int)partnerStateInfo[0];
                        usedMoves[i]->nTried[tempIdA]             = (int)partnerStateInfo[1];
                        usedMoves[i]->nBatches[tempIdA]           = (int)partnerStateInfo[2];
                        usedMoves[i]->nTotAccepted[tempIdA]       = (int)partnerStateInfo[3];
                        usedMoves[i]->nTotTried[tempIdA]          = (int)partnerStateInfo[4];
                        usedMoves[i]->lastAcceptanceRate[tempIdA] = partnerStateInfo[5];
                        if (usedMoves[i]->moveType->numTuningParams > 0)
                            usedMoves[i]->tuningParam[tempIdA][0]     = partnerStateInfo[6];

                        usedMoves[i]->nAccepted[tempIdB]          = 0;
                        usedMoves[i]->nTried[tempIdB]             = 0;
                        usedMoves[i]->nBatches[tempIdB]           = 0;
                        usedMoves[i]->nTotAccepted[tempIdB]       = 0;
                        usedMoves[i]->nTotTried[tempIdB]          = 0;
                        usedMoves[i]->lastAcceptanceRate[tempIdB] = 0.0;
                        if (usedMoves[i]->moveType->numTuningParams > 0)
                            usedMoves[i]->tuningParam[tempIdB][0]     = 0.0;
                        }
                    }
                }

            /*We exchange only if swap successful and (my id is cold or patner id is cold)*/
            if (chainParams.isSS == YES && isSwapSuccessful == YES && (tempIdMy % chainParams.numChains == 0 || (areWeA == YES && chainId[whichElementA] % chainParams.numChains == 0) || (areWeA == NO && chainId[whichElementB] % chainParams.numChains == 0)))
                {
                    run = tempIdMy/chainParams.numChains;

                    myStateInfo[0] = tempIdMy;
                    myStateInfo[1] = marginalLnLSS   [ run ];
                    myStateInfo[2] = stepAcumulatorSS[ run ];
                    myStateInfo[3] = stepScalerSS    [ run ];

                    ierror = MPI_Isend (&myStateInfo, 4, MPI_DOUBLE, procIdPartner, 0, MPI_COMM_WORLD, &request[0]);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }

                    ierror = MPI_Irecv (&partnerStateInfo, 4, MPI_DOUBLE, procIdPartner, 0, MPI_COMM_WORLD, &request[1]);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (2, request, status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }

                    /*we swap chains from the same run*/
                    assert (run == (int)partnerStateInfo[0]/chainParams.numChains);

                    /*If my chain is the cold chain then send current SS values of corresponding run*/
                    if (tempIdMy % chainParams.numChains == 0)
                        {
                        assert ((int)partnerStateInfo[0] % chainParams.numChains != 0);
                        assert (partnerStateInfo[1] == 0.0);
                        marginalLnLSS   [ run ] = (MrBFlt) 0.0;
                        stepAcumulatorSS[ run ] = (MrBFlt) 0.0;
                        stepScalerSS    [ run ] = (MrBFlt) 0.0;
                        }
                    else if ((int)partnerStateInfo[0] % chainParams.numChains == 0)
                        {
                        marginalLnLSS   [ run ] = (MrBFlt) partnerStateInfo[1];
                        stepAcumulatorSS[ run ] = (MrBFlt) partnerStateInfo[2];
                        stepScalerSS    [ run ] = (MrBFlt) partnerStateInfo[3];
                        }
                }

            }
        }
#   else
    if (reweightingChars == YES)
        {
        /* use character reweighting */
        lnLikeStateAonDataB = 0.0;
        for (d=0; d<numCurrentDivisions; d++)
            {
            m = &modelSettings[d];
            tree = GetTree(m->brlens, swapA, state[swapA]);
            lnL = 0.0;
            m->Likelihood (tree->root->left, d, swapA, &lnL, chainId[swapB] % chainParams.numChains);
            lnLikeStateAonDataB += lnL;
            }
        lnLikeStateBonDataA = 0.0;
        for (d=0; d<numCurrentDivisions; d++)
            {
            m = &modelSettings[d];
            tree = GetTree(m->brlens, swapB, state[swapB]);
            lnL = 0.0;
            m->Likelihood (tree->root->left, d, swapB, &lnL, chainId[swapA] % chainParams.numChains);
            lnLikeStateBonDataA += lnL;
            }
        }

    assert (fabs((curLnPr[swapA]-LogPrior(swapA))/curLnPr[swapA]) < 0.0001);
    assert (fabs((curLnPr[swapB]-LogPrior(swapB))/curLnPr[swapB]) < 0.0001);

    tempA = Temperature (chainId[swapA]);
    tempB = Temperature (chainId[swapB]);
    lnLikeA = curLnL[swapA];
    lnLikeB = curLnL[swapB];
    lnPriorA = curLnPr[swapA];
    lnPriorB = curLnPr[swapB];

    if (chainParams.isSS == YES)
        {
        lnLikeA *= powerSS;
        lnLikeB *= powerSS;
        }

    if (reweightingChars == YES)
        {
        if (chainParams.isSS == YES)
            lnR = (tempB * (lnLikeStateAonDataB*powerSS + lnPriorA) + tempA * (lnLikeStateBonDataA*powerSS + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
        else
            lnR = (tempB * (lnLikeStateAonDataB + lnPriorA) + tempA * (lnLikeStateBonDataA + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
        }
    else
        lnR = (tempB * (lnLikeA + lnPriorA) + tempA * (lnLikeB + lnPriorB)) - (tempA * (lnLikeA + lnPriorA) + tempB * (lnLikeB + lnPriorB));
    if (lnR < -100.0)
        r = 0.0;
    else if (lnR > 0.0)
        r =  1.0;
    else
        r =  exp (lnR);

    isSwapSuccessful = NO;
    if (RandomNumber(seed) < r)
        {
        tempX = chainId[swapA];
        chainId[swapA] = chainId[swapB];
        chainId[swapB] = tempX;

        if (reweightingChars == YES)
            {
            curLnL[swapA] = lnLikeStateAonDataB;
            curLnL[swapB] = lnLikeStateBonDataA;
            }
        isSwapSuccessful = YES;
        }
        
    chI = chainId[swapA];
    chJ = chainId[swapB];
    if (chainId[swapB] < chainId[swapA])
        {
        chI = chainId[swapB];
        chJ = chainId[swapA];
        }
    runId = chI / chainParams.numChains;
    chI = chI % chainParams.numChains;
    chJ = chJ % chainParams.numChains;
    swapInfo[runId][chJ][chI]++;
    if (isSwapSuccessful == YES)
        swapInfo[runId][chI][chJ]++;
#   endif
    
    return (NO_ERROR);
}


/* Autotune Dirichlet move */
void AutotuneDirichlet (MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *alphaPi, MrBFlt minTuning, MrBFlt maxTuning)
{
    MrBFlt delta, logTuning, newTuning;

    delta = 1.0 / sqrt(batch);
    delta = 0.01 < delta ? 0.01 : delta;

    logTuning = log(*alphaPi);

    if (acceptanceRate > targetRate)
        logTuning -= delta;
    else
        logTuning += delta;

    newTuning = exp(logTuning);
    if (newTuning > minTuning && newTuning < maxTuning)
        *alphaPi = newTuning;
}


/* Autotune multiplier move */
void AutotuneMultiplier (MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *lambda, MrBFlt minTuning, MrBFlt maxTuning)
{
    MrBFlt delta, logTuning, newTuning;

    delta = 1.0 / sqrt(batch);
    delta = 0.01 < delta ? 0.01 : delta;

    logTuning = log(*lambda);

    if (acceptanceRate > targetRate)
        logTuning += delta;
    else
        logTuning -= delta;

    newTuning = exp(logTuning);
    if (newTuning > minTuning && newTuning < maxTuning)
        *lambda = newTuning;
}


/* Autotune sliding window move */
void AutotuneSlider (MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *width, MrBFlt minTuning, MrBFlt maxTuning)
{
    MrBFlt delta, logTuning, newTuning;

    delta = 1.0 / sqrt(batch);
    delta = 0.01 < delta ? 0.01 : delta;

    logTuning = log(*width);

    if (acceptanceRate > targetRate)
        logTuning += delta;
    else
        logTuning -= delta;

    newTuning = exp(logTuning);
    if (newTuning > minTuning && newTuning < maxTuning)
        *width = newTuning;
}


void BuildExhaustiveSearchTree (Tree *t, int chain, int nTaxInTree, TreeInfo *tInfo)
{
    int         i;
    TreeNode    *p, *q, *r;
        
    if (nTaxInTree == t->nIntNodes + 1) {
        
        /* Get downpass */
        GetDownPass (t);

        /* Calculate cost of this tree and add to counter */
        tInfo->curScore = GetParsimonyLength (t, chain);
        if (tInfo->curScore < tInfo->minScore)
            {
            tInfo->totalScore *= pow ((tInfo->warp/3.0) / (1.0 - tInfo->warp), tInfo->minScore - tInfo->curScore);
            tInfo->totalScore += 1.0;
            tInfo->minScore = tInfo->curScore;
            }
        else
            tInfo->totalScore += pow (tInfo->warp/3.0, tInfo->curScore - tInfo->minScore) * pow (1.0-tInfo->warp, tInfo->minScore - tInfo->curScore);
    }

    else {

        /* find node to connect */
        q=tInfo->leaf[nTaxInTree];

        /* add using this ancestral node */
        p=tInfo->vertex[nTaxInTree-1];
        q->anc=p;
        p->right=q;

        for (i=0;i<2*nTaxInTree-1;i++) {
            /* find node to connect to */
            if (i>=nTaxInTree)
                r=tInfo->vertex[i-nTaxInTree];
            else
                r=tInfo->leaf[i];

            /* add to this node */
            p->left=r;
            if (r->anc==NULL)
                p->anc=NULL;
            else {
                p->anc=r->anc;
                if (r->anc->left==r)
                    r->anc->left=p;
                else
                    r->anc->right=p;
            }
            r->anc=p;

            /* next level */
            BuildExhaustiveSearchTree (t, chain, nTaxInTree+1, tInfo);

            if (tInfo->stopScore > 0.0 && tInfo->totalScore >= tInfo->stopScore)
                return;

            /* restore tree before trying next possibility */
            r->anc=p->anc;
            if (r->anc!=NULL) {
                if (r->anc->left==p)
                    r->anc->left=r;
                else
                    r->anc->right=r;
            }
        }
    }
}


/*------------------------------------------------------------------
|
|   BuildParsTrees: Fill in trees using random add seq with parsimony
|
------------------------------------------------------------------*/
int BuildParsTrees (RandLong *seed, int fromChain, int toChain)
{
    int         k, chn;
    Param       *p, *q;
    Tree        *tree;

    /* Build starting trees for state 0 */
    for (chn=fromChain; chn<toChain; chn++)
        {
        for (k=0; k<numParams; k++)
            {
            p = &params[k];
            if (p->paramType == P_TOPOLOGY)
                {
                assert (p->nSubParams == 1);
                q = p->subParams[0];
                tree = GetTree (q, chn, 0);
                /* fixed topology */
                if (p->paramId == TOPOLOGY_RCL_FIXED ||
                    p->paramId == TOPOLOGY_RCCL_FIXED ||
                    p->paramId == TOPOLOGY_CL_FIXED ||
                    p->paramId == TOPOLOGY_CCL_FIXED ||
                    p->paramId == TOPOLOGY_NCL_FIXED ||
                    p->paramId == TOPOLOGY_PARSIMONY_FIXED)
                    {
                    MrBayesPrint ("%s   Tree %s is fixed so a parsimony-based starting tree is not built\n", spacer, p->name);
                    return (NO_ERROR);
                    }
                /* constrained topology */
                else if (tree->nConstraints > 0)
                    {
                    MrBayesPrint ("%s   Tree %s is constrained and parsimony-based starting trees are not implemented for constrained trees (yet)\n", spacer, p->name);
                    return (NO_ERROR);
                    }
                /* random topology */
                else
                    {
                    if (BuildStepwiseTree (tree, chn, &globalSeed) == ERROR)
                        return (ERROR);
                    }
                if (InitializeTreeCalibrations (tree) == ERROR)
                    return (ERROR);
                FillTopologySubParams(p, chn, 0, seed);
                }
            }
        }

    return (NO_ERROR);
}


/* build (starting) topology stepwise */
int BuildStepwiseTree (Tree *t, int chain, RandLong *seed) {

    int         i, j, nTips;
    TreeNode    *p, *q, *r;
    CLFlt        length;

    // Allocate parsimony matrix if not done already

    /* get the tips */
    for (i=j=0; i<t->nNodes; i++) {
        p =  t->allDownPass[i];
        if (p->left == NULL && p->right == NULL)
            t->allDownPass[j++] = p;
        else if (p->right == NULL && p->anc == NULL && t->isRooted == NO)
            t->allDownPass[j++] = p;
    }
    nTips = j;

    /* order the tips randomly, use last as root */
    for (i=0; i<nTips-1; i++) {
        j = (int) (RandomNumber(seed)*(nTips-1-i));
        j += i;
        p = t->allDownPass[i];
        t->allDownPass[i] = t->allDownPass[j];
        t->allDownPass[j] = p;
    }

    /* build first tree */
    j = 0;
    q = t->allDownPass[0];
    r = t->allDownPass[1];
    p = t->intDownPass[j++];
    q->anc   = p;
    r->anc   = p;
    p->left  = q;
    p->right = r;
    q = t->allDownPass[nTips-1];
    q->anc   = NULL;
    q->right = NULL;
    q->left = p;
    p->anc = q;
    t->root = q;

    /* add nodes one at a time */
    for (i=2; i<nTips-1; i++) {
        r = t->allDownPass[i];
        p = t->intDownPass[j++];
        GetParsDP(t, t->root->left, chain);
        GetParsFP(t, t->root->left, chain);
        q = FindBestNode(t, t->root->left, r, &length, chain);
        p->left = q;
        p->right = r;
        p->anc = q->anc;
        if (q->anc->left == q)
            q->anc->left = p;
        else
            q->anc->right = p;
        q->anc = p;
        r->anc = p;
    }

    /* take care of the root */
    if (t->isRooted == YES) {
        r = t->root;
        q = t->allDownPass[t->nNodes-1];
        p = t->intDownPass[j];
        q->anc = q->right = NULL;
        q->left = p;
        p->anc = q;
        p->left = r->left;
        p->right = r;
        p->left->anc = p;
        r->left = r->right = NULL;
        r->anc = p;
        t->root = q;
    }
    GetDownPass (t);

    return (NO_ERROR);
}


/*------------------------------------------------------------------
|
|   CalcLikeAdgamma: calc likelihood for one adgamma correlation HMM
|
-------------------------------------------------------------------*/
int CalcLikeAdgamma (int d, Param *param, int chain, MrBFlt *lnL)
{
    int             c, i, j, nRates, posit, lastCharId;
    MrBFlt          logScaler, max, prob, *F,
                    *oldF, *tempF, fSpace[2][MAX_RATE_CATS];
    MrBFlt          *rP;
    CLFlt           freq, *lnScaler;
    ModelInfo       *m;
    ModelParams     *mp;
    BitsLong        *inHMM;
    
    /* find nRates for first division in HMM */
    m = &modelSettings[d];
    mp = &modelParams[d];
    nRates = m->numRateCats;

    /* calculate rate category frequencies */
    freq = (CLFlt) ((CLFlt) 1.0 / nRates);

    /* find Markov trans probs */
    F = GetParamSubVals (param,chain, state[chain]);
    for (i=posit=0; i<nRates; i++)
        for (j=0; j<nRates; j++)
            markovTi[0][i][j] = F[posit++];
    
    /* precalculate Markov trans probs up to largest small jump */
    /* but only if needed                                       */
    for (i=1; i<MAX_SMALL_JUMP; i++)
        {
        if (hasMarkovTi[i] == YES)
            {
            if (hasMarkovTi[i-1] == YES || i == 1)
                MultiplyMatrices(nRates, markovTi[i-1], markovTi[0], markovTi[i]);
            else
                MultiplyMatrixNTimes(nRates, markovTi[0], i+1, markovTi[i]);
            }
        }
        
    /* find site scaler for this chain and state */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find rate probs for this chain and state */
    rP = rateProbs[chain] + state[chain] * rateProbRowSize;

    /* set bit vector indicating divisions in this HMM */
    inHMM = (BitsLong *) SafeCalloc (((param->nRelParts/nBitsInALong) + 1), sizeof(BitsLong));
    for (i=0; i<param->nRelParts; i++)
        {
        if (modelSettings[param->relParts[i]].shape ==
            modelSettings[d].shape)
            {
            SetBit(param->relParts[i], inHMM);
            }
        }
    
    /* Perform the so-called forward algorithm of HMMs  */  
    /* set up the space for f(c,i) */
    F = fSpace[0];
    oldF = fSpace[1];

    for (c=0; c<numChar; c++)
        {
        if (IsBitSet(partitionId[c][partitionNum] - 1, inHMM) == YES)
            break;
        }

    /* fill in fi(0) */
    max = 0.0;
    m = &modelSettings[partitionId[c][partitionNum] - 1];
    posit = m->rateProbStart + (compCharPos[c] - m->compCharStart) * m->numRateCats;

    for (i=0; i<nRates; i++)
        {
        F[i] = rP[posit++];
        if (F[i] > max)
            max = F[i];
        }

    for (i=0; i<nRates; i++)
        F[i] /= max;

    /* set logscaler to the value from the first character */
    logScaler = lnScaler[compCharPos[c]] +  log(max);

    /* now step along the sequence to the end */
    lastCharId = charInfo[c].charId;
    for (c++; c<numChar; c++)
        {
        /* skip if excluded */
        if (charInfo[c].isExcluded == YES)
            continue;

        /* skip if part of same codon in translated protein model */
        if ((mp->dataType == DNA || mp->dataType == RNA) && m->dataType == PROTEIN && charInfo[c].charId == lastCharId)
            continue;
        else
            lastCharId = charInfo[c].charId;
        
        /* skip if not in HMM */
        if (IsBitSet(partitionId[c][partitionNum] - 1, inHMM) == NO)
            continue;
        
        /* switch F and oldF, since the previous F is now old */
        tempF = F;
        F = oldF;
        oldF = tempF;

        /* find the position of the rate probs */
        m = &modelSettings[partitionId[c][partitionNum] - 1];
        posit = m->rateProbStart + (compCharPos[c] - m->compCharStart) * m->numRateCats;
                
        /* calculate the HMM forward probs fi(x) at site x in HMM */
        if (siteJump[c] <= MAX_SMALL_JUMP)
            {
            max = 0.0;
            for (i=0; i<nRates; i++)
                {
                prob = 0.0;
                for (j=0; j<nRates; j++)
                    prob += markovTi[siteJump[c]-1][i][j] * oldF[j];
                F[i] = rP[posit++] * prob;
                if (F[i] > max)
                    max = F[i];
                }
            }
        else if (siteJump[c] < BIG_JUMP)    /* intermediate jump, calculate trans probs */
            {
            MultiplyMatrixNTimes(nRates, markovTi[0], siteJump[c], markovTiN);
            max = 0.0;
            for (i=0; i<nRates; i++)
                {
                prob = 0.0;
                for (j=0; j<nRates; j++)
                    prob += markovTiN[i][j] * oldF[j];
                F[i] = rP[posit++] * prob;
                if (F[i] > max)
                    max = F[i];
                }
            }
        else    /* big jump, use stationary freqs */
            {
            max = 0.0;
            for (i=0; i<nRates; i++)
                {
                prob = 0.0;
                for (j=0; j<nRates; j++)
                    prob += (oldF[j] / freq);
                F[i] = rP[posit++] * prob;
                if (F[i] > max)
                    max = F[i];
                }
            }

        /* rescale and adjust total scaler with HMM scaler and site scaler */
        for (i=0; i<nRates; i++)
            F[i] /= max;

        logScaler += lnScaler[compCharPos[c]] +  log(max);

        }
    
    /* now pull the rate probs together at the end, F contains the vals needed */
    prob =  0.0;
    for (i=0; i<nRates; i++)
        prob += (freq * F[i]);

    (*lnL) = logScaler +  log(prob);

    free (inHMM);

    return (NO_ERROR);
}


/* CalcPartFreqStats: Calculate standard deviation of partition frequencies */
void CalcPartFreqStats (PFNODE *p, STATS *stat)
{
    int     i, j, n, min;
    MrBFlt  f, sum, sumsq, stdev;

    n = chainParams.numRuns;
    min = (int)(chainParams.minPartFreq * stat->numSamples);
    if ((MrBFlt)min != chainParams.minPartFreq * stat->numSamples)
        min++;

    /* recursively compute partition frequencies for all subpartitions */
    if (p->left != NULL) 
        CalcPartFreqStats (p->left, stat);
    if (p->right != NULL)
        CalcPartFreqStats (p->right, stat);

    for (i=0; i<n; i++)
        {
        if (p->count[i] >= min)
            break;
        }

    if (i == n)
        return;

    sum = 0.0;
    sumsq = 0.0;
    for (i=0; i<n; i++)
        {
        f = (MrBFlt) (p->count[i]) / (MrBFlt) (stat->numSamples);
        sum += f;
        sumsq += f * f;
        }
    
    f = (sumsq - sum * sum / n) / (n - 1);
    if (f < 0.0)
        stdev = 0.0;
    else
        stdev = sqrt (f);
    
    stat->sum += stdev;
    if (stdev > stat->max)
        stat->max = stdev;

    stat->numPartitions++;

    if (chainParams.allComps == YES)
        {
        for (i=0; i<n; i++)
            {
            for (j=i+1; j<n; j++)
                {
                if (p->count[i] < min && p->count[j] < min)
                    continue;

                sum = 0.0;
                sumsq = 0.0;

                f = (MrBFlt) (p->count[i]) / (MrBFlt) (stat->numSamples);
                sum += f;
                sumsq += f * f;
                
                f = (MrBFlt) (p->count[j]) / (MrBFlt) (stat->numSamples);
                sum += f;
                sumsq += f * f;

                f = (sumsq - sum * sum / 2.0);
                if (f < 0.0)
                    stdev = 0.0;
                else
                    stdev = sqrt (f);
                
                if (chainParams.diagnStat == AVGSTDDEV)
                    stat->pair[i][j] += stdev;
                else if (stdev > stat->pair[i][j])
                    stat->pair[i][j] = stdev;
                stat->pair[j][i]++;
                }
            }
        }
}


/*----------------------------------------------------------------
|
|   CalcTopoConvDiagn: Calculate average and max standard deviation
|                 in clade credibility (partition frequency) values
|
----------------------------------------------------------------*/
void CalcTopoConvDiagn (int numSamples)
{
    int     i, j, n;
    STATS   *stat;
    
    for (n=0; n<numTopologies; n++)
        {
        stat = &chainParams.stat[n];
        stat->numSamples = numSamples;
        stat->numPartitions = 0.0;
        stat->sum = 0.0;
        stat->max = 0.0;

        if (chainParams.allComps == YES)
            {
            for (i=0; i<chainParams.numRuns; i++)
                for (j=0; j<chainParams.numRuns; j++)
                    stat->pair[i][j] = 0.0;
            }
    
        CalcPartFreqStats (partFreqTreeRoot[n], stat);
        
        stat->avgStdDev = stat->sum / stat->numPartitions;
        }
}


/* used in DoCompRefTree */
void PartFreq (PFNODE *p, STATS *stat, int *ntrees)
{
    int     i, n = chainParams.numRuns;
    MrBFlt  f, sum, sumsq, stdev;
    
    /* recursively compute partition frequencies for all subpartitions */
    if (p->left != NULL)
        PartFreq (p->left, stat, ntrees);
    if (p->right != NULL)
        PartFreq (p->right, stat, ntrees);
    
    sum = sumsq = 0.0;
    for (i=0; i<chainParams.numRuns; i++)
        {
        f = (MrBFlt)(p->count[i]) / (MrBFlt)ntrees[i];
        sum += f;
        sumsq += f * f;
        }
    
    f = (sumsq - sum * sum / n) / (n - 1);
    if (f < 0.0)
        stdev = 0.0;
    else
        stdev = sqrt (f);
    
    stat->sum += stdev;
    if (stat->max < stdev)
        stat->max = stdev;
    
    stat->numPartitions++;
}
void CalcTopoConvDiagn2 (int *nTrees)
{
    int     n;
    STATS   *stat;
    
    for (n=0; n<numTopologies; n++)
        {
        stat = &chainParams.stat[n];
        stat->numPartitions = 0.0;
        stat->sum = stat->max = 0.0;
    
        PartFreq (partFreqTreeRoot[n], stat, nTrees);
        
        stat->avgStdDev = stat->sum / stat->numPartitions;
        }
}


int CheckTemperature (void)
{
    if (chainParams.userDefinedTemps == YES)
            {
          if (AreDoublesEqual(chainParams.userTemps[0], 1.0, ETA)==NO)
            {
            MrBayesPrint ("%s   The first user-defined temperature must be 1.0.\n", spacer);
            return (ERROR);
            }
        }

    return (NO_ERROR);
}


void CloseMBPrintFiles (void)
{
    int     i, n;

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return;
#   endif

    for (n=0; n<chainParams.numRuns; n++)
        {
        SafeFclose (&fpParm[n]);
#if defined (PRINT_DUMP)
        SafeFclose (&fpDump[n]);
#endif

        for (i=0; i<numTrees; i++)
            {
            if (fpTree[n][i])
                {
                fprintf (fpTree[n][i], "end;\n");
                SafeFclose (&fpTree[n][i]);
                }
            }
        }

    if (chainParams.mcmcDiagn == YES)
        SafeFclose (&fpMcmc);

    if (chainParams.isSS == YES)
        SafeFclose (&fpSS);
}


/* CompactTree: prune partition tree */
PFNODE *CompactTree (PFNODE *p)
{
    int         i, j;
    PFNODE      *q, *r;

    if (p == NULL)
        return NULL;
    
    i = j = 0;
    if (IsPFNodeEmpty(p) == YES)
        {
        /* steal info from terminal on the way up */
        q = SmallestNonemptyPFNode (p->left, &i, 0);
        r = LargestNonemptyPFNode (p->right, &j, 0);

        if (q != NULL || r != NULL)
            {
            if (i < j)
                q = r;
        
            for (i=0; i<chainParams.numRuns; i++)
                {
                p->count[i] = q->count[i];
                q->count[i] = 0;
                }
            for (i=0; i<nLongsNeeded; i++)
                p->partition[i] = q->partition[i];
            }
        }

    p->left = CompactTree (p->left);
    p->right = CompactTree (p->right);

    /* delete on the way down if empty */
    if (IsPFNodeEmpty(p) == YES)
        {
        Tfree (p);
        return NULL;
        }
    else
        return p;
}


/*-----------------------------------------------------------------
|
|   CopyParams: copy parameters of touched divisions
|
-----------------------------------------------------------------*/
void CopyParams (int chain)
{
    int         i, j, k, fromState, toState, *fromInt, *toInt;
    MrBFlt      *from, *to;
    ModelInfo   *m;
    Param       *p;

    /* copy all params                                               */
    /* now done for all vars, can also be done for only touched vars */
    /* but then m->upDateCl must be kept separate for each chain!    */
    for (i=0; i<numParams; i++)
        {
        p = &params[i];

        from = GetParamVals (p, chain, state[chain]);
        to = GetParamVals (p, chain, (state[chain] ^ 1));

        for (j=0; j<p->nValues; j++)
            to[j] = from[j];

        from = GetParamSubVals (p, chain, state[chain]);
        to = GetParamSubVals (p, chain, (state[chain] ^ 1));

        for (j=0; j<p->nSubValues; j++)
            to[j] = from[j];

        fromInt = GetParamIntVals (p, chain, state[chain]);
        toInt = GetParamIntVals (p, chain, (state[chain] ^ 1));

        for (j=0; j<p->nIntValues; j++)
            toInt[j] = fromInt[j];

        if (p->nStdStateFreqs > 0)
            {
            from = GetParamStdStateFreqs (p, chain, state[chain]);
            to = GetParamStdStateFreqs (p, chain, state[chain] ^ 1);
            for (j=0; j<p->nStdStateFreqs; j++)
                to[j] = from[j];
            }

        if (p->paramType == P_CPPEVENTS)
            {
            fromState = 2*chain + state[chain];
            toState   = 2*chain + (state[chain] ^ 1);
            for (j=0; j<2*numLocalTaxa-2; j++)
                {
                if (p->nEvents[toState][j] != p->nEvents[fromState][j])
                    {
                    if (p->nEvents[fromState][j] == 0)
                        {
                        free (p->position[toState][j]);
                        p->position[toState][j] = NULL;
                        free (p->rateMult[toState][j]);
                        p->rateMult[toState][j] = NULL;
                        }
                    else if (p->nEvents[toState][j] == 0)
                        {
                        p->position[toState][j] = (MrBFlt *) SafeCalloc (p->nEvents[fromState][j], sizeof (MrBFlt));
                        p->rateMult[toState][j] = (MrBFlt *) SafeCalloc (p->nEvents[fromState][j], sizeof (MrBFlt));
                        }
                    else
                        {
                        p->position[toState][j] = (MrBFlt *) SafeRealloc ((void *)p->position[toState][j], p->nEvents[fromState][j] * sizeof (MrBFlt));
                        p->rateMult[toState][j] = (MrBFlt *) SafeRealloc ((void *)p->rateMult[toState][j], p->nEvents[fromState][j] * sizeof (MrBFlt));
                        }
                    p->nEvents[toState][j] = p->nEvents[fromState][j];
                    }
                if (p->nEvents[fromState][j] > 0)
                    {
                    for (k=0; k<p->nEvents[fromState][j]; k++)
                        {
                        p->position[toState][j][k] = p->position[fromState][j][k];
                        p->rateMult[toState][j][k] = p->rateMult[fromState][j][k];
                        }
                    }
                }
            }
        }

    /* copy division params (model settings) for chain */
    /* reset division update flags                     */
    fromState = 2 * chain + state[chain];
    toState   = 2 * chain + (state[chain] ^ 1);
    for (i=0; i<numCurrentDivisions; i++)
        {
        m = &modelSettings[i];
        m->lnLike[toState] = m->lnLike[fromState];
        if (m->parsModelId == YES)
            m->parsTreeLength[toState] = m->parsTreeLength[fromState];
        m->upDateCl = NO;
        m->upDateCijk = NO;
        m->upDateAll = NO;
        }
        
    return;
}


/* CopySiteScalers: Copy site scalers from scratch space into current space */
void CopySiteScalers (ModelInfo *m, int chain)
{
    CLFlt       *from, *to;
#   if defined (BEAGLE_ENABLED)
    int         i, j;
#   endif

#   if defined (BEAGLE_ENABLED)
    if (m->useBeagle == YES)
        {
        j = m->siteScalerScratchIndex;
        for (i=0; i<m->nCijkParts; i++)
            {
            beagleResetScaleFactors (m->beagleInstance,
                                     m->siteScalerIndex[chain] + i);
            beagleAccumulateScaleFactors (m->beagleInstance,
                                          &j,
                                          1,
                                          m->siteScalerIndex[chain] + i);
            j++;
            }
        return;
        }
#   endif
    from = m->scalers[m->siteScalerScratchIndex];
    to   = m->scalers[m->siteScalerIndex[chain]];
    memcpy ((void*) to, (void*) from, (size_t)(m->numChars) * sizeof(CLFlt));
}


/*-----------------------------------------------------------------
|
|   CopyTrees: copies touched trees for chain
|       resets division and node update flags in the process
|       Note: partition information of nodes are not copied if
|       either source or destination tree does not have bitsets allocated
|
-----------------------------------------------------------------*/
void CopyTrees (int chain)
{
    int         i, j, n, nTaxa, nLongsNeeded;
    TreeNode    *p, *q;
    Tree        *from, *to;

    /* reset division update flags */
    for (i=0; i<numCurrentDivisions; i++)
        modelSettings[i].upDateCl = NO;

    for (n=0; n<numTrees; n++)
        {
        from = GetTreeFromIndex (n, chain, state[chain]);       
        to = GetTreeFromIndex (n, chain, (state[chain] ^ 1));
        if (from->bitsets != NULL && to->bitsets != NULL)
            {
            if (from->isRooted == NO)
                nTaxa = from->nNodes - from->nIntNodes;
            else
                nTaxa = from->nNodes - from->nIntNodes - 1;
            nLongsNeeded = (int)((nTaxa - 1) / nBitsInALong) + 1;
            }
        else
            nLongsNeeded = 0;

        /* copy nodes */
        for (j=0; j<from->nNodes; j++)
            {
            /* copy pointers */
            p  = from->nodes + j;
            q  = to->nodes + j;

            if (p->anc != NULL)
                q->anc = to->nodes + p->anc->memoryIndex;
            else
                q->anc = NULL;

            if (p->left != NULL)    
                q->left = to->nodes + p->left->memoryIndex;
            else
                q->left = NULL;

            if (p->right != NULL)   
                q->right = to->nodes + p->right->memoryIndex;
            else
                q->right = NULL;

            CopyTreeNodes (q, p, nLongsNeeded);
            q->upDateCl = q->upDateTi = NO;     /* reset update flags */
            }
        
        for (i=0; i<from->nIntNodes; i++)
            {
            to->intDownPass[i] = to->nodes + from->intDownPass[i]->memoryIndex;
            }
        for (i=0; i<from->nNodes; i++)
            {
            to->allDownPass[i] = to->nodes + from->allDownPass[i]->memoryIndex;
            }

        to->root = to->nodes + from->root->memoryIndex;

        /* rest of tree info is constant and need not be copied */
        }
    
    return;
}


#ifdef VISUAL
BOOL WINAPI CatchInterrupt (DWORD signum)
{
    /* set up signal handler to do the same */
    MrBayesPrint ("\n   Ctrl-C detected\n");
    requestAbortRun = YES;
    return TRUE;
}
#else
void CatchInterrupt (int signum)
{
    /* set up signal handler to do the same */
    signal (signum, CatchInterrupt);
    requestAbortRun = YES;
    MrBayesPrint ("\n   Ctrl-C detected\n");
}
#endif


/*----------------------------------------------------------------
|
|   DebugNodeScalers: Calculate node scalers sum
|
-----------------------------------------------------------------*/
CLFlt DebugNodeScalers (TreeNode *p, int division, int chain)
{
    int             c;
    CLFlt           *scP;
    CLFlt           sum=0.0;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* find scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* remove scalers */
    for (c=0; c<m->numChars; c++)
        sum += scP[c];

    return sum;
}


/*----------------------------------------------------------------
|
|   DebugTreeScalers: Calculate DebugNodeScalers for each node and printit
|
-----------------------------------------------------------------*/
void DebugTreeScalers(int chain, int d)
{
    int i;
    TreeNode        *p;
    ModelInfo       *m;
    Tree            *tree;
    
    m = &modelSettings[d];
    tree = GetTree(m->brlens, chain, state[chain]);
    
    if (m->parsModelId == NO)
        {
        for (i=0; i<tree->nIntNodes; i++)
            {
            p = tree->intDownPass[i];
            
            if (m->unscaledNodes[chain][p->index] == 0)
                {
                printf ("Node:%d Sum scalers:%f\n",p->index,DebugNodeScalers(p, d, chain));
                }
            }
        }
}


int DoMcmc (void)
{
    RandLong    seed;
    int         rc, i, j, run;
    char        c;
    FILE        *tempFile;
    char        temp[20];
    char        *strBuf,*tmpcp;
    double      tmp;

#   if defined (BEST_MPI_ENABLED)
    Tree        *tree;
#   endif

#   if !defined (VISUAL) && !defined (MPI_ENABLED)
    sighandler_t sigint_oldhandler, sigterm_oldhandler;
#   endif

    numPreviousGen = 0;     /* Make sure this is reset */

    /* Check to see that we have a data matrix. Otherwise, the MCMC is rather
       pointless. */
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        goto errorExit;
        }

    if (setUpAnalysisSuccess == NO)
        {
        MrBayesPrint ("%s   The analysis could not be started because there was an error during its setup.\n", spacer);
        MrBayesPrint ("%s   Refer to error messages printed during model set up to adress the problem.\n", spacer);
        goto errorExit;
        }

    /* set file names */
    sumtParams.numRuns  = chainParams.numRuns;
    sumpParams.numRuns  = chainParams.numRuns;
    sumssParams.numRuns = chainParams.numRuns;
    
    if (fileNameChanged == YES)
        {
        SetFileNames();
        fileNameChanged = NO;
        }

    MrBayesPrint ("%s   Running Markov chain\n", spacer);
    
    /* Check the chain temperature parameters */
    if (CheckTemperature () == ERROR)
        goto errorExit;
    
    /* Set the chain random number seeds here. We have two seeds. One
       (called swapSeed) is only used to determine which two chains 
       will swap states in the next trial. The other (called seed) is
       the seed for our work-horse pseudorandom number generator. Note
       that if we are doing MPI, we want the swap seed to be the same
       for every processor. This is taken care of when we initialize
       things in the program. If we are doing MPI, we also want to make
       certain that seed is different for every processor. */
#   if defined (MPI_ENABLED)
    seed = globalSeed + (proc_id + 1);
    if (seed < 0)
        seed = -seed;
#   else
    seed = globalSeed;
#   endif

    /* Get a unique identifier (stamp) for this run. This is used as
       an identifier for each mcmc analysis. It uses runIDSeed to initialize 
       the stamp. All of the processors should have the same seed, so 
       this should be safe. */
    GetStamp ();
    
    MrBayesPrint ("%s   Seed = %d\n", spacer, seed);
    MrBayesPrint ("%s   Swapseed = %d\n", spacer, swapSeed);

    /* Show the model to make sure the user sees it before running the analysis */
    if (ShowModel() == ERROR)
        goto errorExit;
    MrBayesPrint ("\n");

    /* Warn the user or stop analysis in case the model is strange */
    if (CheckModel() == ERROR)
        goto errorExit;
                
    /* Determine the number of local chains and data splits */
    if (SetLocalChainsAndDataSplits() == ERROR)
        goto errorExit;

    /* Set up the moves to be used */
    if (SetUsedMoves () == ERROR)
        goto errorExit;

    /* Check to see that we have at least one move. Otherwise, the MCMC is rather
       pointless. */
    if (numUsedMoves == 0)
        {
        MrBayesPrint ("%s   No move is currently switched on.\n", spacer);
        MrBayesPrint ("%s   There must be at least one move to run an MCMC analysis.\n", spacer);
        MrBayesPrint ("%s   Switch on moves using the 'propset' command.\n", spacer);
        goto errorExit;
        }

    /* Show summary table of moves that will be used */ 
    if (ShowMoveSummary () == ERROR)
        goto errorExit;

    /* Set the likelihood function pointers. */
    if (SetLikeFunctions () == ERROR)
        goto errorExit;

    /* Set up number of characters of each character pattern. */
    if (FillNumSitesOfPat () == ERROR)
        goto errorExit;

    /* Initialize parsimony sets. */
    if (InitParsSets() == ERROR)
        goto errorExit;

    /* Set up a terminal state index matrix for local compression. */
    if (SetUpTermState() == ERROR)
        goto errorExit;

    /* Initialize conditional likelihoods and transition probabilities for chain (the working space). */
    if (InitChainCondLikes () == ERROR)
        goto errorExit;

    /* Initialize adgamma conditional likelihoods */
    if (InitAdGamma () == ERROR)
        goto errorExit;
    
    /* Initialize invariable conditional likelihoods. */
    if (InitInvCondLikes() == ERROR)
        goto errorExit;

    /* Allocate BEST chain variables */
    if (numTopologies > 1 && !strcmp(modelParams[0].topologyPr,"Speciestree"))
        AllocateBestChainVariables();

    /* allocate SS memory for the chains if needed */
    if (chainParams.isSS == YES)
        {
        if (memAllocs[ALLOC_SS] == YES)
            {
            MrBayesPrint ("%s   SS is already allocated\n", spacer);
            goto errorExit;
            }
        else if ((marginalLnLSS = (MrBFlt *) SafeCalloc (chainParams.numRuns, sizeof(MrBFlt))) == NULL)
            {
            MrBayesPrint ("%s   Problem allocating marginalLnLSS\n", spacer);
            goto errorExit;
            }
        else if ((stepScalerSS = (MrBFlt *) SafeCalloc (chainParams.numRuns, sizeof(MrBFlt))) == NULL)
            {
            MrBayesPrint ("%s   Problem allocating stepScalerSS\n", spacer);
            free (marginalLnLSS);
            goto errorExit;
            }
        else if ((stepAcumulatorSS = (MrBFlt *) SafeCalloc (chainParams.numRuns, sizeof(MrBFlt))) == NULL)
            {
            MrBayesPrint ("%s   Problem allocating stepAcumulatorSS\n", spacer);
            free (stepScalerSS);
            free (marginalLnLSS);
            goto errorExit;
            }
        else if ((splitfreqSS = (MrBFlt *) SafeCalloc (chainParams.numStepsSS*numTopologies, sizeof(MrBFlt))) == NULL)
            {
            MrBayesPrint ("%s   Problem allocating splitfreqSS\n", spacer);
            free (stepScalerSS);
            free (marginalLnLSS);
            free (stepAcumulatorSS);
            goto errorExit;
            }
        else
            memAllocs[ALLOC_SS] = YES;
        }

    /* Either append to previous run or deal with starting values */
    if (chainParams.append == YES)
        {
        /* Continue old run */
       
        /* Get starting values from checkpoint file */
        MrBayesPrint ("%s   Getting values from previous run\n", spacer);
        strcpy(inputFileName, chainParams.chainFileName);
        strcat(inputFileName, ".ckp");
        if (OpenTextFileR(inputFileName) == NULL)
            {
            MrBayesPrint ("%s   Could not find the checkpoint file '%s'.\n", spacer, inputFileName);
            MrBayesPrint ("%s   Make sure it is in the working directory.\n", spacer);
            goto errorExit;
            }

        if (DoExecute () == ERROR)
            goto errorExit;

        /* Get number of generations to start from and SS information if needed */
        temp[0] = '\0';
        numPreviousGen = 0;
#   if defined (MPI_ENABLED)
        if (proc_id == 0) {
#   endif
        tempFile = OpenBinaryFileR (inputFileName);
        do { c = fgetc(tempFile);
            } while (c!=':' && c!=EOF);
        if (c!=EOF)
            {
            do { c = fgetc(tempFile);
                } while (c!=':' && c!=EOF);
            }
        if (c!=EOF)
            {
            do { c = fgetc(tempFile);
                } while (!isdigit(c) && c!=EOF);
            }
        if (c!=EOF)
            {
            i=0;
            while (c >= '0' && c <= '9' && i < 18)
                {
                temp[i++] = c;
                c = fgetc(tempFile);
                }
            temp[i] = '\0';
            numPreviousGen = atoi(temp);
            }
        if (chainParams.isSS==YES && c!=EOF)
            {
            do { c = fgetc(tempFile);
                } while (c!=':' && c!=EOF);
            strBuf = (char *) SafeCalloc (chainParams.numRuns*20,sizeof(char));
            if (fgets(strBuf,chainParams.numRuns*20,tempFile)==NULL)
                {
                MrBayesPrint ("%s   Error: Reading SsAcumulators from .ckp file fails.\n", spacer);
                free(strBuf);
                goto errorExit;
                }

            tmpcp=strtok(strBuf," "); 
            for (run=0; run<chainParams.numRuns; run++)
                {
                if (tmpcp == NULL)
                    {
                    MrBayesPrint ("%s   Error: Not enough values in SsAcumulators comment of .ckp file.   \n", spacer);
                    free(strBuf);
                    goto errorExit;
                    }
                    
                tmp=atof(tmpcp);
                stepScalerSS[run]=tmp-10;
                stepAcumulatorSS[run]=exp(10);
                tmpcp=strtok(NULL," ]");
                }
                      
            free(strBuf);
            }
#   if defined (MPI_ENABLED)
        }
        MPI_Bcast (&numPreviousGen, 1, MPI_INT, 0, MPI_COMM_WORLD);
#   endif
        if (numPreviousGen == 0)
            {
            MrBayesPrint ("%s   Could not find the number of generations in previous run.\n", spacer);
            goto errorExit;
            }
        else if (numPreviousGen >= chainParams.numGen)
            {
            MrBayesPrint ("%s   The specified number of generations (%d) was already finished in\n", spacer, chainParams.numGen);
            MrBayesPrint ("%s   the previous run you are trying to append to.\n", spacer);
            goto errorExit;
            }
        else
            MrBayesPrint ("%s   Using samples up to generation %d from previous analysis.\n", spacer, numPreviousGen);
        }
    else
        {
        /* New run */
        
        /* deal with starting param values */
        if (!strcmp(chainParams.startParams,"Reset"))
            {
            MrBayesPrint ("%s   Resetting starting values for substitution model parameters\n", spacer);
            FillNormalParams (&seed, 0, numLocalChains);
            }

        /* deal with starting treeparam values */
        if (!strcmp(chainParams.startTree,"Random"))
            {
            MrBayesPrint ("%s   Resetting starting trees and tree parameters\n", spacer);
            FillTreeParams (&seed, 0, numLocalChains);
            }
        else if (!strcmp(chainParams.startTree,"Parsimony"))
            {
            MrBayesPrint ("%s   Rebuilding starting trees using random addition sequences and parsimony\n", spacer);
            BuildParsTrees (&seed, 0, numLocalChains);
            }

        /* Perturb start trees if requested */
        if (chainParams.numStartPerts > 0)
            {
            MrBayesPrint ("%s   Randomly perturbing starting trees\n", spacer);
            for (i=0; i<numTrees; i++)
                {
                for (j=0; j<numGlobalChains; j++)
                    RandPerturb (GetTreeFromIndex(i, j, 0), chainParams.numStartPerts, &seed);
                }
            }
        }

    /* Set clockRate if we have calibration */
    for (j=0; j<numGlobalChains; j++)
        {
        if (UpdateClockRate(0.0, j) == ERROR) 
            goto errorExit;
        }
    /*
    for (i=0; i<numParams; i++)
        {
        for (j=0; j<numGlobalChains; j++)
            assert (IsTreeConsistent(&params[i], j, 0) == YES);
        }  */

    /* Initialize vectors of print parameters */
    if (InitPrintParams () == ERROR)
        goto errorExit;
    
    /*! setup a signal handler to catch interrupts, ignore failure */
#   ifdef VISUAL
    SetConsoleCtrlHandler (CatchInterrupt, TRUE);
#   else
#       if !defined (MPI_ENABLED)
    /* we do not want to mess with the signal handling in MPI version */
    sigint_oldhandler  = signal(SIGINT, CatchInterrupt);
    sigterm_oldhandler = signal(SIGTERM, CatchInterrupt);
#       endif
#   endif
    requestAbortRun = NO;

    /* Run the Markov chain. */
    rc = RunChain (&seed);
    if (rc == ERROR)
        {
#   ifdef VISUAL
        SetConsoleCtrlHandler (CatchInterrupt, FALSE);
#   else
#       if !defined (MPI_ENABLED)
        signal(SIGINT, sigint_oldhandler);
        signal(SIGTERM, sigterm_oldhandler);
#       endif
#   endif
        goto errorExit;
        }
    else if (rc == ABORT)
        {
        ResetChainIds();
        FreeChainMemory();
#   ifdef VISUAL
        SetConsoleCtrlHandler (CatchInterrupt, FALSE);
#   else
#       if !defined (MPI_ENABLED)
        signal(SIGINT, sigint_oldhandler);
        signal(SIGTERM, sigterm_oldhandler);
#       endif
#   endif
        return ABORT;
        }
        
    /*! restore the default signal handler */
#   ifdef VISUAL
    SetConsoleCtrlHandler (CatchInterrupt, FALSE);
#   else
#       if !defined (MPI_ENABLED)
    signal(SIGINT, sigint_oldhandler);
    signal(SIGTERM, sigterm_oldhandler);
#       endif
#   endif

    /* Reset the global seed at end of chain. We don't want successive
       chains to all start with the same random number seed. */
    globalSeed = seed;

    /* Free up all memory allocated for the chain. */
    FreeChainMemory ();
    
    return (NO_ERROR);
    
    errorExit:
        FreeChainMemory ();
        return (ERROR);
}


int DoMcmcp (void)
{
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }

    sumtParams.numRuns = chainParams.numRuns;
    sumpParams.numRuns = chainParams.numRuns;
    
    if (fileNameChanged == YES)
        {
        SetFileNames();
        fileNameChanged = NO;
        }

    MrBayesPrint ("%s   Successfully set chain parameters\n", spacer);

    return (NO_ERROR);
}


int DoSsParm (char *parmName, char *tkn)
{
    int         tempI;
    MrBFlt      tempD;
    char        tempStr[5];
    static int  negBurninss;

    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        if (!strcmp(parmName, "Burninss"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                negBurninss = NO;
                expecting = Expecting(NUMBER) | Expecting(DASH);
                }
            else if (expecting == Expecting(DASH))
                {
                negBurninss = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (negBurninss == NO)
                    chainParams.burninSS = tempI;
                else
                    chainParams.burninSS = -tempI;
                MrBayesPrint ("%s   Setting burnin for stepping-stone sampling to %d\n", spacer, chainParams.burninSS);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                return (ERROR);
                }         
            }
        else if (!strcmp(parmName, "Nsteps"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                chainParams.numStepsSS = tempI;
                MrBayesPrint ("%s   Setting number of steps in stepping-stone sampling to %d\n", spacer, chainParams.numStepsSS);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                return (ERROR);
                }                
            }
        else if (!strcmp(parmName, "FromPrior"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.startFromPriorSS = YES;
                    else
                        chainParams.startFromPriorSS = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for FromPrior parameter\n", spacer);
                    return (ERROR);
                    }
                MrBayesPrint ("%s   Setting FromPrior=%s\n", spacer, (chainParams.startFromPriorSS==YES)?"Yes":"No");
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                return (ERROR);
                }                
            }
        else if (!strcmp(parmName, "Alpha"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                chainParams.alphaSS = tempD;
                MrBayesPrint ("%s   Setting alpha in stepping-stone sampling to %lf\n", spacer, chainParams.alphaSS);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                return (ERROR);
                }                
            }
        else
            {
            return (ERROR);
            }
        }
    return (NO_ERROR);
}


int DoMcmcParm (char *parmName, char *tkn)
{
    int         tempI;
    MrBFlt      tempD;
    char        *tempStr;
    int         tempStrSize = TEMPSTRSIZE;
    
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }
    *tempStr='\0';

    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Seed (globalSeed) ***************************************************************/
        if (!strcmp(parmName, "Seed"))
            {
                MrBayesPrint ("%s   Error: Setting \"Seed\" in mcmc command is depricated. Use \"set\" command instead.\n", spacer);
                MrBayesPrint ("%s   For more information type \"help set\";\n", spacer);
                free (tempStr);
                return (ERROR);
            /*
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                globalSeed = tempI;
                MrBayesPrint ("%s   Setting seed to %ld\n", spacer, globalSeed);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                free (tempStr);
                return (ERROR);
                }
                */
            }
        /* set Swapseed (global variable swapSeed) ***************************************************************/
        else if (!strcmp(parmName, "Swapseed"))
            {
                MrBayesPrint ("%s   Error: Setting \"Swapseed\" in mcmc command is depricated. Use \"set\" command instead.\n", spacer);
                MrBayesPrint ("%s   For more information type \"help set\";\n", spacer);
                free (tempStr);
                return (ERROR);
                /*
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                swapSeed = tempI;
                MrBayesPrint ("%s   Setting swapseed to %ld\n", spacer, swapSeed);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free (tempStr);
                return (ERROR);
                }
                */
            }
        /* set run ID */
        /* this setting is provided for GRID use only, so that identical runs can be generated */
        else if (!strcmp(parmName, "Runidseed"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                runIDSeed = tempI;
                MrBayesPrint ("%s   Setting run ID [stamp] seed to %ld [for GRID use]\n", spacer, runIDSeed);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free (tempStr);
                return (ERROR);
                }
            }
        /* set Ngen (numGen) ******************************************************************/
        else if (!strcmp(parmName, "Ngen"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   Too few generations\n", spacer);
                    return (ERROR);
                    }
                chainParams.numGen = tempI;
                MrBayesPrint ("%s   Setting number of generations to %d\n", spacer, chainParams.numGen);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free (tempStr);
                return (ERROR);
                }
            }
        /* set Samplefreq (sampleFreq) ********************************************************/
        else if (!strcmp(parmName, "Samplefreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   Sampling chain too infrequently\n", spacer);
                    free (tempStr);
                    return (ERROR);
                    }
                chainParams.sampleFreq = tempI;
                MrBayesPrint ("%s   Setting sample frequency to %d\n", spacer, chainParams.sampleFreq);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free (tempStr);
                return (ERROR);
                }
            }
        /* set Printfreq (printFreq) **********************************************************/
        else if (!strcmp(parmName, "Printfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   Printing to screen too infrequently\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                chainParams.printFreq = tempI;
                MrBayesPrint ("%s   Setting print frequency to %d\n", spacer, chainParams.printFreq);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Printmax (printMax) **********************************************************/
        else if (!strcmp(parmName, "Printmax"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   You need to print at least one chain\n", spacer);
                    return (ERROR);
                    }
                chainParams.printMax = tempI;
                MrBayesPrint ("%s   Setting maximum number of chains to print to screen to %d\n", spacer, chainParams.printMax);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Printall (printAll) ********************************************************/
        else if (!strcmp(parmName, "Printall"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.printAll = YES;
                    else
                        chainParams.printAll = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Printall\n", spacer);
                    free (tempStr);
                    return (ERROR);
                    }
                if (chainParams.allChains == YES)
                    MrBayesPrint ("%s   Printing all chains to screen\n", spacer);
                else
                    MrBayesPrint ("%s   Printing only cold chains to screen\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Swapfreq (swapFreq) ************************************************************/
        else if (!strcmp(parmName, "Swapfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   Swapping states too infrequently\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                chainParams.swapFreq = tempI;
                MrBayesPrint ("%s   Setting swap frequency to %d\n", spacer, chainParams.swapFreq);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free (tempStr);
                return (ERROR);
                }
            }
        /* set Nswaps (numSwaps) ************************************************************/
        else if (!strcmp(parmName, "Nswaps"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   There must be at least one swap per swapping cycle\n", spacer);
                    free (tempStr);
                    return (ERROR);
                    }
                chainParams.numSwaps = tempI;
                MrBayesPrint ("%s   Setting number of swaps per swapping cycle to %d\n", spacer, chainParams.numSwaps);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Allchains (allChains) ********************************************************/
        else if (!strcmp(parmName, "Allchains"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.allChains = YES;
                    else
                        chainParams.allChains = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Allchains\n", spacer);
                    free (tempStr);
                    return (ERROR);
                    }
                if (chainParams.allChains == YES)
                    MrBayesPrint ("%s   Calculating MCMC diagnostics for all chains\n", spacer);
                else
                    MrBayesPrint ("%s   Calculating MCMC diagnostics only for cold chain(s)\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free (tempStr);
                return (ERROR);
                }
            }
        /* set Allcomps (allComps) ************************************************************/
        else if (!strcmp(parmName, "Allcomps"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.allComps = YES;
                    else
                        chainParams.allComps = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Allcomps\n", spacer);
                    free (tempStr);
                    return (ERROR);
                    }
                if (chainParams.allComps == YES)
                    MrBayesPrint ("%s   Calculating MCMC diagnostics for all pairwise run comparisons\n", spacer);
                else
                    MrBayesPrint ("%s   Only calculating overall MCMC diagnostics\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Mcmcdiagn (mcmcDiagn) ********************************************************/
        else if (!strcmp(parmName, "Mcmcdiagn"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.mcmcDiagn = YES;
                    else
                        chainParams.mcmcDiagn = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for mcmc diagnostics\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.mcmcDiagn == YES)
                    MrBayesPrint ("%s   Setting calculation of MCMC diagnostics ('Mcmcdiagn') to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting calculation of MCMC diagnostics ('Mcmcdiagn') to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Diagnfreq (diagnFreq) ************************************************************/
        else if (!strcmp(parmName, "Diagnfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   Diagnosing MCMC behavior too infrequently\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                chainParams.diagnFreq = tempI;
                MrBayesPrint ("%s   Setting diagnosing frequency to %d\n", spacer, chainParams.diagnFreq);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Savetrees (saveTrees) ********************************************************/
        else if (!strcmp(parmName, "Savetrees"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.saveTrees = YES;
                    else
                        chainParams.saveTrees = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Savetrees\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.saveTrees == YES)
                    MrBayesPrint ("%s   Saving trees for MCMC diagnostics in memory (if needed)\n", spacer);
                else
                    MrBayesPrint ("%s   Not saving trees for MCMC diagnostics in memory\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Diagnstat (diagnStat) ********************************************************/
        else if (!strcmp(parmName, "Diagnstat"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Avgstddev"))
                        chainParams.diagnStat = AVGSTDDEV;
                    else /* if (!strcmp(tempStr, "Maxstddev")) */
                        chainParams.diagnStat = MAXSTDDEV;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Savetrees\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.diagnStat == AVGSTDDEV)
                    MrBayesPrint ("%s   Setting diagnostics statistic to Avgstddev\n", spacer);
                else /* if (chainParams.diagnStat == MAXSTDDEV) */
                    MrBayesPrint ("%s   Setting diagnostics statistic to Maxstddev\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Checkpoint (checkPoint) ********************************************************/
        else if (!strcmp(parmName, "Checkpoint"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.checkPoint = YES;
                    else
                        chainParams.checkPoint = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for 'Checkpoint' (check-pointing)\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.checkPoint == YES)
                    MrBayesPrint ("%s   Setting check-pointing ('Checkpoint') to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting check-pointing ('Checkpoint') to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Checkfreq (checkFreq) ************************************************************/
        else if (!strcmp(parmName, "Checkfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 100)
                    {
                    MrBayesPrint ("%s   Check-pointing frequency must be at least 100\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                chainParams.checkFreq = tempI;
                MrBayesPrint ("%s   Setting check-pointing frequency to %d\n", spacer, chainParams.checkFreq);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Minpartfreq (minPartFreq) ************************************************************/
        else if (!strcmp(parmName, "Minpartfreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                if (tempD < 0.01)
                    {
                    MrBayesPrint ("%s   Minimum partition frequency too low (< 0.01)\n", spacer);
                    free (tempStr);
                    return (ERROR);
                    }
                if (tempD > 0.8)
                    {
                    MrBayesPrint ("%s   Minimum partition frequency too high (> 0.8)\n", spacer);
                    free (tempStr);
                    return (ERROR);
                    }
                chainParams.minPartFreq = tempD;
                MrBayesPrint ("%s   Setting minimum partition frequency to %.2f\n", spacer, chainParams.minPartFreq);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Nruns (numRuns) ****************************************************************/
        else if (!strcmp(parmName, "Nruns"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   Too few runs (minimum of 1 run)\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (tempI > MAX_RUNS)
                    {
                    MrBayesPrint ("%s   Too many runs (maximum of %d runs)\n", spacer, MAX_RUNS);
                    free(tempStr);
                    return (ERROR);
                    }
                if (ChangeNumRuns (chainParams.numRuns, tempI) == ERROR)
                    return (ERROR);
                chainParams.numRuns = tempI;
                fileNameChanged = YES;
                MrBayesPrint ("%s   Setting number of runs to %d\n", spacer, chainParams.numRuns);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Nchains (numChains) ************************************************************/
        else if (!strcmp(parmName, "Nchains"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 1)
                    {
                    MrBayesPrint ("%s   Too few chains (minimum of 1 chain)\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (tempI > MAX_CHAINS)
                    {
                    MrBayesPrint ("%s   Too many chains (maximum of %d chains)\n", spacer, MAX_CHAINS);
                    free(tempStr);
                    return (ERROR);
                    }
                if (ChangeNumChains (chainParams.numChains, tempI) == ERROR)
                    return (ERROR);
                chainParams.numChains = tempI;
                MrBayesPrint ("%s   Setting number of chains to %d\n", spacer, chainParams.numChains);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Temp (chainTemp) ***************************************************************/
        else if (!strcmp(parmName, "Temp"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                tempIndex = 0;
                expecting = Expecting(NUMBER) | Expecting(LEFTPAR);
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                chainParams.userDefinedTemps = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                MrBayesPrint ("%s   Setting user-defined temperatures\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else if (expecting == Expecting(COMMA))
                {
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                if (chainParams.userDefinedTemps == NO)
                    {
                    sscanf (tkn, "%lf", &tempD);
                    chainParams.chainTemp = tempD;
                    MrBayesPrint ("%s   Setting heating parameter to %lf\n", spacer, chainParams.chainTemp);
                    expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                    }
                else
                    {
                    if (tempIndex >= MAX_CHAINS)
                        {
                        MrBayesPrint ("%s   Too many user-defined temperatures (%d maximum)\n", spacer, MAX_CHAINS);
                        free(tempStr);
                        return (ERROR);
                        }
                    sscanf (tkn, "%lf", &tempD);
                    chainParams.userTemps[tempIndex++] = tempD;
                    expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
                    }
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Reweight (weightScheme) ********************************************************/
        else if (!strcmp(parmName, "Reweight"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(LEFTPAR);
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting = Expecting(NUMBER);
                whichReweightNum = 0;
                }
            else if (expecting == Expecting(NUMBER))
                {
                if (whichReweightNum < 0 || whichReweightNum > 2)
                    {
                    free(tempStr);
                    return (ERROR);
                    }
                sscanf (tkn, "%lf", &tempD);
                chainParams.weightScheme[whichReweightNum] = tempD;
                if (whichReweightNum < 2)
                    {
                    if (tempD < 0.0 || tempD > 100.0)
                        {
                        MrBayesPrint ("%s   The reweighting parameter must be between 0 and 100\n", spacer);
                        chainParams.weightScheme[0] = chainParams.weightScheme[1] = 0.0;
                        chainParams.weightScheme[2] = 1.0;
                        free(tempStr);
                        return (ERROR);
                        }
                    }
                else
                    {
                    if (tempD <= 0.0 || tempD > 1.0)
                        {
                        MrBayesPrint ("%s   The reweighting increment must be between 0 and 1\n", spacer);
                        chainParams.weightScheme[0] = chainParams.weightScheme[1] = 0.0;
                        chainParams.weightScheme[2] = 1.0;
                        free(tempStr);
                        return (ERROR);
                        }
                    }
                if (whichReweightNum == 0)
                    {
                    expecting = Expecting(COMMA);
                    }
                else if (whichReweightNum == 1)
                    {
                    if (chainParams.weightScheme[0] + chainParams.weightScheme[1] > 100.0)
                        {
                        MrBayesPrint ("%s   The sum of the reweighting parameters cannot exceed 100 %%\n", spacer);
                        chainParams.weightScheme[0] = chainParams.weightScheme[1] = 0.0;
                        chainParams.weightScheme[2] = 1.0;
                        free(tempStr);
                        return (ERROR);
                        }
                    expecting = Expecting(COMMA) | Expecting(RIGHTPAR);
                    }
                else
                    {
                    expecting = Expecting(RIGHTPAR);
                    }
                whichReweightNum++;
                }
            else if ((expecting & Expecting(COMMA)) == Expecting(COMMA))
                expecting = Expecting(NUMBER);
            else if ((expecting & Expecting(RIGHTPAR)) == Expecting(RIGHTPAR))
                {
                if (chainParams.weightScheme[0] >= 100.0)
                    {
                    MrBayesPrint ("%s   Cannot decrease weight of all characters\n", spacer);
                    chainParams.weightScheme[0] = chainParams.weightScheme[1] = 0.0;
                    chainParams.weightScheme[2] = 1.0;
                    free(tempStr);
                    return (ERROR);
                    }
                MrBayesPrint ("%s   Setting reweighting parameter to (%1.2lf v, %1.2lf ^) increment = %1.2lf\n", 
                    spacer, chainParams.weightScheme[0], chainParams.weightScheme[1], chainParams.weightScheme[2]);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Filename (chainFileName) *******************************************************/
        else if (!strcmp(parmName, "Filename"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                sscanf (tkn, "%s", tempStr);
                if (strlen(tempStr)>99)
                    {
                    MrBayesPrint ("%s   Maximum allowed length of chain file name is 99 characters. The given name:\n", spacer);
                    MrBayesPrint ("%s      '%s'\n", spacer,tempStr);
                    MrBayesPrint ("%s   has %d characters.\n", spacer,strlen(tempStr));
                    return (ERROR);
                    }
                strcpy (chainParams.chainFileName, tempStr);
                fileNameChanged = YES;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Relburnin (relativeBurnin) ********************************************************/
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
                    free(tempStr);
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
                free (tempStr);
                return (ERROR);
                }
            }
        /* set Burnin (chainBurnIn) ***********************************************************/
        else if (!strcmp(parmName, "Burnin"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                chainParams.chainBurnIn = tempI;
                MrBayesPrint ("%s   Setting chain burn-in to %d\n", spacer, chainParams.chainBurnIn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Burninfrac (burninFraction) ************************************************************/
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
                chainParams.burninFraction = tempD;
                MrBayesPrint ("%s   Setting burnin fraction to %.2f\n", spacer, chainParams.burninFraction);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Stoprule (stopRule) ********************************************************/
        else if (!strcmp(parmName, "Stoprule"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.stopRule = YES;
                    else
                        chainParams.stopRule = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Stoprule\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.stopRule == YES)
                    MrBayesPrint ("%s   Using stopping rule.\n", spacer);
                else
                    MrBayesPrint ("%s   Not using stopping rule.\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Stopval (stopVal) ************************************************************/
        else if (!strcmp(parmName, "Stopval"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%lf", &tempD);
                if (tempD < 0.000001)
                    {
                    MrBayesPrint ("%s   Stop value too low (< 0.000001)\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (tempD > 0.20)
                    {
                    MrBayesPrint ("%s   Stop value too high (> 0.20)\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                chainParams.stopVal = tempD;
                MrBayesPrint ("%s   Setting burnin fraction to %.2f\n", spacer, chainParams.burninFraction);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Starttree (startTree) **************************************************/
        else if (!strcmp(parmName, "Starttree"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr,"User"))
                        {
                        MrBayesPrint ("%s   The 'user' setting of 'Starttree' is deprecated. Set starting trees using 'Startvals' instead.\n", spacer);
                        free (tempStr);
                        return (ERROR);
                        }
                    else
                        strcpy(chainParams.startTree, tempStr);
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid 'Starttree' argument '%s'. This is a bug -- please report.\n", spacer, tempStr);
                    free(tempStr);
                    return (ERROR);
                    }
                MrBayesPrint ("%s   Setting 'Starttree' to \"%s\"\n", spacer, chainParams.startTree);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Startingtrees (deprecated) **************************************************/
        else if (!strcmp(parmName, "Startingtrees"))
            {
            free (tempStr);
            MrBayesPrint ("%s   Parameter 'Startingtrees' is deprecated. Use the 'Starttree' parameter or the 'Startvals' command instead.\n", spacer);
            return (ERROR);
            }
        /* set Nperts (numStartPerts) *********************************************************/
        else if (!strcmp(parmName, "Nperts"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                chainParams.numStartPerts = tempI;
                MrBayesPrint ("%s   Setting number of perturbations to start tree to %d\n", spacer, chainParams.numStartPerts);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Startparams (startParams) **************************************************/
        else if (!strcmp(parmName, "Startparams"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    strcpy(chainParams.startParams, tempStr);
                else
                    {
                    MrBayesPrint ("%s   Invalid 'Startparams' argument\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                MrBayesPrint ("%s   Setting 'Startparams' to \"%s\"\n", spacer, chainParams.startParams);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Savebrlens (saveBrlens) ********************************************************/
        else if (!strcmp(parmName, "Savebrlens"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        ; /* This is the only option in version 3.2 */
                    else
                        MrBayesPrint ("%s   WARNING: Ignoring savebrlens setting; since version 3.2, branch lengths are always saved\n", spacer);
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for savebrlens\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            }
        /* set Redirect (redirect) ********************************************************/
        else if (!strcmp(parmName, "Redirect"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.redirect = YES;
                    else
                        chainParams.redirect = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for redirecting output\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.redirect == YES)
                    MrBayesPrint ("%s   Setting program to redirect output\n", spacer);
                else
                    MrBayesPrint ("%s   Setting program not to redirect output\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Data (runWithData) ************************************************************/
        else if (!strcmp(parmName, "Data"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.runWithData = YES;
                    else
                        chainParams.runWithData = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Data\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.runWithData == NO)
                    MrBayesPrint ("%s   Running without data (WARNING: use this only for checking priors!)\n", spacer);
                else
                    MrBayesPrint ("%s   Running with data (standard analysis)\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Ordertaxa (chainParams.orderTaxa) *********************************************/
        else if (!strcmp(parmName, "Ordertaxa"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.orderTaxa = YES;
                    else
                        chainParams.orderTaxa = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for ordertaxa\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.orderTaxa == YES)
                    MrBayesPrint ("%s   Setting ordertaxa to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting ordertaxa to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Append (chainParams.append) *********************************************/
        else if (!strcmp(parmName, "Append"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.append = YES;
                    else
                        chainParams.append = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for append\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.append == YES)
                    MrBayesPrint ("%s   Setting append to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting append to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Autotune (chainParams.autotune) *********************************************/
        else if (!strcmp(parmName, "Autotune"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.autotune = YES;
                    else
                        chainParams.autotune = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Autotune\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.autotune == YES)
                    MrBayesPrint ("%s   Setting Autotune to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting Autotune to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Tunefreq (tuneFreq) ************************************************************/
        else if (!strcmp(parmName, "Tunefreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 10)
                    {
                    MrBayesPrint ("%s   Autotuning frequency must be at least 10\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                chainParams.tuneFreq = tempI;
                MrBayesPrint ("%s   Setting autotuning frequency to %d\n", spacer, chainParams.tuneFreq);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                {
                free(tempStr);
                return (ERROR);
                }
            }
        /* set Swapadjacent (swapAdjacentOnly) **************************************************/
        else if (!strcmp(parmName, "Swapadjacent"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        chainParams.swapAdjacentOnly = YES;
                    else
                        chainParams.swapAdjacentOnly = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Swapadjacent\n", spacer);
                    free(tempStr);
                    return (ERROR);
                    }
                if (chainParams.swapAdjacentOnly == YES)
                    MrBayesPrint ("%s   Setting program to attempt swaps only between chains of adjacent temperatures\n", spacer);
                else
                    MrBayesPrint ("%s   Setting program to attempt all possible swaps between chains\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                free(tempStr);
                return (ERROR);
                }
            }
        else
            {
            free(tempStr);
            return (ERROR);
            }
        }
    free(tempStr);
    return (NO_ERROR);
}


int DoSs (void)
{
    int ret, oldBurnin;

    if (chainParams.numGen/chainParams.sampleFreq <= chainParams.burninSS)
        {/*Do not change print out to generations vs samples because of danger of overflow*/
        MrBayesPrint ("%s      ERROR: Burnin %d samples is too large compared with requested total %d samples (%d generations).\n", spacer ,chainParams.burninSS, chainParams.numGen/chainParams.sampleFreq, chainParams.numGen);
        return ERROR;
        }

    oldBurnin = chainParams.burninSS;;
    stepRelativeBurninSS = chainParams.relativeBurnin;

    chainParams.relativeBurnin = YES;
 
    if (chainParams.burninSS < 0)
        chainParams.burninSS =  chainParams.numGen / ((chainParams.numStepsSS-chainParams.burninSS)*chainParams.sampleFreq);
    chainParams.isSS = YES;

    ret=DoMcmc();

    chainParams.isSS = NO;
    chainParams.burninSS = oldBurnin;
    chainParams.relativeBurnin = stepRelativeBurninSS;

    return ret;
}


int DoSsp (void)
{
    return NO_ERROR;
}


int ExhaustiveParsimonySearch (Tree *t, int chain, TreeInfo *tInfo)
{
    int         i, j, k;
    TreeNode    *p;
    
    for (i=j=k=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL || p->right == NULL)
            tInfo->leaf[j++] = p;
        else
            tInfo->vertex[k++] = p;
        }

    tInfo->leaf[0]->anc = tInfo->leaf[1]->anc = tInfo->vertex[0];
    tInfo->vertex[0]->left = tInfo->leaf[0];
    tInfo->vertex[0]->right = tInfo->leaf[1];
    tInfo->leaf[t->nIntNodes+1]->left = tInfo->vertex[0];
    tInfo->vertex[0]->anc = tInfo->leaf[t->nIntNodes+1];

    BuildExhaustiveSearchTree (t, chain, 2, tInfo);

    return (NO_ERROR);
}


int ExtendChainQuery ()
{
    int             extendChain, additionalCycles;
    char            s[100];
    
#   if defined (MPI_ENABLED)
    if (proc_id == 0)
        {
        MrBayesPrint ("\n");
        extendChain = WantTo ("Continue with analysis");
        }
    MPI_Bcast (&extendChain, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (extendChain == YES)
        {
        if (proc_id == 0)
            {
            additionalCycles = 0;
            do
                {
                if (additionalCycles < 0)
                    MrBayesPrint ("%s      Number must be greater than or equal to 0: ", spacer);
                else
                    MrBayesPrint ("%s      Additional number of generations: ", spacer);

                if (fgets (s, 100, stdin) == NULL)
                    {
                    printf ("Error in function: %s at line: %d in file: %s", __FUNCTION__, __LINE__, __FILE__);
                    }
                sscanf (s, "%d", &additionalCycles);

                } while (additionalCycles < 0);
            MrBayesPrint ("\n");
            }
        MPI_Bcast (&additionalCycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
            
        return (additionalCycles);
        }
    else
        return (0);
#   else
    MrBayesPrint ("\n");
    extendChain = WantTo ("Continue with analysis");
    
    if (extendChain == YES)
        {
        additionalCycles = 0;
        do
            {
            if (additionalCycles < 0)
                MrBayesPrint ("%s      Number must be greater than or equal to 0: ", spacer);
            else
                MrBayesPrint ("%s      Additional number of generations: ", spacer);

            if (fgets (s, 20, stdin) == NULL)
                {
                printf ("Error in function: %s at line: %d in file: %s", __FUNCTION__, __LINE__, __FILE__);
                }
            sscanf (s, "%d", &additionalCycles);

            } while (additionalCycles < 0);
        MrBayesPrint ("\n");
        return (additionalCycles);
        }
    else
        return (0);
#   endif

}


int FillNumSitesOfPat (void)
{
    int         i, j, n, *increased, *decreased, nToDecrease, nToIncrease, whichToChange;
    MrBFlt      ran, sum;
    CLFlt       wtIncrement;
    
    wtIncrement = (CLFlt) chainParams.weightScheme[2];
    increased = decreased = NULL;
    
    /* reallocate numSitesOfPat */
    if (memAllocs[ALLOC_NUMSITESOFPAT] == NO)
        {
        MrBayesPrint ("%s   numSitesOfPat is not allocated\n", spacer);
        goto errorExit;
        }
    memAllocs[ALLOC_NUMSITESOFPAT] = NO;
    numSitesOfPat = (CLFlt *) SafeRealloc((void *) numSitesOfPat, numCompressedChars * chainParams.numChains * sizeof(MrBFlt));
    if (!numSitesOfPat)
        {
        MrBayesPrint ("%s   Problem reallocating numSitesOfPat (%d)\n", spacer, numCompressedChars * chainParams.numChains * sizeof(MrBFlt));
        goto errorExit;
        }
    memAllocs[ALLOC_NUMSITESOFPAT] = YES;

    /* copy first numCompressedChars into the remaining bits */
    if (chainParams.numChains > 1)
        {
        for (i=0; i<numCompressedChars; i++)
            {
            for (j=1; j<chainParams.numChains; j++)
                {
                numSitesOfPat[j * numCompressedChars + i] = numSitesOfPat[i];
                }
            }
        }   
        
    /* reweight characters for each chain */
    if (chainParams.numChains > 1)
        {
        if (chainParams.weightScheme[0] + chainParams.weightScheme[1] > 0.0001)
            MrBayesPrint ("%s   Reweighting of characters for chains 1 to %d\n", spacer, chainParams.numChains);

        /* check that we don't have an HMM */
        if (chainHasAdgamma == YES && chainParams.weightScheme[0] + chainParams.weightScheme[1] > 0.0001)
            {
            MrBayesPrint ("%s   Reweighting of characters is not allowed with an autocorrelated gamma model\n", spacer);
            goto errorExit;
            }
        
        /* how many characters */
        n = 0;
        for (i=0; i<numCompressedChars; i++)
            n += (int)numSitesOfPat[0 * numCompressedChars + i];
        nToDecrease = (int)(n * chainParams.weightScheme[0] / 100.0);
        nToIncrease = (int)(n * chainParams.weightScheme[1] / 100.0);
        if (chainParams.weightScheme[0] + chainParams.weightScheme[1] > 0.0001)
            {
            MrBayesPrint ("%s      Decreasing weight of %d characters\n", spacer, nToDecrease);
            MrBayesPrint ("%s      Increasing weight of %d characters\n", spacer, nToIncrease);
            }
        
        /* allocate memory */
        increased = (int *)SafeMalloc(2 * (size_t)numCompressedChars * sizeof(int));
        if (!increased)
            {
            MrBayesPrint ("%s   Problem reallocating increased (%d)\n", spacer, numCompressedChars * chainParams.numChains * sizeof(int));
            goto errorExit;
            }
        decreased = increased + numCompressedChars;

        /* reweight characters for each chain */
        for (j=1; j<chainParams.numChains; j++)
            {
            for (i=0; i<numCompressedChars; i++)
                increased[i] = decreased[i] = 0;

            /* decrease weight of characters */
            for (i=0; i<nToDecrease; i++)
                {
                do
                    {
                    ran = RandomNumber(&swapSeed);
                    sum = 0.0;
                    for (whichToChange=0; whichToChange<numCompressedChars; whichToChange++)
                        {
                        sum += numSitesOfPat[0 * numCompressedChars + whichToChange] / n;
                        if (ran < sum)
                            break;
                        }
                    if (whichToChange < 0 || whichToChange >= numCompressedChars)
                        continue;
                    } while (decreased[whichToChange] >= numSitesOfPat[0 * numCompressedChars + whichToChange]);
                decreased[whichToChange]++;
                numSitesOfPat[j * numCompressedChars + whichToChange] -= wtIncrement;
                if (numSitesOfPat[j * numCompressedChars + whichToChange] < 0)
                    {
                    MrBayesPrint ("%s   Problem reweighting characters\n", spacer);
                    goto errorExit;
                    }
                }

            /* increase weight of characters */
            for (i=0; i<nToDecrease; i++)
                {
                do
                    {
                    ran = RandomNumber(&swapSeed);
                    sum = 0.0;
                    for (whichToChange=0; whichToChange<numCompressedChars; whichToChange++)
                        {
                        sum += numSitesOfPat[0 * numCompressedChars + whichToChange] / n;
                        if (ran < sum)
                            break;
                        }
                    if (whichToChange < 0 || whichToChange >= numCompressedChars)
                        continue;
                    } while ((increased[whichToChange] + decreased[whichToChange]) >= numSitesOfPat[0 * numCompressedChars + whichToChange]);
                increased[whichToChange]++;
                numSitesOfPat[j * numCompressedChars + whichToChange] += wtIncrement;
                if (numSitesOfPat[j * numCompressedChars + whichToChange] < 0)
                    {
                    MrBayesPrint ("%s   Problem reweighting characters\n", spacer);
                    goto errorExit;
                    }
                }

            }
            
        /* free allocated memory */
        free (increased);
        }
        
#   if 0
    /* print site patterns for each chain */
    for (i=0; i<numCompressedChars; i++)
        {
        MrBayesPrint ("%4d -- ", i);
        for (j=0; j<chainParams.numChains; j++)
            {
            MrBayesPrint ("%4.1lf ", numSitesOfPat[j * numCompressedChars + i]);
            }
        MrBayesPrint ("\n");
        }
#   endif

    return (NO_ERROR);
    
    errorExit:
        if (increased)
            free (increased);
        return (ERROR);
}


/* FindBestNode: Recursive function for finding best attachment point */
TreeNode *FindBestNode (Tree *t, TreeNode *p, TreeNode *addNode, CLFlt *minLength, int chain) {

    int         c, n, division;
    TreeNode    *q=NULL, *r=NULL;
    BitsLong    *pA, *pP, *pX;
    CLFlt       *nSitesOfPat, fpLength, length;
    ModelInfo   *m;

    /* Calculate length, looping over divisions */
    fpLength = 0;
    for (n=0; n<t->nRelParts; n++)
        {
        division = t->relParts[n];

        /* Find model settings */
        m = &modelSettings[division];

        /* Find number of site patterns */
        nSitesOfPat = numSitesOfPat + ((1 % chainParams.numChains) * numCompressedChars) + m->compCharStart;

        /* Find final-pass parsimony sets for the node and its ancestor */
        pP    = m->parsSets[p->index      ];
        pA    = m->parsSets[p->anc->index ];
        pX    = m->parsSets[addNode->index];

        for (c=0; c<m->numChars; c++)
            {
            if (((pP[c] | pA[c])&pX[c]) == 0)
                fpLength += nSitesOfPat[c];
            }
        }

    /* If tip, this is the best node and its length is the min length */
    if (p->left == NULL)
        {
        *minLength = fpLength;
        return p;
        }

    /* Find best node and min length above this node */
    if (p->left != NULL) {
        q = FindBestNode(t, p->left,  addNode, minLength, chain);
        r = FindBestNode(t, p->right, addNode, &length,   chain);
        if (length < *minLength) {
            *minLength = length;
            q = r;
        }
    }

    /* Return best node and min length */
    if (*minLength < fpLength)
        return q;
    else /* This node is the best */
        {
        *minLength = fpLength;
        return p;
        }
}


/* FlipCijkSpace: Flip space for cijks with scratch area */
void FlipCijkSpace (ModelInfo* m, int chain)
{
    int         temp;

    temp                = m->cijkIndex[chain];
    m->cijkIndex[chain] = m->cijkScratchIndex;
    m->cijkScratchIndex = temp;
}


/* FlipCondLikeSpace: Flip space for conditional likelihoods with scratch area */
void FlipCondLikeSpace (ModelInfo* m, int chain, int nodeIndex)
{
    int         temp;

    temp                               = m->condLikeIndex[chain][nodeIndex];
    m->condLikeIndex[chain][nodeIndex] = m->condLikeScratchIndex[nodeIndex];
    m->condLikeScratchIndex[nodeIndex] = temp;
}


/* FlipNodeScalerSpace: Flip space for node scalers and scaler flag with scratch area */
void FlipNodeScalerSpace (ModelInfo* m, int chain, int nodeIndex)
{
    int         temp;

    temp                                 = m->nodeScalerIndex[chain][nodeIndex];
    m->nodeScalerIndex[chain][nodeIndex] = m->nodeScalerScratchIndex[nodeIndex];
    m->nodeScalerScratchIndex[nodeIndex] = temp;

    temp                                 = m->unscaledNodes[chain][nodeIndex];
    m->unscaledNodes[chain][nodeIndex]   = m->unscaledNodesScratch[nodeIndex];
    m->unscaledNodesScratch[nodeIndex]   = temp;
}


/* FlipSiteScalerSpace: Flip space for ln site scalers */
void FlipSiteScalerSpace (ModelInfo *m, int chain)
{
    int  temp;

#   if defined (BEAGLE_ENABLED)
    int *tempp;
#   endif

    temp = m->siteScalerIndex[chain];
    m->siteScalerIndex[chain] = m->siteScalerScratchIndex;
    m->siteScalerScratchIndex = temp;

#   if defined (BEAGLE_ENABLED)
    if (m->useBeagle == YES)
        {
        tempp = m->isScalerNode[chain];
        m->isScalerNode[chain] = m->isScalerNodeScratch ;
        m->isScalerNodeScratch = tempp;
        }
#   endif
}


/* FlipTiProbsSpace: Flip space for ti probs with scratch area */
void FlipTiProbsSpace (ModelInfo* m, int chain, int nodeIndex)
{
    int         temp;

    temp                              = m->tiProbsIndex[chain][nodeIndex];
    m->tiProbsIndex[chain][nodeIndex] = m->tiProbsScratchIndex[nodeIndex];
    m->tiProbsScratchIndex[nodeIndex] = temp;
}


void FreeChainMemory (void)
{
    int         i, j, k, nRates;
    ModelInfo   *m;

    /* free model variables for Gibbs gamma */
    for (i=0; i<numCurrentDivisions; i++)
        {
        if (modelSettings[i].gibbsGamma == YES)
            {
            if (modelSettings[i].pInvar != NULL)
                nRates = modelSettings[i].numRateCats + 1;
            else
                nRates = modelSettings[i].numRateCats;
            for (j=0; j<numLocalChains; j++)
                {
                for (k=0; k<nRates; k++)
                    {
                    free(modelSettings[i].catLnScaler[j][k]);
                    free(modelSettings[i].catLike[j][k]);
                    }
                free (modelSettings[i].catLnScaler[j]);
                free (modelSettings[i].catLike[j]);
                }
            free (modelSettings[i].tiIndex);
            free (modelSettings[i].catLike);
            free (modelSettings[i].catLnScaler);
            }
        }

    /* free parsimony sets and node lens */
    for (i=0; i<numCurrentDivisions; i++)
        {
        m = &modelSettings[i];
        if (m->parsSets)
            {
            for (j=0; j<m->numParsSets; j++)
                free (m->parsSets[j]);
            free (m->parsSets);
            m->parsSets = NULL;
            }
        if (m->parsNodeLens)
            {
            free(m->parsNodeLens);
            m->parsNodeLens = NULL;
            }
        }

    /* free model variables for conditional likelihoods */
    for (i=0; i<numCurrentDivisions; i++)
        {
        m = &modelSettings[i];
        if (m->condLikes)
            {
            for (j=0; j<m->numCondLikes; j++)
                {
#   if defined (SSE_ENABLED)
                if (m->useVec != NO)
                    AlignedSafeFree ((void *)(&m->condLikes[j]));
                else
                    free (m->condLikes[j]);
#   else
                free (m->condLikes[j]);
#   endif
                }
            free (m->condLikes);
            m->condLikes = NULL;
            }

        if (m->scalers)
            {
            for (j=0; j<m->numScalers; j++)
#   if defined (SSE_ENABLED)
            if (m->useVec != VEC_NONE)
                AlignedSafeFree ((void *)(&m->scalers[j]));
            else
                free (m->scalers[j]);
#   else
                free (m->scalers[j]);

#   endif
            free (m->scalers);
            m->scalers = NULL;
            }

        if (m->clP)
            {
            free (m->clP);
            m->clP = NULL;
            }
#   if defined (SSE_ENABLED)
        if (m->useVec != VEC_NONE)
            {
            if (m->clP_SSE)
                {
                free (m->clP_SSE);
                m->clP_SSE = NULL;
                }
#if defined (AVX_ENABLED)
            if (m->clP_AVX)
                {
                free (m->clP_AVX);
                m->clP_AVX = NULL;
                }
#endif
            if (m->lnL_Vec)
                AlignedSafeFree ((void *)(&m->lnL_Vec));
            if (m->lnLI_Vec)
                AlignedSafeFree ((void *)(&m->lnLI_Vec));
            }
#   endif

        if (m->tiProbs)
            {
            for (j=0; j<m->numTiProbs; j++)
                free (m->tiProbs[j]);
            free (m->tiProbs);
            m->tiProbs = NULL;
            }
                
        if (m->cijks)
            {
            for (j=0; j<numLocalChains+1; j++)
                free (m->cijks[j]);
            free (m->cijks);
            m->cijks = NULL;
            }

        if (m->condLikeIndex)
            {
            for (j=0; j<numLocalChains; j++)
                free (m->condLikeIndex[j]);
            free (m->condLikeIndex);
            m->condLikeIndex = NULL;
            }

        if (m->condLikeScratchIndex)
            {
            free (m->condLikeScratchIndex);
            m->condLikeScratchIndex=NULL;
            }

        if (m->tiProbsIndex)
            {
            for (j=0; j<numLocalChains; j++)
                free (m->tiProbsIndex[j]);
            free (m->tiProbsIndex);
            m->tiProbsIndex = NULL;
            }
        if (m->tiProbsScratchIndex)
            {
            free (m->tiProbsScratchIndex);
            m->tiProbsScratchIndex = NULL;
            }
        if (m->nodeScalerIndex)
            {
            for (j=0; j<numLocalChains; j++)
                free (m->nodeScalerIndex[j]);
            free (m->nodeScalerIndex);
            m->nodeScalerIndex = NULL;
            }
        if (m->nodeScalerScratchIndex)
            {
            free (m->nodeScalerScratchIndex);
            m->nodeScalerScratchIndex = NULL;
            }
        if (m->unscaledNodes)
            {
            for (j=0; j<numLocalChains; j++)
                free (m->unscaledNodes[j]);
            free (m->unscaledNodes);
            m->unscaledNodes = NULL;
            }
        if (m->unscaledNodesScratch)
            {
            free (m->unscaledNodesScratch);
            m->unscaledNodesScratch = NULL;
            }
        if (m->siteScalerIndex)
            {
            free (m->siteScalerIndex);
            m->siteScalerIndex = NULL;
            }
        if (m->cijkIndex)
            {
            free (m->cijkIndex);
            m->cijkIndex = NULL;
            }
        if (m->ancStateCondLikes)
            {
            free (m->ancStateCondLikes);
            m->ancStateCondLikes = NULL;
            }

#   if defined (BEAGLE_ENABLED)
        if (m->useBeagle == NO)
            continue;

        beagleFinalizeInstance(m->beagleInstance);
        SafeFree((void *)(&m->logLikelihoods));
        SafeFree((void *)(&m->inRates));
        SafeFree((void *)(&m->branchLengths));
        SafeFree((void *)(&m->tiProbIndices));
        SafeFree((void *)(&m->inWeights));
        SafeFree((void *)(&m->bufferIndices));
        SafeFree((void *)(&m->eigenIndices));
        SafeFree((void *)(&m->childBufferIndices));
        SafeFree((void *)(&m->childTiProbIndices));
        SafeFree((void *)(&m->cumulativeScaleIndices));

        m->isScalerNodeScratch += numLocalTaxa;
        SafeFree((void *)&(m->isScalerNodeScratch)); 
        for (j=0; j<numLocalChains; j++)
            {
            m->isScalerNode[j] += numLocalTaxa;
            SafeFree((void *)&(m->isScalerNode[j]));
            }
        SafeFree((void *)(&m->isScalerNode));

        SafeFree((void *)(&m->beagleComputeCount));
        SafeFree((void *)(&m->succesCount));
        SafeFree((void *)(&m->rescaleFreq));

#   endif
        }

    if (memAllocs[ALLOC_CURLNL] == YES) /*alloc in RunChain()*/
        {
        free (maxLnL0);
        free (curLnL);
        memAllocs[ALLOC_CURLNL] = NO;
        }
    if (memAllocs[ALLOC_SS] == YES) /*alloc in mcmc()*/
        {
        free (marginalLnLSS);
        free (stepScalerSS);
        free (stepAcumulatorSS);
        free (splitfreqSS);
        memAllocs[ALLOC_SS] = NO;
        }
    if (memAllocs[ALLOC_CURLNPR] == YES) /*alloc in RunChain()*/
        {
        free (curLnPr); 
        memAllocs[ALLOC_CURLNPR] = NO;
        }
    if (memAllocs[ALLOC_CHAINID] == YES) /*alloc in RunChain()*/
        {
        free (chainId); 
        memAllocs[ALLOC_CHAINID] = NO;
        }
    if (memAllocs[ALLOC_USEDMOVES] == YES) /*alloc in setUsedMoves()*/
        {
        free (usedMoves);
        memAllocs[ALLOC_USEDMOVES] = NO;
        }
    if (memAllocs[ALLOC_TERMSTATE] == YES) /*alloc in SetUpTermState()*/
        {
        free (termState);
        termState = NULL;  
        memAllocs[ALLOC_TERMSTATE] = NO;
        }
    if (memAllocs[ALLOC_ISPARTAMBIG] == YES) /*alloc in SetUpTermState()*/
        {
        free (isPartAmbig);
        isPartAmbig = NULL;
        memAllocs[ALLOC_ISPARTAMBIG] = NO;
        }
    if (memAllocs[ALLOC_PRELIKES] == YES) /*alloc in InitCondLike()*/
        {
        free (preLikeL);
        preLikeL = NULL;
        memAllocs[ALLOC_PRELIKES] = NO;
        }
    if (memAllocs[ALLOC_RATEPROBS] == YES) /*alloc in InitAdGamma() */
        {
        free (rateProbSpace);
        free (rateProbs);
        rateProbs = NULL;
        memAllocs[ALLOC_RATEPROBS] = NO;
        }
    if (memAllocs[ALLOC_SITEJUMP] == YES) /*alloc in InitAdGamma() */
        {
        free (siteJump);
        siteJump = NULL;
        memAllocs[ALLOC_SITEJUMP] = NO;
        }
    if (memAllocs[ALLOC_MARKOVTIS] == YES)  /*alloc in InitAdGamma() */
        {
        for (i=0; i<MAX_SMALL_JUMP; i++)
            if (markovTi[i] != NULL)
                FreeSquareDoubleMatrix(markovTi[i]);
        FreeSquareDoubleMatrix(markovTiN);
        memAllocs[ALLOC_MARKOVTIS] = NO;
        }
    if (memAllocs[ALLOC_SWAPINFO] == YES) /*alloc in RunChain()*/
        {
        for (i=0; i<chainParams.numRuns; i++)
            FreeSquareIntegerMatrix(swapInfo[i]);
        free (swapInfo);
        memAllocs[ALLOC_SWAPINFO] = NO;
        }
    if (memAllocs[ALLOC_POSSELPROBS] == YES) /*alloc in PrintStates() <- RunChain()*/
        {
        free (posSelProbs);
        memAllocs[ALLOC_POSSELPROBS] = NO;
        }
    if (memAllocs[ALLOC_PFCOUNTERS] == YES) /*alloc in SetUpParitionCounters() <- RunChain()*/
        {
        free (partition[0]);
        free (partition);
        for (i=0; i<numTopologies; i++)
            Tfree (partFreqTreeRoot[i]);
        free (partFreqTreeRoot);
        memAllocs[ALLOC_PFCOUNTERS] = NO;
        }
    if (memAllocs[ALLOC_FILEPOINTERS] == YES) /* alloc in (PreparePrintFiles(), ReusePreviousResults()) <- RunChain() */
        {
        CloseMBPrintFiles ();
        if (fpTree != NULL)
            {
            free (fpTree[0]);
            free (fpTree);
            }
        if (fpParm != NULL)
            free (fpParm);
        fpParm = NULL;
        fpTree = NULL;
        fpMcmc = NULL;
        fpSS = NULL;
        memAllocs[ALLOC_FILEPOINTERS] = NO;
        }
    if (memAllocs[ALLOC_STATS] == YES) /* alloc in RunChain() */
        {
        if (chainParams.allComps == YES)
            {
            for (i=0; i<numTopologies; i++)
                FreeSquareDoubleMatrix (chainParams.stat[i].pair);
            }
        free (chainParams.stat);
        memAllocs[ALLOC_STATS] = NO;
        }
    if (memAllocs[ALLOC_DIAGNTREE] == YES)
        {
        FreeTree (chainParams.dtree);
        memAllocs[ALLOC_DIAGNTREE] = NO;
        }
    if (memAllocs[ALLOC_PRINTPARAM] == YES)
        {
        free (printParam);
        free (topologyPrintIndex);
        memAllocs[ALLOC_PRINTPARAM] = NO;
        }
    if (memAllocs[ALLOC_TFILEPOS] == YES)
        {
        free (chainParams.tFilePos);
        chainParams.tFilePos = NULL;
        memAllocs[ALLOC_TFILEPOS] = NO;
        }
    if (memAllocs[ALLOC_TREELIST] == YES)
        {
        for (i=0; i<chainParams.numRuns * numTopologies; i++)
            EraseTreeList (&chainParams.treeList[i]);
        free (chainParams.treeList);
        chainParams.treeList = NULL;
        memAllocs[ALLOC_TREELIST] = NO;
        }
    if (memAllocs[ALLOC_BEST] == YES)
        {
        FreeBestChainVariables();
        memAllocs[ALLOC_BEST] = NO;
        }
}


MrBFlt GetFitchPartials (ModelInfo *m, int chain, int source1, int source2, int destination)
{
    int         c, i;
    BitsLong    x[2], *pS1, *pS2, *pD;
    MrBFlt      length = 0.0;
    CLFlt       *nSitesOfPat;
    
    assert (m->nParsIntsPerSite <= 2 && m->nParsIntsPerSite > 0);
    assert (source1 >= 0 && source1 < m->numParsSets);
    assert (source2 >= 0 && source2 < m->numParsSets);
    assert (destination >= 0 && destination < m->numParsSets);

    /* find parsimony sets for the nodes */
    pS1 = m->parsSets[source1    ];
    pS2 = m->parsSets[source2    ];
    pD  = m->parsSets[destination];
        
    /* Find number of site patterns */
    nSitesOfPat = numSitesOfPat + ((1 % chainParams.numChains) * numCompressedChars) + m->compCharStart;//chainId[chain]

    if (m->nParsIntsPerSite == 1)
        {
        for (c=0; c<m->numChars; c++)
            {
            x[0] = pS1[c] & pS2[c];
            if (x[0] == 0)
                {
                length += nSitesOfPat[c];
                x[0] = pS1[c] | pS2[c];
                }
            pD[c] = x[0];
            }
        }
    else /* if (m->nParsIntsPerSite == 2) */
        {
        assert (m->nParsIntsPerSite == 2);
        for (c=i=0; c<m->numChars; c++)
            {
            x[0] = pS1[i]   & pS2[i];
            x[1] = pS1[i+1] & pS2[i+1];
            if ((x[0] | x[1]) == 0)
                {
                length += nSitesOfPat[c];
                x[0] = pS1[i]   | pS2[i];
                x[1] = pS1[i+1] | pS2[i+1];
                }
            pD[i++] = x[0];
            pD[i++] = x[1];
            }
        }

    return length;
    MrBayesPrint ("%d", chain); /* just because I am tired of seeing the unused parameter error msg */
}


MrBFlt GetParsDP (Tree *t, TreeNode *p, int chain)
{
    int             n, division;
    MrBFlt          length;
    ModelInfo       *m;

    length = 0.0;
    if (p->left != NULL)
        {
        length += GetParsDP (t, p->left, chain);
        length += GetParsDP (t, p->right, chain);

        for (n=0; n<t->nRelParts; n++)
            {
            division = t->relParts[n];
            
            /* Find model settings */
            m = &modelSettings[division];

            /* Get Fitch partials and length */
            length += GetFitchPartials(m,
                                       chain,
                                       p->left->index,
                                       p->right->index,
                                       p->index);
            }
        }

    return length;
}


void GetParsFP (Tree *t, TreeNode *p, int chain)
{
    int             i, c, n, division;
    BitsLong        *pL, *pR, *pP, *pA, x[2];
    ModelInfo       *m;

    if (p->left != NULL)
        {
        for (n=0; n<t->nRelParts; n++)
            {
            division = t->relParts[n];
            
            /* Find model settings */
            m = &modelSettings[division];
            assert (m->nParsIntsPerSite == 1 || m->nParsIntsPerSite == 2);

            /* find parsimony sets for the node and its environment */
            pL   = m->parsSets[p->left->index ];
            pR   = m->parsSets[p->right->index];
            pP   = m->parsSets[p->index       ];
            pA   = m->parsSets[p->anc->index  ];
            
            if (m->nParsIntsPerSite == 1)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    x[0] = pP[c] & pA[c];

                    if (x[0] != pA[c])
                        {   /* means that we allow change of state from p to a */
                        if ((pL[c] & pR[c]) != 0)
                            x[0] = ((pL[c] | pR[c]) & pA[c]) | pP[c] ;
                            /* We still allow only one change from both children of p through p to a.
                               So state from a that belong to one of the children of p can be added to pP[c],
                               if p assume the state then the only change would be on the other child. */
                        else
                            x[0] = pP[c] | pA[c];
                            /* Here we allow two change from both children of p through p to a.
                               Adding pA[c] to pP[c] means that if p assume state exclusive in a then even if
                               both children will be in different state from p we still get optimal parsimony */
                        }
                    pP[c] = x[0];
                    }
                }
            else /* if (m->nParsIntsPerSite == 2) */
                {
                for (c=i=0; c<m->numChars; c++)
                    {
                    x[0] = pP[i]   & pA[i];
                    x[1] = pP[i+1] & pA[i+1];
                    if (x[0] != pA[i] || x[1] != pA[i+1])
                        {
                        x[0] = pL[i] & pR[i];
                        x[1] = pL[i+1] & pR[i+1];
                        if ((x[0] | x[1]) != 0)
                            {
                            x[0] = ((pL[i] | pR[i]) & pA[i]) | pP[i] ;
                            x[1] = ((pL[i+1] | pR[i+1]) & pA[i+1]) | pP[i+1] ;
                            }
                        else
                            {
                            x[0] = pP[i] | pA[i];
                            x[1] = pP[i+1] | pA[i+1];
                            }
                        }
                    pP[i++] = x[0];
                    pP[i++] = x[1];
                    }
                }
            }
        GetParsFP (t, p->left, chain);
        GetParsFP (t, p->right, chain);
        }
}


int GetParsimonyBrlens (Tree *t, int chain, MrBFlt *brlens)
{
    int             c, i, j, n, division;
    BitsLong        *pP, *pA;
    CLFlt           *nSitesOfPat;
    TreeNode        *p;
    ModelInfo       *m;

    /* Reset all brlens */
    for (i=0; i<t->nNodes-1; i++)
        brlens[i] = 0.0;
    
    /* Get final parsimony state sets */
    GetParsDP(t, t->root->left, chain);
    GetParsFP(t, t->root, chain);
    
    /* Get all branch lengths, looping over divisions */
    for (n=0; n<t->nRelParts; n++)      
        {
        division = t->relParts[n];

        /* Find model settings */
        m = &modelSettings[division];

        /* Find number of site patterns */
        nSitesOfPat = numSitesOfPat + ((chainId[chain] % chainParams.numChains) * numCompressedChars) + m->compCharStart;

        /* Record branch lengths in downpass */
        for (i=0; i<t->nNodes-1; i++)
            {
            p = t->allDownPass[i];

            /* Find final-pass parsimony sets for the node and its ancestor */
            pP    = m->parsSets[p->index     ];
            pA    = m->parsSets[p->anc->index];
            
            if (m->nParsIntsPerSite == 1)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    if ((pP[c] & pA[c]) == 0)
                        brlens[i] += nSitesOfPat[c];
                    }
                }
            else
                {
                for (c=j=0; c<m->numChars; c++, j+=2)
                    {
                    if ((pP[j] & pA[j]) == 0 && (pP[j+1] & pA[j+1]) == 0)
                        brlens[i] += nSitesOfPat[c];
                    }
                }
            }
        }

    return (NO_ERROR);
}


MrBFlt GetParsimonyLength (Tree *t, int chain)
{
    int             c, i, n, division;
    BitsLong        *pP, *pA;
    CLFlt           *nSitesOfPat;
    MrBFlt          length;
    TreeNode        *p;
    ModelInfo       *m;

    /* Get length down to internal root */
    length = GetParsDP(t, t->root->left, chain);

    /* return if rooted */
    if (t->isRooted == NO)
        return length;

    /* add in the last terminal if not rooted */
    for (n=0; n<t->nRelParts; n++)      
        {
        division = t->relParts[n];

        /* Find model settings */
        m = &modelSettings[division];

        /* Find number of site patterns */
        nSitesOfPat = numSitesOfPat + ((chainId[chain] % chainParams.numChains) * numCompressedChars) + m->compCharStart;

        /* Deal with last branch */
        p = t->intDownPass[t->nIntNodes-1];

        /* find downpass parsimony sets for the node and its environment */
        pP    = m->parsSets[p->index     ];
        pA    = m->parsSets[p->anc->index];
        
        if (m->nParsIntsPerSite == 1)
            {
            for (c=0; c<m->numChars; c++)
                {
                if ((pP[c] & pA[c]) == 0)
                    {
                    length += nSitesOfPat[c];
                    }
                }
            }
        else /* if (m->nParsIntsPerSite == 2) */
            {
            for (c=i=0; c<m->numChars; c++, i+=2)
                {
                if ((pP[i] & pA[i]) == 0 && (pP[i+1] & pA[i+1]) == 0)
                    {
                    length += nSitesOfPat[c];
                    }
                }
            }
        }

    return length;
}


void GetParsimonySubtreeRootstate (Tree *t, TreeNode *root, int chain)
{
    int             c, i, n, division;
    BitsLong        *pD, *pP, *pA, x[2];
    TreeNode        *p;
    ModelInfo       *m;

    /* Loop over divisions */
    for (n=0; n<t->nRelParts; n++)
        {
        division = t->relParts[n];
            
        /* Find model settings */
        m = &modelSettings[division];
        assert (m->nParsIntsPerSite == 1 || m->nParsIntsPerSite == 2);

        for (i=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            p->marked = NO;
            }

        p = root;
        while (p->anc != NULL)
            {
            p->marked = YES;
            p = p->anc;
            }

        /* Make uppass node by node */
        for (i=t->nIntNodes-1; i>=0; i--)
            {
            p = t->intDownPass[i];

            /* continue if no work needs to be done */
            if (p->marked == NO)
                continue;

            /* find downpass and uppass parsimony sets for the node and its environment */
            pP     = m->parsSets[p->index       ];
            if (p->left->marked == YES)
                pD = m->parsSets[p->right->index];
            else
                pD = m->parsSets[p->left->index ];
            pA     = m->parsSets[p->anc->index  ];
            
            if (m->nParsIntsPerSite == 1)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    x[0] = pD[c] & pA[c];
                    if (x[0] == 0)
                        {
                        x[0] = (pD[c] | pA[c]);
                        }
                    pP[c] = x[0];
                    }
                }
            else if (m->nParsIntsPerSite == 2)
                {
                for (c=i=0; c<m->numChars; c++, i+=2)
                    {
                    x[0] = pD[i  ] & pA[i  ];
                    x[1] = pD[i+1] & pA[i+1];
                    if (x[0] + x[1] == 0)
                        {
                        x[0] = (pD[i  ] | pA[i  ]);
                        x[1] = (pD[i+1] | pA[i+1]);
                        }
                    pP[i  ] = x[0];
                    pP[i+1] = x[1];
                    }
                }
            if (p == root)
                break;
            }
        }

    return;
    MrBayesPrint ("%d", chain); /* just because I am tired of seeing the unused parameter error msg */
}


/* GetRate: retrieve the base rate for the division and chain in current state */
MrBFlt GetRate (int division, int chain)
{
    Param   *p;
    MrBFlt  *values, rate;
    int     i;

    rate = 0.0;

    p = modelSettings[division].rateMult;
    values = GetParamVals (p, chain, state[chain]);
    if (p->nValues == 1)
        rate = values[0];
    else
        {
        for (i=0; i<p->nRelParts; i++)
            {
            if (p->relParts[i] == division)
                {
                rate = values[i];
                break;
                }
            }
        }

    p = modelSettings[division].geneTreeRateMult;
    if (p != NULL)
        {
        values = GetParamVals (p, chain, state[chain]);
        for (i=0; i<p->nRelParts; i++)
            {
            if (p->relParts[i] == division)
                {
                rate *= values[i];
                break;
                }
            }
        }
    
    return rate;
}


void GetStamp (void)
{
    int     i;

    for (i=0; i<10; i++)
        stamp[i] = '0' + (int)(RandomNumber(&runIDSeed) * 10);
    stamp[10] = '\0';

    MrBayesPrint ("%s   MCMC stamp = %s\n", spacer, stamp);
}


void GetSwappers (int *swapA, int *swapB, int run)
{
    int         i;
    
    /* this works for both the serial and parallel versions because the swapSeed is identical for all
        processors, ensuring they all get the same sequence of chainIds to swap */
#   if defined (MPI_ENABLED)

    /* For now, we wonly allow random swaps in the MPI version. Other schemes require
       tagging of messages, or a dedicated server node doing message processing.      */
    (*swapA) = (int) (RandomNumber(&swapSeed) * chainParams.numChains);
    (*swapB) = (int) (RandomNumber(&swapSeed) * (chainParams.numChains - 1));
    if ((*swapB) == (*swapA))
        (*swapB) = chainParams.numChains - 1;

#   else

    if (chainParams.swapAdjacentOnly == NO)
        {
        (*swapA) = (int) (RandomNumber(&swapSeed) * chainParams.numChains);
        (*swapB) = (int) (RandomNumber(&swapSeed) * (chainParams.numChains - 1));
        if ((*swapB) == (*swapA))
            (*swapB) = chainParams.numChains - 1;
        }
    else
        {
        (*swapA) = (int) (RandomNumber(&swapSeed) * (chainParams.numChains - 1));
        (*swapB) = (*swapA) + 1;
        }
#   endif

    i = run * chainParams.numChains;
    (*swapA) += i;
    (*swapB) += i;

    return;
}


void GetTempDownPassSeq (TreeNode *p, int *i, TreeNode **dp)
{
    if (p != NULL)
        {
        GetTempDownPassSeq (p->left,  i, dp);
        GetTempDownPassSeq (p->right, i, dp);
        dp[(*i)++] = p;
        }
}


MrBFlt GibbsSampleGamma (int chain, int division, RandLong *seed)
{
    int             c, i, k, *rateCat, nStates, nRateCats, nPosRateCats, id;
    CLFlt           **catLike, **catLnScaler, *lnScaler, maxLnScaler,
                    *clRoot, f, bs[64], *clInvar, pInvar, freq;
    MrBFlt          ran, lnL, *bsVals, deltaLnL, temp;
    ModelInfo       *m;
    Tree            *t;
    TreeNode        *p;

    m = &modelSettings[division];

    /* find base frequencies */
    bsVals = GetParamSubVals (m->stateFreq, chain, state[chain]);
    nStates = m->numModelStates;
    for (i=0; i<nStates; i++)
        bs[i] = (CLFlt) bsVals[i];

    /* find tree scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find category like array and associated sccaler */
    catLike = m->catLike[chain];
    catLnScaler = m->catLnScaler[chain];
    
    /* find rate category index */
    rateCat = m->tiIndex + chain*m->numChars;
    
    /* find number of rate cats and gamma cats */
    nRateCats = nPosRateCats = m->numRateCats;
    if (m->pInvar != NULL)
        nRateCats++;

    /* find tree */
    t = GetTree (m->brlens, chain, state[chain]);

    /* find invar cond likes (if we have invariable rate category) */
    clInvar = m->invCondLikes;

    /* get pInvar */
    if (m->pInvar == NULL)
        pInvar = 0.0;
    else
        pInvar = (CLFlt) *GetParamVals (m->pInvar, chain, state[chain]);
    freq = ((CLFlt)1.0 - pInvar) / nPosRateCats;

    /* get chain temperature */
    temp = Temperature (chainId[chain]);
    id = chainId[chain] % chainParams.numChains;

    /* calculate rate probs */
    for (k=0; k<nPosRateCats; k++)
        {
        FlipSiteScalerSpace(m, chain);
        ResetSiteScalers(m, chain);
        for (c=0; c<m->numChars; c++)
            rateCat[c] = k;
        for (i=0; i<t->nIntNodes; i++)
            {
            p = t->intDownPass[i];
            if (t->isRooted == NO && p->anc->anc == NULL)
                m->CondLikeRoot (p, division, chain);
            else
                m->CondLikeDown (p, division, chain);
            if (m->unscaledNodes[chain][p->index] == 0)
                m->CondLikeScaler (p, division, chain);
            }
        /* find root conditional likes */
        p = t->root->left;
        clRoot = m->condLikes[m->condLikeIndex[chain][p->index]];
        for (c=0; c<m->numChars; c++)
            {
            catLike[k][c] = 0.0;
            for (i=0; i<nStates; i++)
                catLike[k][c] += bs[i]*clRoot[i];
            catLike[k][c] *= freq;
            catLnScaler[k][c] = lnScaler[c];
            clRoot += nStates;
            }
        FlipSiteScalerSpace(m, chain);
        }

    /* fill in the invar cond likes if needed */
    if (m->pInvar != NULL)
        {
        k = nPosRateCats;
        for (c=0; c<m->numChars; c++)
            {
            catLike[k][c] = 0.0;
            for (i=0; i<nStates; i++)
                catLike[k][c] += bs[i]*clInvar[i];
            catLike[k][c] *= pInvar;
            clInvar += nStates;
            catLnScaler[k][c] = 0.0;
            }
        }

    /* Now Gibbs sample the rate categories */
    for (c=0; c<m->numChars; c++)
        {
        /* find max scaler */
        maxLnScaler = catLnScaler[0][c];
        for (k=1; k<nRateCats; k++)
            {
            if (catLnScaler[k][c] > maxLnScaler && catLike[k][c] > 0.0)
                maxLnScaler = catLnScaler[k][c];
            }
        /* scale values */
        for (k=0; k<nRateCats; k++)
            {
            f = catLnScaler[k][c] - maxLnScaler;
            if (f < -100.0)
                catLike[k][c] = 0.0;
            else
                {
                catLike[k][c] *= (CLFlt) exp (f);
                /* take the temperature into account */
                if (id != 0)
                    catLike[k][c] = (CLFlt) pow(catLike[k][c], temp);
                }
            }
        /* get cumulative sum */
        for (k=1; k<nRateCats; k++)
            catLike[k][c] += catLike[k-1][c];
        /* randomly sample a category; multiply by total to avoid scaling probs */
        ran = RandomNumber(seed) * catLike[nRateCats-1][c];
        for (k=0; k<nRateCats; k++)
            {
            if (ran < catLike[k][c])
                break;
            }
        rateCat[c] = k;
        }

    /* recalculate everything */
    FlipSiteScalerSpace(m, chain);
    ResetSiteScalers(m, chain);
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        if (t->isRooted == NO && p->anc->anc == NULL)
            m->CondLikeRoot (p, division, chain);
        else
            m->CondLikeDown (p, division, chain);
        if (m->unscaledNodes[chain][p->index] == 0)
            m->CondLikeScaler (p, division, chain);
        }
    lnL = 0.0;
    m->Likelihood (t->root->left, division, chain, &lnL, (chainId[chain] % chainParams.numChains));
    
    deltaLnL = lnL - m->lnLike[2*chain + state[chain]];
    m->lnLike[2*chain + state[chain]] =  lnL;

    return (deltaLnL);
}


/*------------------------------------------------------------------------
|
|   InitAdGamma: initialize variables for adgamma model.
|
-------------------------------------------------------------------------*/
int InitAdGamma (void)
{
    int         d, c, i, j, k, maxRates, *corrModel;
    ModelInfo   *m;
    
    /* take care of adgamma model */
    if (chainHasAdgamma == NO)
        return (NO_ERROR);
    
    MrBayesPrint ("%s   Initializing autocorrelated discrete gamma model\n", spacer);
    
    /* allocate corr space */
    corrModel = (int *) SafeCalloc (numCurrentDivisions, sizeof(int));
    if (!corrModel)
        return ERROR;
    
    /* allocate siteJump */
    if (memAllocs[ALLOC_SITEJUMP] == YES)
        {
        MrBayesPrint ("%s   siteJump not free in InitAdGamma\n", spacer);
        free (corrModel);
        return ERROR;
        }
    siteJump = (int *) SafeCalloc (numChar, sizeof(int));
    if (siteJump)
        memAllocs[ALLOC_SITEJUMP] = YES;
    else
        {
        MrBayesPrint ("%s   Problem allocating siteJump in InitAdGamma (%d ints)\n", spacer, numChar);
        free (corrModel);
        return ERROR;
        }
    
    /* reset vector indicating the matrices needed */
    for (i=0; i<MAX_SMALL_JUMP; i++)
        hasMarkovTi[i] = NO;
    
    /* fill in siteJump */
    for (i=0; i<numCurrentDivisions; i++)
        corrModel[i] = 0;

    k = 1;  /* index to corr model, 0 means no corr model */
    maxRates = 0;   /* max no. rates */
    for (d=0; d<numCurrentDivisions; d++)
        modelSettings[d].mark = NO;

    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        
        if (m->correlation == NULL || m->mark == YES)
            continue;
        
        m->mark = YES;
        for (i=0; i<m->correlation->nRelParts; i++)
            {
            if (modelSettings[m->correlation->relParts[i]].shape == 
                modelSettings[d].shape)
                {
                modelSettings[m->correlation->relParts[i]].mark = YES;
                corrModel[m->correlation->relParts[i]] = k;
                }
            }
        k++;

        if (m->numRateCats > maxRates)
            maxRates = m->numRateCats;

        }

    for (c=0; c<numChar; c++)
        {
        if (charInfo[c].isExcluded == YES)
            continue;
        
        if ((k=corrModel[partitionId[c][partitionNum] - 1]) == 0)
            continue;

        /* How far back is last char in this HMM? */
        for (j=c-1; j>=0; j--)
            {
            if (corrModel[partitionId[j][partitionNum] - 1] == k)
                break;
            }

        if (j<0)
            siteJump[c] = 0;
        else if (charInfo[j].bigBreakAfter == YES)
            siteJump[c] = BIG_JUMP;
        else
            {
            siteJump[c] = c - j;
            hasMarkovTi[c-j-1] = YES;
            }
        }

    /* check if any HMM is empty */
    k=0;
    for (i=0; i<numCurrentDivisions; i++)
        {
        if (corrModel[i] > k)
            k = corrModel[i];
        }
    for (i=1; i<=k; i++)
        {
        for (c=j=0; c<numChar; c++)
            {
            if (charInfo[c].isExcluded == NO && corrModel[partitionId[c][partitionNum] - 1] == i)
                j = c;
            }
        if (j == 0)
            {
            MrBayesPrint ("%s   ERROR: HMM model %d is empty.\n",spacer,i);
            free (corrModel);
            return (ERROR);
            }
        }

    /* allocate MarkovTis (space needed for calculations) */
    if (memAllocs[ALLOC_MARKOVTIS] == YES)
        {
        MrBayesPrint ("%s   markovTis not free in InitAdGamma\n", spacer);
        free (corrModel);
        return ERROR;
        }

    for (i=0; i<MAX_SMALL_JUMP; i++)
        {
        if (hasMarkovTi[i] == YES || i == 0)    /* base matrix always needed */
            {
            markovTi[i] = AllocateSquareDoubleMatrix(maxRates);
            if (markovTi[i] == NULL)
                break;
            }
        else
            markovTi[i] = NULL;
        }

    markovTiN = AllocateSquareDoubleMatrix(maxRates);
    if (i >= MAX_SMALL_JUMP && markovTiN)
        memAllocs[ALLOC_MARKOVTIS] = YES;
    else
        {
        MrBayesPrint ("%s   Problem allocating MarkovTis in InitAdGamma (%d MrBFlt)\n", spacer, 2 * MAX_RATE_CATS * MAX_RATE_CATS);
        for (i=0; i<MAX_SMALL_JUMP; i++)
            if (markovTi[i] != NULL) 
                FreeSquareDoubleMatrix (markovTi[i]);
        if (markovTiN != NULL) 
            FreeSquareDoubleMatrix(markovTiN);
        free (corrModel);
        return ERROR;
        }

    /* allocate space for rateProbs needed by adgamma model */

    /* calculate size needed */
    i = 0;
    for (j=0; j<numCurrentDivisions; j++)
        {
        m = &modelSettings[j];
        if (m->correlation != NULL)
            {
            m->rateProbStart = i;
            i += m->numRateCats * m->numChars;
            }
        }
    rateProbRowSize = i;

    /* allocate space */
    if (memAllocs[ALLOC_RATEPROBS] == YES)
        {
        MrBayesPrint ("%s   Space for rate probs not free in InitAdGamma\n", spacer);
        free (corrModel);
        return ERROR;
        }
    rateProbSpace = (MrBFlt *) SafeMalloc (2 * numLocalChains * rateProbRowSize * sizeof(MrBFlt));
    rateProbs = (MrBFlt **) SafeMalloc (numLocalChains * sizeof(MrBFlt *));
    if (!rateProbSpace || !rateProbs)
        {
        MrBayesPrint ("%s   Problem allocating rate probs\n", spacer);
        if (rateProbSpace) 
            free (rateProbSpace);
        if (rateProbs) 
            free (rateProbs);
        free (corrModel);
        return ERROR;
        }
    else
        memAllocs[ALLOC_RATEPROBS] = YES;

    /* set chain rateProbs pointers */
    for (i=j=0; i<numLocalChains; i++)
        {
        rateProbs[i] = rateProbSpace + j;
        j += 2 * rateProbRowSize;
        }

    free (corrModel);

    return (NO_ERROR);
}


/*------------------------------------------------------------------------
|
|   InitAugmentedModels: allocate and initialize space for augmented
|      models
|
-------------------------------------------------------------------------*/
int InitAugmentedModels (void)
{
    int         d, i, j, useAugmentedModels, nRateCats;
    ModelInfo   *m;

    useAugmentedModels = NO;
    for (d=0; d<numCurrentDivisions; d++)
        {
        if (modelSettings[d].gibbsGamma == YES)
                useAugmentedModels = YES;
        }
    
    if (useAugmentedModels == NO)
        return (NO_ERROR);
    
    MrBayesPrint ("%s   Initializing variables for model augmentation\n", spacer);
    
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        if (m->gibbsGamma == NO)
                continue;
        m->tiIndex = (int *) SafeMalloc (numLocalChains * m->numChars * sizeof (int));
        if (!m->tiIndex)
            return ERROR;
        m->catLike = (CLFlt ***) SafeMalloc (numLocalChains * sizeof (CLFlt **));
        if (!m->catLike)
            return ERROR;
        m->catLnScaler = (CLFlt ***) SafeMalloc (numLocalChains * sizeof (CLFlt **));
        if (!m->catLnScaler)
            return ERROR;
        if (m->pInvar == NULL)
            nRateCats = m->numRateCats;
        else
            nRateCats = m->numRateCats + 1;
        for (i=0; i<numLocalChains; i++)
            {
            m->catLike[i] = (CLFlt **) SafeCalloc (nRateCats, sizeof (CLFlt *));
            if (!m->catLike[i])
                return ERROR;
            m->catLnScaler[i] = (CLFlt **) SafeCalloc (nRateCats, sizeof (CLFlt *));
            if (!m->catLnScaler[i])
                return ERROR;
            for (j=0; j<nRateCats; j++)
                {
                m->catLike[i][j] = (CLFlt *) SafeCalloc (m->numChars, sizeof (CLFlt));
                if (!m->catLike[i][j])
                    return ERROR;
                m->catLnScaler[i][j] = (CLFlt *) SafeCalloc (m->numChars, sizeof (CLFlt));
                if (!m->catLnScaler[i][j])
                    return ERROR;
                }
            }
        }

    return NO_ERROR;
}


/*------------------------------------------------------------------------
|
|   InitChainCondLikes: (1) calculate size of cond like arrays, tiprob arrays
|       and scaler arrays; (2) allocate space for cond like, tiprob and
|       scaler arrays; (3) allocate and set node indices pointing to
|       cond like and scaler arrays; (4) initialize tip cond likes;
|       (5) allocate space for precalculated cond likes; (6) allocate
|       space for adgamma probs, if used.
|
-------------------------------------------------------------------------*/
int InitChainCondLikes (void)
{
    int         c, d, i, j, k, s, t, numReps, condLikesUsed, nIntNodes, nNodes, useBeagle,
                clIndex, tiIndex, scalerIndex, indexStep;
    BitsLong    *charBits;
    CLFlt       *cL;
    ModelInfo   *m;
#   if defined (SSE_ENABLED)
    int         j1;
#   endif
#   if defined (BEAGLE_ENABLED)
    double      *nSitesOfPat;
    MrBFlt      freq;
#   endif

    /* Figure out how large cond like array is needed, and how many cond like, scaler and tiprob arrays are needed.
       Also check for possible use of Beagle */
    condLikesUsed = NO;
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];

        MrBayesPrint ("%s   Division %d has %d unique site patterns\n", spacer, d+1, m->numChars);

        /* initialize model settings for chain cond likes */
        m->condLikeLength = 0;
        m->numCondLikes = 0;

        if (m->parsModelId == YES)
            continue;

        condLikesUsed = YES;

        /* figure out length of cond like array */
        if (m->dataType == STANDARD)
            {
#   if defined (BEAGLE_ENABLED)
            m->useBeagle = NO;
#   endif
            for (c=0; c<m->numChars; c++)
                {
                numReps = m->numRateCats;
                if (m->nStates[c] == 2)
                    numReps *= m->numBetaCats;
                m->condLikeLength += m->nStates[c] * numReps;
                }
            }
        else
            {
            if (m->gibbsGamma == YES)
                m->condLikeLength = m->numChars * m->numModelStates;
            else
                m->condLikeLength = m->numChars * m->numRateCats * m->numOmegaCats * m->numModelStates;
#   if defined (BEAGLE_ENABLED)
            /* tentatively decide on whether to use Beagle */
            if (tryToUseBEAGLE == YES)
                {
                if (m->printAncStates == YES || m->printSiteRates == YES || m->printPosSel ==YES || m->printSiteOmegas==YES)
                    {
                    MrBayesPrint ("%s   Non-beagle version of conditional likelihood calculator will be used for division %d due to\n", spacer, d+1);
                    MrBayesPrint ("%s   request of reporting 'ancestral states', 'site rates', 'pos selection' or 'site omegas'.\n", spacer);
                    }                
                else if (m->gibbsGamma == NO)
                    m->useBeagle = YES;
                }
#   endif
            }
        
        /* find size of tree */
        nIntNodes = GetTree(m->brlens, 0, 0)->nIntNodes;
        nNodes = GetTree(m->brlens, 0, 0)->nNodes;

        /* figure out number of cond like arrays */
        m->numCondLikes = (numLocalChains + 1) * (nIntNodes);
        m->numCondLikes += numLocalTaxa;
        /*
#   if !defined (DEBUG_NOSHORTCUTS)
        for (i=0; i<numLocalTaxa; i++)
            {
            if (m->isPartAmbig[i] == NO && m->dataType != STANDARD)
                m->numCondLikes--;
            }
#   endif
        */

        /* figure out number of node and site scalers */
        m->numScalers = (numLocalChains + 1) * (nIntNodes + 1);   /* add 1 for site scalers */

        /* figure out length of ti prob array and number of ti prob arrays */
        m->tiProbLength = 0;
        if (m->dataType == STANDARD)
            {
            m->numTiCats = 0;   /* We do not have repeated similar transition probability matrices */
            if (m->stateFreq->paramId == SYMPI_EQUAL)
                {
                for (k=0; k<9; k++)
                    {
                    if (m->isTiNeeded[k] == YES)
                        m->tiProbLength += (k + 2) * (k + 2) * m->numRateCats;
                    }
                for (k=9; k<13; k++)
                    {
                    if (m->isTiNeeded[k] == YES)
                        m->tiProbLength += (k - 6) * (k - 6) * m->numRateCats;
                    }
                for (k=13; k<18; k++)
                    {
                    if (m->isTiNeeded[k] == YES)
                         m->tiProbLength += (k - 11) * (k - 11) * m->numRateCats;
                    }
                }
            else
                {
                /* deal with unequal state frequencies */
                if (m->isTiNeeded[0] == YES)
                    m->tiProbLength += 4 * m->numRateCats * m->numBetaCats;
                for (c=0; c<m->numChars; c++)
                    {
                    if (m->nStates[c] > 2 && (m->cType[c] == UNORD || m->cType[c] == ORD))
                        {
                        m->tiProbLength += (m->nStates[c] * m->nStates[c]) * m->numRateCats;
                        }
                    }
                }
            }
        else
            {
            m->numTiCats    = m->numRateCats * m->numBetaCats * m->numOmegaCats;   /* A single partition has either gamma, beta or omega categories */
            m->tiProbLength = m->numModelStates * m->numModelStates * m->numTiCats;
            }
        m->numTiProbs = (numLocalChains + 1) * nNodes;
        
        /* set info about eigen systems */
        if (InitEigenSystemInfo (m) == ERROR)
            return (ERROR);
        }

    /* check if conditional likelihoods are needed */
    if (condLikesUsed == YES)
        MrBayesPrint ("%s   Initializing conditional likelihoods\n", spacer);
    else
        return NO_ERROR;

    /* allocate space and fill in info for tips */
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
       
        /* allocate space for conditional likelihoods */
        useBeagle = NO;
#   if defined (BEAGLE_ENABLED)
        if (m->useBeagle == YES)
            {
            if (InitBeagleInstance(m, d) != ERROR)
                useBeagle = YES;
            else
                m->useBeagle = NO;
            }
#   endif
#   if defined (SSE_ENABLED)
        /*if (useBeagle == NO && m->dataType != STANDARD)
            m->useSSE = YES;*/
        if (useBeagle == YES)
            {
            m->useVec = VEC_NONE;
            m->numFloatsPerVec = 0;
            }

#   endif
        if (m->useBeagle == NO && m->useVec == VEC_NONE)
            MrBayesPrint ("%s   Using standard non-SSE likelihood calculator for division %d (%s-precision)\n", spacer, d+1, (sizeof(CLFlt) == 4 ? "single" : "double"));
        else if (m->useBeagle == NO && m->useVec == VEC_SSE)
            MrBayesPrint ("%s   Using standard SSE likelihood calculator for division %d (single-precision)\n", spacer, d+1);
        else if (m->useBeagle == NO && m->useVec == VEC_AVX)
            MrBayesPrint ("%s   Using standard AVX likelihood calculator for division %d (single-precision)\n", spacer, d+1);
        else if (m->useBeagle == NO && m->useVec == VEC_FMA)
            MrBayesPrint ("%s   Using standard FMA likelihood calculator for division %d (single-precision)\n", spacer, d+1);
        else if (m->useBeagle == NO)
            {
            MrBayesPrint ("%s   WARNING! Using unknown SIMD likelihood calculator for division %d (single-precision)\n", spacer, d+1);
            return (ERROR);
            }

        if (useBeagle == NO)
            {
            /* allocate cond like space */
            m->condLikes = (CLFlt**) SafeMalloc(m->numCondLikes * sizeof(CLFlt*));
            if (!m->condLikes)
                return (ERROR);
            for (i=0; i<m->numCondLikes; i++)
                {
#   if defined (SSE_ENABLED)
                if (m->useVec != VEC_NONE)
                    {
                    /* calculate number SSE chars */
                    m->numVecChars = ((m->numChars - 1) / m->numFloatsPerVec) + 1;

                    /* allocate space with padding (m->condLikeLength is without padding) */
                    if (m->gibbsGamma == YES)
                        numReps = 1;
                    else
                        numReps = m->numRateCats * m->numOmegaCats;
                    k = m->numVecChars * m->numFloatsPerVec * m->numModelStates * numReps;
                    
                    if (m->useVec == VEC_AVX)
                        m->condLikes[i] = (CLFlt*) AlignedMalloc (k * sizeof(CLFlt), 32);
                    else
                        m->condLikes[i] = (CLFlt*) AlignedMalloc (k * sizeof(CLFlt), 16);

                    if (!m->condLikes[i])
                        return (ERROR);

                    /* start by filling all with 0.0f; pad when filling in tips */
                    for (j=0; j<k; j++)
                        m->condLikes[i][j] = 0.0f;
                    }
                else
                    {
                    m->condLikes[i] = (CLFlt*) SafeMalloc(m->condLikeLength * sizeof(CLFlt));
                    if (!m->condLikes[i])
                        return (ERROR);
                    }
#   else
                m->condLikes[i] = (CLFlt*) SafeMalloc(m->condLikeLength * sizeof(CLFlt));
                if (!m->condLikes[i])
                    return (ERROR);
#   endif
                }

            /* allocate scaler space and pointers for scaling */
            m->scalers = (CLFlt**) SafeMalloc(m->numScalers * sizeof(CLFlt*));
            if (!m->scalers)
                return (ERROR);
            for (i=0; i<m->numScalers; i++)
                {
#   if defined (SSE_ENABLED)
                if (m->useVec != VEC_NONE)
                    {
                    /* allocate space with padding */
                    if (m->useVec == VEC_SSE)
                        m->scalers[i] = (CLFlt*) AlignedMalloc (m->numVecChars * m->numFloatsPerVec * sizeof(CLFlt), 16);
                    else if (m->useVec == VEC_AVX || m->useVec == VEC_FMA)
                        m->scalers[i] = (CLFlt*) AlignedMalloc (m->numVecChars * m->numFloatsPerVec * sizeof(CLFlt), 32);
                    if (!m->scalers[i])
                        return (ERROR);
                    for (j=0; j<m->numVecChars*m->numFloatsPerVec; j++)
                        m->scalers[i][j] = 0.0f;
                    }
                else
                    {
                    m->scalers[i] = (CLFlt*) SafeMalloc (m->numChars * sizeof(CLFlt));
                    if (!m->scalers[i])
                        return (ERROR);
                    }
#   else
                m->scalers[i] = (CLFlt*) SafeMalloc (m->numChars * sizeof(CLFlt));
                if (!m->scalers[i])
                    return (ERROR);
#   endif
                }

            /* allocate stuff for facilitating scaling and accumulation of cond likes */
            if (m->dataType == STANDARD)
                {
                m->clP = (CLFlt **) SafeMalloc(m->numRateCats * sizeof(CLFlt *));
                if (!m->clP)
                    return (ERROR);
                }
            else
                {
                m->clP = (CLFlt **) SafeMalloc(m->numTiCats * sizeof(CLFlt *));
                if (!m->clP)
                    return (ERROR);
#   if defined (SSE_ENABLED)
                if (m->useVec == VEC_SSE)
                    {
                    m->clP_SSE = (__m128 **) SafeMalloc(m->numTiCats * sizeof(__m128 *));
                    if (!m->clP_SSE)
                        return (ERROR);
                    m->lnL_Vec  = AlignedMalloc (m->numVecChars * m->numFloatsPerVec * sizeof(CLFlt*), 16);
                    m->lnLI_Vec = AlignedMalloc (m->numVecChars * m->numFloatsPerVec * sizeof(CLFlt*), 16);
                    if (!m->lnL_Vec || !m->lnLI_Vec)
                        return (ERROR);
                    }
#if defined (AVX_ENABLED)
                else if (m->useVec == VEC_AVX)
                    {
                    m->clP_AVX = (__m256 **) SafeMalloc(m->numTiCats * sizeof(__m256 *));
                    if (!m->clP_AVX)
                        return (ERROR);
                    m->lnL_Vec  = AlignedMalloc (m->numVecChars * m->numFloatsPerVec * sizeof(CLFlt*), 32);
                    m->lnLI_Vec = AlignedMalloc (m->numVecChars * m->numFloatsPerVec * sizeof(CLFlt*), 32);
                    if (!m->lnL_Vec || !m->lnLI_Vec)
                        return (ERROR);
                    }
#endif
#if defined (FMA_ENABLED)
                else if (m->useVec == VEC_FMA)
                    {
                    m->clP_AVX = (__m256 **) SafeMalloc(m->numTiCats * sizeof(__m256 *));
                    if (!m->clP_AVX)
                        return (ERROR);
                    m->lnL_Vec  = AlignedMalloc (m->numVecChars * m->numFloatsPerVec * sizeof(CLFlt*), 32);
                    m->lnLI_Vec = AlignedMalloc (m->numVecChars * m->numFloatsPerVec * sizeof(CLFlt*), 32);
                    if (!m->lnL_Vec || !m->lnLI_Vec)
                        return (ERROR);
                    }
#endif
#   endif
                }

            /* allocate tiprob space */
            m->tiProbs = (CLFlt**) SafeMalloc(m->numTiProbs * sizeof(CLFlt*));
            if (!m->tiProbs)
                return (ERROR);
            for (i=0; i<m->numTiProbs; i++)
                {
                m->tiProbs[i] = (CLFlt*) SafeMalloc(m->tiProbLength * sizeof(CLFlt));
                if (!m->tiProbs[i])
                    return (ERROR);
                }
            }

        /* allocate eigen system space (needed also for Beagle version */
        if (m->nCijkParts > 0)
            {
            m->cijks = (MrBFlt**) SafeMalloc((numLocalChains + 1) * sizeof(MrBFlt*));
            if (!m->cijks)
                return ERROR;
            for (i=0; i<numLocalChains+1; i++)
                {
                m->cijks[i] = (MrBFlt*) SafeMalloc(m->cijkLength * sizeof(MrBFlt));
                if (!m->cijks[i])
                    return (ERROR);
                }
            }

        /* get size of tree */
        nIntNodes = GetTree(m->brlens,0,0)->nIntNodes;
        nNodes = GetTree(m->brlens,0,0)->nNodes;

            
        /* allocate and set indices from tree nodes to cond like arrays */
        m->condLikeIndex = (int **) SafeMalloc (numLocalChains * sizeof(int *));
        if (!m->condLikeIndex)
            return (ERROR);
        for (i=0; i<numLocalChains; i++)
            {
            m->condLikeIndex[i] = (int *) SafeMalloc (nNodes * sizeof(int));
            if (!m->condLikeIndex[i])
                return (ERROR);
            }
        for (i=0; i<numLocalChains; i++)
            for (j=0; j<nNodes; j++)
                m->condLikeIndex[i][j] = -1;

        /* set up indices for terminal nodes */
        clIndex = 0;
        if (useBeagle == YES)
            indexStep = m->nCijkParts;
        else
            indexStep = 1;
        for (i=0; i<numLocalTaxa; i++)
            {
#   if !defined (DEBUG_NOSHORTCUTS)
            /* TODO: Untill CondLikeRoot_XXX are fixed (case 4 when one of the children is non-ambig) we allocate space for non-ambig tips. if fixed also uncoment down the function */
            /* if (useBeagle == NO && useSSE == NO && m->isPartAmbig[i] == NO && m->dataType != STANDARD)
                continue;
            */
#   endif
            for (j=0; j<numLocalChains; j++)
                m->condLikeIndex[j][i] = clIndex;
            clIndex += 1; /* even for multiple omega cat we need only one set of conditional likelihoods  for terminals for all chains.*/
            }

        /* reserve private space for parsimony-based moves if parsimony model is used */
        if (m->parsModelId == YES && m->parsimonyBasedMove == YES)
            clIndex += nIntNodes;

        /* set up indices for internal nodes */
        for (j=0; j<numLocalChains; j++)
            {
            for (i=0; i<nIntNodes; i++)
                {
                m->condLikeIndex[j][i+numLocalTaxa] = clIndex;
                clIndex += indexStep;
                }
            }

        /* allocate and set up scratch cond like indices */
        m->condLikeScratchIndex = (int *) SafeMalloc (nNodes * sizeof(int));
        if (!m->condLikeScratchIndex)
            return (ERROR);
        for (i=0; i<nNodes; i++)
            m->condLikeScratchIndex[i] = -1;
        for (i=0; i<nIntNodes; i++)
            {
            m->condLikeScratchIndex[i+numLocalTaxa] = clIndex;
            clIndex += indexStep;
            }

        /* allocate and set indices from tree edges to ti prob arrays */
        m->tiProbsIndex = (int **) SafeMalloc (numLocalChains * sizeof(int *));
        if (!m->tiProbsIndex)
            return (ERROR);
        for (i=0; i<numLocalChains; i++)
            {
            m->tiProbsIndex[i] = (int *) SafeMalloc (nNodes * sizeof(int));
            if (!m->tiProbsIndex[i])
                return (ERROR);
            }

        /* set up indices for nodes */
        tiIndex = 0;
        for (i=0; i<numLocalChains; i++)
            {
            for (j=0; j<nNodes; j++)
                {
                m->tiProbsIndex[i][j] = tiIndex;
                tiIndex += indexStep;
                }
            }

        /* allocate and set up scratch transition prob indices */
        m->tiProbsScratchIndex = (int *) SafeMalloc (nNodes * sizeof(int));
        if (!m->tiProbsScratchIndex)
            return (ERROR);
        for (i=0; i<nNodes; i++)
            {
            m->tiProbsScratchIndex[i] = tiIndex;
            tiIndex += indexStep;
            }

        /* allocate and set up rescale frequency */
        m->rescaleFreq = (int*) SafeMalloc((numLocalChains) * sizeof(int));
        for (i=0; i<numLocalChains; ++i)
            {
            if (m->numModelStates == 4 )
                m->rescaleFreq[i] = 1;
            else
                m->rescaleFreq[i] = 1;
            }

        /* allocate and set up number of unscaled nodes + scratch space */
        m->unscaledNodes = (int **) SafeCalloc (numLocalChains, sizeof(int *));
        if (!m->unscaledNodes)
            return (ERROR);
        for (i=0; i<numLocalChains; i++)
            {
            m->unscaledNodes[i] = (int *) SafeCalloc (nNodes, sizeof(int));
            if (!m->unscaledNodes[i])
                return (ERROR);
            }
        m->unscaledNodesScratch = (int *) SafeCalloc (nNodes, sizeof (int));
        if (!m->unscaledNodesScratch)
            return (ERROR);
            
        /* allocate and set up node scaler indices */
        scalerIndex = 0;
        m->nodeScalerIndex = (int **) SafeMalloc (numLocalChains * sizeof(int *));
        if (!m->nodeScalerIndex)
            return (ERROR);
        for (i=0; i<numLocalChains; i++)
            {
            m->nodeScalerIndex[i] = (int *) SafeMalloc (nNodes * sizeof(int));
            if (!m->nodeScalerIndex[i])
                return (ERROR);
            for (j=0; j<nNodes; j++)
                m->nodeScalerIndex[i][j] = -1;
            for (j=0; j<nIntNodes; j++)
                {
                m->nodeScalerIndex[i][j+numLocalTaxa] = scalerIndex;
                scalerIndex += indexStep;
                }
            }
        m->nodeScalerScratchIndex = (int *) SafeMalloc (nNodes * sizeof (int));
        if (!m->nodeScalerScratchIndex)
            return (ERROR);
        for (i=0; i<nNodes; i++)
            m->nodeScalerScratchIndex[i] = -1;
        for (i=0; i<nIntNodes; i++)
            {
            m->nodeScalerScratchIndex[i+numLocalTaxa] = scalerIndex;
            scalerIndex += indexStep;
            }

        /* allocate and set up site scaler indices */
        m->siteScalerIndex = (int *) SafeMalloc (numLocalChains * sizeof(int));
        if (!m->siteScalerIndex)
            return (ERROR);
        for (i=0; i<numLocalChains; i++)
            {
            m->siteScalerIndex[i] = scalerIndex;
            scalerIndex += indexStep;
            }
        m->siteScalerScratchIndex = scalerIndex;
            
#   if defined (BEAGLE_ENABLED)
        /* used only with Beagle advanced dynamic rescaling where we set scaler nodes for each partition  */
        if (m->useBeagle == YES)
            {
            m->succesCount = (int*) SafeMalloc((numLocalChains) * sizeof(int));
            m->beagleComputeCount = (long *) SafeMalloc(sizeof(long) * numLocalChains);
            t=BEAGLE_RESCALE_FREQ/m->numModelStates;
            if (beagleFlags & BEAGLE_FLAG_PRECISION_DOUBLE) /*if double presition is used*/
                t*=BEAGLE_RESCALE_FREQ_DOUBLE;
            for (i=0; i<numLocalChains; i++)
               {
               m->rescaleFreq[i] = t;
               }
            m->isScalerNode = (int**) SafeMalloc((numLocalChains) * sizeof(int*));
            /* we will use m->isScalerNode[chain][node->index] to determine whether the node is scaled or not. We do it only for internal nodes whose indexes start from numLocalTaxa thus we skew the pointer */
            m->isScalerNodeScratch = (int*) SafeMalloc(nIntNodes * sizeof(int)) - numLocalTaxa; 
            assert (NO == 0); /* SafeMalloc set the allocated memmory to 0 while we need to set it to NO */
            for (i=0; i<numLocalChains; i++)
               {
               m->isScalerNode[i] = (int*) SafeMalloc(nIntNodes * sizeof(int)) - numLocalTaxa;
               }
            }
#   endif

        /* allocate and set up cijk indices */
        if (m->nCijkParts > 0)
            {
            m->cijkIndex = (int *) SafeMalloc (numLocalChains * sizeof(int));
            if (!m->cijkIndex)
                return (ERROR);
            for (i=0; i<numLocalChains; i++)
                m->cijkIndex[i] = i*indexStep;
            m->cijkScratchIndex = numLocalChains*indexStep;
            }

#   if defined (BEAGLE_ENABLED)
            /* Set up nSitesOfPat for Beagle */
            if (m->useBeagle == YES)
                {
                nSitesOfPat = (double *) SafeMalloc (m->numChars * sizeof(double));
                for (c=0; c<m->numChars; c++)
                    nSitesOfPat[c] = numSitesOfPat[m->compCharStart + c];
                beagleSetPatternWeights(m->beagleInstance,
                                        nSitesOfPat);
                free (nSitesOfPat);
                nSitesOfPat = NULL;

                /* find category frequencies */
                if (m->pInvar == NO)
                    {
                    freq =  1.0 /  m->numRateCats;
                    
                    /* set category frequencies in beagle instance */
                    if (m->numOmegaCats <= 1)
                        {
                        for (i=0; i<m->numRateCats; i++)
                            m->inWeights[i] = freq;
                        for (i=0; i< (numLocalChains); i++) {
                            beagleSetCategoryWeights(m->beagleInstance,
                                                     m->cijkIndex[i],
                                                     m->inWeights);
                            }
                        beagleSetCategoryWeights(m->beagleInstance,
                                                 m->cijkScratchIndex,
                                                 m->inWeights);
                        }
                    }
                
                /* Set up scalers for Beagle */
                for (i=0; i<m->numScalers*m->nCijkParts; i++)
                    beagleResetScaleFactors(m->beagleInstance, i);
                }
#   endif

        /* fill in tip conditional likelihoods */
        if (m->dataType == STANDARD)
            {
            clIndex = 0;
            for (i=0; i<numLocalTaxa; i++)
                {
                cL = m->condLikes[clIndex++];
                for (t=0; t<m->numRateCats;t++)
                    {
                    charBits = m->parsSets[i];
                    for (c=0; c<m->numChars; c++)
                        {                   
                        if (m->nStates[c] == 2)
                            numReps = m->numBetaCats;
                        else
                            numReps = 1;
                        for (k=0; k<numReps; k++)
                            {
                            for (s=0; s<m->nStates[c]; s++)
                                {
                                if (IsBitSet(s, charBits))
                                    (*cL) = 1.0;
                                cL++;
                                }
                            }
                        charBits += m->nParsIntsPerSite;
                        }
                    }
                }
            }
        else if (useBeagle == NO)
            {
            if (m->gibbsGamma == YES)
                numReps = m->numTiCats / m->numRateCats;
            else
                numReps = m->numTiCats;

            clIndex = 0;
            for (i=0; i<numLocalTaxa; i++)
                {
#   if !defined (DEBUG_NOSHORTCUTS) && !defined (SSE_ENABLED)
                /* TODO: Untill CondLikeRoot_XXX are fixed (case 4 when one of the children is non-ambig) we allocate space for non-ambig tips. if fixed also uncomment up the function */
                /* if (m->isPartAmbig[i] == NO && m->dataType != RESTRICTION)
                    continue;
                */
#   endif
                cL = m->condLikes[clIndex++];
#   if defined (SSE_ENABLED)    /* deal with setup for SIMD code */
                if (m->useVec != VEC_NONE)
                    {
                    for (k=0; k<numReps; k++)
                        {
                        charBits = m->parsSets[i];
                        for (c=0; c<m->numChars/m->numFloatsPerVec; c++)
                            {
                            for (j=0; j<m->numModelStates/m->numStates; j++)
                                {
                                for (s=0; s<m->numStates; s++)
                                    {
                                    for (j1=0; j1<m->numFloatsPerVec; j1++)
                                        {
                                        if (IsBitSet(s, charBits + j1*m->nParsIntsPerSite))
                                            (*cL) = 1.0;
                                        cL++;
                                        }
                                    }
                                }   
                            charBits += m->numFloatsPerVec * m->nParsIntsPerSite;
                            }
                        if (m->numChars % m->numFloatsPerVec > 0)
                            {
                            /* add last characters and padd */
                            for (j=0; j<m->numModelStates/m->numStates; j++)
                                {
                                for (s=0; s<m->numStates; s++)
                                    {
                                    for (j1=0; j1<m->numChars%m->numFloatsPerVec; j1++)
                                        {
                                        if (IsBitSet(s, charBits + j1*m->nParsIntsPerSite))
                                            (*cL) = 1.0;
                                        cL++;
                                        }
                                    for (; j1<m->numFloatsPerVec; j1++)
                                        {
                                        (*cL) = 1.0;
                                        cL++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                else
                    {
                    for (k=0; k<numReps; k++)
                        {
                        charBits = m->parsSets[i];
                        for (c=0; c<m->numChars; c++)
                            {
                            for (j=0; j<m->numModelStates/m->numStates; j++)
                                {
                                for (s=0; s<m->numStates; s++)
                                    {
                                    if (IsBitSet(s, charBits))
                                        (*cL) = 1.0;
                                    cL++;
                                    }
                                }
                                charBits += m->nParsIntsPerSite;
                            }
                        }
                    }

#   else
                for (k=0; k<numReps; k++)
                    {
                    charBits = m->parsSets[i];
                    for (c=0; c<m->numChars; c++)
                        {
                        for (j=0; j<m->numModelStates/m->numStates; j++)
                            {
                            for (s=0; s<m->numStates; s++)
                                {
                                if (IsBitSet(s, charBits))
                                    (*cL) = 1.0;
                                cL++;
                                }
                            }
                            charBits += m->nParsIntsPerSite;
                        }
                    }
#   endif
                }
            }

        if (m->printAncStates == YES)
            {
            m->ancStateCondLikes = (CLFlt *) SafeMalloc (m->condLikeLength * sizeof(CLFlt));
            if (!m->ancStateCondLikes)
                return (ERROR);
            }
        }
    /* allocate space for precalculated likelihoods */
    j = 0;
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        if (m->dataType == STANDARD || m->parsModelId == YES)
            continue;

        i = (m->numModelStates + 1) * m->numModelStates * m->numTiCats;
        if (i > j)
            j = i;
        }
    if (j > 0) /* don't bother allocating precalculated likelihoods if we only have parsimony model or morphological characters */
        {
        if (memAllocs[ALLOC_PRELIKES] == YES)
            {
            MrBayesPrint ("%s   Space for preLikes not free in InitChainCondLikes\n", spacer);
            return ERROR;
            }
        preLikeL = (CLFlt *) SafeMalloc (3 * j * sizeof(CLFlt));
        if (!preLikeL)
            {
            MrBayesPrint ("%s   Problem allocating preLikes\n", spacer);
            return ERROR;
            }
        memAllocs[ALLOC_PRELIKES] = YES;
        preLikeR = preLikeL + j;
        preLikeA = preLikeR + j;
        }

    return NO_ERROR;
}


/*------------------------------------------------------------------------
|
|   InitEigenSystemInfo: set info about eigen decompositions
|
-------------------------------------------------------------------------*/
int InitEigenSystemInfo (ModelInfo *m)
{
    int         ts;
    
    if (m->dataType == STANDARD)
        {
        /* dealt with in ProcessStdChars */
        return (NO_ERROR);
        }

    m->cijkLength = 0;
    m->nCijkParts = 0;
    if (m->dataType == PROTEIN)
        {
        ts = m->numModelStates;
        m->cijkLength = (ts * ts * ts) + (2 * ts);
        m->nCijkParts = 1;
        if (m->switchRates != NULL) /* covarion model */
            {
            m->cijkLength *= m->numRateCats;
            m->nCijkParts = m->numRateCats;
            }
        }
    else if (m->dataType == DNA || m->dataType == RNA)
        {
        if (m->nucModelId == NUCMODEL_4BY4)
            {
            if (m->switchRates==NULL && m->nst != 6 && m->nst != NST_MIXED)
                {
                m->cijkLength = 0;
                m->nCijkParts = 0;
#   if defined (BEAGLE_ENABLED)
                ts = m->numModelStates;
                m->cijkLength = (ts * ts * ts) + (2 * ts);
                m->nCijkParts = 1;
#   endif
                }
            else
                {
                ts = m->numModelStates;
                m->cijkLength = (ts * ts * ts) + (2 * ts);
                m->nCijkParts = 1;
                }
            if (m->switchRates != NULL)
                {
                m->cijkLength *= m->numRateCats;
                m->nCijkParts = m->numRateCats;
                }
            }
        else if (m->nucModelId == NUCMODEL_DOUBLET)
            {
            ts = m->numModelStates;
            m->cijkLength = (ts * ts * ts) + (2 * ts);
            m->nCijkParts = 1;
            }
        else if (m->nucModelId == NUCMODEL_CODON)
            {
            ts = m->numModelStates;
            m->cijkLength = (ts * ts * ts) + (2 * ts);
            m->cijkLength *= m->numOmegaCats;
            m->nCijkParts = m->numOmegaCats;
            }
        else
            {
            MrBayesPrint ("%s   ERROR: Something is wrong if you are here.\n", spacer);
            return ERROR;
            }
        }
#   if defined (BEAGLE_ENABLED)
    else if (m->dataType == RESTRICTION)
        {
                assert (m->numModelStates == 2);
                ts = 2;
                m->cijkLength = (ts * ts * ts) + (2 * ts);
                m->nCijkParts = 1;
        }
#   endif
    return (NO_ERROR);
}


/*------------------------------------------------------------------------
|
|   InitFinalStateCondLikes: allocate space for final conditional
|       likelihoods if needed
|
-------------------------------------------------------------------------*/
int InitFinalStateCondLikes (void)
{
    int         d;
    ModelInfo   *m;
    
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        if (m->printAncStates == YES)
            {
            m->ancStateCondLikes = (CLFlt *) SafeMalloc (m->condLikeLength * sizeof(CLFlt));
            if (!m->ancStateCondLikes)
                return ERROR;
            }
        }
    return (NO_ERROR);
}


/*------------------------------------------------------------------------
|
|   InitInvCondLikes: allocate and initialize invariable conditional
|       likelihoods if needed
|
|       NB! Fills in invariable cond likes for all hidden states; this
|       is convenient although some space is wasted
|
-------------------------------------------------------------------------*/
int InitInvCondLikes (void)
{
    int         c, d, i, s, isConstant, usingInvCondLikes;
    BitsLong    *charBits;
    CLFlt       *cI;
    ModelInfo   *m;

#   if defined (SSE_ENABLED)
    int         c1;
#   endif

    /* allocate space for invariable cond likes */
    usingInvCondLikes = NO;
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];

        if (m->pInvar == NULL)
            continue;

        usingInvCondLikes = YES;
#   if defined (SSE_ENABLED)
        if ( m->useVec != VEC_NONE )
            {
            c1 = m->numVecChars * m->numFloatsPerVec * m->numModelStates;
            if (m->useVec == VEC_AVX)
                m->invCondLikes = (CLFlt *) AlignedMalloc (c1 * sizeof(CLFlt), 32);
            else
                m->invCondLikes = (CLFlt *) AlignedMalloc (c1 * sizeof(CLFlt), 16);
            for (i=0; i<c1; i++)
                m->invCondLikes[i] = 0.0f;
            }
#   else
        m->invCondLikes = (CLFlt *) SafeMalloc (m->numChars * m->numModelStates * sizeof(CLFlt));
#   endif
        if (!m->invCondLikes)
            return ERROR;
        }
    
    if (usingInvCondLikes == NO)
        return NO_ERROR;
    
    MrBayesPrint ("%s   Initializing invariable-site conditional likelihoods\n", spacer);
        
    /* fill in invariable-site conditional likelihoods */
    for (d=0; d<numCurrentDivisions; d++)
        {

        m = &modelSettings[d];
        
        if (m->pInvar == NULL)
            continue;
        
        cI = m->invCondLikes;
        if (m->dataType == STANDARD)
            {
            for (c=0; c<m->numChars; c++)
                {
                for (s=0; s<m->nStates[c]; s++)
                    {
                    isConstant = YES;
                    for (i=0; i<numLocalTaxa; i++)
                        {
                        charBits = &m->parsSets[i][c*m->nParsIntsPerSite];
                        if (IsBitSet(s, charBits) == NO)
                            {
                            isConstant = NO;
                            break;
                            }
                        }
                    if (isConstant == YES)
                        *cI = 1.0;
                    else
                        *cI = 0.0;
                    cI++;
                    }
                }
            }
        else    /* all other models for which pInvar is applicable */
            {
            assert (m->numModelStates == m->numStates);
#   if defined (SSE_ENABLED)
             for (c=0; c<m->numChars/m->numFloatsPerVec; c++)
                {
                for (s=0; s<m->numModelStates; s++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++)
                        {
                        isConstant = YES;
                        //charBits = parsMatrix + m->parsMatrixStart + ((c * m->numFloatsPerVec) + c1) * m->nParsIntsPerSite;
                        for (i=0; i<numLocalTaxa; i++)
                            {
                            charBits = &m->parsSets[i][((c * m->numFloatsPerVec) + c1) *m->nParsIntsPerSite];
                            if (IsBitSet(s, charBits) == NO)
                                {
                                isConstant = NO;
                                break;
                                }
                            //charBits += parsMatrixRowSize;
                            }
                        if (isConstant == YES)
                            *cI = 1.0;
                        cI++;
                        }
                    }
                }
             if (m->numChars % m->numFloatsPerVec != 0)
                {
                for (s=0; s<m->numModelStates; s++)
                    {
                    for (c1=0; c1<m->numChars%m->numFloatsPerVec; c1++)
                        {
                        isConstant = YES;
                        //charBits = parsMatrix + m->parsMatrixStart + (((m->numChars / m->numFloatsPerVec) * m->numFloatsPerVec) + c1) * m->nParsIntsPerSite;
                        for (i=0; i<numLocalTaxa; i++)
                            {
                            charBits = &m->parsSets[i][(((m->numChars / m->numFloatsPerVec) * m->numFloatsPerVec) + c1) *m->nParsIntsPerSite];
                            if (IsBitSet(s, charBits) == NO)
                                {
                                isConstant = NO;
                                break;
                                }
                            //charBits += parsMatrixRowSize;
                            }
                        if (isConstant == YES)
                            *cI = 1.0;
                        cI++;
                        }
                    for (; c1<m->numFloatsPerVec; c1++)
                        {
                        *cI = 1.0;
                        cI++;
                        }
                    }
                }
#   else
            for (c=0; c<m->numChars; c++)
                {
                for (s=0; s<m->numModelStates; s++)
                    {
                    isConstant = YES;
                    for (i=0; i<numLocalTaxa; i++)
                        {
                        charBits = &m->parsSets[i][c*m->nParsIntsPerSite];
                        if (IsBitSet(s, charBits) == NO)
                            {
                            isConstant = NO;
                            break;
                            }
                        }
                    if (isConstant == YES)
                        *cI = 1.0;
                    cI++;
                    }
                }
#   endif
            }
        }   /* next division */

#   if 0
    invCondLikeSize = 0;
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        if (m->pInvar == NULL)
            continue;
        cI = m->invCondLikes;
        if (m->dataType == STANDARD)
            {
            }
        else
            {
            for (c=0; c<m->numChars; c++)
                {
                printf ("%4d -- ", c);
                for (s=0; s<m->numModelStates; s++)
                    {
                    printf ("%1.0lf", *cI);
                    cI++;
                    }
                printf ("\n");
                }
            }
        }
#   endif
        
    return NO_ERROR;
}


/*------------------------------------------------------------------------
|
|   InitParsSets: allocate space for and set parsimony state sets
|
-------------------------------------------------------------------------*/
int InitParsSets (void)
{
    int             c, i, j, k, d, nParsStatesForCont, nIntNodes, nNodes,
                    nuc1, nuc2, nuc3, codingNucCode, allNucCode;
    BitsLong        allAmbig, x, x1, x2, x3, *longPtr, bitsLongOne;
    ModelInfo       *m;
    ModelParams     *mp;

    bitsLongOne = 1;

    /* this variable determines how many parsimony states are used           */
    /* to represent continuous characters (determines weight of these chars) */
    nParsStatesForCont = 3;

    /* find number and size of parsimony sets and node lengths */
    for (d=0; d<numCurrentDivisions; d++)
        {
        m  = &modelSettings[d];
        mp = &modelParams[d];

        /* find how many parsimony ints (BitsLong) are needed for each model site */
        if (mp->dataType == CONTINUOUS)
            {
            /* scale continuous characters down to an ordered parsimony character */
            /* with nParsStatesForCont states, represent this character as a set */
            /* of binary characters by additive binary coding */
            m->nParsIntsPerSite = nParsStatesForCont - 1;
            }
        else
            m->nParsIntsPerSite = 1 + m->numStates / nBitsInALong;

        /* Calculate number of nodes and number of internal nodes */
        nIntNodes = GetTree(m->brlens,0,0)->nIntNodes;
        nNodes    = GetTree(m->brlens,0,0)->nNodes;
        
        /* Calculate number of parsimony sets */
        m->numParsSets = numLocalTaxa;
        if (m->parsimonyBasedMove == YES || !strcmp(chainParams.startTree, "Parsimony"))
            m->numParsSets += nIntNodes;
        if (m->parsModelId == YES)
            m->numParsSets += (numLocalChains + 1) * nIntNodes;

        if (m->parsModelId == YES)
            m->numParsNodeLens = (numLocalChains + 1) * nNodes;
        else
            m->numParsNodeLens = 0;
        }
        
    /* then allocate space for the sets and node lengths */
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        mp = &modelParams[d];

        m->parsSets = (BitsLong **) SafeCalloc (m->numParsSets, sizeof(BitsLong*));
        if (!m->parsSets)
            return (ERROR);
        for (i=0; i<m->numParsSets; i++)
            {
            m->parsSets[i] = (BitsLong *) SafeCalloc (m->numChars*m->nParsIntsPerSite, sizeof(BitsLong));
            if (!m->parsSets[i])
                return (ERROR);
            }

        if (m->numParsNodeLens > 0)
            {
            m->parsNodeLens = (CLFlt *) SafeCalloc (m->numParsNodeLens, sizeof(CLFlt));
            if (!m->parsNodeLens)
                return (ERROR);
            }
        }
    
    /* finally fill in tip parsimony sets */
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        mp = &modelParams[d];

        if (mp->dataType == CONTINUOUS)
            {
            /* Note: This is only a placeholder since continuous characters are not implemented yet.
               Using additive parsimony would be more efficient than using multiple binary chars as here. */
            for (i=0; i<numLocalTaxa; i++)
                {
                for (c=0, j=m->compMatrixStart; j<m->compMatrixStop; j++, c++)
                    {
                    x = compMatrix[pos(i,j,compMatrixRowSize)];

                    for (k=0; k<m->nParsIntsPerSite; k++)
                        {
                        if (x > (unsigned int)(k + 1) * 1000 / (m->nParsIntsPerSite + 1))
                            m->parsSets[i][c*m->nParsIntsPerSite + k] = 1;
                        else
                            m->parsSets[i][c*m->nParsIntsPerSite + k] = 2;
                        }
                    }
                }
            }
        else if (m->nCharsPerSite == 1 && m->nParsIntsPerSite == 1)
            {
            allAmbig = (bitsLongOne<<mp->nStates) - 1UL;
            for (i=0; i<numLocalTaxa; i++)
                {
                for (c=0, j=m->compMatrixStart; j<m->compMatrixStop; j++, c++)
                    {
                    x = compMatrix[pos(i,j,compMatrixRowSize)];

                    if (x == MISSING || x == GAP)
                        m->parsSets[i][c] = allAmbig;
                    else
                        m->parsSets[i][c] = x;
                    }
                }
            }
        else if (!strcmp(mp->nucModel, "Doublet") && (mp->dataType == DNA || mp->dataType == RNA))
            {
            allAmbig = 15;
            for (i=0; i<numLocalTaxa; i++)
                {
                for (c=0, j=m->compMatrixStart; j<m->compMatrixStop; j+=m->nCharsPerSite, c++)
                    {
                    /* fetch the original values x1 and x2 */
                    x1 = compMatrix[pos(i,j,compMatrixRowSize)];
                    if (x1 == MISSING || x1 == GAP)
                        x1 = allAmbig;
                    x2 = compMatrix[pos(i,j+1,compMatrixRowSize)];
                    if (x2 == MISSING || x2 == GAP)
                        x2 = allAmbig;
                    /* squeeze them together in the new value x */
                    x = 0;
                    for (nuc1=0; nuc1<4; nuc1++)
                        {
                        for (nuc2=0; nuc2<4; nuc2++)
                            {
                            if (IsBitSet(nuc1,&x1) == YES && IsBitSet(nuc2, &x2) == YES)
                                x |= (bitsLongOne<<(nuc1*4 + nuc2));
                            }
                        }
                    
                    m->parsSets[i][c] = x;
                    }
                }
            }
        else if (!strcmp(mp->nucModel, "Codon") && (mp->dataType == DNA || mp->dataType == RNA))
            {
            allAmbig = 15;
            for (i=0; i<numLocalTaxa; i++)
                {
                for (c=0, j=m->compMatrixStart; j<m->compMatrixStop; j+=m->nCharsPerSite, c++)
                    {
                    /* fetch the original values x1, x2, and x3*/
                    x1 = compMatrix[pos(i,j,compMatrixRowSize)];
                    if (x1 == MISSING || x1 == GAP)
                        x1 = allAmbig;
                    x2 = compMatrix[pos(i,j+1,compMatrixRowSize)];
                    if (x2 == MISSING || x2 == GAP)
                        x2 = allAmbig;
                    x3 = compMatrix[pos(i,j+2,compMatrixRowSize)];
                    if (x3 == MISSING || x3 == GAP)
                        x3 = allAmbig;

                    /* squeeze them together in the new long string pointed to by longPtr */
                    longPtr = &m->parsSets[i][c*m->nParsIntsPerSite];
                    allNucCode = codingNucCode = 0;
                    for (nuc1=0; nuc1<4; nuc1++)
                        for (nuc2=0; nuc2<4; nuc2++)
                            for (nuc3=0; nuc3<4; nuc3++)
                                {
                                if (mp->codon[allNucCode] != 21)
                                    {
                                    if (IsBitSet(nuc1, &x1) == YES && IsBitSet(nuc2, &x2) == YES && IsBitSet(nuc3, &x3) == YES)
                                        SetBit(codingNucCode, longPtr);
                                    codingNucCode++;
                                    }
                                allNucCode++;
                                }
                    }
                }
            }
        else if (!strcmp(mp->nucModel, "Protein") && (mp->dataType == DNA || mp->dataType == RNA))
            {
            allAmbig = 15;
            for (i=0; i<numLocalTaxa; i++)
                {
                for (c=0, j=m->compMatrixStart; j<m->compMatrixStop; j+=m->nCharsPerSite, c++)
                    {
                    /* fetch the original values x1, x2, and x3*/
                    x1 = compMatrix[pos(i,j,compMatrixRowSize)];
                    if (x1 == MISSING || x1 == GAP)
                        x1 = allAmbig;
                    x2 = compMatrix[pos(i,j+1,compMatrixRowSize)];
                    if (x2 == MISSING || x2 == GAP)
                        x2 = allAmbig;
                    x3 = compMatrix[pos(i,j+2,compMatrixRowSize)];
                    if (x3 == MISSING || x3 == GAP)
                        x3 = allAmbig;

                    /* squeeze them together in the new long string pointed to by longPtr */
                    longPtr = &m->parsSets[i][c*m->nParsIntsPerSite];   /* m->nParsIntsPerSite should be 1 */
                    allNucCode = 0;
                    for (nuc1=0; nuc1<4; nuc1++)
                        for (nuc2=0; nuc2<4; nuc2++)
                            for (nuc3=0; nuc3<4; nuc3++)
                                {
                                if (mp->codon[allNucCode] != 21)
                                    {
                                    if (IsBitSet(nuc1, &x1) == YES && IsBitSet(nuc2, &x2) == YES && IsBitSet(nuc3, &x3) == YES)
                                        SetBit(mp->codon[allNucCode]-1, longPtr);
                                    }
                                allNucCode++;
                                }
                    }
                }
            }
        else
            {
            MrBayesPrint ("%s   Unrecognized data format during bitset compression\n", spacer);
            return ERROR;
            }
        }
    
    return (NO_ERROR);
}


/*------------------------------------------------
|
|   InitPrintParams: Set up arrays of print para-
|      meters and print tree parameters
|
------------------------------------------------*/
int InitPrintParams (void)
{
    int     i, j, k, k1=0;
    Param   *p;

    /* count number of model params to print */
    numPrintParams = 0;
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->printParam == YES &&
            p->paramType != P_TOPOLOGY &&
            p->paramType != P_BRLENS &&
            p->paramType != P_SPECIESTREE &&
            p->paramType != P_CPPEVENTS &&
            p->paramType != P_TK02BRANCHRATES &&
            p->paramType != P_IGRBRANCHRATES &&
            p->paramType != P_MIXEDBRCHRATES)
            numPrintParams++;
        }

    /* count number of tree params to print */
    numPrintTreeParams = 0;
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_TOPOLOGY)
            {
            /* always print parsimony topology (printParam == YES), otherwise */
            /* print topology only if brlens never requested (nPrintSubParams == 0)*/
            if (p->printParam == YES || p->nPrintSubParams == 0)
                numPrintTreeParams++;
            }
        else if (p->paramType == P_BRLENS)
            {
            /* print only if brlens (or events) requested for at least one partition */
            if (p->printParam == YES || p->nPrintSubParams > 0)
                numPrintTreeParams++;
            }
        else if (p->paramType == P_SPECIESTREE)
            {
            /* always print if printParam set to YES */
            if (p->printParam == YES)
                numPrintTreeParams++;
            }
        }

    /* allocate space */
    printParam = (Param **) SafeCalloc (numPrintParams + numPrintTreeParams + numTopologies, sizeof(Param *));
    topologyPrintIndex = (int *) SafeCalloc (numTopologies + numPrintTreeParams, sizeof(int)); 
    if (!printParam || !topologyPrintIndex)
        {
        free (printParam);
        free (topologyPrintIndex);
        MrBayesPrint ("%s   Could not allocate printParam vector in InitPrintParams\n", spacer);
        return (ERROR);
        }
    printTreeParam = printParam + numPrintParams;
    topologyParam = printTreeParam + numPrintTreeParams;
    printTreeTopologyIndex = topologyPrintIndex + numTopologies;
    memAllocs[ALLOC_PRINTPARAM] = YES;

    /* assign normal print params */
    for (i=j=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->printParam == YES &&
            p->paramType != P_TOPOLOGY &&
            p->paramType != P_BRLENS &&
            p->paramType != P_SPECIESTREE &&
            p->paramType != P_CPPEVENTS &&
            p->paramType != P_TK02BRANCHRATES &&
            p->paramType != P_IGRBRANCHRATES &&
            p->paramType != P_MIXEDBRCHRATES)
            printParam[j++] = p;
        }
    
    /* assign tree print params */
    for (i=j=k=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_TOPOLOGY)
            {
            /* always print parsimony topology (printParam == YES), otherwise */
            /* print topology only if brlens never requested (nPrintSubParams == 0)*/
            if (p->printParam == YES || p->nPrintSubParams == 0)
                numPrintTreeParams++;
            }
        else if (p->paramType == P_BRLENS)
            {
            /* print only if brlens (or events) requested for at least one partition */
            if (p->printParam == YES || p->nPrintSubParams > 0)
                printTreeParam[k++] = p;
            }
        else if (p->paramType == P_SPECIESTREE)
            {
            if (p->printParam == YES)
                printTreeParam[k++] = p;
            }
        }

    /* find topologies, topology file index, and printtree topology index */
    for (i=0; i<numPrintTreeParams; i++)
        printTreeTopologyIndex[i] = numTopologies;
    for (i=j=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_SPECIESTREE)
            {
            topologyParam[j] = p;
            for (k=0; k<numPrintTreeParams; k++)
                if (printTreeParam[k] == p)
                    break;
            topologyPrintIndex[j] = k;
            printTreeTopologyIndex[k] = j;
            j++;
            }
        else if (p->paramType == P_TOPOLOGY)
            {
            topologyParam[j] = p;
            for (k=0; k<numPrintTreeParams; k++)
                if (printTreeParam[k] == p)
                    break;
            if (k<numPrintTreeParams)
                {
                topologyPrintIndex[j] = k;
                printTreeTopologyIndex[k] = j;
                }
            else
                {
                for (k=0; k<p->nSubParams; k++)
                    {
                    for (k1=0; k1<numPrintTreeParams; k1++)
                        if (printTreeParam[k1] == p->subParams[k])
                            break;
                    if (k1 < numPrintTreeParams)
                        break;
                    }
                topologyPrintIndex[j] = k1;
                printTreeTopologyIndex[k1] = j;
                }
            j++;
            }
        }

    return (NO_ERROR);
}


int IsPFNodeEmpty (PFNODE *p)
{
    int i;

    for (i=0; i<chainParams.numRuns; i++)
        {
        if (p->count[i] > 0)
            break;
        }
    if (i == chainParams.numRuns)
        return YES;
    else
        return NO;
}


/* LargestNonemptyPFNode: recursive function to largest nonempty node in a subtree */
PFNODE *LargestNonemptyPFNode (PFNODE *p, int *i, int j)
{
    PFNODE *q;

    ++j;
    if (p == NULL)
        return NULL;
    
    q = LargestNonemptyPFNode (p->left, i, j);
    
    if (q != NULL)
        {
        return q;
        }
    else if (IsPFNodeEmpty (p) == NO)
        {
        *i = j;
        return p;
        }
    else
        {
        return LargestNonemptyPFNode (p->right, i, j);
        }
}


/* ln prior ratio for clock trees */
int LogClockTreePriorRatio (Param *param, int chain, MrBFlt *lnPriorRatio)
{
    MrBFlt          oldLnPrior, newLnPrior, theta, N, growth, clockRate, sF, *sR, *eR, *fR;
    char            *sS;
    Model           *mp;
    ModelInfo       *m;
    Tree            *newTree, *oldTree;
    TreeNode        *p, *q=NULL;
    int             i, j;

    (*lnPriorRatio) = 0.0;
    
    mp = &modelParams[param->relParts[0]];
    m  = &modelSettings[param->relParts[0]];
    
    newTree = GetTree (m->brlens, chain, state[chain]);
    oldTree = GetTree (m->brlens, chain, state[chain] ^ 1);

    if (m->clockRate != NULL)
        clockRate = *GetParamVals(m->clockRate, chain, state[chain]);
    else
        clockRate = 1.0;
    
    /* calculate prior ratio on brlens of clock tree */
    if (!strcmp(mp->clockPr,"Coalescence"))
        {
        /* coalescence prior */
        /* first calculate theta as 4*N*mu, 3*N*mu or 2*N*mu */
        N = *(GetParamVals (m->popSize, chain, state[chain]));
        if (!strcmp(mp->ploidy, "Diploid"))
            theta = 4 * N * clockRate;
        else if (!strcmp(mp->ploidy, "Zlinked"))
            theta = 3 * N * clockRate;
        else
            theta = 2 * N * clockRate;
        /* deal with growth */
        if (!strcmp(mp->growthPr, "Fixed"))
            growth = mp->growthFix;
        else
            growth = *(GetParamVals (m->growthRate, chain, state[chain]));
        if (LnCoalescencePriorPr (oldTree, &oldLnPrior, theta, growth) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for coalescence process\n", spacer);
            return (ERROR);
            }
        if (LnCoalescencePriorPr (newTree, &newLnPrior, theta, growth) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for coalescence process\n", spacer);
            return (ERROR);
            }
        (*lnPriorRatio) = (newLnPrior - oldLnPrior);
        }
    else if (!strcmp(mp->clockPr,"Birthdeath"))
        {
        /* birth-death prior */
        sR = GetParamVals (m->speciationRates, chain, state[chain]);
        eR = GetParamVals (m->extinctionRates, chain, state[chain]);
        sS = mp->sampleStrat;
        sF = mp->sampleProb;
        if (LnBirthDeathPriorPr (oldTree, clockRate, &oldLnPrior, *sR, *eR, sS, sF) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
            return (ERROR);
            }
        if (LnBirthDeathPriorPr (newTree, clockRate, &newLnPrior, *sR, *eR, sS, sF) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
            return (ERROR);
            }
        (*lnPriorRatio) = (newLnPrior - oldLnPrior);
        }
    else if (!strcmp(mp->clockPr,"Fossilization"))
        {
        /* fossilized birth-death prior */
        sR = GetParamVals (m->speciationRates, chain, state[chain]);
        eR = GetParamVals (m->extinctionRates, chain, state[chain]);
        sF = mp->sampleProb;
        fR = GetParamVals (m->fossilizationRates, chain, state[chain]);
        sS = mp->sampleStrat;
        if (LnFossilizationPriorPr (oldTree, clockRate, &oldLnPrior, sR, eR, sF, fR, sS) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
            return (ERROR);
            }
        if (LnFossilizationPriorPr (newTree, clockRate, &newLnPrior, sR, eR, sF, fR, sS) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
            return (ERROR);
            }
        (*lnPriorRatio) = (newLnPrior - oldLnPrior);
        }
    else if (!strcmp(mp->clockPr,"Uniform"))
        {
        oldLnPrior = LnUniformPriorPr(oldTree, clockRate);
        newLnPrior = LnUniformPriorPr(newTree, clockRate);
        (*lnPriorRatio) = (newLnPrior - oldLnPrior);
        }
    else if (!strcmp(mp->clockPr,"Speciestreecoalescence"))
        {
        // Defer this calculation to the BEST code
        }

    assert (*lnPriorRatio > NEG_INFINITY);

    /* take care of calibrations */
    if (newTree->isCalibrated == YES)
        {
        for (i=0; i<newTree->nNodes-1; i++)
            {
            p = newTree->allDownPass[i];
            if (p->isDated == NO)
                continue;
            for (j=0; j<oldTree->nNodes-1; j++)
                {
                q = oldTree->allDownPass[j];
                if (p->lockID == q->lockID)
                    break;
                }
            assert (j != oldTree->nNodes-1);
            if (p->isDated == YES && p->calibration->prior != fixed)
                {
                (*lnPriorRatio) += p->calibration->LnPriorRatio(p->age, q->age, p->calibration->priorParams);
                }
            }
        }

    return (NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   LogLike: calculate the log likelihood of the new state of the chain
|
-----------------------------------------------------------------*/
MrBFlt LogLike (int chain)
{
    int             i, d;
    ModelInfo       *m;
    MrBFlt          chainLnLike, lnL;
                        
    /* initialize chain cond like */
    chainLnLike = 0.0;
    
    if (chainParams.runWithData == NO)
        return (chainLnLike);

#   if defined (DEBUG_RUN_WITHOUT_DATA)
    return (chainLnLike);
#   endif

#   if defined (THREADS_ENABLED)
    if (tryToUseThreads && ScheduleLogLikeForAllDivisions()) 
        {
        /* Launch all divisions that require updating simultaneously */
        chainLnLike = LaunchLogLikeForAllDivisionsInParallel(chain);
        } 
    else 
        {
        /* Launch divisions in series */
#   endif
        
    /* Cycle through divisions and recalculate tis and cond likes as necessary. */
    /* Code below does not try to avoid recalculating ti probs for divisions    */
    /* that could share ti probs with other divisions.                          */
    for (d=0; d<numCurrentDivisions; d++)
        {
        
#   if defined (BEST_MPI_ENABLED)
        if (isDivisionActive[d] == NO)
            continue;
#   endif
        m = &modelSettings[d];
        if (m->upDateCl == YES) 
            {   
            /* Work has been delegated to a separate function so we can wrap    */
            /* a thread around it                                               */              
            LaunchLogLikeForDivision(chain, d, &(m->lnLike[2 * chain + state[chain]]));
            }
        if (abortMove == YES)
            return MRBFLT_NEG_MAX;
        chainLnLike += m->lnLike[2*chain + state[chain]];   
        }
        
#   if defined (THREADS_ENABLED)
        }
#   endif

    /* unmark all divisions */
    if (chainHasAdgamma == YES)
        {
        for (d=0; d<numCurrentDivisions; d++)
            {
            m = &modelSettings[d];
            m->mark = NO;
            }
        
        /* update HMM likelihoods if appropriate */
        for (d=0; d<numCurrentDivisions; d++)
            {
#   if defined (BEST_MPI_ENABLED)
            if (isDivisionActive[d] == NO)
                continue;
#   endif
            m = &modelSettings[d];
            
            if (m->upDateCl == YES && m->correlation != NULL && m->mark != YES)
                {
                lnL = 0.0;
                CalcLikeAdgamma(d, m->correlation, chain, &lnL);

                /* store the value for the cases where the HMM is not touched */
                m->lnLike[2*chain + state[chain]] =  lnL;
                
                /* add it to chainLnLike - it was not added above since the division */
                /* lnL was set to zero after the update call to Likelihood_Adgamma */
                chainLnLike += lnL;
                
                /* set mark for other divisions in the HMM
                   (i.e., those with the same correlation parameter AND the same shape parameter) */
                for (i=0; i<m->correlation->nRelParts; i++)
                    {
                    if (modelSettings[m->correlation->relParts[i]].shape == modelSettings[d].shape)
                        {
                        modelSettings[m->correlation->relParts[i]].mark = YES;
                        }
                    }
                }
            }
        }

    return (chainLnLike);   
}


MrBFlt LogOmegaPrior (MrBFlt w1, MrBFlt w2, MrBFlt w3)
{

    /* This function returns the log prior probability of 
       the ratio on three omegas. Here, we have three
       nonsynonymous/synonymous rate ratios, denoted w1, w2,
       and w3. They have the property that w1 < w2 < w3. 
       Remember that w1 = dN1/dS, w2 = dN2/dS, and
       w3 = dN3/dS. We assume that dN1, dN2, dN3, and dS
       are all independent draws from the same exponential
       distribution, and that dN1, dN2, and dN3 are the
       order statistics. The w1, w2, and w3, then, are
       all scaled to the same dS r.v. */
       
    MrBFlt  lnProb;
    
    lnProb = log(36.0) - 4.0 * log(1.0 + w1 + w2 + w3);
     
    return (lnProb);
}

 
/* ------------------------------------------------------------------------------------------------------------- */
/* Joint distribution of branch lengths t under gamma-Dirichlet prior:                                           */
/* (Zhang et al. 2012, Eq. 4; Rannala et al. 2012, Eq. 36):                                                      */
/* ln[f(t|aT,bT,a,c)] =  (aT - a*s - a*c*(s-3)) * ln(T) - bT * T + (a-1) * sum[ln(t_j)] + (a*c-1) * sum[ln(t_k)] */
/*                      + aT * ln(bT) - lnG(aT) - lnB(a,c)                                                       */
/*                                                                                                               */
/* Joint distribution of branch lengths t under invgamma-Dirichlet prior:                                        */
/* (Zhang et al. 2012, Eq. 6; Rannala et al. 2012, Eq. 39):                                                      */
/* ln[f(t|aT,bT,a,c)] = (-aT - a*s - a*c*(s-3)) * ln(T) - bT / T + (a-1) * sum[ln(t_j)] + (a*c-1) * sum[ln(t_k)] */
/*                      + aT * ln(bT) - lnG(aT) - lnB(a,c)                                                       */
/* also see DoCitations()                                                                                        */
/* ------------------------------------------------------------------------------------------------------------- */

/* external (tip): 1, internal: 0 */
#define IsTip(Node) (Node->index < numTaxa || (Node->anc)->index < numTaxa)

MrBFlt LogDirPrior (Tree *t, ModelParams *mp, int PV)
{
    /* ln prior prob. under Dirichlet priors and twoExp prior
     //chi */

    int    i, nb[2] = {0,0};
    MrBFlt lnprior = 0.0, tb[2] = {0,0}, treeL = 0.0;
    MrBFlt aT, bT, a, c;
    TreeNode  *p;
    
    /* Not safe, should define Marcos. YES or NO should never be defined to 2 or 3 or 4! */
    /* PV is 2 or 3: Dirichlet priors */    
    if (PV == 2 || PV == 3)
        {
        /* partially for calculating lnPriorRatio, full part is in LogPrior() */
        aT = mp->brlensDir[0];
        bT = mp->brlensDir[1];
        a  = mp->brlensDir[2];
        c  = mp->brlensDir[3];
    
        for (i = 0; i < t->nNodes; i++)
            {
            p = t->allDownPass[i];
            if (p->anc != NULL)
                {
                treeL += p->length;
                nb[IsTip(p)]++;
                tb[IsTip(p)] += log(p->length);
                }
            }
        lnprior += (a-1)*tb[1] + (a*c -1)*tb[0];
        if (PV == 2)
            lnprior += (aT - a*nb[1] - a*c*nb[0]) * log(treeL) - bT*treeL;
        else
            lnprior += (-aT - a*nb[1] - a*c*nb[0]) * log(treeL) - bT/treeL;
        }
    /* or 4: twoExp prior */
    else if (PV == 4)
        {
        for (i = 0; i < t->nNodes; i++) {
            p = t->allDownPass[i];
            if (p->anc != NULL)
                {
                nb[IsTip(p)]++;
                tb[IsTip(p)] += p->length;
                }
            }
        for (i = 0; i < 2; i++)
            lnprior += nb[i] * log(mp->brlens2Exp[i]) - tb[i] * (mp->brlens2Exp[i]);
        }
    
    return lnprior;
}


MrBFlt LogPrior (int chain)
{
    int             i, j, c, n, nStates, *nEvents, sumEvents, *ist, nRates, nParts[6];
    const int       *rateCat;
    MrBFlt          *st, *sst, lnPrior, sum, x, clockRate, theta, popSize, growth, *alphaDir, newProp[190],
                    sF, *sR, *eR, *fR,  freq, pInvar, lambda, sigma, nu, igrvar, **rateMultiplier;
    char            *sS;
    CLFlt           *nSitesOfPat;
    Param           *p;
    ModelParams     *mp;
    ModelInfo       *m;
    Tree            *t;
    TreeNode        *branch, *q;
    
    /* In the stepping-stone method, constants that appear in prior probability density functions must be fully specified.
       It is not necessary for Bayesian MCMC analyses because such constants cancel out. */

    lnPrior = 0.0;
    for (n=0; n<numParams; n++)
        {
        p = &params[n];
#   if defined (MPI_BEST_ENABLED)
        /* We skip all parameters that are not handled on this processor. The scheme used here
           requires that parameters either be unique to one partition (processor) or that they
           are shared across all partitions and that the first processor has all the relevant
           information about that parameter. */
        if (isDivisionActive[p->relParts[0]] == NO)
            continue;
#   endif
        
        st  = GetParamVals (p, chain, state[chain]);
        sst = GetParamSubVals (p, chain, state[chain]);
        mp = &modelParams[p->relParts[0]];
        m = &modelSettings[p->relParts[0]];

        if (p->paramType == P_TRATIO)
            {
            /* tratio parameter */
            if (p->paramId == TRATIO_DIR)
                {
                alphaDir = mp->tRatioDir;
                /* we convert here from the ratio parameterization used in the parameter
                   struct to the simplex parameterization used for the prior */
                newProp[0] =  st[0] / (st[0] + 1.0);
                newProp[1] =  (1.0 - newProp[0]);
                x = LnGamma(alphaDir[0]+alphaDir[1]) - LnGamma(alphaDir[0]) - LnGamma(alphaDir[1]);
                for (i=0; i<2; i++)
                    x += (alphaDir[i]-1.0)*log(newProp[i]);
                lnPrior += x;
                }
            }
        else if (p->paramType == P_REVMAT)
            {
            /* revmat parameter */
            if (p->paramId == REVMAT_DIR)
                {
                if (p->nValues == 6)
                    alphaDir = mp->revMatDir;
                else /* if (p->nValues == 190) */
                    alphaDir = mp->aaRevMatDir;
                sum = 0.0;
                for (i=0; i<p->nValues; i++)
                    sum += alphaDir[i];
                x = LnGamma(sum);
                for (i=0; i<p->nValues; i++)
                    x -= LnGamma(alphaDir[i]);
                for (i=0; i<p->nValues; i++)
                    x += (alphaDir[i] - 1.0) * log(st[i]);
                lnPrior += x;
                }
            else if (p->paramId == REVMAT_MIX)
                {
                assert (p->nValues == 6);
                alphaDir = &mp->revMatSymDir;
                ist      = GetParamIntVals(p, chain, state[chain]); /* the growth fxn */
                nRates   = GetKFromGrowthFxn(ist);
                /* get the actual rate proportions of the current groups */
                for (i=0; i<nRates; i++)
                    {
                    newProp[i] = 0.0;
                    nParts[i] = 0;
                    }
                for (i=0; i<6; i++)
                    {
                    nParts[ist[i]]++;
                    newProp[ist[i]] += st[i];
                    }
                /* now we calculate probability as usual, with alpha
                   parameter multiplied by number of parts */
                x = LnGamma(6.0 * alphaDir[0]);
                for (i=0; i<nRates; i++)
                    x -= LnGamma(nParts[i] * alphaDir[0]);
                for (i=0; i<nRates; i++)
                    x += (nParts[i] * alphaDir[0] - 1.0) * log(newProp[i]);
                /* finally take model probability into account */
                x += log (1.0 / 203);
                lnPrior += x;
                }
            else
                {
                /* fixed or empirical */
                }
            }
        else if (p->paramType == P_OMEGA)
            {
            /* account for prior on omega proportion if 1 omega category */
            if (p->paramId == OMEGA_DIR)
                {
                alphaDir = mp->omegaDir;
                /* convert from ratio parameterization to simplex representation */
                newProp[0] = st[0] / (st[0] + 1.0);
                newProp[1] = 1.0 - newProp[0];
                x = LnGamma(alphaDir[0]+alphaDir[1]) - LnGamma(alphaDir[0]) - LnGamma(alphaDir[1]);
                for (i=0; i<2; i++)
                    x += (alphaDir[i]-1.0)*log(newProp[i]);
                lnPrior += x;
                }
            /* account for stationary state frequencies of M3 and ny98 omega categories */
            if (p->paramId == OMEGA_BUD || p->paramId == OMEGA_BED ||
                p->paramId == OMEGA_BFD || p->paramId == OMEGA_FUD ||
                p->paramId == OMEGA_FED || p->paramId == OMEGA_FFD ||
                p->paramId == OMEGA_ED  || p->paramId == OMEGA_FD)
                {
                alphaDir = mp->codonCatDir;             
                x = 0.0;
                for (i=0; i<3; i++)
                    x += (alphaDir[i]-1.0)*log(sst[i]);
                lnPrior += x;
                }
            /* account for beta prior on omeganeg in NY98 */
            if (p->paramId == OMEGA_BUD || p->paramId == OMEGA_BUF ||
                p->paramId == OMEGA_BED || p->paramId == OMEGA_BEF ||
                p->paramId == OMEGA_BFD || p->paramId == OMEGA_BFF)
                {
                alphaDir = mp->ny98omega1Beta;
                newProp[0] = st[0] / (st[0] + 1.0);
                newProp[1] = 1.0 - newProp[0];
                x = 0.0;
                for (i=0; i<2; i++)
                    x += (alphaDir[i]-1.0)*log(newProp[i]);
                lnPrior += x;
                }
            /* account for omegapos in NY98 with uniform prior prob */
            if (p->paramId == OMEGA_BUD || p->paramId == OMEGA_BUF ||
                p->paramId == OMEGA_FUD || p->paramId == OMEGA_FUF)
                {
                lnPrior += log(1.0) - log(mp->ny98omega3Uni[1] - mp->ny98omega3Uni[0]);
                }
            /* account for omegapos in NY98 with exponential prior prob */
            if (p->paramId == OMEGA_BED || p->paramId == OMEGA_BEF ||
                p->paramId == OMEGA_FED || p->paramId == OMEGA_FEF)
                {
                lnPrior += (log(mp->ny98omega3Exp) - mp->ny98omega3Exp * st[2]);
                }
            /* account for omegas in M3, which can only be exponential; if fixed, ln prior prob is 0 */
            if (p->paramId == OMEGA_EF || p->paramId == OMEGA_ED)
                {
                lnPrior += LogOmegaPrior (st[0], st[1], st[2]);
                }
            }
        else if (p->paramType == P_PI)
            {
            /* state frequencies parameter */
            if (p->paramId == PI_DIR)
                {
                nStates = p->nSubValues;
                sum = 0.0;
                for (i=0; i<nStates; i++)
                    sum += st[i];
                x = LnGamma(sum);
                for (i=0; i<nStates; i++)
                    x -= LnGamma(st[i]);
                for (i=0; i<nStates; i++)
                    x += (st[i] - 1.0)*log(sst[i]);
                lnPrior += x;
                }
            else if (p->paramId == SYMPI_EXP || p->paramId == SYMPI_EXP_MS)
                {
                lnPrior += - mp->symBetaExp * st[0] + log(mp->symBetaExp);
                }
            else if (p->paramId == SYMPI_UNI || p->paramId == SYMPI_UNI_MS)
                {
                lnPrior += log(1.0) - log(mp->symBetaUni[1] - mp->symBetaUni[0]);
                }                
            if (p->paramId == SYMPI_EXP_MS || p->paramId == SYMPI_UNI_MS || p->paramId == SYMPI_FIX_MS)
                {
                sst = GetParamStdStateFreqs(p, chain, state[chain]);
                sst += 2*m->numBetaCats;
                for (i=0; i<p->nSympi; i++)
                    {
                    nStates = p->sympinStates[i];
                    x = LnGamma(nStates*st[0]) - nStates*LnGamma(st[0]);
                    for (j=0; j<nStates; j++)
                        x += (st[0] - 1.0)*log(sst[j]);
                    lnPrior += x;
                    sst += nStates;
                    }
                }
            }
        else if (p->paramType == P_MIXTURE_RATES)
            {
            /* site rate mixture; first deal with dirichlet prior */
            nStates = p->nSubValues;
            sum = 0.0;
            for (i=0; i<nStates; i++)
                sum += st[i];
            x = LnGamma(sum);
            for (i=0; i<nStates; i++)
                x -= LnGamma(st[i]);
            for (i=0; i<nStates; i++)
                x += (st[i] - 1.0)*log(sst[i]/nStates);
                
            /* now deal with the order statistic; we multiply the density by the number of orderings */
            x += LnFactorial(nStates);

            lnPrior += x;
        }
        else if (p->paramType == P_SHAPE)
            {
            /* gamma/lnorm shape parameter */
            if (p->paramId == SHAPE_UNI)
                {
                lnPrior += log(1.0) - log(mp->shapeUni[1] - mp->shapeUni[0]);
                }
            else if (p->paramId == SHAPE_EXP)
                {
                lnPrior += log(mp->shapeExp) - mp->shapeExp * st[0];
                }
            for (i=0; i<p->nRelParts; i++)
                {
                m = &modelSettings[p->relParts[i]];
                if (m->gibbsGamma == YES)
                    {
                    if (m->pInvar == NULL)
                        lnPrior += log(1.0/m->numRateCats) * m->numUncompressedChars;
                    else
                        {
                        rateCat = m->tiIndex + chain * m->numChars;
                        pInvar = *GetParamVals (m->pInvar, chain, state[chain]);
                        nSitesOfPat = numSitesOfPat + m->compCharStart;
                        freq = (1.0 - pInvar)/m->numRateCats;
                        for (c=0; c<m->numChars; c++)
                            {
                            if (rateCat[c] < m->numRateCats)
                                lnPrior += log(freq) * nSitesOfPat[c];
                            else
                                lnPrior += log(pInvar) * nSitesOfPat[c];
                            }
                        }
                    }
                }
            }
        else if (p->paramType == P_PINVAR)
            {
            /* proportion of invariable sites parameter */
            lnPrior += log(1.0) - log(mp->pInvarUni[1] - mp->pInvarUni[0]);
            }
        else if (p->paramType == P_CORREL)
            {
            /* adGamma model parameter */
            lnPrior += log(1.0) - log(mp->corrUni[1] - mp->corrUni[0]);
            }
        else if (p->paramType == P_SWITCH)
            {
            /* switching rate parameter of covarion model */
            if (p->paramId == SWITCH_UNI)
                {
                lnPrior += log(1.0) - log(mp->covswitchUni[1] - mp->covswitchUni[0]);
                }
            else if (p->paramId == SWITCH_EXP)
                {
                lnPrior += log(mp->covswitchExp) - mp->covswitchExp * st[0];
                }
            }
        else if (p->paramType == P_RATEMULT && p->nValues > 1)
            {
            nStates = p->nValues;
            sum = 0.0;
            for (i=0; i<nStates; i++)
                sum += sst[i+nStates];
            x = LnGamma(sum);
            for (i=0; i<nStates; i++)
                x -= LnGamma(sst[i+nStates]);
            for (i=0; i<nStates; i++)
                x += (sst[i+nStates] - 1.0) * log(st[i]);
            sum = 0.0;  // the constant
            for (i=0; i<nStates; i++)
                sum += sst[i];
            for (i=0; i<nStates; i++) {
                if (i < nStates-1)
                    x += (sst[i+nStates]) * log(sst[i]/sum);
                else
                    x += (sst[i+nStates] - 1.0) * log(sst[i]/sum);
                }
            lnPrior += x;
            }
        else if (p->paramType == P_GENETREERATE && p->nValues > 1)
            {
            nStates = p->nValues;
            sum = 0.0;
            for (i=0; i<nStates; i++)
                sum += sst[i+nStates];
            x = LnGamma(sum);
            for (i=0; i<nStates; i++)
                x -= LnGamma(sst[i+nStates]);
            for (i=0; i<nStates; i++)
                x += (sst[i+nStates] - 1.0) * log(st[i]);
            sum = 0.0;  // the constant
            for (i=0; i<nStates; i++)
                sum += sst[i];
            for (i=0; i<nStates; i++) {
                if (i < nStates-1)
                    x += (sst[i+nStates]) * log(sst[i]/sum);
                else
                    x += (sst[i+nStates] - 1.0) * log(sst[i]/sum);
                }
            lnPrior += x;
            }
        else if (p->paramType == P_TOPOLOGY)
            {
            // Note that a topology can have several unlinked branch length subparameters but only
            // one set of clock branch lengths. To find all the branch length subparameters of a
            // topology, cycle through the p->subParams, which will contain at least one branch length
            // parameter.
            t = GetTree (p, chain, state[chain]);
            if (t->isClock == YES)
                continue;   /* prior probability taken care of in the brlens parameter */
            if (t->nLocks > 0)
                {
                for (i=0; i<t->nNodes-1; i++)
                    {
                    branch = t->allDownPass[i];
                    if (branch->left == NULL)
                        branch->x = 1;
                    else
                        branch->x = branch->left->x + branch->right->x;
                    if (branch->isLocked == YES || branch->anc->anc == NULL)
                        {
                        for (j = 2*branch->x - 3; j>=1; j-=2)
                            {
                            lnPrior -= log ((MrBFlt)j);
                            }
                        branch->x = 1;
                        }
                    }
                }
            else
                {
                for (j = 2*(t->nNodes-t->nIntNodes)-5; j>=1; j-=2)
                    {
                    lnPrior -= log ((MrBFlt)j);
                    }
                }
            }
        else if (p->paramType == P_BRLENS)
            {
            /* branch lengths */
            t = GetTree (p, chain, state[chain]);
            if (t->isClock == YES)
                {
                if (p->paramId == BRLENS_CLOCK_UNI)
                    {
                    /* uniformly distributed branch lengths */
                    clockRate = *(GetParamVals (m->clockRate, chain, state[chain]));
                    lnPrior += LnUniformPriorPr(t, clockRate);
                    }
                else if (p->paramId == BRLENS_CLOCK_COAL)
                    {
                    /* coalescence prior */
                    popSize   = *(GetParamVals (m->popSize, chain, state[chain]));
                    clockRate = *(GetParamVals (m->clockRate, chain, state[chain]));
                    if (strcmp(mp->ploidy, "Diploid") == 0)
                        theta = 4.0 * popSize * clockRate;
                    else if (strcmp(mp->ploidy, "Zlinked") == 0)
                        theta = 3.0 * popSize * clockRate;
                    else
                        theta = 2.0 * popSize * clockRate;
                    if (!strcmp(mp->growthPr, "Fixed"))
                        growth = mp->growthFix;
                    else
                        growth = *(GetParamVals (m->growthRate, chain, state[chain]));
                    if (LnCoalescencePriorPr (t, &x, theta, growth) == ERROR)
                        {
                        MrBayesPrint ("%s   Problem calculating prior for coalescence process\n", spacer);
                        }
                    lnPrior += x;
                    }
                else if (p->paramId == BRLENS_CLOCK_BD)
                    {
                    /* birth-death prior */
                    sR = GetParamVals (m->speciationRates, chain, state[chain]);
                    eR = GetParamVals (m->extinctionRates, chain, state[chain]);
                    sS = mp->sampleStrat;
                    sF = mp->sampleProb;
                    if (m->clockRate != NULL)
                        clockRate = *(GetParamVals (m->clockRate, chain, state[chain]));
                    else
                        clockRate = 1.0;
                    if (LnBirthDeathPriorPr (t, clockRate, &x, *sR, *eR, sS, sF) == ERROR)
                        {
                        MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
                        }
                    lnPrior += x;
                    }
                else if (p->paramId == BRLENS_CLOCK_FOSSIL)
                    {
                    /* fossilized birth-death prior */
                    sR = GetParamVals (m->speciationRates, chain, state[chain]);
                    eR = GetParamVals (m->extinctionRates, chain, state[chain]);
                    sF = mp->sampleProb;
                    fR = GetParamVals (m->fossilizationRates, chain, state[chain]);
                    sS = mp->sampleStrat;
                    if (m->clockRate != NULL)
                        clockRate = *(GetParamVals (m->clockRate, chain, state[chain]));
                    else
                        clockRate = 1.0;
                    if (LnFossilizationPriorPr (t, clockRate, &x, sR, eR, sF, fR, sS) == ERROR)
                        {
                        MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
                        }
                    lnPrior += x;
                    }
                else if (p->paramId == BRLENS_CLOCK_SPCOAL)
                    {
                    /* delegate this calculation to the P_SPECIESTREE parameter */
                    }
                if (t->isCalibrated == YES)
                    {
                    /* take care of calibrations */
                    for (i=0; i<t->nNodes-1; i++)
                        {
                        q = t->allDownPass[i];
                        if (q->isDated == YES && q->calibration->prior != fixed)
                            {
                            lnPrior += q->calibration->LnPriorProb(q->age, q->calibration->priorParams);
                            }
                        }
                    }
                }
            else
                {
                if (p->paramId == BRLENS_UNI)
                    {
                    for (i=0; i<t->nNodes; i++)
                        {
                        branch = t->allDownPass[i];
                        if (branch->anc != NULL)
                            lnPrior += log(1.0) - log(mp->brlensUni[1] - BRLENS_MIN);
                        }
                    }
                else if (p->paramId == BRLENS_EXP)
                    {
                    for (i=0; i<t->nNodes; i++)
                        {
                        branch = t->allDownPass[i];
                        if (branch->anc != NULL)
                            lnPrior += log(mp->brlensExp) - mp->brlensExp * branch->length;
                        }
                    }
                /* Dirichlet priors */
                else if (p->paramId == BRLENS_GamDir)  
                    {
                    lnPrior += LogDirPrior(t, mp, 2);
                    lnPrior += (mp->brlensDir[0]) * log(mp->brlensDir[1]) - LnGamma(mp->brlensDir[0])
                                + LnGamma (mp->brlensDir[2] * numTaxa + mp->brlensDir[2] * mp->brlensDir[3] * (numTaxa-3))
                                - numTaxa * LnGamma(mp->brlensDir[2]) - (numTaxa-3) * LnGamma(mp->brlensDir[2] * mp->brlensDir[3]);
                    } 
                else if (p->paramId == BRLENS_iGmDir)
                    {
                    lnPrior += LogDirPrior(t, mp, 3);
                    lnPrior += (mp->brlensDir[0]) * log(mp->brlensDir[1]) - LnGamma(mp->brlensDir[0])
                                + LnGamma (mp->brlensDir[2] * numTaxa + mp->brlensDir[2] * mp->brlensDir[3] * (numTaxa-3))
                                - numTaxa * LnGamma(mp->brlensDir[2]) - (numTaxa-3) * LnGamma(mp->brlensDir[2] * mp->brlensDir[3]);
                    }
                /* twoExp prior */
                else if (p->paramId == BRLENS_twoExp)
                    lnPrior += LogDirPrior(t, mp, 4);
                }
            }
        else if (p->paramType == P_SPECRATE)
            {
            /* speciation rate parameter */
            if (p->paramId == SPECRATE_UNI)
                {
                for (i=0; i<p->nValues; i++)
                    lnPrior += log(1.0) - log(mp->speciationUni[1] - mp->speciationUni[0]);
                }
            else if (p->paramId == SPECRATE_EXP)
                {
                for (i=0; i<p->nValues; i++)
                    lnPrior += log(mp->speciationExp) - mp->speciationExp * st[i];
                }
            }
        else if (p->paramType == P_EXTRATE)
            {
            /* extinction rate parameter */
            if (p->paramId == EXTRATE_BETA)
                {
                for (i=0; i<p->nValues; i++)
                    {
                    alphaDir = mp->extinctionBeta;
                    newProp[0] =  st[i];
                    newProp[1] =  (1.0 - newProp[0]);
                    lnPrior += LnGamma(alphaDir[0]+alphaDir[1]) - LnGamma(alphaDir[0]) - LnGamma(alphaDir[1]);
                    lnPrior += (alphaDir[0]-1.0)*log(newProp[0]) + (alphaDir[1]-1.0)*log(newProp[1]);
                    }
                }
            }
        else if (p->paramType == P_FOSLRATE)
            {
            /* fossilization rate parameter */
            if (p->paramId == FOSLRATE_BETA)
                {
                for (i=0; i<p->nValues; i++)
                    {
                    alphaDir = mp->fossilizationBeta;
                    newProp[0] =  st[i];
                    newProp[1] =  (1.0 - newProp[0]);
                    // if (newProp[0] > 0.0) /* to avoid psi=0 in [0, x_cut] under diversified sampling */
                    lnPrior += LnGamma(alphaDir[0]+alphaDir[1]) - LnGamma(alphaDir[0]) - LnGamma(alphaDir[1]);
                    lnPrior += (alphaDir[0]-1.0)*log(newProp[0]) + (alphaDir[1]-1.0)*log(newProp[1]);
                    }
                }
            }
        else if (p->paramType == P_POPSIZE)
            {
            /* neutral coalescence population size parameter; one value, or one value per branch of species tree */
            for (i=0; i<p->nValues; i++)
                {
                lnPrior += p->LnPriorProb(st[i], p->priorParams);
                }
            }
        else if (p->paramType == P_AAMODEL)
            {
            lnPrior += sst[(int)st[0]];
            }
        else if (p->paramType == P_BRCORR)
            {

            }
        else if (p->paramType == P_BRSIGMA)
            {

            }
        else if (p->paramType == P_GROWTH)
            {
            /* population growth parameter */
            if (p->paramId == GROWTH_UNI)
                {
                lnPrior += log(1.0) - log(mp->growthUni[1] - mp->growthUni[0]);
                }
            else if (p->paramId == GROWTH_EXP)
                {
                lnPrior += log(mp->growthExp) - mp->growthExp * st[0];
                }
            }
        else if (p->paramType == P_CPPRATE)
            {
            /* rate (lambda) of comp poisson process of relaxed clock */
            if (p->paramId == CPPRATE_EXP)
                {
                lnPrior += log (mp->cppRateExp) - mp->cppRateExp * st[0];
                }
            }
        else if (p->paramType == P_CPPMULTDEV)
            {
            /* standard deviation (log) of lognormal distribution of rate multipliers for cpp relaxed clock */
            /* only fixed value allowed currently */
            }
        else if (p->paramType == P_CPPEVENTS)
            {
            /* events of CPP relaxed clock process */
            lambda = *GetParamVals (m->cppRate, chain, state[chain]);
            sigma = *GetParamVals (m->cppMultDev, chain, state[chain]);
            nEvents = p->nEvents[2*chain+state[chain]];
            rateMultiplier = p->rateMult[2*chain+state[chain]];
            /* cpp events */
            sumEvents = 0;
            for (i=0; i<2*numLocalTaxa-2; i++)
                sumEvents += nEvents[i];
            t = GetTree (p, chain, state[chain]);
            lnPrior += - lambda * TreeLength (p, chain) + (sumEvents * log (lambda));
            /* rate multipliers */
            for (i=0; i<2*numLocalTaxa-2; i++)
                {
                for (j=0; j<nEvents[i]; j++)
                    lnPrior += LnProbLogNormal (0.0, sigma, rateMultiplier[i][j]);
                }
            for (i=0; i<t->nNodes-2; i++)
                {
                branch = t->allDownPass[i];
                assert (fabs(branch->length - (branch->anc->nodeDepth - branch->nodeDepth)) < 0.000001);
                }
            }
        else if (p->paramType == P_TK02VAR)
            {
            /* variance of rates (nu) in Thorne-Kishino model */
            if (p->paramId == TK02VAR_EXP)
                {
                lnPrior += log (mp->tk02varExp) - mp->tk02varExp * st[0];
                }
            else if (p->paramId == TK02VAR_UNI)
                {
                lnPrior += log(1.0) - log (mp->tk02varUni[1] - mp->tk02varUni[0]);
                }
            }
        else if (p->paramType == P_TK02BRANCHRATES || (p->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(p, chain, state[chain]) == RCL_TK02))
            {
            /* branch rates of Thorne-Kishino model */
            t = GetTree (p, chain, state[chain]);
            if (p->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (m->tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (m->mixedvar, chain, state[chain]);
            for (i=0; i<t->nNodes-2; i++)
                {
                branch = t->allDownPass[i];
                if (branch->length > 0.0)  // not ancestral fossil
                    lnPrior += LnProbTK02LogNormal (st[branch->anc->index], nu*branch->length, st[branch->index]);
                }
            }
        else if (p->paramType == P_IGRVAR)
            {
            /* variance of rates in independent gamma rates model */
            if (p->paramId == IGRVAR_EXP)
                {
                lnPrior += log (mp->igrvarExp) - mp->igrvarExp * st[0];
                }
            else if (p->paramId == IGRVAR_UNI)
                {
                lnPrior += log(1.0) - log (mp->igrvarUni[1] - mp->igrvarUni[0]);
                }
            }
        else if (p->paramType == P_IGRBRANCHRATES || (p->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(p, chain, state[chain]) == RCL_IGR))
            {
            /* branch rates of independent branch rate model */
            t = GetTree (p, chain, state[chain]);
            if (p->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (m->igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (m->mixedvar, chain, state[chain]);
            for (i=0; i<t->nNodes-2; i++)
                {
                branch = t->allDownPass[i];
                if (branch->length > 0.0)  // not ancestral fossil
                    lnPrior += LnProbGamma (branch->length/igrvar, branch->length/igrvar, st[branch->index]);
                assert (fabs(sst[branch->index] - branch->length * st[branch->index]) < BRLENS_MIN);
                assert (fabs(branch->length - (branch->anc->nodeDepth - branch->nodeDepth)) < BRLENS_MIN);
                }
            }
        else if (p->paramType == P_MIXEDVAR)
            {
            /* hyper prior of rates in mixed rel clock model */
            if (p->paramId == MIXEDVAR_EXP)
                {
                lnPrior += log (mp->mixedvarExp) - mp->mixedvarExp * st[0];
                }
            else if (p->paramId == MIXEDVAR_UNI)
                {
                lnPrior += log(1.0) - log (mp->mixedvarUni[1] - mp->mixedvarUni[0]);
                }
            }
        else if (p->paramType == P_CLOCKRATE)
            {
            /* base rate of molecular clock */
            lnPrior += p->LnPriorProb(st[0], p->priorParams);
            }
        else if (p->paramType == P_SPECIESTREE)
            {
            /* calculate prior */
            lnPrior += LnSpeciesTreeProb(chain);
            }
        }
    assert (lnPrior == lnPrior);

#   if defined (BEST_MPI_ENABLED)
    /* Assemble prior probabilities across processors */
    myLnPrior = lnPrior;
    MPI_AllReduce (&myLnPrior, &lnPrior, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#   endif

    return (lnPrior);
}


int LnBirthDeathPriorPr (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt sR, MrBFlt eR, char *sS, MrBFlt sF)
{
    if (!strcmp(sS, "Random")) 
        {
        return LnBirthDeathPriorPrRandom (t, clockRate, prob, sR, eR, sF);
        }
    else if (!strcmp(sS, "Diversity")) 
        {
        return LnBirthDeathPriorPrDiversity (t, clockRate, prob, sR, eR, sF);
        }
    else if (!strcmp(sS, "Cluster")) 
        {
        return LnBirthDeathPriorPrCluster (t, clockRate, prob, sR, eR, sF);
        }
    else 
        {
        MrBayesPrint ("%s   ERROR: Sampling strategy for birth-death process not implemented\n", spacer);
        return (ERROR);
        }
}


/*---------------------------------------------------------------------------------
|
|   LnBirthDeathPriorPrRandom
|
|   We assume a rooted tree that satisfies the molecular clock constraint. The
|   tree is labelled as follows:
|
|                                      t_4 (age of tips)
|     \         \         \        /            
|      \         \         \      /              
|       \         \         \    /         
|        \         \         \  /          
|         \         \         \/       t_3 
|          \         \        /            
|           \         \      /             
|            \         \    /              
|             \         \  /               
|              \         \/            t_2  
|               \        /                 
|                \      /                  
|                 \    /                        
|                  \  /                         
|                   \/                 t_1 (age of most recent common ancestor)
|    
|
|   This function calculates the probability of such a tree under the neutral
|   birth death prior with constant birth and death rates, conditioned on
|   a particular time of the first split, t_1, and a particular number of
|   species, n. We assume rho-sampling, that is, a constant sampling pro-
|   bability rho, which is known, across tips of the tree. Variables:
|
|   T:   the unlabeled oriented tree, which is equivalent to a set of unordered
|        speciation times from a point process
|   tau: the labeled unoriented tree
|   b:   birth (speciation) rate
|   d:   death (extintion) rate
|   f:   sampling fraction
|   n:   number of species in the sampled tree
|
|   See:
|   Tanja Stadler (2009) On incomplete sampling under birth-death models and
|   connections to the sampling-based coalescent. Journal of Theoretical Biology
|   261: 58-66.
|
|   We use f(T|n), which is derived from f(T|n,t_or) using Stadler's approach,
|   in which the time of origin of the tree is associated with a uniform prior
|   and integrated out of the density. We then have:
|
|   We have the following distribution for ordered bifurcation times (cf.
|   equation 5 in Stadler, 2009, simplified here using the p0 and p1 functions):
|
|
|              n! * p1(t_1)    n-1
|   f(T|n) = --------------- * prod (b * p1(t_i))
|             (1 - p0(t_1))    i=1
|
|
|   where   t_1   = time of most recent common ancestor
|           p0(t) = prob. of no descendants surviving and being sampled after time t (see LnP0 function below)
|           p1(t) = prob. of one descendant surviving and being sampled after time t (see LnP1 function below)
|
|   To get the distribution on oriented trees, this density needs to be divided by the
|   number of ways of ordering n-1 bifurcation times, (n-1)!, since each way of ordering
|   the bifurcation times corresponds to a distinct oriented tree. The result is the following
|   density on oritented trees:
|
|              n * p1(t_1)     n-1
|   f(T|n) = --------------- * prod (b * p1(t_i))
|             (1 - p0(t_1))    i=1
|
|
|   To translate this to a density on distinct labeled trees, the density needs to be multiplied by
|   (2^(n-1) / n!).
|
|   For the critical process where the speciation and extinction rates are equal, we obtain the
|   following result in the limit (cf. equation 6 in Stadler (2009)):
|
|                  n          n-1          f*b      
|   f(T|n) = -------------- * prod -----------------
|            (1 + f*b*t_1)    i=1   (1 + f*b*t_i)^2
|
---------------------------------------------------------------------------------*/
int LnBirthDeathPriorPrRandom (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt sR, MrBFlt eR, MrBFlt sF)
{
    int             i, nTaxa;
    MrBFlt          *nt, lambda, mu, rho;
    TreeNode        *p;

    /* allocate space for the speciation times */
    nt = (MrBFlt *)SafeMalloc((size_t)(t->nIntNodes) * sizeof(MrBFlt));
    if (!nt)
        {
        MrBayesPrint ("\n   ERROR: Problem allocating nt\n");
        return (ERROR);
        }

    /* transform to standard variables */
    rho    = sF;
    lambda = sR / (1.0 - eR);
    mu     = eR * lambda;

    /* get the node times and put them into a vector */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        nt[i] = p->nodeDepth / clockRate;
        }
    nTaxa = t->nIntNodes + 1;

    /* calculate probability of tree using standard variables */
    if (AreDoublesEqual(lambda,mu,ETA)==NO)
        {
        // birth rate != death rate, see equation (5) in Stadler (2009) and above
        (*prob) = log(nTaxa) + log(lambda - mu) - (lambda - mu) * nt[t->nIntNodes-1];
        (*prob) -= log(rho*lambda + (lambda*(1.0 - rho) - mu)*exp((mu - lambda)*nt[t->nIntNodes-1]));
        for (i=0; i<t->nIntNodes; i++)
            (*prob) += log(rho*lambda) + LnP1Subsample(nt[i], lambda, mu, rho);
        (*prob) += (nTaxa - 1.0) * log(2.0) - LnFactorial(nTaxa);    /* conversion to labeled tree from oriented tree */
        }
    else
        {
        // birth rate == death rate -> the critical branching process
        (*prob) = log(nTaxa/(1.0 + rho*lambda*nt[t->nIntNodes-1]));
        for (i=0; i<t->nIntNodes; i++)
            (*prob) += log(rho*lambda) - 2.0 * log(1.0 + rho*lambda*nt[i]);
        (*prob) += (nTaxa - 1.0) * log(2.0) - LnFactorial(nTaxa);    /* conversion to labeled tree from oriented tree */
        }

    /* free memory */
    free (nt);
    
    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
 |
 |   LnBirthDeathPriorPrDiversity
 |
 |   Eq.5 in Hohna et al. 2011 MBE
 |
 ---------------------------------------------------------------------------------*/
int LnBirthDeathPriorPrDiversity (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt sR, MrBFlt eR, MrBFlt sF)
{
    int             i, nTaxa, n, m;
    MrBFlt          *nt, lambda, mu, nt_min;
    TreeNode        *p;
    
    /* allocate space for the speciation times */
    nt = (MrBFlt *)SafeMalloc((size_t)(t->nIntNodes) * sizeof(MrBFlt));
    if (!nt)
        {
        MrBayesPrint ("\n   ERROR: Problem allocating nt\n");
        return (ERROR);
        }
    
    /* transform to standard variables */
    lambda = sR / (1.0 - eR);
    mu     = eR * lambda;
    
    n      = t->nIntNodes+1;
    m      = (int)floor(n/sF+0.5); /* equal to round(n/sF) plus it is compatible with MS Visual Studio */
    
    /* get the node times and put them into a vector */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        nt[i] = p->nodeDepth / clockRate;
        }
    nTaxa = t->nIntNodes + 1;

    /* find the youngest interal node */
    nt_min = nt[0];
    for (i=0; i<t->nIntNodes-1; i++)
        {
        if (nt_min > nt[i])
            nt_min = nt[i];
        }
    
    /* calculate probability of tree using standard variables */
    if (AreDoublesEqual(lambda,mu,ETA)==NO)
        {
        // birth rate != death rate
        MrBFlt p0_t1;
        p0_t1 = LnP0(nt[t->nIntNodes-1], lambda, mu);
        (*prob) = log(nTaxa); // we need to add here the binomial coefficient
        (*prob) += (m-n) * (LnP0(nt_min, lambda, mu) - p0_t1);
        for (i=0; i<t->nIntNodes-1; i++)
            (*prob) += (LnP1(nt[i], lambda, mu) - p0_t1);
        (*prob) += (nTaxa - 1.0) * log(2.0) - LnFactorial(nTaxa);  /* conversion to labeled tree from oriented tree */
        }
    else
        {
        MrBayesPrint ("\n   ERROR: Critical branchin process for diversity sampling not implemented\n");
        return (ERROR);
        }

    /* condition on tmrca ??? */

    /* free memory */
    free (nt);
    
    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
 |
 |   LnBirthDeathPriorPrCluster
 |
 |   Eq.7 in Hohna et al. 2011 MBE
 |
 ---------------------------------------------------------------------------------*/
int LnBirthDeathPriorPrCluster (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt sR, MrBFlt eR, MrBFlt sF)
{
    int             i, nTaxa, n, m;
    MrBFlt          *nt, lambda, mu, nt_max;
    TreeNode        *p;
    
    /* allocate space for the speciation times */
    nt = (MrBFlt *)SafeMalloc((size_t)(t->nIntNodes) * sizeof(MrBFlt));
    if (!nt)
        {
        MrBayesPrint ("\n   ERROR: Problem allocating nt\n");
        return (ERROR);
        }
    
    /* transform to standard variables */
    lambda = sR / (1.0 - eR);
    mu     = eR * lambda;
    
    n      = t->nIntNodes+1;
    m      = (int)floor(n/sF+0.5); /* equal to round(n/sF) plus it is compatible with MS Visual Studio */
    
    /* get the node times and put them into a vector */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        nt[i] = p->nodeDepth / clockRate;
        }
    nTaxa = t->nIntNodes + 1;
    
    /* find the second oldest interal node */
    nt_max = nt[0];
    for (i=0; i<t->nIntNodes-1; i++)
        {
        if (nt_max < nt[i])
            nt_max = nt[i];
        }

    /* calculate probability of tree using standard variables */
    if (AreDoublesEqual(lambda,mu,ETA)==NO)
        {
        // birth rate != death rate
        MrBFlt p0_t1;
        p0_t1 = LnP0(nt[t->nIntNodes-1], lambda, mu);
        (*prob) = log(nTaxa); // we need to add here the binomial coefficient
        (*prob) += (m-n) * (LnP0(nt_max, lambda, mu) - p0_t1);
        for (i=0; i<t->nIntNodes-1; i++)
            (*prob) += (LnP1(nt[i], lambda, mu) - p0_t1);
        (*prob) += (nTaxa - 1.0) * log(2.0) - LnFactorial(nTaxa);  /* conversion to labeled tree from oriented tree */
        }
    else
        {
        MrBayesPrint ("\n   ERROR: Critical branchin process for cluster sampling not implemented\n");
        return (ERROR);
        }

    /* condition on tmrca ??? */

    /* free memory */
    free (nt);
    
    return (NO_ERROR);
}


/*
 |
 | The probability of having zero lineages remaining after time t in the
 | birth-death process.
 |
 | param: t - speciation time
 | param: b - birth rate
 | param: d - death rate
 | return: log probability of zero remaining lineages
 |
 */
MrBFlt LnP0 (MrBFlt t, MrBFlt b, MrBFlt d)
{
    MrBFlt      p0t;
    
    p0t = d*(1.0-exp((d-b)*t)) / (b -d*exp((d-b)*t));
    
    return (log(p0t));
}

/*
 |
 | The probability of having zero lineages remaining after time t in the
 | birth-death process.
 |
 | param: t - speciation time
 | param: b - birth rate
 | param: d - death rate
 | param: f - sample frequency
 | return: log probability of zero remaining lineages
 |
 */
MrBFlt LnP0Subsample (MrBFlt t, MrBFlt b, MrBFlt d, MrBFlt f)
{
    MrBFlt      p0t;
    
    p0t = (f*d + (b*(1.0-f) - d)*exp((d-b)*t)) / (f*b + (b*(1.0-f)-d)*exp((d-b)*t));
    
    return (log(p0t));
}

/*
 |
 | The probability of having one lineage remaining after time t
 | in the birth-death process.
 |
 | param: t - speciation time
 | param: b - birth rate
 | param: d - death rate
 | return: log probability of one remaining lineage
 |
 */
MrBFlt LnP1 (MrBFlt t, MrBFlt b, MrBFlt d)
{
    MrBFlt      p0t;
    
    p0t = 2.0 * log(b-d) - (b-d)*t;
    
    p0t -= 2.0 * log(b - d*exp((d-b)*t));
    
    return p0t;
}

/*
 |
 | The probability of having one lineage remaining after time t
 | in the birth-death process.
 |
 | param: t - speciation time
 | param: b - birth rate
 | param: d - death rate
 | param: f - sample frequency
 | return: log probability of one remaining lineage
 |
 */
MrBFlt LnP1Subsample (MrBFlt t, MrBFlt b, MrBFlt d, MrBFlt f)
{
    MrBFlt      p0t;
    
    p0t = (b-d) / (f*b + (b*(1.0-f)-d)*exp((d-b)*t));
    
    return (2.0*log(p0t) + (d-b)*t);
}


/* probability that an individual alive at time t before today has
   no sampled extinct or extant descendants
 */
MrBFlt LnP0_fossil (MrBFlt t, MrBFlt lambda, MrBFlt mu, MrBFlt psi, MrBFlt c1, MrBFlt c2)
{
    MrBFlt other;
    
    // c1 = sqrt(pow(lambda-mu-psi, 2) + 4*lambda*psi);
    // c2 = (-lambda + mu + 2*lambda*rho + psi) / c1;
    other = (exp(-c1 *t) * (1 - c2) - (1 + c2)) / (exp(-c1 *t) * (1 - c2) + (1 + c2));
    
    return log(0.5) + log(lambda + mu + psi + c1 * other) - log(lambda);
}

/* probability that an individual alive at time t before today has
   precisely one sampled extant descendant and no sampled extinct descendant
 */
MrBFlt LnP1_fossil (MrBFlt t, MrBFlt rho, MrBFlt c1, MrBFlt c2)
{
    MrBFlt other;
    
    // c1 = sqrt(pow(lambda-mu-psi, 2) + 4*lambda*psi);
    // c2 = (-lambda + mu + 2*lambda*rho + psi) / c1;
    other = 2.0 * (1- c2*c2) * exp(-c1 *t) + pow(1-c2, 2) * exp(-2 *c1 *t) + pow(1+c2, 2);
    
    return log(4.0) + log(rho) - c1 *t - log(other);
}


/* return which time interval t is in */
int Slice_i (MrBFlt t, MrBFlt *t_f, int sl)
{
    int i = 0;
    assert (t > 0.0 && sl >= 0);

    /* we need some tolerance here, for t[i] < t <= t[i-1] to return i */
    while (t < t_f[i] + BRLENS_MIN/5)
        {
        i++;
        if (i > sl)
            return sl;
        }
    return i;
}

/* probability density of an individual at time t giving rise to an edge
   between time t and t_i with q_i(t_i) = 1
 */
MrBFlt  LnQi_fossil (MrBFlt t, MrBFlt *t_f, int sl, MrBFlt *c1, MrBFlt *c2)
{
    MrBFlt lnq;
    int i = Slice_i (t, t_f, sl);
    
    lnq = log(4.0) +c1[i] *(t_f[i] -t);
    lnq -= 2.0 * log(1 +c2[i] +(1 -c2[i]) *exp(c1[i] *(t_f[i] -t)));
    
    return lnq;
}

/* an individual at time t has no sampled descendants when the process is stopped
   (i.e., at time t_s), with t_i < t <= t_{i-1} (i = 1,..,s)
 */
MrBFlt  LnPi_fossil (MrBFlt t, MrBFlt *t_f, int sl, MrBFlt *c1, MrBFlt *c2, MrBFlt *lambda, MrBFlt *mu, MrBFlt *psi)
{
    MrBFlt other;
    int i = Slice_i (t, t_f, sl);
    
    other = lambda[i] +mu[i] +psi[i] -c1[i] * (1 +c2[i] -(1 -c2[i]) *exp(c1[i] *(t_f[i] -t)))
                                            / (1 +c2[i] +(1 -c2[i]) *exp(c1[i] *(t_f[i] -t)));
    return log(other) - log(2 *lambda[i]);
}


int LnFossilizationPriorPr (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt *sR, MrBFlt *eR, MrBFlt sF, MrBFlt *fR, char *sS)
{
    /* fossilization priors 
     //chi */
    
    if (!strcmp(sS, "FossilTip"))
        return LnFossilizedBDPriorFossilTip (t, clockRate, prob, sR, eR, sF, fR);
    else if (!strcmp(sS, "Random"))
        return LnFossilizedBDPriorRandom    (t, clockRate, prob, sR, eR, sF, fR);
    else if (!strcmp(sS, "Diversity"))
        return LnFossilizedBDPriorDiversity (t, clockRate, prob, sR, eR, sF, fR);
    else
        {
        MrBayesPrint ("%s   Sampling strategy %s for fossilized birth-death process not implemented\n", spacer, sS);
        return (ERROR);
        }
}


/*---------------------------------------------------------------------------------
 |
 |                                          time  ___  0
 |     \                            /            |
 |      \                          /             |___  y2
 |       \                   \    /              |
 |        \                   \  /               |
 |         \                   \/         x3  ___|
 |          \                  /                 |___  y1
 |           \         \      /                  |
 |            \         \    /                   |
 |             \         \  /                    |
 |              \         \/              x2  ___|
 |               \        /                      |
 |                \      /                       |
 |                 \    /                        |
 |                  \  /                         |
 |                   \/                   x1  ___|
 |
 |
 |   T:   oriented tree
 |   b:   birth (speciation) rate
 |   d:   death (extintion) rate
 |   p:   extant sampling rate
 |   q:   fossil sampling rate
 |   n:   number of extant taxa
 |   m:   number of fossil tips
 |
 |
 |                     [p1(x1)]^2    n+m-1            m    q
 |   f(T|tmrca) = ---------------- * prod b*p1(xi) * prod -----  .
 |                [1 - ^p0(x1)]^2    i=2             i=1  p1(yi)
 |
 |   f(tmrca) ~ uniform, gamma, etc (see treeAge).
 |
 ---------------------------------------------------------------------------------*/
int LnFossilizedBDPriorFossilTip (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt *sR, MrBFlt *eR, MrBFlt sF, MrBFlt *fR)
{
    /* special case: upon sampling the lineage is dead and won't produce descendants. Each extinct sample is a tip */
    
    int         i, n, m;
    MrBFlt      x, lambda, rho, psi, tmrca, c1, c2, hatP0;
    TreeNode    *p;
    Model       *mp;
    
    /* sR = lambda-mu-psi, eR = (mu+psi)/lambda, fR = psi/(mu+psi) */
    lambda = sR[0] / (1.0 - eR[0]);
    psi    = lambda * eR[0] * fR[0];
    rho    = sF;
    
    tmrca = t->root->left->nodeDepth / clockRate;
    c1 = sqrt(sR[0]*sR[0] + 4 *lambda *psi);
    c2 = (2 *lambda *rho - sR[0]) / c1;

    /* calculate prior prob of the fbd tree */
    (*prob) = 0.0;

    for (n = m = i = 0; i < t->nNodes -1; i++)
        {
        p = t->allDownPass[i];
        x = p->nodeDepth / clockRate;

        if (p->left != NULL && p->right != NULL)  // internal
            {
            if (p != t->root->left)
                (*prob) += log(lambda) + LnP1_fossil(x, rho, c1, c2);
            }
        else if (p->left == NULL && p->right == NULL)  // tip
            {
            if (p->nodeDepth > 0.0)
                {
                (*prob) += log(psi) - LnP1_fossil(x, rho, c1, c2);
                m++;
                }
            else
                n++;
            }
        }

    /* p_0 (t |psi=0, mu->mu+psi) */
    hatP0 = 1.0 - rho*sR[0] / (rho*lambda + (sR[0] - rho*lambda) * exp(-sR[0]*tmrca));
    
    (*prob) += 2.0 * (LnP1_fossil(tmrca, rho, c1, c2) - log(1 - hatP0));
    
    /* condition on tmrca, calibrations are dealt with separately */
    mp = &modelParams[t->relParts[0]];
    if (t->root->left->isDated == NO)
        (*prob) += mp->treeAgePr.LnPriorProb(tmrca, mp->treeAgePr.priorParams);
    
    /* conversion to labeled tree from oriented tree, constant for a given dataset */
    (*prob) += (n + m - 1) * log(2.0) - LnFactorial(n) - LnFactorial(m);

    return (NO_ERROR);
}


int compare_descending (const void *x, const void *y)
{
    if (*((double *)(x)) > *((double *)(y)))
        return -1;
    else if (*((double *)(x)) < *((double *)(y)))
        return 1;
    else
        return 0;
}
int getTotalRateShifts (Model *mp, MrBFlt *shiftTimes)
{
    int    i, j, sLen, sTotal;

    /* get the total number of rate shifts */
    sTotal = mp->fossilSamplingNum;
    for (i = 0; i < mp->fossilSamplingNum; i++) {
        shiftTimes[i] = mp->fossilSamplingTime[i];
    }
    sLen = sTotal;
    for (i = 0; i < mp->birthRateShiftNum; i++) {
        for (j = 0; j < sLen; j++) {
            if (fabs(mp->birthRateShiftTime[i] - shiftTimes[j]) <= BRLENS_MIN/2) break;
        }
        if (fabs(mp->birthRateShiftTime[i] - shiftTimes[j]) > BRLENS_MIN/2) {
            shiftTimes[sTotal++] = mp->birthRateShiftTime[i];
        }
    }
    sLen = sTotal;
    for (i = 0; i < mp->deathRateShiftNum; i++) {
        for (j = 0; j < sLen; j++) {
            if (fabs(mp->deathRateShiftTime[i] - shiftTimes[j]) <= BRLENS_MIN/2) break;
        }
        if (fabs(mp->deathRateShiftTime[i] - shiftTimes[j]) > BRLENS_MIN/2) {
            shiftTimes[sTotal++] = mp->deathRateShiftTime[i];
        }
    }
    
    /* sort shiftTimes[] in ascending order */
    qsort(shiftTimes, sTotal, sizeof(double), compare_descending);

    return sTotal;
}


/*---------------------------------------------------------------------------------
 |
 |   LnFossilizedBDPriorRandom
 |
 |   Gavryushkina, A., D. Welch, T. Stadler, and A. Drummond. 2014. 
 |       Bayesian inference of sampled ancestor trees for epidemiology and fossil calibration. PLoS Comp. Biol.
 |   Zhang C., T. Stadler, S. Klopfstein, T. A. Heath, and F. Ronquist. 2015.
 |       Total-Evidence Dating under the Fossilized Birth-Death Process. Syst. Biol.
 |
 |
 |                                       0  _____________  t3, rho3
 |     \                            /          |
 |      \                          /           |___  y2
 |       \                   \    /   _________|_________  t2, rho2
 |        \                   \  /             |
 |         \                   \/       x3  ___|
 |          \                  /               |___  y1
 |           \         \      /                |
 |           _\         \    /        _________|_________  t1, rho1
 |             \         \  /                  |
 |              \         \/            x2  ___|
 |               \        /                    |
 |                \      /                     |
 |                 \    /                      |
 |                  \  /                       |
 |                   \/                 x1  ___|_________  t_mrca
 |
 |
 |    sl = 2, t1 > t2 > t3 = 0
 |    E = 2, K = 1, M = 2
 |
 ---------------------------------------------------------------------------------*/
int LnFossilizedBDPriorRandom (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt *sR, MrBFlt *eR, MrBFlt sF, MrBFlt *fR)
{
    /* Fossils are sampled with piecewise constant rates in the past.
       Extant taxa are sampled uniformly at random at present. */
    
    int         i, j, i1, i2, i3,  sl,  K, M, E;
    MrBFlt      x, tmrca, *t_f, *lambda, *mu, *psi, *rho, *netDiver, *turnOver, *sampProp, *c1, *c2, *p_t;
    TreeNode    *p;
    Model       *mp;
    
    mp = &modelParams[t->relParts[0]];
    
    /* time of most recent common ancestor */
    tmrca = t->root->left->nodeDepth / clockRate;
    
    /* alloc memory for time of each shift, t_f[sl] = 0 */
    t_f = (MrBFlt *)SafeMalloc((size_t)(mp->fossilSamplingNum + mp->birthRateShiftNum + mp->deathRateShiftNum +1) * sizeof(MrBFlt));
    if (!t_f)
        {
        MrBayesPrint ("%s   ERROR: Problem allocating t_f in LnFossilizedBDPriorRandom\n", spacer);
        return (ERROR);
        }
    /* get the total rate shifts, sl >= 0 */
    sl = getTotalRateShifts(mp, t_f);
    
    t_f[sl] = 0.0;
    assert (t_f[0] < tmrca);

    /* alloc memory for the other parameters */
    lambda   = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    mu       = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    psi      = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    rho      = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    netDiver = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    turnOver = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    sampProp = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    c1       = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    c2       = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    p_t      = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    if (!lambda || !mu || !psi || !rho || !netDiver || !turnOver || !sampProp || !c1 || !c2 || !p_t)
        {
        MrBayesPrint ("%s   ERROR: Problem allocating memory in LnFossilizedBDPriorRandom\n", spacer);
        free(lambda); free(mu); free(psi); free(rho); free(netDiver); free(turnOver); free(sampProp); free(c1); free(c2); free(p_t);
        return (ERROR);
        }
    
    /* initialization */
    i1 = i2 = i3 = 0;
    for (i = 0; i <= sl; i++)
        {
        if (mp->birthRateShiftNum > 0 && mp->birthRateShiftTime[i1] > t_f[i] +BRLENS_MIN/2)
            i1++;
        if (mp->deathRateShiftNum > 0 && mp->deathRateShiftTime[i2] > t_f[i] +BRLENS_MIN/2)
            i2++;
        if (mp->fossilSamplingNum > 0 && mp->fossilSamplingTime[i3] > t_f[i] +BRLENS_MIN/2)
            i3++;
        netDiver[i] = sR[i1];
        turnOver[i] = eR[i2];
        sampProp[i] = fR[i3];
        }
    for (i = 0; i <= sl; i++)
        {
        /* netDiver = lambda-mu, turnOver = mu/lambda, sampProp = psi/(mu+psi) */
        lambda[i] = netDiver[i] / (1.0 - turnOver[i]);
        mu[i] = lambda[i] * turnOver[i];
        psi[i] = mu[i] * sampProp[i] / (1.0 - sampProp[i]);
        rho[i] = 0.0;
        }
    rho[sl] = sF;

    for (i = sl; i >= 0; i--)
        {
        c1[i] = sqrt(pow(lambda[i]-mu[i]-psi[i], 2) + 4*lambda[i]*psi[i]);
        if (i == sl)
            c2[i] = ((1 - 2* (1-rho[i])) *lambda[i] +mu[i] +psi[i]) /c1[i];
        else
            c2[i] = ((1 - 2* (1-rho[i]) *p_t[i+1]) *lambda[i] +mu[i] +psi[i]) /c1[i];
        if (i > 0)
            p_t[i] = (lambda[i] +mu[i] +psi[i] -c1[i] * (1 +c2[i] -(1 -c2[i]) *exp(c1[i] *(t_f[i] -t_f[i-1])))
                                                      / (1 +c2[i] +(1 -c2[i]) *exp(c1[i] *(t_f[i] -t_f[i-1])))) *0.5/lambda[i];
        else
            p_t[i] = (lambda[i] +mu[i] +psi[i] -c1[i] * (1 +c2[i] -(1 -c2[i]) *exp(c1[i] *(t_f[i] -tmrca)))
                                                      / (1 +c2[i] +(1 -c2[i]) *exp(c1[i] *(t_f[i] -tmrca)))) *0.5/lambda[i];
        }
    
#   ifdef DEBUG_FBDPR
    for (i = 0; i <= sl; i++)
        printf ("%d: lambda=%lf mu=%lf psi=%lf d=%lf r=%lf s=%lf t=%lf rho=%lf\n",
                i, lambda[i], mu[i], psi[i], netDiver[i], turnOver[i], sampProp[i], t_f[i], rho[i]);
    for (i = 0; i <= sl; i++)
        printf ("%d: A=%lf B=%lf p%d(t%d)=%lf\n", i, c1[i], c2[i], i+1, i, p_t[i]);
#   endif
    
    /* calculate prior prob of the fbd tree */
    (*prob) = 0.0;
    
    for (K = M = E = 0, i = 0; i < t->nNodes -1; i++)
        {
        p = t->allDownPass[i];
        x = p->nodeDepth / clockRate;
        
        if (p->left != NULL && p->right != NULL)  // internal
            {
            if (p->left->length > 0.0 && p->right->length > 0.0)
                {
                if (p != t->root->left)
                    (*prob) += log(lambda[Slice_i(x, t_f, sl)]) + LnQi_fossil(x, t_f, sl, c1,c2);
                }
            else
                {
                for (j = 0; j < sl; j++)
                    if (AreDoublesEqual(p->nodeDepth, t_f[j]*clockRate, BRLENS_MIN/5) == YES)  break;
                if (j == sl)      /* fossil ancestor between t[j-1] and t[j] */
                    {
                    (*prob) += log(psi[Slice_i(x, t_f, sl)]);
                    }
                else              /* fossil ancestor at silice time t[j] */
                if (rho[j] > 0.0 && rho[j] < 1.0)
                    {
                    (*prob) += log(rho[j]) - log(1 - rho[j]);
                    }
                K++;              /* number of fossil ancestors */
                }
            }
        else if (p->left == NULL && p->length > 0.0)  // tip
            {
            if (p->nodeDepth > 0.0)
                {
                for (j = 0; j < sl; j++)
                    if (AreDoublesEqual(p->nodeDepth, t_f[j]*clockRate, BRLENS_MIN/5) == YES)  break;
                if (j == sl)      /* fossil tip between t[j-1] and t[j] */
                    {
                    (*prob) += LnPi_fossil(x, t_f, sl, c1,c2, lambda,mu,psi)
                             - LnQi_fossil(x, t_f, sl, c1,c2);
                    (*prob) += log(psi[Slice_i(x, t_f, sl)]);
                    }
                else              /* fossil tip at silice time t[j] */
                    {
                    (*prob) += log(p_t[j+1]);
                    if (rho[j] > 0.0)  (*prob) += log(rho[j]);
                    }
                M++;              /* number of fossil tips */
                }
            else
                {
                if (rho[sl] > 0.0)
                    {
                    (*prob) += log(rho[sl]);
                    }
                E++;              /* number of extant taxa */
                }
            }
        
        for (j = 0; j < sl; j++)  /* degree-two vertices at silice time t_j */
            {
            if (p->length > 0.0)
                if ((p->nodeDepth +BRLENS_MIN/10 < t_f[j]*clockRate) && (t_f[j]*clockRate < p->anc->nodeDepth +BRLENS_MIN/10))
                    {
                    (*prob) += LnQi_fossil(t_f[j], t_f, sl, c1,c2);
                    if (rho[j] < 1.0)  (*prob) += log(1 - rho[j]);
                    }
            }
        }
    
    (*prob) += 2.0 * (LnQi_fossil(tmrca, t_f, sl, c1,c2) - log(1- p_t[0]));
    
    /* condition on tmrca, calibrations are dealt with separately */
    if (t->root->left->isDated == NO)
        (*prob) += mp->treeAgePr.LnPriorProb(tmrca, mp->treeAgePr.priorParams);
    
    /* conversion to labeled tree from oriented tree */
    (*prob) += (M + E - 1) * log(2.0);  // - LnFactorial(E + M + K) (# permutation is constant given data)
    
#   ifdef DEBUG_FBDPR
    printf ("K=%d M=%d E=%d\n", K, M, E);
    printf ("prob=%lf\n", *prob);
#   endif
    
    /* free memory */
    free(t_f);
    free(lambda); free(mu); free(psi); free(rho); free(netDiver); free(turnOver); free(sampProp); free(c1); free(c2); free(p_t);
    
    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
 |
 |   LnFossilizedBDPriorDiversity
 | 
 |   Zhang C., T. Stadler, S. Klopfstein, T. A. Heath, and F. Ronquist. 2015.
 |       Total-Evidence Dating under the Fossilized Birth-Death Process. Syst. Biol.
 |
 ---------------------------------------------------------------------------------*/
int LnFossilizedBDPriorDiversity (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt *sR, MrBFlt *eR, MrBFlt sF, MrBFlt *fR)
{
    /* Fossils are sampled with piecewise constant rates in the past.
       Extant taxa are sampled with prop sF to maximize diversity. */
    
    int         i, j, i1, i2, i3,  sl,  K, M, E;
    MrBFlt      x, tmrca, *t_f, x_min, t_min, M_x, *lambda, *mu, *psi, *rho, *netDiver, *turnOver, *sampProp, *c1, *c2, *p_t;
    TreeNode    *p;
    Model       *mp;
    
    mp = &modelParams[t->relParts[0]];
    
    /* time of most recent common ancestor */
    tmrca = t->root->left->nodeDepth / clockRate;
    
    /* alloc memory for time of each shift, t_f[sl] = 0, t_f[sl-1] = x_cut */
    t_f = (MrBFlt *)SafeMalloc((size_t)(mp->fossilSamplingNum + mp->birthRateShiftNum + mp->deathRateShiftNum +2) * sizeof(MrBFlt));
    if (!t_f)
        {
        MrBayesPrint ("%s   ERROR: Problem allocating t_f in LnFossilizedBDPriorDiversity\n", spacer);
        return (ERROR);
        }
    /* get the total rate shifts, plus 1 to shift psi to 0 */
    sl = getTotalRateShifts(mp, t_f) + 1;
    
    /* get time of youngest fossil and internal node */
    t_min = x_min = tmrca;
    for (i = 0; i < t->nNodes -1; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL && p->nodeDepth > 0.0)  //fossil
            if (t_min > p->nodeDepth / clockRate)
                t_min = p->nodeDepth / clockRate;
        if (p->left != NULL && p->right != NULL &&
            p->left->length > 0.0 && p->right->length > 0.0)
            if (x_min > p->nodeDepth / clockRate)
                x_min = p->nodeDepth / clockRate;
        }

    /* lower the cutoff time if not compatible */
    if (x_min > t_min)
        x_min = t_min;
    if (sl > 1 && t_f[sl-2] < x_min)
        x_min = t_f[sl-2];
    
    t_f[sl] = 0.0;
    t_f[sl-1] = x_min * 0.95; // x_cut
    assert (t_f[0] < tmrca);

    /* alloc memory for the other parameters */
    lambda   = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    mu       = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    psi      = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    rho      = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    netDiver = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    turnOver = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    sampProp = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    c1       = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    c2       = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    p_t      = (MrBFlt *)SafeMalloc((size_t)(sl+1) * sizeof(MrBFlt));
    if (!lambda || !mu || !psi || !rho || !netDiver || !turnOver || !sampProp || !c1 || !c2 || !p_t)
        {
        MrBayesPrint ("%s   ERROR: Problem allocating memory in LnFossilizedBDPriorDiversity\n", spacer);
        free(lambda); free(mu); free(psi); free(rho); free(netDiver); free(turnOver); free(sampProp); free(c1); free(c2); free(p_t);
        return (ERROR);
        }
    
    /* initialization */
    i1 = i2 = i3 = 0;
    for (i = 0; i < sl; i++)
        {
        if (mp->birthRateShiftNum > 0 && mp->birthRateShiftTime[i1] > t_f[i] +BRLENS_MIN/2)
            i1++;
        if (mp->deathRateShiftNum > 0 && mp->deathRateShiftTime[i2] > t_f[i] +BRLENS_MIN/2)
            i2++;
        if (mp->fossilSamplingNum > 0 && mp->fossilSamplingTime[i3] > t_f[i] +BRLENS_MIN/2)
            i3++;
        netDiver[i] = sR[i1];
        turnOver[i] = eR[i2];
        sampProp[i] = fR[i3];
        }
    for (i = 0; i < sl; i++)
        {
        /* netDiver = lambda-mu, turnOver = mu/lambda, sampProp = psi/(mu+psi) */
        lambda[i] = netDiver[i] / (1.0 - turnOver[i]);
        mu[i] = lambda[i] * turnOver[i];
        psi[i] = mu[i] * sampProp[i] / (1.0 - sampProp[i]);
        rho[i] = 0.0;
        }
    lambda[sl] = lambda[sl-1];
    mu[sl]  = mu[sl-1];
    psi[sl] = 0.0;   // psi = 0 in [0, x_cut]
    rho[sl] = 1.0;   // not sF

    for (i = sl; i >= 0; i--)
        {
        c1[i] = sqrt(pow(lambda[i]-mu[i]-psi[i], 2) + 4*lambda[i]*psi[i]);
        if (i == sl)
            c2[i] = ((1 - 2* (1-rho[i])) *lambda[i] +mu[i] +psi[i]) /c1[i];
        else
            c2[i] = ((1 - 2* (1-rho[i]) *p_t[i+1]) *lambda[i] +mu[i] +psi[i]) /c1[i];
        if (i > 0)
            p_t[i] = (lambda[i] +mu[i] +psi[i] -c1[i] * (1 +c2[i] -(1 -c2[i]) *exp(c1[i] *(t_f[i] -t_f[i-1])))
                                                      / (1 +c2[i] +(1 -c2[i]) *exp(c1[i] *(t_f[i] -t_f[i-1])))) *0.5/lambda[i];
        else
            p_t[i] = (lambda[i] +mu[i] +psi[i] -c1[i] * (1 +c2[i] -(1 -c2[i]) *exp(c1[i] *(t_f[i] -tmrca)))
                                                      / (1 +c2[i] +(1 -c2[i]) *exp(c1[i] *(t_f[i] -tmrca)))) *0.5/lambda[i];
        }
    
#   ifdef DEBUG_FBDPR
    for (i = 0; i <= sl; i++)
        printf ("%d: lambda=%lf mu=%lf psi=%lf d=%lf r=%lf s=%lf t=%lf rho=%lf\n",
                i, lambda[i], mu[i], psi[i], netDiver[i], turnOver[i], sampProp[i], t_f[i], rho[i]);
    for (i = 0; i <= sl; i++)
        printf ("%d: A=%lf B=%lf p%d(t%d)=%lf\n", i, c1[i], c2[i], i+1, i, p_t[i]);
#   endif
    
    /* first calculate prob of the fbd tree assuming complete sampling */
    (*prob) = 0.0;
    
    for (K = M = E = 0, i = 0; i < t->nNodes -1; i++)
        {
        p = t->allDownPass[i];
        x = p->nodeDepth / clockRate;
        
        if (p->left != NULL && p->right != NULL)  // internal
            {
            if (p->left->length > 0.0 && p->right->length > 0.0)
                {
                if (p != t->root->left)
                    (*prob) += log(lambda[Slice_i(x, t_f, sl)]) + LnQi_fossil(x, t_f, sl, c1,c2);
                }
            else
                {
                for (j = 0; j < sl; j++)
                    if (AreDoublesEqual(p->nodeDepth, t_f[j]*clockRate, BRLENS_MIN/5) == YES)  break;
                if (j == sl)      /* fossil ancestor between t[j-1] and t[j] */
                    {
                    (*prob) += log(psi[Slice_i(x, t_f, sl)]);
                    }
                else              /* fossil ancestor at silice time t[j] */
                if (rho[j] > 0.0 && rho[j] < 1.0)
                    {
                    (*prob) += log(rho[j]) - log(1 - rho[j]);
                    }
                K++;              /* number of fossil ancestors */
                }
            }
        else if (p->left == NULL && p->length > 0.0)  // tip
            {
            if (p->nodeDepth > 0.0)
                {
                for (j = 0; j < sl; j++)
                    if (AreDoublesEqual(p->nodeDepth, t_f[j]*clockRate, BRLENS_MIN/5) == YES)  break;
                if (j == sl)      /* fossil tip between t[j-1] and t[j] */
                    {
                    (*prob) += LnPi_fossil(x, t_f, sl, c1,c2, lambda,mu,psi)
                             - LnQi_fossil(x, t_f, sl, c1,c2);
                    (*prob) += log(psi[Slice_i(x, t_f, sl)]);
                    }
                else              /* fossil tip at silice time t[j] */
                    {
                    (*prob) += log(p_t[j+1]);
                    if (rho[j] > 0.0)  (*prob) += log(rho[j]);
                    }
                M++;              /* number of fossil tips */
                }
            else
                {
                //  (*prob) += log(rho[sl]);  // rho[sl] == 1
                E++;              /* number of extant taxa */
                }
            }
        
        for (j = 0; j < sl; j++)  /* degree-two vertices at silice time t_j */
            {
            if (p->length > 0.0)
                if ((p->nodeDepth +BRLENS_MIN/10 < t_f[j]*clockRate) && (t_f[j]*clockRate < p->anc->nodeDepth +BRLENS_MIN/10))
                    {
                    (*prob) += LnQi_fossil(t_f[j], t_f, sl, c1,c2);
                    if (rho[j] < 1.0)  (*prob) += log(1 - rho[j]);
                    }
            }
        }
    
    (*prob) += 2.0 * (LnQi_fossil(tmrca, t_f, sl, c1,c2) - log(1- p_t[0]));
    
    /* number of extant taxa not sampled */
    M_x = (int)floor(E/sF + 0.5) - E; /* equal to round(E/sF) plus it is compatible with MS Visual Studio */
    
    /* then calculate the prob of the fbd tree assuming diversified sampling of extant */
    (*prob) += M_x * (log(lambda[sl] * (1.0 - exp((mu[sl]-lambda[sl])*t_f[sl-1]))) - log(lambda[sl] - mu[sl] * exp((mu[sl]-lambda[sl])*t_f[sl-1])));
    
    /* condition on tmrca, calibrations are dealt with separately */
    if (t->root->left->isDated == NO)
        (*prob) += mp->treeAgePr.LnPriorProb(tmrca, mp->treeAgePr.priorParams);
    
    /* conversion to labeled tree from oriented tree */
    (*prob) += (M + E - 1) * log(2.0);  // - LnFactorial(E + M + K) (# permutation is constant given data)
    
#   ifdef DEBUG_FBDPR
    printf ("K=%d M=%d E=%d\n", K, M, E);
    printf ("prob=%lf\n", *prob);
#   endif
    
    /* free memory */
    free(t_f);
    free(lambda); free(mu); free(psi); free(rho); free(netDiver); free(turnOver); free(sampProp); free(c1); free(c2); free(p_t);
    
    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
|
|   LnCoalescencePriorPr
|
|   This function calculates the probability of a tree under the neutral
|   coalescence prior with a (potentially) exponentially growing population.
|   We assume a rooted tree that satisfies the molecular clock constraint. The
|   Tree is labelled as follows:
|
|                                      t_4 ___  
|     \         \         \        /            \
|      \         \         \      /             | 
|   I_4 \         \         \    /              | g_4
|        \         \         \  /               |
|         \         \         \/       t_3 ___  /
|          \         \        /                 \
|           \         \      /                  |
|   I_3      \         \    /                   | g_3
|             \         \  /                    |
|              \         \/            t_2 ___  / 
|               \        /                      \
|                \      /                       |
|   I_2           \    /                        | g_2
|                  \  /                         |
|                   \/                 t_1 ___  /
|    
|   Each interval on the tree is specified by successive coalescence events.
|   These intervals are denoted I_2, I_3, I_4, ..., with the subscript denoting
|   how many lineages exist in that interval. The time of each coalescence event
|   is designated t_1, t_2, t_3, ..., where the subscript denotes the number
|   of lineages that exist after the coalescence (t_3, for instance, would be
|   the time of the coalescence that went from four lineages to three lineages).
|   The duration of the i-th interval is designated g_i.
|
|   The probability of the coalescence tree is:
|   
|   prob = (k_C_2 / (N(t_k + t_k))) * exp(-integral_(from x=t_k, to g_k + t_k)  (k_C_2 / N(x)) dx)
|
|   where N(x) = N(0) * exp(-r*x). For the constant population size case,
|   N(0) = N_e. r is the population growth parameter for the exponentially
|   growing population. Here, theta = N(0) * mu when organisms are haploid and
|   theta = 2 * N(0) * mu when the organism is diploid.
|
|   Below, ct holds the n - 1 coalescence times (t_i, above) sorted from the
|   smallest to the largest. Remember that t_4 < t_3 < t_2, etc. 
|
|   2010-03-23
|   The original function (described above) was incorrect in that it used theta = 2 * N * mu
|   in one part of the equation and theta = 4 * N * mu in another. It is now corrected to
|   consistently use theta = 4 * N * mu, which is the standard for diploid populations
|   and theta = 2 * N * mu for haploid populations. The calculations are the same for the
|   diploid and haploid cases, only the interpretation is different. See, e.g., Felsenstein
|   (2004; Inferring Phylogenies). -- Fredrik.
|
---------------------------------------------------------------------------------*/
int LnCoalescencePriorPr (Tree *t, MrBFlt *prob, MrBFlt theta, MrBFlt growth)
{
    int             i, j, k, nNodes;
    MrBFlt          *ct, tempD, lastCoalescenceTime, coalescenceTime, intervalLength;
    TreeNode        *p;

    /* allocate space for the coalescence times */
    ct = (MrBFlt *)SafeMalloc((size_t)(t->nIntNodes) * sizeof(MrBFlt));
    if (!ct)
        {
        MrBayesPrint ("\n   ERROR: Problem allocating ct\n");
        return (ERROR);
        }

    /* get the coalescence times and put them into a vector */
    for (i=j=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        if (p->anc != NULL)
            ct[j++] = p->nodeDepth;  // Don't divide clockRate here, as mu is already in theta
        }
    nNodes = j;

    /* sort the coalescence times */
    SortMrBFlt (ct, 0, nNodes-1);
    
    /*for (i=0, k=numLocalTaxa; i<nNodes; i++)
        {
        printf ("%4d -- %2d %lf\n", i, k, ct[i]);
        k--;
        }*/
        
    /* calculate probability of the tree */
    if (AreDoublesEqual (growth, 0.0, 0.000001) == YES)
        {
        /* use this if there is no population growth */
        tempD = lastCoalescenceTime = 0.0;
        for (i=0, k=numLocalTaxa; i<nNodes; i++)
            {
            coalescenceTime = ct[i];
            intervalLength = coalescenceTime - lastCoalescenceTime;
            lastCoalescenceTime = ct[i];
            tempD += - (k * (k-1) * intervalLength) / (theta);
            k--;
            }
        (*prob) = (numLocalTaxa - 1) * log(2.0 / theta) + tempD;
        }
    else
        {
        /* use this if the population is growing exponentially */
        tempD = lastCoalescenceTime = 0.0;
        for (i=0, k=numLocalTaxa; i<nNodes; i++)
            {
            coalescenceTime = ct[i];
            intervalLength = coalescenceTime - lastCoalescenceTime;
            tempD += growth * coalescenceTime + (((k * (k-1)) / (theta * growth)) * (exp(growth * lastCoalescenceTime) - exp(growth * coalescenceTime)));
            lastCoalescenceTime = ct[i];
            k--;
            }
        (*prob) = (numLocalTaxa - 1) * log(2.0 / theta) + tempD;
        }

    /* printf ("coal pr = %lf theta = %lf, nNodes = %d, nt = %d tempD = %lf\n", *prob, theta, nNodes, numLocalTaxa, tempD); */

    /* free memory */
    free (ct);
    
    return (NO_ERROR);
}


/*---------------------------------------------------------------------------------
|
|   LnUniformPriorPr
|
|   This function calculates the probability of a calibrated clock tree under the
|   uniform prior probability distribution on node depths. The tree is labeled as
|   follows:
|                                                interval between tip dates 
|    t_0                          t_1      time  ____  0 (interval between t_0 and t_1, duration 0) 
|     \                            /            |  
|      \                  t_2     /             |___   1 (interval between t_1 and t_2)
|       \                   \    /              | 
|        \                   \  /               |
|         \                   \/        0.33 ___|      2 (interval between t_2 and t_3)
|          \        t_3       /                 |___ 
|           \         \      /                  |
|            \         \    /                   | 
|             \         \  /                    |
|              \   t_4   \/             0.67 ___| 
|               \     \  /                      |      3 (interval between t_3 and root)
|                \     \/                       |
|                 \    /                        |      Note that t_4 is irrelevant for intervals
|                  \  /                         |      because we need not have a single coalescent
|                   \/                  1.00 ___|___   event beneath t_4 except for the root.
|    
|   The probability depends on the number of nodes falling in each interval between
|   tip dates, the time interval on which each node depth can vary (if it is interior
|   node number i, it can vary on the interval (t_{i+1},1)). Finally, the probability
|   is multiplied by the number of coalescent histories compatible with the node order.
|
---------------------------------------------------------------------------------*/
MrBFlt LnUniformPriorPr (Tree *t, MrBFlt clockRate)
{
    int         i, j, k, *nLineages=NULL, nDatedTips, nLineagesIn, nLineagesOut, nTips;
    MrBFlt      lnProb, treeAge, *nodeDepths=NULL;
    TreeNode    *p, *root;
    Model       *mp;

    lnProb      = 0.0;
    mp          = &modelParams[t->relParts[0]];
    treeAge     = t->root->left->nodeDepth / clockRate;
    assert (t->root->left->isDated == NO || AreDoublesEqual(treeAge, t->root->left->age, 0.000001) == YES);

    /* Calculate number of tips for convenience */
    nTips   = t->nNodes - t->nIntNodes - 1;

    /* First take tree age into account f(t_0) */
    if (t->root->left->isDated == YES)
        {
        lnProb += 0.0;  /* calibrations are dealt with separately in calling function */
        }
    else
        {
        lnProb += mp->treeAgePr.LnPriorProb(treeAge, mp->treeAgePr.priorParams);
        }

    /* If tree is not calibrated or only root is calibrated, it is easy */
    for (i=j=0; i<t->nNodes-2; i++)
        {
        if (t->allDownPass[i]->isDated == YES)
            j++;
        }
    if (j == 0)
        {
        /* Calculate simple probability f(tau|t_0) */
        lnProb += (nTips - 1.0)*log(2.0) - LnFactorial(nTips) - log(nTips-1.0) - (nTips - 2.0)*log(treeAge);
        assert (lnProb > NEG_INFINITY);
        return lnProb;
        }

    /* We have a tree with interior or tip calibrations */

    /* Color subtrees by assigning an index 1,...,k to the x variable of each node,
       where k is the total number of separate subtrees. A subtree is characterized
       by having all interior nodes being unconstrained and all tips being either
       terminals or constrained interior nodes. The root of a subtree may be the
       root of the tree, or an interior node that is constrained and dated, or an
       interior node that is constrained but not dated. */
    i = 0;
    ColorClusters (t->root->left, &i);

    /* Get the probability for each subtree */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        
        /* Skip unless this is the root of a subtree */
        if (p->anc->anc != NULL && p->isDated != YES)
            continue;

        /* Save the root of the subtree */
        root = p;

        /* Create an array containing the sorted times */

        /* Allocate space for node depths */
        nDatedTips = NumDatedTips (root);
        nodeDepths = (MrBFlt *) SafeRealloc ((void *)nodeDepths, (nDatedTips+1)*sizeof(MrBFlt));

        /* Get the dated node depths and sort them. The call to GetDatedNodeDepths also
           returns the root node depth into nodeDepths, which is convenient. For now, this
           only works for dated tips, not for constrained but undated interior nodes. */
        GetDatedNodeDepths (root, nodeDepths);
        SortMrBFlt (nodeDepths, 0, nDatedTips);   /* use index of left and right in call */

        /* Get probability due to the uniform node depths; we do not use first and last tip depth
           for obvious reasons (see figure above) */
        for (j=1; j<nDatedTips-1; j++)
            lnProb -= log ((root->nodeDepth - nodeDepths[j]) / clockRate);

        /* Get probability due to sorting of interior node depths */
        
        /* First get the potential number of lineages leaving each interval j at time nodeDepths[j+1] */
        nLineages = (int *) SafeRealloc ((void *)nLineages, nDatedTips*sizeof(int));
        for (j=0; j<nDatedTips; j++)
            nLineages[j] = j+1;
        
        /* Subtract interior nodes so that we get the real number of lineages leaving each interval.
           The number of lineages entering each interval is always 1 + the number of lineages
           leaving the previous interval. We utilize the fact here that the last node depth is the
           root node depth and is held in nodeDepths[nDatedTips] */
        for (j=0; j<t->nIntNodes; j++)
            {
            p = t->intDownPass[j];
            if (p->x != root->x || p == root || p->isDated == YES)
                continue;
            for (k=0; k<nDatedTips; k++)
                if (p->nodeDepth < nodeDepths[k+1])
                        nLineages[k]--;
            }
        
        /* Now get the density effect of the sorting constraints */
        for (j=1; j<nDatedTips-1; j++)
            {
            nLineagesIn = nLineages[j-1] + 1;
            if (j==nDatedTips-2)
                nLineagesOut = 2;    /* skip the last segment and jump directly to root */
            else
                nLineagesOut = nLineages[j];
            /* only calculate if ln different from 0 */
            if (nLineagesIn > 1 && nLineagesIn - nLineagesOut >= 1)
                {
                lnProb += LnFactorial (nLineagesIn-1) - LnFactorial(nLineagesOut-1);
                }
            }

        /* Finally get the effect of the number of possible coalescent histories */
        for (j=1; j<nDatedTips; j++)
            {
            nLineagesIn = nLineages[j-1] + 1;
            nLineagesOut = nLineages[j];
            if (nLineagesIn != nLineagesOut)
                {
                lnProb += log(2.0) * (nLineagesIn - nLineagesOut);
                lnProb += LnFactorial (nLineagesOut) + LnFactorial (nLineagesOut-1);
                lnProb -= LnFactorial (nLineagesIn)  + LnFactorial (nLineagesIn-1);
                }
            }

        /* Last but not least, change color of root so that it can be a tip in the next
           subtree (if we wanted to use the colors). */
        root->x = root->anc->x;

        }

    free (nodeDepths);
    free (nLineages);

    assert (lnProb > NEG_INFINITY);

    return lnProb;
}


/*------------------------------------------------------------------------
|
|   NewtonRaphsonBrlen: Find one maximum likelihood branch length using
|      the Newton-Raphson method. This function assumes that the tree is
|      a non-clock tree. For clock trees, you have to optimize the node
|      depths instead.
|
------------------------------------------------------------------------*/
int NewtonRaphsonBrlen (Tree *t, TreeNode *p, int chain)
{
    int         c, i, j, s, k, n, d, division, nIterations, maxNumIterations,
                index, *rateCat, r;
    MrBFlt      vOld, *pi, *cijk, *eigenValues, *catRate, baseRate, length, theRate, tolerance,
                expLambdaV[64], pInvar=0.0, likeI=0.0, sum=0.0, sum1, sum2, beta, w, v, x, y, expBetaV,
                expWV, kappa, bigPi, *ptr;
    CLFlt       *nSitesOfPat, *clP, *clA, *tiP, *tiP1, *tiP2, like, like1, like2, CLsum, CLsum1, CLsum2,
                freq, sumLike1, sumLike2, *lnScaler=NULL, *clInvar=NULL;
    ModelInfo   *m;

    /* TODO: Standard model (also check RES for ascertainment bias) */
    tolerance = 0.001;
    maxNumIterations = 5;
    
    nIterations = 0;
    do {
        /* reset f'(v) and f''(v) sums */
        sumLike1 = sumLike2 = 0.0;
        
        /* cycle over character partitions */
        for (d=0; d<t->nRelParts; d++)
            {
            division = t->relParts[d];

            /* get pointer to model */
            m = &modelSettings[division];
            
            /* get number of model states */
            n = m->numModelStates;
            
            /* find conditional likelihoods */
            clP = m->condLikes[m->condLikeIndex[chain][p->index     ]];
            clA = m->condLikes[m->condLikeIndex[chain][p->anc->index]];
            
            /* get state frequencies */
            pi = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
            /* get base rate */
            theRate = 1.0;
            baseRate = GetRate (division, chain);

            /* get category rates */
            if (m->shape == NULL)
                catRate = &theRate;
            else
                catRate = GetParamSubVals (m->shape, chain, state[chain]);
    
            /* find category frequencies and some additional stuff for the invar model */
            if (m->pInvar == NULL)
                freq =  (CLFlt) (1.0 /  m->numRateCats);
            else
                {
                /* invariable sites model */
                pInvar = * GetParamVals(m->pInvar, chain, state[chain]);
                freq = (CLFlt) ((1.0 - pInvar) /  m->numRateCats);
                baseRate /= (1.0 - pInvar);
                clInvar = m->invCondLikes;
                lnScaler = m->scalers[m->siteScalerIndex[chain]];
                }

            /* find the branch lengths times any correction factor to make them
               in terms of expected number of substitutions per character */
            length = p->length;
            if (m->dataType == DNA || m->dataType == RNA)
                {
                if (m->nucModelId == NUCMODEL_DOUBLET)
                    length *= 2.0;
                else if (m->nucModelId == NUCMODEL_CODON)
                    length *= 3.0;
                }

            /* find nSitesOfPat */
            nSitesOfPat = numSitesOfPat + ((chainId[chain] % chainParams.numChains)*numCompressedChars) + m->compCharStart;
            
            /* get tiProbs for current length; use the scratch location */
            FlipTiProbsSpace(m, chain, p->index);
            m->TiProbs (p, division, chain);
            tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];
            FlipTiProbsSpace(m, chain, p->index);

            /* allocate space for first and second derivatives */
            tiP1 = (CLFlt *) SafeCalloc (2*m->tiProbLength, sizeof (CLFlt));
            if (!tiP1)
                return (ERROR);
            tiP2 = tiP1 + m->tiProbLength;
            
            /* calculate first and second derivatives of P(v): P'(v) and P''(v) */
            index = 0;
            if (m->TiProbs == &TiProbs_Fels || m->TiProbs == &TiProbs_Res)
                {
                /* calculate beta */
                sum = 0.0;
                for (i=0; i<m->numModelStates; i++)
                    for (j=i+1; j<m->numModelStates; j++)
                        sum += pi[i]*pi[j];
                beta = 0.5 / sum;

                /* calculate derivatives */
                for (k=0; k<m->numRateCats; k++)
                    {
                    v = length * catRate[k] * baseRate;
                    expBetaV = exp (- beta * v);
                    for (i=0; i<m->numModelStates; i++)
                        {
                        for (j=0; j<m->numModelStates; j++)
                            {
                            if (i == j)
                                {
                                x = - beta * (1.0 - pi[j]) * expBetaV;
                                tiP1[index] = (CLFlt) x;
                                tiP2[index] = (CLFlt) (- beta * x);
                                }
                            else
                                {
                                x = beta * pi[j] * expBetaV;
                                tiP1[index] = (CLFlt) x;
                                tiP2[index] = (CLFlt) (- beta * x);
                                }
                            index++;
                            }
                        }
                    }
                }
            else if (m->TiProbs == &TiProbs_Hky)
                {
                /* get kappa */
                kappa =  *GetParamVals (m->tRatio, chain, state[chain]);
    
                /* calculate beta */
                sum = 0.0;
                for (i=0; i<m->numModelStates; i++)
                    {
                    for (j=i+1; j<m->numModelStates; j++)
                        {
                        if (j - i == 2)
                            sum += kappa * pi[i] * pi[j];
                        else
                            sum += pi[i] * pi[j];
                        }
                    }
                beta = 0.5 / sum;

                /* calculate derivatives */
                for (k=0; k<m->numRateCats; k++)
                    {
                    v = length * catRate[k] * baseRate;
                    expBetaV = exp (- beta * v);
                    for (i=0; i<m->numModelStates; i++)
                        {
                        for (j=0; j<m->numModelStates; j++)
                            {
                            s = (j + 2) % 4;
                            bigPi = pi[j] + pi[s];
                            w = (1.0 + bigPi * (kappa - 1.0)) * beta;
                            expWV = exp (- w * v);
                            if (i == j)
                                {
                                x = - beta * ((pi[j]/bigPi) - pi[j]) * expBetaV;
                                y = - w * (pi[s]/bigPi) * expWV;
                                tiP1[index] = (CLFlt) (x + y);
                                tiP2[index] = (CLFlt) (- beta * x - w * y);
                                }
                            else if (abs(i-j) == 2)
                                {
                                x = - beta * pi[j] * ((pi[j]/bigPi) - pi[j]) * expBetaV;
                                y = w * (pi[j]/bigPi) * expWV;
                                tiP1[index] = (CLFlt) (x + y);
                                tiP2[index] = (CLFlt) (- beta * x - w * y);
                                }
                            else
                                {
                                x = beta * pi[j] * expBetaV;
                                tiP1[index] = (CLFlt) x;
                                tiP2[index] = (CLFlt) (- beta * x);
                                }
                            index++;
                            }
                        }
                    }
                }
            else if (m->TiProbs == &TiProbs_Gen || m->TiProbs == &TiProbs_GenCov)
                {
                /* get eigenvalues and cijk pointers */
                eigenValues = m->cijks[m->cijkIndex[chain]];
                cijk        = eigenValues + (2 * n);

                /* calculate P'(v) and P''(v) */
                for (k=0; k<m->numRateCats; k++)
                    {
                    v = length * catRate[k];
                    for (s=0; s<n; s++)
                        expLambdaV[s] =  exp(eigenValues[s] * v);
                    ptr = cijk;
                    for (i=0; i<n; i++)
                        {
                        for (j=0; j<n; j++)
                            {
                            sum1 = sum2 = 0.0;
                            for (s=0; s<n; s++)
                                {
                                if (fabs(eigenValues[s]) > 0.000001)
                                    {
                                    x = eigenValues[s] * (*ptr++) * expLambdaV[s];
                                    sum1 += x;
                                    sum2 += eigenValues[s] * x;
                                    }
                                else
                                    ptr += n;
                                }
                            tiP1[index] = (CLFlt) ((sum1 < 0.0) ? 0.0 : sum);
                            tiP2[index] = (CLFlt) ((sum2 < 0.0) ? 0.0 : sum);
                            index++;
                            }
                        }
                    }
                }
            else if (m->TiProbs == &TiProbs_GenCov)
                {
                /* get eigenvalues and cijk pointers */
                eigenValues = m->cijks[m->cijkIndex[chain]];
                cijk        = eigenValues + (2 * n);

                /* calculate P'(v) and P''(v) */
                for (k=0; k<m->numTiCats; k++)
                    {
                    if (m->numRateCats > 1)
                        v = length * catRate[k];
                    else
                        v = length;
                    for (s=0; s<n; s++)
                        expLambdaV[s] =  exp(eigenValues[s] * v);
                    for (i=0; i<n; i++)
                        {
                        for (j=0; j<n; j++)
                            {
                            sum1 = sum2 = 0.0;
                            for (s=0; s<n; s++)
                                {
                                if (fabs(eigenValues[s]) > 0.000001)
                                    {
                                    x = eigenValues[s] * (*cijk++) * expLambdaV[s];
                                    sum1 += x;
                                    sum2 += eigenValues[s] * x;
                                    }
                                }
                            tiP1[index] = (CLFlt) ((sum1 < 0.0) ? 0.0 : sum);
                            tiP2[index] = (CLFlt) ((sum2 < 0.0) ? 0.0 : sum);
                            index++;
                            }
                        }
                    /* get new eigenvalues and cijks */
                    eigenValues += m->cijkLength / m->nCijkParts;
                    cijk         = eigenValues + 2 * n;
                    }
                }
            /* calculate f(v), f'(v) and f''(v) for this partition */
            /* note that the scalers are irrelevant because they disappear when
               we take the derivative of the log likelihood */
            if (m->gibbsGamma == YES)
                {
                /* find rate category index */
                rateCat = m->tiIndex + chain*m->numChars;
                for (c=0; c<m->numChars; c++)
                    {
                    like = like1 = like2 = 0.0;
                    r = rateCat[c];
                    if (r < m->numRateCats)
                        {
                        index = r*m->numModelStates*m->numModelStates;
                        for (j=0; j<n; j++)
                            {
                            CLsum = CLsum1 = CLsum2 = 0.0;
                            for (i=0; i<n; i++, index++)
                                {
                                CLsum += clP[i] * tiP[index];
                                CLsum1 += clP[i] * tiP1[index];
                                CLsum2 += clP[i] * tiP2[index];
                                }
                            like += CLsum * (CLFlt) pi[j] * clA[j];
                            like1 += CLsum1 * (CLFlt) pi[j] * clA[j];
                            like2 += CLsum2 * (CLFlt) pi[j] * clA[j];
                            }
                        like *= freq;
                        sumLike1 += nSitesOfPat[c] * (like1 / like);
                        sumLike2 += nSitesOfPat[c] * (like2 * like - like1 * like1) / (like * like);
                        }
                    clP += n;
                    }
                }
            else
                {
                for (c=0; c<m->numChars; c++)
                    {
                    like = like1 = like2 = 0.0;
                    index = 0;
                    for (k=0; k<m->numTiCats; k++)
                        {
                        for (j=0; j<n; j++)
                            {
                            CLsum = CLsum1 = CLsum2 = 0.0;
                            for (i=0; i<n; i++, index++)
                                {
                                CLsum += clP[i] * tiP[index];
                                CLsum1 += clP[i] * tiP1[index];
                                CLsum2 += clP[i] * tiP2[index];
                                }
                            like += CLsum * (CLFlt) pi[j] * clA[j];
                            like1 += CLsum1 * (CLFlt) pi[j] * clA[j];
                            like2 += CLsum2 * (CLFlt) pi[j] * clA[j];
                            }
                        clP += n;
                        }
                    like *= freq;
                    if (m->pInvar != NULL)
                        {
                        /* get right like; like1 and like2 not affected */;
                        for (i=0; i<n; i++)
                            likeI += (*clInvar++) * pi[i];
                        likeI *= pInvar;
                        if (lnScaler[c] < -200.0)
                            {
                            /* we are not going to be able to exponentiate the scaling factor */
                            if (likeI > 1E-70)
                                {
                                /* forget about like; it is going to be insignificant compared to likeI */
                                like = (CLFlt) likeI;
                                }
                            else
                                {
                                /* treat likeI as if 0.0, that is, ignore it completely */
                                like = like + (CLFlt)(0.0);
                                }
                            }
                        else    /* take both likeI and like into account */
                            like = like + (CLFlt) (likeI * exp (-lnScaler[c]));
                        }
                    sumLike1 += nSitesOfPat[c] * (like1 / like);
                    sumLike2 += nSitesOfPat[c] * (like2 * like - like1 * like1) / (like * like);
                    }
                }
            free (tiP1);
            }
        vOld = p->length;
        p->length -= sumLike2 / sumLike1;
        nIterations++;
        } while (fabs(p->length - vOld) > tolerance && nIterations < maxNumIterations);

    return (NO_ERROR);
}


void NodeToNodeDistances (Tree *t, TreeNode *fromNode)
{
    int             i;
    TreeNode        *p;
    
    /* set all distances to 0.0 and also set marks on all nodes to NO */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->x = NO;
        p->d = 0.0;
        }
        
    /* find distances, and mark path, below "fromNode" */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p == fromNode)
            {
            p->x = YES;
            }
        if (p->left != NULL && p->right != NULL && p->anc != NULL)
            {
            if (p->left->x == YES)
                {
                p->x = YES;
                p->d = p->left->d + p->left->length;
                }
            else if (p->right->x == YES)
                {
                p->x = YES;
                p->d = p->right->d + p->right->length;
                }
            }
        }
        
    /* find all other distances */
    for (i=t->nNodes-1; i>=0; i--)
        {
        p = t->allDownPass[i];
        if (p->anc == NULL)
            {
            if (p == fromNode)
                p->d = 0.0;
            else
                p->d = p->left->d + p->left->length;
            }
        else
            {
            if (p->x == NO)
                {
                p->d = p->anc->d + p->length;
                }
            }
        }
}


int NumCppEvents (Param *p, int chain)
{
    int         i, *nEvents, sumEvents;

    nEvents = p->nEvents[2*chain+state[chain]];
    
    sumEvents = 0;
    for (i=0; i<2*numLocalTaxa-2; i++)
        sumEvents += nEvents[i];

    return (sumEvents);
}


/*----------------------------------------------------------------------
|
|   OpenNewMBPrintFile: Open a file the first time for printing
|
------------------------------------------------------------------------*/
FILE *OpenNewMBPrintFile (char *fileName)
{
    int     overWrite;
    FILE    *fp;

    /* Open file, use noWarn to determine if the user should be prompted
       to have the file over-written or appended. */
    if (noWarn == YES)
        {
        /* overwrite or append file, if already present */
        if ((fp = TestOpenTextFileR(fileName)) != NULL)
            {
            SafeFclose (&fp);
            if (autoOverwrite == NO)
                {
                MrBayesPrint ("%s   Appending to file \"%s\"\n", spacer, fileName); 
                return (OpenTextFileA(fileName));
                }
            else 
                MrBayesPrint ("%s   Overwriting file \"%s\"\n", spacer, fileName);      
            }
            return (OpenTextFileW(fileName));
        }   
    else
        {
        /* prompt user if file is already present */
        if ((fp = TestOpenTextFileR(fileName)) != NULL)
            {
            SafeFclose (&fp);
            MrBayesPrint ("\n");
            MrBayesPrint ("%s   File \"%s\" already exists\n", spacer, fileName);
            overWrite = WantTo ("Overwrite information in this file");
            
            if (overWrite == YES)
                {
                MrBayesPrint ("%s   Overwriting file \"%s\"\n", spacer, fileName);
                return (OpenTextFileW(fileName));
                }
            else
                {
                MrBayesPrint ("%s   Appending to file \"%s\"\n", spacer, fileName);
                return (OpenTextFileA(fileName));
                }
            }

        else
            {
            /* file is not present */
            return (OpenTextFileW(fileName));
            }
        }
}


int PickProposal (RandLong *seed, int chainIndex)
{
    MrBFlt      ran;
    int         i;

    ran = RandomNumber(seed);
    
    for (i=0; usedMoves[i]->cumProposalProb[chainIndex] <= ran; i++);
        
    return i;
}


/* Calculate positive selection probabilities */
int PosSelProbs (TreeNode *p, int division, int chain)
{
    int             c, j, k, nStates;
    MrBFlt          catLike, *like, *bs, *omegaCatFreq, *omega,
                    posProb, *ps, sum;
    CLFlt           **clP;
    ModelInfo       *m;
    
    /* find model partition */
    m = &modelSettings[division];

    /* allocate space for conditional likelihood pointer array and site likelihood array */
    clP  = (CLFlt **) calloc (m->numOmegaCats, sizeof(CLFlt *));
    like = (MrBFlt *) calloc (m->numOmegaCats, sizeof(MrBFlt));
    if (!clP || !like)
        {
        MrBayesPrint ("%s   ERROR: Out of memory in PosSelProbs\n", spacer);
        return (ERROR);
        }
    
    /* number of states */
    nStates = m->numModelStates;

    /* find conditional likelihoods */
    clP[0] = m->condLikes[m->condLikeIndex[chain][p->index]];
    for (k=1; k<m->numOmegaCats; k++)
        clP[k] = clP[0] + k*m->numModelStates*m->numChars;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* find category frequencies */
    omegaCatFreq = GetParamSubVals (m->omega, chain, state[chain]);
    
    /* get category omegas */
    omega = GetParamVals (m->omega, chain, state[chain]);

    /* find posSelProbs */
    ps = posSelProbs + m->compCharStart;
    for (c=0; c<m->numChars; c++)
        {
        sum = 0.0;
        for (k=0; k<m->numOmegaCats; k++)
            {
            like[k] = 0.0;
            catLike = 0.0;
            for (j=0; j<nStates; j++)
                catLike += clP[k][j] * bs[j];
            like[k] = catLike * omegaCatFreq[k];
            sum += like[k];
            clP[k] += nStates;
            }
        posProb = 0.0;
        for (k=0; k<m->numOmegaCats; k++)
            {
            if (omega[k] > 1.0)
                posProb += like[k] / sum;
            }
        ps[c] = posProb;
        }

    free (clP);
    free (like);
    
    return NO_ERROR;
}


#if defined (SSE_ENABLED)
/* Calculate positive selection probabilities (SSE version) */
int PosSelProbs_SSE (TreeNode *p, int division, int chain)
{
    int             i, c1, c2, j, k, nStates;
    CLFlt           **catLike, *siteLike;
    MrBFlt          *bs, *omegaCatFreq, *omega,
                    posProb, *ps;
    __m128          m1, m2, *clPtr, **clP, mSiteLike, *mCatLike;
    ModelInfo       *m;
    
    /* find model partition */
    m = &modelSettings[division];
    
    /* number of states */
    nStates = m->numModelStates;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* find category frequencies */
    omegaCatFreq = GetParamSubVals (m->omega, chain, state[chain]);
    
    /* get category omegas */
    omega = GetParamVals (m->omega, chain, state[chain]);
    /* allocate space for category likelihood arrays */
    catLike = (CLFlt **) calloc (m->numOmegaCats, sizeof(CLFlt *));
    mCatLike = (__m128 *) calloc (m->numOmegaCats, sizeof(__m128));
    if (!catLike || !mCatLike)
        {
        MrBayesPrint ("%s   ERROR: Out of memory in PosSelProbs_SSE\n", spacer);
        return (ERROR);
        }

    /* find conditional likelihood pointers */
    clPtr = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP_SSE;
    for (k=0; k<m->numOmegaCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * nStates;
        catLike[k] = (CLFlt *) (&(mCatLike[k]));
        }
    siteLike = (CLFlt *) (&mSiteLike);
    
    /* find posSelProbs */
    ps = posSelProbs + m->compCharStart;
    for (c1=c2=0; c1<m->numVecChars; c1++)
        {
        mSiteLike = _mm_setzero_ps ();
        for (k=0; k<m->numOmegaCats; k++)
            {
            mCatLike[k] = _mm_setzero_ps();
            m1 = _mm_setzero_ps ();
            for (j=0; j<nStates; j++)
                {
                m2 = _mm_mul_ps (clP[k][j], _mm_set1_ps ((CLFlt)bs[j]));
                m1 = _mm_add_ps (m1, m2);
                }
            mCatLike[k] = _mm_mul_ps (m1, _mm_set1_ps ((CLFlt)omegaCatFreq[k]));
            mSiteLike = _mm_add_ps (mSiteLike, mCatLike[k]);
            clP[k] += nStates;
            }

        for (i=0; i<m->numFloatsPerVec && c2 < m->numChars; ++i, ++c2)
            {
            posProb = 0.0;
            for (k=0; k<m->numOmegaCats; k++)
                {
                if (omega[k] > 1.0)
                    {
                    posProb += catLike[k][i] / siteLike[i];
                    }
                }
            ps[c2] = posProb;
            }
        }
    
    free (catLike);
    free (mCatLike);
    
    return NO_ERROR;
}
#endif


/* Calculate omega values for each site */
int SiteOmegas (TreeNode *p, int division, int chain)
{
    int             c, j, k, nStates;
    MrBFlt          catLike, *like, *bs, *omegaCatFreq, *omega,
                    siteOmega, *ps, sum;
    CLFlt           **clP;
    ModelInfo       *m;
    
    /* find model partition */
    m = &modelSettings[division];
    
    /* allocate space for conditional likelihood pointer array and site likelihood array */
    clP  = (CLFlt **) calloc (m->numOmegaCats, sizeof(CLFlt *));
    like = (MrBFlt *) calloc (m->numOmegaCats, sizeof(MrBFlt));
    if (!clP || !like)
        {
        MrBayesPrint ("%s   ERROR: Out of memory in SiteOmegas\n", spacer);
        return (ERROR);
        }
    
    /* number of states */
    nStates = m->numModelStates;

    /* find conditional likelihoods */
    clP[0] = m->condLikes[m->condLikeIndex[chain][p->index]];
    for (k=1; k<m->numOmegaCats; k++)
        clP[k] = clP[0] + k*m->numModelStates*m->numChars;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* find category frequencies */
    omegaCatFreq = GetParamSubVals (m->omega, chain, state[chain]);
    
    /* get category omegas */
    omega = GetParamVals (m->omega, chain, state[chain]);

    /* find site omegas (using posSelProbs space) */
    ps = posSelProbs + m->compCharStart;
    for (c=0; c<m->numChars; c++)
        {
        sum = 0.0;
        for (k=0; k<m->numOmegaCats; k++)
            {
            like[k] = 0.0;
            catLike = 0.0;
            for (j=0; j<nStates; j++)
                catLike += clP[k][j] * bs[j];
            like[k] = catLike * omegaCatFreq[k];
            sum += like[k];
            clP[k] += nStates;
            }
        siteOmega = 0.0;
        for (k=0; k<m->numOmegaCats; k++)
            {
            siteOmega += (like[k]/sum) * omega[k];
            }
        ps[c] = siteOmega;
        }

    free (clP);
    free (like);
    
    return NO_ERROR;
}


#if defined (SSE_ENABLED)
/* Calculate omega values for each site (SSE version) */
int SiteOmegas_SSE (TreeNode *p, int division, int chain)
{
    int             i, c1, c2, j, k, nStates;
    CLFlt           **catLike, *siteLike;
    MrBFlt          *bs, *omegaCatFreq, *omega,
                    siteOmega, *ps;
    __m128          m1, m2, *clPtr, **clP, mSiteLike, *mCatLike;
    ModelInfo       *m;
    
    /* find model partition */
    m = &modelSettings[division];
    
    /* number of states */
    nStates = m->numModelStates;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* find category frequencies */
    omegaCatFreq = GetParamSubVals (m->omega, chain, state[chain]);
    
    /* get category omegas */
    omega = GetParamVals (m->omega, chain, state[chain]);

    /* allocate space for category likelihood arrays */
    catLike = (CLFlt **) calloc (m->numOmegaCats, sizeof(CLFlt *));
    mCatLike = (__m128 *) calloc (m->numOmegaCats, sizeof(__m128));
    if (!catLike || !mCatLike)
        {
        MrBayesPrint ("%s   ERROR: Out of memory in SiteOmegas_SSE\n", spacer);
        return (ERROR);
        }
    
    /* find conditional likelihood pointers */
    clPtr = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP_SSE;
    for (k=0; k<m->numOmegaCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * nStates;
        catLike[k] = (CLFlt *) (&(mCatLike[k]));
        }
    siteLike = (CLFlt *) (&mSiteLike);
    
    /* find site omegas (using posSelProbs space) */
    ps = posSelProbs + m->compCharStart;
    for (c1=c2=0; c1<m->numVecChars; c1++)
        {
        mSiteLike = _mm_setzero_ps ();
        for (k=0; k<m->numOmegaCats; k++)
            {
            mCatLike[k] = _mm_setzero_ps();
            m1 = _mm_setzero_ps ();
            for (j=0; j<nStates; j++)
                {
                m2 = _mm_mul_ps (clP[k][j], _mm_set1_ps ((CLFlt)bs[j]));
                m1 = _mm_add_ps (m1, m2);
                }
            mCatLike[k] = _mm_mul_ps (m1, _mm_set1_ps ((CLFlt)omegaCatFreq[k]));
            mSiteLike = _mm_add_ps (mSiteLike, mCatLike[k]);
            clP[k] += nStates;
            }
        
        for (i=0; i<m->numFloatsPerVec && c2 < m->numChars; ++i, ++c2)
            {
            siteOmega = 0.0;
            for (k=0; k<m->numOmegaCats; k++)
                {
                siteOmega += (catLike[k][i] / siteLike[i]) * omega[k];
                }
            ps[c2] = siteOmega;
            }
        }
    
    free (catLike);
    free (mCatLike);
    
    return NO_ERROR;
}
#endif


/*----------------------------------------------------------------------
|
|   PreparePrintFiles: Prepare .t, .p, and .mcmc files for printing
|
------------------------------------------------------------------------*/
int PreparePrintFiles (void)
{
    int         i, n, previousResults, oldAutoOverwrite, oldNoWarn;
    char        localFileName[100], fileName[220], bkupName[220];
    FILE        *tempFile;

#if defined (MPI_ENABLED)
    if (proc_id != 0)
        return (NO_ERROR);
#endif

    oldNoWarn        = noWarn;
    oldAutoOverwrite = autoOverwrite;

    /* Allocate space for file pointers */
    if (memAllocs[ALLOC_FILEPOINTERS] == YES)
        {
        MrBayesPrint ("%s   File pointers already allocated in PreparePrintFiles\n", spacer);
        return ERROR;
        }
    fpMcmc = NULL;
    fpSS = NULL;
    fpParm = NULL;
    fpTree = NULL;  
    fpParm = (FILE **) SafeCalloc (chainParams.numRuns, sizeof (FILE *));
    if (fpParm == NULL)
        {
        MrBayesPrint ("%s   Could not allocate fpParm in PreparePrintFiles\n", spacer);
        return ERROR;
        }
    memAllocs[ALLOC_FILEPOINTERS] = YES;
    fpTree = (FILE ***) SafeCalloc (chainParams.numRuns, sizeof (FILE **));
    if (fpTree == NULL)
        {
        MrBayesPrint ("%s   Could not allocate fpTree in PreparePrintFiles\n", spacer);
        return ERROR;
        }
    fpTree[0] = (FILE **) SafeCalloc (numTrees*chainParams.numRuns, sizeof (FILE *));
    if (fpTree[0] == NULL)
        {
        MrBayesPrint ("%s   Could not allocate fpTree[0] in PreparePrintFiles\n", spacer);
        return ERROR;
        }
    for (i=1; i<chainParams.numRuns; i++)
        fpTree[i] = fpTree[0] + i*numTrees;

    /* Get root of local file name */
    strcpy (localFileName, chainParams.chainFileName);

    /* Determine whether to overwrite files */
    if (noWarn == NO)
        {
        previousResults = NO;
        if (chainParams.mcmcDiagn == YES)
            {
            sprintf (fileName, "%s.mcmc", localFileName);
            if ((tempFile = TestOpenTextFileR(fileName)) != NULL)
                {
                fclose(tempFile);
                previousResults = YES;
                }
            }
        for (n=0; n<chainParams.numRuns; n++)
            {
            if (chainParams.numRuns == 1)
                sprintf (fileName, "%s.p", localFileName);
            else
                sprintf (fileName, "%s.run%d.p", localFileName, n+1);
            if ((tempFile = TestOpenTextFileR(fileName)) != NULL)
                {
                fclose(tempFile);
                previousResults = YES;
                }

            for (i=0; i<numTrees; i++)
                {
                if (numTrees == 1 && chainParams.numRuns == 1)
                    sprintf (fileName, "%s.t", localFileName);
                else if (numTrees > 1 && chainParams.numRuns == 1)
                    sprintf (fileName, "%s.tree%d.t", localFileName, i+1);
                else if (numTrees == 1 && chainParams.numRuns > 1)
                    sprintf (fileName, "%s.run%d.t", localFileName, n+1);
                else
                    sprintf (fileName, "%s.tree%d.run%d.t", localFileName, i+1, n+1);
                if ((tempFile = TestOpenTextFileR(fileName)) != NULL)
                    {
                    fclose(tempFile);
                    previousResults = YES;
                    }
                }
            }
        if (previousResults == YES)
            {
            MrBayesPrint("\n");
            MrBayesPrint("%s   There are results from a previous run saved using the same filename(s).\n", spacer);
            if (WantTo("Do you want to overwrite these results") == YES)
                {
                MrBayesPrint("\n");
                noWarn = YES;
                autoOverwrite = YES;
                }
            else
                {
                MrBayesPrint("\n");
                MrBayesPrint("%s   Please specify a different file name before running the mcmc analysis.\n", spacer);
                MrBayesPrint("%s      You can do that using 'mcmc filename=<name>'. You can also move or\n", spacer);
                MrBayesPrint("%s      rename the old result files.\n", spacer);
                return ERROR;
                }
            }
        }

    /* Prepare the .mcmc file */
    if (chainParams.mcmcDiagn == YES)
        {
        sprintf (fileName, "%s.mcmc", chainParams.chainFileName);
        if ((fpMcmc = OpenNewMBPrintFile (fileName)) == NULL)
            {
            noWarn = oldNoWarn;
            autoOverwrite = oldAutoOverwrite;
            return (ERROR);
            }
        }
    
    /* Prepare the .p and .t files */
    for (n=0; n<chainParams.numRuns; n++)
        {
        if (chainParams.numRuns == 1)
            sprintf (fileName, "%s.p", localFileName);
        else
            sprintf (fileName, "%s.run%d.p", localFileName, n+1);
        if ((fpParm[n] = OpenNewMBPrintFile (fileName)) == NULL)
            {
            noWarn = oldNoWarn;
            autoOverwrite = oldAutoOverwrite;
            return (ERROR);
            }

        for (i=0; i<numTrees; i++)
            {
            if (numTrees == 1 && chainParams.numRuns == 1)
                sprintf (fileName, "%s.t", localFileName);
            else if (numTrees > 1 && chainParams.numRuns == 1)
                sprintf (fileName, "%s.tree%d.t", localFileName, i+1);
            else if (numTrees == 1 && chainParams.numRuns > 1)
                sprintf (fileName, "%s.run%d.t", localFileName, n+1);
            else
                sprintf (fileName, "%s.tree%d.run%d.t", localFileName, i+1, n+1);
            if ((fpTree[n][i] = OpenNewMBPrintFile (fileName)) == NULL)
                {
                noWarn = oldNoWarn;
                autoOverwrite = oldAutoOverwrite;
                return (ERROR);
                }
            }
        }

    /* Prepare the .ss file */
    if (chainParams.isSS == YES)
        {
        sprintf (fileName, "%s.ss", chainParams.chainFileName);
        if ((fpSS = OpenNewMBPrintFile (fileName)) == NULL)
            {
            noWarn = oldNoWarn;
            autoOverwrite = oldAutoOverwrite;
            return (ERROR);
            }
        }

    /* Remove previous chekpoint file if present */
    sprintf (fileName, "%s%s.ckp", workingDir, chainParams.chainFileName);
    strcpy (bkupName, fileName);
    strcat (bkupName, "~");
    remove (bkupName);
    rename (fileName, bkupName);
    
#   if defined (PRINT_DUMP)
    fpDump = (FILE **) SafeCalloc (chainParams.numRuns, sizeof (FILE *));

    for (n=0; n<chainParams.numRuns; n++)
        {
        if (chainParams.numRuns == 1)
            sprintf (fileName, "%s.dump", localFileName);
        else
            sprintf (fileName, "%s.run%d.dump", localFileName, n+1);
        
        if ((fpDump[n] = OpenNewMBPrintFile (fileName)) == NULL)
            return (ERROR);
        }
#   endif

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   PrintAncStates_Bin: print ancestral states after final pass
|       Binary model with or without rate variation
|
-----------------------------------------------------------------*/
int PrintAncStates_Bin (TreeNode *p, int division, int chain)
{
    int             c, i, k;
    MrBFlt          *bs, freq;
    CLFlt           *clFP, *cL, sum, **clP;
    char            *tempStr;
    int             tempStrSize = TEMPSTRSIZE;
    ModelInfo       *m;
    
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    /* find model settings for this division */
    m = &modelSettings[division];

    /* find state frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* find frequencies of rate categories */
    freq =  1.0 /  m->numRateCats;
    
    /* find the conditional likelihoods from the final pass */
    clFP = m->condLikes[m->condLikeScratchIndex[p->index]];

    clP = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] =  clFP;
        clFP += m->numChars * m->numModelStates;
        }

    /* find the preallocated working space */
    cL = m->ancStateCondLikes;
    
    /* cycle over the compressed characters */
    for (c=0; c<m->numChars; c++)
        {
        cL[0] = cL[1] = 0.0;
        for (k=0; k<m->numRateCats; k++)
            {
            cL[0] += clP[k][0];
            cL[1] += clP[k][1];
            clP[k]+=2;
            }
        cL[0] *= (CLFlt) (bs[0] * freq);
        cL[1] *= (CLFlt) (bs[1] * freq);
        sum = cL[0] + cL[1];
        assert (cL[0]==cL[0]);
        assert (cL[1]==cL[1]);
        assert (sum<9999999999999999999999999999999.0);
        cL[0] /= sum;
        cL[1] /= sum;
        assert (cL[0]==cL[0]);
        assert (cL[1]==cL[1]);
        cL += 2;
        }

    /* print the resulting conditional likelihoods cycling over uncompressed chars */
    for (c=0; c<numChar; c++)
        {
        if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum] != division+1)
            continue;
        i = compCharPos[c] - m->compCharStart;
        cL = m->ancStateCondLikes + (i*2);
        SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(cL[0]));
        if (AddToPrintString (tempStr) == ERROR) return (ERROR);
        SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(cL[1]));
        if (AddToPrintString (tempStr) == ERROR) return (ERROR);
        }

    free (tempStr);
    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   PrintAncStates_Gen: print ancestral states after final pass
|       General model with or without rate variation
|
-----------------------------------------------------------------*/
int PrintAncStates_Gen (TreeNode *p, int division, int chain)
{
    int             c, i, k, nStates, hasPInvar, nGammaCats;
    const int       *rateCat;
    MrBFlt          *bsVals;
    CLFlt           *cL, sum, pInvar=0.0, freq, f, bs[64];
    const CLFlt     *clFP, *clInvar=NULL, *lnScaler, **clP;
    char            *tempStr, *printedChar;
    int             tempStrSize = TEMPSTRSIZE;
    ModelInfo       *m;
    
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    if (!strcmp(modelParams[division].nucModel,"Codon") || !strcmp(modelParams[division].nucModel,"Protein") || !strcmp(modelParams[division].nucModel,"Doublet"))
        {
        assert (modelParams[division].dataType == DNA || modelParams[division].dataType == RNA);
        /* Note that we can have matrix with Protein datatype which is not and should not be covered here */
        printedChar = (char *) SafeMalloc (numChar*sizeof(char));
        }
    else
        {
        printedChar = NULL;
        }

        
    /* find model settings for this division */
    m = &modelSettings[division];

    /* find number of states */
    nStates = m->numModelStates;
    
    /* find state frequencies */
    bsVals = GetParamSubVals (m->stateFreq, chain, state[chain]);
    for (i=0; i<nStates; i++)
        bs[i] = (CLFlt) bsVals[i];

    /* find invar cond likes */
    if (m->pInvar == NULL)
        hasPInvar = NO;
    else
        {
        hasPInvar = YES;
        clInvar = m->invCondLikes;
        pInvar = (CLFlt) *GetParamVals (m->pInvar, chain, state[chain]); 
        }

    /* find number of rate categories */
    nGammaCats = m->numRateCats;

    /* find frequencies of rate categories (only relevant if gibbsGamma == NO) */
    freq = ((CLFlt)1.0 - pInvar) / nGammaCats;
    
    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find rate category index */
    rateCat = m->tiIndex + chain*m->numChars;

    /* find the conditional likelihoods from the final pass */
    clFP = m->condLikes[m->condLikeScratchIndex[p->index]];
    
    /* find the preallocated working space */
    cL = m->ancStateCondLikes;
    
    /* cycle over the compressed characters */
    if (m->gibbsGamma == YES)
        {
        for (c=0; c<m->numChars; c++)
            {
            sum = 0.0;
            if (rateCat[c] < nGammaCats)
                {
                for (i=0; i<nStates; i++)
                    {
                    cL[i] = *(clFP++) * bs[i];
                    sum += cL[i];
                    }
                clInvar += nStates;
                }
            else
                {
                for (i=0; i<nStates; i++)
                    {
                    cL[i] = *(clInvar++) * bs[i];
                    sum += cL[i];
                    }
                clFP += nStates;
                }
            for (i=0; i<nStates; i++)
                cL[i] /= sum;
            cL += nStates;
            }
        }
    else
        {
        /* find conditional likelihood pointers */
        clP = (const CLFlt**)m->clP;
        for (k=0; k<m->numRateCats; k++)
            {
            clP[k] = clFP;
            clFP += m->numChars * m->numModelStates;
            }
        for (c=0; c<m->numChars; c++)
            {
            for (i=0; i<nStates; i++)
                cL[i] = 0.0;
            for (k=0; k<nGammaCats; k++)
                {
                for (i=0; i<nStates; i++)
                    cL[i] += *(clP[k]++);
                }
            for (i=0; i<nStates; i++)
                cL[i] *= bs[i];

            if (hasPInvar == YES && pInvar > 0)
                {
                sum = 0.0;
                for (i=0; i<nStates; i++)
                    sum += clInvar[i];
                if (sum > 0.0)
                    {
                    if (lnScaler[c] < -100.0)
                        {
                        /* only invar cond likes are relevant */
                        for (i=0; i<nStates; i++)
                            cL[i] = clInvar[i] * bs[i] * pInvar;
                        }
                    else
                        {
                        f = (CLFlt) exp(lnScaler[c]) * freq;
                        for (i=0; i<nStates; i++)
                            cL[i] = clInvar[i] * bs[i] * pInvar + cL[i] * f;
                        }
                    }
                clInvar += nStates;
                }
            
            sum = 0.0;
            for (i=0; i<nStates; i++)
                sum += cL[i];
            assert (sum > 0.0);
            for (i=0; i<nStates; i++)
                cL[i] /= sum;
            cL += nStates;
            }
        }

    /* print the resulting conditional likelihoods cycling over uncompressed chars */
    if (printedChar)
        for (c=0; c<numChar; c++)
            printedChar[c] = NO;

    for (c=0; c<numChar; c++)
        {
        if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum] != division+1 || (printedChar &&
            printedChar[c] == YES))
            continue;
        i = compCharPos[c] - m->compCharStart;
        cL = m->ancStateCondLikes + (i*nStates);
        if (printedChar)
            {
            for (i=c+1; i<numChar; i++)
                if (charInfo[c].charId == charInfo[i].charId)
                    printedChar[i] = YES;
            }
        for (i=0; i<nStates; i++)
            {
            SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(cL[i]));
            if (AddToPrintString (tempStr) == ERROR) return (ERROR);
            }
        }
    free (tempStr);
    free (printedChar);
    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   PrintAncStates_NUC4: print ancestral states after final pass
|       4-state nucleotide model with or without rate variation
|
-----------------------------------------------------------------*/
int PrintAncStates_NUC4 (TreeNode *p, int division, int chain)
{
    int             c, i, k, *rateCat, hasPInvar, nGammaCats;
    MrBFlt          *bsVals;
    CLFlt           *cL, sum, pInvar=0.0, bs[4], freq, f;
    const CLFlt     *clFP, *clInvar=NULL, *lnScaler,**clP;
    char            *tempStr;
    int             tempStrSize = TEMPSTRSIZE;
    ModelInfo       *m;
    
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    /* find model settings for this division */
    m = &modelSettings[division];

    /* find state frequencies */
    bsVals = GetParamSubVals (m->stateFreq, chain, state[chain]);
    bs[A] = (CLFlt) bsVals[A];
    bs[C] = (CLFlt) bsVals[C];
    bs[G] = (CLFlt) bsVals[G];
    bs[T] = (CLFlt) bsVals[T];

    /* find invar cond likes */
    if (m->pInvar == NULL)
        hasPInvar = NO;
    else
        {
        hasPInvar = YES;
        clInvar = m->invCondLikes;
        pInvar = (CLFlt) *GetParamVals (m->pInvar, chain, state[chain]); 
        }

    /* find number of rate categories */
    nGammaCats = m->numRateCats;

    /* find frequencies of rate categories (only relevant if gibbsGamma == NO) */
    if (hasPInvar == NO)
        freq =  (CLFlt) 1.0 /  nGammaCats;
    else
        freq = ((CLFlt)1.0 - pInvar) / nGammaCats;
    
    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find rate category index */
    rateCat = m->tiIndex + chain*m->numChars;

    /* find the conditional likelihoods from the final pass */
    clFP = m->condLikes[m->condLikeScratchIndex[p->index]];
    
    /* find the preallocated working space */
    cL = m->ancStateCondLikes;
    
    /* cycle over the compressed characters */
    if (m->gibbsGamma == YES)
        {
        for (c=0; c<m->numChars; c++)
            {
            if (rateCat[c] < nGammaCats)
                {
                cL[A] = clFP[A] * (CLFlt) bs[A];
                cL[C] = clFP[C] * (CLFlt) bs[C];
                cL[G] = clFP[G] * (CLFlt) bs[G];
                cL[T] = clFP[T] * (CLFlt) bs[T];
                }
            else
                {
                cL[A] = clInvar[A] * (CLFlt) bs[A];
                cL[C] = clInvar[C] * (CLFlt) bs[C];
                cL[G] = clInvar[G] * (CLFlt) bs[G];
                cL[T] = clInvar[T] * (CLFlt) bs[T];
                }
            sum = cL[A] + cL[C] + cL[G] + cL[T];
            cL[A] /= sum;
            cL[C] /= sum;
            cL[G] /= sum;
            cL[T] /= sum;
            clInvar += 4;
            clFP += 4;
            cL += 4;
            }
        }
    else
        {
        /* find conditional likelihood pointers */
        clP = (const CLFlt**)m->clP;
        for (k=0; k<m->numRateCats; k++)
            {
            clP[k] = clFP;
            clFP += m->numChars * m->numModelStates;
            }

        for (c=0; c<m->numChars; c++)
            {
            cL[A] = cL[C] = cL[G] = cL[T] = 0.0;
            for (k=0; k<nGammaCats; k++)
                {
                cL[A] += clP[k][A];
                cL[C] += clP[k][C];
                cL[G] += clP[k][G];
                cL[T] += clP[k][T];
                clP[k] += 4;
                }
            cL[A] *= bs[A];
            cL[C] *= bs[C];
            cL[G] *= bs[G];
            cL[T] *= bs[T];
            if (hasPInvar == YES)
                {
                sum = clInvar[A] + clInvar[C] + clInvar[G] + clInvar[T];
                if (sum > 0.0)
                    {
                    if (lnScaler[c] < -100.0)
                        {
                        /* only invar cond likes are relevant */
                        cL[A] = clInvar[A] * bs[A] * pInvar;
                        cL[C] = clInvar[C] * bs[C] * pInvar;
                        cL[G] = clInvar[G] * bs[G] * pInvar;
                        cL[T] = clInvar[T] * bs[T] * pInvar;
                        }
                    else
                        {
                        f = (CLFlt)exp(lnScaler[c]) * freq;
                        cL[A] = clInvar[A] * bs[A] * pInvar + cL[A] * f;
                        cL[C] = clInvar[C] * bs[C] * pInvar + cL[C] * f;
                        cL[G] = clInvar[G] * bs[G] * pInvar + cL[G] * f;
                        cL[T] = clInvar[T] * bs[T] * pInvar + cL[T] * f;
                        }
                    }
                clInvar += 4;
                }
            sum = cL[A] + cL[C] + cL[G] + cL[T];
            cL[A] /= sum;
            cL[C] /= sum;
            cL[G] /= sum;
            cL[T] /= sum;
            cL += 4;
            }
        }

    /* print the resulting conditional likelihoods cycling over uncompressed chars */
    for (c=0; c<numChar; c++)
        {
        if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum] != division+1)
            continue;
        i = compCharPos[c] - m->compCharStart;
        cL = m->ancStateCondLikes + (i*4);
        SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(cL[A]));
        if (AddToPrintString (tempStr) == ERROR) return ERROR;
        SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(cL[C]));
        if (AddToPrintString (tempStr) == ERROR) return ERROR;
        SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(cL[G]));
        if (AddToPrintString (tempStr) == ERROR) return ERROR;
        SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(cL[T]));
        if (AddToPrintString (tempStr) == ERROR) return ERROR;
        }
    free (tempStr);
    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   PrintAncStates_Std: print ancestral states after final pass
|       Standard model with or without rate variation
|
-----------------------------------------------------------------*/
int PrintAncStates_Std (TreeNode *p, int division, int chain)
{
    int             c, i, j, k, s, nStates, numReps;
    MrBFlt          *bsBase, *bs, freq;
    CLFlt           *clFP, *cL, sum,** clP;
    char            *tempStr;
    int             tempStrSize = TEMPSTRSIZE;
    ModelInfo       *m;
    
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    /* find model settings for this division */
    m = &modelSettings[division];

    /* find state frequencies, base index */
    bsBase = GetParamStdStateFreqs (m->stateFreq, chain, state[chain]);
    
    /* find the conditional likelihoods from the final pass */
    clFP = m->condLikes[m->condLikeScratchIndex[p->index]];

    numReps=0;
    for (c=0; c<m->numChars; c++)
        {
        if (m->nStates[c] == 2)
                numReps += m->numBetaCats * 2;
            else
                numReps += m->nStates[c];
        }

    /* find conditional likelihood pointers */
    clP = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clFP;
        clFP += numReps;
        }

    /* find the preallocated working space */
    cL = m->ancStateCondLikes;
    
    /* cycle over the compressed characters */
    for (c=0; c<m->numChars; c++)
        {
        nStates = m->nStates[c];
        bs = bsBase + m->bsIndex[c];

        for (s=0; s<nStates; s++)
            cL[s] = 0.0;

        if (nStates == 2)
            {
            freq = 1.0 / (m->numBetaCats * m->numRateCats);
            for (i=0; i<m->numBetaCats; i++)
                {
                for (k=0; k<m->numRateCats; k++)
                    {
                    for (s=0; s<nStates; s++)
                        cL[s] += clP[k][s] * (CLFlt)(bs[s] * freq);
                    clP[k] += nStates;
                    }
                bs += nStates;
                }
           }
        else
            {
            freq = 1.0 / (m->numRateCats);
            for (k=0; k<m->numRateCats; k++)
                {
                for (s=0; s<nStates; s++)
                    cL[s] += clP[k][s] * (CLFlt)(bs[s] * freq);
                clP[k] += nStates;
                }
            }

        sum = 0.0;
        for (s=0; s<nStates; s++)
            sum += cL[s];

        assert (sum != 0.0);

        for (s=0; s<nStates; s++) 
            cL[s] /= sum;

        cL += nStates;
        }

    /* print the resulting conditional likelihoods cycling over uncompressed chars */
    for (c=0; c<numChar; c++)
        {
        if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum] != division+1)
            continue;
        
        k = compCharPos[c] - m->compCharStart;
        for (i=j=0; i<k; i++)
            j += m->nStates[i];
        cL = m->ancStateCondLikes + j;

        for (i=0; i<m->nStates[k]; i++)
            {
            SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(cL[i]));
            if (AddToPrintString (tempStr) == ERROR) return (ERROR);
            }
        }
    free (tempStr);
    return NO_ERROR;
}


/*-----------------------------------------------------------------------
|
|   PrintCheckPoint: Print checkpoint to file
|
------------------------------------------------------------------------*/
int PrintCheckPoint (int gen)
{
    int         i, j, k, k1, nErrors=0, run, chn, nValues, tempStrSize = TEMPSTRSIZE,
                hasEvents, *intValue, id, oldPrecision;
    char        bkupFileName[220], oldBkupFileName[220], ckpFileName[220], *tempString=NULL;
    MrBFlt      *value, clockRate;
    Param       *p = NULL, *subParm = NULL;
    Tree        *t;
    FILE        *fp = NULL;
    MCMCMove    *mv;
    /*ModelInfo *m = NULL;*/

#if defined (MPI_ENABLED)
    int         sumErrors=0,ierror;
    MrBFlt      r, sum;
#endif

    /* use high precision for checkpointing */
    oldPrecision = precision;
    precision = 15;

    /* allocate tempString */
    if ((tempString = (char *) SafeCalloc (tempStrSize, sizeof(char))) == NULL)
        nErrors++;

#if defined (MPI_ENABLED)
    if (proc_id == 0)
        {
#endif
    if (nErrors == 0)
        {
        /* figure out check-point file names */
        sprintf(ckpFileName, "%s.ckp", chainParams.chainFileName);
        strcpy (bkupFileName, ckpFileName);
        strcat (bkupFileName, "~");
        strcpy (oldBkupFileName, bkupFileName);
        strcat (oldBkupFileName, "~");

        /* shift check-point files */
        remove (oldBkupFileName);
        rename (bkupFileName, oldBkupFileName);
        rename (ckpFileName, bkupFileName);

        /* create new ckp file */
        if ((fp = OpenTextFileW (ckpFileName)) == NULL)
            {
            MrBayesPrint ("%s   Problem opening checkpoint file\n", spacer);
            nErrors++;
            }
        }
    
#if defined (MPI_ENABLED)
        } /* end of if (proc_id == 0)*/
#endif

    ERROR_TEST2("",free(tempString),return(ERROR));
    
    /* write file header */
    MrBayesPrintf (fp, "#NEXUS\n[ID: %s]\n[generation: %d]\n", stamp, gen);

    if (chainParams.isSS == YES)
        {
        /* dump to .ckp file current step contribution */
        MrBayesPrintf (fp, "[SsAcumulators:");
#       if defined (MPI_ENABLED)
        for (j=0; j<chainParams.numRuns ; j++)
            {
            if (stepAcumulatorSS[j]==0)
                r=0;
            else
                r = log (stepAcumulatorSS[j]) + stepScalerSS[j];
            ierror = MPI_Reduce (&r,&sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (ierror != MPI_SUCCESS)
                {
                MrBayesPrint ("%s   Problem with MPI_Reduce\n", spacer);
                return ERROR;
                }
            if (proc_id == 0)
                {
                MrBayesPrintf (fp, " %.4f", sum);
                }
            }
#       else
        for (j=0; j<chainParams.numRuns ; j++)
            {
            MrBayesPrintf (fp, " %.4f", log (stepAcumulatorSS[j]) + stepScalerSS[j]);
            }
#       endif
        MrBayesPrintf (fp, "]\n");
        }

#if defined (MPI_ENABLED)
if (proc_id == 0)
    {
#endif
    MrBayesPrintf (fp, "\nbegin trees;\n\ttranslate\n");

    /* remove old bkup file ('~~') */
    remove (oldBkupFileName);

    /* write translate block */
    for (i=0; i<numLocalTaxa; i++)
        {
        if (i == numLocalTaxa - 1)
            MrBayesPrintf (fp, "      %2d %s;\n", i+1, localTaxonNames[i]);
        else
            MrBayesPrintf (fp, "      %2d %s,\n", i+1, localTaxonNames[i]);
        }
   
#if defined (MPI_ENABLED)
    }
#endif

    /* allocate space for print string */
    printStringSize = tempStrSize;
    printString = (char *) SafeCalloc ((size_t)printStringSize, sizeof(char));
    if (!printString)
        nErrors++;
    else
        strcpy(printString,"");

    ERROR_TEST2("Memory allocation error",free(tempString),return(ERROR));
    /*
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   Memory allocation error on at least one processor\n", spacer);
        free (tempString);
        return ERROR;
        }
#   else
    if (nErrors > 0)
        {
        free (tempString);
        return ERROR;
        }
#   endif
*/
    /* print trees (but not species trees) */
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType != P_BRLENS && p->paramType != P_TOPOLOGY)
            continue;
        if (p->paramType == P_TOPOLOGY && p->subParams[0] != p)
            continue;
        hasEvents = NO;
        for (j=0; j<numLocalChains; j++)
            {
            t = GetTree (p, j, state[j]);
            /* write the tree preamble */
            if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\ttree %s", t->name) == ERROR)
                nErrors++;
            if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                nErrors++;
            if (p->paramType == P_BRLENS && p->nSubParams > 0)
                {
                for (k=0; k<p->nSubParams; k++)
                    {
                    subParm = p->subParams[k];
                    if (subParm->paramType == P_CPPEVENTS)
                        {
                        hasEvents = YES;
                        if (SafeSprintf (&tempString, &tempStrSize, " [&E %s]", subParm->name) == ERROR) nErrors++;
                        if (nErrors == 0 && AddToPrintString (tempString) == ERROR) nErrors++;
                        }
                    if (nErrors == 0 && (subParm->paramType == P_CPPEVENTS || subParm->paramType == P_TK02BRANCHRATES ||
                                         subParm->paramType == P_IGRBRANCHRATES || subParm->paramType == P_MIXEDBRCHRATES))
                        {
                        if (subParm->paramType == P_MIXEDBRCHRATES)
                            {
                            id = *GetParamIntVals(subParm, j, state[j]);
                            if (SafeSprintf (&tempString, &tempStrSize, " [&B %s %d]", subParm->name, id) == ERROR) nErrors++;
                            }
                        else
                            if (SafeSprintf (&tempString, &tempStrSize, " [&B %s]", subParm->name) == ERROR) nErrors++;
                        if (nErrors == 0 && AddToPrintString (tempString) == ERROR) nErrors++;
                        }
                    }
                }

            if (t->isRooted == YES && t->isClock == NO)
                SafeSprintf (&tempString, &tempStrSize, " = ");
            else if (t->isRooted == YES && t->isClock == YES)
                {
                clockRate = *GetParamVals(modelSettings[p->relParts[0]].clockRate, j, state[j]);
                SafeSprintf (&tempString, &tempStrSize, " = [&R] [&clockrate = %s] ", MbPrintNum(clockRate));
                }
            else /* if (t->isRooted == NO) */
                SafeSprintf (&tempString, &tempStrSize, " = ");
            if (nErrors == 0 && AddToPrintString (tempString) == ERROR) nErrors++;
            /* write the tree in (extended) Newick format */
            if (nErrors == 0)
                {
                if (p->paramType == P_TOPOLOGY)
                    WriteNoEvtTreeToPrintString (t->root->left, j, p, NO, t->isRooted);
                else if (hasEvents == NO)
                    WriteNoEvtTreeToPrintString (t->root->left, j, p, YES, t->isRooted);
                else
                    WriteEventTreeToPrintString (t->root->left, j, p, YES);
                if (AddToPrintString (";\n") == ERROR)
                    nErrors++;
                }
            }
        MrBayesPrintf (fp, "%s", printString);
#if defined (MPI_ENABLED)
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Print error on at least one processor\n", spacer);
            goto errorExit;
            }
        if (PrintMPISlaves(fp) == ERROR)
            goto errorExit;
#else
        if (nErrors > 0)
            goto errorExit;
#endif
        strcpy (printString, "");
        strcpy (tempString, "");
        }
    MrBayesPrintf (fp, "end;\n\n");

    /* print species trees */
    if (strcmp(modelParams[0].topologyPr,"Speciestree") == 0)
        {
        /* get the first species tree */
        for (i=0; i<numParams; i++)
            {
            p = &params[i];
            if (p->paramType == P_SPECIESTREE)
                break;
            }
        t = GetTree(p, 0, state[0]);

#if defined (MPI_ENABLED)
if (proc_id == 0)
        {
#endif
        /* write the block header and translate block */
        MrBayesPrintf (fp, "\nbegin trees;\n");
        PrintTranslateBlock (fp, t);
#if defined (MPI_ENABLED)
        }
#endif

        for (j=0; j<numLocalChains; j++)
            {
            t = GetTree (p, j, state[j]);

            /* write the tree preamble */
            if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\ttree %s", t->name) == ERROR)
                nErrors++;
            if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                nErrors++;
                
            subParm = modelSettings[p->relParts[0]].popSize;
            if (subParm->nValues > 1)
                {
                if (SafeSprintf (&tempString, &tempStrSize, " [&N %s]", subParm->name) == ERROR) nErrors++;
                if (nErrors == 0 && AddToPrintString (tempString) == ERROR) nErrors++;
                }

            clockRate = *GetParamVals(modelSettings[p->relParts[0]].clockRate, j, state[j]);
            SafeSprintf (&tempString, &tempStrSize, " = [&R] [&clockrate = %s] ", MbPrintNum(clockRate));
            if (nErrors == 0 && AddToPrintString (tempString) == ERROR) nErrors++;

            /* write the tree in (extended) Newick format */
            if (nErrors == 0)
                {
                WriteNoEvtTreeToPrintString (t->root->left, j, p, YES, t->isRooted);
                if (AddToPrintString (";\n") == ERROR)
                nErrors++;
                }
            }
        MrBayesPrintf (fp, "%s", printString);
#if defined (MPI_ENABLED)
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Print error on at least one processor\n", spacer);
            goto errorExit;
            }
        if (PrintMPISlaves(fp) == ERROR)
            goto errorExit;
#else
        if (nErrors > 0)
            goto errorExit;
#endif
        strcpy (printString, "");
        strcpy (tempString, "");
        MrBayesPrintf (fp, "end;\n\n");
        }
            
    /* start startvals block */
    MrBayesPrintf (fp,"begin mrbayes;\n");
    MrBayesPrintf (fp, "\tstartvals\n");

    /* first print topology values */
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_TOPOLOGY)
            {
            for (j=0; j<numLocalChains; j++)
                {
                t = GetTree (p, j, state[j]);
                run = (chainId[j] / chainParams.numChains) + 1;
                chn = (chainId[j] % chainParams.numChains) + 1;
                if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\t\t%s(%d,%d)=%s\n", p->name, run, chn, t->name) == ERROR)
                    nErrors++;
                if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                    nErrors++;
                }
            MrBayesPrintf (fp, "%s", printString);
#   if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Print error on at least one processor\n", spacer);
                goto errorExit;
                }
            if (PrintMPISlaves(fp) == ERROR)
                goto errorExit;
#   else
            if (nErrors > 0)
                goto errorExit;
#   endif
            strcpy (printString, "");
            strcpy (tempString, "");
            }
        }

    /* now print branch lengths and relaxed clock parameters */
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_BRLENS)
            {
            for (j=0; j<numLocalChains; j++)
                {
                t = GetTree (p, j, state[j]);
                run = (chainId[j] / chainParams.numChains) + 1;
                chn = (chainId[j] % chainParams.numChains) + 1;
                if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\t\t%s(%d,%d)=%s\n", p->name, run, chn, t->name) == ERROR)
                    nErrors++;
                if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                    nErrors++;
                for (k=0; k<p->nSubParams; k++)
                    {
                    if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\t\t%s(%d,%d)=%s\n", p->subParams[k]->name, run, chn, t->name) == ERROR)
                        nErrors++;
                    if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                        nErrors++;
                    }
                }
            MrBayesPrintf (fp, "%s", printString);
#   if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Print error on at least one processor\n", spacer);
                goto errorExit;
                }
            if (PrintMPISlaves(fp) == ERROR)
                goto errorExit;
#   else
            if (nErrors > 0)
                goto errorExit;
#   endif
            strcpy (printString, "");
            strcpy (tempString, "");
            }
        }

    /* now print species tree and population size parameters */
    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_SPECIESTREE)
            {
            for (j=0; j<numLocalChains; j++)
                {
                t = GetTree (p, j, state[j]);
                run = (chainId[j] / chainParams.numChains) + 1;
                chn = (chainId[j] % chainParams.numChains) + 1;
                if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\t\t%s(%d,%d)=%s\n", p->name, run, chn, t->name) == ERROR)
                    nErrors++;
                if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                    nErrors++;
                if (modelSettings[p->relParts[0]].popSize->nValues > 1)
                if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\t\t%s(%d,%d)=%s\n", modelSettings[p->relParts[0]].popSize->name, run, chn, t->name) == ERROR)
                    nErrors++;
                if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                    nErrors++;
                }
            MrBayesPrintf (fp, "%s", printString);
#   if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Print error on at least one processor\n", spacer);
                goto errorExit;
                }
            if (PrintMPISlaves(fp) == ERROR)
                goto errorExit;
#   else
            if (nErrors > 0)
                goto errorExit;
#   endif
            strcpy (printString, "");
            strcpy (tempString, "");
            }
        }

    /* now print param values */
    for (i=0; i<numPrintParams; i++)
        {
        p = printParam[i];
        for (j=0; j<numLocalChains; j++)
            {
            run = (chainId[j] / chainParams.numChains) + 1;
            chn = (chainId[j] % chainParams.numChains) + 1;
            if (p->paramType == P_PI && modelSettings[p->relParts[0]].dataType != STANDARD)
                {
                value = GetParamSubVals (p, j, state[j]);
                nValues = p->nSubValues;
                }
            else
                {
                value = GetParamVals (p, j, state[j]);
                nValues = p->nValues;
                }
            if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\t\t%s(%d,%d)=(%.15le", p->name, run, chn, value[0]) == ERROR)
                nErrors++;
            if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                nErrors++;
            for (k=1; k<nValues; k++)
                {
                if (nErrors==0 && SafeSprintf (&tempString, &tempStrSize, ",%.15le", value[k]) == ERROR)
                    nErrors++;
                if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                    nErrors++;
                }
            /* print int values if present */
            if (p->nIntValues > 0)
                {
                intValue = GetParamIntVals (p, j, state[j]);
                nValues  = p->nIntValues;
                for (k=0; k<nValues; k++)
                    {
                    if (nErrors==0 && SafeSprintf (&tempString, &tempStrSize, ",%d", intValue[k]) == ERROR)
                        nErrors++;
                    if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                        nErrors++;
                    }
                }
            /* print extra params for symdir multistate */
            if (p->nSympi > 0)
                {
                value = GetParamStdStateFreqs (p, j, state[j]);
                if (p->hasBinaryStd == YES)
                    value += 2 * modelSettings[p->relParts[0]].numBetaCats;
                for (k=0; k<p->nSympi; k++)
                    {
                    for (k1=0; k1<p->sympinStates[k]; k1++)
                        {
                        if (nErrors==0 && SafeSprintf (&tempString, &tempStrSize, ",%.15le", *value++) == ERROR)
                            nErrors++;
                        if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                            nErrors++;
                        }
                    }
                }
            /* print extra params for omega */
            if (p->paramType == P_OMEGA)
                {
                value = GetParamSubVals (p, j, state[j]);
                nValues = p->nSubValues/2;
                for (k=0; k<nValues; k++)
                    {
                    if (nErrors==0 && SafeSprintf (&tempString, &tempStrSize, ",%.15le", value[k]) == ERROR)
                        nErrors++;
                    if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                        nErrors++;
                    }
                }
            if (nErrors == 0 && AddToPrintString (")\n") == ERROR)
                nErrors++;
            }
        MrBayesPrintf (fp, "%s", printString);
#if defined (MPI_ENABLED)
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Print error on at least one processor\n", spacer);
            goto errorExit;
            }
        if (PrintMPISlaves(fp) == ERROR)
            goto errorExit;
#else
        if (nErrors > 0)
            goto errorExit;
#endif
        strcpy (printString, "");
        strcpy (tempString, "");
        }
    /* end startvals statement */
    MrBayesPrintf (fp, "\t;\n");

    /* print tuning parameters of moves */
    MrBayesPrintf (fp, "\tpropset\n");
    for (i=0; i<numUsedMoves; i++)
        {
        mv = usedMoves[i];
        if (mv->moveType->Autotune == NULL)
            continue;   /* tuning parameter(s) do not change */
        for (j=0; j<numLocalChains; j++)
            {
            run = (chainId[j] / chainParams.numChains) + 1;
            chn = (chainId[j] % chainParams.numChains) + 1;
            /* format is:   <move_name>$<tuning_param_name>(<run>,<chain>)=<number> */
            if (nErrors == 0 && SafeSprintf (&tempString, &tempStrSize, "\t\t%s$%s(%d,%d)=%.15le\n",
                mv->name, mv->moveType->shortTuningName[0], run, chn, mv->tuningParam[chainId[j]][0]) == ERROR)
                nErrors++;
            if (nErrors == 0 && AddToPrintString (tempString) == ERROR)
                nErrors++;
            }
        MrBayesPrintf (fp, "%s", printString);
#if defined (MPI_ENABLED)
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Print error on at least one processor\n", spacer);
            goto errorExit;
            }
        if (PrintMPISlaves(fp) == ERROR)
            goto errorExit;
#else
        if (nErrors > 0)
            goto errorExit;
#endif
        strcpy (printString, "");
        strcpy (tempString, "");
        }

    /* end propset statement */
    MrBayesPrintf (fp, "\t;\n"); 
    
    /* end mrbayes block */
    MrBayesPrintf (fp, "end;\n\n");

    /* change precision back */
    precision = oldPrecision;
    
    SafeFclose (&fp);
    free (tempString);
    free (printString);
    printStringSize = 0;
    return (NO_ERROR);

errorExit:
    precision = oldPrecision;
    free (tempString);
    free (printString);
    printString = NULL;
    printStringSize = 0;
    SafeFclose (&fp);
    return (ERROR);
}


/*----------------------------------------------------------------------
|
|   PrintMCMCDiagnosticsToFile: Print acceptance ratios, swapping
|      frequencies, and convergence diagnostics to file.
|
------------------------------------------------------------------------*/
int PrintMCMCDiagnosticsToFile (int curGen)
{
    int         i, j, n;
    MCMCMove    *theMove;
    char        *diagnstat;

    /* Simply print header if curGen == 0 */
    if (curGen == 0)
        {
        // MrBayesPrintf (fpMcmc, "[LEGEND:]\n");
        MrBayesPrintf (fpMcmc, "[ID: %s]\n", stamp);
        MrBayesPrintf (fpMcmc, "[   Gen                --  Generation]\n");
        if (chainParams.allChains == YES)
            MrBayesPrintf (fpMcmc, "[   <name1>(<name2>)$acc_run2_chn3 --  Acceptance rate of move <name1> changing parameter <name2> in run 2, chain 3]\n");
        else /* if (chainParams.allChains == NO) */
            MrBayesPrintf (fpMcmc, "[   <name1>(<name2>)$acc_run2      --  Acceptance rate of move <name1> changing parameter 1 in run 2]\n");
        if (chainParams.numChains > 1 && chainParams.numRuns > 1)
            MrBayesPrintf (fpMcmc, "[   Swap(1<>2)$acc_run3            --  Acceptance rate of swaps between chains 1 and 2 in run 3]\n");
        else if (chainParams.numChains > 1 && chainParams.numRuns == 1)
            MrBayesPrintf (fpMcmc, "[   Swap(1<>2)$acc                 --  Acceptance rate of swaps between chains 1 and 2]\n");
        if (chainParams.diagnStat == AVGSTDDEV)
            diagnstat = "Average";
        else
            diagnstat = "Maximum";
        if (chainParams.numRuns > 1 && numTopologies == 1 && chainParams.allComps == NO)
            MrBayesPrintf (fpMcmc, "[   StdDev(s)                      --  %s standard deviation of split frequencies]\n", diagnstat);
        else if (chainParams.numRuns > 1 && numTopologies == 1 && chainParams.allComps == YES)
            MrBayesPrintf (fpMcmc, "[   StdDev(s)(run1-run2)           --  %s standard deviation of split frequencies between runs 1 and 2]\n", diagnstat);
        else if (chainParams.numRuns > 1 && numTopologies > 1 && chainParams.allComps == NO)
            MrBayesPrintf (fpMcmc, "[   StdDev(s.tree1)                --  %s standard deviation of split frequencies for topology 1]\n", diagnstat);
        else if (chainParams.numRuns > 1 && numTopologies > 1 && chainParams.allComps == YES)
            MrBayesPrintf (fpMcmc, "[   StdDev(s.tree1)(run1-run2)     --  %s standard deviation of split frequencies between runs 1 and 2 for topology 1]\n", diagnstat);

        MrBayesPrintf (fpMcmc, "Gen");
        for (n=0; n<chainParams.numRuns; n++)
            {
            if (chainParams.allChains == YES)
                {
                for (i=0; i<chainParams.numChains; i++)
                    {
                    for (j=0; j<numUsedMoves; j++)
                        {
                        theMove = usedMoves[j];
                        MrBayesPrintf (fpMcmc, "\t%s(%s)$acc_run%d_chn%d", theMove->moveType->shortName,
                                theMove->parm->name, n+1, i+1);
                        if (theMove->moveType->Autotune != NULL && chainParams.autotune == YES)
                            MrBayesPrintf (fpMcmc, "\t%s(%s)$%s_run%d_chn%d", theMove->moveType->shortName,
                                theMove->parm->name, theMove->moveType->shortTuningName[0], n+1, i+1);
                        }
                    }
                }
            else
                {
                for (i=0; i<numUsedMoves; i++)
                    {
                    theMove = usedMoves[i];
                    if (chainParams.numRuns == 1)
                        {
                        MrBayesPrintf (fpMcmc, "\t%s(%s)$acc", theMove->moveType->shortName, theMove->parm->name);
                        if (theMove->moveType->Autotune != NULL && chainParams.autotune == YES)
                            MrBayesPrintf (fpMcmc, "\t%s(%s)$%s", theMove->moveType->shortName, theMove->parm->name, theMove->moveType->shortTuningName[0]);
                        }
                    else
                        {
                        MrBayesPrintf (fpMcmc, "\t%s(%s)$acc_run%d", theMove->moveType->shortName, theMove->parm->name, n+1);
                        if (theMove->moveType->Autotune != NULL && chainParams.autotune == YES)
                            MrBayesPrintf (fpMcmc, "\t%s(%s)$%s_run%d", theMove->moveType->shortName, theMove->parm->name, theMove->moveType->shortTuningName[0], n+1);
                        }
                    }
                }
            if (chainParams.numChains > 1)
                {
                for (i=0; i<chainParams.numChains; i++)
                    {
                    for (j=i+1; j<chainParams.numChains; j++)
                        {
                        if (chainParams.numRuns == 1)
                            MrBayesPrintf (fpMcmc, "\tSwap(%d<>%d)$acc", i+1, j+1);
                        else
                            MrBayesPrintf (fpMcmc, "\tSwap(%d<>%d)$acc(%d)", i+1, j+1, n+1);
                        }
                    }
                }
            }

        if (chainParams.numRuns > 1)
            {
            for (n=0; n<numTopologies; n++)
                {
                if (numTopologies == 1)
                    {
                    if (chainParams.diagnStat == AVGSTDDEV)
                        MrBayesPrintf (fpMcmc, "\tAvgStdDev(s)");
                    else
                        MrBayesPrintf (fpMcmc, "\tMaxStdDev(s)");
                    }
                else
                    {
                    if (chainParams.diagnStat == AVGSTDDEV)
                        MrBayesPrintf (fpMcmc, "\tAvgStdDev(s.tree%d)", n+1);
                    else
                        MrBayesPrintf (fpMcmc, "\tMaxStdDev(s.tree%d)", n+1);
                    }

                if (chainParams.allComps == YES)
                    {
                    for (i=0; i<chainParams.numRuns; i++)
                        {
                        for (j=i+1; j<chainParams.numRuns; j++)
                            {
                            if (numTopologies == 1)
                                {
                                if (chainParams.diagnStat == AVGSTDDEV)
                                    MrBayesPrintf (fpMcmc, "\tAvgStdDev(s)(run%d_run%d)", i+1,j+1);
                                else
                                    MrBayesPrintf (fpMcmc, "\tMaxStdDev(s)(run%d_run%d)", i+1, j+1);
                                }
                            else
                                {
                                if (chainParams.diagnStat == AVGSTDDEV)
                                    MrBayesPrintf (fpMcmc, "\tAvgStdDev(s.tree%d)(run%d_run%d)", n+1, i+1, j+1);
                                else
                                    MrBayesPrintf (fpMcmc, "\tMaxStdDev(s.tree%d)(run%d_run%d)", n+1, i+1, j+1);
                                }
                            }
                        }
                    }
                }
            }
        MrBayesPrintf (fpMcmc, "\n");
        fflush (fpMcmc);
        return (NO_ERROR);
        }

#if defined (MPI_ENABLED)
    /* Reassemble info if MPI version */
    if (ReassembleMoveInfo() == ERROR)
        return (ERROR);
    if (chainParams.numChains > 1 && ReassembleSwapInfo() == ERROR)
        return (ERROR);
    if (proc_id != 0)
        return (NO_ERROR);
#endif

    MrBayesPrintf (fpMcmc, "%d", curGen);

    for (n=0; n<chainParams.numRuns; n++)
        {
        if (chainParams.allChains == YES)
            {
            for (j=n*chainParams.numChains; j<(n+1)*chainParams.numChains; j++)
                {
                for (i=0; i<numUsedMoves; i++)
                    {
                    theMove = usedMoves[i];
                    if (theMove->nBatches[j] < 1)
                        MrBayesPrintf (fpMcmc, "\tNA");
                    else
                        MrBayesPrintf (fpMcmc, "\t%.6f", theMove->lastAcceptanceRate[j]);
                    if (theMove->moveType->Autotune != NULL && chainParams.autotune == YES)
                        MrBayesPrintf (fpMcmc, "\t%.6e", theMove->tuningParam[j][0]);
                    }
                }
            }
        else
            {
            j = n*chainParams.numChains;
            for (i=0; i<numUsedMoves; i++)
                {
                theMove = usedMoves[i];
                if (theMove->nBatches[j] < 1)
                    MrBayesPrintf (fpMcmc, "\tNA");
                else
                    MrBayesPrintf (fpMcmc, "\t%.6f", theMove->lastAcceptanceRate[j]);
                if (theMove->moveType->Autotune != NULL && chainParams.autotune == YES)
                    MrBayesPrintf (fpMcmc, "\t%.6e", theMove->tuningParam[j][0]);
                }
            }
        if (chainParams.numChains > 1)
            {
            for (i=0; i<chainParams.numChains; i++)
                {
                for (j=i+1; j<chainParams.numChains; j++)
                    {
                    MrBayesPrintf (fpMcmc, "\t%.6f", (MrBFlt) swapInfo[n][i][j] / (MrBFlt) swapInfo[n][j][i]);
                    }
                }
            }
        }

    if (chainParams.numRuns > 1)
        {
        for (n=0; n<numTopologies; n++)
            {
            if (chainParams.relativeBurnin == NO && curGen < chainParams.chainBurnIn * chainParams.sampleFreq)
                MrBayesPrintf (fpMcmc, "\tNA");
            else
                {
                if (chainParams.diagnStat == AVGSTDDEV)
                    MrBayesPrintf (fpMcmc, "\t%.6f", chainParams.stat[n].avgStdDev);
                else
                    MrBayesPrintf (fpMcmc, "\t%.6f", chainParams.stat[n].max);
                }
            if (chainParams.allComps == YES)
                {
                for (i=0; i<chainParams.numRuns; i++)
                    {
                    for (j=i+1; j<chainParams.numRuns; j++)
                        {
                        if (chainParams.relativeBurnin == NO && curGen < chainParams.chainBurnIn * chainParams.sampleFreq)
                            MrBayesPrintf (fpMcmc, "\tNA");
                        else if (chainParams.diagnStat == AVGSTDDEV)
                            MrBayesPrintf (fpMcmc, "\t%.6f", chainParams.stat[n].pair[i][j] / chainParams.stat[n].pair[j][i]);
                        else /*if (chainParams.diagnStat == MAXSTDDEV) */
                            MrBayesPrintf (fpMcmc, "\t%.6f", chainParams.stat[n].pair[i][j]);
                        }
                    }
                }
            }
        }

    MrBayesPrintf (fpMcmc, "\n");
    fflush (fpMcmc);

#if defined MPI_ENABLED
    /* Redistribute the move info in the parallel version, so that
       swapping occurs correctly; only necessary on processor 0 */
    RedistributeMoveInfo();
#endif

    return (NO_ERROR);
}


/*-----------------------------------------------------------------------
|
|   PrintMPISlaves: Print strings from MPI slave nodes
|
------------------------------------------------------------------------*/
#if defined (MPI_ENABLED)
int PrintMPISlaves (FILE *fp)
{
    char        *s=NULL;
    int         i, len, ierror, nErrors, sumErrors, tag;
    MPI_Status  status;

    nErrors = sumErrors = tag = 0;
    if (proc_id==0)
        {
        s = (char *) SafeCalloc (100, sizeof(char));
        // if (s!=NULL)
        //    lenS = 100;
        // else
        //    lenS = 0;
        }

    for (i=1; i<num_procs; i++)
        {
        /* communicate length */
        if (proc_id == 0)
            {
            /* receive size */
            ierror = MPI_Recv (&len, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            if (ierror != MPI_SUCCESS)
                {
                MrBayesPrint ("%s   Problem receiving string length from proc_id = %d\n", spacer, i);
                nErrors++;
                }
            }
        else if (proc_id == i)
            {
            /* send size */
            len = (int)strlen(printString);
            ierror = MPI_Send (&len, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
            if (ierror != MPI_SUCCESS)
                {
                MrBayesPrint ("%s   Problem sending string length from proc_id = %d\n", spacer, i);
                nErrors++;
                }
            }
        /* reallocate string s on processor 0 if necessary */
        if (nErrors == 0 && proc_id == 0 && len+5 > strlen(s))
            {
            if ((s = (char *) SafeRealloc ((void *)s, ((size_t)len+5)*sizeof(char))) == NULL)
                {
                MrBayesPrint ("%s   Problem reallocating %d chars to string 's' on proc 0 in PrintMPISlaves()\n", spacer, len+5);
                nErrors++;
                }
            }
        /* communicate and print string */
        if (nErrors == 0)
            {
            if (proc_id == 0)
                {
                /* receive string */
                ierror = MPI_Recv (s, len+1, MPI_CHAR, i, tag, MPI_COMM_WORLD, &status);
                if (ierror != MPI_SUCCESS)
                    {
                    MrBayesPrint ("%s   Problem receiving printString from proc_id = %d\n", spacer, i);
                    nErrors++;
                    }
                /* print string */
                if (nErrors == 0)
                    MrBayesPrintf (fp, "%s", s);
                }
            else if (proc_id == i)
                {
                /* send string */
                ierror = MPI_Send (printString, len+1, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
                if (ierror != MPI_SUCCESS)
                    {
                    MrBayesPrint ("%s   Problem sending printString from proc_id = %d\n", spacer, i);
                    nErrors++;
                    }
                }
            }
        if (nErrors > 0)
            break;
        }

    if (proc_id == 0)
        {
        free (s);
        s = NULL;
        }

    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   Problem with the printing in PrintMPISlaves().\n", spacer);
        return (ERROR);
        }

    return (NO_ERROR);
}
#endif


/*----------------------------------------------------------------------
|
|   PrintParamValues: print parameter values and subvalues for param
|
----------------------------------------------------------------------*/
void PrintParamValues (Param *p, int chain, char *s)
{
    int         j;
    MrBFlt      *value0, *value1;
    
    if (p == NULL)
        MrBayesPrint ("%s   %s = NULL\n", spacer, s);
    else
        {
        if (p->nValues > 0)
            {
            value0 = GetParamVals (p, chain, 0);
            value1 = GetParamVals (p, chain, 1);
            for (j=0; j<p->nValues; j++)
                MrBayesPrint ("%s   hyper [%s] = (%lf %lf)\n", spacer, s, value0[j], value1[j]);
            }
        if (p->nSubValues > 0)
            {
            value0 = GetParamSubVals (p, chain, 0);
            value1 = GetParamSubVals (p, chain, 1);
            for (j=0; j<p->nSubValues; j++)
                MrBayesPrint ("%s   %s = (%lf %lf)\n", spacer, s, value0[j], value1[j]);
            }
        }
    MrBayesPrint ("\n\n");

    return;
}


/*----------------------------------------------------------------------
|
|   PrintParsMatrix: Print parsimony (bitset) matrix
|       using hexadecimal representation
|
|
------------------------------------------------------------------------*/
int PrintParsMatrix (void)
{
    int             i, j=0, k, c, d, printWidth, nextColumn, nChars, inputChar;
    BitsLong        x, y, bitsLongOne;
    char            ch;
    ModelInfo       *m;

    bitsLongOne = 1;

    printWidth = 79;

    for (d=0; d<numCurrentDivisions; d++)
        {
        MrBayesPrint ("\nParsimony (bitset) matrix for division %d\n\n", d+1);

        m = &modelSettings[d];

        nChars = 1 + (int) (log((bitsLongOne << m->numStates) - 1) / log(16));
    
        for (c=0; c<m->numChars; c++)
            {
            MrBayesPrint ("Parsimony sets for character %d -- \n", (c / m->nParsIntsPerSite));
            for (i=0; i<numTaxa; i++)
                {
                MrBayesPrint ("%-10.10s   ", taxaNames[i]);
                j = c*m->nParsIntsPerSite;
                for (nextColumn=13; nextColumn < printWidth; nextColumn+=nChars + 1)
                    {
                    if (j >= m->numChars*m->nParsIntsPerSite)
                        break;
                    x = m->parsSets[i][j];
                    for (k=8 - nChars; k<8; k++)
                        {
                        y = (x >> (4* (7 - k))) & 15;
                        // if (y > 16) printf ("y is too big %ld\n",y);
                        if (y < 10)
                            ch = (char) y + '0';
                        else
                            ch = (char) y - 10 + 'A';
                        MrBayesPrint("%c", ch);
                        }
                    MrBayesPrint(" ");
                    j++;
                    }
                MrBayesPrint ("\n");
                }
            MrBayesPrint("\n");
            printf ("Do you want to stop (y/n)?\n");
            inputChar = getchar();
            if (inputChar == 'y' || inputChar == 'Y')
                return NO_ERROR;
            else
                MrBayesPrint ("\n");
            }
        }   /* next division */

    return NO_ERROR;
}


/*------------------------------------------------------------------
|
|   PrintSiteRates_Gen: general n-state models with rate variation
|
-------------------------------------------------------------------*/
int PrintSiteRates_Gen (TreeNode *p, int division, int chain)
{
    int             c, j, k, nStates, hasPInvar;
    MrBFlt          freq, siteLike, invLike, catLike, pInvar=0.0, *bs,
                    *catRate, baseRate;
    MrBFlt          s01, s10, probOn, probOff, *swr, covBF[40];
    CLFlt           *lnScaler, *clP, *siteRates, *clInvar=NULL;
    char            *tempStr;
    int             tempStrSize = TEMPSTRSIZE;
    ModelInfo       *m;
    
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    /* find model settings and nStates, pInvar, invar cond likes */
    m = &modelSettings[division];
    nStates = m->numModelStates;
    if (m->pInvar == NULL)
        {
        hasPInvar = NO;
        }
    else
        {
        hasPInvar = YES;
        pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
        clInvar = m->invCondLikes;
        }

    /* find conditional likelihood pointer */
    clP = m->condLikes[m->condLikeIndex[chain][p->index]];
    
    /* use scratch space for root node for temporary calculations */
    siteRates = m->condLikes[m->condLikeScratchIndex[p->index]];
    
    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);

    /* if covarion model, adjust base frequencies */
    if (m->switchRates != NULL)
        {
        /* find the stationary frequencies */
        swr = GetParamVals(m->switchRates, chain, state[chain]);
        s01 = swr[0];
        s10 = swr[1];
        probOn = s01 / (s01 + s10);
        probOff = 1.0 - probOn;

        /* now adjust the base frequencies; on-state stored first in cond likes */
        for (j=0; j<nStates/2; j++)
            {
            covBF[j] = bs[j] * probOn;
            covBF[j+nStates/2] = bs[j] * probOff;
            }

        /* finally set bs pointer to adjusted values */
        bs = covBF;
        }

    /* find category frequencies */
    if (hasPInvar == NO)
        freq = 1.0 /  m->numRateCats;
    else
        freq = (1.0 - pInvar) /  m->numRateCats;

    /* get rate multipliers (for gamma & partition specific rates) */
    baseRate = GetRate (division, chain);
    
    /* compensate for invariable sites */
    if (hasPInvar == YES)
        baseRate /= (1.0 - pInvar);
        
    /* get category rates */
    catRate = GetParamSubVals (m->shape, chain, state[chain]);

    /* loop over characters */
    if (hasPInvar == NO)
        {
        /* no invariable category */
        for (c=0; c<m->numChars; c++)
            {
            siteLike = 0.0;
            siteRates[c] = 0.0;
            for (k=0; k<m->numRateCats; k++)
                {
                catLike = 0.0;
                for (j=0; j<nStates; j++)
                    catLike += (*(clP++)) * bs[j];
                siteRates[c] += (CLFlt) (catLike * catRate[k]);
                siteLike += catLike;
                }
            siteRates[c] *= (CLFlt) (baseRate / siteLike);  /* category frequencies and site scaler cancel out */
            }
        }
    else
        {
        /* has invariable category */
        for (c=0; c<m->numChars; c++)
            {
            siteLike = invLike = 0.0;
            siteRates[c] = 0.0;
            for (k=0; k<m->numRateCats; k++)
                {
                catLike = 0.0;
                for (j=0; j<nStates; j++)
                    catLike += (*(clP++)) * bs[j];
                siteRates[c] += (CLFlt) (catLike * catRate[k]);
                siteLike += catLike;
                }
            siteLike *= freq;
            siteRates[c] *= (CLFlt) freq;
            for (j=0; j<nStates; j++)
                invLike += (*(clInvar++)) * bs[j];
            siteLike += (invLike /  exp (lnScaler[c]) * pInvar);
            /* we do not need to add the invariable category into siteRates before rescaling because the rate is 0.0 */
            siteRates[c] *= (CLFlt) (baseRate / siteLike);  /* site scaler cancels out; category frequencies dealt with above */
            }
        }
        
    /* print the resulting site rates cycling over uncompressed chars */
    for (c=0; c<numChar; c++)
        {
        if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum] != division+1)
            continue;
        j = compCharPos[c] - m->compCharStart;
        SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(siteRates[j]));
        if (AddToPrintString (tempStr) == ERROR) return (ERROR);
        }

    free (tempStr);
    return NO_ERROR;
}


/*------------------------------------------------------------------
|
|   PrintSiteRates_Std: standard model with rate variation
|
-------------------------------------------------------------------*/
int PrintSiteRates_Std (TreeNode *p, int division, int chain)
{
    int             c, j, k, nStates;
    MrBFlt          siteLike, catLike, *bs, *catRate, baseRate;
    CLFlt           *clP, *siteRates;
    char            *tempStr;
    int             tempStrSize = TEMPSTRSIZE;
    ModelInfo       *m;
    
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    /* find model settings */
    m = &modelSettings[division];

    /* find conditional likelihood pointer */
    clP = m->condLikes[m->condLikeIndex[chain][p->index]];
    
    /* use scratch space for root node for temporary calculations */
    siteRates = m->condLikes[m->condLikeScratchIndex[p->index]];
    
    /* find base frequencies */
    bs = GetParamStdStateFreqs (m->stateFreq, chain, state[chain]);

    /* get rate multiplier */
    baseRate = GetRate (division, chain);
    
    /* get category rates */
    catRate = GetParamSubVals (m->shape, chain, state[chain]);

    /* loop over characters */
    for (c=0; c<m->numChars; c++)
        {
        siteLike = 0.0;
        siteRates[c] = 0.0;
        nStates = m->nStates[c];
        for (k=0; k<m->numRateCats; k++)
            {
            catLike = 0.0;
            for (j=0; j<nStates; j++)
                catLike += (*(clP++)) * bs[j];
            siteRates[c] += (CLFlt) (catLike * catRate[k]);
            siteLike += catLike;
            }
        siteRates[c] *= (CLFlt)(baseRate / siteLike);   /* category frequencies and site scaler cancel out */
        }
        
    /* print the resulting site rates cycling over uncompressed chars */
    for (c=0; c<numChar; c++)
        {
        if (charInfo[c].isExcluded == YES || partitionId[c][partitionNum] != division+1)
            continue;
        j = compCharPos[c] - m->compCharStart;
        SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(siteRates[j]));
        if (AddToPrintString (tempStr) == ERROR) return (ERROR);
        }

    free (tempStr);
    return NO_ERROR;
}


int PrintStates (int curGen, int coldId)
{
    int             d, i, j, k, k1, compressedCharPosition, *printedChar=NULL, origAlignmentChars[3];
    char            *partString=NULL, stateString[4];
    MrBFlt          *st, *sst, sum;
    Param           *p;
    ModelInfo       *m;
    Tree            *tree;
    TreeNode        *node;
    ModelParams     *mp;
    char            *tempStr;
    int             tempStrSize;

    /* allocate the print string */
    printStringSize = tempStrSize = TEMPSTRSIZE;
    printString = (char *)SafeMalloc((size_t)printStringSize * sizeof(char));
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));

    if (!printString)
        {
        MrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
        goto errorExit;
        }
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, (size_t)(tempStrSize) * sizeof(char));
        goto errorExit;
        }
    
    *printString = '\0';
    *tempStr = '\0';

    /* Allocate memory, temporarily, in case we infer positive selection probs, site omegas, or site rates */
    if (inferPosSel == YES || inferSiteOmegas == YES)
        {
        if (memAllocs[ALLOC_POSSELPROBS] == YES)
            {
            MrBayesPrint ("%s   posSelProbs not free in PrintStates\n", spacer);
            goto errorExit;
            }
        posSelProbs = (MrBFlt *)SafeMalloc((size_t)numCompressedChars * sizeof(MrBFlt));
        if (!posSelProbs)
            {
            MrBayesPrint ("%s   Problem allocating posSelProbs (%d)\n", spacer, numCompressedChars * sizeof(MrBFlt));
            goto errorExit;
            }
        for (i=0; i<numCompressedChars; i++)
            posSelProbs[i] =  -10.0;
        memAllocs[ALLOC_POSSELPROBS] = YES;
        }
    if (inferPosSel == YES || inferSiteOmegas == YES || inferSiteRates == YES || inferAncStates == YES)
        {
        printedChar = (int *)SafeMalloc((size_t)numChar * sizeof(int));
        if (!printedChar)
            {
            MrBayesPrint ("%s   Problem allocating printedChar (%d)\n", spacer, numChar * sizeof(int));
            goto errorExit;
            }
        for (i=0; i<numChar; i++)
            printedChar[i] = NO;
        }

    /* Set up the header to the file. */
    if (curGen == 0)
        {
        SafeSprintf (&tempStr, &tempStrSize, "[ID: %s]\n", stamp);
        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
        SafeSprintf (&tempStr, &tempStrSize, "Gen");
        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
        SafeSprintf (&tempStr, &tempStrSize, "\tLnL");
        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
        SafeSprintf (&tempStr, &tempStrSize, "\tLnPr");
        if (AddToPrintString (tempStr) == ERROR) goto errorExit;

        for (i=0; i<numParams; i++)
            {
            p = &params[i];
            if (p->paramType == P_BRLENS)
                {
                /* print tree lengths or heights for all trees */
                tree = GetTree (p, coldId, state[coldId]);
                if (tree->isRooted == YES)
                    {
                    if (FillRelPartsString(p, &partString) == YES)
                        SafeSprintf (&tempStr, &tempStrSize, "\tTH%s\tTL%s", partString, partString);
                    else
                        SafeSprintf (&tempStr, &tempStrSize, "\tTH\tTL");
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                else
                    {
                    if (FillRelPartsString(p, &partString) == YES)
                        SafeSprintf (&tempStr, &tempStrSize, "\tTL%s", partString);
                    else
                        SafeSprintf (&tempStr, &tempStrSize, "\tTL");
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                /* print # cpp events, or relaxed clock model indicator */
                for (j=0; j<p->nSubParams; j++)
                    {
                    if (p->subParams[j]->paramType == P_CPPEVENTS)
                        {
                        if (FillRelPartsString(p->subParams[j], &partString) == YES)
                            SafeSprintf (&tempStr, &tempStrSize, "\tn_CPP%s", partString);
                        else
                            SafeSprintf (&tempStr, &tempStrSize, "\tn_CPP");
                        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                        }
                    else if (p->subParams[j]->paramType == P_MIXEDBRCHRATES)
                        {
                        if (FillRelPartsString(p->subParams[j], &partString) == YES)
                            SafeSprintf (&tempStr, &tempStrSize, "\tm_RCl%s", partString);
                        else
                            SafeSprintf (&tempStr, &tempStrSize, "\tm_RCl");
                        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                        }
                    }
                }
            /* print proportion of ancestral fossils */
            else if (p->paramType == P_FOSLRATE)
                {
                if (FillRelPartsString(p, &partString) == YES)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\tprop_ancfossil%s", partString);
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                else
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\tprop_ancfossil");
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            }

        /* print ordinary parameters header */
        for (i=0; i<numPrintParams; i++)
            {
            p = printParam[i];
            SafeSprintf (&tempStr, &tempStrSize, "\t%s", p->paramHeader);
            if (AddToPrintString (tempStr) == ERROR) goto errorExit;
            }
            
        /* print substitution model parameters header */
        if (inferSiteRates == YES)
            {
            for (i=0; i<numChar; i++)
                printedChar[i] = NO;
    
            for (i=0; i<numChar; i++)
                { 
                if (charInfo[i].isExcluded == YES)
                    continue;
                if (printedChar[i] == YES)
                    continue;
                d = partitionId[i][partitionNum] - 1;
                m = &modelSettings[d];
                mp = &modelParams[d];
                if (m->printSiteRates == YES)
                    {
                    if (m->nCharsPerSite == 1)
                        {
                        SafeSprintf (&tempStr, &tempStrSize, "\tr(%d)", i+1);
                        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                        }
                    else
                        {
                        origAlignmentChars[0] = i;
                        k = 1;
                        for (j=i+1; j<numChar; j++)
                            {
                            if (compCharPos[i] == compCharPos[j])
                                {
                                if (k > m->nCharsPerSite)
                                    return (ERROR);
                                origAlignmentChars[k++] = j;
                                printedChar[j] = YES;
                                }
                            }
                        if (k != m->nCharsPerSite)
                            return (ERROR);
                        SafeSprintf (&tempStr, &tempStrSize, "\tr(%d,", origAlignmentChars[0]);
                        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                        for (j=1; j<k-1; j++)
                            {
                            SafeSprintf (&tempStr, &tempStrSize, "%d,", origAlignmentChars[j]);
                            if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                            }
                        SafeSprintf (&tempStr, &tempStrSize, "%d)", origAlignmentChars[k-1]);
                        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                        }
                    }
                }
            }

        if (inferPosSel == YES)
            {
            for (i=0; i<numChar; i++)
                printedChar[i] = NO;
            for (d=0; d<numCurrentDivisions; d++)
                {
                m = &modelSettings[d];
                tree = GetTree(m->brlens, coldId, state[coldId]);
                if (m->printPosSel == YES)
                    {
                    if (m->PosSelProbs (tree->root->left, d, coldId) == ERROR)
                        goto errorExit;
                    }
                }
            /* for (i=0; i<numChar; i++)
                printf ("%4d -- %3d %3d\n", i, compCharPos[i], compColPos[i]); */
            for (i=0; i<numChar; i++)
                {
                compressedCharPosition = compCharPos[i];
                if (!(posSelProbs[compressedCharPosition] < 0.0) && printedChar[i] == NO && charInfo[i].isExcluded == NO)
                    {
                    for (j=k=0; j<numChar; j++)
                        {
                        if (charInfo[j].charId == charInfo[i].charId)
                            {
                            origAlignmentChars[k++] = j;
                            printedChar[j] = YES;
                            }
                        }
                    SafeSprintf (&tempStr, &tempStrSize, "\tpr+(%d,%d,%d)", origAlignmentChars[0]+1, origAlignmentChars[1]+1, origAlignmentChars[2]+1);
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }   
            for (i=0; i<numChar; i++)
                printedChar[i] = NO;
            }
            
        if (inferSiteOmegas == YES)
            {
            for (i=0; i<numChar; i++)
                printedChar[i] = NO;
            for (d=0; d<numCurrentDivisions; d++)
                {
                m = &modelSettings[d];
                tree = GetTree(m->brlens, coldId, state[coldId]);
                if (m->printSiteOmegas == YES)
                    {
                    if (m->SiteOmegas (tree->root->left, d, coldId) == ERROR)
                        goto errorExit;
                    }
                }
            /* for (i=0; i<numChar; i++)
                printf ("%4d -- %3d %3d\n", i, compCharPos[i], compColPos[i]); */
            for (i=0; i<numChar; i++)
                {
                compressedCharPosition = compCharPos[i];
                if (posSelProbs[compressedCharPosition] >= 0.0 && printedChar[i] == NO && charInfo[i].isExcluded == NO)
                    {
                    for (j=k=0; j<numChar; j++)
                        {
                        if (charInfo[j].charId == charInfo[i].charId)
                            {
                            origAlignmentChars[k++] = j;
                            printedChar[j] = YES;
                            }
                        }
                    SafeSprintf (&tempStr, &tempStrSize, "\tomega(%d,%d,%d)", origAlignmentChars[0]+1, origAlignmentChars[1]+1, origAlignmentChars[2]+1);
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }   
            for (i=0; i<numChar; i++)
                printedChar[i] = NO;
            }
            
        if (inferAncStates == YES)
            {
            for (j=0; j<numChar; j++)
                printedChar[j] = NO;
            for (d=0; d<numCurrentDivisions; d++)
                {
                m = &modelSettings[d];
                mp = &modelParams[d];
                if (m->printAncStates != YES)
                    continue;
                for (i=0; i<numDefinedConstraints; i++)
                    {
                    if (mp->activeConstraints[i] == NO || definedConstraintsType[i] != HARD)
                        continue;
                    for (j=0; j<numChar; j++)
                        {
                        if (partitionId[j][partitionNum] - 1 != d || charInfo[j].isExcluded == YES || printedChar[j] == YES)
                            continue;
                        if (mp->dataType == STANDARD)
                            {
                            for (k=0; k<m->nStates[compCharPos[j] - m->compCharStart]; k++)
                                {
                                SafeSprintf (&tempStr, &tempStrSize, "\tp(%c){%d@%s}", m->StateCode(k), j+1, constraintNames[i]);
                                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                                }
                            }
                        else if ((mp->dataType == DNA || mp->dataType == RNA) && !strcmp(mp->nucModel,"Codon"))
                            {
                            origAlignmentChars[0] = j+1;
                            k1 = 1;
                            for (k=j+1; k<numChar; k++)
                                {
                                if (charInfo[k].charId == charInfo[j].charId)
                                    {
                                    origAlignmentChars[k1++] = k+1;
                                    printedChar[k] = YES;
                                    }
                                }
                            for (k=0; k<m->numStates; k++)
                                {
                                State_CODON(stateString, k, d);
                                SafeSprintf (&tempStr, &tempStrSize, "\tp(%s){%d,%d,%d@%s}",
                                    stateString,
                                    origAlignmentChars[0],
                                    origAlignmentChars[1],
                                    origAlignmentChars[2],
                                    constraintNames[i]);
                                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                                }
                            }
                        else if ((mp->dataType == DNA || mp->dataType == RNA) && !strcmp(mp->nucModel,"Doublet"))
                            {
                            origAlignmentChars[0] = j+1;
                            k1 = 1;
                            for (k=j+1; k<numChar; k++)
                                {
                                if (charInfo[k].charId == charInfo[j].charId)
                                    {
                                    origAlignmentChars[k1++] = k+1;
                                    printedChar[k] = YES;
                                    }
                                }
                            for (k=0; k<m->numStates; k++)
                                {
                                State_DOUBLET(stateString, k);
                                SafeSprintf (&tempStr, &tempStrSize, "\tp(%s){%d,%d@%s}",
                                    stateString,
                                    origAlignmentChars[0],
                                    origAlignmentChars[1],
                                    constraintNames[i]);
                                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                                }
                            }
                        else if ((mp->dataType == DNA || mp->dataType == RNA) && !strcmp(mp->nucModel,"Protein"))
                            {
                            origAlignmentChars[0] = j+1;
                            k1 = 1;
                            for (k=j+1; k<numChar; k++)
                                {
                                if (charInfo[k].charId == charInfo[j].charId)
                                    {
                                    origAlignmentChars[k1++] = k+1;
                                    printedChar[k] = YES;
                                    }
                                }
                            for (k=0; k<m->numStates; k++)
                                {
                                SafeSprintf (&tempStr, &tempStrSize, "\tp(%c){%d,%d,%d@%s}",
                                    m->StateCode(k),
                                    origAlignmentChars[0],
                                    origAlignmentChars[1],
                                    origAlignmentChars[2],
                                    constraintNames[i]);
                                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                                }
                            }
                        else
                            {
                            for (k=0; k<m->numStates; k++)
                                {
                                SafeSprintf (&tempStr, &tempStrSize, "\tp(%c){%d@%s}", m->StateCode(k), j+1, constraintNames[i]);
                                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                                }
                            }
                        }
                    }
                }
            }
            
        SafeSprintf (&tempStr, &tempStrSize, "\n");
        if (AddToPrintString (tempStr) == ERROR) goto errorExit;
        }
        
    /* now print parameter values */
    SafeSprintf (&tempStr, &tempStrSize, "%d", curGen);
    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(curLnL[coldId]));
    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(curLnPr[coldId]));
    if (AddToPrintString (tempStr) == ERROR) goto errorExit;

    for (i=0; i<numParams; i++)
        {
        p = &params[i];
        if (p->paramType == P_BRLENS)
            {
            /* print tree lengths or heights for all trees */
            tree = GetTree (p, coldId, state[coldId]);
            if (tree->isRooted == NO)
                {
                SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(TreeLength(p, coldId)));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                }
            else
                {
                SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(tree->root->left->nodeDepth));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(TreeLength(p, coldId)));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                }
            /* print # cpp events, or relaxed clock model indicator */
            for (j=0; j<p->nSubParams; j++)
                {
                if (p->subParams[j]->paramType == P_CPPEVENTS)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%d", NumCppEvents(p->subParams[j],coldId));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                else if (p->subParams[j]->paramType == P_MIXEDBRCHRATES)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%d", *GetParamIntVals(p->subParams[j],coldId,state[coldId]));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            }
        else if (p->paramType == P_FOSLRATE)
            {
            /* print proportion of ancestral fossils */
            SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(PropAncFossil(p, coldId)));
            if (AddToPrintString (tempStr) == ERROR) goto errorExit;
            }
        }

    /* print ordinary parameters */
    for (i=0; i<numPrintParams; i++)
        {
        p = printParam[i];

        /* get model params and settings */
        mp = &modelParams[p->relParts[0]];
        m  = &modelSettings[p->relParts[0]];
        
        st  = GetParamVals (p, coldId, state[coldId]);
        sst = GetParamSubVals (p, coldId, state[coldId]);

        if (p->paramId == SYMPI_EXP_MS || p->paramId == SYMPI_UNI_MS || p->paramId == SYMPI_FIX_MS)
            {
            /* We print symmetric dirichlet alpha value if not fixed and then multistate character state frequencies */
            if (p->paramId != SYMPI_FIX_MS)
                {
                SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(st[0]));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                }
            sst = GetParamStdStateFreqs (p, coldId, state[coldId]);
            if (p->hasBinaryStd == YES)
                sst += 2 * m->numBetaCats;
            for (j=0; j<p->nSympi; j++)
                {
                for (k=0; k<p->sympinStates[j]; k++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(*sst++));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            }
        else if (p->paramType == P_PI || p->paramType == P_MIXTURE_RATES)
            {
            /* We print the subvalues if we are dealing with state frequencies (state frequencies are held in subvalues) OR
               if we are dealing with rates of a siterate mixture (rates held in subvalues). */
            for (j=0; j<p->nSubValues; j++)
                {
                SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(sst[j]));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                }
            }
        else if (p->paramType == P_TRATIO && !strcmp(mp->tratioFormat,"Dirichlet"))
            {
            SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(st[0] / (1.0 + st[0])));
            if (AddToPrintString (tempStr) == ERROR) goto errorExit;
            SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(1.0 / (1.0 + st[0])));
            if (AddToPrintString (tempStr) == ERROR) goto errorExit;
            }
        else if (p->paramType == P_REVMAT)
            {
            if (!strcmp(mp->revmatFormat,"Ratio"))
                {
                sum = st[p->nValues-1];
                for (j=0; j<p->nValues; j++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(st[j] / sum));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            else
                {
                /* we already have rate proportions */
                for (j=0; j<p->nValues; j++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(st[j]));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            if (p->paramId == REVMAT_MIX)
                {
                /* add model index and k for nst=mixed */
                SafeSprintf (&tempStr, &tempStrSize, "\t%d", FromGrowthFxnToIndex(GetParamIntVals(p, coldId, state[coldId])));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                SafeSprintf (&tempStr, &tempStrSize, "\t%d", GetKFromGrowthFxn(GetParamIntVals(p, coldId, state[coldId])));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                }
            }
        else if (p->paramType == P_RATEMULT)
            {
            if (!strcmp(mp->ratemultFormat,"Ratio"))
                {
                for (j=0; j<p->nValues; j++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(st[j] / st[0]));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            else if (!strcmp(mp->ratemultFormat, "Dirichlet"))
                {
                sum = 0.0;
                for (j=0; j<p->nValues; j++)
                    sum += st[j];
                for (j=0; j<p->nValues; j++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(st[j] / sum));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            else
                {
                for (j=0; j<p->nValues; j++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(st[j]));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            }
        else if (p->paramType == P_AAMODEL)
            {
            for (j=0; j<p->nValues; j++)
                {
                SafeSprintf (&tempStr, &tempStrSize, "\t%d", (int)st[j]);
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                }
            }
        else
            {
            /* run of the mill parameter */
            for (j=0; j<p->nValues; j++)
                {
                SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(st[j]));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                }
            }

        if (p->paramType == P_OMEGA && p->paramId != OMEGA_DIR && p->paramId != OMEGA_FIX && p->paramId != OMEGA_FFF && p->paramId != OMEGA_FF && p->paramId != OMEGA_10FFF)
            {
            /* OK, we also need to print subvalues for the category frequencies in a NY98-like model. */
            if (!strcmp(mp->omegaVar, "M10"))
                {
                for (j=0; j<4; j++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(sst[mp->numM10BetaCats + mp->numM10GammaCats + 4 + j]));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                for (j=0; j<2; j++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(sst[mp->numM10BetaCats + mp->numM10GammaCats + j]));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            else
                {
                for (j=0; j<3; j++)
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(sst[j]));
                    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                    }
                }
            }
        }
        
    /* if user wants site rates, we print those here */
    if (inferSiteRates == YES)
        {
        for (d=0; d<numCurrentDivisions; d++)
            {
            m = &modelSettings[d];
            if (m->printSiteRates == YES)
                {
                mp = &modelParams[d];
                tree = GetTree (m->brlens, coldId, state[coldId]);
                node = tree->root->left;
                m->PrintSiteRates (node, d, coldId);
                }
            }
        }

    /* If the user wants to infer sites that are under positive selection, then we need to print out the posterior
       probability that each site is a positively selected one here. */
    if (inferPosSel == YES)
        {
        /* loop over the divisions, calculating the probability of being in the positively
           selected class for each relevant partition */
        for (d=0; d<numCurrentDivisions; d++)
            {
            m = &modelSettings[d];
            tree = GetTree(m->brlens, coldId, state[coldId]);
            if (m->printPosSel == YES)
                {
                if (m->PosSelProbs (tree->root->left, d, coldId) == ERROR)
                    {
                    goto errorExit;
                    }
                }
            }

        /* print the probabilities for the appropriate sites in the original alignment */
        for (i=0; i<numChar; i++)
            printedChar[i] = NO;
        for (i=0; i<numChar; i++)
            {
            compressedCharPosition = compCharPos[i];
            if (!(posSelProbs[compressedCharPosition] < 0.0) && printedChar[i] == NO && charInfo[i].isExcluded == NO)
                {
                for (j=k=0; j<numChar; j++)
                    {
                    if (charInfo[j].charId == charInfo[i].charId)
                        {
                        origAlignmentChars[k++] = j;
                        printedChar[j] = YES;
                        }
                    }
                SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(posSelProbs[compressedCharPosition]));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                /* printf ("%4d -> (%3d,%3d,%3d) %1.25le\n", i, origAlignmentChars[0]+1, origAlignmentChars[1]+1,
                                                                origAlignmentChars[2]+1, posSelProbs[compressedCharPosition]); */
                }
            }
        }
    
    /* If the user wants omega values for each site, we print those here. */
    if (inferSiteOmegas == YES)
        {
        /* loop over the divisions, calculating the omega value for each site */
        for (d=0; d<numCurrentDivisions; d++)
            {
            m = &modelSettings[d];
            tree = GetTree(m->brlens, coldId, state[coldId]);
            if (m->printSiteOmegas == YES)
                {
                if (m->SiteOmegas (tree->root->left, d, coldId) == ERROR)
                    {
                    goto errorExit;
                    }
                }
            }

        /* print the site omegas for the appropriate sites in the original alignment */
        /* note that we use posSelProbs to pass values between SiteOmegas and this function */
        for (i=0; i<numChar; i++)
            printedChar[i] = NO;
        for (i=0; i<numChar; i++)
            {
            compressedCharPosition = compCharPos[i];
            if (posSelProbs[compressedCharPosition] >= 0.0 && printedChar[i] == NO && charInfo[i].isExcluded == NO)
                {
                for (j=k=0; j<numChar; j++)
                    {
                    if (charInfo[j].charId == charInfo[i].charId)
                        {
                        origAlignmentChars[k++] = j;
                        printedChar[j] = YES;
                        }
                    }
                SafeSprintf (&tempStr, &tempStrSize, "\t%s", MbPrintNum(posSelProbs[compressedCharPosition]));
                if (AddToPrintString (tempStr) == ERROR) goto errorExit;
                /* printf ("%4d -> (%3d,%3d,%3d) %1.25le\n", i, origAlignmentChars[0]+1, origAlignmentChars[1]+1,
                                                                origAlignmentChars[2]+1, posSelProbs[compressedCharPosition]); */
                }
            }
        }
     
    /* free memory for positive selection probs or site omegas */
    if (inferPosSel == YES || inferSiteOmegas == YES)
        {
        if (memAllocs[ALLOC_POSSELPROBS] == YES)
            free (posSelProbs);
        memAllocs[ALLOC_POSSELPROBS] = NO;
        free (printedChar);
        }

    /* if user wants ancestral states for constrained nodes, we obtain and print those here */
    if (inferAncStates == YES)
        {
        for (d=0; d<numCurrentDivisions; d++)
            {
            m = &modelSettings[d];
            if (m->printAncStates == YES)
                {
                mp = &modelParams[d];
                tree = GetTree (m->brlens, coldId, state[coldId]);
                for (i=j=tree->nIntNodes - 1; i>=0; i--)
                    {
                    node = tree->intDownPass[i];
                    m->CondLikeUp (node, d, coldId);
                    }
                for (k=0; k<numDefinedConstraints; k++)
                    {
                    if (mp->activeConstraints[k] == NO || definedConstraintsType[k] != HARD)
                        continue;
                    for (i=tree->nIntNodes-1; i>=0; i--)
                        {
                        node = tree->intDownPass[i];
                        if (node->isLocked == YES && k == node->lockID)
                            m->PrintAncStates (node, d, coldId);
                        }
                    }
                }
            }
        }

    SafeSprintf (&tempStr, &tempStrSize, "\n");
    if (AddToPrintString (tempStr) == ERROR) goto errorExit;
    
    free (tempStr);
    SafeFree ((void *)&partString);
    
    return (NO_ERROR);
    
    errorExit:
        if (printedChar)
            free (printedChar);
        if (memAllocs[ALLOC_POSSELPROBS] == YES)
            free (posSelProbs);
        memAllocs[ALLOC_POSSELPROBS] = NO;
        free (tempStr);
        SafeFree ((void *)&partString);
        return (ERROR);
}


/*----------------------------------------------------------------------
|
|   PrintStatesToFiles: Print trees and model parameters to files. We
|      only come into this function if it is the first cycle of the chain
|      or if we hit a cycle number evenly divisible by the sample frequency,
|      or this is the last cycle of the chain.
|
------------------------------------------------------------------------*/
int PrintStatesToFiles (int curGen)
{
    int             i, j, chn, coldId, runId;
    MrBFlt          clockRate;
    Tree            *tree=NULL;
    Param           *param;
#   if defined (MPI_ENABLED)
    int             id, x, doesThisProcHaveId, procWithChain, ierror, tag, nErrors, sumErrors;
    MPI_Status      status;
#   endif

#   if !defined (MPI_ENABLED)

    /* print parameter values and trees (single-processor version) */
    for (chn=0; chn<numLocalChains; chn++)
        {
        if ((chainId[chn] % chainParams.numChains) == 0)
            {
            coldId = chn;
            runId = chainId[chn] / chainParams.numChains;

            /* print parameter values */
            if (PrintStates (curGen, coldId) == ERROR)
                return (ERROR);
            fprintf (fpParm[runId], "%s", printString);
            fflush (fpParm[runId]);
            free(printString);

            /* print trees */
            for (i=0; i<numPrintTreeParams; i++)
                {
                param = printTreeParam[i];
                tree = GetTree(param, coldId, state[coldId]);
                if (param->paramType == P_TOPOLOGY)
                    {
                    if (tree->isClock == YES)
                        clockRate = *GetParamVals(modelSettings[tree->relParts[0]].clockRate, coldId, state[coldId]);
                    else
                        clockRate = 0.0;
                    if (PrintTree (curGen, param, coldId, NO, clockRate) == ERROR)
                        return (ERROR);
                    }
                else
                    {
                    if (tree->isClock == YES)
                        clockRate = *GetParamVals(modelSettings[tree->relParts[0]].clockRate, coldId, state[coldId]);
                    else
                        clockRate = 0.0;
                    if (PrintTree (curGen, param, coldId, YES, clockRate) == ERROR)
                        return (ERROR);
                    }

                fprintf (fpTree[runId][i], "%s", printString);
                fflush (fpTree[runId][i]);
                free(printString);

                j = printTreeTopologyIndex[i];
                if (j<numTopologies)
                    {
                    if (chainParams.mcmcDiagn == YES && chainParams.numRuns > 1)
                        {
                        if (chainParams.relativeBurnin == YES || curGen >= chainParams.chainBurnIn * chainParams.sampleFreq)
                            {
                            if (AddTreeToPartitionCounters (tree, j, runId) == ERROR)
                                return ERROR;
                            if (chainParams.relativeBurnin == YES && chainParams.saveTrees == YES && (noWarn == NO || curGen <= chainParams.stopTreeGen))
                                {
                                ResetTopologyFromTree (chainParams.dtree, tree);
                                if (AddToTreeList (&chainParams.treeList[numTopologies*runId+j], chainParams.dtree) == ERROR)
                                    return (ERROR);
                                }
                            }
                        }
                    }
                }
            }
        }
#   else
    /* print parameter values and trees (parallel version) */
    
    /* Wait for all of the processors to get to this point before starting the printing. */
    ierror = MPI_Barrier (MPI_COMM_WORLD);
    if (ierror != MPI_SUCCESS)
        {
        MrBayesPrint ("%s   Problem at chain barrier.\n", spacer);
        return ERROR;
        }
    tag = nErrors = 0;
        
    /* Loop over runs. */
    for (runId=0; runId<chainParams.numRuns; runId++)
        {
        /* Get the ID of the chain we want to print. Remember, the ID's should be numbered
           0, 1, 2, ..., numChains X numRuns. Chains numbered 0, numChains, 2 X numChains, ...
           are cold. */
        id = runId * chainParams.numChains;
        
        /* Does this processor have the chain? */
        doesThisProcHaveId = NO;
        coldId = 0;
        for (chn=0; chn<numLocalChains; chn++)
            {
            if (chainId[chn] == id)
                {
                doesThisProcHaveId = YES;
                coldId = chn;
                break;
                }
            }
        
        /* Tell all the processors which has the chain we want to print. We do this using the MPI_AllReduce
           function. If the processor does not have the chain, then it initializes x = 0. If it does
           have the chain, then x = proc_id. When the value of x is summed over all the processors, the sum
           should be the proc_id of the processor with the chain. Possible values are 0, 1, 2, num_procs-1. 
           Note that every processor knows procWithChain because we are using MPI_Allreduce, instead of MPI_Reduce. */
        x = 0;
        if (doesThisProcHaveId == YES)
            x = proc_id;
        ierror = MPI_Allreduce (&x, &procWithChain, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (ierror != MPI_SUCCESS)
            {
            MrBayesPrint ("%s   Problem finding processor with chain to print.\n", spacer);
            return (ERROR);
            }

        /* ****************************************************************************************************/
        /* print parameter values *****************************************************************************/
        
        /* Fill printString with the contents to be printed on proc_id = 0. Note
           that printString is allocated in the function. */
        if (doesThisProcHaveId == YES)
            {
            if (PrintStates (curGen, coldId) == ERROR)
                nErrors++;
            }
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Problem with PrintStates.\n", spacer);
            return ERROR;
            }
        
        /* First communication: Send/receive the length of the printString. */
        if (proc_id == 0 || proc_id == procWithChain)
            {
            if (procWithChain != 0)
                {
                if (proc_id == procWithChain)
                    {
                    /* Find out how large the string is, and send the information to proc_id = 0. */
                    ierror = MPI_Send (&printStringSize, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                    if (ierror != MPI_SUCCESS)
                        nErrors++;
                    }
                else
                    {
                    /* Receive the length of the string from proc_id = procWithChain, and then allocate
                       printString to be that length. */
                    ierror = MPI_Recv (&printStringSize, 1, MPI_LONG, procWithChain, tag, MPI_COMM_WORLD, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        MrBayesPrint ("%s   Problem receiving printStringSize from proc_id = %d\n", spacer, procWithChain);
                        nErrors++;
                        }
                    printString = (char *)SafeMalloc((size_t)printStringSize * sizeof(char));
                    if (!printString)
                        {
                        MrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
                        nErrors++;
                        }
                    strcpy (printString, "");
                    }
                }
            }
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Problem with first communication (states).\n", spacer);
            return ERROR;
            }

        /* Second communication: Send/receive the printString. */
        if (proc_id == 0 || proc_id == procWithChain)
            {
            if (procWithChain != 0)
                {                   
                if (proc_id == procWithChain)
                    {
                    /* Send the printString to proc_id = 0. After we send the string to proc_id = 0, we can
                       free it. */
                    ierror = MPI_Send (&printString[0], (int)printStringSize, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
                    if (ierror != MPI_SUCCESS)
                        nErrors++;
                    free(printString);
                    }
                else
                    {
                    /* Receive the printString from proc_id = procWithChain. */
                    ierror = MPI_Recv (&printString[0], (int)printStringSize, MPI_CHAR, procWithChain, tag, MPI_COMM_WORLD, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        MrBayesPrint ("%s   Problem receiving printString from proc_id = %d\n", spacer, procWithChain);
                        nErrors++;
                        }
                    }
                }
            }
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Problem with second communication (states).\n", spacer);
            return ERROR;
            }

        /* Print the string with the parameter information if we are proc_id = 0. */
        if (proc_id == 0)
            {
            fprintf (fpParm[runId], "%s", printString);
            fflush (fpParm[runId]);
            free(printString);
            }

        /* ****************************************************************************************************/
        /* print trees ****************************************************************************************/

        for (i=0; i<numPrintTreeParams; i++)
            {
            /* Print trees to file. */

            /* Fill printString with the contents to be printed on proc_id = 0. Note
               that printString is allocated in the function. */
            if (doesThisProcHaveId == YES)
                {
                param = printTreeParam[i];
                tree = GetTree(param, coldId, state[coldId]);
                if (param->paramType == P_TOPOLOGY)
                    {
                    if (tree->isClock == YES)
                        clockRate = *GetParamVals(modelSettings[tree->relParts[0]].clockRate, coldId, state[coldId]);
                    else
                        clockRate = 0.0;
                    if (PrintTree (curGen, param, coldId, NO, clockRate) == ERROR)
                        nErrors++;
                    }
                else
                    {
                    if (tree->isClock == YES)
                        clockRate = *GetParamVals(modelSettings[tree->relParts[0]].clockRate, coldId, state[coldId]);
                    else
                        clockRate = 0.0;
                    if (PrintTree (curGen, param, coldId, YES, clockRate) == ERROR)
                        nErrors++;
                    }
                }
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Problem with printing trees.\n", spacer);
                return ERROR;
                }
                
            /* First communication: Send/receive the length of the printString. */
            if (proc_id == 0 || proc_id == procWithChain)
                {
                if (procWithChain != 0)
                    {
                    if (proc_id == procWithChain)
                        {
                        /* Find out how large the string is, and send the information to proc_id = 0. */
                        ierror = MPI_Send (&printStringSize, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                        if (ierror != MPI_SUCCESS)
                            nErrors++;
                        }
                    else
                        {
                        /* Receive the length of the string from proc_id = procWithChain, and then allocate
                           printString to be that length. */
                        ierror = MPI_Recv (&printStringSize, 1, MPI_LONG, procWithChain, tag, MPI_COMM_WORLD, &status);
                        if (ierror != MPI_SUCCESS)
                            {
                            MrBayesPrint ("%s   Problem receiving printStringSize from proc_id = %d\n", spacer, procWithChain);
                            nErrors++;
                            }
                        printString = (char *)SafeMalloc((size_t)printStringSize * sizeof(char));
                        if (!printString)
                            {
                            MrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
                            nErrors++;
                            }
                        strcpy (printString, "");
                        }
                    }
                }
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Problem with first communication (states).\n", spacer);
                return ERROR;
                }

            /* Second communication: Send/receive the printString. */
            if (proc_id == 0 || proc_id == procWithChain)
                {
                if (procWithChain != 0)
                    {                   
                    if (proc_id == procWithChain)
                        {
                        /* Send the printString to proc_id = 0. After we send the string to proc_id = 0, we can
                           free it. */
                        ierror = MPI_Send (&printString[0], (int)printStringSize, MPI_CHAR, 0, tag, MPI_COMM_WORLD);
                        if (ierror != MPI_SUCCESS)
                            nErrors++;
                        free(printString);
                        }
                    else
                        {
                        /* Receive the printString from proc_id = procWithChain. */
                        ierror = MPI_Recv (&printString[0], (int)printStringSize, MPI_CHAR, procWithChain, tag, MPI_COMM_WORLD, &status);
                        if (ierror != MPI_SUCCESS)
                            {
                            MrBayesPrint ("%s   Problem receiving printString from proc_id = %d\n", spacer, procWithChain);
                            nErrors++;
                            }
                        }
                    }
                }
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Problem with second communication (states).\n", spacer);
                return ERROR;
                }

            /* Print the string with the parameter information if we are proc_id = 0. */
            if (proc_id == 0)
                {
                fprintf (fpTree[runId][i], "%s", printString);
                fflush (fpTree[runId][i]);
                j = printTreeTopologyIndex[i];
                if (j < numTopologies)
                    {
                    if (chainParams.numRuns > 1 && chainParams.mcmcDiagn == YES)
                        {
                        if (chainParams.relativeBurnin == YES || curGen >= chainParams.chainBurnIn * chainParams.sampleFreq)
                            {
                            char *s = NULL;
                            StripComments (printString);
                            /* if it is the first tree, we strip out the translate block first (twice) */
                            if (curGen==0) {
                                if (strtok (printString, ";")==NULL) /* get translate lock */
                                    return (ERROR);
                                if (strtok (NULL, ";")==NULL)
                                    return (ERROR);
                                if (strtok (NULL, "\n\t\r ")==NULL) /* get 'tree' */
                                    return (ERROR);
                                }
                            else {
                                if (strtok (printString, "\n\t\r ")==NULL) /* get 'tree' */
                                    return (ERROR);
                                }
                            if (strtok (NULL, " =")==NULL)  /* get 'rep.xxxx' */
                                return (ERROR);
                            if ((s = strtok (NULL, " =;"))==NULL)  /* get Newick string */
                                return (ERROR);
                            ResetTopology (chainParams.dtree, s);
                            if (AddTreeToPartitionCounters (chainParams.dtree, j, runId) == ERROR)
                                return ERROR;
                            if (chainParams.relativeBurnin == YES && chainParams.saveTrees == YES && (noWarn == NO || curGen <= chainParams.stopTreeGen))
                                {
                                if (AddToTreeList (&chainParams.treeList[runId*numTopologies+j], chainParams.dtree) == ERROR)
                                    return (ERROR);
                                }
                            }
                        }
                    }
                free(printString);
                }
            }

        /* Have all of the chains wait here, until the string has been successfully printed on proc_id = 0. */
        ierror = MPI_Barrier (MPI_COMM_WORLD);
        if (ierror != MPI_SUCCESS)
            {
            MrBayesPrint ("%s   Problem at chain barrier.\n", spacer);
            return ERROR;
            }
        }
#   endif
        
    return (NO_ERROR);
}


int PrintSwapInfo (void)
{
    int         i, j, n, maxNumExchanges, len, maxLen, reweightingChars=0;
    char        *tempStr;
    int             tempStrSize;

    if (chainParams.numChains == 1)
        return NO_ERROR;

#   if defined (MPI_ENABLED)
    if (ReassembleSwapInfo() == ERROR)
        return ERROR;
    if (proc_id != 0)
        return NO_ERROR;
#   endif

    tempStrSize = TEMPSTRSIZE;
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    for (n=0; n<chainParams.numRuns; n++)
        {
        maxNumExchanges = 0;
        for (i=0; i<chainParams.numChains; i++)
            for (j=0; j<chainParams.numChains; j++)
                if (i > j && swapInfo[n][i][j] > maxNumExchanges)
                    maxNumExchanges = swapInfo[n][i][j];
        SafeSprintf (&tempStr, &tempStrSize, "%d", maxNumExchanges);
        maxLen = (int) strlen(tempStr);
        if (maxLen < 4)
            maxLen = 4;
            
        reweightingChars = NO;
        if ((chainParams.weightScheme[0] + chainParams.weightScheme[1]) > 0.00001)
            reweightingChars = YES;

        if (chainParams.numRuns == 1)
            MrBayesPrint ("\n%s   Chain swap information:\n\n", spacer);
        else
            MrBayesPrint ("\n%s   Chain swap information for run %d:\n\n", spacer, n+1);

        MrBayesPrint ("%s          ", spacer);
        for (j=0; j<chainParams.numChains; j++)
            {
            SafeSprintf (&tempStr, &tempStrSize, "%d", j+1);
            len = (int) strlen(tempStr);
            MrBayesPrint ("%*c %d ", maxLen-len, ' ', j+1);
            }
        MrBayesPrint ("\n");
        
        MrBayesPrint ("%s        --", spacer);
        for (j=0; j<chainParams.numChains; j++)
            {
            MrBayesPrint ("--");
            for (i=0; i<maxLen; i++)
                MrBayesPrint ("-");
            }
        MrBayesPrint ("\n");
        
        for (i=0; i<chainParams.numChains; i++)
            {
            MrBayesPrint ("%s   %4d | ", spacer, i+1);
            for (j=0; j<chainParams.numChains; j++)
                {
                if (i < j)
                    {
                    if (swapInfo[n][j][i] <= 0)
                        {
                        MrBayesPrint ("%*c%s ", maxLen-3, ' ', " NA ");
                        }
                    else
                        {
                        SafeSprintf (&tempStr, &tempStrSize, "%1.2lf", (MrBFlt)swapInfo[n][i][j]/swapInfo[n][j][i]);
                        len = (int) strlen(tempStr);
                        MrBayesPrint ("%*c%1.2lf ", maxLen-len+1, ' ', (MrBFlt)swapInfo[n][i][j]/swapInfo[n][j][i]);
                        }
                    }
                else if (i == j)
                    {
                    MrBayesPrint ("%*c ", maxLen+1, ' ');
                    }
                else
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "%d", swapInfo[n][i][j]);
                    len = (int) strlen(tempStr);
                    MrBayesPrint ("%*c%d ", maxLen-len+1, ' ', swapInfo[n][i][j]);
                    }
                }
            MrBayesPrint ("\n");
            }
        }

    MrBayesPrint ("\n%s   Upper diagonal: Proportion of successful state exchanges between chains\n", spacer);
    MrBayesPrint ("%s   Lower diagonal: Number of attempted state exchanges between chains\n", spacer);
        
    MrBayesPrint ("\n%s   Chain information:\n\n", spacer);
    MrBayesPrint ("%s     ID -- Heat ", spacer);
    if (reweightingChars == YES)
        MrBayesPrint ("%% Dn %% Up\n");
    else
        MrBayesPrint ("\n");
    
    MrBayesPrint ("%s    -----------", spacer);
    if (reweightingChars == YES)
        MrBayesPrint ("----------\n");
    else
        MrBayesPrint ("\n");
    for (i=0; i<chainParams.numChains; i++)
        {
        MrBayesPrint ("%s   %4d -- %1.2lf ", spacer, i+1, Temperature (i)/*1.0 / (1.0 + chainParams.chainTemp * i)*/);
        if (reweightingChars == YES)
            {
            if (i == 0)
                {
                MrBayesPrint ("  0%%   0%% (cold chain)\n");
                }
            else
                {
                SafeSprintf (&tempStr, &tempStrSize, "%d", (int)chainParams.weightScheme[0]);
                len = (int) strlen(tempStr);
                MrBayesPrint ("%*c%d%% ", 3-len, ' ', (int)chainParams.weightScheme[0]);
                SafeSprintf (&tempStr, &tempStrSize, "%d", (int)chainParams.weightScheme[1]);
                len = (int) strlen(tempStr);
                MrBayesPrint ("%*c%d%% \n", 3-len, ' ', (int)chainParams.weightScheme[1]);
                }
            }
        else
            {
            if (i == 0)
                MrBayesPrint (" (cold chain)\n");
            else
                MrBayesPrint ("\n");
            }
        }
    if (chainParams.userDefinedTemps == NO)
        {
        MrBayesPrint ("\n%s   Heat = 1 / (1 + T * (ID - 1))\n", spacer);
        MrBayesPrint ("%s      (where T = %1.2lf is the temperature and ID is the chain number)\n", spacer, chainParams.chainTemp);
        }
    if (reweightingChars == YES)
        MrBayesPrint ("%s   Reweighting increment = %1.2lf\n", spacer, chainParams.weightScheme[2]);
    MrBayesPrint ("\n");
        
    free (tempStr);
    return (NO_ERROR);
}


/*----------------------------------------------------------------------
|
|   PrintTermState: Print terminal state index matrix
|
------------------------------------------------------------------------*/
int PrintTermState (void)
{
    int             i, j=0, c, d, printWidth, nextColumn, nDigits, nReps;
    ModelInfo       *m;
    ModelParams     *mp;

    printWidth = 79;

    for (d=0; d<numCurrentDivisions; d++)
        {
        MrBayesPrint ("\nTerminal state index matrix for division %d\n\n", d+1);

        m = &modelSettings[d];
        mp = &modelParams[d];

        if (!strcmp(mp->covarionModel, "Yes"))
            nReps = 2;
        else
            nReps = 1;

        nDigits = 1 + (int)(log10(mp->nStates * mp->nStates * nReps));
    
        for (c=m->compCharStart; c<m->compCharStop; c=j)
            {
            for (i=0; i<numTaxa; i++)
                {
                MrBayesPrint ("%-10.10s   ", taxaNames[i]);
                j = c;
                for (nextColumn=13; nextColumn < printWidth; nextColumn+=nDigits + 1)
                    {
                    if (j >= m->compCharStop)
                        break;
                    MrBayesPrint ("%*d ",nDigits, m->termState[i][j-c]);
                    j++;
                    }
                MrBayesPrint ("\n");
                }
            MrBayesPrint("\n");
            }
        }   /* next division */

    return NO_ERROR;
}


/*--------------------------------------------------
|
|   PrintTiProbs: This function is for debugging of
|       tiProbs; it will print a square matrix of
|       tiProbs, check row sums, and check for time
|       reversibility
|
---------------------------------------------------*/
void PrintTiProbs (CLFlt *tP, MrBFlt *bs, int nStates)
{
    int     i, j;
    CLFlt   *tiP, sum;

    tiP = tP;

    printf ("\nTransition matrix\n");
    for (i=0; i<nStates; i++)
        {
        printf ("\t%d", i);
        }
    printf ("\tsum\n");

    for (i=0; i<nStates; i++)
        {
        printf ("%d\t", i);
        sum = 0.0;
        for (j=0; j<nStates; j++)
            {
            printf ("\t%.6f",tP[j]);
            sum += tP[j];
            }
        printf ("\t%.6f\n",sum);
        tP += nStates;
        }

    printf ("\nStationary state frequencies\n");
    for (i=0; i<nStates; i++)
        printf ("%d -- %f\n",i,bs[i]);
    
    printf ("\nTime reversibility\n");

    printf ("State 1\tState 2\tforward\tbackward\tabs diff\n");
    for (i=0; i<nStates; i++)
        {
        for (j=i+1; j<nStates; j++)
            {
            printf ("%d\t%d\t%.6f\t%.6f\t%.6f\n", i, j, tiP[i*nStates+j]*bs[i],
                tiP[j*nStates+i]*bs[j], fabs(tiP[i*nStates+j]*bs[i] - tiP[j*nStates+i]*bs[j]));
            }
        }

    getchar();
    return;
}


int PrintTopConvInfo (void)
{
    int         i, j, n, len, maxLen;
    char        *tempStr;
    int         tempStrSize;
    MrBFlt      maxNumPartitions;
    STATS       *stat;

    if (chainParams.numRuns == 1)
        return NO_ERROR;

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return (NO_ERROR);
#   endif

    tempStrSize = TEMPSTRSIZE;
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    for (n=0; n<numTopologies; n++)
        {
        stat = &(chainParams.stat[n]);
        maxNumPartitions = 0.0;
        for (i=0; i<chainParams.numRuns; i++)
            for (j=0; j<chainParams.numRuns; j++)
                if (i > j && stat->pair[i][j] > maxNumPartitions)
                    maxNumPartitions = stat->pair[i][j];
        SafeSprintf (&tempStr, &tempStrSize, "%d", (int) maxNumPartitions);
        maxLen = (int) strlen(tempStr);
        if (maxLen < 5)
            maxLen = 5;
        
        if (numTopologies == 1)
            {
            if (chainParams.diagnStat == AVGSTDDEV)
                MrBayesPrint ("%s   Pairwise average standard deviation of split frequencies (upper triangle)\n", spacer);
            else
                MrBayesPrint ("%s   Pairwise maximum standard deviation of split frequencies (upper triangle)\n", spacer);
            MrBayesPrint ("%s      and number of qualifying splits for each comparison (lower triangle):\n\n", spacer);
            }
        else
            {
            if (chainParams.diagnStat == AVGSTDDEV)
                MrBayesPrint ("%s   Pairwise average standard deviation of split frequencies in topology %d (upper triangle)\n", spacer, n);
            else
                MrBayesPrint ("%s   Pairwise maximum standard deviation of split frequencies in topology %d (upper triangle)\n", spacer, n);
            MrBayesPrint ("%s      and number of qualifying splits for each comparison (lower triangle):\n\n", spacer);
            }

        MrBayesPrint ("%s          ", spacer);
        for (j=0; j<chainParams.numRuns; j++)
            {
            SafeSprintf (&tempStr, &tempStrSize, "%d", j+1);
            len = (int) strlen(tempStr);
            MrBayesPrint ("%*c %d ", maxLen-len, ' ', j+1);
            }
        MrBayesPrint ("\n");
    
        MrBayesPrint ("%s        --", spacer);
        for (j=0; j<chainParams.numRuns; j++)
            {
            MrBayesPrint ("--");
            for (i=0; i<maxLen; i++)
                MrBayesPrint ("-");
            }
        MrBayesPrint ("\n");
    
        for (i=0; i<chainParams.numRuns; i++)
            {
            MrBayesPrint ("%s   %4d | ", spacer, i+1);
            for (j=0; j<chainParams.numRuns; j++)
                {
                if (i < j)
                    {
                    if (chainParams.diagnStat == AVGSTDDEV)
                        SafeSprintf (&tempStr, &tempStrSize, "%1.3lf", (MrBFlt) (stat->pair[i][j]) / (MrBFlt) (stat->pair[j][i]));
                    else /* if (chainParams.diagnStat == MAXSTDDEV) */
                        SafeSprintf (&tempStr, &tempStrSize, "%1.3lf", stat->pair[i][j]);
                    len = (int) strlen(tempStr);
                    MrBayesPrint ("%*c%s ", maxLen-len+1, ' ', tempStr);
                    }
                else if (i == j)
                    {
                    MrBayesPrint ("%*c ", maxLen+1, ' ');
                    }
                else
                    {
                    SafeSprintf (&tempStr, &tempStrSize, "%d", (int) stat->pair[i][j]);
                    len = (int) strlen(tempStr);
                    MrBayesPrint ("%*c%s ", maxLen-len+1, ' ', tempStr);
                    }
                }
            MrBayesPrint ("\n");
            }
    
        MrBayesPrint ("\n");
        }

    free (tempStr);
    return (NO_ERROR);
}


void PrintToScreen (int curGen, int startGen, time_t endingT, time_t startingT)
{
    int         i, chn, nHours, nMins, nSecs;
    MrBFlt      timePerGen;

#   if defined (MPI_ENABLED)
    int         numLocalColdChains, numFirstAndLastCold;
    
    if (curGen == 0)
        {
        if (chainParams.isSS == NO && chainParams.mcmcDiagn == YES && chainParams.numRuns > 1)
            {
            MrBayesPrint ("\n");
            if (chainParams.relativeBurnin == YES)
                MrBayesPrint ("%s   Using a relative burnin of %.1f %% for diagnostics\n", spacer, 100.0*chainParams.burninFraction);
            else
                MrBayesPrint ("%s   Using an absolute burnin of %d samples for diagnostics\n", spacer, chainParams.chainBurnIn);
            }
        MrBayesPrint ("\n");
        MrBayesPrint ("%s   Chain results (%d generations requested):\n\n", spacer, chainParams.numGen);
        }
    MrBayesPrint ("%s   %4d -- ", spacer, curGen);
    numLocalColdChains = numFirstAndLastCold = 0;
    for (chn=0; chn<numLocalChains; chn++)
        {
        if ((chainId[chn] % chainParams.numChains) == 0)
            {
            numLocalColdChains++;
            if (chn == 0 || chn == numLocalChains - 1)
                numFirstAndLastCold++;
            }
        }

    i = 1;
    for (chn=0; chn<numLocalChains; chn++)
        {
        if (i > chainParams.printMax)   
            {
            if (i == chainParams.printMax +1)
                {
                i++;
                if (numLocalColdChains > 0 && numLocalColdChains > numFirstAndLastCold)
                    MrBayesPrint ("[...%d more local chains...] ", numLocalChains - chainParams.printMax);
                else
                    MrBayesPrint ("(...%d more local chains...) ", numLocalChains - chainParams.printMax);
                continue;
                }
            else
                continue;
            }
        if ((chainId[chn] % chainParams.numChains) == 0)
            {
            i++;
            if (chainParams.printAll == YES)
                MrBayesPrint ("[%1.3lf] ", curLnL[chn]);
            else
                MrBayesPrint ("[%1.3lf] .. ", curLnL[chn]);
            }
        else if (chainParams.printAll == YES)
            {
            i++;
            MrBayesPrint ("(%1.3lf) ", curLnL[chn]);
            }
        if (chn < numLocalChains - 1 && (chainId[chn] / chainParams.numChains != chainId[chn+1] / chainParams.numChains))
            MrBayesPrint ("* ");
        }
        
    if (numLocalColdChains == chainParams.numRuns)
        MrBayesPrint ("(...%d remote chains...) ", (chainParams.numChains*chainParams.numRuns) - numLocalChains);
    else
        MrBayesPrint ("[...%d remote chains...] ", (chainParams.numChains*chainParams.numRuns) - numLocalChains);

    if (curGen > 0)
        {
        timePerGen = (MrBFlt) ((MrBFlt)(endingT-startingT)/(MrBFlt)(curGen-startGen));
        nSecs = (int)((chainParams.numGen - curGen) * timePerGen);
        nHours = nSecs / 3600;
        nSecs  = nSecs % 3600;
        nMins  = nSecs / 60; 
        nSecs  = nSecs % 60;
        MrBayesPrint ("-- %d:%0.2d:%0.2d", nHours, nMins, nSecs);
        }
    MrBayesPrint ("\n");
    fflush (stdout);
    
#   else

    if (curGen == 0)
        {
        if (chainParams.isSS == NO && chainParams.mcmcDiagn == YES && chainParams.numRuns > 1)
            {
            MrBayesPrint ("\n");
            if (chainParams.relativeBurnin == YES)
                MrBayesPrint ("%s   Using a relative burnin of %.1f %% for diagnostics\n", spacer, 100.0*chainParams.burninFraction);
            else
                MrBayesPrint ("%s   Using an absolute burnin of %d samples for diagnostics\n", spacer, chainParams.chainBurnIn);
            }
        MrBayesPrint ("\n");
        MrBayesPrint ("%s   Chain results (%d generations requested):\n\n", spacer, chainParams.numGen);
        }
    MrBayesPrint ("%s   %5d -- ", spacer, curGen);
    if (numLocalChains == 1)
        MrBayesPrint ("%1.3lf ", curLnL[0]);
    else
        {
        i = 0;
        for (chn=0; chn<numLocalChains; chn++)
            {
            if (i >= chainParams.printMax)
                {
                if (i == chainParams.printMax)
                    MrBayesPrint (".. ");
                i++;
                continue;
                }
            if (chainParams.numChains == 1)
                {
                MrBayesPrint ("%1.3lf ", curLnL[chn]);
                i++;
                }
            else if (chainId[chn] % chainParams.numChains == 0)
                {
                if (chainParams.printAll == YES)
                    MrBayesPrint ("[%1.3lf] ", curLnL[chn]);
                else
                    MrBayesPrint ("[%1.3lf][%d] .. ", curLnL[chn], chn % chainParams.numChains + 1);
                i++;
                }
            else if (chainParams.printAll == YES)
                {
                MrBayesPrint ("(%1.3lf) ", curLnL[chn]);
                i++;
                }
            if (chn < numLocalChains - 1 && (chainId[chn] / chainParams.numChains != chainId[chn+1] / chainParams.numChains)
                && i < chainParams.printMax - 1)
                MrBayesPrint ("* ");
            }
        }
        
    if (curGen > 0)
        {
        timePerGen = (MrBFlt) ((MrBFlt)(endingT-startingT)/(MrBFlt)(curGen-startGen));
        nSecs = (int)((chainParams.numGen - curGen) * timePerGen);
        nHours = nSecs / 3600;
        nSecs  = nSecs % 3600;
        nMins  = nSecs / 60; 
        nSecs  = nSecs % 60;
        MrBayesPrint ("-- %d:%0.2d:%0.2d", nHours, nMins, nSecs);
        }
    MrBayesPrint ("\n");
    
    fflush (stdout);
    
#   endif
        
}


int PrintTree (int curGen, Param *treeParam, int chain, int showBrlens, MrBFlt clockRate)
{
    int             i, tempStrSize;
    char            *tempStr;
    Tree            *tree;
    TreeNode        *p=NULL, *q;
    Param           *subParm;

    /* allocate the print string */
    printStringSize = 200;
    printString = (char *)SafeMalloc((size_t)printStringSize * sizeof(char));
    if (!printString)
        {
        MrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
        return (ERROR);
        }
    *printString = '\0';

    tempStrSize = TEMPSTRSIZE;
    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    /* get tree */
    tree = GetTree(treeParam, chain, state[chain]);

    /* order the taxa */
    if (chainParams.orderTaxa == YES)
        {
        for (i=0; i<tree->nNodes-1; i++)
            {
            p = tree->allDownPass[i];
            if (p->left == NULL)
                {
                if (p->index == localOutGroup)
                    p->x = -1;
                else
                    p->x = p->index;
                }
            else if (p->left->x < p->right->x)
                p->x = p->left->x;
            else
                p->x = p->right->x;
            }
        for (i=0; i<tree->nIntNodes; i++)
            {
            if (p->left->x > p->right->x)
                {
                q = p->left;
                p->left = p->right;
                p->right = q;
                }
            }
        }
    
    /* print the translate block information and the top of the file */
    if (curGen == 0)
        {
        /* print #NEXUS and translation block information */
        SafeSprintf (&tempStr, &tempStrSize, "#NEXUS\n");
        if (AddToPrintString (tempStr) == ERROR) return(ERROR);
        SafeSprintf (&tempStr, &tempStrSize, "[ID: %s]\n", stamp);
        if (AddToPrintString (tempStr) == ERROR) return(ERROR);
        SafeSprintf (&tempStr, &tempStrSize, "[Param: tree");
        if (AddToPrintString (tempStr) == ERROR) return(ERROR);
        if (numCurrentDivisions == 1)
            {
            SafeSprintf (&tempStr, &tempStrSize, "]\n");
            if (AddToPrintString (tempStr) == ERROR) return(ERROR);
            }
        else if (numCurrentDivisions == tree->nRelParts)
            {
            SafeSprintf (&tempStr, &tempStrSize, "{all}]\n");
            if (AddToPrintString (tempStr) == ERROR) return(ERROR);
            }
        else
            {
            SafeSprintf (&tempStr, &tempStrSize, "{%d", tree->relParts[0]+1);
            if (AddToPrintString (tempStr) == ERROR) return(ERROR);
            for (i=1; i<tree->nRelParts; i++)
                {
                SafeSprintf (&tempStr, &tempStrSize, ",%d", tree->relParts[i]+1);
                if (AddToPrintString (tempStr) == ERROR) return(ERROR);
                }
            SafeSprintf (&tempStr, &tempStrSize, "}]\n");
            if (AddToPrintString (tempStr) == ERROR) return(ERROR);
            }
        SafeSprintf (&tempStr, &tempStrSize, "begin trees;\n");
        if (AddToPrintString (tempStr) == ERROR) return(ERROR);
        SafeSprintf (&tempStr, &tempStrSize, "   translate\n");
        if (AddToPrintString (tempStr) == ERROR) return(ERROR);
        if (treeParam->paramType == P_SPECIESTREE)
            {
            for (i=0; i<numSpecies; i++)
                {
                if (i != numSpecies - 1)
                    SafeSprintf (&tempStr, &tempStrSize, "      %2d %s,\n", i+1, speciesNameSets[speciespartitionNum].names[i]);
                else
                    SafeSprintf (&tempStr, &tempStrSize, "      %2d %s;\n", i+1, speciesNameSets[speciespartitionNum].names[i]);
                if (AddToPrintString (tempStr) == ERROR) return(ERROR);
                }
            }
        else
            {
            for (i=0; i<numLocalTaxa; i++)
                {
                if (i != numLocalTaxa - 1)
                    SafeSprintf (&tempStr, &tempStrSize, "      %2d %s,\n", i+1, localTaxonNames[i]);
                else
                    SafeSprintf (&tempStr, &tempStrSize, "      %2d %s;\n", i+1, localTaxonNames[i]);
                if (AddToPrintString (tempStr) == ERROR) return(ERROR);
                }
            }
        }
    
    /* write the tree preamble */
    if (SafeSprintf (&tempStr, &tempStrSize, "   tree gen.%d", curGen) == ERROR) return (ERROR);
    if (AddToPrintString (tempStr) == ERROR) return(ERROR);
    if (treeParam->paramType == P_BRLENS && treeParam->nSubParams > 0)
        {
        for (i=0; i<treeParam->nSubParams; i++)
            {
            subParm = treeParam->subParams[i];
            if (subParm->paramType == P_CPPEVENTS)
                {
                if (SafeSprintf (&tempStr, &tempStrSize, " [&E %s]", subParm->name) == ERROR) return (ERROR);
                if (AddToPrintString (tempStr) == ERROR) return(ERROR);
                }
            //  if (subParm->paramType == P_MIXEDBRCHRATES)
            //  {
            //  id = *GetParamIntVals(subParm, chain, state[chain]);
            //  if (SafeSprintf (&tempStr, &tempStrSize, " [&B %s %d]", subParm->name, id) == ERROR) return (ERROR);
            //  }
            else
                if (SafeSprintf (&tempStr, &tempStrSize, " [&B %s]", subParm->name) == ERROR) return (ERROR);
            if (AddToPrintString (tempStr) == ERROR) return(ERROR);
            }
        }
    subParm = modelSettings[treeParam->relParts[0]].popSize;
    if (treeParam->paramType == P_SPECIESTREE && subParm->nValues > 1)
        {
        if (SafeSprintf (&tempStr, &tempStrSize, " [&N %s]", subParm->name) == ERROR) return (ERROR);
        if (AddToPrintString (tempStr) == ERROR) return(ERROR);
        }

    /* write the tree in (extended) Newick format */
    if (tree->isRooted == YES && tree->isCalibrated == NO)
        SafeSprintf (&tempStr, &tempStrSize, " = [&R] ");
    else if (tree->isRooted == YES && tree->isCalibrated == YES)
        SafeSprintf (&tempStr, &tempStrSize, " = [&R] [&clockrate=%s] ", MbPrintNum(clockRate));
    else /* if (tree->isRooted == NO) */
        SafeSprintf (&tempStr, &tempStrSize, " = [&U] ");
    if (AddToPrintString (tempStr) == ERROR) return(ERROR);
    WriteNoEvtTreeToPrintString (tree->root->left, chain, treeParam, showBrlens, tree->isRooted);
    SafeSprintf (&tempStr, &tempStrSize, ";\n");
    if (AddToPrintString (tempStr) == ERROR) return(ERROR);

    free (tempStr); 
    return (NO_ERROR);
}


#if defined (MPI_ENABLED)
int ReassembleMoveInfo (void)
{
    int             i, n, ierror;
    double          x[7], sum[7];
    MCMCMove        *mv;

    for (n=0; n<numGlobalChains; n++)
        {
        for (i=0; i<numUsedMoves; i++)
            {
            mv = usedMoves[i];

            /* collect counts */
            x[0] = mv->nAccepted[n];
            x[1] = mv->nTried[n];
            x[2] = mv->nBatches[n];
            x[3] = mv->nTotAccepted[n];
            x[4] = mv->nTotTried[n];
            x[5] = mv->lastAcceptanceRate[n];
            if (mv->moveType->Autotune != NULL)
                x[6]=mv->tuningParam[n][0];

            ierror = MPI_Allreduce (&x, &sum, 7, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (ierror != MPI_SUCCESS)
                return (ERROR);

            if (proc_id == 0)
                {
                mv->nAccepted[n]          = (int)sum[0];
                mv->nTried[n]             = (int)sum[1];
                mv->nBatches[n]           = (int)sum[2];
                mv->nTotAccepted[n]       = (int)sum[3];
                mv->nTotTried[n]          = (int)sum[4];
                mv->lastAcceptanceRate[n] = (MrBFlt)sum[5];
                if (mv->moveType->Autotune != NULL)
                    mv->tuningParam[n][0]=(MrBFlt)sum[6];
                }
            }
        }

    return (NO_ERROR);
}


int ReassembleParamVals (int *curId)
{
    int             i, j, k, orderLen, nBrlens, lower, upper, numChainsForProc, proc, ierror, *y, *order, *id, *nEvents;
    MrBFlt          *x, *brlens, **position, **rateMult;
    MPI_Status      status;
    MPI_Request     request;
    Tree            *tree;
    Param           *p;

    extern MrBFlt   *paramValues;
    extern int      paramValsRowSize;
    extern int      intValsRowSize;

    for (i=0; i<numLocalChains; i++)
        curId[i] = chainId[i];

    numChainsForProc = numGlobalChains / num_procs;
    if (numGlobalChains % num_procs > 0)
        lower = upper = numChainsForProc+1;
    else
        lower = upper = numChainsForProc;

    for (proc=1; proc<num_procs; proc++, lower=upper)
        {
        if (proc < numGlobalChains % num_procs)
            upper += numChainsForProc+1;
        else
            upper += numChainsForProc;
        
        /* chain ids */
        if (proc_id == 0)
            {
        id = curId + lower;
            ierror = MPI_Irecv (id, upper-lower, MPI_INT, proc, 0, MPI_COMM_WORLD, &request);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            ierror = MPI_Waitall (1, &request, &status);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            }
        else if (proc_id == proc)
            {
            id = curId;
            ierror = MPI_Isend (id, upper-lower, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            ierror = MPI_Waitall (1, &request, &status);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            }
        
        /* chain states */
        if (proc_id == 0)
            {
            ierror = MPI_Irecv (state+lower, upper-lower, MPI_CHAR, proc, 0, MPI_COMM_WORLD, &request);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            ierror = MPI_Waitall (1, &request, &status);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            }
        else if (proc_id == proc)
            {
            ierror = MPI_Isend (state, upper-lower, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &request);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            ierror = MPI_Waitall (1, &request, &status);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            }

        /* normal parameter values */
        if (proc_id == 0)
            {
            x = paramValues + 2*paramValsRowSize*lower;
            ierror = MPI_Irecv (x, paramValsRowSize*2*(upper-lower), MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            ierror = MPI_Waitall (1, &request, &status);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            if (intValsRowSize > 0)
                {
                y = intValues + 2*intValsRowSize*lower;
                ierror = MPI_Irecv (y, intValsRowSize*2*(upper-lower), MPI_INT, proc, 0, MPI_COMM_WORLD, &request);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Waitall (1, &request, &status);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                }
            }
        else if (proc_id == proc)
            {
            x = paramValues;
            ierror = MPI_Isend (x, paramValsRowSize*2*(upper-lower), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            ierror = MPI_Waitall (1, &request, &status);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            if (intValsRowSize > 0)
                {
                y = intValues;
                ierror = MPI_Isend (y, intValsRowSize*2*(upper-lower), MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Waitall (1, &request, &status);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                }
            }

        /* std state frequencies */
        if (stdStateFreqsRowSize > 0)
            {
            if (proc_id == 0)
                {
                x = stdStateFreqs + 2*stdStateFreqsRowSize*lower;
                ierror = MPI_Irecv (x, stdStateFreqsRowSize*2*(upper-lower), MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Waitall (1, &request, &status);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                }
            else if (proc_id == proc)
                {
                x = stdStateFreqs;
                ierror = MPI_Isend (x, stdStateFreqsRowSize*2*(upper-lower), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Waitall (1, &request, &status);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                }
            }
        
        /* mcmc trees */
        brlens = (MrBFlt *) SafeCalloc (2*numLocalTaxa, sizeof(MrBFlt));
        order  = (int *)    SafeCalloc (2*numLocalTaxa,   sizeof(int));
        for (i=lower; i<upper; i++)
            {
            for (j=0; j<numTrees; j++)
                {
                tree = GetTreeFromIndex(j,0,0);
                orderLen = 2*tree->nIntNodes - 1;
                nBrlens = tree->nNodes - 1;
                if (proc_id == 0)
                    {
                    tree = GetTreeFromIndex(j,i,state[i]);
                    ierror = MPI_Irecv (order, orderLen, MPI_INT, proc, 0, MPI_COMM_WORLD, &request);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (1, &request, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Irecv (brlens, nBrlens, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (1, &request, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    if (tree->isRooted == YES)
                        RetrieveRTreeWithIndices(tree, order, brlens);
                    else
                        {
                        RetrieveUTree(tree, order, brlens);
                        if (localOutGroup!=0)
                            MoveCalculationRoot(tree,localOutGroup);
                        }
                    /* since we only transferred some info, there are additional things we need to
                       consider, like constraints and calibrations; tree names are OK on proc 0 */
                    InitializeTreeCalibrations(tree);
                    CheckSetConstraints(tree);
                    SetDatedNodeAges(modelSettings[tree->relParts[0]].brlens, i, state[i]);
                    }
                else if (proc_id == proc)
                    {
                    tree = GetTreeFromIndex(j,i-lower,state[i-lower]);
                    if (tree->isRooted == YES)
                        StoreRTreeWithIndices(tree, order, brlens);
                    else
                        StoreUTree(tree, order, brlens);
                    ierror = MPI_Isend (order, orderLen, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (1, &request, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Isend (brlens, nBrlens, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (1, &request, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    }
                }
            }
        free (brlens);
        free (order);

        /* CPP event parameters */
        for (i=lower; i<upper; i++)
            {
            for (j=0; j<numParams; j++)
                {
                p = &params[j];
                if (p->paramType == P_CPPEVENTS)
                    {
                    if (proc_id == proc)
                        {
                        /* get pointers */
                        nEvents = p->nEvents[2*(i-lower)+state[i-lower]];
                        position = p->position[2*(i-lower)+state[i-lower]];
                        rateMult = p->rateMult[2*(i-lower)+state[i-lower]];

                        /* send number of events */
                        ierror = MPI_Isend (nEvents, 2*numLocalTaxa, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
                        if (ierror != MPI_SUCCESS)
                            return (ERROR);
                        ierror = MPI_Waitall (1, &request, &status);
                        if (ierror != MPI_SUCCESS)
                            return (ERROR);

                        /* send events and clear pointers */
                        for (k=0; k<2*numLocalTaxa; k++)
                            {
                            if (nEvents[k] > 0)
                                {
                                ierror = MPI_Isend (position[k], nEvents[k], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                ierror = MPI_Waitall (1, &request, &status);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);

                                ierror = MPI_Isend (rateMult[k], nEvents[k], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                ierror = MPI_Waitall (1, &request, &status);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);

                                free(position[k]);
                                free(rateMult[k]);
                                position[k] = NULL;
                                rateMult[k] = NULL;
                                nEvents[k] = 0;
                                }
                            }
                        }
                    else if (proc_id == 0)
                        {
                        /* find pointers */
                        nEvents = p->nEvents[2*i];
                        position = p->position[2*i];
                        rateMult = p->rateMult[2*i];

                        /* clear previous events */
                        for (k=0; k<2*numLocalTaxa; k++)
                            {
                            if (nEvents[k] > 0)
                                {
                                free (position[k]);
                                free (rateMult[k]);
                                position[k] = NULL;
                                rateMult[k] = NULL;
                                nEvents[k] = 0;
                                }
                            }

                        /* receive events */
                        ierror = MPI_Irecv (nEvents, 2*numLocalTaxa, MPI_INT, proc, 0, MPI_COMM_WORLD, &request);
                        if (ierror != MPI_SUCCESS)
                            return (ERROR);
                        ierror = MPI_Waitall (1, &request, &status);
                        if (ierror != MPI_SUCCESS)
                            return (ERROR);

                        for (k=0; k<2*numLocalTaxa; k++)
                            {
                            if (nEvents[k] > 0)
                                {
                                position[k] = (MrBFlt *) SafeCalloc (nEvents[k], sizeof(MrBFlt));
                                rateMult[k] = (MrBFlt *) SafeCalloc (nEvents[k], sizeof(MrBFlt));

                                ierror = MPI_Irecv (position[k], nEvents[k], MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                ierror = MPI_Waitall (1, &request, &status);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);

                                ierror = MPI_Irecv (rateMult[k], nEvents[k], MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                ierror = MPI_Waitall (1, &request, &status);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                }
                            }
                        }
                    }
                }
            } 
        }

    return (NO_ERROR);
}


int ReassembleSwapInfo (void)
{
    int i, j, n, x, sum, ierror;
    
    for (n=0; n<chainParams.numRuns; n++)
        {
        for (i=0; i<chainParams.numChains; i++)
            {
            for (j=0; j<chainParams.numChains; j++)
                {
                if (i != j)
                    {
                    if (proc_id == 0)
                        x = 0;
                    else
                        x = swapInfo[n][i][j];
                    ierror = MPI_Allreduce (&x, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                    if (ierror != MPI_SUCCESS)
                        return (ERROR);
                    if (proc_id == 0)
                        swapInfo[n][i][j] += sum;
                    else
                        swapInfo[n][i][j] = 0;
                    }
                }
            }
        }

    return (NO_ERROR);
}


int ReassembleTuningParams (void)
{
    int        i, j, k, lower, ierror;
    MrBFlt     *x, *sum;
    
    x = (MrBFlt *) SafeCalloc (2*numUsedMoves, sizeof(MrBFlt));
    sum = x + numUsedMoves;

    lower = numGlobalChains / num_procs;
    if (numGlobalChains % num_procs != 0)
        lower++;

    for (i=lower; i<numGlobalChains; i++)
        {
        for (j=0; j<numLocalChains; j++)
            {
            if (chainId[j] == i)
                break;
            }

        for (k=0; k<numUsedMoves; k++)
            {
            if (j != numLocalChains && usedMoves[k]->moveType->numTuningParams > 0) /* we have the tuning parameter of interest */
                x[k] = usedMoves[k]->tuningParam[i][0];
            else
                x[k] = 0.0;
            }

        ierror = MPI_Allreduce (x, sum, numUsedMoves, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (ierror != MPI_SUCCESS)
            {
            free (x);
            return (ERROR);
            }

        if (proc_id == 0)
            {
            for (k=0; k<numUsedMoves; k++)
                {
                if (usedMoves[k]->moveType->numTuningParams > 0)
                    usedMoves[k]->tuningParam[i][0] = sum[k];
                }
            }
        }

    free (x);
    return (NO_ERROR);
}


void RedistributeMoveInfo (void)
{
    int         i, j, k;
    MCMCMove    *mv;

    /* Leave if not processor 0, because then we already have the necessary info since
       it was not deleted in ReassembleMoveInfo */
    if (proc_id != 0)
        return;

    /* If we are processor 0, simply delete the unnecessary information */
    for (i=0; i<numGlobalChains; i++)
        {
        for (j=0; j<numLocalChains; j++)
            if (chainId[j] == i)
                break;
        
        if (j == numLocalChains)
            {
            /* we do not have this chain, so delete the move info */
            for (k=0; k<numUsedMoves; k++)
                {
                mv = usedMoves[k];

                /* reset counts */
                mv->nAccepted[i] = 0;
                mv->nTried[i] = 0;
                mv->nBatches[i] = 0;
                mv->nTotAccepted[i] = 0;
                mv->nTotTried[i] = 0;
                mv->lastAcceptanceRate[i] = 0;
        if (mv->moveType->Autotune != NULL)
                    mv->tuningParam[i][0]=0.0;            
                }
            }
        }
}


int RedistributeParamVals (void)
{
    int             i, j, k, orderLen, nBrlens, lower, upper, numChainsForProc, proc, ierror, *y, *order, *nEvents;
    MrBFlt          *x, *brlens, **position, **rateMult;
    MPI_Status      status;
    MPI_Request     request;
    Tree            *tree;
    Param           *p;

    extern MrBFlt   *paramValues;
    extern int      paramValsRowSize;
    extern int      intValsRowSize;

    numChainsForProc = numGlobalChains / num_procs;
    if (numGlobalChains % num_procs > 0)
        lower = upper = numChainsForProc+1;
    else
        lower = upper = numChainsForProc;

    for (proc=1; proc<num_procs; proc++, lower=upper)
        {
        if (proc < numGlobalChains % num_procs)
            upper += numChainsForProc+1;
        else
            upper += numChainsForProc;
        
        /* normal parameter values */
        if (proc_id == 0)
            {
            x = paramValues + 2*paramValsRowSize*lower;
            ierror = MPI_Isend (x, paramValsRowSize*2*(upper-lower), MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            ierror = MPI_Waitall (1, &request, &status);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            if (intValsRowSize > 0)
                {
                y = intValues + 2*intValsRowSize*lower;
                ierror = MPI_Isend (y, intValsRowSize*2*(upper-lower), MPI_INT, proc, 0, MPI_COMM_WORLD, &request);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Waitall (1, &request, &status);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                }
            }
        else if (proc_id == proc)
            {
            x = paramValues;
            ierror = MPI_Irecv (x, paramValsRowSize*2*(upper-lower), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            ierror = MPI_Waitall (1, &request, &status);
            if (ierror != MPI_SUCCESS)
                {
                return (ERROR);
                }
            if (intValsRowSize > 0)
                {
                y = intValues;
                ierror = MPI_Irecv (y, intValsRowSize*2*(upper-lower), MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                ierror = MPI_Waitall (1, &request, &status);
                if (ierror != MPI_SUCCESS)
                    {
                    return (ERROR);
                    }
                }
            }
        
        /* mcmc trees */
        brlens = (MrBFlt *) SafeCalloc (2*numLocalTaxa, sizeof(MrBFlt));
        order  = (int *)    SafeCalloc (2*numLocalTaxa,   sizeof(int));
        for (i=lower; i<upper; i++)
            {
            for (j=0; j<numTrees; j++)
                {
                tree = GetTreeFromIndex(j,0,0);
                orderLen = 2*tree->nIntNodes - 1;
                nBrlens = tree->nNodes - 1;
                if (proc_id == 0)
                    {
                    tree = GetTreeFromIndex(j,i,0);
                    if (tree->isRooted == YES)
                        StoreRTreeWithIndices(tree, order, brlens);
                    else
                        StoreUTree(tree, order, brlens);
                    ierror = MPI_Isend (order, orderLen, MPI_INT, proc, 0, MPI_COMM_WORLD, &request);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (1, &request, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Isend (brlens, nBrlens, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                   ierror = MPI_Waitall (1, &request, &status);
                   if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    }
                else if (proc_id == proc)
                    {
                    tree = GetTreeFromIndex(j,i-lower,0);
                    ierror = MPI_Irecv (order, orderLen, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (1, &request, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Irecv (brlens, nBrlens, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    ierror = MPI_Waitall (1, &request, &status);
                    if (ierror != MPI_SUCCESS)
                        {
                        return (ERROR);
                        }
                    if (tree->isRooted == YES)
                        RetrieveRTreeWithIndices(tree, order, brlens);
                    else
                        {    
                        RetrieveUTree(tree, order, brlens);
                        if (localOutGroup != 0)
                            MoveCalculationRoot(tree,localOutGroup);
                        }
                    /* since we only transferred some info, there are additional things we need to
                       consider, like names, constraints and calibrations */
                    InitializeTreeCalibrations(tree);
                    CheckSetConstraints(tree);
                    SetDatedNodeAges(modelSettings[tree->relParts[0]].brlens, i-lower, 0);
                    strcpy(tree->name, GetTreeFromIndex(j, i, 0)->name);
                    tree = GetTreeFromIndex(j,i-lower,1);
                    strcpy(tree->name, GetTreeFromIndex(j, i, 0)->name);
                    }
                }
            }
        free (brlens);
        free (order);

        /* CPP relaxed clock  parameters */
        for (i=lower; i<upper; i++)
            {
            for (j=0; j<numParams; j++)
                {
                p = &params[j];
                if (p->paramType == P_CPPEVENTS)
                    {
                    if (proc_id == 0)
                        {
                        /* get pointers */
                        nEvents = p->nEvents[2*i];
                        position = p->position[2*i];
                        rateMult = p->rateMult[2*i];

                        /* send number of events */
                        ierror = MPI_Isend (nEvents, 2*numLocalTaxa, MPI_INT, proc, 0, MPI_COMM_WORLD, &request);
                        if (ierror != MPI_SUCCESS)
                            return (ERROR);
                        ierror = MPI_Waitall (1, &request, &status);
                        if (ierror != MPI_SUCCESS)
                            return (ERROR);

                        /* send events and clear pointers */
                        for (k=0; k<2*numLocalTaxa; k++)
                            {
                            if (nEvents[k] > 0)
                                {
                                ierror = MPI_Isend (position[k], nEvents[k], MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                ierror = MPI_Waitall (1, &request, &status);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);

                                ierror = MPI_Isend (rateMult[k], nEvents[k], MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &request);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                ierror = MPI_Waitall (1, &request, &status);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);

                                free(position[k]);
                                free(rateMult[k]);
                                position[k] = NULL;
                                rateMult[k] = NULL;
                                nEvents[k] = 0;
                                }
                            }
                        }
                    else if (proc_id == proc)
                        {
                        /* find pointers */
                        nEvents = p->nEvents[2*(i-lower)];
                        position = p->position[2*(i-lower)];
                        rateMult = p->rateMult[2*(i-lower)];

                        /* clear previous events */
                        for (k=0; k<2*numLocalTaxa; k++)
                            {
                            if (nEvents[k] > 0)
                                {
                                free (position[k]);
                                free (rateMult[k]);
                                position[k] = NULL;
                                rateMult[k] = NULL;
                                nEvents[k] = 0;
                                }
                            }

                        /* receive events */
                        ierror = MPI_Irecv (nEvents, 2*numLocalTaxa, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
                        if (ierror != MPI_SUCCESS)
                            return (ERROR);
                        ierror = MPI_Waitall (1, &request, &status);
                        if (ierror != MPI_SUCCESS)
                            return (ERROR);
                        for (k=0; k<2*numLocalTaxa; k++)
                            {
                            if (nEvents[k] > 0)
                                {
                                position[k] = (MrBFlt *) SafeCalloc (nEvents[k], sizeof(MrBFlt));
                                rateMult[k] = (MrBFlt *) SafeCalloc (nEvents[k], sizeof(MrBFlt));

                                ierror = MPI_Irecv (position[k], nEvents[k], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                ierror = MPI_Waitall (1, &request, &status);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);

                                ierror = MPI_Irecv (rateMult[k], nEvents[k], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                ierror = MPI_Waitall (1, &request, &status);
                                if (ierror != MPI_SUCCESS)
                                    return (ERROR);
                                }
                            }
                       }
                   }
                }
            }

        /* update evolutionary branch lengths or rates (because node indices have changed) */
        if (proc_id == proc)
            {
            for (i=0; i<upper-lower; i++)
                {
                for (j=0; j<numParams; j++)
                    {
                    p = &params[j];
                    if (p->paramType == P_CPPEVENTS)
                        {
                        tree = GetTree(p, i, 0);
                        UpdateCppEvolLengths (p, tree->root->left, i);
                        }
                    else if (p->paramType == P_TK02BRANCHRATES || (p->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(p, i, 0) == RCL_TK02))
                        {
                        tree = GetTree (p, i, 0);
                        UpdateTK02EvolLengths (p, tree, i);
                        }
                    else if (p->paramType == P_IGRBRANCHRATES || (p->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(p, i, 0) == RCL_IGR))
                        {
                        tree = GetTree(p, i, 0);
                        UpdateIgrBrachLengths (p, tree, i);
                        }
                    }
                }
            }
        }

    return (NO_ERROR);
}


int RedistributeTuningParams (void)
{
    int     i, j, k, lower, ierror;
    MrBFlt  *x, *sum;
    
    x = (MrBFlt *) SafeCalloc (2*numUsedMoves, sizeof(MrBFlt));
    sum = x + numUsedMoves;

    lower = numGlobalChains / num_procs;
    if (numGlobalChains % num_procs != 0)
        lower++;

    if (proc_id != 0)
        {
        for (i=0; i<numGlobalChains; i++)
            {
            for (k=0; k<numUsedMoves; k++)
                                {
                                if (usedMoves[k]->moveType->numTuningParams > 0)
                        usedMoves[k]->tuningParam[i][0] = 0.0;
                                }
        }
        }

    for (i=lower; i<numGlobalChains; i++)
        {
        for (j=0; j<numLocalChains; j++)
            {
            if (chainId[j] == i)
                break;
            }

        for (k=0; k<numUsedMoves; k++)
            {
            if (proc_id == 0 && usedMoves[k]->moveType->numTuningParams > 0) /* we have the tuning parameter of interest */
        {
                x[k] = usedMoves[k]->tuningParam[i][0];
        usedMoves[k]->tuningParam[i][0]=0.0;
        }
            else
                x[k] = 0.0;
        }

        ierror = MPI_Allreduce (x, sum, numUsedMoves, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (ierror != MPI_SUCCESS)
            {
            free (x);
            return (ERROR);
            }

        if (j != numLocalChains)   /* we have the chain of interest */
            {
            for (k=0; k<numUsedMoves; k++)
                {
                if (usedMoves[k]->moveType->numTuningParams > 0)
                    usedMoves[k]->tuningParam[i][0] = sum[k];
                }
            }
        }

    free (x);
    return (NO_ERROR);
}

#endif


/* RemovePartition: Remove a partition from the tree keeping track of partition frequencies */
int RemovePartition (PFNODE *r, BitsLong *p, int runId)
{
    int     i, comp;
    
    if (r == NULL)
        {
        return (ERROR);
        }
    else
        {
        for (i=0; i<nLongsNeeded; i++)
            {
            if (r->partition[i] != p[i])
                break;
            }
        
        if (i == nLongsNeeded)
            comp = 0;
        else if (r->partition[i] < p[i])
            comp = -1;
        else
            comp = 1;
        
        if (comp == 0)          /* match */
            {
            if (r->count[runId] == 0)
                return ERROR;
            else
                r->count[runId]--;
            }
        else if (comp < 0)      /* greater than -> into left subtree */
            {
            if ((RemovePartition (r->left, p, runId)) == ERROR)
                return ERROR;
            }
        else
            {
            /* less than -> into right subtree */
            if ((RemovePartition (r->right, p, runId)) == ERROR)
                return ERROR;
            }
        }

    return (NO_ERROR);
}


/* RemoveTreeFromPartitionCounters: Break a tree into partitions and remove those from counters */
int RemoveTreeFromPartitionCounters (Tree *tree, int treeId, int runId)
{
    int         i, j, nTaxa;
    TreeNode    *p;

    if (tree->isRooted == YES)
        nTaxa = tree->nNodes - tree->nIntNodes - 1;
    else
        nTaxa = tree->nNodes - tree->nIntNodes;

    for (i=0; i<nTaxa; i++)
        {
        ClearBits(partition[i], nLongsNeeded);
        SetBit(i, partition[i]);
        }

    for (i=0; i<tree->nIntNodes-1; i++)
        {
        p = tree->intDownPass[i];
        assert (p->index >= tree->nNodes - tree->nIntNodes - (tree->isRooted == YES ? 1 : 0));
        for (j=0; j<nLongsNeeded; j++)
            {
            partition[p->index][j] = partition[p->left->index][j] | partition[p->right->index][j];
            }
        
        if ((RemovePartition (partFreqTreeRoot[treeId], partition[p->index], runId)) == ERROR)
            {
            MrBayesPrint ("%s   Could not remove partition %d in RemoveTreeFromPartitionCounters\n", spacer, p->index);
            ShowParts(stdout,partition[p->index],numLocalTaxa);
            return ERROR;
            }
        }

    return NO_ERROR;
}


/* RemoveTreeSamples: Remove tree samples from partition counters */
int RemoveTreeSamples (int from, int to)
{
    int         i, j, k, longestLine, line;
    char        c, *s, *lineBuf=0;
    FILE        *fp;
    Tree        *t;
    TreeList    *treeList;
    char        *tempStr;
    int         tempStrSize = TEMPSTRSIZE;

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return (NO_ERROR);
#   endif

    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    if (chainParams.saveTrees == YES)
        {
        for (i=0; i<numTopologies; i++)
            {
            t = chainParams.dtree;
            if (topologyParam[i]->tree[0]->isRooted == YES)
                t->isRooted = YES;
            else
                t->isRooted = NO;

            for (j=0; j<chainParams.numRuns; j++)
                {               
                treeList = &chainParams.treeList[j*numTopologies + i];
                
                for (k=from; k<=to; k++)
                    {
                    GetFromTreeList (treeList, t);
                    if (RemoveTreeFromPartitionCounters (t, i, j) == ERROR)
                        {
                        SafeFclose (&fp);
                        free (tempStr);
                        free (lineBuf);
                        return (ERROR);
                        }
                    }
                }
            }
        }
    else
        {
        for (i=0; i<numTopologies; i++)
            {
            t = chainParams.dtree;

            for (j=0; j<chainParams.numRuns; j++)
                {
                if (numPrintTreeParams == 1)
                    {
                    if (chainParams.numRuns == 1)
                        SafeSprintf (&tempStr, &tempStrSize, "%s.t", chainParams.chainFileName);
                    else
                        SafeSprintf (&tempStr, &tempStrSize, "%s.run%d.t", chainParams.chainFileName, j+1);
                    }
                else
                    {
                    if (chainParams.numRuns == 1)
                        SafeSprintf (&tempStr, &tempStrSize, "%s.tree%d.t", chainParams.chainFileName, topologyPrintIndex[i]+1);
                    else
                        SafeSprintf (&tempStr, &tempStrSize, "%s.tree%d.run%d.t", chainParams.chainFileName, topologyPrintIndex[i]+1, j+1);
                    }

                if ((fp = OpenTextFileR (tempStr)) == NULL)
                        {
                        free (tempStr);
                        return (ERROR);
                        }
                
                if (from == 1)
                    {
                    longestLine = LongestLine(fp);
                    lineBuf = (char *) SafeCalloc (longestLine+2,sizeof(char));
                    fseek (fp, LastBlock(fp, lineBuf, longestLine), SEEK_SET);
                    fseek (fp, FirstTree(fp, lineBuf, longestLine), SEEK_SET);
                    fgetpos (fp, &chainParams.tFilePos[j*numTopologies+i]);
                    free(lineBuf);
                    lineBuf = NULL;
                    }
                fsetpos (fp, &chainParams.tFilePos[j*numTopologies+i]);

                longestLine = 0;
                for (k=0; k<=to-from; k++)
                    {
                    line = 0;
                    do {
                        line++;
                        } while ((c = fgetc(fp)) != '\r' && c != '\n');
                    
                    if (line > longestLine)
                        longestLine = line;
                    
                    while ((c = fgetc(fp)) == '\r' || c == '\n')
                        ;
                    }

                lineBuf = (char *) SafeCalloc (longestLine + 10, sizeof (char));
                if (!lineBuf)
                    {
                    SafeFclose (&fp);
                    free (tempStr);
                    return (ERROR);
                    }
                fsetpos (fp, &chainParams.tFilePos[j*numTopologies+i]);
                /* The fsetpos and fgetpos pair are affected by writing to the file,
                    at least on Windows systems. The effect is to put the subsequent
                    fsetpos back a few positions. The following code will deal with this
                    problem without affecting systems where this does not happen. */
                do { c = fgetc(fp);
                } while (c != 't');

                for (k=from; k<=to; k++)
                    {
                    if (fgets (lineBuf, longestLine + 5, fp) == NULL) 
                        {
                        free (tempStr);
                        free (lineBuf);
                        return ERROR;
                        }

                    s = strtok (lineBuf, " ");
                    for (s = strtok (NULL, ";"); *s != '('; s++)
                        ;
                    
                    StripComments (s);
                    if (ResetTopology (t, s) == ERROR)
                        {
                        SafeFclose (&fp);
                        free (tempStr);
                        free (lineBuf);
                        return (ERROR);
                        }
                    
                    if (RemoveTreeFromPartitionCounters (t, i, j) == ERROR)
                        {
                        SafeFclose (&fp);
                        free (tempStr);
                        free (lineBuf);
                        return (ERROR);
                        }
                    }
                fgetpos (fp, &chainParams.tFilePos[j*numTopologies+i]);
                free (lineBuf);
                SafeFclose (&fp);
                }
            }
        }

    /* remove unnecessary nodes from the tree holding partition counters */
    for (i=0; i<numTopologies; i++)
        {
        partFreqTreeRoot[i] = CompactTree (partFreqTreeRoot[i]);
        }

    free (tempStr);
    return (NO_ERROR);
}


int ReopenMBPrintFiles (void)
{
    int     i, n;
    char    fileName[120], localFileName[100];
    
    /* Take care of the mpi procs that do not have a file */
#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return (NO_ERROR);
#   endif

    /* Get root of local file name */
    strcpy (localFileName, chainParams.chainFileName);

    /* Reopen the .p and .t files */
    for (n=0; n<chainParams.numRuns; n++)
        {
        if (chainParams.numRuns == 1)
            sprintf (fileName, "%s.p", localFileName);
        else
            sprintf (fileName, "%s.run%d.p", localFileName, n+1);

        if ((fpParm[n] = OpenTextFileA (fileName)) == NULL)
            return (ERROR);

        for (i=0; i<numTrees; i++)
            {
            if (numTrees == 1 && chainParams.numRuns == 1)
                sprintf (fileName, "%s.t", localFileName);
            else if (numTrees > 1 && chainParams.numRuns == 1)
                sprintf (fileName, "%s.tree%d.t", localFileName, i+1);
            else if (numTrees == 1 && chainParams.numRuns > 1)
                sprintf (fileName, "%s.run%d.t", localFileName, n+1);
            else
                sprintf (fileName, "%s.tree%d.run%d.t", localFileName, i+1, n+1);

            if ((fpTree[n][i] = OpenTextFileA (fileName)) == NULL)
                return (ERROR);
            }
        }

    /* Reopen the .mcmc file */
    if (chainParams.mcmcDiagn == YES)
        {
        sprintf (fileName, "%s.mcmc", localFileName);

        if ((fpMcmc = OpenTextFileA (fileName)) == NULL)
            return (ERROR);
        }
    
#   if defined (PRINT_DUMP)
    for (n=0; n<chainParams.numRuns; n++)
        {
        if (chainParams.numRuns == 1)
            sprintf (fileName, "%s.dump", localFileName);
        else
            sprintf (fileName, "%s.run%d.dump", localFileName, n+1);

        if ((fpDump[n] = OpenTextFileA (fileName)) == NULL)
            return (ERROR);
        }
#   endif

    return (NO_ERROR);
}


int ConfirmAbortRun(void)
{
    char c, line[100];
    int  ret=0, i;

    /* reset requestAbortRun */
    requestAbortRun = NO;

    MrBayesPrint("   Do you really want to stop the run (y/n)?");
    if (fgets (line,98,stdin) == NULL)
        {
        printf ("Error in function: %s at line: %d in file: %s", __FUNCTION__, __LINE__, __FILE__);
        }
    for (i=0; (c=line[i])!='\0' && !isgraph(c); i++)
        ;
    if (c == 'y' || c == 'Y')
        ret=1;
    else 
        {
        MrBayesPrint("   Mcmc run continued ...\n\n");
        ret=0;
        }
    return ret;
}


/*-------------------------------------------------------------------
|
|   ResetChainIds: Make sure parameter values are swapped back 
|       at the end of a Metropolis-coupled MCMC run
|
--------------------------------------------------------------------*/
void ResetChainIds (void)
{
    int     j, k, k1, tempId, *curId, toChn, fromChn, *to, *from, *swap;
    Param   *p;
    MrBFlt  *fromVals, *toVals, *swapVals, **fromPosition, **toPosition,
            **swapPosition, **fromRateMult, **toRateMult, **swapRateMult;
    Tree    *toTree, *fromTree, *swapTree;

    curId = (int *) SafeCalloc (numGlobalChains, sizeof(int));

#if defined (MPI_ENABLED)
    ReassembleParamVals(curId);
    ReassembleTuningParams();
    SetChainIds();
    if (proc_id != 0)
    {
        /* reset state */
        for (j=0; j<numLocalChains; j++)
            state[j] = 0;
        free (curId);
    return;
        }
#else
    for (j=0; j<numGlobalChains; j++)
        curId[j] = chainId[j];
    SetChainIds();
#endif
    
    for (toChn=0; toChn<numGlobalChains; toChn++)
        {
        if (curId[toChn] == toChn)
            {
            if (state[toChn] != 0)
                {
                CopyParams (toChn);
                CopyTrees (toChn);
                state[toChn] ^= 1;
                }
            continue;
            }

        /* we need to swap all values */
        /* first find the chain to swap with */
        for (j=toChn+1; j<numGlobalChains; j++)
            if (curId[j] == toChn)
                break;
        fromChn = j;

        /* normal params */
        CopyParams (toChn);
        CopyTrees (toChn);
        CopyParams (fromChn);
        CopyTrees (fromChn);

        for (j=0; j<numParams; j++)
            {
            p = &params[j];
            toVals = GetParamVals (p, toChn, 0);
            swapVals = GetParamVals (p, toChn, 1);
            fromVals = GetParamVals (p, fromChn, state[fromChn]);
            for (k=0; k<p->nValues; k++)
                {
                toVals[k] = fromVals[k];
                fromVals[k] = swapVals[k];
                }
            toVals = GetParamSubVals (p, toChn, 0);
            swapVals = GetParamSubVals (p, toChn, 1);
            fromVals = GetParamSubVals (p, fromChn, state[fromChn]);
            for (k=0; k<p->nSubValues; k++)
                {
                toVals[k] = fromVals[k];
                fromVals[k] = swapVals[k];
                }
            if (p->nStdStateFreqs > 0)
                {
                toVals = GetParamStdStateFreqs (p, toChn, 0);
                swapVals = GetParamStdStateFreqs (p, toChn, 1);
                fromVals = GetParamStdStateFreqs (p, fromChn, state[fromChn]);
                for (k=0; k<p->nStdStateFreqs; k++)
                    {
                    toVals[k] = fromVals[k];
                    fromVals[k] = swapVals[k];
                    }
                }
            }

        /* mcmc trees */
        for (j=0; j<numTrees; j++)
            {
            toTree = GetTreeFromIndex(j, toChn, 0);
            swapTree = GetTreeFromIndex(j, toChn, 1);
            fromTree = GetTreeFromIndex(j, fromChn, state[fromChn]);
            CopyToTreeFromTree (toTree, fromTree);
            CopyToTreeFromTree (fromTree, swapTree);
            CopyToTreeFromTree (swapTree, toTree);
            swapTree = GetTreeFromIndex(j, fromChn, state[fromChn] ^ 1);
            CopyToTreeFromTree (swapTree, fromTree);
            }
        /* CPP relaxed clock params */
        for (j=0; j<numParams; j++)
            {
            p = &params[j];
            if (p->paramType == P_CPPEVENTS)
                {
                to = p->nEvents[2*toChn];
                swap = p->nEvents[2*toChn+1];
                from = p->nEvents[2*fromChn+state[fromChn]];
                toPosition = p->position[2*toChn];
                swapPosition = p->position[2*toChn+1];
                fromPosition = p->position[2*fromChn+state[fromChn]];
                toRateMult = p->rateMult[2*toChn];
                swapRateMult = p->rateMult[2*toChn+1];
                fromRateMult = p->rateMult[2*fromChn+state[fromChn]];
                for (k=0; k<2*numLocalTaxa; k++)
                    {
                    if (from[k] > 0)
                        {
                        toPosition[k] = (MrBFlt *) SafeRealloc ((void *)toPosition[k], from[k]*sizeof (MrBFlt));
                        toRateMult[k] = (MrBFlt *) SafeRealloc ((void *)toRateMult[k], from[k]*sizeof (MrBFlt));
                        for (k1=0; k1<from[k]; k1++)
                            {
                            toPosition[k][k1] = fromPosition[k][k1];
                            toRateMult[k][k1] = fromRateMult[k][k1];
                            }
                        }
                    else if (to[k] > 0)
                        {
                        free (toPosition[k]);
                        toPosition[k] = NULL;
                        free (toRateMult[k]);
                        toRateMult[k] = NULL;
                        }
                    to[k] = from[k];
                    if (swap[k] > 0)
                        {
                        fromPosition[k] = (MrBFlt *) SafeRealloc ((void *)fromPosition[k], swap[k]*sizeof (MrBFlt));
                        fromRateMult[k] = (MrBFlt *) SafeRealloc ((void *)fromRateMult[k], swap[k]*sizeof (MrBFlt));
                        for (k1=0; k1<swap[k]; k1++)
                            {
                            fromPosition[k][k1] = swapPosition[k][k1];
                            fromRateMult[k][k1] = swapRateMult[k][k1];
                            }
                        }
                    else if (from[k] > 0)
                        {
                        free (fromPosition[k]);
                        fromPosition[k] = NULL;
                        free (fromRateMult[k]);
                        fromRateMult[k] = NULL;
                        }
                    from[k] = swap[k];
                    }
                }   
            }
        /* reset state of chain */
        state[toChn] = 0;

        /* make sure that the id is correct for the from chain */
        tempId = curId[toChn];
        curId[fromChn] = curId[toChn];
        curId[toChn] = tempId;
        }

    free (curId);
}


/* ResetFlips: Reset flipped cond likes etc after rejection */
void ResetFlips (int chain)
{
    int         d, i;
    ModelInfo   *m;
    TreeNode    *p;
    Tree        *tree;
#if defined (BEAGLE_ENABLED)
    int         *isScalerNode=NULL;
#endif    
    
    for (d=0; d<numCurrentDivisions; d++)
    {
        m = &modelSettings[d];
#if defined (BEAGLE_ENABLED)
        if (m->useBeagle == YES)
            isScalerNode = m->isScalerNode[chain];
#endif
        if (m->upDateCl != YES)
            continue;
        
#if defined (BEAGLE_ENABLED)
        if (m->useBeagle == NO || 
            beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS ||
            m->rescaleBeagleAll == YES)
                {
                FlipSiteScalerSpace (m, chain);
                if (m->useBeagle == YES && m->rescaleBeagleAll == YES)
                    m->rescaleFreq[chain] = m->rescaleFreqOld;
                }
#else
        FlipSiteScalerSpace (m, chain);
#endif
            
        
        if (m->upDateCijk == YES && m->nCijkParts > 0)
            FlipCijkSpace (m, chain);
        
        /* cycle over tree */
        tree = GetTree (m->brlens, chain, state[chain]);
        for (i=0; i<tree->nNodes; i++)
            {
            p = tree->allDownPass[i];
            if (p->upDateTi == YES)
                FlipTiProbsSpace (m, chain, p->index);
            if (p->right != NULL)    /* do not flip terminals in case these flags are inappropriately set by moves */
                {
                if (p->upDateCl == YES)
                    {
                    FlipCondLikeSpace (m, chain, p->index);
#if defined (BEAGLE_ENABLED)
                    if (m->useBeagle == NO || 
                        beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS ||
                        (m->rescaleBeagleAll == YES && isScalerNode[p->index] == YES))
                        FlipNodeScalerSpace (m, chain, p->index);
#else
                    FlipNodeScalerSpace (m, chain, p->index);
#endif
                    }
#if defined (BEAGLE_ENABLED)
                else if (m->rescaleBeagleAll == YES)
                    {
                    FlipCondLikeSpace (m, chain, p->index);
                    if (isScalerNode[p->index] == YES)
                        FlipNodeScalerSpace (m, chain, p->index);
                    }
#endif
                }
            }
        
        /* division flag and tree node flags are reset when trees are copied */
    }
}


/*-------------------------------------------------------------------
|
|   ResetScalersPartition: reset scaler nodes of the given tree by appropriately setting isScalerNode array.
| @param isScalerNode   is an array which gets set with information about scaler node. 
|                       For each internal node isScalerNode[node->index] is set to YES if it has to be scaler node.
|                       Note: Only internal nodes can become scaler nodes thus isScalerNode is set only for elemnts in interval [numLocalTaxa, numLocalTaxa+t->nIntNodes]
|
| @param rescaleFreq    effectively represent gaps between rescaling, higher number means more sparse choice of rescaling nodes 
|
--------------------------------------------------------------------*/
int ResetScalersPartition (int *isScalerNode, Tree* t, unsigned rescaleFreq)
{
    int         n;
    TreeNode    *p;
        
    /* set the node depth value of terminal nodes to zero; reset scalerNode */
    for (n=0; n<t->nNodes; n++)
        {
        p = t->allDownPass[n];
        if (p->left == NULL)
            p->x = 0;
        }

    /* loop over interior nodes */
    for (n=0; n<t->nIntNodes; n++)
        {
        p = t->intDownPass[n];
        assert (((p->index - numLocalTaxa) >= 0) && ((p->index - numLocalTaxa) < t->nIntNodes));
        p->x = p->left->x + p->right->x + 1;

        if (p->x > 2 * (int)rescaleFreq)
            {
            assert (p->left->left != NULL && p->right->left != NULL);
            isScalerNode[p->left->index] = YES;
            p->left->x = 0;
            isScalerNode[p->right->index] = YES;
            p->right->x = 0;
            p->x = 1;
            }
        else if (p->x > (int)rescaleFreq)
            {
            if (p->left->x > p->right->x)
                {
                assert (p->left->left != NULL);
                isScalerNode[p->left->index] = YES;
                p->left->x = 0;
                }
            else
                {
                assert (p->right->left != NULL);
                isScalerNode[p->right->index] = YES;
                p->right->x = 0;
                }
            p->x = p->left->x + p->right->x + 1;
            }
        else
            isScalerNode[p->index] = NO;
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------------
 |
 |   ResetSiteScalers: Set log site scalers to 0.0.
 |
 ------------------------------------------------------------------------*/
void ResetSiteScalers (ModelInfo *m, int chain)
{
    int     c;
    CLFlt   *lnScaler;

#if defined (BEAGLE_ENABLED)
    if (m->useBeagle == YES)
        {
        beagleResetScaleFactors(m->beagleInstance, m->siteScalerIndex[chain]);
        return;
        }
#endif
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    for (c=0; c<m->numChars; c++)
        lnScaler[c] = 0.0;
}


/*----------------------------------------------------------------------
|
|   ReusePreviousResults: Save old .p, .t, .ss and .mcmc files with ~ extension,
|      then prepare new print files with the relevant old values added in
|      The number of samples is returned in numSamples
|
------------------------------------------------------------------------*/
int ReusePreviousResults (int *numSamples, int steps)
{
    int         i, n;
    char        localFileName[100], fileName[220], bkupName[220];

    (*numSamples) = 0;

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return (NO_ERROR);
#   endif

    /* Allocate space for file pointers */
    if (memAllocs[ALLOC_FILEPOINTERS] == YES)
        {
        MrBayesPrint ("%s   File pointers already allocated in ReusePreviousResults\n", spacer);
        return ERROR;
        }
    fpMcmc = NULL;
    fpSS = NULL;
    fpParm = NULL;
    fpTree = NULL;  
    fpParm = (FILE **) SafeCalloc (chainParams.numRuns, sizeof (FILE *));
    if (fpParm == NULL)
        {
        MrBayesPrint ("%s   Could not allocate fpParm in ReusePreviousResults\n", spacer);
        return ERROR;
        }
    memAllocs[ALLOC_FILEPOINTERS] = YES;
    fpTree = (FILE ***) SafeCalloc (chainParams.numRuns, sizeof (FILE **));
    if (fpTree == NULL)
        {
        MrBayesPrint ("%s   Could not allocate fpTree in ReusePreviousResults\n", spacer);
        return ERROR;
        }
    fpTree[0] = (FILE **) SafeCalloc (numTrees*chainParams.numRuns, sizeof (FILE *));
    if (fpTree[0] == NULL)
        {
        MrBayesPrint ("%s   Could not allocate fpTree[0] in ReusePreviousResults\n", spacer);
        return ERROR;
        }
    for (i=1; i<chainParams.numRuns; i++)
        fpTree[i] = fpTree[0] + i*numTrees;

    /* Get root of local file name */
    strcpy (localFileName, chainParams.chainFileName);

    /* Store old and prepare new .p and .t files */
    for (n=0; n<chainParams.numRuns; n++)
        {
        if (chainParams.numRuns == 1)
            sprintf (fileName, "%s%s.p", workingDir, localFileName);
        else
            sprintf (fileName, "%s%s.run%d.p", workingDir, localFileName, n+1);
        strcpy(bkupName,fileName);
        strcat(bkupName,"~");
        remove(bkupName);
        if (rename(fileName,bkupName) != 0)
            {
            MrBayesPrint ("%s   Could not rename file %s\n", spacer, fileName);
            return ERROR;
            }

        if ((fpParm[n] = OpenNewMBPrintFile (fileName+strlen(workingDir))) == NULL)
            return (ERROR);
        else if (CopyResults(fpParm[n],bkupName+strlen(workingDir),numPreviousGen) == ERROR)
            return (ERROR);

        for (i=0; i<numTrees; i++)
            {
            if (numTrees == 1 && chainParams.numRuns == 1)
                sprintf (fileName, "%s%s.t", workingDir, localFileName);
            else if (numTrees > 1 && chainParams.numRuns == 1)
                sprintf (fileName, "%s%s.tree%d.t", workingDir, localFileName, i+1);
            else if (numTrees == 1 && chainParams.numRuns > 1)
                sprintf (fileName, "%s%s.run%d.t", workingDir, localFileName, n+1);
            else
                sprintf (fileName, "%s%s.tree%d.run%d.t", workingDir, localFileName, i+1, n+1);
            strcpy(bkupName,fileName);
            strcat(bkupName,"~");
            remove(bkupName);
            if (rename(fileName,bkupName) != 0)
                {
                MrBayesPrint ("%s   Could not rename file %s\n", spacer, fileName);
                return ERROR;
                }
            if ((fpTree[n][i] = OpenNewMBPrintFile (fileName+strlen(workingDir))) == NULL)
                return (ERROR);
            else if (CopyTreeResults(fpTree[n][i],bkupName+strlen(workingDir),numPreviousGen,numSamples) == ERROR)
                return (ERROR);
            }
        }

    /* Store old and prepare new .ss file */
    if (chainParams.isSS == YES)
        {
        sprintf (fileName, "%s%s.ss", workingDir, chainParams.chainFileName);
        strcpy(bkupName,fileName);
        strcat(bkupName,"~");
        remove(bkupName);
        if (rename(fileName,bkupName) != 0)
            {
            MrBayesPrint ("%s   Could not rename file %s\n", spacer, fileName);
            return ERROR;
            }
        if ((fpSS = OpenNewMBPrintFile (fileName+strlen(workingDir))) == NULL)
            return (ERROR);
        else if (CopyProcessSsFile(fpSS,bkupName+strlen(workingDir),steps,marginalLnLSS,splitfreqSS)==ERROR)
            return (ERROR);
        }

    /* Store old and prepare new .mcmc file */
    if (chainParams.mcmcDiagn == YES)
        {
        sprintf (fileName, "%s%s.mcmc", workingDir, chainParams.chainFileName);
        strcpy(bkupName,fileName);
        strcat(bkupName,"~");
        remove(bkupName);
        if (rename(fileName,bkupName) != 0)
            {
            MrBayesPrint ("%s   Could not rename file %s\n", spacer, fileName);
            return ERROR;
            }
        if ((fpMcmc = OpenNewMBPrintFile (fileName+strlen(workingDir))) == NULL)
            return (ERROR);
        else if (CopyResults(fpMcmc,bkupName+strlen(workingDir),numPreviousGen)==ERROR)
            return (ERROR);
        }
    
#   if defined (PRINT_DUMP)
    fpDump = (FILE **) SafeCalloc (chainParams.numRuns, sizeof (FILE *));
    
    for (n=0; n<chainParams.numRuns; n++)
        {
        if (chainParams.numRuns == 1)
            sprintf (fileName, "%s.dump", localFileName);
        else
            sprintf (fileName, "%s.run%d.dump", localFileName, n+1);
        
        if ((fpDump[n] = OpenTextFileA (fileName)) == NULL)
            return (ERROR);
    }
#   endif

    return (NO_ERROR);
}


int RunChain (RandLong *seed)
{
    int         i, j, n, chn, swapA=0, swapB=0, whichMove, acceptMove;
    int         lastDiagnostics;    // the sample no. when last diagnostic was performed
    int         removeFrom, removeTo=0;
    int         stopChain, nErrors;
    MrBFlt      r=0.0, lnLikelihoodRatio, lnPriorRatio, lnProposalRatio, lnLike=0.0, lnPrior=0.0, f=0.0, CPUTime;
    MCMCMove    *theMove, *mv;
    time_t      startingT, endingT, stoppingT1, stoppingT2;
    clock_t     previousCPUTime, currentCPUTime;
    /* Stepping-stone sampling variables */
    int         run, samplesCountSS=0, stepIndexSS=0, numGenInStepSS=0, numGenOld, lastStepEndSS=0, numGenInStepBurninSS=0;
    MrBFlt      stepLengthSS=0, meanSS, varSS, *tempX;
    char        ckpFileName[220], bkupFileName[220];

#   if defined (BEAGLE_ENABLED)
    int         ResetScalersNeeded;  //set to YES if we need to reset node->scalerNode, used in old style rescaling;
#       ifdef DEBUG_BEAGLE
    int         beagleScalingSchemeOld;
#       endif
    ModelInfo   *m;
    ResetScalersNeeded = NO;
    
    for (i=0; i<numCurrentDivisions; i++)
        {
        m = &modelSettings[i];
        if (m->useBeagle == NO || beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS)
            {
            ResetScalersNeeded =YES;
            break;
            }
        }
#   endif

#   if defined (MPI_ENABLED)
    int         ierror, sumErrors;
    MrBFlt      best, sum=0.0;
    MPI_Status  status;
#   endif
#   if defined (DEBUG_RUNCHAIN)
    ModelInfo   *m;
#   endif
    
    /* set nErrors to 0 */
    nErrors = 0;
    if (numLocalTaxa < 4)
        {
        for (i=0; i<numTrees; i++)
            if (GetTreeFromIndex(i, 0, 0)->isRooted == NO)
                break;
        if (i < numTrees && numLocalTaxa < 4)
            {
            MrBayesPrint ("%s   There must be at least four taxa in the analysis\n", spacer);
            return (ERROR);
            }
        else if (i == numTrees && numLocalTaxa < 3)
            {
            MrBayesPrint ("%s   There must be at least three taxa in the analysis\n", spacer);
            return (ERROR);
            }
        }

    /* allocate some memory for the chains */
    if (memAllocs[ALLOC_CURLNL] == YES)
        {
        MrBayesPrint ("%s   curLnL is already allocated\n", spacer);
        nErrors++;
        }
    else if ((curLnL = (MrBFlt *)SafeMalloc((size_t)numLocalChains * sizeof(MrBFlt))) == NULL)
        {
        MrBayesPrint ("%s   Problem allocating curLnL (%d)\n", spacer, numLocalChains * sizeof(MrBFlt));
        nErrors++;
        }
    else if ((maxLnL0 = (MrBFlt *) SafeCalloc ((size_t)(chainParams.numRuns) * (size_t)(chainParams.numChains), sizeof(MrBFlt))) == NULL)
        {
        MrBayesPrint ("%s   Problem allocating maxLnL0\n", spacer, numLocalChains * sizeof(MrBFlt));
        free (curLnL);
        nErrors++;
        }
    else
        memAllocs[ALLOC_CURLNL] = YES;
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   Memory allocation error on at least one processor\n", spacer);
        return ERROR;
        }
#   else
    if (nErrors > 0)
        return ERROR;
#   endif

    if (memAllocs[ALLOC_CURLNPR] == YES)
        {
        MrBayesPrint ("%s   curLnPr is already allocated\n", spacer);
        nErrors++;
        }
    else if ((curLnPr = (MrBFlt *)SafeMalloc((size_t)numLocalChains * sizeof(MrBFlt))) == NULL)
        {
        MrBayesPrint ("%s   Problem allocating curLnPr (%d)\n", spacer, numLocalChains * sizeof(MrBFlt));
        nErrors++;
        }
    else
        memAllocs[ALLOC_CURLNPR] = YES;
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   Memory allocation error on at least one processor\n", spacer);
        return ERROR;
        }
#   else
    if (nErrors > 0)
        return ERROR;
#   endif

    if (memAllocs[ALLOC_CHAINID] == YES)
        {
        MrBayesPrint ("%s   chainId is already allocated\n", spacer);
        nErrors++;
        }
    else if ((chainId = (int *)SafeMalloc((size_t)numLocalChains * sizeof(int))) == NULL)
        {
        MrBayesPrint ("%s   Problem allocating chainId (%d)\n", spacer, numLocalChains * sizeof(int));
        nErrors++;
        }
    else
        memAllocs[ALLOC_CHAINID] = YES;
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   Memory allocation error on at least one processor\n", spacer);
        return ERROR;
        }
#   else
    if (nErrors > 0)
        return ERROR;
#   endif

    if (memAllocs[ALLOC_SWAPINFO] == YES)
        {
        MrBayesPrint ("%s   swapInfo is already allocated\n", spacer);
        nErrors++;
        }
    else if ((swapInfo = (int ***) SafeCalloc (chainParams.numRuns, sizeof (int **))) == NULL)
        {
        MrBayesPrint ("%s   Problem allocating swapInfo\n", spacer);
        nErrors++;
        }
    else
        {
        for (n=0; n<chainParams.numRuns; n++)
            {
            swapInfo[n] = AllocateSquareIntegerMatrix (chainParams.numChains);
            if (!swapInfo[n])
                {
                MrBayesPrint ("%s   Problem allocating swapInfo[%d]\n", spacer, n);
                for (i=0; i<n; i++)
                    free (swapInfo[i]);
                free (swapInfo);
                nErrors++;
                break;
                }
            }
        memAllocs[ALLOC_SWAPINFO] = YES;
        }
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   Memory allocation error on at least one processor\n", spacer);
        return ERROR;
        }
#   else
    if (nErrors > 0)
        return ERROR;
#   endif

    for (n=0; n<chainParams.numRuns; n++)
        for (i=0; i<chainParams.numChains; i++)
            for (j=0; j<chainParams.numChains; j++)
                swapInfo[n][i][j] = 0;

    /* set up counters for topological convergence diagnostics */
    /* allocate tree used for some topological convergence diagnostics */
    if (chainParams.mcmcDiagn == YES && chainParams.numRuns > 1)
        {
        if (SetUpPartitionCounters () == ERROR)
            nErrors++;
#   if defined (MPI_ENABLED)
        if (proc_id == 0)
            {
#   endif
        if (chainParams.relativeBurnin == YES)
            {
            /* we have to remove trees later on */
            if (chainParams.saveTrees == YES)
                {
                if (chainParams.treeList)
                    nErrors++;
                else 
                    {
                    chainParams.treeList = (TreeList *) SafeCalloc (chainParams.numRuns*numTopologies, sizeof (TreeList));
                    if (!chainParams.treeList)
                        nErrors++;
                    }
                if (nErrors == 0)
                    memAllocs[ALLOC_TREELIST] = YES;
                if (noWarn == YES)
                    chainParams.stopTreeGen = (int) (chainParams.numGen * chainParams.burninFraction);
                else
                    chainParams.stopTreeGen = chainParams.numGen;
                }
            else /* if (chainParams.saveTrees == NO) */
                {
                chainParams.tFilePos = (fpos_t*) SafeRealloc ((void *)chainParams.tFilePos, chainParams.numRuns*numTopologies*sizeof (fpos_t));
                if (!chainParams.tFilePos)
                    nErrors++;
                else
                    memAllocs[ALLOC_TFILEPOS] = YES;
                }
            }
#   if defined (MPI_ENABLED)
            }
#   endif
#   if defined (MPI_ENABLED)
        if (proc_id == 0)
            {
            if ((chainParams.stat = (STATS *) SafeCalloc (numTopologies, sizeof (STATS))) == NULL)
                nErrors++;
            else
                {
                memAllocs[ALLOC_STATS] = YES;
                for (i=0; i<numTopologies; i++)
                    chainParams.stat[i].pair = NULL;
                }

            if ((chainParams.dtree = AllocateTree (numLocalTaxa)) == NULL)
                {
                nErrors++;
                }
            else
                memAllocs[ALLOC_DIAGNTREE] = YES;

            if (chainParams.allComps == YES)
                {
                for (i=0; i<numTopologies; i++)
                    {
                    if ((chainParams.stat[i].pair = AllocateSquareDoubleMatrix (chainParams.numRuns)) == NULL)
                        {
                        nErrors++;
                        break;
                        }
                    }
                }
            }
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Memory allocation error on at least one processor\n", spacer);
            return ERROR;
            }
#   else
        if ((chainParams.stat = (STATS *) SafeCalloc (numTopologies, sizeof (STATS))) == NULL)
            return ERROR;
        else
            {
            memAllocs[ALLOC_STATS] = YES;
            for (i=0; i<numTopologies; i++)
                chainParams.stat[i].pair = NULL;
            }
            
        if ((chainParams.dtree = AllocateTree (numLocalTaxa)) == NULL)
            {
            MrBayesPrint ("%s   Could not allocate chainParams.dtree in RunChain\n", spacer);
            return ERROR;
            }
        else
            memAllocs[ALLOC_DIAGNTREE] = YES;
        
        if (chainParams.allComps == YES)
            {
            for (i=0; i<numTopologies; i++)
                {
                if ((chainParams.stat[i].pair = AllocateSquareDoubleMatrix (chainParams.numRuns)) == NULL)
                    {
                    MrBayesPrint ("%s   Could not allocate chainParams.stat.pair in RunChain\n", spacer);
                    return ERROR;
                    }
                }
            }
#   endif
        }

    /* get chain IDs */
    SetChainIds ();
    
    /* distribute parameter starting values and tuning parameters for MPI version */
#   if defined (MPI_ENABLED)
    RedistributeParamVals();
    RedistributeTuningParams();
#   endif

#   if defined (TIMING_ANALIZ)
    CPUCondLikeDown = 0;
    CPUScalers = 0;
    CPUScalersRemove = 0;
    CPUCondLikeRoot = 0;
    CPULilklihood = 0;
#   endif

    /* initialize likelihoods and prior                  */
    /* touch everything and calculate initial cond likes */
    TouchAllPartitions ();
    for (chn=0; chn<numLocalChains; chn++)
        {
        if (chn % chainParams.numChains == 0)
            {
            if (chainParams.numRuns == 1)
                MrBayesPrint ("\n%s   Initial log likelihoods and log prior probs:\n", spacer);
            else
                MrBayesPrint ("\n%s   Initial log likelihoods and log prior probs for run %d:\n", spacer, chn / chainParams.numChains + 1);
            }
        TouchAllTrees (chn);
        TouchAllCijks (chn);
        curLnL[chn] = LogLike(chn);
        curLnPr[chn] = LogPrior(chn);
        for (i=0; i<numCurrentDivisions; i++)
            {
            if (modelSettings[i].gibbsGamma == YES)
                curLnL[chn] += GibbsSampleGamma (chn, i, seed);
            }
        MrBayesPrint ("%s      Chain %d -- %.6lf -- %.6lf\n", spacer, (chn % chainParams.numChains) + 1, curLnL[chn], curLnPr[chn]);
        }
    MrBayesPrint("\n");

#   if defined (MPI_ENABLED)
    if (num_procs > 2)
        MrBayesPrint ("%s   There are %d more chains on other processor(s)\n\n", spacer, numGlobalChains - numLocalChains);
    else if (num_procs==2)
        MrBayesPrint ("%s   There are %d more chains on the other processor\n\n", spacer, numGlobalChains - numLocalChains);
#   endif

    /* All steps are assumed to have the same length. */
    if (chainParams.isSS == YES)
        {
        numGenInStepSS = (chainParams.numGen - chainParams.burninSS*chainParams.sampleFreq)/ chainParams.numStepsSS;
        numGenInStepSS = chainParams.sampleFreq*(numGenInStepSS/chainParams.sampleFreq); /*make muliple of chainParams.sampleFreq*/
        numGenOld = chainParams.numGen;
        chainParams.numGen = (chainParams.burninSS * chainParams.sampleFreq + chainParams.numStepsSS*numGenInStepSS) ; 
        if (stepRelativeBurninSS==YES)
            numGenInStepBurninSS = ((int)(numGenInStepSS*chainParams.burninFraction / chainParams.sampleFreq))*chainParams.sampleFreq;
        else
            numGenInStepBurninSS = chainParams.chainBurnIn * chainParams.sampleFreq;
        MrBayesPrint ("\n");
        MrBayesPrint ("%s   Starting stepping-stone sampling to estimate marginal likelihood.         \n", spacer);
        MrBayesPrint ("%s   %d steps will be used with %d generations (%d samples) within each step.  \n", spacer, chainParams.numStepsSS, numGenInStepSS, numGenInStepSS/chainParams.sampleFreq);
        MrBayesPrint ("%s   Total of %d generations (%d samples) will be collected while first        \n", spacer, chainParams.numGen, chainParams.numGen/chainParams.sampleFreq);
        MrBayesPrint ("%s   %d generations (%d samples) will be discarded as initial burnin.          \n", spacer, chainParams.burninSS*chainParams.sampleFreq, chainParams.burninSS);
        MrBayesPrint ("%s   Additionally at the begining of each step %d generations (%d samples)     \n", spacer, numGenInStepBurninSS, numGenInStepBurninSS/chainParams.sampleFreq);
        MrBayesPrint ("%s   will be discarded as burnin.  \n", spacer);
        if (chainParams.startFromPriorSS==YES)
            MrBayesPrint ("%s   Sampling from prior to posterior, i.e. first step samples from prior.   \n", spacer);
        else
            {
            MrBayesPrint ("%s   Sampling from posterior to prior, i.e. first step samples from close to \n", spacer);
            MrBayesPrint ("%s   posterior.                                                              \n", spacer);
            }
        if (numGenOld != chainParams.numGen)
            {
            MrBayesPrint ("%s   NOTE: Number of generation of each step is reduced to the closest multi-\n", spacer);
            MrBayesPrint ("%s   ple of sampling frequency. That is why, in total it will be taken %d    \n", spacer, chainParams.numGen);
            MrBayesPrint ("%s   generations instead of requested %d.                                    \n", spacer, numGenOld);
            }
        MrBayesPrint ("\n");
        if ((numGenInStepSS-numGenInStepBurninSS)/chainParams.sampleFreq < 1)
            {
            MrBayesPrint ("%s   There is less then one sample in each step of stepping-stone sampling.  \n", spacer);
            MrBayesPrint ("%s   Please adjust burnin, nuber of generations, sampling frequency or       \n", spacer);
            MrBayesPrint ("%s   numnber of step in order to allow at least one sample per step.         \n", spacer);
            return ERROR; /*All MPI run will return here since all of them have the same values*/
            }
        if (numPreviousGen==0 || numPreviousGen < chainParams.burninSS * chainParams.sampleFreq)
            {
            lastStepEndSS = chainParams.burninSS * chainParams.sampleFreq;
            stepIndexSS = chainParams.numStepsSS-1;
            if (numPreviousGen != 0)
                removeTo=(numPreviousGen/chainParams.sampleFreq)+1;
            if (chainParams.startFromPriorSS==YES)
                {
                // powerSS = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-1-stepIndexSS)/(MrBFlt)chainParams.numStepsSS);
                powerSS = 0.0;
                stepLengthSS = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-stepIndexSS)/(MrBFlt)chainParams.numStepsSS)-powerSS;
                }
            else
                {
                powerSS = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)stepIndexSS/(MrBFlt)chainParams.numStepsSS);
                stepLengthSS = 1.0-powerSS;
                }
            samplesCountSS=0;
            }
        else
            {
            stepIndexSS     = (numPreviousGen-chainParams.burninSS * chainParams.sampleFreq)/numGenInStepSS; /* for now it holds number of steps we fully complited*/
            lastStepEndSS   = chainParams.burninSS * chainParams.sampleFreq + stepIndexSS*numGenInStepSS;
            removeTo        = chainParams.burninSS + (stepIndexSS*numGenInStepSS+numGenInStepBurninSS)/chainParams.sampleFreq + 1;
            if (numPreviousGen < (removeTo-1)*chainParams.sampleFreq)
                removeTo=numPreviousGen/chainParams.sampleFreq+1;
            stepIndexSS     = chainParams.numStepsSS-1-stepIndexSS;
            if (chainParams.startFromPriorSS==YES)
                {
                powerSS = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-1-stepIndexSS)/(MrBFlt)chainParams.numStepsSS);
                stepLengthSS = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-stepIndexSS)/(MrBFlt)chainParams.numStepsSS)-powerSS;
                }
            else
                {
                powerSS         = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)stepIndexSS/(MrBFlt)chainParams.numStepsSS);
                stepLengthSS    = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)(stepIndexSS+1)/(MrBFlt)chainParams.numStepsSS)-powerSS;
                }
#   ifdef SAMPLE_ALL_SS
            samplesCountSS  = (numPreviousGen-lastStepEndSS-numGenInStepBurninSS);
#   else
            samplesCountSS  = (numPreviousGen-lastStepEndSS-numGenInStepBurninSS)/chainParams.sampleFreq;
#   endif
            if (samplesCountSS < 0)
                samplesCountSS=0;

            MrBayesPrint("%s   Continue sampling step %d out of %d steps...\n",spacer, chainParams.numStepsSS-stepIndexSS, chainParams.numStepsSS);
            /*marginalLnLSS will be red from file and destributed to other MPI_proc later. stepScalerSS, stepAcumulatorSS are lready red and if (samplesCountSS!=0) they will be redestributed. */
            }

        if (samplesCountSS == 0) /* in appended case it also can happen */
            {
            for (run=0; run<chainParams.numRuns; run++)
                {
                marginalLnLSS[run] = 0.0;
                stepScalerSS[run] = 0.0;
                stepAcumulatorSS[run] = 0.0;
                }
 
            for (chn=0; chn<numLocalChains; chn++)
                {
                if (chainId[chn] % chainParams.numChains == 0)
                    {
                    run = chainId[chn] / chainParams.numChains;
                    stepScalerSS[run] = curLnL[chn]*stepLengthSS;
                    }
                }      
            }
        }

    /* Append to previous analysis if this is requested, otherwise just open new print files */
    if (chainParams.append == YES)
        {
#   if defined (MPI_ENABLED)
    if (proc_id == 0) {
#   endif
        /* We get the number of samples in i */
        if (ReusePreviousResults(&i, chainParams.numStepsSS-stepIndexSS-1) == ERROR || i < 2)
            nErrors++;
        else if (chainParams.numRuns > 1 && chainParams.mcmcDiagn == YES)  /* we potentially need to add tree samples for conv diagn */
            {
            /* Add tree samples to partition counters */
            if (chainParams.relativeBurnin == YES)
                {
                if (numPreviousGen/(i-1) != chainParams.sampleFreq)
                    {
                    MrBayesPrint ("%s   1. Use the same sampling frequency as in the previous run to use relative burnin.\n", spacer);
                    MrBayesPrint ("%s   2. Check (and modify) the number in [generation: number] at line 3 of the .ckp file\n", spacer);
                    MrBayesPrint ("%s      to match the previous number of generations in all the .p and .t files. This may\n", spacer);
                    MrBayesPrint ("%s      happen if checkfreq was smaller than samplefreq.\n", spacer);
                    MrBayesPrint ("%s   3. Rarely, delete the last sample/line in the .p and .t files to achieve 2. above.\n", spacer);
                    MrBayesPrint ("%s      This may happen if ngen was not divisible by samplefreq.\n", spacer);
                    nErrors++;
                    }
                if (chainParams.isSS == NO)
                    {
                    if (noWarn == YES)
                        {
                        /* We definitely know the final number of generations */
                        j = (chainParams.numGen/chainParams.sampleFreq)+1;
                        j = (int) (j*chainParams.burninFraction);
                        }
                    else /* User may extend chain so save all trees if saving trees */
                        j = i;
                    if (j < i)
                        {
                        if (AddTreeSamples(1,j,chainParams.saveTrees) == ERROR) nErrors++;
                        if (AddTreeSamples(j+1,i,NO) == ERROR) nErrors++;
                        /* Since we never need to remove trees from partition counter after total burnin we put NO in the last argument */
                        }
                    else
                        {
                        if (AddTreeSamples(1,i,chainParams.saveTrees) == ERROR) nErrors++;
                        }
                    }
                else
                    {
                    if (SetFilePositions(removeTo) == ERROR) nErrors++;
                    if (AddTreeSamples(removeTo+1,i,chainParams.saveTrees) == ERROR) nErrors++;
                    }
                }
            else if (chainParams.chainBurnIn < i)
                {
                if (AddTreeSamples(chainParams.chainBurnIn,i-1,chainParams.saveTrees) == ERROR) nErrors++;
                }
            }
        if (nErrors == 0)
            {
            if (chainParams.isSS == NO && chainParams.mcmcDiagn == YES && chainParams.numRuns > 1)
                {
                MrBayesPrint ("\n");
                if (chainParams.relativeBurnin == YES)
                    MrBayesPrint ("%s   Using a relative burnin of %.1f %% for diagnostics\n", spacer, 100.0*chainParams.burninFraction);
                else
                    MrBayesPrint ("%s   Using an absolute burnin of %d samples for diagnostics\n", spacer, chainParams.chainBurnIn);
                }
            MrBayesPrint ("\n");
            MrBayesPrint ("%s   Chain results (continued from previous run; %d generations requested):\n\n", spacer, chainParams.numGen);
            }
#   if defined (MPI_ENABLED)
        }
#   endif

        if (chainParams.autotune == YES)
            {
            for (i=0; i<numLocalChains; i++)
                {
                for (j=0; j<numUsedMoves; j++)
                    {
                    if (j==0)
                        f = usedMoves[j]->cumProposalProb[chainId[i]];
                    else
                        f = usedMoves[j]->cumProposalProb[chainId[i]] - usedMoves[j-1]->cumProposalProb[chainId[i]];
                    if (usedMoves[j]->targetRate[chainId[i]] > 0.0)
                        {
                        /* autotuned move; we assume it was perfectly tuned and start tuning from there at appropriate rate */
                        usedMoves[j]->nBatches[chainId[i]] = (int) (f*numPreviousGen/chainParams.tuneFreq);
                        usedMoves[j]->lastAcceptanceRate[chainId[i]] = usedMoves[j]->targetRate[chainId[i]];
                        }
                    else
                        {
                        /* not autotuned move; no previous batches will result in acceptance rate calculated from new samples */
                        usedMoves[j]->nBatches[chainId[i]] = 0;
                        usedMoves[j]->lastAcceptanceRate[chainId[i]] = 0.0;
                        }
                    }
                }
            }
#   if defined (MPI_ENABLED)
        if (chainParams.isSS == YES)
            {
            MPI_Bcast (marginalLnLSS, chainParams.numRuns, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (samplesCountSS != 0)
                {
                MPI_Bcast (stepScalerSS, chainParams.numRuns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Bcast (stepAcumulatorSS, chainParams.numRuns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                }
            /*Set to zero all runs that we are not responsable for*/
            for (run=0; run<chainParams.numRuns; run++)
                {
                for (chn=0; chn<numLocalChains; chn++)
                    {
                    if (chainId[chn] % chainParams.numChains == 0 && run == chainId[chn] / chainParams.numChains)
                        break;
                    }
                if (chn<numLocalChains)
                    continue;

                marginalLnLSS[run] = 0.0;
                stepScalerSS[run] = 0.0;
                stepAcumulatorSS[run] = 0.0;    
                }
            }

        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s    Error appending to previous run\n", spacer);
            return ERROR;
            }            
#   else
        if (nErrors == 1)
            {
            MrBayesPrint ("%s    Error appending to previous run\n", spacer);
            return ERROR;
            }
#   endif
        }
    else
        {
        if (PreparePrintFiles() == ERROR)
            nErrors++;
        }

#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   Error preparing print files on at least one processor\n", spacer);
        CloseMBPrintFiles();
        return ERROR;
        }
#   else
    if (nErrors > 0)
        {
        MrBayesPrint ("%s   Error preparing print files\n", spacer);
        CloseMBPrintFiles();
        return ERROR;
        }
#   endif

    if (chainParams.relativeBurnin == NO)
        lastDiagnostics = chainParams.chainBurnIn;
    else
        lastDiagnostics = 0;
    stopChain = NO;

    for (i=0; i<chainParams.numRuns; i++)
        maxLnL0[i] = -100000000.0;

    startingT=time(0);
    CPUTime = 0.0;
    previousCPUTime = clock();

    /* print headers and starting states */
    if (numPreviousGen==0)
        {
        /* make sure we print headers */
        PrintToScreen(0, 0, time(0), startingT);
        if (PrintStatesToFiles (0) == ERROR)
            {
            MrBayesPrint("%s   Error in printing headers to files\n");
#   if defined (MPI_ENABLED)
            nErrors++;
#   else
            return ERROR;
#   endif
            }
#   if defined (MPI_ENABLED)
        MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (sumErrors > 0)
            {
            MrBayesPrint ("%s   Aborting run.\n");
            return ERROR;
            }
#   endif

        if (chainParams.mcmcDiagn == YES)
            {
            if (PrintMCMCDiagnosticsToFile (0) == ERROR)
                {
                MrBayesPrint ("%s   Problem printing mcmc diagnostics headers to file\n", spacer);
#   if defined (MPI_ENABLED)
                nErrors++;
#   else
                return (ERROR);
#   endif
                }
#   if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Aborting run.\n");
                return ERROR;
                }
#   endif
        
            if (chainParams.isSS == YES && chainParams.burninSS == 0 && chainParams.numRuns > 1)
                {
                /* Remove first sample (generation 0) from diagnostics */
                removeTo=1;
                if (RemoveTreeSamples (1,1) == ERROR)
                    {
                     MrBayesPrint("%s   Problem removing tree samples\n");
#   if defined (MPI_ENABLED)
                       nErrors++;
#   else
                      return ERROR;
#   endif
                    }
#   if defined (MPI_ENABLED)
                MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                if (sumErrors > 0)
                    {
                    MrBayesPrint ("%s   Aborting run.\n");
                    return ERROR;
                    }
#   endif
                }
            }
        if (chainParams.isSS == YES)
            {
            if (chainParams.burninSS == 0)
                MrBayesPrint("%s   Sampling step 1 out of %d steps...\n\n",spacer, chainParams.numStepsSS);

            /*Printing SS header*/
            MrBayesPrintf (fpSS, "[LEGEND: The file contains statistics on the Steppingstone Sampling.]\n");
            MrBayesPrintf (fpSS, "[ID: %s]\n", stamp);
            MrBayesPrintf (fpSS, "[   Step                --  Index of the step ]\n");
            MrBayesPrintf (fpSS, "[   Power               --  At each step we sample from the distribution with density (Likelihood^Power)*Prior ]\n");
            MrBayesPrintf (fpSS, "[   runX                --  Contribution to the marginal log likelihood of run X, i.e. marginal log likelihood for run X is the sum across all steps in column runX.   ]\n");
            if (chainParams.diagnStat == AVGSTDDEV)
                MrBayesPrintf (fpSS, "[   aSplitX             --  Average standard deviation of split frequencies of tree X. -2.0 is printed if no diagnostics was requested. -1.0 is printed if there were no splits with frequency above minimum.]\n");
            else
                MrBayesPrintf (fpSS, "[   mSplitX             --  Maximal standard deviation of split frequencies of tree X. -2.0 is printed if no diagnostics was requested. -1.0 is printed if there were no splits with frequency above minimum.]\n");
            MrBayesPrintf (fpSS, "Step\tPower");
            for (j=0; j<chainParams.numRuns ; j++)
                MrBayesPrintf (fpSS, "\trun%d", j+1);
            if (chainParams.diagnStat == AVGSTDDEV)
                {
                for (j=0; j<numTopologies; j++)
                    MrBayesPrintf (fpSS, "\taSplit%d", j);
                }
            else
                {
                for (j=0; j<numTopologies; j++)
                    MrBayesPrintf (fpSS, "\tmSplit%d", j);
                }
            MrBayesPrintf (fpSS, "\n");
            }
        }

    for (n=numPreviousGen+1; n<=chainParams.numGen; n++) /* begin run chain */
        {
        currentCPUTime = clock();
        if (currentCPUTime - previousCPUTime > 10 * CLOCKS_PER_SEC)
            {
            CPUTime += (currentCPUTime - previousCPUTime) / (MrBFlt) (CLOCKS_PER_SEC);
            previousCPUTime = currentCPUTime;
            }

        /*! requestAbortRun is set by the signal handler when it receives a CTRL-C (serial version only) */
        if (requestAbortRun == YES && ConfirmAbortRun() == 1)
            return ABORT;

        // RandLong oldSeed = *seed;  /* record the old seed for debugging */
        for (chn=0; chn<numLocalChains; chn++)
            {
            /* Do Gibbs resampling of rate categories for current state if time to do so */
            for (i=0; i<numCurrentDivisions; i++)
                {
                if (modelSettings[i].gibbsGamma == YES && n % modelSettings[i].gibbsFreq == 0)
                    curLnL[chn] += GibbsSampleGamma (chn, i, seed);
                }

            /* First copy everything from current state of chain to new state.   */
            /* The global variable state[chain] gives state.                     */

            /* copy all touched trees and reset update flags                     */
            CopyTrees (chn);

            /* copy all model parameters */
            CopyParams (chn);

            /* shift the state of the chain to the new state */
            /* all calculations will be done on this state   */
            state[chn] ^= 1;  /* XORing with 1 switches between 0 and 1 */

            /* decide which move to make */
            whichMove = PickProposal(seed, chainId[chn]);
            theMove = usedMoves[whichMove];
#   if defined SHOW_MOVE
            printf ("Making move '%s'\n", theMove->name);
#   endif

#   if defined (BEST_MPI_ENABLED)
            bestCycleGen = n % (numNonTreeMoves + numTreeMoves + numBestMoves);
            if (bestCycleGen < numNonTreeMoves)
                PickNonTreeProposal(seed, chainId[chn]);
            else if (bestCycleGen < numNonTreeMoves + numTreeMoves)
                PickTreeProposal(seed, chainId[chn]);
            else
                PickBestProposal(chainId[chn]);
#   endif

            /* set prior and proposal ratios */
            lnProposalRatio = 0.0;
            lnPriorRatio = 0.0;

            /* reset abort move flag */
            abortMove = NO;
            
            /* Touch the relevant partitions       */
            /* as a service to the move functions. */
            for (i=0; i<theMove->parm->nRelParts; i++)
                modelSettings[theMove->parm->relParts[i]].upDateCl = YES;

#   ifndef NDEBUG
            if (IsTreeConsistent(theMove->parm, chn, state[chn]) != YES)
                {
                printf ("IsTreeConsistent failed before a move!\n");
                return ERROR;
                }
#   endif
#   if defined (DEBUG_CONSTRAINTS)
            if (theMove->parm->paramType == P_TOPOLOGY && DoesTreeSatisfyConstraints(GetTree (theMove->parm, chn, state[chn]))!=YES)
                {
                printf ("DEBUG ERROR: DoesTreeSatisfyConstraints failed before a move\n");
                return ERROR;
                }
#   endif
            /* make move */
            if ((theMove->moveFxn)(theMove->parm, chn, seed, &lnPriorRatio, &lnProposalRatio, theMove->tuningParam[chainId[chn]]) == ERROR)
                {
                printf ("%s   Error in move %s\n", spacer, theMove->name);
#   if defined (MPI_ENABLED)
                nErrors++;
#   else
                return ERROR;
#   endif
                }

            if (theMove->parm->paramType == P_TOPOLOGY && DoesTreeSatisfyConstraints(GetTree (theMove->parm, chn, state[chn])) != YES)
                {
#   if defined (DEBUG_CONSTRAINTS)
                if (DoesTreeSatisfyConstraints(GetTree (theMove->parm, chn, state[chn]))==ABORT)
                    {
                    printf ("DEBUG ERROR: DoesTreeSatisfyConstraints failed after move '%s'\n", theMove->name);
                    }
#   endif
                abortMove = YES;
                }

            /* abortMove is set to YES if the calculation fails because the likelihood is too small */
            if (abortMove == NO)
                lnLike = LogLike(chn);

            /* calculate acceptance probability */
            if (abortMove == NO)
                {
                lnLikelihoodRatio = lnLike - curLnL[chn];
                lnPrior = curLnPr[chn] + lnPriorRatio;

#   ifndef NDEBUG
                /* We check various aspects of calculations in debug version of code */
                if (IsTreeConsistent(theMove->parm, chn, state[chn]) != YES)
                    {
                    printf ("DEBUG ERROR: IsTreeConsistent failed after move '%s'\n", theMove->name);
                    return ERROR;
                    }
                if (lnPriorRatio != lnPriorRatio)
                    {
                    printf ("DEBUG ERROR: Log prior ratio nan after move '%s'\n", theMove->name);
                    // printf ("Seed: %ld\n", oldSeed);  state[chn] ^= 1;  PrintCheckPoint (n);
                    return ERROR;
                    }
                if (fabs((lnPrior-LogPrior(chn))/lnPrior) > 0.0001)
                    {
                    printf ("DEBUG ERROR: Log prior incorrect after move '%s' :%e :%e\n", theMove->name,lnPrior,LogPrior(chn));
                    // printf ("Seed: %ld\n", oldSeed);  state[chn] ^= 1;  PrintCheckPoint (n);
                    return ERROR;
                    }
                if (lnProposalRatio != lnProposalRatio)
                    {
                    printf ("DEBUG ERROR: Log proposal ratio nan after move '%s'\n", theMove->name);
                    // printf ("Seed: %ld\n", oldSeed);  state[chn] ^= 1;  PrintCheckPoint (n);
                    return ERROR;
                    }
                if (lnLike != lnLike)
                    {
                    printf ("DEBUG ERROR: Log likelihood nan after move '%s'\n", theMove->name);
                    // printf ("Seed: %ld\n", oldSeed);  state[chn] ^= 1;  PrintCheckPoint (n);
                    return ERROR;
                    }
#       if defined (DEBUG_LNLIKELIHOOD) /* slow */
                ResetFlips(chn); /* needed to return flags so they point to old state */
                TouchEverything(chn);
                if (fabs((lnLike-LogLike(chn))/lnLike) > 0.0001)
                    {
                    printf ("DEBUG ERROR: Log likelihood incorrect after move '%s'\n", theMove->name);
                    return ERROR;
                    }
#       endif
                if (theMove->parm->paramType == P_TOPOLOGY && GetTree (theMove->parm, chn, state[chn])->isClock == YES &&
                    IsClockSatisfied (GetTree (theMove->parm, chn, state[chn]),0.001) == NO)
                    {
                    printf ("%s   Branch lengths of the tree do not satisfy the requirements of a clock tree.\n", spacer);
                    ShowNodes(GetTree (theMove->parm, chn, state[chn])->root,0,YES);
                    return (ERROR);
                    }
#   endif

                /* heat */
                lnLikelihoodRatio *= Temperature (chainId[chn]);
                lnPriorRatio      *= Temperature (chainId[chn]);

                if (chainParams.isSS == YES)
                    lnLikelihoodRatio *= powerSS;

                /* calculate the acceptance probability */
                if (lnLikelihoodRatio + lnPriorRatio + lnProposalRatio < -100.0)
                    r = 0.0;
                else if (lnLikelihoodRatio + lnPriorRatio + lnProposalRatio > 0.0)
                    r = 1.0;
                else
                    r = exp(lnLikelihoodRatio + lnPriorRatio + lnProposalRatio);
                }

            /* decide to accept or reject the move */
            acceptMove = NO;
            i = chainId[chn];
            theMove->nTried[i]++;
            theMove->nTotTried[i]++;
            if (abortMove == NO && RandomNumber(seed) < r)
                {
                acceptMove = YES;
                theMove->nAccepted[i]++;
                theMove->nTotAccepted[i]++;
                }

            /* update the chain */
            if (acceptMove == NO)
                {
                /* the new state did not work out so shift chain back */
                if (abortMove == NO)
                    ResetFlips(chn);
                state[chn] ^= 1;
#   if defined (BEAGLE_ENABLED)
                if (recalcScalers == YES)
                    {
                    recalculateScalers(chn);
                    recalcScalers = NO;
                    }
#   endif
                }
            else
                {
                /* if the move is accepted then let the chain stay in the new state */
                /* store the likelihood and prior of the chain */
                curLnL[chn] = lnLike;
                curLnPr[chn] = lnPrior;
                }

            /* check if time to autotune */
            if (theMove->nTried[i] >= chainParams.tuneFreq)
                {
                theMove->lastAcceptanceRate[i] = (MrBFlt) theMove->nAccepted[i] / (MrBFlt) theMove->nTried[i];
                theMove->nTried[i] = 0;
                theMove->nAccepted[i] = 0;
                theMove->nBatches[i]++;                                     /* we only autotune at most 10000 times */
                if (chainParams.autotune == YES && theMove->moveType->Autotune != NULL && theMove->nBatches[i] < MAXTUNINGPARAM)
                    {
                    theMove->moveType->Autotune(theMove->lastAcceptanceRate[i],
                                                theMove->targetRate[i],
                                                theMove->nBatches[i],
                                                &theMove->tuningParam[i][0],
                                                theMove->moveType->minimum[0],
                                                theMove->moveType->maximum[0]);
                    }
                }

            /* ShowValuesForChain (chn); */

            if (curLnL[chn] > maxLnL0[chainId[chn]])
                maxLnL0[chainId[chn]] = curLnL[chn];

            }

        /* attempt swap(s) Non-blocking for MPI if no swap with external process. */
        if (chainParams.numChains > 1 && n % chainParams.swapFreq == 0)
            {
            for (i = 0; i<chainParams.numRuns; i++)
                {
                for (j = 0; j<chainParams.numSwaps; j++)
                    GetSwappers (&swapA, &swapB, i);
                if (AttemptSwap (swapA, swapB, seed) == ERROR)
                    {
                    MrBayesPrint ("%s   Unsuccessful swap of states\n", spacer);
#   if defined (MPI_ENABLED)
                    nErrors++;
#   else
                    return ERROR;
#   endif
                    }
                }
            }

        /* print information to screen. Non-blocking for MPI */
        if (n % chainParams.printFreq == 0)
            {
            PrintToScreen(n, numPreviousGen, time(0), startingT);
#   if defined (TIMING_ANALIZ)
            MrBayesPrint ("%s   Time elapsed:%f CondlikeDownTime:%f CondLikeRoot:%f Likelihood:%f ScalersTime:%f ScalersRemove:%f\n", spacer,
                          CPUTime, CPUCondLikeDown/(MrBFlt)CLOCKS_PER_SEC, CPUCondLikeRoot/(MrBFlt)CLOCKS_PER_SEC, CPULilklihood/(MrBFlt)CLOCKS_PER_SEC,
                          CPUScalers/(MrBFlt)CLOCKS_PER_SEC, CPUScalersRemove/(MrBFlt)CLOCKS_PER_SEC);
#   endif
            }

        /* print information to files */
        /* this will also add tree samples to topological convergence diagnostic counters */
        if (n == chainParams.numGen || n % chainParams.sampleFreq == 0)
            {
#   if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Aborting run.\n");
                return ERROR;
                }
#   endif
            if (PrintStatesToFiles (n) == ERROR)
                {
                MrBayesPrint("%s   Error in printing states to files\n");
#   if defined (MPI_ENABLED)
                nErrors++;
#   else
                return ERROR;
#   endif
                }
#   if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Aborting run.\n");
                return ERROR;
                }
#   endif
            }

        /* print mcmc diagnostics. Blocking for MPI */
        if (chainParams.mcmcDiagn == YES && (n % chainParams.diagnFreq == 0
                                             || n == chainParams.numGen
                                             || (chainParams.isSS == YES && (n-lastStepEndSS) % numGenInStepSS == 0)))
            {
            if (chainParams.numRuns > 1 &&
                ((n > 0 && chainParams.relativeBurnin == YES && (chainParams.isSS == NO || (n > chainParams.burninSS * chainParams.sampleFreq && (n-lastStepEndSS) > numGenInStepBurninSS)))
                 || (n >= chainParams.chainBurnIn * chainParams.sampleFreq && chainParams.relativeBurnin == NO)))
                {
                /* we need some space for coming output */
                MrBayesPrint ("\n");
                /* remove tree samples if using burninpercentage */
                /* the following function returns immediately in MPI if proc_id != 0 */
                if (chainParams.relativeBurnin == YES && chainParams.isSS == NO)
                    {
                    removeFrom = removeTo;
                    removeTo = (int)(chainParams.burninFraction * (n/chainParams.sampleFreq+1)); /* (n/chainParams.sampleFreq+1) is the current number of samples */
                    if (removeFrom < removeTo)
                        {
                        if (RemoveTreeSamples (removeFrom+1, removeTo) == ERROR)
                            {
                            MrBayesPrint("%s   Problem removing tree samples\n");
#   if defined (MPI_ENABLED)
                            nErrors++;
#   else
                            return ERROR;
#   endif
                            }
                        }
#   if defined (MPI_ENABLED)
                    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                    if (sumErrors > 0)
                        {
                        MrBayesPrint ("%s   Aborting run.\n");
                        return ERROR;
                        }
#   endif
                    }

                lastDiagnostics = (n/chainParams.sampleFreq)+1; /* +1 because we always have start tree sampled*/
                if (chainParams.relativeBurnin == YES)
                    {
                    i = lastDiagnostics - removeTo;
                    }
                else
                    i = lastDiagnostics - chainParams.chainBurnIn;
#   if defined (MPI_ENABLED)
                if (proc_id == 0)
                    {
#   endif
                /* calculate statistics */
                CalcTopoConvDiagn (i);
                /* output statistics */
                if (numTopologies == 1)
                    {
                    f = -1.0;
                    if (chainParams.stat[0].numPartitions == 0)
                        {
                        MrBayesPrint ("%s   Average standard deviation of split frequencies: NA (no splits above min. frequency)\n", spacer);
                        }
                    else if (chainParams.diagnStat == AVGSTDDEV)
                        {
                        f = chainParams.stat[0].avgStdDev;
                        MrBayesPrint ("%s   Average standard deviation of split frequencies: %.6f\n", spacer, f);
                        }
                    else
                        {
                        f = chainParams.stat[0].max;
                        MrBayesPrint ("%s   Max standard deviation of split frequencies: %.6f\n", spacer, f);
                        }
                    if (chainParams.isSS == YES)
                        splitfreqSS[chainParams.numStepsSS-stepIndexSS-1] = f;
                    if (chainParams.stat[0].numPartitions > 0 && f <= chainParams.stopVal)
                        stopChain = YES;
                    if (n < chainParams.numGen - chainParams.printFreq && (chainParams.stopRule == NO || stopChain == NO))
                        MrBayesPrint ("\n");
                    }
                else
                    {
                    stopChain = YES;
                    for (i=0; i<numTopologies; i++)
                        {
                        f=-1.0;
                        if (chainParams.stat[i].numPartitions == 0)
                            {
                            if (strcmp(modelParams[0].topologyPr,"Speciestree") == 0)
                                {
                                if (i == numTopologies-1)
                                    MrBayesPrint ("%s   Average standard deviation of split frequencies for species tree (tree %d): NA (no splits above min. frequency)\n", spacer, i+1);
                                else
                                    MrBayesPrint ("%s   Average standard deviation of split frequencies for gene tree %d: NA (no splits above min. frequency)\n", spacer, i+1);
                                }
                            else
                                MrBayesPrint ("%s   Average standard deviation of split frequencies for topology %d: NA (no splits above min. frequency)\n", spacer, i+1);
                            }
                        else if (chainParams.diagnStat == AVGSTDDEV)
                            {
                            f = chainParams.stat[i].avgStdDev;
                            if (strcmp(modelParams[0].topologyPr,"Speciestree") == 0)
                                {
                                if (i == numTopologies-1)
                                    MrBayesPrint ("%s   Average standard deviation of split frequencies for species tree (tree %d): %.6f\n", spacer, i+1, f);
                                else
                                    MrBayesPrint ("%s   Average standard deviation of split frequencies for gene tree %d: %.6f\n", spacer, i+1, f);
                                }
                            else
                                MrBayesPrint ("%s   Average standard deviation of split frequencies for topology %d: %.6f\n", spacer, i+1, f);
                            }
                        else
                            {
                            f = chainParams.stat[i].max;
                            if (strcmp(modelParams[0].topologyPr,"Speciestree") == 0)
                                {
                                if (i == numTopologies-1)
                                    MrBayesPrint ("%s   Max standard deviation of split frequencies for species tree: %.6f\n", spacer, f);
                                else
                                    MrBayesPrint ("%s   Max standard deviation of split frequencies for gene tree %d: %.6f\n", spacer, i+1, f);
                                }
                            else
                                MrBayesPrint ("%s   Max standard deviation of split frequencies for topology %d: %.6f\n", spacer, i+1, f);
                            }
                        if (chainParams.isSS == YES)
                            splitfreqSS[i*chainParams.numStepsSS+chainParams.numStepsSS-stepIndexSS-1] = f;
                        if (chainParams.stat[i].numPartitions == 0 || f > chainParams.stopVal)
                            stopChain = NO;
                        }
                    if (n < chainParams.numGen - chainParams.printFreq && (chainParams.stopRule == NO || stopChain == NO))
                        MrBayesPrint ("\n");
                    }
                if (chainParams.allComps == YES)
                    PrintTopConvInfo ();
#   if defined (MPI_ENABLED)
                    }
                ierror = MPI_Bcast (&stopChain, 1, MPI_INT, 0, MPI_COMM_WORLD);
                if (ierror != MPI_SUCCESS)
                    {
                    MrBayesPrint ("%s   Problem broadcasting stop value\n", spacer);
                    nErrors++;
                    }
#   endif
                }

            /* part of the following function needs to be performed by all MPI processors. Blocking for MPI. */
            if (PrintMCMCDiagnosticsToFile (n) == ERROR)
                {
                MrBayesPrint ("%s   Problem printing mcmc diagnostics to file\n", spacer);
#   if defined (MPI_ENABLED)
                nErrors++;
#   else
                return (ERROR);
#   endif
                }
#   if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   Aborting run.\n");
                return ERROR;
                }
#   endif
            }

        /* check if time to break because stopVal reached */
        if (chainParams.isSS == NO && (chainParams.stopRule == YES && stopChain == YES))
            {
            MrBayesPrint ("\n%s   Analysis stopped because convergence diagnostic hit stop value.\n", spacer);
            break;
            }
            
        /* user may want to extend chain */
        if (chainParams.isSS == NO && (n == chainParams.numGen && autoClose == NO))
            {
            stoppingT1 = time(0);
            currentCPUTime = clock();
            CPUTime += (currentCPUTime - previousCPUTime) / (MrBFlt) CLOCKS_PER_SEC;
            previousCPUTime = currentCPUTime;
            chainParams.numGen += ExtendChainQuery ();
            stoppingT2 = time(0);
            startingT += (stoppingT2-stoppingT1);
            previousCPUTime = clock();
            /* timers should not be increased during the wait for a reply */
            }

        /* Do stepping sampling staf if needed */
        if (chainParams.isSS == YES && n >= chainParams.burninSS*chainParams.sampleFreq)
            {
#   ifndef SAMPLE_ALL_SS
            if ((n-lastStepEndSS) % chainParams.sampleFreq == 0)
#   endif
                {
                if (n > chainParams.burninSS*chainParams.sampleFreq && (n-lastStepEndSS > numGenInStepBurninSS))
                    { /* do sampling */
                    for (chn=0; chn<numLocalChains; chn++)
                        {
                        if (chainId[chn] % chainParams.numChains == 0)
                            {
                            run = chainId[chn] / chainParams.numChains;
                            if (curLnL[chn]*stepLengthSS > stepScalerSS[run] + 200.0)
                                {
                                // adjust scaler;
                                stepAcumulatorSS[run] /= exp (curLnL[chn]*stepLengthSS - 100.0 - stepScalerSS[run]); 
                                stepScalerSS[run]= curLnL[chn]*stepLengthSS - 100.0;
                                }
                            stepAcumulatorSS[run] += exp (curLnL[chn]*stepLengthSS - stepScalerSS[run]);
                            }
                        }
                    samplesCountSS++;
                    }
                }

            if ((n-lastStepEndSS) == numGenInStepBurninSS)
                {
                /* Remove all previouse samples from diagnostics */
                if (chainParams.mcmcDiagn == YES && chainParams.numRuns > 1)
                    {
                    removeFrom = removeTo;
                    removeTo = (int)(n/chainParams.sampleFreq); /* (n/chainParams.sampleFreq+1) is the current number of samples including 0 one*/
                    removeTo++;
                    if (removeFrom < removeTo)
                        {
                        if (RemoveTreeSamples (removeFrom+1, removeTo) == ERROR)
                            {
                            nErrors++;
                            }
                        ERROR_TEST2("Problem removing tree samples",return(ERROR),);
                        }
                    }               
                }

            if ((n-lastStepEndSS) % numGenInStepSS == 0)      /* prepare sample of next step */
                {
                assert (n-lastStepEndSS <= numGenInStepSS);
                lastStepEndSS=n;

                if (n > chainParams.burninSS*chainParams.sampleFreq)
                    {
                    /* dump to file current step contribution */
                    MrBayesPrintf (fpSS, "%3d\t%.4f", chainParams.numStepsSS-stepIndexSS, powerSS);
#   if defined (MPI_ENABLED)
                    for (j=0; j<chainParams.numRuns ; j++)
                        {
                        if (stepAcumulatorSS[j]==0)
                            r=0;
                        else
                            r = log (stepAcumulatorSS[j]/samplesCountSS) + stepScalerSS[j];
                        ierror = MPI_Reduce (&r,&sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                        if (ierror != MPI_SUCCESS)
                            {
                            MrBayesPrint ("%s   Problem with MPI_Reduce\n", spacer);
                            return ERROR;
                            }
                        if (proc_id == 0)
                            {
                            MrBayesPrintf (fpSS, "\t%.6f", sum);
                            }
                        }
#   else
                    for (j=0; j<chainParams.numRuns ; j++)
                        {
                        MrBayesPrintf (fpSS, "\t%.6f", log (stepAcumulatorSS[j]/samplesCountSS) + stepScalerSS[j]);
                        }
#   endif
                    if (chainParams.mcmcDiagn == YES && chainParams.numRuns > 1)
                        {
                        for (i=0; i<numTopologies; i++)
                            {
                            MrBayesPrintf (fpSS, "\t%.6f",splitfreqSS[i*chainParams.numStepsSS+chainParams.numStepsSS-stepIndexSS-1]);
                            }
                        }
                    else{
                        for (i=0; i<numTopologies; i++)
                            {
                            MrBayesPrintf (fpSS, "\t-2.0");
                            }

                        }

                    MrBayesPrintf (fpSS, "\n");
                    fflush (fpSS);

                    stepIndexSS--;

                    if (chainParams.startFromPriorSS==YES)
                        {
                        powerSS = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-1-stepIndexSS)/(MrBFlt)chainParams.numStepsSS);
                        stepLengthSS = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)(chainParams.numStepsSS-stepIndexSS)/(MrBFlt)chainParams.numStepsSS)-powerSS;
                        }
                    else
                        {
                        stepLengthSS = powerSS; 
                        powerSS = BetaQuantile (chainParams.alphaSS, 1.0, (MrBFlt)stepIndexSS/(MrBFlt)chainParams.numStepsSS);
                        stepLengthSS -= powerSS;
                        }
                    if (n != chainParams.numGen)
                        MrBayesPrint("%s   Sampling step %d out of %d steps...\n\n",spacer, chainParams.numStepsSS-stepIndexSS, chainParams.numStepsSS);
                    for (chn=0; chn<numLocalChains; chn++)
                        {
                        if (chainId[chn] % chainParams.numChains == 0)
                            {
                            run = chainId[chn] / chainParams.numChains;
                            if (samplesCountSS != 0)
                                marginalLnLSS[run] += log (stepAcumulatorSS[run]/samplesCountSS) + stepScalerSS[run];
                            stepScalerSS[run] = curLnL[chn]*stepLengthSS;
                            stepAcumulatorSS[run]=0.0;
                            }
                        }
                    samplesCountSS=0;
                    }
                else
                    {
                    MrBayesPrint("\n%s   Sampling step 1 out of %d steps...\n\n",spacer, chainParams.numStepsSS);
                    }

                if (chainParams.backupCheckSS !=0 && (chainParams.numStepsSS-stepIndexSS-1)% chainParams.backupCheckSS == 0)
                    {
                    /* print check-point file. Blocking for MPI */
                    ERROR_TEST2("Error before printing checkpoint",return(ERROR),);
                    if (PrintCheckPoint (n) == ERROR)
                        {
                        nErrors++;
                        }
                    ERROR_TEST2("Error in printing checkpoint",return(ERROR),);
                       
#   if defined (MPI_ENABLED)
                if (proc_id == 0)
                    {
#   endif
                    /* figure out check-point file names */
                    sprintf (ckpFileName, "%s%s.ckp", workingDir, chainParams.chainFileName);
                    sprintf (bkupFileName,"%s.ss%d", ckpFileName,chainParams.numStepsSS-stepIndexSS);
                    if (rename (ckpFileName, bkupFileName)!=0)
                        {
                        MrBayesPrint ("%s   Could not rename file %s to %s\n", spacer, ckpFileName, bkupFileName);
                        return ERROR;
                        }
                    strcpy (bkupFileName, ckpFileName);
                    strcat (bkupFileName, "~");
                    rename (bkupFileName,ckpFileName);
#   if defined (MPI_ENABLED)
                    } /* end of if (proc_id == 0)*/
#   endif
                    }
                }             
            }

        /* print check-point file. Blocking for MPI */
        if (chainParams.checkPoint == YES && (n % chainParams.checkFreq == 0))
            {
            ERROR_TEST2("Error before printing checkpoint",return(ERROR),);
            if (PrintCheckPoint (n) == ERROR)
                {
                nErrors++;
                }
            ERROR_TEST2("Error in printing checkpoint",return(ERROR),);
            }

        } /* end run chain */
    endingT = time(0);
    currentCPUTime = clock();
    CPUTime += (currentCPUTime - previousCPUTime) / (MrBFlt) CLOCKS_PER_SEC;

    CloseMBPrintFiles (); /* redundant because files closed in FreeChainMemory but kept here as a safeguard in case of future changes */

#   if defined (BEAGLE_ENABLED)
#       ifdef DEBUG_BEAGLE
    beagleScalingSchemeOld = beagleScalingScheme;
    beagleScalingScheme = MB_BEAGLE_SCALE_ALWAYS;
    for (chn=0; chn<numLocalChains; chn++)
        {
        if (chn % chainParams.numChains == 0)
            {
            if (chainParams.numRuns == 1)
                MrBayesPrint ("%s   Final log likelihoods and log prior probs (stored and calculated):\n", spacer);
            else
                MrBayesPrint ("%s   Final log likelihoods and log prior probs for run %d (stored and calculated):\n", spacer, chn / chainParams.numChains + 1);
            }
        TouchEverything (chn);
        for (i=0; i<numCurrentDivisions; i++)
            {
            if (modelSettings[i].gibbsGamma == YES)
                curLnL[chn] += GibbsSampleGamma (chn, i, seed);
            }
        MrBayesPrint ("%s      Chain %d -- %.6lf -- %.6lf\n", spacer, (chn % chainParams.numChains) + 1, curLnL[chn], curLnPr[chn]);
        MrBayesPrint ("%s      Chain %d -- %.6lf -- %.6lf\n", spacer, (chn % chainParams.numChains) + 1, LogLike(chn), LogPrior(chn));
        }
    beagleScalingScheme = beagleScalingSchemeOld;
#       endif
#   endif

    /* Make sure current state is reset and values copied back to state 0.
       Note that this can be tricky for Metropolis-coupled chains because
       the chain ids may necessitate some swapping of values among chains. */
    ResetChainIds ();

    MrBayesPrint ("\n");
    if (difftime (endingT, startingT) > 3600.0)
        MrBayesPrint ("%s   Analysis completed in %d hours %d mins %d seconds\n", spacer, 
              (int) (difftime(endingT, startingT)/3600), 
              ((int) (difftime(endingT, startingT))%3600)/60, 
              (int) (difftime(endingT, startingT))%60);
    else if (difftime (endingT, startingT) > 60.0)
        MrBayesPrint ("%s   Analysis completed in %d mins %d seconds\n", spacer, 
              (int) (difftime(endingT, startingT)/60), 
              (int) (difftime(endingT, startingT))%60);
    else if (difftime (endingT, startingT) > 2.0)
        MrBayesPrint ("%s   Analysis completed in %.0f seconds\n", spacer, 
              difftime(endingT, startingT));
    else if (difftime (endingT, startingT) >= 1.0)
        MrBayesPrint ("%s   Analysis completed in 1 second\n", spacer);
    else
        MrBayesPrint ("%s   Analysis completed in less than 1 second\n", spacer);

#   if defined (MPI_ENABLED)
    MrBayesPrint ("%s   Analysis used %1.2f seconds of CPU time on processor 0\n", spacer, (MrBFlt) CPUTime);
#   else
    MrBayesPrint ("%s   Analysis used %1.2f seconds of CPU time\n", spacer, (MrBFlt) CPUTime);
#   endif

#   if defined (MPI_ENABLED)
    /* find the best likelihoods across all of the processors */
    ierror = MPI_Barrier (MPI_COMM_WORLD);
    if (ierror != MPI_SUCCESS)
        {
        MrBayesPrint ("%s   Problem at chain barrier\n", spacer);
        return (ERROR);
        }
    for (j=0; j<numGlobalChains; j++)
        {
        best = maxLnL0[j];
        for (i=1; i<num_procs; i++)
            {
            if (proc_id == 0)
                {
                ierror = MPI_Recv (&maxLnL0[j], 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
                if (ierror != MPI_SUCCESS)
                    { 
                    MrBayesPrint ("%s   Problem with MPI_Recv", spacer);
                    return ERROR;
                    }
                if (maxLnL0[j] > best)
                    best = maxLnL0[j];
                }
            else if (proc_id == i)
                {
                ierror = MPI_Send (&maxLnL0[j], 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
                if (ierror != MPI_SUCCESS)
                    {
                    MrBayesPrint ("%s   Problem with MPI_Send\n", spacer);
                    return ERROR;
                    }
                }
            }
        maxLnL0[j] = best;
        }

    /* Collecting  marginal log likelihoods if SS is used */
    if (chainParams.isSS == YES)
        {
        for (j=0; j<chainParams.numRuns ; j++)
            {
            ierror = MPI_Reduce (&marginalLnLSS[j],&r, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (ierror != MPI_SUCCESS)
                {
                MrBayesPrint ("%s   Problem with MPI_Send\n", spacer);
                return ERROR;
                }
            if (proc_id == 0)
                {
                marginalLnLSS[j]=r;
                }
            }
        }
#   endif

    if (chainParams.numRuns == 1)
        {
        if (chainParams.numChains == 1)
            MrBayesPrint ("%s   Log likelihood of best state was %1.2lf\n", spacer, maxLnL0[0]);
        else
            MrBayesPrint ("%s   Log likelihood of best state for \"cold\" chain was %1.2lf\n", spacer, maxLnL0[0]);
        }
    else
        {
        for (j=0; j<chainParams.numRuns*chainParams.numChains; j++)
            {
            if (j % chainParams.numChains == 0)
                {
                if (chainParams.numChains == 1)
                    MrBayesPrint ("%s   Likelihood of best state for run %d was %1.2lf\n", spacer, j/chainParams.numChains+1, maxLnL0[j/chainParams.numChains]);
                else
                    MrBayesPrint ("%s   Likelihood of best state for \"cold\" chain of run %d was %1.2lf\n", spacer, j/chainParams.numChains+1, maxLnL0[j/chainParams.numChains]);
                }
            }
        }

    if (chainParams.isSS == YES)
        {
        MrBayesPrint ("\n");
        MrBayesPrint ("%s   Marginal likelihood (in natural log units) estimated using stepping-stone sampling based on\n", spacer);
        MrBayesPrint ("%s   %d steps with %d generations (%d samples) within each step. \n\n", spacer, chainParams.numStepsSS, numGenInStepSS, numGenInStepSS/chainParams.sampleFreq);
        MrBayesPrint ("%s       Run   Marginal likelihood (ln)\n",spacer);
        MrBayesPrint ("%s       ------------------------------\n",spacer);
        for (j=0; j<chainParams.numRuns; j++)
            {
            MrBayesPrint ("%s       %3d    %9.2f   \n", spacer, j+1, marginalLnLSS[j]);
            }
        MrBayesPrint ("%s       ------------------------------\n",spacer);
        if (chainParams.numRuns > 1)
            {
            MeanVarianceLog(marginalLnLSS,chainParams.numRuns,&meanSS,&varSS,NULL);
            MrBayesPrint ("%s       Mean:  %9.2f\n\n",spacer,meanSS);
            //MrBayesPrint ("%s       Mean:  %9.2lf  Scaled variance: %.2f of Marginal log likelihood estimates among runs.\n",spacer,meanSS,varSS-2*meanSS);
            //MrBayesPrint ("%s       Note: Scaled variance is given in log units and calculated as \"variance/mean^2\"\n",spacer);     
            }
        MrBayesPrint ("%s   More statistics on stepping-stone sampling is dumped to %s.ss file.\n", spacer, chainParams.chainFileName);

        if (chainParams.mcmcDiagn == YES)
            {
            if ((tempX = (MrBFlt *) SafeCalloc (chainParams.numStepsSS, sizeof(MrBFlt))) == NULL)
                {
                nErrors++;
                }
            ERROR_TEST2("Problem allocating memory", return(ERROR),);

            for (i=0; i<chainParams.numStepsSS; i++)
                {
                tempX[i]=i+1;
                }
            MrBayesPrint ("\n");

            if (numTopologies > 1)
                {
                if (chainParams.diagnStat == AVGSTDDEV)
                    MrBayesPrint ("   Plots of average standard deviation of split frequencies across steps for different topology.");
                else
                    MrBayesPrint ("   Plots of max standard deviation of split frequencies across steps for different topology.");
                }
            else
                {
                if (chainParams.diagnStat == AVGSTDDEV)
                    MrBayesPrint ("   Plot of average standard deviation of split frequencies across steps.");
                else
                    MrBayesPrint ("   Plot of max standard deviation of split frequencies across steps.");
                }
            MrBayesPrint ("\n");
            MrBayesPrint ("   Points at -1.0 (y-axis) indicate that there were no splits\n");
            MrBayesPrint ("   above minimum frequency for corresponding step.");
            if (numTopologies > 1)
                {
                for (i=0; i<numTopologies; i++)
                    {
                    MrBayesPrint ("%s   Topology %d.\n", spacer, i+1);
                    PrintPlot (tempX, splitfreqSS+i*chainParams.numStepsSS, chainParams.numStepsSS);
                    }
                }
            else
                PrintPlot (tempX, splitfreqSS, chainParams.numStepsSS);

            free(tempX);
            }
        }

#   if defined (MPI_ENABLED)
    /* we need to collect the information on the number of accepted moves if
       this is a parallel version */
    if (ReassembleMoveInfo() == ERROR)
        nErrors++;
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   ReassembleMoveInfo failed\n", spacer);
        return ERROR;
        }
#   endif

    /* print acceptance rates for the moves */
    for (j=0; j<chainParams.numChains*chainParams.numRuns; j++)
        {
        if (chainParams.numChains == 1)
            {
            if (chainParams.numRuns == 1)
                MrBayesPrint ("%s   Acceptance rates for the moves:\n", spacer);
            else
                MrBayesPrint ("%s   Acceptance rates for the moves in run %d:\n", spacer, j/chainParams.numChains+1);
            }
        else if (j % chainParams.numChains == 0)
            {
            if (chainParams.numRuns == 1)
                MrBayesPrint ("\n%s   Acceptance rates for the moves in the \"cold\" chain:\n", spacer);
            else
                MrBayesPrint ("\n%s   Acceptance rates for the moves in the \"cold\" chain of run %d:\n", spacer, j/chainParams.numChains+1);
            }
        else if (chainParams.allChains == YES)
            {
            if (chainParams.numRuns == 1)
                MrBayesPrint ("\n%s   Acceptance rates for the moves in chain %d (heated):\n\n", spacer, j+1);
            else
                MrBayesPrint ("\n%s   Acceptance rates for the moves in chain %d of run %d (heated):\n\n", spacer, j%chainParams.numChains+1, j/chainParams.numChains+1);
            }

        if (j % chainParams.numChains == 0 || chainParams.allChains == YES)
            {
            MrBayesPrint ("%s      With prob.   (last %d)   chain accepted proposals by move\n", spacer, chainParams.tuneFreq);

            for (i=0; i<numUsedMoves; i++)
                {
                mv = usedMoves[i];
                if (mv->nBatches[j] < 1)
                    MrBayesPrint ("%s          NA           NA       %s\n",
                    spacer, mv->name);
                else
                    MrBayesPrint ("%s       %6.1f %%     (%3.0f %%)     %s\n", spacer,
                    100.0*mv->nTotAccepted[j]/(MrBFlt)(mv->nTotTried[j]),
                    100.0*mv->lastAcceptanceRate[j],
                    mv->name);
                }
            }
        }

#   if defined MPI_ENABLED
    /* Redistribute move info in case it is needed in a follow-up run */
    RedistributeMoveInfo();
#   endif
    
    /* output information on the success of the chain state swap proposals */
    if (PrintSwapInfo () == ERROR)
        nErrors++;
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   PrintSwapInfo failed\n", spacer);
        return ERROR;
        }
#   else
    if (nErrors > 1)
        {
        MrBayesPrint ("%s   PrintSwapInfo failed\n", spacer);
        return ERROR;
        }
#   endif

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return (NO_ERROR);
#   endif

    if (chainParams.isSS == NO && (chainParams.numRuns > 1 && chainParams.mcmcDiagn == YES))
        {
        f = 0.0;
        for (i=0; i<numTrees; i++)
            {
            if (chainParams.stat[i].avgStdDev > f)
                f = chainParams.stat[i].avgStdDev;
            }
        if (f > 0.10)
            {
            MrBayesPrint ("\n");
            MrBayesPrint ("%s   ************************* WARNING!! ************************************  \n", spacer);
            MrBayesPrint ("%s   MrBayes suspects that your runs have not converged because the tree       \n", spacer);
            MrBayesPrint ("%s   samples are very different (average standard deviation of split frequen-  \n", spacer);
            MrBayesPrint ("%s   cies larger than 0.10 (%1.2lf)). MrBayes suggests that you run the ana-   \n", spacer, f);
            MrBayesPrint ("%s   lysis longer or try to improve the MCMC sampling efficiency by fine-      \n", spacer);
            MrBayesPrint ("%s   tuning MCMC proposal or heating parameters.                               \n", spacer);
            }
        }
    
    return (NO_ERROR);
}


int SafeSprintf (char **target, int *targetLen, char *fmt, ...)
{
    va_list    argp;
    int        retval;

    while (1) {
        /* try to print in the available space */
        va_start(argp, fmt);
#   ifdef VISUAL
        retval = _vsnprintf (*target, *targetLen, fmt, argp);
#   else
        /* With SGI IRIX this function returns 0. Allen Smith suggested using 
           _xpg5_vnsprintf, but this function needs a recent version of IRIX. 
           For now, I just allocate a large buffer for SGI, so there is no need for
           memory reallocation on this kind of machines */
        /* It appears to me that the reason for the bug is that the old version of
           this function relies on nonstandard return values for vsnprintf. With
           the current version of the code, SGI should do fine without special
           treatment.  FR  */
        retval = vsnprintf (*target, *targetLen, fmt, argp);
#   endif
        va_end(argp);

        /* if it worked, retval will be in interval [0,*targetLen-1], else -1 or true length */
        if (retval > -1 && retval < *targetLen)
            return NO_ERROR;

        /* readjust length */
        if (retval > -1)     /* some C compilers will return true length in retval */
            *targetLen = retval + 1;     /* exactly what is needed */
        else                 /* some C compilers will return -1 on buffer overwrite */
            *targetLen += TARGETLENDELTA;

        /* reallocate target */
        *target = SafeRealloc ((void *)*target, *targetLen);
        if (*target == NULL)
            return ERROR;
        }
}


void SetChainIds (void)
{
    /* Fill <chainId[]> with the global chain number.
       Ex. For proc_0, chain[0] = 0;
             chain[1] = 1;
             chain[2] = 2; (numchains = 3)
       For proc_1, chain[0] = 3;
             chain[1] = 4; (numchains = 2)
       etc... 
    */

#   if defined (MPI_ENABLED)

    int     i, proc, numChainsForProc, numGlobalChains;
    int     id;
    int     remainder;

    /* calculate global number of chains */
    numGlobalChains = chainParams.numChains * chainParams.numRuns;
    
    /* there are <remainder> chains left over after
       load balancing the chains */
    remainder = numGlobalChains % num_procs;

    /* get the number of chains handled by this proc */
    numChainsForProc = (int) (numGlobalChains / num_procs);

    /* we must distribute the remaining chains (causing
       the chain load between procs to become unbalanced) */
    if (proc_id < remainder) 
        numChainsForProc++;

    /* NOTE: at this point, procs can have different number of numChainsForProc
       (one more for chains < remainder, one less for procs larger than or equal
       to the remainder) */

    id = 0;
    for (proc=0; proc<num_procs; proc++)
        {
        /* assign or increment chain id */
        if (proc == proc_id)
            {
            for (i=0; i<numChainsForProc; i++)
                chainId[i] = id++;
            }
        else
            {
            /* procs below the remainder have 1 more chain
               than procs above */
            if (proc < remainder) 
                {
                for (i=0; i<(numGlobalChains / num_procs) + 1; i++)
                    id++;
                }
            /* procs above the remainder have one less chain
               than procs below */
            else
                {
                for (i=0; i<(numGlobalChains / num_procs); i++)
                    id++;
                }
            }
        }
#   else

    int chn;
    
    for (chn=0; chn<numLocalChains; chn++)
        chainId[chn] = chn;

#   endif
        
}


/* It sets chainParams.tFilePos[] to point immidiatly after sampled tree in position "samplePos" for all .t files. */
int SetFilePositions (int samplePos)
{
    int i, j, k, longestLine;
    BitsLong    lastBlock;
    char    *lineBuf;
    FILE    *fp;
    char    *tempStr;
    int     tempStrSize = TEMPSTRSIZE;

#   if defined (MPI_ENABLED)
    if (proc_id != 0)
        return (NO_ERROR);
#   endif

    tempStr = (char *) SafeMalloc((size_t)tempStrSize * sizeof(char));
    if (!tempStr)
        {
        MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
        return (ERROR);
        }

    for (i=0; i<numTopologies; i++)
        {
        for (j=0; j<chainParams.numRuns; j++)
            {
            if (numPrintTreeParams == 1)
                {
                if (chainParams.numRuns == 1)
                    SafeSprintf (&tempStr, &tempStrSize, "%s.t", chainParams.chainFileName);
                else
                    SafeSprintf (&tempStr, &tempStrSize, "%s.run%d.t", chainParams.chainFileName, j+1);
                }
            else
                {
                if (chainParams.numRuns == 1)
                    SafeSprintf (&tempStr, &tempStrSize, "%s.tree%d.t", chainParams.chainFileName, i+1);
                else
                    SafeSprintf (&tempStr, &tempStrSize, "%s.tree%d.run%d.t", chainParams.chainFileName, i+1, j+1);
                }

            if ((fp = OpenBinaryFileR (tempStr)) == NULL) 
                {
                MrBayesPrint ("%s   Problem openning file %s.\n", spacer, tempStr);
                free (tempStr);
                return (ERROR);
                }

            longestLine = LongestLine (fp);
            SafeFclose (&fp);

            if ((fp = OpenTextFileR (tempStr)) == NULL) 
                {
                free (tempStr);
                return (ERROR);
                }

            lineBuf = (char *) SafeCalloc (longestLine + 10, sizeof (char));
            if (!lineBuf)
                {
                SafeFclose (&fp);
                free (tempStr);
                return (ERROR);
                }

            lastBlock = LastBlock (fp, lineBuf, longestLine);
            fseek (fp, lastBlock, SEEK_SET);
            fseek (fp, FirstTree(fp, lineBuf, longestLine), SEEK_SET);

            for (k=0; k<samplePos; k++)
                {
                if (fgets (lineBuf, longestLine + 5, fp) == NULL) 
                    {
                    MrBayesPrint ("%s   Not enough records in file %s.\n", spacer, tempStr);
                    SafeFclose (&fp);
                    free (tempStr);
                    free (lineBuf);
                    return ERROR;
                    }
                }
            fgetpos (fp, &chainParams.tFilePos[j*numTopologies+i]);

            SafeFclose (&fp);
            free (lineBuf);
            } /* next run */
        } /* next tree */
    free (tempStr);
    return (NO_ERROR);
}


/* SetFileNames: Set file names */
void SetFileNames (void)
{
    strcpy (sumtParams.sumtFileName, chainParams.chainFileName);
    strcpy (sumtParams.sumtOutfile, chainParams.chainFileName);
    strcpy (sumpParams.sumpFileName, chainParams.chainFileName);
    strcpy (sumpParams.sumpOutfile, chainParams.chainFileName);
    if (chainParams.numRuns == 1)
        {
        sprintf (comptreeParams.comptFileName1, "%s.t", chainParams.chainFileName);
        sprintf (comptreeParams.comptFileName2, "%s.t", chainParams.chainFileName);
        sprintf (plotParams.plotFileName, "%s.p", chainParams.chainFileName);
        MrBayesPrint ("%s   Setting chain output file names to \"%s.<p/t>\"\n", spacer, chainParams.chainFileName);
        }
    else /* if (chainParams.numRuns > 1) */
        {
        sprintf (comptreeParams.comptFileName1, "%s.run1.t", chainParams.chainFileName);
        sprintf (comptreeParams.comptFileName2, "%s.run2.t", chainParams.chainFileName);
        sprintf (plotParams.plotFileName, "%s.run1.p", chainParams.chainFileName);
        MrBayesPrint ("%s   Setting chain output file names to \"%s.run<i>.<p/t>\"\n", spacer, chainParams.chainFileName);
        }
}


/*----------------------------------------------------------------------------
|
|   SetLikeFunctions: This function will set up the pointers from each
|       data division to the relevant likelihood functions. It will also
|       deal with the settings for SIMD code.
|
-----------------------------------------------------------------------------*/
int SetLikeFunctions (void)
{
    int         i;

    ModelInfo   *m;

    /* couple divisions with likelihood functions */
    for (i=0; i<numCurrentDivisions; i++)
        {
        m = &modelSettings[i];
        m->useVec = VEC_NONE;
        
        if (m->dataType == DNA || m->dataType == RNA)
            {
            if (m->parsModelId == YES)
                {
                m->Likelihood = &Likelihood_Pars;
                }
            else
                {
                if (m->nucModelId == NUCMODEL_4BY4)
                    {
                    if (m->numModelStates > 4)
                        {
                        /* covariotide model */
                        /* TODO: allow autocorrelated rates */
                        if (m->gibbsGamma == YES)
                            {
                            m->CondLikeDown = &CondLikeDown_Gen_GibbsGamma;
                            m->CondLikeRoot = &CondLikeRoot_Gen_GibbsGamma;
                            m->CondLikeScaler = &CondLikeScaler_Gen_GibbsGamma;
                            m->Likelihood = &Likelihood_Gen_GibbsGamma;
                            }
                        else
                            {
                            m->CondLikeDown = &CondLikeDown_Gen;
                            m->CondLikeRoot = &CondLikeRoot_Gen;
                            m->CondLikeScaler = &CondLikeScaler_Gen;
                            m->Likelihood = &Likelihood_Gen;
                            }
                        if (m->correlation != NULL)
                            m->Likelihood = &Likelihood_Adgamma;
                        else
                            m->Likelihood = &Likelihood_Gen;
                        if (m->nCijkParts == 1)
                            m->TiProbs = &TiProbs_Gen;
                        else if (m->nCijkParts > 1)
                            m->TiProbs = &TiProbs_GenCov;
                        m->PrintAncStates = &PrintAncStates_NUC4;
                        m->PrintSiteRates = &PrintSiteRates_Gen;
                        }
                    else
                        {
#   if defined (SSE_ENABLED)
                        if (m->printAncStates == YES || m->printSiteRates == YES)
                            {
                            MrBayesPrint ("%s   Non-SSE version of conditional likelihood calculator will be used for division %d\n", spacer, i+1);
                            MrBayesPrint ("%s   due to request of reporting 'ancestral states' or 'site rates'.\n", spacer);
                            }

                        m->CondLikeUp = &CondLikeUp_NUC4;
                        m->PrintAncStates = &PrintAncStates_NUC4;
                        m->PrintSiteRates = &PrintSiteRates_Gen;

                        if (m->gibbsGamma == YES)
                            {
                            m->CondLikeDown = &CondLikeDown_NUC4_GibbsGamma;
                            m->CondLikeRoot = &CondLikeRoot_NUC4_GibbsGamma;
                            m->CondLikeScaler = &CondLikeScaler_NUC4_GibbsGamma;
                            }
                        else if (m->correlation != NULL || m->printAncStates == YES || m->printSiteRates == YES)
                            {
                            m->CondLikeDown = &CondLikeDown_NUC4;
                            m->CondLikeRoot = &CondLikeRoot_NUC4;
                            m->CondLikeScaler = &CondLikeScaler_NUC4;
                            }
                        else
                            {
                            m->useVec = VEC_SSE;
                            m->numFloatsPerVec = 4;
                            m->CondLikeDown = &CondLikeDown_NUC4_SSE;
                            m->CondLikeRoot = &CondLikeRoot_NUC4_SSE;
                            m->CondLikeScaler = &CondLikeScaler_NUC4_SSE;
#if defined (AVX_ENABLED)   // override SSE settings
                            m->useVec = VEC_AVX;
                            m->numFloatsPerVec = 8;
                            m->CondLikeDown = &CondLikeDown_NUC4_AVX;
                            m->CondLikeRoot = &CondLikeRoot_NUC4_AVX;
                            m->CondLikeScaler = &CondLikeScaler_NUC4_AVX;
#if defined (FMA_ENABLED)   // override AVX settings (CondLikeScaler cannot be improved over AVX)
                            m->useVec = VEC_FMA;
                            m->CondLikeDown = &CondLikeDown_NUC4_FMA;
                            m->CondLikeRoot = &CondLikeRoot_NUC4_FMA;
#endif
#endif
                            /* Should be sse versions if we want to handle m->printAncStates == YES || inferSiteRates == YES.
                            For now just set to NULL for early error detection if functions anyway got called by mistake */
                            m->CondLikeUp = NULL;
                            m->PrintAncStates = NULL;
                            m->PrintSiteRates = NULL;
                            }

                        if (m->correlation != NULL)
                            m->Likelihood = &Likelihood_Adgamma;
                        else if (m->gibbsGamma == YES)
                            m->Likelihood = &Likelihood_NUC4_GibbsGamma;
                        else if (m->printAncStates == YES || inferSiteRates == YES)
                            m->Likelihood = &Likelihood_NUC4;
                        else
                            {
                            m->Likelihood = &Likelihood_NUC4_SSE;
#if defined (AVX_ENABLED)
                            m->Likelihood = &Likelihood_NUC4_AVX;
#if defined (FMA_ENABLED)
                            m->Likelihood = &Likelihood_NUC4_FMA;
#endif
#endif
                            }
#   else
                        if (m->gibbsGamma == YES)
                            {
                            m->CondLikeDown = &CondLikeDown_NUC4_GibbsGamma;
                            m->CondLikeRoot = &CondLikeRoot_NUC4_GibbsGamma;
                            m->CondLikeScaler = &CondLikeScaler_NUC4_GibbsGamma;
                            }
                        else
                            {
                            m->CondLikeDown = &CondLikeDown_NUC4;
                            m->CondLikeRoot = &CondLikeRoot_NUC4;
                            m->CondLikeScaler = &CondLikeScaler_NUC4;
                            }

                        if (m->correlation != NULL)
                            m->Likelihood = &Likelihood_Adgamma;
                        else if (m->gibbsGamma == YES)
                            m->Likelihood = &Likelihood_NUC4_GibbsGamma;
                        else
                            m->Likelihood = &Likelihood_NUC4;

                        m->CondLikeUp = &CondLikeUp_NUC4;
                        m->PrintAncStates = &PrintAncStates_NUC4;
                        m->PrintSiteRates = &PrintSiteRates_Gen;
#   endif
                        if (m->nst == 1)
                            m->TiProbs = &TiProbs_Fels;
                        else if (m->nst == 2)
                            m->TiProbs = &TiProbs_Hky;
                        else
                            m->TiProbs = &TiProbs_Gen;
                        m->StateCode = &StateCode_NUC4;
                        }
                    }
                else if (m->nucModelId == NUCMODEL_DOUBLET)
                    {
                    if (m->gibbsGamma == YES)
                        {
                        m->CondLikeDown = &CondLikeDown_Gen_GibbsGamma;
                        m->CondLikeRoot = &CondLikeRoot_Gen_GibbsGamma;
                        m->CondLikeScaler = &CondLikeScaler_Gen_GibbsGamma;
                        m->Likelihood = &Likelihood_Gen_GibbsGamma;
                        }
                    else
                        {
                        m->CondLikeDown = &CondLikeDown_Gen;
                        m->CondLikeRoot = &CondLikeRoot_Gen;
                        m->CondLikeScaler = &CondLikeScaler_Gen;
                        m->Likelihood = &Likelihood_Gen;
                        }
                    if (m->nst == 1)
                        m->TiProbs = &TiProbs_Gen;
                    else if (m->nst == 2)
                        m->TiProbs = &TiProbs_Gen;
                    else
                        m->TiProbs = &TiProbs_Gen;
                    m->CondLikeUp = &CondLikeUp_Gen;
                    m->PrintAncStates = &PrintAncStates_Gen;
                    m->PrintSiteRates = &PrintSiteRates_Gen;
                    }
                else if (m->nucModelId == NUCMODEL_CODON)
                    {
                    /* codon models */
                    if (m->numOmegaCats == 1)
                        {
                        if (m->gibbsGamma == YES)
                            {
                            m->CondLikeDown = &CondLikeDown_Gen_GibbsGamma;
                            m->CondLikeRoot = &CondLikeRoot_Gen_GibbsGamma;
                            m->CondLikeScaler = &CondLikeScaler_Gen_GibbsGamma;
                            m->Likelihood = &Likelihood_Gen_GibbsGamma;
                            }
                        else
                            {
                            m->CondLikeDown = &CondLikeDown_Gen;
                            m->CondLikeRoot = &CondLikeRoot_Gen;
                            m->CondLikeScaler = &CondLikeScaler_Gen;
                            m->Likelihood = &Likelihood_Gen;
#   if defined (SSE_ENABLED)
                            if (m->printAncStates == YES || m->printSiteRates == YES)
                                {
                                MrBayesPrint ("%s   Non-SSE version of conditional likelihood calculator will be used for division %d\n", spacer, i+1);
                                MrBayesPrint ("%s   due to request of reporting 'ancestral states' or 'site rates'.\n", spacer);
                                }
                            else
                                {
                                m->useVec = VEC_SSE;
                                m->numFloatsPerVec = 4;
                                m->CondLikeDown = &CondLikeDown_Gen_SSE;
                                m->CondLikeRoot = &CondLikeRoot_Gen_SSE;
                                m->CondLikeScaler = &CondLikeScaler_Gen_SSE;
                                m->Likelihood = &Likelihood_Gen_SSE;
                                }
#   endif
                            }
                        }
                    else
                        {
                        m->CondLikeDown   = &CondLikeDown_NY98;
                        m->CondLikeRoot   = &CondLikeRoot_NY98;
                        m->CondLikeScaler = &CondLikeScaler_NY98;
                        m->Likelihood     = &Likelihood_NY98;
                        m->PosSelProbs    = &PosSelProbs;
                        m->SiteOmegas     = &SiteOmegas;
#   if defined (SSE_ENABLED)
                        if (m->printAncStates == YES || m->printSiteRates == YES)
                            {
                            MrBayesPrint ("%s   Non-SSE version of conditional likelihood calculator will be used for division %d\n", spacer, i+1);
                            MrBayesPrint ("%s   due to request of reporting 'ancestral states' or 'site rates'.\n", spacer);
                            }
                        else
                            {
                            m->useVec = VEC_SSE;
                            m->numFloatsPerVec = 4;
                            m->CondLikeDown   = &CondLikeDown_NY98_SSE;
                            m->CondLikeRoot   = &CondLikeRoot_NY98_SSE;
                            m->CondLikeScaler = &CondLikeScaler_NY98_SSE;
                            m->Likelihood     = &Likelihood_NY98_SSE;
                            m->PosSelProbs    = &PosSelProbs_SSE;
                            m->SiteOmegas     = &SiteOmegas_SSE;
                            }
#   endif
                        }
                    m->TiProbs = &TiProbs_Gen;
                    if (m->nCijkParts > 1)
                        m->TiProbs = &TiProbs_GenCov;
                    m->CondLikeUp = &CondLikeUp_Gen;
                    m->PrintAncStates = &PrintAncStates_Gen;
                    m->PrintSiteRates = &PrintSiteRates_Gen;
                    }
                else /* if (m->nucModelId == NUCMODEL_AA) */
                    {
                    if (m->gibbsGamma == YES)
                        {
                        m->CondLikeDown = &CondLikeDown_Gen_GibbsGamma;
                        m->CondLikeRoot = &CondLikeRoot_Gen_GibbsGamma;
                        m->CondLikeScaler = &CondLikeScaler_Gen_GibbsGamma;
                        m->Likelihood = &Likelihood_Gen_GibbsGamma;
                        }
                    else
                        {
                        m->CondLikeDown = &CondLikeDown_Gen;
                        m->CondLikeRoot = &CondLikeRoot_Gen;
                        m->CondLikeScaler = &CondLikeScaler_Gen;
                        m->Likelihood = &Likelihood_Gen;
#   if defined (SSE_ENABLED)
                        if (m->printAncStates == YES || m->printSiteRates == YES)
                            {
                            MrBayesPrint ("%s   Non-SSE version of conditional likelihood calculator will be used for division %d\n", spacer, i+1);
                            MrBayesPrint ("%s   due to request of reporting 'ancestral states' or 'site rates'.\n", spacer);
                            }
                        else
                            {
                            m->useVec = VEC_SSE;
                            m->numFloatsPerVec = 4;
                            m->CondLikeDown = &CondLikeDown_Gen_SSE;
                            m->CondLikeRoot = &CondLikeRoot_Gen_SSE;
                            m->CondLikeScaler = &CondLikeScaler_Gen_SSE;
                            m->Likelihood = &Likelihood_Gen_SSE;
                            }
#   endif
                        }
                    if (m->nCijkParts > 1)
                        m->TiProbs = &TiProbs_GenCov;
                    else
                        m->TiProbs = &TiProbs_Gen;
                    m->CondLikeUp = &CondLikeUp_Gen;
                    m->StateCode = &StateCode_AA;
                    m->PrintAncStates = &PrintAncStates_Gen;
                    m->PrintSiteRates = &PrintSiteRates_Gen;
                    }
                }
            }
        else if (m->dataType == PROTEIN)
            {
            if (m->parsModelId == YES)
                {
                m->Likelihood = &Likelihood_Pars;
                }
            else
                {
                /* TODO: allow autocorrelated rates for covarion model */
                if (m->gibbsGamma == YES)
                    {
                    m->CondLikeDown = &CondLikeDown_Gen_GibbsGamma;
                    m->CondLikeRoot = &CondLikeRoot_Gen_GibbsGamma;
                    m->CondLikeScaler = &CondLikeScaler_Gen_GibbsGamma;
                    m->Likelihood = &Likelihood_Gen_GibbsGamma;
                    }
                else
                    {
                    m->CondLikeDown = &CondLikeDown_Gen;
                    m->CondLikeRoot = &CondLikeRoot_Gen;
                    m->CondLikeScaler = &CondLikeScaler_Gen;
                    m->Likelihood = &Likelihood_Gen;
#   if defined (SSE_ENABLED)
                    if (m->printAncStates == YES || m->printSiteRates == YES)
                        {
                        MrBayesPrint ("%s   Non-SSE version of conditional likelihood calculator will be used for division %d\n", spacer, i+1);
                        MrBayesPrint ("%s   due to request of reporting 'ancestral states' or 'site rates'.\n", spacer);
                        }
                    else
                        {
                        m->useVec = VEC_SSE;
                        m->numFloatsPerVec = 4;
                        m->CondLikeDown = &CondLikeDown_Gen_SSE;
                        m->CondLikeRoot = &CondLikeRoot_Gen_SSE;
                        m->CondLikeScaler = &CondLikeScaler_Gen_SSE;
                        m->Likelihood = &Likelihood_Gen_SSE;
                        }
#   endif
                    }
                if (m->correlation != NULL)
                    {
                    if (m->gibbsGamma == YES)
                        {
                        MrBayesPrint ("%s   Adgamma model cannot be used with Gibbs sampling of rate categories\n", spacer);
                        return (ERROR);
                        }
                    else
                        m->Likelihood = &Likelihood_Adgamma;
                    }
                if (m->numModelStates > 20 && m->nCijkParts > 1)
                    m->TiProbs = &TiProbs_GenCov;
                else
                    m->TiProbs = &TiProbs_Gen;
                m->CondLikeUp = &CondLikeUp_Gen;
                m->StateCode = &StateCode_AA;
                m->PrintAncStates = &PrintAncStates_Gen;
                m->PrintSiteRates = &PrintSiteRates_Gen;
                }
            }
        else if (m->dataType == RESTRICTION)
            {
            if (m->parsModelId == YES)
                {
                m->Likelihood = &Likelihood_Pars;
                }
            else
                {
                m->CondLikeDown   = &CondLikeDown_Bin;
                m->CondLikeRoot   = &CondLikeRoot_Bin;
                m->CondLikeScaler = &CondLikeScaler_Gen;
                m->Likelihood     = &Likelihood_Res;
#   if defined (SSE_ENABLED)
                if (m->printAncStates == YES || m->printSiteRates == YES)
                    {
                    MrBayesPrint ("%s   Non-SSE version of conditional likelihood calculator will be used for division %d\n", spacer, i+1);
                    MrBayesPrint ("%s   due to request of reporting 'ancestral states' or 'site rates'.\n", spacer);
                    }
                else
                    {
                    m->useVec = VEC_SSE;
                    m->numFloatsPerVec = 4;
                    m->CondLikeDown   = &CondLikeDown_Bin_SSE;
                    m->CondLikeRoot   = &CondLikeRoot_Bin_SSE;
                    m->CondLikeScaler = &CondLikeScaler_Gen_SSE;
                    m->Likelihood     = &Likelihood_Res_SSE;
                    }
#   endif
                m->TiProbs = &TiProbs_Res;
                m->CondLikeUp = &CondLikeUp_Bin;
                m->StateCode = &StateCode_Std;
                m->PrintAncStates = &PrintAncStates_Bin;
                m->PrintSiteRates = &PrintSiteRates_Gen;
                }
            }
        else if (m->dataType == STANDARD)
            {
            if (m->parsModelId == YES)
                {
                if (m->numModelStates == 2)
                    {
                    m->Likelihood = &Likelihood_Pars; /* this is much faster if number of states do not vary */
                    m->numStates = 2;   /* this is needed for the parsimony calculator */
                    }
                else
                    m->Likelihood = &Likelihood_ParsStd;
                }
            else
                {
                m->CondLikeDown   = &CondLikeDown_Std;
                m->CondLikeRoot   = &CondLikeRoot_Std;
                m->CondLikeScaler = &CondLikeScaler_Std;
                m->Likelihood     = &Likelihood_Std;
                m->TiProbs        = &TiProbs_Std;
                m->CondLikeUp     = &CondLikeUp_Std;
                m->StateCode      = &StateCode_Std;
                m->PrintAncStates = &PrintAncStates_Std;
                m->PrintSiteRates = &PrintSiteRates_Std;
                }
            }       
        else if (m->dataType == CONTINUOUS)
            {
            
            }       
        else
            {
            MrBayesPrint ("%s   ERROR: Data should be one of these types!\n", spacer);
            return ERROR;
            }
        
        }

    return NO_ERROR;
}


/* Determine number of chains and data splits to be handled by MPI processors or threads */
int SetLocalChainsAndDataSplits(void)
{
#   if defined (MPI_ENABLED)
    
    /* tell user how many chains each processor has been assigned */
    if (num_procs <= numGlobalChains)
        {
        if (numGlobalChains % num_procs != 0)
            {
            MrBayesPrint ("%s   The total number of chains (%d) must be evenly divisible by\n", spacer, numGlobalChains);
            MrBayesPrint ("%s   the number of MPI processors (%d).\n", spacer, num_procs);
            /* 
            MrBayesPrint ("%s   The total number of chains (%d) must be evenly divisible by\n", spacer, numGlobalChains);
            MrBayesPrint ("%s   the number of MPI processors (%d), or the number of MPI\n", spacer, num_procs);
            MrBayesPrint ("%s   processors should be a multiple of the number of chains.\n", spacer, num_procs);
            */
            MrBayesPrint ("%s   Please change your MPI settings.\n", spacer, num_procs);
            return ERROR;
            }
        numLocalChains = numGlobalChains / num_procs;
        MrBayesPrint ("%s   Number of chains per MPI processor = %d\n", spacer, numLocalChains);
        }
    else
        {
        MrBayesPrint ("%s   There must be at least as many chains as MPI processors\n", spacer);
        return (ERROR);
        }
        /* 
        {
        if (num_procs % numGlobalChains != 0)
            {
            MrBayesPrint ("%s   The number of MPI processors (%d) must be a multiple of the\n", spacer, num_procs);
            MrBayesPrint ("%s   total number of chains (%d), or the total number of chains\n", spacer, numGlobalChains);
            MrBayesPrint ("%s   should be evenly divisible by the number of MPI processsors.\n", spacer);
            MrBayesPrint ("%s   Please change your MPI settings.\n", spacer);
            return ERROR;
            }
        numLocalChains = 1;
        // numMPIDataSplits = num_procs / numGlobalChains;
        // MrBayesPrint ("%s   Number of MPI data splits per chain = %d\n", spacer, numMPIDataSplits);
        }
        */
#   else

    numLocalChains = numGlobalChains;

#   endif

#   if defined (THREADS_ENABLED)

    if (numLocalChains > 1)
        {
        /* Use pthreads to divide chains and possibly do data splits */
        }
    else
        {
        /* Use pthreads for data splits */
        }
#   endif

    return (NO_ERROR);
}


/*----------------------------------------------------------------------------
|
|   ShowMoveSummary: Show summary of moves that will be used in MCMC sampling
|
-----------------------------------------------------------------------------*/
int ShowMoveSummary (void)
{
    int         i, run, chain, areRunsSame, areChainsSame, chainIndex;
    MCMCMove    *mv;
    MrBFlt      prob;
    
    /* first find out if the probabilities are different in different runs */           
    areRunsSame = YES;
    for (run=1; run<chainParams.numRuns; run++)
        {
        for (chain=0; chain<chainParams.numChains; chain++)
            {
            chainIndex = run*chainParams.numChains + chain;
            for (i=0; i<numUsedMoves; i++)
                {
                mv = usedMoves[i];
                if (AreDoublesEqual (mv->relProposalProb[chainIndex], mv->relProposalProb[chain], 0.000001) == NO)
                    {
                    areRunsSame = NO;
                    break;
                    }
                }
            if (areRunsSame == NO)
                break;
            }
        if (areRunsSame == NO)
            break;
        }

    /* now print values */
    for (run=0; run<chainParams.numRuns; run++)
        {
        if (areRunsSame == YES && run >= 1)
            break;

        /* find out if chains are different within this run */
        areChainsSame = YES;
        for (chain=1; chain<chainParams.numChains; chain++)
            {
            chainIndex = run*chainParams.numChains + chain;
            for (i=0; i<numUsedMoves; i++)
                {
                mv = usedMoves[i];
                if (AreDoublesEqual (mv->relProposalProb[chainIndex], mv->relProposalProb[chainIndex-chain], 0.000001) == NO)
                    {
                    areChainsSame = NO;
                    break;
                    }
                }
            if (areChainsSame == NO)
                break;
            }

        for (chain=0; chain<chainParams.numChains; chain++)
            {
            if (areChainsSame == YES && chain >= 1)
                break;
            
            /* now we can print the values */
            MrBayesPrint("\n");
            if (areRunsSame == YES && areChainsSame == YES)
                MrBayesPrint ("%s   The MCMC sampler will use the following moves:\n", spacer);
            else if (areRunsSame == NO && areChainsSame == YES)
                MrBayesPrint ("%s   The MCMC sampler will use the following moves for run %d:\n", spacer, run+1);
            else if (areRunsSame == YES && areChainsSame == NO)
                MrBayesPrint ("%s   The MCMC sampler will use the following moves for chain %d:\n", spacer, chain+1);
            else if (areRunsSame == NO && areChainsSame == NO)
                MrBayesPrint ("%s   The MCMC sampler will use the following moves for run %d, chain %d:\n", spacer, run+1, chain+1);

            chainIndex = run*chainParams.numChains + chain;
            MrBayesPrint ("%s      With prob.  Chain will use move\n", spacer);
            for (i=0; i<numUsedMoves; i++)
                {
                mv = usedMoves[i];
                prob = mv->cumProposalProb[chainIndex];
                if (i > 0)
                    prob -= usedMoves[i-1]->cumProposalProb[chainIndex];
                if (AreDoublesEqual(prob,0.0,0.000001) == YES)
                    continue;
                MrBayesPrint ("%s       %6.2f %%   %s\n", spacer, 100*prob, mv->name);
                }
            MrBayesPrint ("\n");
            }   /* next chain */
        }   /* next run */

    return (NO_ERROR);
}


/*----------------------------------------------------------------------
|
|   SetUpPartitionCounters: Set up partitions and the root of the
|      partition frequency tree                         
|
|----------------------------------------------------------------------*/
int SetUpPartitionCounters (void)
{
    int     i;
    
#   if defined (MPI_ENABLED)
    /* we only keep partition counters on proc 0 in the MPI version */
    if (proc_id != 0)
        return (NO_ERROR);
#   endif
    nLongsNeeded = 1 + (numLocalTaxa - 1) / nBitsInALong;
    
    if (memAllocs[ALLOC_PFCOUNTERS] == YES)
        {
        MrBayesPrint ("%s   ERROR: pfcounters not free in SetUpPartitionCounters\n", spacer);
        return ERROR;
        }
    partition = (BitsLong **) SafeCalloc (2*numLocalTaxa, sizeof (BitsLong *));
    if (partition == NULL)
        {
        MrBayesPrint ("%s   Failed to allocate partition in SetUpPartitionCounters\n", spacer);
        return ERROR;
        }
    partition[0] = (BitsLong *) SafeCalloc (2*numLocalTaxa * nLongsNeeded, sizeof(BitsLong));
    if (partition[0] == NULL)
        {
        free (partition);
        MrBayesPrint ("%s   Failed to allocate partition[0] in SetUpPartitionCounters\n", spacer);
        return ERROR;
        }
    partFreqTreeRoot = (PFNODE **) SafeCalloc (numTopologies, sizeof (PFNODE *));
    if (partFreqTreeRoot == NULL)
        {
        free (partition);
        free (partition[0]);
        MrBayesPrint ("%s   Failed to allocate partFreqTreeRoot in SetUpPartitionCounters\n", spacer);
        return ERROR;
        }
    memAllocs[ALLOC_PFCOUNTERS] = YES;

    for (i=1; i<2*numLocalTaxa; i++)
        {
        partition[i] = partition[0] + i*nLongsNeeded;
        }
    
    for (i=0; i<numLocalTaxa; i++)
        SetBit (i, partition[i]);

    for (i=0; i<numTopologies; i++)
        partFreqTreeRoot[i] = NULL;

    return NO_ERROR;
}


/*----------------------------------------------------------------------
|
|   SetupTermState: create matrix holding unambiguous states for
|       terminals (used for local compression on terminal branches)
|
-----------------------------------------------------------------------*/
int SetUpTermState (void)
{
    int         i, k, n, c, d, x=0, *termStatePtr;
    BitsLong    *p;
    ModelInfo   *m;
    ModelParams *mp;
    int         numComprChars = 0;

    /* allocate space for termState and isPartAmbig */
    if (memAllocs[ALLOC_TERMSTATE] == YES || memAllocs[ALLOC_ISPARTAMBIG] == YES)
        {
        MrBayesPrint ("%s   termState or isPartAmbig is not free in SetupTermState\n", spacer);
        return ERROR;
        }

#   if defined SSE_ENABLED
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        if ( m->useVec != VEC_NONE )
            {
            m->numVecChars = ((m->numChars - 1) / m->numFloatsPerVec) + 1;
            numComprChars += m->numVecChars * m->numFloatsPerVec;
            }
        else
            numComprChars += m->numChars;
        }
#   else
    numComprChars = numCompressedChars;
#   endif

    termState = (int *) SafeCalloc (numLocalTaxa * numComprChars, sizeof(int));
    if (termState)
        memAllocs[ALLOC_TERMSTATE] = YES;
    else
        {
        MrBayesPrint ("%s   Problem allocating termState\n", spacer);
        return (ERROR);
        }
    isPartAmbig = (int *) SafeCalloc (numLocalTaxa*numCurrentDivisions, sizeof(int));
    if (isPartAmbig)
        memAllocs[ALLOC_ISPARTAMBIG] = YES;
    else
        {
        MrBayesPrint ("%s   Problem allocating isPartAmbig\n", spacer);
        return (ERROR);
        }

    /*initialize isPartAmbig */
    for (i=0; i<numLocalTaxa*numCurrentDivisions; i++)
        isPartAmbig[i] = NO;

    /* loop over divisions */
    termStatePtr = termState;
    for (d=0; d<numCurrentDivisions; d++)
        {
        m = &modelSettings[d];
        mp = &modelParams[d];

        /* don't do anything for continuous data */
        if (mp->dataType == CONTINUOUS)
            continue;
        
        m->termState   = (int **) SafeCalloc (numLocalTaxa, sizeof(int *));
        if (!m->termState)
            {
            MrBayesPrint("%s   Problems allocating termState pointers for division %d\n", spacer, d+1);
            return ERROR;
            }

#   if defined SSE_ENABLED
        if (m->dataType != STANDARD && m->gibbsGamma == NO)
            numComprChars = m->numVecChars * m->numFloatsPerVec;
        else
            numComprChars = m->numChars;

#   else
            numComprChars = m->numChars;
#   endif

        for (i=0; i<numLocalTaxa; i++)
            {
            m->termState[i] = termStatePtr;
            termStatePtr += numComprChars;
            }

        m->isPartAmbig = isPartAmbig + numLocalTaxa * d;

        for (i=0; i<numLocalTaxa; i++)
            {
            p = m->parsSets[i];
            for (c=0; c<m->numChars; c++)
                {
                for (k=n=0; k<m->numStates; k++)
                    {
                    if (IsBitSet(k, p))
                        {
                        x = k;
                        n++;
                        }
                    }
                /* find appropriate index */
                if (n == 1)
                    m->termState[i][c] = x * m->numModelStates;
                else if (n == m->numStates)
                    m->termState[i][c] = m->numStates * m->numModelStates;
                else
                    m->isPartAmbig[i] = YES;

                p += m->nParsIntsPerSite;
                }
            for (; c<numComprChars; c++)
                {
                /* Setting to fully ambig state all padding chars */
                m->termState[i][c] = m->numStates * m->numModelStates;
                }
            }
        }

    /* print the termState matrix */
#   if  defined (DEBUG_SETUPTERMSTATE)
    PrintTermState();
    getchar();
#   endif

    return NO_ERROR;
}


/*----------------------------------------------------------------------------
|
|   SetUsedMoves: This function will extract an array of pointers to the moves
|      that will actually be used in the MCMC sampling. It also makes sure
|      that the parsimonyBased flag is set in the relevant model partitions
|      if there are moves used that are based on parsimony scores
|
-----------------------------------------------------------------------------*/
int SetUsedMoves (void)
{
    int         i, j, moveIndex, numGlobalChains;
    MrBFlt      prob, sum, cumSum;

    /* first count moves */
    numUsedMoves = 0;
    numGlobalChains = chainParams.numChains * chainParams.numRuns;
    for (i=0; i<numApplicableMoves; i++)
        {
        prob = 0.0;
        for (j=0; j<numGlobalChains; j++)
            {
            if (moves[i]->relProposalProb[j] > prob)
                prob = moves[i]->relProposalProb[j];
            }
        if (prob > 0.000001)
            numUsedMoves++;
        }
    
    /* allocate space */
    if (memAllocs[ALLOC_USEDMOVES] == YES)
        {
        MrBayesPrint ("%s   Memory problem: usedMoves not free in SetUsedMoves\n", spacer);
        return (ERROR);
        }
    usedMoves = (MCMCMove **) SafeMalloc (numUsedMoves * sizeof (MCMCMove *));
    if (!usedMoves)
        {
        MrBayesPrint ("%s   Problem allocating usedMoves\n", spacer);
        return (ERROR);
        }
    memAllocs[ALLOC_USEDMOVES] = YES;
        
    /* set move pointers */
    moveIndex = 0;
    for (i=0; i<numApplicableMoves; i++)
        {
        prob = 0.0;
        for (j=0; j<numGlobalChains; j++)
            {
            if (moves[i]->relProposalProb[j] > prob)
                prob = moves[i]->relProposalProb[j];
            }
        if (prob > 0.000001)
            usedMoves[moveIndex++] = moves[i];
        }
    
    if (moveIndex != numUsedMoves)
        {
        MrBayesPrint ("%s   Problem finding the used moves\n", spacer);
        return (ERROR);
        }

    /* set parsimony flag if applicable */
    for (i=0; i<numCurrentDivisions; i++)
        modelSettings[i].parsimonyBasedMove = NO;
    for (i=0; i<numUsedMoves; i++)
        {
        if (usedMoves[i]->moveType->parsimonyBased == YES)
            {
            for (j=0; j<usedMoves[i]->parm->nRelParts; j++)
                modelSettings[usedMoves[i]->parm->relParts[j]].parsimonyBasedMove = YES;
            }       
        }

    /* set cumulative proposal probabilities */
    for (j=0; j<numGlobalChains; j++)
        {
        sum = 0.0;
        for (i=0; i<numUsedMoves; i++)
            {
            sum += usedMoves[i]->relProposalProb[j];
            }
        cumSum = 0.0;
        for (i=0; i<numUsedMoves; i++)
            {
            cumSum += usedMoves[i]->relProposalProb[j];
            usedMoves[i]->cumProposalProb[j] = cumSum / sum;
            }
        }

    /* reset acceptance probability values */
    for (i=0; i<numUsedMoves; i++)
        {
        for (j=0; j<numGlobalChains; j++)
            {
            usedMoves[i]->nAccepted[j] = 0;
            usedMoves[i]->nTried[j] = 0;
            usedMoves[i]->nTotAccepted[j] = 0;
            usedMoves[i]->nTotTried[j] = 0;
            }
        }

    return (NO_ERROR);
}


void ShowValuesForChain (int chn)
{
    int             i;
    char            s[100];
        
    MrBayesPrint ("%s   Chain = %d\n", spacer, chn);
    MrBayesPrint ("%s      numParams = %d\n", spacer, numParams);
    MrBayesPrint ("%s      numTrees  = %d\n", spacer, numTrees);
    MrBayesPrint ("%s      current state: %d\n", spacer, state[chn]);
    
    strcat (spacer, "   ");

    /* tRatio */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "tRatio[%d]", i);
        PrintParamValues (modelSettings[i].tRatio, chn, s);
        }

    /* revMat */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "revMat[%d]", i);
        PrintParamValues (modelSettings[i].revMat, chn, s);
        }

    /* stateFreq */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "stateFreq[%d]", i);
        PrintParamValues (modelSettings[i].stateFreq, chn, s);
        }

    /* omega */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "omega[%d]", i);
        PrintParamValues (modelSettings[i].omega, chn, s);
        }

    /* shape */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "shape[%d]", i);
        PrintParamValues (modelSettings[i].shape, chn, s);
        }

    /* pInvar */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "pInvar[%d]", i);
        PrintParamValues (modelSettings[i].pInvar, chn, s);
        }

    /* correlation */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "correlation[%d]", i);
        PrintParamValues (modelSettings[i].correlation, chn, s);
        }

    /* switchRates */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "switchRates[%d]", i);
        PrintParamValues (modelSettings[i].switchRates, chn, s);
        }

    /* rateMult */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "rateMult[%d]", i);
        PrintParamValues (modelSettings[i].rateMult, chn, s);
        }

    /* speciationRates */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "speciationRates[%d]", i);
        PrintParamValues (modelSettings[i].speciationRates, chn, s);
        }

    /* extinctionRates */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "extinctionRates[%d]", i);
        PrintParamValues (modelSettings[i].extinctionRates, chn, s);
        }

    /* fossilizationRates */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "fossilizationRates[%d]", i);
        PrintParamValues (modelSettings[i].fossilizationRates, chn, s);
        }

    /* popSize */
    for (i=0; i<numCurrentDivisions; i++)
        {
        sprintf (s, "popSize[%d]", i);
        PrintParamValues (modelSettings[i].popSize, chn, s);
        }

    /* topology */
    for (i=0; i<numCurrentDivisions; i++)
        {
        MrBayesPrint ("%s   topology[%d] state 0\n", spacer, i);
        ShowTree(GetTree (modelSettings[i].topology, chn, 0));
        MrBayesPrint ("%s   topology[%d] state 1\n", spacer, i);
        ShowTree(GetTree (modelSettings[i].topology, chn, 1));
        }
        
    /* brlens */
    for (i=0; i<numCurrentDivisions; i++)
        {
        MrBayesPrint ("%s   tree[%d] state 0\n", spacer, i);
        ShowTree(GetTree (modelSettings[i].topology, chn, 0));
        MrBayesPrint ("%s   tree[%d] state 1\n", spacer, i);
        ShowTree(GetTree (modelSettings[i].topology, chn, 1));
        }

    spacer[strlen(spacer) - 3] = '\0';

#   if  0
    for (i=0; i<sizeOfParamValues; i++)
        MrBayesPrint ("%4d -- %lf\n", i, paramValues[i]);
#   endif
}


/* SmallestNonemptyPFNode: recursive function to smallest nonempty node in a subtree */
PFNODE *SmallestNonemptyPFNode (PFNODE *p, int *i, int j)
{
    PFNODE *q;

    ++j;
    if (p == NULL)
        return NULL;
    
    q = SmallestNonemptyPFNode (p->right, i, j);
    
    if (q != NULL)
        {
        return q;
        }
    else if (IsPFNodeEmpty (p) == NO)
        {
        *i = j;
        return p;
        }
    else
        {
        return SmallestNonemptyPFNode (p->left, i, j);
        }
}


/* Talloc: Allocate space for a new node in the tree keeping track of partition frequencies */
PFNODE *Talloc (void)
{
    PFNODE  *temp;

    temp = (PFNODE *) SafeMalloc (sizeof(PFNODE));
    if (temp == NULL)
        return NULL;

    temp->partition = (BitsLong *) SafeCalloc (nLongsNeeded, sizeof (BitsLong));
    if (temp->partition == NULL)
        {
        free (temp);
        return NULL;
        }

    temp->count = (int *) SafeCalloc (chainParams.numRuns, sizeof (int));
    if (temp->count == NULL)
        {
        free (temp->partition);
        free (temp);
        return NULL;
        }

    return temp; 
}


MrBFlt Temperature (int id)
{
    /* let id be number of chain in run */
    id %= chainParams.numChains;
    
    if (chainParams.userDefinedTemps == YES)
        {
        return (chainParams.userTemps[id]);
        }
    else
        {
        return (1.0 / (1.0 + chainParams.chainTemp * id));
        }
}


/* Tfree: Free space for partition frequency counter tree */
void Tfree (PFNODE *r)
{
    if (r != NULL)
        {
        if (r->left != NULL)
            Tfree (r->left);
        if (r->right != NULL)
            Tfree (r->right);

        free (r->partition);
        free (r->count);
        free (r);
        }
}


void TouchAllCijks (int chain)
{
    int i;

    for (i=0; i<numCurrentDivisions; i++)
        {
        if (modelSettings[i].nCijkParts >= 1)
            modelSettings[i].upDateCijk = YES;
        }

    return;
    MrBayesPrint ("%d", chain); /* just because I am tired of seeing the unused parameter error msg */
}


void TouchAllPartitions (void)
{
    int i;

    for (i=0; i<numCurrentDivisions; i++)
        {
        modelSettings[i].upDateCl = YES;
        }

    return;
}


void TouchAllTrees (int chain)
{
    int         i, j;
    Tree        *t;
    TreeNode    *p;

    for (i=0; i<numTrees; i++)
        {
        t = GetTreeFromIndex (i, chain, state[chain]);
        for (j=0; j<t->nNodes; j++)
            {
            p = t->allDownPass[j];
            p->upDateCl = YES;
            p->upDateTi = YES;
            }
        }

    for (i=0; i<numCurrentDivisions; i++)
        modelSettings[i].upDateAll = YES;

    return;
}


/* Touch all update flags so we recalculate likelihood from scratch */
void TouchEverything (int chain)
{
    int         i, j;
    Tree        *t;
    TreeNode    *p;

    for (i=0; i<numCurrentDivisions; i++)
        {
        if (modelSettings[i].nCijkParts >= 1)
            modelSettings[i].upDateCijk = YES;
        modelSettings[i].upDateCl = YES;
        modelSettings[i].upDateAll = YES;
        }

    for (i=0; i<numTrees; i++)
        {
        t = GetTreeFromIndex (i, chain, state[chain]);
        for (j=0; j<t->nNodes; j++)
            {
            p = t->allDownPass[j];
            p->upDateCl = YES;
            p->upDateTi = YES;
            }
        }

    return;
}


/*------------------------------------------------------------------
|
|   TreeLength: Calculates the tree length as the sum of the lengths
|       of all the branches. The tree length is the expected number 
|       of character state transformations per character over the 
|       entire phylogenetic tree.
|
-------------------------------------------------------------------*/
MrBFlt TreeLength (Param *param, int chain)
{
    int             i, j;
    MrBFlt          tl;
    Tree            *t;
    TreeNode        *p;

    if (param->paramId == BRLENS_PARSIMONY)
        {
        tl = 0.0;
        for (j=0; j<param->nRelParts; j++)
            tl += modelSettings[param->relParts[j]].parsTreeLength[2*chain+state[chain]];
        }
    else
        {
        /* get tree */
        t = GetTree (param, chain, state[chain]);
        
        /* loop over all branches of the tree */
        tl = 0.0;
        for (i=0; i<t->nNodes; i++)
            {
            p = t->allDownPass[i];
            if (p->anc != NULL)
                {
                if (p->anc->anc == NULL && t->isRooted == NO)
                    {
                    tl += p->length;
                    }
                else
                    {
                    tl += p->length;
                    }
                }
            }
        }
                
    return (tl);
}


/* proportion of ancestral fossils in a FBD tree */
MrBFlt PropAncFossil (Param *param, int chain)
{
    int             i, m, k;
    Tree            *t;
    TreeNode        *p;
 
    t = GetTree (param, chain, state[chain]);
    
    if (t->isRooted == NO)
        return 0.0;
    
    /* count # of tip and ancestral fossils */
    m = k = 0;
    for (i = 0; i < t->nNodes -2; i++)
        {
        p = t->allDownPass[i];
        if (p->left == NULL && p->right == NULL && p->nodeDepth > 0.0)
            {
            if (p->length > 0.0)
                m++;
            else
                k++;
            }
        }
    
    if (k == 0)
        return 0.0;
    else
        return (MrBFlt)k / (m+k);
}

