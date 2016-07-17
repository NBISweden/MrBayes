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
#include "model.h"
#include "sumpt.h"
#include "utils.h"

#ifdef HAVE_LIBREADLINE
#  if defined(HAVE_READLINE_READLINE_H)
#    include <readline/readline.h>
#  elif defined(HAVE_READLINE_H)
#    include <readline.h>
#  endif /* !defined(HAVE_READLINE_H) */
#endif /* HAVE_LIBREADLINE */
#ifdef HAVE_READLINE_HISTORY
#  if defined(HAVE_READLINE_HISTORY_H)
#    include <readline/history.h>
#  elif defined(HAVE_HISTORY_H)
#    include <history.h>
#  endif /* defined(HAVE_READLINE_HISTORY_H) */
#endif /* HAVE_READLINE_HISTORY */

#ifdef HAVE_LIBREADLINE
static char **readline_completion(const char *, int, int);
#endif

/* local prototypes */
int  CommandLine (int argc, char **argv);
void GetTimeSeed (void);
int  InitializeMrBayes (void);
void PrintHeader (void);
int  ReinitializeMrBayes (void);

/* global variables, declared in this file */
BitsLong    bitsLongWithAllBitsSet;      /* BitsLong with all bits set, for bit ops       */
ModelParams defaultModel;                /* Default model; vals set in InitializeMrBayes  */
int         defTaxa;                     /* flag for whether number of taxa is known      */
int         defChars;                    /* flag for whether number of chars is known     */
int         defMatrix;                   /* flag for whether matrix is successfull read   */
int         defPartition;                /* flag for whether character partition is read  */
int         defPairs;                    /* flag for whether pairs are read               */
Doublet     doublet[16];                 /* holds information on states for doublets      */
int         fileNameChanged;             /* has file name been changed ?                  */
RandLong    globalSeed;                  /* seed that is initialized at start up          */
char        **modelIndicatorParams;      /* model indicator params                        */
char        ***modelElementNames;        /* names for component models                    */
int         nBitsInALong;                /* number of bits in a BitsLong                  */
int         nPThreads;                   /* number of pthreads to use                     */
int         numUserTrees;                /* number of defined user trees                  */
int         readComment;                 /* should we read comment (looking for &) ?      */
int         readWord;                    /* should we read word next ?                    */
RandLong    runIDSeed;                   /* seed used only for determining run ID [stamp] */
RandLong    swapSeed;                    /* seed used only for determining which to swap  */
int         userLevel;                   /* user level                                    */
PolyTree    *userTree[MAX_NUM_USERTREES];/* array of user trees                           */
char        workingDir[100];             /* working directory                             */

#if defined (MPI_ENABLED)
int         proc_id;                     /* process ID (0, 1, ..., num_procs-1)                        */
int         num_procs;                   /* number of active processors                                */
MrBFlt      myStateInfo[7];              /* likelihood/prior/heat/ran/moveInfo vals of me              */
MrBFlt      partnerStateInfo[7];         /* likelihood/prior/heat/ran/moveInfo vals of partner         */
#endif

#if defined (FAST_LOG)
CLFlt       scalerValue[400];
CLFlt       logValue[400];
#endif
/* Define to use a log lookup for 4by4 nucleotide data (actually SLOWER than normal code on intel processors) */


int main (int argc, char *argv[])
{
    int i;

#   if defined (MPI_ENABLED)
    int     ierror;
#   endif

#   if defined (WIN_VERSION)
    HANDLE scbh;
    BOOL ok;
    DWORD lastError;
    COORD largestWindow;
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    int currBottom;
    char poltmp[256];

    scbh = GetStdHandle(STD_OUTPUT_HANDLE);
    GetConsoleScreenBufferInfo(scbh, &csbi);
    currBottom      = csbi.srWindow.Bottom;
    largestWindow   = GetLargestConsoleWindowSize(scbh);

    /* Allow for screen buffer 3000 lines long and 140 characters wide */
    csbi.dwSize.Y = 3000;
    csbi.dwSize.X = 140;

    SetConsoleScreenBufferSize(scbh, csbi.dwSize);
    /* Allow for maximum possible screen height */
    csbi.srWindow.Left      = 0; /* no change relative to current value */
    csbi.srWindow.Top       = 0; /* no change relative to current value */
    csbi.srWindow.Right     = 0; /* no change relative to current value */
    csbi.srWindow.Bottom    = largestWindow.Y - currBottom -10; /**/
    ok = SetConsoleWindowInfo(scbh, FALSE, &csbi.srWindow);
    if (ok == FALSE)
        {
        lastError = GetLastError();
        GetConsoleScreenBufferInfo(scbh, &csbi);
        sprintf(poltmp, "\nlastError = %d", lastError);
        // printf (poltmp);
        }
#   endif

    /* calculate the size of a long - used by bit manipulation functions */
    nBitsInALong = sizeof(BitsLong) * 8;
    for (i=0; i<nBitsInALong; i++)
        SetBit(i, &bitsLongWithAllBitsSet);

#   if defined (__MWERKS__) & defined (MAC_VERSION)
    /* Set up interface when using the Metrowerks compiler. This
       should work for either Macintosh or Windows. */
    SIOUXSetTitle("\pMrBayes v3.2");
    SIOUXSettings.fontface         = 0;  /* plain=0; bold=1 */
    SIOUXSettings.setupmenus       = 0;
    SIOUXSettings.autocloseonquit  = 1;
    SIOUXSettings.asktosaveonclose = 0;
    SIOUXSettings.rows             = 60;
    SIOUXSettings.columns          = 90;
#   endif
    
#   if defined (MPI_ENABLED)
    ierror = MPI_Init(&argc, &argv);
    if (ierror != MPI_SUCCESS)
        {
        MrBayesPrint ("%s   Problem initializing MPI\n", spacer);
        exit (1);
        }
    ierror = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (ierror != MPI_SUCCESS)
        {
        MrBayesPrint ("%s   Problem getting the number of processors\n", spacer);
        exit (1);
        }
    ierror = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    if (ierror != MPI_SUCCESS)
        {
        MrBayesPrint ("%s   Problem getting processors rank\n", spacer);
        exit (1);
        }
#   endif
    
#   ifdef HAVE_LIBREADLINE
    rl_attempted_completion_function = readline_completion;
#   endif
    /* Set up parameter table. */
    SetUpParms ();
    
    /* initialize seed using current time */
    GetTimeSeed ();
    
    /* Initialize the variables of the program. */
    InitializeMrBayes ();
    
    /* Print the nifty header. */
    PrintHeader ();
    
    /* Go to the command line, process any arguments passed to the program
       and then wait for input. */
    i = CommandLine (argc, argv);
    
#   if defined (MPI_ENABLED)
    MPI_Finalize();
#   endif
    
    if (i == ERROR)
        return (1);
    else
        return (0);
}


int CommandLine (int argc, char **argv)
{
    int     i, message, nProcessedArgs;
    char    cmdStr[CMD_STRING_LENGTH];
#   ifdef HAVE_LIBREADLINE
#       ifndef MPI_ENABLED
    char    *cmdStrP;
#       endif
#   endif
#   if defined (MPI_ENABLED)
    int     ierror;
#   endif

    for (i=0;i<CMD_STRING_LENGTH;i++) cmdStr[i]='\0';
    
    /* wait for user-input commands */
    nProcessedArgs = 1; /* first argument is program name and needs not be processed */
    if (nProcessedArgs < argc)
        {
        mode = NONINTERACTIVE;  /* when a command is passed into the program, the default is to exit without listening to stdin */
        autoClose = YES;
        autoOverwrite = YES;
        noWarn = YES;
        quitOnError = YES;
        }
    for (;;)
        {
        if (nProcessedArgs < argc) 
            {
            /* we are here only if a command that has been passed
               into the program remains to be processed */
            if (nProcessedArgs == 1 && (strcmp(argv[1],"-i") == 0 || strcmp(argv[1],"-I") == 0))
                {
                mode = INTERACTIVE;
                autoClose = NO;
                autoOverwrite = YES;
                noWarn = NO;
                quitOnError = NO;
                }
            else
                sprintf (cmdStr, "Execute %s", argv[nProcessedArgs]);
            nProcessedArgs++;
            }
        else
            {
            /* first check if we are in noninteractive mode and quit if so */
            if (mode == NONINTERACTIVE)
                {
                MrBayesPrint ("%s   Tasks completed, exiting program because mode is noninteractive\n", spacer);
                MrBayesPrint ("%s   To return control to the command line after completion of file processing, \n", spacer);
                MrBayesPrint ("%s   set mode to interactive with 'mb -i <filename>' (i is for interactive)\n", spacer);
                MrBayesPrint ("%s   or use 'set mode=interactive'\n\n", spacer);
                return (NO_ERROR);
                }
            /* normally, we simply wait at the prompt for a
               user action */
#   if defined (MPI_ENABLED)
            if (proc_id == 0)
                {
                /* do not use readline because OpenMPI does not handle it */
                MrBayesPrint ("MrBayes > ");
                fflush (stdin);
                if (fgets (cmdStr, CMD_STRING_LENGTH - 2, stdin) == NULL)
                    {
                    if (feof(stdin))
                        MrBayesPrint ("%s   End of File encountered on stdin; quitting\n", spacer);
                    else
                        MrBayesPrint ("%s   Could not read command from stdin; quitting\n", spacer);
                    strcpy (cmdStr,"quit;\n");
                    }
                }
            ierror = MPI_Bcast (&cmdStr, CMD_STRING_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
            if (ierror != MPI_SUCCESS)
                {
                MrBayesPrint ("%s   Problem broadcasting command string\n", spacer);
                }
#   else
#       ifdef HAVE_LIBREADLINE
            cmdStrP = readline("MrBayes > ");
            if (cmdStrP!=NULL) 
                    {
                    strncpy (cmdStr,cmdStrP,CMD_STRING_LENGTH - 2);
                    if (*cmdStrP) 
                        add_history (cmdStrP);
                    free (cmdStrP);
                    }
            else /* fall through to if (feof(stdin))..*/
#       else
            MrBayesPrint ("MrBayes > ");
            fflush (stdin);
            if (fgets (cmdStr, CMD_STRING_LENGTH - 2, stdin) == NULL)
#       endif
                {
                if (feof(stdin))
                    MrBayesPrint ("%s   End of File encountered on stdin; quitting\n", spacer);
                else
                    MrBayesPrint ("%s   Could not read command from stdin; quitting\n", spacer);
                strcpy (cmdStr,"quit;\n");
                }
#   endif
            }
        i = 0;
        while (cmdStr[i] != '\0' && cmdStr[i] != '\n')
            i++;
        cmdStr[i++] = ';';
        cmdStr[i] = '\0';
        MrBayesPrint ("\n");
        if (cmdStr[0] != ';')
            {
            /* check that all characters in the string are valid */
            if (CheckStringValidity (cmdStr) == ERROR)
                {
                MrBayesPrint ("   Unknown character in command string\n\n");
                }
            else
                {
                expecting = Expecting(COMMAND);
                message = ParseCommand (cmdStr);

                if (message == NO_ERROR_QUIT)
                    return (NO_ERROR);

                if (message == ERROR && quitOnError == YES)
                    {
                    MrBayesPrint ("%s   Will exit with signal 1 (error) because quitonerror is set to yes\n", spacer);
                    MrBayesPrint ("%s   If you want control to be returned to the command line on error,\n", spacer);
                    MrBayesPrint ("%s   use 'mb -i <filename>' (i is for interactive) or use 'set quitonerror=no'\n\n", spacer);
                    return (ERROR);
                    }

#   if defined (MPI_ENABLED)
                ierror = MPI_Barrier (MPI_COMM_WORLD);
                if (ierror != MPI_SUCCESS)
                    {
                    MrBayesPrint ("%s   Problem at command barrier\n", spacer);
                    }
#   endif

                MrBayesPrint ("\n");
                }
            }
        }
}


#ifdef HAVE_LIBREADLINE
extern char *command_generator(const char *text, int state);

char **readline_completion (const char *text, int start, int stop)
{
    char **matches = (char **) NULL;

#   ifdef COMPLETIONMATCHES
    if (start == 0)
            matches = rl_completion_matches (text, command_generator);
#   endif

    return (matches);   
}
#endif


void GetTimeSeed (void)
{
    time_t      curTime;

#   if defined (MPI_ENABLED)
    int         ierror;
    
    if (proc_id == 0)
        {
        curTime = time(NULL);
        globalSeed  = (RandLong)curTime;
        if (globalSeed < 0)
            globalSeed = -globalSeed;
        }
    ierror = MPI_Bcast(&globalSeed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (ierror != MPI_SUCCESS)
        {
        MrBayesPrint ("%s   Problem broadcasting seed\n", spacer);
        }
        
    if (proc_id == 0)
        {
        /* Note: swapSeed will often be same as globalSeed */
        curTime = time(NULL);
        swapSeed  = (RandLong)curTime;
        if (swapSeed < 0)
            swapSeed = -swapSeed;
        }
    ierror = MPI_Bcast(&swapSeed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (ierror != MPI_SUCCESS)
        {
        MrBayesPrint ("%s   Problem broadcasting swap seed\n", spacer);
        }

    if (proc_id == 0)
        {
        /* Note: runIDSeed will often be same as globalSeed */
        curTime = time(NULL);
        runIDSeed  = (RandLong)curTime;
        if (runIDSeed < 0)
            runIDSeed = -runIDSeed;
        }
    ierror = MPI_Bcast(&runIDSeed, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    if (ierror != MPI_SUCCESS)
        {
        MrBayesPrint ("%s   Problem broadcasting run ID seed\n", spacer);
        }

#   else
    curTime = time(NULL);
    globalSeed  = (RandLong)curTime;
    if (globalSeed < 0)
        globalSeed = -globalSeed;
        
    /* Note: swapSeed will often be the same as globalSeed */
    curTime = time(NULL);
    swapSeed  = (RandLong)curTime;
    if (swapSeed < 0)
        swapSeed = -swapSeed;

    /* Note: runIDSeed will often be the same as globalSeed */
    curTime = time(NULL);
    runIDSeed  = (RandLong)curTime;
    if (runIDSeed < 0)
        runIDSeed = -runIDSeed;
        
#   endif
}


int InitializeMrBayes (void)
{
    /* this function initializes the program; only call it at the start of execution */
    
    int     i, j, growthFxn[6];

    nBitsInALong         = sizeof(BitsLong) * 8;     /* global variable: number of bits in a BitsLong */
    userLevel            = STANDARD_USER;            /* default user level                            */

    readWord             = NO;                       /* should we read a word next ?                  */
    readComment          = NO;                       /* should we read comments? (used by tree cmd)   */
    fileNameChanged      = NO;                       /* file name changed ? (used by a few commands)  */
    echoMB               = YES;                      /* flag used by Manual to control printing       */

#   if defined (MPI_ENABLED)
    sprintf (manFileName, "commref_mb%sp.txt", VERSION_NUMBER);  /* name of command reference file    */
#   else
    sprintf (manFileName, "commref_mb%s.txt", VERSION_NUMBER);   /* name of command reference file    */
#   endif

    for (i=0; i<NUM_ALLOCS; i++)                     /* set allocated memory to NO                    */
        memAllocs[i] = NO;              
    logToFile = NO;                                  /* should screen output be logged to a file      */
    strcpy(logFileName, "log.out");                  /* name of the log file                          */
    logFileFp = NULL;                                /* file pointer to log file                      */
    replaceLogFile = YES;                            /* should logfile be replace/appended to         */
    autoClose = NO;                                  /* set default autoclose                         */
    autoOverwrite = YES;                             /* set default autoOverwrite                     */
    noWarn = NO;                                     /* set default                                   */
    quitOnError = NO;                                /* set default quitOnError                       */
    inferAncStates = NO;                             /* set default inferAncStates                    */
    inferSiteOmegas = NO;                            /* set default inferSiteOmegas                   */
    inferSiteRates = NO;                             /* set default inferSiteRates                    */
    inferPosSel = NO;                                /* set default inferPosSel                       */
    inComment = NO;                                  /* not in comment                                */
    numComments = 0;                                 /* no comments encountered yet                   */
    mode = INTERACTIVE;                              /* set default mode                              */
    numOpenExeFiles = 0;                             /* no execute files open yet                     */
    scientific = YES;                                /* print to file using scientific format?        */
    precision = 6;                                   /* set default precision                         */
    showmovesParams.allavailable = NO;               /* do not show all available moves               */
    strcpy(workingDir,"");                           /* working directory                             */
#   if defined (BEAGLE_ENABLED)
#       if defined (WIN_VERSION)
    tryToUseBEAGLE = NO;                             /* try to use the BEAGLE library (NO until SSE code works in Win) */
#       else
    tryToUseBEAGLE = NO;                             /* try to use the BEAGLE library if not Win (NO untill SSE single prec. works) */
#       endif
    beagleScalingScheme = MB_BEAGLE_SCALE_ALWAYS;    /* use BEAGLE always scaling                     */
    beagleFlags = BEAGLE_FLAG_PROCESSOR_CPU;         /* default to generic CPU                        */
    beagleResourceNumber = 99;                       /* default to auto-resource selection            */
    // SSE instructions do not work in Windows environment
    // beagleFlags |= BEAGLE_FLAG_VECTOR_SSE;        /* default to SSE code                           */
    beagleResource = NULL;
    beagleResourceCount = 0;                         /* default has no list */
    beagleInstanceCount = 0;                         /* no BEAGLE instances */
    beagleScalingFrequency = 1000;  
#   endif
#   if defined (THREADS_ENABLED)
    tryToUseThreads = NO;                            /* try to use pthread with BEAGLE library        */
#   endif

    /* set the proposal information */
    SetUpMoveTypes ();

    /* Set up rates for standard amino acid models, in case we need them. */
    if (SetAARates () == ERROR)
        return (ERROR);
           
    /* set up doublet information */
    doublet[ 0].first  = 1;   doublet[ 0].second = 1;
    doublet[ 1].first  = 1;   doublet[ 1].second = 2;
    doublet[ 2].first  = 1;   doublet[ 2].second = 4;
    doublet[ 3].first  = 1;   doublet[ 3].second = 8;
    doublet[ 4].first  = 2;   doublet[ 4].second = 1;
    doublet[ 5].first  = 2;   doublet[ 5].second = 2;
    doublet[ 6].first  = 2;   doublet[ 6].second = 4;
    doublet[ 7].first  = 2;   doublet[ 7].second = 8;
    doublet[ 8].first  = 4;   doublet[ 8].second = 1;
    doublet[ 9].first  = 4;   doublet[ 9].second = 2;
    doublet[10].first  = 4;   doublet[10].second = 4;
    doublet[11].first  = 4;   doublet[11].second = 8;
    doublet[12].first  = 8;   doublet[12].second = 1;
    doublet[13].first  = 8;   doublet[13].second = 2;
    doublet[14].first  = 8;   doublet[14].second = 4;
    doublet[15].first  = 8;   doublet[15].second = 8;

#   if defined (FAST_LOG)
    /* set up log table */
    for (i=0; i<400; i++)
        {
        scalerValue[i] = (CLFlt) ldexp (1.0, 1-i);  /* offset 1 needed to deal with scaler == 1.0 */
        logValue[i] = (CLFlt) log (scalerValue[i]);     
        }
#   endif

    /* user trees */
    for (i=0; i<MAX_NUM_USERTREES; i++)
        userTree[i] = NULL;

    /* parameter values */
    paramValues = NULL;
    intValues = NULL;
    
    /* Prior model settings */
    defaultModel.dataType = DNA;                        /* datatype                                     */
    defaultModel.coding = 0;                            /* ascertainment bias                           */
    strcpy(defaultModel.codingString, "All");           /* ascertainment bias string                    */
    strcpy(defaultModel.nucModel, "4by4");              /* nucleotide model                             */
    strcpy(defaultModel.nst, "1");                      /* number of substitution types                 */
    strcpy(defaultModel.aaModelPr, "Fixed");            /* amino acid model prior                       */
    for (i=0; i<10; i++)
        defaultModel.aaModelPrProbs[i] = 0.0;
    strcpy(defaultModel.aaModel, "Poisson");            /* amino acid model                             */
    strcpy(defaultModel.parsModel, "No");               /* do not use parsimony model                   */
    strcpy(defaultModel.geneticCode, "Universal");      /* genetic code                                 */
    strcpy(defaultModel.ploidy, "Diploid");             /* ploidy level                                 */
    strcpy(defaultModel.omegaVar, "Equal");             /* omega variation                              */
    strcpy(defaultModel.ratesModel, "Equal");           /* rates across sites model                     */
    defaultModel.numGammaCats = 4;                      /* number of categories for gamma approximation */
    defaultModel.numLnormCats = 4;                      /* number of categories for lnorm approximation */
    defaultModel.numMixtCats = 4;                       /* number of components in rate mixture         */
    strcpy(defaultModel.useGibbs,"No");                 /* do not use Gibbs sampling of rate cats by default */
    defaultModel.gibbsFreq = 100;                       /* default Gibbs sampling frequency of rate cats*/
    defaultModel.numBetaCats = 5;                       /* number of categories for beta approximation  */
    strcpy(defaultModel.covarionModel, "No");           /* use covarion model? (yes/no)                 */
    strcpy(defaultModel.augmentData, "No");             /* should data be augmented                     */
    strcpy(defaultModel.tRatioPr, "Beta");              /* prior for ti/tv rate ratio                   */
    defaultModel.tRatioFix = 1.0;
    defaultModel.tRatioDir[0] = 1.0;
    defaultModel.tRatioDir[1] = 1.0;
    strcpy(defaultModel.revMatPr, "Dirichlet");         /* prior for GTR model (nucleotides)            */
    for (i=0; i<6; i++)
        {
        defaultModel.revMatFix[i] = 1.0;
        defaultModel.revMatDir[i] = 1.0;
        }
    defaultModel.revMatSymDir = 1.0;                    /* default prior for GTR mixed model            */
    strcpy (defaultModel.aaRevMatPr, "Dirichlet");      /* prior for GTR model (proteins)               */
    for (i=0; i<190; i++)
        {
        defaultModel.aaRevMatFix[i] = 1.0;
        defaultModel.aaRevMatDir[i] = 1.0;
        }
    strcpy(defaultModel.omegaPr, "Dirichlet");          /* prior for omega                              */
    defaultModel.omegaFix = 1.0;
    defaultModel.omegaDir[0] = 1.0;
    defaultModel.omegaDir[1] = 1.0;
    strcpy(defaultModel.ny98omega1pr, "Beta");          /* prior for class 1 omega (Ny98 model)         */
    defaultModel.ny98omega1Fixed = 0.1;
    defaultModel.ny98omega1Beta[0] = 1.0;
    defaultModel.ny98omega1Beta[1] = 1.0;
    strcpy(defaultModel.ny98omega3pr, "Exponential");   /* prior for class 3 omega (Ny98 model)        */
    defaultModel.ny98omega3Fixed = 2.0;
    defaultModel.ny98omega3Uni[0] = 1.0;
    defaultModel.ny98omega3Uni[1] = 50.0;
    defaultModel.ny98omega3Exp = 1.0;
    strcpy(defaultModel.m3omegapr, "Exponential");      /* prior for all three omegas (M3 model)        */
    defaultModel.m3omegaFixed[0] = 0.1;
    defaultModel.m3omegaFixed[1] = 1.0;
    defaultModel.m3omegaFixed[2] = 2.0;
    strcpy(defaultModel.m10betapr, "Uniform");          /* prior for omega variation (M10 model)        */
    strcpy(defaultModel.m10gammapr, "Uniform");
    defaultModel.m10betaUni[0] = 0.0;
    defaultModel.m10betaUni[1] = 20.0;
    defaultModel.m10betaExp = 1.0;
    defaultModel.m10betaFix[0] = 1.0;
    defaultModel.m10betaFix[1] = 1.0;
    defaultModel.m10gammaUni[0] = 0.0;
    defaultModel.m10gammaUni[1] = 20.0;
    defaultModel.m10gammaExp = 1.0;
    defaultModel.m10gammaFix[0] = 1.0;
    defaultModel.m10gammaFix[1] = 1.0;
    defaultModel.numM10GammaCats = 4;
    defaultModel.numM10BetaCats = 4;
    strcpy(defaultModel.codonCatFreqPr, "Dirichlet");   /* prior for selection cat frequencies         */
    defaultModel.codonCatFreqFix[0] = 1.0/3.0;
    defaultModel.codonCatFreqFix[1] = 1.0/3.0;
    defaultModel.codonCatFreqFix[2] = 1.0/3.0;
    defaultModel.codonCatDir[0] = 1.0;
    defaultModel.codonCatDir[1] = 1.0;
    defaultModel.codonCatDir[2] = 1.0;
    strcpy(defaultModel.stateFreqPr, "Dirichlet");      /* prior for character state frequencies        */
    strcpy(defaultModel.stateFreqsFixType, "Equal");
    for (i=0; i<200; i++)
        {
        defaultModel.stateFreqsFix[i] = 0.0;   
        defaultModel.stateFreqsDir[i] = 1.0;
        }    
    defaultModel.numDirParams = 0;
    strcpy(defaultModel.shapePr, "Exponential");        /* prior for gamma/lnorm shape parameter        */
    defaultModel.shapeFix = 0.5;
    defaultModel.shapeUni[0] = MIN_SHAPE_PARAM;
    defaultModel.shapeUni[1] = MAX_SHAPE_PARAM;
    defaultModel.shapeExp = 1.0;
    strcpy(defaultModel.pInvarPr, "Uniform");           /* prior for proportion of invariable sites     */
    defaultModel.pInvarFix = 0.1;
    defaultModel.pInvarUni[0] = 0.0;
    defaultModel.pInvarUni[1] = 1.0;
    strcpy(defaultModel.adGammaCorPr, "Uniform");       /* prior for correlation param of adGamma model */
    defaultModel.corrFix = 0.0;
    defaultModel.corrUni[0] = -1.0;
    defaultModel.corrUni[1] = 1.0;
    strcpy(defaultModel.covSwitchPr, "Uniform");        /* prior for switching rates of covarion model  */
    defaultModel.covswitchFix[0] = 1.0;
    defaultModel.covswitchFix[1] = 1.0;
    defaultModel.covswitchUni[0] = 0.0;
    defaultModel.covswitchUni[1] = 100.0;
    defaultModel.covswitchExp = 1.0;
    strcpy(defaultModel.symPiPr, "Fixed");              /* prior for pi when unidentifiable states used */
    defaultModel.symBetaFix = -1.0;
    defaultModel.symBetaUni[0] = 0.0;
    defaultModel.symBetaUni[1] = 20.0;
    defaultModel.symBetaExp = 2;
    strcpy(defaultModel.brownCorPr, "Fixed");           /* prior on correlation of brownian model       */
    defaultModel.brownCorrFix = 0.0;
    defaultModel.brownCorrUni[0] = -1.0;
    defaultModel.brownCorrUni[1] = 1.0;
    strcpy(defaultModel.brownScalesPr, "Gammamean");    /* prior on scales of brownian model            */
    defaultModel.brownScalesFix = 10.0;
    defaultModel.brownScalesUni[0] = 0.0;
    defaultModel.brownScalesUni[1] = 100.0;
    defaultModel.brownScalesGamma[0] = 1.0;
    defaultModel.brownScalesGamma[1] = 10.0;
    defaultModel.brownScalesGammaMean = 10.0;
    strcpy(defaultModel.topologyPr, "Uniform");         /* prior for tree topology                      */
    defaultModel.topologyFix = -1;                      /* user tree index to use for fixed topology    */
    defaultModel.activeConstraints = NULL;              /* which constraints are active                 */
    strcpy(defaultModel.brlensPr, "Unconstrained");     /* prior on branch lengths                      */
    defaultModel.brlensFix = -1;                        /* user tree index to use for fixed brlens      */
    defaultModel.brlensUni[0] = BRLENS_MIN;
    defaultModel.brlensUni[1] = 10.0;
    defaultModel.brlensExp    = 10.0;
    defaultModel.brlens2Exp[0]= 100.0;                  /* 1st param of twoExp prior (for internal branches) */
    defaultModel.brlens2Exp[1]= 10.0;                   /* 2nd param of twoExp prior (for external branches) */
    defaultModel.brlensDir[0] = 1.0;                    /* 1st param of GammaDir prior   */
//  defaultModel.brlensDir[0] = 3.0;                    /* 1st param of invGamDir prior  */
    defaultModel.brlensDir[1] = 0.1;                    /* 2nd param of GammaDir prior   */
//  defaultModel.brlensDir[1] = 20.0;                   /* 2nd param of invGamDir prior  */
    defaultModel.brlensDir[2] = 1.0;                    /* 3rd param of Dirichlet priors */
    defaultModel.brlensDir[3] = 1.0;                    /* 4th param of Dirichlet priors */
    
    strcpy(defaultModel.unconstrainedPr, "GammaDir");   /* prior on branches if unconstrained           */
    strcpy(defaultModel.clockPr, "Uniform");            /* prior on branch lengths if clock enforced    */
    defaultModel.treeAgePr.prior = standardGamma;       /* calibration prior on tree age */
    strcpy(defaultModel.treeAgePr.name, "Gamma(1.00,1.00)");
    defaultModel.treeAgePr.priorParams[0] = 1.0;
    defaultModel.treeAgePr.priorParams[1] = 1.0;
    defaultModel.treeAgePr.priorParams[2] = -1.0;
    defaultModel.treeAgePr.LnPriorProb = &LnPriorProbGamma_Param_Mean_Sd;
    defaultModel.treeAgePr.LnPriorRatio = &LnProbRatioGamma_Param_Mean_Sd;
    defaultModel.treeAgePr.min = 0.0;
    defaultModel.treeAgePr.max = POS_INFINITY;
    strcpy(defaultModel.clockRatePr, "Fixed");          /* prior on base subst. rate for clock trees    */
    defaultModel.clockRateNormal[0] = 1.0;
    defaultModel.clockRateNormal[1] = 1.0;
    defaultModel.clockRateLognormal[0] = 0.0;           /* mean 0.0 on log scale corresponds to mean rate 1.0 */
    defaultModel.clockRateLognormal[1] = 0.7;           /* double or half the rate in one standard deviation  */
    defaultModel.clockRateGamma[0] = 1.0;
    defaultModel.clockRateGamma[1] = 1.0;
    defaultModel.clockRateExp = 1.0;
    defaultModel.clockRateFix = 1.0;
    strcpy(defaultModel.speciationPr, "Exponential");   /* prior on speciation rate (net diversification) */
    defaultModel.speciationFix = 0.1;
    defaultModel.speciationUni[0] = 0.0;
    defaultModel.speciationUni[1] = 10.0;
    defaultModel.speciationExp = 10.0;
    strcpy(defaultModel.extinctionPr, "Beta");          /* prior on extinction rate (turnover)          */
    defaultModel.extinctionFix = 0.5;
    defaultModel.extinctionBeta[0] = 1;
    defaultModel.extinctionBeta[1] = 1;
    strcpy(defaultModel.fossilizationPr, "Beta");       /* prior on fossilization rate (sampling proportion) */
    defaultModel.fossilizationFix = 0.5;
    defaultModel.fossilizationBeta[0] = 1;
    defaultModel.fossilizationBeta[1] = 1;
    strcpy(defaultModel.sampleStrat, "Random");         /* taxon sampling strategy                      */
    defaultModel.sampleProb = 1.0;                      /* extant taxon sampling fraction               */
    defaultModel.birthRateShiftNum = 0;                 /* number of birth rate shifts                  */
    defaultModel.deathRateShiftNum = 0;                 /* number of death rate shifts                  */
    defaultModel.fossilSamplingNum = 0;                 /* number of fossil sampling rate shifts / slice sampling events */

    strcpy(defaultModel.popSizePr, "Gamma");            /* prior on coalescence population size         */
    defaultModel.popSizeFix = 100.0;                    /* N_e = 100 */
    defaultModel.popSizeUni[0] = 0.0;
    defaultModel.popSizeUni[1] = 1000.0;
    defaultModel.popSizeNormal[0] = 100.0;
    defaultModel.popSizeNormal[1] = 30.0;
    defaultModel.popSizeLognormal[0] = 4.6;             /* mean on log scale corresponds to N_e = 100.0 */
    defaultModel.popSizeLognormal[1] = 0.4;             /* factor 10 in one standard deviation          */
    defaultModel.popSizeGamma[0] = 1.0;                 /* exponential with mean 1/10 = 0.1 */
    defaultModel.popSizeGamma[1] = 10.0;
    strcpy(defaultModel.popVarPr, "Equal");             /* prior on pop. size variation across tree     */
    strcpy(defaultModel.growthPr, "Fixed");             /* prior on coalescence growth rate prior       */
    defaultModel.growthFix = 0.0;
    defaultModel.growthUni[0] = 0.0;
    defaultModel.growthUni[1] = 100.0;
    defaultModel.growthExp = 1.0;
    defaultModel.growthNorm[0] = 0.0;
    defaultModel.growthNorm[1] = 1.0;
    strcpy(defaultModel.nodeAgePr, "Unconstrained");    /* prior on node depths                       */
    strcpy(defaultModel.clockVarPr, "Strict");          /* prior on clock rate variation              */
    strcpy(defaultModel.cppRatePr, "Exponential") ;     /* prior on rate of CPP for relaxed clock     */
    defaultModel.cppRateExp = 0.1;
    defaultModel.cppRateFix = 1.0;
    strcpy(defaultModel.cppMultDevPr, "Fixed");         /* prior on standard dev. of lognormal of rate multipliers of CPP rel clock */
    defaultModel.cppMultDevFix = 0.4;
    strcpy(defaultModel.tk02varPr, "Exponential");      /* prior on nu parameter for BM rel clock     */
    defaultModel.tk02varExp = 1.0;
    defaultModel.tk02varFix = 1.0;
    defaultModel.tk02varUni[0] = 0.0;
    defaultModel.tk02varUni[1] = 5.0;
    strcpy(defaultModel.igrvarPr, "Exponential");       /* prior on variance increase parameter for IGR rel clock */
    defaultModel.igrvarExp = 10.0;
    defaultModel.igrvarFix = 0.1;
    defaultModel.igrvarUni[0] = 0.0;
    defaultModel.igrvarUni[1] = 0.5;
    strcpy(defaultModel.mixedvarPr, "Exponential");     /* prior on var parameter for mixed rel clock */
    defaultModel.mixedvarExp = 10.0;
    defaultModel.mixedvarFix = 0.1;
    defaultModel.mixedvarUni[0] = 0.0;
    defaultModel.mixedvarUni[1] = 5.0;
    strcpy(defaultModel.ratePr, "Fixed");               /* prior on rate for a partition              */
    defaultModel.ratePrDir = 1.0;
    strcpy(defaultModel.generatePr, "Fixed");           /* prior on rate for a gene (multispecies coalescent) */
    defaultModel.generatePrDir = 1.0;

    defaultModel.nStates = 4;                           /* number of states for partition             */

    /* Report format settings */
    strcpy(defaultModel.tratioFormat, "Ratio");         /* default format for tratio                  */
    strcpy(defaultModel.revmatFormat, "Dirichlet");     /* default format for revmat                  */
    strcpy(defaultModel.ratemultFormat, "Scaled");      /* default format for ratemult                */
    strcpy(defaultModel.treeFormat, "Brlens");          /* default format for trees                   */
    strcpy(defaultModel.inferAncStates, "No");          /* do not infer ancestral states              */
    strcpy(defaultModel.inferPosSel, "No");             /* do not infer positive selection            */
    strcpy(defaultModel.inferSiteOmegas, "No");         /* do not infer site omega vals               */
    strcpy(defaultModel.inferSiteRates, "No");          /* do not infer site rates                    */

    /* Allocate and initialize model indicator parameter names */
    modelIndicatorParams = (char **) SafeCalloc (3, sizeof (char *));
    modelIndicatorParams[0] = "aamodel";
    modelIndicatorParams[1] = "gtrsubmodel";
    modelIndicatorParams[2] = "";

    /* Aamodel */
    modelElementNames = (char ***) SafeCalloc (3, sizeof (char **));
    modelElementNames[0] = (char **) SafeCalloc (12, sizeof (char *));
    modelElementNames[0][0]  = "Poisson";
    modelElementNames[0][1]  = "Jones";
    modelElementNames[0][2]  = "Dayhoff";
    modelElementNames[0][3]  = "Mtrev";
    modelElementNames[0][4]  = "Mtmam";
    modelElementNames[0][5]  = "Wag";
    modelElementNames[0][6]  = "Rtrev";
    modelElementNames[0][7]  = "Cprev";
    modelElementNames[0][8]  = "Vt";
    modelElementNames[0][9]  = "Blosum";
    modelElementNames[0][10] = "LG";
    modelElementNames[0][11] = "";

    /* Gtrsubmodel */
    modelElementNames[1] = (char **) SafeCalloc (204, sizeof (char *));
    for (i=0; i<203; i++)
        {
        modelElementNames[1][i]  = (char *) SafeCalloc (7, sizeof (char));
        FromIndexToGrowthFxn(i, growthFxn);
        for (j=0; j<6; j++)
            modelElementNames[1][i][j] = '1' + growthFxn[j];
        modelElementNames[1][i][j] = '\0';
        }
    modelElementNames[1][203]  = "";

    /* Termination */
    modelElementNames[2]    = (char **) SafeCalloc (1, sizeof(char *));
    modelElementNames[2][0] = "";

    /* initialize user trees */
    for (i=0; i<MAX_NUM_USERTREES; i++)
        userTree[i] = NULL;
    numUserTrees = 0;

    /* Reset translate table */
    ResetTranslateTable();

    /* finally reset everything dependent on a matrix being defined */
    return (ReinitializeMrBayes ());
    
}


void PrintHeader (void)
{
    MrBayesPrint ("\n\n");
    MrBayesPrint ("                            MrBayes v%s %s\n\n", VERSION_NUMBER, HOST_CPU);
    MrBayesPrint ("                      (Bayesian Analysis of Phylogeny)\n\n");
#   if defined (MPI_ENABLED)
    MrBayesPrint ("                             (Parallel version)\n");
    MrBayesPrint ("                         (%d processors available)\n\n", num_procs);
#   endif
    MrBayesPrint ("              Distributed under the GNU General Public License\n\n\n");
    MrBayesPrint ("               Type \"help\" or \"help <command>\" for information\n");
    MrBayesPrint ("                     on the commands that are available.\n\n");
    MrBayesPrint ("                   Type \"about\" for authorship and general\n");
    MrBayesPrint ("                       information about the program.\n\n\n");
}


int ReinitializeMrBayes (void)
{
    /* this function resets everything dependent on a matrix */
    
    int             i;
    
    /* reinitialize indentation */
    strcpy (spacer, "");                             /* holds blanks for indentation                  */

    /* reset all taxa flags */
    ResetTaxaFlags();

    /* reset all characters flags */
    ResetCharacterFlags();

    /* chain parameters */
    chainParams.numGen = 1000000;                    /* number of MCMC cycles                         */
    chainParams.sampleFreq = 500;                    /* frequency to sample chain                     */
    chainParams.printFreq = 1000;                    /* frequency to print chain                      */
    chainParams.swapFreq = 1;                        /* frequency of attempting swap of states        */
    chainParams.numSwaps = 1;                        /* number of swaps to try each time              */
    chainParams.isSS = NO;
    chainParams.startFromPriorSS = NO;
    chainParams.numStepsSS = 50;
    chainParams.burninSS = -1;
    chainParams.alphaSS = 0.4;
    chainParams.backupCheckSS = 0;
    chainParams.mcmcDiagn = YES;                     /* write MCMC diagnostics to file ?              */
    chainParams.diagnFreq = 5000;                    /* diagnostics frequency                         */
    chainParams.minPartFreq = 0.10;                  /* min partition frequency for diagnostics       */
    chainParams.allChains = NO;                      /* calculate diagnostics for all chains ?        */
    chainParams.allComps = NO;                       /* do not calc diagn for all run comparisons     */
    chainParams.relativeBurnin = YES;                /* use relative burnin?                          */
    chainParams.burninFraction = 0.25;               /* default burnin fraction                       */
    chainParams.stopRule = NO;                       /* should stopping rule be used?                 */
    chainParams.stopVal = 0.05;                      /* convergence diagnostic value to reach         */
    chainParams.numRuns = 2;                         /* number of runs                                */
    chainParams.numChains = 4;                       /* number of chains                              */
    chainParams.chainTemp = 0.1;                     /* chain temperature                             */
    chainParams.redirect = NO;                       /* should printf be to stdout                    */
    strcpy(chainParams.chainFileName, "temp");       /* chain file name for output                    */
    chainParams.chainBurnIn = 0;                     /* chain burn in length                          */
    chainParams.numStartPerts = 0;                   /* number of perturbations to starting tree      */
    strcpy(chainParams.startTree, "Current");        /* starting tree for chain (random/current)      */
    strcpy(chainParams.startParams, "Current");      /* starting params for chain (reset/current)     */
    chainParams.saveBrlens = YES;                    /* should branch lengths be saved                */
    chainParams.weightScheme[0] = 0.0;               /* percent chars to decrease in weight           */
    chainParams.weightScheme[1] = 0.0;               /* percent chars to increase in weight           */
    chainParams.weightScheme[2] = 1.0;               /* weight increment                              */
    chainParams.calcPbf = NO;                        /* should we calculate the pseudo-BF?            */
    chainParams.pbfInitBurnin = 100000;              /* initial burnin for pseudo BF                  */
    chainParams.pbfSampleFreq = 10;                  /* sample frequency for pseudo BF                */
    chainParams.pbfSampleTime = 2000;                /* how many cycles to calcualate site prob.      */
    chainParams.pbfSampleBurnin = 2000;              /* burnin period for each site for pseudo BF     */
    chainParams.userDefinedTemps = NO;               /* should we use the users temperatures?         */
    for (i=0; i<MAX_CHAINS; i++)
        chainParams.userTemps[i] = 1.0;              /* user-defined chain temperatures               */
    chainParams.swapAdjacentOnly = NO;               /* swap only adjacent temperatures               */
    chainParams.printMax = 8;                        /* maximum number of chains to print to screen   */
    chainParams.printAll = YES;                      /* whether to print heated chains                */
    chainParams.treeList = NULL;                     /* vector of tree lists for saving trees         */
    chainParams.saveTrees = NO;                      /* save tree samples for later removal?          */
    chainParams.tFilePos = NULL;                     /* position for reading trees for removal        */
    chainParams.runWithData = YES;                   /* whether to run with data                      */
    chainParams.orderTaxa = NO;                      /* should taxa be ordered in output trees?       */
    chainParams.append = NO;                         /* append to previous analysis?                  */
    chainParams.autotune = YES;                      /* autotune?                                     */
    chainParams.tuneFreq = 100;                      /* autotuning frequency                          */
    chainParams.checkPoint = YES;                    /* should we checkpoint the run?                 */
    chainParams.checkFreq = 2000;                    /* check-pointing frequency                      */
    chainParams.diagnStat = AVGSTDDEV;               /* mcmc diagnostic to use                        */

    /* sumt parameters */
    strcpy(sumtParams.sumtFileName, "temp");         /* input name for sumt command                   */
    strcpy(sumtParams.sumtConType, "Halfcompat");    /* type of consensus tree output                 */
    sumtParams.calcTreeprobs = YES;                  /* should individual tree probs be calculated    */
    sumtParams.showSumtTrees = NO;                   /* should the individual tree probs be shown     */
    sumtParams.numTrees = 1;                         /* number of trees to summarize                  */
    sumtParams.numRuns = 2;                          /* number of analyses to summarize               */
    sumtParams.orderTaxa = YES;                      /* order taxa in trees ?                         */
    sumtParams.minPartFreq = 0.10;                   /* minimum part. freq. for overall diagnostics   */
    sumtParams.table = YES;                          /* display table of part. freq.?                 */
    sumtParams.summary = YES;                        /* display overall diagnostics?                  */
    sumtParams.showConsensus = YES;                  /* display consensus tree(s)?                    */
    sumtParams.consensusFormat = FIGTREE;            /* format of consensus tree                      */
    strcpy (sumtParams.sumtOutfile, "temp");         /* output name for sumt command                  */
    sumtParams.HPD = YES;                            /* use Highest Posterior Density?                */
    sumtParams.saveBrParams = NO;                    /* save branch/node params?                      */
    sumtParams.minBrParamFreq = 0.50;                /* threshold for printing branch params to file  */

    /* sump parameters */
    strcpy(sumpParams.sumpFileName, "temp");         /* input name for sump command                   */
    strcpy (sumpParams.sumpOutfile, "temp");         /* output name for sump command                  */
    sumpParams.numRuns = 2;                          /* number of analyses to summarize               */
    sumpParams.HPD = YES;                            /* use Highest Posterior Density?                */
    sumpParams.minProb = 0.05;                       /* min. prob. of models to include in summary    */

    /* sumss parameters */
    sumssParams.numRuns= 2;                          /* number of independent analyses to summarize   */
    sumssParams.allRuns = YES;                       /* should data for all runs be printed (yes/no)? */
    sumssParams.stepToPlot = 0;                      /* Which step to plot in the step plot, 0 means burnin */
    sumssParams.askForMorePlots = YES;               /* Should user be asked to plot for different discardfraction (y/n)?  */
    sumssParams.discardFraction = 0.8;               /* Proportion of samples discarded when ploting step plot.*/
    sumssParams.smoothing = 0;                       /* An integer indicating number of neighbors to average over
                                                        when dooing smoothing of curvs on plots */
    /* comparetree parameters */
    strcpy(comptreeParams.comptFileName1, "temp.t"); /* input name for comparetree command            */
    strcpy(comptreeParams.comptFileName2, "temp.t"); /* input name for comparetree command            */
    strcpy(comptreeParams.comptOutfile, "temp.comp");/* output name for comparetree command           */
    comptreeParams.minPartFreq = 0.0;                /* minimum frequency of partitions to include    */

    /* plot parameters */
    strcpy(plotParams.plotFileName, "temp.p");       /* input name for plot command                   */
    strcpy(plotParams.parameter, "lnL");             /* plotted parameter plot command                */
    strcpy(plotParams.match, "Perfect");             /* matching for plot command                     */
    
    return (NO_ERROR);
}

