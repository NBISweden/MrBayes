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
#include "mbbeagle.h"
#include "model.h"
#include "mcmc.h"
#include "sumpt.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

#define NUMCOMMANDS                     62    /* The total number of commands in the program  */
#define NUMPARAMS                       278   /* The total number of parameters  */
#define PARAM(i, s, f, l)               p->string = s;    \
                                        p->fp = f;        \
                                        p->valueList = l; \
                                        p++;
#define HIDE                            0
#define SHOW                            1

/* Debugging options */
#undef SHOW_TOKENS
#undef ECHO_PROCESSED_COMMANDS

/* Local function prototypes */
int      AddNameSet(NameSet **nameSetList, int numNameSets, char **nameSet, int numNames);
int      AddToSet (int i, int j, int k, int id);
int      AllocCharacters (void);
int      AllocMatrix (void);
int      AllocTaxa (void);
char     ChangeCase (char c);
int      CharacterCode (char ch, int *charCode, int chType);
int      CharacterNumber (int charCode, int chType);
int      CheckInitialPartitions (void);
int      Dex (TreeNode *p);
int      DoAbout (void);
int      DoAcknowledgments (void);
int      DoBeginParm (char *parmName, char *tkn);
int      DoBreaks (void);
int      DoBreaksParm (char *parmName, char *tkn);
int      DoCalibrate (void);
int      DoCalibrateParm (char *parmName, char *tkn);
int      DoCharset (void);
int      DoCharsetParm (char *parmName, char *tkn);
int      DoCharStat (void);
int      DoCitations (void);
int      DoConstraint (void);
int      DoConstraintParm (char *parmName, char *tkn);
int      DoCtype (void);
int      DoCtypeParm (char *parmName, char *tkn);
int      DoDelete (void);
int      DoDeleteParm (char *parmName, char *tkn);
int      DoDimensions (void);
int      DoDimensionsParm (char *parmName, char *tkn);
int      DoDisclaimer (void);
int      DoEndBlock (void);
int      DoExecuteParm (char *parmName, char *tkn);
int      DoExclude (void);
int      DoExcludeParm (char *parmName, char *tkn);
int      DoFormat (void);
int      DoFormatParm (char *parmName, char *tkn);
int      DoHelp (void);
int      DoHelpParm (char *parmName, char *tkn);
int      DoInclude (void);
int      DoIncludeParm (char *parmName, char *tkn);
int      DoLog (void);
int      DoLogParm (char *parmName, char *tkn);
int      DoManual (void);
int      DoManualParm (char *parmName, char *tkn);
int      DoMatrix (void);
int      DoMatrixParm (char *parmName, char *tkn);
int      DoNexusParm (char *parmName, char *tkn);
int      DoOutgroup (void);
int      DoOutgroupParm (char *parmName, char *tkn);
int      DoPairs (void);
int      DoPairsParm (char *parmName, char *tkn);
int      DoPartition (void);
int      DoPartitionParm (char *parmName, char *tkn);
int      DoRestore (void);
int      DoRestoreParm (char *parmName, char *tkn);
int      DoSet (void);
int      DoSetParm (char *parmName, char *tkn);
int      DoShowBeagle (void);
int      DoShowMatrix (void);
int      DoShowUserTrees (void);
int      DoSpeciespartition (void);
int      DoSpeciespartitionParm (char *parmName, char *tkn);
int      DoTaxaset (void);
int      DoTaxasetParm (char *parmName, char *tkn);
int      DoTaxaStat (void);
int      DoTaxlabels (void);
int      DoTaxlabelsParm (char *parmName, char *tkn);
int      DoTranslate (void);
int      DoTranslateParm (char *parmName, char *tkn);
int      DoTree (void);
int      DoTreeParm (char *parmName, char *tkn);
int      DoUserTree (void);
int      DoUserTreeParm (char *parmName, char *tkn);
int      DoVersion (void);
int      FindValidParam (char *tk, int *numMatches);
int      FreeCharacters (void);
int      FreeMatrix (void);
int      FreeTaxa (void);
int      GetNumPartDivisions (int n);
int      GetUserHelp (char *helpTkn);
int      IsAmbig (int charCode, int dType);
int      IsMissing (int charCode, int dType);
int      MBResID (char nuc);
int      NucID (char nuc);
void     PrintSettings (char *command);
void     PrintYesNo (int yn, char s[4]);
int      ProtID (char aa);
int      SetPartition (int part);
int      SetSpeciespartition (int part);
int      SetTaxaFromTranslateTable (void);
int      StandID (char nuc);
void     WhatVariableExp (BitsLong exp, char *st);
MrBFlt   WhichCont (int x);

/* globals */
int             autoClose;             /* autoclose                                     */
int             autoOverwrite;         /* Overwrite or append outputfiles when nowarnings=yes */
Calibration     *calibrationPtr;       /* ptr to calibration being set                  */
CharInformation *charInfo;             /* holds critical information about characters   */
BitsLong        **charSet;             /* holds information about defined charsets      */
char            **charSetNames;        /* holds names of character sets                 */
Comptree        comptreeParams;        /* holds parameters for comparetree command      */
char            **constraintNames;     /* holds names of constraints                    */
int             dataType;              /* type of data                                  */
Calibration     defaultCalibration;    /* default calibration                           */
BitsLong        **definedConstraint;          /* bitfields representing taxa sets of defined constraints                                             */
BitsLong        **definedConstraintTwo;       /* bitfields representing second taxa sets of defined constraints (used for PARTIAL constraints)       */
BitsLong        **definedConstraintPruned;    /* bitfields representing taxa sets of defined constraints after delited taxa are removed              */
BitsLong        **definedConstraintTwoPruned; /* bitfields representing second taxa sets of defined constraints for PARTIAL constraints after delited*/
                                              /* taxa are removed and for NEGATIVE constraint it contains complements of definedConstraintPruned     */
int             echoMB;                /* flag used by Manual to prevent echoing        */
BitsLong        expecting;             /* variable denoting expected token type         */
int             foundNewLine;          /* whether a new line has been found             */
int             inComment;             /* flag for whether input stream is commented    */
int             inComparetreeCommand;  /* flag set whenever you enter comparetree cmd   */
int             inferAncStates;        /* should ancestral states be inferred (y/n)     */
int             inferSiteOmegas;       /* should site omegas be inferred (y/n)          */
int             inferSiteRates;        /* should site rates be inferred (y/n)           */
int             inMrbayesBlock;        /* flag for whether we are in a mrbayes block    */
int             inSumtCommand;         /* flag set whenever you enter sumt cmd          */
int             inTreesBlock;          /* flag for whether we are in a trees block      */
int             inValidCommand;        /* a useful flag set whenever you enter a cmd    */
int             isInAmbig, isInPoly;   /* flags whether we are within () or {}          */
int             isTaxsetDef;           /* is a taxon set defined                        */
int             isTranslateDef;        /* is a translation block defined                */
int             isTranslateDiff;       /* is translate different from current taxaset?  */
char            logFileName[100];      /* name of the log file                          */
int             logToFile;             /* should screen output be logged to a file      */
FILE            *logFileFp;            /* file pointer to log file                      */
int             longIntegerSize;       /* size of an unsigned integer                   */
char            manFileName[100];      /* name of the file for the command help info    */
int             *matrix;               /* matrix containing original data               */
int             matrixHasPoly;         /* flag for whether matrix has polymorphisms     */
int             memAllocs[NUM_ALLOCS]; /* allocated memory flags                        */
int             mode;                  /* mode of program (interactive/noninteractive)  */
Calibration     *nodeCalibration;      /* holds information about node calibrations     */
int             noWarn;                /* no warnings on overwriting files              */
int             numChar;               /* number of characters in character matrix      */
int             numCharSets;           /* number of character sets                      */
int             numComments;           /* counts how deeply nested a comment is         */
int             numDefinedConstraints; /* number of constraints defined                 */
int             numDefinedPartitions;  /* number of partitions defined                  */
int             numDefinedSpeciespartitions;  /* number of speciespartitions defined    */
int             numNamedTaxa;          /* number of named taxa during parsing of cmd    */
int             numOpenExeFiles;       /* number of execute files open                  */
int             numSpecies;            /* number of species in current speciespartition */
int             numTaxa;               /* number of taxa in character matrix            */
int             numTaxaSets;           /* number of taxa sets                           */
int             numTranslates;         /* number of taxa in active translate block      */
int             outGroupNum;           /* number of outgroup taxon                      */
ParmInfo        paramTable[NUMPARAMS]; /* information on parameters                     */
char            **partitionNames;      /* hold names of partitions (first is "default") */
int             **partitionId;         /* holds information about defined partitions    */
int             partitionNum;          /* index of current partition                    */
Plot            plotParams;            /* holds parameters for plot command             */
int             precision;             /* precision of samples and summary stats        */
int             quitOnError;           /* quit on error?                                */
int             replaceLogFile;        /* should logfile be replace/appended to         */
int             scientific;            /* use scientific format for samples ?           */
char            spacer[10];            /* holds blanks for printing indentations        */
NameSet         *speciesNameSets;      /* hold species name sets, one for each speciespartition     */
int             **speciespartitionId;  /* holds info about defined speciespartitions    */
char            **speciespartitionNames;    /* hold names of speciespartitions (first is "default") */
int             speciespartitionNum;   /* index of current speciespartition             */
Sump            sumpParams;            /* holds parameters for sump command             */
Sumt            sumtParams;            /* holds parameters for sumt command             */
Sumss           sumssParams;           /* holds parameters for sumss command            */
TaxaInformation *taxaInfo;             /* holds critical information about taxa         */
char            **taxaNames;           /* holds name of taxa                            */
BitsLong        **taxaSet;             /* holds information about defined taxasets      */
char            **taxaSetNames;        /* holds names of taxa sets                      */
int             *tempActiveConstraints;/* temporarily holds active constraints size allcated        */
enum ConstraintType  *definedConstraintsType;  /* Store type of constraint              */
int             *tempSet;              /* temporarily holds defined set                 */
int             *tempSetNeg;           /* holds bitset of negative set of taxa for partial constraint*/
int             theAmbigChar;          /* int containing ambiguous character            */
Calibration     *tipCalibration;       /* holds information about node calibrations     */
char            **transFrom;           /* translation block information                 */
char            **transTo;             /* translation block information                 */
int             userBrlensDef;         /* are the branch lengths on user tree defined   */

#if defined (BEAGLE_ENABLED)
int             tryToUseBEAGLE;        /* try to use the BEAGLE library                 */
int             beagleScalingScheme;   /* BEAGLE dynamic scaling                        */
int             beagleScalingFrequency;/* BEAGLE dynamic scaling frequency              */
long            beagleFlags;           /* BEAGLE required resource flags                */
int             beagleResourceNumber;  /* BEAGLE resource number                        */
int             *beagleResource;       /* BEAGLE resource choice list                   */
int             beagleResourceCount;   /* BEAGLE resource choice list length            */
int             beagleInstanceCount;   /* total number of BEAGLE instances              */
#endif

#if defined (THREADS_ENABLED)
int             tryToUseThreads;       /* try to use pthreads with BEAGLE library       */
#endif

/* local (to this file) */
char            *tokenP, token[CMD_STRING_LENGTH], *cmdStr=NULL;
Calibration     defaultCalibration = {
                    "Unconstrained",      /* name */
                    unconstrained,        /* prior */
                    { -1.0, -1.0, -1.0 }, /* priorParams */
                    NULL,                 /* LnPriorProb */
                    NULL,                 /* LnPriorRatio */
                    -1.0,                 /* min */
                    -1.0                  /* max */
                };

CmdType     commands[] =
            {
            /*  Information on commands initialization:
             
                    1 = Command number (cmdNumber)
                    2 = Command name (string)
                    3 = Special command (YES/NO) (specialCmd) 
                    4 = Pointer to finishing function (fp)
                    5 = Number of valid parameters (numParms)
                    6 = List of valid parameters (parmList) 
                    7 = Expecting (2^TokenType) (expect) (PARAMETER = 4; SEMICOLON = 32; ALPHA = 16384; 
                        ALPHA | QUESTIONMARK | DASH | NUMBER | ASTERISK | EXCLAMATIONMARK | PERCENT | WEIRD | SEMICOLON = 11715360;
                        ALPHA | QUESTIONMARK | DASH | NUMBER | ASTERISK | EXCLAMATIONMARK | PERCENT | WEIRD | VERTICALBAR | SEMICOLON | LEFTPAR | RIGHTPAR | LEFTCURL | RIGHTCURL = 649252640;
                        PARAMETER | SEMICOLON = 36; NUMBER | ALPHA = 49152; ALPHA | SEMICOLON = 16416; EQUALSIGN = 8; NUMBER = 32768)
                    8 = Description of the command (cmdDescription)
                    9 = Where should the command be used (cmdUse) (IN_CMD = used from command line or mrbayes block; IN_FILE = used in data block or in tree block)
                   10 = Should the command be shown when "help" is typed (hiding).
             
              #1                 #2   #3                 #4  #5                                                                                                #6        #7                                                            #8       #9   #10
             -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- */
            {  0,               "#",  NO,              NULL,  1,                                                                                              {0},        4,                                                           "", IN_FILE, HIDE },
            {  1,           "About",  NO,           DoAbout,  0,                                                                                             {-1},       32,                                      "Describes the program",  IN_CMD, SHOW },
            {  2, "Acknowledgments",  NO, DoAcknowledgments,  0,                                                                                             {-1},       32,                              "Shows program acknowledgments",  IN_CMD, SHOW },
            {  3,           "Begin",  NO,              NULL,  6,                                                                              {1,2,3,201,226,227},        4,                         "Denotes beginning of block in file", IN_FILE, SHOW },
            {  4,       "Calibrate",  NO,       DoCalibrate,  1,                                                                                            {119},        4,               "Assigns dates to terminals or interior nodes",  IN_CMD, SHOW },
            {  5,         "Charset",  NO,         DoCharset,  1,                                                                                             {15},        4,                          "Assigns a group of sites to a set",  IN_CMD, SHOW },
            {  6,        "Charstat",  NO,        DoCharStat,  0,                                                                                             {-1},       32,                                 "Shows status of characters",  IN_CMD, SHOW },
            {  7,       "Citations",  NO,       DoCitations,  0,                                                                                             {-1},       32,                   "Citation of program, models, and methods",  IN_CMD, SHOW },
            {  8,     "Comparetree",  NO,     DoCompareTree,  7,                                                                    {127,128,129,130,221,222,223},       36,                     "Compares the trees from two tree files",  IN_CMD, SHOW },
            {  9,      "Constraint",  NO,      DoConstraint,  1,                                                                                             {66},        4,                      "Defines a constraint on tree topology",  IN_CMD, SHOW },
            { 10,           "Ctype",  NO,           DoCtype,  1,                                                                                             {65},        4,                        "Assigns ordering for the characters",  IN_CMD, SHOW },
            { 11,      "Databreaks", YES,          DoBreaks,  1,                                                                                             {93},    32768,           "Defines data breaks for autodiscrete gamma model",  IN_CMD, SHOW },
            { 12,          "Delete", YES,          DoDelete,  1,                                                                                             {47},    49152,                             "Deletes taxa from the analysis",  IN_CMD, SHOW },
            { 13,      "Dimensions",  NO,      DoDimensions,  2,                                                                                            {4,5},        4,                           "Defines size of character matrix", IN_FILE, SHOW },
            { 14,      "Disclaimer",  NO,      DoDisclaimer,  0,                                                                                             {-1},       32,                               "Describes program disclaimer",  IN_CMD, SHOW },
            { 15,             "End",  NO,        DoEndBlock,  0,                                                                                             {-1},       32,                             "Denotes end of a block in file", IN_FILE, SHOW },
            { 16,        "Endblock",  NO,        DoEndBlock,  0,                                                                                             {-1},       32,                 "Alternative way of denoting end of a block", IN_FILE, SHOW },
            { 17,         "Exclude", YES,         DoExclude,  1,                                                                                             {45},    49152,                           "Excludes sites from the analysis",  IN_CMD, SHOW },
            { 18,         "Execute", YES,         DoExecute,  1,                                                                                             {12},    16384,                                            "Executes a file",  IN_CMD, SHOW },
            { 19,          "Format",  NO,          DoFormat,  7,                                                                             {6,7,8,9,10,219,220},        4,                     "Defines character format in data block", IN_FILE, SHOW },
            { 20,            "Help", YES,            DoHelp,  1,                                                                                             {50},    16416,                  "Provides detailed description of commands",  IN_CMD, SHOW },
            { 21,         "Include", YES,         DoInclude,  1,                                                                                             {46},    49152,                                             "Includes sites",  IN_CMD, SHOW },
            { 22,            "Link",  NO,            DoLink, 30,  {55,56,57,58,59,60,61,62,63,72,73,74,75,76,105,118,193,194,195,196,197,242,243,252,253,255,256,
                                                                                                                                                     270,273,274},        4,               "Links parameters across character partitions",  IN_CMD, SHOW },
            { 23,             "Log",  NO,             DoLog,  5,                                                                                 {85,86,87,88,89},        4,                               "Logs screen output to a file",  IN_CMD, SHOW },
            { 24,            "Lset",  NO,            DoLset, 18,                                     {28,29,30,31,32,33,34,40,51,52,53,90,91,131,188,189,276,277},        4,                "Sets the parameters of the likelihood model",  IN_CMD, SHOW },
            { 25,          "Manual",  NO,          DoManual,  1,                                                                                            {126},       36,                  "Prints a command reference to a text file",  IN_CMD, SHOW },
            { 26,          "Matrix", YES,          DoMatrix,  1,                                                                                             {11},649252640,                 "Defines matrix of characters in data block", IN_FILE, SHOW },
            { 27,            "Mcmc",  NO,            DoMcmc, 46,  {17,18,19,20,21,22,23,24,25,26,27,84,98,112,113,114,115,116,132,142,143,144,148,149,150,151,152,
                                                                                     153,154,155,156,157,158,159,160,166,169,190,191,198,199,200,202,213,214,215},       36,                   "Starts Markov chain Monte Carlo analysis",  IN_CMD, SHOW },
            { 28,           "Mcmcp",  NO,           DoMcmcp, 46,  {17,18,19,20,21,22,23,24,25,26,27,84,98,112,113,114,115,116,132,142,143,144,148,149,150,151,152,
                                                                                     153,154,155,156,157,158,159,160,166,169,190,191,198,199,200,202,213,214,215},        4,     "Sets parameters of a chain (without starting analysis)",  IN_CMD, SHOW },
            { 29,        "Outgroup", YES,        DoOutgroup,  1,                                                                                             {78},    49152,                                     "Changes outgroup taxon",  IN_CMD, SHOW },
            { 30,           "Pairs", YES,           DoPairs,  1,                                                                                             {92},    32768,        "Defines nucleotide pairs (doublets) for stem models",  IN_CMD, SHOW },
            { 31,       "Partition",  NO,       DoPartition,  1,                                                                                             {16},        4,                              "Assigns a character partition",  IN_CMD, SHOW },
            { 32,            "Plot",  NO,            DoPlot,  6,                                                                        {106,107,108,109,224,225},       36,                        "Plots parameters from MCMC analysis",  IN_CMD, SHOW },
            { 33,           "Prset",  NO,           DoPrset, 43,  {35,36,37,38,39,41,42,43,44,54,64,67,68,69,70,71,77,100,101,102,103,104,110,111,117,120,121,133,
                                                                                                 168,172,173,174,183,184,185,218,241,246,247,251,254,269,271,272},        4,                         "Sets the priors for the parameters",  IN_CMD, SHOW },
            { 34,         "Propset",  NO,         DoPropset,  1,                                                                                            {186},        4,          "Sets proposal probabilities and tuning parameters",  IN_CMD, SHOW },
            { 35,            "Quit",  NO,            DoQuit,  0,                                                                                             {-1},       32,                                          "Quits the program",  IN_CMD, SHOW },
            { 36,          "Report",  NO,          DoReport,  9,                                                            {122,123,124,125,134,135,136,192,217},        4,                 "Controls how model parameters are reported",  IN_CMD, SHOW },
            { 37,         "Restore", YES,         DoRestore,  1,                                                                                             {48},    49152,                                              "Restores taxa",  IN_CMD, SHOW },
            { 38,             "Set",  NO,             DoSet, 22,           {13,14,94,145,170,171,179,181,182,216,229,233,234,235,236,237,238,239,240,245,268,275},        4,      "Sets run conditions and defines active data partition",  IN_CMD, SHOW },
            { 39,      "Showbeagle",  NO,      DoShowBeagle,  0,                                                                                             {-1},       32,                            "Show available BEAGLE resources",  IN_CMD, SHOW },
            { 40,      "Showmatrix",  NO,      DoShowMatrix,  0,                                                                                             {-1},       32,                             "Shows current character matrix",  IN_CMD, SHOW },
            { 41,   "Showmcmctrees",  NO,   DoShowMcmcTrees,  0,                                                                                             {-1},       32,                          "Shows trees used in mcmc analysis",  IN_CMD, SHOW },
            { 42,       "Showmodel",  NO,       DoShowModel,  0,                                                                                             {-1},       32,                                       "Shows model settings",  IN_CMD, SHOW },
            { 43,       "Showmoves",  NO,       DoShowMoves,  1,                                                                                            {180},       36,                              "Shows moves for current model",  IN_CMD, SHOW },
            { 44,      "Showparams",  NO,      DoShowParams,  0,                                                                                             {-1},       32,                          "Shows parameters in current model",  IN_CMD, SHOW },
            { 45,   "Showusertrees",  NO,   DoShowUserTrees,  0,                                                                                             {-1},       32,                                   "Shows user-defined trees",  IN_CMD, SHOW },
            { 46,"Speciespartition",  NO,DoSpeciespartition,  1,                                                                                            {244},        4,                   "Defines a partition of tips into species",  IN_CMD, SHOW },
            { 47,              "Ss",  NO,              DoSs, 50,  {17,18,19,20,21,22,23,24,25,26,27,84,98,112,113,114,115,116,132,142,143,144,148,149,150,151,152,
                                                                     153,154,155,156,157,158,159,160,166,169,190,191,198,199,200,202,213,214,215,248,249,250,257},       36,                             "Starts stepping-stone sampling",  IN_CMD, SHOW },
            { 48,             "Ssp",  NO,             DoSsp, 50,  {17,18,19,20,21,22,23,24,25,26,27,84,98,112,113,114,115,116,132,142,143,144,148,149,150,151,152,
                                                                     153,154,155,156,157,158,159,160,166,169,190,191,198,199,200,202,213,214,215,248,249,250,257},       36,"Sets parameters of stepping-stone analysis (without starting)",IN_CMD, SHOW },
            { 49,       "Startvals",  NO,       DoStartvals,  1,                                                                                            {187},        4,                         "Sets starting values of parameters",  IN_CMD, SHOW },
            { 50,            "Sump",  NO,            DoSump, 13,                                              {96,97,137,138,139,140,141,161,162,178,211,212,231},       36,                   "Summarizes parameters from MCMC analysis",  IN_CMD, SHOW },
            { 51,           "Sumss",  NO,           DoSumSs, 10,                                                        {258,259,260,261,262,263,264,265,266,267},       36,         "Summarizes parameters from stepping-stone analysis",  IN_CMD, SHOW },
            { 52,            "Sumt",  NO,            DoSumt, 21,                {80,81,82,95,146,147,163,164,165,167,175,177,204,205,206,207,208,209,210,230,232},       36,                        "Summarizes trees from MCMC analysis",  IN_CMD, SHOW },
            { 53,        "Taxastat",  NO,        DoTaxaStat,  0,                                                                                             {-1},       32,                                       "Shows status of taxa",  IN_CMD, SHOW },
            { 54,          "Taxset",  NO,         DoTaxaset,  1,                                                                                             {49},        4,                           "Assigns a group of taxa to a set",  IN_CMD, SHOW },
            { 55,       "Taxlabels", YES,       DoTaxlabels,  1,                                                                                            {228},    49152,                                       "Defines taxon labels", IN_FILE, SHOW },
            { 56,       "Translate", YES,       DoTranslate,  1,                                                                                             {83},    49152,                         "Defines alternative names for taxa", IN_FILE, SHOW },
            { 57,            "Tree",  NO,            DoTree,  1,                                                                                             {79},        4,                                             "Defines a tree", IN_FILE, SHOW },
            { 58,          "Unlink",  NO,          DoUnlink, 30,  {55,56,57,58,59,60,61,62,63,72,73,74,75,76,105,118,193,194,195,196,197,242,243,252,253,255,256,
                                                                                                                                                     270,273,274},        4,             "Unlinks parameters across character partitions",  IN_CMD, SHOW },
            { 59,        "Usertree", YES,        DoUserTree,  1,                                                                                            {203},        8,                                 "Defines a single user tree",  IN_CMD, HIDE },
            { 60,         "Version",  NO,         DoVersion,  0,                                                                                             {-1},       32,                                      "Shows program version",  IN_CMD, SHOW },
            { 61,      "Compareref",  NO,     DoCompRefTree,  7,                                                                    {127,128,129,130,221,222,223},       36,                     "Compares the tree to the reference trees",  IN_CMD, HIDE },
            /* NOTE: If you add a command here, make certain to change NUMCOMMANDS (above, in this file) appropriately! */
            { 999,             NULL,  NO,              NULL,  0,                                                                                             {-1},       32,                                                           "",  IN_CMD, HIDE }  
            };
int                 inDataBlock, inForeignBlock, isInterleaved, isFirstMatrixRead, isFirstInterleavedBlock, 
                    taxonCount, fromI, toJ, everyK, foundDash, foundSlash, foundFirst, isMixed, whichPartition,
                    isNegative, numDivisions, charOrdering, foundExp, foundColon, isFirstNode, nextAvailableNode,
                    pairId, firstPair, inTaxaBlock, inCharactersBlock, foundEqual;
char                gapId, missingId, matchId, tempSetName[100], **tempNames;
CmdType             *commandPtr; /* Points to the commands array entry which corresponds to currently processed command */
ParmInfoPtr         paramPtr;    /* Points to paramTable table array entry which corresponds to currently processed parameter of current command */
TreeNode            *pPtr, *qPtr;

enum ConstraintType constraintType; /* Used only in processing of constraint command to indicate the type of constraint */


int AddToGivenSet (int i, int j, int k, int id, int *Set)
{
    int     m, n;
    
    if (id <= 0)
        {
        MrBayesPrint ("%s   The id for a temporary set should be greater than 0\n", spacer);
        return (ERROR);
        }
    
    if (i < 0 && j < 0)
        return (ERROR);
    else if (i < 0 && j >= 0)
        return (ERROR);
    else if (i >= 0 && j < 0)
        {
        if (k >= 0)
            return (ERROR);
        else
            {
            if (Set[i] != 0)
                {
                MrBayesPrint ("%s   Character %d defined more than once\n", spacer, i+1);
                return (ERROR);
                }
            Set[i] = id;
            }
        }
    else if (i >= 0 && j >= 0)
        {
        if (k < 0)
            {
            for (m=i; m<=j; m++)
                {
                if (Set[m] != 0)
                    {
                    MrBayesPrint ("%s   Character %d defined more than once\n", spacer, m+1);
                    return (ERROR);
                    }
                Set[m] = id;
                }
            }
        else
            {
            n = k;
            for (m=i; m<=j; m++)    
                {
                if (n % k == 0)
                    {
                    if (Set[m] != 0)
                        {
                        MrBayesPrint ("%s   Character %d defined more than once\n", spacer, m+1);
                        return (ERROR);
                        }
                    Set[m] = id;
                    }
                n++;
                }
            }
        }

    return (NO_ERROR);
    
}


int AddToSet (int i, int j, int k, int id)
{
    return AddToGivenSet (i, j, k,id, tempSet);
}


/* AddNameSet: Push a name set onto the end of a list of name sets, with reallocation
      of list to hold the extra element. The calling function needs to keep track of
      the counter holding the length of the list. */
int AddNameSet (NameSet **nameSetList, int numNameSets, char **nameSet, int numNames)
{
    int     i;

    (*nameSetList) = (NameSet*) SafeRealloc ((void*)(*nameSetList), ((size_t)numNameSets+1)*sizeof(NameSet));

    (*nameSetList)[numNameSets].names    = NULL;
    (*nameSetList)[numNameSets].numNames = numNames;

    for (i=0; i<numNames; i++)
        AddString(&((*nameSetList)[numNameSets].names), i, nameSet[i]);
    
    return NO_ERROR;
}


/* AddString: Push a string onto the end of a list, with reallocation of list
      to hold the extra element. The calling function needs to keep track of
      the counter holding the length of the list. */
int AddString (char ***list, int len, char *token)
{
    (*list) = (char **) SafeRealloc ((void *)(*list), ((size_t)len+1)*sizeof(char*));
    if (!(*list))
        return ERROR;

    (*list)[len] = (char *) SafeCalloc ((strlen(token)+1), sizeof(char));
    if (!(*list)[len])
        return ERROR;

    strcpy ((*list)[len], token);
    
    return NO_ERROR;
}


int AllocCharacters (void)
{
    int     i, tempSetSize;

    if (memAllocs[ALLOC_MATRIX] == YES)
        goto errorExit;
    matrix = (int *) SafeMalloc((size_t)numTaxa * (size_t)numChar * sizeof(int));
    if (!matrix)
        {
        MrBayesPrint ("%s   Problem allocating matrix (%d)\n", spacer, numTaxa * numChar * sizeof(int));
        goto errorExit;
        }
    for (i=0; i<numTaxa * numChar; i++)
        matrix[i] = 0;
    memAllocs[ALLOC_MATRIX] = YES;

    if (memAllocs[ALLOC_CHARINFO] == YES)
        goto errorExit;
    charInfo = (CharInformation *) SafeMalloc ((size_t)numChar * sizeof(CharInformation));
    if (!charInfo)
        {
        MrBayesPrint ("%s   Problem allocating charInfo (%d)\n", spacer, numChar * sizeof(CharInformation));
        goto errorExit;
        }
    for (i=0; i<numChar; i++)
        {
        charInfo[i].isExcluded = NO;
        charInfo[i].numStates = 0;
        charInfo[i].charType = 0;
        charInfo[i].isMissAmbig = NO;
        charInfo[i].ctype = UNORD;
        charInfo[i].charId = 0;
        charInfo[i].pairsId = 0;
        charInfo[i].bigBreakAfter = NO;
        }
    memAllocs[ALLOC_CHARINFO] = YES;

    if (memAllocs[ALLOC_CHARSETS] == YES)
        goto errorExit;
    charSetNames = NULL;
    charSet = NULL;
    numCharSets = 0;
    memAllocs[ALLOC_CHARSETS] = YES;    /* safe to do free */

    if (memAllocs[ALLOC_PARTITIONS] == YES)
        goto errorExit;
    partitionNames = NULL;
    partitionId = (int**) SafeMalloc ((size_t)numChar * sizeof(int*));
    for (i=0; i<numChar; i++)
        partitionId[i] = (int *) SafeMalloc (sizeof(int));
    numDefinedPartitions = 0;   /* number of defined partitions */
    memAllocs[ALLOC_PARTITIONS] = YES;  /* safe to do free */

    if (memAllocs[ALLOC_PARTITIONVARS] == YES)
        goto errorExit;
    numVars           = NULL;
    tempLinkUnlinkVec = NULL;
    activeParts       = NULL;
    tempLinkUnlinkVec = NULL;
    tempNum           = NULL;
    linkTable[0]      = NULL;
    tempLinkUnlink[0] = NULL;
    for (i=0; i<NUM_LINKED; i++)
        {
        linkTable[i]      = NULL;
        tempLinkUnlink[i] = NULL;
        activeParams[i]   = NULL;
        }
    memAllocs[ALLOC_PARTITIONVARS] = YES;

    if (memAllocs[ALLOC_TMPSET] == NO)
        goto errorExit;
    if (numChar > numTaxa)
        tempSetSize = numChar;
    else
        tempSetSize = numTaxa;
    tempSet = (int *) SafeRealloc ((void *)tempSet, (size_t)tempSetSize * sizeof(int));
    tempSetNeg = (int *) SafeRealloc ((void *)tempSetNeg, (size_t)tempSetSize * sizeof(int));
    if (!tempSet || !tempSetNeg)
        {
        MrBayesPrint ("%s   Problem reallocating tempSet (%d)\n", spacer, tempSetSize * sizeof(int));
        goto errorExit;
        }

    MrBayesPrint ("%s   Allocated matrix\n", spacer);
    return (NO_ERROR);

    errorExit:
        MrBayesPrint ("%s   Problem allocating matrix\n", spacer);
        FreeMatrix();
        return (ERROR);
}


int AllocMatrix (void)
{
    if (memAllocs[ALLOC_TAXA] == NO && AllocTaxa() == ERROR)
        return ERROR;
    else
        return (AllocCharacters());
}


int AllocTaxa (void)
{
    int i;

    if (defTaxa==NO)
        {
        MrBayesPrint ("%s   Number of taxa not defined\n", spacer);
        return (ERROR);
        }
    if (numTaxa == 0)
        {
        MrBayesPrint ("%s   Number of taxa is 0\n", spacer);
        return (ERROR);
        }

    /* allocate space for taxa */
    if (memAllocs[ALLOC_TAXA] == YES)
        goto errorExit;
    taxaNames = NULL;   /* This variable is allocated in AddString */
    taxaInfo = (TaxaInformation *) SafeMalloc ((size_t)numTaxa * sizeof(TaxaInformation));
    if (!taxaInfo)
        {
        goto errorExit;
        }
    tipCalibration = (Calibration *) SafeMalloc ((size_t)numTaxa * sizeof(Calibration));
    if (!tipCalibration)
        {
        free (taxaInfo);
        taxaInfo = NULL;
        goto errorExit;
        }
    for (i=0; i<numTaxa; i++)
        {
        taxaInfo[i].isDeleted = NO;
        taxaInfo[i].charCount = 0;
        }
    memAllocs[ALLOC_TAXA] = YES;

    /* taxa sets */
    if (memAllocs[ALLOC_TAXASETS] == YES)
        goto errorExit;
    taxaSetNames = NULL;
    taxaSet = NULL;
    numTaxaSets = 0;
    memAllocs[ALLOC_TAXASETS] = YES;    /* safe to free */

    /* species partitions; allocate space and set default species partition */
    if (memAllocs[ALLOC_SPECIESPARTITIONS] == YES)
        goto errorExit;
    speciespartitionNames = NULL;
    speciesNameSets = NULL;
    speciespartitionId = (int**) SafeMalloc ((size_t)numTaxa * sizeof(int*));
    for (i=0; i<numTaxa; i++)
        {
        speciespartitionId[i] = (int *) SafeMalloc (sizeof(int));
        speciespartitionId[i][0] = i + 1;   /* 1-based taxon index, do not ask me why */
        }
    numDefinedSpeciespartitions = 0;   /* number of defined species partitions */
    memAllocs[ALLOC_SPECIESPARTITIONS] = YES;  /* safe to do free */

    /* constraints */
    if (memAllocs[ALLOC_CONSTRAINTS] == YES)
        goto errorExit;
    constraintNames = NULL;
    definedConstraintsType = NULL; 
    definedConstraint = NULL;
    definedConstraintTwo = NULL;
    definedConstraintPruned = NULL;
    definedConstraintTwoPruned = NULL;   
    numDefinedConstraints = 0;
    tempActiveConstraints = NULL;
    memAllocs[ALLOC_CONSTRAINTS] = YES;     /* safe to free */

    /* translate table */
    transFrom = NULL;
    transTo = NULL;
    numTranslates = 0;

    /* tempSet */
    if (memAllocs[ALLOC_TMPSET] == YES)
        goto errorExit;
    tempSet = (int *) SafeMalloc ((size_t)numTaxa * sizeof(int));
    tempSetNeg = (int *) SafeMalloc ((size_t)numTaxa * sizeof(int));
    if (!tempSet || !tempSetNeg)
        goto errorExit;
    memAllocs[ALLOC_TMPSET] = YES;

    /* make sure previous user trees are freed */
    if (numUserTrees > 0)
        {
        MrBayesPrint ("%s   Previous user trees not freed\n", spacer);
        goto errorExit;
        }

    MrBayesPrint ("%s   Allocated taxon set\n", spacer);
    return NO_ERROR;

errorExit:
    MrBayesPrint ("%s   Problem allocating taxon set\n", spacer);
    FreeTaxa();
    return ERROR;
}


char ChangeCase (char c)
{
    int     x;
    
    x = tolower(c);
    return (x);
}


int CharacterCode (char ch, int *charCode, int chType)
{
    if (chType == DNA || chType == RNA)
        {
        if ((*charCode = NucID (ch)) == -1)
            {
            MrBayesPrint ("%s   Unrecognized DNA/RNA character '%c'\n", spacer, ch);
            return (ERROR);
            }
        }
    else if (chType == PROTEIN)
        {
        if ((*charCode = ProtID (ch)) == -1)
            {
            MrBayesPrint ("%s   Unrecognized Protein character '%c'\n", spacer, ch);
            return (ERROR);
            }
        }
    else if (chType == RESTRICTION)
        {
        if ((*charCode = MBResID (ch)) == -1)
            {
            MrBayesPrint ("%s   Unrecognized Restriction character '%c'\n", spacer, ch);
            return (ERROR);
            }
        }
    else if (chType == STANDARD)
        {
        if ((*charCode = StandID (ch)) == -1)
            {
            MrBayesPrint ("%s   Unrecognized Standard character '%c'\n", spacer, ch);
            return (ERROR);
            }
        }
    else if (chType == CONTINUOUS)
        {
        MrBayesPrint ("%s   CharacterCode function cannot check continuous characters\n", spacer);
        }
    else
        {
        MrBayesPrint ("%s   Unrecognized character type (%d)\n", spacer, chType);
        return (ERROR);
        }
        
    return (NO_ERROR);
}


int CharacterNumber (int charCode, int chType)
{
    int i, x = charCode;
    
    if (chType == CONTINUOUS)
        return 0;

    for (i=0; x!=0; i++)
        x >>= 1;

    return (i);
}


int CheckInitialPartitions (void)
{
    int     i;
    
    for (i=0; i<numChar; i++)
        {
        if (partitionId[i][0] <= 0 || partitionId[i][0] > numDivisions)
            {
            MrBayesPrint ("%s   The partition for site %d is incorrect\n", spacer, i+1); 
            return (ERROR);
            }
        }
        
    return (NO_ERROR);
}


int CheckStringValidity (char *s)
{
    int         i, numUnknownChars, tempNumComments, tempInComment;
    char        temp[100];

    i = 0;
    numUnknownChars = 0;
    tempNumComments = numComments;
    tempInComment = inComment;

    while (s[i] != '\0')
        {
        if (tempInComment == NO)
            {
            if (!IsIn(s[i],"=abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.;:,#()[]?-*/'\\'!%\"&~+^$@|{}`>< "))
                {
                if (IsWhite(s[i]) == 1 || IsWhite(s[i]) == 2)
                    {
                    
                    }
                else
                    {
                    if (commandPtr == NULL) 
                        return (ERROR);
                    MrBayesPrint ("%s   Unknown character \"%c\" (ASCII code %d)\n", spacer, s[i], s[i]);
                    if (!strcmp(commandPtr->string,"Matrix"))
                        {
                        if (foundNewLine == NO)
                            {
                            MrBayesPrint ("%s   The error is in character %d for taxon %s\n", spacer, taxaInfo[taxonCount-1].charCount+i+1, "???"); /* bug? */
                            }
                        else
                            {
                            if (taxonCount == 0)
                                MrBayesPrint ("%s   The error is in the first taxon name\n", spacer);
                            else
                                {
                                strcpy(temp, taxaNames[taxonCount]);
                                if (isInterleaved == NO)
                                    MrBayesPrint ("%s   The error is in the name of the taxon following taxon %s\n", spacer, temp);
                                else
                                    {
                                    MrBayesPrint ("%s   The error is in the name of the taxon following taxon %s\n", spacer, temp);
                                    MrBayesPrint ("%s   in one of the interleaved data blocks\n", spacer);
                                    }
                                }
                            }
                        }
                    else if (!strcmp(commandPtr->string,"Execute"))
                        {
                        MrBayesPrint ("%s   Assuming irrelevant characters at beginning of file; processing continues\n", spacer);
                        return (NO_ERROR);
                        }
                    return (ERROR);
                    }
                }
            if (s[i]=='[')
                {
                tempInComment = YES;
                tempNumComments++;
                }
            }
        else if (tempInComment == YES)
            {
            if (s[i]==']')
                {
                tempNumComments--;
                if (tempNumComments == 0)
                    tempInComment = NO;
                }
            }
        i++;
        }
        
    if (numUnknownChars > 0)
        return (ERROR);
    else
        return (NO_ERROR);
}


/* CheckString: This function simply checks a vector of strings for a match against token.
          Upon return, matchIndex contains the index of the matched string. An
          ERROR is returned if there are no matches.  */
int CheckString (char **list, int len, char *token, int *matchIndex)
{
    int         i;      
        
    *matchIndex = -1;
    for (i=0; i<len; i++)
        {
        if (StrCmpCaseInsensitive(token,list[i]) == 0)
            {
            *matchIndex = i;
            return (NO_ERROR);
            }
        }

    return (ERROR); 
}


int Dex (TreeNode *p)
{
    return (p == NULL) ? -1 : p->index;
}


int DoAbout (void)
{
    MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
    MrBayesPrint ("   About the program                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes is a program for the Bayesian estimation of phylogeny. Bayesian       \n");
    MrBayesPrint ("   inference of phylogeny is based upon the posterior probability distribution   \n");
    MrBayesPrint ("   of trees. Trees are labelled T1, T2, ..., Tn, where n is the number of        \n");
    MrBayesPrint ("   possible trees. The posterior probability of the i-th tree is calculated      \n");
    MrBayesPrint ("   using Bayes\'s formula as                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Pr[Ti | X] = Pr[X | Ti] X Pr[Ti] / Pr[X]                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   where X is a character matrix. Here, \"Pr[Ti | X]\" is the posterior          \n");
    MrBayesPrint ("   probability of the i-th tree, \"Pr[X | Ti]\" is the likelihood of the         \n");
    MrBayesPrint ("   i-th tree, and \"Pr[Ti]\" is the prior probability of the i-th tree. The      \n");
    MrBayesPrint ("   denominator of Bayes\'s formula (\"Pr[X]\") is a normalizing constant that    \n");
    MrBayesPrint ("   involves a summation over all possible trees. The likelihood, as described    \n");
    MrBayesPrint ("   above, cannot be calculated with knowledge of only the tree\'s topology. You  \n");
    MrBayesPrint ("   also need to have information on the lenths of the branches and on the        \n");
    MrBayesPrint ("   mechanism of character change. Hence, the likelihood (\"Pr[X | Ti]\")         \n");
    MrBayesPrint ("   involves a multidimensional integral over all possible combinations of        \n");
    MrBayesPrint ("   branch lengths and substitution model parameters.                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   In practice, it is impossible to calculate the posterior probability dist-    \n");
    MrBayesPrint ("   ribution of trees analytically. Instead, the posterior probability            \n");
    MrBayesPrint ("   of trees must be approximated. MrBayes uses a method called Markov chain      \n");
    MrBayesPrint ("   Monte Carlo (MCMC) to approximate the posterior probability of trees.         \n");
    MrBayesPrint ("   The object of MCMC is to construct a Markov chain that has as its state       \n");
    MrBayesPrint ("   space the parameters of the phylogenetic model and a stationary distribution  \n");
    MrBayesPrint ("   that is the posterior probability distribution of trees. MCMC takes valid,    \n");
    MrBayesPrint ("   albeit dependent, samples from the posterior probability distribution of      \n");
    MrBayesPrint ("   trees. The fraction of the time any tree appears in this sample is a          \n");
    MrBayesPrint ("   valid approximation of the posterior probability of the tree. MrBayes keeps   \n");
    MrBayesPrint ("   track of all the parameters of the phylogenetic model. The trees (with branch \n");
    MrBayesPrint ("   lengths) that were sampled by the MCMC procedure are saved in one file        \n");
    MrBayesPrint ("   (a file with a \".t\" extension) whereas the parameters of the model of       \n");
    MrBayesPrint ("   character change are saved in another file (a file with a \".p\" ext-         \n");
    MrBayesPrint ("   ension). You can summarize the results in the \".t\" and \".p\" files         \n");
    MrBayesPrint ("   using the \"sumt\" and \"sump\" commands, respectively.                       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes was originally written by John Huelsenbeck in August of 2000 and was  \n");
    MrBayesPrint ("   intended to be distributed to a small number of people. In March of 2001,     \n");
    MrBayesPrint ("   Fredrik Ronquist started making contributions to the program. The contribu-   \n");
    MrBayesPrint ("   tions were of such a significant nature that he was made a coauthor of the    \n");
    MrBayesPrint ("   program. Version 3 of MrBayes was a fully joint effort, started in the summer \n");
    MrBayesPrint ("   of 2002 when JPH visited Sweden on a grant from the Wenner-Gren Foundations.  \n");
    MrBayesPrint ("   Several others have contributed to the MrBayes code since then, most notably  \n");
    MrBayesPrint ("   Paul van der Mark, Maxim Teslenko and Chi Zhang, all postdocs/programmers in  \n");
    MrBayesPrint ("   Fredrik's lab. A large number of users and students, too many to list here,   \n");
    MrBayesPrint ("   have also contributed importantly to the project (type 'Acknowledgments' for  \n");
    MrBayesPrint ("   a list of some of them).                                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   Since 2003, MrBayes has been distributed from SourceForge. Bugs can be repor- \n");
    MrBayesPrint ("   ted to the MrBayes site on SourceForge.                                       \n");
    MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

    return (NO_ERROR);
}


int DoAcknowledgments (void)
{
    MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
    MrBayesPrint ("   Acknowledgments                                                               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   JPH and FR would like to thank Gautam Altekar, Andrea Betancourt, Jon         \n");
    MrBayesPrint ("   Bollback, Barry Hall, Jimmy McGuire, Rasmus Nielsen, David Swofford,          \n");
    MrBayesPrint ("   Johan Nylander, Mikael Thollesson, and Derrick Zwickl for help during the     \n");
    MrBayesPrint ("   initial development of this program. Gautam Altekar, especially, was instru-  \n");
    MrBayesPrint ("   mental in getting the parallel version of the program working. Important bug- \n");
    MrBayesPrint ("   fixes and additional functionality was contributed by Clemens Lakner, Sebas-  \n");
    MrBayesPrint ("   tian Hoehna, Paul Lewis, Mark Holder, Julian Catchen and Bret Larget. Marc    \n");
    MrBayesPrint ("   Suchard, Daniel Ayres and Aaron Darling got mrbayes working with beagle and   \n");
    MrBayesPrint ("   contributed a lot of related functionality and bug fixes. Aaron Darling was   \n");
    MrBayesPrint ("   instrumental in getting the Windows installer set up. Liu Liang and Dennis    \n");
    MrBayesPrint ("   Pearl helped integrate MrBayes with BEST.                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   Bug fixes and user support was provided by Paul van der Mark (2005-2007),     \n");
    MrBayesPrint ("   Maxim Teslenko (2010-2012) and Chi Zhang (2012-2015).                         \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   Our wives -- Edna Huelsenbeck and Eva Ronquist -- showed extraordinary        \n");
    MrBayesPrint ("   patience with us while we spent many late nights programming.                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   JPH was supported by NSF grants DEB-007540 and MCB-0075404 and a Wenner-      \n");
    MrBayesPrint ("   Gren scholarship while writing this program. FR was supported by grants       \n");
    MrBayesPrint ("   from the Swedish Natural Science Research Council and the Swedish Research    \n");
    MrBayesPrint ("   Council.                                                                      \n");
    MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

    return (NO_ERROR);
}


int DoBeginParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        /* set Data (inDataBlock) *************************************************************/
        if (!strcmp(parmName, "Data"))
            {
            if (FreeModel () == ERROR)
                return (ERROR);
            if (FreeMatrix () == ERROR)
                return (ERROR);
            MrBayesPrint ("   Reading data block\n");
            inDataBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Characters (inCharactersBlock) *************************************************************/
        else if (!strcmp(parmName, "Characters"))
            {
            if (FreeModel () == ERROR)
                return (ERROR);
            if (FreeCharacters () == ERROR)
                return (ERROR);
            MrBayesPrint ("   Reading characters block\n");
            inCharactersBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Taxa (inTaxaBlock) *************************************************************/
        else if (!strcmp(parmName, "Taxa"))
            {
            if (FreeModel () == ERROR)
                return (ERROR);
            if (FreeMatrix () == ERROR)
                return (ERROR);
            MrBayesPrint ("   Reading taxa block\n");
            inTaxaBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Mrbayes (inMrbayesBlock) *******************************************************/
        else if (!strcmp(parmName, "Mrbayes"))
            {
            MrBayesPrint ("   Reading mrbayes block\n");
            inMrbayesBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Trees (inTreesBlock) *******************************************************/
        else if (!strcmp(parmName, "Trees"))
            {
            MrBayesPrint ("   Reading trees block\n");
            inTreesBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "   ");
            }
        /* set Foreign (inForeignBlock) *******************************************************/
        else
            {
            MrBayesPrint ("   Skipping \"%s\" block\n", tkn);
            inForeignBlock = YES;
            expecting = Expecting(SEMICOLON);
            strcpy (spacer, "");
            }
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoBreaks (void)
{
    int         i, numBreaks;
    
    numBreaks = 0;
    for (i=0; i<numChar; i++)
        {
        if (charInfo[i].bigBreakAfter == YES)
            {
            numBreaks++;
            }
        }
        
    if (numBreaks > 0)
        {
        if (numBreaks == 1)
            MrBayesPrint ("%s   One data break found after character ", spacer, numBreaks);
        else
            MrBayesPrint ("%s   %d data breaks found after characters: ", spacer, numBreaks);
        for (i=0; i<numChar; i++)
            {
            if (charInfo[i].bigBreakAfter == YES)
                {
                MrBayesPrint ("%d ", i+1);
                }
            }
        MrBayesPrint ("\n");

        if (numBreaks == 1)
            MrBayesPrint ("%s   Successfully defined one break in data\n", spacer);
        else
            MrBayesPrint ("%s   Successfully defined %d breaks in data\n", spacer, numBreaks);
        }
    else
        {
        MrBayesPrint ("%s   No breaks in data found\n", spacer);
        }
        
    return (NO_ERROR);
}


int DoBreaksParm (char *parmName, char *tkn)
{
    int     i, tempInt;
        
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before you can define breaks in the data\n", spacer);
        return (ERROR);
        }
            
    if (expecting == Expecting(NUMBER))
        {
        sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            for (i=0; i<numChar; i++)
                charInfo[i].bigBreakAfter = NO;
            return (ERROR);
            }
        if (tempInt == numChar)
            {
            MrBayesPrint ("%s   Character number %d is the last character. MrBayes will define the\n", spacer, tempInt);
            MrBayesPrint ("%s   break, even though it doesn't make too much sense.\n", spacer);
            }
        tempInt--;
                    
        charInfo[tempInt].bigBreakAfter = YES;
        
        expecting  = (Expecting(NUMBER) | Expecting(SEMICOLON));
        }
    else
        {
        for (i=0; i<numChar; i++)
            charInfo[i].bigBreakAfter = NO;
        return (ERROR);
        }

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoCalibrate (void)
{
    int         i;

    /* show calibration times (for debugging) */
#   if 0
    MrBayesPrint ("Taxon ages\n");
    for (i=0; i<numTaxa; i++)
        MrBayesPrint ("%4d  --  %s\n", i+1, tipCalibration[i].name);
    MrBayesPrint ("Constraint ages\n");
    for (i=0; i<numDefinedConstraints; i++)
        {
        if (definedConstraintsType[i] != HARD)
            continue;
        MrBayesPrint ("%4d  --  %s\n", i+1, nodeCalibration[i].name);
        }
#   endif

    /* Update model if calibrations enforced */
    for (i=0; i<numCurrentDivisions; i++)
        {
        if (!strcmp(modelParams[i].nodeAgePr,"Calibrated"))
            {
            if (SetUpAnalysis (&globalSeed) == ERROR)
                return (ERROR);
            break;
            }
        }

    return (NO_ERROR);
}


int DoCalibrateParm (char *parmName, char *tkn)
{
    static int              isTaxon, paramIndex;
    static char             nodeName[100], calName[100];
    static MrBFlt           priorParams[3];
    static enum CALPRIOR    calPrior;
    int                     howMany, index;
    char                    s[20], tempStr[100];
    MrBFlt                  tempD;
        
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before you can calibrate nodes\n", spacer);
        return (ERROR);
        }
        
    if (expecting == Expecting(PARAMETER))
        {
        if (strcmp(parmName, "Xxxxxxxxxx") != 0)
            {
            MrBayesPrint ("%s   Unexpected error - Wrong parmName in DoCalibrateParm\n", spacer);
            return (ERROR);
            }

        /* find taxon with this name */
        calibrationPtr = NULL;
        howMany = 0;

        /* first look in constraint names */
        if (CheckString (constraintNames, numDefinedConstraints, tkn, &index) != ERROR && definedConstraintsType[index] == HARD)
            {
            calibrationPtr = &nodeCalibration[index];
            howMany++;
            isTaxon = NO;
            strcpy (nodeName, tkn);
            }
        
        /* then look in terminal taxon names */
        if (CheckString (taxaNames, numTaxa, tkn, &index) != ERROR)
            {
            calibrationPtr = &tipCalibration[index];
            howMany++;
            isTaxon = YES;
            strcpy (nodeName, tkn);
            }

        /* return error if not found or ambiguous */
        if (howMany == 0)
            {
            MrBayesPrint ("%s   No taxon or hard constraint named ""%s"" found. Note that only hard constraint can be calibrated.\n", spacer, tkn);
            return (ERROR);
            }
        else if (howMany > 1)
            {
            MrBayesPrint ("%s   Both a taxon and a constraint named ""%s"" encountered -- please rename one\n", spacer, tkn);
            return (ERROR);
            }

        /* get ready to find the equal sign */
        expecting = Expecting(EQUALSIGN);
        }

    else if (expecting == Expecting(EQUALSIGN))
        {
        /* get ready to find the calibration prior */
        expecting = Expecting(ALPHA);
        }

    else if (expecting == Expecting(ALPHA))
        {
        /* set the calibration prior type */
        if (IsArgValid(tkn,tempStr) == NO_ERROR)
            {
            if (!strcmp (tempStr, "Unconstrained"))
                calPrior = unconstrained;
            else if (!strcmp (tempStr, "Fixed"))
                calPrior = fixed;
            else if (!strcmp (tempStr, "Uniform"))
                calPrior = uniform;
            else if (!strcmp (tempStr, "Offsetexponential"))
                calPrior = offsetExponential;
            else if (!strcmp (tempStr, "Truncatednormal"))
                calPrior = truncatedNormal;
            else if (!strcmp (tempStr, "Lognormal"))
                calPrior = logNormal;
            else if (!strcmp (tempStr, "Offsetlognormal"))
                calPrior = offsetLogNormal;
            else if (!strcmp (tempStr, "Gamma"))
                calPrior = standardGamma;
            else if (!strcmp (tempStr, "Offsetgamma"))
                calPrior = offsetGamma;

            if (calPrior == unconstrained)
                {
                /* reset the values of the calibration */
                MrBayesPrint ("%s   Resetting previous calibration for ""%s""\n", spacer, nodeName);

                calibrationPtr->prior           = defaultCalibration.prior;
                calibrationPtr->priorParams[0]  = defaultCalibration.priorParams[0];
                calibrationPtr->priorParams[1]  = defaultCalibration.priorParams[1];
                calibrationPtr->priorParams[2]  = defaultCalibration.priorParams[2];
                calibrationPtr->LnPriorProb     = defaultCalibration.LnPriorProb;
                calibrationPtr->LnPriorRatio    = defaultCalibration.LnPriorRatio;
                calibrationPtr->min             = defaultCalibration.min;
                calibrationPtr->max             = defaultCalibration.max;
                strcpy(calibrationPtr->name, defaultCalibration.name);
            
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                {
                strcpy (calName, tempStr);
                paramIndex = 0;
                priorParams[0] = priorParams[1] = priorParams[2] =  -1.0;
                expecting = Expecting(LEFTPAR);
                }
            }
        else
            {
            MrBayesPrint ("%s   Invalid calibration prior argument \n", spacer);
            return (ERROR);
            }
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        strcat (calName, "(");
        expecting  = Expecting(NUMBER);
        }
    else if (expecting == Expecting(NUMBER))
        {
        sscanf (tkn, "%lf", &tempD);
        if (paramIndex == 0)
            {
            if (calPrior == logNormal)
                {
                if (tempD < 0.0)
                    {
                    MrBayesPrint ("%s   Mean age must be nonnegative\n", spacer);
                    MrBayesPrint ("%s   Parameters of the lognormal distribution used for dating are mean age\n", spacer);
                    MrBayesPrint ("%s   and standard deviation, both specified on the linear scale, not as log values.\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == standardGamma)
                {
                if (tempD <= 0.0)
                    {
                    MrBayesPrint ("%s   Mean parameter must be positive\n", spacer);
                    MrBayesPrint ("%s   Parameters of the gamma distribution used for dating are mean age and\n", spacer);
                    MrBayesPrint ("%s   standard deviation. In terms of the common shape (alpha) and rate (beta)\n", spacer);
                    MrBayesPrint ("%s   parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    MrBayesPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    return (ERROR);
                    }
                }
            else if (tempD < 0.0)
                {
                if (calPrior == fixed)
                    MrBayesPrint ("%s   Fixed age must be nonnegative\n", spacer);
                else if (calPrior == uniform)
                    {
                    MrBayesPrint ("%s   Minimum age must be nonnegative\n", spacer);
                    MrBayesPrint ("%s   Parameters of the uniform are minimum age and maximum age.\n", spacer);
                    }
                else if (calPrior == truncatedNormal)
                    {
                    MrBayesPrint ("%s   Offset (minimum or truncation) age must be nonnegative.\n", spacer);
                    MrBayesPrint ("%s   Parameters of the truncated normal distribution are offset (minimum\n", spacer);
                    MrBayesPrint ("%s   or truncation) age, mean age and standard deviation.\n", spacer);
                    }
                else if (calPrior == offsetGamma)
                    {
                    MrBayesPrint ("%s   Offset age must be nonnegative\n", spacer);
                    MrBayesPrint ("%s   Parameters of the offset gamma distribution used for dating are offset age,\n", spacer);
                    MrBayesPrint ("%s   mean age, and standard deviation. In terms of the common shape (alpha) and\n", spacer);
                    MrBayesPrint ("%s   rate (beta) parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    MrBayesPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    }
                else if (calPrior == offsetExponential)
                    {
                    MrBayesPrint ("%s   Offset age must be nonnegative\n", spacer);
                    MrBayesPrint ("%s   Parameters of the offset exponential are offset age and mean age.\n", spacer);
                    }
                else if (calPrior == offsetLogNormal)
                    {
                    MrBayesPrint ("%s   Offset age must be nonnegative\n", spacer);
                    MrBayesPrint ("%s   Parameters of the offset lognormal distribution are offset age, mean age,\n", spacer);
                    MrBayesPrint ("%s   and standard deviation. All values are specified on the linear scale, not\n", spacer);
                    MrBayesPrint ("%s   as log values.\n", spacer);
                    }
                return (ERROR);
                }
            priorParams[0] = tempD;
            if (calPrior == fixed)
                expecting = Expecting(RIGHTPAR);
            else
                expecting = Expecting(COMMA);
            }
        else if (paramIndex == 1)
            {
            if (calPrior == uniform)
                {
                if (tempD <= priorParams[0])
                    {
                    MrBayesPrint ("%s   Maximum age of uniform distribution must be larger than minimum age\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == offsetExponential)
                {
                if (tempD <= priorParams[0])
                    {
                    MrBayesPrint ("%s   Mean age must be larger than offset age.\n", spacer);
                    MrBayesPrint ("%s   MrBayes now uses offset and mean rather than offset and rate\n", spacer);
                    MrBayesPrint ("%s   as the parameters for the offset exponential distribution.\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == truncatedNormal)
                {
                if (tempD <= priorParams[0])
                    {
                    MrBayesPrint ("%s   Mean age must be larger than offset (truncation) age.\n", spacer);
                    MrBayesPrint ("%s   Parameters of the truncated normal distribution are offset (minimum\n", spacer);
                    MrBayesPrint ("%s   or truncation) age, mean age and standard deviation\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == logNormal)
                {
                if (tempD <= 0.0)
                    {
                    MrBayesPrint ("%s   Standard deviation must be positive.\n", spacer);
                    MrBayesPrint ("%s   Parameters of the lognormal distribution used for dating are mean age\n", spacer);
                    MrBayesPrint ("%s   and standard deviation, both specified on the linear scale, not as log values.\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == offsetLogNormal)
                {
                if (tempD <= priorParams[0])
                    {
                    MrBayesPrint ("%s   Mean age must be larger than offset age.\n", spacer);
                    MrBayesPrint ("%s   Parameters of the offset lognormal distribution are offset age, mean age,\n", spacer);
                    MrBayesPrint ("%s   and standard deviation. All values are specified on the linear scale, not\n", spacer);
                    MrBayesPrint ("%s   as log values.\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == standardGamma)
                {
                if (tempD <= 0.0)
                    {
                    MrBayesPrint ("%s   Standard deviation must be positive.\n", spacer);
                    MrBayesPrint ("%s   Parameters of the gamma distribution used for dating are mean age and\n", spacer);
                    MrBayesPrint ("%s   standard deviation. In terms of the common shape (alpha) and rate (beta)\n", spacer);
                    MrBayesPrint ("%s   parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    MrBayesPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == offsetGamma)
                {
                if (tempD <= 0.0)
                    {
                    MrBayesPrint ("%s   Mean age must be positive.\n", spacer);
                    MrBayesPrint ("%s   Parameters of the offset gamma distribution used for dating are offset age,\n", spacer);
                    MrBayesPrint ("%s   mean age, and standard deviation. In terms of the common shape (alpha) and\n", spacer);
                    MrBayesPrint ("%s   rate (beta) parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    MrBayesPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    return (ERROR);
                    }
                }

            priorParams[1] = tempD;
            if (calPrior == uniform || calPrior == standardGamma || calPrior == logNormal || calPrior == offsetExponential)
                expecting = Expecting(RIGHTPAR);
            else
                expecting = Expecting(COMMA);
            }
        else /* if (paramIndex == 2) */
            {
            if (calPrior == offsetGamma)
                {
                if (tempD <= 0.0)
                    {
                    MrBayesPrint ("%s   Standard deviation must be positive.\n", spacer);
                    MrBayesPrint ("%s   Parameters of the offset gamma distribution used for dating are offset age,\n", spacer);
                    MrBayesPrint ("%s   mean age, and standard deviation. In terms of the common shape (alpha) and\n", spacer);
                    MrBayesPrint ("%s   rate (beta) parameterization, the expected mean is alpha/beta and the standard\n", spacer);
                    MrBayesPrint ("%s   deviation is the square root of (alpha / beta^2).\n", spacer);
                    return (ERROR);
                    }
                }
            else if (calPrior == offsetLogNormal)
                {
                if (tempD <= 0.0)
                    {
                    MrBayesPrint ("%s   Standard deviation must be positive.\n", spacer);
                    MrBayesPrint ("%s   Parameters of the offset lognormal distribution are offset age, mean age,\n", spacer);
                    MrBayesPrint ("%s   and standard deviation. All values are specified on the linear scale, not\n", spacer);
                    MrBayesPrint ("%s   as log values.\n", spacer);
                    return (ERROR);
                    }
                }
            priorParams[2] = tempD;
            expecting = Expecting(RIGHTPAR);
            }
        sprintf (s, "%1.2lf", tempD);
        strcat (calName, s);
        }
    else if (expecting == Expecting(COMMA))
        {
        strcat (calName, ",");
        paramIndex++;
        expecting  = Expecting(NUMBER);
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        strcat (calName, ")");
        if (isTaxon == YES)
            MrBayesPrint ("%s   Setting age of taxon '%s' to %s\n", spacer, nodeName, calName);
        else
            MrBayesPrint ("%s   Setting age of constraint node '%s' to %s\n", spacer, nodeName, calName);

        /* set calibration based on collected values and settings */
        strcpy(calibrationPtr->name, calName);
        calibrationPtr->priorParams[0]  = priorParams[0];
        calibrationPtr->priorParams[1]  = priorParams[1];
        calibrationPtr->priorParams[2]  = priorParams[2];
        calibrationPtr->prior           = calPrior;
        if (calPrior == fixed)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbFix;
            calibrationPtr->LnPriorRatio    = &LnProbRatioFix;
            calibrationPtr->min             = priorParams[0];
            calibrationPtr->max             = priorParams[0];
            }
        else if (calPrior == uniform)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbUniform;
            calibrationPtr->LnPriorRatio    = &LnProbRatioUniform;
            calibrationPtr->min             = priorParams[0];
            calibrationPtr->max             = priorParams[1];
            }
        else if (calPrior == offsetExponential)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbOffsetExponential_Param_Offset_Mean;
            calibrationPtr->LnPriorRatio    = &LnProbRatioOffsetExponential_Param_Offset_Mean;
            calibrationPtr->min             = priorParams[0];
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == truncatedNormal)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbTruncatedNormal_Param_Trunc_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioTruncatedNormal_Param_Trunc_Mean_Sd;
            calibrationPtr->min             = priorParams[0];
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == logNormal)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbLognormal_Param_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioLognormal_Param_Mean_Sd;
            calibrationPtr->min             = BRLENS_MIN;
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == offsetLogNormal)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbOffsetLognormal_Param_Offset_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioOffsetLognormal_Param_Offset_Mean_Sd;
            calibrationPtr->min             = BRLENS_MIN + priorParams[0];
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == standardGamma)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbGamma_Param_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioGamma_Param_Mean_Sd;
            calibrationPtr->min             = BRLENS_MIN;
            calibrationPtr->max             = POS_INFINITY;
            }
        else if (calPrior == offsetGamma)
            {
            calibrationPtr->LnPriorProb     = &LnPriorProbOffsetGamma_Param_Offset_Mean_Sd;
            calibrationPtr->LnPriorRatio    = &LnProbRatioOffsetGamma_Param_Offset_Mean_Sd;
            calibrationPtr->min             = BRLENS_MIN + priorParams[0];
            calibrationPtr->max             = POS_INFINITY;
            }

        /* get ready to find more calibrated nodes or taxa, if present */
        expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoCharset (void)
{
    /* first add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* add name to charSetNames */
    if (AddString (&charSetNames, numCharSets, tempSetName) == ERROR)
        {
        MrBayesPrint ("%s   Problem adding charset %s to list\n", spacer, tempSetName);
        return (ERROR);
        }

    /* store charset */
    AddBitfield (&charSet, numCharSets, tempSet, numChar);

    /* increment number of char sets */
    numCharSets++;

    return (NO_ERROR);
}


int DoCharsetParm (char *parmName, char *tkn)
{
    int     i, index, tempInt, allDigit;
    
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before charsets can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check that the name of the charset is not a number */
            allDigit = YES;
            for (i=0; i<(int)strlen(tkn); i++)
                {
                if (tkn[i] == '0' || tkn[i] == '1' || tkn[i] == '2' || tkn[i] == '3' || tkn[i] == '4' || 
                    tkn[i] == '5' || tkn[i] == '6' || tkn[i] == '7' || tkn[i] == '8' || tkn[i] == '9' || tkn[i] == '.')
                    {}
                else
                    allDigit = NO;
                }
            if (allDigit == YES)
                {
                MrBayesPrint ("%s   Charset name may not be a number\n", spacer);
                return (ERROR);
                }
            
            /* check size of charset name */
            if (strlen(tkn) > 99)
                {
                MrBayesPrint ("%s   Charset name is too long\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a charset */
            if (numCharSets > 1)
                {
                if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
                    {
                    /* if the charset name has not been used, then we should have an ERROR returned */
                    /* we _want_ to be here */

                    }
                else
                    {
                    MrBayesPrint ("%s   Charset name has been used previously\n", spacer);
                    return (ERROR);
                    }
                }
                
            /* add the name to the character set */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numChar; i++)
                tempSet[i] = 0;
            
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            MrBayesPrint ("%s   Defining charset called '%s'\n", spacer, tkn);
            expecting = Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        /* We are defining a character set in terms of another (called tkn, here). We should be able
           to find tkn in the list of character set names. If we cannot, then we have a problem and
           return an error. */
        if (numCharSets < 1)
            {
            MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
            return (ERROR);
            }
        if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
            {
            MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
            return (ERROR);
            }
        /* add characters from charset "tkn" to new tempset */
        for (i=0; i<numChar; i++)
            {
            if (IsBitSet(i,charSet[index]) == YES)
                tempSet[i] = 1;
            }       
        fromI = toJ = everyK = -1;

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numChar;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
                    {
                    return (ERROR);
                    }
                }
                
            }

        
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoCharStat (void)
{
    int         i, j, numDivs;
    char        tempName[100];
    
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }
            
    if (numDefinedPartitions == 1)
        MrBayesPrint ("%s   1 character partition defined:\n", spacer, numDefinedPartitions);
    else
        MrBayesPrint ("%s   %d character partitions defined:\n", spacer, numDefinedPartitions);
    for (i=0; i<numDefinedPartitions; i++)
        {
        numDivs = GetNumPartDivisions (i);
        if (numDivs == 1)
            MrBayesPrint ("%s      Partition %d (\"%s\") does not divide the characters\n", spacer, i+1, partitionNames[i]);
        else
            MrBayesPrint ("%s      Partition %d (\"%s\") divides the characters into %d parts\n", spacer, i+1, partitionNames[i], numDivs);
        }
    MrBayesPrint ("%s      Current partition is \"%s\"\n", spacer, partitionNames[partitionNum]);
    MrBayesPrint ("\n");

    /* print out list of characters with information about each */
    MrBayesPrint ("%s   Showing character status:\n\n", spacer);
    MrBayesPrint ("%s                                                    Partition(s)\n", spacer);
    MrBayesPrint ("%s      #      Type      In/Out    Ambiguity Order  ", spacer);
    for (i=0; i<numDefinedPartitions; i++)
        MrBayesPrint (" %2d", i+1);
    MrBayesPrint ("\n");
    MrBayesPrint ("%s   -----------------------------------------------", spacer);
    for (i=0; i<numDefinedPartitions; i++)
        MrBayesPrint ("---");
    MrBayesPrint ("\n");
    for (i=0; i<numChar; i++)
        {
        MrBayesPrint ("%s   %4d -- ", spacer, i+1);
                
        if (charInfo[i].charType == DNA)
            MrBayesPrint ("   DNA");
        else if (charInfo[i].charType == RNA)
            MrBayesPrint ("   RNA");
        else if (charInfo[i].charType == PROTEIN)
            MrBayesPrint ("  Prot");
        else if (charInfo[i].charType == RESTRICTION)
            MrBayesPrint ("  Rest");
        else if (charInfo[i].charType == STANDARD)
            MrBayesPrint (" Stand");
        else if (charInfo[i].charType == CONTINUOUS)
            MrBayesPrint ("  Cont");
            
        if (charInfo[i].charType == DNA)
            MrBayesPrint ("   4");
        else if (charInfo[i].charType == RNA)
            MrBayesPrint ("   4");
        else if (charInfo[i].charType == PROTEIN)
            MrBayesPrint ("  20");
        else if (charInfo[i].charType == RESTRICTION)
            MrBayesPrint ("   2");
        else if (charInfo[i].charType == STANDARD)
            MrBayesPrint ("  %2d", charInfo[i].numStates);
        else if (charInfo[i].charType == CONTINUOUS)
            MrBayesPrint (" Inf");
            
        if (charInfo[i].isExcluded == NO)
            MrBayesPrint ("  Included");
        else
            MrBayesPrint ("  Excluded");
            
        if (charInfo[i].isMissAmbig == YES)
            MrBayesPrint ("  MissAmbig");
        else
            MrBayesPrint ("       None");
            
        if (charInfo[i].ctype == UNORD)
            MrBayesPrint (" Unord");
        else if (charInfo[i].ctype == ORD)
            MrBayesPrint ("   Ord");
        else if (charInfo[i].ctype == DOLLO)
            MrBayesPrint (" Dollo");
        else if (charInfo[i].ctype == IRREV)
            MrBayesPrint (" Irrev");

        MrBayesPrint ("  ");
            
        for (j=0; j<numDefinedPartitions; j++)
            MrBayesPrint (" %2d", partitionId[i][j]);

        /* MrBayesPrint ("%4d   ", charSet[i]);*/
        
        if (charInfo[i].pairsId > 0)
            {
            /* find paired character */
            for (j=0; j<numChar; j++)
                {
                if (i != j && charInfo[j].pairsId == charInfo[i].pairsId)
                    {
                    MrBayesPrint (" (coupled with %d)", j+1);
                    break;
                    }
                }
            }
                    
        MrBayesPrint ("\n");
        
        if (charInfo[i].bigBreakAfter == YES)
            {
            MrBayesPrint ("%s   ", spacer);
            MrBayesPrint ("     - - - - - - - - - - - - - - - - - - - -  \n");
            }
        
        /* we may want to pause */
        if (autoClose == NO)
            {
            if ((i+1) % 100 == 0)
                {
                MrBayesPrint ("%s   Hit return key to continue  ", spacer);
                fflush (stdin);
                if (fgets (tempName, 100, stdin) == NULL)
                    {
                    printf ("Error in function: %s at line: %d in file: %s", __func__, __LINE__, __FILE__);
                    }
                }
            }
        }

    return (NO_ERROR);
}


int DoCitations (void)
{
    MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
    MrBayesPrint ("   Citations                                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   If you publish results obtained using MrBayes you may want to cite the        \n");
    MrBayesPrint ("   program using one of these papers:                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Huelsenbeck, J. P. and F. Ronquist. 2001. MRBAYES: Bayesian                \n");
    MrBayesPrint ("         inference of phylogeny. Bioinformatics 17:754-755.                      \n");
    MrBayesPrint ("      Ronquist, F. and J. P. Huelsenbeck. 2003. MRBAYES 3: Bayesian phylogenetic \n");
    MrBayesPrint ("         inference under mixed models. Bioinformatics 19:1572-1574.              \n");
    MrBayesPrint ("      Ronquist, F. et al. 2012. MRBAYES 3.2: Efficient Bayesian phylogenetic     \n");
    MrBayesPrint ("         inference and model selection across a large model space. Systematic    \n");
    MrBayesPrint ("         Biology 61 (in press).                                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   If you use the parallel abilities of the program, you may also want to cite   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Altekar, G., S. Dwarkadas, J. P. Huelsenbeck, and F. Ronquist. 2004.       \n");
    MrBayesPrint ("         Parallel Metropolis-coupled Markov chain Monte Carlo for Bayesian       \n");
    MrBayesPrint ("         phylogenetic inference. Bioinformatics 20:407-415.                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   If you use the BEAGLE library, the appropriate citation is                    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Ayres, D. L., A. Darling, D. J. Zwickl, P. Beerli, M. T. Holder, P. O.     \n");
    MrBayesPrint ("         J. P. Huelsenbeck, F. Ronquist, D. L. Swofford, M. P. Cummings, A.      \n");
    MrBayesPrint ("         Rambaut, and M. A. Suchard. 2012. BEAGLE: an application programming    \n");
    MrBayesPrint ("         interface for statistical phylogenetics. Syst. Biol. 61:170-173.        \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   You should also cite other papers for different ideas that are implemented    \n");
    MrBayesPrint ("   in the program. For example, the program performs Bayesian inference of       \n");
    MrBayesPrint ("   phylogeny, an idea that was first proposed in the following papers:           \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Larget, B., and D. Simon. 1999. Markov chain Monte Carlo                   \n");
    MrBayesPrint ("         algorithms for the Bayesian analysis of phylogenetic trees.             \n");
    MrBayesPrint ("         Mol. Biol. Evol. 16:750-759.                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Li, S. 1996. Phylogenetic tree construction using Markov chain             \n");
    MrBayesPrint ("         Monte carlo. Ph. D. dissertation, Ohio State University, Columbus.      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mau, B. 1996. Bayesian phylogenetic inference via Markov chain             \n");
    MrBayesPrint ("         Monte carlo methods. Ph. D. dissertation, University of                 \n");
    MrBayesPrint ("         Wisconsin, Madison.                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mau, B., and M. Newton. 1997. Phylogenetic inference for binary            \n");
    MrBayesPrint ("         data on dendrograms using Markov chain Monte Carlo. Journal of          \n");
    MrBayesPrint ("         Computational and Graphical Statistics 6:122-131.                       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mau, B., M. Newton, and B. Larget. 1999. Bayesian phylogenetic             \n");
    MrBayesPrint ("         inference via Markov chain Monte carlo methods. Biometrics. 55:1-12.    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Newton, M., B. Mau, and B. Larget. 1999. Markov chain Monte Carlo          \n");
    MrBayesPrint ("         for the Bayesian analysis of evolutionary trees from aligned            \n");
    MrBayesPrint ("         molecular sequences. In Statistics in molecular biology (F. Seillier-   \n");
    MrBayesPrint ("         Moseiwitch, T. P. Speed, and M. Waterman, eds.). Monograph Series       \n");
    MrBayesPrint ("         of the Institute of Mathematical Statistics.                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Rannala, B., and Z. Yang. 1996. Probability distribution of                \n");
    MrBayesPrint ("         molecular evolutionary trees: a new method of phylogenetic              \n");
    MrBayesPrint ("         inference. J. Mol. Evol. 43:304-311.                                    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z., and B. Rannala. 1997. Bayesian phylogenetic inference            \n");
    MrBayesPrint ("         using DNA sequences: a Markov chain Monte carlo method. Molecular       \n");
    MrBayesPrint ("         Biology and Evolution. 14:717-724.                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes uses Markov chain Monte Carlo (MCMC) to approximate the posterior     \n");
    MrBayesPrint ("   probability of trees. MCMC was developed in the following papers:             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Metropolis, N., A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller,          \n");
    MrBayesPrint ("         and E. Teller. 1953. Equations of state calculations by fast            \n");
    MrBayesPrint ("         computing machines. J. Chem. Phys. 21:1087-1091.                        \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Hastings, W. K. 1970. Monte Carlo sampling methods using Markov            \n");
    MrBayesPrint ("         chains and their applications. Biometrika 57:97-109.                    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   In particular, MrBayes implements a variant of MCMC that was described by     \n");
    MrBayesPrint ("   Charles Geyer:                                                                \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Geyer, C. J. 1991. Markov chain Monte Carlo maximum likelihood.            \n");
    MrBayesPrint ("         Pages 156-163 in Computing Science and Statistics: Proceed-             \n");
    MrBayesPrint ("         ings of the 23rd Symposium on the Interface. (E. M. Keramidas,          \n");
    MrBayesPrint ("         ed.). Fairfax Station: Interface Foundation.                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes implements a large number of DNA substitution models. These models    \n");
    MrBayesPrint ("   are of three different structures. The \"4by4\" models are the usual          \n");
    MrBayesPrint ("   flavor of phylogenetic models. The \"Doublet\" model was first proposed       \n");
    MrBayesPrint ("   by                                                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Schoniger, M., and A. von Haeseler. 1994. A stochastic model and the       \n");
    MrBayesPrint ("         evolution of autocorrelated DNA sequences. Molecular Phylogenetics      \n");
    MrBayesPrint ("         and Evolution 3:240-247.                                                \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The program also implements codon models. Two papers, published back-to-back  \n");
    MrBayesPrint ("   were the first to implement a codon model of DNA substitution in which the    \n");
    MrBayesPrint ("   substitution process is modelled on the codon, not on a site-by-site basis:   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Goldman, N., and Z. Yang. 1994. A codon-based model of nucleotide          \n");
    MrBayesPrint ("         substitution for protein coding DNA sequences. Molecular Biology        \n");
    MrBayesPrint ("         and Evolution. 11:725-736.                                              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Muse, S., and B. Gaut. 1994. A likelihood approach for comparing           \n");
    MrBayesPrint ("         synonymous and non-synonymous substitution rates, with application      \n");
    MrBayesPrint ("         to the chloroplast genome. Molecular Biology and Evolution.             \n");
    MrBayesPrint ("         11:715-724.                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The program can be used to detect positively slected amino-acid sites using   \n");
    MrBayesPrint ("   a full hierarchical Bayes analysis. The method is based on the excellent paper\n");
    MrBayesPrint ("   by Nielsen and Yang:                                                          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Nielsen, R., and Z. Yang. 1998. Likelihood models for detecting            \n");
    MrBayesPrint ("         positively selected amino acid sites and applications to the HIV-1      \n");
    MrBayesPrint ("         envelope gene. Genetics. 148:929-936.                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The previous four papers describe three different stuctures for the nuc-      \n");
    MrBayesPrint ("   leotide models implemented in MrBayes--the four-by-four models, the           \n");
    MrBayesPrint ("   16-by-16 (doublet) models and the 64-by-64 (codon) models. The program        \n");
    MrBayesPrint ("   implements three different substitution models within each model structure.   \n");
    MrBayesPrint ("   These include the nst=1 models:                                               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Jukes, T., and C. Cantor. 1969. Evolution of protein molecules.            \n");
    MrBayesPrint ("         Pages 21-132 in Mammalian Protein Metabolism. (H. Munro, ed.).          \n");
    MrBayesPrint ("         Academic Press, New York.                                               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Felsenstein, J. 1981. Evolutionary trees from DNA sequences: A             \n");
    MrBayesPrint ("         maximum likelihood approach. Journal of Molecular Evolution             \n");
    MrBayesPrint ("         17:368-376.                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   the nst=2 models:                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Kimura, M. 1980. A simple method for estimating evolutionary rates         \n");
    MrBayesPrint ("         of base substitutions through comparative studies of nucleotide         \n");
    MrBayesPrint ("         sequences. Journal of Molecular Evolution. 16:111-120.                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Hasegawa, M., T. Yano, and H. Kishino. 1984. A new molecular clock         \n");
    MrBayesPrint ("         of mitochondrial DNA and the evolution of Hominoids. Proc.              \n");
    MrBayesPrint ("         Japan Acad. Ser. B 60:95-98.                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Hasegawa, M., H. Kishino, and T. Yano. 1985. Dating the human-ape          \n");
    MrBayesPrint ("         split by a molecular clock of mitochondrial DNA. Journal of             \n");
    MrBayesPrint ("         Molecular Evolution 22:160-174.                                         \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   and the the nst=6 models:                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Tavare, S. 1986. Some probabilistic and statisical problems on the         \n");
    MrBayesPrint ("         analysis of DNA sequences. Lect. Math. Life Sci. 17:57-86.              \n");
    MrBayesPrint ("         17:368-376.                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes implements a large number of amino-acid models. These include:        \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Poisson --                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Bishop, M.J., and A.E. Friday. 1987. Tetropad relationships: the           \n");
    MrBayesPrint ("         molecular evidence. Pp. 123-139 in Molecules and morphology in          \n");
    MrBayesPrint ("         evolution: conflict or compromise? (C. Patterson, ed.). Cambridge       \n");
    MrBayesPrint ("         University Press, Cambridge, England.                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Jones --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Jones, D.T., W. R. Taylor, and J. M. Thornton. 1992. The rapid generation  \n");
    MrBayesPrint ("         of mutation data matrices from protein sequences. Comput. Appl.         \n");
    MrBayesPrint ("         Biosci. 8:275-282.                                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Dayhoff --                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Dayhoff, M.O., R.M. Schwartz, and B.C. Orcutt. 1978. A model of evol-      \n");
    MrBayesPrint ("         utionary change in proteins. Pp. 345-352 in Atlas of protein sequence   \n");
    MrBayesPrint ("         and structure. Vol. 5, Suppl. 3. National Biomedical Research           \n");
    MrBayesPrint ("          Foundation, Washington, D.C.                                           \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mtrev --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Adachi, J. and M. Hasegawa. 1996. MOLPHY version 2.3: programs for         \n");
    MrBayesPrint ("         molecular phylogenetics based on maximum likelihood.  Computer Science  \n");
    MrBayesPrint ("         Monographs of Institute of Statistical Mathematics 28:1-150.            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Mtmam --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Cao, Y., A. Janke, P.J. Waddell, M. Westerman, O. Takenaka, S. Murata,     \n");
    MrBayesPrint ("         N. Okada, S. Paabo, and M. Hasegawa. 1998. Conflict amongst individual  \n");
    MrBayesPrint ("         mitochondrial proteins in resolving the phylogeny of eutherian orders.  \n");
    MrBayesPrint ("         Journal of Molecular Evolution                                          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z., R. Nielsen, and M. Hasegawa. 1998.  Models of amino acid         \n");
    MrBayesPrint ("         substitution and applications to mitochondrial protein evolution        \n");
    MrBayesPrint ("         Molecular Biology and Evolution 15:1600-1611.                           \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      WAG --                                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Whelan, S. and Goldman, N. 2001. A general empirical model of protein      \n");
    MrBayesPrint ("         evolution derived from multiple protein families using a maximum-       \n");
    MrBayesPrint ("         likelihood approach. Molecular Biology and Evolution 18:691-699.        \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Rtrev --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Dimmic M.W., J.S. Rest, D.P. Mindell, and D. Goldstein. 2002. RArtREV:     \n");
    MrBayesPrint ("         An amino acid substitution matrix for inference of retrovirus and       \n");
    MrBayesPrint ("         reverse transcriptase phylogeny. Journal of Molecular Evolution         \n");
    MrBayesPrint ("         55: 65-73.                                                              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Cprev --                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Adachi, J., P. Waddell, W. Martin, and M. Hasegawa. 2000. Plastid          \n");
    MrBayesPrint ("         genome phylogeny and a model of amino acid substitution for proteins    \n");
    MrBayesPrint ("         encoded by chloroplast DNA. Journal of Molecular Evolution              \n");
    MrBayesPrint ("         50:348-358.                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Blosum --                                                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Henikoff, S., and J. G. Henikoff. 1992. Amino acid substitution            \n");
    MrBayesPrint ("         matrices from protein blocks. Proc. Natl. Acad. Sci., U.S.A.            \n");
    MrBayesPrint ("         89:10915-10919. The matrix implemented in MrBayes is Blosum62.          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Vt --                                                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Muller, T., and M. Vingron. 2000. Modeling amino acid replacement.         \n");
    MrBayesPrint ("         Journal of Computational Biology 7:761-776.                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      LG --                                                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Le, Si Q. & Gascuel, O. 2008 An improved general amino- acid replacement   \n");
    MrBayesPrint ("         matrix. Mol. Biol. Evol. 25, 1307-1320.                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes implements a simple Jukes-Cantor-like model for restriction sites     \n");
    MrBayesPrint ("   and other binary data. A problem with some of these data is that there is a   \n");
    MrBayesPrint ("   coding bias, such that certain characters are missing from any observable     \n");
    MrBayesPrint ("   data matrix. It is impossible, for instance, to observe restriction sites that\n");
    MrBayesPrint ("   are absent in all the studied taxa. However, MrBayes corrects for this coding \n");
    MrBayesPrint ("   bias according to an idea described in                                        \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Felsenstein, J. 1992. Phylogenies from restriction sites: A maximum-       \n");
    MrBayesPrint ("         likelihood approach. Evolution 46:159-173.                              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The model used by MrBayes for 'standard' or morphological data is based on    \n");
    MrBayesPrint ("   the ideas originally presented by                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Lewis, P. O. 2001. A likelihood approach to estimating phylogeny from      \n");
    MrBayesPrint ("         discrete morphological character data. Systematic Biology 50:913-925.   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   For both DNA sequence and amino-acid data, the program allows rates to        \n");
    MrBayesPrint ("   change under a covarion-like model, first described by Tuffley and Steel      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Tuffley, C., and M. Steel. 1998. Modeling the covarion hypothesis          \n");
    MrBayesPrint ("         of nucleotide substitution. Mathematical Biosciences 147:63-91.         \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   and implemented by Huelsenbeck (2002)                                         \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Huelsenbeck, J. P. 2002. Testing a covariotide model of DNA sub-           \n");
    MrBayesPrint ("         stitution. Molecular Biology and Evolution 19(5):698-707.               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   Galtier (2001) implements a different variant of the covarion model in        \n");
    MrBayesPrint ("   a paper that is worth reading:                                                \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Galtier, N. 2001. Maximum-likelihood phylogenetic analysis under a         \n");
    MrBayesPrint ("         covarion-like model. Mol. Biol. Evol. 18:866-873.                       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   A number of models are available that allow rates to vary                     \n");
    MrBayesPrint ("   across the characters. The program implements the proportion                  \n");
    MrBayesPrint ("   of invariable sites model and two variants of gamma distributed               \n");
    MrBayesPrint ("   rate variation. Yang\'s (1993) paper is a good one to cite for                \n");
    MrBayesPrint ("   implementing a gamma-distributed rates model. In the 1994 paper he            \n");
    MrBayesPrint ("   provides a way to approximate the continuous gamma distribution:              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z. 1993. Maximum likelihood estimation of phylogeny from DNA         \n");
    MrBayesPrint ("         sequences when substitution rates differ over sites. Molecular          \n");
    MrBayesPrint ("         Biology and Evolution 10:1396-1401.                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z. 1994. Maximum likelihood phylogenetic estimation from DNA         \n");
    MrBayesPrint ("         sequences with variable rates over sites: Approximate methods.          \n");
    MrBayesPrint ("         Journal of Molecular Evolution 39:306-314.                              \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The program also implements Yang\'s autocorrelated gamma model. In            \n");
    MrBayesPrint ("   this model, the rate at one site depends to some extent on the rate at        \n");
    MrBayesPrint ("   an adjacent site. The appropriate citation for this model is:                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z. 1995. A space-time process model for the evolution of             \n");
    MrBayesPrint ("         DNA sequences. Genetics 139:993-1005.                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The following two papers show how ancestral states on a tree can be recon-    \n");
    MrBayesPrint ("   structed. The Yang et al. paper implements an empirical Bayes approach while  \n");
    MrBayesPrint ("   Huelsenbeck and Bollback use a pure, hierarchical Bayes approach. The method  \n");
    MrBayesPrint ("   used in MrBayes is the latter, since it integrates over uncertainty in model  \n");
    MrBayesPrint ("   parameters.                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z., S. Kumar, and M. Nei. 1995. A new method of inference of         \n");
    MrBayesPrint ("         ancestral nucleotide and amino acid sequences. Genetics 141:1641        \n");
    MrBayesPrint ("         1650.                                                                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Huelsenbeck, J. P., and J. P. Bollback. 2001. Empirical and hier-          \n");
    MrBayesPrint ("         archical Bayesian estimation of ancestral states. Systematic            \n");
    MrBayesPrint ("         Biology 50:351-366.                                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   You may also want to consult a more recent review of Bayesian reconstruction  \n");
    MrBayesPrint ("   of ancestral states and character evolution:                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Ronquist, F. 2004. Bayesian inference of character evolution. Trends in    \n");
    MrBayesPrint ("         Ecology and Evolution 19: 475-481.                                      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes allows you to analyze gene tree - species tree problems using the     \n");
    MrBayesPrint ("   multi-species coalescent approach originally proposed by Edwards et al:       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Edwards, S., L. Liu, and D. Pearl. 2007. High-resolution species trees     \n");
    MrBayesPrint ("         without concatenation. Proc. Natl. Acad. Sci. USA 104: 5936-5941.       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The program implements an incredibly parameter rich model, first described    \n");
    MrBayesPrint ("   by Tuffley and Steel (1997), that orders trees in the same way as the         \n");
    MrBayesPrint ("   so-called parsimony method of phylogenetic inference. The appropriate         \n");
    MrBayesPrint ("   citation is:                                                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Tuffley, C., and M. Steel. 1997. Links between maximum likelihood          \n");
    MrBayesPrint ("         and maximum parsimony under a simple model of site substitution.        \n");
    MrBayesPrint ("         Bull. Math. Bio. 59:581-607.                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   MrBayes implements three relaxed clock models: the Compound Poisson Process   \n");
    MrBayesPrint ("   (CPP), the Thorne-Kishino 2002 (TK02), and the Independent Gamma Rates (IGR)  \n");
    MrBayesPrint ("   models. The CPP model was first described by Huelsenbeck et al. (2000). It    \n");
    MrBayesPrint ("   is an autocorrelated discrete model of rate variation over time. Instead of   \n");
    MrBayesPrint ("   the modified gamma distribution originally proposed for the rate multipliers, \n");
    MrBayesPrint ("   MrBayes uses a lognormal distribution. The extensions necessary to sample over\n");
    MrBayesPrint ("   tree space under this model are original to MrBayes; the original paper only  \n");
    MrBayesPrint ("   considered fixed trees.                                                       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The TK02 model was first described by Thorne and Kishino (2002), and is a     \n");
    MrBayesPrint ("   variant of a model presented by them earlier (Thorne et al., 1998). It is an  \n");
    MrBayesPrint ("   autocorrelated continuous model, in which rates vary according to a lognormal \n");
    MrBayesPrint ("   distribution. Specifically, the rate of a descendant node is assumed to be    \n");
    MrBayesPrint ("   drawn from a lognormal distribution with the mean being the rate of the an-   \n");
    MrBayesPrint ("   cestral node, and the variance being proportional to the length of the branch \n");
    MrBayesPrint ("   separating the nodes (measured in terms of expected substitutions per site at \n");
    MrBayesPrint ("   the base rate of the clock).                                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The final relaxed clock model is the IGR model, in which branch rates are     \n");
    MrBayesPrint ("   modeled as being drawn independently from gamma distributions. The model was  \n");
    MrBayesPrint ("   originally described in the literature as the 'White Noise' model by Lepage   \n");
    MrBayesPrint ("   et al. (2007), but the original MrBayes implementation predates that paper.   \n");
    MrBayesPrint ("   The IGR model is closely related to the uncorrelated gamma model presented    \n");
    MrBayesPrint ("   originally by Drummond et al. (2006), but it is more elegant in that it truly \n");
    MrBayesPrint ("   lacks time structure. See Lepage et al. (2007) for details.                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Huelsenbeck, J. P., B. Larget, and D. Swofford. 2000. A compound Poisson   \n");
    MrBayesPrint ("         process for relaxing the molecular clock. Genetics 154: 1879-1892.      \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Thorne, J. L., H. Kishino, and I. S. Painter. 1998. Estimating the rate    \n");
    MrBayesPrint ("         of evolution of the rate of molecular evolution. Mol. Biol. Evol.       \n");
    MrBayesPrint ("         15: 1647-1657.                                                          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Thorne, J. L., and H. Kishino. 2002. Divergence time and evolutionary      \n");
    MrBayesPrint ("         rate estimation with multilocus data. Syst. Biol. 51: 689-702.          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Drummond, A. J., S. Y. W. Ho, M. J. Phillips, and A. Rambaut. 2006.        \n");
    MrBayesPrint ("         Relaxed phylogenetics and dating with confidence. PLoS Biology          \n");
    MrBayesPrint ("         4: 699-710.                                                             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Lepage, T., D. Bryant, H. Philippe, and N. Lartillot. 2007. A general      \n");
    MrBayesPrint ("         comparison of relaxed molecular clock models. Mol. Biol. Evol.          \n");
    MrBayesPrint ("         24: 2669-2680.                                                          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The standard tree proposals used by MrBayes are described by Lakner et al.    \n");
    MrBayesPrint ("   (2008). The parsimony-biased tree proposals are still undescribed, although   \n");
    MrBayesPrint ("   a rough outline of the idea is presented in the same paper.                   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Lakner, C., P. van der Mark, J. P. Huelsenbeck, B. Larget, and F. Ronquist.\n");
    MrBayesPrint ("         2008. Efficiency of Markov chain Monte Carlo tree proposals in Bayesian \n");
    MrBayesPrint ("         phylogenetics. Syst. Biol. 57: 86-103.                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The topology convergence diagnostic used by MrBayes, the average standard     \n");
    MrBayesPrint ("   deviation of split frequencies, is described by Lakner et al. (2008). The     \n");
    MrBayesPrint ("   potential scale reduction factor, the diagnostic used by MrBayes for contin-  \n");
    MrBayesPrint ("   uous parameters, was first proposed by Gelman and Rubin (1992). The auto-     \n");
    MrBayesPrint ("   tuning mechanism used in MrBayes is based on a paper by Roberts and Rosenthal \n");
    MrBayesPrint ("   (2009).                                                                       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Gelman, A., and D. B. Rubin. 1992. Inference from iterative simulation     \n");
    MrBayesPrint ("         using multiple sequences. Statistical Science 7: 457-472.               \n");
    MrBayesPrint ("         Bull. Math. Bio. 59:581-607.                                            \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Lakner, C., P. van der Mark, J. P. Huelsenbeck, B. Larget, and F. Ronquist.\n");
    MrBayesPrint ("         2008. Efficiency of Markov chain Monte Carlo tree proposals in Bayesian \n");
    MrBayesPrint ("         phylogenetics. Syst. Biol. 57: 86-103.                                  \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Roberts, G. O., and J. S. Rosenthal. 2009. Examples of adaptive MCMC. Jour-\n");
    MrBayesPrint ("         nal of Compuational and Graphical Statistics 18: 349-367.               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   The harmonic mean estimator of model likelihoods, used for Bayes factor tes-  \n");
    MrBayesPrint ("   ting, was discussed by Newton and Raftery (1996). The more accurate stepping- \n");
    MrBayesPrint ("   stone algorithm was first proposed by Xie et al. (2011). The paper by         \n");
    MrBayesPrint ("   Lartillot and Philippe (2006) presents an interesting discussion of the       \n");
    MrBayesPrint ("   shortcomings of the harmonic mean estimator and describes thermodynamic       \n");
    MrBayesPrint ("   integration, a technique that is similar to the stepping-stone algorithm.     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Newton, M. A., and A. E. Raftery. 1994. Approcimate Bayesian inference     \n");
    MrBayesPrint ("         with the weighted likelihood bootstrap. J. R. Stat. Soc. B. 56. 3-48.   \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Lartillot, N., and H. Philippe. 2006. Computing Bayes factors using        \n");
    MrBayesPrint ("         thermodynamic integration. Syst. Biol. 55: 195-207.                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Xie, W., P. O. Lewis, Y. Fan, L. Kuo, and M.-H. Chen. 2011. Improving      \n");
    MrBayesPrint ("         marginal likelihood estimation for Bayesian phylogenetic model          \n");
    MrBayesPrint ("         selection. Syst. Biol. 60: 150-160.                                     \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   For unconstrained branch lengths, MrBayes implements the compound Dirichlet   \n");
    MrBayesPrint ("   priors for branch lengths described by Rannala et al. (2012) and Zhang et al. \n");
    MrBayesPrint ("   (2012). Compared with the i.i.d. exponential and uniform priors for branch    \n");
    MrBayesPrint ("   lengthes in the previous versions of MrBayes, the Dirichlet priors appear more\n");
    MrBayesPrint ("   reasonable and may avoid the problem of extremely long trees, as discussed by \n");
    MrBayesPrint ("   Brown et al. (2010) and Marshall (2010). The two-exponential prior on internal\n");
    MrBayesPrint ("   and external branch lengths described by Yang & Rannala (2005) and Yang (2007)\n");
    MrBayesPrint ("   is also implemented in this version.                                          \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Brown, J. M., S. M. Hedtke, A. R. Lemmon, and E. M. Lemmon. 2010. When     \n");
    MrBayesPrint ("         trees  grow too long: investigating the causes of highly inaccurate     \n");
    MrBayesPrint ("         Bayesian branch-length estimates. Syst. Biol. 59:145-161.               \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Marshall, D. C. 2010. Cryptic failure of partitioned Bayesian phylogenetic \n");
    MrBayesPrint ("         analyses: lost in the land of long trees. Syst. Biol. 59:108-117.       \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Rannala, B., T. Zhu, and Z. Yang. 2012. Tail paradox, partial              \n");
    MrBayesPrint ("         identifiability and influential priors in Bayesian branch length        \n");
    MrBayesPrint ("         inference. Mol. Biol. Evol. 29:325-335.                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Zhang, C., B. Rannala, and Z. Yang. 2012. Robustness of compound Dirichlet \n");
    MrBayesPrint ("         priors for Bayesian inference of branch lengths. Syst. Biol. 61:779-784.\n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z. 2007. Fair-balance paradox, star-tree paradox and Bayesian        \n");
    MrBayesPrint ("         phylogenetics. Mol. Biol. Evol. 24:1639-1655.                           \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      Yang, Z., and B. Rannala. 2005. Branch-length prior influences Bayesian    \n");
    MrBayesPrint ("         posterior probability of phylogeny. Syst. Biol. 54:455-470.             \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

    return (NO_ERROR);
}


int DoConstraint (void)
{
    int         i, howMany;
    int         *tset;

    if (constraintType == PARTIAL)
        tset=tempSetNeg;
    else
        tset=tempSet;

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToGivenSet (fromI, toJ, everyK, 1, tset) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToGivenSet (fromI, toJ, everyK, 1, tset) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToGivenSet (fromI, toJ, everyK, 1, tset) == ERROR)
            return (ERROR);
        }
            
    /* check that this is not a stupid constraint */
    howMany = 0;
    for (i=0; i<numTaxa; i++)
        if (tempSet[i] != 0)
            howMany++;

    if (howMany == 0)
        {
        MrBayesPrint ("%s   This constraint does not include any taxa and will not be defined\n", spacer);
        return (ERROR);
        }

    if (constraintType == HARD)
        {
        if (howMany == numTaxa)
            {
            /* We allow this so we can report states from and calibrate root */
            }
        
        } /*end constraintType == HARD */
    else if (constraintType == PARTIAL)
        {
        if (howMany == 1)
            {
            MrBayesPrint ("%s   This partial constraint includes only one taxon. It is always satisfied and will not be defined.\n", spacer);
            return (ERROR);
            }

        howMany = 0;
        for (i=0; i<numTaxa; i++)
            {
            if (tempSetNeg[i] != 0)
                {
                howMany++;
                if (tempSetNeg[i] == tempSet[i])
                    {
                    MrBayesPrint ("%s   Two sets of taxa in partial constraint are not allowed to intersect. Constraint will not be defined\n", spacer);
                    return (ERROR);
                    }
                }
            }
        if (howMany == 0)
            {
            MrBayesPrint ("%s   This partial constraint does not include any taxa in the second set and will not be defined\n", spacer);
            return (ERROR);
            }
        }
    else if (constraintType == NEGATIVE)
        {
        if (howMany == 1)
            {
            MrBayesPrint ("%s   Negative constraint should include more than one taxon. Constraint will not be defined\n", spacer);
            return (ERROR);
            }
        }

    /* add name to constraintNames */
    if (AddString (&constraintNames, numDefinedConstraints, tempSetName) == ERROR)
        {
        MrBayesPrint ("%s   Problem adding constraint %s to list\n", spacer, tempSetName);
        return (ERROR);
        }

    /* store tempSet */
    AddBitfield (&definedConstraint, numDefinedConstraints, tempSet, numTaxa);
    if (constraintType == PARTIAL)
        {
        AddBitfield (&definedConstraintTwo, numDefinedConstraints, tempSetNeg, numTaxa);
        }
    else
        {
        definedConstraintTwo = (BitsLong **) SafeRealloc ((void *)(definedConstraintTwo), ((size_t)numDefinedConstraints+1)*sizeof(BitsLong *));
        if (definedConstraintTwo==NULL)
            return ERROR;
        definedConstraintTwo[numDefinedConstraints]=NULL;
        }
    
    /* add a default node calibration */
    nodeCalibration = (Calibration *) SafeRealloc ((void *)nodeCalibration, ((size_t)numDefinedConstraints+1)*sizeof(Calibration));
    nodeCalibration[numDefinedConstraints].prior            = defaultCalibration.prior;
    nodeCalibration[numDefinedConstraints].priorParams[0]   = defaultCalibration.priorParams[0];
    nodeCalibration[numDefinedConstraints].priorParams[1]   = defaultCalibration.priorParams[1];
    nodeCalibration[numDefinedConstraints].priorParams[2]   = defaultCalibration.priorParams[2];
    nodeCalibration[numDefinedConstraints].min              = defaultCalibration.min;
    nodeCalibration[numDefinedConstraints].max              = defaultCalibration.max;
    nodeCalibration[numDefinedConstraints].LnPriorProb      = defaultCalibration.LnPriorProb;
    nodeCalibration[numDefinedConstraints].LnPriorRatio     = defaultCalibration.LnPriorRatio;
    strcpy(nodeCalibration[numDefinedConstraints].name, defaultCalibration.name);

    /* increment number of defined constraints */
    numDefinedConstraints++;

    /* reallocate and initialize space for activeConstraints */
    for (i=0; i<numCurrentDivisions; i++)
        {
        modelParams[i].activeConstraints = (int *) SafeRealloc((void *)(modelParams[i].activeConstraints), (size_t)numDefinedConstraints*sizeof(int));
        modelParams[i].activeConstraints[numDefinedConstraints-1] = NO;
        }

    /* reallocate and initialize space for tempActiveConstraints */
    tempActiveConstraints = (int *) SafeRealloc((void *)(tempActiveConstraints), (size_t)numDefinedConstraints*sizeof(int));
    tempActiveConstraints[numDefinedConstraints-1] = NO;

    definedConstraintsType = (enum ConstraintType *) SafeRealloc((void *)(definedConstraintsType), (size_t)numDefinedConstraints*sizeof(enum ConstraintType));
    if (definedConstraintsType==NULL)
        return ERROR;
    definedConstraintsType[numDefinedConstraints-1] = constraintType;

    definedConstraintPruned = (BitsLong **) SafeRealloc ((void *)(definedConstraintPruned), (size_t)numDefinedConstraints*sizeof(BitsLong *));
    if (definedConstraintPruned==NULL)
        return ERROR;
    definedConstraintPruned[numDefinedConstraints-1]=NULL;

    definedConstraintTwoPruned = (BitsLong **) SafeRealloc ((void *)(definedConstraintTwoPruned), (size_t)numDefinedConstraints*sizeof(BitsLong *));
    if (definedConstraintTwoPruned==NULL)
        return ERROR;
    definedConstraintTwoPruned[numDefinedConstraints-1]=NULL;

    /* show taxset (for debugging) */
    // for (i=0; i<numTaxa; i++)
    //     MrBayesPrint ("%4d  %4d\n", i+1, taxaInfo[i].constraints[numDefinedConstraints-1]);

    return (NO_ERROR);
}


int DoConstraintParm (char *parmName, char *tkn)
{
    int         i, index, tempInt;
    MrBFlt      tempD;
    static int  *tempSetCurrent;
    
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before constraints can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check size of constraint name */
            if (strlen(tkn) > 99)
                {
                MrBayesPrint ("%s   Constraint name is too long\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a constraint */
            if (numDefinedConstraints > 0)
                {
                if (CheckString (constraintNames, numDefinedConstraints, tkn, &index) == ERROR)
                    {
                    /* an ERROR returned if the constraint name has not been used. we _want_ to be here */
                    }
                else
                    {
                    MrBayesPrint ("%s   Constraint name '%s' has been used previously\n", spacer, tkn);
                    return (ERROR);
                    }
                }
                
            /* copy the name to the temporary constraint names string */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 0;

            constraintType = HARD; /* set default constrain type */
            tempSetCurrent=tempSet;
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            MrBayesPrint ("%s   Defining constraint called '%s'\n", spacer, tkn);
            foundExp = NO;
            foundFirst = YES;
            foundEqual = NO;
            isNegative = NO;
            foundColon = NO;
            expecting = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(DASH);
            expecting |= Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        foundEqual = YES;
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        isNegative = NO;
        expecting = Expecting(NUMBER);
        expecting |= Expecting(DASH);
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        isNegative = NO;
        foundExp = NO;
        expecting = Expecting(EQUALSIGN);
        }
    else if (expecting == Expecting(DASH))
        {
        if (foundExp == YES)
            isNegative = YES;
        else
            foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (foundFirst == YES && foundEqual == NO)
            {
            /* We are filling in the probability for the constraint. Specifically, we expect exp(number). */
            if (IsSame ("Partial", tkn) == SAME)
                {
                for (i=0; i<numTaxa; i++)
                    tempSetNeg[i] = 0;

                constraintType = PARTIAL;
                expecting = Expecting(EQUALSIGN);
                expecting |= Expecting(ALPHA);
                }
            else if (IsSame ("Hard", tkn) == SAME)
                {
                constraintType = HARD;
                expecting = Expecting(EQUALSIGN);
                expecting |= Expecting(ALPHA);
                }
            else if (IsSame ("Negative", tkn) == SAME)
                {
                constraintType = NEGATIVE;
                expecting = Expecting(EQUALSIGN);
                expecting |= Expecting(ALPHA);
                }
            else if (IsSame ("Exp", tkn) == SAME || IsSame ("Exp", tkn) == CONSISTENT_WITH)
                {
                foundExp = YES;
                foundDash = NO;
                isNegative = NO;
                expecting  = Expecting(LEFTPAR);
                }
            else
                {
                MrBayesPrint ("%s   Do not understand %s\n", spacer, tkn);
                return (ERROR);
                }
            }
        else
            {
            /* We are defining a constraint in terms of a taxon set (called tkn, here) or we are referring to
               the taxon name. We should be able to find tkn in the list of taxon set names or in the list
               of taxon names. If we cannot, then we have a problem and return an error. */
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                if (numTaxaSets < 1)
                    {
                    MrBayesPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                /* add taxa from taxset tkn to new tempSet */
                for (i=0; i<numTaxa; i++)
                    {
                    if (IsBitSet(i, taxaSet[index]) == YES)
                        {
                        tempSetCurrent[i] = 1;
                        }
                    }
                }
            else
                {
                tempSetCurrent[index] = 1;
                }
            fromI = toJ = everyK = -1;

            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            if (constraintType != PARTIAL || foundColon == YES)
                expecting |= Expecting(SEMICOLON);
            else
                expecting |= Expecting(COLON);
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (foundFirst == YES && foundEqual == NO)
            {
            /* We are filling in the probability for the constraint. Specifically, we expect number. */
            sscanf (tkn, "%lf", &tempD);        
            if (foundExp == NO && tempD < 0.0)
                {
                MrBayesPrint ("%s   The probability of a clade cannot be less than zero\n", spacer, tkn);
                return (ERROR);
                }
            if (isNegative == YES || foundDash == YES)
                tempD *= -1.0;
            if (foundExp == YES)
                {
                expecting  = Expecting(RIGHTPAR);
                }
            else
                {
                expecting  = Expecting(EQUALSIGN);
                }
            foundFirst = NO;
            foundDash = NO;
            }
        else
            {       
            if (strlen(tkn) == 1 && !strcmp(tkn, "."))
                {
                tempInt = numTaxa;
                }
            else
                {
                sscanf (tkn, "%d", &tempInt);
                if (tempInt <= 0 || tempInt > numTaxa)
                    {
                    MrBayesPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
                    return (ERROR);
                    }
                }
            tempInt--;
            if (foundDash == YES)
                {
                if (fromI >= 0)
                    toJ = tempInt;
                else
                    {
                    MrBayesPrint ("%s   Improperly formatted constraint\n", spacer);
                    return (ERROR);
                    }
                foundDash = NO;
                }
            else if (foundSlash == YES)
                {
                tempInt++;
                if (tempInt <= 1)
                    {
                    MrBayesPrint ("%s   Improperly formatted constraint\n", spacer);
                    return (ERROR);
                    }
                if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                    everyK = tempInt;
                else
                    {
                    MrBayesPrint ("%s   Improperly formatted constraint\n", spacer);
                    return (ERROR);
                    }
                foundSlash = NO;
                }
            else
                {
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToGivenSet (fromI, toJ, everyK, 1, tempSetCurrent) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    }
                else if (fromI < 0 && toJ < 0)
                    {
                    fromI = tempInt;
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                    {
                    if (AddToGivenSet (fromI, toJ, everyK, 1, tempSetCurrent) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    toJ = everyK = -1;
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                    {
                    if (AddToGivenSet (fromI, toJ, everyK, 1, tempSetCurrent) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    toJ = everyK = -1;
                    }
                else
                    {
                    MrBayesPrint ("%s   Improperly formatted constraint\n", spacer);
                        {
                        return (ERROR);
                        }
                    }
                }

            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(DASH);
            expecting |= Expecting(BACKSLASH);
            if (constraintType != PARTIAL || foundColon == YES)
                expecting |= Expecting(SEMICOLON);
            else
                expecting |= Expecting(COLON);
            }
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(COLON))
        {
        if (foundColon == YES)
            {
            MrBayesPrint ("%s   Improperly formatted constraint: two colon charactors in constraint command.\n", spacer);
            return (ERROR);
            }

        /* add set to tempSet */
        if (fromI >= 0 && toJ < 0)
            {
            if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                return (ERROR);
            }
        else if (fromI >= 0 && toJ >= 0 && everyK < 0)
            {
            if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                return (ERROR);
            }
        else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
            {
            if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                return (ERROR);
            }
        fromI = toJ = everyK = -1;
        foundDash = foundSlash = NO;

        foundColon = YES;
        tempSetCurrent = tempSetNeg;
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoCtype (void)
{
    int         i, foundIllegal, marks[5], numAppliedTo;

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with ctype */
    numAppliedTo = 0;
    for (i=0; i<5; i++)
        marks[i] = NO;
    for (i=0; i<numChar; i++)
        {
        if (tempSet[i] != 0)
            {
            foundIllegal = NO;
            if (charOrdering != UNORD)
                {
                if (charInfo[i].charType == DNA)
                    {
                    foundIllegal = YES;
                    if (marks[0] == NO)
                        MrBayesPrint ("%s   Ctype not applied to DNA states which must be unordered\n", spacer);
                    marks[0] = YES;
                    }
                else if (charInfo[i].charType == RNA)
                    {
                    foundIllegal = YES;
                    if (marks[1] == NO)
                        MrBayesPrint ("%s   Ctype not applied to RNA states which must be unordered\n", spacer);
                    marks[1] = YES;
                    }
                else if (charInfo[i].charType == PROTEIN)
                    {
                    foundIllegal = YES;
                    if (marks[2] == NO)
                        MrBayesPrint ("%s   Ctype not applied to amino acid states which must be unordered\n", spacer);
                    marks[2] = YES;
                    }
                else if (charInfo[i].charType == RESTRICTION)
                    {
                    foundIllegal = YES;
                    if (marks[3] == NO)
                        MrBayesPrint ("%s   Ctype not applied to restriction site states which must be unordered\n", spacer);
                    marks[3] = YES;
                    }
                else if (charInfo[i].charType == CONTINUOUS)
                    {
                    foundIllegal = YES;
                    if (marks[4] == NO)
                        MrBayesPrint ("%s   Ctype not applied to continuous characters\n", spacer);
                    marks[4] = YES;
                    }
                }
            if (foundIllegal == NO)
                {
                charInfo[i].ctype = charOrdering;
                numAppliedTo++;
                }
            }
        }
    if (numAppliedTo > 0)
        {
        MrBayesPrint ("%s   Ctype was applied to %d standard characters\n", spacer, numAppliedTo);
        }
    else
        {
        MrBayesPrint ("%s   No standard characters found to apply ctype to\n", spacer);
        }
    
#   if 0
    for (i=0; i<numChar; i++)
        MrBayesPrint ("%4d -- %d\n", i, ctype[i]);
#   endif

    return (NO_ERROR);
}


int DoCtypeParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
    
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before typesets can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            if (IsSame ("Ordered", tkn) == SAME || IsSame ("Ordered", tkn) == CONSISTENT_WITH)
                charOrdering = ORD;
            else if (IsSame ("Unordered", tkn) == SAME || IsSame ("Unordered", tkn) == CONSISTENT_WITH)
                charOrdering = UNORD;
            else if (IsSame ("Dollo", tkn) == SAME || IsSame ("Dollo", tkn) == CONSISTENT_WITH)
                charOrdering = DOLLO;
            else if (IsSame ("Irreversible", tkn) == SAME || IsSame ("Irreversible", tkn) == CONSISTENT_WITH)
                charOrdering = IRREV;
            else
                {
                MrBayesPrint ("%s   Do not understand delimiter \"%s\"\n", spacer, tkn);
                return (ERROR);
                }
            
            /* clear tempSet */
            for (i=0; i<numChar; i++)
                tempSet[i] = 0;
            
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            MrBayesPrint ("%s   Setting characters to %s\n", spacer, tkn);
            expecting = Expecting(COLON);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting  = Expecting(ALPHA) | Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        /* first, check that we are not trying to put in another character ordering */
        if (IsSame ("Ordered", tkn) == SAME || IsSame ("Ordered", tkn) == CONSISTENT_WITH)
            {
            MrBayesPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
            return (ERROR);
            }
        else if (IsSame ("Unordered", tkn) == SAME || IsSame ("Unordered", tkn) == CONSISTENT_WITH)
            {
            MrBayesPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
            return (ERROR);
            }
        else if (IsSame ("Dollo", tkn) == SAME || IsSame ("Dollo", tkn) == CONSISTENT_WITH)
            {
            MrBayesPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
            return (ERROR);
            }
        else if (IsSame ("Irreversible", tkn) == SAME || IsSame ("Irreversible", tkn) == CONSISTENT_WITH)
            {
            MrBayesPrint ("%s   You cannot specify more than one ordering with a single use of ctype\n", spacer, tkn);
            return (ERROR);
            }
        
        /* We are defining a type set in terms of another (called tkn, here). We should be able
           to find tkn in the list of character set names. If we cannot, then we have a problem and
           return an error. */
        if (IsSame ("All", tkn) == SAME)
            {
            for (i=0; i<numChar; i++)
                tempSet[i] = 1;
            fromI = toJ = everyK = -1;
            expecting = Expecting(SEMICOLON);
            }
        else
            {
            if (numCharSets < 1)
                {
                MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
                {
                MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
                
            /* add characters from charset tkn to new tempset */
            for (i=0; i<numChar; i++)
                {
                if (IsBitSet(i, charSet[index]) == YES)
                    tempSet[i] = 1;
                }
            fromI = toJ = everyK = -1;
            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numChar;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted ctype\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                MrBayesPrint ("%s   Improperly formatted ctype\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted ctype\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                MrBayesPrint ("%s   Improperly formatted ctype\n", spacer);
                    {
                    return (ERROR);
                    }
                }
                
            }

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoDelete (void)
{
    int         i, alreadyDone;

    MrBayesPrint ("%s   Excluding taxa\n", spacer);

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with taxaset */
    alreadyDone = NO;
    for (i=0; i<numTaxa; i++)
        {
        if (tempSet[i] == 1)
            {
            if (taxaInfo[i].isDeleted == YES && alreadyDone == NO)
                {
                MrBayesPrint ("%s   Some taxa already excluded\n", spacer);
                alreadyDone = YES;
                }
            taxaInfo[i].isDeleted = YES;
            }
        }

    SetLocalTaxa ();
    if (SetUpAnalysis(&globalSeed) == ERROR)
        return ERROR;

    /* show tempSet (for debugging) */
#   if 0
    for (i=0; i<numTaxa; i++)
        MrBayesPrint ("%4d  %4d\n", i+1, tempSet[i]);
#   endif

    return (NO_ERROR);
}


int DoDeleteParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
        
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before you can delete taxa\n", spacer);
        return (ERROR);
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        fromI = toJ = everyK = -1;
        foundDash = NO;
        for (i=0; i<numTaxa; i++) /* clear tempSet */
            tempSet[i] = 0;
        foundFirst = YES;
        }

    if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 1;
            }
        else
            {
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                /* we are using a pre-defined taxa set */
                if (numTaxaSets < 1)
                    {
                    MrBayesPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                /* add taxa from taxset tkn to new tempSet */
                for (i=0; i<numTaxa; i++)
                    {
                    if (IsBitSet (i, taxaSet[index]) == YES)
                        tempSet[i] = 1;
                    }
                }
            else
                {
                /* we found the taxon name */
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                        return (ERROR);
                    }
                    
                tempSet[index] = 1;
                }
            }
        foundDash = NO;
        fromI = toJ = everyK = -1;
        
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && !strcmp(tkn, "."))
            tempInt = numTaxa;
        else
            {
            sscanf (tkn, "%d", &tempInt);
            if (tempInt <= 0 || tempInt > numTaxa)
                {
                MrBayesPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
                return (ERROR);
                }
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted delete set\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                MrBayesPrint ("%s   Improperly formatted delete set\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoDimensions (void)
{
    if (inDataBlock == NO && inTaxaBlock == NO && inCharactersBlock == NO)
        {
        MrBayesPrint ("%s   Dimensions can only be defined in a data, characters or taxa block\n", spacer);
        return (ERROR);
        }

    /* other problems are detected already when reading in DoDimensionsParm */
    if (inDataBlock == YES && (defTaxa == NO || defChars == NO))
        {
        MrBayesPrint ("%s   Expecting both Ntax and Nchar to be defined in a data block\n", spacer);
        return (ERROR);
        }

    /* allocate matrix */
    if (inTaxaBlock == YES)
        {
        if (AllocTaxa () == ERROR)
            return ERROR;
        MrBayesPrint ("%s   Defining new set of %d taxa\n", spacer, numTaxa);
        }

    if (inCharactersBlock == YES)
        {
        if (AllocMatrix() == ERROR)
            return (ERROR);
        MrBayesPrint ("%s   Defining new character matrix with %d characters\n", spacer, numChar);
        }

    if (inDataBlock == YES)
        {
        if (AllocMatrix() == ERROR)
            return (ERROR);
        MrBayesPrint ("%s   Defining new matrix with %d taxa and %d characters\n", spacer, numTaxa, numChar);
        }

    return (NO_ERROR);
}


int DoDimensionsParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Ntax (numTaxa) *****************************************************************/
        if (!strcmp(parmName, "Ntax"))
            {
            if (inCharactersBlock == YES)
                {
                MrBayesPrint ("%s   You cannot define ntax in a characters block\n");
                return (ERROR);
                }
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &numTaxa);
                defTaxa = YES;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Nchar (numChar) ****************************************************************/
        else if (!strcmp(parmName, "Nchar"))
            {
            if (inTaxaBlock == YES)
                {
                MrBayesPrint ("%s   You cannot define nchar in a taxa block\n");
                return (ERROR);
                }
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &numChar);
                defChars = YES;
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


int DoDisclaimer (void)
{
    MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
    MrBayesPrint ("   Disclaimer                                                                    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   Copyright 2003 by John P. Huelsenbeck and Fredrik Ronquist                    \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   This software package is provided \"as is\" and without a warranty of any     \n");
    MrBayesPrint ("   kind. In no event shall the authors be held responsible for any damage        \n");
    MrBayesPrint ("   resulting from the use of this software. The program--including source code,  \n");
    MrBayesPrint ("   example data sets, and executables--is distributed free of charge for         \n");
    MrBayesPrint ("   academic use only.                                                            \n");
    MrBayesPrint ("   ---------------------------------------------------------------------------   \n");

    return (NO_ERROR);
}


int DoEndBlock (void)
{
    if (inMrbayesBlock == YES)
        {
        MrBayesPrint ("   Exiting mrbayes block\n");
        inMrbayesBlock = NO;
        }
    else if (inDataBlock == YES)
        {
        MrBayesPrint ("   Exiting data block\n");
        inDataBlock = NO;
        }
    else if (inCharactersBlock == YES)
        {
        MrBayesPrint ("   Exiting characters block\n");
        inCharactersBlock = NO;
        }
    else if (inTaxaBlock == YES)
        {
        MrBayesPrint ("   Exiting taxa block\n");
        if (numNamedTaxa < numTaxa)
            {
            MrBayesPrint ("%s   Leaving taxa block without taxon labels being defined\n", spacer);
            FreeTaxa();
            }
        inTaxaBlock = NO;
        }
    else if (inTreesBlock == YES)
        {
        MrBayesPrint ("   Exiting trees block\n");
        inTreesBlock = NO;
        ResetTranslateTable();
        }
    else if (inForeignBlock == YES)
        {
        MrBayesPrint ("   Exiting foreign block\n");
        inForeignBlock = NO;
        }
    else
        {
        MrBayesPrint ("   Unknown \"end\" statement\n");
        return (ERROR);
        }

    strcpy(spacer,"");  /* reset indentation */
    return (NO_ERROR);
}


int DoExecute (void)
{
    int         rc, cmdLine, lineTerm, longestLineLength, nErrors;
    char        *s, exeFileName[100];
    FILE        *fp;
    CmdType     *oldCommandPtr;
    char        *oldTokenP, oldToken[CMD_STRING_LENGTH];
#   if defined (MPI_ENABLED)
    int         sumErrors;
#   endif
        
    nErrors = 0;
    cmdLine = 0;
    numOpenExeFiles++;
    s = NULL;
    strncpy (exeFileName, inputFileName, 98);
    
    if (numOpenExeFiles > 1)
        MrBayesPrint ("\n%s   Executing file \"%s\"...\n\n", spacer, inputFileName);
    else
        MrBayesPrint ("%s   Executing file \"%s\"\n", spacer, inputFileName);

    /* Save old command ptr, token pointer and token */
    oldCommandPtr = commandPtr;
    oldTokenP     = tokenP;
    strcpy(oldToken, token);

    /* open binary file */
    if ((fp = OpenBinaryFileR(inputFileName)) == NULL)
        nErrors++;

    /* set indentation to 0 */
    strcpy (spacer, "");
    
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif

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
        nErrors++;
        }
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif
            
    /* find length of longest line */
    longestLineLength = LongestLine (fp);
    MrBayesPrint ("%s   Longest line length = %d\n", spacer, longestLineLength);
    longestLineLength += 10;
    
    /* allocate a string long enough to hold a line */
    s = (char *)SafeMalloc((size_t)longestLineLength * sizeof(char));
    if (!s)
        {
        MrBayesPrint ("%s   Problem allocating string for reading file\n", spacer);
        nErrors++;
        }
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif

    /* close binary file */
    SafeFclose (&fp);
    
    /* open text file */
    if ((fp = OpenTextFileR(inputFileName)) == NULL)
        {
        nErrors++;
        }
#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif
    
    /* parse file, reading each line in turn */
    MrBayesPrint ("%s   Parsing file\n", spacer);

    inMrbayesBlock = inDataBlock = inForeignBlock = inTreesBlock = NO;
    foundNewLine = NO;
    expecting = Expecting(COMMAND);
    cmdLine = 0;

    /* read lines into s until end of file */
    while ( fgets(s, longestLineLength, fp) != NULL)
    {
        foundNewLine = YES;
        cmdLine++;

        /* process string if not empty */
        if (strlen(s) > 1)
            {
            /* check that all characters in the string are valid */
            if (CheckStringValidity (s) == ERROR)
                {
                nErrors++;
                }
#           if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
                goto errorExit;
                }
#           else
            if (nErrors > 0)
                goto errorExit;
#           endif
                
            /* interpret commands on line */
            rc = ParseCommand (s);
            if (rc == ERROR)
                nErrors++;
#           if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                {
                MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
                goto errorExit;
                }
#           else
            if (nErrors > 0)
                goto errorExit;
#           endif
            if (rc == NO_ERROR_QUIT)
                nErrors++;
#           if defined (MPI_ENABLED)
            MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            if (sumErrors > 0)
                goto quitExit;
#           else
            if (nErrors > 0)
                goto quitExit;
#           endif
            }
        }
    
    MrBayesPrint ("%s   Reached end of file\n", spacer);

    if (inComment == YES)
        nErrors++;

#   if defined (MPI_ENABLED)
    MPI_Allreduce (&nErrors, &sumErrors, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (sumErrors > 0)
        {
        MrBayesPrint ("%s   There was an error on at least one processor\n", spacer);
        goto errorExit;
        }
#   else
    if (nErrors > 0)
        goto errorExit;
#   endif

    if (s)
        free (s);
    SafeFclose (&fp);
    numOpenExeFiles--;

    if (numOpenExeFiles > 0)
        {
        inMrbayesBlock = YES;
        MrBayesPrint ("\n   Returning execution to calling file ...\n\n");
        strcpy (spacer, "   ");
        }
    else
        strcpy (spacer, "");

    commandPtr = oldCommandPtr;

    return (NO_ERROR);
    
    quitExit:
        if (s)
            free (s);
        SafeFclose (&fp);
        numOpenExeFiles--;
        if (numOpenExeFiles > 0)
            {
            inMrbayesBlock = YES;
            strcpy (spacer, "   ");
            }
        else
            strcpy (spacer, "");

        commandPtr = oldCommandPtr;
        tokenP     = oldTokenP;
        strcpy(token, oldToken);

        return (NO_ERROR_QUIT);
            
    errorExit:
        if (inComment == YES)
            {
            MrBayesPrint ("%s   ERROR: Reached end of file while in comment.\n", spacer);
            inComment = NO;
            numComments = 0;
            }
        if (fp)
            {
            MrBayesPrint ("%s   The error occurred when reading char. %d-%d on line %d\n", spacer,
                (size_t)(tokenP-s)-strlen(token)+1, (size_t)(tokenP-s), cmdLine);
            MrBayesPrint ("%s      in the file '%s'\n", spacer, exeFileName);
            }
        if (s)
            free (s);
        SafeFclose (&fp);

        /* make sure we exit the block we were reading from correctly */
        if (inMrbayesBlock == YES)
            inMrbayesBlock = NO;
        else if (inDataBlock == YES)
            inDataBlock = NO;
        else if (inTreesBlock == YES)
            {
            inTreesBlock = NO;
            ResetTranslateTable();
            }
        else if (inForeignBlock == YES)
            inForeignBlock = NO;

        /* make sure correct return if we came from mrbayes block in another execute file */
        if (numOpenExeFiles > 1)
            {
            inMrbayesBlock = YES;
            MrBayesPrint ("\n   Returning execution to calling file ...\n\n");
            strcpy (spacer, "   ");
            }
        else
            {
            strcpy (spacer, "");
            MrBayesPrint ("\n   Returning execution to command line ...\n\n");
            }

        numOpenExeFiles--;  /* we increase the value above even if no file is successfully opened */

        /* restore state of globals */
        commandPtr = oldCommandPtr;
        tokenP = oldTokenP;
        strcpy(token, oldToken);

        return (ERROR);
}


int DoExecuteParm (char *parmName, char *tkn)
{
    if (strlen(tkn)>99)
        {
        MrBayesPrint ("%s   Maximum allowed length of file name is 99 characters. The given name:\n", spacer);
        MrBayesPrint ("%s      '%s'\n", spacer,tkn);
        MrBayesPrint ("%s   has %d characters.\n", spacer,strlen(tkn));
        return (ERROR);
        }
    strcpy (inputFileName, tkn);
    
    expecting = Expecting (SEMICOLON);

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoExclude (void)
{
    int         i, alreadyDone;

    MrBayesPrint ("%s   Excluding character(s)\n", spacer);

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with charset */
    alreadyDone = NO;
    for (i=0; i<numChar; i++)
        {
        if (tempSet[i] == 1)
            {
            if (charInfo[i].isExcluded == YES && alreadyDone == NO)
                {
                MrBayesPrint ("%s   Some characters already excluded\n", spacer);
                alreadyDone = YES;
                }
            charInfo[i].isExcluded = YES;
            }
        }
        
    foundFirst = NO;

    /* reset analysis to recompress data */
    if (SetUpAnalysis(&globalSeed) == ERROR)
        return ERROR;

    return (NO_ERROR);
}


int DoExcludeParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
        
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before you can exclude characters\n", spacer);
        return (ERROR);
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        fromI = toJ = everyK = -1;
        foundDash = foundSlash = NO;
        for (i=0; i<numChar; i++) /* clear tempSet */
            tempSet[i] = 0;
        foundFirst = YES;
        }

    if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numChar; i++)
                tempSet[i] = 1;
            }
        else if (IsSame ("Missambig", tkn) == SAME || IsSame ("Missambig", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numChar; i++)
                {
                if (charInfo[i].isMissAmbig == YES)
                    tempSet[i] = 1;
                }
            }
        else
            {
            /* we are using a pre-defined character set */
            if (numCharSets < 1)
                {
                MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
                {
                MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            /* add characters from charset tkn to new tempSet */
            for (i=0; i<numChar; i++)
                {
                if (IsBitSet(i, charSet[index]) == YES)
                    tempSet[i] = 1;
                }
            fromI = toJ = everyK = -1;
            }

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numChar;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted exclude set\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                MrBayesPrint ("%s   Improperly formatted exclude set\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted exclude set\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                MrBayesPrint ("%s   Improperly formatted exclude set\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoFormat (void)
{
    if (inDataBlock == NO && inCharactersBlock == NO)
        {
        MrBayesPrint ("%s   Formats can only be defined in a data or characters block\n", spacer);
        return (ERROR);
        }

    return CheckInitialPartitions();
}


int DoFormatParm (char *parmName, char *tkn)
{
    int         i, tempInt;
    char        tempStr[100];
    
    if (inDataBlock == NO && inCharactersBlock == NO)
        {
        MrBayesPrint ("%s   Formats can only be defined in a data or characters block\n", spacer);
        return (ERROR);
        }
    if (defTaxa == NO || defChars == NO)
        {
        MrBayesPrint ("%s   The dimensions of the matrix must be defined before the format\n", spacer);
        return (ERROR);
        }
    
    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        if (!strcmp(parmName, "Interleave"))
            {
            expecting = Expecting(EQUALSIGN) | Expecting(PARAMETER) | Expecting(SEMICOLON);
            isInterleaved = YES;
            }
        }
    else
        {
        /* set Datatype (dataType) ************************************************************/
        if (!strcmp(parmName, "Datatype"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (isMixed == NO)
                        {
                        if (!strcmp(tempStr, "Dna"))
                            dataType = DNA;
                        else if (!strcmp(tempStr, "Rna"))
                            dataType = RNA;
                        else if (!strcmp(tempStr, "Protein"))
                            dataType = PROTEIN;
                        else if (!strcmp(tempStr, "Restriction"))
                            dataType = RESTRICTION;
                        else if (!strcmp(tempStr, "Standard"))
                            dataType = STANDARD;
                        else if (!strcmp(tempStr, "Continuous"))
                            dataType = CONTINUOUS;
                        else if (!strcmp(tempStr, "Mixed"))
                            {
                            dataType = MIXED;
                            isMixed = YES;
                            for (i=0; i<numChar; i++)
                                tempSet[i] = 0;
                            fromI = toJ = everyK = -1;
                            foundDash = foundSlash = NO;
                            numDivisions = 0;
                            MrBayesPrint ("%s   Data is Mixed\n", spacer);
                            }
                        if (dataType == MIXED)
                            expecting = Expecting(LEFTPAR);
                        else
                            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                        }
                    else
                        {
                        if (!strcmp(tempStr, "Dna"))
                            dataType = DNA;
                        else if (!strcmp(tempStr, "Rna"))
                            dataType = RNA;
                        else if (!strcmp(tempStr, "Protein"))
                            dataType = PROTEIN;
                        else if (!strcmp(tempStr, "Restriction"))
                            dataType = RESTRICTION;
                        else if (!strcmp(tempStr, "Standard"))
                            dataType = STANDARD;
                        else if (!strcmp(tempStr, "Continuous"))
                            dataType = CONTINUOUS;
                        else if (!strcmp(tempStr, "Mixed"))
                            {
                            MrBayesPrint ("%s   Cannot have mixed datatype within a mixed datatype\n", spacer);
                            return (ERROR);
                            }
                        expecting = Expecting(COLON);
                        for (i=0; i<numChar; i++)
                            tempSet[i] = 0;
                        fromI = toJ = everyK = -1;
                        foundDash = foundSlash = NO;
                        }
                    if (isMixed == NO)
                        {
                        numDivisions = 1;
                        for (i=0; i<numChar; i++)
                            {
                            charInfo[i].charType = dataType;
                            partitionId[i][0] = numDivisions;
                            }
                        }
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid data type argument\n", spacer);
                    return (ERROR);
                    }
                if (isMixed == NO)
                    MrBayesPrint ("%s   Data is %s\n", spacer, tempStr);
                else if (strcmp(tempStr, "Mixed"))
                    MrBayesPrint ("%s      Data for partition %d is %s\n", spacer, numDivisions+1, tempStr);
                }
            else if (expecting == Expecting(LEFTPAR))
                {
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(COLON))
                {
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(NUMBER))
                {
                if (strlen(tkn) == 1 && tkn[0] == '.')
                    tempInt = numChar;
                else
                    sscanf (tkn, "%d", &tempInt);
                if (tempInt <= 0 || tempInt > numChar)
                    {
                    MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
                    return (ERROR);
                    }
                tempInt--;
                if (foundDash == YES)
                    {
                    if (fromI >= 0)
                        toJ = tempInt;
                    else
                        {
                        MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
                        return (ERROR);
                        }
                    foundDash = NO;
                    }
                else if (foundSlash == YES)
                    {
                    tempInt++;
                    if (tempInt <= 1)
                        {
                        MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
                        return (ERROR);
                        }
                    if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                        everyK = tempInt;
                    else
                        {
                        MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
                        return (ERROR);
                        }
                    foundSlash = NO;
                    }
                else
                    {
                    if (fromI >= 0 && toJ < 0)
                        {
                        if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                            return (ERROR);
                        fromI = tempInt;
                        }
                    else if (fromI < 0 && toJ < 0)
                        {
                        fromI = tempInt;
                        }
                    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                        {
                        if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                            return (ERROR);
                        fromI = tempInt;
                        toJ = everyK = -1;
                        }
                    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                        {
                        if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                            return (ERROR);
                        fromI = tempInt;
                        toJ = everyK = -1;
                        }
                    else
                        {
                        MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
                            {
                            return (ERROR);
                            }
                        }
                        
                    }
                expecting  = Expecting(NUMBER);
                expecting |= Expecting(DASH);
                expecting |= Expecting(BACKSLASH);
                expecting |= Expecting(COMMA);
                expecting |= Expecting(RIGHTPAR);
                }
            else if (expecting == Expecting(DASH))
                {
                foundDash = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(BACKSLASH))
                {
                foundSlash = YES;
                expecting = Expecting(NUMBER);
                }
            else if (expecting == Expecting(COMMA))
                {
                /* add set to tempSet */
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                for (i=0; i<numChar; i++)
                    {
                    if (tempSet[i] == numDivisions)
                        charInfo[i].charType = dataType;
                    }

                /* merge tempSet */
                for (i=0; i<numChar; i++)
                    {
                    if (tempSet[i] != 0)
                        {
                        if (partitionId[i][0] == 0)
                            {
                            charInfo[i].charType = dataType;
                            partitionId[i][0] = numDivisions + 1;
                            }
                        else
                            {
                            MrBayesPrint ("%s   Improperly formatted partition (same character found in multiple partitions)\n", spacer);
                            return (ERROR);
                            }
                        }
                    }
                
                /* increment number of partitions */
                numDivisions++;             
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(RIGHTPAR))
                {
                /* add set to tempSet */
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, numDivisions+1) == ERROR)
                        return (ERROR);
                    }
                    
                /* merge tempSet */
                for (i=0; i<numChar; i++)
                    {
                    if (tempSet[i] != 0)
                        {
                        if (partitionId[i][0] == 0)
                            {
                            charInfo[i].charType = dataType;
                            partitionId[i][0] = numDivisions + 1;
                            }
                        else
                            {
                            MrBayesPrint ("%s   Improperly formatted partition (same character found in multiple partitions)\n", spacer);
                            return (ERROR);
                            }
                        }
                    }
                
                /* increment number of partitions */
                numDivisions++;             
                if (isMixed == YES)
                    dataType = MIXED;
                    
                if (numDivisions > 1)
                    MrBayesPrint ("%s   There are a total of %d default data divisions\n", spacer, numDivisions);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Interleave (isInterleaved) *****************************************************/
        else if (!strcmp(parmName, "Interleave"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        isInterleaved = YES;
                    else
                        isInterleaved = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for interleaved data\n", spacer);
                    return (ERROR);
                    }
                if (isInterleaved == YES)
                    MrBayesPrint ("%s   Data matrix is interleaved\n", spacer);
                else
                    MrBayesPrint ("%s   Data matrix is not interleaved\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Gap (gapId) ********************************************************************/
        else if (!strcmp(parmName, "Gap"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting  = Expecting(ALPHA);
                expecting |= Expecting(QUESTIONMARK);
                expecting |= Expecting(DASH);
                expecting |= Expecting(NUMBER);
                expecting |= Expecting(ASTERISK);
                expecting |= Expecting(EXCLAMATIONMARK);
                expecting |= Expecting(PERCENT);
                expecting |= Expecting(WEIRD);
                expecting |= Expecting(VERTICALBAR);
                }
            else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
                     ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
                     ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
                     ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
                     ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
                     ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
                     ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
                     ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)) ||
                     ((expecting & Expecting(VERTICALBAR)) == Expecting(VERTICALBAR)))
                {
                if (strlen(tkn) == 1)
                    {
                    if (tkn[0] == matchId || tkn[0] == missingId)
                        {
                        MrBayesPrint ("%s   Gap character matches matching or missing characters\n", spacer);
                        return (ERROR);
                        }
                    gapId = tkn[0];
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid gap argument %s\n", spacer, tkn);
                    return (ERROR);
                    }
                MrBayesPrint ("%s   Gaps coded as %s\n", spacer, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Missing (missingId) ************************************************************/
        else if (!strcmp(parmName, "Missing"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting  = Expecting(ALPHA);
                expecting |= Expecting(QUESTIONMARK);
                expecting |= Expecting(DASH);
                expecting |= Expecting(NUMBER);
                expecting |= Expecting(ASTERISK);
                expecting |= Expecting(EXCLAMATIONMARK);
                expecting |= Expecting(PERCENT);
                expecting |= Expecting(WEIRD);
                expecting |= Expecting(VERTICALBAR);
                }
            else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
                     ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
                     ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
                     ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
                     ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
                     ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
                     ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
                     ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)) ||
                     ((expecting & Expecting(VERTICALBAR)) == Expecting(VERTICALBAR)))
                {
                if (strlen(tkn) == 1)
                    {
                    if (tkn[0] == gapId || tkn[0] == matchId)
                        {
                        MrBayesPrint ("%s   Missing character matches matching or gap characters\n", spacer);
                        return (ERROR);
                        }
                    missingId = tkn[0];
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid missing argument %s\n", spacer, tkn);
                    return (ERROR);
                    }
                MrBayesPrint ("%s   Missing data coded as %s\n", spacer, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Matchchar (matchId) ************************************************************/
        else if (!strcmp(parmName, "Matchchar"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting  = Expecting(ALPHA);
                expecting |= Expecting(QUESTIONMARK);
                expecting |= Expecting(DASH);
                expecting |= Expecting(NUMBER);
                expecting |= Expecting(ASTERISK);
                expecting |= Expecting(EXCLAMATIONMARK);
                expecting |= Expecting(PERCENT);
                expecting |= Expecting(WEIRD);
                expecting |= Expecting(VERTICALBAR);
                }
            else if (((expecting & Expecting(ALPHA)) == Expecting(ALPHA)) || 
                     ((expecting & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK)) || 
                     ((expecting & Expecting(DASH)) == Expecting(DASH)) || 
                     ((expecting & Expecting(NUMBER)) == Expecting(NUMBER)) || 
                     ((expecting & Expecting(ASTERISK)) == Expecting(ASTERISK)) || 
                     ((expecting & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK)) || 
                     ((expecting & Expecting(PERCENT)) == Expecting(PERCENT)) || 
                     ((expecting & Expecting(WEIRD)) == Expecting(WEIRD)) ||
                     ((expecting & Expecting(VERTICALBAR)) == Expecting(VERTICALBAR)))
                {
                if (strlen(tkn) == 1)
                    {
                    if (tkn[0] == gapId || tkn[0] == missingId)
                        {
                        MrBayesPrint ("%s   Matching character matches gap or missing characters\n", spacer);
                        return (ERROR);
                        }
                    matchId = tkn[0];
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid matchchar argument %s\n", spacer, tkn);
                    return (ERROR);
                    }
                MrBayesPrint ("%s   Matching characters coded as %s\n", spacer, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* skip Symbols ***************************************************************/
        else if (!strcmp(parmName, "Symbols"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                MrBayesPrint ("%s   WARNING: MrBayes does not support 'symbols' specification; default symbols assumed\n", spacer);
                readWord=YES;
                expecting = Expecting(ALPHA);
                }
            else if (expecting == Expecting(ALPHA))
                {
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* on Equate return ERROR ***************************************************************/
        else if (!strcmp(parmName, "Equate"))
            {
            MrBayesPrint ("%s   ERROR: MrBayes does not support 'Equate' macros; please remove or comment out\n", spacer);
            return (ERROR);
            }
        else
            return (ERROR);
        }

    return (NO_ERROR);
}


int DoHelp (void)
{
    int         i, j, longestDescription;
    CmdType     *p;

    if (foundFirst == NO)
        {
        longestDescription = 0;
        for (i=1; i<NUMCOMMANDS; i++)
            {
            p = commands + i;
            if ((int)strlen(p->string) > longestDescription)
                longestDescription = (int) strlen(p->string);
            }
        
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Commands that are available from the command                                  \n");
        MrBayesPrint ("   line or from a MrBayes block include:                                         \n");
        MrBayesPrint ("                                                                                 \n");
        for (i=1; i<NUMCOMMANDS; i++)
            {
            p = commands + i;
            if (p->cmdUse == IN_CMD && p->hiding == SHOW)
                {
                MrBayesPrint ("   %s", p->string);
                for (j=0; j<longestDescription - (int) strlen(p->string); j++)
                MrBayesPrint (" ");
                MrBayesPrint (" -- %s\n", p->cmdDescription);
                }
            }
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Commands that should be in a NEXUS file (data                                 \n");
        MrBayesPrint ("   block, trees block or taxa block) include:                                                \n");
        MrBayesPrint ("                                                                                 \n");
        for (i=1; i<NUMCOMMANDS; i++)
            {
            p = commands + i;
            if (p->cmdUse == IN_FILE && p->hiding == SHOW)
                {
                MrBayesPrint ("   %s", p->string);
                for (j=0; j<longestDescription - (int) strlen(p->string); j++)
                MrBayesPrint (" ");
                MrBayesPrint (" -- %s\n", p->cmdDescription);
                }
            }
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that this program supports the use of the shortest unambiguous           \n"); 
        MrBayesPrint ("   spelling of the above commands (e.g., \"exe\" instead of \"execute\").        \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    foundFirst = NO;

    return (NO_ERROR);
}


int DoHelpParm (char *parmName, char *tkn)
{
    int         i, j, tkLen, targetLen, numDiff, numMatches;
    CmdType     *p, *q=NULL;

    if (expecting == Expecting(ALPHA))
        {
        p = commands + 0;
        tkLen = (int) strlen(tkn);
        numMatches = 0;
        for (i=0; i<NUMCOMMANDS; i++)
            {
            targetLen = (int) strlen(p->string);
            if (tkLen <= targetLen)
                {
                for (j=0, numDiff=0; j<tkLen; j++)
                    {
                    if (ChangeCase(tkn[j]) != ChangeCase(p->string[j]))
                        numDiff++;
                    }
                if (numDiff == 0)
                    {
                    numMatches++;
                    q = p;
                    if (tkLen == targetLen)
                        break;
                    }
                }       
            p++;
            }
        if (numMatches == 0)
            {
            MrBayesPrint ("%s   Could not find command \"%s\"\n", spacer, tkn);
            return (ERROR);
            }
        else if (numMatches == 1)
            {
            if (GetUserHelp (q->string) == ERROR)
                {
                MrBayesPrint ("%s   Problem getting help for command \"%s\"\n", spacer, q->string);
                }
            }
        else 
            {
            MrBayesPrint ("%s   Ambiguous command \"%s\"\n", spacer, tkn);
            return (ERROR);
            }
            
        expecting = Expecting(SEMICOLON);
        foundFirst = YES;
        }
    else
        return (ERROR);

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoInclude (void)
{
    int         i, alreadyDone;

    MrBayesPrint ("%s   Including character(s)\n", spacer);

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with excludedChars */
    alreadyDone = NO;
    for (i=0; i<numChar; i++)
        {
        if (tempSet[i] == 1)
            {
            if (charInfo[i].isExcluded == NO && alreadyDone == NO)  
                {
                MrBayesPrint ("%s   Some characters already included\n", spacer);
                alreadyDone = YES;
                }
            charInfo[i].isExcluded = NO;
            }
        }

    /* reset analysis to recompress data */
    if (SetUpAnalysis(&globalSeed) == ERROR)
        return ERROR;

    return (NO_ERROR);
}


int DoIncludeParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
        
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before you can include characters\n", spacer);
        return (ERROR);
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        fromI = toJ = everyK = -1;
        foundDash = foundSlash = NO;
        for (i=0; i<numChar; i++) /* clear tempSet */
            tempSet[i] = 0;
        foundFirst = YES;
        }

    if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numChar; i++)
                tempSet[i] = 1;
            }
        else if (IsSame ("Missambig", tkn) == SAME || IsSame ("Missambig", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numChar; i++)
                {
                if (charInfo[i].isMissAmbig == YES)
                    tempSet[i] = 1;
                }
            }
        else
            {
            /* we are using a pre-defined character set */
            if (numCharSets < 1)
                {
                MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
                {
                MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            /* add characters from charset tkn to new tempSet */
            for (i=0; i<numChar; i++)
                {
                if (IsBitSet(i, charSet[index]) == YES)
                    tempSet[i] = 1;
                }
            fromI = toJ = everyK = -1;
            }

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numChar;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted include set\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                MrBayesPrint ("%s   Improperly formatted include set\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted include set\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                MrBayesPrint ("%s   Improperly formatted include set\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoLog (void)
{
    if (logToFile == YES)
        {
        SafeFclose (&logFileFp);
        if (replaceLogFile == YES)
            {
            if ((logFileFp = OpenTextFileW (logFileName)) == NULL)  
                {
                logToFile = NO;
                return (ERROR);
                }
            }
        else
            {
            if ((logFileFp = OpenTextFileA (logFileName)) == NULL)  
                {
                logToFile = NO;
                return (ERROR);
                }
            }
        MrBayesPrint ("%s   Logging screen output to file \"%s\"\n", spacer, logFileName);
        }
    else
        {
        SafeFclose (&logFileFp);
        MrBayesPrint ("%s   Terminating log output\n", spacer);
        }

    return (NO_ERROR);
}


int DoLogParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Start"))
            {
            if (logToFile == YES)
                MrBayesPrint ("%s   Logging to file is already on\n", spacer, logFileName);
            else
                logToFile = YES;
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        else if (!strcmp(parmName, "Stop"))
            {
            if (logToFile == NO)
                MrBayesPrint ("%s   Logging to file is already off\n", spacer, logFileName);
            else
                logToFile = NO;
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        else if (!strcmp(parmName, "Replace"))
            {
            replaceLogFile = YES;
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        else if (!strcmp(parmName, "Append"))
            {
            replaceLogFile = NO;
            expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
            }
        else
            expecting = Expecting(EQUALSIGN);
        }
    else
        {
        if (!strcmp(parmName, "Filename"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (logFileName, tkn);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else
            {
            MrBayesPrint ("%s   Unknown parameter in Log\n", spacer);
            return (ERROR);
            }
        }

    return (NO_ERROR);
}


int DoManual (void)
{
    int     i, j, logSetting;
    char    title[100];
    FILE    *fp, *logfp;
    CmdType *p;
    
    /* try to open file, return error if present */
    if ((fp = OpenTextFileRQuait(manFileName)) != NULL)
        {
        MrBayesPrint ("%s   File \"%s\" already exists \n", spacer, manFileName);
        SafeFclose(&fp);
        return (ERROR);
        }

    /* try to open file for writing, return error if not possible */
    if ((fp = OpenTextFileW(manFileName)) == NULL)
        return (ERROR);

    /* print message */
    MrBayesPrint ("%s   Producing command reference file \"%s\"\n", spacer, manFileName);

    /* temporarily disable normal logging and switch echoing off */
    logSetting = logToFile;
    logfp = logFileFp;
    echoMB = NO;
    logToFile = YES;
    logFileFp = fp;
    
    /* produce command reference file */
    /* header */
    strcpy (title, "Command Reference for MrBayes ver. ");
    strcat (title, VERSION_NUMBER);

    i = (70 - (int) strlen (title)) / 2;
    j = 70 - i - (int) strlen(title);

    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("      %*c%s%*c      \n", i, ' ', title, j, ' ');
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("                   (c) John P. Huelsenbeck, Fredrik Ronquist                     \n");
    MrBayesPrint ("                               and Maxim Teslenko                                \n");
    MrBayesPrint ("                                                                                 \n");

    /* summary */
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   ***************************************************************************   \n");
    MrBayesPrint ("   *                                                                         *   \n");
    MrBayesPrint ("   *  1. Command summary                                                     *   \n");
    MrBayesPrint ("   *                                                                         *   \n");
    MrBayesPrint ("   ***************************************************************************   \n");
    MrBayesPrint ("                                                                                 \n");
    foundFirst = NO;
    if (DoHelp() == ERROR)
        {
        MrBayesPrint ("%s   Could not produce command reference summary\n", spacer);
        goto errorExit;
        }
    
    /* list of MrBayes commands */
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   ***************************************************************************   \n");
    MrBayesPrint ("   *                                                                         *   \n");
    MrBayesPrint ("   *  2. MrBayes commands                                                    *   \n");
    MrBayesPrint ("   *                                                                         *   \n");
    MrBayesPrint ("   ***************************************************************************   \n");
    MrBayesPrint ("                                                                                 \n");
    for (i=1; i<NUMCOMMANDS; i++)
        {
        p = commands + i;
        if (p->cmdUse == IN_CMD && p->hiding == SHOW)
            {
            if (GetUserHelp(p->string)==ERROR)
                goto errorExit;
            }
        }

    /* list of data or tree block commands */
    MrBayesPrint ("                                                                                 \n");
    MrBayesPrint ("   ***************************************************************************   \n");
    MrBayesPrint ("   *                                                                         *   \n");
    MrBayesPrint ("   *  3. 'Data' or 'tree' block commands (in #NEXUS file)                    *   \n");
    MrBayesPrint ("   *                                                                         *   \n");
    MrBayesPrint ("   ***************************************************************************   \n");
    MrBayesPrint ("                                                                                 \n");
    for (i=1; i<NUMCOMMANDS; i++)
        {
        p = commands + i;
        if (p->cmdUse == IN_FILE && p->hiding == SHOW)
            {
            if (GetUserHelp(p->string) == ERROR)
                goto errorExit;
            }
        }

    /* return logging to previous setings and switch echoing on */
    SafeFclose (&fp);
    logToFile = logSetting;
    logFileFp = logfp;
    echoMB = YES;

    MrBayesPrint ("%s   Successfully produced command reference file \"%s\"\n", spacer, manFileName);

    return (NO_ERROR);

    errorExit:
        SafeFclose (&fp);
        logToFile = logSetting;
        logFileFp = logfp;
        echoMB = YES;

        return (ERROR);
}


int DoManualParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        if (!strcmp(parmName, "Filename"))
            {
            if (expecting == Expecting(EQUALSIGN))
                {
                expecting = Expecting(ALPHA);
                readWord = YES;
                }
            else if (expecting == Expecting(ALPHA))
                {
                strcpy (manFileName, tkn);
                expecting = Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else
            {
            MrBayesPrint ("%s   Unknown parameter in Manual\n", spacer);
            return (ERROR);
            }
        }

    return (NO_ERROR);
}


int DoMatrix (void)
{
    int         i, j, hasMissingAmbig;
    
    if (taxonCount != numTaxa)
        {
        MrBayesPrint ("%s   Problem with number of taxa read in (%d taxa read in, while expecting %d)\n", spacer, taxonCount, numTaxa);
        FreeMatrix();
        return (ERROR);
        }
    for (i=0; i<numTaxa; i++)
        {
        if (taxaInfo[i].charCount != numChar)
            {
            MrBayesPrint ("%s   Problem with number of characters read in (%d expected for taxon %d, %d read in)\n", spacer, numChar, i, taxaInfo[i].charCount);
            FreeMatrix();
            return (ERROR);
            }
        }
        
    /* find out which characters have missing or ambiguous states (one time only, so no special function) */
    for (i=0; i<numChar; i++)
        {
        hasMissingAmbig = NO;
        for (j=0; j<numTaxa; j++)
            {
            if (IsMissing (matrix[pos(j,i,numChar)], charInfo[i].charType) == YES)
                hasMissingAmbig = YES;
            if (IsAmbig (matrix[pos(j,i,numChar)], charInfo[i].charType) == YES)
                hasMissingAmbig = YES;
            }
        if (hasMissingAmbig == YES)
            charInfo[i].isMissAmbig = YES;
        }

    MrBayesPrint ("%s   Successfully read matrix\n", spacer);
    if (matrixHasPoly == YES)
        MrBayesPrint ("%s   Matrix contains polymorphisms, interpreted as ambiguity\n", spacer);
    defMatrix = YES;
    isTaxsetDef = YES;

    /* add name of default partition */
    if (AddString (&partitionNames, 0, "Default") == ERROR)
        {
        MrBayesPrint ("%s   Problem adding Default name to partition list\n", spacer);
        return (ERROR);
        }
    numDefinedPartitions = 1;

    if (numDefinedSpeciespartitions == 0)   /* the default species partition could have been added already in DoTaxLabels */
        {
        /* add default speciespartition name to list of valid speciespartitions */
        if (AddString (&speciespartitionNames, 0, "Default") == ERROR)
            {
            MrBayesPrint ("%s   Problem adding Default speciespartition to list\n", spacer);
            return (ERROR);
            }

        /* add default species name set */
        AddNameSet(&speciesNameSets, 0, taxaNames, numTaxa);

        /* set number of defined speciespartitions to 1 */
        numDefinedSpeciespartitions = 1;
        }
        
    if (SetPartition (0) == ERROR)
        return ERROR;

    if (SetSpeciespartition (0) == ERROR)
        return ERROR;

    if (numCurrentDivisions == 1)
        MrBayesPrint ("%s   Setting default partition (does not divide up characters)\n", spacer);
    else
        MrBayesPrint ("%s   Setting default partition, dividing characters into %d parts\n", spacer, numCurrentDivisions);
    
    if (SetModelDefaults () == ERROR)
        return (ERROR);

    if (SetUpAnalysis (&globalSeed) == ERROR)
        return (ERROR);

    /* set default names for some output file names based on processed file */
    strcpy (sumtParams.sumtFileName, inputFileName);
    strcpy (sumtParams.sumtOutfile, inputFileName);
    strcpy (sumpParams.sumpFileName, inputFileName);
    strcpy (sumpParams.sumpOutfile, inputFileName);
    strcpy (comptreeParams.comptOutfile, inputFileName);

    if (chainParams.numRuns == 1)
        {
        sprintf (comptreeParams.comptFileName1, "%s.t", inputFileName);
        sprintf (comptreeParams.comptFileName2, "%s.t", inputFileName);
        }
    else /* if (chainParams.numRuns > 1) */
        {
        sprintf (comptreeParams.comptFileName1, "%s.run1.t", inputFileName);
        sprintf (comptreeParams.comptFileName2, "%s.run2.t", inputFileName);
        }

    if (chainParams.numRuns == 1)
        sprintf (plotParams.plotFileName, "%s.p", inputFileName);
    else /* if (chainParams.numRuns > 1) */
        sprintf (plotParams.plotFileName, "%s.run1.p", inputFileName);

    strcpy (chainParams.chainFileName, inputFileName);

    if (chainParams.numRuns > 1)
        MrBayesPrint ("%s   Setting output file names to \"%s.run<i>.<p|t>\"\n", spacer, chainParams.chainFileName);
    else
        MrBayesPrint ("%s   Setting output file names to \"%s.<p|t>\"\n", spacer, chainParams.chainFileName);

#   if 0
    for (i=0; i<numChar; i++)
        {
        int     j;
        MrBayesPrint ("%4d -- ", i+1);
        for (j=0; j<numTaxa; j++)
            MrBayesPrint ("%2d ", matrix[pos(j,i,numChar)]);
        MrBayesPrint ("\n");
        }
#   endif

    return (NO_ERROR);
}


int DoMatrixParm (char *parmName, char *tkn)
{
    int             i, j, charCode=0, index;
    MrBFlt          charValue;

    expecting  = Expecting(ALPHA);
    expecting |= Expecting(QUESTIONMARK);
    expecting |= Expecting(DASH);
    expecting |= Expecting(NUMBER);
    expecting |= Expecting(ASTERISK);
    expecting |= Expecting(EXCLAMATIONMARK);
    expecting |= Expecting(PERCENT);
    expecting |= Expecting(WEIRD);
    expecting |= Expecting(VERTICALBAR);
    expecting |= Expecting(SEMICOLON);
    expecting |= Expecting(LEFTPAR);
    expecting |= Expecting(RIGHTPAR);
    expecting |= Expecting(LEFTCURL);
    expecting |= Expecting(RIGHTCURL);

    if (defTaxa == NO || defChars == NO)
        {
        MrBayesPrint ("%s   Number of taxa and characters needs to be defined before matrix is read\n", spacer);
        goto errorExit;
        }
    if (inDataBlock == NO && inCharactersBlock == NO)
        {
        MrBayesPrint ("%s   Must be in data or characters block to read in character matrix\n", spacer);
        goto errorExit;
        }

    if (isFirstMatrixRead == YES)
        {
        foundNewLine = YES;
        isFirstInterleavedBlock = YES;
        taxonCount = 0;
        isNegative = NO;
        }
    isFirstMatrixRead = NO;
    
    /* allow line breaks in non-interleaved matrices */
    if (isInterleaved == NO)
        {
        if (foundNewLine == YES && taxonCount > 0)
            {
            if (taxaInfo[taxonCount-1].charCount < numChar)
                foundNewLine = NO;
            }
        }

    if (taxonCount >= numTaxa && foundNewLine == YES)
        {
        if (isInterleaved == YES)
            {
            taxonCount = 0;
            isFirstInterleavedBlock = NO;
            }
        else
            {
            MrBayesPrint ("%s   Too many taxa in matrix\n", spacer);
            goto errorExit;
            }
        }
    
    if (taxaInfo[0].charCount > 4010)
        i = 1;

    if (foundNewLine == YES)
        {
        /* Should be a taxon. */
        if (isFirstInterleavedBlock == YES)
            {
            /* If this is the first interleaved block, then we need to add the taxon
               to the set of taxon names unless there is already a defined taxon set. */
            if (strlen(tkn)>99)
                {
                MrBayesPrint ("%s   Taxon name %s is too long. Maximun 99 characters is allowed.\n", spacer, tkn);
                goto errorExit;
                }
            if (isTaxsetDef == NO && AddString (&taxaNames, taxonCount, tkn) == ERROR)
                {
                MrBayesPrint ("%s   Problem adding taxon %s to taxon set\n", spacer, tkn);
                goto errorExit;
                }
            if (numTaxa < 10)
                MrBayesPrint ("%s   Taxon %d -> %s\n", spacer, taxonCount+1, tkn);
            else if (numTaxa < 100 && numTaxa >= 10)
                MrBayesPrint ("%s   Taxon %2d -> %s\n", spacer, taxonCount+1, tkn);
            else if (numTaxa < 1000 && numTaxa >= 100)
                MrBayesPrint ("%s   Taxon %3d -> %s\n", spacer, taxonCount+1, tkn);
            else
                MrBayesPrint ("%s   Taxon %4d -> %s\n", spacer, taxonCount+1, tkn);
            }
        else
            {
            /* If this is not the first interleaved block, then we need to
               check to see if taxon name is present and in correct place. */
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                MrBayesPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
                goto errorExit;
                }
            if (index != taxonCount)
                {
                MrBayesPrint ("%s   Could not find taxon %s in correct position in list of taxa\n", spacer, tkn);
                goto errorExit;
                }
            }
        foundNewLine = NO;
        isNegative = NO;
        taxonCount++;
        }
    else
        {
        /* Should be a character (either continuous or otherwise). */
        if (charInfo[taxaInfo[taxonCount-1].charCount].charType == CONTINUOUS)
            {
            /* If we have a CONTINUOUS character, then the entire token should either be
               a number or a dash (for a negative sign). */
            if (!strcmp(tkn, "-"))
                {
                /* Dealing with a negative number. We will multiply the next tkn, which
                   had better be a number, by -1. */
                isNegative = YES;
                }
            else
                {
                /* We have a number, we hope. */
                if (tkn[0] == matchId)
                    {
                    /* If the token is a matchchar, then things are simple. */
                    if (taxonCount == 1)
                        {
                        MrBayesPrint ("%s   Matching characters cannot be in first taxon\n", spacer);
                        goto errorExit;
                        }
                    charCode = matrix[pos(0,taxaInfo[taxonCount-1].charCount,numChar)];
                    matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount,numChar)] = charCode;
                    }
                else
                    {
                    /* Otherwise, we have a number. Check that it is a valid number first... */
                    if (!IsIn(tkn[0],"0123456789."))
                        {
                        MrBayesPrint ("%s   Expecting a number for the continuous character\n", spacer);
                        goto errorExit;
                        }
                    /* ... and then put the character into the matrix. Note that matrix
                       is defined as an integer, but we may have floating precision continuous
                       characters. To get around this, we multiply the value of the character
                       by 1000 before putting it into matrix. We will divide by 1000 later on
                       when/if we use the characters. */
                    sscanf (tkn, "%lf", &charValue);
                    charValue *= 1000.0;
                    if (isNegative == YES)
                        {
                        charValue *= -1.0;
                        isNegative = NO;
                        }
                    /*MrBayesPrint ("%d \n", (int)charValue);*/
                    matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = (int)charValue;
                    }
                }
            }
        else
            {
            /* Otherwise, we are dealing with a run-of-the-mill character, and we
               cannot expect the entire token to contain only a single character. We
               must, therefore, go through the token character-by-character. */
            i = 0;
            while (tkn[i] != '\0')
                {
                /*MrBayesPrint ("%c", tkn[i]);*/
                if (tkn[i] == matchId)
                    {
                    if (taxonCount == 1)
                        {
                        MrBayesPrint ("%s   Matching characters cannot be in first taxon\n", spacer);
                        goto errorExit;
                        }
                    charCode = matrix[pos(0,taxaInfo[taxonCount-1].charCount,numChar)];
                    matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
                    }
                else
                    {
                    if ((tkn[i] == ')' && isInAmbig == YES) || (tkn[i] == '}' && isInPoly == YES))
                        {
                        isInAmbig = isInPoly = NO;
                        charCode = theAmbigChar;
                        j = CharacterNumber (charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType);
                        if (j > charInfo[taxaInfo[taxonCount-1].charCount].numStates)
                            charInfo[taxaInfo[taxonCount-1].charCount].numStates = j;
                        matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
                        theAmbigChar = 0;
                        }
                    else if ((tkn[i] == '(' && isInAmbig == YES) || (tkn[i] == '{' && isInPoly == YES))
                        {
                        if (isInAmbig == YES)
                            MrBayesPrint ("%s   Found an inappropriate \"(\"\n", spacer);
                        else
                            MrBayesPrint ("%s   Found an inappropriate \"{\"\n", spacer);
                        goto errorExit;
                        }
                    else if (isInAmbig == YES || isInPoly == YES)
                        {
                        if (tkn[i] == ',')
                            expecting |= Expecting (COMMA);
                        else 
                            {
                            if (CharacterCode(tkn[i], &charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType) == ERROR)
                                goto errorExit;
                            if (charCode == MISSING || charCode == GAP)
                                goto errorExit;
                            theAmbigChar |= charCode;
                            expecting ^= Expecting (COMMA);
                            }
                        }
                    else if (tkn[i] == '{' && isInPoly == NO && isInAmbig == NO)
                        {
                        isInPoly = YES;
                        matrixHasPoly = YES;
                        theAmbigChar = 0;
                        }   
                    else if (tkn[i] == '(' && isInPoly == NO && isInAmbig == NO)
                        {
                        isInAmbig = YES;
                        theAmbigChar = 0;
                        }
                    else if (tkn[i] == '(' && isInPoly == NO && isInAmbig == NO)
                        {
                        isInAmbig = YES;
                        theAmbigChar = 0;
                        }
                    else
                        {
                        if (CharacterCode(tkn[i], &charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType) == ERROR)
                            {
                            MrBayesPrint ("%s   Error while reading character position %d (charCode %d)\n", spacer, taxaInfo[taxonCount-1].charCount+1, charCode);
                            goto errorExit;
                            }
                        if (charCode != MISSING && charCode != GAP)
                            {
                            j = CharacterNumber (charCode, charInfo[taxaInfo[taxonCount-1].charCount].charType);
                            if (j > charInfo[taxaInfo[taxonCount-1].charCount].numStates)
                                charInfo[taxaInfo[taxonCount-1].charCount].numStates = j;
                            }
                        matrix[pos(taxonCount-1,taxaInfo[taxonCount-1].charCount++,numChar)] = charCode;
                        }
                    }
                i++; 
                }
            }
        }

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
    errorExit:
        numTaxa=taxonCount;
        FreeMatrix();
        return (ERROR);
}


int DoNexusParm (char *parmName, char *tkn)
{
    if (!strcmp(parmName, "NEXUS"))
        {
        MrBayesPrint ("%s   Expecting NEXUS formatted file\n", spacer);
        expecting = Expecting(COMMAND);
        }
    else
        {
        MrBayesPrint ("%s   Found %s\n", spacer, tkn);
        return (ERROR);
        }
    
    return (NO_ERROR);
}


int DoOutgroup (void)
{
    MrBayesPrint ("%s   Setting outgroup to taxon \"%s\"\n", spacer, taxaNames[outGroupNum]);
    return (NO_ERROR);
}


int DoOutgroupParm (char *parmName, char *tkn)
{
    int     index, tempInt;

    if (expecting == Expecting(ALPHA))
        {
        if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
            {
            MrBayesPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
            return (ERROR);
            }
        outGroupNum = index;
        
        expecting = Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
            {
            /* OK, as we expect, the taxon is not a digit. So, now we assume that
               the user is assigning the outgroup by its number */
            sscanf (tkn, "%d", &tempInt);
            if (tempInt < 1 || tempInt > numTaxa)
                {
                MrBayesPrint ("%s   Taxon number %d is out of range\n", spacer, tempInt);
                return (ERROR);
                }
            outGroupNum = tempInt - 1;
            }
        else
            {
            outGroupNum = index;
            }
        
        expecting = Expecting(SEMICOLON);
        }
    else
        return (ERROR);

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoPairs (void)
{
    MrBayesPrint ("\n");
    MrBayesPrint ("%s   Successfully defined character pairings\n", spacer);

    defPairs = YES;
    foundFirst = NO;
    
    return (NO_ERROR);
}


int DoPairsParm (char *parmName, char *tkn)
{
    int     i, tempInt;
        
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before you can define pairs of characters\n", spacer);
        return (ERROR);
        }
    
    if (defPairs == YES)
        {
        MrBayesPrint ("%s   Character pairs have been previously defined \n", spacer);
        MrBayesPrint ("%s   Now overwriting old pairings\n", spacer);
        for (i=0; i<numChar; i++)
            charInfo[i].pairsId = 0;
        defPairs = NO;
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        pairId = 1;
        firstPair = YES;
        foundFirst = YES;
        MrBayesPrint ("%s   Defining character pairings:\n\n", spacer);
        MrBayesPrint ("%s      Pair --  First Second \n", spacer);
        }

    if (expecting == Expecting(NUMBER))
        {
        sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numChar)
            {
            MrBayesPrint ("\n");
            MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
            for (i=0; i<numChar; i++)
                charInfo[i].pairsId = 0;
            return (ERROR);
            }
        tempInt--;
        
        if (charInfo[tempInt].pairsId != 0)
            {
            MrBayesPrint ("\n");
            MrBayesPrint ("%s   Character number %d has already been included in a pairing\n", spacer, tempInt+1);
            for (i=0; i<numChar; i++)
                charInfo[i].pairsId = 0;
            return (ERROR);
            }
        if (charInfo[tempInt].charType != DNA && charInfo[tempInt].charType != RNA)
            {
            MrBayesPrint ("\n");
            MrBayesPrint ("%s  Pairings may only include nucleotide data\n", spacer);
            if (charInfo[tempInt].charType == PROTEIN)
                MrBayesPrint ("%s  Character %d is an amino acid character\n", spacer, tempInt+1);
            else if (charInfo[tempInt].charType == RESTRICTION)
                MrBayesPrint ("%s  Character %d is a restriction site character\n", spacer, tempInt+1);
            else if (charInfo[tempInt].charType == STANDARD)
                MrBayesPrint ("%s  Character %d is a \"standard\" character\n", spacer, tempInt+1);
            else if (charInfo[tempInt].charType == CONTINUOUS)
                MrBayesPrint ("%s  Character %d is a continuously varying character\n", spacer, tempInt+1);
            for (i=0; i<numChar; i++)
                charInfo[i].pairsId = 0;
            return (ERROR);
            }
            
        charInfo[tempInt].pairsId = pairId;
        
        if (firstPair == YES)
            {
            MrBayesPrint ("%s      %4d --  %5d  ", spacer, pairId, tempInt+1);
            expecting  = Expecting(COLON);
            firstPair = NO;
            }
        else
            {
            MrBayesPrint ("%5d\n", tempInt+1);
            expecting  = (Expecting(COMMA) | Expecting(SEMICOLON));
            firstPair = YES;
            }
        }
    else if (expecting == Expecting(COMMA))
        {
        pairId++;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting = Expecting(NUMBER);
        }
    else
        {
        for (i=0; i<numChar; i++)
            charInfo[i].pairsId = 0;
        return (ERROR);
        }

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoPartition (void)
{
    int     i, *partTypes;
        
    /* add set to tempSet */
    if (fromI >= 0)
        if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
            return (ERROR);

    /* check that all characters are included */
    for (i=0; i<numChar; i++)
        {
        /* MrBayesPrint ("%4d %4d \n", i, tempSet[i]); */
        if (tempSet[i] == 0)
            {
            MrBayesPrint ("%s   Character %d not included in partition\n", spacer, i+1);
            return (ERROR);
            }
        }

            
    /* check how many partitions were found against how many were expected */
    if (whichPartition != numDivisions - 1)
        {
        MrBayesPrint ("%s   Did not find correct number of partitions (expecting %d, found %d)\n", spacer, numDivisions, whichPartition + 1);
        return (ERROR);
        }

    partTypes = (int *) SafeCalloc (numDivisions, sizeof(int));
    if (!partTypes)
        return ERROR;
    
    /* make certain that the partition labels go from 1 - numDivisions, inclusive */
    for (i=0; i<numChar; i++)
        partTypes[tempSet[i] - 1] = -1; //partTypes is temporary used here not as an indicator of partition type 
    for (i=0; i<numDivisions; i++)
        {
        if (partTypes[i] == 0)
            {
            MrBayesPrint ("%s   Could not find a single character for division %d\n", spacer, i+1);
            return (ERROR);
            }
        }

    /* check if partition overruns data types */
    for (i=0; i<numChar; i++)
        {
        if (partTypes[ tempSet[i]-1 ] == -1)
            partTypes[ tempSet[i]-1 ] = charInfo[i].charType;
        else
            {
            if (partTypes[ tempSet[i]-1 ] != charInfo[i].charType)
                {
                MrBayesPrint ("%s   There are two different data types for partition division %d\n", spacer, tempSet[i]);
                free (partTypes);
                return (ERROR);
                }
            }
        }
    free (partTypes);

    /* add name to list of valid partitions */
    if (AddString (&partitionNames, numDefinedPartitions, tempSetName) == ERROR)
        {
        MrBayesPrint ("%s   Problem adding partition %s to list\n", spacer, tempSetName);
        return (ERROR);
        }
        
    /* add new partition */
    for (i=0; i<numChar; i++) {
        partitionId[i] = (int *) SafeRealloc ((void *)(partitionId[i]), ((size_t)numDefinedPartitions + 1) * sizeof(int));
        if (!partitionId[i])
            return ERROR;
    }

    /* set new partition */
    for (i=0; i<numChar; i++)
        partitionId[i][numDefinedPartitions] = tempSet[i];

    /* increment number of defined partitions */
    numDefinedPartitions++;
    
    return (NO_ERROR);
}


int DoPartitionParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
    
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before partitions can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        /* set Partition () ******************************************************************/
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check size of partition name */
            if (strlen(tkn) > 99)
                {
                MrBayesPrint ("%s   Partition name is too long. Max 100 characters\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a partition */
            if (numDefinedPartitions > 1)
                {
                if (CheckString (partitionNames, numDefinedPartitions, tkn, &index) == ERROR)
                    {
                    /* if the partition name has not been used, then we should have an ERROR returned */
                    /* we _want_ to be here */

                    }
                else
                    {
                    MrBayesPrint ("%s   Partition name '%s' has been used previously\n", spacer, tkn);
                    return (ERROR);
                    }
                }
                
            /* add the name temporarily to tempSetName */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numChar; i++)
                tempSet[i] = 0;
            
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            whichPartition = 0;
            foundFirst = NO;
            numDivisions = 0;
            MrBayesPrint ("%s   Defining partition called '%s'\n", spacer, tkn);
            expecting = Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        /* We are defining a partition in terms of a character set (called tkn, here). We should be able
           to find tkn in the list of character set names. If we cannot, then we have a problem and
           return an error. */
        if (numCharSets < 1)
            {
            MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
            return (ERROR);
            }
        if (CheckString (charSetNames, numCharSets, tkn, &index) == ERROR)
            {
            MrBayesPrint ("%s   Could not find a character set called '%s'\n", spacer, tkn);
            return (ERROR);
            }
        /* add characters from charset tkn to new tempSet */
        for (i=0; i<numChar; i++)
            {
            if (IsBitSet (i, charSet[index]) == YES)
                tempSet[i] = whichPartition + 1;
            }
        fromI = toJ = everyK = -1;

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(COMMA);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (foundFirst == NO)
            {
            sscanf (tkn, "%d", &tempInt);
            numDivisions = tempInt;
            expecting  = Expecting(COLON);
            foundFirst = YES;
            }
        else
            {
            if (strlen(tkn) == 1 && tkn[0] == '.')
                tempInt = numChar;
            else
                sscanf (tkn, "%d", &tempInt);
            if (tempInt <= 0 || tempInt > numChar)
                {
                MrBayesPrint ("%s   Character number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numChar);
                return (ERROR);
                }
            tempInt--;
            if (foundDash == YES)
                {
                if (fromI >= 0)
                    toJ = tempInt;
                else
                    {
                    MrBayesPrint ("%s   Improperly formatted partition\n", spacer);
                    return (ERROR);
                    }
                foundDash = NO;
                }
            else if (foundSlash == YES)
                {
                tempInt++;
                if (tempInt <= 1)
                    {
                    MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
                    return (ERROR);
                    }
                if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                    everyK = tempInt;
                else
                    {
                    MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
                    return (ERROR);
                    }
                foundSlash = NO;
                }
            else
                {
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    }
                else if (fromI < 0 && toJ < 0)
                    {
                    fromI = tempInt;
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    toJ = everyK = -1;
                    }
                else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                        return (ERROR);
                    fromI = tempInt;
                    toJ = everyK = -1;
                    }
                else
                    {
                    MrBayesPrint ("%s   Improperly formatted charset\n", spacer);
                        {
                        return (ERROR);
                        }
                    }
                }
            
            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(SEMICOLON);
            expecting |= Expecting(DASH);
            expecting |= Expecting(BACKSLASH);
            expecting |= Expecting(COMMA);
            }
        }
    else if (expecting == Expecting(COMMA))
        {
        /* add set to tempSet */
        if (fromI >= 0)
            if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                return (ERROR);

        fromI = toJ = everyK = -1;
        foundDash = foundSlash = NO;
        whichPartition++;
        if (whichPartition > numDivisions)
            {
            MrBayesPrint ("%s   Too many partitions of the data (expecting %d)\n", spacer, numDivisions);
            return (ERROR);
            }
        expecting  = Expecting(NUMBER);
        expecting |= Expecting(ALPHA);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting  = Expecting(NUMBER);
        expecting |= Expecting(ALPHA);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoRestore (void)
{
    int         i, alreadyDone;

    MrBayesPrint ("%s   Restore taxa\n", spacer);

    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* merge tempSet with excludedTaxa */
    alreadyDone = NO;
    for (i=0; i<numTaxa; i++)
        {
        if (tempSet[i] == 1)
            {
            if (taxaInfo[i].isDeleted == NO && alreadyDone == NO)
                {
                MrBayesPrint ("%s   Some taxa already included\n", spacer);
                alreadyDone = YES;
                }
            taxaInfo[i].isDeleted = NO;
            }
        }

    SetLocalTaxa();
    if (SetUpAnalysis(&globalSeed) == ERROR)
        return ERROR;

    /* show tempSet (for debugging) */
#   if 0
    for (i=0; i<numTaxa; i++)
        MrBayesPrint ("%4d  %4d\n", i+1, tempSet[i]);
#   endif
        
    return (NO_ERROR);
}


int DoRestoreParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
        
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before you can restore taxa\n", spacer);
        return (ERROR);
        }
        
    if (foundFirst == NO)
        {
        /* this is the first time in */
        fromI = toJ = everyK = -1;
        foundDash = NO;
        for (i=0; i<numTaxa; i++) /* clear tempSet */
            tempSet[i] = 0;
        foundFirst = YES;
        }

    if (expecting == Expecting(ALPHA))
        {
        if (IsSame ("All", tkn) == SAME || IsSame ("All", tkn) == CONSISTENT_WITH)
            {
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 1;
            }
        else
            {
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                /* we are using a pre-defined taxa set */
                if (numTaxaSets < 1)
                    {
                    MrBayesPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                    return (ERROR);
                    }
                /* add taxa from taxset tkn to new tempSet */
                for (i=0; i<numTaxa; i++)
                    {
                    if (IsBitSet (i, taxaSet[index]) == YES)
                        tempSet[i] = 1;
                    }
                }
            else
                {
                /* we found the taxon name */
                if (fromI >= 0 && toJ < 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                        return (ERROR);
                    }
                else if (fromI >= 0 && toJ >= 0)
                    {
                    if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                        return (ERROR);
                    }
                tempSet[index] = 1;
                }
            fromI = toJ = everyK = -1;
            }

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && !strcmp(tkn, "."))
            {
            tempInt = numTaxa;
            }
        else
            {
            sscanf (tkn, "%d", &tempInt);
            if (tempInt <= 0 || tempInt > numTaxa)
                {
                MrBayesPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
                return (ERROR);
                }
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted restore set\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                MrBayesPrint ("%s   Improperly formatted restore set\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoSet (void)
{
    return (NO_ERROR);
}


int DoSetParm (char *parmName, char *tkn)
{
    int         index;
    char        tempStr[100];
    int         tempI;

    if (expecting == Expecting(PARAMETER))
        {
        expecting = Expecting(EQUALSIGN);
        }
    else
        {
        /* set Autoclose (autoClose) **********************************************************/
        if (!strcmp(parmName, "Autoclose"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        autoClose = YES;
                    else
                        autoClose = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for autoclose\n", spacer);
                    return (ERROR);
                    }
                if (autoClose == YES)
                    MrBayesPrint ("%s   Setting autoclose to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting autoclose to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Nowarnings (noWarn) **********************************************************/
        else if (!strcmp(parmName, "Nowarnings"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        noWarn = YES;
                    else
                        noWarn = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for nowarnings\n", spacer);
                    return (ERROR);
                    }
                if (noWarn == YES)
                    MrBayesPrint ("%s   Setting nowarnings to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting nowarnings to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Quitonerror (quitOnError) **************************************************/
        else if (!strcmp(parmName, "Quitonerror"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        quitOnError = YES;
                    else
                        quitOnError = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for quitonerror\n", spacer);
                    return (ERROR);
                    }
                if (quitOnError == YES)
                    MrBayesPrint ("%s   Setting quitonerror to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting quitonerror to no\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Autoreplace (autoOverwrite) **************************************************/
        else if (!strcmp(parmName, "Autoreplace"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        {
                        autoOverwrite = YES;
                        MrBayesPrint ("%s   Setting autoreplace to yes\n", spacer);
                        }
                    else
                        {
                        autoOverwrite = NO;
                        MrBayesPrint ("%s   Setting autoreplace to no\n", spacer);
                        }
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for autoreplace\n", spacer);
                    return (ERROR);
                    }                   
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }           
        /* set Scientific (scientific) *********************************************/
        else if (!strcmp(parmName, "Scientific"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        scientific = YES;
                    else
                        scientific = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for Scientific\n", spacer);
                    return (ERROR);
                    }
                if (scientific == YES)
                    MrBayesPrint ("%s   Setting Scientific to Yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting Scientific to No\n", spacer);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Userlevel (userLevel) **********************************************************/
        else if (!strcmp(parmName, "Userlevel"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Standard"))
                        userLevel = STANDARD_USER;
                    else if (!strcmp (tempStr,"Developer"))
                        userLevel = DEVELOPER;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for userlevel\n", spacer);
                    return (ERROR);
                    }
                MrBayesPrint ("%s   Setting userlevel to %s\n", spacer, tempStr);
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Npthreads (number of pthreads) ****************************************************/
        else if (!strcmp(parmName, "Npthreads"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                nPThreads = tempI;
                MrBayesPrint ("%s   Setting Npthreads to %d\n", spacer, nPThreads);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                return (ERROR);
            }
        /* set Precision (number of decimals) ****************************************************/
        else if (!strcmp(parmName, "Precision"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI < 3 || tempI > 15)
                    {
                    MrBayesPrint ("%s   Precision must be in the range 3 to 15\n", spacer);
                    return ERROR;
                    }
                precision = tempI;
                MrBayesPrint ("%s   Setting Precision to %d\n", spacer, precision);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                return (ERROR);
            }
        /* set Partition (partitionNum) *******************************************************/
        else if (!strcmp(parmName, "Partition"))
            {
            if (defMatrix == NO)
                {
                MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
                return (ERROR);
                }
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA) | Expecting(NUMBER);
            else if (expecting == Expecting(ALPHA))
                {
                /* first check to see if name is there */
                if (CheckString (partitionNames, numDefinedPartitions, tkn, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find \"%s\" as a defined partition\n", spacer, tkn);
                    return (ERROR);
                    }
                if (SetPartition (index) == ERROR)
                    return ERROR;
                if (numCurrentDivisions == 1)
                    MrBayesPrint ("%s   Setting %s as the partition (does not divide up characters).\n", spacer, tkn); 
                else
                    MrBayesPrint ("%s   Setting %s as the partition, dividing characters into %d parts.\n", spacer, tkn, numCurrentDivisions); 
                if (SetModelDefaults () == ERROR)
                    return (ERROR);
                if (SetUpAnalysis (&globalSeed) == ERROR)
                    return (ERROR);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &index);
                if (index > numDefinedPartitions) 
                    {
                    MrBayesPrint ("%s   Partition number %d is not a valid partition. Only %d partitions\n", spacer, index, numDefinedPartitions);
                    MrBayesPrint ("%s   have been defined.\n", spacer);
                    return (ERROR);
                    }
                if (index < 1)
                    {
                    MrBayesPrint ("%s   Partition number %d is not a valid partition. Must be between 1 and %d.\n", spacer, index+1, numDefinedPartitions);
                    return (ERROR);
                    }
                if (SetPartition (index) == ERROR)
                    return ERROR;
                if (numCurrentDivisions == 1)
                    MrBayesPrint ("%s   Setting %s as the partition (does not divide up characters).\n", spacer, partitionNames[index]);
                else
                    MrBayesPrint ("%s   Setting %s as the partition, dividing characters into %d parts.\n", spacer, partitionNames[index], numCurrentDivisions);
                if (SetModelDefaults () == ERROR)
                    return (ERROR);
                if (SetUpAnalysis (&globalSeed) == ERROR)
                    return (ERROR);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Speciespartition (speciespartitionNum) *******************************************************/
        else if (!strcmp(parmName, "Speciespartition"))
            {
            if (defTaxa == NO)
                {
                MrBayesPrint ("%s   A taxaset must be defined first\n", spacer);
                return (ERROR);
                }
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA) | Expecting(NUMBER);
            else if (expecting == Expecting(ALPHA))
                {
                /* first check to see if name is there */
                if (CheckString (speciespartitionNames, numDefinedSpeciespartitions, tkn, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find \"%s\" as a defined speciespartition\n", spacer, tkn);
                    return (ERROR);
                    }
                if (SetSpeciespartition (index) == ERROR)
                    return ERROR;
                MrBayesPrint ("%s   Setting %s as the speciespartition, dividing taxa into %d species.\n", spacer, tkn, numSpecies);
                if (SetModelDefaults () == ERROR)
                    return (ERROR);
                if (SetUpAnalysis (&globalSeed) == ERROR)
                    return (ERROR);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &index);
                if (index > numDefinedSpeciespartitions) 
                    {
                    MrBayesPrint ("%s   Speciespartition number %d is not valid. Only %d speciespartitions\n", spacer, index, numDefinedSpeciespartitions);
                    MrBayesPrint ("%s   have been defined.\n", spacer);
                    return (ERROR);
                    }
                if (index < 1)
                    {
                    MrBayesPrint ("%s   Speciespartition number %d is not valid. Must be between 1 and %d.\n", spacer, index, numDefinedSpeciespartitions);
                    return (ERROR);
                    }
                if (SetSpeciespartition (index-1) == ERROR)
                    return ERROR;
                MrBayesPrint ("%s   Setting %s as the speciespartition, dividing taxa into %d species.\n", spacer, speciespartitionNames[index-1], numSpecies);
                if (SetModelDefaults () == ERROR)
                    return (ERROR);
                if (SetUpAnalysis (&globalSeed) == ERROR)
                    return (ERROR);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Seed (global variable globalSeed) ****************************************************/
        else if (!strcmp(parmName, "Seed"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI == 0 || tempI == 2147483647)
                    {
                    MrBayesPrint ("%s   Error: Seed can be any natural number except 0 and 2147483647\n", spacer);
                    return (ERROR);
                    }
                globalSeed = tempI;
                MrBayesPrint ("%s   Setting seed to %ld\n", spacer, globalSeed);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                return (ERROR);
            }
        /* set Swapseed (global variable swapSeed) ***************************************************************/
        else if (!strcmp(parmName, "Swapseed"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
                sscanf (tkn, "%d", &tempI);
                if (tempI == 0 || tempI == 2147483647)
                    {
                    MrBayesPrint ("%s   Error: Swapseed can be any natural number except 0 and 2147483647\n", spacer);
                    return (ERROR);
                    }
                swapSeed = tempI;
                MrBayesPrint ("%s   Setting swapseed to %ld\n", spacer, swapSeed);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Dir (global variable workingDir) ***************************************************************/
        else if (!strcmp(parmName, "Dir"))
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
                    MrBayesPrint ("%s   Maximum allowed length of working directory name is 99 characters. The given name:\n", spacer);
                    MrBayesPrint ("%s      '%s'\n", spacer,tkn);
                    MrBayesPrint ("%s   has %d characters.\n", spacer,strlen(tkn));
                    return (ERROR);
                    }
                strcpy (workingDir, tkn);
#   if defined (WIN_VERSION)
                /* Reformat to Windows with trailing '\' */
                for (index=0; index<(int)strlen(workingDir); index++)
                    {
                    if (workingDir[index] == '/')
                        workingDir[index] = '\\';
                    }
                if (strlen(workingDir) > 0 && workingDir[strlen(workingDir)-1] != '\\')
                    strcat(workingDir,"\\");
#   else
                /* Reformat to Unix with trailing '/' */
                for (index=0; index<(int)strlen(workingDir); index++)
                    {
                    if (workingDir[index] == '\\')
                        workingDir[index] = '/';
                    }
                if (strlen(workingDir) > 0 && workingDir[strlen(workingDir)-1] != '/')
                    strcat(workingDir,"/");
#   endif
                MrBayesPrint ("%s   Setting working directory to \"%s\"\n", spacer, workingDir);
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Usebeagle (global variable BEAGLE usage) ***************************************************************/    
        else if (!strcmp(parmName, "Usebeagle"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
#   if defined (BEAGLE_ENABLED)
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        tryToUseBEAGLE = YES;
                    else
                        tryToUseBEAGLE = NO;
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for usebeagle\n", spacer);
                    return (ERROR);
                    }
                if (tryToUseBEAGLE == YES)
                    MrBayesPrint ("%s   Setting usebeagle to yes\n", spacer);
                else
                    MrBayesPrint ("%s   Setting usebeagle to no\n", spacer);
#   else
                BeagleNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Beagle resource number (global variable BEAGLE flag) ****************************************/
        else if (!strcmp(parmName, "Beagleresource"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting =  Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
#   if defined (BEAGLE_ENABLED)
                sscanf (tkn, "%d", &tempI);
                if (tempI < 0)
                    {
                    MrBayesPrint ("%s   Beagleresource must be a valid resource number or 99 to disable resource selection\n", spacer);
                    return ERROR;
                    }
                beagleResourceNumber = tempI;
                if (beagleResourceNumber == 99)
                    MrBayesPrint ("%s   Setting Beagleresource to %d (auto)\n", spacer, beagleResourceNumber);
                else
                    MrBayesPrint ("%s   Setting Beagleresource to %d\n", spacer, beagleResourceNumber);
#   else
                BeagleNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        /* set Beagle resources requirements (global variable BEAGLE flag) ****************************************/
        else if (!strcmp(parmName, "Beagledevice"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
#   if defined (BEAGLE_ENABLED)
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    long oldFlags = beagleFlags;
                    if (!strcmp(tempStr, "Gpu"))
                        {
                        beagleFlags &= ~BEAGLE_FLAG_PROCESSOR_CPU;
                        beagleFlags |= BEAGLE_FLAG_PROCESSOR_GPU;
                        BeagleAddGPUDevicesToList(&beagleResource, &beagleResourceCount);                       
                        }
                    else
                        {  
                        beagleFlags &= ~BEAGLE_FLAG_PROCESSOR_GPU;
                        beagleFlags |= BEAGLE_FLAG_PROCESSOR_CPU;
                        BeagleRemoveGPUDevicesFromList(&beagleResource, &beagleResourceCount);
                        }
                    if (BeagleCheckFlagCompatability(beagleFlags) == NO) {
                        beagleFlags = oldFlags;
                        }
                    else {
                        if (beagleFlags & BEAGLE_FLAG_PROCESSOR_GPU)
                            MrBayesPrint ("%s   Setting beagledevice to GPU\n", spacer);
                        else
                            MrBayesPrint ("%s   Setting beagledevice to CPU\n", spacer);
                        }
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for beagledevice\n", spacer);
                    return (ERROR);
                    }
#   else
                BeagleNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beagleprecision"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
#   if defined (BEAGLE_ENABLED)
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    long oldFlags = beagleFlags;
                    if (!strcmp(tempStr, "Single"))
                        {
                        beagleFlags &= ~BEAGLE_FLAG_PRECISION_DOUBLE;
                        beagleFlags |= BEAGLE_FLAG_PRECISION_SINGLE;                       
                        }
                    else
                        {
                        beagleFlags &= ~BEAGLE_FLAG_PRECISION_SINGLE;
                        beagleFlags |= BEAGLE_FLAG_PRECISION_DOUBLE;
                        }
                    if (BeagleCheckFlagCompatability(beagleFlags) == NO) {
                        beagleFlags = oldFlags;
                        }
                    else {
                        if (beagleFlags & BEAGLE_FLAG_PRECISION_DOUBLE)
                            MrBayesPrint ("%s   Setting beagleprecision to double\n", spacer);
                        else
                            MrBayesPrint ("%s   Setting beagleprecision to single\n", spacer);
                        }
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for beagleprecision\n", spacer);
                    return (ERROR);
                    }
#   else
                BeagleNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beagleopenmp"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
#   if defined (BEAGLE_ENABLED)
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    long oldFlags = beagleFlags;
                    if (!strcmp(tempStr, "Yes"))
                        {
                        beagleFlags |= BEAGLE_FLAG_THREADING_OPENMP;
                        }
                    else
                        {
                        beagleFlags &= ~BEAGLE_FLAG_THREADING_OPENMP;                       
                        }
                    if (BeagleCheckFlagCompatability(beagleFlags) == NO) {
                        beagleFlags = oldFlags;
                        }
                    else {
                        if (beagleFlags & BEAGLE_FLAG_THREADING_OPENMP)
                            MrBayesPrint ("%s   Setting beagleopenmp to Yes\n", spacer);
                        else
                            MrBayesPrint ("%s   Setting beagleopenmp to No\n", spacer);
                        }
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for beagleopenmp\n", spacer);
                    return (ERROR);
                    }
#   else
                BeagleNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
        else if (!strcmp(parmName, "Beaglefreq"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(NUMBER);
            else if (expecting == Expecting(NUMBER))
                {
#   if defined (BEAGLE_ENABLED)
                sscanf (tkn, "%d", &tempI);
                if (tempI < 0)
                    {
                    MrBayesPrint ("%s   Beaglefreq must be greater than 0\n", spacer);
                    return ERROR;
                    }
                beagleScalingFrequency= tempI;
                MrBayesPrint ("%s   Setting Beaglefreq to %d\n", spacer, beagleScalingFrequency);
#   else
                BeagleNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else 
                return (ERROR);
            }        
        else if (!strcmp(parmName, "Beaglesse"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
#   if defined (BEAGLE_ENABLED)
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    long oldFlags = beagleFlags;
                    if (!strcmp(tempStr, "Yes"))
                        {
                        beagleFlags |= BEAGLE_FLAG_VECTOR_SSE;
                        }
                    else
                        {
                        beagleFlags &= ~BEAGLE_FLAG_VECTOR_SSE;                     
                        }
                    if (BeagleCheckFlagCompatability(beagleFlags) == NO) {
                        beagleFlags = oldFlags;
                        }
                    else {
                        if (beagleFlags & BEAGLE_FLAG_VECTOR_SSE)
                            MrBayesPrint ("%s   Setting beaglesse to Yes\n", spacer);
                        else
                            MrBayesPrint ("%s   Setting beaglesse to No\n", spacer);
                        }
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for beagleopenmp\n", spacer);
                    return (ERROR);
                    }
#   else
                BeagleNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
#if 0
        else if (!strcmp(parmName, "Beaglevec"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
#   if defined (BEAGLE_ENABLED)
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {                    
                    long oldFlags = beagleFlags;
                    if (!strcmp(tempStr, "Sse"))
                        {                      
                        beagleFlags |= BEAGLE_FLAG_VECTOR_SSE;
                        beagleFlags &= ~BEAGLE_FLAG_VECTOR_AVX;
                        }
                    else if (!strcmp(tempStr, "Avx"))
                        {  
                        beagleFlags |= ~BEAGLE_FLAG_VECTOR_AVX;
                        beagleFlags &= ~BEAGLE_FLAG_VECTOR_SSE;
                        }
                    else if (!strcmp(tempStr, "None"))
                        {
                        beagleFlags &= ~BEAGLE_FLAG_VECTOR_SSE;
                        beagleFlags &= ~BEAGLE_FLAG_VECTOR_AVX;
                        }
                    else
                        {
                        MrBayesPrint("%s   Unrecognized argument for beaglevec\n", spacer);
                        }
                    MrBayesPrint ("%s   Setting beaglevec to %s\n", spacer, tempStr);
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for beagleopenmp\n", spacer);
                    return (ERROR);
                    }
#   else
                BeagleNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            else
                return (ERROR);
            }
#endif
        else if (!strcmp(parmName, "Beaglethreads"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
#   if defined (BEAGLE_ENABLED) && defined (THREADS_ENABLED)
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Yes"))
                        {
                        tryToUseThreads = YES;
                        }
                    else
                        {
                        tryToUseThreads = NO;                       
                        }
                    
                    if (tryToUseThreads == YES)
                        MrBayesPrint ("%s   Setting beaglethreads to Yes\n", spacer);
                    else
                        MrBayesPrint ("%s   Setting beaglethreads to No\n", spacer);                    
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for beaglethreads\n", spacer);
                    return (ERROR);
                    }
#   else
                BeagleThreadsNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
                expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
                }
            }
        else if (!strcmp(parmName, "Beaglescaling"))
            {
            if (expecting == Expecting(EQUALSIGN))
                expecting = Expecting(ALPHA);
            else if (expecting == Expecting(ALPHA))
                {
#   if defined (BEAGLE_ENABLED)
                if (IsArgValid(tkn, tempStr) == NO_ERROR)
                    {
                    if (!strcmp(tempStr, "Always"))
                        {
                        beagleScalingScheme = MB_BEAGLE_SCALE_ALWAYS;
                        }
                    else
                        {
                        beagleScalingScheme = MB_BEAGLE_SCALE_DYNAMIC;                      
                        }
                    
                    if (beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS)
                        MrBayesPrint ("%s   Setting beaglescaling to Always\n", spacer);
                    else
                        MrBayesPrint ("%s   Setting beaglescaling to Dynamic\n", spacer);                    
                    }
                else
                    {
                    MrBayesPrint ("%s   Invalid argument for beaglescaling\n", spacer);
                    return (ERROR);
                    }
#   else
                BeagleThreadsNotLinked();
#   endif
                if (defMatrix == YES && SetUpAnalysis(&globalSeed) == ERROR)
                    return ERROR;
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


int DoShowMatrix (void)
{
    int         i, j, nameLen, start, finish, ct, longestName;
    char        tempStr[100], stride;
    
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }
            
    longestName = 0;
    for (i=0; i<numTaxa; i++)
        {
        nameLen = (int) strlen(taxaNames[i]);
        if (nameLen > longestName)
            longestName = nameLen;
        }
            
    stride = 50;
    start = finish = 0;
    do
        {
        finish += stride;
        if (finish > numChar)
            finish = numChar;

        MrBayesPrint ("%s   ", spacer);
        for (j=0; j<longestName; j++)
            MrBayesPrint (" ");
        MrBayesPrint ("  ");
        MrBayesPrint ("%d\n", start+1);

        for (i=0; i<numTaxa; i++)
            {
            strcpy (tempStr, taxaNames[i]);
            nameLen = (int) strlen(tempStr);
            
            MrBayesPrint ("%s   ", spacer);
            if (nameLen >= longestName)
                {
                for (j=0; j<longestName; j++)
                    MrBayesPrint ("%c", tempStr[j]);
                }
            else
                {
                MrBayesPrint ("%s", tempStr);
                for (j=0; j<longestName-nameLen; j++)
                    MrBayesPrint (" ");
                }
            MrBayesPrint ("  ");

            for (j=start; j<finish; j++)
                {
                ct = charInfo[j].charType;
                if (ct == DNA || ct == RNA)
                    MrBayesPrint ("%c", WhichNuc(matrix[pos(i,j,numChar)]));
                else if (ct == PROTEIN)
                    MrBayesPrint ("%c", WhichAA(matrix[pos(i,j,numChar)]));
                else if (ct == STANDARD)
                    MrBayesPrint ("%c", WhichStand(matrix[pos(i,j,numChar)]));
                else if (ct == RESTRICTION)
                    MrBayesPrint ("%c", WhichRes(matrix[pos(i,j,numChar)]));
                else if (ct == CONTINUOUS)
                    {
                    if (WhichCont(matrix[pos(i,j,numChar)]) < 0.0)
                        MrBayesPrint (" %2.2lf", WhichCont(matrix[pos(i,j,numChar)]));
                    else
                        MrBayesPrint ("  %2.2lf", WhichCont(matrix[pos(i,j,numChar)]));
                    }
                else
                    {
                    MrBayesPrint ("%s   Unknown data type\n", spacer);
                    return (ERROR);
                    }
                
                }
            MrBayesPrint ("\n");
            }
        MrBayesPrint ("\n");
        start = finish;
        } while (finish != numChar);

    return (NO_ERROR);
}


int DoShowUserTrees (void)
{
    int         i;

    if (numUserTrees == 0)
        {
        MrBayesPrint ("%s   No user trees have been defined\n", spacer);
        }
    else
        {
        for (i=0; i<numUserTrees; i++)
            {
            MrBayesPrint ("\n   Tree #%d -- '%s':\n\n", i+1, userTree[i]->name);
            ShowConTree (stdout, userTree[i], 70, NO);
            MrBayesPrint ("\n");
            }
        }

    return (NO_ERROR);
}


int DoShowBeagle (void)
{
#   if defined (BEAGLE_ENABLED)
    BeaglePrintResources();
#   else
    BeagleNotLinked();
#   endif
    return (NO_ERROR);
}


int DoTaxlabels (void)
{
    isTaxsetDef = YES;

    /* add default speciespartition name to list of valid speciespartitions */
    if (AddString (&speciespartitionNames, 0, "Default") == ERROR)
        {
        MrBayesPrint ("%s   Problem adding Default speciespartition to list\n", spacer);
        return (ERROR);
        }

    /* add default species name set */
    AddNameSet(&speciesNameSets, 0, taxaNames, numTaxa);

    /* set number of defined speciespartitions to 1 */
    numDefinedSpeciespartitions = 1;
        
    return (NO_ERROR);
}


int DoTaxlabelsParm (char *parmName, char *tkn)
{
    int         index;

    if (inTaxaBlock == NO)
        {
        MrBayesPrint ("%s   You must be in a taxa block to read a taxlabels command\n", spacer);
        return (ERROR);
        }

    if (defTaxa == NO)
        {
        MrBayesPrint ("%s   The number of taxa must be given before a set of taxon labels can be read\n", spacer);
        return ERROR;
        }

    if (isTaxsetDef == YES)
        {
        MrBayesPrint ("%s   A set of taxon labels has already been defined\n", spacer);
        if (defMatrix == NO)
            if (WantTo ("Do you want to delete the current set of taxon labels") == NO)
                return (SKIP_COMMAND);
            else
                FreeTaxa();
        else
            if (WantTo ("Do you want to delete the current character matrix") == NO)
                return (SKIP_COMMAND);
            else
                FreeMatrix();
        }

    if (expecting == Expecting(ALPHA) ||
        expecting == Expecting(NUMBER))
        {
        if (CheckString (taxaNames, numNamedTaxa, tkn, &index) == ERROR)
            {
            if (strlen(tkn)>99)
                {
                MrBayesPrint ("%s   Taxon name %s is too long. Maximun 99 characters is allowed.\n", spacer, tkn);
                return (ERROR);
                }
            if (AddString (&taxaNames, numNamedTaxa, tkn) == ERROR)
                {
                MrBayesPrint ("%s   Problem adding label %s to list of taxon labels\n", spacer, tkn);
                return (ERROR);
                }
            numNamedTaxa++;
            }
        else
            {
            MrBayesPrint ("%s   Taxon label '%s' is included twice in list of taxon labels\n", spacer, tkn);
            return (ERROR);
            }
        if (numNamedTaxa < numTaxa)
            {
            expecting = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            }
        else
            expecting |= Expecting(SEMICOLON);
        }

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
    MrBayesPrint ("%s", tkn);
}


int DoSpeciespartition (void)
{
    int     i, *partCount;
        
    /* add set to tempSet */
    if (fromI >= 0)
        if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
            {
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }

    /* set numDivisions; not set while reading the speciespartition */
    numDivisions = whichPartition + 1;
    
    /* check that all species are included */
    for (i=0; i<numTaxa; i++)
        {
        if (tempSet[i] == 0)
            {
            MrBayesPrint ("%s   Tip %d not included in speciespartition\n", spacer, i+1);
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        /*MrBayesPrint ("%4d %4d \n", i, tempSet[i]);*/
        }

    partCount = (int *) SafeCalloc (numDivisions, sizeof(int));
    if (!partCount)
        {
        for (i=0; i<numDivisions; i++)
            free(tempNames[i]);
        free (tempNames);
        tempNames = NULL;
        return ERROR;
        }

    /* make certain that the partition labels go from 1 - numTaxa, inclusive */
    for (i=0; i<numTaxa; i++)
        {
        if (tempSet[i] < 1 || tempSet[i] > numTaxa)
            {
            MrBayesPrint ("%s   Speciespartition index for tip %d out of bound (%d)\n", spacer, i+1, tempSet[i]);
            free (partCount);
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        partCount[tempSet[i] - 1]++;
        }
    for (i=0; i<numDivisions; i++)
        {
        if (partCount[i] == 0)
            {
            MrBayesPrint ("%s   Could not find a single tip for species %d\n", spacer, i+1);
            free (partCount);
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        }
    free (partCount);

    /* add name to list of valid partitions */
    if (AddString (&speciespartitionNames, numDefinedSpeciespartitions, tempSetName) == ERROR)
        {
        MrBayesPrint ("%s   Problem adding speciespartition %s to list\n", spacer, tempSetName);
        for (i=0; i<numDivisions; i++)
            free(tempNames[i]);
        free (tempNames);
        tempNames = NULL;
        return (ERROR);
        }

    /* add new partition */
    for (i=0; i<numTaxa; i++)
        {
        speciespartitionId[i] = (int *) SafeRealloc ((void *)(speciespartitionId[i]), ((size_t)numDefinedSpeciespartitions + 1) * sizeof(int));
        if (!speciespartitionId[i])
            {
            for (i=0; i<numDivisions; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return ERROR;
            }
        }

    /* set new partition */
    for (i=0; i<numTaxa; i++)
        speciespartitionId[i][numDefinedSpeciespartitions] = tempSet[i];

    /* add new set of species names */
    AddNameSet(&speciesNameSets, numDefinedSpeciespartitions, tempNames, numDivisions);

    /* free species names */
    for (i=0; i<numDivisions; i++)
        free(tempNames[i]);
    free (tempNames);
    tempNames = NULL;

    /* increment number of defined partitions */
    numDefinedSpeciespartitions++;
    
    return (NO_ERROR);
}


int DoSpeciespartitionParm (char *parmName, char *tkn)
{
    int             i, index, tempInt;
    
    if (defTaxa == NO || numTaxa == 0)
        {
        MrBayesPrint ("%s   A matrix or taxaset must be specified before partitions can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        /* set Speciespartition name ******************************************************************/
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check size of partition name */
            if (strlen(tkn) > 99)
                {
                MrBayesPrint ("%s   Partition name is too long. Max 100 characters\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a partition */
            if (numDefinedSpeciespartitions > 0)
                {
                if (CheckString (speciespartitionNames, numDefinedSpeciespartitions, tkn, &index) == ERROR)
                    {
                    /* if the partition name has not been used, then we should have an ERROR returned */
                    /* we _want_ to be here */

                    }
                else
                    {
                    MrBayesPrint ("%s   Speciespartition name '%s' has been used previously\n", spacer, tkn);
                    return (ERROR);
                    }
                }
                
            /* add the name temporarily to tempSetName */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 0;
    
            /* make sure tempNames is NULL */
            assert (tempNames == NULL);

            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            whichPartition = 0;
            foundFirst = NO;
            numDivisions = 0;
            MrBayesPrint ("%s   Defining speciespartition called '%s'\n", spacer, tkn);
            expecting = Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(ALPHA);
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (foundFirst == NO)
            {
            AddString(&tempNames, whichPartition, tkn);
            foundFirst = YES;
            expecting = Expecting(COLON);
            }
        else
            {
            /* We are defining a species partition in terms of a tip name (called tkn, here). We should be able
               to find tkn in the list of tip names. If we cannot, then we have a problem and
               return an error. */
            if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                {
                MrBayesPrint ("%s   Could not find a tip called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            /* add index of the tip named tkn to new tempSet */
            tempSet[index] = whichPartition + 1;
            fromI = toJ = everyK = -1;

            expecting  = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(SEMICOLON);
            expecting |= Expecting(COMMA);
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && tkn[0] == '.')
            tempInt = numTaxa;
        else
            sscanf (tkn, "%d", &tempInt);
        if (tempInt <= 0 || tempInt > numTaxa)
            {
            MrBayesPrint ("%s   Tip number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
            for (i=0; i<whichPartition; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted speciespartition\n", spacer);
                for (i=0; i<whichPartition; i++)
                    free(tempNames[i]);
                free (tempNames);
                tempNames = NULL;
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                MrBayesPrint ("%s   Improperly formatted speciespartition\n", spacer);
                for (i=0; i<whichPartition; i++)
                    free(tempNames[i]);
                free (tempNames);
                tempNames = NULL;
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted speciespartition\n", spacer);
                for (i=0; i<whichPartition; i++)
                    free(tempNames[i]);
                free (tempNames);
                tempNames = NULL;
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                    {
                    for (i=0; i<whichPartition; i++)
                        free(tempNames[i]);
                    free (tempNames);
                    tempNames = NULL;
                    return (ERROR);
                    }
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                MrBayesPrint ("%s   Improperly formatted speciespartition\n", spacer);
                    {
                    for (i=0; i<whichPartition; i++)
                        free(tempNames[i]);
                    free (tempNames);
                    tempNames = NULL;
                    return (ERROR);
                    }
                }
            }
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        expecting |= Expecting(COMMA);
        }
    else if (expecting == Expecting(COMMA))
        {
        /* add set to tempSet */
        if (fromI >= 0)
            if (AddToSet (fromI, toJ, everyK, whichPartition+1) == ERROR)
                {
                for (i=0; i<whichPartition; i++)
                    free(tempNames[i]);
                free (tempNames);
                tempNames = NULL;
                return (ERROR);
                }

        fromI = toJ = everyK = -1;
        foundDash = foundSlash = foundFirst = NO;
        whichPartition++;
        if (whichPartition > numTaxa)
            {
            MrBayesPrint ("%s   Too many speciespartitions (expecting maximum %d speciespartitions)\n", spacer, numTaxa);
            for (i=0; i<whichPartition; i++)
                free(tempNames[i]);
            free (tempNames);
            tempNames = NULL;
            return (ERROR);
            }
        expecting  = Expecting(ALPHA);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting  = Expecting(NUMBER);
        expecting |= Expecting(ALPHA);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        {
        for (i=0; i<whichPartition; i++)
            free(tempNames[i]);
        free (tempNames);
        tempNames = NULL;
        return (ERROR);
        }

    return (NO_ERROR);
}


int DoTaxaset (void)
{
    /* add set to tempSet */
    if (fromI >= 0 && toJ < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK < 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
    else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
        {
        if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
            return (ERROR);
        }
        
    /* add name to taxaSetNames */
    if (AddString (&taxaSetNames, numTaxaSets, tempSetName) == ERROR)
        {
        MrBayesPrint ("%s   Problem adding taxset %s to list\n", spacer, tempSetName);
        return (ERROR);
        }

    /* merge tempSet with taxaSet */
    AddBitfield (&taxaSet, numTaxaSets, tempSet, numTaxa);
    
    /* increment number of char sets */
    numTaxaSets++;
    
    /* show taxset (for debugging) */
#   if 0
    for (i=0; i<numTaxa; i++)
        MrBayesPrint ("%4d  %4d\n", i+1, tempSet[i]);
#   endif

    return (NO_ERROR);
}


int DoTaxasetParm (char *parmName, char *tkn)
{
    int     i, index, tempInt;
    
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A matrix must be specified before taxsets can be defined\n", spacer);
        return (ERROR);
        }

    if (expecting == Expecting(PARAMETER))
        {
        if (!strcmp(parmName, "Xxxxxxxxxx"))
            {
            /* check size of taxset name */
            if (strlen(tkn) > 99)
                {
                MrBayesPrint ("%s   Taxset name is too long\n", spacer);
                return (ERROR);
                }
                
            /* check to see if the name has already been used as a taxset */
            if (numTaxaSets > 0)
                {
                if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                    {
                    /* if the taxset name has not been used, then we should have an ERROR returned */
                    /* we _want_ to be here */

                    }
                else
                    {
                    MrBayesPrint ("%s   Taxset name has been used previously\n", spacer);
                    return (ERROR);
                    }
                }
            else if (numTaxaSets > 30)
                {
                MrBayesPrint ("%s   You cannot define more than 30 taxsets\n", spacer);
                return (ERROR);
                }
                
            /* add the name to the taxa set */
            strcpy (tempSetName, tkn);
            
            /* clear tempSet */
            for (i=0; i<numTaxa; i++)
                tempSet[i] = 0;
            
            fromI = toJ = everyK = -1;
            foundDash = foundSlash = NO;
            MrBayesPrint ("%s   Defining taxset called '%s'\n", spacer, tkn);
            expecting = Expecting(EQUALSIGN);
            }
        else
            return (ERROR);
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        /* We are defining a taxon set in terms of another (called tkn, here) or we are referring to
           the taxon name. We should be able to find tkn in the list of character set names or in the list
           of taxon names. If we cannot, then we have a problem and return an error. */
        if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
            {
            if (numTaxaSets < 1)
                {
                MrBayesPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            if (CheckString (taxaSetNames, numTaxaSets, tkn, &index) == ERROR)
                {
                MrBayesPrint ("%s   Could not find a taxset called '%s'\n", spacer, tkn);
                return (ERROR);
                }
            /* add taxa from taxset tkn to new tempSet */
            for (i=0; i<numTaxa; i++)
                {
                if (IsBitSet (i, taxaSet[index]) == YES)
                    tempSet[i] = 1;
                }
            }
        else
            {
            tempSet[index] = 1;
            }
        fromI = toJ = everyK = -1;

        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (strlen(tkn) == 1 && !strcmp(tkn, "."))
            {
            tempInt = numTaxa;
            }
        else
            {
            sscanf (tkn, "%d", &tempInt);
            if (tempInt <= 0 || tempInt > numTaxa)
                {
                MrBayesPrint ("%s   Taxon number %d is out of range (should be between %d and %d)\n", spacer, tempInt, 1, numTaxa);
                return (ERROR);
                }
            }
        tempInt--;
        if (foundDash == YES)
            {
            if (fromI >= 0)
                toJ = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted taxset\n", spacer);
                return (ERROR);
                }
            foundDash = NO;
            }
        else if (foundSlash == YES)
            {
            tempInt++;
            if (tempInt <= 1)
                {
                MrBayesPrint ("%s   Improperly formatted taxset\n", spacer);
                return (ERROR);
                }
            if (fromI >= 0 && toJ >= 0 && fromI < toJ)
                everyK = tempInt;
            else
                {
                MrBayesPrint ("%s   Improperly formatted taxset\n", spacer);
                return (ERROR);
                }
            foundSlash = NO;
            }
        else
            {
            if (fromI >= 0 && toJ < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                }
            else if (fromI < 0 && toJ < 0)
                {
                fromI = tempInt;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK < 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else if (fromI >= 0 && toJ >= 0 && everyK >= 0)
                {
                if (AddToSet (fromI, toJ, everyK, 1) == ERROR)
                    return (ERROR);
                fromI = tempInt;
                toJ = everyK = -1;
                }
            else
                {
                MrBayesPrint ("%s   Improperly formatted taxset\n", spacer);
                    {
                    return (ERROR);
                    }
                }
            }
        
        expecting  = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(SEMICOLON);
        expecting |= Expecting(DASH);
        expecting |= Expecting(BACKSLASH);
        }
    else if (expecting == Expecting(DASH))
        {
        foundDash = YES;
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(BACKSLASH))
        {
        foundSlash = YES;
        expecting = Expecting(NUMBER);
        }
    else
        return (ERROR);

    return (NO_ERROR);
}


int DoTaxaStat (void)
{
    int         i, j, maxLen, nameLen, nIncludedTaxa;
    char        tempName[100];
    
    if (defMatrix == NO)
        {
        MrBayesPrint ("%s   A character matrix must be defined first\n", spacer);
        return (ERROR);
        }
        
    /* find maximum length of taxon name */
    maxLen = nIncludedTaxa = 0;
    for (i=0; i<numTaxa; i++)
        {
        strcpy (tempName, taxaNames[i]);
        if ((int)strlen(tempName) > maxLen)
            maxLen = (int) strlen(tempName);
        if (taxaInfo[i].isDeleted == NO)
            nIncludedTaxa++;
        }
            
    MrBayesPrint ("%s   Showing taxon status:\n\n", spacer);
    if (nIncludedTaxa == numTaxa)
        MrBayesPrint ("%s     Number of taxa        = %d (all of which are included)\n", spacer, numTaxa);
    else
        MrBayesPrint ("%s     Number of taxa        = %d (of which %d are included)\n", spacer, numTaxa, nIncludedTaxa);
    MrBayesPrint ("%s     Number of constraints = %d\n\n", spacer, numDefinedConstraints);
    
    if (numDefinedConstraints > 0)
        {
        for (j=0; j<numDefinedConstraints; j++)
            {
            strcpy (tempName, constraintNames[j]);

            /* for now, ignore the probability */
            if (definedConstraintsType[j] == HARD)
                MrBayesPrint ("%s     %2d -- Trees with 'hard' constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
            else if (definedConstraintsType[j] == PARTIAL)
                MrBayesPrint ("%s     %2d -- Trees with 'partial' constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
            else
                MrBayesPrint ("%s     %2d -- Trees with 'negative' constraint \"%s\" are infinitely\n", spacer, j+1, tempName);
            MrBayesPrint ("%s           more probable than those without \n", spacer);
            }
        MrBayesPrint ("\n");
        for (j=0; j<maxLen; j++)
            MrBayesPrint (" ");
        MrBayesPrint ("                             Constraints\n");
        }
    MrBayesPrint ("%s     Taxon  ", spacer);
    for (j=0; j<maxLen; j++)
        MrBayesPrint (" ");
    MrBayesPrint ("   Inclusion");
    MrBayesPrint ("   ");
    for (j=0; j<numDefinedConstraints; j++)
        MrBayesPrint (" %2d", j+1);
    MrBayesPrint ("\n");
    MrBayesPrint ("%s   -------", spacer);
    for (j=0; j<maxLen; j++)
        MrBayesPrint ("-");
    MrBayesPrint ("--------------");
    
    if (numDefinedConstraints > 0)
        {
        MrBayesPrint ("----");
        for (j=0; j<numDefinedConstraints; j++)
            MrBayesPrint ("---");
        }
    MrBayesPrint ("\n");
    for (i=0; i<numTaxa; i++)
        {
        strcpy (tempName, taxaNames[i]);
        nameLen = (int) strlen(tempName);
        
        if (i == outGroupNum)
            MrBayesPrint ("%s ->%4d (%s) ", spacer, i+1, tempName);
        else
            MrBayesPrint ("%s   %4d (%s) ", spacer, i+1, tempName);
        for (j=0; j<(maxLen-nameLen); j++)
            MrBayesPrint (" ");
        MrBayesPrint (" -- ");
        
        if (taxaInfo[i].isDeleted == YES)
            MrBayesPrint ("Deleted ");
        else
            MrBayesPrint ("Included");
            
        MrBayesPrint ("    ");
            
        for (j=0; j<numDefinedConstraints; j++)
            {
            if (definedConstraintsType[j] == HARD)
                {
                if (IsBitSet(i, definedConstraint[j]) == NO)
                    MrBayesPrint ("  .");
                else
                    MrBayesPrint ("  *");
                }
            else if (definedConstraintsType[j] == PARTIAL)
                {
                if (IsBitSet(i, definedConstraint[j]) == YES)
                    MrBayesPrint ("  +");
                else if (IsBitSet(i, definedConstraintTwo[j]) == YES)
                    MrBayesPrint ("  -");
                else
                    MrBayesPrint ("  .");
                }
            else if (definedConstraintsType[j] == NEGATIVE)
                {
                if (IsBitSet(i, definedConstraint[j]) == NO)
                    MrBayesPrint ("  .");
                else
                    MrBayesPrint ("  #");
                }
            }
        MrBayesPrint ("\n");
        }
        
    MrBayesPrint ("\n");
    MrBayesPrint ("%s   '.' indicate that the taxon is not present in the constraint. \n", spacer);
    MrBayesPrint ("%s   '*' indicate that the taxon is present in the 'hard' constraint. \n", spacer);
    MrBayesPrint ("%s   '+' indicate that the taxon is present in the first groupe of 'partial' constraint. \n", spacer);
    MrBayesPrint ("%s   '-' indicate that the taxon is present in the second groupe of 'partial' constraint. \n", spacer);
    MrBayesPrint ("%s   '#' indicate that the taxon is present in the 'negative' constraint. \n", spacer);
    MrBayesPrint ("%s   Arrow indicates current outgroup. \n", spacer);

    return (NO_ERROR);
}


int DoTranslate (void)
{
    int     i, j;

    if (inTreesBlock == NO)
        {
        MrBayesPrint ("%s   You must be in a trees block to read a translate command\n", spacer);
        return (ERROR);
        }
    numTranslates++;    /* number of taxa in translate table */
    isTranslateDef = YES;

    isTranslateDiff = NO;
    if (isTaxsetDef == NO)
        SetTaxaFromTranslateTable();
    else
        {
        for (i=0; i<numTranslates; i++)
            {
            strcpy (token, transFrom[i]);
            if (CheckString (taxaNames, numTaxa, token, &j) == ERROR)
                {
                isTranslateDiff = YES;
                }
            }
        if (numTranslates != numTaxa)
            isTranslateDiff = YES;
        }

    return (NO_ERROR);
}


int DoTranslateParm (char *parmName, char *tkn)
{
    int         index;
    static int  whichTranslate;

    if (inTreesBlock == NO)
        {
        MrBayesPrint ("%s   You must be in a trees block to read a translate command\n", spacer);
        return (ERROR);
        }

    if (isTranslateDef == YES)
        {
        MrBayesPrint ("%s   A translation has already been defined for this tree block\n", spacer);
        return (ERROR);
        }
        
    if (expecting == Expecting(ALPHA) ||
        expecting == Expecting(NUMBER))
        {
        if (numTaxa == 0)
            {
            MrBayesPrint ("%s   Data matrix should be defined before translation table could be set.\n", spacer);
            return (ERROR);
            }
        if (numTranslates == numTaxa)
            {
            MrBayesPrint ("%s   Too many entries in translation table. Maximum number of taxon names to translate is %d\n", spacer,numTaxa);
            return (ERROR);
            }
        if (whichTranslate == 0)
            {
            if (CheckString (transTo, numTranslates, tkn, &index) == ERROR)
                {
                if (AddString (&transTo, numTranslates, tkn) == ERROR)
                    {
                    MrBayesPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
                    return (ERROR);
                    }
                }
            else
                {
                MrBayesPrint ("%s   Already found name (%s) in list\n", spacer, tkn);
                return (ERROR);
                }           
            whichTranslate++;
            expecting = Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            }
        else 
            {
            if (CheckString (transFrom, numTranslates, tkn, &index) == ERROR)
                {
                if (AddString (&transFrom, numTranslates, tkn) == ERROR)
                    {
                    MrBayesPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
                    return (ERROR);
                    }
                }
            else
                {
                MrBayesPrint ("%s   Already found name (%s) in list\n", spacer, tkn);
                return (ERROR);
                }           
            whichTranslate = 0;
            expecting = Expecting(COMMA);
            expecting |= Expecting(SEMICOLON);
            }
        }
    else if (expecting == Expecting(COMMA))
        {
        numTranslates++;
        expecting = Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
    MrBayesPrint ("%s", tkn);
}


int DoTree (void)
{
    readComment = NO;

    if (inSumtCommand == YES || inComparetreeCommand == YES)
        return (DoSumtTree ());

    return (NO_ERROR);
}


int DoTreeParm (char *parmName, char *tkn)
{
    int                 i, tempInt, index;
    MrBFlt              tempD;
    char                tempName[100];
    static BitsLong     lastExpecting; /* keep track of what we expected before a comment, in case we want to skip a comment */
    static char         *tempNameString=NULL; /* Contains multiple tokens which form name string of param set*/
    static int          foundAmpersand, foundColon, foundComment, foundE, foundB, foundN, foundFirst,
                        foundCurly, /* is set to YES when we are between two curly bracets ONLY while processing CppEvent name */
                        foundClockrate, 
                        foundName, /*is set to YES when param set name token is found and set to NO once full param set name is processed*/
                        eSetIndex, /* is set in the begining of reading CppEvent for a node/branch to the index of currently processed CppEvent set */
                        bSetIndex, eventIndex, treeIndex, nextIntNodeIndex;
    static PolyNode     *pp, *qq;
    static PolyTree     *t;
    
    /* This function will read in components of a tree description. We expect one of the following formats:
    
          tree <name> = [&R] <newick-description>;
          tree <name> = [&U] <newick-description>;
          tree <name> [&E CppEvents]  = [&R] [&clockrate = 1.23] ((1:0.021[&E CppEvents 2: (0.10 1.11,0.83 3.17)],...
          tree <name> [&B TK02Brlens] = [&R] [&clockrate = 1.23] ((1:0.021[&B TK02Brlens 0.019],...
          tree <name> [&B IgrBrlens]  = [&R] [&clockrate = 1.23] ((1:0.021[&B IgrBrlens 0.019],...
     
       Values will be stored in event sets that go with the tree and that are used to initialize the relaxed clock
       parameters before a run is started. Note that several sets of events can be stored with each tree.
    */

    if (isTaxsetDef == NO)
        {
        MrBayesPrint ("%s   Taxon labels must be specified before a tree could be red in\n", spacer);
        return (ERROR);
        }
    if (inTreesBlock == NO)
        {
        MrBayesPrint ("%s   You must be in a trees block to read a tree\n", spacer);
        return (ERROR);
        }
    
    if (expecting == Expecting(PARAMETER))
        {
        /* this is the name of the tree */
        if (inSumtCommand==YES || inComparetreeCommand == YES)
            {
            /* we are reading in a tree to sumt or comparetree counters */
            t = sumtParams.tree;
            ResetPolyTree (t);
            }
        else
            {
            /* we are reading in a user tree */
            /* check if the tree exists */
            treeIndex = 0;
            for (i=0; i<numUserTrees; i++)
                if (strcmp(tkn,userTree[i]->name) == 0)
                    break;
            treeIndex = i;
            if (treeIndex < numUserTrees)
                {
                MrBayesPrint ("%s   Overwriting tree '%s'.\n", spacer, userTree[treeIndex]);
                FreePolyTree (userTree[treeIndex]);
                }
            if ((userTree[treeIndex] = AllocatePolyTree (numTaxa)) == NULL)
                return (ERROR);
            t = userTree[treeIndex];
            }
        strncpy (t->name, tkn, 99);
        foundColon = foundAmpersand = foundEqual = foundComment = NO;
        foundE = foundB = foundN = foundFirst = foundClockrate = foundName = NO;
        eSetIndex = bSetIndex = eventIndex = 0;
        nextAvailableNode = 0;
        if (isTranslateDef == YES && isTranslateDiff == YES)
            nextIntNodeIndex = numTranslates;
        else
            nextIntNodeIndex = numTaxa;
        pp = &t->nodes[nextAvailableNode++];
        t->root = pp;
        t->isRooted = NO;  /* expect unrooted tree */
        t->isClock = NO;   /* expect nonclock tree */
        t->isCalibrated = NO;  /* expect uncalibrated tree */
        t->isRelaxed = NO;    /* expect strict clock if clock tree */
        t->clockRate = 0.0;     /* expect no clock rate */
        t->popSizeSet = NO;     
        readComment = YES;
        expecting = Expecting(EQUALSIGN) | Expecting(LEFTCOMMENT);
        lastExpecting = expecting;
        }
    else if (expecting == Expecting(EQUALSIGN))
        {
        if (foundClockrate == YES)
            expecting = Expecting(NUMBER);
        else
            {
            for (i=0; i<numTaxa; i++)
                tempSet[i] = NO;
            foundEqual = YES;
            expecting = Expecting(LEFTPAR) | Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        if (foundE == YES)
            {
            expecting = Expecting(NUMBER);
            }
        else
            {
            if (nextAvailableNode >= 2*numTaxa)
                {
                MrBayesPrint ("%s   Too many nodes on tree '%s'\n", spacer, t->name);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            qq = &t->nodes[nextAvailableNode++];
            qq->anc = pp;
            pp->left = qq;
            pp->index = nextIntNodeIndex++;
            pp = qq;
            expecting = Expecting(LEFTPAR);
            expecting |= Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        }
    else if (expecting == Expecting(ALPHA))
        {
        if (foundAmpersand == YES)
            {
            if (strcmp(tkn,"E") == 0)
                {
                foundE = YES;
                expecting = Expecting(ALPHA);
                }
            else if (strcmp(tkn,"B") == 0)
                {
                foundB = YES;
                expecting = Expecting(ALPHA);
                }
            else if (strcmp(tkn,"N") == 0)
                {
                foundN = YES;
                expecting = Expecting(ALPHA);
                }
            else if (strcmp(tkn, "R") == 0)
                {
                t->isRooted = YES;
                t->isClock = YES;   /* assume clock if rooted */
                expecting = Expecting(RIGHTCOMMENT);
                }
            else if (strcmp(tkn, "U") == 0)
                {
                t->isRooted = NO;
                expecting = Expecting(RIGHTCOMMENT);
                }
            else if (strcmp(tkn, "clockrate") == 0)
                {
                t->isCalibrated = YES;
                foundClockrate = YES;
                expecting = Expecting(EQUALSIGN);
                }
            else
                {
                inComment = YES;
                numComments++;
                expecting = lastExpecting;
                }
            foundAmpersand = NO;
            }
        else if (foundName == YES && foundCurly == YES)
            {
            if (strcmp("all",tkn) == 0)
                {
                SafeStrcat (&tempNameString,tkn);
                expecting = Expecting(RIGHTCURL);
                }
            else
                {
                MrBayesPrint ("%s   Urecognized argument '%s'\n", spacer, tkn);
                return (ERROR);
                }
            }
        else if (foundE == YES) /* We have seen &E */
            {
            if (foundEqual == NO) /* We have not seen name before and we are in header */
                {
                t->nESets++;
                t->isRelaxed = YES;
                t->nEvents  = (int **) SafeRealloc ((void *)t->nEvents, t->nESets*sizeof(int *));
                t->position = (MrBFlt ***) SafeRealloc ((void *)t->position, t->nESets*sizeof(MrBFlt **));
                t->rateMult = (MrBFlt ***) SafeRealloc ((void *)t->rateMult, t->nESets*sizeof(MrBFlt **));
                t->nEvents[t->nESets-1]  = (int *) SafeCalloc (2*(size_t)numTaxa, sizeof(int));
                t->position[t->nESets-1] = (MrBFlt **) SafeCalloc (2*(size_t)numTaxa, sizeof(MrBFlt *));
                t->rateMult[t->nESets-1] = (MrBFlt **) SafeCalloc (2*(size_t)numTaxa, sizeof(MrBFlt *));
                t->eSetName = (char **) SafeRealloc ((void *)t->eSetName, t->nESets*sizeof(char **));
                }
            SafeStrcpy (&tempNameString,tkn);
            foundName = YES;
            expecting = Expecting(LEFTCURL);
            if (foundEqual == YES)
                expecting |= Expecting(NUMBER);
            else
                expecting |= Expecting(RIGHTCOMMENT);
            }
        else if (foundB == YES)
            {
            if (foundEqual == NO)
                {
                t->nBSets++;
                t->isRelaxed = YES;
                t->effectiveBrLen = (MrBFlt **) SafeRealloc ((void *)t->effectiveBrLen, (size_t)(t->nBSets)*sizeof(MrBFlt *));
                t->effectiveBrLen[t->nBSets-1] = (MrBFlt *) SafeCalloc (2*(size_t)numTaxa, sizeof(MrBFlt));
                for (i=0; i<2*numTaxa; i++)
                    t->effectiveBrLen[t->nBSets-1][i] = 1.0;
                t->bSetName = (char **) SafeRealloc ((void *)t->bSetName, (size_t)(t->nBSets)*sizeof(char *));
                t->bSetName[t->nBSets-1] = (char *) SafeCalloc (strlen(tkn)+1, sizeof(char));
                }
            SafeStrcpy (&tempNameString,tkn);
            foundName = YES;
            expecting = Expecting(LEFTCURL);
            if (foundEqual == YES)
                expecting |= Expecting(NUMBER);
            else
                expecting |= Expecting(RIGHTCOMMENT);
            }
        else if (foundN == YES)
            {
            if (foundEqual == NO)
                {
                if (t->popSizeSet == YES)
                    {
                    MrBayesPrint ("%s   Cannot hold more than one population size set\n", spacer);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                t->popSizeSet = YES;
                if (isTranslateDef == YES && isTranslateDiff == YES)
                    t->popSize = (MrBFlt *) SafeCalloc (2*numTranslates, sizeof(MrBFlt));
                else
                    t->popSize = (MrBFlt *) SafeCalloc (2*numLocalTaxa, sizeof(MrBFlt));
                }
            SafeStrcpy (&tempNameString,tkn);
            foundName = YES;
            expecting = Expecting(LEFTCURL);
            if (foundEqual == YES)
                expecting |= Expecting(NUMBER);
            else
                expecting |= Expecting(RIGHTCOMMENT);
            }
        else   /* taxon name */
            {
            if (isTranslateDef == YES)
                {
                /* we are using the translation table */
                if (CheckString (transTo, numTranslates, tkn, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find token '%s' in taxon translation table\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                strcpy (tempName, transFrom[index]);
                if (isTranslateDiff == NO && CheckString (taxaNames, numTaxa, tempName, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find taxon '%s' in list of taxa\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                if (tempSet[index] == YES)
                    {
                    MrBayesPrint ("%s   Taxon name '%s' already used in tree\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                tempSet[index] = YES;
                strcpy (pp->label, tempName);
                pp->index = index;
                }
            else
                {
                /* Check to see if the name is in the list of taxon names. */
                if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find taxon '%s' in list of taxa\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                if (tempSet[index] == YES)
                    {
                    MrBayesPrint ("%s   Taxon name '%s' already used in tree\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                tempSet[index] = YES;
                strcpy (pp->label, tkn);
                pp->index = index;
                }
            expecting  = Expecting(COMMA);
            expecting |= Expecting(COLON);
            expecting |= Expecting(RIGHTPAR);
            }
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        if (foundE == YES)
            expecting = Expecting(RIGHTCOMMENT);
        else
            {
            if (pp->anc == NULL)
                {
                MrBayesPrint ("%s   Incorrect tree format: cannot go down\n", spacer);//, tkn
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            if (pp->anc->left == pp)
                {
                MrBayesPrint ("%s   Incorrect tree format: all nodes except tips should have more then one child. Either a single\n", spacer);
                MrBayesPrint ("%s   taxon is surrounded with brackets or there is a clade surrounded by double brackets.\n", spacer);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            pp = pp->anc;
            if (pp->anc == NULL)
                {
                /* finish up tree */
                t->nNodes = nextAvailableNode;
                t->nIntNodes = t->nNodes;
                for (i=0; i<t->nNodes; i++)
                    {
                    if (t->nodes[i].left == NULL)
                        t->nIntNodes--;
                    }
                GetPolyDownPass(t);
                
                /* check that number of taxa is correct */
                if (t->isRooted == NO && t->nNodes-t->nIntNodes == t->nIntNodes + 1)
                    t->isRooted = YES;
                if ((t->isRooted == YES && t->nNodes-t->nIntNodes != t->nIntNodes + 1) ||
                    (t->isRooted == NO  && t->nNodes-t->nIntNodes != t->nIntNodes + 2))
                    {
                    /* we are protected from adding too many taxa by taxon-matching code above */
                    if (t->isRooted == YES && t->nNodes-t->nIntNodes == t->nIntNodes + 2)
                        {
                        MrBayesPrint ("%s   The tree is declared as rooted (by comment [&R]) but\n", spacer);
                        MrBayesPrint ("%s   the given tree has unrooted structure.\n", spacer);
                        }
                    else
                        MrBayesPrint ("%s   Taxa missing in tree, or NOT a binary tree\n", spacer);

                    return (ERROR);
                    }

                /* check other properties */
                if (t->isClock == YES && t->isRooted == NO)
                    {
                    MrBayesPrint ("%s   Tree has clock rate but is not rooted\n", spacer);
                    return (ERROR);
                    /* Note: any deviation from an ultrametric tree must be assumed to be due to dated
                       tips at this point */
                    }
                if (t->isRelaxed == YES && t->isClock == NO)
                    {
                    MrBayesPrint ("%s   Tree has relaxed clock rates but is not a clock tree\n", spacer);
                    return (ERROR);
                    }
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    {
                    if (treeIndex == numUserTrees)
                        numUserTrees++;
                    MrBayesPrint ("%s   Successfully read tree '%s'\n", spacer, userTree[treeIndex]->name);
                    }
                if (t->popSize == NULL)
                    {
                    readComment = NO;
                    expecting = Expecting(SEMICOLON);
                    }
                else
                    {
                    readComment = YES;
                    expecting = Expecting(LEFTCOMMENT);
                    lastExpecting = expecting;
                    }
                }
            else
                {
                expecting = Expecting(COMMA);
                expecting |= Expecting(COLON);
                expecting |= Expecting(RIGHTPAR);
                }
            }
        }
    else if (expecting == Expecting(COLON))
        {
        foundColon = YES;
        if (foundE == YES)
            expecting = Expecting(LEFTPAR);
        else
            expecting  = Expecting(NUMBER);
        expecting |= Expecting(LEFTCOMMENT);
        lastExpecting = expecting;
        }
    else if (expecting == Expecting(COMMA))
        {
        if (foundName == YES)
            {
            SafeStrcat (&tempNameString,",");
            expecting = Expecting(NUMBER);
            }
        else if (foundE == YES)
            {
            expecting = Expecting(NUMBER);
            }
        else
            {
            if (nextAvailableNode >= 2*numTaxa)
                {
                MrBayesPrint ("%s   Too many nodes on tree '%s'\n", spacer, t->name);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            qq = &t->nodes[nextAvailableNode++];
            pp->sib = qq;
            qq->anc = pp->anc;
            pp = qq;
            expecting = Expecting(LEFTPAR);
            expecting |= Expecting(ALPHA);
            expecting |= Expecting(NUMBER);
            expecting |= Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        }
    else if (expecting == Expecting(NUMBER))
        {
        if (foundClockrate == YES)
            {
            sscanf (tkn, "%lf", &tempD);
            t->clockRate = tempD;
            foundClockrate = NO;
            expecting = Expecting(RIGHTCOMMENT);
            }
        else if (foundName == YES && foundCurly == YES)
            {
            /* still assembling name of a param set */
            SafeStrcat (&tempNameString,tkn);       
            expecting = Expecting(RIGHTCURL) | Expecting(COMMA);
            }
        else if (foundN == YES)
            {
            /* we only know now that name is complete if it does not have curlies in it */
            foundName = NO;

            if (strcmp(tempNameString,t->popSizeSetName) != 0)
                {
                MrBayesPrint ("%s   Could not find population size set '%s'\n", spacer, tempNameString);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }

            sscanf (tkn, "%lf", &tempD);
            t->popSize[pp->index] = tempD;
            foundN = NO;
            expecting = Expecting(RIGHTCOMMENT);
            }
        else if (foundB == YES)
            {
            /* we only know now that name is complete if it does not have curlies in it */
            foundName = NO;

            /* find the right effective branch length set */
            for (i=0; i<t->nBSets; i++)
                if (strcmp(t->bSetName[i],tempNameString) == 0)
                    break;
            if (i == t->nBSets)
                {
                MrBayesPrint ("%s   Could not find effective branch length set '%s'\n", spacer, tempNameString);
                if (inSumtCommand == NO && inComparetreeCommand == NO)
                    FreePolyTree (userTree[treeIndex]);
                return (ERROR);
                }
            bSetIndex = i;

            sscanf (tkn, "%lf", &tempD);
            t->effectiveBrLen[bSetIndex][pp->index] = tempD;
            foundB = NO;
            expecting = Expecting(RIGHTCOMMENT);
            }
        else if (foundE == YES)
            {
            if (foundColon == NO)
                {
                /* we only know now that name is complete if it does not have curlies in it */
                foundName = NO;

                /* find the right event set */
                for (i=0; i<t->nESets; i++)
                    if (strcmp(t->eSetName[i],tempNameString) == 0)
                        break;
                if (i == t->nESets)
                    {
                    MrBayesPrint ("%s   Could not find event set '%s'\n", spacer, tempNameString);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                eSetIndex = i;

                sscanf (tkn, "%d", &tempInt);
                if (tempInt < 0)
                    {
                    MrBayesPrint ("%s   Wrong number of events (%d) for event set '%s'\n", spacer, tempInt, t->eSetName[eSetIndex]);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                t->nEvents[eSetIndex][pp->index]  = tempInt;
                if (tempInt > 0)
                    {
                    t->position[eSetIndex][pp->index] = (MrBFlt *) SafeCalloc (tempInt, sizeof(MrBFlt));
                    t->rateMult[eSetIndex][pp->index] = (MrBFlt *) SafeCalloc (tempInt, sizeof(MrBFlt));
                    expecting = Expecting (COLON);
                    if (inSumtCommand == YES || inComparetreeCommand == YES)
                        expecting |= Expecting (RIGHTCOMMENT);  /* we allow empty event specifications in sumt and comparetree */
                    }
                else
                    expecting = Expecting (RIGHTCOMMENT);
                eventIndex = 0;
                }
            else if (foundFirst == NO)
                {
                /* processing the first number in the cpp event pair <position rate> */
                sscanf (tkn, "%lf", &tempD);
                t->position[eSetIndex][pp->index][eventIndex] = tempD;
                expecting = Expecting(NUMBER);
                foundFirst = YES;
                }
            else
                {
                /* processing the second number in the cpp event pair <position rate> */
                foundFirst = NO;
                sscanf (tkn, "%lf", &tempD);
                t->rateMult[eSetIndex][pp->index][eventIndex] = tempD;
                eventIndex++;
                if (eventIndex == t->nEvents[eSetIndex][pp->index])
                    {
                    expecting = Expecting(RIGHTPAR);
                    foundColon = NO;
                    }
                else
                    expecting = Expecting(COMMA);
                }
            }
        else if (foundColon == YES)
            {
            /* branch length */
            sscanf (tkn, "%lf", &tempD);
            pp->length = tempD;
            foundColon = NO;
            t->brlensDef = YES;
            expecting  = Expecting(COMMA);
            expecting |= Expecting(RIGHTPAR);
            expecting |= Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        else    /* taxon identifier */
            {
            if (isTranslateDef == YES)
                {
                /* we are using the translation table */
                if (CheckString (transTo, numTranslates, tkn, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find token '%s' in taxon translation table\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                strcpy (tempName, transFrom[index]);
                if (isTranslateDiff == NO && CheckString (taxaNames, numTaxa, tempName, &index) == ERROR)
                    {
                    MrBayesPrint ("%s   Could not find taxon '%s' in list of taxa\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                if (tempSet[index] == YES)
                    {
                    MrBayesPrint ("%s   Taxon name '%s' already used in tree\n", spacer, tkn);
                    if (inSumtCommand == NO && inComparetreeCommand == NO)
                        FreePolyTree (userTree[treeIndex]);
                    return (ERROR);
                    }
                tempSet[index] = YES;
                strcpy (pp->label, tempName);
                pp->index = index;
                }
            else
                {
                /* Simply use taxon number; first check to see if the name is in the list of taxon names. */
                if (CheckString (taxaNames, numTaxa, tkn, &index) == ERROR)
                    {
                    /* The number could not be found as a taxon name in the list of taxon names. We will
                        assume that the user has then input taxa as numbers and not the names. */
                    sscanf (tkn, "%d", &index);
                    if (index < 1 || index > numTaxa)
                        {
                        MrBayesPrint ("%s   Taxon number %d is out of range\n", spacer, index);
                        if (inSumtCommand == NO && inComparetreeCommand == NO)
                            FreePolyTree (userTree[treeIndex]);
                        return (ERROR);
                        }
                    index--;
                    if (tempSet[index] == YES)
                        {
                        MrBayesPrint ("%s   Taxon name %d has already been used in tree '%s'\n", spacer, index+1, t->name);
                        if (inSumtCommand == NO && inComparetreeCommand == NO)
                            FreePolyTree (userTree[treeIndex]);
                        return (ERROR);
                        }
                    }
                else
                    {
                    /* The number is in the list of taxon names */
                    if (index < 0 || index >= numTaxa)
                        {
                        MrBayesPrint ("%s   Taxon name %s could not be found\n", spacer, tkn);
                        if (inSumtCommand == NO && inComparetreeCommand == NO)
                            FreePolyTree (userTree[treeIndex]);
                        return (ERROR);
                        }
                    if (tempSet[index] == YES)
                        {
                        MrBayesPrint ("%s   Taxon %d has already been used in tree '%s'\n", spacer, index+1, t->name);
                        if (inSumtCommand == NO && inComparetreeCommand == NO)
                            FreePolyTree (userTree[treeIndex]);
                        return (ERROR);
                        }
                    }
                tempSet[index] = YES;
                strcpy (pp->label, taxaNames[index]);
                pp->index = index;
                }
            expecting  = Expecting(COMMA);
            expecting |= Expecting(COLON);
            expecting |= Expecting(RIGHTPAR);
            expecting |= Expecting(LEFTCOMMENT);
            lastExpecting = expecting;
            }
        }
    else if (expecting == Expecting(LEFTCOMMENT))
        {
        expecting = Expecting(AMPERSAND);
        foundComment = YES;
        }
    else if (expecting == Expecting(RIGHTCOMMENT))
        {
        if (foundEqual == NO)
            {
            /* We may have a complete name of a set of branch parameters, which needs to be recorded */
            if (foundName == YES)
                {
                if (foundE == YES)
                    {
                    t->eSetName[t->nESets-1] = (char *) SafeCalloc (strlen(tempNameString)+1,sizeof(char));
                    strcat(t->eSetName[t->nESets-1],tempNameString);
                    }
                else if (foundB == YES)
                    {
                    t->bSetName[t->nBSets-1] = (char *) SafeCalloc (strlen(tempNameString)+1,sizeof(char));
                    strcat(t->bSetName[t->nBSets-1],tempNameString);
                    }
                else if (foundN == YES)
                    {
                    t->popSizeSetName = (char *) SafeCalloc (strlen(tempNameString)+1,sizeof(char));
                    strcpy(t->popSizeSetName,tempNameString);
                    }
                foundName = NO;
                }
            expecting = Expecting(EQUALSIGN);
            }
        else
            {
            if (pp->anc == NULL)
                {
                if (pp->left == NULL)
                    expecting = Expecting(LEFTPAR);
                else
                    expecting = Expecting(SEMICOLON);
                }
            else if (pp == pp->anc->left)
                expecting = Expecting(COMMA);
            else
                expecting = Expecting(RIGHTPAR);
            }
        foundE = foundB = foundN = NO;
        expecting |= Expecting(LEFTCOMMENT);
        }
    else if (expecting == Expecting(AMPERSAND))
        {
        foundAmpersand = YES;
        foundComment = NO;
        expecting = Expecting (ALPHA);
        }
    else if (foundComment == YES)
        {
        numComments++;
        foundComment = NO;
        }
    else if (expecting == Expecting(LEFTCURL))
        {
        if (foundName == YES)
            {
            foundCurly=YES;
            SafeStrcat (&tempNameString,"{");               
            expecting = Expecting(NUMBER) | Expecting(ALPHA);
            }
        else
            return(ERROR);
        }
    else if (expecting == Expecting(RIGHTCURL))
        {
        if (foundName == YES)
            {
            SafeStrcat (&tempNameString,"}");
            foundCurly=NO;
            if (foundEqual == NO)
                {
                /* We are processing a name of a set of branch params in the header of a tree.  */
                expecting = Expecting(RIGHTCOMMENT);
                }
            else
                {
                /* We are processing a param value of a branch param set  */
                expecting = Expecting(NUMBER);
                }
            }
        else
            return(ERROR);
        }

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
}


int DoUserTree (void)
{
    MrBayesPrint ("%s   Usertree command deprecated. Define the tree in a treeblock and use 'Startvals' instead.\n", spacer);
    return (ERROR);
}


int DoUserTreeParm (char *parmName, char *tkn)
{
    if (expecting == Expecting(EQUALSIGN))
        {
        expecting = Expecting(LEFTPAR);
        expecting |= Expecting(RIGHTPAR);
        expecting |= Expecting(COLON);
        expecting |= Expecting(NUMBER);
        expecting |= Expecting(ALPHA);
        expecting |= Expecting(SEMICOLON);
        }
    else if (expecting == Expecting(LEFTPAR))
        {
        expecting = Expecting(LEFTPAR);
        expecting |= Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(ALPHA))
        {
        expecting = Expecting(COLON);
        expecting |= Expecting(COMMA);
        expecting |= Expecting(RIGHTPAR);
        }
    else if (expecting == Expecting(NUMBER))
        {
        expecting = Expecting(COLON);
        expecting |= Expecting(COMMA);
        expecting |= Expecting(RIGHTPAR);
        }
    else if (expecting == Expecting(COLON))
        {
        expecting = Expecting(NUMBER);
        }
    else if (expecting == Expecting(COMMA))
        {
        expecting = Expecting(LEFTPAR);
        expecting |= Expecting(ALPHA);
        expecting |= Expecting(NUMBER);
        }
    else if (expecting == Expecting(RIGHTPAR))
        {
        expecting = Expecting(RIGHTPAR);
        expecting |= Expecting(COMMA);
        expecting |= Expecting(COLON);
        expecting |= Expecting(SEMICOLON);
        }
    else
        return (ERROR);

    return (NO_ERROR);
    MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */
    MrBayesPrint ("%s", tkn);
}

int DoVersion (void)
{
    MrBayesPrint ("   ---------------------------------------------------------------------------\n");
    MrBayesPrint ("   Version\n");
    MrBayesPrint ("\n");
    MrBayesPrint ("   MrBayes v%s\n", VERSION_NUMBER);
    MrBayesPrint("\n");
    MrBayesPrint("   Features: ");
#ifdef SSE_ENABLED
    MrBayesPrint(" SSE");
#endif
#ifdef AVX_ENABLED
    MrBayesPrint(" AVX");
#endif
#ifdef FMA_ENABLED
    MrBayesPrint(" FMA");
#endif
#ifdef BEAGLE_ENABLED
    MrBayesPrint(" Beagle");
#endif
#ifdef MPI_ENABLED
    MrBayesPrint(" MPI");
#endif
#ifdef HAVE_LIBREADLINE
    MrBayesPrint(" readline");
#endif
    MrBayesPrint("\n");
#if defined(HOST_TYPE) && defined(HOST_CPU)
    MrBayesPrint("   Host type: %s (CPU: %s)\n", HOST_TYPE, HOST_CPU);
#endif
#if defined(COMPILER_VENDOR) && defined(COMPILER_VERSION)
    MrBayesPrint("   Compiler:  %s %s\n", COMPILER_VENDOR, COMPILER_VERSION);
#endif
    MrBayesPrint ("   ---------------------------------------------------------------------------\n");

    return (NO_ERROR);
}


BitsLong Expecting (int y)
{
    BitsLong x;
    
    x = (BitsLong) pow (2.0, (MrBFlt)y);
    
    return (x);
}


#ifdef HAVE_LIBREADLINE
/* This function is for commandline substitution: first word is always a command */
char *command_generator (const char *text, int state)
{
    static int      list_index, len;
    char           *command;
    char           *dupstring;

    if (state == 0)
        {
        list_index = 0;
        len = (int) strlen (text);
        }

    while ((command = commands[list_index].string) != NULL)
        {
        list_index++;

        if (StrCmpCaseInsensitiveLen (command, text, len) == 0)
            {
            /* memory is freed by the readline library so we need a strdup here */
            dupstring = SafeMalloc (strlen (command) + 1);
            strcpy (dupstring, command);
            return dupstring;
            }
        }

    return NULL;
}
#endif


int FindValidCommand (char *tk, int *numMatches)
{
    int             i, j, tkLen, targetLen, numDiff;
    CmdType         *p;

    p = commands + 0;
    tkLen = (int) strlen(tk);

    (*numMatches) = 0;
    for (i=0; i<NUMCOMMANDS; i++)
        {
        targetLen = (int) strlen(p->string);
        if (tkLen <= targetLen)
            {
            for (j=0, numDiff=0; j<tkLen; j++)
                {
                if (ChangeCase(tk[j]) != ChangeCase(p->string[j]))
                    numDiff++;
                }
            if (numDiff == 0)
                {
                (*numMatches)++;
                commandPtr = p;
                if (tkLen == targetLen)
                    break;
                }
            }
        p++;
        }

    inValidCommand = NO;
    if (*numMatches == 1)
        {
        inValidCommand = YES;
        return (NO_ERROR);
        }
    else
        return (ERROR);
}


int FindValidParam (char *tk, int *numMatches)
{
    int         i, j, tkLen, targetLen, numDiff;
    CmdType     *p;
    ParmInfoPtr q;

    if (commandPtr)
        p = commandPtr;
    else
        {
        MrBayesPrint ("%s   Command pointer is NULL\n", spacer);
        return (ERROR);
        }
    tkLen = (int) strlen(tk);

    *numMatches = 0;
    for (i=0; i<p->numParms; i++)
        {
        q = paramTable + (p->parmList[i]);
        targetLen = (int) strlen(q->string);
        /* printf ("%s %d (%s %d)\n", q->string, targetLen, tk, p->numParms); */
        if (!strcmp(q->string, "Xxxxxxxxxx"))
            {
            (*numMatches)++;
            paramPtr = q;
            }
        else if (tkLen <= targetLen)
            {
            for (j=0, numDiff=0; j<tkLen; j++)
                {
                if (ChangeCase(tk[j]) != ChangeCase(q->string[j]))
                    numDiff++;
                }
            if (numDiff == 0)
                {
                (*numMatches)++;
                paramPtr = q;
                if (tkLen == targetLen)
                    break;
                }
            }   
        }
    
    if (*numMatches == 1)
        return (NO_ERROR);
    else
        return (ERROR);
}


int FreeCharacters (void)
{
    int     i, memoryLetFree;
    
    memoryLetFree = NO;

    if (memAllocs[ALLOC_TMPSET] == YES)
        {
        if (numChar > numTaxa)
            tempSet = (int *) SafeRealloc ((void *)tempSet, (size_t)numTaxa*sizeof(int));
            tempSetNeg = (int *) SafeRealloc ((void *)tempSetNeg, (size_t)numTaxa*sizeof(int));
        }
    if (memAllocs[ALLOC_MATRIX] == YES)
        {
        free (matrix);
        matrix = NULL;
        defMatrix = NO;
        memAllocs[ALLOC_MATRIX] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_CHARINFO] == YES)
        {
        free (charInfo);
        charInfo = NULL;
        memAllocs[ALLOC_CHARINFO] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_CHARSETS] == YES)
        {
        for (i=0; i<numCharSets; i++)
            {
            free (charSetNames[i]);
            free (charSet[i]);
            }
        free (charSetNames);
        free (charSet);
        charSetNames = NULL;
        charSet = NULL;
        numCharSets = 0;
        memAllocs[ALLOC_CHARSETS] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_PARTITIONS] == YES)
        {
        for (i=0; i<numDefinedPartitions; i++)
            free (partitionNames[i]);
        free (partitionNames);
        partitionNames = NULL;
        for (i=0; i<numChar; i++)
            free (partitionId[i]);
        free (partitionId);
        numDefinedPartitions = 0;
        memAllocs[ALLOC_PARTITIONS] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_PARTITIONVARS] == YES)
        {
        free (numVars);
        numVars = NULL;
        free (tempNum);
        tempNum = NULL;
        free (activeParams[0]);
        activeParams[0] = NULL;
        free (linkTable[0]);
        linkTable[0] = NULL;
        tempLinkUnlinkVec = NULL;
        activeParts = NULL;
        tempLinkUnlinkVec = NULL;
        for (i=0; i<NUM_LINKED; i++)
            {
            linkTable[i] = NULL;
            activeParams[i] = NULL;
            }
        memAllocs[ALLOC_PARTITIONVARS] = NO;
        memoryLetFree = YES;
        }

    ResetCharacterFlags();

    if (memoryLetFree == YES)
        MrBayesPrint ("%s   Deleting previously defined characters\n", spacer);

    return (NO_ERROR);
}


int FreeMatrix (void)
{
    if (FreeCharacters() == ERROR)
        return ERROR;

    return (FreeTaxa());
}


int FreeTaxa (void)
{
    int i, memoryLetFree;

    memoryLetFree = NO;
    if (memAllocs[ALLOC_TAXA] == YES)
        {
        if (taxaNames)
            {
            for (i=0; i<taxonCount; i++)
                free (taxaNames[i]);
            }
        free (taxaNames);
        taxaNames = NULL;
        free (taxaInfo);
        taxaInfo = NULL;
        free (tipCalibration);
        tipCalibration = NULL;
        numTaxa = 0;
        memAllocs[ALLOC_TAXA] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_TMPSET] == YES)
        {
        free (tempSet);
        tempSet = NULL;
        free (tempSetNeg);
        tempSetNeg = NULL;
        memAllocs[ALLOC_TMPSET] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_TAXASETS] == YES)
        {
        for (i=0; i<numTaxaSets; i++)
            {
            free (taxaSetNames[i]);
            free (taxaSet[i]);
            }
        free (taxaSetNames);
        taxaSetNames = NULL;
        free (taxaSet);
        taxaSet = NULL;
        numTaxaSets = 0;
        memAllocs[ALLOC_TAXASETS] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_SPECIESPARTITIONS] == YES)
        {
        for (i=0; i<numDefinedSpeciespartitions; i++)
            free (speciespartitionNames[i]);
        free (speciespartitionNames);
        speciespartitionNames = NULL;
        for (i=0; i<numTaxa; i++)
            free (speciespartitionId[i]);
        free (speciespartitionId);
        speciespartitionId = NULL;
        numDefinedSpeciespartitions = 0;
        memAllocs[ALLOC_SPECIESPARTITIONS] = NO;
        memoryLetFree = YES;
        }
    if (memAllocs[ALLOC_CONSTRAINTS] == YES)
        {
        for (i=0; i<numDefinedConstraints; i++)
            {
            free(definedConstraint[i]);
            free(definedConstraintTwo[i]);
            free(definedConstraintPruned[i]);
            free(definedConstraintTwoPruned[i]);
            free (constraintNames[i]);
            }
        free (definedConstraint);
        definedConstraint = NULL;
        free (definedConstraintTwo);
        definedConstraintTwo = NULL;
        free (definedConstraintsType);
        definedConstraintsType = NULL;
        free (constraintNames);
        constraintNames = NULL;
        free (nodeCalibration);
        nodeCalibration = NULL;
        numDefinedConstraints = 0;
        free (tempActiveConstraints);
        tempActiveConstraints = NULL;
        memAllocs[ALLOC_CONSTRAINTS] = NO;
        memoryLetFree = YES;
        }
    if (numUserTrees > 0)
        {
        MrBayesPrint ("%s   Deleting user trees\n", spacer);
        for (i=0; i<numUserTrees; i++)
            {
            FreePolyTree(userTree[i]);
            userTree[i] = NULL;
            }
        numUserTrees = 0;
        }

    FreeCharacters();

    if (memoryLetFree == YES)
        MrBayesPrint ("%s   Deleting previously defined taxa\n", spacer);

    /* reinitialize taxa variables */
    ResetTaxaFlags();

    return NO_ERROR;
}


int GetNumPartDivisions (int n)
{
    int         i, maxDiv, numDivs, *divFound;
    
    maxDiv = 0;
    for (i=0; i<numChar; i++)
        if (partitionId[i][n] > maxDiv)
            maxDiv = partitionId[i][n];

    divFound = (int *) SafeCalloc (maxDiv, sizeof(int));
    
    for (i=0; i<maxDiv; i++)
        divFound[i] = NO;
    
    for (i=0; i<numChar; i++)
        divFound[partitionId[i][n]] = YES;
        
    numDivs = 0;
    for (i=0; i<maxDiv; i++)
        if (divFound[i] == YES)
            numDivs++;
    
    free (divFound);

    return (numDivs + 1);
}


int GetToken (char *token, int *tokenType, char **sourceH)
{
    int             allNumbers, foundExp, foundExpSign;
    register char   *temp;
    char            *tempMax;
    
    (*tokenType) = 0;
    temp = token;
    tempMax = temp + CMD_STRING_LENGTH - 10;
    
    while (IsWhite(**sourceH) == 1 || IsWhite(**sourceH) == 2)
        {
        if (IsWhite(**sourceH) == 2)
            {
            *tokenType = RETURNSYMBOL;
            /* foundNewLine = YES;  Why is this commented out?? */
            /* MrBayesPrint ("RETURN\n"); */
            }
        ++(*sourceH);
        }
    
    if (readWord == YES && **sourceH != '"')
        {
        if (**sourceH==';')
            {
            *temp++ = ';';
            *tokenType = SEMICOLON;
            }
        else
            {
            while (isgraph(**sourceH) && **sourceH!=';')
                {
                if (temp > tempMax)
                    {
                    *tokenType = NOTHING;
                    token[20]='\0';
                    MrBayesPrint ("%s   Error while parsing a string. Token \"%s...[followed by at least %d  more charectors]\" is too long.\n", spacer,token,tempMax-token-20);
                    MrBayesPrint ("%s   Maximum allowed lenght of a token is %d\n", spacer,tempMax-token);
                    return (ERROR);
                    }
                *temp++ = *(*sourceH)++;
                }
            *tokenType = ALPHA;
            }
        *temp = '\0';
        readWord = NO;
        return (NO_ERROR);;
        }

    *tokenType = UNKNOWN_TOKEN_TYPE;
    if (IsIn(**sourceH,"="))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = EQUALSIGN;
        }
    else if (IsIn(**sourceH,";"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = SEMICOLON;
        }
    else if (IsIn(**sourceH,":"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = COLON;
        }
    else if (IsIn(**sourceH,","))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = COMMA;
        }
    else if (IsIn(**sourceH,"#"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = POUNDSIGN;
        }
    else if (IsIn(**sourceH,"("))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = LEFTPAR;
        }
    else if (IsIn(**sourceH,")"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = RIGHTPAR;
        }
    else if (IsIn(**sourceH,"{"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = LEFTCURL;
        }
    else if (IsIn(**sourceH,"}"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = RIGHTCURL;
        }
    else if (IsIn(**sourceH,"["))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = LEFTCOMMENT;
        }
    else if (IsIn(**sourceH,"]"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = RIGHTCOMMENT;
        }
    else if (IsIn(**sourceH,"?"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = QUESTIONMARK;
        }
    else if (IsIn(**sourceH,"-"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = DASH;
        }
    else if (IsIn(**sourceH,"$"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = DOLLAR;
        }
    else if (IsIn(**sourceH,"\"") && readWord == YES)
        {
        (*sourceH)++;
        while (**sourceH != '"' && **sourceH != '\0')
            {
            if (temp > tempMax)
                {
                *tokenType = NOTHING;
                token[20]='\0';
                MrBayesPrint ("%s   Error while parsing a string. Token \"%s...[followed by at least %d  more charectors]\" is too long.\n", spacer,token,tempMax-token-20);
                MrBayesPrint ("%s   Maximum allowed lenght of a token is %d\n", spacer,tempMax-token);
                return (ERROR);
                }
            *temp++ = *((*sourceH)++);
            }
        *temp='\0';
        *tokenType = ALPHA;
        (*sourceH)++;
        readWord = NO;
        }
    else if (IsIn(**sourceH,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789."))
        {
        if (IsIn(**sourceH,"0123456789."))
            allNumbers = TRUE;
        else
            allNumbers = FALSE;
        foundExp = foundExpSign = FALSE;
        *temp++ = *(*sourceH)++;
        while (IsIn(**sourceH,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.-+"))
            {
            if (temp > tempMax)
                {
                *tokenType = NOTHING;
                token[20]='\0';
                MrBayesPrint ("%s   Error while parsing a string. Token \"%s...[followed by at least %d  more charectors]\" is too long.\n", spacer,token,tempMax-token-20);
                MrBayesPrint ("%s   Maximum allowed lenght of a token is %d\n", spacer,tempMax-token);
                return (ERROR);
                }
            if (allNumbers == TRUE && !IsIn((*sourceH)[-1],"Ee") && **sourceH=='-')
                break;
            else if (allNumbers == TRUE && IsIn(**sourceH,"Ee") && foundExp == NO)
                foundExp = TRUE;
            else if (allNumbers == TRUE && IsIn(**sourceH,"+-") && IsIn((*sourceH)[-1],"Ee"))
                foundExpSign = TRUE;
            else if (!IsIn(**sourceH,"0123456789."))
                allNumbers = FALSE;
            *temp++ = *(*sourceH)++;
            }
        if (allNumbers == TRUE)
            *tokenType = NUMBER;
        else
            *tokenType = ALPHA;
        }
    else if (IsIn(**sourceH,"*"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = ASTERISK;
        }
    else if (IsIn(**sourceH,"/"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = FORWARDSLASH;
        }
    else if (IsIn(**sourceH,"'\\'"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = BACKSLASH;
        }
    else if (IsIn(**sourceH,"!"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = EXCLAMATIONMARK;
        }
    else if (IsIn(**sourceH,"%"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = PERCENT;
        }
    else if (IsIn(**sourceH,"\""))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = QUOTATIONMARK;
        }
    else if (IsIn(**sourceH,"&"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = AMPERSAND;
        }
    else if (IsIn(**sourceH,"~+^@{}`><"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = WEIRD;
        }
    else if (IsIn(**sourceH,"|"))
        {
        *temp++ = *(*sourceH)++;
        *tokenType = VERTICALBAR;
        }

    *temp = '\0';
    return (NO_ERROR);
}


int GetUserHelp (char *helpTkn)
{
    int         i, j, k, tempInt;
    char        tempString[100];
    Model       *mp;
    
    if (!strcmp(helpTkn, "Begin"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Begin                                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used to format data or commands in the program. The correct   \n");
        MrBayesPrint ("   usage is                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin <data or mrbayes>;                                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The two valid uses of the \"begin\" command, then, are                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin data;                                                                \n");
        MrBayesPrint ("      begin mrbayes;                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The \"data\" specifier is used to specify the beginning of a data block; your \n");
        MrBayesPrint ("   character data should follow. For example, the following is an example of     \n");
        MrBayesPrint ("   a data block for four taxa and ten DNA sites:                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin data;                                                                \n");
        MrBayesPrint ("         dimensions ntax=4 nchar=10;                                             \n");
        MrBayesPrint ("         format datatype=dna;                                                    \n");
        MrBayesPrint ("         matrix                                                                  \n");
        MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
        MrBayesPrint ("         taxon_2  AAGGATTCCA                                                     \n");
        MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
        MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        MrBayesPrint ("         ;                                                                       \n");
        MrBayesPrint ("      end;                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The other commands -- dimensions, format, and matrix -- are discussed         \n");
        MrBayesPrint ("   in the appropriate help menu. The only thing to note here is that the         \n");
        MrBayesPrint ("   block begins with a \"begin data\" command. The \"mrbayes\" command is        \n");
        MrBayesPrint ("   used to enter commands specific to the MrBayes program into the file.         \n");
        MrBayesPrint ("   This allows you to automatically process commands on execution of the         \n");
        MrBayesPrint ("   program. The following is a simple mrbayes block:                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin mrbayes;                                                             \n");
        MrBayesPrint ("         charset first  = 1-10\\3;                                               \n");
        MrBayesPrint ("         charset second = 2-10\\3;                                               \n");
        MrBayesPrint ("         charset third  = 3-10\\3;                                               \n");
        MrBayesPrint ("      end;                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This mrbayes block sets off the three \"charset\" commands, used to           \n");
        MrBayesPrint ("   predefine some blocks of characters. The mrbayes block can be very useful.    \n");
        MrBayesPrint ("   For example, in this case, it would save you the time of typing the char-     \n");
        MrBayesPrint ("   acter sets each time you executed the file. Also, note that every             \n");
        MrBayesPrint ("   \"begin <data or mrbayes>\" command ends with an \"end\". Finally, you can    \n");
        MrBayesPrint ("   have so-called foreign blocks in the file. An example of a foreign block      \n");
        MrBayesPrint ("   would be \"begin paup\". The program will simply skip this block. This is     \n");
        MrBayesPrint ("   useful because it means that you can use the same file for MrBayes, PAUP*     \n");
        MrBayesPrint ("   or MacClade (although it isn't clear why you would want to use those other    \n");
        MrBayesPrint ("   programs).                                                                    \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "End"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   End                                                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used to terminate a data or mrbayes block. The correct        \n");
        MrBayesPrint ("   usage is                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      end;                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For more information on this, check the help for the \"begin\" command.       \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Endblock"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Endblock                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This is an older, deprecated version of \"End\", see that command.            \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Plot"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Plot                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command plots specified parameters in the .p file or one of the .p files \n");
        MrBayesPrint ("   created during an MCMC analysis. An x-y graph of the parameter over the course\n");
        MrBayesPrint ("   of the chain is created. The command can be useful for visually diagnosing    \n");
        MrBayesPrint ("   convergence for many of the parameters of the phylogenetic model. The para-   \n");
        MrBayesPrint ("   meter to be plotted is specified by the \"parameter\" option. Several para-   \n");
        MrBayesPrint ("   meters can be plotted at once by using the \"match\" option, which has a      \n");
        MrBayesPrint ("   default value of \"perfect\". For example, if you were to set \"parameter = pi\"\n");
        MrBayesPrint ("   and \"match = consistentwith\", then all of the state frequency parameters    \n");
        MrBayesPrint ("   would be plotted. You can also set \"match=all\", in which case all of the    \n");
        MrBayesPrint ("   parameters are plotted.                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that the \"Sump\" command provides a different set of convergence diag-  \n");
        MrBayesPrint ("   nostics tools that you may also want to explore. Unlike \"Plot\", \"Sump\" can\n");
        MrBayesPrint ("   compare two or more parameter samples and will calculate convergence diagnos- \n");
        MrBayesPrint ("   tics as wel as parameter summaries for the pooled sample.                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Relburnin     -- If this option is set to 'Yes', then a proportion of the     \n");
        MrBayesPrint ("                    samples will be discarded as burnin when creating the plot.  \n");
        MrBayesPrint ("                    The proportion to be discarded is set with Burninfrac (see   \n");
        MrBayesPrint ("                    Burninfrac below). When the Relburnin option is set to 'No', \n");
        MrBayesPrint ("                    then a specific number of samples is discarded instead. This \n");
        MrBayesPrint ("                    number is set by Burnin (see below). Note that the burnin    \n");
        MrBayesPrint ("                    setting is shared across the 'comparetree', 'sump' and 'sumt'\n");
        MrBayesPrint ("                    commands.                                                    \n");
        MrBayesPrint ("   Burnin        -- Determines the number of samples (not generations) that will \n");
        MrBayesPrint ("                    be discarded when summary statistics are calculated. The     \n");
        MrBayesPrint ("                    value of this option is only relevant when Relburnin is set  \n");
        MrBayesPrint ("                    to 'No'.                                                     \n");
        MrBayesPrint ("   Burninfrac    -- Determines the fraction of samples that will be discarded    \n");
        MrBayesPrint ("                    when creating a plot. The value of this parameter is only    \n");
        MrBayesPrint ("                    relevant when Relburnin is set to 'Yes'. Example: A value of \n");
        MrBayesPrint ("                    this option of 0.25 means that 25%% of the samples will be   \n");
        MrBayesPrint ("                    discarded.                                                   \n");
        MrBayesPrint ("   Filename      -- The name of the file to plot.                                \n");
        MrBayesPrint ("   Parameter     -- Specification of parameters to be plotted. See above for     \n");
        MrBayesPrint ("                    details.                                                     \n");
        MrBayesPrint ("   Match         -- Specifies how to match parameter names to the Parameter      \n");
        MrBayesPrint ("                    specification. See above for details.                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Current settings:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter       Options                      Current Setting                  \n");
        MrBayesPrint ("   ------------------------------------------------------------                  \n");
        MrBayesPrint ("   Relburnin       Yes/No                       %s                               \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        MrBayesPrint ("   Burnin          <number>                     %d                               \n", chainParams.chainBurnIn);
        MrBayesPrint ("   Burninfrac      <number>                     %1.2lf                           \n", chainParams.burninFraction);
        MrBayesPrint ("   Filename        <name>                       %s                               \n", plotParams.plotFileName);
        MrBayesPrint ("   Parameter       <name>                       %s                               \n", plotParams.parameter);
        MrBayesPrint ("   Match           Perfect/Consistentwith/All   %s                               \n", plotParams.match);
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Dimensions"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Dimensions                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used in a data block to define the number of taxa and         \n");
        MrBayesPrint ("   characters. The correct usage is                                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      dimensions ntax=<number> nchar=<number>                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The dimensions must be the first command in a data block. The following       \n");
        MrBayesPrint ("   provides an example of the proper use of this command:                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin data;                                                                \n");
        MrBayesPrint ("         dimensions ntax=4 nchar=10;                                             \n");
        MrBayesPrint ("         format datatype=dna;                                                    \n");
        MrBayesPrint ("         matrix                                                                  \n");
        MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
        MrBayesPrint ("         taxon_2  AAGGATTCCA                                                     \n");
        MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
        MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        MrBayesPrint ("         ;                                                                       \n");
        MrBayesPrint ("      end;                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Here, the dimensions command tells MrBayes to expect a matrix with four       \n");
        MrBayesPrint ("   taxa and 10 characters.                                                       \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Format"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Format                                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used in a data block to define the format of the char-        \n");
        MrBayesPrint ("   acter matrix. The correct usage is                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      format datatype=<name> ... <parameter>=<option>                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The format command must be the second command in a data block. The following  \n");
        MrBayesPrint ("   provides an example of the proper use of this command:                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin data;                                                                \n");
        MrBayesPrint ("         dimensions ntax=4 nchar=10;                                             \n");
        MrBayesPrint ("         format datatype=dna gap=-;                                              \n");
        MrBayesPrint ("         matrix                                                                  \n");
        MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
        MrBayesPrint ("         taxon_2  AAGGAT--CA                                                     \n");
        MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
        MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        MrBayesPrint ("         ;                                                                       \n");
        MrBayesPrint ("      end;                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Here, the format command tells MrBayes to expect a matrix with DNA char-      \n");
        MrBayesPrint ("   acters and with gaps coded as \"-\".                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The following are valid options for format:                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Datatype   -- This parameter MUST BE INCLUDED in the format command. More-    \n");
        MrBayesPrint ("                 over, it must be the first parameter in the line. The           \n");
        MrBayesPrint ("                 datatype command specifies what type of characters are          \n");
        MrBayesPrint ("                 in the matrix. The following are valid options:                 \n");
        MrBayesPrint ("                    Datatype = Dna: DNA states (A,C,G,T,R,Y,M,K,S,W,H,B,         \n");
        MrBayesPrint ("                               V,D,N)                                            \n");
        MrBayesPrint ("                    Datatype = Rna: DNA states (A,C,G,U,R,Y,M,K,S,W,H,B,         \n");
        MrBayesPrint ("                               V,D,N)                                            \n");
        MrBayesPrint ("                    Datatype = Protein: Amino acid states (A,R,N,D,C,Q,E,        \n");
        MrBayesPrint ("                               G,H,I,L,K,M,F,P,S,T,W,Y,V)                        \n");
        MrBayesPrint ("                    Datatype = Restriction: Restriction site (0,1) states        \n");
        MrBayesPrint ("                    Datatype = Standard: Morphological (0,1) states              \n");
        MrBayesPrint ("                    Datatype = Continuous: Real number valued states             \n");
        MrBayesPrint ("                    Datatype = Mixed(<type>:<range>,...,<type>:<range>): A       \n");
        MrBayesPrint ("                               mixture of the above datatypes. For example,      \n");
        MrBayesPrint ("                               \"datatype=mixed(dna:1-100,protein:101-200)\"     \n");
        MrBayesPrint ("                               would specify a mixture of DNA and amino acid     \n");
        MrBayesPrint ("                               characters with the DNA characters occupying      \n");
        MrBayesPrint ("                               the first 100 sites and the amino acid char-      \n");
        MrBayesPrint ("                               acters occupying the last 100 sites.              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Interleave -- This parameter specifies whether the data matrix is in          \n");
        MrBayesPrint ("                 interleave format. The valid options are \"Yes\" or \"No\",     \n");
        MrBayesPrint ("                 with \"No\" as the default. An interleaved matrix looks like    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    format datatype=dna gap=- interleave=yes;                    \n");
        MrBayesPrint ("                    matrix                                                       \n");
        MrBayesPrint ("                    taxon_1  AACGATTCGT                                          \n");
        MrBayesPrint ("                    taxon_2  AAGGAT--CA                                          \n");
        MrBayesPrint ("                    taxon_3  AACGACTCCT                                          \n");
        MrBayesPrint ("                    taxon_4  AAGGATTCCT                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    taxon_1  CCTGGTAC                                            \n");
        MrBayesPrint ("                    taxon_2  CCTGGTAC                                            \n");
        MrBayesPrint ("                    taxon_3  ---GGTAG                                            \n");
        MrBayesPrint ("                    taxon_4  ---GGTAG                                            \n");
        MrBayesPrint ("                    ;                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Gap        -- This parameter specifies the format for gaps. Note that         \n");
        MrBayesPrint ("                 gap character can only be a single character and that it        \n");
        MrBayesPrint ("                 cannot correspond to a standard state (e.g., A,C,G,T,R,Y,       \n");
        MrBayesPrint ("                 M,K,S,W,H,B,V,D,N for nucleotide data).                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Missing    -- This parameter specifies the format for missing data. Note      \n");
        MrBayesPrint ("                 that the missing character can only be a single character and   \n");
        MrBayesPrint ("                 cannot correspond to a standard state (e.g., A,C,G,T,R,Y,       \n");
        MrBayesPrint ("                 M,K,S,W,H,B,V,D,N for nucleotide data). This is often an        \n");
        MrBayesPrint ("                 unnecessary parameter to set because many data types, such      \n");
        MrBayesPrint ("                 as nucleotide or amino acid, already have a missing char-       \n");
        MrBayesPrint ("                 acter specified. However, for morphological or restriction      \n");
        MrBayesPrint ("                 site data, \"missing=?\" is often used to specify ambiguity     \n");
        MrBayesPrint ("                 or unobserved data.                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Matchchar  -- This parameter specifies the matching character for the         \n");
        MrBayesPrint ("                 matrix. For example,                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    format datatype=dna gap=- matchchar=.;                       \n");
        MrBayesPrint ("                    matrix                                                       \n");
        MrBayesPrint ("                    taxon_1  AACGATTCGT                                          \n");
        MrBayesPrint ("                    taxon_2  ..G...--CA                                          \n");
        MrBayesPrint ("                    taxon_3  .....C..C.                                          \n");
        MrBayesPrint ("                    taxon_4  ..G.....C.                                          \n");
        MrBayesPrint ("                    ;                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                 is equivalent to                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    format datatype=dna gap=-;                                   \n");
        MrBayesPrint ("                    matrix                                                       \n");
        MrBayesPrint ("                    taxon_1  AACGATTCGT                                          \n");
        MrBayesPrint ("                    taxon_2  AAGGAT--CA                                          \n");
        MrBayesPrint ("                    taxon_3  AACGACTCCT                                          \n");
        MrBayesPrint ("                    taxon_4  AAGGATTCCT                                          \n");
        MrBayesPrint ("                    ;                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The only non-standard NEXUS format option is the use of the \"mixed\",        \n");
        MrBayesPrint ("   \"restriction\", \"standard\" and \"continuous\" datatypes. Hence, if         \n");
        MrBayesPrint ("   you use any of these datatype specifiers, a program like PAUP* or             \n");
        MrBayesPrint ("   MacClade will report an error (as they should because MrBayes is not          \n");
        MrBayesPrint ("   strictly NEXUS compliant).                                                    \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Matrix"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Matrix                                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command specifies the actual data for the phylogenetic analysis.         \n");
        MrBayesPrint ("   The character matrix should follow the dimensions and format commands         \n");
        MrBayesPrint ("   in a data block. The matrix can have all of the characters for a taxon        \n");
        MrBayesPrint ("   on a single line:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin data;                                                                \n");
        MrBayesPrint ("         dimensions ntax=4 nchar=10;                                             \n");
        MrBayesPrint ("         format datatype=dna gap=-;                                              \n");
        MrBayesPrint ("         matrix                                                                  \n");
        MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
        MrBayesPrint ("         taxon_2  AAGGAT--CA                                                     \n");
        MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
        MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        MrBayesPrint ("         ;                                                                       \n");
        MrBayesPrint ("      end;                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or be in \"interleaved\" format:                                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin data;                                                                \n");
        MrBayesPrint ("         dimensions ntax=4 nchar=20;                                             \n");
        MrBayesPrint ("         format datatype=dna gap=- interleave=yes;                               \n");
        MrBayesPrint ("         matrix                                                                  \n");
        MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
        MrBayesPrint ("         taxon_2  AAGGAT--CA                                                     \n");
        MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
        MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("         taxon_1  TTTTCGAAGC                                                     \n");
        MrBayesPrint ("         taxon_2  TTTTCGGAGC                                                     \n");
        MrBayesPrint ("         taxon_3  TTTTTGATGC                                                     \n");
        MrBayesPrint ("         taxon_4  TTTTCGGAGC                                                     \n");
        MrBayesPrint ("         ;                                                                       \n");
        MrBayesPrint ("      end;                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that the taxon names must not have spaces. If you really want to         \n");
        MrBayesPrint ("   indicate a space in a taxon name (perhaps between a genus and species         \n");
        MrBayesPrint ("   name), then you might use an underline (\"_\"). There should be at            \n");
        MrBayesPrint ("   least a single space after the taxon name, separating the name from           \n");
        MrBayesPrint ("   the actual data on that line. There can be spaces between the char-           \n");
        MrBayesPrint ("   acters.                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   If you have mixed data, then you specify all of the data in the same          \n");
        MrBayesPrint ("   matrix. Here is an example that includes two different data types:            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      begin data;                                                                \n");
        MrBayesPrint ("         dimensions ntax=4 nchar=20;                                             \n");
        MrBayesPrint ("         format datatype=mixed(dna:1-10,standard:21-30) interleave=yes;          \n");
        MrBayesPrint ("         matrix                                                                  \n");
        MrBayesPrint ("         taxon_1  AACGATTCGT                                                     \n");
        MrBayesPrint ("         taxon_2  AAGGAT--CA                                                     \n");
        MrBayesPrint ("         taxon_3  AACGACTCCT                                                     \n");
        MrBayesPrint ("         taxon_4  AAGGATTCCT                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("         taxon_1  0001111111                                                     \n");
        MrBayesPrint ("         taxon_2  0111110000                                                     \n");
        MrBayesPrint ("         taxon_3  1110000000                                                     \n");
        MrBayesPrint ("         taxon_4  1000001111                                                     \n");
        MrBayesPrint ("         ;                                                                       \n");
        MrBayesPrint ("      end;                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The matrix command is terminated by a semicolon.                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Finally, just a note on data presentation. It is much easier for others       \n");
        MrBayesPrint ("   to (1) understand your data and (2) repeat your analyses if you make          \n");
        MrBayesPrint ("   your data clean, comment it liberally (using the square brackets), and        \n");
        MrBayesPrint ("   embed the commands you used in a publication in the mrbayes block.            \n");
        MrBayesPrint ("   Remember that the data took a long time for you to collect. You might         \n");
        MrBayesPrint ("   as well spend a little time making the data file look nice and clear to       \n");
        MrBayesPrint ("   any that may later request the data for further analysis.                     \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Pairs"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Pairs                                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used to specify pairs of nucleotides. For example, your       \n");
        MrBayesPrint ("   data may be RNA sequences with a known secondary structure of stems and       \n");
        MrBayesPrint ("   loops. Substitutions in nucleotides involved in a Watson-Crick pairing        \n");
        MrBayesPrint ("   in stems are not strictly independent; a change in one changes the prob-      \n");
        MrBayesPrint ("   ability of a change in the partner. A solution to this problem is to          \n");
        MrBayesPrint ("   expand the model around the pair of nucleotides in the stem. This             \n");
        MrBayesPrint ("   command allows you to do this. The correct usage is:                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      pairs <NUC1>:<NUC2>, <NUC1>:<NUC2>,..., <NUC1>:<NUC2>;                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example,                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      pairs 30:56, 31:55, 32:54, 33:53, 34:52, 35:51, 36:50;                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   specifies pairings between nucleotides 30 and 56, 31 and 55, etc. Only        \n");
        MrBayesPrint ("   nucleotide data (DNA or RNA) may be paired using this command. Note that      \n");
        MrBayesPrint ("   in order for the program to actually implement a \"doublet\" model            \n");
        MrBayesPrint ("   involving a 16 X 16 rate matrix, you must specify that the structure of       \n");
        MrBayesPrint ("   the model is 16 X 16 using \"lset nucmodel=doublet\".                         \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Databreaks"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Databreaks                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used to specify breaks in your input data matrix. Your        \n");
        MrBayesPrint ("   data may be a mixture of genes or a mixture of different types of data.       \n");
        MrBayesPrint ("   Some of the models implemented by MrBayes account for nonindependence at      \n");
        MrBayesPrint ("   adjacent characters. The autocorrelated gamma model, for example, allows      \n");
        MrBayesPrint ("   rates at adjacent sites to be correlated. However, there is no way for        \n");
        MrBayesPrint ("   such a model to tell whether two sites, adjacent in the matrix, are           \n");
        MrBayesPrint ("   actually separated by many kilobases or megabases in the genome. The          \n");
        MrBayesPrint ("   databreaks command allows you to specify such breaks. The correct             \n");
        MrBayesPrint ("   usage is:                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      databreaks <break 1> <break 2> <break 3> ...                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example, say you have a data matrix of 3204 characters that include       \n");
        MrBayesPrint ("   nucleotide data from three genes. The first gene covers characters 1 to       \n");
        MrBayesPrint ("   970, the second gene covers characters 971 to 2567, and the third gene        \n");
        MrBayesPrint ("   covers characters 2568 to 3204. Also, let's assume that the genes are         \n");
        MrBayesPrint ("   not directly adjacent to one another in the genome, as might be likely        \n");
        MrBayesPrint ("   if you have mitochondrial sequences. In this case, you can specify            \n");
        MrBayesPrint ("   breaks between the genes using:                                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      databreaks 970 2567;                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The first break, between genes one and two, is after character 970 and        \n");
        MrBayesPrint ("   the second break, between genes two and three, is after character 2567.       \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Acknowledgments"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Acknowledgments                                                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the authors' acknowledgments.                              \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "About"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   About                                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command provides some general information about the program.             \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Version"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Version                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the release version of the program.                        \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Citations"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Citations                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows a thorough list of citations you may consider using        \n");
        MrBayesPrint ("   when publishing the results of a MrBayes analysis.                            \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showmatrix"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Showmatrix                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the character matrix currently in memory.                  \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showbeagle"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Showbeagle                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows available BEAGLE resources.                                \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Speciespartition"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Speciespartition                                                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Defines a partition of tips into species. The format for the speciespartition \n");
        MrBayesPrint ("   command is                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      Speciespartition <name> = <species name>:<taxon list> ,...,<sp nm>:<tx lst>\n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The command enumerates comma separated list of pairs consisting of 'species   \n");
        MrBayesPrint ("   name' and 'taxon list'. The 'taxon list' is a standard taxon list, as used by \n");
        MrBayesPrint ("   the 'Taxset' command. This means that you can use either the index or the name\n");
        MrBayesPrint ("   of a sequence ('taxon'). Ranges are specified using a dash, and a period can  \n");
        MrBayesPrint ("   be used as a synonym of the last sequence in the matrix.                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For exammple: speciespartition species = SpeciesA: 1, SpeciesB: 2-.           \n");
        MrBayesPrint ("   Here, we name two species. SpeciesA is represented by a single sequence while \n");
        MrBayesPrint ("   SpeciesB is represented by all remaining sequences in the matrix.             \n");
        MrBayesPrint ("   Each sequence is specified by its row index in the data matrix.               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   As with ordinary partitioning you may define multiple species partitioning    \n");
        MrBayesPrint ("   scheme. You have to use command 'set speciespartition' to enable use of one of\n");
        MrBayesPrint ("   them.                                                                         \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Currently defined Speciespartitions:                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Number  Speciespartition name        Number of species                        \n");
        MrBayesPrint ("   --------------------------------------------------------------------------    \n");
        for (i=0; i<numDefinedSpeciespartitions; i++)
            {
            tempInt=0;
            for (j=0; j<numTaxa; j++)
                {
                if (tempInt < speciespartitionId[j][i])
                    tempInt = speciespartitionId[j][i];
                }
            MrBayesPrint ("   %4d    %-24.24s   %4d",i+1, speciespartitionNames[i], tempInt);
            MrBayesPrint ("\n");
            }
        MrBayesPrint ("                                                                                \n");
        MrBayesPrint ("   --------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Constraint"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Constraint                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command defines a tree constraint. The format for the constraint         \n");
        MrBayesPrint ("   command is                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      constraint <name> [hard|negative|partial] = <taxon list> [:<taxon list>]   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   There are three types of constraint implemented in MrBayes. The type of the   \n");
        MrBayesPrint ("   constraint is specified by using one of the three keywords 'hard', 'negative',\n");
        MrBayesPrint ("   or 'partial' right after the name of the constraint. If no type is specified, \n");
        MrBayesPrint ("   then the constraint is assumed to be 'hard'.                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   In a rooted tree, a 'hard' constraint forces the taxa in the list to form a   \n");
        MrBayesPrint ("   monophyletic group. In an unrooted tree, the taxon split that separates the   \n");
        MrBayesPrint ("   taxa in the list from other taxa is forced to be present. The interpretation  \n");
        MrBayesPrint ("   of this depends on whether the tree is rooted on a taxon outside the list or  \n");
        MrBayesPrint ("   a taxon in the list. If the outgroup is excluded , the taxa in the list are   \n");
        MrBayesPrint ("   assumed to form a monophyletic group, but if the outgroup is included, the    \n");
        MrBayesPrint ("   taxa that are not in the list are forced together.                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   A 'negative' constraint bans all the trees that have the listed taxa in the   \n");
        MrBayesPrint ("   same subtree. In other words, it is the opposite of a hard constraint.        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   A 'partial' or backbone constraint is defined in terms of two sets of taxa    \n");
        MrBayesPrint ("   separated by a colon character. The constraint forces all taxa in the first   \n");
        MrBayesPrint ("   list to form a monophyletic group that does not include any taxon in the      \n");
        MrBayesPrint ("   second list. Taxa that are not included in either list can be placed in any   \n");
        MrBayesPrint ("   position on the tree, either inside or outside the constrained group. In an   \n");
        MrBayesPrint ("   unrooted tree, the two taxon lists can be switched with each other with no    \n");
        MrBayesPrint ("   effect. For a rooted tree, it is the taxa in the first list that have to be   \n");
        MrBayesPrint ("   monophyletic, that is, these taxa must share a common ancestor not shared with\n");
        MrBayesPrint ("   any taxon in the second list. The taxa in the second list may or may not fall \n");
        MrBayesPrint ("   in a monophyletic group depending on the rooting of the tree.                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   A list of taxa can be specified using a taxset, taxon names, taxon numbers, or\n");
        MrBayesPrint ("   any combination of the above, sepatated by spaces. The constraint is treated  \n");
        MrBayesPrint ("   as an absolute requirement of trees, that is, trees that are not compatible   \n");
        MrBayesPrint ("   with the constraint have zero prior (and hence zero posterior) probabilty.    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   If you are interested in inferring ancestral states for a particular node,    \n");
        MrBayesPrint ("   you need to 'hard' constrain that node first using the 'constraint' command.  \n");
        MrBayesPrint ("   The same applies if you wish to calibrate an interior node in a dated         \n");
        MrBayesPrint ("   analysis. For more information on how to infer ancestral states, see the help \n");
        MrBayesPrint ("   for the 'report' command. For more on dating, see the 'calibrate' command.    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   It is important to note that simply defining a constraint using this          \n");
        MrBayesPrint ("   command is not sufficient for the program to actually implement the           \n");
        MrBayesPrint ("   constraint in an analysis. You must also enforce the constraints using        \n");
        MrBayesPrint ("   'prset topologypr = constraints (<list of constraints>)'. For more infor-     \n");
        MrBayesPrint ("   mation on this, see the help on the 'prset' command.                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Examples:                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      constraint myclade = Homo Pan Gorilla                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Defines a hard constraint forcing Homo, Pan, and Gorilla to form a mono-      \n");
        MrBayesPrint ("   phyletic group or a split that does not include any other taxa.               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      constraint forbiddenclade negative = Homo Pan Gorilla                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Defines a negative constraint that associates all trees where Homon, Pan, and \n");
        MrBayesPrint ("   Gorilla form a monophyletic group with zero posterior probability. In other   \n");
        MrBayesPrint ("   words, such trees will not be sampled during MCMC.                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      constraint backbone partial = Homo Gorilla : Mus                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Defines a partial constraint that keeps Mus outside of the clade defined by   \n");
        MrBayesPrint ("   the most recent common ancestor of Homo and Gorilla. Other taxa are allowed to\n");
        MrBayesPrint ("   sit anywhere in the tree. Note that this particular constraint is meaningless \n");
        MrBayesPrint ("   in unrooted trees. MrBayes does not assume anything about the position of the \n");
        MrBayesPrint ("   outgroup unless it is explicitly included in the partial constraint. Therefore\n");
        MrBayesPrint ("   a partial constraint must have at least two taxa on each side of the ':' to be\n");
        MrBayesPrint ("   useful in analyses of unrooted trees. The case is different for rooted trees, \n");
        MrBayesPrint ("   where it is sufficient for a partial constraint to have more than one taxon   \n");
        MrBayesPrint ("   before the ':', as in the example given above, to constrain tree space.       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   To define a more complex constraint tree, simply combine constraints into a   \n");
        MrBayesPrint ("   list when issuing the 'prset topologypr' command.                             \n");
        MrBayesPrint ("                                                                                 \n");
        if (numDefinedConstraints > 0)
            {
            MrBayesPrint ("   Currently defined constraints:                                                \n");
            MrBayesPrint ("                                                                                 \n");
            MrBayesPrint ("   Number  Constraint name          type      Number of taxa in[:out]            \n");
            MrBayesPrint ("   --------------------------------------------------------------------------    \n");       
            }
        for (i=0; i<numDefinedConstraints; i++)
            {
            strncpy (tempString, constraintNames[i], 22);
            MrBayesPrint ("   %4d    %-22.22s   ",i+1, tempString);
            if (definedConstraintsType[i] == HARD)
                MrBayesPrint ("hard      ");
            else if (definedConstraintsType[i] == PARTIAL)
                MrBayesPrint ("partial   ");
            else
                {
                assert (definedConstraintsType[i] == NEGATIVE);
                MrBayesPrint ("negative  ");
                }
            k = NumBits (definedConstraint[i], numTaxa/nBitsInALong + 1);
            MrBayesPrint ("%d", k);
            if (definedConstraintsType[i] == PARTIAL)
                {
                k = NumBits (definedConstraintTwo[i], numTaxa/nBitsInALong + 1);
                MrBayesPrint (":%d", k);
                }
        MrBayesPrint ("\n");
            }
        MrBayesPrint ("                                                                                \n");
        MrBayesPrint ("   --------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Calibrate"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Calibrate                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command dates a terminal or interior node in the tree. The format is     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      calibrate <node_name> = <age_prior>                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   where <node_name> is the name of a defined interior constraint node or the    \n");
        MrBayesPrint ("   name of a terminal node (tip) and <age_prior> is a prior probability distribu-\n");
        MrBayesPrint ("   tion on the age of the node. The latter can either be a fixed date or a date  \n");
        MrBayesPrint ("   drawn from one of the available prior probability distributions. In general,  \n");
        MrBayesPrint ("   the available prior probability distributions are parameterized in terms of   \n");
        MrBayesPrint ("   the expected mean age of the distribution to facilitate for users. Some dis-  \n");
        MrBayesPrint ("   tributions put a positive probability on all ages above 0.0, while others in- \n");
        MrBayesPrint ("   clude a minimum-age constraint and sometimes a maximum-age constraint. The    \n");
        MrBayesPrint ("   available distributions and their parameters are:                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      calibrate <node_name> = fixed(<age>)                                       \n");
        MrBayesPrint ("      calibrate <node_name> = uniform(<min_age>,<max_age>)                       \n");
        MrBayesPrint ("      calibrate <node_name> = offsetexponential(<min_age>,<mean_age>)            \n");
        MrBayesPrint ("      calibrate <node_name> = truncatednormal(<min_age>,<mean_age>,<stdev>)      \n");
        MrBayesPrint ("      calibrate <node_name> = lognormal(<mean_age>,<stdev>)                      \n");
        MrBayesPrint ("      calibrate <node_name> = offsetlognormal(<min_age>,<mean_age>,<stdev>)      \n");
        MrBayesPrint ("      calibrate <node_name> = gamma(<mean_age>,<stdev>)                          \n");
        MrBayesPrint ("      calibrate <node_name> = offsetgamma(<min_age>,<mean_age>,<stdev>)          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that mean_age is always the mean age and stdev the standard deviation of \n");
        MrBayesPrint ("   the distribution measured in user-defined time units. This way of specifying  \n");
        MrBayesPrint ("   the distribution parameters is often different from the parameterization used \n");
        MrBayesPrint ("   elsewhere in the program. For instance, the standard parameters of the gamma  \n");
        MrBayesPrint ("   distribution used by MrBayes are shape (alpha) and rate (beta). If you want   \n");
        MrBayesPrint ("   to use the standard parameterization, the conversions are as follows:         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      exponential distributon: mean    = 1 / rate                                \n");
        MrBayesPrint ("      gamma distributon:       mean    = alpha / beta                            \n");
        MrBayesPrint ("                               st.dev. = square_root (alpha / beta^2)            \n");
        MrBayesPrint ("      lognormal distributon:   mean    = exp (mean_log + st.dev._log^2/2)        \n");
        MrBayesPrint ("                               st.dev. = square_root ((exp (st.dev._log^2) - 1)  \n");
        MrBayesPrint ("                                         * (exp (2*mean_log + st.dev._log^2))    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The truncated normal distribution is an exception in that the mean_age and    \n");
        MrBayesPrint ("   stdev parameters are the mean and standard deviation of the underlying non-   \n");
        MrBayesPrint ("   truncated normal distribution. The truncation will cause the modified distri- \n");
        MrBayesPrint ("   bution to have a higher mean and lower standard deviation. The magnitude of   \n");
        MrBayesPrint ("   that effect depends on how much of the tail of the distribution is removed.   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that previous to version 3.2.2, MrBayes used the standard rate parameter-\n");
        MrBayesPrint ("   ization of the offset exponential. This should not cause a problem in most    \n");
        MrBayesPrint ("   cases because the old parameterization will result in an error in more recent \n");
        MrBayesPrint ("   versions of MrBayes, and the likely source of the error is given in the error \n");
        MrBayesPrint ("   message.                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For a practical example, assume that we had three fossil terminals named      \n");
        MrBayesPrint ("   'FossilA', 'FossilB', and 'FossilC'. Assume further that we want to fix the   \n");
        MrBayesPrint ("   age of FossilA to 100.0 million years, we think that FossilB is somewhere     \n");
        MrBayesPrint ("   between 100.0 and 200.0 million years old, and that FossilC is at least 300.0 \n");
        MrBayesPrint ("   million years old, possibly older but relatively unlikely to be more than     \n");
        MrBayesPrint ("   400.0 million years old. Then we might use the commands:                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      calibrate FossilA = fixed(100) FossilB = uniform(100,200)                  \n");
        MrBayesPrint ("      calibrate FossilC = offsetexponential(300,400)                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that it is possible to give more than one calibration for each           \n");
        MrBayesPrint ("   'calibrate' statement. Thus, 'calibrate FossilA=<setting> FossilB=<setting>'  \n");
        MrBayesPrint ("   would be a valid statement.                                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   To actually use the calibrations to obtain dated trees, you also need to set  \n");
        MrBayesPrint ("   a clock model using relevant 'brlenspr' and 'nodeagepr' options of the 'prset'\n");
        MrBayesPrint ("   command. You may also want to examine the 'clockvarpr' and 'clockratepr' op-  \n");
        MrBayesPrint ("   tions. Furthermore, you need to activate the relevant constraint(s) using     \n");
        MrBayesPrint ("   'topologypr', if you use any dated interior nodes in the tree.                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   You may wish to remove a calibration from an interior or terminal node, which \n");
        MrBayesPrint ("   has previously been calibrated. You can do that using                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      calibrate <node_name> = unconstrained                                      \n");
        MrBayesPrint ("                                                                                 \n");
        j = 0;
        for (i=0; i<numTaxa; i++)
            if (tipCalibration[i].prior != unconstrained)
                j++;
        for (i=0; i<numDefinedConstraints; i++)
            if (nodeCalibration[i].prior != unconstrained)
                j++;
        if (j > 0)
            {
            MrBayesPrint ("                                                                                 \n");
            MrBayesPrint ("   Currently defined calibrations:                                               \n");
            MrBayesPrint ("                                                                                 \n");
            MrBayesPrint ("   Node name                Type       Calibration                               \n");
            MrBayesPrint ("   ------------------------------------------------------------------            \n");       
            for (i=0; i<numTaxa+numDefinedConstraints; i++)
                {
                if (i<numTaxa)
                    calibrationPtr = &tipCalibration[i];
                else
                    calibrationPtr = &nodeCalibration[i-numTaxa];
                if (calibrationPtr != NULL && calibrationPtr->prior != unconstrained)
                    {
                    if (i<numTaxa)
                        strncpy (tempString, taxaNames[i], 22);
                    else
                        strncpy (tempString, constraintNames[i-numTaxa], 22);
                    if (i<numTaxa)
                        MrBayesPrint ("   %-22.22s   Terminal   %s\n", tempString, calibrationPtr->name);
                    else
                        MrBayesPrint ("   %-22.22s   Interior   %s\n", tempString, calibrationPtr->name);
                    }
                }
            }
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showmodel"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Showmodel                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the current model settings. The correct usage is           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      showmodel                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   After typing \"showmodel\", the modelling assumptions are shown on a          \n");
        MrBayesPrint ("   partition-by-partition basis.                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Execute"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Execute                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command executes a file called <file name>. The correct usage is:        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      execute <file name>                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example,                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      execute replicase.nex                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   would execute the file named \"replicase.nex\". This file must be in the      \n");
        MrBayesPrint ("   same directory as the executable.                                             \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Lset"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Lset                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command sets the parameters of the likelihood model. The likelihood      \n");
        MrBayesPrint ("   function is the probability of observing the data conditional on the phylo-   \n");
        MrBayesPrint ("   genetic model. In order to calculate the likelihood, you must assume a        \n");
        MrBayesPrint ("   model of character change. This command lets you tailor the biological        \n");
        MrBayesPrint ("   assumptions made in the phylogenetic model. The correct usage is              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      lset <parameter>=<option> ... <parameter>=<option>                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example, \"lset nst=6 rates=gamma\" would set the model to a general      \n");
        MrBayesPrint ("   model of DNA substition (the GTR) with gamma-distributed rate variation       \n");
        MrBayesPrint ("   across sites.                                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Applyto   -- This option allows you to apply the lset commands to specific    \n");
        MrBayesPrint ("                partitions. This command should be the first in the list of      \n");
        MrBayesPrint ("                commands specified in lset. Moreover, it only makes sense to     \n");
        MrBayesPrint ("                be using this command if the data have been partitioned. A       \n");
        MrBayesPrint ("                default partition is set on execution of a matrix. If the data   \n");
        MrBayesPrint ("                are homogeneous (i.e., all of the same data type), then this     \n");
        MrBayesPrint ("                partition will not subdivide the characters. Up to 30 other      \n");
        MrBayesPrint ("                partitions can be defined, and you can switch among them using   \n");
        MrBayesPrint ("                \"set partition=<partition name>\". Now, you may want to         \n");
        MrBayesPrint ("                specify different models to different partitions of the data.    \n");
        MrBayesPrint ("                Applyto allows you to do this. For example, say you have         \n");
        MrBayesPrint ("                partitioned the data by codon position, and you want to apply    \n");
        MrBayesPrint ("                a nst=2 model to the first two partitions and nst=6 to the       \n");
        MrBayesPrint ("                last. This could be implemented in two uses of lset:             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                   lset applyto=(1,2) nst=2                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                   lset applyto=(3) nst=6                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                The first applies the parameters after \"applyto\" to the        \n");
        MrBayesPrint ("                first and second partitions. The second lset applies nst=6       \n");
        MrBayesPrint ("                to the third partition. You can also use applyto=(all), which    \n");
        MrBayesPrint ("                attempts to apply the parameter settings to all of the data      \n");
        MrBayesPrint ("                partitions. Importantly, if the option is not consistent with    \n");
        MrBayesPrint ("                the data in the partition, the program will not apply the        \n");
        MrBayesPrint ("                lset option to that partition.                                   \n");
        MrBayesPrint ("   Nucmodel  -- This specifies the general form of the nucleotide substitution   \n");
        MrBayesPrint ("                model. The options are \"4by4\" [the standard model of DNA       \n");
        MrBayesPrint ("                substitution in which there are only four states (A,C,G,T/U)],   \n");
        MrBayesPrint ("                \"doublet\" (a model appropriate for modelling the stem regions  \n");
        MrBayesPrint ("                of ribosomal genes where the state space is the 16 doublets of   \n");
        MrBayesPrint ("                nucleotides), \"codon\" (the substitution model is expanded      \n");
        MrBayesPrint ("                around triplets of nucleotides--a codon), and \"Protein\"        \n");
        MrBayesPrint ("                (triplets of nucleotides are translated to amino acids, which    \n");
        MrBayesPrint ("                form the basis of the substitution model).                       \n");
        MrBayesPrint ("   Nst       -- Sets the number of substitution types: \"1\" constrains all of   \n");
        MrBayesPrint ("                the rates to be the same (e.g., a JC69 or F81 model); \"2\" all- \n");
        MrBayesPrint ("                ows transitions and transversions to have potentially different  \n");
        MrBayesPrint ("                rates (e.g., a K80 or HKY85 model); \"6\" allows all rates to    \n");
        MrBayesPrint ("                be different, subject to the constraint of time-reversibility    \n");
        MrBayesPrint ("                (e.g., a GTR model). Finally, 'nst' can be set to 'mixed', which \n");
        MrBayesPrint ("                results in the Markov chain sampling over the space of all poss- \n");
        MrBayesPrint ("                ible reversible substitution models, including the GTR model and \n");
        MrBayesPrint ("                all models that can be derived from it model by grouping the six \n");
        MrBayesPrint ("                rates in various combinations. This includes all the named models\n");
        MrBayesPrint ("                above and a large number of others, with or without name.        \n");
        MrBayesPrint ("   Code      -- Enforces the use of a particular genetic code. The default       \n");
        MrBayesPrint ("                is the universal code. Other options include \"vertmt\" for      \n");
        MrBayesPrint ("                vertebrate mitocondrial, \"invermt\", \"mycoplasma\", \"yeast\", \n");
        MrBayesPrint ("                \"ciliate\", \"echinoderm\", \"euplotid\", and \"metmt\" (for    \n");
        MrBayesPrint ("                metazoan mitochondrial except vertebrates).                      \n");
        MrBayesPrint ("   Ploidy    -- Specifies the ploidy of the organism. Options are \"Haploid\",   \n");
        MrBayesPrint ("                \"Diploid\" or \"Zlinked\". This option is used when a coalescent\n");
        MrBayesPrint ("                prior is used on trees.                                          \n");
        MrBayesPrint ("   Rates     -- Sets the model for among-site rate variation. In general, the    \n");
        MrBayesPrint ("                rate at a site is considered to be an unknown random variable.   \n");
        MrBayesPrint ("                The valid options are:                                           \n");
        MrBayesPrint ("                * equal    -- No rate variation across sites.                    \n");
        MrBayesPrint ("                * gamma    -- Gamma-distributed rates across sites. The rate     \n");
        MrBayesPrint ("                              at a site is drawn from a gamma distribution.      \n");
        MrBayesPrint ("                              The gamma distribution has a single parameter      \n");
        MrBayesPrint ("                              that describes how much rates vary.                \n");
        MrBayesPrint ("                * lnorm    -- Lognormal-distributed rates across sites. The      \n");
        MrBayesPrint ("                              rate at a site is drawn from a lognormal           \n");
        MrBayesPrint ("                              distribution. the lognormal distribiton has a      \n");
        MrBayesPrint ("                              single parameter, sigma (SD) that describes how    \n");
        MrBayesPrint ("                              much rates vary (mean fixed to log(1.0) == 0.0.    \n");
        MrBayesPrint ("                * adgamma  -- Autocorrelated rates across sites. The marg-       \n");
        MrBayesPrint ("                              inal rate distribution is gamma, but adjacent      \n");
        MrBayesPrint ("                              sites have correlated rates.                       \n");
        MrBayesPrint ("                * propinv  -- A proportion of the sites are invariable.          \n");
        MrBayesPrint ("                * invgamma -- A proportion of the sites are invariable while     \n");
        MrBayesPrint ("                              the rate for the remaining sites are drawn from    \n");
        MrBayesPrint ("                              a gamma distribution.                              \n");
        MrBayesPrint ("                * kmixture -- Site rates come from a mixture with k categories.  \n");
        MrBayesPrint ("                              Category rates are drawn from an ordered flat      \n");
        MrBayesPrint ("                              Dirichlet distribution with mean rather than sum   \n");
        MrBayesPrint ("                              equal to 1.0.                                      \n");
        MrBayesPrint ("                Note that MrBayes versions 2.0 and earlier supported options     \n");
        MrBayesPrint ("                that allowed site specific rates (e.g., ssgamma). In versions    \n");
        MrBayesPrint ("                3.0 and later, site specific rates are allowed, but set using    \n");
        MrBayesPrint ("                the 'prset ratepr' command for each partition.                   \n");
        MrBayesPrint ("   Ngammacat -- Sets the number of rate categories for the gamma distribution.   \n");
        MrBayesPrint ("                The gamma distribution is continuous. However, it is virtually   \n");
        MrBayesPrint ("                impossible to calculate likelihoods under the continuous gamma   \n");
        MrBayesPrint ("                distribution. Hence, an approximation to the continuous gamma    \n");
        MrBayesPrint ("                is used; the gamma distribution is broken into ncat categories   \n");
        MrBayesPrint ("                of equal weight (1/ncat). The mean rate for each category rep-   \n");
        MrBayesPrint ("                resents the rate for the entire cateogry. This option allows     \n");
        MrBayesPrint ("                you to specify how many rate categories to use when approx-      \n");
        MrBayesPrint ("                imating the gamma. The approximation is better as ncat is inc-   \n");
        MrBayesPrint ("                reased. In practice, \"ncat=4\" does a reasonable job of         \n");
        MrBayesPrint ("                approximating the continuous gamma.                              \n");
        MrBayesPrint ("   Nlnormcat -- Used to set the number of discrete categories used for the ap-   \n");
        MrBayesPrint ("                proximation of the lognormal distribution, in the same way as    \n");
        MrBayesPrint ("                the Ngammacat setting for the discrete gamma approximation.      \n");
        MrBayesPrint ("                Default value is 4.                                              \n");
        MrBayesPrint ("   Nmixtcat  -- Used to set the number of components in the k-mixture model of   \n");
        MrBayesPrint ("                rate variation across sites. Default value is 4.                 \n");
#if 0
        /* Temporarily disable this because of conflict with likelihood calculators. It should be renamed to samplerates when reintroduced. */
        MrBayesPrint ("   Usegibbs  -- Specifies whether site probabilities under the discrete gamma    \n");
        MrBayesPrint ("                model of rate variation across sites will be summed across rate  \n");
        MrBayesPrint ("                categories ('Usegibbs=No') or sampled using a Gibbs sampler      \n");
        MrBayesPrint ("                ('Usegibbs=Yes'). The Gibbs sampling approach is much faster and \n");
        MrBayesPrint ("                requires less memory but the likelihood of the sampled points    \n");
        MrBayesPrint ("                will be considerably higher than with the standard approach of   \n");
        MrBayesPrint ("                summing probabilities, so you need to be aware of this when com- \n");
        MrBayesPrint ("                paring your results with those you obtain with other programs.   \n");
        MrBayesPrint ("                Assume that you are using n rate categories in your discrete     \n");
        MrBayesPrint ("                gamma distribution. Then the Gibbs approach is up to n times     \n");
        MrBayesPrint ("                faster and requires 1/n as much memory as the standard method.   \n");
        MrBayesPrint ("                Unfortunately, the state space also becomes larger so the chain  \n");
        MrBayesPrint ("                may need more time to converge. The approach should work best    \n");
        MrBayesPrint ("                for large trees, where the uncertainty concerning the best rate  \n");
        MrBayesPrint ("                category for each site is negligible. Gibbs sampling cannot be   \n");
        MrBayesPrint ("                used for the autocorrelated discrete gamma model, for standard   \n");
        MrBayesPrint ("                data, or for restriction data. Also, MrBayes will not use Gibbs  \n");
        MrBayesPrint ("                sampling when you want to infer site rates.                      \n");
        MrBayesPrint ("   Gibbsfreq -- Sets the frequency with which the rate categories of the discrete\n");
        MrBayesPrint ("                gamma will be Gibbs sampled. In practice, we have found that a   \n");
        MrBayesPrint ("                resampling frequency of every 100 MCMC generations works well for\n");
        MrBayesPrint ("                reasonably long runs. The more frequent the Gibbs sampling, the  \n");
        MrBayesPrint ("                slower the Gibbs sampling approach will be. If you have k rate   \n");
        MrBayesPrint ("                categories and Gibbs sample them every n generations, then the   \n");
        MrBayesPrint ("                time it takes to complete n generations will roughly be propor-  \n");
        MrBayesPrint ("                tional to n+k. Compare this with the traditional approach of     \n");
        MrBayesPrint ("                summing across the n rate categories in every generation, which  \n");
        MrBayesPrint ("                requires time proportional to n*k. In practice, however, the     \n");
        MrBayesPrint ("                speed difference is not quite as large as this.                  \n");
#endif
        MrBayesPrint ("   Nbetacat  -- Sets the number of rate categories for the beta distribution.    \n");
        MrBayesPrint ("                A symmetric beta distribution is used to model the stationary    \n");
        MrBayesPrint ("                frequencies when morphological data are used. This option        \n");
        MrBayesPrint ("                specifies how well the beta distribution will be approximated.   \n");
        MrBayesPrint ("   Omegavar  -- Allows the nonsynonymous/synonymous rate ratio (omega) to vary   \n");
        MrBayesPrint ("                across codons. Ny98 assumes that there are three classes, with   \n");
        MrBayesPrint ("                potentially different omega values (omega1, omega2, omega3):     \n");
        MrBayesPrint ("                omega2 = 1; 0 < omega1 < 1; and omega3 > 1. Like the Ny98 model, \n");
        MrBayesPrint ("                the M3 model has three omega classes. However, their values are  \n");
        MrBayesPrint ("                less constrained, with omega1 < omega2 < omega3. The default     \n");
        MrBayesPrint ("                (omegavar = equal) has no variation on omega across sites.       \n");
        MrBayesPrint ("   Covarion  -- This forces the use of a covarion-like model of substitution     \n");
        MrBayesPrint ("                for nucleotide or amino acid data. The valid options are \"yes\" \n");
        MrBayesPrint ("                and \"no\". The covarion model allows the rate at a site to      \n");
        MrBayesPrint ("                change over its evolutionary history. Specifically, the site     \n");
        MrBayesPrint ("                is either on or off. When it is off, no substitutions are poss-  \n");
        MrBayesPrint ("                ible. When the process is on, substitutions occur according to   \n");
        MrBayesPrint ("                a specified substitution model (specified using the other        \n");
        MrBayesPrint ("                lset options).                                                   \n");
        MrBayesPrint ("   Coding    -- This specifies how characters were sampled. If all site patterns \n");
        MrBayesPrint ("                had the possibility of being sampled, then \"All\" should be     \n");
        MrBayesPrint ("                specified (the default). Otherwise \"Variable\" (only variable   \n");
        MrBayesPrint ("                characters had the possibility of being sampled), \"Informative\"\n");
        MrBayesPrint ("                (only parsimony informative characters has the possibility of    \n");
        MrBayesPrint ("                being sampled), \"Nosingletons\" (characters which are constant  \n");
        MrBayesPrint ("                in all but one taxon were not sampled), \"Noabsencesites\" (char-\n");
        MrBayesPrint ("                acters for which all taxa were coded as absent were not sampled),\n");
        MrBayesPrint ("                \"Nopresencesites\" (characters for which all taxa were coded as \n");
        MrBayesPrint ("                present were not sampled). \"All\" works for all data types.     \n");
        MrBayesPrint ("                However, the others only work for morphological (All/Variable/   \n");
        MrBayesPrint ("                Informative/Nosingletons) or restriction site (All/Variable/     \n");
        MrBayesPrint ("                Informative/Nosingletons/Noabsencesites/Nopresencesites/         \n");
        MrBayesPrint ("                Nosingletonpresence/Nosingletonabsence) data.                    \n");
        MrBayesPrint ("   Parsmodel -- This forces calculation under the so-called parsimony model      \n");
        MrBayesPrint ("                described by Tuffley and Steel (1998). The options are \"yes\"   \n");
        MrBayesPrint ("                or \"no\". Note that the biological assumptions of this model    \n");
        MrBayesPrint ("                are anything but parsimonious. In fact, this model assumes many  \n");
        MrBayesPrint ("                more parameters than the next most complicated model implemented \n");
        MrBayesPrint ("                in this program. If you really believe that the parsimony model  \n");
        MrBayesPrint ("                makes the biological assumptions described by Tuffley and Steel, \n");
        MrBayesPrint ("                then the parsimony method is miss-named.                         \n");
    /*  MrBayesPrint ("   Augment   -- This allows the chain to consider the missing entries of         \n");
        MrBayesPrint ("                the data matrix as random variables. A Gibbs sampler is          \n");
        MrBayesPrint ("                used to sample states.                                           \n"); */
        MrBayesPrint ("                                                                                 \n");
        if (numCurrentDivisions == 0)
            tempInt = 1;
        else
            tempInt = numCurrentDivisions;
        for (i=0; i<tempInt; i++)
            {
            if (numCurrentDivisions == 0)
                {
                MrBayesPrint ("   Default model settings:                                                       \n");
                mp = &defaultModel;
                }
            else
                {
                MrBayesPrint ("   Model settings for partition %d:                                              \n", i+1);
                mp = &modelParams[i];
                }
            MrBayesPrint ("                                                                                 \n");
            MrBayesPrint ("   Parameter    Options                               Current Setting            \n");
            MrBayesPrint ("   ------------------------------------------------------------------            \n");       
            MrBayesPrint ("   Nucmodel     4by4/Doublet/Codon/Protein              %s                       \n", mp->nucModel);
            MrBayesPrint ("   Nst          1/2/6/Mixed                             %s                       \n", mp->nst);
            MrBayesPrint ("   Code         Universal/Vertmt/Invermt/Yeast/Mycoplasma/                       \n");
            MrBayesPrint ("                Ciliate/Echinoderm/Euplotid/Metmt       %s                       \n", mp->geneticCode);
            MrBayesPrint ("   Ploidy       Haploid/Diploid/Zlinked                 %s                       \n", mp->ploidy);
            MrBayesPrint ("   Rates        Equal/Gamma/LNorm/Propinv/                                       \n");
            MrBayesPrint ("                Invgamma/Adgamma/Kmixture               %s                       \n", mp->ratesModel);
            MrBayesPrint ("   Ngammacat    <number>                                %d                       \n", mp->numGammaCats);
            MrBayesPrint ("   Nlnormcat    <number>                                %d                       \n", mp->numLnormCats);
            MrBayesPrint ("   Nmixtcat     <number>                                %d                       \n", mp->numMixtCats);
#if 0
/* Temporarily disable this because of conflict with likelihood calculators. It should be renamed to samplerates when reintroduced. */
            MrBayesPrint ("   Usegibbs     Yes/No                                  %s                       \n", mp->useGibbs);
            MrBayesPrint ("   Gibbsfreq    <number>                                %d                       \n", mp->gibbsFreq);
#endif
            MrBayesPrint ("   Nbetacat     <number>                                %d                       \n", mp->numBetaCats);
            MrBayesPrint ("   Omegavar     Equal/Ny98/M3                           %s                       \n", mp->omegaVar);
            MrBayesPrint ("   Covarion     No/Yes                                  %s                       \n", mp->covarionModel);
            MrBayesPrint ("   Coding       All/Variable/Informative/Nosingletons                            \n");
            MrBayesPrint ("                Noabsencesites/Nopresencesites/                                  \n");
            MrBayesPrint ("                Nosingletonabsence/Nosingletonpresence  %s                       \n", mp->codingString);
            MrBayesPrint ("   Parsmodel    No/Yes                                  %s                       \n", mp->parsModel);
        /*  MrBayesPrint ("   Augment      No/Yes                                  %s                       \n", mp->augmentData); */
            MrBayesPrint ("   ------------------------------------------------------------------            \n");       
            MrBayesPrint ("                                                                                 \n");
            }
        }
    else if (!strcmp(helpTkn, "Prset"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Prset                                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command sets the priors for the phylogenetic model. Remember that        \n");
        MrBayesPrint ("   in a Bayesian analysis, you must specify a prior probability distribution     \n");
        MrBayesPrint ("   for the parameters of the likelihood model. The prior distribution rep-       \n");
        MrBayesPrint ("   resents your prior beliefs about the parameter before observation of the      \n");
        MrBayesPrint ("   data. This command allows you to tailor your prior assumptions to a large     \n");
        MrBayesPrint ("   extent.                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Applyto       -- This option allows you to apply the prset commands to        \n");
        MrBayesPrint ("                    specific partitions. This command should be the first        \n");
        MrBayesPrint ("                    in the list of commands specified in prset. Moreover, it     \n");
        MrBayesPrint ("                    only makes sense to be using this command if the data        \n");
        MrBayesPrint ("                    have been partitioned. A default partition is set on         \n");
        MrBayesPrint ("                    execution of a matrix. If the data are homogeneous           \n");
        MrBayesPrint ("                    (i.e., all of the same data type), then this partition       \n");
        MrBayesPrint ("                    will not subdivide the characters. Up to 30 other part-      \n");
        MrBayesPrint ("                    itions can be defined, and you can switch among them using   \n");
        MrBayesPrint ("                    \"set partition=<partition name>\". Now, you may want to     \n");
        MrBayesPrint ("                    specify different priors to different partitions of the      \n");
        MrBayesPrint ("                    data. Applyto allows you to do this. For example, say        \n");
        MrBayesPrint ("                    you have partitioned the data by codon position, and         \n");
        MrBayesPrint ("                    you want to fix the statefreqs to equal for the first two    \n");
        MrBayesPrint ("                    partitions but apply a flat Dirichlet prior to the state-    \n");
        MrBayesPrint ("                    freqs of the last. This could be implemented in two uses of  \n");
        MrBayesPrint ("                    prset:                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset applyto=(1,2) statefreqs=fixed(equal)               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset applyto=(3) statefreqs=dirichlet(1,1,1,1)           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    The first applies the parameters after \"applyto\"           \n");
        MrBayesPrint ("                    to the first and second partitions. The second prset         \n");
        MrBayesPrint ("                    applies a flat Dirichlet to the third partition. You can     \n");
        MrBayesPrint ("                    also use applyto=(all), which attempts to apply the para-    \n");
        MrBayesPrint ("                    meter settings to all of the data partitions. Importantly,   \n");
        MrBayesPrint ("                    if the option is not consistent with the data in the part-   \n");
        MrBayesPrint ("                    ition, the program will not apply the prset option to        \n");
        MrBayesPrint ("                    that partition.                                              \n");
        MrBayesPrint ("   Tratiopr      -- This parameter sets the prior for the transition/trans-      \n");
        MrBayesPrint ("                    version rate ratio (tratio). The options are:                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset tratiopr = beta(<number>, <number>)                 \n");
        MrBayesPrint ("                       prset tratiopr = fixed(<number>)                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    The program assumes that the transition and transversion     \n");
        MrBayesPrint ("                    rates are independent gamma-distributed random variables     \n");
        MrBayesPrint ("                    with the same scale parameter when beta is selected. If you  \n");
        MrBayesPrint ("                    want a diffuse prior that puts equal emphasis on transition/ \n");
        MrBayesPrint ("                    transversion rate ratios above 1.0 and below 1.0, then use a \n");
        MrBayesPrint ("                    flat Beta, beta(1,1), which is the default. If you wish to   \n");
        MrBayesPrint ("                    concentrate this distribution more in the equal-rates region,\n");
        MrBayesPrint ("                    then use a prior of the type beta(x,x), where the magnitude  \n");
        MrBayesPrint ("                    of x determines how much the prior is concentrated in the    \n");
        MrBayesPrint ("                    equal rates region. For instance, a beta(20,20) puts more    \n");
        MrBayesPrint ("                    probability on rate ratios close to 1.0 than a beta(1,1). If \n");
        MrBayesPrint ("                    you think it is likely that the transition/transversion rate \n");
        MrBayesPrint ("                    ratio is 2.0, you can use a prior of the type beta(2x,x),    \n");
        MrBayesPrint ("                    where x determines how strongly the prior is concentrated on \n");
        MrBayesPrint ("                    tratio values near 2.0. For instance, a beta(2,1) is much    \n");
        MrBayesPrint ("                    more diffuse than a beta(80,40) but both have the expected   \n");
        MrBayesPrint ("                    tratio 2.0 in the absence of data. The parameters of the     \n");
        MrBayesPrint ("                    Beta can be interpreted as counts: if you have observed x    \n");
        MrBayesPrint ("                    transitions and y transversions, then a beta(x+1,y+1) is a   \n");
        MrBayesPrint ("                    good representation of this information. The fixed option    \n");
        MrBayesPrint ("                    allows you to fix the tratio to a particular value.          \n");
        MrBayesPrint ("   Revmatpr      -- This parameter sets the prior for the substitution rates     \n");
        MrBayesPrint ("                    of the GTR model for nucleotide data. The options are:       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset revmatpr = dirichlet(<number>,<number>,...,<number>)\n");
        MrBayesPrint ("                       prset revmatpr = fixed(<number>,<number>,...,<number>)    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    The program assumes that the six substitution rates          \n");
        MrBayesPrint ("                    are independent gamma-distributed random variables with the  \n");
        MrBayesPrint ("                    same scale parameter when dirichlet is selected. The six     \n");
        MrBayesPrint ("                    numbers in brackets each corresponds to a particular substi- \n");
        MrBayesPrint ("                    tution type. Together, they determine the shape of the prior.\n");
        MrBayesPrint ("                    The six rates are in the order A<->C, A<->G, A<->T, C<->G,   \n");
        MrBayesPrint ("                    C<->T, and G<->T. If you want an uninformative prior you can \n");
        MrBayesPrint ("                    use dirichlet(1,1,1,1,1,1), also referred to as a 'flat'     \n");
        MrBayesPrint ("                    Dirichlet. This is the default setting. If you wish a prior  \n");
        MrBayesPrint ("                    where the C<->T rate is 5 times and the A<->G rate 2 times   \n");
        MrBayesPrint ("                    higher, on average, than the transversion rates, which are   \n");
        MrBayesPrint ("                    all the same, then you should use a prior of the form        \n");
        MrBayesPrint ("                    dirichlet(x,2x,x,x,5x,x), where x determines how much the    \n");
        MrBayesPrint ("                    prior is focused on these particular rates. For more info,   \n");
        MrBayesPrint ("                    see tratiopr. The fixed option allows you to fix the substi- \n");
        MrBayesPrint ("                    tution rates to particular values.                           \n");
        MrBayesPrint ("   Revratepr     -- This parameter sets the prior for each substitution rate of  \n");
        MrBayesPrint ("                    the GTR model subspace when 'nst' is set to 'mixed' (see the \n");
        MrBayesPrint ("                    'lset' command). The only option is                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset revratepr = symdir(<number>)                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    which will associate each independent rate in the rate matrix\n");
        MrBayesPrint ("                    with a modified symmetric Dirichlet prior, where a singleton \n");
        MrBayesPrint ("                    rate has the specified alpha parameter, while a rate that    \n");
        MrBayesPrint ("                    applies to n pairwise substitution types has an alpha that is\n");
        MrBayesPrint ("                    n times the specified number. The higher the specified num-  \n");
        MrBayesPrint ("                    ber, the more focused the prior will be on equal rates. The  \n");
        MrBayesPrint ("                    default value is 1, which gives an effect similar to a flat  \n");
        MrBayesPrint ("                    Dirichlet.                                                   \n");
        MrBayesPrint ("   Aamodelpr     -- This parameter sets the rate matrix for amino acid data.     \n");
        MrBayesPrint ("                    You can either fix the model by specifying aamodelpr=fixed   \n");
        MrBayesPrint ("                    (<model name>), where <model name> is 'poisson' (a glorified \n");
        MrBayesPrint ("                    Jukes-Cantor model), 'jones', 'dayhoff', 'mtrev', 'mtmam',   \n");
        MrBayesPrint ("                    'wag', 'rtrev', 'cprev', 'vt', 'blosum', 'lg', 'equalin'     \n");
        MrBayesPrint ("                    (a glorified Felsenstein 1981 model), or 'gtr'. You can also \n");
        MrBayesPrint ("                    average over the first ten models by specifying aamodelpr=   \n");
        MrBayesPrint ("                    mixed. If you do so, the Markov chain will sample each model \n");
        MrBayesPrint ("                    according to its probability. The sampled model is reported  \n");
        MrBayesPrint ("                    as an index: poisson(0), jones(1), dayhoff(2), mtrev(3),     \n");
        MrBayesPrint ("                    mtmam(4), wag(5), rtrev(6), cprev(7), vt(8), or blosum(9).   \n");
        MrBayesPrint ("                    The 'Sump' command summarizes the MCMC samples and calculates\n");
        MrBayesPrint ("                    the posterior probability estimate for each of these models. \n");
        MrBayesPrint ("   Aarevmatpr    -- This parameter sets the prior for the substitution rates     \n");
        MrBayesPrint ("                    of the GTR model for amino acid data. The options are:       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset aarevmatpr = dirichlet(<number>,<number>,...,<number>)\n");
        MrBayesPrint ("                       prset aarevmatpr = fixed(<number>,<number>,...,<number>)  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    The options are the same as those for 'Revmatpr' except that \n");
        MrBayesPrint ("                    they are defined over the 190 rates of the time-reversible   \n");
        MrBayesPrint ("                    GTR model for amino acids instead of over the 6 rates of the \n");
        MrBayesPrint ("                    GTR model for nucleotides. The rates are in the order A<->R, \n");
        MrBayesPrint ("                    A<->N, etc to Y<->V. In other words, amino acids are listed  \n");
        MrBayesPrint ("                    in alphabetic order based on their full name. The first amino\n");
        MrBayesPrint ("                    acid (Alanine) is then combined in turn with all amino acids \n");
        MrBayesPrint ("                    following it in the list, starting with amino acid 2 (Argi-  \n");
        MrBayesPrint ("                    nine) and finishing with amino acid 20 (Valine). The second  \n");
        MrBayesPrint ("                    amino acid (Arginine) is then combined in turn with all amino\n");
        MrBayesPrint ("                    acids following it, starting with amino acid 3 (Asparagine)  \n");
        MrBayesPrint ("                    and finishing with amino acid 20 (Valine), and so on.        \n");
        MrBayesPrint ("   Omegapr       -- This parameter specifies the prior on the nonsynonymous/     \n");
        MrBayesPrint ("                    synonymous rate ratio. The options are:                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset omegapr = dirichlet(<number>,<number>)              \n");
        MrBayesPrint ("                       prset omegapr = fixed(<number>)                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
        MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        MrBayesPrint ("                    case when there is no variation in omega across sites (i.e., \n");
        MrBayesPrint ("                    \"lset omegavar=equal\").                                    \n");
        MrBayesPrint ("   Ny98omega1pr  -- This parameter specifies the prior on the nonsynonymous/     \n");
        MrBayesPrint ("                    synonymous rate ratio for sites under purifying selection.   \n");
        MrBayesPrint ("                    The options are:                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset Ny98omega1pr = beta(<number>,<number>)              \n");
        MrBayesPrint ("                       prset Ny98omega1pr = fixed(<number>)                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
        MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        MrBayesPrint ("                    case where omega varies across sites using the model of      \n");
        MrBayesPrint ("                    Nielsen and Yang (1998) (i.e., \"lset omegavar=ny98\"). If   \n");
        MrBayesPrint ("                    fixing the parameter, you must specify a number between      \n");
        MrBayesPrint ("                    0 and 1.                                                     \n");
        MrBayesPrint ("   Ny98omega3pr  -- This parameter specifies the prior on the nonsynonymous/     \n");
        MrBayesPrint ("                    synonymous rate ratio for positively selected sites. The     \n");
        MrBayesPrint ("                    options are:                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset Ny98omega3pr = uniform(<number>,<number>)           \n");
        MrBayesPrint ("                       prset Ny98omega3pr = exponential(<number>)                \n");
        MrBayesPrint ("                       prset Ny98omega3pr = fixed(<number>)                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
        MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        MrBayesPrint ("                    case where omega varies across sites according to the        \n");
        MrBayesPrint ("                    NY98 model. Note that if the NY98 model is specified         \n");
        MrBayesPrint ("                    that this parameter must be greater than 1, so you should    \n");
        MrBayesPrint ("                    not specify a uniform(0,10) prior, for example.              \n");
        MrBayesPrint ("   M3omegapr     -- This parameter specifies the prior on the nonsynonymous/     \n");
        MrBayesPrint ("                    synonymous rate ratios for all three classes of sites for    \n");
        MrBayesPrint ("                    the M3 model. The options are:                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset M3omegapr = exponential                             \n");
        MrBayesPrint ("                       prset M3omegapr = fixed(<number>,<number>,<number>)       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
        MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        MrBayesPrint ("                    case where omega varies across sites using the M3 model of   \n");
        MrBayesPrint ("                    Yang et al. (2000) (i.e., \"lset omegavar=M3\"). Under the   \n");
        MrBayesPrint ("                    exponential prior, the four rates (dN1, dN2, dN3, and dS)    \n");
        MrBayesPrint ("                    are all considered to be independent draws from the same     \n");
        MrBayesPrint ("                    exponential distribution (the parameter of the exponential   \n");
        MrBayesPrint ("                    does not matter, and so you don't need to specify it). The   \n");
        MrBayesPrint ("                    rates dN1, dN2, and dN3 are taken to be the order statistics \n");
        MrBayesPrint ("                    with dN1 < dN2 < dN3. These three rates are all scaled to    \n");
        MrBayesPrint ("                    the same synonymous rate, dS. The other option is to simply  \n");
        MrBayesPrint ("                    fix the three rate ratios to some values.                    \n");
        MrBayesPrint ("   Codoncatfreqs -- This parameter specifies the prior on frequencies of sites   \n");
        MrBayesPrint ("                    under purifying, neutral, and positive selection. The        \n");
        MrBayesPrint ("                    options are:                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset codoncatfreqs = dirichlet(<num>,<num>,<num>)        \n");
        MrBayesPrint ("                       prset codoncatfreqs = fixed(<number>,<number>,<number>)   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only in effect if the nucleotide sub-      \n");
        MrBayesPrint ("                    stitution model is set to codon using the lset command       \n");
        MrBayesPrint ("                    (lset nucmodel=codon). Moreover, it only applies to the      \n");
        MrBayesPrint ("                    case where omega varies across sites using the models of     \n");
        MrBayesPrint ("                    Nielsen and Yang (1998) (i.e., \"lset omegavar=ny98\")       \n");
        MrBayesPrint ("                    or Yang et al. (2000) (i.e., \"lset omegavar=M3\")           \n");
        MrBayesPrint ("                    Note that the sum of the three frequencies must be 1.        \n");
        MrBayesPrint ("   Statefreqpr   -- This parameter specifies the prior on the state freq-        \n");
        MrBayesPrint ("                    uencies. The options are:                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset statefreqpr = dirichlet(<number>)                   \n");
        MrBayesPrint ("                       prset statefreqpr = dirichlet(<number>,...,<number>)      \n");
        MrBayesPrint ("                       prset statefreqpr = fixed(equal)                          \n");
        MrBayesPrint ("                       prset statefreqpr = fixed(empirical)                      \n");
        MrBayesPrint ("                       prset statefreqpr = fixed(<number>,...,<number>)          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    For the dirichlet, you can specify either a single number    \n");
        MrBayesPrint ("                    or as many numbers as there are states. If you specify a     \n");
        MrBayesPrint ("                    single number, then the prior has all states equally         \n");
        MrBayesPrint ("                    probable with a variance related to the single parameter     \n");
        MrBayesPrint ("                    passed in.                                                   \n");
        MrBayesPrint ("   Shapepr       -- This parameter specifies the prior for the gamma/lnorm shape \n");
        MrBayesPrint ("                    parameter for among-site rate variation. The options are:    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset shapepr = uniform(<number>,<number>)                \n");
        MrBayesPrint ("                       prset shapepr = exponential(<number>)                     \n");
        MrBayesPrint ("                       prset shapepr = fixed(<number>)                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Pinvarpr      -- This parameter specifies the prior for the proportion of     \n");
        MrBayesPrint ("                    invariable sites. The options are:                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset pinvarpr = uniform(<number>,<number>)               \n");
        MrBayesPrint ("                       prset pinvarpr = fixed(<number>)                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    Note that the valid range for the parameter is between 0     \n");
        MrBayesPrint ("                    and 1. Hence, \"prset pinvarpr=uniform(0,0.8)\" is valid     \n");
        MrBayesPrint ("                    while \"prset pinvarpr=uniform(0,10)\" is not. The def-      \n");
        MrBayesPrint ("                    ault setting is \"prset pinvarpr=uniform(0,1)\".             \n");
        MrBayesPrint ("   Ratecorrpr    -- This parameter specifies the prior for the autocorrelation   \n");
        MrBayesPrint ("                    parameter of the autocorrelated gamma distribution for       \n");
        MrBayesPrint ("                    among-site rate variation. The options are:                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset ratecorrpr = uniform(<number>,<number>)             \n");
        MrBayesPrint ("                       prset ratecorrpr = fixed(<number>)                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    Note that the valid range for the parameter is between -1    \n");
        MrBayesPrint ("                    and 1. Hence, \"prset ratecorrpr=uniform(-1,1)\" is valid    \n");
        MrBayesPrint ("                    while \"prset ratecorrpr=uniform(-11,10)\" is not. The       \n");
        MrBayesPrint ("                    default setting is \"prset ratecorrpr=uniform(-1,1)\".       \n");
        MrBayesPrint ("   Covswitchpr   -- This option sets the prior for the covarion switching        \n");
        MrBayesPrint ("                    rates. The options are:                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset covswitchpr = uniform(<number>,<number>)            \n");
        MrBayesPrint ("                       prset covswitchpr = exponential(<number>)                 \n");
        MrBayesPrint ("                       prset covswitchpr = fixed(<number>,<number>)              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    The covarion model has two rates: a rate from on to off      \n");
        MrBayesPrint ("                    and a rate from off to on. The rates are assumed to have     \n");
        MrBayesPrint ("                    independent priors that individually are either uniformly    \n");
        MrBayesPrint ("                    or exponentially distributed. The other option is to         \n");
        MrBayesPrint ("                    fix the switching rates, in which case you must specify      \n");
        MrBayesPrint ("                    both rates. (The first number is off->on and the second      \n");
        MrBayesPrint ("                    is on->off).                                                 \n");
        MrBayesPrint ("   Symdirihyperpr - This option sets the prior for the stationary frequencies    \n");
        MrBayesPrint ("                    of the states for morphological (standard) data. There can   \n");
        MrBayesPrint ("                    be as many as 10 states for standard data. However, the      \n");
        MrBayesPrint ("                    labelling of the states is somewhat arbitrary. For example,  \n");
        MrBayesPrint ("                    the state \"1\" for different characters does not have the   \n");
        MrBayesPrint ("                    same meaning. This is not true for DNA characters, for ex-   \n");
        MrBayesPrint ("                    ample, where a \"G\" has the same meaning across characters. \n");
        MrBayesPrint ("                    The fact that the labelling of morphological characters is   \n");
        MrBayesPrint ("                    arbitrary makes it difficult to allow unequal character-     \n");
        MrBayesPrint ("                    state frequencies. MrBayes gets around this problem by       \n");
        MrBayesPrint ("                    assuming that the states have a symmetric Dirichlet prior    \n");
        MrBayesPrint ("                    (i.e. all Dirichlet parameters are equal). The variation in  \n");
        MrBayesPrint ("                    the Dirichlet can be controlled by this parameter.           \n");
        MrBayesPrint ("                    Symdirihyperpr specifies the distribution on the parameter   \n");
        MrBayesPrint ("                    of the symmetric Dirichlet. The valid options are:           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset Symdirihyperpr = uniform(<number>,<number>)         \n");
        MrBayesPrint ("                       prset Symdirihyperpr = exponential(<number>)              \n");
        MrBayesPrint ("                       prset Symdirihyperpr = fixed(<number>)                    \n");
        MrBayesPrint ("                       prset Symdirihyperpr = fixed(infinity)                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    If \"fixed(infinity)\" is chosen, the Dirichlet prior is     \n");
        MrBayesPrint ("                    fixed such that all character states have equal frequency.   \n");
        MrBayesPrint ("   Topologypr    -- This parameter specifies the prior probabilities of          \n");
        MrBayesPrint ("                    phylogenies. The options are:                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset topologypr = uniform                                \n");
        MrBayesPrint ("                       prset topologypr = speciestree                            \n");
        MrBayesPrint ("                       prset topologypr = constraints(<list>)                    \n");
        MrBayesPrint ("                       prset topologypr = fixed(<treename>)                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    If the prior is selected to be \"uniform\", the default,     \n");
        MrBayesPrint ("                    then all possible trees are considered a priori equally      \n");
        MrBayesPrint ("                    probable. The 'speciestree' option is used when the topology \n");
        MrBayesPrint ("                    is constrained to fold inside a species tree together with   \n");
        MrBayesPrint ("                    other (gene) trees. The constraints option allows you to     \n");
        MrBayesPrint ("                    specify complicated prior probabilities on trees (constraints\n");
        MrBayesPrint ("                    are discussed more fully in \"help constraint\"). Note that  \n");
        MrBayesPrint ("                    you must specify a list of constraints that you wish to be   \n");
        MrBayesPrint ("                    obeyed. The list can be either the constraints' name or      \n");
        MrBayesPrint ("                    number. Finally, you can fix the topology to that of a user  \n");
        MrBayesPrint ("                    tree defined in a trees block. Branch lengths will still be  \n");
        MrBayesPrint ("                    sampled as usual on the fixed topology.                      \n");
        MrBayesPrint ("   Brlenspr      -- This parameter specifies the prior probability dist-         \n");
        MrBayesPrint ("                    ribution on branch lengths. The options are specified using: \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset brlenspr = <setting>                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    where <setting> is one of                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       unconstrained:uniform(<num>,<num>)                        \n");
        MrBayesPrint ("                       unconstrained:exponential(<number>)                       \n");
        MrBayesPrint ("                       unconstrained:twoexp(<num>,<num>)                         \n");
        MrBayesPrint ("                       unconstrained:gammadir(<num>,<num>,<num>,<num>)           \n");
        MrBayesPrint ("                       unconstrained:invgamdir(<num>,<num>,<num>,<num>)          \n");
        MrBayesPrint ("                       clock:uniform                                             \n");
        MrBayesPrint ("                       clock:birthdeath                                          \n");
        MrBayesPrint ("                       clock:coalescence                                         \n");
        MrBayesPrint ("                       clock:fossilization                                       \n");
        MrBayesPrint ("                       clock:speciestree                                         \n");
        MrBayesPrint ("                       fixed(<treename>)                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    Trees with unconstrained branch lengths are unrooted         \n");
        MrBayesPrint ("                    whereas clock-constrained trees are rooted. The option       \n");
        MrBayesPrint ("                    after the colon specifies the details of the probability     \n");
        MrBayesPrint ("                    density of branch lengths. If you choose a birth-death       \n");
        MrBayesPrint ("                    or coalescence prior, you may want to modify the details     \n");
        MrBayesPrint ("                    of the parameters of those processes (speciation rate,       \n");
        MrBayesPrint ("                    extinction rate and sample probability for the birth-death   \n");
        MrBayesPrint ("                    prior; population size and clock rate parameter for the      \n");
        MrBayesPrint ("                    coalescence prior). When gene trees are constrained to fold  \n");
        MrBayesPrint ("                    inside species trees, the appropriate branch length prior is \n");
        MrBayesPrint ("                    'clock:speciestree'. Under this model, it is possible to     \n");
        MrBayesPrint ("                    control whether the population size is constant or variable  \n");
        MrBayesPrint ("                    across the species tree using the 'popvarpr' setting.        \n");
        MrBayesPrint ("                    Branch lengths can also be fixed but only if the topology is \n");
        MrBayesPrint ("                    fixed.                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    For unconstrained branch lengths, MrBayes offers five alter- \n");
        MrBayesPrint ("                    native prior distributions. The first two are the simple     \n");
        MrBayesPrint ("                    'uniform' and 'exponential' priors. The 'uniform' prior takes\n");
        MrBayesPrint ("                    two parameters, the lower and upper bound of the uniform dis-\n");
        MrBayesPrint ("                    tribution, respectively. The 'exponential' prior takes a sin-\n");
        MrBayesPrint ("                    gle parameter, the rate of the exponential distribution. The \n");
        MrBayesPrint ("                    mean of the exponential distribution is the inverse of the   \n");
        MrBayesPrint ("                    rate. For instance, an 'exp(10)' distribution has an expected\n");
        MrBayesPrint ("                    mean of 0.1.                                                 \n");
        MrBayesPrint ("                    MrBayes also offers three more complex prior distributions   \n");
        MrBayesPrint ("                    on unconstrained branch lengths. The two-exponential prior   \n");
        MrBayesPrint ("                    (Yang and Rannala 2005; Yang 2007) uses two different expo-  \n"); 
        MrBayesPrint ("                    nential distributions, one for internal and one for external \n");
        MrBayesPrint ("                    branch lengths. The two-exponential prior is invoked using   \n");
        MrBayesPrint ("                    'twoexp(<r_I>,<r_E>)', where '<r_I>' is a number specifying  \n");
        MrBayesPrint ("                    the rate of the exponential distribution on internal branch  \n");
        MrBayesPrint ("                    lengths, while '<r_E>' is the rate for external branch       \n");
        MrBayesPrint ("                    lengths. The prior mean for internal branch lengths is then  \n");
        MrBayesPrint ("                    1/r_I, and for external ones is 1/r_E. For instance, to set  \n");
        MrBayesPrint ("                    prior mean of internal branch lengths to 0.01, and external  \n");
        MrBayesPrint ("                    ones to 0.1, use 'twoexp(100,10)'.                           \n");
        MrBayesPrint ("                    The setting 'twoexp(10,10)' is equivalent to 'exp(10)'.      \n");
        MrBayesPrint ("                    The compound Dirichlet priors 'gammadir(<a_T>,<b_T>,<a>,<c>)'\n");
        MrBayesPrint ("                    and 'invgamdir(<a_T>,<b_T>,<a>,<c>)' specify a fairly diffuse\n");
        MrBayesPrint ("                    prior on tree length 'T', and then partition the tree length \n");
        MrBayesPrint ("                    into branch lengths according to a Dirichlet distribution    \n");
        MrBayesPrint ("                    (Rannala et al. 2012). If 'T' is considered drawn from a     \n");
        MrBayesPrint ("                    gamma distribution with parameters a_T and b_T, and with mean\n");
        MrBayesPrint ("                    a_T/b_T, we recommend setting a_T = 1; if it is instead con- \n");
        MrBayesPrint ("                    sidered drawn from an inverse gamma (invgamma) distribution  \n");
        MrBayesPrint ("                    with parameters a_T and b_T, and with mean b_T/(a_T -1), then\n");
        MrBayesPrint ("                    we reccommend setting a_T = 3. In the latter case, b_T should\n");
        MrBayesPrint ("                    be chosen so that the prior mean of T is reasonable for the  \n");
        MrBayesPrint ("                    data. In the former case, setting b_T = 0.1 (corresponding to\n");
        MrBayesPrint ("                    a mean tree length of 10) should be appropriate for a wide   \n");
        MrBayesPrint ("                    range of tree lengths (at least in the interval 1 to 100).   \n");
        MrBayesPrint ("                    The concentration parameter a of the Dirichlet distribution  \n");
        MrBayesPrint ("                    is inversely related to the variance of the branch lengths,  \n");
        MrBayesPrint ("                    while c is the ratio of the prior means for the internal and \n");
        MrBayesPrint ("                    external branch lengths. The default setting, a = c = 1,     \n");
        MrBayesPrint ("                    specifies a uniform Dirichlet distribution of branch lengths \n");
        MrBayesPrint ("                    given the tree length. For instance, 'gammadir(1,0.1,1,1)'   \n");
        MrBayesPrint ("                    specifies a compound Dirichlet prior on branch lengths, where\n");
        MrBayesPrint ("                    tree length is associated with a gamma distribution with mean\n");
        MrBayesPrint ("                    10, and branch length proportions are associated with a uni- \n");
        MrBayesPrint ("                    form Dirichlet distribution (default).                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    For clock trees with calibrated external nodes (fossils),    \n");
        MrBayesPrint ("                    MrBayes also offers the fossilized birth-death prior:        \n");
        MrBayesPrint ("                    'clock:fossilization'.                                       \n");
        MrBayesPrint ("                    If 'SampleStrat' is set to 'fossiltip', it assumes that upon \n");
        MrBayesPrint ("                    sampling the lineage is dead and won't produce descendants,  \n");
        MrBayesPrint ("                    meaning each fossil sample is a tip. If 'SampleStrat' is set \n");
        MrBayesPrint ("                    to 'random' (default), fossils are sampled serially along the\n");
        MrBayesPrint ("                    birth-death tree (Stadler 2010), so they can be tips or an-  \n");
        MrBayesPrint ("                    cestors. See 'Speciationpr', 'Extinctionpr', 'SampleStrat',  \n");
        MrBayesPrint ("                    'Fossilizationpr' for more information.                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Treeagepr     -- This parameter specifies the prior probability distribution  \n");
        MrBayesPrint ("                    on the tree age when a uniform or fossilization prior is used\n");
        MrBayesPrint ("                    on the branch lengths of a clock tree.                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    The options are:                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset treeagepr = <setting>                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    where <setting> is one of                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       fixed(<age>)                                              \n");
        MrBayesPrint ("                       uniform(<min_age>,<max_age>)                              \n");
        MrBayesPrint ("                       offsetexponential(<min_age>,<mean_age>)                   \n");
        MrBayesPrint ("                       truncatednormal(<min_age>,<mean_age>,<st.dev.>)           \n");
        MrBayesPrint ("                       lognormal(<mean_age>,<st.dev.>)                           \n");
        MrBayesPrint ("                       offsetlognormal(<min_age>,<mean_age>,<st.dev.>)           \n");
        MrBayesPrint ("                       gamma(<mean_age>,<st.dev.>)                               \n");
        MrBayesPrint ("                       offsetgamma(<min_age>,<mean_age>,<st.dev.>)               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    These are the same options used for the 'Calibrate' command. \n");
        MrBayesPrint ("                    Note that, unlike elsewhere in MrMayes, we always use the    \n");
        MrBayesPrint ("                    mean and standard deviation of the resulting age distribution\n");
        MrBayesPrint ("                    rather than the standard parameterization, if different. This\n");
        MrBayesPrint ("                    is to facilitate for the users who want to focus on the in-  \n");
        MrBayesPrint ("                    formation conveyed about the age. For those who wish to use  \n");
        MrBayesPrint ("                    the standard parameterization, there are simple conversions  \n");
        MrBayesPrint ("                    between the two. See the 'Calibrate' command for more infor- \n");
        MrBayesPrint ("                    mation.                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    The tree age is simply the age of the most recent common     \n");
        MrBayesPrint ("                    ancestor of the tree. If the clock rate is fixed to 1.0,     \n");
        MrBayesPrint ("                    which is the default, the tree age is equivalent to the      \n");
        MrBayesPrint ("                    expected number of substitutions from the root to the tip of \n");
        MrBayesPrint ("                    the tree, that is, tree height. The tree age prior ensures   \n");
        MrBayesPrint ("                    that the joint probability for the uniform prior (or fossil- \n");
        MrBayesPrint ("                    ization prior) model of branch lengths on a clock tree is    \n");
        MrBayesPrint ("                    proper. The default setting is 'gamma(1,1)'. If the root node\n");
        MrBayesPrint ("                    in the tree is calibrated, the root calibration replaces the \n");
        MrBayesPrint ("                    tree age prior.                                              \n");
        MrBayesPrint ("   Speciationpr  -- This parameter sets the prior on the net speciation rate (net\n");
        MrBayesPrint ("                    diversification), that is, (lambda - mu) in the birth-death  \n");
        MrBayesPrint ("                    model and the general case of fossilized birth-death model.  \n");
        MrBayesPrint ("                    Or, (lambda - mu - psi) in the special case of f-b-d model   \n");
        MrBayesPrint ("                    (fossiltip). Values of this parameter are > 0. Prior options:\n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset speciationpr = uniform(<number>,<number>)           \n");
        MrBayesPrint ("                       prset speciationpr = exponential(<number>)                \n");
        MrBayesPrint ("                       prset speciationpr = fixed(<number>)                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only relevant if the (fossil) birth-death  \n");
        MrBayesPrint ("                    process is selected as the prior on branch lengths.          \n");
        MrBayesPrint ("   Extinctionpr  -- This parameter sets the prior on the relative extinction rate\n");
        MrBayesPrint ("                    (turnover), that is, (mu / lambda) in the birth-death model  \n");
        MrBayesPrint ("                    and the general case of fossilized birth-death model.        \n");
        MrBayesPrint ("                    Or, (mu + psi) / lambda in the special case of f-b-d model   \n");
        MrBayesPrint ("                    (fossiltip). Values of this parameter are in range (0,1).    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset extinctionpr = beta(<number>,<number>)              \n");
        MrBayesPrint ("                       prset extinctionpr = fixed(<number>)                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only relevant if the (fossil) birth-death  \n");
        MrBayesPrint ("                    process is selected as the prior on branch lengths.          \n");
        MrBayesPrint (" Fossilizationpr -- This parameter sets the prior on the relative fossilization  \n");
        MrBayesPrint ("                    rate (sampling proportion), psi/(mu+psi), in the fossilized  \n");
        MrBayesPrint ("                    b-d model. Values of this parameter are in range (0,1).      \n");
        MrBayesPrint ("                    If SampleStrat is used to divide up time intervals, it sets  \n");
        MrBayesPrint ("                    the prior for the fossilization parameter in each interval.  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset fossilizationpr = beta(<number>,<number>)           \n");
        MrBayesPrint ("                       prset fossilizationpr = fixed(<number>)                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only relevant if the fossilized birth-death\n");
        MrBayesPrint ("                    process is selected as the prior on branch lengths.          \n");
        MrBayesPrint ("   SampleStrat   -- This parameter sets the strategy under which species were    \n");
        MrBayesPrint ("                    sampled in the analysis. For the birth-death prior, 'birth-  \n");
        MrBayesPrint ("                    death' (Hohna et al. 2011), three strategies: 'random',      \n");
        MrBayesPrint ("                    'diversity' and 'cluster' sampling can be used for extant    \n");
        MrBayesPrint ("                    taxa. No extinct sample (fossil) is allowed in this prior.   \n");
        MrBayesPrint ("                    For data with extant and extinct samples, use 'prset brlenspr\n");
        MrBayesPrint ("                    =clock:fossilization'. (Stadler 2010; Zhang et al. 2015)     \n");
        MrBayesPrint ("                    For the fossilized birth-death prior, 'fossiltip' assumes    \n");
        MrBayesPrint ("                    extant taxa are sampled randomly, and extinct taxa (fossils) \n");
        MrBayesPrint ("                    are sampled with constant rate and upon sampling the lineage \n");
        MrBayesPrint ("                    is dead and won't produce any descendant, so fossils are all \n");
        MrBayesPrint ("                    at tips. Except 'fossiltip', the following strategies allow  \n");
        MrBayesPrint ("                    fossils also being ancestors of other samples.               \n");
        MrBayesPrint ("                    'random' (default) assumes extant taxa are sampled uniformly \n");
        MrBayesPrint ("                    at random with prob rho, 'diversity' assumes extant taxa are \n");
        MrBayesPrint ("                    sampled to maximize diversity with prop (set in sampleprob). \n");
        MrBayesPrint ("                    Fossils are sampled on the birth-death tree with piecewise   \n");
        MrBayesPrint ("                    constant rates, psi_i (i = 1,...,s+1). Time is divided by <s>\n");
        MrBayesPrint ("                    slices in the past, each at time <t_i> (s >= 0).             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset samplestrat = random                                \n");
        MrBayesPrint ("                       prset samplestrat = diversity                             \n");
        MrBayesPrint ("                       prset samplestrat = cluster                               \n");
        MrBayesPrint ("                       prset samplestrat = fossiltip                             \n");
        MrBayesPrint ("                       prset samplestrat = random    <s>: ... <t_i> ...          \n");
        MrBayesPrint ("                       prset samplestrat = diversity <s>: ... <t_i> ...          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Sampleprob    -- This parameter sets the fraction of extant species that are  \n");
        MrBayesPrint ("                    sampled in the analysis. This is used with the birth-death   \n");
        MrBayesPrint ("                    prior on trees (Yang and Rannala 1997; Stadler 2009; Hohna   \n");
        MrBayesPrint ("                    et al. 2011), and the fossilized birth-death prior (Stadler  \n");
        MrBayesPrint ("                    2010, Zhang et al. 2015).                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset sampleprob = <number>                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Popsizepr     -- This parameter sets the prior on the population size compo-  \n");
        MrBayesPrint ("                    nent of the coalescent parameter. The options are:           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset popsizepr = uniform(<number>,<number>)              \n");
        MrBayesPrint ("                       prset popsizepr = lognormal(<number>,<number>)            \n");
        MrBayesPrint ("                       prset popsizepr = normal(<number>,<number>)               \n");
        MrBayesPrint ("                       prset popsizepr = gamma(<number>,<number>)                \n");
        MrBayesPrint ("                       prset popsizepr = fixed(<number>)                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only relevant if the coalescence process is\n");
        MrBayesPrint ("                    selected as the prior on branch lengths. Note that the set-  \n");
        MrBayesPrint ("                    ting of 'ploidy' in 'lset' is important for how this para-   \n");
        MrBayesPrint ("                    meter is interpreted.                                        \n");
        MrBayesPrint ("   Popvarpr      -- In a gene tree - species tree model, this parameter deter-   \n");
        MrBayesPrint ("                    mines whether the population size is the same for the entire \n");
        MrBayesPrint ("                    species tree ('popvarpr = equal', the default), or varies    \n");
        MrBayesPrint ("                    across branches of the species tree ('popvarpr=variable').   \n");
/*      MrBayesPrint ("   Growthpr      -- This parameter sets the prior on the exponential growth      \n");
        MrBayesPrint ("                    parameter of the coalescence process. The options are:       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset growthpr = uniform(<number>,<number>)               \n");
        MrBayesPrint ("                       prset growthpr = exponential(<number>)                    \n");
        MrBayesPrint ("                       prset growthpr = fixed(<number>)                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    This parameter is only relevant if the coalescence           \n");
        MrBayesPrint ("                    process is selected as the prior on branch lengths.          \n"); */
        MrBayesPrint ("   Nodeagepr     -- This parameter specifies the assumptions concerning the age  \n");
        MrBayesPrint ("                    of the terminal and interior nodes in the tree. The default  \n");
        MrBayesPrint ("                    model ('nodeagepr = unconstrained') assumes that all terminal\n");
        MrBayesPrint ("                    nodes are of the same age while the age of interior nodes is \n");
        MrBayesPrint ("                    unconstrained. The alternative ('nodeagepr = calibrated')    \n");
        MrBayesPrint ("                    option derives a prior probability distribution on terminal  \n");
        MrBayesPrint ("                    and interior node ages from the calibration settings (see    \n");
        MrBayesPrint ("                    the 'calibrate' command). The 'nodeagepr' parameter is only  \n");
        MrBayesPrint ("                    relevant for clock trees.                                    \n");
        MrBayesPrint ("   Clockratepr   -- This parameter specifies the prior assumptions concerning the\n");
        MrBayesPrint ("                    base substitution rate of the tree, measured in expected num-\n");
        MrBayesPrint ("                    ber of substitutions per site per time unit. The default set-\n");
        MrBayesPrint ("                    ting is 'Fixed(1.0)', which effectively means that the time  \n");
        MrBayesPrint ("                    unit is the number of expected substitutions per site.       \n");
/*      MrBayesPrint ("                    If you apply age constraints to the tree, the default setting\n");
        MrBayesPrint ("                    changes automatically to 'Exponential(<x>)', where '<x>' (the\n");
        MrBayesPrint ("                    rate of exponential) is ten times the age of the maximum age \n");
        MrBayesPrint ("                    constraint. This will give you a very vague prior, which may \n");
        MrBayesPrint ("                    or may not be adequate for your particular problem.          \n"); */
        MrBayesPrint ("                    If you do not have any age calibrations in the tree, you can \n");
        MrBayesPrint ("                    still calibrate the tree using 'Clockratepr'. For instance,  \n");
        MrBayesPrint ("                    if you know that your sequence data evolve at a rate of 0.20 \n");
        MrBayesPrint ("                    substitutions per million years, you might calibrate the tree\n");
        MrBayesPrint ("                    by fixing the substitution rate to 0.20 using                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset clockratepr = fixed(0.20)                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    after which the tree will be calibrated using millions of    \n");
        MrBayesPrint ("                    years as the unit.                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    You can also assign a prior probability distribution to the  \n");
        MrBayesPrint ("                    substitution rate, accommodating the uncertainty of it.      \n");
        MrBayesPrint ("                    When you calibrate the nodes, you should properly set this   \n");
        MrBayesPrint ("                    prior to match the time unit of the calibrations.            \n");
        MrBayesPrint ("                    You can choose among normal, lognormal, exponential and gamma\n");
        MrBayesPrint ("                    distributions for this purpose. For instance, to assign a    \n");
        MrBayesPrint ("                    normal distribution truncated at 0, so that only positive    \n");
        MrBayesPrint ("                    values are allowed, and with mean 0.20 and standard deviation\n");
        MrBayesPrint ("                    of 0.02, you would use                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset clockratepr = normal(0.20,0.02)                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    The lognormal distribution is parameterized in terms of the  \n");
        MrBayesPrint ("                    mean and standard deviation on the log scale (natural logs). \n");
        MrBayesPrint ("                    For instance,                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset clockratepr = lognormal(-1.61,0.10)                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    specifies a lognormal distribution with a mean of log values \n");
        MrBayesPrint ("                    of -1.61 and a standard deviation of log values of 0.10. In  \n");
        MrBayesPrint ("                    such a case, the mean value of the lognormal distribution is \n");
        MrBayesPrint ("                    equal to e^(-1.61 + 0.10^2/2) = 0.20.                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    Note that the 'Clockratepr' parameter has no effect on non-  \n");
        MrBayesPrint ("                    clock trees.                                                 \n");
        MrBayesPrint ("   Clockvarpr    -- This parameter allows you to specify the type of clock you   \n");
        MrBayesPrint ("                    are assuming. The default is 'strict', which corresponds to  \n");
        MrBayesPrint ("                    the standard clock model where the evolutionary rate is      \n");
        MrBayesPrint ("                    constant throughout the tree. For relaxed clock models, you  \n");
        MrBayesPrint ("                    can use 'cpp', 'tk02', 'igr'. ('mixed' is not working)       \n");
        MrBayesPrint ("                    'cpp' invokes a relaxed clock model where the rate evolves   \n");
        MrBayesPrint ("                    according to a Compound Poisson Process (CPP) (Huelsenbeck   \n");
        MrBayesPrint ("                    et al., 2000).                                               \n");
        MrBayesPrint ("                    'tk02' invokes the Brownian Motion model described by Thorne \n");
        MrBayesPrint ("                    and Kishino (2002). [autocorrelated lognormal distributions] \n");
        MrBayesPrint ("                    'igr' invokes the Independent Gamma Rate (IGR) model where   \n");
        MrBayesPrint ("                    each branch has an independent rate drawn from a gamma       \n");
        MrBayesPrint ("                    distribution (LePage et al., 2007).                          \n");
        MrBayesPrint ("                    Each of the relaxed clock models has additional parameters   \n");
        MrBayesPrint ("                    with priors. For the CPP model, it is 'cppratepr' and        \n");
        MrBayesPrint ("                    'cppmultdevpr'; for the TK02 model, it is 'tk02varpr'; for   \n");
        MrBayesPrint ("                    the IGR  model, it is 'igrvarpr'.                            \n");
        MrBayesPrint ("                    The 'clockvarpr' parameter is only relevant for clock trees. \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    For backward compatibility, 'bm' is allowed as a synonym of  \n");
        MrBayesPrint ("                    'tk02', and 'ibr' as a synonym of 'igr'.                     \n");
        MrBayesPrint ("   Cppratepr     -- This parameter allows you to specify a prior probability     \n");
        MrBayesPrint ("                    distribution on the rate of the Poisson process generating   \n");
        MrBayesPrint ("                    changes in the evolutionary rate in the CPP relaxed clock    \n");
        MrBayesPrint ("                    model. You can either fix the rate or associate it with an   \n");
        MrBayesPrint ("                    exponential prior using                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset cppratepr = fixed(<number>)                         \n");
        MrBayesPrint ("                       prset cppratepr = exponential(<number>)                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    For instance, if you fix the rate to 2, then on a branch     \n");
        MrBayesPrint ("                    with the length equual to one expresed in terms of average   \n");
        MrBayesPrint ("                    expected number of substitution per site, you expect to see, \n"); 
        MrBayesPrint ("                    on average, two rate-modifying events.                       \n");
        MrBayesPrint ("                    If you put an exponential(0.1) on the rate, you will be      \n");
        MrBayesPrint ("                    estimating the rate against a prior probability distribution \n");
        MrBayesPrint ("                    where the expected rate is 10 (= 1/0.1).                     \n");
        MrBayesPrint ("   Cppmultdevpr  -- This parameter allows you to specify the standard deviation  \n");
        MrBayesPrint ("                    of the log-normal distribution from which the rate multi-    \n");
        MrBayesPrint ("                    pliers of the CPP relaxed clock model are drawn. The standard\n");
        MrBayesPrint ("                    deviation is given on the log scale. The default value of 1.0\n");
        MrBayesPrint ("                    thus corresponds to rate multipliers varying from 0.37 (1/e) \n");
        MrBayesPrint ("                    to 2.7 (e) when they are +/- one standard deviation from the \n");
        MrBayesPrint ("                    expected mean. The expected mean of the logarithm of the mul-\n");
        MrBayesPrint ("                    pliers is fixed to 0, ensuring that the expected mean rate is\n");
        MrBayesPrint ("                    1.0. You can change the default value by using               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset cppmultdevpr = fixed(<number>)                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    where <number> is the standard deviation on the log scale.   \n");
        MrBayesPrint ("   TK02varpr     -- This parameter allows you to specify the prior probability   \n");
        MrBayesPrint ("                    distribution for the variance of the rate multiplier in the  \n");
        MrBayesPrint ("                    Thorne-Kishino ('Brownian motion') relaxed clock model.      \n");
        MrBayesPrint ("                    Specifically, the parameter specifies the rate at which the  \n");
        MrBayesPrint ("                    variance increases with respect to the base rate of the      \n");
        MrBayesPrint ("                    clock. If you have a branch of a length corresponding to 0.4 \n");
        MrBayesPrint ("                    expected changes per site according to the base rate of the  \n");
        MrBayesPrint ("                    clock, and the tk02var parameter has a value of 2.0, then the\n");
        MrBayesPrint ("                    rate multiplier at the end of the branch will be drawn from a\n");
        MrBayesPrint ("                    lognormal distribution with a variance of 0.4*2.0 (on the    \n");
        MrBayesPrint ("                    linear, not the logarithm scale). The mean is the same as the\n");
        MrBayesPrint ("                    rate multiplier at the start of the branch (again on the     \n");
        MrBayesPrint ("                    linear scale).                                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    You can set the parameter to a fixed value, or specify that  \n");
        MrBayesPrint ("                    it is drawn from an exponential or uniform distribution:     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset tk02varpr = fixed(<number>)                         \n");
        MrBayesPrint ("                       prset tk02varpr = exponential(<number>)                   \n");
        MrBayesPrint ("                       prset tk02varpr = uniform(<number>,<number>)              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    For backward compatibility, 'bmvarpr' is allowed as a synonym\n");
        MrBayesPrint ("                    of 'tko2varpr'.                                              \n");
        MrBayesPrint ("   Igrvarpr      -- This parameter allows you to specify a prior on the variance \n");
        MrBayesPrint ("                    of the gamma distribution from which the branch lengths are  \n");
        MrBayesPrint ("                    drawn in the independent branch rate (IGR) relaxed clock     \n");
        MrBayesPrint ("                    model. Specifically, the parameter specifies the rate at     \n");
        MrBayesPrint ("                    which the variance increases with respect to the base rate of\n");
        MrBayesPrint ("                    the clock. If you have a branch of a length corresponding to \n");
        MrBayesPrint ("                    0.4 expected changes per site according to the base rate of  \n");
        MrBayesPrint ("                    the clock, and the igrvar parameter has a value of 2.0, then \n");
        MrBayesPrint ("                    the effective branch length will be drawn from a distribution\n");
        MrBayesPrint ("                    with a variance of 0.4*2.0.                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    You can set the parameter to a fixed value, or specify that  \n");
        MrBayesPrint ("                    it is drawn from an exponential or uniform distribution:     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset igrvarpr = fixed(<number>)                          \n");
        MrBayesPrint ("                       prset igrvarpr = exponential(<number>)                    \n");
        MrBayesPrint ("                       prset igrvarpr = uniform(<number>,<number>)               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    For backward compatibility, 'ibrvarpr' is allowed as a syn-  \n");
        MrBayesPrint ("                    onym of 'igrvarpr'.                                          \n");
        MrBayesPrint ("   Ratepr        -- This parameter allows you to specify the site specific rates \n");
        MrBayesPrint ("                    model or any other model that allows different partitions to \n");
        MrBayesPrint ("                    evolve at different rates. First, you must have defined a    \n");
        MrBayesPrint ("                    partition of the characters. For example, you may define a   \n");
        MrBayesPrint ("                    partition that divides the characters by codon position, if  \n");
        MrBayesPrint ("                    you have DNA data. You can also divide your data using a     \n");
        MrBayesPrint ("                    partition that separates different genes from each other.    \n");
        MrBayesPrint ("                    The next step is to make the desired partition the active one\n");
        MrBayesPrint ("                    using the set command. For example, if your partition is     \n");
        MrBayesPrint ("                    called \"by_codon\", then you make that the active partition \n");
        MrBayesPrint ("                    using \"set partition=by_codon\". Now that you have defined  \n");
        MrBayesPrint ("                    and activated a partition, you can specify the rate multi-   \n");
        MrBayesPrint ("                    pliers for the various partitions. The options are:          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                       prset ratepr = fixed                                      \n");
        MrBayesPrint ("                       prset ratepr = variable                                   \n");
        MrBayesPrint ("                       prset ratepr = dirichlet(<number>,<number>,...,<number>)  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                    If you specify \"fixed\", then the rate multiplier for       \n");
        MrBayesPrint ("                    that partition is set to 1 (i.e., the rate is fixed to       \n");
        MrBayesPrint ("                    the average rate across partitions). On the other hand,      \n");
        MrBayesPrint ("                    if you specify \"variable\", then the rate is allowed to     \n");
        MrBayesPrint ("                    vary across partitions subject to the constraint that the    \n");
        MrBayesPrint ("                    average rate of substitution across the partitions is 1.     \n");
        MrBayesPrint ("                    You must specify a variable rate prior for at least two      \n");
        MrBayesPrint ("                    partitions, otherwise the option is not activated when       \n");
        MrBayesPrint ("                    calculating likelihoods. The variable option automatically   \n");
        MrBayesPrint ("                    associates the partition rates with a dirichlet(1,...,1)     \n");
        MrBayesPrint ("                    prior. The dirichlet option is an alternative way of setting \n");
        MrBayesPrint ("                    a partition rate to be variable, and also gives accurate     \n");
        MrBayesPrint ("                    control of the shape of the prior. The parameters of the     \n");
        MrBayesPrint ("                    Dirichlet are listed in the order of the partitions that the \n");
        MrBayesPrint ("                    ratepr is applied to. For instance, \"prset applyto=(1,3,4)  \n");
        MrBayesPrint ("                    ratepr = dirichlet(10,40,15)\" would set the Dirichlet para- \n");
        MrBayesPrint ("                    meter 10 to partition 1, 40 to partition 3, and 15 to parti- \n");
        MrBayesPrint ("                    tion 4. The Dirichlet distribution is applied to the weighted\n");
        MrBayesPrint ("                    rates; that is, it weights the partition rates according to  \n");
        MrBayesPrint ("                    the number of included characters in each partition.         \n");
        MrBayesPrint ("   Generatepr    -- This parameter is similar to 'Ratepr' but applies to gene    \n");
        MrBayesPrint ("                    trees in the multispecies coalescent, whereas 'Ratepr' app-  \n");
        MrBayesPrint ("                    lies to partitions within genes.                             \n");
        MrBayesPrint ("                                                                                 \n");
        if (numCurrentDivisions == 0)
            tempInt = 1;
        else
            tempInt = numCurrentDivisions;
        for (i=0; i<tempInt; i++)
            {
            if (numCurrentDivisions == 0)
                {
                MrBayesPrint ("   Default model settings:                                                       \n");
                mp = &defaultModel;
                }
            else
                {
                MrBayesPrint ("   Model settings for partition %d:                                              \n", i+1);
                mp = &modelParams[i];
                }
            MrBayesPrint ("                                                                                 \n");
            MrBayesPrint ("   Parameter        Options                      Current Setting                 \n");
            MrBayesPrint ("   ------------------------------------------------------------------            \n");

            MrBayesPrint ("   Tratiopr         Beta/Fixed                   %s", mp->tRatioPr);
            if (!strcmp(mp->tRatioPr, "Beta"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->tRatioDir[0], mp->tRatioDir[1]);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->tRatioFix);

            MrBayesPrint ("   Revmatpr         Dirichlet/Fixed              %s", mp->revMatPr);
            if (!strcmp(mp->revMatPr, "Dirichlet"))
                MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", mp->revMatDir[0],
                mp->revMatDir[1], mp->revMatDir[2], mp->revMatDir[3],
                mp->revMatDir[4], mp->revMatDir[5]);
            else
                MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", mp->revMatFix[0],
                mp->revMatFix[1], mp->revMatFix[2], mp->revMatFix[3],
                mp->revMatFix[4], mp->revMatFix[5]);

            MrBayesPrint ("   Aamodelpr        Fixed/Mixed                  %s", mp->aaModelPr);
            if (!strcmp(mp->aaModelPr, "Fixed"))
                MrBayesPrint ("(%s)\n", mp->aaModel);
            else
                MrBayesPrint ("\n");

            MrBayesPrint ("   Aarevmatpr       Dirichlet/Fixed              %s", mp->aaRevMatPr);
            if (!strcmp(mp->aaRevMatPr, "Dirichlet"))
                {
                for (j=1; j<190; j++)
                    if (AreDoublesEqual (mp->aaRevMatDir[0], mp->aaRevMatDir[j], 0.00001) == NO)
                        break;
                if (j==190)
                    MrBayesPrint ("(%1.1lf,%1.1lf,...)\n", mp->aaRevMatDir[0], mp->aaRevMatDir[0]);
                else
                    MrBayesPrint (" (use 'Showmodel' to see values set by user)\n");
                }
            else
                {
                for (j=1; j<190; j++)
                    if (AreDoublesEqual (mp->aaRevMatFix[0], mp->aaRevMatFix[j], 0.00001) == NO)
                        break;
                if (j==190)
                    MrBayesPrint ("(%1.1lf,%1.1lf,...)\n", mp->aaRevMatFix[0], mp->aaRevMatFix[0]);
                else
                    MrBayesPrint (" (use 'Showmodel' to see values set by user)\n");
                }

            MrBayesPrint ("   Omegapr          Dirichlet/Fixed              %s", mp->omegaPr);
            if (!strcmp(mp->omegaPr, "Dirichlet"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->omegaDir[0], mp->omegaDir[1]);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->omegaFix);

            MrBayesPrint ("   Ny98omega1pr     Beta/Fixed                   %s", mp->ny98omega1pr);
            if (!strcmp(mp->ny98omega1pr, "Beta"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->ny98omega1Beta[0], mp->ny98omega1Beta[1]);
            else if (!strcmp(mp->ny98omega1pr, "Fixed"))
                MrBayesPrint ("(%1.1lf)\n", mp->ny98omega1Fixed);
                
            MrBayesPrint ("   Ny98omega3pr     Uniform/Exponential/Fixed    %s", mp->ny98omega3pr);
            if (!strcmp(mp->ny98omega3pr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->ny98omega3Uni[0], mp->ny98omega3Uni[1]);
            else if (!strcmp(mp->ny98omega3pr, "Exponential"))
                MrBayesPrint ("(%1.1lf)\n", mp->ny98omega3Exp);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->ny98omega3Fixed);

            MrBayesPrint ("   M3omegapr        Exponential/Fixed            %s", mp->m3omegapr);
            if (!strcmp(mp->m3omegapr, "Exponential"))
                MrBayesPrint ("\n");
            else if (!strcmp(mp->m3omegapr, "Fixed"))
                MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", mp->m3omegaFixed[0], mp->m3omegaFixed[1], mp->m3omegaFixed[2]);
                
            MrBayesPrint ("   Codoncatfreqs    Dirichlet/Fixed              %s", mp->codonCatFreqPr);
            if (!strcmp(mp->codonCatFreqPr, "Dirichlet"))
                MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", mp->codonCatDir[0], mp->codonCatDir[1], mp->codonCatDir[2]);
            else
                MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf)\n", mp->codonCatFreqFix[0], mp->codonCatFreqFix[1], mp->codonCatFreqFix[2]);

            MrBayesPrint ("   Statefreqpr      Dirichlet/Fixed              %s", mp->stateFreqPr);
            if (!strcmp(mp->stateFreqPr, "Dirichlet"))
                {
                if (mp->dataType == DNA || mp->dataType == RNA)
                    {
                    if (!strcmp(mp->nucModel, "4by4"))
                        MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", mp->stateFreqsDir[0], mp->stateFreqsDir[1],
                            mp->stateFreqsDir[2], mp->stateFreqsDir[3]);
                    else
                        MrBayesPrint ("\n");
                    }
                else if (mp->dataType == RESTRICTION)
                    {
                    MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->stateFreqsDir[0], mp->stateFreqsDir[1]);
                    }
                else
                    MrBayesPrint ("\n");
                }
            else if (!strcmp(mp->stateFreqPr, "Fixed"))
                {
                if (mp->dataType == DNA || mp->dataType == RNA)
                    {
                    if (!strcmp(mp->nucModel, "4by4"))
                        MrBayesPrint ("(%1.1lf,%1.1lf,%1.1lf,%1.1lf)\n", mp->stateFreqsFix[0], mp->stateFreqsFix[1],
                            mp->stateFreqsFix[2], mp->stateFreqsFix[3]);
                    else
                        MrBayesPrint ("\n");
                    }
                else if (mp->dataType == RESTRICTION)
                    {
                    MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->stateFreqsFix[0], mp->stateFreqsFix[1]);
                    }
                else
                    MrBayesPrint ("\n");
                }

            MrBayesPrint ("   Shapepr          Uniform/Exponential/Fixed    %s", mp->shapePr);
            if (!strcmp(mp->shapePr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->shapeUni[0], mp->shapeUni[1]);
            else if (!strcmp(mp->shapePr, "Exponential"))
                MrBayesPrint ("(%1.1lf)\n", mp->shapeExp);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->shapeFix);

            MrBayesPrint ("   Ratecorrpr       Uniform/Fixed                %s", mp->adGammaCorPr);
            if (!strcmp(mp->adGammaCorPr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->corrUni[0], mp->corrUni[1]);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->corrFix);

            MrBayesPrint ("   Pinvarpr         Uniform/Fixed                %s", mp->pInvarPr);
            if (!strcmp(mp->pInvarPr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->pInvarUni[0], mp->pInvarUni[1]);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->pInvarFix);

            MrBayesPrint ("   Covswitchpr      Uniform/Exponential/Fixed    %s", mp->covSwitchPr);
            if (!strcmp(mp->covSwitchPr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->covswitchUni[0], mp->covswitchUni[1]);
            else if (!strcmp(mp->covSwitchPr, "Exponential"))
                MrBayesPrint ("(%1.1lf)\n", mp->covswitchExp);
            else
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->covswitchFix[0], mp->covswitchFix[1]);

            MrBayesPrint ("   Symdirihyperpr   Uniform/Exponential/Fixed    %s", mp->symPiPr);
            if (!strcmp(mp->symPiPr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->symBetaUni[0], mp->symBetaUni[1]);
            else if (!strcmp(mp->covSwitchPr, "Exponential"))
                MrBayesPrint ("(%1.1lf)\n", mp->symBetaExp);
            else
                {
                if (mp->symBetaFix < 0)
                    MrBayesPrint ("(Infinity)\n");
                else
                    MrBayesPrint ("(%1.1lf)\n", mp->symBetaFix);
                }
            
            MrBayesPrint ("   Topologypr       Uniform/Constraints/Fixed/   %s", mp->topologyPr);
            if (!strcmp(mp->topologyPr, "Constraints"))
                {
                MrBayesPrint ("(");
                for (j=0; j<numDefinedConstraints; j++)
                    {
                    if (mp->activeConstraints[j] == YES)
                        {
                        MrBayesPrint ("%d", j+1);
                        break;
                        }
                    }
               for (j++; j<numDefinedConstraints; j++)
                    {
                    if (mp->activeConstraints[j] == YES)
                        {
                        MrBayesPrint (",%d", j+1);
                        }
                    }
                MrBayesPrint (")\n");
                }
            else if (!strcmp(mp->topologyPr, "Fixed"))
                MrBayesPrint("(%s)\n", userTree[mp->topologyFix]->name);
            else
                MrBayesPrint ("\n");
            MrBayesPrint ("                    Speciestree                  \n");
            
            MrBayesPrint ("   Brlenspr         Unconstrained/Clock/Fixed    %s", mp->brlensPr);
            if (!strcmp(mp->brlensPr, "Unconstrained"))
                {
                if (!strcmp(mp->unconstrainedPr, "Uniform"))
                    MrBayesPrint (":Uni(%1.1lf,%1.1lf)\n", mp->brlensUni[0], mp->brlensUni[1]);
                else if (!strcmp(mp->unconstrainedPr, "GammaDir"))
                    MrBayesPrint (":GammaDir(%1.1lf,%1.3lf,%1.1lf,%1.1lf)\n",
                                mp->brlensDir[0], mp->brlensDir[1], mp->brlensDir[2], mp->brlensDir[3]);
                else if (!strcmp(mp->unconstrainedPr, "invGamDir"))
                    MrBayesPrint (":invGamDir(%1.1lf,%1.3lf,%1.1lf,%1.1lf)\n",
                                mp->brlensDir[0], mp->brlensDir[1], mp->brlensDir[2], mp->brlensDir[3]);
                else if (!strcmp(mp->unconstrainedPr, "twoExp"))
                    MrBayesPrint (":twoExp(%1.1lf,%1.1lf)\n", mp->brlens2Exp[0], mp->brlens2Exp[1]);
                else
                    MrBayesPrint (":Exp(%1.1lf)\n", mp->brlensExp);
                }
            else if (!strcmp(mp->brlensPr, "Clock"))
                {
                if (!strcmp(mp->clockPr,"Fixed"))
                    MrBayesPrint (":%s(%s)\n", mp->clockPr, userTree[mp->brlensFix]->name);
                else
                    MrBayesPrint (":%s\n", mp->clockPr);
                }
            else if (!strcmp(mp->brlensPr, "Fixed"))
                MrBayesPrint("(%s)\n", userTree[mp->brlensFix]->name);
            
            MrBayesPrint ("   Treeagepr        Gamma/Uniform/Fixed/         %s\n", mp->treeAgePr.name);
            MrBayesPrint ("                    Truncatednormal/Lognormal/   \n");
            MrBayesPrint ("                    Offsetlognormal/Offsetgamma/ \n");
            MrBayesPrint ("                    Offsetexponential            \n");
            
            MrBayesPrint ("   Speciationpr     Uniform/Exponential/Fixed    %s", mp->speciationPr);
            if (!strcmp(mp->speciationPr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->speciationUni[0], mp->speciationUni[1]);
            else if (!strcmp(mp->speciationPr, "Exponential"))
                MrBayesPrint ("(%1.1lf)\n", mp->speciationExp);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->speciationFix);
            
            MrBayesPrint ("   Extinctionpr     Beta/Fixed                   %s", mp->extinctionPr);
            if (!strcmp(mp->extinctionPr, "Beta"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->extinctionBeta[0], mp->extinctionBeta[1]);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->extinctionFix);
            
            MrBayesPrint ("   Fossilizationpr  Beta/Fixed                   %s", mp->fossilizationPr);
            if (!strcmp(mp->fossilizationPr, "Beta"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->fossilizationBeta[0], mp->fossilizationBeta[1]);
            else
                MrBayesPrint ("(%1.2lf)\n", mp->fossilizationFix);
            
            MrBayesPrint ("   SampleStrat      Random/Diversity/Cluster/    %s\n", mp->sampleStrat);
            MrBayesPrint ("                    FossilTip                    \n");
            // if (!strcmp(mp->sampleStrat, "Random") || !strcmp(mp->sampleStrat, "Diversity"))
            
            MrBayesPrint ("   Sampleprob       <number>                     %1.8lf\n", mp->sampleProb);
            
            MrBayesPrint ("   Popsizepr        Lognormal/Gamma/Uniform/     %s", mp->popSizePr);
            if (!strcmp(mp->popSizePr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->popSizeUni[0], mp->popSizeUni[1]);
            else if (!strcmp(mp->popSizePr, "Lognormal"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->popSizeLognormal[0], mp->popSizeLognormal[1]);
            else if (!strcmp(mp->popSizePr, "Normal"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->popSizeNormal[0], mp->popSizeNormal[1]);
            else if (!strcmp(mp->popSizePr, "Gamma"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->popSizeGamma[0], mp->popSizeGamma[1]);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->popSizeFix);
            MrBayesPrint ("                    Normal/Fixed                 \n");

            MrBayesPrint ("   Popvarpr         Equal/Variable               %s\n", mp->popVarPr);

            /*
            MrBayesPrint ("   Growthpr         Uniform/Exponential/         \n");
            MrBayesPrint ("                    Fixed/Normal                 %s", mp->growthPr);
            if (!strcmp(mp->growthPr, "Uniform"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->growthUni[0], mp->growthUni[1]);
            else if (!strcmp(mp->growthPr, "Exponential"))
                MrBayesPrint ("(%1.1lf)\n", mp->growthExp);
            else if (!strcmp(mp->growthPr, "Normal"))
                MrBayesPrint ("(%1.1lf,%1.1lf)\n", mp->growthNorm[0], mp->growthNorm[1]);
            else
                MrBayesPrint ("(%1.1lf)\n", mp->growthFix); 
            */

            MrBayesPrint ("   Nodeagepr        Unconstrained/Calibrated     %s\n", mp->nodeAgePr);

            MrBayesPrint ("   Clockratepr      Fixed/Normal/Lognormal/      %s", mp->clockRatePr);
            if (!strcmp(mp->clockRatePr, "Fixed"))
                MrBayesPrint ("(%1.2lf)\n", mp->clockRateFix);
            else if (!strcmp(mp->clockRatePr,"Exponential"))
                MrBayesPrint ("(%1.2lf)\n", mp->clockRateExp);
            else if (!strcmp(mp->clockRatePr,"Normal"))
                MrBayesPrint ("(%1.2lf,%1.2lf)\n", mp->clockRateNormal[0], mp->clockRateNormal[1]);
            else if (!strcmp(mp->clockRatePr,"Lognormal"))
                MrBayesPrint ("(%1.2lf,%1.2lf)\n", mp->clockRateLognormal[0], mp->clockRateLognormal[1]);
            else
                {
                assert (!strcmp(mp->clockRatePr,"Gamma"));
                MrBayesPrint ("(%1.2lf,%1.2lf)\n", mp->clockRateGamma[0], mp->clockRateGamma[1]);
                }
            MrBayesPrint ("                    Exponential/Gamma            \n");

            MrBayesPrint ("   Clockvarpr       Strict/Cpp/TK02/Igr/Mixed    %s\n", mp->clockVarPr);

            MrBayesPrint ("   Cppratepr        Fixed/Exponential            %s", mp->cppRatePr);
            if (!strcmp(mp->cppRatePr, "Fixed"))
                MrBayesPrint ("(%1.2lf)\n", mp->cppRateFix);
            else /* if (!strcmp(mp->cppRatePr,"Exponential")) */
                MrBayesPrint ("(%1.2lf)\n", mp->cppRateExp);

            MrBayesPrint ("   Cppmultdevpr     Fixed                        %s", mp->cppMultDevPr);
            MrBayesPrint ("(%1.2lf)\n", mp->cppMultDevFix);

            MrBayesPrint ("   TK02varpr        Fixed/Exponential/Uniform    %s", mp->tk02varPr);
            if (!strcmp(mp->tk02varPr, "Fixed"))
                MrBayesPrint ("(%1.2lf)\n", mp->tk02varFix);
            else if (!strcmp(mp->tk02varPr,"Exponential"))
                MrBayesPrint ("(%1.2lf)\n", mp->tk02varExp);
            else
                {
                assert (!strcmp(mp->tk02varPr,"Uniform"));
                MrBayesPrint ("(%1.2lf,%1.2lf)\n", mp->tk02varUni[0], mp->tk02varUni[1]);
                }

            MrBayesPrint ("   Igrvarpr         Fixed/Exponential/Uniform    %s", mp->igrvarPr);
            if (!strcmp(mp->igrvarPr, "Fixed"))
                MrBayesPrint ("(%1.2lf)\n", mp->igrvarFix);
            else if (!strcmp(mp->igrvarPr,"Exponential"))
                MrBayesPrint ("(%1.2lf)\n", mp->igrvarExp);
            else
                {
                assert (!strcmp(mp->igrvarPr,"Uniform"));
                MrBayesPrint ("(%1.2lf,%1.2lf)\n", mp->igrvarUni[0], mp->igrvarUni[1]);
                }
            
            /*  MrBayesPrint ("   Mixedvarpr       Fixed/Exponential/Uniform    %s", mp->mixedvarPr);
            if (!strcmp(mp->mixedvarPr, "Fixed"))
                MrBayesPrint ("(%1.2lf)\n", mp->mixedvarFix);
            else if (!strcmp(mp->mixedvarPr,"Exponential"))
                MrBayesPrint ("(%1.2lf)\n", mp->mixedvarExp);
            else
                {
                assert (!strcmp(mp->mixedvarPr,"Uniform"));
                MrBayesPrint ("(%1.2lf,%1.2lf)\n", mp->mixedvarUni[0], mp->mixedvarUni[1]);
                }  */

            MrBayesPrint ("   Ratepr           Fixed/Variable=Dirichlet     %s", mp->ratePr);
            if (!strcmp(mp->ratePr, "Dirichlet"))
                MrBayesPrint ("(...,%1.1lf,...)\n", mp->ratePrDir);
            else
                MrBayesPrint ("\n");

            MrBayesPrint ("   Generatepr       Fixed/Variable=Dirichlet     %s", mp->generatePr);
            if (!strcmp(mp->generatePr, "Dirichlet"))
                MrBayesPrint ("(...,%1.1lf,...)\n", mp->generatePrDir);
            else
                MrBayesPrint ("\n");

            MrBayesPrint ("   ------------------------------------------------------------------            \n");
            MrBayesPrint ("                                                                                 \n");
            }
        }
    else if (!strcmp(helpTkn, "Ctype"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Ctype                                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command sets the character ordering for standard-type data. The          \n");
        MrBayesPrint ("   correct usage is:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      ctype <ordering>:<characters>                                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The available options for the <ordering> specifier are:                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("     unordered    -- Movement directly from one state to another is              \n");
        MrBayesPrint ("                     allowed in an instant of time.                              \n");
        MrBayesPrint ("     ordered      -- Movement is only allowed between adjacent characters.       \n");
        MrBayesPrint ("                     For example, perhaps only between 0 <-> 1 and 1 <-> 2       \n");
        MrBayesPrint ("                     for a three state character ordered as 0 - 1 - 2.           \n");
        MrBayesPrint ("     irreversible -- Rates of change for losses are 0.                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The characters to which the ordering is applied is specified in manner        \n");
        MrBayesPrint ("   that is identical to commands such as \"include\" or \"exclude\". For         \n");
        MrBayesPrint ("   example,                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      ctype ordered: 10 23 45                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   defines charactes 10, 23, and 45 to be of type ordered. Similarly,            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      ctype irreversible: 54 - 67  71-92                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   defines characters 54 to 67 and characters 71 to 92 to be of type             \n");
        MrBayesPrint ("   irreversible. You can use the \".\" to denote the last character, and         \n");
        MrBayesPrint ("   \"all\" to denote all of the characters. Finally, you can use the             \n");
        MrBayesPrint ("   specifier \"\\\" to apply the ordering to every n-th character or             \n");
        MrBayesPrint ("   you can use predefined charsets to specify the character.                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Only one ordering can be used on any specific application of ctype.           \n");
        MrBayesPrint ("   If you want to apply different orderings to different characters, then        \n");
        MrBayesPrint ("   you need to use ctype multiple times. For example,                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      ctype ordered: 1-50                                                        \n");
        MrBayesPrint ("      ctype irreversible: 51-100                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   sets characters 1 to 50 to be ordered and characters 51 to 100 to be          \n");
        MrBayesPrint ("   irreversible.                                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The ctype command is only sensible with morphological (here called            \n");
        MrBayesPrint ("   \"standard\") characters. The program ignores attempts to apply char-         \n");
        MrBayesPrint ("   acter orderings to other types of characters, such as DNA characters.         \n");

        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Propset"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Propset                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command allows the user to change the details of the MCMC samplers       \n");
        MrBayesPrint ("   (moves) that update the state of the chain. The useage is:                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      propset  <move_name>$<tuning-parameter>=<value>                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Assume we have a topology parameter called 'Tau{all}', which is sampled by    \n");
        MrBayesPrint ("   the move 'ExtTBR(Tau{all})' (note that the parameter name is included in the  \n");
        MrBayesPrint ("   move name). This move has three tuning parameters: (1) 'prob', the relative   \n");
        MrBayesPrint ("   proposal probability (a weight defining its probability relative to other     \n");
        MrBayesPrint ("   moves); (2) 'p_ext', the extension probability; and (3) 'lambda', the tuning  \n");
        MrBayesPrint ("   parameter of the branch length multiplier. A list of the tuning parameters is \n");
        MrBayesPrint ("   available by using 'Showmoves' (see below). To change the relative proposal   \n");
        MrBayesPrint ("   probability to 20 and the extension probability to 0.7, use:                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      propset etbr(tau{all})$prob=20 etbr(tau{all})$p_ext=0.7                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This change would apply to all chains in all runs. It is also possible to set \n");
        MrBayesPrint ("   the tuning parameters of individual runs and chains using the format:         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      propset  <move_name>$<tuning-parameter>(<run>,<chain>)=<value>             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   where <run> and <chain> are the index numbers of the run and chain for which  \n");
        MrBayesPrint ("   you want to change the value. If you leave out the index of the run, the      \n");
        MrBayesPrint ("   change will apply to all runs; if you leave out the index of the chain, the   \n");
        MrBayesPrint ("   change will similarly apply to all chains. To switch off the exttbr(tau{all}) \n");
        MrBayesPrint ("   move in chain 2 of all runs, use:                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      propset  etbr(tau{all})$prob(,2)=0                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   It is important to note that all moves are not available until the model has  \n");
        MrBayesPrint ("   been completely defined. Any change to the model will cause all proposal      \n");
        MrBayesPrint ("   tuning parameters to return to their default values. To see a list of all the \n");
        MrBayesPrint ("   moves that are currently switched on for the model, use 'showmoves'. You can  \n");
        MrBayesPrint ("   also see other available moves by using 'showmoves allavailable=yes'. A list  \n");
        MrBayesPrint ("   of the moves for each parameter in the model is available by using the command\n");
        MrBayesPrint ("   'Showparams'. If you change proposal probabilities, make sure that all        \n");
        MrBayesPrint ("   parameters that are not fixed in your model have at least one move switched   \n");
        MrBayesPrint ("   on.                                                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   One word of warning: You should be extremely careful when modifying any       \n");
        MrBayesPrint ("   of the chain parameters using 'propset'. It is quite possible to completely   \n");
        MrBayesPrint ("   wreck any hope of achieving convergence by inappropriately setting the        \n");
        MrBayesPrint ("   tuning parameters. In general, you want to set move tuning parameters such    \n");
        MrBayesPrint ("   that the acceptance rate of the move is intermediate (we suggest targeting    \n");
        MrBayesPrint ("   the range 10%% to 70%% acceptance, if possible). If the acceptance rate is    \n");
        MrBayesPrint ("   outside of this range, the MCMC chain will probably not sample that parameter \n");
        MrBayesPrint ("   very efficiently. The acceptance rates for all moves in the cold chain(s) are \n");
        MrBayesPrint ("   summarized at the end of each run in the screen output. The acceptance rates  \n");
        MrBayesPrint ("   (potentially for all chains, cold and heated) are also printed to the .mcmc   \n");
        MrBayesPrint ("   file if Mcmc convergence diagnostics are turned on (using 'Mcmc' or 'Mcmcp'). \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Log"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Log                                                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command allows output to the screen to also be output to a file.         \n");
        MrBayesPrint ("   The useage is:                                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      log start/stop filename=<name> append/replace                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The options are:                                                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Start/Stop     -- Starts or stops logging of output to file.                  \n");
        MrBayesPrint ("   Append/Replace -- Either append to or replace existing file.                  \n");
        MrBayesPrint ("   Filename       -- Name of log file (currently, the name of the log            \n");
        MrBayesPrint ("                     file is \"%s\").\n", logFileName);
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Translate"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Translate                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used by MrBayes to specify the mapping between taxon names    \n");
        MrBayesPrint ("   and taxon numbers in a Nexus tree file. For instance,                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      translate                                                                  \n");
        MrBayesPrint ("         1 Homo,                                                                 \n");
        MrBayesPrint ("         2 Pan,                                                                  \n");
        MrBayesPrint ("         3 Gorilla,                                                              \n");
        MrBayesPrint ("         4 Hylobates;                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   establishes that the taxon labeled 1 in the trees that follow is Homo, the    \n");
        MrBayesPrint ("   taxon labeled 2 is Pan, etc.                                                  \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Usertree"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Usertree                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command allows you to specify a user tree. The user tree can then be     \n");
        MrBayesPrint ("   used as a starting tree for a MCMC analysis. The format for the command is    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      usertree = <tree in Newick format>                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example,                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      usertree = (A,B,(C,D))                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   specifies an unrooted tree of four species. Note that the program re-         \n");
        MrBayesPrint ("   quires that trees are binary (i.e., strictly bifurcating). Hence, there       \n");
        MrBayesPrint ("   can be only one three-way split, as shown in the example. If the tree         \n");
        MrBayesPrint ("   is not binary, the program will return an error.                              \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Mcmc"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Mcmc                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command starts the Markov chain Monte Carlo (MCMC) analysis. The         \n");
        MrBayesPrint ("   posterior probability of phylogenetic trees (and other parameters of the      \n");
        MrBayesPrint ("   substitution model) cannot be determined analytically. Instead, MCMC is       \n");
        MrBayesPrint ("   used to approximate the posterior probabilities of trees by drawing           \n");
        MrBayesPrint ("   (dependent) samples from the posterior distribution. This program can         \n");
        MrBayesPrint ("   implement a variant of MCMC called \"Metropolis-coupled Markov chain Monte    \n");
        MrBayesPrint ("   Carlo\", or MCMCMC for short. Basically, \"Nchains\" are run, with            \n");
        MrBayesPrint ("   Nchains - 1 of them heated. The chains are labelled 1, 2, ..., Nchains.       \n");
        MrBayesPrint ("   The heat that is applied to the i-th chain is B = 1 / (1 + temp X i). B       \n");
        MrBayesPrint ("   is the power to which the posterior probability is raised. When B = 0, all    \n");
        MrBayesPrint ("   trees have equal probability and the chain freely visits trees. B = 1 is      \n");
        MrBayesPrint ("   the \"cold\" chain (or the distribution of interest). MCMCMC can mix          \n");
        MrBayesPrint ("   better than ordinary MCMC; after all of the chains have gone through          \n");
        MrBayesPrint ("   one cycle, two chains are chosen at random and an attempt is made to          \n");
        MrBayesPrint ("   swap the states (with the probability of a swap being determined by the       \n");
        MrBayesPrint ("   Metropolis et al. equation). This allows the chain to potentially jump        \n");
        MrBayesPrint ("   a valley in a single bound. The correct usage is                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      mcmc <parameter> = <value> ... <parameter> = <value>                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example,                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      mcmc ngen=100000 nchains=4 temp=0.5                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   performs a MCMCMC analysis with four chains with the temperature set to       \n");
        MrBayesPrint ("   0.5. The chains would be run for 100,000 cycles.                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Ngen         -- This option sets the number of cycles for the MCMC alg-       \n");
        MrBayesPrint ("                   orithm. This should be a big number as you want the chain     \n");
        MrBayesPrint ("                   to first reach stationarity, and then remain there for        \n");
        MrBayesPrint ("                   enough time to take lots of samples.                          \n");
        MrBayesPrint ("   Nruns        -- How many independent analyses are started simultaneously.     \n");
        MrBayesPrint ("   Nchains      -- How many chains are run for each analysis for the MCMCMC      \n");
        MrBayesPrint ("                   variant. The default is 4: 1 cold chain and 3 heated chains.  \n");
        MrBayesPrint ("                   If Nchains is set to 1, MrBayes will use regular MCMC sam-    \n");
        MrBayesPrint ("                   pling, without heating.                                       \n");
        MrBayesPrint ("   Temp         -- The temperature parameter for heating the chains. The higher  \n");
        MrBayesPrint ("                   the temperature, the more likely the heated chains are to     \n");
        MrBayesPrint ("                   move between isolated peaks in the posterior distribution.    \n");
        MrBayesPrint ("                   However, excessive heating may lead to very low acceptance    \n");
        MrBayesPrint ("                   rates for swaps between different chains. Before changing the \n");
        MrBayesPrint ("                   default setting, however, note that the acceptance rates of   \n");
        MrBayesPrint ("                   swaps tend to fluctuate during the burn-in phase of the run.  \n");
        MrBayesPrint ("   Reweight     -- Here, you specify three numbers, that respectively represent  \n");
        MrBayesPrint ("                   the percentage of characters to decrease in weight, the       \n");
        MrBayesPrint ("                   percentage of characters to increase in weight, and the       \n");
        MrBayesPrint ("                   increment. An increase/decrease in weight is acheived by      \n");
        MrBayesPrint ("                   replicating/removing a character in the matrix. This is       \n");
        MrBayesPrint ("                   only done to non-cold chains. The format for this parameter   \n");
        MrBayesPrint ("                   is \"reweight=(<number>,<number>)\" or \"reweight=(<number>,  \n");
        MrBayesPrint ("                   <number>,<number>)\".                                         \n");
        MrBayesPrint ("   Swapfreq     -- This specifies how often swaps of states between chains are   \n");
        MrBayesPrint ("                   attempted. You must be running at least two chains for this   \n");
        MrBayesPrint ("                   option to be relevant. The default is Swapfreq=1, resulting   \n");
        MrBayesPrint ("                   in Nswaps (see below) swaps being tried each generation of    \n");
        MrBayesPrint ("                   the run. If Swapfreq is set to 10, then Nswaps swaps will be  \n");
        MrBayesPrint ("                   tried every tenth generation of the run.                      \n");
        MrBayesPrint ("   Nswaps       -- The number of swaps tried for each swapping generation of the \n");
        MrBayesPrint ("                   chain (see also Swapfreq).                                    \n");
        MrBayesPrint ("   Samplefreq   -- This specifies how often the Markov chain is sampled. You     \n");
        MrBayesPrint ("                   can sample the chain every cycle, but this results in very    \n");
        MrBayesPrint ("                   large output files. Thinning the chain is a way of making     \n");
        MrBayesPrint ("                   these files smaller and making the samples more independent.  \n");
        MrBayesPrint ("   Printfreq    -- This specifies how often information about the chain is       \n");
        MrBayesPrint ("                   printed to the screen.                                        \n");
        MrBayesPrint ("   Printall     -- If set to NO, only cold chains in a MCMC analysis are printed \n");
        MrBayesPrint ("                   to screen. If set to YES, both cold and heated chains will be \n");
        MrBayesPrint ("                   output. This setting only affects the printing to screen, it  \n");
        MrBayesPrint ("                   does not change the way values are written to file.           \n");
        MrBayesPrint ("   Printmax     -- The maximum number of chains to print to screen.              \n");
        MrBayesPrint ("   Mcmcdiagn    -- Determines whether acceptance ratios of moves and swaps will  \n");
        MrBayesPrint ("                   be printed to file. The file will be named similarly to the   \n");
        MrBayesPrint ("                   '.p' and '.t' files, but will have the ending '.mcmc'. If     \n");
        MrBayesPrint ("                   more than one independent analysis is run simultaneously (see \n");
        MrBayesPrint ("                   Nruns below), convergence diagnostics for tree topology will  \n");
        MrBayesPrint ("                   also be printed to this file. The convergence diagnostic used \n");
        MrBayesPrint ("                   is the average standard deviation in partition frequency      \n");
        MrBayesPrint ("                   values across independent analyses. The Burnin setting (see   \n");
        MrBayesPrint ("                   below) determines how many samples will be discarded as burnin\n");
        MrBayesPrint ("                   before calculating the partition frequencies. The Minpartfreq \n");
        MrBayesPrint ("                   setting (see below) determines the minimum partition frequency\n");
        MrBayesPrint ("                   required for a partition to be included in the calculation. As\n");
        MrBayesPrint ("                   the independent analyses approach stationarity (converge), the\n");
        MrBayesPrint ("                   value of the diagnostic is expected to approach zero.         \n");
        MrBayesPrint ("   Diagnfreq    -- The number of generations between the calculation of MCMC     \n");
        MrBayesPrint ("                   diagnostics (see Mcmcdiagn above).                            \n");
        MrBayesPrint ("   Diagnstat    -- The statistic to use for run-time convergence diagnostics.    \n");
        MrBayesPrint ("                   Choices are 'Avgstddev' for average standard deviation of     \n");
        MrBayesPrint ("                   split frequencies and 'Maxstddev' for maximum standard devia- \n");
        MrBayesPrint ("                   tion of split frequencies.                                    \n");
        MrBayesPrint ("   Savetrees    -- If you are using a relative burnin for run-time convergence   \n");
        MrBayesPrint ("                   diagnostics, tree samples need to be deleted from split       \n");
        MrBayesPrint ("                   frequency counters as the cut-off point for the burnin moves  \n");
        MrBayesPrint ("                   during the run. If 'Savetrees' is set to 'No', tree samples   \n");
        MrBayesPrint ("                   to be discarded are read back in from file. If 'Savetrees' is \n");
        MrBayesPrint ("                   set to 'Yes', the tree samples to be removed will be stored   \n");
        MrBayesPrint ("                   in the internal memory instead. This can use up a lot of      \n");
        MrBayesPrint ("                   memory in large analyses.                                     \n");
        MrBayesPrint ("   Minpartfreq  -- The minimum frequency required for a partition to be included \n");
        MrBayesPrint ("                   in the calculation of the topology convergence diagnostic. The\n");
        MrBayesPrint ("                   partition is included if the minimum frequency is reached in  \n");
        MrBayesPrint ("                   at least one of the independent tree samples that are com-    \n");
        MrBayesPrint ("                   pared.                                                        \n");
        MrBayesPrint ("   Allchains    -- If this option is set to YES, acceptance ratios for moves are \n");
        MrBayesPrint ("                   recorded for all chains, cold or heated. By default, only the \n");
        MrBayesPrint ("                   acceptance ratios for the cold chain are recorded.            \n");
        MrBayesPrint ("   Allcomps     -- If this option is set to YES, topological convergence diag-   \n");
        MrBayesPrint ("                   nostics are calculated over all pairwise comparisons of runs. \n");
        MrBayesPrint ("                   If it is set to NO, only the overall value is reported.       \n");
        MrBayesPrint ("   Relburnin    -- If this option is set to YES, then a proportion of the sampled\n");
        MrBayesPrint ("                   values will be discarded as burnin when calculating the con-  \n");
        MrBayesPrint ("                   vergence diagnostic. The proportion to be discarded is set    \n");
        MrBayesPrint ("                   with Burninfrac (see below). When the Relburnin option is set \n");
        MrBayesPrint ("                   to NO, then a specific number of samples will be discarded    \n");
        MrBayesPrint ("                   instead. This number is set by Burnin (see below).            \n");
        MrBayesPrint ("   Burnin       -- Determines the number of samples (not generations) that will  \n");
        MrBayesPrint ("                   be discarded when convergence diagnostics are calculated.     \n");
        MrBayesPrint ("                   The value of this option is only relevant when Relburnin is   \n");
        MrBayesPrint ("                   set to NO.                                                    \n");
        MrBayesPrint ("   BurninFrac   -- Determines the fraction of samples that will be discarded     \n");
        MrBayesPrint ("                   when convergence diagnostics are calculated. The value of     \n");
        MrBayesPrint ("                   this option is only relevant when Relburnin is set to YES.    \n");
        MrBayesPrint ("                   Example: A value for this option of 0.25 means that 25%% of   \n");
        MrBayesPrint ("                   the samples will be discarded.                                \n");
        MrBayesPrint ("   Stoprule     -- If this option is set to NO, then the chain is run the number \n");
        MrBayesPrint ("                   of generations determined by Ngen. If it is set to YES, and   \n");
        MrBayesPrint ("                   topological convergence diagnostics are calculated (Mcmcdiagn \n");
        MrBayesPrint ("                   is set to YES), then the chain will be stopped before the pre-\n");
        MrBayesPrint ("                   determined number of generations if the convergence diagnostic\n");
        MrBayesPrint ("                   falls below the stop value.                                   \n");
        MrBayesPrint ("   Stopval      -- The critical value for the topological convergence diagnostic.\n");
        MrBayesPrint ("                   Only used when Stoprule and Mcmcdiagn are set to yes, and     \n");
        MrBayesPrint ("                   more than one analysis is run simultaneously (Nruns > 1).     \n");
        MrBayesPrint ("   Checkpoint   -- If this parameter is set to 'Yes', all the current parameter  \n");
        MrBayesPrint ("                   values of all chains will be printed to a check-pointing file \n");
        MrBayesPrint ("                   every 'Checkfreq' generation of the analysis. The file will be\n");
        MrBayesPrint ("                   named <Filename>.ckp and allows you to restart the analysis   \n");
        MrBayesPrint ("                   from the last check point. This can be handy if you are       \n");
        MrBayesPrint ("                   running a long analysis and want to extend it, or if there is \n");
        MrBayesPrint ("                   a risk that a long analysis will be inadvertently interupted  \n");
        MrBayesPrint ("                   by hardware failure or other factors that are out of your     \n");
        MrBayesPrint ("                   control.                                                      \n");
        MrBayesPrint ("   Checkfreq    -- The number of generations between check-pointing. See the     \n");
        MrBayesPrint ("                   'Checkpoint' parameter above for more information.            \n");
        MrBayesPrint ("   Filename     -- The name of the files that will be generated. Two files       \n");
        MrBayesPrint ("                   are generated: \"<Filename>.t\" and \"<Filename>.p\".         \n");
        MrBayesPrint ("                   The .t file contains the trees whereas the .p file con-       \n");
        MrBayesPrint ("                   tains the sampled values of the parameters.                   \n");
        MrBayesPrint ("   Startparams  -- The starting values for the model parameters are set to       \n");
        MrBayesPrint ("                   arbitrary or random values when the parameters are created.   \n");
        MrBayesPrint ("                   These starting values can be altered using the 'Startvals'    \n");
        MrBayesPrint ("                   command. The 'Startparams=reset' option allows you to reset   \n");
        MrBayesPrint ("                   the starting values to the default at the start of the ana-   \n");
        MrBayesPrint ("                   lysis, overriding any previous user-defined starting values.  \n");
        MrBayesPrint ("                   Under the default option, 'current', the chains will use the  \n");
        MrBayesPrint ("                   current starting values.                                      \n");
        MrBayesPrint ("   Starttree    -- The starting tree(s) for the chain can either be randomly     \n");
        MrBayesPrint ("                   selected or user-defined. It might be a good idea to          \n");
        MrBayesPrint ("                   start from randomly chosen trees; convergence seems           \n");
        MrBayesPrint ("                   likely if independently run chains, each of which             \n");
        MrBayesPrint ("                   started from different random trees, converge to the same     \n");
        MrBayesPrint ("                   answer. If you want the chain to start from user-defined      \n");
        MrBayesPrint ("                   trees instead, you first need to read in your tree(s) from a  \n");
        MrBayesPrint ("                   Nexus file with a 'trees' block, and then you need to set the \n");
        MrBayesPrint ("                   starting tree(s) using the 'Startvals' command. Finally, you  \n");
        MrBayesPrint ("                   need to make sure that 'Starttree' is set to 'current'. If    \n");
        MrBayesPrint ("                   you do not set the starting tree(s), the chains will start    \n");
        MrBayesPrint ("                   with random trees. Setting 'Starttree' to 'random' causes     \n");
        MrBayesPrint ("                   new starting trees to be drawn randomly at the start of the   \n");
        MrBayesPrint ("                   run, overwriting any previous user-defined starting trees.    \n");
        MrBayesPrint ("   Nperts       -- This is the number of random perturbations to apply to the    \n");
        MrBayesPrint ("                   user starting tree. This allows you to have something         \n");
        MrBayesPrint ("                   between completely random and user-defined trees start        \n");
        MrBayesPrint ("                   the chain.                                                    \n");
        MrBayesPrint ("   Data         -- When Data is set to NO, the chain is run without data. This   \n");
        MrBayesPrint ("                   should be used only for examining induced priors. DO NOT SET  \n");
        MrBayesPrint ("                   'DATA' TO 'NO' UNLESS YOU KNOW WHAT YOU ARE DOING!            \n");
        MrBayesPrint ("   Ordertaxa    -- Determines whether taxa should be ordered before trees are    \n");
        MrBayesPrint ("                   printed to file. If set to 'Yes', terminals in the sampled    \n");
        MrBayesPrint ("                   trees will be reordered to match the order of the taxa in the \n");
        MrBayesPrint ("                   data matrix as closely as possible. By default, trees will be \n");
        MrBayesPrint ("                   printed without reordering of taxa.                           \n");
        MrBayesPrint ("   Append       -- Set this to 'Yes' to append the results of the current run to \n");
        MrBayesPrint ("                   a previous run. MrBayes will first read in the results of the \n");
        MrBayesPrint ("                   previous run (number of generations and sampled splits) and   \n");
        MrBayesPrint ("                   will then continue that run where you left it off. Make sure  \n");
        MrBayesPrint ("                   that the output file names used in the previous run are the   \n");
        MrBayesPrint ("                   same as those in the current run.                             \n");
        MrBayesPrint ("   Autotune     -- Set this to 'Yes' to autotune the proposals that change       \n");
        MrBayesPrint ("                   substitution model parameters. When set to 'No', the tuning   \n");
        MrBayesPrint ("                   parameters are fixed to their starting values. Note that the  \n");
        MrBayesPrint ("                   autotuning occurs independently for each chain. The target    \n");
        MrBayesPrint ("                   acceptance rate for each move can be changed using the        \n");
        MrBayesPrint ("                   'Propset' command.                                            \n");
        MrBayesPrint ("   Tunefreq     -- When a proposal has been tried 'Tunefreq' times, its tuning   \n");
        MrBayesPrint ("                   parameter is adjusted to reach the target acceptance rate     \n");
        MrBayesPrint ("                   if 'Autotune' is set to 'Yes'.                                \n");
        MrBayesPrint ("                                                                                 \n");
        PrintSettings ("Mcmc");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Mcmcp"))
        {
        // PrintYesNo (chainParams.saveBrlens, yesNoStr);
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Mcmcp                                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command sets the parameters of the Markov chain Monte Carlo (MCMC)       \n");
        MrBayesPrint ("   analysis without actually starting the chain. This command is identical       \n");
        MrBayesPrint ("   in all respects to Mcmc, except that the analysis will not start after        \n");
        MrBayesPrint ("   this command is issued. For more details on the options, check the help       \n");
        MrBayesPrint ("   menu for Mcmc.\n");
        MrBayesPrint ("                                                                                 \n");
        PrintSettings ("Mcmc");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Ss"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Ss                                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used to start stepping-stone sampling, which is an efficient  \n");
        MrBayesPrint ("   and accurate method for estimating the marginal likelihood of the currently   \n");
        MrBayesPrint ("   specified model. It is considerably more accurate than the harmonic mean of   \n");
        MrBayesPrint ("   the likelihoods from a standard MCMC run on the model (calculated by the      \n");
        MrBayesPrint ("   'Sump' command) but it requires a separate MCMC-like run. To be more specific,\n");
        MrBayesPrint ("   stepping-stone sampling uses importance sampling to estimate each ratio in a  \n");
        MrBayesPrint ("   series of discrete steps bridging the posterior and prior distributions.      \n");
        MrBayesPrint ("   The importance distributions that are used are called power posterior distri- \n");
        MrBayesPrint ("   butions, and are defined as prior*(likelihood^beta). By varying beta from 1 to\n");
        MrBayesPrint ("   0, we get a series of distributions that connect the posterior (beta = 1) to  \n");
        MrBayesPrint ("   the prior (beta = 0).                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The power posterior distributions are sampled using MCMC. First, we start a   \n");
        MrBayesPrint ("   standard MCMC chain on the posterior distribution, and let it run until we    \n");
        MrBayesPrint ("   have reached the criterion specified by the 'Burninss' option. After this, we \n");
        MrBayesPrint ("   step through the power posterior distributions until we reach the prior dis-  \n");
        MrBayesPrint ("   tribution. In each of the 'Nsteps' steps, we sample from a new power poster-  \n");
        MrBayesPrint ("   ior distribution with a distinct beta value. The beta values correspond to    \n");
        MrBayesPrint ("   'Nsteps' evenly spaced quantiles in a Beta distribution with the parameters   \n");
        MrBayesPrint ("   'Alpha' and 1.0. For the first sampling step, the beta value is equal to the  \n");
        MrBayesPrint ("   last quantile, i.e., it is close to 1.0. For each successive step, the beta   \n");
        MrBayesPrint ("   value takes on the value of the next quantile, in decreasing order, until it  \n");
        MrBayesPrint ("   reaches the value of 0.0. If you change value of 'FromPrior' from default 'No'\n");
        MrBayesPrint ("   to 'Yes' then the direction of power posterior change during SS analizes is   \n");
        MrBayesPrint ("   opposite to the one described above, i.e. we start from sampling prior and    \n");
        MrBayesPrint ("   finish close to posterior.                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The 'Ss' procedure uses the same machinery as the standard 'Mcmc' algorithm,  \n");
        MrBayesPrint ("   and shares most of its parameters with the 'Mcmc' and 'Mcmcp' commands. All   \n");
        MrBayesPrint ("   'Mcmc' parameters, except those related to burnin, have the same meaning and  \n");
        MrBayesPrint ("   usage in the 'Ss' command as they have in the 'Mcmc' command. The 'Mcmc'      \n");
        MrBayesPrint ("   burnin parameters are used to set up burnin within each step. The 'Ss' command\n");
        MrBayesPrint ("   also uses its own burnin parameter, 'Burninss' (see below for details). The   \n");
        MrBayesPrint ("   'Ss' command also has its own parameters for specifying the number of steps   \n");
        MrBayesPrint ("   and the shape of the Beta distribution from which the beta values are computed\n");
        MrBayesPrint ("   (see below).                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that the 'Ngen' parameter of 'Mcmc' is used to set the maximum number of \n");
        MrBayesPrint ("   generations processed, including both the burnin and the following steps in   \n");
        MrBayesPrint ("   the stepping-stone sampling phase. For instance, assume that 'Burninss' is set\n");
        MrBayesPrint ("   to '-1', 'Nsteps' to '49', 'Ngen' to '1000000' and 'Samplefreq' to '1000'.    \n");
        MrBayesPrint ("   We will then get 1,000 samples in total (1,000,000 / 1,000). These will fall  \n");
        MrBayesPrint ("   into 50 bins, one of which represents the burnin and is discarded. Each step  \n");
        MrBayesPrint ("   in the algorithm will thus be represented by 20 samples.                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   More information on 'Mcmc' parameters is available in the help for the 'Mcmc' \n");
        MrBayesPrint ("   and 'Mcmcp' commands. Only the exclusive 'Ss' parameters are listed below.    \n");
        MrBayesPrint ("   These can only be set up using the 'Ss' command, while the parameters shared  \n");
        MrBayesPrint ("   with 'Mcmc' and 'Mcmcp' can also be set up using those commands.              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The correct usage is                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      ss <parameter>=<value> ... <parameter>=<value>                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that a command:                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      ss <setting parameters shared with mcmc> <setting exclusive ss parameters> \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   would be equivalent to executing two commands:                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("     mcmcp <setting parameters shared with mcmc>;                                \n");
        MrBayesPrint ("     ss <setting exclusive ss parameters>;                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For more information on the stepping-stone algorithm, see:                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Xie, W., P. O. Lewis, Y. Fan, L. Kuo, and M.-H. Chen. 2011. Improving marginal\n");
        MrBayesPrint ("      likelihood estimation for Bayesian phylogenetic model selection. Systematic\n");
        MrBayesPrint ("      Biology 60:150-160.                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Available options:                                                            \n");
        MrBayesPrint ("   (NB: Only exclusive ss parameters listed here. For additional parameters, see \n");
        MrBayesPrint ("        help on 'mcmc' or 'mcmcp'.                                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Alpha        -- The beta values used in the stepping-stone sampling procedure \n");
        MrBayesPrint ("                   correspond to evenly spaced quantiles from a Beta('Alpha',1.0)\n");
        MrBayesPrint ("                   distribution. The parameter 'Alpha' determines the skewness of\n");
        MrBayesPrint ("                   the beta values. If 'Alpha' is set to '1.0', the beta values  \n");
        MrBayesPrint ("                   would be spaced uniformly on the interval (0.0,1.0). However, \n");
        MrBayesPrint ("                   better results are obtained if the beta values are skewed.    \n");
        MrBayesPrint ("                   Empirically, it was observed that 'Alpha' values in the range \n");
        MrBayesPrint ("                   of 0.3 to 0.5 produce the most accurate results.              \n");
        MrBayesPrint ("   Burninss     -- Fixed number of samples discarded before sampling of the first\n");
        MrBayesPrint ("                   step starts. 'Burninss' can be specified using either a pos-  \n");
        MrBayesPrint ("                   itive or a negative number. If the number is positive, it is  \n");
        MrBayesPrint ("                   interpreted as the number of samples to discard as burnin. If \n");
        MrBayesPrint ("                   the number is negative, its absolute value is interpreted as  \n");
        MrBayesPrint ("                   the length of the burnin in terms of the length of each of the\n");
        MrBayesPrint ("                   following steps in the stepping-stone algorithm. For instance,\n");
        MrBayesPrint ("                   a value of '-1' means that the length of the burnin is the    \n");
        MrBayesPrint ("                   same as the length of each of the subsequent steps.           \n");
        MrBayesPrint ("   Nsteps       -- Number of steps in the stepping-stone algorithm. Typically, a \n");
        MrBayesPrint ("                   number above 30 is sufficient for accurate results.           \n");
        MrBayesPrint ("   FromPrior    -- If it is set to 'Yes', it indicates that in the first step we \n"); 
        MrBayesPrint ("                   sample from the prior, with each consequtive step we sample   \n");
        MrBayesPrint ("                   closer to the posterior. 'No' indicates the opposite direction\n");
        MrBayesPrint ("                   of power posterior change, i.e. in the first step we sample   \n");
        MrBayesPrint ("                   close to the posterior, and with each consequtive step we     \n");
        MrBayesPrint ("                   sample closer to the prior.                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Current settings:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter          Options               Current Setting                      \n");
        MrBayesPrint ("   --------------------------------------------------------                      \n");
        MrBayesPrint ("   Alpha              <number>              %1.2lf\n", chainParams.alphaSS);
        MrBayesPrint ("   BurninSS           <number>              %d\n", chainParams.burninSS);
        MrBayesPrint ("   Nsteps             <number>              %d\n", chainParams.numStepsSS);
        MrBayesPrint ("   FromPrior           Yes/No               %s                                   \n", chainParams.startFromPriorSS == YES ? "Yes" : "No");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
else if (!strcmp(helpTkn, "Ssp"))
        {
        // PrintYesNo (chainParams.saveBrlens, yesNoStr);
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Ssp                                                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command sets the parameters of the stepping-stone sampling               \n");
        MrBayesPrint ("   analysis without actually starting the chain. This command is identical       \n");
        MrBayesPrint ("   in all respects to Ss, except that the analysis will not start after          \n");
        MrBayesPrint ("   this command is issued. For more details on the options, check the help       \n");
        MrBayesPrint ("   menu for Ss.\n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Current settings:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter          Options               Current Setting                      \n");
        MrBayesPrint ("   --------------------------------------------------------                      \n");
        MrBayesPrint ("   Alpha              <number>              %1.2lf\n", chainParams.alphaSS);
        MrBayesPrint ("   BurninSS           <number>              %d\n", chainParams.burninSS);
        MrBayesPrint ("   Nsteps             <number>              %d\n", chainParams.numStepsSS);
        MrBayesPrint ("   FromPrior           Yes/No               %s                                   \n", chainParams.startFromPriorSS == YES ? "Yes" : "No");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
else if (!strcmp(helpTkn, "Set"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Set                                                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used to set some general features of the model or program     \n");
        MrBayesPrint ("   behavior. The correct usage is                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      set <parameter>=<value> ... <parameter>=<value>                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Available options:                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Seed         -- Sets the seed number for the random number generator. The     \n");
        MrBayesPrint ("                   random number seed is initialized haphazardly at the beg-     \n");
        MrBayesPrint ("                   inning of each MrBayes session. This option allows you to     \n");
        MrBayesPrint ("                   set the seed to some specific value, thereby allowing you     \n");
        MrBayesPrint ("                   to exactly repeat an analysis. If the analysis uses swapping  \n");
        MrBayesPrint ("                   between cold and heated chains, you must also set the swap    \n");
        MrBayesPrint ("                   seed (see below) to exactly repeat the analysis.              \n");
        MrBayesPrint ("   Swapseed     -- Sets the seed used for generating the swapping sequence       \n");
        MrBayesPrint ("                   when Metropolis-coupled heated chains are used. This seed     \n");
        MrBayesPrint ("                   is initialized haphazardly at the beginning of each MrBayes   \n");
        MrBayesPrint ("                   session. This option allows you to set the seed to some       \n");
        MrBayesPrint ("                   specific value, thereby allowing you to exactly repeat a      \n");
        MrBayesPrint ("                   swap sequence. See also the 'Seed' option.                    \n");
        MrBayesPrint ("   Dir          -- The working directory. Specifies the absolute or relative path\n");
        MrBayesPrint ("                   to the working directory. If left empty, the working directory\n");
        MrBayesPrint ("                   is the current directory.                                     \n");
        MrBayesPrint ("   Partition    -- Set this option to a valid partition id, either the number or \n");
        MrBayesPrint ("                   name of a defined partition, to enforce a specific partition- \n");
        MrBayesPrint ("                   ing of the data. When a data matrix is read in, a partition   \n");
        MrBayesPrint ("                   called \"Default\" is automatically created. It divides the   \n");
        MrBayesPrint ("                   data into one part for each data type. If you only have one   \n");
        MrBayesPrint ("                   data type, DNA for instance, the default partition will not   \n");
        MrBayesPrint ("                   divide up the data at all. The default partition is always    \n");
        MrBayesPrint ("                   the first partition, so 'set partition=1' is the same as      \n");
        MrBayesPrint ("                   'set partition=default'.                                      \n");
        MrBayesPrint ("   Speciespartition -- Set this option to a valid speciespartition id, either the\n");
        MrBayesPrint ("                   number or name of a defined speciespartition, to enforce a    \n");
        MrBayesPrint ("                   specific partitioning of taxa to species. When a data matrix  \n");
        MrBayesPrint ("                   is read in, a speciespartition called \"Default\" is auto-    \n");
        MrBayesPrint ("                   matically created. It assigns one taxon for each species. The \n"); 
        MrBayesPrint ("                   default speciespartition is always the first speciespartition,\n");
        MrBayesPrint ("                   so 'set speciespartition=1' is the same as                    \n");
        MrBayesPrint ("                   'set speciespartition=default'.                               \n");
        MrBayesPrint ("   Autoclose    -- If autoclose is set to 'yes', then the program will not prompt\n");
        MrBayesPrint ("                   you during the course of executing a file. This is particular-\n");
        MrBayesPrint ("                   ly useful when you run MrBayes in batch mode.                 \n");
        MrBayesPrint ("   Nowarnings   -- If nowarnings is set to yes, then the program will not prompt \n");
        MrBayesPrint ("                   you when overwriting or appending an ouput file that is al-   \n");
        MrBayesPrint ("                   ready present. If 'nowarnings=no' (the default setting), then \n");
        MrBayesPrint ("                   the program propts the user before overwriting output files.  \n");
        MrBayesPrint ("   Autoreplace  -- When nowarnings is set to yes, then MrBayes will by default   \n");
        MrBayesPrint ("                   overwrite output files that already exists. This may cause    \n");
        MrBayesPrint ("                   irrecoverable loss of previous results if you have not removed\n");
        MrBayesPrint ("                   or renamed the files from previous runs. To override this be- \n");
        MrBayesPrint ("                   havior, set autoreplace to no, in which case new output will  \n");
        MrBayesPrint ("                   be appended to existing files instead.                        \n");
        MrBayesPrint ("   Quitonerror  -- If quitonerror is set to yes, then the program will quit when \n");
        MrBayesPrint ("                   an error is encountered, after printing an error message. If  \n");
        MrBayesPrint ("                   quitonerror is set to no (the default setting), then the      \n");
        MrBayesPrint ("                   program will wait for additional commands from the command    \n");
        MrBayesPrint ("                   line after the error message is printed.                      \n");
        MrBayesPrint ("   Scientific   -- Set this option to 'Yes' to write sampled values to file in   \n");
        MrBayesPrint ("                   scientific format and to 'No' to write them in fixed format.  \n");
        MrBayesPrint ("                   Fixed format is easier for humans to read but you risk losing \n");
        MrBayesPrint ("                   precision for small numbers. For instance, sampled values that\n");
        MrBayesPrint ("                   are less than 1E-6 will print to file as '0.000000' if fixed  \n");
        MrBayesPrint ("                   format is used and 'precision' is set to 6.                   \n");
        MrBayesPrint ("   Precision    -- Precision allows you to set the number of decimals to be prin-\n");
        MrBayesPrint ("                   ted when sampled values are written to file. Precision must be\n");
        MrBayesPrint ("                   in the range 3 to 15.                                         \n");
#   if defined (BEAGLE_ENABLED)
        MrBayesPrint ("   Usebeagle    -- Set this option to 'Yes' to attempt to use the BEAGLE library \n");
        MrBayesPrint ("                   to compute the phylogenetic likelihood on a variety of high-  \n");
        MrBayesPrint ("                   performance hardware including multicore CPUs and GPUs. Some  \n"); 
        MrBayesPrint ("                   models in MrBayes are not yet supported by BEAGLE.            \n");               
        MrBayesPrint ("   Beagleresource -- Set this option to the number of a specific resource you    \n");
        MrBayesPrint ("                   wish to use with BEAGLE (use 'Showbeagle' to see the list of  \n");
        MrBayesPrint ("                   available resources). Set to '99' for auto-resource selection.\n");
        MrBayesPrint ("   Beagledevice -- Set this option to 'GPU' or 'CPU' to select processor.        \n"); 
        MrBayesPrint ("   Beagleprecision -- Selection 'Single' or 'Double' precision computation.      \n");
        MrBayesPrint ("   Beaglescaling -- 'Always' rescales partial likelihoods at each evaluation.    \n");
        MrBayesPrint ("                    'Dynamic' rescales less frequently and should run faster.    \n");
        MrBayesPrint ("   Beaglesse    -- Use SSE instructions on Intel CPU processors.                 \n");
        MrBayesPrint ("   Beagleopenmp -- Use OpenMP to parallelize across multi-core CPU processors.   \n");
#   endif
#   if defined (THREADS_ENABLED)
        MrBayesPrint ("   Beaglethreads -- Set this option to 'Yes' to employ multiple threads to drive \n");
        MrBayesPrint ("                   multiple BEAGLE resource simultaneously. This is highly       \n");
        MrBayesPrint ("                   recommended for more than one GPU, and for sufficiently large \n");
        MrBayesPrint ("                   data partitions, multi-core CPUs should also demonstrate      \n");
        MrBayesPrint ("                   speed-ups.                                                    \n");
#   endif
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Current settings:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter          Options               Current Setting                      \n");
        MrBayesPrint ("   --------------------------------------------------------                      \n");
        MrBayesPrint ("   Seed               <number>              %ld                                  \n", globalSeed);
        MrBayesPrint ("   Swapseed           <number>              %ld                                  \n", swapSeed);
        MrBayesPrint ("   Dir                <name>                \"%s\"\n", workingDir);
        if (defMatrix == YES)
            MrBayesPrint ("   Partition          <name>                %s\n", partitionNames[partitionNum]);
        else
            MrBayesPrint ("   Partition          <name>                \"\"\n");
        if (defTaxa == YES)
            MrBayesPrint ("   Speciespartition   <name>                %s\n", speciespartitionNames[speciespartitionNum]);
        else
            MrBayesPrint ("   Speciespartition   <name>                \"\"\n");
        MrBayesPrint ("   Autoclose          Yes/No                %s                                   \n", autoClose == YES ? "Yes" : "No");
        MrBayesPrint ("   Nowarnings         Yes/No                %s                                   \n", noWarn == YES ? "Yes" : "No");
        MrBayesPrint ("   Autoreplace        Yes/No                %s                                   \n", autoOverwrite == YES ? "Yes" : "No");
        MrBayesPrint ("   Quitonerror        Yes/No                %s                                   \n", quitOnError == YES ? "Yes" : "No");
        MrBayesPrint ("   Scientific         Yes/No                %s                                   \n", scientific == YES ? "Yes" : "No");
        MrBayesPrint ("   Precision          <number>              %d                                   \n", precision);
#   if defined (BEAGLE_ENABLED)
        MrBayesPrint ("   Usebeagle          Yes/No                %s                                   \n", tryToUseBEAGLE == YES ? "Yes" : "No");
        MrBayesPrint ("   Beagleresource     <number>              %d                                   \n", beagleResourceNumber);
        MrBayesPrint ("   Beagledevice       CPU/GPU               %s                                   \n", beagleFlags & BEAGLE_FLAG_PROCESSOR_GPU ? "GPU" : "CPU");
        MrBayesPrint ("   Beagleprecision    Single/Double         %s                                   \n", beagleFlags & BEAGLE_FLAG_PRECISION_SINGLE ? "Single" : "Double");
        MrBayesPrint ("   Beaglescaling      Always/Dynamic        %s                                   \n", beagleScalingScheme == MB_BEAGLE_SCALE_ALWAYS ? "Always" : "Dynamic");
        MrBayesPrint ("   Beaglesse          Yes/No                %s                                   \n", beagleFlags & BEAGLE_FLAG_VECTOR_SSE ? "Yes" : "No");
        MrBayesPrint ("   Beagleopenmp       Yes/No                %s                                   \n", beagleFlags & BEAGLE_FLAG_THREADING_OPENMP ? "Yes" : "No");        
#   endif
#   if defined (THREADS_ENABLED)
        MrBayesPrint ("   Beaglethreads      Yes/No                %s                                   \n", tryToUseThreads == YES ? "Yes" : "No");
#   endif
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Charset"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Charset                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command defines a character set. The format for the charset command      \n"); 
        MrBayesPrint ("   is                                                                            \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      charset <name> = <character numbers>                                       \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example, \"charset first_pos = 1-720\\3\" defines a character set         \n");
        MrBayesPrint ("   called \"first_pos\" that includes every third site from 1 to 720.            \n");
        MrBayesPrint ("   The character set name cannot have any spaces in it. The slash (\\)           \n");
        MrBayesPrint ("   is a nifty way of telling the program to assign every third (or               \n");
        MrBayesPrint ("   second, or fifth, or whatever) character to the character set.                \n");
        MrBayesPrint ("   This option is best used not from the command line, but rather as a           \n");
        MrBayesPrint ("   line in the mrbayes block of a file. Note that you can use \".\" to           \n");
        MrBayesPrint ("   stand in for the last character (e.g., charset 1-.\\3).                       \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Outgroup"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Outgroup                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command assigns a taxon to the outgroup. The correct usage is:           \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      outgroup <number>/<taxon name>                                             \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example, \"outgroup 3\" assigns the third taxon in the matrix to be       \n");
        MrBayesPrint ("   the outgroup. Similarly, \"outgroup Homo_sapiens\" assings the taxon          \n");
        MrBayesPrint ("   \"Homo_sapiens\" to be the outgroup (assuming that there is a taxon named     \n");
        MrBayesPrint ("   \"Homo_sapiens\" in the matrix). Only a single taxon can be assigned to       \n");
        MrBayesPrint ("   be the outgroup.                                                              \n");
        MrBayesPrint ("                                                                                 \n");
        if (defTaxa == YES)
            MrBayesPrint ("   Current outgroup: %s (taxon no. %d)\n", taxaNames[outGroupNum], outGroupNum+1);
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showusertrees"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Showusertrees                                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the currently defined user trees. The correct usage        \n");
        MrBayesPrint ("   is \"showusertrees\".                                                         \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showmcmctrees"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Showmcmctrees                                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the current trees used by the Markov chains.               \n");
        MrBayesPrint ("   is \"showmcmctrees\".                                                         \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Deroot"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Deroot                                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command deroots the user tree. If the tree is already unrooted, a        \n");
        MrBayesPrint ("   warning is issued. The correct usage is \"deroot\".                           \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Root"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Root                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command roots the tree. If the tree is already rooted, a warning         \n");
        MrBayesPrint ("   is issued. The tree is rooted at the midpoint between the outgroup species    \n");
        MrBayesPrint ("   and the ingroup species. The correct usage is \"root\".                       \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Taxset"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Taxset                                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command defines a taxon set. The format for the taxset command           \n"); 
        MrBayesPrint ("   is                                                                            \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      taxset <name> = <taxon names or numbers>                                   \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example, \"taxset apes = Homo Pan Gorilla Orang gibbon\" defines a        \n");
        MrBayesPrint ("   taxon set called \"apes\" that includes five taxa (namely, apes).             \n");
        MrBayesPrint ("   You can assign up to 30 taxon sets. This option is best used                  \n");
        MrBayesPrint ("   not from the command line but rather as a line in the mrbayes block           \n");
        MrBayesPrint ("   of a file.                                                                    \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Taxlabels"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Taxlabels                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command defines taxon labels. It could be used within taxa block.        \n"); 
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Charstat"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Charstat                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the status of all the characters. The correct usage        \n");
        MrBayesPrint ("   is                                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      charstat                                                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   After typing \"charstat\", the character number, whether it is excluded       \n");
        MrBayesPrint ("   or included, and the partition identity are shown. The output is paused       \n");
        MrBayesPrint ("   every 100 characters. This pause can be turned off by setting autoclose       \n");
        MrBayesPrint ("   to \"yes\" (set autoclose=yes).                                               \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Taxastat"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Taxastat                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the status of all the taxa. The correct usage is           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      taxastat                                                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   After typing \"taxastat\", the taxon number, name, and whether it is          \n");
        MrBayesPrint ("   excluded or included are shown.                                               \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Partition"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Partition                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command allows you to specify a character partition. The format for      \n"); 
        MrBayesPrint ("   this command is                                                               \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      partition <name> = <num parts>:<chars in first>, ...,<chars in last>       \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example, \"partition by_codon = 3:1st_pos,2nd_pos,3rd_pos\" specifies     \n"); 
        MrBayesPrint ("   a partition called \"by_codon\" which consists of three parts (first,         \n"); 
        MrBayesPrint ("   second, and third codon positions). Here, we are assuming that the sites      \n"); 
        MrBayesPrint ("   in each partition were defined using the charset command. You can specify     \n"); 
        MrBayesPrint ("   a partition without using charset as follows:                                 \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      partition by_codon = 3:1 4 6 9 12,2 5 7 10 13,3 6 8 11 14                  \n"); 
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   However, we recommend that you use the charsets to define a set of char-      \n"); 
        MrBayesPrint ("   acters and then use these predefined sets when defining the partition.        \n"); 
        MrBayesPrint ("   Also, it makes more sense to define a partition as a line in the mrbayes      \n"); 
        MrBayesPrint ("   block than to issue the command from the command line (then again, you        \n"); 
        MrBayesPrint ("   may be a masochist, and want to do extra work).                               \n"); 
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Exclude"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Exclude                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command excludes characters from the analysis. The correct usage is      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      exclude <number> <number> <number>                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or                                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      exclude <number> - <number>                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or                                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      exclude <charset>                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or some combination thereof. Moreover, you can use the specifier \"\\\" to    \n");
        MrBayesPrint ("   exclude every nth character. For example, the following                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      exclude 1-100\\3                                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   would exclude every third character. As a specific example,                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      exclude 2 3 10-14 22                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   excludes sites 2, 3, 10, 11, 12, 13, 14, and 22 from the analysis. Also,      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      exclude all                                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   excludes all of the characters from the analysis. Excluding all characters    \n");
        MrBayesPrint ("   does not leave you much information for inferring phylogeny.                  \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Include"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Include                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command includes characters that were previously excluded from the       \n");
        MrBayesPrint ("   analysis. The correct usage is                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      include <number> <number> <number>                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or                                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      include <number> - <number>                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or                                                                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      include <charset>                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   or some combination thereof. Moreover, you can use the specifier \"\\\" to    \n");
        MrBayesPrint ("   include every nth character. For example, the following                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      include 1-100\\3                                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   would include every third character. As a specific example,                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      include 2 3 10-14 22                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   includes sites 2, 3, 10, 11, 12, 13, 14, and 22 from the analysis. Also,      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      include all                                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   includes all of the characters in the analysis. Including all of the          \n");
        MrBayesPrint ("   characters (even if many of them are bad) is a very total-evidence-like       \n");
        MrBayesPrint ("   thing to do. Doing this will make a certain group of people very happy.       \n");
        MrBayesPrint ("   On the other hand, simply using this program would make those same people     \n");
        MrBayesPrint ("   unhappy.                                                                      \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Delete"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Delete                                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command deletes taxa from the analysis. The correct usage is:            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      delete <name and/or number and/or taxset> ...                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   A list of the taxon names or taxon numbers (labelled 1 to ntax in the order   \n");
        MrBayesPrint ("   in the matrix) or taxset(s) can be used.  For example, the following:         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      delete 1 2 Homo_sapiens                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   deletes taxa 1, 2, and the taxon labelled Homo_sapiens from the analysis.     \n");
        MrBayesPrint ("   You can also use \"all\" to delete all of the taxa. For example,              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      delete all                                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   deletes all of the taxa from the analysis. Of course, a phylogenetic anal-    \n");
        MrBayesPrint ("   ysis that does not include any taxa is fairly uninteresting.                  \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Restore"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Restore                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command restores taxa to the analysis. The correct usage is:             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      restore <name and/or number and/or taxset> ...                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   A list of the taxon names or taxon numbers (labelled 1 to ntax in the order   \n");
        MrBayesPrint ("   in the matrix) or taxset(s) can be used.  For example, the following:         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      restore 1 2 Homo_sapiens                                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   restores taxa 1, 2, and the taxon labelled Homo_sapiens to the analysis.      \n");
        MrBayesPrint ("   You can also use \"all\" to restore all of the taxa. For example,             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      restore all                                                                \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   restores all of the taxa to the analysis.                                     \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Quit"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Quit                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command quits the program. The correct usage is:                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      quit                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   It is a very easy command to use properly.                                    \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Disclaimer"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Disclaimer                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the disclaimer for the program. In short, the disclaimer   \n");
        MrBayesPrint ("   states that the authors are not responsible for any silly things you may do   \n");
        MrBayesPrint ("   to your computer or any unforseen but possibly nasty things the computer      \n");
        MrBayesPrint ("   program may inadvertently do to you.                                          \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Unlink"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Unlink                                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command unlinks model parameters across partitions of the data. The      \n");
        MrBayesPrint ("   correct usage is:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      unlink <parameter name> = (<all> or <partition list>)                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   A little background is necessary to understand this command. Upon exe-        \n");
        MrBayesPrint ("   cution of a file, a default partition is set up. This partition refer-        \n");
        MrBayesPrint ("   enced either by its name (\"default\") or number (0). If your data are        \n");
        MrBayesPrint ("   all of one type, then this default partition does not actually divide up      \n");
        MrBayesPrint ("   your characters. However, if your datatype is mixed, then the default         \n");
        MrBayesPrint ("   partition contains as many divisions as there are datatypes in your           \n");
        MrBayesPrint ("   character matrix. Of course, you can also define other partitions, and        \n");
        MrBayesPrint ("   switch among them using the set command (\"set partition=<name/number>\").    \n");
        MrBayesPrint ("   Importantly, you can also assign model parameters to individual part-         \n");
        MrBayesPrint ("   itions or to groups of them using the \"applyto\" option in lset and          \n");
        MrBayesPrint ("   prset. When the program attempts to perform an analysis, the model is         \n");
        MrBayesPrint ("   set for individual partitions. If the same parameter applies to differ-       \n");
        MrBayesPrint ("   partitions and if that parameter has the same prior, then the program         \n");
        MrBayesPrint ("   will link the parameters: that is, it will use a single value for the         \n");
        MrBayesPrint ("   parameter. The program's default, then, is to strive for parsimony.           \n");
        MrBayesPrint ("   However, there are lots of cases where you may want unlink a parameter        \n");
        MrBayesPrint ("   across partitions. For example, you may want a different transition/          \n");
        MrBayesPrint ("   transversion rate ratio to apply to different partitions. This command        \n");
        MrBayesPrint ("   allows you to unlink the parameters, or to make them different across         \n");
        MrBayesPrint ("   partitions. The converse of this command is \"link\", which links to-         \n");
        MrBayesPrint ("   gether parameters that were previously told to be different. The list         \n");
        MrBayesPrint ("   of parameters that can be unlinked includes:                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      Tratio          -- Transition/transversion rate ratio                      \n");
        MrBayesPrint ("      Revmat          -- Substitution rates of GTR model                         \n");
        MrBayesPrint ("      Omega           -- Nonsynonymous/synonymous rate ratio                     \n");
        MrBayesPrint ("      Statefreq       -- Character state frequencies                             \n");
        MrBayesPrint ("      Shape           -- Gamma/LNorm shape parameter                             \n");
        MrBayesPrint ("      Pinvar          -- Proportion of invariable sites                          \n");
        MrBayesPrint ("      Correlation     -- Correlation parameter of autodiscrete gamma             \n");
        MrBayesPrint ("      Ratemultiplier  -- Rate multiplier for partitions                          \n");
        MrBayesPrint ("      Switchrates     -- Switching rates for covarion model                      \n");
        MrBayesPrint ("      Topology        -- Topology of tree                                        \n");
        MrBayesPrint ("      Brlens          -- Branch lengths of tree                                  \n");
        MrBayesPrint ("      Speciationrate  -- Speciation rates for birth-death process                \n");
        MrBayesPrint ("      Extinctionrate  -- Extinction rates for birth-death process                \n");
    //  MrBayesPrint ("   Fossilizationrate  -- Fossilization rates for fossilized birth-death process  \n");
        MrBayesPrint ("      Popsize         -- Population size for coalescence process                 \n");
        MrBayesPrint ("      Growthrate      -- Growth rate of coalescence process                      \n"); 
        MrBayesPrint ("      Aamodel         -- Aminoacid rate matrix                                   \n"); 
        MrBayesPrint ("      Cpprate         -- Rate of Compound Poisson Process (CPP)                  \n"); 
        MrBayesPrint ("      Cppmultdev      -- Standard dev. of CPP rate multipliers (log scale)       \n"); 
        MrBayesPrint ("      Cppevents       -- CPP events                                              \n"); 
        MrBayesPrint ("      TK02var         -- Variance increase in TK02 relaxed clock model           \n"); 
        MrBayesPrint ("      Igrvar          -- Variance increase in IGR relaxed clock model            \n");
        MrBayesPrint ("      Mixedvar        -- Variance increase in Mixed relaxed clock model          \n");
    //  MrBayesPrint ("      TK02branchrates -- Branch rates of TK02  relaxed clock model               \n");
    //  MrBayesPrint ("      Igrbranchrates  -- Branch rates of IGR   relaxed clock model               \n");
    //  MrBayesPrint ("      Mixedbrchrates  -- Branch rates of Mixed relaxed clock model               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example,                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      unlink shape=(all)                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   unlinks the gamma/lnorm shape parameter across all partitions of the data.    \n");
        MrBayesPrint ("   You can use \"showmodel\" to see the current linking status of the            \n");
        MrBayesPrint ("   characters.                                                                   \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Link"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Link                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command links model parameters across partitions of the data. The        \n");
        MrBayesPrint ("   correct usage is:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      link <parameter name> = (<all> or <partition list>)                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The list of parameters that can be linked includes:                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      Tratio          -- Transition/transversion rate ratio                      \n");
        MrBayesPrint ("      Revmat          -- Substitution rates of GTR model                         \n");
        MrBayesPrint ("      Omega           -- Nonsynonymous/synonymous rate ratio                     \n");
        MrBayesPrint ("      Statefreq       -- Character state frequencies                             \n");
        MrBayesPrint ("      Shape           -- Gamma/LNorm shape parameter                             \n");
        MrBayesPrint ("      Pinvar          -- Proportion of invariable sites                          \n");
        MrBayesPrint ("      Correlation     -- Correlation parameter of autodiscrete gamma             \n");
        MrBayesPrint ("      Ratemultiplier  -- Rate multiplier for partitions                          \n");
        MrBayesPrint ("      Switchrates     -- Switching rates for covarion model                      \n");
        MrBayesPrint ("      Topology        -- Topology of tree                                        \n");
        MrBayesPrint ("      Brlens          -- Branch lengths of tree                                  \n");
        MrBayesPrint ("      Speciationrate  -- Speciation rates for birth-death process                \n");
        MrBayesPrint ("      Extinctionrate  -- Extinction rates for birth-death process                \n");
    //  MrBayesPrint ("   Fossilizationrate  -- Fossilization rates for fossilized birth-death process  \n");
        MrBayesPrint ("      Popsize         -- Population size for coalescence process                 \n");
        MrBayesPrint ("      Growthrate      -- Growth rate of coalescence process                      \n");
        MrBayesPrint ("      Aamodel         -- Aminoacid rate matrix                                   \n");
        MrBayesPrint ("      Cpprate         -- Rate of Compound Poisson Process (CPP)                  \n");
        MrBayesPrint ("      Cppmultdev      -- Standard dev. of CPP rate multipliers (log scale)       \n");
        MrBayesPrint ("      Cppevents       -- CPP events                                              \n");
        MrBayesPrint ("      TK02var         -- Variance increase in TK02 relaxed clock model           \n");
        MrBayesPrint ("      Igrvar          -- Variance increase in IGR relaxed clock model            \n");
        MrBayesPrint ("      Mixedvar        -- Variance increase in Mixed relaxed clock model          \n");
    //  MrBayesPrint ("      TK02branchrates -- Branch rates of TK02  relaxed clock model               \n");
    //  MrBayesPrint ("      Igrbranchrates  -- Branch rates of IGR   relaxed clock model               \n");
    //  MrBayesPrint ("      Mixedbrchrates  -- Branch rates of Mixed relaxed clock model               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   For example,                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      link shape=(all)                                                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   links the gamma/lnorm shape parameter across all partitions of the data.      \n");
        MrBayesPrint ("   You can use \"showmodel\" to see the current linking status of the            \n");
        MrBayesPrint ("   characters. For more information on this command, see the help menu           \n");
        MrBayesPrint ("   for link's converse, unlink (\"help unlink\");                                \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Help"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Help                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command provides useful information on the use of this program. The      \n");
        MrBayesPrint ("   correct usage is                                                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      help                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   which gives a list of all available commands with a brief description of      \n");
        MrBayesPrint ("   each or                                                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      help <command>                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   which gives detailed information on the use of <command>.                     \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Sump"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Sump                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   During an MCMC analysis, MrBayes prints the sampled parameter values to one or\n");
        MrBayesPrint ("   more tab-delimited text files, one for each independent run in your analysis. \n");
        MrBayesPrint ("   The command 'Sump' summarizes the information in this parameter file or these \n");
        MrBayesPrint ("   parameter files. By default, the root of the parameter file name(s) is assumed\n");
        MrBayesPrint ("   to be the name of the last matrix-containing nexus file. MrBayes also remem-  \n");
        MrBayesPrint ("   bers the number of independent runs in the last analysis that you set up, re- \n");
        MrBayesPrint ("   gardless of whether you actually ran it. For instance, if there were two in-  \n");
        MrBayesPrint ("   dependent runs, which is the initial setting when you read in a new matrix,   \n");
        MrBayesPrint ("   MrBayes will assume that there are two parameter files with the endings       \n");
        MrBayesPrint ("   '.run1.p' and '.run2.p'. You can change the root of the file names and the    \n");
        MrBayesPrint ("   number of runs using the 'Filename' and 'Nruns' settings.                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   When you invoke the 'Sump' command, three items are output: (1) a generation  \n");
        MrBayesPrint ("   plot of the likelihood values; (2) estimates of the marginal likelihood of    \n");
        MrBayesPrint ("   the model; and (3) a table with the mean, variance, and 95 percent credible   \n");
        MrBayesPrint ("   interval for the sampled parameters. All three items are output to screen.    \n");
        MrBayesPrint ("   The table of marginal likelihoods is also printed to a file with the ending   \n");
        MrBayesPrint ("   '.lstat' and the parameter table to a file with the ending '.pstat'. For some \n");
        MrBayesPrint ("   model parameters, there may also be a '.mstat' file.                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   When running 'Sump' you typically want to discard a specified number or       \n");
        MrBayesPrint ("   fraction of samples from the beginning of the chain as the burn in. This is   \n");
        MrBayesPrint ("   done using the same mechanism used by the 'mcmc' command. That is, if you     \n");
        MrBayesPrint ("   run an mcmc analysis with a relative burn in of 25 %% of samples for con-     \n");
        MrBayesPrint ("   vergence diagnostics, then the same burn in will be used for a subsequent     \n");
        MrBayesPrint ("   sump command, unless a different burn in is specified. That is, issuing       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   sump                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   immediately after 'mcmc', will result in using the same burn in settings as   \n");
        MrBayesPrint ("   for the 'mcmc' command. All burnin settings are reset to default values every \n");
        MrBayesPrint ("   time a new matrix is read in, namely relative burnin ('relburnin=yes') with   \n");
        MrBayesPrint ("   25 %% of samples discarded ('burninfrac = 0.25').                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Relburnin    -- If this option is set to 'Yes', then a proportion of the      \n");
        MrBayesPrint ("                   samples will be discarded as burnin when calculating summary  \n");
        MrBayesPrint ("                   statistics. The proportion to be discarded is set with        \n");
        MrBayesPrint ("                   'Burninfrac' (see below). When the 'Relburnin' option is set  \n");
        MrBayesPrint ("                   to 'No', then a specific number of samples is discarded       \n");
        MrBayesPrint ("                   instead. This number is set by 'Burnin' (see below). Note that\n");
        MrBayesPrint ("                   the burnin setting is shared across the 'sumt', 'sump', and   \n");
        MrBayesPrint ("                   'mcmc' commands.                                              \n");
        MrBayesPrint ("   Burnin       -- Determines the number of samples (not generations) that will  \n");
        MrBayesPrint ("                   be discarded when summary statistics are calculated. The      \n");
        MrBayesPrint ("                   value of this option is only applicable when 'Relburnin' is   \n");
        MrBayesPrint ("                   set to 'No'.                                                  \n");
        MrBayesPrint ("   Burninfrac   -- Determines the fraction of samples that will be discarded when\n");
        MrBayesPrint ("                   summary statistics are calculated. The setting only takes     \n");
        MrBayesPrint ("                   effect if 'Relburnin' is set to 'Yes'.                        \n");
        MrBayesPrint ("   Nruns        -- Determines how many '.p' files from independent analyses that \n");
        MrBayesPrint ("                   will be summarized. If Nruns > 1 then the names of the files  \n");
        MrBayesPrint ("                   are derived from 'Filename' by adding '.run1.p', '.run2.p',   \n");
        MrBayesPrint ("                   etc. If Nruns=1, then the single file name is obtained by     \n");
        MrBayesPrint ("                   adding '.p' to 'Filename'.                                    \n");
        MrBayesPrint ("   Filename     -- The name of the file to be summarized. This is the base of the\n");
        MrBayesPrint ("                   file name to which endings are added according to the current \n");
        MrBayesPrint ("                   setting of the 'Nruns' parameter. If 'Nruns' is 1, then only  \n");
        MrBayesPrint ("                   '.p' is added to the file name. Otherwise, the endings will   \n");
        MrBayesPrint ("                   be '.run1.p', '.run2.p', etc.                                 \n");
        MrBayesPrint ("   Outputname   -- Base name of the file(s) to which 'Sump' results will be      \n");
        MrBayesPrint ("                   printed.                                                      \n");
        MrBayesPrint ("   Hpd          -- Determines whether credibility intervals will be given as the \n");
        MrBayesPrint ("                   region of Highest Posterior Density ('Yes') or as the interval\n");
        MrBayesPrint ("                   containing the median 95 %% of sampled values ('No').         \n");
        MrBayesPrint ("   Minprob      -- Determines the minimum probability of submodels to be included\n");
        MrBayesPrint ("                   in summary statistics. Only applicable to models that explore \n");
        MrBayesPrint ("                   submodel spaces, like 'nst=mixed' and 'aamodelpr=mixed'.      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Current settings:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
        MrBayesPrint ("   --------------------------------------------------------                      \n");
        MrBayesPrint ("   Relburnin       Yes/No                   %s                                   \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        MrBayesPrint ("   Burnin          <number>                 %d                                   \n", chainParams.chainBurnIn);
        MrBayesPrint ("   Burninfrac      <number>                 %1.2lf                               \n", chainParams.burninFraction);
        MrBayesPrint ("   Nruns           <number>                 %d                                   \n", sumpParams.numRuns);
        if (sumpParams.numRuns == 1)
            MrBayesPrint ("   Filename        <name>                   %s<.p>\n", sumpParams.sumpFileName);
        else
            MrBayesPrint ("   Filename        <name>                   %s<.run<i>.p>\n", sumpParams.sumpFileName);
        MrBayesPrint ("   Outputname      <name>                   %s<.pstat etc>\n", sumpParams.sumpOutfile);
        MrBayesPrint ("   Hpd             Yes/No                   %s                                   \n", sumpParams.HPD == YES ? "Yes" : "No");
        MrBayesPrint ("   Minprob         <number>                 %1.3lf                               \n", sumpParams.minProb);
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Sumss"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Sumss                                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command summarizes results of stepping stone analyses. It is a tool to   \n");
        MrBayesPrint ("   investigate the obtained results, and to help find the proper step burn-in.   \n");
        MrBayesPrint ("   To get more help information on stepping-stone analyses, use 'help ss'.       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   During stepping-stone analysis, MrBayes collects the sampled likelihoods in   \n");
        MrBayesPrint ("   order to estimate the marginal likelihood at the end. It also prints the sam- \n");
        MrBayesPrint ("   pled parameter values to one or more tab-delimited text files, one for each   \n");
        MrBayesPrint ("   independent run in your analysis. The command 'Sumss' summarizes likelihood   \n");
        MrBayesPrint ("   values stored in these parameter files and calculates marginal likelihood es- \n");
        MrBayesPrint ("   timates. The names of the files that are summarized are exactly the same as   \n");
        MrBayesPrint ("   the names of the files used for the 'sump' command. In fact, the 'filename'   \n");
        MrBayesPrint ("   setting is a shared setting for the 'sump' and 'sumss' commands. That is, if  \n");
        MrBayesPrint ("   you change the setting in one of the commands, it would change the setting in \n");
        MrBayesPrint ("   the other command as well.                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   When you invoke the 'Sumss' command, three items are output: (1) 'Step contri-\n");
        MrBayesPrint ("   bution table' - summarizes the contribution of each step to the overall esti- \n");
        MrBayesPrint ("   mate; (2) 'Step plot' - plot of the likelihood values for the initial burn-in \n");
        MrBayesPrint ("   phase or a chosen step in the stepping-stone algorithm; (3) 'Joined plot' -   \n");
        MrBayesPrint ("   summarizes sampling across all steps in the algorithm.                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Step contribution table                                                       \n");
        MrBayesPrint ("   The printed table is similar to the one output to the .ss file. The main pur- \n");
        MrBayesPrint ("   pose of the table is to summarize marginal likelihood for different values of \n");
        MrBayesPrint ("   the step burn-in after the stepping stone  analysis has finished. The burn-in \n");
        MrBayesPrint ("   is controlled by the 'Relburnin', 'Burnin' and 'Burninfrac' settings.         \n");
        MrBayesPrint ("   Note that during stepping-stone analyses, step contributions to marginal      \n");
        MrBayesPrint ("   likelihood are calculated based on all generations excluding burn-in. 'Sumss' \n");
        MrBayesPrint ("   on the other hand makes estimates based only on the sampled generations. This \n");
        MrBayesPrint ("   may lead to slight difference in results compared to the one printed to the   \n");
        MrBayesPrint ("   .ss file.                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Step plot                                                                     \n");
        MrBayesPrint ("   The main objective of the plot is to provide a close look at a given step in  \n");
        MrBayesPrint ("   the analysis. Which step is printed here is defined by the 'Steptoplot' set-  \n");
        MrBayesPrint ("   ting. The plot could be used to inspect if the chosen step burn-in is appro-  \n");
        MrBayesPrint ("   priate for the given step. It could also be used to check if the initial burn-\n");
        MrBayesPrint ("   in phase has converged. Note that the amount of discarded samples is controled\n");
        MrBayesPrint ("   by the 'Discardfrac' setting, and not by the ordinary burn-in settings.       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Joined plot                                                                   \n");
        MrBayesPrint ("   Different steps sample from different power posterior distributions. When we  \n");
        MrBayesPrint ("   switch from one distribution to another, it takes some number of generations  \n");
        MrBayesPrint ("   before the chain settles at the correct stationary distribution. This lag is  \n");
        MrBayesPrint ("   called a 'temperature lag' and if the corresponding samples are not removed,  \n");
        MrBayesPrint ("   it will result in a biased estimate. It is difficult to determine the lag be- \n");
        MrBayesPrint ("   forehand, but MrBayes allows you to explore different step burn-in settings   \n");
        MrBayesPrint ("   after you have finished the stepping-stone algorithm, without having to rerun \n");
        MrBayesPrint ("   the whole analysis. The 'Joined plot' helps to facilitate the choice of the   \n");
        MrBayesPrint ("   right step burn-in. The plot summarizes samples across all steps and gives you\n");
        MrBayesPrint ("   a quick overview of the whole analysis.                                       \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Specifically, the following procedure is used to obtain the joined plot. Each \n");
        MrBayesPrint ("   step has the same number N of samples taken. We number each sample 1 to N     \n");
        MrBayesPrint ("   within steps according to the order in which the samples are taken. The first \n"); 
        MrBayesPrint ("   sample in each step is numbered 1, and the last sample is N. For each number i\n");
        MrBayesPrint ("   in [1,..., N], we sum up log likelihoods for all samples numbered i across all\n");
        MrBayesPrint ("   steps. The joined plot is a graph of the step number versus the normalized    \n");
        MrBayesPrint ("   sums we get in the procedure describe above. This directly visualizes the tem-\n");
        MrBayesPrint ("   perature lag and allows you to select the appropriate step burn-in.           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Ideally, after you discard the appropriate step burn-in, the graph should     \n");
        MrBayesPrint ("   appear as white noise around the estimated value. If you see an increasing or \n");
        MrBayesPrint ("   decreasing tendency in the beginning of the graph, you should increase the    \n");
        MrBayesPrint ("   step burn-in. If you see an increasing or decreasing tendency across the whole\n");
        MrBayesPrint ("   graph, then the initial burn-in phase was not long enough. In this case, you  \n");
        MrBayesPrint ("   need to rerun the analysis with a longer initial burn-in.                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   To make it easier to observe tendencies in the plotted graph you can choose   \n");
        MrBayesPrint ("   different levels of curve smoothing. If 'Smoothing' is set to k, it means that\n");
        MrBayesPrint ("   for each step i we take an average over step i and k neighboring samples in   \n");
        MrBayesPrint ("   both directions, i.e., the k-smoothed estimate for step i is an average over  \n");
        MrBayesPrint ("   values for steps [i-k,...,i+k].                                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Allruns      -- If set to 'Yes', it forces all runs to be printed on the same \n");
        MrBayesPrint ("                   graph when drawing joined and step plots. If set to 'No', each\n");
        MrBayesPrint ("                   run is printed on a separat plot.                             \n");
        MrBayesPrint ("   Askmore      -- Long analyses may produce huge .p files. Reading in them may  \n");
        MrBayesPrint ("                   take several minutes. If you want to investigate different    \n");
        MrBayesPrint ("                   aspects of your analyses, it could be very inconvenient to    \n");
        MrBayesPrint ("                   wait for several minutes each time you want to get a new sum- \n");
        MrBayesPrint ("                   mary for different settings. If you set 'Askmore' to 'YES',   \n");
        MrBayesPrint ("                   sumss will read .p files only once. After responding to the   \n");
        MrBayesPrint ("                   original query, it will interactivaly ask you if you wish to  \n");
        MrBayesPrint ("                   produce more tables and plots for different settings of       \n");
        MrBayesPrint ("                   'Burnin' or 'Smoothing' (see below).                          \n");
        MrBayesPrint ("   Relburnin    -- If this option is set to 'Yes', then a proportion of the      \n");
        MrBayesPrint ("                   samples from each step will be discarded as burnin when calcu-\n");
        MrBayesPrint ("                   lsting summary statistics. The proportion to be discarded is  \n");
        MrBayesPrint ("                   set with 'Burninfrac' (see below). When the 'Relburnin' option\n");
        MrBayesPrint ("                   is set to 'No', then a specific number of samples is discarded\n");
        MrBayesPrint ("                   instead. This number is set by 'Burnin'. Note that the burnin \n");
        MrBayesPrint ("                   settings --- 'Relburnin', 'Burnin', and 'Burninfrac' --- are  \n");
        MrBayesPrint ("                   shared across the 'sumt', 'sump', 'sumss' and 'mcmc' commands.\n");
        MrBayesPrint ("   Burnin       -- Determines the number of samples (not generations) that will  \n");
        MrBayesPrint ("                   be discarded from each step when summary statistics are calcu-\n");
        MrBayesPrint ("                   lated. The value of this option is only applicable when       \n");
        MrBayesPrint ("                   'Relburnin' is set to 'No'.                                   \n");
        MrBayesPrint ("   Burninfrac   -- Determines the fraction of samples that will be discarded from\n");
        MrBayesPrint ("                   each step when summary statistics are calculated. The setting \n");
        MrBayesPrint ("                   only takes effect if 'Relburnin' is set to 'Yes'.             \n");
        MrBayesPrint ("   Discardfrac  -- Determines the fraction of samples that will be discarded when\n");
        MrBayesPrint ("                   a step plot is printed. It is similar to the 'Burninfrac' set-\n");
        MrBayesPrint ("                   ting, but unlike 'Burninfrac' it is used only for better vis- \n");
        MrBayesPrint ("                   ualization of the step plot. It has no effect on the number of\n");
        MrBayesPrint ("                   samples discarded during marginal likelihood computation.     \n");
        MrBayesPrint ("   Filename     -- The name of the file to be summarized. This is the base of the\n");
        MrBayesPrint ("                   file name to which endings are added according to the current \n");
        MrBayesPrint ("                   setting of the 'Nruns' parameter. If 'Nruns' is 1, then only  \n");
        MrBayesPrint ("                   '.p' is added to the file name. Otherwise, the endings will   \n");
        MrBayesPrint ("                   be '.run1.p', '.run2.p', etc. Note that the 'Filename' setting\n");
        MrBayesPrint ("                   is shared with 'sump' command.                                \n");
        MrBayesPrint ("   Nruns        -- Determines how many '.p' files from independent analyses that \n");
        MrBayesPrint ("                   will be summarized. If Nruns > 1 then the names of the files  \n");
        MrBayesPrint ("                   are derived from 'Filename' by adding '.run1.p', '.run2.p',   \n");
        MrBayesPrint ("                   etc. If Nruns=1, then the single file name is obtained by     \n");
        MrBayesPrint ("                   adding '.p' to 'Filename'.                                    \n");
        MrBayesPrint ("   Steptoplot   -- Defines which step will be printed in the step plot.If the    \n");
        MrBayesPrint ("                   value is set to 0, then the initial sample from the posterior \n");
        MrBayesPrint ("                   will be used.                                                 \n");
        MrBayesPrint ("   Smoothing    -- Determines smoothing of the joined plot (see above). A value  \n");
        MrBayesPrint ("                   equal to 0 results in no smoothing.                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Current settings:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
        MrBayesPrint ("   --------------------------------------------------------                      \n");
        MrBayesPrint ("   Allruns         Yes/No                   %s                                   \n", sumssParams.allRuns == YES ? "Yes" : "No");
        MrBayesPrint ("   Askmore         Yes/No                   %s                                   \n", sumssParams.askForMorePlots == YES ? "Yes" : "No");
        MrBayesPrint ("   Relburnin       Yes/No                   %s                                   \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        MrBayesPrint ("   Burnin          <number>                 %d                                   \n", chainParams.chainBurnIn);
        MrBayesPrint ("   Burninfrac      <number>                 %1.2lf                               \n", chainParams.burninFraction);
        MrBayesPrint ("   Discardfrac     <number>                 %1.2lf                               \n", sumssParams.discardFraction);
        if (sumpParams.numRuns == 1)
            MrBayesPrint ("   Filename        <name>                   %s<.p>\n", sumpParams.sumpFileName);
        else
            MrBayesPrint ("   Filename        <name>                   %s<.run<i>.p>\n", sumpParams.sumpFileName);        
        MrBayesPrint ("   Nruns           <number>                 %d                                   \n", sumpParams.numRuns);
        MrBayesPrint ("   Steptoplot      <number>                 %d                                   \n", sumssParams.stepToPlot);
        MrBayesPrint ("   Smoothing       <number>                 %d                                   \n", sumssParams.smoothing);
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Comparetree"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Comparetree                                                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command compares the trees in two files, called \"filename1\" and        \n");
        MrBayesPrint ("   \"filename2\". It will output a bivariate plot of the split frequencies       \n");
        MrBayesPrint ("   as well as plots of the tree distance as a function of the generation. The    \n");
        MrBayesPrint ("   plots can be used to get a quick indication of whether two runs have con-     \n");
        MrBayesPrint ("   verged onto the same set of trees. The \"Comparetree\" command will also      \n");
        MrBayesPrint ("   produce a \".pairs\" file and a \".dists\" file (these file endings are added \n");
        MrBayesPrint ("   to the end of the \"Outputname\"). The \".pairs\" file contains the paired    \n");
        MrBayesPrint ("   split frequencies from the two tree samples; the \".dists\" file contains the \n");
        MrBayesPrint ("   tree distance values.                                                         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Note that the \"Sumt\" command provides a different set of convergence diag-  \n");
        MrBayesPrint ("   nostics tools that you may also want to explore. Unlike \"Comparetree\",      \n");
        MrBayesPrint ("   \"Sumt\" can compare more than two tree samples and will calculate consensus  \n");
        MrBayesPrint ("   trees and split frequencies from the pooled samples.                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Relburnin     -- If this option is set to 'Yes', then a proportion of the     \n");
        MrBayesPrint ("                    samples will be discarded as burnin when calculating summary \n");
        MrBayesPrint ("                    statistics. The proportion to be discarded is set with       \n");
        MrBayesPrint ("                    Burninfrac (see below). When the Relburnin option is set to  \n");
        MrBayesPrint ("                    'No', then a specific number of samples is discarded instead.\n");
        MrBayesPrint ("                    This number is set by Burnin (see below). Note that the      \n");
        MrBayesPrint ("                    burnin setting is shared with the 'mcmc', 'sumt', 'sump' and \n");
        MrBayesPrint ("                    'plot' commands.                                             \n");
        MrBayesPrint ("   Burnin        -- Determines the number of samples (not generations) that will \n");
        MrBayesPrint ("                    be discarded when summary statistics are calculated. The     \n");
        MrBayesPrint ("                    value of this option is only relevant when Relburnin is set  \n");
        MrBayesPrint ("                    to 'No'.                                                     \n");
        MrBayesPrint ("   BurninFrac    -- Determines the fraction of samples that will be discarded    \n");
        MrBayesPrint ("                    when summary statistics are calculated. The value of this    \n");
        MrBayesPrint ("                    option is only relevant when Relburnin is set to 'Yes'.      \n");
        MrBayesPrint ("                    Example: A value for this option of 0.25 means that 25%% of  \n");
        MrBayesPrint ("                    the samples will be discarded.                               \n");
        MrBayesPrint ("   Minpartfreq   -- The minimum probability of partitions to include in summary  \n");
        MrBayesPrint ("                    statistics.                                                  \n");
        MrBayesPrint ("   Filename1     -- The name of the first tree file to compare.                  \n");
        MrBayesPrint ("   Filename2     -- The name of the second tree file to compare.                 \n");
        MrBayesPrint ("   Outputname    -- Name of the file to which 'comparetree' results will be      \n");
        MrBayesPrint ("                    printed.                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Current settings:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
        MrBayesPrint ("   --------------------------------------------------------                      \n");
        MrBayesPrint ("   Relburnin       Yes/No                   %s                                   \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        MrBayesPrint ("   Burnin          <number>                 %d                                   \n", chainParams.chainBurnIn);
        MrBayesPrint ("   Burninfrac      <number>                 %1.2lf                               \n", chainParams.burninFraction);
        MrBayesPrint ("   Minpartfreq     <number>                 %1.2lf                               \n", comptreeParams.minPartFreq);
        MrBayesPrint ("   Filename1       <name>                   %s                                   \n", comptreeParams.comptFileName1);
        MrBayesPrint ("   Filename2       <name>                   %s                                   \n", comptreeParams.comptFileName2);
        MrBayesPrint ("   Outputname      <name>                   %s                                   \n", comptreeParams.comptOutfile);
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Sumt"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Sumt                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used to produce summary statistics for trees sampled during   \n");
        MrBayesPrint ("   a Bayesian MCMC analysis. You can either summarize trees from one individual  \n");
        MrBayesPrint ("   analysis, or trees coming from several independent analyses. In either case,  \n");
        MrBayesPrint ("   all the sampled trees are read in and the proportion of the time any single   \n");
        MrBayesPrint ("   taxon bipartition (split) is found is counted. The proportion of the time that\n");
        MrBayesPrint ("   the bipartition is found is an approximation of the posterior probability of  \n");
        MrBayesPrint ("   the bipartition. (Remember that a taxon bipartition is defined by removing a  \n");
        MrBayesPrint ("   branch on the tree, dividing the tree into those taxa to the left and right   \n");
        MrBayesPrint ("   of the removed branch. This set is called a taxon bipartition.) The branch    \n");
        MrBayesPrint ("   length of the bipartition is also recorded, if branch lengths have been saved \n");
        MrBayesPrint ("   to file. The result is a list of the taxon bipartitions found, the frequency  \n");
        MrBayesPrint ("   with which they were found, the posterior probability of the bipartition      \n");
        MrBayesPrint ("   and, the mean and variance of the branch lengths or node depths, and various  \n");
        MrBayesPrint ("   other statistics.                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The key to the partitions is output to a file with the suffix '.parts'. The   \n");
        MrBayesPrint ("   summary statistics pertaining to bipartition probabilities are output to a    \n");
        MrBayesPrint ("   file with the suffix '.tstat', and the statistics pertaining to branch or node\n");
        MrBayesPrint ("   parameters are output to a file with the suffix '.vstat'.                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   A consensus tree is also printed to a file with the suffix '.con.tre' and     \n");
        MrBayesPrint ("   printed to the screen as a cladogram, and as a phylogram if branch lengths    \n");
        MrBayesPrint ("   have been saved. The consensus tree is either a 50 percent majority rule tree \n");
        MrBayesPrint ("   or a majority rule tree showing all compatible partitions. If branch lengths  \n");
        MrBayesPrint ("   have been recorded during the run, the '.con.tre' file will contain a consen- \n");
        MrBayesPrint ("   sus tree with branch lengths and interior nodes labelled with support values. \n");
        MrBayesPrint ("   By default, the consensus tree will also contain other summary information in \n");
        MrBayesPrint ("   a format understood by the program 'FigTree'. To use a simpler format under-  \n");
        MrBayesPrint ("   stood by other tree-drawing programs, such as 'TreeView', set 'Conformat' to  \n");
        MrBayesPrint ("   'Simple'.                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   MrBayes alo produces a file with the ending \".trprobs\" that contains a list \n");
        MrBayesPrint ("   of all the trees that were found during the MCMC analysis, sorted by their    \n");
        MrBayesPrint ("   probabilities. This list of trees can be used to construct a credible set of  \n");
        MrBayesPrint ("   trees. For example, if you want to construct a 95 percent credible set of     \n");
        MrBayesPrint ("   trees, you include all of those trees whose cumulative probability is less    \n");
        MrBayesPrint ("   than or equal to 0.95. You have the option of displaying the trees to the     \n");
        MrBayesPrint ("   screen using the \"Showtreeprobs\" option. The default is to not display the  \n");
        MrBayesPrint ("   trees to the screen; the number of different trees sampled by the chain can   \n");
        MrBayesPrint ("   be quite large. If you are analyzing a large set of taxa, you may actually    \n");
        MrBayesPrint ("   want to skip the calculation of tree probabilities entirely by setting        \n");
        MrBayesPrint ("   'Calctreeprobs' to 'No'.                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   When calculating summary statistics you probably want to skip those trees that\n");
        MrBayesPrint ("   were sampled in the initial part of the run, the so-called burn-in period. The\n");
        MrBayesPrint ("   number of skipped samples is controlled by the 'Relburnin', 'Burnin', and     \n");
        MrBayesPrint ("   'Burninfrac' settings, just as for the 'Mcmc' command. Since version 3.2.0,   \n");
        MrBayesPrint ("   the burn-in settings are shared across the 'Sumt', 'Sump' and 'Mcmc' commands.\n");
        MrBayesPrint ("   That is, changing the burn-in setting for one command will change the settings\n");
        MrBayesPrint ("   for subsequent calls to any of the other commands.                            \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   If you are summarizing the trees sampled in several independent analyses,     \n");
        MrBayesPrint ("   such as those resulting from setting the 'Nruns' option of the 'Mcmc' command \n");
        MrBayesPrint ("   to a value larger than 1, MrBayes will also calculate convergence diagnostics \n");
        MrBayesPrint ("   for the sampled topologies and branch lengths. These values can help you      \n");
        MrBayesPrint ("   determine whether it is likely that your chains have converged.               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The 'Sumt' command expands the 'Filename' according to the current values of  \n");
        MrBayesPrint ("   the 'Nruns' and 'Ntrees' options. For instance, if both 'Nruns' and 'Ntrees'  \n");
        MrBayesPrint ("   are set to 1, 'Sumt' will try to open a file named '<Filename>.t'. If 'Nruns' \n");
        MrBayesPrint ("   is set to 2 and 'Ntrees' to 1, then 'Sumt' will open two files, the first     \n");
        MrBayesPrint ("   named '<Filename>.run1.t' and the second '<Filename>.run2.t', etc. By default,\n");
        MrBayesPrint ("   the 'Filename' option is set such that 'Sumt' automatically summarizes all the\n");
        MrBayesPrint ("   results from your immediately preceding 'Mcmc' command. You can also use the  \n");
        MrBayesPrint ("   'Sumt' command to summarize tree samples in older analyses. If you want to do \n");
        MrBayesPrint ("   that, remember to first read in a matrix so that MrBayes knows what taxon     \n");
        MrBayesPrint ("   names to expect in the trees. Then set the 'Nruns', 'Ntrees' and 'Filename'   \n");
        MrBayesPrint ("   options appropriately if they differ from the MrBayes defaults.               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Relburnin     -- If this option is set to YES, then a proportion of the       \n");
        MrBayesPrint ("                    samples will be discarded as burnin when calculating summary \n");
        MrBayesPrint ("                    statistics. The proportion to be discarded is set with       \n");
        MrBayesPrint ("                    Burninfrac (see below). When the Relburnin option is set to  \n");
        MrBayesPrint ("                    NO, then a specific number of samples is discarded instead.  \n");
        MrBayesPrint ("                    This number is set by Burnin (see below). Note that the      \n");
        MrBayesPrint ("                    burnin setting is shared across the 'sumt', 'sump', and      \n");
        MrBayesPrint ("                    'mcmc' commands.                                             \n");
        MrBayesPrint ("   Burnin        -- Determines the number of samples (not generations) that will \n");
        MrBayesPrint ("                    be discarded when summary statistics are calculated. The     \n");
        MrBayesPrint ("                    value of this option is only relevant when Relburnin is set  \n");
        MrBayesPrint ("                    to NO.                                                       \n");
        MrBayesPrint ("   BurninFrac    -- Determines the fraction of samples that will be discarded    \n");
        MrBayesPrint ("                    when summary statistics are calculated. The value of this    \n");
        MrBayesPrint ("                    option is only relevant when Relburnin is set to YES.        \n");
        MrBayesPrint ("                    Example: A value for this option of 0.25 means that 25%% of  \n");
        MrBayesPrint ("                    the samples will be discarded.                               \n");
        MrBayesPrint ("   Nruns         -- Determines how many '.t' files from independent analyses that\n");
        MrBayesPrint ("                    will be summarized. If Nruns > 1 then the names of the files \n");
        MrBayesPrint ("                    are derived from 'Filename' by adding '.run1.t', '.run2.t',  \n");
        MrBayesPrint ("                    etc. If Nruns=1 and Ntrees=1 (see below), then only '.t' is  \n");
        MrBayesPrint ("                    added to 'Filename'.                                         \n");
        MrBayesPrint ("   Ntrees        -- Determines how many trees there are in the sampled model. If \n");
        MrBayesPrint ("                    'Ntrees' > 1 then the names of the files are derived from    \n");
        MrBayesPrint ("                    'Filename' by adding '.tree1.t', '.tree2.t', etc. If there   \n");
        MrBayesPrint ("                    are both multiple trees and multiple runs, the filenames will\n");
        MrBayesPrint ("                    be '<Filename>.tree1.run1.t', '<Filename>.tree1.run2.t', etc.\n");
        MrBayesPrint ("   Filename      -- The name of the file(s) to be summarized. This is the base of\n");
        MrBayesPrint ("                    the file name, to which endings are added according to the   \n");
        MrBayesPrint ("                    current settings of the 'Nruns' and 'Ntrees' options.        \n");
        MrBayesPrint ("   Minpartfreq   -- The minimum probability of partitions to include in summary  \n");
        MrBayesPrint ("                    statistics.                                                  \n");
        MrBayesPrint ("   Contype       -- Type of consensus tree. 'Halfcompat' results in a 50%% major-\n");
        MrBayesPrint ("                    ity rule tree, 'Allcompat' adds all compatible groups to such\n");
        MrBayesPrint ("                    a tree.                                                      \n");
        MrBayesPrint ("   Conformat     -- Format of consensus tree. The 'Figtree' setting results in a \n");
        MrBayesPrint ("                    consensus tree formatted for the program FigTree, with rich  \n");
        MrBayesPrint ("                    summary statistics. The 'Simple' setting results in a simple \n");
        MrBayesPrint ("                    consensus tree written in a format read by a variety of pro- \n");
        MrBayesPrint ("                    grams.                                                       \n");
        MrBayesPrint ("   Outputname    -- Base name of the file(s) to which 'sumt' results will be     \n");
        MrBayesPrint ("                    printed. The default is the same as 'Filename'.              \n");
        MrBayesPrint ("   Calctreeprobs -- Determines whether tree probabilities should be calculated.  \n");
        MrBayesPrint ("   Showtreeprobs -- Determines whether tree probabilities should be displayed on \n");
        MrBayesPrint ("                    screen.                                                      \n");
        MrBayesPrint ("   Hpd           -- Determines whether credibility intervals will be given as the\n");
        MrBayesPrint ("                    region of Highest Posterior Density ('Yes') or as the inter- \n");
        MrBayesPrint ("                    val containing the median 95 %% of sampled values ('No').    \n");
        MrBayesPrint ("   Savebrparams  -- Set this option to 'yes' to save all sampled branch and node \n");
        MrBayesPrint ("                    parameter values to a separate file with the filename ending \n");
        MrBayesPrint ("                    in '.brparams'. All partitions with a posterior probability  \n");
        MrBayesPrint ("                    larger than Minbrparamfreq will be included.                 \n");
        MrBayesPrint ("   Minbrparamfreq -- The minimum probability of partitions for which to save     \n");
        MrBayesPrint ("                     parameter values to file if 'Savebrparams' is set to 'yes'. \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Current settings:                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
        MrBayesPrint ("   --------------------------------------------------------                      \n");
        MrBayesPrint ("   Relburnin       Yes/No                   %s                                   \n", chainParams.relativeBurnin == YES ? "Yes" : "No");
        MrBayesPrint ("   Burnin          <number>                 %d                                   \n", chainParams.chainBurnIn);
        MrBayesPrint ("   Burninfrac      <number>                 %1.2lf                               \n", chainParams.burninFraction);
        MrBayesPrint ("   Nruns           <number>                 %d                                   \n", sumtParams.numRuns);
        MrBayesPrint ("   Ntrees          <number>                 %d                                   \n", sumtParams.numTrees);
        if (sumtParams.numRuns == 1 && sumtParams.numTrees == 1)
            MrBayesPrint ("   Filename        <name>                   %s<.t>\n", sumtParams.sumtFileName);
        else if (sumtParams.numRuns == 1 && sumtParams.numTrees > 1)
            MrBayesPrint ("   Filename        <name>                   %s<.tree<i>.t>\n", sumtParams.sumtFileName);
        else if (sumtParams.numRuns > 1 && sumtParams.numTrees == 1)
            MrBayesPrint ("   Filename        <name>                   %s<.run<i>.t>\n", sumtParams.sumtFileName);
        else if (sumtParams.numRuns > 1 && sumtParams.numTrees > 1)
            MrBayesPrint ("   Filename        <name>                   %s<.tree<i>.run<i>.t>\n", sumtParams.sumtFileName);
        MrBayesPrint ("   Minpartfreq     <number>                 %1.2lf                               \n", sumtParams.minPartFreq);
        MrBayesPrint ("   Contype         Halfcompat/Allcompat     %s\n", sumtParams.sumtConType);
        MrBayesPrint ("   Conformat       Figtree/Simple           %s                                   \n", sumtParams.consensusFormat == SIMPLE ? "Simple" : "Figtree");
        MrBayesPrint ("   Outputname      <name>                   %s<.parts etc>\n", sumtParams.sumtOutfile);
        MrBayesPrint ("   Calctreeprobs   Yes/No                   %s                                   \n", sumtParams.calcTreeprobs == YES ? "Yes" : "No");
        MrBayesPrint ("   Showtreeprobs   Yes/No                   %s                                   \n", sumtParams.showSumtTrees == YES ? "Yes" : "No");
        MrBayesPrint ("   Hpd             Yes/No                   %s                                   \n", sumtParams.HPD == YES ? "Yes" : "No");
        MrBayesPrint ("   Savebrparams    Yes/No                   %s                                   \n", sumtParams.saveBrParams == YES ? "Yes" : "No");
        MrBayesPrint ("   Minbrparamfreq  <number>                 %1.2lf                               \n", sumtParams.minBrParamFreq);
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Tree"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Tree                                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command is used by MrBayes to write trees to a nexus tree file. Trees    \n");
        MrBayesPrint ("   are written in the Newick format. For instance,                               \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      tree ((1,2),3,4);                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   describes an unrooted tree with taxa 1 and 2 being more closely related to    \n");
        MrBayesPrint ("   each other than to taxa 3 and 4. If branch lengths are saved to file, they    \n");
        MrBayesPrint ("   are given after a colon sign immediately following the terminal taxon or the  \n");
        MrBayesPrint ("   interior node they refer to. An example of an unrooted tree with branch       \n");
        MrBayesPrint ("   lengths is:                                                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      tree ((1:0.064573,2:0.029042):0.041239,3:0.203988,4:0.187654);             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Trees that are rooted (clock trees) are written with a basal dichotomy        \n");
        MrBayesPrint ("   instead of a basal trichotomy. If the tree described above had been rooted    \n");
        MrBayesPrint ("   on the branch leading to taxon 4, it would have been represented as:          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      tree (((1,2),3),4);                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Report"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Report                                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command allows you to control how the posterior distribution is          \n");
        MrBayesPrint ("   reported. For rate parameters, it allows you to choose among several popular  \n");
        MrBayesPrint ("   parameterizations. The report command also allows you to request printing of  \n");
        MrBayesPrint ("   some model aspects that are usually not reported. For instance, if a node is  \n");
        MrBayesPrint ("   constrained in the analysis, MrBayes can print the probabilities of the       \n");
        MrBayesPrint ("   ancestral states at that node. Similarly, if there is rate variation in the   \n");
        MrBayesPrint ("   model, MrBayes can print the inferred site rates, and if there is omega varia-\n");
        MrBayesPrint ("   tion, MrBayes can print the inferred omega (positive selection) values for    \n");
        MrBayesPrint ("   each codon. In a complex model with several partitions, each partition is     \n");
        MrBayesPrint ("   controlled separately using the same 'Applyto' mechanism as in the 'Lset' and \n");
        MrBayesPrint ("   'Prset' commands.                                                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Options:                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Applyto   -- This option allows you to apply the report commands to specific  \n");
        MrBayesPrint ("                partitions. This command should be the first in the list of      \n");
        MrBayesPrint ("                commands specified in 'report'.                                  \n");
        MrBayesPrint ("                For example,                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                   report applyto=(1,2) tratio=ratio                             \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                   report applyto=(3) tratio=dirichlet                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("                would result in the transition and transversion rates of the     \n");
        MrBayesPrint ("                first and second partitions in the model being reported as a     \n");
        MrBayesPrint ("                ratio and the transition and transversion rates of the third     \n");
        MrBayesPrint ("                partition being reported as proportions of the rate sum (the     \n");
        MrBayesPrint ("                Dirichlet parameterization).                                     \n");
        MrBayesPrint ("   Tratio    -- This specifies the report format for the transition and trans-   \n");
        MrBayesPrint ("                version rates of a nucleotide substituion model with nst=2.      \n");
        MrBayesPrint ("                If 'ratio' is selected, the rates will be reported as a ratio    \n");
        MrBayesPrint ("                (transition rate/transversion rate). If 'dirichlet' is selected, \n");
        MrBayesPrint ("                the transition and transversion rates will instead be reported   \n");
        MrBayesPrint ("                as proportions of the rate sum. For example, if the transition   \n");
        MrBayesPrint ("                rate is three times the transversion rate and 'ratio' is selec-  \n");
        MrBayesPrint ("                ted, this will reported as a single value, '3.0'. If 'dirichlet' \n");
        MrBayesPrint ("                is selected instead, the same rates will be reported using two   \n");
        MrBayesPrint ("                values, '0.75 0.25'. The sum of the Dirichlet values is always 1.\n");
        MrBayesPrint ("                Although the Dirichlet format may be unfamiliar to some users,   \n");
        MrBayesPrint ("                it is more convenient for specifying priors than the ratio       \n");
        MrBayesPrint ("                format.                                                          \n");
        MrBayesPrint ("   Revmat    -- This specifies the report format for the substitution rates of   \n");
        MrBayesPrint ("                a GTR substitution model for nucleotide or amino acid data. If   \n");
        MrBayesPrint ("                'ratio' is selected, the rates will be reported scaled to the    \n");
        MrBayesPrint ("                G-T rate (for nucleotides) or the Y-V rate (for amino acids). If \n");
        MrBayesPrint ("                'dirichlet' is specified instead, the rates are reported as pro- \n");
        MrBayesPrint ("                portions of the rate sum. For instance, assume that the C-T rate \n");
        MrBayesPrint ("                is twice the A-G rate and four times the transversion rates,     \n");
        MrBayesPrint ("                which are equal. If the report format is set to 'ratio', this    \n");
        MrBayesPrint ("                would be reported as '1.0 2.0 1.0 1.0 4.0 1.0' since the rates   \n");
        MrBayesPrint ("                are reported in the order rAC, rAG, rAT, rCG, rCT, rGT and scaled\n");
        MrBayesPrint ("                relative to the last rate, the G-T rate. If 'dirichlet' is selec-\n");
        MrBayesPrint ("                ted instead, the same rates would have been reported as '0.1 0.2 \n");
        MrBayesPrint ("                0.1 0.1 0.4 0.1' since the rates are now scaled so that they sum \n");
        MrBayesPrint ("                to 1.0. The Dirichlet format is the parameterization used for    \n");
        MrBayesPrint ("                formulating priors on the rates.                                 \n");
        MrBayesPrint ("   Ratemult  -- This specifies the report format used for the rate multiplier of \n");
        MrBayesPrint ("                different model partitions. Three formats are available. If      \n");
        MrBayesPrint ("                'scaled' is selected, then rates are scaled such that the mean   \n");
        MrBayesPrint ("                rate per site across partitions is 1.0. If 'ratio' is chosen,    \n");
        MrBayesPrint ("                the rates are scaled relative to the rate of the first parti-    \n");
        MrBayesPrint ("                tion. Finally, if 'dirichlet' is chosen, the rates are given as  \n");
        MrBayesPrint ("                proportions of the rate sum. The latter is the format used       \n");
        MrBayesPrint ("                when formulating priors on the rate multiplier.                  \n");
        MrBayesPrint ("   Tree      -- This specifies the report format used for the tree(s). Two op-   \n");
        MrBayesPrint ("                tions are available. 'Topology' results in only the topology     \n");
        MrBayesPrint ("                being printed to file, whereas 'brlens' causes branch lengths to \n");
        MrBayesPrint ("                to be printed as well.                                           \n");
        MrBayesPrint ("   Ancstates -- If this option is set to 'yes', MrBayes will print the pro-      \n");
        MrBayesPrint ("                bability of the ancestral states at all constrained nodes. Typ-  \n");
        MrBayesPrint ("                ically, you are interested in the ancestral states of only a few \n");
        MrBayesPrint ("                characters and only at one node in the tree. To perform such     \n");
        MrBayesPrint ("                an analysis, first define and enforce a topology constraint      \n");
        MrBayesPrint ("                using 'constraint' and 'prset topologypr = constraints (...)'.   \n");
        MrBayesPrint ("                Then put the character(s) of interest in a separate partition and\n");
        MrBayesPrint ("                set MrBayes to report the ancestral states for that partition.   \n");
        MrBayesPrint ("                For instance, if the characters of interest are in partition 2,  \n");
        MrBayesPrint ("                use 'report applyto=(2) ancstates=yes' to force MrBayes to print \n");
        MrBayesPrint ("                the probability of the ancestral states of those characters at   \n");
        MrBayesPrint ("                the constrained node to the '.p' file.                           \n");
        MrBayesPrint ("   Siterates -- If this option is set to 'yes' and the relevant model has rate   \n");
        MrBayesPrint ("                variation across sites, then the site rates, weighted over rate  \n");
        MrBayesPrint ("                categories, will be reported to the '.p' file.                   \n");
        MrBayesPrint ("   Possel    -- If this option is set to 'yes' and the relevant model has omega  \n");
        MrBayesPrint ("                variation across sites, the probability that each model site     \n");
        MrBayesPrint ("                (codon in this case) is positively selected will be written to   \n");
        MrBayesPrint ("                file.                                                            \n");
        MrBayesPrint ("   Siteomega -- If this option is set to 'yes' and the relevant model has omega  \n");
        MrBayesPrint ("                variation across sites, the weighted omega value (over omega     \n");
        MrBayesPrint ("                categories) for each model site will be reported to file.        \n");
        MrBayesPrint ("                                                                                 \n");
        if (numCurrentDivisions == 0)
            tempInt = 1;
        else
            tempInt = numCurrentDivisions;
        for (i=0; i<tempInt; i++)
            {
            if (numCurrentDivisions == 0)
                {
                MrBayesPrint ("   Default report settings:                                                       \n");
                mp = &defaultModel;
                }
            else
                {
                MrBayesPrint ("   Current report settings for partition %d:                                              \n", i+1);
                mp = &modelParams[i];
                }
            MrBayesPrint ("                                                                                 \n");
            MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
            MrBayesPrint ("   --------------------------------------------------------                      \n");
            MrBayesPrint ("   Tratio          Ratio/Dirichlet          %s                                   \n", mp->tratioFormat);
            MrBayesPrint ("   Revmat          Ratio/Dirichlet          %s                                   \n", mp->revmatFormat);
            MrBayesPrint ("   Ratemult        Scaled/Ratio/Dirichlet   %s                                   \n", mp->ratemultFormat);
            MrBayesPrint ("   Tree            Brlens/Topology          %s                                   \n", mp->treeFormat);
            MrBayesPrint ("   Ancstates       Yes/No                   %s                                   \n", mp->inferAncStates);
            MrBayesPrint ("   Siterates       Yes/No                   %s                                   \n", mp->inferSiteRates);
            MrBayesPrint ("   Possel          Yes/No                   %s                                   \n", mp->inferPosSel);
            MrBayesPrint ("   Siteomega       Yes/No                   %s                                   \n", mp->inferSiteOmegas);
            MrBayesPrint ("                                                                                 \n");
            MrBayesPrint ("   ------------------------------------------------------------------            \n");       
            }
        }
    else if (!strcmp(helpTkn, "Manual"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Manual                                                                        \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command allows you to generate a text file containing help information   \n");
        MrBayesPrint ("   on all the available commands. This text file can be used as an up-to-date    \n");
        MrBayesPrint ("   command reference. You can set the name of the text file using the            \n");
        MrBayesPrint ("   \"filename\" option; the default is \"commref_mb<version>.txt\".              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Parameter       Options                  Current Setting                      \n");
        MrBayesPrint ("   --------------------------------------------------------                      \n");
        MrBayesPrint ("   Filename        <name>                   %s                                   \n", manFileName);
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showmoves"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Showmoves                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows the MCMC samplers (moves) that are switched on for the     \n");
        MrBayesPrint ("   parameters in the current model. The basic usage is                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      showmoves                                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   If you want to see all available moves, use                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      showmoves allavailable=yes                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   If you want to change any of the tuning parameters for the moves, use the     \n");
        MrBayesPrint ("   'propset' command.                                                            \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Showparams"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Showparams                                                                    \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This command shows all of the parameters in the current model. The basic      \n");
        MrBayesPrint ("   usage is                                                                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      showparams                                                                 \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The parameters are listed together with their priors, the available moves,    \n");
        MrBayesPrint ("   and the current value(s), which will be used as the starting values in the    \n");
        MrBayesPrint ("   next mcmc analysis.                                                           \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else if (!strcmp(helpTkn, "Startvals"))
        {
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        MrBayesPrint ("   Startvals                                                                     \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   Use this command to change the current values for parameters in your model.   \n");
        MrBayesPrint ("   These values will be used as the starting values in the next mcmc analysis.   \n");
        MrBayesPrint ("   The basic format is:                                                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      startvals <param>=(<value_1>,<value_2>,...,<value_n>)                      \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   for all substitution model parameters. The format is slightly different for   \n");
        MrBayesPrint ("   parameters that are written to a tree file:                                   \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      startvals <param>=<tree_name>                                              \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   This version of the command will look for a tree with the specified name      \n");
        MrBayesPrint ("   among the trees read in previously when parsing a tree block. The information \n");
        MrBayesPrint ("   stored in that tree will be used to set the starting value of the parameter.  \n");
        MrBayesPrint ("   The parameters that are set using this mechanism include topology and branch  \n");
        MrBayesPrint ("   length parameters, as well as relaxed clock branch rates, cpp events and      \n");
        MrBayesPrint ("   cpp branch rate multipliers.                                                  \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   The above versions of the command will set the value for all runs and chains. \n");
        MrBayesPrint ("   You can also set the value for an individual run and chain by using the format\n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      startvals <param>(<run>,<chain>)=(<value_1>,...)                           \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   where <run> is the index of the run and <chain> the index of the chain. If    \n");
        MrBayesPrint ("   the run index is omitted, the values will be changed for all runs. Similarly, \n");
        MrBayesPrint ("   if the chain index is omitted, all chains will be set to the specified value. \n");
        MrBayesPrint ("   For example, if we wanted to set the values of the stationary frequency       \n");
        MrBayesPrint ("   parameter pi{1} to (0.1,0.1,0.4,0.4) for all chains in run 1, and to          \n");
        MrBayesPrint ("   (0.3,0.3,0.2,0.2) for chain 3 of run 2, we would use                          \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("      startvals pi{1}(1,)=(0.1,0.1,0.4,0.4) pi{1}(2,3)=(0.3,0.3,0.2,0.2)         \n");
        MrBayesPrint ("                                                                                 \n");
        MrBayesPrint ("   ---------------------------------------------------------------------------   \n");
        }
    else
        {
        return (ERROR);
        }
        
    return (NO_ERROR);
}


/* IsAmbig: This function returns YES if character is set as ambiguous
   either by using parenthetic notation or by ambiguity codes. It returns
   NO if character is unambiguous, missing or gapped */ 
int IsAmbig (int charCode, int dType)
{
    if (dType == DNA || dType == RNA || dType == STANDARD || dType == RESTRICTION || dType == PROTEIN)
        {
        if (charCode != MISSING && charCode != GAP)
            if (NBits(charCode) > 1)
                return (YES);
        }
    else if (dType == CONTINUOUS)
        {
        /* do nothing, these cannot be partly ambiguous */
        }
    else
        {
        MrBayesPrint ("Unknown datatype in \"IsAmbig\"\n", spacer);
        }

    return (NO);
}


int IsArgValid (char *tk, char *validArg)
{
    int         i, j, k, tkLen, targetLen, numDiff, numStrMatches;
    char        tempStr[100];
    ParmInfoPtr p;

    p = paramPtr;
    tkLen = (int) strlen(tk);

    numStrMatches = i = j = 0;
    do
        {
        if (p->valueList[i] == '|' || p->valueList[i] == '\0')
            {
            tempStr[j++] = '\0';
            targetLen = (int) strlen(tempStr);
            if (tkLen <= targetLen)
                {
                numDiff = 0;
                for (k=0; k<tkLen; k++)
                    if (ChangeCase(tk[k]) != ChangeCase(tempStr[k]))
                        numDiff++;
                if (numDiff == 0)
                    {
                    numStrMatches++;
                    strcpy (validArg, tempStr);
                    }
                }
            j = 0;
            }
        else
            tempStr[j++] = p->valueList[i];
        i++;
        }
    while (p->valueList[i] != '\0');
        
    if (numStrMatches == 0)
        {
        MrBayesPrint ("%s   No valid match for argument \"%s\"\n", spacer, tk);
        return (ERROR);
        }
    else if (numStrMatches == 1)
        {
        return (NO_ERROR);
        }
    else
        {
        MrBayesPrint ("%s   Argument \"%s\" is ambiguous\n", spacer, tk);
        return (ERROR);
        }
}


int IsIn (char ch, char *s)
{
    while (*s)
        {
        if (*s++ == ch)
            return 1;
        }
    return 0;
}


int IsMissing (int charCode, int dType)
{
    if (dType == DNA || dType == RNA)
        {
        if (charCode == 15 || charCode == 16)
            return (YES);
        }
    else if (dType == STANDARD || dType == PROTEIN)
        {
        if (charCode == MISSING || charCode == GAP)
            return (YES);
        }
    else if (dType == RESTRICTION)
        {
        if (charCode == 3 || charCode == 4)
            return (YES);
        }
    else if (dType == CONTINUOUS)
        {

        }
    else
        {
        MrBayesPrint ("Unknown datatype in \"IsMissing\"\n", spacer);
        }
    return (NO);
}


int IsSame (char *s1, char *s2)
{
    int         i, nDiff, isIdentical, len;
    
    isIdentical = YES;
    if (strlen(s1) != strlen(s2))
        isIdentical = NO; /* strings cannot be identical because they are different lengths */
    
    /* now, we go through both strings, one character at a time, to see if
       any are different */
    if (strlen(s1) > strlen(s2))
        len = (int) strlen(s2);
    else
        len = (int) strlen(s1);
    i = nDiff = 0;
    while (i < len)
        {
        if (tolower(s1[i]) != tolower(s2[i]))
            nDiff++;
        i++;
        }
    if (nDiff == 0 && isIdentical == YES)
        return (SAME);
    else if (nDiff == 0 && isIdentical == NO)
        return (CONSISTENT_WITH);
    else
        return (DIFFERENT);
}


int IsWhite (char c)
{
    if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
        {
        if (c == '\n' || c == '\r')
            return 2;
        return 1;
        }
    return 0;
}


int NucID (char nuc)
{
    char        n;
    
    if (nuc == 'U' || nuc == 'u')
        n = 'T';
    else
        n = nuc;

    if (n == 'A' || n == 'a')
        {
        return 1;
        }
    else if (n == 'C' || n == 'c')
        {
        return 2;
        }
    else if (n == 'G' || n == 'g')
        {
        return 4;
        }
    else if (n == 'T' || n == 't')
        {
        return 8;
        }
    else if (n == 'R' || n == 'r')
        {
        return 5;
        }
    else if (n == 'Y' || n == 'y')
        {
        return 10;
        }
    else if (n == 'M' || n == 'm')
        {
        return 3;
        }
    else if (n == 'K' || n == 'k')
        {
        return 12;
        }
    else if (n == 'S' || n == 's')
        {
        return 6;
        }
    else if (n == 'W' || n == 'w')
        {
        return 9;
        }
    else if (n == 'H' || n == 'h')
        {
        return 11;
        }
    else if (n == 'B' || n == 'b')
        {
        return 14;
        }
    else if (n == 'V' || n == 'v')
        {
        return 7;
        }
    else if (n == 'D' || n == 'd')
        {
        return 13;
        }
    else if (n == 'N' || n == 'n')
        {
        return 15;
        }
    else if (n == gapId)
        {
        return GAP;
        }
    else if (n == missingId)
        {
        return MISSING;
        }
    else
        return -1;
}


/*-------| ParseCommand |------------------------------------------------
|
|   This function is used to parse a file. The expected format is:
|   
|      command parameter=value parameter=value ... ;
|
|   For example, the following is a valid line for this parser:
|
|      lset nst=2;
|
|   In some cases, however, the format is:
|
|      command stuff more stuff ... ;
|
|   For example, when reading a data file, the matrix command might be:
|
|      matrix
|         taxon_1 data
|         taxon_2 data
|         taxon_3 data
|         ;
|
|   Like before, the command and all of the stuff for that command are
|   terminated by a semicolon.
|
*/
int ParseCommand (char *s)
{
    int             rc, tokenType, inError, numMatches, skipCmd;
    char            errStr[100];

    numMatches = 0;     /* Avoid gcc warnings (actually set in call to FindValidCommand) */
    cmdStr = s;
    tokenP = &s[0];

#   if defined (ECHO_PROCESSED_COMMANDS)
        MrBayesPrint ("Currently processing command: %s\n", s);
#   endif
    
    inError = skipCmd = NO;
    do
        {
        /* Get the next token. A token is a valid word in a line. Token type is defined in "bayes.h". */
        if (GetToken (token, &tokenType, &tokenP))
            {
            inError = YES; 
            break;
            }
        if (strlen(token) > 0 || tokenType == ALPHA)
            {
#           if defined (SHOW_TOKENS)
            MrBayesPrint ("%s\n", token);
#           endif
            if (tokenType == LEFTCOMMENT)
                {
                /* If the token is a left comment "[", then we don't want to
                   actually process commands until we find a right comment.  */
                /* The exception is if readComment is set to YES, in which case
                   we will leave it to the parser functions to decide on whether
                   they want to read the comment or not */
                if (readComment == NO || inComment == YES)
                    {
                    inComment = YES;
                    numComments++;
                    }
                }
            if (inComment == NO && inForeignBlock == NO)
                {
                if (tokenType != SEMICOLON)
                    {
                    /* If the token is not a semicolon, then we will be processing 
                       either a command or a parameter. */
                    if (expecting == Expecting(COMMAND))
                        {
                        /* We are expecting to find a command (defined above in "commands[]"). Find the 
                           correct command and set a pointer to that command. */
                        commandPtr = NULL;
                        if (FindValidCommand (token, &numMatches) == ERROR)
                            {
                            /* We couldn't find the command or the user did not specify enough letters
                               to unambiguously determine the command. The command pointer (commandPtr)
                               is NULL. */
                            if (numMatches == 0)    
                                MrBayesPrint ("%s   Could not find command \"%s\"\n", spacer, token);
                            else 
                                MrBayesPrint ("%s   Ambiguous command \"%s\"\n", spacer, token);
                            inError = YES;
                            }
                        else
                            {
                            /* We did find a valid command. Set what we are expecting to see next. */
                            expecting = commandPtr->expect;
                            
                            /* Check to see if we have one of the so-called special cases in which a 
                               command is not necessarily followed by a parameter (e.g., matrix). If we
                               do have a special case, then we want to set the parameter pointer (paramPtr)
                               appropriately. In this case, simply go to the first parameter in the parmList. */
                            if (commandPtr->specialCmd == YES)
                                {
                                isFirstMatrixRead = YES;
                                foundFirst = NO;
                                paramPtr = paramTable + commandPtr->parmList[0];
                                }
                            if (strcmp(commandPtr->string, "Execute")==0)
                                {
                                /* set the tokenizer to recognize quoted strings */
                                readWord = YES;
                                }
                            }
                        }
                    else 
                        {
                        /* We are expecting to find a parameter or a value for the parameter, not a command. */
                        if ((expecting & Expecting(PARAMETER)) == Expecting(PARAMETER) && 
                            (expecting & Expecting(tokenType)) != Expecting(tokenType))
                            {
                            /* Specifically, if we are here, we need to go through the parameter list,
                               checking to see if the token is a valid parameter. */
                            expecting = (expecting & Expecting(PARAMETER));
                            if (FindValidParam (token, &numMatches) == ERROR)
                                {
                                /* The token is not a valid parameter. */
                                if (numMatches == 0)
                                    MrBayesPrint ("%s   Could not find parameter \"%s\"\n", spacer, token);
                                else 
                                    MrBayesPrint ("%s   Ambiguous parameter \"%s\"\n", spacer, token);
                                inError = YES;
                                }
                            else
                                {
                                /* The token is a valid parameter. Call the appropriate function ("DoXxxxParm"). */
                                if ((paramPtr->fp)(paramPtr->string, token) == ERROR)
                                    {
                                    if (strcmp("Xxxxxxxxxx", paramPtr->string))
                                        MrBayesPrint ("%s   Error when setting parameter \"%s\" (1)\n", spacer, paramPtr->string);
                                    inError = YES;
                                    }
                                }
                            }
                        else
                            {
                            /* Otherwise, we are expecting a value for the parameter. Call the appropriate function ("DoXxxxParm"). */
                            if ((expecting & Expecting(tokenType)) != 0)
                                expecting = (expecting & Expecting(tokenType));
                            if ((Expecting(tokenType) & expecting) == Expecting(tokenType))
                                {
                                if ((paramPtr->fp)(paramPtr->string, token) == ERROR)
                                    {
                                    if (strcmp("Xxxxxxxxxx", paramPtr->string))
                                        MrBayesPrint ("%s   Error when setting parameter \"%s\" (2)\n", spacer, paramPtr->string);
                                    inError = YES;
                                    }
                                }
                            else
                                {
                                inError = YES;
                                WhatVariableExp (expecting, errStr);
                                MrBayesPrint ("%s   Expecting '%s'\n", spacer, errStr+1);  /* there will be an initial space in errStr so print from pos 1 */
                                if (numOpenExeFiles > 0)
                                    MrBayesPrint ("%s   Instead found '%s' in command '%s'\n",
                                        spacer, token, commandPtr->string);
                                else
                                    MrBayesPrint ("%s   Instead found '%s' in command '%s' at position %d\n",
                                        spacer, token, commandPtr->string, tokenP - cmdStr - strlen(token)+1);
                                }
                            }
                        }
                    }
                else
                    {
                    /* The token is a semicolon. This means that we are at the end of processing one command. We
                       need to clean things up. We do this by calling the finishing function ("DoXxxx"). */
                    if ((Expecting(SEMICOLON) & expecting) == Expecting(SEMICOLON))
                        {
                        if (commandPtr->cmdFxnPtr != NULL)
                            {
                            /* Finish up the command here. */
                            rc = (commandPtr->cmdFxnPtr) ();
                            if (rc  == ERROR || rc == ABORT)
                                {
                                if (rc == ABORT)
                                    {
                                    MrBayesPrint ("   Mcmc run aborted\n");
                                    }
                                else if (rc == SKIP_COMMAND)
                                    {
                                    MrBayesPrint ("   Cancelled execution of command\n");
                                    skipCmd = YES;
                                    }
                                else
                                    {
                                    MrBayesPrint ("%s   Error in command \"%s\"\n", spacer, commandPtr->string);
                                    inError = YES;
                                    }
                                }
                            }
                        /* if the user typed "quit", then we want to bail out of this loop, with a NO_ERROR_QUIT */
                        if (!strcmp(commandPtr->string, "Quit"))
                            return (NO_ERROR_QUIT);
                        expecting = Expecting(COMMAND);
                        }
                    else
                        {
                        inError = YES;
                        WhatVariableExp (expecting, errStr);
                        MrBayesPrint ("%s   Expecting %s\n", spacer, errStr);
                        }
                    }
                }
            /* Check to see if a comment is terminated. A comment can either be a right comment "]" or, if we were in a foreign nexus block
               (e.g., a "paup" block) the terminating comment will be "end". */
            if (tokenType == RIGHTCOMMENT)
                {
                if (inComment == NO && readComment == NO)
                    {
                    MrBayesPrint ("%s   Found \"]\", without having previously found \"[\"\n", spacer);
                    inError = YES; 
                    }
                else if (inComment == NO && readComment == YES)
                    {
                    /* This is OK, we just pass through and rely on the command to handle the RIGHTCOMMENT */
                    }
                else
                    {
                    numComments--;
                    if (numComments == 0)
                        inComment = NO;
                    }
                }
            if ((IsSame(token, "end") == SAME || IsSame(token, "endblock") == SAME) && inForeignBlock == YES)
                {
                strcpy (spacer, "");
                inForeignBlock = NO;
                }
            }
        
        } while ((*token || tokenType == ALPHA) && inError == NO && skipCmd == NO);
        
    if (inError == YES)
        {
        readComment = NO;   /* reset this in case it is set to YES in command and we get an error exit */
        return (ERROR);
        }
    else
        return (NO_ERROR);
}


void PrintSettings (char *command)
{
    char yesNoStr[20];

    if (!strcmp(command,"Mcmc"))
        {
        MrBayesPrint ("   Parameter       Options               Current Setting                         \n");
        MrBayesPrint ("   -----------------------------------------------------                         \n");
        MrBayesPrint ("   Ngen            <number>              %d                                      \n", chainParams.numGen);
        MrBayesPrint ("   Nruns           <number>              %d                                      \n", chainParams.numRuns);
        MrBayesPrint ("   Nchains         <number>              %d                                      \n", chainParams.numChains);
        MrBayesPrint ("   Temp            <number>              %lf                                     \n", chainParams.chainTemp);
        MrBayesPrint ("   Reweight        <number>,<number>     %1.2lf v %1.2lf ^                       \n", chainParams.weightScheme[0], chainParams.weightScheme[1]);
        MrBayesPrint ("   Swapfreq        <number>              %d                                      \n", chainParams.swapFreq);
        MrBayesPrint ("   Nswaps          <number>              %d                                      \n", chainParams.numSwaps);
        MrBayesPrint ("   Samplefreq      <number>              %d                                      \n", chainParams.sampleFreq);
        MrBayesPrint ("   Printfreq       <number>              %d                                      \n", chainParams.printFreq);
        PrintYesNo (chainParams.printAll, yesNoStr);
        MrBayesPrint ("   Printall        Yes/No                %s                                      \n", yesNoStr);
        MrBayesPrint ("   Printmax        <number>              %d                                      \n", chainParams.printMax);
        PrintYesNo (chainParams.mcmcDiagn, yesNoStr);
        MrBayesPrint ("   Mcmcdiagn       Yes/No                %s                                      \n", yesNoStr);
        MrBayesPrint ("   Diagnfreq       <number>              %d                                      \n", chainParams.diagnFreq);
        if (chainParams.diagnStat == AVGSTDDEV)
            strcpy (yesNoStr, "Avgstddev");
        else
            strcpy (yesNoStr, "Maxstddev");
        MrBayesPrint ("   Diagnstat       Avgstddev/Maxstddev   %s                                     \n", yesNoStr);
        MrBayesPrint ("   Minpartfreq     <number>              %1.2lf                                 \n", chainParams.minPartFreq);
        PrintYesNo (chainParams.allChains, yesNoStr);
        MrBayesPrint ("   Allchains       Yes/No                %s                                     \n", yesNoStr);
        PrintYesNo (chainParams.allComps, yesNoStr);
        MrBayesPrint ("   Allcomps        Yes/No                %s                                     \n", yesNoStr);
        PrintYesNo (chainParams.relativeBurnin, yesNoStr);
        MrBayesPrint ("   Relburnin       Yes/No                %s                                     \n", yesNoStr);
        MrBayesPrint ("   Burnin          <number>              %d                                     \n", chainParams.chainBurnIn);
        MrBayesPrint ("   Burninfrac      <number>              %1.2lf                                 \n", chainParams.burninFraction);
        PrintYesNo (chainParams.stopRule, yesNoStr);
        MrBayesPrint ("   Stoprule        Yes/No                %s                                     \n", yesNoStr);
        MrBayesPrint ("   Stopval         <number>              %1.2lf                                 \n", chainParams.stopVal);
        PrintYesNo (chainParams.saveTrees, yesNoStr);
        MrBayesPrint ("   Savetrees       Yes/No                %s                                     \n", yesNoStr);
        PrintYesNo (chainParams.checkPoint, yesNoStr);
        MrBayesPrint ("   Checkpoint      Yes/No                %s                                     \n", yesNoStr);
        MrBayesPrint ("   Checkfreq       <number>              %d                                     \n", chainParams.checkFreq);
        MrBayesPrint ("   Filename        <name>                %s.<p/t>\n", chainParams.chainFileName);
        MrBayesPrint ("   Startparams     Current/Reset         %s                                     \n", chainParams.startParams);
        MrBayesPrint ("   Starttree       Current/Random/       %s                                     \n", chainParams.startTree);
        MrBayesPrint ("                   Parsimony                                                    \n");
        MrBayesPrint ("   Nperts          <number>              %d                                     \n", chainParams.numStartPerts);
        PrintYesNo (chainParams.runWithData, yesNoStr);
        MrBayesPrint ("   Data            Yes/No                %s                                     \n", yesNoStr);
        MrBayesPrint ("   Ordertaxa       Yes/No                %s                                     \n", chainParams.orderTaxa == YES? "Yes" : "No");
        MrBayesPrint ("   Append          Yes/No                %s                                     \n", chainParams.append == YES? "Yes" : "No");
        MrBayesPrint ("   Autotune        Yes/No                %s                                     \n", chainParams.autotune == YES? "Yes" : "No");
        MrBayesPrint ("   Tunefreq        <number>              %d                                     \n", chainParams.tuneFreq);
        MrBayesPrint ("                                                                                \n");
        }
}


void PrintYesNo (int yn, char s[4])
{
    if (yn == YES)
        strcpy (s, "Yes");
    else
        strcpy (s, "No");
}


int ProtID (char aa)
{
    if (aa == 'A' || aa == 'a')      /* Ala */
        {
        return 1;
        }
    else if (aa == 'R' || aa == 'r') /* Arg */
        {
        return 2;
        }
    else if (aa == 'N' || aa == 'n') /* Asn */
        {
        return 4;
        }
    else if (aa == 'D' || aa == 'd') /* Asp */
        {
        return 8;
        }
    else if (aa == 'C' || aa == 'c') /* Cys */
        {
        return 16;
        }
    else if (aa == 'Q' || aa == 'q') /* Gln */
        {
        return 32;
        }
    else if (aa == 'E' || aa == 'e') /* Glu */
        {
        return 64;
        }
    else if (aa == 'G' || aa == 'g') /* Gly */
        {
        return 128;
        }
    else if (aa == 'H' || aa == 'h') /* His */
        {
        return 256;
        }
    else if (aa == 'I' || aa == 'i') /* Ile */
        {
        return 512;
        }
    else if (aa == 'L' || aa == 'l') /* Leu */
        {
        return 1024;
        }
    else if (aa == 'K' || aa == 'k') /* Lys */
        {
        return 2048;
        }
    else if (aa == 'M' || aa == 'm') /* Met */
        {
        return 4096;
        }
    else if (aa == 'F' || aa == 'f') /* Phe */
        {
        return 8192;
        }
    else if (aa == 'P' || aa == 'p') /* Pro */
        {
        return 16384;
        }
    else if (aa == 'S' || aa == 's') /* Ser */
        {
        return 32768;
        }
    else if (aa == 'T' || aa == 't') /* Thr */
        {
        return 65536;
        }
    else if (aa == 'W' || aa == 'w') /* Trp */
        {
        return 131072;
        }
    else if (aa == 'Y' || aa == 'y') /* Tyr */
        {
        return 262144;
        }
    else if (aa == 'V' || aa == 'v') /* Val */
        {
        return 524288;
        }
    else if (aa == 'X' || aa == 'x') /* Nonidentified */
        {
        return MISSING;
        }
    else if (aa == gapId)
        {
        return GAP;
        }
    else if (aa == missingId)
        {
        return MISSING;
        }
    else
        return -1;
}


int RemoveLastFromString (char *s1)
{
    int     i, j, numPrev, numRemoved;
    
    /* We remove the last name from the string simply by deleting the last "|". */
       
    i = numPrev = 0;
    while (s1[i] != '\0')
        {
        if (s1[i] == '|')
            numPrev++;
        i++;
        }
        
    i = j = numRemoved = 0;
    while (s1[i] != '\0')
        {
        if (s1[i] == '|')
            j++;
        if (numPrev == j)
            {
            s1[i] = ' ';
            numRemoved++;
            break;
            }
        i++;
        }

    if (numRemoved != 1)
        {
        MrBayesPrint ("%s   Could not find name to remove\n", spacer);
        return (ERROR);
        }

    return (NO_ERROR);
}


int MBResID (char nuc)
{
    char        n;
    
    n = nuc;

    if (n == '0' || n == 'a' || n == 'A')
        {
        return 1;
        }
    else if (n == '1' || n == 'b' || n == 'B')
        {
        return 2;
        }
    else if (n == gapId)
        {
        return GAP;
        }
    else if (n == missingId)
        {
        return MISSING;
        }
    else
        return -1;
}


/* Reset character flags */
void ResetCharacterFlags (void)
{
    /* reset all characters flags */
    numChar              = 0;                        /* number of defined characters                  */
    defChars             = NO;                       /* flag for whether number of characters is known*/
    defMatrix            = NO;                       /* flag for whether matrix is successfull read   */
    matrixHasPoly        = NO;                       /* flag for whether matrix has polymorphisms     */
    isInAmbig            = NO;                       /* flag for whether the parser is within ()      */
    isInPoly             = NO;                       /* flag for whether the parser is within {}      */
    defPartition         = NO;                       /* flag for whether character partition is read  */
    defPairs             = NO;                       /* flag indicating whether pairs have been defnd */
    numDefinedPartitions = 0;                        /* number of defined partitions                  */
    partitionNum         = 0;                        /* partition number currently enforced           */
    numCurrentDivisions  = 0;                        /* number of partitions of data                  */
    numCharSets          = 0;                        /* holds number of character sets                */
    numDivisions         = 1;                        /* holds number of partitions                    */
    isMixed              = NO;                       /* are data mixed ?                              */
    dataType             = NONE;                     /* holds datatype                                */
    matchId              = '\0';                     /* no default for match character                */
    gapId                = '\0';                     /* no default for gap character                  */
    missingId            = '\0';                     /* no default for missing characters             */
}


/* Reset taxa flags */
void ResetTaxaFlags (void)
{
    numTaxa                 = 0;                         /* number of taxa                                */
    numNamedTaxa            = 0;                         /* number of named taxa                          */
    defTaxa                 = NO;                        /* flag for whether number of taxa is known      */
    isTaxsetDef             = NO;                        /* is a taxlabels set defined                    */
    numDefinedConstraints   = 0;                         /* holds number of defined constraints           */
    definedConstraint       = NULL;
    definedConstraintTwo    = NULL;
    definedConstraintPruned       = NULL;
    definedConstraintTwoPruned    = NULL;
    constraintNames         = NULL;
    nodeCalibration         = NULL;
    tempActiveConstraints   = NULL;                      /* holds temp info on active constraints         */
    outGroupNum             = 0;                         /* default outgroup                              */
    numTaxaSets             = 0;                         /* holds number of taxa sets                     */
}


/* SetPartition: Set model partition */
int SetPartition (int part)
{
    int     i, j;
    
    /* Free space for modelParams and modelSettings */
    if (memAllocs[ALLOC_MODEL] == YES)
        {
        for (i=0; i<numCurrentDivisions; i++)
          free (modelParams[i].activeConstraints);
        free (modelParams);
        free (modelSettings);
        modelParams = NULL;
        modelSettings = NULL;
        memAllocs[ALLOC_MODEL] = NO;
        }

    /* Set model partition */
    partitionNum = part;
    numCurrentDivisions = 0;

    /* Set numCurrentDivisions to maximum division a character belongs to in partition part */
    for (i=0; i<numChar; i++)
        {
        j = partitionId[i][part];
        if (j > numCurrentDivisions)
            numCurrentDivisions = j;
        }

    /* Allocate space for partition models */
    modelParams = (Model *) SafeCalloc (numCurrentDivisions, sizeof (Model));
    modelSettings = (ModelInfo *) SafeCalloc (numCurrentDivisions, sizeof (ModelInfo));
    if (!modelParams || !modelSettings)
        {
        MrBayesPrint ("%s   Could not allocate modelParams or modelSettings\n", spacer);
        if (modelParams)
            free (modelParams);
        if (modelSettings)
            free (modelSettings);
        return (ERROR);
        }
    memAllocs[ALLOC_MODEL] = YES;

    numVars = (int *) SafeRealloc ((void *) numVars, 3 * (size_t)numCurrentDivisions * sizeof(int));
    tempLinkUnlinkVec = numVars + numCurrentDivisions;
    activeParts       = numVars + 2*numCurrentDivisions;

    tempNum = (MrBFlt *) SafeRealloc ((void *) tempNum, 6 * sizeof(MrBFlt));

    activeParams[0] = (int *) SafeRealloc ((void *) (activeParams[0]), (size_t)NUM_LINKED * (size_t)numCurrentDivisions * sizeof(int));
    for (i=1; i<NUM_LINKED; i++)
        activeParams[i] = activeParams[0] + i*numCurrentDivisions;
 
    linkTable[0] = (int *) SafeRealloc ((void *) (linkTable[0]), 3 * (size_t)NUM_LINKED * (size_t)numCurrentDivisions * sizeof(int));
    tempLinkUnlink[0] = linkTable[0] + NUM_LINKED*numCurrentDivisions;
    for (i=1; i<NUM_LINKED; i++)
        {
        linkTable[i]      = linkTable[0] + i*numCurrentDivisions;
        tempLinkUnlink[i] = tempLinkUnlink[0] + i*numCurrentDivisions;
        }

    return (NO_ERROR);
}


/* SetSpeciespartition: Set speciespartition */
int SetSpeciespartition (int part)
{
    int     i, j;
    
    /* Set model partition */
    speciespartitionNum = part;
    numSpecies = 0;

    /* Set numSpecies to maximum species a taxon belongs to in partition part */
    for (i=0; i<numTaxa; i++)
        {
        j = speciespartitionId[i][part];
        if (j > numSpecies)
            numSpecies = j;
        }

    return (NO_ERROR);
}


int SetTaxaFromTranslateTable (void)
{
    int     i;

    if (numTaxa != 0)
        return ERROR;

    for (i=0; i<numTranslates; i++)
        {
        if (strlen(transFrom[i])>99)
            {
            MrBayesPrint ("%s   Taxon name %s is too long. Maximun 99 characters is allowed.\n", spacer, transFrom[i]);
            return (ERROR);
            }
        AddString(&taxaNames, numTaxa, transFrom[i]);
        numTaxa++;
        }
    
    return NO_ERROR;
}


void SetUpParms (void)
{
    ParmInfoPtr p = paramTable;

    PARAM   (0, "NEXUS",          DoNexusParm,       "NEXUS|\0");
    PARAM   (1, "Data",           DoBeginParm,       "\0");
    PARAM   (2, "Mrbayes",        DoBeginParm,       "\0");
    PARAM   (3, "Trees",          DoBeginParm,       "\0");
    PARAM   (4, "Ntax",           DoDimensionsParm,  "\0");
    PARAM   (5, "Nchar",          DoDimensionsParm,  "\0");
    PARAM   (6, "Interleave",     DoFormatParm,      "Yes|No|\0");
    PARAM   (7, "Datatype",       DoFormatParm,      "Dna|Rna|Protein|Restriction|Standard|Continuous|Mixed|\0");
    PARAM   (8, "Gap",            DoFormatParm,      "\0");
    PARAM   (9, "Missing",        DoFormatParm,      "\0");
    PARAM  (10, "Matchchar",      DoFormatParm,      "\0");
    PARAM  (11, "MatrixInfo",     DoMatrixParm,      "\0");
    PARAM  (12, "Filename",       DoExecuteParm,     "\0");
    PARAM  (13, "Autoclose",      DoSetParm,         "Yes|No|\0");
    PARAM  (14, "Partition",      DoSetParm,         "\0");
    PARAM  (15, "Xxxxxxxxxx",     DoCharsetParm,     "\0");
    PARAM  (16, "Xxxxxxxxxx",     DoPartitionParm,   "\0");
    PARAM  (17, "Seed",           DoMcmcParm,        "\0");
    PARAM  (18, "Ngen",           DoMcmcParm,        "\0");
    PARAM  (19, "Samplefreq",     DoMcmcParm,        "\0");
    PARAM  (20, "Printfreq",      DoMcmcParm,        "\0");
    PARAM  (21, "Nchains",        DoMcmcParm,        "\0");
    PARAM  (22, "Temp",           DoMcmcParm,        "\0");
    PARAM  (23, "Filename",       DoMcmcParm,        "\0");
    PARAM  (24, "Burnin",         DoMcmcParm,        "\0");
    PARAM  (25, "Starttree",      DoMcmcParm,        "Random|Current|User|Parsimony|NJ|\0");
    PARAM  (26, "Nperts",         DoMcmcParm,        "\0");
    PARAM  (27, "Savebrlens",     DoMcmcParm,        "Yes|No|\0");
    PARAM  (28, "Nucmodel",       DoLsetParm,        "4by4|Doublet|Codon|Protein|\0");
    PARAM  (29, "Nst",            DoLsetParm,        "1|2|6|Mixed|\0");
    PARAM  (30, "Aamodel",        DoLsetParm,        "Poisson|Equalin|Jones|Dayhoff|Mtrev|Mtmam|Wag|Rtrev|Cprev|Vt|Blosum|Blossum|LG|\0");
    PARAM  (31, "Parsmodel",      DoLsetParm,        "Yes|No|\0");
    PARAM  (32, "Omegavar",       DoLsetParm,        "Equal|Ny98|M3|M10|\0");
    PARAM  (33, "Code",           DoLsetParm,        "Universal|Vertmt|Invermt|Mycoplasma|Yeast|Ciliate|Echinoderm|Euplotid|Metmt|\0");
    PARAM  (34, "Coding",         DoLsetParm,        "All|Variable|Informative|Nosingletons|Noabsencesites|Nopresencesites|Nosingletonpresence|Nosingletonabsence|\0");
    PARAM  (35, "Seqerror",       DoPrsetParm,       "\0");
    PARAM  (36, "Tratiopr",       DoPrsetParm,       "Beta|Fixed|\0");
    PARAM  (37, "Revmatpr",       DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM  (38, "Omegapr",        DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM  (39, "Statefreqpr",    DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM  (40, "Ngammacat",      DoLsetParm,        "\0");
    PARAM  (41, "Shapepr",        DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM  (42, "Ratecorrpr",     DoPrsetParm,       "Uniform|Fixed|\0");
    PARAM  (43, "Pinvarpr",       DoPrsetParm,       "Uniform|Fixed|\0");
    PARAM  (44, "Covswitchpr",    DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM  (45, "Xxxxxxxxxx",     DoExcludeParm,     "\0");
    PARAM  (46, "Xxxxxxxxxx",     DoIncludeParm,     "\0");
    PARAM  (47, "Xxxxxxxxxx",     DoDeleteParm,      "\0");
    PARAM  (48, "Xxxxxxxxxx",     DoRestoreParm,     "\0");
    PARAM  (49, "Xxxxxxxxxx",     DoTaxasetParm,     "\0");
    PARAM  (50, "Xxxxxxxxxx",     DoHelpParm,        "\0");
    PARAM  (51, "Applyto",        DoLsetParm,        "\0");
    PARAM  (52, "Rates",          DoLsetParm,        "Equal|Gamma|LNorm|Propinv|Invgamma|Adgamma|Kmixture|\0");
    PARAM  (53, "Covarion",       DoLsetParm,        "Yes|No|\0");
    PARAM  (54, "Applyto",        DoPrsetParm,       "\0");
    PARAM  (55, "Tratio",         DoLinkParm,        "\0");
    PARAM  (56, "Revmat",         DoLinkParm,        "\0");
    PARAM  (57, "Omega",          DoLinkParm,        "\0");
    PARAM  (58, "Statefreq",      DoLinkParm,        "\0");
    PARAM  (59, "Shape",          DoLinkParm,        "\0");
    PARAM  (60, "Pinvar",         DoLinkParm,        "\0");
    PARAM  (61, "Correlation",    DoLinkParm,        "\0");
    PARAM  (62, "Ratemultiplier", DoLinkParm,        "\0");
    PARAM  (63, "Switchrates",    DoLinkParm,        "\0");
    PARAM  (64, "Symdirihyperpr", DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM  (65, "Xxxxxxxxxx",     DoCtypeParm,       "\0");
    PARAM  (66, "Xxxxxxxxxx",     DoConstraintParm,  "\0");
    PARAM  (67, "Topologypr",     DoPrsetParm,       "Uniform|Constraints|Fixed|Speciestree|\0");
    PARAM  (68, "Brlenspr",       DoPrsetParm,       "Unconstrained|Clock|Relaxedclock|Fixed|\0");
    PARAM  (69, "Speciationpr",   DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM  (70, "Extinctionpr",   DoPrsetParm,       "Beta|Fixed|\0");
    PARAM  (71, "Popsizepr",      DoPrsetParm,       "Lognormal|Uniform|Gamma|Normal|Fixed|\0");
    PARAM  (72, "Topology",       DoLinkParm,        "\0");
    PARAM  (73, "Brlens",         DoLinkParm,        "\0");
    PARAM  (74, "Speciationrate", DoLinkParm,        "\0");
    PARAM  (75, "Extinctionrate", DoLinkParm,        "\0");
    PARAM  (76, "Popsize",        DoLinkParm,        "\0");
    PARAM  (77, "Ratepr",         DoPrsetParm,       "Variable|Dirichlet|Fixed|\0");
    PARAM  (78, "Xxxxxxxxxx",     DoOutgroupParm,    "\0");
    PARAM  (79, "Xxxxxxxxxx",     DoTreeParm,        "\0");
    PARAM  (80, "Filename",       DoSumtParm,        "\0");
    PARAM  (81, "Burnin",         DoSumtParm,        "\0");
    PARAM  (82, "Contype",        DoSumtParm,        "Halfcompat|Allcompat|\0");
    PARAM  (83, "Xxxxxxxxxx",     DoTranslateParm,   "\0");
    PARAM  (84, "Swapfreq",       DoMcmcParm,        "\0");
    PARAM  (85, "Start",          DoLogParm,         "\0");
    PARAM  (86, "Stop",           DoLogParm,         "\0");
    PARAM  (87, "Filename",       DoLogParm,         "\0");
    PARAM  (88, "Append",         DoLogParm,         "\0");
    PARAM  (89, "Replace",        DoLogParm,         "\0");
    PARAM  (90, "Nbetacat",       DoLsetParm,        "\0");
    PARAM  (91, "Augment",        DoLsetParm,        "Yes|No|\0");
    PARAM  (92, "Xxxxxxxxxx",     DoPairsParm,       "\0");
    PARAM  (93, "Xxxxxxxxxx",     DoBreaksParm,      "\0");
    PARAM  (94, "Nowarnings",     DoSetParm,         "Yes|No|\0");
    PARAM  (95, "Showtreeprobs",  DoSumtParm,        "Yes|No|\0");
    PARAM  (96, "Filename",       DoSumpParm,        "\0");
    PARAM  (97, "Burnin",         DoSumpParm,        "\0");
    PARAM  (98, "Reweight",       DoMcmcParm,        "\0");
    PARAM  (99, "Noop",           DoMcmcParm,        "\0");
    PARAM (100, "Ny98omega1pr",   DoPrsetParm,       "Beta|Fixed|\0");
    PARAM (101, "Ny98omega3pr",   DoPrsetParm,       "Uniform|Exponential|Fixed|\0");
    PARAM (102, "Codoncatfreqs",  DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM (103, "Sampleprob",     DoPrsetParm,       "\0");
    PARAM (104, "Aamodelpr",      DoPrsetParm,       "Fixed|Mixed|\0");
    PARAM (105, "Aamodel",        DoLinkParm,        "\0");
    PARAM (106, "Filename",       DoPlotParm,        "\0");
    PARAM (107, "Parameter",      DoPlotParm,        "\0");
    PARAM (108, "Match",          DoPlotParm,        "Perfect|Consistentwith|All|\0");
    PARAM (109, "Burnin",         DoPlotParm,        "\0");
    PARAM (110, "Brownscalepr",   DoPrsetParm,       "Uniform|Gamma|Gammamean|Fixed|\0");
    PARAM (111, "Browncorrpr",    DoPrsetParm,       "Uniform|Fixed|\0");
    PARAM (112, "Pbf",            DoMcmcParm,        "Yes|No|\0");
    PARAM (113, "Pbfinitburnin",  DoMcmcParm,        "\0");
    PARAM (114, "Pbfsamplefreq",  DoMcmcParm,        "\0");
    PARAM (115, "Pbfsampletime",  DoMcmcParm,        "\0");
    PARAM (116, "Pbfsampleburnin",DoMcmcParm,        "\0");
    PARAM (117, "Growthpr",       DoPrsetParm,       "Uniform|Exponential|Fixed|Normal|\0");
    PARAM (118, "Growthrate",     DoLinkParm,        "\0");
    PARAM (119, "Xxxxxxxxxx",     DoCalibrateParm,   "Unconstrained|Fixed|Uniform|Offsetexponential|Truncatednormal|Lognormal|Offsetlognormal|Gamma|Offsetgamma|\0");
    PARAM (120, "Calwaitpr",      DoPrsetParm,       "Exponential|Fixed|\0");     /* not used but leave it in to not destroy mapping to commands */
    PARAM (121, "M3omegapr",      DoPrsetParm,       "Exponential|Fixed|\0");
    PARAM (122, "Applyto",        DoReportParm,      "\0");
    PARAM (123, "Tratio",         DoReportParm,      "Dirichlet|Ratio|\0");
    PARAM (124, "Revmat",         DoReportParm,      "Dirichlet|Ratio|\0");
    PARAM (125, "Ratemult",       DoReportParm,      "Dirichlet|Scaled|Ratio|\0");
    PARAM (126, "Filename",       DoManualParm,      "\0");
    PARAM (127, "Filename1",      DoCompareTreeParm, "\0");
    PARAM (128, "Filename2",      DoCompareTreeParm, "\0");
    PARAM (129, "Outputname",     DoCompareTreeParm, "\0");
    PARAM (130, "Burnin",         DoCompareTreeParm, "\0");
    PARAM (131, "Ploidy",         DoLsetParm,        "Haploid|Diploid|Zlinked|\0");
    PARAM (132, "Swapadjacent",   DoMcmcParm,        "Yes|No|\0");
    PARAM (133, "Treeagepr",      DoPrsetParm,       "Fixed|Uniform|Offsetexponential|Truncatednormal|Lognormal|Offsetlognormal|Gamma|Offsetgamma|\0");
    PARAM (134, "Ancstates",      DoReportParm,      "Yes|No|\0");
    PARAM (135, "Siterates",      DoReportParm,      "Yes|No|\0");
    PARAM (136, "Possel",         DoReportParm,      "Yes|No|\0");
    PARAM (137, "Plot",           DoSumpParm,        "Yes|No|\0");
    PARAM (138, "Table",          DoSumpParm,        "Yes|No|\0");
    PARAM (139, "Minprob",        DoSumpParm,        "\0");
    PARAM (140, "Printtofile",    DoSumpParm,        "Yes|No|\0");
    PARAM (141, "Outputname",     DoSumpParm,        "\0");
    PARAM (142, "Redirect",       DoMcmcParm,        "Yes|No|\0");
    PARAM (143, "Swapseed",       DoMcmcParm,        "\0");
    PARAM (144, "Runidseed",      DoMcmcParm,        "\0");
    PARAM (145, "Quitonerror",    DoSetParm,         "Yes|No|\0");
    PARAM (146, "Savebrparams",   DoSumtParm,        "Yes|No|\0");
    PARAM (147, "Minbrparamfreq", DoSumtParm,        "\0");
    PARAM (148, "Minpartfreq",    DoMcmcParm,        "\0");
    PARAM (149, "Allchains",      DoMcmcParm,        "Yes|No|\0");
    PARAM (150, "Mcmcdiagn",      DoMcmcParm,        "Yes|No|\0");
    PARAM (151, "Diagnfreq",      DoMcmcParm,        "\0");
    PARAM (152, "Nruns",          DoMcmcParm,        "\0");
    PARAM (153, "Stoprule",       DoMcmcParm,        "Yes|No|\0");
    PARAM (154, "Stopval",        DoMcmcParm,        "\0");
    PARAM (155, "Relburnin",      DoMcmcParm,        "Yes|No|\0");
    PARAM (156, "Burninfrac",     DoMcmcParm,        "\0");
    PARAM (157, "Allcomps",       DoMcmcParm,        "Yes|No|\0");
    PARAM (158, "Printall",       DoMcmcParm,        "Yes|No|\0");
    PARAM (159, "Printmax",       DoMcmcParm,        "\0");
    PARAM (160, "Data",           DoMcmcParm,        "Yes|No|\0");
    PARAM (161, "Nruns",          DoSumpParm,        "\0");
    PARAM (162, "Allruns",        DoSumpParm,        "Yes|No|\0");
    PARAM (163, "Nruns",          DoSumtParm,        "\0");
    PARAM (164, "Ntrees",         DoSumtParm,        "\0");
    PARAM (165, "Calctreeprobs",  DoSumtParm,        "Yes|No|\0");
    PARAM (166, "Ordertaxa",      DoMcmcParm,        "Yes|No|\0");
    PARAM (167, "Ordertaxa",      DoSumtParm,        "Yes|No|\0");
    PARAM (168, "Aarevmatpr",     DoPrsetParm,       "Dirichlet|Fixed|\0");
    PARAM (169, "Nswaps",         DoMcmcParm,        "\0");
    PARAM (170, "Autoreplace",    DoSetParm,         "Yes|No|\0");
    PARAM (171, "Npthreads",      DoSetParm,         "\0");
    PARAM (172, "Cppratepr",      DoPrsetParm,       "Fixed|Exponential|\0");
    PARAM (173, "Cppmultdevpr",   DoPrsetParm,       "Fixed|\0");
    PARAM (174, "TK02varpr",      DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (175, "Pfile",          DoSumtParm,        "\0");
    PARAM (176, "Pfile",          DoSumtParm,        "\0");
    PARAM (177, "Autocomplete",   DoSumtParm,        "Yes|No|\0");
    PARAM (178, "Autocomplete",   DoSumpParm,        "Yes|No|\0");
    PARAM (179, "Userlevel",      DoSetParm,         "Standard|Developer|\0");
    PARAM (180, "Allavailable",   DoShowmovesParm,   "Yes|No|\0");
    PARAM (181, "Seed",           DoSetParm,         "\0");
    PARAM (182, "Swapseed",       DoSetParm,         "\0");
    PARAM (183, "Clockratepr",    DoPrsetParm,       "Fixed|Normal|Lognormal|Exponential|Gamma|\0");
    PARAM (184, "Nodeagepr",      DoPrsetParm,       "Unconstrained|Calibrated|\0");
    PARAM (185, "Clockvarpr",     DoPrsetParm,       "Strict|Cpp|TK02|Igr|Bm|Ibr|Mixed|\0");
    PARAM (186, "Xxxxxxxxxx",     DoPropsetParm,     "\0");
    PARAM (187, "Xxxxxxxxxx",     DoStartvalsParm,   "\0");
    PARAM (188, "Usegibbs",       DoLsetParm,        "Yes|No|\0");
    PARAM (189, "Gibbsfreq",      DoLsetParm,        "\0");
    PARAM (190, "Checkpoint",     DoMcmcParm,        "Yes|No|\0");
    PARAM (191, "Checkfreq",      DoMcmcParm,        "\0");
    PARAM (192, "Tree",           DoReportParm,      "Topology|Brlens|\0");
    PARAM (193, "Cpprate",        DoLinkParm,        "\0");
    PARAM (194, "Cppmultdev",     DoLinkParm,        "\0");
    PARAM (195, "Cppevents",      DoLinkParm,        "\0");
    PARAM (196, "TK02var",        DoLinkParm,        "\0");
    PARAM (197, "TK02branchrates",DoLinkParm,        "\0");
    PARAM (198, "Savetrees",      DoMcmcParm,        "Yes|No|\0");
    PARAM (199, "Diagnstat",      DoMcmcParm,        "Avgstddev|Maxstddev|\0");
    PARAM (200, "Startparams",    DoMcmcParm,        "Reset|Current|\0");
    PARAM (201, "Characters",     DoBeginParm,       "\0");
    PARAM (202, "Startingtrees",  DoMcmcParm,        "\0");
    PARAM (203, "Xxxxxxxxxx",     DoUserTreeParm,    "\0");
    PARAM (204, "Outputname",     DoSumtParm,        "\0");
    PARAM (205, "Table",          DoSumtParm,        "Yes|No|\0");
    PARAM (206, "Summary",        DoSumtParm,        "Yes|No|\0");
    PARAM (207, "Consensus",      DoSumtParm,        "Yes|No|\0");
    PARAM (208, "Minpartfreq",    DoSumtParm,        "\0");
    PARAM (209, "Relburnin",      DoSumtParm,        "Yes|No|\0");
    PARAM (210, "Burninfrac",     DoSumtParm,        "\0");
    PARAM (211, "Relburnin",      DoSumpParm,        "Yes|No|\0");
    PARAM (212, "Burninfrac",     DoSumpParm,        "\0");
    PARAM (213, "Append",         DoMcmcParm,        "Yes|No|\0");
    PARAM (214, "Autotune",       DoMcmcParm,        "Yes|No|\0");
    PARAM (215, "Tunefreq",       DoMcmcParm,        "\0");
    PARAM (216, "Scientific",     DoSetParm,         "Yes|No|\0");
    PARAM (217, "Siteomega",      DoReportParm,      "Yes|No|\0");
    PARAM (218, "Igrvarpr",       DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (219, "Symbols",        DoFormatParm,      "\0");
    PARAM (220, "Equate",         DoFormatParm,      "\0");
    PARAM (221, "Relburnin",      DoCompareTreeParm, "Yes|No|\0");
    PARAM (222, "Burninfrac",     DoCompareTreeParm, "\0");
    PARAM (223, "Minpartfreq",    DoCompareTreeParm, "\0");
    PARAM (224, "Relburnin",      DoPlotParm,        "Yes|No|\0");
    PARAM (225, "Burninfrac",     DoPlotParm,        "\0");
    PARAM (226, "Taxa",           DoBeginParm,       "\0");
    PARAM (227, "Xxxxxxxxxx",     DoBeginParm,       "\0");
    PARAM (228, "Xxxxxxxxxx",     DoTaxlabelsParm,   "\0");
    PARAM (229, "Dir",            DoSetParm,         "\0");
    PARAM (230, "Conformat",      DoSumtParm,        "Figtree|Simple|\0");
    PARAM (231, "Hpd",            DoSumpParm,        "Yes|No|\0");
    PARAM (232, "Hpd",            DoSumtParm,        "Yes|No|\0");
    PARAM (233, "Usebeagle",      DoSetParm,         "Yes|No|\0");
    PARAM (234, "Beagledevice",   DoSetParm,         "Cpu|Gpu|\0");
    PARAM (235, "Beagleprecision",DoSetParm,         "Single|Double|\0");
    PARAM (236, "Beaglesse",      DoSetParm,         "Yes|No|\0");
    PARAM (237, "Beagleopenmp",   DoSetParm,         "Yes|No|\0");
    PARAM (238, "Beaglethreads",  DoSetParm,         "Yes|No|\0");
    PARAM (239, "Beaglescaling",  DoSetParm,         "Always|Dynamic|\0");
    PARAM (240, "Beaglefreq",     DoSetParm,         "\0");
    PARAM (241, "Popvarpr",       DoPrsetParm,       "Equal|Variable|\0");
    PARAM (242, "Igrvar",         DoLinkParm,        "\0");
    PARAM (243, "Igrbranchrates", DoLinkParm,        "\0");
    PARAM (244, "Xxxxxxxxxx",     DoSpeciespartitionParm,   "\0");
    PARAM (245, "Speciespartition",  DoSetParm,      "\0");
    PARAM (246, "Revratepr",      DoPrsetParm,       "Symdir|\0");
    PARAM (247, "Samplestrat",    DoPrsetParm,       "Random|Diversity|Cluster|FossilTip|\0");
    PARAM (248, "Burninss",       DoSsParm,          "\0");
    PARAM (249, "Nsteps",         DoSsParm,          "\0");
    PARAM (250, "Alpha",          DoSsParm,          "\0");
    PARAM (251, "Bmvarpr",        DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (252, "Bmvar",          DoLinkParm,        "\0");
    PARAM (253, "Bmbranchrates",  DoLinkParm,        "\0");
    PARAM (254, "Ibrvarpr",       DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (255, "Ibrvar",         DoLinkParm,        "\0");
    PARAM (256, "Ibrbranchlens",  DoLinkParm,        "\0");
    PARAM (257, "FromPrior",      DoSsParm,          "Yes|No|\0");
    PARAM (258, "Filename",       DoSumSsParm,       "\0");
    PARAM (259, "Burnin",         DoSumSsParm,       "\0");
    PARAM (260, "Nruns",          DoSumSsParm,       "\0");
    PARAM (261, "Allruns",        DoSumSsParm,       "Yes|No|\0");
    PARAM (262, "Askmore",        DoSumSsParm,       "Yes|No|\0");
    PARAM (263, "Relburnin",      DoSumSsParm,       "Yes|No|\0");
    PARAM (264, "Burninfrac",     DoSumSsParm,       "\0");
    PARAM (265, "Discardfrac",    DoSumSsParm,       "\0");
    PARAM (266, "Smoothing",      DoSumSsParm,       "\0");
    PARAM (267, "Steptoplot",     DoSumSsParm,       "\0");
    PARAM (268, "Precision",      DoSetParm,         "\0");
    PARAM (269, "Fossilizationpr",   DoPrsetParm,    "Beta|Fixed|\0");
    PARAM (270, "Fossilizationrate", DoLinkParm,     "\0");
    PARAM (271, "Generatepr",     DoPrsetParm,       "Variable|Fixed|\0");
    PARAM (272, "Mixedvarpr",     DoPrsetParm,       "Fixed|Exponential|Uniform|\0");
    PARAM (273, "Mixedvar",       DoLinkParm,        "\0");
    PARAM (274, "Mixedbrchrates", DoLinkParm,        "\0");
    PARAM (275, "Beagleresource", DoSetParm,         "\0");
    PARAM (276, "Nlnormcat",      DoLsetParm,        "\0");
    PARAM (277, "Nmixtcat",       DoLsetParm,        "\0");

    /* NOTE: If a change is made to the parameter table, make certain you change
            NUMPARAMS (now 278; one more than last index) at the top of this file. */
    /* CmdType commands[] */
}


void ShowNodes (TreeNode *p, int indent, int isThisTreeRooted)
{
    if (p != NULL)
        {
        printf ("   ");
        if (p->left == NULL && p->right == NULL && p->anc != NULL)
            {
            printf ("%*cN %d (l=%d r=%d a=%d) %1.15lf (%s) isDated=%d ",
            indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->length, p->label, p->isDated);
            }
        else if (p->left != NULL && p->right == NULL && p->anc == NULL)
            {
            if (isThisTreeRooted == NO)
                {
                if (p->label[0] == '\0' || p->label[0] == '\n' || p->label[0] == ' ')
                    printf ("%*cN %d (l=%d r=%d a=%d) (---) ",
                    indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
                else
                    printf ("%*cN %d (l=%d r=%d a=%d) (%s) ",
                    indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->label);
                }
            else
                {
                printf ("%*cN %d (l=%d r=%d a=%d) X.XXXXXX ",
                indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
                }
            }
        else
            {
            if (p->anc != NULL)
                {
                if (p->anc->anc == NULL && isThisTreeRooted == YES)
                    printf ("%*cN %d (l=%d r=%d a=%d) X.XXXXXX ",
                    indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc));
                else    
                    printf ("%*cN %d (l=%d r=%d a=%d) %1.15lf ",
                    indent, ' ', Dex(p), Dex(p->left), Dex(p->right), Dex(p->anc), p->length);
                }
            }
        if (isThisTreeRooted == YES)
            printf ("depth=%1.15lf\n", p->nodeDepth);
        else
            printf ("\n");
        ShowNodes (p->left,  indent + 2, isThisTreeRooted);
        ShowNodes (p->right, indent + 2, isThisTreeRooted);
        }
}


int StandID (char nuc)
{
    char        n;
    
    /* Note that if you change how many states are recognized, you need 
       to look at IsMissing */
    n = nuc;

    if (n == '0')
        {
        return 1;
        }
    else if (n == '1')
        {
        return 2;
        }
    else if (n == '2')
        {
        return 4;
        }
    else if (n == '3')
        {
        return 8;
        }
    else if (n == '4')
        {
        return 16;
        }
    else if (n == '5')
        {
        return 32;
        }
    else if (n == '6')
        {
        return 64;
        }
    else if (n == '7')
        {
        return 128;
        }
    else if (n == '8')
        {
        return 256;
        }
    else if (n == '9')
        {
        return 512;
        }
    else if (n == missingId)
        {
        return MISSING;
        }
    else if (n == gapId)
        {
        return GAP;
        }
    else
        return -1;
}


void State_CODON (char *state, int code, int division)
{
    state[0] = StateCode_NUC4(modelParams[division].codonNucs[code][0]);
    state[1] = StateCode_NUC4(modelParams[division].codonNucs[code][1]);
    state[2] = StateCode_NUC4(modelParams[division].codonNucs[code][2]);
    state[3] = '\0';
}


void State_DOUBLET (char *state, int code)
{
    state[0] = code/4 + 'A';
    state[1] = code%4 + 'A';
    state[2] = '\0';
}


int StateCode_AA (int n)
{
    if (n == 0)
        return 'A';      /* Ala */
    else if (n == 1)
        return 'R';      /* Arg */
    else if (n == 2)
        return 'N';      /* Asn */
    else if (n == 3)
        return 'D';      /* Asp */
    else if (n == 4)
        return 'C';      /* Cys */
    else if (n == 5)
        return 'Q';      /* Gln */
    else if (n == 6)
        return 'E';      /* Glu */
    else if (n == 7)
        return 'G';      /* Gly */
    else if (n == 8)
        return 'H';      /* His */
    else if (n == 9)
        return 'I';      /* Ile */
    else if (n == 10)
        return 'L';      /* Leu */
    else if (n == 11)
        return 'K';      /* Lys */
    else if (n == 12)
        return 'M';      /* Met */
    else if (n == 13)
        return 'F';      /* Phe */
    else if (n == 14)
        return 'P';      /* Pro */
    else if (n == 15)
        return 'S';      /* Ser */
    else if (n == 16)
        return 'T';      /* Thr */
    else if (n == 17)
        return 'W';      /* Trp */
    else if (n == 18)
        return 'Y';      /* Tyr */
    else if (n == 19)
        return 'V';      /* Val */
    else
        return '?';
}


int StateCode_NUC4 (int n)
{
    if (n == 0)
        return 'A';
    else if (n == 1)
        return 'C';
    else if (n == 2)
        return 'G';
    else if (n == 3)
        return 'T';
    else return '?';
}


int StateCode_Std (int n)
{
    if (n <= 9 && n >= 0)
        return '0' + n;
    else return '?';
}


void WhatVariableExp (BitsLong exp, char *st)
{
    int         n;
    
    strcpy (st, "");
    n = 0;
    if (exp == 0)
        strcat(st, " nothing");
    else
        {
        if ((exp & Expecting(COMMAND)) == Expecting(COMMAND))
            {
            strcat(st, " command");
            n++;
            }
        if ((exp & Expecting(PARAMETER)) == Expecting(PARAMETER))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " parameter");
            n++;
            }
        if ((exp & Expecting(EQUALSIGN)) == Expecting(EQUALSIGN))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " =");
            n++;
            }
        if ((exp & Expecting(COLON)) == Expecting(COLON))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " :");
            n++;
            }
        if ((exp & Expecting(SEMICOLON)) == Expecting(SEMICOLON))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " ;");
            n++;
            }
        if ((exp & Expecting(COMMA)) == Expecting(COMMA))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " ,");
            n++;
            }
        if ((exp & Expecting(POUNDSIGN)) == Expecting(POUNDSIGN))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " #");
            n++;
            }
        if ((exp & Expecting(QUESTIONMARK)) == Expecting(QUESTIONMARK))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " ?");
            n++;
            }
        if ((exp & Expecting(DASH)) == Expecting(DASH))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " -");
            n++;
            }
        if ((exp & Expecting(LEFTPAR)) == Expecting(LEFTPAR))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " (");
            n++;
            }
        if ((exp & Expecting(RIGHTPAR)) == Expecting(RIGHTPAR))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " )");
            n++;
            }
        if ((exp & Expecting(LEFTCOMMENT)) == Expecting(LEFTCOMMENT))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " [");
            n++;
            }
        if ((exp & Expecting(RIGHTCOMMENT)) == Expecting(RIGHTCOMMENT))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " ]");
            n++;
            }
        if ((exp & Expecting(ALPHA)) == Expecting(ALPHA))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " <name>");
            n++;
            }
        if ((exp & Expecting(NUMBER)) == Expecting(NUMBER))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " <number>");
            n++;
            }
        if ((exp & Expecting(RETURNSYMBOL)) == Expecting(RETURNSYMBOL))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " return");
            n++;
            }
        if ((exp & Expecting(ASTERISK)) == Expecting(ASTERISK))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " *");
            n++;
            }
        if ((exp & Expecting(BACKSLASH)) == Expecting(BACKSLASH))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " /");
            n++;
            }
        if ((exp & Expecting(BACKSLASH)) == Expecting(BACKSLASH))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " \\");
            n++;
            }
        if ((exp & Expecting(EXCLAMATIONMARK)) == Expecting(EXCLAMATIONMARK))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " !");
            n++;
            }
        if ((exp & Expecting(PERCENT)) == Expecting(PERCENT))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " %");
            n++;
            }
        if ((exp & Expecting(LEFTCURL)) == Expecting(LEFTCURL))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " {");
            n++;
            }
        if ((exp & Expecting(RIGHTCURL)) == Expecting(RIGHTCURL))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " }");
            n++;
            }
        if ((exp & Expecting(WEIRD)) == Expecting(WEIRD))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " <whatever>");
            n++;
            }
        if ((exp & Expecting(VERTICALBAR)) == Expecting(VERTICALBAR))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " |");
            n++;
            }
        if ((exp & Expecting(UNKNOWN_TOKEN_TYPE)) == Expecting(UNKNOWN_TOKEN_TYPE))
            {
            if (n > 0)
                strcat(st, " or");
            strcat(st, " no clue");
            n++;
            }
        }
}


char WhichAA (int x)
{
    if (x == 1)
        return ('A');
    else if (x == 2)
        return ('R');
    else if (x == 4)
        return ('N');
    else if (x == 8)
        return ('D');
    else if (x == 16)
        return ('C');
    else if (x == 32)
        return ('Q');
    else if (x == 64)
        return ('E');
    else if (x == 128)
        return ('G');
    else if (x == 256)
        return ('H');
    else if (x == 512)
        return ('I');
    else if (x == 1024)
        return ('L');
    else if (x == 2048)
        return ('K');
    else if (x == 4096)
        return ('M');
    else if (x == 8192)
        return ('F');
    else if (x == 16384)
        return ('P');
    else if (x == 32768)
        return ('S');
    else if (x == 65536)
        return ('T');
    else if (x == 131072)
        return ('W');
    else if (x == 262144)
        return ('Y');
    else if (x == 524288)
        return ('V');
    else if (x > 0 && x < 524288)
        return ('*');
    else if (x == MISSING)
        return ('?');
    else if (x == GAP)
        return ('-');
    else 
        return (' ');
}


MrBFlt WhichCont (int x)
{
    return ((MrBFlt)(x / 1000.0));
}


char WhichNuc (int x)
{
    if (x == 1)
        return ('A');
    else if (x == 2)
        return ('C');
    else if (x == 3)
        return ('M');
    else if (x == 4)
        return ('G');
    else if (x == 5)
        return ('R');
    else if (x == 6)
        return ('S');
    else if (x == 7)
        return ('V');
    else if (x == 8)
        return ('T');
    else if (x == 9)
        return ('W');
    else if (x == 10)
        return ('Y');
    else if (x == 11)
        return ('H');
    else if (x == 12)
        return ('K');
    else if (x == 13)
        return ('D');
    else if (x == 14)
        return ('B');
    else if (x == 15)
        return ('N');
    else if (x == MISSING)
        return ('?');
    else if (x == GAP)
        return ('-');
    else 
        return (' ');
}


char WhichRes (int x)
{
    if (x == 1)
        return ('0');
    else if (x == 2)
        return ('1');
    else if (x == 3)
        return ('*');
    else if (x == MISSING)
        return ('N');
    else if (x == GAP)
        return ('-');
    else 
        return (' ');
}


char WhichStand (int x)
{
    if (x == 1)
        return ('0');
    else if (x == 2)
        return ('1');
    else if (x == 4)
        return ('2');
    else if (x == 8)
        return ('3');
    else if (x == 16)
        return ('4');
    else if (x == 32)
        return ('5');
    else if (x == 64)
        return ('6');
    else if (x == 128)
        return ('7');
    else if (x == 256)
        return ('8');
    else if (x == 512)
        return ('9');
    else if (x > 0 && x < 512)
        return ('*');
    else if (x == MISSING)
        return ('N');
    else if (x == GAP)
        return ('-');
    else 
        return (' ');
}

