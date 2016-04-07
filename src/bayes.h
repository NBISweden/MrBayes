#ifndef __BAYES_H__
#define __BAYES_H__

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>

#ifdef HAVE_CONFIG_H
#   include "config.h"
#define VERSION_NUMBER  PACKAGE_VERSION
#elif !defined (XCODE_VERSION) /* some defaults that would otherwise be guessed by configure */
#   define PACKAGE_NAME "MrBayes"
#   define PACKAGE_VERSION "3.2.7-dev"
#   define VERSION_NUMBER  PACKAGE_VERSION
#   undef  HAVE_LIBREADLINE
#   define UNIX_VERSION 1
#   define SSE_ENABLED  1
#   undef  AVX_ENABLED
#   undef  FMA_ENABLED
#   undef  MPI_ENABLED
#   undef  BEAGLE_ENABLED
#   undef  FAST_LOG
#endif

/* Set SSE_ENABLED if SSE SIMD extensions available. */
#ifdef HAVE_SSE
#define SSE_ENABLED
#endif

/* Set AVX_ENABLED if AVX SIMD extensions available. */
#ifdef HAVE_AVX
#define AVX_ENABLED
#endif

/* Set FMA_ENABLED if FMA SIMD extensions available. */
#if defined(HAVE_FMA3) || defined(HAVE_FMA4)
#define FMA_ENABLED
#endif

/* Set COMPLETIONMATCHES if we have tho readline library */
#ifdef HAVE_LIBREADLINE
#define COMPLETIONMATCHES
#endif

#if defined (MPI_ENABLED)
#include "mpi.h"
#endif

#if defined (BEAGLE_ENABLED)
#include "libhmsbeagle/beagle.h"
#endif

#if !defined (UNIX_VERSION) && !defined (WIN_VERSION) && !defined (MAC_VERSION)
#  ifdef __MWERKS__
#    define MAC_VERSION
#  elif defined __APPLE__
#    define MAC_VERSION
#  elif defined __unix__
#    define UNIX_VERISON
#  else
#    define WIN_VERSION
#  endif
#endif

#if defined (WIN_VERSION)
#   include <windows.h>
#   include <winbase.h>
#endif

/* Previous problems with bitfield operations may have been caused by several things. One potential
   problem has been that MrBayes has used signed ints or signed longs for the bit operations, and
   not all compilers perform bitfield operations on signed ints and longs in the expected manner
   for unsigned ints and longs. Actually, the ANSI standard does not specify how to implement
   bit operations on signed operands.
   
   Another problem is that some bitshift operations in MrBayes have been performed on an integer
   constant "1". This is probably interpreted by most compilers as a signed int and the result of
   the bitshift operation is therefore also a signed int. When a signed int has had a different
   number of bits than the type specified as SafeLong, this may also have resulted in errors on
   some systems.
 
   Both of these problems are fixed now by using unsigned longs for the bitfield operations and
   by setting the "1" value needed by some bitfield functions using a variable defined as the same
   type as the bitfield containers themselves. Furthermore, I have separated out the cases where
   a signed long is required or is more appropriate than the container type of bitfields.

   FR 2013-07-06
 */
typedef unsigned long BitsLong;
typedef long RandLong;

#define MRBFLT_MAX DBL_MAX  /* maximum possible value that can be stored in MrBFlt */
#define MRBFLT_MIN DBL_MIN  /* minimum possible value that can be stored in MrBFlt */
#define MRBFLT_NEG_MAX (-DBL_MAX)  /* maximum possible negative value that can be stored in MrBFlt */
typedef double MrBFlt;      /* double used for parameter values and generally for floating point values, 
                               if set to float MPI would not work becouse of use MPI_DOUBLE */
typedef float CLFlt;        /* single-precision float used for cond likes (CLFlt) to increase speed and reduce memory requirement */
                            /* set CLFlt to double if you want increased precision */
                            /* NOTE: CLFlt = double not compatible with SSE_ENABLED, AVX_ENABLED or FMA_ENABLED */

/*
 * Make sure we define SIMD instruction flags in a stepwise manner. That is, if we have FMA, make sure we have AVX;
 * if we have AVX, make sure we have SSE.
 */
#if defined (FMA_ENABLED)
#  if !defined (AVX_ENABLED)
#    define AVX_ENABLED
#  endif
#endif

#if defined (AVX_ENABLED)
#  if !defined (SSE_ENABLED)
#    define SSE_ENABLED
#  endif
#endif

/* Define compiler for appropriate SIMD vector data alignment and free operations */
#if defined (SSE_ENABLED)
#  if defined (WIN_VERSION)
#    define MS_VCPP_SIMD
#  else
#    define GCC_SIMD
#    undef ICC_SIMD
#  endif
#endif

/* Include header for SSE code */
#if defined (SSE_ENABLED)
#  include <xmmintrin.h>
#endif

/* Include header for AVX code */
#if defined (AVX_ENABLED)
#  include <immintrin.h>
#endif

#if defined (FMA_ENABLED)
#  include <immintrin.h>
#endif

/* Define vector code constants */
#define VEC_NONE    0
#define VEC_SSE     1
#define VEC_AVX     2
#define VEC_FMA     3

/* For comparing floating points: two values are the same if the absolute difference is less then
   this value.
*/
#ifndef ETA
#define ETA (1E-30)
#endif

/* This defined DEBUG() is not used anywhere
#if defined (DEBUGOUTPUT)
#define DEBUG(fmt, arg) printf("%s:%d ",__FILE__,__LINE__); printf(fmt,arg);
#endif
*/

/* TEMPSTRSIZE determines size of temporary sprintf buffer (for SafeSprintf) */
/* A patch was sent in by Allen Smith for SafeSprintf, but I could not get
   it compiled on SGI IRIX 6.5 (too old?) with _xpg5_vsnprintf undefined.
   The code below is a hack so SafeSprintf never has to reallocate memory.
*/
#ifdef __sgi
#define TEMPSTRSIZE 1000
#else
#define TEMPSTRSIZE 200
#endif

/* NO_ERROR is defined in bayes.h (as 0) and also in WinError.h (as 0L)
 ERROR is defined in bayes.h (as 1) and also in WinGDI.h (as 0). we use the bayes.h value */
#undef  NO_ERROR
#undef  ERROR
#define NO_ERROR                0
#define ERROR                   1
#define NO_ERROR_QUIT           2
#define ABORT                   3
#define SKIP_COMMAND            4

#undef  FALSE
#undef  TRUE
#define FALSE                   0
#define TRUE                    1

#define NO                      0
#define YES                     1

#define UP                      0
#define DOWN                    1

#define UPPER                   0
#define MIDDLE                  1
#define LOWER                   2

#define NONINTERACTIVE          0
#define INTERACTIVE             1

#define STANDARD_USER           1
#define DEVELOPER               3

#define DIFFERENT               0
#define SAME                    1
#define CONSISTENT_WITH         2

#define LINETERM_UNIX           0
#define LINETERM_MAC            1
#define LINETERM_DOS            2

#define SCREENWIDTH             60
#define SCREENWIDTH2            61

#define AVGSTDDEV               0
#define MAXSTDDEV               1

#define NONE                    0
#define DNA                     1
#define RNA                     2
#define PROTEIN                 3
#define RESTRICTION             4
#define STANDARD                5
#define MIXED                   6
#define CONTINUOUS              7

#define AAMODEL_POISSON         0
#define AAMODEL_JONES           1
#define AAMODEL_DAY             2
#define AAMODEL_MTREV           3
#define AAMODEL_MTMAM           4
#define AAMODEL_WAG             5
#define AAMODEL_RTREV           6
#define AAMODEL_CPREV           7
#define AAMODEL_VT              8
#define AAMODEL_BLOSUM          9
#define AAMODEL_LG              10
#define AAMODEL_EQ              11
#define AAMODEL_GTR             12 /* aa models with free parameters must be listed last */

#define NUCMODEL_4BY4           0
#define NUCMODEL_DOUBLET        1
#define NUCMODEL_CODON          2
#define NUCMODEL_AA             3

#define NST_MIXED              -1  /* anything other than 1, 2, or 6 */

#define MISSING                 10000000
#define GAP                     10000001

#define UNORD                   0
#define ORD                     1
#define DOLLO                   2
#define IRREV                   3

#define IN_CMD                  0
#define IN_FILE                 1

#define NOTHING                 0
#define COMMAND                 1
#define PARAMETER               2
#define EQUALSIGN               3
#define COLON                   4
#define SEMICOLON               5
#define COMMA                   6
#define POUNDSIGN               7
#define QUESTIONMARK            8
#define DASH                    9
#define LEFTPAR                 10
#define RIGHTPAR                11
#define LEFTCOMMENT             12
#define RIGHTCOMMENT            13
#define ALPHA                   14
#define NUMBER                  15
#define RETURNSYMBOL            16
#define ASTERISK                17
#define BACKSLASH               18
#define FORWARDSLASH            19
#define EXCLAMATIONMARK         20
#define PERCENT                 21
#define QUOTATIONMARK           22
#define WEIRD                   23
#define UNKNOWN_TOKEN_TYPE      24
#define LEFTCURL                25
#define RIGHTCURL               26
#define DOLLAR                  27
#define AMPERSAND               28
#define VERTICALBAR             29

#define MAX_Q_RATE              100.0f
#define MIN_SHAPE_PARAM         0.00001f
#define MAX_SHAPE_PARAM         100.0f
#define MAX_SITE_RATE           10.0f
#define MAX_GAMMA_CATS          20
#define MAX_GAMMA_CATS_SQUARED  400
#define BRLENS_MIN              0.00000001f  // 1E-8f
#define BRLENS_MAX              100.0f
/* BRLENS_MIN must be bigger than TIME_MIN */
#define TIME_MIN                1.0E-11f
#define TIME_MAX                100.0f
#define RELBRLENS_MIN           0.00000001f  // 1E-8f
#define RELBRLENS_MAX           100.0f
#define KAPPA_MIN               0.001f
#define KAPPA_MAX               1000.0f
#define GROWTH_MIN              0.000001f
#define GROWTH_MAX              1000000.0f
#define RATE_MIN                0.000001f
#define RATE_MAX                100.0f
#define CPPRATEMULTIPLIER_MIN   0.001f
#define CPPRATEMULTIPLIER_MAX   1000.0f
#define SYMPI_MIN               0.000001f
#define SYMPI_MAX               100.0f
#define ALPHA_MIN               0.0001f
#define ALPHA_MAX               10000.0f
#define DIR_MIN                 0.000001f
#define PI_MIN                  0.000001f
#define OFFSETEXPLAMBDA_MIN     0.000001f
#define OFFSETEXPLAMBDA_MAX     100000.0f
#define TREEHEIGHT_MIN          0.00000001f
#define TREEHEIGHT_MAX          1000.0f
#define TREEAGE_MIN             0.00000001f
#define TREEAGE_MAX             1000000.0f
#define CPPLAMBDA_MIN           0.00001f
#define CPPLAMBDA_MAX           100.0f
#define TK02VAR_MIN             0.000001f
#define TK02VAR_MAX             10000.0f
#define IGRVAR_MIN              0.000001f
#define IGRVAR_MAX              10000.0f
#define MIXEDVAR_MIN            0.000001f
#define MIXEDVAR_MAX            10000.0f
#define OMEGA_MIN               0.001f
#define OMEGA_MAX               1000.0f

#define POS_MIN                 1E-25f
#define POS_MAX                 1E25f
#define POS_INFINITY            1E25f
#define NEG_INFINITY            -1000000.0f
#define POSREAL_MIN             1E-25f
#define POSREAL_MAX             1E25f

#define CMD_STRING_LENGTH       100000

#define pos(i,j,n)              ((i)*(n)+(j))

#define NUM_ALLOCS               100

#define ALLOC_MATRIX             0
#define ALLOC_CHARINFO           2
#define ALLOC_CHARSETS           3
#define ALLOC_TAXA               4
#define ALLOC_TMPSET             5
#define ALLOC_PARTITIONS         6
#define ALLOC_PARTITIONVARS      7
#define ALLOC_TAXASETS           8
#define ALLOC_CONSTRAINTS        9
#define ALLOC_USERTREE           10
#define ALLOC_SUMTPARAMS         11
#define ALLOC_TERMSTATE          12
#define ALLOC_ISPARTAMBIG        13
#define ALLOC_AVAILNODES         25
#define ALLOC_AVAILINDICES       26
#define ALLOC_CURLNL             28
#define ALLOC_CURLNPR            29
#define ALLOC_CHAINID            30
#define ALLOC_PARAMS             31
#define ALLOC_TREE               32
#define ALLOC_NODES              33
#define ALLOC_LOCTAXANAMES       34
#define ALLOC_COMPMATRIX         39
#define ALLOC_NUMSITESOFPAT      40
#define ALLOC_COMPCOLPOS         41
#define ALLOC_COMPCHARPOS        42
#define ALLOC_ORIGCHAR           43
#define ALLOC_PARAMVALUES        46
#define ALLOC_MCMCTREES          47
#define ALLOC_MOVES              48
#define ALLOC_PRELIKES           52
#define ALLOC_SITEJUMP           54
#define ALLOC_MARKOVTIS          55
#define ALLOC_RATEPROBS          56
#define ALLOC_STDTYPE            57
#define ALLOC_PACKEDTREES        58
#define ALLOC_SUMPSTRING         62
#define ALLOC_SUMPINFO           63
#define ALLOC_SWAPINFO           64
#define ALLOC_SYMPIINDEX         65
#define ALLOC_POSSELPROBS        66
#define ALLOC_PBF                68
#define ALLOC_LOCALTAXONCALIBRATION  69
#define ALLOC_SPR_PARSSETS       72
#define ALLOC_PFCOUNTERS         74
#define ALLOC_FILEPOINTERS       75
#define ALLOC_STATS              76
#define ALLOC_DIAGNTREE          77
#define ALLOC_USEDMOVES          82
#define ALLOC_MODEL              83
#define ALLOC_STDSTATEFREQS      84
#define ALLOC_PRINTPARAM         85
#define ALLOC_TREELIST           86
#define ALLOC_TFILEPOS           87
#define ALLOC_BEST               88
#define ALLOC_SPECIESPARTITIONS  89
#define ALLOC_SS                 90
#define ALLOC_SAMPLEFOSSILSLICE  91

#define LINKED                  0
#define UNLINKED                1

/*paramType*/
#define NUM_LINKED              31
#define P_TRATIO                0
#define P_REVMAT                1
#define P_OMEGA                 2
#define P_PI                    3
#define P_SHAPE                 4
#define P_PINVAR                5
#define P_CORREL                6
#define P_SWITCH                7
#define P_RATEMULT              8
#define P_TOPOLOGY              9
#define P_BRLENS                10
#define P_SPECRATE              11
#define P_EXTRATE               12
#define P_FOSLRATE              13
#define P_POPSIZE               14
#define P_AAMODEL               15
#define P_BRCORR                16
#define P_BRSIGMA               17
#define P_GROWTH                18
#define P_CPPMULTDEV            19
#define P_CPPRATE               20
#define P_CPPEVENTS             21
#define P_TK02VAR               22
#define P_TK02BRANCHRATES       23
#define P_IGRVAR                24
#define P_IGRBRANCHRATES        25
#define P_CLOCKRATE             26
#define P_SPECIESTREE           27
#define P_GENETREERATE          28
#define P_MIXEDVAR              29
#define P_MIXEDBRCHRATES        30
/* NOTE: If you add another parameter, change NUM_LINKED */

// #define CPPm                 0       /* CPP rate multipliers */
// #define CPPi                 1       /* CPP independent rates */
#define RCL_TK02                0
#define RCL_IGR                 1       /* type of mixed relaxed clock model */

#define MAX_NUM_USERTREES       200     /* maximum number of user trees MrBayes will read */
#define MAX_CHAINS              256     /* maximum numbder of chains you can run actually only half of it becouse of m->lnLike[MAX_CHAINS] */

// #define PARAM_NAME_SIZE      400

typedef void * VoidPtr;
typedef int (*CmdFxn)(void);
typedef int (*ParmFxn)(char *, char *);

/* typedef for a ln prior prob fxn */
typedef MrBFlt (*LnPriorProbFxn)(MrBFlt val, MrBFlt *priorParams);

/* typedef for a ln prior prob ratio fxn */
typedef MrBFlt (*LnPriorRatioFxn)(MrBFlt newVal, MrBFlt oldVal, MrBFlt *priorParams);

typedef struct
    {
    MrBFlt          sum;            /* sum of standard deviations */
    MrBFlt          max;            /* maximum standard deviation */
    MrBFlt          numPartitions;
    MrBFlt          numSamples;
    MrBFlt          avgStdDev;
    MrBFlt          **pair;
    } STATS;

/* enumeration for calibration prior type */
enum CALPRIOR
    {
    unconstrained,
    fixed,
    uniform,
    offsetExponential,
    truncatedNormal,
    logNormal,
    offsetLogNormal,
    standardGamma,
    offsetGamma 
    };

enum ConstraintType
    {
    PARTIAL,
    NEGATIVE,
    HARD
    };

enum CodingType
    {
    ALL                 = 0,
    NOABSENCESITES      = 1,
    NOPRESENCESITES     = 2,
    VARIABLE            = 3,
    NOSINGLETONPRESENCE = 4,
    NOSINGLETONABSENCE  = 8,
    NOSINGLETONS        = 12,
    INFORMATIVE         = 15
    };

/* typedef for calibration */
typedef struct calibration
    {
    char                name[100];
    enum CALPRIOR       prior;
    MrBFlt              priorParams[3];
    LnPriorProbFxn      LnPriorProb;
    LnPriorRatioFxn     LnPriorRatio;
    MrBFlt              min;
    MrBFlt              max;
    }
    Calibration;

/* typedef for tree (topology) list element */
typedef struct element
    {
    struct element *next;
    int             *order;
    } TreeListElement;

/* typedef for list of trees (topologies) */
typedef struct
    {
    TreeListElement *first;
    TreeListElement *last;
    } TreeList;

/* typedef for packed tree */
typedef struct
    {
    int     *order;
    MrBFlt  *brlens;
    } PackedTree;

/* typedef for binary tree node */
/* NOTE: Any variable added here must also be copied in CopyTrees */
typedef struct node
    {
    char            *label;                 /*!< name of node if tip                        */
    struct node     *left, *right, *anc;    /*!< pointers to adjacent nodes                 */
    int             memoryIndex;            /*!< memory index (do not change)               */
    int             index;                  /*!< index to node (0 to numLocalTaxa for tips) */
    int             upDateCl;               /*!< cond likes need update?                    */
    int             upDateTi;               /*!< transition probs need update?              */
    int             isLocked;               /*!< is node locked?                            */
    int             lockID;                 /*!< id of lock                                 */
    int             isDated;                /*!< is node dated (calibrated)?                */
    int             marked, x, y;           /*!< scratch variables                          */
    MrBFlt          d;                      /*!< scratch variable                           */
    BitsLong        *partition;             /*!< pointer to bitfield describing splits      */
    MrBFlt          length;                 /*!< length of pending branch                   */
    MrBFlt          nodeDepth;              /*!< node depth (height)                        */
    MrBFlt          age;                    /*!< age of node                                */
    Calibration     *calibration;           /*!< pointer to calibration data                */
    }
    TreeNode;

/* typedef for binary tree */
typedef struct 
    {
    char            name[100];          /*!< name of tree                                 */
    int             memNodes;           /*!< number of allocated nodes (do not exceed!)   */
    int             nNodes;             /*!< number of nodes in tree (including lower root in rooted trees) */
    int             nIntNodes;          /*!< number of interior nodes in tree (excluding lower root in rooted trees) */  
    int             isRooted;           /*!< is tree rooted?                              */
    int             isClock;            /*!< is tree clock?                               */
    int             isCalibrated;       /*!< is tree calibrated?                          */
    int             nRelParts;          /*!< number of relevant partitions                */
    int             *relParts;          /*!< pointer to relevant partitions               */
    int             checkConstraints;   /*!< does tree have constraints?                  */
    int             nConstraints;       /*!< number of constraints                        */
    int             *constraints;       /*!< pointer to constraints                       */
    int             nLocks;             /*!< number of constrained (locked) nodes         */
    TreeNode        **allDownPass;      /*!< downpass array of all nodes                  */
    TreeNode        **intDownPass;      /*!< downpass array of interior nodes (including upper but excluding lower root in rooted trees) */
    TreeNode        *root;              /*!< pointer to root (lower root in rooted trees) */
    TreeNode        *nodes;             /*!< array containing the nodes                   */
    BitsLong        *bitsets;           /*!< pointer to bitsets describing splits         */
    BitsLong        *flags;             /*!< pointer to cond like flags                   */
    int             fromUserTree;       /*!< YES is set for the trees whoes branch lengthes are set from user tree(as start tree or fix branch length prior), NO otherwise */       
    }
    Tree;

/* typedef for node in polytomous tree */
typedef struct pNode
    {
    char            label[100];         /*!< name of node if terminal                     */
    struct pNode    *left, *sib, *anc;  /*!< pointers to adjacent nodes                   */
    int             x, y, mark;         /*!< scratch variables                            */
    int             partitionIndex;     /*!< partition index in sumt (scratch)            */
    int             index;              /*!< index of node (if < numLocalTaxa = local taxon index) */
    int             memoryIndex;        /*!< immutable index of memory position           */
    int             isLocked;           /*!< is the node locked?                          */
    int             lockID;             /*!< id of lock                                   */
    int             isDated;            /*!< is node dated?                               */
    MrBFlt          length;             /*!< age of node                                  */
    MrBFlt          depth;              /*!< depth (height) of node                       */
    MrBFlt          age;                /*!< age of node                                  */
    MrBFlt          support, f;         /*!< scratch variables                            */
    BitsLong        *partition;         /*!< pointer to partition (split) bitset          */
    Calibration     *calibration;       /*!< pointer to dating of node                    */
    }
    PolyNode;

/* typedef for polytomous tree */
typedef struct 
    {
    char            name[100];           /*!< name of tree                                */
    int             memNodes;            /*!< number of allocated nodes; do not exceed!   */
    int             nNodes;              /*!< number of nodes in tree                     */
    int             nIntNodes;           /*!< number of interior nodes in tree            */
    PolyNode        **allDownPass;       /*!< downpass array over all nodes               */
    PolyNode        **intDownPass;       /*!< downpass array over interior nodes          */
    PolyNode        *root;               /*!< pointer to root (lower for rooted trees     */
    PolyNode        *nodes;              /*!< array holding the tree nodes                */  
    BitsLong        *bitsets;            /*!< bits describing partitions (splits)         */
    int             nBSets;              /*!< number of effective branch length sets      */
    int             nESets;              /*!< number of breakpoint rate sets              */
    char            **bSetName;          /*!< names of effective branch length sets       */
    char            **eSetName;          /*!< names of breakpoint rate sets               */
    int             **nEvents;           /*!< number of branch events of bp rate set      */
    MrBFlt          ***position;         /*!< position of branch events                   */
    MrBFlt          ***rateMult;         /*!< parameter of branch events                  */
    MrBFlt          **effectiveBrLen;    /*!< effective branch lengths of ebl set         */
    int             brlensDef;           /*!< are brlens defined ?                        */
    int             isRooted;            /*!< is tree rooted?                             */
    int             isClock;             /*!< is tree clock?                              */
    int             isCalibrated;        /*!< is tree calibrated?                         */
    int             isRelaxed;           /*!< is tree relaxed?                            */
    MrBFlt          clockRate;           /*!< clock rate                                  */
    int             popSizeSet;          /*!< does tree have a population size set?       */
    MrBFlt          *popSize;            /*!< the population size                         */
    char            *popSizeSetName;     /*!< name of the population size set             */
    }
    PolyTree;

/* struct for holding model parameter info for the mcmc run */
typedef struct param
    {
    int             index;              /* index to the parameter (0, 1, 2, ...)          */
    int             paramType;          /* the type of the parameter                      */
    int             paramId;            /* unique ID for parameter x prior combination    */
    MrBFlt          *values;            /* main values of parameter                       */
    MrBFlt          *subValues;         /* subvalues of parameter                         */
    int             *intValues;         /* integer values (model index/growth fxn)        */
    int             nValues;            /* number of values                               */
    int             nSubValues;         /* number of subvalues                            */
    int             nIntValues;         /* number of intvalues                            */
    MrBFlt          min;                /* minimum value of parameter                     */
    MrBFlt          max;                /* maximum value of parameter                     */
    int             *relParts;          /* pointer to relevant divisions                  */
    int             nRelParts;          /* number of relevant divisions                   */
    int             upDate;             /* update flag (for copying)                      */
    struct param    **subParams;        /* pointers to subparams (for topology)           */
    int             nSubParams;         /* number of subparams                            */
    Tree            **tree;             /* pointer to tree ptrs (for brlens & topology)   */
    int             treeIndex;          /* index to first tree in mcmcTree                */
    int             hasBinaryStd;       /* has binary standard chars                      */
    int             *sympiBsIndex;      /* pointer to sympi bsIndex (std chars)           */
    int             *sympinStates;      /* pointer to sympi nStates (std chars)           */
    int             *sympiCType;        /* pointer to sympi cType (std chars)             */
    int             nSympi;             /* number of sympis                               */
    int             printParam;         /* whether parameter should be printed            */
    int             nPrintSubParams;    /* number of subparams that should be printed     */
    char            *paramHeader;       /* a string holding header for param values       */
    char            *name;              /* string holding name of parameter               */
    char            *paramTypeName;     /* pointer to description of parameter type       */
    int             checkConstraints;   /* is tree parameter constrained?                 */
    int             fill;               /* flags whether the parameter should be filled   */
    int             nStdStateFreqs;     /* number of std state frequencies                */
    MrBFlt          *stdStateFreqs;     /* pointer to std state frequencies               */
    int             **nEvents;          /* number of branch events for Cpp model          */
                                        /* nEvents[0..2*numCains][0..numNodes=2*numTaxa]  */
    MrBFlt          ***position;        /* event positions for Cpp relaxed clock model    */
    MrBFlt          ***rateMult;        /* rate multipliers for Cpp relaxed clock model   */
    int             affectsLikelihood;  /* does parameter directly influence likelihood?  */
    MrBFlt*         priorParams;        /* pointer to the prior parameters                */
    LnPriorProbFxn  LnPriorProb;        /* ln prior prob function                         */
    LnPriorRatioFxn LnPriorRatio;       /* ln prior prob ratio function                   */
    } Param;

#if defined(THREADS_ENABLED)
#include <pthread.h>

typedef struct s_launch_struct 
    {
    int chain;
    int division;
    MrBFlt* lnL;                    
    } LaunchStruct; 
#endif

/* parameter ID values */
/* identifies unique model parameter x prior combinations */
#define TRATIO_DIR                      1
#define TRATIO_FIX                      2
#define REVMAT_DIR                      3
#define REVMAT_FIX                      4
#define OMEGA_DIR                       5
#define OMEGA_FIX                       6
#define SYMPI_UNI                       7
#define SYMPI_UNI_MS                    8
#define SYMPI_EXP                       9
#define SYMPI_EXP_MS                    10
#define SYMPI_FIX                       11
#define SYMPI_FIX_MS                    12
#define SYMPI_EQUAL                     13
#define PI_DIR                          14
#define PI_USER                         15
#define PI_EMPIRICAL                    16
#define PI_EQUAL                        17
#define PI_FIXED                        18
#define SHAPE_UNI                       19
#define SHAPE_EXP                       20
#define SHAPE_FIX                       21
#define PINVAR_UNI                      22
#define PINVAR_FIX                      23
#define CORREL_UNI                      24
#define CORREL_FIX                      25
#define SWITCH_UNI                      26
#define SWITCH_EXP                      27
#define SWITCH_FIX                      28
#define RATEMULT_DIR                    29
#define RATEMULT_FIX                    30
#define TOPOLOGY_NCL_UNIFORM            31
#define TOPOLOGY_NCL_CONSTRAINED        32
#define TOPOLOGY_NCL_FIXED              33
#define TOPOLOGY_NCL_UNIFORM_HOMO       34
#define TOPOLOGY_NCL_UNIFORM_HETERO     35
#define TOPOLOGY_NCL_CONSTRAINED_HOMO   36
#define TOPOLOGY_NCL_CONSTRAINED_HETERO 37
#define TOPOLOGY_NCL_FIXED_HOMO         38
#define TOPOLOGY_NCL_FIXED_HETERO       39
#define TOPOLOGY_CL_UNIFORM             40
#define TOPOLOGY_CL_CONSTRAINED         41
#define TOPOLOGY_CL_FIXED               42
#define TOPOLOGY_CCL_UNIFORM            43
#define TOPOLOGY_CCL_CONSTRAINED        44
#define TOPOLOGY_CCL_FIXED              45
#define TOPOLOGY_PARSIMONY_UNIFORM      46
#define TOPOLOGY_PARSIMONY_CONSTRAINED  47
#define TOPOLOGY_PARSIMONY_FIXED        48
#define BRLENS_UNI                      49
#define BRLENS_EXP                      50
#define BRLENS_GamDir                   51
#define BRLENS_iGmDir                   52
#define BRLENS_twoExp                   53
#define BRLENS_FIXED                    54
#define BRLENS_CLOCK_UNI                55
#define BRLENS_CLOCK_COAL               56
#define BRLENS_CLOCK_BD                 57
#define BRLENS_CLOCK_FIXED              58
#define BRLENS_CLOCK_SPCOAL             59
#define BRLENS_CLOCK_FOSSIL             60
#define BRLENS_PARSIMONY                61
#define SPECRATE_UNI                    62
#define SPECRATE_EXP                    63
#define SPECRATE_FIX                    64
#define EXTRATE_BETA                    65
#define EXTRATE_FIX                     66
#define FOSLRATE_BETA                   67
#define FOSLRATE_FIX                    68
#define POPSIZE_UNI                     69
#define POPSIZE_GAMMA                   70
#define POPSIZE_FIX                     71
#define POPSIZE_NORMAL                  72
#define POPSIZE_LOGNORMAL               73
#define AAMODEL_FIX                     74
#define AAMODEL_MIX                     75
#define GROWTH_UNI                      76
#define GROWTH_EXP                      77
#define GROWTH_FIX                      78
#define GROWTH_NORMAL                   79
#define OMEGA_BUD                       80
#define OMEGA_BUF                       81
#define OMEGA_BED                       82
#define OMEGA_BEF                       83
#define OMEGA_BFD                       84
#define OMEGA_BFF                       85
#define OMEGA_FUD                       86
#define OMEGA_FUF                       87
#define OMEGA_FED                       88
#define OMEGA_FEF                       89
#define OMEGA_FFD                       90
#define OMEGA_FFF                       91
#define OMEGA_ED                        92
#define OMEGA_EF                        93
#define OMEGA_FD                        94
#define OMEGA_FF                        95
#define OMEGA_10UUB                     96
#define OMEGA_10UUF                     97
#define OMEGA_10UEB                     98
#define OMEGA_10UEF                     99
#define OMEGA_10UFB                     100
#define OMEGA_10UFF                     101
#define OMEGA_10EUB                     102
#define OMEGA_10EUF                     103
#define OMEGA_10EEB                     104
#define OMEGA_10EEF                     105
#define OMEGA_10EFB                     106
#define OMEGA_10EFF                     107
#define OMEGA_10FUB                     108
#define OMEGA_10FUF                     109
#define OMEGA_10FEB                     110
#define OMEGA_10FEF                     111
#define OMEGA_10FFB                     112
#define OMEGA_10FFF                     113
#define CPPRATE_FIX                     114
#define CPPRATE_EXP                     115
#define CPPMULTDEV_FIX                  116
#define TK02VAR_FIX                     117
#define TK02VAR_EXP                     118
#define TK02VAR_UNI                     119
#define TOPOLOGY_RCL_UNIFORM            120
#define TOPOLOGY_RCL_CONSTRAINED        121
#define TOPOLOGY_RCL_FIXED              122
#define TOPOLOGY_RCCL_UNIFORM           123
#define TOPOLOGY_RCCL_CONSTRAINED       124
#define TOPOLOGY_RCCL_FIXED             125
#define TOPOLOGY_SPECIESTREE            126
#define CPPEVENTS                       127
#define TK02BRANCHRATES                 128
#define TOPOLOGY_FIXED                  129
#define IGRVAR_FIX                      130
#define IGRVAR_EXP                      131
#define IGRVAR_UNI                      132
#define IGRBRANCHRATES                  133
#define CLOCKRATE_FIX                   134
#define CLOCKRATE_NORMAL                135
#define CLOCKRATE_LOGNORMAL             136
#define CLOCKRATE_GAMMA                 137
#define CLOCKRATE_EXP                   138
#define SPECIESTREE_UNIFORM             139
#define GENETREERATEMULT_DIR            140
#define GENETREERATEMULT_FIX            141
#define REVMAT_MIX                      142
#define MIXEDVAR_FIX                    143
#define MIXEDVAR_EXP                    144
#define MIXEDVAR_UNI                    145
#define MIXEDBRCHRATES                  146

#if defined (BEAGLE_ENABLED)
#define MB_BEAGLE_SCALE_ALWAYS          0
#define MB_BEAGLE_SCALE_DYNAMIC         1
#if defined (_DEBUG)
#define MB_PRINT_DYNAMIC_RESCALE_FAIL_STAT
#endif
#endif

/* typedef for a MoveFxn */
typedef int (MoveFxn)(Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);

/* typedef for an ApplicFxn */
typedef int (ApplicFxn)(Param *param);

/* typedef for an AutotuneFxn */
typedef void (AutotuneFxn)(MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *tuningParameter, MrBFlt minTuning, MrBFlt maxTuning);

/* struct holding info on each move type that the program handles */
typedef struct
    {
    MoveFxn      *moveFxn;           /* pointer to the move function                 */
    ApplicFxn    *isApplicable;      /* pointer to function determining whether move is applicable to a parameter */
    int          nApplicable;        /* number of relevant params                    */
    int          applicableTo[40];   /* pointer to ID of relevant params             */
    char         *name;              /* name of the move type                        */
    char         *shortName;         /* abbreviated name of the move type            */
    char         *paramName;         /* name of subparameter if complex parameter    */
    int          subParams;          /* are we changing subparams (brlens of topol.) */
    char         *tuningName[5];     /* name of tuning params                        */
    char         *shortTuningName[5];/* short name of tuning params                  */
    MrBFlt       relProposalProb;    /* default relative proposal probability        */
    int          numTuningParams;    /* number of tuning parameters                  */
    MrBFlt       tuningParam[5];     /* default tuning parameters for the proposal   */
    MrBFlt       minimum[5];         /* minimum values for tuning params             */
    MrBFlt       maximum[5];         /* maximum values for tuning params             */
    int          parsimonyBased;     /* this move is based on parsimony (YES/NO)     */
    int          level;              /* user level of this move                      */
    AutotuneFxn  *Autotune;          /* pointer to the autotuning function           */
    MrBFlt       targetRate;         /* default target acceptance rate for autotuning */
    } MoveType;

/* max number of move types */
#define NUM_MOVE_TYPES 101

/* struct holding info on each move */
/* Note: This allows several different moves to affect the same parameter */
/* It also allows the same move to affect different parameters as before */
/* This is also a good place to keep the proposal probs */
typedef struct
    {
    char        *name;              /* name of the move                             */
    MoveType    *moveType;          /* pointer to the move type                     */
    MoveFxn     *moveFxn;           /* pointer to the move function                 */
    Param       *parm;              /* ptr to parameter the move applies to         */
    MrBFlt      *relProposalProb;   /* the relative proposal probability            */
    MrBFlt      *cumProposalProb;   /* the cumulative proposal probability          */
    int         *nAccepted;         /* number of accepted moves                     */
    int         *nTried;            /* number of tried moves                        */
    int         *nBatches;          /* counter for autotuning rounds                */
    int         *nTotAccepted;      /* total number of accepted moves               */
    int         *nTotTried;         /* total number of tried moves                  */
    MrBFlt      *targetRate;        /* target acceptance rate for autotuning        */
    MrBFlt      *lastAcceptanceRate;/* acceptance rate in last complete batch       */
    MrBFlt      **tuningParam;      /* tuning parameters for the move               */
    } MCMCMove;

typedef int (*LikeDownFxn)(TreeNode *, int, int);
typedef int (*LikeRootFxn)(TreeNode *, int, int);
typedef int (*LikeScalerFxn)(TreeNode *, int, int);
typedef int (*LikeFxn)(TreeNode *, int, int, MrBFlt *, int);
typedef int (*TiProbFxn)(TreeNode *, int, int);
typedef int (*LikeUpFxn)(TreeNode *, int, int);
typedef int (*PrintAncStFxn)(TreeNode *, int, int);
typedef int (*StateCodeFxn) (int);
typedef int (*PrintSiteRateFxn) (TreeNode *, int, int);
typedef int (*PosSelProbsFxn) (TreeNode *, int, int);
typedef int (*SiteOmegasFxn) (TreeNode *, int, int);

typedef struct cmdtyp           
    {
    int         cmdNumber;
    char        *string;
    int         specialCmd;
    CmdFxn      cmdFxnPtr;
    short       numParms;
    short       parmList[50];
    int         expect;
    char        *cmdDescription;
    int         cmdUse;
    int         hiding;
    }
    CmdType;
    
typedef struct parm
    {
    char        *string;    /* parameter name */
    char        *valueList; /* list of values that could be input */
    ParmFxn     fp;         /* function pointer */
    }
    ParmInfo, *ParmInfoPtr;

typedef struct model
    {
    int         dataType;          /* data type for partition                      */
    int         nStates;           /* number of states for this type of data       */
    int         codon[64];         /* gives protein ID for each codon              */
    int         codonNucs[64][3];  /* gives triplet for each codon                 */
    int         codonAAs[64];      /* gives protein ID for implemented code        */
    
    char        nucModel[100];     /* nucleotide model used                        */
    char        nst[100];          /* number of substitution types                 */
    char        parsModel[100];    /* use the (so-called) parsimony model          */
    char        geneticCode[100];  /* genetic code used                            */
    int         coding;            /* type of patterns encoded                     */
    char        codingString[100]; /* string describing type of patterns encoded   */
    char        ploidy[100];       /* ploidy level                                 */
    char        omegaVar[100];     /* type of omega variation model                */
    char        ratesModel[100];   /* rates across sites model                     */
    int         numGammaCats;      /* number of categories for gamma approximation */
    char        useGibbs[100];     /* flags whether Gibbs sampling of discrete gamma is used */
    int         gibbsFreq;         /* frequency of Gibbs resampling of discrete gamma */

    int         numBetaCats;       /* number of categories for beta approximation  */
    int         numM10GammaCats;   /* number of cats for gamma approx (M10 model)  */
    int         numM10BetaCats;    /* number of cats for beta approx (M10 model)   */
    char        covarionModel[100];/* use covarion model? (yes/no)                 */
    char        augmentData[100];  /* should data be augmented                     */

    char        tRatioPr[100];     /* prior for ti/tv rate ratio                   */
    MrBFlt      tRatioFix;   
    MrBFlt      tRatioDir[2];      
    char        revMatPr[100];     /* prior for GTR model                          */
    MrBFlt      revMatFix[6];
    MrBFlt      revMatDir[6];
    MrBFlt      revMatSymDir;      /* prior for mixed GTR subspace model           */
    char        aaModelPr[100];    /* prior for amino acid model                   */
    char        aaModel[100];
    MrBFlt      aaModelPrProbs[10];
    char        aaRevMatPr[100];   /* prior for aa GTR model                       */
    MrBFlt      aaRevMatFix[190];
    MrBFlt      aaRevMatDir[190];
    char        omegaPr[100];      /* prior for omega                              */
    MrBFlt      omegaFix;
    MrBFlt      omegaDir[2];
    char        ny98omega1pr[100]; /* prior for class 1 omega (Ny98 model)         */
    MrBFlt      ny98omega1Fixed;
    MrBFlt      ny98omega1Beta[2];
    char        ny98omega3pr[100]; /* prior for class 3 omega (Ny98 model)         */
    MrBFlt      ny98omega3Fixed;
    MrBFlt      ny98omega3Uni[2];
    MrBFlt      ny98omega3Exp;
    char        m3omegapr[100];    /* prior for all three omegas (M3 model)        */
    MrBFlt      m3omegaFixed[3];
    char        m10betapr[100];    /* prior for omega variation (M10 model)        */
    char        m10gammapr[100];
    MrBFlt      m10betaExp;
    MrBFlt      m10betaUni[2];
    MrBFlt      m10betaFix[2];
    MrBFlt      m10gammaExp;
    MrBFlt      m10gammaUni[2];
    MrBFlt      m10gammaFix[2];
    char        codonCatFreqPr[100];  /* prior for selection cat frequencies       */
    MrBFlt      codonCatFreqFix[3];
    MrBFlt      codonCatDir[3];
    char        stateFreqPr[100];  /* prior for character state frequencies        */
    MrBFlt      stateFreqsFix[200];
    MrBFlt      stateFreqsDir[200];
    char        stateFreqsFixType[100];
    int         numDirParams;
    char        shapePr[100];      /* prior for gamma/lnorm shape parameter        */
    MrBFlt      shapeFix;
    MrBFlt      shapeUni[2];
    MrBFlt      shapeExp;
    char        pInvarPr[100];     /* prior for proportion of invariable sites     */
    MrBFlt      pInvarFix;
    MrBFlt      pInvarUni[2];
    char        adGammaCorPr[100]; /* prior for correlation param of adGamma model */
    MrBFlt      corrFix;
    MrBFlt      corrUni[2];
    char        covSwitchPr[100];  /* prior for switching rates of covarion model  */
    MrBFlt      covswitchFix[2];
    MrBFlt      covswitchUni[2];
    MrBFlt      covswitchExp;
    char        symPiPr[100];      /* prior for pi when unidentifiable states used */
    MrBFlt      symBetaFix;
    MrBFlt      symBetaUni[2];
    MrBFlt      symBetaExp;
    char        ratePr[100];       /* prior on rate for a partition                */
    MrBFlt      ratePrDir;
    char        generatePr[100];   /* prior on rate for a gene (one or more partitions) */
    MrBFlt      generatePrDir;
    char        brownCorPr[100];   /* prior for correlation of Brownian model      */
    MrBFlt      brownCorrFix;
    MrBFlt      brownCorrUni[2];
    char        brownScalesPr[100];  /* prior for scales of Brownian model         */
    MrBFlt      brownScalesFix;
    MrBFlt      brownScalesUni[2];
    MrBFlt      brownScalesGamma[2];
    MrBFlt      brownScalesGammaMean;

    char        topologyPr[100];   /* prior for tree topology                      */
    int         topologyFix;       /* user tree index for fixed topology           */
    int         *activeConstraints;  /* which constraints are active?              */
    int         numActiveConstraints;
    int         numActiveLocks;
    char        brlensPr[100];     /* prior on branch lengths                      */
    int         brlensFix;         /* user tree index for fixed brlens             */
    MrBFlt      brlensUni[2];
    MrBFlt      brlensExp;
    MrBFlt      brlens2Exp[2];
    MrBFlt      brlensDir[4];
    MrBFlt      brlensGamma[2];
    char        speciesTreeBrlensPr[100];  /* prior on branch lengths of species tree       */
    char        unconstrainedPr[100];  /* prior on branch lengths if unconstrained          */
    char        clockPr[100];          /* prior on branch if clock enforced                 */
    char        clockVarPr[100];       /* prior on clock rate variation (strict, cpp, tk02, igr, ...) */
    char        nodeAgePr[100];        /* prior on node depths (unconstrained, constraints) */
    char        speciationPr[100];     /* prior on speciation rate (net diversification)    */
    MrBFlt      speciationFix;
    MrBFlt      speciationUni[2];
    MrBFlt      speciationExp;
    char        extinctionPr[100];     /* prior on relative extinction rate (turnover)      */
    MrBFlt      extinctionFix;
    MrBFlt      extinctionBeta[2];
    char        fossilizationPr[100];  /* prior on fossilization rate (sampling proportion) */
    MrBFlt      fossilizationFix;
    MrBFlt      fossilizationBeta[2];
    char        sampleStrat[100];      /* taxon sampling strategy (for b-d process)         */
    int         sampleFSNum;           /* number of fossil slice sampling events (s)        */
    MrBFlt     *sampleFSTime;          /* fossil slice sampling times (t_i,   i=1,..,s)     */
    MrBFlt     *sampleFSProb;          /* fossil slice sampling probs (rho_i, i=1,..,s)     */
    MrBFlt      sampleProb;            /* extant taxon sampling fraction (rho)              */
    Calibration treeAgePr;             /* prior on tree age for uniform clock trees         */
    char        clockRatePr[100];      /* prior on base substitution rate of tree for clock trees */
    MrBFlt      clockRateNormal[2];
    MrBFlt      clockRateLognormal[2];
    MrBFlt      clockRateGamma[2];
    MrBFlt      clockRateExp;
    MrBFlt      clockRateFix;
    char        popSizePr[100];        /* prior on population size                      */
    MrBFlt      popSizeFix;
    MrBFlt      popSizeUni[2];
    MrBFlt      popSizeLognormal[2];
    MrBFlt      popSizeGamma[2];
    MrBFlt      popSizeNormal[2];
    char        popVarPr[100];         /* prior on pop. size variation across tree      */
    char        growthPr[100];         /* prior on coalescence growth rate              */
    MrBFlt      growthFix;
    MrBFlt      growthUni[2];
    MrBFlt      growthExp;
    MrBFlt      growthNorm[2];
    char        cppRatePr[100];        /* prior on CPP rate                             */
    MrBFlt      cppRateFix;
    MrBFlt      cppRateExp;
    char        cppMultDevPr[100];     /* prior on CPP rate multiplier Lognormal variance */
    MrBFlt      cppMultDevFix;
    char        tk02varPr[100];        /* prior on TK02 lognormal rate variance         */
    MrBFlt      tk02varFix;
    MrBFlt      tk02varUni[2];
    MrBFlt      tk02varExp;
    char        igrvarPr[100];         /* prior on IGR gamma distribution variance      */
    MrBFlt      igrvarFix;
    MrBFlt      igrvarUni[2];
    MrBFlt      igrvarExp;
    char        mixedvarPr[100];       /* prior on mixed relaxed clock rate variance    */
    MrBFlt      mixedvarFix;
    MrBFlt      mixedvarUni[2];
    MrBFlt      mixedvarExp;

    char        tratioFormat[30];      /* format used to report tratio                  */
    char        revmatFormat[30];      /* format used to report revmat                  */
    char        ratemultFormat[30];    /* format used to report ratemult                */
    char        treeFormat[30];        /* format used to report trees/topologies        */
    char        inferAncStates[5];     /* should ancestral states be inferred (Yes/No)? */
    char        inferSiteOmegas[5];    /* should site omega vals be inferred (Yes/No)?  */
    char        inferSiteRates[5];     /* should site rates be inferred (Yes/No)?       */
    char        inferPosSel[5];        /* should site selection be inferred (Yes/No)?   */
    } Model, ModelParams;

typedef struct chain
    {
    int         numGen;                /* number of MCMC cycles                         */
    int         sampleFreq;            /* frequency to sample chain                     */
    int         printFreq;             /* frequency to print chain                      */
    int         swapFreq;              /* frequency to attempt swap of states           */
    int         numRuns;               /* number of runs                                */
    int         numChains;             /* number of chains                              */
    int         isSS;                  /* do we do Steppingstone Sampling               */
    int         startFromPriorSS;      /* If Yes SS is moving from Prior to Posterior   */
    int         numStepsSS;            /* Number of steps in SS                         */
    int         burninSS;              /* Fixed burnin for SS                           */
    MrBFlt      alphaSS;               /* Beta values are distributed according to quantiles of Beta(alphaSS,1.0) distribution */
    int         backupCheckSS;         /* Frequency of checkpoints backup               */
    MrBFlt      chainTemp;             /* chain temperature                             */
    int         userDefinedTemps;      /* should we use the users temperatures?         */
    MrBFlt      userTemps[MAX_CHAINS]; /* user-defined chain temperatures               */
    char        chainFileName[100];    /* chain file name for output                    */
    int         chainBurnIn;           /* chain burn in length                          */
    int         numStartPerts;         /* number of perturbations to starting tree      */
    char        startTree[100];        /* starting tree for chain (current/random)      */
    char        startParams[100];      /* starting values for chain (current/reset)     */
    int         saveBrlens;            /* should branch lengths be saved                */
    MrBFlt      weightScheme[3];       /* percent chars to increase/decrease in weight  */
    int         calcPbf;               /* should we calculate the pseudo Bayes factor   */
    int         pbfInitBurnin;         /* initial burnin when calculating pseudo BF     */
    int         pbfSampleFreq;         /* sample frequency for pseudo BF                */
    int         pbfSampleTime;         /* how many cycles to calcualate site prob.      */
    int         pbfSampleBurnin;       /* burnin period for each site for pseudo BF     */
    int         swapAdjacentOnly;      /* whether we only swap adjacent temperatures    */
    int         redirect;              /* should output go to stdout                    */
    int         allChains;             /* should stats be output for all chains?        */
    int         numSwaps;              /* number of swaps to try each time              */
    int         mcmcDiagn;             /* should mcmc diagnostics be output?            */
    int         diagnFreq;             /* mcmc diagnostics frequency                    */
    int         diagnStat;             /* statistic to use for mcmc diagnostics         */
    int         relativeBurnin;        /* should a relative burnin be used ?            */
    MrBFlt      burninFraction;        /* the sample fraction to discard as burnin      */
    int         allComps;              /* top conv diagnosis for all pairs?             */
    MrBFlt      minPartFreq;           /* minimum partition frequency for conv diagn    */
    MrBFlt      stopVal;               /* top conv diagn value to reach before stopping */
    int         stopRule;              /* use stop rule?                                */
    STATS       *stat;                 /* ptr to structs with mcmc diagnostics info     */
    Tree        *dtree;                /* pointing to tree used for conv diagnostics    */
    TreeList    *treeList;             /* vector of tree lists for saving trees         */
    int         saveTrees;             /* save tree samples for later removal?          */
    int         stopTreeGen;           /* generation after which no trees need be saved */
    fpos_t      *tFilePos;             /* position for reading trees for removal        */
    int         printMax;              /* maximum number of chains to print             */
    int         printAll;              /* whether to print all or only cold chains      */
    int         checkPoint;            /* should we use check-pointing?                 */
    int         checkFreq;             /* check-pointing frequency                      */
    int         runWithData;           /* should we run with or without data?           */
    int         orderTaxa;             /* order taxa before printing tree to file?      */
    int         append;                /* order taxa before printing tree to file?      */
    int         autotune;              /* autotune tuning parameters of proposals ?     */
    int         tuneFreq;              /* autotuning frequency                          */
    } Chain;

typedef struct modelinfo
    {
    /* General model information */
    int         dataType;                   /* data type for partition                  */
    int         nucModelId;                 /* structure of nucleotide model            */
    int         nst;                        /* # substitution types                     */
    int         aaModelId;                  /* amino acid model type                    */
    int         parsModelId;                /* is parsimony model used YES/NO           */

    /* Specific model information */
    int         numGammaCats;               /* number of gamma cats (1 if inapplic.)    */
    int         numBetaCats;                /* number of beta cats (1 if inapplic.)     */
    int         numOmegaCats;               /* number of omega cats (1 if inapplic.)    */
    int         numTiCats;                  /* number of cats needing different tis     */
    int         numModelStates;             /* number of states including hidden ones   */
    int         numStates;                  /* number of observable discrete states     */

    /* Model parameter pointers */
    Param       *tRatio;                    /* ptr to tRatio used in model              */
    Param       *revMat;                    /* ptr to revMat used in model              */
    Param       *omega;                     /* ptr to omega used in model               */
    Param       *stateFreq;                 /* ptr to statFreq used in model            */
    Param       *shape;                     /* ptr to shape used in model               */
    Param       *pInvar;                    /* ptr to pInvar used in model              */
    Param       *correlation;               /* ptr to correlation used in model         */
    Param       *switchRates;               /* ptr to switchRates (off<->on)            */
    Param       *rateMult;                  /* ptr to parition rateMult used in model   */
    Param       *geneTreeRateMult;          /* ptr to gene tree rateMult used in model  */
    Param       *speciationRates;           /* ptr to speciationRates used in model     */
    Param       *extinctionRates;           /* ptr to extinctionRates used in model     */
    Param       *fossilizationRates;        /* ptr to fossilizationRates                */
    Param       *popSize;                   /* ptr to population size used in model     */
    Param       *growthRate;                /* ptr to growth rate used in model         */
    Param       *topology;                  /* ptr to topology used in model            */
    Param       *brlens;                    /* ptr to brlens (and tree) used in model   */
    Param       *speciesTree;               /* ptr to species tree used in model        */
    Param       *aaModel;                   /* ptr to amino acid matrix used            */
    Param       *cppMultDev;                /* ptr to cpp ratemult lognormal variance   */
    Param       *cppRate;                   /* ptr to CPP rate used in model            */
    Param       *cppEvents;                 /* ptr to CPP events                        */
    Param       *tk02var;                   /* ptr to variance for TK02 relaxed clock   */
    Param       *tk02BranchRates;           /* ptr to branch rates for TK02 relaxed clock */
    Param       *igrvar;                    /* ptr to gamma var for IGR relaxed clock   */
    Param       *igrBranchRates;            /* ptr to branch rates for IGR relaxed clock*/
    Param       *mixedvar;                  /* ptr to var for mixed relaxed clock       */
    Param       *mixedBrchRates;            /* ptr to branch rates for mixed relaxed clock */
    Param       *clockRate;                 /* ptr to clock rate parameter              */

    /* Information about characters and transformations */
    int         numChars;                   /* number of compressed characters          */
    int         numUncompressedChars;       /* number of uncompressed characters        */
    int         numDummyChars;              /* number of dummy characters               */
    int         compMatrixStart;            /* start column in compressed matrix        */
    int         compMatrixStop;             /* stop column in compressed matrix         */
    int         compCharStart;              /* start char among compressed chars        */
    int         compCharStop;               /* stop char among compressed chars         */
    int         parsMatrixStart;            /* start column in parsimony matrix         */
    int         parsMatrixStop;             /* stop collumn in parsimony matrix         */
    int         nParsIntsPerSite;           /* # parsimony ints per character           */  
    int         nCharsPerSite;              /* number chars per site (eg 3 for codon)   */
    int         rateProbStart;              /* start of rate probs (for adgamma)        */
                
     /* Variables for eigen decompositions */
    int         cijkLength;                 /* stores length of cijk vector                 */
    int         nCijkParts;                 /* stores number of cijk partitions (1 except for omega/covarion models) */
    int         upDateCijk;                 /* whether cijk vector needs to be updated      */

    /* Variables for standard model */
    int         *tiIndex;                   /* index to trans probs for each compressed char*/
    int         *bsIndex;                   /* index to stat freqs for each compressed char */
    int         *nStates;                   /* # states of each compressed char             */
    int         *cType;                     /* whether char is ord, unord or irrev          */
    int         *weight;                    /* prior weight of each compressed char         */
    int         isTiNeeded[20];             /* marks whether a trans prob matrix is needed  */

    /* Gibbs sampling of gamma site rate parameters */
    CLFlt       ***catLike;                 /* likelihood for Gibbs sampling of gamma       */
    CLFlt       ***catLnScaler;             /* scaler for Gibbs sampling of gamma           */
    int         gibbsGamma;                 /* flags whether Gibbs sampling of discrete gamma is used */
    int         gibbsFreq;                  /* frequency of Gibbs resampling of discrete gamma */
    
    /* Variables for parsimony sets and parsimony calculations */
    MrBFlt      parsTreeLength[MAX_CHAINS*2];   /* parsimony tree lengths for chains        */
    BitsLong    **parsSets;                 /* parsimony sets                               */
    int         numParsSets;                /* number of parsimony sets                     */
    CLFlt       *parsNodeLens;              /* parsimony node lengths                       */
    int         numParsNodeLens;            /* number of parsimony node lengths             */

    /* Miscellaneous parameters */
    int         mark;                       /* scratch parameter                            */
    int         parsimonyBasedMove;         /* is parsimony-based move used (YES/NO)        */

    /* Variables for conditional likelihoods */
    int         upDateCl;                   /* flags whether update of cond likes needed    */
    int         upDateAll;                  /* flags whether update of entire tree is needed*/
    int         *isPartAmbig;               /* is tip partially ambiguous?                  */
    int         **termState;                /* index arrays for terminal branch shortcuts   */
    CLFlt       *invCondLikes;              /* space for the invariable cond likes          */
    CLFlt       **condLikes;                /* space for the cond likes                     */
    CLFlt       **tiProbs;                  /* space for the ti probs                       */
    CLFlt       **scalers;                  /* space for the node and site scalers          */
    CLFlt       **clP;                      /* handy pointers to cond likes for ti cats     */
#if defined (SSE_ENABLED)
    __m128      **clP_SSE;                  /* handy pointers to cond likes, SSE version    */
    int         numVecChars;                /* number of compact SIMD vectors               */
    int         numFloatsPerVec;            /* number of floats per SIMD vector             */
    CLFlt       *lnL_Vec;                   /* temp storage for log site likes              */
    CLFlt       *lnLI_Vec;                  /* temp storage for log site invariable likes   */
#if defined (AVX_ENABLED)
    __m256      **clP_AVX;                  /* handy pointers to cond likes, AVX version    */
#endif
#endif
    MrBFlt      **cijks;                    /* space for cijks                              */
    int         **condLikeIndex;            /* index to cond like space for nodes & chains  */
    int         *condLikeScratchIndex;      /* index to scratch space for node cond likes   */
    int         **unscaledNodes;            /* # unscaled nodes at a node                   */
    int         *unscaledNodesScratch;      /* scratch # unscaled nodes at a node           */
    int         **nodeScalerIndex;          /* index to scaler space for nodes & chains     */
    int         *nodeScalerScratchIndex;    /* index to scratch space for node scalers      */
    int         *siteScalerIndex;           /* index to site scaler space for chains        */
    int         siteScalerScratchIndex;     /* index to scratch space for site scalers      */
    int         **tiProbsIndex;             /* index to ti probs for branches & chains      */
    int         *tiProbsScratchIndex;       /* index to scratch space for branch ti probs   */
    int         *cijkIndex;                 /* index to cijks for chains                    */
    int         cijkScratchIndex;           /* index to scratch space for cijks             */
    int         numCondLikes;               /* number of cond like arrays                   */
    int         numScalers;                 /* number of scaler arrays                      */
    int         numTiProbs;                 /* number of ti prob arrays                     */
    int         condLikeLength;             /* length of cond like array (incl. ti cats)    */
    int         tiProbLength;               /* length of ti prob array                      */
    MrBFlt      lnLike[MAX_CHAINS];         /* log like for chain                           */
    CLFlt       *ancStateCondLikes;         /* ancestral state cond like array              */

    /* Likelihood function pointers */
    LikeDownFxn         CondLikeDown;       /* function for calculating partials            */
    LikeRootFxn         CondLikeRoot;       /* function for calculating partials at root    */
    LikeScalerFxn       CondLikeScaler;     /* function for scaling partials                */
    LikeFxn             Likelihood;         /* function for getting cond likes for tree     */
    TiProbFxn           TiProbs;            /* function for calculating transition probs    */
    LikeUpFxn           CondLikeUp;         /* final-pass calculation of cond likes         */
    PrintAncStFxn       PrintAncStates;     /* function for sampling ancestral states       */
    StateCodeFxn        StateCode;          /* function for getting states from codes       */
    PrintSiteRateFxn    PrintSiteRates;     /* function for samling site rates              */
    PosSelProbsFxn      PosSelProbs;        /* function for sampling pos. selection probs   */
    SiteOmegasFxn       SiteOmegas;         /* function for sampling site omega values      */

    /* Report variables */
    int         printAncStates;             /* should ancestral states be printed (YES/NO)  */
    int         printSiteRates;             /* should site rates be printed (YES/NO)        */
    int         printPosSel;                /* should selection be printed (YES/NO)         */
    int         printSiteOmegas;            /* should site omegas be printed (YES/NO)       */

    /* likelihood calculator flags and variables */
    int         useBeagle;                  /* use Beagle for this partition?               */
    int         useVec;                     /* use SSE for this partition?                  */
    int*        rescaleFreq;                /* rescale frequency for each chain             */

#if defined (BEAGLE_ENABLED)
    /* Beagle variables */
    int         useBeagleResource;          /* try to use this BEAGLE resource number       */
    MrBFlt*     branchLengths;              /* array of branch lengths for Beagle           */
    MrBFlt*     inRates;                    /* array of category rates for Beagle           */
    int*        tiProbIndices;              /* array of trans prob indices for Beagle       */
    MrBFlt*     logLikelihoods;             /* array of log likelihoods from Beagle         */
    int         beagleInstance;             /* beagle instance for division                 */
    MrBFlt*     inWeights;                  /* array of weights for Beagle root likelihood  */
    int*        bufferIndices;              /* array of partial indices for root likelihood */
    int*        eigenIndices;               /* array of eigen indices for root likelihood   */
    int*        childBufferIndices;         /* array of child partial indices (unrooted)    */
    int*        childTiProbIndices;         /* array of child ti prob indices (unrooted)    */
    int*        cumulativeScaleIndices;     /* array of cumulative scale indices            */
    int         rescaleBeagleAll;           /* set to rescale all nodes                     */
    int         rescaleFreqOld;             /* holds rescale frequency of current state     */
    int         recalculateScalers;         /* shoud we recalculate scalers for current state YES/NO */
    int*        succesCount;                /* count number of succesful computation since last reset of scalers */
    int**       isScalerNode;               /* for each node and chain set to YES if scaled node */
    int*        isScalerNodeScratch;        /* scratch space to hold isScalerNode of proposed state*/
    long*       beagleComputeCount;         /* count of number of calls to likelihood       */
#endif

    } ModelInfo;

typedef struct sumt
    {
    int        *absentTaxa;            /* information on absent taxa                    */
    int         brlensDef;             /* branch lengths defined?                       */
    char        sumtFileName[100];     /* name of input file                            */
    char        sumtOutfile[120];      /* name of output file                           */
    char        curFileName[120];      /* name of file being processed                  */
    int         burnin;                /* actual burnin when parsing tree files         */
    char        sumtConType[100];      /* consensus tree type                           */
    int         calcTreeprobs;         /* should the individual tree probs be calculated*/
    int         showSumtTrees;         /* should the individual tree probs be shown     */
    int         numRuns;               /* number of independent analyses to summarize   */
    int         numTrees;              /* number of tree params to summarize            */
    int         orderTaxa;             /* order taxa in trees?                          */
    MrBFlt      minPartFreq;           /* minimum part. freq. for overall diagnostics   */
    int         table;                 /* show table of partition frequencies?          */
    int         summary;               /* show summary diagnostics ?                    */
    int         showConsensus;         /* show consensus trees ?                        */
    int         consensusFormat;       /* format of consensus tree                      */
    PolyTree   *tree;                  /* for storing tree read from file               */
    int        *order;                 /* for storing topology read from file           */
    int         orderLen;              /* length of order array                         */
    int         numTreesInLastBlock;   /* number of trees in last block                 */
    int         numTreesEncountered;   /* number of trees encounted in total            */
    int         numTreesSampled;       /* number of sampled trees in total              */
    int         isRooted;              /* is sumt tree rooted ?                         */
    int         isRelaxed;             /* is sumt tree a relaxed clock tree ?           */
    int         isClock;               /* is sumt tree a clock tree ?                   */
    int         isCalibrated;          /* is sumt tree calibrated ?                     */
    int         nESets;                /* number of event sets                          */
    int         nBSets;                /* number of branch rate sets                    */
    char      **bSetName;              /* name of effective branch length sets          */
    char      **eSetName;              /* name of event sets                            */
    int         popSizeSet;            /* do sumt trees have population size set?       */
    char       *popSizeSetName;        /* name of population size set                   */
    int         BitsLongsNeeded;       /* number of safe longs needed for taxon bits    */
    int         runId;                 /* id of run being processed                     */
    int         numTaxa;               /* number of sumt taxa                           */
    char      **taxaNames;             /* names of sumt taxa                            */
    int        *numFileTrees;          /* number of trees per file                      */
    int        *numFileTreesSampled;   /* number of trees sampled per file              */
    int         HPD;                   /* use highest posterior density?                */
    int         saveBrParams;          /* save node/branch parameters?                  */
    MrBFlt      minBrParamFreq;        /* threshold for printing branch params to file  */
    } Sumt;

/* formats for consensus trees */
#define SIMPLE      0
#define FIGTREE     1

typedef struct comptree
    {
    char        comptFileName1[120];   /* name of first input file                      */
    char        comptFileName2[120];   /* name of second input file                     */
    char        comptOutfile[120];     /* name of output file                           */
    int         burnin;                /* actual burnin used when parsing tree files    */
    MrBFlt      minPartFreq;           /* use partitions with frequency >= minPartFreq  */
    } Comptree;

typedef struct sump
    {
    char        sumpFileName[100];     /* name of input file                            */
    char        sumpOutfile[120];      /* name of output file                           */
    //int       plot;                  /* output plot (y/n)?                            */
    int         table;                 /* output table (y/n)?                           */
    int         margLike;              /* output marginal likelihood (y/n)?             */
    int         numRuns;               /* number of independent analyses to summarize   */
    int         allRuns;               /* should data for all runs be printed (yes/no)? */
    int         HPD;                   /* use highest posterior density?                */
    MrBFlt      minProb;               /* cut-off for model probabilities to show       */
    } Sump;

typedef struct sumss
    {
    //int       plot;                  /* output plot (y/n)?                            */
    int         numRuns;               /* number of independent analyses to summarize   */
    int         allRuns;               /* should data for all runs be printed (yes/no)? */
    int         stepToPlot;            /* Which step to plot in the step plot           */
    int         askForMorePlots;       /* Should user be asked to plot for different discardfraction (y/n)?  */
    int         smoothing;             /* An integer indicating number of neighbors to average over when dooing smoothing of curvs on plots */
    MrBFlt      discardFraction;       /* Proportion of samples discarded when ploting step plot.*/
    } Sumss;

typedef struct plot
    {
    char        plotFileName[120];     /* name of input file                            */
    char        parameter[100];        /* parameter(s) to be plotted                    */
    char        match[100];            /* whether the match needs to be perfect         */
    } Plot;

typedef struct
    {
    int         numTrees;              /* number of trees to reassemble                 */
    int         numRuns;               /* number of runs to reassemble                  */
    } ReassembleInfo;

typedef struct doublet
    {
    BitsLong    first, second;
    } Doublet;

typedef struct matrix
    {
    BitsLong *origin;
    int rowSize;
    int nRows;
    int column;
    int row;
    } Matrix;

typedef struct charinfo
    {
    int dType;
    int cType;
    int nStates;
    int constant[10];
    int singleton[10];
    int variable;
    int informative;
    } CharInfo;
    
typedef struct 
    {
    int         isExcluded;            /* is the character excluded                     */
    int         numStates;             /* number of observed states for the character   */
    int         charType;              /* type of character                             */
    int         isMissAmbig;           /* is the character missing or ambiguous         */
    int         ctype;                 /* ordering of character                         */
    int         charId;                /* char ID index for doublet and codon models    */
    int         pairsId;               /* char ID for doublets                          */
    int         bigBreakAfter;         /* is there a large break after this character   */
    }
    CharInformation;

typedef struct 
    {
    int         isDeleted;             /* is the taxon deleted                          */
    int         charCount;             /* count holder                                  */
    }
    TaxaInformation;

typedef struct
    {
    MrBFlt      curScore;
    MrBFlt      minScore;
    MrBFlt      totalScore;
    MrBFlt      stopScore;
    MrBFlt      warp;
    TreeNode    **leaf;
    TreeNode    **vertex;
    }
    TreeInfo;

typedef struct
    {
    int     allavailable;
    }
    ShowmovesParams;

typedef struct
    {
    int     numNames;
    char    **names;
    }
    NameSet;

/* global variables */
extern int              abortMove;                              /* flag determining whether to abort move        */
extern int              *activeParams[NUM_LINKED];              /* a table holding the parameter status          */
extern int              *activeParts;                           /* partitions changes should apply to            */
extern int              autoClose;                              /* autoclose                                     */
extern int              autoOverwrite;                          /* Overwrite or append outputfiles when nowarnings=yes */
extern int              chainHasAdgamma;                        /* indicates if chain has adgamma HMMs           */
extern Chain            chainParams;                            /* holds parameters for Markov chain             */
extern CharInformation  *charInfo;                              /* holds critical information about characters   */
extern char             **charSetNames;                         /* holds names of character sets                 */
extern int              *compCharPos;                           /* char position in compressed matrix            */
extern int              *compColPos;                            /* column position in compressed matrix          */
extern BitsLong         *compMatrix;                            /* compressed character matrix                   */
extern int              compMatrixRowSize;                      /* row size of compressed matrix                 */
extern Comptree         comptreeParams;                         /* holds parameters for comparetree command      */
extern char             **constraintNames;                      /* holds names of constraints                    */
extern BitsLong         **definedConstraint;                    /* holds information about defined constraints   */
extern BitsLong         **definedConstraintTwo;                 /* bitfields representing second taxa sets of defined constraints (for PARTIAL constraints) */
extern BitsLong         **definedConstraintPruned;              /* bitfields representing taxa sets of defined constraints after delited taxa are removed */
extern BitsLong         **definedConstraintTwoPruned;           /* bitfields representing second taxa sets of defined constraints after delited taxa are removed (for PARTIAL constraints) */
extern int              dataType;                               /* type of data                                  */
extern Calibration      defaultCalibration;                     /* default model settings                        */
extern ModelParams      defaultModel;                           /* default model settings                        */
extern int              defChars;                               /* flag for whether number of characters is known*/
extern int              defMatrix;                              /* flag for whether matrix is successfull read   */
extern int              defPairs;                               /* flag for whether constraints on tree are read */
extern int              defPartition;                           /* flag for whether character partition is read  */
extern int              defTaxa;                                /* are taxon labels defined ?                    */
extern Doublet          doublet[16];                            /* holds information on states for doublets      */
extern int              echoMB;                                 /* flag used by Manual to prevent echoing        */
extern BitsLong         expecting;                              /* variable denoting expected token type         */
extern int              fileNameChanged;                        /* has file name been changed?                   */
extern int              foundNewLine;                           /* whether a new line has been found             */
extern char             gapId;                                  /* gap character Id                              */
extern RandLong         globalSeed;                             /* seed that is initialized at start up          */
extern char             **headerNames;                          /* string to hold headers in sump and plot       */
extern int              inComment;                              /* flag for whether input stream is commented    */
extern int              inferAncStates;                         /* should ancestral states be inferred (y/n)     */
extern int              inferSiteOmegas;                        /* should site omega values be inferred (y/n)    */
extern int              inferSiteRates;                         /* should site rates be inferred (y/n)           */
extern int              inferPosSel;                            /* should positive selection be inferred (y/n)   */
extern char             inputFileName[100];                     /* input (NEXUS) file name                       */
extern int              inTreesBlock;                           /* are we in the sumt block                      */
extern int              inValidCommand;                         /* a useful flag set whenever you enter a cmd    */
extern int              isInAmbig, isInPoly;                    /* flags whether we are within () or {}          */
extern int              isMixed;                                /* flags whether dataset is mixed                */
extern int              inMrbayesBlock;                         /* flag for whether we are in a mrbayes block    */
extern int              *intValues;                             /* integer values of parameters                  */
extern int              isTaxsetDef;                            /* is a taxon set defined                        */
extern int              isTranslateDef;                         /* is a translation block defined                */
extern int              isTranslateDiff;                        /* is translate different from current taxaset?  */
extern int              *linkTable[NUM_LINKED];                 /* how parameters are linked across parts        */
extern int              localOutGroup;                          /* outgroup for non-excluded taxa                */
extern char             **localTaxonNames;                      /* points to names of non-excluded taxa          */
extern FILE             *logFileFp;                             /* file pointer to log file                      */
extern char             logFileName[100];                       /* name of the log file                          */
extern int              logToFile;                              /* should screen output be logged to a file      */
extern char             manFileName[100];                       /* name of man file                              */
extern char             matchId;                                /* mach character Id                             */
extern int              *matrix;                                /* matrix containing original data               */
extern int              matrixHasPoly;                          /* flag for whether matrix has polymorphisms     */
extern int              memAllocs[NUM_ALLOCS];                  /* allocated memory flags                        */
extern int              mode;                                   /* mode of program (interactive/noninteractive)  */
extern char             **modelIndicatorParams;                 /* model indicator params                        */
extern char             ***modelElementNames;                   /* names for component models                    */
extern MCMCMove         **moves;                                /* vector of applicable moves                    */
extern MoveType         moveTypes[NUM_MOVE_TYPES];              /* holds information on the move types           */
extern char             missingId;                              /* missing character Id                          */
extern Tree             **mcmcTree;                             /* pointers to mcmc trees                        */
extern Model            *modelParams;                           /* holds model params for partitions             */
extern ModelInfo        *modelSettings;                         /* stores important info on model params         */
extern int              nBitsInALong;                           /* number of bits in a BitsLong                  */
extern Calibration      *nodeCalibration;                       /* holds information about node calibrations     */
extern int              noWarn;                                 /* no warnings on overwriting files              */
extern int              nPThreads;                              /* number of pthreads to use                     */
extern int              numActiveLocks;                         /* number of active, locked nodes                */
extern int              numApplicableMoves;                     /* number of moves applicable to parameters      */
extern int              numChar;                                /* number of characters in character matrix      */
extern int              numCharSets;                            /* holds number of character sets                */
extern int              numComments;                            /* number of nested comments                     */
extern int              numCompressedChars;                     /* number of compressed characters               */
extern int              numCurrentDivisions;                    /* number of partitions of data                  */
extern int              numDefinedConstraints;                  /* number of constraints defined                 */
extern enum ConstraintType  *definedConstraintsType;            /* Store type of constraint                      */
extern int              numDefinedPartitions;                   /* number of partitions defined                  */
extern int              numDefinedSpeciespartitions;            /* number of species partitions defined          */
extern int              numGlobalChains;                        /* number of global chains                       */
extern int              numLocalTaxa;                           /* number of non-excluded taxa                   */
extern int              numLocalChar;                           /* number of non-excluded characters             */
extern int              numMoveTypes;                           /* the number of move types                      */
extern int              numOpenExeFiles;                        /* number of execute files open                  */
extern int              numParams;                              /* number of parameters in model                 */
extern int              numDivisions;                           /* number of current divisions                   */
extern int              numPrintParams;                         /* number of substitution model parameters to print */
extern int              numPrintTreeParams;                     /* number of tree model parameters to print      */
extern CLFlt            *numSitesOfPat;                         /* no. sites of each pattern                     */
extern int              numSpecies;                             /* number of species in current speciespartition */
extern int              numTaxa;                                /* number of taxa in character matrix            */
extern int              numTaxaSets;                            /* holds number of taxa sets                     */
extern int              numTopologies;                          /* number of topologies for one chain and state  */
extern int              numTranslates;                          /* number of taxa in active translate block      */
extern int              numTrees;                               /* number of trees for one chain and state       */
extern int              numUserTrees;                           /* number of defined user trees                  */
extern int              *numVars;                               /* number of variables in setting arrays         */
extern int              *origChar;                              /* index from compressed char to original char   */
extern int              outGroupNum;                            /* number of outgroup taxon                      */
extern ParmInfo         paramTable[];                           /* information on parameters                     */
extern MrBFlt           *paramValues;                           /* values of parameters                          */
extern int              **partitionId;                          /* holds information about defined partitions    */
extern char             **partitionNames;                       /* hold names of partitions (first is "default") */
extern MrBFlt           *parameterValues;                       /* vector holding sump or plot parameters        */
extern Param            *params;                                /* vector of parameters in model                 */
extern int              partitionNum;                           /* index of current partition                    */
extern Plot             plotParams;                             /* holds parameters for plot command             */
extern int              precision;                              /* precision of samples and summary stats        */
extern int              *printAncStates;                        /* divisions to print anc states for             */
extern int              quitOnError;                            /* quit on error?                                */
extern int              readComment;                            /* should we read comment (looking for &)?       */
extern int              readWord;                               /* should we read a word next?                   */
extern ReassembleInfo   reassembleParams;                       /* holds parameters for reassemble command       */
extern int              replaceLogFile;                         /* should logfile be replace/appended to         */
extern RandLong         runIDSeed;                              /* seed used only for generating run ID [stamp]  */
extern BitsLong         bitsLongWithAllBitsSet;                 /* a BitsLong with all bits set, for bit ops     */
extern int              setUpAnalysisSuccess;                   /* Set to YES if analysis is set without error   */
extern int              scientific;                             /* use scientific format for samples ?           */
extern ShowmovesParams  showmovesParams;                        /* holds parameters for Showmoves command        */
extern char             spacer[10];                             /* holds blanks for printing indentations        */
extern NameSet          *speciesNameSets;                       /* hold species name sets, one for each speciespartition */
extern int              **speciespartitionId;                   /* holds info about defined speciespartitions    */
extern char             **speciespartitionNames;                /* hold names of speciespartitions (first is "default") */
extern int              speciespartitionNum;                    /* index of current species partition            */
extern char             stamp[11];                              /* holds a unique identifier for each analysis   */
extern RandLong         swapSeed;                               /* seed used only for determining which to swap  */
extern int              state[MAX_CHAINS];                      /* state of chain                                */
extern MrBFlt           *stdStateFreqs;                         /* std char state frequencies                    */
extern int              *stdType;                               /* compressed std char type: ord, unord, irrev   */
extern Sump             sumpParams;                             /* holds parameters for sump command             */
extern char             sumpToken[];                            /* string holding a .p file token                */
extern char             *sumpTokenP;                            /* pointer to a .p file token                    */
extern Sumt             sumtParams;                             /* holds parameters for sumt command             */
extern Sumss            sumssParams;                            /* holds parameters for sumss command            */
extern int              stdStateFreqsRowSize;                   /* row size for stdStateFreqs                    */
extern int              *sympiIndex;                            /* sympi state freq index for multistate chars   */
extern TaxaInformation  *taxaInfo;                              /* holds critical information about taxa         */
extern char             **taxaNames;                            /* holds name of taxa                            */
extern char             **taxaSetNames;                         /* holds names of taxa sets                      */
extern BitsLong         **taxaSet;                              /* holds information about defined taxasets      */
extern int              *tempActiveConstraints;                 /* info on the active constraints in prset       */
extern int              *tempLinkUnlink[NUM_LINKED];            /* for changing parameter linkage                */
extern int              *tempLinkUnlinkVec;                     /* for changing parameter linkage                */
extern MrBFlt           *tempNum;                               /* vector of numbers used for setting arrays     */
extern int              *tempSet;                               /* temporarily holds defined character set       */
extern int              theAmbigChar;                           /* int containing ambiguous character            */
extern int              *tiIndex;                               /* compressed std char ti index                  */
extern Calibration      *tipCalibration;                        /* holds tip calibrations                        */
extern char             **transFrom;                            /* translation block information                 */
extern char             **transTo;                              /* translation block information                 */
extern int              userBrlensDef;                          /* are the branch lengths on user tree defined   */
extern int              userLevel;                              /* the level of the user                         */
extern PolyTree         *userTree[];                            /* array of user trees                           */
extern char             workingDir[100];                        /* working directory                             */
#if defined (BEAGLE_ENABLED)
extern int              tryToUseBEAGLE;                         /* try to use the BEAGLE library                 */
extern long             beagleFlags;                            /* BEAGLE requirement flags                      */
extern int              beagleResourceNumber;                   /* BEAGLE resource number                        */
extern int*             beagleResource;                         /* BEAGLE resource list                          */
extern int              beagleResourceCount;                    /* BEAGLE resource list length                   */
extern int              beagleInstanceCount;                    /* total number of BEAGLE instances              */
extern int              beagleScalingScheme;                    /* BEAGLE dynamic scaling                        */
extern int              beagleScalingFrequency;                 /* BEAGLE rescaling frequency                    */
extern int              recalcScalers;                      /* shoud we recalculate scalers for one of divisions for current state YES/NO */
#endif
#if defined (THREADS_ENABLED)
extern int              tryToUseThreads;                        /* try to use pthreads with BEAGLE library       */
#endif

/* Aamodel parameters */
extern MrBFlt           aaJones[20][20];             /* rates for Jones model                        */
extern MrBFlt           aaDayhoff[20][20];           /* rates for Dayhoff model                      */
extern MrBFlt           aaMtrev24[20][20];           /* rates for mtrev24 model                      */
extern MrBFlt           aaMtmam[20][20];             /* rates for mtmam model                        */
extern MrBFlt           aartREV[20][20];             /* rates for rtREV model                        */
extern MrBFlt           aaWAG[20][20];               /* rates for WAG model                          */
extern MrBFlt           aacpREV[20][20];             /* rates for aacpREV model                      */
extern MrBFlt           aaVt[20][20];                /* rates for VT model                           */
extern MrBFlt           aaBlosum[20][20];            /* rates for Blosum62 model                     */
extern MrBFlt           aaLG[20][20];                /* rates for LG model                           */
extern MrBFlt           jonesPi[20];                 /* stationary frequencies for Jones model       */
extern MrBFlt           dayhoffPi[20];               /* stationary frequencies for Dayhoff model     */
extern MrBFlt           mtrev24Pi[20];               /* stationary frequencies for mtrev24 model     */
extern MrBFlt           mtmamPi[20];                 /* stationary frequencies for mtmam model       */
extern MrBFlt           rtrevPi[20];                 /* stationary frequencies for rtREV model       */
extern MrBFlt           wagPi[20];                   /* stationary frequencies for WAG model         */
extern MrBFlt           cprevPi[20];                 /* stationary frequencies for aacpREV model     */
extern MrBFlt           vtPi[20];                    /* stationary frequencies for VT model          */
extern MrBFlt           blosPi[20];                  /* stationary frequencies for Blosum62 model    */
extern MrBFlt           lgPi[20];                    /* stationary frequencies for LG model          */

#if defined (PRINT_DUMP)
FILE                    *dumpFile;                   /* for debugging logs */
#endif

#if defined (MPI_ENABLED)
extern int              proc_id;                                /* process ID (0, 1, ..., num_procs-1)                        */
extern int              num_procs;                              /* number of active processors                                */
extern MrBFlt           myStateInfo[7];                         /* likelihood/prior/heat/ran/moveInfo vals of me              */
extern MrBFlt           partnerStateInfo[7];                    /* likelihood/prior/heat/ran/moveInfo vals of partner         */
#endif

#if defined (FAST_LOG)
extern CLFlt            scalerValue[];
extern CLFlt            logValue[];
#endif

#endif  /* __BAYES_H__ */
