/* global variables */
extern int				abortMove;					            /* flag determining whether to abort move        */
extern int				*activeParams[NUM_LINKED];              /* a table holding the parameter status          */
extern int              *activeParts;                           /* partitions changes should apply to            */
extern int				autoClose;                              /* autoclose                                     */
extern int 				autoOverwrite;                          /* Overwrite or append outputfiles when nowarnings=yes */
extern int				chainHasAdgamma;						/* indicates if chain has adgamma HMMs			 */
extern Chain			chainParams;                            /* holds parameters for Markov chain             */
extern CharInformation	*charInfo;								/* holds critical information about characters   */
extern char				**charSetNames;                         /* holds names of character sets                 */
extern int				*compCharPos;		                    /* char position in compressed matrix            */
extern int				*compColPos;		                    /* column position in compressed matrix		     */
extern SafeLong		    *compMatrix;		                    /* compressed character matrix					 */
extern int				compMatrixRowSize;	                    /* row size of compressed matrix				 */
extern Comptree			comptreeParams;                         /* holds parameters for comparetree command      */
extern char				**constraintNames;					    /* holds names of constraints                    */
extern SafeLong         **definedConstraint;                    /* holds information about defined constraints   */
extern SafeLong        **definedConstraintTwo;                  /* bitfields representing second taxa sets of defined constraints (used for PARTIAL constraints)                                  */
extern SafeLong        **definedConstraintPruned;               /* bitfields representing taxa sets of defined constraints after delited taxa are removed                                         */
extern SafeLong        **definedConstraintTwoPruned;            /* bitfields representing second taxa sets of defined constraints  after delited taxa are removed(used for PARTIAL constraints)   */
extern int				dataType;                               /* type of data                                  */
extern Calibration      defaultCalibration;                     /* default model settings                        */
extern ModelParams      defaultModel;                           /* default model settings                        */
extern int				defChars;                               /* flag for whether number of characters is known*/
extern int				defMatrix;                              /* flag for whether matrix is successfull read   */
extern int				defPairs;                               /* flag for whether constraints on tree are read */
extern int				defPartition;                           /* flag for whether character partition is read  */
extern int				defTaxa;                                /* are taxon labels defined ?                    */
extern Doublet			doublet[16];                            /* holds information on states for doublets      */
extern int				echoMB;							   	    /* flag used by Manual to prevent echoing        */
extern SafeLong         expecting;								/* variable denoting expected token type         */
extern int				fileNameChanged;					    /* has file name been changed?                   */
extern int				foundNewLine;                           /* whether a new line has been found             */
extern char				gapId;                                  /* gap character Id                              */
extern SafeLong			globalSeed;                             /* seed that is initialized at start up          */
extern char				**headerNames;                          /* string to hold headers in sump and plot       */
extern int 				inComment;                              /* flag for whether input stream is commented    */
extern int				inferAncStates;					   	    /* should ancestral states be inferred (y/n)     */
extern int				inferSiteOmegas;					   	/* should site omega values be inferred (y/n)    */
extern int				inferSiteRates;					   	    /* should site rates be inferred (y/n)           */
extern int				inferPosSel;					   	    /* should positive selection be inferred (y/n)   */
extern char				inputFileName[100];                     /* input (NEXUS) file name                       */
extern int				inTreesBlock;                           /* are we in the sumt block                      */
extern int				inValidCommand;                         /* a useful flag set whenever you enter a cmd    */
extern int  			isInAmbig, isInPoly;                    /* flags whether we are within () or {}          */
extern int  			isMixed;			                    /* flags whether dataset is mixed                */
extern int				inMrbayesBlock;                         /* flag for whether we are in a mrbayes block    */
extern int				*intValues;                             /* integer values of parameters                  */
extern int				isTaxsetDef;							/* is a taxon set defined                        */
extern int				isTranslateDef;							/* is a translation block defined                */
extern int				isTranslateDiff;					    /* is translate different from current taxaset?  */
extern int              *linkTable[NUM_LINKED];                 /* how parameters are linked across parts        */
extern int				localOutGroup;				            /* outgroup for non-excluded taxa                */
extern char				**localTaxonNames;						/* points to names of non-excluded taxa          */
extern FILE				*logFileFp;                             /* file pointer to log file                      */
extern char				logFileName[100];                       /* name of the log file                          */
extern int				logToFile;                              /* should screen output be logged to a file      */
extern char				manFileName[100];						/* name of man file								 */
extern char				matchId;                                /* mach character Id                             */
extern int				*matrix;                                /* matrix containing original data               */
extern int  			matrixHasPoly;                          /* flag for whether matrix has polymorphisms     */
extern int				memAllocs[NUM_ALLOCS];                  /* allocated memory flags                        */
extern int				mode;				                    /* mode of program (interactive/noninteractive)  */
extern char             **modelIndicatorParams;                 /* model indicator params                        */
extern char             ***modelElementNames;                   /* names for component models                    */
extern MCMCMove		    **moves;								/* vector of applicable moves					 */
extern MoveType		    moveTypes[NUM_MOVE_TYPES];              /* holds information on the move types           */
extern char				missingId;                              /* missing character Id                          */
extern Tree				**mcmcTree;								/* pointers to mcmc trees						 */
extern Model			*modelParams;							/* holds model params for partitions             */
extern ModelInfo		*modelSettings;							/* stores important info on model params         */
extern int              nBitsInALong;                           /* number of bits in a SafeLong                  */
extern Calibration      *nodeCalibration;                       /* holds information about node calibrations     */
extern int				noWarn;                					/* no warnings on overwriting files              */
extern int              nPThreads;                              /* number of pthreads to use                     */
extern int				numActiveLocks;					   	    /* number of active, locked nodes                */
extern int				numApplicableMoves;						/* number of moves applicable to parameters      */
extern int				numChar;                                /* number of characters in character matrix      */
extern int				numCharSets;                            /* holds number of character sets                */
extern int				numComments;                            /* number of nested comments				     */
extern int				numCompressedChars;                     /* number of compressed characters				 */
extern int				numCurrentDivisions;                    /* number of partitions of data                  */
extern int				numDefinedConstraints;                  /* number of constraints defined                 */
extern enum ConstraintType  *definedConstraintsType;            /* Store type of constraint                      */
extern int				numDefinedPartitions;                   /* number of partitions defined                  */
extern int				numDefinedSpeciespartitions;            /* number of species partitions defined          */
extern int				numGlobalChains;						/* number of global chains						 */
extern int				numLocalTaxa;						    /* number of non-excluded taxa                   */
extern int				numLocalChar;							/* number of non-excluded characters             */
extern int				numMoveTypes;		                    /* the number of move types                      */
extern int				numOpenExeFiles;					    /* number of execute files open                  */
extern int				numParams;								/* number of parameters in model				 */
extern int				numDivisions;                           /* number of current divisions                   */
extern int				numPrintParams;						    /* number of substitution model parameters to print */
extern int				numPrintTreeParams;						/* number of tree model parameters to print      */
extern CLFlt			*numSitesOfPat;		                    /* no. sites of each pattern					 */
extern int				numSpecies;                             /* number of species in current speciespartition */
extern int				numTaxa;                                /* number of taxa in character matrix            */
extern int				numTaxaSets;                            /* holds number of taxa sets                     */
extern int				numTopologies;						    /* number of topologies for one chain and state	 */
extern int				numTranslates;                          /* number of taxa in active translate block      */
extern int				numTrees;						        /* number of trees for one chain and state	     */
extern int				numUserTrees;						    /* number of defined user trees				     */
extern int              *numVars;                               /* number of variables in setting arrays         */
extern int				*origChar;			                    /* index from compressed char to original char   */
extern int				outGroupNum;                            /* number of outgroup taxon                      */
extern ParmInfo			paramTable[];						    /* information on parameters                     */
extern MrBFlt			*paramValues;						    /* values of parameters                          */
extern int              **partitionId;                          /* holds information about defined partitions    */
extern char				**partitionNames;                       /* hold names of partitions (first is "default") */
extern MrBFlt			*parameterValues;                       /* vector holding sump or plot parameters        */
extern Param			*params;								/* vector of parameters in model				 */
extern int				partitionNum;                           /* index of current partition                    */
extern Plot				plotParams;                             /* holds parameters for plot command             */
extern int 				precision;                              /* precision of samples and summary stats        */
extern int				*printAncStates;                        /* divisions to print anc states for             */
extern int				quitOnError;							/* quit on error?					             */
extern int				readComment;							/* should we read comment (looking for &)?       */
extern int				readWord;							    /* should we read a word next?                   */
extern ReassembleInfo	reassembleParams;		                /* holds parameters for reassemble command       */
extern int				replaceLogFile;                         /* should logfile be replace/appended to         */
extern SafeLong			runIDSeed;                              /* seed used only for generating run ID [stamp]  */
extern SafeLong			safeLongWithAllBitsSet;                 /* a SafeLong with all bits set, for bit ops     */
extern int 				setUpAnalysisSuccess;                   /* Set to YES if analysis is set without error   */
extern int 				scientific;                             /* use scientific format for samples ?           */
extern ShowmovesParams	showmovesParams;					    /* holds parameters for Showmoves command        */
extern char				spacer[10];                             /* holds blanks for printing indentations        */
extern NameSet		    *speciesNameSets;                       /* hold species name sets, one for each speciespartition */
extern int              **speciespartitionId;                   /* holds info about defined speciespartitions    */
extern char				**speciespartitionNames;                /* hold names of speciespartitions (first is "default") */
extern int				speciespartitionNum;                    /* index of current species partition            */
extern char			    stamp[11];                              /* holds a unique identifier for each analysis   */
extern SafeLong			swapSeed;                               /* seed used only for determining which to swap  */
extern int				state[MAX_CHAINS];						/* state of chain								 */
extern MrBFlt			*stdStateFreqs;				            /* std char state frequencies					 */
extern int				*stdType;				                /* compressed std char type: ord, unord, irrev  */
extern Sump				sumpParams;                             /* holds parameters for sump command             */
extern char				sumpToken[];							/* string holding a .p file token                */
extern char				*sumpTokenP;							/* pointer to a .p file token					 */
extern Sumt				sumtParams;                             /* holds parameters for sumt command             */
extern Sumss			sumssParams;                            /* holds parameters for sumss command            */
extern char				stamp[11];                   			/* holds a unique identifier for each analysis   */
extern int				stdStateFreqsRowSize;					/* row size for stdStateFreqs					 */
extern int				*sympiIndex;							/* sympi state freq index for multistate chars   */
extern TaxaInformation	*taxaInfo;								/* holds critical information about taxa         */
extern char				**taxaNames;                            /* holds name of taxa                            */
extern char				**taxaSetNames;                         /* holds names of taxa sets                      */
extern SafeLong         **taxaSet;                              /* holds information about defined taxasets      */
extern int              *tempActiveConstraints;                 /* info on the active constraints in prset       */
extern int              *tempLinkUnlink[NUM_LINKED];            /* for changing parameter linkage                */
extern int              *tempLinkUnlinkVec;                     /* for changing parameter linkage                */
extern MrBFlt           *tempNum;                               /* vector of numbers used for setting arrays     */
extern int				*tempSet;                               /* temporarily holds defined character set       */
extern int  			theAmbigChar;                           /* int containing ambiguous character            */
extern int				*tiIndex;				                /* compressed std char ti index                  */
extern Calibration		*tipCalibration;                        /* holds tip calibrations                        */
extern char				**transFrom;                            /* translation block information                 */
extern char				**transTo;                              /* translation block information                 */
extern int				userBrlensDef;                          /* are the branch lengths on user tree defined   */
extern int              userLevel;                              /* the level of the user                         */ 	
extern PolyTree			*userTree[];						    /* array of user trees							 */
extern char			    workingDir[100];                        /* working directory                             */
#if defined (BEAGLE_ENABLED)
extern int              tryToUseBEAGLE;                         /* try to use the BEAGLE library                 */
extern long             beagleFlags;                            /* BEAGLE requirement flags                      */
extern int*				beagleResource;						    /* BEAGLE resource list                          */
extern int				beagleResourceCount;					/* BEAGLE resource list length                   */
extern int              beagleInstanceCount;                    /* total number of BEAGLE instances              */
extern int              beagleScalingScheme;                    /* BEAGLE dynamic scaling                        */
extern int				beagleScalingFrequency;					/* BEAGLE rescaling frequency                    */
extern int				recalcScalers;						/* shoud we recalculate scalers for one of divisions for current state YES/NO */
#endif
#if defined (THREADS_ENABLED)
extern int			    tryToUseThreads;						/* try to use pthreads with BEAGLE library       */
#endif

/* Aamodel parameters */
extern MrBFlt			aaJones[20][20];	         /* rates for Jones model                        */
extern MrBFlt			aaDayhoff[20][20];           /* rates for Dayhoff model                      */
extern MrBFlt			aaMtrev24[20][20];	         /* rates for mtrev24 model                      */
extern MrBFlt			aaMtmam[20][20];	         /* rates for mtmam model                        */
extern MrBFlt			aartREV[20][20];             /* rates for rtREV model                        */
extern MrBFlt			aaWAG[20][20];               /* rates for WAG model                          */
extern MrBFlt			aacpREV[20][20];             /* rates for aacpREV model                      */
extern MrBFlt			aaVt[20][20];                /* rates for VT model                           */
extern MrBFlt			aaBlosum[20][20];            /* rates for Blosum62 model                     */
extern MrBFlt			jonesPi[20];                 /* stationary frequencies for Jones model       */
extern MrBFlt			dayhoffPi[20];               /* stationary frequencies for Dayhoff model     */
extern MrBFlt			mtrev24Pi[20];               /* stationary frequencies for mtrev24 model     */
extern MrBFlt			mtmamPi[20];                 /* stationary frequencies for mtmam model       */
extern MrBFlt			rtrevPi[20];                 /* stationary frequencies for rtREV model       */
extern MrBFlt			wagPi[20];                   /* stationary frequencies for WAG model         */
extern MrBFlt			cprevPi[20];                 /* stationary frequencies for aacpREV model     */
extern MrBFlt			vtPi[20];                    /* stationary frequencies for VT model          */
extern MrBFlt			blosPi[20];                  /* stationary frequencies for Blosum62 model    */

#if defined (PRINT_RATEMULTIPLIERS_CPP)
extern FILE             *rateMultfp;
#endif

#if defined (MPI_ENABLED)
extern int 				proc_id;                                /* process ID (0, 1, ..., num_procs-1)                        */
extern int 				num_procs;                              /* number of active processors                                */
extern MrBFlt			myStateInfo[7];                         /* likelihood/prior/heat/ran/moveInfo vals of me              */
extern MrBFlt			partnerStateInfo[7];                    /* likelihood/prior/heat/ran/moveInfo vals of partner         */
#endif

#if defined (FAST_LOG)
extern CLFlt			scalerValue[];
extern CLFlt			logValue[];
#endif
