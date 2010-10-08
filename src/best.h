#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <stdarg.h>
#include "mb.h"

/*#define DEBUG*/
#define LSPNAME  30
#define ERROR 1
#define NO_ERROR 0
#define YES 1
#define NO 0
#define NA -1

#define NMOVE 4
#define MOVETREE 1
#define MOVETHETA 2
#define MOVESR  3
#define MOVEER  4

/* parameter ID values */
#define P_TOP   0
#define P_POP   1
#define P_SR 	2
#define P_ER	3
#define NSPECIES 200
#define NGENE 1010

#define FPN(file) fputc('\n', file)
#define F0 stdout
#define FOR(i,n) for(i=0; i<n; i++)
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGamma(alpha))





typedef struct spmodel
	{
       double thetaprior[2];
	   int thetainvgamma;
	int sRprior[2];
	int eRprior[2];
        double sF;
	double treeHeightExp;
	}
	ModelParam;

typedef struct MCMCPARAMETERS 
	{
	int nchain;
   	int nsptree;
       	char chainFileName[100];
	}  
	McmcPara; 

typedef struct PARAM
	{
	int			paramType;		/* the type of the parameter					*/
	double			*values;		/* main values of parameter						*/
	int			nSubParams;		/* number of subparams							*/
	}
	SParam; 

typedef struct COALPOP 
	{
   	int nin[NSPECIES], ncoal[NSPECIES], nout[NSPECIES],nodes[NSPECIES];
   	double tj[NSPECIES][NSPECIES];
	}  
	CoalTime; 
/* typedef for a MoveFxn */
typedef int (SPMoveFxn)(SParam *param, int chain, long int *seed, double *lnLikeRatio, double *lnPriorRatio, double *lnProposalRatio, double *mvp);

typedef struct
	{
	SPMoveFxn		*moveFxn;			/* pointer to the move function					*/
	int			nApplicable;		/* number of relevant params					*/
	int			applicableTo[40];	/* pointer to ID of relevant params				*/
	char		*name;				/* name of the move type						*/
	char		shortName[10];		/* abbreviated name of the move type            */
	double		relProposalProb;	/* this holds the set proposal probability      */
	double		proposalParam[2];	/* parameters for the proposal mechanism        */
	double		cumProposalProb;
       int		nparam;
	SParam        *para;
	int         parsimonyBased;     /* this move is based on parsimony (YES/NO)     */
	} SPMoveType;
/* max number of move types */

typedef struct
	{
	char		*name;				/* pointer to the name of the move type         */
	char		*shortName;	        /* pointer to the short name of the move        */
	MoveFxn		*moveFxn;			/* pointer to the move function					*/
	SParam		*parm;				/* ptr to parameter the move applies to			*/
	double		relProposalProb;	/* the actual proposal probability used			*/
	double		cumProposalProb;	/* the cumulative proposal probability			*/
	int			*nAccepted;			/* number of accepted moves						*/
	int			*nTried;			/* number of tried moves						*/
	double		proposalParam[2];	/* parameters for the proposal mechanism        */
	} SPMCMCMove;

/* tool functions*/
FILE *gfopen(char *filename, char *mode);
void SetSeed (unsigned int seed);
double Lngamma (double x);
double rndu (void);

