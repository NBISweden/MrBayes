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


/**************** typedefs used by Fredrik's code, derived from BEST code *******************/

/* struct for constraints (similar to Distance in BEST code, except
   that a bitfield is used to hold info on the taxon pair)           */
typedef struct {
	double      depth;
    SafeLong*   pairSet;
} Depth;


/**************** Declaration of functions that are called from MrBayes **************/
void    AllocateBestChainVariables(void);
int     FillSpeciesTreeParams (SafeLong* seed, int from, int to);
void    FreeBestChainVariables(void);
int     IsSpeciesTreeConsistent (Tree *speciesTree, int chain);
double  LnSpeciesTreeProb(int chain);
double  LnJointGeneTreeSpeciesTreePr(Tree **geneTrees, int numGeneTrees, Tree *speciesTree, int chain);
int     Move_GeneTree1 (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_GeneTree2 (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_GeneTree3 (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_NodeSliderGeneTree (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
int     Move_SpeciesTree (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp);
void    ShowUpperTriangMatrix(double *values, int squareSize);

/* NOTE: To add and set up more move functions, a struct needs to be added to SetUpMoveTypes in model.c */
