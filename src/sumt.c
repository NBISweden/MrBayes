/*
 *  MrBayes 3.1.2
 *
 *  copyright 2002-2005
 *
 *  John P. Huelsenbeck
 *  Section of Ecology, Behavior and Evolution
 *  Division of Biological Sciences
 *  University of California, San Diego
 *  La Jolla, CA 92093-0116
 *
 *  johnh@biomail.ucsd.edu
 *
 *	Fredrik Ronquist
 *  Paul van der Mark
 *  School of Computational Science
 *  Florida State University
 *  Tallahassee, FL 32306-4120
 *
 *  ronquist@scs.fsu.edu
 *  paulvdm@scs.fsu.edu
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

/* id-string for ident, do not edit: cvs will update this string */
const char sumtID[]="$Id: sumt.c,v 3.49 2009/02/04 13:19:49 ronquist Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "mb.h"
#include "globals.h"
#include "command.h"
#include "bayes.h"
#include "mbmath.h"
#include "sumt.h"
#include "sump.h"
#include "mcmc.h"
#include "model.h"
#include "utils.h"
#if defined(__MWERKS__)
#include "SIOUX.h"
#endif

typedef struct stnode
	{
	struct stnode	*left, *right, *anc;
	int				memoryIndex, index, upDateCl, upDateTi, marked, x, y,
					scalerNode, taxonName;
	safeLong        scalersSet, clSpace, tiSpace;
	char			label[100];
	MrBFlt			length, age;
	}
	SumtNode;


#define	MAX_PARTITIONS					10000
#define	MAX_TREES						1000

#undef	DEBUG_CONTREE

/* local prototypes */
int      AddTreeToList (int whichList);
int      AllocBits (int n);
void     AssignIntPart (safeLong *lft, safeLong *rht, safeLong *p);
void     AssignTipPart (int n, safeLong *p);
int		 BrlensVals (int treeNo, char *s, int longestLineLength, int lastTreeBlockBegin, int lastTreeBlockEnd);
void     CalculateTreeToTreeDistance (int *lst[2], MrBFlt *lngs[2], int nnds, int hasBrlens, MrBFlt *d1, MrBFlt *d2, MrBFlt *d3);
int      CheckSumtSpecies (void);
int      ConTree (void);
int      DerootSumtTree (SumtNode *p, int n, int outGrp);
int      FindParts (int nodeDepthConTree);
int      FindTree (void);
void     FinishSumtTree (SumtNode *p, int *i, int isThisTreeRooted);
int      FirstTaxonInPartition (safeLong *partition, int length);
int      FreeBits (void);
int		 GetClockRates (char *rateName, char *fileName);
void     GetConDownPass (PolyNode **downPass, PolyNode *p, int *i);
int      GetPartitions (void);
void     GetSumtDownPass (SumtNode *p, SumtNode **dp, int *i);
void     GetSumtToken (int *tokenType);
int      Label (PolyNode *p, int addIndex, char *label, int maxLength);
int		 OpenBrlensFile (int treeNo);
int      OpenComptFiles (void);
int      OpenSumtFiles (int treeNo);
int      PartFinder (safeLong *p, MrBFlt bl, int *partID);
int		 PrintBrlensToFile (void);
void     PrintParts (FILE *fp, safeLong *p, int nTaxaToShow);
int      PruneSumt (void);
int      ReallocateBits (void);
int      ReallocateFullCompTrees (int whichList);
int      ReallocateFullTrees (void);
int      ReorderParts (void);
int      RootSumtTree (SumtNode *p, int n, int out);
void     ShowBits (safeLong *p, int nBitsToShow);
void	 ShowConNodes (int nNodes, PolyNode *root);
int		 ShowConPhylogram (FILE *fp, int nNodes, PolyNode *root, int screenWidth, int isCalibrated);
void     ShowParts (FILE *fp, safeLong *p, int nTaxaToShow);
void     ShowSumtNodes (SumtNode *p, int indent, int isThisTreeRooted);
void     SortIndParts (int *item, MrBFlt *assoc, int count, int descendingOrder);
void     SortIndParts2 (int *item, MrBFlt *assoc, int left, int right, int descendingOrder);
int      SortParts (int *item, int count);
void     SortParts2 (int *item, int left, int right);
int      SumtDex (SumtNode *p);
int      TreeProb (void);
void     WriteConTree (PolyNode *p, FILE *fp, int showSupport);
void     WriteTree (PolyNode *p, FILE *fp);

extern int DoUserTree (void);
extern int DoUserTreeParm (char *parmName, char *tkn);
extern int SafeFclose(FILE **);

extern int inSumtCommand;
extern int inComparetreeCommand;

/* local (to this file) */
char		*sumtTokenP, sumtToken[CMD_STRING_LENGTH];
int			taxonLongsNeeded, numPartsAllocated, *numFoundOfThisPart, *numFoundOfThisPart1, *numFoundOfThisPart2, numTreesInLastBlock,
			**numFoundInRunOfPart, numTreePartsFound, numSumtTrees, numSumTreesSampled, numTranslates,
			whichTranslate, sumtBrlensDef, nextAvailableSumtNode, numSumtTaxa,
			isFirstSumtNode, foundSumtColon, *sumTaxaFound, numIncludedTaxa,
			isSumtTreeDefined, isSumtTreeRooted, *partOrigOrder, numAsterices,
			numTreeParts, *treePartNums, *fullTreePartIds, *numOfThisFullTree, numFullTreesAllocated, numFullTreesFound,
			*prunedTaxa, *absentTaxa, *firstPrunedTaxa, *firstAbsentTaxa, comparingFiles, fileNum, numCompTrees[2], numCompTreesSampled[2],
			numFullCompTreesFound[2], numFullCompTreesAllocated[2], *fullCompTreePartIds1, *fullCompTreePartIds2,
			numBrlens, printingBrlens, runIndex;
safeLong	*treeBits, *treePartsFound, *taxonMask;
MrBFlt		*aBrlens, *sBrlens, *treePartLengths, *fullCompTreePartLengths1, *fullCompTreePartLengths2,
			*brlens, *aWithinBrlens, *sWithinBrlens, *sumB, *sumsqB;
SumtNode	*pSumtPtr, *qSumtPtr, *sumtRoot, *sumtNodes;
PolyNode	*conNodes, *conRoot;
FILE		*fpParts, *fpCon, *fpTrees, *fpCompParts, *fpCompDists, *fpBrlens;
int			nodeDepthConTree = NO;
int			treeIndex;
int			numClockRates;
MrBFlt		*clockRate;



int AddTreeToList (int whichList)

{

	int			i, *x;
	MrBFlt		*y;
	
	if (numTreeParts == 0)
		{
		MrBayesPrint ("%s   Too few tree partitions\n", spacer);
		return (ERROR);
		}
		
	if (numFullCompTreesFound[whichList] + 1 > numFullCompTreesAllocated[whichList])
		{
		numFullCompTreesAllocated[whichList] += 500;
		if (ReallocateFullCompTrees (whichList) == ERROR)
			return (ERROR);
		}
		
	if (whichList == 0)
		{
		x = &fullCompTreePartIds1[numFullCompTreesFound[0] * 2 * numTaxa];
		for (i=0; i<numTreeParts; i++)
			x[i] = treePartNums[i];
		y = &fullCompTreePartLengths1[numFullCompTreesFound[0] * 2 * numTaxa];
		for (i=0; i<numTreeParts; i++)
			y[i] = treePartLengths[i];
		numFullCompTreesFound[0]++;
		}
	else
		{
		x = &fullCompTreePartIds2[numFullCompTreesFound[1] * 2 * numTaxa];
		for (i=0; i<numTreeParts; i++)
			x[i] = treePartNums[i];
		y = &fullCompTreePartLengths2[numFullCompTreesFound[1] * 2 * numTaxa];
		for (i=0; i<numTreeParts; i++)
			y[i] = treePartLengths[i];
		numFullCompTreesFound[1]++;
		}

	return (NO_ERROR);
	
}





int AllocBits (int n)

{

	int					i, j, offSet;
	safeLong				x, y;
	
	/* decide how many unsigned ints are going to be needed to
	   represent a taxon number */
	taxonLongsNeeded = (n / (sizeof(safeLong)*8)) + 1;
	
	/* how many partitions have been allocated */
	numPartsAllocated = MAX_PARTITIONS;
	
	/* how many trees (partition number information) have been allocated */
	numFullTreesAllocated = MAX_TREES;
	numFullCompTreesAllocated[0] = MAX_TREES;
	numFullCompTreesAllocated[1] = MAX_TREES;
	
	/* allocate memory */
	if (memAllocs[ALLOC_TREEBITS] == YES)
		{
		MrBayesPrint ("%s   treeBits not free in AllocBits\n", spacer);
		goto errorExit;
		}
	treeBits = NULL;
	treePartNums = NULL;
	treePartLengths = NULL;
	treeBits = (safeLong *)SafeMalloc((size_t) (2 * n * taxonLongsNeeded * sizeof(safeLong)));
	treePartNums = (int *)SafeMalloc((size_t) (2 * n * sizeof(int)));
	if (!treeBits || !treePartNums)
		{
		MrBayesPrint ("%s   Problem allocating treeBits (%d)\n", spacer, 2 * n * taxonLongsNeeded * sizeof(safeLong));
		goto errorExit;
		}
	treePartLengths = (MrBFlt *)SafeMalloc((size_t) (2 * n * sizeof(MrBFlt)));
	if (!treePartLengths)
		{
		MrBayesPrint ("%s   Problem allocating treePartLengths (%d)\n", spacer, 2 * n * sizeof(MrBFlt));
		goto errorExit;
		}
	memAllocs[ALLOC_TREEBITS] = YES;
	
	if (memAllocs[ALLOC_FULLTREEINFO] == YES)
		{
		MrBayesPrint ("%s   fullTreePartIds not free in AllocBits\n", spacer);
		goto errorExit;
		}
	fullTreePartIds = numOfThisFullTree = NULL;
	fullTreePartIds = (int *)SafeMalloc((size_t) (2 * numTaxa * MAX_TREES * sizeof(int)));
	if (!fullTreePartIds)
		{
		MrBayesPrint ("%s   Problem allocating fullTreePartIds (%d)\n", spacer, 2 * n * MAX_TREES * sizeof(int));
		goto errorExit;
		}
	numOfThisFullTree = (int *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(int)));
	if (!numOfThisFullTree)
		{
		MrBayesPrint ("%s   Problem allocating numOfThisFullTree (%d)\n", spacer, MAX_TREES * sizeof(int));
		goto errorExit;
		}
	memAllocs[ALLOC_FULLTREEINFO] = YES;

	if (memAllocs[ALLOC_TREEPARTS] == YES)
		{
		MrBayesPrint ("%s   treePartsFound not free in AllocBits\n", spacer);
		goto errorExit;
		}
	treePartsFound = NULL;
	treePartsFound = (safeLong *)SafeMalloc((size_t) (taxonLongsNeeded * MAX_PARTITIONS * sizeof(safeLong)));
	if (!treePartsFound)
		{
		MrBayesPrint ("%s   Problem allocating treePartsFound (%d)\n", spacer, taxonLongsNeeded * MAX_PARTITIONS * sizeof(safeLong));
		goto errorExit;
		}
	memAllocs[ALLOC_TREEPARTS] = YES;

	if (memAllocs[ALLOC_NUMOFPART] == YES)
		{
		MrBayesPrint ("%s   numFoundOfThisPart not free in AllocBits\n", spacer);
		goto errorExit;
		}
	numFoundOfThisPart = NULL;
	numFoundOfThisPart = (int *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(int)));
	if (!numFoundOfThisPart)
		{
		MrBayesPrint ("%s   Problem allocating numFoundOfThisPart (%d)\n", spacer, MAX_PARTITIONS * sizeof(int));
		goto errorExit;
		}

	if (sumtParams.numRuns > 1)
		{
		if (memAllocs[ALLOC_NUMINRUNOFPART] == YES)
			{
			MrBayesPrint ("%s   numFoundOfThisPart not free in AllocBits\n", spacer);
			goto errorExit;
			}
		numFoundInRunOfPart = (int **) calloc ((size_t)(sumtParams.numRuns), sizeof (int *));
		if (!numFoundInRunOfPart)
			{
			MrBayesPrint ("%s   Problem allocating numFoundInRunOfPart (%d bytes)\n", spacer, sumtParams.numRuns * sizeof (int *));
			goto errorExit;
			}
		for (i=0; i<sumtParams.numRuns; i++)
			{
			numFoundInRunOfPart[i] = (int *) SafeMalloc ((size_t) (MAX_PARTITIONS * sizeof(int)));
			if (!numFoundInRunOfPart[i])
				{
				MrBayesPrint ("%s   Problem allocating numFoundInRunOfPart[%d] (%d bytes)\n", spacer, i, MAX_PARTITIONS * sizeof (int *));
				goto errorExit;
				}
			}
		memAllocs[ALLOC_NUMINRUNOFPART] = YES;
		}

	if (comparingFiles == YES)
		{
		numFoundOfThisPart1 = numFoundOfThisPart2 = NULL;
		numFoundOfThisPart1 = (int *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(int)));
		if (!numFoundOfThisPart1)
			{
			MrBayesPrint ("%s   Problem allocating numFoundOfThisPart1 (%d)\n", spacer, MAX_PARTITIONS * sizeof(int));
			goto errorExit;
			}
		numFoundOfThisPart2 = (int *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(int)));
		if (!numFoundOfThisPart2)
			{
			MrBayesPrint ("%s   Problem allocating numFoundOfThisPart2 (%d)\n", spacer, MAX_PARTITIONS * sizeof(int));
			goto errorExit;
			}
		}
	memAllocs[ALLOC_NUMOFPART] = YES;
	
	if (comparingFiles == YES)
		{
		if (memAllocs[ALLOC_FULLCOMPTREEINFO] == YES)
			{
			MrBayesPrint ("%s   fullTreePartIds1 not free in AllocBits\n", spacer);
			goto errorExit;
			}
		fullCompTreePartIds1 = fullCompTreePartIds2 = NULL;
		fullCompTreePartLengths1 = fullCompTreePartLengths2 = NULL;
		fullCompTreePartIds1 = (int *)SafeMalloc((size_t) (2 * numTaxa * MAX_TREES * sizeof(int)));
		if (!fullCompTreePartIds1)
			{
			MrBayesPrint ("%s   Problem allocating fullCompTreePartIds1 (%d)\n", spacer, 2 * n * MAX_TREES * sizeof(int));
			goto errorExit;
			}
		fullCompTreePartIds2 = (int *)SafeMalloc((size_t) (2 * numTaxa * MAX_TREES * sizeof(int)));
		if (!fullCompTreePartIds2)
			{
			MrBayesPrint ("%s   Problem allocating fullCompTreePartIds2 (%d)\n", spacer, 2 * n * MAX_TREES * sizeof(int));
			goto errorExit;
			}
		fullCompTreePartLengths1 = (MrBFlt *)SafeMalloc((size_t) (2 * numTaxa * MAX_TREES * sizeof(MrBFlt)));
		if (!fullCompTreePartLengths1)
			{
			MrBayesPrint ("%s   Problem allocating fullCompTreePartLengths1 (%d)\n", spacer, 2 * n * MAX_TREES * sizeof(MrBFlt));
			goto errorExit;
			}
		fullCompTreePartLengths2 = (MrBFlt *)SafeMalloc((size_t) (2 * numTaxa * MAX_TREES * sizeof(MrBFlt)));
		if (!fullCompTreePartLengths2)
			{
			MrBayesPrint ("%s   Problem allocating fullCompTreePartLengths2 (%d)\n", spacer, 2 * n * MAX_TREES * sizeof(MrBFlt));
			goto errorExit;
			}
		memAllocs[ALLOC_FULLCOMPTREEINFO] = YES;
		}	

	if (sumtBrlensDef == YES)
		{
		if (memAllocs[ALLOC_ABRLENS] == YES)
			{
			MrBayesPrint ("%s   aBrlens not free in AllocBits\n", spacer);
			goto errorExit;
			}
		aBrlens = NULL;
		aBrlens = (MrBFlt *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(MrBFlt)));
		if (!aBrlens)
			{
			MrBayesPrint ("%s   Problem allocating aBrlens (%d)\n", spacer, MAX_PARTITIONS * sizeof(MrBFlt));
			goto errorExit;
			}
		memAllocs[ALLOC_ABRLENS] = YES;
		
		if (memAllocs[ALLOC_SBRLENS] == YES)
			{
			MrBayesPrint ("%s   sBrlens not free in AllocBits\n", spacer);
			goto errorExit;
			}
		sBrlens = NULL;
		sBrlens = (MrBFlt *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(MrBFlt)));
		if (!sBrlens)
			{
			MrBayesPrint ("%s   Problem allocating sBrlens (%d)\n", spacer, MAX_PARTITIONS * sizeof(MrBFlt));
			goto errorExit;
			}
		memAllocs[ALLOC_SBRLENS] = YES;
		}

	if (sumtParams.numRuns > 1 && sumtBrlensDef == YES)
		{
		if (memAllocs[ALLOC_A_WITHIN_BRLENS] == YES)
			{
			MrBayesPrint ("%s   aWithinBrlens not free in AllocBits\n", spacer);
			goto errorExit;
			}
		aWithinBrlens = (MrBFlt *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(MrBFlt)));
		if (!aWithinBrlens)
			{
			MrBayesPrint ("%s   Problem allocating aWithinBrlens (%d)\n", spacer, MAX_PARTITIONS * sizeof(MrBFlt));
			goto errorExit;
			}
		memAllocs[ALLOC_A_WITHIN_BRLENS] = YES;
		
		if (memAllocs[ALLOC_S_WITHIN_BRLENS] == YES)
			{
			MrBayesPrint ("%s   sWithinBrlens not free in AllocBits\n", spacer);
			goto errorExit;
			}
		sWithinBrlens = (MrBFlt *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(MrBFlt)));
		if (!sWithinBrlens)
			{
			MrBayesPrint ("%s   Problem allocating sWithinBrlens (%d)\n", spacer, MAX_PARTITIONS * sizeof(MrBFlt));
			goto errorExit;
			}
		memAllocs[ALLOC_S_WITHIN_BRLENS] = YES;
		
		if (memAllocs[ALLOC_SUMB] == YES)
			{
			MrBayesPrint ("%s   sumB not free in AllocBits\n", spacer);
			goto errorExit;
			}
		sumB = (MrBFlt *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(MrBFlt)));
		if (!sumB)
			{
			MrBayesPrint ("%s   Problem allocating sumB (%d bytes)\n", spacer, MAX_PARTITIONS * sizeof(MrBFlt));
			goto errorExit;
			}
		memAllocs[ALLOC_SUMB] = YES;

		if (memAllocs[ALLOC_SUMSQB] == YES)
			{
			MrBayesPrint ("%s   sumsqB not free in AllocBits\n", spacer);
			goto errorExit;
			}
		sumsqB = (MrBFlt *)SafeMalloc((size_t) (MAX_PARTITIONS * sizeof(MrBFlt)));
		if (!sumsqB)
			{
			MrBayesPrint ("%s   Problem allocating sumsqB (%d bytes)\n", spacer, MAX_PARTITIONS * sizeof(MrBFlt));
			goto errorExit;
			}
		memAllocs[ALLOC_SUMSQB] = YES;
		}

	if (memAllocs[ALLOC_TAXAFOUND] == YES)
		{
		MrBayesPrint ("%s   sumTaxaFound not free in AllocBits\n", spacer);
		goto errorExit;
		}
	sumTaxaFound = NULL;
	sumTaxaFound = (int *)SafeMalloc((size_t) (numTaxa * sizeof(int)));
	if (!sumTaxaFound)
		{
		MrBayesPrint ("%s   Problem allocating sumTaxaFound (%d)\n", spacer, numTaxa * sizeof(int));
		goto errorExit;
		}
	for (i=0; i<numTaxa; i++)
		sumTaxaFound[i] = NO;
	memAllocs[ALLOC_TAXAFOUND] = YES;

	if (memAllocs[ALLOC_TAXONMASK] == YES)
		{
		MrBayesPrint ("%s   taxonMask not free in AllocBits\n", spacer);
		goto errorExit;
		}
	taxonMask = NULL;
	taxonMask = (safeLong *)SafeMalloc((size_t) (taxonLongsNeeded * sizeof(safeLong)));
	if (!taxonMask)
		{
		MrBayesPrint ("%s   Problem allocating taxonMask (%d)\n", spacer, taxonLongsNeeded * sizeof(safeLong));
		goto errorExit;
		}
	memAllocs[ALLOC_TAXONMASK] = YES;
	
	if (memAllocs[ALLOC_PRUNEINFO] == YES)
		{
		MrBayesPrint ("%s   prunedTaxa not free in AllocBits\n", spacer);
		goto errorExit;
		}
	prunedTaxa = NULL;
	if (sumtParams.numRuns > 1)
		i = 4 * numTaxa;
	else
		i = 2 * numTaxa;
	prunedTaxa = (int *)SafeMalloc((size_t) (i * sizeof(int)));
	if (!prunedTaxa)
		{
		MrBayesPrint ("%s   Problem allocating prunedTaxa (%d)\n", spacer, i * sizeof(safeLong));
		goto errorExit;
		}
	absentTaxa = prunedTaxa + numTaxa;
	if (sumtParams.numRuns > 1)
		{
		firstPrunedTaxa = prunedTaxa + 2 * numTaxa;
		firstAbsentTaxa = prunedTaxa + 3 * numTaxa;
		}
	memAllocs[ALLOC_PRUNEINFO] = YES;
		
	if (memAllocs[ALLOC_SUMTTREE] == YES)
		{
		MrBayesPrint ("%s   sumtNodes not free in AllocBits\n", spacer);
		goto errorExit;
		}
	sumtNodes = NULL;
	sumtNodes = (SumtNode *)SafeMalloc((size_t) (2 * numTaxa * sizeof(SumtNode)));
	if (!sumtNodes)
		{
		MrBayesPrint ("%s   Problem allocating sumtNodes (%d)\n", spacer, 2 * numTaxa * sizeof(SumtNode));
		goto errorExit;
		}
	for (i=0; i<2*numTaxa; i++)
		{
		sumtNodes[i].left = sumtNodes[i].right = sumtNodes[i].anc = NULL;
		sumtNodes[i].memoryIndex = i;
		sumtNodes[i].length = 0.0;
		sumtNodes[i].marked = NO;
		sumtNodes[i].index = 0;
		}
	sumtRoot = NULL;
	isSumtTreeDefined = NO;
	memAllocs[ALLOC_SUMTTREE] = YES;

	/* initialize */
	for (i=0; i<2*n*taxonLongsNeeded; i++)
		treeBits[i] = 0;
	for (i=0; i<taxonLongsNeeded*MAX_PARTITIONS; i++)
		treePartsFound[i] = 0;
	for (i=0; i<MAX_PARTITIONS; i++)
		numFoundOfThisPart[i] = 0;
	if (sumtParams.numRuns > 1)
		{
		for (i=0; i<sumtParams.numRuns; i++)
			for (j=0; j<MAX_PARTITIONS; j++)
				numFoundInRunOfPart[i][j] = 0;
		}
	
	if (comparingFiles == YES)
		{
		for (i=0; i<MAX_PARTITIONS; i++)
			numFoundOfThisPart1[i] = numFoundOfThisPart2[i] = 0;
		}
	
	if (sumtBrlensDef == YES)
		{
		for (i=0; i<MAX_PARTITIONS; i++)
			aBrlens[i] = sBrlens[i] = 0.0;
		if (sumtParams.numRuns > 1)
			{
			for (i=0; i<MAX_PARTITIONS; i++)
				aWithinBrlens[i] = sWithinBrlens[i] = sumB[i] = sumsqB[i] = 0.0;
			}
		}
	numTreePartsFound = 0;
	numFullTreesFound = 0;
	for (i=0; i<taxonLongsNeeded; i++)
		taxonMask[i] = 0;

	for (i=0; i<n; i++)
		{
		offSet = 0;
		while ((i+1) > nBitsInALong*(offSet+1))
			offSet++;
		x = 1 << (i - offSet*nBitsInALong);
		y = taxonMask[offSet];
		taxonMask[offSet] = x | y;
		}

	for (i=0; i<numTaxa; i++)
		absentTaxa[i] = prunedTaxa[i] = NO;
	if (sumtParams.numRuns > 1)
		{
		for (i=0; i<numTaxa; i++)
			firstAbsentTaxa[i] = firstPrunedTaxa[i] = NO;
		}

	return (NO_ERROR);

	errorExit:
		FreeBits();
		return (ERROR);
	
}






void AssignTipPart (int n, safeLong *p)

{

	safeLong		x;

	p += n / nBitsInALong;
	x = 1 << (n % nBitsInALong);
	(*p) ^= x;
	
}





void AssignIntPart (safeLong *lft, safeLong *rht, safeLong *p)

{

	int			i;
	
	for (i=0; i<taxonLongsNeeded; i++)
		{
		(*p) = ((*lft) | (*rht));
		p++;
		lft++;
		rht++;
		}
	
}





int BrlensVals (int treeNo, char *s, int longestLineLength, int lastTreeBlockBegin, int lastTreeBlockEnd)

{
	int		i, runNo, foundBegin, inTreeBlock, inSumtComment, lineNum, tokenType;
	FILE	*fp;
	char	fileName[100];

	/* Calculate number of branch lengths to print to file */
	for (i=0; (MrBFlt)numFoundOfThisPart[i]/(sumtParams.numRuns * numSumTreesSampled) >= sumtParams.brlensFreqDisplay; i++)
		;
	numBrlens = i;
	if (numBrlens < 1)
		{
		MrBayesPrint ("%s   No branch lengths above the frequency to display (%lf).", spacer, sumtParams.brlensFreqDisplay);
		return ERROR;
		}

	/* Open brlens file */
	if (OpenBrlensFile(treeNo) == ERROR)
		return ERROR;

	/* allocate space for brlens */
	brlens = (MrBFlt *) SafeMalloc (numBrlens * sizeof(MrBFlt));
	if (brlens == NULL)
		{
		SafeFclose (&fpBrlens);
		return ERROR;
		}
	for (i=0; i<numBrlens; i++)
		brlens[i] = -1.0;
		
	/* Change global setting so that DoTree does the right thing */
	printingBrlens = YES;

	for (runNo = 0; runNo < sumtParams.numRuns; runNo++)
		{
		/* Open tree file */
		if (sumtParams.numRuns == 1 && sumtParams.numTrees == 1)
			sprintf (fileName, "%s.t", sumtParams.sumtFileName);
		else if (sumtParams.numRuns > 1 && sumtParams.numTrees == 1)
			sprintf (fileName, "%s.run%d.t", sumtParams.sumtFileName, runNo+1);
		else if (sumtParams.numRuns == 1 && sumtParams.numTrees > 1)
			sprintf (fileName, "%s.tree%d.t", sumtParams.sumtFileName, treeNo+1);
		else if (sumtParams.numRuns > 1 && sumtParams.numTrees > 1)
			sprintf (fileName, "%s.tree%d.run%d.t", sumtParams.sumtFileName, treeNo+1, runNo+1);
		
		/* open binary file */
		if ((fp = OpenBinaryFileR (fileName)) == NULL)
			{
			SafeFclose (&fpBrlens);
			free (brlens);
			return ERROR;
			}
		
		/* find length of longest line */
		longestLineLength = LongestLine (fp);
		longestLineLength += 10;
	
		/* allocate a string long enough to hold a line */
		if (runNo == 0)
			s = (char *) SafeMalloc (sizeof (char) * longestLineLength);
		else
			{
			free (s);
			s = (char *) SafeMalloc (sizeof (char) * longestLineLength);
			}

		if (!s)
			{
			free (brlens);
			SafeFclose (&fpBrlens);
			return ERROR;
			}

		/* close binary file */
		SafeFclose (&fp);
	
		/* open text file */
		if ((fp = OpenTextFileR (fileName)) == NULL)
			{
			SafeFclose (&fpBrlens);
			free (brlens);
			free (s);
			return ERROR;
			}

		/* Check file for appropriate blocks. We want to find the last tree block
			in the file and start from there. */
		foundBegin = inTreeBlock = inSumtComment = NO;
		lineNum = lastTreeBlockBegin = lastTreeBlockEnd = 0;
		while (fgets (s, longestLineLength, fp) != NULL)
			{
			sumtTokenP = &s[0];
			do
				{
				GetSumtToken (&tokenType);
				if (IsSame("[", sumtToken) == SAME)
					inSumtComment = YES;
				if (IsSame("]", sumtToken) == SAME)
					inSumtComment = NO;
					
				if (inSumtComment == NO)
					{
					if (foundBegin == YES)
						{
						if (IsSame("Trees", sumtToken) == SAME)
							{
							inTreeBlock = YES;
							foundBegin = NO;
							lastTreeBlockBegin = lineNum;
							}
						}
					else
						{
						if (IsSame("Begin", sumtToken) == SAME)
							{
							foundBegin = YES;
							}
						else if (IsSame("End", sumtToken) == SAME)
							{
							if (inTreeBlock == YES)
								{
								inTreeBlock = NO;
								lastTreeBlockEnd = lineNum;
								}
							}
						}
					}
					
				} while (*sumtToken);
			lineNum++;
			}
				
		/* Now fast rewind tree file */
		(void)fseek(fp, 0L, 0);	
	
		/* ...and fast forward to beginning of last tree block. */
		for (i=0; i<lastTreeBlockBegin+1; i++)
			fgets (s, longestLineLength, fp);	
	
		/* Set up cheap status bar. */
		if (runNo ==0)
			{
			MrBayesPrint ("\n%s   Rereading trees to process branch length values. Reading status:\n\n", spacer);
			MrBayesPrint ("%s   0      10      20      30      40      50      60      70      80      90     100\n", spacer);
			MrBayesPrint ("%s   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v\n", spacer);
			MrBayesPrint ("%s   *", spacer);
			numAsterices = 0;
			}

		/* Parse file, tree-by-tree. We are only parsing lines between the "begin trees" and "end" statements. */
		expecting = Expecting(COMMAND);
		numSumtTrees = numSumTreesSampled = 0;
		inTreesBlock = YES;
		ResetTranslateTable();
		for (i=0; i<lastTreeBlockEnd - lastTreeBlockBegin - 1; i++)
			{
			fgets (s, longestLineLength, fp);
			/*MrBayesPrint ("%s", s);*/
			if (ParseCommand (s) == ERROR)
				{
				free (brlens);
				SafeFclose (&fpBrlens);
				return ERROR;
				}
			}
		inTreesBlock = NO;
		
		/* Finish cheap status bar. */
		if (runNo == sumtParams.numRuns - 1)
			{
			if (numAsterices < 80)
				for (i=0; i<80 - numAsterices; i++)
					MrBayesPrint ("*");
			MrBayesPrint ("\n\n");
			}

		SafeFclose (&fp);
		}	/* next file */
		
	/* reset status variable */
	printingBrlens = NO;

	/* close brlens file */
	SafeFclose (&fpBrlens);
	MrBayesPrint ("%s   Branch length values printed to file.\n", spacer);

	return NO_ERROR;
}





void CalculateTreeToTreeDistance (int *lst[2], MrBFlt *lngs[2], int nnds, int hasBrlens, MrBFlt *d1, MrBFlt *d2, MrBFlt *d3)

{

	int			i, j, foundPart, *identifiedParts, *list1, *list2;
	MrBFlt		*lengths1, *lengths2, tl1=0.0, tl2=0.0, diff1, diff2;
	
	(*d1) = (*d2) = (*d3) = 0.0;

	list1 = lst[0];
	list2 = lst[1];
	lengths1 = lngs[0];
	lengths2 = lngs[1];

	identifiedParts = (int *)SafeMalloc((size_t) (nnds * sizeof(int)));
	if (!identifiedParts)
		{
		MrBayesPrint ("%s   Could not allocate identifiedParts\n", spacer);
		}
	for (i=0; i<nnds; i++)
		identifiedParts[i] = 0;
		
	if (hasBrlens == YES)
		{
		tl1 = tl2 = 0.0;
		for (i=0; i<nnds; i++)
			{
			tl1 += lengths1[i];
			tl2 += lengths2[i];
			}
		}
	
	for (i=0; i<nnds; i++)
		{
		foundPart = NO;
		for (j=0; j<nnds; j++)
			{
			if (list2[j] == list1[i])
				{
				foundPart = YES;
				break;
				}
			}
		if (foundPart == YES)
			{
			if (hasBrlens == YES)
				{
				diff1 = lengths1[i] - lengths2[j];
				if (diff1 < 0.0)
					diff1 = -diff1;
				diff2 = lengths1[i]/tl1 - lengths2[j]/tl2;
				if (diff2 < 0.0)
					diff2 = -diff2;
				(*d2) += diff1;
				(*d3) += diff2;
				}
			identifiedParts[j] = 1;
			}
		else
			{
			if (hasBrlens == YES)
				{
				diff1 = lengths1[i];
				if (diff1 < 0.0)
					diff1 = -diff1;
				diff2 = lengths1[i]/tl1;
				if (diff2 < 0.0)
					diff2 = -diff2;
				(*d2) += diff1;
				(*d3) += diff2;
				}
			(*d1) += 2.0;
			}
		}
		
	for (i=0; i<nnds; i++)
		{
		if (identifiedParts[i] == 0)
			{
			if (hasBrlens == YES)
				{
				diff1 = lengths2[i];
				if (diff1 < 0.0)
					diff1 = -diff1;
				diff2 = lengths2[i]/tl2;
				if (diff2 < 0.0)
					diff2 = -diff2;
				(*d2) += diff1;
				(*d3) += diff2;
				}
			}
		}

#	if 0		
	printf ("DISTANCES: %lf %lf %lf (%lf %lf)\n", *d1, *d2, *d3, tl1, tl2);
	for (i=0; i<nnds; i++)
		{
		printf ("%4d -- %4d (%lf) %4d (%lf)\n", i, list1[i], lengths1[i], list2[i], lengths2[i]);
		}
#	endif
						
	free (identifiedParts);

}





int CheckSumtSpecies (void)

{

	int 			i, nNodes, whichTaxon, numNotFound;
	SumtNode		**downPass, *q;

	/* allocate memory for downpass */
	downPass = (SumtNode **)SafeMalloc((size_t) (2 * numTaxa * sizeof(SumtNode *)));
	if (!downPass)
		{
		MrBayesPrint ("%s   Could not allocate downPass\n", spacer);
		goto errorExit;
		}
		
	/* get the downpass sequence */
	i = 0;
	GetSumtDownPass (sumtRoot, downPass, &i);
	nNodes = i;
	
	/* find which taxa are in the tree */
	for (i=0; i<nNodes; i++)
		{
		q = downPass[i];
		if ((q->left == NULL && q->right == NULL) || (q->anc == NULL && isSumtTreeRooted == NO))
			{
			if (CheckString (q->label, taxaNames, &whichTaxon) == ERROR)
				{
				MrBayesPrint ("%s   Could not find taxon %s in original list of taxa\n", spacer, q->label);
				goto errorExit;
				}
			whichTaxon--;
			sumTaxaFound[whichTaxon] = YES;
			}
		}

	/* now, print out which taxa are not in list */
	numNotFound = 0;
	for (i=0; i<numTaxa; i++)
		if (sumTaxaFound[i] == NO)
			absentTaxa[i] = YES;
		
	/* clean up on way out of function */
	free (downPass);

	return (NO_ERROR);
	
	errorExit:
		if (downPass)
			free (downPass);
		return(ERROR);

}





int ConTree (void)

{

	int			i, j, targetNode, nBits, nextConNode, isCompat, localOutgroupNum, numTerminalsEncountered;
	safeLong		x, *mask = NULL, *partition = NULL, *ingroupPartition = NULL, *outgroupPartition = NULL;
	MrBFlt		freq;
	char		tempName[100];
	PolyNode	*cp, *q, *r, *ql, *rl, *pl, **downPass = NULL;
	
	/* check that we have at least three species */
	j = 0;
	for (i=0; i<numTaxa; i++)
		if (sumTaxaFound[i] == YES)
			j++;
	if (j < 3)
		{
		MrBayesPrint ("%s   Too few taxa included to show consensus trees\n", spacer);
		goto errorExit;
		}
	
	/* Set the outgroup. Remember that the outgroup number goes from 0 to numTaxa-1.
	   The outgroup may have been deleted, so we should probably set localOutgroupNum
	   to reflect this. */
	j = 0;
	localOutgroupNum = 0;
	for (i=0; i<numTaxa; i++)
		{
		if (sumTaxaFound[i] == YES)
			{
			if (i == outGroupNum)
				{
				localOutgroupNum = j;
				break;
				}
			j++;
			}
		}
	
	/* now, make a consensus tree */
	
	/* note that numIncludedTaxa is initialized in ReorderParts */

	/* first allocate some stuff for the consensus tree */
	if (memAllocs[ALLOC_CONNODES] == YES)
		{
		MrBayesPrint ("%s   conNodes is already allocated\n", spacer);
		goto errorExit;
		}
	conNodes = (PolyNode *)SafeMalloc((size_t) (2 * numTaxa * sizeof(PolyNode)));
	if (!conNodes)
		{
		MrBayesPrint ("%s   Could not allocate conNodes\n", spacer);
		goto errorExit;
		}
	memAllocs[ALLOC_CONNODES] = YES;
	for (i=0; i<2*numTaxa; i++)
		{
		conNodes[i].left = conNodes[i].sib = conNodes[i].anc = NULL;
		conNodes[i].x = conNodes[i].y = conNodes[i].index = conNodes[i].mark = 0;
		conNodes[i].length = conNodes[i].support = conNodes[i].f = 0.0;
		}
	if (memAllocs[ALLOC_OUTPART] == YES)
		{
		MrBayesPrint ("%s   outgroupPartition is already allocated\n", spacer);
		goto errorExit;
		}
	outgroupPartition = (safeLong *) calloc (3 * taxonLongsNeeded, sizeof(safeLong));
	if (!outgroupPartition)
		{
		MrBayesPrint ("%s   Could not allocate outgroupPartition\n", spacer);
		goto errorExit;
		}
	ingroupPartition = outgroupPartition + taxonLongsNeeded;
	mask = ingroupPartition + taxonLongsNeeded;
	memAllocs[ALLOC_OUTPART] = YES;

	/* initialize terminal consensus nodes */
	j = 0;
	for (i=0; i<numTaxa; i++)
		{
		if (sumTaxaFound[i] == YES)
			{
			if (GetNameFromString (taxaNames, tempName, i+1) == ERROR)
				{
				MrBayesPrint ("%s   Error getting taxon names \n", spacer);
				return (ERROR);
				}
			conNodes[j].left = NULL;
			conNodes[j].sib = NULL;
			conNodes[j].index = j;
			strcpy(conNodes[j].label,tempName);
			j++;
			}
		}
	for (i=numIncludedTaxa; i<2*numIncludedTaxa; i++)
		{
		conNodes[j].left = NULL;
		conNodes[j].sib = NULL;
		conNodes[j].index = j;
		strcpy (conNodes[j].label, "");
		j++;
		}
		
	/* Set mask - needed to trim last bits in the partition bit field.
	   This could be done when a new matrix is read in
	   and adjusted when taxa are deleted or restored. */
	for (i=0; i<numIncludedTaxa; i++)
		SetBit (i, mask);

#	if defined (DEBUG_CONTREE)
	for (i=0; i<taxonLongsNeeded; i++)
		ShowBits (&mask[i], nBitsInALong);
	MrBayesPrint (" <- lastMask \n");
#	endif

	/* Set ingroup and outgroup partitions.
	   This could be done when a new matrix is read in
	   and adjusted when an outgroup command is issued.
	   This mechanism allows multiple taxa in outgroup 
	   but only one outgroup taxon is used here. */
	x = 1;
	x <<= (localOutgroupNum) % nBitsInALong;
	i = (localOutgroupNum) / nBitsInALong;
	outgroupPartition[i] = x;
	for (i = 0; i < taxonLongsNeeded; i++)
		ingroupPartition[i] = outgroupPartition[i];
	FlipBits (ingroupPartition, taxonLongsNeeded, mask);
#	if defined (DEBUG_CONTREE)
	ShowBits (&outgroupPartition[0], numIncludedTaxa);
	MrBayesPrint (" <- outgroupPartition \n");
	ShowBits (&ingroupPartition[0], numIncludedTaxa);
	MrBayesPrint (" <- ingroupPartition \n");
	MrBayesPrint ("%s   Hold: ", spacer);
	fgets (tempName, 100, stdin);
#	endif
	
	/* create bush 
	   ->x counts number of subtended terminals 
	   make sure conRoot->left is in outgroup */
	conRoot = &conNodes[numIncludedTaxa];
	conRoot->anc = conRoot->sib = NULL;
	conRoot->x = numIncludedTaxa;
	j = FirstTaxonInPartition (outgroupPartition, taxonLongsNeeded);
	conRoot->left = cp = &conNodes[j];
	cp->anc = conRoot;
	cp->x = 1;
	for (i=0; i<numIncludedTaxa; i++)
		{
		if (i != j)
			{
			cp->sib = &conNodes[i];
			cp = cp->sib;
			cp->anc = conRoot;
			cp->x = 1;
			}
		}
	cp->sib = NULL;

	/* Resolve bush according to partitions.
	   Partitions may include incompatible ones.
	   Partitions should be sorted from most frequent to least frequent 
	   for quit test to work when a 50% majority rule tree is requested. */
	nextConNode = numIncludedTaxa + 1;
	if (isSumtTreeRooted == YES)
		targetNode = 2 * numIncludedTaxa - 2;
	else
		targetNode = 2 * numIncludedTaxa - 3;

	numTerminalsEncountered = 0;
	for (i=0; i<numTreePartsFound; i++)
		{
		/* calculate frequency and test if time to quit */
		if (nextConNode > targetNode && numTerminalsEncountered == numIncludedTaxa)
			break;
		freq = (MrBFlt)numFoundOfThisPart[i]/ (MrBFlt)(sumtParams.numRuns * numSumTreesSampled);
		if (freq < 0.50 && !strcmp(sumtParams.sumtConType, "Halfcompat"))
			break;
		
		/* get partition */
		partition = &treePartsFound[i*taxonLongsNeeded];

		/* flip bits if necessary */
		/* This code is needed if single outgroup is indexed incorrectly or if the partition
		   defines a clade in a multispecies outgroup but the indexing is reversed. Note that
		   bits should not be flipped for rooted trees. */
		if (isSumtTreeRooted == NO)
			{
			if (!IsPartNested(partition, ingroupPartition, taxonLongsNeeded) && !IsPartCompatible(partition, ingroupPartition, taxonLongsNeeded))
				FlipBits(partition, taxonLongsNeeded, mask);
			}
		
		/* count bits in this partition */
		for (j=nBits=0; j<taxonLongsNeeded; j++)
			{
			x = partition[j];
			for (x = partition[j]; x != 0; x &= (x - 1))
				nBits++;
			}
			
		/* flip this partition if it leaves single outgroup outside and tree is unrooted */
		if (nBits == numIncludedTaxa - 1  && isSumtTreeRooted == NO)
			{
			nBits = 1;
			FlipBits(partition, taxonLongsNeeded, mask);
			}

		/* find out if this is an informative partition */
		if (nBits == numIncludedTaxa)
			{
			/* this is the root (for seeting age of root node when nodeDepthConTree == YES */
			q = conRoot;
			q->age = aBrlens[i];
			q->support = freq * 100.0;
			}
		else if (nBits > 1)
			{
			/* this is an informative partition */
			/* find anc of partition */
			j = FirstTaxonInPartition (partition, taxonLongsNeeded);
			for (cp = &conNodes[j]; cp!=NULL; cp = cp->anc)
				if (cp->x > nBits)
					break;
					
			/* do not include if incompatible with ancestor or any of descendants
			   do not check terminals or root because it is
			   redundant and partitions have not been set for those */
			isCompat = YES;
			if (cp->anc != NULL && !IsPartCompatible(partition, cp->partition, taxonLongsNeeded))
				isCompat = NO;
			for (q=cp->left; q!=NULL; q=q->sib)
				{
				if (q->x > 1 && !IsPartCompatible(q->partition, partition, taxonLongsNeeded))
					isCompat = NO;
				if (isCompat == NO)
					break;
				}
			if (isCompat == NO)
				continue;

			/* set new node */
			q = &conNodes[nextConNode++];
			if (sumtBrlensDef == YES)
				{
				if (nodeDepthConTree == YES)
					q->age = aBrlens[i];
				else
					q->length = aBrlens[i];
				}
			else
				q->length = 0.0;
			q->support = freq * 100;
			q->x = nBits;
			q->partition = partition;

			/* go through descendants of anc */
			ql = pl = NULL;
			for (r=cp->left; r!=NULL; r=r ->sib)
				{
				/* test if r is in the new partition or not */
				if ((r->x > 1 && IsPartNested(r->partition, partition, taxonLongsNeeded)) || (r->x == 1 && (partition[r->index / nBitsInALong] & (1 << (r->index % nBitsInALong))) != 0))
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
						cp->left = r;
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
			q->anc = cp;
			}
		else
			/* singleton partition */
			{
			j = FirstTaxonInPartition(partition, taxonLongsNeeded);
			q = &conNodes[j];
			if (sumtBrlensDef == YES)
				{
				if (nodeDepthConTree == YES)
					q->age = aBrlens[i];
				else
					q->length = aBrlens[i];
				}
			else
				q->length = 0.0;
			q->support = freq * 100;
			numTerminalsEncountered++;
			}
		}

	if (sumtParams.orderTaxa == YES)
		{
		/* rearrange tree so that terminals are in order */
		/* first allocate space for downPass */
		downPass = (PolyNode **) calloc (nextConNode, sizeof (PolyNode *));	
		if (!downPass)
			return ERROR;
		i = 0;
		GetConDownPass (downPass, conRoot, &i);

		/* label by minimum index */
		for (i=0; i<nextConNode; i++)
			{
			cp = downPass[i];
			if (cp->left == NULL)
				{
				if (cp->index == localOutgroupNum)
					cp->x = -1;
				else
					cp->x = cp->index;
				}
			else
				{
				j = nextConNode;
				for (q=cp->left; q!=NULL; q=q->sib)
					{
					if (q->x < j)
						j = q->x;
					}
				cp->x = j;
				}
			}
		/* and rearrange */
		for (i=0; i<nextConNode; i++)
			{
			cp = downPass[i];
			if (cp->left == NULL || cp->anc == NULL)
				continue;
			for (ql=NULL, q=cp->left; q->sib!=NULL; ql=q, q=q->sib)
				{
				for (rl=q, r=q->sib; r!=NULL; rl=r, r=r->sib)
					{
					if (r->x < q->x)
						{
						if (ql == NULL)
							cp->left = r;
						if (r == q->sib) /* swap adjacent q and r */
							{
							if (ql != NULL)
								ql->sib = r;
							pl = r->sib;
							r->sib = q;
							q->sib = pl;
							}
						else	/* swap separated q and r */
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
		}

	/* take care of phylogram based on node depths */
	if (nodeDepthConTree == YES)
		{
		for (i=0; i<nextConNode; i++)
			{
			cp = downPass[i];
			if (cp->anc == NULL)
				cp->length = 0.0;
			else
				cp->length = cp->anc->age - cp->age;
			}
		}
		
	/* draw tree to stdout and fp */
	MrBayesPrint ("\n%s   Clade credibility values:\n\n", spacer);
	ShowConTree (stdout, nextConNode, conRoot, 80, YES);
	if (logToFile == YES)
		ShowConTree (logFileFp, nextConNode, conRoot, 80, YES);
	if (sumtBrlensDef == YES)
		{
		MrBayesPrint ("\n");
		if (nodeDepthConTree == YES)
			MrBayesPrint ("%s   Phylogram (based on average node depths):\n\n", spacer);
		else
			MrBayesPrint ("%s   Phylogram (based on average branch lengths):\n\n", spacer);
		if (nodeDepthConTree == YES && AreDoublesEqual(clockRate[0], 1.0, 0.000001) == NO)
			ShowConPhylogram (stdout, nextConNode, conRoot, 80, YES);
		else
			ShowConPhylogram (stdout, nextConNode, conRoot, 80, NO);
		if (logToFile == YES)
			{
			if (nodeDepthConTree == YES && AreDoublesEqual (clockRate[0], 1.0, 0.000001) == NO)
				ShowConPhylogram (stdout, nextConNode, conRoot, 80, YES);
			else
				ShowConPhylogram (stdout, nextConNode, conRoot, 80, NO);
			}
		}

    /* print taxa block */
    MrBayesPrintf (fpCon, "begin taxa;\n");
    MrBayesPrintf (fpCon, "\tdimensions ntax=%d;\n", numSumtTaxa);
    MrBayesPrintf (fpCon, "\ttaxlabels\n", numSumtTaxa);
    for (i=0; i<nextConNode; i++)
        {
        if (conNodes[i].left == NULL)
            MrBayesPrintf (fpCon, "\t\t%s\n", conNodes[i].label);
        }
    MrBayesPrintf (fpCon, "\t\t;end;\n");
    
	MrBayesPrintf (fpCon, "begin trees;\n");
    MrBayesPrintf (fpCon, "\ttranslate\n");
    for (i=0; i<numSumtTaxa; i++)
        {
        for (j=0; j<nextConNode; j++)
            if (conNodes[j].index == i) break;
        if (i == numSumtTaxa-1)
            MrBayesPrintf (fpCon, "\t\t%d\t%s\n", conNodes[i].index+1, conNodes[i].label);
        else
            MrBayesPrintf (fpCon, "\t\t%d\t%s,\n", conNodes[i].index+1, conNodes[i].label);
        }
    MrBayesPrintf (fpCon, "\t\t;\n");

    MrBayesPrintf (fpCon, "   [Note: This tree contains information on the topology, \n");
	MrBayesPrintf (fpCon, "          branch lengths (if present), and the probability\n");
	MrBayesPrintf (fpCon, "          of the partition indicated by the branch.]\n");
	if (!strcmp(sumtParams.sumtConType, "Halfcompat"))
		MrBayesPrintf (fpCon, "   tree con_50_majrule = ");
	else
		MrBayesPrintf (fpCon, "   tree con_all_compat = ");
	WriteConTree (conRoot, fpCon, YES);
	MrBayesPrintf (fpCon, ";\n");
	if (sumtBrlensDef == YES)
		{
		MrBayesPrintf (fpCon, "\n");
		MrBayesPrintf (fpCon, "   [Note: This tree contains information only on the topology\n");
		MrBayesPrintf (fpCon, "          and branch lengths (mean of the posterior probability density).]\n");
		if (!strcmp(sumtParams.sumtConType, "Halfcompat"))
			MrBayesPrintf (fpCon, "   tree con_50_majrule = ");
		else
			MrBayesPrintf (fpCon, "   tree con_all_compat = ");
		WriteConTree (conRoot, fpCon, NO);
		MrBayesPrintf (fpCon, ";\n");
		}
	MrBayesPrintf (fpCon, "end;\n");
	
	/* free memory and file pointers */
	if (memAllocs[ALLOC_CONNODES] == YES)
		{
		free (conNodes);
		memAllocs[ALLOC_CONNODES] = NO;
		}
	if (memAllocs[ALLOC_OUTPART] == YES)
		{
		free (outgroupPartition);
		memAllocs[ALLOC_OUTPART] = NO;
		}
		
	if (downPass)
		free (downPass);

	return (NO_ERROR);
	
	errorExit:
		if (memAllocs[ALLOC_CONNODES] == YES)
			{
			free (conNodes);
			memAllocs[ALLOC_CONNODES] = NO;
			}
		if (memAllocs[ALLOC_OUTPART] == YES)
			{
			free (outgroupPartition);
			memAllocs[ALLOC_OUTPART] = NO;
			}
		if (downPass)
			free (downPass);

		return (ERROR);
}





int DerootSumtTree (SumtNode *p, int n, int outGrp)

{

	int 			i, nNodes, isMarked, *localTaxaFound=NULL, sumtOut;
	MrBFlt			tempBrLen;
	SumtNode		**downPass=NULL, *lft, *m1, *m2, *um1, *um2, *out;

#	if 0
	ShowTree (sumtRoot, YES, n);
#	endif

	/* check that we are not already unrooted */
	if (isSumtTreeRooted == NO)
		{
		MrBayesPrint ("%s   Tree is already unrooted\n", spacer);
		goto errorExit;
		}
	
	/* Find the outgroup number (0, 1, ..., n-1) and the number of nodes on the tree. The
	   number of nodes may change later, if the user deleted taxa. */
	sumtOut = outGrp;
	nNodes = 2 * n;

	/* allocate space for derooting tree */
	downPass = (SumtNode **)SafeMalloc((size_t) (2 * numTaxa * sizeof(SumtNode *)));
	if (!downPass)
		{
		MrBayesPrint ("%s   Could not allocate downPass\n", spacer);
		goto errorExit;
		}
	localTaxaFound = (int *)SafeMalloc((size_t) (numTaxa * sizeof(int)));
	if (!localTaxaFound)
		{
		MrBayesPrint ("%s   Could not allocate localTaxaFound\n", spacer);
		goto errorExit;
		}
	for (i=0; i<numTaxa; i++)
		localTaxaFound[i] = NO;
		
	/* get the downpass sequence for the tree */
	i = 0;
	GetSumtDownPass (sumtRoot, downPass, &i);
	nNodes = i;

	/* bring the outgroup around to the first right position */

	/* first, mark the outgroup tip and the path from the outgroup to the root */
	isMarked = NO;
	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		p->marked = NO;
		if (p->left == NULL && p->right == NULL && p->anc != NULL)
			{
			localTaxaFound[p->index] = YES;
			if (p->index == sumtOut)
				{
				p->marked = YES;
				isMarked = YES;
				}
			}
		}	
	
	/* If we don't find the outgroup, it was deleted by the user. Instead, we
	   will use the first undeleted taxon in the matrix that is found in the
	   tree. */
	if (isMarked == NO)
		{
		/* here we find the first undeleted taxon and designate it as the outgroup */
		for (i=0; i<numTaxa; i++)
			{
			if (localTaxaFound[i] == YES)
				{
				sumtOut = i;
				break;
				}
			}
		/* and then mark it */
		isMarked = NO;
		for (i=0; i<nNodes; i++)
			{
			p = downPass[i];
			p->marked = NO;
			if (p->left == NULL && p->right == NULL && p->anc != NULL)
				{
				if (p->index == sumtOut)
					{
					p->marked = YES;
					isMarked = YES;
					}
				}
			}	
		/* if we still have not marked an outgroup, we have trouble */
		if (isMarked == NO)
			{
			MrBayesPrint ("%s   Could not find outgroup taxon\n", spacer);
			goto errorExit;
			}
		}
		
	/* now we mark the path from the outgroup tip to the root */
	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		if (p->left != NULL && p->right != NULL)
			if (p->left->marked == YES || p->right->marked == YES)
				p->marked = YES;
		}	
		
	/* now we rotate the tree until the outgroup is to the left or to the right of the root */
	lft = sumtRoot->left;
	while (lft->left->index != sumtOut && lft->right->index != sumtOut)
		{
		if (lft->left->marked == YES && lft->right->marked == NO)
			{
			m1 = lft->left;
			um1 = lft->right;
			if (m1->left != NULL && m1->right != NULL)
				{
				if (m1->left->marked == YES)
					{
					m2 = m1->left;
					um2 = m1->right;
					}
				else
					{
					m2 = m1->right;
					um2 = m1->left;
					}
				lft->left = m2;
				lft->right = m1;
				m2->anc = m1->anc = lft;
				m1->left = um2;
				m1->right = um1;
				um1->anc = um2->anc = m1;
				m1->marked = NO;
				um1->length += m1->length;
				m2->length *= 0.5;
				m1->length = m2->length;
				}
			else
				{
				MrBayesPrint ("%s   Rooting routine is lost (1)\n", spacer);
				goto errorExit;
				}
			}
		else if (lft->left->marked == NO && lft->right->marked == YES)
			{
			m1 = lft->right;
			um1 = lft->left;
			if (m1->left != NULL && m1->right != NULL)
				{
				if (m1->left->marked == YES)
					{
					m2 = m1->left;
					um2 = m1->right;
					}
				else
					{
					m2 = m1->right;
					um2 = m1->left;
					}
				lft->left = m1;
				lft->right = m2;
				m2->anc = m1->anc = lft;
				m1->left = um1;
				m1->right = um2;
				um1->anc = um2->anc = m1;
				m1->marked = NO;
				um1->length += m1->length;
				m2->length *= 0.5;
				m1->length = m2->length;
				}
			else
				{
				MrBayesPrint ("%s   Rooting routine is lost (2)\n", spacer);
				goto errorExit;
				}
			}
		else
			{
			MrBayesPrint ("%s   Rooting routine is lost (3)\n", spacer);
			goto errorExit;
			}
		}

	/* make certain outgroup is to the right of the root */
	if (sumtRoot->left->left->index == sumtOut)
		{
		m1 = sumtRoot->left->left;
		m2 = sumtRoot->left->right;
		lft = sumtRoot->left;
		lft->left = m2;
		lft->right = m1;
		}
		
	/* now, take outgroup and make it point down */
	m1 = sumtRoot;
	m2 = sumtRoot->left;
	lft = sumtRoot->left->left;
	out = sumtRoot->left->right;
	tempBrLen = out->length;
	if (tempBrLen < 0.0)
		tempBrLen = 0.0;
	lft->anc = out;
	out->left = lft;
	out->right = out->anc = NULL;
	lft->length += tempBrLen;
	m1->left = m1->right = m1->anc = NULL;
	m2->left = m2->right = m2->anc = NULL;
	sumtRoot = out;
	
	/* the tree is now unrooted */
	isSumtTreeRooted = NO;

	/* reindex internal nodes of tree */
	i = numTaxa;
	FinishSumtTree (sumtRoot, &i, NO);

#	if 0
	ShowNodes (sumtRoot, 3, isSumtTreeRooted);
	ShowTree (sumtRoot, isSumtTreeRooted, n);
#	endif

	/* free memory */
	free (downPass);
	free (localTaxaFound);

	return (NO_ERROR);
	
	errorExit:
		if (downPass)
			free (downPass);
		if (localTaxaFound)
			free (localTaxaFound);
		return(ERROR);

}





int DoCompareTree (void)

{


	int			i, j, k, n, lineTerm, longestLineLength, tokenType, foundBegin, inTreeBlock, lineNum, 
				numTreesInBlock, numTreeBlocks, lastTreeBlockBegin, lastTreeBlockEnd,
				lineWidth, blockErrors, inSumtComment, numExcludedTaxa,
				longestLineLength1, longestLineLength2, xaxis, yaxis, starHolder[80], *treePartList[2], minNumTrees,
				screenWidth, screenHeigth, numY[60], nSamples, nCompPartitions, numIncludedTaxa1, numIncludedTaxa2,
				oldSumtBrlensDef, brlensDef;
	safeLong	temporarySeed, len;
	safeLong		*x;
	MrBFlt		xProb, yProb, xInc, yInc, xUpper, xLower, yUpper, yLower, *treeLengthList[2], *dT1=NULL, *dT2=NULL, *dT3=NULL, d1, d2, d3, 
				meanY[60], xVal, yVal, minX, minY, maxX, maxY, sums[3];
	char		*s=NULL, tempStr[100], tempName[100], prCh;
	FILE		*fp[2];
	time_t		curTime;
	
#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
#	endif

    /* Make sure we read trees using sumt code instead of with the user tree code */
    inComparetreeCommand = YES;

    /* Are we comparing two files, or using sumt? */
	comparingFiles = YES;
	fileNum = 0;
	oldSumtBrlensDef = sumtBrlensDef;
	
	/* we do not want to print brlens to file */
	printingBrlens = NO;
	
	/* set file pointers to NULL */
	fp[0] = fp[1] = NULL;

	/* Check that a data set has been read in. We check taxon names against
	   those read in. */
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before comparetree can be used\n", spacer);
		goto errorExit;
		}

	/* open binary file */
	if ((fp[0] = OpenBinaryFileR(comptreeParams.comptFileName1)) == NULL)
		goto errorExit;
	if ((fp[1] = OpenBinaryFileR(comptreeParams.comptFileName2)) == NULL)
		goto errorExit;
		
	/* find out what type of line termination is used for file 1 */
	lineTerm = LineTermType (fp[0]);
	if (lineTerm == LINETERM_MAC)
		MrBayesPrint ("%s   Macintosh line termination for file 1\n", spacer);
	else if (lineTerm == LINETERM_DOS)
		MrBayesPrint ("%s   DOS line termination for file 1\n", spacer);
	else if (lineTerm == LINETERM_UNIX)
		MrBayesPrint ("%s   UNIX line termination for file 1\n", spacer);
	else
		{
		MrBayesPrint ("%s   Unknown line termination for file 1\n", spacer);
		goto errorExit;
		}

	/* find out what type of line termination is used for file 2 */
	lineTerm = LineTermType (fp[1]);
	if (lineTerm == LINETERM_MAC)
		MrBayesPrint ("%s   Macintosh line termination for file 2\n", spacer);
	else if (lineTerm == LINETERM_DOS)
		MrBayesPrint ("%s   DOS line termination for file 2\n", spacer);
	else if (lineTerm == LINETERM_UNIX)
		MrBayesPrint ("%s   UNIX line termination for file 2\n", spacer);
	else
		{
		MrBayesPrint ("%s   Unknown line termination for file 2\n", spacer);
		goto errorExit;
		}

	/* find length of longest line in either file */
	longestLineLength1 = LongestLine (fp[0]);
	longestLineLength2 = LongestLine (fp[1]);
	if (longestLineLength1 > longestLineLength2)
		longestLineLength = longestLineLength1;
	else
		longestLineLength = longestLineLength2;
	MrBayesPrint ("%s   Longest line length = %d\n", spacer, longestLineLength);
	longestLineLength += 10;
	
	/* allocate a string long enough to hold a line */
	if (memAllocs[ALLOC_SUMTSTRING] == YES)
		{
		MrBayesPrint ("%s   Comparetree string is already allocated\n", spacer);
		goto errorExit;
		}
	s = (char *)SafeMalloc((size_t) (longestLineLength * sizeof(char)));
	if (!s)
		{
		MrBayesPrint ("%s   Problem allocating string for reading comparetree file\n", spacer);
		goto errorExit;
		}
	memAllocs[ALLOC_SUMTSTRING] = YES;
		
	/* close binary file */
	SafeFclose (&fp[0]);
	SafeFclose (&fp[1]);
	
	/* read in data file 1 ***************************************************************************/

	/* tell user we are ready to go */
	MrBayesPrint ("%s   Summarizing trees in file %s\n", spacer, comptreeParams.comptFileName1);
	
	/* open text file */
	if ((fp[0] = OpenTextFileR(comptreeParams.comptFileName1)) == NULL)
		goto errorExit;
	
	/* Check file for appropriate blocks. We want to find the last tree block
	   in the file and start from there. */
	foundBegin = inTreeBlock = blockErrors = inSumtComment = NO;
	lineNum = numTreesInBlock = lastTreeBlockBegin = lastTreeBlockEnd = numTreeBlocks = numTreesInLastBlock = 0;
	while (fgets (s, longestLineLength, fp[0]) != NULL)
		{
		sumtTokenP = &s[0];
		do
			{
			GetSumtToken (&tokenType);
			if (IsSame("[", sumtToken) == SAME)
				inSumtComment = YES;
			if (IsSame("]", sumtToken) == SAME)
				inSumtComment = NO;
				
			if (inSumtComment == NO)
				{
				if (foundBegin == YES)
					{
					if (IsSame("Trees", sumtToken) == SAME)
						{
						numTreesInBlock = 0;
						inTreeBlock = YES;
						foundBegin = NO;
						lastTreeBlockBegin = lineNum;
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
							numTreeBlocks++;
							inTreeBlock = NO;
							lastTreeBlockEnd = lineNum;
							}
						else
							{
							MrBayesPrint ("%s   Found inappropriate \"End\" statement in file\n", spacer);
							blockErrors = YES;
							}
						numTreesInLastBlock = numTreesInBlock;
						}
					else if (IsSame("Tree", sumtToken) == SAME)
						{
						if (inTreeBlock == YES)
							{
							brlensDef = NO;
							for (j=0; s[j]!='\0'; j++)
								{
								if (s[j] == ':')
									{
									brlensDef = YES;
									break;
									}
								else if (s[j] == ',')
									break;
								}
							sumtBrlensDef = brlensDef;
							numTreesInBlock++;
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
		MrBayesPrint ("%s   Unterminated tree block in file %s. You probably need to\n", spacer, comptreeParams.comptFileName1);
		MrBayesPrint ("%s   add a new line to the end of the file with \"End;\" on it.\n", spacer);
		goto errorExit;
		}
	if (inSumtComment == YES)
		{
		MrBayesPrint ("%s   Unterminated comment in file %s\n", spacer, comptreeParams.comptFileName1);
		goto errorExit;
		}
	if (blockErrors == YES)
		{
		MrBayesPrint ("%s   Found formatting errors in file %s\n", spacer, comptreeParams.comptFileName1);
		goto errorExit;
		}
	if (lastTreeBlockEnd < lastTreeBlockBegin)
		{
		MrBayesPrint ("%s   Problem reading tree file %s\n", spacer, comptreeParams.comptFileName1);
		goto errorExit;
		}
	if (numTreesInLastBlock <= 0)
		{
		MrBayesPrint ("%s   No trees were found in last tree block of file %s\n", spacer, comptreeParams.comptFileName1);
		goto errorExit;
		}
    if (sumtParams.relativeBurnin == NO && sumtParams.sumtBurnIn > numTreesInLastBlock)
		{
		MrBayesPrint ("%s   No trees are sampled as the burnin exceeds the number of trees in last block\n", spacer);
		MrBayesPrint ("%s   Try setting burnin to a number less than %d\n", spacer, numTreesInLastBlock);
		goto errorExit;
		}
		
	/* tell the user that everything is fine */
	if (numTreeBlocks == 1)
		MrBayesPrint ("%s   Found one tree block in file \"%s\" with %d trees in last block\n", spacer, comptreeParams.comptFileName1, numTreesInLastBlock);
	else
		{
		MrBayesPrint ("%s   Found %d tree blocks in file \"%s\" with %d trees in last block\n", spacer, numTreeBlocks, comptreeParams.comptFileName1, numTreesInLastBlock);
		MrBayesPrint ("%s   Only the %d trees in last tree block will be summarized\n", spacer, numTreesInLastBlock);
		}
		
	/* Now we read the file for real. First, rewind file pointer to beginning of file... */
	(void)fseek(fp[0], 0L, 0);	
	
	/* ...and fast forward to beginning of last tree block. */
	for (i=0; i<lastTreeBlockBegin+1; i++)
		fgets (s, longestLineLength, fp[0]);
		
	/* Allocate things we will need for trees... */
	if (AllocBits (numTaxa) == ERROR)
		goto errorExit;
		
	/* Set up cheap status bar. */
	MrBayesPrint ("\n%s   Tree reading status:\n\n", spacer);
	MrBayesPrint ("%s   0      10      20      30      40      50      60      70      80      90     100\n", spacer);
	MrBayesPrint ("%s   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v\n", spacer);
	MrBayesPrint ("%s   *", spacer);
	numAsterices = 0;
		
	/* ...and parse file, tree-by-tree. We are only parsing lines between the "begin trees" and "end" statements.
	   We don't actually get those lines, however, but rather the lines between those statements. */
	expecting = Expecting(COMMAND);
	numSumtTrees = numSumTreesSampled = 0;
	numCompTrees[0] = numCompTreesSampled[0] = 0;
	numFullCompTreesFound[0] = 0;
	inTreesBlock = YES;
    ResetTranslateTable();
	for (i=0; i<lastTreeBlockEnd - lastTreeBlockBegin - 1; i++)
		{
		fgets (s, longestLineLength, fp[0]);
		/*MrBayesPrint ("%s", s);*/
		if (ParseCommand (s) == ERROR)
			goto errorExit;
		}
	inTreesBlock = NO;
	
	/* Finish cheap status bar. */
	if (numAsterices < 80)
		for (i=0; i<80 - numAsterices; i++)
			MrBayesPrint ("*");
	MrBayesPrint ("\n\n");
	
	/* how many taxa were included */
	numIncludedTaxa1 = 0;
	for (i=0; i<numTaxa; i++)
		if (absentTaxa[i] == NO && prunedTaxa[i] == NO)
			numIncludedTaxa1++;
			
	/* print out information on absent taxa */
	numExcludedTaxa = 0;
	for (i=0; i<numTaxa; i++)
		if (absentTaxa[i] == YES)
			numExcludedTaxa++;
	if (numExcludedTaxa > 0)
		{
		if (numExcludedTaxa == 1)
			MrBayesPrint ("%s   The following species was absent from trees:\n", spacer);
		else
			MrBayesPrint ("%s   The following %d species were absent from trees:\n", spacer, numExcludedTaxa);
		MrBayesPrint ("%s      ", spacer);
		j = lineWidth = 0;
		for (i=0; i<numTaxa; i++)
			{
			if (absentTaxa[i] == YES)
				{
				j++;
				if (GetNameFromString (taxaNames, tempStr, i+1) == ERROR)
					{
					MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
					return (ERROR);
					}
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
		if (prunedTaxa[i] == YES && absentTaxa[i] == NO)
			numExcludedTaxa++;
	
	if (numExcludedTaxa > 0)
		{
		if (numExcludedTaxa == 1)
			MrBayesPrint ("%s   The following species was pruned from trees:\n", spacer);
		else
			MrBayesPrint ("%s   The following %d species were pruned from trees:\n", spacer, numExcludedTaxa);
		MrBayesPrint ("%s      ", spacer);
		j = lineWidth = 0;
		for (i=0; i<numTaxa; i++)
			{
			if (prunedTaxa[i] == YES && absentTaxa[i] == NO)
				{
				j++;
				if (GetNameFromString (taxaNames, tempStr, i+1) == ERROR)
					{
					MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
					return (ERROR);
					}
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
	
	/* tell user how many trees were successfully read */
	MrBayesPrint ("%s   Read %d trees from last tree block (sampling %d of them)\n", spacer, numSumtTrees, numSumTreesSampled);
	
	/* Check that at least one tree was read in. */
	if (numSumTreesSampled <= 0)
		{
		MrBayesPrint ("%s   No trees read in\n", spacer);
		goto errorExit;
		}
		
		
	/* read in data file 2 ***************************************************************************/

	fileNum = 1;

	/* tell user we are ready to go */
	MrBayesPrint ("%s   Summarizing trees in file %s\n", spacer, comptreeParams.comptFileName2);
		
	/* open text file */
	if ((fp[1] = OpenTextFileR(comptreeParams.comptFileName2)) == NULL)
		goto errorExit;
	
	/* Check file for appropriate blocks. We want to find the last tree block
	   in the file and start from there. */
	foundBegin = inTreeBlock = blockErrors = inSumtComment = NO;
	lineNum = numTreesInBlock = lastTreeBlockBegin = lastTreeBlockEnd = numTreeBlocks = numTreesInLastBlock = 0;
	while (fgets (s, longestLineLength, fp[1]) != NULL)
		{
		sumtTokenP = &s[0];
		do
			{
			GetSumtToken (&tokenType);
			if (IsSame("[", sumtToken) == SAME)
				inSumtComment = YES;
			if (IsSame("]", sumtToken) == SAME)
				inSumtComment = NO;
				
			if (inSumtComment == NO)
				{
				if (foundBegin == YES)
					{
					if (IsSame("Trees", sumtToken) == SAME)
						{
						numTreesInBlock = 0;
						inTreeBlock = YES;
						foundBegin = NO;
						lastTreeBlockBegin = lineNum;
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
							numTreeBlocks++;
							inTreeBlock = NO;
							lastTreeBlockEnd = lineNum;
							}
						else
							{
							MrBayesPrint ("%s   Found inappropriate \"End\" statement in file\n", spacer);
							blockErrors = YES;
							}
						numTreesInLastBlock = numTreesInBlock;
						}
					else if (IsSame("Tree", sumtToken) == SAME)
						{
						if (inTreeBlock == YES)
							numTreesInBlock++;
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
		MrBayesPrint ("%s   Unterminated tree block in file %s. You probably need to\n", spacer, comptreeParams.comptFileName2);
		MrBayesPrint ("%s   add a new line to the end of the file with \"End;\" on it.\n", spacer);
		goto errorExit;
		}
	if (inSumtComment == YES)
		{
		MrBayesPrint ("%s   Unterminated comment in file %s\n", spacer, comptreeParams.comptFileName2);
		goto errorExit;
		}
	if (blockErrors == YES)
		{
		MrBayesPrint ("%s   Found formatting errors in file %s\n", spacer, comptreeParams.comptFileName2);
		goto errorExit;
		}
	if (lastTreeBlockEnd < lastTreeBlockBegin)
		{
		MrBayesPrint ("%s   Problem reading tree file %s\n", spacer, comptreeParams.comptFileName2);
		goto errorExit;
		}
	if (numTreesInLastBlock <= 0)
		{
		MrBayesPrint ("%s   No trees were found in last tree block of file %s\n", spacer, comptreeParams.comptFileName2);
		goto errorExit;
		}
    if (sumtParams.relativeBurnin == NO && sumtParams.sumtBurnIn > numTreesInLastBlock)
		{
		MrBayesPrint ("%s   No trees are sampled as the burnin exceeds the number of trees in last block\n", spacer);
		MrBayesPrint ("%s   Try setting burnin to a number less than %d\n", spacer, numTreesInLastBlock);
		goto errorExit;
		}
		
	/* tell the user that everything is fine */
	if (numTreeBlocks == 1)
		MrBayesPrint ("%s   Found one tree block in file \"%s\" with %d trees in last block\n", spacer, comptreeParams.comptFileName2, numTreesInLastBlock);
	else
		{
		MrBayesPrint ("%s   Found %d tree blocks in file \"%s\" with %d trees in last block\n", spacer, numTreeBlocks, comptreeParams.comptFileName2, numTreesInLastBlock);
		MrBayesPrint ("%s   Only the %d trees in last tree block will be summarized\n", spacer, numTreesInLastBlock);
		}
		
	/* Now we read the file for real. First, rewind file pointer to beginning of file... */
	(void)fseek(fp[1], 0L, 0);	
	
	/* ...and fast forward to beginning of last tree block. */
	for (i=0; i<lastTreeBlockBegin+1; i++)
		fgets (s, longestLineLength, fp[1]);
		
	/* Set up cheap status bar. */
	MrBayesPrint ("\n%s   Tree reading status:\n\n", spacer);
	MrBayesPrint ("%s   0      10      20      30      40      50      60      70      80      90     100\n", spacer);
	MrBayesPrint ("%s   v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v\n", spacer);
	MrBayesPrint ("%s   *", spacer);
	numAsterices = 0;
		
	/* ...and parse file, tree-by-tree. We are only parsing lines between the "begin trees" and "end" statements.
	   We don't actually get those lines, however, but rather the lines between those statements. */
	expecting = Expecting(COMMAND);
	/*numSumtTrees = numSumTreesSampled = 0;*/
	numCompTrees[1] = numCompTreesSampled[1] = 0;
	numFullCompTreesFound[1] = 0;
	inTreesBlock = YES;
    ResetTranslateTable();
	for (i=0; i<lastTreeBlockEnd - lastTreeBlockBegin - 1; i++)
		{
		fgets (s, longestLineLength, fp[1]);
		/*MrBayesPrint ("%s", s);*/
		if (ParseCommand (s) == ERROR)
			goto errorExit;
		}
	inTreesBlock = NO;
	
	/* Finish cheap status bar. */
	if (numAsterices < 80)
		for (i=0; i<80 - numAsterices; i++)
			MrBayesPrint ("*");
	MrBayesPrint ("\n\n");
	
	/* how many taxa were included */
	numIncludedTaxa2 = 0;
	for (i=0; i<numTaxa; i++)
		if (absentTaxa[i] == NO && prunedTaxa[i] == NO)
			numIncludedTaxa2++;

	/* print out information on absent taxa */
	numExcludedTaxa = 0;
	for (i=0; i<numTaxa; i++)
		if (absentTaxa[i] == YES)
			numExcludedTaxa++;
	if (numExcludedTaxa > 0)
		{
		if (numExcludedTaxa == 1)
			MrBayesPrint ("%s   The following species was absent from trees:\n", spacer);
		else
			MrBayesPrint ("%s   The following %d species were absent from trees:\n", spacer, numExcludedTaxa);
		MrBayesPrint ("%s      ", spacer);
		j = lineWidth = 0;
		for (i=0; i<numTaxa; i++)
			{
			if (absentTaxa[i] == YES)
				{
				j++;
				if (GetNameFromString (taxaNames, tempStr, i+1) == ERROR)
					{
					MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
					return (ERROR);
					}
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
		if (prunedTaxa[i] == YES && absentTaxa[i] == NO)
			numExcludedTaxa++;
	
	if (numExcludedTaxa > 0)
		{
		if (numExcludedTaxa == 1)
			MrBayesPrint ("%s   The following species was pruned from trees:\n", spacer);
		else
			MrBayesPrint ("%s   The following %d species were pruned from trees:\n", spacer, numExcludedTaxa);
		MrBayesPrint ("%s      ", spacer);
		j = lineWidth = 0;
		for (i=0; i<numTaxa; i++)
			{
			if (prunedTaxa[i] == YES && absentTaxa[i] == NO)
				{
				j++;
				if (GetNameFromString (taxaNames, tempStr, i+1) == ERROR)
					{
					MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
					return (ERROR);
					}
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
	
	/* tell user how many trees were successfully read */
	MrBayesPrint ("%s   Read %d trees from last tree block (sampling %d of them)\n", spacer, numCompTrees[1], numCompTreesSampled[1]);
	
	/* summarize information ***************************************************************************/
	
	if (numIncludedTaxa1 != numIncludedTaxa2)
		{
		MrBayesPrint ("%s   The number of taxa to be compared in each file is not the same\n", spacer);
		goto errorExit;
		}
	
	/* how many taxon bipartitions are there for each tree */
	if (isSumtTreeRooted == YES)
		nCompPartitions = 2 * numIncludedTaxa1 - 2;
	else
		nCompPartitions = 2 * numIncludedTaxa1 - 3;

	/* Check that at least one tree was read in. */
	if (numSumTreesSampled <= 0)
		{
		MrBayesPrint ("%s   No trees read in\n", spacer);
		goto errorExit;
		}
		
	/* Sort partitions... */
	if (memAllocs[ALLOC_PARTORIGORDER] == YES)
		{
		MrBayesPrint ("%s   numFoundOfThisPart not free in AllocBits\n", spacer);
		goto errorExit;
		}
	partOrigOrder = (int *)SafeMalloc((size_t) (numTreePartsFound * sizeof(int)));
	if (!partOrigOrder)
		{
		MrBayesPrint ("%s   Problem allocating partOrigOrder (%d)\n", spacer, numTreePartsFound * sizeof(int));
		goto errorExit;
		}
	memAllocs[ALLOC_PARTORIGORDER] = YES;
	for (i=0; i<numTreePartsFound; i++)
		partOrigOrder[i] = i;
	if (SortParts (numFoundOfThisPart, numTreePartsFound) == ERROR)
		goto errorExit;
	
	/* ...and reorder bits if some taxa were not included. */
	numIncludedTaxa = numTaxa;
	if (ReorderParts () == ERROR)
		goto errorExit;
		
	/* open output files for summary information (two files) */
	if (OpenComptFiles () == ERROR)
		goto errorExit;
		
	/* print to screen */
    MrBayesPrint ("                                                                                   \n");
	MrBayesPrint ("%s   General explanation:                                                          \n", spacer);
    MrBayesPrint ("                                                                                   \n");
    MrBayesPrint ("%s   A taxon bibartition is specified by removing a branch, thereby divid-         \n", spacer);
    MrBayesPrint ("%s   ing the species into those to the left and those to the right of the          \n", spacer);
    MrBayesPrint ("%s   branch. Here, taxa to one side of the removed branch are denoted \".\"        \n", spacer);
    MrBayesPrint ("%s   and those to the other side are denoted \"*\". The output includes the        \n", spacer);
    MrBayesPrint ("%s   bipartition number (sorted from highest to lowest probability), bi-           \n", spacer);
    MrBayesPrint ("%s   partition (e.g., ...**..), number of times the bipartition was ob-            \n", spacer);
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
		if (GetNameFromString (taxaNames, tempName, k+1) == ERROR)
			{
			MrBayesPrint ("%s   Error getting taxon names \n", spacer);
			return (ERROR);
			}
		if (sumTaxaFound[k] == YES)
			{
			MrBayesPrint ("%s   %4d -- %s\n", spacer, j++, tempName);
			}
		}
    MrBayesPrint ("                                                                                   \n");
	MrBayesPrint ("%s   List of taxon bipartitions found in tree file:                                \n\n", spacer);

	x = &treePartsFound[0];
	for (i=0; i<numTreePartsFound; i++)
		{
		if ((MrBFlt)numFoundOfThisPart[i]/numSumTreesSampled >= sumtParams.freqDisplay)
			{
			MrBayesPrint ("%s   %4d -- ", spacer, i+1);
			ShowParts (stdout, &x[0], numIncludedTaxa);

			MrBayesPrint ("   %4d %4d %1.3lf %1.3lf\n", 
			numFoundOfThisPart1[i], numFoundOfThisPart2[i], 
			(MrBFlt)numFoundOfThisPart1[i]/numCompTreesSampled[0], 
			(MrBFlt)numFoundOfThisPart2[i]/numCompTreesSampled[1]);
			
			MrBayesPrintf (fpCompParts, "%d\t%d\t%d\t%1.3lf\t%1.3lf\n", 
			i+1, numFoundOfThisPart1[i], numFoundOfThisPart2[i], 
			(MrBFlt)numFoundOfThisPart1[i]/numCompTreesSampled[0], 
			(MrBFlt)numFoundOfThisPart2[i]/numCompTreesSampled[1]);
			}
		x += taxonLongsNeeded;
		}
		
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
			for (i=0; i<numTreePartsFound; i++)
				{
				xProb = (MrBFlt)numFoundOfThisPart1[i]/numCompTreesSampled[0];
				yProb = (MrBFlt)numFoundOfThisPart2[i]/numCompTreesSampled[1];
				if (xProb > xLower && xProb <= xUpper && yProb > yLower && yProb <= yUpper)
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

	MrBayesPrint ("   ^                                                                              ^\n");
	MrBayesPrint ("  0.00                                                                          1.00\n");
		
	/* get tree-to-tree distances */
	minNumTrees = numFullCompTreesFound[0];
	if (numFullCompTreesFound[1] < minNumTrees)
		minNumTrees = numFullCompTreesFound[1];
	if (memAllocs[ALLOC_TOPO_DISTANCES] == YES)
		{
		MrBayesPrint ("%s   Topological distances all ready allocated\n", spacer);
		goto errorExit;
		}
	dT1 = (MrBFlt *)SafeMalloc((size_t) (minNumTrees * sizeof(MrBFlt)));
	if (!dT1)
		{
		MrBayesPrint ("%s   Problem allocating topological distances\n", spacer);
		goto errorExit;
		}
	dT2 = (MrBFlt *)SafeMalloc((size_t) (minNumTrees * sizeof(MrBFlt)));
	if (!dT2)
		{
		MrBayesPrint ("%s   Problem allocating topological distances\n", spacer);
		goto errorExit;
		}
	dT3 = (MrBFlt *)SafeMalloc((size_t) (minNumTrees * sizeof(MrBFlt)));
	if (!dT3)
		{
		MrBayesPrint ("%s   Problem allocating topological distances\n", spacer);
		goto errorExit;
		}
	memAllocs[ALLOC_TOPO_DISTANCES] = YES;
		
	for (i=0; i<minNumTrees; i++)
		{
		treePartList[0] = &fullCompTreePartIds1[i * 2 * numTaxa];
		treeLengthList[0] = &fullCompTreePartLengths1[i * 2 * numTaxa];
		treePartList[1] = &fullCompTreePartIds2[i * 2 * numTaxa];
		treeLengthList[1] = &fullCompTreePartLengths2[i * 2 * numTaxa];
		CalculateTreeToTreeDistance (treePartList, treeLengthList, nCompPartitions, sumtBrlensDef, &d1, &d2, &d3);
		dT1[i] = d1;
		dT2[i] = d2;
		dT3[i] = d3;
		}
		
	for (i=0; i<minNumTrees; i++)
		{
		/*MrBayesPrint ("%s   %4d -- %lf %lf %lf\n", spacer, i+1, dT1[i], dT2[i], dT3[i]);*/	
		if (sumtBrlensDef == YES)
			MrBayesPrintf (fpCompDists, "%d\t%lf\t%lf\t%lf\n", i+1, dT1[i], dT2[i], dT3[i]);	
		else
			MrBayesPrintf (fpCompDists, "%d\t%lf\n", i+1, dT1[i]);	
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
		xVal = (MrBFlt) (i + comptreeParams.comptBurnIn);
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
		xVal = (MrBFlt) (i + comptreeParams.comptBurnIn);
		yVal = dT1[i];
		k = (int)(((xVal - minX) / (maxX - minX)) * screenWidth);
		if (k >= screenWidth)
			k = screenWidth - 1;
		meanY[k] += yVal;
		numY[k]++;
		}
	MrBayesPrint ("\n   +");
	for (i=0; i<screenWidth; i++)
		MrBayesPrint ("-");
	MrBayesPrint ("+ %1.2lf\n", maxY);
	for (j=screenHeigth-1; j>=0; j--)
		{
		MrBayesPrint ("   |");
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
	MrBayesPrint ("   +");
	for (i=0; i<screenWidth; i++)
		{
		if (i % (screenWidth/10) == 0 && i != 0)
			MrBayesPrint ("+");
		else
			MrBayesPrint ("-");
		}
	MrBayesPrint ("+ %1.2lf\n", minY);
	MrBayesPrint ("   ^");
	for (i=0; i<screenWidth; i++)
		MrBayesPrint (" ");
	MrBayesPrint ("^\n");
	MrBayesPrint ("   %1.0lf", minX);
	for (i=0; i<screenWidth; i++)
		MrBayesPrint (" ");
	MrBayesPrint ("%1.0lf\n\n", maxX);

	if (sumtBrlensDef == YES)
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
				xVal = (MrBFlt) (i + comptreeParams.comptBurnIn);
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
				xVal = (MrBFlt) (i + comptreeParams.comptBurnIn);
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
			MrBayesPrint ("\n   +");
			for (i=0; i<screenWidth; i++)
				MrBayesPrint ("-");
			MrBayesPrint ("+ %1.2lf\n", maxY);
			for (j=screenHeigth-1; j>=0; j--)
				{
				MrBayesPrint ("   |");
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
			MrBayesPrint ("   +");
			for (i=0; i<screenWidth; i++)
				{
				if (i % (screenWidth/10) == 0 && i != 0)
					MrBayesPrint ("+");
				else
					MrBayesPrint ("-");
				}
			MrBayesPrint ("+ %1.2lf\n", minY);
			MrBayesPrint ("   ^");
			for (i=0; i<screenWidth; i++)
				MrBayesPrint (" ");
			MrBayesPrint ("^\n");
			MrBayesPrint ("   %1.0lf", minX);
			for (i=0; i<screenWidth; i++)
				MrBayesPrint (" ");
			MrBayesPrint ("%1.0lf\n\n", maxX);
			}
		}

	/* calculate average tree-to-tree distances */
	curTime = time(NULL);
	temporarySeed  = (safeLong)curTime;
	if (temporarySeed < 0)
		temporarySeed = -temporarySeed;
	sums[0] = sums[1] = sums[2] = 0.0;
	nSamples = 1000;
	for (n=0; n<nSamples; n++)
		{
		i = (int) RandomNumber(&temporarySeed) * minNumTrees;
		j = (int) RandomNumber(&temporarySeed) * minNumTrees;
		treePartList[0] = &fullCompTreePartIds1[i * 2 * numTaxa];
		treeLengthList[0] = &fullCompTreePartLengths1[i * 2 * numTaxa];
		treePartList[1] = &fullCompTreePartIds2[j * 2 * numTaxa];
		treeLengthList[1] = &fullCompTreePartLengths2[j * 2 * numTaxa];
		CalculateTreeToTreeDistance (treePartList, treeLengthList, nCompPartitions, sumtBrlensDef, &d1, &d2, &d3);
		sums[0] += d1;
		sums[1] += d2;
		sums[2] += d3;
		}
	MrBayesPrint ("%s   Mean tree-to-tree distances, based on %d trees randomly sampled from both files:\n\n", spacer, nSamples);
	MrBayesPrint ("%s                                 Mean(Robinson-Foulds) = %1.3lf\n", spacer, sums[0]/nSamples);
	if (sumtBrlensDef == YES)
		{
		MrBayesPrint ("%s             Mean(Robinson-Foulds with branch lengths) = %1.3lf\n", spacer, sums[1]/nSamples);
		MrBayesPrint ("%s      Mean(Robinson-Foulds with scaled branch lengths) = %1.3lf\n", spacer, sums[2]/nSamples);
		}

	/* free memory and file pointers */
	if (memAllocs[ALLOC_SUMTSTRING] == YES)
		{
		free (s);
		memAllocs[ALLOC_SUMTSTRING] = NO;
		}
	if (memAllocs[ALLOC_TOPO_DISTANCES] == YES)
		{
		free (dT1);
		free (dT2);
		free (dT3);
		memAllocs[ALLOC_TOPO_DISTANCES] = NO;
		}
	if (memAllocs[ALLOC_PARTORIGORDER] == YES)
		{
		free (partOrigOrder);
		memAllocs[ALLOC_PARTORIGORDER] = NO;
		}
	FreeBits ();
	sumtBrlensDef = oldSumtBrlensDef;
	expecting = Expecting(COMMAND);

	/* close files */
	SafeFclose (&fp[0]);
	SafeFclose (&fp[1]);
	SafeFclose (&fpCompParts);
	SafeFclose (&fpCompDists);
	
#	if defined (MPI_ENABLED)
		}
#	endif

    inComparetreeCommand = NO;
	return (NO_ERROR);
	
	/* error exit */			
	errorExit:
		sumtBrlensDef = oldSumtBrlensDef;
		expecting = Expecting(COMMAND);
		if (memAllocs[ALLOC_SUMTSTRING] == YES)
			{
			free (s);
			memAllocs[ALLOC_SUMTSTRING] = NO;
			}
		if (memAllocs[ALLOC_TOPO_DISTANCES] == YES)
			{
			free (dT1);
			free (dT2);
			free (dT3);
			memAllocs[ALLOC_TOPO_DISTANCES] = NO;
			}
		if (memAllocs[ALLOC_PARTORIGORDER] == YES)
			{
			free (partOrigOrder);
			memAllocs[ALLOC_PARTORIGORDER] = NO;
			}
		FreeBits ();
		SafeFclose (&fp[0]);
		SafeFclose (&fp[1]);
		SafeFclose (&fpCompParts);
		SafeFclose (&fpCompDists);
		strcpy (spacer, "");
		strcpy (sumtToken, "Comparetree");
		i = 0;
		if (FindValidCommand (sumtToken, &i) == ERROR)
			MrBayesPrint ("%s   Could not find comparetree\n", spacer);
        inComparetreeCommand = NO;
		return (ERROR);	
	
}





int DoCompareTreeParm (char *parmName, char *tkn)

{

	int			tempI;

	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before comparetree can be used\n", spacer);
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
		/* set Burnin (comptreeParams.comptBurnIn) *******************************************************/
		else if (!strcmp(parmName, "Burnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
				comptreeParams.comptBurnIn = tempI;
				MrBayesPrint ("%s   Setting comparetree burnin to %ld\n", spacer, comptreeParams.comptBurnIn);
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





int DoSumt (void)

{


	int		i, j=0, k, n, lineTerm=0, longestLineLength=0, tokenType, foundBegin, inTreeBlock, lineNum, 
			numTreesInBlock, numTreeBlocks, lastTreeBlockBegin=0, lastTreeBlockEnd=0,
			len, lineWidth, longestName, blockErrors, inSumtComment, numExcludedTaxa, treeNo,
			runNo, firstFileNumTrees=0, firstFileNumSumTreesSampled=0, numTreePartsToPrint, maxWidthID,
			maxWidthNumberPartitions, tableWidth=0, unreliable, oneUnreliable, firstSumtBrlensDef=0,
			numTotalTreesSampled, precision, nSummarySamples=0, nPSRFSamples=0, includeInStat=NO;
	MrBFlt		var=0.0, f, var_s, sum_s, stddev_s=0.0, sumsq_s, varB=0.0, varW=0.0, sqrt_R=0.0,
                sumStdDev=0.0, maxStdDev=0.0, sumPSRF=0.0, maxPSRF=0.0, avgStdDev=0.0, avgPSRF=0.0;
	safeLong		*x;
	char		*s=NULL, tempName[100], tempStr[100], fileName[100], treeName[100], rateName[100], clockFile[100],
				*tempStringP;
	FILE		*fp, *fpParam;
	
#	if defined (MPI_ENABLED)
	if (proc_id == 0)
		{
#	endif

	/* Ensure that we read trees with sumt code and not user tree code */
    inSumtCommand = YES;

	/* Are we comparing two files, or using sumt? */
	comparingFiles = NO;

	/* Initially, we do not want to print brlens to file */
	printingBrlens = NO;

	/* Are we focusing on brlens or node depths ? */
	if (!strcmp(sumtParams.phylogramType, "Nodedepths"))
		nodeDepthConTree = YES;
	else
		nodeDepthConTree = NO;
	
	/* set file pointers to NULL */
	fp = fpParam = fpParts = fpCon = fpTrees = fpBrlens = NULL;

	/* Check that a data set has been read in. We check taxon names against
	   those read in. */
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before sumt can be used\n", spacer);
		goto errorExit;
		}

	/* Check that if there is anything to do */
    if (sumtParams.table == NO && sumtParams.summary == NO && sumtParams.showConsensus == NO)
		{
		MrBayesPrint ("%s   Nothing to do, all output parameters (Table, Summary, Consensus) set to 'NO'\n", spacer);
		goto errorExit;
		}

    for (treeNo = 0; treeNo < sumtParams.numTrees; treeNo++)
		{
		/* initialize across-file tree counter */
		numTotalTreesSampled = 0;

		/* initialize oneUnreliable && unreliable */
		oneUnreliable = unreliable = NO;
		
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
			MrBayesPrint ("%s   Summarizing trees in files \"%s.run1.t\", \"%s.run2.t\", etc\n", spacer, fileName, fileName);

		for (runNo=0; runNo < sumtParams.numRuns; runNo++)
			{
			/* initialize tree counters */
			numSumtTrees = numSumTreesSampled = 0;

			/* open binary file */
			if (sumtParams.numRuns == 1)
				sprintf (tempName, "%s.t", fileName);
			else
				sprintf (tempName, "%s.run%d.t", fileName, runNo+1);

			if ((fp = OpenBinaryFileR(tempName)) == NULL)
				{
				if (strcmp (fileName+strlen(fileName)-2, ".t") == 0)
					{
					MrBayesPrint ("%s   You probably need to remove '.t' from 'Filename'\n", spacer);
					MrBayesPrint ("%s   Also make sure that 'Nruns' and 'Ntrees' are set correctly\n", spacer);
					}
				else
					MrBayesPrint ("%s   Make sure that 'Nruns' and 'Ntrees' are set correctly\n", spacer);
				goto errorExit;
				}
		
			/* find out what type of line termination is used */
			if (runNo == 0 && treeNo == 0)
				{
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
					goto errorExit;
					}
				if (sumtParams.numRuns > 1)
					MrBayesPrint ("%s   Examining first file ...\n", spacer);
				else
					MrBayesPrint ("%s   Examining file ...\n", spacer);
				}
			else if (LineTermType (fp) != lineTerm)
				{
				MrBayesPrint ("%s   Inconsistent line termination for file %d\n", spacer, runNo + 1);
				goto errorExit;
				}

			/* find length of longest line */
			longestLineLength = LongestLine (fp);
			longestLineLength += 10;
	
			/* allocate a string long enough to hold a line */
			if (runNo == 0 && treeNo == 0)
				{
				if (memAllocs[ALLOC_SUMTSTRING] == YES)
					{
					MrBayesPrint ("%s   Sumt string is already allocated\n", spacer);
					goto errorExit;
					}
				s = (char *)SafeMalloc((size_t) (longestLineLength * sizeof(char)));
				if (!s)
					{
					MrBayesPrint ("%s   Problem allocating string for reading sumt file\n", spacer);
					goto errorExit;
					}
				memAllocs[ALLOC_SUMTSTRING] = YES;
				}
			else
				{
				free (s);
				s = (char *) SafeMalloc (sizeof (char) * longestLineLength);
				}
		
			/* close binary file */
			SafeFclose (&fp);
	
			/* open text file */
			if ((fp = OpenTextFileR(tempName)) == NULL)
				goto errorExit;
	
			/* Check file for appropriate blocks. We want to find the last tree block
			   in the file and start from there. */
			foundBegin = inTreeBlock = blockErrors = inSumtComment = NO;
			lineNum = numTreesInBlock = lastTreeBlockBegin = lastTreeBlockEnd = numTreeBlocks = numTreesInLastBlock = 0;
			while (fgets (s, longestLineLength, fp) != NULL)
				{
				sumtTokenP = &s[0];
				do
					{
					GetSumtToken (&tokenType);
					if (IsSame("[", sumtToken) == SAME)
						inSumtComment = YES;
					if (IsSame("]", sumtToken) == SAME)
						inSumtComment = NO;
						
					if (inSumtComment == YES)
						{
						if (IsSame ("Param", sumtToken) == SAME)
							{
							/* extract the tree name */
							GetSumtToken (&tokenType);	/* get the colon */
							GetSumtToken (&tokenType);	/* get the tree name */
							strcpy (treeName, sumtToken);
							GetSumtToken (&tokenType);
							while (IsSame("]", sumtToken) != SAME)
								{
								strcat (treeName, sumtToken);
								GetSumtToken (&tokenType);
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
								lastTreeBlockBegin = lineNum;
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
									numTreeBlocks++;
									inTreeBlock = NO;
									lastTreeBlockEnd = lineNum;
									}
								else
									{
									MrBayesPrint ("%s   Found inappropriate \"End\" statement in file\n", spacer);
									blockErrors = YES;
									}
								numTreesInLastBlock = numTreesInBlock;
								}
							else if (IsSame("Tree", sumtToken) == SAME)
								{
								if (inTreeBlock == YES)
									{
									numTreesInBlock++;
									if (numTreesInBlock == 1)
										{
										sumtBrlensDef = NO;
										for (i=0; s[i]!='\0'; i++)
											{
											if (s[i] == ':')
												{
												sumtBrlensDef = YES;
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

			if (runNo == 0)
				firstSumtBrlensDef = sumtBrlensDef;
			else if (firstSumtBrlensDef != sumtBrlensDef)
				{
				MrBayesPrint ("%s   Tree files with and without brlens mixed\n");
				goto errorExit;
				}
				
			/* Now, check some aspects of the tree file, such as the number of tree blocks and whether they are properly terminated. */
			if (inTreeBlock == YES)
				{
				MrBayesPrint ("%s   Unterminated tree block in file %s. You probably need to\n", spacer, tempName);
				MrBayesPrint ("%s   add a new line to the end of the file with \"End;\" on it.\n", spacer);
				goto errorExit;
				}
			if (inSumtComment == YES)
				{
				MrBayesPrint ("%s   Unterminated comment in file %s\n", spacer, tempName);
				goto errorExit;
				}
			if (blockErrors == YES)
				{
				MrBayesPrint ("%s   Found formatting errors in file %s\n", spacer, tempName);
				goto errorExit;
				}
			if (lastTreeBlockEnd < lastTreeBlockBegin)
				{
				MrBayesPrint ("%s   Problem reading tree file %s\n", spacer, tempName);
				goto errorExit;
				}
			if (numTreesInLastBlock <= 0)
				{
				MrBayesPrint ("%s   No trees were found in last tree block of file %s\n", spacer, tempName);
				goto errorExit;
				}
            if (sumtParams.relativeBurnin == NO && sumtParams.sumtBurnIn > numTreesInLastBlock)
				{
				MrBayesPrint ("%s   No trees are sampled as the burnin exceeds the number of trees in last block\n", spacer);
				MrBayesPrint ("%s   Try setting burnin to a number less than %d\n", spacer, numTreesInLastBlock);
				goto errorExit;
				}
				
			/* tell the user that everything is fine */
			if (runNo == 0)
				{
				if (numTreeBlocks == 1)
					MrBayesPrint ("%s   Found one tree block in file \"%s\" with %d trees in last block\n", spacer, tempName, numTreesInLastBlock);
				else
					{
					MrBayesPrint ("%s   Found %d tree blocks in file \"%s\" with %d trees in last block\n", spacer, numTreeBlocks, tempName, numTreesInLastBlock);
					MrBayesPrint ("%s   Only the %d trees in last tree block will be summarized\n", spacer, numTreesInLastBlock);
					}
				if (sumtParams.numRuns > 1)
					MrBayesPrint ("%s   Expecting the same number of trees in the last tree block of all files\n", spacer);
				firstFileNumTrees = numTreesInLastBlock;
				}
			else
				{
				if (numTreesInLastBlock != firstFileNumTrees)
					{
					MrBayesPrint ("%s   Found %d trees in first file but %d trees in file \"%s\"\n", spacer, firstFileNumTrees, numTreesInLastBlock, tempName);
					goto errorExit;
					}
				}
		
			/* get the tree clock rates if we are focusing on node depths */
			if (nodeDepthConTree == YES)
				{
				strcpy (clockFile, tempName);
				clockFile[strlen(clockFile)-1] = 'p';
				tempStringP = treeName;
				for (tempStringP=treeName; (*tempStringP)!='\0'; tempStringP++)
					if (*tempStringP == '{')
						break;
				strcpy (rateName, "clockRate");
				strcat (rateName, tempStringP);
				if (GetClockRates (rateName, clockFile) == ERROR)
					{
					MrBayesPrint ("%s   Could not get clock rates for the tree\n", spacer);
					goto errorExit;
					}
				treeIndex = 0;	/* reset treeIndex before reading trees */
				}

			/* Now we read the file for real. First, rewind file pointer to beginning of file... */
			(void)fseek(fp, 0L, 0);	
	
			/* ...and fast forward to beginning of last tree block. */
			for (i=0; i<lastTreeBlockBegin+1; i++)
				fgets (s, longestLineLength, fp);

			/* Allocate things we will need for trees... */
			if (treeNo == 0 && runNo == 0)
				{
				if (AllocBits (numTaxa) == ERROR)
					goto errorExit;
				}

			/* Set up cheap status bar. */
			if (runNo == 0)
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
			inTreesBlock = YES;
            ResetTranslateTable();
			for (i=0; i<lastTreeBlockEnd - lastTreeBlockBegin - 1; i++)
				{
				fgets (s, longestLineLength, fp);
				/*MrBayesPrint ("%s", s);*/
				runIndex = runNo;	/* not an elegant solution to use a global (for this file), but we need the run number in DoTree */
				if (ParseCommand (s) == ERROR)
					goto errorExit;
				}
			inTreesBlock = NO;
	
			/* Finish cheap status bar. */
			if (runNo == sumtParams.numRuns - 1)
				{
				if (numAsterices < 80)
					{
					for (i=0; i<80 - numAsterices; i++)
						MrBayesPrint ("*");
					}
				MrBayesPrint ("\n\n");
				}
	
			/* check or print out information on absent and pruned taxa */
			if (runNo == 0 && treeNo == 0 && sumtParams.numRuns > 1)
				{
				for (i=0; i<numTaxa; i++)
					{
					firstPrunedTaxa[i] = prunedTaxa[i];
					firstAbsentTaxa[i] = absentTaxa[i];
					}
				}
			else if (runNo != sumtParams.numRuns - 1 || treeNo != 0)
				{
				for (i=0; i<numTaxa; i++)
					{
					if (prunedTaxa[i] != firstPrunedTaxa[i])
						break;
					if (absentTaxa[i] != firstAbsentTaxa[i])
						break;
					}
				if (i != numTaxa)
					{
					if (runNo != sumtParams.numRuns - 1) 
						MrBayesPrint ("\n\n%s   Mismatch in the absent or pruned taxa\n", spacer);
					else
						MrBayesPrint ("\n\n%s   Mismatch in the absent or pruned taxa\n", spacer);
					goto errorExit;
					}
				}
			else
				{
				/* print out information on absent taxa */
				numExcludedTaxa = 0;
				for (i=0; i<numTaxa; i++)
					if (absentTaxa[i] == YES)
						numExcludedTaxa++;
				if (numExcludedTaxa > 0)
					{
					if (numExcludedTaxa == 1)
						MrBayesPrint ("%s   The following species was absent from trees:\n", spacer);
					else
						MrBayesPrint ("%s   The following %d species were absent from trees:\n", spacer, numExcludedTaxa);
					MrBayesPrint ("%s      ", spacer);
					j = lineWidth = 0;
					for (i=0; i<numTaxa; i++)
						{
						if (absentTaxa[i] == YES)
							{
							j++;
							if (GetNameFromString (taxaNames, tempStr, i+1) == ERROR)
								{
								MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
								return (ERROR);
								}
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
					if (prunedTaxa[i] == YES && absentTaxa[i] == NO)
						numExcludedTaxa++;
				
				if (numExcludedTaxa > 0)
					{
					if (numExcludedTaxa == 1)
						MrBayesPrint ("%s   The following species was pruned from trees:\n", spacer);
					else
						MrBayesPrint ("%s   The following %d species were pruned from trees:\n", spacer, numExcludedTaxa);
					MrBayesPrint ("%s      ", spacer);
					j = lineWidth = 0;
					for (i=0; i<numTaxa; i++)
						{
						if (prunedTaxa[i] == YES && absentTaxa[i] == NO)
							{
							j++;
							if (GetNameFromString (taxaNames, tempStr, i+1) == ERROR)
								{
								MrBayesPrint ("%s   Could not find taxon %d\n", spacer, i+1);
								return (ERROR);
								}
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

			/* Update total number of sampled trees */
			numTotalTreesSampled += numSumTreesSampled;

			/* tell user how many trees were successfully read */
			if (sumtParams.numRuns == 1)
				MrBayesPrint ("%s   Read %d trees from last tree block (sampling %d of them)\n", spacer, numSumtTrees, numSumTreesSampled);
			else if (sumtParams.numRuns > 1)
				{
				if (runNo == 0)
					firstFileNumSumTreesSampled = numSumTreesSampled;
				else if (numSumTreesSampled != firstFileNumSumTreesSampled)
					{
					if (runNo != sumtParams.numRuns - 1)
						MrBayesPrint ("%s   Found %d post-burnin trees in the first file but %d post-burnin trees in file %d\n", spacer,
							firstFileNumSumTreesSampled, numSumTreesSampled, runNo+1);
					else
						MrBayesPrint ("\n\n%s   Found %d post-burnin trees in the first file but %d post-burnin trees in file %d\n", spacer,
							firstFileNumSumTreesSampled, numSumTreesSampled, runNo+1);
					goto errorExit;
					}
				if (runNo == sumtParams.numRuns - 1)
					{
					MrBayesPrint ("%s   Read a total of %d trees in %d files (sampling %d of them)\n", spacer, sumtParams.numRuns * numSumtTrees,
						sumtParams.numRuns, numTotalTreesSampled);
					MrBayesPrint ("%s      (Each file contained %d trees of which %d were sampled)\n", spacer, numSumtTrees,
						numSumTreesSampled);
					}
				}

			/* Check that at least one tree was read in. */
			if (numSumTreesSampled <= 0)
				{
				MrBayesPrint ("%s   No trees read in\n", spacer);
				goto errorExit;
				}

			SafeFclose (&fp);
			}	/* next run for this tree */
				
		/* Sort partitions... */
		if (treeNo == 0)
			{
			if (memAllocs[ALLOC_PARTORIGORDER] == YES)
				{
				MrBayesPrint ("%s   partOrigOrder not free in DoSumt\n", spacer);
				goto errorExit;
				}
			partOrigOrder = (int *)SafeMalloc((size_t) (numTreePartsFound * sizeof(int)));
			if (!partOrigOrder)
				{
				MrBayesPrint ("%s   Problem allocating partOrigOrder (%d)\n", spacer, numTreePartsFound * sizeof(int));
				goto errorExit;
				}
			memAllocs[ALLOC_PARTORIGORDER] = YES;
			}
		else
			{
			free (partOrigOrder);
			partOrigOrder = (int *)SafeMalloc((size_t) (numTreePartsFound * sizeof(int)));
			if (!partOrigOrder)
				{
				MrBayesPrint ("%s   Problem allocating partOrigOrder (%d)\n", spacer, numTreePartsFound * sizeof(int));
				goto errorExit;
				}
			}

		for (i=0; i<numTreePartsFound; i++)
			partOrigOrder[i] = i;
		if (SortParts (numFoundOfThisPart, numTreePartsFound) == ERROR)
			goto errorExit;

		/* ...and reorder bits if some taxa were not included. */
		numIncludedTaxa = numTaxa;
		if (ReorderParts () == ERROR)
			goto errorExit;
		
		/* open output files for summary information (three files) */
		if (OpenSumtFiles (treeNo) == ERROR)
			goto errorExit;

        /* Print partitions to screen. */
		if (treeNo == 0)
			{
			longestName = 0;
			for (k=0; k<numTaxa; k++)
				{
				if (GetNameFromString (taxaNames, tempName, k+1) == ERROR)
					{
					MrBayesPrint ("%s   Error getting taxon names \n", spacer);
					return (ERROR);
					}
				len = (int) strlen (tempName);
				if (len > longestName)
					longestName = len;
				}
            if (sumtParams.table == YES)
                {
                MrBayesPrint ("                                                                                   \n");
			    MrBayesPrint ("%s   General explanation:                                                          \n", spacer);
			    MrBayesPrint ("                                                                                   \n");
			    MrBayesPrint ("%s   A taxon bibartition is specified by removing a branch, thereby dividing the   \n", spacer);
			    MrBayesPrint ("%s   species into those to the left and those to the right of the branch. Here,    \n", spacer);
			    MrBayesPrint ("%s   taxa to one side of the removed branch are denoted \".\" and those to the     \n", spacer);
			    MrBayesPrint ("%s   other side are denoted \"*\". The output includes the bipartition number      \n", spacer);
			    MrBayesPrint ("%s   (ID; sorted from highest to lowest probability), bipartition (e.g., ...**..), \n", spacer);
			    MrBayesPrint ("%s   number of times the bipartition was observed (#obs), the posterior probabil-  \n", spacer);
			    MrBayesPrint ("%s   ity of the bipartition, and, if branch lengths were recorded on the trees in  \n", spacer);
			    MrBayesPrint ("%s   the file, either the average (Mean(v)) and variance (Var(v) of the lengths    \n", spacer);
			    MrBayesPrint ("%s   (if the 'phylogramtype' option is set to 'brlens'), or the average (Mean(d))  \n", spacer);
			    MrBayesPrint ("%s   and variance (Var(d)) of the node depths (if the 'phylogramtype' option is    \n", spacer);
			    MrBayesPrint ("%s   set to 'nodedepths'). The latter option is preferable for clock trees and     \n", spacer);
			    MrBayesPrint ("%s   will take calibration points into account, if available.                      \n", spacer);	
			    MrBayesPrint ("\n");
			    MrBayesPrint ("%s   Each \".\" or \"*\" in the bipartition represents a taxon that is to the left \n", spacer);
			    MrBayesPrint ("%s   or right of the removed branch. A list of the taxa in the bipartition is given\n", spacer);
			    MrBayesPrint ("%s   before the list of bipartitions. If you summarize several independent analy-  \n", spacer);
			    MrBayesPrint ("%s   ses, convergence diagnostics are presented for both the posterior probabil-   \n", spacer);
			    MrBayesPrint ("%s   ities of bipartitions (bipartition or split frequencies) and branch lengths or\n", spacer);
			    MrBayesPrint ("%s   node depths (if recorded on the trees in the files). In the former case, the  \n", spacer);
			    MrBayesPrint ("%s   diagnostic is the standard deviation of the partition frequencies (Stdev(s)), \n", spacer);
			    MrBayesPrint ("%s   in the second case it is the potential scale reduction factor (PSRF) of Gelman\n", spacer);
			    MrBayesPrint ("%s   and Rubin (1992). Stdev(s) is expected to approach 0 and PSRF is expected to  \n", spacer);
			    MrBayesPrint ("%s   approach 1 as runs converge onto the posterior probability distribution. Note \n", spacer);
			    MrBayesPrint ("%s   that these values may be unreliable if the partition is not present in all    \n", spacer);
			    MrBayesPrint ("%s   runs (the last column indicates the number of runs that sampled the partition \n", spacer);
			    MrBayesPrint ("%s   if more than one run is summarized). The PSRF is also sensitive to small      \n", spacer);
			    MrBayesPrint ("%s   sample sizes and it should only be considered a rough guide to convergence    \n", spacer);
			    MrBayesPrint ("%s   since some of the assumptions allowing one to interpret it as a true potential\n", spacer);
			    MrBayesPrint ("%s   scale reduction factor are violated in the phylogenetic context.              \n", spacer);
			    MrBayesPrint ("%s                                                                                 \n", spacer);
			    MrBayesPrint ("%s   List of taxa in bipartitions:                                                 \n", spacer);
			    MrBayesPrint ("                                                                                   \n");
			    j = 1;
			    for (k=0; k<numTaxa; k++)
				    {
				    if (GetNameFromString (taxaNames, tempName, k+1) == ERROR)
					    {
					    MrBayesPrint ("%s   Error getting taxon names \n", spacer);
					    return (ERROR);
					    }
				    if (sumTaxaFound[k] == YES)
					    {
					    MrBayesPrint ("%s   %4d -- %s\n", spacer, j++, tempName);
					    }
				    }
                }
			}
	
        if (sumtParams.numTrees > 1 && (sumtParams.table == YES || sumtParams.summary == YES))
			{
			MrBayesPrint ("%s   Results for tree number %d\n", spacer, treeNo+1);
			MrBayesPrint ("%s   ==========================\n\n", spacer);
			}

		if (sumtParams.table == YES)
            {
            MrBayesPrint ("                                                                                   \n");
		    MrBayesPrint ("%s   Summary statistics for taxon bipartitions:                                \n\n", spacer);
            }

		/* print header with legend to .parts file */
		if (treeNo == 0 && sumtParams.table == YES)
            {
            MrBayesPrintf (fpParts, "[   ID        = Partition ID number]\n");
		    MrBayesPrintf (fpParts, "[   Partition = Description of partition in .* format]\n");
            MrBayesPrintf (fpParts, "[   #obs      = Number of trees sampled with the partition]\n");
		    MrBayesPrintf (fpParts, "[   Probab.   = Posterior probability of the partition]\n");
		    if (sumtParams.numRuns > 1)
		    	{
			    MrBayesPrintf (fpParts, "[   Stddev(s) = Standard deviation of partition probabilities across partitions]\n");
			    }
		    if (sumtBrlensDef == YES)
			    {
                if (nodeDepthConTree == YES)
                    {
			        MrBayesPrintf (fpParts, "[   Mean(v)   = Mean branch length]\n");
			        MrBayesPrintf (fpParts, "[   Var(v)    = Variance of branch length]\n");
			        if (sumtParams.numRuns > 1)
				        {
				        MrBayesPrintf (fpParts, "[   PSRF   = Potential scale reduction factor of branch length]\n");
				        }
                    }
                else
                    {
			        MrBayesPrintf (fpParts, "[   Mean(d)   = Mean node depth]\n");
			        MrBayesPrintf (fpParts, "[   Var(d)    = Variance of node depth]\n");
			        if (sumtParams.numRuns > 1)
				        {
				        MrBayesPrintf (fpParts, "[   PSRF   = Potential scale reduction factor of node depth]\n");
				        }
                    }
			    }
		    if (sumtParams.numRuns > 1)
			    {
			    MrBayesPrintf (fpParts, "[   Nruns  = Number of runs with this partition]\n");
			    }
            MrBayesPrintf (fpParts, "\n\n");
            }

		/* calculate a couple of numbers that are handy to have */
		numTreePartsToPrint = 0;
		for (i=0; i<numTreePartsFound; i++)
			{
			if ((MrBFlt)numFoundOfThisPart[i]/(MrBFlt)numTotalTreesSampled >= sumtParams.freqDisplay)
				numTreePartsToPrint++;
			}
		maxWidthID = (int) (log10 (numTreePartsToPrint)) + 1;
		if (maxWidthID < 2)
			maxWidthID = 2;
		maxWidthNumberPartitions = (int) (log10 (numFoundOfThisPart[0])) + 1;
		if (maxWidthNumberPartitions < 4)
			maxWidthNumberPartitions = 4;

		/* print header to screen and to parts file simultaneously */
		if (sumtParams.table == YES)
            {
            MrBayesPrint ("%s   ", spacer);
		    MrBayesPrintf (fpParts, "%s   ", spacer);
		    MrBayesPrint ("%*s -- Partition", maxWidthID, "ID");
		    MrBayesPrintf (fpParts, "%*s -- Partition", maxWidthID, "ID");
		    tableWidth = maxWidthID + 13;
		    for (i=9; i<numIncludedTaxa; i++)
			    {
			    MrBayesPrint (" ");
			    MrBayesPrintf (fpParts, " ");
			    tableWidth++;
			    }
		    MrBayesPrint ("  #obs");
		    MrBayesPrintf (fpParts, "  #obs");
		    tableWidth += 6;
		    for (i=4; i<maxWidthNumberPartitions; i++)
			    {
			    MrBayesPrint (" ");
			    MrBayesPrintf (fpParts, " ");
			    tableWidth++;
			    }
		    MrBayesPrint ("   Probab.");
		    MrBayesPrintf (fpParts, "   Probab.");
		    tableWidth += 10;
		    if (sumtParams.numRuns > 1)
			    {
			    MrBayesPrint (" Stddev(s)");
			    MrBayesPrintf (fpParts, " Stddev(s)");
			    tableWidth += 10;
			    }
		    if (sumtBrlensDef == YES)
			    {
			    if (nodeDepthConTree == YES)
                    {
				    MrBayesPrint ("   Mean(d)   Var(d) ");
				    MrBayesPrintf (fpParts, "   Mean(d)   Var(d) ");
                    }
			    else
                    {
				    MrBayesPrint ("   Mean(v)   Var(v) ");
				    MrBayesPrintf (fpParts, "   Mean(v)   Var(v) ");
                    }
			    tableWidth += 20;
			    if (sumtParams.numRuns > 1)
				    {
				    MrBayesPrint ("   PSRF");
				    MrBayesPrintf (fpParts, "   PSRF");
				    tableWidth += 7;
				    }
			    }
		    if (sumtParams.numRuns > 1)
			    {
			    MrBayesPrint ("  Nruns");
			    MrBayesPrintf (fpParts, "  Nruns");
			    tableWidth += 7;
			    }
		    MrBayesPrint ("\n%s   ", spacer);
		    MrBayesPrintf (fpParts, "\n%s   ", spacer);
		    for (i=0; i<tableWidth; i++)
                {
			    MrBayesPrint ("-");
			    MrBayesPrintf (fpParts, "-");
                }
		    MrBayesPrint ("\n");
		    MrBayesPrintf (fpParts, "\n");
            }

        /* now, show partitions that were found on screen; print to .parts file simultaneously */
		x = &treePartsFound[0];
		for (i=0; i<numTreePartsToPrint; i++)
			{
			if (sumtParams.table == YES)
                {
                MrBayesPrint ("%s   %*d -- ", spacer, maxWidthID, i+1);
			    MrBayesPrintf (fpParts, "%s   %*d -- ", spacer, maxWidthID, i+1);
			    ShowParts (stdout, &x[0], numIncludedTaxa);
			    fflush(stdout);
                ShowParts (fpParts, &x[0], numIncludedTaxa);
			    for (j=numIncludedTaxa; j<9; j++)
                    {
				    MrBayesPrint (" ");
				    MrBayesPrintf (fpParts, " ");
                    }
                }
			if (sumtBrlensDef == YES)
				{
				if (numFoundOfThisPart[i] == 1)
					var = 0.0;
				else
					var = sBrlens[i] / (numFoundOfThisPart[i]-1);
				}
			if (sumtParams.numRuns > 1)
				{
				sum_s = 0.0;
				sumsq_s = 0.0;
				for (n=j=0; n<sumtParams.numRuns; n++)
					{
					if (numFoundInRunOfPart[n][i] > 0)
						j++;
					f = (MrBFlt) numFoundInRunOfPart[n][i] / (MrBFlt) firstFileNumSumTreesSampled;
                    sum_s += f;
					sumsq_s += f * f;
					}
				var_s = sumsq_s - sum_s * sum_s / (MrBFlt) sumtParams.numRuns;
				var_s /= (sumtParams.numRuns - 1);
				if (var_s > 0.0)
					stddev_s = sqrt (var_s);
				else
					stddev_s = 0.0;
				if (sumtBrlensDef == YES)
					{
					if (numFoundOfThisPart[i] - j == 0)
						varW = 0.0;
					else
						varW = sWithinBrlens[i] / (MrBFlt) (numFoundOfThisPart[i] - j);
					if (j == 1)
						varB = 0.0;
					else
						varB = (sumsqB[i] - sumB[i]*sumB[i]/(MrBFlt) j) / (MrBFlt) (j - 1);
					if (varW > 0.0)
						sqrt_R = sqrt ((MrBFlt)(firstFileNumSumTreesSampled - 1) / (MrBFlt) (firstFileNumSumTreesSampled)
							+ ((MrBFlt) (j + 1) / (MrBFlt) (j)) * (varB / varW));
                    else
                        sqrt_R = -1.0;      /* we are not going to write this value anyway but avoid math errors */
					}
				if (j == sumtParams.numRuns)
					unreliable = NO;
				else
					{
					unreliable = YES;
					oneUnreliable = YES;
					}
				}
			f = (MrBFlt) numFoundOfThisPart[i] / (MrBFlt) numTotalTreesSampled;
			if (sumtParams.table == YES)
                {
                MrBayesPrint ("  %*d  %1.6lf", maxWidthNumberPartitions, numFoundOfThisPart[i], f);
			    MrBayesPrintf (fpParts, "  %*d  %1.6lf", maxWidthNumberPartitions, numFoundOfThisPart[i], f);
			    if (sumtParams.numRuns > 1)
                    {
				    MrBayesPrint ("  %1.6lf", stddev_s);
				    MrBayesPrintf (fpParts, "  %1.6lf", stddev_s);
                    }
			    if (sumtBrlensDef == YES)
				    {
				    if (nodeDepthConTree == YES)
					    {
					    precision = 6 - (int) log10 (sumtRoot->left->age);
					    if (precision < 0)
						    precision = 0;
					    else if (precision > 3)
						    precision = 3;
					    MrBayesPrint ("  %8.*lf  %8.*lf", precision, aBrlens[i], precision, var);
					    MrBayesPrintf (fpParts, "  %8.*lf  %8.*lf", precision, aBrlens[i], precision, var);
					    }
				    else
                        {
					    MrBayesPrint ("  %1.6lf  %1.6lf", aBrlens[i], var);
					    MrBayesPrintf (fpParts, "  %1.6lf  %1.6lf", aBrlens[i], var);
                        }
				    if (sumtParams.numRuns > 1)
					    {
					    if (varW <= 0.0 && varB > 0.0)
                            {
						    MrBayesPrint ("  > 2.0");
						    MrBayesPrintf (fpParts, "  > 2.0");
                            }
					    else if ((varW <= 0.0 && varB <= 0.0) || j == 1)
                            {
						    MrBayesPrint ("   N/A ");
						    MrBayesPrintf (fpParts, "   N/A ");
                            }
					    else
                            {
						    MrBayesPrint ("  %1.3lf", sqrt_R);
						    MrBayesPrintf (fpParts, "  %1.3lf", sqrt_R);
                            }
					    }
				    }
			    if (sumtParams.numRuns > 1)
				    {
				    MrBayesPrint ("  %3d", j);
				    MrBayesPrintf (fpParts, "  %3d", j);
			    	}
			    if (unreliable == YES)
                    {
				    MrBayesPrint (" *\n");
				    MrBayesPrintf (fpParts, " *\n");
                    }
			    else
                    {
				    MrBayesPrint ("\n");
				    MrBayesPrintf (fpParts, "\n");
                    }
                }
			x += taxonLongsNeeded;
			}
	    if (sumtParams.table == YES)
            {
            MrBayesPrint ("%s   ", spacer);
	        MrBayesPrintf (fpParts, "%s   ", spacer);
	        for (i=0; i<tableWidth; i++)
                {
		        MrBayesPrint ("-");
		        MrBayesPrintf (fpParts, "-");
                }
	        MrBayesPrint ("\n");
	        MrBayesPrintf (fpParts, "\n");
	        if (oneUnreliable == YES)
                {
		        MrBayesPrint ("%s   * The partition was not found in all runs so the values are unreliable\n\n", spacer);
		        MrBayesPrintf (fpParts, "%s   * The partition was not found in all runs so the values are unreliable\n\n", spacer);
                }
		    else
                {
			    MrBayesPrint ("\n");
			    MrBayesPrintf (fpParts, "\n");
                }
            }

		/* calculate summary statistics */
		if (sumtParams.summary == YES && sumtParams.numRuns > 1)
            {
            x = &treePartsFound[0];
		    nSummarySamples = 0;
            sumStdDev = 0.0;
            maxStdDev = 0.0;
            sumPSRF = 0.0;
            maxPSRF = 0.0;
            for (i=0; i<numTreePartsFound; i++)
                {
				sum_s = 0.0;
				sumsq_s = 0.0;
				includeInStat = NO;
                for (n=j=0; n<sumtParams.numRuns; n++)
					{
					if (numFoundInRunOfPart[n][i] > 0)
						j++;
					f = (MrBFlt) numFoundInRunOfPart[n][i] / (MrBFlt) firstFileNumSumTreesSampled;
                    if (f > sumtParams.minPartFreq)
                        includeInStat = YES;
                    sum_s += f;
					sumsq_s += f * f;
					}
				var_s = sumsq_s - sum_s * sum_s / (MrBFlt) sumtParams.numRuns;
				var_s /= (sumtParams.numRuns - 1);
				if (var_s > 0.0)
					stddev_s = sqrt (var_s);
				else
					stddev_s = 0.0;
                if (includeInStat == YES)
                    {
                    nSummarySamples++;
                    sumStdDev += stddev_s;
                    if (stddev_s > maxStdDev)
                        maxStdDev = stddev_s;
				    if (sumtBrlensDef == YES)
					    {
					    if (numFoundOfThisPart[i] - j == 0)
						    varW = 0.0;
					    else
						    varW = sWithinBrlens[i] / (MrBFlt) (numFoundOfThisPart[i] - j);
					    if (j == 1)
						    varB = 0.0;
					    else
						    varB = (sumsqB[i] - sumB[i]*sumB[i]/(MrBFlt) j) / (MrBFlt) (j - 1);
					    if (varW > 0.0)
						    sqrt_R = sqrt ((MrBFlt)(firstFileNumSumTreesSampled - 1) / (MrBFlt) (firstFileNumSumTreesSampled)
							    + ((MrBFlt) (j + 1) / (MrBFlt) (j)) * (varB / varW));
					    else
						    sqrt_R = -1.0;
                        if (sqrt_R > maxPSRF && j > 1)
                            maxPSRF = sqrt_R;
                        if (sqrt_R > 0.0 && j > 1)
                            {
                            sumPSRF += sqrt_R;
                            nPSRFSamples++;
                            }
					    }
                    }
                x += taxonLongsNeeded;
                }
            
            /* Exclude trivial splits when calculating average standard deviation of split frequencies.
               Note that when node depths are requested, the depth of the entire tree is also a trivial
               split included in the table. */
            if (nodeDepthConTree == YES)
                avgStdDev = sumStdDev / (nSummarySamples-numIncludedTaxa-1);
            else
                avgStdDev = sumStdDev / (nSummarySamples-numIncludedTaxa);
            avgPSRF   = sumPSRF / nPSRFSamples;

            MrBayesPrint ("%s   Summary statistics for partitions with frequency > %1.2lf in at least one run:\n", spacer, sumtParams.minPartFreq);
            MrBayesPrint ("%s       Average standard deviation of split frequencies = %1.6lf\n", spacer, avgStdDev);
            MrBayesPrint ("%s       Maximum standard deviation of split frequencies = %1.6lf\n", spacer, maxStdDev);
            if (sumtBrlensDef == YES)
                {
                MrBayesPrint ("%s       Average potential scale reduction factor = %1.3lf\n", spacer, avgPSRF);
                MrBayesPrint ("%s       Maximum potential scale reduction factor = %1.3lf\n", spacer, maxPSRF);
                }
            MrBayesPrint ("\n");

            MrBayesPrintf (fpParts, "%s   Summary statistics for partitions with frequency > %1.3lf in at least one run:\n", spacer, sumtParams.minPartFreq);
            MrBayesPrintf (fpParts, "%s       Average standard deviation of split frequencies = %1.3lf\n", spacer, avgStdDev);
            MrBayesPrintf (fpParts, "%s       Maximum standard deviation of split frequencies = %1.3lf\n", spacer, maxStdDev);
            if (sumtBrlensDef == YES)
                {
                MrBayesPrintf (fpParts, "%s       Average potential scale reduction factor = %1.3lf\n", spacer, avgPSRF);
                MrBayesPrintf (fpParts, "%s       Maximum potential scale reduction factor = %1.3lf\n", spacer, maxPSRF);
                }
            }
            MrBayesPrintf (fpParts, "\n");

        /* get branch lengths and print to file if appropriate */
		if (sumtBrlensDef == YES && sumtParams.printBrlensToFile == YES)
			{
			if (BrlensVals (treeNo, s, longestLineLength, lastTreeBlockBegin, lastTreeBlockEnd) == ERROR)
				goto errorExit;
			}

		/* make the majority rule consensus tree */
        if (sumtParams.showConsensus == YES && ConTree () == ERROR)
			goto errorExit;
			
		/* get probabilities of individual trees */
		if (TreeProb () == ERROR)
			goto errorExit;
		
		SafeFclose (&fpParts);
		SafeFclose (&fpCon);
		SafeFclose (&fpTrees);
		SafeFclose (&fpBrlens);
		} /* next tree */

	/* free memory and file pointers */
	if (memAllocs[ALLOC_SUMTSTRING] == YES)
		{
		free (s);
		memAllocs[ALLOC_SUMTSTRING] = NO;
		}
	if (memAllocs[ALLOC_PARTORIGORDER] == YES)
		{
		free (partOrigOrder);
		memAllocs[ALLOC_PARTORIGORDER] = NO;
		}
	FreeBits ();
	expecting = Expecting(COMMAND);
	nodeDepthConTree = NO;

#	if defined (MPI_ENABLED)
		}
#	endif

	inSumtCommand = NO;
	return (NO_ERROR);
	
	/* error exit */			
	errorExit:
		expecting = Expecting(COMMAND);
		if (memAllocs[ALLOC_SUMTSTRING] == YES)
			{
			free (s);
			memAllocs[ALLOC_SUMTSTRING] = NO;
			}
		if (memAllocs[ALLOC_PARTORIGORDER] == YES)
			{
			free (partOrigOrder);
			memAllocs[ALLOC_PARTORIGORDER] = NO;
			}
		FreeBits ();
		SafeFclose (&fp);
		SafeFclose (&fpParts);
		SafeFclose (&fpCon);
		SafeFclose (&fpTrees);
		SafeFclose (&fpBrlens);
		strcpy (spacer, "");
		strcpy (sumtToken, "Sumt");
		i = 0;
		if (FindValidCommand (sumtToken, &i) == ERROR)
			MrBayesPrint ("%s   Could not find sumt\n", spacer);
		nodeDepthConTree = NO;
		inSumtCommand = NO;
		return (ERROR);	
	
}





int DoSumtParm (char *parmName, char *tkn)

{

	int			tempI;
	MrBFlt		tempD;
	char		tempStr[100];
	
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
				strcpy (sumtParams.sumtFileName, tkn);
				MrBayesPrint ("%s   Setting sumt filename to %s\n", spacer, sumtParams.sumtFileName);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Pfile (sumtParams.pFile) ***************************************************/
		else if (!strcmp(parmName, "Pfile"))
			{
			if (expecting == Expecting(EQUALSIGN))
				{
				expecting = Expecting(ALPHA);
				readWord = YES;
				}
			else if (expecting == Expecting(ALPHA))
				{
				strcpy (sumtParams.pFile, tkn);
				MrBayesPrint ("%s   Setting sumt pfile to %s\n", spacer, sumtParams.pFile);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Displaygeq (sumtParams.freqDisplay) *******************************************************/
		else if (!strcmp(parmName, "Displaygeq"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);
				sumtParams.freqDisplay = tempD;
				MrBayesPrint ("%s   Showing partitions with probability greater than or equal to %lf\n", spacer, sumtParams.freqDisplay);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Relburnin (sumtParams.relativeBurnin) ********************************************************/
		else if (!strcmp(parmName, "Relburnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumtParams.relativeBurnin = YES;
					else
						sumtParams.relativeBurnin = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for Relburnin\n", spacer);
					free(tempStr);
					return (ERROR);
					}
				if (sumtParams.relativeBurnin == YES)
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
		/* set Burnin (sumtParams.sumtBurnIn) ***********************************************************/
		else if (!strcmp(parmName, "Burnin"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%d", &tempI);
                sumtParams.sumtBurnIn = tempI;
				MrBayesPrint ("%s   Setting sumt burn-in to %d\n", spacer, sumtParams.sumtBurnIn);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				{
				free(tempStr);
				return (ERROR);
				}
			}
		/* set Burninfrac (sumtParams.sumtBurnInFraction) ************************************************************/
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
                sumtParams.sumtBurnInFraction = tempD;
				MrBayesPrint ("%s   Setting burnin fraction to %.2f\n", spacer, sumtParams.sumtBurnInFraction);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else 
				{
				free(tempStr);
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
					MrBayesPrint ("%s   Setting sumt nruns to %ld\n", spacer, sumtParams.numRuns);
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
					MrBayesPrint ("%s   Setting sumt ntrees to %ld\n", spacer, sumtParams.numTrees);
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
		/* set Phylogramtype (sumtParams.phylogramType) **************************************************/
		else if (!strcmp(parmName, "Phylogramtype"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					strcpy (sumtParams.phylogramType, tempStr);
					MrBayesPrint ("%s   Setting sumt phylogramtype to '%s'\n", spacer, sumtParams.phylogramType);
					}
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Calctrprobs (sumtParams.calcTrprobs) *********************************************/
		else if (!strcmp(parmName, "Calctrprobs"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumtParams.calcTrprobs = YES;
					else
						sumtParams.calcTrprobs = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for calctrprobs\n", spacer);
					return (ERROR);
					}
				if (sumtParams.calcTrprobs == YES)
					MrBayesPrint ("%s   Setting calctrprobs to yes\n", spacer);
				else
					MrBayesPrint ("%s   Setting calctrprobs to no\n", spacer);
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
		/* set Printbrlens (sumtParams.printBrlensToFile) *********************************************/
		else if (!strcmp(parmName, "Printbrlens"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(ALPHA);
			else if (expecting == Expecting(ALPHA))
				{
				if (IsArgValid(tkn, tempStr) == NO_ERROR)
					{
					if (!strcmp(tempStr, "Yes"))
						sumtParams.printBrlensToFile = YES;
					else
						sumtParams.printBrlensToFile = NO;
					}
				else
					{
					MrBayesPrint ("%s   Invalid argument for printbrlens\n", spacer);
					return (ERROR);
					}
				if (sumtParams.printBrlensToFile == YES)
					MrBayesPrint ("%s   Setting printbrlens to yes\n", spacer);
				else
					MrBayesPrint ("%s   Setting printbrlens to no\n", spacer);
				expecting = Expecting(PARAMETER) | Expecting(SEMICOLON);
				}
			else
				return (ERROR);
			}
		/* set Brlensgeq (sumtParams.brlensFreqDisplay) *******************************************************/
		else if (!strcmp(parmName, "Brlensgeq"))
			{
			if (expecting == Expecting(EQUALSIGN))
				expecting = Expecting(NUMBER);
			else if (expecting == Expecting(NUMBER))
				{
				sscanf (tkn, "%lf", &tempD);
				sumtParams.brlensFreqDisplay = tempD;
				MrBayesPrint ("%s   Printing branch lengths to file for partitions with probability >= %lf\n", spacer, sumtParams.brlensFreqDisplay);
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
				if (sumpParams.plot == YES)
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
				if (sumpParams.plot == YES)
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





int DoTranslate (void)

{

	if (inTreesBlock == NO)
		{
		MrBayesPrint ("%s   You must be in a trees block to read a translate command\n", spacer);
		return (ERROR);
		}
	numTranslates++;
	isTranslateDef = YES;
		
#	if 0
	MrBayesPrint ("%s   Defining a translation table\n", spacer);
	
	/* print out translate table */
	{
	int		i, j, len, longestLen;
	char		tempName[100];
	/*MrBayesPrint ("%s\n", transFrom);
	MrBayesPrint ("%s\n", transTo);*/
	longestLen = 0;
	for (i=0; i<numTranslates; i++)
		{
		if (GetNameFromString (transFrom, tempName, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Error getting translate names \n", spacer);
			return (ERROR);
			}
		len = strlen(tempName);
		if (len > longestLen)
			longestLen = len;
		}	
	
	for (i=0; i<numTranslates; i++)
		{
		if (GetNameFromString (transFrom, tempName, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Error getting translate names \n", spacer);
			return (ERROR);
			}
		len = strlen(tempName);		
		MrBayesPrint ("%s      %s", spacer, tempName);
		for (j=0; j<longestLen - len; j++)
			MrBayesPrint (" ");
		MrBayesPrint (" is translated to ");
		if (GetNameFromString (transTo, tempName, i+1) == ERROR)
			{
			MrBayesPrint ("%s   Error getting translate names \n", spacer);
			return (ERROR);
			}
		MrBayesPrint ("%s\n", tempName);
		}
	}
#	endif

	return (NO_ERROR);

}





int DoTranslateParm (char *parmName, char *tkn)

{

	int			howMany;

	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before a translate command can be used\n", spacer);
		return (ERROR);
		}
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
		
	if (expecting == Expecting(ALPHA))
		{
		if (numTranslates == numTaxa)
			{
			MrBayesPrint ("%s   Too many entries in translation table\n", spacer);
			return (ERROR);
			}
		if (whichTranslate == 0)
			{
			if (CheckString (tkn, transTo, &howMany) == ERROR)
				{
				if (AddToString (tkn, transTo, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
					return (ERROR);
					}
				if (howMany - 1 != numTranslates)
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
			if (CheckString (tkn, transFrom, &howMany) == ERROR)
				{
				if (AddToString (tkn, transFrom, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
					return (ERROR);
					}
				if (howMany - 1 != numTranslates)
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
	else if (expecting == Expecting(NUMBER))
		{
		if (numTranslates == numTaxa)
			{
			MrBayesPrint ("%s   Too many entries in translation table\n", spacer);
			return (ERROR);
			}
		if (whichTranslate == 0)
			{
			if (CheckString (tkn, transTo, &howMany) == ERROR)
				{
				if (AddToString (tkn, transTo, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
					return (ERROR);
					}
				if (howMany - 1 != numTranslates)
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
			if (CheckString (tkn, transFrom, &howMany) == ERROR)
				{
				if (AddToString (tkn, transFrom, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Problem adding taxon %s to list\n", spacer, tkn);
					return (ERROR);
					}
				if (howMany - 1 != numTranslates)
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





int DoSumtTree (void)

{

	int			i, z, printEvery, nAstPerPrint, numTreesInThisFile, burnin;
	MrBFlt		x, y;
	
	/* check that we are in a trees block */
	if (inTreesBlock == NO)
		{
		MrBayesPrint ("%s   You must be in a trees block to read a tree\n", spacer);
		goto errorExit;
		}
		
	/* check that the tree is not too big */
	if (numSumtTaxa > numTaxa)
		{
		MrBayesPrint ("%s   Too many taxa in tree\n", spacer);
		goto errorExit;
		}
	
	/* check that the tree is rooted if we calculate dated consensus tree */
	if (nodeDepthConTree == YES && isSumtTreeRooted == NO)
		{
		MrBayesPrint ("%s   Tree is not rooted\n", spacer);
		return (ERROR);
		}

	/* Increment number of trees read in. */
	numSumtTrees++;
	if (comparingFiles == YES)
		numCompTrees[fileNum]++;

	/*  update status bar */
	if (comparingFiles == YES)
		numTreesInThisFile = numCompTrees[fileNum];
	else
		numTreesInThisFile = numSumtTrees;
	if (numTreesInLastBlock * sumtParams.numRuns < 80)
		{
		printEvery = 1;
		nAstPerPrint = 80 / (numTreesInLastBlock * sumtParams.numRuns);
		if (numSumtTrees % printEvery == 0)
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
		x = (MrBFlt)(numTreesInLastBlock * sumtParams.numRuns) / (MrBFlt) (80);
		y = (MrBFlt)(numTreesInThisFile + numTreesInLastBlock * runIndex) / x;
		z = (int)y;
		if (numAsterices < z)
			{
			MrBayesPrint ("*");
			numAsterices++;
			}
		}
	
	/* prepare tree */
    if (sumtParams.relativeBurnin == NO)
        burnin = sumtParams.sumtBurnIn;
    else
        burnin = (int) (numTreesInLastBlock * sumtParams.sumtBurnInFraction);

	if ((comparingFiles == NO && numSumtTrees > burnin) || (comparingFiles == YES && numCompTrees[fileNum] > comptreeParams.comptBurnIn))
		{
		/* Increment the number of trees sampled. */
		numSumTreesSampled++;
		if (comparingFiles == YES)
			numCompTreesSampled[fileNum]++;

		/* Add a node to the root of the tree if
		   this is a rooted tree. */
		if (isSumtTreeRooted == YES)
			{
			if (pSumtPtr->anc == NULL)
				{
				qSumtPtr = &sumtNodes[nextAvailableSumtNode];
				nextAvailableSumtNode++;
				qSumtPtr->left = pSumtPtr;
				pSumtPtr->anc = qSumtPtr;
				pSumtPtr = qSumtPtr;
				sumtRoot = pSumtPtr;
				}
			else
				{
				MrBayesPrint ("%s   Tree is not rooted correctly\n", spacer);
				goto errorExit;
				}
			}
			
		/* Check to see if all of the species in the character matrix (already
		   read in) are present in the trees. We only do this for the first
		   tree read in. */
		if (numSumTreesSampled == 1)
			{
			if (CheckSumtSpecies () == ERROR)
				goto errorExit;
			}
			
		/* Prune deleted taxa from list of trees, if necessary */
		if (PruneSumt () == ERROR)
			goto errorExit;
		if (numSumTreesSampled == 1)
			{
			for (i=0; i<numTaxa; i++)
				{
				if (taxaInfo[i].isDeleted == YES)
					{
					sumTaxaFound[i] = NO;
					}
				}
			}		

		/* Reroot tree, if necessary, on the outgroup or on 
		   the first undeleted taxon. We only do this if the
		   tree is unrooted. If the tree is rooted, then we
		   assume that the root is OK. This seems like extra
		   work, but it is very convenient if all of the trees
		   that are read in are consistently rooted. */
		if (isSumtTreeRooted == NO)
			{
			if (RootSumtTree (sumtRoot, numSumtTaxa, outGroupNum) == ERROR)
				goto errorExit;
			if (DerootSumtTree (sumtRoot, numSumtTaxa, outGroupNum) == ERROR)
				goto errorExit;
			}	
			
		/* Relabel some interior nodes */	
		i = numTaxa;
		if (isSumtTreeRooted == YES)
			FinishSumtTree (sumtRoot, &i, YES);
		else
			FinishSumtTree (sumtRoot, &i, NO);
		
		/* reset brlens if printing brlen vals to file */
		if (printingBrlens == YES)
			{
			for (i=0; i<numBrlens; i++)
				brlens[i] = -1.0;
			}

		/* get partitions for tree */
		if (GetPartitions () == ERROR)
			goto errorExit;
			
		/* find the partitions in the table and increment
		   the appropriate ones */
		if (FindParts (nodeDepthConTree) == ERROR)
			goto errorExit;
			
		/* print brlens to file and bail out if appropriate */
		if (printingBrlens == YES)
			{
			if (PrintBrlensToFile () == ERROR)
				return ERROR;
			else
				return NO_ERROR;
			}

		/* find the tree in the list of trees */
		if (FindTree () == ERROR)
			goto errorExit;
			
		/* add the tree to the list of trees */
		if (comparingFiles == YES)
			{
			if (AddTreeToList (fileNum) == ERROR)
				goto errorExit;
			}

		/* Display the tree nodes. */
#		if 0
		if (isSumtTreeRooted == YES)
			MrBayesPrint ("%s   Rooted tree %d:\n", spacer, numSumtTrees);
		else
			MrBayesPrint ("%s   Unrooted tree %d:\n", spacer, numSumtTrees);
		ShowSumtNodes (sumtRoot, 3, isSumtTreeRooted);
#		endif
		}
	
	return (NO_ERROR);

	errorExit:
		return (ERROR);
	
}





int DoSumtTreeParm (char *parmName, char *tkn)

{

	int			i, tempInt, howMany;
	MrBFlt		tempD;
	char		tempName[100];
	
	if (defMatrix == NO)
		{
		MrBayesPrint ("%s   A matrix must be specified before a tree can be read in\n", spacer);
		goto errorExit;
		}
	if (inTreesBlock == NO)
		{
		MrBayesPrint ("%s   You must be in a trees block to read a tree\n", spacer);
		goto errorExit;
		}
	
	if (expecting == Expecting(PARAMETER))
		{
		/* this should be the name of the tree, but we don't need to do anything with it */
		expecting = Expecting(EQUALSIGN);
		}
	else if (expecting == Expecting(EQUALSIGN))
		{
		for (i=0; i<2*numTaxa; i++)
			{
			sumtNodes[i].left = sumtNodes[i].right = sumtNodes[i].anc = NULL;
			sumtNodes[i].memoryIndex = i;
			sumtNodes[i].length = 0.0;
			sumtNodes[i].marked = NO;
			sumtNodes[i].index = 0;
			sumtNodes[i].taxonName = NO;
			}
		sumtRoot = NULL;
		isSumtTreeDefined = NO;
		isFirstSumtNode = YES;
		pSumtPtr = qSumtPtr = &sumtNodes[0];
		nextAvailableSumtNode = 0;
		numSumtTaxa = 0;
		foundSumtColon = NO;
		isSumtTreeRooted = YES;
		for (i=0; i<numTaxa; i++)
			tempSet[i] = NO;
		expecting  = Expecting(LEFTPAR);
		}
	else if (expecting == Expecting(LEFTPAR))
		{
		if (isFirstSumtNode == YES)
			{
			pSumtPtr = &sumtNodes[nextAvailableSumtNode];
			nextAvailableSumtNode++;
			isFirstSumtNode = NO;
			sumtRoot = pSumtPtr;
			}
		else
			{
			if (nextAvailableSumtNode+1 >= 2*numTaxa)
				{
				MrBayesPrint ("%s   Too many nodes on sumt tree\n", spacer);
				goto errorExit;
				}
			if (pSumtPtr->left == NULL)
				{
				pSumtPtr = &sumtNodes[nextAvailableSumtNode];
				nextAvailableSumtNode++;
				qSumtPtr->left = pSumtPtr;
				pSumtPtr->anc = qSumtPtr;
				qSumtPtr = pSumtPtr;
				}
			else if (pSumtPtr->right == NULL)
				{
				pSumtPtr = &sumtNodes[nextAvailableSumtNode];
				nextAvailableSumtNode++;
				qSumtPtr->right = pSumtPtr;
				pSumtPtr->anc = qSumtPtr;
				qSumtPtr = pSumtPtr;
				}
			else if (pSumtPtr->anc == NULL)
				{
				pSumtPtr = &sumtNodes[nextAvailableSumtNode];
				nextAvailableSumtNode++;
				qSumtPtr->anc = pSumtPtr;
				pSumtPtr->left = qSumtPtr;
				qSumtPtr = pSumtPtr;
				sumtRoot = pSumtPtr;
				pSumtPtr->marked = YES;
				isSumtTreeRooted = NO;
				}
			else
				{
				MrBayesPrint ("\n   ERROR: Tree is not bifurcating\n");
				goto errorExit;
				}
			}
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(LEFTPAR);
		}
	else if (expecting == Expecting(ALPHA))
		{
		if (nextAvailableSumtNode+1 >= 2*numTaxa)
			{
			MrBayesPrint ("%s   Too many nodes on sumt tree\n", spacer);
			return (ERROR);
			}
			
		if (isTranslateDef == YES)
			{
			/* we are using the translation table */
			if (CheckString (tkn, transTo, &howMany) == ERROR)
				{
				MrBayesPrint ("%s   Could not find taxon %s in list of translation taxa\n", spacer, tkn);
				goto errorExit;
				}
			else
				{
				if (GetNameFromString (transFrom, tempName, howMany) == ERROR)
					{
					MrBayesPrint ("%s   Error getting taxon names \n", spacer);
					goto errorExit;
					}
				else
					{
					if (CheckString (tempName, taxaNames, &howMany) == ERROR)
						{
						MrBayesPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
						goto errorExit;
						}
					if (tempSet[howMany-1] == YES)
						{
						MrBayesPrint ("%s   Taxon name %s already used in tree\n", spacer, tkn);
						goto errorExit;
						}
					else
						tempSet[howMany-1] = YES;
					}
				}
			}
		else
			{
			/* Check to see if the name is in the list of taxon names. */
			if (CheckString (tkn, taxaNames, &howMany) == ERROR)
				{
				MrBayesPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
				goto errorExit;
				}
			if (tempSet[howMany-1] == YES)
				{
				MrBayesPrint ("%s   Taxon name %s already used in tree\n", spacer, tkn);
				goto errorExit;
				}
			else
				tempSet[howMany-1] = YES;
			}
		if (pSumtPtr->left == NULL)
			{
			pSumtPtr = &sumtNodes[nextAvailableSumtNode];
			strcpy (pSumtPtr->label, tkn);
			pSumtPtr->taxonName = YES;
			pSumtPtr->index = howMany - 1;
			nextAvailableSumtNode++;
			qSumtPtr->left = pSumtPtr;
			pSumtPtr->anc = qSumtPtr;
			qSumtPtr = pSumtPtr;
			}
		else if (pSumtPtr->right == NULL)
			{
			pSumtPtr = &sumtNodes[nextAvailableSumtNode];
			strcpy (pSumtPtr->label, tkn);
			pSumtPtr->taxonName = YES;
			pSumtPtr->index = howMany - 1;
			nextAvailableSumtNode++;
			qSumtPtr->right = pSumtPtr;
			pSumtPtr->anc = qSumtPtr;
			qSumtPtr = pSumtPtr;
			}
		else if (pSumtPtr->anc == NULL)
			{
			pSumtPtr = &sumtNodes[nextAvailableSumtNode];
			strcpy (pSumtPtr->label, tkn);
			pSumtPtr->taxonName = YES;
			pSumtPtr->index = howMany - 1;
			nextAvailableSumtNode++;
			qSumtPtr->anc = pSumtPtr;
			pSumtPtr->left = qSumtPtr;
			qSumtPtr = pSumtPtr;
			sumtRoot = pSumtPtr;
			pSumtPtr->marked = YES;
			isSumtTreeRooted = NO;
			}
		else
			{
			MrBayesPrint ("%s   Tree is not bifurcating\n", spacer);
			goto errorExit;
			}
		numSumtTaxa++;
		expecting  = Expecting(COMMA);
		if (sumtBrlensDef == YES)
			expecting |= Expecting(COLON);
		expecting |= Expecting(RIGHTPAR);
		}
	else if (expecting == Expecting(RIGHTPAR))
		{
		if (pSumtPtr->marked == NO)
			{
			if (pSumtPtr->anc != NULL)
				{
				pSumtPtr = pSumtPtr->anc;
				qSumtPtr = pSumtPtr;
				}
			else
				{
				MrBayesPrint ("%s   Cannot go down\n", spacer);
				goto errorExit;
				}
			}
		else
			{
			if (pSumtPtr->left != NULL)
				{
				pSumtPtr = pSumtPtr->left;
				qSumtPtr = pSumtPtr;
				}
			else
				{
				MrBayesPrint ("%s   Cannot go down\n", spacer);
				goto errorExit;
				}
			}
		expecting  = Expecting(COMMA);
		if (sumtBrlensDef == YES)
			expecting |= Expecting(COLON);
		expecting |= Expecting(RIGHTPAR);
		expecting |= Expecting(SEMICOLON);
		}
	else if (expecting == Expecting(COLON))
		{
		foundSumtColon = YES;
		expecting  = Expecting(NUMBER);
		}
	else if (expecting == Expecting(COMMA))
		{
		if (pSumtPtr->marked == NO)
			{
			if (pSumtPtr->anc != NULL)
				{
				pSumtPtr = pSumtPtr->anc;
				qSumtPtr = pSumtPtr;
				}
			else
				{
				MrBayesPrint ("%s   Cannot go down\n", spacer);
				goto errorExit;
				}
			}
		else
			{
			if (pSumtPtr->left != NULL)
				{
				pSumtPtr = pSumtPtr->left;
				qSumtPtr = pSumtPtr;
				}
			else
				{
				MrBayesPrint ("%s   Cannot go down\n", spacer);
				goto errorExit;
				}
			}
		expecting  = Expecting(ALPHA);
		expecting |= Expecting(NUMBER);
		expecting |= Expecting(LEFTPAR);
		}
	else if (expecting == Expecting(NUMBER))
		{
		if (foundSumtColon == YES)
			{
			/* branch length */
			sscanf (tkn, "%lf", &tempD);
			if (pSumtPtr->marked == NO)
				pSumtPtr->length = tempD;
			else
				{
				if (pSumtPtr->left != NULL)
					pSumtPtr->left->length = tempD;
				else
					{
					MrBayesPrint ("%s   Cannot assign branch length to left node\n", spacer);
					goto errorExit;
					}
				}
			foundSumtColon = NO;
			expecting  = Expecting(COMMA);
			expecting |= Expecting(RIGHTPAR);
			}
		else
			{
			if (isTranslateDef == YES)
				{
				/* we are using the translation table */
				if (CheckString (tkn, transTo, &howMany) == ERROR)
					{
					MrBayesPrint ("%s   Could not find taxon %s in list of translation taxa\n", spacer, tkn);
					goto errorExit;
					}
				else
					{
					if (GetNameFromString (transFrom, tempName, howMany) == ERROR)
						{
						MrBayesPrint ("%s   Error getting partition names \n", spacer);
						goto errorExit;
						}
					else
						{
						if (CheckString (tempName, taxaNames, &howMany) == ERROR)
							{
							MrBayesPrint ("%s   Could not find taxon %s in list of taxa\n", spacer, tkn);
							goto errorExit;
							}
						if (tempSet[howMany-1] == YES)
							{
							MrBayesPrint ("%s   Taxon name %s already used in tree\n", spacer, tkn);
							goto errorExit;
							}
						else
							tempSet[howMany-1] = YES;
						}
					}
				tempInt = howMany - 1;
				}
			else
				{
				/* simply use taxon number */
				sscanf (tkn, "%d", &tempInt);

				if (nextAvailableSumtNode+1 >= 2 * numTaxa)
					{
					MrBayesPrint ("%s   Too many nodes on sumt tree\n", spacer);
					goto errorExit;
					}
				/* Check to see if the name is in the list of taxon names. */
				if (CheckString (tkn, taxaNames, &howMany) == ERROR)
					{
					/* The number could not be found as a taxon name in the list of taxon names. We will
					   assume that the user has then input taxa as numbers and not the names. */
					if (tempSet[tempInt-1] == YES)
						{
						MrBayesPrint ("%s   Taxon name %d has already been used in tree\n", spacer, tempInt);
						goto errorExit;
						}
					else
						tempSet[tempInt-1] = YES;
					tempInt--;
					}
				else
					{
					/* The taxon name is in the list of taxon names */
					howMany--;
					if (howMany < 0 || howMany >= numTaxa)
						{
						MrBayesPrint ("%s   Taxon number is out of range\n", spacer);
						goto errorExit;
						}
					if (tempSet[howMany] == YES)
						{
						MrBayesPrint ("%s   Taxon %d has already been used in tree\n", spacer, howMany+1);
						goto errorExit;
						}
					else
						tempSet[howMany] = YES;
					tempInt = howMany;
					}
							
				if (GetNameFromString (taxaNames, tempName, tempInt+1) == ERROR)
					{
					MrBayesPrint ("%s   Error getting partition names \n", spacer);
					goto errorExit;
					}
				}
						
			if (pSumtPtr->left == NULL)
				{
				pSumtPtr = &sumtNodes[nextAvailableSumtNode];
				strcpy (pSumtPtr->label, tempName);
				pSumtPtr->taxonName = YES;
				pSumtPtr->index = tempInt;
				nextAvailableSumtNode++;
				qSumtPtr->left = pSumtPtr;
				pSumtPtr->anc = qSumtPtr;
				qSumtPtr = pSumtPtr;
				}
			else if (pSumtPtr->right == NULL)
				{
				pSumtPtr = &sumtNodes[nextAvailableSumtNode];
				strcpy (pSumtPtr->label, tempName);
				pSumtPtr->taxonName = YES;
				pSumtPtr->index = tempInt;
				nextAvailableSumtNode++;
				qSumtPtr->right = pSumtPtr;
				pSumtPtr->anc = qSumtPtr;
				qSumtPtr = pSumtPtr;
				}
			else if (pSumtPtr->anc == NULL)
				{
				pSumtPtr = &sumtNodes[nextAvailableSumtNode];
				strcpy (pSumtPtr->label, tempName);
				pSumtPtr->taxonName = YES;
				pSumtPtr->index = tempInt;
				nextAvailableSumtNode++;
				qSumtPtr->anc = pSumtPtr;
				pSumtPtr->left = qSumtPtr;
				qSumtPtr = pSumtPtr;
				sumtRoot = pSumtPtr;
				pSumtPtr->marked = YES;
				isSumtTreeRooted = NO;
				}
			else
				{
				MrBayesPrint ("%s   Tree is not bifurcating\n", spacer);
				return (ERROR);
				}
			numSumtTaxa++;
			expecting  = Expecting(COMMA);
			expecting |= Expecting(COLON);
			expecting |= Expecting(RIGHTPAR);
			}
		}

	return (NO_ERROR);
	
	errorExit:
		inTreesBlock = NO;
		return (ERROR);
	MrBayesPrint ("%s", parmName); /* just because I am tired of seeing the unused parameter error msg */

}





int FindParts (int nodeDepthConTree)


{

	int				i, j, nNodes, partNum;
	MrBFlt			bl, rate, x, y;
	SumtNode		**downPass, *p;

	/* allocate memory for downpass */
	downPass = (SumtNode **)SafeMalloc((size_t) (2 * numTaxa * sizeof(SumtNode *)));
	if (!downPass)
		{
		MrBayesPrint ("%s   Could not allocate downPass\n", spacer);
		goto errorExit;
		}
				
	/* get the downpass sequence */
	i = 0;
	GetSumtDownPass (sumtRoot, downPass, &i);
	nNodes = i;
		
	/* calculate the node ages (depths) if appropriate */
	if (nodeDepthConTree == YES)
		{
		if (treeIndex >= numClockRates)
			{
			MrBayesPrint ("%s   More trees than clock rate samples\n", spacer);
			return (ERROR);
			}
		rate = clockRate[treeIndex++];
		for (i=0; i<nNodes-1; i++)
			{
			p = downPass[i];
			if (p->left == NULL && p->right == NULL)
				p->age = 0.0;
			else
				{
				x = p->left->age + p->left->length / rate;
				y = p->right->age + p->right->length / rate;
				if (x > y)
					p->age = x;
				else
					p->age = y;
				}
			}
		for (i=nNodes-3; i>=0; i--)
			{
			p = downPass[i];
			p->age = p->anc->age - p->length / rate;
			}
		}

	j = 0;
	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		if (p->anc != NULL)
			{
			if (sumtBrlensDef == YES)
				{
				if (nodeDepthConTree == YES)
                    {
					bl = p->age;
                    /* correct age of terminals that are likely to be extant */
                    if (p->left == NULL && bl < (downPass[nNodes-2]->age/(100.0*nNodes)))
                        bl = 0.0;
                    }
				else
					bl = p->length;
				}
			else
				bl = -1.0;
			if (p->anc->anc != NULL || (p->anc->anc == NULL && nodeDepthConTree == YES))
				{
				if (PartFinder (&treeBits[(p->index) * taxonLongsNeeded], bl, &partNum) == ERROR)
					goto errorExit;
				treePartNums[j] = partNum;
				treePartLengths[j] = bl;
				j++;
				}
			else if (p->anc->anc == NULL && isSumtTreeRooted == NO)
				{
				if (PartFinder (&treeBits[(p->index) * taxonLongsNeeded], bl, &partNum) == ERROR)
					goto errorExit;
				treePartNums[j] = partNum;
				treePartLengths[j] = bl;
				j++;
				}
			}
		}	
		
	/* count and sort all of the tree partitions IDs */		
	numTreeParts = j;
	SortIndParts (&treePartNums[0], &treePartLengths[0], numTreeParts, YES);

	free (downPass);
	
#	if 0
	MrBayesPrint ("Partition IDs (%d): \n", numTreeParts);
	for (i=0; i<numTreeParts; i++)
		MrBayesPrint ("      %3d %lf\n", treePartNums[i], treePartLengths[i]);
#	endif

	return (NO_ERROR);
	
	errorExit:
		if (downPass)
			free (downPass);
		return (ERROR);
	
}





int FindTree (void)

{

	int			i, n, foundTree, nDiff, whichTree, *x;
	
	if (numTreeParts == 0)
		{
		MrBayesPrint ("%s   Too few tree partitions\n", spacer);
		return (ERROR);
		}
		
	if (numFullTreesFound == 0)
		{
		x = &fullTreePartIds[0];
		for (i=0; i<numTreeParts; i++)
			x[i] = treePartNums[i];
		numOfThisFullTree[0] = 1;
		numFullTreesFound++;
		whichTree = 0;
		}
	else
		{
		foundTree = NO;
		x = &fullTreePartIds[0];
		for (n=0; n<numFullTreesFound; n++)
			{
			nDiff = 0;
			for (i=0; i<numTreeParts; i++)
				{
				if (x[i] != treePartNums[i])
					{
					nDiff++;
					break;
					}
				}
			if (nDiff == 0)
				{
				foundTree = YES;
				break;
				}
			x += (2 * numTaxa);
			}
		
		if (foundTree == YES)
			{
			numOfThisFullTree[n]++;
			whichTree = n;
			}
		else
			{
			if (numFullTreesFound+1 > numFullTreesAllocated)
				{
				numFullTreesAllocated += 500;
				if (ReallocateFullTrees () == ERROR)
					return (ERROR);
				}
			
			x = &fullTreePartIds[numFullTreesFound * 2 * numTaxa];
			for (i=0; i<numTreeParts; i++)
				x[i] = treePartNums[i];
			numOfThisFullTree[numFullTreesFound] = 1;
			whichTree = numFullTreesFound;
			numFullTreesFound++;
			}
		}

	return (NO_ERROR);
	
}





void FinishSumtTree (SumtNode *p, int *i, int isThisTreeRooted)

{

	/* We only reindex the internal nodes of the tree. We
	   assume that the tip nodes have already been indexed
	   0, 1, 2, ..., numTaxa-1. */
	   
	if (p != NULL)
		{
		FinishSumtTree (p->left,  i, isThisTreeRooted);
		FinishSumtTree (p->right, i, isThisTreeRooted);
		p->marked = NO;
		if (p->left == NULL && p->right == NULL && p->anc != NULL)
			{
			}
		else if (p->left != NULL && p->right == NULL && p->anc == NULL)
			{
			if (isThisTreeRooted == YES)
				p->index = (*i)++;
			}
		else
			{
			p->index = (*i)++;
			}
		}
		
}





int FirstTaxonInPartition (safeLong *partition, int length)

{

	int				i, j, nBits, taxon;
	safeLong			x;

	nBits = sizeof(safeLong) * 8;

	for (i=taxon=0; i<length; i++)
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






int FreeBits (void)

{
	int		i;

	if (memAllocs[ALLOC_TREEBITS] == YES)
		{
		free (treeBits);
		free (treePartNums);
		free (treePartLengths);
		memAllocs[ALLOC_TREEBITS] = NO;
		}
	if (memAllocs[ALLOC_TREEPARTS] == YES)
		{
		free (treePartsFound);
		memAllocs[ALLOC_TREEPARTS] = NO;
		}
	if (memAllocs[ALLOC_NUMOFPART] == YES)
		{
		free (numFoundOfThisPart);
		if (comparingFiles == YES)
			{
			free (numFoundOfThisPart1);
			free (numFoundOfThisPart2);
			}
		memAllocs[ALLOC_NUMOFPART] = NO;
		}
	if (memAllocs[ALLOC_NUMINRUNOFPART] == YES)
		{
		for (i=0; i<sumtParams.numRuns; i++)
			{
			if (numFoundInRunOfPart[i])
				free (numFoundInRunOfPart[i]);
			}
		free (numFoundInRunOfPart);
		memAllocs[ALLOC_NUMINRUNOFPART] = NO;
		}
	if (memAllocs[ALLOC_ABRLENS] == YES)
		{
		free (aBrlens);
		memAllocs[ALLOC_ABRLENS] = NO;
		}
	if (memAllocs[ALLOC_SBRLENS] == YES)
		{
		free (sBrlens);
		memAllocs[ALLOC_SBRLENS] = NO;
		}
	if (memAllocs[ALLOC_A_WITHIN_BRLENS] == YES)
		{
		free (aWithinBrlens);
		memAllocs[ALLOC_A_WITHIN_BRLENS] = NO;
		}
	if (memAllocs[ALLOC_S_WITHIN_BRLENS] == YES)
		{
		free (sWithinBrlens);
		memAllocs[ALLOC_S_WITHIN_BRLENS] = NO;
		}
	if (memAllocs[ALLOC_SUMB] == YES)
		{
		free (sumB);
		memAllocs[ALLOC_SUMB] = NO;
		}
	if (memAllocs[ALLOC_SUMSQB] == YES)
		{
		free (sumsqB);
		memAllocs[ALLOC_SUMSQB] = NO;
		}
	if (memAllocs[ALLOC_TAXONMASK] == YES)
		{
		free (taxonMask);
		memAllocs[ALLOC_TAXONMASK] = NO;
		}
	if (memAllocs[ALLOC_TAXAFOUND] == YES)
		{
		free (sumTaxaFound);
		memAllocs[ALLOC_TAXAFOUND] = NO;
		}
	if (memAllocs[ALLOC_SUMTTREE] == YES)
		{
		free (sumtNodes);
		memAllocs[ALLOC_SUMTTREE] = NO;
		}
	if (memAllocs[ALLOC_FULLTREEINFO] == YES)
		{
		free (fullTreePartIds);
		free (numOfThisFullTree);
		memAllocs[ALLOC_FULLTREEINFO] = NO;
		}
	if (memAllocs[ALLOC_PRUNEINFO] == YES)
		{
		free (prunedTaxa);
		absentTaxa = NULL;
		memAllocs[ALLOC_PRUNEINFO] = NO;
		}
	if (memAllocs[ALLOC_FULLCOMPTREEINFO] == YES)
		{
		free (fullCompTreePartIds1);
		free (fullCompTreePartIds2);
		free (fullCompTreePartLengths1);
		free (fullCompTreePartLengths2);
		memAllocs[ALLOC_FULLCOMPTREEINFO] = NO;
		}
	if (memAllocs[ALLOC_CLOCKRATE] == YES)
		{
		free (clockRate);
		numClockRates = 0;
		memAllocs[ALLOC_CLOCKRATE] = NO;
		}
	return (NO_ERROR);

}





/* GetClockRates: Extract clock rates from .p file name 'fileName' */
int GetClockRates (char *rateName, char *fileName)

{

	int		i, j, index, lineLength, inSumpComment=NO, lineNum, lastNonDigitLine, numParamLines,
			nHeaders, allDigitLine, lastTokenWasDash, nNumbersOnThisLine, tokenType, nLines,
			numColumns, numLinesToRead, numLinesRead, burnin;
	FILE		*fp;
	char		*inputLine, *headerLine, temp[100];
	MrBFlt		tempD;

	if ((fp = OpenBinaryFileR (fileName)) == NULL)
		{
		MrBayesPrint ("%s   Could not open file '%s' to read clock rates\n", spacer, fileName);
		return (ERROR);
		}
	
	lineLength = LongestLine (fp) + 10;
	
	/* close binary file */
	SafeFclose (&fp);
	
	inputLine = (char *) calloc (3 * (lineLength), sizeof (char));
	if (inputLine == NULL)
		{
		MrBayesPrint ("%s   Could not allocate inputLine in GetClockRates\n", spacer);
		return (ERROR);
		}
	headerLine = inputLine + lineLength;
	headerNames = headerLine + lineLength;
	headerLine[0] = '\0';
	for (i=0; i<lineLength-1; i++)
		headerNames[i] = ' ';
	headerNames[i] = '\0';

	if ((fp = OpenTextFileR(fileName)) == NULL)
		{
		MrBayesPrint ("%s   Could not open text file in GetClockRates\n", spacer);
		free (inputLine);
		return (ERROR);
		}

	/* find last block */
	inComment = NO;
	lineNum = lastNonDigitLine = numParamLines = 0;
	while (fgets (inputLine, lineLength, fp) != NULL)
		{
		sumpTokenP = inputLine;
		allDigitLine = YES;
		lastTokenWasDash = NO;
		nNumbersOnThisLine = 0;
		do
			{
			GetSumpToken (&tokenType, &sumpTokenP, sumpToken);
			/*printf ("%s (%d)\n", sumpToken, tokenType);*/
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
			strcpy (headerLine, inputLine);
			}
		else
			{
			if (nNumbersOnThisLine > 0)
				numParamLines++;
			}
		}
		
	/* Now, check some aspects of the file. */
	if (inSumpComment == YES)
		{
		MrBayesPrint ("%s   Unterminated comment in file '%s'\n", spacer, fileName);
		free (inputLine);
		SafeFclose (&fp);
		return (ERROR);
		}
	if (numParamLines <= 0)
		{
		MrBayesPrint ("%s   No parameters were found in file '%s'\n", spacer, fileName);
		free (inputLine);
		SafeFclose (&fp);
		return (ERROR);
		}

    if (sumtParams.relativeBurnin == NO)
        burnin = sumtParams.sumtBurnIn;
    else
        burnin = (int) (numParamLines * sumtParams.sumtBurnInFraction);

    if (burnin > numParamLines)
		{
		MrBayesPrint ("%s   No clock rates retrieved as there are too few samples in last block of file '%s'\n", spacer, fileName);
		MrBayesPrint ("%s   Try setting burnin to a number less than %d\n", spacer, numParamLines);
		SafeFclose (&fp);
		free (inputLine);
		return (ERROR);
		}

	/* Fast forward to last block */
	(void)fseek(fp, 0L, 0);	
	for (lineNum=0; lineNum<lastNonDigitLine+burnin; lineNum++)
		fgets (inputLine, lineLength, fp);

	/* Count number of rows and check number of columns */
	inSumpComment = NO;
	nLines = 0;
	numColumns = 0;
	while (fgets (inputLine, lineLength, fp) != NULL)
		{
		sumpTokenP = inputLine;
		allDigitLine = YES;
		lastTokenWasDash = NO;
		nNumbersOnThisLine = 0;
		do
			{
			GetSumpToken (&tokenType, &sumpTokenP, sumpToken);
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
			MrBayesPrint ("%s   Found a line with non-digit characters (line %d)\n", spacer, lineNum);
			free (inputLine);
			SafeFclose (&fp);
			return (ERROR);
			}
		else
			{
			if (nNumbersOnThisLine > 0)
				{
				nLines++;
				if (nLines == 1)
					numColumns = nNumbersOnThisLine;
				else
					{
					if (nNumbersOnThisLine != numColumns)
						{
						MrBayesPrint ("%s   Number of columns is not even (%d in first line and %d in %d line of file '%s')\n", spacer, numColumns, nNumbersOnThisLine, lineNum, fileName);
						free (inputLine);
						SafeFclose (&fp);
						return (ERROR);
						}
					}
				}
			}
		}
	numRows = numClockRates = nLines;
	if (numClockRates == 0 || numColumns == 0)
		{
		MrBayesPrint ("%s   No postburnin rows or no columns in file '%s'\n", spacer, fileName);
		free (inputLine);
		SafeFclose (&fp);
		return (ERROR);
		}

	/* allocate space to hold clock rates */
	if (memAllocs[ALLOC_CLOCKRATE] == YES)
		{
		free (clockRate);
		memAllocs[ALLOC_CLOCKRATE] = NO;
		}
	clockRate = (MrBFlt *) calloc (numClockRates, sizeof (MrBFlt));
	if (!clockRate)
		{
		MrBayesPrint ("%s   Could not allocate clockRate '%s'\n", spacer);
		free (inputLine);
		SafeFclose (&fp);
		return (ERROR);
		}
	memAllocs[ALLOC_CLOCKRATE] = YES;
	for (i=0; i<numClockRates; i++)
		clockRate[i] = 1.0;

	/* separate header line into titles for each column */
	if (GetHeaders (headerLine, &nHeaders) == ERROR)
		{
		free (inputLine);
		SafeFclose (&fp);
		return (ERROR);
		}

	/* find index of column with clock rates for the tree */
	index = -1;
	for (i=0; i<nHeaders; i++)
		{
		if (GetNameFromString (headerNames, temp, i+1) == ERROR)
			{
			free (inputLine);
			SafeFclose (&fp);
			return (ERROR);
			}
		if (strcmp (temp,rateName) == 0)
			index = i;
		}

	if (index != -1)
		{
		/* Only read clock rates if they have been recorded, otherwise assume 1.0 */
		/* First, rewind file pointer to beginning of file... */
		(void)fseek(fp, 0L, 0);	
		
		/* ...and fast forward to beginning of last unburned parameter line. */
		for (lineNum=0; lineNum<lastNonDigitLine+burnin; lineNum++)
			if(fgets (inputLine, lineLength, fp)==0) 

		/* ...and parse file, line-by-line. We are only parsing lines that have digits that should be read. */
		inSumpComment = NO;
		numLinesToRead = numParamLines - burnin;
		numLinesRead = j = 0;
		while (fgets (inputLine, lineLength, fp) != NULL)
			{
			sumpTokenP = inputLine;
			allDigitLine = YES;
			lastTokenWasDash = NO;
			nNumbersOnThisLine = 0;
			do
				{
				GetSumpToken (&tokenType, &sumpTokenP, sumpToken);
				if (IsSame("[", sumpToken) == SAME)
					inSumpComment = YES;
				if (IsSame("]", sumpToken) == SAME)
					inSumpComment = NO;
				if (inSumpComment == NO)
					{
					if (tokenType == NUMBER)
						{
						/* read the information from this line */
						if (j >= numRows * numColumns)
							{
							MrBayesPrint ("%s   Too many parameter values read in (%d)\n", spacer, j);
							free (inputLine);
							SafeFclose (&fp);
							return (ERROR);
							}
						sscanf (sumpToken, "%lf", &tempD);
						if (lastTokenWasDash == YES)
							tempD *= -1.0;
						if (nNumbersOnThisLine == index)
							clockRate[numLinesRead] = tempD;
						j++;
						nNumbersOnThisLine++;
						lastTokenWasDash = NO;
						}
					else if (tokenType == DASH)
						{
						lastTokenWasDash = YES;
						}
					else if (tokenType != UNKNOWN_TOKEN_TYPE)
						{
						/* we have a problem */
						MrBayesPrint ("%s   Found a line with non-digit characters (line %d) in file '%s'\n", spacer, lineNum, fileNum);
						free (inputLine);
						SafeFclose (&fp);
						return (ERROR);
						}
					}
				} while (*sumpToken);
			lineNum++;
			if (nNumbersOnThisLine > 0)
				numLinesRead++;
			}
		/* Check how many parameter line was read in. */
		if (numLinesRead != numParamLines - burnin)
			{
			MrBayesPrint ("%s   Unable to read all lines that should contain clock rates\n", spacer);
			SafeFclose (&fp);
			free (inputLine);
			return (ERROR);
			}
		}

	SafeFclose (&fp);
	free (inputLine);

	return (NO_ERROR);
}





void GetConDownPass (PolyNode **downPass, PolyNode *p, int *i)

{

	PolyNode	*q;
	
	if (p->left != NULL)
		{
		for (q=p->left; q!=NULL; q=q->sib)
			GetConDownPass(downPass, q, i);
		}

	downPass[(*i)++] = p;

}





/* get the actual down pass sequences */
void GetSumtDownPass (SumtNode *p, SumtNode **dp, int *i)

{
	
	if (p != NULL )
		{
		GetSumtDownPass (p->left,  dp, i);
		GetSumtDownPass (p->right, dp, i);
		if (p->left != NULL && p->right != NULL && p->anc != NULL)
			{
			dp[(*i)++] = p;
			}
		else if (p->left == NULL && p->right == NULL && p->anc != NULL)
			{
			dp[(*i)++] = p;
			}
		else if (p->left != NULL && p->right == NULL && p->anc == NULL)
			{
			dp[(*i)++] = p;
			}
		}
		
}





int GetPartitions (void)

{

	int				i, nNodes;
	SumtNode		**downPass, *p;

	/* allocate memory for downpass */
	downPass = (SumtNode **)SafeMalloc((size_t) (2 * numTaxa * sizeof(SumtNode *)));
	if (!downPass)
		{
		MrBayesPrint ("%s   Could not allocate downPass\n", spacer);
		goto errorExit;
		}
		
	/* set all partitions to 0 */
	for (i=0; i<2*numTaxa*taxonLongsNeeded; i++)
		treeBits[i] = 0;
		
	/* get the downpass sequence */
	i = 0;
	GetSumtDownPass (sumtRoot, downPass, &i);
	nNodes = i;
		
	/* find which taxa are in the tree */
	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		if (p->left == NULL && p->right == NULL)
			{
			AssignTipPart (p->index, &treeBits[(p->index) * taxonLongsNeeded]);
			}
		else if (p->anc != NULL)
			{
			AssignIntPart (&treeBits[(p->left->index) * taxonLongsNeeded], &treeBits[(p->right->index) * taxonLongsNeeded], &treeBits[(p->index) * taxonLongsNeeded]);
			}
		}
				
#	if 0
	MrBayesPrint ("Partitions:\n");
	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		if (p->left == NULL && p->right == NULL)
			{
			ShowParts (stdout, &treeBits[(p->index) * taxonLongsNeeded], numTaxa);
			MrBayesPrint ("\n");
			}
		else if (p->anc != NULL)
			{
			if (p->anc->anc != NULL)
				{
				ShowParts (stdout, &treeBits[(p->index) * taxonLongsNeeded], numTaxa);
				MrBayesPrint ("\n");
				}
			else if (p->anc->anc == NULL && isSumtTreeRooted == NO)
				{
				ShowParts (stdout, &treeBits[(p->index) * taxonLongsNeeded], numTaxa);
				MrBayesPrint ("\n");
				}
			}
		}
#	endif

	free (downPass);

	return (NO_ERROR);
	
	errorExit:
		if (downPass)
			free (downPass);
		return (ERROR);
	
}





void GetSumtToken (int *tokenType)

{
		
	int				allNumbers, foundExp, foundExpSign;
	register char	*temp;
	
	(*tokenType) = 0;
	temp = sumtToken;
	
	while (IsWhite(*sumtTokenP) == 1 || IsWhite(*sumtTokenP) == 2)
		{
		if (IsWhite(*sumtTokenP) == 2)
			{
			*tokenType = RETURNSYMBOL;
			foundNewLine = YES;
			/* MrBayesPrint ("RETURN\n"); */
			}
		++sumtTokenP;
		}
	
	*tokenType = UNKNOWN_TOKEN_TYPE;
	if (IsIn(*sumtTokenP,"="))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = EQUALSIGN;
		}
	else if (IsIn(*sumtTokenP,";"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = SEMICOLON;
		}
	else if (IsIn(*sumtTokenP,":"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = COLON;
		}
	else if (IsIn(*sumtTokenP,","))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = COMMA;
		}
	else if (IsIn(*sumtTokenP,"#"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = POUNDSIGN;
		}
	else if (IsIn(*sumtTokenP,"("))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = LEFTPAR;
		}
	else if (IsIn(*sumtTokenP,")"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = RIGHTPAR;
		}
	else if (IsIn(*sumtTokenP,"{"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = LEFTCURL;
		}
	else if (IsIn(*sumtTokenP,"}"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = RIGHTCURL;
		}
	else if (IsIn(*sumtTokenP,"["))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = LEFTCOMMENT;
		}
	else if (IsIn(*sumtTokenP,"]"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = RIGHTCOMMENT;
		}
	else if (IsIn(*sumtTokenP,"?"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = QUESTIONMARK;
		}
	else if (IsIn(*sumtTokenP,"-"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = DASH;
		}
	else if (IsIn(*sumtTokenP,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789."))
		{
		allNumbers = TRUE;
		if (!IsIn(*sumtTokenP,"0123456789."))
			allNumbers = FALSE;
        foundExp = foundExpSign = FALSE;
		*temp++ = *sumtTokenP++;
		while(IsIn(*sumtTokenP,"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.-+"))
			{
			if(allNumbers == TRUE && !IsIn(sumtTokenP[-1],"Ee") && *sumtTokenP=='-')
                break;
            else if (allNumbers == TRUE && IsIn(*sumtTokenP,"Ee") && foundExp == NO)
                foundExp = TRUE;
            else if (allNumbers == TRUE && IsIn(*sumtTokenP,"+-") && IsIn(sumtTokenP[-1],"Ee"))
                foundExpSign = TRUE;
            else if (!IsIn(*sumtTokenP,"0123456789."))
				allNumbers = FALSE;
            *temp++ = *sumtTokenP++;
			}
		if (allNumbers == TRUE)
			*tokenType = NUMBER;
		else
			*tokenType = ALPHA;
		}
	else if (IsIn(*sumtTokenP,"*"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = ASTERISK;
		}
	else if (IsIn(*sumtTokenP,"/"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = FORWARDSLASH;
		}
	else if (IsIn(*sumtTokenP,"'\\'"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = BACKSLASH;
		}
	else if (IsIn(*sumtTokenP,"!"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = EXCLAMATIONMARK;
		}
	else if (IsIn(*sumtTokenP,"%"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = PERCENT;
		}
	else if (IsIn(*sumtTokenP,"\""))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = QUOTATIONMARK;
		}
	else if (IsIn(*sumtTokenP,"&~+^$@|{}`><"))
		{
		*temp++ = *sumtTokenP++;
		*tokenType = WEIRD;
		}

	*temp = '\0';
	
}





/* Label: Calculate length of label and fill in char *label if not NULL */
int Label (PolyNode *p, int addIndex, char *label, int maxLength)
{
    int     i, j0, j1, k, n, length, nameLength;

    if (p == NULL)
        return 0;

    /* first calculate length */
    if (addIndex != NO)
        length = (int)(strlen(p->label)) + 4 + (int)(log10(p->index+1));
    else
        length = (int)(strlen(p->label));
    length = (length > maxLength ? maxLength : length);

    /* fill in label if label != NULL */
    if (label != NULL)
        {
        if (addIndex != NO)
            nameLength = length - 4 - (int)(log10(p->index+1));
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
            n = p->index + 1;
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





int OpenBrlensFile (int treeNo)

{

	int			i, len;
	char		fileName[100];

	/* set file name */
	if (sumtParams.numTrees > 1)
		sprintf (fileName, "%s.tree%d.brlens", sumtParams.sumtFileName, treeNo+1);
	else
		sprintf (fileName, "%s.brlens", sumtParams.sumtFileName);

	/* open file checking for over-write as appropriate */
	if ((fpBrlens = OpenNewMBPrintFile(fileName)) == NULL)
		return ERROR;

	/* print unique identifier to the file */
	len = (int) strlen (stamp);
	if (len <= 1)
		{
		MrBayesPrintf (fpBrlens, "[ID: None Available]\n");
		}
	else
		{
		MrBayesPrintf (fpBrlens, "[ID: %s]\n", stamp);
		}

	/* print header */
	for (i=0; i<numBrlens; i++)
		{
		if (i == numBrlens - 1)	
			MrBayesPrintf (fpBrlens, "v%d\n", i+1);
		else
			MrBayesPrintf (fpBrlens, "v%d\t", i+1);
		}

	return (NO_ERROR);
}





int OpenComptFiles (void)

{

	int			len;
	char		pFilename[100], dFilename[100];

	/* set file names */
	strcpy (pFilename, comptreeParams.comptOutfile);
	strcpy (dFilename, comptreeParams.comptOutfile);
	strcat (pFilename, ".parts");
	strcat (dFilename, ".dists");

	/* Check to see if files are present. We don't want to
	   inadvertantly over-write files unless the user specifies
	   that we do so. */
	if ((fpCompParts = OpenNewMBPrintFile (pFilename)) == NULL)
		return ERROR;
	if ((fpCompDists = OpenNewMBPrintFile (dFilename)) == NULL)
		return ERROR;
		
	/* print unique identifiers to each file */
	len = (int) strlen (stamp);
	if (len <= 1)
		{
		fprintf (fpCompParts, "[ID: None Available]\n");
		fprintf (fpCompDists, "[ID: None Available]\n");
		}
	else
		{
		fprintf (fpCompParts, "[ID: %s]\n", stamp);
		fprintf (fpCompDists, "[ID: %s]\n", stamp);
		}

	return (NO_ERROR);
}





int OpenSumtFiles (int treeNo)

{

	int			len;
	char		pFilename[100], cFilename[100], tFilename[100];

	/* set file names */
	if (sumtParams.numTrees > 1)
		{
        sprintf (pFilename, "%s.tree%d.parts", sumtParams.sumtOutfile, treeNo+1);
		sprintf (cFilename, "%s.tree%d.con", sumtParams.sumtOutfile, treeNo+1);
		sprintf (tFilename, "%s.tree%d.trprobs", sumtParams.sumtOutfile, treeNo+1);
		}
	else
		{
		sprintf (pFilename, "%s.parts", sumtParams.sumtFileName);
		sprintf (cFilename, "%s.con", sumtParams.sumtFileName);
		sprintf (tFilename, "%s.trprobs", sumtParams.sumtFileName);
		}
	
	/* Open files checking for over-write as appropriate */
	if ((fpParts = OpenNewMBPrintFile(pFilename)) == NULL)
		return ERROR;
	if ((fpCon = OpenNewMBPrintFile(cFilename)) == NULL)
		{
		SafeFclose (&fpParts);
		return ERROR;
		}
	if (sumtParams.calcTrprobs == YES)
		{
		if ((fpTrees = OpenNewMBPrintFile(tFilename)) == NULL)
			{
			SafeFclose (&fpParts);
			SafeFclose (&fpCon);
			return ERROR;
			}
		}

	/* print #NEXUS if appropriate */
	fprintf (fpCon,   "#NEXUS\n\n");
	if (sumtParams.calcTrprobs == YES)
		fprintf (fpTrees, "#NEXUS\n\n");

	/* print unique identifiers to each file */
	len = (int) strlen (stamp);
	if (len <= 1)
		{
		if (sumtParams.calcTrprobs == YES)
			fprintf (fpParts, "[ID: None Available]\n");
		fprintf (fpCon,   "[ID: None Available]\n");
		fprintf (fpTrees, "[ID: None Available]\n");
		}
	else
		{
		if (sumtParams.calcTrprobs == YES)
			fprintf (fpParts, "[ID: %s]\n", stamp);
		fprintf (fpCon,   "[ID: %s]\n", stamp);
		fprintf (fpTrees, "[ID: %s]\n", stamp);
		}

	return (NO_ERROR);		
}





int PartFinder (safeLong *p, MrBFlt bl, int *partID)

{
			
	int			i, n, foundPart, nDiff, whichPart;
	safeLong		*x;
	MrBFlt		f;
	
	if (printingBrlens == YES)
		{
		foundPart = NO;
		x = &treePartsFound[0];
		for (n=0; n<numBrlens; n++)
			{
			for (i=0; i<taxonLongsNeeded; i++)
				if (x[i] != p[i])
					break;
			if (i == taxonLongsNeeded)
				{
				foundPart = YES;
				break;
				}
			x += taxonLongsNeeded;
			}
		
		if (foundPart == YES)
			{
			brlens[n] = bl;
			}
		return (NO_ERROR);

		}

	foundPart = NO;
	x = &treePartsFound[0];
	for (n=0; n<numTreePartsFound; n++)
		{
		nDiff = 0;
		for (i=0; i<taxonLongsNeeded; i++)
			if (x[i] != p[i])
				nDiff++;
		if (nDiff == 0)
			{
			foundPart = YES;
			break;
			}
		x += taxonLongsNeeded;
		}
	
	if (foundPart == YES)
		{
		numFoundOfThisPart[n]++;
		if (sumtParams.numRuns > 1)
			numFoundInRunOfPart[runIndex][n]++;
		if (comparingFiles == YES)
			{
			if (fileNum == 0)
				numFoundOfThisPart1[n]++;
			else
				numFoundOfThisPart2[n]++;
			}
		whichPart = n;
		}
	else
		{
		if (numTreePartsFound+1 > numPartsAllocated)
			{
			numPartsAllocated += 500;
			if (ReallocateBits () == ERROR)
				return (ERROR);
			}
		
		x = &treePartsFound[numTreePartsFound * taxonLongsNeeded];
		for (i=0; i<taxonLongsNeeded; i++)
			x[i] = p[i];
		numFoundOfThisPart[numTreePartsFound] = 1;
		if (sumtParams.numRuns > 1)
			{
			for (i=0; i<sumtParams.numRuns; i++)
				numFoundInRunOfPart[i][numTreePartsFound] = 0;
			numFoundInRunOfPart[runIndex][numTreePartsFound] = 1;
			}
		if (comparingFiles == YES)
			{
			if (fileNum == 0)
				numFoundOfThisPart1[numTreePartsFound] = 1;
			else
				numFoundOfThisPart2[numTreePartsFound] = 1;
			}
		whichPart = numTreePartsFound;
		numTreePartsFound++;
		}

	(*partID) = whichPart;

	if (bl >= 0.0)
		{
		if (numFoundOfThisPart[whichPart] == 1)
			{
			aBrlens[whichPart] = bl;
			sBrlens[whichPart] = 0.0;
			}
		else
			{
			f = aBrlens[whichPart];
			aBrlens[whichPart] += (bl - aBrlens[whichPart]) / (MrBFlt) (numFoundOfThisPart[whichPart]);
			sBrlens[whichPart] += (bl - aBrlens[whichPart]) * (bl - f);
			}
		if (sumtParams.numRuns > 1)
			{
			if (numFoundOfThisPart[whichPart] == 1)
				{
				aWithinBrlens[whichPart] = bl;
				sWithinBrlens[whichPart] = 0.0;
				sumB[whichPart] = bl;
				sumsqB[whichPart] = bl * bl;
				}
			else if (numFoundInRunOfPart[runIndex][whichPart] == 1)
				{
				aWithinBrlens[whichPart] = bl;
				sumB[whichPart] += bl;
				sumsqB[whichPart] += bl * bl;
				}
			else
				{
				f = aWithinBrlens[whichPart];
				aWithinBrlens[whichPart] += (bl - aWithinBrlens[whichPart]) / (MrBFlt) (numFoundInRunOfPart[runIndex][whichPart]);
				sWithinBrlens[whichPart] += (bl - aWithinBrlens[whichPart]) * (bl - f);
				sumB[whichPart] += (aWithinBrlens[whichPart] - f);
				sumsqB[whichPart] += ((aWithinBrlens[whichPart] * aWithinBrlens[whichPart]) - (f * f));
				}
			}
		}
	
	return (NO_ERROR);

}





int PrintBrlensToFile (void)

{
	int		i;

	/* print header */
	for (i=0; i<numBrlens; i++)
		{
		if (brlens[i] < 0.0)
			MrBayesPrintf (fpBrlens, "N/A");
		else
			MrBayesPrintf (fpBrlens, "%.6f", brlens[i]);

		if (i == numBrlens - 1)
			MrBayesPrintf (fpBrlens, "\n");
		else
			MrBayesPrintf (fpBrlens, "\t");
		}

	return NO_ERROR;
}





void PrintParts (FILE *fp, safeLong *p, int nTaxaToShow)

{

	int			i, flipBits;
	safeLong		x, y;
	
	flipBits = NO;
	if (isSumtTreeRooted == NO)
		{
		y = p[0];
		x = 1 << (0 % nBitsInALong);
		if (x & y)
			flipBits = YES;
		}

	for (i=0; i<nTaxaToShow; i++)
		{
		x = 0;
		y = p[i / nBitsInALong];
		x = 1 << (i % nBitsInALong);
		if (flipBits == NO)
			{
			if ((x & y) == 0)
				MrBayesPrintf (fp, ".");
			else
				MrBayesPrintf (fp, "*");
			}
		else
			{
			if ((x & y) == 0)
				MrBayesPrintf (fp, "*");
			else
				MrBayesPrintf (fp, ".");
			}
		}

}





int PruneSumt (void)

{

	int 			i, nNodes, wasDerooted, whichTaxon, deletedOne;
	SumtNode		**downPass, *p, *q, *sis, *qAnc;

	/* allocate memory for downpass */
	downPass = (SumtNode **)SafeMalloc((size_t) (2 * numTaxa * sizeof(SumtNode *)));
	if (!downPass)
		{
		MrBayesPrint ("%s   Could not allocate downPass\n", spacer);
		goto errorExit;
		}
		
	/* root tree, if it was previously unrooted */
	wasDerooted = NO;
	if (isSumtTreeRooted == NO)
		{
		if (RootSumtTree (sumtRoot, numSumtTaxa, outGroupNum) == ERROR)
			goto errorExit;
		wasDerooted = YES;
		}
	
	/* get the downpass sequence */
	i = 0;
	GetSumtDownPass (sumtRoot, downPass, &i);
	nNodes = i;
	
	/* find which taxa are in the tree */
	do
		{
		deletedOne = NO;
		for (i=0; i<nNodes; i++)
			{
			p = downPass[i];
			if (p->left == NULL && p->right == NULL)
				{
				if (CheckString (p->label, taxaNames, &whichTaxon) == ERROR)
					{
					MrBayesPrint ("%s   Could not find taxon %s in original list of taxa\n", spacer, p->label);
					goto errorExit;
					}
				whichTaxon--;
				if (taxaInfo[whichTaxon].isDeleted == YES)
					{
					prunedTaxa[whichTaxon] = YES;
					q = p->anc;
					if (q->left == p)
						sis = q->right;
					else
						sis = q->left;
					sis->length += q->length;
					if (q->anc == NULL)
						{
						MrBayesPrint ("%s   Could not find root of tree\n", spacer);
						goto errorExit;
						}
					else
						qAnc = q->anc;
					if (qAnc->left == q)
						{
						qAnc->left = sis;
						sis->anc = qAnc;
						}
					else
						{
						qAnc->right = sis;
						sis->anc = qAnc;
						}
					p->left = p->right = p->anc = NULL;
					q->left = q->right = q->anc = NULL;
					numSumtTaxa--;
					deletedOne = YES;
					break;
					}
				}
			}
		if (deletedOne == YES)
			{
			i = 0;
			GetSumtDownPass (sumtRoot, downPass, &i);
			nNodes = i;
			}
		} while (deletedOne == YES);
		
	/* clean up on way out of function */
	free (downPass);
	
	/* deroot tree, if we rooted it earlier */
	if (wasDerooted == YES)
		{
		if (DerootSumtTree (sumtRoot, numSumtTaxa, outGroupNum) == ERROR)
			goto errorExit;
		}

#	if 0
	if (isSumtTreeRooted == YES)
		MrBayesPrint ("%s   Rooted tree %d:\n", spacer, numSumtTrees);
	else
		MrBayesPrint ("%s   Unrooted tree %d:\n", spacer, numSumtTrees);
	ShowNodes (sumtRoot, 3, isSumtTreeRooted);
#	endif

	return (NO_ERROR);
	
	errorExit:
		if (downPass)
			free (downPass);
		return(ERROR);

}





int ReallocateBits (void)

{
	int		i;

	numPartsAllocated += 100;
	
	if (memAllocs[ALLOC_TREEPARTS] == NO)
		goto errorExit;
	treePartsFound = (safeLong *)realloc((void *)treePartsFound, (size_t) (taxonLongsNeeded * numPartsAllocated * sizeof(safeLong)));
	if (!treePartsFound)
		{
		MrBayesPrint ("%s   Problem reallocating treePartsFound (%d)\n", spacer, taxonLongsNeeded * numPartsAllocated * sizeof(safeLong));
		goto errorExit;
		}
		
	if (memAllocs[ALLOC_NUMOFPART] == NO)
		goto errorExit;
	numFoundOfThisPart = (int *)realloc((void *)numFoundOfThisPart, (size_t) (numPartsAllocated * sizeof(int)));
	if (!numFoundOfThisPart)
		{
		MrBayesPrint ("%s   Problem reallocating numFoundOfThisPart (%d)\n", spacer, numPartsAllocated * sizeof(int));
		goto errorExit;
		}

	if (sumtParams.numRuns > 1)
		{
		if (memAllocs[ALLOC_NUMINRUNOFPART] == NO)
			goto errorExit;
		for (i=0; i<sumtParams.numRuns; i++) {
			numFoundInRunOfPart[i] = (int *)realloc((void *)numFoundInRunOfPart[i], (size_t) (numPartsAllocated * sizeof(int)));
			if (!numFoundInRunOfPart[i])
				{
				MrBayesPrint ("%s   Problem reallocating numFoundInRunOfPart (%d)\n", spacer, numPartsAllocated * sizeof(int));
				goto errorExit;
				}
			}
		}
	if (comparingFiles == YES)
		{
		numFoundOfThisPart1 = (int *)realloc((void *)numFoundOfThisPart1, (size_t) (numPartsAllocated * sizeof(int)));
		if (!numFoundOfThisPart1)
			{
			MrBayesPrint ("%s   Problem reallocating numFoundOfThisPart1 (%d)\n", spacer, numPartsAllocated * sizeof(int));
			goto errorExit;
			}
		numFoundOfThisPart2 = (int *)realloc((void *)numFoundOfThisPart2, (size_t) (numPartsAllocated * sizeof(int)));
		if (!numFoundOfThisPart)
			{
			MrBayesPrint ("%s   Problem reallocating numFoundOfThisPart2 (%d)\n", spacer, numPartsAllocated * sizeof(int));
			goto errorExit;
			}
		}
		
	if (sumtBrlensDef == YES)
		{
		if (memAllocs[ALLOC_ABRLENS] == NO)
			goto errorExit;
		aBrlens = (MrBFlt *)realloc((void *)aBrlens, (size_t) (numPartsAllocated * sizeof(MrBFlt)));
		if (!aBrlens)
			{
			MrBayesPrint ("%s   Problem reallocating aBrlens (%d)\n", spacer, numPartsAllocated * sizeof(MrBFlt));
			goto errorExit;
			}
			
		if (memAllocs[ALLOC_SBRLENS] == NO)
			goto errorExit;
		sBrlens = (MrBFlt *)realloc((void *)sBrlens, (size_t) (numPartsAllocated * sizeof(MrBFlt)));
		if (!sBrlens)
			{
			MrBayesPrint ("%s   Problem reallocating sBrlens (%d)\n", spacer, numPartsAllocated * sizeof(MrBFlt));
			goto errorExit;
			}
		}
		
	if (sumtParams.numRuns > 1 && sumtBrlensDef == YES)
		{
		if (memAllocs[ALLOC_A_WITHIN_BRLENS] == NO)
			goto errorExit;
		aWithinBrlens = (MrBFlt *)realloc((void *)aWithinBrlens, (size_t) (numPartsAllocated * sizeof(MrBFlt)));
		if (!aWithinBrlens)
			{
			MrBayesPrint ("%s   Problem reallocating aWithinBrlens (%d bytes)\n", spacer, numPartsAllocated * sizeof(MrBFlt));
			goto errorExit;
			}
		
		if (memAllocs[ALLOC_S_WITHIN_BRLENS] == NO)
			goto errorExit;
		sWithinBrlens = (MrBFlt *)realloc((void *)sWithinBrlens, (size_t) (numPartsAllocated * sizeof(MrBFlt)));
		if (!sWithinBrlens)
			{
			MrBayesPrint ("%s   Problem reallocating sWithinBrlens (%d bytes)\n", spacer, numPartsAllocated * sizeof(MrBFlt));
			goto errorExit;
			}

		if (memAllocs[ALLOC_SUMB] == NO)
			goto errorExit;
		sumB = (MrBFlt *)realloc((void *)sumB, (size_t) (numPartsAllocated * sizeof(MrBFlt)));
		if (!sumB)
			{
			MrBayesPrint ("%s   Problem reallocating sumB (%d bytes)\n", spacer, numPartsAllocated * sizeof(MrBFlt));
			goto errorExit;
			}

		if (memAllocs[ALLOC_SUMSQB] == NO)
			goto errorExit;
		sumsqB = (MrBFlt *)realloc((void *)sumsqB, (size_t) (numPartsAllocated * sizeof(MrBFlt)));
		if (!sumsqB)
			{
			MrBayesPrint ("%s   Problem reallocating sumsqB (%d bytes)\n", spacer, numPartsAllocated * sizeof(MrBFlt));
			goto errorExit;
			}			
		}

	return (NO_ERROR);
	
	errorExit:
		FreeBits();
		return (ERROR);
	
}





int ReallocateFullCompTrees (int whichList)

{

	numFullCompTreesAllocated[whichList] += 100;
	
	if (memAllocs[ALLOC_FULLCOMPTREEINFO] == NO)
		goto errorExit;
		
	if (whichList == 0)
		{
		fullCompTreePartIds1 = (int *)realloc((void *)fullCompTreePartIds1, (size_t) (2 * numTaxa * numFullCompTreesAllocated[whichList] * sizeof(int)));
		if (!fullCompTreePartIds1)
			{
			MrBayesPrint ("%s   Problem reallocating fullCompTreePartIds1 (%d)\n", spacer, 2 * numTaxa * numFullCompTreesAllocated[whichList] * sizeof(int));
			goto errorExit;
			}
		fullCompTreePartLengths1 = (MrBFlt *)realloc((void *)fullCompTreePartLengths1, (size_t) (2 * numTaxa * numFullCompTreesAllocated[whichList] * sizeof(MrBFlt)));
		if (!fullCompTreePartLengths1)
			{
			MrBayesPrint ("%s   Problem reallocating fullCompTreePartLengths1 (%d)\n", spacer, 2 * numTaxa * numFullCompTreesAllocated[whichList] * sizeof(MrBFlt));
			goto errorExit;
			}
		}
	else
		{
		fullCompTreePartIds2 = (int *)realloc((void *)fullCompTreePartIds2, (size_t) (2 * numTaxa * numFullCompTreesAllocated[whichList] * sizeof(int)));
		if (!fullCompTreePartIds2)
			{
			MrBayesPrint ("%s   Problem reallocating fullCompTreePartIds2 (%d)\n", spacer, 2 * numTaxa * numFullCompTreesAllocated[whichList] * sizeof(int));
			goto errorExit;
			}
		fullCompTreePartLengths2 = (MrBFlt *)realloc((void *)fullCompTreePartLengths2, (size_t) (2 * numTaxa * numFullCompTreesAllocated[whichList] * sizeof(MrBFlt)));
		if (!fullCompTreePartLengths2)
			{
			MrBayesPrint ("%s   Problem reallocating fullCompTreePartLengths2 (%d)\n", spacer, 2 * numTaxa * numFullCompTreesAllocated[whichList] * sizeof(MrBFlt));
			goto errorExit;
			}
		}

	return (NO_ERROR);
	
	errorExit:
		FreeBits();
		return (ERROR);
	
}





int ReallocateFullTrees (void)

{


	numFullTreesAllocated += 100;
	
	if (memAllocs[ALLOC_FULLTREEINFO] == NO)
		goto errorExit;
		
	numOfThisFullTree = (int *)realloc((void *)numOfThisFullTree, (size_t) (numFullTreesAllocated * sizeof(int)));
	if (!numOfThisFullTree)
		{
		MrBayesPrint ("%s   Problem reallocating numOfThisFullTree (%d)\n", spacer, numFullTreesAllocated * sizeof(int));
		goto errorExit;
		}
		
	fullTreePartIds = (int *)realloc((void *)fullTreePartIds, (size_t) (2 * numTaxa * numFullTreesAllocated * sizeof(int)));
	if (!fullTreePartIds)
		{
		MrBayesPrint ("%s   Problem reallocating fullTreePartIds (%d)\n", spacer, 2 * numTaxa * numFullTreesAllocated * sizeof(int));
		goto errorExit;
		}

	return (NO_ERROR);
	
	errorExit:
		FreeBits();
		return (ERROR);
	
}





int ReorderParts (void)

{

	int			n, i, j;
	safeLong		*x, y, z, *newBits;

	newBits = (safeLong *)SafeMalloc((size_t) (taxonLongsNeeded * sizeof(safeLong)));
	if (!newBits)
		{
		MrBayesPrint ("%s   Could not allocate newBits\n", spacer);
		return (ERROR);
		}
		
	numIncludedTaxa = 0;
	for (i=0; i<numTaxa; i++)
		if (sumTaxaFound[i] == YES)
			numIncludedTaxa++;

	x = &treePartsFound[0];
	for (n=0; n<numTreePartsFound; n++)
		{
		for (i=0; i<taxonLongsNeeded; i++)
			newBits[i] = 0;

		j = 0;
		for (i=0; i<numTaxa; i++)
			{
			if (sumTaxaFound[i] == YES)
				{
				y = x[i / nBitsInALong];
				z = 1 << (i % nBitsInALong);
				if (y & z)
					{
					z = 1 << (j % nBitsInALong);
					newBits[j/nBitsInALong] |= z;
					}
				j++;
				}
			}
			
		for (i=0; i<taxonLongsNeeded; i++)
			x[i] = newBits[i];
		x += taxonLongsNeeded;
		}
		
	free (newBits);
		
	return (NO_ERROR);
		
}





void ResetTranslateTable (void) {

    int i;

    isTranslateDef = NO;
	numTranslates = whichTranslate = 0;
	for (i=0; i<numTaxa*100; i++)
		{
		transFrom[i] = ' ';
		transTo[i] = ' ';
		if (i == numTaxa*100 - 1)
			{
			transFrom[i] = '\0';
			transTo[i] = '\0';
			}
		}

}





int RootSumtTree (SumtNode *p, int n, int out)

{

	int 			i, j, nNodes, *usedMemIndex=NULL, availableMemIndex[2], isMarked, sumtOut=0, 
					*localTaxaFound=NULL, localOutgroupNum;
	MrBFlt			tempBrLen;
	SumtNode		**downPass=NULL, *first, *second, *lft, *rht, *m1, *m2, *um1, *um2;

	/* get down pass sequence */
	if (isSumtTreeRooted == YES)
		{
		MrBayesPrint ("%s   Tree is already rooted\n", spacer);
		goto errorExit;
		}
	else
		{
		nNodes = 2 * n - 2;
		}
		
	/* allocate memory */
	downPass = (SumtNode **)SafeMalloc((size_t) (2 * numTaxa * sizeof(SumtNode *)));
	if (!downPass)
		{
		MrBayesPrint ("%s   Could not allocate downPass\n", spacer);
		goto errorExit;
		}
	i = 0;
	GetSumtDownPass (sumtRoot, downPass, &i);
	nNodes = i;
	localTaxaFound = (int *)SafeMalloc((size_t) (numTaxa * sizeof(int)));
	if (!localTaxaFound)
		{
		MrBayesPrint ("%s   Could not allocate localTaxaFound\n", spacer);
		goto errorExit;
		}
	for (i=0; i<numTaxa; i++)
		localTaxaFound[i] = NO;
	usedMemIndex = (int *)SafeMalloc((size_t) (2 * numTaxa * sizeof(int)));
	if (!usedMemIndex)
		{
		MrBayesPrint ("%s   Could not allocate usedMemIndex\n", spacer);
		goto errorExit;
		}
	for (i=0; i<2*numTaxa; i++)
		usedMemIndex[i] = NO;

	/* set local outgroup number */
	localOutgroupNum = outGroupNum;

	/* find available nodes */
	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		usedMemIndex[p->memoryIndex] = YES;
		}	
	j = 0;
	for (i=0; i<2*numTaxa; i++)
		{
		if (usedMemIndex[i] == NO)
			{
			if (j <= 1)
				{
				availableMemIndex[j] = i;
				j++;
				}
			else
				{
				/* some taxa are not included in the trees */
				}
			}
		}
	
	first  = &sumtNodes[availableMemIndex[0]];
	second = &sumtNodes[availableMemIndex[1]];
	
	/* root tree with previously down taxon as first right */
	lft = sumtRoot->left;
	rht = sumtRoot;
	tempBrLen = lft->length;
	if (tempBrLen <= 0.0)
		tempBrLen = 0.0;
	lft->anc = rht->anc = first;
	rht->left = rht->right = NULL;
	first->left = lft;
	first->right = rht;
	first->anc = second;
	second->left = first;
	second->right = second->anc = NULL;
	lft->length = rht->length = tempBrLen * (MrBFlt) 0.5;
	first->length = second->length = 0.0;
	sumtRoot = second;
	isSumtTreeRooted = YES;
	
	/* update downpass sequence */
	i = 0;
	GetSumtDownPass (sumtRoot, downPass, &i);
	nNodes = i;
	
	/* now, bring the outgroup around to the first right position */
	isMarked = NO;
	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		if (p->left == NULL && p->right == NULL)
			{
			localTaxaFound[p->index] = YES;
			if (p->index == out)
				{
				p->marked = YES;
				isMarked = YES;
				}
			else	
				p->marked = NO;
			}
		}	

	if (isMarked == NO)
		{
		/* We will arbitrarily root it by the first taxon in the matrix
		   that was actually found on the tree */
		for (i=0; i<numTaxa; i++)
			{
			if (localTaxaFound[i] == YES)
				{
				sumtOut = i;
				break;
				}
			}
		isMarked = NO;
		for (i=0; i<2*n; i++)
			{
			p = downPass[i];
			if (p->left == NULL && p->right == NULL)
				{
				localTaxaFound[p->index] = YES;
				if (p->index == sumtOut)
					{
					p->marked = YES;
					isMarked = YES;
					}
				else	
					p->marked = NO;
				}
			}	
		if (isMarked == NO)
			{
			MrBayesPrint ("%s   Could not find outgroup taxon\n", spacer);
			goto errorExit;
			}
		localOutgroupNum = sumtOut;
		}

	for (i=0; i<nNodes; i++)
		{
		p = downPass[i];
		if (p->left != NULL && p->right != NULL)
			if (p->left->marked == YES || p->right->marked == YES)
				p->marked = YES;
		}	

	lft = sumtRoot->left;
	while (lft->left->index != localOutgroupNum && lft->right->index != localOutgroupNum)
		{
		if (lft->left->marked == YES && lft->right->marked == NO)
			{
			m1 = lft->left;
			um1 = lft->right;
			if (m1->left != NULL && m1->right != NULL)
				{
				if (m1->left->marked == YES)
					{
					m2 = m1->left;
					um2 = m1->right;
					}
				else
					{
					m2 = m1->right;
					um2 = m1->left;
					}
				lft->left = m2;
				lft->right = m1;
				m2->anc = m1->anc = lft;
				m1->left = um2;
				m1->right = um1;
				um1->anc = um2->anc = m1;
				m1->marked = NO;
				um1->length += m1->length;
				m2->length *= 0.5;
				m1->length = m2->length;
				}
			else
				{
				MrBayesPrint ("%s   Rooting routine is lost (4)\n", spacer);
				goto errorExit;
				}
			}
		else if (lft->left->marked == NO && lft->right->marked == YES)
			{
			m1 = lft->right;
			um1 = lft->left;
			if (m1->left != NULL && m1->right != NULL)
				{
				if (m1->left->marked == YES)
					{
					m2 = m1->left;
					um2 = m1->right;
					}
				else
					{
					m2 = m1->right;
					um2 = m1->left;
					}
				lft->left = m1;
				lft->right = m2;
				m2->anc = m1->anc = lft;
				m1->left = um1;
				m1->right = um2;
				um1->anc = um2->anc = m1;
				m1->marked = NO;
				um1->length += m1->length;
				m2->length *= 0.5;
				m1->length = m2->length;
				}
			else
				{
				MrBayesPrint ("%s   Rooting routine is lost (5)\n", spacer);
				goto errorExit;
				}
			}
		else
			{
			MrBayesPrint ("%s   Rooting routine is lost (6)\n", spacer);
			goto errorExit;
			}
		}

	/* make certain outgroup is to the right of the root */
	if (sumtRoot->left->left->index == localOutgroupNum)
		{
		m1 = sumtRoot->left->left;
		m2 = sumtRoot->left->right;
		lft = sumtRoot->left;
		lft->left = m2;
		lft->right = m1;
		}

	/* reindex internal nodes of tree */
	i = numTaxa;
	FinishSumtTree (sumtRoot, &i, YES);

	/* free memory */
	free (downPass);
	free (usedMemIndex);
	free (localTaxaFound);
	
	return (NO_ERROR);
	
	errorExit:
		if (downPass)
			free (downPass);
		if (usedMemIndex)
			free (usedMemIndex);
		if (localTaxaFound)
			free (localTaxaFound);

		return (ERROR);

}





void ShowBits (safeLong *p, int nBitsToShow)

{

	int			i;
	safeLong		x, y;
	
	for (i=0; i<nBitsToShow; i++)
		{
		x = 0;
		y = p[i / nBitsInALong];
		x = 1 << (i % nBitsInALong);
		if ((x & y) == 0)
			MrBayesPrint ("0");
		else
			MrBayesPrint ("1");
		if ((i+1) % nBitsInALong == 0)
			MrBayesPrint (" ");
		}

}




int ShowConPhylogram (FILE *fp, int nNodes, PolyNode *root, int screenWidth, int isCalibrated)

{

	int 			i, j, k, nLines, from, to, treeWidth=0, barLength, printExponential,
                    precision, width, newPos, curPos, nTimes, numSpaces, maxLabelLength;
	char			*printLine, *markLine, temp[30], *label;
	MrBFlt			scale, f, scaleBar;
	PolyNode		*p, *q, **allDownPass;

    /* set max label length */
    maxLabelLength = 20;

	/* allocate space for label, printLine and markLine */
	printLine = (char *) calloc ((2*screenWidth+2),sizeof(char));
    label = (char *) calloc (maxLabelLength+1, sizeof(char));
	if (!printLine || !label)
		return ERROR;
	markLine = printLine + screenWidth + 1;

	/* allocate allDownPass */
	allDownPass = (PolyNode **) SafeMalloc (nNodes * sizeof(PolyNode *));
	if (!allDownPass)
		{
		free (printLine);
		return ERROR;
		}

	/* get fresh downpass sequence */
	i = 0;
	GetConDownPass(allDownPass, root, &i);
	
	/* calculate scale */
	scale = 0.0;
	root->f = 0.0;
	for (i=nNodes-2; i>=0; i--)
		{
		p = allDownPass[i];
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
	for (i=0; i<nNodes; i++)
		{
		p = allDownPass[i];
		p->x = (int) (0.5 + (p->f / scale));
		}

	/* calculate y coordinates and lines to print */
	for (i=nLines=0; i<nNodes; i++)
		{
		p = allDownPass[i];
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

		for (j=0; j<nNodes; j++)
			{
			p = allDownPass[j];
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
					sprintf(printLine+to+2,"%s", label);
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

	if (nodeDepthConTree == YES)
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
			
		if (isCalibrated == YES)
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

	free (allDownPass);
	free (printLine);

	return NO_ERROR;

}





void ShowConNodes (int nNodes, PolyNode *root)

{

	int 			i;
	PolyNode		*p;

	/* this is the tree, on a node-by-node basis */
	printf ("root = %d\n", root->index);
	for (i=0; i<2*numTaxa; i++)
		{
		p = &conNodes[i];
		if (!(p->left == NULL && p->sib == NULL && p->anc == NULL))
			{
			printf ("%4d -- %2d ", i, p->index);
			if (p->sib != NULL)
				printf ("(%2d ", p->sib->index);
			else
				printf ("(-- ");
				
			if (p->left != NULL)
				printf ("%2d ", p->left->index);
			else
				printf ("-- ");

			if (p->anc != NULL)
				printf ("%2d)", p->anc->index);
			else
				printf ("--)");
			
			if (p->left == NULL && p->anc != NULL)
				printf ("  %s (%d)\n", p->label, p->x);
			else
				printf (" (%d)\n", p->x);
			}
		}

	return;
}

		
		
		
		
int ShowConTree (FILE *fp, int nNodes, PolyNode *root, int screenWidth, int showSupport)

{

	int 			i, j, k, treeWidth, minBranchLength, maxWidth, isTreeDivided,
					printWidth, nLines, nodesToBePrinted, from, to, maxLabelLength,
                    maxLength;
	char			*printLine, *markLine, temp[20], *label;
	PolyNode		*p=NULL, *q, **allDownPass;

    maxLength = 20;         /* max length of label */
	minBranchLength = 5;    /* min length of branch in tree */
	isTreeDivided = NO;
	
	/* allocate space for printLine, markLine and label */
	printLine = (char *) calloc (maxLength+1+(2*screenWidth+2),sizeof(char));
	if (!printLine)
		return ERROR;
	markLine = printLine + screenWidth + 1;
    label = markLine + screenWidth + 1;

	/* allocate allDownPass */
	allDownPass = (PolyNode **) SafeMalloc (nNodes * sizeof(PolyNode *));
	if (!allDownPass)
		{
		free (printLine);
		return ERROR;
		}

	/* get fresh downpass sequence and internal node indices */
	i = 0;
	GetConDownPass(allDownPass, root, &i);
	for (i=k=0; i<nNodes; i++)
		if (allDownPass[i]->left == NULL)
			k++;
	for (i=0; i<nNodes; i++)
		{
		p = allDownPass[i];
		if (p->left != NULL)
			p->index = k++;
		}
	
	/* calculate max length of labels including taxon index number */
    maxLabelLength = 0;
    for (i=0; i<nNodes; i++)
		{
		p = allDownPass[i];
		if (p->left == NULL)
            {
            j = Label(p,YES,NULL,maxLength);
            if (j > maxLabelLength)
                maxLabelLength = j;
            }
        }

    /* make sure label can hold an interior node index number */
    j = (int) (3.0 + log10((MrBFlt)nNodes)); 
    maxLabelLength = (maxLabelLength > j ? maxLabelLength : j);

	/* calculate remaining screen width for tree
	   and maxWidth in terms of branches */
	treeWidth = screenWidth - maxLabelLength - 1;
	maxWidth = treeWidth / minBranchLength;
	
	/* unmark whole tree */
	for (i=0; i<nNodes; i++)
		allDownPass[i]->mark = 0;
	nodesToBePrinted = nNodes;

	while (nodesToBePrinted > 0)
		{
		/* count depth of nodes in unprinted tree */
		for (i=0; i<nNodes; i++)
			{
			p = allDownPass[i];
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
		for (i=nNodes-2; i>=0; i--)
			{
			p = allDownPass[i];
			if (p->mark == 0 && p->anc->mark == 1)
				{	
				p->mark = 1;
				p->x = (int) (treeWidth - 0.5 - ((treeWidth - 1) * (p->x / (MrBFlt) printWidth)));
				}
			}

		/* calculate y coordinates of nodes to be printed and lines to print */
		for (i=nLines=0; i<nNodes; i++)
			{
			p = allDownPass[i];
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

			for (j=0; j<nNodes; j++)
				{
				p = allDownPass[j];
				if (p->mark != 1 || p->y != i)
					continue;

				/* this branch should be printed
				   add label if the branch is terminal in tree to be printed */
				if (p->left == NULL)
                    {
                    Label (p,YES,label,maxLength);
					sprintf(printLine+treeWidth+1,"%s", label);
                    }
				else if (p->left->mark == 2)
					sprintf(printLine+treeWidth+1,"(%d)", p->index);
				
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
						sprintf(temp, "%d", (int) (p->support + 0.5));
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
							sprintf(temp, "%d", (int) (p->support + 0.5));
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
		
			MrBayesPrintf(fp, "%s\n",printLine);
			}

		/* mark printed branches */
		for (i=0; i<nNodes; i++)
			{
			p = allDownPass[i];
			if (p->mark == 1)
				{
				if (p->anc == NULL)
					p->mark = 2;
				else if (p->anc->mark == 0)
					p->mark = 0;	/* this branch will have to be printed again */
				else
					p->mark = 2;
				}
			}

		}	/* next subtree */
	
	free (allDownPass);
	free (printLine);

	return NO_ERROR;
	
}





void ShowParts (FILE *fp, safeLong *p, int nTaxaToShow)

{

	int			i, flipBits;
	safeLong		x, y;
	
	flipBits = NO;
//	if (isSumtTreeRooted == NO)
//		{
//		y = p[0];
//		x = 1 << (0 % nBitsInALong);
//		if (x & y)
//			flipBits = YES;
//		}

	for (i=0; i<nTaxaToShow; i++)
		{
		x = 0;
		y = p[i / nBitsInALong];
		x = 1 << (i % nBitsInALong);
		if (flipBits == NO)
			{
			if ((x & y) == 0)
				MrBayesPrintf (fp, ".");
			else
				MrBayesPrintf (fp, "*");
			}
		else
			{
			if ((x & y) == 0)
				MrBayesPrintf (fp, "*");
			else
				MrBayesPrintf (fp, ".");
			}
		}

}





void ShowSumtNodes (SumtNode *p, int indent, int isThisTreeRooted)

{

	if (p != NULL)
		{
		MrBayesPrint ("   ");
		if (p->left == NULL && p->right == NULL && p->anc != NULL)
			{
			MrBayesPrint("%*cN %d (l=%d r=%d a=%d) %lf (%s) ", 
			indent, ' ', SumtDex(p), SumtDex(p->left), SumtDex(p->right), SumtDex(p->anc), p->length, p->label);
			}
		else if (p->left != NULL && p->right == NULL && p->anc == NULL)
			{
			if (isThisTreeRooted == NO)
				{
				if (p->label[0] == '\0' || p->label[0] == '\n' || p->label[0] == ' ')
					MrBayesPrint("%*cN %d (l=%d r=%d a=%d) (---) ", 
					indent, ' ', SumtDex(p), SumtDex(p->left), SumtDex(p->right), SumtDex(p->anc));
				else
					MrBayesPrint("%*cN %d (l=%d r=%d a=%d) (%s) ", 
					indent, ' ', SumtDex(p), SumtDex(p->left), SumtDex(p->right), SumtDex(p->anc), p->label);
				}
			else
				{
				MrBayesPrint("%*cN %d (l=%d r=%d a=%d) X.XXXXXX ", 
				indent, ' ', SumtDex(p), SumtDex(p->left), SumtDex(p->right), SumtDex(p->anc));
				}
			}
		else
			{
			if (p->anc != NULL)
				{
				if (p->anc->anc == NULL && isThisTreeRooted == YES)
					MrBayesPrint("%*cN %d (l=%d r=%d a=%d) X.XXXXXX ", 
					indent, ' ', SumtDex(p), SumtDex(p->left), SumtDex(p->right), SumtDex(p->anc));
				else	
					MrBayesPrint("%*cN %d (l=%d r=%d a=%d) %lf ", 
					indent, ' ', SumtDex(p), SumtDex(p->left), SumtDex(p->right), SumtDex(p->anc), p->length);
				}
			}
		MrBayesPrint ("\n");
		ShowSumtNodes (p->left,  indent + 2, isThisTreeRooted);
		ShowSumtNodes (p->right, indent + 2, isThisTreeRooted);
		}
   
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





void SortIndParts (int *item, MrBFlt *assoc, int count, int descendingOrder)

{

	SortIndParts2 (item, assoc, 0, count-1, descendingOrder);

}





void SortIndParts2 (int *item, MrBFlt *assoc, int left, int right, int descendingOrder)

{

	register int	i, j, x;
	MrBFlt			y;

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
				y = (MrBFlt) item[i];
				item[i] = item[j];
				item[j] = (int) y;
				
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
			SortIndParts2 (item, assoc, left, j, descendingOrder);
		if (i < right)
			SortIndParts2 (item, assoc, i, right, descendingOrder);
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
				y = (MrBFlt) item[i];
				item[i] = item[j];
				item[j] = (int) y;
				
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
			SortIndParts2 (item, assoc, left, j, descendingOrder);
		if (i < right)
			SortIndParts2 (item, assoc, i, right, descendingOrder);
		}

}





int SortParts (int *item, int count)

{

	int				i, *tempVect;
	
	SortParts2 (item, 0, count-1);
		
	tempVect = (int *)SafeMalloc((size_t) (numTreePartsFound * sizeof(int)));
	if (!tempVect)
		{
		MrBayesPrint ("%s   Problem allocating tempVect (%d)\n", spacer, numTreePartsFound * sizeof(int));
		goto errorExit;
		}
		
	for (i=0; i<numTreePartsFound; i++)
		{
		tempVect[i] = partOrigOrder[i];
		partOrigOrder[i] = 0;
		}

	for (i=0; i<numTreePartsFound; i++)
		partOrigOrder[tempVect[i]] = i;
		
	/*for (i=0; i<numTreePartsFound; i++)
		printf ("%d -> %d\n", partOrigOrder[i], i);*/

	free (tempVect);
	
	return (NO_ERROR);
	
	errorExit:
		if (tempVect)
			free (tempVect);
		return (ERROR);

}





void SortParts2 (int *item, int left, int right)

{

	register int	i, j, k;
	int				yI;
	safeLong			yL;
	MrBFlt			x, y;

	i = left;
	j = right;
	x = (MrBFlt) item[(left+right)/2];
	do 
		{
		while (item[i] > x && i < right)
			i++;
		while (x > item[j] && j > left)
			j--;
		if (i <= j)
			{
			yI = item[i];
			item[i] = item[j];
			item[j] = yI;
			
			yI = partOrigOrder[i];
			partOrigOrder[i] = partOrigOrder[j];
			partOrigOrder[j] = yI;
			
			if (sumtBrlensDef == YES)
				{
				y = aBrlens[i];
				aBrlens[i] = aBrlens[j];
				aBrlens[j] = y;
				
				y = sBrlens[i];
				sBrlens[i] = sBrlens[j];
				sBrlens[j] = y;
				}

			if (sumtParams.numRuns > 1 && sumtBrlensDef == YES)
				{
				y = aWithinBrlens[i];
				aWithinBrlens[i] = aWithinBrlens[j];
				aWithinBrlens[j] = y;
				
				y = sWithinBrlens[i];
				sWithinBrlens[i] = sWithinBrlens[j];
				sWithinBrlens[j] = y;

				y = sumB[i];
				sumB[i] = sumB[j];
				sumB[j] = y;

				y = sumsqB[i];
				sumsqB[i] = sumsqB[j];
				sumsqB[j] = y;
				}

			for (k = 0; k<taxonLongsNeeded; k++)
				{
				yL = treePartsFound[i*taxonLongsNeeded + k];
				treePartsFound[i*taxonLongsNeeded + k] = treePartsFound[j*taxonLongsNeeded + k];
				treePartsFound[j*taxonLongsNeeded + k] = yL;
				}
				
			if (comparingFiles == YES)
				{
				yI = numFoundOfThisPart1[i];
				numFoundOfThisPart1[i] = numFoundOfThisPart1[j];
				numFoundOfThisPart1[j] = yI;
				
				yI = numFoundOfThisPart2[i];
				numFoundOfThisPart2[i] = numFoundOfThisPart2[j];
				numFoundOfThisPart2[j] = yI;
				}

			if (sumtParams.numRuns > 1)
				{
				for (k=0; k<sumtParams.numRuns; k++)
					{
					yI = numFoundInRunOfPart[k][i];
					numFoundInRunOfPart[k][i] = numFoundInRunOfPart[k][j];
					numFoundInRunOfPart[k][j] = yI;
					}
				}
				
			i++;
			j--;
			}
		} while (i <= j);
	if (left < j)
		SortParts2 (item, left, j);
	if (i < right)
		SortParts2 (item, i, right);

}





int SumtDex (SumtNode *p)

{

	return (p == NULL) ? -1 : p->index;

}





int TreeProb (void)

{

	int			i, j, n, num, targetNode, nBits, nextConNode, isCompat, 
				localOutgroupNum, origPartNum, reorderedPartNum, *tempTreeNum=NULL, *tempNumOfTree=NULL,
				nInSets[5];
	safeLong		x, *mask, *partition, *ingroupPartition, *outgroupPartition=NULL;
	MrBFlt		treeProb, cumTreeProb;
	char		tempName[100];
	PolyNode	*cp, *q, *r, *ql, *pl;
	
	/* check if we need to do this */
	if (sumtParams.calcTrprobs == NO)
		return (NO_ERROR);

	MrBayesPrint ("%s   Calculating tree probabilities...\n\n", spacer);

	/* check that we have at least three species */
	j = 0;
	for (i=0; i<numTaxa; i++)
		if (sumTaxaFound[i] == YES)
			j++;
	if (j < 3)
		{
		MrBayesPrint ("%s   Too few taxa included to show tree probabilities\n", spacer);
		goto errorExit;
		}
		
	/* sort trees, from most probable to least probable */
	tempTreeNum = (int *)SafeMalloc((size_t) (numFullTreesFound * sizeof(int)));
	if (!tempTreeNum)
		{
		MrBayesPrint ("%s   Problem allocating tempTreeNum (%d)\n", spacer, numFullTreesFound * sizeof(int));
		goto errorExit;
		}
	tempNumOfTree = (int *)SafeMalloc((size_t) (numFullTreesFound * sizeof(int)));
	if (!tempNumOfTree)
		{
		MrBayesPrint ("%s   Problem allocating tempNumOfTree (%d)\n", spacer, numFullTreesFound * sizeof(int));
		goto errorExit;
		}
	for (i=0; i<numFullTreesFound; i++)
		{
		tempTreeNum[i] = i;
		tempNumOfTree[i] = numOfThisFullTree[i];
		}
	SortInts (tempNumOfTree, tempTreeNum, numFullTreesFound, YES);
	
	/* Set the outgroup. Remember that the outgroup number goes from 0 to numTaxa-1.
	   The outgroup may have been deleted, so we should probably set localOutgroupNum
	   to reflect this. */
	j = 0;
	localOutgroupNum = 0;
	for (i=0; i<numTaxa; i++)
		{
		if (sumTaxaFound[i] == YES)
			{
			if (i == outGroupNum)
				{
				localOutgroupNum = j;
				break;
				}
			j++;
			}
		}
	
	/* note that numIncludedTaxa is initialized in ReorderParts */

	/* First allocate some stuff for the trees. We use the same routines that
	   we used when making consensus trees. However, all of the trees should
	   be strictly bifurcating. */
	if (memAllocs[ALLOC_CONNODES] == YES)
		{
		MrBayesPrint ("%s   conNodes is already allocated\n", spacer);
		goto errorExit;
		}
	conNodes = (PolyNode *)SafeMalloc((size_t) (2 * numTaxa * sizeof(PolyNode)));
	if (!conNodes)
		{
		MrBayesPrint ("%s   Could not allocate conNodes\n", spacer);
		goto errorExit;
		}
	memAllocs[ALLOC_CONNODES] = YES;
	for (i=0; i<2*numTaxa; i++)
		{
		conNodes[i].left = conNodes[i].sib = conNodes[i].anc = NULL;
		conNodes[i].x = conNodes[i].y = conNodes[i].index = conNodes[i].mark = 0;
		conNodes[i].length = conNodes[i].support = conNodes[i].f = 0.0;
		}

	if (memAllocs[ALLOC_OUTPART] == YES)
		{
		MrBayesPrint ("%s   outgroupPartition is already allocated\n", spacer);
		goto errorExit;
		}
	outgroupPartition = (safeLong *) calloc (3 * taxonLongsNeeded, sizeof(safeLong));
	if (!outgroupPartition)
		{
		MrBayesPrint ("%s   Could not allocate outgroupPartition\n", spacer);
		goto errorExit;
		}
	ingroupPartition = outgroupPartition + taxonLongsNeeded;
	mask = ingroupPartition + taxonLongsNeeded;
	memAllocs[ALLOC_OUTPART] = YES;

	/* Set mask - needed to trim last element in partition.
	   This could be done when a new matrix is read in
	   and adjusted when taxa are deleted or restored. */
	for (i=0; i<numIncludedTaxa; i++)
		SetBit (i, mask);

	/* Set ingroup and outgroup partitions.
	   This could be done when a new matrix is read in
	   and adjusted when an outgroup command is issued.
	   This mechanism allows multiple taxa in outgroup. */
	x = 1;
	x <<= (localOutgroupNum) % nBitsInALong;
	i = (localOutgroupNum) / nBitsInALong;
	outgroupPartition[i] = x;
	for (i = 0; i < taxonLongsNeeded; i++)
		ingroupPartition[i] = outgroupPartition[i];
	FlipBits (ingroupPartition, taxonLongsNeeded, mask);	
	/*ShowBits (&outgroupPartition[0], numIncludedTaxa);
	MrBayesPrint (" <- outgroupPartition\n");
	ShowBits (&ingroupPartition[0], numIncludedTaxa);
	MrBayesPrint (" <- ingroupPartition\n");*/
	
	/* now, resolve each tree in the list of trees */
	cumTreeProb = 0.0;
	nInSets[0] = nInSets[1] = nInSets[2] = nInSets[3] = nInSets[4] = 0;
	for (num=0; num<numFullTreesFound; num++)   /* loop over all of the trees that were found */
		{
		
		/* figure out which tree we want */
		n = tempTreeNum[num];

		/* initialize terminal consensus nodes */
		j = 0;
		for (i=0; i<numTaxa; i++)
			{
			if (sumTaxaFound[i] == YES)
				{
				if (GetNameFromString (taxaNames, tempName, i+1) == ERROR)
					{
					MrBayesPrint ("%s   Error getting taxon names \n", spacer);
					return (ERROR);
					}
				conNodes[j].left = NULL;
				conNodes[j].sib = NULL;
				conNodes[j].index = j;
				strcpy (conNodes[j].label, tempName);
				j++;
				}
			}
		for (i=numIncludedTaxa; i<2*numIncludedTaxa; i++)
			{
			conNodes[j].left = NULL;
			conNodes[j].sib = NULL;
			conNodes[j].index = j;
			j++;
			}
			
		/* create bush 
		   ->x counts number of subtended terminals 
		   make sure conRoot->left is in outgroup */
		conRoot = &conNodes[numIncludedTaxa];
		conRoot->anc = conRoot->sib = NULL;
		conRoot->x = numIncludedTaxa;
		j = FirstTaxonInPartition(outgroupPartition, taxonLongsNeeded);
		conRoot->left = cp = &conNodes[j];
		cp->anc = conRoot;
		cp->x = 1;
		for (i=0; i<numIncludedTaxa; i++)
			{
			if (i != j)
				{
				cp->sib = &conNodes[i];
				cp = cp->sib;
				cp->anc = conRoot;
				cp->x = 1;
				}
			}
		cp->sib = NULL;

		/* Resolve bush according to partitions.
		   Partitions may include incompatible ones. */
		nextConNode = numIncludedTaxa + 1;
		if (isSumtTreeRooted == YES)
			targetNode = 2 * numIncludedTaxa - 2;
		else
			targetNode = 2 * numIncludedTaxa - 3;
			
		for (i=0; i<numTreeParts; i++) /* loop over partitions for this tree */
			{
			/* get partition */
			origPartNum = fullTreePartIds[n * 2 * numTaxa + i];
			reorderedPartNum = partOrigOrder[origPartNum];
			partition = &treePartsFound[reorderedPartNum*taxonLongsNeeded];

			/* flip bits if necessary */
			if (isSumtTreeRooted == NO)
				{
				if (!IsPartNested(partition, ingroupPartition, taxonLongsNeeded) && !IsPartCompatible(partition, ingroupPartition, taxonLongsNeeded))
					FlipBits(partition, taxonLongsNeeded, mask);
				}
			
			/* count bits in this partition */
			for (j=nBits=0; j<taxonLongsNeeded; j++)
				{
				x = partition[j];
				for (x = partition[j]; x != 0; x &= (x - 1))
					nBits++;
				}
				
			/* flip this partition if it leaves single outgroup outside */
			if (nBits == numIncludedTaxa - 1  && isSumtTreeRooted == NO)
				{
				nBits = 1;
				FlipBits(partition, taxonLongsNeeded, mask);
				}

			/*ShowBits (&partition[0], numIncludedTaxa);
			MrBayesPrint (" <- partition (%d)\n", nBits);*/

			if (nBits > 1 && nBits < numIncludedTaxa)  /* this is an informative partition */
				{
				/* find anc of partition */
				j = FirstTaxonInPartition (partition, taxonLongsNeeded);
				for (cp = &conNodes[j]; cp!=NULL; cp = cp->anc)
					if (cp->x > nBits)
						break;

				/* Do not include if incompatible with ancestor
				   or any of descendants.
				   Do not check terminals or root because it is
				   redundant and partitions have not been set for those. */
				isCompat = YES;
				if (cp->anc != NULL && !IsPartCompatible(partition, cp->partition, taxonLongsNeeded))
					isCompat = NO;
				for (q=cp->left; q!=NULL; q=q->sib)
					{
					if (q->x > 1 && !IsPartCompatible(q->partition, partition, taxonLongsNeeded))
						isCompat = NO;
					if (isCompat == NO)
						{
						MrBayesPrint ("%s   Found an incompatible partition in the tree (1 %d %d)\n", spacer, num, n);
						goto errorExit;
						}
					}
				if (isCompat == NO)
					{
					MrBayesPrint ("%s   Found an incompatible partition in the tree (2)\n", spacer);
					goto errorExit;
					}

				/* check for number of nodes in tree */
				if (nextConNode > targetNode)
					{
					MrBayesPrint ("%s   Too many nodes in tree\n", spacer);
					goto errorExit;
					}
			
				/* set new node */
				q = &conNodes[nextConNode++];
				q->x = nBits;
				q->partition = partition;

				/* go through descendants of anc */
				ql = pl = NULL;
				for (r=cp->left; r!=NULL; r=r ->sib)
					{
					/* test if r is in the new partition or not */
					if ((r->x > 1 && IsPartNested(r->partition, partition, taxonLongsNeeded)) || (r->x == 1 && (partition[r->index / nBitsInALong] & (1 << (r->index % nBitsInALong))) != 0))
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
							cp->left = r;
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
				q->anc = cp;
				}
			else
				/* singleton partition */
				{
				j = FirstTaxonInPartition(partition, taxonLongsNeeded);
				q = &conNodes[j];
				}
			}

		/* get probability of tree */
		treeProb = (MrBFlt)tempNumOfTree[num] / (MrBFlt) (sumtParams.numRuns * numSumTreesSampled);
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
		if (sumtParams.showSumtTrees == YES)
			{
			MrBayesPrint ("\n%s   Tree %d (p = %1.3lf, P = %1.3lf):\n\n", spacer, num+1, treeProb, cumTreeProb);
			ShowConTree (stdout, nextConNode, conRoot, 80, NO);
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
			j = 0;
			for (i=0; i<numTaxa; i++)
				{
				if (sumTaxaFound[i] == YES)
					{
					if (GetNameFromString (taxaNames, tempName, i+1) == ERROR)
						{
						MrBayesPrint ("%s   Error getting taxon names \n", spacer);
						return (ERROR);
						}
					if (j+1 == numIncludedTaxa)
						MrBayesPrintf (fpTrees, "   %2d %s;\n", j+1, tempName);
					else
						MrBayesPrintf (fpTrees, "   %2d %s,\n", j+1, tempName);
					j++;
					}
				}
			
			}
		MrBayesPrintf (fpTrees, "   tree tree_%d [p = %1.3lf, P = %1.3lf] = [&W %1.6lf] ", num+1, treeProb, cumTreeProb, treeProb);
		WriteTree (conRoot, fpTrees);
		MrBayesPrintf (fpTrees, ";\n");
		if (num == numFullTreesFound - 1)
			MrBayesPrintf (fpTrees, "end;\n");
	
		}	
		
	/* print out general information on credible sets of trees */
	MrBayesPrint ("%s   Credible sets of trees (%d trees sampled):\n", spacer, nInSets[0] + nInSets[1] + nInSets[2] + nInSets[3] + nInSets[4]);
	if (nInSets[0] != 0)
		MrBayesPrint ("%s      50 %% credible set contains %d trees\n", spacer, nInSets[0] + 1);
	if (nInSets[0] + nInSets[1] != 0)
		MrBayesPrint ("%s      90 %% credible set contains %d trees\n", spacer, nInSets[0] + nInSets[1] + 1);
	if (nInSets[0] + nInSets[1] + nInSets[2] != 0)
		MrBayesPrint ("%s      95 %% credible set contains %d trees\n", spacer, nInSets[0] + nInSets[1] + nInSets[2] + 1);
	MrBayesPrint ("%s      99 %% credible set contains %d trees\n", spacer, nInSets[0] + nInSets[1] + nInSets[2] + nInSets[3] + 1);

	
	/* free memory and file pointers */
	if (memAllocs[ALLOC_CONNODES] == YES)
		{
		free (conNodes);
		memAllocs[ALLOC_CONNODES] = NO;
		}
	if (memAllocs[ALLOC_OUTPART] == YES)
		{
		free (outgroupPartition);
		memAllocs[ALLOC_OUTPART] = NO;
		}
	free (tempTreeNum);
	free (tempNumOfTree);
		
	return (NO_ERROR);
	
	errorExit:
		if (memAllocs[ALLOC_CONNODES] == YES)
			{
			free (conNodes);
			memAllocs[ALLOC_CONNODES] = NO;
			}
		if (memAllocs[ALLOC_OUTPART] == YES)
			{
			free (outgroupPartition);
			memAllocs[ALLOC_OUTPART] = NO;
			}
		if (tempTreeNum)
			free (tempTreeNum);
		if (tempNumOfTree)
			free (tempNumOfTree);
		return (ERROR);
	
}





void WriteConTree (PolyNode *p, FILE *fp, int showSupport)

{

	PolyNode		*q;

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
		if (sumtBrlensDef == YES)
			fprintf (fp, "%d:%lf", p->index+1, p->length);
		else
			fprintf (fp, "%d", p->index+1);
		}
		
	if (p->sib == NULL && p->anc != NULL)
		{
		if (p->anc->anc != NULL)
			{
			if (sumtBrlensDef == YES && showSupport == NO)
				fprintf (fp, "):%lf", p->anc->length); 
			else if (sumtBrlensDef == NO && showSupport == YES)
				fprintf (fp, ")%1.2lf", p->anc->support/100.0); 
			else if (sumtBrlensDef == YES && showSupport == YES)
				fprintf (fp, ")%1.2lf:%lf", p->anc->support/100.0, p->anc->length);
			else
				fprintf (fp, ")");
			}
		else
			fprintf (fp, ")");
		}

}





void WriteTree (PolyNode *p, FILE *fp)

{

	PolyNode		*q;

	if (p->anc != NULL)
		if (p->anc->left == p)
			fprintf (fp, "(");

	for (q = p->left; q != NULL; q = q->sib)
		{
		if (q->anc->left != q)  /* Note that q->anc always exists (it is p) */
			fprintf (fp, ",");
		WriteTree (q, fp);
		}
	if (p->left == NULL)
		{
		fprintf (fp, "%d", p->index+1);
		}
		
	if (p->sib == NULL && p->anc != NULL)
		{
		fprintf (fp, ")");
		}

}
