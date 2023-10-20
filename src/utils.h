#ifndef UTILS_H_
#define UTILS_H_

/* M_PI and M_PI_2 not part of standard C */
#ifndef M_PI
#    define M_PI    3.14159265358979323846264338327950288
#endif
#ifndef M_PI_2
#    define M_PI_2  1.57079632679489661923132169163975144
#endif

struct MrBComplex
{
    MrBFlt re;
    MrBFlt im;
};
typedef struct MrBComplex MrBComplex;

typedef struct
    {
    MrBFlt     mean;
    MrBFlt     median;
    MrBFlt     lower;
    MrBFlt     upper;
    MrBFlt     var;
    MrBFlt     PSRF;
    MrBFlt     avrESS;
    MrBFlt     minESS;
    }
    Stat;


/* For explanation why the following two macros exists, see
 * http://stackoverflow.com/questions/38569628/calling-a-free-wrapper-dereferencing-type-punned-pointer-will-break-strict-al
 */

#define SAFEFREE(ptr) (ptr = SafeFree(ptr))
#define ALIGNEDSAFEFREE(ptr) (ptr = AlignedSafeFree(ptr))

int      AddBitfield (BitsLong ***list, int listLen, int *set, int setLen);
#if defined (SSE_ENABLED)
void    *AlignedMalloc (size_t size, size_t alignment);
void    *AlignedSafeFree (void *ptr);
#endif
int      AreBitfieldsEqual (BitsLong *p, BitsLong *q, int length);
int      Bit (int n, BitsLong *p);
void     ClearBit (int i, BitsLong *bits);
void     ClearBits (BitsLong *bits, int nLongs);
void     CopyBits (BitsLong *dest, BitsLong *source, int nLongs);
int      CopyResults (FILE *toFile, char *fromFileName, long long lastGen);
int      CopyProcessSsFile (FILE *toFile, char *fromFileName, int lastStep, MrBFlt *marginalLnLSS, MrBFlt *splitfreqSS);
int      CopyTreeResults (FILE *toFile, char *fromFileName, long long lastGen, int *treeNum);
int      FirstTaxonInPartition (BitsLong *partition, int length);
long     FirstTree (FILE *fp, char *lineBuf, int longestLine);
int      Flip01 (int x);
void     FlipBits (BitsLong *partition, int length, BitsLong *mask);
void     FlipOneBit (int n, BitsLong *p);
int      FromGrowthFxnToIndex (int *growthFxn);
void     FromIndexToGrowthFxn (int index, int *growthFxn);
void     GetIntSummary (int **vals, int nRows, int *rowCount, Stat *theStats, int HPD);
int      GetKFromGrowthFxn (int *growthFxn);
void     GetSummary (MrBFlt **vals, int nRows, int *rowCount, Stat *theStats, int HPD);
int      HarmonicArithmeticMeanOnLogs (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *harm_mean);
int      IsBitSet (int i, BitsLong *bits);
int      IsConsistentWith (const char *token, const char *expected);
int      IsPartNested (BitsLong *smaller, BitsLong *larger, int length);
int      IsPartCompatible (BitsLong *smaller, BitsLong *larger, int length);
int      IsSectionEmpty (BitsLong *bitField1, BitsLong *bitField2, int length);
int      IsUnionEqThird (BitsLong *bitField1, BitsLong *bitField2, BitsLong *bitField3, int length);
long     LastBlock (FILE *fp, char *lineBuf, int longestLine);
int      LineTermType (FILE *fp);
MrBFlt   LnDirichlet (MrBFlt *alphai, MrBFlt *xi, int lengthi);
int      LongestLine (FILE *fp);
void     LowerUpperMedian (MrBFlt *vals, int nVals, MrBFlt *lower, MrBFlt *upper, MrBFlt *median);
void     LowerUpperMedianHPD (MrBFlt *vals, int nVals, MrBFlt *lower, MrBFlt *upper, MrBFlt *median);
void     MeanVariance (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *var);
void     MeanVarianceLog (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *var, MrBFlt *varEst);
int      NextTaxonInPartition (int currentTaxon, BitsLong *partition, int length);
int      NBits (int x);
int      NumBits (BitsLong *x, int len);
char    *MbPrintNum (MrBFlt num);
void     MrBayesPrint (char *format, ...);
void     MrBayesPrintf (FILE *f, char *format, ...);
FILE    *OpenBinaryFileR (char *name);
FILE    *OpenTextFileA (char *name);
FILE    *OpenTextFileR (char *name);
FILE    *OpenTextFileRQuait (char *name);
FILE    *OpenTextFileW (char *name);
MrBFlt   PotentialScaleReduction (MrBFlt **vals, int nRows, int *count);
void     EstimatedSampleSize (MrBFlt **vals, int nRuns, int *count, MrBFlt *returnESS);
void    *SafeCalloc (size_t n, size_t s);
int      SafeFclose (FILE **fp);
void    *SafeFree (void *ptr);
void    *SafeMalloc (size_t s);
void    *SafeRealloc (void *ptr, size_t s);
char    *SafeStrcat (char **target, const char *source);
char    *SafeStrcpy (char **target, const char *source);
void     SetBit (int i, BitsLong *bits);
void     SortInts (int *item, int *assoc, int count, int descendingOrder);
void     SortInts2 (int *item, int *assoc, int left, int right, int descendingOrder);
void     SortMrBFlt_Asc (MrBFlt *item, int left, int right);
void     SortMrBFlt_Des (MrBFlt *item, int left, int right);
int      StrCmpCaseInsensitiveLen (const char *s, const char *t, size_t len);
int      StrCmpCaseInsensitive (char *s, char *t);
void     StripComments (char *s);
FILE    *TestOpenTextFileR (char *name);
void     UpdateGrowthFxn (int *growthFxn);
int      UpperTriangIndex (int i, int j, int size);
int      WantTo (const char *msg);

/* tree utility functions */
int       AddToTreeList (TreeList *treeList, Tree *tree);
Tree     *AllocateTree (int numTaxa);
Tree     *AllocateFixedTree (int numTaxa, int isRooted);
int       AllocateTreePartitions (Tree *t);
PolyTree *AllocatePolyTree (int numTaxa);
int       AllocatePolyTreePartitions (PolyTree *pt);
int       AllocatePolyTreeRelClockParams (PolyTree *pt, int nBSets, int nESets);
int       AreTopologiesSame (Tree *t1, Tree *t2);
int       AreTreesSame (Tree *t1, Tree *t2);
int       BuildConstraintTree (Tree *t, PolyTree *pt, char **localTaxonNames);
int       BuildRandomRTopology (Tree *t, RandLong *seed);
int       BuildRandomUTopology (Tree *t, RandLong *seed);
int       CheckConstraints (Tree *t);
int       CheckSetConstraints (Tree *t);
void      ColorClusters (TreeNode *p, int *index);
void      CopySubtreeToTree (Tree *subtree, Tree *t);
int       CopyToPolyTreeFromPolyTree (PolyTree *to, PolyTree *from);
int       CopyToSpeciesTreeFromPolyTree (Tree *to, PolyTree *from);
int       CopyToTreeFromPolyTree (Tree *to, PolyTree *from);
void      CopyPolyNodes (PolyNode *p, PolyNode *q, int nLongsNeeded);
int       CopyToTreeFromTree (Tree *to, Tree *from);
void      CopyTreeNodes (TreeNode *p, TreeNode *q, int nLongsNeeded);
void      CopyTreeToSubtree (Tree *t, Tree *subtree);
int       Deroot (PolyTree *pt);
void      EraseTreeList (TreeList *treeList);
void      findAllowedClockrate (Tree *t, MrBFlt *minClockRate, MrBFlt *maxClockRate);
void      FreePolyTree (PolyTree *pt);
void      FreePolyTreePartitions (PolyTree *pt);
void      FreePolyTreePopSizeParams (PolyTree *pt);
void      FreePolyTreeRelClockParams (PolyTree *pt);
void      FreeTree (Tree *t);
void      FreeTreePartitions (Tree *pt);
void      GetDatedNodeDepths (TreeNode *p, MrBFlt *nodeDepths);
void      GetDatedNodes (TreeNode *p, TreeNode **datedNodes);
void      GetDownPass (Tree *t);
void      GetNodeDownPass (Tree *t, TreeNode *p, int *i, int *j);
void      GetPolyAges (PolyTree *t);
void      GetPolyDepths (PolyTree *t);
void      GetPolyDownPass (PolyTree *t);
void      GetPolyNodeDownPass (PolyTree *t, PolyNode *p, int *i, int *j);
int       GetRandomEmbeddedSubtree (Tree *t, int nTerminals, RandLong *seed, int *nEmbeddedTrees);
int       GetFromTreeList (TreeList *treeList, Tree *tree);
int       InitBrlens (Tree *t, MrBFlt v);
int       InitCalibratedBrlens (Tree *t, MrBFlt minLength, RandLong *seed);
int       InitClockBrlens (Tree *t);
int       IsCalibratedClockSatisfied (Tree *t,MrBFlt *minClockRate,MrBFlt *maxClockRate , MrBFlt tol);
int       IsClockSatisfied (Tree *t, MrBFlt tol);
int       IsTreeConsistent (Param *param, int chain, int state);
int       LabelTree (Tree *t, char **taxonNames);
void      Mark (TreeNode *p);
void      MarkDistance (TreeNode *p, int YESorNO, int dist, int *n);
void      MarkUnconstrained (TreeNode *p);
int       MoveCalculationRoot (Tree *t, int outgroup);
int       MovePolyCalculationRoot (PolyTree *t, int outgroup);
int       NumConstrainedTips (TreeNode *p);
int       NumDatedTips (TreeNode *p);
void      OrderTips (PolyTree *t);
void      PrintNewick (char **s, int *len, Tree *t);
void      PrintNodes (Tree *t);
void      PrintPolyNodes (PolyTree *pt);
void      PrintTranslateBlock (FILE *fp, Tree *t);
int       PrunePolyTree (PolyTree *pt);
int       RandPerturb (Tree *t, int nPert, RandLong *seed);
int       RandResolve (Tree *tt, PolyTree *t, RandLong *seed, int destinationIsRooted);
int       ResetBrlensFromTree (Tree *tree, Tree *vTree);
void      ResetIntNodeIndices (PolyTree *t);
void      ResetPolyTree (PolyTree *t);
void      ResetPolyTreePartitions (PolyTree *pt);
void      ResetPolyTreeRelClockParams (PolyTree *pt);
int       ResetRootHeight (Tree *t, MrBFlt rootHeight);
void      ResetTipIndices (PolyTree *pt);
int       ResetTopology (Tree *t, char *s);
int       ResetTopologyFromTree (Tree *tree, Tree *top);
int       ResetTopologyFromPolyTree (Tree *tree, PolyTree *top);
void      ResetTreePartitions (Tree *t);
int       RetrieveRTopology (Tree *t, int *order);
int       RetrieveRTree (Tree *t, int *order, MrBFlt *brlens);
int       RetrieveRTreeWithIndices (Tree *t, int *order, MrBFlt *brlens);
int       RetrieveUTopology (Tree *t, int *order);
int       RetrieveUTree (Tree *t, int *order, MrBFlt *brlens);
void      SetDatedNodeAges (Param* param, int chain, int state);
void      SetNodeDepths (Tree *t);
int       SetTreeNodeAges (Param *param, int chain, int state);
int       ShowPolyNodes (PolyTree *pt);
int       ShowTree (Tree *t);
int       StoreRPolyTopology (PolyTree *t, int *order);
int       StoreRPolyTree (PolyTree *t, int *order, MrBFlt *brlens);
int       StoreRTopology (Tree *t, int *order);
int       StoreRTree (Tree *t, int *order, MrBFlt *brlens);
int       StoreRTreeWithIndices (Tree *t, int *order, MrBFlt *brlens);
int       StoreUPolyTopology (PolyTree *t, int *order);
int       StoreUPolyTree (PolyTree *t, int *order, MrBFlt *brlens);
int       StoreUTopology (Tree *t, int *order);
int       StoreUTree (Tree *t, int *order, MrBFlt *brlens);
MrBFlt    TreeLen (Tree *t);
void      Unmark (TreeNode *p);
void      UpdateTreeWithClockrate (Tree *t, MrBFlt clockRate);
void      WriteEventTree (TreeNode *p, int chain, Param *param);
void      WriteEventTreeToPrintString (TreeNode *p, int chain, Param *param, int printAll);
void      WriteNoEvtTreeToPrintString (TreeNode *p, int chain, Param *param, int showBrlens, int isRooted);
void      WriteEvolTree (TreeNode *p, int chain, Param *param);
void      WriteTopologyToFile (FILE *fp, TreeNode *p, int isRooted);

/* math utility functions */
MrBComplex **AllocateSquareComplexMatrix (int dim);
MrBFlt  **AllocateSquareDoubleMatrix (int dim);
int     **AllocateSquareIntegerMatrix (int dim);
int       AutodGamma (MrBFlt *M, MrBFlt rho, int K);
void      BetaBreaks (MrBFlt alpha, MrBFlt beta, MrBFlt *values, int K);
MrBFlt    BetaQuantile (MrBFlt alpha, MrBFlt beta, MrBFlt x);
void      CalcCijk (int dim, MrBFlt *c_ijk, MrBFlt **u, MrBFlt **v);
void      CopyComplexMatrices (int dim, MrBComplex **from, MrBComplex **to);
void      CopyDoubleMatrices (int dim, MrBFlt **from, MrBFlt **to);
void      DirichletRandomVariable (MrBFlt *alp, MrBFlt *z, int n, RandLong *seed);
int       DiscreteGamma (MrBFlt *rK, MrBFlt alfa, MrBFlt beta, int K, int median);
void      FreeSquareComplexMatrix (MrBComplex **m);
void      FreeSquareDoubleMatrix (MrBFlt **m);
void      FreeSquareIntegerMatrix (int **m);
int       GetEigens (int dim, MrBFlt **q, MrBFlt *eigenValues, MrBFlt *eigvalsImag, MrBFlt **eigvecs, MrBFlt **inverseEigvecs, MrBComplex **Ceigvecs, MrBComplex **CinverseEigvecs);
MrBFlt    LnFactorial (int value);
MrBFlt    LnGamma (MrBFlt alp);
MrBFlt    LnPriorProbExponential (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbExponential_Param_Mean (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbFix (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbGamma (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbGamma_Param_Mean_Sd (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbLognormal (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbLognormal_Param_Mean_Sd (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbNormal (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbOffsetExponential (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbOffsetExponential_Param_Offset_Mean (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbOffsetGamma (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbOffsetGamma_Param_Offset_Mean_Sd (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbOffsetLognormal (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbOffsetLognormal_Param_Offset_Mean_Sd (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbTruncatedNormal (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbTruncatedNormal_Param_Trunc_Mean_Sd (MrBFlt val, MrBFlt *params);
MrBFlt    LnPriorProbUniform (MrBFlt val, MrBFlt *params);
MrBFlt    LnProbRatioExponential (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioExponential_Param_Mean (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioFix (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioGamma (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioGamma_Param_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioLognormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioLognormal_Param_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioNormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioOffsetExponential (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioOffsetExponential_Param_Offset_Mean (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioOffsetGamma (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioOffsetGamma_Param_Offset_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioOffsetLognormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioOffsetLognormal_Param_Offset_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioTruncatedNormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioTruncatedNormal_Param_Trunc_Mean_Sd (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbRatioUniform (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt    LnProbGamma (MrBFlt alpha, MrBFlt beta, MrBFlt x);
MrBFlt    LnProbTruncGamma (MrBFlt alpha, MrBFlt beta, MrBFlt x, MrBFlt min, MrBFlt max);
MrBFlt    LnProbLogNormal (MrBFlt mu, MrBFlt sigma, MrBFlt x);
MrBFlt    LnProbLogNormal_Mean_LogVar (MrBFlt mean, MrBFlt sigma2, MrBFlt x);
MrBFlt    LnProbLogNormal_Mean_Var (MrBFlt mean, MrBFlt var, MrBFlt x);
MrBFlt    LnProbNormal (MrBFlt mu, MrBFlt sigma, MrBFlt x);
MrBFlt    LnRatioLogNormal (MrBFlt mu, MrBFlt sigma, MrBFlt xNew, MrBFlt xOld);
MrBFlt    LogNormalRandomVariable (MrBFlt mean, MrBFlt var, RandLong *seed);
MrBFlt    MaximumValue (MrBFlt x, MrBFlt y);
MrBFlt    MinimumValue (MrBFlt x, MrBFlt y);
void      MultiplyMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result);
int       MultiplyMatrixNTimes (int dim, MrBFlt **Mat, int power, MrBFlt **Result);
MrBFlt    PointNormal (MrBFlt prob);
MrBFlt    PsiGammaLnProb (MrBFlt alpha, MrBFlt value);
MrBFlt    PsiGammaLnRatio (MrBFlt alpha, MrBFlt numerator, MrBFlt denominator);
MrBFlt    PsiGammaRandomVariable (MrBFlt alpha, RandLong *seed);
MrBFlt    QuantileGamma (MrBFlt x, MrBFlt alfa, MrBFlt beta);
MrBFlt    RandomNumber (RandLong *seed);
MrBFlt    QuantileLogNormal (MrBFlt prob, MrBFlt mu, MrBFlt sigma);
int       DiscreteLogNormal (MrBFlt *rK, MrBFlt sigma, int K, int median);
MrBFlt    LogNormalPoint (MrBFlt x, MrBFlt mu, MrBFlt sigma);

/* qsort utility function */
int       cmpMrBFlt(const void *a, const void *b);

#endif  /* UTILS_H_ */
