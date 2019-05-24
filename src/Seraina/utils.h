#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>


typedef struct
	{
    MrBFlt          mean;
    MrBFlt          median;
    MrBFlt          lower;
    MrBFlt          upper;
    MrBFlt          var;
    MrBFlt          PSRF;
    MrBFlt          avrESS;
    MrBFlt          minESS;
	}
	Stat;


int      AddBitfield (SafeLong ***list, int listLen, int *set, int setLen);
#if defined (SSE_ENABLED)
void     AlignedSafeFree(void **ptr);
#endif
void     ClearBit (int i, SafeLong *bits);
void     ClearBits (SafeLong *bits, int nLongs);
int      CopyResults (FILE *toFile, char *fromFileName, int lastGen);
int      CopyProcessSsFile (FILE *toFile, char *fromFileName, int lastStep, MrBFlt *marginalLnLSS, MrBFlt *splitfreqSS);
int      CopyTreeResults (FILE *toFile, char *fromFileName, int lastGen, int *treeNum);
int      FirstTaxonInPartition (SafeLong *partition, int length);
SafeLong FirstTree (FILE *fp, char *lineBuf, int longestLine);
int      Flip01 (int x);
void     FlipBits (SafeLong *partition, int length, SafeLong *mask);
void     FlipOneBit (int n, SafeLong *p);
int      FromGrowthFxnToIndex(int *growthFxn);
void     FromIndexToGrowthFxn(int index, int *growthFxn);
void     GetIntSummary (int **vals, int nRows, int *rowCount, Stat *theStats, int HPD);
int      GetKFromGrowthFxn(int *growthFxn);
void     GetSummary (MrBFlt **vals, int nRows, int *rowCount, Stat *theStats, int HPD);
int      HarmonicArithmeticMeanOnLogs (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *harm_mean);
int      IsBitSet (int i, SafeLong *bits);
int      IsConsistentWith (const char *token, const char *expected);
int      IsPartNested (SafeLong *smaller, SafeLong *larger, int length);
int      IsPartCompatible (SafeLong *smaller, SafeLong *larger, int length);
int      IsSectionEmpty (SafeLong *bitField1, SafeLong *bitField2, int length);
int      IsUnionEqThird (SafeLong *bitField1, SafeLong *bitField2, SafeLong *bitField3, int length);
SafeLong LastBlock (FILE *fp, char *lineBuf, int longestLine);
int		 LineTermType (FILE *fp);
int      LongestLine (FILE *fp);
void     LowerUpperMedian (MrBFlt *vals, int nVals, MrBFlt *lower, MrBFlt *upper, MrBFlt *median);
void     LowerUpperMedianHPD (MrBFlt *vals, int nVals, MrBFlt *lower, MrBFlt *upper, MrBFlt *median);
void     MeanVariance (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *var);
void MeanVarianceLog (MrBFlt *vals, int nVals, MrBFlt *mean, MrBFlt *var, MrBFlt *varEst );
int      NextTaxonInPartition (int currentTaxon, SafeLong *partition, int length);
int      NumBits (SafeLong *x, int len);
char    *MbPrintNum (MrBFlt num);
void     MrBayesPrint (char *format, ...);
void     MrBayesPrintf (FILE *f, char *format, ...);
FILE    *OpenBinaryFileR (char *name);
FILE 	*OpenTextFileA (char *name);
FILE    *OpenTextFileR (char *name);
FILE    *OpenTextFileRQuait (char *name);
FILE 	*OpenTextFileW (char *name);
MrBFlt   PotentialScaleReduction (MrBFlt **vals, int nRows, int *count);
void     EstimatedSampleSize (MrBFlt **vals, int nRuns, int *count, MrBFlt *returnESS);
void    *SafeCalloc(size_t n, size_t s);
int      SafeFclose(FILE **fp);
void     SafeFree(void **ptr);
void    *SafeMalloc(size_t s);
void    *SafeRealloc(void *ptr, size_t s);
char    *SafeStrcat(char **target, const char *source);
char    *SafeStrcpy (char **target, const char *source);
void	 SetBit (int i, SafeLong *bits);
void     SortInts(int *item, int *assoc, int count, int descendingOrder);
void     SortInts2(int *item, int *assoc, int left, int right, int descendingOrder);
void     SortMrBFlt (MrBFlt *item, int left, int right);
int      StrCmpCaseInsensitive (char *s, char *t);
void     StripComments (char *s);
FILE    *TestOpenTextFileR (char *name);
void     UpdateGrowthFxn(int *growthFxn);
int      UpperTriangIndex(int i, int j, int size);
int      WantTo (const char *msg);
