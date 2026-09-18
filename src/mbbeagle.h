#ifndef MBBEAGLE_H_
#define MBBEAGLE_H_

void   BeaglePrintResources (void);
void   BeaglePrintFlags (long inFlags);
void   BeagleNotLinked (void);
void   BeagleThreadsNotAvailable (void);
int    BeagleCheckFlagCompatability (long inFlags);
void   BeagleAddGPUDevicesToList (int **beagleResource, int *beagleResourceCount);
void   BeagleRemoveGPUDevicesFromList (int **beagleResource, int *beagleResourceCount);
int    ScheduleLogLikeForAllDivisions (void);

void   LaunchBEAGLELogLikeForDivision (int chain, int d, ModelInfo* m, Tree* tree, MrBFlt* lnL);
void  *LaunchThreadLogLikeForDivision (void *arguments);
MrBFlt LaunchLogLikeForAllDivisionsInParallel (int chain);
void   recalculateScalers (int chain);

int    InitBeagleInstance (ModelInfo *m, int division);

int    createBeagleInstance(ModelInfo *m, int nCijkParts, int numGammaCats, int numModelStates, int numCondLikes, int numScalers, int numChars, int numTiProbs, int numPartAmbigTips, int division);

int    TreeCondLikes_Beagle_Always_Rescale (Tree *t, int division, int chain);
int    TreeLikelihood_Beagle (Tree *t, int division, int chain, MrBFlt *lnL, int whichSitePats);
int    TreeTiProbs_Beagle (Tree *t, int division, int chain);
int    TreeCondLikes_Beagle_No_Rescale (Tree *t, int division, int chain);
int    TreeCondLikes_Beagle_Rescale_All (Tree *t, int division, int chain);

int    InitBeagleMultiPartitionInstance (void);
void   LaunchBEAGLELogLikeMultiPartition(int* divisions, int divisionCount, int chain, MrBFlt* lnL);
int    TreeTiProbs_BeagleMultiPartition (int* divisions, int divisionCount, int chain);
int    TreeCondLikes_BeagleMultiPartition_No_Rescale (int* divisions, int divisionCount, int chain);
int    TreeCondLikes_BeagleMultiPartition_Rescale_All (int* divisions, int divisionCount, int chain);
int    TreeCondLikes_BeagleMultiPartition_Always_Rescale (int* divisions, int divisionCount, int chain);
int    TreeLikelihood_BeagleMultiPartition (int* divisions, int divisionCount, int chain, MrBFlt *lnL, int whichSitePats);


//extern char *beagleGetVersion (void);

#endif  /* MBBEAGLE_H_ */
