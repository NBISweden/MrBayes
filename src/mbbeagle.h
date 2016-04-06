#ifndef __MBBEAGLE_H__
#define __MBBEAGLE_H__

void   BeaglePrintResources (void);
void   BeaglePrintFlags (long inFlags);
void   BeagleNotLinked (void);
void   BeagleThreadsNotLinked (void);
int    BeagleCheckFlagCompatability (long inFlags);
void   BeagleAddGPUDevicesToList (int **beagleResource, int *beagleResourceCount);
void   BeagleRemoveGPUDevicesFromList (int **beagleResource, int *beagleResourceCount);
int    ScheduleLogLikeForAllDivisions (void);

void   LaunchBEAGLELogLikeForDivision (int chain, int d, ModelInfo* m, Tree* tree, MrBFlt* lnL);
void  *LaunchThreadLogLikeForDivision (void *arguments);
MrBFlt LaunchLogLikeForAllDivisionsInParallel (int chain);
void   recalculateScalers (int chain);

int    InitBeagleInstance (ModelInfo *m, int division);

int    TreeCondLikes_Beagle (Tree *t, int division, int chain);
int    TreeLikelihood_Beagle (Tree *t, int division, int chain, MrBFlt *lnL, int whichSitePats);
int    TreeTiProbs_Beagle (Tree *t, int division, int chain);

#endif  /* __MBBEAGLE_H__ */
