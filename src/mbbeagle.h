void BeaglePrintResources(void);
void BeaglePrintFlags(long inFlags);
void BeagleNotLinked(void);
void BeagleThreadsNotLinked(void);
int BeagleCheckFlagCompatability(long inFlags);
void BeagleAddGPUDevicesToList(int **beagleResource, int *beagleResourceCount);
void BeagleRemoveGPUDevicesFromList(int **beagleResource, int *beagleResourceCount);
int ScheduleLogLikeForAllDivisions(void);

void LaunchBEAGLELogLikeForDivision(int chain, int d, ModelInfo* m, Tree* tree, MrBFlt* lnL);
void *LaunchThreadLogLikeForDivision(void *arguments);
MrBFlt LaunchLogLikeForAllDivisionsInParallel(int chain);
