#ifndef __MCMC_H__
#define __MCMC_H__

int     AddToPrintString (char *tempStr);
void    AutotuneDirichlet (MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *alphaPi, MrBFlt minTuning, MrBFlt maxTuning);
void    AutotuneMultiplier (MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *lambda, MrBFlt minTuning, MrBFlt maxTuning);
void    AutotuneSlider (MrBFlt acceptanceRate, MrBFlt targetRate, int batch, MrBFlt *width, MrBFlt minTuning, MrBFlt maxTuning);
int     DoMcmc (void);
int     DoMcmcp (void);
int     DoMcmcParm (char *parmName, char *tkn);
int     DoSs (void);
int     DoSsp (void);
int     DoSsParm (char *parmName, char *tkn);
int     ExhaustiveParsimonySearch (Tree *t, int chain, TreeInfo *tInfo);
MrBFlt  GetParsDP (Tree *t, TreeNode *p, int chain);
void    GetParsFP (Tree *t, TreeNode *p, int chain);
int     GetParsimonyBrlens (Tree *t, int chain, MrBFlt *brlens);
MrBFlt  GetParsimonyLength (Tree *t, int chain);
void    GetParsimonySubtreeRootstate (Tree *t, TreeNode *root, int chain);
MrBFlt  GetRate (int division, int chain);
int     LnBirthDeathPriorPr (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt sR, MrBFlt eR, char *sS, MrBFlt sF);
int     LnCoalescencePriorPr (Tree *t, MrBFlt *prob, MrBFlt theta, MrBFlt growth);
MrBFlt  LnUniformPriorPr (Tree *t, MrBFlt clockRate);
int     LnFossilizationPriorPr (Tree *t, MrBFlt clockRate, MrBFlt *prob, MrBFlt *sR, MrBFlt *eR, MrBFlt sF, MrBFlt *fR, char *sS);
int     LogClockTreePriorRatio (Param *param, int chain, MrBFlt *lnPriorRatio);
MrBFlt  LogDirPrior (Tree *t, ModelParams *mp, int PV);
MrBFlt  LogOmegaPrior (MrBFlt w1, MrBFlt w2, MrBFlt w3);
FILE    *OpenNewMBPrintFile (char *fileName);
int     ResetScalersPartition (int *isScalerNode, Tree* t, unsigned rescaleFreq);
int     SafeSprintf (char **target, int *targetLen, char *fmt, ...);
int     SetFilePositions (int samplePos);
MrBFlt  TreeLength (Param *param, int chain);

#endif  /* __MCMC_H__ */
