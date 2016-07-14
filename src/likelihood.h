#ifndef __LIKELIHOOD_H__
#define __LIKELIHOOD_H__

//#define TIMING_ANALIZ
#if defined (TIMING_ANALIZ)
    static clock_t         CPUCondLikeDown;
    static clock_t         CPUScalers;
    static clock_t         CPUScalersRemove;
    static clock_t         CPUCondLikeRoot;
    static clock_t         CPULilklihood;

    #define TIME(X1,CPUtime)\
        {CPUTimeStart = clock();\
        X1;\
        CPUtime += (clock()-CPUTimeStart);}
#else
    #define TIME(X1,CPUtime)\
        X1;
#endif

#define A                           0
#define C                           1
#define G                           2
#define T                           3
#define AA                          0
#define AC                          1
#define AG                          2
#define AT                          3
#define CA                          4
#define CC                          5
#define CG                          6
#define CT                          7
#define GA                          8
#define GC                          9
#define GG                          10
#define GT                          11
#define TA                          12
#define TC                          13
#define TG                          14
#define TT                          15

int       CondLikeDown_Bin (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       CondLikeDown_Bin_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeDown_Gen (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       CondLikeDown_Gen_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeDown_Gen_GibbsGamma (TreeNode *p, int division, int chain);
int       CondLikeDown_NUC4 (TreeNode *p, int division, int chain);
#if defined (FMA_ENABLED)
int       CondLikeDown_NUC4_FMA (TreeNode *p, int division, int chain);
#endif
#if defined (AVX_ENABLED)
int       CondLikeDown_NUC4_AVX (TreeNode *p, int division, int chain);
#endif
#if defined (SSE_ENABLED)
int       CondLikeDown_NUC4_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeDown_NUC4_GibbsGamma (TreeNode *p, int division, int chain);
int       CondLikeDown_NY98 (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       CondLikeDown_NY98_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeDown_Std (TreeNode *p, int division, int chain);
int       CondLikeRoot_Bin (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       CondLikeRoot_Bin_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeRoot_Gen (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       CondLikeRoot_Gen_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeRoot_Gen_GibbsGamma (TreeNode *p, int division, int chain);
int       CondLikeRoot_NUC4 (TreeNode *p, int division, int chain);
#if defined (FMA_ENABLED)
int       CondLikeRoot_NUC4_FMA (TreeNode *p, int division, int chain);
#endif
#if defined (AVX_ENABLED)
int       CondLikeRoot_NUC4_AVX (TreeNode *p, int division, int chain);
#endif
#if defined (SSE_ENABLED)
int       CondLikeRoot_NUC4_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeRoot_NUC4_GibbsGamma (TreeNode *p, int division, int chain);
int       CondLikeRoot_NY98 (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       CondLikeRoot_NY98_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeRoot_Std (TreeNode *p, int division, int chain);
int       CondLikeScaler_Gen (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       CondLikeScaler_Gen_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeScaler_Gen_GibbsGamma (TreeNode *p, int division, int chain);
int       CondLikeScaler_NUC4 (TreeNode *p, int division, int chain);
#if defined (AVX_ENABLED)
int       CondLikeScaler_NUC4_AVX (TreeNode *p, int division, int chain);
#endif
#if defined (SSE_ENABLED)
int       CondLikeScaler_NUC4_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeScaler_NUC4_GibbsGamma (TreeNode *p, int division, int chain);
int       CondLikeScaler_NY98 (TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       CondLikeScaler_NY98_SSE (TreeNode *p, int division, int chain);
#endif
int       CondLikeScaler_Std (TreeNode *p, int division, int chain);
int       CondLikeUp_Bin (TreeNode *p, int division, int chain);
int       CondLikeUp_Gen (TreeNode *p, int division, int chain);
int       CondLikeUp_NUC4 (TreeNode *p, int division, int chain);
int       CondLikeUp_Std (TreeNode *p, int division, int chain);
void      LaunchLogLikeForDivision (int chain, int d, MrBFlt* lnL);
int       Likelihood_Adgamma (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
int       Likelihood_Gen (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#if defined (SSE_ENABLED)
int       Likelihood_Gen_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#endif
int       Likelihood_Gen_GibbsGamma (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
int       Likelihood_NUC4 (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#if defined (FMA_ENABLED)
int       Likelihood_NUC4_FMA (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#endif
#if defined (AVX_ENABLED)
int       Likelihood_NUC4_AVX (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#endif
#if defined (SSE_ENABLED)
int       Likelihood_NUC4_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#endif
int       Likelihood_NUC4_GibbsGamma (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
int       Likelihood_NY98 (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#if defined (SSE_ENABLED)
int       Likelihood_NY98_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#endif
int       Likelihood_Pars (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
int       Likelihood_ParsStd (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
int       Likelihood_Res (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#if defined (SSE_ENABLED)
int       Likelihood_Res_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
#endif
int       Likelihood_Std (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats);
int       TiProbs_Fels (TreeNode *p, int division, int chain);
int       TiProbs_Gen (TreeNode *p, int division, int chain);
int       TiProbs_GenCov (TreeNode *p, int division, int chain);
int       TiProbs_Hky (TreeNode *p, int division, int chain);
int       TiProbs_JukesCantor (TreeNode *p, int division, int chain);
int       TiProbs_Std (TreeNode *p, int division, int chain);
int       TiProbs_Res (TreeNode *p, int division, int chain);

#endif  /* __LIKELIHOOD_H__ */
