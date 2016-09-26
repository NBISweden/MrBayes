/*
 *  MrBayes 3
 *
 *  (c) 2002-2013
 *
 *  John P. Huelsenbeck
 *  Dept. Integrative Biology
 *  University of California, Berkeley
 *  Berkeley, CA 94720-3140
 *  johnh@berkeley.edu
 *
 *  Fredrik Ronquist
 *  Swedish Museum of Natural History
 *  Box 50007
 *  SE-10405 Stockholm, SWEDEN
 *  fredrik.ronquist@nrm.se
 *
 *  With important contributions by
 *
 *  Paul van der Mark (paulvdm@sc.fsu.edu)
 *  Maxim Teslenko (maxkth@gmail.com)
 *  Chi Zhang (zhangchicool@gmail.com)
 *
 *  and by many users (run 'acknowledgments' to see more info)
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

#include "bayes.h"
#include "likelihood.h"
#include "mbbeagle.h"
#include "model.h"
#include "utils.h"

#define LIKE_EPSILON                1.0e-300

/* global variables declared here */
CLFlt     *preLikeL;                  /* precalculated cond likes for left descendant */
CLFlt     *preLikeR;                  /* precalculated cond likes for right descendant*/
CLFlt     *preLikeA;                  /* precalculated cond likes for ancestor        */

/* global variables used here but declared elsewhere */
extern int      *chainId;
extern int      numLocalChains;
extern int      rateProbRowSize;            /* size of rate probs for one chain one state   */
extern MrBFlt   **rateProbs;                /* pointers to rate probs used by adgamma model */

/* local prototypes */
void      CopySiteScalers (ModelInfo *m, int chain);
void      FlipCondLikeSpace (ModelInfo *m, int chain, int nodeIndex);
void      FlipCijkSpace (ModelInfo *m, int chain);
void      FlipNodeScalerSpace (ModelInfo *m, int chain, int nodeIndex);
void      FlipSiteScalerSpace (ModelInfo *m, int chain);
void      FlipTiProbsSpace (ModelInfo *m, int chain, int nodeIndex);
MrBFlt    GetRate (int division, int chain);
int       RemoveNodeScalers(TreeNode *p, int division, int chain);
#if defined (SSE_ENABLED)
int       RemoveNodeScalers_SSE(TreeNode *p, int division, int chain);
#endif
#if defined (AVX_ENABLED)
int       RemoveNodeScalers_AVX(TreeNode *p, int division, int chain);
#endif
void      ResetSiteScalers (ModelInfo *m, int chain);
int       SetBinaryQMatrix (MrBFlt **a, int whichChain, int division);
int       SetNucQMatrix (MrBFlt **a, int n, int whichChain, int division, MrBFlt rateMult, MrBFlt *rA, MrBFlt *rS);
int       SetStdQMatrix (MrBFlt **a, int nStates, MrBFlt *bs, int cType);
int       SetProteinQMatrix (MrBFlt **a, int n, int whichChain, int division, MrBFlt rateMult);
int       UpDateCijk (int whichPart, int whichChain);


#if !defined (SSE_ENABLED) || 1
/*----------------------------------------------------------------
|
|   CondLikeDown_Bin: binary model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeDown_Bin (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *clL, *clR, *clP, *pL, *pR, *tiPL, *tiPR;
    ModelInfo       *m;
    
    /* find model settings for this division */
    m = &modelSettings[division];

    /* Flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    tiPL = pL;
    tiPR = pR;
    for (k=0; k<m->numRateCats; k++)
        {
        for (c=0; c<m->numChars; c++)
            {
            *(clP++) = (tiPL[0]*clL[0] + tiPL[1]*clL[1])
                      *(tiPR[0]*clR[0] + tiPR[1]*clR[1]);
            *(clP++) = (tiPL[2]*clL[0] + tiPL[3]*clL[1])
                      *(tiPR[2]*clR[0] + tiPR[3]*clR[1]);

            clL += 2;
            clR += 2;
            }
        tiPL += 4;
        tiPR += 4;
        }

    return NO_ERROR;
    
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeDown_Bin_SSE: binary model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeDown_Bin_SSE (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *pL, *pR, *tiPL, *tiPR;
    __m128          *clL, *clR, *clP;
    __m128          m1, m2, m3, m4, m5, m6;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    tiPL = pL;
    tiPR = pR;
    for (k=0; k<m->numRateCats; k++)
        {
        for (c=0; c<m->numVecChars; c++)
            {
            m1 = _mm_load1_ps (&tiPL[0]);
            m2 = _mm_load1_ps (&tiPR[0]);
            m5 = _mm_mul_ps (m1, clL[0]);
            m6 = _mm_mul_ps (m2, clR[0]);

            m1 = _mm_load1_ps (&tiPL[1]);
            m2 = _mm_load1_ps (&tiPR[1]);
            m3 = _mm_mul_ps (m1, clL[1]);
            m4 = _mm_mul_ps (m2, clR[1]);

            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            *clP++ = _mm_mul_ps (m5, m6);

            m1 = _mm_load1_ps (&tiPL[2]);
            m2 = _mm_load1_ps (&tiPR[2]);
            m5 = _mm_mul_ps (m1, clL[0]);
            m6 = _mm_mul_ps (m2, clR[0]);

            m1 = _mm_load1_ps (&tiPL[3]);
            m2 = _mm_load1_ps (&tiPR[3]);
            m3 = _mm_mul_ps (m1, clL[1]);
            m4 = _mm_mul_ps (m2, clR[1]);

            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);
           
            *clP++ = _mm_mul_ps (m5, m6);
            clL += 2;
            clR += 2;
            }
        tiPL += 4;
        tiPR += 4;
        }

    return NO_ERROR;
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeDown_Gen: general n-state model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeDown_Gen (TreeNode *p, int division, int chain)
{
    int             a, b, c, h, i, k, j, shortCut, *lState=NULL, *rState=NULL,
                    nObsStates, nStates, nStatesSquared, preLikeJump;
    CLFlt           likeL, likeR, *pL, *pR, *tiPL, *tiPR, *clL, *clR, *clP;
    ModelInfo       *m;
#   if !defined (DEBUG_NOSHORTCUTS)
    int catStart;
#   endif
    
    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nObsStates = m->numStates;
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;
    preLikeJump = nObsStates * nStates;

    /* flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeL[a++] += tiPL[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeR[a++] += tiPR[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }
#   endif
    switch (shortCut)
        {
        case 0:
            tiPL = pL;
            tiPR = pR;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    for (i=h=0; i<nStates; i++)
                        {
                        likeL = likeR = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += tiPL[h]*clL[j];
                            likeR += tiPR[h++]*clR[j];
                            }
                        *(clP++) = likeL * likeR;
                        }
                    clL += nStates;
                    clR += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = lState[c] + k*(preLikeJump+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeR = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeR += tiPR[h++]*clR[j];
                            }
                        *(clP++) = preLikeL[a++] * likeR;
                        }
                    clR += nStates;
                    }
                tiPR += nStatesSquared;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = rState[c] + k*(preLikeJump+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeL = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += tiPL[h++]*clL[j];
                            }
                        *(clP++) = preLikeR[a++] * likeL;
                        }
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                }
            break;
        case 3:
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = rState[c] + k*(preLikeJump+nStates);
                    b = lState[c] + k*(preLikeJump+nStates);
                    for (i=0; i<nStates; i++)
                        {
                        *(clP++) = preLikeR[a++] * preLikeL[b++];
                        }
                    }
                }
            break;
        }

    return NO_ERROR;
}


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeDown_Gen_SSE: general n-state model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeDown_Gen_SSE (TreeNode *p, int division, int chain)
{
    int             c, c1, h, i, j, k, t, shortCut, *lState=NULL, *rState=NULL, nStates, nStatesSquared, nObsStates, preLikeJump;
    CLFlt           *pL, *pR, *tiPL, *tiPR;
    __m128          *clL, *clR, *clP;
    __m128          mTiPL, mTiPR, mL, mR, mAcumL, mAcumR;
    ModelInfo       *m;
    CLFlt           *preLikeRV[4];
    CLFlt           *preLikeLV[4];

#   if !defined (DEBUG_NOSHORTCUTS)
    int             a, b, catStart;
#   endif
    
    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nObsStates = m->numStates;
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;
    preLikeJump = nObsStates * nStates;

    /* Flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeL[a++] += tiPL[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeR[a++] += tiPR[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }
#   endif

    switch (shortCut)
        {
        case 0:
            tiPL = pL;
            tiPR = pR;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numVecChars; c++)
                    {
                    for (i=h=0; i<nStates; i++)
                        {
                        mAcumL = _mm_setzero_ps();
                        mAcumR = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h]);
                            mTiPR  = _mm_load1_ps (&tiPR[h++]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            }
                        *(clP++) = _mm_mul_ps (mAcumL, mAcumR);
                        }
                    clL += nStates;
                    clR += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeLV[c1] = &preLikeL[lState[t] + k*(preLikeJump+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        mAcumL = _mm_set_ps (*(preLikeLV[3]++), *(preLikeLV[2]++), *(preLikeLV[1]++), *(preLikeLV[0]++));
                        mAcumR = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPR  = _mm_load1_ps (&tiPR[h++]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            }
                        *(clP++) = _mm_mul_ps (mAcumL,mAcumR);
                        }
                    clR += nStates;
                    }
                tiPR += nStatesSquared;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeRV[c1] = &preLikeR[rState[t] + k*(preLikeJump+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        mAcumR = _mm_set_ps (*(preLikeRV[3]++), *(preLikeRV[2]++), *(preLikeRV[1]++), *(preLikeRV[0]++));
                        mAcumL = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h++]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            }
                        *(clP++) = _mm_mul_ps (mAcumL,mAcumR);
                        }
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                }
            break;
        case 3:
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeRV[c1] = &preLikeR[rState[t] + k*(preLikeJump+nStates)];
                        preLikeLV[c1] = &preLikeL[lState[t] + k*(preLikeJump+nStates)];
                        }
                    for (i=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following 2 statments we assume that SSE register can hold exactly 4 ClFlts. */
                        mL = _mm_set_ps (*(preLikeLV[3]++), *(preLikeLV[2]++), *(preLikeLV[1]++), *(preLikeLV[0]++));
                        mR = _mm_set_ps (*(preLikeRV[3]++), *(preLikeRV[2]++), *(preLikeRV[1]++), *(preLikeRV[0]++));
                        *(clP++) = _mm_mul_ps (mL,mR);
                        }
                    }
                }
            break;
        }
    return NO_ERROR;
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeDown_Gen_GibbsGamma: general n-state model with rate
|       variation modeled using discrete gamma with Gibbs resampling
|
-----------------------------------------------------------------*/
int CondLikeDown_Gen_GibbsGamma (TreeNode *p, int division, int chain)
{
    int             a, b, c, i, j, r, *rateCat, shortCut, *lState=NULL, *rState=NULL,
                    nObsStates, nStates, nStatesSquared, nGammaCats;
    CLFlt           likeL, likeR, *pL, *pR, *tiPL, *tiPR, *clL, *clR, *clP;
    ModelInfo       *m;
#   if !defined (DEBUG_NOSHORTCUTS)
    int k, catStart;
#   endif
    
    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nObsStates = m->numStates;
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;

    /* flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* find rate category index and number of gamma categories */
    rateCat = m->tiIndex + chain * m->numChars;
    nGammaCats = m->numRateCats;

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<nGammaCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeL[a++] += tiPL[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<nGammaCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeR[a++] += tiPR[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }
#   endif

    switch (shortCut)
        {
        case 0:
            for (c=0; c<m->numChars; c++)
                {
                r = (*rateCat++);
                if (r < nGammaCats)
                    {
                    tiPL = pL + r*nStatesSquared;
                    tiPR = pR + r*nStatesSquared;
                    for (i=0; i<nStates; i++)
                        {
                        likeL = likeR = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += (*tiPL++) * clL[j];
                            likeR += (*tiPR++) * clR[j];
                            }
                        *(clP++) = likeL * likeR;
                        }
                    }
                else
                    clP += nStates;
                clL += nStates;
                clR += nStates;
                }
            break;
        case 1:
            for (c=0; c<m->numChars; c++)
                {
                r = (*rateCat++);
                if (r < nGammaCats)
                    {
                    tiPR = pR + r*nStatesSquared;
                    a = lState[c] + r*(nStatesSquared+nStates);
                    for (i=0; i<nStates; i++)
                        {
                        likeR = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeR += (*tiPR++)*clR[j];
                            }
                        *(clP++) = preLikeL[a++] * likeR;
                        }
                    }
                else
                    clP += nStates;
                clR += nStates;
                }
            break;
        case 2:
            for (c=0; c<m->numChars; c++)
                {
                r = (*rateCat++);
                if (r < nGammaCats)
                    {
                    tiPL = pL + r*nStatesSquared;
                    a = rState[c] + r*(nStatesSquared+nStates);
                    for (i=0; i<nStates; i++)
                        {
                        likeL = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += (*tiPL++)*clL[j];
                            }
                        *(clP++) = preLikeR[a++] * likeL;
                        }
                    }
                else
                    clP += nStates;
                clL += nStates;
                }
            break;
        case 3:
            for (c=0; c<m->numChars; c++)
                {
                r = (*rateCat++);
                if (r < nGammaCats)
                    {
                    a = lState[c] + r*(nStatesSquared+nStates);
                    b = rState[c] + r*(nStatesSquared+nStates);
                    for (i=0; i<nStates; i++)
                        *(clP++) = preLikeL[a++]*preLikeR[b++];
                    }
                else
                    clP += nStates;
                }
            break;
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeDown_NUC4: 4by4 nucleotide model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeDown_NUC4 (TreeNode *p, int division, int chain)
{
    int             c, h, i, j, k, shortCut, *lState=NULL, *rState=NULL;
    CLFlt           *clL, *clR, *clP, *pL, *pR, *tiPL, *tiPR;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* flip space so that we do not overwrite old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=j=0; k<m->numRateCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeL[j++] = tiPL[0];
                preLikeL[j++] = tiPL[4];
                preLikeL[j++] = tiPL[8];
                preLikeL[j++] = tiPL[12];
                tiPL++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeL[j++] = 1.0;
            tiPL += 12;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=j=0; k<m->numRateCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeR[j++] = tiPR[0];
                preLikeR[j++] = tiPR[4];
                preLikeR[j++] = tiPR[8];
                preLikeR[j++] = tiPR[12];
                tiPR++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeR[j++] = 1.0;
            tiPR += 12;
            }
        }
#   endif

    switch (shortCut)
        {
        case 0:
            tiPL = pL;
            tiPR = pR;
            for (k=h=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                                *(tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T]);
                    clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                                *(tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T]);
                    clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                                *(tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T]);
                    clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                                *(tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T]);
                    clL += 4;
                    clR += 4;
                    }
                tiPL += 16;
                tiPR += 16;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=h=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    i = lState[c] + k*20;
                    clP[h++] =   preLikeL[i++]
                                *(tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T]);
                    clP[h++] =   preLikeL[i++]
                                *(tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T]);
                    clP[h++] =   preLikeL[i++]
                                *(tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T]);
                    clP[h++] =   preLikeL[i++]
                                *(tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T]);
                    clR += 4;
                    }
                tiPR += 16;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=h=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    i = rState[c] + k*20;
                    clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                                *preLikeR[i++];
                    clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                                *preLikeR[i++];
                    clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                                *preLikeR[i++];
                    clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                                *preLikeR[i++];
                    clL += 4;
                    }
                tiPL += 16;
                }
            break;
        case 3:
            for (k=h=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    i = j = k*20;
                    i += lState[c];
                    j += rState[c];
                    clP[h++] =   preLikeL[i++]*preLikeR[j++];
                    clP[h++] =   preLikeL[i++]*preLikeR[j++];
                    clP[h++] =   preLikeL[i++]*preLikeR[j++];
                    clP[h++] =   preLikeL[i++]*preLikeR[j++];
                    }
                }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeDown_NUC4_GibbsGamma: 4by4 nucleotide model with rate
|       variation approximated using Gibbs sampling of gamma
|
-----------------------------------------------------------------*/
int CondLikeDown_NUC4_GibbsGamma (TreeNode *p, int division, int chain)
{
    int             c, h, i, j, r, *rateCat, shortCut, *lState=NULL, *rState=NULL,
                    nGammaCats;
    CLFlt           *clL, *clR, *clP, *pL, *pR, *tiPL, *tiPR;
    ModelInfo       *m;
#   if !defined (DEBUG_NOSHORTCUTS)
    int k;
#   endif
    
    m = &modelSettings[division];

    /* flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* find rate category index  and number of gamma categories */
    rateCat = m->tiIndex + chain * m->numChars;
    nGammaCats = m->numRateCats;

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=j=0; k<nGammaCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeL[j++] = tiPL[0];
                preLikeL[j++] = tiPL[4];
                preLikeL[j++] = tiPL[8];
                preLikeL[j++] = tiPL[12];
                tiPL++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeL[j++] = 1.0;
            tiPL += 12;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState =  m->termState[p->right->index];
        tiPR = pR;
        for (k=j=0; k<nGammaCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeR[j++] = tiPR[0];
                preLikeR[j++] = tiPR[4];
                preLikeR[j++] = tiPR[8];
                preLikeR[j++] = tiPR[12];
                tiPR++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeR[j++] = 1.0;
            tiPR += 12;
            }
        }
#   endif

    switch (shortCut)
        {
        case 0:
            for (c=h=0; c<m->numChars; c++)
                {
                r = rateCat[c];
                if (r < nGammaCats)
                    {
                    tiPL = pL + r * 16;
                    tiPR = pR + r * 16;
                    clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                                *(tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T]);
                    clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                                *(tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T]);
                    clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                                *(tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T]);
                    clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                                *(tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T]);
                    }
                else
                    h += 4;
                clL += 4;
                clR += 4;
                }
            break;
        case 1:
            for (c=h=0; c<m->numChars; c++)
                {
                r = rateCat[c];
                if (r < nGammaCats)
                    {
                    tiPR = pR + r * 16;
                    i = lState[c] + r * 20;
                    clP[h++] =   preLikeL[i++]
                                *(tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T]);
                    clP[h++] =   preLikeL[i++]
                                *(tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T]);
                    clP[h++] =   preLikeL[i++]
                                *(tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T]);
                    clP[h++] =   preLikeL[i++]
                                *(tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T]);
                    }
                else
                    h += 4;
                clR += 4;
                }
            break;
        case 2:
            for (c=h=0; c<m->numChars; c++)
                {
                r = rateCat[c];
                if (r < nGammaCats)
                    {
                    tiPL = pL + r * 16;
                    i = rState[c] + r * 20;
                    clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                                *preLikeR[i++];
                    clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                                *preLikeR[i++];
                    clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                                *preLikeR[i++];
                    clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                                *preLikeR[i++];
                    }
                else
                    h += 4;
                clL += 4;
                }
            break;
        case 3:
            for (c=h=0; c<m->numChars; c++)
                {
                r = rateCat[c];
                if (r < nGammaCats)
                    {
                    i = lState[c] + r * 20;
                    j = rState[c] + r * 20;
                    clP[h++] =   preLikeL[i++]*preLikeR[j++];
                    clP[h++] =   preLikeL[i++]*preLikeR[j++];
                    clP[h++] =   preLikeL[i++]*preLikeR[j++];
                    clP[h++] =   preLikeL[i++]*preLikeR[j++];
                    }
                else
                    h += 4;
                }
            break;
        }

    return NO_ERROR;
}


#if defined (FMA_ENABLED)
/*----------------------------------------------------------------
 |
 |   CondLikeDown_NUC4_FMA: 4by4 nucleotide model with or without rate
 |       variation, using AVX + FMA instructions
 |
 -----------------------------------------------------------------*/
int CondLikeDown_NUC4_FMA (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *pL, *pR, *tiPL, *tiPR;
    __m256          *clL, *clR, *clP;
    __m256          m1, m2, m3, m4;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    
    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    
    tiPL = pL;
    tiPR = pR;
    for (k=0; k<m->numRateCats; k++)
    {
        for (c=0; c<m->numVecChars; c++)
        {
            m1 = _mm256_broadcast_ss (&tiPL[AA]);
            m2 = _mm256_broadcast_ss (&tiPR[AA]);
            m3 = _mm256_mul_ps (m1, clL[A]);
            m4 = _mm256_mul_ps (m2, clR[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[AC]);
            m2 = _mm256_broadcast_ss (&tiPR[AC]);
            m3 = _mm256_fmadd_ps (m1, clL[C], m3);
            m4 = _mm256_fmadd_ps (m2, clR[C], m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[AG]);
            m2 = _mm256_broadcast_ss (&tiPR[AG]);
            m3 = _mm256_fmadd_ps (m1, clL[G], m3);
            m4 = _mm256_fmadd_ps (m2, clR[G], m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[AT]);
            m2 = _mm256_broadcast_ss (&tiPR[AT]);
            m3 = _mm256_fmadd_ps (m1, clL[T], m3);
            m4 = _mm256_fmadd_ps (m2, clR[T], m4);
            
            *clP++ = _mm256_mul_ps (m3, m4);

            m1 = _mm256_broadcast_ss (&tiPL[CA]);
            m2 = _mm256_broadcast_ss (&tiPR[CA]);
            m3 = _mm256_mul_ps (m1, clL[A]);
            m4 = _mm256_mul_ps (m2, clR[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[CC]);
            m2 = _mm256_broadcast_ss (&tiPR[CC]);
            m3 = _mm256_fmadd_ps (m1, clL[C], m3);
            m4 = _mm256_fmadd_ps (m2, clR[C], m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[CG]);
            m2 = _mm256_broadcast_ss (&tiPR[CG]);
            m3 = _mm256_fmadd_ps (m1, clL[G], m3);
            m4 = _mm256_fmadd_ps (m2, clR[G], m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[CT]);
            m2 = _mm256_broadcast_ss (&tiPR[CT]);
            m3 = _mm256_fmadd_ps (m1, clL[T], m3);
            m4 = _mm256_fmadd_ps (m2, clR[T], m4);
            
            *clP++ = _mm256_mul_ps (m3, m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[GA]);
            m2 = _mm256_broadcast_ss (&tiPR[GA]);
            m3 = _mm256_mul_ps (m1, clL[A]);
            m4 = _mm256_mul_ps (m2, clR[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[GC]);
            m2 = _mm256_broadcast_ss (&tiPR[GC]);
            m3 = _mm256_fmadd_ps (m1, clL[C], m3);
            m4 = _mm256_fmadd_ps (m2, clR[C], m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[GG]);
            m2 = _mm256_broadcast_ss (&tiPR[GG]);
            m3 = _mm256_fmadd_ps (m1, clL[G], m3);
            m4 = _mm256_fmadd_ps (m2, clR[G], m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[GT]);
            m2 = _mm256_broadcast_ss (&tiPR[GT]);
            m3 = _mm256_fmadd_ps (m1, clL[T], m3);
            m4 = _mm256_fmadd_ps (m2, clR[T], m4);
            
            *clP++ = _mm256_mul_ps (m3, m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[TA]);
            m2 = _mm256_broadcast_ss (&tiPR[TA]);
            m3 = _mm256_mul_ps (m1, clL[A]);
            m4 = _mm256_mul_ps (m2, clR[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[TC]);
            m2 = _mm256_broadcast_ss (&tiPR[TC]);
            m3 = _mm256_fmadd_ps (m1, clL[C], m3);
            m4 = _mm256_fmadd_ps (m2, clR[C], m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[TG]);
            m2 = _mm256_broadcast_ss (&tiPR[TG]);
            m3 = _mm256_fmadd_ps (m1, clL[G], m3);
            m4 = _mm256_fmadd_ps (m2, clR[G], m4);
            
            m1 = _mm256_broadcast_ss (&tiPL[TT]);
            m2 = _mm256_broadcast_ss (&tiPR[TT]);
            m3 = _mm256_fmadd_ps (m1, clL[T], m3);
            m4 = _mm256_fmadd_ps (m2, clR[T], m4);
            
            *clP++ = _mm256_mul_ps (m3, m4);
            
            clL += 4;
            clR += 4;
        }
        tiPL += 16;
        tiPR += 16;
    }
    
    return NO_ERROR;
    
}
#endif


#if defined (AVX_ENABLED)
/*----------------------------------------------------------------
 |
 |   CondLikeDown_NUC4_AVX: 4by4 nucleotide model with or without rate
 |       variation, using AVX instructions
 |
 -----------------------------------------------------------------*/
int CondLikeDown_NUC4_AVX (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *pL, *pR, *tiPL, *tiPR;
    __m256          *clL, *clR, *clP;
    __m256          m1, m2, m3, m4, m5, m6;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    
    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    
    tiPL = pL;
    tiPR = pR;
    for (k=0; k<m->numRateCats; k++)
    {
        for (c=0; c<m->numVecChars; c++)
        {
            m1 = _mm256_broadcast_ss (&tiPL[AA]);
            m2 = _mm256_broadcast_ss (&tiPR[AA]);
            m5 = _mm256_mul_ps (m1, clL[A]);
            m6 = _mm256_mul_ps (m2, clR[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[AC]);
            m2 = _mm256_broadcast_ss (&tiPR[AC]);
            m3 = _mm256_mul_ps (m1, clL[C]);
            m4 = _mm256_mul_ps (m2, clR[C]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[AG]);
            m2 = _mm256_broadcast_ss (&tiPR[AG]);
            m3 = _mm256_mul_ps (m1, clL[G]);
            m4 = _mm256_mul_ps (m2, clR[G]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[AT]);
            m2 = _mm256_broadcast_ss (&tiPR[AT]);
            m3 = _mm256_mul_ps (m1, clL[T]);
            m4 = _mm256_mul_ps (m2, clR[T]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            *clP++ = _mm256_mul_ps (m5, m6);

            m1 = _mm256_broadcast_ss (&tiPL[CA]);
            m2 = _mm256_broadcast_ss (&tiPR[CA]);
            m5 = _mm256_mul_ps (m1, clL[A]);
            m6 = _mm256_mul_ps (m2, clR[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[CC]);
            m2 = _mm256_broadcast_ss (&tiPR[CC]);
            m3 = _mm256_mul_ps (m1, clL[C]);
            m4 = _mm256_mul_ps (m2, clR[C]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[CG]);
            m2 = _mm256_broadcast_ss (&tiPR[CG]);
            m3 = _mm256_mul_ps (m1, clL[G]);
            m4 = _mm256_mul_ps (m2, clR[G]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[CT]);
            m2 = _mm256_broadcast_ss (&tiPR[CT]);
            m3 = _mm256_mul_ps (m1, clL[T]);
            m4 = _mm256_mul_ps (m2, clR[T]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            *clP++ = _mm256_mul_ps (m5, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[GA]);
            m2 = _mm256_broadcast_ss (&tiPR[GA]);
            m5 = _mm256_mul_ps (m1, clL[A]);
            m6 = _mm256_mul_ps (m2, clR[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[GC]);
            m2 = _mm256_broadcast_ss (&tiPR[GC]);
            m3 = _mm256_mul_ps (m1, clL[C]);
            m4 = _mm256_mul_ps (m2, clR[C]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[GG]);
            m2 = _mm256_broadcast_ss (&tiPR[GG]);
            m3 = _mm256_mul_ps (m1, clL[G]);
            m4 = _mm256_mul_ps (m2, clR[G]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[GT]);
            m2 = _mm256_broadcast_ss (&tiPR[GT]);
            m3 = _mm256_mul_ps (m1, clL[T]);
            m4 = _mm256_mul_ps (m2, clR[T]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            *clP++ = _mm256_mul_ps (m5, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[TA]);
            m2 = _mm256_broadcast_ss (&tiPR[TA]);
            m5 = _mm256_mul_ps (m1, clL[A]);
            m6 = _mm256_mul_ps (m2, clR[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[TC]);
            m2 = _mm256_broadcast_ss (&tiPR[TC]);
            m3 = _mm256_mul_ps (m1, clL[C]);
            m4 = _mm256_mul_ps (m2, clR[C]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[TG]);
            m2 = _mm256_broadcast_ss (&tiPR[TG]);
            m3 = _mm256_mul_ps (m1, clL[G]);
            m4 = _mm256_mul_ps (m2, clR[G]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[TT]);
            m2 = _mm256_broadcast_ss (&tiPR[TT]);
            m3 = _mm256_mul_ps (m1, clL[T]);
            m4 = _mm256_mul_ps (m2, clR[T]);
            m5 = _mm256_add_ps (m3, m5);
            m6 = _mm256_add_ps (m4, m6);
            
            *clP++ = _mm256_mul_ps (m5, m6);
            
            clL += 4;
            clR += 4;
        }
        tiPL += 16;
        tiPR += 16;
    }
    
    return NO_ERROR;
    
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeDown_NUC4_SSE: 4by4 nucleotide model with or without rate
|       variation, using SSE instructions
|
-----------------------------------------------------------------*/
int CondLikeDown_NUC4_SSE (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *pL, *pR, *tiPL, *tiPR;
    __m128          *clL, *clR, *clP;
    __m128          m1, m2, m3, m4, m5, m6;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    tiPL = pL;
    tiPR = pR;
    for (k=0; k<m->numRateCats; k++)
        {
        for (c=0; c<m->numVecChars; c++)
            {
            m1 = _mm_load1_ps (&tiPL[AA]);
            m2 = _mm_load1_ps (&tiPR[AA]);
            m5 = _mm_mul_ps (m1, clL[A]);
            m6 = _mm_mul_ps (m2, clR[A]);

            m1 = _mm_load1_ps (&tiPL[AC]);
            m2 = _mm_load1_ps (&tiPR[AC]);
            m3 = _mm_mul_ps (m1, clL[C]);
            m4 = _mm_mul_ps (m2, clR[C]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            m1 = _mm_load1_ps (&tiPL[AG]);
            m2 = _mm_load1_ps (&tiPR[AG]);
            m3 = _mm_mul_ps (m1, clL[G]);
            m4 = _mm_mul_ps (m2, clR[G]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            m1 = _mm_load1_ps (&tiPL[AT]);
            m2 = _mm_load1_ps (&tiPR[AT]);
            m3 = _mm_mul_ps (m1, clL[T]);
            m4 = _mm_mul_ps (m2, clR[T]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            *clP++ = _mm_mul_ps (m5, m6);

            m1 = _mm_load1_ps (&tiPL[CA]);
            m2 = _mm_load1_ps (&tiPR[CA]);
            m5 = _mm_mul_ps (m1, clL[A]);
            m6 = _mm_mul_ps (m2, clR[A]);

            m1 = _mm_load1_ps (&tiPL[CC]);
            m2 = _mm_load1_ps (&tiPR[CC]);
            m3 = _mm_mul_ps (m1, clL[C]);
            m4 = _mm_mul_ps (m2, clR[C]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            m1 = _mm_load1_ps (&tiPL[CG]);
            m2 = _mm_load1_ps (&tiPR[CG]);
            m3 = _mm_mul_ps (m1, clL[G]);
            m4 = _mm_mul_ps (m2, clR[G]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            m1 = _mm_load1_ps (&tiPL[CT]);
            m2 = _mm_load1_ps (&tiPR[CT]);
            m3 = _mm_mul_ps (m1, clL[T]);
            m4 = _mm_mul_ps (m2, clR[T]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            *clP++ = _mm_mul_ps (m5, m6);

            m1 = _mm_load1_ps (&tiPL[GA]);
            m2 = _mm_load1_ps (&tiPR[GA]);
            m5 = _mm_mul_ps (m1, clL[A]);
            m6 = _mm_mul_ps (m2, clR[A]);

            m1 = _mm_load1_ps (&tiPL[GC]);
            m2 = _mm_load1_ps (&tiPR[GC]);
            m3 = _mm_mul_ps (m1, clL[C]);
            m4 = _mm_mul_ps (m2, clR[C]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            m1 = _mm_load1_ps (&tiPL[GG]);
            m2 = _mm_load1_ps (&tiPR[GG]);
            m3 = _mm_mul_ps (m1, clL[G]);
            m4 = _mm_mul_ps (m2, clR[G]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            m1 = _mm_load1_ps (&tiPL[GT]);
            m2 = _mm_load1_ps (&tiPR[GT]);
            m3 = _mm_mul_ps (m1, clL[T]);
            m4 = _mm_mul_ps (m2, clR[T]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            *clP++ = _mm_mul_ps (m5, m6);

            m1 = _mm_load1_ps (&tiPL[TA]);
            m2 = _mm_load1_ps (&tiPR[TA]);
            m5 = _mm_mul_ps (m1, clL[A]);
            m6 = _mm_mul_ps (m2, clR[A]);

            m1 = _mm_load1_ps (&tiPL[TC]);
            m2 = _mm_load1_ps (&tiPR[TC]);
            m3 = _mm_mul_ps (m1, clL[C]);
            m4 = _mm_mul_ps (m2, clR[C]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            m1 = _mm_load1_ps (&tiPL[TG]);
            m2 = _mm_load1_ps (&tiPR[TG]);
            m3 = _mm_mul_ps (m1, clL[G]);
            m4 = _mm_mul_ps (m2, clR[G]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            m1 = _mm_load1_ps (&tiPL[TT]);
            m2 = _mm_load1_ps (&tiPR[TT]);
            m3 = _mm_mul_ps (m1, clL[T]);
            m4 = _mm_mul_ps (m2, clR[T]);
            m5 = _mm_add_ps (m3, m5);
            m6 = _mm_add_ps (m4, m6);

            *clP++ = _mm_mul_ps (m5, m6);

            clL += 4;
            clR += 4;
            }
        tiPL += 16;
        tiPR += 16;
        }

    return NO_ERROR;
    
}
#endif


#if !defined (SSE_ENABLED) || 1
/*----------------------------------------------------------------
|
|   CondLikeDown_NY98: codon model with omega variation
|
-----------------------------------------------------------------*/
int CondLikeDown_NY98 (TreeNode *p, int division, int chain)
{
    int             a, b, c, h, i, j, k, shortCut, *lState=NULL, *rState=NULL, nStates, nStatesSquared;
    CLFlt           likeL, likeR, *pL, *pR, *tiPL, *tiPR, *clL, *clR, *clP;
    ModelInfo       *m;
    
    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;

    /* Flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }
#   endif

    switch (shortCut)
        {
        case 0:
            tiPL = pL;
            tiPR = pR;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    for (i=h=0; i<nStates; i++)
                        {
                        likeL = likeR = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += tiPL[h]*clL[j];
                            likeR += tiPR[h++]*clR[j];
                            }
                        *(clP++) = likeL * likeR;
                        }
                    clL += nStates;
                    clR += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = lState[c] + k*(nStatesSquared+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeR = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeR += tiPR[h++]*clR[j];
                            }
                        *(clP++) = preLikeL[a++] * likeR;
                        }
                    clR += nStates;
                    }
                tiPR += nStatesSquared;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = rState[c] + k*(nStatesSquared+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeL = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += tiPL[h++]*clL[j];
                            }
                        *(clP++) = preLikeR[a++] * likeL;
                        }
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                }
            break;
        case 3:
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = rState[c] + k*(nStatesSquared+nStates);
                    b = lState[c] + k*(nStatesSquared+nStates);
                    for (i=0; i<nStates; i++)
                        {
                        *(clP++) = preLikeR[a++] * preLikeL[b++];
                        }
                    }
                }
            break;
        }

    return NO_ERROR;
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeDown_NY98_SSE: codon model with omega variation
|
-----------------------------------------------------------------*/
int CondLikeDown_NY98_SSE (TreeNode *p, int division, int chain)
{
    int             c, c1, h, i, j, k, t, shortCut, *lState=NULL, *rState=NULL, nStates, nStatesSquared;
    CLFlt           *pL, *pR, *tiPL, *tiPR;
    __m128          *clL, *clR, *clP;
    __m128          mTiPL, mTiPR, mL, mR, mAcumL, mAcumR;
    ModelInfo       *m;
    CLFlt           *preLikeRV[4];
    CLFlt           *preLikeLV[4];
#   if !defined (DEBUG_NOSHORTCUTS)
    int             a;
#   endif
    
    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;

    /* Flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }
#   endif

    switch (shortCut)
        {
        case 0:
            tiPL = pL;
            tiPR = pR;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numVecChars; c++)
                    {
                    for (i=h=0; i<nStates; i++)
                        {
                        mAcumL = _mm_setzero_ps();
                        mAcumR = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h]);
                            mTiPR  = _mm_load1_ps (&tiPR[h++]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            }
                        *(clP++) = _mm_mul_ps (mAcumL, mAcumR);
                        }
                    clL += nStates;
                    clR += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeLV[c1] = &preLikeL[lState[t] + k*(nStatesSquared+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following statment we assume that SSE register can hold exactly 4 ClFlts. */
                        mAcumL = _mm_set_ps (*(preLikeLV[3]++), *(preLikeLV[2]++), *(preLikeLV[1]++), *(preLikeLV[0]++));
                        mAcumR = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPR  = _mm_load1_ps (&tiPR[h++]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            }
                        *(clP++) = _mm_mul_ps (mAcumL,mAcumR);
                        }
                    clR += nStates;
                    }
                tiPR += nStatesSquared;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeRV[c1] = &preLikeR[rState[t] + k*(nStatesSquared+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following statment we assume that SSE register can hold exactly 4 ClFlts. */
                        mAcumR = _mm_set_ps (*(preLikeRV[3]++), *(preLikeRV[2]++), *(preLikeRV[1]++), *(preLikeRV[0]++));
                        mAcumL = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h++]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            }
                        *(clP++) = _mm_mul_ps (mAcumL,mAcumR);
                        }
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                }
            break;
        case 3:
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeRV[c1] = &preLikeR[rState[t] + k*(nStatesSquared+nStates)];
                        preLikeLV[c1] = &preLikeL[lState[t] + k*(nStatesSquared+nStates)];
                        }
                    for (i=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following 2 statments we assume that SSE register can hold exactly 4 ClFlts. */
                        mL = _mm_set_ps (*(preLikeLV[3]++), *(preLikeLV[2]++), *(preLikeLV[1]++), *(preLikeLV[0]++));
                        mR = _mm_set_ps (*(preLikeRV[3]++), *(preLikeRV[2]++), *(preLikeRV[1]++), *(preLikeRV[0]++));
                        *(clP++) = _mm_mul_ps (mL,mR);
                        }
                    }
                }
            break;
        }

    return NO_ERROR;
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeDown_Std: variable number of states model
|       with or without rate variation
|
-----------------------------------------------------------------*/
int CondLikeDown_Std (TreeNode *p, int division, int chain)
{
    int             a, c, h, i, j, k, nStates, nCats, tmp;
    CLFlt           *clL, *clR, *clP, *pL, *pR, *tiPL, *tiPR, likeL, likeR;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* Flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];

    /* Conditional likelihood space is assumed to be arranged in numGammaCats blocks of data. Each block contains all data for one gamma category.
    Each gamma cat block consist of numChars sequences of data, each of this sequences corresponds to a character of data matrix. 
    A sequence consists of nStates for all non-binary data, otherwise length of sequence is nStates*numBetaCats (i.e. 2*numBetaCats) */

    /* calculate ancestral probabilities */
    for (k=h=0; k<m->numRateCats; k++)
        {
        /* calculate ancestral probabilities */
        for (c=0; c<m->numChars; c++)
            {
            nStates = m->nStates[c];
        
            /* the following lines ensure that nCats is 1 unless */
            /* the character is binary and beta categories are used  */
            if (nStates == 2)
                nCats = m->numBetaCats;
            else
                nCats = 1;

            tmp = k*nStates*nStates; /* tmp contains offset to skip rate cats that already processed*/
            tiPL = pL + m->tiIndex[c] + tmp;
            tiPR = pR + m->tiIndex[c] + tmp;
            tmp = (m->numRateCats-1)*2*2; /* tmp contains size of block of tpi matrices across all rate cats (minus one) for single beta category. Further used only if character is binary to jump to next beta category */
                
            for (j=0; j<nCats;j++)
                {
                for (a=0; a<nStates; a++)
                    {
                    likeL = likeR = 0.0;
                    for (i=0; i<nStates; i++)
                        {
                        likeL += *(tiPL++) * clL[i];
                        likeR += *(tiPR++) * clR[i];
                        }
                    clP[h++] = likeL * likeR;
                    }
                clL += nStates;
                clR += nStates;
        
                tiPL += tmp;
                tiPR += tmp;
                }
            }
        }

    return NO_ERROR;
}


#if !defined (SSE_ENABLED) || 1
/*----------------------------------------------------------------
|
|   CondLikeRoot_Bin: binary model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeRoot_Bin (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *clL, *clR, *clP, *clA, *pL, *pR, *pA, *tiPL, *tiPR, *tiPA;
    ModelInfo       *m;

    /* find model settings for this division */
    m = &modelSettings[division];
    
    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    tiPL = pL;
    tiPR = pR;
    tiPA = pA;
    for (k=0; k<m->numRateCats; k++)
        {
        for (c=0; c<m->numChars; c++)
            {
            *(clP++) = (tiPL[0]*clL[0] + tiPL[1]*clL[1])
                      *(tiPR[0]*clR[0] + tiPR[1]*clR[1])
                      *(tiPA[0]*clA[0] + tiPA[1]*clA[1]);
            *(clP++) = (tiPL[2]*clL[0] + tiPL[3]*clL[1])
                      *(tiPR[2]*clR[0] + tiPR[3]*clR[1])
                      *(tiPA[2]*clA[0] + tiPA[3]*clA[1]);

            clA += 2;
            clL += 2;
            clR += 2;
            }
        tiPA += 4;
        tiPL += 4;
        tiPR += 4;
        }

    return NO_ERROR;
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeRoot_Bin_SSE:binary model with or without rate
|       variation 
|
-----------------------------------------------------------------*/
int CondLikeRoot_Bin_SSE (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *pL, *pR, *pA, *tiPL, *tiPR, *tiPA;
    __m128          *clL, *clR, *clP, *clA;
    __m128          m1, m2, m3, m4, m5, m6, m7;
    ModelInfo       *m;

    m = &modelSettings[division];

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    tiPL = pL;
    tiPR = pR;
    tiPA = pA;
    for (k=0; k<m->numRateCats; k++)
        {
        for (c=0; c<m->numVecChars; c++)
            {
            m1 = _mm_load1_ps (&tiPL[0]);
            m5 = *clL++;
            m2 = _mm_mul_ps (m1, m5);
            m1 = _mm_load1_ps (&tiPL[2]);
            m6 = _mm_mul_ps (m1, m5);

            m1 = _mm_load1_ps (&tiPL[1]);
            m5 = *clL++;
            m3 = _mm_mul_ps (m1, m5);
            m1 = _mm_load1_ps (&tiPL[3]);
            m5 = _mm_mul_ps (m1, m5);

            m4 = _mm_add_ps (m2, m3); /* in m4 we get (tiPL[0]*clL[0] + tiPL[1]*clL[1]) */
            m6 = _mm_add_ps (m5, m6); /* in m6 we get (tiPL[2]*clL[0] + tiPL[3]*clL[1]) */

            m1 = _mm_load1_ps (&tiPR[0]);
            m5 = *clR++;
            m2 = _mm_mul_ps (m1, m5);
            m1 = _mm_load1_ps (&tiPR[2]);
            m7 = _mm_mul_ps (m1, m5);

            m1 = _mm_load1_ps (&tiPR[1]);
            m5 = *clR++;
            m3 = _mm_mul_ps (m1, m5);
            m1 = _mm_load1_ps (&tiPR[3]);
            m5 = _mm_mul_ps (m1, m5);

            m1 = _mm_add_ps (m2, m3); /* in m1 we get (tiPR[0]*clR[0] + tiPR[1]*clR[1]) */
            m7 = _mm_add_ps (m5, m7); /* in m7 we get (tiPR[2]*clR[0] + tiPR[3]*clR[1]) */

            m4 = _mm_mul_ps (m1, m4); /* in m4 we get (tiPL[0]*clL[0] + tiPL[1]*clL[1])*(tiPR[0]*clR[0] + tiPR[1]*clR[1]) */
            m7 = _mm_mul_ps (m6, m7); /* in m7 we get (tiPL[2]*clL[0] + tiPL[3]*clL[1])*(tiPR[2]*clR[0] + tiPR[3]*clR[1]) */

            m1 = _mm_load1_ps (&tiPA[0]);
            m5 = *clA++;
            m2 = _mm_mul_ps (m1, m5);
            m1 = _mm_load1_ps (&tiPA[2]);
            m6 = _mm_mul_ps (m1, m5);

            m1 = _mm_load1_ps (&tiPA[1]);
            m5 = *clA++;
            m3 = _mm_mul_ps (m1, m5);
            m1 = _mm_load1_ps (&tiPA[3]);
            m1 = _mm_mul_ps (m1, m5);

            m2 = _mm_add_ps (m2, m3); /* in m1 we get (tiPA[0]*clA[0] + tiPA[1]*clA[1]) */
            m1 = _mm_add_ps (m1, m6); /* in m1 we get (tiPA[2]*clA[0] + tiPA[3]*clA[1]) */

            *clP++ = _mm_mul_ps (m2, m4);
            *clP++ = _mm_mul_ps (m1, m7);

            }
        tiPL += 4;
        tiPR += 4;
        tiPA += 4;
        }

    return NO_ERROR;
    
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeRoot_Gen: general n-state model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeRoot_Gen (TreeNode *p, int division, int chain)
{
    int             a, b, c, d, h, i, j, k, shortCut, *lState=NULL, *rState=NULL, *aState=NULL,
                    nObsStates, nStates, nStatesSquared, preLikeJump;
    CLFlt           likeL, likeR, likeA, *clL, *clR, *clP, *clA, *pL, *pR, *pA,
                    *tiPL, *tiPR, *tiPA;
    ModelInfo       *m;
#   if !defined (DEBUG_NOSHORTCUTS)
    int catStart;
#   endif
    
    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nObsStates = m->numStates;
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;
    preLikeJump = nObsStates * nStates;

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeL[a++] += tiPL[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeR[a++] += tiPR[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for anc branch, always terminal */
    if (m->isPartAmbig[p->anc->index] == YES)
        {
        shortCut = 4;
        }
    else 
        {
        aState = m->termState[p->anc->index];
        tiPA = pA;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeA[a++] = tiPA[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeA[a++] += tiPA[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeA[a++] = 1.0;
            tiPA += nStatesSquared;
            }
        }
#   else
    shortCut = 4;
#   endif

    //shortCut = 4;
    switch (shortCut)
        {
        case 4:
            tiPL = pL;
            tiPR = pR;
            tiPA = pA;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    for (i=h=0; i<nStates; i++)
                        {
                        likeL = likeR = likeA = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += tiPL[h]*clL[j];
                            likeR += tiPR[h]*clR[j];
                            likeA += tiPA[h++]*clA[j];
                            }
                        *(clP++) = likeL * likeR * likeA;
                        }
                    clL += nStates;
                    clR += nStates;
                    clA += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                tiPA += nStatesSquared;
                }
            break;
        case 0:
            tiPR = pR;
            tiPL = pL;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = aState[c] + k*(preLikeJump+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeR = likeL = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeR += tiPR[h]*clR[j];
                            likeL += tiPL[h++]*clL[j];
                            }
                        *(clP++) = preLikeA[a++] * likeR * likeL;
                        }
                    clR += nStates;
                    clL += nStates;
                    }
                tiPR += nStatesSquared;
                tiPL += nStatesSquared;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = lState[c] + k*(preLikeJump+nStates);
                    b = aState[c] + k*(preLikeJump+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeR = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeR += tiPR[h++]*clR[j];
                            }
                        *(clP++) = preLikeL[a++] * preLikeA[b++] * likeR;
                        }
                    clR += nStates;
                    }
                tiPR += nStatesSquared;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = rState[c] + k*(preLikeJump+nStates);
                    b = aState[c] + k*(preLikeJump+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeL = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += tiPL[h++]*clL[j];
                            }
                        *(clP++) = preLikeR[a++] * preLikeA[b++] * likeL;
                        }
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                }
            break;  
        case 3:
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = rState[c] + k*(preLikeJump+nStates);
                    b = lState[c] + k*(preLikeJump+nStates);
                    d = aState[c] + k*(preLikeJump+nStates);
                    for (i=0; i<nStates; i++)
                        {
                        *(clP++) = preLikeR[a++] * preLikeL[b++] * preLikeA[d++];
                        }
                    }
                }
            break;
        }

    return NO_ERROR;
}


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeRoot_Gen_SSE:general n-state model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeRoot_Gen_SSE (TreeNode *p, int division, int chain)
{
    int             c, c1, t, h, i, j, k, shortCut, *lState=NULL, *rState=NULL, *aState=NULL, nObsStates, preLikeJump,
                    nStates, nStatesSquared;
    CLFlt           *pL, *pR, *pA,
                    *tiPL, *tiPR, *tiPA;
    __m128          *clL, *clR, *clP, *clA;
    __m128          mTiPL, mTiPR, mTiPA, mL, mR, mA, mAcumL, mAcumR, mAcumA;
    ModelInfo       *m;
    CLFlt           *preLikeRV[4];
    CLFlt           *preLikeLV[4];
    CLFlt           *preLikeAV[4];

#   if !defined (DEBUG_NOSHORTCUTS)
    int a, b, catStart;
#   endif

    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nObsStates = m->numStates;
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;
    preLikeJump = nObsStates * nStates;

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeL[a++] += tiPL[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeR[a++] += tiPR[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for anc branch, always terminal */
    if (m->isPartAmbig[p->anc->index] == YES)
        {
        shortCut = 4;
        }
    else 
        {
        aState = m->termState[p->anc->index];
        tiPA = pA;
        for (k=a=0; k<m->numRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeA[a++] = tiPA[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeA[a++] += tiPA[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeA[a++] = 1.0;
            tiPA += nStatesSquared;
            }
        }
#   else
    shortCut = 4;
#   endif

        switch (shortCut)
        {
        case 4:
            tiPL = pL;
            tiPR = pR;
            tiPA = pA;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=0; c<m->numVecChars; c++)
                    {
                    for (i=h=0; i<nStates; i++)
                        {
                        mAcumL = _mm_setzero_ps();
                        mAcumR = _mm_setzero_ps();
                        mAcumA = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h]);
                            mTiPR  = _mm_load1_ps (&tiPR[h]);
                            mTiPA  = _mm_load1_ps (&tiPA[h++]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mA     = _mm_mul_ps (mTiPA, clA[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            mAcumA = _mm_add_ps (mA, mAcumA);
                            }
                        mAcumL = _mm_mul_ps (mAcumL, mAcumR);
                        *(clP++) = _mm_mul_ps (mAcumL, mAcumA);
                        }
                    clL += nStates;
                    clR += nStates;
                    clA += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                tiPA += nStatesSquared;
                }
            break;
        case 0:
            tiPL =pL;
            tiPR =pR;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeAV[c1] = &preLikeA[aState[t] + k*(preLikeJump+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following statment we assume that SSE register can hold exactly 4 ClFlts. */
                        mAcumA = _mm_set_ps (*(preLikeAV[3]++), *(preLikeAV[2]++), *(preLikeAV[1]++), *(preLikeAV[0]++));
                        mAcumL = _mm_setzero_ps();
                        mAcumR = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            mTiPR  = _mm_load1_ps (&tiPR[h++]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            }
                        mAcumL = _mm_mul_ps (mAcumL, mAcumR);
                        *(clP++) = _mm_mul_ps (mAcumL, mAcumA);
                        }
                    clR += nStates;
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeLV[c1] = &preLikeL[lState[t] + k*(preLikeJump+nStates)];
                        preLikeAV[c1] = &preLikeA[aState[t] + k*(preLikeJump+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following statment we assume that SSE register can hold exactly 4 ClFlts. */
                        mAcumL = _mm_set_ps (*(preLikeLV[3]++), *(preLikeLV[2]++), *(preLikeLV[1]++), *(preLikeLV[0]++));
                        mAcumA = _mm_set_ps (*(preLikeAV[3]++), *(preLikeAV[2]++), *(preLikeAV[1]++), *(preLikeAV[0]++));
                        mAcumR = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPR  = _mm_load1_ps (&tiPR[h++]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            }
                        mAcumL = _mm_mul_ps (mAcumL, mAcumR);
                        *(clP++) = _mm_mul_ps (mAcumL, mAcumA);
                        }
                    clR += nStates;
                    }
                tiPR += nStatesSquared;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeRV[c1] = &preLikeR[rState[t] + k*(preLikeJump+nStates)];
                        preLikeAV[c1] = &preLikeA[aState[t] + k*(preLikeJump+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following statment we assume that SSE register can hold exactly 4 ClFlts. */
                        mAcumR = _mm_set_ps (*(preLikeRV[3]++), *(preLikeRV[2]++), *(preLikeRV[1]++), *(preLikeRV[0]++));
                        mAcumA = _mm_set_ps (*(preLikeAV[3]++), *(preLikeAV[2]++), *(preLikeAV[1]++), *(preLikeAV[0]++));
                        mAcumL = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h++]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            }
                        mAcumL = _mm_mul_ps (mAcumL, mAcumR);
                        *(clP++) = _mm_mul_ps (mAcumL,mAcumA);
                        }
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                }
            break;
        case 3:
            for (k=0; k<m->numRateCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeRV[c1] = &preLikeR[rState[t] + k*(preLikeJump+nStates)];
                        preLikeLV[c1] = &preLikeL[lState[t] + k*(preLikeJump+nStates)];
                        preLikeAV[c1] = &preLikeA[aState[t] + k*(preLikeJump+nStates)];
                        }
                    for (i=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following 2 statments we assume that SSE register can hold exactly 4 ClFlts. */
                        mL = _mm_set_ps (*(preLikeLV[3]++), *(preLikeLV[2]++), *(preLikeLV[1]++), *(preLikeLV[0]++));
                        mR = _mm_set_ps (*(preLikeRV[3]++), *(preLikeRV[2]++), *(preLikeRV[1]++), *(preLikeRV[0]++));
                        mA = _mm_set_ps (*(preLikeAV[3]++), *(preLikeAV[2]++), *(preLikeAV[1]++), *(preLikeAV[0]++));
                        mL = _mm_mul_ps (mL,mR);
                        *(clP++) = _mm_mul_ps (mL,mA);
                        }
                    }
                }
            break;
        }

    return NO_ERROR;
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeRoot_Gen_GibbsGamma: general n-state model with rate
|       variation modeled using a discrete gamma distribution with
|       Gibbs resampling of rate categories
|
-----------------------------------------------------------------*/
int CondLikeRoot_Gen_GibbsGamma (TreeNode *p, int division, int chain)
{
    int             a, b, c, i, j, r, *rateCat, shortCut, *lState=NULL,
                    *rState=NULL, *aState=NULL, nObsStates, nStates,
                    nStatesSquared, nRateCats;
    CLFlt           likeL, likeR, likeA, *clL, *clR, *clP, *clA, *pL, *pR, *pA,
                    *tiPL, *tiPR, *tiPA;
    ModelInfo       *m;
#   if !defined (DEBUG_NOSHORTCUTS)
    int k, catStart;
#endif
    
    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nObsStates = m->numStates;
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;

    /* flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);

    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* find rate category index and number of rate categories */
    rateCat = m->tiIndex + chain * m->numChars;
    nRateCats = m->numRateCats;

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<nRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeL[a++] += tiPL[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<nRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeR[a++] += tiPR[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for anc branch, always terminal */
    if (m->isPartAmbig[p->anc->index] == YES)
        {
        shortCut = 4;
        }
    else 
        {
        aState = m->termState[p->anc->index];
        tiPA = pA;
        for (k=a=0; k<nRateCats; k++)
            {
            catStart = a;
            for (i=0; i<nObsStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeA[a++] = tiPA[j];
            for (b=1; b<nStates/nObsStates; b++)
                {
                a = catStart;
                for (i=0; i<nObsStates; i++)
                    {
                    for (j=i+b*nObsStates; j<nStatesSquared; j+=nStates)
                        preLikeA[a++] += tiPA[j];
                    }
                }
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeA[a++] = 1.0;
            tiPA += nStatesSquared;
            }
        }
#   else
    shortCut = 4;
#   endif

    switch (shortCut)
        {
    case 4:
        for (c=0; c<m->numChars; c++)
            {
            r = (*rateCat++);
            if (r < nRateCats)
                {
                tiPL = pL + r*nStatesSquared;
                tiPR = pR + r*nStatesSquared;
                tiPA = pA + r*nStatesSquared;
                for (i=0; i<nStates; i++)
                    {
                    likeL = likeR = likeA = 0.0;
                    for (j=0; j<nStates; j++)
                        {
                        likeL += (*tiPL++) * clL[j];
                        likeR += (*tiPR++) * clR[j];
                        likeA += (*tiPA++) * clA[j];
                        }
                    *(clP++) = likeL * likeR * likeA;
                    }
                }
            else
                clP += nStates;
            clL += nStates;
            clR += nStates;
            clA += nStates;
            }
        break;
    case 0:
    case 3:
        for (c=0; c<m->numChars; c++)
            {
            r = (*rateCat++);
            if (r < nRateCats)
                {
                tiPL = pL + r*nStatesSquared;
                tiPR = pR + r*nStatesSquared;
                a = aState[c] + r*(nStatesSquared+nStates);
                for (i=0; i<nStates; i++)
                    {
                    likeL = likeR = 0.0;
                    for (j=0; j<nStates; j++)
                        {
                        likeL += (*tiPL++) * clL[j];
                        likeR += (*tiPR++) * clR[j];
                        }
                    *(clP++) = likeL * likeR * preLikeA[a++];
                    }
                }
            else
                clP += nStates;
            clL += nStates;
            clR += nStates;
            }
        break;
    case 1:
        for (c=0; c<m->numChars; c++)
            {
            r = (*rateCat++);
            if (r < nRateCats)
                {
                tiPR = pR + r*nStatesSquared;
                a = lState[c] + r*(nStatesSquared+nStates);
                b = aState[c] + r*(nStatesSquared+nStates);
                for (i=0; i<nStates; i++)
                    {
                    likeR = 0.0;
                    for (j=0; j<nStates; j++)
                        {
                        likeR += (*tiPR++) * clR[j];
                        }
                    *(clP++) = preLikeL[a++] * likeR * preLikeA[b++];
                    }
                }
            else
                clP += nStates;
            clR += nStates;
            }
        break;
    case 2:
        for (c=0; c<m->numChars; c++)
            {
            r = (*rateCat++);
            if (r < nRateCats)
                {
                tiPL = pL + r*nStatesSquared;
                a = rState[c] + r*(nStatesSquared+nStates);
                b = aState[c] + r*(nStatesSquared+nStates);
                for (i=0; i<nStates; i++)
                    {
                    likeL = 0.0;
                    for (j=0; j<nStates; j++)
                        {
                        likeL += (*tiPL++) * clL[j];
                        }
                    *(clP++) = likeL * preLikeR[a++] * preLikeA[b++];
                    }
                }
            else
                clP += nStates;
            clL += nStates;
            }
        break;
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeRoot_NUC4: 4by4 nucleotide model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeRoot_NUC4 (TreeNode *p, int division, int chain)
{
    int             a, c, h, i, j, k, shortCut, *lState=NULL, *rState=NULL, *aState=NULL;
    CLFlt           *clL, *clR, *clP, *clA, *pL, *pR, *pA, *tiPL, *tiPR, *tiPA;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=j=0; k<m->numRateCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeL[j++] = tiPL[0];
                preLikeL[j++] = tiPL[4];
                preLikeL[j++] = tiPL[8];
                preLikeL[j++] = tiPL[12];
                tiPL++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeL[j++] = 1.0;
            tiPL += 12;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=j=0; k<m->numRateCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeR[j++] = tiPR[0];
                preLikeR[j++] = tiPR[4];
                preLikeR[j++] = tiPR[8];
                preLikeR[j++] = tiPR[12];
                tiPR++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeR[j++] = 1.0;
            tiPR += 12;
            }
        }

    /* find likelihoods of site patterns for anc branch, always terminal */
    if (m->isPartAmbig[p->anc->index] == YES)
        {
        shortCut = 4;
        }
    else 
        {
        aState = m->termState[p->anc->index];
        tiPA = pA;
        for (k=j=0; k<m->numRateCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeA[j++] = tiPA[0];
                preLikeA[j++] = tiPA[4];
                preLikeA[j++] = tiPA[8];
                preLikeA[j++] = tiPA[12];
                tiPA++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeA[j++] = 1.0;
            tiPA += 12;
            }
        }
#   else
    shortCut = 4;
#   endif

    switch (shortCut)
        {
    case 4:
        tiPL = pL;
        tiPR = pR;
        tiPA = pA;
        for (k=h=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                            *(tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T])
                            *(tiPA[AA]*clA[A] + tiPA[AC]*clA[C] + tiPA[AG]*clA[G] + tiPA[AT]*clA[T]);
                clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                            *(tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T])
                            *(tiPA[CA]*clA[A] + tiPA[CC]*clA[C] + tiPA[CG]*clA[G] + tiPA[CT]*clA[T]);
                clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                            *(tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T])
                            *(tiPA[GA]*clA[A] + tiPA[GC]*clA[C] + tiPA[GG]*clA[G] + tiPA[GT]*clA[T]);
                clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                            *(tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T])
                            *(tiPA[TA]*clA[A] + tiPA[TC]*clA[C] + tiPA[TG]*clA[G] + tiPA[TT]*clA[T]);
                clL += 4;
                clR += 4;
                clA += 4;
                }
            tiPL += 16;
            tiPR += 16;
            tiPA += 16;
            }
        break;

    case 0:
        tiPL = pL;
        tiPR = pR;
        for (k=h=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                i = aState[c] + k*20;
                clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                            *(tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T])
                            *preLikeA[i++];
                clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                            *(tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T])
                            *preLikeA[i++];
                clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                            *(tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T])
                            *preLikeA[i++];
                clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                            *(tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T])
                            *preLikeA[i++];
                clL += 4;
                clR += 4;
                }
            tiPL += 16;
            tiPR += 16;
            }
        break;

    case 1:
        tiPR = pR;
        for (k=h=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                i = lState[c] + k*20;
                j = aState[c] + k*20;
                clP[h++] =   (tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T])
                            *preLikeL[i++]*preLikeA[j++];
                clP[h++] =   (tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T])
                            *preLikeL[i++]*preLikeA[j++];
                clP[h++] =   (tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T])
                            *preLikeL[i++]*preLikeA[j++];
                clP[h++] =   (tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T])
                            *preLikeL[i++]*preLikeA[j++];
                clR += 4;
                }
            tiPR += 16;
            }
        break;

    case 2:
        tiPL = pL;
        for (k=h=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                i = rState[c] + k*20;
                j = aState[c] + k*20;
                clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                            *preLikeR[i++]*preLikeA[j++];
                clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                            *preLikeR[i++]*preLikeA[j++];
                clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                            *preLikeR[i++]*preLikeA[j++];
                clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                            *preLikeR[i++]*preLikeA[j++];
                clL += 4;
                }
            tiPL += 16;
            }
        break;

    case 3:
        for (k=h=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                a = lState[c] + k*20;
                i = rState[c] + k*20;
                j = aState[c] + k*20;
                clP[h++] =   preLikeL[a++]*preLikeR[i++]*preLikeA[j++];
                clP[h++] =   preLikeL[a++]*preLikeR[i++]*preLikeA[j++];
                clP[h++] =   preLikeL[a++]*preLikeR[i++]*preLikeA[j++];
                clP[h++] =   preLikeL[a++]*preLikeR[i++]*preLikeA[j++];
                }
            }
        break;
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeRoot_NUC4_GibbsGamma: 4by4 nucleotide model with rate
|       variation approimated by Gibbs sampling from gamma
|
-----------------------------------------------------------------*/
int CondLikeRoot_NUC4_GibbsGamma (TreeNode *p, int division, int chain)
{
    int             c, h, i, j, r, *rateCat, shortCut, *lState=NULL, *rState=NULL, *aState=NULL,
                    nRateCats;
    CLFlt           *clL, *clR, *clP, *clA, *pL, *pR, *pA, *tiPL, *tiPR, *tiPA;
    ModelInfo       *m;
#   if !defined (DEBUG_NOSHORTCUTS)
    int k;
#   endif
    
    m = &modelSettings[division];

    /* flip conditional likelihood space */
    FlipCondLikeSpace (m, chain, p->index);

        /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* find rate category index and number of gamma categories */
    rateCat = m->tiIndex + chain * m->numChars;
    nRateCats = m->numRateCats;

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=j=0; k<nRateCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeL[j++] = tiPL[0];
                preLikeL[j++] = tiPL[4];
                preLikeL[j++] = tiPL[8];
                preLikeL[j++] = tiPL[12];
                tiPL++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeL[j++] = 1.0;
            tiPL += 12;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=j=0; k<nRateCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeR[j++] = tiPR[0];
                preLikeR[j++] = tiPR[4];
                preLikeR[j++] = tiPR[8];
                preLikeR[j++] = tiPR[12];
                tiPR++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeR[j++] = 1.0;
            tiPR += 12;
            }
        }

    /* find likelihoods of site patterns for anc branch, always terminal */
    if (m->isPartAmbig[p->anc->index] == YES)
        {
        shortCut = 4;
        }
    else 
        {
        aState = m->termState[p->anc->index];
        tiPA = pA;
        for (k=j=0; k<nRateCats; k++)
            {
            for (i=0; i<4; i++)
                {
                preLikeA[j++] = tiPA[0];
                preLikeA[j++] = tiPA[4];
                preLikeA[j++] = tiPA[8];
                preLikeA[j++] = tiPA[12];
                tiPA++;
                }
            /* for ambiguous */
            for (i=0; i<4; i++)
                preLikeA[j++] = 1.0;
            tiPA += 12;
            }
        }
#   else
    shortCut = 4;
#   endif

    switch (shortCut)
        {
    case 4:
        for (c=h=0; c<m->numChars; c++)
            {
            r = rateCat[c];
            if (r < nRateCats)
                {
                tiPL = pL + r * 16;
                tiPR = pR + r * 16;
                tiPA = pA + r * 16;
                clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                            *(tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T])
                            *(tiPA[AA]*clA[A] + tiPA[AC]*clA[C] + tiPA[AG]*clA[G] + tiPA[AT]*clA[T]);
                clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                            *(tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T])
                            *(tiPA[CA]*clA[A] + tiPA[CC]*clA[C] + tiPA[CG]*clA[G] + tiPA[CT]*clA[T]);
                clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                            *(tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T])
                            *(tiPA[GA]*clA[A] + tiPA[GC]*clA[C] + tiPA[GG]*clA[G] + tiPA[GT]*clA[T]);
                clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                            *(tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T])
                            *(tiPA[TA]*clA[A] + tiPA[TC]*clA[C] + tiPA[TG]*clA[G] + tiPA[TT]*clA[T]);
                }
            else
                h += 4;
            clL += 4;
            clR += 4;
            clA += 4;
            }
        break;

    case 0:
    case 3:
        for (c=h=0; c<m->numChars; c++)
            {
            r = rateCat[c];
            if (r < nRateCats)
                {
                tiPL = pL + r * 16;
                tiPR = pR + r * 16;
                i = aState[c] + r * 20;
                clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                            *(tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T])
                            *preLikeA[i++];
                clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                            *(tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T])
                            *preLikeA[i++];
                clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                            *(tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T])
                            *preLikeA[i++];
                clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                            *(tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T])
                            *preLikeA[i++];
                }
            else
                h += 4;
            clL += 4;
            clR += 4;
            }
        break;

    case 1:
        for (c=h=0; c<m->numChars; c++)
            {
            r = rateCat[c];
            if (r < nRateCats)
                {
                tiPR = pR + r * 16;
                i = lState[c] + r * 20;
                j = aState[c] + r * 20;
                clP[h++] =   (tiPR[AA]*clR[A] + tiPR[AC]*clR[C] + tiPR[AG]*clR[G] + tiPR[AT]*clR[T])
                            *preLikeL[i++]*preLikeA[j++];
                clP[h++] =   (tiPR[CA]*clR[A] + tiPR[CC]*clR[C] + tiPR[CG]*clR[G] + tiPR[CT]*clR[T])
                            *preLikeL[i++]*preLikeA[j++];
                clP[h++] =   (tiPR[GA]*clR[A] + tiPR[GC]*clR[C] + tiPR[GG]*clR[G] + tiPR[GT]*clR[T])
                            *preLikeL[i++]*preLikeA[j++];
                clP[h++] =   (tiPR[TA]*clR[A] + tiPR[TC]*clR[C] + tiPR[TG]*clR[G] + tiPR[TT]*clR[T])
                            *preLikeL[i++]*preLikeA[j++];
                }
            else
                h += 4;
            clR += 4;
            }
        break;

    case 2:
        for (c=h=0; c<m->numChars; c++)
            {
            r = rateCat[c];
            if (r < nRateCats)
                {
                tiPL = pL + r * 16;
                i = rState[c] + r * 20;
                j = aState[c] + r * 20;
                clP[h++] =   (tiPL[AA]*clL[A] + tiPL[AC]*clL[C] + tiPL[AG]*clL[G] + tiPL[AT]*clL[T])
                            *preLikeR[i++]*preLikeA[j++];
                clP[h++] =   (tiPL[CA]*clL[A] + tiPL[CC]*clL[C] + tiPL[CG]*clL[G] + tiPL[CT]*clL[T])
                            *preLikeR[i++]*preLikeA[j++];
                clP[h++] =   (tiPL[GA]*clL[A] + tiPL[GC]*clL[C] + tiPL[GG]*clL[G] + tiPL[GT]*clL[T])
                            *preLikeR[i++]*preLikeA[j++];
                clP[h++] =   (tiPL[TA]*clL[A] + tiPL[TC]*clL[C] + tiPL[TG]*clL[G] + tiPL[TT]*clL[T])
                            *preLikeR[i++]*preLikeA[j++];
                }
            else
                h += 4;
            clL += 4;
            }
        break;
        }

    return NO_ERROR;
}


#if defined (FMA_ENABLED)
/*----------------------------------------------------------------
 |
 |   CondLikeRoot_NUC4_FMA: 4by4 nucleotide model with or without rate
 |       variation using AVX + FMA instructions
 |
 -----------------------------------------------------------------*/
int CondLikeRoot_NUC4_FMA (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *pL, *pR, *pA, *tiPL, *tiPR, *tiPA;
    __m256          *clL, *clR, *clP, *clA;
    __m256          m1, m2, m3, m4, m5, m6;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    
    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];
    
    tiPL = pL;
    tiPR = pR;
    tiPA = pA;
    for (k=0; k<m->numRateCats; k++)
    {
        for (c=0; c<m->numVecChars; c++)
        {
            m1 = _mm256_broadcast_ss (&tiPL[AA]);
            m2 = _mm256_broadcast_ss (&tiPR[AA]);
            m3 = _mm256_broadcast_ss (&tiPA[AA]);
            m4 = _mm256_mul_ps (m1, clL[A]);
            m5 = _mm256_mul_ps (m2, clR[A]);
            m6 = _mm256_mul_ps (m3, clA[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[AC]);
            m2 = _mm256_broadcast_ss (&tiPR[AC]);
            m3 = _mm256_broadcast_ss (&tiPA[AC]);
            m4 = _mm256_fmadd_ps (m1, clL[C], m4);
            m5 = _mm256_fmadd_ps (m2, clR[C], m5);
            m6 = _mm256_fmadd_ps (m3, clA[C], m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[AG]);
            m2 = _mm256_broadcast_ss (&tiPR[AG]);
            m3 = _mm256_broadcast_ss (&tiPA[AG]);
            m4 = _mm256_fmadd_ps (m1, clL[G], m4);
            m5 = _mm256_fmadd_ps (m2, clR[G], m5);
            m6 = _mm256_fmadd_ps (m3, clA[G], m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[AT]);
            m2 = _mm256_broadcast_ss (&tiPR[AT]);
            m3 = _mm256_broadcast_ss (&tiPA[AT]);
            m4 = _mm256_fmadd_ps (m1, clL[T], m4);
            m5 = _mm256_fmadd_ps (m2, clR[T], m5);
            m6 = _mm256_fmadd_ps (m3, clA[T], m6);
            
            m4 = _mm256_mul_ps (m4, m5);
            *clP++ = _mm256_mul_ps (m4, m6);
           
            m1 = _mm256_broadcast_ss (&tiPL[CA]);
            m2 = _mm256_broadcast_ss (&tiPR[CA]);
            m3 = _mm256_broadcast_ss (&tiPA[CA]);
            m4 = _mm256_mul_ps (m1, clL[A]);
            m5 = _mm256_mul_ps (m2, clR[A]);
            m6 = _mm256_mul_ps (m3, clA[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[CC]);
            m2 = _mm256_broadcast_ss (&tiPR[CC]);
            m3 = _mm256_broadcast_ss (&tiPA[CC]);
            m4 = _mm256_fmadd_ps (m1, clL[C], m4);
            m5 = _mm256_fmadd_ps (m2, clR[C], m5);
            m6 = _mm256_fmadd_ps (m3, clA[C], m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[CG]);
            m2 = _mm256_broadcast_ss (&tiPR[CG]);
            m3 = _mm256_broadcast_ss (&tiPA[CG]);
            m4 = _mm256_fmadd_ps (m1, clL[G], m4);
            m5 = _mm256_fmadd_ps (m2, clR[G], m5);
            m6 = _mm256_fmadd_ps (m3, clA[G], m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[CT]);
            m2 = _mm256_broadcast_ss (&tiPR[CT]);
            m3 = _mm256_broadcast_ss (&tiPA[CT]);
            m4 = _mm256_fmadd_ps (m1, clL[T], m4);
            m5 = _mm256_fmadd_ps (m2, clR[T], m5);
            m6 = _mm256_fmadd_ps (m3, clA[T], m6);
            
            m4 = _mm256_mul_ps (m4, m5);
            *clP++ = _mm256_mul_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[GA]);
            m2 = _mm256_broadcast_ss (&tiPR[GA]);
            m3 = _mm256_broadcast_ss (&tiPA[GA]);
            m4 = _mm256_mul_ps (m1, clL[A]);
            m5 = _mm256_mul_ps (m2, clR[A]);
            m6 = _mm256_mul_ps (m3, clA[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[GC]);
            m2 = _mm256_broadcast_ss (&tiPR[GC]);
            m3 = _mm256_broadcast_ss (&tiPA[GC]);
            m4 = _mm256_fmadd_ps (m1, clL[C], m4);
            m5 = _mm256_fmadd_ps (m2, clR[C], m5);
            m6 = _mm256_fmadd_ps (m3, clA[C], m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[GG]);
            m2 = _mm256_broadcast_ss (&tiPR[GG]);
            m3 = _mm256_broadcast_ss (&tiPA[GG]);
            m4 = _mm256_fmadd_ps (m1, clL[G], m4);
            m5 = _mm256_fmadd_ps (m2, clR[G], m5);
            m6 = _mm256_fmadd_ps (m3, clA[G], m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[GT]);
            m2 = _mm256_broadcast_ss (&tiPR[GT]);
            m3 = _mm256_broadcast_ss (&tiPA[GT]);
            m4 = _mm256_fmadd_ps (m1, clL[T], m4);
            m5 = _mm256_fmadd_ps (m2, clR[T], m5);
            m6 = _mm256_fmadd_ps (m3, clA[T], m6);
            
            m4 = _mm256_mul_ps (m4, m5);
            *clP++ = _mm256_mul_ps (m4, m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[TA]);
            m2 = _mm256_broadcast_ss (&tiPR[TA]);
            m3 = _mm256_broadcast_ss (&tiPA[TA]);
            m4 = _mm256_mul_ps (m1, clL[A]);
            m5 = _mm256_mul_ps (m2, clR[A]);
            m6 = _mm256_mul_ps (m3, clA[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[TC]);
            m2 = _mm256_broadcast_ss (&tiPR[TC]);
            m3 = _mm256_broadcast_ss (&tiPA[TC]);
            m4 = _mm256_fmadd_ps (m1, clL[C], m4);
            m5 = _mm256_fmadd_ps (m2, clR[C], m5);
            m6 = _mm256_fmadd_ps (m3, clA[C], m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[TG]);
            m2 = _mm256_broadcast_ss (&tiPR[TG]);
            m3 = _mm256_broadcast_ss (&tiPA[TG]);
            m4 = _mm256_fmadd_ps (m1, clL[G], m4);
            m5 = _mm256_fmadd_ps (m2, clR[G], m5);
            m6 = _mm256_fmadd_ps (m3, clA[G], m6);
            
            m1 = _mm256_broadcast_ss (&tiPL[TT]);
            m2 = _mm256_broadcast_ss (&tiPR[TT]);
            m3 = _mm256_broadcast_ss (&tiPA[TT]);
            m4 = _mm256_fmadd_ps (m1, clL[T], m4);
            m5 = _mm256_fmadd_ps (m2, clR[T], m5);
            m6 = _mm256_fmadd_ps (m3, clA[T], m6);
            
            m4 = _mm256_mul_ps (m4, m5);
            *clP++ = _mm256_mul_ps (m4, m6);

            clL += 4;
            clR += 4;
            clA += 4;
        }
        tiPL += 16;
        tiPR += 16;
        tiPA += 16;
    }
    
    return NO_ERROR;
}
#endif


#if defined (AVX_ENABLED)
/*----------------------------------------------------------------
 |
 |   CondLikeRoot_NUC4_AVX: 4by4 nucleotide model with or without rate
 |       variation using AVX instructions
 |
 -----------------------------------------------------------------*/
int CondLikeRoot_NUC4_AVX (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *pL, *pR, *pA, *tiPL, *tiPR, *tiPA;
    __m256          *clL, *clR, *clP, *clA;
    __m256          m1, m2, m3, m4, m5, m6, m7, m8, m9;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    
    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];
    
    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];
    
    tiPL = pL;
    tiPR = pR;
    tiPA = pA;
    for (k=0; k<m->numRateCats; k++)
    {
        for (c=0; c<m->numVecChars; c++)
        {
            m1 = _mm256_broadcast_ss (&tiPL[AA]);
            m2 = _mm256_broadcast_ss (&tiPR[AA]);
            m3 = _mm256_broadcast_ss (&tiPA[AA]);
            m7 = _mm256_mul_ps (m1, clL[A]);
            m8 = _mm256_mul_ps (m2, clR[A]);
            m9 = _mm256_mul_ps (m3, clA[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[AC]);
            m2 = _mm256_broadcast_ss (&tiPR[AC]);
            m3 = _mm256_broadcast_ss (&tiPA[AC]);
            m4 = _mm256_mul_ps (m1, clL[C]);
            m5 = _mm256_mul_ps (m2, clR[C]);
            m6 = _mm256_mul_ps (m3, clA[C]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[AG]);
            m2 = _mm256_broadcast_ss (&tiPR[AG]);
            m3 = _mm256_broadcast_ss (&tiPA[AG]);
            m4 = _mm256_mul_ps (m1, clL[G]);
            m5 = _mm256_mul_ps (m2, clR[G]);
            m6 = _mm256_mul_ps (m3, clA[G]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[AT]);
            m2 = _mm256_broadcast_ss (&tiPR[AT]);
            m3 = _mm256_broadcast_ss (&tiPA[AT]);
            m4 = _mm256_mul_ps (m1, clL[T]);
            m5 = _mm256_mul_ps (m2, clR[T]);
            m6 = _mm256_mul_ps (m3, clA[T]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m7 = _mm256_mul_ps (m7, m8);
            *clP++ = _mm256_mul_ps (m7, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[CA]);
            m2 = _mm256_broadcast_ss (&tiPR[CA]);
            m3 = _mm256_broadcast_ss (&tiPA[CA]);
            m7 = _mm256_mul_ps (m1, clL[A]);
            m8 = _mm256_mul_ps (m2, clR[A]);
            m9 = _mm256_mul_ps (m3, clA[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[CC]);
            m2 = _mm256_broadcast_ss (&tiPR[CC]);
            m3 = _mm256_broadcast_ss (&tiPA[CC]);
            m4 = _mm256_mul_ps (m1, clL[C]);
            m5 = _mm256_mul_ps (m2, clR[C]);
            m6 = _mm256_mul_ps (m3, clA[C]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[CG]);
            m2 = _mm256_broadcast_ss (&tiPR[CG]);
            m3 = _mm256_broadcast_ss (&tiPA[CG]);
            m4 = _mm256_mul_ps (m1, clL[G]);
            m5 = _mm256_mul_ps (m2, clR[G]);
            m6 = _mm256_mul_ps (m3, clA[G]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[CT]);
            m2 = _mm256_broadcast_ss (&tiPR[CT]);
            m3 = _mm256_broadcast_ss (&tiPA[CT]);
            m4 = _mm256_mul_ps (m1, clL[T]);
            m5 = _mm256_mul_ps (m2, clR[T]);
            m6 = _mm256_mul_ps (m3, clA[T]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m7 = _mm256_mul_ps (m7, m8);
            *clP++ = _mm256_mul_ps (m7, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[GA]);
            m2 = _mm256_broadcast_ss (&tiPR[GA]);
            m3 = _mm256_broadcast_ss (&tiPA[GA]);
            m7 = _mm256_mul_ps (m1, clL[A]);
            m8 = _mm256_mul_ps (m2, clR[A]);
            m9 = _mm256_mul_ps (m3, clA[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[GC]);
            m2 = _mm256_broadcast_ss (&tiPR[GC]);
            m3 = _mm256_broadcast_ss (&tiPA[GC]);
            m4 = _mm256_mul_ps (m1, clL[C]);
            m5 = _mm256_mul_ps (m2, clR[C]);
            m6 = _mm256_mul_ps (m3, clA[C]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[GG]);
            m2 = _mm256_broadcast_ss (&tiPR[GG]);
            m3 = _mm256_broadcast_ss (&tiPA[GG]);
            m4 = _mm256_mul_ps (m1, clL[G]);
            m5 = _mm256_mul_ps (m2, clR[G]);
            m6 = _mm256_mul_ps (m3, clA[G]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[GT]);
            m2 = _mm256_broadcast_ss (&tiPR[GT]);
            m3 = _mm256_broadcast_ss (&tiPA[GT]);
            m4 = _mm256_mul_ps (m1, clL[T]);
            m5 = _mm256_mul_ps (m2, clR[T]);
            m6 = _mm256_mul_ps (m3, clA[T]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m7 = _mm256_mul_ps (m7, m8);
            *clP++ = _mm256_mul_ps (m7, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[TA]);
            m2 = _mm256_broadcast_ss (&tiPR[TA]);
            m3 = _mm256_broadcast_ss (&tiPA[TA]);
            m7 = _mm256_mul_ps (m1, clL[A]);
            m8 = _mm256_mul_ps (m2, clR[A]);
            m9 = _mm256_mul_ps (m3, clA[A]);
            
            m1 = _mm256_broadcast_ss (&tiPL[TC]);
            m2 = _mm256_broadcast_ss (&tiPR[TC]);
            m3 = _mm256_broadcast_ss (&tiPA[TC]);
            m4 = _mm256_mul_ps (m1, clL[C]);
            m5 = _mm256_mul_ps (m2, clR[C]);
            m6 = _mm256_mul_ps (m3, clA[C]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[TG]);
            m2 = _mm256_broadcast_ss (&tiPR[TG]);
            m3 = _mm256_broadcast_ss (&tiPA[TG]);
            m4 = _mm256_mul_ps (m1, clL[G]);
            m5 = _mm256_mul_ps (m2, clR[G]);
            m6 = _mm256_mul_ps (m3, clA[G]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m1 = _mm256_broadcast_ss (&tiPL[TT]);
            m2 = _mm256_broadcast_ss (&tiPR[TT]);
            m3 = _mm256_broadcast_ss (&tiPA[TT]);
            m4 = _mm256_mul_ps (m1, clL[T]);
            m5 = _mm256_mul_ps (m2, clR[T]);
            m6 = _mm256_mul_ps (m3, clA[T]);
            m7 = _mm256_add_ps (m4, m7);
            m8 = _mm256_add_ps (m5, m8);
            m9 = _mm256_add_ps (m6, m9);
            
            m7 = _mm256_mul_ps (m7, m8);
            *clP++ = _mm256_mul_ps (m7, m9);
            
            clL += 4;
            clR += 4;
            clA += 4;
        }
        tiPL += 16;
        tiPR += 16;
        tiPA += 16;
    }
    
    return NO_ERROR;
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeRoot_NUC4_SSE: 4by4 nucleotide model with or without rate
|       variation using SSE instructions
|
-----------------------------------------------------------------*/
int CondLikeRoot_NUC4_SSE (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *pL, *pR, *pA, *tiPL, *tiPR, *tiPA;
    __m128          *clL, *clR, *clP, *clA;
    __m128          m1, m2, m3, m4, m5, m6, m7, m8, m9;
    ModelInfo       *m;

    m = &modelSettings[division];

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    tiPL = pL;
    tiPR = pR;
    tiPA = pA;
    for (k=0; k<m->numRateCats; k++)
        {
        for (c=0; c<m->numVecChars; c++)
            {
            m1 = _mm_load1_ps (&tiPL[AA]);
            m2 = _mm_load1_ps (&tiPR[AA]);
            m3 = _mm_load1_ps (&tiPA[AA]);
            m7 = _mm_mul_ps (m1, clL[A]);
            m8 = _mm_mul_ps (m2, clR[A]);
            m9 = _mm_mul_ps (m3, clA[A]);

            m1 = _mm_load1_ps (&tiPL[AC]);
            m2 = _mm_load1_ps (&tiPR[AC]);
            m3 = _mm_load1_ps (&tiPA[AC]);
            m4 = _mm_mul_ps (m1, clL[C]);
            m5 = _mm_mul_ps (m2, clR[C]);
            m6 = _mm_mul_ps (m3, clA[C]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m1 = _mm_load1_ps (&tiPL[AG]);
            m2 = _mm_load1_ps (&tiPR[AG]);
            m3 = _mm_load1_ps (&tiPA[AG]);
            m4 = _mm_mul_ps (m1, clL[G]);
            m5 = _mm_mul_ps (m2, clR[G]);
            m6 = _mm_mul_ps (m3, clA[G]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m1 = _mm_load1_ps (&tiPL[AT]);
            m2 = _mm_load1_ps (&tiPR[AT]);
            m3 = _mm_load1_ps (&tiPA[AT]);
            m4 = _mm_mul_ps (m1, clL[T]);
            m5 = _mm_mul_ps (m2, clR[T]);
            m6 = _mm_mul_ps (m3, clA[T]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m7 = _mm_mul_ps (m7, m8);
            *clP++ = _mm_mul_ps (m7, m9);

            m1 = _mm_load1_ps (&tiPL[CA]);
            m2 = _mm_load1_ps (&tiPR[CA]);
            m3 = _mm_load1_ps (&tiPA[CA]);
            m7 = _mm_mul_ps (m1, clL[A]);
            m8 = _mm_mul_ps (m2, clR[A]);
            m9 = _mm_mul_ps (m3, clA[A]);

            m1 = _mm_load1_ps (&tiPL[CC]);
            m2 = _mm_load1_ps (&tiPR[CC]);
            m3 = _mm_load1_ps (&tiPA[CC]);
            m4 = _mm_mul_ps (m1, clL[C]);
            m5 = _mm_mul_ps (m2, clR[C]);
            m6 = _mm_mul_ps (m3, clA[C]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m1 = _mm_load1_ps (&tiPL[CG]);
            m2 = _mm_load1_ps (&tiPR[CG]);
            m3 = _mm_load1_ps (&tiPA[CG]);
            m4 = _mm_mul_ps (m1, clL[G]);
            m5 = _mm_mul_ps (m2, clR[G]);
            m6 = _mm_mul_ps (m3, clA[G]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m1 = _mm_load1_ps (&tiPL[CT]);
            m2 = _mm_load1_ps (&tiPR[CT]);
            m3 = _mm_load1_ps (&tiPA[CT]);
            m4 = _mm_mul_ps (m1, clL[T]);
            m5 = _mm_mul_ps (m2, clR[T]);
            m6 = _mm_mul_ps (m3, clA[T]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m7 = _mm_mul_ps (m7, m8);
            *clP++ = _mm_mul_ps (m7, m9);

            m1 = _mm_load1_ps (&tiPL[GA]);
            m2 = _mm_load1_ps (&tiPR[GA]);
            m3 = _mm_load1_ps (&tiPA[GA]);
            m7 = _mm_mul_ps (m1, clL[A]);
            m8 = _mm_mul_ps (m2, clR[A]);
            m9 = _mm_mul_ps (m3, clA[A]);

            m1 = _mm_load1_ps (&tiPL[GC]);
            m2 = _mm_load1_ps (&tiPR[GC]);
            m3 = _mm_load1_ps (&tiPA[GC]);
            m4 = _mm_mul_ps (m1, clL[C]);
            m5 = _mm_mul_ps (m2, clR[C]);
            m6 = _mm_mul_ps (m3, clA[C]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m1 = _mm_load1_ps (&tiPL[GG]);
            m2 = _mm_load1_ps (&tiPR[GG]);
            m3 = _mm_load1_ps (&tiPA[GG]);
            m4 = _mm_mul_ps (m1, clL[G]);
            m5 = _mm_mul_ps (m2, clR[G]);
            m6 = _mm_mul_ps (m3, clA[G]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m1 = _mm_load1_ps (&tiPL[GT]);
            m2 = _mm_load1_ps (&tiPR[GT]);
            m3 = _mm_load1_ps (&tiPA[GT]);
            m4 = _mm_mul_ps (m1, clL[T]);
            m5 = _mm_mul_ps (m2, clR[T]);
            m6 = _mm_mul_ps (m3, clA[T]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m7 = _mm_mul_ps (m7, m8);
            *clP++ = _mm_mul_ps (m7, m9);

            m1 = _mm_load1_ps (&tiPL[TA]);
            m2 = _mm_load1_ps (&tiPR[TA]);
            m3 = _mm_load1_ps (&tiPA[TA]);
            m7 = _mm_mul_ps (m1, clL[A]);
            m8 = _mm_mul_ps (m2, clR[A]);
            m9 = _mm_mul_ps (m3, clA[A]);

            m1 = _mm_load1_ps (&tiPL[TC]);
            m2 = _mm_load1_ps (&tiPR[TC]);
            m3 = _mm_load1_ps (&tiPA[TC]);
            m4 = _mm_mul_ps (m1, clL[C]);
            m5 = _mm_mul_ps (m2, clR[C]);
            m6 = _mm_mul_ps (m3, clA[C]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m1 = _mm_load1_ps (&tiPL[TG]);
            m2 = _mm_load1_ps (&tiPR[TG]);
            m3 = _mm_load1_ps (&tiPA[TG]);
            m4 = _mm_mul_ps (m1, clL[G]);
            m5 = _mm_mul_ps (m2, clR[G]);
            m6 = _mm_mul_ps (m3, clA[G]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m1 = _mm_load1_ps (&tiPL[TT]);
            m2 = _mm_load1_ps (&tiPR[TT]);
            m3 = _mm_load1_ps (&tiPA[TT]);
            m4 = _mm_mul_ps (m1, clL[T]);
            m5 = _mm_mul_ps (m2, clR[T]);
            m6 = _mm_mul_ps (m3, clA[T]);
            m7 = _mm_add_ps (m4, m7);
            m8 = _mm_add_ps (m5, m8);
            m9 = _mm_add_ps (m6, m9);

            m7 = _mm_mul_ps (m7, m8);
            *clP++ = _mm_mul_ps (m7, m9);

            clL += 4;
            clR += 4;
            clA += 4;
            }
        tiPL += 16;
        tiPR += 16;
        tiPA += 16;
        }

    return NO_ERROR;
}
#endif


#if !defined (SSE_ENABLED) || 1
/*----------------------------------------------------------------
|
|   CondLikeRoot_NY98: codon model with omega variation
|
-----------------------------------------------------------------*/
int CondLikeRoot_NY98 (TreeNode *p, int division, int chain)
{
    int             a, b, c, d, h, i, j, k, shortCut, *lState=NULL, *rState=NULL, *aState=NULL,
                    nStates, nStatesSquared;
    CLFlt           likeL, likeR, likeA, *clL, *clR, *clP, *clA, *pL, *pR, *pA,
                    *tiPL, *tiPR, *tiPA;
    ModelInfo       *m;
    
    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for anc branch, always terminal */
    if (m->isPartAmbig[p->anc->index] == YES)
        {
        shortCut = 4;
        }
    else 
        {
        aState = m->termState[p->anc->index];
        tiPA = pA;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeA[a++] = tiPA[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeA[a++] = 1.0;
            tiPA += nStatesSquared;
            }
        }
#   else
    shortCut = 4;
#   endif

        switch (shortCut)
        {
        case 4:
            tiPL = pL;
            tiPR = pR;
            tiPA = pA;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    for (i=h=0; i<nStates; i++)
                        {
                        likeL = likeR = likeA = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeA += tiPA[h]*clA[j];
                            likeL += tiPL[h]*clL[j];
                            likeR += tiPR[h++]*clR[j];
                            }
                        *(clP++) = likeL * likeR * likeA;
                        }
                    clL += nStates;
                    clR += nStates;
                    clA += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                tiPA += nStatesSquared;
                }
            break;
        case 0:
            tiPR = pR;
            tiPL = pL;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    b = aState[c] + k*(nStatesSquared+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeR = likeL = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeR += tiPR[h]*clR[j];
                            likeL += tiPL[h++]*clL[j];
                            }
                        *(clP++) =  preLikeA[b++] * likeL * likeR;
                        }
                    clR += nStates;
                    clL += nStates;
                    }
                tiPR += nStatesSquared;
                tiPL += nStatesSquared;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = lState[c] + k*(nStatesSquared+nStates);
                    b = aState[c] + k*(nStatesSquared+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeR = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeR += tiPR[h++]*clR[j];
                            }
                        *(clP++) = preLikeL[a++] * preLikeA[b++] * likeR;
                        }
                    clR += nStates;
                    }
                tiPR += nStatesSquared;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = rState[c] + k*(nStatesSquared+nStates);
                    b = aState[c] + k*(nStatesSquared+nStates);
                    for (i=h=0; i<nStates; i++)
                        {
                        likeL = 0.0;
                        for (j=0; j<nStates; j++)
                            {
                            likeL += tiPL[h++]*clL[j];
                            }
                        *(clP++) = preLikeR[a++] * preLikeA[b++] * likeL;
                        }
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                }
            break;
        case 3:
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numChars; c++)
                    {
                    a = rState[c] + k*(nStatesSquared+nStates);
                    b = lState[c] + k*(nStatesSquared+nStates);
                    d = aState[c] + k*(nStatesSquared+nStates);
                    for (i=0; i<nStates; i++)
                        {
                        *(clP++) = preLikeR[a++] * preLikeL[b++] * preLikeA[d++];
                        }
                    }
                }
            break;
        }

    return NO_ERROR;
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeRoot_NY98_SSE: codon model with omega variation
|
-----------------------------------------------------------------*/
int CondLikeRoot_NY98_SSE (TreeNode *p, int division, int chain)
{
    int             c, c1, t, h, i, j, k, shortCut, *lState=NULL, *rState=NULL, *aState=NULL,
                    nStates, nStatesSquared;
    CLFlt           *pL, *pR, *pA,
                    *tiPL, *tiPR, *tiPA;
    __m128          *clL, *clR, *clP, *clA;
    __m128          mTiPL, mTiPR, mTiPA, mL, mR, mA, mAcumL, mAcumR, mAcumA;
    ModelInfo       *m;
    CLFlt           *preLikeRV[4];
    CLFlt           *preLikeLV[4];
    CLFlt           *preLikeAV[4];

#   if !defined (DEBUG_NOSHORTCUTS)
    int             a;

#   endif

    /* find model settings for this division and nStates, nStatesSquared */
    m = &modelSettings[division];
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* find likelihoods of site patterns for left branch if terminal */
    shortCut = 0;
#   if !defined (DEBUG_NOSHORTCUTS)
    if (p->left->left == NULL && m->isPartAmbig[p->left->index] == NO)
        {
        shortCut |= 1;
        lState = m->termState[p->left->index];
        tiPL = pL;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeL[a++] = tiPL[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeL[a++] = 1.0;
            tiPL += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for right branch if terminal */
    if (p->right->left == NULL && m->isPartAmbig[p->right->index] == NO)
        {
        shortCut |= 2;
        rState = m->termState[p->right->index];
        tiPR = pR;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeR[a++] = tiPR[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeR[a++] = 1.0;
            tiPR += nStatesSquared;
            }
        }

    /* find likelihoods of site patterns for anc branch, always terminal */
    if (m->isPartAmbig[p->anc->index] == YES)
        {
        shortCut = 4;
        }
    else 
        {
        aState = m->termState[p->anc->index];
        tiPA = pA;
        for (k=a=0; k<m->numOmegaCats; k++)
            {
            for (i=0; i<nStates; i++)
                for (j=i; j<nStatesSquared; j+=nStates)
                    preLikeA[a++] = tiPA[j];
            /* for ambiguous */
            for (i=0; i<nStates; i++)
                preLikeA[a++] = 1.0;
            tiPA += nStatesSquared;
            }
        }
#   else
    shortCut = 4;
#   endif
        switch (shortCut)
        {
        case 4:
            tiPL = pL;
            tiPR = pR;
            tiPA = pA;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=0; c<m->numVecChars; c++)
                    {
                    for (i=h=0; i<nStates; i++)
                        {
                        mAcumL = _mm_setzero_ps();
                        mAcumR = _mm_setzero_ps();
                        mAcumA = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h]);
                            mTiPR  = _mm_load1_ps (&tiPR[h]);
                            mTiPA  = _mm_load1_ps (&tiPA[h++]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mA     = _mm_mul_ps (mTiPA, clA[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            mAcumA = _mm_add_ps (mA, mAcumA);
                            }
                        mAcumL = _mm_mul_ps (mAcumL, mAcumR);
                        *(clP++) = _mm_mul_ps (mAcumL, mAcumA);
                        }
                    clL += nStates;
                    clR += nStates;
                    clA += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                tiPA += nStatesSquared;
                }
            break;
        case 0:
            tiPL =pL;
            tiPR =pR;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeAV[c1] = &preLikeA[aState[t] + k*(nStatesSquared+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following statment we assume that SSE register can hold exactly 4 ClFlts. */
                        mAcumA = _mm_set_ps (*(preLikeAV[3]++), *(preLikeAV[2]++), *(preLikeAV[1]++), *(preLikeAV[0]++));
                        mAcumL = _mm_setzero_ps();
                        mAcumR = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            mTiPR  = _mm_load1_ps (&tiPR[h++]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            }
                        mAcumL = _mm_mul_ps (mAcumL, mAcumR);
                        *(clP++) = _mm_mul_ps (mAcumL, mAcumA);
                        }
                    clR += nStates;
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                tiPR += nStatesSquared;
                }
            break;
        case 1:
            tiPR = pR;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeLV[c1] = &preLikeL[lState[t] + k*(nStatesSquared+nStates)];
                        preLikeAV[c1] = &preLikeA[aState[t] + k*(nStatesSquared+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following statment we assume that SSE register can hold exactly 4 ClFlts. */
                        mAcumL = _mm_set_ps (*(preLikeLV[3]++), *(preLikeLV[2]++), *(preLikeLV[1]++), *(preLikeLV[0]++));
                        mAcumA = _mm_set_ps (*(preLikeAV[3]++), *(preLikeAV[2]++), *(preLikeAV[1]++), *(preLikeAV[0]++));
                        mAcumR = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPR  = _mm_load1_ps (&tiPR[h++]);
                            mR     = _mm_mul_ps (mTiPR, clR[j]);
                            mAcumR = _mm_add_ps (mR, mAcumR);
                            }
                        mAcumL = _mm_mul_ps (mAcumL, mAcumR);
                        *(clP++) = _mm_mul_ps (mAcumL, mAcumA);
                        }
                    clR += nStates;
                    }
                tiPR += nStatesSquared;
                }
            break;
        case 2:
            tiPL = pL;
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeRV[c1] = &preLikeR[rState[t] + k*(nStatesSquared+nStates)];
                        preLikeAV[c1] = &preLikeA[aState[t] + k*(nStatesSquared+nStates)];
                        }
                    for (i=h=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following statment we assume that SSE register can hold exactly 4 ClFlts. */
                        mAcumR = _mm_set_ps (*(preLikeRV[3]++), *(preLikeRV[2]++), *(preLikeRV[1]++), *(preLikeRV[0]++));
                        mAcumA = _mm_set_ps (*(preLikeAV[3]++), *(preLikeAV[2]++), *(preLikeAV[1]++), *(preLikeAV[0]++));
                        mAcumL = _mm_setzero_ps();
                        for (j=0; j<nStates; j++)
                            {
                            mTiPL  = _mm_load1_ps (&tiPL[h++]);
                            mL     = _mm_mul_ps (mTiPL, clL[j]);
                            mAcumL = _mm_add_ps (mL, mAcumL);
                            }
                        mAcumL = _mm_mul_ps (mAcumL, mAcumR);
                        *(clP++) = _mm_mul_ps (mAcumL,mAcumA);
                        }
                    clL += nStates;
                    }
                tiPL += nStatesSquared;
                }
            break;
        case 3:
            for (k=0; k<m->numOmegaCats; k++)
                {
                for (c=t=0; c<m->numVecChars; c++)
                    {
                    for (c1=0; c1<m->numFloatsPerVec; c1++,t++)
                        {
                        preLikeRV[c1] = &preLikeR[rState[t] + k*(nStatesSquared+nStates)];
                        preLikeLV[c1] = &preLikeL[lState[t] + k*(nStatesSquared+nStates)];
                        preLikeAV[c1] = &preLikeA[aState[t] + k*(nStatesSquared+nStates)];
                        }
                    for (i=0; i<nStates; i++)
                        {
                        assert (m->numFloatsPerVec == 4); /* In the following 2 statments we assume that SSE register can hold exactly 4 ClFlts. */
                        mL = _mm_set_ps (*(preLikeLV[3]++), *(preLikeLV[2]++), *(preLikeLV[1]++), *(preLikeLV[0]++));
                        mR = _mm_set_ps (*(preLikeRV[3]++), *(preLikeRV[2]++), *(preLikeRV[1]++), *(preLikeRV[0]++));
                        mA = _mm_set_ps (*(preLikeAV[3]++), *(preLikeAV[2]++), *(preLikeAV[1]++), *(preLikeAV[0]++));
                        mL = _mm_mul_ps (mL,mR);
                        *(clP++) = _mm_mul_ps (mL,mA);
                        }
                    }
                }
            break;
        }

    return NO_ERROR;
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeRoot_Std: variable number of states model
|       with or without rate variation
|
-----------------------------------------------------------------*/
int CondLikeRoot_Std (TreeNode *p, int division, int chain)
{
    int             a, c, h, i, j, k, nStates=0, nCats=0, tmp;
    CLFlt           *clL, *clR, *clP, *clA, *pL, *pR, *pA, *tiPL, *tiPR, *tiPA,
                    likeL, likeR, likeA;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* flip state of node so that we are not overwriting old cond likes */
    FlipCondLikeSpace (m, chain, p->index);
    
    /* find conditional likelihood pointers */
    clL = m->condLikes[m->condLikeIndex[chain][p->left->index ]];
    clR = m->condLikes[m->condLikeIndex[chain][p->right->index]];
    clP = m->condLikes[m->condLikeIndex[chain][p->index       ]];
    clA = m->condLikes[m->condLikeIndex[chain][p->anc->index  ]];

    /* find transition probabilities (or calculate instead) */
    pL = m->tiProbs[m->tiProbsIndex[chain][p->left->index ]];
    pR = m->tiProbs[m->tiProbsIndex[chain][p->right->index]];
    pA = m->tiProbs[m->tiProbsIndex[chain][p->index       ]];

    /* calculate ancestral probabilities */
    for (k=h=0; k<m->numRateCats; k++)
        {
        /* calculate ancestral probabilities */
        for (c=0; c<m->numChars; c++)
            {
            nStates = m->nStates[c];
        
            /* the following lines ensure that nCats is 1 unless */
            /* the character is binary and beta categories are used  */
            if (nStates == 2)
                nCats = m->numBetaCats;
            else
                nCats = 1;

            tmp = k*nStates*nStates; /* tmp contains offset to skip gamma cats that already processed*/
            tiPL = pL + m->tiIndex[c] + tmp;
            tiPR = pR + m->tiIndex[c] + tmp;
            tiPA = pA + m->tiIndex[c] + tmp;
            tmp = (m->numRateCats-1)*2*2; /* tmp contains size of block of tpi matrices across all rate cats (minus one) for single beta category. Further used only if character is binary to jump to next beta category */
                
            for (j=0; j<nCats;j++)
                {
                for (a=0; a<nStates; a++)
                    {
                    likeL = likeR = likeA = 0.0;
                    for (i=0; i<nStates; i++)
                        {
                        likeL += *(tiPL++) * clL[i];
                        likeR += *(tiPR++) * clR[i];
                        likeA += *(tiPA++) * clA[i];
                        }
                    clP[h++] = likeL * likeR * likeA;
                    }
                clL += nStates;
                clR += nStates;
                clA += nStates;
        
                tiPL += tmp;
                tiPR += tmp;
                tiPA += tmp;
                }
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeUp_Bin: pull likelihoods up and calculate scaled
|       finals, binary model with or without rate variation
|
-----------------------------------------------------------------*/
int CondLikeUp_Bin (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *clFA, *clFP, *clDP, *tiP, condLikeUp[2], sum[2];
    ModelInfo       *m;
    
    /* find model settings for this division */
    m = &modelSettings[division];

    if (p->anc->anc == NULL)
        {
        /* this is the root node */
        /* find conditional likelihood pointers = down cond likes */
        /* use conditional likelihood scratch space for final cond likes */
        clDP = m->condLikes[m->condLikeIndex[chain][p->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index]];

        for (k=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                *(clFP++) = *(clDP++);
                *(clFP++) = *(clDP++);
                }
            }
        }
    else
        {
        /* find conditional likelihood pointers */
        /* use conditional likelihood scratch space for final cond likes */
        clFA = m->condLikes[m->condLikeScratchIndex[p->anc->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index     ]];
        clDP = m->condLikes[m->condLikeIndex[chain][p->index     ]];
        
        /* find transition probabilities */
        tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];

        for (k=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                condLikeUp[0] = condLikeUp[1] = 0.0;
                
                sum[0] = tiP[0]*clDP[0] + tiP[1]*clDP[1];
                sum[1] = tiP[2]*clDP[0] + tiP[3]*clDP[1];

                if (sum[0] != 0.0) condLikeUp[0] = clFA[0] / sum[0];
                if (sum[1] != 0.0) condLikeUp[1] = clFA[1] / sum[1];
                
                *(clFP++) = (condLikeUp[0]*tiP[0] + condLikeUp[1]*tiP[1])*clDP[0];
                *(clFP++) = (condLikeUp[0]*tiP[2] + condLikeUp[1]*tiP[3])*clDP[1];
                
                clFA += 2;
                clDP += 2;
                }
            tiP += 4;
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeUp_Gen: pull likelihoods up and calculate scaled
|       finals for an interior node
|
-----------------------------------------------------------------*/
int CondLikeUp_Gen (TreeNode *p, int division, int chain)
{
    int             a, c, i, j, k, nStates, nStatesSquared, nRateCats;
    CLFlt           *clFA, *clFP, *clDP, *tiP, *condLikeUp, sum;
    ModelInfo       *m;
    
    /* find model settings for this division */
    m = &modelSettings[division];

    /* find number of states in the model */
    nStates = m->numModelStates;
    nStatesSquared = nStates * nStates;

    /* find number of gamma cats */
    nRateCats = m->numRateCats;
    if (m->gibbsGamma == YES)
        nRateCats = 1;

    /* use preallocated scratch space */
    condLikeUp = m->ancStateCondLikes;

    /* calculate final states */
    if (p->anc->anc == NULL)
        {
        /* this is the root node */
        /* find conditional likelihood pointers = down cond likes */
        /* use conditional likelihood scratch space for final cond likes */
        clDP = m->condLikes[m->condLikeIndex[chain][p->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index]];
        
        /* final cond likes = downpass cond likes */
        for (k=0; k<nRateCats; k++)
            {
            /* copy cond likes */ 
            for (c=0; c<m->numChars*nStates; c++)
                *(clFP++) = *(clDP++);
            }
        }
    else
        {
        /* find conditional likelihood pointers */
        /* use conditional likelihood scratch space for final cond likes */
        clFA = m->condLikes[m->condLikeScratchIndex[p->anc->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index     ]];
        clDP = m->condLikes[m->condLikeIndex[chain][p->index     ]];
        
        /* find transition probabilities */
        tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];
        
        for (k=0; k<nRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                for (a=j=0; a<nStates; a++)
                    {
                    sum = 0.0;
                    for (i=0; i<nStates; i++)
                        sum += tiP[j++]*clDP[i];
                    if (sum != 0.0) condLikeUp[a] = clFA[a] / sum;
                    }
                    
                for (a=j=0; a<nStates; a++)
                    {
                    sum = 0.0;
                    for (i=0; i<nStates; i++)
                        {
                        sum += condLikeUp[i] * tiP[j++];
                        }
                    *(clFP++) = sum * clDP[a];
                    }

                clFA += nStates;
                clDP += nStates;
                }
            tiP += nStatesSquared;
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeUp_NUC4: pull likelihoods up and calculate scaled
|       finals for an interior node
|
-----------------------------------------------------------------*/
int     CondLikeUp_NUC4 (TreeNode *p, int division, int chain)
{
    int             c, k, nRateCats;
    CLFlt           *clFA, *clFP, *clDP, *tiP, condLikeUp[4], sum[4];
    ModelInfo       *m;
    
    /* find model settings for this division */
    m = &modelSettings[division];

    /* find number of rate cats */
    nRateCats = m->numRateCats;
    if (m->gibbsGamma == YES)
        nRateCats = 1;

    /* calculate final states */
    if (p->anc->anc == NULL)
        {
        /* this is the root node */
        /* find conditional likelihood pointers = down cond likes */
        /* use conditional likelihood scratch space for final cond likes */
        clDP = m->condLikes[m->condLikeIndex[chain][p->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index]];
        
        /* final cond likes = downpass cond likes */
        for (k=0; k<nRateCats; k++)
            {
            /* copy cond likes */ 
            for (c=0; c<m->numChars; c++)
                {
                *(clFP++) = *(clDP++);
                *(clFP++) = *(clDP++);
                *(clFP++) = *(clDP++);
                *(clFP++) = *(clDP++);
                }
            }
        }
    else
        {
        /* find conditional likelihood pointers */
        /* use conditional likelihood scratch space for final cond likes */
        clFA = m->condLikes[m->condLikeScratchIndex[p->anc->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index     ]];
        clDP = m->condLikes[m->condLikeIndex[chain][p->index     ]];
        
        /* find transition probabilities */
        tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];
        
        for (k=0; k<nRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {
                condLikeUp[A] = condLikeUp[C] = condLikeUp[G] = condLikeUp[T] = 0.0;

                sum[A] = (tiP[AA]*clDP[A] + tiP[AC]*clDP[C] + tiP[AG]*clDP[G] + tiP[AT]*clDP[T]);
                sum[C] = (tiP[CA]*clDP[A] + tiP[CC]*clDP[C] + tiP[CG]*clDP[G] + tiP[CT]*clDP[T]);
                sum[G] = (tiP[GA]*clDP[A] + tiP[GC]*clDP[C] + tiP[GG]*clDP[G] + tiP[GT]*clDP[T]);
                sum[T] = (tiP[TA]*clDP[A] + tiP[TC]*clDP[C] + tiP[TG]*clDP[G] + tiP[TT]*clDP[T]);

                if (sum[A] != 0.0) condLikeUp[A] = clFA[A] / sum[A];
                if (sum[C] != 0.0) condLikeUp[C] = clFA[C] / sum[C];
                if (sum[G] != 0.0) condLikeUp[G] = clFA[G] / sum[G];
                if (sum[T] != 0.0) condLikeUp[T] = clFA[T] / sum[T];

/*
                clFP[A] = (condLikeUp[A]*tiP[AA] + condLikeUp[C]*tiP[CA] + condLikeUp[G]*tiP[GA] + condLikeUp[T]*tiP[TA])*clDP[A];
                clFP[C] = (condLikeUp[A]*tiP[AC] + condLikeUp[C]*tiP[CC] + condLikeUp[G]*tiP[GC] + condLikeUp[T]*tiP[TC])*clDP[C];
                clFP[G] = (condLikeUp[A]*tiP[AG] + condLikeUp[C]*tiP[CG] + condLikeUp[G]*tiP[GG] + condLikeUp[T]*tiP[TG])*clDP[G];
                clFP[T] = (condLikeUp[A]*tiP[AT] + condLikeUp[C]*tiP[CT] + condLikeUp[G]*tiP[GT] + condLikeUp[T]*tiP[TT])*clDP[T];
*/

                clFP[A] = (condLikeUp[A]*tiP[AA] + condLikeUp[C]*tiP[AC] + condLikeUp[G]*tiP[AG] + condLikeUp[T]*tiP[AT])*clDP[A];
                clFP[C] = (condLikeUp[A]*tiP[CA] + condLikeUp[C]*tiP[CC] + condLikeUp[G]*tiP[CG] + condLikeUp[T]*tiP[CT])*clDP[C];
                clFP[G] = (condLikeUp[A]*tiP[GA] + condLikeUp[C]*tiP[GC] + condLikeUp[G]*tiP[GG] + condLikeUp[T]*tiP[GT])*clDP[G];
                clFP[T] = (condLikeUp[A]*tiP[TA] + condLikeUp[C]*tiP[TC] + condLikeUp[G]*tiP[TG] + condLikeUp[T]*tiP[TT])*clDP[T];

                clFA += 4;
                clFP += 4;
                clDP += 4;
                }
            tiP += 16;
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeUp_Std: pull likelihoods up and calculate scaled
|       finals for an interior node
|
-----------------------------------------------------------------*/
int     CondLikeUp_Std (TreeNode *p, int division, int chain)
{
    int             a, c, i, j, k, t, nStates, nCats, coppySize,tmp;
    CLFlt           *clFA, *clFP, *clDP, *pA, *tiP, condLikeUp[10], sum;
    ModelInfo       *m;
    
    /* find model settings for this division */
    m = &modelSettings[division];
    
    /* calculate final states */
    if (p->anc->anc == NULL)
        {
        /* this is the root node */
        /* find conditional likelihood pointers = down cond likes */
        /* use conditional likelihood scratch space for final cond likes */
        clDP = m->condLikes[m->condLikeIndex[chain][p->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index]];
        
        coppySize=0;
        /* final cond likes = downpass cond likes */
        for (c=0; c<m->numChars; c++)
            {
            /* calculate nStates and nCats */
            nStates = m->nStates[c];
            
            /* the following lines ensure that nCats is 1 unless */
            /* the character is binary and beta categories are used  */
            if (nStates == 2)
                nCats = m->numBetaCats;
            else
                nCats = 1;

            coppySize+=nCats*nStates;
            }

        /* finally multiply with the rate cats */
        coppySize *= m->numRateCats;

        /* copy cond likes */ 
        for (k=0; k<coppySize; k++)
            *(clFP++) = *(clDP++);
        }
    else
        {
        /* find conditional likelihood pointers */
        /* use conditional likelihood scratch space for final cond likes */
        clFA = m->condLikes[m->condLikeScratchIndex[p->anc->index]];
        clFP = m->condLikes[m->condLikeScratchIndex[p->index     ]];
        clDP = m->condLikes[m->condLikeIndex[chain][p->index     ]];

        /* find transition probabilities */
        pA = m->tiProbs[m->tiProbsIndex[chain][p->index]];
        
        for (k=0; k<m->numRateCats; k++)
            {
            for (c=0; c<m->numChars; c++)
                {

                /* calculate nStates and nCats */
                nStates = m->nStates[c];
                
                /* the following lines ensure that nCats is 1 unless */
                /* the character is binary and beta categories are used  */
                if (nStates == 2)
                    nCats = m->numBetaCats;
                else
                    nCats = 1;

                tmp = k*nStates*nStates; /* tmp contains offset to skip rate cats that already processed*/
                tiP = pA + m->tiIndex[c] + tmp;
                tmp = (m->numRateCats-1)*2*2; /* tmp contains size of block of tpi matrices across all rate cats (minus one) for single beta category. Further used only if character is binary to jump to next beta category */

                /* now calculate the final cond likes */
                for (t=0; t<nCats; t++)
                    {
                    for (a=j=0; a<nStates; a++)
                        {
                        sum = 0.0;
                        for (i=0; i<nStates; i++)
                            sum += tiP[j++]*clDP[i];
                        if (sum == 0.0)
                            condLikeUp[a] = 0.0;    /* we lost the conditional likelihood in the downpass (can occur in gamma model) */
                        else
                            condLikeUp[a] = clFA[a] / sum;
                        }
                        
                    for (a=j=0; a<nStates; a++)
                        {
                        sum = 0.0;
                        for (i=0; i<nStates; i++)
                            {
                            sum += condLikeUp[i] * tiP[j++];
                            }
                        clFP[a] = sum * clDP[a];
                        }

                    clFP += nStates;
                    clFA += nStates;
                    clDP += nStates;
                    tiP += tmp;
                    }
                }
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   CondLikeScaler_Gen: general n-state model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeScaler_Gen (TreeNode *p, int division, int chain)
{
    int             c, k, n, nStates;
    CLFlt           scaler, **clP, *clPtr, *scP, *lnScaler;
    ModelInfo       *m;

    m = &modelSettings[division];
    nStates = m->numModelStates;

    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numChars * m->numModelStates;
        }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* rescale */
    for (c=0; c<m->numChars; c++)
        {
        scaler = 0.0;
        for (k=0; k<m->numRateCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                if (clP[k][n] > scaler)
                    scaler = clP[k][n];
                }
            }

        for (k=0; k<m->numRateCats; k++)
            {
            for (n=0; n<nStates; n++)
                clP[k][n] /= scaler;
            clP[k] += n;
            }

        scP[c]       = (CLFlt) log (scaler);    /* store node scaler */
        lnScaler[c] += scP[c];  /* add into tree scaler  */
        }

    m->unscaledNodes[chain][p->index] = 0;

    return (NO_ERROR);
}


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeScaler_Gen_SSE: general n-state model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeScaler_Gen_SSE (TreeNode *p, int division, int chain)
{
    int             c, k, n, nStates;
    CLFlt           *scP, *lnScaler;
    __m128          *clPtr, **clP, m1;
    ModelInfo       *m;

    m = &modelSettings[division];
    nStates = m->numModelStates;

    /* find conditional likelihood pointers */
    clPtr = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP_SSE;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
        }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];
    //scP_SSE = (__m128 *) scP;

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* rescale */
    for (c=0; c<m->numVecChars; c++)
        {
        //scaler = 0.0;
        m1 = _mm_setzero_ps ();
        for (k=0; k<m->numRateCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                m1 = _mm_max_ps (m1, clP[k][n]);
                }
            }
        _mm_store_ps (scP,  m1);
        scP += m->numFloatsPerVec;

        for (k=0; k<m->numRateCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                *clP[k] = _mm_div_ps (*clP[k], m1);
                clP[k]++;
                }
            }
        }
    
    /* Reset scP to original position*/
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];
    for (c=0; c<m->numChars; c++)
        {
        scP[c]       = (CLFlt) log (scP[c]);    /* store node scaler */
        lnScaler[c] += scP[c];                  /* add into tree scaler  */
        }

    m->unscaledNodes[chain][p->index] = 0;

    return (NO_ERROR);
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeScaler_Gen_GibbsGamma: general n-state model with Gibbs
|       sampling of rate categories in discrete gamma
|
-----------------------------------------------------------------*/
int CondLikeScaler_Gen_GibbsGamma (TreeNode *p, int division, int chain)
{
    int             c, i, j, n, nStates, *rateCat, nRateCats;
    CLFlt           scaler, *clP, *scP, *lnScaler;
    ModelInfo       *m;

    m = &modelSettings[division];
    nStates = m->numModelStates;

    /* find conditional likelihood pointer */
    clP = m->condLikes[m->condLikeIndex[chain][p->index]];

    /* flip node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* find rate category index and number of rate categories */
    rateCat = m->tiIndex + chain * m->numChars;
    nRateCats = m->numRateCats;

    /* scale */
    i = j = 0;
    for (c=0; c<m->numChars; c++)
        {
        if (rateCat[c] < nRateCats)
            {
            scaler = 0.0;
            for (n=0; n<nStates; n++)
                {
                if (clP[i] > scaler)
                    scaler = clP[i];
                i++;
                }


            for (n=0; n<nStates; n++)
                clP[j++] /= scaler;

            scP[c]       = (CLFlt) log (scaler);    /* store node scaler */
            lnScaler[c] += scP[c];                  /* add into tree scaler  */

            }
        else
            {
            scP[c] = 0.0;
            /* no need to add it to the lnScaler */
            i += nStates;
            j += nStates;
            }
        }

    m->unscaledNodes[chain][p->index] = 0;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   CondLikeScaler_NUC4: 4by4 nucleotide model with or without rate
|       variation
|
-----------------------------------------------------------------*/
int CondLikeScaler_NUC4 (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           scaler, *scP, *lnScaler, *clPtr, **clP;
    ModelInfo       *m;

    m = &modelSettings[division];

    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numChars * m->numModelStates;
        }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* rescale values */
    for (c=0; c<m->numChars; c++)
        {
        scaler = 0.0;
        for (k=0; k<m->numRateCats; k++)
            {
            if (clP[k][A] > scaler)
                scaler = clP[k][A];
            if (clP[k][C] > scaler)
                scaler = clP[k][C];
            if (clP[k][G] > scaler)
                scaler = clP[k][G];
            if (clP[k][T] > scaler)
                scaler = clP[k][T];
            }

        for (k=0; k<m->numRateCats; k++)
            {
            clP[k][A] /= scaler;
            clP[k][C] /= scaler;
            clP[k][G] /= scaler;
            clP[k][T] /= scaler;
            clP[k] += 4;
            }

        scP[c]       = (CLFlt) log(scaler); /* store node scaler */
        lnScaler[c] += scP[c];  /* add into tree scaler  */
        }

    m->unscaledNodes[chain][p->index] = 0;   /* set unscaled nodes to 0 */

    return NO_ERROR;
}


#if defined (AVX_ENABLED)
/*----------------------------------------------------------------
 |
 |   CondLikeScaler_NUC4_AVX: 4by4 nucleotide model with or without rate
 |       variation using AVX (or AVX + FMA) code
 |
 -----------------------------------------------------------------*/
int CondLikeScaler_NUC4_AVX (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *scP, *lnScaler;
    __m256          *clPtr, **clP, *scP_AVX, m1;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    
    /* find conditional likelihood pointers */
    clPtr = (__m256 *) m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP_AVX;
    for (k=0; k<m->numRateCats; k++)
    {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
    }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];
    scP_AVX = (__m256 *) scP;
    
    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* rescale */
    for (c=0; c<m->numVecChars; c++)
    {
        m1 = _mm256_setzero_ps ();

        for (k=0; k<m->numRateCats; k++)
        {
            m1 = _mm256_max_ps (m1, clP[k][A]);
            m1 = _mm256_max_ps (m1, clP[k][C]);
            m1 = _mm256_max_ps (m1, clP[k][G]);
            m1 = _mm256_max_ps (m1, clP[k][T]);
        }
        
        for (k=0; k<m->numRateCats; k++)
        {
            *clP[k] = _mm256_div_ps (*clP[k], m1);
            clP[k]++;
            *clP[k] = _mm256_div_ps (*clP[k], m1);
            clP[k]++;
            *clP[k] = _mm256_div_ps (*clP[k], m1);
            clP[k]++;
            *clP[k] = _mm256_div_ps (*clP[k], m1);
            clP[k]++;
        }
        
        (*scP_AVX++) = m1;
    }

    /* update site scalers */
    for (c=0; c<m->numChars; c++)
        lnScaler[c] += (scP[c] = logf (scP[c]));    /* add log of new scaler into tree scaler  */
    
    m->unscaledNodes[chain][p->index] = 0;   /* set unscaled nodes to 0 */

    return NO_ERROR;

}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeScaler_NUC4_SSE: 4by4 nucleotide model with or without rate
|       variation using SSE code
|
-----------------------------------------------------------------*/
int CondLikeScaler_NUC4_SSE (TreeNode *p, int division, int chain)
{
    int             c, k;
    CLFlt           *scP, *lnScaler;
    __m128          *clPtr, **clP, *scP_SSE, m1;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* find conditional likelihood pointers */
    clPtr = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP_SSE;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
        }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];
    scP_SSE = (__m128 *) scP;

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* rescale */
    for (c=0; c<m->numVecChars; c++)
        {
        m1 = _mm_setzero_ps ();
        for (k=0; k<m->numRateCats; k++)
            {
            m1 = _mm_max_ps (m1, clP[k][A]);
            m1 = _mm_max_ps (m1, clP[k][C]);
            m1 = _mm_max_ps (m1, clP[k][G]);
            m1 = _mm_max_ps (m1, clP[k][T]);
            }

        for (k=0; k<m->numRateCats; k++)
            {
            *clP[k] = _mm_div_ps (*clP[k], m1);
            clP[k]++;
            *clP[k] = _mm_div_ps (*clP[k], m1);
            clP[k]++;
            *clP[k] = _mm_div_ps (*clP[k], m1);
            clP[k]++;
            *clP[k] = _mm_div_ps (*clP[k], m1);
            clP[k]++;
            }

        (*scP_SSE++) = m1;
        }

    /* update site scalers */
    for (c=0; c<m->numChars; c++)
        lnScaler[c] += (scP[c] = (CLFlt)(log (scP[c])));    /* add log of new scaler into tree scaler  */

    m->unscaledNodes[chain][p->index] = 0;   /* number of unscaled nodes is 0 */

    return NO_ERROR;
    
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeScaler_NUC4_GibbsGamma: 4by4 nucleotide model with rate
|       variation approximated by Gibbs sampling from gamma
|
-----------------------------------------------------------------*/
int CondLikeScaler_NUC4_GibbsGamma (TreeNode *p, int division, int chain)
{
    int             c, i, j, nRateCats, *rateCat;
    CLFlt           scaler, *clP, *scP, *lnScaler;
    ModelInfo       *m;

    m = &modelSettings[division];

    /* find conditional likelihood pointer */
    clP = m->condLikes[m->condLikeIndex[chain][p->index]];

    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* find rate category index and number of gamma categories */
    rateCat = m->tiIndex + chain * m->numChars;
    nRateCats = m->numRateCats;

    /* scale */
    i = j = 0;
    for (c=0; c<m->numChars; c++)
        {
        if (rateCat[c] < nRateCats)
            {
            scaler = 0.0;
            if (clP[i] > scaler)
                scaler = clP[i];
            i++;
            if (clP[i] > scaler)
                scaler = clP[i];
            i++;
            if (clP[i] > scaler)
                scaler = clP[i];
            i++;
            if (clP[i] > scaler)
                scaler = clP[i];
            i++;

            clP[j++] /= scaler;
            clP[j++] /= scaler;
            clP[j++] /= scaler;
            clP[j++] /= scaler;

            scP[c]       = (CLFlt) log (scaler);    /* store node scaler */
            lnScaler[c] += scP[c];                  /* add into tree scaler  */
            }
        else
            {
            scP[c] = 0.0;   /* store node scaler */
            /* no need to add it to the lnScaler */
            i += 4;
            j += 4;
            }
        }

    m->unscaledNodes[chain][p->index] = 0;

    return NO_ERROR;
}


#if !defined (SSE_ENABLED) || 1
/*----------------------------------------------------------------
|
|   CondLikeScaler_NY98: codon model with omega variation
|
-----------------------------------------------------------------*/
int CondLikeScaler_NY98 (TreeNode *p, int division, int chain)
{
    int             c, k, n, nStates;
    CLFlt           scaler, **clP, *clPtr, *scP, *lnScaler;
    ModelInfo       *m;

    m = &modelSettings[division];
    nStates = m->numModelStates;

    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP;
    for (k=0; k<m->numOmegaCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numChars * m->numModelStates;
        }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* rescale */
    for (c=0; c<m->numChars; c++)
        {
        scaler = 0.0;
        for (k=0; k<m->numOmegaCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                if (clP[k][n] > scaler)
                    scaler = clP[k][n];
                }
            }

        for (k=0; k<m->numOmegaCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                clP[k][n] /= scaler;
                }
            clP[k] += n;
            }

        scP[c]       = (CLFlt) log (scaler);    /* store node scaler */
        lnScaler[c] += scP[c];                  /* add into tree scaler  */
        }

    m->unscaledNodes[chain][p->index] = 0;

    return (NO_ERROR);
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   CondLikeScaler_NY98_SSE: codon model with omega variation
|
-----------------------------------------------------------------*/
int CondLikeScaler_NY98_SSE (TreeNode *p, int division, int chain)
{
    int             c, k, n, nStates;
    CLFlt           *scP, *lnScaler;
    __m128          *clPtr, **clP, m1;
    ModelInfo       *m;

    m = &modelSettings[division];
    nStates = m->numModelStates;

    /* find conditional likelihood pointers */
    clPtr = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP_SSE;
    for (k=0; k<m->numOmegaCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
        }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];
    //scP_SSE = (__m128 *) scP;

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* rescale */
    for (c=0; c<m->numVecChars; c++)
        {
        //scaler = 0.0;
        m1 = _mm_setzero_ps ();
        for (k=0; k<m->numOmegaCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                m1 = _mm_max_ps (m1, clP[k][n]);
                }
            }
        _mm_store_ps (scP,  m1);
        scP += m->numFloatsPerVec;

        for (k=0; k<m->numOmegaCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                *clP[k] = _mm_div_ps (*clP[k], m1);
                clP[k]++;
                }
            }
        }
    
    /* Reset scP to original position*/
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];
    for (c=0; c<m->numChars; c++)
        {
        scP[c]       = (CLFlt) log (scP[c]);    /* store node scaler */
        lnScaler[c] += scP[c];                  /* add into tree scaler  */
        }

    m->unscaledNodes[chain][p->index] = 0;

    return (NO_ERROR);
}
#endif


/*----------------------------------------------------------------
|
|   CondLikeScaler_Std: variable states model with or without
|       rate variation
|
-----------------------------------------------------------------*/
int CondLikeScaler_Std (TreeNode *p, int division, int chain)
{
    int             c, n, k, nStates, numReps;
    CLFlt           scaler, *clPtr, **clP, *scP, *lnScaler;
    ModelInfo       *m;

    m = &modelSettings[division];

    numReps=0;
    for (c=0; c<m->numChars; c++)
        {
        if (m->nStates[c] == 2)
            numReps += m->numBetaCats * 2;
        else
            numReps += m->nStates[c];
        }

    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += numReps;
        }
    
    /* find node scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* rescale */
    for (c=0; c<m->numChars; c++)
        {
        scaler = 0.0;
        nStates = m->nStates[c];
        if (nStates == 2)
            nStates = m->numBetaCats * 2;

        for (k=0; k<m->numRateCats; k++)
            {
            for (n=0; n<nStates; n++)
                {
                if (clP[k][n] > scaler)
                    scaler = clP[k][n];
                }
            }

        for (k=0; k<m->numRateCats; k++)
            {
            for (n=0; n<nStates; n++)
                clP[k][n] /= scaler;
            clP[k] += nStates;
            }

        scP[c]       = (CLFlt) log (scaler);    /* store node scaler */
        lnScaler[c] += scP[c];                  /* add into tree scaler  */
        }

    m->unscaledNodes[chain][p->index] = 0;
        
    return NO_ERROR;
}


/*------------------------------------------------------------------
|
|   Likelihood_Adgamma: all n-state models with autocorrelated
|        discrete gamma rate variation, NOT morph, restriction,
|        codon or doublet models; just fill in rateProbs
|
-------------------------------------------------------------------*/
int Likelihood_Adgamma (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, j, k, i, nStates, nStatesDiv2;
    MrBFlt          *bs, *swr, s01, s10, probOn, probOff, covBF[40];
    MrBFlt          like, *rP;
    CLFlt           *clP;
    ModelInfo       *m;
    
    /* NOTE: whichSitePats offsets numSitesOfPat by whichSitePats X numCompressedChars.
       This is done so we can use the character reweighting scheme for "heating" chains. This was easy to
       accomplish for all of the models except this one, which doesn't use numSitesOfPat when calculating
       likelihoods. Either we disallow autocorrelated rates when using MCMC with character reweighting, or
       we properly calculate likelihoods when some site patterns have increased or decreased weight. For
       now, we do not allow MCMCMC with character reweighting with this HMM; we bail out in the function
       FillNumSitesOfPat if we have Adgamma rate variation and reweighting. */
    k = whichSitePats;  /* FIXME: Not used (from clang static analyzer) */
    
    /* find model settings */
    m = &modelSettings[division];
    
    /* get the number of states */
    nStates = m->numModelStates;
    nStatesDiv2 = nStates / 2;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);

    /* find conditional likelihood pointer */
    clP = m->condLikes[m->condLikeIndex[chain][p->index]];

    /* find pointer to rate probabilities */
    rP = rateProbs[chain] + state[chain] * rateProbRowSize + m->rateProbStart;

    /* loop over characters and calculate rate probs */
    if (m->switchRates != NULL)
        {
        swr = GetParamVals (m->switchRates, chain, state[chain]);
        s01 = swr[0];
        s10 = swr[1];
        probOn = s01 / (s01 + s10);
        probOff =  1.0 - probOn;
        for (j=0; j<nStatesDiv2; j++)
            {
            covBF[j] = bs[j] * probOn;
            covBF[j+nStatesDiv2] = bs[j] * probOff;
            }
        bs = covBF;
        }

    for (c=i=0; c<m->numChars; c++)
        {
        for (k=0; k<m->numRateCats; k++)
            {
            like =  0.0;
            for (j=0; j<nStates; j++)
                like += (*(clP++)) *  bs[j];
            rP[i++] = like;
            }
        }

    /* reset lnL, likelihood calculated later for this model */
    *lnL =  0.0;

    return (NO_ERROR);
}


/*------------------------------------------------------------------
|
|   Likelihood_Gen: general n-state models with or without rate
|       variation
|
-------------------------------------------------------------------*/
int Likelihood_Gen (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, j, k, nStates, hasPInvar;
    MrBFlt          s01, s10, probOn, probOff, *swr;
    MrBFlt          covBF[40], freq, *bs, like, likeI, pInvar=0.0, lnLike;
    CLFlt           *clPtr, **clP, *lnScaler, *nSitesOfPat, *clInvar=NULL;
    ModelInfo       *m;
    
    /* find model settings and nStates, pInvar, invar cond likes */
    m = &modelSettings[division];
    nStates = m->numModelStates;
    if (m->pInvar == NULL)
        {
        hasPInvar = NO;
        }
    else
        {
        hasPInvar = YES;
        pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
        clInvar = m->invCondLikes;
        }

    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numChars * m->numModelStates;
        }

    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);

    /* if covarion model, adjust base frequencies */
    if (m->switchRates != NULL)
        {
        /* find the stationary frequencies */
        swr = GetParamVals(m->switchRates, chain, state[chain]);
        s01 = swr[0];
        s10 = swr[1];
        probOn = s01 / (s01 + s10);
        probOff =  1.0 - probOn;

        /* now adjust the base frequencies; on-state stored first in cond likes */
        for (j=0; j<nStates/2; j++)
            {
            covBF[j] = bs[j] * probOn;
            covBF[j+nStates/2] = bs[j] * probOff;
            }

        /* finally set bs pointer to adjusted values */
        bs = covBF;
        }

    /* find category frequencies */
    if (hasPInvar == NO)
        freq =  1.0 /  m->numRateCats;
    else
        freq = (1.0 - pInvar) /  m->numRateCats;

    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* reset lnL */
    *lnL = 0.0;

    /* loop over characters */
    if (hasPInvar == NO)
        {
        for (c=0; c<m->numChars; c++)
            {
            like = 0.0;
            for (k=0; k<m->numRateCats; k++)
                for (j=0; j<nStates; j++)
                    {
                    like += (*(clP[k]++)) * bs[j];
#   ifdef DEBUG_LIKELIHOOD
                    // printf ("char=%d cat=%d j=%d like %E\n",c, k,j,like);
#   endif
                    }
            like *= freq;

            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }
    else
        {
        /* has invariable category */
        for (c=0; c<m->numChars; c++)
            {
            likeI = like = 0.0;
            for (k=0; k<m->numRateCats; k++)
                for (j=0; j<nStates; j++)
                    {
                    like += (*(clP[k]++)) * bs[j];
                    }
            like *= freq;
            for (j=0; j<nStates; j++)
                likeI += (*(clInvar++)) * bs[j] * pInvar;
            if (lnScaler[c] < -200.0)
                {
                /* we are not going to be able to exponentiate the scaling factor */
                if (likeI > 1E-70)
                    {
                    /* forget about like; it is going to be insignificant compared to likeI */
                    lnLike = log(likeI);
                    }
                else
                    {
                    /* treat likeI as if 0.0, that is, ignore it completely */
                    lnLike = log(like) + lnScaler[c];
                    }
                }
            else
                lnLike = log (like + (likeI / exp (lnScaler[c]))) + lnScaler[c];

            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += lnLike * nSitesOfPat[c];
                }
            }       
        }
        
    return NO_ERROR;
}


#if defined (SSE_ENABLED)
//#   if 0
//CLFlt DeleteME[1000];
//int PrintOld_SSE (TreeNode *p, int division, int chain){
//
//    int             c, c1, j, k, nStates;
//    //MrBFlt            *swr, likeI, pInvar=0.0, lnLike;
//    CLFlt           *temp_vector;
//    __m128          *clPtr, **clP;
//    ModelInfo       *m;
//
//    m = &modelSettings[division];
//    nStates = m->numModelStates;
//    /* find conditional likelihood pointers */
//
//    temp_vector =  DeleteME;
//
//    clPtr = (__m128 *) (m->condLikes[m->condLikeIndex[chain][p->index]]);
//    clP = m->clP_SSE;
//    for (k=0; k<m->numRateCats; k++)
//        {
//        clP[k] = clPtr;
//        clPtr += m->numVecChars * m->numModelStates;
//        }
//
//    for (c=0; c<m->numChars; c++)
//        {
//        c1 = c / FLOATS_PER_VEC;
//        for (k=0; k<m->numRateCats; k++)
//            {
//            for (j=0; j<nStates; j++)
//                {
//                *temp_vector++ = *(((CLFlt*)&clP[k][c1*nStates+j])+c % FLOATS_PER_VEC);
//                }
//            }
//        }
//    temp_vector=DeleteME;
//
//    return 1;
//}
//#   endif


/*------------------------------------------------------------------
|
|   Likelihood_Gen_SSE: general n-state model with or without rate
|       variation
|
-------------------------------------------------------------------*/
int Likelihood_Gen_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{   
    int             c, j, k, nStates, hasPInvar;
    MrBFlt          like, *bs;
    MrBFlt          s01, s10, probOn, probOff, *swr, covBF[40], freq, likeI, pInvar=0.0, lnLike;
    CLFlt           *lnScaler, *nSitesOfPat, *lnL_SSE, *lnLI_SSE;
    __m128          *clPtr, **clP, *clInvar=NULL;
    __m128          m1, mCatLike, mLike, mFreq;
    ModelInfo       *m;

    /* find model settings and nStates, pInvar, invar cond likes */
    m = &modelSettings[division];
    nStates = m->numModelStates;
    if (m->pInvar == NULL)
        {
        hasPInvar = NO;
        }
    else
        {
        hasPInvar = YES;
        pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
        clInvar = (__m128 *) (m->invCondLikes);
        }

    /* find conditional likelihood pointers */
    clPtr = (__m128 *) (m->condLikes[m->condLikeIndex[chain][p->index]]);
    clP = m->clP_SSE;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
        }
    lnL_SSE  = m->lnL_Vec;
    lnLI_SSE = m->lnLI_Vec;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);

    /* if covarion model, adjust base frequencies */
    if (m->switchRates != NULL)
        {
        /* find the stationary frequencies */
        swr = GetParamVals(m->switchRates, chain, state[chain]);
        s01 = swr[0];
        s10 = swr[1];
        probOn = s01 / (s01 + s10);
        probOff =  1.0 - probOn;

        /* now adjust the base frequencies; on-state stored first in cond likes */
        for (j=0; j<nStates/2; j++)
            {
            covBF[j] = bs[j] * probOn;
            covBF[j+nStates/2] = bs[j] * probOff;
            }

        /* finally set bs pointer to adjusted values */
        bs = covBF;
        }

    /* find category frequencies */
    if (hasPInvar == NO)
        freq =  1.0 /  m->numRateCats;
    else
        freq = (1.0 - pInvar) /  m->numRateCats;

    mFreq = _mm_set1_ps ((CLFlt)(freq));

    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* reset lnL */
    *lnL = 0.0;

    for (c=0; c<m->numVecChars; c++)
        {
        mLike = _mm_setzero_ps ();
        for (k=0; k<m->numRateCats; k++)
            {
            mCatLike = _mm_setzero_ps ();
            for (j=0; j<nStates; j++)
                {
                m1 = _mm_mul_ps (clP[k][j], _mm_set1_ps ((CLFlt)bs[j]));
                mCatLike = _mm_add_ps (mCatLike, m1);
                }
            m1 = _mm_mul_ps (mCatLike, mFreq);
            mLike = _mm_add_ps (mLike, m1);
            clP[k] += nStates;
            }
        _mm_store_ps (lnL_SSE, mLike);
        lnL_SSE += m->numFloatsPerVec;
        }

    /* loop over characters */
    if (hasPInvar == NO)
        {
        for (c=0; c<m->numChars; c++)
            {
            like = m->lnL_Vec[c];
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }
    else
        {
        /* has invariable category */
        for (c=0; c<m->numVecChars; c++)
            {
            mCatLike = _mm_setzero_ps ();
            for (j=0; j<nStates; j++)
                {
                m1 = _mm_mul_ps (clInvar[j], _mm_set1_ps ((CLFlt)bs[j]));
                mCatLike = _mm_add_ps (mCatLike, m1);
                }
            clInvar += nStates;
            _mm_store_ps (lnL_SSE, mCatLike);
            lnLI_SSE += m->numFloatsPerVec;
            }

        for (c=0; c<m->numChars; c++)
            {
            like  = m->lnL_Vec[c];
            likeI = m->lnLI_Vec[c];
            if (lnScaler[c] < -200.0)
                {
                /* we are not going to be able to exponentiate the scaling factor */
                if (likeI > 1E-70)
                    {
                    /* forget about like; it is going to be insignificant compared to likeI */
                    lnLike = log(likeI);
                    }
                else
                    {
                    /* treat likeI as if 0.0, that is, ignore it completely */
                    lnLike = log(like) + lnScaler[c];
                    }
                }
            else
                lnLike = log (like + (likeI / exp (lnScaler[c]))) + lnScaler[c];

            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += lnLike * nSitesOfPat[c];
                }
            }       
        }
        
    return NO_ERROR;
    
}
#endif


/*------------------------------------------------------------------
|
|   Likelihood_Gen_GibbsGamma: general n-state models using
|       Gibbs resampling of discrete gamma rate categories
|
-------------------------------------------------------------------*/
int Likelihood_Gen_GibbsGamma (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, j, nStates, nGammaCats, *rateCat;
    MrBFlt          s01, s10, probOn, probOff, *swr;
    MrBFlt          covBF[40], *bs, like;
    CLFlt           *clP, *lnScaler, *nSitesOfPat, *clInvar=NULL;
    ModelInfo       *m;
    
    /* find model settings, nStates and invar cond likes */
    m = &modelSettings[division];
    nStates = m->numModelStates;
    clInvar = m->invCondLikes;

    /* find conditional likelihood pointer */
    clP = m->condLikes[m->condLikeIndex[chain][p->index]];
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);

    /* if covarion model, adjust base frequencies */
    if (m->switchRates != NULL)
        {
        /* find the stationary frequencies */
        swr = GetParamVals(m->switchRates, chain, state[chain]);
        s01 = swr[0];
        s10 = swr[1];
        probOn = s01 / (s01 + s10);
        probOff =  1.0 - probOn;

        /* now adjust the base frequencies; on-state stored first in cond likes */
        for (j=0; j<nStates/2; j++)
            {
            covBF[j] = bs[j] * probOn;
            covBF[j+nStates/2] = bs[j] * probOff;
            }

        /* finally set bs pointer to adjusted values */
        bs = covBF;
        }

    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* find rate category index and number of gamma categories */
    rateCat = m->tiIndex + chain * m->numChars;
    nGammaCats = m->numRateCats;

    /* reset lnL */
    *lnL = 0.0;

    /* loop over characters */
    if (m->pInvar == NULL)
        {
        for (c=0; c<m->numChars; c++)
            {
            like = 0.0;
            for (j=0; j<nStates; j++)
                {
                like += (*(clP++)) * bs[j];
#   ifdef DEBUG_LIKELIHOOD
                // printf ("char=%d cat=%d j=%d like %E\n",c, k,j,like);
#   endif
                }

            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }
    else
        {
        /* has invariable category */
        for (c=0; c<m->numChars; c++)
            {
            like = 0.0;
            if (rateCat[c] < nGammaCats)
                {
                for (j=0; j<nStates; j++)
                    like += (*(clP++)) * bs[j];
                clInvar += nStates;
                }
            else
                {
                for (j=0; j<nStates; j++)
                    like += (*(clInvar++)) * bs[j];
                clP += nStates;
                }

            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (log(like) + lnScaler[c]) * nSitesOfPat[c];
                }
            }       
        }
        
    return NO_ERROR;
}


/*------------------------------------------------------------------
|
|   Likelihood_NUC4: 4by4 nucleotide models with or without rate
|       variation
|
-------------------------------------------------------------------*/
int Likelihood_NUC4 (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, k, hasPInvar;
    MrBFlt          freq, likeI, *bs, like, pInvar=0.0;
    CLFlt           *clPtr, **clP, *lnScaler, *nSitesOfPat, *clInvar=NULL;
    ModelInfo       *m;

    /* find model settings and pInvar, invar cond likes */
    m = &modelSettings[division];
    if (m->pInvar == NULL)
        {
        hasPInvar = NO;
        }
    else
        {
        hasPInvar = YES;
        pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
        clInvar = m->invCondLikes;
        }

    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numChars * m->numModelStates;
        }
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);

    /* find category frequencies */
    if (hasPInvar == NO)
        freq =  1.0 /  m->numRateCats;
    else
        freq =  (1.0 - pInvar) /  m->numRateCats;

    /* find tree scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* reset lnL */
    *lnL = 0.0;

    /* loop over characters */
    if (hasPInvar == NO)
        {
        for (c=0; c<m->numChars; c++)
            {
            like = 0.0;
            for (k=0; k<m->numRateCats; k++)
                {
                like += (clP[k][A] * bs[A] + clP[k][C] * bs[C] + clP[k][G] * bs[G] + clP[k][T] * bs[T]);
                clP[k] += 4;
                }
            like *= freq;
            
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }
    else
        {
        /* has invariable category */
        for (c=0; c<m->numChars; c++)
            {
            like = 0.0;
            for (k=0; k<m->numRateCats; k++)
                {
                like += (clP[k][A] * bs[A] + clP[k][C] * bs[C] + clP[k][G] * bs[G] + clP[k][T] * bs[T]);
                clP[k] += 4;
                }
            like *= freq;
            likeI = (clInvar[A] * bs[A] + clInvar[C] * bs[C] + clInvar[G] * bs[G] + clInvar[T] * bs[T]) * pInvar;
            if (lnScaler[c] < -200)
                {
                /* we are not going to be able to exponentiate the scaling factor */
                if (likeI > 1E-70)
                    {
                    /* forget about like; it is going to be insignificant compared to likeI */
                    like = likeI;
                    }
                else
                    {
                    /* treat likeI as if 0.0, that is, ignore it completely */
                    }
                }
            else
                like = like + (likeI / exp (lnScaler[c]));

            clInvar += 4;

            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }       
        }
        

    return NO_ERROR;
}


/*------------------------------------------------------------------
|
|   Likelihood_NUC4_GibbsGamma: 4by4 nucleotide models with rate
|       variation using Gibbs sampling from gamma rate categories
|
-------------------------------------------------------------------*/
int Likelihood_NUC4_GibbsGamma (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, i, r, nGammaCats, *rateCat;
    MrBFlt          *bs, like;
    CLFlt           *clP, *lnScaler, *nSitesOfPat, *clInvar;
    ModelInfo       *m;

    /* find model settings and invar cond likes */
    m = &modelSettings[division];
    clInvar = m->invCondLikes;

    /* find conditional likelihood pointer */
    clP = m->condLikes[m->condLikeIndex[chain][p->index]];
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);

    /* find tree scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* find rate category index  and number of gamma categories */
    rateCat = m->tiIndex + chain * m->numChars;
    nGammaCats = m->numRateCats;

    /* reset lnL */
    *lnL = 0.0;

    /* loop over characters */
    if (m->pInvar == NULL)
        {
        for (c=i=0; c<m->numChars; c++)
            {
            like = (clP[A] * bs[A] + clP[C] * bs[C] + clP[G] * bs[G] + clP[T] * bs[T]);
            clP += 4;
            
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }
    else
        {
        /* has invariable category */
        for (c=i=0; c<m->numChars; c++)
            {
            r = rateCat[c];
            if (r < nGammaCats)
                like = (clP[A] * bs[A] + clP[C] * bs[C] + clP[G] * bs[G] + clP[T] * bs[T]);
            else
                like = (clInvar[A] * bs[A] + clInvar[C] * bs[C] + clInvar[G] * bs[G] + clInvar[T] * bs[T]);
            clInvar += 4;
            clP += 4;

            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (log (like) + lnScaler[c]) * nSitesOfPat[c];
                }
            }       
        }

    return NO_ERROR;
}


//#if defined (SSE_ENABLED)
///*------------------------------------------------------------------
// |
// | Likelihood_NUC4_GibbsGamma: 4by4 nucleotide models with rate
// |     variation using Gibbs sampling from gamma rate categories
// |
// -------------------------------------------------------------------*/
//int Likelihood_NUC4_GibbsGamma_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
//{
//    int             c, i, r, nRateCats, *rateCat;
//    MrBFlt          *bs, like;
//    CLFlt           *lnScaler, *nSitesOfPat, *lnL_SSE, *lnLI_SSE;
//    __m128          *clP, *clInvar=NULL;
//    __m128          m1, mA, mC, mG, mT, mFreq, mPInvar, mLike;
//    ModelInfo       *m;
//    
//#if defined (FAST_LOG)
//    int             k, index;
//    MrBFlt          likeAdjust = 1.0, f;
//#endif
//    
//    /* find model settings and invar cond likes */
//    m = &modelSettings[division];
//    clInvar = (__m128 *)m->invCondLikes;
//    /* find conditional likelihood pointer */
//    clP = (__m128 *)m->condLikes[m->condLikeIndex[chain][p->index]];
//    
//    lnL_SSE  = m->lnL_SSE;
//    lnLI_SSE = m->lnLI_SSE;
//    
//    /* find base frequencies */
//    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
//    
//    /* find tree scaler */
//    lnScaler = m->scalers[m->siteScalerIndex[chain]];
//    
//    /* find nSitesOfPat */
//    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
//    
//    /* find rate category index  and number of rate categories */
//    rateCat = m->tiIndex + chain * m->numChars;
//    nRateCats = m->numRateCats;
//    
//    /* reset lnL */
//    *lnL = 0.0;
//    
//    /* calculate variable likelihood */
//    for (c=0; c<m->numVecChars; c++)
//    {
//        mLike = _mm_mul_ps (clP[A], mA);
//        m1    = _mm_mul_ps (clP[C], mC);
//        mLike = _mm_add_ps (mLike, m1);
//        m1    = _mm_mul_ps (clP[G], mG);
//        mLike = _mm_add_ps (mLike, m1);
//        m1    = _mm_mul_ps (clP[T], mT);
//        mLike = _mm_add_ps (mLike, m1);
//        
//        clP += 4;
//        _mm_store_ps (lnL_SSE, mLike);
//        lnL_SSE += FLOATS_PER_VEC;
//    }
//    
//    /* calculate invariable likelihood */
//    if (hasPInvar == YES)
//    {
//        for (c=0; c<m->numVecChars; c++)
//        {
//            mLike = _mm_mul_ps (clInvar[A], mA);
//            m1    = _mm_mul_ps (clInvar[C], mC);
//            mLike = _mm_add_ps (mLike, m1);
//            m1    = _mm_mul_ps (clInvar[G], mG);
//            mLike = _mm_add_ps (mLike, m1);
//            m1    = _mm_mul_ps (clInvar[T], mT);
//            mLike = _mm_add_ps (mLike, m1);
//            mLike = _mm_mul_ps (mLike, mPInvar);
//            
//            _mm_store_ps (lnLI_SSE, mLike);
//            clInvar += 4;
//            lnLI_SSE += FLOATS_PER_VEC;
//        }
//    }
//    
//    
//    /* loop over characters */
//    if (m->pInvar == NULL)
//    {
//        for (c=i=0; c<m->numChars; c++)
//        {
//            like = m->lnL_SSE[c];
//            /* check against LIKE_EPSILON (values close to zero are problematic) */
//            if (like < LIKE_EPSILON)
//            {
//                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30lf\n", spacer, division, c, like);
//                (*lnL) = MRBFLT_NEG_MAX;
//                return ERROR;
//            }
//            else
//            {
//#if defined (FAST_LOG)
//                f = frexp (like, &index);
//                index = 1-index;
//                (*lnL) += (lnScaler[c] +  logValue[index]) * nSitesOfPat[c];
//                for (k=0; k<(int)nSitesOfPat[c]; k++)
//                    likeAdjust *= f;
//#else
//                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
//#endif
//            }
//        }
//    }
//    else
//    {
//        /* has invariable category */
//        for (c=i=0; c<m->numChars; c++)
//        {
//            r = rateCat[c];
//            if (r < nRateCats)
//                like = m->lnL_SSE[c];
//            else
//                like = m->lnLI_SSE[c];
//            
//            /* check against LIKE_EPSILON (values close to zero are problematic) */
//            if (like < LIKE_EPSILON)
//            {
//                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30lf\n", spacer, division, c, like);
//                (*lnL) = MRBFLT_NEG_MAX;
//                return ERROR;
//            }
//            else
//            {
//                (*lnL) += (log (like) + lnScaler[c]) * nSitesOfPat[c];
//            }
//        }       
//    }
//    
//#if defined (FAST_LOG)
//    (*lnL) += log (likeAdjust);
//#endif
//    
//    return NO_ERROR;
//}
//#endif


#if defined (FMA_ENABLED)
/*------------------------------------------------------------------
 |
 |   Likelihood_NUC4_FMA: 4by4 nucleotide models with or without rate
 |       variation using AVX + FMA code
 |
 -------------------------------------------------------------------*/
int Likelihood_NUC4_FMA (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, k, hasPInvar;
    MrBFlt          freq, *bs, pInvar=0.0, like, likeI;
    CLFlt           *lnScaler, *nSitesOfPat, *lnL_Vec, *lnLI_Vec;
    __m256          *clPtr, **clP, *clInvar=NULL;
    __m256          mA, mC, mG, mT, mFreq, mPInvar=_mm256_set1_ps(0.0f), mLike;
    ModelInfo       *m;

    /* find model settings and pInvar, invar cond likes */
    m = &modelSettings[division];
    if (m->pInvar == NULL)
    {
        hasPInvar = NO;
    }
    else
    {
        hasPInvar = YES;
        pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
        mPInvar = _mm256_set1_ps ((CLFlt)(pInvar));
        clInvar = (__m256 *) (m->invCondLikes);
    }
    
    /* find conditional likelihood pointers */
    clPtr = (__m256 *) (m->condLikes[m->condLikeIndex[chain][p->index]]);
    clP = m->clP_AVX;
    for (k=0; k<m->numRateCats; k++)
    {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
    }
    lnL_Vec  = m->lnL_Vec;
    lnLI_Vec = m->lnLI_Vec;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    mA = _mm256_set1_ps ((CLFlt)(bs[A]));
    mC = _mm256_set1_ps ((CLFlt)(bs[C]));
    mG = _mm256_set1_ps ((CLFlt)(bs[G]));
    mT = _mm256_set1_ps ((CLFlt)(bs[T]));
    
    /* find category frequencies */
    if (hasPInvar == NO)
        freq =  1.0 / m->numRateCats;
    else
        freq =  (1.0 - pInvar) / m->numRateCats;
    mFreq = _mm256_set1_ps ((CLFlt)(freq));
    
    /* find tree scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* reset lnL */
    *lnL = 0.0;
    
    /* calculate variable likelihood */
    for (c=0; c<m->numVecChars; c++)
    {
        mLike = _mm256_setzero_ps ();
        for (k=0; k<m->numRateCats; k++)
        {
            mLike = _mm256_fmadd_ps (clP[k][A], mA, mLike);
            mLike = _mm256_fmadd_ps (clP[k][C], mC, mLike);
            mLike = _mm256_fmadd_ps (clP[k][G], mG, mLike);
            mLike = _mm256_fmadd_ps (clP[k][T], mT, mLike);
            clP[k] += 4;
        }
        mLike = _mm256_mul_ps (mLike, mFreq);
        _mm256_store_ps (lnL_Vec, mLike);
        lnL_Vec += m->numFloatsPerVec;
    }
    
    /* calculate invariable likelihood */
    if (hasPInvar == YES)
    {
        for (c=0; c<m->numVecChars; c++)
        {
            mLike = _mm256_mul_ps (clInvar[A], mA);
            mLike = _mm256_fmadd_ps (clInvar[C], mC, mLike);
            mLike = _mm256_fmadd_ps (clInvar[G], mG, mLike);
            mLike = _mm256_fmadd_ps (clInvar[T], mT, mLike);
            mLike = _mm256_mul_ps (mLike, mPInvar);
            _mm256_store_ps (lnLI_Vec, mLike);
            clInvar += 4;
            lnLI_Vec += m->numFloatsPerVec;
        }
    }
    
    /* accumulate results */
    if (hasPInvar == NO)
    {
        for (c=0; c<m->numChars; c++)
        {
            like = m->lnL_Vec[c];
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
            }
            else
            {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
            }
        }
    }
    else
    {
        /* has invariable category */
        for (c=0; c<m->numChars; c++)
        {
            like  = m->lnL_Vec[c];
            likeI = m->lnLI_Vec[c];
            if (lnScaler[c] < -200)
            {
                /* we are not going to be able to exponentiate the scaling factor */
                if (likeI > 1E-70)
                {
                    /* forget about like; it is going to be insignificant compared to likeI */
                    like = likeI;
                }
                else
                {
                    /* treat likeI as if 0.0, that is, ignore it completely */
                }
            }
            else
                like = like + (likeI / exp (lnScaler[c]));
            
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
            }
            else
            {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
            }
        }
    }
    
    return NO_ERROR;
}
#endif


#if defined (AVX_ENABLED)
/*------------------------------------------------------------------
 |
 |   Likelihood_NUC4_AVX: 4by4 nucleotide models with or without rate
 |       variation using AVX code
 |
 -------------------------------------------------------------------*/
int Likelihood_NUC4_AVX (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, k, hasPInvar;
    MrBFlt          freq, *bs, pInvar=0.0, like, likeI;
    CLFlt           *lnScaler, *nSitesOfPat, *lnL_Vec, *lnLI_Vec;
    __m256          *clPtr, **clP, *clInvar=NULL;
    __m256          m1, mA, mC, mG, mT, mFreq, mPInvar=_mm256_set1_ps(0.0f), mLike;
    ModelInfo       *m;
    
    /* find model settings and pInvar, invar cond likes */
    m = &modelSettings[division];
    if (m->pInvar == NULL)
    {
        hasPInvar = NO;
    }
    else
    {
        hasPInvar = YES;
        pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
        mPInvar = _mm256_set1_ps ((CLFlt)(pInvar));
        clInvar = (__m256 *) (m->invCondLikes);
    }
    
    /* find conditional likelihood pointers */
    clPtr = (__m256 *) (m->condLikes[m->condLikeIndex[chain][p->index]]);
    clP = m->clP_AVX;
    for (k=0; k<m->numRateCats; k++)
    {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
    }
    lnL_Vec  = m->lnL_Vec;
    lnLI_Vec = m->lnLI_Vec;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    mA = _mm256_set1_ps ((CLFlt)(bs[A]));
    mC = _mm256_set1_ps ((CLFlt)(bs[C]));
    mG = _mm256_set1_ps ((CLFlt)(bs[G]));
    mT = _mm256_set1_ps ((CLFlt)(bs[T]));
    
    /* find category frequencies */
    if (hasPInvar == NO)
        freq =  1.0 / m->numRateCats;
    else
        freq =  (1.0 - pInvar) / m->numRateCats;
    mFreq = _mm256_set1_ps ((CLFlt)(freq));
    
    /* find tree scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* reset lnL */
    *lnL = 0.0;
    
    /* calculate variable likelihood */
    for (c=0; c<m->numVecChars; c++)
    {
        mLike = _mm256_setzero_ps ();
        for (k=0; k<m->numRateCats; k++)
        {
            m1    = _mm256_mul_ps (clP[k][A], mA);
            mLike = _mm256_add_ps (mLike, m1);
            m1    = _mm256_mul_ps (clP[k][C], mC);
            mLike = _mm256_add_ps (mLike, m1);
            m1    = _mm256_mul_ps (clP[k][G], mG);
            mLike = _mm256_add_ps (mLike, m1);
            m1    = _mm256_mul_ps (clP[k][T], mT);
            mLike = _mm256_add_ps (mLike, m1);
            clP[k] += 4;
        }
        mLike = _mm256_mul_ps (mLike, mFreq);
        _mm256_store_ps (lnL_Vec, mLike);
        lnL_Vec += m->numFloatsPerVec;
    }
    
    /* calculate invariable likelihood */
    if (hasPInvar == YES)
    {
        for (c=0; c<m->numVecChars; c++)
        {
            mLike = _mm256_mul_ps (clInvar[A], mA);
            m1    = _mm256_mul_ps (clInvar[C], mC);
            mLike = _mm256_add_ps (mLike, m1);
            m1    = _mm256_mul_ps (clInvar[G], mG);
            mLike = _mm256_add_ps (mLike, m1);
            m1    = _mm256_mul_ps (clInvar[T], mT);
            mLike = _mm256_add_ps (mLike, m1);
            mLike = _mm256_mul_ps (mLike, mPInvar);
            
            _mm256_store_ps (lnLI_Vec, mLike);
            clInvar += 4;
            lnLI_Vec += m->numFloatsPerVec;
        }
    }
    
    /* accumulate results */
    if (hasPInvar == NO)
    {
        for (c=0; c<m->numChars; c++)
        {
            like = m->lnL_Vec[c];
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
            }
            else
            {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
            }
        }
    }
    else
    {
        /* has invariable category */
        for (c=0; c<m->numChars; c++)
        {
            like  = m->lnL_Vec[c];
            likeI = m->lnLI_Vec[c];
            if (lnScaler[c] < -200)
            {
                /* we are not going to be able to exponentiate the scaling factor */
                if (likeI > 1E-70)
                {
                    /* forget about like; it is going to be insignificant compared to likeI */
                    like = likeI;
                }
                else
                {
                    /* treat likeI as if 0.0, that is, ignore it completely */
                }
            }
            else
                like = like + (likeI / exp (lnScaler[c]));
            
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
            }
            else
            {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
            }
        }
    }
    
    return NO_ERROR;
}
#endif


#if defined (SSE_ENABLED)
/*------------------------------------------------------------------
|
|   Likelihood_NUC4_SSE: 4by4 nucleotide models with or without rate
|       variation
|
-------------------------------------------------------------------*/
int Likelihood_NUC4_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, k, hasPInvar;
    MrBFlt          freq, *bs, pInvar=0.0, like, likeI;
    CLFlt           *lnScaler, *nSitesOfPat, *lnL_Vec, *lnLI_Vec;
    __m128          *clPtr, **clP, *clInvar=NULL;
    __m128          m1, mA, mC, mG, mT, mFreq, mPInvar=_mm_set1_ps(0.0f), mLike;
    ModelInfo       *m;

    /* find model settings and pInvar, invar cond likes */
    m = &modelSettings[division];
    if (m->pInvar == NULL)
        {
        hasPInvar = NO;
        }
    else
        {
        hasPInvar = YES;
        pInvar =  *(GetParamVals (m->pInvar, chain, state[chain]));
        mPInvar = _mm_set1_ps ((CLFlt)(pInvar));
        clInvar = (__m128 *) (m->invCondLikes);
        }

    /* find conditional likelihood pointers */
    clPtr = (__m128 *) (m->condLikes[m->condLikeIndex[chain][p->index]]);
    clP = m->clP_SSE;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
        }
    lnL_Vec  = m->lnL_Vec;
    lnLI_Vec = m->lnLI_Vec;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    mA = _mm_set1_ps ((CLFlt)(bs[A]));
    mC = _mm_set1_ps ((CLFlt)(bs[C]));
    mG = _mm_set1_ps ((CLFlt)(bs[G]));
    mT = _mm_set1_ps ((CLFlt)(bs[T]));

    /* find category frequencies */
    if (hasPInvar == NO)
        freq =  1.0 / m->numRateCats;
    else
        freq =  (1.0 - pInvar) / m->numRateCats;
    mFreq = _mm_set1_ps ((CLFlt)(freq));

    /* find tree scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* reset lnL */
    *lnL = 0.0;

    /* calculate variable likelihood */
    for (c=0; c<m->numVecChars; c++)
        {
        mLike = _mm_setzero_ps ();
        for (k=0; k<m->numRateCats; k++)
            {
            m1    = _mm_mul_ps (clP[k][A], mA);
            mLike = _mm_add_ps (mLike, m1);
            m1    = _mm_mul_ps (clP[k][C], mC);
            mLike = _mm_add_ps (mLike, m1);
            m1    = _mm_mul_ps (clP[k][G], mG);
            mLike = _mm_add_ps (mLike, m1);
            m1    = _mm_mul_ps (clP[k][T], mT);
            mLike = _mm_add_ps (mLike, m1);
            clP[k] += 4;
            }
        mLike = _mm_mul_ps (mLike, mFreq);
        _mm_store_ps (lnL_Vec, mLike);
        lnL_Vec += m->numFloatsPerVec;
        }
    
    /* calculate invariable likelihood */
    if (hasPInvar == YES)
        {
        for (c=0; c<m->numVecChars; c++)
            {
            mLike = _mm_mul_ps (clInvar[A], mA);
            m1    = _mm_mul_ps (clInvar[C], mC);
            mLike = _mm_add_ps (mLike, m1);
            m1    = _mm_mul_ps (clInvar[G], mG);
            mLike = _mm_add_ps (mLike, m1);
            m1    = _mm_mul_ps (clInvar[T], mT);
            mLike = _mm_add_ps (mLike, m1);
            mLike = _mm_mul_ps (mLike, mPInvar);

            _mm_store_ps (lnLI_Vec, mLike);
            clInvar += 4;
            lnLI_Vec += m->numFloatsPerVec;
            }
        }

    /* accumulate results */
    if (hasPInvar == NO)
        {
        for (c=0; c<m->numChars; c++)
            {
            like = m->lnL_Vec[c];
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }
    else
        {
        /* has invariable category */
        for (c=0; c<m->numChars; c++)
            {
            like  = m->lnL_Vec[c];
            likeI = m->lnLI_Vec[c];
            if (lnScaler[c] < -200)
                {
                /* we are not going to be able to exponentiate the scaling factor */
                if (likeI > 1E-70)
                    {
                    /* forget about like; it is going to be insignificant compared to likeI */
                    like = likeI;
                    }
                else
                    {
                    /* treat likeI as if 0.0, that is, ignore it completely */
                    }
                }
            else
                like = like + (likeI / exp (lnScaler[c]));

            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }

    return NO_ERROR;
}
#endif


/*------------------------------------------------------------------
|
|   Likelihood_NY98: Codon model with three selection categories,
|       after Nielsen and Yang (1998).
|
-------------------------------------------------------------------*/
int Likelihood_NY98 (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, j, k, nStates;
    MrBFlt          catLike, like, *bs, *omegaCatFreq;
    CLFlt           **clP,*clPtr, *lnScaler, *nSitesOfPat;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* number of states */
    nStates = m->numModelStates;

    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP;
    for (k=0; k<m->numOmegaCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numChars * m->numModelStates;
        }
    
    /* find codon frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* find category frequencies */
    omegaCatFreq = GetParamSubVals (m->omega, chain, state[chain]);

    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    *lnL = 0.0; /* reset lnL */

    for (c=m->numDummyChars; c<m->numChars; c++)
        {
        like = 0.0;
        for (k=0; k<m->numOmegaCats; k++)
            {
            catLike = 0.0;
            for (j=0; j<nStates; j++)
                catLike += clP[k][j] * bs[j];
            like += catLike * omegaCatFreq[k];
            clP[k] += nStates;
            }
        /* check against LIKE_EPSILON (values close to zero are problematic) */
        if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
            MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
            (*lnL) = MRBFLT_NEG_MAX;
            abortMove = YES;
            return ERROR;
            }
        else
            {
            (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
            }
        }

    return NO_ERROR;
}


#if defined (SSE_ENABLED)
/*------------------------------------------------------------------
|
|   Likelihood_NY98_SSE: Codon model with three selection categories,
|       after Nielsen and Yang (1998).
|
-------------------------------------------------------------------*/
int Likelihood_NY98_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, j, k, nStates;
    MrBFlt          like, *bs, *omegaCatFreq;
    CLFlt           *lnScaler, *nSitesOfPat, *lnL_Vec;
    __m128          *clPtr, **clP;
    __m128          m1, mCatLike, mLike;
    ModelInfo       *m;
    
    m = &modelSettings[division];

    /* number of states */
    nStates = m->numModelStates;

    /* find conditional likelihood pointers */
    clPtr = (__m128 *) m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP_SSE;
    for (k=0; k<m->numOmegaCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * nStates;
        }
    
    /* find codon frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* find category frequencies */
    omegaCatFreq = GetParamSubVals (m->omega, chain, state[chain]);

    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    *lnL = 0.0; /* reset lnL */

    lnL_Vec  = m->lnL_Vec;
    for (c=0; c<m->numVecChars; c++)
        {
        mLike = _mm_setzero_ps ();
        for (k=0; k<m->numOmegaCats; k++)
            {
            mCatLike = _mm_setzero_ps ();
            for (j=0; j<nStates; j++)
                {
                m1 = _mm_mul_ps (clP[k][j], _mm_set1_ps ((CLFlt)bs[j]));
                mCatLike = _mm_add_ps (mCatLike, m1);
                }
            m1 = _mm_mul_ps (mCatLike, _mm_set1_ps ((CLFlt)omegaCatFreq[k]));
            mLike = _mm_add_ps (mLike, m1);
            clP[k] += nStates;
            }
        _mm_store_ps (lnL_Vec, mLike);
        lnL_Vec += m->numFloatsPerVec;
        }
    for (c=m->numDummyChars; c<m->numChars; c++)
        {
        like = m->lnL_Vec[c];
        /* check against LIKE_EPSILON (values close to zero are problematic) */
        if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
            MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
            (*lnL) = MRBFLT_NEG_MAX;
            abortMove = YES;
            return ERROR;
            }
        else    
            {
            (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
            }
        }

    return NO_ERROR;
}
#endif


/*------------------------------------------------------------------
|
|   Likelihood_Res: restriction site model with or without rate
|       variation
|
-------------------------------------------------------------------*/
int Likelihood_Res (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, k;
    MrBFlt          *bs, freq, like, pUnobserved, pObserved;
    CLFlt           *clPtr, **clP, *lnScaler, *nSitesOfPat;
    ModelInfo       *m;

    
    m = &modelSettings[division];

    /* find conditional likelihood pointer */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numChars * m->numModelStates;
        }

    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);

    /* find category frequencies */
    freq =  1.0 /  m->numRateCats;

    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    *lnL = 0.0; /* reset lnL */

    pUnobserved = 0.0;
    for (c=0; c<m->numDummyChars; c++)
        {
        like = 0.0;
        for (k=0; k<m->numRateCats; k++)
            {
            like += (clP[k][0]*bs[0] + clP[k][1]*bs[1]) * freq;
            clP[k] += 2;
            }
        pUnobserved += like *  exp(lnScaler[c]);
        }

    pObserved =  1.0 - pUnobserved;
    if (pObserved < LIKE_EPSILON)
        {
#   ifdef DEBUG_LIKELIHOOD
        MrBayesPrint ("%s   WARNING: p(Observed) < LIKE_EPSILON - for division %d p(Observed) = %1.30le\n", spacer, division+1, pObserved);
#   endif
        (*lnL) = MRBFLT_NEG_MAX;
        abortMove = YES;
        return ERROR;
        }

    for (c=m->numDummyChars; c<m->numChars; c++)
        {
        like = 0.0;
        for (k=0; k<m->numRateCats; k++)
            {
            like += (clP[k][0]*bs[0] + clP[k][1]*bs[1]) * freq;
            clP[k] += 2;
            }
        /* check against LIKE_EPSILON (values close to zero are problematic) */
        if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
            MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
            (*lnL) = MRBFLT_NEG_MAX;
            abortMove = YES;
            return ERROR;
            }
        else    
            {
            (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
            }
        }

    /* correct for absent characters */
    (*lnL) -=  log(pObserved) * (m->numUncompressedChars);

    return NO_ERROR;
}


#if defined (SSE_ENABLED)
/*------------------------------------------------------------------
|
|   Likelihood_Res_SSE: 4by4 nucleotide models with or without rate
|       variation
|
-------------------------------------------------------------------*/
int Likelihood_Res_SSE (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, k;
    MrBFlt          freq, *bs, like, pUnobserved, pObserved;
    CLFlt           *lnScaler, *nSitesOfPat, *lnL_Vec;
    __m128          *clPtr, **clP;
    __m128          m1, mA, mB, mFreq, mLike;
    ModelInfo       *m;

    /* find model settings and pInvar, invar cond likes */
    m = &modelSettings[division];

    /* find conditional likelihood pointers */
    clPtr = (__m128 *) (m->condLikes[m->condLikeIndex[chain][p->index]]);
    clP = m->clP_SSE;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += m->numVecChars * m->numModelStates;
        }
    lnL_Vec  = m->lnL_Vec;
    
    /* find base frequencies */
    bs = GetParamSubVals (m->stateFreq, chain, state[chain]);
    mA = _mm_set1_ps ((CLFlt)(bs[0]));
    mB = _mm_set1_ps ((CLFlt)(bs[1]));

    freq =  1.0 / m->numRateCats;
    mFreq = _mm_set1_ps ((CLFlt)(freq));

    /* find tree scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    /* reset lnL */
    *lnL = 0.0;

    /* calculate variable likelihood */
    for (c=0; c<m->numVecChars; c++)
        {
        mLike = _mm_setzero_ps ();
        for (k=0; k<m->numRateCats; k++)
            {
            m1    = _mm_mul_ps (clP[k][0], mA);
            mLike = _mm_add_ps (mLike, m1);
            m1    = _mm_mul_ps (clP[k][1], mB);
            mLike = _mm_add_ps (mLike, m1);
            clP[k] += 2;
            }
        mLike = _mm_mul_ps (mLike, mFreq);
        _mm_store_ps (lnL_Vec, mLike);
        lnL_Vec += m->numFloatsPerVec;
        }

    pUnobserved = 0.0;
    for (c=0; c<m->numDummyChars; c++)
        {
        like  = m->lnL_Vec[c];
        pUnobserved += like *  exp(lnScaler[c]);
        }

    pObserved =  1.0 - pUnobserved;
    if (pObserved < LIKE_EPSILON)
        {
#   ifdef DEBUG_LIKELIHOOD
        MrBayesPrint ("%s   WARNING: p(Observed) < LIKE_EPSILON - for division %d p(Observed) = %1.30le\n", spacer, division+1, pObserved);
#   endif
        (*lnL) = MRBFLT_NEG_MAX;
        abortMove = YES;
        return ERROR;
        }

    for (c=m->numDummyChars; c<m->numChars; c++)
        {
        like  = m->lnL_Vec[c];
        /* check against LIKE_EPSILON (values close to zero are problematic) */
        if (like < LIKE_EPSILON)
            {
#   ifdef DEBUG_LIKELIHOOD
            MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
            (*lnL) = MRBFLT_NEG_MAX;
            abortMove = YES;
            return ERROR;
            }
        else    
            {
            (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
            }
        }

    /* correct for absent characters */
    (*lnL) -=  log(pObserved) * (m->numUncompressedChars);

    return NO_ERROR;
}
#endif


/*------------------------------------------------------------------
|
|   Likelihood_Std: variable states model with or without rate
|       variation
|
-------------------------------------------------------------------*/
int Likelihood_Std (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             b, c, j, k, nBetaCats, nRateCats, nStates, numReps;
    MrBFlt          catLike, catFreq, rateFreq, like, *bs, *bsBase,
                    pUnobserved, pObserved;
    CLFlt           *clPtr, **clP, *lnScaler, *nSitesOfPat;
    ModelInfo       *m;

    m = &modelSettings[division];

    numReps=0;
    for (c=0; c<m->numChars; c++)
        {
        if (m->nStates[c] == 2)
            numReps += m->numBetaCats * 2;
        else
            numReps += m->nStates[c];
        }
    /* find conditional likelihood pointers */
    clPtr = m->condLikes[m->condLikeIndex[chain][p->index]];
    clP   = m->clP;
    for (k=0; k<m->numRateCats; k++)
        {
        clP[k] = clPtr;
        clPtr += numReps;
        }
    
    /* find base frequencies */
    bsBase = GetParamStdStateFreqs (m->stateFreq, chain, state[chain]);

    /* find rate category number and frequencies */
    nRateCats = m->numRateCats;
    rateFreq = 1.0 / nRateCats;

    /* find site scaler */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];
    
    /* find nSitesOfPat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    
    *lnL = 0.0; /* reset lnL */

    if (m->numBetaCats == 1)
        {
        pUnobserved = 0.0;
        catFreq = rateFreq;
        for (c=j=0; c<m->numDummyChars; c++)
            {
            like = 0.0;
            nStates = m->nStates[c];
            bs = bsBase + m->bsIndex[c];
            for (k=0; k<nRateCats; k++)
                {
                catLike = 0.0;
                for (j=0; j<nStates; j++)
                    catLike += clP[k][j] * bs[j];
                like += catLike * catFreq;
                clP[k] += nStates;
                }
            pUnobserved += like *  exp(lnScaler[c]);
            }

        pObserved =  1.0 - pUnobserved;
        if (pObserved < LIKE_EPSILON)
            pObserved = LIKE_EPSILON;

        for (c=m->numDummyChars; c<m->numChars; c++)
            {
            like = 0.0;
            nStates = m->nStates[c];
            bs = bsBase + m->bsIndex[c];

            for (k=0; k<nRateCats; k++)
                {
                catLike = 0.0;
                for (j=0; j<nStates; j++)
                    catLike += clP[k][j] * bs[j];
                like += catLike * catFreq;
                clP[k] += nStates;
                }
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }
    else
        {
        pUnobserved = 0.0;
        for (c=j=0; c<m->numDummyChars; c++)
            {
            like = 0.0;
            nStates = m->nStates[c];
            bs = bsBase + m->bsIndex[c];
            if (nStates == 2)
                {
                nBetaCats = m->numBetaCats;
                catFreq = rateFreq / nBetaCats;
                }
            else
                {
                nBetaCats = 1;
                catFreq = rateFreq;
                }
            for (b=0; b<nBetaCats; b++)
                {
                for (k=0; k<nRateCats; k++)
                    {
                    catLike = 0.0;
                    for (j=0; j<nStates; j++)
                        catLike += clP[k][j] * bs[j];
                    like += catLike * catFreq;
                    clP[k] += nStates;
                    }
                bs += nStates;
                }
            pUnobserved += like *  exp(lnScaler[c]);
            }

        pObserved =  1.0 - pUnobserved;
        if (pObserved < LIKE_EPSILON)
            pObserved = LIKE_EPSILON;

        for (c=m->numDummyChars; c<m->numChars; c++)
            {
            like = 0.0;
            nStates = m->nStates[c];
            bs = bsBase + m->bsIndex[c];
            if (nStates == 2)
                {
                nBetaCats = m->numBetaCats;
                catFreq = rateFreq / nBetaCats;
                }
            else
                {
                nBetaCats = 1;
                catFreq = rateFreq;
                }
            for (b=0; b<nBetaCats; b++)
                {
                for (k=0; k<nRateCats; k++)
                    {
                    catLike = 0.0;
                    for (j=0; j<nStates; j++)
                        catLike += clP[k][j] * bs[j];
                    like += catLike * catFreq;
                    clP[k] += nStates;
                    }
                bs += nStates;
                }
            /* check against LIKE_EPSILON (values close to zero are problematic) */
            if (like < LIKE_EPSILON)
                {
#   ifdef DEBUG_LIKELIHOOD
                MrBayesPrint ("%s   WARNING: In LIKE_EPSILON - for division %d char %d has like = %1.30le\n", spacer, division+1, c+1, like);
#   endif
                (*lnL) = MRBFLT_NEG_MAX;
                abortMove = YES;
                return ERROR;
                }
            else    
                {
                (*lnL) += (lnScaler[c] +  log(like)) * nSitesOfPat[c];
                }
            }
        }

    /* correct for absent characters */
    (*lnL) -=  log(pObserved) * (m->numUncompressedChars);

    return NO_ERROR;
}


/*------------------------------------------------------------------
|
|   Likelihood_Pars: likelihood under the Tuffley and Steel (1997)
|       model for characters with constant number of states. The idea
|       is described in:
|
|       Tuffley, C., and M. Steel. 1997. Links between maximum likelihood
|          and maximum parsimony under a simple model of site substitution.
|          Bull. Math. Bio. 59:581-607.
|
|       The likelihood under the Tuffley and Steel (1997) model is:
|       
|       L = k^[-(T + n)]
|      
|       where L is the likelihood
|             k is the number of character states
|             T is the parsimony tree length
|             n is the number of characters 
|
|   The parsimony calculator does not use character packing; this is
|       to enable reweighting of characters 
|
|   Note that this is an empirical Bayes approach in that it uses the
|       maximum likelihood branch length.
|
-------------------------------------------------------------------*/
int Likelihood_Pars (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, i, nStates;
    BitsLong        done, *pL, *pR, *pP, *pA, *oldpP, x;
    CLFlt           nParsChars, treeLength;
    CLFlt           length, *nSitesOfPat, *newNodeLength, oldNodeLength;
    Tree            *t;
    ModelInfo       *m;

    /* Find model settings */
    m = &modelSettings[division];

    /* Get tree */
    t = GetTree(m->brlens,chain,state[chain]);
    
    /* Get parsimony tree length */
    treeLength = (CLFlt) m->parsTreeLength[2 * chain + state[chain]];
    
    /* Get number of states */
    nStates = m->numStates;

    /* Get number of sites of pat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;

    /* Mark the nodes that can be stop nodes                 */
    /* (there must not be any touched side nodes below them) */
    p = t->root;
    p->marked = YES;
    for (i=t->nIntNodes-1; i>=0; i--)
        {
        p = t->intDownPass[i];
        p->marked = NO;
        if (p->upDateCl == YES && p->anc->marked == YES)
            {
            if (p->left->upDateCl == NO || p->right->upDateCl == NO)
                p->marked = YES;
            }
        }

    /* Now make downpass node by node */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];

        /* continue if no work needs to be done */
        if (p->upDateCl == NO)
            continue;

        /* flip space */
        FlipCondLikeSpace(m, chain, p->index);
        
        /* find parsimony sets for the node and its environment */
        pL    = m->parsSets[m->condLikeIndex[chain][p->left->index ]];
        pR    = m->parsSets[m->condLikeIndex[chain][p->right->index]];
        oldpP = m->parsSets[m->condLikeScratchIndex[p->index       ]];
        pP    = m->parsSets[m->condLikeIndex[chain][p->index       ]];

        /* find old and new node lengths */
        oldNodeLength =  m->parsNodeLens[m->condLikeScratchIndex[p->index]];
        newNodeLength = &m->parsNodeLens[m->condLikeIndex[chain][p->index]];
        
        if (t->isRooted == NO && p->anc->anc == NULL)
            {
            pA = m->parsSets[m->condLikeIndex[chain][p->anc->index]];
            length = 0.0;
            for (c=0; c<m->numChars; c++)
                {
                x = pL[c] & pR[c];
                if (x == 0)
                    {
                    x = pL[c] | pR[c];
                    length += nSitesOfPat[c];
                    }
                if ((x & pA[c]) == 0)
                    length += nSitesOfPat[c];
                pP[c] = x;
                }
            treeLength += (length - oldNodeLength);
            newNodeLength[0] = length;
            }
        else
            {
            length = 0.0;
            done = 0;
            for (c=0; c<m->numChars; c++)
                {
                x = pL[c] & pR[c];
                if (x == 0)
                    {
                    x = pL[c] | pR[c];
                    length += nSitesOfPat[c];
                    }
                pP[c] = x;
                done |= (x^oldpP[c]);
                }
            treeLength += (length - oldNodeLength);
            newNodeLength[0] = length;
            if (p->marked == YES && done == 0)
                break;
            }
        }

    /* Count number of characters in the partition. It is calculated
       on the fly because this number is going to differ for
       different chains if character reweighting is used. */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;
    nParsChars = 0.0;
    for (c=0; c<m->numChars; c++)
        nParsChars += nSitesOfPat[c];

    /* Calculate likelihood from parsimony tree length */
    *lnL = - ((treeLength + nParsChars) *  log (nStates));

    /* Store current parsimony tree length */
    m->parsTreeLength[2 * chain + state[chain]] = treeLength;

    return (NO_ERROR);
}


#if 0
int Likelihood_ParsCodon (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             x, y;
    TreeNode        *q;
    
    /* no warnings */
    q = p;
    x = division;
    y = chain;
    *lnL = 0.0;
    x = whichSitePats;

    MrBayesPrint ("%s   Parsimony calculator for codons not yet implemented\n", spacer);
    
    return ERROR;
}
#   endif


/*------------------------------------------------------------------
|
|   Likelihood_Pars: likelihood under the Tuffley and Steel (1997)
|       model for characters with constant number of states. The idea
|       is described in:
|
|       Tuffley, C., and M. Steel. 1997. Links between maximum likelihood
|          and maximum parsimony under a simple model of site substitution.
|          Bull. Math. Bio. 59:581-607.
|
|       The likelihood under the Tuffley and Steel (1997) model is:
|       
|       L = k^[-(T + n)]
|      
|       where L is the likelihood
|             k is the number of character states
|             T is the parsimony tree length
|             n is the number of characters 
|
|   The parsimony calculator does not use character packing; this is
|       to enable reweighting of characters 
|
|   Note that this is an empirical Bayes approach in that it uses the
|       maximum likelihood branch length.
|
|   This variant of the calculator assumes that the number of states
|       is variable. It does not take state order into account.
|
-------------------------------------------------------------------*/
int Likelihood_ParsStd (TreeNode *p, int division, int chain, MrBFlt *lnL, int whichSitePats)
{
    int             c, i, *nStates;
    BitsLong        *pL, *pR, *pP, *pA, x;
    CLFlt           *treeLength;
    CLFlt           *nSitesOfPat;
    Tree            *t;
    ModelInfo       *m;

    /* Find model settings */
    m = &modelSettings[division];

    /* Get tree */
    t = GetTree(m->brlens,chain,state[chain]);
    
    /* Allocate space for parsimony tree length */
    treeLength = (CLFlt *) SafeCalloc (m->numChars, sizeof (CLFlt));
    
    /* Get number of states */
    nStates = m->nStates;

    /* Get number of sites of pat */
    nSitesOfPat = numSitesOfPat + (whichSitePats*numCompressedChars) + m->compCharStart;

    /* Make downpass node by node; do not skip any nodes */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];

        /* flip space */
        FlipCondLikeSpace(m, chain, p->index);
        
        /* find parsimony sets for the node and its environment */
        pL    = m->parsSets[m->condLikeIndex[chain][p->left->index ]];
        pR    = m->parsSets[m->condLikeIndex[chain][p->right->index]];
        pP    = m->parsSets[m->condLikeIndex[chain][p->index       ]];

        if (t->isRooted == NO && p->anc->anc == NULL)
            {
            pA = m->parsSets[m->condLikeIndex[chain][p->anc->index]];
            for (c=0; c<m->numChars; c++)
                {
                x = pL[c] & pR[c];
                if (x == 0)
                    {
                    x = pL[c] | pR[c];
                    treeLength[c] += nSitesOfPat[c];
                    }
                if ((x & pA[c]) == 0)
                    treeLength[c] += nSitesOfPat[c];
                pP[c] = x;
                }
            }
        else
            {
            for (c=0; c<m->numChars; c++)
                {
                x = pL[c] & pR[c];
                if (x == 0)
                    {
                    x = pL[c] | pR[c];
                    treeLength[c] += nSitesOfPat[c];
                    }
                pP[c] = x;
                }
            }
        }

    /* Calculate the likelihood one character at a time */
    *lnL = 0.0;
    for (c=0; c<m->numChars; c++)
        {
        *lnL -= ((treeLength[c] + nSitesOfPat[c]) * log (nStates[c]));
        }

    /* Free space for parsimony character states */
    free (treeLength);

    return (NO_ERROR);
}


/*-----------------------------------------------------------------
|
|   LaunchLogLikeForDivision: calculate the log likelihood of the 
|       new state of the chain for a single division
|
-----------------------------------------------------------------*/
void LaunchLogLikeForDivision(int chain, int d, MrBFlt* lnL)
{
    int i;
    TreeNode        *p;
    ModelInfo       *m;
    Tree            *tree;
#   if defined (TIMING_ANALIZ)
    clock_t         CPUTimeStart;
#   endif
    
    m = &modelSettings[d];
    tree = GetTree(m->brlens, chain, state[chain]);
    
    if (m->upDateCijk == YES)
        {
        if (UpDateCijk(d, chain)== ERROR)
            {
            (*lnL) = MRBFLT_NEG_MAX; /* effectively abort the move */
            return;
            }
        m->upDateAll = YES;
        }
    
#   if defined (BEAGLE_ENABLED)
    if (m->useBeagle == YES)
        {
        LaunchBEAGLELogLikeForDivision(chain, d, m, tree, lnL);
        return;
        }
#   endif
        
    /* Flip and copy or reset site scalers */
    FlipSiteScalerSpace(m, chain);
    if (m->upDateAll == YES)
        ResetSiteScalers(m, chain);
    else
        CopySiteScalers(m, chain);
    
    if (m->parsModelId == NO)
        {
        for (i=0; i<tree->nIntNodes; i++)
            {
            p = tree->intDownPass[i];
            
            if (p->left->upDateTi == YES)
                {
                /* shift state of ti probs for node */
                FlipTiProbsSpace (m, chain, p->left->index);
                m->TiProbs (p->left, d, chain);
                }
            
            if (p->right->upDateTi == YES)
                {
                /* shift state of ti probs for node */
                FlipTiProbsSpace (m, chain, p->right->index);
                m->TiProbs (p->right, d, chain);
                }
            
            if (tree->isRooted == NO)
                {
                if (p->anc->anc == NULL /* && p->upDateTi == YES */)
                    {
                    /* shift state of ti probs for node */
                    FlipTiProbsSpace (m, chain, p->index);
                    m->TiProbs (p, d, chain);
                    }
                }
            
            if (p->upDateCl == YES)
                {
                if (tree->isRooted == NO)
                    {
                    if (p->anc->anc == NULL)
                        {
                        TIME(m->CondLikeRoot (p, d, chain),CPUCondLikeRoot);
                        }
                    else
                        {
                        TIME(m->CondLikeDown (p, d, chain),CPUCondLikeDown);                        
                        }
                    }
                else
                    {
                    TIME(m->CondLikeDown (p, d, chain),CPUCondLikeDown);
                    }

                if (m->unscaledNodes[chain][p->index] == 0 && m->upDateAll == NO)
                    {
#if defined (SSE_ENABLED)
                    if (m->useVec == VEC_SSE)
                        {
                        TIME(RemoveNodeScalers_SSE (p, d, chain),CPUScalersRemove);
                        }
#if defined (AVX_ENABLED)
                    else if (m->useVec == VEC_AVX)
                        {
                        TIME(RemoveNodeScalers_AVX (p, d, chain),CPUScalersRemove);
                        }
#endif
                    else
                        {
                        TIME(RemoveNodeScalers (p, d, chain),CPUScalersRemove);
                        }
#   else
                    TIME(RemoveNodeScalers (p, d, chain),CPUScalersRemove);
#   endif
                    }
                FlipNodeScalerSpace (m, chain, p->index);
                m->unscaledNodes[chain][p->index] = 1 + m->unscaledNodes[chain][p->left->index] + m->unscaledNodes[chain][p->right->index];
                
                if (m->unscaledNodes[chain][p->index] >= m->rescaleFreq[chain] && p->anc->anc != NULL)
                    {
                    TIME(m->CondLikeScaler (p, d, chain),CPUScalers);
                    }
                }
            }
        }
    TIME(m->Likelihood (tree->root->left, d, chain, lnL, (chainId[chain] % chainParams.numChains)),CPULilklihood);
    return;
}


/*----------------------------------------------------------------
|
|   RemoveNodeScalers: Remove node scalers
|
-----------------------------------------------------------------*/
int RemoveNodeScalers (TreeNode *p, int division, int chain)
{
    int             c;
    CLFlt           *scP, *lnScaler;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    assert (m->unscaledNodes[chain][p->index] == 0);

    /* find scalers */
    scP = m->scalers[m->nodeScalerIndex[chain][p->index]];

    /* find site scalers */
    lnScaler = m->scalers[m->siteScalerIndex[chain]];

    /* remove scalers */
    for (c=0; c<m->numChars; c++)
        lnScaler[c] -= scP[c];

    return NO_ERROR;
}


#if defined (AVX_ENABLED)
/*----------------------------------------------------------------
 |
 |   RemoveNodeScalers_AVX: Remove node scalers, AVX code
 |
 -----------------------------------------------------------------*/
int RemoveNodeScalers_AVX (TreeNode *p, int division, int chain)
{
    int             c;
    __m256          *scP_AVX, *lnScaler_AVX;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    assert (m->unscaledNodes[chain][p->index] == 0);
    
    /* find scalers */
    scP_AVX = (__m256*)(m->scalers[m->nodeScalerIndex[chain][p->index]]);
    
    /* find site scalers */
    lnScaler_AVX = (__m256*)(m->scalers[m->siteScalerIndex[chain]]);
    
    /* remove scalers */
    for (c=0; c<m->numVecChars; c++)
    {
        lnScaler_AVX[c] = _mm256_sub_ps(lnScaler_AVX[c], scP_AVX[c]);
    }
    
    return NO_ERROR;
    
}
#endif


#if defined (SSE_ENABLED)
/*----------------------------------------------------------------
|
|   RemoveNodeScalers_SSE: Remove node scalers, SSE code
|
-----------------------------------------------------------------*/
int RemoveNodeScalers_SSE (TreeNode *p, int division, int chain)
{
    int             c;
    __m128          *scP_SSE, *lnScaler_SSE;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    assert (m->unscaledNodes[chain][p->index] == 0);

    /* find scalers */
    scP_SSE = (__m128*)(m->scalers[m->nodeScalerIndex[chain][p->index]]);

    /* find site scalers */
    lnScaler_SSE = (__m128*)(m->scalers[m->siteScalerIndex[chain]]);

    /* remove scalers */
    for (c=0; c<m->numVecChars; c++)
        {
        lnScaler_SSE[c] = _mm_sub_ps(lnScaler_SSE[c], scP_SSE[c]);
        }

    return NO_ERROR;
    
}
#endif


int SetBinaryQMatrix (MrBFlt **a, int whichChain, int division)
{
    MrBFlt          scaler, *bs;
    ModelInfo       *m;
        
    /* set up pointers to the appropriate model information */
    m = &modelSettings[division];
    assert (m->numModelStates == 2);

    bs = GetParamSubVals (m->stateFreq, whichChain, state[whichChain]);
    scaler = 1.0 / (2*bs[0]*bs[1]);
    a[0][0]= -bs[1]*scaler;
    a[0][1]=  bs[1]*scaler;
    a[1][0]=  bs[0]*scaler;
    a[1][1]= -bs[0]*scaler;

    return (NO_ERROR);
}


int SetNucQMatrix (MrBFlt **a, int n, int whichChain, int division, MrBFlt rateMult, MrBFlt *rA, MrBFlt *rS)
{
    register int    i, j, k;
    int             isTransition=0, nDiff, rtNum=0;
    MrBFlt          scaler, mult=0.0, probOn, sum, *swr, s01, s10, s[4][4], nonsyn, *rateValues=NULL, *bs, dN, dS;
    ModelInfo       *m;
    ModelParams     *mp;
#   if defined BEAGLE_ENABLED
    MrBFlt          trans;
#   endif

    /* set up pointers to the appropriate model information */
    mp = &modelParams[division];
    m = &modelSettings[division];
    assert (m->numModelStates == n);

    /* All of the models that are set up in this function require the frequencies
       of the nucleotides (or doublets or codons). They will also require either
       a transition/transversion rate ratio or the GTR rate parameters. The 
       "rateValues" will either be
       
          rateValues[0] = transtion/transversion rate (kappa)
       
       for nst=2 models or
       
          rateValues[0] = A <-> C rate
          rateValues[1] = A <-> G rate
          rateValues[2] = A <-> T rate
          rateValues[3] = C <-> G rate
          rateValues[4] = C <-> T rate
          rateValues[5] = G <-> T rate
          
       for nst=6 models. */
    bs = GetParamSubVals (m->stateFreq, whichChain, state[whichChain]);
    if (m->nst == 2)
        {
        rateValues = GetParamVals(m->tRatio, whichChain, state[whichChain]);
#   if defined (BEAGLE_ENABLED)
        /* transversions assumed to have rate 1.0; */
        trans = rateValues[0];
        if (m->numModelStates == 4)   /* code to satisfy Beagle */
            {
            rateValues = (MrBFlt *) SafeCalloc (6, sizeof(MrBFlt));
            rateValues[0] = rateValues[2] = rateValues[3] = rateValues[5] =1.0; /* Setting transversions */
            rateValues[1] = rateValues[4] = trans; /* Setting transitions */
            }
#   endif
        }

    else if (m->nst == 6 || m->nst == NST_MIXED)
        rateValues = GetParamVals(m->revMat, whichChain, state[whichChain]);
#   if defined (BEAGLE_ENABLED)
    else if (m->nst == 1 && m->numModelStates == 4)   /* code to satisfy Beagle */
        {
        rateValues = (MrBFlt *) SafeCalloc (6, sizeof(MrBFlt));
        for (i=0; i<6; i++)
            rateValues[i] = 1.0;
        }
#   endif

    if (n == 4) 
        {
        /* 4 X 4 model:
        
           Here, we set the rate matrix for the GTR model (Tavare, 1986). We
           need not only the 6 rates for this model (rateValues), but also the 
           base frequencies (bs). */
            
        /* set diagonal of Q matrix to 0 */
        for (i=0; i<4; i++)
            a[i][i] = 0.0;

        /* initialize Q matrix */
        scaler = 0.0;
        for (i=0; i<4; i++)
            {
            for (j=i+1; j<4; j++)
                {
                if (i == 0 && j == 1)
                    mult = rateValues[0];
                else if (i == 0 && j == 2)
                    mult = rateValues[1];
                else if (i == 0 && j == 3)
                    mult = rateValues[2];
                else if (i == 1 && j == 2)
                    mult = rateValues[3];
                else if (i == 1 && j == 3)
                    mult = rateValues[4];
                else if (i == 2 && j == 3)
                    mult = rateValues[5];
                a[i][i] -= (a[i][j] = bs[j] * mult);
                a[j][j] -= (a[j][i] = bs[i] * mult);
                scaler += bs[i] * a[i][j];
                scaler += bs[j] * a[j][i];
                }
            }
            
        /* rescale Q matrix */
        scaler = 1.0 / scaler;
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
                a[i][j] *= scaler;
        }
    else if (n == 8) /* we have a 4 X 4 covarion model */
        {
        /* 8 X 8 covarion model:
        
           Here, we set the rate matrix for the covarion model (Tuffley and
           Steel, 1997). We need the rate parameters of the model 
           (contained in rateValues), the frequencies of the four nucleotides,
           and the switching rates to completely specify the rate matrix. We
           first set up the 4 X 4 submatrix that represents changes (the upper
           left portion of the 8 X 8 matrix). Note that if we have rate
           variation across sites, that we need to deal with the multiplication
           in the rate matrix (i.e., we cannot simply deal with rate variation
           by multiplying the branch length by a rate multiplier as we can
           with other models). Instead, we multiply the scaled rate matrix
           by the rate multiplier. */

        /* Get the switching rates. The rate of off->on is s01 and the rate
           of on->off is s10. The stationary probability of the switch process
           is prob1 = s01/(s01+s10) and prob0 = s10/(s01+s10). */
        swr = GetParamVals (m->switchRates, whichChain, state[whichChain]);
        s01 = swr[0];
        s10 = swr[1];
        probOn = s01 / (s01 + s10);
        
        /* set matrix a to 0 */
        for (i=0; i<8; i++)
            for (j=0; j<8; j++)
                a[i][j] = 0.0;

        /* set up the 4 X 4 matrix representing substitutions (s[][]; upper left) */
        if (m->nst == 1)
            {
            scaler = 0.0;
            for (i=0; i<4; i++)
                {
                for (j=i+1; j<4; j++)
                    {
                    s[i][j] = bs[j];
                    s[j][i] = bs[i];
                    scaler += bs[i] * s[i][j] * probOn;
                    scaler += bs[j] * s[j][i] * probOn;
                    }
                }
            }
        else if (m->nst == 2)
            {
            scaler = 0.0;
            for (i=0; i<4; i++)
                {
                for (j=i+1; j<4; j++)
                    {
                    if ((i == 0 && j == 2) || (i == 2 && j == 0) || (i == 1 && j == 3) || (i == 3 && j == 1))
                        mult = rateValues[0];
                    else
                        mult = 1.0;
                    s[i][j] = bs[j] * mult;
                    s[j][i] = bs[i] * mult;
                    scaler += bs[i] * s[i][j] * probOn;
                    scaler += bs[j] * s[j][i] * probOn;
                    }
                }
            }
        else
            {
            scaler = 0.0;
            for (i=0; i<4; i++)
                {
                for (j=i+1; j<4; j++)
                    {
                    if (i == 0 && j == 1)
                        mult = rateValues[0];
                    else if (i == 0 && j == 2)
                        mult = rateValues[1];
                    else if (i == 0 && j == 3)
                        mult = rateValues[2];
                    else if (i == 1 && j == 2)
                        mult = rateValues[3];
                    else if (i == 1 && j == 3)
                        mult = rateValues[4];
                    else if (i == 2 && j == 3)
                        mult = rateValues[5];

                    s[i][j] = bs[j] * mult;
                    s[j][i] = bs[i] * mult;
                    scaler += bs[i] * s[i][j] * probOn;
                    scaler += bs[j] * s[j][i] * probOn;
                    }
                }
            }

        /* rescale off diagonal elements of s[][] matrix */
        scaler = 1.0 / scaler;
        for (i=0; i<4; i++)
            {
            for (j=0; j<4; j++)
                {
                if (i != j)
                    s[i][j] *= scaler;
                }
            }
            
        /* now, scale s[][] by rate factor */
        for (i=0; i<4; i++)
            {
            for (j=0; j<4; j++)
                {
                if (i != j)
                    s[i][j] *= rateMult;
                }
            }
            
        /* put in diagonal elements of s[][] */
        for (i=0; i<4; i++)
            {
            sum = 0.0;
            for (j=0; j<4; j++)
                {
                if (i != j)
                    sum += s[i][j];
                }
            s[i][i] = -(sum + s10);
            }
                
        /* Now, put s[][] into top left portion of a matrix and fill in the
           other parts of the matrix with the appropriate switching rates. */
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
                a[i][j] = s[i][j];
        for (i=4; i<8; i++)
            a[i][i] = -s01;
        a[0][4] = s10;
        a[1][5] = s10;
        a[2][6] = s10;
        a[3][7] = s10;
        a[4][0] = s01;
        a[5][1] = s01;
        a[6][2] = s01;
        a[7][3] = s01;
        
#       if 0
        for (i=0; i<8; i++)
            {
            for (j=0; j<8; j++)
                printf ("%1.10lf ", a[i][j]);
            printf ("\n");
            }
        for (i=0; i<4; i++)
            printf ("%lf ", bs[i]);
        printf ("\n");
        printf ("s01 = %lf s10 = %lf pi1 = %lf pi0 = %lf\n", s01, s10, probOn, 1-probOn);
#       endif
        }
    else if (n == 16) 
        {
        /* 16 X 16 doublet model:
        
           We have a doublet model. The states are in the order AA, AC, AG, AT, CA, CC
           CG, CT, GA, GC, GG, GT, TA, TC, TG, TT. The rate matrix is straight-forward
           to set up. We simply multiply the rate parameter (e.g., the ti/tv rate
           ratio) by the doublet frequencies. */
           
        /* set diagonal of Q matrix to 0 */
        for (i=0; i<16; i++)
            a[i][i] = 0.0;

        if (m->nst == 1) /* F81-like doublet model */
            {
            scaler = 0.0;
            for (i=0; i<16; i++)
                {
                for (j=i+1; j<16; j++)
                    {
                    if (((doublet[i].first & doublet[j].first) == 0) && ((doublet[i].second & doublet[j].second) == 0))
                        mult = 0.0;
                    else
                        mult = 1.0;                 
                    a[i][i] -= (a[i][j] = bs[j] * mult);
                    a[j][j] -= (a[j][i] = bs[i] * mult);
                    scaler += bs[i] * a[i][j];
                    scaler += bs[j] * a[j][i];
                    }
                }
            }
        else if (m->nst == 2) /* HKY-like doublet model */
            {
            scaler = 0.0;
            for (i=0; i<16; i++)
                {
                for (j=i+1; j<16; j++)
                    {
                    if (((doublet[i].first & doublet[j].first) == 0) && ((doublet[i].second & doublet[j].second) == 0))
                        mult = 0.0;
                    else
                        {
                        if ((doublet[i].first & doublet[j].first) == 0)
                            {
                            if ((doublet[i].first + doublet[j].first) == 5 || (doublet[i].first + doublet[j].first) == 10)
                                mult = rateValues[0];
                            else
                                mult = 1.0;
                            }
                        else
                            {
                            if ((doublet[i].second + doublet[j].second) == 5 || (doublet[i].second + doublet[j].second) == 10)
                                mult = rateValues[0];
                            else
                                mult = 1.0;
                            }
                        }               
                    a[i][i] -= (a[i][j] = bs[j] * mult);
                    a[j][j] -= (a[j][i] = bs[i] * mult);
                    scaler += bs[i] * a[i][j];
                    scaler += bs[j] * a[j][i];
                    }
                }
            }
        else /* GTR-like doublet model */
            {
            scaler = 0.0;
            for (i=0; i<16; i++)
                {
                for (j=i+1; j<16; j++)
                    {
                    if (((doublet[i].first & doublet[j].first) == 0) && ((doublet[i].second & doublet[j].second) == 0))
                        mult = 0.0;
                    else
                        {
                        if ((doublet[i].first & doublet[j].first) == 0)
                            {
                            if ((doublet[i].first + doublet[j].first) == 3)
                                mult = rateValues[0];
                            else if ((doublet[i].first + doublet[j].first) == 5)
                                mult = rateValues[1];
                            else if ((doublet[i].first + doublet[j].first) == 9)
                                mult = rateValues[2];
                            else if ((doublet[i].first + doublet[j].first) == 6)
                                mult = rateValues[3];
                            else if ((doublet[i].first + doublet[j].first) == 10)
                                mult = rateValues[4];
                            else
                                mult = rateValues[5];
                            }
                        else
                            {
                            if ((doublet[i].second + doublet[j].second) == 3)
                                mult = rateValues[0];
                            else if ((doublet[i].second + doublet[j].second) == 5)
                                mult = rateValues[1];
                            else if ((doublet[i].second + doublet[j].second) == 9)
                                mult = rateValues[2];
                            else if ((doublet[i].second + doublet[j].second) == 6)
                                mult = rateValues[3];
                            else if ((doublet[i].second + doublet[j].second) == 10)
                                mult = rateValues[4];
                            else
                                mult = rateValues[5];
                            }
                        }               
                    a[i][i] -= (a[i][j] = bs[j] * mult);
                    a[j][j] -= (a[j][i] = bs[i] * mult);
                    scaler += bs[i] * a[i][j];
                    scaler += bs[j] * a[j][i];
                    }
                }
            }
                    
            
        /* rescale Q matrix */
        scaler = 1.0 / scaler;
        for (i=0; i<16; i++)
            for (j=0; j<16; j++)
                a[i][j] *= scaler;
        }
    else
        {
        /* 64(ish) X 64(ish) codon model:
        
           Here, we set the rate matrix for the codon model (see Goldman and
           Yang, 1994). Note that we can specifiy any general type of codon
           model, with these constraints:
           
            a[i][j] = 0                      -> if i and j differ at 2 or 3 nucleotides
            a[i][j] = rateValues[0] * bs[j]  -> if synonymous A <-> C change
            a[i][j] = rateValues[1] * bs[j]  -> if synonymous A <-> G change
            a[i][j] = rateValues[2] * bs[j]  -> if synonymous A <-> T change
            a[i][j] = rateValues[3] * bs[j]  -> if synonymous C <-> G change
            a[i][j] = rateValues[4] * bs[j]  -> if synonymous C <-> T change
            a[i][j] = rateValues[5] * bs[j]  -> if synonymous G <-> T change
            
            a[i][j] = rateValues[0] * nonsyn * bs[j]  -> if nonsynonymous A <-> C change
            a[i][j] = rateValues[1] * nonsyn * bs[j]  -> if nonsynonymous A <-> G change
            a[i][j] = rateValues[2] * nonsyn * bs[j]  -> if nonsynonymous A <-> T change
            a[i][j] = rateValues[3] * nonsyn * bs[j]  -> if nonsynonymous C <-> G change
            a[i][j] = rateValues[4] * nonsyn * bs[j]  -> if nonsynonymous C <-> T change
            a[i][j] = rateValues[5] * nonsyn * bs[j]  -> if nonsynonymous G <-> T change
            
          Other models, such as the one used by Nielsen & Yang (1998) can be obtained
          from this model by restricing transitions and transversions to have the same rate.
          nonsyn is the nonsynonymous/synonymous rate ratio (often called the
          dN/dS ratio). If we are in this part of the function, then we rely on it
          being called with the "rateMult" parameter specifying the dN/dS ratio. Note
          that the size of the matrix will never be 64 X 64 as we only consider changes
          among coding triplets (i.e., we exclude the stop codons). */
          
        /* get the nonsynonymous/synonymous rate ratio */
        nonsyn = rateMult; 
        
        /* set diagonal of Q matrix to 0 */
        for (i=0; i<n; i++)
            a[i][i] = 0.0;
            
        /* set dN and dS rates to zero */
        dN = dS = 0.0;

        if (m->nst == 1) /* F81-like codon model */
            {
            scaler = 0.0;
            for (i=0; i<n; i++)
                {
                for (j=i+1; j<n; j++)
                    {
                    nDiff = 0;
                    for (k=0; k<3; k++)
                        {
                        if (mp->codonNucs[i][k] != mp->codonNucs[j][k])
                            nDiff++;
                        }
                    if (nDiff > 1)
                        {
                        mult = 0.0;
                        }
                    else
                        {
                        if (mp->codonAAs[i] == mp->codonAAs[j])
                            mult = 1.0;
                        else
                            mult = nonsyn;
                        }
                    
                    a[i][i] -= (a[i][j] = bs[j] * mult);
                    a[j][j] -= (a[j][i] = bs[i] * mult);
                    if (mp->codonAAs[i] == mp->codonAAs[j])
                        dS += (bs[i] * a[i][j] + bs[j] * a[j][i]);
                    else
                        dN += (bs[i] * a[i][j] + bs[j] * a[j][i]);
                    scaler += bs[i] * a[i][j];
                    scaler += bs[j] * a[j][i];
                    }
                }
            }
        else if (m->nst == 2) /* HKY-like codon model */
            {
            scaler = 0.0;
            for (i=0; i<n; i++)
                {
                for (j=i+1; j<n; j++)
                    {
                    nDiff = 0;
                    for (k=0; k<3; k++)
                        {
                        if (mp->codonNucs[i][k] != mp->codonNucs[j][k])
                            {
                            nDiff++;
                            if ((mp->codonNucs[i][k] == 0 && mp->codonNucs[j][k] == 2) || (mp->codonNucs[i][k] == 2 && mp->codonNucs[j][k] == 0) ||
                                (mp->codonNucs[i][k] == 1 && mp->codonNucs[j][k] == 3) || (mp->codonNucs[i][k] == 3 && mp->codonNucs[j][k] == 1))
                                isTransition = YES;
                            else
                                isTransition = NO;
                            }
                        }
                    if (nDiff > 1)
                        {
                        mult = 0.0;
                        }
                    else
                        {
                        if (mp->codonAAs[i] == mp->codonAAs[j])
                            mult = 1.0;
                        else
                            mult = nonsyn;
                        if (isTransition == YES)
                            mult *= rateValues[0];
                        }
                    
                    a[i][i] -= (a[i][j] = bs[j] * mult);
                    a[j][j] -= (a[j][i] = bs[i] * mult);
                    if (mp->codonAAs[i] == mp->codonAAs[j])
                        dS += (bs[i] * a[i][j] + bs[j] * a[j][i]);
                    else
                        dN += (bs[i] * a[i][j] + bs[j] * a[j][i]);
                    scaler += bs[i] * a[i][j];
                    scaler += bs[j] * a[j][i];
                    }
                }
            }
        else /* GTR-like codon model */
            {
            scaler = 0.0;
            for (i=0; i<n; i++)
                {
                for (j=i+1; j<n; j++)
                    {
                    nDiff = 0;
                    for (k=0; k<3; k++)
                        {
                        if (mp->codonNucs[i][k] != mp->codonNucs[j][k])
                            {
                            nDiff++;
                            if ((mp->codonNucs[i][k] == 0 && mp->codonNucs[j][k] == 1) || (mp->codonNucs[i][k] == 1 && mp->codonNucs[j][k] == 0))
                                rtNum = 0;
                            else if ((mp->codonNucs[i][k] == 0 && mp->codonNucs[j][k] == 2) || (mp->codonNucs[i][k] == 2 && mp->codonNucs[j][k] == 0))
                                rtNum = 1;
                            else if ((mp->codonNucs[i][k] == 0 && mp->codonNucs[j][k] == 3) || (mp->codonNucs[i][k] == 3 && mp->codonNucs[j][k] == 0))
                                rtNum = 2;
                            else if ((mp->codonNucs[i][k] == 1 && mp->codonNucs[j][k] == 2) || (mp->codonNucs[i][k] == 2 && mp->codonNucs[j][k] == 1))
                                rtNum = 3;
                            else if ((mp->codonNucs[i][k] == 1 && mp->codonNucs[j][k] == 3) || (mp->codonNucs[i][k] == 3 && mp->codonNucs[j][k] == 1))
                                rtNum = 4;
                            else
                                rtNum = 5;
                            }
                        }
                    if (nDiff > 1)
                        {
                        mult = 0.0;
                        }
                    else
                        {
                        if (mp->codonAAs[i] == mp->codonAAs[j])
                            mult = 1.0;
                        else
                            mult = nonsyn;
                        if (rtNum == 0)
                            mult *= rateValues[0];
                        else if (rtNum == 1)
                            mult *= rateValues[1];
                        else if (rtNum == 2)
                            mult *= rateValues[2];
                        else if (rtNum == 3)
                            mult *= rateValues[3];
                        else if (rtNum == 4)
                            mult *= rateValues[4];
                        else
                            mult *= rateValues[5];
                        }
                    
                    a[i][i] -= (a[i][j] = bs[j] * mult);
                    a[j][j] -= (a[j][i] = bs[i] * mult);
                    if (mp->codonAAs[i] == mp->codonAAs[j])
                        dS += (bs[i] * a[i][j] + bs[j] * a[j][i]);
                    else
                        dN += (bs[i] * a[i][j] + bs[j] * a[j][i]);
                    scaler += bs[i] * a[i][j];
                    scaler += bs[j] * a[j][i];
                    }
                }
            }

        /* rescale Q matrix */
        if (m->nucModelId == NUCMODEL_CODON && m->numOmegaCats > 1)
            {
            /* If we have a positive selection model with multiple categories, then
               we do not rescale the rate matrix until we have finished generating
               all of the necessary rate matrices. The rescaling occurrs in 
               UpDateCijk. */
            (*rA) = dN;
            (*rS) = dS;
            }
        else
            {
            scaler = 1.0 / scaler;
            for (i=0; i<n; i++)
                for (j=0; j<n; j++)
                    a[i][j] *= scaler;
            (*rA) = (*rS) = 1.0;
            }           
        }

#   if 0
    for (i=0; i<n; i++)
        {
        for (j=0; j<n; j++)
            printf ("%0.5lf ", a[i][j]);
        printf ("\n");
        }
#   endif

#   if defined (BEAGLE_ENABLED)
    if ((m->nst == 1 || m->nst == 2) && m->numModelStates == 4)
        free (rateValues);
#   endif

    return (NO_ERROR);
}


int SetProteinQMatrix (MrBFlt **a, int n, int whichChain, int division, MrBFlt rateMult)
{
    register int    i, j, k;
    int             aaModelID;
    MrBFlt          scaler, probOn, sum, *swr, s01, s10, *bs, *rt;
    ModelInfo       *m;
        
    /* set up pointers to the appropriate model information */
    m = &modelSettings[division];

    /* get amino acid model ID 
        AAMODEL_POISSON         0
        AAMODEL_JONES           1
        AAMODEL_DAY             2
        AAMODEL_MTREV           3
        AAMODEL_MTMAM           4
        AAMODEL_WAG             5
        AAMODEL_RTREV           6
        AAMODEL_CPREV           7
        AAMODEL_VT              8
        AAMODEL_BLOSUM          9
        AAMODEL_LG             10
        AAMODEL_EQ             11
        AAMODEL_GTR            12 */
        
    if (m->aaModelId >= 0)
        aaModelID = m->aaModelId;
    else
        aaModelID = (int)*GetParamVals(m->aaModel, whichChain, state[whichChain]);
    
    /* Make certain that we have either 20 or 40 states. Anything
       else means we have a real problem. */
    if (n != 20 && n != 40)
        {
        MrBayesPrint ("%s   ERROR: There should be 20 or 40 states for the aa model\n");
        return (ERROR);
        }

    if (n == 20)
        {
        /* We have a run-of-the-mill amino acid model (i.e., 20 X 20). */
        if (aaModelID == AAMODEL_POISSON)
            {
            scaler = 1.0 / 19.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = scaler;
                    a[j][i] = scaler;
                    }
                }
            for (i=0; i<20; i++)
                a[i][i] = -1.0;
            }
        else if (aaModelID == AAMODEL_EQ)
            {
            bs = GetParamSubVals (m->stateFreq, whichChain, state[whichChain]);
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = 0.0;
            scaler = 0.0;   
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][i] -= (a[i][j] = bs[j]);
                    a[j][j] -= (a[j][i] = bs[i]);
                    scaler += bs[i] * a[i][j];
                    scaler += bs[j] * a[j][i];
                    }
                }
            scaler = 1.0 / scaler;
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] *= scaler;
            }
        else if (aaModelID == AAMODEL_GTR)
            {
            bs = GetParamSubVals (m->stateFreq, whichChain, state[whichChain]);
            rt = GetParamVals (m->revMat, whichChain, state[whichChain]);
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = 0.0;
            scaler = 0.0;
            for (i=k=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][i] -= (a[i][j] = bs[j] * rt[k]);
                    a[j][j] -= (a[j][i] = bs[i] * rt[k]);
                    k++;
                    }
                }
            for (i=0; i<20; i++)
                scaler += -(bs[i] * a[i][i]);
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] /= scaler;
            }
        else if (aaModelID == AAMODEL_JONES)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aaJones[i][j];
            }
        else if (aaModelID == AAMODEL_DAY)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aaDayhoff[i][j];
            }
        else if (aaModelID == AAMODEL_MTREV)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aaMtrev24[i][j];
            }
        else if (aaModelID == AAMODEL_MTMAM)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aaMtmam[i][j];
            }
        else if (aaModelID == AAMODEL_RTREV)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aartREV[i][j];
            }
        else if (aaModelID == AAMODEL_WAG)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aaWAG[i][j];
            }
        else if (aaModelID == AAMODEL_CPREV)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aacpREV[i][j];
            }
        else if (aaModelID == AAMODEL_VT)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aaVt[i][j];
            }
        else if (aaModelID == AAMODEL_BLOSUM)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aaBlosum[i][j];
            }
        else if (aaModelID == AAMODEL_LG)
            {
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = aaLG[i][j];
            }
        else
            {
            MrBayesPrint ("%s   ERROR: Don't understand which amino acid model is needed\n");
            return (ERROR);
            }
#       if 0
        for (i=0; i<20; i++)
            {
            for (j=0; j<20; j++)
                printf ("%1.3lf ", a[i][j]);
            printf ("\n");
            }
#       endif
        }
    else
        {
        /* 40 X 40 covarion model:
        
           We have a covarion model, and must set up the other quadrants. Note that if
           we are at this point in the code, that we have already set up the upper left
           portion of the 40 X 40 rate matrix. Note that if we have rate
           variation across sites, that we need to deal with the multiplication
           in the rate matrix (i.e., we cannot simply deal with rate variation
           by multiplying the branch length by a rate multiplier as we can
           with other models). Instead, we multiply the scaled rate matrix
           by the rate multiplier. */

        /* Get the switching rates. The rate of off->on is s01 and the rate
           of on->off is s10. The stationary probability of the switch process
           is prob1 = s01/(s01+s10) and prob0 = s10/(s01+s10). */
        swr = GetParamVals (m->switchRates, whichChain, state[whichChain]);
        s01 = swr[0];
        s10 = swr[1];
        probOn = s01 / (s01 + s10);
        
        /* set matrix a[][] to 0 */
        for (i=0; i<40; i++)
            for (j=0; j<40; j++)
                a[i][j] = 0.0;  
                
        /* fill in upper-left sub matrix (where substitutions occur */
        if (aaModelID == AAMODEL_POISSON)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = 0.05;
                    a[j][i] = 0.05;
                    scaler += 0.05 * a[i][j] * probOn;
                    scaler += 0.05 * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_EQ)
            {
            bs = GetParamSubVals (m->stateFreq, whichChain, state[whichChain]);
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = bs[j];
                    a[j][i] = bs[i];
                    scaler += bs[i] * a[i][j] * probOn;
                    scaler += bs[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_GTR)
            {
            bs = GetParamSubVals (m->stateFreq, whichChain, state[whichChain]);
            rt = GetParamVals (m->revMat, whichChain, state[whichChain]);
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] = 0.0;
            scaler = 0.0;
            for (i=k=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][i] -= (a[i][j] = bs[j] * rt[k]);
                    a[j][j] -= (a[j][i] = bs[i] * rt[k]);
                    k++;
                    }
                }
            for (i=0; i<20; i++)
                scaler += -(bs[i] * a[i][i]);
            for (i=0; i<20; i++)
                for (j=0; j<20; j++)
                    a[i][j] /= scaler;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = bs[j];
                    a[j][i] = bs[i];
                    scaler += bs[i] * a[i][j] * probOn;
                    scaler += bs[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_JONES)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aaJones[i][j];
                    a[j][i] = aaJones[j][i];
                    scaler += jonesPi[i] * a[i][j] * probOn;
                    scaler += jonesPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_DAY)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aaDayhoff[i][j];
                    a[j][i] = aaDayhoff[j][i];
                    scaler += dayhoffPi[i] * a[i][j] * probOn;
                    scaler += dayhoffPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_MTREV)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aaMtrev24[i][j];
                    a[j][i] = aaMtrev24[j][i];
                    scaler += mtrev24Pi[i] * a[i][j] * probOn;
                    scaler += mtrev24Pi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_MTMAM)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aaMtmam[i][j];
                    a[j][i] = aaMtmam[j][i];
                    scaler += mtmamPi[i] * a[i][j] * probOn;
                    scaler += mtmamPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_RTREV)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aartREV[i][j];
                    a[j][i] = aartREV[j][i];
                    scaler += rtrevPi[i] * a[i][j] * probOn;
                    scaler += rtrevPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_WAG)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aaWAG[i][j];
                    a[j][i] = aaWAG[j][i];
                    scaler += wagPi[i] * a[i][j] * probOn;
                    scaler += wagPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_CPREV)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aacpREV[i][j];
                    a[j][i] = aacpREV[j][i];
                    scaler += cprevPi[i] * a[i][j] * probOn;
                    scaler += cprevPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_VT)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aaVt[i][j];
                    a[j][i] = aaVt[j][i];
                    scaler += vtPi[i] * a[i][j] * probOn;
                    scaler += vtPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_BLOSUM)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aaBlosum[i][j];
                    a[j][i] = aaBlosum[j][i];
                    scaler += blosPi[i] * a[i][j] * probOn;
                    scaler += blosPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else if (aaModelID == AAMODEL_LG)
            {
            scaler = 0.0;
            for (i=0; i<20; i++)
                {
                for (j=i+1; j<20; j++)
                    {
                    a[i][j] = aaLG[i][j];
                    a[j][i] = aaLG[j][i];
                    scaler += lgPi[i] * a[i][j] * probOn;
                    scaler += lgPi[j] * a[j][i] * probOn;
                    }
                }
            }
        else
            {
            MrBayesPrint ("%s   ERROR: Don't understand which amino acid model is needed\n");
            return (ERROR);
            }

        /* rescale off diagonal elements of Q matrix */
        scaler = 1.0 / scaler;
        for (i=0; i<20; i++)
            {
            for (j=0; j<20; j++)
                {
                if (i != j)
                    a[i][j] *= scaler;
                }
            }
            
        /* now, scale by rate factor */
        for (i=0; i<20; i++)
            {
            for (j=0; j<20; j++)
                {
                if (i != j)
                    a[i][j] *= rateMult;
                }
            }
            
        /* put in diagonal elements */
        for (i=0; i<20; i++)
            {
            sum = 0.0;
            for (j=0; j<20; j++)
                {
                if (i != j)
                    sum += a[i][j];
                a[i][i] = -(sum + s10);
                }
            }
                
        /* fill in the other three submatrices */
        for (i=20; i<40; i++)
            a[i][i] = -s01;
        for (i=0; i<20; i++)
            {
            a[i][20+i] = s10;
            a[20+i][i] = s01;
            }
                       
        }

    return (NO_ERROR);
}


int SetStdQMatrix (MrBFlt **a, int nStates, MrBFlt *bs, int cType)
{
    register int    i, j;
    MrBFlt          scaler;

    /* This function sets up ordered or unordered models for standard characters
       with unequal stationary state frequencies. It requires the stationary
       frequencies of the states (passed when calling the function). It also
       needs to know the number of states and the type (ordered or unordered) 
       of the character. */

    /* set Q matrix to 0 */
    for (i=0; i<nStates; i++)
        for (j=0; j<nStates; j++)
            a[i][j] = 0.0;

    /* initialize Q matrix */
    scaler = 0.0;
    if (cType == UNORD)
        {
        /* unordered characters */
        for (i=0; i<nStates; i++)
            {
            for (j=0; j<nStates; j++)
                {
                if (i != j)
                    {
                    a[i][i] -= (a[i][j] = bs[j]);
                    scaler += bs[i] * a[i][j];
                    }
                }
            }
        }
    else
        {
        /* ordered characters */
        for (i=0; i<nStates; i++)
            {
            for (j=0; j<nStates; j++)
                {
                if (abs(i - j) == 1)
                    {
                    a[i][i] -= (a[i][j] = bs[j]);
                    scaler += bs[i] * a[i][j];
                    }
                }
            }
        }
        
    /* rescale Q matrix */
    for (i=0; i<nStates; i++)
        for (j=0; j<nStates; j++)
            a[i][j] /= scaler;

#   if defined DEBUG_SETSTDQMATRIX
    for (i=0; i<nStates; i++)
        {
        for (j=0; j<nStates; j++)
            printf ("%0.5lf ", a[i][j]);
        printf ("\n");
        }
#   endif

    return (NO_ERROR);
}


int TiProbs_Fels (TreeNode *p, int division, int chain)
{
    int         i, j, k, index;
    MrBFlt      t, u, x, z, beta, bigPi_j[4], pij, bigPij,
                *catRate, baseRate, theRate, *pis, length;
    CLFlt       *tiP;
    ModelInfo   *m;

    m = &modelSettings[division];

    /* find transition probabilities */
    tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];

    /* get base frequencies */
    pis = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
    
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;
    
    /* rescale beta */
    beta =  (0.5 / ((pis[0] + pis[2])*(pis[1] + pis[3]) + ((pis[0]*pis[2]) + (pis[1]*pis[3]))));

    bigPi_j[0] =  (pis[0] + pis[2]);
    bigPi_j[1] =  (pis[1] + pis[3]);
    bigPi_j[2] =  (pis[0] + pis[2]);
    bigPi_j[3] =  (pis[1] + pis[3]);

    /* find length */
    if (m->cppEvents != NULL)
        {
        length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
        }
    else if (m->tk02BranchRates != NULL)
        {
        length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
        }
    else if (m->igrBranchRates != NULL)
        {
        length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
        }
    else if (m->mixedBrchRates != NULL)
        {
        length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
        }
    else
        length = p->length;

    /* numerical errors will ensue if we allow very large or very small branch lengths,
       which might occur in relaxed clock models */

    /* fill in values */
    for (k=index=0; k<m->numRateCats; k++)
        {
        t =  length * baseRate * catRate[k];

        if (t < TIME_MIN)
            {
            /* Fill in identity matrix */
            for (i=0; i<4; i++)
                {
                for (j=0; j<4; j++)
                    {
                    if (i == j)
                        tiP[index++] = 1.0;
                    else
                        tiP[index++] = 0.0;
                    }
                }
            }
        else if (t > TIME_MAX)
            {
            /* Fill in stationary matrix */
            for (i=0; i<4; i++)
                for (j=0; j<4; j++)
                    tiP[index++] = (CLFlt) pis[j];
            }
        else
            {
            /* calculate probabilities */
            for (i=0; i<4; i++)
                {
                for (j=0; j<4; j++)
                    {
                    bigPij = bigPi_j[j];
                    pij =  pis[j];
                    u =  1.0/bigPij -  1.0;
                    x =  exp(-beta * t);
                    z = (bigPij - pij) / bigPij;
                    
                    if (i == j)
                        tiP[index++] = (CLFlt) (pij + pij * u * x + z * x);
                    else
                        tiP[index++] = (CLFlt) (pij + pij * u * x - (pij/bigPij) * x);
                    }
                }
            }
        }

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   TiProbs_Gen: Calculates transition probabilities for general
|       models with or without rate variation. This function does
|       not work with:
|      
|       1. codon models with omega variation or
|       2. covarion models with rate variation
|
|   In either of these cases, TiProbs_GenCov is used
|
-----------------------------------------------------------------*/
int TiProbs_Gen (TreeNode *p, int division, int chain)
{
    register int    i, j, k, n, s, index;
    MrBFlt          t, *catRate, baseRate, *eigenValues, *cijk, *bs,
                    EigValexp[64], sum, *ptr, theRate, correctionFactor,
                    length;
    CLFlt           *tiP;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    n = m->numModelStates;
    
    /* find the correction factor to make branch lengths
       in terms of expected number of substitutions per character */
    correctionFactor = 1.0;
    if (m->dataType == DNA || m->dataType == RNA)
        {
        if (m->nucModelId == NUCMODEL_DOUBLET)
            correctionFactor = 2.0;
        else if (m->nucModelId == NUCMODEL_CODON)
            correctionFactor = 3.0;
        }

    /* find transition probabilities */
    tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];
    
    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
        
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;

    /* get eigenvalues and cijk pointers */
    eigenValues = m->cijks[m->cijkIndex[chain]];
    cijk        = eigenValues + (2 * n);

    /* find length */
    if (m->cppEvents != NULL)
        {
        length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
        }
    else if (m->tk02BranchRates != NULL)
        {
        length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
        }
    else if (m->igrBranchRates != NULL)
        {
        length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
        }
    else if (m->mixedBrchRates != NULL)
        {
        length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
        }
    else
        length = p->length;

    /* fill in values */
    for (k=index=0; k<m->numRateCats; k++)
        {
        t =  length * baseRate * catRate[k] * correctionFactor;

        if (t < TIME_MIN)
            {
            /* Fill in identity matrix */
            for (i=0; i<n; i++)
                {
                for (j=0; j<n; j++)
                    {
                    if (i == j)
                        tiP[index++] = 1.0;
                    else
                        tiP[index++] = 0.0;
                    }
                }
            }
        else if (t > TIME_MAX)
            {
            /* Get base freq */
            bs = GetParamSubVals(m->stateFreq, chain, state[chain]);
            /* Fill in stationary matrix */
            for (i=0; i<n; i++)
                for (j=0; j<n; j++)
                    tiP[index++] = (CLFlt) bs[j];
            }
        else
            {
            /* We actually need to do some work... */
            for (s=0; s<n; s++)
                EigValexp[s] =  exp(eigenValues[s] * t);

            ptr = cijk;
            for (i=0; i<n; i++)
                {
                for (j=0; j<n; j++)
                    {
                    sum = 0.0;
                    for (s=0; s<n; s++)
                        sum += (*ptr++) * EigValexp[s];
                    tiP[index++] = (CLFlt) ((sum < 0.0) ? 0.0 : sum);
                    }
                }
            }
        }

#   if 0
    printf ("v = %lf (%d)\n", t, p->index);
    for (i=index=0; i<n; i++)
        {
        for (j=0; j<n; j++)
            printf ("%1.4lf ", tiP[index++]);
        printf ("\n");
        }
    printf ("\n");
#   endif

    return NO_ERROR;
}


/*----------------------------------------------------------------
|
|   TiProbs_GenCov: Calculates transition probabilities for codon
|       models with omega variation or covarion models with
|       rate variation.
|
-----------------------------------------------------------------*/
int TiProbs_GenCov (TreeNode *p, int division, int chain)
{
    register int    i, j, k, n, s, index;
    int             sizeOfSingleCijk;
    MrBFlt          t, *eigenValues, *cijk, EigValexp[64], sum, *ptr, correctionFactor,
                    length, *bs;
    CLFlt           *tiP;
    ModelInfo       *m;
    
    m = &modelSettings[division];
    n = m->numModelStates;
    
    /* find the correction factor to make branch lengths
       in terms of expected number of substitutions per character */
    correctionFactor = 1.0;
    if (m->dataType == DNA || m->dataType == RNA)
        {
        if (m->nucModelId == NUCMODEL_DOUBLET)
            correctionFactor = 2.0;
        else if (m->nucModelId == NUCMODEL_CODON)
            correctionFactor = 3.0;
        }

    /* find transition probabilities */
    tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];
            
    /* get eigenvalues and cijk pointers */
    eigenValues = m->cijks[m->cijkIndex[chain]];
    cijk        = eigenValues + (2 * n);
    
    /* get offset size (we need to move the pointers to the appropriate
       cijk information for these models) */
    sizeOfSingleCijk = m->cijkLength / m->nCijkParts;

    /* find length */
    if (m->cppEvents != NULL)
        {
        length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
        }
    else if (m->tk02BranchRates != NULL)
        {
        length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
        }
    else if (m->igrBranchRates != NULL)
        {
        length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
        }
    else if (m->mixedBrchRates != NULL)
        {
        length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
        }
    else
        length = p->length;

    /* numerical errors will ensue if we allow very large or very small branch lengths,
       which might occur in relaxed clock models */

    /* fill in values */
    for (k=index=0; k<m->nCijkParts; k++)
        {
        t =  length * correctionFactor;
            
        if (t < TIME_MIN)
            {
            /* Fill in identity matrix */
            for (i=0; i<n; i++)
                {
                for (j=0; j<n; j++)
                    {
                    if (i == j)
                        tiP[index++] = 1.0;
                    else
                        tiP[index++] = 0.0;
                    }
                }
            }
        else if (t > TIME_MAX)
            {
            /* Get base freq */
            bs = GetParamSubVals(m->stateFreq, chain, state[chain]);
            /* Fill in stationary matrix */
            for (i=0; i<n; i++)
                for (j=0; j<n; j++)
                    tiP[index++] = (CLFlt) bs[j];
            }
        else
            {
            /* We actually need to do some work... */
            for (s=0; s<n; s++)
                EigValexp[s] =  exp(eigenValues[s] * t);

            ptr = cijk;
            for (i=0; i<n; i++)
                {
                for (j=0; j<n; j++)
                    {
                    sum = 0.0;
                    for (s=0; s<n; s++)
                        sum += (*ptr++) * EigValexp[s];
                    tiP[index++] = (CLFlt) ((sum < 0.0) ? 0.0 : sum);
                    }
                }
                
            /* increment pointers by m->cijkLength */
            if (k+1 < m->nCijkParts)
                {
                /* shift pointers */
                eigenValues += sizeOfSingleCijk;
                cijk        += sizeOfSingleCijk;
                }
            }
        }
        
#   if 0
    for (i=index=0; i<n; i++)
        {
        for (j=0; j<n; j++)
            printf ("%1.4lf ", tiP[index++]);
        printf ("\n");
        }
#   endif

    return NO_ERROR;
}


/*-----------------------------------------------------------------
|
|   TiProbs_Hky: update transition probabilities for 4by4
|       nucleotide model with nst == 2 (K80/HKY85)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TiProbs_Hky (TreeNode *p, int division, int chain)
{
    int         i, j, k, index;
    MrBFlt      t, kap, u, w, x, y, z, beta, bigPi_j[4], pij, bigPij, *pis,
                *catRate, baseRate, theRate, length;
    CLFlt       *tiP;
    ModelInfo   *m;
    
    m = &modelSettings[division];

    /* find transition probabilities */
    tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];

    /* get kappa */
    kap =  *GetParamVals (m->tRatio, chain, state[chain]);
    
    /* get base frequencies */
    pis = GetParamSubVals (m->stateFreq, chain, state[chain]);
    
    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
    
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;
    
    /* rescale beta */
    beta =  0.5 / ((pis[0] + pis[2])*(pis[1] + pis[3]) + kap*((pis[0]*pis[2]) + (pis[1]*pis[3])));

    bigPi_j[0] = pis[0] + pis[2];
    bigPi_j[1] = pis[1] + pis[3];
    bigPi_j[2] = pis[0] + pis[2];
    bigPi_j[3] = pis[1] + pis[3];

    /* find length */
    if (m->cppEvents != NULL)
        {
        length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
        }
    else if (m->tk02BranchRates != NULL)
        {
        length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
        }
    else if (m->igrBranchRates != NULL)
        {
        length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
        }
    else if (m->mixedBrchRates != NULL)
        {
        length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
        }
    else
        length = p->length;

    /* numerical errors will ensue if we allow very large or very small branch lengths,
       which might occur in relaxed clock models */

    /* fill in values */
    for (k=index=0; k<m->numRateCats; k++)
        {
        t =  length * baseRate * catRate[k];

        if (t < TIME_MIN)
            {
            /* Fill in identity matrix */
            for (i=0; i<4; i++)
                {
                for (j=0; j<4; j++)
                    {
                    if (i == j)
                        tiP[index++] = 1.0;
                    else
                        tiP[index++] = 0.0;
                    }
                }
            }
        else if (t > TIME_MAX)
            {
            /* Fill in stationary matrix */
            for (i=0; i<4; i++)
                for (j=0; j<4; j++)
                    tiP[index++] = (CLFlt) pis[j];
            }
        else
            {
            /* calculate probabilities */
            for (i=0; i<4; i++)
                {
                for (j=0; j<4; j++)
                    {
                    bigPij = bigPi_j[j];
                    pij = pis[j];
                    u =  1.0/bigPij -  1.0;
                    w = -beta * (1.0 + bigPij * (kap -  1.0));
                    x =  exp(-beta * t);
                    y =  exp(w * t);
                    z = (bigPij - pij) / bigPij;
                    
                    if (i == j)
                        tiP[index++] = (CLFlt) (pij + pij * u * x + z * y);
                    else if ((i == 0 && j == 2) || (i == 2 && j == 0) || (i == 1 && j == 3) || (i == 3 && j == 1))
                        tiP[index++] = (CLFlt) (pij + pij * u * x - (pij/bigPij) * y);
                    else
                        tiP[index++] = (CLFlt) (pij * (1.0 - x));
                    }
                }
            }
        }
        
    return NO_ERROR;
}


/*-----------------------------------------------------------------
|
|   TiProbs_JukesCantor: update transition probabilities for 4by4
|       nucleotide model with nst == 1 (Jukes-Cantor)
|       with or without rate variation
|
------------------------------------------------------------------*/
int TiProbs_JukesCantor (TreeNode *p, int division, int chain)
{
    /* calculate Jukes Cantor transition probabilities */
    
    int         i, j, k, index;
    MrBFlt      t, *catRate, baseRate, theRate, length;
    CLFlt       pNoChange, pChange, *tiP;
    ModelInfo   *m;
    
    m = &modelSettings[division];

    /* find transition probabilities */
    tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];

    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* compensate for invariable sites if appropriate */
    if (m->pInvar != NULL)
        baseRate /= (1.0 - (*GetParamVals(m->pInvar, chain, state[chain])));
    
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;

    /* find length */
    if (m->cppEvents != NULL)
        {
        length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
        }
    else if (m->tk02BranchRates != NULL)
        {
        length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
        }
    else if (m->igrBranchRates != NULL)
        {
        length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
        }
    else if (m->mixedBrchRates != NULL)
        {
        length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
        }
    else
        length = p->length;

    /* numerical errors will ensue if we allow very large or very small branch lengths,
       which might occur in relaxed clock models */

    /* fill in values */
    for (k=index=0; k<m->numRateCats; k++)
        {
        t = length * baseRate * catRate[k];
            
        if (t < TIME_MIN)
            {
            /* Fill in identity matrix */
            for (i=0; i<4; i++)
                {
                for (j=0; j<4; j++)
                    {
                    if (i == j)
                        tiP[index++] = 1.0;
                    else
                        tiP[index++] = 0.0;
                    }
                }
            }
        else if (t > TIME_MAX)
            {
            /* Fill in stationary matrix */
            for (i=0; i<4; i++)
                for (j=0; j<4; j++)
                    tiP[index++] = 0.25;
            }
        else
            {
            /* calculate probabilities */
            pChange   = (CLFlt) (0.25 - 0.25 * exp(-(4.0/3.0)*t));
            pNoChange = (CLFlt) (0.25 + 0.75 * exp(-(4.0/3.0)*t));
            for (i=0; i<4; i++)
                {
                for (j=0; j<4; j++)
                    {
                    if (i == j)
                        tiP[index++] = pNoChange;
                    else
                        tiP[index++] = pChange;
                    }
                }
            }
        }

    return NO_ERROR;
}


/*-----------------------------------------------------------------
|
|   TiProbs_Res: update transition probabilities for binary
|       restriction site model with or without rate variation
|
------------------------------------------------------------------*/
int TiProbs_Res (TreeNode *p, int division, int chain)
{
    int         k, index;
    MrBFlt      baseRate, eV, mu, theRate, v,
                *bs, *catRate, length;
    CLFlt       *tiP;
    ModelInfo   *m;
    
    /* find model settings for the division */
    m = &modelSettings[division];

    /* find transition probabilities */
    tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];

    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;

    /* find base frequencies */
    bs = GetParamSubVals(m->stateFreq, chain, state[chain]);

    /* calculate scaling factor */
    mu =  1.0 / (2.0 * bs[0] * bs[1]);
    
    /* find length */
    if (m->cppEvents != NULL)
        {
        length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
        }
    else if (m->tk02BranchRates != NULL)
        {
        length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
        }
    else if (m->igrBranchRates != NULL)
        {
        length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
        }
    else if (m->mixedBrchRates != NULL)
        {
        length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
        }
    else
        length = p->length;

    /* numerical errors will ensue if we allow very large or very small branch lengths,
       which might occur in relaxed clock models */

    /* fill in values */
    for (k=index=0; k<m->numRateCats; k++)
        {       
        v =  length * baseRate * catRate[k];
            
        if (v < TIME_MIN)
            {
            /* Fill in identity matrix */
            tiP[index++] = (CLFlt) (bs[0] + bs[1]);
            tiP[index++] = (CLFlt) (bs[1] - bs[1]);
            tiP[index++] = (CLFlt) (bs[0] - bs[0]);
            tiP[index++] = (CLFlt) (bs[1] + bs[0]);
            }
        else if (v > TIME_MAX)
            {
            /* Fill in stationary matrix */
            tiP[index++] = (CLFlt) bs[0];
            tiP[index++] = (CLFlt) bs[1];
            tiP[index++] = (CLFlt) bs[0];
            tiP[index++] = (CLFlt) bs[1];
            }
        else
            {
            /* calculate probabilities */
            eV =  exp(-mu * v);
            tiP[index++] = (CLFlt) (bs[0] + bs[1] * eV);
            tiP[index++] = (CLFlt) (bs[1] - bs[1] * eV);
            tiP[index++] = (CLFlt) (bs[0] - bs[0] * eV);
            tiP[index++] = (CLFlt) (bs[1] + bs[0] * eV);
            }
        }

    return NO_ERROR;
}


/*-----------------------------------------------------------------
|
|   TiProbs_Std: update transition probabilities for
|       variable states model with or without rate variation
|
------------------------------------------------------------------*/
int TiProbs_Std (TreeNode *p, int division, int chain)
{
    int         b, c, i, j, k, n, s, nStates, index=0, index2;
    MrBFlt      v, eV1, eV2, eV3, eV4, eV5, *catRate,
                baseRate, theRate, pi, f1, f2, f3, f4, f5, f6, f7, root,
                *eigenValues, *cijk, sum, *bs, mu, length;
    CLFlt       pNoChange, pChange, *tiP;
    ModelInfo   *m;
#   if defined (DEBUG_TIPROBS_STD)
    int         index3;
#   endif

    m = &modelSettings[division];

    /* find transition probabilities */
    tiP = m->tiProbs[m->tiProbsIndex[chain][p->index]];
    
    /* get base rate */
    baseRate = GetRate (division, chain);
    
    /* get category rates */
    theRate = 1.0;
    if (m->shape != NULL)
        catRate = GetParamSubVals (m->shape, chain, state[chain]);
    else if (m->mixtureRates != NULL)
        catRate = GetParamSubVals (m->mixtureRates, chain, state[chain]);
    else
        catRate = &theRate;
    
#   if defined (DEBUG_TIPROBS_STD)
    /* find base frequencies */
    bs = GetParamStdStateFreqs (m->stateFreq, chain, state[chain]);
#   endif

    /* find length */
    if (m->cppEvents != NULL)
        {
        length = GetParamSubVals (m->cppEvents, chain, state[chain])[p->index];
        }
    else if (m->tk02BranchRates != NULL)
        {
        length = GetParamSubVals (m->tk02BranchRates, chain, state[chain])[p->index];
        }
    else if (m->igrBranchRates != NULL)
        {
        length = GetParamSubVals (m->igrBranchRates, chain, state[chain])[p->index];
        }
    else if (m->mixedBrchRates != NULL)
        {
        length = GetParamSubVals (m->mixedBrchRates, chain, state[chain])[p->index];
        }
    else
        length = p->length;

    /* numerical errors will ensue if we allow very large or very small branch lengths, which might
       occur in relaxed clock models; an elegant solution would be to substitute the stationary
       probs and initial probs but for now we truncate lengths at small or large values TODO */
    if (length > BRLENS_MAX)
        length = BRLENS_MAX;
    else if (length < BRLENS_MIN)
        length = BRLENS_MIN;

    /* fill in values; this has to be done differently if state freqs are not equal */
    if (m->stateFreq->paramId == SYMPI_EQUAL)
        {
        /* equal state frequencies */
        /* fill in values for unordered characters */
        index = 0;
#   if defined (DEBUG_TIPROBS_STD)
        index3 = 0;
#   endif
        for (nStates=2; nStates<=10; nStates++)
            {
            if (m->isTiNeeded[nStates-2] == NO)
                continue;
            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                eV1 =  exp(-(nStates / (nStates -  1.0)) * v);
                pChange   = (CLFlt) ((1.0 / nStates) - ((1.0 / nStates) * eV1));
                pNoChange = (CLFlt) ((1.0 / nStates) + ((nStates - 1.0) / nStates) * eV1);
                if (pChange<0.0)
                    pChange = (CLFlt) 0.0;
                for (i=0; i<nStates; i++)
                    {
                    for (j=0; j<nStates; j++)
                        {
                        if (i == j)
                            tiP[index++] = pNoChange;
                        else
                            tiP[index++] = pChange;
                        }
                    }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }
#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }

        /* fill in values for 3-state ordered character */
        if (m->isTiNeeded[9] == YES)
            {
            nStates = 3;
            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                eV1 =  exp (-(3.0 / 4.0) * v);
                eV2 =  exp (-(9.0 / 4.0) * v);
                
                /* pij(0,0) */
                tiP[index] = (CLFlt) ((1.0 / 3.0) + (eV1 / 2.0) + (eV2 / 6.0));
                /* pij(0,1) = pij(1,0) */
                tiP[index+1] = tiP[index+3] = (CLFlt) ((1.0 / 3.0) - (eV2 / 3.0));
                /* pij(0,2) */
                tiP[index+2] = (CLFlt) ((1.0 / 3.0) - (eV1 / 2.0) + (eV2 / 6.0));
                /* pij(1,1) */
                tiP[index+4] = (CLFlt) ((1.0 / 3.0) + (2.0 * eV2 / 3.0));
                
                /* fill in mirror part of matrix */
                index += 5;
                index2 = index - 2;
                for (i=0; i<4; i++)
                    tiP[index++] = tiP[index2--];

                /* make sure no value is negative */
                for (i=index-(nStates*nStates); i<index; i++) {
                    if (tiP[i] < 0.0)
                        tiP[i] = (CLFlt) 0.0;
                }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }

#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }

        /* 4-state ordered character */
        if (m->isTiNeeded[10] == YES)
            {
            nStates = 4;
            pi = 1.0 / 4.0;
            root =  sqrt (2.0);
            f1 = root +  1.0;
            f2 = root -  1.0;

            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                eV1 =  1.0 / (exp ((4.0 * v) / 3.0));
                eV2 =  exp ((2.0 * (root - 2.0) * v) / 3.0) / root;
                eV3 =  1.0 / (root *  exp ((2.0 * (root + 2.0) * v) / 3.0));
                
                /* pij(0,0) */
                tiP[index] = (CLFlt) (pi * (1.0 + eV1 + (f1*eV2) + (f2*eV3)));
                /* pij(0,1) = pij(1,0) */
                tiP[index+1] = tiP[index+4] = (CLFlt) (pi * (1.0 - eV1 + eV2 - eV3));
                /* pij(0,2) = tiP(1,3) */
                tiP[index+2] = tiP[index+7] = (CLFlt) (pi * (1.0 - eV1 - eV2 + eV3));
                /* pij(0,3) */
                tiP[index+3] = (CLFlt) (pi * (1.0 + eV1 - (f1*eV2) - (f2*eV3)));
                /* pij(1,1) */
                tiP[index+5] = (CLFlt) (pi * (1.0 + eV1 + (f2*eV2) + (f1*eV3)));
                /* pij(1,2) */
                tiP[index+6] = (CLFlt) (pi * (1.0 + eV1 - (f2*eV2) - (f1*eV3)));

                /* fill in mirror part of matrix */
                index += 8;
                index2 = index - 1;
                for (i=0; i<8; i++)
                    tiP[index++] = tiP[index2--];
        
                /* make sure no value is negative */
                for (i=index-(nStates*nStates); i<index; i++) {
                    if (tiP[i] < 0.0)
                        tiP[i] = (CLFlt) 0.0;
                }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }
#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }

        /* 5-state ordered character */
        if (m->isTiNeeded[11] == YES)
            {
            nStates = 5;
            pi = 1.0 / 5.0;
            root =  sqrt (5.0);

            f5 = root /  4.0;
            f1 =  0.75 + f5;;
            f2 =  1.25 + f5;
            f3 =  1.25 - f5;
            f4 =  0.75 - f5;
            f5 = f5 *  2.0;
            f6 = f5 +  0.5;
            f7 = f5 -  0.5;

            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                v *=  5.0 /  16.0;

                eV1 =  exp ((root -  3.0) * v);
                eV2 =  exp (-(root +  3.0) * v);
                eV3 =  exp ((root -  5.0) * v);
                eV4 =  exp (-(root +  5.0) * v);

                /* pij(0,0) */
                tiP[index] = (CLFlt) (pi* (1.0 + (f1*eV3) + (f2*eV1) + (f3*eV2) + (f4*eV4)));
                /* pij(0,1) = pij(1,0) */
                tiP[index+1] = tiP[index+5] =
                    (CLFlt) (pi*(1.0 - (eV3/2.0) + (f5*eV1) - (f5*eV2) - (eV4/2.0)));
                /* pij(0,2) = pij(2,0) */
                tiP[index+2] = tiP[index+10] = (CLFlt) (pi*(1.0 - (f6*eV3) + (f7*eV4)));
                /* pij(0,3) = pij(1,4) */
                tiP[index+3] = tiP[index+9] =
                    (CLFlt) (pi*(1.0 - (eV3/2.0) - (f5*eV1) + (f5*eV2) - (eV4/2.0)));
                /* pij(0,4) */
                tiP[index+4] = (CLFlt) (pi*(1.0 + (f1*eV3) - (f2*eV1) - (f3*eV2) + (f4*eV4)));
                /* pij(1,1) */
                tiP[index+6] = (CLFlt) (pi*(1.0 + (f4*eV3) + (f3*eV1) + (f2*eV2) + (f1*eV4)));
                /* pij(1,2) = pij(2,1) */
                tiP[index+7] = tiP[index+11] = (CLFlt) (pi*(1.0 + (f7*eV3) - (f6*eV4)));
                /* pij(1,3) */
                tiP[index+8] = (CLFlt) (pi*(1.0 + (f4*eV3) - (f3*eV1) - (f2*eV2) + (f1*eV4)));
                /* pij(2,2) */
                tiP[index+12] = (CLFlt) (pi*(1.0 + (2.0*eV3) + (2.0*eV4)));

                /* fill in mirror part of matrix */
                index += 13;
                index2 = index - 2;
                for (i=0; i<12; i++)
                    tiP[index++] = tiP[index2--];

                /* make sure no value is negative */
                for (i=index-(nStates*nStates); i<index; i++) {
                    if (tiP[i] < 0.0)
                        tiP[i] = (CLFlt) 0.0;
                }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }
#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }

        /* 6-state ordered character */
        if (m->isTiNeeded[12] == YES)
            {
            nStates = 6;
            pi =  1.0 /  6.0;
            root =  sqrt (3.0);

            f4 = (3.0 / (2.0 * root));
            f1 =  1.0 + f4;
            f2 =  1.0 - f4;
            f3 =  0.5 + f4;
            f4 =  0.5 - f4;

            for (k=0; k<m->numRateCats; k++)
                {
                /* calculate probabilities */
                v =  length * catRate[k] * baseRate;
                v /=  5.0;

                eV1 =  exp (-9 * v);
                eV2 =  exp (-6 * v);
                eV3 =  exp (-3 * v);
                eV4 =  exp (3.0 * (root - 2.0) * v);
                eV5 =  exp (-3.0 * (root + 2.0) * v);

                /* pij(0,0) */
                tiP[index] = (CLFlt) (pi* (1.0 + (0.5*eV1) + eV2 + (1.5*eV3) + (f1*eV4) + (f2*eV5)));
                /* pij(0,1) = pij(1,0) */
                tiP[index+1] = tiP[index+6] = (CLFlt) (pi*(1.0 - eV1 - eV2 + (f3*eV4) + (f4*eV5)));
                /* pij(0,2) = pij(2,0) */
                tiP[index+2] = tiP[index+12] = 
                    (CLFlt) (pi*(1.0 + (0.5*eV1) - eV2 - (1.5*eV3) + (0.5*eV4) + (0.5*eV5)));
                /* pij(0,3) = pij(2,5) */
                tiP[index+3] = tiP[index+17] = 
                    (CLFlt) (pi*(1.0 + (0.5*eV1) + eV2 - (1.5*eV3) - (0.5*eV4) - (0.5*eV5)));
                /* pij(0,4) = pij(1,5) */
                tiP[index+4] = tiP[index+11] = (CLFlt) (pi*(1.0 - eV1 + eV2 - (f3*eV4) - (f4*eV5)));
                /* pij(0,5) */
                tiP[index+5] = (CLFlt) (pi*(1.0 + (0.5*eV1) - eV2 + (1.5*eV3) - (f1*eV4) - (f2*eV5)));
                /* pij(1,1) */
                tiP[index+7] = (CLFlt) (pi*(1.0 + (2.0*eV1) + eV2 + eV4 + eV5));
                /* pij(1,2) = pij(2,1) */
                tiP[index+8] = tiP[index+13] = (CLFlt) (pi*(1.0 - eV1 + eV2 - (f4*eV4) - (f3*eV5)));
                /* pij(1,3) = pij(2,4) */
                tiP[index+9] = tiP[index+16] = (CLFlt) (pi*(1.0 - eV1 - eV2 + (f4*eV4) + (f3*eV5)));
                /* pij(1,4) */
                tiP[index+10] = (CLFlt) (pi*(1.0 + (2.0*eV1) - eV2 - eV4 - eV5));
                /* pij(2,2) */
                tiP[index+14] = (CLFlt) (pi*(1.0 + (0.5*eV1) + eV2 + (1.5*eV3) + (f2*eV4) + (f1*eV5)));
                /* pij(2,3) */
                tiP[index+15] = (CLFlt) (pi*(1.0 + (0.5*eV1) - eV2 + (1.5*eV3) - (f2*eV4) - (f1*eV5)));

                /* fill in mirror part of matrix */
                index += 18;
                index2 = index - 1;
                for (i=0; i<18; i++)
                    tiP[index++] = tiP[index2--];

                /* make sure no value is negative */
                for (i=index-(nStates*nStates); i<index; i++) {
                    if (tiP[i] < 0.0)
                        tiP[i] = (CLFlt) 0.0;
                }
#   if defined (DEBUG_TIPROBS_STD)
                PrintTiProbs (tiP+index-(nStates*nStates), bs+index3, nStates);
#   endif
                }
#   if defined (DEBUG_TIPROBS_STD)
            index3 += nStates;
#   endif
            }
        }
    else
        {
        /* unequal state frequencies */
        index = 0;

        /* first fill in for binary characters using beta categories if needed */
        if (m->isTiNeeded[0] == YES)
            {
            /* find base frequencies */
            bs = GetParamStdStateFreqs (m->stateFreq, chain, state[chain]);

            /* cycle through beta and gamma cats */
            for (b=0; b<m->numBetaCats; b++)
                {
                mu =  1.0 / (2.0 * bs[0] * bs[1]);
                for (k=0; k<m->numRateCats; k++)
                    {
                    /* calculate probabilities */
                    v =  length * catRate[k] * baseRate;
                    eV1 =  exp(- mu * v);
                    tiP[index++] = (CLFlt) (bs[0] + (bs[1] * eV1));
                    tiP[index++] = (CLFlt) (bs[1] - (bs[1] * eV1));
                    tiP[index++] = (CLFlt) (bs[0] - (bs[0] * eV1));
                    tiP[index++] = (CLFlt) (bs[1] + (bs[0] * eV1));
                    }
                /* update stationary state frequency pointer */
                bs += 2;
                }
            }

        /* now use general algorithm for the other cases */
        if (m->cijkLength > 0)
            {
            /* first update cijk if necessary */
            if (m->cijkLength > 0 && m->upDateCijk == YES)
                {
                if (UpDateCijk (division, chain) == ERROR)
                    return (ERROR);
                }

            /* then get first set of eigenvalues */
            eigenValues = m->cijks[m->cijkIndex[chain]];

            /* and cycle through the relevant characters */
            for (c=0; c<m->stateFreq->nSympi; c++)
                {
                n = m->stateFreq->sympinStates[c];

                /* fill in values */
                for (k=0; k<m->numRateCats; k++)
                    {
                    v =  length * baseRate * catRate[k];
                    cijk = eigenValues + (2 * n);

                    for (i=0; i<n; i++)
                        {
                        for (j=0; j<n; j++)
                            {
                            sum = 0.0;
                            for (s=0; s<n; s++)
                                sum += (*cijk++) * exp(eigenValues[s] * v);
                            tiP[index++] = (CLFlt) ((sum <  0.0) ?  0.0 : sum);
                            }
                        }
                    }

                /* update eigenValues pointer */
                eigenValues += (n * n * n) + (2 * n);
                }
            }
        }

    return NO_ERROR;
}


int UpDateCijk (int whichPart, int whichChain)
{
    int         c, i, j, k, n, n3, isComplex, sizeOfSingleCijk, cType, numQAllocated;
    MrBFlt      **q[100], **eigvecs, **inverseEigvecs;
    MrBFlt      *eigenValues, *eigvalsImag, *cijk;
    MrBFlt      *bs, *bsBase, *rateOmegaValues=NULL, rA=0.0, rS=0.0, posScaler, *omegaCatFreq=NULL;
    MrBComplex     **Ceigvecs, **CinverseEigvecs;
    ModelInfo   *m;
    Param       *p;
#   if defined (BEAGLE_ENABLED)
    int         u;
    double      *beagleEigvecs=NULL, *beagleInverseEigvecs=NULL;
#   endif

    /* get a pointer to the model settings for this partition */
    m = &modelSettings[whichPart];
    assert (m->upDateCijk == YES);
    
    /* we should only go through here if we have cijk information available for the partition */
    if (m->cijkLength > 0) 
        {
        /* flip cijk space */
        FlipCijkSpace(m, whichChain);
        
        /* figure out information on either omega values or rate values, if necessary */
        if (m->dataType == DNA || m->dataType == RNA)
            {
            if (m->nucModelId == NUCMODEL_CODON)                                                    /* we have a NY98 model     */
                {
                rateOmegaValues = GetParamVals(m->omega, whichChain, state[whichChain]);
                if (m->numOmegaCats > 1)
                    omegaCatFreq = GetParamSubVals (m->omega, whichChain, state[whichChain]);
                }
            else if (m->nCijkParts > 1 && m->nucModelId == NUCMODEL_4BY4 && m->numModelStates == 8)
                {
                /* we have a covarion (covariotide) model with rate variation */
                if (m->shape != NULL)
                    rateOmegaValues = GetParamSubVals (m->shape, whichChain, state[whichChain]);
                else if (m->mixtureRates != NULL)
                    rateOmegaValues = GetParamSubVals (m->mixtureRates, whichChain, state[whichChain]);
                }
            }
        else if (m->dataType == PROTEIN)
            {
            if (m->nCijkParts > 1)
                {
                /* we have a covarion model with rate variation */
                if (m->shape != NULL)
                    rateOmegaValues = GetParamSubVals (m->shape, whichChain, state[whichChain]);
                else if (m->mixtureRates != NULL)
                    rateOmegaValues = GetParamSubVals (m->mixtureRates, whichChain, state[whichChain]);
                }
            }
#   if defined (BEAGLE_ENABLED)
        else if (m->dataType == RESTRICTION){}
#   endif
        else if (m->dataType != STANDARD)
            {
            MrBayesPrint ("%s   ERROR: Should not be updating cijks!\n", spacer);
            return (ERROR);
            }
        
        if (m->dataType == STANDARD)
            {
            /* set pointers and other stuff needed */
            numQAllocated = 1;
            p = m->stateFreq;
            eigenValues = m->cijks[m->cijkIndex[whichChain]];
            q[0] = AllocateSquareDoubleMatrix (10);
            eigvecs = AllocateSquareDoubleMatrix (10);
            inverseEigvecs = AllocateSquareDoubleMatrix (10);
            Ceigvecs = AllocateSquareComplexMatrix (10);
            CinverseEigvecs = AllocateSquareComplexMatrix (10);
            bsBase = GetParamStdStateFreqs (m->stateFreq, whichChain, state[whichChain]);
            
            /* cycle over characters needing cijks */
            for (c=0; c<p->nSympi; c++)
                {
                n = p->sympinStates[c];
                bs = bsBase + p->sympiBsIndex[c];
                cType = p->sympiCType[c];
                n3 = n * n * n;
                eigvalsImag = eigenValues + n;
                cijk = eigenValues + (2 * n);
                if (SetStdQMatrix (q[0], n, bs, cType) == ERROR)
                    return (ERROR);
                isComplex = GetEigens (n, q[0], eigenValues, eigvalsImag, eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
                if (isComplex == NO)
                    {
                    CalcCijk (n, cijk, eigvecs, inverseEigvecs);
                    }
                else
                    {
                    if (isComplex == YES)
                        MrBayesPrint ("%s   ERROR: Complex eigenvalues found!\n", spacer);
                    else
                        MrBayesPrint ("%s   ERROR: Computing eigenvalues problem!\n", spacer);
                    goto errorExit;
                    }
                eigenValues += (n3 + (2 * n));
                }
            }
        else
            {
            /* all other data types */
            numQAllocated = m->nCijkParts;
            sizeOfSingleCijk = m->cijkLength / m->nCijkParts;
            n = m->numModelStates;
#   if defined (BEAGLE_ENABLED)
            if (m->useBeagle == YES)
                eigenValues = m->cijks[m->cijkIndex[whichChain]/m->nCijkParts];
            else
                eigenValues = m->cijks[m->cijkIndex[whichChain]];
#   else
            eigenValues = m->cijks[m->cijkIndex[whichChain]];
#   endif
            eigvalsImag = eigenValues + n;
            cijk        = eigenValues + (2 * n);
            for (k=0; k<numQAllocated; k++)
                q[k] = AllocateSquareDoubleMatrix (n);
            eigvecs = AllocateSquareDoubleMatrix (n);
            inverseEigvecs = AllocateSquareDoubleMatrix (n);
            Ceigvecs = AllocateSquareComplexMatrix (n);
            CinverseEigvecs = AllocateSquareComplexMatrix (n);
            
            if (m->nCijkParts == 1)
                {
                if (m->dataType == DNA || m->dataType == RNA)
                    {
                    if (m->nucModelId == NUCMODEL_CODON)
                        {
                        if (SetNucQMatrix (q[0], n, whichChain, whichPart, rateOmegaValues[0], &rA, &rS) == ERROR)
                            goto errorExit;
                        }
                    else
                        {
                        if (SetNucQMatrix (q[0], n, whichChain, whichPart, 1.0, &rA, &rS) == ERROR)
                            goto errorExit;
                        }
                    }
#   if defined (BEAGLE_ENABLED)
                else if (m->dataType == RESTRICTION)
                    {
                    SetBinaryQMatrix (q[0], whichChain, whichPart);
                    }
#   endif
                else
                    {
                    if (SetProteinQMatrix (q[0], n, whichChain, whichPart, 1.0) == ERROR)
                        goto errorExit;
                    }
                isComplex = GetEigens (n, q[0], eigenValues, eigvalsImag, eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
#   if defined (BEAGLE_ENABLED)
                if (isComplex == YES)
                    {
                    if (isComplex == YES)
                        MrBayesPrint ("%s   ERROR: Complex eigenvalues found!\n", spacer);
                    else
                        MrBayesPrint ("%s   ERROR: Computing eigenvalues problem!\n", spacer);
                    goto errorExit;
                    }
                if (m->useBeagle == YES)
                    {
                    /* TODO: only allocate this space once at initialization */
                    beagleEigvecs = (double*) SafeCalloc (2*n*n, sizeof(double));
                    beagleInverseEigvecs = beagleEigvecs + n*n;
                    for (i=k=0; i<n; i++)
                        {
                        // eigenValues[i] = 0.1;
                        for (j=0; j<n; j++)
                            {
                            beagleEigvecs[k] = eigvecs[i][j];
                            beagleInverseEigvecs[k] = inverseEigvecs[i][j];
                            k++;
                            }
                        }
                    beagleSetEigenDecomposition(m->beagleInstance,
                                                m->cijkIndex[whichChain],
                                                beagleEigvecs,
                                                beagleInverseEigvecs,
                                                eigenValues);
                    free(beagleEigvecs);
                    }
                else
                    {
                    CalcCijk (n, cijk, eigvecs, inverseEigvecs);
                    }
#   else
                if (isComplex == NO)
                    {
                    CalcCijk (n, cijk, eigvecs, inverseEigvecs);
                    }
                else
                    {
                    MrBayesPrint ("%s   ERROR: Complex eigenvalues found!\n", spacer);
                    goto errorExit;
                    }
#   endif
                }
            else
                {
                /* Here, we calculate the rate matrices (Q) for various nucleotide and amino acid
                   data models. Usually, when the rate matrix is set in SetNucQMatrix, it is scaled
                   such that the average substitution rate is one. However, there is a complication
                   for positive selection models using codon rate matrices. First, we have more than
                   one matrix; in fact, we have as many rate matrices as there are omega values. Second,
                   the mean substitution rate still has to be one. And third, we want the synonymous
                   rate to be the same across the rate matrices. For positive selection models, the Q
                   matrix comes out of SetNucQMatrix unscaled. Once we have all m->nCijkParts rate 
                   matrices, we then scale again, this time to ensure that the mean substitution rate is one. */

                /* First, calculate rate matrices for each category: */
                posScaler = 0.0;
                for (k=0; k<m->nCijkParts; k++)
                    {
                    if (m->dataType == DNA || m->dataType == RNA)
                        {
                        if (SetNucQMatrix (q[k], n, whichChain, whichPart, rateOmegaValues[k], &rA, &rS) == ERROR)
                            goto errorExit;
                        }
                    else
                        {
                        if (SetProteinQMatrix (q[k], n, whichChain, whichPart, rateOmegaValues[k]) == ERROR)
                            goto errorExit;
                        }
                    if (m->nucModelId == NUCMODEL_CODON && m->numOmegaCats > 1)
                        posScaler += omegaCatFreq[k] * (rS + rA);
                    }
                    
                /* Then rescale the rate matrices, if this is a positive selection model: */
                if (m->nucModelId == NUCMODEL_CODON && m->numOmegaCats > 1)
                    {
                    posScaler = 1.0 / posScaler;
                    for (k=0; k<m->nCijkParts; k++)
                        {
                        for (i=0; i<n; i++)
                            for (j=0; j<n; j++)
                                q[k][i][j] *= posScaler;
                        }
                    }

                /* Finally, calculate eigenvalues, etc.: */
#   if defined (BEAGLE_ENABLED)
                if (m->useBeagle == YES)
                    {
                    /* TODO: only allocate this space once at initialization */
                    beagleEigvecs = (double*) SafeCalloc (2*n*n, sizeof(double));
                    beagleInverseEigvecs = beagleEigvecs + n*n;
                    }
#   endif
                for (k=0; k<m->nCijkParts; k++)
                    {
                    isComplex = GetEigens (n, q[k], eigenValues, eigvalsImag, eigvecs, inverseEigvecs, Ceigvecs, CinverseEigvecs);
#   if defined (BEAGLE_ENABLED)
                    if (isComplex == YES)
                        {
                        if (isComplex == YES)
                            MrBayesPrint ("%s   ERROR: Complex eigenvalues found!\n", spacer);
                        else
                            MrBayesPrint ("%s   ERROR: Computing eigenvalues problem!\n", spacer);
                        goto errorExit;
                        }
                    if (m->useBeagle == YES)
                        {
                        for (i=u=0; i<n; i++)
                            {
                            for (j=0; j<n; j++)
                                {
                                beagleEigvecs[u] = eigvecs[i][j];
                                beagleInverseEigvecs[u] = inverseEigvecs[i][j];
                                u++;
                                }
                            }

                        beagleSetEigenDecomposition(m->beagleInstance,
                                                    m->cijkIndex[whichChain] + k,
                                                    beagleEigvecs,
                                                    beagleInverseEigvecs,
                                                    eigenValues);
                        }
                    else
                        {
                        CalcCijk (n, cijk, eigvecs, inverseEigvecs);
                        }
#   else
                    if (isComplex == NO)
                        {
                        CalcCijk (n, cijk, eigvecs, inverseEigvecs);
                        }
                    else
                        {
                        MrBayesPrint ("%s   ERROR: Complex eigenvalues found!\n", spacer);
                        goto errorExit;
                        }
#   endif
                    /* shift pointers */
                    eigenValues += sizeOfSingleCijk;
                    eigvalsImag += sizeOfSingleCijk;
                    cijk        += sizeOfSingleCijk;
                    }
#   if defined (BEAGLE_ENABLED)
                free(beagleEigvecs);
#   endif
                }
            }
            
        for (k=0; k<numQAllocated; k++)
            FreeSquareDoubleMatrix (q[k]);
        FreeSquareDoubleMatrix (eigvecs);
        FreeSquareDoubleMatrix (inverseEigvecs);
        FreeSquareComplexMatrix (Ceigvecs);
        FreeSquareComplexMatrix (CinverseEigvecs);
        }
        
    return (NO_ERROR);

    errorExit:      
        for (k=0; k<numQAllocated; k++)
            FreeSquareDoubleMatrix (q[k]);
        FreeSquareDoubleMatrix (eigvecs);
        FreeSquareDoubleMatrix (inverseEigvecs);
        FreeSquareComplexMatrix (Ceigvecs);
        FreeSquareComplexMatrix (CinverseEigvecs);

        return ERROR;
}

