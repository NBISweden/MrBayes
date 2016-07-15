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
#include "best.h"
#include "command.h"
#include "mcmc.h"
#include "model.h"
#include "proposal.h"
#include "utils.h"

/* debugging compiler statements */
#undef  DEBUG_LOCAL
#undef  DEBUG_UNROOTED_SLIDER
#undef  DEBUG_ParsSPR
#undef  DEBUG_ExtSS
#undef  DEBUG_CSLIDER
#undef  DEBUG_ExtSPRClock
#undef  DEBUG_ParsSPRClock
#undef  DEBUG_ExtTBR
#undef  DEBUG_NNIClock
#undef  DEBUG_SPLITMERGE
#undef  DEBUG_FBDPR


extern int *chainId;

void TouchAllTreeNodes (ModelInfo *m, int chain);


int Move_Aamodel (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change amino acid model for model mixing 
       amino acid model ID's
        AAMODEL_POISSON         0
        AAMODEL_JONES           1
        AAMODEL_DAY             2
        AAMODEL_MTREV           3
        AAMODEL_MTMAM           4
        AAMODEL_WAG             5
        AAMODEL_RTREV           6 
        AAMODEL_CPREV           7 
        AAMODEL_VT              8
        AAMODEL_BLOSUM          9 */

    int         i, oldM, newM;
    MrBFlt      *bs, *subValue;
    ModelParams *mp;
    
    /* get model params */
    mp = &modelParams[param->relParts[0]];

    subValue = GetParamSubVals(param, chain, state[chain]);

    /* get old value of model */
    newM = oldM = (int)*GetParamVals(param, chain, state[chain]);
    
    /* get a new model ID */
    do
        {
        newM = (int)(RandomNumber(seed) * 10);
        } while (newM == oldM);

    /* set proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* set prior ratio */
    *lnPriorRatio = subValue[newM] - subValue[oldM];
    
    /* copy new amino acid model ID back */
    *GetParamVals(param, chain, state[chain]) = (MrBFlt)newM;
    
    /* set amino acid frequencies */
    bs = GetParamSubVals (modelSettings[param->relParts[0]].stateFreq, chain, state[chain]);
    if (newM == AAMODEL_POISSON)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = 1.0 / 20.0;
        }
    else if (newM == AAMODEL_JONES)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = jonesPi[i];
        }
    else if (newM == AAMODEL_DAY)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = dayhoffPi[i];
        }
    else if (newM == AAMODEL_MTREV)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = mtrev24Pi[i];
        }
    else if (newM == AAMODEL_MTMAM)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = mtmamPi[i];
        }
    else if (newM == AAMODEL_WAG)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = wagPi[i];
        }
    else if (newM == AAMODEL_RTREV)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = rtrevPi[i];
        }
    else if (newM == AAMODEL_CPREV)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = cprevPi[i];
        }
    else if (newM == AAMODEL_VT)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = vtPi[i];
        }
    else if (newM == AAMODEL_BLOSUM)
        {
        for (i=0; i<mp->nStates; i++)
            bs[i] = blosPi[i];
        }

    /* Set update flags for all partitions that share this amino acid model. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


int Move_AddDeleteCPPEvent (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* add or delete one Poisson process event */
    
    int         i, k, addEvent, *nEvents, numEvents;
    MrBFlt      sigma, m, lognormalLnProb, **position, **rateMultiplier, length, pos, rate;
    TreeNode    *p, *q;
    ModelInfo   *model;
    Tree        *t;

    /* get the model settings */
    model = &modelSettings[param->relParts[0]];

    /* get cpp rate */
    rate = *GetParamVals (model->cppRate, chain, state[chain]);
    
    /* get sigma of lognormal of rate multipliers */
    sigma = *GetParamVals (model->cppMultDev, chain, state[chain]);

    /* get the cpp event data */
    nEvents = param->nEvents[2*chain+state[chain]];
    position = param->position[2*chain+state[chain]];
    rateMultiplier = param->rateMult[2*chain+state[chain]];

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* pick a branch */
    do
        {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes -2))];
        } while (p->anc == NULL || (p->anc->anc == NULL));

    /* get number of events for convenience */
    numEvents = nEvents[p->index];

    /* add or delete ? */
    addEvent = NO;
    if (numEvents == 0)
        addEvent = YES;
    else if (RandomNumber(seed) < 0.5)
        addEvent = YES;
    
    if (addEvent == NO)
        {
        /* delete event */

        /* choose random event */
        k = (int) (RandomNumber(seed) * numEvents);
        
        /* save multiplier to be deleted */
        m = rateMultiplier[p->index][k];

        /* rearrange and reduce */
        for (i=k; i<numEvents-1; i++)
            {
            position[p->index][i] = position[p->index][i+1];
            rateMultiplier[p->index][i] = rateMultiplier[p->index][i+1];
            }
        if (numEvents-1 > 0)
            {
            position[p->index] = (MrBFlt *) SafeRealloc ((void *) position[p->index], (numEvents-1)*sizeof(MrBFlt));
            rateMultiplier[p->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[p->index], (numEvents-1)*sizeof(MrBFlt));
            assert (position[p->index] != NULL && rateMultiplier[p->index] != NULL);
            }
        else
            {
            free (position[p->index]);
            free (rateMultiplier[p->index]);
            position[p->index] = rateMultiplier[p->index] = NULL;
            }
        /* update number of events */
        nEvents[p->index]--;
        }
    else /* if (addEvent == YES) */
        {
        /* add event */

        /* generate new multiplier */
        m = LogNormalRandomVariable (0.0, sigma, seed);

        /* generate new position */
        pos = RandomNumber(seed);

        /* find place in current array */
        for (k=0; k<numEvents; k++)
            {
            if (position[p->index][k] > pos)
                break;
            }

        /* rearrange and insert */
        position[p->index] = (MrBFlt *) SafeRealloc ((void *)position[p->index], (numEvents+1)*sizeof(MrBFlt));
        rateMultiplier[p->index] = (MrBFlt *) SafeRealloc ((void *)rateMultiplier[p->index], (numEvents+1)*sizeof(MrBFlt));
        assert (position[p->index] != NULL && rateMultiplier[p->index] != NULL);
        for (i=numEvents; i>k; i--)
            {
            position[p->index][i] = position[p->index][i-1];
            rateMultiplier[p->index][i] = rateMultiplier[p->index][i-1];
            }
        position[p->index][k] = pos;
        rateMultiplier[p->index][k] = m;

        /* update number of events */
        nEvents[p->index]++;
        }
    
    /* the CPP process is relative to expected substitutions */
    length = p->length;
    
    lognormalLnProb = LnProbLogNormal(0.0, sigma, m);
    if (addEvent == YES)
        (*lnPriorRatio) = lognormalLnProb + log (rate); 
    else
        (*lnPriorRatio) = -(lognormalLnProb + log(rate));

    if (addEvent == YES)
        /* note that nEvents[p->index] now contains k+1 after addition */
        (*lnProposalRatio) = log (length / ((double) nEvents[p->index])) - lognormalLnProb;
    else
        /* note that nEvents[p->index] contains k after deletion */
        (*lnProposalRatio) = log ((double)(nEvents[p->index]+1) / length) + lognormalLnProb;

    /* take care of asymmetric add and delete probabilities around 0 and 1 events */
    if (addEvent == YES && nEvents[p->index] == 1)
        (*lnProposalRatio) += log (0.5);
    else if (addEvent == NO && nEvents[p->index] == 0)
        (*lnProposalRatio) += log (2.0);

    /* update evolLengths in subtree above new event */
    if (UpdateCppEvolLengths (param, p, chain)==ERROR)
        {
        abortMove=YES;
        return (NO_ERROR);
        }
    
    /* set update of cond likes down to root */
    /* crown tree update flags set in UpdateCppEvolLengths */
    q = p->anc;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }

    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


int Move_Adgamma (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change correlation parameter (-1, 1) of adgamma model */

    int         i, isValidP;
    MrBFlt      oldP, newP, window, minP, maxP, ran, *markovTiValues;
    ModelParams *mp;

    /* get size of window, centered on current rho */
    window = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* get minimum and maximum values for rho */
    minP = mp->corrUni[0];
    maxP = mp->corrUni[1];

    /* get address of markovTi */
    markovTiValues = GetParamSubVals (param, chain, state[chain]);

    /* get old value of rho */
    newP = oldP = *GetParamVals(param, chain, state[chain]);

    /* change value for rho */
    ran = RandomNumber(seed);
     if (maxP-minP < window)
        {
        window = maxP-minP;
        }
    newP = oldP + window * (ran - 0.5);
    
    /* check that new value is valid */
    isValidP = NO;
    do
        {
        if (newP < minP)
            newP = 2* minP - newP;
        else if (newP > maxP)
            newP = 2 * maxP - newP;
        else
            isValidP = YES;
        } while (isValidP == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* get prior ratio */
    *lnPriorRatio = 0.0;
    
    /* copy new rho value back */
    *GetParamVals(param, chain, state[chain]) = newP;

    /* fill in new Markov trans probs */
    AutodGamma (markovTiValues, newP, mp->numGammaCats);
        
    /* Set update flags for all partitions that share this rho. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Update flags for divisions already set */

    return (NO_ERROR);
}


int Move_Beta (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change symmetric Dirichlet variance using multiplier */

    int         i, j, k, isValidB, isPriorExp, nStates;
    MrBFlt      oldB, newB, minB, maxB, priorExp=0.0, *bs, ran, factor, tuning,
                x, y;
    ModelParams *mp;

    /* get tuning parameter */
    tuning = mvp[0];    /* multiplier tuning parameter lambda */

    /* get model paramaters */
    mp = &modelParams[param->relParts[0]];

    /* get prior, minimum and maximum values for rate     */
    if (!strcmp(mp->symPiPr,"Uniform"))
        {
        isPriorExp = NO;
        minB = mp->symBetaUni[0];
        maxB = mp->symBetaUni[1];
        }
    else
        {
        isPriorExp = YES;
        priorExp = mp->symBetaExp;
        minB = SYMPI_MIN;
        maxB = SYMPI_MAX;
        }

    /* get old value of symDir */
    oldB = *GetParamVals(param, chain, state[chain]);

    /* change value */
    ran = RandomNumber(seed);
    factor = exp(tuning * (ran - 0.5));
    newB = oldB * factor;

    /* check validity */
    isValidB = NO;
    do
        {
        if (newB < minB)
            newB = minB * minB / newB;
        else if (newB > maxB)
            newB = maxB * maxB / newB;
        else
            isValidB = YES;
        } while (isValidB == NO);

    /* set new value of symDir */
    *GetParamVals(param, chain, state[chain]) = newB;

    /* get proposal ratio */
    *lnProposalRatio = log (newB / oldB);

    /* get prior ratio */
    if (isPriorExp == YES)
        {
        *lnPriorRatio = priorExp * (oldB - newB);
        }
    else
        *lnPriorRatio = 0.0;

    /* fill in the new betacat frequencies */
    bs = GetParamStdStateFreqs(param, chain, state[chain]);
    k = mp->numBetaCats;
    BetaBreaks (newB, newB, bs, k);
    k *= 2;
    for (i=k-2; i>0; i-=2)
        {
        bs[i] = bs[i/2];
        }
    for (i=1; i<k; i+=2)
        {
        bs[i] = 1.0 - bs[i-1];
        }
        
    /* if there are multistate characters, update prior probability of current pis */
    bs += 2 * mp->numBetaCats;
    for (i=0; i<param->nSympi; i++)
        {
        /* get number of states */
        nStates = param->sympinStates[i];

        /* get prior ratio update */
        x = LnGamma(newB*nStates) - nStates*LnGamma(newB);
        y = LnGamma(oldB*nStates) - nStates*LnGamma(oldB);
        for (j=0; j<nStates; j++)
            {
            x += (newB-1.0)*log(bs[j]);
            y += (oldB-1.0)*log(bs[j]);
            }
        (*lnPriorRatio) += x - y;

        /* update pointer to state freqs */
        bs += nStates;
        }
        
    /* Set update flags for all tree nodes. Note that the conditional
       likelihood update flags have been set for the relevant partitions
       before we even call the move function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* may need to hit update flag for cijks if we have multistate characters */
    for (i=0; i<param->nRelParts; i++)
        {
        if (modelSettings[param->relParts[i]].nCijkParts > 0)
            modelSettings[param->relParts[i]].upDateCijk = YES;
        }

    return (NO_ERROR);
}


int Move_BrLen (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change one branch length */

    MrBFlt      tuning, maxV, minV, m, newM, brlensPrExp=0.0;
    TreeNode    *p;
    ModelParams *mp;
    Tree        *t;
    int isVPriorExp;
    
    tuning = mvp[0]; /* Larget & Simon's tuning parameter lambda */

    mp = &modelParams[param->relParts[0]];

    /* max and min brlen */
    if (param->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0];
        maxV = mp->brlensUni[1];
        isVPriorExp = NO;
        }
    else if (param->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensPrExp = mp->brlensExp;
        isVPriorExp = YES;  
        }

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);

    /* pick a branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * t->nNodes)];
        }
    while (p->anc == NULL || (t->isRooted == YES && p->anc->anc == NULL));

    /* determine new length */
    m = p->length;
    newM = m * exp(tuning * (RandomNumber(seed) - 0.5));

    /* reflect new length if necessary */
    while (newM < minV || newM > maxV)
        {
        if (newM < minV)
            newM = minV * minV / newM;
        else if (newM > maxV)
            newM = maxV * maxV / newM;
        }
    p->length = newM;

    /* calculate proposal ratio */
    /* must be based on new length after reflection */
    (*lnProposalRatio) = log(newM / m);

    /* set flags for update of transition probabilities at p */
    p->upDateTi = YES;

    /* set the update flag for cond likes if p is connected to root in unrooted */
    /* tree, if this is not done, cond likes are not updated in this case       */
    if (t->isRooted == NO && p->anc->anc == NULL)  
        p->upDateCl = YES;

    /* set flags for update of cond likes from p->anc and down to root */
    while (p->anc->anc != NULL)
        {
        p = p->anc;
        p->upDateCl = YES; 
        }

    /* update prior if exponential prior on branch lengths */
    if (param->paramId == BRLENS_EXP)
        (*lnPriorRatio) = brlensPrExp * (m - newM);
    /* Dirichlet or twoExp prior */
    else if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);

    return (NO_ERROR);
}


int Move_ClockRate_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change clock rate using multiplier */
    
    int             i, j, k, *nEvents;
    MrBFlt          oldRate, newRate, factor, lambda, nu, igrvar, *brlens, *igrRate, *tk02Rate,
                    N, newTheta, oldTheta, growth, newLnPrior, oldLnPrior;
    Tree            *t, *oldT;
    TreeNode        *p, *q;
    Param           *treeParam, *subParm;
    ModelParams     *mp;
    ModelInfo       *m;
    
    /* get old value of clock rate */
    oldRate = *GetParamVals(param, chain, state[chain]);

    /* Rely on general algorithm to change the value */
    Move_PosRealMultiplier(param, chain, seed, lnPriorRatio, lnProposalRatio, mvp);
    if (abortMove == YES)
        return NO_ERROR;

    /* get new value of clock rate */
    newRate = *GetParamVals(param, chain, state[chain]);

    /* calculate factor */
    factor = newRate / oldRate;

    /* clock rate applies to all clock trees */
    for (i = 0; i < numTrees; i++)
        {
        t = GetTreeFromIndex(i, chain, state[chain]);
        if (t->isClock == NO)
            continue;
        if (!strcmp(modelParams[t->relParts[0]].clockPr, "Fixed"))
            continue;

        oldT = GetTreeFromIndex(i, chain, 1^state[chain]);
        treeParam = modelSettings[t->relParts[0]].brlens;
        
        /* adjust the node depths and lengths */
        for (j = 0; j < t->nNodes-1; j++)
            {
            p = t->allDownPass[j];
            q = oldT->allDownPass[j];
            p->nodeDepth *= factor; /* no harm done if nodeDepth==0.0 (undated tips) */
            p->length *= factor;    /* no harm done if length==0.0 (root or fossil ancestors)*/
            if (p->length < 0.0 || p->length > BRLENS_MAX ||
                (q->length > BRLENS_MIN && p->length < BRLENS_MIN) ||
                (q->length < TIME_MIN   && p->length > TIME_MIN))
                {  /* consider ancestral fossil (brl=0) in fossilized bd tree */
                abortMove = YES;
                return (NO_ERROR);
                }
            }

        /* prior ratio for coalecent tree, as theta is changed */
        mp = &modelParams[t->relParts[0]];
        if (!strcmp(mp->clockPr,"Coalescence"))
            {
            m = &modelSettings[t->relParts[0]];
            N = *GetParamVals(m->popSize, chain, state[chain]);
            if (!strcmp(mp->ploidy, "Diploid"))
                N *= 4.0;
            else if (!strcmp(mp->ploidy, "Zlinked"))
                N *= 3.0;
            else
                N *= 2.0;
            oldTheta = N * oldRate;
            newTheta = N * newRate;
            if (!strcmp(mp->growthPr, "Fixed"))
                growth = mp->growthFix;
            else
                growth = *GetParamVals(m->growthRate, chain, state[chain]);

            if (LnCoalescencePriorPr (oldT, &oldLnPrior, oldTheta, growth) == ERROR)
                {
                MrBayesPrint ("%s   Problem calculating prior for coalescent process\n", spacer);
                return (ERROR);
                }
            if (LnCoalescencePriorPr (t, &newLnPrior, newTheta, growth) == ERROR)
                {
                MrBayesPrint ("%s   Problem calculating prior for coalescent process\n", spacer);
                return (ERROR);
                }
            (*lnPriorRatio) += newLnPrior - oldLnPrior;
            }
            
        /* adjust proposal and prior ratio for relaxed clock models */
        for (k = 0; k < treeParam->nSubParams; k++)
            {
            subParm = treeParam->subParams[k];
            if (subParm->paramType == P_CPPEVENTS)
                {
                nEvents = subParm->nEvents[2*chain+state[chain]];
                lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);
                /* proposal ratio */
                for (j=0; j<t->nNodes-2; j++)
                    {
                    p = t->allDownPass[j];
                    q = oldT->allDownPass[j];
                    (*lnProposalRatio) += nEvents[p->index ] * log (p->length  / q->length);
                    }
                /* prior ratio */
                (*lnPriorRatio) += lambda * (TreeLen(oldT) - TreeLen(t));
                /* update effective evolutionary lengths */
                if (UpdateCppEvolLengths (subParm, t->root->left, chain) == ERROR)
                    {
                    abortMove = YES;
                    return (NO_ERROR);
                    }
                }
            else if ( subParm->paramType == P_TK02BRANCHRATES ||
                     (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
                {
                if (subParm->paramType == P_TK02BRANCHRATES)
                    nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
                else
                    nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
                tk02Rate = GetParamVals (subParm, chain, state[chain]);
                brlens = GetParamSubVals (subParm, chain, state[chain]);

                /* no proposal ratio effect */

                /* prior ratio and update of brlens */
                for (j = 0; j < t->nNodes-2; j++)
                    {
                    p = t->allDownPass[j];
                    q = oldT->allDownPass[j];
                    if (p->length > 0.0)  // not ancestral fossil
                        {
                        (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[q->anc->index], nu*q->length, tk02Rate[q->index]);
                        (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->anc->index], nu*p->length, tk02Rate[p->index]);
                        
                        brlens[p->index] = p->length * (tk02Rate[p->anc->index]+tk02Rate[p->index])/2.0;
                        }
                    }
                }
            else if ( subParm->paramType == P_IGRBRANCHRATES ||
                     (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
                {
                if (subParm->paramType == P_IGRBRANCHRATES)
                    igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
                else
                    igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
                igrRate = GetParamVals (subParm, chain, state[chain]);
                brlens = GetParamSubVals (subParm, chain, state[chain]);
            
                /* prior ratio and update of brlens */
                for (j = 0; j < t->nNodes-2; j++)
                    {
                    p = t->allDownPass[j];
                    q = oldT->allDownPass[j];
                    if (p->length > 0.0)  // not ancestral fossil
                        {
                        (*lnPriorRatio) -= LnProbGamma (q->length/igrvar, q->length/igrvar, igrRate[q->index]);
                        (*lnPriorRatio) += LnProbGamma (p->length/igrvar, p->length/igrvar, igrRate[p->index]);

                        brlens[p->index] = igrRate[p->index] * p->length;
                        }
                    }
                }
            }
        }
 
    return (NO_ERROR);
}


int Move_CPPEventPosition (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move the position of one CPP event */

    int         i, j, k, *nEvents;
    MrBFlt      pos, temp, **position, **rateMultiplier;
    TreeNode    *p=NULL, *q;
    Tree        *t;

    /* get the cpp event data */
    nEvents = param->nEvents[2*chain+state[chain]];
    position = param->position[2*chain+state[chain]];
    rateMultiplier = param->rateMult[2*chain+state[chain]];

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* pick a branch and an event */
    for (i=j=0; i<t->nNodes - 2; i++)
        {
        p = t->allDownPass[i];
        j += nEvents[p->index];
        }       
    if (j == 0)
        {
        abortMove = YES;
        return (NO_ERROR);
        }
    k = (int) (RandomNumber(seed) * j);
    for (i=j=0; i<t->nNodes - 2; i++)
        {
        p = t->allDownPass[i];
        j += nEvents[p->index];
        if (j > k)
            break;
        }       
    if (position[p->index] == NULL)
        getchar();

    /* find local index */
    k = k - (j - nEvents[p->index]);

    /* find new position */
    pos = RandomNumber(seed);
    if (pos < POS_MIN || 1.0 - pos < POS_MIN)
        {
        abortMove = YES;
        return (NO_ERROR);
        }
    position[p->index][k] = pos;

    /* sort events; bubble sort for now */
    for (i=0; i<nEvents[p->index]; i++)
        {
        for (j=i+1; j<nEvents[p->index]; j++)
            {
            if (position[p->index][j] < position[p->index][i])
                {
                temp = position[p->index][i];
                position[p->index][i] = position[p->index][j];
                position[p->index][j] = temp;
                temp = rateMultiplier[p->index][i];
                rateMultiplier[p->index][i] = rateMultiplier[p->index][j];
                rateMultiplier[p->index][j] = temp;
                }
            }
        }

    /* calculate prior and proposal ratio */
    (*lnPriorRatio) = (*lnProposalRatio) = 0.0;

    /* update branch evolution lengths */
    if (UpdateCppEvolLengths (param, p, chain) == ERROR)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* set update of cond likes down to root */
    /* update of crowntree set in UpdateCppEvolLengths */
    q = p->anc;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }

    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


int Move_CPPRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move the CPP rate (lambda) using multiplier */

    int         i, j, *nEvents, sumEvents;
    MrBFlt      oldLambda, newLambda, treeLength, tuning;
    Model       *mp;
    TreeNode    *p;
    Tree        *t;

    /* get tuning parameter */
    tuning = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get the CPP rate */
    oldLambda = *GetParamVals (param, chain, state[chain]);

    /* set new value */
    newLambda = oldLambda * exp ((0.5 - RandomNumber(seed))*tuning);
    
    /* reflect if necessary */
    while (newLambda < CPPLAMBDA_MIN || newLambda > CPPLAMBDA_MAX)
        {
        if (newLambda < CPPLAMBDA_MIN)
            newLambda = CPPLAMBDA_MIN * CPPLAMBDA_MIN / newLambda;
        if (newLambda > CPPLAMBDA_MAX)
            newLambda = CPPLAMBDA_MAX * CPPLAMBDA_MAX / newLambda;
        }
    
    /* store new value */
    (*GetParamVals (param, chain, state[chain])) = newLambda;

    /* calculate prior ratio */
    (*lnPriorRatio) = 0.0;
    for (i=0; i<param->nSubParams; i++)
        {
        nEvents = param->subParams[i]->nEvents[2*chain+state[chain]];
        sumEvents = 0;
        t = GetTree (param->subParams[i], chain, state[chain]); 
        treeLength = 0.0;
        for (j=0; j<t->nNodes-2; j++)
            {
            p = t->allDownPass[j];
            sumEvents += nEvents[p->index];
            treeLength += p->length;
            }
        (*lnPriorRatio) += (oldLambda - newLambda) * treeLength;
        (*lnPriorRatio) += sumEvents * log (newLambda / oldLambda);
        }

    /* adjust for prior on cppRate */
    if (!strcmp(mp->cppRatePr,"Exponential"))
        (*lnPriorRatio) +=  mp->cppRateExp * (oldLambda - newLambda);

    /* calculate proposal ratio */
    (*lnProposalRatio) = log (newLambda / oldLambda);

    /* we do not need to update likelihoods */
    for (i=0; i<param->nRelParts; i++)
        {
        modelSettings[param->relParts[i]].upDateCl = NO;
        }

    return (NO_ERROR);
}


int Move_CPPRateMultiplier_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move one CPP rate multiplier using multiplier */

    int         i, j, k, *nEvents;
    MrBFlt      newRateMultiplier, oldRateMultiplier, tuning, minM, maxM, sigma, **rateMultiplier;
    TreeNode    *p = NULL;
    ModelInfo   *m;
    Tree        *t;
    TreeNode    *q;

    /* get the tuning parameter */
    tuning = mvp[0];
    
    /* get the model settings */
    m = &modelSettings[param->relParts[0]];

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get CPP event data */
    nEvents = param->nEvents[2*chain+state[chain]];
    rateMultiplier = param->rateMult[2*chain+state[chain]];

    /* get minimum and maximum of CPP rate multiplier */
    minM = param->min;
    maxM = param->max;
    
    /* pick a branch and a rateMultiplier */
    for (i=j=0; i<t->nNodes - 2; i++)
        {
        p = t->allDownPass[i];
        j += nEvents[p->index];
        }       
    if (j == 0)
        {
        abortMove = YES;
        return (NO_ERROR);
        }
    k = (int) (RandomNumber(seed) * j);
    for (i=j=0; i<t->nNodes - 2; i++)
        {
        p = t->allDownPass[i];
        j += nEvents[p->index];
        if (j > k)
            break;
        }

    /* find local index */
    k = nEvents[p->index] - (j - k);

    /* find new rateMultiplier */
    oldRateMultiplier = rateMultiplier[p->index][k];
    newRateMultiplier = oldRateMultiplier * (exp (0.5 - RandomNumber(seed) * tuning));

    /* reflect if necessary */
    while (newRateMultiplier < minM || newRateMultiplier > maxM)
        {
        if (newRateMultiplier < minM)
            newRateMultiplier = minM * minM / newRateMultiplier;
        if (newRateMultiplier > maxM)
            newRateMultiplier = maxM * maxM / newRateMultiplier;
        }

    rateMultiplier[p->index][k] = newRateMultiplier;
    
    /* calculate prior ratio */
    sigma = *GetParamVals (m->cppMultDev, chain, state[chain]);
    (*lnPriorRatio) = LnRatioLogNormal (0.0, sigma, newRateMultiplier, oldRateMultiplier);

    /* calculate proposal ratio */
    (*lnProposalRatio) = log (newRateMultiplier / oldRateMultiplier);

    /* update branch evolution lengths */
    if (UpdateCppEvolLengths (param, p, chain)==ERROR)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* set update of cond likes down to root */
    /* update of crowntree set in UpdateCppEvolLengths */
    q = p->anc;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }

    return (NO_ERROR);
}


int Move_CPPRateMultiplierRnd (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move one CPP rate multiplier by redrawing from prior */

    int         i, j, k, *nEvents;
    MrBFlt      sigma, newRateMultiplier, oldRateMultiplier, **rateMultiplier;
    TreeNode    *p=NULL, *q;
    ModelInfo   *m;
    Tree        *t;

    /* get the model settings */
    m = &modelSettings[param->relParts[0]];

    /* get the CPP event data */
    nEvents = param->nEvents[2*chain+state[chain]];
    rateMultiplier = param->rateMult[2*chain+state[chain]];

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* pick a branch and a rateMultiplier */
    for (i=j=0; i<t->nNodes - 2; i++)
        {
        p = t->allDownPass[i];
        j += nEvents[p->index];
        }       
    if (j == 0)
        {
        abortMove = YES;
        return (NO_ERROR);
        }
    k = (int) (RandomNumber(seed) * j);
    for (i=j=0; i<t->nNodes - 2; i++)
        {
        p = t->allDownPass[i];
        j += nEvents[p->index];
        if (j > k)
            break;
        }

    /* find local index */
    k = nEvents[p->index] - (j - k);

    /* record old rate multiplier */
    oldRateMultiplier = rateMultiplier[p->index][k];

    /* find stdev of lognormal */
    sigma = *GetParamVals (m->cppMultDev, chain, state[chain]);

    /* set new value */
    do {
        newRateMultiplier = LogNormalRandomVariable (0.0, sigma, seed);
        } while (newRateMultiplier < param->min || newRateMultiplier > param->max);
    rateMultiplier[p->index][k] = newRateMultiplier;

    /* calculate prior ratio */
    (*lnPriorRatio) = LnRatioLogNormal(0.0, sigma, newRateMultiplier, oldRateMultiplier);
        
    /* calculate proposal ratio */
    (*lnProposalRatio) += LnRatioLogNormal (0.0, sigma, oldRateMultiplier, newRateMultiplier);

    /* update branch evolution lengths */
    if (UpdateCppEvolLengths (param, p, chain) == ERROR)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* set update of cond likes down to root */
    /* update of crowntree set in UpdateCppEvolLengths */
    q = p->anc;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }

    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


int Move_AddBranch (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Move an ancestral fossil (brl = 0) to fossil tip (brl > 0)

       __|__              __|__
      |     |            |     |
            |       -->       q|___
            |       -->        |   |
           q|___p              |   |p
           r|                 r|
     
       1. Pich a fossil among those with brl = 0 (prob = 1/k)
       2. Propose brl from a uniform(0, ?) distribution
     */
    
    int    i, j, k, mFossil, kFossil;
    MrBFlt minDepth, maxDepth, newLength, clockRate, x, oldQLength, oldRLength,
           *brlens=NULL, nu=0.0, *tk02Rate=NULL, igrvar=0.0, *igrRate=NULL;
    TreeNode    *p=NULL, *q=NULL, *r;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m;
    Param       *subParm;
    Calibration *calibrationPtr = NULL;
    
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* get number of ancestral and tip fossils */
    kFossil = mFossil = 0;
    for (i = 0; i < t->nNodes -1; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;  // reset marked node
        if (p->left == NULL && p->right == NULL && p->nodeDepth > 0.0)
            {
            if (p->length > 0.0)
                {
                mFossil++;        // count tip fossil
                }
            else
                {
                p->marked = YES;  // mark  anc fossil
                kFossil++;        // count anc fossil
                }
            }
        }
    if (kFossil == 0)  // no ancestral fossil, nothing to do
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get clock rate */
    clockRate = *(GetParamVals(m->clockRate, chain, state[chain]));
    
    /* pick an ancestral fossil randomly */
    j = (int) (RandomNumber(seed) * kFossil);
    for (i = k = 0; i < t->nNodes -1; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES)
            k++;
        if (k > j)
            break;
        }
    /* now p is pointing to the ancestral fossil
       whose brl needs to be changed to >0. let's do it! */
    q = p->anc;
    if (q->left == p)
        r = q->right;
    else
        r = q->left;

    /* determine lower and upper bound of forward move, abort if impossible */
    minDepth = p->nodeDepth + BRLENS_MIN;
    if (q->anc->anc == NULL)
        maxDepth = TREEHEIGHT_MAX;
    else
        maxDepth = q->anc->nodeDepth - BRLENS_MIN;
    
    if (q->isDated == YES)
        calibrationPtr = q->calibration;
    else if (q->anc->anc == NULL)  // q is root but not dated
        calibrationPtr = &mp->treeAgePr;
    
    if (calibrationPtr != NULL)
        {
        if (calibrationPtr->prior == fixed || calibrationPtr->min * clockRate > minDepth)
            {
            abortMove = YES;
            return (NO_ERROR);
            }
        if (calibrationPtr->max * clockRate < maxDepth)
            maxDepth = calibrationPtr->max * clockRate;
        }
    if (minDepth >= maxDepth)
        {
        abortMove = YES;
        return (NO_ERROR);
        }
    
    /* record old lengths and depths */
    oldQLength = q->length;
    oldRLength = r->length;
    // oldDepth = q->nodeDepth;
    
    /* propose the branch length leading to the fossil */
    newLength = (RandomNumber(seed)) * (maxDepth - minDepth);
    
    /* adjust brls and depths, set flags for update of trans probs */
    p->length   = newLength;
    p->upDateTi = YES;
    q->nodeDepth += newLength;
    if (q->anc->anc != NULL)
        {
        q->length  -= newLength;
        q->upDateTi = YES;
        }
    r->length += newLength;
    r->upDateTi = YES;
    // newDepth = q->nodeDepth;

    /* adjust age of q if dated */
    if (calibrationPtr != NULL)
        {
        q->age = q->nodeDepth / clockRate;
        }
    
    /* set flags for update of cond likes from p/r to root */
    r->upDateCl = YES;
    q = p;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }
    q = p->anc;

    /* calculate prior ratio */
    if (LogClockTreePriorRatio(param, chain, &x) == ERROR)
        return (ERROR);
    (*lnPriorRatio) += x;
    
    /* calculate proposal ratio; because we abort in the extreme conditions (k == 0 and m == 0), no correction for this needed */
    /* note that we assume here that the relative proposal probabilities are the same for the add move and the delete move */
    (*lnProposalRatio) = log(kFossil) - log(mFossil +1);
    
    /* add the Jacobian term */
    (*lnProposalRatio) += log((maxDepth - minDepth) / clockRate);
            
    /* adjust proposal and prior ratio for relaxed clock models */
    for (i=0; i<param->nSubParams; i++)
        {
        subParm = param->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            /* CPP is not compatible with FBD ancestral fossils until we have a good way to implement it !! */
            MrBayesPrint ("%s   CPP clock is not compatible with Fossilization prior currently\n", spacer);
            return (ERROR);
            }
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            
            /* prior ratio */
            tk02Rate[p->index] = tk02Rate[q->index];
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[q->index], nu*oldRLength, tk02Rate[r->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[q->index], nu* p->length, tk02Rate[p->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[q->index], nu* r->length, tk02Rate[r->index]);
            if (q->anc->anc != NULL)
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[q->anc->index], nu*oldQLength, tk02Rate[q->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[q->anc->index], nu* q->length, tk02Rate[q->index]);
                }
            
            /* update effective evolutionary lengths */
            brlens[p->index] = p->length * (tk02Rate[p->index]+tk02Rate[q->index])/2.0;
            brlens[r->index] = r->length * (tk02Rate[r->index]+tk02Rate[q->index])/2.0;
            if (q->anc->anc != NULL)
                {
                brlens[q->index] = q->length * (tk02Rate[q->index]+tk02Rate[q->anc->index])/2.0;
                }
            }
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            
            /* prior ratio */
            igrRate[p->index] = igrRate[q->index];
            (*lnPriorRatio) -= LnProbGamma (oldRLength/igrvar, oldRLength/igrvar, igrRate[r->index]);
            (*lnPriorRatio) += LnProbGamma (p->length /igrvar, p->length /igrvar, igrRate[p->index]);
            (*lnPriorRatio) += LnProbGamma (r->length /igrvar, r->length /igrvar, igrRate[r->index]);
            if (q->anc->anc != NULL)
                {
                (*lnPriorRatio) -= LnProbGamma (oldQLength/igrvar, oldQLength/igrvar, igrRate[q->index]);
                (*lnPriorRatio) += LnProbGamma (q->length /igrvar, q->length /igrvar, igrRate[q->index]);
                }
            
            /* update effective evolutionary lengths */
            brlens[p->index] = igrRate[p->index] * p->length;
            brlens[r->index] = igrRate[r->index] * r->length;
            if (q->anc->anc != NULL)
                {
                brlens[q->index] = igrRate[q->index] * q->length;
                }
            }
        }
    
    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


int Move_DelBranch (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Move a fossil tip (brl > 0) to be ancestral (brl =0)

       __|__              __|__
      |     |            |     |
           q|___     -->       |
            |   |    -->       |
            |   |p            q|___p
           r|                 r|
         
       1. Pich a fossil among those with brl > 0 (prob = 1/m)
       2. Set brl = 0
     */
    
    int    i, j, k, mFossil, kFossil;
    MrBFlt minDepth, maxDepth, clockRate, x, oldPLength, oldQLength, oldRLength,
           *brlens=NULL, nu=0.0, *tk02Rate=NULL, igrvar=0.0, *igrRate=NULL;
    TreeNode    *p=NULL, *q=NULL, *r;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m;
    Param       *subParm;
    Calibration *calibrationPtr = NULL;

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* get number of ancestral and tip fossils */
    kFossil = mFossil = 0;
    for (i = 0; i < t->nNodes -1; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;  // reset marked node
        if (p->left == NULL && p->right == NULL && p->nodeDepth > 0.0)
            {
            if (p->length > 0.0)
                {
                p->marked = YES;  // mark  tip fossil
                mFossil++;        // count tip fossil
                }
            else
                {
                kFossil++;        // count anc fossil
                }
            }
        }
    if (mFossil == 0)  // no tip fossil, nothing to do
        {
        abortMove = YES;
        return (NO_ERROR);
        }
    
    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get clock rate */
    clockRate = *(GetParamVals(m->clockRate, chain, state[chain]));

    /* pick a tip fossil randomly */
    j = (int) (RandomNumber(seed) * mFossil);
    for (i = k = 0; i < t->nNodes -1; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES)
            k++;
        if (k > j)
            break;
        }
    /* now p is pointing to the fossil tip
       whose brl needs to be changed to 0. let's do it */
    q = p->anc;
    if (q->left == p)
        r = q->right;
    else
        r = q->left;

    /* determine lower and upper bound of backward move, abort if impossible */
    minDepth = p->nodeDepth + BRLENS_MIN;
    if (q->anc->anc == NULL)
        maxDepth = TREEHEIGHT_MAX;
    else
        maxDepth = q->anc->nodeDepth - BRLENS_MIN;
    
    if (q->isDated == YES)
        calibrationPtr = q->calibration;
    else if (q->anc->anc == NULL)  // q is root but not dated
        calibrationPtr = &mp->treeAgePr;
    
    if (calibrationPtr != NULL)
        {
        if (calibrationPtr->prior == fixed || calibrationPtr->min * clockRate > minDepth)
            {
            abortMove = YES;
            return (NO_ERROR);
            }
        if (calibrationPtr->max * clockRate < maxDepth)
            maxDepth = calibrationPtr->max * clockRate;
        }
    if (r->nodeDepth > p->nodeDepth -BRLENS_MIN || minDepth >= maxDepth)
        {  /* the sister node (another fossil) is older than the current fossil */
        abortMove = YES;
        return (NO_ERROR);
        }

    /* record old lengths and depths */
    oldPLength = p->length;
    oldQLength = q->length;
    oldRLength = r->length;
    // oldDepth = q->nodeDepth;

    /* set the brl to 0 for the fossil tip, it becomes an ancestral fossil */
    /* set flags for update of transition probabilities too */
    q->nodeDepth = p->nodeDepth;
    if (q->anc->anc != NULL)
        {
        q->length += p->length;
        q->upDateTi = YES;
        }
    r->length  -= p->length;
    r->upDateTi = YES;
    p->length   = 0.0;
    p->upDateTi = YES;
    // newDepth = q->nodeDepth;
    
    /* adjust age of q if dated */
    if (calibrationPtr != NULL)
        {
        q->age = q->nodeDepth / clockRate;
        }
    
    /* set flags for update of cond likes from p/r to root */
    r->upDateCl = YES;
    q = p;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }
    q = p->anc;

    /* calculate prior ratio */
    if (LogClockTreePriorRatio(param, chain, &x) == ERROR)
        return (ERROR);
    (*lnPriorRatio) += x;
    
    /* calculate proposal ratio; since we abort in the extreme cases (m == 0 and k == 0), no need to correct for these */
    /* note that we assume that the add and del move have the same relative proposal prob */
    (*lnProposalRatio) = log(mFossil) - log(kFossil +1);

    /* add the Jacobian term */
    (*lnProposalRatio) -= log((maxDepth - minDepth) / clockRate);

    /* adjust proposal and prior ratio for relaxed clock models */
    for (i=0; i<param->nSubParams; i++)
        {
        subParm = param->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            /* CPP is not compatible with FBD ancestral fossils until we have a good way to implement it !! */
            MrBayesPrint ("%s   CPP clock is not compatible with Fossilization prior currently\n", spacer);
            return (ERROR);
            }
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
           {
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            
            /* prior ratio */
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[q->index], nu*oldPLength, tk02Rate[p->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[q->index], nu*oldRLength, tk02Rate[r->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[q->index], nu* r->length, tk02Rate[r->index]);
            if (q->anc->anc != NULL)
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[q->anc->index], nu*oldQLength, tk02Rate[q->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[q->anc->index], nu* q->length, tk02Rate[q->index]);
                }
            
            /* update effective evolutionary lengths */
            brlens[p->index] = 0.0;  // tk02Rate[p->index] = tk02Rate[q->index];
            brlens[r->index] = r->length * (tk02Rate[r->index]+tk02Rate[q->index])/2.0;
            if (q->anc->anc != NULL)
                {
                brlens[q->index] = q->length * (tk02Rate[q->index]+tk02Rate[q->anc->index])/2.0;
                }
            }
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            
            (*lnPriorRatio) -= LnProbGamma (oldPLength/igrvar, oldPLength/igrvar, igrRate[p->index]);
            (*lnPriorRatio) -= LnProbGamma (oldRLength/igrvar, oldRLength/igrvar, igrRate[r->index]);
            (*lnPriorRatio) += LnProbGamma (r->length /igrvar, r->length /igrvar, igrRate[r->index]);
            if (q->anc->anc != NULL)
                {
                (*lnPriorRatio) -= LnProbGamma (oldQLength/igrvar, oldQLength/igrvar, igrRate[q->index]);
                (*lnPriorRatio) += LnProbGamma (q->length /igrvar, q->length /igrvar, igrRate[q->index]);
                }
            
            /* update effective evolutionary lengths */
            brlens[p->index] = 0.0;  // igrRate[p->index] = igrRate[q->index];
            brlens[r->index] = igrRate[r->index] * r->length;
            if (q->anc->anc != NULL)
                {
                brlens[q->index] = igrRate[q->index] * q->length;
                }
            }
        }

    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


int Move_Extinction (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change relative extinction rate using sliding window */
    
    int         i, isValidM, valIndex;
    MrBFlt      *valPtr, oldM, newM, window, minM, maxM, *sR, *eR, sF, *fR, oldLnPrior, newLnPrior,
                oldProp[2], newProp[2], x, y, *alphaDir, clockRate;
    char        *sS;
    ModelParams *mp;
    ModelInfo   *m;
    Tree        *t;

    /* get size of window, centered on current value */
    window = mvp[0];

    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get minimum and maximum values */
    minM = 0.0;
    maxM = 0.99999;

    /* get pointer to value to be changed */
    valIndex = (int)(RandomNumber(seed) * param->nValues);
    valPtr = GetParamVals(param, chain, state[chain]) + valIndex;
    
    /* get old value */
    oldM = *valPtr;

    /* change value */
    if (window > maxM-minM)
        window = maxM-minM;
    newM = oldM + window * (RandomNumber(seed) - 0.5);
    
    /* check that new value is valid */
    isValidM = NO;
    do  {
        if (newM < minM)
            newM = 2 * minM - newM;
        else if (newM > maxM)
            newM = 2 * maxM - newM;
        else
            isValidM = YES;
        } while (isValidM == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* calculate prior ratio */
    t  = GetTree(modelSettings[param->relParts[0]].brlens,chain,state[chain]);
    sR = GetParamVals (m->speciationRates, chain, state[chain]);
    eR = GetParamVals (param, chain, state[chain]);
    sF = mp->sampleProb;
    sS = mp->sampleStrat;
    clockRate = *GetParamVals (m->clockRate, chain, state[chain]);
    
    if (!strcmp(mp->clockPr,"Birthdeath"))
        {
        if (LnBirthDeathPriorPr (t, clockRate, &oldLnPrior, *sR, *eR, sS, sF) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
            return (ERROR);
            }
        *valPtr = newM;  // update with new value
        if (LnBirthDeathPriorPr (t, clockRate, &newLnPrior, *sR, *eR, sS, sF) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
            return (ERROR);
            }
        }
    else if (!strcmp(mp->clockPr,"Fossilization"))
        {
        fR = GetParamVals (m->fossilizationRates, chain, state[chain]);
        if (LnFossilizationPriorPr (t, clockRate, &oldLnPrior, sR, eR, sF, fR, sS) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
            return (ERROR);
            }
        *valPtr = newM;  // update with new value
        // for (i=0; i<param->nValues; i++)  *(GetParamVals(param, chain, state[chain]) + i) = newM;
        if (LnFossilizationPriorPr (t, clockRate, &newLnPrior, sR, eR, sF, fR, sS) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
            return (ERROR);
            }
        }
    else {
        MrBayesPrint ("%s   Move_Extinction not applicable\n", spacer);
        return (ERROR);
        }
    
    /* get proportions */
    oldProp[0] = oldM;
    oldProp[1] = 1.0 - oldM;
    newProp[0] = newM;
    newProp[1] = 1.0 - newM;
    
    /* adjust prior ratio according to beta distribution */
    alphaDir = mp->extinctionBeta;
    x = y = 0.0;
    for (i=0; i<2; i++)
        x += (alphaDir[i]-1.0)*log(newProp[i]);
    for (i=0; i<2; i++)
        y += (alphaDir[i]-1.0)*log(oldProp[i]);
    (*lnPriorRatio) = x - y + newLnPrior - oldLnPrior;
    
    return (NO_ERROR);
}


int Move_Fossilization (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change fossilization rate using sliding window */
    
    int         i, isValidM, valIndex;
    MrBFlt      *valPtr, oldM, newM, window, minM, maxM, *sR, *eR, sF, *fR, oldLnPrior, newLnPrior,
                oldProp[2], newProp[2], x, y, *alphaDir, clockRate;
    char        *sS;
    ModelParams *mp;
    ModelInfo   *m;
    Tree        *t;
    
    /* get size of window, centered on current value */
    window = mvp[0];
    
    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get minimum and maximum values */
    minM = 0.00001;
    maxM = 1.0;
    
    /* get pointer to value to be changed */
    valIndex = (int)(RandomNumber(seed) * param->nValues);
    valPtr = GetParamVals(param, chain, state[chain]) + valIndex;

    /* get old value */
    oldM = *valPtr;
    
    /* change value */
    if (window > maxM-minM)
        window = maxM-minM;
    newM = oldM + window * (RandomNumber(seed) - 0.5);
    
    /* check that new value is valid */
    isValidM = NO;
    do  {
        if (newM < minM)
            newM = 2 * minM - newM;
        else if (newM > maxM)
            newM = 2 * maxM - newM;
        else
            isValidM = YES;
        } while (isValidM == NO);
    
    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* calculate prior ratio */
    t  = GetTree(modelSettings[param->relParts[0]].brlens,chain,state[chain]);
    sR = GetParamVals (m->speciationRates, chain, state[chain]);
    eR = GetParamVals (m->extinctionRates, chain, state[chain]);
    fR = GetParamVals (param, chain, state[chain]);
    sF = mp->sampleProb;
    sS = mp->sampleStrat;
    clockRate = *GetParamVals(m->clockRate, chain, state[chain]);

    if (LnFossilizationPriorPr (t, clockRate, &oldLnPrior, sR, eR, sF, fR, sS) == ERROR)
        {
        MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
        return (ERROR);
        }
    *valPtr = newM;  // update with new value
    if (LnFossilizationPriorPr (t, clockRate, &newLnPrior, sR, eR, sF, fR, sS) == ERROR)
        {
        MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
        return (ERROR);
        }
    
    /* get proportions */
    oldProp[0] = oldM;
    oldProp[1] = 1.0 - oldM;
    newProp[0] = newM;
    newProp[1] = 1.0 - newM;

    /* adjust prior ratio according to beta distribution */
    alphaDir = mp->fossilizationBeta;
    x = y = 0.0;
    for (i=0; i<2; i++)
        x += (alphaDir[i]-1.0)*log(newProp[i]);
    for (i=0; i<2; i++)
        y += (alphaDir[i]-1.0)*log(oldProp[i]);
    (*lnPriorRatio) = x - y + newLnPrior - oldLnPrior;            
    
    return (NO_ERROR);
}


int Move_ExtFossilSPRClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* This move is identical to the Move_ExtSPRClock move except that it only moves fossil subtrees. */
    
    int         i, j, topologyHasChanged=NO, isStartLocked=NO, isStopLocked=NO, nRootNodes, directionUp,
                n1=0, n2=0, n3=0, n4=0, n5=0, *nEvents, numMovableNodesOld, numMovableNodesNew;
    MrBFlt      x, y, oldBrlen=0.0, newBrlen=0.0, extensionProb, igrvar, *igrRate=NULL,
    v1=0.0, v2=0.0, v3=0.0, v4=0.0, v5=0.0, v3new=0.0, lambda, *tk02Rate=NULL,
    **position=NULL, **rateMultiplier=NULL, *brlens, nu, minV, clockRate;
    TreeNode    *p, *a, *b, *u, *v, *oldA;
    Tree        *t;
    ModelInfo   *m;
    Param       *subParm;
    
    extensionProb = mvp[0]; /* extension probability */
    
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* get model params and model info */
    m = &modelSettings[param->relParts[0]];
    
    /* get clock rate */
    clockRate = *GetParamVals (m->clockRate, chain, state[chain]);
    
    /* get min and max branch lengths in relative time and substitution units */
    minV = BRLENS_MIN;
    
#   if defined (DEBUG_ExtSPRClock)
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
#   endif
    
    /* mark all nodes that only have fossil children with YES and count number movable nodes in current tree */
    numMovableNodesOld=0;
    for (i=0; i<t->nNodes-2; ++i)
    {
        p = t->allDownPass[i];
        if (p->left == NULL)
        {
            if (p->calibration == NULL)
                p->x = NO;
            else
            {
                p->x = YES;
            }
        }
        else
        {
            if (p->left->x == YES && p->right->x == YES)
            {
                p->x = YES;
            }
            else
                p->x = NO;
        }
        a = p->anc->left;
        b = p->anc->right;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL
            || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN) || p->x == NO)
            numMovableNodesOld++;
    }
    
    if (numMovableNodesOld==0)
        return (NO_ERROR);
    
    /* pick a branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes - 2))];
        a = p->anc->left;
        b = p->anc->right;
    }
    while (p->anc->isLocked == YES || p->anc->anc->anc == NULL
           || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN)
           || (p->length < TIME_MIN && p->calibration->prior == fixed)
           || p->x == NO);
    /* skip constraints, siblings of root (and root); and consider ancestral fossils in fbd tree;
       skip all nodes that subtend extant terminals */
    
    /* set up pointers for nodes around the picked branch */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    oldA = a;
    
    /* record branch length for insertion in back move */
    if (v->length > 0.0)  /* side branch, not anc fossil */
    {
        if (v->nodeDepth > a->nodeDepth)
            oldBrlen = b->nodeDepth - v->nodeDepth - 2.0*minV;
        else
            oldBrlen = b->nodeDepth - a->nodeDepth - 2.0*minV;
    }
    else  /* ancestral fossil */
    {
        y = (b->nodeDepth - minV > v->calibration->max * clockRate) ? (v->calibration->max * clockRate) : (b->nodeDepth - minV);
        x = (a->nodeDepth + minV < v->calibration->min * clockRate) ? (v->calibration->min * clockRate) : (a->nodeDepth + minV);
        oldBrlen = y - x;
    }
    v1 = a->length;
    v2 = u->length;
    v3 = v->length;
    
    /* reassign events for CPP and adjust prior and proposal ratios for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
    {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
        {
            /* get pointers to CPP events */
            nEvents = subParm->nEvents[2*chain+state[chain]];
            position = subParm->position[2*chain+state[chain]];
            rateMultiplier = subParm->rateMult[2*chain+state[chain]];
            n1 = nEvents[a->index];
            n2 = nEvents[u->index];
            n3 = nEvents[v->index];
            if (n2 > 0)
            {
                position[a->index] = (MrBFlt *) SafeRealloc ((void *) position[a->index], (n1+n2) * sizeof (MrBFlt));
                rateMultiplier[a->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[a->index], (n1+n2) * sizeof (MrBFlt));
            }
            for (j=0; j<n1; j++)
                position[a->index][j] *= v1 / (v1+v2);
            for (j=n1; j<n1+n2; j++)
            {
                position[a->index][j] = (position[u->index][j-n1] * v2 + v1) / (v1+v2);
                rateMultiplier[a->index][j] = rateMultiplier[u->index][j-n1];
            }
            nEvents[a->index] = n1+n2;
            nEvents[u->index] = 0;
            if (n2 > 0)
            {
                free (position[u->index]);
                free (rateMultiplier[u->index]);
                position[u->index] = rateMultiplier[u->index] = NULL;
            }
        }   /* end CPP events parm */
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
        {
            /* adjust prior ratio */
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            if (v->length > 0.0)
                (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[v->anc->index], nu*v->length, tk02Rate[v->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[a->anc->index], nu*a->length, tk02Rate[a->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*u->length, tk02Rate[u->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*(a->length+u->length), tk02Rate[a->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = (tk02Rate[a->index] + tk02Rate[b->index]) / 2.0 * (a->length + u->length);
        }   /* end tk02 branch rate parameter */
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
        {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            
            /* adjust prior ratio for old branches */
            if (v->length > 0.0)
                (*lnPriorRatio) -= LnProbGamma(v->length/igrvar, v->length/igrvar, igrRate[v->index]);
            (*lnPriorRatio) -= LnProbGamma(a->length/igrvar, a->length/igrvar, igrRate[a->index]);
            (*lnPriorRatio) -= LnProbGamma(u->length/igrvar, u->length/igrvar, igrRate[u->index]);
            (*lnPriorRatio) += LnProbGamma((a->length+u->length)/igrvar, (a->length+u->length)/igrvar, igrRate[a->index]);
            
            /* adjust effective branch lengths and rates */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = igrRate[a->index] * (a->length + u->length);
        }
    }   /* next subparameter */
    
    /* cut tree */
    a->anc = b;
    if (b->left == u)
        b->left = a;
    else
        b->right = a;
    a->length += u->length;
    a->upDateTi = YES;
    
    /* determine initial direction of move and whether the reverse move would be stopped by constraints */
    if (a->left == NULL || a->isLocked == YES || a->nodeDepth < v->nodeDepth + minV)
    {
        isStartLocked = YES;
        directionUp = NO;
    }
    else
    {
        isStartLocked = NO;
        if (RandomNumber(seed) < 0.5)
            directionUp = YES;
        else
            directionUp = NO;
    }
    
    /* move around in root subtree */
    for (nRootNodes=0; nRootNodes==0 || RandomNumber(seed)<extensionProb; nRootNodes++)
    {
        if (directionUp == YES)
        {   /* going up tree */
            if (a->left == NULL || a->isLocked == YES || a->nodeDepth < v->nodeDepth + minV)
                break;      /* can't go farther */
            topologyHasChanged = YES;
            b = a;
            if (a->left->length < TIME_MIN)
                a = a->right;
            else if (a->right->length < TIME_MIN)
                a = a->left;
            else if (RandomNumber(seed) < 0.5)
                a = a->left;
            else
                a = a->right;
        }
        else
        {   /* going down tree */
            topologyHasChanged = YES;
            if (RandomNumber(seed) < 0.5 || b->anc->anc == NULL || b->isLocked == YES)
            {
                directionUp = YES; /* switch direction */
                /* find sister of a */
                if (b->left == a)
                    a = b->right;
                else
                    a = b->left;
                /* as long as we are moving upwards
                 the cond likes to update will be
                 flagged by the last pass from u to the root */
            }
            else
            {   /* continue down */
                a = b;
                b = b->anc;
                a->upDateCl = YES;
            }
        }
    }
    
    /* determine whether the forward move was or would have been stopped by constraints */
    isStopLocked = NO;
    if (directionUp == YES)
    {
        if (a->left == NULL || a->isLocked == YES || a->nodeDepth < v->nodeDepth + minV)
            isStopLocked = YES;
    }
    
    /* reattach u */
    if (u->left == v)
        u->right = a;
    else
        u->left = a;
    a->anc = u;
    u->anc = b;
    if (b->left == a)
        b->left = u;
    else
        b->right = u;
    
    if (v->length > 0.0)  /* side branch, not anc fossil */
    {
        if (a->nodeDepth > v->nodeDepth)
            newBrlen = b->nodeDepth - a->nodeDepth - 2.0*minV;
        else
            newBrlen = b->nodeDepth - v->nodeDepth - 2.0*minV;
    }
    else  /* ancestral fossil */
    {
        y = (b->nodeDepth - minV > v->calibration->max * clockRate) ? (v->calibration->max * clockRate) : (b->nodeDepth - minV);
        x = (a->nodeDepth + minV < v->calibration->min * clockRate) ? (v->calibration->min * clockRate) : (a->nodeDepth + minV);
        newBrlen = y - x;
    }
    if (newBrlen <= 0.0)
    {
        abortMove = YES;
        return (NO_ERROR);
    }
    
    /* adjust lengths */
    if (v->length > 0.0)  /* side branch, not anc fossil */
    {
        u->nodeDepth = b->nodeDepth - minV - RandomNumber(seed) * newBrlen;
        v->length = u->nodeDepth - v->nodeDepth;
    }
    else  /* ancestral fossil */
    {
        u->nodeDepth = y - RandomNumber(seed) * newBrlen;
        v->nodeDepth = u->nodeDepth;
        v->age = u->age = u->nodeDepth / clockRate;
    }
    u->length = b->nodeDepth - u->nodeDepth;
    a->length = u->nodeDepth - a->nodeDepth;
    
    v3new = v->length;
    v4 = a->length;
    v5 = u->length;
    
    /* adjust events, prior ratio and proposal ratio for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
    {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
        {
            /* reassign events for CPP */
            nEvents = subParm->nEvents[2*chain+state[chain]];
            position = subParm->position[2*chain+state[chain]];
            rateMultiplier = subParm->rateMult[2*chain+state[chain]];
            for (j=0; j<nEvents[a->index]; j++)
            {
                if (position[a->index][j] > v4 / (v4+v5))
                    break;
            }
            n4 = j;
            n5 = nEvents[a->index] - j;
            nEvents[u->index] = n5;
            if (n5 > 0)
            {
                position[u->index] = (MrBFlt *) SafeRealloc ((void *) position[u->index], n5 * sizeof (MrBFlt));
                rateMultiplier[u->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[u->index], n5 * sizeof (MrBFlt));
                for (j=n4; j<nEvents[a->index]; j++)
                {
                    position[u->index][j-n4] = (position[a->index][j] * (v4+v5) - v4) / v5;
                    rateMultiplier[u->index][j-n4] = rateMultiplier[a->index][j];
                }
                if (n4 > 0)
                {
                    position[a->index] = (MrBFlt *) SafeRealloc ((void *) position[a->index], n4 * sizeof (MrBFlt));
                    rateMultiplier[a->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[a->index], n4 * sizeof (MrBFlt));
                    for (j=0; j<n4; j++)
                        position[a->index][j] *= ((v4+v5) / v4);
                }
                else
                {
                    free (position[a->index]);
                    free (rateMultiplier[a->index]);
                    position[a->index] = rateMultiplier[a->index] = NULL;
                }
                nEvents[a->index] = n4;
            }
            else
            {
                for (j=0; j<nEvents[a->index]; j++)
                    position[a->index][j] *= ((v4+v5) / v4);
            }
            
            /* adjust proposal ratio for length change in v branch*/
            (*lnProposalRatio) += n3 * log (v3new / v3);
            
            /* adjust prior ratio for length change */
            lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);
            (*lnPriorRatio) += lambda * (v3 - v3new);
            
            /* update effective branch lengths */
            if (UpdateCppEvolLengths (subParm, oldA, chain) == ERROR)
            {
                abortMove = YES;
                return (NO_ERROR);
            }
            if (UpdateCppEvolLengths (subParm, u, chain) == ERROR)
            {
                abortMove = YES;
                return (NO_ERROR);
            }
        }   /* end cpp events parameter */
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
        {
            /* adjust prior ratio */
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*(a->length+u->length), tk02Rate[a->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[a->anc->index], nu*a->length, tk02Rate[a->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*u->length, tk02Rate[u->index]);
            if (v->length > 0.0)
                (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[v->anc->index], nu*v->length, tk02Rate[v->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = a->length * (tk02Rate[a->index] + tk02Rate[a->anc->index]) / 2.0;
            brlens[v->index] = v->length * (tk02Rate[v->index] + tk02Rate[v->anc->index]) / 2.0;
            brlens[u->index] = u->length * (tk02Rate[u->index] + tk02Rate[u->anc->index]) / 2.0;
        }   /* end tk02 branch rate parameter */
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
        {
            /* adjust prior ratio */
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            (*lnPriorRatio) -= LnProbGamma ((a->length+u->length)/igrvar, (a->length+u->length)/igrvar, igrRate[a->index]);
            (*lnPriorRatio) += LnProbGamma (a->length/igrvar, a->length/igrvar, igrRate[a->index]);
            (*lnPriorRatio) += LnProbGamma (u->length/igrvar, u->length/igrvar, igrRate[u->index]);
            if (v->length > 0.0)
                (*lnPriorRatio) += LnProbGamma (v->length/igrvar, v->length/igrvar, igrRate[v->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[v->index] = igrRate[v->index] * v->length;
            brlens[u->index] = igrRate[u->index] * u->length;
            brlens[a->index] = igrRate[a->index] * a->length;
        }   /* end igr branch rate parameter */
    }   /* next subparameter */
    
    /* set tiprobs update flags */
    a->upDateTi = YES;
    u->upDateTi = YES;
    v->upDateTi = YES;
    
    /* set flags for update of cond likes from u and down to root */
    p = u;
    while (p->anc != NULL)
    {
        p->upDateCl = YES;
        p = p->anc;
    }
    
    /* adjust prior ratio for clock tree */
    if (LogClockTreePriorRatio(param, chain, &x) == ERROR)
        return (ERROR);
    (*lnPriorRatio) += x;
    
    if (topologyHasChanged == YES)
    {
        /* get down pass sequence if tree topology has changed */
        GetDownPass (t);
        /* calculate proposal ratio for tree change */
        (*lnProposalRatio) += log (newBrlen / oldBrlen);
        if (isStartLocked == NO && isStopLocked == YES)
            (*lnProposalRatio) += log (2.0 * (1.0 - extensionProb));
        else if (isStartLocked == YES && isStopLocked == NO)
            (*lnProposalRatio) -= log (2.0 * (1.0 - extensionProb));
    }
    
    /* adjust proposal prob for number movable nodes in new tree */
    numMovableNodesNew=0;
    for (i=0; i<t->nNodes-2; ++i)
    {
        p = t->allDownPass[i];
        if (p->left == NULL)
        {
            if (p->calibration == NULL)
                p->x = NO;
            else
            {
                p->x = YES;
            }
        }
        else
        {
            if (p->left->x == YES && p->right->x == YES)
            {
                p->x = YES;
            }
            else
                p->x = NO;
        }
        a = p->anc->left;
        b = p->anc->right;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL
            || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN) || p->x == NO)
            numMovableNodesNew++;
    }
    
    
    if (numMovableNodesNew!=numMovableNodesOld)
    {
        (*lnProposalRatio) += log (numMovableNodesOld / numMovableNodesNew);
    }
    
#   if defined (DEBUG_ExtSPRClock)
    ShowNodes (t->root, 2, YES);
    printf ("After\nProposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  a: %d  b: %d\n",v->index, u->index, a->index, b->index);
    printf ("No. nodes moved in root subtree: %d\n",nRootNodes);
    printf ("Has topology changed? %d\n",topologyHasChanged);
#   endif
    
    return (NO_ERROR);
}


int Move_ExtSPR (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change topology (and branch lengths) using SPR (unrooted) with extension probability.
       Pick either external or internal branches instead of just internal branches. */
    
    int         i, j, topologyHasChanged, nCrownNodes, nRootNodes, directionLeft, directionUp, 
                isVPriorExp, moveInRoot, isStartConstrained, isStopConstrained;
    MrBFlt      m, x, y, tuning, maxV, minV, extensionProb, brlensExp=0.0;
    TreeNode    *p, *q, *a, *b, *c, *d, *u, *v;
    Tree        *t;
    ModelParams *mp;
    
    /* these parameters should be possible to set by user */
    extensionProb = mvp[0]; /* extension probability */
    tuning = mvp[1];        /* Larget & Simon's tuning parameter lambda */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1] < BRLENS_MAX ? mp->brlensUni[1] : BRLENS_MAX;
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else  /* (param->subParams[0]->paramId == BRLENS_EXP) */
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);
    
    topologyHasChanged = NO;

#   if defined (DEBUG_ExtSPR)
    printf ("Before:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
#   endif
    
    /* pick a random branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes -1))];
        q = p->anc->right;  if (q == p) q = p->anc->left;
        i = j = 0;
        if (p->left == NULL)
            j = 2;
        if (p->anc->anc == NULL)
            i = 2;
        if (p->anc->anc != NULL && (p->anc->isLocked == YES || p->anc->anc->anc == NULL))
            i++;
        if (p->anc->anc != NULL && (q->isLocked == YES || q->left == NULL))
            i++;
        if (p->left != NULL && (p->left->isLocked == YES || p->left->left == NULL))
            j++;
        if (p->left != NULL && (p->right->isLocked == YES || p->right->left == NULL))
            j++;
        } while (i == 2 && j == 2);
    
    /* change in root tree ? */
    if (j == 2)
        moveInRoot = YES;
    else if (i == 2)
        moveInRoot = NO;
    else if (RandomNumber(seed) < 0.5)
        moveInRoot = YES;
    else
        moveInRoot = NO;

    /* determine whether start is constrained on backward move */
    isStartConstrained = isStopConstrained = NO;
    if (moveInRoot == YES && i == 1)
        isStartConstrained = YES;
    else if (moveInRoot == NO && j == 1)
        isStartConstrained = YES;

    /* set up pointers for nodes around the picked branch */
    /* cut the tree into crown, root and attachment part  */
    /* change the relevant lengths in the attachment part */
    v = p;
    u = p->anc;

    /* modify length of middle branch */
    m = v->length;
    x = m * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV) x = minV * minV / x;
        if (x > maxV) x = maxV * maxV / x;
        }
    v->length = x;
    v->upDateTi = YES;

    /* update proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / m);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (m - x);

    /* move around in root subtree */
    nRootNodes = 0;
    if (moveInRoot == YES)
        {
        /* mark nodes in root part */
        /* also determine direction of move in root part */
        if (u->left == v)
            a = u->right;
        else
            a = u->left;
        b = u->anc;
        if (u->anc->anc == NULL || u->isLocked == YES)
            directionUp = YES;
        else if (a->left == NULL || a->isLocked == YES)
            directionUp = NO;
        else if (RandomNumber(seed) < 0.5)
            directionUp = YES;
        else
            directionUp = NO;
        /* cut root part*/
        if (directionUp == NO)
            {
            b = a;  /* switch a and b */
            a = u->anc;
            b->anc = a;
            if (a->left == u)
                a->left = b;
            else
                a->right = b;
            }
        else  // if (directionUp == YES)
            {
            a->anc = b;
            if (b->left == u)
                b->left = a;
            else
                b->right = a;
            y = a->length;
            a->length = u->length;
            u->length = y;
            a->upDateTi = YES;
            u->upDateTi = YES;
            }

        for (nRootNodes=0; RandomNumber(seed)<extensionProb || nRootNodes==0; nRootNodes++) 
            {
            if (directionUp == YES) 
                {   /* going up tree */
                if (a->left == NULL || a->isLocked == YES)
                    break;      /* can't go further */
                b = a;
                if (RandomNumber(seed) < 0.5)
                    a = a->left;
                else
                    a = a->right;
                if (u->isLocked == YES)
                    {
                    b->isLocked = YES;
                    u->isLocked = NO;
                    b->lockID = u->lockID;
                    u->lockID = 0;
                    }
                }
            else  // directionUp == NO
                {   /* going down tree */
                if (a->anc == NULL || u->isLocked == YES)
                    break;      /* can't go further */
                if (RandomNumber(seed) < 0.5)
                    {
                    directionUp = YES; /* switch direction */
                    /* find sister of a */
                    if (a->left == b) 
                        {
                        b = a;
                        a = a->right;
                        }
                    else 
                        {  
                        b = a;
                        a = a->left;
                        }
                    /* as long as we are moving upwards, the cond likes to update
                       will be flagged by the last pass from u to the root */
                    }   
                else 
                    {   /* continue down */
                    b = a;
                    a = a->anc;
                    b->upDateCl = YES; 
                    if (b->isLocked == YES)
                        {
                        u->isLocked = YES;
                        b->isLocked = NO;
                        u->lockID = b->lockID;
                        b->lockID = 0;
                        }
                    }
                }
            }
            
        topologyHasChanged = YES;
        /* check whether stop is constrained */
        if (directionUp == YES && (a->left == NULL || a->isLocked == YES))
            isStopConstrained = YES;
        if (directionUp == NO  && (a->anc  == NULL || u->isLocked == YES))
            isStopConstrained = YES;
        
        /* modify branch length */
        m = u->length;
        x = m * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        u->length = x;
        u->upDateTi = YES;

        /* update proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / m);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (m - x);

        /* combine the subtrees */
        if (directionUp == YES)
            {
            u->anc = b;
            if (u->left == v)
                u->right = a;
            else 
                u->left = a;
            a->anc = u;
            if (b->left == a)
                b->left = u;
            else
                b->right = u;
            }
        else  // if (directionUp == NO)
            {
            u->anc = a;
            if (u->left == v)
                u->right = b;
            else 
                u->left = b;
            b->anc = u;
            if (a->left == b)
                a->left = u;
            else
                a->right = u;
            /* the modified branch contained in u->length will have to be moved to b->length to enable back move
               BUT if we haven't moved, it is better to keep it in place (necessary for rooted trees) */
            y = u->length;
            u->length = b->length;
            b->length = y;
            b->upDateTi = YES;
            u->upDateTi = YES;
            }
        }

    /* move around in crown subtree */
    nCrownNodes = 0;
    if (moveInRoot == NO)
        {
        /* set up pointers for crown part */
        /* also determine direction of move in crown part */
        if (v->right->left == NULL || v->right->isLocked == YES)
            directionLeft = YES;
        else if (v->left->left == NULL || v->left->isLocked == YES)
            directionLeft = NO;
        else if (RandomNumber(seed) < 0.5)
            directionLeft = YES;
        else
            directionLeft = NO;
        if (directionLeft == YES)
            {
            c = v->left;
            d = v->right;
            }
        else
            {
            c = v->right;
            d = v->left;
            }

        /* store brlen nodes and brlen to move */
        x = c->length;

        /* cut and reconnect crown part */
        c->anc = d;
        d->anc = c;
    
        for (nCrownNodes=0; RandomNumber(seed)<extensionProb || nCrownNodes==0; nCrownNodes++) 
            {
            if (c->left == NULL || c->isLocked == YES)
                break;  /* can't go further */
            if (RandomNumber(seed) < 0.5)
                {
                /* rotate c anticlockwise - prepare pointers for move left */
                c->anc = c->left;  /* the root will be in the direction we are heading */
                c->left = c->right;
                c->right = d;
                }
            else 
                {
                /* rotate c clockwise - prepare pointers for move right */
                c->anc = c->right;  /* the root will be in the direction we are heading */
                c->right = c->left;
                c->left = d;  
                }
            /* OK - let's move! c->anc points in the right direction
               don't forget to move the branch lengths as well */
            d = c;
            c = c->anc;
            d->length = c->length;  /* also rotate other info ?? */
            d->upDateCl = YES; 
            d->upDateTi = YES;
            }
            
        topologyHasChanged = YES;
        /* check if stop constrained */
        if (c->left == NULL || c->isLocked == YES)
            isStopConstrained = YES;
        
        /* combine the subtrees */
        c->anc = v;
        d->anc = v;
        if (directionLeft == YES)
            {
            v->left = c;
            v->right = d;
            }
        else
            {
            v->left = d;
            v->right = c;
            }
        
        /* the dangling branch is inserted in reverted position such that the back move will be possible
           if we have moved around in crown subtree otherwise it is left in its original position */
        d->length = x;
        
        /* modify branch length */
        m = d->length;
        x = m * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        d->length = x;
        d->upDateTi = YES;

        /* update proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / m);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (m - x);
        }

    /* adjust proposal ratio for constraints */
    if (isStartConstrained == NO && isStopConstrained == YES)
        (*lnProposalRatio) += log (2.0 * (1.0 - extensionProb));
    else if (isStartConstrained == YES && isStopConstrained == NO)
        (*lnProposalRatio) -= log (2.0 * (1.0 - extensionProb));

    /* set flags for update of cond likes from v and down to root */
    p = v;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);

#   if defined (DEBUG_ExtSPR)
    printf ("After:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
    printf ("Proposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  c: %d  d: %d  a: %d  b: %d\n",v->index, u->index,
            c->index, d->index, a->index, b->index);
    printf ("No. nodes moved in root subtree: %d\n",nRootNodes);
    printf ("No. nodes moved in crown subtree: %d\n",nCrownNodes);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    getchar();
#   endif

    return (NO_ERROR);
}


int Move_ExtSPR1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change topology (and branch lengths) using SPR (unrooted) with extension probability. */
    /* Pick only internal branches. For a description, see Lakner et al. (2008). */
    
    int         i, j, topologyHasChanged, nCrownNodes, nRootNodes, directionLeft, directionUp, 
                isVPriorExp, moveInRoot, isStartConstrained, isStopConstrained;
    MrBFlt      m, x, y, tuning, maxV, minV, extensionProb, brlensExp=0.0;
    TreeNode    *p, *a, *b, *c, *d, *u, *v, *brlenNode[7];
    Tree        *t;
    ModelParams *mp;
    
    /* these parameters should be possible to set by user */
    extensionProb = mvp[0]; /* extension probability */
    tuning = mvp[1];        /* Larget & Simon's tuning parameter lambda */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1] < BRLENS_MAX ? mp->brlensUni[1] : BRLENS_MAX;
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);
    
    topologyHasChanged = NO;

#   if defined (DEBUG_ExtSPR)
    printf ("Before:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
#   endif
    
    /* pick an internal branch that is free to move in either end
       (i and j keep track of number of locked directions) */
    do  {
        p = t->intDownPass[(int)(RandomNumber(seed) * (t->nIntNodes-1))];
        if (p->anc->left == p)
            a = p->anc->right;
        else
            a = p->anc->left;
        i = j = 0;
        if (a->isLocked == YES || a->left == NULL)
            i++;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL)
            i++;
        if (p->left->isLocked == YES || p->left->left == NULL)
            j++;
        if (p->right->isLocked == YES || p->right->left == NULL)
            j++;
        } while (i == 2 && j == 2);
        
    /* set up pointers for nodes around the picked branch */
    /* cut the tree into crown, root and attachment part */
    /* change the relevant lengths in the attachment part */
    /* the lengths of a and v are automatically contained in the */
    /* "attachment" part but the length of c has to be stored in x */
    v = p;
    u = p->anc;

    /* store brlen node */
    brlenNode[3] = v;

    /* change in root tree ? */
    if (j == 2)
        moveInRoot = YES;
    else if (i == 2)
        moveInRoot = NO;
    else if (RandomNumber(seed) < 0.5)
        moveInRoot = YES;
    else
        moveInRoot = NO;

    /* determine whether start is constrained on backward move */
    isStartConstrained = isStopConstrained = NO;
    if (moveInRoot == YES && i == 1)
        isStartConstrained = YES;
    else if (moveInRoot == NO && j == 1)
        isStartConstrained = YES;

    /* set up pointers for crown part */
    /* also determine direction of move in crown part */
    if (v->right->left == NULL || v->right->isLocked == YES)
        directionLeft = YES;
    else if (v->left->left == NULL || v->left->isLocked == YES)
        directionLeft = NO;
    else if (RandomNumber(seed) < 0.5)
        directionLeft = YES;
    else
        directionLeft = NO;
    if (directionLeft == YES)
        {
        c = v->left;
        d = v->right;
        }
    else
        {
        c = v->right;
        d = v->left;
        }

    /* store brlen nodes and brlen to move */
    brlenNode[0] = d;
    brlenNode[1] = c;
    x = c->length;

    /* cut and reconnect crown part */
    c->anc = d;
    d->anc = c;
    
    /* mark nodes in root part */
    /* also determine direction of move in root part */
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    if (u->anc->anc == NULL || u->isLocked == YES)
        directionUp = YES;
    else if (a->left == NULL || a->isLocked == YES)
        directionUp = NO;
    else if (RandomNumber(seed) < 0.5)
        directionUp = YES;
    else
        directionUp = NO;
    if (directionUp == NO)
        {
        /* switch a and b */
        b = a;
        a = u->anc;
        }

    /* store brlen nodes */
    if (directionUp == YES)
        {
        brlenNode[4] = u;
        brlenNode[5] = a;
        }
    else
        {
        brlenNode[4] = b;
        brlenNode[5] = u;
        }

    /* cut root part*/
    if (directionUp == NO)
        {
        b->anc = a;
        if (a->left == u)
            a->left = b;
        else
            a->right = b;
        }
    else 
        {
        a->anc = b;
        if (b->left == u)
            b->left = a;
        else
            b->right = a;
        y = a->length;
        a->length = u->length;
        u->length = y;
        a->upDateTi = YES;
        u->upDateTi = YES;
        }

    /* move around in root subtree */
    nRootNodes = 0;
    if (moveInRoot == YES)
        {
        for (nRootNodes=0; RandomNumber(seed)<extensionProb || nRootNodes==0; nRootNodes++) 
            {
            if (directionUp == YES) 
                {   /* going up tree */
                if (a->left == NULL || a->isLocked == YES)
                    break;      /* can't go further */
                topologyHasChanged = YES;
                b = a;
                if (RandomNumber(seed) < 0.5)
                    a = a->left;
                else
                    a = a->right;
                if (u->isLocked == YES)
                    {
                    b->isLocked = YES;
                    u->isLocked = NO;
                    b->lockID = u->lockID;
                    u->lockID = 0;
                    }
                }
            else 
                {   /* going down tree */
                if (a->anc == NULL || u->isLocked == YES)
                    break;      /* can't go further */
                topologyHasChanged = YES;
                if (RandomNumber(seed)<0.5) 
                    {
                    directionUp = YES; /* switch direction */
                    /* find sister of a */
                    if (a->left == b) 
                        {
                        b = a;
                        a = a->right;
                        }
                    else 
                        {  
                        b = a;
                        a = a->left;
                        }
                    /* as long as we are moving upwards
                    the cond likes to update will be
                    flagged by the last pass from u to the root */
                    }   
                else 
                    {   /* continue down */
                    b = a;
                    a = a->anc;
                    b->upDateCl = YES; 
                    if (b->isLocked == YES)
                        {
                        u->isLocked = YES;
                        b->isLocked = NO;
                        u->lockID = b->lockID;
                        b->lockID = 0;
                        }
                    }
                }
            }
        /* check whether stop is constrained */
        if (directionUp == YES) 
            {
            if (a->left == NULL || a->isLocked == YES) 
                isStopConstrained = YES;
            }
        else 
            {
            if (a->anc  == NULL || u->isLocked == YES)
                isStopConstrained = YES;
            }
        }

    /* store brlen nodes */
    if (nRootNodes > 0)
        {
        if (directionUp == YES)
            {
            brlenNode[6] = a;
            brlenNode[5] = u;
            }
        else
            {
            brlenNode[6] = u;
            brlenNode[5] = b;
            }
        }

    /* move around in crown subtree */
    nCrownNodes = 0;
    if (moveInRoot == NO)       
        {
        for (nCrownNodes=0; RandomNumber(seed)<extensionProb || nCrownNodes==0; nCrownNodes++) 
            {
            if (c->left == NULL || c->isLocked == YES)
                break;  /* can't go further */
            topologyHasChanged = YES;
            if (RandomNumber(seed) < 0.5) 
                {
                /* rotate c anticlockwise - prepare pointers for move left */
                c->anc = c->left;  /* the root will be in the direction we are heading */
                c->left = c->right;
                c->right = d;
                }
            else 
                {
                /* rotate c clockwise - prepare pointers for move right */
                c->anc = c->right;  /* the root will be in the direction we are heading */
                c->right = c->left;
                c->left = d;  
                }
            /* OK - let's move!; c->anc points in the right direction
            don't forget to move the branch lengths as well */
            d = c;
            c = c->anc;
            d->length = c->length;
            d->upDateCl = YES; 
            d->upDateTi = YES;
            }
        /* check if stop constrained */
        if (c->left == NULL || c->isLocked == YES)
            isStopConstrained = YES;
        }

    /* store brlen nodes */
    if (nCrownNodes > 0)
        {
        brlenNode[2] = c;
        brlenNode[1] = d;
        }

    /* adjust proposal ratio for constraints */
    if (isStartConstrained == NO && isStopConstrained == YES)
        (*lnProposalRatio) += log (2.0 * (1.0 - extensionProb));
    else if (isStartConstrained == YES && isStopConstrained == NO)
        (*lnProposalRatio) -= log (2.0 * (1.0 - extensionProb));

    /* combine the subtrees */
    c->anc = v;
    d->anc = v;
    if (directionLeft == YES) 
        {
        v->left = c;
        v->right = d;
        }
    else 
        {
        v->left = d;
        v->right = c;
        }

    /* the dangling branch is inserted in reverted position
       such that the back move will be possible
       if we have moved around in crown subtree
       otherwise it is left in its original position */
    if (nCrownNodes > 0)
        {
        d->length = x;
        d->upDateTi = YES;
        }
    else
        {
        c->length = x;
        }

    if (directionUp == YES) 
        {
        u->anc = b;
        if (u->left == v)
            u->right = a;
        else 
            u->left = a;
        a->anc = u;
        if (b->left == a)
            b->left = u;
        else
            b->right = u;
        /* the dangling branch is contained in u->length
           and will automatically be inserted in the right position
           to enable the back move regardless of whether it was
           initially directed upwards or downwards
           BUT if we haven't moved in root subtree, it is advantageous (necessary
           for rooted trees) to avoid switching branches, which occurs otherwise
           if directionUp == YES */
        if (nRootNodes == 0) 
            {
            x = u->length;
            u->length = a->length;
            a->length = x;
            a->upDateTi = NO;
            u->upDateTi = NO;
            }
        }
    else 
        {
        u->anc = a;
        if (u->left == v)
            u->right = b;
        else 
            u->left = b;
        b->anc = u;
        if (a->left == b)
            a->left = u;
        else
            a->right = u;
        /* the modified branch contained in u->length will have
           to be moved to b->length to enable back move
           BUT if we haven't moved, it is better to keep it in place
           (necessary for rooted trees) */
        if (nRootNodes > 0) 
            {
            x = u->length;
            u->length = b->length;
            b->length = x;
            b->upDateTi = YES;
            u->upDateTi = YES;
            }
        }
    
    /* modify branch lengths */
    /* first modify length of middle branch */
    m = brlenNode[3]->length;
    x = m * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV)
            x = minV * minV / x;
        else if (x > maxV)
            x = maxV * maxV / x;
        }
    brlenNode[3]->length = x;
    brlenNode[3]->upDateTi = YES;

    /* update proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / m);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (m - x);

    if (moveInRoot == NO)
        {
        /* if no move in crown, then select randomly, otherwise always the moved branch */
        if (nCrownNodes == 0 && RandomNumber(seed) < 0.5)
            p = brlenNode[0];
        else
            p = brlenNode[1];

        /* modify branch length */
        m = p->length;
        x = m * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV)
                x = minV * minV / x;
            else if (x > maxV)
                x = maxV * maxV / x;
            }
        p->length = x;
        p->upDateTi = YES;

        /* update proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / m);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (m - x);
        }
            
    if (moveInRoot == YES)
        {
        /* if no move in root, then select randomly, otherwise always the moved branch */
        if (nRootNodes == 0 && RandomNumber(seed) < 0.5)
            p = brlenNode[4];
        else
            p = brlenNode[5];
        
        /* modify branch length but not if 'root' branch in rooted tree */
        if (t->isRooted == NO || p->anc->anc != NULL)
            {
            m = p->length;
            x = m * exp(tuning * (RandomNumber(seed) - 0.5));
            while (x < minV || x > maxV)
                {
                if (x < minV)
                    x = minV * minV / x;
                else if (x > maxV)
                    x = maxV * maxV / x;
                }
            p->length = x;
            p->upDateTi = YES;

            /* update proposal and prior ratio based on length modification */
            (*lnProposalRatio) += log (x / m);
            if (isVPriorExp == YES)
                (*lnPriorRatio) += brlensExp * (m - x); 
            }
        }

    /* set flags for update of cond likes from v and down to root */
    p = v;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);

#   if defined (DEBUG_ExtSPR)
    printf ("After:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
    printf ("Proposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  c: %d  d: %d  a: %d  b: %d\n",v->index, u->index,
            c->index, d->index, a->index, b->index);
    printf ("No. nodes moved in root subtree: %d\n",nRootNodes);
    printf ("No. nodes moved in crown subtree: %d\n",nCrownNodes);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    getchar();
#   endif

    return (NO_ERROR);
}


int Move_ExtSPRClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change branch lengths and topology (potentially) using SPR-type move 
       with extension probability (rather than window, attachment rate or similar).
       The move is Metropolized, which should improve mixing. However, this means 
       that it must be combined with a node slider move to be efficient.

       The move picks a branch and then moves its lower attachment point 
       from its original position, one node at a time, with
       a probability determined by the extensionProb parameter. This is
       done in a way consistent with the clock constraints and any locked
       nodes there might be in the tree. The lower attachment point is 
       minimally moved one node away.
       
       On the ending branch, the attachment point is reinserted randomly
       along the branch (below the minimum age of the node). */
    
    int         i, j, topologyHasChanged=NO, isStartLocked=NO, isStopLocked=NO, nRootNodes, directionUp,
                n1=0, n2=0, n3=0, n4=0, n5=0, *nEvents;
    MrBFlt      x, y, oldBrlen=0.0, newBrlen=0.0, extensionProb, igrvar, *igrRate=NULL,
                v1=0.0, v2=0.0, v3=0.0, v4=0.0, v5=0.0, v3new=0.0, lambda, *tk02Rate=NULL,
                **position=NULL, **rateMultiplier=NULL, *brlens, nu, minV, clockRate;
    TreeNode    *p, *a, *b, *u, *v, *oldA;
    Tree        *t;
    ModelInfo   *m;
    Param       *subParm;

    extensionProb = mvp[0]; /* extension probability */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);
 
    /* get model params and model info */
    m = &modelSettings[param->relParts[0]];
    
    /* get clock rate */
    clockRate = *GetParamVals (m->clockRate, chain, state[chain]);

    /* get min and max branch lengths in relative time and substitution units */
    minV = BRLENS_MIN;

#   if defined (DEBUG_ExtSPRClock)
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
#   endif
    
    /* pick a branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes - 2))];
        a = p->anc->left;
        b = p->anc->right;
        }
    while (p->anc->isLocked == YES || p->anc->anc->anc == NULL
           || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN)
           || (p->length < TIME_MIN && p->calibration->prior == fixed));
    /* skip constraints, siblings of root (and root); and consider ancestral fossils in fbd tree */
    
    /* set up pointers for nodes around the picked branch */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    oldA = a;

    /* record branch length for insertion in back move */
    if (v->length > 0.0)  /* side branch, not anc fossil */
        {
        if (v->nodeDepth > a->nodeDepth)
            oldBrlen = b->nodeDepth - v->nodeDepth - 2.0*minV;
        else
            oldBrlen = b->nodeDepth - a->nodeDepth - 2.0*minV;
        }
    else  /* ancestral fossil */
        {
        y = (b->nodeDepth - minV > v->calibration->max * clockRate) ? (v->calibration->max * clockRate) : (b->nodeDepth - minV);
        x = (a->nodeDepth + minV < v->calibration->min * clockRate) ? (v->calibration->min * clockRate) : (a->nodeDepth + minV);
        oldBrlen = y - x;
        }
    v1 = a->length;
    v2 = u->length;
    v3 = v->length;

    /* reassign events for CPP and adjust prior and proposal ratios for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
        {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            /* get pointers to CPP events */
            nEvents = subParm->nEvents[2*chain+state[chain]];
            position = subParm->position[2*chain+state[chain]];
            rateMultiplier = subParm->rateMult[2*chain+state[chain]];
            n1 = nEvents[a->index];
            n2 = nEvents[u->index];
            n3 = nEvents[v->index];
            if (n2 > 0)
                {
                position[a->index] = (MrBFlt *) SafeRealloc ((void *) position[a->index], (n1+n2) * sizeof (MrBFlt));
                rateMultiplier[a->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[a->index], (n1+n2) * sizeof (MrBFlt));
                }
            for (j=0; j<n1; j++)
                position[a->index][j] *= v1 / (v1+v2);
            for (j=n1; j<n1+n2; j++)
                {
                position[a->index][j] = (position[u->index][j-n1] * v2 + v1) / (v1+v2);
                rateMultiplier[a->index][j] = rateMultiplier[u->index][j-n1];
                }
            nEvents[a->index] = n1+n2;
            nEvents[u->index] = 0;
            if (n2 > 0)
                {
                free (position[u->index]);
                free (rateMultiplier[u->index]);
                position[u->index] = rateMultiplier[u->index] = NULL;
                }
            }   /* end CPP events parm */
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            /* adjust prior ratio */
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            if (v->length > 0.0)
                (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[v->anc->index], nu*v->length, tk02Rate[v->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[a->anc->index], nu*a->length, tk02Rate[a->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*u->length, tk02Rate[u->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*(a->length+u->length), tk02Rate[a->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = (tk02Rate[a->index] + tk02Rate[b->index]) / 2.0 * (a->length + u->length);
            }   /* end tk02 branch rate parameter */
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);

             /* adjust prior ratio for old branches */
            if (v->length > 0.0)
                (*lnPriorRatio) -= LnProbGamma(v->length/igrvar, v->length/igrvar, igrRate[v->index]);
            (*lnPriorRatio) -= LnProbGamma(a->length/igrvar, a->length/igrvar, igrRate[a->index]);
            (*lnPriorRatio) -= LnProbGamma(u->length/igrvar, u->length/igrvar, igrRate[u->index]);
            (*lnPriorRatio) += LnProbGamma((a->length+u->length)/igrvar, (a->length+u->length)/igrvar, igrRate[a->index]);

            /* adjust effective branch lengths and rates */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = igrRate[a->index] * (a->length + u->length);
            }
        }   /* next subparameter */

    /* cut tree */
    a->anc = b;
    if (b->left == u)
        b->left = a;
    else
        b->right = a;
    a->length += u->length;
    a->upDateTi = YES;

    /* determine initial direction of move and whether the reverse move would be stopped by constraints */
    if (a->left == NULL || a->isLocked == YES || a->nodeDepth < v->nodeDepth + minV)
        {
        isStartLocked = YES;
        directionUp = NO;
        }
    else
        {
        isStartLocked = NO;
        if (RandomNumber(seed) < 0.5)
            directionUp = YES;
        else
            directionUp = NO;
        }
        
    /* move around in root subtree */
    for (nRootNodes=0; nRootNodes==0 || RandomNumber(seed)<extensionProb; nRootNodes++) 
        {
        if (directionUp == YES) 
            {   /* going up tree */
            if (a->left == NULL || a->isLocked == YES || a->nodeDepth < v->nodeDepth + minV)
                break;      /* can't go farther */
            topologyHasChanged = YES;
            b = a;
            if (a->left->length < TIME_MIN)
                a = a->right;
            else if (a->right->length < TIME_MIN)
                a = a->left;
            else if (RandomNumber(seed) < 0.5)
                a = a->left;
            else
                a = a->right;
            }
        else 
            {   /* going down tree */
            topologyHasChanged = YES;
            if (RandomNumber(seed) < 0.5 || b->anc->anc == NULL || b->isLocked == YES)
                {
                directionUp = YES; /* switch direction */
                /* find sister of a */
                if (b->left == a)
                    a = b->right;
                else
                    a = b->left;
                /* as long as we are moving upwards
                the cond likes to update will be
                flagged by the last pass from u to the root */
                }   
            else 
                {   /* continue down */
                a = b;
                b = b->anc;
                a->upDateCl = YES;
                }
            }
        }
        
    /* determine whether the forward move was or would have been stopped by constraints */
    isStopLocked = NO;
    if (directionUp == YES)
        {
        if (a->left == NULL || a->isLocked == YES || a->nodeDepth < v->nodeDepth + minV)
            isStopLocked = YES;
        }

    /* reattach u */
    if (u->left == v)
        u->right = a;
    else
        u->left = a;
    a->anc = u;
    u->anc = b;
    if (b->left == a)
        b->left = u;
    else
        b->right = u;

    if (v->length > 0.0)  /* side branch, not anc fossil */
        {
        if (a->nodeDepth > v->nodeDepth)
            newBrlen = b->nodeDepth - a->nodeDepth - 2.0*minV;
        else
            newBrlen = b->nodeDepth - v->nodeDepth - 2.0*minV;
        }
    else  /* ancestral fossil */
        {
        y = (b->nodeDepth - minV > v->calibration->max * clockRate) ? (v->calibration->max * clockRate) : (b->nodeDepth - minV);
        x = (a->nodeDepth + minV < v->calibration->min * clockRate) ? (v->calibration->min * clockRate) : (a->nodeDepth + minV);
        newBrlen = y - x;
        }
    if (newBrlen <= 0.0)
        {
        abortMove = YES;
        return (NO_ERROR);
        }
    
    /* adjust lengths */
    if (v->length > 0.0)  /* side branch, not anc fossil */
        {
        u->nodeDepth = b->nodeDepth - minV - RandomNumber(seed) * newBrlen;
        v->length = u->nodeDepth - v->nodeDepth;
        }
    else  /* ancestral fossil */
        {
        u->nodeDepth = y - RandomNumber(seed) * newBrlen;
        v->nodeDepth = u->nodeDepth;
        v->age = u->age = u->nodeDepth / clockRate;
        }
    u->length = b->nodeDepth - u->nodeDepth;
    a->length = u->nodeDepth - a->nodeDepth;
    
    v3new = v->length;
    v4 = a->length;
    v5 = u->length;

    /* adjust events, prior ratio and proposal ratio for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
        {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            /* reassign events for CPP */
            nEvents = subParm->nEvents[2*chain+state[chain]];
            position = subParm->position[2*chain+state[chain]];
            rateMultiplier = subParm->rateMult[2*chain+state[chain]];
            for (j=0; j<nEvents[a->index]; j++)
                {
                if (position[a->index][j] > v4 / (v4+v5))
                    break;
                }
            n4 = j;
            n5 = nEvents[a->index] - j;
            nEvents[u->index] = n5;
            if (n5 > 0)
                {
                position[u->index] = (MrBFlt *) SafeRealloc ((void *) position[u->index], n5 * sizeof (MrBFlt));
                rateMultiplier[u->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[u->index], n5 * sizeof (MrBFlt));            
                for (j=n4; j<nEvents[a->index]; j++)
                    {
                    position[u->index][j-n4] = (position[a->index][j] * (v4+v5) - v4) / v5;
                    rateMultiplier[u->index][j-n4] = rateMultiplier[a->index][j];
                    }
                if (n4 > 0)
                    {
                    position[a->index] = (MrBFlt *) SafeRealloc ((void *) position[a->index], n4 * sizeof (MrBFlt));
                    rateMultiplier[a->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[a->index], n4 * sizeof (MrBFlt));
                    for (j=0; j<n4; j++)
                        position[a->index][j] *= ((v4+v5) / v4);
                    }
                else
                    {
                    free (position[a->index]);
                    free (rateMultiplier[a->index]);
                    position[a->index] = rateMultiplier[a->index] = NULL;
                    }
                nEvents[a->index] = n4;
                }
            else
                {
                for (j=0; j<nEvents[a->index]; j++)
                    position[a->index][j] *= ((v4+v5) / v4);
                }

            /* adjust proposal ratio for length change in v branch*/
            (*lnProposalRatio) += n3 * log (v3new / v3);

            /* adjust prior ratio for length change */
            lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);
            (*lnPriorRatio) += lambda * (v3 - v3new);

            /* update effective branch lengths */
            if (UpdateCppEvolLengths (subParm, oldA, chain) == ERROR)
                {
                abortMove = YES;
                return (NO_ERROR);
                }
            if (UpdateCppEvolLengths (subParm, u, chain) == ERROR)
                {
                abortMove = YES;
                return (NO_ERROR);
                }
            }   /* end cpp events parameter */
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            /* adjust prior ratio */
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*(a->length+u->length), tk02Rate[a->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[a->anc->index], nu*a->length, tk02Rate[a->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*u->length, tk02Rate[u->index]);
            if (v->length > 0.0)
                (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[v->anc->index], nu*v->length, tk02Rate[v->index]);

            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = a->length * (tk02Rate[a->index] + tk02Rate[a->anc->index]) / 2.0;
            brlens[v->index] = v->length * (tk02Rate[v->index] + tk02Rate[v->anc->index]) / 2.0;
            brlens[u->index] = u->length * (tk02Rate[u->index] + tk02Rate[u->anc->index]) / 2.0;
            }   /* end tk02 branch rate parameter */
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            /* adjust prior ratio */
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            (*lnPriorRatio) -= LnProbGamma ((a->length+u->length)/igrvar, (a->length+u->length)/igrvar, igrRate[a->index]);
            (*lnPriorRatio) += LnProbGamma (a->length/igrvar, a->length/igrvar, igrRate[a->index]);
            (*lnPriorRatio) += LnProbGamma (u->length/igrvar, u->length/igrvar, igrRate[u->index]);
            if (v->length > 0.0)
                (*lnPriorRatio) += LnProbGamma (v->length/igrvar, v->length/igrvar, igrRate[v->index]);

            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[v->index] = igrRate[v->index] * v->length;
            brlens[u->index] = igrRate[u->index] * u->length;
            brlens[a->index] = igrRate[a->index] * a->length;
            }   /* end igr branch rate parameter */
        }   /* next subparameter */

    /* set tiprobs update flags */
    a->upDateTi = YES;
    u->upDateTi = YES;
    v->upDateTi = YES;

    /* set flags for update of cond likes from u and down to root */
    p = u;
    while (p->anc != NULL)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }

    /* adjust prior ratio for clock tree */
    if (LogClockTreePriorRatio(param, chain, &x) == ERROR)
        return (ERROR);
    (*lnPriorRatio) += x;

    if (topologyHasChanged == YES)
        {
        /* get down pass sequence if tree topology has changed */
        GetDownPass (t);
        /* calculate proposal ratio for tree change */
        (*lnProposalRatio) += log (newBrlen / oldBrlen);
        if (isStartLocked == NO && isStopLocked == YES)
            (*lnProposalRatio) += log (2.0 * (1.0 - extensionProb));
        else if (isStartLocked == YES && isStopLocked == NO)
            (*lnProposalRatio) -= log (2.0 * (1.0 - extensionProb));
        }

#   if defined (DEBUG_ExtSPRClock)
    ShowNodes (t->root, 2, YES);
    printf ("After\nProposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  a: %d  b: %d\n",v->index, u->index, a->index, b->index);
    printf ("No. nodes moved in root subtree: %d\n",nRootNodes);
    printf ("Has topology changed? %d\n",topologyHasChanged);
#   endif

    return (NO_ERROR);
}


int Move_ExtSS (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change branch lengths and topology (potentially) using Subtree Swapping (unrooted) 
       with extension probability.

       This move type picks two subtrees and swaps their position. Like the SPR and TBR,
       it is a superset of the NNI but there are some interesting differences. With the
       SPR and TBR, it is not possible to go between all five-tip trees in a single
       step. For instance, going from ((1,2),3,(4,5)) to ((1,5),3,(4,2)) requires two
       steps. The SS move can go between all pairs of five-tip trees in a single step.
       Some six-tip tree pairs will require two steps though.

       Unlike the published version of the move (Lakner et al, Syst Bio), this version
       does _not_ multiply all branch lengths between the subtrees.
       
       */
    
    int         i, numFree, topologyHasChanged, nCrownNodes, nRootNodes, directionLeft, directionUp, 
                isVPriorExp, moveInRoot;
    MrBFlt      m, x, tuning, maxV, minV, extensionProb, brlensExp=0.0;
    TreeNode    *p, *q, *a, *b, *c, *d, *u, *v;
    Tree        *t;
    ModelParams *mp;

    (*lnPriorRatio) = (*lnProposalRatio) = 0.0;

    /* these parameters should be possible to set by user */
    extensionProb = mvp[0]; /* extension probability */
    tuning = mvp[1];        /* Larget & Simon's tuning parameter lambda */
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1];
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);

    topologyHasChanged = NO;

    /* unmark all tree */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;
        }

    /* pick a branch */
    do
        {
        p = t->allDownPass[(int)(RandomNumber(seed) * t->nNodes)];
        } while (p->anc == NULL);
        
    /* set up pointers for nodes around the picked branch */
    v = p;
    u = p->anc;

    /* check the possible move directions */
    numFree = 0;
    if (v->left != NULL && v->left->isLocked == NO)
        numFree ++;
    if (v->right != NULL && v->right->isLocked == NO)
        numFree++;
    if (u->anc != NULL && u->isLocked == NO)
        numFree++;
    if (u->left == v)
        {
        if (u->right != NULL && u->right->isLocked == NO)
            numFree++;
        }
    else
        {
        if (u->left != NULL && u->left->isLocked == NO)
            numFree++;
        }

    /* select one of them randomly */
    i = (int) (RandomNumber(seed) * numFree) + 1;
    numFree = 0;
    a = b = c = d = p;
    directionLeft = directionUp = moveInRoot = NO;
    if (v->left != NULL && v->left->isLocked == NO)
        {
        numFree ++;
        if (i == numFree)
            {
            moveInRoot = NO;
            directionLeft = YES;
            c = v->left;
            d = v;
            }
        }
    if (v->right != NULL && v->right->isLocked == NO)
        {
        numFree ++;
        if (i == numFree)
            {
            moveInRoot = NO;
            directionLeft = NO;
            c = v->right;
            d = v;
            }
        }
    if (u->anc != NULL && u->isLocked == NO)
        {
        numFree ++;
        if (i == numFree)
            {
            moveInRoot = YES;
            directionUp = NO;
            a = u->anc;
            b = u;
            }
        }
    if (u->left == v)
        {
        if (u->right != NULL && u->right->isLocked == NO)
            {
            numFree ++;
            if (i == numFree)
                {
                moveInRoot = YES;
                directionUp = YES;
                a = u->right;
                b = u;
                }
            }
        }
    else
        {
        if (u->left != NULL && u->left->isLocked == NO)
            {
            numFree ++;
            if (i == numFree)
                {
                moveInRoot = YES;
                directionUp = YES;
                a = u->left;
                b = u;
                }
            }
        }

#   if defined (DEBUG_ExtSS)
    printf ("Before:\n");
    ShowNodes (t->root, 2, NO);
    printf ("v: %d  u: %d  c: %d  d: %d  a: %d  b: %d\n",v->index, u->index, 
        c->index, d->index, a->index, b->index);
    printf ("directionUp = %d -- directionLeft = %d -- moveInRoot = %d\n", directionUp, directionLeft, moveInRoot);
    getchar();
#   endif
    
    /* move around and potentially swap in root subtree */
    nRootNodes = 0;
    if (moveInRoot == YES)
        {
        for (nRootNodes=0; nRootNodes==0 || RandomNumber(seed)<extensionProb; nRootNodes++) 
            {
            if (directionUp == YES) 
                {   /* going up tree */
                if (a->left == NULL || a->isLocked == YES)
                    break;      /* can't go further */
                topologyHasChanged = YES;
                b = a;
                if (RandomNumber(seed) < 0.5)
                    a = a->left;
                else
                    a = a->right;
                }
            else 
                {   /* going down tree */
                if (a->anc == NULL || a->isLocked == YES)
                    break;      /* can't go further */
                topologyHasChanged = YES;
                b->marked = YES;

                if (RandomNumber(seed) < 0.5) 
                    {
                    directionUp = YES; /* switch direction */
                    /* find sister of a */
                    if (a->left == b) 
                        {
                        b = a;
                        a = a->right;
                        }
                    else 
                        {  
                        b = a;
                        a = a->left;
                        }
                    }   
                else 
                    {   /* continue down */
                    b = a;
                    a = a->anc;
                    }
                }
            }
        /* swap the root subtrees */
        if (nRootNodes > 0)
            {
            if (directionUp == YES)
                {
                v->anc = b;
                a->anc = u;
                if (b->left == a)
                    b->left = v;
                else if (b->right == a)
                    b->right = v;
                if (u->left == v)
                    u->left = a;
                else
                    u->right = a;
                }
            else
                {
                /* rotate the nodes from b to u*/
                p = b;
                q = a;
                x = b->length;
                while (p->left->marked == YES || p->right->marked == YES)
                    {                   
                    if (p->left->marked == YES) 
                        {
                        /* rotate p anticlockwise - prepare pointers for move left */
                        p->anc = p->left;  /* the root will be in the direction we are heading */
                        p->left = p->right;
                        p->right = q;
                        }
                    else 
                        {
                        /* rotate p clockwise - prepare pointers for move right */
                        p->anc = p->right;  /* the root will be in the direction we are heading */
                        p->right = p->left;
                        p->left = q;  
                        }
                    /* OK - let's move!; p->anc points in the right direction
                    don't forget to move the branch lengths as well */
                    q = p;
                    p = p->anc;
                    q->length = p->length;
                    q->upDateTi = YES;
                    }
                /* rotations finished, take care of u */
                if (u->left == v)
                    u->left = u->anc;
                else
                    u->right = u->anc;
                u->length = x;
                /* now swap the subtrees of u and b */
                if (a->left == b)
                    a->left = u;
                else
                    a->right = u;
                u->anc = a;
                v->anc = b;
                if (b->left == a)
                    b->left = v;
                else
                    b->right = v;
                }
            }
        }

    /* move around and potentially swap in crown subtree */
    nCrownNodes = 0;
    if (moveInRoot == NO)       
        {
        x = v->length;  /* save v length in case there is a move */
        for (nCrownNodes=0; nCrownNodes==0 || RandomNumber(seed)<extensionProb; nCrownNodes++) 
            {
            if (c->left == NULL || c->isLocked == YES)
                break;  /* can't go further */

            topologyHasChanged = YES;
            
            /* prepare d for move */
            d->anc = c;
            d->length = c->length;
            d->upDateTi = YES;
            d->upDateCl = YES;
            if (d->isLocked == YES)
                {
                c->isLocked = YES;
                d->isLocked = NO;
                c->lockID = d->lockID;
                d->lockID = -1;
                }
            
            /* go left or right with equal probability */
            if (RandomNumber(seed) < 0.5) 
                {
                /* rotate c anticlockwise - prepare pointers for move left */
                c->anc = c->left;  /* the root will be in the direction we are heading */
                c->left = c->right;
                c->right = d;
                }
            else 
                {
                /* rotate c clockwise - prepare pointers for move right */
                c->anc = c->right;  /* the root will be in the direction we are heading */
                c->right = c->left;
                c->left = d;
                }
            /* OK - let's move!; c->anc points in the right direction */
            d = c;
            c = c->anc;
            }

        /* swap the crown subtrees */
        if (nCrownNodes > 0)
            {
            d->anc = u;
            d->length = x;
            if (u->left == v)
                u->left = d;
            else
                u->right = d;

            c->anc = v;
            if (directionLeft == YES)
                v->left = c;
            else
                v->right = c;
            }
        }

    /* modify branch lengths */
    if (nCrownNodes > 0)
        {
        p = c;
        q = d;
        }
    else if (nRootNodes > 0)
        {
        if (directionUp == YES)
            {
            p = v;
            q = a;
            }
        else
            {
            p = v;
            q = u;
            }
        }
    else
        {
        p = v;
        if (RandomNumber(seed) < 0.5)
            {
            if (RandomNumber(seed) < 0.5)
                q = u;
            else
                {
                if (u->left == v)
                    q = u->right;
                else
                    q = u->left;
                }
            }
        else
            {
            if (RandomNumber(seed) < 0.5)
                q = v->left;
            else
                q = v->right;
            }
        }

    if (p != NULL)
        {
        m = p->length;
        x = m * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV)
                x = minV * minV / x;
            else if (x > maxV)
                x = maxV * maxV / x;
            }
        p->length = x;
        p->upDateTi = YES;

        /* update proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / m);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (m - x);
        }

    if (q != NULL && q->anc != NULL)
        {
        m = q->length;
        x = m * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV)
                x = minV * minV / x;
            else if (x > maxV)
                x = maxV * maxV / x;
            }
        q->length = x;
        q->upDateTi = YES;

        /* update proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / m);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (m - x);
        }

    /* set flags for update of cond likes from v and down to root */
    p = v;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }

    if (topologyHasChanged == YES)
        {
        /* set flags for update of cond likes from u and down to root */
        p = u;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        }

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);
    
#   if defined (DEBUG_ExtSS)
    printf ("After:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
    printf ("Proposal ratio: %f\n",exp(*lnProposalRatio));
    printf ("v: %d  u: %d  c: %d  d: %d  a: %d  b: %d\n",v->index, u->index, 
        c->index, d->index, a->index, b->index);
    printf ("No. nodes moved in root subtree: %d\n",nRootNodes);
    printf ("No. nodes moved in crown subtree: %d\n",nCrownNodes);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    printf ("directionUp = %d -- directionLeft = %d\n", directionUp, directionLeft);
    getchar();
#   endif

    return (NO_ERROR);
}


int Move_ExtSSClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change branch lengths and topology (potentially) using SS-type move 
       with extension probability (rather than window, attachment rate or similar). */

    /* This move picks a branch at random. It then moves away from this branch, one 
       one node at a time, with a probability determined by the extensionProb parameter.
       The process stops when a tip is reached or when a move further upwards would break
       the clock assumption. When the extension process stops, the subtrees supported by
       the two chosen branches are swapped. Since 2010-11-01, the move is Metropolized for
       increased efficiency. */
    /* Note: this move is not compatible with fossilized birth-death model with ancestral fossils */
    
    int         i, *nEvents, numFreeOld, numFreeNew;
    MrBFlt      x, oldALength, oldCLength, extensionProb, igrvar, *igrRate,
                *tk02Rate, *brlens, nu, ran, cumulativeProb, forwardProb,
                backwardProb, minV;
    TreeNode    *p, *q, *a, *c;
    Tree        *t;
    Param       *subParm;

    extensionProb = mvp[0]; /* extension probability */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get min and max brlens in relative time and subst units */
    minV = BRLENS_MIN;

    /* calculate the number of free nodes */
    numFreeOld = t->nNodes-2;
    if (t->nConstraints > 1)
        {
        numFreeOld = 0;
        for (i=0; i<t->nNodes-2; i++)
            {
            p = t->allDownPass[i];
            if (p->anc->left == p)
                q = p->anc->right;
            else
                q = p->anc->left;
            if (p->anc->isLocked == NO || q->isLocked == NO)
                numFreeOld++;
            }
        }

    /* pick a branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes -2))];
        if (p->anc->left == p)
            q = p->anc->right;
        else
            q = p->anc->left;
        }
    while ((p->anc->isLocked == YES && q->isLocked == YES) || p->length < TIME_MIN || q->length < TIME_MIN);
    /* choose subtree that can be swapped */

    /* set up pointers for nodes around the picked branch */
    a = p;
    if (p->anc->left == p)
        q = p->anc->right;
    else
        q = p->anc->left;
    if (p->anc->anc->left == p->anc)
        c = p->anc->anc->right;
    else
        c = p->anc->anc->left;

    /* record branch length */
    oldALength = a->length;

    /* reset scratch variables */
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        p->x = -1;
        p->y = NO;
        }

    /* calculate distance from picked node */
    p = a->anc;
    p->x = 0;
    while (p->isLocked == NO && p->anc != NULL)
        {
        p->anc->x = p->x + 1;
        p = p->anc;
        }
    for (i=t->nIntNodes-2; i>=0; i--)
        {
        p = t->intDownPass[i];
        if (p->x < 0 && p->anc->x >= 0 && p != a && p->isLocked == NO)
            p->x = p->anc->x + 1;
        }

    /* mark the free nodes and calculate the total score */
    cumulativeProb = 0.0; 
    for (i=0; i<t->nNodes-2; i++)
        {
        p = t->allDownPass[i];
        if (p != a && p->anc->x > 0 && a->anc->nodeDepth > p->nodeDepth + minV && p->anc->nodeDepth > a->nodeDepth + minV)
            {
            p->y = YES;
            p->d = pow(0.5 * extensionProb, p->anc->x);
            cumulativeProb += p->d;
            }
        else
            p->d = 0.0;
        }

    /* find the target node */
    ran = RandomNumber(seed) * cumulativeProb;
    x = 0.0;
    for (i=0; i<t->nNodes-2; i++)
        {
        p = t->allDownPass[i];
        if (p->y == YES)
            {
            x += p->d;
            if (x > ran)
                break;
            }
        }
    if (i == t->nNodes - 2)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* record first forward prob */
    forwardProb = p->d / cumulativeProb;

    /* record partner swap branch */
    c = p;
    oldCLength = c->length;

    /* calculate second forward prob */

    /* reset scratch variables */
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        p->x = -1;
        p->y = NO;
        }

    /* calculate distance from picked node */
    p = c->anc;
    p->x = 0;
    while (p->isLocked == NO && p->anc != NULL)
        {
        p->anc->x = p->x + 1;
        p = p->anc;
        }
    for (i=t->nIntNodes-1; i>=0; i--)
        {
        p = t->intDownPass[i];
        if (p->x < 0 && p != c && p->anc->x >= 0 && p->isLocked == NO)
            p->x = p->anc->x + 1;
        }

    /* mark the free nodes and calculate the total score */
    cumulativeProb = 0.0; 
    for (i=0; i<t->nNodes-2; i++)
        {
        p = t->allDownPass[i];
        if (p != c && p->anc->x > 0 && c->anc->nodeDepth > p->nodeDepth + minV && p->anc->nodeDepth > c->nodeDepth + minV)
            {
            p->y = YES;
            p->d = pow(0.5 * extensionProb, p->anc->x);
            cumulativeProb += p->d;
            }
        else
            p->d = 0.0;
        }

    /* now we can calculate second forward prob */
    forwardProb += a->d / cumulativeProb;

    /* swap subtrees */
    if (a->anc->left == a)
        a->anc->left = c;
    else
        a->anc->right = c;
    if (c->anc->left == c)
        c->anc->left = a;
    else
        c->anc->right = a;
    p = a->anc;
    a->anc = c->anc;
    c->anc = p;
    a->length = a->anc->nodeDepth - a->nodeDepth;
    c->length = c->anc->nodeDepth - c->nodeDepth;

    /* get down pass sequence */
    GetDownPass (t);

    /* set tiprobs update flags */
    a->upDateTi = YES;
    c->upDateTi = YES;

    /* set flags for update of cond likes from a->anc and down to root */
    p = a->anc;
    while (p->anc != NULL)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }

    /* set flags for update of cond likes from c->anc and down to root */
    p = c->anc;
    while (p->anc != NULL)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }

    /* adjust prior ratio for clock tree */
    if (LogClockTreePriorRatio(param, chain, &x) == ERROR)
        return (ERROR);
    (*lnPriorRatio) += x;

    /* calculate first backward prob */

    /* reset scratch variables */
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        p->x = -1;
        p->y = NO;
        }

    /* calculate distance from picked node */
    p = a->anc;
    p->x = 0;
    while (p->isLocked == NO && p->anc != NULL)
        {
        p->anc->x = p->x + 1;
        p = p->anc;
        }
    for (i=t->nIntNodes-1; i>=0; i--)
        {
        p = t->intDownPass[i];
        if (p->x < 0 && p != a && p->anc->x >= 0 && p->isLocked == NO)
            p->x = p->anc->x + 1;
        }

    /* mark the free nodes and calculate the total score */
    cumulativeProb = 0.0; 
    for (i=0; i<t->nNodes-2; i++)
        {
        p = t->allDownPass[i];
        if (p != a && p->anc->x > 0 && a->anc->nodeDepth > p->nodeDepth + minV && p->anc->nodeDepth > a->nodeDepth + minV)
            {
            p->y = YES;
            p->d = pow(0.5 * extensionProb, p->anc->x);
            cumulativeProb += p->d;
            }
        else
            p->d = 0.0;
        }

    /* calculate first backward prob */
    backwardProb = c->d / cumulativeProb;

    /* calculate second backward prob */

    /* reset scratch variables */
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        p->x = -1;
        p->y = NO;
        }

    /* calculate distance from picked node */
    p = c->anc;
    p->x = 0;
    while (p->isLocked == NO && p->anc != NULL)
        {
        p->anc->x = p->x + 1;
        p = p->anc;
        }
    for (i=t->nIntNodes-1; i>=0; i--)
        {
        p = t->intDownPass[i];
        if (p->x < 0 && p != c && p->anc->x >= 0 && p->isLocked == NO)
            p->x = p->anc->x + 1;
        }

    /* mark the free nodes and calculate the total score */
    cumulativeProb = 0.0; 
    for (i=0; i<t->nNodes-2; i++)
        {
        p = t->allDownPass[i];
        if (p != c && p->anc->x > 0 && c->anc->nodeDepth > p->nodeDepth + minV && p->anc->nodeDepth > c->nodeDepth + minV)
            {
            p->y = YES;
            p->d = pow(0.5 * extensionProb, p->anc->x);
            cumulativeProb += p->d;
            }
        else
            p->d = 0.0;
        }

    /* calculate second backward prob */
    backwardProb += a->d / cumulativeProb;

    /* now we can calculate proposal ratio */
    (*lnProposalRatio) += log (backwardProb / forwardProb);

    /* adjust for number of free nodes */
    numFreeNew = t->nNodes-2;
    if (t->nConstraints > 1)
        {
        numFreeNew = 0;
        for (i=0; i<t->nNodes-2; i++)
            {
            p = t->allDownPass[i];
            if (p->anc->left == p)
                q = p->anc->right;
            else
                q = p->anc->left;
            if (p->anc->isLocked == NO || q->isLocked == NO)
                numFreeNew++;
            }
        (*lnProposalRatio) += log(numFreeOld / numFreeNew);
        }

    /* adjust proposal and prior ratio for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
        {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            nEvents = subParm->nEvents[2*chain+state[chain]];

            /* proposal ratio */
            (*lnProposalRatio) += nEvents[a->index] * log (a->length / oldALength);
            (*lnProposalRatio) += nEvents[c->index] * log (c->length / oldCLength);

            /* prior ratio: no effect because tree length is the same */

            /* update effective evolutionary lengths */
            if (UpdateCppEvolLengths (subParm, a, chain) == ERROR || UpdateCppEvolLengths (subParm, c, chain) == ERROR)
                {
                abortMove = YES;
                return (NO_ERROR);
                }
            }
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);

            /* no proposal ratio effect */

            /* prior ratio and update of effective evolutionary lengths */
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[c->anc->index], nu*oldALength, tk02Rate[a->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[a->anc->index], nu* a->length, tk02Rate[a->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[a->anc->index], nu*oldCLength, tk02Rate[c->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[c->anc->index], nu* c->length, tk02Rate[c->index]);
            brlens[a->index] = a->length * (tk02Rate[a->index] + tk02Rate[a->anc->index])/2.0;
            brlens[c->index] = c->length * (tk02Rate[c->index] + tk02Rate[c->anc->index])/2.0;
            }
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            /* get relevant parameters */
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);

            /* prior ratio and update of effective evolutionary lengths */
            (*lnPriorRatio) -= LnProbGamma (oldALength/igrvar, oldALength/igrvar, igrRate[a->index]);
            (*lnPriorRatio) -= LnProbGamma (oldCLength/igrvar, oldCLength/igrvar, igrRate[c->index]);
            (*lnPriorRatio) += LnProbGamma (a->length /igrvar, a->length /igrvar, igrRate[a->index]);
            (*lnPriorRatio) += LnProbGamma (c->length /igrvar, c->length /igrvar, igrRate[c->index]);
            brlens[a->index] = igrRate[a->index] * a->length;
            brlens[c->index] = igrRate[c->index] * c->length;
            }
        }
    
    assert (*lnPriorRatio == *lnPriorRatio);

    return (NO_ERROR);
}


int Move_ExtTBR (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change topology (and branch lengths) using TBR with extension probability. */

    /* This move picks an internal branch and two "danglies", modifies their lengths
       independently according to the method of Larget & Simon (1999: MBE); it then
       moves the danglies away from their original position one node at a time with
       a probability determined by the extensionProb parameter. When the danglies are moved,
       their direction is changed -- "reflection" necessary to enable the back move.

       This move type has been tested on all combinations of rooted and unrooted,
       constrained and unconstrained trees */
    
    int         i, j, topologyHasChanged, nCrownNodes, nRootNodes, directionLeft, directionUp, 
                isVPriorExp, alwaysMoveRoot, isCrownStartConstrained, isRootStartConstrained, isStopConstrained;
    MrBFlt      m, x, y, tuning, maxV, minV, extensionProb, brlensExp=0.0;
    TreeNode    *p, *a, *b, *c, *d, *u, *v;
    Tree        *t;
    ModelParams *mp;

    /* these parameters should be possible to set by user */
    extensionProb = mvp[0]; /* extension probability */
    tuning = mvp[1];        /* Larget & Simon's tuning parameter lambda */
    
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1] < BRLENS_MAX ? mp->brlensUni[1] : BRLENS_MAX;
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);

    topologyHasChanged = NO;

#   if defined (DEBUG_ExtTBR)
    printf ("Before:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
#   endif
    
    /* pick an internal branch */
    do  {
        p = t->intDownPass[(int)(RandomNumber(seed) * (t->nIntNodes-1))];
        if (p->anc->left == p)
            a = p->anc->right;
        else
            a = p->anc->left;
        i = j = 0;
        if (a->isLocked == YES || a->left == NULL)
            i++;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL)
            i++;
        if (p->left->isLocked == YES || p->left->left == NULL)
            j++;
        if (p->right->isLocked == YES || p->right->left == NULL)
            j++;
        } while (i == 2 && j == 2);

    /* determine whether to move first step unconditionally in root or in crown */
    if (j == 2)
        alwaysMoveRoot = YES;
    else if (i == 2)
        alwaysMoveRoot = NO;
    else if (RandomNumber(seed) < 0.5)
        alwaysMoveRoot = YES;
    else
        alwaysMoveRoot = NO;

    /* determine any starting constraints */
    isCrownStartConstrained = isRootStartConstrained = NO;
    if (i >= 1)
        isRootStartConstrained = YES;
    if (j >= 1)
        isCrownStartConstrained = YES;

    /* set up pointers for nodes around the picked branch */
    /* cut the tree into crown, root and attachment part */
    /* change the relevant lengths in the attachment part */
    /* the lengths of a and v are automatically contained in the */
    /* "attachment" part but the length of c has to be stored in x */
    v = p;
    u = p->anc;

    /* set up pointers for crown part */
    /* also determine direction of move in crown part */
    if (v->right->left == NULL || v->right->isLocked == YES)
        directionLeft = YES;
    else if (v->left->left == NULL || v->left->isLocked == YES)
        directionLeft = NO;
    else if (RandomNumber(seed) < 0.5)
        directionLeft = YES;
    else
        directionLeft = NO;
    if (directionLeft == YES)
        {
        c = v->left;
        d = v->right;
        }
    else
        {
        c = v->right;
        d = v->left;
        }

    /* cut and reconnect crown part */
    c->anc = d;
    d->anc = c;
    
    /* record c length and adjust with multiplier using reflection */
    m = c->length;
    x = c->length * exp(tuning * (RandomNumber(seed) - 0.5));       /* save the modified dangling branch for later use */
    while (x < minV || x > maxV)
        {
        if (x < minV) x = minV * minV / x;
        if (x > maxV) x = maxV * maxV / x;
        }
    
    /* calculate proposal and prior ratio based on length modification */
    (*lnProposalRatio) = log (x / m);
    if (isVPriorExp == YES)
        (*lnPriorRatio) = brlensExp * (m - x);

    /* record v length and adjust with multiplier using reflection*/
    m = v->length;
    v->length *= exp(tuning * (RandomNumber(seed) - 0.5));
    while (v->length < minV || v->length > maxV)
        {
        if (v->length < minV)
            v->length = minV * minV / v->length;
        else if (v->length > maxV)
            v->length = maxV * maxV / v->length;
        }
    v->upDateTi = YES;

    /* adjust proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (v->length / m);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (m - v->length);

    /* mark nodes in root part */
    /* also determine direction of move in root part */
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    if (u->anc->anc == NULL || u->isLocked == YES)
        directionUp = YES;
    else if (a->left == NULL || a->isLocked == YES)
        directionUp = NO;
    else if (RandomNumber(seed) < 0.5)
        directionUp = YES;
    else
        directionUp = NO;
    if (directionUp == NO)
        {
        /* switch a and b */
        b = a;
        a = u->anc;
        }

    /* cut root part */
    if (directionUp == NO)
        {
        b->anc = a;
        if (a->left == u)
            a->left = b;
        else
            a->right = b;
        }
    else 
        {
        a->anc = b;
        if (b->left == u)
            b->left = a;
        else
            b->right = a;
        y = a->length;
        a->length = u->length;
        u->length = y;
        a->upDateTi = YES;
        }

    /* adjust length of branch to be modified */
    /* if it is not the root branch of a rooted tree */
    if (t->isRooted == NO || u->anc->anc != NULL) 
        {
        m = u->length;
        u->length *= exp(tuning * (RandomNumber(seed) - 0.5));
        while (u->length < minV || u->length > maxV)
            {
            if (u->length < minV)
                u->length = minV * minV / u->length;
            else if (u->length > maxV)
                u->length = maxV * maxV / u->length;
            }

        /* adjust proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (u->length / m);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (m - u->length);
        }
    u->upDateTi = YES;

    /* move around in root subtree */
    for (nRootNodes=0; (alwaysMoveRoot == YES && nRootNodes == 0) || RandomNumber(seed) < extensionProb; nRootNodes++)
        {
        if (directionUp == YES) 
            {   /* going up tree */
            if (a->left == NULL || a->isLocked == YES)
                break;      /* can't go further */
            topologyHasChanged = YES;
            b = a;
            if (RandomNumber(seed) < 0.5)
                a = a->left;
            else
                a = a->right;
            if (u->isLocked == YES)
                {
                b->isLocked = YES;
                u->isLocked = NO;
                b->lockID = u->lockID;
                u->lockID = 0;
                }
            }
        else 
            {   /* going down tree */
            if (a->anc == NULL || u->isLocked == YES)
                break;      /* can't go further */
            topologyHasChanged = YES;
            if (RandomNumber(seed)<0.5) 
                {
                directionUp = YES; /* switch direction */
                /* find sister of a */
                if (a->left == b) 
                    {
                    b = a;
                    a = a->right;
                    }
                else 
                    {  
                    b = a;
                    a = a->left;
                    }
                /* as long as we are moving upwards
                the cond likes to update will be
                flagged by the last pass from u to the root */
                }   
            else 
                {   /* continue down */
                b = a;
                a = a->anc;
                b->upDateCl = YES; 
                if (b->isLocked == YES)
                    {
                    u->isLocked = YES;
                    b->isLocked = NO;
                    u->lockID = b->lockID;
                    b->lockID = 0;
                    }
                }
            }
        }

    /* adjust proposal ratio for root move if unbalanced */
    isStopConstrained = NO;
    if (directionUp == YES && (a->left == NULL || a->isLocked == YES))
        isStopConstrained = YES;
    if (directionUp == NO && (a->anc  == NULL || u->isLocked == YES))
        isStopConstrained = YES;
    if (nRootNodes > 0)
        {
        if (isRootStartConstrained == YES && isStopConstrained == NO)
            (*lnProposalRatio) -= log (2.0 * (1.0 - extensionProb));
        else if (isRootStartConstrained == NO && isStopConstrained == YES)
            (*lnProposalRatio) += log (2.0 * (1.0 - extensionProb));
        }

    /* move around in crown subtree */
    for (nCrownNodes=0; (alwaysMoveRoot == NO && nCrownNodes == 0) || RandomNumber(seed) < extensionProb; nCrownNodes++)
        {
        if (c->left == NULL || c->isLocked == YES)
            break;  /* can't go further */
        topologyHasChanged = YES;
        if (RandomNumber(seed) < 0.5) 
            {
            /* rotate c anticlockwise - prepare pointers for move left */
            c->anc = c->left;  /* the root will be in the direction we are heading */
            c->left = c->right;
            c->right = d;
            }
        else 
            {
            /* rotate c clockwise - prepare pointers for move right */
            c->anc = c->right;  /* the root will be in the direction we are heading */
            c->right = c->left;
            c->left = d;  
            }
        /* OK - let's move!; c->anc points in the right direction
        don't forget to move the branch lengths as well */
        d = c;
        c = c->anc;
        d->length = c->length;
        d->upDateCl = YES; 
        d->upDateTi = YES;
        }

    /* adjust proposal ratio for crown move if unbalanced */
    isStopConstrained = NO;
    if (c->left == NULL || c->isLocked == YES)
        isStopConstrained = YES;
    if (nCrownNodes > 0)
        {
        if (isCrownStartConstrained == YES && isStopConstrained == NO)
            (*lnProposalRatio) -= log (2.0 * (1.0 - extensionProb));
        else if (isCrownStartConstrained == NO && isStopConstrained == YES)
            (*lnProposalRatio) += log (2.0 * (1.0 - extensionProb));
        }

    /* combine the subtrees */
    c->anc = v;
    d->anc = v;
    if (directionLeft == YES) 
        {
        v->left = c;
        v->right = d;
        }
    else 
        {
        v->left = d;
        v->right = c;
        }

    /* the dangling branch is inserted in reverted position
       such that the back move will be possible
       if we have moved around in crown subtree
       otherwise it is left in its original position */
    if (nCrownNodes > 0)
        {
        d->length = x;
        d->upDateTi = YES;
        }
    else
        {
        c->length = x;
        c->upDateTi = YES;
        }

    if (directionUp == YES) 
        {
        u->anc = b;
        if (u->left == v)
            u->right = a;
        else 
            u->left = a;
        a->anc = u;
        if (b->left == a)
            b->left = u;
        else
            b->right = u;
        /* the dangling branch is contained in u->length
           and will automatically be inserted in the right position
           to enable the back move regardless of whether it was
           initially directed upwards or downwards
           BUT if we haven't moved in root subtree, it is advantageous (necessary
           for rooted trees) to avoid switching branches, which occurs otherwise
           if directionUp == YES */
        if (nRootNodes == 0) 
            {
            y = u->length;
            u->length = a->length;
            a->length = y;
            a->upDateTi = YES;
            u->upDateTi = NO;   /* u retains its old length */
            }
        }
    else 
        {
        u->anc = a;
        if (u->left == v)
            u->right = b;
        else
            u->left = b;
        b->anc = u;
        if (a->left == b)
            a->left = u;
        else
            a->right = u;
        /* the modified branch contained in u->length will have
           to be moved to b->length to enable back move
           BUT if we haven't moved, it is better to keep it in place
           (necessary for rooted trees) */
        if (nRootNodes > 0) 
            {
            y = u->length;
            u->length = b->length;
            b->length = y;
            b->upDateTi = YES;
            }
        }
        
    /* set flags for update of cond likes from v and down to root */
    p = v;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }
    
    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);

#   if defined (DEBUG_ExtTBR)
    printf ("After:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
    printf ("Proposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  c: %d  d: %d  a: %d  b: %d\n",v->index, u->index, 
            c->index, d->index, a->index, b->index);
    printf ("No. nodes moved in root subtree: %d\n",nRootNodes);
    printf ("No. nodes moved in crown subtree: %d\n",nCrownNodes);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    getchar();
#   endif

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_GeneRate_Dir: Change gene rate multiplier using Dirichlet
|      proposal.
|
----------------------------------------------------------------*/
int Move_GeneRate_Dir (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, nRates, isValid;
    MrBFlt      alphaPi, *value, *subValue, numSites, *alphaDir, x, y, sum,
                rate_pot, *dirParm, *oldRate, *newRate;

    /* allocate memory */
    dirParm = (MrBFlt *) SafeCalloc (3*(numTopologies-1), sizeof(MrBFlt));
    oldRate = dirParm + numCurrentDivisions;
    newRate = dirParm + 2*numCurrentDivisions;

    /* get number of rates */
    nRates = param->nValues;

    /* get pointer to rates and number of uncompressed chars */
    value = GetParamVals(param, chain, state[chain]);
    subValue = GetParamSubVals(param, chain, state[chain]);

    /* get Dirichlet parameters */
    alphaDir = subValue + nRates;

    /* calculate old ratesum proportions */
    numSites = 0.0;
    for (i=0; i<nRates; i++)
        numSites += subValue[i];  /* numSites should be equal to the number of sites */
    for (i=0; i<nRates; i++)
        oldRate[i] = value[i] * subValue[i] / numSites;
    
    /* get so called alphaPi parameter */
    alphaPi = mvp[0] * nRates;
    
    /* multiply old ratesum proportions with some large number to get new values close to the old ones */
    for (i=0; i<nRates; i++)
        dirParm[i] = oldRate[i] * alphaPi;
    
    /* get new values */
    DirichletRandomVariable (dirParm, newRate, nRates, seed);

    /* check new values. we rely on newRate be already normalized  */
    while (1)
        {
        sum = 0.0;
        rate_pot = 1.0;
        isValid=1;
        for (i=0; i<nRates; i++)
            {
            if (newRate[i] <= DIR_MIN)
                {
                if (newRate[i] < DIR_MIN)
                    {
                    newRate[i] = DIR_MIN;
                    isValid=0;
                    }
                rate_pot -= DIR_MIN;
                }
            else
                sum += newRate[i];
            }
        if (isValid==1) break;
        for (i=0; i<nRates; i++)
            {
            if (newRate[i]!=DIR_MIN)
                newRate[i] = rate_pot * newRate[i] / sum;
            }
        }

    /* calculate and copy new rate ratio values back */
    for (i=0; i<nRates; i++)
        value[i] = newRate[i] * (numSites / subValue[i]);
    
    /* get proposal ratio */
    sum = 0.0;
    for (i=0; i<nRates; i++)
        sum += newRate[i]*alphaPi;
    x = LnGamma(sum);
    for (i=0; i<nRates; i++)
        x -= LnGamma(newRate[i]*alphaPi);
    for (i=0; i<nRates; i++)
        x += (newRate[i]*alphaPi-1.0)*log(oldRate[i]);
    sum = 0.0;
    for (i=0; i<nRates; i++)
        sum += oldRate[i]*alphaPi;
    y = LnGamma(sum);
    for (i=0; i<nRates; i++)
        y -= LnGamma(oldRate[i]*alphaPi);
    for (i=0; i<nRates; i++)
        y += (oldRate[i]*alphaPi-1.0)*log(newRate[i]);
    (*lnProposalRatio) = x - y;

    /* get prior ratio */
    x = y = 0.0;
    for (i=0; i<nRates; i++)
        x += (alphaDir[i]-1.0)*log(newRate[i]);
    for (i=0; i<nRates; i++)
        y += (alphaDir[i]-1.0)*log(oldRate[i]);
    (*lnPriorRatio) = x - y;

    /* Set update flags for all partitions that share the rate multiplier. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* may need to hit update flag for cijks when you have a covarion model */
    for (i=0; i<param->nRelParts; i++)
        if (modelSettings[param->relParts[i]].nCijkParts > 1)
            modelSettings[param->relParts[i]].upDateCijk = YES;

    free (dirParm);

    return (NO_ERROR);
}


int Move_RateShape_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change gamma/lnorm shape parameter using multiplier */
    
    int         i, isAPriorExp, isValidA;
    MrBFlt      oldA, newA, minA, maxA, alphaExp=0.0, ran, factor, tuning, *rates;
    ModelParams *mp;

    /* get tuning parameter */
    tuning = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* get minimum and maximum values for alpha */
    if (param->paramId == SHAPE_UNI)
        {
        minA = mp->shapeUni[0];
        maxA = mp->shapeUni[1];
        if (minA < MIN_SHAPE_PARAM)
            minA = MIN_SHAPE_PARAM;
        if (maxA > MAX_SHAPE_PARAM)
            maxA = MAX_SHAPE_PARAM;
        isAPriorExp = NO;
        }
    else
        {
        minA = MIN_SHAPE_PARAM;
        maxA = MAX_SHAPE_PARAM;
        alphaExp = mp->shapeExp;
        isAPriorExp = YES;
        }

    /* get old value of alpha */
    oldA = *GetParamVals(param, chain, state[chain]);

    /* change value for alpha */
    ran = RandomNumber(seed);
    factor = exp(tuning * (ran - 0.5));
    newA = oldA * factor;

    /* check validity */
    isValidA = NO;
    do  {
        if (newA < minA)
            newA = minA * minA / newA;
        else if (newA > maxA)
            newA = maxA * maxA / newA;
        else
            isValidA = YES;
        } while (isValidA == NO);

    /* get proposal ratio */
    *lnProposalRatio = log(newA / oldA);
    
    /* get prior ratio */
    if (isAPriorExp == NO)
        *lnPriorRatio = 0.0;
    else
        *lnPriorRatio = -alphaExp * (newA - oldA);
    
    /* copy new alpha value back */
    *GetParamVals(param, chain, state[chain]) = newA;
    
    /* now, update rate category information */
    rates = GetParamSubVals (param, chain, state[chain]);
    
    if(!strcmp(mp->ratesModel, "LNorm"))
        {
        if (DiscreteLogNormal (rates, newA, mp->numGammaCats, 1) == ERROR)
            return (ERROR);
        }
    else  /* gamma rate */
        {
        if (DiscreteGamma (rates, newA, newA, mp->numGammaCats, 0) == ERROR)
            return (ERROR);
        }

    /* Set update flags for all partitions that share this alpha. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* We need to update flags when we have a covarion model */
    for (i=0; i<param->nRelParts; i++)
        if (modelSettings[param->relParts[i]].nCijkParts > 1)
            modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


int Move_Growth_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    MrBFlt          oldG, newG, lambda, minG, maxG, ran, oldLnPrior, newLnPrior, curTheta;
    ModelParams     *mp;
    ModelInfo       *m;
    Tree            *t;

    /* get tuning parameter */
    lambda = mvp[0];

    /* get model params */
    m = &modelSettings[param->relParts[0]];
    mp = &modelParams[param->relParts[0]];
    curTheta = (*GetParamVals(m->popSize, chain, state[chain])) * (*GetParamVals(m->clockRate, chain, state[chain]));
    if (!strcmp(mp->ploidy, "Diploid"))
        curTheta *= 4.0;
    else if (!strcmp(mp->ploidy, "Zlinked"))
        curTheta *= 3.0;
    else
        curTheta *= 2.0;
    
    /* get minimum and maximum values for growth */
    minG = param->min;
    maxG = param->max;
    
    /* get old value of theta */
    newG = oldG = *GetParamVals(param, chain, state[chain]);
    if (newG < minG)
        newG = oldG = minG;

    /* change value of growth */
    ran = RandomNumber(seed);
    newG = oldG * exp (lambda * (ran - 0.5));
    
    /* check that new value is valid */
    while (newG < minG || newG > maxG)
        {
        if (newG < minG)
            newG = minG * minG / newG;
        else if (newG > maxG)
            newG = maxG * maxG / newG;
        }
    
    /* get proposal ratio */
    (*lnProposalRatio) = log (newG / oldG);
    
    /* get prior ratio */
    t = GetTree(modelSettings[param->relParts[0]].brlens,chain,state[chain]);
    if (LnCoalescencePriorPr (t, &oldLnPrior, curTheta, oldG) == ERROR)
        {
        MrBayesPrint ("%s   Problem calculating prior for coalescent process\n", spacer);
        return (ERROR);
        }
    if (LnCoalescencePriorPr (t, &newLnPrior, curTheta, newG) == ERROR)
        {
        MrBayesPrint ("%s   Problem calculating prior for coalescent process\n", spacer);
        return (ERROR);
        }
    (*lnPriorRatio) = newLnPrior - oldLnPrior + param->LnPriorRatio(newG, oldG, param->priorParams);

    /* copy new growth value back */
    *GetParamVals(param, chain, state[chain]) = newG;

    return (NO_ERROR);
}


int Move_IgrBranchRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move one IGR relaxed clock branch rate using multiplier */

    int         i;
    MrBFlt      newRate, oldRate, tuning, minR, maxR, igrvar, *igrRate, *brlens;
    TreeNode    *p = NULL;
    ModelInfo   *m;
    Tree        *t;
    TreeNode    *q;

    /* get the tuning parameter */
    tuning = mvp[0];
    
    /* get the model settings */
    m = &modelSettings[param->relParts[0]];

    /* get the IGR branch rate and effective branch length data */
    igrRate = GetParamVals (param, chain, state[chain]);
    brlens = GetParamSubVals (param, chain, state[chain]);

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get minimum and maximum rate */
    minR = RATE_MIN;
    maxR = RATE_MAX;
    
    /* randomly pick a branch */
    do  {
        i = (int) (RandomNumber(seed) * (t->nNodes -2));
        p = t->allDownPass[i];
        }
    while (p->length < TIME_MIN);  // not ancestral fossil
    
    /* find new rate using multiplier */
    oldRate = igrRate[p->index];
    newRate = oldRate * exp ((0.5 - RandomNumber(seed)) * tuning);
    
    /* reflect if necessary */
    while (newRate < minR || newRate > maxR)
        {
        if (newRate < minR)
            newRate = minR * minR / newRate;
        if (newRate > maxR)
            newRate = maxR * maxR / newRate;
        }
    
    igrRate[p->index] = newRate;

    /* calculate prior ratio */
    igrvar = *GetParamVals (m->igrvar, chain, state[chain]);
    (*lnPriorRatio) = LnProbGamma (p->length/igrvar, p->length/igrvar, newRate)
                    - LnProbGamma (p->length/igrvar, p->length/igrvar, oldRate);

    /* calculate proposal ratio */
    (*lnProposalRatio) = log (newRate / oldRate);
    
    /* update branch evolution lengths */
    brlens[p->index] = newRate * p->length;

    /* set update of transition probability */
    p->upDateTi = YES;

    /* set update of cond likes down to root */
    q = p->anc;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }

    return (NO_ERROR);
}


int Move_IgrBranchRate2 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move one IGR relaxed clock branch rate using sliding window */
    
    int         i;
    MrBFlt      newRate, oldRate, window, minR, maxR, igrvar, *igrRate, *brlens;
    TreeNode    *p = NULL;
    ModelInfo   *m;
    Tree        *t;
    TreeNode    *q;
    
    /* get the tuning parameter */
    window = mvp[0];
    
    /* get the model settings */
    m = &modelSettings[param->relParts[0]];
    
    /* get the IGR branch rate and effective branch length data */
    igrRate = GetParamVals (param, chain, state[chain]);
    brlens = GetParamSubVals (param, chain, state[chain]);
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* get minimum and maximum rate */
    minR = RATE_MIN;
    maxR = RATE_MAX;
    
    /* randomly pick a branch */
    do  {
        i = (int) (RandomNumber(seed) * (t->nNodes -2));
        p = t->allDownPass[i];
        }
    while (p->length < TIME_MIN);  // not ancestral fossil
    
    /* find new rate using multiplier */
    oldRate = igrRate[p->index];
    newRate = oldRate + window * (RandomNumber(seed) - 0.5);
    
    /* reflect if necessary */
    while (newRate < minR || newRate > maxR)
        {
        if (newRate < minR)
            newRate = 2 * minR - newRate;
        if (newRate > maxR)
            newRate = 2 * maxR - newRate;
        }
    
    igrRate[p->index] = newRate;
    
    /* calculate prior ratio */
    igrvar = *GetParamVals (m->igrvar, chain, state[chain]);
    (*lnPriorRatio) = LnProbGamma (p->length/igrvar, p->length/igrvar, newRate)
                    - LnProbGamma (p->length/igrvar, p->length/igrvar, oldRate);
    
    /* calculate proposal ratio */
    (*lnProposalRatio) = 0.0;
    
    /* update branch evolution lengths */
    brlens[p->index] = newRate * p->length;
    
    /* set update of transition probability */
    p->upDateTi = YES;
    
    /* set update of cond likes down to root */
    q = p->anc;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }
    
    return (NO_ERROR);
}


int Move_IgrVar (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move the variance of the IGR relaxed clock model using multiplier */

    int         i, j;
    MrBFlt      oldIgrvar, newIgrvar, minIgrvar, maxIgrvar, tuning, *igrRate;
    Model       *mp;
    TreeNode    *p;
    Tree        *t;

    /* get tuning parameter */
    tuning = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get the min and max values */
    minIgrvar = IGRVAR_MIN;
    maxIgrvar = IGRVAR_MAX;
    if (!strcmp(mp->igrvarPr,"Uniform"))
        {
        minIgrvar = (mp->igrvarUni[0] < IGRVAR_MIN) ? IGRVAR_MIN : mp->igrvarUni[0];
        maxIgrvar = (mp->igrvarUni[1] > IGRVAR_MAX) ? IGRVAR_MAX : mp->igrvarUni[1];
        }
    
    /* get the igr variance */
    oldIgrvar = *GetParamVals (param, chain, state[chain]);

    /* set new value */
    newIgrvar = oldIgrvar * exp ((0.5 - RandomNumber(seed))*tuning);
    
    /* reflect if necessary */
    while (newIgrvar < minIgrvar || newIgrvar > maxIgrvar)
        {
        if (newIgrvar < minIgrvar)
            newIgrvar = minIgrvar * minIgrvar / newIgrvar;
        if (newIgrvar > maxIgrvar)
            newIgrvar = maxIgrvar * maxIgrvar / newIgrvar;
        }
    
    /* store new value */
    (*GetParamVals (param, chain, state[chain])) = newIgrvar;

    /* calculate prior ratio */
    for (i=0; i<param->nSubParams; i++)
        {
        igrRate = GetParamVals (param->subParams[i], chain, state[chain]);
        t = GetTree (param->subParams[i], chain, state[chain]);
        for (j=0; j<t->nNodes-2; j++)
            {
            p = t->allDownPass[j];
            if (p->length > 0.0)  // not ancestral fossil
                {
                (*lnPriorRatio) -= LnProbGamma (p->length/oldIgrvar, p->length/oldIgrvar, igrRate[p->index]);
                (*lnPriorRatio) += LnProbGamma (p->length/newIgrvar, p->length/newIgrvar, igrRate[p->index]);
                }
            }
        }

    /* take prior on Igrvar into account */
    if (!strcmp(mp->igrvarPr,"Exponential"))
        (*lnPriorRatio) += mp->igrvarExp * (oldIgrvar - newIgrvar);
    
    /* calculate proposal ratio */
    (*lnProposalRatio) = log (newIgrvar / oldIgrvar);

    /* we do not need to update likelihoods */
    for (i=0; i<param->nRelParts; i++)
        {
        modelSettings[param->relParts[i]].upDateCl = NO;
        }

    return (NO_ERROR);
}


int Move_MixedBranchRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move one relaxed clock branch rate using multiplier */
    
    int         i, *rclModel=NULL;
    MrBFlt      newRate, oldRate, tuning, minR, maxR, mxvar, *mxRate, *brlens;
    TreeNode    *p = NULL, *q;
    ModelInfo   *m;
    Tree        *t;
    
    /* get the tuning parameter */
    tuning = mvp[0];
    
    /* get the model settings */
    m = &modelSettings[param->relParts[0]];
    
    /* get the branch rate and effective branch length data */
    mxRate = GetParamVals (param, chain, state[chain]);
    brlens = GetParamSubVals (param, chain, state[chain]);
    rclModel = GetParamIntVals(param, chain, state[chain]);

    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* get minimum and maximum rate */
    minR = RATE_MIN;
    maxR = RATE_MAX;
    
    /* randomly pick a length */
    do  {
        i = (int) (RandomNumber(seed) * (t->nNodes -2));
        p = t->allDownPass[i];
        }
    while (p->length < TIME_MIN);  // not ancestral fossil
    
    /* find new rate using multiplier */
    oldRate = mxRate[p->index];
    newRate = oldRate * exp ((0.5 - RandomNumber(seed)) * tuning);
    
    /* reflect if necessary */
    while (newRate < minR || newRate > maxR)
        {
        if (newRate < minR)
            newRate = minR * minR / newRate;
        if (newRate > maxR)
            newRate = maxR * maxR / newRate;
        }
    
    mxRate[p->index] = newRate;

    /* calculate prior ratio */
    mxvar = *GetParamVals (m->mixedvar, chain, state[chain]);

    if (*rclModel == RCL_TK02)
        {
        (*lnPriorRatio) += LnRatioTK02LogNormal (mxRate[p->anc->index], mxvar*p->length, newRate, oldRate);
        if (p->left != NULL)
            {
            if (p->left->length > 0.0)
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (oldRate, mxvar*p->left->length,  mxRate[p->left->index ]);
                (*lnPriorRatio) += LnProbTK02LogNormal (newRate, mxvar*p->left->length,  mxRate[p->left->index ]);
                }
            if (p->right->length > 0.0)
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (oldRate, mxvar*p->right->length, mxRate[p->right->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (newRate, mxvar*p->right->length, mxRate[p->right->index]);
                }
            }
        
        /* update branch evolution lengths */
        brlens[p->index] = p->length * (newRate + mxRate[p->anc->index]) / 2.0;
        if (p->left != NULL)
            {
            brlens[p->left->index ] = p->left->length  * (mxRate[p->left->index ] + newRate) / 2.0;
            brlens[p->right->index] = p->right->length * (mxRate[p->right->index] + newRate) / 2.0;
            }
        
        /* set update of ti probs */
        p->upDateTi = YES;
        if (p->left != NULL)
            {
            p->left ->upDateTi = YES;
            p->right->upDateTi = YES;
            }
        }
    else if (*rclModel == RCL_IGR)
        {
        (*lnPriorRatio) -= LnProbGamma (p->length/mxvar, p->length/mxvar, oldRate);
        (*lnPriorRatio) += LnProbGamma (p->length/mxvar, p->length/mxvar, newRate);
        
        brlens[p->index] = newRate * p->length;

        /* set update of transition probability */
        p->upDateTi = YES;
        }
    
    /* calculate proposal ratio */
    (*lnProposalRatio) = log (newRate / oldRate);
    
    /* set update of cond likes down to root */
    p->upDateCl = YES;
    q = p->anc;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }
    
    return (NO_ERROR);
}


int Move_MixedVar (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move the variance of the mixed relaxed clock models using multiplier */

    int         i, j, *rclModel=NULL;
    MrBFlt      oldVar, newVar, minVar, maxVar, tuning, *igrRate, *tk02Rate;
    Model       *mp;
    TreeNode    *p;
    Tree        *t;
    
    /* get tuning parameter */
    tuning = mvp[0];
    
    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get the min and max values */
    minVar = MIXEDVAR_MIN;
    maxVar = MIXEDVAR_MAX;
    if (!strcmp(mp->mixedvarPr,"Uniform"))
        {
        minVar = (mp->mixedvarUni[0] < MIXEDVAR_MIN) ? MIXEDVAR_MIN : mp->mixedvarUni[0];
        maxVar = (mp->mixedvarUni[1] > MIXEDVAR_MAX) ? MIXEDVAR_MAX : mp->mixedvarUni[1];
        }
    
    /* get the variance */
    oldVar = *GetParamVals (param, chain, state[chain]);
    
    /* set new value */
    newVar = oldVar * exp ((0.5 - RandomNumber(seed))*tuning);
    
    /* reflect if necessary */
    while (newVar < minVar || newVar > maxVar)
        {
        if (newVar < minVar)
            newVar = minVar * minVar / newVar;
        if (newVar > maxVar)
            newVar = maxVar * maxVar / newVar;
        }
    
    /* store new value */
    (*GetParamVals (param, chain, state[chain])) = newVar;
    
    /* calculate prior ratio */
    for (i=0; i<param->nSubParams; i++)
        {
        rclModel = GetParamIntVals (param->subParams[i], chain, state[chain]);

        if (*rclModel == RCL_TK02)
            {
            tk02Rate = GetParamVals (param->subParams[i], chain, state[chain]);
            t = GetTree (param->subParams[i], chain, state[chain]);
            for (j=0; j<t->nNodes-2; j++)
                {
                p = t->allDownPass[j];
                if (p->length > 0.0)  // not ancestral fossil
                    {
                    (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[p->anc->index], oldVar*p->length, tk02Rate[p->index]);
                    (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->anc->index], newVar*p->length, tk02Rate[p->index]);
                    }
                }
            }
        else if (*rclModel == RCL_IGR)
            {
            igrRate = GetParamVals (param->subParams[i], chain, state[chain]);
            t = GetTree (param->subParams[i], chain, state[chain]);
            for (j=0; j<t->nNodes-2; j++)
                {
                p = t->allDownPass[j];
                if (p->length > 0.0)  // not ancestral fossil
                    {
                    (*lnPriorRatio) -= LnProbGamma (p->length/oldVar, p->length/oldVar, igrRate[p->index]);
                    (*lnPriorRatio) += LnProbGamma (p->length/newVar, p->length/newVar, igrRate[p->index]);
                    }
                }
            }
        }

    /* take prior on Mixedvar into account */
    if (!strcmp(mp->mixedvarPr,"Exponential"))
        (*lnPriorRatio) += mp->mixedvarExp * (oldVar - newVar);
    
    /* calculate proposal ratio */
    (*lnProposalRatio) = log (newVar / oldVar);
    
    /* we do not need to update likelihoods */
    for (i=0; i<param->nRelParts; i++)
        {
        modelSettings[param->relParts[i]].upDateCl = NO;
        }
    
    return (NO_ERROR);
}


int Move_RelaxedClockModel (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* rjMCMC move between TK02 (correlated lognormal) and IGR (independent gamma rate)  
     //chi */
    
    int         i, *rclModel;
    MrBFlt      *mxvar, *mxRate, *brlens, ratio, tk02var, igrvar;
    Tree        *t;
    TreeNode    *p = NULL;
    ModelInfo   *m;
    Model       *mp;

    /* get model settings and parameters */
    m = &modelSettings[param->relParts[0]];
    mp = &modelParams[param->relParts[0]];
    mxvar  = GetParamVals (m->mixedvar, chain, state[chain]);
    mxRate = GetParamVals (param, chain, state[chain]);
    brlens = GetParamSubVals (param, chain, state[chain]);
    t = GetTree (param, chain, state[chain]);
    
    /* get current value of model indicator */
    rclModel = GetParamIntVals(param, chain, state[chain]);

    /* get tk02/igr var ratio */
    ratio = mvp[0];

    (*lnPriorRatio) = (*lnProposalRatio) = 0.0;
    
    /* rjMCMC between models: Pr(TK02) = Pr(IGR) = 1/2 */
    /* the current model is TK02, move to IGR */
    if ((*rclModel) == RCL_TK02)
        {
        /* move the var parameter */
        tk02var = (*mxvar);
     // ratio  *= RandomNumber(seed);
        igrvar  = tk02var / ratio;
        if (igrvar < IGRVAR_MIN || igrvar > IGRVAR_MAX)
            {
            abortMove = YES;
            return (NO_ERROR);
            }
        
        /* take prior on Mixedvar into account */
        if (!strcmp(mp->mixedvarPr,"Exponential"))
            (*lnPriorRatio) += mp->mixedvarExp * (tk02var - igrvar);

        /* match the rates and change the effective branch lengths */
        for (i = 0; i < t->nNodes -2; i++)
            {
            p = t->allDownPass[i];
            if (p->length > 0.0)  // not ancestral fossil
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (mxRate[p->anc->index], tk02var*p->length, mxRate[p->index]);
                (*lnPriorRatio) += LnProbGamma (p->length/igrvar, p->length/igrvar, mxRate[p->index]);

                brlens[p->index] = mxRate[p->index] * p->length;
                }
            }
        
        /* In this move, we simply match the parameters in each model,
           the dimension is same, the Jacobian is 1/ratio. */
        (*lnProposalRatio) -= log(ratio);
            
        /* switch model */
        (*rclModel) = RCL_IGR;
        (*mxvar) = igrvar;
        }
    /* the current model is IGR, move to TK02 */
    else
        {
        /* move the var parameter */
        igrvar  = (*mxvar);
     // ratio  *= RandomNumber(seed);
        tk02var = igrvar * ratio;
        if (tk02var < TK02VAR_MIN || tk02var > TK02VAR_MAX)
            {
            abortMove = YES;
            return (NO_ERROR);
            }

        /* take prior on Mixedvar into account */
        if (!strcmp(mp->mixedvarPr,"Exponential"))
            (*lnPriorRatio) += mp->mixedvarExp * (igrvar - tk02var);
    
        /* match the rates and change the effective branch lengths */
        for (i = 0; i < t->nNodes -2; i++)
            {
            p = t->allDownPass[i];
            if (p->length > 0.0)  // not ancestral fossil
                {
                (*lnPriorRatio) -= LnProbGamma (p->length/igrvar, p->length/igrvar, mxRate[p->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (mxRate[p->anc->index], tk02var*p->length, mxRate[p->index]);

                brlens[p->index] = p->length * (mxRate[p->index] + mxRate[p->anc->index]) /2.0;
                }
            }
            
        /* In this move, we simply match the parameters in each model,
           the dimension is same, the Jacobian is ratio. */
        (*lnProposalRatio) += log(ratio);
            
        /* switch model */
        (*rclModel) = RCL_TK02;
        (*mxvar) = tk02var;
        }

    /* since effective branch lengths are updated, we need to update likelihood calculation */
    TouchAllTreeNodes(m, chain);

    return (NO_ERROR);
    MrBayesPrint ("%lf", *seed); /* just because I am tired of seeing the unused parameter error msg */
}


/*----------------------------------------------------------------
|
|   Move_Local: This proposal mechanism changes the topology and
|      branch lengths of an unrooted tree using the LOCAL mech-
|      anism described by Larget & Simon (1999):
|
|      Larget, B. L., and D. L. Simon. 1999. Markov chain 
|         Monte Carlo algorithms for the Bayesian analysis 
|         of phylogenetic trees. Molecular Biology and 
|         Evolution 16:750-759.
|
|      Programmed by FR 2001-10-14 and partly rewritten 2002-02-21
|      for more elegance and the ability to deal with rooted trees.
|      Support for locked nodes added 2004-01-12 based on mb v2.01.
|      Calculation of the Hastings ratio corrected 2004-07-01.
|      Boundary conditions correctly taken care of 2004-09-29.
|      NB! An alternative to reflection is to skip moves, which might
|          be better for the LOCAL given the complexity of taking
|          the boundary conditions into account
|
----------------------------------------------------------------*/
int Move_Local (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         topologyHasChanged, isVPriorExp, directionUp, moveX;
    MrBFlt      oldM, newM, x, y, newX, newY,
                tuning, minV, maxV, brlensExp=0.0;
    TreeNode    *v, *u, *a, *b, *c, *d;
    Tree        *t;
    ModelParams *mp;

    tuning = mvp[0]; /* Larget & Simon's tuning parameter lambda */
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1];
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);
    
    topologyHasChanged = NO;

#   if defined (DEBUG_LOCAL)
    printf ("Before:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
#   endif
    
    /* pick an internal branch */
    do
        {
        v = t->intDownPass[(int)(RandomNumber(seed)*t->nIntNodes)];
        } while (v->anc->anc == NULL);
        
    /* set up pointers for crown part */
    if (RandomNumber(seed) < 0.5)
        {
        c = v->left;
        d = v->right;
        }
    else
        {
        c = v->right;
        d = v->left;
        }

    /* set up pointers for root part */
    u = v->anc;
    if (RandomNumber(seed) < 0.5 || (t->isRooted == YES && u->anc->anc == NULL))
        {
        directionUp = YES;
        if (u->left == v)
            a = u->right;
        else
            a = u->left;
        b = u->anc;
        }
    else
        {
        directionUp = NO;
        if (u->left == v)
            b = u->right;
        else
            b = u->left;
        a = u->anc;
        }

    /* store old and new path length as well as old x and y */
    oldM = c->length + v->length;
    if (directionUp == YES)
        {
        oldM += a->length;
        x = a->length;
        }
    else
        {
        oldM += u->length;
        x = u->length;
        }

    y = x + v->length;

    /* pick dangly to move */
    if (RandomNumber(seed) < 0.5)
        moveX = YES;
    else
        moveX = NO;

    /* find new m value */
    newM = exp(tuning * (RandomNumber(seed) - 0.5)) * oldM;

    /* pick dangly to move and
       pick new attachment point */
    if (moveX == YES)
        {
        /* choose new x */

        /* first update y */
        newY = y * (newM / oldM);

        /* find reinsertion point */
        if (v->isLocked == YES) 
            {
            newX = RandomNumber(seed) *  newY;
            }
        else
            {
            newX = RandomNumber(seed) * newM;
            }
        }
    else
        {
        /* choose new y */

        /* first update x */
        newX = x * (newM / oldM);

        /* find reinsertion point */
        if (v->isLocked == YES)
            {
            newY = RandomNumber(seed) * (newM - newX) + newX;
            }
        else
            {
            newY = RandomNumber(seed) * newM;
            }
        }

    /* adjust proposal and prior ratio based on length modification */
    /* and insertion mechanism */   
    (*lnProposalRatio) += 3.0 * log (newM / oldM);
    if (isVPriorExp == YES)
        (*lnPriorRatio) = brlensExp * (oldM - newM);

    /* make topology move if necessary and then set branch lengths */
    if (newX > newY)
        {
        /* check if we need to abort */
        if (newY < minV || newY > maxV || (newX-newY) < minV || (newX-newY) > maxV || (newM-newX) < minV || (newM-newX) > maxV)
            {
            abortMove = YES;
            return NO_ERROR;
            }

        /* topology has changed */
        topologyHasChanged = YES;
        /* detach v and d */
        /* this scheme differs from that used by Larget and Simon but is more
           convenient because it avoids tree rotations */
        if (u->left == v)
            u->left = c;
        else
            u->right = c;
        c->anc = u;
        if (directionUp == YES)
            {
            /* place v and d below a */
            if (v->left == d)
                v->right = a;
            else
                v->left = a;
            a->anc = v;
            if (u->left == a)
                u->left = v;
            else
                u->right = v;
            /* v->anc is already u */
            /* adjust lengths */
            c->length = newM - newX;
            v->length = newX - newY;
            a->length = newY;
            }
        else
            {
            /* place v and d below u */
            if (u->isLocked == YES)
                {
                v->isLocked = YES;
                u->isLocked = NO;
                v->lockID = u->lockID;
                u->lockID = 0;
                }
            if (v->left == d)
                v->right = u;
            else
                v->left = u;
            u->anc = v;
            v->anc = a;
            if (a->left == u)
                a->left = v;
            else
                a->right = v;
            /* adjust lengths */
            c->length = newM - newX;
            u->length = newX - newY;
            v->length = newY;
            }
        }
    else
        {
        /* check if we need to abort */
        if (newX < minV || newX > maxV || (newY-newX) < minV || (newY-newX) > maxV || (newM-newY) < minV || (newM-newY) > maxV)
            {
            abortMove = YES;
            return NO_ERROR;
            }

        /* topology has not changed */
        c->length = newM - newY;
        v->length = newY - newX;
        if (directionUp == YES)
            a->length = newX;
        else
            u->length = newX;
        }
                
    /* set update of transition probs */
    c->upDateTi = YES;
    v->upDateTi = YES;
    if (directionUp == YES)
        a->upDateTi = YES;
    else
        u->upDateTi = YES;
        
    /* set flags for update of cond likes from v and u down to root */
    v->upDateCl = YES; 
    u->upDateCl = YES; 
    if (directionUp == YES)
        v = b;
    else
        v = a;
    while (v->anc != NULL)
        {
        v->upDateCl = YES; 
        v = v->anc;
        }

    /* get downpass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }
        
    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);
    
#   if defined (DEBUG_LOCAL)
    printf ("After:\n");
    ShowNodes (t->root, 2, NO);
    getchar();
    printf ("Proposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  c: %d  d: %d  a: %d  b: %d\n",v->index, u->index, 
            c->index, d->index, a->index, b->index);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    getchar();
#   endif

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_LocalClock: This proposal mechanism changes the topology and
|      branch lengths of a rooted tree using the LOCAL (clock) mech-
|      anism described by Larget & Simon (1999):
|
|      Larget, B. L., and D. L. Simon. 1999. Markov chain 
|         Monte Carlo algorithms for the Bayesian analysis 
|         of phylogenetic trees. Molecular Biology and 
|         Evolution 16:750-759.
|
|      Programmed by JH 2002-07-07
|      Modified by FR 2004-05-22 to handle locked and dated trees
|      Modified by FR 2005-11-09 to take care of erroneous
|           Hastings ratio. The fix implemented here results in
|           a move that does not change tree height.
|
----------------------------------------------------------------*/
int Move_LocalClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             i, topologyHasChanged, vIsRoot, aSister, bSister, cSister, *nEvents;
    MrBFlt          x, y, h1, h2, h3, h[3], tempD, ran, distUv, distCv,
                    oldALength, oldBLength, oldCLength, oldULength, oldVLength, lambda, nu,
                    *tk02Rate=NULL, *brlens, oldTreeLength, newTreeLength;
 // MrBFlt          newDepth, oldDepth, factor, diff;
    TreeNode        *u, *v, *w=NULL, *a, *b, *c, *deepestChild, *p;
    Tree            *t;
    Param           *subParm;

    /* get tree */
    t = GetTree (param, chain, state[chain]);

#if defined (DEBUG_LOCAL)
    /* check branch lengths and node depths */
    for (i=0; i<t->nNodes-2; i++) {
        p = t->allDownPass[i];
        /* the two checkings don't consider ancestral fossil (brl=0) in fossilized bd tree */
        if (p->length < minV) {
            printf ("%s   ERROR when entering LocalClock: node %d has length %lf", spacer, p->index, p->length);
            return ERROR;
        }
        if (p->nodeDepth >= p->anc->nodeDepth) {
            printf ("%s   ERROR when entering LocalClock: node %d has depth %lf larger than its ancestor %d depth %lf", spacer, p->index, p->nodeDepth, p->anc->index, p->anc->nodeDepth);
            return ERROR;
        }
    }
#endif

    topologyHasChanged = NO;

#   if defined (DEBUG_LOCAL)
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
#   endif

    /* set up pointers */
    do
        {
        u = t->intDownPass[(int)(RandomNumber(seed)*(t->nIntNodes-1))];
        } while (u->anc->anc == NULL);
    v = u->anc;
    a = u->left;
    b = u->right;
    if (v->left == u)
        c = v->right;
    else
        c = v->left;
    vIsRoot = NO;
    if (v->anc->anc != NULL)
        w = v->anc;
    else
        vIsRoot = YES;
    
    oldALength = a->length;
    oldBLength = b->length;
    oldCLength = c->length;
    oldVLength = v->length;
    oldULength = u->length;
    oldTreeLength = TreeLength (param, chain);
        
    /* get distances from root of move (w or v) to nodes a, b, and c */
    if (vIsRoot == NO)
        h1 = h2 = h3 = v->length;
    else
        h1 = h2 = h3 = 0.0;
    h1 += u->length + a->length;
    h2 += u->length + b->length;
    h3 += c->length;
    h[0] = h1;
    h[1] = h2;
    h[2] = h3;
    
    /* we also need the distances between u <-> v and c <-> v to calculate the hastings' term */
    distUv = u->length;
    distCv = c->length;
        
    /* sort distances (simply make three comparisons and swap values, if necessary) */
    if (h[0] > h[1])
        {
        tempD = h[1];
        h[1] = h[0];
        h[0] = tempD;
        }
    if (h[0] > h[2])
        {
        tempD = h[2];
        h[2] = h[0];
        h[0] = tempD;
        }
    if (h[1] > h[2])
        {
        tempD = h[2];
        h[2] = h[1];
        h[1] = tempD;
        }
        
    /* Find the child node (a, b, or c) that is closest to the root (i.e., has smallest h_i; i=1,2,3). This
       part deals with the possibility that some of the nodes are at the same nodeDepth and randomly assigns
       a node to be deepest in case of ties. */
    if (AreDoublesEqual (h1, h2, 0.00000001) == YES && AreDoublesEqual (h1, h3, 0.00000001) == YES && AreDoublesEqual (h2, h3, 0.00000001) == YES)
        {
        ran = RandomNumber(seed);
        if (ran < 0.33333333)
            deepestChild = a;
        else if (ran > 0.66666666)
            deepestChild = b;
        else
            deepestChild = c;
        }
    else if (AreDoublesEqual (h1, h2, 0.00000001) == YES && AreDoublesEqual (h1, h3, 0.00000001) == NO && AreDoublesEqual (h2, h3, 0.00000001) == NO)
        {
        if (h1 < h3)
            {
            ran = RandomNumber(seed);
            if (ran < 0.5)
                deepestChild = a;
            else
                deepestChild = b;
            }
        else
            deepestChild = c;
        }
    else if (AreDoublesEqual (h1, h2, 0.00000001) == NO && AreDoublesEqual (h1, h3, 0.00000001) == YES && AreDoublesEqual (h2, h3, 0.00000001) == NO)
        {
        if (h1 < h2)
            {
            ran = RandomNumber(seed);
            if (ran < 0.5)
                deepestChild = a;
            else
                deepestChild = c;
            }
        else
            deepestChild = b;
        }
    else if (AreDoublesEqual (h1, h2, 0.00000001) == NO && AreDoublesEqual (h1, h3, 0.00000001) == NO && AreDoublesEqual (h2, h3, 0.00000001) == YES)
        {
        if (h2 < h1)
            {
            ran = RandomNumber(seed);
            if (ran < 0.5)
                deepestChild = b;
            else
                deepestChild = c;
            }
        else
            deepestChild = a;
        }
    else
        {
        if (h1 < h2 && h1 < h3)
            deepestChild = a;
        else if (h2 < h1 && h2 < h3)
            deepestChild = b;
        else
            deepestChild = c;
        }
    
    /* get x and y */
    /* for most of the branches, the proposal ratio is 0.0 so it makes sense to set this first */
    (*lnProposalRatio) = 0.0;
    if (u->isDated == YES && v->isDated == YES)
        {
        /* this proposal is wasted, change nothing */
        if (vIsRoot == NO)
            {
            y = v->length;
            x = y + u->length;
            }
        else
            {
            y = 0.0;
            x = u->length;
            }
        }
    else if (u->isDated == YES && v->isDated == NO)
        {
        /* we can only change the position of v */
        if (vIsRoot == NO)
            {
            /* the upper limit of v's height is determined either by u-length + v->length or by c->length + v->length (h[0]) */
            x = v->length + u->length;
            if (x > h[0])
                x = h[0];
            y = RandomNumber(seed) * x;
            }
        else
            {
            /* v is root: we leave tree height unchanged so we cannot change anything */
            x = u->length;
            y = 0.0;
            }
        }
    else if (u->isDated == NO && v->isDated == YES)
        {
        /* we can only change the position of u */
        if (vIsRoot == NO)
            y = v->length;
        else
            y = 0.0;
        if (u->isLocked == YES)
            {
            if (h1 > h2)
                {
                x = y + RandomNumber(seed) * (h2 - y);
                }
            else
                {
                x = y + RandomNumber(seed) * (h1 - y);
                }
            }
        else
            {
            x = y + RandomNumber(seed) * (h[1] - y);
            }
        }
    /* if we reach the statements down here, neither u nor v is dated */
    else if (u->isLocked == YES)
        {
        if (h1 > h2)
            {
            y = RandomNumber(seed) * h[0];
            x = y + RandomNumber(seed) * (h2 - y);
            }
        else
            {
            y = RandomNumber(seed) * h[0];
            x = y + RandomNumber(seed) * (h1 - y);
            }
        }
    else if (vIsRoot == NO)
        {
        /* this is the standard variant for nonroot v */
        x = RandomNumber(seed) * h[1];
        y = RandomNumber(seed) * h[0];
        }
    else
        {
        /* this is the standard variant when v is the root */
        /*oldDepth = t->root->left->nodeDepth;
          factor = exp((RandomNumber(seed) - 0.5) * 2.0 * log(1.2));
          t->root->left->nodeDepth = newDepth =  factor * h[0] - h[0] + oldDepth;
          adjust h[0], h[1], and h[2] 
          diff = newDepth - oldDepth;
          h[0] += diff;
          h[1] += diff;
          h[2] += diff;*/
        /* set y to 0.0 and select new x */
        y = 0.0;
        x = RandomNumber(seed) * h[1];
        /* Adjust proposal ratio. We deal with topology bias below. Note that this
           proposal ratio is very different from what appeared in Larget and Simon */
        /*(*lnProposalRatio) += (t->nIntNodes-1) * log(oldDepth / newDepth);*/
        /*(*lnProposalRatio) += 2.0 * log (factor);*/
        }
        
    /* decide which topology we will construct (cSister is what we started with) */
    aSister = bSister = cSister = NO;
    /* if u is locked then we cannot change topology */
    if (u->isLocked == YES)
        cSister = YES;
    else if (MaximumValue (x, y) < h[0])
        {
        ran = RandomNumber(seed);
        if (ran < 0.33333333)
            aSister = YES;
        else if (ran > 0.66666666)
            bSister = YES;
        else 
            cSister = YES;
        }
    else
        {
        if (deepestChild == a)
            aSister = YES;
        else if (deepestChild == b)
            bSister = YES;
        else 
            cSister = YES;
        }
    
    /* adjust lengths of nodes u and v */
    u->length = MaximumValue (x, y) - MinimumValue (x, y);
    v->length = MinimumValue (x, y);
    if (vIsRoot == NO)
        v->nodeDepth = w->nodeDepth - v->length;
    u->nodeDepth = v->nodeDepth - u->length;
    
    /* adjust pointers and lengths of nodes a, b, and c */
    topologyHasChanged = NO;
    if (cSister == YES)
        {
        if (v->left == u)
            v->right = c;
        else
            v->left = c;
        u->left = a;
        u->right = b;
        a->anc = b->anc = u;
        c->anc = v;
        a->length = u->nodeDepth - a->nodeDepth;
        b->length = u->nodeDepth - b->nodeDepth;
        c->length = v->nodeDepth - c->nodeDepth;
        }
    else if (bSister == YES)
        {
        if (v->left == u)
            v->right = b;
        else
            v->left = b;
        u->left = a;
        u->right = c;
        a->anc = c->anc = u;
        b->anc = v;
        a->length = u->nodeDepth - a->nodeDepth;
        b->length = v->nodeDepth - b->nodeDepth;
        c->length = u->nodeDepth - c->nodeDepth;
        topologyHasChanged = YES;
        }
    else if (aSister == YES)
        {
        if (v->left == u)
            v->right = a;
        else
            v->left = a;
        u->left = b;
        u->right = c;
        b->anc = c->anc = u;
        a->anc = v;
        a->length = v->nodeDepth - a->nodeDepth;
        b->length = u->nodeDepth - b->nodeDepth;
        c->length = u->nodeDepth - c->nodeDepth;
        topologyHasChanged = YES;
        }

    /* check that all branch lengths are good */
    if (a->length < 0.0 && b->length < 0.0 && c->length < 0.0 && u->length < 0.0 && v->length < 0.0)
        {
        abortMove = YES;
        return NO_ERROR;
        }

    /* calculate the proposal ratio due to asymmetric topology changes */
    if (u->isLocked == NO)
        {
        if (v->isDated == YES || vIsRoot == YES)
            {
            if (distUv > distCv && MaximumValue (x, y) < h[0])
                (*lnProposalRatio) += log(3.0);
            else if (distUv < distCv && MaximumValue (x, y) > h[0])
                (*lnProposalRatio) += log(1.0 / 3.0);
            }
        else
            {
            /* note that Larget and Simon did not have the correct Hastings ratio
               for this case */
            if (distUv > distCv && MaximumValue (x, y) < h[0])
                (*lnProposalRatio) += log(3.0 / 2.0);
            else if (distUv < distCv && MaximumValue (x, y) > h[0])
                (*lnProposalRatio) += log(2.0 / 3.0);
            }
        }

    /* set update of transition probs */
    a->upDateTi = b->upDateTi = c->upDateTi = u->upDateTi = YES;
    if (vIsRoot == NO)
        v->upDateTi = YES;

    /* set flags for update of cond likes from u down to root */
    p = u;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
        
    /* get downpass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        GetDownPass (t);
        
    /* adjust proposal and prior ratio for relaxed clock models */
    newTreeLength = TreeLength(param, chain);
    for (i=0; i<param->subParams[0]->nSubParams; i++)
        {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);
            /* proposal ratio */
            (*lnProposalRatio) += nEvents[a->index] * log (a->length / oldALength);
            (*lnProposalRatio) += nEvents[b->index] * log (b->length / oldBLength);
            (*lnProposalRatio) += nEvents[c->index] * log (c->length / oldCLength);
            (*lnProposalRatio) += nEvents[u->index] * log (u->length / oldULength);
            if (v->anc->anc != NULL && v->isDated == NO)
                (*lnProposalRatio) += nEvents[v->index] * log (v->length / oldVLength);
            /* prior ratio */
            (*lnPriorRatio) += lambda * ((oldTreeLength - newTreeLength)/t->root->left->nodeDepth);
            /* update effective evolutionary lengths */
            if (v->anc->anc == NULL || v->isDated == YES)
                {
                if (UpdateCppEvolLengths (subParm, v->left, chain) == ERROR ||
                    UpdateCppEvolLengths (subParm, v->right, chain) == ERROR)
                    {
                    abortMove = YES;
                    return (NO_ERROR);
                    }
                }
            else
                {
                if (UpdateCppEvolLengths (subParm, v, chain) == ERROR)
                    {
                    abortMove = YES;
                    return (NO_ERROR);
                    }
                }
            }
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            nu /= t->root->left->nodeDepth;     /* variance increase measured relative to tree height */
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            /* no proposal ratio effect */
            /* prior ratio and update of effective evolutionary lengths */
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[u->index], nu*oldALength, tk02Rate[a->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[u->index], nu*u->left->length, tk02Rate[u->left->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[u->index], nu*oldBLength, tk02Rate[b->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[u->index], nu*u->right->length, tk02Rate[u->right->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[v->index], nu*oldCLength, tk02Rate[c->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[v->index], nu*oldULength, tk02Rate[u->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[v->index], nu*v->left->length, tk02Rate[v->left->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[v->index], nu*v->right->length, tk02Rate[v->right->index]);
            brlens[a->index] = a->length * (tk02Rate[a->index] + tk02Rate[a->anc->index])/2.0;
            brlens[b->index] = a->length * (tk02Rate[a->index] + tk02Rate[a->anc->index])/2.0;
            brlens[c->index] = c->length * (tk02Rate[c->index] + tk02Rate[c->anc->index])/2.0;
            brlens[u->index] = u->length * (tk02Rate[u->index] + tk02Rate[u->anc->index])/2.0;
            if (v->anc->anc != NULL && v->isDated == NO)
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[w->index], nu*oldVLength, tk02Rate[v->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[w->index], nu*v->length, tk02Rate[v->index]);
                brlens[v->index] = v->length * (tk02Rate[v->index] + tk02Rate[v->anc->index])/2.0;
                }
            }
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            /* to do */
            }
        }

    /* calculate and adjust prior ratio for clock tree */
    if (LogClockTreePriorRatio (param, chain, &x) == ERROR)
        return (ERROR);
    (*lnPriorRatio) += x;
    
#   if defined (DEBUG_LOCAL)
    printf ("After:\n");
    ShowNodes (t->root, 2, YES);
    printf ("Has topology changed? %d\n",topologyHasChanged);

    /* check branch lengths and node depths */
    for (i=0; i<t->nNodes-2; i++) {
        p = t->allDownPass[i];
        if (p->length < minV) {
            printf ("%s   ERROR when leaving LocalClock: node %d has length %lf", spacer, p->index, p->length);
            return ERROR;
        }
        if (p->nodeDepth >= p->anc->nodeDepth) {
            printf ("%s   ERROR when leaving LocalClock: node %d has depth %lf larger than its ancestor %d depth %lf", spacer, p->index, p->nodeDepth, p->anc->index, p->anc->nodeDepth);
            return ERROR;
        }
    }
#endif

    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


#if 0
/*--------------------------------------------------------------------
|
|   Move_LSPR: Change topology using move based on likelihood scores
|
|--------------------------------------------------------------------*/
int Move_LSPR (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change branch lengths and topology (potentially) using SPR-type move 
       biased according to likelihood scores. NOT work for constrained trees. */
    
    int         i, j, n, division, topologyHasChanged, isVPriorExp, nNodes;
    BitsLong    *pA, *pV, *pP;
    MrBFlt      x, minV, maxV, brlensExp=0.0, minLength=0.0, curLength=0.0, length=0.0,
                cumulativeProb, warpFactor, sum, ran, tuning, increaseProb, decreaseProb,
                divFactor, nStates, rateMult, temp;
    CLFlt       *nSitesOfPat, *globalNSitesOfPat, *tempCondLikes, **tempCondLikePtr;
    TreeNode    *p, *q, *a, *b, *u, *v, *c=NULL, *d, *candidateNodes[20], *vLeft, *vRight;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m = NULL;

    temp = mvp[0];      /* tuning parameter determining how heavily to weight according to likelihood scores */
    var = mvp[1];       /* variance of lognormal for proposing branch lengths */
    increaseProb = decreaseProb = mvp[2];   /* reweighting probabilities */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1];
        isVPriorExp = NO;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

#   if defined (DEBUG_MLSPR)
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
#   endif
    
    /* set topologyHasChanged to NO */
    topologyHasChanged = NO;
    
    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;
        p->x = 0;
        p->d = 0.0;
        }

    /* pick a random branch */
    do
        {
        p = t->allDownPass[(int)(RandomNumber(seed)*(t->nNodes - 1))];
        } while (p->anc->anc == NULL || p->anc->isLocked == YES);
        
    /* set up pointers for nodes around the picked branch */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    vLeft = v->left;
    vRight = vRight;

    /* store the branch lengths */
    aLength = a->length;
    uLength = u->length;
    vLength = v->length;
    if (v->left != NULL)
        {
        vLeftLength = v->left->length;
        vRightLength = v->right->length;
        }
    else
        vLeftLength = vRightLength = 0.0;

    /* get the ML branch lengths */
    /* set initial branch lengths */
    /* cycle through using Newton Raphson and reoptimization a fixed number of iterations */
    for (i=0; i<5; i++)
        {
        }
    
    /* get variance of lognormal */
    
    /* clip tree */
    a->anc = b;
    if (b->left == u)
        b->left = a;
    else
        b->right = a;

    /* count distance to root */
    q = b;
    nNodes = 0;
    while (q->anc != NULL)
        {
        nNodes++;
        q = q->anc;
        }
    
    /* allocate space for temporary cond likes and condlike pointers */
    tempCondLikes = (CLFlt *) SafeCalloc (nNodes*m->numChars*m->numModelStates, sizeof (CLFlt));
    tempCondLikePtr = (CLFlt **) SafeCalloc (nNodes, sizeof (CLFlt *));
    if (!tempCondLikes || !tempCondLikePtr)
        {
        free (tempCondLikes);
        free (tempCondLikePtr);
        return (ERROR);
        }

    /* shift pointers over */
    q = b;
    j = 0;
    while (q->anc != NULL)
        {
        tempCondLikePtr[j] = m->condLike[chain][q->index][q->clSpace];
        m->condLike[chain][q->index][q->clSpace] = tempCondLikes + j*m->numChars*m->numModelStates; 
        j++;
        q = q->anc;
        }

    /* set length to 0.1 for now; test ML brlen later */
    aLength = a->length;
    a->length = 0.1;
    uLength = u->length;
    u->length = 0.1;
    vLength = v->length;
    v->length = 0.1;

    /* get downpass cond likes for the root part */
    q = b;
    while (q->anc != NULL)
        m->condLikeDown (q, division, chain);

    /* get final pass cond likes for the root part */
    GetLikeFPRootPath (a);

    /* get downpass parsimony states for the crown part */
    GetParsDP (t, v, chain);

    /* mark all nodes in the root part of the tree */
    t->root->left->marked = YES;
    for (i=t->nNodes-3; i>=0; i--)
        {
        p = t->allDownPass[i];
        if (p->anc->marked == YES && p != u)
            p->marked = YES;
        }

    /* find number of site patterns and modify randomly */
    globalNSitesOfPat = numSitesOfPat + (chainId[chain] % chainParams.numChains) * numCompressedChars + m->compCharStart;
    nSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!nSitesOfPat)
        {
        MrBayesPrint ("%s   Problem allocating nSitesOfPat in Move_LSPR\n", spacer);
        free (tempCondLikes);
        free (tempCondLikePtr);
        return (ERROR);
        }
    for (i=0; i<numCompressedChars; i++)
        {
        nSitesOfPat[i] = globalNSitesOfPat[i];
        for (j=0; j<globalNSitesOfPat[i]; j++)
            {
            ran = RandomNumber(seed);
            if (ran < decreaseProb)
                nSitesOfPat[i]--;
            else if (ran > 1.0 - increaseProb)
                nSitesOfPat[i]++;
            }
        }

    /* cycle through the possibilities and record ln likelihood of each in p->d */
    minLength = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO)
            continue;
        /* find the parsimony length */
        p->d = 0.0;
        for (n=0; n<t->nRelParts; n++)
            {
            division = t->relParts[n];
            
            /* Find model settings */
            m = &modelSettings[division];

            nStates = m->numModelStates;
            if (m->dataType == STANDARD)
                nStates = 2;
            rateMult = GetRate(division, chain);

            divFactor = warpFactor + log(nStates-1) - log (3) - log(rateMult);

            /* find downpass parsimony sets for the node and its environment */
            pP   = parsPtr[chain][p->index]      + m->parsMatrixStart + Bit(division, p->clSpace)      * parsMatrixRowSize;
            pA   = parsPtr[chain][p->anc->index] + m->parsMatrixStart + Bit(division, p->anc->clSpace) * parsMatrixRowSize;
            pV   = parsPtr[chain][v->index]      + m->parsMatrixStart + Bit(division, v->clSpace)      * parsMatrixRowSize;
        
            length = 0.0;
            for (j=0; j<m->numChars; j++)
                {
                x = (pP[j] | pA[j]) & pV[j];
                if (x == 0)
                    length += nSitesOfPat[j];
                }
            p->d += divFactor * length;
            }
        if (i == 0)
            minLength = p->d;
        else if (p->d < minLength)
            minLength = p->d;
        if (p == a)
            curLength = p->d;
        }

    /* find the sum given the warp factor */
    sum = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES)
            {
            p->d = exp (minLength - p->d);
            sum += p->d;
            }
        }

    /* generate a random uniform */
    ran = RandomNumber(seed);

    /* select the appropriate reattachment point */
    cumulativeProb = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES)
            {
            c = p;
            cumulativeProb += p->d / sum;
            if (cumulativeProb > ran)
                break;
            }
        }
    if (c->marked != YES)
        {
        printf ("Could not select node\n");
        getchar();
        }

    /* calculate the proposal ratio */
    if (c == a)
        (*lnProposalRatio) = 0.0;
    else
        (*lnProposalRatio) = c->d - curLength;

    /* reattach */
    d = c->anc;
    c->anc = u;
    if (u->left == v)
        u->right = c;
    else
        u->left = c;
    if (d->left == c)
        d->left = u;
    else
        d->right = u;
    u->anc = d;

    /* reassign branch lengths */
    if (c != a)
        {
        topologyHasChanged = YES;
        if (RandomNumber(seed) < 0.5)
            {
            x = u->length;
            u->length = a->length;
            a->length = x;
            }
        if (RandomNumber(seed) < 0.5)
            {
            x = c->length;
            c->length = u->length;
            u->length = x;
            }
        /* hit c length with multiplier (a and u dealt with below) */
        x = c->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV)
                x = minV * minV / x;
            else if (x > maxV)
                x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / c->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (c->length - x);
        c->length = x;
        }
    
    /* hit a length with multiplier (even if no topology change was made) */
    x = a->length * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV)
            x = minV * minV / x;
        else if (x > maxV)
            x = maxV * maxV / x;
        }

    /* calculate proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / a->length);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (a->length - x);
    a->length = x;

    /* hit u length with multiplier (even if no topology change was made) */
    x = u->length * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV)
            x = minV * minV / x;
        else if (x > maxV)
            x = maxV * maxV / x;
        }

    /* calculate proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / u->length);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (u->length - x);
    u->length = x;

    /* set tiprobs update flags */
    a->upDateTi = YES;
    u->upDateTi = YES;
    c->upDateTi = YES;  /* could be same as a but that does not matter */

    /* set flags for update of cond likes from u and down to root */
    p = u;
    while (p->anc != NULL)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }

    /* set flags for update of cond likes from b and down to root */
    p = b;
    while (p->anc != NULL && p->upDateCl == NO)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }

    free (nSitesOfPat);

#   if defined (DEBUG_MLSPR)
    printf ("After:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
    printf ("Proposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  a: %d  b: %d\n",v->index, u->index, a->index, b->index);
    printf ("No. nodes moved in root subtree: %d\n",nRootNodes);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    getchar();
#   endif

    return (NO_ERROR);
}


/*--------------------------------------------------------------------
|
|  Move_LSPR1: Change topology using move based on likelihood scores
|
|--------------------------------------------------------------------*/
int Move_LSPR1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change branch lengths and topology (potentially) using SPR-type move 
       biased according to likelihood scores. NOT work for constrained trees. */
    
    int         i, j, n, division, topologyHasChanged, isVPriorExp, nNodes;
    BitsLong    *pA, *pV, *pP;
    MrBFlt      x, minV, maxV, brlensExp=0.0, minLength=0.0, curLength=0.0, length=0.0,
                cumulativeProb, warpFactor, sum, ran, tuning, increaseProb, decreaseProb,
                divFactor, nStates, rateMult, temp;
    CLFlt       *nSitesOfPat, *globalNSitesOfPat, *tempCondLikes, **tempCondLikePtr;
    TreeNode    *p, *q, *a, *b, *u, *v, *c=NULL, *d, *candidateNodes[20], *vLeft, *vRight;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m = NULL;

    temp = mvp[0];      /* tuning parameter determining how heavily to weight according to likelihood scores */
    var = mvp[1];       /* variance of lognormal for proposing branch lengths */
    increaseProb = decreaseProb = mvp[2];   /* reweighting probabilities */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1];
        isVPriorExp = NO;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

#   if defined (DEBUG_MLSPR)
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
#   endif
    
    /* set topologyHasChanged to NO */
    topologyHasChanged = NO;
    
    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;
        p->x = 0;
        p->d = 0.0;
        }

    /* pick a random branch */
    do
        {
        p = t->allDownPass[(int)(RandomNumber(seed)*(t->nNodes - 1))];
        } while (p->anc->anc == NULL || p->anc->isLocked == YES);
        
    /* set up pointers for nodes around the picked branch */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    vLeft = v->left;
    vRight = vRight;

    /* store the branch lengths */
    aLength = a->length;
    uLength = u->length;
    vLength = v->length;
    if (v->left != NULL)
        {
        vLeftLength = v->left->length;
        vRightLength = v->right->length;
        }
    else
        vLeftLength = vRightLength = 0.0;

    /* save DP cond likes */
    /* count distance to root */
    q = b;
    nNodes = 0;
    while (q->anc != NULL)
        {
        nNodes++;
        q = q->anc;
        }
    
    /* allocate space for temporary cond likes and condlike pointers */
    tempCondLikes = (CLFlt *) SafeCalloc (nNodes*m->numChars*m->numModelStates, sizeof (CLFlt));
    tempCondLikePtr = (CLFlt **) SafeCalloc (nNodes, sizeof (CLFlt *));
    if (!tempCondLikes || !tempCondLikePtr)
        {
        free (tempCondLikes);
        free (tempCondLikePtr);
        return (ERROR);
        }

    /* shift pointers over */
    q = b;
    j = 0;
    while (q->anc != NULL)
        {
        tempCondLikePtr[j] = m->condLike[chain][q->index][q->clSpace];
        m->condLike[chain][q->index][q->clSpace] = tempCondLikes + j*m->numChars*m->numModelStates; 
        j++;
        q = q->anc;
        }
    
    /* get cond like uppass up to b */
    getLikeUPRootPath (t, b);

    /* get ML branch lengths */
    NRBrlenOptimizer (t, v, 5, 3);

    /* cycle through using Newton Raphson and reoptimization a fixed number of iterations */
    for (i=0; i<numIterations; i++)
        {
        if (v->left != NULL)
            {
            getBaseLikeUpLeft (t, v);   /* store instead of DP */
            NewtonRaphsonBrlen (t, v->left, chain);
            getBaseLikeUpRight (t, v);
            GetNewtonRaphsonBrlen (t, v->right, chain);
            m->CondLikeDown (v);
            }
        if (u->left == v)
            getBaseLikeUpLeft (t, u);
        else
            getBaseLikeUpRight (t, u);
        NewtonRaphsonBrlen (t, v, chain);
        if (u->left == v)
            getBaseLikeUpRight (t, u);
        else
            getBaseLikeUpLeft (t, u);
        NewtonRaphsonBrlen (t, a->length, chain);
        m->CondLikeDown (t, u);
        if (b->left == u)
            getBaseLikeUpLeft (t, b);
        else
            getBaseLikeUpRight (t, b);
        NewtonRaphsonBrlen (t, u->length, chain);
        getLikeUp(t, u);
        getLikeUp(t, v);
        }
    
    /* get variance of lognormal for forward move */
    f = log (a->length) - log (aLength);
    fvar = f*f;
    f = log (v->length) - log (vLength);
    fvar += f*f;
    f = log (u->length) - log (uLength);
    fvar += f*f;
    if (v->left != NULL)
        {
        f = log (v->left->length) - log (vLeftLength);
        fvar += f*f;
        f = log (v->right->length) - log (vRightLength);
        fvar += f*f;
        fvar /= 5.0;
        }
    else
        fvar /= 3.0;

    /* clip tree */
    a->anc = b;
    if (b->left == u)
        b->left = a;
    else
        b->right = a;

    /* get ML branch length for a */
    NewtonRaphsonBrlen (t, a, chain, 3);

    /* propose new length for a */
    f = PointNormal(RandomNumber(seed));
    f *= fvar;
    f += log (a->length);
    a->length = f;

    /* get downpass cond likes for the root part */
    q = b;
    while (q->anc != NULL)
        m->condLikeDown (q, division, chain);

    /* get uppass cond likes for the root part */
    GetLikeUp (t, t->root->left);

    /* cycle through the possibilities and record ln likelihood of each in p->d */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO)
            continue;
        /* attach crown tree here */
        pLength = p->length;
        /* find ml branch lengths */
        NewtonRaphsonBrlens5 (t, v, chain, 5, 3);
        /* find score */
        m->CondLikeDown (t, v);
        m->CondLikeRoot (t, u);
        m->Likelihood (t, u, &lnL);
        p->d = lnL * warp;
        if (i == 0)
            maxLnL = p->d;
        else if (p->d > maxLnL)
            maxLnL = p->d;
        if (p == a)
            curLnL = p->d;
        /* detach crown tree */
        /* restore p->length */
        p->length = pLength;
        }

    /* find the sum given the warp factor */
    sum = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES)
            {
            p->d = exp (maxLnL - p->d);
            sum += p->d;
            }
        }

    /* generate a random uniform */
    ran = RandomNumber(seed);

    /* select the appropriate reattachment point */
    cumulativeProb = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES)
            {
            c = p;
            cumulativeProb += p->d / sum;
            if (cumulativeProb > ran)
                break;
            }
        }
    if (c->marked != YES)
        {
        printf ("Could not select node\n");
        getchar();
        }

    /* calculate the proposal ratio based on biased reattachment */
    if (c == a)
        (*lnProposalRatio) = 0.0;
    else
        (*lnProposalRatio) = (maxLnL - log(c->d)) - curLnL;

    /* reattach */
    if (c != a)
        topologyHasChanged = YES;
    d = c->anc;
    c->anc = u;
    if (u->left == v)
        u->right = c;
    else
        u->left = c;
    if (d->left == c)
        d->left = u;
    else
        d->right = u;
    u->anc = d;

    /* optimize branch lengths */
    NewtonRaphsonBrlens5 (t, v, chain, 5, 5);

    /* calculate variance of lognormal for back move */
    f = log (a->length) - log (aLength);
    fvarNew = f*f;
    f = log (v->length) - log (vLength);
    fvarNew += f*f;
    f = log (u->length) - log (uLength);
    fvarNew += f*f;
    if (v->left != NULL)
        {
        f = log (v->left->length) - log (vLeftLength);
        fvarNew += f*f;
        f = log (v->right->length) - log (vRightLength);
        fvarNew += f*f;
        fvarNew /= 5.0;
        }
    else
        fvarNew /= 3.0;
    
    /* draw new branch lengths */
    c->length = fvar * PointNormal(RandomNumber(seed)) + log(c->length);
    u->length = fvar * PointNormal(RandomNumber(seed)) + log(u->length);
    v->length = fvar * PointNormal(RandomNumber(seed)) + log(v->length);
    if (v->left != NULL)
        {
        v->left->length = fvar * PointNormal(RandomNumber(seed)) + log(v->left->length);
        v->right->length = fvar * PointNormal(RandomNumber(seed)) + log(v->right->length);
        }

    /* calculate proposal ratio for branch lengths */

    /* set tiprobs update flags */
    a->upDateTi = YES;
    u->upDateTi = YES;
    c->upDateTi = YES;  /* could be same as a but that does not matter */
    v->upDateTi = YES;
    if (v->left != NULL)
        {
        v->left->upDateTi = YES;
        v->right->upDateTi = YES;
        }

    /* set flags for update of cond likes from v and down to root */
    q = v;
    while (q->anc != NULL)
        {
        q->upDateCl = YES; 
        q = q->anc;
        }

    /* set flags for update of cond likes from b and down to root */
    q = b;
    while (q->anc != NULL && q->upDateCl == NO)
        {
        q->upDateCl = YES; 
        q = q->anc;
        }

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }

    /* restore old conditional likelihoods */

#   if defined (DEBUG_MLSPR)
    printf ("After:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
    printf ("Proposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  a: %d  b: %d\n",v->index, u->index, a->index, b->index);
    printf ("No. nodes moved in root subtree: %d\n",nRootNodes);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    getchar();
#   endif

    return (NO_ERROR);
}
#endif


/* Move_NNI, change topology using NNI move */
int Move_NNI (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    TreeNode    *p, *u, *v, *a, *b, *c;
    Tree        *t;
    
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* pick an internal branch */
    do  {
        p = t->intDownPass[(int)(RandomNumber(seed) * t->nIntNodes)];
        }
    while (p->anc->anc == NULL || p->isLocked == YES);
        
    /* set up area of rearrangement */
    u = p;
    v = u->anc;
    a = u->left;
    b = u->right;
    if (v->left == u)
        c = v->right;
    else
        c = v->left;
        
    /* change topology */
    if (RandomNumber(seed) < 0.5)
        {
        if (v->left == u)
            v->right = b;
        else
            v->left = b;
        u->left = a;
        u->right = c;
        a->anc = c->anc = u;
        b->anc = v;
        }
    else
        {
        if (v->left == u)
            v->right = a;
        else
            v->left = a;
        u->left = b;
        u->right = c;
        b->anc = c->anc = u;
        a->anc = v;
        }

    /* set update of cond likes */
    while (p->anc != NULL)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }
    
    GetDownPass (t);
    
    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


int Move_NNIClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change clock tree using NNI move */
    
    int         i, *nEvents, numFreeOld, numFreeNew;
    MrBFlt      x, *tk02Rate=NULL, *brlens, *igrRate=NULL, igrvar=0.0, nu=0.0, oldALength, oldCLength;
    TreeNode    *p, *q, *a, *c, *u, *v;
    Tree        *t;
    Param       *subParm;

    /* no tuning parameter */

    /* make absolutely sure the proposal ratio and prior ratio are reset */
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);

#   if defined (DEBUG_NNIClock)
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
#   endif
    
    /* count number of free interior branches */
    numFreeOld = 0;
    for (i=0; i<t->nIntNodes-1; i++)
        {
        p = t->intDownPass[i];
        if (p->anc->left == p)
            q = p->anc->right;
        else
            q = p->anc->left;
        if (p->isLocked == NO && p->nodeDepth >= q->nodeDepth + BRLENS_MIN)
            numFreeOld++;
        }

    /* In extremely constrained trees, it might be impossible to change the tree before nodes have changed in position */
    if (numFreeOld == 0)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* pick an interior branch, around which it is possible to make an NNI */
    do  {
        p = t->intDownPass[(int)(RandomNumber(seed) * (t->nIntNodes-1))];
        if (p->anc->left == p)
            q = p->anc->right;
        else
            q = p->anc->left;
        }
    while (p->isLocked == YES || p->nodeDepth < q->nodeDepth + BRLENS_MIN);
        
    /* set up pointers for nodes around the picked branch */
    /* consider ancestral fossil (brl=0) in fossilized bd tree */
    if (p->left->length < TIME_MIN)
        a = p->right;
    else if (p->right->length < TIME_MIN)
        a = p->left;
    else if (RandomNumber(seed) < 0.5)
        a = p->left;
    else
        a = p->right;
    v = p;
    u = p->anc;
    if (u->left == v)
        c = u->right;
    else
        c = u->left;

    /* record branch lengths */
    oldALength = a->length;
    oldCLength = c->length;
    
    /* make topology change */
    a->anc = u;
    c->anc = v;
    if (v->left == a)
        v->left = c;
    else
        v->right = c;
    if (u->left == c)
        u->left = a;
    else
        u->right = a;

    /* adjust branch lengths */
    a->length = u->nodeDepth - a->nodeDepth;
    c->length = v->nodeDepth - c->nodeDepth;
    assert (a->length > BRLENS_MIN);
    assert (c->length > BRLENS_MIN);

    /* no reassignment of CPP events or branch rates necessary */

    /* set tiprobs update flags */
    a->upDateTi = YES;
    c->upDateTi = YES;

    /* set flags for update of cond likes from v and down to root */
    p = v;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }

    /* get down pass sequence */
    GetDownPass (t);

    /* count number of free interior branches after the move */
    numFreeNew = 0;
    for (i=0; i<t->nIntNodes-1; i++)
        {
        p = t->intDownPass[i];
        if (p->anc->left == p)
            q = p->anc->right;
        else
            q = p->anc->left;
        if (p->isLocked == NO && p->nodeDepth >= q->nodeDepth + BRLENS_MIN)
            numFreeNew++;
        }
        
    /* get proposal ratio if number of free branches has changed */
    if (numFreeNew != numFreeOld)
        (*lnProposalRatio) = log((MrBFlt)numFreeOld / (MrBFlt)numFreeNew);

    /* calculate and adjust prior ratio for clock trees */
    if (LogClockTreePriorRatio(param, chain, &x) == ERROR)
        return (ERROR);
    (*lnPriorRatio) += x;

    /* adjust proposal and prior ratio for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
        {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            /* proposal ratio */
            (*lnProposalRatio) += nEvents[a->index] * log (a->length / oldALength);
            (*lnProposalRatio) += nEvents[c->index] * log (c->length / oldCLength);
            /* prior ratio: no effect because tree length is the same */
            /* update effective evolutionary lengths */
            if (UpdateCppEvolLengths (subParm, a, chain) == ERROR || UpdateCppEvolLengths (subParm, c, chain) == ERROR)
                {
                abortMove = YES;
                return (NO_ERROR);
                }
            }
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            /* prior ratio and update of effective evolutionary lengths */
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[v->index], nu*oldALength, tk02Rate[a->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[u->index], nu* a->length, tk02Rate[a->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[u->index], nu*oldCLength, tk02Rate[c->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[v->index], nu* c->length, tk02Rate[c->index]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = a->length * (tk02Rate[a->index] + tk02Rate[a->anc->index])/2.0;
            brlens[c->index] = c->length * (tk02Rate[c->index] + tk02Rate[c->anc->index])/2.0;
            }
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            /* prior ratio and update of effective evolutionary lengths */
            (*lnPriorRatio) -= LnProbGamma (oldALength/igrvar, oldALength/igrvar, igrRate[a->index]);
            (*lnPriorRatio) -= LnProbGamma (oldCLength/igrvar, oldCLength/igrvar, igrRate[c->index]);
            (*lnPriorRatio) += LnProbGamma (a->length /igrvar, a->length /igrvar, igrRate[a->index]);
            (*lnPriorRatio) += LnProbGamma (c->length /igrvar, c->length /igrvar, igrRate[c->index]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = igrRate[a->index] * a->length;
            brlens[c->index] = igrRate[c->index] * c->length;
            }
        }
    
#   if defined (DEBUG_NNIClock)
    printf ("After:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
    printf ("Proposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  a: %d  b: %d\n",v->index, u->index, a->index, b->index);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    getchar();
#   endif

    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


/* Move_NNI_Hetero, change topology with unlinked brlens using NNI */
int Move_NNI_Hetero (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, brIndex, moveType;
    TreeNode    *p, *u, *v, *a, *b, *c;
    Tree        *t;
    
    (*lnPriorRatio) = (*lnProposalRatio) = 0.0;
        
    /* get first tree */
    t = GetTree (param, chain, state[chain]);

    /* pick an internal branch */
    do
        {
        brIndex = (int) (RandomNumber(seed) * t->nIntNodes);
        p = t->intDownPass[brIndex];
        } while (p->anc->anc == NULL || p->isLocked == YES);
        
    /* decide on how to change the tree */
    if (RandomNumber(seed) < 0.5)
        moveType = 0;
    else
        moveType = 1;
    
    /* cycle through trees */
    for (i=0; i<param->nSubParams; i++)
        {
        /* get tree */
        t = GetTree (param->subParams[i], chain, state[chain]);
            
        /* find p */
        p = t->intDownPass[brIndex];

        /* set up area of rearrangement */
        u = p;
        v = u->anc;
        a = u->left;
        b = u->right;
        if (v->left == u)
            c = v->right;
        else
            c = v->left;

        /* change topology */
        if (moveType == 0)
            {
            if (v->left == u)
                v->right = b;
            else
                v->left = b;
            u->left = a;
            u->right = c;
            a->anc = c->anc = u;
            b->anc = v;
            }
        else if (moveType == 1)
            {
            if (v->left == u)
                v->right = a;
            else
                v->left = a;
            u->left = b;
            u->right = c;
            b->anc = c->anc = u;
            a->anc = v;
            }

        /* set update of ti probs */
        a->upDateTi = YES;
        b->upDateTi = YES;
        c->upDateTi = YES;
        u->upDateTi = YES;
        v->upDateTi = YES;
        
        /* set update of conditional likelihoods */
        while (p->anc != NULL)
            {
            p->upDateCl = YES; 
            p = p->anc;
            }

        /* reset tree downpass sequences */
        GetDownPass (t);
        
        }
    
    return (NO_ERROR);
    MrBayesPrint ("%lf", *mvp); /* just because I am tired of seeing the unused parameter error msg */
}


/*-----------------------------------------------------------------------------------
|
|   Move_NodeSlider: move the position of one node without changing topology
|
-------------------------------------------------------------------------------------*/
int Move_NodeSlider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    MrBFlt      tuning, maxV, minV, oldM, newM, brlensPrExp=0.0, newMin, newMax, oldMin, oldMax;
    TreeNode    *p, *q;
    ModelParams *mp;
    Tree        *t;
    int isVPriorExp;
    
    tuning = mvp[0]; /* Larget & Simon's tuning parameter lambda */

    mp = &modelParams[param->relParts[0]];

    /* max and min brlen (time) */
    if (param->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1];
        isVPriorExp = NO;
        }
    else if (param->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensPrExp = mp->brlensExp;
        isVPriorExp = YES;
        }
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);

    /* pick an interior branch */
    do  {
        p = t->intDownPass[(int)(RandomNumber(seed) * t->nIntNodes)];
        }
    while (p->anc == NULL || (t->isRooted == YES && p->anc->anc == NULL));

    /* pick one descendant branch */
    if (RandomNumber(seed) < 0.5)
        q = p->left;
    else
        q = p->right;
    
    /* determine new length */
    oldM = (q->length + p->length);
    newM = oldM * exp(tuning * (RandomNumber(seed) - 0.5));
    while (newM < 2.0 * minV || newM > 2.0 * maxV)
        {
        if (newM < 2.0 * minV)
            newM = 4.0 * minV * minV / newM;
        else if (newM > 2.0 * maxV)
            newM = 4.0 * maxV * maxV / newM;
        }

    /* determine new lengths of p and q */
    newMin = minV > newM - maxV ? minV : newM - maxV;
    newMax = maxV < newM - minV ? maxV : newM - minV;
    oldMin = minV > oldM - maxV ? minV : oldM - maxV;
    oldMax = maxV < oldM - minV ? maxV : oldM - minV;

    q->length = newMin + RandomNumber(seed) * (newMax - newMin);
    p->length = newM - q->length;

    /* the proposal ratio for two sliding windows */
    (*lnProposalRatio) = log ((newMax - newMin) / (oldMax - oldMin));
    
    /* The proposal ratio for shrinking/expanding two variables (x1 = p->length, x2 = q->length)
       by the same factor c = newM/oldM is c^2. This can be derived by variable transformation:
       y1 = x1, y2 = x2/x1. The proposal ratio in the transformed variables is c, the Jacobian is y1,
       so the proposal ratio in the original variables is c*c = c^2.
       (see Move_TreeLen for m variables and Yang 2006 CME P171 S5.4.4 for details) */
    (*lnProposalRatio) += 2.0 * log(newM / oldM);

    /* set flags for update of transition probabilities at p and q */
    p->upDateTi = YES;
    q->upDateTi = YES;
    p->upDateCl = YES;

    /* set flags for update of cond likes from p->anc and down to root */
    while (p->anc->anc != NULL)
        {
        p = p->anc;
        p->upDateCl = YES;
        }

    /* update prior if exponential prior on branch lengths */
    if (param->paramId == BRLENS_EXP)
        (*lnPriorRatio) = brlensPrExp * (oldM - newM);
    /* Dirichlet or twoExp prior */
    else if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);

    return (NO_ERROR);
    
}


/*-----------------------------------------------------------------------------------
|
|   Move_NodeSliderClock: Move the position of one (root or nonroot) node in clock tree.
|      In calibrated trees, we need to move also calibrated terminal nodes.
|
-------------------------------------------------------------------------------------*/
int Move_NodeSliderClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, *nEvents;
    MrBFlt      window, minDepth, maxDepth, oldDepth, newDepth, minL, minR,
                oldLeftLength=0.0, oldRightLength=0.0, oldPLength=0.0, x, clockRate,
                lambda=0.0, nu=0.0, igrvar=0.0, *brlens=NULL, *tk02Rate=NULL, *igrRate=NULL;
    TreeNode    *p, *q;
    ModelParams *mp;
    ModelInfo   *m;
    Tree        *t;
    Param       *subParm;
    Calibration *calibrationPtr;

    window = mvp[0]; /* window size */
 
    m = &modelSettings[param->relParts[0]];
    mp = &modelParams[param->relParts[0]];

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get clock rate */
    if (m->clockRate == NULL)
        clockRate = 1.0;
    else
        clockRate = *GetParamVals (m->clockRate, chain, state[chain]);

    /* check whether or not we can change root */
    if ((t->root->left->isDated == YES && t->root->left->calibration->prior == fixed) ||
        ((!strcmp(mp->clockPr, "Uniform") || !strcmp(mp->clockPr, "Fossilization")) && mp->treeAgePr.prior == fixed))
        i = t->nNodes - 2;
    else
        i = t->nNodes - 1;

    /* pick a node that can be changed in position */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * i)];
        }
    while ((p->left == NULL && p->isDated == NO) ||
           (p->left == NULL && p->length < TIME_MIN) ||
           (p->isDated == YES && p->calibration->prior == fixed));

#if defined (DEBUG_CSLIDER)
    printf ("Before node slider (clock):\n");
    printf ("Picked branch with index %d and depth %f\n", p->index, p->nodeDepth);
    if (p->anc->anc == NULL)
        printf ("Old clock rate: %f\n", clockRate);
    ShowNodes (t->root, 0, t->isRooted);
    getchar();
#endif

    /* store values needed later for prior calculation (relaxed clocks) */
    oldPLength = p->length;
    if (p->left != NULL)
        {
        oldLeftLength = p->left->length;
        oldRightLength = p->right->length;
        }
    else
        oldLeftLength = oldRightLength = 0.0;

    /* determine lower and upper bound */
    if (p->left == NULL)
        minDepth = 0.0;
    else      // internal node
        {
        if (p->left->length > 0.0)
            minL = p->left->nodeDepth + BRLENS_MIN;
        else  // ancestral fossil
            {
            assert (p->left->calibration != NULL);
            minL = p->left->calibration->min * clockRate;
            }
        if (p->right->length > 0.0)
            minR = p->right->nodeDepth + BRLENS_MIN;
        else  // ancestral fossil
            {
            assert (p->right->calibration != NULL);
            minR = p->right->calibration->min * clockRate;
            }
        if (minL > minR)
            minDepth = minL;
        else
            minDepth = minR;
        }
    
    if (p->anc->anc == NULL)
        maxDepth = TREEHEIGHT_MAX;
    else
        maxDepth = p->anc->nodeDepth - BRLENS_MIN;
    if (p->left != NULL && p->left->length < TIME_MIN)
        {
        assert (p->left->calibration != NULL);
        if (maxDepth > p->left->calibration->max * clockRate)
            maxDepth = p->left->calibration->max * clockRate;
        }
    if (p->right != NULL && p->right->length < TIME_MIN)
        {
        assert (p->right->calibration != NULL);
        if (maxDepth > p->right->calibration->max * clockRate)
            maxDepth = p->right->calibration->max * clockRate;
        }
    
    if (p->isDated == YES)
        calibrationPtr = p->calibration;
    else if (p->anc->anc == NULL && (!strcmp(mp->clockPr,"Uniform") || !strcmp(mp->clockPr, "Fossilization")))
        calibrationPtr = &mp->treeAgePr;
    else
        calibrationPtr = NULL;
    if (calibrationPtr != NULL)
        {
        if (maxDepth > calibrationPtr->max * clockRate)
            maxDepth = calibrationPtr->max * clockRate;
        if (minDepth < calibrationPtr->min * clockRate)
            minDepth = calibrationPtr->min * clockRate;
        }

    /* abort if impossible */
    if (minDepth > maxDepth -BRLENS_MIN)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* save some reflection time */
    if (maxDepth-minDepth < window)
        {
        window = maxDepth-minDepth;
        }

    /* pick the new node depth */
    oldDepth = p->nodeDepth;
    newDepth = oldDepth + (RandomNumber(seed) - 0.5) * window;
 
    /* reflect the new node depth */
    while (newDepth < minDepth || newDepth > maxDepth)
        {
        if (newDepth < minDepth)
            newDepth = 2.0 * minDepth - newDepth;
        if (newDepth > maxDepth)
            newDepth = 2.0 * maxDepth - newDepth;
        }

    p->nodeDepth = newDepth;

    /* determine new branch lengths around p and set update of transition probabilities */
    if (p->left != NULL)
        {
        if (p->left->length > 0.0) {
            p->left->length = p->nodeDepth - p->left->nodeDepth;
            p->left->upDateTi = YES;
            }
        else
            p->left->nodeDepth = p->nodeDepth;
        if (p->right->length > 0.0) {
            p->right->length = p->nodeDepth - p->right->nodeDepth;
            p->right->upDateTi = YES;
            }
        else
            p->right->nodeDepth = p->nodeDepth;
        }
    if (p->anc->anc != NULL)
        {
        p->length = p->anc->nodeDepth - p->nodeDepth;
        p->upDateTi = YES;
        }

    /* adjust age of p if dated */
    if (calibrationPtr != NULL)
        p->age = p->nodeDepth / clockRate;
    if ((p->left != NULL) && (p->left->length < TIME_MIN))
        p->left->age = p->nodeDepth / clockRate;
    if ((p->right != NULL) && (p->right->length < TIME_MIN))
        p->right->age = p->nodeDepth / clockRate;

    /* set flags for update of cond likes from p and down to root */
    q = p;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }

    /* calculate proposal ratio */
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* calculate and adjust prior ratio for clock tree */
    if (LogClockTreePriorRatio (param, chain, &x) == ERROR)
        return (ERROR);
    (*lnPriorRatio) += x;

    /* adjust proposal and prior ratio for relaxed clock models */
    for (i=0; i<param->nSubParams; i++)
        {
        subParm = param->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);
            /* proposal ratio */
            if (p->left != NULL)
                {
                (*lnProposalRatio) += nEvents[p->left->index ] * log (p->left->length  / oldLeftLength);
                (*lnProposalRatio) += nEvents[p->right->index] * log (p->right->length / oldRightLength);
                }
            if (p->anc->anc != NULL)
                (*lnProposalRatio) += nEvents[p->index] * log (p->length / oldPLength);

            /* prior ratio */
            if (p->anc->anc == NULL) // two branches changed in same direction
                (*lnPriorRatio) += lambda * (2.0 * (oldDepth - newDepth));
            else if (p->left != NULL) // two branches changed in one direction, one branch in the other direction
                (*lnPriorRatio) += lambda * (oldDepth - newDepth);
            else /* if (p->left == NULL) */ // one branch changed
                (*lnPriorRatio) += lambda * (newDepth - oldDepth);

            /* update effective evolutionary lengths */
            if (UpdateCppEvolLengths (subParm, p, chain) == ERROR)
                {
                abortMove = YES;
                return (NO_ERROR);
                }
            }
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);

            /* prior ratio & update effective evolutionary lengths */
            if (p->left != NULL)
                {
                if (p->left->length > 0.0)
                    {
                    (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[p->index], nu*oldLeftLength, tk02Rate[p->left->index]);
                    (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->index], nu*p->left->length, tk02Rate[p->left->index]);
                    brlens[p->left->index] = p->left->length * (tk02Rate[p->left->index]+tk02Rate[p->index])/2.0;
                    }
                if (p->right->length > 0.0)
                    {
                    (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[p->index], nu*oldRightLength, tk02Rate[p->right->index]);
                    (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->index], nu*p->right->length, tk02Rate[p->right->index]);
                    brlens[p->right->index] = p->right->length * (tk02Rate[p->right->index]+tk02Rate[p->index])/2.0;
                    }
                }
            if (p->anc->anc != NULL)
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[p->anc->index], nu*oldPLength, tk02Rate[p->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->anc->index], nu*p->length, tk02Rate[p->index]);
                brlens[p->index] = p->length * (tk02Rate[p->index]+tk02Rate[p->anc->index])/2.0;
                }
            }
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
                
            /* prior ratio & update effective evolutionary lengths */
            if (p->left != NULL)
                {
                if (p->left->length > 0.0)
                    {
                    (*lnPriorRatio) -= LnProbGamma (oldLeftLength/igrvar, oldLeftLength/igrvar, igrRate[p->left->index ]);
                    (*lnPriorRatio) += LnProbGamma (p->left->length/igrvar, p->left->length/igrvar, igrRate[p->left->index ]);
                    brlens[p->left->index ] = igrRate[p->left->index ] * p->left->length;
                    }
                if (p->right->length > 0.0)
                    {
                    (*lnPriorRatio) -= LnProbGamma (oldRightLength/igrvar, oldRightLength/igrvar, igrRate[p->right->index]);
                    (*lnPriorRatio) += LnProbGamma (p->right->length/igrvar, p->right->length/igrvar, igrRate[p->right->index]);
                    brlens[p->right->index] = igrRate[p->right->index] * p->right->length;
                    }
                }
            if (p->anc->anc != NULL)
                {
                (*lnPriorRatio) -= LnProbGamma (oldPLength/igrvar, oldPLength/igrvar, igrRate[p->index]);
                (*lnPriorRatio) += LnProbGamma (p->length /igrvar, p->length /igrvar, igrRate[p->index]);
                brlens[p->index] = igrRate[p->index] * p->length;
                }
            }
        }

#if defined (DEBUG_CSLIDER)
    printf ("After node slider (clock):\n");
    printf ("Old depth: %f -- New depth: %f -- LnPriorRatio %f -- LnProposalRatio %f\n",
        oldDepth, newDepth, (*lnPriorRatio), (*lnProposalRatio));
    ShowNodes (t->root, 0, t->isRooted);
    getchar();
#endif

    return (NO_ERROR);
}


int Move_Nu (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move the variance of the TK02 relaxed clock lognormal using multiplier */

    int         i, j;
    MrBFlt      oldNu, newNu, minNu, maxNu, tuning, *tk02Rate;
    Model       *mp;
    TreeNode    *p;
    Tree        *t;

    /* get tuning parameter */
    tuning = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get the min and max values */
    minNu = TK02VAR_MIN;
    maxNu = TK02VAR_MAX;
    if (!strcmp(mp->tk02varPr,"Uniform"))
        {
        minNu = (mp->tk02varUni[0] < TK02VAR_MIN) ? TK02VAR_MIN : mp->tk02varUni[0];
        maxNu = (mp->tk02varUni[1] > TK02VAR_MAX) ? TK02VAR_MAX : mp->tk02varUni[1];
        }
    
    /* get the TK02 lognormal variance */
    oldNu = *GetParamVals (param, chain, state[chain]);

    /* set new value */
    newNu = oldNu * exp ((0.5 - RandomNumber(seed))*tuning);
    
    /* reflect if necessary */
    while (newNu < minNu || newNu > maxNu)
        {
        if (newNu < minNu)
            newNu = minNu * minNu / newNu;
        if (newNu > maxNu)
            newNu = maxNu * maxNu / newNu;
        }
    
    /* store new value */
    (*GetParamVals (param, chain, state[chain])) = newNu;

    /* calculate prior ratio */
    for (i=0; i<param->nSubParams; i++)
        {
        tk02Rate = GetParamVals (param->subParams[i], chain, state[chain]);
        t = GetTree (param->subParams[i], chain, state[chain]);
        for (j=0; j<t->nNodes-2; j++)
            {
            p = t->allDownPass[j];
            if (p->length > 0.0)  // not ancestral fossil
                {
                (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[p->anc->index], oldNu*p->length, tk02Rate[p->index]);
                (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->anc->index], newNu*p->length, tk02Rate[p->index]);
                }
            }
        }

    /* take prior on nu into account */
    if (!strcmp(mp->tk02varPr,"Exponential"))
        (*lnPriorRatio) += mp->tk02varExp * (oldNu - newNu);
    
    /* calculate proposal ratio */
    (*lnProposalRatio) = log (newNu / oldNu);

    /* we do not need to update likelihoods */
    for (i=0; i<param->nRelParts; i++)
        {
        modelSettings[param->relParts[i]].upDateCl = NO;
        }

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_Omega: Change the nonysnonymous/synonymous rate ratio
|      Note that this is appropriate when omegavar=equal
|
----------------------------------------------------------------*/
int Move_Omega (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change omega using sliding window */
    
    int         i, isValidO;
    MrBFlt      oldO, newO, window, minO, maxO, ran, *alphaDir, oldPropRatio, newPropRatio, x, y;
    ModelParams *mp;

    /* get size of window, centered on current omega value */
    window = mvp[0];
    
    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* get minimum and maximum values for omega */
    minO = OMEGA_MIN;
    maxO = OMEGA_MAX;

    /* get old value of omega */
    oldO = *GetParamVals(param, chain, state[chain]);

    /* get Dirichlet parameters */
    alphaDir = mp->tRatioDir;

    /* change value for omega */
    ran = RandomNumber(seed);
    if (maxO-minO < window)
        {
        window = maxO-minO;
        }
    newO = oldO + window * (ran - 0.5);
    
    /* check that new value is valid */
    isValidO = NO;
    do  {
        if (newO < minO)
            newO = 2.0 * minO - newO;
        else if (newO > maxO)
            newO = 2.0 * maxO - newO;
        else
            isValidO = YES;
        }
    while (isValidO == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* get prior ratio from Dirichlet */
    oldPropRatio = oldO / (oldO + 1.0);
    newPropRatio = newO / (newO + 1.0);
    x = ((alphaDir[0] - 1.0) * log(newPropRatio)) + ((alphaDir[1] - 1.0) * log (1.0 - newPropRatio));
    y = ((alphaDir[0] - 1.0) * log(oldPropRatio)) + ((alphaDir[1] - 1.0) * log (1.0 - oldPropRatio));
    (*lnPriorRatio) = x - y;
    
    /* copy new omega value back */
    *GetParamVals(param, chain, state[chain]) = newO;

    /* Set update flags for all partitions that share this kappa. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for covarion, doublet, and codon models, we do want to update
       the cijk flag. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_Omega_M: Change the nonysnonymous/synonymous rate ratio
|      using multiplier. Note that this is appropriate when
|      omegavar=equal
|
----------------------------------------------------------------*/
int Move_Omega_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change omega using multiplier */
    
    int         i, isValidO;
    MrBFlt      oldO, newO, minO, maxO, tuning, ran, factor, *alphaDir, oldPropRatio, newPropRatio, x, y;
    ModelParams *mp;

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get tuning parameter */
    tuning = mvp[0];

    /* get minimum and maximum values for omega */
    minO = OMEGA_MIN;
    maxO = OMEGA_MAX;

    /* get old value of omega */
    oldO = *GetParamVals(param, chain, state[chain]);

    /* get Dirichlet parameters */
    alphaDir = mp->omegaDir;

    /* change value for omega */
    ran = RandomNumber(seed);
    factor = exp(tuning * (ran - 0.5));
    newO = oldO * factor;
    
    /* check that new value is valid */
    isValidO = NO;
    do
        {
        if (newO < minO)
            newO = minO * minO / newO;
        else if (newO > maxO)
            newO = maxO * maxO / newO;
        else
            isValidO = YES;
        } while (isValidO == NO);

    /* get proposal ratio */
    *lnProposalRatio = log(newO / oldO);
    
    /* get prior ratio from Dirichlet */
    oldPropRatio = oldO / (oldO + 1.0);
    newPropRatio = newO / (newO + 1.0);
    x = ((alphaDir[0] - 1.0) * log(newPropRatio)) + ((alphaDir[1] - 1.0) * log (1.0 - newPropRatio));
    y = ((alphaDir[0] - 1.0) * log(oldPropRatio)) + ((alphaDir[1] - 1.0) * log (1.0 - oldPropRatio));
    (*lnPriorRatio) = x - y;
    
    /* copy new omega value back */
    *GetParamVals(param, chain, state[chain]) = newO;

    /* Set update flags for all partitions that share this omega. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for covarion, doublet, and codon models, we do want to update
       the cijk flag. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_OmegaBeta_M: Change parameters of the beta distribution
|      using multiplier for the M10 model. Note that this is 
|      appropriate when omegavar=M10
|
----------------------------------------------------------------*/
int Move_OmegaBeta_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, isValidVal, whichParam;
    MrBFlt      oldVal, newVal, minVal, maxVal, *vals, *subVals, tuning, ran, factor;
    ModelParams *mp;
    
    /* do we pick alpha or beta of the Beta distribution to change */
    if (RandomNumber(seed) < 0.5)
        whichParam = 0;
    else
        whichParam = 1;

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get tuning parameter */
    tuning = mvp[0];

    /* get minimum and maximum values for omega */
    minVal = 0.05;
    maxVal = 100.0;

    /* get old value of omega */
    vals = GetParamVals(param, chain, state[chain]);
    subVals = GetParamSubVals(param, chain, state[chain]);
    oldVal = subVals[mp->numM10BetaCats + mp->numM10GammaCats + 4 + whichParam];

    /* change value for alpha/beta */
    ran = RandomNumber(seed);
    factor = exp(tuning * (ran - 0.5));
    newVal = oldVal * factor;
    
    /* check that new value is valid */
    isValidVal = NO;
    do
        {
        if (newVal < minVal)
            newVal = minVal * minVal / newVal;
        else if (newVal > maxVal)
            newVal = maxVal * maxVal / newVal;
        else
            isValidVal = YES;
        } while (isValidVal == NO);

    /* get proposal ratio */
    *lnProposalRatio = log(newVal / oldVal);
    
    /* get prior ratio */
    if (!strcmp(mp->m10betapr, "Exponential"))
        (*lnPriorRatio) = (log(mp->m10betaExp) - newVal * mp->m10betaExp) - (log(mp->m10betaExp) - oldVal * mp->m10betaExp);
    else
        (*lnPriorRatio) = 0.0;
    
    /* copy new omega value back */
    subVals[mp->numM10BetaCats + mp->numM10GammaCats + 4 + whichParam] = newVal;
    
    /* update the omega values */
    BetaBreaks (subVals[mp->numM10BetaCats + mp->numM10GammaCats + 4], subVals[mp->numM10BetaCats + mp->numM10GammaCats + 5], &vals[0], mp->numM10BetaCats);

    /* Set update flags for all partitions that share this kappa. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for covarion, doublet, and codon models, we do want to update
       the cijk flag. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_OmegaGamma_M: Change parameters of the gamma distribution
|      using multiplier for the M10 model. Note that this is 
|      appropriate whenomegavar=M10
|
----------------------------------------------------------------*/
int Move_OmegaGamma_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, isValidVal, whichParam;
    MrBFlt      oldVal, newVal, minVal, maxVal, *vals, *subVals, tuning, ran, factor, quantile95;
    ModelParams *mp;
    
    /* do we pick alpha or beta of the Gamma distribution to change */
    if (RandomNumber(seed) < 0.5)
        whichParam = 0;
    else
        whichParam = 1;

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get tuning parameter */
    tuning = mvp[0];

    /* get minimum and maximum values for omega */
    minVal = 0.05;
    maxVal = 100.0;

    /* get values */
    vals = GetParamVals(param, chain, state[chain]);
    subVals = GetParamSubVals(param, chain, state[chain]);
    oldVal = subVals[mp->numM10BetaCats + mp->numM10GammaCats + 6 + whichParam];

    /* change value for alpha/beta */
    do
        {
        ran = RandomNumber(seed);
        factor = exp(tuning * (ran - 0.5));
        newVal = oldVal * factor;
        
        /* check that new value is valid */
        isValidVal = NO;
        do
            {
            if (newVal < minVal)
                newVal = minVal * minVal / newVal;
            else if (newVal > maxVal)
                newVal = maxVal * maxVal / newVal;
            else
                isValidVal = YES;
            } while (isValidVal == NO);

        /* check that the distribution does not go too far to the right */
        if (whichParam == 0)
            quantile95 = QuantileGamma (0.95, newVal, subVals[mp->numM10BetaCats + mp->numM10GammaCats + 7]);
        else
            quantile95 = QuantileGamma (0.95, subVals[mp->numM10BetaCats + mp->numM10GammaCats + 6], newVal);

        } while (quantile95 > 100.0);
        
    /* get proposal ratio */
    *lnProposalRatio = log(newVal / oldVal);
    
    /* get prior ratio */
    if (!strcmp(mp->m10gammapr, "Exponential"))
        (*lnPriorRatio) = (log(mp->m10gammaExp) - newVal * mp->m10gammaExp) - (log(mp->m10gammaExp) - oldVal * mp->m10gammaExp);
    else
        (*lnPriorRatio) = 0.0;
    
    /* copy new value back */
    subVals[mp->numM10BetaCats + mp->numM10GammaCats + 6 + whichParam] = newVal;
    
    /* update the omega values */
    if (DiscreteGamma (&vals[mp->numM10BetaCats], subVals[mp->numM10BetaCats + mp->numM10GammaCats + 6],
                       subVals[mp->numM10BetaCats + mp->numM10GammaCats + 7], mp->numM10GammaCats, 0) == ERROR)
        return (ERROR);
    for (i=0; i<mp->numM10GammaCats; i++)
        vals[mp->numM10BetaCats + i] += 1.0;

    /* Set update flags for all partitions that share this kappa. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for covarion, doublet, and codon models, we do want to update
       the cijk flag. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


int Move_OmegaCat (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, localNumCats, numBetaAndGammaCats;
    MrBFlt      dirichletParameters[3], *newSubVals, *oldSubVals, *newFreqs, *oldFreqs, *priorParams, sum, alpha, x, y;
    ModelParams *mp;

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* how many categories are there */
    localNumCats = 3;
    numBetaAndGammaCats = 0;
    if (!strcmp(mp->omegaVar, "M10"))
        {
        localNumCats = 2;
        numBetaAndGammaCats = mp->numM10BetaCats + mp->numM10GammaCats;
        }
        
    /* get the values we need */
    newSubVals = GetParamSubVals (param, chain, state[chain]);
    oldSubVals = GetParamSubVals (param, chain, state[chain] ^ 1);
    if (!strcmp(mp->omegaVar, "M10"))
        {
        newFreqs = newSubVals + numBetaAndGammaCats;
        oldFreqs = oldSubVals + numBetaAndGammaCats;
        priorParams = newSubVals + (numBetaAndGammaCats + 2);
        }
    else
        {
        newFreqs = newSubVals + 0;
        oldFreqs = oldSubVals + 0;
        priorParams = newFreqs + 3;
        }

    /* get parameter of proposal mechanism */
    alpha = mvp[0];

    /* multiply old values with some large number to get new values close to the old ones */
    for (i=0; i<localNumCats; i++)
        dirichletParameters[i] = oldFreqs[i] * alpha;

    /* get the new category frequencies */
    DirichletRandomVariable (dirichletParameters, newFreqs, localNumCats, seed);
    sum = 0.0;
    for (i=0; i<localNumCats; i++)
        {
        if (newFreqs[i] < 0.0001)
            newFreqs[i] = 0.0001;
        sum += newFreqs[i];
        }
    for (i=0; i<localNumCats; i++)
        newFreqs[i] /= sum;
        
    /* and get the new frequencies of the omega values, if we have another
       distribution for omega too */
    if (!strcmp(mp->omegaVar, "M10"))
        {
        for (i=0; i<mp->numM10BetaCats; i++)
            newSubVals[i] = newFreqs[0] / mp->numM10BetaCats;
        for (i=mp->numM10BetaCats; i<mp->numM10BetaCats+mp->numM10GammaCats; i++)
            newSubVals[i] = newFreqs[1] / mp->numM10GammaCats;
        }

    /* get proposal ratio */
    sum = 0.0;
    for (i=0; i<localNumCats; i++)
        sum += newFreqs[i]*alpha;
    x = LnGamma(sum);
    for (i=0; i<localNumCats; i++)
        x -= LnGamma(newFreqs[i]*alpha);
    for (i=0; i<localNumCats; i++)
        x += (newFreqs[i]*alpha-1.0)*log(oldFreqs[i]);
    sum = 0.0;
    for (i=0; i<localNumCats; i++)
        sum += oldFreqs[i]*alpha;
    y = LnGamma(sum);
    for (i=0; i<localNumCats; i++)
        y -= LnGamma(oldFreqs[i]*alpha);
    for (i=0; i<localNumCats; i++)
        y += (oldFreqs[i]*alpha-1.0)*log(newFreqs[i]);
    (*lnProposalRatio) = x - y;

    /* get prior ratio */
    x = y = 0.0;        /* ignore the gamma part, it is identical */
    for (i=0; i<localNumCats; i++)
        x += (priorParams[i]-1.0)*log(newFreqs[i]);
    for (i=0; i<localNumCats; i++)
        y += (priorParams[i]-1.0)*log(oldFreqs[i]);
    (*lnPriorRatio) = x - y;
        
    /* Set update flags for all partitions that share this omega. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_OmegaM3: Change the nonysnonymous/synonymous rate ratio
|      of one class of the M3 model
|
----------------------------------------------------------------*/
int Move_OmegaM3 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, isValidO, whichOmega;
    MrBFlt      newO, window, minO, maxO, ran, *value, *oldValue, x, y;

    /* get size of window, centered on current omega value */
    window = mvp[0];

    /* get old value of omega */
    value = GetParamVals(param, chain, state[chain]);
    oldValue = GetParamVals(param, chain, state[chain] ^ 1);
    whichOmega = (int) (RandomNumber(seed)*3.0);
    
    /* get minimum and maximum values for omega */
    if (whichOmega == 0)
        minO = 0.0;
    else
        minO = value[whichOmega-1];
    if (whichOmega == 2)
        maxO = OMEGA_MAX;
    else
        maxO = value[whichOmega+1];

    /* change value for omega */
    ran = RandomNumber(seed);
     if (maxO-minO < window)
        {
        window = maxO-minO;
        }
    newO = oldValue[whichOmega] + window * (ran - 0.5);
    
    /* check that new value is valid */
    isValidO = NO;
    do
        {
        if (newO < minO)
            newO = 2* minO - newO;
        else if (newO > maxO)
            newO = 2 * maxO - newO;
        else
            isValidO = YES;
        } while (isValidO == NO);

    /* copy new omega value back */
    value[whichOmega] = newO;

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* get prior ratio */
    x = LogOmegaPrior (value[0], value[1], value[2]);
    y = LogOmegaPrior (oldValue[0], oldValue[1], oldValue[2]);
    *lnPriorRatio = x - y;

    /* Set update flags for all partitions that share this omega. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_OmegaNeu: Change the nonysnonymous/synonymous rate ratio
|      for neutral sites
|
----------------------------------------------------------------*/
int Move_OmegaNeu (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, isOPriorExp, isValidO;
    MrBFlt      oldO, newO, window, minO, maxO, ran, *value, x, y;

    /* get size of window, centered on current omega value */
    window = mvp[0];

    /* get old value of omega */
    value = GetParamVals(param, chain, state[chain]);
    newO = oldO = value[1];

    /* get minimum and maximum values for omega */
    minO = value[0];
    maxO = value[2];
    
    /* the only way we would be updating the middle category (omega2) is
       if we have an exponential prior on all three omegas */
    isOPriorExp = YES;

    /* change value for omega */
    ran = RandomNumber(seed);
    if (maxO-minO < window)
        {
        window = maxO-minO;
        }
    newO = oldO + window * (ran - 0.5);
    
    /* check that new value is valid */
    isValidO = NO;
    do
        {
        if (newO < minO)
            newO = 2 * minO - newO;
        else if (newO > maxO)
            newO = 2 * maxO - newO;
        else
            isValidO = YES;
        } while (isValidO == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* copy new omega value back */
    value[1] = newO;

    /* get prior ratio */
    if (isOPriorExp == NO)
        {
        *lnPriorRatio = 0.0;
        }
    else
        {
        x = LogOmegaPrior (value[0], newO, value[2]);
        y = LogOmegaPrior (value[0], oldO, value[2]);
        *lnPriorRatio = x - y;
        }

    /* Set update flags for all partitions that share this omega. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_OmegaPos: Change the nonysnonymous/synonymous rate ratio
|      for positively selected sites
|
----------------------------------------------------------------*/
int Move_OmegaPos (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, isValidO, omegaUni, omegaExp1, omegaExp2;
    MrBFlt      oldO, newO, window, minO=0.0, maxO=0.0, ran, *value, x, y;
    ModelParams *mp;

    /* get size of window, centered on current omega value */
    window = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get old value of omega */
    value = GetParamVals(param, chain, state[chain]);
    newO = oldO = value[2];
    
    /* determine prior for omega */
    omegaUni = omegaExp1 = omegaExp2 = NO;
    if (param->paramId == OMEGA_BUD || param->paramId == OMEGA_BUF || param->paramId == OMEGA_FUD || param->paramId == OMEGA_FUF)
        omegaUni = YES;
    else if (param->paramId == OMEGA_BED || param->paramId == OMEGA_BEF || param->paramId == OMEGA_FED || param->paramId == OMEGA_FEF)
        omegaExp1 = YES;
    else if (param->paramId == OMEGA_ED || param->paramId == OMEGA_EF)
        omegaExp2 = YES;
        
    /* get minimum and maximum values for omega */
    if (omegaUni == YES)
        {
        minO = mp->ny98omega3Uni[0];
        if (minO < value[1])
            minO = value[1];
        maxO = mp->ny98omega3Uni[1];
        if (maxO > KAPPA_MAX)
            maxO = KAPPA_MAX;
        }
    else if (omegaExp1 == YES || omegaExp2 == YES)
        {
        minO = value[1];
        maxO = KAPPA_MAX;
        }

    /* change value for omega */
    ran = RandomNumber(seed);
    if (maxO-minO < window)
        {
        window = maxO-minO;
        }
    newO = oldO + window * (ran - 0.5);
    
    /* check that new value is valid */
    isValidO = NO;
    do
        {
        if (newO < minO)
            newO = 2* minO - newO;
        else if (newO > maxO)
            newO = 2 * maxO - newO;
        else
            isValidO = YES;
        } while (isValidO == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* copy new omega value back */
    value[2] = newO;

    /* get prior ratio (part 2) */
    if (omegaUni == YES)
        {
        *lnPriorRatio = 0.0;
        }
    else if (omegaExp1 == YES)
        {
        *lnPriorRatio = mp->ny98omega3Exp * (oldO - newO);
        }
    else if (omegaExp2 == YES)
        {
        x = LogOmegaPrior (value[0], value[1], newO);
        y = LogOmegaPrior (value[0], value[1], oldO);
        *lnPriorRatio = x - y;
        }

    /* Set update flags for all partitions that share this omega. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_OmegaPur: Change the nonysnonymous/synonymous rate ratio
|      for purifying selection sites
|
----------------------------------------------------------------*/
int Move_OmegaPur (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, isOPriorExp, isValidO;
    MrBFlt      oldO, newO, window, minO, maxO, ran, *value, x, y;

    /* get size of window, centered on current omega value */
    window = mvp[0];

    /* get old value of omega */
    value = GetParamVals(param, chain, state[chain]);
    newO = oldO = value[0];
    
    /* get minimum and maximum values for omega */
    minO = 0.0;
    maxO = value[1];
    
    /* get prior for omega */
    if (param->paramId == OMEGA_BUD || param->paramId == OMEGA_BUF || param->paramId == OMEGA_BED || 
        param->paramId == OMEGA_BEF || param->paramId == OMEGA_BFD || param->paramId == OMEGA_BFF) 
        isOPriorExp = NO;
    else
        isOPriorExp = YES;

    /* change value for omega */
    ran = RandomNumber(seed);
    if (maxO-minO < window)
        {
        window = maxO-minO;
        }
    newO = oldO + window * (ran - 0.5);
    
    /* check that new value is valid */
    isValidO = NO;
    do
        {
        if (newO < minO)
            newO = 2* minO - newO;
        else if (newO > maxO)
            newO = 2 * maxO - newO;
        else
            isValidO = YES;
        } while (isValidO == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* copy new omega value back */
    value[0] = newO;

    /* get prior ratio (part 2) */
    if (isOPriorExp == NO)
        {
        *lnPriorRatio = 0.0;
        }
    else
        {
        x = LogOmegaPrior (newO, value[1], value[2]);
        y = LogOmegaPrior (oldO, value[1], value[2]);
        *lnPriorRatio = x - y;
        }

    /* Set update flags for all partitions that share this omega. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_ParsEraser1: This proposal mechanism changes the topology and
|      branch lengths of an unrooted tree. A randomly chosen region of
|      the tree is erased. Parsimony is used to guide the selection of
|      a new topology for the erased part of the tree. The parsimony
|      branch lengths are used to guide the proposal of new branch
|      lengths. This variant (1) uses exhaustive enumeration.
|
|      Programmed by FR 2004-10-23--
|
----------------------------------------------------------------*/
int Move_ParsEraser1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, j, isVPriorExp, nSubTerminals, nEmbeddedTrees;
    MrBFlt      alphaPi, warp, minV, maxV, minP, maxP, brlensExp=0.0, newM, oldM, maxLen,
                *brlensCur, *brlensNew, *parslensCur, *parslensNew,
                curLength, newLength, lnJacobian, lnRandomRatio, alpha[2], prob[2],
                minLenCur, minLenNew, f;
    TreeNode    *p=NULL;
    Tree        *t, *subtree, *subtree1, memTree[2];
    ModelParams *mp;
    ModelInfo   *m;
    TreeInfo    tInfo;

    /* set pointers to NULL */
    subtree = subtree1 = NULL;
    brlensCur = NULL;
    for (i=0; i<2; i++)
        {
        memTree[i].allDownPass = NULL;
        memTree[i].intDownPass = NULL;
        memTree[i].nodes = NULL;
        }
    tInfo.leaf = NULL;

    /* Set alpha Pi for Dirichlet p generator */
    alphaPi = mvp[0];
    alphaPi = 0.05;
    
    /* Set the parsimony warp factor */
    warp = mvp[1];
    warp = 0.2;
    
    /* Set the number of terminals (nSubTerminals, column 3) in erased tree */
    /* Erased Nodes => Leaves => Terminals => Embedded trees => Embedded histories => New trees
                  2 => 3      => 4         => 2              => 2 = 2!             => 3 = 1*3
                  3 => 4      => 5         => 5              => 6 = 3!             => 15 = 1*3*5
                  4 => 5      => 6         => 14             => 24 = 4!            => 105 = 1*3*5*7
                  5 => 6      => 7         => 42             => 120 = 5!           => 945 = 1*3*5*7*9
                  etc               */  
    nSubTerminals = (int) (RandomNumber(seed) * 4) + 4;
    nSubTerminals = 7;

    /* initialize log prior and log proposal probabilities */
    *lnPriorRatio = *lnProposalRatio = 0.0;
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1];
        isVPriorExp = NO;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }
    minP = 3.0 * ((1.0 / 4.0) - ((1.0 / 4.0) * exp (-4.0 * minV / 3.0)));
    maxP = 3.0 * ((1.0 / 4.0) - ((1.0 / 4.0) * exp (-4.0 * maxV / 3.0)));

    /* allocate some memory for this move */
    brlensCur = (MrBFlt *) SafeMalloc (8 * nSubTerminals * sizeof (MrBFlt));
    if (!brlensCur)
        {
        MrBayesPrint ("%s   ERROR: Could not allocate brlensCur\n", spacer);
        goto errorExit;
        }
    brlensNew = brlensCur + 2*nSubTerminals;
    parslensCur = brlensCur + 4 * nSubTerminals;
    parslensNew = brlensCur + 6 * nSubTerminals;

    subtree = &memTree[0];
    subtree->nNodes = 2 * nSubTerminals - 2;
    subtree->nIntNodes = nSubTerminals - 2;
    subtree->nodes = (TreeNode *) SafeCalloc (subtree->nNodes, sizeof (TreeNode));
    subtree->allDownPass = (TreeNode **) SafeCalloc (subtree->nNodes, sizeof (TreeNode *));
    subtree->intDownPass = (TreeNode **) SafeCalloc (subtree->nIntNodes, sizeof (TreeNode *));
    if (!subtree->nodes || !subtree->intDownPass || !subtree->allDownPass)
        {
        MrBayesPrint ("%s   ERROR: Could not allocate subtree\n", spacer);
        goto errorExit;
        }

    subtree1 = &memTree[1];
    subtree1->nNodes = 2 * nSubTerminals - 2;
    subtree1->nIntNodes = nSubTerminals - 2;
    subtree1->nodes = (TreeNode *) SafeCalloc (subtree1->nNodes, sizeof (TreeNode));
    subtree1->allDownPass = (TreeNode **) SafeCalloc (subtree1->nNodes, sizeof (TreeNode *));
    subtree1->intDownPass = (TreeNode **) SafeCalloc (subtree1->nIntNodes, sizeof (TreeNode *));
    if (!subtree1->nodes || !subtree1->intDownPass || !subtree1->allDownPass)
        {
        MrBayesPrint ("%s   ERROR: Could not allocate subtree1\n", spacer);
        goto errorExit;
        }

    tInfo.leaf = (TreeNode **) SafeCalloc (t->nNodes, sizeof(TreeNode *));
    if (!tInfo.leaf)
        {
        MrBayesPrint ("%s   ERROR: Could not allocate tInfo.leaf\n", spacer);
        goto errorExit;
        }
    tInfo.vertex = tInfo.leaf + t->nNodes - t->nIntNodes;

    /* Select a random embedded subtree with nSubTerminals terminals */
    if (GetRandomEmbeddedSubtree (t, nSubTerminals, seed, &nEmbeddedTrees) == ERROR)
        {
        MrBayesPrint ("%s   ERROR: Could not get subtree\n", spacer);
        goto errorExit;
        }

    /* Set update flags (We'd better do it before the marked nodes disappear) */
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        if (p->marked == YES)
            {
            p->upDateCl = YES; 
            p->upDateTi = YES;
            }
        else if (p->left->upDateCl == YES || p->right->upDateCl == YES)
                p->upDateCl = YES; 
        }

    /* Fill in subtrees */
    CopyTreeToSubtree (t, subtree); 
    CopyTreeToSubtree (t, subtree1);

    /* Calculate downstates and upstate of root node of subtree */
    GetParsDP (t, t->root->left, chain);
    for (i=0; i<t->nIntNodes; i++)
        {
        p = t->intDownPass[i];
        if (p->marked == YES && p->anc->marked == NO)
            break;
        }
    GetParsimonySubtreeRootstate (t, p->anc, chain);

    /* Get parsimony length of current tree */
    curLength = GetParsimonyLength (subtree, chain);
    
    /* Get the Markov and parsimony branch lengths of the current subtree */
    GetParsimonyBrlens (subtree, chain, parslensCur);
    for (i=0; i<subtree->nNodes-1; i++)
        brlensCur[i] = subtree->allDownPass[i]->length;

    /* Calculate parsimony score of all trees relative to shortest tree (1.0) */
    tInfo.totalScore = 0.0;
    tInfo.stopScore = -1.0;
    tInfo.minScore = curLength;
    tInfo.warp = warp;
    ExhaustiveParsimonySearch (subtree, chain, &tInfo);
        
    /* Choose one of these trees randomly based on its score */
    tInfo.stopScore = RandomNumber(seed) * tInfo.totalScore;
    tInfo.totalScore = 0.0;
    ExhaustiveParsimonySearch (subtree1, chain, &tInfo);
    /* ShowNodes (subtree1->root, 0 , NO); */
    /* getchar(); */

    /* Get length of that tree */

    newLength = GetParsimonyLength (subtree1, chain);

    /* Get the parsimony branch lengths of the new subtree */
    GetParsimonyBrlens (subtree1, chain, parslensNew);

    /* Find the maximum length of a branch */
    maxLen = 0.0;
    for (i=0; i<t->nRelParts; i++)
        {
        j = t->relParts[i];
        m = &modelSettings[j];
        maxLen += m->numUncompressedChars;
        }
    
    /* Find the Markov branch lengths of the new subtree */
    /* Calculate Jacobian and prob ratio for the Dirichlet random number generator */
    lnJacobian = lnRandomRatio = 0.0;
    minLenCur = minLenNew = 0.0;
    for (i=0; i<subtree1->nNodes-1; i++)
        {
        minLenCur += parslensCur[i];
        minLenNew += parslensNew[i];
        }
    for (i=0; i<subtree1->nNodes-1; i++)
        {
        p = subtree1->allDownPass[i];
        f = newLength / minLenNew;
        alpha[0] = parslensNew[i] * f * alphaPi + 1.0;
        alpha[1] = (maxLen - parslensNew[i] * f) * alphaPi + 1.0;
        DirichletRandomVariable (alpha, prob, 2, seed);
        if (prob[0] >= maxP || prob[0] <= minP)
            {
            abortMove = YES;
            return NO_ERROR;
            }

        p->length = (-3.0 / 4.0) * log (1.0 - 4.0 * prob[0] / 3.0);
        lnJacobian += (-4.0 * brlensCur[i] / 3.0) - log (1.0 - 4.0 * prob[0] / 3.0);
        lnRandomRatio -= log (pow (prob[0], alpha[0] - 1.0) * pow (prob[1], alpha[1] - 1.0));
        f = curLength / minLenNew;
        alpha[0] = parslensCur[i] * f * alphaPi + 1.0;
        alpha[1] = (maxLen - parslensCur[i] * f) * alphaPi + 1.0;
        prob[0] = 3.0 * ((1.0 / 4.0) - ((1.0 / 4.0) * exp (-4.0 * brlensCur[i] / 3.0)));
        prob[1] = 1.0 - prob[0];
        lnRandomRatio += log (pow (prob[0], alpha[0] - 1.0) * pow (prob[1], alpha[1] - 1.0));
        }

    /* Store the new Markov branch lengths */
    for (i=0; i<subtree1->nNodes-1; i++)
        brlensNew[i] = subtree1->allDownPass[i]->length;

    /* Calculate the proposal ratio */
    (*lnProposalRatio) = lnJacobian + lnRandomRatio + log (warp/3.0) * (curLength - newLength) + log (1.0-warp) * (newLength - curLength);

    /* Calculate the prior ratio */
    if (isVPriorExp == YES)
        {
        newM = oldM = 0.0;
        for (i=0; i<subtree->nNodes-1; i++)
            {
            oldM += brlensCur[i];
            newM += brlensNew[i];
            }
        (*lnPriorRatio) += brlensExp * (oldM - newM);
        }

    /* Copy subtree into tree */
    CopySubtreeToTree (subtree1, t);
    /* ShowNodes (subtree1->root, 0, NO); */
    /* ShowNodes (t->root, 0, NO); */

    /* Update node sequences */
    GetDownPass (t);
    
    /* correct for difference in number of embedded subtrees */
    if (GetRandomEmbeddedSubtree (t, nSubTerminals, seed, &i) == ERROR)
        {
        MrBayesPrint ("%s   Could not count number of subtrees in Move_ParsEraser1\n", spacer);
        goto errorExit;
        }
    if (i != nEmbeddedTrees)
        (*lnProposalRatio) += log ((MrBFlt) nEmbeddedTrees / (MrBFlt) i);

    /* Free memory allocated for this move */
    free (subtree->allDownPass);
    free (subtree->intDownPass);
    free (subtree->nodes);
    free (subtree1->allDownPass);
    free (subtree1->intDownPass);
    free (subtree1->nodes);
    free (brlensCur);
    free (tInfo.leaf);

    return (NO_ERROR);

errorExit:

    free (subtree->allDownPass);
    free (subtree->intDownPass);
    free (subtree->nodes);
    free (subtree1->allDownPass);
    free (subtree1->intDownPass);
    free (subtree1->nodes);
    free (brlensCur);
    free (tInfo.leaf);

    return (ERROR);
}


int Move_ParsFossilSPRClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change branch lengths and topology (potentially) using SPR-type move, parsimony-biased */
    
    /* This move is identical to ParsSPRClock except that it uses s/n weighting and it only picks fossil subtrees. */
    
    int         i, j, n, division, n1=0, n2=0, n3=0, n4=0, n5=0, *nEvents, numMovableNodesOld, numMovableNodesNew;
    BitsLong    *pA, *pV, *pP, y[2];
    MrBFlt      x, oldBrlen=0.0, newBrlen=0.0, v1=0.0, v2=0.0, v3=0.0, v4=0.0, v5=0.0,
    v3new=0.0, lambda, **position=NULL, **rateMultiplier=NULL, *brlens,
    igrvar, *igrRate=NULL, nu, *tk02Rate=NULL, minLength=0.0, length=0.0,
    cumulativeProb, warpFactor, sum1, sum2, ran, increaseProb, decreaseProb,
    divFactor, nStates, rateMult, v_approx, minV;
    CLFlt       *nSitesOfPat, *nSites, *globalNSitesOfPat;
    TreeNode    *p, *a, *b, *u, *v, *c=NULL, *d;
    Tree        *t;
    ModelInfo   *m=NULL;
    Param       *subParm;
    
    warpFactor = mvp[0];                  /* tuning parameter determining how heavily to weight according to parsimony scores */
    increaseProb = decreaseProb = mvp[1]; /* reweighting probabilities */
    
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* get model params and model info */
    m = &modelSettings[param->relParts[0]];
    
    /* get min and max brlen in relative time and subst units */
    minV = BRLENS_MIN;
    
#   if defined (DEBUG_ParsSPRClock)
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
#   endif
    
    /* mark all nodes that only have fossil children with YES and count number movable nodes in current tree */
    numMovableNodesOld=0;
    for (i=0; i<t->nNodes-2; ++i)
    {
        p = t->allDownPass[i];
        if (p->left == NULL)
        {
            if (p->calibration == NULL)
                p->x = NO;
            else
            {
                p->x = YES;
            }
        }
        else
        {
            if (p->left->x == YES && p->right->x == YES)
            {
                p->x = YES;
            }
            else
                p->x = NO;
        }
        a = p->anc->left;
        b = p->anc->right;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL
            || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN) || p->x == NO)
            numMovableNodesOld++;
    }
    
    if (numMovableNodesOld==0)
        return (NO_ERROR);
    
    /* pick a branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes - 2))];
        a = p->anc->left;
        b = p->anc->right;
    }
    while (p->anc->isLocked == YES || p->anc->anc->anc == NULL
           || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN) || p->x == NO);
    /* skip constraints, siblings of root (and root); and consider ancestral fossils in fbd tree;
     skip all nodes that subtend extant terminals */
    
    /* set up pointers for nodes around the picked branch */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    
    /* record branch length for insertion in back move */
    if (v->length > 0.0)  /* side branch, not anc fossil */
    {
        if (v->nodeDepth > a->nodeDepth)
            oldBrlen = b->nodeDepth - v->nodeDepth - 2.0*minV;
        else
            oldBrlen = b->nodeDepth - a->nodeDepth - 2.0*minV;
    }
    v1 = a->length;
    v2 = u->length;
    v3 = v->length;
    
    /* reassign events for CPP and adjust prior and proposal ratios for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
    {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
        {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            position = subParm->position[2*chain+state[chain]];
            rateMultiplier = subParm->rateMult[2*chain+state[chain]];
            n1 = nEvents[a->index];
            n2 = nEvents[u->index];
            n3 = nEvents[v->index];
            if (n2 > 0)
            {
                position[a->index] = (MrBFlt *) SafeRealloc ((void *) position[a->index], (n1+n2) * sizeof (MrBFlt));
                rateMultiplier[a->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[a->index], (n1+n2) * sizeof (MrBFlt));
            }
            for (j=0; j<n1; j++)
                position[a->index][j] *= v1 / (v1+v2);
            for (j=n1; j<n1+n2; j++)
            {
                position[a->index][j] = (position[u->index][j-n1] * v2 + v1) / (v1+v2);
                rateMultiplier[a->index][j] = rateMultiplier[u->index][j-n1];
            }
            nEvents[a->index] = n1+n2;
            nEvents[u->index] = 0;
            if (n2 > 0)
            {
                free (position[u->index]);
                free (rateMultiplier[u->index]);
                position[u->index] = rateMultiplier[u->index] = NULL;
            }
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] += brlens[u->index];   /* only change in effective branch lengths so far */
        }   /* end CPP events parm */
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
        {
            /* adjust prior ratio */
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            if (v->length > 0.0)
                (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[v->anc->index], nu*v->length, tk02Rate[v->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[a->anc->index], nu*a->length, tk02Rate[a->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*u->length, tk02Rate[u->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*(a->length+u->length), tk02Rate[a->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = (tk02Rate[a->index] + tk02Rate[b->index]) / 2.0 * (a->length + u->length);
        }   /* end tk02 branch rate parameter */
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
        {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            
            /* adjust prior ratio for old branches */
            if (v->length > 0.0)
                (*lnPriorRatio) -= LnProbGamma(v->length/igrvar, v->length/igrvar, igrRate[v->index]);
            (*lnPriorRatio) -= LnProbGamma(a->length/igrvar, a->length/igrvar, igrRate[a->index]);
            (*lnPriorRatio) -= LnProbGamma(u->length/igrvar, u->length/igrvar, igrRate[u->index]);
            (*lnPriorRatio) += LnProbGamma((a->length+u->length)/igrvar, (a->length+u->length)/igrvar, igrRate[a->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = igrRate[a->index] * (a->length + u->length);
        }   /* end igr branch rate parameter */
    }   /* next subparameter */
    
    /* cut tree */
    a->anc = b;
    if (b->left == u)
        b->left = a;
    else
        b->right = a;
    a->length += u->length;
    a->upDateTi = YES;
    
    /* get final parsimony states for the root part */
    GetParsDP (t, t->root->left, chain);
    GetParsFP (t, t->root->left->left, chain);
    GetParsFP (t, t->root->left->right, chain);
    
    /* get downpass parsimony states for the crown part */
    GetParsDP (t, v, chain);
    
    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
    {
        p = t->allDownPass[i];
        p->marked = NO;
        p->d = 0.0;
    }
    
    /* mark nodes in the root part of the tree, first mark a */
    a->marked = YES;
    /* then move down towards root taking constraints into account */
    p = b;
    while (p->isLocked == NO && p->anc->anc != NULL)
    {
        p->marked = YES;
        p = p->anc;
    }
    /* make sure sisters of last node are marked otherwise it will not be marked in the uppass */
    p->left->marked = YES;
    p->right->marked = YES;
    /* finally move up, skip constraints and ancestral fossil */
    for (i=t->nNodes-2; i>=0; i--)
    {
        p = t->allDownPass[i];
        if (p != u && p->marked == NO && p->anc->marked == YES && p->anc->isLocked == NO
            && p->anc->nodeDepth > v->nodeDepth + minV && p->length > 0.0)
            p->marked = YES;
    }
    
    /* unmark nodes if the picked branch is 0 (ancestral fossil) */
    if (v->length < TIME_MIN)
    {
        n = 0;
        for (i=0; i<t->nNodes-1; i++)
        {
            p = t->allDownPass[i];
            if (p->nodeDepth > v->nodeDepth - minV || p->anc->nodeDepth < v->nodeDepth + minV)
                p->marked = NO;
            if (p->marked == YES)
                n++;
        }
        if (n < 2)  /* no new position to move */
        {
            abortMove = YES;
            return (NO_ERROR);
        }
    }
    
    /* find number of site patterns and modify randomly */
    globalNSitesOfPat = numSitesOfPat + ((chainId[chain] % chainParams.numChains) * numCompressedChars) + m->compCharStart;
    nSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!nSitesOfPat)
    {
        MrBayesPrint ("%s   Problem allocating nSitesOfPat in Move_ParsSPRClock\n", spacer);
        return (ERROR);
    }
    for (i=0; i<numCompressedChars; i++)
    {
        nSitesOfPat[i] = globalNSitesOfPat[i];
        for (j=0; j<globalNSitesOfPat[i]; j++)
        {
            ran = RandomNumber(seed);
            if (ran < decreaseProb)
                nSitesOfPat[i]--;
            else if (ran > 1.0 - increaseProb)
                nSitesOfPat[i]++;
        }
    }
    
    /* cycle through the possibilities and record the parsimony length */
    for (i=0; i<t->nNodes; i++)
    {
        p = t->allDownPass[i];
        if (p->marked == NO)
            continue;
        /* find the parsimony length */
        p->d = 0.0;
        for (n=0; n<t->nRelParts; n++)
        {
            division = t->relParts[n];
            
            /* Find model settings */
            m = &modelSettings[division];
            
            /* find nStates and ratemult */
            nStates = m->numModelStates;
            if (m->dataType == STANDARD)
                nStates = 2;
            rateMult = GetRate(division, chain);
            
            /* find nSitesOfPat */
            nSites = nSitesOfPat + m->compCharStart;
            
            /* find downpass parsimony sets for the node and its environment */
            pP   = m->parsSets[p->index     ];
            pA   = m->parsSets[p->anc->index];
            pV   = m->parsSets[v->index     ];
            
            length = 0.0;
            if (m->nParsIntsPerSite == 1)
            {
                for (j=0; j<m->numChars; j++)
                {
                    y[0] = (pP[j] | pA[j]) & pV[j];
                    if (y[0] == 0)
                        length += nSites[j];
                }
            }
            else /* if (m->nParsIntsPerSite == 2) */
            {
                for (j=0; j<2*m->numChars; j+=2)
                {
                    y[0] = (pP[j] | pA[j]) & pV[j];
                    y[1] = (pP[j+1] | pA[j+1]) & pV[j+1];
                    if ((y[0] | y[1]) == 0)
                        length += nSites[j/2];
                }
            }
            
            /* find nStates and v approximation using parsimony-based s/n approximation */
            nStates = m->numModelStates;
            if (m->dataType == STANDARD)
                nStates = 2;
            v_approx = length/m->numUncompressedChars + 0.0001;
            
            /* get division warp factor (prop. to prob. of change) */
            divFactor = - warpFactor * log(1.0/nStates - exp(-nStates/(nStates-1)*v_approx)/nStates);
            
            p->d += divFactor * length;
        }
    }
    
    /* find the min length and the sum for the forward move */
    minLength = -1.0;
    for (i=0; i<t->nNodes; i++)
    {
        p = t->allDownPass[i];
        if (p->marked == NO || p == a)
            continue;
        if (minLength < 0.0 || p->d < minLength)
            minLength = p->d;
    }
    sum1 = 0.0;
    for (i=0; i<t->nNodes; i++)
    {
        p = t->allDownPass[i];
        if (p->marked == YES && p != a)
            sum1 += exp (minLength - p->d);
    }
    
    /* generate a random uniform */
    ran = RandomNumber(seed) * sum1;
    
    /* select the appropriate reattachment point (not a!) */
    cumulativeProb = 0.0;
    for (i=0; i<t->nNodes; i++)
    {
        p = t->allDownPass[i];
        if (p->marked == YES && p != a)
        {
            c = p;
            cumulativeProb += exp (minLength - p->d);
            if (cumulativeProb > ran)
                break;
        }
    }
    
    /* calculate the proposal ratio */
    (*lnProposalRatio) = c->d - minLength + log(sum1);
    
    /* find the min length and the sum for the backward move */
    minLength = -1.0;
    for (i=0; i<t->nNodes; i++)
    {
        p = t->allDownPass[i];
        if (p->marked == NO || p == c)
            continue;
        if (minLength < 0.0 || p->d < minLength)
            minLength = p->d;
    }
    sum2 = 0.0;
    for (i=0; i<t->nNodes; i++)
    {
        p = t->allDownPass[i];
        if (p->marked == YES && p != c)
            sum2 += exp (minLength - p->d);
    }
    
    /* calculate the proposal ratio */
    (*lnProposalRatio) += minLength - a->d - log(sum2);
    
    /* reattach u */
    d = c->anc;
    c->anc = u;
    if (u->left == v)
        u->right = c;
    else
        u->left = c;
    u->anc = d;
    if (d->left == c)
        d->left = u;
    else
        d->right = u;
    
    if (v->length > 0.0)  /* side branch, not anc fossil */
    {
        if (c->nodeDepth > v->nodeDepth)
            newBrlen = d->nodeDepth - c->nodeDepth - 2.0*minV;
        else
            newBrlen = d->nodeDepth - v->nodeDepth - 2.0*minV;
        if (newBrlen <= 0.0)
        {
            abortMove = YES;
            free (nSitesOfPat);
            return (NO_ERROR);
        }
        
        /* adjust lengths */
        u->nodeDepth = d->nodeDepth - minV - RandomNumber(seed) * newBrlen;
        v->length = u->nodeDepth - v->nodeDepth;
        
        /* calculate proposal ratio for tree change */
        (*lnProposalRatio) += log (newBrlen / oldBrlen);
    }
    u->length = d->nodeDepth - u->nodeDepth;
    c->length = u->nodeDepth - c->nodeDepth;
    
    v3new = v->length;
    v4 = c->length;
    v5 = u->length;
    
    /* reassign events for CPP and adjust prior and proposal ratios for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
    {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
        {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            position = subParm->position[2*chain+state[chain]];
            rateMultiplier = subParm->rateMult[2*chain+state[chain]];
            for (j=0; j<nEvents[c->index]; j++)
            {
                if (position[c->index][j] > v4 / (v4+v5))
                    break;
            }
            n4 = j;
            n5 = nEvents[c->index] - j;
            nEvents[u->index] = n5;
            if (n5 > 0)
            {
                position[u->index] = (MrBFlt *) SafeRealloc ((void *) position[u->index], n5 * sizeof (MrBFlt));
                rateMultiplier[u->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[u->index], n5 * sizeof (MrBFlt));
                for (j=n4; j<nEvents[c->index]; j++)
                {
                    position[u->index][j-n4] = (position[c->index][j] * (v4+v5) - v4) / v5;
                    rateMultiplier[u->index][j-n4] = rateMultiplier[c->index][j];
                }
                if (n4 > 0)
                {
                    position[c->index] = (MrBFlt *) SafeRealloc ((void *) position[c->index], n4 * sizeof (MrBFlt));
                    rateMultiplier[c->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[c->index], n4 * sizeof (MrBFlt));
                    for (j=0; j<n4; j++)
                        position[c->index][j] *= ((v4+v5) / v4);
                }
                else
                {
                    free (position[c->index]);
                    free (rateMultiplier[c->index]);
                    position[c->index] = rateMultiplier[c->index] = NULL;
                }
                nEvents[c->index] = n4;
            }
            else
            {
                for (j=0; j<nEvents[c->index]; j++)
                    position[c->index][j] *= ((v4+v5) / v4);
            }
            
            /* adjust proposal ratio */
            (*lnProposalRatio) += n3 * log (v3new / v3);
            
            /* adjust prior ratio */
            lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);
            (*lnPriorRatio) += lambda * (v3 - v3new);
            
            /* update effective branch lengths */
            if (UpdateCppEvolLengths (subParm, a, chain) == ERROR)
            {
                abortMove = YES;
                free (nSitesOfPat);
                return (NO_ERROR);
            }
            
            if (UpdateCppEvolLengths (subParm, u, chain) == ERROR)
            {
                abortMove = YES;
                free (nSitesOfPat);
                return (NO_ERROR);
            }
        }   /* end cpp events parameter */
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
        {
            /* adjust prior ratio */
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*(c->length+u->length), tk02Rate[c->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[c->anc->index], nu*c->length, tk02Rate[c->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*u->length, tk02Rate[u->index]);
            if (v->length > 0.0)
                (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[v->anc->index], nu*v->length, tk02Rate[v->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[c->index] = c->length * (tk02Rate[c->index] + tk02Rate[c->anc->index]) / 2.0;
            brlens[v->index] = v->length * (tk02Rate[v->index] + tk02Rate[v->anc->index]) / 2.0;
            brlens[u->index] = u->length * (tk02Rate[u->index] + tk02Rate[u->anc->index]) / 2.0;
        }   /* end tk02 branch rate parameter */
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
        {
            /* adjust prior ratio */
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            
            (*lnPriorRatio) -= LnProbGamma ((c->length+u->length)/igrvar, (c->length+u->length)/igrvar, igrRate[c->index]);
            (*lnPriorRatio) += LnProbGamma (c->length/igrvar, c->length/igrvar, igrRate[c->index]);
            (*lnPriorRatio) += LnProbGamma (u->length/igrvar, u->length/igrvar, igrRate[u->index]);
            if (v->length > 0.0)
                (*lnPriorRatio) += LnProbGamma (v->length/igrvar, v->length/igrvar, igrRate[v->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[v->index] = igrRate[v->index] * v->length;
            brlens[u->index] = igrRate[u->index] * u->length;
            brlens[c->index] = igrRate[c->index] * c->length;
        }   /* end igr branch rate parameter */
    }   /* next subparameter */

    /* set tiprobs update flags */
    c->upDateTi = YES;
    u->upDateTi = YES;
    v->upDateTi = YES;
    
    /* set flags for update of cond likes down to root */
    p = u;
    while (p->anc != NULL)
    {
        p->upDateCl = YES;
        p = p->anc;
    }
    p = b;
    while (p->anc != NULL)
    {
        p->upDateCl = YES;
        p = p->anc;
    }
    
    /* get down pass sequence */
    GetDownPass (t);
    
    /* adjust prior ratio for clock tree */
    if (LogClockTreePriorRatio (param, chain, &x) == ERROR)
    {
        free (nSitesOfPat);
        return (ERROR);
    }
    (*lnPriorRatio) += x;
    
    /* adjust proposal prob for number movable nodes in new tree */
    numMovableNodesNew=0;
    for (i=0; i<t->nNodes-2; ++i)
    {
        p = t->allDownPass[i];
        if (p->left == NULL)
        {
            if (p->calibration == NULL)
                p->x = NO;
            else
            {
                p->x = YES;
            }
        }
        else
        {
            if (p->left->x == YES && p->right->x == YES)
            {
                p->x = YES;
            }
            else
                p->x = NO;
        }
        a = p->anc->left;
        b = p->anc->right;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL
            || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN) || p->x == NO)
            numMovableNodesNew++;
    }
    
    if (numMovableNodesNew!=numMovableNodesOld)
    {
        (*lnProposalRatio) += log (numMovableNodesOld / numMovableNodesNew);
    }
    
#   if defined (DEBUG_ParsSPRClock)
    ShowNodes (t->root, 2, YES);
    printf ("After\nProposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  a: %d  b: %d c: %d\n",v->index, u->index, a->index, b->index, c->index);
    getchar();
#   endif
    
    free (nSitesOfPat);
    return (NO_ERROR);
}


int Move_ParsSPR (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change topology (and branch lengths) using SPR (asymmetric) biased according to parsimony scores. */

    int         i, j, n, division, topologyHasChanged, isVPriorExp;
    BitsLong    *pA, *pV, *pP, y[2];
    MrBFlt      x, minV, maxV, brlensExp=0.0, minLength=0.0, length=0.0,
                cumulativeProb, warpFactor, ran, tuning, increaseProb, decreaseProb,
                divFactor, nStates, rateMult, v_typical, sum1, sum2, tempsum, tempc, tempy;
    CLFlt       *nSitesOfPat, *nSites, *globalNSitesOfPat;
    TreeNode    *p, *q, *a, *b, *u, *v, *c=NULL, *d;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m = NULL;

    warpFactor = mvp[0];                  /* tuning parameter determining how heavily to weight according to parsimony scores */
    increaseProb = decreaseProb = mvp[1]; /* reweighting probabilities */
    v_typical = mvp[2];                   /* typical branch length for conversion of parsimony score to log prob ratio */
    tuning = mvp[3];                      /* multiplier tuning parameter */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1] < BRLENS_MAX ? mp->brlensUni[1] : BRLENS_MAX;
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);

#   if defined (DEBUG_ParsSPR)
    // WriteTopologyToFile (stdout, t->root->left, t->isRooted);  fprintf (stdout, ";\t");
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
#   endif
    
    /* set topologyHasChanged to NO */
    topologyHasChanged = NO;
    
    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;
        p->d = 0;
        }

    /* pick a random branch */
    p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes -2))];
    q = p->anc->right;
    if (q == p)
        q = p->anc->left;
    if ((p->anc->anc->anc == NULL || p->anc->isLocked == YES) && (q->left == NULL || q->isLocked == YES))
        {
        abortMove = YES;
        return (NO_ERROR);
        }
        
    /* set up pointers for nodes around the picked branch */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;

    /* clip tree */
    a->anc = b;
    if (b->left == u)
        b->left = a;
    else
        b->right = a;

    /* get final parsimony states for the root part */
    GetParsDP (t, t->root->left, chain);
    GetParsFP (t, t->root->left, chain);

    /* get downpass parsimony states for the crown part */
    GetParsDP (t, v, chain);

    /* mark all nodes in the root part of the tree, taking constraints into account */
    /* first mark a */
    a->marked = YES;
    /* then move down towards root */
    if (u->isLocked == NO)
        {
        p = a->anc;
        while (p->anc != NULL)
            {
            p->marked = YES;
            if (p->isLocked == YES)
                break;
            p = p->anc;
            }
        }

    /* finally move up */
    for (i=t->nNodes-2; i>=0; i--)
        {
        p = t->allDownPass[i];
        if (p->marked == NO && p->anc->marked == YES && p->anc->isLocked == NO && p != u)
            p->marked = YES;
        }

    /* find number of site patterns and modify randomly */
    globalNSitesOfPat = numSitesOfPat + (chainId[chain] % chainParams.numChains) * numCompressedChars + m->compCharStart;
    nSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!nSitesOfPat)
        {
        MrBayesPrint ("%s   Problem allocating nSitesOfPat in Move_ParsSPR\n", spacer);
        return (ERROR);
        }
    for (i=0; i<numCompressedChars; i++)
        {
        nSitesOfPat[i] = globalNSitesOfPat[i];
        for (j=0; j<globalNSitesOfPat[i]; j++)
            {
            ran = RandomNumber(seed);
            if (ran < decreaseProb)
                nSitesOfPat[i]--;
            else if (ran > 1.0 - increaseProb)
                nSitesOfPat[i]++;
            }
        }

    /* cycle through the possibilities and record the parsimony length */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO)
            continue;
        /* find the parsimony length */
        p->d = 0.0;
        for (n=0; n<t->nRelParts; n++)
            {
            division = t->relParts[n];
            
            /* Find model settings */
            m = &modelSettings[division];

            /* find nSitesOfPat */
            nSites = nSitesOfPat + m->compCharStart;

            /* find downpass parsimony sets for the node and its environment */
            pP   = m->parsSets[p->index];
            pA   = m->parsSets[p->anc->index];
            pV   = m->parsSets[v->index];
        
            length = 0.0;
            if (m->nParsIntsPerSite == 1)
                {
                for (j=0; j<m->numChars; j++)
                    {
                    y[0] = (pP[j] | pA[j]) & pV[j];
                    if (y[0] == 0)
                        length += nSites[j];
                    }
                }
            else /* if (m->nParsIntsPerSite == 2) */
                {
                for (j=0; j<2*m->numChars; j+=2)
                    {
                    y[0] = (pP[j] | pA[j]) & pV[j];
                    y[1] = (pP[j+1] | pA[j+1]) & pV[j+1];
                    if ((y[0] | y[1]) == 0)
                        length += nSites[j/2];
                    }
                }

            /* find nStates and ratemult */
            nStates = m->numModelStates;
            if (m->dataType == STANDARD)
                nStates = 2;
            rateMult = GetRate(division, chain);

            /* get division warp factor */
            divFactor = - warpFactor * log((1.0/nStates) - exp(-nStates/(nStates-1)*v_typical*rateMult)/nStates);

            p->d += divFactor * length;
            }
        }

    /* find the min length and the sum for the forward move */    
    minLength = -1.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO || p == a)
            continue;
        if (minLength < 0.0 || p->d < minLength)
            minLength = p->d;
        }
    sum1 = 0.0; tempc = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES && p != a)
            {
            /* Kahan summation to reduce numerical error */
            tempy = exp (minLength - p->d) - tempc;
            tempsum = sum1 + tempy;
            tempc = (tempsum - sum1) - tempy;
            sum1 = tempsum;
            }
        }
    
    /* generate a random uniform */
    ran = RandomNumber(seed) * sum1;

    /* select the appropriate reattachment point */
    cumulativeProb = 0.0; tempc = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES && p != a)
            {
            c = p;
            /* Kahan summation to reduce numerical error */
            tempy = exp (minLength - p->d) - tempc;
            tempsum = cumulativeProb + tempy;
            tempc = (tempsum - cumulativeProb) - tempy;
            cumulativeProb = tempsum;
            if (cumulativeProb > ran)
                break;
            }
        }

    /* calculate the proposal ratio */
    (*lnProposalRatio) = c->d - minLength + log(sum1);

    /* find the min length and the sum for the backward move */
    minLength = -1.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO || p == c)
            continue;
        if (minLength < 0.0 || p->d < minLength)
            minLength = p->d;
        }
    sum2 = 0.0; tempc = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES && p != c)
            {
            /* Kahan summation to reduce numerical error */
            tempy = exp (minLength - p->d) - tempc;
            tempsum = sum2 + tempy;
            tempc = (tempsum - sum2) - tempy;
            sum2 = tempsum;
            }
        }

    /* calculate the proposal ratio */
    (*lnProposalRatio) += minLength - a->d - log(sum2);

    /* reattach */
    d = c->anc;
    c->anc = u;
    if (u->left == v)
        u->right = c;
    else
        u->left = c;
    if (d->left == c)
        d->left = u;
    else
        d->right = u;
    u->anc = d;

    /* c cannot be a, as a is skiped in the selection for reattachment point */
    assert (c != a);
    /* transfer lock if necessary */
    /* if u is locked, then we have moved upwards and need to leave the u lock behind */
    if (u->isLocked == YES)
        {
        u->isLocked = NO;
        a->isLocked = YES;
        a->lockID = u->lockID;
        u->lockID = -1;
        }
    /* if c is on root path and locked, we need to transfer lock to u */
    p = b;
    while (p != NULL)
        {
        if (p == c)
            break;
        p = p->anc;
        }
    if (p == c && c->isLocked == YES)
        {
        u->isLocked = YES;
        u->lockID = c->lockID;
        c->isLocked = NO;
        c->lockID = -1;
        }

    /* reassign branch lengths */
    p = c;
    while (p->anc->anc != NULL)
        {
        if (p == a)
            break;
        p = p->anc;
        }
    if (p == a)
        {
        /* c is descendant to a so move a->length and not u->length */
        x = u->length;
        u->length = a->length;
        a->length = x;
        }
    p = a;
    while (p->anc->anc != NULL)
        {
        if (p == c)
            break;
        p = p->anc;
        }
    if (p == c)
        {
        /* c is ancestor to a so insert above instead of below */
        x = c->length;
        c->length = u->length;
        u->length = x;
        }

    topologyHasChanged = YES;

    /* hit c length with multiplier (a and u dealt with below) */
    x = c->length * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV)
            x = minV * minV / x;
        else if (x > maxV)
            x = maxV * maxV / x;
        }
    /* calculate proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / c->length);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (c->length - x);
    c->length = x;
    
    /* hit a length with multiplier */
    x = a->length * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV)
            x = minV * minV / x;
        else if (x > maxV)
            x = maxV * maxV / x;
        }
    /* calculate proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / a->length);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (a->length - x);
    a->length = x;

    /* hit u length with multiplier */
    x = u->length * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV)
            x = minV * minV / x;
        else if (x > maxV)
            x = maxV * maxV / x;
        }
    /* calculate proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / u->length);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (u->length - x);
    u->length = x;

    /* set tiprobs update flags */
    a->upDateTi = YES;
    u->upDateTi = YES;
    c->upDateTi = YES;  /* could be same as a but that does not matter */

    /* set flags for update of cond likes from u and down to root */
    p = u;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }

    /* set flags for update of cond likes from b and down to root */
    p = b;
    while (p->anc != NULL && p->upDateCl == NO)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);
    
#   if defined (DEBUG_ParsSPR)
    // WriteTopologyToFile (stdout, t->root->left, t->isRooted);
    // fprintf (stdout, ";\t");  fprintf (stdout, "%lf\n", *lnProposalRatio);
    printf ("After:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
    printf ("Proposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  a: %d  b: %d\n",v->index, u->index, a->index, b->index);
    printf ("Has topology changed? %d\n",topologyHasChanged);
    getchar();
#   endif

    free (nSitesOfPat);

    return (NO_ERROR);
}


int Move_ParsSPR1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change topology (and branch lengths) using SPR (symmetric) biased according to parsimony scores,
       controlled by a window defined by a certain node distance radius. Note: w = e^{-S} */
    
    int         i, j, k, n, division, topologyHasChanged, moveInRoot, nNeighbor, nRoot, nCrown, iA, jC, isVPriorExp;
    BitsLong    *pA, *pB, *pP, *pC, *pD, y[2];
    MrBFlt      x, minV, maxV, brlensExp=0.0, minLength=0.0, length=0.0, *parLength=NULL, prob, ran, tuning, warpFactor,
                sum1, sum2, tempsum, tempc, tempy;
    CLFlt       *nSites, *nSitesOfPat=NULL, *globalNSitesOfPat;
    TreeNode    *p, *q, *r, *a, *b, *u, *v, *c, *d, *newB, *newA, *newC, **pRoot=NULL, **pCrown=NULL;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m=NULL;
    
    warpFactor = mvp[0];                  /* tuning parameter determining how heavily to weight according to parsimony scores */
//  increaseProb = decreaseProb = mvp[1]; /* reweighting probabilities */
//  v_typical = mvp[2];                   /* typical branch length for conversion of parsimony score to log prob ratio */
    tuning = mvp[3];                      /* multiplier tuning parameter */
    nNeighbor = (int)mvp[4];              /* distance to move picked branch in root and crown part */
    
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;
    
    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1] < BRLENS_MAX ? mp->brlensUni[1] : BRLENS_MAX;
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);
    
    /* set topologyHasChanged to NO */
    topologyHasChanged = NO;
    
    /* pick a random branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes -1))];
        q = p->anc->right;  if (q == p) q = p->anc->left;
        i = j = 0;
        if (p->left == NULL)
            j = 2;
        if (p->anc->anc == NULL)
            i = 2;
        if (p->anc->anc != NULL && (p->anc->isLocked == YES || p->anc->anc->anc == NULL))
            i++;
        if (p->anc->anc != NULL && (q->isLocked == YES || q->left == NULL))
            i++;
        if (p->left != NULL && (p->left->isLocked == YES || p->left->left == NULL))
            j++;
        if (p->left != NULL && (p->right->isLocked == YES || p->right->left == NULL))
            j++;
        } while (i == 2 && j == 2);
    
    /* determine moving direction */
    if (j == 2)
        moveInRoot = YES;
    else if (i == 2)
        moveInRoot = NO;
    else if (RandomNumber(seed) < 0.5)
        moveInRoot = YES;
    else
        moveInRoot = NO;

    /* set up pointers for nodes around the picked branch */
    /* should never change u, v, a, b, c, d pointers */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    c = v->left;
    d = v->right;

    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->x = 0;
        p->marked = NO;
        }

    /* mark nodes nNeighbor away in root (negative) or crown (positive) respecting constraints */
    nRoot = nCrown = 1;
    if (moveInRoot == YES)
        {
        /* clip root part of tree */
        a->anc = b;
        if (b->left == u)
            b->left = a;
        else
            b->right = a;
    
        /* mark the root part */
        if (u->isLocked == NO )
            {
            p = a; q = b; n = 0;
            while (q->anc != NULL)
                {
                q->marked = YES;
                q->x = n;    // temporary, for MarkDistance below
                if (q->left == p)
                    MarkDistance(q->right, YES, nNeighbor, &nRoot);
                else
                    MarkDistance(q->left,  YES, nNeighbor, &nRoot);
                q->x = --n;  // final
                nRoot++;
                if (q->isLocked == YES || abs(q->x) >= nNeighbor)
                    break;
                p = q; q = q->anc;
                }
            }
        if (a->isLocked == NO)
            {
            MarkDistance(a->left,  YES, nNeighbor, &nRoot);
            MarkDistance(a->right, YES, nNeighbor, &nRoot);
            }
            
        /* get final parsimony states for the root part */
        GetParsDP (t, t->root->left, chain);
        GetParsFP (t, t->root->left, chain);
        /* get final parsimony states for the crown part */
        GetParsDP (t, v, chain);
        GetParsFP (t, v, chain);
        }
    else  /* moveInRoot == NO */
        {
        /* clip crown part of tree */
        c->anc = d;
        d->anc = c;

        /* mark the crown part */
        if (c->isLocked == NO)
            {
            MarkDistance(c->left,  NO, nNeighbor, &nCrown);
            MarkDistance(c->right, NO, nNeighbor, &nCrown);
            }
        if (d->isLocked == NO)
            {
            MarkDistance(d->left,  NO, nNeighbor, &nCrown);
            MarkDistance(d->right, NO, nNeighbor, &nCrown);
            }
        
        /* get final parsimony states for the root part */
        if (u->anc != NULL) {
            a->anc = b;  /* clip */
            if (b->left == u) b->left = a;
            else             b->right = a;
            GetParsDP (t, t->root->left, chain);
            GetParsFP (t, t->root->left, chain);
            a->anc = u;  /* change back */
            if (b->left == a) b->left = u;
            else             b->right = u;
            }
        /* get final parsimony states for the crown part */
        GetParsDP (t, c, chain);
        GetParsDP (t, d, chain);
        GetParsFP (t, c, chain);
        GetParsFP (t, d, chain);
        }
    
    /* find number of site patterns and modify randomly */
    globalNSitesOfPat = numSitesOfPat + (chainId[chain] % chainParams.numChains) * numCompressedChars + m->compCharStart;
    nSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!nSitesOfPat)  goto errorExit;
    for (i=0; i<numCompressedChars; i++)
        {
        nSitesOfPat[i] = globalNSitesOfPat[i];
    /*  for (j=0; j<globalNSitesOfPat[i]; j++)
            {
            ran = RandomNumber(seed);
            if (ran < decreaseProb)
                nSitesOfPat[i]--;
            else if (ran > 1.0 - increaseProb)
                nSitesOfPat[i]++;
            }  // this is slow at the moment */
        }

    /* need to alloc a matrix for parsimony lengths, an array of pointers to crown part,
       and an array of pointers to root part. */
    parLength = (MrBFlt *) SafeCalloc ((size_t)nRoot * (size_t)nCrown, sizeof(MrBFlt));
    pRoot  = (TreeNode **) SafeCalloc(nRoot,  sizeof(TreeNode *));
    pCrown = (TreeNode **) SafeCalloc(nCrown, sizeof(TreeNode *));
    if (!parLength || !pRoot || !pCrown)  goto errorExit;
    
    /* starting position */
    pRoot[0] = a; pCrown[0] = c;
    for (i=j=1, n=t->nNodes-2; n>=0; n--)
        {  /* and the rest */
        p = t->allDownPass[n];
        if (p->marked == YES && p->x < 0)
            pRoot[i++] = p;
        if (p->marked == YES && p->x > 0)
            pCrown[j++] = p;
        }
    assert (i==nRoot && j==nCrown);
    
    /* cycle through the possibilities and record the parsimony length */
    for (j=0; j<nCrown; j++)
        {
        for (i=0; i<nRoot; i++)
            {
            parLength[i+j*nRoot] = 0.0;
            for (n=0; n<t->nRelParts; n++)
                {
                division = t->relParts[n];
                
                /* Find model settings */
                m = &modelSettings[division];
                
                /* find nSitesOfPat */
                nSites = nSitesOfPat + m->compCharStart;
                
                /* find parsimony length for each candidate position */
                length = 0.0;
                if (moveInRoot == YES)
                    {
                    pA = m->parsSets[pRoot[i]->index];
                    pB = m->parsSets[pRoot[i]->anc->index];
                    pP = m->parsSets[v->index];
                    
                    if (m->nParsIntsPerSite == 1)
                        {
                        for (k=0; k<m->numChars; k++)
                            {
                            y[0] = (pA[k] | pB[k]) & pP[k];
                            if (y[0] == 0)
                                length += nSites[k];
                            }
                        }
                    else /* if (m->nParsIntsPerSite == 2) */
                        {
                        for (k=0; k<2*m->numChars; k+=2)
                            {
                            y[0] = (pA[k] | pB[k]) & pP[k];
                            y[1] = (pA[k+1] | pB[k+1]) & pP[k+1];
                            if ((y[0] | y[1]) == 0)
                                length += nSites[k/2];
                            }
                        }
                    }
                else if (u->anc == NULL)
                    {
                    pP = m->parsSets[u->index];
                    pC = m->parsSets[pCrown[j]->index];
                    pD = m->parsSets[pCrown[j]->anc->index];
                    
                    if (m->nParsIntsPerSite == 1)
                        {
                        for (k=0; k<m->numChars; k++)
                            {
                            y[0] = pP[k] & (pC[k] | pD[k]);
                            if (y[0] == 0)
                                length += nSites[k];
                            }
                        }
                    else /* if (m->nParsIntsPerSite == 2) */
                        {
                        for (k=0; k<2*m->numChars; k+=2)
                            {
                            y[0] = pP[k] & (pC[k] | pD[k]);
                            y[1] = pP[k+1] & (pC[k+1] | pD[k+1]);
                            if ((y[0] | y[1]) == 0)
                                length += nSites[k/2];
                            }
                        }
                    }
                else
                    {
                    pA = m->parsSets[a->index];
                    pB = m->parsSets[b->index];
                    pC = m->parsSets[pCrown[j]->index];
                    pD = m->parsSets[pCrown[j]->anc->index];
                    
                    if (m->nParsIntsPerSite == 1)
                        {
                        for (k=0; k<m->numChars; k++)
                            {
                            y[0] = (pA[k] | pB[k]) & (pC[k] | pD[k]);
                            if (y[0] == 0)
                                length += nSites[k];
                            }
                        }
                    else /* if (m->nParsIntsPerSite == 2) */
                        {
                        for (k=0; k<2*m->numChars; k+=2)
                            {
                            y[0] = (pA[k] | pB[k]) & (pC[k] | pD[k]);
                            y[1] = (pA[k+1] | pB[k+1]) & (pC[k+1] | pD[k+1]);
                            if ((y[0] | y[1]) == 0)
                                length += nSites[k/2];
                            }
                        }
                    }
                
                /* get division warp factor */
                parLength[i+j*nRoot] += warpFactor * length;
                }
            }
        }
    
    /* find the min length and the sum for the forward move */
    minLength = -1.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            if (minLength > parLength[i+j*nRoot] || minLength < 0.0)
                minLength = parLength[i+j*nRoot];
            }
    sum1 = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = sum1 + tempy;  tempc = (tempsum - sum1) - tempy;
            sum1 = tempsum;
            // sum1 += exp(minLength - parLength[i+j*nRoot]);
            // printf("%d %d %lf\n", i, j, exp(minLength - parLength[i+j*nRoot]));
            }

    /* generate a random uniform */
    ran = RandomNumber(seed) * sum1;
    
    /* select the appropriate reattachment point */
    newA = a; newC = c;
    prob = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            // prob += exp (minLength - parLength[i+j*nRoot]);
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = prob + tempy;  tempc = (tempsum - prob) - tempy;
            prob = tempsum;
            if (prob > ran) {
                /* proposed new attaching position */
                newA = pRoot[i];
                newC = pCrown[j];
                goto outLoop;
                }
            }
outLoop:;
    iA = i; jC = j;
    
    /* calculate the proposal ratio */
    (*lnProposalRatio) = parLength[i+j*nRoot] - minLength + log(sum1);
    
    /* find the min length and the sum for the backward move */
    minLength = -1.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == iA && j == jC)  // exclude new position
                continue;
            if (minLength > parLength[i+j*nRoot] || minLength < 0.0)
                minLength = parLength[i+j*nRoot];
            }
    sum2 = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == iA && j == jC)  // exclude new position
                continue;
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = sum2 + tempy;  tempc = (tempsum - sum2) - tempy;
            sum2 = tempsum;
            // sum2 += exp (minLength - parLength[i+j*nRoot]);
            // printf("%d %d %lf\n", i, j, exp(minLength - parLength[i+j*nRoot]));
            }
    
    /* calculate the proposal ratio */
    (*lnProposalRatio) += minLength - parLength[0] - log(sum2);
    
    if (moveInRoot == YES)  /* root part has changed */
        {
        /* reattach the root part */
        newB = newA->anc;
        newA->anc = u;
        if (u->left == v)
            u->right = newA;
        else
            u->left = newA;
        u->anc = newB;
        if (newB->left == newA)
            newB->left = u;
        else
            newB->right = u;

        /* if u is locked, then we have moved upwards and need to leave the u lock behind */
        if (u->isLocked == YES)
            {
            u->isLocked = NO;
            a->isLocked = YES;
            a->lockID = u->lockID;
            u->lockID = -1;
            }
        
        p = newA;
        while (p->anc != NULL)
            {
            if (p == a) break;
            p = p->anc;
            }
        if (p == a)
            {
            /* newA is descendant to a so move a->length not u->length */
            x = u->length;
            u->length = a->length;
            a->length = x;
            }

        p = b;
        while (p->anc != NULL)
            {
            if (p == newA) break;
            p = p->anc;
            }
        if (p == newA)
            {
            /* newA is ancestor to a so insert above instead of below */
            x = newA->length;
            newA->length = u->length;
            u->length = x;
            /* newA is on root path and locked, we need to transfer lock to u */
            if (newA->isLocked == YES) {
                u->isLocked = YES;
                u->lockID = newA->lockID;
                newA->isLocked = NO;
                newA->lockID = -1;
                }
            }
        
        /* hit a length with multiplier */
        x = a->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / a->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (a->length - x);
        a->length = x;

        /* hit u length with multiplier */
        x = u->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / u->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (u->length - x);
        u->length = x;

        /* hit newA length with multiplier */
        x = newA->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / newA->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (newA->length - x);
        newA->length = x;
         
        /* set tiprobs update flags */
        newA->upDateTi = YES;
        a->upDateTi = YES;
        u->upDateTi = YES;
            
        /* set flags for update of cond likes */
        p = u;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        p = b;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        }
    
    if (moveInRoot == NO)  /* crown part has changed */
        {
        r = newC;
        q = newB = newC->anc;
        /* rotate nodes from newC to c or d (whichever is closest) */
        tempc = r->length;
        while (r != c && r != d)
            {
            p = q->anc;
            /* rotate pointers of q */
            if (q->left == r)
                q->left = p;
            else
                q->right = p;
            q->anc = r;
            /* swap q and old */
            tempy = q->length;
            q->length = tempc;
            q->upDateTi = YES;
            tempc = tempy;
            /* make sure we get q and r initialized for next round */
            r = q;
            q = p;
            }
        newB->length = tempc;
            
        /* hit q length with multiplier while we are at it */
        x = q->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }        
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / q->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (q->length - x);
        q->length = x;
        q->upDateTi = YES;

        /* hit newB length with multiplier */
        x = newB->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / newB->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (newB->length - x);
        newB->length = x;
        newB->upDateTi = YES;

        /* hit newC length with multiplier */
        x = newC->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / newC->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (newC->length - x);
        newC->length = x;
        newC->upDateTi = YES;
 
        /* reattach the crown part */
        v->left = newC;
        v->right = newB;
        newC->anc = newB->anc = v;
        
        /* set flags for update of cond likes */
        p = newC;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        p = r;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        }
    
    topologyHasChanged = YES;

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }
    
    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);
    
    /* free up local memory */
    free (parLength); free (pRoot); free (pCrown); free (nSitesOfPat);
    
    return (NO_ERROR);
    
errorExit:
    MrBayesPrint ("%s   Problem allocating memory in Move_ParsSPR\n", spacer);
    free (parLength); free (pRoot); free (pCrown); free (nSitesOfPat);
    
    return (ERROR);
}


int Move_ParsSPR2 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change topology (and branch lengths) using SPR (symmetric) biased according to parsimony scores,
       controlled by a window defined by a certain node distance radius. Note: S/N */
    
    int         i, j, k, n, division, topologyHasChanged, moveInRoot, nNeighbor, nRoot, nCrown, iA, jC, isVPriorExp;
    BitsLong    *pA, *pB, *pP, *pC, *pD, y[2];
    MrBFlt      x, minV, maxV, brlensExp=0.0, minLength=0.0, length=0.0, *parLength=NULL, prob, ran, tuning, warpFactor,
                v_typical, divFactor, nStates, sum1, sum2, tempsum, tempc, tempy;
    CLFlt       *nSites, *nSitesOfPat=NULL, *globalNSitesOfPat;
    TreeNode    *p, *q, *r, *a, *b, *u, *v, *c, *d, *newB, *newA, *newC, **pRoot=NULL, **pCrown=NULL;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m=NULL;
    
    warpFactor = mvp[0];                  /* tuning parameter determining how heavily to weight according to parsimony scores */
//  increaseProb = decreaseProb = mvp[1]; /* reweighting probabilities */
    v_typical = mvp[2];                   /* typical branch length for conversion of parsimony score to log prob ratio */
    tuning = mvp[3];                      /* multiplier tuning parameter */
    nNeighbor = (int)mvp[4];              /* distance to move picked branch in root and crown part */
    
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;
    
    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1] < BRLENS_MAX ? mp->brlensUni[1] : BRLENS_MAX;
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);
    
    /* set topologyHasChanged to NO */
    topologyHasChanged = NO;
    
    /* pick a random branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes -1))];
        q = p->anc->right;  if (q == p) q = p->anc->left;
        i = j = 0;
        if (p->left == NULL)
            j = 2;
        if (p->anc->anc == NULL)
            i = 2;
        if (p->anc->anc != NULL && (p->anc->isLocked == YES || p->anc->anc->anc == NULL))
            i++;
        if (p->anc->anc != NULL && (q->isLocked == YES || q->left == NULL))
            i++;
        if (p->left != NULL && (p->left->isLocked == YES || p->left->left == NULL))
            j++;
        if (p->left != NULL && (p->right->isLocked == YES || p->right->left == NULL))
            j++;
        } while (i == 2 && j == 2);
    
    /* pick an internal branch
    do  {
        p = t->intDownPass[(int)(RandomNumber(seed)*(t->nIntNodes-1))];
        q = p->anc->left;  if (q == p)  q = p->anc->right;
        i = j = 0;
        if (q->isLocked == YES || q->left == NULL)
            i++;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL)
            i++;
        if (p->left->isLocked == YES || p->left->left == NULL)
            j++;
        if (p->right->isLocked == YES || p->right->left == NULL)
            j++;
        } while (i == 2 && j == 2);
    */

    /* determine moving direction */
    if (j == 2)
        moveInRoot = YES;
    else if (i == 2)
        moveInRoot = NO;
    else if (RandomNumber(seed) < 0.5)
        moveInRoot = YES;
    else
        moveInRoot = NO;

    /* set up pointers for nodes around the picked branch */
    /* should never change u, v, a, b, c, d pointers */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;
    c = v->left;
    d = v->right;

    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->x = 0;
        p->marked = NO;
        }

    /* mark nodes nNeighbor away in root (negative) or crown (positive) respecting constraints */
    nRoot = nCrown = 1;
    if (moveInRoot == YES)
        {
        /* clip root part of tree */
        a->anc = b;
        if (b->left == u)
            b->left = a;
        else
            b->right = a;
    
        /* mark the root part */
        if (u->isLocked == NO )
            {
            p = a; q = b; n = 0;
            while (q->anc != NULL)
                {
                q->marked = YES;
                q->x = n;    // temporary, for MarkDistance below
                if (q->left == p)
                    MarkDistance(q->right, YES, nNeighbor, &nRoot);
                else
                    MarkDistance(q->left,  YES, nNeighbor, &nRoot);
                q->x = --n;  // final
                nRoot++;
                if (q->isLocked == YES || abs(q->x) >= nNeighbor)
                    break;
                p = q; q = q->anc;
                }
            }
        if (a->isLocked == NO)
            {
            MarkDistance(a->left,  YES, nNeighbor, &nRoot);
            MarkDistance(a->right, YES, nNeighbor, &nRoot);
            }
            
        /* get final parsimony states for the root part */
        GetParsDP (t, t->root->left, chain);
        GetParsFP (t, t->root->left, chain);
        /* get final parsimony states for the crown part */
        GetParsDP (t, v, chain);
        GetParsFP (t, v, chain);
        }
    else  /* moveInRoot == NO */
        {
        /* clip crown part of tree */
        c->anc = d;
        d->anc = c;

        /* mark the crown part */
        if (c->isLocked == NO)
            {
            MarkDistance(c->left,  NO, nNeighbor, &nCrown);
            MarkDistance(c->right, NO, nNeighbor, &nCrown);
            }
        if (d->isLocked == NO)
            {
            MarkDistance(d->left,  NO, nNeighbor, &nCrown);
            MarkDistance(d->right, NO, nNeighbor, &nCrown);
            }
        
        /* get final parsimony states for the root part */
        if (u->anc != NULL) {
            a->anc = b;  /* clip */
            if (b->left == u) b->left = a;
            else             b->right = a;
            GetParsDP (t, t->root->left, chain);
            GetParsFP (t, t->root->left, chain);
            a->anc = u;  /* change back */
            if (b->left == a) b->left = u;
            else             b->right = u;
            }
        /* get final parsimony states for the crown part */
        GetParsDP (t, c, chain);
        GetParsDP (t, d, chain);
        GetParsFP (t, c, chain);
        GetParsFP (t, d, chain);
        }
    
    /* find number of site patterns and modify randomly */
    globalNSitesOfPat = numSitesOfPat + (chainId[chain] % chainParams.numChains) * numCompressedChars + m->compCharStart;
    nSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!nSitesOfPat)  goto errorExit;
    for (i=0; i<numCompressedChars; i++)
        {
        nSitesOfPat[i] = globalNSitesOfPat[i];
    /*  for (j=0; j<globalNSitesOfPat[i]; j++)
            {
            ran = RandomNumber(seed);
            if (ran < decreaseProb)
                nSitesOfPat[i]--;
            else if (ran > 1.0 - increaseProb)
                nSitesOfPat[i]++;
            }  // this is slow at the moment */
        }

    /* need to alloc a matrix for parsimony lengths, an array of pointers to crown part,
       and an array of pointers to root part. */
    parLength = (MrBFlt *) SafeCalloc ((size_t)nRoot * (size_t)nCrown, sizeof(MrBFlt));
    pRoot  = (TreeNode **) SafeCalloc(nRoot,  sizeof(TreeNode *));
    pCrown = (TreeNode **) SafeCalloc(nCrown, sizeof(TreeNode *));
    if (!parLength || !pRoot || !pCrown)  goto errorExit;
    
    /* starting position */
    pRoot[0] = a; pCrown[0] = c;
    for (i=j=1, n=t->nNodes-2; n>=0; n--)
        {  /* and the rest */
        p = t->allDownPass[n];
        if (p->marked == YES && p->x < 0)
            pRoot[i++] = p;
        if (p->marked == YES && p->x > 0)
            pCrown[j++] = p;
        }
    assert (i==nRoot && j==nCrown);
    
    /* cycle through the possibilities and record the parsimony length */
    for (j=0; j<nCrown; j++)
        {
        for (i=0; i<nRoot; i++)
            {
            parLength[i+j*nRoot] = 0.0;
            for (n=0; n<t->nRelParts; n++)
                {
                division = t->relParts[n];
                
                /* Find model settings */
                m = &modelSettings[division];
                
                /* find nSitesOfPat */
                nSites = nSitesOfPat + m->compCharStart;
                
                /* find parsimony length for each candidate position */
                length = 0.0;
                if (moveInRoot == YES)
                    {
                    pA = m->parsSets[pRoot[i]->index];
                    pB = m->parsSets[pRoot[i]->anc->index];
                    pP = m->parsSets[v->index];
                    
                    if (m->nParsIntsPerSite == 1)
                        {
                        for (k=0; k<m->numChars; k++)
                            {
                            y[0] = (pA[k] | pB[k]) & pP[k];
                            if (y[0] == 0)
                                length += nSites[k];
                            }
                        }
                    else /* if (m->nParsIntsPerSite == 2) */
                        {
                        for (k=0; k<2*m->numChars; k+=2)
                            {
                            y[0] = (pA[k] | pB[k]) & pP[k];
                            y[1] = (pA[k+1] | pB[k+1]) & pP[k+1];
                            if ((y[0] | y[1]) == 0)
                                length += nSites[k/2];
                            }
                        }
                    }
                else if (u->anc == NULL)
                    {
                    pP = m->parsSets[u->index];
                    pC = m->parsSets[pCrown[j]->index];
                    pD = m->parsSets[pCrown[j]->anc->index];
                    
                    if (m->nParsIntsPerSite == 1)
                        {
                        for (k=0; k<m->numChars; k++)
                            {
                            y[0] = pP[k] & (pC[k] | pD[k]);
                            if (y[0] == 0)
                                length += nSites[k];
                            }
                        }
                    else /* if (m->nParsIntsPerSite == 2) */
                        {
                        for (k=0; k<2*m->numChars; k+=2)
                            {
                            y[0] = pP[k] & (pC[k] | pD[k]);
                            y[1] = pP[k+1] & (pC[k+1] | pD[k+1]);
                            if ((y[0] | y[1]) == 0)
                                length += nSites[k/2];
                            }
                        }
                    }
                else
                    {
                    pA = m->parsSets[a->index];
                    pB = m->parsSets[b->index];
                    pC = m->parsSets[pCrown[j]->index];
                    pD = m->parsSets[pCrown[j]->anc->index];
                    
                    if (m->nParsIntsPerSite == 1)
                        {
                        for (k=0; k<m->numChars; k++)
                            {
                            y[0] = (pA[k] | pB[k]) & (pC[k] | pD[k]);
                            if (y[0] == 0)
                                length += nSites[k];
                            }
                        }
                    else /* if (m->nParsIntsPerSite == 2) */
                        {
                        for (k=0; k<2*m->numChars; k+=2)
                            {
                            y[0] = (pA[k] | pB[k]) & (pC[k] | pD[k]);
                            y[1] = (pA[k+1] | pB[k+1]) & (pC[k+1] | pD[k+1]);
                            if ((y[0] | y[1]) == 0)
                                length += nSites[k/2];
                            }
                        }
                    }
                
                /* find nStates and ratemult */
                nStates = m->numModelStates;
                if (m->dataType == STANDARD)
                    nStates = 2;
                v_typical = length/m->numUncompressedChars + 0.0001;
                
                /* get division warp factor (prop. to prob. of change) */
                divFactor = - warpFactor * log(1.0/nStates - exp(-nStates/(nStates-1)*v_typical)/nStates);
                parLength[i+j*nRoot] += divFactor * length;
                }
            }
        }
    
    /* find the min length and the sum for the forward move */
    minLength = -1.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            if (minLength > parLength[i+j*nRoot] || minLength < 0.0)
                minLength = parLength[i+j*nRoot];
            }
    sum1 = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = sum1 + tempy;  tempc = (tempsum - sum1) - tempy;
            sum1 = tempsum;
            // sum1 += exp(minLength - parLength[i+j*nRoot]);
            // printf("%d %d %lf\n", i, j, exp(minLength - parLength[i+j*nRoot]));
            }

    /* generate a random uniform */
    ran = RandomNumber(seed) * sum1;
    
    /* select the appropriate reattachment point */
    newA = a; newC = c;
    prob = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            // prob += exp (minLength - parLength[i+j*nRoot]);
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = prob + tempy;  tempc = (tempsum - prob) - tempy;
            prob = tempsum;
            if (prob > ran) {
                /* proposed new attaching position */
                newA = pRoot[i];
                newC = pCrown[j];
                goto outLoop;
                }
            }
outLoop:;
    iA = i; jC = j;
    
    /* calculate the proposal ratio */
    (*lnProposalRatio) = parLength[i+j*nRoot] - minLength + log(sum1);
    
    /* find the min length and the sum for the backward move */
    minLength = -1.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == iA && j == jC)  // exclude new position
                continue;
            if (minLength > parLength[i+j*nRoot] || minLength < 0.0)
                minLength = parLength[i+j*nRoot];
            }
    sum2 = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == iA && j == jC)  // exclude new position
                continue;
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = sum2 + tempy;  tempc = (tempsum - sum2) - tempy;
            sum2 = tempsum;
            // sum2 += exp (minLength - parLength[i+j*nRoot]);
            // printf("%d %d %lf\n", i, j, exp(minLength - parLength[i+j*nRoot]));
            }
    
    /* calculate the proposal ratio */
    (*lnProposalRatio) += minLength - parLength[0] - log(sum2);
    
    if (moveInRoot == YES)  /* root part has changed */
        {
        /* reattach the root part */
        newB = newA->anc;
        newA->anc = u;
        if (u->left == v)
            u->right = newA;
        else
            u->left = newA;
        u->anc = newB;
        if (newB->left == newA)
            newB->left = u;
        else
            newB->right = u;

        /* if u is locked, then we have moved upwards and need to leave the u lock behind */
        if (u->isLocked == YES)
            {
            u->isLocked = NO;
            a->isLocked = YES;
            a->lockID = u->lockID;
            u->lockID = -1;
            }
        
        p = newA;
        while (p->anc != NULL)
            {
            if (p == a) break;
            p = p->anc;
            }
        if (p == a)
            {
            /* newA is descendant to a so move a->length not u->length */
            x = u->length;
            u->length = a->length;
            a->length = x;
            }

        p = b;
        while (p->anc != NULL)
            {
            if (p == newA) break;
            p = p->anc;
            }
        if (p == newA)
            {
            /* newA is ancestor to a so insert above instead of below */
            x = newA->length;
            newA->length = u->length;
            u->length = x;
            /* newA is on root path and locked, we need to transfer lock to u */
            if (newA->isLocked == YES) {
                u->isLocked = YES;
                u->lockID = newA->lockID;
                newA->isLocked = NO;
                newA->lockID = -1;
                }
            }
        
        /* hit a length with multiplier */
        x = a->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / a->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (a->length - x);
        a->length = x;

        /* hit u length with multiplier */
        x = u->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / u->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (u->length - x);
        u->length = x;

        /* hit newA length with multiplier */
        x = newA->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / newA->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (newA->length - x);
        newA->length = x;
         
        /* set tiprobs update flags */
        newA->upDateTi = YES;
        a->upDateTi = YES;
        u->upDateTi = YES;
            
        /* set flags for update of cond likes */
        p = u;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        p = b;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        }
    
    if (moveInRoot == NO)  /* crown part has changed */
        {
        r = newC;
        q = newB = newC->anc;
        /* rotate nodes from newC to c or d (whichever is closest) */
        tempc = r->length;
        while (r != c && r != d)
            {
            p = q->anc;
            /* rotate pointers of q */
            if (q->left == r)
                q->left = p;
            else
                q->right = p;
            q->anc = r;
            /* swap q and old */
            tempy = q->length;
            q->length = tempc;
            q->upDateTi = YES;
            tempc = tempy;
            /* make sure we get q and r initialized for next round */
            r = q;
            q = p;
            }
        newB->length = tempc;
            
        /* hit q length with multiplier while we are at it */
        x = q->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }        
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / q->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (q->length - x);
        q->length = x;
        q->upDateTi = YES;

        /* hit newB length with multiplier */
        x = newB->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / newB->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (newB->length - x);
        newB->length = x;
        newB->upDateTi = YES;

        /* hit newC length with multiplier */
        x = newC->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / newC->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (newC->length - x);
        newC->length = x;
        newC->upDateTi = YES;
 
        /* reattach the crown part */
        v->left = newC;
        v->right = newB;
        newC->anc = newB->anc = v;
        
        /* set flags for update of cond likes */
        p = newC;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        p = r;
        while (p->anc != NULL)
            {
            p->upDateCl = YES;
            p = p->anc;
            }
        }
    
    topologyHasChanged = YES;

    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }
    
    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);
    
    /* free up local memory */
    free (parLength); free (pRoot); free (pCrown); free (nSitesOfPat);
    
    return (NO_ERROR);
    
errorExit:
    MrBayesPrint ("%s   Problem allocating memory in Move_ParsSPR\n", spacer);
    free (parLength); free (pRoot); free (pCrown); free (nSitesOfPat);
    
    return (ERROR);
}


int Move_ParsSPRClock (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change branch lengths and topology (potentially) using SPR-type move, parsimony-biased */

    /* This move picks a branch and then chooses a reattachment point based on
       the parsimony score. On the ending branch, the attachment point is reinserted
       randomly along the branch (below the minimum age of the node). Since 2010-11-02
       the move is Metropolized to improve efficiency. */
    
    /* Since 2015-11-06, this move uses the s/n brlen approximation */
    
    int         i, j, n, division, n1=0, n2=0, n3=0, n4=0, n5=0, *nEvents, numMovableNodesOld, numMovableNodesNew;
    BitsLong    *pA, *pV, *pP, y[2];
    MrBFlt      x, oldBrlen=0.0, newBrlen=0.0, v1=0.0, v2=0.0, v3=0.0, v4=0.0, v5=0.0,
                v3new=0.0, lambda, **position=NULL, **rateMultiplier=NULL, *brlens,
                igrvar, *igrRate=NULL, nu, *tk02Rate=NULL, minLength=0.0, length=0.0,
                cumulativeProb, warpFactor, sum1, sum2, ran, divFactor, nStates, rateMult, v_approx, minV;
    CLFlt       *nSitesOfPat, *nSites, *globalNSitesOfPat;
    TreeNode    *p, *a, *b, *u, *v, *c=NULL, *d;
    Tree        *t;
    ModelInfo   *m=NULL;
    Param       *subParm;

    warpFactor = mvp[0];                  /* tuning parameter determining how heavily to weight according to parsimony scores */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    /* get model params and model info */
    m = &modelSettings[param->relParts[0]];
    
    /* get min and max brlen in relative time and subst units */
    minV = BRLENS_MIN;

#   if defined (DEBUG_ParsSPRClock)
    printf ("Before:\n");
    ShowNodes (t->root, 2, YES);
    getchar();
#   endif
    
    numMovableNodesOld=0;
    for (i=0; i<t->nNodes-2; ++i)
        {
        p = t->allDownPass[i];
        a = p->anc->left;
        b = p->anc->right;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL
            || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN))
            ++numMovableNodesOld;
        }
    
    /* pick a branch */
    do  {
        p = t->allDownPass[(int)(RandomNumber(seed) * (t->nNodes - 2))];
        a = p->anc->left;
        b = p->anc->right;
        }
    while (p->anc->isLocked == YES || p->anc->anc->anc == NULL
           || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN));
    /* skip constraints, siblings of root (and root); and consider ancestral fossils in fbd tree */

    /* set up pointers for nodes around the picked branch */
    v = p;
    u = p->anc;
    if (u->left == v)
        a = u->right;
    else
        a = u->left;
    b = u->anc;

    /* record branch length for insertion in back move */
    if (v->length > 0.0)  /* side branch, not anc fossil */
        {
        if (v->nodeDepth > a->nodeDepth)
            oldBrlen = b->nodeDepth - v->nodeDepth - 2.0*minV;
        else
            oldBrlen = b->nodeDepth - a->nodeDepth - 2.0*minV;
        }
    v1 = a->length;
    v2 = u->length;
    v3 = v->length;

    /* reassign events for CPP and adjust prior and proposal ratios for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
        {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            position = subParm->position[2*chain+state[chain]];
            rateMultiplier = subParm->rateMult[2*chain+state[chain]];
            n1 = nEvents[a->index];
            n2 = nEvents[u->index];
            n3 = nEvents[v->index];
            if (n2 > 0)
                {
                position[a->index] = (MrBFlt *) SafeRealloc ((void *) position[a->index], (n1+n2) * sizeof (MrBFlt));
                rateMultiplier[a->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[a->index], (n1+n2) * sizeof (MrBFlt));
                }
            for (j=0; j<n1; j++)
                position[a->index][j] *= v1 / (v1+v2);
            for (j=n1; j<n1+n2; j++)
                {
                position[a->index][j] = (position[u->index][j-n1] * v2 + v1) / (v1+v2);
                rateMultiplier[a->index][j] = rateMultiplier[u->index][j-n1];
                }
            nEvents[a->index] = n1+n2;
            nEvents[u->index] = 0;
            if (n2 > 0)
                {
                free (position[u->index]);
                free (rateMultiplier[u->index]);
                position[u->index] = rateMultiplier[u->index] = NULL;
                }
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] += brlens[u->index];   /* only change in effective branch lengths so far */
            }   /* end CPP events parm */
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            /* adjust prior ratio */
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            if (v->length > 0.0)
                (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[v->anc->index], nu*v->length, tk02Rate[v->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[a->anc->index], nu*a->length, tk02Rate[a->index]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*u->length, tk02Rate[u->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*(a->length+u->length), tk02Rate[a->index]);
            
            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = (tk02Rate[a->index] + tk02Rate[b->index]) / 2.0 * (a->length + u->length);
            }   /* end tk02 branch rate parameter */
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);

            /* adjust prior ratio for old branches */
            if (v->length > 0.0)
                (*lnPriorRatio) -= LnProbGamma(v->length/igrvar, v->length/igrvar, igrRate[v->index]);
            (*lnPriorRatio) -= LnProbGamma(a->length/igrvar, a->length/igrvar, igrRate[a->index]);
            (*lnPriorRatio) -= LnProbGamma(u->length/igrvar, u->length/igrvar, igrRate[u->index]);
            (*lnPriorRatio) += LnProbGamma((a->length+u->length)/igrvar, (a->length+u->length)/igrvar, igrRate[a->index]);

            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[a->index] = igrRate[a->index] * (a->length + u->length);
            } /* end igr branch rate parameter */
        }   /* next subparameter */

    /* cut tree */
    a->anc = b;
    if (b->left == u)
        b->left = a;
    else
        b->right = a;
    a->length += u->length;
    a->upDateTi = YES;

    /* get final parsimony states for the root part */
    GetParsDP (t, t->root->left, chain);
    GetParsFP (t, t->root->left->left, chain);
    GetParsFP (t, t->root->left->right, chain);

    /* get downpass parsimony states for the crown part */
    GetParsDP (t, v, chain);

    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->marked = NO;
        p->d = 0.0;
        }

    /* mark nodes in the root part of the tree, first mark a */
    a->marked = YES;
    /* then move down towards root taking constraints into account */
    p = b;
    while (p->isLocked == NO && p->anc->anc != NULL)
        {
        p->marked = YES;
        p = p->anc;
        }
    /* make sure sisters of last node are marked otherwise it will not be marked in the uppass */
    p->left->marked = YES;
    p->right->marked = YES;
    /* finally move up, skip constraints and ancestral fossil */
    for (i=t->nNodes-2; i>=0; i--)
        {
        p = t->allDownPass[i];
        if (p != u && p->marked == NO && p->anc->marked == YES && p->anc->isLocked == NO
            && p->anc->nodeDepth > v->nodeDepth + minV && p->length > 0.0)
            p->marked = YES;
        }

    /* unmark nodes if the picked branch is 0 (ancestral fossil) */
    if (v->length < TIME_MIN)
        {
        n = 0;
        for (i=0; i<t->nNodes-1; i++)
            {
            p = t->allDownPass[i];
            if (p->nodeDepth > v->nodeDepth - minV || p->anc->nodeDepth < v->nodeDepth + minV)
                p->marked = NO;
            if (p->marked == YES)
                n++;
            }
        if (n < 2)  /* no new position to move */
            {
            abortMove = YES;
            return (NO_ERROR);
            }
        }
    
    /* find number of site patterns and modify randomly */
    globalNSitesOfPat = numSitesOfPat + (chainId[chain] % chainParams.numChains) * numCompressedChars + m->compCharStart;
    nSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!nSitesOfPat)
        {
        MrBayesPrint ("%s   Problem allocating nSitesOfPat in Move_ParsSPRClock\n", spacer);
        return (ERROR);
        }
    for (i=0; i<numCompressedChars; i++)
        {
        nSitesOfPat[i] = globalNSitesOfPat[i];
        // for (j=0; j<globalNSitesOfPat[i]; j++)
        //    {
        //    ran = RandomNumber(seed);
        //    if (ran < decreaseProb)
        //        nSitesOfPat[i]--;
        //    else if (ran > 1.0 - increaseProb)
        //        nSitesOfPat[i]++;
        //    }
        }

    /* cycle through the possibilities and record the parsimony length */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO)
            continue;
        /* find the parsimony length */
        p->d = 0.0;
        for (n=0; n<t->nRelParts; n++)
            {
            division = t->relParts[n];
            
            /* Find model settings */
            m = &modelSettings[division];

            /* find nStates and ratemult */
            nStates = m->numModelStates;
            if (m->dataType == STANDARD)
                nStates = 2;
            rateMult = GetRate(division, chain);

            /* find nSitesOfPat */
            nSites = nSitesOfPat + m->compCharStart;

            /* find downpass parsimony sets for the node and its environment */
            pP   = m->parsSets[p->index     ];
            pA   = m->parsSets[p->anc->index];
            pV   = m->parsSets[v->index     ];
            
            length = 0.0;
            if (m->nParsIntsPerSite == 1)
                {
                for (j=0; j<m->numChars; j++)
                    {
                    y[0] = (pP[j] | pA[j]) & pV[j];
                    if (y[0] == 0)
                        length += nSites[j];
                    }
                }
            else /* if (m->nParsIntsPerSite == 2) */
                {
                for (j=0; j<2*m->numChars; j+=2)
                    {
                    y[0] = (pP[j] | pA[j]) & pV[j];
                    y[1] = (pP[j+1] | pA[j+1]) & pV[j+1];
                    if ((y[0] | y[1]) == 0)
                        length += nSites[j/2];
                    }
                }

            /* find nStates and v approximation using parsimony-based s/n approximation */
            nStates = m->numModelStates;
            if (m->dataType == STANDARD)
                nStates = 2;
            v_approx = length/m->numUncompressedChars + 0.0001;
            
            /* get division warp factor (prop. to prob. of change) */
            divFactor = - warpFactor * log(1.0/nStates - exp(-nStates/(nStates-1)*v_approx)/nStates);
            
            p->d += divFactor * length;
            }
        }

    /* find the min length and the sum for the forward move */
    minLength = -1.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO || p == a)
            continue;
        if (minLength < 0.0 || p->d < minLength)
            minLength = p->d;
        }
    sum1 = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES && p != a)
            sum1 += exp (minLength - p->d);
        }

    /* generate a random uniform */
    ran = RandomNumber(seed) * sum1;

    /* select the appropriate reattachment point (not a!) */
    cumulativeProb = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES && p != a)
            {
            c = p;
            cumulativeProb += exp (minLength - p->d);
            if (cumulativeProb > ran)
                break;
            }
        }

    /* calculate the proposal ratio */
    (*lnProposalRatio) = c->d - minLength + log(sum1);

    /* find the min length and the sum for the backward move */
    minLength = -1.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == NO || p == c)
            continue;
        if (minLength < 0.0 || p->d < minLength)
            minLength = p->d;
        }
    sum2 = 0.0;
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->marked == YES && p != c)
            sum2 += exp (minLength - p->d);
        }
    
    /* calculate the proposal ratio */
    (*lnProposalRatio) += minLength - a->d - log(sum2);

    /* reattach u */
    d = c->anc;
    c->anc = u;
    if (u->left == v)
        u->right = c;
    else
        u->left = c;
    u->anc = d;
    if (d->left == c)
        d->left = u;
    else
        d->right = u;

    if (v->length > 0.0)  /* side branch, not anc fossil */
        {
        if (c->nodeDepth > v->nodeDepth)
            newBrlen = d->nodeDepth - c->nodeDepth - 2.0*minV;
        else
            newBrlen = d->nodeDepth - v->nodeDepth - 2.0*minV;
        if (newBrlen <= 0.0)
            {
            abortMove = YES;
            free (nSitesOfPat);
            return (NO_ERROR);
            }

        /* adjust lengths */
        u->nodeDepth = d->nodeDepth - minV - RandomNumber(seed) * newBrlen;
        v->length = u->nodeDepth - v->nodeDepth;
            
        /* calculate proposal ratio for tree change */
        (*lnProposalRatio) += log (newBrlen / oldBrlen);
        }
    u->length = d->nodeDepth - u->nodeDepth;
    c->length = u->nodeDepth - c->nodeDepth;
    
    v3new = v->length;
    v4 = c->length;
    v5 = u->length;

    /* reassign events for CPP and adjust prior and proposal ratios for relaxed clock models */
    for (i=0; i<param->subParams[0]->nSubParams; i++)
        {
        subParm = param->subParams[0]->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            position = subParm->position[2*chain+state[chain]];
            rateMultiplier = subParm->rateMult[2*chain+state[chain]];
            for (j=0; j<nEvents[c->index]; j++)
                {
                if (position[c->index][j] > v4 / (v4+v5))
                    break;
                }
            n4 = j;
            n5 = nEvents[c->index] - j;
            nEvents[u->index] = n5;
            if (n5 > 0)
                {
                position[u->index] = (MrBFlt *) SafeRealloc ((void *) position[u->index], n5 * sizeof (MrBFlt));
                rateMultiplier[u->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[u->index], n5 * sizeof (MrBFlt));            
                for (j=n4; j<nEvents[c->index]; j++)
                    {
                    position[u->index][j-n4] = (position[c->index][j] * (v4+v5) - v4) / v5;
                    rateMultiplier[u->index][j-n4] = rateMultiplier[c->index][j];
                    }
                if (n4 > 0)
                    {
                    position[c->index] = (MrBFlt *) SafeRealloc ((void *) position[c->index], n4 * sizeof (MrBFlt));
                    rateMultiplier[c->index] = (MrBFlt *) SafeRealloc ((void *) rateMultiplier[c->index], n4 * sizeof (MrBFlt));
                    for (j=0; j<n4; j++)
                        position[c->index][j] *= ((v4+v5) / v4);
                    }
                else
                    {
                    free (position[c->index]);
                    free (rateMultiplier[c->index]);
                    position[c->index] = rateMultiplier[c->index] = NULL;
                    }
                nEvents[c->index] = n4;
                }
            else
                {
                for (j=0; j<nEvents[c->index]; j++)
                    position[c->index][j] *= ((v4+v5) / v4);
                }

            /* adjust proposal ratio */
            (*lnProposalRatio) += n3 * log (v3new / v3);

            /* adjust prior ratio */
            lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);
            (*lnPriorRatio) += lambda * (v3 - v3new);

            /* update effective branch lengths */
            if (UpdateCppEvolLengths (subParm, a, chain) == ERROR)
                {
                abortMove = YES;
                free (nSitesOfPat);
                return (NO_ERROR);
                }

            if (UpdateCppEvolLengths (subParm, u, chain) == ERROR)
                {
                abortMove = YES;
                free (nSitesOfPat);
                return (NO_ERROR);
                }
            }   /* end cpp events parameter */
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            /* adjust prior ratio */
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            (*lnPriorRatio) -= LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*(c->length+u->length), tk02Rate[c->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[c->anc->index], nu*c->length, tk02Rate[c->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[u->anc->index], nu*u->length, tk02Rate[u->index]);
            if (v->length > 0.0)
                (*lnPriorRatio) += LnProbTK02LogNormal(tk02Rate[v->anc->index], nu*v->length, tk02Rate[v->index]);

            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[c->index] = c->length * (tk02Rate[c->index] + tk02Rate[c->anc->index]) / 2.0;
            brlens[v->index] = v->length * (tk02Rate[v->index] + tk02Rate[v->anc->index]) / 2.0;
            brlens[u->index] = u->length * (tk02Rate[u->index] + tk02Rate[u->anc->index]) / 2.0;
            }   /* end tk02 branch rate parameter */
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            /* adjust prior ratio */
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            (*lnPriorRatio) -= LnProbGamma ((c->length+u->length)/igrvar, (c->length+u->length)/igrvar, igrRate[c->index]);
            (*lnPriorRatio) += LnProbGamma (c->length/igrvar, c->length/igrvar, igrRate[c->index]);
            (*lnPriorRatio) += LnProbGamma (u->length/igrvar, u->length/igrvar, igrRate[u->index]);
            if (v->length > 0.0)
                (*lnPriorRatio) += LnProbGamma (v->length/igrvar, v->length/igrvar, igrRate[v->index]);

            /* adjust effective branch lengths */
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            brlens[v->index] = igrRate[v->index] * v->length;
            brlens[u->index] = igrRate[u->index] * u->length;
            brlens[c->index] = igrRate[c->index] * c->length;
            }   /* end igr branch rate parameter */
        }   /* next subparameter */
    
    /* set tiprobs update flags */
    c->upDateTi = YES;
    u->upDateTi = YES;
    v->upDateTi = YES;

    /* set flags for update of cond likes down to root */
    p = u;
    while (p->anc != NULL)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }
    p = b;
    while (p->anc != NULL)
        {
        p->upDateCl = YES; 
        p = p->anc;
        }

    /* get down pass sequence */
    GetDownPass (t);

    /* adjust proposal ratio for number of movable nodes */
    numMovableNodesNew=0;
    for (i=0; i<t->nNodes-2; ++i)
        {
        p = t->allDownPass[i];
        a = p->anc->left;
        b = p->anc->right;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL
            || (p == b && a->length < TIME_MIN) || (p == a && b->length < TIME_MIN))
            ++numMovableNodesNew;
        }
    if (numMovableNodesOld != numMovableNodesNew)
        (*lnProposalRatio) += log(numMovableNodesOld/numMovableNodesNew);

    /* adjust prior ratio for clock tree */
    if (LogClockTreePriorRatio (param, chain, &x) == ERROR)
        {
        free (nSitesOfPat);
        return (ERROR);
        }
    (*lnPriorRatio) += x;

#   if defined (DEBUG_ParsSPRClock)
    ShowNodes (t->root, 2, YES);
    printf ("After\nProposal ratio: %f\n",(*lnProposalRatio));
    printf ("v: %d  u: %d  a: %d  b: %d c: %d\n",v->index, u->index, a->index, b->index, c->index);
    getchar();
#   endif

    free (nSitesOfPat);
    return (NO_ERROR);
}


int Move_ParsTBR1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change topology and map branch lengths using TBR-type move biased according to parsimony scores,
       controlled by a window defined by a certain node distance radius. */
    
    int         i, j, k, n, division, topologyHasChanged, nNeighbor, nRoot, nCrown, iA, jC, isVPriorExp;
    BitsLong    *pA, *pB, *pC, *pD, y[2];
    MrBFlt      x, minV, maxV, brlensExp=0.0, minLength=0.0, length=0.0, *parLength=NULL, prob, ran, tuning, warpFactor,
                sum1, sum2, tempsum, tempc, tempy;
    CLFlt       *nSites, *nSitesOfPat=NULL, *globalNSitesOfPat;
    TreeNode    *p, *q, *r, *a, *b, *u, *v, *c, *d, *newB, *newA, *newC, **pRoot=NULL, **pCrown=NULL;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m=NULL;
    
    warpFactor = mvp[0];                  /* tuning parameter determining how heavily to weight according to parsimony scores */
//  increaseProb = decreaseProb = mvp[1]; /* reweighting probabilities */
//  v_typical = mvp[2];                   /* typical branch length for conversion of parsimony score to log prob ratio */
    tuning = mvp[3];                      /* multiplier tuning parameter */
    nNeighbor = (int)mvp[4];              /* distance to move picked branch in root and crown part */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;
    
    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1] < BRLENS_MAX ? mp->brlensUni[1] : BRLENS_MAX;
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);
    
    /* set topologyHasChanged to NO */
    topologyHasChanged = NO;
    
    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->x = 0;
        p->marked = NO;
        }

    /* pick an internal branch */
    do  {
        p = t->intDownPass[(int)(RandomNumber(seed)*(t->nIntNodes-1))];
        q = p->anc->left;
        if (q == p)  q = p->anc->right;
        i = j = 0;
        if (q->isLocked == YES || q->left == NULL)
            i++;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL)
            i++;
        if (p->left->isLocked == YES || p->left->left == NULL)
            j++;
        if (p->right->isLocked == YES || p->right->left == NULL)
            j++;
        }
    while (i == 2 && j == 2);
    
    /* set up pointers for nodes around the picked branch */
    v = p;            u = p->anc;
    c = p->left;      d = p->right;
    if (u->left == v) a = u->right;
    else              a = u->left;
    b = u->anc;
    /* clip root part of tree */
    a->anc = b;
    if (b->left == u) b->left = a;
    else              b->right = a;
    /* clip crown part of tree */
    c->anc = d;
    d->anc = c;
    /* should never change u, v, a, b, c, d pointers */

    /* mark nodes nNeighbor away in root (negative) and crown (positive) respecting constraints */
    /* first move down towards root */
    nRoot = nCrown = 0;
    if (u->isLocked == NO)
        {
        p = a; q = b; n = 0;
        while (q->anc != NULL)
            {
            q->marked = YES;
            q->x = n;    // temporary, for MarkDistance below
            if (q->left == p)
                MarkDistance(q->right, YES, nNeighbor, &nRoot);
            else
                MarkDistance(q->left,  YES, nNeighbor, &nRoot);
            q->x = --n;  // final
            nRoot++;
            if (q->isLocked == YES || abs(q->x) >= nNeighbor)
                break;
            p = q; q = q->anc;
            }
        }
    /* then move up in root part */
    a->marked = YES; nRoot++;
    if (a->isLocked == NO)
        {
        MarkDistance(a->left,  YES, nNeighbor, &nRoot);
        MarkDistance(a->right, YES, nNeighbor, &nRoot);
        }
    /* finally in crown part */
    c->marked = YES; nCrown++;
    if (c->isLocked == NO)
        {
        MarkDistance(c->left,  NO, nNeighbor, &nCrown);
        MarkDistance(c->right, NO, nNeighbor, &nCrown);
        }
    if (d->isLocked == NO)
        {
        MarkDistance(d->left,  NO, nNeighbor, &nCrown);
        MarkDistance(d->right, NO, nNeighbor, &nCrown);
        }

    /* need to alloc a matrix for parsimony lengths, an array of pointers to crown part,
       and an array of pointers to root part. */
    parLength = (MrBFlt *) SafeCalloc ((size_t)nRoot * (size_t)nCrown, sizeof(MrBFlt));
    pRoot  = (TreeNode **) SafeCalloc(nRoot,  sizeof(TreeNode *));
    pCrown = (TreeNode **) SafeCalloc(nCrown, sizeof(TreeNode *));
    if (!parLength || !pRoot || !pCrown)  goto errorExit;
    /* starting position */
    pRoot[0] = a; pCrown[0] = c;
    for (i=j=1, n=t->nNodes-2; n>=0; n--)
        {  /* and the rest */
        p = t->allDownPass[n];
        if (p->marked == YES && p->x < 0)
            pRoot[i++] = p;
        if (p->marked == YES && p->x > 0)
            pCrown[j++] = p;
        }
    assert (i==nRoot && j==nCrown);

    /* get final parsimony state sets for the root part */
    GetParsDP (t, t->root->left, chain);
    GetParsFP (t, t->root->left, chain);
    /* get final parsimony state sets for the crown part */
    GetParsDP (t, c, chain);
    GetParsDP (t, d, chain);
    GetParsFP (t, c, chain);
    GetParsFP (t, d, chain);

    /* find number of site patterns and modify randomly */
    globalNSitesOfPat = numSitesOfPat + (chainId[chain] % chainParams.numChains) * numCompressedChars + m->compCharStart;
    nSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!nSitesOfPat)  goto errorExit;
    for (i=0; i<numCompressedChars; i++)
        {
        nSitesOfPat[i] = globalNSitesOfPat[i];
    /*  for (j=0; j<globalNSitesOfPat[i]; j++)
            {
            ran = RandomNumber(seed);
            if (ran < decreaseProb)
                nSitesOfPat[i]--;
            else if (ran > 1.0 - increaseProb)
                nSitesOfPat[i]++;
            }  // this is not used at the moment */
        }

    /* cycle through the possibilities and record the parsimony length */
    for (j=0; j<nCrown; j++)
        {
        for (i=0; i<nRoot; i++)
            {
            parLength[i+j*nRoot] = 0.0;
            for (n=0; n<t->nRelParts; n++)
                {
                division = t->relParts[n];
                
                /* Find model settings */
                m = &modelSettings[division];
                
                /* find nSitesOfPat */
                nSites = nSitesOfPat + m->compCharStart;
                    
                /* find downpass parsimony sets for the potential new connection nodes and their environment */
                pA = m->parsSets[pRoot[i]->index];
                pB = m->parsSets[pRoot[i]->anc->index];
                pC = m->parsSets[pCrown[j]->index];
                pD = m->parsSets[pCrown[j]->anc->index];
                    
                length = 0.0;
                if (m->nParsIntsPerSite == 1)
                    {
                    for (k=0; k<m->numChars; k++)
                        {
                        y[0] = (pC[k] | pD[k]) & (pA[k] | pB[k]);
                        if (y[0] == 0)
                            length += nSites[k];
                        }
                    }
                else /* if (m->nParsIntsPerSite == 2) */
                    {
                    for (k=0; k<2*m->numChars; k+=2)
                        {
                        y[0] = (pC[k] | pD[k]) & (pA[k] | pB[k]);
                        y[1] = (pC[k+1] | pD[k+1]) & (pA[k+1] | pB[k+1]);;
                        if ((y[0] | y[1]) == 0)
                            length += nSites[k/2];
                        }
                    }
                parLength[i+j*nRoot] += warpFactor * length;
                }
            }
        }

    /* find the min length and the sum for the forward move */
    minLength = -1.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            if (minLength > parLength[i+j*nRoot] || minLength < 0.0)
                minLength = parLength[i+j*nRoot];
            }
    sum1 = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = sum1 + tempy;  tempc = (tempsum - sum1) - tempy;
            sum1 = tempsum;
            // sum1 += exp(minLength - parLength[i+j*nRoot]);
            // printf("%d %d %lf\n", i, j, exp(minLength - parLength[i+j*nRoot]));
            }

    /* generate a random uniform */
    ran = RandomNumber(seed) * sum1;

    /* select the appropriate reattachment point */
    newA = a; newC = c;
    prob = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            // prob += exp (minLength - parLength[i+j*nRoot]);
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = prob + tempy;  tempc = (tempsum - prob) - tempy;
            prob = tempsum;
            if (prob > ran) {
                /* proposed new attaching position */
                newA = pRoot[i];
                newC = pCrown[j];
                goto outLoop;
                }
            }
outLoop:;
    iA = i; jC = j;

    /* calculate the proposal ratio */
    (*lnProposalRatio) = parLength[i+j*nRoot] - minLength + log(sum1);

    /* find the min length and the sum for the backward move */
    minLength = -1.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == iA && j == jC)  // exclude new position
                continue;
            if (minLength > parLength[i+j*nRoot] || minLength < 0.0)
                minLength = parLength[i+j*nRoot];
            }
    sum2 = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == iA && j == jC)  // exclude new position
                continue;
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = sum2 + tempy;  tempc = (tempsum - sum2) - tempy;
            sum2 = tempsum;
            // sum2 += exp (minLength - parLength[i+j*nRoot]);
            // printf("%d %d %lf\n", i, j, exp(minLength - parLength[i+j*nRoot]));
            }

    /* calculate the proposal ratio */
    (*lnProposalRatio) += minLength - parLength[0] - log(sum2);

    /* reattach the root part */
    newB = newA->anc;
    newA->anc = u;
    if (u->left == v)
        u->right = newA;
    else
        u->left = newA;
    u->anc = newB;
    if (newB->left == newA)
        newB->left = u;
    else
        newB->right = u;

    /* transfer lock and reassign branch lengths, if necessary */
    if (newA != a)
        {
        /* if u is locked, then we have moved upwards and need to leave the u lock behind */
        if (u->isLocked == YES)
            {
            u->isLocked = NO;
            a->isLocked = YES;
            a->lockID = u->lockID;
            u->lockID = -1;
            }

        p = newA;
        while (p->anc != NULL)
            {
            if (p == a) break;
            p = p->anc;
            }
        if (p == a)
            {
            /* newA is descendant to a so move a->length not u->length */
            x = u->length;
            u->length = a->length;
            a->length = x;
            }

        p = b;
        while (p->anc != NULL)
            {
            if (p == newA) break;
            p = p->anc;
            }
        if (p == newA)
            {
            /* newA is ancestor to a so insert above instead of below */
            x = newA->length;
            newA->length = u->length;
            u->length = x;
            /* newA is on root path and locked, we need to transfer lock to u */
            if (newA->isLocked == YES) {
                u->isLocked = YES;
                u->lockID = newA->lockID;
                newA->isLocked = NO;
                newA->lockID = -1;
                }
            }
            
        /* hit u length with multiplier */
        x = u->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / u->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (u->length - x);
        u->length = x;

        /* set tiprobs update flags */
        u->upDateTi = YES;
        newA->upDateTi = YES;
        a->upDateTi = YES;
        }
    
    r = newC;
    q = newB = newC->anc;
    if (newC != c)  // crown part has changed
        {
        /* rotate nodes from newC to c or d (whichever is closest) */
        tempc = r->length;
        while (r != c && r != d)
            {
            p = q->anc;
            /* rotate pointers of q */
            if (q->left == r)
                q->left = p;
            else
                q->right = p;
            q->anc = r;
            /* swap q and old */
            tempy = q->length;
            q->length = tempc;
            q->upDateTi = YES;
            tempc = tempy;
            /* make sure we get q and r initialized for next round */
            r = q;
            q = p;
            }
        newB->length = tempc;
        
        /* hit newB length with multiplier */
        x = newB->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / newB->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (newB->length - x);
        newB->length = x;
        newB->upDateTi = YES;
        }
    
    /* reattach the crown part */
    v->left = newC;
    v->right = newB;
    newC->anc = newB->anc = v;
    
    topologyHasChanged = YES;
    
    /* hit v length with multiplier */
    x = v->length * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV) x = minV * minV / x;
        if (x > maxV) x = maxV * maxV / x;
        }
    /* calculate proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / v->length);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (v->length - x);
    v->length = x;
    v->upDateTi = YES;
 
    /* set flags for update of cond likes */
    p = u;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
    p = b;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
    p = newC;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
    p = r;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
    
    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);
    
    /* free up local memory */
    free (parLength); free (pRoot); free (pCrown); free (nSitesOfPat);

    return (NO_ERROR);

errorExit:
    MrBayesPrint ("%s   Problem allocating memory in Move_ParsTBR\n", spacer);
    free (parLength); free (pRoot); free (pCrown); free (nSitesOfPat);

    return (ERROR);
}


int Move_ParsTBR2 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* Change topology and map branch lengths using TBR-type move biased according to parsimony scores,
       controlled by a window defined by a certain node distance radius. */
    
    int         i, j, k, n, division, topologyHasChanged, nNeighbor, nRoot, nCrown, iA, jC, isVPriorExp;
    BitsLong    *pA, *pB, *pC, *pD, y[2];
    MrBFlt      x, minV, maxV, brlensExp=0.0, minLength=0.0, length=0.0, *parLength=NULL, prob, ran, tuning, warpFactor,
                v_typical, divFactor, nStates, sum1, sum2, tempsum, tempc, tempy;
    CLFlt       *nSites, *nSitesOfPat=NULL, *globalNSitesOfPat;
    TreeNode    *p, *q, *r, *a, *b, *u, *v, *c, *d, *newB, *newA, *newC, **pRoot=NULL, **pCrown=NULL;
    Tree        *t;
    ModelParams *mp;
    ModelInfo   *m=NULL;
    
    warpFactor = mvp[0];                  /* tuning parameter determining how heavily to weight according to parsimony scores */
//  increaseProb = decreaseProb = mvp[1]; /* reweighting probabilities */
    v_typical = mvp[2];                   /* typical branch length for conversion of parsimony score to log prob ratio */
    tuning = mvp[3];                      /* multiplier tuning parameter */
    nNeighbor = (int)mvp[4];              /* distance to move picked branch in root and crown part */

    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;
    
    /* get model params and model info */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* max and min brlen */
    if (param->subParams[0]->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0] > BRLENS_MIN ? mp->brlensUni[0] : BRLENS_MIN;
        maxV = mp->brlensUni[1] < BRLENS_MAX ? mp->brlensUni[1] : BRLENS_MAX;
        isVPriorExp = NO;
        }
    else if (param->subParams[0]->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->subParams[0]->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->subParams[0]->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);
    
    /* set topologyHasChanged to NO */
    topologyHasChanged = NO;
    
    /* reset node variables that will be used */
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->x = 0;
        p->marked = NO;
        }

    /* pick an internal branch */
    do  {
        p = t->intDownPass[(int)(RandomNumber(seed)*(t->nIntNodes-1))];
        q = p->anc->left;
        if (q == p)  q = p->anc->right;
        i = j = 0;
        if (q->isLocked == YES || q->left == NULL)
            i++;
        if (p->anc->isLocked == YES || p->anc->anc->anc == NULL)
            i++;
        if (p->left->isLocked == YES || p->left->left == NULL)
            j++;
        if (p->right->isLocked == YES || p->right->left == NULL)
            j++;
        }
    while (i == 2 && j == 2);
    
    /* set up pointers for nodes around the picked branch */
    v = p;            u = p->anc;
    c = p->left;      d = p->right;
    if (u->left == v) a = u->right;
    else              a = u->left;
    b = u->anc;
    /* clip root part of tree */
    a->anc = b;
    if (b->left == u) b->left = a;
    else              b->right = a;
    /* clip crown part of tree */
    c->anc = d;
    d->anc = c;
    /* should never change u, v, a, b, c, d pointers */

    /* mark nodes nNeighbor away in root (negative) and crown (positive) respecting constraints */
    /* first move down towards root */
    nRoot = nCrown = 0;
    if (u->isLocked == NO)
        {
        p = a; q = b; n = 0;
        while (q->anc != NULL)
            {
            q->marked = YES;
            q->x = n;    // temporary, for MarkDistance below
            if (q->left == p)
                MarkDistance(q->right, YES, nNeighbor, &nRoot);
            else
                MarkDistance(q->left,  YES, nNeighbor, &nRoot);
            q->x = --n;  // final
            nRoot++;
            if (q->isLocked == YES || abs(q->x) >= nNeighbor)
                break;
            p = q; q = q->anc;
            }
        }
    /* then move up in root part */
    a->marked = YES; nRoot++;
    if (a->isLocked == NO)
        {
        MarkDistance(a->left,  YES, nNeighbor, &nRoot);
        MarkDistance(a->right, YES, nNeighbor, &nRoot);
        }
    /* finally in crown part */
    c->marked = YES; nCrown++;
    if (c->isLocked == NO)
        {
        MarkDistance(c->left,  NO, nNeighbor, &nCrown);
        MarkDistance(c->right, NO, nNeighbor, &nCrown);
        }
    if (d->isLocked == NO)
        {
        MarkDistance(d->left,  NO, nNeighbor, &nCrown);
        MarkDistance(d->right, NO, nNeighbor, &nCrown);
        }

    /* need to alloc a matrix for parsimony lengths, an array of pointers to crown part,
       and an array of pointers to root part. */
    parLength = (MrBFlt *) SafeCalloc ((size_t)nRoot * (size_t)nCrown, sizeof(MrBFlt));
    pRoot  = (TreeNode **) SafeCalloc(nRoot,  sizeof(TreeNode *));
    pCrown = (TreeNode **) SafeCalloc(nCrown, sizeof(TreeNode *));
    if (!parLength || !pRoot || !pCrown)  goto errorExit;
    /* starting position */
    pRoot[0] = a; pCrown[0] = c;
    for (i=j=1, n=t->nNodes-2; n>=0; n--)
        {  /* and the rest */
        p = t->allDownPass[n];
        if (p->marked == YES && p->x < 0)
            pRoot[i++] = p;
        if (p->marked == YES && p->x > 0)
            pCrown[j++] = p;
        }
    assert (i==nRoot && j==nCrown);

    /* get final parsimony state sets for the root part */
    GetParsDP (t, t->root->left, chain);
    GetParsFP (t, t->root->left, chain);
    /* get final parsimony state sets for the crown part */
    GetParsDP (t, c, chain);
    GetParsDP (t, d, chain);
    GetParsFP (t, c, chain);
    GetParsFP (t, d, chain);

    /* find number of site patterns and modify randomly */
    globalNSitesOfPat = numSitesOfPat + (chainId[chain] % chainParams.numChains) * numCompressedChars + m->compCharStart;
    nSitesOfPat = (CLFlt *) SafeCalloc (numCompressedChars, sizeof(CLFlt));
    if (!nSitesOfPat)  goto errorExit;
    for (i=0; i<numCompressedChars; i++)
        {
        nSitesOfPat[i] = globalNSitesOfPat[i];
    /*  for (j=0; j<globalNSitesOfPat[i]; j++)
            {
            ran = RandomNumber(seed);
            if (ran < decreaseProb)
                nSitesOfPat[i]--;
            else if (ran > 1.0 - increaseProb)
                nSitesOfPat[i]++;
            }  // this is not used at the moment */
        }

    /* cycle through the possibilities and record the parsimony length */
    for (j=0; j<nCrown; j++)
        {
        for (i=0; i<nRoot; i++)
            {
            parLength[i+j*nRoot] = 0.0;
            for (n=0; n<t->nRelParts; n++)
                {
                division = t->relParts[n];
                
                /* Find model settings */
                m = &modelSettings[division];
                
                /* find nSitesOfPat */
                nSites = nSitesOfPat + m->compCharStart;
                
                /* find downpass parsimony sets for the potential new connection nodes and their environment */
                pA = m->parsSets[pRoot[i]->index];
                pB = m->parsSets[pRoot[i]->anc->index];
                pC = m->parsSets[pCrown[j]->index];
                pD = m->parsSets[pCrown[j]->anc->index];
                    
                length = 0.0;
                if (m->nParsIntsPerSite == 1)
                    {
                    for (k=0; k<m->numChars; k++)
                        {
                        y[0] = (pC[k] | pD[k]) & (pA[k] | pB[k]);
                        if (y[0] == 0)
                            length += nSites[k];
                        }
                    }
                else /* if (m->nParsIntsPerSite == 2) */
                    {
                    for (k=0; k<2*m->numChars; k+=2)
                        {
                        y[0] = (pC[k] | pD[k]) & (pA[k] | pB[k]);
                        y[1] = (pC[k+1] | pD[k+1]) & (pA[k+1] | pB[k+1]);;
                        if ((y[0] | y[1]) == 0)
                            length += nSites[k/2];
                        }
                    }
                    
                /* find nStates and ratemult */
                nStates = m->numModelStates;
                if (m->dataType == STANDARD)
                    nStates = 2;
                v_typical = length/m->numUncompressedChars + 0.0001;
                
                /* get division warp factor (prop. to prob. of change) */
                divFactor = - warpFactor * log(1.0/nStates - exp(-nStates/(nStates-1)*v_typical)/nStates);
                parLength[i+j*nRoot] += divFactor * length;
                }
            }
        }

    /* find the min length and the sum for the forward move */
    minLength = -1.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            if (minLength > parLength[i+j*nRoot] || minLength < 0.0)
                minLength = parLength[i+j*nRoot];
            }
    sum1 = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = sum1 + tempy;  tempc = (tempsum - sum1) - tempy;
            sum1 = tempsum;
            // sum1 += exp(minLength - parLength[i+j*nRoot]);
            // printf("%d %d %lf\n", i, j, exp(minLength - parLength[i+j*nRoot]));
            }

    /* generate a random uniform */
    ran = RandomNumber(seed) * sum1;

    /* select the appropriate reattachment point */
    newA = a; newC = c;
    prob = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == 0 && j == 0)  // exclude original position
                continue;
            // prob += exp (minLength - parLength[i+j*nRoot]);
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = prob + tempy;  tempc = (tempsum - prob) - tempy;
            prob = tempsum;
            if (prob > ran) {
                /* proposed new attaching position */
                newA = pRoot[i];
                newC = pCrown[j];
                goto outLoop;
                }
            }
outLoop:;
    iA = i; jC = j;

    /* calculate the proposal ratio */
    (*lnProposalRatio) = parLength[i+j*nRoot] - minLength + log(sum1);

    /* find the min length and the sum for the backward move */
    minLength = -1.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == iA && j == jC)  // exclude new position
                continue;
            if (minLength > parLength[i+j*nRoot] || minLength < 0.0)
                minLength = parLength[i+j*nRoot];
            }
    sum2 = 0.0; tempc = 0.0;
    for (j=0; j<nCrown; j++)
        for (i=0; i<nRoot; i++)
            {
            if (i == iA && j == jC)  // exclude new position
                continue;
            /* Kahan summation to reduce numerical error */
            tempy = exp(minLength - parLength[i+j*nRoot]) - tempc;
            tempsum = sum2 + tempy;  tempc = (tempsum - sum2) - tempy;
            sum2 = tempsum;
            // sum2 += exp (minLength - parLength[i+j*nRoot]);
            // printf("%d %d %lf\n", i, j, exp(minLength - parLength[i+j*nRoot]));
            }

    /* calculate the proposal ratio */
    (*lnProposalRatio) += minLength - parLength[0] - log(sum2);

    /* reattach the root part */
    newB = newA->anc;
    newA->anc = u;
    if (u->left == v)
        u->right = newA;
    else
        u->left = newA;
    u->anc = newB;
    if (newB->left == newA)
        newB->left = u;
    else
        newB->right = u;

    /* transfer lock and reassign branch lengths, if necessary */
    if (newA != a)
        {
        /* if u is locked, then we have moved upwards and need to leave the u lock behind */
        if (u->isLocked == YES)
            {
            u->isLocked = NO;
            a->isLocked = YES;
            a->lockID = u->lockID;
            u->lockID = -1;
            }

        p = newA;
        while (p->anc != NULL)
            {
            if (p == a) break;
            p = p->anc;
            }
        if (p == a)
            {
            /* newA is descendant to a so move a->length not u->length */
            x = u->length;
            u->length = a->length;
            a->length = x;
            }

        p = b;
        while (p->anc != NULL)
            {
            if (p == newA) break;
            p = p->anc;
            }
        if (p == newA)
            {
            /* newA is ancestor to a so insert above instead of below */
            x = newA->length;
            newA->length = u->length;
            u->length = x;
            /* newA is on root path and locked, we need to transfer lock to u */
            if (newA->isLocked == YES) {
                u->isLocked = YES;
                u->lockID = newA->lockID;
                newA->isLocked = NO;
                newA->lockID = -1;
                }
            }
            
        /* hit u length with multiplier */
        x = u->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / u->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (u->length - x);
        u->length = x;

        /* set tiprobs update flags */
        u->upDateTi = YES;
        newA->upDateTi = YES;
        a->upDateTi = YES;
        }
    
    r = newC;
    q = newB = newC->anc;
    if (newC != c)  // crown part has changed
        {
        /* rotate nodes from newC to c or d (whichever is closest) */
        tempc = r->length;
        while (r != c && r != d)
            {
            p = q->anc;
            /* rotate pointers of q */
            if (q->left == r)
                q->left = p;
            else
                q->right = p;
            q->anc = r;
            /* swap q and old */
            tempy = q->length;
            q->length = tempc;
            q->upDateTi = YES;
            tempc = tempy;
            /* make sure we get q and r initialized for next round */
            r = q;
            q = p;
            }
        newB->length = tempc;
        
        /* hit newB length with multiplier */
        x = newB->length * exp(tuning * (RandomNumber(seed) - 0.5));
        while (x < minV || x > maxV)
            {
            if (x < minV) x = minV * minV / x;
            if (x > maxV) x = maxV * maxV / x;
            }
        /* calculate proposal and prior ratio based on length modification */
        (*lnProposalRatio) += log (x / newB->length);
        if (isVPriorExp == YES)
            (*lnPriorRatio) += brlensExp * (newB->length - x);
        newB->length = x;
        newB->upDateTi = YES;
        }
    
    /* reattach the crown part */
    v->left = newC;
    v->right = newB;
    newC->anc = newB->anc = v;
    
    topologyHasChanged = YES;
    
    /* hit v length with multiplier */
    x = v->length * exp(tuning * (RandomNumber(seed) - 0.5));
    while (x < minV || x > maxV)
        {
        if (x < minV) x = minV * minV / x;
        if (x > maxV) x = maxV * maxV / x;
        }
    /* calculate proposal and prior ratio based on length modification */
    (*lnProposalRatio) += log (x / v->length);
    if (isVPriorExp == YES)
        (*lnPriorRatio) += brlensExp * (v->length - x);
    v->length = x;
    v->upDateTi = YES;
 
    /* set flags for update of cond likes */
    p = u;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
    p = b;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
    p = newC;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
    p = r;
    while (p->anc != NULL)
        {
        p->upDateCl = YES;
        p = p->anc;
        }
    
    /* get down pass sequence if tree topology has changed */
    if (topologyHasChanged == YES)
        {
        GetDownPass (t);
        }

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);
    
    /* free up local memory */
    free (parLength); free (pRoot); free (pCrown); free (nSitesOfPat);

    return (NO_ERROR);

errorExit:
    MrBayesPrint ("%s   Problem allocating memory in Move_ParsTBR\n", spacer);
    free (parLength); free (pRoot); free (pCrown); free (nSitesOfPat);

    return (ERROR);
}


int Move_Pinvar (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change proportion of invariable sites (pInvar) */

    int             i, c, isValidP, *rateCat, nGammaCats;
    MrBFlt          oldP, newP, window, minP, maxP, ran, lnInvarRatio, lnVarRatio;
    CLFlt           *nSitesOfPat;
    ModelParams     *mp;
    ModelInfo       *m;

    /* get size of window, centered on current pInvar value */
    window = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* get minimum and maximum values for pInvar */
    minP = mp->pInvarUni[0];
    maxP = mp->pInvarUni[1];

    /* get old value of pInvar */
    newP = oldP = *GetParamVals(param, chain, state[chain]);

    /* change value for pInvar */
    ran = RandomNumber(seed);
    if (maxP-minP < window)
        {
        window = maxP-minP;
        }

    newP = oldP + window * (ran - 0.5);

    /* check validity */
    isValidP = NO;
    do
        {
        if (newP < minP)
            newP = 2* minP - newP;
        else if (newP > maxP)
            newP = 2 * maxP - newP;
        else
            isValidP = YES;
        } while (isValidP == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* get prior ratio */
    *lnPriorRatio = 0.0;
    lnInvarRatio = log(newP) - log(oldP);
    lnVarRatio = log(1.0-newP) - log(1.0-oldP);
    for (i=0; i<param->nRelParts; i++)
        {
        m = &modelSettings[param->relParts[i]];
        if (m->gibbsGamma == YES)
            {
            /* find rate category index and number of gamma categories */
            rateCat = m->tiIndex + chain * m->numChars;
            nGammaCats = m->numRateCats;

            /* find nSitesOfPat */
            nSitesOfPat = numSitesOfPat + (chainId[chain] % chainParams.numChains)*numCompressedChars + m->compCharStart;
            
            /* loop over characters */
            for (c=0; c<m->numChars; c++)
                {
                if (rateCat[c] < nGammaCats)
                    *lnPriorRatio += lnVarRatio * nSitesOfPat[c];
                else
                    *lnPriorRatio += lnInvarRatio * nSitesOfPat[c];
                }
            }
        }
    
    /* copy new pInvar value back */
    *GetParamVals(param, chain, state[chain]) = newP;

    /* Set update flags for all partitions that share this pInvar. Note that the conditional
       likelihood update flags for divisions have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
    
    /* However, you do need to update cijk flags if this is a covarion model */
    /* TO DO */
    
    return (NO_ERROR);
}


int Move_PopSize_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             isValidN, valIndex;
    MrBFlt          *valPtr, oldN, newN, tuning, minN, maxN, ran, oldLnPrior, newLnPrior, growth,
                    oldT, newT, clockRate;
    ModelParams     *mp;
    ModelInfo       *m;
    Tree            *t;

    /* get multiplier tuning parameter */
    tuning = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* get model settings */
    m = &modelSettings[param->relParts[0]];

    /* get minimum and maximum values for population size */
    if (param->paramId == POPSIZE_UNI)
        {
        minN = mp->popSizeUni[0];
        maxN = mp->popSizeUni[1];
        }
    else
        {
        minN = 0.00001;
        maxN = 10000000;
        }

    /* get pointer to value to be changed */
    valIndex = (int)(RandomNumber(seed) * param->nValues);
    valPtr = GetParamVals(param, chain, state[chain]) + valIndex;

    /* get old value of population size */
    oldN = *valPtr;

    /* get old prior for species tree coalescence */
    if (m->brlens->paramId == BRLENS_CLOCK_SPCOAL)
        {
        oldLnPrior = LnSpeciesTreeProb(chain);
        }
    
    /* change value for theta */
    ran = RandomNumber(seed);
    newN = oldN * exp(tuning * (ran - 0.5));
    
    /* check that new value is valid */
    isValidN = NO;
    do {
        if (newN < minN)
            newN = 2* minN - newN;
        else if (newN > maxN)
            newN = 2 * maxN - newN;
        else
            isValidN = YES;
        }
    while (isValidN == NO);

    /* copy new population size value back */
    (*valPtr) = newN;

    /* get proposal ratio */
    *lnProposalRatio = log (newN / oldN);
    
    /* get prior ratio */
    if (m->brlens->paramId == BRLENS_CLOCK_SPCOAL)
        {
        newLnPrior = LnSpeciesTreeProb(chain);
        }
    else
        {
        t = GetTree(modelSettings[param->relParts[0]].brlens,chain,state[chain]);
        m = &modelSettings[param->relParts[0]];
        clockRate = *GetParamVals(m->clockRate, chain, state[chain]);
        if (!strcmp(mp->ploidy, "Diploid"))
            clockRate *= 4.0;
        else if (!strcmp(mp->ploidy, "Zlinked"))
            clockRate *= 3.0;
        else
            clockRate *= 2.0;
        newT = oldN * clockRate;
        oldT = newN * clockRate;
        if (!strcmp(mp->growthPr, "Fixed"))
            growth = mp->growthFix;
        else
            growth = *(GetParamVals (m->growthRate, chain, state[chain]));
        if (LnCoalescencePriorPr (t, &oldLnPrior, oldT, growth) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for coalescent process\n", spacer);
            return (ERROR);
            }
        if (LnCoalescencePriorPr (t, &newLnPrior, newT, growth) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for coalescent process\n", spacer);
            return (ERROR);
            }
        }

    (*lnPriorRatio) = param->LnPriorRatio(newN, oldN, param->priorParams);
    (*lnPriorRatio) += newLnPrior - oldLnPrior;

    return (NO_ERROR);
}


/* Generalized lognormal move for positive real random variables */
int Move_PosRealLognormal (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i;
    MrBFlt      oldX, newX, minX, maxX, tuning, u, z;

    /* get tuning parameter */
    tuning = mvp[0];

    /* get minimum and maximum values for X */
    minX = param->min;
    maxX = param->max;

    /* get old value of X */
    newX = oldX = *GetParamVals(param, chain, state[chain]);

    /* change value of X */
    u = RandomNumber(seed);
    z = PointNormal(u);

    newX = exp (log(oldX) + z * tuning);
    
    /* check that new value is valid */
    if (newX < minX || newX > maxX) {
        abortMove = YES;
        return (NO_ERROR);
    }
    
    /* get proposal ratio */
    (*lnProposalRatio) = log (newX / oldX);
    
    /* get prior ratio */
    (*lnPriorRatio) = param->LnPriorRatio(newX, oldX, param->priorParams);
    
    /* copy new value back */
    (*GetParamVals(param, chain, state[chain])) = newX;

    /* Set update flags for tree nodes if relevant */
    if (param->affectsLikelihood == YES)
        {
        for (i=0; i<param->nRelParts; i++)
            TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        }

    return (NO_ERROR);
}


/* Generalized multiplier move for positive real random variables */
int Move_PosRealMultiplier (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, isValid;
    MrBFlt      oldX, newX, minX, maxX, tuning, ran, factor;

    /* get tuning parameter */
    tuning = mvp[0];

    /* get minimum and maximum values for X */
    minX = param->min;
    maxX = param->max;

    /* get old value of X */
    newX = oldX = *GetParamVals(param, chain, state[chain]);

    /* change value of X */
    ran = RandomNumber(seed);
    factor = exp(tuning * (ran - 0.5));
    newX = oldX * factor;
    
    /* check that new value is valid */
    isValid = NO;
    do
        {
        if (newX < minX)
            newX = minX * minX / newX;
        else if (newX > maxX)
            newX = maxX * maxX / newX;
        else
            isValid = YES;
        } while (isValid == NO);

    /* get proposal ratio */
    (*lnProposalRatio) = log (newX / oldX);
    
    /* get prior ratio */
    (*lnPriorRatio) = param->LnPriorRatio(newX, oldX, param->priorParams);
    
    /* copy new value back */
    *(GetParamVals(param, chain, state[chain])) = newX;

    /* Set update flags for tree nodes if relevant */
    if (param->affectsLikelihood == YES)
        {
        for (i=0; i<param->nRelParts; i++)
            TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        }

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_RateMult_Dir: Change rate multiplier using Dirichlet
|      proposal.
|
----------------------------------------------------------------*/
int Move_RateMult_Dir (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, nRates, isValid;
    MrBFlt      alphaPi, *value, *subValue, numSites, *alphaDir, x, y, sum,
                rate_pot, *dirParm, *oldRate, *newRate;

    /* allocate memory */
    dirParm = (MrBFlt *) SafeCalloc (3*numCurrentDivisions, sizeof(MrBFlt));
    oldRate = dirParm + numCurrentDivisions;
    newRate = dirParm + 2*numCurrentDivisions;

    /* get number of rates */
    nRates = param->nValues;

    /* get pointer to rates and number of uncompressed chars */
    value = GetParamVals(param, chain, state[chain]);
    subValue = GetParamSubVals(param, chain, state[chain]);

    /* get Dirichlet parameters */
    alphaDir = subValue + nRates;

    /* calculate old ratesum proportions */
    numSites = 0.0;
    for (i=0; i<nRates; i++)
        numSites += subValue[i];  /* numSites should be equal to the number of sites */
    for (i=0; i<nRates; i++)
        oldRate[i] = value[i] * subValue[i] / numSites;
    
    /* get alphaPi tuning parameter */
    alphaPi = mvp[0] * nRates;

    /* multiply old ratesum proportions with some large number to get new values close to the old ones */
    for (i=0; i<nRates; i++)
        dirParm[i] = oldRate[i] * alphaPi;
    
    /* get new values */
    DirichletRandomVariable (dirParm, newRate, nRates, seed);
    
    /* check new values. we rely on newRate be already normalized  */
    while (1)
        {
        sum = 0.0;
        rate_pot = 1.0;
        isValid=1;
        for (i=0; i<nRates; i++)
            {
            if (newRate[i] <= DIR_MIN)
                {
                if (newRate[i] < DIR_MIN)
                    {
                    newRate[i] = DIR_MIN;
                    isValid=0;
                    }
                rate_pot -= DIR_MIN;
                }
            else
                sum += newRate[i];
            }
        if (isValid==1) break;
        for (i=0; i<nRates; i++)
            {
            if (newRate[i]!=DIR_MIN)
                newRate[i] = rate_pot * newRate[i] / sum;
            }
        }

    /* calculate and copy new rate ratio values back */
    for (i=0; i<nRates; i++)
        value[i] = newRate[i] * (numSites / subValue[i]);
    
    /* get proposal ratio */
    sum = 0.0;
    for (i=0; i<nRates; i++)
        sum += newRate[i]*alphaPi;
    x = LnGamma(sum);
    for (i=0; i<nRates; i++)
        x -= LnGamma(newRate[i]*alphaPi);
    for (i=0; i<nRates; i++)
        x += (newRate[i]*alphaPi-1.0)*log(oldRate[i]);
    sum = 0.0;
    for (i=0; i<nRates; i++)
        sum += oldRate[i]*alphaPi;
    y = LnGamma(sum);
    for (i=0; i<nRates; i++)
        y -= LnGamma(oldRate[i]*alphaPi);
    for (i=0; i<nRates; i++)
        y += (oldRate[i]*alphaPi-1.0)*log(newRate[i]);
    (*lnProposalRatio) = x - y;

    /* get prior ratio */
    x = y = 0.0;
    for (i=0; i<nRates; i++)
        x += (alphaDir[i]-1.0)*log(newRate[i]);
    for (i=0; i<nRates; i++)
        y += (alphaDir[i]-1.0)*log(oldRate[i]);
    (*lnPriorRatio) = x - y;

    /* Set update flags for all partitions that share the rate multiplier. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* may need to hit update flag for cijks when you have a covarion model */
    for (i=0; i<param->nRelParts; i++)
        if (modelSettings[param->relParts[i]].nCijkParts > 1)
            modelSettings[param->relParts[i]].upDateCijk = YES;

    free (dirParm);

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_RateMult_Slider: Change rate multiplier using slider
|      proposal.
|
----------------------------------------------------------------*/
int Move_RateMult_Slider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, indexI, indexJ, nRates;
    MrBFlt      delta, *value, *subValue, sum, numSites, *alphaDir, x, y,
                oldRate[2], newRate[2], min, max;

    /* get number of rates */
    nRates = param->nValues;

    /* get pointer to rates and number of uncompressed chars */
    value = GetParamVals(param, chain, state[chain]);
    subValue = GetParamSubVals(param, chain, state[chain]);

    /* get Dirichlet prior parameters */
    alphaDir = subValue + nRates;

    /* randomly select two rates */
    indexI = (int) (RandomNumber(seed) * nRates);
    indexJ = (int) (RandomNumber(seed) * (nRates - 1));
    if (indexJ == indexI)
        indexJ = nRates - 1;

    /* calculate old ratesum proportions */
    numSites = 0.0;
    for (i=0; i<nRates; i++)
        numSites += subValue[i];  /* numSites should be equal to the number of sites */
    oldRate[0] = value[indexI] * subValue[indexI] / numSites;
    oldRate[1] = value[indexJ] * subValue[indexJ] / numSites;
    sum = oldRate[0] + oldRate[1];
    
    /* get delta tuning parameter */
    delta = mvp[0];

    /* reflect */
    min = DIR_MIN / sum;
    max = 1.0 - min;
    if (delta > max-min) /* we do it to avoid following long while loop in case if delta is high */
        delta = max-min;

    x = oldRate[0] / sum;
    y = x + delta * (RandomNumber(seed) - 0.5);
    while (y < min || y > max)
        {
        if (y < min)
            y = 2.0 * min - y;
        if (y > max)
            y = 2.0 * max - y;
        }
    
    /* set the new values */
    newRate[0] = y * sum;
    newRate[1] = sum - newRate[0];
    value[indexI] = newRate[0] * numSites / subValue[indexI];
    value[indexJ] = newRate[1] * numSites / subValue[indexJ];

    /* get proposal ratio */
    (*lnProposalRatio) = 0.0;

    /* get prior ratio */
    (*lnPriorRatio)  = (alphaDir[indexI]-1.0) * (log(newRate[0]) - log(oldRate[0]));
    (*lnPriorRatio) += (alphaDir[indexJ]-1.0) * (log(newRate[1]) - log(oldRate[1]));

    /* Set update flags for all partitions that share the rate multiplier. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* may need to hit update flag for cijks when you have a covarion model */
    for (i=0; i<param->nRelParts; i++)
        if (modelSettings[param->relParts[i]].nCijkParts > 1)
            modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_Revmat_Dir: Change rate matrix using Dirichlet proposal
|      mechanism.
|
----------------------------------------------------------------*/
int Move_Revmat_Dir (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change revMat using Dirichlet proposal */
    
    int             i, nRates,isValid;
    MrBFlt          oldRate[200], newRate[200], dirParm[200], *value, sum, x, y, rate_pot, *alphaDir, alphaPi;
    ModelParams     *mp;
    ModelInfo       *m;

    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];
    m  = &modelSettings[param->relParts[0]];

    /* get rates and nRates */
    value = GetParamVals(param, chain, state[chain]);
    nRates = param->nValues;

    /* get so called alpha_pi parameter and adjust for number of components */
    alphaPi = mvp[0] * nRates;

    /* get Dirichlet parameters */
    if (m->dataType == PROTEIN)
        alphaDir = mp->aaRevMatDir;
    else
        alphaDir = mp->revMatDir;

    /* copy old rates */
    for (i=0; i<nRates; i++)
        oldRate[i] = value[i];
    
    /* multiply old ratesum props with some large number to get new values close to the old ones */
    for (i=0; i<nRates; i++)
        dirParm[i] = oldRate[i] * alphaPi;
    
    /* get new values */
    DirichletRandomVariable (dirParm, newRate, nRates, seed);

    /* check new values. we rely on newRate be already normalized  */
    while (1)
        {
        sum = 0.0;
        rate_pot = 1.0;
        isValid=1;
        for (i=0; i<nRates; i++)
            {
            if (newRate[i] <= RATE_MIN)
                {
                if (newRate[i] < RATE_MIN)
                    {
                    newRate[i] = RATE_MIN;
                    isValid=0;
                    }
                rate_pot -= RATE_MIN;
                }
            else
                sum += newRate[i];
            }
        if (isValid==1) break;
        for (i=0; i<nRates; i++)
            {
            if (newRate[i]!=RATE_MIN)
                newRate[i] = rate_pot * newRate[i] / sum;
            }
        }

    /* copy new rate ratio values back */
    for (i=0; i<nRates; i++)
        value[i] = newRate[i];
    
    /* get proposal ratio */
    sum = 0.0;
    for (i=0; i<nRates; i++)
        sum += newRate[i]*alphaPi;
    x = LnGamma(sum);
    for (i=0; i<nRates; i++)
        x -= LnGamma(newRate[i]*alphaPi);
    for (i=0; i<nRates; i++)
        x += (newRate[i]*alphaPi-1.0)*log(oldRate[i]);
    sum = 0.0;
    for (i=0; i<nRates; i++)
        sum += oldRate[i]*alphaPi;
    y = LnGamma(sum);
    for (i=0; i<nRates; i++)
        y -= LnGamma(oldRate[i]*alphaPi);
    for (i=0; i<nRates; i++)
        y += (oldRate[i]*alphaPi-1.0)*log(newRate[i]);
    (*lnProposalRatio) = x - y;

    /* get prior ratio */
    x = y = 0.0;
    for (i=0; i<nRates; i++)
        x += (alphaDir[i]-1.0)*log(newRate[i]);
    for (i=0; i<nRates; i++)
        y += (alphaDir[i]-1.0)*log(oldRate[i]);
    (*lnPriorRatio) = x - y;

    /* Set update flags for all partitions that share this revmat. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* Set update flags for cijks for all affected partitions */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_Revmat_DirMix: Dirichlet proposal for REVMAT_MIX. From
|      Huelsenbeck et al. (2004), but note that the prior density
|      is different in that paper because they set the rate sum
|      to 6, not to 1.
|
----------------------------------------------------------------*/
int Move_Revmat_DirMix (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, j, k, isValid, *growthFxn, nRates, groupSize[6];
    MrBFlt      *value, dirParm[6], newRate[6], oldRate[6], alphaPi, symDir, sum, rate_pot, x, y;
    ModelParams *mp;

    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];

    /* get growthFunction and nRates */
    value     = GetParamVals (param, chain, state[chain]);
    growthFxn = GetParamIntVals (param, chain, state[chain]);
    nRates    = GetKFromGrowthFxn(growthFxn);

    /* we can't do anything if there is only one rate */
    if (nRates == 1)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* extract unique rates from value vector */
    for (i=0; i<nRates; i++)
        oldRate[i] = 0.0;
    for (i=0; i<6; i++)
        oldRate[growthFxn[i]] += value[i];

    /* get so called alpha_pi parameter and adjust for number of components */
    alphaPi = mvp[0] * nRates;

    /* get symmetric dirichlet parameter */
    symDir  = mp->revMatSymDir;

    /* multiply old ratesum props with some large number to get new values close to the old ones */
    for (i=0; i<nRates; i++)
        dirParm[i] = oldRate[i] * alphaPi;

    /* get new values */
    DirichletRandomVariable (dirParm, newRate, nRates, seed);

    /* check new values. we rely on newRate be already normalized  */
    while (1)
        {
        sum = 0.0;
        rate_pot = 1.0;
        isValid=1;
        for (i=0; i<nRates; i++)
            {
            if (newRate[i] <= RATE_MIN)
                {
                if (newRate[i] < RATE_MIN)
                    {
                    newRate[i] = RATE_MIN;
                    isValid=0;
                    }
                rate_pot -= RATE_MIN;
                }
            else
                sum += newRate[i];
            }
        if (isValid==1) break;
        for (i=0; i<nRates; i++)
            {
            if (newRate[i]!=RATE_MIN)
                newRate[i] = rate_pot * newRate[i] / sum;
            }
        }

    /* copy new unique rate ratio values back into the value array */
    for (i=0; i<nRates; i++)
        {
        k = 0;
        for (j=i; j<6; j++)
            {
            if (growthFxn[j] == i)
                k++;
            }
        for (j=i; j<6; j++)
            {
            if (growthFxn[j] == i)
                value[j] = newRate[i] / (MrBFlt) k;
            }
        }
    
    /* get proposal ratio */
    sum = 0.0;
    for (i=0; i<nRates; i++)
        sum += newRate[i]*alphaPi;
    x = LnGamma(sum);
    for (i=0; i<nRates; i++)
        x -= LnGamma(newRate[i]*alphaPi);
    for (i=0; i<nRates; i++)
        x += (newRate[i]*alphaPi-1.0)*log(oldRate[i]);
    sum = 0.0;
    for (i=0; i<nRates; i++)
        sum += oldRate[i]*alphaPi;
    y = LnGamma(sum);
    for (i=0; i<nRates; i++)
        y -= LnGamma(oldRate[i]*alphaPi);
    for (i=0; i<nRates; i++)
        y += (oldRate[i]*alphaPi-1.0)*log(newRate[i]);
    (*lnProposalRatio) = x - y;

    /* get group sizes, needed for prior ratio */
    for (i=0; i<nRates; i++)
        groupSize[i] = 0;
    for (i=0; i<6; i++)
        groupSize[growthFxn[i]]++;

    /* get prior ratio */
    x = y = 0.0;
    for (i=0; i<nRates; i++)
        x += (groupSize[i]*symDir-1.0)*log(newRate[i]);
    for (i=0; i<nRates; i++)
        y += (groupSize[i]*symDir-1.0)*log(oldRate[i]);
    (*lnPriorRatio) = x - y;

    /* Set update flags for all partitions that share this revmat. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* Set update flags for cijks for all affected partitions */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_Revmat_Slider: Change rate matrix using sliding window
|       move. Choose a pair of rates (e.g. r(A<>C), and r(A<>G)) at
|       random and denote them rA, and rB. Let oldProp = rA/(rA + rB)
|       and newProp = oldProp + delta(U - 0.5), where U is a uniform
|       random variable on the interval (0, 1] and delta is a tuning
|       parameter. Values that fall outside the boundaries are reflected
|       back in. Then set new_rA = newProp*(rA+rB) and new_rB =
|       (1-newProp)*(piA+piB). The Hastings ratio of this move is 1.0.
|
----------------------------------------------------------------*/
int Move_Revmat_Slider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, j, nRates;
    MrBFlt      delta, *newRate, *oldRate, *priorAlpha, x, y, sum, min, max;
    ModelParams *mp;
    ModelInfo   *m;

    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];
    m  = &modelSettings[param->relParts[0]];

    /* get Dirichlet parameters */
    if (m->dataType == PROTEIN)
        priorAlpha = mp->aaRevMatDir;
    else
        priorAlpha = mp->revMatDir;

    /* get the values we need */
    nRates = param->nValues;
    newRate = GetParamVals (param, chain, state[chain]);
    oldRate = GetParamVals (param, chain, state[chain] ^ 1);

    /* choose a pair to change */
    i = (int) (RandomNumber(seed) * nRates);
    j = (int) (RandomNumber(seed) * (nRates-1));
    if (i == j)
        j = nRates-1;
    
    /* find new proportion */
    sum = oldRate[i] + oldRate[j];

    /* get window size */
    delta = mvp[0];

    /* reflect */
    min = RATE_MIN / sum;
    max = 1.0 - min;
    if (delta > max-min) /* we do it to avoid following long while loop in case if delta is high */
        delta = max-min;

    x = oldRate[i] / sum;
    y = x + delta * (RandomNumber(seed) - 0.5);
    while (y < min || y > max)
        {
        if (y < min)
            y = 2.0 * min - y;
        if (y > max)
            y = 2.0 * max - y;
        }
    
    /* set the new values */
    newRate[i] = y * sum;
    newRate[j] = sum - newRate[i];

    /* get proposal ratio */
    (*lnProposalRatio) = 0.0;

    /* get prior ratio */
    (*lnPriorRatio)  = (priorAlpha[i]-1.0) * (log(newRate[i]) - log(oldRate[i]));
    (*lnPriorRatio) += (priorAlpha[j]-1.0) * (log(newRate[j]) - log(oldRate[j]));

    /* Set update for entire tree */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for many models we do want to update the cijk flag, as the transition
       probability matrices require diagonalizing the rate matrix. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_Revmat_SplitMerge1: Split or merge rates of rate matrix.
|      See Huelsenbeck et al. (2004). Note that the prior used
|      here is different from theirs. Also, a Beta proposal is
|      used instead of a uniform to propose new rate proportions.
|
----------------------------------------------------------------*/
int Move_Revmat_SplitMerge1 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, j, k, index_i, index_j, n_i, n_j, foundFirstI, foundFirstJ,
                *newGrowthFxn, *oldGrowthFxn, nOldRates, nNewRates, merge,
                groupSize[6], nCompositeRates;
    MrBFlt      R, R_i, R_j, *newValue, *oldValue, newRate[6], oldRate[6], symDir,
                prob_split, prob_merge, dirParm[2], rateProps[2], x, alphaPi;
    ModelParams *mp;

    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];

    /* get the values we need */
    oldValue     = GetParamVals(param, chain, state[chain] ^ 1);
    newValue     = GetParamVals(param, chain, state[chain]);
    oldGrowthFxn = GetParamIntVals(param, chain, state[chain] ^ 1);
    newGrowthFxn = GetParamIntVals (param, chain, state[chain]);
    nOldRates    = GetKFromGrowthFxn(oldGrowthFxn);
    symDir       = mp->revMatSymDir;
    alphaPi      = mvp[0];      /* tuning parameter alpha */

    /* get the old rates */
    for (i=0; i<nOldRates; i++)
        oldRate[i] = 0.0;
    for (i=0; i<6; i++)
        oldRate[oldGrowthFxn[i]] += oldValue[i];

    /* decide whether to split or merge */
    if (nOldRates == 1)
        merge = NO;
    else if (nOldRates == 6)
        merge = YES;
    else if (RandomNumber(seed) < 0.5)
        merge = YES;
    else
        merge = NO;

    /* now split or merge */
    R = R_i = R_j = 0.0;
    if (merge == YES)
        {
        /* merge two rates */
        nNewRates = nOldRates - 1;

        /* determine split and merge probs */
        if (nNewRates == 1)
            prob_split = 1.0;
        else
            prob_split = 0.5;
        if (nOldRates == 6)
            prob_merge = 1.0;
        else
            prob_merge = 0.5;

        /* select two rates randomly */
        index_i = (int) (RandomNumber(seed) * nOldRates);
        index_j = (int) (RandomNumber(seed) * (nOldRates - 1));
        if (index_j == index_i)
            index_j = nOldRates - 1;

        /* make sure index_i is lower index */
        if (index_i > index_j)
            {
            i = index_i;
            index_i = index_j;
            index_j = i;
            }

        /* find group sizes */
        n_i = n_j = 0;
        for (i=0; i<6; i++)
            {
            if (oldGrowthFxn[i] == index_i)
                n_i++;
            else if (oldGrowthFxn[i] == index_j)
                n_j++;
            }

        /* adjust growth function */
        for (i=0; i<6; i++)
            {
            if (oldGrowthFxn[i] == index_j)
                newGrowthFxn[i] = index_i;
            else if (oldGrowthFxn[i] > index_j)
                newGrowthFxn[i] = oldGrowthFxn[i] - 1;
            else
                newGrowthFxn[i] = oldGrowthFxn[i];
            }

        /* find the new rates */
        for (i=0; i<nNewRates; i++)
            {
            if (i == index_i)
                newRate[i] = (oldRate[index_i] + oldRate[index_j]);
            else if (i < index_j)
                newRate[i] = oldRate[i];
            else if (i >= index_j)
                newRate[i] = oldRate[i+1];
            }

        /* copy new unique rate values back into the value array */
        for (i=0; i<nNewRates; i++)
            {
            k = 0;
            for (j=i; j<6; j++)
                {
                if (newGrowthFxn[j] == i)
                    k++;
                }
            for (j=i; j<6; j++)
                {
                if (newGrowthFxn[j] == i)
                    newValue[j] = newRate[i] / (MrBFlt) k;
                }
            }

        /* get the new and old rates (sum over parts) */
        R_i = oldRate[index_i];
        R_j = oldRate[index_j];
        R   = R_i + R_j;

        /* check group sizes after merge (before split in back move) */
        for (i=0; i<nNewRates; i++)
            groupSize[i] = 0;
        for (i=0; i<6; i++)
            groupSize[newGrowthFxn[i]]++;
        nCompositeRates = 0;
        for (i=0; i<nNewRates; i++)
            {
            if (groupSize[i] > 1)
                nCompositeRates++;
            }

        /* calculate prior ratio (different in the paper) */
        (*lnPriorRatio) = LnGamma(n_i * symDir) + LnGamma(n_j * symDir) - LnGamma ((n_i + n_j) * symDir);
        (*lnPriorRatio) += ((n_i + n_j) * symDir - 1.0) * log(R) - (n_i * symDir - 1.0) * log(R_i) - (n_j * symDir - 1.0) * log(R_j);

        /* calculate proposal ratio */
        (*lnProposalRatio) = log ((prob_split / prob_merge) * ((nOldRates * (nOldRates - 1)) / (2.0 * nCompositeRates)) * (1.0 / ((pow(2, n_i + n_j - 1) - 1))));

        /* adjust for Beta proposal in back move */
        dirParm[0] = alphaPi * n_i;
        dirParm[1] = alphaPi * n_j;
        rateProps[0] = R_i / R;
        rateProps[1] = R_j / R;
        x  = LnGamma(dirParm[0] + dirParm[1]);
        x -= LnGamma(dirParm[0]);
        x -= LnGamma(dirParm[1]);
        x += (dirParm[0] - 1.0) * log(rateProps[0]);
        x += (dirParm[1] - 1.0) * log(rateProps[1]);
        (*lnProposalRatio) += x;

        /* Jacobian for the rate proportion */
        (*lnProposalRatio) -= log(R);
        }
    else
        {
        /* split two rates */
        nNewRates = nOldRates + 1;

        /* determine split and merge probs */
        if (nNewRates == 6)
            prob_merge = 1.0;
        else
            prob_merge = 0.5;
        if (nOldRates == 1)
            prob_split = 1.0;
        else
            prob_split = 0.5;

        /* check group sizes before split */
        for (i=0; i<nOldRates; i++)
            groupSize[i] = 0;
        for (i=0; i<6; i++)
            groupSize[oldGrowthFxn[i]]++;
        nCompositeRates = 0;
        for (i=0; i<nOldRates; i++)
            {
            if (groupSize[i] > 1)
                nCompositeRates++;
            }

        /* randomly select a rate with two or more components to split */
        k = (int) (RandomNumber(seed) * nCompositeRates);

        for (i=j=0; i<nOldRates; i++)
            {
            if (groupSize[i] > 1)
                {
                if (k == j)
                    break;
                j++;
                }
            }
        assert (i < nOldRates && groupSize[i] > 1);
        index_i = i;

        /* adjust growth function */
        do {
            foundFirstI = foundFirstJ = NO;
            k = 0;
            index_j = -1;
            for (i=0; i<6; i++)
                {
                if (oldGrowthFxn[i] == index_i)
                    {
                    if (foundFirstI == NO)
                        {
                        newGrowthFxn[i] = index_i;
                        foundFirstI = YES;
                        }
                    else
                        {
                        if (RandomNumber(seed) < 0.5)
                            {
                            if (foundFirstJ == NO)
                                {
                                foundFirstJ = YES;
                                index_j = k + 1;    /* one more than previous max */
                                newGrowthFxn[i] = index_j;
                                }
                            else
                                {
                                newGrowthFxn[i] = index_j;
                                }
                            }
                        else
                            newGrowthFxn[i] = index_i;
                        }
                    }
                else if (foundFirstJ == YES && oldGrowthFxn[i] >= index_j)
                    newGrowthFxn[i] = oldGrowthFxn[i] + 1;
                else
                    newGrowthFxn[i] = oldGrowthFxn[i];
                if (foundFirstJ == NO && oldGrowthFxn[i] > k)
                    k = oldGrowthFxn[i];
                }
            } while (foundFirstJ == NO);

        /* find group sizes */
        n_i = n_j = 0;
        for (i=0; i<6; i++)
            {
            if (newGrowthFxn[i] == index_i)
                n_i++;
            else if (newGrowthFxn[i] == index_j)
                n_j++;
            }

        /* find old rate */
        R = oldRate[index_i];

        /* propose new rates */
        dirParm[0] = alphaPi * n_i;
        dirParm[1] = alphaPi * n_j;

        DirichletRandomVariable(dirParm, rateProps, 2, seed);
        R_i = rateProps[0] * R;
        R_j = rateProps[1] * R;

        if (R_i/n_i < RATE_MIN)
            {
            R_i = RATE_MIN*n_i;
            rateProps[0] = R_i/R;
            rateProps[1] = 1-rateProps[0];
            R_j = rateProps[1] * R;
            assert (R_j/n_j < RATE_MIN);
            }
        else if (R_j/n_j < RATE_MIN)
            {
            R_j = RATE_MIN*n_j;
            rateProps[1] = R_j/R;
            rateProps[0] = 1-rateProps[1];
            R_i = rateProps[0] * R;
            assert (R_i/n_i < RATE_MIN);
            }

        /* set the new rates */
        for (i=0; i<nNewRates; i++)
            {
            if (i == index_i)
                newRate[i] = R_i;
            else if (i == index_j)
                newRate[i] = R_j;
            else if (i > index_j)
                newRate[i] = oldRate[i-1];
            else
                newRate[i] = oldRate[i];
            }

        /* copy new unique rate values back into the value array */
        for (i=0; i<nNewRates; i++)
            {
            k = 0;
            for (j=i; j<6; j++)
                {
                if (newGrowthFxn[j] == i)
                    k++;
                }
            for (j=i; j<6; j++)
                {
                if (newGrowthFxn[j] == i)
                    newValue[j] = newRate[i] / (MrBFlt) k;
                }
            }

        /* calculate prior ratio (different in the paper) */
        (*lnPriorRatio) = LnGamma((n_i + n_j) * symDir) - LnGamma(n_i * symDir) - LnGamma(n_j * symDir);
        (*lnPriorRatio) += (n_i * symDir - 1.0) * log(R_i) + (n_j * symDir - 1.0) * log(R_j) - ((n_i + n_j) * symDir - 1.0) * log(R);;

        /* calculate proposal ratio */
        (*lnProposalRatio) = log ((prob_merge / prob_split) * ((2.0 * nCompositeRates) / (nNewRates * (nNewRates - 1))) * ((pow(2, n_i + n_j - 1) - 1)));
        
        /* adjust for Beta proposal */
        x  = LnGamma(dirParm[0] + dirParm[1]);
        x -= LnGamma(dirParm[0]);
        x -= LnGamma(dirParm[1]);
        x += (dirParm[0] - 1.0) * log(rateProps[0]);
        x += (dirParm[1] - 1.0) * log(rateProps[1]);
        (*lnProposalRatio) -= x;

        /* Jacobian for rate proportion */
        (*lnProposalRatio) += log (R);
        }

#if defined (DEBUG_SPLITMERGE)
    if (*lnPriorRatio != *lnPriorRatio)
        {
        printf ("prob_merge=%f prob_split=%f nCompositeRates=%d nOldRates=%d nNewRates=%d\n", prob_merge, prob_split, nCompositeRates, nOldRates, nNewRates);
        printf ("merge=%s n_i=%d n_j=%d rateProps[0]=%f R=%f R_i=%f R_j=%f\n", merge == NO ? "NO" : "YES", n_i, n_j, rateProps[0], R, R_i, R_j);
        printf ("Old rates={%f,%f,%f,%f,%f,%f}\n", oldValue[0], oldValue[1], oldValue[2], oldValue[3], oldValue[4], oldValue[5]);
        printf ("Old growth fxn={%d,%d,%d,%d,%d,%d}\n", oldGrowthFxn[0], oldGrowthFxn[1], oldGrowthFxn[2], oldGrowthFxn[3], oldGrowthFxn[4], oldGrowthFxn[5]);
        printf ("New rates={%f,%f,%f,%f,%f,%f}\n", newValue[0], newValue[1], newValue[2], newValue[3], newValue[4], newValue[5]);
        printf ("New growth fxn={%d,%d,%d,%d,%d,%d}\n", newGrowthFxn[0], newGrowthFxn[1], newGrowthFxn[2], newGrowthFxn[3], newGrowthFxn[4], newGrowthFxn[5]);
        printf ("lnPriorRatio=%f  lnProposalRatio=%f\n", *lnPriorRatio, *lnProposalRatio);
        getchar();
        }
#endif

    /* Set update flags for all partitions that share this revmat. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* Set update flags for cijks for all affected partitions */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_Revmat_SplitMerge2: Componentwise split or merge move.
|
----------------------------------------------------------------*/
int Move_Revmat_SplitMerge2 (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, k, n_i, n_j, index_i, index_j, groupIndex_i, groupIndex_j,
                *newGrowthFxn, *oldGrowthFxn;
    MrBFlt      R_i, R_j, r_j, alphaPi, *newValue, *oldValue, symDir,
                dirParm[2], rateProps[2], x;
    ModelParams *mp;

    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];

    /* get the values we need */
    oldValue     = GetParamVals(param, chain, state[chain] ^ 1);
    newValue     = GetParamVals(param, chain, state[chain]);
    oldGrowthFxn = GetParamIntVals(param, chain, state[chain] ^ 1);
    newGrowthFxn = GetParamIntVals (param, chain, state[chain]);
    symDir       = mp->revMatSymDir;
    alphaPi      = mvp[0];      /* tuning parameter */

    /* pick two component rates at random without replacement */
    index_i = (int) (RandomNumber(seed) * 6);
    index_j = (int) (RandomNumber(seed) * 5);
    if (index_j == index_i)
        index_j = 5;
    groupIndex_i = oldGrowthFxn[index_i];
    groupIndex_j = oldGrowthFxn[index_j];

    if (oldGrowthFxn[index_i] != oldGrowthFxn[index_j])
        {
        /* the rates are different, so merge them */

        /* calculate n_i, n_j, R_i and R_j before merge */
        n_i = n_j = 0;
        R_i = R_j = 0.0;
        for (i=0; i<6; i++)
            {
            if (oldGrowthFxn[i] == groupIndex_i)
                {
                n_i++;
                R_i += oldValue[i];
                }
            if (oldGrowthFxn[i] == groupIndex_j)
                {
                n_j++;
                R_j += oldValue[i];
                }
            }

        /* merge component rates by adding j to i */
        newGrowthFxn[index_j] = oldGrowthFxn[index_i];

        /* select a new rate for r_j */
        if (n_j == 1)
            r_j = oldValue[index_j];
        else
            {
            dirParm[0] = alphaPi * 1;
            dirParm[1] = alphaPi * (n_j - 1);

            DirichletRandomVariable(dirParm, rateProps, 2, seed);
            r_j = rateProps[0] * R_j;

            if (R_j - r_j < RATE_MIN)
                {
                r_j = R_j - RATE_MIN;
                rateProps[0] = r_j/R_j;
                rateProps[1] = 1 - rateProps[0];
                }
            }

        /* update new growth function */
        UpdateGrowthFxn(newGrowthFxn);

        /* we divide R_i + r_j equally among components of merged group,
           and R_j - r_j equally among split group */
        for (i=0; i<6; i++)
            {
            if (oldGrowthFxn[i] == oldGrowthFxn[index_i] || i == index_j)
                newValue[i] = (R_i + r_j) / (MrBFlt)(n_i + 1);
            else if (oldGrowthFxn[i] == oldGrowthFxn[index_j])
                newValue[i] = (R_j - r_j) / (MrBFlt)(n_j - 1);
            else
                newValue[i] = oldValue[i];
            }

        /* calculate prior ratio */
        if (n_j > 1)
            {
            /* no category disappeared */
            (*lnPriorRatio) += LnGamma (n_i * symDir) + LnGamma(n_j * symDir);
            (*lnPriorRatio) -= LnGamma((n_i +1)* symDir) + LnGamma((n_j-1) * symDir);
            (*lnPriorRatio) += ((n_i + 1) * symDir - 1.0) * log(R_i + r_j) + ((n_j - 1) * symDir - 1.0) * log(R_j - r_j);
            (*lnPriorRatio) -= (n_i * symDir - 1.0) * log(R_i) + (n_j * symDir - 1.0) * log(R_j);
            }
        else
            {
            /* j category disappeared */
            (*lnPriorRatio) += LnGamma (n_i * symDir) + LnGamma(n_j * symDir);
            (*lnPriorRatio) -= LnGamma((n_i +1)* symDir);
            (*lnPriorRatio) += ((n_i + 1) * symDir - 1.0) * log(R_i + r_j);
            (*lnPriorRatio) -= (n_i * symDir - 1.0) * log(R_i) + (n_j * symDir - 1.0) * log(R_j);
            }

        /* calculate proposal ratio; this is the probability of choosing the right category for rate j when splitting */
        k = GetKFromGrowthFxn(newGrowthFxn);
        (*lnProposalRatio) = log (1.0 / k);

        /* adjust for Beta proposal in back move */
        dirParm[0] = alphaPi * 1;
        dirParm[1] = alphaPi * n_i;
        rateProps[0] = r_j / (R_i + r_j);
        rateProps[1] = 1.0 - rateProps[0];
        x  = LnGamma(dirParm[0] + dirParm[1]);
        x -= LnGamma(dirParm[0]);
        x -= LnGamma(dirParm[1]);
        x += (dirParm[0] - 1.0) * log(rateProps[0]);
        x += (dirParm[1] - 1.0) * log(rateProps[1]);
        (*lnProposalRatio) += x;

        /* adjust for Beta proposal in forward move */
        if (n_j > 1)
            {
            dirParm[0] = alphaPi * 1;
            dirParm[1] = alphaPi * n_j;
            rateProps[0] = r_j / R_j;
            rateProps[1] = 1.0 - rateProps[0];
            x  = LnGamma(dirParm[0] + dirParm[1]);
            x -= LnGamma(dirParm[0]);
            x -= LnGamma(dirParm[1]);
            x += (dirParm[0] - 1.0) * log(rateProps[0]);
            x += (dirParm[1] - 1.0) * log(rateProps[1]);
            (*lnProposalRatio) -= x;
            }

        /* Jacobian */
        (*lnProposalRatio) -= log (R_i + r_j);
        if (n_j > 1)
            (*lnProposalRatio) += log (R_j);
        }
    else
        {
        /* split component rates because they are the same */

        /* split component rates by selecting new group for j from (0,K), with j starting a new group if index becomes the same */
        k = GetKFromGrowthFxn(oldGrowthFxn);
        newGrowthFxn[index_j] = (int) (RandomNumber(seed) * k);
        if (newGrowthFxn[index_j] == oldGrowthFxn[index_j])
            newGrowthFxn[index_j] = k + 1;

        /* update growth function and group indices */
        UpdateGrowthFxn(newGrowthFxn);
        groupIndex_i = newGrowthFxn[index_i];
        groupIndex_j = newGrowthFxn[index_j];

        /* calculate n_i, n_j, R_i and R_j after split */
        n_i = n_j = 0;
        R_i = R_j = 0.0;
        for (i=0; i<6; i++)
            {
            if (i == index_j)
                {
                R_i += oldValue[i];
                n_i++;
                }
            else if (newGrowthFxn[i] == groupIndex_i)
                {
                n_i++;
                R_i += oldValue[i];
                }
            else if (newGrowthFxn[i] == groupIndex_j)
                {
                n_j++;
                R_j += oldValue[i];
                }
            }

        /* select a new rate for r_j */
        dirParm[0] = alphaPi * 1;
        dirParm[1] = alphaPi * (n_i - 1);

        DirichletRandomVariable(dirParm, rateProps, 2, seed);
        r_j = rateProps[0] * R_i;

        if (R_i-r_j < RATE_MIN)
            {
            r_j = R_i - RATE_MIN;
            rateProps[0] = r_j/R_i;
            rateProps[1] = 1 - rateProps[0];
            }

        /* update n_i, n_j, R_i and R_j after split */
        n_i -= 1;
        n_j += 1;
        R_i -= r_j;
        R_j += r_j;

        /* we divide R_i equally among remaining components of split group,
           and R_j equally among new or expanded group */
        for (i=0; i<6; i++)
            {
            if (newGrowthFxn[i] == groupIndex_i)
                newValue[i] = R_i / (MrBFlt)(n_i);
            else if (newGrowthFxn[i] == groupIndex_j)
                newValue[i] = R_j / (MrBFlt)(n_j);
            else
                newValue[i] = oldValue[i];
            }

        /* calculate prior ratio */
        if (n_j > 1)
            {
            /* no new category created by split */
            (*lnPriorRatio) += LnGamma((n_i +1)* symDir) + LnGamma((n_j-1) * symDir);
            (*lnPriorRatio) -= LnGamma (n_i * symDir) + LnGamma(n_j * symDir);
            (*lnPriorRatio) += (n_i * symDir - 1.0) * log(R_i) + (n_j * symDir - 1.0) * log(R_j);
            (*lnPriorRatio) -= ((n_i + 1) * symDir - 1.0) * log(R_i + r_j) + ((n_j - 1) * symDir - 1.0) * log(R_j - r_j);
            }
        else
            {
            /* new category created by split */
            (*lnPriorRatio) += LnGamma((n_i +1)* symDir);
            (*lnPriorRatio) -= LnGamma (n_i * symDir) + LnGamma(n_j * symDir);
            (*lnPriorRatio) += (n_i * symDir - 1.0) * log(R_i) + (n_j * symDir - 1.0) * log(R_j);
            (*lnPriorRatio) -= ((n_i + 1) * symDir - 1.0) * log(R_i + r_j);
            }

        /* calculate proposal ratio; this is one over the probability of choosing the right category for rate j when splitting */
        k = GetKFromGrowthFxn(oldGrowthFxn);
        (*lnProposalRatio) = log (k);

        /* adjust for Beta proposal in back move */
        if (n_j > 1)
            {
            dirParm[0] = alphaPi * 1;
            dirParm[1] = alphaPi * (n_j - 1);
            rateProps[0] = r_j / R_j;
            rateProps[1] = 1.0 - rateProps[0];
            x  = LnGamma(dirParm[0] + dirParm[1]);
            x -= LnGamma(dirParm[0]);
            x -= LnGamma(dirParm[1]);
            x += (dirParm[0] - 1.0) * log(rateProps[0]);
            x += (dirParm[1] - 1.0) * log(rateProps[1]);
            (*lnProposalRatio) += x;
            }

        /* adjust for Beta proposal in forward move */
        dirParm[0] = alphaPi * 1;
        dirParm[1] = alphaPi * n_i;
        rateProps[0] = r_j / (R_i + r_j);
        rateProps[1] = 1.0 - rateProps[0];
        x  = LnGamma(dirParm[0] + dirParm[1]);
        x -= LnGamma(dirParm[0]);
        x -= LnGamma(dirParm[1]);
        x += (dirParm[0] - 1.0) * log(rateProps[0]);
        x += (dirParm[1] - 1.0) * log(rateProps[1]);
        (*lnProposalRatio) -= x;

        /* Jacobian */
        (*lnProposalRatio) += log (R_i + r_j);
        if (n_j > 1)
            (*lnProposalRatio) -= log (R_j);
        }

    /* Set update flags for all partitions that share this revmat. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* Set update flags for cijks for all affected partitions */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


int Move_Speciation (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change speciation rate using sliding window */
    
    int         isLPriorExp, isValidL, valIndex;
    MrBFlt      *valPtr, oldL, newL, minL, maxL, lambdaExp=0.0, *sR, *eR, sF, *fR, oldLnPrior, newLnPrior,
                window, clockRate;
    char        *sS;
    ModelParams *mp;
    ModelInfo   *m;
    Tree        *t;

    /* get size of window, centered on current value */
    window = mvp[0];

    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get minimum and maximum values */
    if (param->paramId == SPECRATE_UNI)
        {
        minL = mp->speciationUni[0];
        maxL = mp->speciationUni[1];
        isLPriorExp = NO;
        }
    else
        {
        minL = 0.000001;
        maxL = 1000.0;
        lambdaExp = mp->speciationExp;
        isLPriorExp = YES;
        }

    /* get pointer to value to be changed */
    valIndex = (int)(RandomNumber(seed) * param->nValues);
    valPtr = GetParamVals(param, chain, state[chain]) + valIndex;
    
    /* get old value */
    oldL = *valPtr;

    /* change value */
    if (maxL-minL < window)
        window = maxL-minL;
    newL = oldL + window * (RandomNumber(seed) - 0.5);
    
    /* check that new value is valid */
    isValidL = NO;
    do  {
        if (newL < minL)
            newL = 2 * minL - newL;
        else if (newL > maxL)
            newL = 2 * maxL - newL;
        else
            isValidL = YES;
        } while (isValidL == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* calculate prior ratio */
    t  = GetTree(modelSettings[param->relParts[0]].brlens,chain,state[chain]);
    sR = GetParamVals (param, chain, state[chain]);
    eR = GetParamVals (m->extinctionRates, chain, state[chain]);
    sF = mp->sampleProb;
    sS = mp->sampleStrat;
    clockRate = *GetParamVals (m->clockRate, chain, state[chain]);
    
    if (!strcmp(mp->clockPr,"Birthdeath"))
        {
        if (LnBirthDeathPriorPr (t, clockRate, &oldLnPrior, *sR, *eR, sS, sF) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
            return (ERROR);
            }
        *valPtr = newL;  // update with new value
        if (LnBirthDeathPriorPr (t, clockRate, &newLnPrior, *sR, *eR, sS, sF) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
            return (ERROR);
            }
        }
    else if (!strcmp(mp->clockPr,"Fossilization"))
        {
        fR = GetParamVals (m->fossilizationRates, chain, state[chain]);
        if (LnFossilizationPriorPr (t, clockRate, &oldLnPrior, sR, eR, sF, fR, sS) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
            return (ERROR);
            }
        *valPtr = newL;  // update with new value
        // for (i=0; i<param->nValues; i++)  *(GetParamVals(param, chain, state[chain]) + i) = newL;
        if (LnFossilizationPriorPr (t, clockRate, &newLnPrior, sR, eR, sF, fR, sS) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
            return (ERROR);
            }
        }
    else {
        MrBayesPrint ("%s   Move_Speciation not applicable\n", spacer);
        return (ERROR);
        }

    if (isLPriorExp == NO)
        *lnPriorRatio = newLnPrior - oldLnPrior;
    else
        *lnPriorRatio = -lambdaExp * (newL - oldL) + (newLnPrior - oldLnPrior);
    
    return (NO_ERROR);
}


int Move_Speciation_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change speciation rate using multiplier */
    
    int         isLPriorExp, isValidL, valIndex;
    MrBFlt      *valPtr, oldL, newL, minL, maxL, lambdaExp=0.0, *sR, *eR, sF, *fR, oldLnPrior, newLnPrior,
                tuning, clockRate;
    char        *sS;
    ModelParams *mp;
    ModelInfo   *m;
    Tree        *t;

    /* get tuning parameter */
    tuning = mvp[0];

    /* get model params and settings */
    mp = &modelParams[param->relParts[0]];
    m = &modelSettings[param->relParts[0]];
    
    /* get minimum and maximum values */
    if (param->paramId == SPECRATE_UNI)
        {
        minL = mp->speciationUni[0];
        maxL = mp->speciationUni[1];
        isLPriorExp = NO;
        }
    else
        {
        minL = 0.000001;
        maxL = 1000.0;
        lambdaExp = mp->speciationExp;
        isLPriorExp = YES;
        }

    /* get pointer to value to be changed */
    valIndex = (int)(RandomNumber(seed) * param->nValues);
    valPtr = GetParamVals(param, chain, state[chain]) + valIndex;
    
    /* get old value */
    oldL = *valPtr;

    /* change value */
    newL = oldL * exp(tuning * (RandomNumber(seed) - 0.5));
    
    /* check that new value is valid */
    isValidL = NO;
    do  {
        if (newL < minL)
            newL = minL * minL / newL;
        else if (newL > maxL)
            newL = maxL * maxL / newL;
        else
            isValidL = YES;
        } while (isValidL == NO);

    /* get proposal ratio */
    *lnProposalRatio = log (newL / oldL);
    
    /* calculate prior ratio */
    t  = GetTree(modelSettings[param->relParts[0]].brlens,chain,state[chain]);
    sR = GetParamVals (param, chain, state[chain]);
    eR = GetParamVals (m->extinctionRates, chain, state[chain]);
    sF = mp->sampleProb;
    sS = mp->sampleStrat;
    clockRate = *GetParamVals(m->clockRate, chain, state[chain]);
    
    if (!strcmp(mp->clockPr,"Birthdeath"))
        {
        if (LnBirthDeathPriorPr (t, clockRate, &oldLnPrior, *sR, *eR, sS, sF) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
            return (ERROR);
            }
        *valPtr = newL;  // update with new value
        if (LnBirthDeathPriorPr (t, clockRate, &newLnPrior, *sR, *eR, sS, sF) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for birth-death process\n", spacer);
            return (ERROR);
            }
        }
    else if (!strcmp(mp->clockPr,"Fossilization"))
        {
        fR = GetParamVals (m->fossilizationRates, chain, state[chain]);
        if (LnFossilizationPriorPr (t, clockRate, &oldLnPrior, sR, eR, sF, fR, sS) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
            return (ERROR);
            }
        *valPtr = newL;  // update with new value
        // for (i=0; i<param->nValues; i++)  *(GetParamVals(param, chain, state[chain]) + i) = newL;
        if (LnFossilizationPriorPr (t, clockRate, &newLnPrior, sR, eR, sF, fR, sS) == ERROR)
            {
            MrBayesPrint ("%s   Problem calculating prior for fossilized birth-death process\n", spacer);
            return (ERROR);
            }
        }
    else {
        MrBayesPrint ("%s   Move_Speciation_M not applicable\n", spacer);
        return (ERROR);
        }

    if (isLPriorExp == NO)
        *lnPriorRatio = newLnPrior - oldLnPrior;
    else
        *lnPriorRatio = -lambdaExp * (newL - oldL) + (newLnPrior - oldLnPrior);
    
    return (NO_ERROR);
}


int Move_Statefreqs (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change pi */
    int         i, nStates, isValid;
    MrBFlt      dirichletParameters[64], *newPi, *oldPi, *priorAlpha, sum, alphaPi, x, y;

    /* get the values we need */
    nStates = param->nSubValues;
    priorAlpha = GetParamVals(param, chain, state[chain]);
    newPi = GetParamSubVals (param, chain, state[chain]);
    oldPi = GetParamSubVals (param, chain, state[chain] ^ 1);
        
    /* tuning parameter */
    alphaPi = mvp[0]*nStates;

    /* multiply old values with some large number to get new values close to the old ones */
    for (i=0; i<nStates; i++)
        dirichletParameters[i] = oldPi[i] * alphaPi;

    do  {
        DirichletRandomVariable (dirichletParameters, newPi, nStates, seed);
        isValid = YES;
        for (i=0; i<nStates; i++)
            {
            if (newPi[i] < PI_MIN)
                {
                isValid = NO;
                break;
                }
            }
        } while (isValid == NO);

    /* get proposal ratio */
    sum = 0.0;
    for (i=0; i<nStates; i++)
        sum += newPi[i]*alphaPi;
    x = LnGamma(sum);
    for (i=0; i<nStates; i++)
        x -= LnGamma(newPi[i]*alphaPi);
    for (i=0; i<nStates; i++)
        x += (newPi[i]*alphaPi-1.0)*log(oldPi[i]);
    sum = 0.0;
    for (i=0; i<nStates; i++)
        sum += oldPi[i]*alphaPi;
    y = LnGamma(sum);
    for (i=0; i<nStates; i++)
        y -= LnGamma(oldPi[i]*alphaPi);
    for (i=0; i<nStates; i++)
        y += (oldPi[i]*alphaPi-1.0)*log(newPi[i]);
    (*lnProposalRatio) = x - y;

    /* get prior ratio */
    y = x = 0.0;                    /* the Gamma part of the prior is the same */
    for (i=0; i<nStates; i++)
        x += (priorAlpha[i]-1.0)*log(newPi[i]);
    for (i=0; i<nStates; i++)
        y += (priorAlpha[i]-1.0)*log(oldPi[i]);
    (*lnPriorRatio) = x - y;
        
    /* Touch the entire tree */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for many models we do want to update the cijk flag, as the transition
       probability matrices require diagonalizing the rate matrix. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/*----------------------------------------------------------------
|
|   Move_Statefreqs_Slider: Change state frequencies using Slider proposal
|       mechanism.
|       Choose pairs of the parameter values (e.g. pi(A), and pi(G)) at
|       random and denote them piA, and piB. Let oldProp = piA/(piA + piB)
|       and newProp = oldProp + delta(U - 0.5), where U is a uniform random variable
|       on the interval (0, 1] and delta is a tuning parameter. Values
|       that fall outside the boundaries are reflected back in. Then
|       set newPiA = newProp*(piA+piB) and newPiB = (1-newProp)*(piA+piB).
|       The Hastings ratio of this move is 1.0.
|
----------------------------------------------------------------*/
int Move_Statefreqs_Slider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, j, nStates, isValid;
    MrBFlt      delta, *newPi, *oldPi, *priorAlpha, x, y, sum, min, max;

    /* get the values we need */
    nStates = param->nSubValues;
    priorAlpha = GetParamVals(param, chain, state[chain]);
    newPi = GetParamSubVals (param, chain, state[chain]);
    oldPi = GetParamSubVals (param, chain, state[chain] ^ 1);

    /* get window size */
    delta = mvp[0];

    /* choose a pair to change */
    i = (int) (RandomNumber(seed) * nStates);
    j = (int) (RandomNumber(seed) * (nStates-1));
    if (i == j)
        j = nStates-1;
    
    /* find new proportion */
    sum = oldPi[i] + oldPi[j];

    /* reflect */
    isValid = NO;
    min = PI_MIN / sum;
    max = 1.0 - min;

    x   = oldPi[i] / sum;
    if (delta > max-min) /* we do it to avoid following long while loop in case if delta is high */
        {
        delta = max-min;
        }
    y = x + delta * (RandomNumber(seed) - 0.5);

    do {
        if (y < min)
            y = 2.0 * min - y;
        else if (y > max)
            y = 2.0 * max - y;
        else
            isValid = YES;
        } while (isValid == NO);

    /* set the new values */
    newPi[i] = y * sum;
    newPi[j] = sum - newPi[i];

    /* get proposal ratio */
    *lnProposalRatio = 0.0;

    /* get prior ratio */
    /* (the Gamma part of the prior is the same) */
    x = (priorAlpha[i]-1.0)*log(newPi[i]);
    x += (priorAlpha[j]-1.0)*log(newPi[j]);
    y = (priorAlpha[i]-1.0)*log(oldPi[i]);
    y += (priorAlpha[j]-1.0)*log(oldPi[j]);
    (*lnPriorRatio) = x - y;

    /* Set update for entire tree */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for many models we do want to update the cijk flag, as the transition
       probability matrices require diagonalizing the rate matrix. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


int Move_StatefreqsSymDirMultistate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change state freqs of multistate characters */
    /* ideally, we would let the likelihood calculator deal with only the affected character
       but we do not have the mechanism for doing that in the current version of mrbayes, so
       take the hit of updating all chars of the morph partition(s). */
    int     i, nStates, charIndex;
    MrBFlt  dirichletParameters[10], symDirAlphai, *newPi, *oldPi, sum, alphaPi, x, y;
    Model   *mp;

    /* tuning parameters */
    alphaPi = mvp[0];

    /* get model paramaters */
    mp = &modelParams[param->relParts[0]];

    /* select one character at random */
    charIndex = (int) (RandomNumber(seed) * param->nSympi);
    
    /* get the values we need */
    symDirAlphai = *GetParamVals(param, chain, state[chain]);
    newPi = GetParamStdStateFreqs (param, chain, state[chain]);
    oldPi = GetParamStdStateFreqs (param, chain, state[chain] ^ 1);
    newPi += 2 * mp->numBetaCats;
    oldPi += 2 * mp->numBetaCats;
    for (i=0; i<charIndex; i++)
        {
        oldPi += param->sympinStates[i];
        newPi += param->sympinStates[i];
        }
    nStates = param->sympinStates[charIndex];
    
    /* multiply old values with some large number to get new values close to the old ones */
    for (i=0; i<nStates; i++)
        dirichletParameters[i] = oldPi[i] * alphaPi;

    DirichletRandomVariable (dirichletParameters, newPi, nStates, seed);

    sum = 0.0;
    for (i=0; i<nStates; i++)
        {
        if (newPi[i] < 0.0001)
            newPi[i] = 0.0001;
        sum += newPi[i];
        }
    for (i=0; i<nStates; i++)
        newPi[i] /= sum;

    /* get proposal ratio */
    sum = 0.0;
    for (i=0; i<nStates; i++)
        sum += newPi[i]*alphaPi;
    x = LnGamma(sum);
    for (i=0; i<nStates; i++)
        x -= LnGamma(newPi[i]*alphaPi);
    for (i=0; i<nStates; i++)
        x += (newPi[i]*alphaPi-1.0)*log(oldPi[i]);
    sum = 0.0;
    for (i=0; i<nStates; i++)
        sum += oldPi[i]*alphaPi;
    y = LnGamma(sum);
    for (i=0; i<nStates; i++)
        y -= LnGamma(oldPi[i]*alphaPi);
    for (i=0; i<nStates; i++)
        y += (oldPi[i]*alphaPi-1.0)*log(newPi[i]);
    (*lnProposalRatio) = x - y;

    /* get prior ratio */
    y = x = 0.0;    /* the Gamma part of the prior is the same */
    for (i=0; i<nStates; i++)
        x += (symDirAlphai-1.0)*log(newPi[i]);
    for (i=0; i<nStates; i++)
        y += (symDirAlphai-1.0)*log(oldPi[i]);
    (*lnPriorRatio) = x - y;
        
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        
    /* Set update flags for cijks for all affected partitions. Only cijks for the changed character 
       actually need to be updated but we can't do that in the current version of the program. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


int Move_SwitchRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change switch rate of covarion model using sliding window */
    
    int         i, isSPriorExp, isValidS, whichRate;
    MrBFlt      oldS, newS, window, minS, maxS, sExp=0.0, ran, *value;
    ModelParams *mp;

    /* decide which switching rate to change */
    if (RandomNumber(seed) < 0.5)
        whichRate = 0;
    else
        whichRate = 1;
        
    /* get size of window, centered on current switching rates value */
    window = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* get minimum and maximum values for switching rate */
    if (param->paramId == SWITCH_UNI)
        {
        minS = mp->covswitchUni[0];
        maxS = mp->covswitchUni[1];
        isSPriorExp = NO;
        }
    else
        {
        minS = 0.01;
        maxS = KAPPA_MAX;
        sExp = mp->covswitchExp;
        isSPriorExp = YES;
        }

    /* get old value of switching rate */
    value = GetParamVals(param, chain, state[chain]);
    newS = oldS = value[whichRate];

    /* change value for switching rate */
    ran = RandomNumber(seed);
    if (maxS-minS < window)
        {
        window = maxS-minS;
        }
    newS = oldS + window * (ran - 0.5);
    
    /* check that new value is valid */
    isValidS = NO;
    do
        {
        if (newS < minS)
            newS = 2* minS - newS;
        else if (newS > maxS)
            newS = 2 * maxS - newS;
        else
            isValidS = YES;
        } while (isValidS == NO);

    /* get proposal ratio */
    *lnProposalRatio = 0.0;
    
    /* get prior ratio */
    if (isSPriorExp == NO)
        *lnPriorRatio = 0.0;
    else
        *lnPriorRatio = -sExp * (newS - oldS);
    
    /* copy new switching rate value back */
    value[whichRate] = newS;

    /* Set update flags for all partitions that share this switching rate. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for covarion, doublet, and codon models, we do want to update
       the cijk flag. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


int Move_SwitchRate_M (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change switch rate of covarion model using multiplier */

    int         i, isSPriorExp, isValidS, whichRate;
    MrBFlt      oldS, newS, minS, maxS, sExp=0.0, tuning, ran, factor, *value;
    ModelParams *mp;

    /* decide which switching rate to change */
    if (RandomNumber(seed) < 0.5)
        whichRate = 0;
    else
        whichRate = 1;
        
    /* get tuning parameter */
    tuning = mvp[0];

    /* get model params */
    mp = &modelParams[param->relParts[0]];
    
    /* get minimum and maximum values for switching rate */
    if (param->paramId == SWITCH_UNI)
        {
        minS = mp->covswitchUni[0];
        maxS = mp->covswitchUni[1];
        isSPriorExp = NO;
        }
    else
        {
        minS = 0.01;
        maxS = KAPPA_MAX;
        sExp = mp->covswitchExp;
        isSPriorExp = YES;
        }

    /* get old value of switching rate */
    value = GetParamVals(param, chain, state[chain]);
    newS = oldS = value[whichRate];

    /* change value for switching rate */
    ran = RandomNumber(seed);
    factor = exp(tuning * (ran - 0.5));
    newS = oldS * factor;
    
    /* check that new value is valid */
    isValidS = NO;
    do
        {
        if (newS < minS)
            newS = minS * minS / newS;
        else if (newS > maxS)
            newS = maxS * maxS / newS;
        else
            isValidS = YES;
        } while (isValidS == NO);

    /* get proposal ratio */
    *lnProposalRatio = log (newS / oldS);
    
    /* get prior ratio */
    if (isSPriorExp == NO)
        *lnPriorRatio = 0.0;
    else
        *lnPriorRatio = -sExp * (newS - oldS);
    
    /* copy new switching rate value back */
    value[whichRate] = newS;

    /* Set update flags for all partitions that share this switching rate. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for covarion, doublet, and codon models, we do want to update
       the cijk flag. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


int Move_TK02BranchRate (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* move one TK02 relaxed clock branch rate using multiplier */
    
    int         i;
    MrBFlt      newRate, oldRate, tuning, minR, maxR, nu, *tk02Rate, *brlens;
    TreeNode    *p = NULL;
    ModelInfo   *m;
    Tree        *t;
    TreeNode    *q;
    
    /* get model settings */
    m = &modelSettings[param->relParts[0]];
    
    /* get the tuning parameter */
    tuning = mvp[0];
    
    /* get the TK02 branch rate and effective branch length data */
    tk02Rate = GetParamVals (param, chain, state[chain]);
    brlens   = GetParamSubVals (param, chain, state[chain]);
    
    /* get tree */
    t = GetTree (param, chain, state[chain]);
    
    /* get minimum and maximum rate */
    minR = RATE_MIN;
    maxR = RATE_MAX;
    
    /* randomly pick a rate */
    do  {
        i = (int) (RandomNumber(seed) * (t->nNodes -2));
        p = t->allDownPass[i];
        }
    while (p->length < TIME_MIN);  // not ancestral fossil
    
    /* find new rate */
    oldRate = tk02Rate[p->index];
    newRate = oldRate * exp ((0.5 - RandomNumber(seed)) * tuning);
    
    /* reflect if necessary */
    while (newRate < minR || newRate > maxR)
        {
        if (newRate < minR)
            newRate = minR * minR / newRate;
        if (newRate > maxR)
            newRate = maxR * maxR / newRate;
        }
    
    tk02Rate[p->index] = newRate;
    
    /* calculate prior ratio */
    nu = *GetParamVals (m->tk02var, chain, state[chain]);
    (*lnPriorRatio) = LnRatioTK02LogNormal (tk02Rate[p->anc->index], nu*p->length, newRate, oldRate);
    if (p->left != NULL)
        {
        if (p->left->length > 0.0)
            {
            (*lnPriorRatio) -= LnProbTK02LogNormal (oldRate, nu*p->left->length,  tk02Rate[p->left->index ]);
            (*lnPriorRatio) += LnProbTK02LogNormal (newRate, nu*p->left->length,  tk02Rate[p->left->index ]);
            }
        if (p->right->length > 0.0)
            {
            (*lnPriorRatio) -= LnProbTK02LogNormal (oldRate, nu*p->right->length, tk02Rate[p->right->index]);
            (*lnPriorRatio) += LnProbTK02LogNormal (newRate, nu*p->right->length, tk02Rate[p->right->index]);
            }
        }
    
    /* calculate proposal ratio */
    (*lnProposalRatio) = log (newRate / oldRate);
    
    /* update branch evolution lengths */
    brlens[p->index] = p->length * (newRate + tk02Rate[p->anc->index]) / 2.0;
    if (p->left != NULL)
        {
        if (p->left->length > 0.0)
            {
            brlens[p->left->index] = p->left->length  * (tk02Rate[p->left->index] + newRate) / 2.0;
            }
        if (p->right->length > 0.0)
            {
            brlens[p->right->index] = p->right->length * (tk02Rate[p->right->index] + newRate) / 2.0;
            }
        }
    
    /* set update of ti probs */
    p->upDateTi = YES;
    if (p->left != NULL)
        {
        p->left ->upDateTi = YES;
        p->right->upDateTi = YES;
        }
    
    /* set update of cond likes down to root */
    /* update of crowntree set in UpdateCppEvolLengths */
    p->upDateCl = YES;
    q = p->anc;
    while (q->anc != NULL)
        {
        q->upDateCl = YES;
        q = q->anc;
        }
    
    return (NO_ERROR);
}


int Move_Tratio_Dir (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change tratio using Dirichlet proposal */
    
    int         i;
    MrBFlt      oldK, alphaPi, *alphaDir, oldProp[2], newProp[2], dirParm[2], sum, x, y;
    ModelParams *mp;

    /* get model params */
    mp = &modelParams[param->relParts[0]];

    /* get so called alphaPi parameter */
    alphaPi = mvp[0];

    /* get old value of kappa */
    oldK = *GetParamVals(param, chain, state[chain]);

    /* get Dirichlet parameters */
    alphaDir = mp->tRatioDir;

    /* calculate old ratesum proportions */
    oldProp[0] = oldK / (oldK + 1.0);
    oldProp[1] = 1.0 - oldProp[0];
    
    /* multiply old ratesum props with some large number to get new values close to the old ones */
    dirParm[0] = oldProp[0] * alphaPi;
    dirParm[1] = oldProp[1] * alphaPi;
    
    /* get new values */
    DirichletRandomVariable (dirParm, newProp, 2, seed);

    if (newProp[0] < DIR_MIN)
        {
        newProp[0] = DIR_MIN;
        newProp[1] = 1.0-DIR_MIN;
        }
    else if (newProp[1] < DIR_MIN)
        {
        newProp[1] = DIR_MIN;
        newProp[0] = 1.0-DIR_MIN;
        }

    /* calculate and copy new kappa value back */
    *GetParamVals(param, chain, state[chain]) = newProp[0] / newProp[1];

    /* get proposal ratio */
    sum = 0.0;
    for (i=0; i<2; i++)
        sum += newProp[i]*alphaPi;
    x = LnGamma(sum);
    for (i=0; i<2; i++)
        x -= LnGamma(newProp[i]*alphaPi);
    for (i=0; i<2; i++)
        x += (newProp[i]*alphaPi-1.0)*log(oldProp[i]);
    sum = 0.0;
    for (i=0; i<2; i++)
        sum += oldProp[i]*alphaPi;
    y = LnGamma(sum);
    for (i=0; i<2; i++)
        y -= LnGamma(oldProp[i]*alphaPi);
    for (i=0; i<2; i++)
        y += (oldProp[i]*alphaPi-1.0)*log(newProp[i]);
    (*lnProposalRatio) = x - y;

    /* get prior ratio */
    x = y = 0.0;
    for (i=0; i<2; i++)
        x += (alphaDir[i]-1.0)*log(newProp[i]);
    for (i=0; i<2; i++)
        y += (alphaDir[i]-1.0)*log(oldProp[i]);
    (*lnPriorRatio) = x - y;
        
    /* Set update flags for all partitions that share this kappa. Note that the conditional
       likelihood update flags have been set before we even call this function. */
    for (i=0; i<param->nRelParts; i++)
        TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);

    /* Set update flags for cijks for all affected partitions. If this is a simple 4 X 4 model,
       we don't take any hit, because we will never go into a general transition probability
       calculator. However, for covarion, doublet, and codon models, we do want to update
       the cijk flag. */
    for (i=0; i<param->nRelParts; i++)
        modelSettings[param->relParts[i]].upDateCijk = YES;

    return (NO_ERROR);
}


/* Code added by Jeremy Brown and modified by Maxim Teslenko */
int Move_TreeLen (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    /* change all branch lengths */

    MrBFlt      begin_tl, treescaler, tuning, maxV, minV, brlensPrExp=0.0;
    TreeNode    *p;
    ModelParams *mp;
    Tree        *t;
    int i,branch_counter,  isVPriorExp;

    tuning = mvp[0]; /* Larget & Simon's tuning parameter lambda */

    mp = &modelParams[param->relParts[0]];

    /* max and min brlen */
    if (param->paramId == BRLENS_UNI)
        {
        minV = mp->brlensUni[0];
        maxV = mp->brlensUni[1];
        isVPriorExp = NO;
        }
    else if (param->paramId == BRLENS_GamDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 2;
        }
    else if (param->paramId == BRLENS_iGmDir)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 3;
        }
    else if (param->paramId == BRLENS_twoExp)
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        isVPriorExp = 4;
        }
    else
        {
        minV = BRLENS_MIN;
        maxV = BRLENS_MAX;
        brlensPrExp = mp->brlensExp;
        isVPriorExp = YES;
        }

    /* get tree */
    t = GetTree (param, chain, state[chain]);

    assert (t->isRooted == NO);

    /* Dirichlet or twoExp prior */
    if (isVPriorExp > 1)
        (*lnPriorRatio) = -LogDirPrior(t, mp, isVPriorExp);

    treescaler = exp(tuning * (RandomNumber(seed) - 0.5));
    
    begin_tl = 0.0;
    branch_counter=0;

    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->anc != NULL)
            {
            if (p->length*treescaler < minV || p->length*treescaler > maxV)
                {
                abortMove = YES;
                return NO_ERROR;
                }
            begin_tl += p->length;
            branch_counter++;               
            }
        }
    assert (branch_counter==t->nNodes-1);
    
    /* iterate scaling over all branches */
    for (i=0; i < t->nNodes; i++)
        {
        p = t->allDownPass[i];
        if (p->anc != NULL)
            {
            /* set new length */
            p->length *= treescaler;

            /* set flags for update of transition probabilities at p */
            p->upDateTi = YES;
            p->anc->upDateCl = YES; 
            }
        }

    /* calculate proposal ratio */
    (*lnProposalRatio) = branch_counter * log(treescaler);

    /* update prior if exponential prior on branch lengths */
    if (param->paramId == BRLENS_EXP)
        (*lnPriorRatio) = brlensPrExp * (begin_tl* (1 - treescaler));
    /* Dirichlet or twoExp prior */
    else if (isVPriorExp > 1)
        (*lnPriorRatio) += LogDirPrior(t, mp, isVPriorExp);

    return (NO_ERROR);
}


/*-----------------------------------------------------------------------------------
|
|   Move_TreeStretch: Shrink or grow a clock tree
|
-------------------------------------------------------------------------------------*/
int Move_TreeStretch (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int         i, j, *nEvents, numChangedNodes;
    MrBFlt      minV, maxV, tuning, factor, lambda=0.0, x,
                *brlens=NULL, nu=0.0, igrvar=0.0, *tk02Rate=NULL, *igrRate=NULL;
    TreeNode    *p, *q;
    ModelParams *mp;
    ModelInfo   *m;
    Tree        *t, *oldT;
    Param       *subParm;
    Calibration *calibrationPtr;

    tuning = mvp[0]; /* Larget & Simon tuning parameter lambda */
 
    (*lnProposalRatio) = (*lnPriorRatio) = 0.0;

    m = &modelSettings[param->relParts[0]];
    mp = &modelParams[param->relParts[0]];

    /* get trees */
    t = GetTree (param, chain, state[chain]);
    oldT = GetTree (param, chain, 1^state[chain]);

    /* min and max branch lengths in relative time and substitution units */
    minV = BRLENS_MIN;
    maxV = BRLENS_MAX;

    /* determine multiplication factor */
    factor = exp(tuning * (RandomNumber(seed) - 0.5));

    /* multiply all changeable ages and node depths by this factor */
    numChangedNodes = 0;
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
            
        /* skip extant tip and fixed calibration */
        if (p->left == NULL && p->isDated == NO)
            continue;
        if (p->isDated == YES)
            calibrationPtr = p->calibration;
        else if (p->anc->anc == NULL && (!strcmp(mp->clockPr,"Uniform") || !strcmp(mp->clockPr,"Fossilization")))
            calibrationPtr = &mp->treeAgePr;
        else
            calibrationPtr = NULL;
        if (calibrationPtr != NULL && calibrationPtr->prior == fixed)
            continue;
        
        /* now stretch the node */
        if (calibrationPtr != NULL)
            {
            p->age *= factor;
            if (p->age < calibrationPtr->min || p->age > calibrationPtr->max)
                {
                abortMove = YES;
                return (NO_ERROR);
                }
            }
        p->nodeDepth *= factor;
        numChangedNodes++;
        
        /* deal with ancestral fossils */
        if (p->left != NULL)
            {
            if (p->left->length < TIME_MIN)
                {
                p->left->length = 0.0;
                p->nodeDepth = p->left->nodeDepth;
                if (calibrationPtr != NULL)
                    {
                    assert (p->left->calibration != NULL);
                    p->age = p->left->age;
                    if (p->age < calibrationPtr->min || p->age > calibrationPtr->max)
                        {
                        abortMove = YES;
                        return (NO_ERROR);
                        }
                    }
                numChangedNodes--;
                }
            if (p->right->length < TIME_MIN)
                {
                p->right->length = 0.0;
                p->nodeDepth = p->right->nodeDepth;
                if (calibrationPtr != NULL)
                    {
                    assert (p->right->calibration != NULL);
                    p->age = p->right->age;
                    if (p->age < calibrationPtr->min || p->age > calibrationPtr->max)
                        {
                        abortMove = YES;
                        return (NO_ERROR);
                        }
                    }
                numChangedNodes--;
                }
            assert (!(p->left->length == 0.0 && p->right->length == 0.0));
            }
        }
    
    /* update brls */
    for (i=0; i<t->nNodes-1; i++)
        {
        p = t->allDownPass[i];
        if (p->left != NULL)
            {
            if (p->left->length > 0.0)
                p->left->length = p->nodeDepth - p->left->nodeDepth;
            if (p->right->length > 0.0)
                p->right->length = p->nodeDepth - p->right->nodeDepth;
            }
        }

    /* check that all branch lengths are proper, which need not be the case */
    for (i = 0; i < t->nNodes -2; i++)
        {
        p = t->allDownPass[i];
        q = oldT->allDownPass[i];
        if (p->length < 0.0 || p->length > maxV || (q->length > minV && p->length < minV) || (q->length < TIME_MIN && p->length > TIME_MIN))
            {  /* consider ancestral fossil (brl=0) in fossilized bd tree */
            abortMove = YES;
            return NO_ERROR;
            }
        }

    /* calculate proposal ratio */
    (*lnProposalRatio) = numChangedNodes * log(factor);

    /* calculate prior ratio */
    if (LogClockTreePriorRatio(param, chain, &x) == ERROR)
        return ERROR;
    (*lnPriorRatio) += x;

    /* adjust proposal and prior ratio for relaxed clock models */
    for (i=0; i<param->nSubParams; i++)
        {
        subParm = param->subParams[i];
        if (subParm->paramType == P_CPPEVENTS)
            {
            nEvents = subParm->nEvents[2*chain+state[chain]];
            lambda = *GetParamVals (modelSettings[subParm->relParts[0]].cppRate, chain, state[chain]);
            /* proposal ratio */
            for (j=0; j<t->nNodes-2; j++)
                {
                p = t->allDownPass[j];
                q = oldT->allDownPass[j];
                (*lnProposalRatio) += nEvents[p->index ] * log (p->length  / q->length);
                }
            /* prior ratio */
            (*lnPriorRatio) += lambda * (TreeLen(oldT) - TreeLen(t));
            /* update effective evolutionary lengths */
            if (UpdateCppEvolLengths (subParm, t->root->left, chain) == ERROR)
                {
                abortMove = YES;
                return (NO_ERROR);
                }
            }
        else if ( subParm->paramType == P_TK02BRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_TK02))
            {
            if (subParm->paramType == P_TK02BRANCHRATES)
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].tk02var, chain, state[chain]);
            else
                nu = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            tk02Rate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);

            /* prior ratio and update of brlens */
            for (j=0; j<t->nNodes-2; j++)
                {
                p = t->allDownPass[j];
                q = oldT->allDownPass[j];
                if (p->length > 0.0)  // not ancestral fossil
                    {
                    (*lnPriorRatio) -= LnProbTK02LogNormal (tk02Rate[q->anc->index], nu*q->length, tk02Rate[q->index]);
                    (*lnPriorRatio) += LnProbTK02LogNormal (tk02Rate[p->anc->index], nu*p->length, tk02Rate[p->index]);
                    brlens[p->index] = p->length * (tk02Rate[p->anc->index]+tk02Rate[p->index])/2.0;
                    if (brlens[p->index] < RELBRLENS_MIN || brlens[p->index] > RELBRLENS_MAX)
                        {
                        abortMove = YES;
                        return (NO_ERROR);
                        }
                    }
                }
            }
        else if ( subParm->paramType == P_IGRBRANCHRATES ||
                 (subParm->paramType == P_MIXEDBRCHRATES && *GetParamIntVals(subParm, chain, state[chain]) == RCL_IGR))
            {
            if (subParm->paramType == P_IGRBRANCHRATES)
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].igrvar, chain, state[chain]);
            else
                igrvar = *GetParamVals (modelSettings[subParm->relParts[0]].mixedvar, chain, state[chain]);
            igrRate = GetParamVals (subParm, chain, state[chain]);
            brlens = GetParamSubVals (subParm, chain, state[chain]);
            
            /* prior ratio and update of igr branch lengths and rates (stretched in the same way as tree) */
            for (j=0; j<t->nNodes-2; j++)
                {
                p = t->allDownPass[j];
                q = oldT->allDownPass[j];
                if (p->length > 0.0)  // not ancestral fossil
                    {
                    (*lnPriorRatio) -= LnProbGamma (q->length/igrvar, q->length/igrvar, igrRate[q->index]);
                    (*lnPriorRatio) += LnProbGamma (p->length/igrvar, p->length/igrvar, igrRate[p->index]);
                    brlens[p->index] = p->length * igrRate[p->index];
                    if (brlens[p->index] < RELBRLENS_MIN || brlens[p->index] > RELBRLENS_MAX)
                        {
                        abortMove = YES;
                        return (NO_ERROR);
                        }
                    }
                }
            }
        }

    TouchAllTreeNodes(m, chain);

#if defined (DEBUG_TREESTRETCH)
    printf ("After treestretch:\n");
    printf ("Old tree height: %f -- New tree height: %f -- lnPriorRatio = %f -- lnProposalRatio = %f\n",
        oldT->root->left->nodeDepth, t->root->left->nodeDepth, (*lnPriorRatio), (*lnProposalRatio));
#endif

    return (NO_ERROR);
}


/* Generalized normal move for real random variables */
int Move_RealNormal (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             i;
    MrBFlt          oldX, newX, tuning, minX, maxX, u, z;

    /* get tuning parameter */
    tuning = mvp[0];

    /* get minimum and maximum values for X */
    minX = param->min;
    maxX = param->max;

    /* get old value of X */
    newX = oldX = *GetParamVals(param, chain, state[chain]);

    /* change value */
    u = RandomNumber(seed);
    z = PointNormal(u);
    newX = oldX + z * tuning;
    
    /* check that new value is valid */
    if (newX < minX || newX > maxX)
        {
        abortMove = YES;
        return (NO_ERROR);
        }

    /* get proposal ratio */
    (*lnProposalRatio) = 0.0;
    
    /* get prior ratio */
    (*lnPriorRatio) = param->LnPriorRatio(newX, oldX, param->priorParams);

    /* copy new value back */
    *GetParamVals(param, chain, state[chain]) = newX;

    /* Set update flags for tree nodes if relevant */
    if (param->affectsLikelihood == YES)
        {
        for (i=0; i<param->nRelParts; i++)
            TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        }

    return (NO_ERROR);
}


/* Generalized slider move for real random variables */
int Move_RealSlider (Param *param, int chain, RandLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)
{
    int             i, isValid;
    MrBFlt          oldX, newX, window, minX, maxX, u;

    /* get size of window, centered on current value */
    window = mvp[0];

    /* get minimum and maximum values for X */
    minX = param->min;
    maxX = param->max;

    /* get old value of X */
    newX = oldX = *GetParamVals(param, chain, state[chain]);

    /* change value */
    u = RandomNumber(seed);
    newX = oldX + window * (u - 0.5);
    
    /* check that new value is valid */
    isValid = NO;
    do
        {
        if (newX < minX)
            newX = 2* minX - newX;
        else if (newX > maxX)
            newX = 2 * maxX - newX;
        else
            isValid = YES;
        } while (isValid == NO);

    /* get proposal ratio */
    (*lnProposalRatio) = 0.0;
    
    /* get prior ratio */
    (*lnPriorRatio) = param->LnPriorRatio(newX, oldX, param->priorParams);

    /* copy new value back */
    *GetParamVals(param, chain, state[chain]) = newX;

    /* Set update flags for tree nodes if relevant */
    if (param->affectsLikelihood == YES)
        {
        for (i=0; i<param->nRelParts; i++)
            TouchAllTreeNodes(&modelSettings[param->relParts[i]],chain);
        }

    return (NO_ERROR);
}


void TouchAllTreeNodes (ModelInfo *m, int chain)
{
    int         i;
    Tree        *t;
    TreeNode    *p;
    
    t = GetTree(m->brlens, chain, state[chain]);
    for (i=0; i<t->nNodes; i++)
        {
        p = t->allDownPass[i];
        p->upDateCl = YES;
        p->upDateTi = YES;
        }
    m->upDateAll = YES;
}

