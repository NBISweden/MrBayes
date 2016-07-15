# by Paul O. Lewis
# paul.lewis@uconn.edu

import sys
from math import log,exp,lgamma,pow
import scipy.stats, scipy.optimize

def Pr(from_state, to_state, pi0, exp_minus_mu_t):
    # compute transition probability from state from_state to state to_state
    # with exp_minus_mu_t = exp(-mu*t), where mu*t is computed from the branch
    # length v and relative frequencies of state 0 (pi0) and state 1 (pi1) by
    # solving the following equation for mu*t:
    #     v = 2*pi0*pi1*mu*t
    pi1 = 1.0 - pi0
    if from_state == to_state:
        if from_state == 0:
            return pi0 + pi1*exp_minus_mu_t
        else:
            return pi1 + pi0*exp_minus_mu_t
    else:
        if from_state == 0:
            return pi1 - pi1*exp_minus_mu_t
        else:
            return pi0 - pi0*exp_minus_mu_t

def logLikelihood(pi0, v, state1, state2, state3, state4):
    # computes log-likelihood assuming tree (1,2,(3,4)) with all 5 branch lengths equal to v
    # and assuming relative equilibrium frequency of state 0 is pi0
    pi1 = 1.0 - pi0
    mu_t = v/(2.0*pi0*pi1)
    eterm = exp(-mu_t)
    condlike0 = Pr(0, state3, pi0, eterm) * Pr(0, state4, pi0, eterm)
    condlike1 = Pr(1, state3, pi0, eterm) * Pr(1, state4, pi0, eterm)
    condlikesum0 = Pr(0, 0, pi0, eterm) * condlike0 + Pr(0, 1, pi0, eterm) * condlike1
    condlikesum1 = Pr(1, 0, pi0, eterm) * condlike0 + Pr(1, 1, pi0, eterm) * condlike1
    root0 = Pr(0, state1, pi0, eterm) * Pr(0, state2, pi0, eterm) * condlikesum0
    root1 = Pr(1, state1, pi0, eterm) * Pr(1, state2, pi0, eterm) * condlikesum1
    loglike = log(pi0*root0 + pi1*root1)
    return loglike

if __name__ == '__main__':
    ncateg = 5      # number of beta categories used for discrete beta frequency distribution

    # if prset Symdirihyperpr = Fixed(x), then alpha = beta = x
    alpha  = 1.0
    beta   = 1.0

    # starting tree = (1:0.1,2:0.1,(3:0.1,4:0.1):0.1)
    v      = 0.1    # all branch lengths equal v

    # states for one character are hard-coded in this example
    state1 = 0      # state for tip 1
    state2 = 0      # state for tip 2
    state3 = 1      # state for tip 3
    state4 = 1      # state for tip 4

    print 'pi0 ~ Beta(%.5f, %.5f) divided into %d equal-volume categories' % (alpha, beta, ncateg)

    likelihood = 0.0
    for i in range(ncateg):
        # Find point on x-axis (pi0) that represents the median of the (i+1)th category
        category_median = (0.5+i)/ncateg
        pi0 = scipy.stats.beta.ppf(category_median, alpha, beta)

        loglike = logLikelihood(pi0, v, state1, state2, state3, state4)
        likelihood += exp(loglike)

        print '  categ %d: pi0 = %.5f, pi1 = %.5f, logLike = %.5f' % (i+1, pi0, 1.0-pi0, loglike)

    # the likelihood for each category (likelihood) is multiplied by the
    # probability of being in that category (1/ncateg)
    loglikelihood = log(likelihood) - log(ncateg)

    print 'log-likelihood =',loglikelihood

# [corresponding mb block]
#
# #nexus
#
# begin data;
#   dimensions ntax=4 nchar=1;
#   format datatype=standard;
#   matrix
#     A 0
#     B 0
#     C 1
#     D 1
#   ;
# end;
# 
# begin trees;
#   translate
#     1 A,
#     2 B,
#     3 C,
#     4 D;
#   tree only = (1:0.1,2:0.1,(3:0.1,4:0.1):0.1);
# end;
# 
# begin mrbayes;
#   lset nbetacat=5 coding=all;
#   prset brlenspr = unconstrained:exponential(1.0);
#   prset Symdirihyperpr = Exponential(0.5);
#   mcmcp ngen=1 nchain=1 nrun=1;
#   startvals Tau=only V=only;
#   mcmc;
# end;

