/*
 *  MrBayes 3
 *
 *  copyright 2002--
 *
 *  John P. Huelsenbeck
 *  Section of Ecology, Behavior and Evolution
 *  Division of Biological Sciences
 *  University of California, San Diego
 *  La Jolla, CA 92093-0116
 *
 *  johnh@biomail.ucsd.edu
 *
 *	Fredrik Ronquist
 *  Paul van der Mark
 *  School of Computational Science
 *  Florida State University
 *  Tallahassee, FL 32306-4120
 *
 *  fredrik.ronquist@nrm.se
 *  paulvdm@scs.fsu.edu
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

#define		TREEBUFINCREASE   200

/* tree utility functions */
int       AddToTreeList (TreeList *treeList, Tree *tree);
Tree     *AllocateTree (int numTaxa);
Tree     *AllocateFixedTree (int numTaxa, int isRooted);
int       AllocateTreeFlags (Tree *t);
int       AllocateTreePartitions (Tree *t);
PolyTree *AllocatePolyTree (int numTaxa);
int       AllocatePolyTreePartitions (PolyTree *pt);
int       AllocatePolyTreeRelClockParams (PolyTree *pt, int nBSets, int nESets);
int		  AreTopologiesSame(Tree *t1, Tree *t2);
int		  AreTreesSame (Tree *t1, Tree *t2);
int       BuildConstraintTree (Tree *t, PolyTree *pt, char **localTaxonNames);
int       BuildRandomRTopology (Tree *t, safeLong *seed);
int       BuildRandomUTopology (Tree *t, safeLong *seed);
int       CheckConstraints (Tree *t);
int       CheckSetConstraints (Tree *t);
void      ColorClusters (TreeNode *p, int *index);
void      CopySubtreeToTree (Tree *subtree, Tree *t);
int       CopyToPolyTreeFromPolyTree (PolyTree *to, PolyTree *from);
int       CopyToTreeFromPolyTree (Tree *to, PolyTree *from);
void      CopyPolyNodes (PolyNode *p, PolyNode *q);
int       CopyToTreeFromTree (Tree *to, Tree *from);
void      CopyTreeNodes (TreeNode *p, TreeNode *q);
void      CopyTreeToSubtree (Tree *t, Tree *subtree);
int       Deroot(PolyTree *pt);
void      EraseTreeList (TreeList *treeList);
void      FreePolyTree (PolyTree *pt);
void      FreePolyTreePartitions (PolyTree *pt);
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
int       GetRandomEmbeddedSubtree (Tree *t, int nTerminals, safeLong *seed, int *nEmbeddedTrees);
int       GetFromTreeList (TreeList *treeList, Tree *tree);
int       InitBrlens (Tree *t, MrBFlt v);
int       InitCalibratedBrlens (Tree *t, MrBFlt minLength, safeLong *seed);
int       InitClockBrlens (Tree *t);
int	      IsCalibratedClockSatisfied (Tree *t, MrBFlt tol);
int       IsClockSatisfied (Tree *t, MrBFlt tol);
int       LabelTree (Tree *t, char **taxonNames);
void      Mark (TreeNode *p);
void      MarkDatedSubtree (TreeNode *p);
int	      MoveCalculationRoot (Tree *t, int outgroup);
int	      MovePolyCalculationRoot (PolyTree *t, int outgroup);
int       NumConstrainedTips (TreeNode *p);
int       NumDatedTips (TreeNode *p);
void      OrderTips (PolyTree *t);
void      PrintNewick (char **s, int *len, Tree *t);
void      PrintNodes (Tree *t);
void      PrintPolyNodes (PolyTree *pt);
int       PrunePolyTree (PolyTree *pt);
int       RandPerturb (Tree *t, int nPert, safeLong *seed);
int       RandResolve (PolyTree *t, safeLong *seed, int destinationIsRooted);
int       ResetBrlensFromTree (Tree *tree, Tree *vTree);
void      ResetPolyTree (PolyTree *t);
void      ResetPolyTreePartitions (PolyTree *pt);
void      ResetPolyTreeRelClockParams (PolyTree *pt);
int       ResetRootHeight (Tree *t, MrBFlt rootHeight);
int       ResetTopology (Tree *t, char *s);
int       ResetTopologyFromTree (Tree *tree, Tree *top);
int       ResetTopologyFromPolyTree (Tree *tree, PolyTree *top);
void      ResetTreePartitions (Tree *t);
int       RetrieveRTopology (Tree *t, int *order);
int       RetrieveRTree (Tree *t, int *order, MrBFlt *brlens);
int       RetrieveUTopology (Tree *t, int *order);
int       RetrieveUTree (Tree *t, int *order, MrBFlt *brlens);
int       ShowTree (Tree *t);
int       StoreRPolyTopology (PolyTree *t, int *order);
int       StoreRPolyTree (PolyTree *t, int *order, MrBFlt *brlens);
int       StoreRTopology (Tree *t, int *order);
int       StoreRTree (Tree *t, int *order, MrBFlt *brlens);
int       StoreUPolyTopology (PolyTree *t, int *order);
int       StoreUPolyTree (PolyTree *t, int *order, MrBFlt *brlens);
int       StoreUTopology (Tree *t, int *order);
int       StoreUTree (Tree *t, int *order, MrBFlt *brlens);
MrBFlt    TreeLen (Tree *t);
void      Unmark (TreeNode *p);
void      WriteEventTree (TreeNode *p, int chain, Param *param);
void      WriteEventTreeToPrintString (TreeNode *p, int chain, Param *param, int printAll);
void      WriteEvolTree (TreeNode *p, int chain, Param *param);
void      WriteTopologyToFile (FILE *fp, TreeNode *p, int isRooted);
void      WriteTreeToPrintString (TreeNode *p, int showBrlens, int isRooted);
