/*
 *  Best 2.2
 *
 *  This file contains the functions 
 *  for calculating the probability of 
 *  gene trees given the species tree 
 *  and the prior probability of the 
 *  species tree
 *
 *  Liang Liu
 *  Department of Statistics
 *  The Ohio State University
 *  Columbus, Ohio
 *  
 *  liuliang@stat.ohio-state.edu
 */

#include	"best.h"
#include	"command.h"
#include "globals.h"
#include "mb.h"
#include "model.h"
#include "tree.h"

double      NodeDistance(SPTree *tree, int inode, int jnode);
int	        ReadaTree (FILE *fTree,SPTree *tree);
void        ReadGeneTree(FILE *fTree);
void        SPWriteTreeToFile (SPTree *tree, int inode, int showBrlens, int showTheta, int showMu, int isRooted);
int         Constraint(SPTree *genetree, int numgenetree, SPTree *speciestree, Distance *constraint);
double      CalNodeAge(int node, SPTree *tree);
double      PopNodeDistance(int inode, int jnode, SPTree *tree);
int         FindaPosition(int nodenumber,int root,double currentdistance,double *position,SPTree *speciestree);  
int         StartSptree(SPTree *speciestree, int numchange);
void        quick_struct(Distance *item,int count);
void        qs_struct(Distance *item,int left,int right);
int         GetCoaltime(SPTree *genetree,SPTree *speciestree, CoalTime *coal);
int         GetNcoal(int inode, CoalTime *coal, SPTree *genetree, SPTree *speciestree, CoalTime *coal_t,int *index);
void        SPMrBayesPrint (char *format, ...);
void        PrintInf(FILE *fp);
int         SPLogLike(SPTree *genetree, SPTree *speciestree, double *lnl);
int         LnLikehood1Tree(SPTree *genetree, SPTree *speciestree, double *lnp);
void        SPPrintToScreen (int curGen, time_t endingT, time_t startingT);
int         SPLogPrior(SPTree *speciestree, double *lnprior);
int         Move_SPSpeciation (SParam *param, int chain, long *seed, double *lnLikeRatio, double *lnPriorRatio, double *lnProposalRatio, double *mvp);
int         Move_SPExtinction (SParam *param, int chain, long *seed, double *lnLikeRatio, double *lnPriorRatio, double *lnProposalRatio, double *mvp);
double 	    SPLnP1 (double t, double l, double m, double r);
double 	    SPLnVt (double t, double l, double m, double r);
void        CopyParam(int chn);
int         SPPickProposal (void);
void        SPPreparePrintFiles (void);
void        DelAllmissingspecies(int nummissing,int *missnode, SPTree *speciestree);
void        DeleteaSpecies(int inode,SPTree *speciestree);
int         SPAddToPrintString (char *tempStr);
double 	    Toclocktree(SPTree *t, int node);
double 	    TreeL(Tree *t);
void        ToGenetree(Tree *file[],int *updatedtreeid, int nfile,double *GeneMu);
MrBFlt 	    Prob_sptree(Distance *tau,int ntime);
int         LnLikehood1Tree_invgamma(SPTree *genetree, SPTree *speciestree, double *a, double *b);
int         SPLogLike_invgamma(SPTree *genetree, SPTree *speciestree, double *lnl);
int         poisson(double x);
double 	    ChangeBrlen(SPTree *speciestree, int spnode, Tree *genetree, TreeNode *p);
int         FindSpnodeDownGenenode(SPTree *speciestree, int spnode, TreeNode *p);
int         FindspNodeBelow(SPTree *speciestree, int spnode, double dis_gene);
int         SPTreeConstraint(Distance *minimumdistance, Distance *distance, long int nconstraints, int nspecies);
int         MaximumTree(Distance *dist, SPTree *speciestree, int nconstraints);
int         ChangeConstraint(Distance *dist, int nconstraints);
void        ToSingleGenetree(Tree *file,int i,MrBFlt GeneMu);
void        FindDescendantTaxa(SPTree *tree, int inode, int *taxa, int *index);
long int    ClockTreeNodeDist(SPTree *clocktree, int ngene, Distance *dist);
int 	    CheckConstraint(SPTree *genetrees, int ngene, Distance *constraint, int nconstraints);
long int    OneClockTreeNodeDist(SPTree *clocktree, Distance *dist);
int 		populationMutation (Tree *genetree, SPTree *speciestree, MrBFlt genemu);
int			SPSaveSprintf(char **target, int *targetLen, char *fmt, ...);
int         GetConstraints(SPTree *s, Distance *constr);
void        GetMinDists(SPTree *clocktree, double **md);

McmcPara 	mcmc;
SPTree 	sptree;
SPTree	*gtree;
ModelParam	modelParam;
int		nGene;
char		spacer[10]="  ";
FILE		*fptree;     		/* Output tree file */
FILE		*fpparm;     		/* Output parameter file */          
int		*spnode;       /*Vector of taxaIDs for each tip */
long int 	seed;
int		curGeneration=0;  /*should be defined in mcmc.c and used to print species tree*/
double	speciationR, extinctionR, sampleF; 

/*global variables*/
int		jointGenePr=0 ; 		/*best or mrbayes*/
int		spMupr = 0; 		/*non-clock species tree or clock species tree*/
int		spMupralpha = 5; 		/*for non-clock species tree model, gamma(spMupralpha,r/spMupralpha) for mutation rate ratio r*/
int		speciestreePr = 0; 	/*birth-death prior or uniform prior*/
double 	poissonMean = 1.0; 	/*the number of nodes changed for proposing a new species tree*/



/**************************** Functions that interact directly with MrBayes ****************************/

/*------------------------------------------------------------------
|
|	FillSpeciesTreeParams: Fill in species trees (start value)
|
------------------------------------------------------------------*/
int FillSpeciesTreeParams (SafeLong *seed, int fromChain, int toChain)

{

    // TODO: BEST Build an appropriate starting tree

    int			k, chn;
	Param		*p;
	Tree		*tree;

	/* Build species trees for state 0 */
	for (chn=fromChain; chn<toChain; chn++)
		{
		for (k=0; k<numParams; k++)
			{
			p = &params[k];
			if (p->paramType == P_SPECIESTREE && p->fill == YES)
				{
                tree = GetTree(p, chn, 0);
                BuildRandomRTopology(tree, seed);
                InitClockBrlens(tree);
				if (LabelTree (tree, localTaxonNames) == ERROR)
					return (ERROR);
				if (chn == toChain-1)	/* last chain to fill */
					p->fill = NO;
                }
            }
        }

	return (NO_ERROR);
}





/*------------------------------------------------------------------
|
|	Move_GeneTree: Propose a new gene tree
|
------------------------------------------------------------------*/
int Move_GeneTree (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)

{
    Tree			*t, *spt;
    ModelInfo       *m;
    ModelParams     *mp;

	/* get model params */
	mp = &modelParams[param->relParts[0]];
	
	/* get model settings */
    m = &modelSettings[param->relParts[0]];

    /* get gene tree */
    t = GetTree (param, chain, state[chain]);

    /* get species tree */
    spt = GetTree (m->speciestree, chain, state[chain]);

    // TODO: BEST Modify gene tree.
    //printf ("Proposing new gene tree, index = %d\n", param->treeIndex);

    // TODO: BEST Calculate proposal ratio
    (*lnProposalRatio) = 0.0;

    // TODO: BEST Calculate prior ratio taking species tree into account
    (*lnPriorRatio) = 0.0;
    
    return (NO_ERROR);
}





/*------------------------------------------------------------------
|
|	Move_SpeciesTree: Propose a new species tree
|
------------------------------------------------------------------*/
int Move_SpeciesTree (Param *param, int chain, SafeLong *seed, MrBFlt *lnPriorRatio, MrBFlt *lnProposalRatio, MrBFlt *mvp)

{

    /* Move species tree */

    int             i;
    Tree			*t, *genetree;
    ModelInfo       *m;
    ModelParams     *mp;

	/* get model params */
	mp = &modelParams[param->relParts[0]];
	
	/* get model settings */
    m = &modelSettings[param->relParts[0]];

    /* get species tree */
    t = GetTree (param, chain, state[chain]);

    /* cycle over gene trees */
    for (i=0; i<param->nSubParams; i++)
        genetree = GetTree(param->subParams[i], chain, state[chain]);

    // TODO: BEST code needed here::
    // printf ("Modifying species tree...\n");
    
    // Modify the species tree, given info on the gene trees

    // Calculate proposal ratio
    (*lnProposalRatio) = 0.0;

#if defined (BEST_MPI_ENABLED)
    // Broadcast the proposed species tree to all processors if MPI version
#endif

    // Calculate the ln prior probability ratio of the new to old species trees from hyperpriors
    (*lnPriorRatio) = 0.0;

#if defined (BEST_MPI_ENABLED)
    // Let each processor calculate the ln probability ratio of its current gene tree(s)
    //    given the new and old species tree in the MPI version

    // Assemble the ln probability ratios across the processors and to lnPriorRatio
#else
    // Calculate the ln probability ratio of the current gene trees
    //    given the new and old species trees

#endif

    // Add ln probability ratio to (*lnPriorRatio)
    (*lnPriorRatio) += 0.0;
    
    return (NO_ERROR);

}





/*************************** best functions ***********************************/

//recursive function that loads the IDs of the (index) tips descending from inode into taxa
void FindDescendantTaxa(SPTree *tree, int inode, int *taxa, int *index) {
	if(inode < tree->nTaxa) {
		taxa[(*index)++] = inode;
	} else {
		FindDescendantTaxa(tree, tree->nodes[inode].sons[0], taxa, index);
		FindDescendantTaxa(tree, tree->nodes[inode].sons[1], taxa, index);
	}
}

int poisson(double x)
{
	int    poi_value;             
  	double t_sum;                 

  	poi_value = 0;
  	t_sum = 0.0;

  	while(1)
  	{
    		t_sum = t_sum - x * log(rndu());
    		if (t_sum >= 1.0) break;
    			poi_value++;
  	}
 	return(poi_value);
}

void ToSingleGenetree(Tree *file,int i,MrBFlt GeneMu)
{
  	TreeNode *p;
  	int j,nspecies;
  	MrBFlt tlg,tl;
  
  	nspecies = gtree[i].nTaxa;
      
   if(file->isRooted == 1)
   {
	gtree[i].root = file->root->index;
	gtree[i].nodes[gtree[i].root].sons[0] = file->root->left->index;
	gtree[i].nodes[gtree[i].root].sons[1] = file->root->right->index;
   }
   else
   {
   	gtree[i].root = nspecies*2-2;
   	gtree[i].nodes[gtree[i].root].sons[0] = file->root->left->index;
   	gtree[i].nodes[gtree[i].root].sons[1] = file->root->index;
   }
   gtree[i].nodes[gtree[i].root].nson = 2;
   gtree[i].nodes[gtree[i].root].father = -1;
   gtree[i].nodes[gtree[i].root].brlens = 0.0;
   for(j=0;j<file->nNodes;j++)
     { 
       p=file->allDownPass[j];
       if(p != file->root)
          {
              
              if(p->anc != NULL)
                     {
                      
                      if(p != file->root->left) 
                           {
                              gtree[i].nodes[p->index].brlens = p->length;
                              gtree[i].nodes[p->index].father = p->anc->index;
                           }
                      else {
                            gtree[i].nodes[p->index].brlens = p->length/2;
                            gtree[i].nodes[p->index].father = gtree[i].root;
                           }
                     }
              else gtree[i].nodes[p->index].father = -1;
              if(p->left != NULL)
                    { 
                     gtree[i].nodes[p->index].nson = 2;
                     gtree[i].nodes[p->index].sons[0] = p->left->index;
                     gtree[i].nodes[p->index].sons[1] = p->right->index;
                    }
              else
                    { 
                       gtree[i].nodes[p->index].nson = 0;
                       gtree[i].nodes[p->index].sons[0] = -2;
                       gtree[i].nodes[p->index].sons[1] = -2;
                    }
          }
       else
          {
                       gtree[i].nodes[p->index].nson = 0;
                       gtree[i].nodes[p->index].sons[0] = -2;
                       gtree[i].nodes[p->index].sons[1] = -2;
                       gtree[i].nodes[p->index].brlens = p->left->length/2;
                       gtree[i].nodes[p->index].father = gtree[i].root;
          }
             
     }
  

 
  Toclocktree(&gtree[i],gtree[i].root);
 
    tl = TreeL(file);
    tlg = 0.0;
    for(j=0;j<nspecies*2-1;j++) tlg += gtree[i].nodes[j].brlens;
 
  for(j=0;j<nspecies*2-1;j++)
	{
		gtree[i].nodes[j].brlens /= (GeneMu*tlg/tl);
            gtree[i].nodes[j].age /= (GeneMu*tlg/tl);
	} 
  
}


void ToGenetree(Tree *file[],int nfile,double *GeneMu)
{
  	TreeNode *p;
  	int i,j;
  	double tlg,tl;

  	for(i=0;i<nfile;i++)
  	{
		if(file[i]->isRooted != 1)
		{
			gtree[i].root = gtree[i].nTaxa*2-2;
   			gtree[i].nodes[gtree[i].root].sons[0] = file[i]->root->left->index;
   			gtree[i].nodes[gtree[i].root].sons[1] = file[i]->root->index;
   			gtree[i].nodes[gtree[i].root].nson = 2;
   			gtree[i].nodes[gtree[i].root].father = -1;
   			gtree[i].nodes[gtree[i].root].brlens = 0.0;

   			for(j=0;j<file[i]->nNodes;j++)
     			{ 
       				p=file[i]->allDownPass[j];
       				if(p != file[i]->root)
          			{
              				if(p->anc != NULL)
                     			{
                       				if(p != file[i]->root->left) 
                           			{
                              			gtree[i].nodes[p->index].brlens = p->length;
                              			gtree[i].nodes[p->index].father = p->anc->index;
                           			}
                      				else 
						{
                            			gtree[i].nodes[p->index].brlens = p->length/2;
                            			gtree[i].nodes[p->index].father = gtree[i].root;
                           			}
                     			}	
              				else 
						gtree[i].nodes[p->index].father = -1;
              				if(p->left != NULL)
                    			{ 
                      			gtree[i].nodes[p->index].nson = 2;
                     			gtree[i].nodes[p->index].sons[0] = p->left->index;
                     			gtree[i].nodes[p->index].sons[1] = p->right->index;
                    			}
               				else
                    			{ 
                       			gtree[i].nodes[p->index].nson = 0;
                       			gtree[i].nodes[p->index].sons[0] = -2;
                       			gtree[i].nodes[p->index].sons[1] = -2;
                    			}
          			}
       				else
          			{
                       		gtree[i].nodes[p->index].nson = 0;
                       		gtree[i].nodes[p->index].sons[0] = -2;
                       		gtree[i].nodes[p->index].sons[1] = -2;
                       		gtree[i].nodes[p->index].brlens = p->left->length/2;
                       		gtree[i].nodes[p->index].father = gtree[i].root;
          			}		
         		}
			/*for(j=0;j<file[i]->nNodes;j++)
                                printf("node %d %d %d %d %lf\n",j,gtree[i].nodes[j].sons[0],gtree[i].nodes[j].sons[1],gtree[i].nodes[j].father,gtree[i].nodes[j].brlens);*/
		}
		else
		{
			for(j=0;j<file[i]->nNodes-1;j++)
                	{
 				p=file[i]->allDownPass[j];
                        	gtree[i].nodes[p->index].brlens = p->length;
				gtree[i].nodes[p->index].father = p->anc->index;

                        	if(p->left != NULL)
				{
                                gtree[i].nodes[p->index].sons[0] = p->left->index;
				gtree[i].nodes[p->index].sons[1] = p->right->index;
				gtree[i].nodes[p->index].nson = 2;
				}
				else
                                {
                                gtree[i].nodes[p->index].nson = 0;
                                gtree[i].nodes[p->index].sons[0] = -2;
                                gtree[i].nodes[p->index].sons[1] = -2;
                                }
			}
			p=file[i]->allDownPass[file[i]->nNodes-1];
			gtree[i].root = p->left->index;
			gtree[i].nodes[gtree[i].root].father = -1;

			/*for(j=0;j<file[i]->nNodes-1;j++)
				printf("node %d %d %d %d %lf\n",j,gtree[i].nodes[j].sons[0],gtree[i].nodes[j].sons[1],gtree[i].nodes[j].father,gtree[i].nodes[j].brlens);*/
		}  	
	}
  	
	for(i=0;i<nfile;i++) Toclocktree(&gtree[i],gtree[i].root);
  
	/*normalized by total branch length*/
 	if (!strcmp(modelParams[0].brlensPr,"Clock"))
      	{
		for(i=0;i<nfile;i++)
                {
                                for(j=0;j<gtree[i].nTaxa*2-1;j++)
                                {
                                gtree[i].nodes[j].brlens /= GeneMu[i];
                                gtree[i].nodes[j].age    /= GeneMu[i];
                                }
		}
      	}
	else
	{
  		for(i=0;i<nfile;i++)
   		{
    			tl = TreeL(file[i]);
    			tlg = 0.0;
    			for(j=0;j<gtree[i].nTaxa*2-1;j++) tlg += gtree[i].nodes[j].brlens;
    			for(j=0;j<gtree[i].nTaxa*2-1;j++)
			{
			gtree[i].nodes[j].brlens /= (GeneMu[i]*tlg/tl);
			gtree[i].nodes[j].age   /= (GeneMu[i]*tlg/tl);
			}
   		}
	}

}

long int ClockTreeNodeDist(SPTree *clocktree, int ngene, Distance *dist) {
	int i, j, k, w, indexa, indexb, son0, son1;
	long int indexc=0;
	int taxa0[clocktree[0].nTaxa], taxa1[clocktree[0].nTaxa];

	for(w=0; w<ngene; w++) {
		for(j=clocktree[w].nTaxa; j<2*clocktree[w].nTaxa-1;j++) {
			son0 = clocktree[w].nodes[j].sons[0];
			son1 = clocktree[w].nodes[j].sons[1];
			indexa = 0; 
			FindDescendantTaxa(&clocktree[w], son0, taxa0, &indexa);
			indexb = 0;
			FindDescendantTaxa(&clocktree[w], son1, taxa1, &indexb);
		
			for(i=0; i<indexa; i++) for(k=0; k<indexb; k++)
				if(spnode[taxa0[i]] != spnode[taxa1[k]]) { //if from different species
					dist[indexc].nodes[0] = spnode[taxa0[i]];
					dist[indexc].nodes[1] = spnode[taxa1[k]];
					dist[indexc].dist = 2*(clocktree[w].nodes[j].age);
					indexc++;
				}
		}
	}
	return (indexc);
}

long int OneClockTreeNodeDist(SPTree *clocktree, Distance *dist)
{
	int i, j, k, indexa, indexb, son0, son1;
	long int indexc=0;
	int taxa0[NSPECIES], taxa1[NSPECIES]; 

	for(j=clocktree->nTaxa; j<2*clocktree->nTaxa-1;j++)
		{	
			son0 = clocktree->nodes[j].sons[0];
			son1 = clocktree->nodes[j].sons[1];
			indexa = 0; 
			FindDescendantTaxa(clocktree, son0, taxa0, &indexa);
			indexb = 0;
			FindDescendantTaxa(clocktree, son1, taxa1, &indexb);
		
			for(i=0; i<indexa; i++)
				for(k=0; k<indexb; k++)
				{
					if(spnode[taxa0[i]] != spnode[taxa1[k]])
					{
						dist[indexc].nodes[0] = spnode[taxa0[i]];
						dist[indexc].nodes[1] = spnode[taxa1[k]];
						dist[indexc].dist = 2*(clocktree->nodes[j].age);
						indexc++; 
					}
				}
		}
	return (indexc);
}

void OutNodes(treenode n,int id,int lvl) {
	for(int i=0; i<lvl; i++) printf("  ");
	printf("NODE %d: Age=%f, %d sons, species %d\n",id,n.age,n.nson,spnode[id]);
}
void PrintSPTree(SPTree *s,int nid,int lvl=0) {
	OutNodes(s->nodes[nid],nid,lvl);
	if(s->nodes[nid].nson>0) PrintSPTree(s,s->nodes[nid].sons[0],lvl+1);
	if(s->nodes[nid].nson>1) PrintSPTree(s,s->nodes[nid].sons[1],lvl+1);
}

void PrintMinMat(double **md, int nsp) {
	for(int i=0; i<nsp; i++) { printf("\n");
	for(int j=0; j<nsp; j++) printf("%05f ",md[0][i*nsp+j]);
	}
	printf("\n");
}

/*Scans across the ngene clocktrees to find mindist between all species*/
long int GetMinDists(SPTree *clocktree, int ngene, double **md) {
	bool trace=0;
	int i, j, k, w, nR, nL, son0, son1, nsp=sptree.nSpecies;
	long int indexc=0;
	int taxa0[clocktree[0].nTaxa], taxa1[clocktree[0].nTaxa];
	for(i=0; i<nsp;i++) for(j=0;j<nsp;j++) md[0][i*nsp+j]=0.0;

	for(w=0; w<ngene; w++) {
		if(trace) { printf("\nGene %d\n",w); PrintSPTree(&(clocktree[w]),clocktree[0].root); }
		for(j=clocktree[w].nTaxa; j<2*clocktree[w].nTaxa-1;j++) { //for each internode
			son0 = clocktree[w].nodes[j].sons[0];
			son1 = clocktree[w].nodes[j].sons[1];
			nR = 0;
			FindDescendantTaxa(&clocktree[w], son0, taxa0, &nR);
			nL = 0;
			FindDescendantTaxa(&clocktree[w], son1, taxa1, &nL);
			if(trace) for(i=0; i<nR; i++) for(k=0; k<nL; k++) printf("%dx%d ",taxa0[i],taxa1[k]);
			
			for(i=0; i<nR; i++) for(k=0; k<nL; k++)
				if(spnode[taxa0[i]] != spnode[taxa1[k]]) { //if from different species
				   if(spnode[taxa0[i]]<spnode[taxa1[k]]) indexc=spnode[taxa0[i]]*nsp+spnode[taxa1[k]];
				   else indexc=spnode[taxa1[k]]*nsp+spnode[taxa0[i]]; //figure out where to put this min
				   if((w==0 && md[0][indexc]==0.0) || 2*(clocktree[w].nodes[j].age)<md[0][indexc]) {
						md[0][indexc] = 2*(clocktree[w].nodes[j].age);
						if(trace) PrintMinMat(md,nsp);
						if(md[0][indexc]<.000001) printf("\n0 BL ERROR: Gene %d %d.%d",w,indexc/nsp,indexc%nsp);
					}
				}
		}
	}
	if(trace) for(i=0;i<nsp;i++) { printf("\n");
		for(j=i+1;j<nsp;j++) {
			printf("%05f ",md[0][i*nsp+j]);
			if(md[0][i*nsp+j]<.000001) system("PAUSE");
		}
	}
}




double TreeL(Tree *t)
{
	int			i;
	double		tl;
	TreeNode		*p;

	tl = 0.0;
	for (i=0; i<t->nNodes; i++)
	{
		p = t->allDownPass[i];
		if (p->anc != NULL)
		{
			if (p->anc->anc == NULL)
			{
				if (t->isRooted == NO)
					tl += p->length;
			}
			else
			{
				tl += p->length;
			}
		}
	}				
	return (tl);	
}

double Toclocktree(SPTree *t, int node)
{
  	double a,b;
  
  
  	if(t->nodes[node].nson==0)
     	{
         t->nodes[node].age = 0.0;
         return(0.0);
     	}
  	else
     	{   
         a = Toclocktree(t,t->nodes[node].sons[0]) + t->nodes[t->nodes[node].sons[0]].brlens;
         b = Toclocktree(t,t->nodes[node].sons[1]) + t->nodes[t->nodes[node].sons[1]].brlens;
         
         if(a>b)
            { 
               t->nodes[t->nodes[node].sons[1]].brlens += (a-b);
               t->nodes[node].age = a;
            }
         else if(a<b)
             { 
               t->nodes[t->nodes[node].sons[0]].brlens += (b-a);
               t->nodes[node].age = b;
             }
         else
             {
               t->nodes[node].age = a;
             }
         return(t->nodes[node].age);
     	}
}

int populationMutation (Tree *genetree,SPTree *speciestree, MrBFlt genemu)
{  
	int inode, i, k, inode_gene, stop, index=0;
	int *genetreenodes;
	TreeNode *p=NULL;

	genetreenodes = (int*)malloc(2*sptree.nTaxa*sizeof(int));
	if(!genetreenodes)
		{
		printf("allocating problem for genetreenodes\n");
		return(ERROR);
		}
	FOR(i,sptree.nSpecies)
      { 
        	FOR(inode,sptree.speciesIndex[i][0])
        	{        
        		inode_gene = sptree.speciesIndex[i][inode+1];
			FOR(k,(genetree->nNodes-1))
				{
				if(genetree->allDownPass[k]->index == inode_gene)
					{
					p = genetree->allDownPass[k];
					break;
					}
				}
   
        		stop=0;
        		do{
				/*check if the node is already taken care of*/
				FOR(k,index)
              			if(p->index == genetreenodes[k]) {stop=1;break;}
				if(stop == 1)
					break;
				/*change the branch length of node p*/				
				p->length = ChangeBrlen(speciestree, i, genetree, p);

				/*copy p to genetreenode*/
				genetreenodes[index] = p->index;
                 		index++;

				/*reset p*/
				p = p->anc;
			}while(p->index != genetree->root->left->index);
		}
	} 

	/*genes may have different mutation rates affecting all branches*/
	FOR(k,(genetree->nNodes-1))
	{
		p = genetree->allDownPass[k];
		p->length *= genemu;
	}
   	
    	if(index != 2*sptree.nTaxa-2)
    	{
      	printf("Error for the coalescence time calculation\n");
       	return(ERROR);
    	}
      free(genetreenodes);
	return(NO_ERROR);
}
				
double ChangeBrlen(SPTree *speciestree, int spnode, Tree *genetree, TreeNode *p)
{
	double length;
	int 	inode_sp, jnode_sp, father;

	/*find the node in the species tree that is immidiately down the gene node p*/
    	inode_sp = FindSpnodeDownGenenode(speciestree,spnode,p);
	jnode_sp = FindSpnodeDownGenenode(speciestree,spnode,p->anc);

	if(inode_sp == jnode_sp)
	{
		length = (p->length) * (speciestree->nodes[inode_sp].mu);
	}
	else
	{
		father = speciestree->nodes[inode_sp].father;
		length = (speciestree->nodes[father].age - p->nodeDepth)*(speciestree->nodes[inode_sp].mu);
		while(father != jnode_sp)
		{
			inode_sp = father;
			father = speciestree->nodes[father].father;
			length += (speciestree->nodes[father].age - speciestree->nodes[inode_sp].age)*(speciestree->nodes[inode_sp].mu);
		}
		length += (p->anc->nodeDepth - speciestree->nodes[father].age)*(speciestree->nodes[father].mu);
  	}
	return(length);		
}

int FindSpnodeDownGenenode(SPTree *speciestree, int spnode, TreeNode *p)
{
	int father, findnode=spnode;
	double depth;

	depth = p->nodeDepth;
	father = speciestree->nodes[spnode].father;

	if(p->index < sptree.nTaxa)
		return findnode;

	while(speciestree->nodes[father].age <= depth)
	{
		if(father == speciestree->root)
		{
			findnode = father;
			break;
		}
		else
		{
			findnode = father;
			father = speciestree->nodes[father].father;
		}
	}
	
	return (findnode);

}
			
int LnJointGenetreePr(Tree *t[],int num_tree, double *lncoalPrior, double *GeneMu, SPTree *speciestree)
{ 	
    	double	lnLike, lnPrior; 
 	int 		k, numchange;

	/*convert non-clock gene trees to clock gene trees*/
    	ToGenetree(t,num_tree,GeneMu);

	/*generate a random species tree*/
	numchange = poisson(1/poissonMean);
	if(StartSptree(&sptree, numchange)==ERROR)
	   {
		   printf(" Problem with startSptree \n");
		   return ERROR;

	   }  

   /*Calculate likelihood*/
	if(modelParam.thetainvgamma == 1) {
		if(SPLogLike_invgamma(gtree,&sptree,&lnLike) == ERROR) {
     		printf("Error for Lnlike\n");
     		return ERROR;
      }
	} else {
		if(SPLogLike(gtree,&sptree,&lnLike) == ERROR)
          	{
           		printf("Error for Lnlike\n");
           		return ERROR;
      	}
	}

	/*Calculate prior probability*/
    	if(SPLogPrior(&sptree,&lnPrior) == ERROR)
    	{
     		printf("Error for LnPrior\n");
           	return ERROR;
    	}                                
 
	/*Copy the species tree*/
	FOR(k, (2*sptree.nSpecies-1))
   	{
         	speciestree->nodes[k].nson = sptree.nodes[k].nson;
          	speciestree->nodes[k].sons[0] = sptree.nodes[k].sons[0];
          	speciestree->nodes[k].sons[1] = sptree.nodes[k].sons[1];
		speciestree->nodes[k].father = sptree.nodes[k].father;
          	speciestree->nodes[k].brlens = sptree.nodes[k].brlens;
          	speciestree->nodes[k].age = sptree.nodes[k].age;
		speciestree->nodes[k].theta = sptree.nodes[k].theta;
		speciestree->root = sptree.root;    
		speciestree->nSpecies = sptree.nSpecies;
		/*if(k == speciestree->root) speciestree->nodes[k].theta = 0.0;*/
   	}

	/*lnproposal = (-poissonMean + numchange * log(poissonMean));
	for(k = 1; k <= numchange; k++)
		lnproposal -= log(k); */

	*lncoalPrior = (lnLike + lnPrior);
    	return NO_ERROR;		       
}

int SPPrintTreeTitle (int curGen, FILE *fout)
{

	int				i;
	char			name[200], spfilename[100];
	FILE			*sumt;

	/* print the translate block information and the top of the file */
	if (curGen == 1)
		{
		/* print #NEXUS and translation block information */
		fprintf (fout, "#NEXUS\n");
		fprintf (fout, "[ID: %ld]\n", seed);
		fprintf (fout, "begin trees;\n");
		fprintf (fout, "   translate\n");
		for (i=0; i<sptree.nSpecies-1; i++)
			{
			GetNameFromString (taxaSetNames, name, i+1);			
			fprintf (fout, "      %2d %s,\n", i+1, name);
			}
		GetNameFromString (taxaSetNames, name, i+1);
		fprintf (fout, "      %2d %s;\n", i+1, name);
		}

	/* write the tree in Newick format */
	fprintf (fout, "   tree rep.%d = ", curGen);
	
	/* create a file to summarize species trees in the .sptree files */
	if (curGen == 1)
		{
		sprintf(spfilename, "%s.sumt", chainParams.chainFileName);
		sumt = fopen(spfilename, "w");
		sprintf(spfilename, "%s.sptree", chainParams.chainFileName);
	
		/* print #NEXUS and translation block information */
		fprintf (sumt, "#NEXUS\n");
		fprintf (sumt, "Begin data;\n");
		fprintf (sumt, "Dimensions ntax = %d nchar = 1;\n", sptree.nSpecies);
		fprintf (sumt, "FORMAT DATATYPE = dna gap=-  MISSING=?  interleave=yes;\nMatrix\n");
		for (i=0; i<sptree.nSpecies; i++)
			{
			GetNameFromString (taxaSetNames, name, i+1);			
			fprintf (sumt, " %s A\n", name);
			}
		fprintf (sumt, "\n;\nEND;\n");
		fprintf (sumt, "Begin mrbayes;\n  prset best = 1;\n");
		fprintf (sumt, "  sumt burnin = %d nruns = %d filename = %s contype = allcompat;\nEND;\n", chainParams.numGen/chainParams.sampleFreq/2, chainParams.numRuns, spfilename);
		
		fclose (sumt);
		}
 		   		
   	return (NO_ERROR);		
}


int SPAddToPrintString (char *tempStr)
{
	size_t			len1, len2;
	
	len1 = (int) strlen(printString);
	len2 = (int) strlen(tempStr);
	if (len1 + len2 + 5 > printStringSize)
		{
		printStringSize += len1 + len2 - printStringSize + 200;
		printString = (char*)realloc((void *)printString, printStringSize * sizeof(char));
		if (!printString)
			{
			SPMrBayesPrint ("%s   Problem reallocating printString (%d)\n", spacer, printStringSize * sizeof(char));
			goto errorExit;
			}
		}
	strcat(printString, tempStr);	
#	if 0
	printf ("printString(%d) -> \"%s\"\n", printStringSize, printString);
#	endif	
	return (NO_ERROR);
	
	errorExit:
		return (ERROR);
}

void SPPreparePrintFiles (void)
{
	char		localFileName[100], fileName[100];

	/* Get root of local file name */
	strcpy (localFileName, mcmc.chainFileName);

	/* Prepare the .p, .t */
	sprintf (fileName, "%s.p", localFileName);
	fpparm =(FILE*)gfopen(fileName,"w");
       sprintf (fileName, "%s.t", localFileName);
	fptree =(FILE*)gfopen(fileName,"w");
}

void InitiateParam(void)
{ 
 	modelParam.sRprior[0]=0;
  	modelParam.sRprior[1]=10;
  	modelParam.eRprior[0]=0;
  	modelParam.eRprior[1]=10;
  	modelParam.sF=1.0;
  	modelParam.treeHeightExp = 1.0;
}

int SPLogPrior(SPTree *speciestree, double *lnprior)
{
   	/* prior for the species tree*/
   	if(speciestreePr == 0)
		*lnprior = 0.0;
   	else
 	{
  		if (SPLnBirthDeathPriorPr (lnprior, speciationR, extinctionR, sampleF,speciestree) == ERROR)
   		{	
      		printf ("   Problem calculating prior for birth-death process\n");
      		return(ERROR);
   		}
		*lnprior += log(modelParam.treeHeightExp) - modelParam.treeHeightExp*speciestree->nodes[speciestree->root].age;
	}
   	return(NO_ERROR);
} 

void SPMrBayesPrint (char *format, ...)
{

	va_list                 ptr;

	va_start (ptr, format);

	vprintf (format, ptr);
	fflush (stdout);

	va_end(ptr);
}


int ReadaTree (FILE *fTree,SPTree *tree)
{
/* Read a tree from fTree, using the parenthesis node representation of trees.
   Branch lengths are read in nodes[].branch, and branch (node) labels 
   (integers) are preceeded by # and read in nodes[].label.  If the clade label
   $ is used, the label is read into nodes[].divtime first and then moved into
   nodes[].label in the routine DownTreeCladeLabel().

   Both names and numbers for species are accepted.  
   Species names are considered case-sensitive, with trailing blanks ignored.

   copyname = 0: species numbers and names are both accepted, but names have to match 
                 the names in taxaName[], which are from the sequence data file.
              1: species names are copied into taxaName[], but species numbers
                 are accepted.  
              2: the tree can have species names only, which are copied into 
                 taxaName[].
*/
   int cnode, cfather=-1;  /* current node and father */
   int inodeb=0;  /* node number that will have the next branch length */
   int i, level=0, ch=' ';
   char  skips[]="\"\'";
   int nnode;   

   nnode=tree->nTaxa;  
   FOR(i,2*tree->nTaxa-1) {
      tree->nodes[i].father=-1;
      tree->nodes[i].nson=0; 
   }


   while(isspace(ch)) 
     {
      ch=fgetc(fTree);  
     }
   ungetc(ch,fTree);

   for (;;) {
      ch = fgetc (fTree);
      if (ch==EOF) return(-1);
      else if (!isgraph(ch) || ch==skips[0] || ch==skips[1]) continue;
      else if (ch=='(') {
         level++;
         cnode=nnode++;
   
         if(nnode>2*tree->nTaxa-1)
              {
                  printf("check tree: perhaps too many '('s");
                  exit(-1);
              }
         if (cfather>=0) {
            tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = cnode;
            tree->nodes[cnode].father=cfather;
         }
         else
            tree->root=cnode;
         cfather=cnode;
      }
      else if (ch==')') { level--;  inodeb=cfather; cfather=tree->nodes[cfather].father; }
      else if (ch==':') fscanf(fTree,"%lf",&tree->nodes[inodeb].brlens); 
      else if (ch==',') ;
      else if (ch==';' && level!=0) 
         {
            printf("; in treefile");
            exit(-1);
         }
      else if (isdigit(ch))
      { 
         ungetc(ch, fTree); 
         fscanf(fTree,"%d",&inodeb); 
         inodeb--;
         tree->nodes[inodeb].father=cfather;
         tree->nodes[cfather].sons[tree->nodes[cfather].nson++]=inodeb;
      }
      if (level<=0) break;
   }
   
   for ( ; ; ) {
      while(isspace(ch=fgetc(fTree)) && ch!=';' ) 
         ;
      if (ch==':')       fscanf(fTree, "%lf", &tree->nodes[tree->root].brlens);
      else if (ch==';')  break;
      else  { ungetc(ch,fTree);  break; }
   }

   if(nnode!=2*tree->nTaxa-1) { printf(" # of nodes != %d\n",2*tree->nTaxa-1); exit(-1);}

   return (0);
}


int ReadControlfile(FILE *fdata)
{ 

  	int i,j,k,m,*species,max;

	 
  	/* parameters for MCMC */
	seed=swapSeed;
   	SetSeed(seed);

   	/* parameters for prior */
	if(!strcmp(modelParams[0].thetaPr,"Gamma"))
  		modelParam.thetainvgamma = 0;
	else
		modelParam.thetainvgamma = 1;
	modelParam.thetaprior[0] = modelParams[0].thetaGamma[0];
	modelParam.thetaprior[1] = modelParams[0].thetaGamma[1];	

	/* get information for gene trees  */
  	nGene = numCurrentDivisions;  
  	gtree = (SPTree*)malloc(nGene*sizeof(SPTree)); 
  	FOR(i,nGene) 
  	{
	  	gtree[i].nTaxa = numTaxa;     
	  	gtree[i].nodes =(Treenode*)malloc((2*gtree[i].nTaxa-1)*sizeof(Treenode));
      		if(!gtree[i].nodes)
      		{ 
		   	printf(" allocating problem for gtree.nodes\n");
           		return(ERROR);
      		}
	  	gtree[i].taxaName = (char **)calloc(gtree[i].nTaxa , sizeof(char *));
	  	gtree[i].taxaName[0]=(char *)calloc(gtree[i].nTaxa*LSPNAME,sizeof(char));

	  	for(j = 0; j < gtree[i].nTaxa; j++) 
			gtree[i].taxaName[j] = gtree[i].taxaName[0]+j*LSPNAME;

	  	if(!gtree[i].taxaName)
	  	{
	   		printf(" allocating problem for spname\n");
	   		return(ERROR);
	  	} 
     
  	}

  	FOR(i,nGene)
      		FOR(j,gtree[i].nTaxa) 
			GetNameFromString (taxaNames, gtree[i].taxaName[j], j+1);

  	FOR(i,nGene) 
	{
		if (!strcmp(modelParams[i].ploidy,"Haploid"))
			gtree[i].haploid = 1;
		else if (!strcmp(modelParams[i].ploidy,"Zlink"))
			gtree[i].haploid = -1;
		else
			gtree[i].haploid = 0;
	}

  	/* allocate space for the species tree */
   	sptree.nSpecies = numTaxaSets;
   	sptree.nodes = (Treenode*)malloc((2*sptree.nSpecies-1)*sizeof(Treenode));
   	if(!sptree.nodes)
   	{
       		printf("allocating problem for sptree.nodes\n");
       		return(ERROR);
   	}

   	species = (int*)malloc(sptree.nSpecies*sizeof(int));
	sptree.nTaxa = 0;
 	FOR(i,sptree.nSpecies)  
        {
	   species[i] = 0;
	   FOR(k,numTaxa)
		{
		if(taxaInfo[k].taxaSet[i] == 1)
			species[i]++;
		}
	  
	   sptree.nTaxa += species[i];
	}

	if(sptree.nTaxa != gtree[0].nTaxa)
	{
		printf("Wrong Taxaset!\n");
                return(ERROR);
        }

   	spnode = (int*)malloc(sptree.nTaxa*sizeof(int));

   	sptree.taxaName = (char **)calloc(sptree.nTaxa,sizeof(char *));
   	sptree.taxaName[0]=(char *)calloc(sptree.nTaxa*LSPNAME,sizeof(char));


	for(j = 0; j < sptree.nTaxa; j++) 
		sptree.taxaName[j] = sptree.taxaName[0]+j*LSPNAME;
	if(!sptree.taxaName)
	  {
	   	printf(" allocating problem for spname\n");
	   	return(ERROR);
	  } 

   	max = 0;
   	FOR(i,sptree.nSpecies) 
         	if(max<species[i]) 
			max = species[i];
   	max++;

   	sptree.speciesIndex = (int**)calloc(sptree.nSpecies,sizeof(int*));
   	if(!sptree.speciesIndex)
	   {
	   printf(" allocating problem for sptree.speciesIndex\n");
	   return(ERROR);
	   } 
   	sptree.speciesIndex[0] = (int*)calloc(sptree.nSpecies*max,sizeof(int));
   	FOR(i,sptree.nSpecies)
        	sptree.speciesIndex[i] = sptree.speciesIndex[0] + i*max;

   	FOR(i,sptree.nSpecies)
   	{
      		sptree.speciesIndex[i][0] = species[i];
		m = 1;
		FOR(k,numTaxa)
		{
			if(taxaInfo[k].taxaSet[i] == 1)
			{
				spnode[k] = i;
				sptree.speciesIndex[i][m] = k;
				m++;
			}
		}

		/*FOR(k,sptree.speciesIndex[i][0])
			printf("index %d ",sptree.speciesIndex[i][k+1]);
		printf("\n");*/
   	}  
   	FOR(i,sptree.nTaxa)   
		GetNameFromString (taxaNames, sptree.taxaName[i], i+1);		

	sptree.nPop = 0;
   	FOR(i,sptree.nSpecies)
       		if(species[i]>1) 
			sptree.nPop++;
   	sptree.nPop += (sptree.nSpecies-1);

   	sptree.popIndex = (int*)malloc(sptree.nPop*sizeof(int));
   	j = 0;
   	FOR(i,sptree.nSpecies)
       		if(species[i]>1) 
          	{
               		sptree.popIndex[j] = i;
               		j++;
          	}

   	for(i=1;i<=sptree.nSpecies-1;i++)
       		sptree.popIndex[j++] = sptree.nSpecies-1+i;
   
   	sptree.nconstraint = sptree.nSpecies*(sptree.nSpecies-1)/2;

   	sptree.constraint = (double*)malloc(sptree.nconstraint*sizeof(double));
   	
	if(!sptree.constraint)
       	{
          	printf("allocating problem for sptree.constraint\n");
          	return(ERROR);
       	}

   	sptree.treeConstraint = (Distance*)malloc((sptree.nSpecies-1)*sizeof(Distance));
   	
	if(!sptree.treeConstraint)
	{
		printf("allocating problem for treecontraint\n");
		return(ERROR);
	}

  	/* allocating space for gene tree */

       	FOR(i,nGene)
        {


           gtree[i].taxaIndex = (int*)malloc(sptree.nTaxa*sizeof(int));
           if(!gtree[i].taxaIndex)
            {
                printf("allocating problem for gtree.taxaIndex\n");
                return(ERROR);
             }


           gtree[i].speciesIndex = (int**)calloc(sptree.nSpecies,sizeof(int*));
           if(!gtree[i].speciesIndex)
	   {
	   printf(" allocating problem for gtree[i].speciesIndex\n");
	   return(ERROR);
	   } 
           gtree[i].speciesIndex[0] = (int*)calloc(sptree.nSpecies*max,sizeof(int));
           FOR(j,sptree.nSpecies)
                gtree[i].speciesIndex[j] = gtree[i].speciesIndex[0] + j*max;


           gtree[i].popIndex = (int*)malloc(sptree.nPop*sizeof(int));
           if(!gtree[i].popIndex)
	   {
	   printf(" allocating problem for gtree[i].popIndex\n");
	   return(ERROR);
	   } 
        }

	FOR(i,nGene)
        	FOR(j,sptree.nTaxa)
         		gtree[i].taxaIndex[j] = j;

	FOR(i,nGene)
        {
 		gtree[i].nSpecies = sptree.nSpecies;
           	gtree[i].nPop = sptree.nPop;

		/*copy species index*/
           	FOR(j,sptree.nSpecies)
           	{
             		gtree[i].speciesIndex[j][0] = sptree.speciesIndex[j][0];
             		for(m=1;m<sptree.speciesIndex[j][0]+1;m++)           
             			gtree[i].speciesIndex[j][m] = sptree.speciesIndex[j][m];

             	}
		/*copy population index*/  
        	for(j=0;j<sptree.nPop;j++)
              		gtree[i].popIndex[j] = sptree.popIndex[j];      
     
        }
  	free(species); 
  	return(NO_ERROR);  
}

void PrintInf(FILE *fp)
{
  int i,j,k;
    
  /* write information into log.txt */
  fprintf(fp,"%s      Setting MCMC parameters\n",spacer);
  fprintf(fp,"%s      seed = %ld;\n",spacer, seed); 
  fprintf(fp,"%s      Theta = gamma(%lf, %lf); \n",spacer, modelParam.thetaprior[0],modelParam.thetaprior[1]);
  fprintf(fp,"\n\n\n");

  /* print Taxa information */
  fprintf(fp,"%s      Taxa Matrix\n",spacer);
  fprintf(fp,"%s          \t",spacer);
  FOR(i,sptree.nTaxa) fprintf(fp,"%s\t",sptree.taxaName[i]);
  fprintf(fp,"\n");

  FOR(j,nGene)
  { 
      fprintf(fp,"%s      Gene%d\t",spacer, j+1);
      FOR(i,sptree.nTaxa)
        fprintf(fp,"%d\t",gtree[j].taxaIndex[i]+1);
      fprintf(fp,"\n");
  }
  fprintf(fp,"\n\n");



  /*print species information */
  fprintf(fp,"%s      Species Matrix in Species SPTree\n",spacer);
  fprintf(fp,"%s      Species\tnTaxa\tList\n",spacer);
  FOR(i,sptree.nSpecies) 
      {
         fprintf(fp,"%s      %d.\t%d\t", spacer, i+1, sptree.speciesIndex[i][0]);
         FOR(j,sptree.speciesIndex[i][0])
          fprintf(fp,"%d\t",sptree.speciesIndex[i][j+1]+1);
         fprintf(fp,"\n");
      }
  fprintf(fp,"\n");
 
  FOR(j,nGene)
  {
      fprintf(fp,"%s      Species Matrix in Gene Tree%d\n",spacer, j+1);
      fprintf(fp,"%s      Species\tnTaxa\tList\n",spacer);
      FOR(i,sptree.nSpecies) 
       {
            fprintf(fp,"%s      %d.\t%d\t",spacer, i+1,gtree[j].speciesIndex[i][0]);
         FOR(k,gtree[j].speciesIndex[i][0])
          fprintf(fp,"%d\t",gtree[j].speciesIndex[i][k+1]+1);
         fprintf(fp,"\n");

        }
     fprintf(fp,"\n");
  }
  fprintf(fp,"\n\n");

  /*print population information*/
  fprintf(fp,"%s      Population Matrix in Species Tree\n",spacer);
  fprintf(fp,"%s      Species tree Populations:\t",spacer);
  FOR(i,sptree.nPop) 
      {
         fprintf(fp,"%d\t", sptree.popIndex[i]);
      }
  fprintf(fp,"\n");
   
  fprintf(fp,"%s      Gene tree Populations\n",spacer);
  FOR(j,nGene)
  {FOR(i,gtree[j].nPop) 
      {
         fprintf(fp,"%d\t", gtree[j].popIndex[i]);
      }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n\n");



}
  

/*
void FindMissingname(Tree *genetree)
{ 
  
	int i,j,k,index,stop,missingindex,totalname;
	char **name;

  
	totalname = 0;
	FOR(i, nGene) totalname += gtree[i].nTaxa;

	name = (char**)calloc(totalname*sizeof(char*));
  
	if(!name)
	{
	  printf("allocating problem for name\n");
	  return(ERROR);
	}

 	For(i,gtree[0].nTaxa) name[i] = gtree[0].taxaName[i];
 
	totalname = gtree[0].nTaxa;

	for(i=1;i<nGene;i++)
	{  
	 
		index = 0;
		FOR(j,gtree[i].nTaxa)
		{ 
		  find = 0;
          	for(k=0;k<totalname;k++)
		  {
                if(!strcmp(gtree[i].taxaName[j],name[k])) 
				{
					gtree[i].taxaIndex[index] = k;
					find = 1;
					index++;
					break;
				}
		  }               
		  if(find == 0) 
		  {			  
			  name[totalname] = gtree[i].taxaName[j]);
			  gtree[i].taxaIndex[index] = totalname;
			  index++;
			  totalname++;	             
		  }
		}
	}

	sptree.nTaxa = totalname;
	sptree.taxaName = (char **)malloc(sptree.nTaxa * sizeof(char *));
    	if(!sptree.taxaName)
	{
	  printf(" allocating problem for sptree.taxaName\n");
	  return(ERROR);
	} 
	for(i=0;i<sptree.nTaxa;i++)
		sptree.taxaName[i] = name[i];

	for(i=0;i<nGene;i++)
	{
		index = 0;
		gtree[i].nmissTaxa = sptree.nTaxa - gtree[i].nTaxa;
		gtree[i].mtaxaIndex = (int*)malloc(gtree[i].nmissTaxa*(int));
		for(j=0;j<sptree.nTaxa-1;j++)
		{
			k=0;
			if(j == gtree[i].taxaIndex[k])
			{
				k++;
				continue;
			}
			else
			{
				gtree[i].mtaxaIndex[index];
				index++;
			}
	}

        free(name);
	return(NO_ERROR);  
}
*/

void ReadGeneTree(FILE *fTree)
{
  int i;  

  /* Read in Gene Trees */
  FOR(i,nGene) ReadaTree(fTree,&gtree[i]);

 
  /* Calculate the age of nodes */
  FOR(i,nGene) CalNodeAge(gtree[i].root, &gtree[i]);


}
  
int SPPrintTree(int curGen, SPTree *tree, int showBrlens, int showTheta, int showMu, int isRooted)
{
    char			*tempStr;
	int             tempStrSize, i;
	FILE			*sumt;
	char			name[200], spfilename[100];

	/* allocate the print string */
	printStringSize = 200;
	printString = (char *)SafeMalloc((size_t) (printStringSize * sizeof(char)));
	if (!printString)
		{
		SPMrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
		return (ERROR);
		}
	*printString = '\0';

	tempStrSize = 200;
	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
		{
		SPMrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
		return (ERROR);
		}

	/*print tree file header*/
	if (curGen == 1)
	{
		/* print #NEXUS and translation block information */
		sprintf (tempStr, "#NEXUS\n");
		if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
		sprintf (tempStr, "[ID: %ld]\n", seed);
		if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
		sprintf (tempStr, "begin trees;\n");
		if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
		sprintf (tempStr, "   translate\n");
		if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
		for (i=0; i<sptree.nSpecies-1; i++)
			{
			GetNameFromString (taxaSetNames, name, i+1);			
			sprintf (tempStr, "      %2d %s,\n", i+1, name);
			if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);
			}
		GetNameFromString (taxaSetNames, name, i+1);
		sprintf (tempStr, "      %2d %s;\n", i+1, name);
        if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);

		/*create .sumt file*/
# if defined (MPI_ENABLED)
		if(proc_id == 0)
		{
# endif
		sprintf(spfilename, "%s.sumt", chainParams.chainFileName);
		sumt = fopen(spfilename, "w");
		sprintf(spfilename, "%s.sptree", chainParams.chainFileName);
	
		/* print #NEXUS and translation block information */
		sprintf(spfilename, "%s.sumt", chainParams.chainFileName);
		sumt = fopen(spfilename, "w");
		sprintf(spfilename, "%s.sptree", chainParams.chainFileName);
		fprintf (sumt, "#NEXUS\n");
		fprintf (sumt, "Begin data;\n");
		fprintf (sumt, "Dimensions ntax = %d nchar = 1;\n", sptree.nSpecies);
		fprintf (sumt, "FORMAT DATATYPE = dna gap=-  MISSING=?  interleave=yes;\nMatrix\n");
		for (i=0; i<sptree.nSpecies; i++)
			{
			GetNameFromString (taxaSetNames, name, i+1);			
			fprintf (sumt, " %s A\n", name);
			}
		fprintf (sumt, "\n;\nEND;\n");
		fprintf (sumt, "Begin mrbayes;\n  prset best = 1;\n");
		fprintf (sumt, "  sumt burnin = %d nruns = %d filename = %s contype = allcompat;\nEND;\n", chainParams.numGen/chainParams.sampleFreq/2, chainParams.numRuns, spfilename);
		fclose (sumt);
#if defined (MPI_ENABLED)
		}
# endif
	}

	/* write the tree in Newick format */
	sprintf (tempStr, "   tree rep.%d = ", curGen);
	if (SPAddToPrintString (tempStr) == ERROR) return(ERROR);

    SPSaveSprintf (&tempStr, &tempStrSize,"(");
	SPAddToPrintString (tempStr);
					
    SPWriteTreeToFile (tree, tree->root, showBrlens, showTheta, showMu, isRooted);

    if(showTheta == YES) 
		SPSaveSprintf (&tempStr, &tempStrSize,")[#%lf];\n",tree->nodes[tree->root].theta);
    else 
		SPSaveSprintf (&tempStr, &tempStrSize,");\n");
	SPAddToPrintString (tempStr);
	free (tempStr); 

	return (NO_ERROR);					
}

void SPWriteTreeToFile (SPTree *tree, int inode, int showBrlens, int showTheta, int showMu, int isRooted)

{
		char			*tempStr;
		int                      tempStrSize = 200;

		tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
		if (!tempStr)
			SPMrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));


			
		if (tree->nodes[inode].nson == 0)
			{
				if (showBrlens == YES)
				{
    				SPSaveSprintf (&tempStr, &tempStrSize, "%d:%lf", inode+1, tree->nodes[inode].brlens);
					SPAddToPrintString (tempStr);
					if((tree->nodes[inode].theta>0) && showTheta == YES) 
					{
						SPSaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[inode].theta);
						SPAddToPrintString (tempStr);
					}
				}
				else
				{
					SPSaveSprintf (&tempStr, &tempStrSize, "%d", inode+1);
					SPAddToPrintString (tempStr);
				}
			}
		else
			{
				if (inode != tree->root)
				{
					SPSaveSprintf (&tempStr, &tempStrSize, "(");
					SPAddToPrintString (tempStr);
				}
				SPWriteTreeToFile (tree,tree->nodes[inode].sons[0],  showBrlens, showTheta, showMu, isRooted);
				SPSaveSprintf (&tempStr, &tempStrSize, ",");
				SPAddToPrintString (tempStr);
				SPWriteTreeToFile (tree,tree->nodes[inode].sons[1], showBrlens, showTheta, showMu, isRooted);	
				if (inode != tree->root)
				{
					if (tree->nodes[inode].father == tree->root && isRooted == NO)
					{
						if (showBrlens == YES)
						{
							SPSaveSprintf (&tempStr, &tempStrSize, ",%d:%lf", tree->nodes[inode].father + 1, tree->nodes[tree->nodes[inode].father].brlens);
							SPAddToPrintString (tempStr);
							if((tree->nodes[tree->nodes[inode].father].theta>0) && showTheta == YES) 
							{
								SPSaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[tree->nodes[inode].father].theta);
								SPAddToPrintString (tempStr);
							}
						}
						else
						{
							SPSaveSprintf (&tempStr, &tempStrSize, ",%d", tree->nodes[inode].father + 1);
							SPAddToPrintString (tempStr);
						}
					}
				
					if (showBrlens == YES && isRooted == YES) /*tree->nodes[inode].father != tree->root)*/
					{
						SPSaveSprintf (&tempStr, &tempStrSize,"):%lf", tree->nodes[inode].brlens);
						SPAddToPrintString (tempStr);
						if((tree->nodes[inode].theta > 0) && showTheta == YES)
						{
							SPSaveSprintf (&tempStr, &tempStrSize, "[#%lf]", tree->nodes[inode].theta);
							SPAddToPrintString (tempStr);
						}
					}
					else
					{
						SPSaveSprintf (&tempStr, &tempStrSize, ")");
						SPAddToPrintString (tempStr);
					}					
				}
			}
	free (tempStr);
		
}


double SpeciesDistance(int inode, int jnode, SPTree *genetree)
{
  int i,j;
  double distance;
  double con = 0.0;

  if(genetree->speciesIndex[inode][0] == 0 || genetree->speciesIndex[jnode][0] == 0)
     return(NA);

  FOR(i,genetree->speciesIndex[inode][0])
     FOR(j,genetree->speciesIndex[jnode][0])
          {      
                 distance = NodeDistance(genetree,genetree->speciesIndex[inode][i+1],genetree->speciesIndex[jnode][j+1]);
                 if(j==0&&i==0) 
                 {
                     con = distance;
                     continue;
                 }
                 if(con>distance) con = distance;
           }
  return(con);
}

double NodeDistance(SPTree *tree, int inode, int jnode)
{
  int i, *ancestor,father,nancestor,stop=0;

  ancestor=(int*)malloc(tree->nTaxa*sizeof(int));

  i=0;
  father = inode;
  do{
      father = tree->nodes[father].father;
      ancestor[i] = father;
      i++;
    }while( father != tree->root);

  nancestor = i;

  father = jnode;
  do{
      father = tree->nodes[father].father;
      FOR(i,nancestor) 
         if(ancestor[i] == father) 
             {
                stop = 1;
                break;
              }
     }while(stop == 0);
  free(ancestor);
  return(2*tree->nodes[father].age);
}
 
int Constraint(SPTree *genetree, int numgenetree, SPTree *speciestree, Distance *constraint)
{
  	int i, j, w, k;
  	double distance;

  	w = 0;
  	FOR(j,(speciestree->nSpecies-1))
    		for(k=j+1;k<speciestree->nSpecies;k++)
     		{
			constraint[w].nodes[0] = j;
             	constraint[w].nodes[1] = k;
		 	constraint[w].dist = 100000.0;
      		FOR(i,numgenetree)
          		{  
             		distance = SpeciesDistance(j,k,&genetree[i]);
             		if(constraint[w].dist > distance) constraint[w].dist = distance;
          		}
      		if(constraint[w].dist <= 0.0 || constraint[w].dist >= 100000.0)
          		{
             		printf("constraint%d is wrong\n",w);
             		return(ERROR);
           		}
       		w++;
     		}
  	return(NO_ERROR);
}

double CalNodeAge(int node, SPTree *tree)
{ 
  if(tree->nodes[node].nson == 0) 
      {
         tree->nodes[node].age = 0.0;
         return(0.0);
      }
  else
      {
         tree->nodes[node].age = 0.5*(CalNodeAge(tree->nodes[node].sons[0],tree) + tree->nodes[tree->nodes[node].sons[0]].brlens + CalNodeAge(tree->nodes[node].sons[1],tree) + tree->nodes[tree->nodes[node].sons[1]].brlens);
         return( tree->nodes[node].age );
      }
}
       
int StartSptree(SPTree *speciestree, int numchange) { //*speciestree is the address of the global sptree
   int i, j;
   Distance onetreeConstraint[speciestree->nSpecies];

	//CNKA 10/10 we don't want all cophenetic distances for all genes, just the min across genes
	speciestree->mindist = (double*)malloc(speciestree->nSpecies*speciestree->nSpecies*sizeof(double));
   GetMinDists(gtree, nGene, &speciestree->mindist);
   GetConstraints(speciestree, onetreeConstraint);

	/*jiggle polytomys*/
  	for(i=0;i<sptree.nSpecies-2;i++)
		for(j=i+1;j<sptree.nSpecies-1;j++) {
			if(onetreeConstraint[i].dist == onetreeConstraint[j].dist) {
				onetreeConstraint[j].dist -= 1e-10;
				//if(onetreeConstraint[j].dist < 0.0) onetreeConstraint[j].dist = 0.0;
			}
		}

	//change constraints
	for(i=0; i<numchange; i++)
		ChangeConstraint(onetreeConstraint, speciestree->nSpecies-1);

  	//propose a species tree
  	MaximumTree(onetreeConstraint, speciestree, speciestree->nSpecies-1);

  	// initial theta
  	FOR(i, 2*speciestree->nSpecies-1) speciestree->nodes[i].theta = -1.0;
  	FOR(i, speciestree->nPop) speciestree->nodes[speciestree->popIndex[i]].theta = 0.1;
  	return(NO_ERROR); 
}

int CheckConstraint(SPTree *genetrees, int ngene, Distance *constraint, int nconstraints) {
	int i, j, k, l, node0, node1;
	double dist;

	FOR(i, nconstraints)
	{
		node0 = constraint[i].nodes[0];
		node1 = constraint[i].nodes[1];
		FOR(j, ngene)
		{
			FOR(k, sptree.speciesIndex[node0][0])
				FOR(l,sptree.speciesIndex[node1][0])
				{
					dist = NodeDistance(&genetrees[j],sptree.speciesIndex[node0][k+1],sptree.speciesIndex[node1][l+1]);
					if(dist <  constraint[i].dist)
					{
						printf("constrain %d %d %f genedist %d %d %d %f\n",node0,node1,constraint[i].dist,j,sptree.speciesIndex[node0][k+1],sptree.speciesIndex[node1][l+1],dist);
						return ERROR;
					}
				}
		}
	}
	return (NO_ERROR);
}

int ChangeConstraint(Distance *dist, int nconstraints)
{
	int i;
	double window=0.0001;

	i = (int)(rndu()*nconstraints);
	dist[i].dist -= rndu()*window;
	
  	if(dist[i].dist < 1e-10) 
     	{ 
      	dist[i].dist = 1e-10;
     	}
   	return(NO_ERROR);
}


int MaximumTree(Distance *dist, SPTree *speciestree, int nconstraints)
{
	int i, j, k, father, son;

	/*sort the distances*/	
	quick_struct(dist, nconstraints);

	for(i=0; i<2*(nconstraints+1)-1; i++)
	{
		speciestree->nodes[i].father = -2;
	}

	k = nconstraints + 1;
	for(i=0; i<nconstraints; i++)
	{
		speciestree->nodes[k].nson = 2;
		speciestree->nodes[k].age = dist[i].dist/2;
		for(j=0; j<2; j++)
		{
			speciestree->nodes[dist[i].nodes[j]].age = 0.0;
			speciestree->nodes[dist[i].nodes[j]].nson = 0;
			speciestree->nodes[dist[i].nodes[j]].sons[0] = -2;
			speciestree->nodes[dist[i].nodes[j]].sons[1] = -2;

			father = speciestree->nodes[dist[i].nodes[j]].father;
			son = dist[i].nodes[j];
			while(father != -2)
			{
				son = father;
				father = speciestree->nodes[father].father;
			} 
			speciestree->nodes[k].sons[j] = son;
			speciestree->nodes[son].father = k;
			speciestree->nodes[son].brlens = speciestree->nodes[k].age - speciestree->nodes[son].age; 
		}
		k++;
	}
	speciestree->root = k-1;
	speciestree->nodes[k-1].father = -1;
	return(NO_ERROR);
}

int SPTreeConstraint(Distance *minimumdistance, Distance *distance, long int nconstraints, int nspecies) {
	int i, k;
	int node[NSPECIES], index[2];
   long int j;
	
	quick_struct(minimumdistance, nconstraints);

  	distance[0].dist = minimumdistance[0].dist;
  	distance[0].nodes[0]=node[0]=minimumdistance[0].nodes[0];
  	distance[0].nodes[1]=node[1]=minimumdistance[0].nodes[1];
  
  	for(i=2; i<nspecies; i++)
  	{ 
	  	FOR(j, nconstraints)  
	  	{ 
		  	index[0]=0;index[1]=0;
		  	FOR(k,i)
		  	{
             		if (minimumdistance[j].nodes[0] == node[k]) {index[0]++;}
             		if (minimumdistance[j].nodes[1] == node[k]) {index[1]++;}
		  	}

		  	if((index[0]==0)&&(index[1]==1)) 
		  	{
			   	distance[i-1].nodes[0]=minimumdistance[j].nodes[1];
               distance[i-1].nodes[1]=minimumdistance[j].nodes[0];
               distance[i-1].dist=minimumdistance[j].dist;
               node[i]=minimumdistance[j].nodes[0];
               break;
		  	}        
		  	if((index[0]==1)&&(index[1]==0))
		  	{
			   	distance[i-1].nodes[0]=minimumdistance[j].nodes[0];
               distance[i-1].nodes[1]=minimumdistance[j].nodes[1];
               distance[i-1].dist = minimumdistance[j].dist;
               node[i] = minimumdistance[j].nodes[1];
               break;
		  	}         
	  	}  
  	}
	return(NO_ERROR);
}

//loads constraints on the species tree given that s->mindist is already set
int GetConstraints(SPTree *s, Distance *constr) {
	bool trace=0;
	int i,j,k=0,nsp=s->nSpecies;
	int node[nsp], index[2];

	//convert s->mindist to a distance structure
	Distance minimumdistance[nsp*(nsp-1)/2];
	for(i=0; i<nsp; i++) for(j=i+1;j<nsp;j++) {
		minimumdistance[k].dist=s->mindist[i*nsp+j];
		minimumdistance[k].nodes[0]=i;
		minimumdistance[k++].nodes[1]=j;
	}
	quick_struct(minimumdistance, nsp*(nsp-1)/2);

  	constr[0].dist = minimumdistance[0].dist;
  	constr[0].nodes[0]=node[0]=minimumdistance[0].nodes[0];
  	constr[0].nodes[1]=node[1]=minimumdistance[0].nodes[1];

  	for(i=2; i<nsp; i++)
  	{
	  	FOR(j, nsp*(nsp-1)/2)
	  	{
		  	index[0]=index[1]=0;
		  	FOR(k,i)
		  	{
       		if (minimumdistance[j].nodes[0] == node[k]) {index[0]++;}
       		if (minimumdistance[j].nodes[1] == node[k]) {index[1]++;}
		  	}

		  	if((index[0]==0)&&(index[1]==1)) {
		   	constr[i-1].nodes[0]=minimumdistance[j].nodes[1];
	         constr[i-1].nodes[1]=minimumdistance[j].nodes[0];
	         constr[i-1].dist=minimumdistance[j].dist;
	         node[i]=minimumdistance[j].nodes[0];
	         break;
		  	}
		  	if((index[0]==1)&&(index[1]==0)) {
		   	constr[i-1].nodes[0]=minimumdistance[j].nodes[0];
            constr[i-1].nodes[1]=minimumdistance[j].nodes[1];
            constr[i-1].dist = minimumdistance[j].dist;
            node[i] = minimumdistance[j].nodes[1];
            break;
		  	}
	  	}
  	}
  	if(trace) for(i=0;i<nsp-1;i++)
	  printf("\n%d. (%d.%d:%f)",i+1,constr[i].nodes[0],constr[i].nodes[1],constr[i].dist);
	return(NO_ERROR);
}


int FindaPosition(int nodenumber,int root,double currentdistance,double *position,SPTree *speciestree)  /*ancestor contains nodes and distance */
{
  

  position[0]=0.0;
  position[1]=0.0;

  while(nodenumber!=root)
  { 
	  position[1] += speciestree->nodes[nodenumber].brlens;
	  if (currentdistance<=position[1])
	  {
		  position[0]=nodenumber;
		  position[1]=currentdistance-(position[1]-speciestree->nodes[nodenumber].brlens);
		  break;
	  } 
	  if(currentdistance == position[1])
	  {
		  printf("Equal distance for position %f\n",currentdistance);
		  return(ERROR);
	  }
	  nodenumber=speciestree->nodes[nodenumber].father;
  }
  if(nodenumber==root)
  {
	  position[0]=root;
	  position[1]=currentdistance-position[1];
  }
  
  return(NO_ERROR);
}

void quick_struct(Distance item[],int count)
         {qs_struct(item,0,count-1);}

void qs_struct(Distance item[],int left,int right)
{

   register int i, j;
   double x;
   Distance temp;

  i=left; j=right;
  x=item[(left+right)/2].dist;


  do {
       while(item[i].dist<x && i<right) i++;
       while(item[j].dist>x && j>left) j--;

       if(i<=j)  {
          temp=item[i];
          item[i]=item[j];
          item[j]=temp;
          i++; j--;
          }
     } while(i<=j);


    if(left<j)  qs_struct(item,left, j);
    if(i<right)  qs_struct(item,i, right);
}

int SPLogLike(SPTree *genetree, SPTree *speciestree, double *lnl)
{
  	int k;
  	double lnp;

  	*lnl = 0.0;

  	FOR(k,nGene) 
    	{
       	if( LnLikehood1Tree(&genetree[k],speciestree, &lnp) == ERROR)
          	{
            	printf("Error for calculating the log likelhood of gene trees\n");
               	return(ERROR);
          	}
       	*lnl += lnp;
      }
  	return(NO_ERROR);
}

int LnLikehood1Tree(SPTree *genetree, SPTree *speciestree, double *lnp)
{ 
   	int i,j, ipop,k=0, nt;
   	double *t,y,theta;
   	CoalTime coal,coal_t;

  	if( GetCoaltime(genetree,speciestree,&coal) == ERROR)
     	{
         printf(" bugs in GetCoaltime \n");
         return(ERROR);
     	}
 
  	i=0;
  	GetNcoal(speciestree->root, &coal, genetree, speciestree, &coal_t, &i);
   
  	for(ipop=0; ipop<genetree->nPop; ipop++)
  	{      
	  nt=coal_t.ncoal[ipop]; 
	  t=coal_t.tj[ipop];   
	  for(i=0; i<nt-1; i++)  for (j=i+1; j<nt; j++)
	  {
            if (t[j]<t[i])  
			{ 
				y=t[i]; 
				t[i]=t[j];
				t[j]=y;
			}
	  }
  	}

  	*lnp = 0.0;
 	for(i=0; i<genetree->nPop; i++) 
	{
      	ipop=i; t=coal_t.tj[ipop];
         	/* handle the haploid gene*/
         	if(genetree->haploid == 1) /*haploid*/
              	theta = 0.25*speciestree->nodes[coal_t.nodes[ipop]].theta;
         	else if (genetree->haploid == -1) /*zlink*/
				theta = 0.75*speciestree->nodes[coal_t.nodes[ipop]].theta;
			else  /*diploid*/
              	theta = speciestree->nodes[coal_t.nodes[ipop]].theta;

#	if defined (DEBUG)
         	if(theta<1.e-6)
            {
               printf("theta is wrong for node %d  theta %lf\n",coal_t.nodes[ipop],theta);
               return(ERROR);
            }
#	endif
 
         	if(coal_t.nin[ipop]>coal_t.nout[ipop]) 
		{
         		*lnp+=(coal_t.nin[ipop]-coal_t.nout[ipop])*log(2/theta);
         		for(j=coal_t.nin[ipop],k=0; j>coal_t.nout[ipop]; j--,k++) 
			{
            	y = t[k]- (k==0? 0:t[k-1]);       
            	*lnp -= j*(j-1)/theta*y;
         		}
      	}
      	if(coal_t.nout[ipop]>1) 
		{
         	y = (coal_t.nout[ipop]==coal_t.nin[ipop] ? speciestree->nodes[coal_t.nodes[ipop]].brlens : speciestree->nodes[coal_t.nodes[ipop]].brlens-t[k-1]);
         	*lnp -= coal_t.nout[ipop]*(coal_t.nout[ipop]-1)/theta*y;
         
      	}
   	}  
   	return(NO_ERROR);
}

int SPLogLike_invgamma(SPTree *genetree, SPTree *speciestree, double *lnl)
{
  	int i,j,k;
  	double a[NSPECIES],b[NSPECIES], totalcoal=0.0, alpha=modelParam.thetaprior[0],beta=modelParam.thetaprior[1];
  

  	/*initial a and b*/
  	for(i=0; i<speciestree->nPop;i++)
	{
		a[i] = 0.0;
 		b[i] = 0.0;
	}

  	*lnl = 0.0;

  	FOR(k,nGene) {
		if( LnLikehood1Tree_invgamma(&genetree[k],speciestree, a,b) == ERROR) {
      	printf("Error for calculating the log likelhood of gene trees\n");
      	return(ERROR);
    	}
   }

   for(i=0; i<speciestree->nPop;i++){
     	*lnl += alpha*log(beta)-(a[i]+alpha)*log(b[i]+beta);
		totalcoal += a[i];
		for(j=0;j<a[i];j++) *lnl += log(alpha+j);
		speciestree->nodes[sptree.popIndex[i]].theta = (b[i]+beta)/(alpha+a[i]-1); 
	}
  	if(totalcoal != (speciestree->nTaxa-1)*nGene)
  	{
		printf("the number of coalescence is %f\n",totalcoal);
		return(ERROR);
  	}
  	return(NO_ERROR);
}

int LnLikehood1Tree_invgamma(SPTree *genetree, SPTree *speciestree, double *a, double *b)
{ 
   	int i,j, ipop, k=0,nt, popnode=0;
   	double *t,y;
   	CoalTime coal,coal_t;

  	if( GetCoaltime(genetree,speciestree,&coal) == ERROR)
     	{
         printf(" bugs in GetCoaltime \n");
         return(ERROR);
     	}
  
  	i=0;
  	GetNcoal(speciestree->root, &coal, genetree, speciestree, &coal_t, &i);
  
  	for(ipop=0; ipop<genetree->nPop; ipop++)
 	{      
	  nt=coal_t.ncoal[ipop]; 
	  t=coal_t.tj[ipop];   
	  for(i=0; i<nt-1; i++)  for (j=i+1; j<nt; j++)
	  {
            if (t[j]<t[i])  
			{ 
				y=t[i]; 
				t[i]=t[j];
				t[j]=y;
			}
	  }
  	}

 	for(i=0; i<genetree->nPop; i++) 
	{
      	ipop=i; t=coal_t.tj[ipop];
    
	   	/* find the population node*/
	   	for(j=0;j<sptree.nPop;j++)
		{
			if(sptree.popIndex[j] == coal_t.nodes[i])
			{
				popnode = j;
				break;
			}
		}

         if(coal_t.nin[ipop]>coal_t.nout[ipop]) {
			 a[popnode] += (coal_t.nin[ipop]-coal_t.nout[ipop]);

         for(j=coal_t.nin[ipop],k=0; j>coal_t.nout[ipop]; j--,k++) {
            y = t[k]- (k==0? 0:t[k-1]);         
		if(genetree->haploid == 1) /*haploid*/
		  	 b[popnode] += 4*j*(j-1)*y;
		else if (genetree->haploid == -1) /*zlink*/
			b[popnode] += 1.3333333*j*(j-1)*y;
		else /*diploid*/
			 b[popnode] += j*(j-1)*y;

         }
      }
      if(coal_t.nout[ipop]>1) {
         y = (coal_t.nout[ipop]==coal_t.nin[ipop] ? speciestree->nodes[coal_t.nodes[ipop]].brlens : speciestree->nodes[coal_t.nodes[ipop]].brlens-t[k-1]);
         if(genetree->haploid == 1) /*haploid*/
		  	 b[popnode] += 4*coal_t.nout[ipop]*(coal_t.nout[ipop]-1)*y;
		 else if(genetree->haploid == -1) /*zlink*/
		  	 b[popnode] += 1.3333333*coal_t.nout[ipop]*(coal_t.nout[ipop]-1)*y;
		else /*diploid*/
			 b[popnode] += coal_t.nout[ipop]*(coal_t.nout[ipop]-1)*y;
         
      }
   }
   return(NO_ERROR);
}

int GetCoaltime(SPTree *genetree,SPTree *speciestree, CoalTime *coal)
{  
	int inode,i,k,inode_gene,inode_sp,stop,index=0;
	int *genetreenodes;
	double distance_gene;

	genetreenodes = (int*)malloc(2*genetree->nTaxa*sizeof(int));
	if(!genetreenodes)
	{
		printf("allocating problem for genetreenodes\n");
		return(ERROR);
	}
	FOR(i,genetree->nSpecies)
      { 
        	FOR(inode,genetree->speciesIndex[i][0])
        	{        
        		inode_gene=genetree->speciesIndex[i][inode+1];      
        		stop=0;        
        		do{
				inode_gene = genetree->nodes[inode_gene].father; 
				FOR(k,index)
              			if(inode_gene==genetreenodes[k]) 
					{
						stop=1;
						break;
					}
           			if(stop==1) break;

				/*find the node in the species tree that is right below the inode_gene*/
				distance_gene = genetree->nodes[inode_gene].age;
            		inode_sp = FindspNodeBelow(speciestree,i,distance_gene);
				/*record the coalescence*/
				coal->ncoal[index]=1;
                 		coal->nodes[index]=inode_sp;
                 		coal->tj[index][0]= distance_gene - speciestree->nodes[inode_sp].age;
                 		genetreenodes[index]=inode_gene;
                 		index++;					
           		} while(inode_gene!=genetree->root);
		}
    	}
	
    	/* check if all the Taxa coalesce */
    	if(index != genetree->nTaxa-1)
    	{
       	printf("Error for the coalescence time calculation\n");
       	return(ERROR);
    	}
     	free(genetreenodes);
	return(NO_ERROR);
}

int FindspNodeBelow(SPTree *speciestree, int spnode, double dis_gene)
{
	int inode=spnode,father;
	double dis_sp;

	father = speciestree->nodes[inode].father;
	dis_sp = speciestree->nodes[father].age;
	while(dis_gene >= dis_sp)
	{
		if(father == speciestree->root)
		{
			inode = speciestree->root;
			break;
		}
		inode = speciestree->nodes[inode].father;
		father = speciestree->nodes[inode].father;
		dis_sp = speciestree->nodes[father].age;
	}
	return(inode);
}


int GetNcoal(int inode, CoalTime *coal, SPTree *genetree, SPTree *speciestree, CoalTime *coal_t,int *index)
{
  	int k,i,kson,w;

  	if(inode<speciestree->nSpecies)   
  	{
     		if(genetree->speciesIndex[inode][0] == 1 )    
        		return(1);
    		else
       	{
         		coal_t->nodes[*index] = inode;
         		coal_t->nin[*index] = genetree->speciesIndex[inode][0];
         		w=0;
         		FOR(i,genetree->nTaxa-1)
         		{
           			if(coal->nodes[i]==inode) 
				{  
					coal_t->tj[*index][w]=coal->tj[i][0];
					w++;
				}
	  		}
  
	  		coal_t->ncoal[*index]=w;
         		w=*index;
	  		(*index)++;
         		coal_t->nout[w] = coal_t->nin[w] - coal_t->ncoal[w];
         		return(coal_t->nout[w]);
      	}
  	}
  	else
  	{
	 	w=0;
	 	FOR(i,genetree->nTaxa-1) 
	 	{
        		if(coal->nodes[i]==inode) 
			{  
				coal_t->tj[*index][w]=coal->tj[i][0];
				w++;
			}
	 	}
  
	 	coal_t->ncoal[*index]=w;
	 	coal_t->nodes[*index]=inode;
	 	w=*index;
	 	(*index)++;
	 	coal_t->nin[w]=0;
	 	FOR(k,speciestree->nodes[inode].nson)
	 	{
		 	kson=speciestree->nodes[inode].sons[k];
		 	coal_t->nin[w]+=GetNcoal(kson,coal,genetree,speciestree,coal_t,index);
	 	}
	 	coal_t->nout[w]=coal_t->nin[w]-coal_t->ncoal[w];
 
	 	return(coal_t->nout[w]);
  	}
}

double SPLnP1 (double t, double l, double m, double r) {
	double		p0t;
	p0t = r*(l-m) / (r*l + (l*(1.0-r)-m)*exp((m-l)*t) );
	return (log(1.0/r) + 2.0*log(p0t) + (m-l)*t);
}

double SPLnVt (double t, double l, double m, double r) {
	double		p0t;
	p0t = r*(l-m) / (r*l + (l*(1.0-r)-m)*exp((m-l)*t) );
	return (log(1.0 - (1.0/r) * p0t * exp((m-l)*t)));
}

int SPLnBirthDeathPriorPr (double *prob, double sR, double eR, double sF,SPTree *speciestree) {
	int		i, j;
	double	*nt,rootTime;
   nt = (double*)malloc((speciestree->nSpecies-1)*sizeof(double));
   if(!nt) {
   	printf("allocating problem for nt\n");
      return(ERROR);
   }
	// rescale all of the node times on the tree
	rootTime = speciestree->nodes[speciestree->root].age;
   j = 0;
   for (i=speciestree->nSpecies; i<2*speciestree->nSpecies-1; i++)
   	if(i != speciestree->root) nt[j++] = speciestree->nodes[i].age;
   nt[speciestree->nSpecies-2] = rootTime;
	for (i=0; i<speciestree->nSpecies-1; i++) {
		nt[i] /= rootTime;
      if(nt[i]>1.0) {
			printf("nt time error\n");
			return(ERROR);
		}
	}
	/* I think this is correct. It looks as if Yang and Rannala (1997)
	   have the root time constrained to be 1.0. */
	rootTime = 1.0;
							
	/* calculate probabilities of tree */
	if ((sR - eR)> 1e-8) {
		(*prob) = (speciestree->nSpecies - 2.0) * log(sR);
		for (i=0; i<speciestree->nSpecies-2; i++)
			(*prob) += SPLnP1 (nt[i], sR, eR, sF) - SPLnVt (rootTime, sR, eR, sF);
	} else {
		(*prob) = 0.0;
		for (i=0; i<speciestree->nSpecies-2; i++)
			(*prob) += log (1.0 + sF * eR) - (2.0 * log(1.0 + sF * eR * nt[i]));
	}

	/*plus the probability for the root*/
	*prob -= modelParam.treeHeightExp*speciestree->nodes[speciestree->root].age;
	free(nt);
	return (NO_ERROR);
}

int FindConstraint(int node, int *anode)
{
	int i;

	if(node == sptree.treeConstraint[0].nodes[0])
	{
		*anode = sptree.treeConstraint[0].nodes[1];
		return(0);
	}
	
	for(i=0;i<sptree.nSpecies;i++)
		if(node == sptree.treeConstraint[i].nodes[1])
			break;
	*anode = sptree.treeConstraint[i].nodes[0];
	return(i);
}

MrBFlt Prob_sptree(Distance *tau,int ntime)
{
  int i,j,round=100;
  MrBFlt time,*lnp,lamda,roottime,randomroottime,max,lnprior2=0.0;
 
  lnp = (MrBFlt *)malloc((size_t) ((round) * sizeof(MrBFlt)));
		
  if (!lnp)
			{
			printf ("%s   Problem allocating lnp \n", spacer);
			exit(-1);
			}

 
  /* find the root depth and save it into tau[ntime-1]*/
  for(i=(ntime-2);i>=0;i--)
  {
	  if(tau[ntime-1].dist < tau[i].dist)
	  {
		  time = tau[ntime-1].dist;
		  tau[ntime-1].dist = tau[i].dist;
		  tau[i].dist = time;
	  }
  }

 
  for(i=0;i<round;i++)
   {
     lnp[i] = 0.0;
     lamda = (double)(i+1)/10;
    
     /*generate a random roottime */
     roottime = tau[ntime-1].dist*0.5;
     randomroottime = (1-exp(-roottime))*rndu();
     
     /* given the randomroottime, calculate the density function */
     for(j=0;j<ntime-1;j++)
       {
         time = tau[j].dist/2;
         if(time < randomroottime) lnp[i] += log(1-exp(-lamda*time))-log((1-exp(-lamda*randomroottime)));
        
       }  
   }
   
   /*find the maximum*/
   max = lnp[0];
 
 
   for(i=1;i<round;i++) if(max<lnp[i]) max = lnp[i];
   for(i=0;i<round;i++) if(max-lnp[i]<=10) 
   {
	   lnprior2 += exp(lnp[i]-max);
   }   

   /*log prior is equal to the log of the average*/
   lnprior2 = -log(round) + max + log(lnprior2);

   /* plus the probability of the roottime*/    
   lnprior2 += log(1-exp(-roottime)); 
   
   free(lnp);
  
   return(lnprior2);
}

#define TARGETLENDELTA (100)

int SPSaveSprintf(char **target, int *targetLen, char *fmt, ...) {
  va_list    argp;
  int        len,retval;

  va_start(argp, fmt);
#ifdef VISUAL
  len = _vsnprintf(NULL, 0, fmt, argp);
#else
  len = vsnprintf(NULL, 0, fmt, argp);
#endif

  va_end(argp);

  if(len>*targetLen)
        {
/*        fprintf(stderr, "* increasing buffer from %d to %d bytes\n", *targetLen, len+TARGETLENDELTA); */
        *targetLen = len+TARGETLENDELTA; /* make it a little bigger */
            *target = (char*)realloc(*target, *targetLen);
        }

  va_start(argp,fmt);
  retval=vsprintf(*target, fmt, argp);
  va_end(argp);

/*   fprintf(stderr, "* savesprintf /%s/\n",*target); */
  return retval;
}

