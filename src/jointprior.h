int		LnJointGenetreePr(Tree *t[],int *updatedtreeid, int num_tree, MrBFlt *lnprior, MrBFlt *GeneMu, SPTree *speciestree);
int		ReadControlfile(FILE *fdata);
void 	InitiateParam(void);
int		SPPrintTreeTitle (int curGen, FILE *fout);
int		SPLnBirthDeathPriorPr (double *prob, double sR, double eR, double sF,SPTree *speciestree);
int		SPPrintTree(int curGen, SPTree *tree, int showBrlens, int showTheta, int showMu, int isRooted);
/*species tree mutation rate*/
int		populationMutation (Tree *genetree,SPTree *speciestree, MrBFlt genemu);


