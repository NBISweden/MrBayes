#undef complex
struct complex
	{
	MrBFlt re;
	MrBFlt im;
	};

typedef struct complex complex;


complex **AllocateSquareComplexMatrix (int dim);
MrBFlt  **AllocateSquareDoubleMatrix (int dim);
int     **AllocateSquareIntegerMatrix (int dim);
int     AutodGamma (MrBFlt *M, MrBFlt rho, int K);
void    BetaBreaks (MrBFlt alpha, MrBFlt beta, MrBFlt *values, int K);
void    CalcCijk (int dim, MrBFlt *c_ijk, MrBFlt **u, MrBFlt **v);
void    CopyComplexMatrices (int dim, complex **from, complex **to);
void    CopyDoubleMatrices (int dim, MrBFlt **from, MrBFlt **to);
void    DirichletRandomVariable (MrBFlt *alp, MrBFlt *z, int n, safeLong *seed);
int     DiscreteGamma (MrBFlt *rK, MrBFlt alfa, MrBFlt beta, int K, int median);
void    FreeSquareComplexMatrix (complex **m);
void    FreeSquareDoubleMatrix (MrBFlt **m);
void    FreeSquareIntegerMatrix (int **m);
int     GetEigens (int dim, MrBFlt **q, MrBFlt *eigenValues, MrBFlt *eigvalsImag, MrBFlt **eigvecs, MrBFlt **inverseEigvecs, complex **Ceigvecs, complex **CinverseEigvecs);
MrBFlt  LnFactorial (int value);
MrBFlt  LnGamma (MrBFlt alp);
MrBFlt  LnProbGamma (MrBFlt alpha, MrBFlt beta, MrBFlt x);
MrBFlt  LnProbLogNormal (MrBFlt exp, MrBFlt var, MrBFlt x);
MrBFlt  LnProbScaledGamma (MrBFlt alpha, MrBFlt x);
MrBFlt  LnRatioLogNormal (MrBFlt exp, MrBFlt var, MrBFlt xNew, MrBFlt xOld);
MrBFlt  LogNormalRandomVariable (MrBFlt mean, MrBFlt var, safeLong *seed);
void    MultiplyMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result);
int     MultiplyMatrixNTimes (int dim, MrBFlt **Mat, int power, MrBFlt **Result);
MrBFlt  PointNormal (MrBFlt prob);
MrBFlt  PsiGammaLnProb (MrBFlt alpha, MrBFlt value);
MrBFlt  PsiGammaLnRatio (MrBFlt alpha, MrBFlt numerator, MrBFlt denominator);
MrBFlt  PsiGammaRandomVariable (MrBFlt alpha, safeLong *seed);
MrBFlt  QuantileGamma (MrBFlt x, MrBFlt alfa, MrBFlt beta);
MrBFlt  RandomNumber (safeLong *seed);
