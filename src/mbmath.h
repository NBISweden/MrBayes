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
MrBFlt  BetaQuantile (MrBFlt alpha, MrBFlt beta, MrBFlt x);
void    CalcCijk (int dim, MrBFlt *c_ijk, MrBFlt **u, MrBFlt **v);
void    CopyComplexMatrices (int dim, complex **from, complex **to);
void    CopyDoubleMatrices (int dim, MrBFlt **from, MrBFlt **to);
void    DirichletRandomVariable (MrBFlt *alp, MrBFlt *z, int n, SafeLong *seed);
int     DiscreteGamma (MrBFlt *rK, MrBFlt alfa, MrBFlt beta, int K, int median);
void    FreeSquareComplexMatrix (complex **m);
void    FreeSquareDoubleMatrix (MrBFlt **m);
void    FreeSquareIntegerMatrix (int **m);
int     GetEigens (int dim, MrBFlt **q, MrBFlt *eigenValues, MrBFlt *eigvalsImag, MrBFlt **eigvecs, MrBFlt **inverseEigvecs, complex **Ceigvecs, complex **CinverseEigvecs);
MrBFlt  LnFactorial (int value);
MrBFlt  LnGamma (MrBFlt alp);
MrBFlt  LnPriorProbExponential(MrBFlt val, MrBFlt *params);
MrBFlt  LnPriorProbFix(MrBFlt val, MrBFlt *params);
MrBFlt  LnPriorProbGamma(MrBFlt val, MrBFlt *params);
MrBFlt  LnPriorProbLognormal(MrBFlt val, MrBFlt *params);
MrBFlt  LnPriorProbNormal(MrBFlt val, MrBFlt *params);
MrBFlt  LnPriorProbTruncatedNormal(MrBFlt val, MrBFlt *params);
MrBFlt  LnPriorProbUniform(MrBFlt val, MrBFlt *params);
MrBFlt  LnProbRatioExponential (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt  LnProbRatioGamma (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt  LnProbRatioLognormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt  LnProbRatioNormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt  LnProbRatioTruncatedNormal (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt  LnProbRatioUniform (MrBFlt newX, MrBFlt oldX, MrBFlt *params);
MrBFlt  LnProbGamma (MrBFlt alpha, MrBFlt beta, MrBFlt x);
MrBFlt  LnProbTruncGamma (MrBFlt alpha, MrBFlt beta, MrBFlt x, MrBFlt min, MrBFlt max);
MrBFlt  LnProbLogNormal (MrBFlt exp, MrBFlt sd, MrBFlt x);
MrBFlt  LnProbTK02LogNormal (MrBFlt mean, MrBFlt var, MrBFlt x);
MrBFlt  LnProbScaledGamma (MrBFlt alpha, MrBFlt x);
MrBFlt  LnRatioTK02LogNormal (MrBFlt exp, MrBFlt sd, MrBFlt xNew, MrBFlt xOld);
MrBFlt  LnRatioLogNormal (MrBFlt exp, MrBFlt sd, MrBFlt xNew, MrBFlt xOld);
MrBFlt  LogNormalRandomVariable (MrBFlt mean, MrBFlt var, SafeLong *seed);
void    MultiplyMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result);
int     MultiplyMatrixNTimes (int dim, MrBFlt **Mat, int power, MrBFlt **Result);
MrBFlt  PointNormal (MrBFlt prob);
MrBFlt  PsiGammaLnProb (MrBFlt alpha, MrBFlt value);
MrBFlt  PsiGammaLnRatio (MrBFlt alpha, MrBFlt numerator, MrBFlt denominator);
MrBFlt  PsiGammaRandomVariable (MrBFlt alpha, SafeLong *seed);
MrBFlt  QuantileGamma (MrBFlt x, MrBFlt alfa, MrBFlt beta);
MrBFlt  RandomNumber (SafeLong *seed);
