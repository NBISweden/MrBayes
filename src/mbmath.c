/*
 *  MrBayes 3
 *
 *  (c) 2002-2010
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
 *  Maxim Teslenko (maxim.teslenko@nrm.se)
 *
 *  and by many users (run 'acknowledgements' to see more info)
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdarg.h>
#include "mb.h"
#include "globals.h"
#include "mbmath.h"
#include "bayes.h"
#include "model.h"
#include "utils.h"

#define	MAX_GAMMA_CATS						20
#define PI                                  3.14159265358979324
#define PIOVER2 							1.57079632679489662
#define POINTGAMMA(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define PAI2								6.283185307
#define TINY								1.0e-20
#define EVALUATE_COMPLEX_NUMBERS 			2
#if !defined(MAX)
#define MAX(a,b)							(((a) > (b)) ? (a) : (b))
#endif
#if !defined(MIN)
#define MIN(a,b)							(((a) < (b)) ? (a) : (b))
#endif
#define SQUARE(a)							((a)*(a))

/* prototypes */
void    AddTwoMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result);
void    BackSubstitutionRow (int dim, MrBFlt **u, MrBFlt *b);
void    Balanc (int dim, MrBFlt **a, int *low, int *high, MrBFlt *scale);
void    BalBak (int dim, int low, int high, MrBFlt *scale, int m, MrBFlt **z);
MrBFlt  BetaCf (MrBFlt a, MrBFlt b, MrBFlt x);
MrBFlt  BetaQuantile (MrBFlt alpha, MrBFlt beta, MrBFlt x);
MrBFlt  CdfBinormal (MrBFlt h1, MrBFlt h2, MrBFlt r);
MrBFlt  CdfNormal (MrBFlt x);
complex Complex (MrBFlt a, MrBFlt b);
MrBFlt  ComplexAbsoluteValue (complex a);
complex ComplexAddition (complex a, complex b);
complex ComplexConjugate (complex a);
complex ComplexDivision (complex a, complex b);
void    ComplexDivision2 (MrBFlt ar, MrBFlt ai, MrBFlt br, MrBFlt bi, MrBFlt *cr, MrBFlt *ci);
complex ComplexExponentiation (complex a);
int     ComplexInvertMatrix (int dim, complex **a, MrBFlt *dwork, int *indx, complex **aInverse, complex *col);
complex ComplexLog (complex a);
void    ComplexLUBackSubstitution (int dim, complex **a, int *indx, complex *b);
int     ComplexLUDecompose (int dim, complex **a, MrBFlt *vv, int *indx, MrBFlt *pd);
complex ComplexMultiplication (complex a, complex b);
complex ComplexSquareRoot (complex a);
complex ComplexSubtraction (complex a, complex b);
int     ComputeEigenSystem (int dim, MrBFlt **a, MrBFlt *v, MrBFlt *vi, MrBFlt **u, int *iwork, MrBFlt *dwork);
void    ComputeLandU (int dim, MrBFlt **aMat, MrBFlt **lMat, MrBFlt **uMat);
void    ComputeMatrixExponential (int dim, MrBFlt **a, int qValue, MrBFlt **f);
void    DivideByTwos (int dim, MrBFlt **a, int power);
MrBFlt  D_sign (MrBFlt a, MrBFlt b);
int     EigensForRealMatrix (int dim, MrBFlt **a, MrBFlt *wr, MrBFlt *wi, MrBFlt **z, int *iv1, MrBFlt *fv1);
void    ElmHes (int dim, int low, int high, MrBFlt **a, int *interchanged);
void    ElTran (int dim, int low, int high, MrBFlt **a, int *interchanged, MrBFlt **z);
void    Exchange (int j, int k, int l, int m, int n, MrBFlt **a, MrBFlt *scale);
MrBFlt  Factorial (int x);
void    ForwardSubstitutionRow (int dim, MrBFlt **L, MrBFlt *b);
MrBFlt  GammaRandomVariable (MrBFlt a, MrBFlt b, safeLong *seed);
void    GaussianElimination (int dim, MrBFlt **a, MrBFlt **bMat, MrBFlt **xMat);
int     Hqr2 (int dim, int low, int high, MrBFlt **h, MrBFlt *wr, MrBFlt *wi, MrBFlt **z);
MrBFlt  IncompleteBetaFunction (MrBFlt alpha, MrBFlt beta, MrBFlt x);
MrBFlt  IncompleteGamma (MrBFlt x, MrBFlt alpha, MrBFlt LnGamma_alpha);
int     InvertMatrix (int dim, MrBFlt **a, MrBFlt *col, int *indx, MrBFlt **aInv);
MrBFlt  LBinormal (MrBFlt h1, MrBFlt h2, MrBFlt r);
int     LogBase2Plus1 (MrBFlt x);
void    LUBackSubstitution (int dim, MrBFlt **a, int *indx, MrBFlt *b);
int     LUDecompose (int dim, MrBFlt **a, MrBFlt *vv, int *indx, MrBFlt *pd);
void    MultiplyMatrixByScalar (int dim, MrBFlt **a, MrBFlt scalar, MrBFlt **result);
MrBFlt  PointChi2 (MrBFlt prob, MrBFlt v);
MrBFlt  PointNormal (MrBFlt prob);
void    PrintComplexVector (int dim, complex *vec);
void    PrintSquareComplexMatrix (int dim, complex **m);
void    PrintSquareDoubleMatrix (int dim, MrBFlt **matrix);
void    PrintSquareIntegerMatrix (int dim, int **matrix);
complex ProductOfRealAndComplex (MrBFlt a, complex b);
MrBFlt  RndGamma (MrBFlt s, safeLong *seed);
MrBFlt  RndGamma1 (MrBFlt s, safeLong *seed);
MrBFlt  RndGamma2 (MrBFlt s, safeLong *seed);
int     SetQvalue (MrBFlt tol);
void    SetToIdentity (int dim, MrBFlt **matrix);
MrBFlt  Tha (MrBFlt h1, MrBFlt h2, MrBFlt a1, MrBFlt a2);
void    TiProbsUsingEigens (int dim, MrBFlt *cijk, MrBFlt *eigenVals, MrBFlt v, MrBFlt r, MrBFlt **tMat, MrBFlt **fMat, MrBFlt **sMat);
void    TiProbsUsingPadeApprox (int dim, MrBFlt **qMat, MrBFlt v, MrBFlt r, MrBFlt **tMat, MrBFlt **fMat, MrBFlt **sMat);





/*---------------------------------------------------------------------------------
|
|   AddTwoMatrices
|
|   Takes the sum of two matrices, "a" and "b", and puts the results in a matrix
|   called "result".
|
---------------------------------------------------------------------------------*/
void AddTwoMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result)

{

	int			row, col;

	for (row=0; row<dim; row++)
		{
		for (col=0; col<dim; col++) 
			{
			result[row][col] = a[row][col] + b[row][col];
			}
		}

}





/*---------------------------------------------------------------------------------
|
|   AllocateSquareComplexMatrix
|
|   Allocate memory for a square (dim X dim) complex matrix.
|
---------------------------------------------------------------------------------*/
complex **AllocateSquareComplexMatrix (int dim)

{

	int 		i;
	complex 	**m;

	m = (complex **) SafeMalloc((size_t)((dim)*sizeof(complex*)));
	if (!m) 
		{
		MrBayesPrint ("%s   Error: Problem allocating a square complex matrix.\n", spacer);
		exit (0);
		}
	m[0]=(complex *) SafeMalloc((size_t)((dim*dim)*sizeof(complex)));
	if (!m[0]) 
		{
		MrBayesPrint ("%s   Error: Problem allocating a square complex matrix.\n", spacer);
		exit (0);
		}
	for(i=1;i<dim;i++) 
		{
		m[i] = m[i-1] + dim;
		}
		
	return (m);
	
}






/*---------------------------------------------------------------------------------
|
|   AllocateSquareDoubleMatrix
|
|   Allocate memory for a square (dim X dim) matrix of doubles.
|
---------------------------------------------------------------------------------*/
MrBFlt **AllocateSquareDoubleMatrix (int dim)

{

	int			i;
	MrBFlt		**m;
	
	m = (MrBFlt **)SafeMalloc((size_t)((dim)*sizeof(MrBFlt*)));
	if (!m)
		{
		MrBayesPrint ("%s   Error: Problem allocating a square matrix of doubles.\n", spacer);
		exit(1);
		}
	m[0] = (MrBFlt *)SafeMalloc((size_t)((dim*dim)*sizeof(MrBFlt)));
	if (!m[0])
		{
		MrBayesPrint ("%s   Error: Problem allocating a square matrix of doubles.\n", spacer);
		exit(1);
		}
	for(i=1; i<dim; i++)
		{
		m[i] = m[i-1] + dim;
		}

	return (m);
	
}





/*---------------------------------------------------------------------------------
|
|   AllocateSquareIntegerMatrix
|
|   Allocate memory for a square (dim X dim) matrix of integers.
|
---------------------------------------------------------------------------------*/
int **AllocateSquareIntegerMatrix (int dim)

{

	int		i, **m;
	
	m = (int **)SafeMalloc((size_t)((dim)*sizeof(int*)));
	if (!m)
		{
		MrBayesPrint ("%s   Error: Problem allocating a square matrix of integers.\n", spacer);
		exit(1);
		}
	m[0] = (int *)SafeMalloc((size_t)((dim*dim)*sizeof(int)));
	if (!m[0])
		{
		MrBayesPrint ("%s   Error: Problem allocating a square matrix of integers.\n", spacer);
		exit(1);
		}
	for(i=1; i<dim; i++)
		{
		m[i] = m[i-1] + dim;
		}

	return (m);
	
}





/*---------------------------------------------------------------------------------
|
|   AutodGamma
|
|   Auto-discrete-gamma distribution of rates over sites, K equal-probable
|   categories, with the mean for each category used.         
|   This routine calculates M[], using rho and K (numGammaCats)     
|
---------------------------------------------------------------------------------*/
int AutodGamma (MrBFlt *M, MrBFlt rho, int K)

{

	int			i, j, i1, i2;
	MrBFlt		point[MAX_GAMMA_CATS], x, y, large = 20.0, sum;

	for (i=0; i<K-1; i++) 
		point[i] = PointNormal ((i + 1.0) / K);
	for (i=0; i<K; i++) 
		{
		for (j=0; j<K; j++) 
			{
			x = (i < K-1 ? point[i]:large);
			y = (j < K-1 ? point[j]:large);
        	M[i * K + j] = CdfBinormal (x, y, rho);
			}
		}
	for (i1=0; i1<2*K-1; i1++) 
		{
		for (i2=0; i2<K*K; i2++) 
			{
			i = i2 / K; 
			j = i2 % K;
			if (AreDoublesEqual(i+j, 2*(K-1.0)-i1, ETA)==NO)
				continue;
			y = 0;
			if (i > 0) 
				y -= M[(i-1)*K+j];
			if (j > 0) 
				y -= M[i*K+(j-1)];
			if (i > 0 && j > 0) 
				y += M[(i-1)*K+(j-1)];
     	    M[i*K+j] = (M[i*K+j] + y) * K;
			}
		}
	for (i=0; i<K; i++)
		{
		sum = 0.0;
		for (j=0; j<K; j++)
			{
			if (M[i*K+j] < 0.0)
				M[i*K+j] = 0.0;
			sum += M[i*K+j];
			}
		for (j=0; j<K; j++)
			M[i*K+j] /= sum;
		}
		
#	if 0
	MrBayesPrint ("rho = %lf\n", rho);
	for (i=0; i<K; i++)
		{
		for (j=0; j<K; j++)
			MrBayesPrint ("%lf ", M[i*K + j]);
		MrBayesPrint ("\n");
		}
#	endif
	
	return (NO_ERROR);
	
}





/*---------------------------------------------------------------------------------
|
|   BackSubstitutionRow
|
---------------------------------------------------------------------------------*/
void BackSubstitutionRow (int dim, MrBFlt **u, MrBFlt *b)

{

	int				i, j;
	MrBFlt			dotProduct;

	b[dim-1] /= u[dim-1][dim-1];
	for (i=dim-2; i>=0; i--) 
		{
		dotProduct = 0.0;
		for (j=i+1; j<dim; j++)
			dotProduct += u[i][j] * b[j];
		b[i] = (b[i] - dotProduct) / u[i][i];
		}

}





/*---------------------------------------------------------------------------------
|
|   Balanc
|
|   This subroutine balances a real matrix and isolates
|   eigenvalues whenever possible.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * a contains the input matrix to be balanced
|
|   On output:
|
|    * a contains the balanced matrix.
|
|    * low and high are two integers such that a(i,j)
|      is equal to zero if
|         (1) i is greater than j and
|         (2) j=1,...,low-1 or i=igh+1,...,n.
|
|    * scale contains information determining the
|      permutations and scaling factors used.
|
|   Suppose that the principal submatrix in rows pLow through pHigh
|   has been balanced, that p(j) denotes the index interchanged
|   with j during the permutation step, and that the elements
|   of the diagonal matrix used are denoted by d(i,j). Then
|      scale(j) = p(j),    for j = 1,...,pLow-1
|               = d(j,j),      j = pLow,...,pHigh
|               = p(j)         j = pHigh+1,...,dim.
|   The order in which the interchanges are made is dim to pHigh+1,
|   then 1 to pLow-1.
|
|   Note that 1 is returned for pHigh if pHigh is zero formally.
|
|   The algol procedure exc contained in balance appears in
|   balanc in line.  (Note that the algol roles of identifiers
|   k,l have been reversed.)
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|
|   This function was converted from FORTRAN by D. L. Swofford.
|   
---------------------------------------------------------------------------------*/
void Balanc (int dim, MrBFlt **a, int *low, int *high, MrBFlt *scale)

{

	int			i, j, k, l, m, noconv;
	MrBFlt		c, f, g, r, s, b2;

	b2 = FLT_RADIX * FLT_RADIX;
	k = 0;
	l = dim - 1;
   	
	for (j=l; j>=0; j--)
		{
		for (i=0; i<=l; i++)
			{
			if (i != j)
				{
				  if (AreDoublesEqual(a[j][i],0.0, ETA)==NO)
					goto next_j1;
				}
			}
			
		/* bug that DLS caught */
		/*m = l;
		Exchange(j, k, l, m, dim, a, scale);
		if (l < 0)
			goto leave;
		else
			j = --l;*/
		m = l;
		Exchange(j, k, l, m, dim, a, scale);
		if (--l < 0)
			goto leave;
		next_j1:
			;
		}

	for (j=k; j<=l; j++)
		{
		for (i=k; i<=l; i++)
			{
			if (i != j)
				{
				if (AreDoublesEqual(a[i][j], 0.0, ETA)==NO)
					goto next_j;
				}
			}
		m = k;
		Exchange(j, k, l, m, dim, a, scale);
		k++;
		next_j:
			;
		}

	for (i=k; i<=l; i++)
		scale[i] = 1.0;
	
	do	{
		noconv = FALSE;
		for (i=k; i<=l; i++)
			{
			c = 0.0;
			r = 0.0;
			for (j=k; j<=l; j++)
				{
				if (j != i)
					{
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
					}
				}
			if (AreDoublesEqual(c,0.0,ETA)==NO && AreDoublesEqual(r,0.0,ETA)==NO)
				{
				g = r / FLT_RADIX;
				f = 1.0;
				s = c + r;
				while (c < g)
					{
					f *= FLT_RADIX;
					c *= b2;
					}
				g = r * FLT_RADIX;
				while (c >= g)
					{
					f /= FLT_RADIX;
					c /= b2;
					}
				if ((c + r) / f < s * .95)
					{
					g = 1.0 / f;
					scale[i] *= f;
					noconv = TRUE;				
					for (j=k; j<dim; j++)
						a[i][j] *= g;
					for (j=0; j<=l; j++)
						a[j][i] *= f;
					}
				}
			}	
		}
		while (noconv);
	leave:
		*low = k;
		*high = l;
	
#	if 0 
/* begin f2c version of code:
   balanc.f -- translated by f2c (version 19971204) */
int balanc (int *nm, int *n, MrBFlt *a, int *low, int *igh, MrBFlt *scale)

{

	/* System generated locals */
	int a_dim1, a_offset, i__1, i__2;
	MrBFlt d__1;

	/* Local variables */
	static MrBFlt iexc;
	static MrBFlt c__, f, g;
	static MrBFlt i__, j, k, l, m;
	static MrBFlt r__, s, radix, b2;
	static MrBFlt jj;
	static logical noconv;

	/* parameter adjustments */
	--scale;
	a_dim1 = *nm;
	a_offset = a_dim1 + 1;
	a -= a_offset;

	/* function Body */
	radix = 16.0;

	b2 = radix * radix;
	k = 1;
	l = *n;
	goto L100;
	
	/* .......... in-line procedure for row and column exchange .......... */
	L20:
	scale[m] = (MrBFlt) j;
	if (j == m) 
		goto L50;

	i__1 = l;
	for (i__ = 1; i__ <= i__1; ++i__) 
		{
		f = a[i__ + j * a_dim1];
		a[i__ + j * a_dim1] = a[i__ + m * a_dim1];
		a[i__ + m * a_dim1] = f;
		/* L30: */
		}

	i__1 = *n;
	for (i__ = k; i__ <= i__1; ++i__) 
		{
		f = a[j + i__ * a_dim1];
		a[j + i__ * a_dim1] = a[m + i__ * a_dim1];
		a[m + i__ * a_dim1] = f;
		/* L40: */
		}

	L50:
	switch (iexc) 
		{
		case 1:  
			goto L80;
		case 2:  
			goto L130;
		}
		
	/* .......... search for rows isolating an eigenvalue and push them down .......... */
	L80:
	if (l == 1) 
		goto L280;
	--l;
	
	/* .......... for j=l step -1 until 1 do -- .......... */
	L100:
	i__1 = l;
	for (jj = 1; jj <= i__1; ++jj) 
		{
		j = l + 1 - jj;
		i__2 = l;
		for (i__ = 1; i__ <= i__2; ++i__) 
			{
			if (i__ == j) 
				goto L110;
			if (a[j + i__ * a_dim1] != 0.) 
				goto L120;
			L110:
			;
			}
		m = l;
		iexc = 1;
		goto L20;
		L120:
		;
		}

	goto L140;
	/* .......... search for columns isolating an eigenvalue and push them left .......... */
	L130:
	++k;

	L140:
	i__1 = l;
	for (j = k; j <= i__1; ++j) 
		{
		i__2 = l;
		for (i__ = k; i__ <= i__2; ++i__) 
			{
			if (i__ == j) 
				goto L150;
			if (a[i__ + j * a_dim1] != 0.) 
				goto L170;
			L150:
			;
			}
		m = k;
		iexc = 2;
		goto L20;
		L170:
		;
		}
		
	/* .......... now balance the submatrix in rows k to l .......... */
	i__1 = l;
	for (i__ = k; i__ <= i__1; ++i__) 
		{
		/* L180: */
		scale[i__] = 1.0;
		}
	/* .......... iterative loop for norm reduction .......... */
	L190:
	noconv = FALSE;

	i__1 = l;
	for (i__ = k; i__ <= i__1; ++i__) 
		{
		c__ = 0.0;
		r__ = 0.0;
		i__2 = l;
		for (j = k; j <= i__2; ++j) 
			{
			if (j == i__) 
				goto L200;
			c__ += (d__1 = a[j + i__ * a_dim1], abs(d__1));
			r__ += (d__1 = a[i__ + j * a_dim1], abs(d__1));
			L200:
			;
			}
		
		/* .......... guard against zero c or r due to underflow .......... */
		if (c__ == 0. || r__ == 0.) 
			goto L270;
		g = r__ / radix;
		f = 1.0;
		s = c__ + r__;
		L210:
		if (c__ >= g) 
			goto L220;
		f *= radix;
		c__ *= b2;
		goto L210;
		L220:
		g = r__ * radix;
		L230:
		if (c__ < g) 
			goto L240;
		f /= radix;
		c__ /= b2;
		goto L230;
		
		/*     .......... now balance .......... */
		L240:
		if ((c__ + r__) / f >= s * .95) 
			goto L270;
		g = 1.0 / f;
		scale[i__] *= f;
		noconv = TRUE;
		
		i__2 = *n;
		for (j = k; j <= i__2; ++j) 
			{
			/* L250: */
			a[i__ + j * a_dim1] *= g;
			}

		i__2 = l;
		for (j = 1; j <= i__2; ++j) 
			{
			/* L260: */
			a[j + i__ * a_dim1] *= f;
			}

		L270:
		;
		}

	if (noconv) 
		goto L190;

	L280:
	*low = k;
	*igh = l;
	return 0;
	
} 
/* end f2c version of code */
#	endif
	
}





/*---------------------------------------------------------------------------------
|
|   BalBak
|
|   This subroutine forms the eigenvectors of a real general 
|   matrix by back transforming those of the corresponding 
|   balanced matrix determined by  balance. 
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * low and high are integers determined by  balance
|
|    * scale contains information determining the permutations 
|      and scaling factors used by balance
|
|    * m is the number of columns of z to be back transformed
|
|    * z contains the real and imaginary parts of the eigen-
|      vectors to be back transformed in its first m columns
|
|   On output:
|
|    * z contains the real and imaginary parts of the
|      transformed eigenvectors in its first m columns
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|   
---------------------------------------------------------------------------------*/
void BalBak (int dim, int low, int high, MrBFlt *scale, int m, MrBFlt **z)

{

	int			i, j, k, ii;
	MrBFlt		s;

	if (m != 0) /* change "==" to "!=" to eliminate a goto statement */
		{
		if (high != low) /* change "==" to "!=" to eliminate a goto statement */
			{
			for (i=low; i<=high; i++)
				{
				s = scale[i];
				for (j=0; j<m; j++)
					z[i][j] *= s;
				}
			}
		for (ii=0; ii<dim; ii++)
			{
			i = ii;
			if ((i < low) || (i > high)) /* was (i >= lo) && (i<= hi) but this */
				{                        /* eliminates a goto statement        */
				if (i < low)
					i = low - ii;
				k = (int)scale[i];
				if (k != i) /* change "==" to "!=" to eliminate a goto statement */
					{
					for (j = 0; j < m; j++)
						{
						s = z[i][j];
						z[i][j] = z[k][j];
						z[k][j] = s;
						}
					}
				}
			}
		}

#if 0
/* begin f2c version of code:
   balbak.f -- translated by f2c (version 19971204) */
int balbak (int *nm, int *n, int *low, int *igh, MrBFlt *scale, int *m, MrBFlt *z__)

{

	/* system generated locals */
	int z_dim1, z_offset, i__1, i__2;

	/* Local variables */
	static int i__, j, k;
	static MrBFlt s;
	static int ii;

	/* parameter adjustments */
	--scale;
	z_dim1 = *nm;
	z_offset = z_dim1 + 1;
	z__ -= z_offset;

	/* function Body */
	if (*m == 0) 
		goto L200;
	if (*igh == *low) 
		goto L120;

	i__1 = *igh;
	for (i__ = *low; i__ <= i__1; ++i__) 
		{
		s = scale[i__];
		/* .......... left hand eigenvectors are back transformed */
		/*            if the foregoing statement is replaced by */
		/*            s=1.0d0/scale(i) ........... */
		i__2 = *m;
		for (j = 1; j <= i__2; ++j) 
			{
			/* L100: */
			z__[i__ + j * z_dim1] *= s;
			}

		/* L110: */
		}
		
	/* .........for i=low-1 step -1 until 1, igh+1 step 1 until n do -- .......... */
	L120:
	i__1 = *n;
	for (ii = 1; ii <= i__1; ++ii) 
		{
		i__ = ii;
		if (i__ >= *low && i__ <= *igh) 
			goto L140;
	if (i__ < *low) 
		i__ = *low - ii;
	k = (integer) scale[i__];
	if (k == i__) 
		goto L140;

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) 
		{
		s = z__[i__ + j * z_dim1];
		z__[i__ + j * z_dim1] = z__[k + j * z_dim1];
		z__[k + j * z_dim1] = s;
		/* L130: */
		}
	L140:
	;
	}

	L200:
	return 0;
	
}
/* end f2c version of code */
#endif
		
}





void BetaBreaks (MrBFlt alpha, MrBFlt beta, MrBFlt *values, int K)

{

	int				i;
	MrBFlt			r, quantile, lower, upper;
			
	r = (1.0 / K) * 0.5;
	lower = 0.0;
	upper = (1.0 / K);
	r = (upper - lower) * 0.5 + lower;
	for (i=0; i<K; i++)
		{
		quantile = BetaQuantile (alpha, beta, r);
		values[i] = quantile;
		lower += (1.0/K);
		upper += (1.0/K);
		r += (1.0/K);
		}
		
#	if 0
	for (i=0; i<K; i++)
		{
		MrBayesPrint ("%4d %lf %lf\n", i, values[i]);
		}
#	endif

}





MrBFlt BetaCf (MrBFlt a, MrBFlt b, MrBFlt x)

{

	int			m, m2;
	MrBFlt		aa, c, d, del, h, qab, qam, qap;
	
	qab = a + b;
	qap = a + 1.0;
	qam = a - 1.0;
	c = 1.0;
	d = 1.0 - qab * x / qap;
	if (fabs(d) < (1.0e-30))
		d = (1.0e-30);
	d = 1.0 / d;
	h = d;
	for (m=1; m<=100; m++)
		{
		m2 = 2 * m;
		aa = m * (b-m) * x / ((qam+m2) * (a+m2));
		d = 1.0 + aa * d;
		if (fabs(d) < (1.0e-30))
			d = (1.0e-30);
		c = 1.0 + aa / c;
		if (fabs(c) < (1.0e-30))
			c = (1.0e-30);
		d = 1.0 / d;
		h *= d * c;
		aa = -(a+m) * (qab+m) * x / ((a+m2) * (qap+m2));
		d = 1.0 + aa * d;
		if (fabs(d) < (1.0e-30))
			d = (1.0e-30);
		c = 1.0 + aa / c;
		if (fabs(c) < (1.0e-30))
			c = (1.0e-30);
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.0) < (3.0e-7))
			break;
		}
	if (m > 100)
		{
		MrBayesPrint ("%s   Error in BetaCf.\n", spacer);
		exit(0);
		}
	return (h);
	
}





MrBFlt BetaQuantile (MrBFlt alpha, MrBFlt beta, MrBFlt x)

{

	int		i, stopIter, direction, nswitches;
	MrBFlt	curPos, curFraction, increment;
	
	i = nswitches = 0;
	curPos = 0.5;
	stopIter = NO;
	increment = 0.25;
	curFraction = IncompleteBetaFunction (alpha, beta, curPos);
	if (curFraction > x)
		direction = DOWN;
	else
		direction = UP;

	while (stopIter == NO)
		{
		curFraction = IncompleteBetaFunction (alpha, beta, curPos);
		if (curFraction > x && direction == DOWN)
			{
			/* continue going down */
			while (curPos - increment <= 0.0)
				{
				increment /= 2.0;
				}
			curPos -= increment;
			}
		else if (curFraction > x && direction == UP)
			{
			/* switch directions, and go down */
			nswitches++;
			direction = DOWN;
			while (curPos - increment <= 0.0)
				{
				increment /= 2.0;
				}
			increment /= 2.0;
			curPos -= increment;
			}
		else if (curFraction < x && direction == UP)
			{
			/* continue going up */
			while (curPos + increment >= 1.0)
				{
				increment /= 2.0;
				}
			curPos += increment;
			}
		else if (curFraction < x && direction == DOWN)
			{
			/* switch directions, and go up */
			nswitches++;
			direction = UP;
			while (curPos + increment >= 1.0)
				{
				increment /= 2.0;
				}
			increment /= 2.0;
			curPos += increment;
			}
		else
			{
			stopIter = YES;
			}
		if (i > 1000 || nswitches > 20)
			stopIter = YES;
		i++;
		}
		
	return (curPos);
		
}





/*---------------------------------------------------------------------------------
|
|   CalcCijk
|
|   This function precalculates the product of the eigenvectors and their
|   inverse for faster calculation of transition probabilities. The output
|   is a vector of precalculated values. The input is the eigenvectors (u) and
|   the inverse of the eigenvector matrix (v).
|
---------------------------------------------------------------------------------*/
void CalcCijk (int dim, MrBFlt *c_ijk, MrBFlt **u, MrBFlt **v)

{

	register int 	i, j, k;
	MrBFlt 			*pc;

	pc = c_ijk;
	for (i=0; i<dim; i++)
		for (j=0; j<dim; j++)
			for (k=0; k<dim; k++)
			 	*pc++ = u[i][k] * v[k][j];

}





/*---------------------------------------------------------------------------------
|
|   CdfBinormal
|
|   F(h1,h2,r) = prob(x<h1, y<h2), where x and y are standard binormal.
|
---------------------------------------------------------------------------------*/
MrBFlt CdfBinormal (MrBFlt h1, MrBFlt h2, MrBFlt r)

{

	return (LBinormal(h1, h2, r) + CdfNormal(h1) + CdfNormal(h2) - 1.0);
   
}    





/*---------------------------------------------------------------------------------
|
|   CdfNormal
|
|   Calculates the cumulative density distribution (CDF) for the normal using:
|
|   Hill, I. D.  1973.  The normal integral.  Applied Statistics, 22:424-427.
|      (AS66)                                                  
|
---------------------------------------------------------------------------------*/
MrBFlt CdfNormal (MrBFlt x)

{

	int 			invers = 0;
	MrBFlt 			p, limit = 10.0, t = 1.28, y = x*x/2.0;

	if (x < 0.0) 
		{  
		invers = 1;  
		x  *= -1.0; 
		}
	if (x > limit)  
		return (invers ? 0 : 1);
	if (x < t)  
		p = 0.5 - x * (0.398942280444 - 0.399903438504 * y /
			(y + 5.75885480458 - 29.8213557808 /
			(y + 2.62433121679 + 48.6959930692 /
			(y + 5.92885724438))));
	else 
		p = 0.398942280385 * exp(-y) /
			(x - 3.8052e-8 + 1.00000615302 /
			(x + 3.98064794e-4 + 1.98615381364 /
			(x - 0.151679116635 + 5.29330324926 /
			(x + 4.8385912808 - 15.1508972451 /
			(x + 0.742380924027 + 30.789933034 /
			(x + 3.99019417011))))));
			
	return (invers ? p : 1-p);

}





/*---------------------------------------------------------------------------------
|
|   Complex
|
|   Returns a complex number with specified real and imaginary parts.
|
---------------------------------------------------------------------------------*/
complex Complex (MrBFlt a, MrBFlt b)

{

    complex c;
    
    c.re = a;
    c.im = b;
    
    return (c);
    
}





/*---------------------------------------------------------------------------------
|
|   ComplexAbsoluteValue
|
|   Returns the complex absolute value (modulus) of a complex number.
|
---------------------------------------------------------------------------------*/
MrBFlt ComplexAbsoluteValue (complex a)

{

	MrBFlt 		x, y, answer, temp;
	
	x = fabs(a.re);
	y = fabs(a.im);
	if(AreDoublesEqual(x, 0.0, ETA)==YES)  /* x == 0.0 */
		answer = y;
	else if (AreDoublesEqual(y, 0.0, ETA)==YES) /* y == 0.0 */
		answer = x;
	else if (x > y) 
		{
		temp = y / x;
		answer = x * sqrt(1.0 + temp * temp);
		}
	else
		{
		temp = x / y;
		answer = y * sqrt(1.0 + temp * temp);
		}

	return (answer);
	
}





/*---------------------------------------------------------------------------------
|
|   ComplexAddition
|
|   Returns the complex sum of two complex numbers.
|
---------------------------------------------------------------------------------*/
complex ComplexAddition (complex a, complex b)

{

    complex 	c;
    
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    
    return (c);
    
}





/*---------------------------------------------------------------------------------
|
|   ComplexConjugate
|
|   Returns the complex conjugate of a complex number.
|
---------------------------------------------------------------------------------*/
complex ComplexConjugate (complex a)

{

    complex 	c;
    
    c.re = a.re;
    c.im = -a.im;
    
    return (c);
    
}





/*---------------------------------------------------------------------------------
|
|   ComplexDivision
|
|   Returns the complex quotient of two complex numbers.
|
---------------------------------------------------------------------------------*/
complex ComplexDivision (complex a, complex b)

{

    complex		c;
    MrBFlt 		r, den;
    
	if(fabs(b.re) >= fabs(b.im)) 
		{
		r = b.im / b.re;
		den = b.re + r * b.im;
		c.re = (a.re + r * a.im) / den;
		c.im = (a.im - r * a.re) / den;
		} 
	else
		{
		r = b.re / b.im;
		den = b.im + r * b.re;
		c.re = (a.re * r + a.im) / den;
		c.im = (a.im * r - a.re) / den;
		}
	
	return (c);
	
}





/*---------------------------------------------------------------------------------
|
|   ComplexDivision2
|
|   Returns the complex quotient of two complex numbers. It does not require that
|   the numbers be in a complex structure.
|
---------------------------------------------------------------------------------*/
void ComplexDivision2 (MrBFlt ar, MrBFlt ai, MrBFlt br, MrBFlt bi, MrBFlt *cr, MrBFlt *ci)

{

	MrBFlt		s, ais, bis, ars, brs;

	s = fabs(br) + fabs(bi);
	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;
	s = brs*brs + bis*bis;
	*cr = (ars*brs + ais*bis) / s;
	*ci = (ais*brs - ars*bis) / s;
	
}





/*---------------------------------------------------------------------------------
|
|   ComplexExponentiation
|
|   Returns the complex exponential of a complex number.
|
---------------------------------------------------------------------------------*/
complex ComplexExponentiation (complex a)

{

	complex		c;

	c.re = exp(a.re);
	if (AreDoublesEqual(a.im,0.0, ETA)==YES) /* == 0 */
		c.im = 0; 
	else
		{ 
		c.im = c.re*sin(a.im); 
		c.re *= cos(a.im); 
		}

	return (c);
   
}




/*---------------------------------------------------------------------------------
|
|   ComplexInvertMatrix
|
|   Inverts a matrix of complex numbers using the LU-decomposition method. 
|   The program has the following variables:
|
|      a        -- the matrix to be inverted
|      aInverse -- the results of the matrix inversion
|      dim      -- the dimension of the square matrix a and its inverse
|      dwork    -- a work vector of doubles
|      indx     -- a work vector of integers
|      col      -- carries the results of the back substitution
|      
|   The function returns YES (1) or NO (0) if the results are singular.
|
---------------------------------------------------------------------------------*/
int ComplexInvertMatrix (int dim, complex **a, MrBFlt *dwork, int *indx, complex **aInverse, complex *col)

{

	int             isSingular, i, j;

	isSingular = ComplexLUDecompose (dim, a, dwork, indx, (MrBFlt *)NULL);

	if (isSingular == 0) 
		{
		for (j=0; j<dim; j++) 
			{
			for (i=0; i<dim; i++)
				col[i] = Complex (0.0, 0.0);
			col[j] = Complex (1.0, 0.0);
			ComplexLUBackSubstitution (dim, a, indx, col);
			for (i=0; i<dim; i++)
				aInverse[i][j] = col[i];
			}
		}

	return (isSingular);
    
}





/*---------------------------------------------------------------------------------
|
|   ComplexExponentiation
|
|   Returns the complex exponential of a complex number.
|
---------------------------------------------------------------------------------*/
complex ComplexLog (complex a)

{

	complex 	c;
	
	c.re = log(ComplexAbsoluteValue(a));
	if (AreDoublesEqual(a.re,0.0,ETA)==YES) /* == 0.0 */ 
		{
		c.im = PIOVER2;
		} 
	else 
		{
		c.im = atan2(a.im, a.re);
		}
		
	return (c);
    
}






/*---------------------------------------------------------------------------------
|
|   ComplexLUBackSubstitution
|
|   Perform back-substitution into a LU-decomposed matrix to obtain
|   the inverse.
|      
---------------------------------------------------------------------------------*/
void ComplexLUBackSubstitution (int dim, complex **a, int *indx, complex *b)

{

	int             i, ip, j, ii = -1;
	complex     	sum;

	for (i = 0; i < dim; i++) 
		{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0) 
			{
			for (j = ii; j <= i - 1; j++)
				sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][j], b[j]));
			} 
		else if (AreDoublesEqual(sum.re,0.0,ETA)==NO || AreDoublesEqual(sum.im, 0.0, ETA)==NO) /* 2x != 0.0 */
			ii = i;
		b[i] = sum;
		}
	for (i = dim - 1; i >= 0; i--) 
		{
		sum = b[i];
		for (j = i + 1; j < dim; j++)
			sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][j], b[j]));
		b[i] = ComplexDivision (sum, a[i][i]);
		}
		
}






/*---------------------------------------------------------------------------------
|
|   ComplexLUDecompose
|
|   Replaces the matrix a with its LU-decomposition. 
|   The program has the following variables:
|
|      a        -- the matrix
|      dim      -- the dimension of the square matrix a and its inverse
|      vv       -- a work vector of doubles
|      indx     -- row permutation according to partitial pivoting sequence
|      pd       -- 1 if number of row interchanges was even, -1 if number of
|                  row interchanges was odd. Can be NULL.
|      
|   The function returns YES (1) or NO (0) if the results are singular.
|
---------------------------------------------------------------------------------*/
int ComplexLUDecompose (int dim, complex **a, MrBFlt *vv, int *indx, MrBFlt *pd)

{

	int             i, imax, j, k;
	MrBFlt      	big, dum, temp, d;
	complex			sum, cdum;

	d = 1.0;
	imax = 0;

	for (i = 0; i < dim; i++) 
		{
		big = 0.0;
		for (j = 0; j < dim; j++) 
			{
			if ((temp = ComplexAbsoluteValue (a[i][j])) > big)
				big = temp;
			}
		if (AreDoublesEqual(big, 0.0, ETA)==YES) /* == 0.0 */
			{
			MrBayesPrint ("%s   Error: Problem in ComplexLUDecompose\n", spacer);
			return (1);
			}
		vv[i] = 1.0 / big;
		}

	for (j = 0; j < dim; j++) 
		{
		for (i = 0; i < j; i++) 
			{
			sum = a[i][j];
			for (k = 0; k < i; k++) 
				sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][k], a[k][j]));
			a[i][j] = sum;
			}
		big = 0.0;
		for (i = j; i < dim; i++) 
			{
			sum = a[i][j];
			for (k = 0; k < j; k++)
			sum = ComplexSubtraction (sum, ComplexMultiplication (a[i][k], a[k][j]));
			a[i][j] = sum;
			dum = vv[i] * ComplexAbsoluteValue (sum);
			if (dum >= big) 
				{
				big = dum;
				imax = i;
				}
			}
		if (j != imax) 
			{
			for (k = 0; k < dim; k++) 
				{
				cdum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = cdum;
				}       
			d = -d;
			vv[imax] = vv[j];
			}
		indx[j] = imax;
		if (AreDoublesEqual(a[j][j].re, 0.0, ETA)==YES && AreDoublesEqual(a[j][j].im, 0.0, ETA)==YES) /* 2x == 0.0 */
			a[j][j] = Complex (1.0e-20, 1.0e-20);
		if (j != dim - 1)
			{
			cdum = ComplexDivision (Complex(1.0, 0.0), a[j][j]);
			for (i = j + 1; i < dim; i++)
			a[i][j] = ComplexMultiplication (a[i][j], cdum);
			}
		}

	if (pd != NULL)
		*pd = d;
		
	return (0);
	
}





/*---------------------------------------------------------------------------------
|
|   ComplexMultiplication
|
|   Returns the complex product of two complex numbers.
|
---------------------------------------------------------------------------------*/
complex ComplexMultiplication (complex a, complex b)

{

    complex 	c;
    
    c.re = a.re * b.re - a.im * b.im;
    c.im = a.im * b.re + a.re * b.im;
    
    return (c);
    
}





/*---------------------------------------------------------------------------------
|
|   ComplexSquareRoot
|
|   Returns the complex square root of a complex number.
|
---------------------------------------------------------------------------------*/
complex ComplexSquareRoot (complex a)

{

    complex 		c;
    MrBFlt 			x, y, w, r;
    
    if (AreDoublesEqual(a.re, 0.0, ETA)==YES && AreDoublesEqual(a.im, 0.0, ETA)==YES) /* 2x == 0.0 */
    	{
        c.re = 0.0;
        c.im = 0.0;
        return (c);
    	}
    else
    	{
        x = fabs(a.re);
        y = fabs(a.im);
        if (x >= y)
        	{
            r = y / x;
            w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
        	}
        else
        	{
            r = x / y;
            w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
        	}
        if (a.re >= 0.0)
        	{
            c.re = w;
            c.im = a.im / (2.0 * w);
        	}
        else
        	{
            c.im = (a.im >= 0.0) ? w : -w;
            c.re = a.im / (2.0 * c.im);
        	}
        return (c);
		}
    
}





/*---------------------------------------------------------------------------------
|
|   ComplexSubtraction
|
|   Returns the complex difference of two complex numbers.
|
---------------------------------------------------------------------------------*/
complex ComplexSubtraction (complex a, complex b)

{

    complex 	c;
    
    c.re = a.re - b.re;
    c.im = a.im - b.im;
    
    return (c);
    
}





/*---------------------------------------------------------------------------------
|
|   ComputeEigenSystem
|
|   Calculates the eigenvalues, eigenvectors, and the inverse of the eigenvectors
|   for a matrix of real numbers.
|
---------------------------------------------------------------------------------*/
int ComputeEigenSystem (int dim, MrBFlt **a, MrBFlt *v, MrBFlt *vi, MrBFlt **u, int *iwork, MrBFlt *dwork)

{

	int			i, rc;

	rc = EigensForRealMatrix (dim, a, v, vi, u, iwork, dwork);
	if (rc != NO_ERROR)
		{
		MrBayesPrint ("%s   Error in ComputeEigenSystem.\n", spacer);
		return (ERROR);
		}
	for (i=0; i<dim; i++)
		{
		if (AreDoublesEqual(vi[i], 0.0, ETA)==NO) /* != 0.0 */
			return (EVALUATE_COMPLEX_NUMBERS);
		}

	return (NO_ERROR);
	
}





/*---------------------------------------------------------------------------------
|
|   ComputeLandU
|
|   This function computes the L and U decomposition of a matrix. Basically,
|   we find matrices lMat and uMat such that
|
|      lMat * uMat = aMat
|
---------------------------------------------------------------------------------*/
void ComputeLandU (int dim, MrBFlt **aMat, MrBFlt **lMat, MrBFlt **uMat)

{

	int			i, j, k, m, row, col;

	for (j=0; j<dim; j++) 
		{
		for (k=0; k<j; k++)
			for (i=k+1; i<j; i++)
				aMat[i][j] = aMat[i][j] - aMat[i][k] * aMat[k][j];

		for (k=0; k<j; k++)
			for (i=j; i<dim; i++)
				aMat[i][j] = aMat[i][j] - aMat[i][k]*aMat[k][j];

		for (m=j+1; m<dim; m++)
  		    aMat[m][j] /= aMat[j][j]; 
		}

	for (row=0; row<dim; row++)
		{
		for (col=0; col<dim; col++) 
			{
			if (row <= col) 
				{
				uMat[row][col] = aMat[row][col];
				lMat[row][col] = (row == col ? 1.0 : 0.0);
				}
			else 
				{
				lMat[row][col] = aMat[row][col];
				uMat[row][col] = 0.0;
				}
			}
		}

}





/*---------------------------------------------------------------------------------
|
|   ComputeMatrixExponential
|
|   The method approximates the matrix exponential, f = e^a, using
|   the algorithm 11.3.1, described in:
|  
|   Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
|      The Johns Hopkins University Press, Baltimore, Maryland.
|
|   The method has the advantage of error control. The error is controlled by
|   setting qValue appropriately (using the function SetQValue).
|
---------------------------------------------------------------------------------*/
void ComputeMatrixExponential (int dim, MrBFlt **a, int qValue, MrBFlt **f)

{

	int			i, j, k, negativeFactor;
	MrBFlt 		maxAValue, c, **d, **n, **x, **cX;

	d  = AllocateSquareDoubleMatrix (dim);
	n  = AllocateSquareDoubleMatrix (dim);
	x  = AllocateSquareDoubleMatrix (dim);
	cX = AllocateSquareDoubleMatrix (dim);

	SetToIdentity (dim, d);
	SetToIdentity (dim, n);
	SetToIdentity (dim, x);

	maxAValue = 0;
	for (i=0; i<dim; i++)
		maxAValue = MAX (maxAValue, a[i][i]);

	j = MAX (0, LogBase2Plus1 (maxAValue));

	DivideByTwos (dim, a, j);
	
	c = 1;
	for (k=1; k<=qValue; k++) 
		{
		c = c * (qValue - k + 1.0) / ((2.0 * qValue - k + 1.0) * k);

		/* X = AX */
		MultiplyMatrices (dim, a, x, x);

		/* N = N + cX */
		MultiplyMatrixByScalar (dim, x, c, cX);
		AddTwoMatrices (dim, n, cX, n);

		/* D = D + (-1)^k*cX */
		negativeFactor = (k % 2 == 0 ? 1 : -1);
		if (negativeFactor == -1)
			MultiplyMatrixByScalar (dim, cX, negativeFactor, cX);
		AddTwoMatrices (dim, d, cX, d);      
		}

	GaussianElimination (dim, d, n, f);

	for (k = 0; k < j; k++)
		MultiplyMatrices (dim, f, f, f);
	
	for (i=0; i<dim; i++)
		{
		for (j=0; j<dim; j++)
			{
			if (f[i][j] < 0.0)
				f[i][j] = fabs(f[i][j]);
			}
		}
		
	FreeSquareDoubleMatrix (d);
	FreeSquareDoubleMatrix (n);
	FreeSquareDoubleMatrix (x);
	FreeSquareDoubleMatrix (cX);
	
}





/*---------------------------------------------------------------------------------
|
|   CopyComplexMatrices
|
|   Copies the contents of one matrix of complex numbers to another matrix.
|
---------------------------------------------------------------------------------*/
void CopyComplexMatrices (int dim, complex **from, complex **to)

{

	int			i, j;
	
	for (i=0; i<dim; i++)
		{
		for (j=0; j<dim; j++) 
			{
			to[i][j].re = from[i][j].re;
			to[i][j].im = from[i][j].im;
			}
		}

}





/*---------------------------------------------------------------------------------
|
|   CopyDoubleMatrices
|
|   Copies the contents of one matrix of doubles to another matrix.
|
---------------------------------------------------------------------------------*/
void CopyDoubleMatrices (int dim, MrBFlt **from, MrBFlt **to)

{

	int			i, j;
	
	for (i=0; i<dim; i++)
		{
		for (j=0; j<dim; j++) 
			{
			to[i][j] = from[i][j];
			}
		}
		
}





/*---------------------------------------------------------------------------------
|
|   DirichletRandomVariable
|
|   Generate a Dirichlet-distributed random variable. The parameter of the
|   Dirichlet is contained in the vector alp. The random variable is contained
|   in the vector z.
|      
---------------------------------------------------------------------------------*/
void DirichletRandomVariable (MrBFlt *alp, MrBFlt *z, int n, safeLong *seed)

{

	int		i;
	MrBFlt	sum;

	sum = 0.0;
	for(i=0; i<n; i++)
		{
		z[i] = RndGamma (alp[i], seed) / 1.0;
		sum += z[i];
		}
	for(i=0; i<n; i++)
		z[i] /= sum;
	
}





/*---------------------------------------------------------------------------------
|
|   DiscreteGamma
|
|   Discretization of gamma distribution with equal proportions in each
|   category.
|
---------------------------------------------------------------------------------*/
int DiscreteGamma (MrBFlt *rK, MrBFlt alfa, MrBFlt beta, int K, int median)

{

	int 			i;
	MrBFlt 			gap05 = 1.0/(2.0*K), t, factor = alfa/beta*K, lnga1;

	if (median) 
		{
		for (i=0; i<K; i++) 
			rK[i] = POINTGAMMA((i*2.0+1.0)*gap05, alfa, beta);
		for (i=0,t=0; i<K; i++) 
			t += rK[i];
		for (i=0; i<K; i++)     
			rK[i] *= factor / t;
		}
	else 
		{
		lnga1 = LnGamma(alfa+1);
		/* calculate the points in the gamma distribution */
		for (i=0; i<K-1; i++) 
			rK[i] = POINTGAMMA((i+1.0)/K, alfa, beta);
		/* calculate the cumulative values */
		for (i=0; i<K-1; i++) 
			rK[i] = IncompleteGamma(rK[i] * beta, alfa + 1.0, lnga1);
		rK[K-1] = 1.0;
		/* calculate the relative values and rescale */
		for (i=K-1; i>0; i--)
			{
			rK[i] -= rK[i-1];
			rK[i] *= factor;
			}
		rK[0] *= factor;
		}

	return (NO_ERROR);
	
}





/*---------------------------------------------------------------------------------
|
|   DivideByTwos
|
|   Divides all of the elements of the matrix a by 2^power.
|      
---------------------------------------------------------------------------------*/
void DivideByTwos (int dim, MrBFlt **a, int power)

{

	int			divisor = 1, i, row, col;

	for (i=0; i<power; i++)
		divisor = divisor * 2;

	for (row=0; row<dim; row++)
		for (col=0; col<dim; col++)
  		    a[row][col] /= divisor;

}





/*---------------------------------------------------------------------------------
|
|   D_sign
|
|   This function is called from "Hqr2".
|
---------------------------------------------------------------------------------*/
MrBFlt D_sign (MrBFlt a, MrBFlt b)

{

	MrBFlt		x;

	x = (a >= 0 ? a : -a);
	
	return (b >= 0 ? x : -x);
	
}





/*---------------------------------------------------------------------------------
|
|   Eigens
|
|   The matrix of interest is a. The ouptut is the real and imaginary parts of the 
|   eigenvalues (wr and wi). z contains the real and imaginary parts of the 
|   eigenvectors. iv2 and fv1 are working vectors.
|      
---------------------------------------------------------------------------------*/
int EigensForRealMatrix (int dim, MrBFlt **a, MrBFlt *wr, MrBFlt *wi, MrBFlt **z, int *iv1, MrBFlt *fv1)

{

	static int	is1, is2;
	int			ierr;

	Balanc (dim, a, &is1, &is2, fv1);
	ElmHes (dim, is1, is2, a, iv1);
	ElTran (dim, is1, is2, a, iv1, z);
	ierr = Hqr2 (dim, is1, is2, a, wr, wi, z);
	if (ierr == 0)
		BalBak (dim, is1, is2, fv1, dim, z);

	return (ierr);
	
}





/*---------------------------------------------------------------------------------
|
|   ElmHes
|
|   Given a real general matrix, this subroutine
|   reduces a submatrix situated in rows and columns
|   low through high to upper Hessenberg form by
|   stabilized elementary similarity transformations.
|
|   On input:
|
|    * dim is the order of the matrix
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc.  if  balanc  has not been used,
|      set low=1, high=dim.
|
|    * a contains the input matrix.
|
|   On output:
|
|    * a contains the hessenberg matrix.  The multipliers
|      which were used in the reduction are stored in the
|      remaining triangle under the hessenberg matrix.
|
|    * interchanged contains information on the rows and columns
|      interchanged in the reduction.
|
|   Only elements low through high are used.
|
---------------------------------------------------------------------------------*/
void ElmHes (int dim, int low, int high, MrBFlt **a, int *interchanged)

{
	int			i, j, m, la, mm1, kp1, mp1;
	MrBFlt		x, y;
	
	la = high - 1;
	kp1 = low + 1;
	if (la < kp1)
		return; /* remove goto statement, which exits at bottom of function */

	for (m=kp1; m<=la; m++)
		{
		mm1 = m - 1;
		x = 0.0;
		i = m;
	
		for (j=m; j<=high; j++)
			{
			if (fabs(a[j][mm1]) > fabs(x)) /* change direction of inequality */
				{                          /* remove goto statement          */
				x = a[j][mm1];
				i = j;
				}
			}
	
		interchanged[m] = i;
		if (i != m) /* change "==" to "!=", eliminating goto statement */
			{
			/* interchange rows and columns of a */
			for (j=mm1; j<dim; j++)
				{
				y = a[i][j];
				a[i][j] = a[m][j];
				a[m][j] = y;
				}
			for (j=0; j<=high; j++)
				{
				y = a[j][i];
				a[j][i] = a[j][m];
				a[j][m] = y;
				}
			}

		if (AreDoublesEqual(x, 0.0, ETA)==NO) /* change "==" to "!=", eliminating goto statement */
			{
			mp1 = m + 1;
		
			for (i=mp1; i<=high; i++)
				{
				y = a[i][mm1];
				if (AreDoublesEqual(y, 0.0, ETA)==NO) /* != 0.0 */
					{
					y /= x;
					a[i][mm1] = y;
					for (j = m; j < dim; j++)
						a[i][j] -= y * a[m][j];
					for (j = 0; j <= high; j++)
						a[j][m] += y * a[j][i];
					}
				}
			}
		}

#if 0
/* begin f2c version of code:
   elmhes.f -- translated by f2c (version 19971204) */
int elmhes (int *nm, int *n, int *low, int *igh, MrBFlt *a, int *int__)

{

    /*system generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    MrBFlt d__1;

    /* local variables */
    static int i__, j, m;
    static MrBFlt x, y;
    static int la, mm1, kp1, mp1;

    /* parameter adjustments */
    a_dim1 = *nm;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --int__;

    /* function body */
    la = *igh - 1;
    kp1 = *low + 1;
    if (la < kp1) 
		goto L200;

    i__1 = la;
    for (m = kp1; m <= i__1; ++m) 
	    {
		mm1 = m - 1;
		x = 0.;
		i__ = m;
		i__2 = *igh;
		for (j = m; j <= i__2; ++j) 
			{
		    if ((d__1 = a[j + mm1 * a_dim1], abs(d__1)) <= abs(x)) 
				goto L100;
		    x = a[j + mm1 * a_dim1];
		    i__ = j;
			L100:
		    ;
		}

	int__[m] = i__;
	if (i__ == m) 
	    goto L130;

	/* .......... interchange rows and columns of a.......... */
	i__2 = *n;
	for (j = mm1; j <= i__2; ++j) 
		{
	    y = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = a[m + j * a_dim1];
	    a[m + j * a_dim1] = y;
		/* L110: */
		}

	i__2 = *igh;
	for (j = 1; j <= i__2; ++j) 
		{
	    y = a[j + i__ * a_dim1];
	    a[j + i__ * a_dim1] = a[j + m * a_dim1];
	    a[j + m * a_dim1] = y;
		/* L120: */
		}
		
	/* .......... end interchange .......... */
	L130:
	if (x == 0.) 
	    goto L180;
	mp1 = m + 1;

	i__2 = *igh;
	for (i__ = mp1; i__ <= i__2; ++i__) 
		{
	    y = a[i__ + mm1 * a_dim1];
	    if (y == 0.) 
			goto L160;
	    y /= x;
	    a[i__ + mm1 * a_dim1] = y;

	    i__3 = *n;
	    for (j = m; j <= i__3; ++j) 
	    	{
			/* L140: */
			a[i__ + j * a_dim1] -= y * a[m + j * a_dim1];
			}

		i__3 = *igh;
		for (j = 1; j <= i__3; ++j) 
			{
			/* L150: */
			a[j + m * a_dim1] += y * a[j + i__ * a_dim1];
			}

		L160:
			;
		}

	L180:
		;
	}

	L200:
	return 0;
	
}
/* end f2c version of code */
#endif
		
}





/*---------------------------------------------------------------------------------
|
|   ElTran
|
|   This subroutine accumulates the stabilized elementary
|   similarity transformations used in the reduction of a
|   real general matrix to upper Hessenberg form by ElmHes.
|
|   On input:
|
|    * dim is the order of the matrix.
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc. If Balanc has not been used,
|      set low=0, high=dim-1.
|
|    * a contains the multipliers which were used in the
|      reduction by  ElmHes in its lower triangle
|      below the subdiagonal.
|
|    * interchanged contains information on the rows and columns
|      interchanged in the reduction by ElmHes.
|      only elements low through high are used.
|
|   On output:
|
|    * z contains the transformation matrix produced in the
|      reduction by ElmHes.
|
|   This routine is a translation of the Algol procedure from
|   Handbook for Automatic Computation, vol. II, Linear Algebra,
|   by Wilkinson and Reinsch, Springer-Verlag.
|   
---------------------------------------------------------------------------------*/
void ElTran (int dim, int low, int high, MrBFlt **a, int *interchanged, MrBFlt **z)

{

	int			i, j, mp;

	/* initialize z to identity matrix */
	for (j=0; j<dim; j++)
		{
		for (i=0; i<dim; i++)
			z[i][j] = 0.0;
		z[j][j] = 1.0;
		}
	for (mp=high-1; mp>=low+1; mp--) /* there were a number of additional    */
		{                            /* variables (kl, la, m, mm, mp1) that  */
		for (i=mp+1; i<=high; i++)   /* have been eliminated here simply by  */
			z[i][mp] = a[i][mp-1];   /* initializing variables appropriately */
		i = interchanged[mp];        /* in the loops                         */
		if (i != mp) /* change "==" to "!=" to eliminate a goto statement */
			{
			for (j=mp; j<=high; j++)
				{
				z[mp][j] = z[i][j];
				z[i][j] = 0.0;
				}
			z[i][mp] = 1.0;
			}
		}
	
#if 0
/* begin f2c version of code:
   eltran.f -- translated by f2c (version 19971204) */
int eltran (int *nm, int *n, int *low, int *igh, MrBFlt *a, int *int__, MrBFlt *z__)

{

	/* system generated locals */
	int a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;

	/* local variables */
	static int i__, j, kl, mm, mp, mp1;

	/*     .......... initialize z to identity matrix .......... */
	
	/* parameter adjustments */
	z_dim1 = *nm;
	z_offset = z_dim1 + 1;
	z__ -= z_offset;
	--int__;
	a_dim1 = *nm;
	a_offset = a_dim1 + 1;
	a -= a_offset;

	/* function Body */
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) 
		{
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) 
			{
			/* L60: */
			z__[i__ + j * z_dim1] = 0.0;
			}
		z__[j + j * z_dim1] = 1.0;
		/* L80: */
		}

	kl = *igh - *low - 1;
	if (kl < 1) 
		goto L200;

	/* .......... for mp=igh-1 step -1 until low+1 do -- .......... */
	i__1 = kl;
	for (mm = 1; mm <= i__1; ++mm) 
		{
		mp = *igh - mm;
		mp1 = mp + 1;
		i__2 = *igh;
		for (i__ = mp1; i__ <= i__2; ++i__) 
			{
			/* L100: */
			z__[i__ + mp * z_dim1] = a[i__ + (mp - 1) * a_dim1];
			}
		i__ = int__[mp];
		if (i__ == mp) 
			goto L140;
		i__2 = *igh;
		for (j = mp; j <= i__2; ++j) 
			{
			z__[mp + j * z_dim1] = z__[i__ + j * z_dim1];
			z__[i__ + j * z_dim1] = 0.;
			/* L130: */
			}
		z__[i__ + mp * z_dim1] = 1.;
		L140:
			;
		}

	L200:
	return 0;

}
/* end f2c version of code */
#endif
	
}





/*---------------------------------------------------------------------------------
|
|   Exchange
|
---------------------------------------------------------------------------------*/
void Exchange (int j, int k, int l, int m, int n, MrBFlt **a, MrBFlt *scale)

{

	int			i;
	MrBFlt		f;

	scale[m] = (MrBFlt)j;
	if (j != m)
		{
		for (i = 0; i <= l; i++)
			{
			f = a[i][j];
			a[i][j] = a[i][m];
			a[i][m] = f;
			}	
		for (i = k; i < n; i++)
			{
			f = a[j][i];
			a[j][i] = a[m][i];
			a[m][i] = f;
			}
		}
		
}





/*---------------------------------------------------------------------------------
|
|   Factorial
|
|   Returns x!
|      
---------------------------------------------------------------------------------*/
MrBFlt Factorial (int x)

{

	int 		i;
	MrBFlt		fac;
	
	fac = 1.0;
	for (i=0; i<x; i++)
		{
		fac *= (i+1);
		}
		
	return (fac);
	
}





/*---------------------------------------------------------------------------------
|
|   ForwardSubstitutionRow
|
---------------------------------------------------------------------------------*/
void ForwardSubstitutionRow (int dim, MrBFlt **L, MrBFlt *b)

{

	int			i, j;
	MrBFlt		dotProduct;

	b[0] = b[0] / L[0][0];
	for (i=1; i<dim; i++) 
		{
		dotProduct = 0.0;
		for (j=0; j<i; j++)
      		dotProduct += L[i][j] * b[j];
		b[i] = (b[i] - dotProduct) / L[i][i];
		}
	
}





/*---------------------------------------------------------------------------------
|
|   FreeSquareComplexMatrix
|
|   Frees a matrix of complex numbers.
|      
---------------------------------------------------------------------------------*/
void FreeSquareComplexMatrix (complex **m)

{

    free((char *) (m[0]));
    free((char *) (m));
    
}





/*---------------------------------------------------------------------------------
|
|   FreeSquareDoubleMatrix
|
|   Frees a matrix of doubles.
|      
---------------------------------------------------------------------------------*/
void FreeSquareDoubleMatrix (MrBFlt **m)

{

	free((char *) (m[0]));
	free((char *) (m));
	
}




/*---------------------------------------------------------------------------------
|
|   FreeSquareIntegerMatrix
|
|   Frees a matrix of integers.
|      
---------------------------------------------------------------------------------*/
void FreeSquareIntegerMatrix (int **m)

{

	free((char *) (m[0]));
	free((char *) (m));
	
}





/*---------------------------------------------------------------------------------
|
|   GammaRandomVariable
|
|   This function generates a gamma-distributed random variable with parameters
|   a and b. The mean is E(X) = a / b and the variance is Var(X) = a / b^2.
|      
---------------------------------------------------------------------------------*/
MrBFlt GammaRandomVariable (MrBFlt a, MrBFlt b, safeLong *seed)

{

	return (RndGamma (a, seed) / b);

}





/*---------------------------------------------------------------------------------
|
|   GaussianElimination
|      
---------------------------------------------------------------------------------*/
void GaussianElimination (int dim, MrBFlt **a, MrBFlt **bMat, MrBFlt **xMat)

{

	int			i, k;
	MrBFlt		*bVec, **lMat, **uMat;

	lMat = AllocateSquareDoubleMatrix (dim);
	uMat = AllocateSquareDoubleMatrix (dim);
	bVec = (MrBFlt *)SafeMalloc((size_t) ((dim) * sizeof(MrBFlt)));
	if (!bVec)
		{
		MrBayesPrint ("%s   Error: Problem allocating bVec\n", spacer);
		exit (0);
		}

	ComputeLandU (dim, a, lMat, uMat);

	for (k=0; k<dim; k++) 
		{
		
		for (i=0; i<dim; i++)
			bVec[i] = bMat[i][k];

		/* Answer of Ly = b (which is solving for y) is copied into b. */
		ForwardSubstitutionRow (dim, lMat, bVec);

		/* Answer of Ux = y (solving for x and the y was copied into b above) 
		   is also copied into b. */
		BackSubstitutionRow (dim, uMat, bVec);

		for (i=0; i<dim; i++)
			xMat[i][k] = bVec[i];

		}
	
	FreeSquareDoubleMatrix (lMat);
	FreeSquareDoubleMatrix (uMat);
	free (bVec);

}





/*---------------------------------------------------------------------------------
|
|   GetEigens
|
---------------------------------------------------------------------------------*/
int GetEigens (int dim, MrBFlt **q, MrBFlt *eigenValues, MrBFlt *eigvalsImag, MrBFlt **eigvecs, MrBFlt **inverseEigvecs, complex **Ceigvecs, complex **CinverseEigvecs)

{

	int			i, j, rc, *iWork, isComplex;
	MrBFlt		**tempWork, *dWork;
	complex 	**cWork, *Ccol;

	/* allocate memory */
	dWork = (MrBFlt *)SafeMalloc((size_t) (dim * sizeof(MrBFlt)));
	iWork = (int *)SafeMalloc((size_t) (dim * sizeof(int)));
	if (!dWork || !iWork)
		{
		MrBayesPrint ("%s   Error: Problem in GetEigens\n", spacer);
		exit (0);
		}

	/* calculate eigenvalues and eigenvectors */
	isComplex = NO;
	rc = ComputeEigenSystem (dim, q, eigenValues, eigvalsImag, eigvecs, iWork, dWork);
	if (rc != NO_ERROR)
		{
		if (rc == EVALUATE_COMPLEX_NUMBERS)
			isComplex = YES;
		}

	/* invert eigenvectors */
	if (isComplex == NO)
		{
		tempWork = AllocateSquareDoubleMatrix (dim);
		CopyDoubleMatrices (dim, eigvecs, tempWork);
		InvertMatrix (dim, tempWork, dWork, iWork, inverseEigvecs);
		FreeSquareDoubleMatrix (tempWork);
		}
	else
		{
		for(i=0; i<dim; i++)
			{
			  if (fabs(eigvalsImag[i])<1E-20) /* == 0.0 */
				{ 
				for(j=0; j<dim; j++)
					{
					Ceigvecs[j][i].re = eigvecs[j][i];
					Ceigvecs[j][i].im = 0.0;
					}
				}
			else if (eigvalsImag[i] > 0)
				{ 
				for (j=0; j<dim; j++)
					{
					Ceigvecs[j][i].re = eigvecs[j][i];
					Ceigvecs[j][i].im = eigvecs[j][i + 1];
					}
				}
			else if (eigvalsImag[i] < 0)
				{ 
				for (j=0; j<dim; j++)
					{
					Ceigvecs[j][i].re =  eigvecs[j][i-1];
					Ceigvecs[j][i].im = -eigvecs[j][i];
					}
				}
			}
		Ccol = (complex *)SafeMalloc((size_t) (dim * sizeof(complex)));
		if (!Ccol)
			{
			MrBayesPrint ("%s   Error: Problem in GetEigens\n", spacer);
			exit (0);
			}
		cWork = AllocateSquareComplexMatrix (dim);
		CopyComplexMatrices (dim, Ceigvecs, cWork);
		ComplexInvertMatrix (dim, cWork, dWork, iWork, CinverseEigvecs, Ccol);
		free (Ccol);
		FreeSquareComplexMatrix (cWork);
		}

	free (dWork);
	free (iWork);
	
	return (isComplex);

}





/*---------------------------------------------------------------------------------
|
|   Hqr2
|
|   This subroutine finds the eigenvalues and eigenvectors
|   of a real upper Hessenberg matrix by the QR method. The
|   eigenvectors of a real general matrix can also be found
|   if ElmHes  and ElTran or OrtHes and OrTran have
|   been used to reduce this general matrix to Hessenberg form
|   and to accumulate the similarity transformations.
|
|   On input:
|
|    * dim is the order of the matrix.
|
|    * low and high are integers determined by the balancing
|      subroutine  balanc. If  balanc has not been used,
|      set low=0, high=dim-1.
|
|    * h contains the upper hessenberg matrix. Information about
|      the transformations used in the reduction to Hessenberg
|      form by  ElmHes  or OrtHes, if performed, is stored
|      in the remaining triangle under the Hessenberg matrix.
|
|   On output:
|
|    * h has been destroyed.
|
|    * wr and wi contain the real and imaginary parts,
|      respectively, of the eigenvalues. The eigenvalues
|      are unordered except that complex conjugate pairs
|      of values appear consecutively with the eigenvalue
|      having the positive imaginary part first. If an
|      error exit is made, the eigenvalues should be correct
|      for indices j,...,dim-1.
|
|    * z contains the transformation matrix produced by ElTran
|      after the reduction by ElmHes, or by OrTran after the
|      reduction by OrtHes, if performed. If the eigenvectors
|      of the Hessenberg matrix are desired, z must contain the
|      identity matrix.
|
|   Calls ComplexDivision2 for complex division.
|
|   This function returns:
|      zero       for normal return,
|      j          if the limit of 30*n iterations is exhausted
|                 while the j-th eigenvalue is being sought.
|
|   This subroutine is a translation of the ALGOL procedure HQR2,
|   Num. Math. 14, 219,231(1970) by Martin, Peters, and Wilkinson.
|   Handbook for Automatic Computation, vol. II - Linear Algebra,
|   pp. 357-391 (1971).
|   
---------------------------------------------------------------------------------*/
int Hqr2 (int dim, int low, int high, MrBFlt **h, MrBFlt *wr, MrBFlt *wi, MrBFlt **z)

{

	int			i, j, k, l, m, na, en, notlas, mp2, itn, its, enm2, twoRoots;
	MrBFlt		norm, p=0.0, q=0.0, r=0.0, s=0.0, t, w=0.0, x, y=0.0, ra, sa, vi, vr, zz=0.0, tst1, tst2;

	norm = 0.0;
	k = 0;  /* used for array indexing. FORTRAN version: k = 1 */
	
	/* store roots isolated by balance, and compute matrix norm */
	for (i=0; i<dim; i++)
		{
		for (j=k; j<dim; j++)
			norm += fabs(h[i][j]);

		k = i;
		if ((i < low) || (i > high))
			{
			wr[i] = h[i][i];
			wi[i] = 0.0;
			}
		}
	en = high;
	t = 0.0;
	itn = dim * 30;

	/* search for next eigenvalues */
	while (en >= low) /* changed from an "if(en < lo)" to eliminate a goto statement */
		{
		its = 0;
		na = en - 1;
		enm2 = na - 1;
		twoRoots = FALSE;

		for (;;)
			{
			for (l=en; l>low; l--) /* changed indexing, got rid of lo, ll */
				{
				s = fabs(h[l-1][l-1]) + fabs(h[l][l]);
				if (AreDoublesEqual(s, 0.0, ETA)==YES) /* == 0.0 */
					s = norm;
				tst1 = s;
				tst2 = tst1 + fabs(h[l][l-1]);
				if (fabs(tst2 - tst1) < ETA) /* tst2 == tst1 */
					break; /* changed to break to remove a goto statement */
				}
	
			/* form shift */
			x = h[en][en];
			if (l == en) /* changed to break to remove a goto statement */
				break;
			y = h[na][na];
			w = h[en][na] * h[na][en];
			if (l == na)         /* used to return to other parts of the code */
				{
				twoRoots = TRUE;
				break;
				}
			if (itn == 0)
				return (en);
				
			/* form exceptional shift */
			if ((its == 10) || (its == 20)) /* changed to remove a goto statement */
				{
				t += x;
				for (i = low; i <= en; i++)
					h[i][i] -= x;
				s = fabs(h[en][na]) + fabs(h[na][enm2]);
				x = 0.75 * s;
				y = x;
				w = -0.4375 * s * s;
				}
			its++;
			itn--;
			
			/* look for two consecutive small sub-diagonal elements */
			for (m=enm2; m>=l; m--)
				{
				/* removed m = enm2 + l - mm and above loop to remove variables */
				zz = h[m][m];
				r = x - zz;
				s = y - zz;
				p = (r * s - w) / h[m+1][m] + h[m][m+1];
				q = h[m+1][m+1] - zz - r - s;
				r = h[m+2][m+1];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l)
					break; /* changed to break to remove a goto statement */
				tst1 = fabs(p) * (fabs(h[m-1][m-1]) + fabs(zz) + fabs(h[m+1][m+1]));
				tst2 = tst1 + fabs(h[m][m-1]) * (fabs(q) + fabs(r));
				if (fabs(tst2 - tst1) < ETA) /* tst2 == tst1 */
					break; /* changed to break to remove a goto statement */
				}
		
			mp2 = m + 2;
			for (i = mp2; i <= en; i++)
				{
				h[i][i-2] = 0.0;
				if (i != mp2) /* changed "==" to "!=" to remove a goto statement */
					h[i][i-3] = 0.0;
				}
	
			/* MrBFlt QR step involving rows l to en and columns m to en */
			for (k=m; k<=na; k++)
				{
				notlas = (k != na);
				if (k != m) /* changed "==" to "!=" to remove a goto statement */
					{
					p = h[k][k-1];
					q = h[k+1][k-1];
					r = 0.0;
					if (notlas)
						r = h[k+2][k-1];
					x = fabs(p) + fabs(q) + fabs(r);
					if (x < ETA) /* == 0.0 */
						continue; /* changed to continue remove a goto statement */
					p /= x;
					q /= x;
					r /= x;
					}
	
		        /*s = sqrt(p*p+q*q+r*r);
		        sgn = (p<0)?-1:(p>0);
		        s = sgn*sqrt(p*p+q*q+r*r);*/
				s = D_sign(sqrt(p*p + q*q + r*r), p);
				if (k != m) /* changed "==" to "!=" to remove a goto statement */
					h[k][k-1] = -s * x;
				else if (l != m) /* else if gets rid of another goto statement */
					h[k][k-1] = -h[k][k-1];
				p += s;
				x = p / s;
				y = q / s;
				zz = r / s;
				q /= p;
				r /= p;
				if (!notlas) /* changed to !notlas to remove goto statement (see **) */
					{
					/* row modification */
					for (j=k; j<dim; j++)
						{
						p = h[k][j] + q * h[k+1][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						} 
					j = MIN(en, k + 3);
					
					/* column modification */
					for (i=0; i<=j; i++)
						{
						p = x * h[i][k] + y * h[i][k+1];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						}
						
					/* accumulate transformations */
					for (i=low; i<=high; i++)
						{
						p = x * z[i][k] + y * z[i][k+1];
						z[i][k] -= p;
						z[i][k+1] -= p * q;
						}
					}
				else /* (**) also put in else */
					{
					/* row modification */
					for (j=k; j<dim; j++)
						{
						p = h[k][j] + q * h[k+1][j] + r * h[k+2][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						h[k+2][j] -= p * zz;
						}
					j = MIN(en, k + 3);
					
					/* column modification */
					for (i = 0; i <= j; i++)
						{
						p = x * h[i][k] + y * h[i][k+1] + zz * h[i][k+2];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						h[i][k+2] -= p * r;
						}
						
					/* accumulate transformations */
					for (i = low; i <= high; i++)
						{
						p = x * z[i][k] + y * z[i][k+1] + zz * z[i][k+2];
						z[i][k] -= p;
						z[i][k+1] -= p * q;
						z[i][k+2] -= p * r;
						}
					}
				}
			}

		if (twoRoots)
			{
			/* two roots found */
			p = (y - x) / 2.0;
			q = p * p + w;
			zz = sqrt(fabs(q));
			h[en][en] = x + t;
			x = h[en][en];
			h[na][na] = y + t;
			if (q >= -1e-12) /* change "<" to ">=", and also change "0.0" to */
				{            /* a small number (Swofford's change)           */
				/* real pair */
				zz = p + D_sign(zz, p);
				wr[na] = x + zz;
				wr[en] = wr[na];
				if (fabs(zz) > ETA) /* != 0.0 */
					wr[en] = x - w/zz;
				wi[na] = 0.0;
				wi[en] = 0.0;
				x = h[en][na];
				s = fabs(x) + fabs(zz);
				p = x / s;
				q = zz / s;
				r = sqrt(p*p + q*q);
				p /= r;
				q /= r;
				
				/* row modification */
				for (j=na; j<dim; j++)
					{
					zz = h[na][j];
					h[na][j] = q * zz + p * h[en][j];
					h[en][j] = q * h[en][j] - p * zz;
					}
					
				/* column modification */
				for (i = 0; i <= en; i++)
					{
					zz = h[i][na];
					h[i][na] = q * zz + p * h[i][en];
					h[i][en] = q * h[i][en] - p * zz;
					}
					
				/* accumulate transformations */
				for (i = low; i <= high; i++)
					{
					zz = z[i][na];
					z[i][na] = q * zz + p * z[i][en];
					z[i][en] = q * z[i][en] - p * zz;
					}
				}
			else
				{
				/* complex pair */
				wr[na] = x + p;
				wr[en] = x + p;
				wi[na] = zz;
				wi[en] = -zz;
				}
			en = enm2;
			}
		else
			{
			/* one root found */
			h[en][en] = x + t;
			wr[en] = h[en][en];
			wi[en] = 0.0;
			en = na;
			}
		}
	
	if (fabs(norm) < ETA) /* == 0.0 */
		return (0); /* was a goto end of function */

	for (en=dim-1; en>=0; en--)
		{
		/*en = n - nn - 1; and change for loop */
		p = wr[en];
		q = wi[en];
		na = en - 1;

		if (q < -1e-12)
			{
			/* last vector component chosen imaginary so that eigenvector
			   matrix is triangular */
			m = na;
			if (fabs(h[en][na]) > fabs(h[na][en]))
				{
				h[na][na] = q / h[en][na];
				h[na][en] = -(h[en][en] - p) / h[en][na];
				}
			else
				ComplexDivision2 (0.0, -h[na][en], h[na][na] - p, q, &h[na][na], &h[na][en]);

			h[en][na] = 0.0;
			h[en][en] = 1.0;
			enm2 = na - 1;
			if (enm2 >= 0) /* changed direction to remove goto statement */
				{
				for (i=enm2; i>=0; i--)
					{
					w = h[i][i] - p;
					ra = 0.0;
					sa = 0.0;
			
					for (j=m; j<=en; j++)
						{
						ra += h[i][j] * h[j][na];
						sa += h[i][j] * h[j][en];
						}
			
					if (wi[i] < 0.0) /* changed direction to remove goto statement */
						{
						zz = w;
						r = ra;
						s = sa;
						}
					else
						{
						m = i;
						if (fabs(wi[i])<ETA) /* == 0.0 */ /* changed direction to remove goto statement */
							ComplexDivision2 (-ra, -sa, w, q, &h[i][na], &h[i][en]);
						else
							{
							/* solve complex equations */
							x = h[i][i+1];
							y = h[i+1][i];
							vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
							vi = (wr[i] - p) * 2.0 * q;
							if ((fabs(vr)<ETA) && (fabs(vi)<ETA))
								{
								tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
								vr = tst1;
								do	{
									vr *= .01;
									tst2 = tst1 + vr;
									}
									while (tst2 > tst1); /* made into a do/while loop */
								}
							ComplexDivision2 (x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, &h[i][na], &h[i][en]);
							if (fabs(x) > fabs(zz) + fabs(q)) /* changed direction to remove goto statement */
								{
								h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;
								h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;
								}
							else
								ComplexDivision2 (-r - y * h[i][na], -s - y * h[i][en], zz, q, &h[i+1][na], &h[i+1][en]);
							}
							
						/* overflow control */
						tst1 = fabs(h[i][na]);
						tst2 = fabs(h[i][en]);
						t = MAX(tst1, tst2);
						if (t > ETA) /* t != 0.0 */
							{
							tst1 = t;
							tst2 = tst1 + 1.0 / tst1;
							if (tst2 <= tst1)
								{
								for (j = i; j <= en; j++)
									{
									h[j][na] /= t;
									h[j][en] /= t;
									}
								}
							}
						}
					}
				}
			}
		else if (fabs(q)<ETA)
			{
			/* real vector */
			m = en;
			h[en][en] = 1.0;
			if (na >= 0)
				{
				for (i=na; i>=0; i--)
					{
					w = h[i][i] - p;
					r = 0.0;
					for (j = m; j <= en; j++)
						r += h[i][j] * h[j][en];
					if (wi[i] < 0.0) /* changed direction to remove goto statement */
						{
						zz = w;
						s = r;
						continue;  /* changed to continue to remove goto statement */
						}
					else
						{
						m = i;
						if (fabs(wi[i])<ETA) /* changed to remove goto statement */
							{
							t = w;
							if (fabs(t)<ETA)  /* changed to remove goto statement */
								{
								tst1 = norm;
								t = tst1;
								do	{
									t *= .01;
									tst2 = norm + t;
									}
									while (tst2 > tst1);
								}			
							h[i][en] = -r / t;
							}
						else
							{
							/* solve real equations */
							x = h[i][i+1];
							y = h[i+1][i];
							q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
							t = (x * s - zz * r) / q;
							h[i][en] = t;
							if (fabs(x) > fabs(zz))  /* changed direction to remove goto statement */
								h[i+1][en] = (-r - w * t) / x;
							else
								h[i+1][en] = (-s - y * t) / zz;
							}
				
						/* overflow control */
						t = fabs(h[i][en]);
						if (t > ETA)
							{
							tst1 = t;
							tst2 = tst1 + 1. / tst1;
							if (tst2 <= tst1)
								{
								for (j = i; j <= en; j++)
									h[j][en] /= t;
								}
							}
						}
					}
				}
			}
		}
	
	for (i=0; i<dim; i++)
		{
		if ((i < low) || (i > high)) /* changed to rid goto statement */
			{
			for (j=i; j<dim; j++)
				z[i][j] = h[i][j];
			}
		}

	/* multiply by transformation matrix to give vectors of original
	   full matrix */
	for (j=dim-1; j>=low; j--)
		{
		m = MIN(j, high);
		for (i=low; i<=high; i++)
			{
			zz = 0.0;
			for (k = low; k <= m; k++)
				zz += z[i][k] * h[k][j];
			z[i][j] = zz;
			}
		}

	return (0);
	
#if 0
int hqr2 (int *nm, int *n, int *low, int *igh, MrBFlt *h__, MrBFlt *wr, MrBFlt *wi, MrBFlt *z__, int *ierr)
	
{

	/* system generated locals */
	int h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3;
	MrBFlt d__1, d__2, d__3, d__4;

	/* builtin functions */
	MrBFlt sqrt(doublereal), d_sign(doublereal *, doublereal *);

	/* Local variables */
	static MrBFlt norm;
	static int i__, j, k, l, m;
	static MrBFlt p, q, r__, s, t, w, x, y;
	static int na, ii, en, jj;
	static MrBFlt ra, sa;
	static int ll, mm, nn;
	static MrBFlt vi, vr, zz;
	static logical notlas;
	static int mp2, itn, its, enm2;
	static MrBFlt tst1, tst2;

	/* parameter adjustments */
	z_dim1 = *nm;
	z_offset = z_dim1 + 1;
	z__ -= z_offset;
	--wi;
	--wr;
	h_dim1 = *nm;
	h_offset = h_dim1 + 1;
	h__ -= h_offset;

	/* function Body */
	*ierr = 0;
	norm = 0.;
	k = 1;
	
	/* .......... store roots isolated by balanc and compute matrix norm .......... */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) 
		{
		i__2 = *n;
		for (j = k; j <= i__2; ++j) 
			{
			/* L40: */
			norm += (d__1 = h__[i__ + j * h_dim1], abs(d__1));
			}
		k = i__;
		if (i__ >= *low && i__ <= *igh) 
			goto L50;
		wr[i__] = h__[i__ + i__ * h_dim1];
		wi[i__] = 0.;
		L50:
			;
		}

	en = *igh;
	t = 0.;
	itn = *n * 30;
	
	/* ..........search for next eigenvalues.......... */
	L60:
	if (en < *low) 
		goto L340;
	its = 0;
	na = en - 1;
	enm2 = na - 1;
	
	/* ..........look for single small sub-diagonal element for l=en step -1 until low do -- .......... */
	L70:
	i__1 = en;
	for (ll = *low; ll <= i__1; ++ll) 
		{
		l = en + *low - ll;
		if (l == *low) 
			goto L100;
		s = (d__1 = h__[l - 1 + (l - 1) * h_dim1], abs(d__1)) + (d__2 = h__[l + l * h_dim1], abs(d__2));
		if (s == 0.0) 
			s = norm;
		tst1 = s;
		tst2 = tst1 + (d__1 = h__[l + (l - 1) * h_dim1], abs(d__1));
		if (tst2 == tst1) 
			goto L100;
		/* L80: */
		}
		
	/* .......... form shift .......... */
	L100:
	x = h__[en + en * h_dim1];
	if (l == en) 
		goto L270;
	y = h__[na + na * h_dim1];
	w = h__[en + na * h_dim1] * h__[na + en * h_dim1];
	if (l == na) 
		goto L280;
	if (itn == 0) 
		goto L1000;
	if (its != 10 && its != 20) 
		goto L130;

	/* .......... form exceptional shift .......... */
	t += x;

	i__1 = en;
	for (i__ = *low; i__ <= i__1; ++i__) 
		{
		/* L120: */
		h__[i__ + i__ * h_dim1] -= x;
		}

	s = (d__1 = h__[en + na * h_dim1], abs(d__1)) + (d__2 = h__[na + enm2 * h_dim1], abs(d__2));
	x = s * 0.75;
	y = x;
	w = s * -0.4375 * s;
	L130:
	++its;
	--itn;
	
	/* .......... look for two consecutive small sub-diagonal elements for m=en-2 step -1 until l do -- .......... */
	i__1 = enm2;
	for (mm = l; mm <= i__1; ++mm) 
		{
		m = enm2 + l - mm;
		zz = h__[m + m * h_dim1];
		r__ = x - zz;
		s = y - zz;
		p = (r__ * s - w) / h__[m + 1 + m * h_dim1] + h__[m + (m + 1) * h_dim1];
		q = h__[m + 1 + (m + 1) * h_dim1] - zz - r__ - s;
		r__ = h__[m + 2 + (m + 1) * h_dim1];
		s = abs(p) + abs(q) + abs(r__);
		p /= s;
		q /= s;
		r__ /= s;
		if (m == l) 
			goto L150;
		tst1 = abs(p) * ((d__1 = h__[m - 1 + (m - 1) * h_dim1], abs(d__1)) + 
		abs(zz) + (d__2 = h__[m + 1 + (m + 1) * h_dim1], abs(d__2)));
		tst2 = tst1 + (d__1 = h__[m + (m - 1) * h_dim1], abs(d__1)) * (abs(q) + abs(r__));
		if (tst2 == tst1) 
			goto L150;
		/* L140: */
		}
	L150:
	mp2 = m + 2;

	i__1 = en;
	for (i__ = mp2; i__ <= i__1; ++i__) 
		{
		h__[i__ + (i__ - 2) * h_dim1] = 0.0;
		if (i__ == mp2)
			goto L160;
		h__[i__ + (i__ - 3) * h_dim1] = 0.;
		L160:
			;
		}
		
	/*     .......... MrBFlt qr step involving rows l to en and columns m to en .......... */
	i__1 = na;
	for (k = m; k <= i__1; ++k) 
		{
		notlas = k != na;
		if (k == m) 
			goto L170;
		p = h__[k + (k - 1) * h_dim1];
		q = h__[k + 1 + (k - 1) * h_dim1];
		r__ = 0.;
		if (notlas) 
			r__ = h__[k + 2 + (k - 1) * h_dim1];
		x = abs(p) + abs(q) + abs(r__);
		if (x == 0.) 
			goto L260;
		p /= x;
		q /= x;
		r__ /= x;
		L170:
		d__1 = sqrt(p * p + q * q + r__ * r__);
		s = d_sign(&d__1, &p);
		if (k == m) 
			goto L180;
		h__[k + (k - 1) * h_dim1] = -s * x;
		goto L190;
		L180:
		if (l != m) 
			{
			h__[k + (k - 1) * h_dim1] = -h__[k + (k - 1) * h_dim1];
			}
		L190:
		p += s;
		x = p / s;
		y = q / s;
		zz = r__ / s;
		q /= p;
		r__ /= p;
		if (notlas) 
			goto L225;
		
		/* .......... row modification .......... */
		i__2 = *n;
		for (j = k; j <= i__2; ++j) 
			{
			p = h__[k + j * h_dim1] + q * h__[k + 1 + j * h_dim1];
			h__[k + j * h_dim1] -= p * x;
			h__[k + 1 + j * h_dim1] -= p * y;
			/* L200: */
			}

		/* computing MIN */
		i__2 = en, i__3 = k + 3;
		j = min(i__2,i__3);
		
		/* .......... column modification .......... */
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) 
			{
			p = x * h__[i__ + k * h_dim1] + y * h__[i__ + (k + 1) * h_dim1];
			h__[i__ + k * h_dim1] -= p;
			h__[i__ + (k + 1) * h_dim1] -= p * q;
			/* L210: */
			}
			
		/* .......... accumulate transformations .......... */
		i__2 = *igh;
		for (i__ = *low; i__ <= i__2; ++i__) 
			{
			p = x * z__[i__ + k * z_dim1] + y * z__[i__ + (k + 1) * z_dim1];
			z__[i__ + k * z_dim1] -= p;
			z__[i__ + (k + 1) * z_dim1] -= p * q;
			/* L220: */
			}
		goto L255;
		L225:
		
		/* .......... row modification .......... */
		i__2 = *n;
		for (j = k; j <= i__2; ++j) 
			{
			p = h__[k + j * h_dim1] + q * h__[k + 1 + j * h_dim1] + r__ * h__[k + 2 + j * h_dim1];
			h__[k + j * h_dim1] -= p * x;
			h__[k + 1 + j * h_dim1] -= p * y;
			h__[k + 2 + j * h_dim1] -= p * zz;
			/* L230: */
			}

		/* computing MIN */
		i__2 = en, i__3 = k + 3;
		j = min(i__2,i__3);
		
		/* .......... column modification .......... */
		i__2 = j;
		for (i__ = 1; i__ <= i__2; ++i__) 
			{
			p = x * h__[i__ + k * h_dim1] + y * h__[i__ + (k + 1) * h_dim1] + 
			zz * h__[i__ + (k + 2) * h_dim1];
			h__[i__ + k * h_dim1] -= p;
			h__[i__ + (k + 1) * h_dim1] -= p * q;
			h__[i__ + (k + 2) * h_dim1] -= p * r__;
			/* L240: */
			}
		
		/* .......... accumulate transformations .......... */
		i__2 = *igh;
		for (i__ = *low; i__ <= i__2; ++i__) 
			{
			p = x * z__[i__ + k * z_dim1] + y * z__[i__ + (k + 1) * z_dim1] + zz * z__[i__ + (k + 2) * z_dim1];
			z__[i__ + k * z_dim1] -= p;
			z__[i__ + (k + 1) * z_dim1] -= p * q;
			z__[i__ + (k + 2) * z_dim1] -= p * r__;
			/* L250: */
			}
		L255:
		L260:
			;
		}
	goto L70;
	
	/* .......... one root found .......... */
	L270:
	h__[en + en * h_dim1] = x + t;
	wr[en] = h__[en + en * h_dim1];
	wi[en] = 0.;
	en = na;
	goto L60;
	
	/* .......... two roots found .......... */
	L280:
	p = (y - x) / 2.;
	q = p * p + w;
	zz = sqrt((abs(q)));
	h__[en + en * h_dim1] = x + t;
	x = h__[en + en * h_dim1];
	h__[na + na * h_dim1] = y + t;
	if (q < 0.) 
		goto L320;
	
	/* .......... real pair .......... */
	zz = p + d_sign(&zz, &p);
	wr[na] = x + zz;
	wr[en] = wr[na];
	if (zz != 0.) 
		{
		wr[en] = x - w / zz;
		}
	wi[na] = 0.0;
	wi[en] = 0.0;
	x = h__[en + na * h_dim1];
	s = abs(x) + abs(zz);
	p = x / s;
	q = zz / s;
	r__ = sqrt(p * p + q * q);
	p /= r__;
	q /= r__;
	
	/* .......... row modification .......... */
	i__1 = *n;
	for (j = na; j <= i__1; ++j) 
		{
		zz = h__[na + j * h_dim1];
		h__[na + j * h_dim1] = q * zz + p * h__[en + j * h_dim1];
		h__[en + j * h_dim1] = q * h__[en + j * h_dim1] - p * zz;
		/* L290: */
		}
	
	/* .......... column modification .......... */
	i__1 = en;
	for (i__ = 1; i__ <= i__1; ++i__) 
		{
		zz = h__[i__ + na * h_dim1];
		h__[i__ + na * h_dim1] = q * zz + p * h__[i__ + en * h_dim1];
		h__[i__ + en * h_dim1] = q * h__[i__ + en * h_dim1] - p * zz;
		/* L300: */
		}
		
	/* .......... accumulate transformations .......... */
	i__1 = *igh;
	for (i__ = *low; i__ <= i__1; ++i__) 
		{
		zz = z__[i__ + na * z_dim1];
		z__[i__ + na * z_dim1] = q * zz + p * z__[i__ + en * z_dim1];
		z__[i__ + en * z_dim1] = q * z__[i__ + en * z_dim1] - p * zz;
		/* L310: */
		}
	goto L330;
	
	/* .......... complex pair .......... */
	L320:
	wr[na] = x + p;
	wr[en] = x + p;
	wi[na] = zz;
	wi[en] = -zz;
	L330:
	en = enm2;
	goto L60;
	
	/* .......... all roots found.  backsubstitute to find vectors of upper triangular form .......... */
	L340:
	if (norm == 0.0) 
		goto L1001;

	/* .......... for en=n step -1 until 1 do -- .......... */
	i__1 = *n;
	for (nn = 1; nn <= i__1; ++nn) 
		{
		en = *n + 1 - nn;
		p = wr[en];
		q = wi[en];
		na = en - 1;
		if (q < 0.) 
			goto L710;
		else if (q == 0) 
			goto L600;
		else 
			goto L800;
			
		/* .......... real vector .......... */
		L600:
		m = en;
		h__[en + en * h_dim1] = 1.0;
		if (na == 0) 
			goto L800;
		
		/*     .......... for i=en-1 step -1 until 1 do -- .......... */
		i__2 = na;
		for (ii = 1; ii <= i__2; ++ii) 
			{
			i__ = en - ii;
			w = h__[i__ + i__ * h_dim1] - p;
			r__ = 0.0;

			i__3 = en;
			for (j = m; j <= i__3; ++j) 
				{
				/* L610: */
				r__ += h__[i__ + j * h_dim1] * h__[j + en * h_dim1];
				}

			if (wi[i__] >= 0.0) 
				goto L630;
			zz = w;
			s = r__;
			goto L700;
			L630:
			m = i__;
			if (wi[i__] != 0.0) 
				goto L640;
			t = w;
			if (t != 0.0)
				goto L635;
			tst1 = norm;
			t = tst1;
			L632:
			t *= 0.01;
			tst2 = norm + t;
			if (tst2 > tst1) 
				goto L632;
			L635:
			h__[i__ + en * h_dim1] = -r__ / t;
			goto L680;
			
			/* .......... solve real equations .......... */
			L640:
			x = h__[i__ + (i__ + 1) * h_dim1];
			y = h__[i__ + 1 + i__ * h_dim1];
			q = (wr[i__] - p) * (wr[i__] - p) + wi[i__] * wi[i__];
			t = (x * s - zz * r__) / q;
			h__[i__ + en * h_dim1] = t;
			if (abs(x) <= abs(zz)) 
				goto L650;
			h__[i__ + 1 + en * h_dim1] = (-r__ - w * t) / x;
			goto L680;
			L650:
			h__[i__ + 1 + en * h_dim1] = (-s - y * t) / zz;

			/*     .......... overflow control .......... */
			L680:
			t = (d__1 = h__[i__ + en * h_dim1], abs(d__1));
			if (t == 0.0) 
				goto L700;
			tst1 = t;
			tst2 = tst1 + 1.0 / tst1;
			if (tst2 > tst1) 
				goto L700;
			i__3 = en;
			for (j = i__; j <= i__3; ++j) 
				{
				h__[j + en * h_dim1] /= t;
				/* L690: */
				}

			L700:
				;
			}
			
		/* .......... end real vector .......... */
		goto L800;
		
		/* .......... complex vector .......... */
		L710:
		m = na;
		
		/* .......... last vector component chosen imaginary so that eigenvector matrix is triangular .......... */
		if ((d__1 = h__[en + na * h_dim1], abs(d__1)) <= (d__2 = h__[na + en *
		h_dim1], abs(d__2))) 
			goto L720;
		h__[na + na * h_dim1] = q / h__[en + na * h_dim1];
		h__[na + en * h_dim1] = -(h__[en + en * h_dim1] - p) / h__[en + na * h_dim1];
		goto L730;
		L720:
		d__1 = -h__[na + en * h_dim1];
		d__2 = h__[na + na * h_dim1] - p;
		cdiv_(&c_b49, &d__1, &d__2, &q, &h__[na + na * h_dim1], &h__[na + en *
		h_dim1]);
		L730:
		h__[en + na * h_dim1] = 0.0;
		h__[en + en * h_dim1] = 1.0;
		enm2 = na - 1;
		if (enm2 == 0) 
			goto L800;

		/*     .......... for i=en-2 step -1 until 1 do -- .......... */
		i__2 = enm2;
		for (ii = 1; ii <= i__2; ++ii) 
			{
			i__ = na - ii;
			w = h__[i__ + i__ * h_dim1] - p;
			ra = 0.0;
			sa = 0.0;

			i__3 = en;
			for (j = m; j <= i__3; ++j) 
				{
				ra += h__[i__ + j * h_dim1] * h__[j + na * h_dim1];
				sa += h__[i__ + j * h_dim1] * h__[j + en * h_dim1];
				/* L760: */
				}

			if (wi[i__] >= 0.0) 
				goto L770;
			zz = w;
			r__ = ra;
			s = sa;
			goto L795;
			L770:
			m = i__;
			if (wi[i__] != 0.0) 
				goto L780;
			d__1 = -ra;
			d__2 = -sa;
			cdiv_(&d__1, &d__2, &w, &q, &h__[i__ + na * h_dim1], &h__[i__ + en * h_dim1]);
			goto L790;
			
			/*     .......... solve complex equations .......... */
			L780:
			x = h__[i__ + (i__ + 1) * h_dim1];
			y = h__[i__ + 1 + i__ * h_dim1];
			vr = (wr[i__] - p) * (wr[i__] - p) + wi[i__] * wi[i__] - q * q;
			vi = (wr[i__] - p) * 2.0 * q;
			if (vr != 0.0 || vi != 0.0) 
				goto L784;
			tst1 = norm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(zz));
			vr = tst1;
			L783:
			vr *= 0.01;
			tst2 = tst1 + vr;
			if (tst2 > tst1) 
				goto L783;
			L784:
			d__1 = x * r__ - zz * ra + q * sa;
			d__2 = x * s - zz * sa - q * ra;
			cdiv_(&d__1, &d__2, &vr, &vi, &h__[i__ + na * h_dim1], &h__[i__ + en * h_dim1]);
			if (abs(x) <= abs(zz) + abs(q)) 
				goto L785;
			h__[i__ + 1 + na * h_dim1] = (-ra - w * h__[i__ + na * h_dim1] + q * h__[i__ + en * h_dim1]) / x;
			h__[i__ + 1 + en * h_dim1] = (-sa - w * h__[i__ + en * h_dim1] - q * h__[i__ + na * h_dim1]) / x;
			goto L790;
			L785:
			d__1 = -r__ - y * h__[i__ + na * h_dim1];
			d__2 = -s - y * h__[i__ + en * h_dim1];
			cdiv_(&d__1, &d__2, &zz, &q, &h__[i__ + 1 + na * h_dim1], &h__[i__ + 1 + en * h_dim1]);

			/*     .......... overflow control .......... */
			L790:
			/* Computing MAX */
			d__3 = (d__1 = h__[i__ + na * h_dim1], abs(d__1)), d__4 = (d__2 = h__[i__ + en * h_dim1], abs(d__2));
			t = max(d__3,d__4);
			if (t == 0.0) 
				goto L795;
			tst1 = t;
			tst2 = tst1 + 1.0 / tst1;
			if (tst2 > tst1) 
				goto L795;
			i__3 = en;
			for (j = i__; j <= i__3; ++j) 
				{
				h__[j + na * h_dim1] /= t;
				h__[j + en * h_dim1] /= t;
				/* L792: */
				}
			L795:
				;
			}
		/*     .......... end complex vector .......... */
		L800:
			;
		}
	/*     .......... end back substitution vectors of isolated roots .......... */
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) 
		{
		if (i__ >= *low && i__ <= *igh) 
			goto L840;
		i__2 = *n;
		for (j = i__; j <= i__2; ++j) 
			{
			/* L820: */
			z__[i__ + j * z_dim1] = h__[i__ + j * h_dim1];
			}
		L840:
		;
		}
		
	/* .......... multiply by transformation matrix to give vectors of original full matrix. */
	/*            for j=n step -1 until low do -- .......... */
	i__1 = *n;
	for (jj = *low; jj <= i__1; ++jj) 
		{
		j = *n + *low - jj;
		m = min(j,*igh);

		i__2 = *igh;
		for (i__ = *low; i__ <= i__2; ++i__) 
			{
			zz = 0.0;
			i__3 = m;
			for (k = *low; k <= i__3; ++k) 
				{
				/* L860: */
				zz += z__[i__ + k * z_dim1] * h__[k + j * h_dim1];
				}

			z__[i__ + j * z_dim1] = zz;
			/* L880: */
			}
		}

	goto L1001;
	/* .......... set error -- all eigenvalues have not converged after 30*n iterations .......... */
	L1000:
	*ierr = en;
	L1001:
	return 0;
	
}
/* end f2c version of code */
#endif

}





MrBFlt IncompleteBetaFunction (MrBFlt alpha, MrBFlt beta, MrBFlt x)

{

	MrBFlt		bt, gm1, gm2, gm3, temp;
	
	if (x < 0.0 || x > 1.0) 
		{
		MrBayesPrint ("%s   Error: Problem in IncompleteBetaFunction.\n", spacer);
		exit (0);
		}
	if (fabs(x) < ETA || fabs(x-1.0)<ETA) /* x == 0.0 || x == 1.0 */
		{
		bt = 0.0;
		}
	else
		{
		gm1 = LnGamma (alpha + beta);
		gm2 = LnGamma (alpha);
		gm3 = LnGamma (beta);
		temp = gm1 - gm2 - gm3 + (alpha) * log(x) + (beta) * log(1.0 - x);
		bt = exp(temp);
		}
	if (x < (alpha + 1.0)/(alpha + beta + 2.0))
		return (bt * BetaCf(alpha, beta, x) / alpha);
	else
		return (1.0 - bt * BetaCf(beta, alpha, 1.0-x) / beta);

}





/*---------------------------------------------------------------------------------
|
|   IncompleteGamma
|
|   Returns the incomplete gamma ratio I(x,alpha) where x is the upper
|   limit of the integration and alpha is the shape parameter.  Returns (-1)
|   if in error.   
|
|   Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied
|      Statistics, 19:285-287 (AS32)
|
---------------------------------------------------------------------------------*/
MrBFlt IncompleteGamma (MrBFlt x, MrBFlt alpha, MrBFlt LnGamma_alpha)

{

	int 			i;
	MrBFlt 		p = alpha, g = LnGamma_alpha,
					accurate = 1e-8, overflow = 1e30,
					factor, gin = 0.0, rn = 0.0, a = 0.0, b = 0.0, an = 0.0, 
					dif = 0.0, term = 0.0, pn[6];

	if (fabs(x) < ETA) 
		return (0.0);
	if (x < 0 || p <= 0) 
		return (-1.0);

	factor = exp(p*log(x)-x-g);   
	if (x>1 && x>=p) 
		goto l30;
	gin = 1.0;  
	term = 1.0;  
	rn = p;
	l20:
		rn++;
		term *= x/rn;   
		gin += term;
		if (term > accurate) 
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = 1.0-p;   
		b = a+x+1.0;  
		term = 0.0;
		pn[0] = 1.0;  
		pn[1] = x;  
		pn[2] = x+1;  
		pn[3] = x*b;
		gin = pn[2]/pn[3];
	l32:
		a++;  
		b += 2.0;  
		term++;   
		an = a*term;
		for (i=0; i<2; i++) 
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (fabs(pn[5]) < ETA) 
			goto l35;
		rn = pn[4]/pn[5];   
		dif = fabs(gin-rn);
		if (dif>accurate) 
			goto l34;
		if (dif<=accurate*rn) 
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++) 
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow) 
			goto l32;
		for (i=0; i<4; i++) 
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = 1.0-factor*gin;
	l50:
		return (gin);

}





/*---------------------------------------------------------------------------------
|
|   InvertMatrix
|
|   Calculates aInv = a^{-1} using LU-decomposition. The input matrix a is
|   destroyed in the process. The program returns an error if the matrix is
|   singular. col and indx are work vectors.
|      
---------------------------------------------------------------------------------*/
int InvertMatrix (int dim, MrBFlt **a, MrBFlt *col, int *indx, MrBFlt **aInv)

{

	int			rc, i, j;
	
	rc = LUDecompose (dim, a, col, indx, (MrBFlt *)NULL);
	if (rc == FALSE)
		{
		for (j = 0; j < dim; j++)
			{
			for (i = 0; i < dim; i++)
				col[i] = 0.0;
			col[j] = 1.0;
			LUBackSubstitution (dim, a, indx, col);
			for (i = 0; i < dim; i++)
				aInv[i][j] = col[i];
			}
		}
		
	return (rc);
	
}





/*---------------------------------------------------------------------------------
|
|   LBinormal
|
|   L(h1,h2,r) = prob(x>h1, y>h2), where x and y are standard binormal, 
|   with r=corr(x,y),  error < 2e-7.
|
|   Drezner Z., and G.O. Wesolowsky (1990) On the computation of the
|      bivariate normal integral.  J. Statist. Comput. Simul. 35:101-107.
|
---------------------------------------------------------------------------------*/
MrBFlt LBinormal (MrBFlt h1, MrBFlt h2, MrBFlt r)

{
      
	int i;
	MrBFlt 		x[]={0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992};
	MrBFlt 		w[]={0.018854042, 0.038088059, 0.0452707394,0.038088059,0.018854042};
	MrBFlt 		Lh=0.0, r1, r2, r3, rr, aa, ab, h3, h5, h6, h7, h12, temp1, temp2, exp1, exp2;

	h12 = (h1 * h1 + h2 * h2) / 2.0;
	if (fabs(r) >= 0.7) 
		{
		r2 = 1.0 - r * r;   
		r3 = sqrt(r2);
		if (r < 0) 
			h2 *= -1;
		h3 = h1 * h2;   
		h7 = exp(-h3 / 2.0);
		if (fabs(r-1.0)>ETA)  /* fabs(r) != 1.0 */
			{
			h6 = fabs(h1-h2);  
			h5 = h6 * h6 / 2.0; 
			h6 /= r3; 
			aa = 0.5 - h3 / 8;  
			ab = 3.0 - 2.0 * aa * h5;
			temp1 = -h5 / r2;
			if (temp1 < -100.0)
				exp1 = 0.0;
			else
				exp1 = exp(temp1);
			Lh = 0.13298076 * h6 * ab * (1.0 - CdfNormal(h6)) - exp1 * (ab + aa * r2) * 0.053051647;
			for (i=0; i<5; i++) 
				{
				r1 = r3 * x[i];
				rr = r1 * r1;   
				r2 = sqrt(1.0 - rr);
				temp1 = -h5 / rr;
				if (temp1 < -100.0)
					exp1 = 0.0;
				else
					exp1 = exp(temp1);
				temp2 = -h3 / (1.0 + r2);
				if (temp2 < -100.0)
					exp2 = 0.0;
				else
					exp2 = exp(temp2);
				Lh -= w[i] * exp1 * (exp2 / r2 / h7 - 1.0 - aa * rr);
				}
			}
		if (r > 0) 
			Lh = Lh * r3 * h7 + (1.0 - CdfNormal(MAX(h1, h2)));
		else if (r<0) 
			Lh = (h1 < h2 ? CdfNormal(h2) - CdfNormal(h1) : 0) - Lh * r3 * h7;
		}
	else 
		{
		h3 = h1 * h2;
		if (fabs(r)>ETA) 
			{
			for (i=0; i<5; i++) 
				{
				r1 = r * x[i]; 
				r2 = 1.0 - r1 * r1;
				temp1 = (r1 * h3 - h12) / r2;
				if (temp1 < -100.0)
					exp1 = 0.0;
				else
					exp1 = exp(temp1);
				Lh += w[i] * exp1 / sqrt(r2);
				}
			}
		Lh = (1.0 - CdfNormal(h1)) * (1.0 - CdfNormal(h2)) + r * Lh;
		}
	return (Lh);
   
}    





/*---------------------------------------------------------------------------------
|
|   LnFactorial: Calculates the log of the factorial for an integer
|
---------------------------------------------------------------------------------*/
MrBFlt  LnFactorial (int value)
{
	int		i;
	MrBFlt  result;

	result = 0.0;

	for (i = 2; i<=value; i++)
		result += log(i);

	return result;
}





/*---------------------------------------------------------------------------------
|
|   LnGamma
|
|   Calculates the log of the gamma function. The Gamma function is equal
|   to:
|
|      Gamma(alp) = {integral from 0 to infinity} t^{alp-1} e^-t dt
|
|   The result is accurate to 10 decimal places. Stirling's formula is used
|   for the central polynomial part of the procedure.
|
|   Pike, M. C. and I. D. Hill.  1966.  Algorithm 291: Logarithm of the gamma
|      function.  Communications of the Association for Computing
|      Machinery, 9:684.
|      
---------------------------------------------------------------------------------*/
MrBFlt LnGamma (MrBFlt alp)

{

	MrBFlt 		x = alp, f = 0.0, z;
	
	if (x < 7) 
		{
		f = 1.0;
		z = x-1.0;
		while (++z < 7.0)  
			f *= z;
		x = z;   
		f = -log(f);
		}
	z = 1.0 / (x*x);
	return  (f + (x-0.5)*log(x) - x + 0.918938533204673 + 
			(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
			0.083333333333333)/x);  

}





/* Log probability for a value drawn from a gamma distribution */
MrBFlt LnProbGamma (MrBFlt alpha, MrBFlt beta, MrBFlt x)

{
    MrBFlt lnProb;

    lnProb = (alpha-1.0)*log(x) + alpha*log(beta) - x*beta - LnGamma(alpha);

    return lnProb;
}





/* Log probability for a value drawn from a lognormal distribution */
MrBFlt LnProbLogNormal (MrBFlt exp, MrBFlt var, MrBFlt x)

{
	MrBFlt		mu, sqdev, lnProb;

	mu = log (exp) - (var * var / 2.0);
	sqdev = log(x) - mu;
	sqdev *= sqdev;

	lnProb = - log (x * var * sqrt (2.0 * PI)) - (sqdev / (2.0*var*var));

	return lnProb;
}





/* Log probability for a value drawn from a scaled gamma distribution */
MrBFlt LnProbScaledGamma (MrBFlt alpha, MrBFlt x)

{
    MrBFlt lnProb;

    lnProb = (alpha - 1.0) * log(x) - LnGamma(alpha) + alpha*log(alpha) - x*alpha;

    return lnProb;
}





/* Log ratio for two values drawn from a lognormal distribution */
MrBFlt LnRatioLogNormal (MrBFlt exp, MrBFlt var, MrBFlt xNew, MrBFlt xOld)

{
	MrBFlt		mu, sqdev, sqdevNew, lnRatio;

	mu = log (exp) - (var * var / 2.0);
	sqdev = log (xOld) - mu;
	sqdev *= sqdev;
	sqdevNew = log (xNew) - mu;
	sqdevNew *= sqdevNew;

	lnRatio = log (xOld) - log (xNew);
	lnRatio += (sqdev - sqdevNew) / (2.0 * var * var);

	return lnRatio;
}





/*---------------------------------------------------------------------------------
|
|   LogBase2Plus1
|
|   This function is called from ComputeMatrixExponential.
|      
---------------------------------------------------------------------------------*/
int LogBase2Plus1 (MrBFlt x)

{

	int		j = 0;

	while(x > 1.0 - 1.0e-07) 
		{
		x /= 2.0;
		j++;
		}
		
	return (j);

}





/*---------------------------------------------------------------------------------
|
|   LogNormalRandomVariable
|
|   Draw a random variable from a lognormal distribution.
|      
---------------------------------------------------------------------------------*/
MrBFlt LogNormalRandomVariable (MrBFlt mean, MrBFlt var, safeLong *seed)

{

	MrBFlt      x;
    
    x = PointNormal(RandomNumber(seed));
    
    x*= var;
    x += log(mean);
    
    return exp(x);
}





/*---------------------------------------------------------------------------------
|
|   LUBackSubstitution
|
|   Back substitute into an LU-decomposed matrix.
|      
---------------------------------------------------------------------------------*/
void LUBackSubstitution (int dim, MrBFlt **a, int *indx, MrBFlt *b)

{

	int			i, ip, j, ii = -1;
	MrBFlt		sum;

	for (i=0; i<dim; i++)
		{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0)
			{
			for (j=ii; j<=i-1; j++)
				sum -= a[i][j] * b[j];
			}
		else if (fabs(sum)>ETA)
			ii = i;
		b[i] = sum;
		}
	for (i=dim-1; i>=0; i--)
		{
		sum = b[i];
		for (j=i+1; j<dim; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
		}
		
}





/*---------------------------------------------------------------------------------
|
|   LUDecompose
|
|   Calculate the LU-decomposition of the matrix a. The matrix a is replaced.
|      
---------------------------------------------------------------------------------*/
int LUDecompose (int dim, MrBFlt **a, MrBFlt *vv, int *indx, MrBFlt *pd)

{

	int			i, imax=0, j, k;
	MrBFlt		big, dum, sum, temp, d;

	d = 1.0;
	for (i=0; i<dim; i++)
		{
		big = 0.0;
		for (j = 0; j < dim; j++)
			{
			if ((temp = fabs(a[i][j])) > big)
				big = temp;
			}
		if (fabs(big)<ETA)
			{
			MrBayesPrint ("%s   Error: Problem in LUDecompose\n", spacer);
			return (ERROR);
			}
		vv[i] = 1.0 / big;
		}
	for (j=0; j<dim; j++)
		{
		for (i = 0; i < j; i++)
			{
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			}
		big = 0.0;
		for (i=j; i<dim; i++)
			{
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			dum = vv[i] * fabs(sum);
			if (dum >= big)
				{
				big = dum;
				imax = i;
				}
			}
		if (j != imax)
			{
			for (k=0; k<dim; k++)
				{
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
				}	
			d = -d;
			vv[imax] = vv[j];
			}
		indx[j] = imax;
		if (fabs(a[j][j])<ETA)
			a[j][j] = TINY;
		if (j != dim - 1)
			{
			dum = 1.0 / (a[j][j]);
			for (i=j+1; i<dim; i++)
				a[i][j] *= dum;
			}
		}
	if (pd != NULL)
		*pd = d;
		
	return (NO_ERROR);
	
}





/*---------------------------------------------------------------------------------
|
|   MultiplyMatrices
|
|   Multiply matrix a by matrix b and put the results in matrix result.
|
---------------------------------------------------------------------------------*/
void MultiplyMatrices (int dim, MrBFlt **a, MrBFlt **b, MrBFlt **result)

{

	register int	i, j, k;
	MrBFlt			**temp;

	temp = AllocateSquareDoubleMatrix (dim);

	for (i=0; i<dim; i++)
		{
		for (j=0; j<dim; j++) 
			{
			temp[i][j] = 0.0;
			for (k=0; k<dim; k++) 
				{
				temp[i][j] += a[i][k] * b[k][j];
				}
			}
		}
	for (i=0; i<dim; i++)
		{
		for (j=0; j<dim; j++) 
			{
			result[i][j] = temp[i][j];
			}
		}
		
	FreeSquareDoubleMatrix (temp);

}





/*---------------------------------------------------------------------------------
|
|   MultiplyMatrixByScalar
|
|   Multiply the elements of matrix a by a scalar.
|
---------------------------------------------------------------------------------*/
void MultiplyMatrixByScalar (int dim, MrBFlt **a, MrBFlt scalar, MrBFlt **result)

{

	int			row, col;

	for (row=0; row<dim; row++)
		for (col=0; col<dim; col++)
     		 result[row][col] = a[row][col] * scalar;

}





/*---------------------------------------------------------------------------------
|
|   MultiplyMatrixNTimes
|
---------------------------------------------------------------------------------*/
int MultiplyMatrixNTimes (int dim, MrBFlt **Mat, int power, MrBFlt **Result)

{

	register int	i, j;
	int				k, numSquares, numRemaining;
	MrBFlt			**TempIn, **TempOut;

	if (power < 0)
		{
		MrBayesPrint ("%s   Error: Power cannot be a negative number.\n", spacer);
		return (ERROR);
		}
	else if (power == 0)
		{
		for (i=0; i<dim; i++)
			for (j=0; j<dim; j++)
				Result[i][j] = 1.0;
		}
	else
		{
		TempIn  = AllocateSquareDoubleMatrix (dim);
		TempOut = AllocateSquareDoubleMatrix (dim);

		/* how many times can I multiply the matrices together */
		numSquares = 0;
		while (pow (2.0, (MrBFlt)(numSquares)) < power)
			numSquares++;
		numRemaining = power - (int)(pow(2.0, (MrBFlt)(numSquares)));
		
		/* now, multiply matrix by power of 2's */
		CopyDoubleMatrices (dim, Mat, TempIn);
		for (k=0; k<numSquares; k++)
			{
			MultiplyMatrices (dim, TempIn, TempIn, TempOut);
			CopyDoubleMatrices (dim, TempOut, TempIn);
			}
			
		/* TempIn is Mat^numSquares. Now, multiply it by Mat numRemaining times */
		for (k=0; k<numSquares; k++)
			{
			MultiplyMatrices (dim, TempIn, Mat, TempOut);
			CopyDoubleMatrices (dim, TempOut, TempIn);
			}
			
		/* copy result */
		CopyDoubleMatrices (dim, TempIn, Result);
		
		FreeSquareDoubleMatrix (TempIn);
		FreeSquareDoubleMatrix (TempOut);
		}

	return (NO_ERROR);
	
}





/*---------------------------------------------------------------------------------
|
|   PointChi2
|
|   Returns z so that Prob(x < z) = prob where x is Chi2 distributed with df=v. 
|   Returns -1 if in error.   0.000002 < prob < 0.999998.
|
---------------------------------------------------------------------------------*/
MrBFlt PointChi2 (MrBFlt prob, MrBFlt v)

{

	MrBFlt 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0, 
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0) 
		return (-1.0);
	g = LnGamma (v/2.0);
	xx = v/2.0;   
	c = xx - 1.0;
	if (v >= -1.24*log(p)) 
		goto l1;
	ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (ch-e<0) 
		return (ch);
	goto l4;
	l1:
		if (v > 0.32) 
			goto l3;
		ch = 0.4;   
		a = log(1.0-p);
	l2:
		q = ch;  
		p1 = 1.0+ch*(4.67+ch);  
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0) 
			goto l4;
		else                       
			goto l2;
	l3: 
		x = PointNormal (p);
		p1 = 0.222222/v;   
		ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)  
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;   
		p1 = 0.5*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < 0.0) 
			{
			MrBayesPrint ("%s   Error: Problem in PointChi2", spacer);
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));   
		b = t/ch;  
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a)) / 2520.0;
		s6 = (120.0+c*(346.0+127.0*c)) / 5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e) 
			goto l4;
		return (ch);

}





/*---------------------------------------------------------------------------------
|
|   PointNormal
|
|   Returns z so That Prob{x<z} = prob where x ~ N(0,1) and
|   (1e-12) < prob < 1-(1e-12).  Returns (-9999) if in error.
|
|   Odeh, R. E. and J. O. Evans.  1974.  The percentage points of the normal
|     distribution.  Applied Statistics, 22:96-97 (AS70).
|
|   Newer methods:
|
|   Wichura, M. J.  1988.  Algorithm AS 241: The percentage points of the
|      normal distribution.  37:477-484.
|   Beasley, JD & S. G. Springer.  1977.  Algorithm AS 111: The percentage
|      points of the normal distribution.  26:118-121.
|
---------------------------------------------------------------------------------*/
MrBFlt PointNormal (MrBFlt prob)

{

	MrBFlt 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
 					a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
 					b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
 					y, z = 0, p = prob, p1;

	p1 = (p<0.5 ? p : 1-p);
	if (p1<1e-20) 
	   return (-9999);
	y = sqrt (log(1/(p1*p1)));   
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	
	return (p<0.5 ? -z : z);

}





/*---------------------------------------------------------------------------------
|
|   PrintComplexVector
|
|   Prints a vector of dim complex numbers.
|
---------------------------------------------------------------------------------*/
void PrintComplexVector (int dim, complex *vec)

{

	int     i;

	MrBayesPrint ("{");
	for (i = 0; i < (dim - 1); i++) 
		{
		MrBayesPrint ("%lf + %lfi, ", vec[i].re, vec[i].im);
		if(i == 1) 
			MrBayesPrint("\n    ");
		}
	MrBayesPrint ("%lf + %lfi}\n", vec[dim - 1].re, vec[dim - 1].im);
    
}





/*---------------------------------------------------------------------------------
|
|   PrintSquareComplexMatrix
|
|   Prints a square matrix of complex numbers.
|
---------------------------------------------------------------------------------*/
void PrintSquareComplexMatrix (int dim, complex **m)

{

	int     row, col;

	MrBayesPrint ("{");
	for (row = 0; row < (dim - 1); row++) 
		{
		MrBayesPrint ("{");
		for(col = 0; col < (dim - 1); col++) 
			{
			MrBayesPrint ("%lf + %lfi, ", m[row][col].re, m[row][col].im);
			if(col == 1) 
				MrBayesPrint ("\n    ");
			}
		MrBayesPrint ("%lf + %lfi},\n", 
		m[row][dim - 1].re, m[row][dim - 1].im);
		}
	MrBayesPrint ("{");
	for (col = 0; col < (dim - 1); col++) 
		{
		MrBayesPrint ("%lf + %lfi, ", m[dim - 1][col].re, m[dim - 1][col].im);
		if(col == 1) 
			MrBayesPrint ("\n    ");
		}
	MrBayesPrint ("%lf + %lfi}}", m[dim - 1][dim - 1].re, m[dim - 1][dim - 1].im);
	MrBayesPrint ("\n");
	
}





/*---------------------------------------------------------------------------------
|
|   PrintSquareDoubleMatrix
|
|   Prints a square matrix of doubles.
|
---------------------------------------------------------------------------------*/
void PrintSquareDoubleMatrix (int dim, MrBFlt **matrix)

{

	int			i, j;
	
	for (i=0; i<dim; i++) 
		{
		for(j=0; j<dim; j++)
			MrBayesPrint ("%1.6lf ", matrix[i][j]);
		MrBayesPrint ("\n");
		}

}





/*---------------------------------------------------------------------------------
|
|   PrintSquareIntegerMatrix
|
|   Prints a square matrix of integers.
|
---------------------------------------------------------------------------------*/
void PrintSquareIntegerMatrix (int dim, int **matrix)

{

	int			i, j;
	
	for (i=0; i<dim; i++) 
		{
		for(j=0; j<dim; j++)
			MrBayesPrint ("%d ", matrix[i][j]);
		MrBayesPrint ("\n");
		}

}





/*---------------------------------------------------------------------------------
|
|   ProductOfRealAndComplex
|
|   Returns the complex product of a real and complex number.
|
---------------------------------------------------------------------------------*/
complex ProductOfRealAndComplex (MrBFlt a, complex b)

{

    complex 	c;
    
    c.re = a * b.re;
    c.im = a * b.im;
    
    return (c);
    
}





/*---------------------------------------------------------------------------------
|
|   PsiExp: Returns psi (also called digamma) exponentiated
|		Algorithm from http://lib.stat.cmu.edu/apstat/103
|
---------------------------------------------------------------------------------*/
MrBFlt  PsiExp (MrBFlt alpha)

{
	MrBFlt		digamma, y, r, s, c, s3, s4, s5, d1;
	
	s = 1.0e-05;
	c = 8.5;
	s3 = 8.333333333333333333333333e-02;
	s4 = 8.333333333333333333333333e-03;
	s5 = 3.968253968e-03;
	d1 = -0.577215664901532860606512;	/* negative of Euler's constant */

	digamma = 0.0;
	y = alpha;
	if (y <= 0.0)
		return (0.0);
	
	if (y <= s)
		{
		digamma = d1 - 1.0 / y;
		return (exp (digamma));
		}
	
	while (y < c)
		{
		digamma -= 1.0 / y;
		y += 1.0;
		}

	r = 1.0 / y;
	digamma += (log (y) - 0.5 * r);
	r *= r;
	digamma -= r * (s3 - r * (s4 - r * s5));
	
	return (exp (digamma));

}





/*---------------------------------------------------------------------------------
|
|   PsiGammaLnProb: Calculates the log probability of a PsiGamma distributed
|      variable
|
---------------------------------------------------------------------------------*/
MrBFlt  PsiGammaLnProb (MrBFlt alpha, MrBFlt value)
{
	MrBFlt	beta, lnProb;

	beta = PsiExp (alpha);

	lnProb = alpha * log (beta) - LnGamma (alpha) + (alpha - 1.0) * log (value) - beta * value;

	return lnProb;
}





/*---------------------------------------------------------------------------------
|
|   PsiGammaLnRatio: Calculates the log prob ratio of two PsiGamma distributed
|      variables
|
---------------------------------------------------------------------------------*/
MrBFlt  PsiGammaLnRatio (MrBFlt alpha, MrBFlt numerator, MrBFlt denominator)

{
	MrBFlt beta, lnRatio;

	beta = PsiExp (alpha);

	lnRatio = (alpha - 1.0) * (log (numerator) - log (denominator)) - beta * (numerator - denominator);
	
	return (lnRatio);
}





/*---------------------------------------------------------------------------------
|
|   PsiGammaRandomVariable: Returns a random draw from the PsiGamma
|
---------------------------------------------------------------------------------*/
MrBFlt  PsiGammaRandomVariable (MrBFlt alpha, safeLong *seed)
{
	return GammaRandomVariable (alpha, PsiExp(alpha), seed);
}





/*---------------------------------------------------------------------------------
|
|   QuantileGamma
|
---------------------------------------------------------------------------------*/
MrBFlt QuantileGamma (MrBFlt x, MrBFlt alfa, MrBFlt beta)

{

	MrBFlt		lnga1, quantile;

	lnga1 = LnGamma(alfa + 1.0);
	quantile = POINTGAMMA(x, alfa, beta);
	
	return (quantile);

}





/*---------------------------------------------------------------------------------
|
|   RandomNumber
|
|   This pseudorandom number generator is described in:
|   Park, S. K. and K. W. Miller.  1988.  Random number generators: good
|      ones are hard to find.  Communications of the ACM, 31(10):1192-1201.
|
---------------------------------------------------------------------------------*/
MrBFlt RandomNumber (safeLong *seed)

{

	safeLong	lo, hi, test;
	
	hi = (*seed) / 127773;
	lo = (*seed) % 127773;
	test = 16807 * lo - 2836 * hi;
	if (test > 0)
		*seed = test;
	else
		*seed = test + 2147483647;
	return ((MrBFlt)(*seed) / (MrBFlt)2147483647);

}





/*---------------------------------------------------------------------------------
|
|   RndGamma
|
---------------------------------------------------------------------------------*/
MrBFlt RndGamma (MrBFlt s, safeLong *seed)

{

	MrBFlt  r=0.0;
	
	if (s <= 0.0)      
		puts ("Gamma parameter less than zero\n");
	else if (s < 1.0)  
		r = RndGamma1 (s, seed);
	else if (s > 1.0)  
		r = RndGamma2 (s, seed);
	else    /* 0-log() == -1 * log(), but =- looks confusing */
		r -= log(RandomNumber(seed));
		
	return (r);
   
}





/*---------------------------------------------------------------------------------
|
|   RndGamma1
|
---------------------------------------------------------------------------------*/
MrBFlt RndGamma1 (MrBFlt s, safeLong *seed)

{

	MrBFlt			r, x=0.0, small=1e-37, w;
	static MrBFlt   a, p, uf, ss=10.0, d;
	
	if (fabs(s-ss)>ETA) /* s != ss */ 
		{
		a  = 1.0 - s;
		p  = a / (a + s * exp(-a));
		uf = p * pow(small / a, s);
		d  = a * log(a);
		ss = s;
		}
	for (;;) 
		{
		r = RandomNumber (seed);
		if (r > p)        
			x = a - log((1.0 - r) / (1.0 - p)), w = a * log(x) - d;
		else if (r>uf)  
			x = a * pow(r / p, 1.0 / s), w = x;
		else            
			return (0.0);
		r = RandomNumber (seed);
		if (1.0 - r <= w && r > 0.0)
		if (r*(w + 1.0) >= 1.0 || -log(r) <= w)  
			continue;
		break;
		}
		
	return (x);
   
}





/*---------------------------------------------------------------------------------
|
|   RndGamma2
|
---------------------------------------------------------------------------------*/
MrBFlt RndGamma2 (MrBFlt s, safeLong *seed)

{

	MrBFlt			r , d, f, g, x;
	static MrBFlt	b, h, ss=0.0;
	
	if (fabs(s-ss)>ETA) /* s != ss */
		{
		b  = s - 1.0;
		h  = sqrt(3.0 * s - 0.75);
		ss = s;
		}
	for (;;) 
		{
		r = RandomNumber (seed);
		g = r - r * r;
		f = (r - 0.5) * h / sqrt(g);
		x = b + f;
		if (x <= 0.0) 
			continue;
		r = RandomNumber (seed);
		d = 64 * r * r * g * g * g;
		if (d * x < x - 2.0 * f * f || log(d) < 2.0 * (b * log(x / b) - f))  
			break;
		}
		
	return (x);
   
}





/*---------------------------------------------------------------------------------
|
|   SetQvalue
|
|   The Pade method for calculating the matrix exponential, tMat = e^{qMat * v}, 
|   has an error, e(p,q), that can be controlled by setting p and q to appropriate
|   values. The error is:
|
|      e(p,q) = 2^(3-(p+q)) * ((p!*q!) / (p+q)! * (p+q+1)!)
|
|   Setting p = q will minimize the error for a given amount of work. This function 
|   assumes that p = q. The function takes in as a parameter the desired tolerance
|   for the accuracy of the matrix exponentiation, and returns qV = p = q, that
|   will achieve the tolerance. The Pade approximation method is described in:
|  
|   Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
|      The Johns Hopkins University Press, Baltimore, Maryland.
|
|   The function is called from TiProbsUsingPadeApprox.
|
---------------------------------------------------------------------------------*/
int SetQvalue (MrBFlt tol)

{

	int			qV;
	MrBFlt 		x;
	
	x = pow(2.0, 3.0 - (0 + 0)) * Factorial(0) * Factorial (0) / (Factorial(0+0) * Factorial (0+0+1));
	qV = 0;
	while (x > tol)
		{
		qV++;
		x = pow(2.0, 3.0 - (qV + qV)) * Factorial(qV) * Factorial (qV) / (Factorial(qV+qV) * Factorial (qV+qV+1));
		}
		
	return (qV);

}





/*---------------------------------------------------------------------------------
|
|   SetToIdentity
|
|   Make a dim X dim identity matrix.
|
---------------------------------------------------------------------------------*/
void SetToIdentity (int dim, MrBFlt **matrix)

{

	int			row, col;

	for (row=0; row<dim; row++)
		for (col=0; col<dim; col++)
      		matrix[row][col] = (row == col ? 1.0 : 0.0);
	
}





/*---------------------------------------------------------------------------------
|
|   Tha
|
|   Calculate Owen's (1956) T(h,a) function, -inf <= h, a <= inf, 
|   where h = h1/h2, a = a1/a2, from the program of: 
|
|   Young, J. C. and C. E. Minder.  1974.  Algorithm AS 76.  An integral  
|      useful in calculating non-central t and bivariate normal  
|      probabilities.  Appl. Statist., 23:455-457.  [Correction: Appl.  
|      Statist., 28:113 (1979).  Remarks: Appl. Statist. 27:379 (1978),
|      28: 113 (1979), 34:100-101 (1985), 38:580-582 (1988)]  
|
|   See also: 
|
|   Johnson, N. L.  and S. Kotz.  1972.  Distributions in statistics:
|      multivariate distributions.  Wiley and Sons.  New York.  pp. 93-100.
|
---------------------------------------------------------------------------------*/
MrBFlt Tha (MrBFlt h1, MrBFlt h2, MrBFlt a1, MrBFlt a2)

{

	int 			ng = 5, i;
	MrBFlt	 		U[] = {0.0744372, 0.2166977, 0.3397048, 0.4325317, 0.4869533},
					R[] = {0.1477621, 0.1346334, 0.1095432, 0.0747257, 0.0333357},
					pai2 = 6.283185307, tv1 = 1e-35, tv2 = 15.0, tv3 = 15.0, tv4 = 1e-5,
					a, h, rt, t, x1, x2, r1, r2, s, k, sign = 1.0;

	if (fabs(h2) < tv1) 
		return (0.0);
	h = h1 / h2;
	if (fabs(a2) < tv1) 
		{
		t = CdfNormal(h);
		if (h >= 0.0) 
			t = (1.0 - t) / 2.0;
		else      
			t /= 2.0;
		return (t*(a1 >= 0.0 ? 1.0 : -1.0));
		}
	a = a1 / a2;
	if (a < 0.0) 
		sign = -1.0;  
	a = fabs(a);  
	h = fabs(h);   
	k = h*a;
	if (h > tv2 || a < tv1) 
		return (0.0);
	if (h < tv1) 
		return (atan(a)/pai2*sign);
	if (h < 0.3 && a > 7.0) /* (Boys RJ, 1989) */
		{             
		x1 = exp(-k*k/2.0)/k;
		x2 = (CdfNormal(k)-0.5)*sqrt(pai2);
		t = 0.25 - (x1+x2)/pai2*h + ((1.0+2.0/(k*k))*x1+x2)/(6.0*pai2)*h*h*h;
		return (MAX(t,0)*sign);
		}
	t = -h*h / 2.0;  
	x2 = a;  
	s = a*a;
	if (log(1.0+s)-t*s >= tv3) 
		{
		x1 = a/2;  
		s /= 4.0;
	for (;;) /* truncation point by Newton iteration */
		{        
		x2 = x1 + (t*s+tv3-log(s+1.0)) / (2.0*x1*(1.0/(s+1.0)-t));
		s = x2*x2;
		if (fabs(x2-x1) < tv4) 
			break;
		x1 = x2;
		}
	}
	for (i=0,rt=0; i<ng; i++) /* Gauss quadrature */
		{          
		r1 = 1.0+s*SQUARE(0.5+U[i]);
		r2 = 1.0+s*SQUARE(0.5-U[i]);
		rt+= R[i]*(exp(t*r1)/r1 + exp(t*r2)/r2);
		}
		
	return (MAX(rt*x2/pai2,0)*sign);

}





/*---------------------------------------------------------------------------------
|
|   TiProbsUsingEigens
|
---------------------------------------------------------------------------------*/
void TiProbsUsingEigens (int dim, MrBFlt *cijk, MrBFlt *eigenVals, MrBFlt v, MrBFlt r, MrBFlt **tMat, MrBFlt **fMat, MrBFlt **sMat)

{

	int				i, j, s;
	MrBFlt			sum, sumF, sumS, *ptr, EigValexp[192];

	for (s=0; s<dim; s++)
		EigValexp[s] = exp(eigenVals[s] * v * r);

	ptr = cijk;
	for (i=0; i<dim; i++)
		{
		for (j=0; j<dim; j++)
			{
			sum = 0.0;
			for(s=0; s<dim; s++)
				sum += (*ptr++) * EigValexp[s];
			tMat[i][j] = (sum < 0.0) ? 0.0 : sum;
			}
		}
		
#	if 0
	for (i=0; i<dim; i++)
		{
		sum = 0.0;
		for (j=0; j<dim; j++)
			{
			sum += tMat[i][j];
			}
		if (sum > 1.0001 || sum < 0.9999)
			{
			MrBayesPrint ("%s   Warning: Transition probabilities do not sum to 1.0 (%lf)\n", spacer, sum);
			}
		}
#	endif
	
	if (fMat != NULL && sMat != NULL)
		{
		ptr = cijk;
		for (i=0; i<dim; i++)
			{
			for (j=0; j<dim; j++)
				{
				sumF = sumS = 0.0;
				for(s=0; s<dim; s++)
					{
					sumF += (*ptr  ) * eigenVals[s] *                r *     EigValexp[s];
					sumS += (*ptr++) * eigenVals[s] * eigenVals[s] * r * r * EigValexp[s];
					}
				fMat[i][j] = sumF;
				sMat[i][j] = sumS;
				}
			}
		}
		
}





/*---------------------------------------------------------------------------------
|
|   TiProbsUsingPadeApprox
|
|   The method approximates the matrix exponential, tMat = e^{qMat * v}, using
|   the Pade approximation method, described in:
|  
|   Golub, G. H., and C. F. Van Loan. 1996. Matrix Computations, Third Edition.
|      The Johns Hopkins University Press, Baltimore, Maryland.
|
|   The method approximates the matrix exponential with accuracy tol.
|
---------------------------------------------------------------------------------*/
void TiProbsUsingPadeApprox (int dim, MrBFlt **qMat, MrBFlt v, MrBFlt r, MrBFlt **tMat, MrBFlt **fMat, MrBFlt **sMat)

{

	int			qValue;
	MrBFlt		**a, tol;
	
	tol = 0.0000001;
	
	a = AllocateSquareDoubleMatrix (dim);
	
	MultiplyMatrixByScalar (dim, qMat, v * r, a);

	qValue = SetQvalue (tol);

	ComputeMatrixExponential (dim, a, qValue, tMat);
	
	FreeSquareDoubleMatrix (a);

	if (fMat != NULL && sMat != NULL)
		{
		MultiplyMatrices (dim, qMat, tMat, fMat);
		MultiplyMatrices (dim, qMat, fMat, sMat);
		}
				
}

