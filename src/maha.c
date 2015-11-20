/* This file is part of teigen.
   Copyright (C) 2012-2015  Jeffrey L. Andrews and Paul D. McNicholas

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA  02110-1301, USA. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "teigen.h"

/* cholesky computes the Cholesky decomposition of a symmetric
   positive definite (square) matrix 'A' with dimensions 'n' by 'n'
   and stores the result in 'L'.

   Source: http://rosettacode.org/wiki/Cholesky_decomposition#C */
static void cholesky(double *A, int n, double* L) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < (i+1); j++) {
            double s = 0;
	    int k;

            for (k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
		sqrt(A[i * n + i] - s) :
		(1.0 / L[j * n + j] * (A[i * n + j] - s));
        }
    }
}

/* sum_and_sq adds the square of each element in the vector and return
   the result */
static double sum_and_sq(double *vec, int size) {
    double sum = 0;
    int i;

    for (i = 0; i < size; i++) {
	sum += (vec[i]*vec[i]);
    }

    return sum;
}

/* mymaha calculates the Mahalanobis distance of all rows in x and vector mu
   with respect to cov. Col and row are the number of columns and rows in x respectively
   and nsq is the number of rows and columns in the square matrix cov. The out vector
   stores the resulting calculations.
   Note: matrices x and cov have been converted into vectors by taking elements
   column-wise.
   For example, the C array {1, 2, 3, 4, 5, 6, 7, 8, 9} represents the matrix
   {1 4 7}
   {2 5 8}
   {3 6 9}

   Source of algorithm and various comments: the mvnfast package available on CRAN,
   translated from C++ in mahaCpp.cpp file.

   #Equivalent R function:

   .fastMahalanobis <- function(X, mean, mcov)
   {
   dec <- chol(mcov)
   tmp <- forwardsolve(t(dec), t(X) - mean )
   colSums( tmp ^ 2 )
   }
*/
void mymaha(double *x, int *col, int *row, double *mu, double *cov, int *nsq, double *out) {
    int xcol = *col; int xrow = *row; int n = *nsq;

    double *buff = (double *)calloc(n*n, sizeof(double));
    double *diag = (double *)calloc(n, sizeof(double));
    double *tmp = (double *)calloc(xcol, sizeof(double));

    int i, icol, irow;

    // calculate transposed cholesky decomposition
    // buffer array used to hold data between calculations
    cholesky(cov, n, buff);
    transpose(buff, &n, &n, cov);

    // get diagonal entries of transposed chol. decomp.
    for (i = 0; i < n; i++)
	diag[i] = cov[i*n + i];

    // for each of the "n" random vectors, forward solve the
    // corresponding linear system since we are using the lower
    // triangle Cholesky
    for (icol = 0; icol < xrow; icol++) {
        for (irow = 0; irow < xcol; irow++) {
	    double acc = 0.0;
	    int ii;

	    for (ii = 0; ii < irow; ii++)
		    
		acc += tmp[ii]*cov[RC2IDX(irow, ii, n)];

	    tmp[irow] = (x[RC2IDX(icol, irow, xrow)] - mu[irow] - acc) / diag[irow];

	}

	out[icol] = sum_and_sq(tmp, xcol);
    }

    free(tmp);
    free(diag);
    free(buff);
}

#ifdef DO_MAIN
//print an nxn matrix to the console nicely
void show_matrix(double *A, int n) {
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
	    printf("%2.5f ", A[j * n + i]);
	}
        printf("\n");
    }
}

int main(void) {
    double x[] = {-0.58,-1.06,0.08,0.28,-1.58,-0.83,-1.83,-0.45,-0.15,
		  -0.28,0.29,-1.15,-0.21,-1.37,-0.85,-1.68,-0.17,1.16,-0.12,1.01,0.16,
		  -1.05,-1.56,-0.96,0.88,1.74,-0.95,1.09,-0.66,-0.52,1.72,0.29,0.65,-0.18,
		  -2.64,-1.04,-0.25,0.1,-1.24,-0.97,-2.7,-1.36,-0.35,1.12,-0.32,1.76,0.26,
		  0.69,0.57,1.4,0.05,0.09,0.25,0.66,-1.1,1.54,0.52,-0.56,1.57,-0.34};

    int xcol = 3; int xrow = 20;
    double mu[] = {-0.47, -0.27, 0.19};

    double cov[] = {0.69,0.02,-0.06,0.02,1.25,0.17,-0.06,0.17,1.22};
    int n = 3; // size of square matrix cov (row, column)

    show_matrix(cov, 3);

    double *out = (double *)calloc(xrow, sizeof(double));

    int i;

    mymaha(x, &xcol, &xrow, mu, cov, &n, out);

    for (i = 0; i<xrow; i++)
	printf("%f\n", out[i]);

    free(out);
    return 0;
}
#endif /* DO_MAIN */
