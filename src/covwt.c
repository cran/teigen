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
#include <math.h>
#include <stdio.h>

#include "teigen.h"

/* 	mycovwt calculates the estimated weighted covariance matrix based on the
	implementation in R. The algorithm roughly follows these steps:
		1) wt <- wt/sum(wt)
		2) x <- sqrt(wt) * (x - t(array(center, dim(x)[2:1])))
		3) t(x) %*% x
	'x' is a matrix with 'm' rows and 'n' columns in vector form, while 'wt'
	is a vector if weights for each obervation with length 'm'. 'center' is
	a vector of length 'size' specifying the centers to be used when
	computing covariances. Resulting covariance matrix stored in 'rv'. */
void mycovwt(double *x, int *m, int *n, double *wt, double *center, int *size, double *rv){
	double sum = 0;
	int row = *m; 
	int col = *n;
	int i;
	
	for (i = 0; i < row; i++){
		sum += wt[i];
	}
	
	for (i = 0; i < row; i++){
		wt[i] /= sum;
		wt[i] = sqrt(wt[i]);
	}

	// same operation in R as rep_len(center, row*col)
	double arr[row*col];
	for (i = 0; i < row*col; i++){
		arr[i] = center[i % *size];
	}
	double arr_t[row*col];
	
	// tranpose with 'col' as the row length and 'row' as the col length
	// because the sweep function in R (for teigen situations) uses
	// 'array(center, dim(x)[2:1])' which gives a resulting matrix with 
	// the row and column dimensions of x reversed.
	transpose(arr, &col, &row, arr_t);
	for (i = 0; i < row*col; i++){
		x[i] = (x[i] - arr_t[i]) * wt[i%row];
	}
	
	//use arr_t to hold transpose of x
	transpose(x, &row, &col, arr_t);
	
	//multiply the transpose of x by x
	matrix_mult(arr_t, x, col, row, col, rv);
	
}

#ifdef DO_MAIN
/* show_matrix is a handy function for displaying a matrix to the console
   in proper visual form when given a matrix and its dimensions */
void show_matrix(double *A, int m, int n) {
	int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
           printf("%2.5f ", A[j * m + i]);
        printf("\n");
    }
}

int main(void) {
	// example numbers for calculating mycovwt
	double x[] = {1,2,3,4,5,6,7,8};
	//dimensions of x
	int r = 4; 
	int c = 2; 
	double wt[] = {56,57,58,59};
	double center[] = {1.23, 4.56};
	
	double* rv = (double*)calloc(c*c, sizeof(double));
	
	mycovwt(x, &r, &c, wt, center, &c, rv);
	show_matrix(rv, c, c);
	
	// should be (calculated in R):
	double correct[] = {2.91812 , 3.78358 , 3.78358 , 5.09795};
	printf("\n\nThe correct matrix should be:\n");
	show_matrix(correct, c, c);
	
	free(rv);
	return 0;
}
#endif /* DO_MAIN */
