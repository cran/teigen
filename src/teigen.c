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

#include "teigen.h"

void matrix_mult(double *mat1, double *mat2, int m, int n, int col, double *mult_mat){
	double total = 0;
	int i, j, k;
	
	for (k = 0; k < col; k++){
		for (i = 0; i < m; i++){
			for (j = 0; j < n; j++){
				total += mat1[RC2IDX(i, j, m)] * mat2[RC2IDX(j, k, n)];
			}
			mult_mat[RC2IDX(i, k, m)] = total;
			total = 0;
		}
	}
}

void transpose(double *mat, int *row, int *col, double *rv) {
    int i, j;

    for (i = 0; i < *row; i++) {
		for (j = 0; j < *col; j++) {
			rv[j + i*(*col)] = mat[i +j*(*row)];
		}
		/* to help understand rv[j + i*(*col)] = mat[i +j*(*row)], consider:
	       row = 2, col = 4:
	       j + i*col = {0, 1, 2, 3, 4, 5, 6, 7}
	       i + j*row = {0, 2, 4, 6, 1, 3, 5, 7}

	       row = 4, col = 2:
	       j + i*col = {0, 1, 2, 3, 4, 5, 6, 7}
	       i + j*row = {0, 4, 1, 5, 2, 6, 3, 7} */
    }
}
