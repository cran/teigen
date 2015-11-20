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

#ifndef _TEIGEN_H_
#define _TEIGEN_H_

/* convert matrix indexing to vector indexing
   r: row number, c: column number, nr: number of rows in matrix */
#define RC2IDX(r,c,nr) ((r) + ((nr)*(c)))

/* transpose computes the transpose the matrix given a matrix 'mat'
   with 'row' rows and 'col' columns and dimensions removed (i.e., the
   vector form of a matrix) and stores the result in 'rv' */
void transpose(double *mat, int *row, int *col, double *rv);

/* matrix_mult multiplies matrix 'mat1' with 'm' rows by 'n' columns with
   matrix 'mat2' with 'n' rows and 'col' columns. Result is stored in the
   'mult_mat' array. */
void matrix_mult(double *mat1, double *mat2, int m, int n, int col, double *mult_mat);

#endif /* _TEIGEN_H_ */
