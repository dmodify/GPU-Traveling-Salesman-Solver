/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 Jeffrey Scholz
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef UTIL_H
#define UTIL_H
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
int read_euc2d_file( char *filename, int **coords, int *dimension );
//__device__ int euclid_edgelen( int i, int j, int *coords );

int tour_len( int *tour, int *coords, int N, int before, int after ); 
int tour_len(int **tour, int *coords, int N);
int next_tour (int** array, int N, int fixed_idx, int skip);
float get_ticks();
void swap( int **array, int i, int j, int N);
void reverse( int **array, int from, int N );
int next_tour_M (int** array, int* coords, int N, int *record);
int tour_len_M (int *tour, int *coords, int N , int part, int before, int after, int record);
void m_reverse(int **array, int from, int to);
__device__ void swap_S(int array[], int i, int j, int N);
__device__ int tour_len_simple_S(int *tour, int *coords, int N, int record);
#endif

