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

