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

#include "util.h"
#include <stdbool.h>
#include <cuda.h>

int euclid_edgelen_h(int i, int j, int *coords);
int tour_len(int **tour, int *coords, int N)
{
	if(N == 0)
		return 0;
	int len = 0;
	for (int i=0;i<N;i++)
	{
		len += euclid_edgelen_h(*(*tour+i), *(*tour+(i+1)%N), coords);
	}
	return len;
}

__device__ void reverse_S(int array[], int from, int N)
{
	int head = from;
	int tail = N - 1;

	while(head < tail)
	{
		array[head] = array[tail] + array[head];
		array[tail] = array[head] - array[tail];
		array[head] = array[head] - array[tail];
		++head;
		--tail;
	}
}
__device__ void swap_S(int array[], int i, int j, int N)
{
	if (i == j)
		return;
	if(i<N && j<N)
	{
		array[i] = array[i] + array[j];
		array[j] = array[i] - array[j];
		array[i] = array[i] - array[j];
	}
	else
	{
		printf("swap fail: out of bound\n");
	}
}

int euclid_edgelen_h(int i, int j, int *coords) 
{
	int dist;
	float dx, dy;
	dx = coords[2 * i] - coords[2 * j];
	dy = coords[2 * i + 1] - coords[2 * j + 1];
	dist = (int)(sqrtf(dx*dx + dy*dy) + 0.5);
	return dist;
}


int read_euc2d_file(char *filename, int **coords, int *dimension) {
	FILE *fp;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("file read error (1)\n");
		exit(-1);
	}

	while (!feof(fp)) {
		char buffer[256];
		if(!fgets(buffer, 256, fp)) {
			if (!feof(fp)) {
				fprintf(stderr, "file read failure\n");
				exit(-1);
			}
		}
		if (strstr(buffer, "DIMENSION")) {
			sscanf(buffer, "DIMENSION : %d", dimension);
		}
	}
	rewind(fp);
	*coords = (int*)malloc(*dimension * 2 * sizeof(int));

	if (coords == NULL) {
		fprintf(stderr, "out of heap memory");
		exit(-1);
	}

	bool ncs = false;
	while (!feof(fp)) {
		char buffer[256];
		int i, x, y;
		if(!fgets(buffer, 256, fp)) {
			if (!feof(fp)) {
				fprintf(stderr, "file read failure\n");
				exit(-1);
			}
		}
       
		if (!ncs && !strstr(buffer, "NODE_COORD_SECTION"))
			continue;
		else if (!ncs && strstr(buffer, "NODE_COORD_SECTION")) {
			ncs = true;
			continue;
		}
		else {
			sscanf(buffer, "%d %d %d", &i, &x, &y);
			i--;
			(*coords)[i * 2] = x;
			(*coords)[i * 2 + 1] = y;
		}

	}
	return 0;
}
