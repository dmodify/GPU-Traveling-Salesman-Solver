#include "util.h"
#include "assert.h"
#include <cuda.h>
#include <stdlib.h>

// 0..N = tour
// N..2N = coords
// 2N..3N = distances
#define TOUR_IDX 0
#define COORDS_IDX 1024
#define DISTS_IDX 3*1024
#define OLD_IDX 4*1024 
#define BEST_IDX 5*1024
#define DELTA_COST 6*1024
#define FINAL_IDX 7*1024

#define GPU_CHECKERROR(err)(gpuCheckError(err, __FILE__, __LINE__))
static void gpuCheckError(cudaError_t err, const char* file, int line){
	if(err != cudaSuccess){
		printf("%s in %s at line %d\n", cudaGetErrorString(err), file,
			 line);
		exit(EXIT_FAILURE);
	}
}
__device__ __forceinline__ int euclid_edgelen(int i, int j, int *coords) {
	int dist;
	float dx, dy;
	dx = coords[2 * i] - coords[2 * j];
	dy = coords[2 * i + 1] - coords[2 * j + 1];
	dist = (int)(sqrtf(dx*dx + dy*dy) + 0.5);
	return dist;
}
__global__ void wipe_shared_mem()
{
	int s = 1024;
	int tidx = threadIdx.x;
	extern __shared__ int tour_coords_dists[];
	tour_coords_dists[tidx] = 0;
	tour_coords_dists[s*1+tidx] = 0;
	tour_coords_dists[s*2+tidx] = 0;
	tour_coords_dists[s*3+tidx] = 0;
	tour_coords_dists[s*4+tidx] = 0;
	tour_coords_dists[s*5+tidx] = 0;
	tour_coords_dists[s*6+tidx] = 0;
	tour_coords_dists[s*7+tidx] = 0;
}

__global__ void two_opt(int *coords, int *record_paths, int N, int seed)
{
	int tidx = threadIdx.x;
	int bidx = blockIdx.x;
	//if (tidx == 18)
	//	printf("hi I'm bidx %d\n", bidx);
	extern __shared__ int tour_coords_dists[];
	
	// erase the memory
	tour_coords_dists[DISTS_IDX + tidx] = 0;
	tour_coords_dists[OLD_IDX + tidx] = 0;
	tour_coords_dists[BEST_IDX + tidx] = 0;
	tour_coords_dists[DELTA_COST + tidx] = 0;

	// copy from global
	tour_coords_dists[COORDS_IDX + 2*tidx] = coords[2*tidx];
	tour_coords_dists[COORDS_IDX + 2*tidx+1] = coords[2*tidx+1];

	tour_coords_dists[TOUR_IDX + tidx] = tidx;
	//tour_coords_dists[TOUR_IDX + tidx] = tour[tidx];

//	int p = 2;
//	while (p < N) {
//		p = p*2;
//	}
//	int pow_N = p;
	int can_improve = 1;
//	int iterations = INT_MAX;

	__syncthreads();
	if (tidx == 0) {
		unsigned long a = 16807;
		unsigned long m = 2147483647;
		unsigned long x = 13 + bidx + seed;

		// "Any one who considers arithmetical methods of producing
		// random digits is, of course, in a state of sin."
		//					-- John von Neumann

		// generate initial tour
		//for (int i = 0; i < N; i++)
		//	tour_coords_dists[i] = i;

		// `randomly' shuffle the tour
		for (int i = 0; i < N/2; i++) {
			x = (a * x) % m;
			unsigned int idx_1 = x % N;
			x = (a * x) % m;
			unsigned int idx_2 = x % N;
	
			int temp = tour_coords_dists[idx_1];
			tour_coords_dists[idx_1] = tour_coords_dists[idx_2];
			tour_coords_dists[idx_2] = temp;
		}
	}

	__syncthreads();

	do {
		for (int edge = 0; edge < N; edge++) {
	//		__syncthreads();
			// everyone gets a copy of the edge we are thinking about swapping
			int delta_cost;
			if (tidx > edge) {
				int edge_len_a = euclid_edgelen(tour_coords_dists[TOUR_IDX + edge], 
							tour_coords_dists[TOUR_IDX + (edge + 1)%N], 
							tour_coords_dists+COORDS_IDX);
		//		__syncthreads();
				// everyone looks at all the other edges that could be swapped
				int edge_len_b = euclid_edgelen(tour_coords_dists[TOUR_IDX + tidx],
							tour_coords_dists[TOUR_IDX + (tidx + 1)%N],
							tour_coords_dists+COORDS_IDX);

		//		__syncthreads();
				// consider two replacements for edge_a and edge_b
				int poss_edge_1 = euclid_edgelen(tour_coords_dists[TOUR_IDX+edge],
								tour_coords_dists[TOUR_IDX+tidx],
								tour_coords_dists+COORDS_IDX);

		//		__syncthreads();
				int poss_edge_2 = euclid_edgelen(tour_coords_dists[TOUR_IDX+(edge+1)%N]
								,tour_coords_dists[TOUR_IDX+(tidx+1)%N]
								,tour_coords_dists+COORDS_IDX);
				
		//		__syncthreads();
				delta_cost = poss_edge_1 + poss_edge_2 - edge_len_a - edge_len_b;

		//		__syncthreads();
				// cannot swap neigbor edges
				if (tidx == edge|| tidx == (edge-1)%N 
						|| tidx == (edge+1)%N || (tidx+1)%N == edge)
					delta_cost = INT_MAX;
			}
				
			__syncthreads();

			// write deltas to shared memory
			//printf("%d tidx %d delta cost \n", tidx, delta_cost);
			tour_coords_dists[DISTS_IDX + tidx] = delta_cost;
//			__syncthreads();
			tour_coords_dists[OLD_IDX + tidx] = tidx; 
			__syncthreads();

/*
			// reduce and get the index at which the optimum swap exists
			for (int stride = pow_N/2; stride > 32; stride /= 2) {
				if (tidx < stride) {
					if (tour_coords_dists[DISTS_IDX+tidx] >
						tour_coords_dists[DISTS_IDX+tidx+stride]) {

						tour_coords_dists[OLD_IDX+tidx] =
						tour_coords_dists[OLD_IDX+tidx+stride];

						tour_coords_dists[DISTS_IDX+tidx] = 
						tour_coords_dists[DISTS_IDX+tidx+stride];
					}

				}
				__syncthreads();
			}
*/
			if (tidx < 512) {
				if (tour_coords_dists[DISTS_IDX+tidx] >
						tour_coords_dists[DISTS_IDX+tidx+512]) {

						tour_coords_dists[OLD_IDX+tidx] =
						tour_coords_dists[OLD_IDX+tidx+512];

						tour_coords_dists[DISTS_IDX+tidx] = 
						tour_coords_dists[DISTS_IDX+tidx+512];


				}
			}
			__syncthreads();

			if (tidx < 256) {
				if (tour_coords_dists[DISTS_IDX+tidx] >
						tour_coords_dists[DISTS_IDX+tidx+256]) {

						tour_coords_dists[OLD_IDX+tidx] =
						tour_coords_dists[OLD_IDX+tidx+256];

						tour_coords_dists[DISTS_IDX+tidx] = 
						tour_coords_dists[DISTS_IDX+tidx+256];
				}
			}
			__syncthreads();

			if (tidx < 128) {
				if (tour_coords_dists[DISTS_IDX+tidx] >
						tour_coords_dists[DISTS_IDX+tidx+128]) {

						tour_coords_dists[OLD_IDX+tidx] =
						tour_coords_dists[OLD_IDX+tidx+128];

						tour_coords_dists[DISTS_IDX+tidx] = 
						tour_coords_dists[DISTS_IDX+tidx+128];
				}
			}
			__syncthreads();

			if (tidx < 64) {
				if (tour_coords_dists[DISTS_IDX+tidx] >
						tour_coords_dists[DISTS_IDX+tidx+64]) {

						tour_coords_dists[OLD_IDX+tidx] =
						tour_coords_dists[OLD_IDX+tidx+64];

						tour_coords_dists[DISTS_IDX+tidx] = 
						tour_coords_dists[DISTS_IDX+tidx+64];
				}
			}
			__syncthreads();

			// unroll the last warp
			if (tidx < 32) {
				if (tour_coords_dists[DISTS_IDX+tidx] > 
					tour_coords_dists[DISTS_IDX+tidx+32]) {

					tour_coords_dists[OLD_IDX+tidx] =
					tour_coords_dists[OLD_IDX+tidx+32];

					tour_coords_dists[DISTS_IDX+tidx] = 
					tour_coords_dists[DISTS_IDX+tidx+32];
				}
				if (tour_coords_dists[DISTS_IDX+tidx] > 
					tour_coords_dists[DISTS_IDX+tidx+16]) {

					tour_coords_dists[OLD_IDX+tidx] =
					tour_coords_dists[OLD_IDX+tidx+16];

					tour_coords_dists[DISTS_IDX+tidx] = 
					tour_coords_dists[DISTS_IDX+tidx+16];
				}
				if (tour_coords_dists[DISTS_IDX+tidx] > 
					tour_coords_dists[DISTS_IDX+tidx+8]) {

					tour_coords_dists[OLD_IDX+tidx] =
					tour_coords_dists[OLD_IDX+tidx+8];

					tour_coords_dists[DISTS_IDX+tidx] = 
					tour_coords_dists[DISTS_IDX+tidx+8];
				}
				if (tour_coords_dists[DISTS_IDX+tidx] > 
					tour_coords_dists[DISTS_IDX+tidx+4]) {

					tour_coords_dists[OLD_IDX+tidx] =
					tour_coords_dists[OLD_IDX+tidx+4];

					tour_coords_dists[DISTS_IDX+tidx] = 
					tour_coords_dists[DISTS_IDX+tidx+4];
				}
				if (tour_coords_dists[DISTS_IDX+tidx] > 
					tour_coords_dists[DISTS_IDX+tidx+2]) {

					tour_coords_dists[OLD_IDX+tidx] =
					tour_coords_dists[OLD_IDX+tidx+2];

					tour_coords_dists[DISTS_IDX+tidx] = 
					tour_coords_dists[DISTS_IDX+tidx+2];
				}
				if (tour_coords_dists[DISTS_IDX+tidx] > 
					tour_coords_dists[DISTS_IDX+tidx+1]) {

					tour_coords_dists[OLD_IDX+tidx] =
					tour_coords_dists[OLD_IDX+tidx+1];

					tour_coords_dists[DISTS_IDX+tidx] = 
					tour_coords_dists[DISTS_IDX+tidx+1];
				}
			}


			__syncthreads();
			// store the optimum swap in the index
			if (tidx == edge) {
				tour_coords_dists[BEST_IDX + edge] =
					tour_coords_dists[OLD_IDX];
				tour_coords_dists[DELTA_COST + edge] = 
					tour_coords_dists[DISTS_IDX];

			}
			__syncthreads();
			
			assert(tour_coords_dists[tidx] >= 0);	
			assert(tour_coords_dists[OLD_IDX+tidx] >= 0);
		
		}
		// reduce. now that we have all the candidates from examining each edge,
		// we want to find the best. Wipe dists and reuse them to get indexes
		tour_coords_dists[FINAL_IDX + tidx] = tidx;
		//__syncthreads();
/*
		for (int stride = pow_N/2; stride > 32; stride /= 2) {
			if (tidx < stride) {
				if (tour_coords_dists[DELTA_COST+tidx] >
					tour_coords_dists[DELTA_COST+tidx+stride]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+stride];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+stride];
				}

			}
			__syncthreads();

		}
*/
		__syncthreads();
		if (tidx < 512) {
			if (tour_coords_dists[DELTA_COST+tidx] >
					tour_coords_dists[DELTA_COST+tidx+512]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+512];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+512];

			}
		}
		__syncthreads();

		if (tidx < 256) {
			if (tour_coords_dists[DELTA_COST+tidx] >
					tour_coords_dists[DELTA_COST+tidx+256]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+256];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+256];
			}
		}
		__syncthreads();

		if (tidx < 128) {
			if (tour_coords_dists[DELTA_COST+tidx] >
					tour_coords_dists[DELTA_COST+tidx+128]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+128];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+128];
			}
		}
		__syncthreads();

		if (tidx < 64) {
			if (tour_coords_dists[DELTA_COST+tidx] >
					tour_coords_dists[DELTA_COST+tidx+64]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+64];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+64];
			}
		}
		__syncthreads();

		if (tidx < 32) {
				if (tour_coords_dists[DELTA_COST+tidx] > 
					tour_coords_dists[DELTA_COST+tidx+32]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+32];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+32];
				}
				if (tour_coords_dists[DELTA_COST+tidx] > 
					tour_coords_dists[DELTA_COST+tidx+16]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+16];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+16];
				}
				if (tour_coords_dists[DELTA_COST+tidx] > 
					tour_coords_dists[DELTA_COST+tidx+8]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+8];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+8];
				}
				if (tour_coords_dists[DELTA_COST+tidx] > 
					tour_coords_dists[DELTA_COST+tidx+4]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+4];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+4];
				}
				if (tour_coords_dists[DELTA_COST+tidx] > 
					tour_coords_dists[DELTA_COST+tidx+2]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+2];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+2];
				}
				if (tour_coords_dists[DELTA_COST+tidx] > 
					tour_coords_dists[DELTA_COST+tidx+1]) {

					tour_coords_dists[FINAL_IDX+tidx] =
					tour_coords_dists[FINAL_IDX+tidx+1];

					tour_coords_dists[DELTA_COST+tidx] = 
					tour_coords_dists[DELTA_COST+tidx+1];
				}
			}
		__syncthreads();
		//assert(tour_coords_dists[tidx] >= 0);	
		//assert(tour_coords_dists[OLD_IDX+tidx] >= 0);
		//assert(tour_coords_dists[OLD_IDX+tidx] >= 0);
		//assert(tour_coords_dists[FINAL_IDX+tidx] >= 0);

		//__syncthreads();
		// reverse the path
		int best_edge = tour_coords_dists[FINAL_IDX];
		int swap_edge = tour_coords_dists[BEST_IDX + best_edge];
		int start, end;
		if (best_edge > swap_edge) {
			start = swap_edge;
			end = best_edge;
		}
		else {
			start = best_edge;
			end = swap_edge;
		}
		
		__syncthreads();

		int n_swap = end - start;
		if (tidx < n_swap/2) {
			
			int temp1 = tour_coords_dists[TOUR_IDX+start+n_swap-tidx];
			int temp2 = tour_coords_dists[TOUR_IDX+start+1+tidx];

			tour_coords_dists[TOUR_IDX+start+n_swap-tidx] = temp2;
			tour_coords_dists[TOUR_IDX+start+1+tidx] = temp1;
		}
		//__syncthreads();

		if (tour_coords_dists[DELTA_COST] >= 0) {
			can_improve = 0;
			//if (tidx == 0)
			//	printf("can't improve\n");
			
		}
		
		__syncthreads();
		
	} while(can_improve);

	// copy back
	record_paths[bidx*N + tidx] = tour_coords_dists[tidx];
}

int main(int argc, char **argv) 
{

	if (argc < 2) {
		fprintf(stderr, "Usage: ./%s tsp_file\n", argv[0]);
		exit(-1);
	}

	int *tour;
	int *coords; // this gets malloc'ed, must free later
	int N;
//	int record = 1000000000;
	int B = 12;
	read_euc2d_file(argv[1], &coords, &N);
	printf("N : %d\n", N);

	int *record_paths = (int*)malloc(B*N*sizeof(int));// store the acutal path 
        if (!record_paths) {
                fprintf(stderr, "malloc failure\n");
                exit(-1);
        }

	tour = (int*)malloc(N*sizeof(int));
	if (!tour) {
		fprintf(stderr, "heap error");
		exit(-1);
	}

	for (int i = 0; i < N; i++)
		tour[i] = i;

	int *d_coords, *d_record, *d_record_paths;
	GPU_CHECKERROR(cudaMalloc((void**)&d_coords, 2*N*sizeof(int)));
	GPU_CHECKERROR(cudaMalloc((void**)&d_record, sizeof(int)));
	GPU_CHECKERROR(cudaMalloc((void**)&d_record_paths, B*N*sizeof(int)));

	GPU_CHECKERROR(cudaMemcpy((void*)d_coords, (void*)coords,
				2*N*sizeof(int), cudaMemcpyHostToDevice));

	
	srand(time(NULL));
	int seed = rand()%10000000;	

	cudaEvent_t start, stop;
	float time;
	GPU_CHECKERROR(cudaEventCreate(&start));
	GPU_CHECKERROR(cudaEventCreate(&stop));

	GPU_CHECKERROR(cudaEventRecord(start, 0));

	int size = (1024 + 2*1024 + 1024 + 6*1024) * sizeof(int);	
	//clock_t start = clock();

	wipe_shared_mem<<<B, 1024, size>>>();
	
	GPU_CHECKERROR(cudaDeviceSynchronize());
	two_opt<<<B, N, size>>>(d_coords, d_record_paths, N, seed);	

	GPU_CHECKERROR(cudaEventRecord(stop, 0));
	GPU_CHECKERROR(cudaEventSynchronize(stop));
	GPU_CHECKERROR(cudaDeviceSynchronize());       

	GPU_CHECKERROR(cudaMemcpy((void*)record_paths, (void*)d_record_paths, 
				B*N*sizeof(int), cudaMemcpyDeviceToHost));

	cudaEventElapsedTime(&time, start, stop);
	printf("*************************************\n");
	printf("Time for the kernel: %'f ms\n", time);
	printf("*************************************\n");


	// find shortest tour
	int best_len = INT_MAX;
	int block_num = 0;
	int *best_tour = (int*)malloc(N*sizeof(int));
	if (!best_tour) {
		fprintf(stderr, "malloc failure");
		exit(-1);
	}

	// loop around one extra time to compare the last tour
	for (int i = 0; i < B*N+1; i++) {
		if (i != 0 && i%N == 0) {
			int temp_len = tour_len(&tour, coords, N);
			if (temp_len < best_len) {
				best_len = temp_len;
				for (int j = 0; j < N; j++) 
					best_tour[j] = tour[j];
				block_num = i/N;
			}
		}
		if (i < B*N) {
			int t = record_paths[i];
			tour[i%N] = t;
		}
	}
	printf("best tour is %d\n", block_num); 
	printf("best lenght is %d\n", best_len);
	int offset = block_num*N;

	printf("\n");

	FILE *fp;
	fp = fopen("out", "w");
	fprintf(fp, "%d\n", N);
	for (int i = 0; i < N; i++)
		fprintf(fp, "%d ", *(record_paths+i+offset));
	fprintf(fp, "\n");
	fclose(fp);
	printf("tour written to out\n");	
	
	free(tour);
	free(record_paths);
	free(coords);

	return 0;
}
