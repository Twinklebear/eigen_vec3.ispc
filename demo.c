#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include "eigen_vec3.h"

// Generate a random 3x3 matrix with values between [-100, 100]
void random_matrix(float *m){
	float vals[6] = { 0 };
	for (int i = 0; i < 6; ++i){
		vals[i] = rand() / (float)RAND_MAX;
		vals[i] *= 200;
		vals[i] -= 100;
	}
#define AT(m, r, c) ((m)[(r) * 3 + (c)])
	AT(m, 0, 0) = vals[0];
	AT(m, 0, 1) = vals[1];
	AT(m, 1, 0) = vals[1];
	AT(m, 0, 2) = vals[2];
	AT(m, 2, 0) = vals[2];

	AT(m, 1, 1) = vals[3];
	AT(m, 1, 2) = vals[4];
	AT(m, 2, 1) = vals[4];

	AT(m, 2, 2) = vals[5];
#undef AT
}

int main(int argc, char **argv){
	if (argc < 2){
		fprintf(stderr, "Usage: ./<exe> num_trials\nProgram will run num_trials of solving 100000 random matrices\n");
		return 1;
	}
	srand(time(0));
	const int N = 100000;
	const int TRIALS = atoi(argv[1]);

	// Allocate enough room for all our matrices and results and fill them
	float *mats_ispc = malloc(sizeof(float) * 9 * N);
	float *vecs_ispc = malloc(sizeof(float) * 9 * N);
	float *vals_ispc = malloc(sizeof(float) * 3 * N);

	unsigned long long ispc_elapsed = 0;
	for (int trial = 0; trial < TRIALS; ++trial){
		for (int i = 0; i < N; ++i){
			random_matrix(&mats_ispc[i * 9]);
		}

		// Time the ispc solving of all matrices in parallel
		struct timeval start, end;
		gettimeofday(&start, NULL);
		eigen_vec3(mats_ispc, vecs_ispc, vals_ispc, N);
		gettimeofday(&end, NULL);
		ispc_elapsed += 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000;
	}
	printf("ISPC Kopp impl computing %d matrices took %llums\n", N * TRIALS, ispc_elapsed);

	free(mats_ispc);
	free(vecs_ispc);
	free(vals_ispc);
	return 0;
}

