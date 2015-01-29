#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <sys/time.h>
#include "eigen_vec3.h"
#include "dsyevh3.h"

/*
 * Compare our performance both speed and error-wise against Kopp's
 * C implementation. To run this demo you'll need to download their
 * source code for the hybrid method here: www.mpi-hd.mpg.de/personalhomes/globes/3x3/
 * and drop it into the project's source code.
 * These benchmarks likely aren't very rigourous.
 */

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
	const float TOL = 0.00005;
	const int N = 100000;
	const int TRIALS = atoi(argv[1]);
	// Count number of times vectors or values differed
	int vec_differing = 0;
	int val_differing = 0;

	// Allocate enough room for all our matrices and results and fill them
	float *mats_ispc = malloc(sizeof(float) * 9 * N);
	float *vecs_ispc = malloc(sizeof(float) * 9 * N);
	float *vals_ispc = malloc(sizeof(float) * 3 * N);

	double *mats = malloc(sizeof(double) * 9 * N);
	double *vecs = malloc(sizeof(double) * 9 * N);
	double *vals = malloc(sizeof(double) * 3 * N);

	unsigned long long ispc_elapsed = 0, baseline_elapsed = 0;
	for (int trial = 0; trial < TRIALS; ++trial){
		for (int i = 0; i < N; ++i){
			random_matrix(&mats_ispc[i * 9]);
		}
		// Copy our random matrices over to double precision
		for (int i = 0; i < N; ++i){
			for (int k = 0; k < 3; ++k){
				for (int j = 0; j < 3; ++j){
					mats[i * 9 + k * 3 + j] = mats_ispc[i * 9 + k * 3 + j];
				}
			}
		}

		// Time the ispc solving of all matrices in parallel
		struct timeval start, end;
		gettimeofday(&start, NULL);
		eig_vec_test(mats_ispc, vecs_ispc, vals_ispc, N);
		gettimeofday(&end, NULL);
		ispc_elapsed += 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000;

		// Compare against serial C implementation, we should definitely be faster since we're solving in parallel
		// but we need to run this anyway for an accuracy benchmark, so why not compare performance as well
		gettimeofday(&start, NULL);
		for (int i = 0; i < N; ++i){
			dsyevh3(&mats[i * 9], &vecs[i * 9], &vals[i * 3]);
		}
		gettimeofday(&end, NULL);
		baseline_elapsed += 1000 * (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000;

		// Compare accuracy of results for this trial using Kopp's implementation as our baseline
		for (int i = 0; i < N; ++i){
			for (int k = 0; k < 3; ++k){
				float diff[3] = { fabs(vecs[i * 9 + k] - vecs_ispc[i * 9 + k]),
					fabs(vecs[i * 9 + 3 + k] - vecs_ispc[i * 9 + 3 + k]),
					fabs(vecs[i * 9 + 6 + k] - vecs_ispc[i * 9 + 6 + k]) };
				float val_diff = fabs(vals[i * 3 + k] - vals_ispc[i * 3 + k]);
				if (diff[0] > TOL || diff[1] > TOL || diff[2] > TOL){
					++vec_differing;
				}
				if (val_diff > TOL){
					++val_differing;
				}
			}
		}
	}
	printf("ISPC Kopp impl computing %d matrices took %llums\n", N * TRIALS, ispc_elapsed);
	printf("C Kopp impl computing %d matrices took %llums\n", N * TRIALS, baseline_elapsed);

	free(mats_ispc);
	free(vecs_ispc);
	free(vals_ispc);
	free(mats);
	free(vecs);
	free(vals);

	printf("%% of eigenvectors within %f of baseline = %.2f%%\n", TOL, 100.f - (vec_differing * 100.f) / (3 * N * TRIALS));
	printf("%% of eigenvalues within %f of baseline = %.2f%%\n", TOL, 100.f - (val_differing * 100.f) / (3 * N * TRIALS));
	return 0;
}

