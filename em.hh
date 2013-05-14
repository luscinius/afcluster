#pragma once

/**
 * Hard EM clustering
 * K = number of clusters
 * N = number of samples
 * row_length = row length of the data matrix
 * data = data matrix
 * assignment = vector of cluster assignment; must be pre-allocated
 * Z = flat array for storing assignment confidence; must be pre-allocated
 * num_trials = number of times to repeat clustering before chosing the result
 * 	with the minimal distortion
 * dist = distance: e -- Euclidean, k -- KL, s -- symmetrized KL,
 * 		r --- KL using raw counts
 * verbose = verbosity level
 * Returns the number of clusters formed
 */
int hard_em(int K, int N, int row_length, double** data, int num_nt,
double** data_1, int* assignment, double* Z, int num_trials=1,
char dist_type='e', int verbose=0);


/**
 * Consensus EM clustering
 * Arguments shared with the regular kmeans utility have the same meanings.
 * Regular k-means clustering is repeated num_trials times, each time choosing
 * boot_frac*N samples for the bootstrap dataset. Distance matrix is formed
 * as D[i,j] = 1 - #(i and j were clustered together) / #(i and j were in the 
 * same bootstrap dataset). Then hierarchical clustering is performed using 
 * this distance matrix. This procedure mitigates the randomness of the 
 * initialization in k-means algotithm and reduces the effect of the outliers.
 */
int cons_em(int K, int N, int row_length, double** data, int num_nt,
double** data_1, int* assignment, int num_trials, char dist_type='e', 
int verbose=0, double boot_frac=0.80, double cons_threshold=10.);


/**
 * Soft EM clustering using KL distance
 * Args shared with hard_em routine have the same meaning;
 * Z must be pre-allocated as well.
 * threshold = threshold for centroid shift which causes to stop iterations
 */
void soft_em(int K, int N, int row_length, double** data, int num_nt,
double** data_1, double* Z, double threshold, int verbose);
