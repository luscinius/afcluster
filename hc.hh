# pragma once

/** hierarchical clustering routine using min-heap
 * K = target number of clusters
 * N = number of data points
 * assignment = vector to store cluster assignment in
 * dist_cutoff = cutoff for maximal distance; stop clustering should we reach
 * 	it
 * distmatrix = upper triangular distance matrix
 * heap = upper triangular matrix for the heap of distance indices
 * 
 * It is assumed that assignment, distmatrix and heap are pre-allocated
 * and that those will be de-allocated elsewhere
 */
int hcluster(int K, int N, int* assignment, double dist_cutoff, 
	double** distmatrix, int** heap, int verbose);
