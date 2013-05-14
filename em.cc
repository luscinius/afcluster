#include "em.hh"
#include "freqs.hh"
#include "dists.hh"
#include "centroid.hh"
#include "hc.hh"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;


/**
 * Evaluate confidence of assignment to each centroid
 */
static void eval_confidence(int K, int N, int row_length, double **data,
double **centroids, double *Z)
{
	static const double KL_DIST_CUTOFF = 10; // set exponential to zero
		// if we exceed this negative exponent
	for (int n=0; n<N; n++){
		double* row = Z + n*K;
		for(int k=0; k<K; k++){
			row[k] = kl_distance(centroids[k], data[n], row_length, NULL);
		}
		double min_dist = *row;
		for(int k=1; k<K; k++){
			if (min_dist > row[k]){
				min_dist = row[k];
			}
		}
		double S = 0;
		for(int k=0; k<K; k++){
			row[k] -= min_dist;
			S += row[k] = row[k] > KL_DIST_CUTOFF ? 0 : exp(-row[k]);
		}
		for(int k=0; k<K; k++){
			row[k] /= S;
		}
	}
	return;
}


// Auxiliary function to select procedures for evaluation of 
// distances and centroids
void select_dists_cent(char dist_type, double (**distf)(double*, double*, int,
double*), void (**eval_centroid)(int, int, double**, int, double**, double*,
double*))
{
	switch (dist_type){
		case 'a':
			*distf = &d2ast_distance;
			*eval_centroid = &d2ast_centroid;
			break;
		case 'c':
			*distf = &chi2_distance;
			*eval_centroid = &mm_centroid;
			break;
		case 'd':
			*distf = &d2_distance;
			*eval_centroid = &d2_centroid;
			break;
		case 'k':
			*distf = &kl_distance;
			*eval_centroid = &kl_centroid;
			break;
		case 's':
			*distf = &symkl_distance;
			*eval_centroid = &kl_centroid;
			break;
		case 'e': // fall through
		default:
			*distf = &euclidean_distance;
			*eval_centroid = &mean_centroid;
	}
	return;
}


/** EM clustering routine
 * This is an auxiliary routine, and it assumes that the memory for the 
 * necessary computations is already pre-allocated. This routine is called
 * from other routines. Calling routines can do a single clustering run
 *  or multiple runs with consensus clustering.
 * K = number of clusters
 * N = number of samples
 * row_length = row length
 * data = data matrix
 * assignment = assignment of samples to clusters
 * numMembers = number of members in each cluster
 * centroids = matrix for centroids
 * tmp_data = array of double* of length N for storing elements in each cluster;
 * 		used in centroid calculation
 * distf = function for evaluating distance
 * eval_centroid = function for evaluating centroids
 * 
 * Return value is the distortion (sum of the distances from each data 
 * point to its corresponding centroid; in the case of the squared 
 * Euclidean distance it gives intra-cluster variance).
 */
static double em_routine(int K, int N, int row_length, double** data,
int num_nt, double** data_1, int* assignment, int* numMembers,
double** centroids, double** centroids_tilde, double** tmp_data,
double**tmp_data_1, char dist_type,
int verbose=0, int allow_empty_clusters = 0)
{
	// choose auxiliary functions
	double (*distf)(double*, double*, int, double*);
	void (*eval_centroid)(int, int, double**, int, double**, double*, double*);
	select_dists_cent(dist_type, &distf, &eval_centroid);

	while (true){ // repeat clustering attempts till we get a result
restart:
		if (verbose > 0) {
			cerr<<"\tEM routine: starting clustering attempt"<<endl;
		}
		memset(assignment, 0, N*sizeof(*assignment));
		memset(numMembers, -1, K*sizeof(*numMembers));
		// initial centroid assignment --- initialize with randomly chosen 
		// data items
		for (int k=0; k<K; ++k) {
			int point_number = rand()% N;
			eval_centroid(1, row_length, data + point_number, num_nt, data_1 +
			point_number, centroids[k], centroids_tilde[k]);
		}

		// start iterations of EM clustering
		while (true){
			// compute the distances to the new centroids 
			// and compute the new assignments
			bool assignmentChanged = false;  // set to true if 
			// at least one element changes its cluster assignment
			for(int n=0; n<N; n++){ // update assignment of each element
				int new_assignment=-1;  // new assignment
				double min_dist=-1;  // distance to the closest centroid
				// initial value set here is irrelevant as it 
				// gets changed later
				for(int k=0; k<K; k++){ // loop through centroids to 
				// find the closest one
					if (!allow_empty_clusters || numMembers[k]){
						// evaluate the distance; note the order of 
						// arguments: centroid followed by data --- it's 
						// important for non-symmetric distances
						double dist = (*distf)(centroids[k], data[n], 
							row_length, centroids_tilde[k]);
						if (new_assignment==-1 || dist < min_dist) {
							min_dist = dist;
							new_assignment = k;
						}
					}
				}
				// check for assignment change
				assignmentChanged = assignmentChanged || (assignment[n] != 
					new_assignment);
				// update assignment
				assignment[n] = new_assignment;
			}

			if (!assignmentChanged){
				// clustering succeded; compute distortion and return
				double distortion = 0;
				for (int n=0; n<N; n++) {
					distortion += distf(centroids[assignment[n]], data[n],
					row_length, centroids_tilde[assignment[n]]);
				}
				if (verbose > 0) {
					cerr<<"\tEM routine: attempt succeded"<<endl;
				}
				return distortion;
			}

			// evaluate the number of members in each cluster
			memset(numMembers, 0, K*sizeof(*numMembers));
			for(int n=0; n<N; n++){
				numMembers[assignment[n]] += 1;
			}

			// check for empty clusters
			if (!allow_empty_clusters){
				for(int k=0; k<K; k++){
					if(!numMembers[k]){
						// got an empty cluster; continue from the beginning
						goto restart;
					}
				}
			}

			// evaluate the new centroids
			for(int k=0; k<K; k++){
				if(!allow_empty_clusters || numMembers[k]){
					int elems_in_cluster = 0;
					for(int n=0; n<N; ++n){
						if(assignment[n] == k){
							tmp_data[elems_in_cluster] = data[n];
							tmp_data_1[elems_in_cluster] = data_1[n];
							elems_in_cluster++;
						}
					}
					eval_centroid(elems_in_cluster, row_length, tmp_data,
					num_nt, tmp_data_1, centroids[k], centroids_tilde[k]);
				}
			}
		}
	}
}


// Count number of distinct clusters from assignment vector.
// Cluster numbers are non-negative. They can be non-contiguous:
// e. g., 0, 2, 4, 5, 7
static int count_num_clusters(int const N, const int * assignment)
{
	int max_id = 0;
	for (int n=0; n<N; n++){
		if (assignment[n] > max_id){
			max_id = assignment[n];
		}
	}
	max_id++; // increment to get the size of the index vector
	int *ix = new int[max_id];
	memset(ix, 0, max_id*sizeof(*ix));
	for (int n=0; n<N; n++) {
		ix[assignment[n]] = 1;
	}
	int  num_clusters=0;
	for (int k=0; k<max_id; k++) {
		num_clusters += ix[k];
	}
	delete[] ix;

	return num_clusters; 
}


// implementation of a publicly accessible function
int hard_em(int K, int N, int row_length, double** data, int num_nt,
double** data_1, int* assignment, double *Z, int num_trials, char dist_type,
int verbose)
{
	// Allocate matrices for centroids, distances, assignment and the 
	// vector for the number of members
	double**  centroids = new double*[K];
	centroids[0] = new double[row_length*K];
	for(int k=1; k<K; ++k) {
		centroids[k] = centroids[0] + k*row_length;
	} // K*row_length matrix for centroid locations

	double**  centroids_tilde = new double*[K];
	centroids_tilde[0] = new double[row_length*K];
	for(int k=1; k<K; ++k) {
		centroids_tilde[k] = centroids_tilde[0] + k*row_length;
	} // K*row_length matrix for X_tilde vars for centroids

	int* tmp_assignment = new int[N];
	int *numMembers = new int[K];  // number of members in each cluster
	double** tmp_data = new double*[N];
	double** tmp_data_1 = new double*[N];

	// call clustering routine
	double min_distortion = 0;
	for (int t=0; t<num_trials; t++) {
		if (verbose > 0) {
			cerr<<"Calling EM routine, attempt "<<t+1<<endl;
		}
		double distortion =	em_routine(K, N, row_length, data, num_nt, data_1,
			tmp_assignment, numMembers, centroids, centroids_tilde,
			tmp_data, tmp_data_1, dist_type, verbose-1);
		if (verbose > 0) {
			cerr<<"Resulting distortion: "<<distortion<<endl;
		}
		if (!t || (distortion < min_distortion)) {
			if (verbose > 0) {
				cerr<<"Updating best partitioning and minimal distortion"<<endl;
			}
			min_distortion = distortion;
			memcpy(assignment, tmp_assignment, N*sizeof(*assignment));
		}
	}
	if (verbose > 0) {
		cerr<<"Resulting minimal distortion: "<<min_distortion<<endl;
	}

	//evaluate confidence if using KL with raw counts
	if (dist_type == 'k') {
		eval_confidence(K, N, row_length, data, centroids, Z);
	}

	// delete allocated arrays
	delete[] centroids[0];
	delete[] centroids;
	delete[] centroids_tilde[0];
	delete[] centroids_tilde;
	delete[] numMembers;
	delete[] tmp_data;
	delete[] tmp_data_1;
	delete[] tmp_assignment;

	// return the number of clusters
	return count_num_clusters(N, assignment);
}


// implementation of a publicly accessible function
int cons_em(int K, int N, int row_length, double** data, int num_nt,
double** data_1, int* assignment, int num_trials, char dist_type,
int verbose, double boot_frac, double cons_threshold)
{
	// Allocate matrices for centroids, distances, assignment and the 
	// vector for the number of members
	double**  centroids = new double*[K];
	centroids[0] = new double[row_length*K];
	for(int k=1; k<K; ++k) {
		centroids[k] = centroids[0] + k*row_length;
	} // K*row_length matrix for centroid locations

	double**  centroids_tilde = new double*[K];
	centroids_tilde[0] = new double[row_length*K];
	for(int k=1; k<K; ++k) {
		centroids_tilde[k] = centroids_tilde[0] + k*row_length;
	} // K*row_length matrix for X_tilde vars for centroids

	int *numMembers = new int[K];  // number of members in each cluster
	double** tmp_data = new double*[N];
	double** tmp_data_1 = new double*[N];
	double** boot_data = new double*[N];
	double** boot_data_1 = new double*[N];
	int* boot_list = new int[N];
	// CC = number of times the two samples appeared in the same bootstrap dataset
	// (upper triangular matrix)
	// it will be also used for heap later
	int** CC = new int*[N-1];
	for(int m=0; m<N-1; m++) {
		CC[m] = new int[N-m-1];
	}
	// CM = number of times the two samples were clustered together
	// (upper triangular matrix)	
	// we decalre it as double since it will be used 
	// as the distance matrix after post-processing
	double** CM = new double*[N-1];
	for(int m=0; m<N-1; m++) {
		CM[m] = new double[N-1-m];
	}
	// initialize the correlation matrices
	for(int m=0; m<N-1; m++) {
		memset(CC[m], 0, (N-1-m) * sizeof(**CC));
		memset(CM[m], 0, (N-1-m) * sizeof(**CM));
	}

	// perform EM clustering num_trials times
	for(int t=0; t<num_trials; t++) {
		// choose bootstrap samples
		if (verbose>0){
			cerr<<"Bootstrap dataset "<<t+1<<endl;
		}
		int B=0;  // number of samples in a bootstrap dataset
		for(int n=0; n<N; n++) {
			if (rand() <= boot_frac * RAND_MAX) {
				boot_list[B] = n;
				boot_data[B] = data[n];
				boot_data_1[B] = data_1[n];
				++B;
			}
		}
		if (verbose > 0) {
			cerr<<B<<" samples chosen for bootstrap dataset"<<endl;
		}
		if (B<2*K){ // not enough samples
			if (verbose>0){
				cerr<<"Not enough samples; skipping the dataset"<<endl;
			}
			continue;
		}
		// run EM routine on this bootstrap dataset
		if (verbose > 0) {
			cerr<<"Calling EM routine, attempt "<<t+1<<endl;
		}
		em_routine(K, B, row_length, boot_data, num_nt, boot_data_1,
		assignment, numMembers, centroids, centroids_tilde, tmp_data,
		tmp_data_1, dist_type, verbose-1);
		// update correlation matrices CC and CM
		for(int m=0; m<B; m++) {
			int p1 = boot_list[m];
			for(int n=m+1; n<B; n++) { 
				int p2 = boot_list[n] - p1 - 1; // boot_list is sorted => p2>=0
				CC[p1][p2] += 1;
				if (assignment[m]==assignment[n]) // m and n are the positions
				// in the boot_list --- em_routine operates with them, not
				// p1 and p2
				{
					CM[p1][p2] += 1;
				}
			}
		}
	}
	// evaluate the distance matrix
	for(int m=0; m<N-1; m++) {
		for(int n=0; n<N-1-m; n++) {
			// make sure denominator is non-zero
			double r = CM[m][n] 
				/ (double) ( (CC[m][n]>0) ? CC[m][n] : 1);
			// correlation distance
			CM[m][n] = 1-r;
		}
	}

	// now perform hierarchical clustering with the distance matrix
	hcluster(K, N, assignment, cons_threshold, CM, CC, verbose -1);

	// de-allocate memory
	delete[] centroids[0];
	delete[] centroids;
	delete[] centroids_tilde[0];
	delete[] centroids_tilde;

	delete[] numMembers;
	delete[] tmp_data;
	delete[] tmp_data_1;
	delete[] boot_data;
	delete[] boot_data_1;
	delete[] boot_list;
	for(int m=0; m<N-1; m++)
	{
		delete[] CC[m];
		delete[] CM[m];
	} 
	delete[] CC;	
	delete[] CM;

	// return the number of clusters
	return count_num_clusters(N, assignment);
}


// implementation of a publicly accessible function
void soft_em(int K, int N, int row_length, double** data, int num_nt,
double** data_1, double* Z, double threshold, int verbose)
{
	// total count of words in each sequence
	double* count = new double[N];
	memset(count, 0, N * sizeof(*count));
	for (int n=0; n<N; n++){
		for (int l=0; l<row_length; l++) {
			count[n] += data[n][l];
		}
	}
	
	// K*row_length matrices for centroid locations and K*N for probabilities
	double** centroids = new double*[K];
	double** old_centroids = new double*[K];
	for(int k=0; k<K; ++k) {
		*(centroids + k) = new double[row_length];
		*(old_centroids + k) = new double[row_length];
	}

	for(int k=0; k<K; ++k){
	}

	if (verbose > 0) {
		cerr<<"Soft EM routine: starting clustering"<<endl;
	}
	// initial centroid assignment
	for (int k=0; k<K; ++k) {
		memcpy(centroids[k], data[rand()% N], row_length * 
		sizeof(**centroids));
		// need to normalize the vector for KL divergence
		normalize_row(centroids[k], row_length);
	}

	// start iterations of soft EM clustering
	bool flag = false;  // this initial value is unused
	do {
		// retain centroid locations
		for(int k=0; k<K; ++k) {
			memcpy(old_centroids[k], centroids[k], row_length * 
			sizeof(**centroids));
		}
		// evaluate the new assignment probabilities
		eval_confidence(K, N, row_length, data, centroids, Z);

		// evaluate new centroids
		for (int k=0; k<K; k++) {
			memset(centroids[k], 0, row_length * sizeof(**centroids));
			for (int l=0; l<row_length; l++) {
				for (int n=0; n<N; n++) {
					centroids[k][l] += Z[k + n*K] * data[n][l];
				}
			}
			normalize_row(centroids[k], row_length);
		}

		// check that the shift of at least one of centroids 
		// is sufficiently large to warrant another iteration
		flag = false;
		for (int k=0; k<K; k++) {
			if (sqrt(euclidean_distance(centroids[k], old_centroids[k],
			row_length, NULL)) > threshold) {
				flag = true;
				if (verbose>1){
					cerr<<"\tsufficient centroid shift detected; "
					"need another iteration"<<endl;
				}
				break;
			}
		}
	} while (flag);

	// delete allocated objects
	delete[] count;
	for(int k=0; k<K; ++k) {
		delete[] centroids[k];
		delete[] old_centroids[k];
	}
	delete[] centroids;
	delete[] old_centroids;
	
	return;
}
