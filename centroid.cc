#include <string.h> // memset
#include <math.h>

// mean centroid
void mean_centroid(int N, int row_length, double** data, int num_nt,
double** data_1, double* centroid, double* centroid_tilde)
{
	memset(centroid, 0, row_length * sizeof(*centroid));
	for(int n=0; n<N; n++){
		for(int l=0; l<row_length; l++){
			centroid[l] += data[n][l];
		}
	}
	for(int l=0; l<row_length; l++){
		centroid[l] /= N;
	}
	return;
}


// d2 centroid
void d2_centroid(int N, int row_length, double** data, int num_nt,
double** data_1, double* centroid, double* centroid_tilde)
{
	// add the total counts
	memset(centroid, 0, row_length * sizeof(*centroid));
	for(int n=0; n<N; ++n){
		for(int l=0; l<row_length; ++l){
			centroid[l] += data[n][l];
		}
	}
	// now normalize the entries
	double total_sq_count = 0;
	for(int l=0; l<row_length; ++l){
		total_sq_count += centroid[l] * centroid[l];
	}
	double norm = sqrt(total_sq_count);
	for(int l=0; l<row_length; ++l){
		centroid[l] /= norm;
	}
	return;
}


// KL centroid --- operates on raw counts which do not have to be normalized
void kl_centroid(int N, int row_length, double** data, int num_nt,
double** data_1, double* centroid, double* centroid_tilde)
{
	// add the total counts
	memset(centroid, 0, row_length * sizeof(*centroid));
	for(int n=0; n<N; ++n){
		for(int l=0; l<row_length; ++l){
			centroid[l] += data[n][l];
		}
	}
	// now normalize the entries
	double total_count = 0;
	for(int l=0; l<row_length; ++l){
		total_count += centroid[l];
	}
	for(int l=0; l<row_length; ++l){
		centroid[l] /= total_count;
	}
	return;
}


// MM centroid --- evaluates frequencies of single nucleotides
// and coumputes frequencies of words using zero order Markov model;
// i.e., as product of single nucleotide frequencies:
// f_{w1...w_k} = f_w1 f_w2 ... f_wk
void mm_centroid(int N, int row_length, double** data, const int num_nt, 
double** data_1, double* centroid, double* centroid_tilde)
{
	double nt_freq[num_nt];
	memset(nt_freq, 0, num_nt * sizeof(nt_freq[0]));
	double S=0;
	// compute frequencies of individual nucleotides
	for (int n=0; n<N; n++){
		for (int l=0; l<num_nt; l++){
			nt_freq[l] += data_1[n][l];
		}
	}
	for (int l=0; l<num_nt; l++){
		S += nt_freq[l];
	}
	for (int l=0; l<num_nt; l++){
		nt_freq[l] /= S;
	}
	// compute frequencies of words using Markov model
	for(int l=0; l<row_length; l++){
		double p = 1;
		for(int M=1; M<row_length; M*=num_nt){
			int digit = l/M % num_nt;
			p *= nt_freq[digit];
		}
		centroid[l] = p;
	}
}


// d2* centroid
void d2ast_centroid(int N, int row_length, double** data, int num_nt,
double** data_1, double* centroid, double* centroid_tilde)
{
	// compute frequencies from MM
	mm_centroid(N, row_length, data, num_nt, data_1, centroid, NULL);
	// compute X_tilde
	memset(centroid_tilde, 0, row_length * sizeof(*centroid_tilde));
	for(int n=0; n<N; ++n){
		for(int l=0; l<row_length; ++l){
			centroid_tilde[l] += data[n][l];
		}
	}
	double total_count = 0;
	for(int l=0; l<row_length; ++l){
		total_count += centroid_tilde[l];
	}
	double S = 0;
	for(int l=0; l<row_length; ++l){
		centroid_tilde[l] -= total_count * centroid[l];
		centroid_tilde[l] /= sqrt(centroid[l]);
		S += centroid_tilde[l] * centroid_tilde[l];
	}
	S = sqrt(S);
	for(int l=0; l<row_length; ++l){
		centroid_tilde[l] /= S;
	}
	return;
}
