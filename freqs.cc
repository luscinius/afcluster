#include <math.h>
#include "freqs.hh"
#include <iostream>


static bool is_valid_nt(char c)
{
	c = toupper(c);
	if (c == 'A') {
		return true;
	} 
	if (c == 'C') {
		return true;
	} 
	if (c == 'G') {
		return true;
	} 
	if (c == 'T') {
		return true;
	}
	return false;
}


static int nt2int(char c)
{
	switch (c) {
		case 'A':
		case 'a': 
			return 0;
		case 'C':
		case 'c': 
			return 1;
		case 'G':
		case 'g': 
			return 2;
		case 'T':
		case 't': 
			return 3;
	}
	// should never get here
	return -1;
}

void normalize_row(double* fv, int L)
{
	double total = 0;
	for (int l=0; l<L; l++){
		total += fv[l];
	}
	for (int l=0; l<L; l++){
		fv[l] /= total;
	}
	return;
}


// calculate vector of counts of different polynucleotides with overlaps
// pc=pseudocount; set it to 0 to disable
void fill_overlap_count_vector(string seq, int K, double* fv, int normalize,
	double pc)
{
	int L = seq.length();
	for(int i=0; i<L; ++i){
		seq[i] = toupper(seq[i]);
	}
	int N = 1;  // number of K-mers
	for(int k=0; k<K; ++k){
		N *= NUM_NT;
	}
	int MODULUS = N/NUM_NT;
	for(int i=0; i<N; ++i){
		*(fv+i) = pc;  // initialize with the pseudocount
	}
	// fill the vector
	int pos = 0;
	bool in_progress = false; // set to true when extending a valid word
	int ix = 0; // integer corresponding to the current word
	while (pos<L){
		if (!in_progress){
			if (pos>L-K){
				break; // not enough characters left to form a K-mer
			}
			ix = 0; 
			in_progress = true;
			for(int i=pos; (i<L) && (i<pos+K); ++i){
				if (!is_valid_nt(seq[i])){
					// invalid character
					pos = i+1;
					in_progress = false;
					break;
				}	
				ix = NUM_NT*ix + nt2int(seq[i]);
			}
			if (in_progress){
				*(fv+ix) += 1;
				pos += K;
			}
		} else { // if (!in_progress)
			if (!is_valid_nt(seq[pos])){
				in_progress = false;
				pos += 1;
				continue;
			}
			ix = NUM_NT * (ix % MODULUS) + nt2int(seq[pos]);
			*(fv+ix) += 1;
			pos += 1;
		}
	}
	// normalize the counts to obtain frequencies
	if (normalize){
		normalize_row(fv, N);
	}
	return;
}


// calculate vector of counts of different polynucleotides without overlaps
// pc=pseudocount; set it to 0 to disable
void fill_count_vector(string seq, int K, double* fv, int normalize,
	double pc)
{
	int L = seq.length();
	for(int i=0; i<L; ++i){
		seq[i] = toupper(seq[i]);
	}
	int N = 1;  // number of K-mers
	for(int k=0; k<K; ++k){
		N *= NUM_NT;
	}
	for(int i=0; i<N; ++i){
		*(fv+i) = pc;  // initialize with the pseudocount
	}
	// fill the vector
	for(int i=0; i<=L-K; i+=K){
		int j;
		int ix=0;
		for(j=0; j<K; ++j){
			if (!is_valid_nt(seq[i+j])){
				break;
			}
			ix = ix*NUM_NT + nt2int(seq[i+j]);
		}
		if (j==K){
			*(fv+ix) += 1;
		}
	}
	// normalize if requested
	if (normalize){
		normalize_row(fv, N);
	}
	return;
}



// normalize frequency matrix to make its columns univariant
void normalize_freq_matrix(double** data, int N, int row_length)
{
	double tmp, SQ, S, V;
	for(int l=0; l<row_length; ++l){
		SQ=0;
		S=0;
		for(int n=0; n<N; ++n){
			tmp = *(*(data+n)+l);
			SQ += tmp*tmp;
			S += tmp;
		}
		V = sqrt( SQ/N - (S*S)/(N*N) );  // variance
		if (V>0){
			for(int n=0; n<N; ++n){
				*(*(data+n)+l) /= V;
			}
		}
	}
	return;
}
