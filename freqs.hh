#pragma once

#include <string>
#include <vector>

using namespace std;

// number of characters
enum {NUM_NT = 4};

// normalize row dividing by the total counts
void normalize_row(double* fv, int L);

// build vector of frequencies of overlapping K-mers 
// from string seq
// normalize = whether to calculate raw or normalized 
// counts (default=normalized)
// pc = pseudocount
void fill_overlap_count_vector(string seq, int K, double* freq_vector, 
	int normalize=1, double pc=1);

// build vector of frequencies of non-overlapping K-mers
// from string seq
// normalize = whether to calculate raw or normalized 
// counts (default=normalized)
// pc = pseudocount
void fill_count_vector(string seq, int K, double* freq_vector, 
	int normalize=1, double pc=1);

// normalize frequency matrix to make its columns univariant
void normalize_freq_matrix(double** data, int N, int row_length);

