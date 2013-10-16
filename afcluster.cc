#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include "seqio.hh"
#include "freqs.hh"
#include "em.hh"

enum {PC=1};  // default for pseudocount
enum {TRIALS=100};  // number of bootstrap datasets for consensus clustering
const double SHIFT = 1e-5;  // threshold for centroid shift

void print_usage(char* progname)
{
	cout<<"Centroid based (k-means-like) clustering of sequences "
		"in n-mer frequency space\n\n"
	"usage: "<<progname<<" [-h] [-C] [-c num_clusters] [-d dist_type] "
	"[-f fraction] [-k nmer_length] [-m threshold] [-n] [-o] "
	"[-p pseudocount] [-q threshold] [-r] [-s] [-S seed] "
	"[-t num_trials] [-u] [-v] [-w] fasta_file\n\n"
	"\t-C: do consensus clustering\n"
	"\t-c num_clusters (5 clusters by default)\n"
	"\t-d dist_type:\n"
		"\t\ta: d2* distance\n"
		"\t\tc: chi square statistic\n"
		"\t\td: d2 distance\n"
		"\t\te: regular euclidean (L2) distance; default\n"
		"\t\tk: Kullback-Leibler divergence\n"
		"\t\ts: symmetrized Kullback-Leibler divergence\n"
	"\t-f fraction: perform consensus clustering using bootstrap "
		"datasets with the specified fraction of samples; implies "
		"consensus clustering\n"
	"\t-h print this message and exit\n"
	"\t-k nmer_length: length of word (2-mers by default)\n"
	"\t-m threshold: for consensus hierarchical clustering\n"
	"\t-n: normalize frequency matrix to make "
		"each column univariant; implies L2 distance\n"
	"\t-o: count overlapping k-mers (non-overlapping by default)\n"
	"\t-p pseudocount: (default: "<<PC<<"); can be fractional\n"
	"\t-q threshold: threshold for centroid shift (default: "<<SHIFT<<
		"); implies soft EM\n"
	"\t-r: reverse complement and stack together\n"
	"\t-s: stack frequencies of k-mers with different k together\n"
	"\t-S seed: initial seed for random number generator\n"
	"\t-t num_trials: repeat clustering num_trials; defaults: "
		"1 for regular clustering, "<<TRIALS<<" for consensus clustering. "
		"For regular clustering choose the best partitioning, for consensus "
		"clustering build the consensus partitioning. "
		"Option ignored for soft EM\n"
	"\t-u: use soft EM algorithm; resets -d to 'k' and turns off -n and -w\n"
	"\t-v: output progress messages (repeat for increased verbosity)\n"
	"\t-w: write sequences from each cluster to a file"<<endl;
	return;
}


int main(int argc, char **argv)
{
	bool consensus_flag = false; // consensus clustering; off by default
	int num_clusters = 5;  // target number of clusters
	char dist_type = 'e'; // default to L2
	double boot_frac = 0.80; // fraction of samples for bootstrap dataset
	double cons_threshold = 1.5; // threshold for consensus agreement
		// this value is a priori higher than 1; effectively it means
		// no threshold
	int K = 2;  // k-mer size
	bool normalize_flag = false; // make columns univariant ("whiten");
		// off by default
	void (*fill_freq_vector)(string, int, double*, int, double) 
		= &fill_count_vector; // function to count k-mers;
		// do not use overlap count by default
	double pseudocount = PC; // pseudocount for k-mers
	double shift_threshold = SHIFT; // threshold for centroid shift
	bool rc_flag = false; // append reverse complement to each sequence;
		// off by default
	bool stack_all_flag = false; // stack frequencies of all k-mers up to K;
		// off by default
	unsigned int random_seed = time(NULL); // initial seed
	int num_trials = 0; // later we check if it was specified using "-t" arg
		// and set the correct default if it was not; note that the default 
		// depends on whether consensus clustering is performed
	int verbose_level = 0;
	bool soft_em_flag = false; // use soft EM; off by default
	bool write_flag = false; // write sequences by cluster; off by default

	int opt;
	while((opt = getopt (argc, argv, "Cc:d:f:k:m:nop:q:rsS:t:uvwh")) != -1 ){
		switch(opt){
			case 'C':
				consensus_flag = true;
				break;
            case 'c':
                num_clusters = atoi(optarg);
                break;
			case 'd':
                dist_type = *optarg;
                break;
			case 'f':
				consensus_flag = true;
				boot_frac = atof(optarg);
				break;
            case 'k':
                K = atoi(optarg);
                break;
			case 'm':
				consensus_flag = true;
				cons_threshold = atof(optarg);
				break;
            case 'n':
                normalize_flag = true;
                break;
            case 'o':
                fill_freq_vector = &fill_overlap_count_vector;
                break;
			case 'p':
				pseudocount = atof(optarg);
				break;
			case 'q':
				soft_em_flag = true;
				shift_threshold = atof(optarg);
				break;
            case 'r':
                rc_flag = true;
                break;
            case 's':
                stack_all_flag = true;
                break;
			case 'S':
				random_seed = atoi(optarg);
				break;
            case 't':
                num_trials = atoi(optarg);
                break;
			case 'u':
				soft_em_flag = true;
				break;
            case 'v':
                ++verbose_level;
                break;
            case 'w':
                write_flag = true;
                break;
            case 'h':
				print_usage(argv[0]);
				return EXIT_SUCCESS;
            case '?':
                print_usage(argv[0]);
                return EXIT_FAILURE;
            default:
                // You won't actually get here
                break;
        }
	}

	// cannot request consensus and soft EM sumultaneously
	if (consensus_flag && soft_em_flag){
		cout<<"Cannot specify '-C' and '-u' simultaneously"<<endl;
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	// if "-t" argument was absent, use the appropriate default
	if (!num_trials) {
		if (consensus_flag) {
			num_trials = TRIALS;
		} else {
			num_trials = 1;
		}
	}
	
	// turn off row and column normalization for soft EM
	if (soft_em_flag){
		if (dist_type != 'k' && dist_type != 'c') {
			cout<<"Can only use soft EM with chi squared or KL distance"<<endl;
			print_usage(argv[0]);
			return EXIT_FAILURE;
		}
		normalize_flag = false;
	}

	if(argc != optind + 1){
		cout<<"Missing/extra input file name"<<endl;
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}
	char* fasta_file_name = *(argv + optind );

	// initialize random number generator
	srand(random_seed);

	// initialize constants
	int L = 1; // length of the freq. vector L = NUM_NT**K
	for(int k=0; k<K; ++k){
		L *= NUM_NT;
	}
	// when stacking frequencies of different k-mers one needs to adjust 
	// L so that L = NUM_NT + ... + NUM_NT**k 
	// = (NUM_NT**(k+1) -NUM_NT) / (NUM_NT -1).
	if (stack_all_flag){
		L = (L - 1) * NUM_NT / (NUM_NT -1);
	}

	// count sequences
	int N = CountReads(fasta_file_name);  // number of sequences
	if (verbose_level>0){
		cerr<<"Found "<<N<<" sequences\n";
	}

	// read sequences and build count matrix
	RecordGenerator rec_gen(fasta_file_name);
	double **freq = new double*[N];
	freq[0] = new double[N * L];
	double **freq_1  = new double*[N];
	freq_1[0] = new double[N * NUM_NT];
	for(int i=0; i<N; ++i){
		SeqRecord rec = rec_gen.next();
		string seq = rc_flag ? rec.seq() + rec.rc().seq() : rec.seq();
		freq[i] = freq[0] + i*L;
		if (stack_all_flag){
			// stack normalized counts 
			int start_pos = 0;
			for(int k=0; k<K; ++k){
				fill_freq_vector(seq, k+1,
					freq[i] + start_pos, (dist_type!='r'), pseudocount);
				start_pos = (start_pos + 1) * NUM_NT;
			}
		} else{
			// compute regular frequencies; normalize row unless dist_type=='r'
			fill_freq_vector(seq, K, *(freq+i), 
				(dist_type=='e'), pseudocount);
		}
		freq_1[i] = freq_1[0] + i*NUM_NT;
		fill_freq_vector(seq, 1, freq_1[i], false, pseudocount);
	}

	// whiten if requested; in this case reset the distance to euclidean L2
	if (normalize_flag){
		dist_type = 'e';
		normalize_freq_matrix(freq, N, L);
	}
	if (verbose_level > 0){
		cerr<<"Word frequencies calculated\n";
	}

	// soft EM clustering
	if (soft_em_flag){
		// allocate memory for probabilities
		double* Z = new double[num_clusters * N];
		// run clustering
		soft_em(num_clusters, N, L, freq, NUM_NT, freq_1, Z, 
		shift_threshold, verbose_level);
		// output the result
		RecordGenerator rec_gen(fasta_file_name);
		for (int n=0; n<N; n++) {
			cout<<rec_gen.next().id();
			for (int k=0; k<num_clusters; k++){
				cout.setf(ios_base::fixed);
				cout.precision(4);
				cout<<"\t"<<Z[k + n * num_clusters];
			}
			cout<<endl;
		}

		// de-allocate memory
		delete[] Z;

		return EXIT_SUCCESS;
	}

	// hard EM clustering
	int* assignment = new int[N];
	double* Z = new double[num_clusters * N];
	if (consensus_flag){
		if (verbose_level > 0) {
			cerr<<"Running EM clustering "<<num_trials<<" times "
			"to build the consensus partitioning"<<endl;
		}
		num_clusters = cons_em(num_clusters, N, L, freq, NUM_NT, freq_1,
		assignment, num_trials, dist_type, verbose_level-1, boot_frac,
		cons_threshold);
	} else{
		if (verbose_level > 0) {
			cerr<<"Running regular EM clustering "<<num_trials<<" times "
			"to chose the best partitioning"<<endl;
		}
		num_clusters = hard_em(num_clusters, N, L, freq, NUM_NT, freq_1,
		assignment, Z, num_trials, dist_type, verbose_level-1);
	}
	if (verbose_level > 0){
		cerr<<"Clustering done\n";
	}

	// release memory previously allocated for frequency data
	delete[] freq[0];
	delete[] freq;
	delete[] freq_1[0];
	delete[] freq_1;

	// reassign cluster names so that they are sorted 
	// in the descending order of the number of members
	int* num_members = new int[num_clusters];
	int* ix = new int[num_clusters];
	int* rix = new int[num_clusters];
	for (int i=0; i<num_clusters; ++i){
		*(num_members + i) = 0;
		*(ix + i) = i;
	}
	for (int i=0; i<N; ++i){
		num_members[assignment[i]]++;
	}
	for (int i=0; i<num_clusters; ++i){
		int max_pos = i;
		for (int j=i+1; j<num_clusters; ++j){
			if (num_members[ix[j]] > num_members[ix[max_pos]]) {
				max_pos = j;
			}
		}
		int tmp = ix[i];
		ix[i] = ix[max_pos];
		ix[max_pos] = tmp;
	}
	for (int i=0; i<num_clusters; ++i){
		rix[ix[i]] = i;
	}
	for (int i=0; i<N; ++i) {
		*(assignment + i) = *(rix + *(assignment + i));
	}
	delete[] num_members;
	//delete[] ix and rix later

	cout.setf(ios_base::fixed);
	cout.precision(4);

	// output assignment to stdout
	RecordGenerator rec_gen_2(fasta_file_name);
	for(int n=0; n<N; ++n){
		cout<<rec_gen_2.next().id()<<"\t"<<*(assignment+n);
		if (!consensus_flag && dist_type=='r' && num_trials==1){
			for (int k=0; k<num_clusters; k++){
				cout<<"\t"<<Z[ix[k] + n * num_clusters];
			}
		}
		cout<<endl;
	}
	delete[] ix;
	delete[] rix;

	// write sequences to individual files if requested
	if (write_flag) {
		vector<SeqRecord> rec_vec;
		FastaRead(fasta_file_name, rec_vec);
		int N = rec_vec.size();  
		char fname[256];
		for (int i=0; i<num_clusters; ++i){
			sprintf(fname, "part_%i.fna", i+1);
			ofstream ofs(fname);
			for (int j=0; j<N; ++j) {
				if (*(assignment + j) == i) {
					ofs<<">"<<rec_vec[j].desc()<<endl;
					ofs<<rec_vec[j].seq()<<endl;
				}
			}
		}
	}

	// free allocated memory
	delete[] Z;
	delete[] assignment;
	return EXIT_SUCCESS;
}

