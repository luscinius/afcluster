#include <assert.h>
#include <iostream>

using namespace std;

/** struct to represent the node in the hierarchical clustering
 * a node contains number of elements (1 in the beginning of agglomerative
 * clustering) and the one it was joined to (itself in the beginning)
 */

struct HCNode
{
    int num_members; // number of elements in the cluster
    int joined_to;  // node to which the current one was merged
};


/** Functions to work with the heap represented as an array.
 * In the array notation descendatns of node with index i 
 * are those with index 2i+1 and 2i+2; ancestor of node i is (i-1)/2.
 * Heap is supplemented with the index which allows one to find any 
 * element within the heap in constant time:
 * idx[n] = position of element with value n in the array.
 * This is needed for the logarithmic time removal of any element
 * (not only the minimal one).
 * Elements are compared using the corresponding distance matrix row.
 */

// function for swapping elements of an array and updating index
static inline void swap_items(int pos1, int pos2, int* heap, int* idx)
{
	int tmp = heap[pos1];
	heap[pos1] = heap[pos2];
	heap[pos2] = tmp;
	idx[heap[pos1]] = pos1;
	idx[heap[pos2]] = pos2;
	return;
}


// move the element at pos to the correct location in the heap,
// updating heap and index
static void heapify(int pos, const int heap_len, double * const data, int *heap,
int *idx)
{
	int anc = (pos-1)/2;
	if (pos>0 && data[heap[pos]] < data[heap[anc]]){ // move up as needed
		while(true){
			swap_items(pos, anc, heap, idx);
			pos = anc;
			anc = (pos-1)/2; // ancestor
			if (!pos || data[heap[pos]] >= data[heap[anc]]){
				break;
			}
		}
	} else { // move down as needed
		while (true) {
			int left = 2*pos + 1;
			int right = 2*pos + 2;
			int min = pos;
			if (left < heap_len && data[heap[left]] < data[heap[min]]){
				min = left;
			}
			if  (right < heap_len && data[heap[right]] < data[heap[min]]){
				min = right;
			}
			if (min == pos){ // no need to move further down
				break;
			} else{
				swap_items(pos, min, heap, idx);
				pos=min;
			}
		}
	}
    return;
}


// build min-heap of distances in data[], store it in heap[]; build index
// we need index for logarithmic time removal of arbitrary element
static int build_heap(const int N, const int n, HCNode * const nodes, 
double * const data, int *heap, int* idx)
{
	int heap_length=0;
	// add elements one by one and move them to the correct position
	for(int m=0; m<N-n-1; m++){
		if ((nodes+n+m+1)->joined_to==(n+m+1)){
			heap[heap_length] = m;
			idx[m] = heap_length;
			heapify(heap_length, heap_length+1, data, heap, idx);
			++heap_length;
		}
	}
	return heap_length;
}


// delete an element from the heap
// pos = position of element to be deleted
static void delete_element(int pos, int& heap_len, double * const data, 
int* heap, int* idx)
{
	// decrement heap array length and move the last element to pos; 
	// then heapify
	heap[pos] = heap[--heap_len]; 
	idx[heap[pos]] = pos; 
	heapify(pos, heap_len, data, heap, idx);
	return;
}


/** Function for distance martix update
 * pos1, pos2 = nodes to be joined
 * N1, N2 = number of members in their clusters
 * distmatrix = distance matrix to be updated
 * nodes = nodes for hierarchical clustering
 * N = original size of the distance matrix
 * 
 * Note that the node nd2 is discarded, therefore it requires no update
 */
static void update_distmatrix(int nd1, int nd2, int N1, int N2, 
double **distmatrix, HCNode *nodes, int N)
{
	assert(nd1 < nd2);
	for(int n=0; n<nd1; n++){ // update D[n, nd1-n-1]; n<nd1<nd2
		if((nodes+n)->joined_to == n){
			distmatrix[n][nd1-1-n] = (N1*distmatrix[n][nd1-1-n] 
			+ N2*distmatrix[n][nd2-1-n]) / (N1 + N2);
		}
	}
	for(int n=nd1+1; n<nd2; n++){ // update D[nd1, n-nd1-1]; nd1<n<nd2
		if((nodes+n)->joined_to == n){ 
			distmatrix[nd1][n-nd1-1] = (N1 * distmatrix[nd1][n-nd1-1]
			+ N2 * distmatrix[n][nd2-n-1]) / (N1 + N2);
		}
	}
	for(int n=nd2+1; n<N; n++){ // update D[nd1, n-nd1-1]; n>nd2
		if((nodes+n)->joined_to == n){ 
			distmatrix[nd1][n-nd1-1] = (N1 * distmatrix[nd1][n-nd1-1]
			+ N2 * distmatrix[nd2][n-nd2-1]) / (N1 + N2);
		}
	}
	return;
}


// publicly accessible function
int hcluster(int K, int N, int* assignment, double dist_cutoff, 
double** distmatrix, int** heap, int verbose)
{
	// initialize the array of nodes
	HCNode *nodes = new HCNode[N];
	for (int n=0; n<N; n++){
		(nodes + n)->num_members = 1;
		(nodes + n)->joined_to = n;
	}

	// initialize the heap for the distances
	int *heap_len = new int[N-1]; // length of distance heap for each element
	int **idx = new int*[N-1];
	for (int n=0; n<N-1; n++){
		idx[n] = new int[N-1-n];
	}
	for (int n=0; n<N-1; n++){
		heap_len[n] = build_heap(N, n, nodes, distmatrix[n], heap[n], idx[n]);
	}
	if(verbose>0){
		cerr<<"Starting hierarchical clustering"<<endl;
	}

	// perform hierarchical clustering
	int num_clusters;
	for(num_clusters=N; num_clusters>K; num_clusters--){
		// find the closest pair of nodes
		int nd1 = 0;
		int nd2 = nd1 + 1 + heap[nd1][0];
		double Dmin = distmatrix[nd1][heap[nd1][0]];
		for (int n=1; n<N-1; n++){
			if((nodes+n)->joined_to == n && heap_len[n] && 
			distmatrix[n][heap[n][0]] < Dmin){
				Dmin = distmatrix[n][heap[n][0]];
				nd1 = n;
				nd2 = nd1 + 1 + heap[nd1][0];
			}
		}
		if (Dmin > dist_cutoff){
			break;
		}
		// merge nodes
		if (verbose>0){
			cerr<<"Merging node "<<nd2<<" into "<<nd1<<" with distance "
			<<Dmin<<" at step "<<N+1-num_clusters<<endl;
		}
		// update node info
		int N1 = (nodes + nd1)->num_members;
		int N2 = (nodes + nd2)->num_members;
		(nodes+nd1)->num_members = N1 + N2;
		(nodes+nd2)->joined_to = nd1;
		// remove entries for nd2 (being merged into nd1) from the heap
		for(int n=0; n<nd2; n++){
			if((nodes+n)->joined_to==n){ // only unmerged nodes are relevant
				delete_element(idx[n][nd2-1-n], heap_len[n], distmatrix[n], 
				heap[n], idx[n]);
			}
		}
		// update the distance matrix --- do it only AFTER the removal of the 
		// entries for nd2 to avoid breaking the heap property in > 1 place
		update_distmatrix(nd1, nd2, N1, N2, distmatrix, nodes, N);
		// update the position for the entries for nd1 within the heap
		// do this after the removal of the entries for nd2 and 
		// the distance matrix update
		for(int n=0; n<nd1; n++){
			if((nodes+n)->joined_to==n){
				heapify(idx[n][nd1-1-n], heap_len[n], distmatrix[n], heap[n], 
				idx[n]);
			}
		}
		//rebuild the heap for nd1
		heap_len[nd1] = build_heap(N, nd1, nodes, distmatrix[nd1], heap[nd1], 
		idx[nd1]);
	}

	// post-process the data: identify assignment
	int largest_assignment = 0;
	for (int n=0; n<N; n++){
		int anc = (nodes+n)->joined_to;
		if(anc==n){
			assignment[n] = largest_assignment++;
		} else {
			assignment[n] = assignment[anc];
		}
	}

	// destroy the array of nodes
	delete[] nodes;
	// destroy heap related data
	for (int n=0; n<N-1; n++){
		delete[] idx[n];
	}
	delete[] idx;
	delete[] heap_len;

	return num_clusters;
}
