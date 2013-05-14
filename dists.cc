#include <string.h> // memcpy
#include <math.h>

double euclidean_distance(double *x, double *y, int L, double* xt)
{
	double S = 0;
	double z;
	for(int i=0; i<L; i++){
		z = x[i] - y[i];
		S += z*z;
	}
	return S;
}


double kl_distance(double *q, double *p, int L, double* qt)
{
	double total_count = 0;
	for(int i=0; i<L; ++i){
		total_count += p[i];
	}
	double S = 0;
	for(int i=0; i<L; i++){
		S += p[i] * log(p[i]/(total_count * q[i]));
	}
	return S;
}


double symkl_distance(double *q, double *p, int L, double* qt)
{
	double total_count = 0;
	for(int i=0; i<L; ++i){
		total_count += p[i];
	}
	double S = 0;
	for(int i=0; i<L; i++){
		S += (p[i]-total_count*q[i]) * log(p[i]/(total_count*q[i]));
	}
	return S;
}


double d2_distance(double *q, double *p, int L, double* qt)
{
	double norm = 0;
	double S = 0;
	for(int i=0; i<L; i++){
		norm += p[i] * p[i];
	}
	norm = sqrt(norm);
	for(int i=0; i<L; i++){
		S += q[i] * p[i] / norm;
	}
	return 1 - S;
}


double chi2_distance(double *q, double *p, int L, double* qt)
{
	double chi2 = 0;
	double total_count = 0;
	for(int i=0; i<L; i++){
		total_count += p[i];
	}
	for(int i=0; i<L; i++){
		double exp_count = q[i] * total_count;
		chi2 += (p[i] - exp_count) * (p[i] - exp_count) / exp_count;
	}
	//return (chi2 - (L-3) * log(chi2))/2;
	return chi2;
}


double d2ast_distance(double *q, double *p, int L, double* qt)
{
	double x[L];
	memcpy(x, p, L*sizeof(*p));
	double total_count = 0;
	for(int i=0; i<L; i++){
		total_count += p[i];
	}
	double S = 0;
	for(int l=0; l<L; ++l){
		x[l] -= total_count * q[l];
		x[l] /= sqrt(q[l]);
		S += x[l] * x[l];
	}
	S = sqrt(S);
	for(int l=0; l<L; ++l){
		x[l] /= S;
	}	
	S = 0;
	for(int l=0; l<L; ++l){
		S += x[l] * qt[l];
	}
	return (1-S)/2;
}
