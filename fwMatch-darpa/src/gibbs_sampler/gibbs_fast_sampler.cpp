// X = gibbs_fast_sampler(node_potentials, edge_potentials, node_count, sample_count, burnin, sampint);
// node_potentials: column vector with potentials of the nodes
// edge_potentials: full/dense matrix with the potentials of the edges (must be symmetric)
// node_count: number of nodes
// sample_count: number of samples desired (suggestion: 20000)
// burnin: number of burn in samples (suggestion: 10000)
// sampint: sample interval (1 for no gap) (suggestion : 100
// randomseed: just some int number so the algorithm is deterministic


# include <iostream>
# include <cmath>
# include "mex.h"
# include <ctime>

#define matrix_index(x,t,M) (((t)*(M)) + x)
using namespace std;

double flip(){
	double flipv = static_cast< double >( rand() ) / RAND_MAX;
	return(flipv);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	//START OF MAT INPUT
	if(nlhs < 1 || nlhs > 2) {
		cout<<"Function needs to return 1 or 2 arguments\n";
		return;
	}

	if(nrhs != 7) {
		cout<<"Function takes 7 arguments\n";
		return;
	}	
    
    // unpack parameters
	int node_count = (int)(mxGetScalar(prhs[2]));
	int sample_count = (int)(mxGetScalar(prhs[3]));
	int burnin = (int)(mxGetScalar(prhs[4]));
	int sampling_interval = (int)(mxGetScalar(prhs[5]));
    int random_seed = (int)(mxGetScalar(prhs[6]));
    
    // unpack nodes (potentials)
	double* nodes_ptr = mxGetPr(prhs[0]);
    double* node_potentials = new double[node_count];
    for(int i=0; i<node_count; i++)
        node_potentials[i] = nodes_ptr[i];
    
    // unpack edges (potentials)
	double* edges_ptr = mxGetPr(prhs[1]);
    double** edge_potentials = new double*[node_count];
    for(int i=0; i<node_count; i++)
		edge_potentials[i] = new double[node_count];
    for(int i=0; i<node_count; i++)
		for(int j=0; j<node_count; j++)
			edge_potentials[i][j] = (double) edges_ptr[matrix_index(i,j,node_count)];
	
	cout << "Running fast Gibbs Sampler" << endl;
	
	//END OF MAT INPUT
	
	// allocating sample matrix
	bool** X = new bool*[sample_count];
	for(int i=0; i<sample_count ;i++){
		X[i] = new bool[node_count];
	}
    
    // allocating pseudo log like vector
    unsigned int total_samples_taken = burnin + sample_count*sampling_interval;
    double* logLike = new double[total_samples_taken];
	
    // random initial value
    srand(random_seed);
    bool* x = new bool[node_count];
    for(int j=0; j<node_count; j++){
		x[j] = (flip() > 0.5);
	}
    
	int sampleSkip = 0;
	for(int i=0; i<sample_count; i++){
		if(i==0){
			sampleSkip = burnin + sampling_interval;
		} else{
			sampleSkip = sampling_interval;
		}
	
		for(int ctr=0; ctr<sampleSkip; ctr++){
			double currentSampleLogLike = 0;
            
            for(int j=0; j<node_count; j++){
				double log_odds = node_potentials[j];
				for(int k=0; k<node_count; k++){
					if(k != j){
                        log_odds = log_odds + edge_potentials[j][k] * x[k];
					}
				}
				double prob = exp(log_odds)/(exp(log_odds)+1);
                
				x[j] = (flip() <= prob);
                X[i][j] = x[j];
                
                // compute pseudo loglike (must accumulate for all nodes)
                if (x[j]) {
                    currentSampleLogLike += log(prob);
                } else {
                    currentSampleLogLike += log(1-prob);
                }
			}
            
            // save log like (must do some math to figure the index)
            if (i==0) {
                logLike[ctr] = currentSampleLogLike;
            } else {
                logLike[burnin + ctr + sampling_interval*i] = currentSampleLogLike;
            }
		}
	
	}
	
	plhs[0] = mxCreateDoubleMatrix(sample_count,node_count,mxREAL);
	double* X_ret = mxGetPr(plhs[0]);
	for(int l=0; l<sample_count; l++){
		for(int s=0; s<node_count; s++){
			X_ret[matrix_index(l,s,sample_count)] = X[l][s];
		}
	}
    
    plhs[1] = mxCreateDoubleMatrix(1, total_samples_taken,mxREAL);
	double* LogLike_ret = mxGetPr(plhs[1]);
	for(int i=0; i<total_samples_taken; i++){
        LogLike_ret[i] = logLike[i];
	}
}

