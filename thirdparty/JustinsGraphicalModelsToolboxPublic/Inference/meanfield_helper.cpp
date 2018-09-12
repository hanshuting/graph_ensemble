#include <stdio.h>
#include <map>
#include "mex.h"   //--This one is required
//#define nanoMatDebug 1 // 1 for safety/error checking; 0 for faster and unsafe(er)
//#include "nanoMatT.cpp"
//#include "etreeC.cpp"
#include "math.h"
//#define mmin(X,Y) ((X) < (Y) ? (X) : (Y))
//#define mmax(X,Y) ((X) > (Y) ? (X) : (Y))
#include "../eigenstuff.cpp"

double randd(){
 return double(rand())/RAND_MAX;
}

double randn(){
    // box-muller transform
    double u1     = randd();
    double u2     = randd();
    double R      = sqrt(-2*log(u1));
    double theta  = 6.283185307179586*u2;
    return R*cos(theta);
}


// double log_sum_exp(dMat x){
//     double C = -1*x.max(1)(0);
//     return log((x+C).expp().sum(1)(0))-C;
// }

// mcmc_helper(im,im_prev,buf,int(xwho),int(xwho2));
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    // trw_simple_helper(model,c_ij, psi_i, rho, n, m1, m2);
    
          
    int i=0;
    int nnodes     = mapdouble(mxGetField(prhs[i],0,"nnodes"))(0);
    int ncliques   = mapdouble(mxGetField(prhs[i],0,"ncliques"))(0);
    int nvals      = mapdouble(mxGetField(prhs[i],0,"nvals"))(0);
    MatrixXi pairs = mapdouble(mxGetField(prhs[i],0,"pairs")).cast<int>();
    pairs.array() -= 1;
    MatrixXi N1   = mapdouble(mxGetField(prhs[i],0,"N1")).cast<int>();
    N1.array() -= 1;
    MatrixXi N2   = mapdouble(mxGetField(prhs[i],0,"N2")).cast<int>();
    N2.array() -= 1;
    //vector<MatrixXi> N1;
    //vector<MatrixXi> N2;
    //mxArray *mxN1 = mxGetField(prhs[i],0,"N1");
    //mxArray *mxN2 = mxGetField(prhs[i],0,"N2");
    //for(int node=0; node<nnodes; node++){
    //    N1.push_back( mapdouble(mxGetCell(mxN1,node)).cast<int>() );
    //    N2.push_back( mapdouble(mxGetCell(mxN2,node)).cast<int>() );
    //    N1[node].cwise() -= 1;
    //    N2[node].cwise() -= 1;
    //}
    
    i++;
        
    MatrixMd psi_ij   = mapdouble(prhs[i++]);
    MatrixMd psi_i    = mapdouble(prhs[i++]);
    int maxiter       = mapdouble(prhs[i++])(0);
    MatrixMd b_i      = mapdouble(prhs[i++]);
    MatrixMd b_ij     = mapdouble(prhs[i++]);
    double convthresh = mapdouble(prhs[i++])(0);
    MatrixMd reps_out = mapdouble(prhs[i++]);
    
    MatrixXd b_i_save = b_i;
    
    int reps;
    double conv;
    
    for(reps=0; reps<maxiter; reps++){
        b_i_save = b_i;
    
        //for(int j=0; j<nnodes; j++){
        for(int j0=0; j0<nnodes*2-1; j0++){
            int j, mode;
            if(j0<nnodes)
                j    = j0;
            else
                j    = nnodes - 2 - (j0-nnodes);

            for(int yj=0;yj<nvals; yj++){
                b_i(yj,j) = psi_i(yj,j);
                for(int k=0; k<N1.cols(); k++){
                    int d = N1(j, k);
                    if(d==-2) continue;
                    i = pairs(d, 1);
                    for(int yi=0; yi<nvals; yi++){
                        int index = yj + yi*nvals;
                        b_i(yj, j) *= pow(psi_ij(index, d), b_i(yi, i));
                    }
                }
                for(int k=0; k<N1.cols(); k++){
                    int d = N2(j, k);
                    if(d==-2) continue;            
                    i = pairs(d,0);
                    for(int yi=0; yi<nvals; yi++){
                        int index = yi + yj*nvals;
                        b_i(yj,j) = b_i(yj,j)*pow(psi_ij(index,d),b_i(yi,i));
                    }
                }
            }
           b_i.col(j) = b_i.col(j) / b_i.col(j).sum();
        }
        conv = (b_i-b_i_save).array().abs().maxCoeff();
        //fprintf('iter: %d  conv: %f \n', iter, conv);
        if(conv < convthresh){
            break;
        }
    }
    reps_out(0) = reps;
}

