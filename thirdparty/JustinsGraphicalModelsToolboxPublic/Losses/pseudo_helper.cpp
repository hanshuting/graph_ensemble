#include <stdio.h>
#include <map>
#include "mex.h"   //--This one is required
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
    
    // need from model struct:
    //model =
    //
    //  nnodes: 10000
    //ncliques: 20000
    //   nvals: 2
    //   pairs: [20000x2 double]
    //      N1: {10000x1 cell}
    //      N2: {10000x1 cell}
          
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
    i++;
        
    MatrixMd psi_ij  = mapdouble(prhs[i++]);
    MatrixMd psi_i   = mapdouble(prhs[i++]);
    MatrixXi x       = mapdouble(prhs[i++]).cast<int>();
    x.array() -= 1;
    MatrixMd L       = mapdouble(prhs[i++]);
    MatrixMd dpsi_ij = mapdouble(prhs[i++]);
    MatrixMd dpsi_i  = mapdouble(prhs[i++]);
    
    MatrixXd logb = MatrixXd::Zero(nvals,1);
    MatrixXd b    = MatrixXd::Zero(nvals,1);
    for(int i=0; i<nnodes; i++){
        if(x(i)==-1) continue;
        bool bad = false;
        for(int xi=0; xi<nvals; xi++){
            logb(xi) = psi_i(xi,i);
            for(int cwhere=0; cwhere<N1.cols(); cwhere++){
                int c = N1(i,cwhere);
                if(c==-2) break;
                int j = pairs(c,1);
                int xj = x(j);
                if(xj==-1){ bad=true; continue;}
                int where = xi + xj*nvals;
                logb(xi) += psi_ij(where,c);
            }
            for(int cwhere=0; cwhere<N2.cols(); cwhere++){
                int c = N2(i,cwhere);
                if(c==-2) break;
                int j = pairs(c,0);
                int xj = x(j);
                if(xj==-1){ bad=true; continue;}
                int where = xj + xi*nvals;
                logb(xi) += psi_ij(where,c);
            }
        }
        if(bad) continue;
        double A = log_sum_exp(logb,nvals);
        L(0) += -logb(x(i)) + A;
        
        for(int xi=0; xi<nvals; xi++){
            b(xi) = exp(logb(xi) - A);
            //cout << b(xi) << endl;
            
            double mult = double(xi==x(i)) - b(xi);
            
            dpsi_i(xi, i) -= mult;
            for(int cwhere=0; cwhere<N1.cols(); cwhere++){
                int c = N1(i, cwhere);
                if(c==-2) break;
                int j = pairs(c,1);
                int xj = x(j);
                int where = xi + xj*nvals;
                
                dpsi_ij(where, c) -= mult;
            }
            for(int cwhere=0; cwhere<N2.cols(); cwhere++){
                int c = N2(i, cwhere);
                if(c==-2) break;
                int j = pairs(c,0);
                int xj = x(j);
                int where = xj + xi*nvals;
                
                dpsi_ij(where, c) -= mult;
            }
        }
    }    
}

