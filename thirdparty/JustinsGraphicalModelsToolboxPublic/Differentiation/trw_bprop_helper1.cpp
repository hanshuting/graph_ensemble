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
    // trw_bprop_helper1(model, psi_ij_rho, psi_i, rho, maxiter, n, m1, m2, b_i);
          
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
        
    MatrixMd psi_ij = mapdouble(prhs[i++]);
    MatrixMd psi_i  = mapdouble(prhs[i++]);
    double rho      = mapdouble(prhs[i++])(0);
    int maxiter     = mapdouble(prhs[i++])(0);
    //double damp     = mapdouble(prhs[i++])(0);
    double convthresh = mapdouble(prhs[i++])(0);
    MatrixMd n      = mapdouble(prhs[i++]);
    MatrixMd m1     = mapdouble(prhs[i++]);
    MatrixMd m2     = mapdouble(prhs[i++]);
    MatrixMd b_i    = mapdouble(prhs[i++]);
    //MatrixMd b_ij   = mapdouble(prhs[i++]);
    //MatrixMd Z_ij   = mapdouble(prhs[i++]);
    MatrixMd mstor  = mapdouble(prhs[i++]);
    MatrixMd b_ij0  = mapdouble(prhs[i++]);
    int w  = 0;
    int dorec       = mapdouble(prhs[i++]).cast<int>()(0);
    MatrixMd actualiters = mapdouble(prhs[i++]);
    MatrixMd w_out       = mapdouble(prhs[i++]);
    
    MatrixXd S (nvals,nvals);
    MatrixXd m0(nvals,1);
    
    MatrixXd b_i_save = b_i;
    
    int reps;
    double conv;
    for(reps=0; reps<maxiter; reps++){
        for(int c0=0; c0<ncliques*2; c0++){
            int c, mode;
            if(c0<ncliques){
                c    = c0;
                mode = 1;
            } else{
                c    = ncliques - 1 - (c0-ncliques);
                mode = 2;
            }
            int i = pairs(c,0);
            int j = pairs(c,1);

            if( mode==1 ){
                
                for(int yi=0; yi<nvals; yi++)
                    n(yi, i) = 1;
                for(int k=0; k<N1.cols(); k++){
                    int d = N1(i, k);
                    if(d==-2) continue;
                    for(int yi=0; yi<nvals; yi++)
                        n(yi, i) *= m1(yi, d);
                }
                for(int k=0; k<N2.cols(); k++){
                    int d = N2(i, k);
                    if(d==-2) continue;
                    for(int yi=0; yi<nvals; yi++)
                        n(yi, i) *= m2(yi, d);
                }
                //for(int yi=0; yi<nvals; yi++)
                //    n(yi, i) = pow(n(yi, i), rho);
                if(rho==1){
                    //nothing
                } else if(rho==.5){
                    n.col(i) = n.col(i).array().sqrt();
                } else
                    n.col(i) = n.col(i).array().pow(rho);
                
                // compute m(y_j)
                for(int yj=0; yj<nvals; yj++){
                    m0(yj) = 0;
                    for(int yi=0; yi<nvals; yi++){
                        int index = yi + yj*nvals;
                        S(yi,yj)  = psi_ij(index, c)*psi_i(yi, i)*n(yi, i)/m1(yi, c);
                        m0(yj)   += S(yi,yj);
                    }
                }
                if( dorec==1 )
                    for(int yj=0; yj<nvals; yj++){
                        mstor(w) = m2(yj,c);
                        w=w+1;
                    }
                double k = S.sum();
                m2.col(c) = m0/k;
                
            } else{
                for( int yj=0; yj<nvals; yj++)
                    n(yj, j) = 1;
                for(int k=0; k<N1.cols(); k++){
                    int d = N1(j, k);
                    if(d==-2) continue;
                    for( int yj=0; yj<nvals; yj++)
                        n(yj, j) *= m1(yj, d);
                }
                for(int k=0; k<N2.cols(); k++){
                    int d = N2(j, k);
                    if(d==-2) continue;
                    for( int yj=0; yj<nvals; yj++)
                        n(yj, j) *= m2(yj, d);
                }
                //for( int yj=0; yj<nvals; yj++)
                //    n(yj, j) = pow(n(yj, j), rho);
                if(rho==1){
                    //nothing
                } else if(rho==.5){
                    n.col(j) = n.col(j).array().sqrt();
                } else
                    n.col(j) = n.col(j).array().pow(rho);
                
                for(int yi=0; yi<nvals; yi++){
                    m0(yi) = 0;
                    for(int yj=0; yj<nvals; yj++){
                        int index = yi + yj*nvals;
                        S(yi,yj)  = psi_ij(index, c)*psi_i(yj, j)*n(yj, j)/m2(yj, c);
                        m0(yi)   += S(yi,yj);
                    }
                }
                if( dorec==1 )
                    for(int yi=0; yi<nvals; yi++){
                        mstor(w) = m1(yi,c);
                        w=w+1;
                    }
                double k = S.sum();
                m1.col(c) = m0/k;
            }
        }

        b_i = n.array()*psi_i.array();
        norm_cols(b_i);
        
        
        conv = (b_i-b_i_save).array().abs().maxCoeff();
        if(conv < convthresh){
            break;        
        }
        b_i_save = b_i;
    }
    actualiters(0) = reps+1;
    if(actualiters(0) > maxiter) actualiters(0)=maxiter;
    w_out(0) = w;
    
//     b_ij = psi_ij;
//     for(int c=0; c<ncliques; c++){
//         int i = pairs(c,0);
//         int j = pairs(c,1);
//         for(int yi=0; yi<nvals; yi++){
//             for(int yj=0; yj<nvals; yj++){
//                 int index = yi + yj*nvals;
//                 b_ij(index,c) = b_ij(index,c)*psi_i(yi,i)*psi_i(yj,j)*n(yi,i)*n(yj,j)/m1(yi,c)/m2(yj,c);
//             }
//         }
//     }
//     
//     Z_ij = b_ij.colwise().sum();
//     for(int c=0; c<ncliques; c++)
//         b_ij.col(c) = b_ij.col(c)/Z_ij(c);

    // recompute all inwards messages
    for(int i=0; i<b_i.cols(); i++){
        for(int yi=0; yi<nvals; yi++){
            n(yi, i) = 1;
            for(int k=0; k<N1.cols(); k++){
                int d = N1(i, k);
                if(d==-2) continue;
                n(yi, i) *= m1(yi, d);
            }
            for(int k=0; k<N2.cols(); k++){
                int d = N2(i, k);
                if(d==-2) continue;
                n(yi, i) *= m2(yi, d);
            }
        }
        
        if(rho==1){
            //nothing
        } else if(rho==.5){
            n.col(i) = n.col(i).array().sqrt();
        } else
            n.col(i) = n.col(i).array().pow(rho);
    }

    for(int c=0; c<ncliques; c++){
        int i = pairs(c, 0);
        int j = pairs(c, 1);
        for(int yi=0; yi<nvals; yi++){
            for(int yj=0; yj<nvals; yj++){
                int index = yi + yj*nvals;
                b_ij0(index, c) = b_ij0(index, c)*psi_i(yi, i)*psi_i(yj, j)*n(yi, i)*n(yj, j)/m1(yi, c)/m2(yj, c);
            }
        }
    }
    
}