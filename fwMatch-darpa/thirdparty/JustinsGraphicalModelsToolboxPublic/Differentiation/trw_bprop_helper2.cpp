#include <stdio.h>
#include <map>
#include "mex.h"   //--This one is required
#include "math.h"
//#define mmin(X,Y) ((X) < (Y) ? (X) : (Y))
//#define mmax(X,Y) ((X) > (Y) ? (X) : (Y))
#include "../eigenstuff.cpp"
#include <iostream>

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
    // trw_bprop_helper2(model, psi_ij_rho, psi_i, rho, maxiter, n, m1, m2, b_i,...
    //                   dm1, dm2, dn, dpsi_i, dpsi_ij);
          
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
    double rho       = mapdouble(prhs[i++])(0);
    int maxiter      = mapdouble(prhs[i++])(0);
    MatrixMd n       = mapdouble(prhs[i++]);
    MatrixMd m1      = mapdouble(prhs[i++]);
    MatrixMd m2      = mapdouble(prhs[i++]);
    MatrixMd b_i     = mapdouble(prhs[i++]);
    MatrixMd mstor   = mapdouble(prhs[i++]);
    MatrixMd dm1     = mapdouble(prhs[i++]);
    MatrixMd dm2     = mapdouble(prhs[i++]);
    MatrixMd dn      = mapdouble(prhs[i++]);
    MatrixMd dpsi_i  = mapdouble(prhs[i++]);
    MatrixMd dpsi_ij = mapdouble(prhs[i++]);
    MatrixMd b_ij0   = mapdouble(prhs[i++]);
    MatrixMd db_ij0  = mapdouble(prhs[i++]);
    int dorec        = mapdouble(prhs[i++]).cast<int>()(0);
    
    int w  = mstor.rows()-1;
    
    MatrixXd S (nvals,nvals);
    MatrixXd m0(nvals,1);

    for(int c=0; c<ncliques; c++){
        int i = pairs(c, 0);
        int j = pairs(c, 1);
        for(int yi=0; yi<nvals; yi++){
            for(int yj=0; yj<nvals; yj++){
                int index = yi + yj*nvals;
                //b_ij0(index,c) = b_ij0(index,c)*psi_i(yi,i)*psi_i(yj,j)*n(yi,i)*n(yj,j)/m1(yi,c)/m2(yj,c);
                dpsi_i(yi, i) = dpsi_i(yi, i) + db_ij0(index, c)*b_ij0(index, c)/psi_i(yi, i);
                dpsi_i(yj, j) = dpsi_i(yj, j) + db_ij0(index, c)*b_ij0(index, c)/psi_i(yj, j);
                dn(yi, i)     = dn(yi, i)     + db_ij0(index, c)*b_ij0(index, c)/n(yi, i);
                dn(yj, j)     = dn(yj, j)     + db_ij0(index, c)*b_ij0(index, c)/n(yj, j);
                dm1(yi, c)    = dm1(yi, c)    - db_ij0(index, c)*b_ij0(index, c)/m1(yi, c);
                dm2(yj, c)    = dm2(yj, c)    - db_ij0(index, c)*b_ij0(index, c)/m2(yj, c);
            }
        }
    }

    for(int i=0; i<b_i.cols(); i++){
        for(int yi=0; yi<nvals; yi++){
            for(int k=0; k<N1.cols(); k++){
                int d = N1(i, k);
                if(d==-2) continue;
                //n(yi, i) *= m1(yi, d);
                dm1(yi,d) += rho * dn(yi,i)*n(yi,i)/m1(yi,d);
            }
            for(int k=0; k<N2.cols(); k++){
                int d = N2(i, k);
                if(d==-2) continue;
                //n(yi, i) *= m2(yi, d);
                dm2(yi,d) += rho * dn(yi,i)*n(yi,i)/m2(yi,d);
            }
        }
    }
    
    int reps;
    double conv;
    for(reps=0; reps<maxiter; reps++){
        
        for(int c0=2*ncliques-1; c0>=0; c0--){
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
                // n.col(i) = n.col(i).array().pow(rho);
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
            
                double k = S.sum();
                
                MatrixXd dm0 = dm2.col(c)/k; 
                dm0.array() -= (dm2.col(c).array()*m0.array()).sum()/(k*k);
                
                MatrixXd dn  = 0*n.col(i);
                
                //#pragma omp parallel for num_threads(2)
                for(int yj=0; yj<nvals; yj++){
                    for(int yi=0; yi<nvals; yi++){
                        int index = yi + yj*nvals;
                        dpsi_ij(index,c) += dm0(yj)*S(yi,yj)/psi_ij(index,c);
                        dpsi_i(yi,i)     += dm0(yj)*S(yi,yj)/psi_i(yi,i);
                        dn(yi)           += dm0(yj)*S(yi,yj)/n(yi,i);
                        dm1(yi,c)        -= dm0(yj)*S(yi,yj)/m1(yi,c);
                    }
                }
                
               for(int yi=0; yi<nvals; yi++){
                    for(int k=0; k<N1.cols(); k++){
                        int d = N1(i, k);
                        if(d==-2) continue;
                        dm1(yi,d) += rho*dn(yi)*n(yi,i)/m1(yi,d);
                    }
                    for(int k=0; k<N2.cols(); k++){
                        int d = N2(i, k);
                        if(d==-2) continue;
                        dm2(yi,d) += rho*dn(yi)*n(yi,i)/m2(yi,d);
                    }
               }
               dm2.col(c) *= 0.0;
               if( dorec )
                   for(int yj=nvals-1; yj>=0; yj--){ m2(yj,c)=mstor(w); w--;}
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
                
//                n.col(j) = n.col(j).array().pow(rho);
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
                double k = S.sum();

                MatrixXd dm0 = dm1.col(c)/k; 
                dm0.array() -= (dm1.col(c).array()*m0.array()).sum()/(k*k);
                MatrixXd dn  = 0*n.col(i);

                for(int yj=0; yj<nvals; yj++){
                    for(int yi=0; yi<nvals; yi++){
                        int index = yi + yj*nvals;
                        dpsi_ij(index,c) += dm0(yi)*S(yi,yj)/psi_ij(index,c);
                        dpsi_i(yj,j)     += dm0(yi)*S(yi,yj)/psi_i(yj,j);
                        dn(yj)           += dm0(yi)*S(yi,yj)/n(yj,j);
                        dm2(yj,c)        -= dm0(yi)*S(yi,yj)/m2(yj,c);
                    }
                }
                
                for(int yj=0; yj<nvals; yj++){
                    for(int k=0; k<N1.cols(); k++){
                        int d = N1(j, k);
                        if(d==-2) continue;
                        dm1(yj, d) += rho*dn(yj)*n(yj, j)/m1(yj, d);
                    }
                    for(int k=0; k<N2.cols(); k++){
                        int d = N2(j, k);
                        if(d==-2) continue;
                        dm2(yj, d) += rho*dn(yj)*n(yj, j)/m2(yj, d);
                    }
                }
                
               dm1.col(c) *= 0.0;
               if( dorec )
                   for(int yi=nvals-1; yi>=0; yi--){ m1(yi,c)=mstor(w); w--;}
            }
        }
    }
}

