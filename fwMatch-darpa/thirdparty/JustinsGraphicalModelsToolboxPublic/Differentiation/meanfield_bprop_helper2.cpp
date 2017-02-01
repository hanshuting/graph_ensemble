#include <stdio.h>
#include <map>
#include "mex.h"   //--This one is required
#include "math.h"
//#define mmin(X,Y) ((X) < (Y) ? (X) : (Y))
//#define mmax(X,Y) ((X) > (Y) ? (X) : (Y))
#include "../eigenstuff.cpp"
#include <iostream>

//meanfield_bprop_helper2(model, psi_ij, psi_i, maxiter, convthresh, ...
//    b_ij, b_i, mstor, dorec, actualiters, w, dpsi_ij, dpsi_i, db_ij, db_i);


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
    int maxiter     = mapdouble(prhs[i++])(0);
    double convthresh = mapdouble(prhs[i++])(0);
    MatrixMd b_ij   = mapdouble(prhs[i++]);
    MatrixMd b_i    = mapdouble(prhs[i++]);
    MatrixMd mstor  = mapdouble(prhs[i++]);
    int dorec       = mapdouble(prhs[i++]).cast<int>()(0);
    
    MatrixMd w_out       = mapdouble(prhs[i++]);
    int w = w_out(0);
    
    MatrixMd dpsi_ij     = mapdouble(prhs[i++]);
    MatrixMd dpsi_i      = mapdouble(prhs[i++]);
    MatrixMd db_ij       = mapdouble(prhs[i++]);
    MatrixMd db_i        = mapdouble(prhs[i++]);
    
    MatrixXd b_i0(nvals,1);
    MatrixXd db_i0(nvals,1);
    
    MatrixXd b_i_save = b_i;
        
    for(int c=0; c<ncliques; c++){
        for(int yi=0; yi<nvals; yi++){
            for(int yj=0; yj<nvals; yj++){
                int i = pairs(c,0);
                int j = pairs(c,1);
                int index = yi + yj*nvals;
                db_i(yi,i) += db_ij(index,c)*b_ij(index,c)/b_i(yi,i);
                db_i(yj,j) += db_ij(index,c)*b_ij(index,c)/b_i(yj,j);
            }
        }
    }
        
    
// % for c=1:model.ncliques
// %     for yi=1:model.nvals
// %         for yj=1:model.nvals
// %             i = model.pairs(c,1);
// %             j = model.pairs(c,2);
// %             index = yi + (yj-1)*model.nvals;
// %             db_i(yi,i) = db_i(yi,i) + db_ij(index,c)*b_ij(index,c)/b_i(yi,i);
// %             db_i(yj,j) = db_i(yj,j) + db_ij(index,c)*b_ij(index,c)/b_i(yj,j);
// %         end
// %     end
// % end

    
    //cout << "w " << w << endl;
    //return;
    
    int reps;
    double conv;
    for(reps=0; reps<maxiter; reps++){
        for(int j0=nnodes*2-2; j0>=0; j0--){
            
            int j, mode;
            if(j0<nnodes)
                j    = j0;
            else
                j    = nnodes - 2 - (j0-nnodes);
            
            // recover non-normalized messages
            for(int yj=0; yj<nvals; yj++){
                //cout << b_i0 << "--\n" << endl;
                b_i0(yj) = psi_i(yj,j);
                
                for(int k=0; k<N1.cols(); k++){
                    int d = N1(j, k);
                    if(d==-2) continue;
                    int i = pairs(d,1);
                    for(int yi=0; yi<nvals; yi++){
                        int index = yj + yi*nvals;
                        b_i0(yj) *= pow(psi_ij(index,d),b_i(yi,i));
                    }
                }
                for(int k=0; k<N2.cols(); k++){
                    int d = N2(j, k);
                    if(d==-2) continue;
                    int i = pairs(d,0);
                    for(int yi=0; yi<nvals; yi++){
                        int index = yi + yj*nvals;
                        b_i0(yj) *= pow(psi_ij(index,d),b_i(yi,i));
                    }
                }
            }
            //cout << "b_i0\n" << b_i0 << endl;
            //cout << "db_i\n" << db_i.col(j) << endl;
            
            double s = 0;
            for(int yj=0; yj<nvals; yj++)
                s += b_i0(yj);
            
            double inner = 0;
            for(int yj=0; yj<nvals; yj++)
                inner += db_i(yj,j)*b_i0(yj);
            for(int yj=0; yj<nvals; yj++)
                db_i0(yj) = db_i(yj,j)/s - inner/s/s;
    
            //cout << "db_i0\n\n" << db_i0 << endl;
            
            for(int yj=0; yj<nvals; yj++){
                for(int k=0; k<N1.cols(); k++){
                    int d = N1(j, k);
                    if(d==-2) continue;
                    int i = pairs(d,1);
                    for(int yi=0; yi<nvals; yi++){
                        int index = yj + yi*nvals;
                        //b_i0(yj) *= pow(psi_ij(index,d),b_i(yi,i));
                        db_i(yi,i) += db_i0(yj)*b_i0(yj)*log(psi_ij(index,d));
                        dpsi_ij(index,d) += db_i0(yj)*b_i0(yj)/psi_ij(index,d)*b_i(yi,i);
                    }
                }
                for(int k=0; k<N2.cols(); k++){
                    int d = N2(j, k);
                    if(d==-2) continue;
                    int i = pairs(d,0);
                    for(int yi=0; yi<nvals; yi++){
                        int index = yi + yj*nvals;
                        //b_i0(yj) *= pow(psi_ij(index,d),b_i(yi,i));
                        db_i(yi,i) += db_i0(yj)*b_i0(yj)*log(psi_ij(index,d));
                        dpsi_ij(index,d) += db_i0(yj)*b_i0(yj)/psi_ij(index,d)*b_i(yi,i);
                    }
                }
            }
        
            for(int yj=0; yj<nvals; yj++)
                dpsi_i(yj,j) += db_i0(yj)*b_i0(yj)/psi_i(yj,j);
        
            if(dorec==1){
                for(int yj=nvals-1; yj>=0; yj--){
                    w--;
                    if(w<0){
                        cout << "W ERROR!" << endl;
                        return;
                    }
                    b_i(yj,j) = mstor(w);
                } 
            }
            for(int yj=0; yj<nvals; yj++)
                db_i(yj,j) = 0;
        }
    }
}