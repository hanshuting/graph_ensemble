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

//meanfield_bprop_helper1(model, psi_ij, psi_i, maxiter, convthresh, ...
//    b_i, b_ij, mstor, dorec, actualiters, w);


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
    int maxiter     = mapdouble(prhs[i++])(0);
    double convthresh = mapdouble(prhs[i++])(0);
    MatrixMd b_ij   = mapdouble(prhs[i++]);
    MatrixMd b_i    = mapdouble(prhs[i++]);
    MatrixMd mstor  = mapdouble(prhs[i++]);
    int w  = 0;
    int dorec       = mapdouble(prhs[i++]).cast<int>()(0);
    MatrixMd actualiters = mapdouble(prhs[i++]);
    MatrixMd w_out       = mapdouble(prhs[i++]);
    
    MatrixXd b_i0(nvals,1);
    MatrixXd db_i0(nvals,1);
    
    MatrixXd b_i_save = b_i;
        
// % b_i0 = ones(model.nvals,1);
// % for reps=1:maxiter    
// %     b_i_old = b_i;
// %     
// %     for j=[1:model.nnodes model.nnodes-1:-1:2]
// %         if dorec, for yj=1:model.nvals, w=w+1; mstor(w) = b_i(yj,j); end; end
// %         
// %         for yj=1:model.nvals
// %             b_i0(yj) = psi_i(yj,j);
// %             for d=model.N1(j,:)
// %                 if d==-1, continue; end
// %                 i = model.pairs(d,2);
// %                 for yi=1:model.nvals
// %                     index = yj + (yi-1)*model.nvals;
// %                     %b_i(yj,j) = b_i(yj,j)*psi_ij(index,d)^b_i(yi,i);
// %                     b_i0(yj) = b_i0(yj)*psi_ij(index,d)^b_i(yi,i);
// %                 end
// %             end
// %             for d=model.N2(j,:)
// %                 if d==-1, continue; end
// %                 i = model.pairs(d,1);
// %                 for yi=1:model.nvals
// %                     index = yi + (yj-1)*model.nvals;
// %                     %b_i(yj,j) = b_i(yj,j)*psi_ij(index,d)^b_i(yi,i);
// %                     b_i0(yj) = b_i0(yj)*psi_ij(index,d)^b_i(yi,i);
// %                 end
// %             end
// %         end
// %         %b_i(:,j) = b_i0(:,j) / sum( b_i0(:,j) );
// %         b_i(:,j) = norm(b_i0);
// %     end
// %     
// %     conv = max(abs(b_i_old(:)-b_i(:)));
// %     if conv < convthresh
// %         % do backprop by same number of iterations
// %         maxiter = reps;
// %         break
// %     end
// % end
    
    int reps;
    double conv;
    for(reps=0; reps<maxiter; reps++){
        for(int j0=0; j0<nnodes*2-1; j0++){
            int j, mode;
            if(j0<nnodes)
                j    = j0;
            else
                j    = nnodes - 2 - (j0-nnodes);
            
            if( dorec==1 ){
                for(int yj=0; yj<nvals; yj++){
                    if(w >= mstor.rows()*mstor.cols()){
                        cout << "ERROR: " << w << " vs. " << mstor.rows()*mstor.cols() << endl;
                    }else{
                        mstor(w) = b_i(yj, j);
                        w=w+1;
                    }
                }
            }
            
            for(int yj=0; yj<nvals; yj++){
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
            double mysum = 0;
            for(int yj=0; yj<nvals; yj++)
                mysum += b_i0(yj);
            for(int yj=0; yj<nvals; yj++)
                b_i(yj,j) = b_i0(yj)/mysum;
        }
        
        conv = (b_i-b_i_save).array().abs().maxCoeff();
        if(conv < convthresh){
            break;        
        }
        b_i_save = b_i;
    }
    actualiters(0) = reps+1;
    if(actualiters(0) > maxiter) actualiters(0)=maxiter;
    w_out(0) = w;    
}