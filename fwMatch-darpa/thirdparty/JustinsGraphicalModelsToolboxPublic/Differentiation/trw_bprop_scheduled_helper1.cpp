#include <stdio.h>
#include <map>
#include "mex.h"   //--This one is required
#include "math.h"
#include "../eigenstuff.cpp"

#define NTHREAD 6

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    // trw_bprop_helper1(model, psi_ij_rho, psi_i, rho, maxiter, n, m1, m2, b_i);
          
    int i=0;
    int nnodes     = mapdouble(mxGetField(prhs[i],0,"nnodes"))(0);
    int ncliques   = mapdouble(mxGetField(prhs[i],0,"ncliques"))(0);
    int nvals      = mapdouble(mxGetField(prhs[i],0,"nvals"))(0);
    MatrixXi pairs = mapdouble(mxGetField(prhs[i],0,"pairs")).cast<int>(); pairs.array() -= 1;
    MatrixXi N1    = mapdouble(mxGetField(prhs[i],0,"N1")).cast<int>(); N1.array() -= 1;
    MatrixXi N2    = mapdouble(mxGetField(prhs[i],0,"N2")).cast<int>(); N2.array() -= 1;
    MatrixXi tree2clique  = mapdouble(mxGetField(prhs[i],0,"tree2clique")).cast<int>(); tree2clique.array() -= 1;
    MatrixXi treeschedule = mapdouble(mxGetField(prhs[i],0,"treeschedule")).cast<int>(); treeschedule.array() -= 1;
    
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
    
    int dorec       = mapdouble(prhs[i++]).cast<int>()(0);
    MatrixMd actualiters = mapdouble(prhs[i++]);
    MatrixMi w      = mapint32(prhs[i++]);
    
    MatrixXd b_i_save = b_i;
    
    //tree_ncliques = sum(double(model.tree2clique>0),1);
    //ntree = length(tree_ncliques);
    int ntree = tree2clique.cols();
    MatrixXi tree_ncliques = MatrixXi::Zero(ntree,1);
    for(int tree=0; tree<tree2clique.cols(); tree++)
        for(i=0; i<tree2clique.rows(); i++)
            if(tree2clique(i,tree) != -1)
                tree_ncliques(tree)++;

    int reps;
    double conv;
    for(reps=0; reps<maxiter; reps++){
        for(int block=0; block<treeschedule.cols(); block++){
            #pragma omp parallel for schedule(dynamic) num_threads(NTHREAD)
            for(int treenum=0; treenum<treeschedule.rows(); treenum++){
                MatrixXd S (nvals,nvals);
                MatrixXd m0(nvals,1);
                
                int tree = treeschedule(treenum,block);
                if(tree==-1)
                    continue;
                int ncliques = tree_ncliques(tree);
                for(int c0=0; c0<ncliques*2; c0++){
                    int c, mode;
                    if(c0<ncliques){
                        c    = tree2clique(c0,tree);
                        mode = 1;
                    } else{
                        c    = tree2clique(ncliques - 1 - (c0-ncliques),tree);
                        mode = 2;
                    }
                    
                    //cout << "FW c " << c << " mode " << mode << endl;
                    
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
                                mstor(w(tree),tree) = m2(yj,c);
                                w(tree)=w(tree)+1;
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
                        
                        //#pragma omp parallel for num_threads(10)
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
                                mstor(w(tree),tree) = m1(yi,c);
                                w(tree)=w(tree)+1;
                            }
                        double k = S.sum();
                        m1.col(c) = m0/k;
                    }
                }
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
    #pragma omp parallel for num_threads(NTHREAD)
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

    #pragma omp parallel for num_threads(NTHREAD)
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