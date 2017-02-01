#include <stdio.h>
#include <map>
#include "mex.h"   //--This one is required
#include "math.h"
#include "../eigenstuff.cpp"
#include <iostream>

#define NTHREAD 6

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    // trw_bprop_helper2(model, psi_ij_rho, psi_i, rho, maxiter, n, m1, m2, b_i,...
    //                   dm1, dm2, dn, dpsi_i, dpsi_ij);
          
    int i=0;
    int nnodes     = mapdouble(mxGetField(prhs[i],0,"nnodes"))(0);
    int ncliques   = mapdouble(mxGetField(prhs[i],0,"ncliques"))(0);
    int nvals      = mapdouble(mxGetField(prhs[i],0,"nvals"))(0);
    MatrixXi pairs = mapdouble(mxGetField(prhs[i],0,"pairs")).cast<int>(); pairs.array() -= 1;
    MatrixXi N1   = mapdouble(mxGetField(prhs[i],0,"N1")).cast<int>(); N1.array() -= 1;
    MatrixXi N2   = mapdouble(mxGetField(prhs[i],0,"N2")).cast<int>(); N2.array() -= 1;    
    MatrixXi tree2clique  = mapdouble(mxGetField(prhs[i],0,"tree2clique")).cast<int>();  tree2clique.array() -= 1;
    MatrixXi treeschedule = mapdouble(mxGetField(prhs[i],0,"treeschedule")).cast<int>(); treeschedule.array() -= 1;
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
    
    MatrixMi w      = mapint32(prhs[i++]);

    #pragma omp parallel for num_threads(NTHREAD)
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

    #pragma omp parallel for num_threads(NTHREAD)
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
        // must re-order blocks
        //for(int block=0; block<treeschedule.cols(); block++){
        for(int block=treeschedule.cols()-1; block>=0; block--){
            // need not re-order trees, but why not...
            // helps if someone specifies not parallel trees
            // for(int treenum=treeschedule.rows()-1; treenum>=0; treenum--){
            #pragma omp parallel for schedule(dynamic) num_threads(NTHREAD)
            for(int treenum=0; treenum<treeschedule.rows(); treenum++){
                MatrixXd S (nvals,nvals);
                MatrixXd m0(nvals,1);
                
                int tree = treeschedule(treenum,block);
                if(tree==-1)
                    continue;
                int ncliques = tree_ncliques(tree);
                
                for(int c0=2*ncliques-1; c0>=0; c0--){
                    int c, mode;
                    if(c0<ncliques){
                        c    = tree2clique(c0,tree);
                        mode = 1;
                    } else{
                        c    = tree2clique(ncliques - 1 - (c0-ncliques),tree);
                        mode = 2;
                    }
    
                    //cout << "BW c " << c << " mode " << mode << endl;
                    
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
                            for(int yj=nvals-1; yj>=0; yj--){
                                w(tree)=w(tree)-1;
                                m2(yj,c)=mstor(w(tree),tree);
                            }
                    } else if( mode==2 ){
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
                            for(int yi=nvals-1; yi>=0; yi--){
                                w(tree)=w(tree)-1;
                                m1(yi,c)=mstor(w(tree),tree);
                            }
                    }
                }
            }
        }
    }
}

