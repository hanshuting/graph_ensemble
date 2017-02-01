/* Copyright (c) 2006, National ICT Australia
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * Authors      : Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 * Created      : 20/11/2007
 * Last Updated :
 */

#ifndef _L2N2_BMRMINNERSOLVER_CPP_
#define _L2N2_BMRMINNERSOLVER_CPP_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "l2n2_bmrmdualinnersolver.hpp"
#include "configuration.hpp"

// which type of gradient to store in gradientSet?
#ifdef REAL_NUMBER_IS_FLOAT
        #define DENSE_TO_SPARSE_THRESHOLD  0.4    // actual threshold = 0.5
#else
        #define DENSE_TO_SPARSE_THRESHOLD  0.6    // actual threshold = 0.75
#endif


// remove ALL of the gradients with idle age over gradIdleAge
//#define REMOVE_ALL_IDLE_GRADIENTS
// otherwise, remove the oldest gradient with idle age over gradIdleAge
//chteo: in general, BMRM needs more number of iterations to converge if we remove all IDLE gradients at onece

CL2N2_BMRMDualInnerSolver::CL2N2_BMRMDualInnerSolver(double lambda)
   : CBMRMInnerSolver(lambda),
     f(0),
     Q(0),
     a(0),
     b(0),
     l(0),
     u(0),
     tol(1e-6)
{
   // inherited from CInnerSolver
   iter = 0;
   dim = 0;
   numOfConstraint = 1;
   
   // set (default) private member values
   gradIdleAge = 10;
   gradType = SML::DENSE;
   maxGradSetSize = 2048;  // max num of gradients

   // synchronize members with configuration file
   Configuration &config = Configuration::GetInstance();

   if(config.IsSet("L2N2_BMRMDualInnerSolver.maxGradSetSize"))
      maxGradSetSize = config.GetInt("L2N2_BMRMDualInnerSolver.maxGradSetSize");

   if(config.IsSet("L2N2_BMRMDualInnerSolver.gradIdleAge")) 
   {
     gradIdleAge = config.GetInt("L2N2_BMRMDualInnerSolver.gradIdleAge");
     gradIdleAge = std::max(gradIdleAge, 2);
   }

   // pre-allocate memory for gradient, offset and time-stamps
   gradientSet.resize(maxGradSetSize,0);                  
   offsetSet.resize(maxGradSetSize,0);                  
   timeStamp.resize(maxGradSetSize,-1);   
   
   // allocate maximum memory needed
   x = (double*)calloc(maxGradSetSize, sizeof(double));
   Q = (double*)calloc(maxGradSetSize*maxGradSetSize, sizeof(double));
   f = (double*)calloc(maxGradSetSize, sizeof(double));
   l = (double*)calloc(maxGradSetSize, sizeof(double));
   u = (double*)calloc(maxGradSetSize, sizeof(double));
   a = (double*)calloc(maxGradSetSize, sizeof(double));
   b = new double(1); 
   
   // initializations
   Reset();

}

CL2N2_BMRMDualInnerSolver::~CL2N2_BMRMDualInnerSolver()
{
        // free for parent
        if(x) free(x);
        if(Q) free(Q);
        if(f) free(f);
        if(l) free(l);
        if(u) free(u);
        if(a) free(a);
        if(b) delete b;
   
        for(int i=0; i < dim; i++)
                if(gradientSet[i])
                        delete gradientSet[i];   
}


void CL2N2_BMRMDualInnerSolver::Reset()
{
        // QP mem and gradient related mem
        for(int i=0; i<maxGradSetSize; i++) 
        {
                x[i] = 0;
                f[i] = 0;                
                l[i] = 0;
                u[i] = SML::INFTY;
                a[i] = 1;      
                
                if(gradientSet[i]) 
                {
                        delete gradientSet[i];
                        gradientSet[i] = 0;
                }
                offsetSet[i] = 0;
                timeStamp[i] = -1;
        }                
        memset(Q, 0, sizeof(double)*maxGradSetSize*maxGradSetSize); 
        
        // variables
        iter = 0;
        dim = 0;
}


#ifdef REMOVE_ALL_IDLE_GRADIENTS
/**  Update matrix and vectors needed in optimisation round.
 *
 *   \param a [read] Gradient
 *   \param b [read] Offset 
 */
void CL2N2_BMRMDualInnerSolver::Update(TheMatrix& a, Scalar b)
{
        iter++;

        // remove idle gradients
        int first_idle = 0;
        int last_idle = dim-1;
        int removeCnt = 0;
        
        // those to-keep gradients float to the top
        while(first_idle < last_idle)
        {                
                // look for the bottom most to-keep gradient and remove bottom most to-remove gradients
                while((last_idle > 0) && (iter-timeStamp[last_idle] >= gradIdleAge)) 
                {
                        if(gradientSet[last_idle]) 
                        {
                                delete gradientSet[last_idle];
                                gradientSet[last_idle] = 0;                                
                                offsetSet[last_idle] = 0;
                                x[last_idle] = 0;
                                f[last_idle] = 0; 
                                timeStamp[last_idle] = -1;
                        }                        
                        removeCnt++; 
                        last_idle--;                              
                }
                                
                // look for top most to-remove gradient
                while((first_idle < dim) && (iter-timeStamp[first_idle] < gradIdleAge)) 
                        first_idle++;

                // replace the top most to-remove gradient wiht the bottom most to-keep gradient
                if(first_idle < last_idle)
                {       
                        // 1. remove/replace the elements in gradientSet, offsetSet, x, f, timeStamp                        
                        if(gradientSet[first_idle]) 
                                delete gradientSet[first_idle];
                        gradientSet[first_idle] = gradientSet[last_idle];
                        gradientSet[last_idle] = 0;
                        
                        offsetSet[first_idle] = offsetSet[last_idle];
                        offsetSet[last_idle] = 0;
                        
                        x[first_idle] = x[last_idle];
                        x[last_idle] = 0;
                        
                        f[first_idle] = f[last_idle];                                                     
                        f[last_idle] = 0;
                        
                        timeStamp[first_idle] = timeStamp[last_idle];
                        timeStamp[last_idle] = -1;
                        
                        // 2. memmove row                
                        memmove(Q+first_idle*dim, Q+last_idle*dim, sizeof(double)*dim);
                        
                        // 3. memmove column
                        for(int i=0; i<last_idle; i++)
                                Q[i*dim + first_idle] = Q[i*dim + last_idle];                                
                        
                        // 4. one more to-remove gradient sink to the bottom
                        removeCnt++;
                }
        }
        
        // make the hessian matrix compact
        if(removeCnt > 0)
        {
                int prevdim = dim;
                dim -= removeCnt;
                assert(dim >= 0);
                
                for(int i=1; i<dim; i++)
                        memmove(Q+i*dim, Q+i*prevdim, sizeof(double)*dim);
        }
        
        if(verbosity > 0)
                std::cout << "removed gradient: " << removeCnt << std::endl;
        
        // new gradient just comes in... increase dim of vector of lagrangian multipliers.
        dim++;
        int idx = dim-1; 
        
        // store new column at #idx#      
        if(gradientSet[idx]) delete gradientSet[idx];  // clean up slot
        if(a.Density() < DENSE_TO_SPARSE_THRESHOLD)
                gradientSet[idx] = new TheMatrix(a, SML::SPARSE);    // store the gradient                
        else
                gradientSet[idx] = new TheMatrix(a, SML::DENSE);    // store the gradient
        offsetSet[idx] = b;      // store the offset  
   
        // use new slot
        timeStamp[idx] = iter;  
      
        // step 1: insert new row into Q
        int lastRow = (dim-1)*dim;
        for(int col=0; col < dim; col++) 
        {
                Scalar value = 0;
                gradientSet[idx]->Dot(*gradientSet[col],value);
                Q[lastRow+col] = 1.0/bmrmLambda*value;
        }

        // step 2: adjust memory area to fit a new column
        int prevDim = dim-1;
        int size = sizeof(double)*prevDim;
        for(int row=dim-2; row >= 1; row--)
        {
                double *src  = Q + (row*prevDim);
                double *dest = Q + (row*prevDim)+ row;
                memmove(dest, src, size);
        }

        // step 3: insert new column (at the right-most-side of the matrix)
        for(int i=0; i < dim-1; i++)
                Q[i*dim + dim-1] = Q[lastRow + i];

        // update vectors
        x[dim-1] = 0;
        f[dim-1] = -offsetSet[dim-1];
        
#ifdef DEBUG
        // check Q matrix !    
        MatrixCorrectnessCheck();
#endif  
}

#else
/**  Update matrix and vectors needed in optimisation round.
 *
 *   \param a [read] Gradient
 *   \param b [read] Offset 
 */
void CL2N2_BMRMDualInnerSolver::Update(TheMatrix& a, Scalar b)
{
        int idx = -1;                      // default: look for i such that x[i] is smallest
        int first_new_slot = -1;
        int oldest_slot = 0;
        int smallest_slot = 0;
   
        iter++;                            // another column generated
   
        // deciding which slot to store the latest gradient
        // steps:
        // 1. if there is any gradient qualified for retirement, then replace it.
        // 2. else if this is the first iteration, use the first slot
        // 3. otherwise, occupy a new slot
        
        for(int i=0; i < maxGradSetSize; i++) 
        {
                if(timeStamp[i] == -1 && first_new_slot == -1) first_new_slot = i;
                if(timeStamp[i] < timeStamp[oldest_slot] && timeStamp[i] >= 0) oldest_slot = i; 
                if(x[i] < x[smallest_slot]) smallest_slot = i;
        }
        
        if(timeStamp[oldest_slot] < iter - gradIdleAge) 
                idx = oldest_slot;
        else if(first_new_slot != -1) 
        {
                idx = first_new_slot;
                dim++;
        }
        else 
                idx = smallest_slot;
   
#ifdef DEBUG   
        std::cout << "new column idx: " << idx << "   gradientSet size:" << gradientSet.size() << std::endl;
#endif      
   
        // store new column at #idx#      
        if(gradientSet[idx]) delete gradientSet[idx];  // clean up slot
        if(a.Density() < DENSE_TO_SPARSE_THRESHOLD)
                gradientSet[idx] = new TheMatrix(a, SML::SPARSE);    // store the gradient                
        else
                gradientSet[idx] = new TheMatrix(a, SML::DENSE);    // store the gradient
        offsetSet[idx] = b;      // store the offset  
   
        // use new slot
        if(timeStamp[idx] == -1)
        {
                // This gradient is new!
                timeStamp[idx] = iter;  
      
                // step 1: insert new row into Q
                int lastRow = (dim-1)*dim;
                for(int col=0; col < dim; col++) 
                {
                        Scalar value = 0;
                        gradientSet[idx]->Dot(*gradientSet[col],value);
                        Q[lastRow+col] = 1.0/bmrmLambda*value;
                }
      
                // step 2: adjust memory area to fit a new column
                int prevDim = dim-1;
                int size = sizeof(double)*prevDim;
                for(int row=dim-2; row >= 1; row--)
                {
                        double *src  = Q + (row*prevDim);
                        double *dest = Q + (row*prevDim)+ row;
                        memmove(dest, src, size);
                }
      
                // step 3: insert new column (at the right-most-side of the matrix)
                for(int i=0; i < dim-1; i++)
                        Q[i*dim + dim-1] = Q[lastRow + i];
      
                // update vectors
                x[dim-1] = 0;
                f[dim-1] = -offsetSet[dim-1];
        }
        // reuse slot
        else 
        { 
                // This gradient is new!
                timeStamp[idx] = iter;  
      
                // insert new row into Q[idx,:]
                int theRow = idx*dim;
                for(int col=0; col < dim; col++) 
                {
                        Scalar value = 0;
                        gradientSet[idx]->Dot(*gradientSet[col],value);
                        Q[theRow+col] = 1.0/bmrmLambda*value;
                }
      
                // insert new column into Q[:,idx]
                for(int i=0; i < dim; i++)
                        Q[i*dim + idx] = Q[theRow + i];
      
                // update vectors      
                x[idx] = 0;
                f[idx] = -offsetSet[idx];
        }
   
#ifdef DEBUG
        // check Q matrix !    
        MatrixCorrectnessCheck();
#endif  
}
#endif


/** Return solution.
 *
 *  \param w      [write] solution vector (preallocated)
 *  \param xi     [write] lower bound of empirical risk, R_emp[w]
 *  \param regval [write] value of regularization term
 *  \param objval [write] objective value
 */
void CL2N2_BMRMDualInnerSolver::GetSolution(TheMatrix& w, Scalar &xi, Scalar &regval, Scalar &objval)
{
        assert(x != 0);
        double factor = 1.0/bmrmLambda;  // bmrmLambda is regularization constant
   
   
        // compute objective value (explicitly)
        // n := dimensionality of x (the solution of QP)
        // Q := hessian matrix
        // f := linear part of obj

        objval = 0;
        double tmp = 0;
        double fx = 0;
        for(int i=0; i < dim; i++) 
        {
                tmp = 0;
                for(int j=0; j < dim; j++)
                        tmp += Q[j + i*dim]*x[j];
                objval += x[i]*tmp;
                fx += f[i]*x[i];
        }
   
        // since the dual of minimization problem is maximization
        // and we solved the  -maximization version,
        // the objective value should be -objval
        objval = -0.5*objval - fx;
   
   
        // Compute new w
        w.Zero();
        for(int i=0; i < dim; i++)
                if(x[i] > L2N2_BMRMDualInnerSolver::ZERO_EPS)
                        w.ScaleAdd(-x[i], *gradientSet[i]);
        w.Scale(factor);
   
        // compute xi
        w.Norm2(regval);
        regval = regval * regval * bmrmLambda * 0.5; 
        xi = objval - regval;
   
   
        // update time-stamp of constraints
        for(int i=0; i < dim; i++)
                if(x[i] > L2N2_BMRMDualInnerSolver::ZERO_EPS) 
                        timeStamp[i] = iter;         
   
        // count nonzero terms
        if(verbosity > 0) 
        {
                double smallest = 1000000;
                int howmany=0;
                //std::cout << "nz: ";
                for(int i=0; i < dim; i++) 
                {
                        if(SML::abs(x[i])> L2N2_BMRMDualInnerSolver::ZERO_EPS) 
                        {
                                //std::cout << i << "..";
                                howmany++;
                                smallest = std::min(x[i],smallest);
                        }
                }
                std::cout << "\ngradient nnz: " << howmany << std::endl;
                std::cout << "qp size: " << dim << std::endl;
                std::cout << "smallest |x_i| > ZERO_EPS(" << L2N2_BMRMDualInnerSolver::ZERO_EPS << "): " << smallest << std::endl;
        }
   
}


void CL2N2_BMRMDualInnerSolver::Solve(TheMatrix& w, TheMatrix& a, Scalar loss, Scalar &xi, Scalar &regval, Scalar &objval)
{
        Scalar w_dot_a = 0.0;
        w.Dot(a, w_dot_a);
        Update(a, loss - w_dot_a);

        SolveQP();

        GetSolution(w, xi, regval, objval);
}


#ifdef DEBUG
/** Check if the Q matrix update is correct
 */
void CL2N2_BMRMDualInnerSolver::MatrixCorrectnessCheck() 
{
        std::cout << "in matrix correctness check! " << std::endl;
        Scalar dot = 0;
        double *correctmat = new double[dim*dim];
        memset(correctmat, 0, sizeof(double)*dim*dim);
   
        for(int i=0; i < dim; i++)
        {
                for(int j=i; j < dim; j++)
                {
                        //std::cout << "i:" << i << "   j:" << j << std::endl;
                        dot = 0;
                        gradientSet[i]->Dot(*gradientSet[j], dot);
                        correctmat[i*dim + j] = 1.0/bmrmLambda*dot;
                        correctmat[j*dim + i] = correctmat[i*dim+j];
                }
        }
   
        for(int i=0; i < dim*dim; i++)
        {
                if(fabs(correctmat[i] - Q[i]) > 1e-15) 
                {
                        std::cout << "residual : " << fabs(correctmat[i] - Q[i]) << std::endl;
                        std::cout << "matrix update correctness check : i = " << i << std::endl;
                        std::cout << "ERROR: rank-1 updated matrix is incorrect ! " << std::endl;
                        exit(EXIT_FAILURE);
                }
        }
   
        if(correctmat) delete[] correctmat;
}
#endif

#endif
