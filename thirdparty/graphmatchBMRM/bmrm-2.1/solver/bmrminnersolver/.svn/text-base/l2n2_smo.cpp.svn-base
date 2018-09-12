/**
 *
 *
 *   Authors : Ivor Tsang (ivo@cse.ust.hk)
 *             Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *   Code adapted from CVM
 */

#ifndef _L2N2_SMO_CPP_
#define _L2N2_SMO_CPP_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "l2n2_smo.hpp"
#include "configuration.hpp"

CL2N2_SMO::CL2N2_SMO(double lambda)
   : CL2N2_BMRMDualInnerSolver(lambda),
     maxSMOIter(10000),
     smo_eps(1e-6),
     smo_eps_sq(smo_eps*smo_eps),
     G(0),
     QD(0),
     x_status(0)
{
   // synchronize members with configuration file
   Configuration &config = Configuration::GetInstance();
   
   if(config.IsSet("L2N2_SMO.maxIter"))
      maxSMOIter = config.GetInt("L2N2_SMO.maxIter");
   
   if(config.IsSet("L2N2_SMO.tolerance")) {
      smo_eps = config.GetDouble("L2N2_SMO.tolerance");
      smo_eps_sq = smo_eps*smo_eps;
   }
   
   
   // pre-allocate memory for gradient, offset and time-stamps
   gradientSet.resize(maxGradSetSize);                  
   offsetSet.resize(maxGradSetSize);                  
   timeStamp.resize(maxGradSetSize, -1);
   
   G = (double*)calloc(maxGradSetSize, sizeof(double));
   QD = (double*)calloc(maxGradSetSize, sizeof(double));
   x_status = (int*)calloc(maxGradSetSize, sizeof(int));
   
   // initialize x_status
   for(int i=0; i<maxGradSetSize; i++)
      x_status[i] = SMO::LOWER_BOUND;


   // make the initial solution feasible
   x[0] = 1.0;
   x_status[0] = SMO::FREE;
   
}

CL2N2_SMO::~CL2N2_SMO()
{  
   if(G) free(G);
   if(QD) free(QD);
   if(x_status) free(x_status);  
}



/**  Solve:
 *   
 *   min  0.5*x'Qx -b'x
 *   s.t. x[i] >= 0 \forall i
 *        1'x = 1
 *
 */

void CL2N2_SMO::SolveQP()
{
   int i=-1, j=-1;
   int smo_iter = 0;

 
   // optimization step
   for(smo_iter=0; smo_iter < maxSMOIter; smo_iter++)
   {
      //std::cout << "before wss" << std::endl;
      
      int lbcnt = 0;
      for(int kk=0; kk < dim; kk++)
         if(x_status[kk] == SMO::LOWER_BOUND) lbcnt++;

      if(smo_select_working_set(i,j)!=0)
         break;
      
      

      // update x[i] and x[j], handle bounds carefully
      const double *Q_i = &Q[i*dim]; // get the i-th row of the Q matrix
      const double *Q_j = &Q[j*dim]; // get the j-th row of the Q matrix
      const double old_x_i = x[i];
      const double old_x_j = x[j];
      const double sum = x[i] + x[j];		// original sum (must maintain during update)
      const double P_val = std::max((double)(Q_i[i]+Q_j[j]-2*Q_i[j]), SMO::SMO_ZERO_EPS);	// it should be non-negative
      const double Q_val = (G[i]-G[j]) - old_x_i * P_val;				
      
      // this handle the case where P_val is -ve (in non PSD kernel matrix)
      // if you really want to speed things up then you can try to get rid of this bit
      if ( P_val <= SMO::SMO_ZERO_EPS )	 // non quadratic problem
      {
         if ( Q_val >= 0 )
            x[i] = 0;
         else 
            x[i] = sum;	// unbounded linear problem
      }
      else
      {
         x[i] = (-Q_val/P_val);
         if ( x[i] < 0 )
            x[i] = 0;
         else if ( x[i] > sum )
            x[i] = sum;
      }
      x[j] = (sum - x[i]);
      
      // update the x's status
      x_status[i] = (x[i] <= 0) ? SMO::LOWER_BOUND : SMO::FREE;
      x_status[j] = (x[j] <= 0) ? SMO::LOWER_BOUND : SMO::FREE;
      
      // update G
      double delta_x_i = x[i] - old_x_i;
      double delta_x_j = x[j] - old_x_j;
      for(int k=0; k < dim; k++)
         G[k] += Q_i[k]*delta_x_i + Q_j[k]*delta_x_j;
      
      // stop if little improvement (eps^2 is a very very small number)
      if ( SML::abs(delta_x_i) + SML::abs(delta_x_j) < smo_eps*smo_eps ) {
         //smo_eps = std::min(1e-15, smo_eps * 0.5);
         break;
      }
   }
   
   std::cout << "smo iter: " << smo_iter << std::endl;
   if(smo_iter >= maxSMOIter) {
      std::cout << "*************WARNING: maximum SMO iteration exceeded!" << std::endl;
   }
   if(smo_iter == 0) {
      std::cout << "-----------------WARNING: 0 SMO step! -------------------" << std::endl;
      smo_eps *= 0.5;
   }

   
   double gnorm2 = 0;
   for(int k=0; k < dim; k++) {
      gnorm2 += G[k]*G[i];
   }
   std::cout << "gnorm2: " << sqrt(gnorm2);
   
}


/**  SMO working set selection
 *   Steepest descent direction
 */
int CL2N2_SMO::smo_select_working_set(int &out_i, int &out_j)
{ 
   //double Gmax  = -SMO::INF;
   int Gmax_idx = -1;
   int Gmin_idx = -1;
   double obj_diff_min = SMO::INF;
   
   for(int i=0; i < dim; i++) 
   {
      const double *Q_i = 0;
      Q_i = &Q[i*dim];       
      
      for(int j=0; j < dim; j++)
       {
          if(j == i) continue;
          //if(x_status[j] != SMO::LOWER_BOUND)
          {
             double grad_diff = -G[i] + G[j];  //Gmax+G[j];
             double obj_diff; 
             double quad_coef=max(Q_i[i]+QD[j]-2*Q_i[j],SMO::TAU);
             if (grad_diff*grad_diff/quad_coef >= smo_eps)
             {
                if (quad_coef > 0)
                   obj_diff = -(grad_diff*grad_diff)/quad_coef;
                else
                   obj_diff = -(grad_diff*grad_diff)/SMO::TAU;
                if (obj_diff <= obj_diff_min)
                {
                   Gmax_idx = i;
                   Gmin_idx = j;
                   obj_diff_min = obj_diff;
                }
             }
          }
       }
    }


   if(Gmin_idx == -1)
      return 1;
   
   out_i = Gmax_idx;
   out_j = Gmin_idx;
   
   return 0;
}



/**  Update matrix and vectors needed in optimisation round.
 *
 *   \param a [read] Gradient
 *   \param b [read] Offset 
 */
void CL2N2_SMO::Update(TheMatrix& a, Scalar b)
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
   for(int i=0; i < maxGradSetSize; i++) {
      if(timeStamp[i] == -1 && first_new_slot == -1) first_new_slot = i;
      if(timeStamp[i] < timeStamp[oldest_slot] && timeStamp[i] >= 0) oldest_slot = i; 
      if(x[i] < x[smallest_slot]) smallest_slot = i;
   }
   if(timeStamp[oldest_slot] < iter - gradIdleAge) 
      idx = oldest_slot;
   else if(first_new_slot != -1) {
      idx = first_new_slot;
      dim++;
   }
   else 
      idx = smallest_slot;
   
   if(verbosity > 1) 
      std::cout << "new column idx: " << idx << "   gradientSet size:" << gradientSet.size() << std::endl;
   
   // store new column at #idx#   
   if(gradientSet[idx]) delete gradientSet[idx];
   gradientSet[idx] = new TheMatrix(a, gradType);    // store the gradient
   offsetSet[idx] = b;                   // store the offset  
   
   // use new slot
   if(timeStamp[idx] == -1)
   {
      // chteo: if this is not active after immediate optimization, we still make to make it current since it is new!
      timeStamp[idx] = iter;  
      
      // step 1: insert new row into Q
      int lastRow = (dim-1)*dim;
      for(int col=0; col < dim; col++) {
         Scalar value = 0;
         gradientSet[idx]->Dot(*gradientSet[col],value);
         Q[lastRow+col] = value/bmrmLambda;
      }
      
      // step 2: adjust memory area to fit a new column
      int prevDim = dim-1;
      int size    = sizeof(double)*prevDim;
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
      //x[dim-1] = 0;  // because i want to initialize x[0] to 1 in the very first iteration
      //x_status[dim-1] = SMO::LOWER_BOUND;
      f[dim-1] = offsetSet[dim-1];
   }
   // reuse slot
   else 
   { 
      // This gradient is new!
      timeStamp[idx] = iter;  
      
      // insert new row into Q[idx,:]
      int theRow = idx*dim;
      for(int col=0; col < dim; col++) {
         Scalar value = 0;
         gradientSet[idx]->Dot(*gradientSet[col],value);
         Q[theRow+col] = value/bmrmLambda;
      }
      
      // insert new column into Q[:,idx]
      for(int i=0; i < dim; i++)
         Q[i*dim + idx] = Q[theRow + i];
      
      // update vectors      
      x[idx] = 0;
      x_status[dim-1] = SMO::LOWER_BOUND;
      f[idx] = offsetSet[idx];
   }
   

   // compute G and QD
   for(int i=0; i < dim; i++) {
      G[i] = 0;
      double q_dot_x = 0;
      const double *q = &Q[i*dim];
      for(int j=0; j < dim; j++) {
         q_dot_x += q[j]*x[j];
      }
      G[i] = q_dot_x - f[i];
      QD[i] = q[i];
   }
   
#ifdef DEBUG
   // check Q matrix !    
   MatrixCorrectnessCheck();
#endif
   
}


/** Return solution.
 *
 *  \param w      [write] solution vector (preallocated)
 *  \param xi     [write] lower bound of empirical risk, R_emp[w]
 *  \param regval [write] value of regularization term
 *  \param objval [write] objective value
 */
void CL2N2_SMO::GetSolution(TheMatrix& w, Scalar &xi, Scalar &regval, Scalar &objval)
{
   assert(x != 0);
   double factor = 1.0/bmrmLambda;  // lambda is regularization constant
   
   
   // compute objective value (explicitly)
   // n := dimensionality of x (the solution of QP)
   // G := gradient of the QP objective
   // f := linear part of obj

   objval = 0;
   for(int i=0; i < dim; i++)
      objval += x[i]*(G[i]-f[i]);
   objval /= -2.0;

   if(verbosity > 1)
      std::cout << "inner qp obj: " << objval << std::endl;
   
   
   // Compute new w
   w.Zero();
   for(int i=0; i < dim; i++)
      if(x[i] > SML::ZERO_EPS)
	 w.ScaleAdd(-x[i], *gradientSet[i]);
   w.Scale(factor);
   
   // compute xi
   w.Norm2(regval);
   regval = regval * regval * bmrmLambda * 0.5; 
   xi = objval - regval;
   
   
   // update time-stamp of constraints
   for(int i=0; i < dim; i++)
      if(x[i] > 0) timeStamp[i] = iter;
   
   // constraint set size
   if(verbosity > 1)
      std::cout << "qp_size: " << dim << std::endl;
   
   
   // display nonzero terms
   if(verbosity > 2) {
      int howmany=0;
      std::cout << "nz: ";
      for(int i=0; i < dim; i++) {
         if(SML::abs(x[i])> SML::ZERO_EPS) {
            std::cout << i << "..";
            howmany++;
         }
      }
      std::cout << "   ("<<howmany<<")"<<std::endl;
   }
   
}

#endif
