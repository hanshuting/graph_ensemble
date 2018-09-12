/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * Library of solvers for QP task required in StructSVM learning.
 *
 * Written (W) 1999-2007 Vojtech Franc, xfrancv@cmp.felk.cvut.cz
 * Copyright (C) 1999-2007 Center for Machine Perception, CTU FEL Prague 
-------------------------------------------------------------------- 
Synopsis:

  exitflag = qpssvm_solver( &get_col, diag_H, f, b, I, x, n, tmax, 
             tolabs, tolrel, &t, &History, verb );   

  exitflag = qpssvm_solver( &get_col, diag_H, f, b, I, x, n, tmax, 
             tolabs, tolrel, &QP, &QD, verb );   
Description:
 
 It solves the following QP task:
  
   min 0.5*x'*H*x + f'*x
    x

 subject to 
 
   sum(x(find(I==k))) <= b   for all k=1:max(I)
   x >= 0

 where I is a set of positive indices from (1 to max(I)).

 A precision of the found solution is given by the parameters tmax, 
 tolabs and tolrel which define the stopping conditions:
 
 UB-LB <= tolabs      ->  exitflag = 1   Abs. tolerance.
 UB-LB <= UB*tolrel   ->  exitflag = 2   Relative tolerance.
 t >= tmax            ->  exitflag = 0   Number of iterations.

 UB ... Upper bound on the optimal solution, i.e., Q_P.
 LB ... Lower bound on the optimal solution, i.e., Q_D.
 t  ... Number of iterations.


Inputs/Outputs:

 const void* (*get_col)(unsigned int) retunr pointer to i-th column of H
 diag_H [double n x n] diagonal of H.
 f [double n x 1] is an arbitrary vector.
 b [double 1 x 1] scalar
 I [unsigned int n x 1] Indices (1..max(I)); max(I) <= n
 x [double n x 1] solution vector (inital solution).
 n [unsigned int 1 x 1] dimension of H.
 tmax [unsigned int 1 x 1] Max number of steps.
 tolrel [double 1 x 1] Relative tolerance.
 tolabs [double 1 x 1] Absolute tolerance.
 t [unsigned int 1 x 1] Number of iterations.
 History [double 2 x t] Value of LB and UB wrt. number of iterations.
 verb [int 1 x 1] if > 0 then prints info every verb-th iteation.

 For more info refer to TBA

 Modifications:
 01-Oct-2007, VF
 20-Feb-2006, VF
 18-feb-2006, VF

-------------------------------------------------------------------- */

// #include <math.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include <string.h>
// #include <stdint.h>
// #include <limits.h>

// #include "classifier/svm/libocas_common.h"
// #include "classifier/svm/qpssvmlib.h"




// chteo: 281207:1649
// fitting Vojtech Franc's GMNP into bmrm...

#ifndef _L2N2_GMN_CPP_
#define _L2N2_GMN_CPP_

#include "configuration.hpp"
#include "l2n2_gmn.hpp"

#define INDEX2(ROW,COL,NUM_ROWS) ((COL)*(NUM_ROWS)+(ROW))
#define MIN(A,B) ((A) > (B) ? (B) : (A))
#define MAX(A,B) ((A) < (B) ? (B) : (A))
#define ABS(A) ((A) < 0 ? -(A) : (A))



CL2N2_GMN::CL2N2_GMN(double lambda)
   : CL2N2_BMRMDualInnerSolver(lambda)
{   
   bmrmC = 1.0;
   tmax = 10000000;
   qpsolvertolrel = 1e-3;
   
   Configuration &config = Configuration::GetInstance();
   
   if(config.IsSet("L2N2_GMN.maxIter"))
   {
      tmax = config.GetInt("L2N2_GMN.maxIter");
   }
   
   if(config.IsSet("BMRM.epsilonTol"))
   {
      double epstol = config.GetDouble("BMRM.epsilonTol");
      if(epstol > 0.0)
	 qpsolvertolrel = epstol*0.1;
   }
   else if(config.IsSet("BMRM.gammaTol"))
   {
      double gamtol = config.GetDouble("BMRM.gammaTol");
      if(gamtol > 0.0)
	 qpsolvertolrel = gamtol*0.1;
   }
   
   assert(maxGradSetSize > 0 and maxGradSetSize < 100000);
   diag_Q = (double*)calloc(maxGradSetSize, sizeof(double));
   I = (unsigned int*)calloc(maxGradSetSize, sizeof(unsigned int));
   
   for(int i=0; i<maxGradSetSize; i++)
      I[i] = 1;  // because we have only one type of constraint matrix

   if(verbosity > 0)
   {
      std::cout << "CL2N2+GMN instantiated!" << std::endl;
      std::cout << "  tmax: " << tmax << std::endl;

   }
   std::cout << "  qpsolvertolrel: " << qpsolvertolrel << std::endl;
}


CL2N2_GMN::~CL2N2_GMN()
{
   if(diag_Q) free(diag_Q);
   if(I) free(I);
}


/**  Update matrix and vectors needed in optimisation round.
 *
 *   \param a [read] Gradient
 *   \param b [read] Offset 
 */
void CL2N2_GMN::Update(TheMatrix& a, Scalar b)
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
   //for(int i=0; i < maxGradSetSize; i++) {
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
   
   if(verbosity > 2) 
      std::cout << "new column idx: " << idx << "   gradientSet size:" << gradientSet.size() << std::endl;
   
   // store new column at #idx#   
   if(gradientSet[idx]) delete gradientSet[idx];
   gradientSet[idx] = new TheMatrix(a, gradType);    // store the gradient
   offsetSet[idx] = b;                   // store the offset  


//    for(int i=0; i<dim; i++)
//    {
//       for(int j=i; j<dim; j++)
//       {
// 	 Scalar value = 0.0;
// 	 gradientSet[i]->Dot(*gradientSet[j],value);
// 	 Q[i*dim+j] = 1.0/bmrmLambda*value;
// 	 Q[i+dim*j] = 1.0/bmrmLambda*value;
//       }
//    }
   
//    timeStamp[idx] = iter;  
//    x[idx] = 0;
//    f[idx] = -offsetSet[idx];




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
      for(int col=0; col < dim; col++) {
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
      
}


/** Return solution.
 *
 *  \param w      [write] solution vector (preallocated)
 *  \param xi     [write] lower bound of empirical risk, R_emp[w]
 *  \param regval [write] value of regularization term
 *  \param objval [write] objective value
 */
void CL2N2_GMN::GetSolution(TheMatrix& w, Scalar &xi, Scalar &regval, Scalar &objval)
{
   assert(x != 0);
   double factor = 1.0/bmrmLambda;  // lambda is regularization constant
   
   
   // compute objective value (explicitly)
   // n := dimensionality of x (the solution of QP)
   // Q := hessian matrix
   // f := linear part of obj

   objval = 0;
   double tmp = 0;
   double fx = 0;
   for(int i=0; i < dim; i++) {
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
      if(x[i] > SML::ZERO_EPS)
         w.ScaleAdd(-x[i], *gradientSet[i]);
   w.Scale(factor);
   
   // compute xi
   w.Norm2(regval);
   regval = regval * regval * bmrmLambda * 0.5; 
   xi = objval - regval;
   
   
   // update time-stamp of constraints
   for(int i=0; i < dim; i++)
      if(x[i] > SML::ZERO_EPS) 
         timeStamp[i] = iter;
   
   // constraint set size
   if(verbosity >= 1)
   {
      int howmany = 0;
      std::cout << "qp_size: " << dim << std::endl;
      for(int i=0; i < dim; i++) {
         if(SML::abs(x[i])> SML::ZERO_EPS) {
            howmany++;
         }
      }
      std::cout << "gradient nnz: " << howmany << std::endl;
   }
   
   
   // display nonzero terms
   if(verbosity >= 3) {
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



void CL2N2_GMN::SolveQP()
{      
   for(int i=0; i<dim; i++)
      diag_Q[i] = Q[i*dim+i];
   double RHS = bmrmC; // RHS of linear constraint
   double tolabs = 0.0;  // following ocas
   qp_objval = 0.0;
   double dummy = 0.0;
   int exitflag = 0;
   exitflag = qpssvm_solver(diag_Q, f, RHS, I, x, dim, tmax, tolabs, tol*0.5, &qp_objval, &dummy, verbosity);

   if(verbosity >= 1)
      printf("qpssvm exitflag=%d\t qp_bojval=%f\n",exitflag, qp_objval);
      
}

int CL2N2_GMN::qpssvm_solver(double *diag_H,
			     double *f,
			     double b,
			     unsigned int *I,
			     double *x,
			     unsigned int n,
			     unsigned int tmax,
			     double tolabs,
			     double tolrel,
			     double *QP,
			     double *QD,
			     unsigned int verb)
{
  double *x_nequ;
  double *d;
  double *col_u, *col_v;
  double LB;
  double UB;
  double tmp;
  double improv;
  double tmp_num;
  double tmp_den=0;
  double tau=0;
  double delta;
  double yu;
  unsigned int *inx;
  unsigned int *nk;
  unsigned int m;
  unsigned int t;
  unsigned int u=0;
  unsigned int v=0;
  unsigned int k;
  unsigned int i, j;
  int exitflag;

  
  /* ------------------------------------------------------------ 
    Initialization                                               
  ------------------------------------------------------------ */

  x_nequ=NULL;
  inx=NULL;
  nk=NULL;
  d=NULL;

  /* count cumber of constraints */
  for( i=0, m=0; i < n; i++ ) m = MAX(m,I[i]);

  /* alloc and initialize x_nequ */
  x_nequ = (double*) calloc(m, sizeof(double));
  if( x_nequ == NULL )
  {
	  exitflag=-2;
	  goto cleanup;
  }

  /* alloc Inx */
  inx = (unsigned int*) calloc(m*n, sizeof(unsigned int));
  if( inx == NULL )
  {
	  exitflag=-2;
	  goto cleanup;
  }

  nk = (unsigned int*) calloc(m, sizeof(unsigned int));
  if( nk == NULL )
  {
	  exitflag=-2;
	  goto cleanup;
  }

  for( i=0; i < m; i++ ) x_nequ[i] = b;
  for( i=0; i < n; i++ ) {
     k = I[i]-1;
     x_nequ[k] -= x[i];
     inx[INDEX2(nk[k],k,n)] = i;
     nk[k]++;
  }
    
  /* alloc d [n x 1] */
  d = (double*) calloc(n, sizeof(double));
  if( d == NULL )
  {
	  exitflag=-2;
	  goto cleanup;
  }
 
  /* d = H*x + f; */
  for( i=0; i < n; i++ ) {
    if( x[i] > 0 ) {

       //chteo:
       //col_u = (double*)get_col(i);
       col_u = &Q[i*dim];

      for( j=0; j < n; j++ ) {
          d[j] += col_u[j]*x[i];
      }
    }
  }
  for( i=0; i < n; i++ ) d[i] += f[i];
  
  /* UB = 0.5*x'*(f+d); */
  /* LB = 0.5*x'*(f-d); */
  for( i=0, UB = 0, LB=0; i < n; i++) {
    UB += x[i]*(f[i]+d[i]);
    LB += x[i]*(f[i]-d[i]);
  }
  UB = 0.5*UB;
  LB = 0.5*LB;

  /*
  for k=1:m,
    tmp = min(d(find(I==k)));
    if tmp < 0, LB = LB + b*tmp; end
  end
  */
  
  for( i=0; i < m; i++ ) {
    for( j=0, tmp = GMN::INF; j < nk[i]; j++ ) {
      tmp = MIN(tmp, d[inx[INDEX2(j,i,n)]]);
    }
    if( tmp < 0) LB += b*tmp;
  }
  
  exitflag = 0;
  t = 0;

  /* -- Main loop ---------------------------------------- */
  while( (exitflag == 0) && (t < tmax)) 
  {
    t++;

    exitflag = 1;
    for( k=0; k < m; k++ ) 
    {       
      /*
      inx = find(I==k);
      [tmp,u] = min(d(inx)); u = inx(u);
      */
        
     for( j=0, tmp = GMN::INF, delta = 0; j < nk[k]; j++ ) {
        i = inx[INDEX2(j,k,n)];
        delta += x[i]*d[i];
        if( tmp > d[i] ) {
          tmp = d[i];
          u = i;
        }
      }

      /* if d(u) < 0, yu = b; else yu = 0; end  */
      if( d[u] < 0) yu = b; else yu = 0;
     
      /* delta = x(inx)'*d(inx) - yu*d(u); */
      delta -= yu*d[u];
            
      if( delta > tolabs/m && delta > tolrel*ABS(UB)/m) 
      {
         exitflag = 0;
         
         if( yu > 0 ) 
         {
	    //chteo:
	    //col_u = (double*)get_col(u);      
	    col_u = &Q[u*dim];

           improv = -GMN::INF;
           for( j=0; j < nk[k]; j++ ) {
             i = inx[INDEX2(j,k,n)];
           
/*           for(i = 0; i < n; i++ ) {
             if( (I[i]-1 == k) && (i != u) && (x[i] > 0)) {              */
             if(x[i] > 0) {             
               
               tmp_num = x[i]*(d[i] - d[u]); 
               tmp_den = x[i]*x[i]*(diag_H[u] - 2*col_u[i] + diag_H[i]);
               if( tmp_den > 0 ) {
                 if( tmp_num < tmp_den ) {
                    tmp = tmp_num*tmp_num / tmp_den;
                 } else {
                    tmp = tmp_num - 0.5 * tmp_den;
                 }
               }
               if( tmp > improv ) {
                 improv = tmp;
                 tau = MIN(1,tmp_num/tmp_den);
                 v = i;
               }
             }
           }

           tmp_num = -x_nequ[k]*d[u];
           if( tmp_num > 0 ) {
             tmp_den = x_nequ[k]*x_nequ[k]*diag_H[u];
             if( tmp_den > 0 ) {
               if( tmp_num < tmp_den ) {
                 tmp = tmp_num*tmp_num / tmp_den;
               } else {
                   tmp = tmp_num - 0.5 * tmp_den;
               }
             }
           } else {
             tmp = -GMN::INF; 
           }
           
           if( tmp > improv ) {
              tau = MIN(1,tmp_num/tmp_den);
              for( i = 0; i < n; i++ ) {             
                d[i] += x_nequ[k]*tau*col_u[i];
              }
             x[u] += tau*x_nequ[k];
             x_nequ[k] -= tau*x_nequ[k];
               
           } else {
            
             /* updating with the best line segment */
	      // chteo:
	      //col_v = (double*)get_col(v);
	      col_v = &Q[v*dim];

             for( i = 0; i < n; i++ ) {             
               d[i] += x[v]*tau*(col_u[i]-col_v[i]);
             }

             x[u] += tau*x[v];
             x[v] -= tau*x[v];
           }
         }
         else
         {
           improv = -GMN::INF;
           for( j=0; j < nk[k]; j++ ) {
             i = inx[INDEX2(j,k,n)];
           
/*           for(i = 0; i < n; i++ ) {
             if( (I[i]-1 == k) && (x[i] > 0)) {*/
             if( x[i] > 0 && d[i] > 0) {
                
               tmp_num = x[i]*d[i]; 
               tmp_den = x[i]*x[i]*diag_H[i];
               if( tmp_den > 0 ) {
                 if( tmp_num < tmp_den ) {
                    tmp = tmp_num*tmp_num / tmp_den;
                 } else {
                    tmp = tmp_num - 0.5 * tmp_den;
                 }
               }
               if( tmp > improv ) {
                 improv = tmp;
                 tau = MIN(1,tmp_num/tmp_den);
                 v = i;
               }
             }    
           }

           /* updating with the best line segment */
           //chteo:
	   //col_v = (double*)get_col(v);
	   col_v = &Q[v*dim];
	   
           for( i = 0; i < n; i++ ) {             
             d[i] -= x[v]*tau*col_v[i];
           }

           x_nequ[k] += tau*x[v];
           x[v] -= tau*x[v];         
         }

         UB = UB - improv;
      }
                   
    }

    /* -- Computing LB --------------------------------------*/

    /*
    LB = 0.5*x'*(f-d);   
    for k=1:n,
      LB = LB + b*min(d(find(I==k)));
    end */
    
    for( i=0, UB = 0, LB=0; i < n; i++) {
       UB += x[i]*(f[i]+d[i]);
       LB += x[i]*(f[i]-d[i]);
    }
    UB = 0.5*UB;
    LB = 0.5*LB;

    for( k=0; k < m; k++ ) { 
      for( j=0,tmp = GMN::INF; j < nk[k]; j++ ) {
        i = inx[INDEX2(j,k,n)];

        tmp = MIN(tmp, d[i]);
      }
      if( tmp < 0) LB += b*tmp;
    }

    //if( verb > 0 && (exitflag > 0 || (t % verb)==0 )) {
    if(verbosity >= 2) {
       printf("%d: UB=%.10f, LB=%.10f, UB-LB=%.10f, (UB-LB)/|UB|=%.10f \n",
        t, UB, LB, UB-LB, (UB!=0) ? (UB-LB)/ABS(UB) : 0);      
    }    

  }

  /* -- Find which stopping consition has been used -------- */
  if( UB-LB < tolabs ) exitflag = 1;
  else if(UB-LB < ABS(UB)*tolrel ) exitflag = 2;
  else exitflag = 0;

  /*----------------------------------------------------------   
    Set up outputs                                          
  ---------------------------------------------------------- */
  *QP = UB;
  *QD = LB;

  /*----------------------------------------------------------
    Clean up
  ---------------------------------------------------------- */
cleanup:
  free( d );
  free( inx );
  free( nk );
  free( x_nequ );  
  
  return( exitflag ); 

}

#endif
