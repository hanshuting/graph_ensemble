/* 
 * Purpose      : L1-norm Linear Program solver using COIN-OR Clp 
 *
 */

#ifndef _L1N1_CLP_CPP_
#define _L1N1_CLP_CPP_

#include "l1n1_clp.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"

// Clp headers
#include "CoinBuild.hpp"


/*
    L1-Norm linear program:

    minimize   :  d'x_bar
     x_bar
     
    subject to : | -1_cm  A_mn  -A_mn | |\xi|    | 0|
                                        | v | <= | 0| 
                                        | u |    |-B|


    where      : d := [1; lambda; -lambda] (length: [1,n,n])
                 C is a positive scalar
                 lambda is a column vector of size n
                 x_bar := [\xi; v; u]
                 u,v column vectors of size n
                 1_c# is a column vector of 1 of size #,
                 A_#* is a dense matrix of size #-by-*, and
                 m << n (m starts from 1 and increase by 1 after 
                 each lp optimization round.)

*/


CL1N1_Clp::CL1N1_Clp(double lambda, const int &thedim)
   : CBMRMInnerSolver(lambda),
     newRowIndices(0),
     newRowElements(0)
{
  // parameters
  iter = 0;
  dim = thedim;
  numOfConstraint = 0;

  // sanity check
  assert(dim > 0);
  
  // build simplex model
  sim = new ClpSimplex();
  sim->setLogLevel(0);
  sim->resize(0, 1+2*dim);

  // set objective coefficients    
  sim->setObjCoeff(0, 1.0);
  for(int i=1; i<=dim; i++) 
     sim->setObjCoeff(i, bmrmLambda);
  for(int i=dim+1; i<=2*dim; i++)
     sim->setObjCoeff(i, bmrmLambda);

  // set bounds for variables (i.e., columns as in COIN Clp terminology)
  for (int i=0; i<2*dim+1; i++)
     sim->setColumnBounds(i, 0.0, COIN_DBL_MAX);
  
  
  // something useful in constraint matrix update phase
  newRowIndices = (int*)calloc(2*dim+1, sizeof(int));
  newRowElements = (double*)calloc(2*dim+1, sizeof(double));    
  
  for(int i=0; i<2*dim+1; i++)
     newRowIndices[i] = i;
  
  newRowElements[0] = -1.0;
  // Update() set the rest of the elements in newRowElements[]
    
}


CL1N1_Clp::~CL1N1_Clp()
{        
        if(sim) delete sim;
        if(newRowIndices) free(newRowIndices);
        if(newRowElements) free(newRowElements);  
}


void CL1N1_Clp::Solve(TheMatrix& w, TheMatrix& a, Scalar loss, Scalar &xi, Scalar &regval, Scalar &objval)
{
   Scalar w_dot_a = 0.0;
   w.Dot(a, w_dot_a);
   Update(a, loss - w_dot_a);

   //sim->dual();
   sim->primal();

   GetSolution(w, xi, regval, objval);
}


/** Clear constraint matrix
 */
void CL1N1_Clp::Reset()
{
        // delete rows in constraint matrix
        int nrows = sim->numberRows();
        int *rows = new int[nrows];
        for(int i=0; i<nrows; i++) rows[i] = i;
        sim->deleteRows(nrows, rows);      
        delete rows;
}


/**  Update matrix and vectors needed in optimisation round.
 *
 *   \param a [read] Gradient
 *   \param b [read] Offset 
 */
void CL1N1_Clp::Update(TheMatrix& a, Scalar b)
{
  iter++;
 
  // addRow: set values explicitly to accommodate different float types of gradient

  int length = a.Length();
  assert(length == dim);

  
  // resetting the content of newRowElements and newRowIndices
  memset(newRowElements, 0, sizeof(double)*(1+dim+dim));
  memset(newRowIndices, 0, sizeof(int)*(1+dim+dim));
  newRowElements[0] = -1;
  newRowIndices[0] = 0;

  Scalar val = 0;
  int nnz = 0;

  for(int i=0; i<length; i++)  {
     a.Get(i,val);
     if(SML::abs(val) > SML::ZERO_EPS) {
        newRowElements[1+nnz] = val;
        newRowIndices[1+nnz] = i+1;  // gradient starts from second column
        nnz++;
     }
  }

  if(verbosity > 0)
     std::cout << "nnz: " << nnz << "   length:" << dim << std::endl;

  for(int i=0; i<nnz; i++)  {
     newRowElements[1+nnz+i] = -newRowElements[1+i];
     newRowIndices[1+nnz+i]  = newRowIndices[1+i]+dim;
  }
  
  sim->addRow(1+nnz+nnz, newRowIndices, newRowElements, -COIN_DBL_MAX, -b);
  numOfConstraint++;
  
  if(verbosity > 0) {
     printf("num cols : %d\n", sim->numberColumns());
     printf("num rows : %d\n", sim->numberRows());
  }
}



/** Return solution.
 *
 *  \param w      [write] solution vector (preallocated)
 *  \param xi     [write] lower bound of empirical risk, R_emp[w]
 *  \param regval [write] value of regularization term
 *  \param objval [write] objective value
 */
void CL1N1_Clp::GetSolution(TheMatrix& w, Scalar &xi, Scalar &regval, Scalar &objval)
{
   double *solution = 0;
   solution = sim->primalColumnSolution(); 
   assert(solution);
  
   // construct new w
   int length = w.Length();
   for(int i=0; i<length; i++) {
      w.Set(i, solution[1+i] - solution[1+dim+i]);   // w := u-v
  }

  // compute xi, objective value and regularizer value
  xi = solution[0];
  objval = sim->objectiveValue();
  w.Norm1(regval);
  regval *= bmrmLambda;
}

#endif
