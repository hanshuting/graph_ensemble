/* Notes :
 *
 *  Code adapted from CVM
 */

#ifndef _L2N2_SMO_HPP_
#define _L2N2_SMO_HPP_

#include "sml.hpp"
#include "model.hpp"
#include "l2n2_bmrmdualinnersolver.hpp"
#include <iostream>
#include <cassert>


namespace SMO {
   const double INF = 1e30;
   const double SMO_ZERO_EPS = 1e-2;
   enum {LOWER_BOUND, FREE};
   const double TAU = 1e-10;
}

class CL2N2_SMO : public CL2N2_BMRMDualInnerSolver 
{
   protected: 
      int maxSMOIter;
      double smo_eps;
      double smo_eps_sq;
      double *G, *QD;  // G is the gradient of the QP obj; QD is the diag of Q matrix
      int *x_status;

      int smo_select_working_set(int &i, int &j);

      void SolveQP();
      void Update(TheMatrix& a, Scalar b);
      void GetSolution(TheMatrix& w, Scalar &xi, Scalar &regval, Scalar &objval);

   public:      
      CL2N2_SMO(double lambda);
      
      virtual ~CL2N2_SMO();
};

#endif
