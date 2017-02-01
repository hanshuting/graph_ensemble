
#ifndef _L2N2_GMN_HPP_
#define _L2N2_GMN_HPP_

#include "sml.hpp"
#include "model.hpp"
#include "l2n2_bmrmdualinnersolver.hpp"
#include <iostream>
#include <cassert>



namespace GMN
{
   const double INF = 1e30;
}

class CL2N2_GMN : public CL2N2_BMRMDualInnerSolver 
{
   public:      
      CL2N2_GMN(double lambda);      
      virtual ~CL2N2_GMN();
      
      virtual void SolveQP();
      virtual void Update(TheMatrix& a, Scalar b);
      virtual void GetSolution(TheMatrix& w, Scalar &xi, Scalar &regval, Scalar &objval);
      
   private: 
      unsigned int tmax;
      double qpsolvertolrel;
      double qp_objval;
      double *diag_Q;
      unsigned int *I;
      Scalar bmrmC;
      
      int qpssvm_solver(double *diag_H,
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
			unsigned int verb);
      
};

#endif
