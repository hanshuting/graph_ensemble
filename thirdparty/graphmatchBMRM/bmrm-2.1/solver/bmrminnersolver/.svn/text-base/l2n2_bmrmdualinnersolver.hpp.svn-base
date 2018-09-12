
#ifndef _L2N2_BMRMDUALINNERSOLVER_HPP_
#define _L2N2_BMRMDUALINNERSOLVER_HPP_

#include "model.hpp"
#include "sml.hpp"
#include "bmrminnersolver.hpp"


namespace L2N2_BMRMDualInnerSolver
{
        const double ZERO_EPS = 1e-20;
}


/** 
 *   when \Omega = 0.5|w|_2^2, the dual of the problem can be solved instead of the primal.
 */
class CL2N2_BMRMDualInnerSolver : public CBMRMInnerSolver 
{
   protected:

      /** Variables
       */
      double *x;

      /** Linear part of the objective function
       */
      double *f;

      /** Quadratic matrix of the objective function
       */
      double *Q;
      
      /** Constraint matrix
       */
      double *a;

      /** RHS of the (in)equality constraints
       */
      double *b;

      /** Lower bound vector for the variables
       */
      double *l;

      /** Upper bound vector for the variables
       */
      double *u;

      /** Tolerance for optimization error
       */
      double tol;
      
      /** Gradient set
       */
      vector<TheMatrix*> gradientSet;
      
      /** Offsets set
       */
      vector<double> offsetSet;

      /** Type of gradient to be store in gradientSet (e.g., SPARSE or DENSE)
       */
      int gradType;

      /** Time stamp for elements in aSet
       */
      vector<int> timeStamp;

      /** Number of consecutive iterations a gradient has a zero lagrange multiplier
       */
      int gradIdleAge;

      /** Maximum number of gradients to be keep in aSet at a particular time
       */
      int maxGradSetSize;



      /** Solve QP
       */
      virtual void SolveQP()=0;


      /** Update the solver with new gradient
       */
      virtual void Update(TheMatrix& a, Scalar b);

      /** Get solution and lower bound of empirical risk
       */
      virtual void GetSolution(TheMatrix& w, Scalar &xi, Scalar &regval, Scalar &objval);


#ifdef DEBUG
      /** Routine to check the correctness of the quadratic matrix of the objective function
       */
      void MatrixCorrectnessCheck();
#endif

   public:     

      /** Constructor
       */
      CL2N2_BMRMDualInnerSolver(double lambda);      

      /** Destructor
       */
      virtual ~CL2N2_BMRMDualInnerSolver();
      
      /** Solve the problem
       */
      virtual void Solve(TheMatrix& w, TheMatrix& a, Scalar loss, Scalar &xi, Scalar &regval, Scalar &objval);

      /** Reset the gradientSet, offsetSet, and mem for QP
       */
      virtual void Reset();
      
      /** With good QP tolerance annealing heuristic 
       *  the whole problem can be solved in much less number of iterations
       */
      virtual void SetTolerance(const double &theTolerance) {tol = theTolerance;}
};

#endif
