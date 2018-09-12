
/* 
 * Purpose      : L1-norm Linear Program solver using COIN-OR Clp 
 *
 */

#ifndef _L1N1_CLP_HPP_
#define _L1N1_CLP_HPP_

#include "sml.hpp"
#include "bmrminnersolver.hpp"

// from COIN-OR Clp
#include "ClpSimplex.hpp"

/* 
 * solve the 1-norm linear program
 *
 * minimize    c' * x + |x|_1
 * subject to  r <= A*x <= r+b
 *             l <= x <= u
 *
 */
class CL1N1_Clp : public CBMRMInnerSolver 
{
protected:

        /** CLP simplex solver
         */
        ClpSimplex *sim;
      
        /** Pre-allocated indices array for new constaint
         */
        int *newRowIndices;
      
        /** Pre-allocated values array for new constraint
         */
        double *newRowElements;

        /** Update the constraint matrix
        */
        virtual void Update(TheMatrix& a, Scalar b);

        /** Return the updated solution
         */
        virtual void GetSolution(TheMatrix& w, Scalar &xi, Scalar &regval, Scalar &objval);

public:      

        /** Constructor
         */
        CL1N1_Clp(double lambda, const int &thedim);
      
        /** Destructor
         */
        virtual ~CL1N1_Clp();
      
        /** Solve the problem
         */
        virtual void Solve(TheMatrix& w, TheMatrix& a, Scalar loss, Scalar &xi, Scalar &regval, Scalar &objval);
        
        /** Clear constraint matrix
         */
        virtual void Reset();
};

#endif
