#ifndef _BMRM_HPP_
#define _BMRM_HPP_

#include "common.hpp"
#include "solver.hpp"
#include "bmrminnersolver.hpp"
#include "model.hpp"
#include "loss.hpp"


/**   Class for BMRM solver.
 *    This type of solver iteratively builds up a convex lower-bound of the 
 *      objective function, and performs minimization on the lower-bound.
 */
class CBMRM : public CSolver
{
   public:
      // Constructors
      CBMRM(CModel* model, CLoss* loss);

      // Destructor
      virtual ~CBMRM();

      // Methods
      virtual void Train();

   protected:
      // program parameters
      int    verbosity;         // level of program verbosity
      int    convergenceLog;    // keep convergence log or not? (1:yes, 0:no)
      unsigned int maxNumOfIter;      // maximum number of iteration (column generation)
      double epsilonTol;        // tolerance for epsilon termination criterion
      double gammaTol;          // tolerance for gamma termination criterion
      double lambda;            // regularization constant
      string checkpointPrefix;  // prefix for intermediate model files
      unsigned int checkpointInterval;
      unsigned int checkpointMode;   
      
      enum CHECKPOINT_MODE {KEEP_ALL, KEEP_LATEST};

      // BMRM inner solver
      CBMRMInnerSolver* innerSolver;  // pointer to inner solver object
      
      // Methods
      virtual void ConfirmProgramParameters();

};

#endif
