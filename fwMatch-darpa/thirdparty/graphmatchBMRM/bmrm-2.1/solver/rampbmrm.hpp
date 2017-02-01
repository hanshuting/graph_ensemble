// Created: 08/01/2008

#ifndef _RAMPBMRM_HPP_
#define _RAMPBMRM_HPP_

#include "common.hpp"
#include "solver.hpp"
#include "bmrminnersolver.hpp"
#include "model.hpp"
#include "loss.hpp"


/**   Class for BMRM for ramp losses
 *    This type of solver iteratively solves a Convex-Concave Procedure
 *    for nonconvex ramp losses.
 */
class CRampBMRM : public CSolver
{
protected:
        // program parameters
        int verbosity;           // level of program verbosity
        int convergenceLog;      // keep convergence log or not? (1:yes, 0:no)
        unsigned int maxNumOfOuterIter;   // maximum number of outer iteration
        unsigned int maxNumOfInnerIter;   // maximum number of inner iteration (column generation)
        Scalar epsilonTol;          // tolerance for epsilon termination criterion
        Scalar gammaTol;            // tolerance for gamma termination criterion
	Scalar tauTol;              // tolerance for CCCP stopping criterion i.e., |w_t - w_{t+1}|_2 < tau
        Scalar lambda;              // regularization constant
        string checkpointPrefix;    // prefix for intermediate model files
        unsigned int checkpointInterval;
        unsigned int checkpointMode;   
      
        enum CHECKPOINT_MODE {KEEP_ALL, KEEP_LATEST};

        // BMRM inner solver
        CBMRMInnerSolver* innerSolver;  // pointer to inner solver object
      
        // convex loss
        CLoss* _loss_vex;
      
        // linearization of concave loss
        CLoss* _loss_cav;
              
        // Methods
        virtual void ConfirmProgramParameters();
        
public:
        // Constructors
        CRampBMRM(CModel* model, CLoss* loss_vex, CLoss* loss_cav);

        // Destructor
        virtual ~CRampBMRM();

        // Methods
        virtual void Train();

};

#endif
