//created: 08/01/2008
#ifndef _RAMPBMRM_CPP_
#define _RAMPBMRM_CPP_

#include "common.hpp"
#include "rampbmrm.hpp"
#include "timer.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "bmrminnersolver.hpp"
#include "loss.hpp"
#include "bmrminnersolverfactory.hpp"

#include <fstream>
#include <sstream>


/** Constructor
 *
 *  @param model [read] pointer to loss model object
 *  @param loss_vex [read] pointer to convex loss
 *  @param loss_cav [read] pointer to linearization of concave loss
 */
CRampBMRM::CRampBMRM(CModel *model, CLoss* loss_vex, CLoss* loss_cav)
        : CSolver(model, 0),
          _loss_vex(loss_vex),
          _loss_cav(loss_cav)
{
        // set private members (default) values
        verbosity         = 0;
        convergenceLog    = 0;
        maxNumOfOuterIter = 500;
        maxNumOfInnerIter = 10000;
        epsilonTol        = 1e-4;
        gammaTol         = -1e15;  // negative value means this criteria is off
	tauTol           = 1e-2;
        lambda           = 1e-3;
        checkpointPrefix = "model.checkpoint";
        checkpointInterval = 1000000;  // no checkpoint by default
        checkpointMode = KEEP_LATEST;

        // check and confirm program parameters
        ConfirmProgramParameters();

        // instantiate inner solver
        innerSolver = CBMRMInnerSolverFactory::GetBMRMInnerSolver(*_model, lambda);
}


/**  Destructor
 */
CRampBMRM::~CRampBMRM()
{
        // destroy inner solver
        if(innerSolver) delete innerSolver;
}


/**  Start training/learning a model w.r.t. the loss object (and the data supplied to it).
 */
void CRampBMRM::Train()
{

        // Timers (CPU and wall-clock)
        CTimer totalTime;             // total runtime of the training
        CTimer innerSolverTime;       // time for inner optimization (e.g., QP or LP)
        CTimer lossAndGradientTime;   // time for loss and gradient computation

        unsigned int outerIter = 0;        
        unsigned int innerIter = 0;         // iteration count                
        unsigned int totalInnerIter = 0;
        Scalar xi = 0.0;              // convex lower-bound (approximate) of loss value
        Scalar exactObjVal = 0.0;     // (exact) objective function value
        Scalar approxObjVal = 0.0;    // convex lower-bound (approximate) of objective function value
        Scalar minExactObjVal = SML::INFTY; // minimum of all previously evaluated (exact) objective function value
        Scalar regVal = 0.0;          // value of the regularizer term e.g., 0.5*w'*w
        Scalar epsilon = 0.0;         // (duality) gap := exactObjVal - approxObjVal
        Scalar gamma = 0.0;           // := minExactObjVal - approxObjVal
        Scalar tau = 0.0;             // := rampLossVal(t) - rampLossVal(t+1)
        Scalar w_dot_a = 0.0;         // temp for <w,a> where 'w' is weight vector and 'a' is gradient
        double innerSolverTol = 1.0;  // optimization tolerance for inner solver

        unsigned int row = 0; 
        unsigned int col = 0;
        TheMatrix &w = _model->GetW();
        w.Shape(row, col);   
        
        TheMatrix w_old(row, col, SML::DENSE);    // keep w at previous outer iter
        TheMatrix a_vex(row, col, SML::DENSE);    // gradient vector for convex loss
        TheMatrix a_cav(row, col, SML::DENSE);    // gradient vector for linearization of concave loss
                
        Scalar lossVal_vex = 0.0;                 // convex loss function value        
        Scalar lossVal_cav = 0.0;                 // value of linearization of concave loss function
        Scalar rampLossVal = 0.0;                 // := lossVal_vex - lossVal_cav
        Scalar prevRampLossVal = SML::INFTY;
        
        // start training
        totalTime.Start();

        // loop for CCCP 
        Scalar curEpsilonTol = epsilonTol;//*10;  // decrease epsilon tolerance gradually
        Scalar curGammaTol = gammaTol;//*10;                
        w_old.Zero();
        a_cav.Zero();
        
        while(1)
        {
                outerIter++;
                std::cout << "\nOuter Iter: " << outerIter << std::endl;
                
                // loop for solving convex + linear loss
                innerIter = 0;         // inner iteration count                
                xi = 0.0;              // convex lower-bound (approximate) of loss value
                exactObjVal = 0.0;     // (exact) objective function value
                approxObjVal = 0.0;    // convex lower-bound (approximate) of objective function value
                minExactObjVal = 1e99; // minimum of all previously evaluated (exact) objective function value
                regVal = 0.0;          // value of the regularizer term e.g., 0.5*w'*w
                epsilon = 0.0;         // (duality) gap := exactObjVal - approxObjVal
                gamma = 0.0;           // := minExactObjVal - approxObjVal
                w_dot_a = 0.0;         // temp for <w,a> where 'w' is weight vector and 'a' is gradient
                innerSolverTol = 1.0;  // optimization tolerance for inner solver
                
                // reset inner solver and weight vector for newly updated "convex + linear" function
                innerSolver->Reset();
                w.Zero();
                
                while(1)
                {
                        innerIter++;
      
                        // column generation
                        lossAndGradientTime.Start();
                        _loss_vex->ComputeLossAndGradient(lossVal_vex, a_vex);
                        w.Dot(a_cav,lossVal_cav);
                        rampLossVal = lossVal_vex - lossVal_cav;                        
                        a_vex.Minus(a_cav);       
                        lossAndGradientTime.Stop();

                        // update convergence monitor                        
			exactObjVal = regVal + rampLossVal;
                        minExactObjVal = std::min(minExactObjVal, exactObjVal);
                        epsilon = exactObjVal - approxObjVal;
                        gamma = minExactObjVal - approxObjVal;
                                                
                        // dump convergence statistics on stdout
                        if(verbosity < 1) 
                        {
                                printf(".");
                                if(innerIter%100 == 0) 
                                        printf("%d",innerIter);
                                fflush(stdout);
                        }
                        else if(verbosity == 1)
                                printf("#%d   eps %.6e   loss %.6e   xi %.6e   reg %.6e\n",innerIter, epsilon, rampLossVal, xi, regVal);      
                        else if(verbosity > 1)
                                printf("#%d   pobj %.6e (%6e)   aobj %.6e   eps %.6e   gam %.6e   loss %.6e   xi % .6e   reg %.6e\n",
                                        innerIter, exactObjVal, minExactObjVal, approxObjVal, epsilon, gamma, rampLossVal, xi, regVal);      

                        // stopping criteria for solving convex + linear loss
                        if((innerIter >= 2) && ((gamma < curGammaTol) || (epsilon < curEpsilonTol)))
                                break;
                        if(innerIter >= maxNumOfInnerIter)
                        { 
                                printf("\nWARNING: inner loop exceeded maximum number of iterations (%d) !\n", maxNumOfInnerIter);
                                break;
                        }
                        
                        // adjust inner solver optimization tolerance
                        innerSolverTol = std::min(innerSolverTol, std::max((double)epsilon, (double)curEpsilonTol));
                        innerSolver->SetTolerance(innerSolverTol*0.5);

                        // run inner solver
                        innerSolverTime.Start();
                        innerSolver->Solve(w, a_vex, rampLossVal, xi, regVal, approxObjVal);
                        innerSolverTime.Stop();
                }  //end of inner loop
                totalInnerIter += innerIter;                
                
                Scalar diff_w_norm2 = 0.0;
                w_old.Minus(w);
                w_old.Norm2(diff_w_norm2);
                w_old.Assign(w);
                                                                
                // compute the linearization of concave loss
                lossAndGradientTime.Start();
                Scalar dummy = 0.0;
                _loss_cav->ComputeLossAndGradient(dummy, a_cav);
                lossAndGradientTime.Stop();
                                			
                // check point
                if(outerIter%checkpointInterval == 0) 
                {
                        if(checkpointMode == KEEP_LATEST)
                                _model->Save(checkpointPrefix);
                        else 
                        {
                                ostringstream oss;
                                oss << checkpointPrefix << "." << outerIter;
                                _model->Save(oss.str());
                        }
                }
                
                // decrease epsilon/gamma tolerance
                curEpsilonTol = std::max((Scalar)(curEpsilonTol*0.5), epsilonTol);
                curGammaTol = std::max((Scalar)(curGammaTol*0.5), gammaTol);
                
                // update CCCP convergence monitor
                tau = fabs(prevRampLossVal - rampLossVal);
                prevRampLossVal = rampLossVal;
                
                if(verbosity < 1) std::cout << std::endl;
                std::cout << "rampLossVal: " << rampLossVal << std::endl;
                std::cout << "epsilon: " << epsilon << std::endl;                
                std::cout << "gamma: " << gamma << std::endl;                
                std::cout << "tau: " << tau << std::endl;
                std::cout << "|w-w_old|_2: " << diff_w_norm2 << std::endl;
                
                // stopping criteria for CCCP
                if((outerIter >= 2) && (tau < tauTol))
                        break;
                if(outerIter >= maxNumOfOuterIter) 
                {
                        printf("\nWARNING: outer loop exceeded maximum number of iterations (%d) !\n", maxNumOfOuterIter);
                        break;
                }                
        } 
            
        // legends
        if(verbosity >= 1) 
        {
                printf("\nLegends::\n");
                if(verbosity > 1)
                        printf("pobj: primal objective function value\n"
							   "aobj: approximate objective function value\n");
                printf("eps: epsilon (approximation error) \n"
					   "gam: lower bound on eps \n"
					   "loss: loss function value \n"
					   "xi: approximation to loss \n"
					   "reg: regularizer value\n");
        }
      

        Scalar norm1 = 0, norm2 = 0, norminf = 0;
        w.Norm1(norm1);
        w.Norm2(norm2);
        w.NormInf(norminf);
        printf("\n");
        printf("No. of outer ite.:  %d\n", outerIter);
        printf("No. of inner ite.:  %d\n", totalInnerIter);
        printf("Primal obj. val. : %.6e\n",exactObjVal);
        printf("Approx obj. val. : %.6e\n",approxObjVal);
        printf("Primal - Approx. : %.6e\n",exactObjVal-approxObjVal);
        printf("Loss:              %.6e\n",lossVal_vex);
        printf("|w|_1:             %.6e\n",norm1);
        printf("|w|_2:             %.6e\n",norm2);
        printf("|w|_oo:            %.6e\n",norminf);
   
        totalTime.Stop();
        // end of training

        // display timing profile
        printf("\nCPU seconds in:\n");
        printf("1. loss and gradient: %8.2f\n", lossAndGradientTime.CPUTotal());
        printf("2. solver:            %8.2f\n", innerSolverTime.CPUTotal()); 
        printf("               Total: %8.2f\n", totalTime.CPUTotal());

}


/**   Validate program parameters set in Configuration.
 */
void CRampBMRM::ConfirmProgramParameters()
{
        Configuration &config = Configuration::GetInstance();  // make sure configuration file is read before this!

        if(config.IsSet("RampBMRM.verbosity")) 
        {
                verbosity = config.GetInt("RampBMRM.verbosity");
        }
   
        if(config.IsSet("RampBMRM.convergenceLog")) 
        {
                convergenceLog = config.GetInt("RampBMRM.convergenceLog");
        }
   
        if(config.IsSet("RampBMRM.maxNumOfOuterIter")) 
        {
                maxNumOfOuterIter = config.GetInt("RampBMRM.maxNumOfOuterIter");
                if(maxNumOfInnerIter < 0)
                throw CBMRMException("RampBMRM.maxNumOfOuterIter must be > 0\n","CRampBMRM::ConfirmProgramParameters()");                
        }
   
        if(config.IsSet("RampBMRM.maxNumOfInnerIter")) 
        {
                maxNumOfInnerIter = config.GetInt("RampBMRM.maxNumOfInnerIter");
                if(maxNumOfInnerIter < 0)
                throw CBMRMException("RampBMRM.maxNumOfInnerIter must be > 0\n","CRampBMRM::ConfirmProgramParameters()");
        }
        
        // One of the BMRM stopping criterion
        if(config.IsSet("RampBMRM.epsilonTol")) 
        {
                epsilonTol = config.GetDouble("RampBMRM.epsilonTol");
        }

        // One of the BMRM stopping criterion
        if(config.IsSet("RampBMRM.gammaTol")) 
        {
                gammaTol = config.GetDouble("RampBMRM.gammaTol");
        }
		
        // The CCCP stopping criteria
	if(config.IsSet("RampBMRM.tauTol")) 
        {
			tauTol = config.GetDouble("RampBMRM.tauTol");
	}
        
        if(config.IsSet("RampBMRM.lambda"))           
        {
                lambda = config.GetDouble("RampBMRM.lambda");
                if(lambda <= 0)
                        throw CBMRMException("RampBMRM.lambda must be > 0\n","CRampBMRM::ConfirmProgramParameters()");
        }
        else 
        {
                config.SetDouble("RampBMRM.lambda", lambda);
        }

        if(config.IsSet("RampBMRM.checkpointInterval")) 
        {
                checkpointInterval = config.GetInt("RampBMRM.checkpointInterval");
                if(checkpointInterval < 1)
                        throw CBMRMException("RampBMRM.checkpointInterval must be a positive integer!\n","CRampBMRM::ConfirmProgramParameters()");
        }

        if(config.IsSet("RampBMRM.checkpointPrefix")) 
        {
                checkpointPrefix = config.GetString("RampBMRM.checkpointPrefix");
        }

        if(config.IsSet("RampBMRM.checkpointMode")) 
        {
                string mode = config.GetString("RampBMRM.checkpointMode");
                if(mode == "LATEST")
                        checkpointMode = KEEP_LATEST;
                if(mode == "ALL")
                        checkpointMode = KEEP_ALL;
        }   
}

#endif
