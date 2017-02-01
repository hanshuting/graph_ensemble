#ifndef _BMRM_CPP_
#define _BMRM_CPP_

#include "common.hpp"
#include "bmrm.hpp"
#include "timer.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "bmrminnersolver.hpp"
#include "loss.hpp"
#include "bmrminnersolverfactory.hpp"

#include <fstream>
#include <sstream>


/**  
 *  Constructor
 *
 *  @param model [read] pointer to loss model object
 */
CBMRM::CBMRM(CModel *model, CLoss* loss)
   : CSolver(model, loss)
{
   // set private members (default) values
   verbosity        = 0;
   convergenceLog   = 0;
   maxNumOfIter     = 10000;
   epsilonTol       = 1e-4;
   gammaTol         = 0.0;;
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
CBMRM::~CBMRM()
{
   // destroy inner solver
   if(innerSolver) { delete innerSolver; innerSolver = 0; }
}



/**  Start training/learning a model w.r.t. the loss object (and the data supplied to it).
 *
 */
void CBMRM::Train()
{

   // Timers (CPU and wall-clock)
   CTimer totalTime;             // total runtime of the training
   CTimer innerSolverTime;       // time for inner optimization (e.g., QP or LP)
   CTimer lossAndGradientTime;   // time for loss and gradient computation
   
   unsigned int iter = 0;              // iteration count
   Scalar loss = 0.0;            // loss function value        
   Scalar xi = 0.0;              // convex lower-bound (approximate) of loss value
   Scalar exactObjVal = 0.0;     // (exact) objective function value
   Scalar approxObjVal = 0.0;    // convex lower-bound (approximate) of objective function value
   Scalar minExactObjVal = 1e99; // minimum of all previously evaluated (exact) objective function value
   Scalar regVal = 0.0;          // value of the regularizer term e.g., 0.5*w'*w
   Scalar epsilon = 0.0;         // (duality) gap := exactObjVal - approxObjVal
   Scalar gamma = 0.0;           // := minExactObjVal - approxObjVal
   //Scalar w_dot_a = 0.0;         // temp for <w,a> where 'w' is weight vector and 'a' is gradient
   double innerSolverTol = 1.0;  // optimization tolerance for inner solver

   ofstream lossFp;              // keep loss values
   ofstream xiFp;                // keep keep xi values
   ofstream exactObjValFp;       // keep exactObjVal values
   ofstream approxObjValFp;      // keep approxObjVal values
   ofstream regValFp;            // keep regVal values
   ofstream epsilonFp;           // keep epsilon values
   ofstream gammaFp;             // keep gamma values

   unsigned int row = 0; 
   unsigned int col = 0;
   TheMatrix &w = _model->GetW();
   w.Shape(row, col);   
   TheMatrix a(row, col, SML::DENSE);   // gradient vector
  
   
   // convergence log files
   if(convergenceLog) 
   {
      lossFp.open("loss.dump");    
      xiFp.open("xi.dump");
      exactObjValFp.open("exactobj.dump");
      approxObjValFp.open("approxobj.dump");
      regValFp.open("regval.dump");
      epsilonFp.open("epsilon.dump");
      gammaFp.open("gamma.dump");
   }

   // start training
   totalTime.Start();


   while(1)
   {
      iter++;
      
      // column generation
      lossAndGradientTime.Start();
      _loss->ComputeLossAndGradient(loss, a);
      lossAndGradientTime.Stop();
      
      // update convergence monitor
      exactObjVal = loss + regVal;
      minExactObjVal = std::min(minExactObjVal, exactObjVal);
      epsilon = exactObjVal - approxObjVal;
      gamma   = minExactObjVal - approxObjVal;

      // dump convergence statistics into files
      if(convergenceLog) {
         lossFp         << loss         << endl;
         xiFp           << xi           << endl;
         exactObjValFp  << exactObjVal  << endl;   
         approxObjValFp << approxObjVal << endl;
         regValFp       << regVal       << endl;    
         epsilonFp      << epsilon      << endl;
         gammaFp        << gamma        << endl;
      }
      
      // dump convergence statistics on stdout
      if(verbosity < 1) {
         printf(".");
         if(iter%100 == 0) 
            printf("%d",iter);
         fflush(stdout);
      }
      else if(verbosity == 1)
         printf("#%d   eps %.6e   loss %.6e   xi %.6e   reg %.6e\n",iter, epsilon, loss, xi, regVal);      
      else if(verbosity > 1)
         printf("#%d   pobj %.6e   aobj %.6e   eps %.6e   gam %.6e   loss %.6e   xi % .6e   reg %.6e\n",
		iter, exactObjVal, approxObjVal, epsilon, gamma, loss, xi, regVal);      


      // check point
      if(iter%checkpointInterval == 0) 
      {
         if(checkpointMode == KEEP_LATEST)
            _model->Save(checkpointPrefix);
         else {
            ostringstream oss;
            oss << checkpointPrefix << "." << iter;
            _model->Save(oss.str());
         }
      }
        
      // stopping criteria
      if((iter >= 2) && ((gamma < gammaTol) || (epsilon < epsilonTol)))
         break;
      if(iter >= maxNumOfIter)
      { 
         printf("\nWARNING: program exceeded maximum number of iterations (%d) !\n", maxNumOfIter);
         break;
      }

      
      // adjust inner solver optimization tolerance
      innerSolverTol = std::min(innerSolverTol, std::max((double)epsilon, (double)epsilonTol));
      innerSolver->SetTolerance(innerSolverTol*0.5);        

      // run inner solver
      innerSolverTime.Start();
      innerSolver->Solve(w, a, loss, xi, regVal, approxObjVal);
      innerSolverTime.Stop();
     
   }

   // legends
   if(verbosity >= 1) {
      printf("\nLegends::\n");
      if(verbosity > 1)
	 printf("pobj: primal objective function value\naobj: approximate objective function value\n");
      printf("eps: epsilon (approximation error) \ngam: lower bound on eps \nloss: loss function value \nxi: approximation to loss \nreg: regularizer value\n");
   }
      

   Scalar norm1 = 0, norm2 = 0, norminf = 0;
   w.Norm1(norm1);
   w.Norm2(norm2);
   w.NormInf(norminf);
   printf("\n");
   printf("No. of iterations:  %d\n",iter);
   printf("Primal obj. val.: %.6e\n",exactObjVal);
   printf("Approx obj. val.: %.6e\n",approxObjVal);
   printf("Primal - Approx.: %.6e\n",exactObjVal-approxObjVal);
   printf("Loss:             %.6e\n",loss);
   printf("|w|_1:            %.6e\n",norm1);
   printf("|w|_2:            %.6e\n",norm2);
   printf("|w|_oo:           %.6e\n",norminf);
   
   totalTime.Stop();
   // end of training

   // display timing profile
   printf("\nCPU seconds in:\n");
   printf("1. loss and gradient: %8.2f\n", lossAndGradientTime.CPUTotal());
   printf("2. solver:            %8.2f\n",innerSolverTime.CPUTotal()); 
   printf("               Total: %8.2f\n", totalTime.CPUTotal());

   // clean up
   if(convergenceLog)
   {
      lossFp.close();
      xiFp.close();
      exactObjValFp.close();
      approxObjValFp.close();
      regValFp.close();
      epsilonFp.close();
      gammaFp.close();
   }   
}


/**   Validate program parameters set in Configuration.
 */
void CBMRM::ConfirmProgramParameters()
{
   Configuration &config = Configuration::GetInstance();  // make sure configuration file is read before this!

   if(config.IsSet("BMRM.verbosity")) {
      verbosity = config.GetInt("BMRM.verbosity");
   }
   
   if(config.IsSet("BMRM.convergenceLog")) {
      convergenceLog = config.GetInt("BMRM.convergenceLog");
   }
   
   if(config.IsSet("BMRM.maxNumOfIter")) {
      maxNumOfIter = config.GetInt("BMRM.maxNumOfIter");
      if(maxNumOfIter < 0)
	 throw CBMRMException("BMRM.maxNumOfIter must be > 0\n","CBMRM::ConfirmProgramParameters()");
   }
   
   if(config.IsSet("BMRM.epsilonTol")) {
      epsilonTol = config.GetDouble("BMRM.epsilonTol");
   }

   if(config.IsSet("BMRM.gammaTol")) {
      gammaTol = config.GetDouble("BMRM.gammaTol");
   }
   
   if(config.IsSet("BMRM.lambda"))           {
      lambda = config.GetDouble("BMRM.lambda");
      if(lambda <= 0)
         throw CBMRMException("BMRM.lambda must be > 0\n","CBMRM::ConfirmProgramParameters()");
   }
   else {
      config.SetDouble("BMRM.lambda", lambda);
   }


   if(config.IsSet("BMRM.checkpointInterval")) {
      checkpointInterval = config.GetInt("BMRM.checkpointInterval");
      if(checkpointInterval < 1)
         throw CBMRMException("BMRM.checkpointInterval must be a positive integer!\n","CBMRM::ConfirmProgramParameters()");
   }

   if(config.IsSet("BMRM.checkpointPrefix")) {
      checkpointPrefix = config.GetString("BMRM.checkpointPrefix");
   }

   if(config.IsSet("BMRM.checkpointMode")) {
      string mode = config.GetString("BMRM.checkpointMode");
      if(mode == "LATEST")
         checkpointMode = KEEP_LATEST;
      if(mode == "ALL")
         checkpointMode = KEEP_ALL;
   }   
}

#endif
