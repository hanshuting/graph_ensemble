/* Copyright (c) 2006, National ICT Australia 
 * All rights reserved. 
 * 
 * The contents of this file are subject to the Mozilla Public License 
 * Version 1.1 (the "License"); you may not use this file except in 
 * compliance with the License. You may obtain a copy of the License at 
 * http://www.mozilla.org/MPL/ 
 * 
 * Software distributed under the License is distributed on an "AS IS" 
 * basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the 
 * License for the specific language governing rights and limitations 
 * under the License. 
 * 
 * Authors: Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *          
 *
 * Created: (29/11/2007) 
 *
 * Last Updated:
 */

#ifndef _SOLVERFACTORY_HPP_
#define _SOLVERFACTORY_HPP_

#include <string>

#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "model.hpp"
#include "loss.hpp"

// solver headers
//#include "bmrm.hpp"
//#include "sublbfgs.hpp"
#include "smd.hpp"

#include "lbfgs.hpp"

#include "olbfgs.hpp"


/**  
 * Factory class for creating new solvers 
 *
 * When you subclass the CSolver class you must also add an entry into
 * this factory class (grep YOUR_SOLVER for an example). 
 */
class CSolverFactory
{
   public:
      
      /**  Return an instance of solver based on user's argument in configuration file
       *
       *   @param model [read] Pointer to CModel object
       *   @return loss object
       */
      static CSolver* GetSolver(CModel* &model, CLoss* &loss)
      {         
         CSolver* solver = 0;
         Configuration &config = Configuration::GetInstance();
         
        
         
         // select the loss function specified by user (in configuration file)
         if(config.IsSet("Solver.type"))
         {
            string solverType = config.GetString("Solver.type");
            //if(solverType == "BMRM")
            //{ 
            //   solver = new CBMRM(model, loss);
	    //}
//             else if(solverType == "SUBLBFGS") 
//             {
//               solver = new CSUBLBFGS(model);
//             }
            else if(solverType == "SMD") 
            {
        	      solver = new CSMD(model, loss);
            }
            else if(solverType == "OLBFGS") 
            {
        	      solver = new COLBFGS(model, loss);
            }
            else if(solverType == "LBFGS") 
            {
        	      solver = new CLBFGS(model, loss);
            }
            //else if(solverType == "YOUR_SOLVER") 
            //{
            //  solver = new CYourSolver(model);
            //}
            else
            {
               throw CBMRMException("ERROR: unrecognised solver ("+solverType+")\n", "CSolverFactory::GetSolver()");
            }
         } 
         else
         {
            throw CBMRMException("ERROR: no solver specified!\n", "CSolverFactory::GetSolver()");
         }
         
         return solver;  
      }
      
};

#endif
