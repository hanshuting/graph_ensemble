/* Copyright (c) 2006, National ICT Australia
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * Authors      : Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 * Created      : 28/11/2007
 * Last Updated :
 */

#ifndef _BMRMINNERSOLVER_CPP_
#define _BMRMINNERSOLVER_CPP_


#include "bmrminnersolver.hpp"
#include "configuration.hpp"


CBMRMInnerSolver::CBMRMInnerSolver(double lambda)
   : iter(0), 
     verbosity(0), 
     bmrmLambda(lambda), 
     dim(0), 
     numOfConstraint(0)
{
        // read configuration file
        Configuration &config = Configuration::GetInstance();
        
        //bmrmLambda = config.GetDouble("BMRM.lambda");
        assert(lambda > 0.0);
        
        if(config.IsSet("InnerSolver.verbosity"))  
                verbosity = config.GetInt("InnerSolver.verbosity");

}

#endif
