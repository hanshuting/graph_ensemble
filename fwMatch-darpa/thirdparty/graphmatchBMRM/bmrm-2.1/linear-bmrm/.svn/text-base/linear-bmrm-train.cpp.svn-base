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
 * Created: (02/11/2007) 
 *
 * Last Updated: (10/01/2008)   
 *
 * Note: [chteo:100108:1743] now with support for model selection over lambda parameters.
 *       The model files will be appended with "_lambda_XXX" where XXX is the lambda used 
 *       to produce that model.
 */

#include "configuration.hpp"
#include "sml.hpp"
#include "model.hpp"
#include "modelfactory.hpp"
#include "solver.hpp"
#include "bmrm.hpp"
#include "loss.hpp"
#include "lossfactory.hpp"
#include "data.hpp"
#include "datafactory.hpp"
#include <fstream>

int main(int argc, char** argv)
{	
	// sanity check
	if(argc < 2) 
	{
		std::cout << "Usage: ./linear-bmrm-train config.file" << std::endl;
		std::cout << "Check the configfiles directory for examples" << std::endl;
		std::cout << "ERROR: No configuration file given!" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	CData* data = 0;
	CLoss* loss = 0;
	CSolver* solver = 0;
	CModel* model = 0;
 	               
	try {
                // Note that the config object is also a central information 
                // storage object. Individual functions/classes can access common 
                // information they need by calling Configuration::GetInstance()
                Configuration &config = Configuration::GetInstance();
                config.ReadFromFile(argv[1]);
                                
		// in learning mode
		config.SetString("Program.mode", "LEARNING");
		
		// serial computation with centralised data
		config.SetString("Computation.mode", "SERIAL_CENTRALISED_DS");
		
                model = CModelFactory::GetModel();     
		if(config.IsSet("Model.hotStartModel")) 
			model->Initialize(config.GetString("Model.hotStartModel"));		
		data = CDataFactory::GetData();		
		loss = CLossFactory::GetLoss(model, data); // loss will initialize model if model is not hotstarted                                             
                                
                vector<double> lambdas;
                if(config.IsSet("BMRM.lambdas"))
                {
                        lambdas = config.GetDoubleVector("BMRM.lambdas");
                        std::cout << "Main(): Multiple lambda parameters detected! Training with each of them..." << std::endl;
                        
                        for(size_t i=0; i < lambdas.size(); i++)
                        {
                                config.SetDouble("BMRM.lambda", lambdas[i]);
                                std::cout << "\n[Learning using lambda: " << lambdas[i] << "]" << std::endl;
                                                                   		
                                if(solver) delete solver;
                                //solver = CSolverFactory::GetSolver(model, loss);		
				solver = new CBMRM(model, loss);
                                solver->Train();		
                                loss->Evaluate(model);
		
                                ostringstream oss;
                                oss << lambdas[i];
                                std::string lambda_str = oss.str();
                                std::string modelFn = "model_lambda_" + lambda_str;
                                if(config.IsSet("Model.modelFile"))
                                        modelFn = config.GetString("Model.modelFile");
                                modelFn = modelFn + "_lambda_" + lambda_str;
                                model->Save(modelFn);
                        }
                }
                else
                {
		        //solver = CSolverFactory::GetSolver(model, loss);		
		        solver = new CBMRM(model, loss);		
                        solver->Train();		
                        loss->Evaluate(model);
		
                        std::string modelFn = "model";
                        if(config.IsSet("Model.modelFile"))
                                modelFn = config.GetString("Model.modelFile");
                        model->Save(modelFn);
                }                
                
		// cleaning up 
		if(solver) delete solver;
		if(model) delete model;
		if(loss) delete loss;
		if(data) delete data;		
	} 
	catch(CBMRMException e) {
		std::cerr << e.Report() << std::endl;
	}
	
	return EXIT_SUCCESS;	
}


/*
  Note:
  1.  For kernel case, SML::TheMatrix will be a structure which is no longer a matrix class.
      I.e., it could be a structure with the same method interfaces as in the original 
      TheMatrix class but with different function definitions.
      
 */
