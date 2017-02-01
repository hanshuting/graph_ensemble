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
 * Authors:  Simon Guenter (simon.guenter@nicta.com.au)
 *
 * Created: (02/01/2008) 
 *
 * Last Updated: (02/01/2008)   
 */

#include "configuration.hpp"
#include "sml.hpp"
#include "model.hpp"
#include "modelfactory.hpp"
#include "solver.hpp"
//#include "solverfactory.hpp"
#include "sublbfgs.hpp"
#include "loss.hpp"
#include "softmarginloss.hpp"
#include "data.hpp"
#include "datafactory.hpp"
#include <fstream>


int main(int argc, char** argv)
{
    srand(12345678);
	
    // sanity check
    if(argc < 2) 
	{
	    std::cout << "Usage: ./sublbfgs-train config.file" << std::endl;
	    std::cout << "Check the configfiles directory for examples" << std::endl;
	    std::cout << "ERROR: No configuration file given!" << std::endl;
	    exit(EXIT_FAILURE);
	}
    
    // Note that the config object is also a central information 
    // storage object. Individual functions/classes can access common 
    // information they need by calling Configuration::GetInstance()
    Configuration &config = Configuration::GetInstance();
    config.ReadFromFile(argv[1]);
    
    CVecData* vecdata = 0;
    CData* data = 0;
    CLoss* loss = 0;
    CSUBLBFGS* solver = 0;
    CModel* model = 0;
   
    try {
	// in learning mode
	config.SetString("Program.mode", "LEARNING");
	
	// serial computation with centralised data
	config.SetString("Computation.mode", "SERIAL_CENTRALISED_DS");


	
	model = CModelFactory::GetModel();     
	if(config.IsSet("Model.hotStartModel")) 
	    model->Initialize(config.GetString("Model.hotStartModel"));
	
	data = CDataFactory::GetData();     
        
        if (dynamic_cast<CVecData *>(data))
            vecdata = (CVecData *) data;
        else
            throw CBMRMException("Data must be of type CVecData\n","CSUBLBFGS::CSUBLFBGS()");
	
        vecdata->CreateFeatureMatrixRowViews();
        vecdata->CreateLabelMatrixRowViews();
	
	loss = new CSoftMarginLoss(model, vecdata); 
	
	solver = new CSUBLBFGS(model,vecdata);
	double C = config.GetDouble("SUBLBFGS.C");
	solver->setC(C);
	solver->Train();
	
	loss->Evaluate(model);
	
	std::string modelFn = "model";
	if(config.IsSet("Model.modelFile"))
	    modelFn = config.GetString("Model.modelFile");
	model->Save(modelFn);
	
	// cleaning up 
	delete solver; solver = 0;
	delete model; model = 0;
	delete loss; loss = 0;
	delete data; data = 0;
	
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
