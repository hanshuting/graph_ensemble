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
 * Created: (08/01/2008) 
 *
 * Last Updated:
 */

#include "configuration.hpp"
#include "sml.hpp"
#include "model.hpp"
#include "modelfactory.hpp"
#include "rampbmrm.hpp"
#include "loss.hpp"
#include "ramplossfactory.hpp"
#include "data.hpp"
#include "datafactory.hpp"
#include <fstream>


/** BMRM for ramp losses
 *  note:
 *  1. this type of bmrm takes two loss functions, one for convex part and
 *     the other for linearization of concave part
 *  2. users must supply the linearization loss for the convex loss function
 *     they want to work with
 */

int main(int argc, char** argv)
{
        /*  /0. decide how the configuration file is written. [user know which loss_cav is for loss_vex, so no configuration arg needed for the loss_cav]
         *  1. CRampBMRM
         *  /2. CRampLossFactory and GetRampLoss()
         *  /3. modify multiclass such that it takes a parameter to switch off the additive labelloss term
         */
	
	// sanity check
	if(argc < 2) 
	{
		std::cout << "Usage: ./ramp-bmrm-train config.file" << std::endl;
		std::cout << "Check the configfiles directory for examples" << std::endl;
		std::cout << "ERROR: No configuration file given!" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// Note that the config object is also a central information 
	// storage object. Individual functions/classes can access common 
	// information they need by calling Configuration::GetInstance()
	Configuration &config = Configuration::GetInstance();
	config.ReadFromFile(argv[1]);
	
        CModel* model = 0;
	CData* data = 0;          // loss_vec and loss_cav use the same dataset
	CLoss* loss_vex = 0;      // convex loss
        CLoss* loss_cav = 0;      // linearization of concave loss
	CRampBMRM* rampbmrm = 0;  // loss_vex and loss_cav use the same model
	
	try {
		// in learning mode
		config.SetString("Program.mode", "LEARNING");
		
		// serial computation with centralised data
		config.SetString("Computation.mode", "SERIAL_CENTRALISED_DS");
		
		model = CModelFactory::GetModel();     
		if(config.IsSet("Model.hotStartModel")) 
			model->Initialize(config.GetString("Model.hotStartModel"));
		
		data = CDataFactory::GetData();     
		
		CRampLossFactory::GetRampLoss(model, data, loss_vex, loss_cav); // loss will initialize model if model is not hotstarted
		
		rampbmrm = new CRampBMRM(model, loss_vex, loss_cav);
		
		rampbmrm->Train();
		
		loss_vex->Evaluate(model);
		
		std::string modelFn = "model";
		if(config.IsSet("Model.modelFile"))
			modelFn = config.GetString("Model.modelFile");
		model->Save(modelFn);
		
		// cleaning up 
		delete rampbmrm;
		delete model;
		delete loss_vex;
                delete loss_cav;
		delete data;
		
	} 
	catch(CBMRMException e) 
        {
		std::cout << e.Report() << std::endl;
	}
	
	return EXIT_SUCCESS;	
}
