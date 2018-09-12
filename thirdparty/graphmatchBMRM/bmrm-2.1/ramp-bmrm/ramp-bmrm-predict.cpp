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
 * Created: (09/01/2008) 
 *
 * Last Updated:
 */

#include "configuration.hpp"
#include "bmrm.hpp"
#include "sml.hpp"
#include "loss.hpp"
#include "ramplossfactory.hpp"
#include "data.hpp"
#include "datafactory.hpp"
#include "model.hpp"
#include "modelfactory.hpp"

int main(int argc, char** argv)
{
   
	// sanity check
	if(argc < 2) 
	{
		std::cout << "Usage: ramp-bmrm-predict config.file" << std::endl;
		std::cout << "Check the configfiles directory for examples" << std::endl;
		std::cout << "ERROR: No configuration file given!" << std::endl;
		exit(EXIT_FAILURE);
	}

	// the very first thing to do!
	Configuration &config = Configuration::GetInstance();
	config.ReadFromFile(argv[1]);
	
	CData* data = 0;
	CLoss* loss_vex = 0;
	CLoss* loss_cav = 0;
	CModel* model = 0;
	
	try {
		// serial computation with centralised data
		config.SetString("Computation.mode", "SERIAL_CENTRALISED_DS");
		
		std::string modelFilename = config.GetString("Model.modelFile");
		std::string programMode = config.GetString("Program.mode");

		data = CDataFactory::GetData();     		
		model = CModelFactory::GetModel();
		model->Initialize(modelFilename, data->dim());
		CRampLossFactory::GetRampLoss(model,data, loss_vex, loss_cav);
		
		if(programMode == "PREDICTION")
			loss_vex->Predict(model);
		else if(programMode == "EVALUATION")
			loss_vex->Evaluate(model);
		else
			throw CBMRMException("unknown program mode <" + programMode +">","main()");
		
		// compute ramp loss function value
		Scalar lossVal_vex = 0.0;
		Scalar lossVal_cav = 0.0;
		
		loss_vex->ComputeLoss(lossVal_vex);
		loss_cav->ComputeLoss(lossVal_cav);

		std::cout << "a) Convex loss function value: " << lossVal_vex << std::endl;
		std::cout << "b) Concave loss function linearization value: " << lossVal_cav << std::endl;
		std::cout << "c) Ramp loss function value (a-b) : " << lossVal_vex - lossVal_cav << std::endl;

		// cleaning up
		delete model;
		delete loss_vex;
		delete loss_cav;
		delete data;
		
	}
	catch(CBMRMException e) {
		cout << e.Report() << endl;
	}
	
	return EXIT_SUCCESS;
}
