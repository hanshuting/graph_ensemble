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
 * Last Updated: (02/11/2007)   
 */

#include "configuration.hpp"
#include "bmrm.hpp"
#include "sml.hpp"
#include "loss.hpp"
#include "lossfactory.hpp"
#include "data.hpp"
#include "datafactory.hpp"
#include "model.hpp"
#include "modelfactory.hpp"

int main(int argc, char** argv)
{
   
	// sanity check
	if(argc < 2) 
	{
		std::cout << "Usage: linear-bmrm-predict config.file" << std::endl;
		std::cout << "Check the configfiles directory for examples" << std::endl;
		std::cout << "ERROR: No configuration file given!" << std::endl;
		exit(EXIT_FAILURE);
	}

	CData* data = 0;
	CLoss* loss = 0;
	CModel* model = 0;
	
	try {
                // the very first thing to do!
                Configuration &config = Configuration::GetInstance();
                config.ReadFromFile(argv[1]);
	

		// serial computation with centralised data
		config.SetString("Computation.mode", "SERIAL_CENTRALISED_DS");
		
		std::string modelFilename = config.GetString("Model.modelFile");
		std::string programMode = config.GetString("Program.mode");
		
		data = CDataFactory::GetData();     
		model = CModelFactory::GetModel();
		model->Initialize(modelFilename, data->dim());
		loss = CLossFactory::GetLoss(model,data);
		
		if(programMode == "PREDICTION")
			loss->Predict(model);
		else if(programMode == "EVALUATION")
			loss->Evaluate(model);
		else
			throw CBMRMException("unknown program mode <" + programMode +">","main()");
		
		
		// cleaning up
		delete model; model = 0;
		delete loss; loss = 0;
		delete data; data = 0;
		
	}
	catch(CBMRMException e) {
		cout << e.Report() << endl;
	}
	
	return EXIT_SUCCESS;
}
