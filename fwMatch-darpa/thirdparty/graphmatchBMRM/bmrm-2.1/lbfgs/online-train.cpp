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
 * Authors: Simon Guenter
 *
 * Created: (04/01/2008) 
 */

#include "configuration.hpp"
#include "sml.hpp"
#include "model.hpp"
#include "modelfactory.hpp"
#include "solver.hpp"
//#include "solverfactory.hpp"
#include "loss.hpp"
#include "softmarginloss.hpp"
//#include "lossfactory.hpp"
#include "data.hpp"
//#include "datafactory.hpp"
#include "onlinevecdata.hpp"
#include <fstream>
#include "lbfgs.hpp"
#include "olbfgs.hpp"
#include "smd.hpp"



/** main for online training 
    data is restricted to be of type vecdata
*/
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
	
	// Note that the config object is also a central information 
	// storage object. Individual functions/classes can access common 
	// information they need by calling Configuration::GetInstance()
	Configuration &config = Configuration::GetInstance();
	config.ReadFromFile(argv[1]);
	
	CVecData* data = 0;
	CLoss* loss = 0;
	//CSoftMarginLoss* loss  = 0;
	CSolver* solver = 0;
	CModel* model = 0;
	try {
		// in learning mode
		config.SetString("Program.mode", "LEARNING");
		
		// serial computation with centralised data
		config.SetString("Computation.mode", "SERIAL_CENTRALISED_DS");
		
		model = CModelFactory::GetModel();     
		if(config.IsSet("Model.hotStartModel")) 
			model->Initialize(config.GetString("Model.hotStartModel"));
        
        // HACK
        int usenum =1; // number of calls per mini batch
        int batchsize=1; 

       
        if(config.IsSet("Data.batchsize")) {
            batchsize = config.GetInt("Data.batchsize");
        }

        
        // numbers of calls depend on the solver 
        if(config.IsSet("Solver.type"))
         {
             
             string solverType = config.GetString("Solver.type");
             if(solverType == "SMD") {
                 usenum = 1;
                 if(config.IsSet("SMD.smdmu"))           {
                     double smdmu = config.GetDouble("SMD.smdmu");
                     if(smdmu>0)
                         usenum=2;
                 }
                 
             }
             if(solverType == "OLBFGS") {
                 usenum = 2;
                 if(config.IsSet("OLBFGS.smdmu"))           {
                     double smdmu = config.GetDouble("OLBFGS.smdmu");
                     if(smdmu>0)
                         usenum=4;
                 }
             }
             
            
         }
        else
            throw CBMRMException("ERROR: no solver specified!\n", "");
      
        
        // data MUST be of type vecdata
	data =  new COnlineVecData(batchsize, usenum);
        
        Scalar lambda = 0.0; // regularizer used for mini batch
        Scalar origlambda = 0.0; // original regularizer
        Scalar multlambda = ((Scalar) data->slice_size()) / ((Scalar)data->size());
        //multlambda = 1.0;
        string solverType = config.GetString("Solver.type");
       
        // setting number of online iterations (depends on mini batchsize
        //  and number of passes through data set)
        if(solverType == "SMD") {
	    solver = new CSMD(model,loss);
            if(config.IsSet("SMD.maxNumOfIter"))           {
                int iter = config.GetInt("SMD.maxNumOfIter");
                iter= (int) ceil(((double) iter) / multlambda);
                config.SetInt("SMD.maxNumOfIter",iter);
            }
            if(config.IsSet("SMD.lambda"))           {
                lambda = config.GetDouble("SMD.lambda");
                origlambda = lambda;
                lambda*=multlambda;
                config.SetDouble("SMD.lambda",lambda);
            }
        }   
        if(solverType == "OLBFGS") {
             solver = new COLBFGS(model,loss);
            if(config.IsSet("OLBFGS.maxNumOfIter"))           {
                int iter = config.GetInt("OLBFGS.maxNumOfIter");
                iter= (int) ceil(((double) iter) / multlambda);
                config.SetInt("OLBFGS.maxNumOfIter",iter);
            }
            if(config.IsSet("OLBFGS.lambda"))           {
                lambda = config.GetDouble("OLBFGS.lambda");
                origlambda = lambda;
                lambda*=multlambda;
                config.SetDouble("OLBFGS.lambda",lambda);
            }
        }
        // HACK: communication with solver
        config.SetDouble("Online.minibatches",1.0 / multlambda);
        
		
	//loss = CLossFactory::GetLoss(model, data); // loss will initialize model if model is not hotstarted
	loss = new CSoftMarginLoss(model, data);

	//solver = CSolverFactory::GetSolver(model, loss);
		
	solver->Train();
		
        // do ADDITIONAL LBFGS (if requested)
        int iter = 0;
        if(config.IsSet("LBFGS.maxNumOfIter"))           {
           iter = config.GetInt("LBFGS.maxNumOfIter");
        }
        
        if(iter!=0) { // apply batch LBFGS 
            
            config.SetDouble("LBFGS.lambda",origlambda);
            delete solver;
            //delete loss;
            
            data = new CVecData();  // reload of data 
            //loss = CLossFactory::GetLoss(model, data); 
	    loss = new CSoftMarginLoss(model, data);
                
            solver = new CLBFGS(model,loss);
            solver->Train();
            loss->Evaluate(model);
        }
		
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
