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
 * Authors: S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
 *          Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 *
 * Created: (08/01/2008) 
 *
 * Last Updated:
 */

#ifndef _RAMPLOSSFACTORY_HPP_
#define _RAMPLOSSFACTORY_HPP_

#include <string>

#include "common.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "data.hpp"
#include "vecdata.hpp"
#include "model.hpp"

// loss function headers
#include "loss.hpp"
//#include "softmarginloss.hpp"
//#include "squaredsoftmarginloss.hpp"
//#include "epsiloninsensitiveloss.hpp"
#include "wtamulticlassloss.hpp"
#include "rocscoreloss.hpp"
#include "ndcgrankloss.hpp"
#include "rampndcgrankloss.hpp"
#include "smmsegmentationloss.hpp"
//#include "preferencerankloss.hpp"
//#include "ordinalregressionloss.hpp"


/**  
 * Factory class for creating new loss instances. 
 *
 * When you subclass the CLoss class you must also add an entry into
 * this factory class (grep YOUR_LOSS_FUNCTION for an example). 
 * (dynamic) cast the generic data object handler to your specific data object type, if
 * you have a special subclass of CData.
 */
class CRampLossFactory
{
public:
      
        /**  Instantiate one convex loss function and another linearization of the concave loss 
         *   based on user's argument in configuration file
         *
         *   @param model [read] Pointer to model object
         *   @param data [read] Pointer to data object  
         *   @param loss_vex [write] Convex loss function
         *   @param loss_cav [write] Linearization of concave loss function corresponding to loss_vex     
         */
        static void GetRampLoss(CModel* &model, CData* &data, CLoss* &loss_vex, CLoss* &loss_cav)
        {                         
                Configuration &config = Configuration::GetInstance();
	 
                // select the loss function specified by user (in configuration file)
                if(config.IsSet("Loss.lossFunctionType"))
                {       
                        std::string lossFunctionType = config.GetString("Loss.lossFunctionType");
                        if(lossFunctionType == "WTA_MULTICLASS")
                        {
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss_vex = new CWTAMulticlassLoss(model, vecdata);
                                
                                // in the ramp losses, users must know which loss_cav to use for their specific loss_vex!                                
                                loss_cav = new CWTAMulticlassLoss(model, vecdata, false); // with additive label loss switched off
                        }
                        else if(lossFunctionType == "ROC_SCORE")
                        {
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss_vex = new CROCScoreLoss(model, vecdata);
                                
                                // in the ramp losses, users must know which loss_cav to use for their specific loss_vex!                                
                                loss_cav = new CROCScoreLoss(model, vecdata, false); // with additive label loss switched off
                        }
                        else if(lossFunctionType == "NDCG_RANK")
                        {
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss_vex = new CNDCGRankLoss(model, vecdata);
                                
                                // in the ramp losses, users must know which loss_cav to use for their specific loss_vex!                                
                                loss_cav = new CRampNDCGRankLoss(model, vecdata); 
                        }
                        else if(lossFunctionType == "SMM_SEGMENTATION")
                        {
                                CSeqData *seqdata = 0;
                                if(! (seqdata = dynamic_cast<CSeqData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CSeqData","CLossFactory::GetLoss()");
                                }
                                loss_vex = new CSMMSegmentationLoss(model, seqdata);
                                
                                // in the ramp losses, users must know which loss_cav to use for their specific loss_vex!                                
                                loss_cav = new CSMMSegmentationLoss(model, seqdata, false); 
                        }

                        //else if(lossFunctionType == "YOUR_LOSS_FUNCTION")
                        //{
                        //      loss = new CYourLoss(w, data);
                        //}
                        else
                        {
                                throw CBMRMException("ERROR: unrecognised loss function ("+lossFunctionType+")\n", 
                                                     "CLossFactory::GetLoss()");
                        }
                }
                else
                {
                        throw CBMRMException("ERROR: no loss function specified!\n", "CLossFactory::GetLoss()");
                }
	 
        }      
};

#endif
