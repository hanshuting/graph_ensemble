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
 * Created: (02/11/2007) 
 *
 * Last Updated: (07/11/2007)   
 */

#ifndef _LOSSFACTORY_HPP_
#define _LOSSFACTORY_HPP_

#include <string>

#include "common.hpp"
#include "configuration.hpp"
#include "bmrmexception.hpp"
#include "data.hpp"
#include "vecdata.hpp"
#include "model.hpp"

// loss function headers
#include "loss.hpp"
#include "softmarginloss.hpp"
#include "squaredsoftmarginloss.hpp"
#include "rocscoreloss.hpp"
#include "epsiloninsensitiveloss.hpp"
#include "wtamulticlassloss.hpp"
#include "expmulticlassloss.hpp"
#include "preferencerankloss.hpp"
#include "ordinalregressionloss.hpp"
#include "faceloss.hpp"
#include "multilabelloss.hpp"
#include "smmsegmentationloss.hpp"
#include "graphmatchloss.hpp"

/**  
 * Factory class for creating new Loss instances. 
 *
 * When you subclass the CLoss class you must also add an entry into
 * this factory class (grep YOUR_LOSS_FUNCTION for an example). 
 * (dynamic) cast the generic data object handler to your specific data object type, if
 * you have a special subclass of CData.
 */
class CLossFactory
{
public:
      
        /**  Return an instance of loss function based on user's argument in configuration file
         *
         *   @param data [read] Pointer to data object
         *   @return loss object
         */
        static CLoss* GetLoss(CModel* &model, CData* &data)
        {         
                CLoss* loss = 0;
                Configuration &config = Configuration::GetInstance();
	 
                // select the loss function specified by user (in configuration file)
                if(config.IsSet("Loss.lossFunctionType"))
                {
                        string lossFunctionType = config.GetString("Loss.lossFunctionType");
	    
                        if(lossFunctionType == "SOFT_MARGIN")
                        { 
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data)))
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CSoftMarginLoss(model, vecdata);
                        } 
                        else if(lossFunctionType == "SQUARED_SOFT_MARGIN")
                        { 
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CSquaredSoftMarginLoss(model, vecdata);
                        }
                        else if(lossFunctionType == "ROC_SCORE")
                        { 
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CROCScoreLoss(model, vecdata);
                        }
                        else if(lossFunctionType == "EPSILON_INSENSITIVE")
                        { 
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data)))
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CEpsilonInsensitiveLoss(model, vecdata);
                        }
                        else if(lossFunctionType == "WTA_MULTICLASS")
                        {
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CWTAMulticlassLoss(model, vecdata);
                        }
                        else if(lossFunctionType == "EXP_MULTICLASS")
                        {
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CExpMulticlassLoss(model, vecdata);
                        } 
                        else if(lossFunctionType == "PREFERENCE_RANK")
                        {
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CPreferenceRankLoss(model, vecdata);
                        } 
                        else if(lossFunctionType == "ORDINAL_REGRESSION")
                        {
                                CVecData *vecdata = 0;
                                if(! (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new COrdinalRegressionLoss(model, vecdata);
                        }
                        else if(lossFunctionType == "MULTI_LABEL_CLASSIFICATION")
                        {
                                CMultilabelVecData *mlvecdata = 0;
                                if(! (mlvecdata = dynamic_cast<CMultilabelVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CMultilabelVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CMultilabelLoss(model, mlvecdata);
                        }
                        else if(lossFunctionType == "SMM_SEGMENTATION")
                        {
                                CSeqData *seqdata = 0;
                                if(! (seqdata = dynamic_cast<CSeqData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CMultilabelVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CSMMSegmentationLoss(model, seqdata);
                        }
                        
                        else if(lossFunctionType == "FACE_LOSS")
                        {
                                CVecData *vecdata = 0;
                                if(not (vecdata = dynamic_cast<CVecData*>(data))) 
                                {
                                        throw CBMRMException("unable to cast data into CVecData",
                                                             "CLossFactory::GetLoss()");
                                }
                                loss = new CFaceLoss(model, vecdata);
                        } 
                        else if(lossFunctionType == "GRAPHMATCH") 
                        {
                          CGraphData *graphdata = 0;
                          if(not (graphdata = dynamic_cast<CGraphData*>(data))) {
                              throw CBMRMException("unable to cast data into CGraphData","CLossFactory::GetLoss()");
                          }
                          loss = new CGraphMatchLoss(model, graphdata);
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
                return loss;  
        }      
};

#endif
