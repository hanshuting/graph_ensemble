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
 *          Qinfeng Shi (shiqinfeng@gmail.com)
 *
 * Created: 30/01/2008
 *
 * Last Updated: 4/2/8
 */

#ifndef _SMMSEGMENTATIONLOSS_HPP_
#define _SMMSEGMENTATIONLOSS_HPP_

#include <string>
#include "sml.hpp"
#include "data.hpp"
#include "seqdata.hpp"
#include "loss.hpp"
#include "model.hpp"

/** Class for semi markov model for sequence segmentation
 *  
 *  Reference:
 *  1. Q. Shi, Y. Altun, A.J. Smola, and S.V.N. Vishwanathan,
 *     "Semi-Markov Models for Sequence Segmentation",
 *     in EMNLP, 2007.
 */
class CSMMSegmentationLoss : public CLoss 
{    
protected:
        /** Pointer to data
         */
        CSeqData* _data;
    
        /** whether to include label loss \Delta or not
         */
        bool _additiveLabelLoss;
        
        /** Number of examples in data
         */
        unsigned int m;                

        /** Minimum segment duration
         */
        unsigned int minDuration;

        /** Maximum segment duration
         */
        unsigned int maxDuration;

        /** Maximum sequence length
         */
        unsigned int maxSequenceLength;
        
        /** The margin value vector used in dynamic programming
         */
        vector<double> M;
        
        /** The label loss value vector used in dynamic programming
         */
        vector<double> L;
        
        /** The back pointers vector used in dynamic programming to retrieve the optimal path
         */
        vector<int> A;
        
        /** find best label WITH label loss (for training)
         */
        void find_best_label(const vector<unsigned int> &y, const CSeqFeature::seqfeature_struct &x, const TheMatrix &w, vector<unsigned int> &ybar, Scalar &marginloss, Scalar &labelloss);
        
        /** find best label WITHOUT label loss (for testing)
         */
        void find_best_label(const CSeqFeature::seqfeature_struct &x, const TheMatrix &w, vector<unsigned int> &ybar, Scalar &marginloss);
        
        /** compute precision and recall
	 */
	void PrecRec(const vector<unsigned int> &y, const vector<unsigned int> &ybar, double &prec, double& rec);

        /** "Partial" and decomposed label loss propsed by ref[1]
         */
        double PartialDelta(const unsigned int ybar_i, const vector<unsigned int> &y);
        
        /** Full label loss propsed by ref[1]
         */
        double Delta(const vector<unsigned int> &ybar, const vector<unsigned int> &y);
        
        /** Compute decision function values and predict labels
         *
         *  @param model [read] Model object containing a weight vector
         *  @param ybar [write] Predicted labels
         *  @param f_ybar [write] Decision function values
         */
        virtual void DecisionAndPrediction(CModel* model, int* ybar, Scalar* f_ybar);
    
public:
        CSMMSegmentationLoss() {}
        CSMMSegmentationLoss(CModel* &model, CSeqData* &data, bool additiveLabelLoss=true);
        virtual ~CSMMSegmentationLoss();
    
        // Methods
        virtual void Usage();
        virtual void ComputeLoss(Scalar& loss);
        virtual void ComputeLossAndGradient(Scalar& loss, TheMatrix& grad);
        virtual void Predict(CModel* model);
        virtual void Evaluate(CModel* model);

};

#endif
