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
 * Created: (30/01/2008) 
 *
 * Last Updated:
 */

#ifndef _SMMSEGMENTATIONLOSS_CPP_
#define _SMMSEGMENTATIONLOSS_CPP_

#include "common.hpp"  // WriteFile()
#include "smmsegmentationloss.hpp"
#include "configuration.hpp"
#include <sstream>
#include <vector>


/** Constructor
 *
 *  @param w    [r/w] external pointer to weight vector
 *  @param data [r/w] pointer to data container
 *  @param additiveLabelLoss [read] flag to determine whether to "add" label loss to loss function or not (default: true)
 */
CSMMSegmentationLoss::CSMMSegmentationLoss(CModel* &model, CSeqData* &data, bool additiveLabelLoss)
        : _data(data),
          _additiveLabelLoss(additiveLabelLoss),
          m(0),
          minDuration(0),
          maxDuration(0),
          maxSequenceLength(0)
{                
	// initialize model
	if(not model->IsInitialized())
	{
		model->Initialize(1, _data->dim(), _data->bias());
	}

	// keep a pointer to model
	_model = model;

        // get/set problem parameters
        Configuration &config = Configuration::GetInstance();
        
        if(config.IsSet("SMM_SEGMENTATION.minDuration"))
                minDuration = config.GetInt("SMM_SEGMENTATION.minDuration");
        else
                minDuration = _data->MinDuration();
        
        if(config.IsSet("SMM_SEGMENTATION.maxDuration"))
                maxDuration = config.GetInt("SMM_SEGMENTATION.maxDuration");
        else
                maxDuration = _data->MaxDuration();                
        
        m = _data->slice_size();        
        scalingFactor = 1;
        maxSequenceLength = _data->MaxSequenceLength(); 
        
        M.resize(maxSequenceLength+1,0);
        L.resize(maxSequenceLength+1,0);
        A.resize(maxSequenceLength+1,0);
        
//        if(additiveLabelLoss)
//        {
//                FullDelta = &CSMMSegmentationLoss::FullSymmSetDiff;
//                PartialDelta = &CSMMSegmentationLoss::PartialSymmSetDiff;
//        }
//        else
//        {
//                FullDelta = &CSMMSegmentationLoss::FullZero;
//                PartialDelta = &CSMMSegmentationLoss::PartialZero;
//
//        }
		
	if(verbosity > 0)
	{
                std::cout << "In CSMMSegmentationLoss::CSMMSegmentationLoss()" << std::endl;
                if(!_additiveLabelLoss)
                        std::cout << "  with modification for ramp loss" << std::endl;
	}
        
        assert(m > 0);
        assert(minDuration <= maxDuration);
        printf("min duration: %d\n", minDuration);
        printf("max duration: %d\n", maxDuration);
}

CSMMSegmentationLoss::~CSMMSegmentationLoss()
{}

/** 
 * Display (usage) message. 
 */
void CSMMSegmentationLoss::Usage()
{
	std::cout << "Semi-Markov Model for Sequence Segmentation loss function:" << std::endl;
        std::cout << "1. Min segment duration (int SMM_SEGMENTATION.minDuration XX)" << std::endl;
        std::cout << "2. Max segment duration (int SMM_SEGMENTATION.maxDuration XX)" << std::endl;
}


/**  Compute loss
 */
void CSMMSegmentationLoss::ComputeLoss(Scalar& loss)
{
	    TheMatrix &w = _model->GetW();
        loss = 0;
        
        const vector<CSeqLabel::seqlabel_struct> &Y = _data->labels();
        const vector<CSeqFeature::seqfeature_struct> &X = _data->features();

        
        for(unsigned int i=0; i < m; i++)
        {
                vector<unsigned int> ybar;
                double labelloss = 0;
                double marginloss = 0;
                
                // find best label y' and return the score wrt to y'
                find_best_label(Y[i].pos, X[i], w, ybar, marginloss, labelloss);
                
                // accumulate loss
                loss += marginloss + labelloss;
        }
        
        loss *= scalingFactor;  
}


/**   Compute loss and gradient
 */
void CSMMSegmentationLoss::ComputeLossAndGradient(Scalar& loss, TheMatrix& grad)
{
	TheMatrix &w = _model->GetW();
        loss = 0;
        grad.Zero();
        TheMatrix g(grad, SML::DENSE);

        const vector<CSeqLabel::seqlabel_struct> &Y = _data->labels();
        const vector<CSeqFeature::seqfeature_struct> &X = _data->features();

        
        for(unsigned int i=0; i < m; i++)
        //unsigned int i=0;
        {
                vector<unsigned int> ybar(X[i].len,0);
                double labelloss = 0;
                double marginloss = 0;
                Scalar w_dot_g = 0.0;;

                // find best label y' and return the score wrt to y'
                if(_additiveLabelLoss)
                        find_best_label(Y[i].pos, X[i], w, ybar, marginloss, labelloss);
                else
                        find_best_label(X[i], w, ybar, marginloss);
                
                // construct the gradient vector for the part of true y
                const vector<unsigned int> &y = Y[i].pos;
		g.Zero();                
                for(unsigned int j=1; j < y.size(); j++)
                {
                        //std::cout << "y   j=" << j << std::endl; 
                        
                        g.Add(*(X[i].phi_1[y[j]]));
                        g.Add(*(X[i].phi_2[y[j-1]][y[j]-y[j-1]-1])); ////
                }
                if(y.size() > 0)
                {
                        g.Add(*(X[i].phi_2[y[y.size()-1]][X[i].len-1 - y[y.size()-1]-1]));////
                }

                // for predicted y'
                for(unsigned int j=1; j < ybar.size(); j++)
                {                        
                        grad.Add(*(X[i].phi_1[ybar[j]]));                         
                        grad.Add(*(X[i].phi_2[ybar[j-1]][ybar[j]-ybar[j-1]-1])); ////                          
                }
                if(ybar.size() > 0)
                {
                        grad.Add(*(X[i].phi_2[ybar[ybar.size()-1]][X[i].len-1 - ybar[ybar.size()-1]-1]));////                 
                }
                grad.Minus(g);

      		// accumulate the loss
		w.Dot(g, w_dot_g);
                loss = loss - w_dot_g + marginloss + labelloss;                
        }

	grad.Scale(scalingFactor);	
	loss *= scalingFactor;        
        
        if(verbosity)
        {
                Scalar gnorm = 0.0;
                grad.Norm2(gnorm);
                std::cout << "gradient norm=" << gnorm << std::endl;
        }
        //Evaluate(_model);
}


/** find best label (with label loss): g(w) := max_y' <w,\phi(x,y')> + Delta(y', y)
 *
 *  @param x [read] sequence
 *  @param y [read] actual label for x
 *  @param w [read] weight vector
 *  @param ybar [write] found best label
 *  @param marginloss [write] margin loss <w,\Phi(x,y')> w.r.t to best y'
 *  @param labelloss [write] label loss \Delta(y',y) w.r.t. to best y'
 *
 */
void CSMMSegmentationLoss::find_best_label(const vector<unsigned int> &y, const CSeqFeature::seqfeature_struct &x, const TheMatrix &w, vector<unsigned int> &ybar, Scalar &marginloss, Scalar &labelloss)
{
        using namespace std;
        
        // reset return values
        marginloss = 0;
        labelloss = 0;
        ybar.clear();
        
        // reset temporaries
        fill(M.begin(), M.end(), 0);
        fill(L.begin(), L.end(), 0);
        fill(A.begin(), A.end(),-1);
		        
	double maxval = -SML::INFTY;
        Scalar w_dot_phi1 = 0;
        Scalar w_dot_phi2 = 0;
        double marginval = 0;
        double labelval = 0;
        unsigned int right = 0;
        unsigned int left = 0;
        unsigned int start = 0;
        unsigned int end = 0;
        double sum = 0;
        
	// compute DP statistics for positions 1 to len-1
        L[0] += y.size()-2;
        A[1] = 0;

        for(right=1; right < x.len; right++)
	{
	        // \Phi = (phi1, phi2[left,right])
                // <w, \Phi> = <w,phi1> + <w,phi[left,right]>                
                maxval = -SML::INFTY;
                w_dot_phi1 = 0.0;
                w.Dot(*(x.phi_1[right]), w_dot_phi1);

                start = max(0,int(right-maxDuration));
                end = right;//-minDuration+1;
                for(left=start; left < end; left++)
                {
                        w.Dot(*(x.phi_2[left][right-left-1]), w_dot_phi2); 
                        marginval = w_dot_phi1 + w_dot_phi2;                       
                        labelval = PartialDelta(right,y);                        
                        sum = M[left]+marginval + L[left]+labelval;
                        if(sum > maxval)
                        {
                                A[right] = left;
                                M[right] = M[left] + marginval;
                                L[right] = L[left] + labelval;
                                maxval = sum;
                        }
                }	  
        }
        
        // get optimal path (i.e. segmentation)        
        int i = x.len-1;
        while(A[i] >= 0)
        {
                ybar.push_back(A[i]);
                i = A[i];
        }
        
        marginloss = M[x.len-1];
        labelloss = L[x.len-1];
        reverse(ybar.begin(), ybar.end());
}


/** find best label (without label loss): g(w) := max_y' <w,\phi(x,y')>
 *
 *  @param x [read] sequence
 *  @param w [read] weight vector
 *  @param ybar [write] found best label
 *  @param marginloss [write] margin loss <w,\Phi(x,y')> w.r.t to best y'
 */
void CSMMSegmentationLoss::find_best_label(const CSeqFeature::seqfeature_struct &x, const TheMatrix &w, vector<unsigned int> &ybar, Scalar &marginloss)
{
        using namespace std;
        
        // reset return values
        marginloss = 0;
        ybar.clear();
        
        // reset temporaries
        fill(M.begin(), M.end(), 0);
        fill(A.begin(), A.end(),-1);
		        
	double maxval = -SML::INFTY;
        Scalar w_dot_phi1 = 0;
        Scalar w_dot_phi2 = 0;
        double marginval = 0;
        unsigned int right = 0;
        unsigned int left = 0;
        unsigned int start = 0;
        unsigned int end = 0;
        double sum = 0;
        
	// compute DP statistics for positions 1 to len-1
        A[1] = 0;

        for(right=1; right < x.len; right++)
	{
	        // \Phi = (phi1, phi2[left,right])
                // <w, \Phi> = <w,phi1> + <w,phi[left,right]>                
                maxval = -SML::INFTY;
                w_dot_phi1 = 0.0;
                w.Dot(*(x.phi_1[right]), w_dot_phi1);

                start = max(0,int(right-maxDuration));
                end = right;//-minDuration+1;
                for(left=start; left < end; left++)
                {
                        w.Dot(*(x.phi_2[left][right-left-1]), w_dot_phi2); 
                        marginval = w_dot_phi1 + w_dot_phi2;                       
                        sum = M[left]+marginval;
                        if(sum > maxval)
                        {
                                A[right] = left;
                                M[right] = M[left] + marginval;
                                maxval = sum;
                        }
                }	  
        }
        
        // get optimal path (i.e. segmentation)        
        int i = x.len-1;
        while(A[i] >= 0)
        {
                ybar.push_back(A[i]);
                i = A[i];
        }
        
        marginloss = M[x.len-1];
        reverse(ybar.begin(), ybar.end());
}



/** Partial label loss propsed by ref[1]
 *
 *  @param ybar_i [read] the newly selected starting position of a segment
 *  @param y [read] the full label i.e., all segment starting positions of a sequence
 *
 *  @return partial label loss value
 */
double CSMMSegmentationLoss::PartialDelta(const unsigned int ybar_i, const vector<unsigned int> &y)
{
        if(binary_search(y.begin(), y.end(), ybar_i))
                return -1;
        else
                return 1;
}


double CSMMSegmentationLoss::Delta(const vector<unsigned int> &ybar, const vector<unsigned int> &y)
{
        unsigned int same = 0;
        unsigned int i=0, j=0;
        
	while(i<y.size() && j<ybar.size())
	{
		if(y[i] < ybar[j] && i<y.size()) 
			i++;
		else if(y[i] > ybar[j] && j<ybar.size()) 
			j++;
		else if(y[i] == ybar[j])
		{
			same += 1;
			i++;
			j++;
		}
	}
        
        return (y.size()+ybar.size() - 2*same);
}

        


void CSMMSegmentationLoss::Predict(CModel* model)
{
//        // sanity check
//        assert(model);
//        assert(_data);
//        assert(m > 0);
//               
//        // more defensive precaution
//        assert(_data->dim() == model->GetDimOfW());
//        
//        if(numOfClass != model->GetNumOfW())
//        {
//                if(matW) 
//                {
//                        delete matW;                        
//                        matW = new TheMatrix(model->GetNumOfW(), model->GetDimOfW(), SML::DENSE);
//                }                        
//        }
//        
//        numOfClass = std::min(numOfClass, model->GetNumOfW());                      
//        assert(matW);
//        
//        int* ybar = new int[m];
//	Scalar* f_ybar = new Scalar[m];
//        
//        DecisionAndPrediction(model, ybar, f_ybar);
// 
//        Configuration &config = Configuration::GetInstance();
//        if(config.GetString("Program.mode") != "LEARNING")
//	{
//                bool success = false;
//                std::string predictedLabelsFn = "predictedLabels";
//                if(config.IsSet("Prediction.predictedLabelsFile"))
//                        predictedLabelsFn = config.GetString("Prediction.predictedLabelsFile");
//                
//                success = WriteFile<int*>(predictedLabelsFn, m, ybar);
//                if(success) std::cout << "Predicted labels file written." << std::endl;
//                
//                std::string decisionFunctionValuesFn = "decisionFunctionValuesFile";
//                if(config.IsSet("Prediction.decisionFunctionValuesFile"))
//                        decisionFunctionValuesFn = config.GetString("Prediction.decisionFunctionValuesFile");
//                
//                success = WriteFile<Scalar*>(decisionFunctionValuesFn, m, f_ybar);
//                if(success) std::cout << "Decision function values file written." << std::endl;
//        }
//        
//	// clean up
//	if(ybar) delete [] ybar;
//	if(f_ybar) delete [] f_ybar;
}
 
/**  compute precision and recall
 *
 *   @param y [read] actual label
 *   @param ybar [read] predicted label
 *   @param prec [write] precision
 *   @param rec [write] recall
 */
void CSMMSegmentationLoss::PrecRec(const vector<unsigned int> &y, const vector<unsigned int> &ybar, double &prec, double &rec)
{
        double same = 0;
	prec = 0.0;
	rec = 0.0;
	

	//F_alpha = (1+alpha)*prec*rec/(alpha*prec+rec)
	unsigned int i=0, j=0;
	while(i<y.size() && j<ybar.size())
	{
		if(y[i] < ybar[j] && i<y.size()) 
			i++;
		else if(y[i] > ybar[j] && j<ybar.size()) 
			j++;
		else if(y[i] == ybar[j])
		{
			same += 1;
			i++;
			j++;
		}
	}
	
	prec = same/ybar.size();
	rec = same/y.size();
}

/**  Evaluate the performance of the model on the examples in _data
 *
 *   @param model [read] model object containing a weight vector
 */

void CSMMSegmentationLoss::Evaluate(CModel* model)
{
       // sanity check
       assert(model);
       assert(_data);
       assert(m > 0);
       
       if(not _data->HasLabel())
       {
               throw CBMRMException("Data object does not contain labels","CSMMSegmentationLoss::Evaluate()");
       }
       
       // more defensive precaution
       assert(_data->HasLabel());
       assert(_data->dim() == model->GetDimOfW());
       
	    TheMatrix &w = _model->GetW();
        const vector<CSeqLabel::seqlabel_struct> &Y = _data->labels();
        const vector<CSeqFeature::seqfeature_struct> &X = _data->features();

        double unweighted_avgf1 = 0;
        double weighted_prec = 0;
        double weighted_rec = 0;
        double weighted_avgf1 = 0;
        double totalsegment = 0;
        
        for(unsigned int i=0; i < m; i++)
        {
                vector<unsigned int> ybar;
                double marginloss = 0;
                
                // find best label y' and return the score wrt to y'
                find_best_label(X[i], w, ybar, marginloss);
                
                // accumulate loss                
                //double loss = marginloss + labelloss;                
		double prec = 0;
                double rec = 0;
                double f1 = 0;
                
                PrecRec(Y[i].pos,ybar,prec,rec);                			
                f1 = (2*prec*rec)/(prec+rec);
                unweighted_avgf1 += f1;
                weighted_prec += Y[i].pos.size()*prec;
                weighted_rec += Y[i].pos.size()*rec;
                totalsegment += Y[i].pos.size();
                
		printf("#%2d   F1: %2.4f\n",i,f1);
        }
        weighted_prec /= totalsegment;
        weighted_rec /= totalsegment;
        weighted_avgf1 = 2*weighted_prec*weighted_rec/(weighted_prec+weighted_rec);
        unweighted_avgf1 /= m;
                
        std::cout << "1. unweighted macro-average F1: " << unweighted_avgf1 << std::endl;
        std::cout << "2. weighted macro-average F1  : " << weighted_avgf1 << std::endl;
        
		// micro or macro averaging f-measure
	   
}


void CSMMSegmentationLoss::DecisionAndPrediction(CModel* model, int *ybar, Scalar *f_ybar)
{
//        // sanity check	
//        assert(ybar);
//        assert(f_ybar);
//        
//        TheMatrix &w = model->GetW();        
//        matW->Assign(w);        	
//	        
//	for(unsigned int i=0; i < m; i++)
//	{
//		_data->ComputeFi(*matW, *f, i);
//		Scalar max_f_ybar_i = -SML::INFTY;
//		
//		for(unsigned int j = 0; j < numOfClass; j++)
//		{
//			f->Get(j, f_ybar[i]);
//			if(f_ybar[i] > max_f_ybar_i)
//			{
//				max_f_ybar_i = f_ybar[i];
//				ybar[i] = j+1;
//			}
//		} 		
//	}                        
}

#endif
