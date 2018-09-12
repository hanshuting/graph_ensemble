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
 * Created: (14/01/2008) 
 *
 * Last Updated:
 */

#ifndef _ROCSCORELOSS_CPP_
#define _ROCSCORELOSS_CPP_

#include "rocscoreloss.hpp"
#include "sml.hpp"


CROCScoreLoss::CROCScoreLoss(CModel* &model, CVecData* &data, bool additiveLabelLoss) 
        : CBinaryClassificationLoss(model, data),
          _additiveLabelLoss(additiveLabelLoss),
          m_pos(0),
          m_neg(0),
          orig_idx(0),
          idx(0)
{
        const TheMatrix &Y = _data->labels();
        Scalar y_i = 0.0;
        for(unsigned int i=0; i < _data->size(); i++)
        {
                Y.Get(i,y_i);
                if(y_i > 0.0) 
                        m_pos++;
                else
                        m_neg++;
        }
        assert(m == (m_pos + m_neg));
        
        orig_idx = new int[m];        
        for(int i=0; i<m; i++) 
                orig_idx[i] = i;
        idx = new int[m];        
        scalingFactor = 1.0/m_pos/m_neg;
        
        if(verbosity > 0)
        {
                std::cout << "In CROCScoreLoss" << std::endl;
                std::cout << "Num. of +ve labels: " << m_pos << std::endl;
                std::cout << "Num. of -ve labels: " << m_neg << std::endl;
        }
}
        
        
CROCScoreLoss::~CROCScoreLoss()
{
        if(orig_idx) delete [] orig_idx;
        if(idx) delete [] idx;        
}



/**  
 * Display (usage) message.
 */
void CROCScoreLoss::Usage()
{
        std::cout << "ROCScoreLoss: No special parameters!" << std::endl; 
}



/**  
 *  Compute ROCScore loss. CAUTION: f is passed by reference and is
 *  changed within this function. This is done for efficiency reasons,
 *  otherwise we would have had to create a new copy of f. 
 *   
 *  @param loss [write] loss value computed.
 *  @param f [read/write] prediction vector. 
 */
void CROCScoreLoss::Loss(Scalar& loss, TheMatrix& f)
{
        LossAndGrad(loss, f, *l);
}

/**  
 *  Compute loss and partial derivative of ROCScore loss w.r.t f
 *   
 *  @param loss [write] loss value computed.
 *  @param f [r/w] = X*w
 *  @param l [write] partial derivative of loss w.r.t. f
 */
void CROCScoreLoss::LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l)
{
        const TheMatrix &Y = _data->labels();
        Scalar *Y_array = Y.Data();
        Scalar *f_array = f.Data();
        
        if(_additiveLabelLoss)
        {
                for(int i=0; i<m; i++) 
                        f_array[i] -= 0.5*Y_array[i];
        }
        
        indirect_less_than<Scalar> ilt(f_array);
        memcpy(idx, orig_idx, sizeof(int)*m);
        sort(idx, idx+m, ilt);
        
        int l_neg = m_neg;
        int l_pos = 0;
        //Scalar l_i = 0.0;
        Scalar y_i = 0.0;
                
        for(int i=0; i<m; i++)
        {
                y_i = Y_array[idx[i]];
                if(y_i < 0.0)
                {
                        l.Set(idx[i], l_pos);
                        l_neg--;
                }
                else
                {
                        l.Set(idx[i], -l_neg);
                        l_pos++;
                }                
        }        
        l.Dot(f, loss);                           
}

#endif
