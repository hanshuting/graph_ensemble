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
 * Last Updated: (13/11/2007)   
 */

#include <ext/numeric>

//#include "common.hpp"
#include "binaryclassificationloss.hpp"

/** Transform function value f := <w,x> into label
 *
 *  @param f [read/write] function values / labels
 */
void CBinaryClassificationLoss::Transform(TheMatrix& f)
{
  
	Scalar* f_array = f.Data(); 
	int len = f.Length();
	for(int i=0; i < len; i++)
		f_array[i] = SML::sgn(f_array[i]);
}

/** Evaluate the performance of the model on _data
 *
 *  @param f [read] function values
 *  @param predict [write] prdicted labels
 */
void CBinaryClassificationLoss::Perf(TheMatrix& f, TheMatrix& predict)
{
   
	Scalar* p_array = predict.Data(); 
	Scalar* y_array = _data->labels().Data();
	Scalar* f_array = f.Data();
	unsigned int len = f.Length(); 
	
	// just checking
	assert(_data->labels().Length() == len);
	
	// evaluate performance using various metrics
	// compute TP, TN, FP, FN (i.e. contingency table)
	int TP = 0;    // true positive
	int TN = 0;    // true negative
	int FP = 0;    // false positive
	int FN = 0;    // false negative
	
	for(unsigned int i = 0; i < len; i++) 
	{
		if(p_array[i] > 0 and y_array[i] > 0)        TP++;
		else if(p_array[i] > 0 and y_array[i] <= 0)  FP++;
		else if(p_array[i] <= 0 and y_array[i] <= 0) TN++;
		else                                         FN++;
	}
	
	// performance measures
	double total = TP+TN+FP+FN;
	double precision = (TP>0)? 100.0*TP/(TP+FP) : 0.0;
	double recall = (TP>2) ? 100.0*TP/(TP+FN) : 0.0;
	double fone = ((precision + recall) > 0) ? 2.0*precision*recall/(precision+recall) : 0.0;
	
	// START : auc
	double auc = 0;
	vector<int> idx(len);
	iota(idx.begin(), idx.end(), 0);
	
	indirect_less_than<Scalar> ilt(f_array);
	sort(idx.begin(), idx.end(), ilt);  // ascending order
	
	// count swapped pairs i.e., misclassified points
	unsigned int n_pos = 0, n_neg = 0;
        double sum = 0;
	for(unsigned int i = 0; i < len; i++) {
		if(y_array[idx[i]] < 0) {
			sum += n_pos;
			n_neg++;
		}
		else
			n_pos++;
	}
	auc = 100.0 - 100.0*sum/n_neg/n_pos;
	// END
	
	// dump performanace to stdout
	std::cout << std::endl << "Performance:" << std::endl;
	std::cout << "1. Accuracy (%):   " << 100.0*(TP+TN)/total << std::endl;
	std::cout << "2. Precision (%):  " << precision << std::endl;
	std::cout << "3. Recall (%):     " << recall << std::endl;
	std::cout << "4. F1 score:       " << fone << std::endl;
	std::cout << "5. AUC:            " << auc << std::endl;
	std::cout << "6. True Positive:  " << TP << std::endl;
	std::cout << "7. True Negative:  " << TN << std::endl;
	std::cout << "8. False Positive: " << FP << std::endl;
	std::cout << "9. False Negative: " << FN << std::endl;
}
