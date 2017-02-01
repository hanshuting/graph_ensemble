#ifndef _EPSILONINSENSITIVELOSS_HPP_
#define _EPSILONINSENSITIVELOSS_HPP_

#include "sml.hpp"
#include "univariateregressionloss.hpp"
#include "model.hpp"

/**  
 * Class to represent epsilon insensitive regression loss:
 * loss = max(0, |y - f(x)| - epsilon)
 * where f(x) := <w, x> and epsilon \geq 0. By default margin = 0.1.
 */
class CEpsilonInsensitiveLoss : public CUnivariateRegressionLoss 
{
protected:
	/** The epsilon parameter in the \epsilon-insensitve loss function
	 */
	Scalar epsilon;
	
	void Loss(Scalar& loss, TheMatrix& f);
	void LossAndGrad(Scalar& loss, TheMatrix& f, TheMatrix& l);
	
public:    
	CEpsilonInsensitiveLoss(CModel* &model, CVecData* &data);
	virtual ~CEpsilonInsensitiveLoss() {}
	virtual void Usage();
    
      
};

#endif
