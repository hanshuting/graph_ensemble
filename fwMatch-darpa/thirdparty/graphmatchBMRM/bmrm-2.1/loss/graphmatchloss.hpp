#ifndef _GRAPHMATCHLOSS_HPP_
#define _GRAPHMATCHLOSS_HPP_

#include "common.hpp"
#include "sml.hpp"
#include "info.hpp"
#include "data.hpp"
#include "graphdata.hpp"
#include "loss.hpp"
#include <string>
#include "timer.hpp"


/** Class for graphmatch classification loss.
 */
class CGraphMatchLoss : public CLoss
{
  public:
    CGraphMatchLoss(CModel* &model, CGraphData* &data);
    virtual ~CGraphMatchLoss();

    // Interfaces
    void Usage();
    void ComputeLoss(Scalar& loss)
    {
      throw CBMRMException("ERROR: not implemented!\n", "CGraphMatchLoss::ComputeLoss()");
    }
    void ComputeLossAndGradient(Scalar& loss, TheMatrix& grad);
    void Predict(CModel *model);
    void Evaluate(CModel *model);

    Scalar LabelLoss(int* y, int* ybar);

    void LoadModel(std::string modelFilename="");
    void SaveModel(std::string modelFilename="");

  protected:
    void Phi(int n, int* y, Scalar* res);

    CGraphData* _data;

    int quad;
};

#endif
