#include <algorithm>
#include <cstdint>
#include <mexcpp.h>
#include "QPBO.h"

using namespace mexcpp;


// I should make it a lambda but I don't want to use C++11 for the
// platforms that don't support it.
#define IROUND(x) (boost::math::iround(scale * ((x) - lo)))

// Separated out in its own function to allow measuring overhead.
void fillQPBOModel(QPBO<COST_TYPE> &q,
                   int nNodes,
                   int nEdges,
                   int *QPBOEdges,
                   const Mat<COST_TYPE> &nodeCost,
                   const Mat<COST_TYPE> &edgeCost,
                   const Mat<int> &edgeList) {
  q.AddNode(nNodes);

  for (int n = 0; n < nNodes; n++) {
    /*
    mexPrintf("q.AddUnaryTerm(%d, %g, %g)\n",
        n,
        nodeCost(0, n),
        nodeCost(1, n));
        */

    q.AddUnaryTerm(n, nodeCost(0, n), nodeCost(1, n));
  }

  for (int e = 0; e < nEdges; e++) {
    /*
    mexPrintf("Edge (%d, %d) costs (0/00=%g, 1/01=%g, 2/10=%g, 3/11=%g)\n",
        int(edgeList(0, e)) - 1,
        int(edgeList(1, e)) - 1,
        edgeCost(0, e),
        edgeCost(1, e),
        edgeCost(2, e),
        edgeCost(3, e));
        */

    QPBOEdges[e] = q.AddPairwiseTerm(edgeList(0, e) - 1,
                                     edgeList(1, e) - 1,
                                     edgeCost(0, e),
                                     edgeCost(1, e),
                                     edgeCost(2, e),
                                     edgeCost(3, e));
  }
}

static COST_TYPE vals[] = {0, 0.5, 1.0};

// Lookup as exactIdxs[xj][xi] since x indexes columns, xi indexes rows.
static const int exactIdxs[][2] = { {0, 1}, {2, 3} };
// If xi = 0.5 but xj is 0 or 1, xiHalfIdxs[xj][0] and xiHalfIdxs[xj][1] give the two indices to set to 0.5
static const int xiHalfIdxs[][2] = { {0, 1}, {2, 3} };
static const int xjHalfIdxs[][2] = { {0, 2}, {1, 3} };

void recoverEdgeBeliefs(Mat<COST_TYPE> &edgeBelief,
                        int nEdges,
                        QPBO<COST_TYPE> &q,
                        const Mat<COST_TYPE> &edgeCost,
                        const Mat<int> &edgeList) {
  // Recover edge beliefs (harder).
  //
  // Down the rows of edgeCost (and thus also edgeBelief), we store the beliefs for
  // the (x_i, x_j) = [0] 00, [1] 01, [2] 10, [3] 11 configurations.
  //
  // The consistency constraints imply
  // edgeBelief(0,e) + edgeBelief(1,e) == nodeBelief(0,i)
  // edgeBelief(2,e) + edgeBelief(3,e) == nodeBelief(1,i)
  // edgeBelief(0,e) + edgeBelief(2,e) == nodeBelief(0,j)
  // edgeBelief(1,e) + edgeBelief(3,e) == nodeBelief(1,j)
  //
  // Due to some theory that I don't understand, the edge beliefs are also constrained
  // to lie in {0, 1/2, 1}.

  COST_TYPE Ediag, Eoffdiag;
  int xi, xj;
  for (int i, j, e = 0; e < nEdges; e++) {
    i = edgeList(0,e) - 1;
    j = edgeList(1,e) - 1;

    xi = q.GetLabel(i);
    xj = q.GetLabel(j);

    if ((xi > 0 && xi != 1) || (xj > 0 && xj != 1)) {
      mexErrMsgTxt("GetLabel returned something wrong.");
    }

    if (xi >= 0 && xj >= 0) {
      edgeBelief(exactIdxs[xj][xi],e) = 1;
    } else if (xi < 0 && xj >= 0) {
      edgeBelief(xiHalfIdxs[xj][0],e) = 0.5;
      edgeBelief(xiHalfIdxs[xj][1],e) = 0.5;
    } else if (xi >= 0 && xj < 0) {
      edgeBelief(xjHalfIdxs[xi][0],e) = 0.5;
      edgeBelief(xjHalfIdxs[xi][1],e) = 0.5;
    } else if (xi < 0 && xj < 0) {
      Ediag    = edgeCost(0,e) + edgeCost(3,e);
      Eoffdiag = edgeCost(1,e) + edgeCost(2,e);

      if (Ediag < Eoffdiag) {
        edgeBelief(0,e) = 0.5;
        edgeBelief(3,e) = 0.5;
      } else {
        edgeBelief(1,e) = 0.5;
        edgeBelief(2,e) = 0.5;
      }
    } else {
      mexErrMsgTxt("Unreachable clause in recoverEdgeBeliefs");
    }
  }
}

void mexFunction (int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  /////////////////////////////////////
  // Process arguments
  /////////////////////////////////////

  if (nOut < 2 || nIn != 3) {
    mexErrMsgTxt("[nodeBelief, edgeBelief, [energyLB]] = QPBO_mex(nodeCost, edgeCost, edgeList)");
  }

  const Mat<COST_TYPE> nodeCost(pIn[0]);
  const Mat<COST_TYPE> edgeCost(pIn[1]);
  const Mat<int> edgeList(pIn[2]);

  const int nNodes = nodeCost.N;
  const int nEdges = edgeCost.N;

  if (nodeCost.M != 2) { mexErrMsgTxt("nodeCost must be 2 x nNodes"); }
  if (edgeCost.M != 4) { mexErrMsgTxt("edgeCost must be 4 x nEdges"); }
  if (edgeList.M != 2) { mexErrMsgTxt("edgeList must be 2 x nEdges"); }
  if (edgeCost.N != edgeList.N) { mexErrMsgTxt("edgeCost and edgeList must have same number of columns"); }

  /////////////////////////////////////
  // Prepare to round
  /////////////////////////////////////

  COST_TYPE lo = std::min(nodeCost.min(), edgeCost.min());
  COST_TYPE hi = std::max(nodeCost.max(), edgeCost.max());

  /////////////////////////////////////
  // Create the model.
  /////////////////////////////////////
  // TOOD: Profile how long this construction takes relative to the solver.

  QPBO<COST_TYPE> q(nNodes, nEdges);

  int QPBOEdges[nEdges];
  fillQPBOModel(q, nNodes, nEdges, QPBOEdges, nodeCost, edgeCost, edgeList);

  /////////////////////////////////////
  // Solve the model.
  /////////////////////////////////////
  q.Solve();
//  q.ComputeWeakPersistencies();

  Mat<COST_TYPE> nodeBelief(2, nNodes);
  Mat<COST_TYPE> edgeBelief(4, nEdges);

  /////////////////////////////////////
  // Recover the fractional labels (beliefs).
  /////////////////////////////////////

  // Recover node beliefs (easy)
  for (int label, n = 0; n < nNodes; n++) {
    label = q.GetLabel(n);
//    mexPrintf("label[%d] = %d\n", n, label);
    if (label >= 0) {
      // Node n was correctly (and fully) labeled.
      nodeBelief(label,n) = 1;
    } else {
      // Node n was unlabled; we have a 1/2 fractional solution.
      nodeBelief(0,n) = 0.5;
      nodeBelief(1,n) = 0.5;
    }
  }

  recoverEdgeBeliefs(edgeBelief, nEdges, q, edgeCost, edgeList);

  pOut[0] = nodeBelief;
  pOut[1] = edgeBelief;
  if (nOut == 3) {
    pOut[2] = scalar(q.ComputeTwiceLowerBound() / 2);
  }
}

