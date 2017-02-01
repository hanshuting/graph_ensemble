#include <algorithm>
#include <tr1/cstdint>
#include <mexcpp.h>
#include "csa.hh"
#define BOOST_ALL_DYN_LINK
#include <boost/math/special_functions/round.hpp>

using namespace mexcpp;

void mexFunction (int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  if (nOut != 1 || nIn != 1) {
    mexErrMsgTxt("permMat = csaAssignPermMat(cost)");
  }

  const Mat<double> cost(pIn[0]);
  const int N = cost.M;
  const int Nsq = N * N;

  mxAssert(N == cost.N, "cost must be square.");

  Mat<double> permMat(N, N);

  // Get ready to round
  double lo = *std::min_element(cost.re, cost.re + Nsq);
  double hi = *std::max_element(cost.re ,cost.re + Nsq);
//  mexPrintf("lo = %g, hi = %g\n", lo, hi);

  double range = hi - lo;

  double scale;
  if (range >= 1.0) {
    scale = (INT32_MAX - 1) / range;
  } else {
    scale = double(INT32_MAX - 1);
  }

  // Make things smaller coz costs will be added and whatnot
  scale *= 1e-4;

  int graph[3 * Nsq];
  int i = 0;
  int graphCost;
  int graphMaxC = 0;
  for (int c = 1; c <= N; c++) {
    for (int r = 1; r <= N; r++) {
      // Sparsify inlined. Edge labels must be 1-indexed.
      graph[3*i+0] = r;
      graph[3*i+1] = N + c;
      //mexPrintf("pre_round: cost[%d] - lo = %g\n", i, cost[i] - lo);
      graphCost = boost::math::iround(scale * (cost[i] - lo));
      graph[3*i+2] = graphCost;
      graphMaxC = std::max(graphMaxC, graphCost);

      // DEBUG
      //mexPrintf("graph[%d, %d, %d] = %d, %d, %d\n", 3*i+0, 3*i+1, 3*i+2, r, N + c, graphCost);
      i++;
    }
  }

  // The CSA package segfaults if all the edge weights are zero.
  // In this case, set all the weights to one, and then later
  // remember to set the returned graph weights back to zero.
  if (graphMaxC == 0) {
    for (i = 0; i < Nsq; i++) {
        graph[3*i+2] = 1;
    }
  }

  // KT: Note that we aren't even returning the cost, so no need to
  // check for maxc == 0 later.

  // Run CSA.  It will either run successfully or segfault or loop
  // forever or return garbage.  But it claims to always return a
  // valid result if the input is valid.  The checks above try to
  // ensure that the input is ok, but I don't check that a perfect
  // match is present (which CSA requires but does not check for,
  // grumble, grumble).  In that case, you're on your own, since
  // I'm not sure how to quickly check for that condition.
  CSA csa (2*N, Nsq, graph);
  int e = csa.edges();
//  mexPrintf("e = %d\n", e);

  int r, c, junkCost;
  for (i = 0; i < e; i++) {
    csa.edge(i, r, c, junkCost);
    // Column index c has an extra N added, and also 1-indexed.
    //mexPrintf("i = %d, r = %d, c - N - 1 = %d\n", i, r, c - N - 1);
    permMat(r - 1,c - N - 1) = 1;
  }

  /*
  mexPrintf("permMat: \n");
  for (int r = 0; r < N; r++) {
    for (int c = 0; c < N; c++) {
      mexPrintf("%d ", permMat(r,c));
    }
    mexPrintf("\n");
  }
  */

  pOut[0] = permMat;
}

