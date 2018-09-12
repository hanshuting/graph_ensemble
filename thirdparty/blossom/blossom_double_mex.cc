#include <tr1/cstdint>
#include "mexcpp.h"
#include "print_options.h"
#include "PerfectMatching.h"

using namespace mexcpp;

enum {
  inNodes,
  iiVec,
  ijVec,
  iwVec,
  nI,
};

enum {
  oeVec,
  nO,
};

const char *usage = "Usage: eVec = blossom_double_mex(nNodes, iVec, jVec, wVec)";

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  if (nIn != 4) {
    mexErrMsgIdAndTxt("blossom_double_mex:args", usage);
  }

  // casting...
  int nNodes = scalar<double>(pIn[inNodes]);
  // MATLAB defaults to making everything double. We want double weights
  // (Kolmogorov's code can solve them), but iVec and jVec are just for
  // convenience.
  Mat<double> iVec(pIn[iiVec]);
  Mat<double> jVec(pIn[ijVec]);
  Mat<double> wVec(pIn[iwVec]);

  mxAssert(iVec.length == jVec.length && jVec.length == wVec.length,
           "Input vectors must be the same length.");

  size_t nEdges = iVec.length;

  // We really want it on the heap, don't we?
  PerfectMatching pm(nNodes, nEdges);

  // Structure for the checking call.
  // (Maybe it's so big it needs to be allocated on the heap... yeah, that
  // must be it. Damn C++.)
  //
  // This means your sparse version could also work, by allocating weights
  // and other structures on the heap.
  int *edges = new int[2* nEdges];
  int32_t ii, jj;
  for (int i = 0; i < nEdges; i++) {
    ii = iVec[i] - 1;
    jj = jVec[i] - 1;
    pm.AddEdge(ii, jj, wVec[i]);
    edges[2*i]     = ii;
    edges[2*i + 1] = jj;

    //mexPrintf("%s:%d -- Added (%d, %d, %d)\n", __FILE__, __LINE__, iVec[i] - 1, jVec[i] - 1, wVec[i]);
  }

  PerfectMatching::Options opts;

  //mexPrintf("%s:%d -- options\n", __FILE__, __LINE__);
  //mexPrintf(print_options(opts).c_str());
  pm.options = opts;

  //mexPrintf("%s:%d -- interrogation\n", __FILE__, __LINE__);
  //mexPrintf(interrogate(pm).c_str());

  //mexPrintf("%s:%d -- %d nodes, %d edges; solving\n", __FILE__, __LINE__, nNodes, nEdges);
  pm.Solve();
  int opt = CheckPerfectMatchingOptimality(nNodes, nEdges, edges, wVec.re, &pm);
  mxAssert(opt == 0, "Blossom did not solve to optimality.");
  delete edges;

//  mexPrintf("%s:%d -- checking\n", __FILE__, __LINE__, nEdges);
//  int res = ComputePerfectMatchingOptimality(nNodes, nEdges,

  Mat<bool> eVec(1, nEdges);
  for (int i = 0; i < nEdges; i++) {
    eVec[i] = (bool) pm.GetSolution(i);
  }

  pOut[oeVec] = eVec;
}
