#include <tr1/cstdint>
//#include <cstdio>
#include "mexcpp.h"
//#include "print_options.h"

#include "PerfectMatching.h"

using namespace mexcpp;

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  if (nIn != 1) {
    mexErrMsgIdAndTxt("blossom_sparse_mex:args", "Should only be called from blossom.m");
  }

  // A IS ASSUMED TO BE SYMMETRIC. Further, we assume nnz = nzMax.
  // In particular, no extraneous zeros.
  SparseMat A(const_cast<mxArray *>(pIn[0]));

  mxAssert(A.N == A.M, "A must be square and symmetric.");
  int nNodes = A.N;
  mxAssert(nNodes % 2 == 0,  "A must be even rows/cols.");
  mxAssert(A.nzMax % 2 == 0, "A must have even number of rows.");
  int maxEdges = A.nzMax;

  PerfectMatching pm(nNodes, maxEdges);

  mwIndex i;
  PerfectMatching::NodeId ii, jj;

  // These must be allocated on the heap.
  int *edges = new int[2 * maxEdges];

  // For some reason we crash without this copy. Whatever, keep it safe.
  //double weights[maxEdges];
  int iEdge = 0;

  mwIndex aIdx;
  for (mwIndex j = 0; j < nNodes; j++) {
    for (aIdx = A.jc[j]; aIdx < A.jc[j+1]; aIdx++) {
      i = A.ir[aIdx];

      // Upper triangular only
      if (j > i) {
        // Translate 1 to 0 indexing.
        // weights[iEdge] = A.pr[aIdx];
        pm.AddEdge(i, j, A.pr[aIdx]);
//        mexPrintf("%s:%d -- (i=%d, j=%d), iEdge=%d, aIdx=%d; weight = %g\n", __FILE__, __LINE__, i, j, iEdge, aIdx, weights[iEdge]);
        // Checking structure
        edges[2*aIdx]     = i;
        edges[2*aIdx + 1] = j;
        //iEdge++;
      }
    }
  }

//  int nEdges = iEdge - 1;
  int nEdges = aIdx - 1;

  PerfectMatching::Options opts;

  //mexPrintf("%s:%d -- options\n", __FILE__, __LINE__);
  //mexPrintf(print_options(opts).c_str());
  pm.options = opts;

  //mexPrintf("%s:%d -- interrogation\n", __FILE__, __LINE__);
  //mexPrintf(interrogate(pm).c_str());

  //mexPrintf("%s:%d -- %d nodes, %d edges; solving\n", __FILE__, __LINE__, nNodes, maxEdges);
  pm.Solve();
  int opt = CheckPerfectMatchingOptimality(nNodes, nEdges, edges, A.pr, &pm);
  mxAssert(opt == 0, "Blossom did not solve to optimality.");

  std::vector<double> is;
  std::vector<double> js;

  for (int iEdge = 0; iEdge < nEdges; iEdge++) {
    if (pm.GetSolution(iEdge)) {
      // i is part of the matching
      // Translate 0 to 1 indexing
      is.push_back(edges[2*iEdge] + 1);
      js.push_back(edges[2*iEdge + 1] + 1);
    }
  }

  // Prepare a sparse call (true = I own the memory)
  mxAssert(is.size() == js.size(), "iVec and jVec sizes differed!");
  int matchingEdges = is.size();
  //fprintf(stderr, "%s:%d matchingEdges = %d\n", __FILE__, __LINE__, matchingEdges);
  Mat<double> iVec(is, true, true);
  Mat<double> jVec(js, true, true);
  Mat<double> wVec(matchingEdges, 1, true);
  std::fill(wVec.re, wVec.re + matchingEdges, 1.0);

  SparseMat ret(iVec, jVec, wVec, nNodes, nNodes);
  pOut[0] = ret;
  if (nOut == 2) {
    pOut[1] = scalar(ComputePerfectMatchingCost(nNodes, nEdges, edges, A.pr, &pm));
  }

}

