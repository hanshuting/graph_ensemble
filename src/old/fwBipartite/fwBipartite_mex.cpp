#include "mexcpp.h"

// faster (inlined) rounding
#include <boost/math/special_functions/round.hpp>
#include <ctime>
#include <dlib/optimization.h>
#include <dlib/matrix.h>
#include <algorithm>
#include <vector>

using namespace mexcpp;

void updateLogBeliefs(const Mat<double> &beliefs, double *logBeliefs) {
  for (int i = 0; i < beliefs.length; i++) {
    if (beliefs[i] == 0 || beliefs[i] == 1) {
      logBeliefs[2*i]     = 0;
      logBeliefs[2*i + 1] = 0;
    } else {
#ifdef FASTLOG
      logBeliefs[2*i]     = fastlog(beliefs[i]);
      logBeliefs[2*i + 1] = fastlog(1 - beliefs[i]);
#else
      logBeliefs[2*i]     = log(beliefs[i]);
      logBeliefs[2*i + 1] = log(1 - beliefs[i]);
#endif
    }
  }
}

void computeRiPlusRjMinus1(int N, const double *rho, double *riPRjM1) {
  int Nsq = N * N;
  int i, j;
  for (int k = 0; k < Nsq; k++) {
    i = k / N;
    j = k % N;
    riPRjM1[k] = rho[i] + rho[j] - 1;
    //mexPrintf("riPRjM1(%d,%d) = %g\n", i, j, riPRjM1[k]);
  }
}

double evalObj(const Mat<double> &A,
               const Mat<double> &beliefs,
               const double *logBeliefs,
               const double *riPRjM1) {
  double obj = 0;
  for (int i = 0; i < beliefs.length; i++) {
    obj += - beliefs[i] * A[i]
           + beliefs[i] * logBeliefs[2*i]
           - riPRjM1[i]*(1 - beliefs[i])*logBeliefs[2*i + 1];
  }

  return obj;
}

bool absLess(double a, double b) {
  return fabs(a) < fabs(b);
}

double dlibLap(const Mat<double> &C, int *map, int64_t &cost) {
  // Normalize and round the cost function
  double absMaxC = fabs(*std::max_element(C.re, C.re + C.length, absLess));

  dlib::matrix<int64_t, 0, 0> dC(C.M, C.N);

  // Prepare variables for lap
//  mexPrintf("dlibLap C:\n");
  for (int i = 0; i < C.M; i++) {
    for (int j = 0; j < C.M; j++) {
      // Maximum cost problem.
      // Possibly, transpose?
      dC(i,j) = -boost::math::llround(C(j,i) / absMaxC * (INT64_MAX - 1)); 
//      mexPrintf("%lld ", dC(i,j));
    }
//    mexPrintf("\n");
  }

  clock_t start = clock();
  std::vector<long> mapVec(dlib::max_cost_assignment(dC));
  clock_t end = clock();

  cost = dlib::assignment_cost(dC, mapVec);

  std::copy(mapVec.begin(), mapVec.end(), map);

  return (end - start) / (CLOCKS_PER_SEC + 0.0);
}


void printMat(Mat<double> X) {
  for (int i = 0; i < X.M; i++) {
    for (int j = 0; j < X.N; j++) {
      mexPrintf("%g ", X(i,j));
    }
    mexPrintf("\n");
  }
}

void printArray(int N, int *a) {
  for (int i = 0; i < N; i++) {
    mexPrintf("%d ", a[i]);
  }
  mexPrintf("\n");
}

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  /*
  mexPrintf("pOut[0] = %x , pIn[0] = %x\n", pOut[0], pIn[0]);
  */
  SigintHandler sa;

  // Make a local copy
  Mat<double> beliefs(mxDuplicateArray(pIn[0]));
  Mat<double> A(pIn[1]);

  double *rho;
  if (nIn >= 3) {
    Mat<double> tmpRho(pIn[2]);
    if (tmpRho.length == A.M) {
      rho = tmpRho.re;

//      mexPrintf("fwBipartite: rho: ");
//      for (int i = 0; i < A.M; i++) {
//        mexPrintf("%g ", rho[i]);
//      }
//      mexPrintf("\n");

    } else {
      goto fillBethe;
    }
  } else {
fillBethe:
    rho = new double[A.M];
    std::fill(rho, rho + A.M, 1.0); // Bethe weights
  }

  double TolFun = 1e-6;
  if (nIn >= 4) {
    TolFun = scalar<double>(pIn[3]);
  }

  bool isDebug = false;
#ifndef NDEBUG
  if (nIn >= 5) {
    isDebug = scalar<bool>(pIn[4]);
  }
#endif

  int maxIter = INT_MAX;
  if (nIn >= 6) {
    maxIter = scalar<double>(pIn[5]);
  }

  if (A.M != A.N) {
    mexErrMsgIdAndTxt("fwBipartite_mex:A", "A must be square");
  }

  if ((A.M != beliefs.M) || (A.N != beliefs.N)) {
    mexErrMsgIdAndTxt("fwBipartite_mex:beliefs", "A and oldBeliefs must have the same shape");
  }

  // Own this array, so we destroy it ourselves when we're done.
  Mat<double> grad(A.M, A.M, true);

  int t = 0;
  std::vector<double> objs;
  std::vector<double> lapTimes;

#ifndef NDEBUG
  std::vector<mxArray *> grads;
  std::vector<Mat<int> > maps;
#endif

  double step;
  double obj;
//  int *map;
  int dlibMap[A.M];

  // Precompute log beliefs. They are expensive. Also, store log(b[i]) and log(1 - b[i])
  // consecutively.
  double logBeliefs[2 * A.length];
  double riPRjM1[A.length];
  computeRiPlusRjMinus1(A.M, rho, riPRjM1);

  int64_t lapCost, dlibLapCost;

  // NOTE: My signal handler thing might be bullshit.
  while (!(sa.interrupted) && (t < maxIter) && (t <= 1 || fabs(objs[objs.size() - 2] - objs[objs.size() - 1]) > TolFun)) {
//  while ((t < maxIter) && (t <= 1 || fabs(objs[objs.size() - 2] - objs[objs.size() - 1]) > TolFun)) {
    // Precompute logarithms and re-use for ent and dent; they are expensive.
    updateLogBeliefs(beliefs, logBeliefs);

    for (int k = 0; k < A.length; k++) {
      grad[k] = -A[k] + (riPRjM1[k] + 1) + logBeliefs[2*k] + riPRjM1[k]*logBeliefs[2*k + 1];
    }

    // Linear assignment
//    lapTimes.push_back(lap(grad, workArrays, lapCost));
//    map = workArrays.colsol;

    lapTimes.push_back(dlibLap(grad, dlibMap, dlibLapCost));
//    mxAssert(-lap
//    if (-lapCost != dlibLapCost) {
//      mexPrintf("lap cost = %lld and dlibLap cost = %lld\n", lapCost, dlibLapCost);
//      mexPrintf("lap cost = %lld and dlibLap cost = %lld not same\n", lapCost, dlibLapCost);
//      for (int i = 0; i < A.M; i++) {
//        mexPrintf("map[%d] = %d ; dlibMap[%d] = %d\n", i, map[i], i, dlibMap[i]);
//      }
//    }

#ifndef NDEBUG
    if (isDebug) {
      // Add grad and map to the debug log
      grads.push_back(mxDuplicateArray(grad.pm));

      Mat<int> mapMat(1, A.M, dlibMap);
      // Translate to 1-index
      for (int i = 0; i < mapMat.length; i++) {
        mapMat[i] += 1;
      }

      maps.push_back(mapMat);

      mexPrintf("grad: \n");
      printMat(grad);
      mexPrintf("map: \n");
      printArray(A.M, dlibMap);
    }
#endif

    objs.push_back(evalObj(A, beliefs, logBeliefs, riPRjM1));

    t++;
    step = 2.0 / (2.0 + t);

    // Update beliefs
    for (int i = 0; i < A.length; i++) {
      beliefs[i] *= 1 - step;
    }

#ifndef NDEBUG
    if (isDebug) {
      mexPrintf("post de-step, pre step\n");
      printMat(beliefs);
    }
#endif

    int iSel;
    for (int j = 0; j < A.M; j++) {
      iSel = dlibMap[j];
      beliefs(iSel,j) += step;
    }

#ifndef NDEBUG
    if (isDebug) {
      mexPrintf("post step\n");
      printMat(beliefs);

      if (t > 1) {
        mexPrintf("[MEX t=%d] objectives(end-1) = %g; objectives(end) = %g; diff = %g, TolFun = %g\n", t, objs[objs.size() - 2], objs[objs.size() - 1], fabs(objs[objs.size() - 2] - objs[objs.size() - 1]), TolFun);
      }
    }
#endif
  }

  // Evaluate last objective
  updateLogBeliefs(beliefs, logBeliefs);
  objs.push_back(evalObj(A, beliefs, logBeliefs, riPRjM1));

  //pOut[0] = scalar(0.0);
#ifndef NDEBUG
  if (isDebug) {
    mexPrintf("End beliefs\n");
    printMat(beliefs);
  }
#endif

  mexPrintf("nOut = %d\n", nOut);
  if (nOut >= 1) {
    pOut[0] = beliefs;
  }
  if (nOut >= 2) {
    pOut[1] = Mat<double>(objs, false);
  }
  if (nOut >= 3) {
    pOut[2] = Mat<double>(lapTimes, false);
  }

#ifndef NDEBUG
  if (isDebug && nOut >= 3) {
    int nIters = grads.size();
    StructMat debugLog(1, nIters, { "grad", "map" });
    for (int i = 0; i < nIters; i++) {
      Entry e = debugLog[i];
      e.set("grad", grads[i]);
      e.set("map",  maps[i]);
    }

    pOut[2] = debugLog;
  }
#endif
  /*
  mexPrintf("pOut[0] = %lx, beliefs.pm = %lx, pIn[0] = %lx, oldBeliefs.pm = %lx\n",
      pOut[0], beliefs.pm, pIn[0], oldBeliefs.pm);
      */
}

