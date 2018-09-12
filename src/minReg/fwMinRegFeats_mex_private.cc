// [theta, tau, objs] = fwMinRegFeats_mex_private(N, Gty, G, lamba, rho, params)
// N      - Number of items in the permutation
// Gty    - Kx1 vector of inner productsr; precomputed.
// G      - (MN^2)xK matrix of features, again accessed column-wise (contiguously)
// lambda - scalar double regularization parameter
// params - structure matrix of param-value options.
//
#include <mexcpp.h>
#include <csa.hh>

// TODO: Alternative sink for boost logging.

// faster (inlined) rounding
#include <algorithm>
#include <Eigen/Dense>
#include <cstdint>
#include <cmath>
#include <ctime>
#include <sstream>
#include <vector>

#define BOOST_ALL_DYN_LINK
#include <boost/math/special_functions/round.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>

using namespace mexcpp;
using namespace Eigen;

// Correct order is ROW MAJOR
void initCsaGraph(int N, int *csaGraph) {
  int k = 0;
  for (int c = 0; c < N; c++) {
    for (int r = 0; r < N; r++) {
      /* other order
      csaGraph[3*k + 0] = r + 1;
      csaGraph[3*k + 1] = N + c + 1;
      */
      csaGraph[3*k + 0] = c + 1;
      csaGraph[3*k + 1] = N + r + 1;
      k++;
    }
  }
}

void csaAssignPerm(int N, ArrayXd &grad, int *csaGraph, int *map) {
  // ROUND IT!!!
  grad -= grad.minCoeff();
  double hi = grad.maxCoeff();
  double scale;
  if (hi >= 1.0) {
    scale = (INT32_MAX - 1) / hi;
  } else {
    scale = double(INT32_MAX - 1);
  }

  BOOST_LOG_TRIVIAL(debug) << "INT32_MAX - 1 = " << INT32_MAX - 1 << " hi = " << hi << " scale = " << scale;

  // For safety... sometimes the edge numbers change.
  initCsaGraph(N, csaGraph);

  int c;
  for (int i = 0; i < N*N; i++) {
    BOOST_LOG_TRIVIAL(debug) << "i = " << i << " scale = " << scale << " grad(i) = " << grad(i);
    c = boost::math::iround(scale * grad(i));
    mxAssert(c >= 0, "We subtracted grad.minCoeff(), so all elements should be nonnegative");
    csaGraph[3*i + 2] = c;

  }

  for (int i = 0; i < N*N; i++) {
    BOOST_LOG_TRIVIAL(debug) << "csaGraph i = " << i
                             << ": col = " << csaGraph[3*i]
                             << " row = " << csaGraph[3*i + 1]
                             << " cost = " << csaGraph[3*i + 2];
  }

  BOOST_LOG_TRIVIAL(debug) << "Right before calling CSA";
  // CSA n = number of nodes, which are twice our N (bipartite)
  // CSA m = number of edges, e.g. N2.
  CSA csa(2*N, N * N, csaGraph);

  // Sentinel
  std::fill(map, map + N, -1);

  int r, a, b, cost;
  // We extract a column permutation. We initialized the graph as column major,
  // where row indices are 1:N and column indices are (N+1):2N. We want
  // map[c] = r to denote the edge (r,c)
  //
  // TODO: EXPLAIN MORE CLEARLY; THIS DOESN'T MAKE SENSE TO ANYONE BUT YOU.
  for (int i = 0; i < csa.edges(); i++) {
    // mxPrintF blah blah blah
    csa.edge(i, a, b, cost);
    r = a - 1;
    c = b - N - 1;
    mxAssert(r >= 0 && r < N, "Row out of bounds from CSA assignment");
    mxAssert(c >= 0 && c < N, "Col out of bounds from CSA assignment");
    map[c] = r;
  }

  mxAssert(std::find(map, map + N, -1) == map + N, "map was not fully assigned.");
}

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  boost::log::core::get()->set_filter(
      boost::log::trivial::severity >= boost::log::trivial::warning
      );

  SigintHandler sa;

  //////////////////////////////////////////////////////
  // INPUTS
  //////////////////////////////////////////////////////

  const int N = scalar<double>(pIn[0]);
  mexPrintf("N\n");
  const Mat<double>::EigenMatrixMap Gty(Mat<double>(pIn[1]).asEigenMatrixMap());
  mexPrintf("Gty\n");
  const Mat<double>::EigenMatrixMap G(Mat<double>(pIn[2]).asEigenMatrixMap());
  mexPrintf("G\n");

  double lambda = scalar<double>(pIn[3]);
  mexPrintf("lambda\n");
  const Mat<double>::EigenMatrixMap rho(Mat<double>(pIn[4]).asEigenMatrixMap());
  mexPrintf("rho\n");
  StructMat       params(pIn[5]);
  mexPrintf("params\n");

  // Parse dimensions
  const int K   = Gty.size();
  const int MN2 = G.rows();
  const int N2  = N * N;
  const int M   = MN2 / N2;

  // Dimensionality consistency checks
  mxAssert(M * N2 == MN2, "MN2 is not M * N2");

  // Additional parameters
  double  TolGap  = params.getS<double>("TolGap");
  int32_t maxIter = params.getS<double>("maxIter");
  bool    debug   = params.getS<bool>("debug");

  BOOST_LOG_TRIVIAL(debug) << "N = " << N
                           << ", K = " << K
                           << "MN2 = " << MN2
                           << "N2 = " << N2
                           << "M = " << M
                           << "debug = " << debug;

  //////////////////////////////////////////////////////
  // OUTPUTS
  //////////////////////////////////////////////////////
  Mat<double> mlabTheta(K, 1);
  Mat<double> mlabTau(MN2, 1);

  pOut[0] = mlabTheta;
  pOut[1] = mlabTau;

  Mat<double>::EigenMatrixMap theta(mlabTheta.asEigenMatrixMap());
  Map<ArrayXd> tau(mlabTau.re, MN2, 1);

  // This is dynamic, so we don't allocate the Mat until end.
  std::vector<double> objsVec;

  //////////////////////////////////////////////////////
  // INITIALIZE FRANK-WOLFE
  //////////////////////////////////////////////////////
  ArrayXd riPRjM1(N2, 1);
  int ind;
  // TODO: Make this line not segfault.
  // ArrayXd riPRjM1 = rho.replicate(1, N) + rho.transpose().replicate(N, 1);
  for (int c = 0; c < N; c++) {
    for (int r = 0; r < N; r++) {
      ind = N*c + r;
      riPRjM1(ind) = rho(r) + rho(c) - 1;
    }
  }

  /*
  std::stringstream ss;
  ss << riPRjM1;
  BOOST_LOG_TRIVIAL(debug) << "riPRjM1 = " << ss.str();
  */

  // Initialize tau and logTau
  tau.fill(1.0 / N);

  ArrayXd logTau(MN2, 1);
  ArrayXd log1MTau(MN2, 1);

  logTau = tau.log();
  log1MTau = (1 - tau).log();

  ArrayXd grad(N2, 1);
  MatrixXd momentDiffs(K, 1);
  double dualityGap = INFINITY;
  double sGrad;
  double xGrad;
  double entropy;
  double step;
  double obj;
  int t = 0;
  int map[N];
  int csaGraph[3*N];
  initCsaGraph(N, csaGraph);

  int GStartInd;
  int tauStartInd;
  const double *GStartPtr;

  // Initialize momentDiffs
  momentDiffs = Gty - G.transpose()*tau.matrix();

  // NOTE: My signal handler thing might be bullshit.
  while (!(sa.interrupted) &&
         (t < maxIter) &&
         (dualityGap > TolGap)) {
    dualityGap = 0;

    // Deterministic step-size. The literature writes 2/(2+t), but that's for
    // 1-based indexing.
    step = 2.0 / (3.0 + t);

    // Shrink all beliefs by (1 - step)
    tau *= (1 - step);

    // Warning: Henceforth, tau has an invalid value.

    //////////////////////////////////////////////////////
    // COMPUTE PER-SAMPLE FRANK-WOLFE STEPS
    //////////////////////////////////////////////////////
    for (int m = 0; m < M; m++) {
      grad.setZero(N2, 1);

      // Compute gradient of data term.
      for (int k = 0; k < K; k++) {
        // MN^2 rows x K cols. Each of M samples gets a block of N^2 rows.
        // P.block(i, j, rows, cols)          // P(i+1 : i+rows, j+1 : j+cols)
        grad = grad - (momentDiffs(k) * G.array().block(N2*m, k, N2, 1)) / lambda;
      }

      // vectorize
      grad = grad + (riPRjM1 + 1) + logTau.segment(N2*m, N2) +
              riPRjM1 * (log1MTau.segment(N2*m, N2));

      // Solve the LP.
      csaAssignPerm(N, grad, csaGraph, map);
      for (int c = 0; c < N; c++) {
        BOOST_LOG_TRIVIAL(debug) << "map[" << c << "] = " << map[c];
      }

      // Update beliefs. Note the we already shrank it before the m-loop, so
      // we need only add step in the map direction.
      int r;
      sGrad = 0;
      for (int c = 0; c < N; c++) {
        r = map[c];
        tau(N2*m + N*c + r) += step;
        sGrad += grad(N*c + r);
      }

      BOOST_LOG_TRIVIAL(debug) << "tau.segment(N2*m, N2).rows() = " <<
                                  tau.segment(N2*m, N2).rows() <<
                                  " tau.segment(N2*m, N2).cols() = " <<
                                  tau.segment(N2*m, N2).cols() <<
                                  " grad = " <<
                                  grad;

      xGrad = (tau.segment(N2*m, N2) * grad).sum();

      dualityGap += xGrad - sGrad;
    }

    // Now all \tau are valid again. Update logTau.
    logTau = tau.log();
    log1MTau = (1 - tau).log();

    entropy = (tau * logTau + (riPRjM1.replicate(M, 1)) * (1 - tau) * log1MTau).sum();

    momentDiffs = Gty - G.transpose()*tau.matrix();

    BOOST_LOG_TRIVIAL(warning) << " tau' = " << tau.transpose();
    BOOST_LOG_TRIVIAL(warning) << " entropy " << entropy;

    obj = (0.5 * lambda) * momentDiffs.array().square().sum() + entropy;

    objsVec.push_back(obj);

    if (debug) {
      BOOST_LOG_TRIVIAL(debug) << "obj = " << objsVec.back() << " dualityGap = " << dualityGap << " TolGap = " << TolGap;
      mexPrintf("[.MEX t=%d] objectives(end) = %g; dualityGap = %g, TolGap = %g\n", t, objsVec.back(), dualityGap, TolGap);
    }

    t++;
  }

  pOut[2] = Mat<double>(objsVec);

  // Recover \theta from the LAST iterate (not necessarily the best but easier this way)
  theta = momentDiffs / lambda;
}

