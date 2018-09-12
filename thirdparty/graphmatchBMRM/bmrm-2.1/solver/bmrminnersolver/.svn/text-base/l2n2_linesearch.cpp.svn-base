
#ifndef _L2N2_LINESEARCH_CPP_
#define _L2N2_LINESEARCH_CPP_

#include "l2n2_linesearch.hpp"


void CL2N2_LineSearch::Solve(TheMatrix& w, TheMatrix& a, Scalar loss, Scalar &xi, Scalar &regval, Scalar &objval)
{
   iter++;

   Scalar slope = 0.0;
   Scalar gradientNorm2 = 0.0;
   Scalar w_dot_a = 0.0;
   Scalar H = 0.0;
   Scalar eta = 0.0;
   
   slope = loss + bmrmLambda*prevWNorm2Sq - r;
   a.Norm2(gradientNorm2);
   w.Dot(a,w_dot_a);
   H = bmrmLambda*bmrmLambda*prevWNorm2Sq + 2.0*bmrmLambda*w_dot_a + gradientNorm2*gradientNorm2;
   eta = (Scalar)std::min(1.0, (double)bmrmLambda*slope/H);

   //std::cout << "slope=" << slope << "   H=" << H << "   eta=" << eta << "   r=" << r << std::endl;

   w.Scale(1.0-eta);
   w.ScaleAdd(-eta/bmrmLambda, a);
   r = (1.0-eta)*r - eta*(loss - w_dot_a);
   w.Norm2(prevWNorm2Sq);
   prevWNorm2Sq = prevWNorm2Sq*prevWNorm2Sq;
   regval = 0.5*bmrmLambda*prevWNorm2Sq;
   xi = 0; //slope/2.0*eta;
   objval = 0; //xi + regval;
}

#endif
