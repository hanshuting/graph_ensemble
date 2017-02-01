#include <stdio.h>
#include <map>
//#include "etreeC.cpp"
//#include <adolc.h>
//#include <drivers/drivers.h>    // use of "Easy to Use" drivers
//#include <taping.h>
//#include "rad.hpp"
#include <Eigen/Eigen>
//#include <Eigen/Array>
#include <ctime>
#include <iostream>

using namespace std;
using namespace Eigen;

// here is a better version of ADOL-C (removes memory horror)
// http://www.sc.rwth-aachen.de/willkomm/adolcsc.html

// compilation problems!
// 1) needed to change permissions for stuff in autoconf
// 2) needed to copy malloc.h from source directory to include directory
// 3) compile like this:
// mex gradcomp_multiscale_helper_adol.cpp -L$HOME/adolc_base/lib -Ieigen/ -I$HOME/adolc_base/include/adolc/ -ladolc

/* MAKE ETREES WORK **************************
namespace Eigen {
    template<> struct NumTraits<etree> {
        typedef etree Real;
        typedef etree FloatingPoint;
        enum {
            IsComplex = 0,
            HasFloatingPoint = 1,
            ReadCost = 1,
            AddCost = 1,
            MulCost = 1
        };
    };
}
typedef Matrix<etree  ,Dynamic,Dynamic> MatrixXe;
/**/

/* MAKE ADVARS WORK **************************
namespace Eigen {
    template<> struct NumTraits<ADvar> {
        typedef ADvar Real;
        typedef ADvar FloatingPoint;
        enum {
            IsComplex = 0,
            HasFloatingPoint = 1,
            ReadCost = 1,
            AddCost = 1,
            MulCost = 1
        };
    };
}
// with with Eigen
ostream& operator<<(ostream& out, const ADvar& r){
    out << r.val();
    return out;
}
inline const ADvar& ei_conj(const ADvar& x)  { return x; }
inline const ADvar& ei_real(const ADvar& x)  { return x; }
inline ADvar ei_abs2(const ADvar& x)  { return x*x; }
inline ADvar ei_exp(const ADvar&  x)  { return exp(x); }
inline ADvar ei_log(const ADvar&  x)  { return log(x); }
inline ADvar ei_imag(const ADvar&)    { return 0.; }
//inline ADvar ei_abs(const ADvar&  x)  { return fabs(x); } //broken
inline ADvar ei_sqrt(const ADvar& x)  { return sqrt(x); }
inline ADvar ei_sin(const ADvar&  x)  { return sin(x); }
inline ADvar ei_cos(const ADvar&  x)  { return cos(x); }
inline ADvar ei_pow(const ADvar& x, ADvar y)  { return pow(x, y); }
// unary operator to get adjoints by "p.unaryExpr(ADadj())" for more info see:
// http://eigen.tuxfamily.org/dox/classEigen_1_1MatrixBase.html#1d1835f182c0141dc073a1045c8a2e9e
struct ADadj {
  const double operator()(const ADvar& x) const { return x.adj(); }
};
typedef Matrix<ADvar  ,Dynamic,Dynamic> MatrixXadvar;
/**/

/* MAKE adol-c WORK *******************
namespace Eigen {
template<> struct NumTraits<adouble>
{
  typedef adouble Real;
  typedef adouble FloatingPoint;
  enum {
    IsComplex = 0,
    HasFloatingPoint = 1,
    ReadCost = 1,
    AddCost = 1,
    MulCost = 1
  };
};
}
// the Adolc's type adouble is defined in the adtl namespace
// therefore, the following ei_* functions *must* be defined
// in the same namespace
namespace adtl {
  inline const adouble& ei_conj(const adouble& x)  { return x; }
  inline const adouble& ei_real(const adouble& x)  { return x; }
  inline adouble ei_imag(const adouble&)    { return 0.; }
  inline adouble ei_abs(const adouble&  x)  { return fabs(x); }
  inline adouble ei_abs2(const adouble& x)  { return x*x; }
  inline adouble ei_sqrt(const adouble& x)  { return sqrt(x); }
  inline adouble ei_exp(const adouble&  x)  { return exp(x); }
  inline adouble ei_log(const adouble&  x)  { return log(x); }
  inline adouble ei_sin(const adouble&  x)  { return sin(x); }
  inline adouble ei_cos(const adouble&  x)  { return cos(x); }
  inline adouble ei_pow(const adouble& x, adouble y)  { return pow(x, y); }
}
typedef Matrix<adouble,Dynamic,Dynamic> MatrixXa;
/**/


typedef Map<MatrixXd> MatrixMd;
typedef Map<MatrixXi> MatrixMi;

Eigen::Map<Eigen::MatrixXd> mapdouble(const mxArray* matlab_array){
    // error checking!  amazing!
    if(mxGetClassID(matlab_array)!=mxDOUBLE_CLASS){
        cout << "Trying to wrap non-double array" << endl;
        throw "Trying to wrap non-double array";
    }
    double* data = (double*)mxGetData(matlab_array);
    int rows  = mxGetM(matlab_array);
    int cols  = mxGetN(matlab_array);
    return Eigen::Map < Eigen::MatrixXd >(data, rows, cols);
}

Eigen::Map<Eigen::MatrixXi> mapint32(const mxArray* matlab_array){
    if(mxGetClassID(matlab_array)!=mxINT32_CLASS){
        cout << "Trying to wrap non-int32 array" << endl;
        throw "Trying to wrap non-int array";
    }
    int* data = (int*)mxGetData(matlab_array);
    int rows  = mxGetM(matlab_array);
    int cols  = mxGetN(matlab_array);
    return Eigen::Map < Eigen::MatrixXi >(data, rows, cols);
}

// intended for eigen matrix types
template <typename T>
void norm_cols(T& A){
//     cout << "A:\n" << A << endl;
//     for(int i=0; i<A.cols(); i++){
//         double mysum = 0;
//         for(int j=0; j<A.rows(); j++){
//             cout << "i: " << i << " j: " << j << " A[j,i]: " << A[j,i] << endl;
//             mysum += A[j,i];
//         }
//         cout << "mysum: " << mysum << endl;
//         for(int j=0; j<A.rows(); j++)
//             A[j,i] /= mysum;
//     }
    for(int i=0; i<A.cols(); i++)
        A.col(i) = A.col(i) / A.col(i).sum();
}

// actually, this should work for any matrix type!
template <typename T>
T imresize_bigger(const T& A){
    int M  = A.rows();
    int N  = A.cols();
    int M2 = M*2;
    int N2 = N*2;
    
        // want f(M2-1)=M-1-r
        //      f(0)   = r
        //      f(x)   = a*x+b
        //      b      = r;
        //      a*(M2-1)+r = M-1-r
        //      a          = (M-1-2*r)/(M2-1)
    
    T C = T(M2,N2);
    for(int y2=0; y2<M2; y2++){
        for(int x2=0; x2<N2; x2++){
            double r = -.25;
            double a = (M-1-2*r)/(M2-1);
            double b = r;
            double y = a*y2+b;
            a = (N-1-2*r)/(N2-1);
            b = r;
            double x = a*x2+b;
            if(y<0) y=0;
            if(x<0) x=0;
            if(y>M-1) y=M-1;
            if(x>N-1) x=N-1;
            int y_lo = floor(y);
            int x_lo = floor(x);
            int y_hi = ceil(y);
            int x_hi = ceil(x);
            double yfrac = 1-(y-y_lo);
            double xfrac = 1-(x-x_lo);
            
            C(y2, x2) =    yfrac *   xfrac *A(y_lo,x_lo) + 
                        (1-yfrac)*   xfrac *A(y_hi,x_lo)  + 
                           yfrac *(1-xfrac)*A(y_lo,x_hi)  + 
                        (1-yfrac)*(1-xfrac)*A(y_hi,x_hi);
        }
    }
    
    return C;
}

template <typename T>
T imresize_smaller(T A){
    int M = A.rows();
    int N = A.cols();
    int M2 = M/2;
    int N2 = N/2;
    
    //cout << "M2: " << M2 << "  N2: " << N2 << endl;
    
    T C = T(M2,N2);
    
    for(int y2=0; y2<M2; y2++){
        for(int x2=0; x2<N2; x2++){
            int y_lo = y2*2;
            int x_lo = x2*2;
            int y_hi = y_lo+1;
            int x_hi = x_lo+1;
            C(y2, x2) = .25*A(y_lo,x_lo) + 
                        .25*A(y_hi,x_lo) + 
                        .25*A(y_lo,x_hi) + 
                        .25*A(y_hi,x_hi);
        }
    }
    return C;
}

template <typename T>
T rescale(T A, int diff){
    if(diff==0)
        return A;
    if(diff<0) // make bigger
        return rescale(imresize_bigger(A),diff+1);
    else
        return rescale(imresize_smaller(A),diff-1);
}

// templated function to print an arbitrary list of crap
template <typename T>
void printStuff(T first, T last){
     for(; first != last; ++first)
         cout << *first << endl;
} //use like: printStuff(ingredients.begin(), ingredients.end()

// // this is done in a somewhat stupid way returning a vector just because
// // I don't feel like remembering how to use templates properly
// template <typename T>
// T log_sum_exp(T A){
//     T rez = T(1,1);
//     T damin = T(1,1);
//     damin(0) = A.minCoeff();
//     rez(0) = 0;
//     for(int i=0; i<A.rows()*A.cols(); i++)
//         rez(0) += exp(A(i)-damin(0));
//     rez(0) = log(rez(0))+damin(0);
//     return rez;
// }

inline double log_sum_exp(const MatrixXd& A){
    int nvals = A.rows()*A.cols();
    
    double rez   = 0;
    double damin = A.maxCoeff();
    
    for(int i=0; i<nvals; i++)
        rez += exp(A(i)-damin);
    rez = log(rez)+damin;
    return rez;    
}

inline double log_sum_exp(const MatrixXd& A, int nvals){
    double rez   = 0;
    double damin = A.maxCoeff();
    
    for(int i=0; i<nvals; i++)
        rez += exp(A(i)-damin);
    rez = log(rez)+damin;
    return rez;    
}


// void print_tapestats(){
//     int counts[11];
//     tapestats(1, counts);
//     cout << "#ind:  " << counts[0] << endl;
//     cout << "#dep:  " << counts[1] << endl;
//     cout << "#live: " << counts[2] << endl;
//     cout << "#vstk: " << counts[3] << endl;
//     cout << "#bsiz: " << counts[4] << endl;
//     cout << "#ops:  " << counts[5] << endl;
//     cout << "#6:    " << counts[6] << endl;
//     cout << "#7:    " << counts[7] << endl;
//     cout << "#8:    " << counts[8] << endl;
//     cout << "#9:    " << counts[9] << endl;
//     cout << "#10:   " << counts[10] << endl;    
// }

// template <typename T>
// vector<MatrixXd> niceGrad(T first, T last){
//     // first, count the total number of elements
//     int numels = 0;
//     for(T it=first; it != last; ++it)
//         numels += it->rows()*it->cols();
//     //cout << "numels: " << numels << endl;
//     // declare a big vector to stick everything in
//     double x[numels];
//     double g[numels];
//     int el=0;
//     for(T it=first; it != last; ++it)
//         for(int i=0; i<it->rows()*it->cols(); i++)
//             x[el++] = (*it)(i);
//     // call gradient function
//     //double f;
//     //function(1,1,numels,x,&f);  
//     clock_t start = clock();
//     gradient(1  ,numels,x,g);
//     cout << "time3: " << (double(clock())-double(start))/CLOCKS_PER_SEC << endl;
//     //for(int i=0; i<numels; i++)
//     //    cout << x[i] << " " << g[i] << endl;
//     // create output matrices holding gradient values
//     vector<MatrixXd> vec;
//     el = 0;
//     for(T it=first; it != last; ++it){
//         // new matrix of same size
//         MatrixXd M(it->rows(),it->cols());
//         for(int i=0; i<it->rows()*it->cols(); i++)
//             M(i) = g[el++];
//         vec.push_back(M);
//     }
//     return vec;
// }
// 
// template <typename T>
// MatrixXa adMat(T& x0){
//     //MatrixXa x = x0.cast<adouble>(); // not sure why doesn't work
//     MatrixXa x(x0.rows(),x0.cols());
//     for(int i=0; i<x.rows()*x.cols(); i++){
//         x[i] <<= x0[i];
//         //x[i] = x0[i];
//         //x[i].declareIndependent();
//     }
//     return x;
// }