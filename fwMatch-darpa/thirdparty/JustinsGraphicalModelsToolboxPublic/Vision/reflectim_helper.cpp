/* Simple test program for the RAD package. */

/* Compute f(n,x) = the product of x[0], x[1], ..., x[n-1] */
/* and the gradient of f for two x values, and print the result. */

//#include "rad.hpp"
#include <stdio.h>
#include "mex.h"   //--This one is required
//#define nanoMatDebug 1 // 1 for safety/error checking; 0 for faster and unsafe(er)
//#include "nanoMatT.cpp"
#include "math.h"
//#define mmin(X,Y) ((X) < (Y) ? (X) : (Y))
//#define mmax(X,Y) ((X) > (Y) ? (X) : (Y))
#include "../eigenstuff.cpp"

/*
typedef nanoMat<double> dMat;
typedef nanoMat<int>    iMat;
typedef nanoMat<ADvar> adMat;
typedef ADvar var;
*/

// typedef nanoMat<double> dMat;
// typedef nanoMat<int>    iMat;
// typedef nanoMat<double> adMat;
// typedef double var;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    int i, j, k, l, m, n, x, y, x0, y0;
    
    MatrixMd im  = mapdouble(prhs[0]);
    int siz = (int)mapdouble(prhs[1])(0);
    
    plhs[0] = mxCreateDoubleMatrix(im.rows()+2*siz,im.cols()+2*siz,mxREAL);
    MatrixMd im2 = mapdouble(plhs[0]);
    
    for(y0=-siz; y0<im.rows()+siz; y0++){
        for(x0=-siz; x0<im.cols()+siz; x0++){
            x = x0;
            y = y0;
            // seemingly pointless while loop deals with points
            // where, e.g. siz < - im.M 
            while(1){
                if( x<0 )
                    x = -x-1;
                else if( x>= im.cols())
                    x = im.cols()*2-x-1;
                else
                    break;
            }
            while(1){    
                if( y<0 )
                    y = -y-1;
                else if( y>= im.rows())
                    y = im.rows()*2-y-1;
                else
                    break;
            }
            im2(y0+siz,x0+siz) = im(y,x);
        }
    }
    
    
    //printf("buf_size: %d  ly: %d  lx: %d  num_colors: %d \n", buf_size, ly, lx, num_colors);
    
    return;
    
}