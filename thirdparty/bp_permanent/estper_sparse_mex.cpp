#include "mexcpp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

using namespace mexcpp;

#define NEG_INF -9999
#define MAX_ITER 1000000000
#define CONVERGENCE_INSURANCE 1
#define THRESHOLD 1e-6

double damp = .5;

typedef struct node_t_struct {
  int id;
  int neighborhood; //number of neighbors (not counting itself)
  struct node_t_struct **neighbors; //vector of pointers to neighbors
  double **message; //vector of messages from this node to all other nodes
  double *f;
  double *g;
  double **prevmessage;
  double *phi;   //potential
} node_t;

double safelog(double x) {
  /* returns log of x + smallest possible double */
    if (x==0) {
     return log(x+DBL_MIN);
    }
  return log(x);
}

void printMatrix(Mat<double> &W, int n)
     /*prints matrix for debugging purposes*/
{
  int i,j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      mexPrintf("%f ", W(i,j));
    }
    mexPrintf("\n");
  }
}

node_t * initNode(int n, int id, int total)
{
  int i,j;
  node_t *node;
  node = (node_t *)malloc(sizeof(node_t));
  node->id = id;
  node->neighborhood = n;
  node->phi = (double*)malloc(n*sizeof(double));
  node->f = (double*)malloc(n*sizeof(double));
  node->g = (double*)malloc(n*sizeof(double));
  node->neighbors = (node_t**)malloc(n*sizeof(node_t*));
  node->message = (double**)malloc(n*sizeof(double*));
  node->prevmessage = (double**)malloc(n*sizeof(double*));

  for (i=0; i<total; i++) {
    node->message[i] = (double*)malloc(n*sizeof(double));
    node->prevmessage[i] = (double*)malloc(n*sizeof(double));
    for (j=0; j<total; j++) {
      node->message[i][j]=1;
      node->prevmessage[i][j]=1;
    }
  }
  for (i=0; i<n; i++) {
    node->phi[i]=1;
    node->f[i]=1;
    node->g[i]=0;
  }

  return node;
}

void freeNode(node_t *node, int n)
{
  int i;
  free(node->phi);
  free(node->f);
  free(node->g);
  free(node->neighbors);

  for (i=0; i<n; i++) {
    free(node->message[i]);
    free(node->prevmessage[i]);
  }
  free(node->message);
  free(node->prevmessage);
}

void updateF(node_t *node, int n) {
  /*updates f function*/
  int i,j;
  //compute full products

  for (i=0; i<n; i++) {
    node->f[i] = node->phi[i];
    for (j=0; j<n; j++) {
      node->f[i] *= node->neighbors[j]->message[node->id][i];
    }
  }
}

void updateG(node_t *node, int n) {
  /*update g value. Assumes f is up to date */
  int i,j;
  //compute full summation \sum_{x_i} \frac{f(x_i)}{m_{y_j}(x_i)}

  for (i=0; i<n; i++) {
    node->g[i] = 0;
    for (j=0; j<n; j++) {
      if (node->neighbors[i]->f[j]>DBL_MIN) {
	node->g[i] +=
	  node->neighbors[i]->f[j]/node->message[node->neighbors[i]->id][j];
      }
    }
  }
}

void saveMessages(node_t *node, int n)
     /*saves current messages in prevmessages matrix*/
{
  int i, j;
  //save previous messages
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      node->prevmessage[i][j] = node->message[i][j];
    }
  }
}

void updateMessages(node_t *node, int n)
     /*updates the messages stored in the given node
       by using messages from neighbors in opposite
     Also compute change from previous message*/
{
  int i,j;
  double sum, fraction;

  for (i=0; i<n; i++) {

    //check for zero/zero, set to zero;

    /*if (node->neighbors[i]->prevmessage[node->id][node->neighbors[i]->id]<
	THRESHOLD) {
      mexPrintf("%f/%f\n",
	     node->f[node->neighbors[i]->id],
	     node->neighbors[i]->prevmessage[node->id][node->neighbors[i]->id]);mexPrintf("numerator below threshold? %d\n", node->f[node->neighbors[i]->id]<THRESHOLD);
	     }*/

    if (node->f[node->neighbors[i]->id]< DBL_MIN) {
      fraction = 0;
      //mexPrintf("set to zero");
    } else {
      fraction = node->f[node->neighbors[i]->id]/
	node->neighbors[i]->prevmessage[node->id][node->neighbors[i]->id];
      //mexPrintf("computed");
    }

    //update message from current node to its ith neighbor
    for (j=0; j<n; j++) {
      if (j==node->id) {
	node->message[node->neighbors[i]->id][j] = fraction;
      } else {
	node->message[node->neighbors[i]->id][j] =
	  node->neighbors[i]->g[node->id] - fraction;
      }
    }

    //normalize message:
    sum = 0;
    for (j=0; j<n; j++) {
      sum+=node->message[node->neighbors[i]->id][j];
    }
    for (j=0; j<n; j++) {
      if (sum>0) {
	node->message[node->neighbors[i]->id][j]/=sum;
      }
    }

    if (isnan(sum)) {
      mexPrintf("Warning. sum of message is %f\n", sum);
      mexPrintf("Message from %d to %d\n",
	     node->id, node->neighbors[i]->id);
      //while(1){}
    }

    //dampen message:
    for (j=0; j<n; j++) {
      if (~isnan(node->message[node->neighbors[i]->id][j])) {
      //mexPrintf("M: %f\n", node->message[node->neighbors[i]->id][j]);

	node->message[node->neighbors[i]->id][j]=
	  exp(safelog(node->prevmessage[node->neighbors[i]->id][j])+
	      damp*(safelog(node->message[node->neighbors[i]->id][j])
		    -safelog(node->prevmessage[node->neighbors[i]->id][j])));
	if (isnan(node->message[node->neighbors[i]->id][j])) {
	  node->message[node->neighbors[i]->id][j] = 0;
	}
      } else {
	node->message[node->neighbors[i]->id][j] =
	  node->prevmessage[node->neighbors[i]->id][j];
      }
    }
  }
}

int permCheck(node_t **nodes, int n)
     /*checks if messages did not change in the last iterations*/
{
  int i,j,k;
  double change=0;

  for (i=0; i<n; i++) {
    //printMatrix(nodes[i]->message, n);
    for (j=0; j<n; j++) {
      for (k=0; k<n; k++) {
	change+=
	  fabs(nodes[i]->message[j][k] - nodes[i]->prevmessage[j][k]);
      }
    }
  }
  //mexPrintf("%5f, ", change);
  if (change < THRESHOLD)
   return 1;
  return 0;
}

//void readMatrix(double** W, char* filename, int n)
//     /*reads a matrix in from an ascii text file*/
//{
//  int i,j;
//  FILE *fin;
//
//  fin = fopen(filename, "r");
//  for (i=0; i<n; i++) {
//    for (j=0; j<n; j++) {
//      fscanf(fin, "%lf", &W(i,j));
//      W(i,j) = sqrt(W(i,j));
//    }
//  }
//  fclose(fin);
//}

void writeOut(double bethe, int iters, char*filename)
     /*writes a double to ascii file*/
{
  FILE *fout;

  fout = fopen(filename, "wt");
  fprintf(fout, "%f\t%d", bethe, iters);
  fclose(fout);
}

double computeBethe(node_t **X, node_t **Y, int n) {
  /*computes bethe free energy from converged state*/
  double bethe = 0, sumx, sumy, sumxy;
  int i,j,k,l;

  double **bx;
  double **by;
  double**** bxy;
  bx = (double**)malloc(sizeof(double*)*n);
  by = (double**)malloc(sizeof(double*)*n);
  bxy = (double****)malloc(sizeof(double***)*n);
  //init and compute beliefs
  for (i=0; i<n; i++) {
    bxy[i] = (double***)malloc(sizeof(double**)*n);
    bx[i] = (double*)malloc(sizeof(double)*n);
    by[i] = (double*)malloc(sizeof(double)*n);
    for (j=0; j<n; j++) {
      bx[i][j] = X[i]->f[j];
      by[i][j] = Y[i]->f[j];
      bxy[i][j] = (double**)malloc(sizeof(double*)*n);
      for (k=0; k<n; k++) {
	bxy[i][j][k] = (double*)malloc(sizeof(double)*n);
	for (l=0; l<n; l++) {
	  if ((i==l && j!=k) || (i!=l && j==k)) {
	    bxy[i][j][k][l]=0;
	  } else {
	    bxy[i][j][k][l] = //belief b(x_i = k, y_j = l)
	      X[i]->f[k]*Y[j]->f[l]/
	      (X[i]->message[j][l]*Y[j]->message[i][k]);
	  }
	}
      }
    }
  }

  //normalize all beliefs
  for (i=0; i<n; i++) {
    sumx = 0;
    sumy = 0;
    for (j=0; j<n; j++) {
      sumx+=bx[i][j];
      sumy+=by[i][j];
    }
    for (j=0; j<n; j++) {
      bx[i][j]/=sumx;
      by[i][j]/=sumy;
    }
    for (j=0; j<n; j++) {
      sumxy = 0;
      for (k=0; k<n; k++) {
	for (l=0; l<n; l++) {
	  sumxy+=bxy[i][j][k][l];
	}
      }
      for (k=0; k<n; k++) {
	for (l=0; l<n; l++) {
	  bxy[i][j][k][l]/=sumxy;
	}
      }
    }
  }

  //now compute bethe

  for (i=0; i<n; i++) {
    //compute singleton values
    for (j=0; j<n; j++) {
      if (bx[i][j]>THRESHOLD) {
	bethe+=(n-1)*bx[i][j] * safelog(X[i]->phi[j]/bx[i][j]);
      }
      if (by[i][j]>THRESHOLD) {
	bethe+=(n-1)*by[i][j] * safelog(Y[i]->phi[j]/by[i][j]);
      }

      //compute pairwise values


      for (l=0; l<n; l++) {
	for (k=0; k<n; k++) {
	  if (bxy[i][j][k][l]>THRESHOLD) {
	    bethe+=bxy[i][j][k][l] * safelog(bxy[i][j][k][l]/
	    			 (X[i]->phi[k]*Y[j]->phi[l]));
	  }
	}
      }

    }
  }

  sumy = 0;
  for (i=0; i<n; i++) {
    sumx=0;
    for (j=0; j<n; j++) {
      sumx+=bxy[0][0][i][j];
    }
    sumy += fabs(sumx-bx[0][i]);
  }
  //mexPrintf("%f\n", sumy); //this should check if the first pseudo marginals are consistent.

  for (i=0; i<n; i++) {
    free(bx[i]);
    free(by[i]);
    for (j=0; j<n; j++) {
      for (k=0; k<n; k++) {
	free(bxy[i][j][k]);
      }
      free(bxy[i][j]);
    }
    free(bxy[i]);
  }

  return bethe;
}

void mexFunction(int nOut, mxArray *pOut[], int nIn, const mxArray *pIn[]) {
  //assume n objects, nxn matrix

  if (nIn < 1) {
    mexErrMsgIdAndTxt("estper_mex:args", "Usage: [FBethe, B] = estper_mex(A)");
  }

  int i,j,n,counter,converged=CONVERGENCE_INSURANCE,iters=0,links=0;
  double bethe;
  node_t **alpha, **beta;

  Mat<double> W(pIn[0]);
  mxAssert(W.N == W.M, "Input matrix W must be square.");
  n = W.N;

  printMatrix(W, n);

  /*********************************
    allocate and initialize
  *************************************/

  alpha = (node_t**)malloc(n*sizeof(node_t*));
  beta = (node_t**)malloc(n*sizeof(node_t*));

  for (i=0; i<n; i++) { //initialize nodes
    alpha[i] = initNode(n,i,n);
    beta[i] = initNode(n,i,n);
  }

    //mexPrintf("Initialized nodes\n");

    //now connect the nodes

    for (i=0; i<n; i++) {
      counter = 0;
      for (j=0; j<n; j++) {
	if (W(i,j)>NEG_INF) {
	  alpha[i]->neighbors[counter] = beta[j];
	  alpha[i]->phi[counter] = W(i,j);
	  counter++;
	  links++;
	}
      }
    }
    for (i=0; i<n; i++) {
      counter = 0;
      for (j=0; j<n; j++) {
	if (W(j,i)>NEG_INF) {
	  beta[i]->neighbors[counter] = alpha[j];
	  beta[i]->phi[counter] = W(j,i);
	  counter++;
	  links++;
	}
      }
    }

    /*********************************************
			belief propagation
    *********************************************/
    while(converged>0) {
      for (i=0; i<n; i++) {
        saveMessages(alpha[i],n);
	saveMessages(beta[i],n);
      }
      for (i=0; i<n; i++) {
	updateF(alpha[i],n);
	updateF(beta[i],n);
      }
      //mexPrintf("Updated fs\n");
      for (i=0; i<n; i++) {
	updateG(alpha[i],n);
	updateG(beta[i],n);
      }
      //mexPrintf("Updated gs\n");
      for (i=0; i<n; i++) {
	updateMessages(alpha[i],n);
	updateMessages(beta[i],n);
      }

      //check for convergence
      converged -= permCheck(alpha,n);
      //mexPrintf("Checked for convergence == %d\n", converged);

      if (++iters>MAX_ITER) {
	mexPrintf("\nReached maximum iterations without converging. Failing\n");
	//return -1;
	converged=0;
      }
      if (iters>1000 && iters%1000==0) {
	//mexPrintf("Bethe at iteration %d is %f\n", iters, computeBethe(alpha, beta,n));
	damp = damp*.99;
      }

      //mexPrintf("Message x_1->y_2 %f %f %f %f %f\n",
      //     alpha[0]->message[0][0],alpha[0]->message[0][1],
      //     alpha[0]->message[0][2],alpha[0]->message[0][3],
      //     alpha[0]->message[0][4]);
    }
    /********************************************
			output answer
    *********************************************/
    //mexPrintf("\nCompleted after %d iterations\n", iters);

    //bethe = exp(-computeBethe(alpha,beta,n));

    pOut[0] = scalar<double>(computeBethe(alpha,beta,n));
    Mat <double> B(n,n);
    // TODO: FILL UP B
    pOut[1] = B;

    //clean up
    for (i=0; i<n; i++) {
      freeNode(alpha[i],n);
      freeNode(beta[i],n);
    }
    free(alpha);
    free(beta);
}
