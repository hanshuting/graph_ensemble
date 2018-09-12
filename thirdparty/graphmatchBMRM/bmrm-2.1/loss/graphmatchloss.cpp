#ifndef _GRAPHMATCHLOSS_CPP_
#define _GRAPHMATCHLOSS_CPP_

#include "graphmatchloss.hpp"
#include "configuration.hpp"
#include <fstream>
#include <sstream>

#include "math.h"

#include "lap.hpp"
#include "gap.h"

#include "ctime"

#include <map>

//#include <omp.h>

#include "CImg.h"
using namespace cimg_library;

using namespace std;

/**
 * Difference between node features.
 * Here we just use the coordinatewise exponential decay.
 */
inline void Phi1(int* a, int* b, Scalar* res)
{
  for (int i = 0; i < WEIGHTS; i ++) res[i] = -(a[i] - b[i])*(a[i] - b[i]);
}

/**
 * Difference between edge features.
 * Here we just use (a && b).
 */
inline int Phi2(int a, int b)
{
  return a*b;
}

int instance = 0;

/**
 * Find the most violated constraint (if y is nonzero), or perform evaluation.
 *
 * G, Gp -- node features for the two graphs.
 * Gadj, Gpadj -- adjacency matrices for the two graphs.
 * w -- current weight vector.
 * y -- the true assignment (if available).
 * res -- space to store the chosen assignment.
 *
 * weight -- should we perform unweighted or weighted assignment?
 * quad -- should the quadratic term be included?
 */
void assignment(int** G, int** Gp, int** Gadj, int** Gpadj, Scalar* w, int* y, int* res, int weight, int quad, double delta)
{
  double** C;
  Scalar* phi;

  long* row;
  long* col;

  C = new double* [NODES];
  for (int i = 0; i < NODES; i ++) C[i] = new double [NODES];

  phi = new Scalar [WEIGHTS];
  row = new long [NODES];
  col = new long [NODES];

  // Don't bother doing normalisation if we are doing learning...
  if (weight == 1) delta = -1;

  // Randomly permute the correct labelling (since the solver may return an identity matrix on overflow/underflow).
  int* maparray = new int [NODES];
  for (int i = 0; i < NODES; i ++) maparray[i] = i;

  for (int i = 1; i < NODES; i ++)
  {
    int swap = rand() % (i+1);
    int v1 = maparray[i];
    int v2 = maparray[swap];

    maparray[swap] = v1;
    maparray[i] = v2;
  }

  map<int,int> mapf;

  for (int i = 0; i < NODES; i ++) mapf[i] = maparray[i];


  // Compute C_ii' = <Phi1(G_i, G'_i'), w_1>, and add the loss if we are doing column generation.
  // (actually we compute -C, since the lap solver we are using minimises).
  for (int i = 0; i < NODES; i ++)
    for (int ip = 0; ip < NODES; ip ++)
    {
      C[i][ip] = 0;
      Phi1(G[i], Gp[mapf[ip]], phi);
      if (weight == 0) for (int k = 0; k < WEIGHTS; k ++) C[i][ip] -= phi[k];
      else             for (int k = 0; k < WEIGHTS; k ++) C[i][ip] -= w[k]*phi[k];

      if (y != NULL)
        if (y[i] != mapf[ip]) C[i][ip] -= 1.0/NODES;
    }

  // Compute D_(i,i',j,j'),
  // where d_ij = <Phi2(G_i, G'_j), w_2>, C is defined as above.
  if (quad)
  { // Quadratic assignment.
    long double* D = new long double[NODES*NODES*NODES*NODES];

    for (int i = 0; i < NODES; i ++)
      for (int j = 0; j < NODES; j ++)
        {
          int gf = Gadj[i][j];
          for (int ip = 0; ip < NODES; ip ++)
            for (int jp = 0; jp < NODES; jp ++)
              {
                int gpf = Gpadj[mapf[ip]][mapf[jp]];
                if (weight == 0) D[i*NODES*NODES*NODES + ip*NODES*NODES + j*NODES + jp] = Phi2(gf, gpf);
                else             D[i*NODES*NODES*NODES + ip*NODES*NODES + j*NODES + jp] = w[WEIGHTS]*Phi2(gf, gpf);
              }
        }

    // Add the linear term along the diagonal.
    for (int i = 0; i < NODES; i ++)
      for (int ip = 0; ip < NODES; ip ++)
      {
        D[i*NODES*NODES*NODES + ip*NODES*NODES + i*NODES + ip] = -C[i][ip];
      }

//Output to a file if we want to use Matlab...
//    if (delta > 0 && weight == 0)
//     {
//       char* filename = new char [100];
//       sprintf(filename, "matrix%d", instance++);
//       FILE* f = fopen(filename, "w");
//       delete [] filename;
//       for (int i = 0; i < NODES; i ++)
//         for (int j = 0; j < NODES; j ++)
//         {
//           for (int ip = 0; ip < NODES; ip ++)
//             for (int jp = 0; jp < NODES; jp ++)
//               fprintf(f, "%Lf ", D[i*NODES*NODES*NODES + ip*NODES*NODES + j*NODES + jp]);
//           fprintf(f, "\n");
//         }
//       fclose(f);
//     }

    if (delta > 0)
      normalize_matrix(D,NODES,NODES,delta);

    // The solver is quite sensitive to the scale of the terms (it can easily overflow or underflow).
    // Hence we scale the entire matrix.
    Scalar scale = 0;
    for (int i = 0; i < NODES*NODES*NODES*NODES; i ++) scale += fabs(D[i]);
    scale /= NODES*NODES*NODES*NODES;
    for (int i = 0; i < NODES*NODES*NODES*NODES; i ++) D[i] /= 100*scale;

    gap(row, D, NODES, NODES, delta);

    delete [] D;
  }
  else
  { // Linear assignment.

    // The lap solver likes ints rather than longs...
    int* row2 = new int [NODES];
    int* col2 = new int [NODES];

    lap(NODES, C, row2, col2);

    for (int i = 0; i < NODES; i ++)
    {
      row[i] = row2[i];
      col[i] = col2[i];
    }

    delete [] row2;
    delete [] col2;
  }

  // Permute the results back to the original ordering.
  for (int i = 0; i < NODES; i ++)
    res[i] = mapf[row[i]];

  for (int i = 0; i < NODES; i ++) delete [] C[i];
  delete [] C;

  delete [] phi;
  delete [] row;
  delete [] col;
}

/**
 * Compute the feature vector.
 *
 * Phi(G, G', y) = [ sum_ij y_ii'*Phi1(G_i, G'_i'),
                     sum_(i,j,i',j') y_ij*y_i'j'Phi2(G_ij, G'_i'j') ]
 */
void CGraphMatchLoss::Phi(int n, int* y, Scalar* res)
{
  Scalar* temp = new Scalar [WEIGHTS];

  for (int i = 0; i < WEIGHTS + quad; i ++) res[i] = 0;

  for (int i = 0; i < NODES; i ++)
    for (int ip = 0; ip < NODES; ip ++)
    {
      Phi1(_data->_G[n][i], _data->_Gp[n][ip], temp);
      if (y[i] == ip)
        for (int k = 0; k < WEIGHTS; k ++) res[k] += temp[k];

      if (quad)
      { // Include the quadratic term.
        if (y[i] == ip)
          for (int j = 0; j < NODES; j ++)
            for (int jp = 0; jp < NODES; jp ++)
            {
              int gf = _data->_Gadj[n]->mat[i][j];
              int gpf = _data->_Gpadj[n]->mat[ip][jp];
              if (y[j] == jp) res[WEIGHTS] += Phi2(gf, gpf);
            }
      }
    }

  delete [] temp;
}

/** Constructor. */
CGraphMatchLoss::CGraphMatchLoss(CModel* &model, CGraphData* &data) : CLoss(model, 1, WEIGHTS + Configuration::GetInstance().GetInt("quadratic"), data->bias()), _data(data)
{
  Configuration &config = Configuration::GetInstance();

  quad = config.GetInt("quadratic");
  verbosity = config.GetInt("Loss.verbosity");
}

/** Destructor. */
CGraphMatchLoss::~CGraphMatchLoss() {}

void CGraphMatchLoss::Usage() {}

void CGraphMatchLoss::ComputeLossAndGradient(Scalar& loss, TheMatrix& grad)
{
  TheMatrix &w = _model->GetW();
  loss = 0;
  grad.Zero();


  Scalar* dat = w.Data();
  Scalar* raw_g = grad.Data();

  Configuration &config = Configuration::GetInstance();

  Scalar lambda = config.GetDouble("BMRM.lambda");
  Scalar delta = config.GetDouble("delta");

  for (int i = 0; i < WEIGHTS + quad; i ++) raw_g[i] = 0;//lambda*dat[i];

  Scalar* lossi = new Scalar [_data->_N];
  Scalar** raw_gi = new Scalar* [_data->_N];

  int i;
//  #pragma omp parallel for
  for (i = 0; i < _data->_N; i ++)
  {
    Scalar* resy;
    Scalar* resybar;

    int* ybar = new int [NODES];

    resy = new Scalar [WEIGHTS + quad];
    resybar = new Scalar [WEIGHTS + quad];

    int** mat = NULL;
    int** matp = NULL;
    if (quad)
    {
      mat = _data->_Gadj[i]->mat;
      matp = _data->_Gpadj[i]->mat;
    }
    assignment(_data->_G[i], _data->_Gp[i], mat, matp, dat, _data->_Y[i], ybar, 1, quad, delta);
    Phi(i, _data->_Y[i], resy);
    Phi(i, ybar, resybar);

    Scalar inp = 0;
    for (int j = 0; j < WEIGHTS + quad; j ++) inp += (resybar[j]-resy[j])*dat[j];

    Scalar labloss = LabelLoss(_data->_Y[i], ybar);
    //if (inp + labloss > 0)
      lossi[i] = labloss;

    raw_gi[i] = new Scalar [WEIGHTS + quad];
    for (int j = 0; j < WEIGHTS + quad; j ++)
    {
      //if (inp + labloss > 0)
      {
        lossi[i] += dat[j]*(resybar[j]-resy[j]);
        raw_gi[i][j] = (1.0/_data->_N)*(resybar[j]-resy[j]);
      }
    }

    delete [] ybar;

    delete [] resy;
    delete [] resybar;
  }

  for (int j = 0; j < _data->_N; j ++)
  {
    loss += lossi[j];
    for (int k = 0; k < WEIGHTS + quad; k ++)
      raw_g[k] += raw_gi[j][k];

    delete [] raw_gi[j];
  }
  delete [] raw_gi;
  delete [] lossi;

  loss = loss/_data->_N;
}

void CGraphMatchLoss::Predict(CModel *model)
{
  Evaluate(model);
}

Scalar mean(Scalar* data, int n)
{
  Scalar total = 0;
  for (int i = 0; i < n; i ++) total += data[i];
  return total/n;
}

Scalar sde(Scalar* data, int n)
{
  Scalar avg = 0;
  Scalar avg2 = 0;

  for (int i = 0; i < n; i ++)
  {
    avg += data[i];
    avg2 += data[i]*data[i];
  }
  avg /= n;
  avg2 /= n;

  Scalar sd = sqrt(avg2 - avg*avg);
  return sd / sqrt((float) n);
}

void savequad(char* filename, char* im, pair<Scalar,Scalar>* corn, int** adjmat)
{
  CImg<unsigned char> image_bw(im);

  CImg<unsigned char> image_out(image_bw.width(), image_bw.height(), 1, 3, 0);
  const unsigned char blue[] = { 0,0,255 };

  for (int y = 0; y < image_out.height(); y ++)
    for (int x = 0; x < image_out.width(); x ++)
      for (int c = 0; c < 3; c ++)
        image_out(x,y,c) = image_bw(x,y);

  for (int i = 0; i < NODES; i ++)
    for (int j = 0; j < NODES; j ++)
    {
      if (adjmat[i][j])
      {
        int p1x = (int) corn[i].first;
        int p1y = (int) corn[i].second;
        int p2x = (int) corn[j].first;
        int p2y = (int) corn[j].second;
        image_out.draw_line(p1x,p1y,p2x,p2y,blue);
      }
    }
  image_out.save(filename);
}

void saveimage(char* filename, char* im1, char* im2, pair<Scalar,Scalar>* corn1, pair<Scalar,Scalar>* corn2, int* match, int* correct)
{

  CImg<unsigned char> image_lbw(im1);
  CImg<unsigned char> image_rbw(im2);

  CImg<unsigned char> image_l(image_lbw.width(), image_lbw.height(), 1, 3, 0);
  CImg<unsigned char> image_r(image_rbw.width(), image_rbw.height(), 1, 3, 0);

  for (int y = 0; y < image_l.height(); y ++)
    for (int x = 0; x < image_l.width(); x ++)
      for (int c = 0; c < 3; c ++)
      {
        image_l(x,y,c) = image_lbw(x,y);
        image_r(x,y,c) = image_rbw(x,y);
      }

  CImg<unsigned char> image_out(2*image_l.width(), image_l.height(), 1, 3, 0);
  const unsigned char red[] = { 255,0,0 };
  const unsigned char blue[] = { 0,0,255 };

  for (int i = 0; i < NODES; i ++)
  {
    int plx = (int) corn1[i].first;
    int ply = (int) corn1[i].second;
    int prx = (int) corn2[i].first;
    int pry = (int) corn2[i].second;

    image_l.draw_triangle(plx, ply-2, plx+2, ply+2, plx-2, ply+2,blue);
    image_r.draw_triangle(prx, pry-2, prx+2, pry+2, prx-2, pry+2,blue);
  }
  for (int y = 0; y < image_l.height(); y ++)
    for (int x = 0; x < image_l.width(); x ++)
      for (int c = 0; c < 3; c ++)
      {
        image_out(x,y,c) = image_l(x,y,c);
        image_out(x+image_l.width(),y,c) = image_r(x,y,c);
      }

  for (int i = 0; i < NODES; i ++)
  {
    int plx = (int) corn1[i].first;
    int ply = (int) corn1[i].second;
    int prx = image_l.width() + (int) corn2[match[i]].first;
    int pry = (int) corn2[match[i]].second;

    if (match[i] == correct[i])
      image_out.draw_line(plx,ply,prx,pry,blue);
    else
      image_out.draw_line(plx,ply,prx,pry,red);
  }
  image_out.save(filename);
}

/**
 * Compare the performance of the model after learning to the model before learning.
 * The performance before learning is found using a constant weight vector.
 */
void CGraphMatchLoss::Evaluate(CModel *model)
{
  TheMatrix &w = _model->GetW();
  Scalar loss_weight = 0;
  Scalar loss_noweight = 0;
  Scalar* time_weight = new Scalar[_data->_N];
  Scalar* time_noweight = new Scalar[_data->_N];

  Scalar* dat = w.Data();
  int* ybar;

  Scalar ll = 0;
  Configuration &config = Configuration::GetInstance();
  Scalar delta = config.GetDouble("delta");

  ybar = new int [NODES];

  Scalar* errorw = new Scalar[_data->_N];
  Scalar* errorn = new Scalar[_data->_N];

  for (int i = 0; i < _data->_N; i ++)
  {
    int** mat = NULL;
    int** matp = NULL;
    if (quad)
    {
      mat = _data->_Gadj[i]->mat;
      matp = _data->_Gpadj[i]->mat;
    }
    clock_t start = clock();
    assignment(_data->_G[i], _data->_Gp[i], mat, matp, dat, NULL, ybar, 1, quad, delta);
    time_weight[i] = (clock() - start)/ (double) CLOCKS_PER_SEC;

    ll = LabelLoss(_data->_Y[i], ybar); loss_weight += ll;
    errorw[i] = ll;

    // Print the performance for this instance.
    //printf("Weight: %f\n", ll);

    // Save the matching images in the folder `matches'.
    //int a1 = _data->pairs[i].first;
    //int a2 = _data->pairs[i].second;
    //char* save1_loc = new char [100];
    //sprintf(save1_loc, "matches/weight_%d_%d_%d.jpg", quad, a1, a2);
    //saveimage(save1_loc, _data->imnames[a1], _data->imnames[a2], _data->corners[a1], _data->corners[a2], ybar, _data->_Y[i]);
    //delete [] save1_loc;

    start = clock();
    assignment(_data->_G[i], _data->_Gp[i], mat, matp, dat, NULL, ybar, 0, quad, delta);

    time_noweight[i] = (1.0*clock() - start)/ (double) CLOCKS_PER_SEC;

    ll = LabelLoss(_data->_Y[i], ybar); loss_noweight += ll;
    errorn[i] = ll;

    // Print the performance for this instance.
    //printf("No weight: %f\n", ll);

    // Save the matching images in the folder `matches'.
    //char* save2_loc = new char [100];
    //sprintf(save2_loc, "matches/noweight_%d_%d_%d.jpg", quad, a1, a2);
    //saveimage(save2_loc, _data->imnames[a1], _data->imnames[a2], _data->corners[a1], _data->corners[a2], ybar, _data->_Y[i]);
    //delete [] save2_loc;

    if (quad)
    {
      // Save the adjacency images in the folder `matches'.
      //char* adj_loc = new char [100];
      //sprintf(adj_loc, "matches/adj_%d.jpg", a1);
      //savequad(adj_loc, _data->imnames[a1], _data->corners[a1], _data->_Gadj[i]->mat);
      //delete [] adj_loc;
    }
  }

  loss_weight = loss_weight/_data->_N;
  loss_noweight = loss_noweight/_data->_N;

  cout << "Error with weights: " << loss_weight << endl;
  cout << "Error without weights: " << loss_noweight << endl;

// Output the results to a temporary file.
//   FILE* res_temp = fopen("res_temp.txt", "w");
//   fprintf(res_temp, "%f %f %f %f %f %f %f %f\n", loss_weight, loss_noweight, sde(errorw, _data->_N), sde(errorn, _data->_N), mean(time_weight, _data->_N), mean(time_noweight, _data->_N), sde(time_weight, _data->_N), sde(time_noweight, _data->_N));
//   fclose(res_temp);

  delete [] errorw;
  delete [] errorn;

  delete [] time_weight;
  delete [] time_noweight;
}

/** Loss function. */
Scalar CGraphMatchLoss::LabelLoss(int* y, int* ybar)
{
   // y = actual, ybar = predicted.
   float mat = 0;
   for (int i = 0; i < NODES; i ++)
     if (y[i] == ybar[i]) mat ++;

   return 1 - mat/NODES;
}

#endif

