/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    NonLocalMeansDenoise.cxx
  Language:  C++
  Website:   itksnap.org/c3d
  Copyright (c) 2014 Paul A. Yushkevich
  
  This file is part of C3D, a command-line companion tool to ITK-SNAP

  C3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

// THIS CODE IS A DERIVATIVE OF THE ORIGINAL CODE BY Jose V. Manjon and Pierrick Coupe
// WITH MODIFICATIONS BY Long Xie

/* Jose V. Manjon - jmanjon@fis.upv.es                                     */
/* Pierrick Coupe - pierrick.coupe@gmail.com                               */
/* Brain Imaging Center, Montreal Neurological Institute.                  */
/* Mc Gill University                                                      */
/*                                                                         */
/* Copyright (C) 2010 Jose V. Manjon and Pierrick Coupe                    */

/***************************************************************************
*              Adaptive Non-Local Means Denoising of MR Images            *
%              With Spatially Varying Noise Levels                        *
%                                                                         *
%              Jose V. Manjon, Pierrick Coupe, Luis Marti-Bonmati,        *
%              D. Louis Collins and Montserrat Robles                      *
***************************************************************************/

#include "NonLocalMeansDenoise.h"
#include <iostream>
#include <string>
#include "math.h"
#include <cstdlib>
#include <itkMultiThreaderBase.h>


#define pi 3.1415926535

using namespace std;


struct DenoiseParameters
{
  int v;
  int f;
  bool rician;

  DenoiseParameters() :
    v(3), f(1), rician(false){}
};

template <class TFloat>
struct ThreadArgument
{
  int rows;
  int cols;
  int slices;
  TFloat * in_image;
  TFloat * means_image;
  TFloat * var_image;
  TFloat * estimate;
  TFloat * label;
  TFloat * bias;
  int ini;
  int fin;
  int radioB;
  int radioS;
  bool rician;
  TFloat max;
};

template <class TFloat>
class NLMDenoiseProblem
{
public:
  typedef itk::OrientedRASImage<TFloat, 3> ImageType;
  typedef typename ImageType::Pointer      ImagePointer;

  /*Returns the modified Bessel function I0(x) for any real x.*/
  static TFloat
  bessi0(TFloat x);

  /*Returns the modified Bessel function I1(x) for any real x.*/
  static TFloat
  bessi1(TFloat x);

  static TFloat
  Epsi(TFloat snr);

  /* Function which compute the weighted average for one block */
  static void
  Average_block(TFloat * ima,
                int      x,
                int      y,
                int      z,
                int      neighborhoodsize,
                TFloat * average,
                TFloat   weight,
                int      sx,
                int      sy,
                int      sz,
                bool     rician);

  /* Function which computes the value assigned to each voxel */
  static void
  Value_block(TFloat * Estimate,
              TFloat * Label,
              int      x,
              int      y,
              int      z,
              int      neighborhoodsize,
              TFloat * average,
              TFloat   global_sum,
              int      sx,
              int      sy,
              int      sz);

  static TFloat
  distance(TFloat * ima, int x, int y, int z, int nx, int ny, int nz, int f, int sx, int sy, int sz);

  static TFloat
  distance2(TFloat * ima, TFloat * medias, int x, int y, int z, int nx, int ny, int nz, int f, int sx, int sy, int sz);

  static void
  Regularize(TFloat * in, TFloat * out, int r, int sx, int sy, int sz);

  static void *
  ThreadFunc(void * pArguments);

  static ImagePointer
  NLMdenoise(ImageType * image, DenoiseParameters & param, std::ostream & sout = std::cout);

private:
};

//bool rician=false;
//double max;

/*Returns the modified Bessel function I0(x) for any real x.*/
template <class TFloat>
TFloat
NLMDenoiseProblem<TFloat>::bessi0(TFloat x)
{
  TFloat ax, ans, a;
  TFloat y;
  if ((ax = fabs(x)) < 3.75)
  {
    y = x / 3.75;
    y *= y;
    ans =
      1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
  }
  else
  {
    y = 3.75 / ax;
    ans = (exp(ax) / sqrt(ax));
    a = y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))));
    ans = ans * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + a))));
  }
  return ans;
}

/*Returns the modified Bessel function I1(x) for any real x.*/
template <class TFloat>
TFloat
NLMDenoiseProblem<TFloat>::bessi1(TFloat x)
{
  TFloat ax, ans;
  TFloat y;
  if ((ax = fabs(x)) < 3.75)
  {
    y = x / 3.75;
    y *= y;
    ans =
      ax * (0.5 + y * (0.87890594 +
                       y * (0.51498869 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
  }
  else
  {
    y = 3.75 / ax;
    ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
    ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
    ans *= (exp(ax) / sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}

template <class TFloat>
TFloat
NLMDenoiseProblem<TFloat>::Epsi(TFloat snr)
{
  TFloat val;
  val = 2 + snr * snr -
        (pi / 8) * exp(-(snr * snr) / 2) *
          ((2 + snr * snr) * bessi0((snr * snr) / 4) + (snr * snr) * bessi1((snr * snr) / 4)) *
          ((2 + snr * snr) * bessi0((snr * snr) / 4) + (snr * snr) * bessi1((snr * snr) / 4));
  if (val < 0.001)
    val = 1;
  if (val > 10)
    val = 1;
  return val;
}

/* Function which compute the weighted average for one block */
template <class TFloat>
void
NLMDenoiseProblem<TFloat>::Average_block(TFloat * ima,
                                         int      x,
                                         int      y,
                                         int      z,
                                         int      neighborhoodsize,
                                         TFloat * average,
                                         TFloat   weight,
                                         int      sx,
                                         int      sy,
                                         int      sz,
                                         bool     rician)
{
  int  x_pos, y_pos, z_pos;
  bool is_outside;
  int  a, b, c, ns, sxy, count;

  // extern bool rician;

  ns = 2 * neighborhoodsize + 1;
  sxy = sx * sy;

  count = 0;

  for (c = 0; c < ns; c++)
  {
    for (b = 0; b < ns; b++)
    {
      for (a = 0; a < ns; a++)
      {
        is_outside = false;
        x_pos = x + a - neighborhoodsize;
        y_pos = y + b - neighborhoodsize;
        z_pos = z + c - neighborhoodsize;

        if ((z_pos < 0) || (z_pos > sz - 1))
          is_outside = true;
        if ((y_pos < 0) || (y_pos > sy - 1))
          is_outside = true;
        if ((x_pos < 0) || (x_pos > sx - 1))
          is_outside = true;

        if (rician)
        {
          if (is_outside)
            average[count] = average[count] + ima[z * (sxy) + (y * sx) + x] * ima[z * (sxy) + (y * sx) + x] * weight;
          else
            average[count] = average[count] + ima[z_pos * (sxy) + (y_pos * sx) + x_pos] *
                                                ima[z_pos * (sxy) + (y_pos * sx) + x_pos] * weight;
        }
        else
        {
          if (is_outside)
            average[count] = average[count] + ima[z * (sxy) + (y * sx) + x] * weight;
          else
            average[count] = average[count] + ima[z_pos * (sxy) + (y_pos * sx) + x_pos] * weight;
        }
        count++;
      }
    }
  }
}

/* Function which computes the value assigned to each voxel */
template <class TFloat>
void
NLMDenoiseProblem<TFloat>::Value_block(TFloat * Estimate,
                                       TFloat * Label,
                                       int      x,
                                       int      y,
                                       int      z,
                                       int      neighborhoodsize,
                                       TFloat * average,
                                       TFloat   global_sum,
                                       int      sx,
                                       int      sy,
                                       int      sz)
{
  int    x_pos, y_pos, z_pos;
  bool   is_outside;
  TFloat value = 0.0;
  TFloat label = 0.0;
  TFloat denoised_value = 0.0;
  int    count = 0;
  int    a, b, c, ns, sxy;

  ns = 2 * neighborhoodsize + 1;
  sxy = sx * sy;


  for (c = 0; c < ns; c++)
  {
    for (b = 0; b < ns; b++)
    {
      for (a = 0; a < ns; a++)
      {
        is_outside = false;
        x_pos = x + a - neighborhoodsize;
        y_pos = y + b - neighborhoodsize;
        z_pos = z + c - neighborhoodsize;

        if ((z_pos < 0) || (z_pos > sz - 1))
          is_outside = true;
        if ((y_pos < 0) || (y_pos > sy - 1))
          is_outside = true;
        if ((x_pos < 0) || (x_pos > sx - 1))
          is_outside = true;
        if (!is_outside)
        {
          value = Estimate[z_pos * (sxy) + (y_pos * sx) + x_pos];
          value = value + (average[count] / global_sum);

          label = Label[(x_pos + y_pos * sx + z_pos * sxy)];
          Estimate[z_pos * (sxy) + (y_pos * sx) + x_pos] = value;
          Label[(x_pos + y_pos * sx + z_pos * sxy)] = label + 1;
        }
        count++;
      }
    }
  }
}

template <class TFloat>
TFloat
NLMDenoiseProblem<
  TFloat>::distance(TFloat * ima, int x, int y, int z, int nx, int ny, int nz, int f, int sx, int sy, int sz)
{
  TFloat d, acu, distancetotal;
  int    i, j, k, ni1, nj1, ni2, nj2, nk1, nk2;

  distancetotal = 0;
  for (k = -f; k <= f; k++)
  {
    nk1 = z + k;
    nk2 = nz + k;
    if (nk1 < 0)
      nk1 = -nk1;
    if (nk2 < 0)
      nk2 = -nk2;
    if (nk1 >= sz)
      nk1 = 2 * sz - nk1 - 1;
    if (nk2 >= sz)
      nk2 = 2 * sz - nk2 - 1;

    for (j = -f; j <= f; j++)
    {
      nj1 = y + j;
      nj2 = ny + j;
      if (nj1 < 0)
        nj1 = -nj1;
      if (nj2 < 0)
        nj2 = -nj2;
      if (nj1 >= sy)
        nj1 = 2 * sy - nj1 - 1;
      if (nj2 >= sy)
        nj2 = 2 * sy - nj2 - 1;

      for (i = -f; i <= f; i++)
      {
        ni1 = x + i;
        ni2 = nx + i;
        if (ni1 < 0)
          ni1 = -ni1;
        if (ni2 < 0)
          ni2 = -ni2;
        if (ni1 >= sx)
          ni1 = 2 * sx - ni1 - 1;
        if (ni2 >= sx)
          ni2 = 2 * sx - ni2 - 1;

        distancetotal =
          distancetotal + ((ima[nk1 * (sx * sy) + (nj1 * sx) + ni1] - ima[nk2 * (sx * sy) + (nj2 * sx) + ni2]) *
                           (ima[nk1 * (sx * sy) + (nj1 * sx) + ni1] - ima[nk2 * (sx * sy) + (nj2 * sx) + ni2]));
      }
    }
  }

  acu = (2 * f + 1) * (2 * f + 1) * (2 * f + 1);
  d = distancetotal / acu;

  return d;
}

template <class TFloat>
TFloat
NLMDenoiseProblem<TFloat>::distance2(TFloat * ima,
                                     TFloat * medias,
                                     int      x,
                                     int      y,
                                     int      z,
                                     int      nx,
                                     int      ny,
                                     int      nz,
                                     int      f,
                                     int      sx,
                                     int      sy,
                                     int      sz)
{
  TFloat d, acu, distancetotal;
  int    i, j, k, ni1, nj1, ni2, nj2, nk1, nk2;

  acu = 0;
  distancetotal = 0;

  for (k = -f; k <= f; k++)
  {
    nk1 = z + k;
    nk2 = nz + k;
    if (nk1 < 0)
      nk1 = -nk1;
    if (nk2 < 0)
      nk2 = -nk2;
    if (nk1 >= sz)
      nk1 = 2 * sz - nk1 - 1;
    if (nk2 >= sz)
      nk2 = 2 * sz - nk2 - 1;
    for (j = -f; j <= f; j++)
    {
      nj1 = y + j;
      nj2 = ny + j;
      if (nj1 < 0)
        nj1 = -nj1;
      if (nj2 < 0)
        nj2 = -nj2;
      if (nj1 >= sy)
        nj1 = 2 * sy - nj1 - 1;
      if (nj2 >= sy)
        nj2 = 2 * sy - nj2 - 1;
      for (i = -f; i <= f; i++)
      {
        ni1 = x + i;
        ni2 = nx + i;
        if (ni1 < 0)
          ni1 = -ni1;
        if (ni2 < 0)
          ni2 = -ni2;
        if (ni1 >= sx)
          ni1 = 2 * sx - ni1 - 1;
        if (ni2 >= sx)
          ni2 = 2 * sx - ni2 - 1;

        d = (ima[nk1 * (sx * sy) + (nj1 * sx) + ni1] - medias[nk1 * (sx * sy) + (nj1 * sx) + ni1]) -
            (ima[nk2 * (sx * sy) + (nj2 * sx) + ni2] - medias[nk2 * (sx * sy) + (nj2 * sx) + ni2]);
        distancetotal = distancetotal + d * d;
      }
    }
  }

  acu = (2 * f + 1) * (2 * f + 1) * (2 * f + 1);
  d = distancetotal / acu;

  return d;
}

template <class TFloat>
void
NLMDenoiseProblem<TFloat>::Regularize(TFloat * in, TFloat * out, int r, int sx, int sy, int sz)
{
  TFloat acu;
  int    ind, i, j, k, ni, nj, nk, ii, jj, kk;

  // TFloat * temp=mxCalloc(sx*sy*sz,sizeof(TFloat));
  TFloat * temp = new TFloat[sx * sy * sz];
  // TFloat *temp = (TFloat*)malloc(sx*sy*sz*sizeof(TFloat));


  /* separable convolution */
  for (k = 0; k < sz; k++)
    for (j = 0; j < sy; j++)
      for (i = 0; i < sx; i++)
      {
        if (in[k * (sx * sy) + (j * sx) + i] == 0)
          continue;

        acu = 0;
        ind = 0;
        for (ii = -r; ii <= r; ii++)
        {
          ni = i + ii;
          if (ni < 0)
            ni = -ni;
          if (ni >= sx)
            ni = 2 * sx - ni - 1;
          if (in[k * (sx * sy) + (j * sx) + ni] > 0)
          {
            acu += in[k * (sx * sy) + (j * sx) + ni];
            ind++;
          }
        }
        if (ind == 0)
          ind = 1;
        out[k * (sx * sy) + (j * sx) + i] = acu / ind;
      }
  for (k = 0; k < sz; k++)
    for (j = 0; j < sy; j++)
      for (i = 0; i < sx; i++)
      {
        if (out[k * (sx * sy) + (j * sx) + i] == 0)
          continue;

        acu = 0;
        ind = 0;
        for (jj = -r; jj <= r; jj++)
        {
          nj = j + jj;
          if (nj < 0)
            nj = -nj;
          if (nj >= sy)
            nj = 2 * sy - nj - 1;

          if (out[k * (sx * sy) + (nj * sx) + i] > 0)
          {
            acu += out[k * (sx * sy) + (nj * sx) + i];
            ind++;
          }
        }
        if (ind == 0)
          ind = 1;
        temp[k * (sx * sy) + (j * sx) + i] = acu / ind;
      }
  for (k = 0; k < sz; k++)
    for (j = 0; j < sy; j++)
      for (i = 0; i < sx; i++)
      {
        if (temp[k * (sx * sy) + (j * sx) + i] == 0)
          continue;

        acu = 0;
        ind = 0;
        for (kk = -r; kk <= r; kk++)
        {
          nk = k + kk;
          if (nk < 0)
            nk = -nk;
          if (nk >= sz)
            nk = 2 * sz - nk - 1;
          if (temp[nk * (sx * sy) + (j * sx) + i] > 0)
          {
            acu += temp[nk * (sx * sy) + (j * sx) + i];
            ind++;
          }
        }
        if (ind == 0)
          ind = 1;
        out[k * (sx * sy) + (j * sx) + i] = acu / ind;
      }
  // mxFree(temp);
  delete [] temp;
  // free(temp);

  return;
}

template <class TFloat>
void *
NLMDenoiseProblem<TFloat>::ThreadFunc(void * pArguments)
{
  TFloat *bias, *Estimate, *Label, *ima, *means, *variances, *average, epsilon, mu1, var1, totalweight, wmax, t1, t1i,
    t2, d, w, distanciaminima, max;
  int  rows, cols, slices, ini, fin, v, f, init, i, j, k, rc, ii, jj, kk, ni, nj, nk, Ndims;
  bool rician;

  // extern bool rician;
  // extern TFloat max;

  ThreadArgument<TFloat> arg;
  arg = *(ThreadArgument<TFloat> *) pArguments;

  rows = arg.rows;
  cols = arg.cols;
  slices = arg.slices;
  ini = arg.ini;
  fin = arg.fin;
  ima = arg.in_image;
  means = arg.means_image;
  variances = arg.var_image;
  Estimate = arg.estimate;
  bias = arg.bias;
  Label = arg.label;
  v = arg.radioB;
  f = arg.radioS;
  rician = arg.rician;
  max = arg.max;

  /* filter */
  epsilon = 0.00001;
  mu1 = 0.95;
  var1 = 0.5;
  init = 0;
  rc = rows * cols;

  Ndims = (2 * f + 1) * (2 * f + 1) * (2 * f + 1);

  average = (TFloat *)malloc(Ndims * sizeof(TFloat));

  wmax = 0.0;

  for (k = ini; k < fin; k += 2)
    for (j = 0; j < rows; j += 2)
      for (i = 0; i < cols; i += 2)
      {
        /* init  */
        for (init = 0; init < Ndims; init++)
          average[init] = 0.0;
        totalweight = 0.0;
        distanciaminima = 100000000000000;

        if (ima[k * rc + (j * cols) + i] > 0 && (means[k * rc + (j * cols) + i]) > epsilon &&
            (variances[k * rc + (j * cols) + i] > epsilon))
        {
          /* calculate minimum distance */
          for (kk = -v; kk <= v; kk++)
          {
            nk = k + kk;
            for (jj = -v; jj <= v; jj++)
            {
              nj = j + jj;
              for (ii = -v; ii <= v; ii++)
              {
                ni = i + ii;
                if (ii == 0 && jj == 0 && kk == 0)
                  continue;
                if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
                {
                  if (ima[nk * rc + (nj * cols) + ni] > 0 && (means[nk * (rc) + (nj * cols) + ni]) > epsilon &&
                      (variances[nk * rc + (nj * cols) + ni] > epsilon))
                  {
                    t1 = (means[k * rc + (j * cols) + i]) / (means[nk * rc + (nj * cols) + ni]);
                    t1i = (max - means[k * (rc) + (j * cols) + i]) / (max - means[nk * (rc) + (nj * cols) + ni]);
                    t2 = (variances[k * rc + (j * cols) + i]) / (variances[nk * rc + (nj * cols) + ni]);

                    if (((t1 > mu1 && t1 < (1 / mu1)) || (t1i > mu1 && t1i < (1 / mu1))) && t2 > var1 &&
                        t2 < (1 / var1))
                    {
                      d = distance2(ima, means, i, j, k, ni, nj, nk, f, cols, rows, slices);
                      if (d < distanciaminima)
                        distanciaminima = d;
                    }
                  }
                }
              }
            }
          }
          if (distanciaminima == 0)
            distanciaminima = 1;

          /* rician correction */
          if (rician)
          {
            for (kk = -f; kk <= f; kk++)
            {
              nk = k + kk;
              for (ii = -f; ii <= f; ii++)
              {
                ni = i + ii;
                for (jj = -f; jj <= f; jj++)
                {
                  nj = j + jj;
                  if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
                  {
                    if (distanciaminima == 100000000000000)
                      bias[nk * (rc) + (nj * cols) + ni] = 0;
                    else
                      bias[nk * (rc) + (nj * cols) + ni] = (distanciaminima);
                  }
                }
              }
            }
          }

          /* block filtering */
          for (kk = -v; kk <= v; kk++)
          {
            nk = k + kk;
            for (jj = -v; jj <= v; jj++)
            {
              nj = j + jj;
              for (ii = -v; ii <= v; ii++)
              {
                ni = i + ii;
                if (ii == 0 && jj == 0 && kk == 0)
                  continue;

                if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
                {
                  if (ima[nk * rc + (nj * cols) + ni] > 0 && (means[nk * (rc) + (nj * cols) + ni]) > epsilon &&
                      (variances[nk * rc + (nj * cols) + ni] > epsilon))
                  {
                    t1 = (means[k * rc + (j * cols) + i]) / (means[nk * rc + (nj * cols) + ni]);
                    t1i = (max - means[k * (rc) + (j * cols) + i]) / (max - means[nk * (rc) + (nj * cols) + ni]);
                    t2 = (variances[k * rc + (j * cols) + i]) / (variances[nk * rc + (nj * cols) + ni]);

                    if (((t1 > mu1 && t1 < (1 / mu1)) || (t1i > mu1 && t1i < (1 / mu1))) && t2 > var1 &&
                        t2 < (1 / var1))
                    {
                      d = distance(ima, i, j, k, ni, nj, nk, f, cols, rows, slices);

                      if (d > 3 * distanciaminima)
                        w = 0;
                      else
                        w = exp(-d / distanciaminima);

                      if (w > wmax)
                        wmax = w;

                      if (w > 0)
                      {
                        Average_block(ima, ni, nj, nk, f, average, w, cols, rows, slices, rician);
                        totalweight = totalweight + w;
                      }
                    }
                  }
                }
              }
            }
          }

          if (wmax == 0.0)
            wmax = 1.0;
          Average_block(ima, i, j, k, f, average, wmax, cols, rows, slices, rician);
          totalweight = totalweight + wmax;
          Value_block(Estimate, Label, i, j, k, f, average, totalweight, cols, rows, slices);
        }
        else
        {
          wmax = 1.0;
          Average_block(ima, i, j, k, f, average, wmax, cols, rows, slices, rician);
          totalweight = totalweight + wmax;
          Value_block(Estimate, Label, i, j, k, f, average, totalweight, cols, rows, slices);
        }
      }

  return 0;
}

//TFloat* NLMdenoise(TFloat *ima, int v, int f, bool rician, int ndim, const int *dims)
template <class TFloat>
typename NLMDenoiseProblem<TFloat>::ImagePointer
  NLMDenoiseProblem<TFloat>
  ::NLMdenoise(ImageType *input, DenoiseParameters &param, std::ostream &sout)
{
  // dimension
  const unsigned int ndim = 3;

  // read inputs
  int v = param.v;
  int f = param.f;
  bool rician = param.rician;

  // read image
  ImagePointer ITKima = input;
  TFloat *ima = ITKima->GetBufferPointer();

  // get image size
  itk::Size<ndim> dims = ITKima->GetBufferedRegion().GetSize();
  int nelement = dims[2] * dims[1] * dims[0];

  // allocate memory for output image
  ImagePointer ITKfima = ImageType::New();
  ITKfima->SetRegions(ITKima->GetBufferedRegion());
  ITKfima->SetSpacing(ITKima->GetSpacing());
  ITKfima->SetOrigin(ITKima->GetOrigin());
  ITKfima->SetDirection(ITKima->GetDirection());
  ITKfima->SetMetaDataDictionary(ITKima->GetMetaDataDictionary());
  ITKfima->Allocate();
  TFloat *fima = ITKfima->GetBufferPointer();

  // allocate memory for intermediate matrices
  TFloat *means = new TFloat [nelement];
  TFloat *variances = new TFloat [nelement];
  TFloat *Estimate = new TFloat [nelement];
  TFloat *Label = new TFloat [nelement];
  TFloat *bias;
  if (rician) bias = new TFloat [nelement];

  // declare temporary variables
  TFloat SNR,mean,var,label,estimate;
  int i,j,k,ii,jj,kk,ni,nj,nk,indice,Nthreads,ini,fin,r;

  int Ndims = pow((2*f+1),ndim);
  TFloat *average = new TFloat [Ndims];

  // initialize all to 0
  for (i = 0; i < dims[2] *dims[1] * dims[0];i++)
  {
    Estimate[i] = 0.0;
    Label[i] = 0.0;
    fima[i] = 0.0;
    if(rician) bias[i]=0.0;
  }

  // get max of the image and image with the mean of neighborhood 3x3x3 patch
  TFloat max=0;
  for(k=0;k<dims[2];k++)
  {
    for(j=0;j<dims[1];j++)
    {
      for(i=0;i<dims[0];i++)
      {
        // get the maximum of the image array
        if(ima[k*(dims[0]*dims[1])+(j*dims[0])+i]>max)
          max=ima[k*(dims[0]*dims[1])+(j*dims[0])+i];

        // get the mean of the neighborhood 3x3x3 patch
        mean=0;
        indice=0;
        for(ii=-1;ii<=1;ii++)
        {
          for(jj=-1;jj<=1;jj++)
          {
            for(kk=-1;kk<=1;kk++)
            {
              ni=i+ii;
              nj=j+jj;
              nk=k+kk;

              if(ni<0) ni=-ni;
              if(nj<0) nj=-nj;
              if(nk<0) nk=-nk;
              if(ni>=dims[0]) ni=2*dims[0]-ni-1;
              if(nj>=dims[1]) nj=2*dims[1]-nj-1;
              if(nk>=dims[2]) nk=2*dims[2]-nk-1;

              mean = mean + ima[nk*(dims[0]*dims[1])+(nj*dims[0])+ni];
              indice=indice+1;

            }
          }
        }
        mean=mean/indice;
        means[k*(dims[0]*dims[1])+(j*dims[0])+i]=mean;
      }
    }
  }

  // image with the variance of neighborhood 3x3x3 patch
  for(k=0;k<dims[2];k++)
  {
    for(j=0;j<dims[1];j++)
    {
      for(i=0;i<dims[0];i++)
      {
        var=0;
        indice=0;
        for(ii=-1;ii<=1;ii++)
        {
          for(jj=-1;jj<=1;jj++)
          {
            for(kk=-1;kk<=1;kk++)
            {
              ni=i+ii;
              nj=j+jj;
              nk=k+kk;
              if(ni>=0 && nj>=0 && nk>0 && ni<dims[0] && nj<dims[1] && nk<dims[2])
              {
                var = var + (ima[nk*(dims[0]*dims[1])+(nj*dims[0])+ni]-means[k*(dims[0]*dims[1])+(j*dims[0])+i])*(ima[nk*(dims[0]*dims[1])+(nj*dims[0])+ni]-means[k*(dims[0]*dims[1])+(j*dims[0])+i]);
                indice=indice+1;
              }
            }
          }
        }
        var=var/(indice-1);
        variances[k*(dims[0]*dims[1])+(j*dims[0])+i]=var;
      }
    }
  }

  // Threaded portion of the code
  itk::MultiThreaderBase::Pointer mt = itk::MultiThreaderBase::New();
  itk::ImageRegion<3>             full_region;
  for (int d = 0; d < 3; d++)
    full_region.SetSize(d, dims[d]);

  mt->ParallelizeImageRegion<3>(
    full_region,
    [dims, ima, fima, variances, means, Estimate, bias, Label, f, v, rician, max](
      const itk::ImageRegion<3> & thread_region) {
      ThreadArgument<TFloat> ta;
      ta.cols = dims[0];
      ta.rows = dims[1];
      ta.slices = dims[2];
      ta.in_image = ima;
      ta.var_image = variances;
      ta.means_image = means;
      ta.estimate = Estimate;
      ta.bias = bias;
      ta.label = Label;
      ta.ini = thread_region.GetIndex(2);
      ta.fin = thread_region.GetIndex(2) + thread_region.GetSize(2);
      ta.radioB = v;
      ta.radioS = f;
      ta.rician = rician;
      ta.max = max;
      ThreadFunc(&ta);
    },
    nullptr);

  // correction?
  if(rician)
  {
    r=5;
    Regularize(bias,variances,r,dims[0],dims[1],dims[2]);
    for(i=0;i<dims[0]*dims[1]*dims[2];i++)
    {
      if(variances[i]>0)
      {
        SNR=means[i]/sqrt(variances[i]);
        bias[i]=2*(variances[i]/Epsi(SNR));
        if(std::isnan(bias[i])) bias[i]=0;
      }
    }
  }

  /* Aggregation of the estimators (i.e. means computation) */
  label = 0.0;
  estimate = 0.0;
  for(i=0;i<dims[0]*dims[1]*dims[2];i++)
  {
    label = Label[i];
    if (label == 0.0) fima[i] = ima[i];
    else
    {
      estimate = Estimate[i];
      estimate = (estimate/label);
      if(rician)
      {
        estimate = (estimate-bias[i])<0?0:(estimate-bias[i]);
        fima[i]=sqrt(estimate);
      }
      else fima[i]=estimate;
    }
  }


  // free memory
  delete [] means;
  delete [] variances;
  delete [] Estimate;
  delete [] Label;
  if(rician) delete bias;
  delete [] average;

  // write output image
  return ITKfima;
}

template <class TPixel, unsigned int VDim>
class NLMDenoiseProblemTemplated
{
public:
  using ImageType = itk::OrientedRASImage<TPixel, VDim>;
  using ImagePointer = typename ImageType::Pointer;
  static ImagePointer execute(ImageType *image, DenoiseParameters param, std::ostream &sout)
  {
    throw ConvertException("NonLocalMeansDenoise only supported for 3D images");
  }
};

template <class TPixel>
class NLMDenoiseProblemTemplated<TPixel, 3>
{
public:
  using ImageType = itk::OrientedRASImage<TPixel, 3>;
  using ImagePointer = typename ImageType::Pointer;
  static ImagePointer execute(ImageType *image, DenoiseParameters param, std::ostream &sout)
  {
    return NLMDenoiseProblem<TPixel>::NLMdenoise(image, param, sout);
  }
};

/*

int main(int argc, char *argv[])
{
  if(argc < 2)
    return usage();

  DenoiseParameters param;

         // Process parameters
  CommandLineHelper cl(argc, argv);

  while(!cl.is_at_end())
  {
    // Read the next command
    std::string arg = cl.read_command();

    if(arg == "-i")
    {
      param.fnImage = cl.read_existing_filename();
    }
    else if(arg == "-o")
    {
      param.fnOutput = cl.read_output_filename();
    }
    else if(arg == "-v")
    {
      param.v = cl.read_integer();
    }
    else if(arg == "-f")
    {
      param.f = cl.read_integer();
    }
    else if(arg == "-r")
    {
      int rint = cl.read_integer();
      if( rint > 0 )
        param.rician = true;
      else
        param.rician = false;
    }
    else
    {
      cerr << "Unknown option " << arg << endl;
      return -1;
    }
  }

         // Check parameters
  check(param.v >= 1, "Search radius (v) must a positive integer.");
  check(param.f >= 1, "Patch radius (f) must a positive integer.");

         // perform denoising
  NLMDenoiseProblem<TFloat, 3>::NLMdenoise(param);
  cout << "Done NLM denoising." << endl;
}

*/

template <class TPixel, unsigned int VDim>
void
NonLocalMeansDenoise<TPixel, VDim>
::operator() (int search_radius, int patch_radius, bool rician)
{
  // This filter is only supported for 3D images
  if(VDim != 3)
    throw ConvertException("NonLocalMeansDenoise only supported for 3D images");

  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Populate the parameters
  DenoiseParameters param;
  param.v = search_radius;
  param.f = patch_radius;
  param.rician = rician;

  // Perform NLM
  *c->verbose << "Applying Manjon et al. (2009) Non-Local Means Denoising to #" << c->m_ImageStack.size() << endl;
  ImagePointer result = NLMDenoiseProblemTemplated<TPixel, VDim>::execute(img, param, *c->verbose);

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class NonLocalMeansDenoise<double, 2>;
template class NonLocalMeansDenoise<double, 3>;
template class NonLocalMeansDenoise<double, 4>;
