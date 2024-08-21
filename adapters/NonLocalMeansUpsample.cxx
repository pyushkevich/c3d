/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    NonLocalMeansUpsample.cxx
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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%
//%   Jose V. Manjon - jmanjon@fis.upv.es
//%   Universidad Politecinca de Valencia, Spain
//%   Pierrick Coupe - pierrick.coupe@gmail.com
//%   Brain Imaging Center, Montreal Neurological Institute.
//%   Mc Gill University
//%
//%   Copyright (C) 2010 Jose V. Manjon and Pierrick Coupe
//%
//%    Usage of NLMUpsample:
//%
//%    fima=NLMUpsample(ima,f)
//%
//%    ima: LR volume
//%    f: Magnification factor in each dimension (for example [2,2,2])
//%    fima: HR upsampled volume
//%
//%**************************************************************************

#include "NonLocalMeansUpsample.h"

#include <iostream>
#include <string>
#include "math.h"
#include <cstdlib>
#include "itkResampleImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkNoiseImageFilter.h"

/* Multithreading stuff*/
#ifdef _WIN32
#  include <windows.h>
#  include <process.h>
#else
#  include <pthread.h>
#endif

#define pi 3.1415926535

using namespace std;

struct UpsampleParameters
{
  int    v;
  int    f;
  int    lf_x;
  int    lf_y;
  int    lf_z;

  UpsampleParameters()
    : v(3)
    , f(1)
    , lf_x(2)
    , lf_y(2)
    , lf_z(2)
  {}
};

template <class TFloat>
struct ThreadArgument
{
  int      rows;
  int      cols;
  int      slices;
  double * in_image;
  double * out_image;
  double * mean_image;
  double * pesos;
  int      ini;
  int      fin;
  int      radio;
  int      f;
  /*int th;*/
  double * sigma;
};


template <class TFloat>
class NLMUpsampleProblem
{
public:
  typedef itk::OrientedRASImage<TFloat, 3>          ImageType;
  typedef typename ImageType::Pointer       ImagePointer;
  using myargument = ThreadArgument<TFloat>;

  static TFloat
  distancia(TFloat * ima, int x, int y, int z, int nx, int ny, int nz, int f, int sx, int sy, int sz);

  static void *
  ThreadFunc(void * pArguments);

  static void
  PatchCorrection(TFloat * ima, TFloat * fima, int * dims, TFloat * h, UpsampleParameters param);

  static ImagePointer
  Resample(ImagePointer ITKima, UpsampleParameters param);

  static void
  MeanCorrection(TFloat * bima, TFloat * ref, TFloat * fima, int * dims, UpsampleParameters param);

  static ImagePointer
  ComputeLevel(ImagePointer input, int radius);

  static ImagePointer
  main_loop(ImageType *input, UpsampleParameters param, std::ostream &sout);

private:
};


template <class TFloat>
TFloat
NLMUpsampleProblem<TFloat>::distancia(TFloat * ima,
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
  TFloat d, acu, distancetotal, inc;
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
void *
NLMUpsampleProblem<TFloat>::ThreadFunc(void * pArguments)
{
  TFloat *ima, *fima, *medias, *pesos, w, d, *h, t1;
  int     ii, jj, kk, ni, nj, nk, i, j, k, ini, fin, rows, cols, slices, v, p, p1, f, rc;

  myargument arg;
  arg = *(myargument *)pArguments;

  rows = arg.rows;
  cols = arg.cols;
  slices = arg.slices;
  ini = arg.ini;
  fin = arg.fin;
  ima = arg.in_image;
  fima = arg.out_image;
  medias = arg.mean_image;
  pesos = arg.pesos;
  v = arg.radio;
  f = arg.f;
  /*th=arg.th;*/
  h = arg.sigma;
  rc = rows * cols;

  /* filter*/
  for (k = ini; k < fin; k++)
  {
    for (j = 0; j < rows; j++)
    {
      for (i = 0; i < cols; i++)
      {
        p = k * rc + j * cols + i;
        if (h[p] < 1)
          continue;

        for (kk = 0; kk <= v; kk++)
        {
          nk = k + kk;
          for (ii = -v; ii <= v; ii++)
          {
            ni = i + ii;
            for (jj = -v; jj <= v; jj++)
            {
              nj = j + jj;

              if (kk == 0 && jj < 0)
                continue;
              if (kk == 0 && jj == 0 && ii <= 0)
                continue;

              if (ni >= 0 && nj >= 0 && nk >= 0 && ni < cols && nj < rows && nk < slices)
              {
                p1 = nk * rc + nj * cols + ni;

                t1 = abs(medias[p] - medias[p1]);

                if (t1 > 0.6 * h[p])
                  continue;

                d = distancia(ima, i, j, k, ni, nj, nk, f, cols, rows, slices);

                if (h[p] > 0)
                  d = d / (2 * h[p] * h[p]) - 1;
                if (d < 0)
                  d = 0;

                w = exp(-d);

                fima[p] = fima[p] + w * ima[p1];
                pesos[p] = pesos[p] + w;

                fima[p1] = fima[p1] + w * ima[p];
                pesos[p1] = pesos[p1] + w;
              }
            }
          }
        }
      }
    }
  }

#ifdef _WIN32
  _endthreadex(0);
#else
  pthread_exit(0);
#endif

  return 0;
}

template <class TFloat>
void
NLMUpsampleProblem<TFloat>::PatchCorrection(TFloat *           ima,
                                            TFloat *           fima,
                                            int *              dims,
                                            TFloat *           h,
                                            UpsampleParameters param)
{
  // read inputs
  int v = param.v;
  int f = param.f;
  int fac[3];
  fac[0] = param.lf_x;
  fac[1] = param.lf_y;
  fac[2] = param.lf_z;
  int nelement = dims[2] * dims[1] * dims[0];

  // allocate memory for intermediate matrices
  TFloat * tmp = new TFloat[nelement];
  TFloat * medias = new TFloat[nelement];
  TFloat * pesos = new TFloat[nelement];

  /*Declarations*/
  TFloat off, media;
  int    ini, fin, i, j, k, ii, jj, kk, ni, nj, nk, indice, Nthreads, rc, ft;
  bool   salir;

  // start computation
  rc = dims[0] * dims[1];
  for (k = 0; k < dims[2]; k++)
  {
    for (j = 0; j < dims[1]; j++)
    {
      for (i = 0; i < dims[0]; i++)
      {
        media = 0;
        for (ii = -1; ii <= 1; ii++)
        {
          ni = i + ii;
          if (ni < 0)
            ni = -ni;
          if (ni >= dims[0])
            ni = 2 * dims[0] - ni - 1;
          for (jj = -1; jj <= 1; jj++)
          {
            nj = j + jj;
            if (nj < 0)
              nj = -nj;
            if (nj >= dims[1])
              nj = 2 * dims[1] - nj - 1;
            for (kk = -1; kk <= 1; kk++)
            {
              nk = k + kk;
              if (nk < 0)
                nk = -nk;
              if (nk >= dims[2])
                nk = 2 * dims[2] - nk - 1;

              media = media + ima[nk * rc + nj * dims[0] + ni];
            }
          }
        }
        medias[k * rc + j * dims[0] + i] = media / 27;
      }
    }
  }

  ft = fac[2] * fac[1] * fac[0];
  for (k = 0; k < dims[2] / fac[2]; k++)
    for (j = 0; j < dims[1] / fac[1]; j++)
      for (i = 0; i < dims[0] / fac[0]; i++)
      {
        media = 0;
        for (kk = 0; kk < fac[2]; kk++)
          for (jj = 0; jj < fac[1]; jj++)
            for (ii = 0; ii < fac[0]; ii++)
              media += ima[(k * fac[2] + kk) * rc + (j * fac[1] + jj) * dims[0] + i * fac[0] + ii];
        tmp[k * rc + (j * dims[0]) + i] = media / ft;
      }

  for (k = 0; k < dims[2] * dims[1] * dims[0]; k++)
  {
    fima[k] = ima[k];
    pesos[k] = 1.0;
  }
  /*th=0.6*h;  //3/sqrt(27) */

  // start multi-thread computation
  itk::MultiThreaderBase::Pointer mt = itk::MultiThreaderBase::New();
  itk::ImageRegion<3>             full_region;
  for (int d = 0; d < 3; d++)
    full_region.SetSize(d, dims[d]);
  
  mt->ParallelizeImageRegion<3>(
    full_region,
    [dims, ima, fima, medias, pesos, v, f, h](const itk::ImageRegion<3> & thread_region) 
      { 
      myargument ta;
      ta.cols = dims[0];
      ta.rows = dims[1];
      ta.slices = dims[2];
      ta.in_image = ima;
      ta.out_image = fima;
      ta.mean_image = medias;
      ta.pesos = pesos;
      ta.ini = thread_region.GetIndex(2);
      ta.fin = thread_region.GetIndex(2) + thread_region.GetSize(2);
      ta.radio = v;
      ta.f = f;
      ta.sigma = h;
      ThreadFunc(&ta);
    },
    nullptr);



  for (i = 0; i < dims[0] * dims[1] * dims[2]; i++)
    fima[i] /= pesos[i];


  /* apply mean preservation constraint*/
  for (k = 0; k < dims[2]; k = k + fac[2])
    for (j = 0; j < dims[1]; j = j + fac[1])
      for (i = 0; i < dims[0]; i = i + fac[0])
      {
        salir = false;

        media = 0;
        for (kk = 0; kk < fac[2]; kk++)
          for (jj = 0; jj < fac[1]; jj++)
            for (ii = 0; ii < fac[0]; ii++)
              media += fima[(k + kk) * rc + (j + jj) * dims[0] + i + ii];
        media = media / ft;

        off = tmp[(k / fac[2]) * rc + (j / fac[1]) * dims[0] + (i / fac[0])] - media;

        for (kk = 0; kk < fac[2]; kk++)
          for (jj = 0; jj < fac[1]; jj++)
            for (ii = 0; ii < fac[0]; ii++)
            {
              fima[(k + kk) * rc + (j + jj) * dims[0] + (i + ii)] += off;
            }
      }

  // free memory
  delete[] tmp;
  delete[] medias;
  delete[] pesos;

  // return;
}

template <class TFloat>
typename NLMUpsampleProblem<TFloat>::ImagePointer
NLMUpsampleProblem<TFloat>::Resample(ImagePointer ITKima, UpsampleParameters param)
{
  // code from Paul's c3d

  // Get the size of the image in voxels
  typename ImageType::SizeType sz;
  sz[0] = (unsigned long)(ITKima->GetBufferedRegion().GetSize(0) * param.lf_x + 0.5);
  sz[1] = (unsigned long)(ITKima->GetBufferedRegion().GetSize(1) * param.lf_y + 0.5);
  sz[2] = (unsigned long)(ITKima->GetBufferedRegion().GetSize(2) * param.lf_z + 0.5);

  // Build the resampling filter
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer                   fltSample = ResampleFilterType::New();

  // Specify interpolator
  typedef itk::BSplineInterpolateImageFunction<ImageType, TFloat> CubicInterpolatorType;

  // Initialize the resampling filter with an identity transform
  fltSample->SetInput(ITKima);
  fltSample->SetTransform(itk::IdentityTransform<TFloat, 3>::New());
  fltSample->SetInterpolator(CubicInterpolatorType::New());

  // Compute the spacing of the new image
  typename ImageType::SpacingType spc_pre = ITKima->GetSpacing();
  typename ImageType::SpacingType spc_post = spc_pre;
  for (size_t i = 0; i < 3; i++)
    spc_post[i] *= ITKima->GetBufferedRegion().GetSize()[i] * 1.0 / sz[i];

  // Get the bounding box of the input image
  typename ImageType::PointType origin_pre = ITKima->GetOrigin();

  // Recalculate the origin. The origin describes the center of voxel 0,0,0
  // so that as the voxel size changes, the origin will change as well.
  typename ImageType::SpacingType off_pre = (ITKima->GetDirection() * spc_pre) * 0.5;
  typename ImageType::SpacingType off_post = (ITKima->GetDirection() * spc_post) * 0.5;
  typename ImageType::PointType   origin_post = origin_pre - off_pre + off_post;

  // Set the image sizes and spacing.
  fltSample->SetSize(sz);
  fltSample->SetOutputSpacing(spc_post);
  fltSample->SetOutputOrigin(origin_post);
  fltSample->SetOutputDirection(ITKima->GetDirection());

  // Perform resampling
  fltSample->UpdateLargestPossibleRegion();
  return fltSample->GetOutput();
}

template <class TFloat>
void
NLMUpsampleProblem<TFloat>::MeanCorrection(TFloat *           bima,
                                           TFloat *           ref,
                                           TFloat *           fima,
                                           int *              dims,
                                           UpsampleParameters param)
{
  /* Declarations */
  TFloat media, off;
  int    i, j, k, ii, jj, kk;

  // compute size
  int lf[3];
  lf[0] = param.lf_x;
  lf[1] = param.lf_y;
  lf[2] = param.lf_z;
  int sx = dims[0] / lf[0];
  int sy = dims[1] / lf[1];
  int sz = dims[2] / lf[2];

  /*  mean correction */
  for (k = 0; k < dims[2]; k = k + lf[2])
    for (j = 0; j < dims[1]; j = j + lf[1])
      for (i = 0; i < dims[0]; i = i + lf[0])
      {
        media = 0;
        for (kk = 0; kk < lf[2]; kk++)
          for (jj = 0; jj < lf[1]; jj++)
            for (ii = 0; ii < lf[0]; ii++)
              media = media + bima[(k + kk) * dims[0] * dims[1] + (j + jj) * dims[0] + i + ii];
        media = media / (lf[0] * lf[1] * lf[2]);

        off = ref[(int)((k / lf[2]) * sx * sy + (j / lf[1]) * sx + (i / lf[0]))] - media;

        for (kk = 0; kk < lf[2]; kk++)
          for (jj = 0; jj < lf[1]; jj++)
            for (ii = 0; ii < lf[0]; ii++)
              fima[(k + kk) * dims[0] * dims[1] + (j + jj) * dims[0] + i + ii] =
                bima[(k + kk) * dims[0] * dims[1] + (j + jj) * dims[0] + i + ii] + off;
      }
}

template <class TFloat>
typename NLMUpsampleProblem<TFloat>::ImagePointer
NLMUpsampleProblem<TFloat>::ComputeLevel(ImagePointer input, int radius)
{
  // compute local std map
  typedef itk::NoiseImageFilter<ImageType, ImageType> NoiseImageFilterType;
  typename NoiseImageFilterType::Pointer              noiseImageFilter = NoiseImageFilterType::New();
  noiseImageFilter->SetInput(input);
  noiseImageFilter->SetRadius(radius);
  noiseImageFilter->Update();

  // perform mean filter
  typedef itk::MeanImageFilter<ImageType, ImageType> MeanImageFilterType;
  typename MeanImageFilterType::Pointer              meanImageFilter = MeanImageFilterType::New();
  meanImageFilter->SetInput(noiseImageFilter->GetOutput());
  meanImageFilter->SetRadius(radius);
  meanImageFilter->Update();

  return meanImageFilter->GetOutput();
}


template <class TFloat>
typename NLMUpsampleProblem<TFloat>::ImagePointer
NLMUpsampleProblem<TFloat>::main_loop(ImageType *input, UpsampleParameters param, std::ostream &sout)
{
  ImagePointer ITKima = input;
  TFloat *     ima = ITKima->GetBufferPointer();

  // get image size
  itk::Size<3> ITKdims = ITKima->GetBufferedRegion().GetSize();
  int dims[3], n = 1;
  for(unsigned int i = 0; i < 3; i++)
  {
    dims[i] = ITKdims[i];
    n *= dims[i];
  }

  // normalize signal to [0 256]
  TFloat iMin = ima[0], iMax = ima[0];
  for (size_t i = 0; i < n; i++)
  {
    iMax = (iMax > ima[i]) ? iMax : ima[i];
    iMin = (iMin < ima[i]) ? iMin : ima[i];
  }
  if (iMin == iMax)
    throw ConvertException("NLMUpsampleProblem: the image has no intensity range");

  for (size_t i = 0; i < n; i++)
    ima[i] = 256 * (ima[i] - iMin) / (iMax - iMin);

  // perform initial interpolation
  ImagePointer ITKbima = Resample(ITKima, param);
  TFloat *     bima = ITKbima->GetBufferPointer();

  // update dimension info
  ITKdims = ITKbima->GetBufferedRegion().GetSize();
  n = 1;
  for(unsigned int i = 0; i < 3; i++)
  {
    dims[i] = ITKdims[i];
    n *= dims[i];
  }

  // allocate memory for output image
  ImagePointer ITKfima = ImageType::New();
  ITKfima->SetRegions(ITKbima->GetBufferedRegion());
  ITKfima->SetSpacing(ITKbima->GetSpacing());
  ITKfima->SetOrigin(ITKbima->GetOrigin());
  ITKfima->SetDirection(ITKbima->GetDirection());
  ITKfima->SetMetaDataDictionary(ITKbima->GetMetaDataDictionary());
  ITKfima->Allocate();
  TFloat* fima = ITKfima->GetBufferPointer();

  // perform mean correction
  MeanCorrection(bima, ima, fima, dims, param);

  // compute level (sigma
  int          radius = 1;
  ImagePointer ITKlevel = ComputeLevel(ITKfima, radius);
  TFloat *     level = ITKlevel->GetBufferPointer();

  // perform iteration
  int      maxiter = 1000;
  std::vector<TFloat> mean(maxiter), mean_down(maxiter);
  TFloat   tol = 1.2;
  int      down = 0;
  int      ii = 1, iii = 1;
  TFloat * tmpima = new TFloat[n];
  for (size_t i = 0; i < n; i++)
  {
    bima[i] = fima[i];
    tmpima[i] = fima[i];
  }

  while(true)
  {
    // perform upsample
    PatchCorrection(bima, fima, dims, level, param);

    // compute average differences between ima and fima
    mean[ii] = 0;
    for (size_t i = 0; i < n; i++)
      mean[ii] += abs(bima[i] - fima[i]);
    mean[ii] /= n;
    sout << "  Iter " << ii << ": down fac = " << iii << "; abs mean diff = " << mean[ii] << "." << endl;

    // determine continue or not
    if (ii > 1)
    {
      if (mean[ii - 1] / mean[ii] < tol && down == 0)
      {
        down = 1;
        for (size_t i = 0; i < n; i++)
        {
          level[i] = level[i] / 2;
          mean_down[iii] += abs(tmpima[i] - fima[i]);
          tmpima[i] = fima[i];
        }
        if (iii > 1)
          if (mean_down[iii - 1] / mean_down[iii] < tol)
            break;
        iii += 1;
      }
      else
        down = 0;

      if (mean[ii] <= 0.001)
        break;

      if (ii >= maxiter)
        break;
    }

    // prepare for the next iteration
    for (size_t i = 0; i < n; i++)
      bima[i] = fima[i];
    ii += 1;
  }

  // rescale it back to original range
  for (size_t i = 0; i < n; i++)
    fima[i] = iMin + (iMax-iMin) * (fima[i] / 256);

  // Return the upsampled image
  return ITKfima;
}

template <class TPixel, unsigned int VDim>
class NLMUpsampleProblemTemplated
{
public:
  using ImageType = itk::OrientedRASImage<TPixel, VDim>;
  using ImagePointer = typename ImageType::Pointer;
  static ImagePointer execute(ImageType *image, UpsampleParameters param, std::ostream &sout)
  {
    throw ConvertException("NonLocalMeansUpsample only supported for 3D images");
  }
};

template <class TPixel>
class NLMUpsampleProblemTemplated<TPixel, 3>
{
public:
  using ImageType = itk::OrientedRASImage<TPixel, 3>;
  using ImagePointer = typename ImageType::Pointer;
  static ImagePointer execute(ImageType *image, UpsampleParameters param, std::ostream &sout)
  {
    return NLMUpsampleProblem<TPixel>::main_loop(image, param, sout);
  }
};


template <class TPixel, unsigned int VDim>
void
NonLocalMeansUpsample<TPixel, VDim>::operator()(SizeType factor, int search_radius, int patch_radius)
{
  // This filter is only supported for 3D images
  if(VDim != 3)
    throw ConvertException("NonLocalMeansUpsample only supported for 3D images");

  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Populate the parameters
  UpsampleParameters param;
  param.lf_x = factor[0];
  param.lf_y = factor[1];
  param.lf_z = factor[2];
  param.v = search_radius;
  param.f = patch_radius;

  // Perform NLM
  *c->verbose << "Applying Manjon et al. (2010) Non-Local Means Upsampling to #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Upsampling factor: " << factor << endl;
  ImagePointer result = NLMUpsampleProblemTemplated<TPixel, VDim>::execute(img, param, *c->verbose);

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class NonLocalMeansUpsample<double, 2>;
template class NonLocalMeansUpsample<double, 3>;
template class NonLocalMeansUpsample<double, 4>;
