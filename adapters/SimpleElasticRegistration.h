/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SimpleElasticRegistration.h
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

#ifndef __SimpleElasticRegistration_h_
#define __SimpleElasticRegistration_h_

#include "ConvertAdapter.h"
#include "itkCovariantVector.h"
#include "itkBSplineInterpolateImageFunction.h"
#include <fftw3.h>

typedef double VecType;

#define DOUBLE_FFT
#ifdef DOUBLE_FFT

#define myfftw_complex fftw_complex
#define myfftw_plan fftw_plan
#define myfftw_malloc fftw_malloc
#define myfftw_execute fftw_execute
#define myfftw_plan_dft_r2c fftw_plan_dft_r2c
#define myfftw_plan_dft_c2r fftw_plan_dft_c2r
#define myfftw_destroy_plan fftw_destroy_plan

#else

#define myfftw_complex fftwf_complex
#define myfftw_plan fftwf_plan
#define myfftw_malloc fftwf_malloc
#define myfftw_execute fftwf_execute
#define myfftw_plan_dft_r2c fftwf_plan_dft_r2c
#define myfftw_plan_dft_c2r fftwf_plan_dft_c2r
#define myfftw_destroy_plan fftwf_destroy_plan

#endif

template<class TPixel, unsigned int VDim>
class SimpleElasticRegistration : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  // Other typedefs
  typedef itk::CovariantVector<VecType, VDim>  VectorType;
  typedef itk::OrientedRASImage<VectorType, VDim> VectorImageType;
  typedef typename VectorImageType::Pointer VectorImagePointer;
  typedef typename itk::ImageRegionIteratorWithIndex<VectorImageType> VectorIterator;
  typedef itk::OrientedRASImage<VecType, VDim> FloatImage;

  // typedef itk::GaussianInterpolateImageFunction<ImageType, double> GaussInterpolator;
  typedef itk::BSplineInterpolateImageFunction<ImageType, double> Interpolator;

  SimpleElasticRegistration(Converter *c) : c(c) {}

  void operator() ();



  Converter *c;

  double ComputeObjectiveAndGradient(bool need_gradient);
  static SimpleElasticRegistration<TPixel, VDim> *glob_ref;

  void SetX(double *x);
  void TestGradientComputation();
  void TestGradient2();
  void RandomVectorField();
  void SaveRaw(VecType *, const char *);
  void TestFourier();

  // Meaning of 'forward' flag:
  // true:  solve LL x = y for unknown x
  // false: compute y = LL x for unknown y
  void SolvePDE(myfftw_plan &plan, bool forward); 

  // Stuff for optimization
  ImagePointer iref, imsk, imov;
  typename FloatImage::Pointer ikernel;

  // The vector field can't be stored as an image of vectors because of FFT. Instead,
  // we store it as a bunch of floating point images
  typename FloatImage::Pointer v[VDim], g[VDim];

  // Interpolators
  typename Interpolator::Pointer giref, gimov;

  // FFTW plans
  myfftw_plan f_fft_v[VDim], f_fft_g[VDim], i_fft; 

  // Sizes of arrays for fft, etc
  int fft_sz[VDim], fft_szc[VDim];
  size_t fft_n, fft_nc;
  
  // Scaling arrays for laplacian operator computations
  vnl_vector<double> kx[VDim], cx[VDim];

  // Complex array
  myfftw_complex *fcmp;
  VecType *vraw[VDim], *graw[VDim], *vtmp;

  // Alpha, gamma factors
  double alpha, gamma, sigma;


  // Mask scaling factor
  double mask_scale;
};

#endif

