/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SimpleElasticRegistration.cxx
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

#include "SimpleElasticRegistration.h"
#include "itkFFTWRealToComplexConjugateImageFilter.h"
#include "itkFFTWComplexConjugateToRealImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageFileWriter.h"
#include <vnl/vnl_random.h>

extern "C" {
#include "cg_user.h"
}

template <class TPixel, unsigned int VDim>
double
cg_func(double *x, long int n)
{
  SimpleElasticRegistration<TPixel, VDim> *o = SimpleElasticRegistration<TPixel, VDim>::glob_ref;
  o->SetX(x);
  return o->ComputeObjectiveAndGradient(false);
}

template <class TPixel, unsigned int VDim>
double
cg_valgrad(double *g, double *x, long int n)
{
  // Compute the objective
  SimpleElasticRegistration<TPixel, VDim> *o = SimpleElasticRegistration<TPixel, VDim>::glob_ref;
  o->SetX(x);
  double rc = o->ComputeObjectiveAndGradient(true);

  // Fill in the gradient
  for(size_t d = 0; d < VDim; d++)
    for(size_t i = 0; i < o->fft_n; i++, g++)
      *g = o->graw[d][i];

  return rc;
}

template <class TPixel, unsigned int VDim>
void 
cg_grad(double *g, double *x, long int n)
{
  cg_valgrad<TPixel, VDim>(g, x, n);
}


template <class TPixel, unsigned int VDim>
SimpleElasticRegistration<TPixel, VDim> *
SimpleElasticRegistration<TPixel, VDim>::glob_ref = NULL;

template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::operator() ()
{
  // Get reference and moving images from stack
  this->imov = c->m_ImageStack.back(); c->m_ImageStack.pop_back();
  this->imsk = c->m_ImageStack.back(); c->m_ImageStack.pop_back();
  this->iref = c->m_ImageStack.back(); c->m_ImageStack.pop_back();

  // Save a pointer to self for optimization
  glob_ref = this;

  // Initialize parameters
  alpha = 40.96; gamma = 1; sigma = 2.0;
  
  // Initialize the Fourier transform code
  // We have some trouble. The FFTW format is row-major and ITK is column-major
  // We simply reverse the size and spacing arrays to deal with it
  double fft_sp[VDim];
  for(size_t d = 0; d < VDim; d++)
    { 
    fft_szc[d] = fft_sz[d] = iref->GetBufferedRegion().GetSize()[VDim - (d+1)];
    fft_sp[d] = iref->GetSpacing()[VDim - (d+1)];
    }

  // Adjust the size along the last dimension (per FFTW)
  fft_szc[VDim-1] = fft_sz[VDim-1] / 2 + 1;

  // Compute size of real/complex arrays
  fft_n = iref->GetBufferedRegion().GetNumberOfPixels();
  fft_nc = fft_n * fft_szc[VDim-1] / fft_sz[VDim - 1];
  
  // Create the output complex array
  fcmp = (myfftw_complex *) myfftw_malloc(fft_nc * sizeof(myfftw_complex));

  // Initialize all the vector fields
  for(size_t d = 0; d < VDim; d++)
    {
    // The displacement field that's optimized over
    vraw[d] = (VecType *) myfftw_malloc(fft_n * sizeof(VecType));
    v[d] = FloatImage::New();
    v[d]->SetRegions(iref->GetBufferedRegion());
    v[d]->GetPixelContainer()->SetImportPointer(vraw[d], fft_n, true);
    v[d]->FillBuffer(0.0f);

    // The gradient vector field
    graw[d] = (VecType *) myfftw_malloc(fft_n * sizeof(VecType));
    g[d] = FloatImage::New();
    g[d]->SetRegions(iref->GetBufferedRegion());
    g[d]->GetPixelContainer()->SetImportPointer(graw[d], fft_n, true);
    g[d]->FillBuffer(0.0f);

    // Create the FFT plans for each FFT we will perform
    f_fft_v[d] = myfftw_plan_dft_r2c(VDim, fft_sz, vraw[d], fcmp, FFTW_ESTIMATE);
    f_fft_g[d] = myfftw_plan_dft_r2c(VDim, fft_sz, graw[d], fcmp, FFTW_ESTIMATE);

    // Create frequency arrays and cosine arrays for laplacian computation
    size_t N = fft_sz[d];
    double c = 6.283185307179586 / (N * fft_sp[d]);
    kx[d].set_size(N);
    cx[d].set_size(N);
    for(size_t j = 0; j <= N - N/2; j++)
      kx[d][j] = c * j;
    for(size_t j = 1; j <= N/2; j++)
      kx[d][N-j] = -kx[d][j];
    for(size_t j = 0; j < N; j++)
      {
      cx[d][j] = 2.0 * (1.0 - cos(kx[d][j]));
      }
    }
  
  // Create the temporary float array for FFT
  vtmp = (VecType *) myfftw_malloc(fft_n * sizeof(VecType));

  // Create the inverse fft
  i_fft = myfftw_plan_dft_c2r(VDim, fft_sz, fcmp, vtmp, FFTW_ESTIMATE);

  // Create the sampling function
  giref = Interpolator::New();
  giref->SetSplineOrder(3);
  giref->SetInputImage(iref);

  gimov = Interpolator::New();
  gimov->SetSplineOrder(3);
  gimov->SetInputImage(imov);

  // Integrate the mask to save time
  double mask_sum = 0.0;
  for(Iterator itmsk(imsk, imsk->GetBufferedRegion()); !itmsk.IsAtEnd(); ++itmsk)
    mask_sum += itmsk.Value();
  if(mask_sum <= 0.0)
    throw(ConvertException("Mask is empty or negative in SimpleElasticRegistration"));
  this->mask_scale = 1.0 / mask_sum;

  // Test interpolator derivative
  itk::ContinuousIndex<double, VDim> q;
  for(size_t d = 0; d < VDim; d++)
    q[d] = iref->GetBufferedRegion().GetIndex()[d] + 0.5 * iref->GetBufferedRegion().GetSize()[d];
  cout << "F[q] " << gimov->EvaluateAtContinuousIndex(q) << endl;
  cout << "G[q] " << gimov->EvaluateDerivativeAtContinuousIndex(q) << endl;
  for(size_t d = 0; d < VDim; d++)
    {
    double eps = 1.0e-6;
    itk::ContinuousIndex<double, VDim> dq = q;
    dq[d] = q[d] + eps / imov->GetSpacing()[d];
    double f1 = gimov->EvaluateAtContinuousIndex(dq);
    dq[d] = q[d] - eps / imov->GetSpacing()[d];
    double f2 = gimov->EvaluateAtContinuousIndex(dq);
    cout << "Dx = " << (f1-f2) / (2 * eps) << endl;
    }

  // Another test
    {
    // Define point x in reference space
    itk::Index<VDim> idx;
    for(size_t d = 0; d < VDim; d++)
      idx[d] = iref->GetBufferedRegion().GetIndex()[d] + iref->GetBufferedRegion().GetSize()[d] / 2;
    itk::Point<double, VDim> x;
    itk::ContinuousIndex<double, VDim> qx, qy;
    iref->TransformIndexToRASPhysicalPoint(idx, x);
    iref->TransformRASPhysicalPointToContinuousIndex(x, qx);
    cout << "X = " << x << endl;
    cout << "QX = " << qx << endl;

    // Define the vector v
    double w[] = {4.0, 4.0, 4.0};
    itk::Point<double, VDim> y;
    for(size_t d = 0; d < VDim; d++)
      y[d] = x[d] + w[d];
    iref->TransformRASPhysicalPointToContinuousIndex(y, qy);
    cout << "Y = " << y << endl;
    cout << "QY = " << qy << endl;

    double f0 = pow(giref->EvaluateAtContinuousIndex(qx) - gimov->EvaluateAtContinuousIndex(qy),2);

    // Define an offset vector for gradient computation
    double off[] = {0.5, -1.0, 2.0}, eps = 1.0e-6;

    // Compute the gradient
    itk::CovariantVector<double, VDim> gmov = gimov->EvaluateDerivativeAtContinuousIndex(qy);
    cout << "f0 " << f0 << endl;
    cout << "G[y] " << gmov << endl;
    
    double ad = 0.0;
    for(size_t d = 0; d < 3; d++)
      {
      ad += off[d] * -2.0 * (giref->EvaluateAtContinuousIndex(qx) - gimov->EvaluateAtContinuousIndex(qy)) * gmov[d];
      }

    // Shift by off
    for(size_t d = 0; d < VDim; d++)
      y[d] = x[d] + w[d] + eps * off[d];
    iref->TransformRASPhysicalPointToContinuousIndex(y, qy);
    double f1 = pow(giref->EvaluateAtContinuousIndex(qx) - gimov->EvaluateAtContinuousIndex(qy),2);
    cout << "qy " << qy << endl;
    cout << "f1 - f0 " << f1 - f0 << endl;

    // Shift by off
    for(size_t d = 0; d < VDim; d++)
      y[d] = x[d] + w[d] - eps * off[d];
    iref->TransformRASPhysicalPointToContinuousIndex(y, qy);
    double f2 = pow(giref->EvaluateAtContinuousIndex(qx) - gimov->EvaluateAtContinuousIndex(qy),2);
    cout << "qy " << qy << endl;
    cout << "f2 - f0 " << f2 - f0 << endl;

    double nd = (f1 - f2) / (eps * 2.0);

    cout << "AD = " << ad << ", CD = " << nd << " DIFF " << ad - nd << endl;

    }

  // Test Fourier
  // TestFourier();
  // TestGradientComputation();
  // TestGradient2();

  // Reset the vector field to zero
  for(size_t i = 0; i < fft_n; i++)
    for(size_t d = 0; d < VDim; d++)
      vraw[d][i] = 0.0;

  // Compute the objective function
  *c->verbose << "Performing elastic registration" << endl;
  *c->verbose << "  Initial Objective " << ComputeObjectiveAndGradient(true) << endl;

  // Perform optimization
  double *xopt = new double[fft_n * VDim];
  memset(xopt, 0, sizeof(double) * fft_n * VDim);

  cg_parameter cg_parm;
  cg_default(&cg_parm);
  cg_parm.PrintLevel=0;
  cg_parm.maxit_fac = 100.0 / (VDim * fft_n);
  cg_descent(xopt, fft_n * VDim, NULL, &cg_parm, 1.e-8, 
    cg_func<TPixel,VDim>, 
    cg_grad<TPixel,VDim>,
    cg_valgrad<TPixel,VDim>, NULL);

  // Do simple gradient descent
  double dt = 0.01;
  for(size_t i = 0; i < 40; i++)
    {
    for(size_t d = 0; d < 3; d++)
      {
      for(size_t k = 0; k < fft_n; k++)
        {
        vraw[d][k] -= graw[d][k] * dt;
        }
      }

    double f = ComputeObjectiveAndGradient(true);
    *c->verbose << "  Iteration " << i << "\t Objective " << f << endl;
    }

  // Save transformed image
  for(size_t d = 0; d < 3; d++)
    {
    Iterator itref(iref, iref->GetBufferedRegion());
    size_t k = 0;
    for(; !itref.IsAtEnd(); ++itref, ++k)
      {
      // Check if we are inside the mask
      double vmask = imsk->GetBufferPointer()[k];
      if(vmask == 0.0)
        {
        itref.Set(0.0);
        continue;
        }

      // Get the index in reference space
      itk::Point<double, VDim> p;
      itk::ContinuousIndex<double, VDim> q, qr;
      IndexType idx_ref = itref.GetIndex();
      iref->TransformIndexToRASPhysicalPoint(idx_ref, p);
      iref->TransformRASPhysicalPointToContinuousIndex(p, qr);
      for(size_t d = 0; d < VDim; d++)
        p[d] += vraw[d][k];
      imov->TransformRASPhysicalPointToContinuousIndex(p, q);

      // Get the two image values
      if(gimov->IsInsideBuffer(q))
        {
        // Update the objective
        double I_mov = gimov->EvaluateAtContinuousIndex(q);
        itref.Set(I_mov);
        }
      else
        {
        itref.Set(0.0);
        }
      }
    }

  // Do some processing ...
  // ImagePointer result = ...;

  string name[] = {"warpxvec.nii", "warpyvec.nii", "warpzvec.nii"}; 
  for(size_t d = 0; d < 3; d++)
    SaveRaw(vraw[d], name[d].c_str());
  
  // Put result on stack
  // c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(iref);
}

template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::SetX(double *x)
{
  // X is a vector over which we optimize. The actual deformation vector
  // has the form V = K X, where K is the inverse of L'L, i.e., X is the 
  // solution of the PDE LL V = X
  
  // Repeat over each dimension
  for(size_t d = 0; d < VDim; d++)
    {
    // Copy x[d] into v[d]
    for(size_t i = 0; i < fft_n; i++, x++)
      vraw[d][i] = *x;

    // Solve the PDE
    SolvePDE(f_fft_v[d], true);

    // Put the solution into v
    for(size_t i = 0; i < fft_n; i++)
      vraw[d][i] = vtmp[i] / fft_n;
    }
}


template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::SolvePDE(myfftw_plan &plan, bool forward)
{
  // Take the forward FFT of v[d] into fcmp
  myfftw_execute(plan);

  // Create index
  int index[VDim];
  for(size_t d=0;d<VDim;d++)
    index[d]=0;

  // Apply operator LL to v, where L = (gamma I - alpha del)
  for(size_t j = 0; j < fft_nc; j++)
    {
    double w = gamma;
    for(size_t d = 0; d < VDim; d++)
      w += alpha * cx[d][index[d]];
    double wtotal = w * w;

    if(forward)
      wtotal = 1.0 / wtotal;

    fcmp[j][0] *= wtotal;
    fcmp[j][1] *= wtotal;

    // Increment index
    for(size_t d = VDim-1; d >= 0; d--)
      {
      if(++index[d] == fft_szc[d])
        index[d] = 0;
      else
        break;
      }
    }

  // Take the inverse fft, put it into vtmp
  myfftw_execute(i_fft);
}

template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::SaveRaw(VecType *v, const char *fn)
{
  typedef itk::Image<VecType, VDim> IType;
  typedef itk::ImageFileWriter<IType> WType;
  typename IType::Pointer i = IType::New();
  itk::ImageRegion<VDim> r;
  for(size_t d = 0; d < VDim; d++)
    {
    r.SetIndex(d,0);
    r.SetSize(d, iref->GetBufferedRegion().GetSize(d));
    }
  i->SetRegions(r);
  i->GetPixelContainer()->SetImportPointer(v, false);
  typename WType::Pointer w = WType::New();
  w->SetFileName(fn);
  w->SetInput(i);
  w->Update();
}


template <class TPixel, unsigned int VDim>
double
SimpleElasticRegistration<TPixel, VDim>
::ComputeObjectiveAndGradient(bool need_gradient)
{
  // Compute the gradient of the image objective
  double f_img = 0.0;
  double sfactor = 2.0 / (sigma * sigma);
  Iterator itref(iref, iref->GetBufferedRegion());
  size_t k = 0;

  // Clear the gradient
  if(need_gradient)
    {
    for(size_t d = 0; d < VDim; d++)
      for(size_t i = 0; i < fft_n; i++)
        graw[d][i] = 0.0;
    }

  for(; !itref.IsAtEnd(); ++itref, ++k)
    {
    // Check if we are inside the mask
    double vmask = imsk->GetBufferPointer()[k];
    if(vmask == 0.0)
      continue;

    // Get the index in reference space
    itk::Point<double, VDim> p;
    itk::ContinuousIndex<double, VDim> q, qr;
    IndexType idx_ref = itref.GetIndex();
    iref->TransformIndexToRASPhysicalPoint(idx_ref, p);
    iref->TransformRASPhysicalPointToContinuousIndex(p, qr);
    for(size_t d = 0; d < VDim; d++)
      p[d] += vraw[d][k];
    imov->TransformRASPhysicalPointToContinuousIndex(p, q);

    // Get the two image values
    if(gimov->IsInsideBuffer(q))
      {
      // Update the objective
      double I_ref = giref->EvaluateAtContinuousIndex(qr);
      double I_mov = gimov->EvaluateAtContinuousIndex(q);
      f_img += vmask * (I_ref - I_mov) * (I_ref - I_mov);

      // Update the gradient of the objective
      if(need_gradient) 
        {
        itk::CovariantVector<double, VDim> gidx, gmov = gimov->EvaluateDerivativeAtContinuousIndex(q);
        // imov->TransformLocalVectorToPhysicalVector(gmov, gidx);
        for(size_t d = 0; d < VDim; d++)
          {
          // Sfactor = 2.0 / sigma^2
          graw[d][k] = - sfactor * vmask * (I_ref - I_mov) * gmov[d];

          }
        }
      }
    } 

  // Compute the regularization component of the energy
  double f_reg = 0.0;

  for(size_t d = 0; d < VDim; d++)
    { 
    // Take the forward FFT of v[d] into fcmp
    SolvePDE(f_fft_v[d], false);

    // Compute the sum of the norm over all voxels
    double normv = 0.0;
    for(size_t i = 0; i < fft_n; i++)
      normv += vtmp[i] * v[d]->GetBufferPointer()[i];

    // Correct for scaling introduced by the FFT
    normv /= fft_n;

    // Add the norm to the total objective
    f_reg += normv;

    // If gradient is needed, solve the PDE with the previously computed gradient
    if(need_gradient)
      {
      // Take the gradient vector in L2 space and map it to Sobolev space
      SolvePDE(f_fft_g[d], true);

      // Update the gradient as g = 2 v + K (g_int)
      for(size_t i = 0; i < fft_n; i++)
        graw[d][i] = 2.0 * vraw[d][i] + vtmp[i] / fft_n;

      // Update the gradient as g = 2 L L v - g_int
      // for(size_t i = 0; i < fft_n; i++)
      //  graw[d][i] = 2.0 * vtmp[i] / fft_n + graw[d][i];


      }
    }

  double f_total = f_img / (sigma * sigma) + f_reg;

  printf("RMSE(mask) = %8.4f    VMAG(vol) = %8.4f    TOTAL(vol) = %8.4f \n",
    sqrt(f_img * mask_scale), 
    sqrt(f_reg / fft_n),
    f_total / fft_n);

  return f_total;
}

template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::RandomVectorField()
{
  vnl_random randy;

  // Fill the vectors
  for(size_t d = 0; d < VDim; d++)
    {
    // Create index
    int index[VDim];
    for(size_t q=0;q<VDim;q++)
      index[q]=0;

    for(size_t i = 0; i < fft_nc; i++)
      {
      bool lowfrq = true;
      for(size_t q = 0; q < VDim; q++)
        if(index[q] > 5) lowfrq = false;
      if(lowfrq)
        {
        fcmp[i][0] = randy.drand32(-1.0, 1.0);
        fcmp[i][1] = randy.drand32(-1.0, 1.0);
        }
      else
        {
        fcmp[i][0] = 0.0;
        fcmp[i][1] = 0.0;
        }

      // Increment index
      for(size_t d = VDim-1; d >= 0; d--)
        {
        if(++index[d] == fft_szc[d])
          index[d] = 0;
        else
          break;
        }
      }

    myfftw_execute(i_fft);
    for(size_t i = 0; i < fft_n; i++)
      vraw[d][i] = 4 * vtmp[i] / 125;
    }

  typename itk::ImageFileWriter<FloatImage>::Pointer writer = itk::ImageFileWriter<FloatImage>::New();
  writer->SetInput(v[0]);
  writer->SetFileName("dfield.nii");
  writer->Update();
}

template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::TestGradientComputation()
{
  // Create a random vector field
  RandomVectorField();

  // Evaluate the function and the gradient
  double fc = this->ComputeObjectiveAndGradient(true);

  // Random number generator
  vnl_random randy;
  double eps = 1.0e-4;

  // Create a random offset in V and compute directional derivative in the
  // direction of the offset
  VecType *offraw[VDim];
  for(size_t d = 0; d < VDim; d++)
    offraw[d] = (VecType *) myfftw_malloc(sizeof(VecType) * fft_n);

  // Compute forward and backward finite differences
  for(size_t k = 0; k < 10; k++)
    {
    double dirderiv = 0.0;
    for(size_t d = 0; d < VDim; d++)
      {
      for(size_t i = 0; i < fft_n; i++)
        {
        offraw[d][i] = randy.drand32(-1.0, 1.0);
        dirderiv += graw[d][i] * offraw[d][i];
        }

      // Now solve the PDE, i.e., replace offset in L2 space with offset in Sobolev
      myfftw_plan poff = myfftw_plan_dft_r2c(VDim, fft_sz, offraw[d], fcmp, FFTW_ESTIMATE);
      SolvePDE(poff, true);
      myfftw_destroy_plan(poff);
      for(size_t i = 0; i < fft_n; i++)
        offraw[d][i] = vtmp[i] / fft_n;
      }

    double f[2], doff[] = {1.0, -1.0};
    for(size_t q = 0; q < 2; q++)
      {
      // Apply offset to v (in Sobolev space)
      for(size_t d = 0; d < VDim; d++)
        for(size_t i = 0; i < fft_n; i++)
          vraw[d][i] += doff[q] * eps * offraw[d][i];

      // Compute offset
      f[q] = ComputeObjectiveAndGradient(false);

      // Unapply offset to v (in Sobolev space)
      for(size_t d = 0; d < VDim; d++)
        for(size_t i = 0; i < fft_n; i++)
          vraw[d][i] -= doff[q] * eps * offraw[d][i];
      }

    // Compute the deriv
    double cdiffder = (f[0] - f[1]) / (2 * eps);

    cout << "FUNCTION " << fc << "\t";
    cout << "ANALYTIC " << dirderiv << "\t";
    cout << "NUMERIC  " << cdiffder << endl;
    }
}

template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::TestGradient2()
{
  vnl_random randy;

  // Create a random array X
  double *X = new double[VDim * fft_n];
  for(size_t i = 0; i < VDim * fft_n; i++)
    X[i] = randy.drand32(-40.0, 40.0);

  // Create an array to store the gradient
  double *G = new double[VDim * fft_n];
  double *O = new double[VDim * fft_n];

  // Call the same function the optimizer calls
  double fc = cg_valgrad<TPixel,VDim>(G, X, VDim * fft_n);
  cout << "TestGradient2 Solution " << fc << endl;
  SaveRaw(vraw[0], "gradient2field.nii");

  for(size_t q = 0; q < 10; q++)
    {
    // Create a random offset vector
    for(size_t i = 0; i < VDim * fft_n; i++)
      O[i] = randy.drand32(-1.0, 1.0);

    // Compute the variational derivative in the offset direction
    double ad= 0.0;
    for(size_t i = 0; i < VDim * fft_n; i++)
      ad += O[i] * G[i];

    // Compute central difference approximation
    double f[2], dx[] = {-1.0, 1.0}, eps = 1.0e-4;
    for(size_t s = 0; s < 2; s++)
      {
      // add the offset
      for(size_t i = 0; i < VDim * fft_n; i++)
        X[i] += dx[s] * eps * O[i];

      // compute objective
      f[s] = cg_func<TPixel,VDim>(X, VDim * fft_n);

      // subtract the offset
      for(size_t i = 0; i < VDim * fft_n; i++)
        X[i] -= dx[s] * eps * O[i];
      }

    double nd = (f[1] - f[0]) / (2 * eps);
    printf("ANALYTIC %12.8f \t NUMERIC %12.8f \t DIFFERENCE %12.8f\n", ad, nd, ad-nd);
    }

  delete X;
  delete O;
  delete G;
}

template <class TPixel, unsigned int VDim>
void
SimpleElasticRegistration<TPixel, VDim>
::TestFourier()
{
  // Create a plan for the reference image
  myfftw_plan plan = myfftw_plan_dft_r2c(VDim, fft_sz, iref->GetBufferPointer(), fcmp, FFTW_ESTIMATE);
  SolvePDE(plan, true);
  for(size_t i = 0; i < fft_n; i++) vtmp[i] /= fft_n;
  SaveRaw(vtmp, "iref_pde_forward.nii");
  SolvePDE(plan, false);
  for(size_t i = 0; i < fft_n; i++) vtmp[i] /= fft_n;
  SaveRaw(vtmp, "iref_pde_inverse.nii");
}


// Invocations
template class SimpleElasticRegistration<double, 2>;
template class SimpleElasticRegistration<double, 3>;
template class SimpleElasticRegistration<double, 4>;
