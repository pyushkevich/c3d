/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ComputeFFT.cxx
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

#include "ComputeFFT.h"
#include "itkRealToHalfHermitianForwardFFTImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkComplexToImaginaryImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ComputeFFT<TPixel, VDim>
::operator() ()
{
  // Compute the fourier transform of the image. This will take one image
  // off the stack and place two images (real, imag) on the stack
  ImagePointer image = c->m_ImageStack.back();
  
  // Create the fourier transform filter
  typedef itk::Image<ComplexPixel, VDim> UnorientedComplexImageType;
  typedef itk::RealToHalfHermitianForwardFFTImageFilter<ImageType> FFTFilter;
  typename FFTFilter::Pointer fltFourier = FFTFilter::New();

  // Get the real and imaginary components
  typedef itk::ComplexToRealImageFilter<UnorientedComplexImageType,ImageType> RealFilter;
  typedef itk::ComplexToImaginaryImageFilter<UnorientedComplexImageType,ImageType> ImagFilter;
  typename RealFilter::Pointer fltReal = RealFilter::New();
  typename ImagFilter::Pointer fltImag = ImagFilter::New();

  // Set inputs and outputs
  cout << "DOING FFT" << endl;
  try 
    {
    fltFourier->SetInput(image);
    fltFourier->Update();
    }
  catch(itk::ExceptionObject &exc) 
    {
    cerr << "Exception caught : " << exc;
    cerr << endl;
    return;
    }
  cout << "DID MAIN PART" << endl;

  fltReal->SetInput(fltFourier->GetOutput());
  fltImag->SetInput(fltFourier->GetOutput());

  // Compute the transforms
  fltReal->Update();
  fltImag->Update();
  cout << "DID FFT" << endl;

  // Pop the last guy from the stack and push the real and imag
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltReal->GetOutput());
  c->m_ImageStack.push_back(fltImag->GetOutput());
  cout << "FINISHED STACK" << endl;

}

// Invocations
template class ComputeFFT<double, 2>;
template class ComputeFFT<double, 3>;
template class ComputeFFT<double, 4>;
