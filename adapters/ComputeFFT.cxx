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
