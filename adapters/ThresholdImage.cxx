#include "ThresholdImage.h"
#include "itkBinaryThresholdImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ThresholdImage<TPixel, VDim>
::operator() (double u1, double u2, double vIn, double vOut)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();
  
  // Say what we are doing
  *c->verbose << "Thresholding #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Mapping range [" << u1 << ", " << u2 << "] to " << vIn << endl;
  *c->verbose << "  Values outside are mapped to " << vOut << endl;

  // Do the thresholding
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(input);
  filter->SetLowerThreshold(u1);
  filter->SetUpperThreshold(u2);
  filter->SetInsideValue(vIn);
  filter->SetOutsideValue(vOut);
  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class ThresholdImage<double, 2>;
template class ThresholdImage<double, 3>;
