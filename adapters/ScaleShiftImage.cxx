#include "ScaleShiftImage.h"
#include "itkShiftScaleImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ScaleShiftImage<TPixel, VDim>
::operator() (double a, double b)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();
  
  // Say what we are doing
  *c->verbose << "Scaling #" << c->m_ImageStack.size() << " by " << a << " and adding " << b << endl;

  // If a=0, this means setting the image to a constant
  if(a == 0.0)
    {
    c->CopyImage();
    c->m_ImageStack.back()->FillBuffer(b);
    return;
    }

  // Create and run filter
  typedef itk::ShiftScaleImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(input);
  filter->SetScale(a);
  filter->SetShift(b / a);
  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class ScaleShiftImage<double, 2>;
template class ScaleShiftImage<double, 3>;
