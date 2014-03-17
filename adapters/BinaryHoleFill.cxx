#include "BinaryHoleFill.h"
#include "itkBinaryFillholeImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BinaryHoleFill<TPixel, VDim>
::operator() (double foreground, bool full_conn)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Apply the fill hole algorithm
  typedef itk::BinaryFillholeImageFilter<ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(img);
  filter->SetForegroundValue(foreground);
  filter->SetFullyConnected(full_conn);

  // Some debug message
  *c->verbose << "Performing binary hole fill for intensity value " << foreground 
    << " in # " << c->m_ImageStack.size() << endl;

  filter->Update();
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class BinaryHoleFill<double, 2>;
template class BinaryHoleFill<double, 3>;
