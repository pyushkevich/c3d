#include "MultiplyImages.h"
#include "itkMultiplyImageFilter.h"

template <class TPixel, unsigned int VDim>
void
MultiplyImages<TPixel, VDim>
::operator() ()
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Binary operations require two images on the stack" << endl;
    throw -1;
    }

  // Get the last two images
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Write something
  *c->verbose << "Multiplying #" << c->m_ImageStack.size() - 1 
    << " by #" << c->m_ImageStack.size() << endl;

  // Perform the multiplication
  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType> FilterType;
  typename FilterType::Pointer flt = FilterType::New();
  flt->SetInput1(i1);
  flt->SetInput2(i2);
  flt->Update();

  // Replace the images with the product
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(flt->GetOutput());

}

// Invocations
template class MultiplyImages<double, 2>;
template class MultiplyImages<double, 3>;
