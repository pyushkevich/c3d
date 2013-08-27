#include "AddImages.h"
#include "itkAddImageFilter.h"

template <class TPixel, unsigned int VDim>
void
AddImages<TPixel, VDim>
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
  *c->verbose << "Adding #" << c->m_ImageStack.size() - 1 << " and "  
    << c->m_ImageStack.size() - 2 << endl;

  // Perform the multiplication
  typedef itk::AddImageFilter<ImageType, ImageType, ImageType> FilterType;
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
template class AddImages<double, 2>;
template class AddImages<double, 3>;
