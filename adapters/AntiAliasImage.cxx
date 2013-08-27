#include "AntiAliasImage.h"
#include "itkAntiAliasBinaryImageFilter.h"

template <class TPixel, unsigned int VDim>
void
AntiAliasImage<TPixel, VDim>
::operator() (double xIsoSurface)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Report what the filter is doing
  *c->verbose << "Anti-aliasing #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Root Mean Square error: " << c->m_AntiAliasRMS << endl;
  *c->verbose << "  Iterations: "; 
  if(c->m_Iterations == 0) 
    *c->verbose << "Unlimited" << endl; 
  else 
    *c->verbose << c->m_Iterations << endl;

  // Apply antialiasing to the image
  typedef itk::AntiAliasBinaryImageFilter<ImageType,ImageType> AntiFilterType;
  typename AntiFilterType::Pointer fltAnti = AntiFilterType::New();
  fltAnti->SetInput(input);
  fltAnti->SetMaximumRMSError(c->m_AntiAliasRMS);
  if(c->m_Iterations > 0)
    fltAnti->SetNumberOfIterations(c->m_Iterations);
  fltAnti->SetIsoSurfaceValue(xIsoSurface);
  // fltAnti->AddObserver(itk::ProgressEvent(),command);
  fltAnti->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltAnti->GetOutput());
}

// Invocations
template class AntiAliasImage<double, 2>;
template class AntiAliasImage<double, 3>;
