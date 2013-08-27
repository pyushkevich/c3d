#include "PeronaMalik.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

template <class TPixel, unsigned int VDim>
void
PeronaMalik<TPixel, VDim>
::operator() (double conductance, size_t nIter)
{
  // Get an image off the stack
  ImagePointer image = c->m_ImageStack.back();

  // Create a filter
  typedef itk::GradientAnisotropicDiffusionImageFilter<
    ImageType,ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  *c->verbose << "Performing anisotropic diffusion on #" << c->m_ImageStack.size() << endl;
  filter->SetInput(image);
  filter->SetConductanceParameter(conductance);
  filter->SetNumberOfIterations(nIter);
  filter->SetTimeStep(0.0125);
  filter->UseImageSpacingOn();
  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class PeronaMalik<double, 2>;
template class PeronaMalik<double, 3>;
