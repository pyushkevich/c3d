#include "ReciprocalImage.h"
#include "itkUnaryFunctorImageFilter.h"

template <class TIn, class TOut>
class ReciprocalFunctor
{
public:
  TOut operator() (const TIn &a)
    { return (TOut) (1.0 / a); }
};

template <class TPixel, unsigned int VDim>
void
ReciprocalImage<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Verbose
  *c->verbose << "Taking the reciprocal of #" << c->m_ImageStack.size() << endl;

  // Simply go through and divide
  typedef ReciprocalFunctor<TPixel, TPixel> FunctorType;
  typedef itk::UnaryFunctorImageFilter<ImageType, ImageType, FunctorType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(img);
  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class ReciprocalImage<double, 2>;
template class ReciprocalImage<double, 3>;
