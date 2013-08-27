#include "StapleAlgorithm.h"
#include "itkSTAPLEImageFilter.h"

template <class TPixel, unsigned int VDim>
void
StapleAlgorithm<TPixel, VDim>
::operator() (double ival)
{
  size_t i;
  
  // Create a STAPLE filter
  typedef itk::STAPLEImageFilter<ImageType, ImageType> Filter;
  typename Filter::Pointer fltStaple = Filter::New();

  // Add each of the images on the stack to the filter
  for(i = 0; i < c->m_ImageStack.size(); i++)
    fltStaple->SetInput(i, c->m_ImageStack[i]);

  // Configure the STAPLE filter
  fltStaple->SetForegroundValue(ival);

  // Describe what we are doing
  *c->verbose << "Executing STAPLE EM Algorithm on " << c->m_ImageStack.size() << " images." << endl;

  // Plug in the filter's components
  fltStaple->Update();

  // Dump sensitivity/specificity values
  *c->verbose << "  Elapsed Iterations: " << fltStaple->GetElapsedIterations() << endl;
  for(i = 0; i < c->m_ImageStack.size(); i++)
    {
    *c->verbose << "  Rater " << i 
      << ": Sensitivity = " << fltStaple->GetSensitivity(i) 
      << "; Specificity = " << fltStaple->GetSpecificity(i) << endl;
    }

  // Store the output
  c->m_ImageStack.clear();
  c->m_ImageStack.push_back(fltStaple->GetOutput());

}

// Invocations
template class StapleAlgorithm<double, 2>;
template class StapleAlgorithm<double, 3>;
