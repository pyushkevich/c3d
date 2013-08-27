#include "HistogramMatch.h"
#include "itkHistogramMatchingImageFilter.h"

template <class TPixel, unsigned int VDim>
void
HistogramMatch<TPixel, VDim>
::operator() (int nmatch)
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Reslice operation requires two images on the stack" << endl;
    throw -1;
    }

  // Get the reference and source images
  ImagePointer iref = c->m_ImageStack[c->m_ImageStack.size() - 2];
  ImagePointer isrc = c->m_ImageStack.back();

  // Create the filter
  typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> HistogramFilter;
  typename HistogramFilter::Pointer filter = HistogramFilter::New();

  filter->SetReferenceImage(iref);
  filter->SetSourceImage(isrc);
  filter->SetNumberOfMatchPoints(nmatch);
  filter->ThresholdAtMeanIntensityOff();

  *c->verbose << "Histogram matching #" << c->m_ImageStack.size() 
    << " to reference" << c->m_ImageStack.size() - 1 << endl;
  *c->verbose << "  Number of match points: " << filter->GetNumberOfMatchPoints() << endl;
  *c->verbose << "  Number of histogram levels: " << filter->GetNumberOfHistogramLevels() << endl;

  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class HistogramMatch<double, 2>;
template class HistogramMatch<double, 3>;
