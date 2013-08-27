#include "ExtractRegion.h"
#include "itkRegionOfInterestImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ExtractRegion<TPixel, VDim>
::operator() (RegionType bbox)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Make sure the bounding box is within the contents of the image
  bbox.Crop(input->GetBufferedRegion());

  // Report the bounding box size
  *c->verbose << "  Extracting bounding box " << bbox.GetIndex() << " " << bbox.GetSize() << endl;

  // Chop off the region
  typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> TrimFilter;
  typename TrimFilter::Pointer fltTrim = TrimFilter::New();
  fltTrim->SetInput(input);
  fltTrim->SetRegionOfInterest(bbox);
  fltTrim->Update();

  // What happened to the origin of the image?
  ImagePointer output = fltTrim->GetOutput();

  // Update the image stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(output);
}

// Invocations
template class ExtractRegion<double, 2>;
template class ExtractRegion<double, 3>;
