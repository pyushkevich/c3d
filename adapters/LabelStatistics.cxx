#include "LabelStatistics.h"
#include "itkLabelStatisticsImageFilter.h"
#include <set>

template <class TPixel, unsigned int VDim>
void
LabelStatistics<TPixel, VDim>
::operator() ()
{
  // Get images from stack
  size_t n = c->m_ImageStack.size();
  if(n < 2)
    throw ConvertException("Label statistics require two image inputs");
  ImagePointer label = c->m_ImageStack[n-1];
  ImagePointer image = c->m_ImageStack[n-2];

  // Create a short image for the labels
  typedef itk::Image<short, VDim> LabelImageType;
  typename LabelImageType::Pointer slab = LabelImageType::New();

  // Allocate the image
  slab->CopyInformation(label);
  slab->SetRegions(label->GetBufferedRegion());
  slab->Allocate();

  // Accumulator of label values    
  std::set<short> sval;

  // Round off doubles to create labels
  size_t nv = label->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < nv; i++)
    {
    slab->GetBufferPointer()[i] = (short) (label->GetBufferPointer()[i] + 0.5);
    sval.insert(slab->GetBufferPointer()[i]);
    }

  // Create the label statistics filter
  typedef itk::LabelStatisticsImageFilter<ImageType, LabelImageType> StatFilter;
  typename StatFilter::Pointer fltStat = StatFilter::New();
  
  // Set the inputs
  fltStat->SetInput(image);
  fltStat->SetLabelInput(slab);

  // Update the filter
  fltStat->Update();

  // Get the voxel dimensions
  double dim = 1.0;
  for(size_t i = 0; i < VDim; i++)
    dim *= label->GetSpacing()[i];

  // Get the number of labels                                                             .
  printf("LabelID        Mean        StdD         Max         Min       Count     Vol(mm^%1d)        Extent(Vox)\n", VDim);
  for(set<short>::iterator it = sval.begin(); it != sval.end(); ++it)
    {
    // printf("xxxxx    xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx  xxxxxxxxxx");
    printf("%5i    %10.5f  %10.5f  %10.5f  %10.5f  %10lu  %12.3f ",
      (int) *it, 
      fltStat->GetMean(*it), 
      fltStat->GetSigma(*it), 
      fltStat->GetMaximum(*it), 
      fltStat->GetMinimum(*it), 
      (long unsigned) fltStat->GetCount(*it),
      fltStat->GetCount(*it) * dim);

    // Bounding box
    typename StatFilter::BoundingBoxType bbox = fltStat->GetBoundingBox(*it);
    for(size_t i = 0; i < VDim; i++)
      {
      printf(" %5lu", 1 + bbox[i*2+1] - bbox[i*2]);
      }
    printf("\n");
    }
}

// Invocations
template class LabelStatistics<double, 2>;
template class LabelStatistics<double, 3>;

