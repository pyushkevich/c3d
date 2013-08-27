#include "BinaryImageCentroid.h"

template <class TPixel, unsigned int VDim>
void
BinaryImageCentroid<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  *c->verbose << "Computing centroid of #" << c->m_ImageStack.size() << endl;

  // Find all pixels that do not match background
  itk::ContinuousIndex<double, VDim> center_idx;
  center_idx.Fill(0.0);
  size_t n = 0;
  for(Iterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(it.Get() != c->m_Background)
      {
      n++;
      for(size_t d = 0; d < VDim; d++)
        center_idx[d] += it.GetIndex()[d];
      }      
    }

  // Find the centroid in voxel units
  for(size_t d = 0; d < VDim; d++)
    center_idx[d] /= n;

  // Find the centroid in Nifti units
  itk::Point<double, VDim> center_pt;
  img->TransformContinuousIndexToRASPhysicalPoint(center_idx, center_pt);

  // Print the resuts
  cout << "CENTROID_VOX " << center_idx << endl;
  cout << "CENTROID_MM " << center_pt << endl;
}

// Invocations
template class BinaryImageCentroid<double, 2>;
template class BinaryImageCentroid<double, 3>;
