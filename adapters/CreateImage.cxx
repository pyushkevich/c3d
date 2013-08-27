#include "CreateImage.h"
#include "itkImage.h"

template <class TPixel, unsigned int VDim>
void
CreateImage<TPixel, VDim>
::operator()(SizeType dims, RealVector voxelSize)
{
  // Create a new region
  RegionType region;
  region.SetSize(dims);

  // Create the image
  ImagePointer img = ImageType::New();
  img->SetRegions(region);
  img->Allocate();
  img->FillBuffer(c->m_Background);
  
  // Set the voxel size
  img->SetSpacing(voxelSize.data_block());

  // Report
  *c->verbose << "Creating #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Dimensions: " << dims << endl;
  *c->verbose << "  Spacing: " << voxelSize << endl;
  
  // Push the image into the stack
  c->m_ImageStack.push_back(img);
}

// Invocations
template class CreateImage<double, 2>;
template class CreateImage<double, 3>;
