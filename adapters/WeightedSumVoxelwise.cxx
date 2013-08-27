#include "WeightedSumVoxelwise.h"

template <class TPixel, unsigned int VDim>
void
WeightedSumVoxelwise<TPixel, VDim>
::operator() ()
{
  // Check input validity 
  if((c->m_ImageStack.size() % 2) == 1)
    throw ConvertException("Number of images on stack incompatible with"
      " weighted averaging (must be even");

  if(!c->CheckStackSameDimensions(0))
    throw ConvertException("Images on the stack don't have same dimensions");

  // Make a copy of the last image on the stack
  ImagePointer isum = ImageType::New();
  isum->SetRegions(c->m_ImageStack.back()->GetBufferedRegion());
  isum->CopyInformation(c->m_ImageStack.back());
  isum->Allocate();
  isum->FillBuffer(0);

  *c->verbose << "Taking weighted sum of " << c->m_ImageStack.size() / 2 << " images.";

  // Compute the weighted sum, two images at a time
  while(c->m_ImageStack.size())
    {
    // Get images from the stack
    ImagePointer img = c->m_ImageStack.back(); 
    c->m_ImageStack.pop_back();
    ImagePointer wgt = c->m_ImageStack.back(); 
    c->m_ImageStack.pop_back();

    // Compute the weighted sum
    Iterator is(img, img->GetBufferedRegion());
    Iterator iw(wgt, wgt->GetBufferedRegion());
    Iterator it(isum, isum->GetBufferedRegion());
    for(; !is.IsAtEnd(); ++is, ++it, ++iw)
      it.Set(it.Get() + is.Get() * iw.Get());
    }

  // Put result on stack
  c->m_ImageStack.push_back(isum);
}

// Invocations
template class WeightedSumVoxelwise<double, 2>;
template class WeightedSumVoxelwise<double, 3>;
