#include "CopyTransform.h"

template <class TPixel, unsigned int VDim>
void
CopyTransform<TPixel, VDim>
::operator() ()
{
  // Pop off the last image
  size_t n = c->m_ImageStack.size();
  if(n < 2)
    {
    throw string("Two images must be on the stack");
    }

  // Get data source and transform source
  ImagePointer dsrc = c->m_ImageStack[n-1];
  ImagePointer tsrc = c->m_ImageStack[n-2];

  // Make sure dimensions match
  if(dsrc->GetBufferedRegion().GetSize() != tsrc->GetBufferedRegion().GetSize())
    {
    throw string("Dimensions of images must match");
    }

  // Print out what is being done
  *c->verbose << "Copying transform from #" << n-2 << " to #" << n-1 << endl;

  // Update the data source
  dsrc->SetOrigin(tsrc->GetOrigin());
  dsrc->SetSpacing(tsrc->GetSpacing());
  dsrc->SetDirection(tsrc->GetDirection());

  // Push and pop
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(dsrc);
}

// Invocations
template class CopyTransform<double, 2>;
template class CopyTransform<double, 3>;
