#include "CoordinateMap.h"

template <class TPixel, unsigned int VDim>
void
CoordinateMap<TPixel, VDim>
::operator() (bool physical)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  *c->verbose << "Replacing #" << c->m_ImageStack.size()
    << " with " << VDim << "coordinate maps" << endl;

  // Create three new images
  ImagePointer coord[VDim];
  Iterator it[VDim];

  for(size_t i = 0; i < VDim; i++)
    {
    coord[i] = ImageType::New();
    coord[i]->SetRegions(img->GetBufferedRegion());
    coord[i]->CopyInformation(img);
    coord[i]->Allocate();
    it[i] = Iterator(coord[i], img->GetBufferedRegion());
    }

  while(!it[0].IsAtEnd())
    {
    IndexType idx = it[0].GetIndex();
    if(physical)
      {
      typename ImageType::PointType p;
      img->TransformIndexToRASPhysicalPoint(idx, p);
      for(size_t i = 0; i < VDim; i++)
        it[i].Set(p[i]);
      }
    else
      {
      for(size_t i = 0; i < VDim; i++)
        it[i].Set(idx[i]);
      }
    for(size_t i = 0; i < VDim; i++)
      ++it[i];
    }

  // Do some processing ...
  // ImagePointer result = ...;
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  for(size_t i = 0; i < VDim; i++)
    c->m_ImageStack.push_back(coord[i]);
}

// Invocations
template class CoordinateMap<double, 2>;
template class CoordinateMap<double, 3>;
