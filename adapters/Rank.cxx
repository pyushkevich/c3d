#include "Rank.h"
#include "itkAddImageFilter.h"
#include "itkMetaDataObject.h"


template <class TPixel, unsigned int VDim>
void
Rank<TPixel, VDim>
::operator() ()
{
  // Create a maximum image
  ImagePointer i0 = c->m_ImageStack[0];

  // This image replaces each voxel in the input image
  // with the rank (from 0 to n-1). 

  // Say something
  *c->verbose << "Ranking " << c->m_ImageStack.size() << 
    " images. " << endl;


  // For each of the images, update the vote
  size_t n = c->m_ImageStack.size();
  for(size_t j = 1; j < n; j++)
    {
    // Get the next image pointer
    ImagePointer ij = c->m_ImageStack[j];

    // Check the image dimensions
    if(ij->GetBufferedRegion() != c->m_ImageStack.back()->GetBufferedRegion())
      throw ConvertException("All images must have same dimensions");
    }

  size_t nvox = c->m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
  std::vector< std::pair<double, size_t> > sortvec(n);

  for(size_t k = 0; k < nvox; k++)
    {
    for(size_t j = 0; j < n; j++)
      {
      sortvec[j] = std::pair<double,size_t>(c->m_ImageStack[j]->GetBufferPointer()[k], j);
      }

    if(k == 189552)
      for(size_t j = 0; j < n; j++)
        printf("J = %li, v = %f, r = %li\n", j, sortvec[j].first, sortvec[j].second);
      

    std::sort(sortvec.begin(), sortvec.end());

    for(size_t j = 0; j < n; j++)
      {
      c->m_ImageStack[sortvec[j].second]->GetBufferPointer()[k] = n - j;
      }
    }
}

// Invocations
template class Rank<double, 2>;
template class Rank<double, 3>;
