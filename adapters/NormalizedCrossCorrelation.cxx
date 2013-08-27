#include "NormalizedCrossCorrelation.h"
#include "itkConstNeighborhoodIterator.h"

template <class TPixel, unsigned int VDim>
void
NormalizedCrossCorrelation<TPixel, VDim>
::operator() (itk::Size<VDim> radius)
{
  // Get two images from stack
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size()-1];
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size()-2];

  // Iterate over both images with a neighborhood iterator
  typedef typename  itk::ConstNeighborhoodIterator<ImageType> Iter;
  Iter q1(radius, i1,i1->GetBufferedRegion());
  Iter q2(radius, i2,i2->GetBufferedRegion());

  typename ImageType::Pointer iout = ImageType::New();
  iout->CopyInformation(i1);
  iout->SetRegions(i1->GetBufferedRegion());
  iout->Allocate();

  Iterator it(iout, iout->GetBufferedRegion());

  for(; !q1.IsAtEnd() && !q2.IsAtEnd(); ++q1, ++q2, ++it)
    {
    // Compute correlation at this location
    double sumx=0.0, ssqx = 0.0; 
    double sumy=0.0, ssqy = 0.0; 
    double sumxy=0.0;
    size_t n = 0;

    for(size_t i = 0; i < q1.Size(); i++)
      {
      bool inb1,inb2;
      double x = q1.GetPixel(i, inb1);
      double y = q2.GetPixel(i, inb2);
      sumx += x;
      sumy += y;
      sumxy += x * y;
      ssqx += x * x;
      ssqy += y * y;
      n++;
      }

    double cc = (sumxy - sumx * sumy / n) / 
      (sqrt((ssqx - sumx * sumx / n) * (ssqy - sumy * sumy / n)));

    it.Set(cc);
    }

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(iout);
}

// Invocations
template class NormalizedCrossCorrelation<double, 2>;
template class NormalizedCrossCorrelation<double, 3>;
