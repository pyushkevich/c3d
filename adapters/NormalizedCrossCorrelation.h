#ifndef __NormalizedCrossCorrelation_h_
#define __NormalizedCrossCorrelation_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class NormalizedCrossCorrelation : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  NormalizedCrossCorrelation(Converter *c) : c(c) {}

  void operator() (itk::Size<VDim> radius);

private:
  Converter *c;

};

#endif

