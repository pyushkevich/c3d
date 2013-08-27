#ifndef __ResampleImage_h_
#define __ResampleImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ResampleImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ResampleImage(Converter *c) : c(c) {}

  void operator() (SizeType &size);

private:
  Converter *c;

};

#endif

