#ifndef __BinaryImageCentroid_h_
#define __BinaryImageCentroid_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class BinaryImageCentroid : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  BinaryImageCentroid(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

