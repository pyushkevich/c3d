#ifndef __ImageERF_h_
#define __ImageERF_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ImageERF : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ImageERF(Converter *c) : c(c) {}

  void operator() (double thresh, double scale);

private:
  Converter *c;

};

#endif

