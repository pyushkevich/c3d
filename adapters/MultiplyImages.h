#ifndef __MultiplyImages_h_
#define __MultiplyImages_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class MultiplyImages : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  MultiplyImages(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

