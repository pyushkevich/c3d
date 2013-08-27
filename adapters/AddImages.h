#ifndef __AddImages_h_
#define __AddImages_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class AddImages : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  AddImages(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

