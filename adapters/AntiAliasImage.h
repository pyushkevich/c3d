#ifndef __AntiAliasImage_h_
#define __AntiAliasImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class AntiAliasImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  AntiAliasImage(Converter *c) : c(c) {}

  void operator() (double xIsoSurface);

private:
  Converter *c;

};

#endif

