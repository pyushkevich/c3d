#ifndef __ScalarToRGB_h_
#define __ScalarToRGB_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ScalarToRGB : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ScalarToRGB(Converter *c) : c(c) {}

  void operator() (const std::string &colormap);

private:
  Converter *c;

};

#endif

