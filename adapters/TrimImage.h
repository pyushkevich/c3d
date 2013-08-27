#ifndef __TrimImage_h_
#define __TrimImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class TrimImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  TrimImage(Converter *c) : c(c) {}

  enum TrimMode {
    SPECIFY_MARGIN,
    SPECIFY_FINALSIZE };

  void operator() (const RealVector &vec, TrimMode mode);

private:
  Converter *c;

};

#endif

