#ifndef __BiasFieldCorrection_h_
#define __BiasFieldCorrection_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class BiasFieldCorrection : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  BiasFieldCorrection(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

