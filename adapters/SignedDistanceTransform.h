#ifndef __SignedDistanceTransform_h_
#define __SignedDistanceTransform_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class SignedDistanceTransform : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  SignedDistanceTransform(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

