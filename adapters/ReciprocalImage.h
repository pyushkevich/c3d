#ifndef __ReciprocalImage_h_
#define __ReciprocalImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ReciprocalImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ReciprocalImage(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

