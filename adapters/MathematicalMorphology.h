#ifndef __MathematicalMorphology_h_
#define __MathematicalMorphology_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class MathematicalMorphology : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  MathematicalMorphology(Converter *c) : c(c) {}

  void operator() (bool erode, TPixel value, SizeType size);

private:
  Converter *c;

};

#endif

