#ifndef __SplitMultilabelImage_h_
#define __SplitMultilabelImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class SplitMultilabelImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  SplitMultilabelImage(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

