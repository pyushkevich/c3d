#ifndef __NormalizeLocalWindow_h_
#define __NormalizeLocalWindow_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class NormalizeLocalWindow : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  NormalizeLocalWindow(Converter *c) : c(c) {}

  void operator() (SizeType radius);

private:
  Converter *c;

};

#endif

