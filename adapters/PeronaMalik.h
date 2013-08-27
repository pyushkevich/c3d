#ifndef __PeronaMalik_h_
#define __PeronaMalik_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class PeronaMalik : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  PeronaMalik(Converter *c) : c(c) {}

  void operator() (double conductance, size_t nIter);

private:
  Converter *c;

};

#endif

