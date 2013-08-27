#ifndef __ReadImage_h_
#define __ReadImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ReadImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ReadImage(Converter *c) : c(c) {}

  void operator() (const char *file);

private:
  Converter *c;

};

#endif

