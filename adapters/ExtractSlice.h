
#ifndef __ExtractSlice_h_
#define __ExtractSlice_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ExtractSlice : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ExtractSlice(Converter *c) : c(c) {}

  void operator() (string axis, char* pos);

private:
  Converter *c;

};

#endif

