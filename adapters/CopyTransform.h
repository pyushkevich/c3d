#ifndef __CopyTransform_h_
#define __CopyTransform_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class CopyTransform : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  CopyTransform(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

