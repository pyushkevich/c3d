#ifndef __CreateImage_h_
#define __CreateImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class CreateImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  CreateImage(Converter *c) : c(c) {}

  void operator() (SizeType dims, RealVector voxelSize);

private:
  Converter *c;

};

#endif

