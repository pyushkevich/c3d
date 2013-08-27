#ifndef __VoxelwiseRegression_h_
#define __VoxelwiseRegression_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class VoxelwiseRegression : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  VoxelwiseRegression(Converter *c) : c(c) {}

  void operator() (size_t order);

private:
  Converter *c;

};

#endif

