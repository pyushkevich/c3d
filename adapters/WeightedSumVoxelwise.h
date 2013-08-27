#ifndef __WeightedSumVoxelwise_h_
#define __WeightedSumVoxelwise_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class WeightedSumVoxelwise : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  WeightedSumVoxelwise(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

