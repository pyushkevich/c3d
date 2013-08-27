#ifndef __WeightedSum_h_
#define __WeightedSum_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class WeightedSum : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  WeightedSum(Converter *c) : c(c) {}

  void operator() (std::vector<double> weights);

private:
  Converter *c;

};

#endif

