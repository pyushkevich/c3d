#ifndef __ReplaceIntensities_h_
#define __ReplaceIntensities_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ReplaceIntensities : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ReplaceIntensities(Converter *c) : c(c) {}

  void operator() (vector<double> &xRule);

private:
  Converter *c;

};

#endif

