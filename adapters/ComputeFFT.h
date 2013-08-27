#ifndef __ComputeFFT_h_
#define __ComputeFFT_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ComputeFFT : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ComputeFFT(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

