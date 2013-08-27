#ifndef __UnaryMathOperation_h_
#define __UnaryMathOperation_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class UnaryMathOperation : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  UnaryMathOperation(Converter *c) : c(c) {}

  void operator() (double (*func)(double));

private:
  Converter *c;

};

#endif

