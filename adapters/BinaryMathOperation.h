#ifndef __BinaryMathOperation_h_
#define __BinaryMathOperation_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class BinaryMathOperation : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  // Operations supported by this class
  enum Operation 
    {
    ADD = 0,
    SUBTRACT,
    MULTIPLY,
    DIVIDE,
    MAXIMUM,
    MINIMUM 
    };

  BinaryMathOperation(Converter *c) : c(c) {}

  void operator() (Operation op);

private:
  Converter *c;

};

#endif

