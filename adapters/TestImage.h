#ifndef __TestImage_h_
#define __TestImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class TestImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  TestImage(Converter *c) : c(c) {}

  void operator() (bool test_header, bool test_image, double tol);

private:
  Converter *c;

};

#endif

