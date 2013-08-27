#ifndef __ResliceImage_h_
#define __ResliceImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ResliceImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ResliceImage(Converter *c) : c(c) {}

  void operator() (string format, string fn_tran);

private:
  Converter *c;

};

#endif

