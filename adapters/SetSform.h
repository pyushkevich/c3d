#ifndef __SetSform_h_
#define __SetSform_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class SetSform : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  SetSform(Converter *c) : c(c) {}

  void operator() (string fn_tran);

private:
  Converter *c;

};

#endif

