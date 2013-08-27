#ifndef __SetOrientation_h_
#define __SetOrientation_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class SetOrientation : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  SetOrientation(Converter *c) : c(c) {}

  void operator() (std::string rai);

private:
  Converter *c;

};

#endif

