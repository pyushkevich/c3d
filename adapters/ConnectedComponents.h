#ifndef __ConnectedComponents_h_
#define __ConnectedComponents_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ConnectedComponents : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ConnectedComponents(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

