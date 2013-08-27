#ifndef __Rank_h_
#define __Rank_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class Rank : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  Rank(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

