#ifndef __TileImages_h_
#define __TileImages_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class TileImages : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  TileImages(Converter *c) : c(c) {}

  void operator() (const std::string &tileParam);

private:
  Converter *c;

};

#endif

