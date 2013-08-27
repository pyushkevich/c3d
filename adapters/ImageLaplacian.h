#ifndef __ImageLaplacian_h_
#define __ImageLaplacian_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ImageLaplacian : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  ImageLaplacian(Converter *c) : c(c) {}

  void operator() ();

private:
  Converter *c;

};

#endif

