#ifndef __LandmarksToSpheres_h_
#define __LandmarksToSpheres_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class LandmarksToSpheres : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  LandmarksToSpheres(Converter *c) : c(c) {}

  void operator() (const char *fnland, double radius);

private:
  Converter *c;

};

#endif

