#ifndef __LevelSetSegmentation_h_
#define __LevelSetSegmentation_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class LevelSetSegmentation : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  LevelSetSegmentation(Converter *c) : c(c) {}

  void operator() (int nIter);

private:
  Converter *c;

};

#endif

