#ifndef __UpdateMetadataKey_h_
#define __UpdateMetadataKey_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class UpdateMetadataKey : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  UpdateMetadataKey(Converter *c) : c(c) {}

  void operator() (const char *key, const char *value);

private:
  Converter *c;

};

#endif

