#ifndef __GeneralLinearModel_h_
#define __GeneralLinearModel_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class GeneralLinearModel : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  GeneralLinearModel(Converter *c) : c(c) {}

  void operator() (string fn_matrix, string fn_contrast);

private:
  Converter *c;

};

#endif

