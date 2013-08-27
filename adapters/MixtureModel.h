#ifndef __MixtureModel_h_
#define __MixtureModel_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class MixtureModel : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  MixtureModel(Converter *c) : c(c) {}

  void operator() (std::vector<double> mu, std::vector<double> sigma);

private:
  Converter *c;

};

#endif

