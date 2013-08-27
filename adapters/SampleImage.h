#ifndef __SampleImage_h_
#define __SampleImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class SampleImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  SampleImage(Converter *c) : c(c) {}

  void operator() (const RealVector &x);

  // After the adapter executes, we can obtain the result
  itkGetMacro(Result, double);

private:
  Converter *c;

  double m_Result;
};

#endif

