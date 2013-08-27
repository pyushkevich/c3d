#ifndef __CreateInterpolator_h_
#define __CreateInterpolator_h_

#include "ConvertAdapter.h"
#include "itkSmartPointer.h"

namespace itk {
  template <class TImage, class TCoordRep> class InterpolateImageFunction;
}

template<class TPixel, unsigned int VDim>
class CreateInterpolator : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  CreateInterpolator (Converter *c) : c(c) {}

  // Helper function: get interpolator based on current flag values
  typedef itk::InterpolateImageFunction<ImageType, double> InterpolatorType;

  void CreateNN();
  void CreateLinear();
  void CreateCubic();
  void CreateSinc();
  void CreateGaussian(RealVector sigma);
  void CreateMultiLabel(RealVector sigma);

private:
  itk::SmartPointer<InterpolatorType> m_Interp;
  Converter *c;
};

#endif

