#ifndef __AlignByLandmarks_h_
#define __AlignByLandmarks_h_

#include "ConvertAdapter.h"
#include <map>
#include "itkPoint.h"

template<class TPixel, unsigned int VDim>
class AlignByLandmarks : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  AlignByLandmarks(Converter *c) : c(c) {}

  void operator() (int dof, std::string fn_output);

private:
  Converter *c;

  typedef std::map<unsigned long, itk::Point<double, VDim> > CentroidMap;

  CentroidMap ExtractCentroids(ImageType *image);

};

#endif

