#ifndef __ConvertAdapter_h_
#define __ConvertAdapter_h_

#include "ConvertImageND.h"

// Common typedefs for all child classes
#define CONVERTER_STANDARD_TYPEDEFS \
  typedef ImageConverter<TPixel, VDim> Converter; \
  typedef typename Converter::ImageType ImageType; \
  typedef typename Converter::UnorientedImageType UnorientedImageType; \
  typedef typename Converter::ImagePointer ImagePointer; \
  typedef typename Converter::SizeType SizeType; \
  typedef typename Converter::IndexType IndexType; \
  typedef typename Converter::IntegerVector IntegerVector; \
  typedef typename Converter::RealVector RealVector; \
  typedef typename Converter::RegionType RegionType; \
  typedef typename Converter::ComplexPixel ComplexPixel; \
  typedef typename Converter::ComplexImageType ComplexImageType; \
  typedef typename Converter::Iterator Iterator; \
  typedef typename Converter::ConstIterator ConstIterator;

/**
 * Parent class for all adapters
 */
template<class TPixel, unsigned int VDim>
class ConvertAdapter
{
public:

  CONVERTER_STANDARD_TYPEDEFS

  ConvertAdapter() {}

protected:

};

#endif
