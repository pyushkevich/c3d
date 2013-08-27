#ifndef __ApplyMetric_h_
#define __ApplyMetric_h_

#include "ConvertAdapter.h"
#include <itkAffineTransform.h>
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"


template<class TPixel, unsigned int VDim>
class ApplyMetric : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  typedef vnl_matrix_fixed<double, VDim+1, VDim+1> MatrixType;
  typedef typename itk::AffineTransform<double, VDim> TransformType;
  typedef typename TransformType::Pointer TransformPointer;
  // Interpolators
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double> NNInterpolatorType;
  typedef itk::LinearInterpolateImageFunction< ImageType, double > LinInterpolatorType;

  ApplyMetric(Converter *c) : c(c) {}

  void operator() (const char *metric, const char *ftran_fn, const char *mtran_fn);
  void ReadMatrix(const char *fname, itk::Matrix<double,VDim+1,VDim+1> &mat);
  void CreateHalfwayImageSpace( ImagePointer fixed, ImagePointer moving, ImagePointer halfway);
  double GetValueInternalSymmetric( ImagePointer fixed, ImagePointer moving, ImagePointer halfway,
                                  TransformPointer ftran,                              
                                  TransformPointer mtran, const char * metric_name );
  void Flip_LPS_to_RAS(  itk::Matrix<double,VDim+1,VDim+1> &matrix,
                      itk::Matrix<double,VDim,VDim> &amat,
                      itk::Vector<double, VDim> &aoff);
  void Flip_RAS_to_LPS(  itk::Matrix<double,VDim+1,VDim+1> &matrix,
                      itk::Matrix<double,VDim,VDim> &amat,
                      itk::Vector<double, VDim> &aoff);




private:
  Converter *c;

};

#endif

