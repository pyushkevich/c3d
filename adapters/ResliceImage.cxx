#include "ResliceImage.h"
#include <string>
#include <iostream>
#include <itkAffineTransform.h>
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include "itkResampleImageFilter.h"
#include "gsGSAffine3DTransform.h"

#include <itkTransformFactory.h>

template<unsigned int VDim>
void ReadMatrix(const char *fname, itk::Matrix<double,VDim+1,VDim+1> &mat)
  {
  ifstream fin(fname);
  for(size_t i = 0; i < VDim+1; i++)
    for(size_t j = 0; j < VDim+1; j++)
      if(fin.good())
        {
        fin >> mat[i][j];
        }
      else
        {
        throw ConvertException("Unable to read matrix %s", fname);
        }
  fin.close();
  }

template <class TPixel, unsigned int VDim>
void
ResliceImage<TPixel, VDim>
::operator() (string format, string fn_tran)
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Reslice operation requires two images on the stack" << endl;
    throw -1;
    }

  // Get the reference and source images
  ImagePointer iref = c->m_ImageStack[c->m_ImageStack.size() - 2];
  ImagePointer isrc = c->m_ImageStack.back();

  // Create an initial identity transform
  typedef itk::AffineTransform<double, VDim> TranType;
  typename TranType::Pointer atran = TranType::New();
  atran->SetIdentity();

  // Read transform based on type
  if(format=="itk")
    {
    typedef itk::MatrixOffsetTransformBase<double, VDim, VDim> MOTBType;
    typedef itk::AffineTransform<double, VDim> AffTran;
    itk::TransformFactory<MOTBType>::RegisterTransform();
    itk::TransformFactory<AffTran>::RegisterTransform();

    itk::TransformFileReader::Pointer fltReader = itk::TransformFileReader::New();
    fltReader->SetFileName(fn_tran);
    fltReader->Update();

    itk::TransformBase *base = fltReader->GetTransformList()->front();
    MOTBType *motb = dynamic_cast<MOTBType *>(base);

    if(motb)
      {
      atran->SetMatrix(motb->GetMatrix());
      atran->SetOffset(motb->GetOffset());
      }
    }
  else if(format=="matrix")
    {
    // Read the matrix
    itk::Matrix<double,VDim+1,VDim+1> matrix;
    itk::Matrix<double,VDim,VDim> amat;
    itk::Vector<double, VDim> aoff;

    ReadMatrix<VDim>(fn_tran.c_str(), matrix);
    amat.GetVnlMatrix().update(
      matrix.GetVnlMatrix().extract(VDim, VDim));
    aoff.GetVnlVector().update(
      matrix.GetVnlMatrix().get_column(VDim).extract(VDim));

    // Extrernal matrices are assumed to be RAS to RAS, so we must convert to LPS to LPS
    vnl_vector<double> v_lps_to_ras(VDim, 1.0);
    v_lps_to_ras[0] = v_lps_to_ras[1] = -1.0;
    vnl_diag_matrix<double> m_lps_to_ras(v_lps_to_ras);
    vnl_matrix<double> mold = amat.GetVnlMatrix();
    amat.GetVnlMatrix().update(m_lps_to_ras * mold * m_lps_to_ras);
    aoff.GetVnlVector().update(m_lps_to_ras * aoff.GetVnlVector());

    // Put the values in the transform
    atran->SetMatrix(amat);
    atran->SetOffset(aoff);
    }

  // Build the resampling filter
  typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleFilterType;
  typename ResampleFilterType::Pointer fltSample = ResampleFilterType::New();

  // Initialize the resampling filter with an identity transform
  fltSample->SetInput(isrc);
  fltSample->SetTransform(atran);
  
  // Set the unknown intensity to positive value
  fltSample->SetDefaultPixelValue(c->m_Background);

  // Set the interpolator
  fltSample->SetInterpolator(c->GetInterpolator());

  // Calculate where the transform is taking things
  itk::ContinuousIndex<double, VDim> idx[3];
  for(size_t i = 0; i < VDim; i++)
    {
    idx[0][i] = 0.0;
    idx[1][i] = iref->GetBufferedRegion().GetSize(i) / 2.0;
    idx[2][i] = iref->GetBufferedRegion().GetSize(i) - 1.0;
    }
  for(size_t j = 0; j < VDim; j++)
    {
    itk::ContinuousIndex<double, VDim> idxmov;
    itk::Point<double, VDim> pref, pmov;
    iref->TransformContinuousIndexToPhysicalPoint(idx[j], pref);
    pmov = atran->TransformPoint(pref);
    isrc->TransformPhysicalPointToContinuousIndex(pmov, idxmov);
    *c->verbose << "  Reference voxel " << idx[j] << " => moving voxel " << idxmov << endl;
    }

  // Describe what we are doing
  *c->verbose << "Reslicing #" << c->m_ImageStack.size() 
    << " with reference" << c->m_ImageStack.size() - 1 << endl;
  *c->verbose << "  Interpolation method: " << c->m_Interpolation << endl;
  *c->verbose << "  Background intensity: " << c->m_Background << endl;
  *c->verbose << "  Affine Transform: " << endl;
  vnl_matrix_fixed<double, VDim+1, VDim+1> amat(0.0);
  vnl_vector_fixed<double, VDim+1> atmp(1.0);
  amat.update(atran->GetMatrix().GetVnlMatrix(), 0, 0);
  atmp.update(atran->GetOffset().GetVnlVector(), 0);
  amat.set_column(VDim, atmp);
  c->PrintMatrix(*c->verbose, amat, "%12.5f ", "    ");

  // Set the spacing, origin, direction of the output
  fltSample->UseReferenceImageOn();
  fltSample->SetReferenceImage(iref);
  fltSample->Update();
    
  // Change the source to the output 
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltSample->GetOutput());
}

// Invocations
template class ResliceImage<double, 2>;
template class ResliceImage<double, 3>;
