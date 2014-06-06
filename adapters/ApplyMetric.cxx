/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ApplyMetric.cxx
  Language:  C++
  Website:   itksnap.org/c3d
  Copyright (c) 2014 Paul A. Yushkevich
  
  This file is part of C3D, a command-line companion tool to ITK-SNAP

  C3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#include <string>
#include <iostream>
#include "ApplyMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMutualInformationHistogramImageToImageMetric.h"
#include "itkNormalizedMutualInformationHistogramImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkAffineTransform.h"
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include "itkResampleImageFilter.h"
#include "gsGSAffine3DTransform.h"
#include <vnl/vnl_inverse.h>
#include <vnl/algo/vnl_matrix_inverse.h>


template <class TPixel, unsigned int VDim>
void
ApplyMetric<TPixel, VDim>
::ReadMatrix(const char *fname, itk::Matrix<double,VDim+1,VDim+1> &mat)
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
ApplyMetric<TPixel, VDim>
::Flip_LPS_to_RAS(  itk::Matrix<double,VDim+1,VDim+1> &matrix,
                    itk::Matrix<double,VDim,VDim> &amat,
                    itk::Vector<double, VDim> &aoff) 
{   

    // Convert lps to ras
    vnl_vector<double> v_ras_to_lps(VDim, 1.0);
    v_ras_to_lps[0] = v_ras_to_lps[1] = -1.0;
    vnl_diag_matrix<double> m_ras_to_lps(v_ras_to_lps);

    vnl_matrix<double> amatvnl = amat.GetVnlMatrix();
    amatvnl = m_ras_to_lps * amatvnl * m_ras_to_lps;
    vnl_vector_fixed<double, VDim > aoffs ;
    vnl_vector_fixed<double, VDim + 1> aoffl ;
    aoffs = m_ras_to_lps * aoff.GetVnlVector();
    aoffl.fill(1.0);
    for (size_t i=0; i<VDim; i++)
      aoffl(i) = aoffs(i);

    matrix.GetVnlMatrix().set_identity();
    matrix.GetVnlMatrix().update( amatvnl, 0, 0);
    matrix.GetVnlMatrix().set_column(VDim, aoffl);


}

template <class TPixel, unsigned int VDim>
void
ApplyMetric<TPixel, VDim>
::Flip_RAS_to_LPS(  itk::Matrix<double,VDim+1,VDim+1> &matrix,
                    itk::Matrix<double,VDim,VDim> &amat,
                    itk::Vector<double, VDim> &aoff)
{
  // Get the transform matrix and the offset vector
  vnl_matrix<double> A_ras = matrix.GetVnlMatrix().extract(VDim, VDim); 
  vnl_vector<double> b_ras = matrix.GetVnlMatrix().extract(VDim, 1, 0, VDim).get_column(0);

  // Extrernal matrices are assumed to be RAS to RAS, so we must convert to LPS to LPS
  vnl_vector<double> v_lps_to_ras(VDim, 1.0);
  v_lps_to_ras[0] = v_lps_to_ras[1] = -1.0;
  vnl_diag_matrix<double> m_lps_to_ras(v_lps_to_ras);
  vnl_matrix<double> A_lps = m_lps_to_ras * A_ras * m_lps_to_ras;
  vnl_vector<double> b_lps = m_lps_to_ras * b_ras;

  amat = A_lps;
  aoff.SetVnlVector(b_lps);
}



template <class TPixel, unsigned int VDim>
void
ApplyMetric<TPixel, VDim>
::operator() (const char *metric_name, const char *fn_ftran, const char *fn_mtran)
{
  //typedef itk::AffineTransform<double, VDim> TransformType;
  // Two images must be on a stack
  if(c->m_ImageStack.size() < 2)
    throw ConvertException("Two images required for metric computation");

  

  // Get the images
  ImagePointer iref = c->m_ImageStack[c->m_ImageStack.size() - 2];
  ImagePointer imov = c->m_ImageStack[c->m_ImageStack.size() - 1];

  *c->verbose << "Fixed  Image Transform: " << fn_ftran << endl;
  *c->verbose << "Moving Image Transform: " << fn_mtran << endl;

  // Create the appropriate metric
  typedef itk::ImageToImageMetric<ImageType, ImageType> MetricType;
  typename MetricType::Pointer metric;
  if(!strcmp(metric_name,"MI"))
    {
    metric = 
      itk::MutualInformationHistogramImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else if(!strcmp(metric_name,"MMI"))
    {
    metric = 
      itk::MattesMutualInformationImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else if(!strcmp(metric_name,"MSQ"))
    {
    metric = 
      itk::MeanSquaresImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else if(!strcmp(metric_name,"NCOR"))
    {
    metric = 
      itk::NormalizedCorrelationImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else if(!strcmp(metric_name,"NMI"))
    {
    metric = 
      itk::NormalizedMutualInformationHistogramImageToImageMetric<
        ImageType,ImageType>::New();
    }
  else throw ConvertException("Unknown metric %s", metric_name);

  // Configure the identity transform
  typename TransformType::Pointer tran = TransformType::New();


  double mvalue;

  if (!strcmp(fn_mtran,"none"))
    {
    tran->SetIdentity();
    metric->SetInterpolator(NNInterpolatorType::New());
    }
  else if (!strcmp(fn_ftran,"none"))
    {
    // Read the matrix
    itk::Matrix<double,VDim+1,VDim+1> matrix;
    itk::Matrix<double,VDim,VDim> amat;
    itk::Vector<double, VDim> aoff;

    ReadMatrix(fn_mtran, matrix);
    Flip_RAS_to_LPS(matrix, amat, aoff);

    // Put the values in the transform
    tran->SetMatrix(amat);
    tran->SetOffset(aoff);
    metric->SetInterpolator(LinInterpolatorType::New());
    }

  if (!strcmp(fn_ftran,"none"))
    {  
    //metric->DebugOn();

    // Configure the metric
    metric->SetMovingImage(imov);
    metric->SetFixedImage(iref);
    metric->SetTransform(tran);
    metric->SetFixedImageRegion(iref->GetBufferedRegion());
    metric->Initialize();

    //std::cout << "metric transform parameters: " << metric->GetTransform()->GetParameters() << endl;
    mvalue = metric->GetValue(tran->GetParameters());
    }
  else
    {

    // Create the image
    ImagePointer halfway = ImageType::New();
    CreateHalfwayImageSpace( iref, imov, halfway );

    typename TransformType::Pointer ftran = TransformType::New();
    typename TransformType::Pointer mtran = TransformType::New();

    // Read the matrix
    itk::Matrix<double,VDim+1,VDim+1> fmatrix;
    itk::Matrix<double,VDim,VDim> famat;
    itk::Vector<double, VDim> faoff;

    ReadMatrix(fn_ftran, fmatrix);
    Flip_RAS_to_LPS(fmatrix, famat, faoff);

    // Put the values in the transform
    ftran->SetMatrix(famat);
    ftran->SetOffset(faoff);
    
    // Read the matrix
    itk::Matrix<double,VDim+1,VDim+1> mmatrix;
    itk::Matrix<double,VDim,VDim> mamat;
    itk::Vector<double, VDim> maoff;

    ReadMatrix(fn_mtran, mmatrix);
    Flip_RAS_to_LPS(mmatrix, mamat, maoff);

    // Put the values in the transform
    mtran->SetMatrix(mamat);
    mtran->SetOffset(maoff);
    
    mvalue = GetValueInternalSymmetric( iref, imov, halfway, 
                                        ftran, mtran, metric_name );


    }

    // Print the value
    cout << metric_name << " = " << mvalue << endl;

    // Get image from stack
    // ImagePointer img = c->m_ImageStack.back();

    // Do some processing ...
    // ImagePointer result = ...;
  
    // Put result on stack
    // c->m_ImageStack.pop_back();
    // c->m_ImageStack.push_back(result);
}

template <class TPixel, unsigned int VDim>
void
ApplyMetric<TPixel, VDim>
::CreateHalfwayImageSpace( ImagePointer fixed, ImagePointer moving, ImagePointer halfway)
{
//  c3d_affine_tool -sform $TP1 -sform $TP0 -inv -mult -sqrt -sform $TP0 -mult -o $WDIR/hwspace.mat

  MatrixType mfixed  = fixed->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  MatrixType mmoving = moving->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  
  MatrixType mcomb = mmoving * vnl_matrix_inverse<double>(mfixed);
  // Peform Denman-Beavers iteration
  MatrixType Z, Y = mcomb;
  Z.set_identity();

  for(size_t i = 0; i < 16; i++) 
    {    
    MatrixType Ynext = 0.5 * (Y + vnl_matrix_inverse<double>(Z));
    MatrixType Znext = 0.5 * (Z + vnl_matrix_inverse<double>(Y));
    Y = Ynext;
    Z = Znext;
    }    

  MatrixType mhalf = Y * mfixed;


  halfway->SetBufferedRegion(fixed->GetBufferedRegion());
  // Set the voxel size
  halfway->SetSpacing(fixed->GetSpacing());
  halfway->Allocate();
  halfway->FillBuffer( 0.0 );

  // Set the matrix
  halfway->SetVoxelSpaceToRASPhysicalSpaceMatrix( mhalf );

}

template <class TPixel, unsigned int VDim>
double
ApplyMetric<TPixel, VDim>
::GetValueInternalSymmetric( ImagePointer fixed, ImagePointer moving, ImagePointer halfway, 
                                  TransformPointer ftran, 
                                  TransformPointer mtran, const char * metric_name )
{

  typedef  itk::ImageRegionConstIteratorWithIndex< ImageType > HalfwayIteratorType;
  typedef  itk::Point<double, VDim> InputPointType;
  typedef  itk::Point<double, VDim> OutputPointType;
  bool subtractmean = false;
  double measure = 0.0;
  int nPixels = 0;
  typename LinInterpolatorType::Pointer movinginterpolator = LinInterpolatorType::New();
  typename LinInterpolatorType::Pointer fixedinterpolator = LinInterpolatorType::New();
  movinginterpolator->SetInputImage( moving );
  fixedinterpolator->SetInputImage( fixed );

  if(!strcmp(metric_name,"MSQ"))
    {
    //std::cerr << "Internal method: " << parameters << std::endl;

    HalfwayIteratorType ti( halfway , halfway->GetBufferedRegion() );


    typename ImageType::IndexType index;



    while(!ti.IsAtEnd())
      {

      index = ti.GetIndex();

      InputPointType inputPoint;
      halfway->TransformIndexToPhysicalPoint( index, inputPoint );

      OutputPointType transformedFixedPoint = ftran->TransformPoint( inputPoint );

      if( !fixedinterpolator->IsInsideBuffer( transformedFixedPoint ))
        {
        ++ti;
        continue;
        }

      OutputPointType transformedMovingPoint = mtran->TransformPoint( inputPoint );

      if( !movinginterpolator->IsInsideBuffer( transformedMovingPoint ) )
        {
        ++ti;
        continue;
        }

      if( movinginterpolator->IsInsideBuffer( transformedMovingPoint ) &&
          fixedinterpolator->IsInsideBuffer( transformedFixedPoint ))
        {
        const double movingValue  = movinginterpolator->Evaluate( transformedMovingPoint );
        const double fixedValue   = fixedinterpolator->Evaluate( transformedFixedPoint );
        nPixels++;
        const double diff = movingValue - fixedValue;
        measure += diff * diff;
        }

      ++ti;
      }

    if( !nPixels )
      {
      throw ConvertException("All the points mapped to outside of the moving image");
      }
    else
      {
      measure /= nPixels;
      }
    }
  else if(!strcmp(metric_name,"NCOR"))
    {

    typedef  itk::ImageRegionConstIteratorWithIndex< ImageType > HalfwayIteratorType;


    HalfwayIteratorType ti( halfway , halfway->GetBufferedRegion() );

    typename ImageType::IndexType index;


    double sff = 0.0;
    double smm = 0.0;
    double sfm = 0.0;
    double sf  = 0.0;
    double sm  = 0.0;

    while(!ti.IsAtEnd())
      {

      index = ti.GetIndex();

      InputPointType inputPoint;
      halfway->TransformIndexToPhysicalPoint( index, inputPoint );
      OutputPointType transformedFixedPoint = ftran->TransformPoint( inputPoint );

      if( !fixedinterpolator->IsInsideBuffer( transformedFixedPoint ))
        {
        ++ti;
        continue;
        }

      OutputPointType transformedMovingPoint = mtran->TransformPoint( inputPoint );

      if( !movinginterpolator->IsInsideBuffer( transformedMovingPoint ) )
        {
        ++ti;
        continue;
        }

      if( movinginterpolator->IsInsideBuffer( transformedMovingPoint ) &&
          fixedinterpolator->IsInsideBuffer( transformedFixedPoint ))
        {
        const double movingValue  = movinginterpolator->Evaluate( transformedMovingPoint );
        const double fixedValue   = fixedinterpolator->Evaluate( transformedFixedPoint );

        sff += fixedValue  * fixedValue;
        smm += movingValue * movingValue;
        sfm += fixedValue  * movingValue;
        if ( subtractmean )
          {
          sf += fixedValue;
          sm += movingValue;
          }
        nPixels++;
        }

      ++ti;
      }

    if ( subtractmean && nPixels > 0 )
      {
      sff -= ( sf * sf / nPixels );
      smm -= ( sm * sm / nPixels );
      sfm -= ( sf * sm / nPixels );
      }

    const double denom = -1.0 * vcl_sqrt(sff * smm );

    if( nPixels > 0 && denom != 0.0)
      {
      measure = sfm / denom;
      }
    else
      {
      measure = 0.0;
      }

    }
    else
      throw ConvertException("Metric %s not supported for symmetric computation", metric_name);
      
  return measure;

}



// Invocations
template class ApplyMetric<double, 2>;
template class ApplyMetric<double, 3>;
template class ApplyMetric<double, 4>;
