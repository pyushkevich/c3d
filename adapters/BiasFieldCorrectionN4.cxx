/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    BiasFieldCorrectionN4.cxx
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

#include "BiasFieldCorrectionN4.h"

#include "itkBSplineControlPointImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShrinkImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BiasFieldCorrectionN4<TPixel, VDim>
::operator() ()
{

  // Get image from stack
  ImagePointer mri = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();

  // Bias filter
  typedef itk::N4BiasFieldCorrectionImageFilter<ImageType, ImageType, ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();

  // distance (in mm) of the mesh resolution at the base level
  float splineDistance = 100;

  typename CorrecterType::ArrayType numberOfControlPoints;

  typename ImageType::IndexType inputImageIndex =
    mri->GetLargestPossibleRegion().GetIndex();
  typename ImageType::SizeType inputImageSize =
    mri->GetLargestPossibleRegion().GetSize();

  typename ImageType::PointType newOrigin = mri->GetOrigin();

  itk::SizeValueType lowerBound[VDim];
  itk::SizeValueType upperBound[VDim];

  for( unsigned int d = 0; d < VDim; d++ )
    {
    float domain = static_cast<float>( mri->
      GetLargestPossibleRegion().GetSize()[d] - 1 ) * mri->GetSpacing()[d];
    unsigned int numberOfSpans = static_cast<unsigned int>(
      vcl_ceil( domain / splineDistance ) );
    unsigned long extraPadding = static_cast<unsigned long>( ( numberOfSpans *
      splineDistance - domain ) / mri->GetSpacing()[d] + 0.5 );
    lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
    upperBound[d] = extraPadding - lowerBound[d];
    newOrigin[d] -= ( static_cast<float>( lowerBound[d] ) *
      mri->GetSpacing()[d] );

    numberOfControlPoints[d] = numberOfSpans + correcter->GetSplineOrder();
    }
  correcter->SetNumberOfControlPoints( numberOfControlPoints );

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
  typename PadderType::Pointer padder = PadderType::New();
  padder->SetInput( mri );
  padder->SetPadLowerBound( lowerBound );
  padder->SetPadUpperBound( upperBound );
  padder->SetConstant( 0 );
  padder->Update();

  typename PadderType::Pointer padder2 = PadderType::New();
  padder2->SetInput( padder->GetOutput() );
  padder2->SetPadLowerBound( lowerBound );
  padder2->SetPadUpperBound( upperBound );
  padder2->SetConstant( 0 );
  padder2->Update();

  // Set up a filter to shrink image by a factor
  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput( padder2->GetOutput() );
  shrinker->SetShrinkFactors( 4 );
  shrinker->Update();

  // Compute mask using Otsu threshold
  typedef itk::OtsuThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer otsu = ThresholderType::New();
  otsu->SetInput( padder->GetOutput() );
  otsu->SetNumberOfHistogramBins( 200 );
  otsu->SetInsideValue( 0 );
  otsu->SetOutsideValue( 1 );
  otsu->Update();
  ImagePointer mask = otsu->GetOutput();

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> MaskPadderType;
  typename MaskPadderType::Pointer maskPadder = MaskPadderType::New();
  maskPadder->SetInput( otsu->GetOutput() );
  maskPadder->SetPadLowerBound( lowerBound );
  maskPadder->SetPadUpperBound( upperBound );
  maskPadder->SetConstant( 0 );
  maskPadder->Update();

  // Shrink the mask
  typename ShrinkerType::Pointer maskshrinker = ShrinkerType::New();
  maskshrinker->SetInput( maskPadder->GetOutput() );
  maskshrinker->SetShrinkFactors( 4 );
  maskshrinker->Update();

  correcter->SetInput( shrinker->GetOutput() );
  correcter->SetMaskImage( maskshrinker->GetOutput() );

  // These parameters are pretty standard
  correcter->SetSplineOrder( 3 );
  correcter->SetNumberOfHistogramBins( 200 );
  correcter->SetBiasFieldFullWidthAtHalfMaximum( 0.15 );
  correcter->SetConvergenceThreshold( 0.001 );
  correcter->SetWienerFilterNoise( 0.01 );
  correcter->SetBiasFieldFullWidthAtHalfMaximum( 0.15 );

  // You will probably want to have an option for the maximum number of
  //  iterations at each level, the shrink factor, and the spline distance.
  typename CorrecterType::ArrayType numberOfFittingLevels;
  numberOfFittingLevels.Fill( 3 );
		correcter->SetNumberOfFittingLevels( numberOfFittingLevels );
  typename CorrecterType::VariableSizeArrayType maximumNumberOfIterations;
  maximumNumberOfIterations.SetSize( 3 );
  maximumNumberOfIterations[0] = 100;
  maximumNumberOfIterations[1] = 50;
  maximumNumberOfIterations[2] = 50;
  correcter->SetMaximumNumberOfIterations( maximumNumberOfIterations );

  // Progress meter
  // typedef CommandIterationUpdate<CorrecterType> CommandType;
  // typename CommandType::Pointer observer = CommandType::New();
  // correcter->AddObserver( itk::IterationEvent(), observer );
  correcter->Update();

  /**
   * Reconstruct the bias field at full image resolution.  Divide
   * the original input image by the bias field to get the final
   * corrected image.
   */
  typedef itk::BSplineControlPointImageFilter<typename
    CorrecterType::BiasFieldControlPointLatticeType, typename
    CorrecterType::ScalarImageType> BSplinerType;

  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
  bspliner->SetSplineOrder( correcter->GetSplineOrder() );
  bspliner->SetSize( mri->GetLargestPossibleRegion().GetSize() );
  bspliner->SetOrigin( mri->GetOrigin() );
  bspliner->SetDirection( mri->GetDirection() );
  bspliner->SetSpacing( mri->GetSpacing() );
  bspliner->Update();

  typename ImageType::Pointer logField = ImageType::New();
  logField->SetOrigin( bspliner->GetOutput()->GetOrigin() );
  logField->SetSpacing( bspliner->GetOutput()->GetSpacing() );
  logField->SetRegions( bspliner->GetOutput()->GetLargestPossibleRegion().GetSize() );
  logField->SetDirection( bspliner->GetOutput()->GetDirection() );
  logField->Allocate();

  itk::ImageRegionIterator<typename CorrecterType::ScalarImageType> ItB(
    bspliner->GetOutput(),
    bspliner->GetOutput()->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> ItF( logField,
    logField->GetLargestPossibleRegion() );
  for( ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF )
    {
    ItF.Set( ItB.Get()[0] );
    }

  typedef itk::ExpImageFilter<ImageType, ImageType> ExpFilterType;
  typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
  expFilter->SetInput( logField );
  expFilter->Update();

  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DividerType;
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput1( mri );
  divider->SetInput2( expFilter->GetOutput() );
  divider->Update();

  typename ImageType::RegionType inputRegion;
  inputRegion.SetIndex( inputImageIndex );
  inputRegion.SetSize( inputImageSize );

  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  typename CropperType::Pointer cropper = CropperType::New();
  cropper->SetInput( divider->GetOutput() );
  cropper->SetExtractionRegion( inputRegion );
  cropper->Update();

  typename CropperType::Pointer biasFieldCropper = CropperType::New();
  biasFieldCropper->SetInput( expFilter->GetOutput() );
  biasFieldCropper->SetExtractionRegion( inputRegion );
  biasFieldCropper->Update();

  // Update
  c->m_ImageStack.push_back( cropper->GetOutput() );
//   c->m_ImageStack.push_back( biasFieldCropper->GetOutput() );
}

// Invocations
template class BiasFieldCorrectionN4<double, 2>;
template class BiasFieldCorrectionN4<double, 3>;
template class BiasFieldCorrectionN4<double, 4>;
