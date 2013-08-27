#include "BiasFieldCorrection.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BiasFieldCorrection<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer mri = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();

  // Set up a filter to shrink image by a factor
  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput(mri);
  shrinker->SetShrinkFactors(4);
  shrinker->Update();

  // Compute mask using Otsu threshold
  typedef itk::OtsuThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer otsu = ThresholderType::New();
  otsu->SetInput(mri);
  otsu->SetNumberOfHistogramBins( 200 );
  otsu->SetInsideValue( 0 );
  otsu->SetOutsideValue( 1 );
  otsu->Update();
  ImagePointer mask = otsu->GetOutput();

  // Shrink the mask
  typename ShrinkerType::Pointer maskshrinker = ShrinkerType::New();
  maskshrinker->SetInput( mask);
  maskshrinker->SetShrinkFactors(4);
  maskshrinker->Update();

  // Bias filter
  typedef itk::N3MRIBiasFieldCorrectionImageFilter<ImageType, ImageType, ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();
  correcter->SetInput( shrinker->GetOutput() );
  correcter->SetMaskImage( maskshrinker->GetOutput() );

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

  // Update
  c->m_ImageStack.push_back(divider->GetOutput());
}

// Invocations
template class BiasFieldCorrection<double, 2>;
template class BiasFieldCorrection<double, 3>;
