#include "LabelOverlapMeasures.h"
#include "itkLabelOverlapMeasuresImageFilter.h"

#include <iomanip>

template <class TPixel, unsigned int VDim>
void
LabelOverlapMeasures<TPixel, VDim>
::operator() ()
{
  // Get images from stack
  size_t n = c->m_ImageStack.size();
  if( n < 2 )
    throw ConvertException( "Label overlap measures require two image inputs" );
  ImagePointer target = c->m_ImageStack[n-1];
  ImagePointer source = c->m_ImageStack[n-2];

  // Create a short image for the labels
  typedef itk::Image<short, VDim> LabelImageType;
  typename LabelImageType::Pointer slabSource = LabelImageType::New();
  typename LabelImageType::Pointer slabTarget = LabelImageType::New();

  // Allocate the images
  slabSource->SetRegions( source->GetBufferedRegion() );
  slabSource->Allocate();

  slabTarget->SetRegions( target->GetBufferedRegion() );
  slabTarget->Allocate();

  // Round off doubles to create labels
  size_t nvS = slabSource->GetBufferedRegion().GetNumberOfPixels();
  for( size_t i = 0; i < nvS; i++ )
    {
    slabSource->GetBufferPointer()[i]
      = (short) ( source->GetBufferPointer()[i] + 0.5 );
    }

  size_t nvT = slabTarget->GetBufferedRegion().GetNumberOfPixels();
  for( size_t i = 0; i < nvT; i++ )
    {
    slabTarget->GetBufferPointer()[i]
      = (short) ( target->GetBufferPointer()[i] + 0.5 );
    }

  // Create the label statistics fltOverlap
  typedef itk::LabelOverlapMeasuresImageFilter<LabelImageType> OverlapFilter;
  typename OverlapFilter::Pointer filter = OverlapFilter::New();

  // Set the inputs
  filter->SetSourceImage( slabSource );
  filter->SetTargetImage( slabTarget );

  // Update the fltOverlap
  filter->Update();

  std::cout << "                                          "
            << "************ All Labels *************" << std::endl;
  std::cout << std::setw( 10 ) << "   "
    << std::setw( 17 ) << "Total"
    << std::setw( 17 ) << "Union (jaccard)"
    << std::setw( 17 ) << "Mean (dice)"
    << std::setw( 17 ) << "Volume sim."
    << std::setw( 17 ) << "False negative"
    << std::setw( 17 ) << "False positive" << std::endl;
  std::cout << std::setw( 10 ) << "   ";
  std::cout << std::setw( 17 ) << filter->GetTotalOverlap();
  std::cout << std::setw( 17 ) << filter->GetUnionOverlap();
  std::cout << std::setw( 17 ) << filter->GetMeanOverlap();
  std::cout << std::setw( 17 ) << filter->GetVolumeSimilarity();
  std::cout << std::setw( 17 ) << filter->GetFalseNegativeError();
  std::cout << std::setw( 17 ) << filter->GetFalsePositiveError();
  std::cout << std::endl;

  std::cout << "                                       "
            << "************ Individual Labels *************" << std::endl;
  std::cout << std::setw( 10 ) << "Label"
            << std::setw( 17 ) << "Target"
            << std::setw( 17 ) << "Union (jaccard)"
            << std::setw( 17 ) << "Mean (dice)"
            << std::setw( 17 ) << "Volume sim."
            << std::setw( 17 ) << "False negative"
            << std::setw( 17 ) << "False positive" << std::endl;

  typename OverlapFilter::MapType labelMap = filter->GetLabelSetMeasures();
  typename OverlapFilter::MapType::const_iterator it;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
    {
    if( (*it).first == 0 )
      {
      continue;
      }

    int label = (*it).first;

    std::cout << std::setw( 10 ) << label;
    std::cout << std::setw( 17 ) << filter->GetTargetOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetUnionOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetMeanOverlap( label );
    std::cout << std::setw( 17 ) << filter->GetVolumeSimilarity( label );
    std::cout << std::setw( 17 ) << filter->GetFalseNegativeError( label );
    std::cout << std::setw( 17 ) << filter->GetFalsePositiveError( label );
    std::cout << std::endl;
    }

}

// Invocations
template class LabelOverlapMeasures<double, 2>;
template class LabelOverlapMeasures<double, 3>;
