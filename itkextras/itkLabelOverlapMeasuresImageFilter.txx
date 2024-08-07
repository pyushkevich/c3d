/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLabelOverlapMeasuresImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2010/04/24 00:03:21 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkLabelOverlapMeasuresImageFilter_txx
#define _itkLabelOverlapMeasuresImageFilter_txx

#include "itkLabelOverlapMeasuresImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

namespace itk {

#if defined(__GNUC__) && (__GNUC__ <= 2) //NOTE: This class needs a mutex for gnu 2.95
/** Used for mutex locking */
#define LOCK_HASHMAP this->m_Mutex.lock()
#define UNLOCK_HASHMAP this->m_Mutex.unlock()
#else
#define LOCK_HASHMAP
#define UNLOCK_HASHMAP
#endif

template<class TLabelImage>
LabelOverlapMeasuresImageFilter<TLabelImage>
::LabelOverlapMeasuresImageFilter()
{
  // this filter requires two input images
  this->SetNumberOfRequiredInputs( 2 );
  this->DynamicMultiThreadingOff();
}

template<class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if( this->GetSourceImage() )
    {
    LabelImagePointer source = const_cast
      <LabelImageType *>( this->GetSourceImage() );
    source->SetRequestedRegionToLargestPossibleRegion();
    }
  if( this->GetTargetImage() )
    {
    LabelImagePointer target = const_cast
      <LabelImageType *>( this->GetTargetImage() );
    target->SetRequestedRegionToLargestPossibleRegion();
    }
}

template<class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::EnlargeOutputRequestedRegion( DataObject *data )
{
  Superclass::EnlargeOutputRequestedRegion( data );
  data->SetRequestedRegionToLargestPossibleRegion();
}


template<class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::AllocateOutputs()
{
  // Pass the source through as the output
  LabelImagePointer image =
    const_cast<TLabelImage *>( this->GetSourceImage() );
  this->SetNthOutput( 0, image );

  // Nothing that needs to be allocated for the remaining outputs
}

template<class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::BeforeThreadedGenerateData()
{
  int numberOfThreads = this->GetNumberOfWorkUnits();

  // Resize the thread temporaries
  this->m_LabelSetMeasuresPerThread.resize( numberOfThreads );

  // Initialize the temporaries
  for( int n = 0; n < numberOfThreads; n++ )
    {
    this->m_LabelSetMeasuresPerThread[n].clear();
    }

  // Initialize the final map
  this->m_LabelSetMeasures.clear();
}

template<class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::AfterThreadedGenerateData()
{
  // Run through the map for each thread and accumulate the set measures.
  for( int n = 0; n < this->GetNumberOfWorkUnits(); n++ )
    {
    // iterate over the map for this thread
    for( MapConstIterator threadIt = this->m_LabelSetMeasuresPerThread[n].begin();
      threadIt != this->m_LabelSetMeasuresPerThread[n].end();
      ++threadIt )
      {
      // does this label exist in the cumulative stucture yet?
      MapIterator mapIt = this->m_LabelSetMeasures.find( ( *threadIt ).first );
      if( mapIt == this->m_LabelSetMeasures.end() )
        {
        // create a new entry
        typedef typename MapType::value_type MapValueType;
        mapIt = this->m_LabelSetMeasures.insert( MapValueType(
          (*threadIt).first, LabelSetMeasures() ) ).first;
        }

      // accumulate the information from this thread
      (*mapIt).second.m_Source += (*threadIt).second.m_Source;
      (*mapIt).second.m_Target += (*threadIt).second.m_Target;
      (*mapIt).second.m_Union += (*threadIt).second.m_Union;
      (*mapIt).second.m_Intersection +=
        (*threadIt).second.m_Intersection;
      (*mapIt).second.m_SourceComplement +=
        (*threadIt).second.m_SourceComplement;
      (*mapIt).second.m_TargetComplement +=
        (*threadIt).second.m_TargetComplement;
      } // end of thread map iterator loop
    } // end of thread loop
}

template<class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::ThreadedGenerateData( const RegionType& outputRegionForThread,
  ThreadIdType threadId )
{
  ImageRegionConstIterator<LabelImageType> ItS( this->GetSourceImage(),
    outputRegionForThread );
  ImageRegionConstIterator<LabelImageType> ItT( this->GetTargetImage(),
    outputRegionForThread );

  // support progress methods/callbacks
  ProgressReporter progress( this, threadId,
    2*outputRegionForThread.GetNumberOfPixels() );

  for( ItS.GoToBegin(), ItT.GoToBegin(); !ItS.IsAtEnd(); ++ItS, ++ItT )
    {
    LabelType sourceLabel = ItS.Get();
    LabelType targetLabel = ItT.Get();

    // is the label already in this thread?
    MapIterator mapItS =
      this->m_LabelSetMeasuresPerThread[threadId].find( sourceLabel );
    MapIterator mapItT =
      this->m_LabelSetMeasuresPerThread[threadId].find( targetLabel );

    if( mapItS == this->m_LabelSetMeasuresPerThread[threadId].end() )
      {
      // create a new label set measures object
      typedef typename MapType::value_type MapValueType;
      LOCK_HASHMAP;
      mapItS = this->m_LabelSetMeasuresPerThread[threadId].insert(
        MapValueType( sourceLabel, LabelSetMeasures() ) ).first;
      UNLOCK_HASHMAP;
      }

    if( mapItT == this->m_LabelSetMeasuresPerThread[threadId].end() )
      {
      // create a new label set measures object
      typedef typename MapType::value_type MapValueType;
      LOCK_HASHMAP;
      mapItT = this->m_LabelSetMeasuresPerThread[threadId].insert(
        MapValueType( targetLabel, LabelSetMeasures() ) ).first;
      UNLOCK_HASHMAP;
      }

    (*mapItS).second.m_Source++;
    (*mapItT).second.m_Target++;

    if( sourceLabel == targetLabel )
      {
      (*mapItS).second.m_Intersection++;
      (*mapItS).second.m_Union++;
      }
    else
      {
      (*mapItS).second.m_Union++;
      (*mapItT).second.m_Union++;

      (*mapItS).second.m_SourceComplement++;
      (*mapItT).second.m_TargetComplement++;
      }

    progress.CompletedPixel();
    }
}

/**
 *  measures
 */
template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetTotalOverlap()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_Intersection );
    denominator += static_cast<RealType>( (*mapIt).second.m_Target );
    }
  return ( numerator / denominator );
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetTargetOverlap( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_Intersection ) /
    static_cast<RealType>( (*mapIt).second.m_Target );
  return value;
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetUnionOverlap()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_Intersection );
    denominator += static_cast<RealType>( (*mapIt).second.m_Union );
    }
  return ( numerator / denominator );
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetUnionOverlap( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_Intersection ) /
    static_cast<RealType>( (*mapIt).second.m_Union );
  return value;
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetMeanOverlap()
{
  RealType uo = this->GetUnionOverlap();
  return ( 2.0 * uo / ( 1.0 + uo ) );
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetMeanOverlap( LabelType label )
{
  RealType uo = this->GetUnionOverlap( label );
  return ( 2.0 * uo / ( 1.0 + uo ) );
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeSimilarity()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += ( ( static_cast<RealType>( (*mapIt).second.m_Source ) -
      static_cast<RealType>( (*mapIt).second.m_Target ) ) );
    denominator += ( ( static_cast<RealType>( (*mapIt).second.m_Source ) +
      static_cast<RealType>( (*mapIt).second.m_Target ) ) );
    }
  return ( 2.0 * numerator / denominator );
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetVolumeSimilarity( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value = 2.0 *
    ( static_cast<RealType>( (*mapIt).second.m_Source ) -
      static_cast<RealType>( (*mapIt).second.m_Target ) ) /
    ( static_cast<RealType>( (*mapIt).second.m_Source ) +
      static_cast<RealType>( (*mapIt).second.m_Target ) );
  return value;
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetFalseNegativeError()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_TargetComplement );
    denominator += static_cast<RealType>( (*mapIt).second.m_Target );
    }
  return ( numerator / denominator );
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetFalseNegativeError( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_TargetComplement ) /
    static_cast<RealType>( (*mapIt).second.m_Target );
  return value;
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetFalsePositiveError()
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType>::Zero )
      {
      continue;
      }
    numerator += static_cast<RealType>( (*mapIt).second.m_SourceComplement );
    denominator += static_cast<RealType>( (*mapIt).second.m_Source );
    }
  return ( numerator / denominator );
}

template<class TLabelImage>
typename LabelOverlapMeasuresImageFilter<TLabelImage>::RealType
LabelOverlapMeasuresImageFilter<TLabelImage>
::GetFalsePositiveError( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast<RealType>( (*mapIt).second.m_SourceComplement ) /
    static_cast<RealType>( (*mapIt).second.m_Source );
  return value;
}

template<class TLabelImage>
void
LabelOverlapMeasuresImageFilter<TLabelImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

}


}// end namespace itk
#endif
