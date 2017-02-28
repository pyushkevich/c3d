/*=========================================================================

  Program:   ALFABIS fast medical image registration programs
  Language:  C++
  Website:   github.com/pyushkevich/greedy
  Copyright (c) Paul Yushkevich, University of Pennsylvania. All rights reserved.

  This program is part of ALFABIS: Adaptive Large-Scale Framework for
  Automatic Biomedical Image Segmentation.

  ALFABIS development is funded by the NIH grant R01 EB017255.

  ALFABIS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ALFABIS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ALFABIS.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/
#ifndef __UnaryFunctorVectorImageFilter_txx_
#define __UnaryFunctorVectorImageFilter_txx_

#include <itkImageLinearIteratorWithIndex.h>
#include "ImageRegionConstIteratorWithIndexOverride.h"
#include "itkProgressReporter.h"

/**
 * Constructor
 */
template< typename TInputImage, typename TOutputImage, typename TFunction  >
UnaryFunctorVectorImageFilter< TInputImage, TOutputImage, TFunction >
::UnaryFunctorVectorImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->InPlaceOff();
}

template< typename TInputImage, typename TOutputImage, typename TFunction  >
void
UnaryFunctorVectorImageFilter< TInputImage, TOutputImage, TFunction >
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();
  OutputImageType *outPtr = this->GetOutput();
  if(outPtr)
    outPtr->SetNumberOfComponentsPerPixel(this->GetFunctor().GetNumberOfComponentsPerPixel());
}

template< typename TInputImage, typename TOutputImage, typename TFunction  >
void
UnaryFunctorVectorImageFilter< TInputImage, TOutputImage, TFunction >
::ThreadedGenerateData(
    const OutputImageRegionType & outputRegionForThread,
    itk::ThreadIdType threadId)
{
  // Get the image
  InputImageType *input = const_cast<InputImageType *>(this->GetInput());
  OutputImageType *output = this->GetOutput();

  // Set up the iterator that will go through all the lines in the
  // output region. We assume that the lines span the whole length of
  // the input, i.e., the threading direction does not interfere
  typedef itk::ImageLinearIteratorWithIndex<OutputImageType> OutputIteratorBaseType;
  typedef IteratorExtenderWithOffset<OutputIteratorBaseType> OutputIteratorType;

  // This is the line iterator, although for even greater speed we operate
  // directly on pointers, so we only use it's NextLine functionality()
  OutputIteratorType itOutput(output, outputRegionForThread);

  // Create the input iterator that similarly iterates over line
  typedef itk::ImageLinearIteratorWithIndex<InputImageType> InputIteratorType;
  InputIteratorType itInput(input, outputRegionForThread);

  // Get the number of components
  int nc = output->GetNumberOfComponentsPerPixel();

  // Get the line length
  int line_length = outputRegionForThread.GetSize(0);

  // Progress
  const size_t numberOfLinesToProcess = outputRegionForThread.GetNumberOfPixels() / line_length;
  itk::ProgressReporter progress( this, threadId, numberOfLinesToProcess );

  // Start iterating over lines
  while(!itOutput.IsAtEnd())
    {
    // Pointer to the beginning of the scan line
    long offset_in_pixels = itOutput.GetPosition() - output->GetBufferPointer();
    long offset_in_comp = offset_in_pixels * nc;

    // Where we are scanning from
    OutputImageComponentType *p_scan_pixel = output->GetBufferPointer() + offset_in_comp;
    
    // Iterate over this line
    for(int p = 0; p < line_length; p++, p_scan_pixel+=nc, ++itInput)
      {
      m_Functor(itInput.Get(), p_scan_pixel);
      }

    itInput.NextLine();
    itOutput.NextLine();
    progress.CompletedPixel();
    }
}



#endif

