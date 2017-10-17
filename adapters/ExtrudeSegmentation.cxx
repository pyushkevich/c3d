/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ExtrudeSegmentation.cxx
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

#include "ExtrudeSegmentation.h"

#include "itkInPlaceImageFilter.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkProgressReporter.h"

template <class TInputImage, class TOutputImage, class TFunction>
class LineFunctorImageFilter 
  : public itk::InPlaceImageFilter<TInputImage, TOutputImage>
{
public:
  typedef LineFunctorImageFilter<TInputImage,TOutputImage,TFunction> Self;
  typedef itk::InPlaceImageFilter<TInputImage, TInputImage>    Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkTypeMacro(LineFunctorImageFilter, itk::InPlaceImageFilter)

  itkNewMacro(Self)

  /** Some typedefs. */
  typedef TFunction                                        FunctorType;

  typedef TInputImage                                      InputImageType;
  typedef typename    InputImageType::ConstPointer         InputImagePointer;
  typedef typename    InputImageType::RegionType           InputImageRegionType;
  typedef typename    InputImageType::PixelType            InputImagePixelType;

  typedef TOutputImage                                     OutputImageType;
  typedef typename     OutputImageType::Pointer            OutputImagePointer;
  typedef typename     OutputImageType::RegionType         OutputImageRegionType;
  typedef typename     OutputImageType::PixelType          OutputImagePixelType;

  /** Get the functor object.  The functor is returned by reference.
   * (Functors do not have to derive from itk::LightObject, so they do
   * not necessarily have a reference count. So we cannot return a
   * SmartPointer.) */
  FunctorType &       GetFunctor() { return m_Functor; }
  const FunctorType & GetFunctor() const { return m_Functor; }

  /** Set the functor object.  This replaces the current Functor with a
   * copy of the specified Functor. This allows the user to specify a
   * functor that has ivars set differently than the default functor.
   * This method requires an operator!=() be defined on the functor
   * (or the compiler's default implementation of operator!=() being
   * appropriate). */
  void SetFunctor(const FunctorType & functor)
  {
    if ( m_Functor != functor )
      {
      m_Functor = functor;
      this->Modified();
      }
  }

protected:
  LineFunctorImageFilter();
  ~LineFunctorImageFilter() {}

  void ThreadedGenerateData(
    const OutputImageRegionType & outputRegionForThread,
    itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
  LineFunctorImageFilter(const Self &);         // purposely not implemented
  void operator=(const Self &);                 // purposely not implemented

  FunctorType m_Functor;
};


/**
 * Constructor
 */
template< typename TInputImage, typename TOutputImage, typename TFunction  >
LineFunctorImageFilter< TInputImage, TOutputImage, TFunction >
::LineFunctorImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->InPlaceOff();
}

template< typename TInputImage, typename TOutputImage, typename TFunction  >
void
LineFunctorImageFilter< TInputImage, TOutputImage, TFunction >
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
  typedef itk::ImageLinearIteratorWithIndex<OutputImageType> OutputIteratorType;

  // This is the line iterator, although for even greater speed we operate
  // directly on pointers, so we only use it's NextLine functionality()
  OutputIteratorType itOutput(output, outputRegionForThread);

  // Create the input iterator that similarly iterates over line
  typedef itk::ImageLinearIteratorWithIndex<InputImageType> InputIteratorType;
  InputIteratorType itInput(input, outputRegionForThread);

  // Get the line length
  int line_length = outputRegionForThread.GetSize(0);

  // Progress
  const size_t numberOfLinesToProcess = outputRegionForThread.GetNumberOfPixels() / line_length;
  itk::ProgressReporter progress( this, threadId, numberOfLinesToProcess );

  // Create a functor for this thread
  FunctorType thread_functor = m_Functor;

  // Start iterating over lines
  while(!itOutput.IsAtEnd())
    {
    // Iterate over this line
    for(int p = 0; p < line_length; p++, ++itInput, ++itOutput)
      {
      itOutput.Set(m_Functor(itInput.Get(), p, line_length));
      }

    itInput.NextLine();
    itOutput.NextLine();
    progress.CompletedPixel();
    }
}



// --------

template <class TPixel, unsigned int VDim>
class ExtrudeSegmentationLineFunctor
{
public:
  ExtrudeSegmentationLineFunctor() { m_Background = 0; }

  void SetBackground(TPixel value) { m_Background = value; }

  TPixel operator() (const TPixel &value, int line_pos, int line_len)
    {
    if(line_pos == 0)
      m_Activated = false;

    if(value == m_Background)
      m_Activated = true;

    return m_Activated ? m_Background : value;
    }

  bool operator != (const ExtrudeSegmentationLineFunctor<TPixel, VDim> &other)
    {
    return m_Background != other.m_Background;
    }

protected:
  TPixel m_Background;
  bool m_Activated;
};


template <class TPixel, unsigned int VDim>
class MinimumIntensityProjectionFunctor
{
public:
  MinimumIntensityProjectionFunctor() : m_LineMinimum(0) {}

  TPixel operator() (const TPixel &value, int line_pos, int line_len)
    {
    if(line_pos == 0)
      {
      m_LineMinimum = value;
      return value;
      }

    if(value < m_LineMinimum)
      m_LineMinimum = value;

    return m_LineMinimum;
    }

  bool operator != (const ExtrudeSegmentationLineFunctor<TPixel, VDim> &other)
    {
    return m_LineMinimum != other.m_LineMinimum;
    }

protected:
  TPixel m_LineMinimum;
};


template <class TPixel, unsigned int VDim>
void
ExtrudeSegmentation<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer image = c->PopImage();

  // Create the filter
  typedef MinimumIntensityProjectionFunctor<TPixel, VDim> FunctorType;
  typedef LineFunctorImageFilter<ImageType, ImageType, FunctorType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);
  
  // Get the background
  /*
  FunctorType functor;
  functor.SetBackground(c->m_Background);
  filter->SetFunctor(functor);
  */

  // Run the filter
  filter->Update();

  // Do some processing ...
  c->PushImage(filter->GetOutput());
}

// Invocations
template class ExtrudeSegmentation<double, 2>;
template class ExtrudeSegmentation<double, 3>;
template class ExtrudeSegmentation<double, 4>;
