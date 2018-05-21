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
#ifndef __UnaryFunctorVectorImageFilter_h_
#define __UnaryFunctorVectorImageFilter_h_

#include "itkInPlaceImageFilter.h"

/**
 * This is a filter similar to UnaryFunctorImageFilter but allows output to
 * be a VectorImage. The main difference is that the functor must implement
 * GetNumberOfComponentsPerPixel() method. Additionally, the functor's 
 * operator has a different signature to avoid unnecessary memory allocation.
 *
 *  
 */
template <class TInputImage, class TOutputImage, class TFunction>
class UnaryFunctorVectorImageFilter : 
  public itk::InPlaceImageFilter<TInputImage, TOutputImage>
{
public:
  typedef UnaryFunctorVectorImageFilter<TInputImage,TOutputImage,TFunction> Self;
  typedef itk::InPlaceImageFilter<TInputImage, TOutputImage> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkTypeMacro(UnaryFunctorVectorImageFilter, itk::InPlaceImageFilter)

  itkNewMacro(Self)

  /** Some typedefs. */
  typedef TFunction FunctorType;

  typedef TInputImage                                      InputImageType;
  typedef typename    InputImageType::ConstPointer         InputImagePointer;
  typedef typename    InputImageType::RegionType           InputImageRegionType;
  typedef typename    InputImageType::PixelType            InputImagePixelType;

  typedef TOutputImage                                     OutputImageType;
  typedef typename     OutputImageType::Pointer            OutputImagePointer;
  typedef typename     OutputImageType::RegionType         OutputImageRegionType;
  typedef typename     OutputImageType::PixelType          OutputImagePixelType;
  typedef typename     OutputImageType::InternalPixelType  OutputImageComponentType;

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

  UnaryFunctorVectorImageFilter();
  ~UnaryFunctorVectorImageFilter() {}

  void GenerateOutputInformation() ITK_OVERRIDE;
  
  void ThreadedGenerateData(
    const OutputImageRegionType & outputRegionForThread,
    itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
  UnaryFunctorVectorImageFilter(const Self &);  //purposely not implemented
  void operator=(const Self &);                 //purposely not implemented

  FunctorType m_Functor;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "UnaryFunctorVectorImageFilter.txx"
#endif

#endif

