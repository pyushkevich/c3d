/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    NormalizedCrossCorrelation.cxx
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

#include "NormalizedCrossCorrelation.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkUnaryFunctorImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMeanImageFilter.h"

#include "itkImageRegionSplitterDirection.h"
#include "itkTimeProbe.h"

#include "itkImageFileWriter.h"
#include "OneDimensionalInPlaceAccumulateFilter.h"
#include "VectorImageTools.h"


template <class TInputPixel, class TOutputPixel>
class ImageToStatsFunctor
{
public:
  TOutputPixel operator() (const TInputPixel &a, const TInputPixel &b)
    {
    TOutputPixel c;
    c[0] = a;
    c[1] = b;
    c[2] = a*b;
    c[3] = a*a;
    c[4] = b*b;
    c[5] = 1.0;
    return c;
    }

  bool operator != (const ImageToStatsFunctor<TInputPixel, TOutputPixel> &other)
    { return false; }
};

template <class TInputPixel, class TOutputPixel>
class StatsToNCCFunctor
{
public:
  StatsToNCCFunctor(int n = 0) { this->m_N = n; }

  TOutputPixel operator() (const TInputPixel &c)
    {
    double a = c[0], b = c[1], ab = c[2], a2 = c[3], b2 = c[4], n = c[5];
    double Sab = (ab - a * b / n); 
    double Saa = (a2 - a * a / n); 
    double Sbb = (b2 - b * b / n); 
    return Sab / sqrt(Saa * Sbb);
    }

  bool operator != (const StatsToNCCFunctor<TInputPixel, TOutputPixel> &other)
    { return m_N != other.m_N; }

private:
  int m_N;
};

template <class TPixel, unsigned int VDim>
void
NormalizedCrossCorrelation<TPixel, VDim>
::operator() (itk::Size<VDim> radius)
{
  // Get two images from stack
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size()-1];
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size()-2];

  // Alternative approach
  typedef itk::Vector<TPixel, 6> StatsVector;
  typedef itk::Image<StatsVector, VDim> StatsImage;
  typedef ImageToStatsFunctor<TPixel,StatsVector> StatsFunctor;
  typedef itk::BinaryFunctorImageFilter<ImageType, ImageType, StatsImage, StatsFunctor> StatsFilter;
  typename StatsFilter::Pointer fltStats = StatsFilter::New();
  fltStats->SetInput1(i1);
  fltStats->SetInput2(i2);
  fltStats->Update();

  // Masquerade as a vector image
  typedef itk::VectorImage<TPixel, VDim> VecImageType;
  typename VecImageType::Pointer vecImage = WrapImageOfVectorAsVectorImage(fltStats->GetOutput());

  // Perform the accumulation
  typename VecImageType::Pointer imgAccum = AccumulateNeighborhoodSumsInPlace(vecImage.GetPointer(), radius);
    
  typedef StatsToNCCFunctor<typename VecImageType::PixelType, TPixel> NCCFunctor;
  typedef itk::UnaryFunctorImageFilter<VecImageType, ImageType, NCCFunctor> NCCFilter;
  typename NCCFilter::Pointer fltNCC = NCCFilter::New();

  int n = 1;
  for(int i = 0; i < VDim; i++)
   n *= 2 * radius[i] + 1;

  NCCFunctor funk(n);
  fltNCC->SetInput(imgAccum);
  fltNCC->SetFunctor(funk);
  fltNCC->Update();

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltNCC->GetOutput());
}

// Invocations
template class NormalizedCrossCorrelation<double, 2>;
template class NormalizedCrossCorrelation<double, 3>;
template class NormalizedCrossCorrelation<double, 4>;
