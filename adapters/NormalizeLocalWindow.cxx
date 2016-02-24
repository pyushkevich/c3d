/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    NormalizeLocalWindow.cxx
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

#include "NormalizeLocalWindow.h"
#include "OneDimensionalInPlaceAccumulateFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkVectorImage.h"
#include "VectorImageTools.h"

template <class TInputPixel, class TOutputPixel>
class NormalizeLocalWindowImageToStatsFunctor
{
public:
  // The first parameter is the image, and the second parameter is the mask
  TOutputPixel operator() (const TInputPixel &x, const TInputPixel &m)
    {
    TOutputPixel c;
    c[0] = m;
    c[1] = x * m;
    c[2] = x * x * m;
    return c;
    }

  bool operator != (const NormalizeLocalWindowImageToStatsFunctor<TInputPixel, TOutputPixel> &other)
    { return false; }
};

template <class TInputPixel1, class TInputPixel2, class TOutputPixel>
class NormalizedLocalWindowStatsToResultFunctor
{
public:
  typedef NormalizedLocalWindowStatsToResultFunctor<TInputPixel1, TInputPixel2, TOutputPixel> Self;

  TOutputPixel operator() (const TInputPixel1 &c, const TInputPixel2 &x)
    {
    // Get the raw data
    double sum_m = c[0];
    double sum_mx = c[1];
    double sum_mxx = c[2];

    // If no mask, return 0
    if(sum_m == 0)
      return 0;

    // Compute the standard deviation
    double var_xx = (sum_mxx - sum_mx * sum_mx / sum_m) / sum_m;

    // Compute the difference
    return (x - sum_mx / sum_m) / sqrt(var_xx);
    }

  bool operator != (const Self &other)
    { return false; }
};


template <class TPixel, unsigned int VDim>
void
NormalizeLocalWindow<TPixel, VDim>
::operator() (SizeType radius)
{
  // This filter replaces the intensity g(x) at each voxel by (g-mu)/sigma, where mu and sigma
  // are the mean and standard deviation of intensity in a neighborhood. To compute the filter,
  // we need to compute these sliding window statistics
  //
  // The first input is the image and the second input is a mask to which the computation of the
  // statistics is restricted. The statistics are only computed inside of the mask
  // Get two images from stack
  ImagePointer image = c->m_ImageStack[c->m_ImageStack.size()-2];
  ImagePointer mask = c->m_ImageStack[c->m_ImageStack.size()-1];

  // Alternative approach
  typedef itk::Vector<TPixel, 3> StatsVector;
  typedef itk::Image<StatsVector, VDim> StatsImage;
  typedef NormalizeLocalWindowImageToStatsFunctor<TPixel,StatsVector> StatsFunctor;
  typedef itk::BinaryFunctorImageFilter<ImageType, ImageType, StatsImage, StatsFunctor> StatsFilter;
  typename StatsFilter::Pointer fltStats = StatsFilter::New();
  fltStats->SetInput1(image);
  fltStats->SetInput2(mask);
  fltStats->Update();

  // Masquerade as a vector image
  typedef itk::VectorImage<TPixel, VDim> VecImageType;
  typename VecImageType::Pointer vecImage = WrapImageOfVectorAsVectorImage(fltStats->GetOutput());

  // Perform the accumulation
  typename VecImageType::Pointer imgAccum = AccumulateNeighborhoodSumsInPlace(vecImage.GetPointer(), radius);
    
  typedef NormalizedLocalWindowStatsToResultFunctor<typename VecImageType::PixelType, TPixel, TPixel> NormalizeFunctor;
  typedef itk::BinaryFunctorImageFilter<VecImageType, ImageType, ImageType, NormalizeFunctor> NormalizeFilter;
  typename NormalizeFilter::Pointer fltNorm = NormalizeFilter::New();

  fltNorm->SetInput1(imgAccum);
  fltNorm->SetInput2(image);
  fltNorm->Update();

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltNorm->GetOutput());
}

// Invocations
template class NormalizeLocalWindow<double, 2>;
template class NormalizeLocalWindow<double, 3>;
template class NormalizeLocalWindow<double, 4>;
