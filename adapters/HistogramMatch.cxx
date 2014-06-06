/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    HistogramMatch.cxx
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

#include "HistogramMatch.h"
#include "itkHistogramMatchingImageFilter.h"

template <class TPixel, unsigned int VDim>
void
HistogramMatch<TPixel, VDim>
::operator() (int nmatch)
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Reslice operation requires two images on the stack" << endl;
    throw -1;
    }

  // Get the reference and source images
  ImagePointer iref = c->m_ImageStack[c->m_ImageStack.size() - 2];
  ImagePointer isrc = c->m_ImageStack.back();

  // Create the filter
  typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> HistogramFilter;
  typename HistogramFilter::Pointer filter = HistogramFilter::New();

  filter->SetReferenceImage(iref);
  filter->SetSourceImage(isrc);
  filter->SetNumberOfMatchPoints(nmatch);
  filter->ThresholdAtMeanIntensityOff();

  *c->verbose << "Histogram matching #" << c->m_ImageStack.size() 
    << " to reference" << c->m_ImageStack.size() - 1 << endl;
  *c->verbose << "  Number of match points: " << filter->GetNumberOfMatchPoints() << endl;
  *c->verbose << "  Number of histogram levels: " << filter->GetNumberOfHistogramLevels() << endl;

  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class HistogramMatch<double, 2>;
template class HistogramMatch<double, 3>;
template class HistogramMatch<double, 4>;
