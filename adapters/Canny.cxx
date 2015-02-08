/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    Canny.cxx
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

#include "Canny.h"

#include "itkCannyEdgeDetectionImageFilter.h"

template <class TPixel, unsigned int VDim>
void
Canny<TPixel, VDim>
::operator() (const RealVector &sigma, double tLower, double tUpper)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Create the filter
  typedef itk::CannyEdgeDetectionImageFilter<ImageType, ImageType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(img);

  // Set up the variance
  typename FilterType::ArrayType var;
  for(int d = 0; d < VDim; d++)
      var[d] = sigma[d] * sigma[d];
  filter->SetVariance(var);

  // Set the thresholds
  filter->SetLowerThreshold(tLower);
  filter->SetUpperThreshold(tUpper);

  *c->verbose << "Performing Canny edge detection on #" << c->m_ImageStack.size() << std::endl;
  *c->verbose << "  Variance        : " << var    << std::endl;
  *c->verbose << "  Lower Threshold : " << tLower << std::endl;
  *c->verbose << "  Upper Threshold : " << tUpper << std::endl;

  // Do some processing ...
  filter->Update();

  // Iterate through the output points and save them
  FILE *fcanny = fopen("canny.obj","wt");
  for(Iterator it(filter->GetOutput(), filter->GetOutput()->GetBufferedRegion());
    !it.IsAtEnd(); ++it)
    {
    if(it.Get())
      {
      itk::Point<double, VDim> p;
      filter->GetOutput()->TransformIndexToRASPhysicalPoint(it.GetIndex(), p);
      fprintf(fcanny, "v ");
      for(int i = 0; i < VDim; i++) fprintf(fcanny, "%f ", p[i]);
      for(int i = 0; i < VDim; i++) fprintf(fcanny, "255 ");
      fprintf(fcanny, "\n");
      }
    }
  fclose(fcanny);
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class Canny<double, 2>;
template class Canny<double, 3>;
template class Canny<double, 4>;
