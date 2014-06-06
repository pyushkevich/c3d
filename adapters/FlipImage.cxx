/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    FlipImage.cxx
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

#include "FlipImage.h"
#include "itkFlipImageFilter.h"

template <class TPixel, unsigned int VDim>
void
FlipImage<TPixel, VDim>
::operator() (string axis)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Create a flip filter
  typedef itk::FlipImageFilter<ImageType> FlipType;
  typename FlipType::Pointer flipper = FlipType::New();
  typename FlipType::FlipAxesArrayType flipax;

  // Set up the axes to flip
  for(size_t i = 0; i < VDim; i++)
    if(axis.find('x'+i) != string::npos || axis.find('X'+i) != string::npos)
      flipax[i] = true;
    else
      flipax[i] = false;

  // Say what we are doing
  *c->verbose << "Flipping #" << c->m_ImageStack.size()-1 << " about " << flipax << endl;

  // Do some processing ...
  flipper->SetInput(img);
  flipper->SetFlipAxes(flipax);
  flipper->Update();

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(flipper->GetOutput());
}

// Invocations
template class FlipImage<double, 2>;
template class FlipImage<double, 3>;
template class FlipImage<double, 4>;
