/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    FastMarching.cxx
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

#include "FastMarching.h"
#include "itkFastMarchingImageFilter.h"

template <class TPixel, unsigned int VDim>
void
FastMarching<TPixel, VDim>
::operator() (double stopping_value)
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    throw ConvertException("Level set segmentation requires two images on the stack!");

  // Get the last two images
  ImagePointer init = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer speed = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Report what the filter is doing
  *c->verbose << "Running fast marching filter (";
  *c->verbose << "#" << c->m_ImageStack.size() - 1 << " is speed, ";
  *c->verbose << "#" << c->m_ImageStack.size() << " is init)" << endl;

  // Create filter
  using FilterType = itk::FastMarchingImageFilter<ImageType, ImageType>;
  typename FilterType::Pointer filter = FilterType::New();

  // All the non-zero pixels in the initial image are provided as trial points
  using NodeContainer = typename FilterType::NodeContainer;
  using NodeType = typename FilterType::NodeType;
  auto seeds = NodeContainer::New();
  seeds->Initialize();
  for(Iterator it(init, init->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(it.Get() > 0)
      {
      NodeType node;
      node.SetValue(1.0);
      node.SetIndex(it.GetIndex());
      seeds->push_back(node);
      }
    }

  // Set up the filter
  filter->SetTrialPoints(seeds);
  filter->SetInput(speed);
  filter->SetStoppingValue(stopping_value);
  filter->Update();

  // Do some processing ...
  ImagePointer result = filter->GetOutput();
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class FastMarching<double, 2>;
template class FastMarching<double, 3>;
template class FastMarching<double, 4>;
