/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    CreateImage.cxx
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

#include "CreateImage.h"
#include "itkImage.h"

template <class TPixel, unsigned int VDim>
void
CreateImage<TPixel, VDim>
::operator()(SizeType dims, RealVector voxelSize)
{
  // Create a new region
  RegionType region;
  region.SetSize(dims);

  // Create the image
  ImagePointer img = ImageType::New();
  img->SetRegions(region);
  img->Allocate();
  img->FillBuffer(c->m_Background);
  
  // Set the voxel size
  img->SetSpacing(voxelSize.data_block());

  // Report
  *c->verbose << "Creating #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Dimensions: " << dims << endl;
  *c->verbose << "  Spacing: " << voxelSize << endl;
  
  // Push the image into the stack
  c->m_ImageStack.push_back(img);
}

// Invocations
template class CreateImage<double, 2>;
template class CreateImage<double, 3>;
template class CreateImage<double, 4>;
