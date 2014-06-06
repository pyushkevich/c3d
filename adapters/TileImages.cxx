/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    TileImages.cxx
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

#include "TileImages.h"

#include <itkTileImageFilter.h>

template <class TPixel, unsigned int VDim>
void
TileImages<TPixel, VDim>
::operator() (const std::string &tileParam)
{
  // Create the tile filter
  typedef typename itk::TileImageFilter<ImageType, ImageType> TileFilterType;
  typename TileFilterType::Pointer fltTile = TileFilterType::New();

  // Add all of the images as input to the tile filter
  for(int i = 0; i < c->m_ImageStack.size(); i++)
    {
    fltTile->SetInput(i, c->m_ImageStack[i]);
    }

  // Set the layout. The layout can either be a letter (x, y, z) which means tiling
  // in one dimension (i.e., tile a bunch of pngs) or it can be an integer vector
  // 2x2x0 which matches the input of the tile filter
  typename TileFilterType::LayoutArrayType loArray;
  
  if(tileParam == "x" || tileParam == "X" || tileParam == "0")
    {
    loArray.Fill(1);
    loArray[0] = c->m_ImageStack.size();
    }
  else if(tileParam == "y" || tileParam == "Y" || tileParam == "1")
    {
    loArray.Fill(1);
    loArray[1] = c->m_ImageStack.size();
    }
  else if(tileParam == "z" || tileParam == "Z" || tileParam == "2")
    {
    if(VDim < 3) throw ConvertException("Can not tile in z-dimension using c2d, use c3d");
    loArray.Fill(1);
    loArray[2] = c->m_ImageStack.size();
    }
  else if(tileParam == "w" || tileParam == "W" || tileParam == "t" || tileParam == "T" || tileParam == "3")
    {
    if(VDim < 4) throw ConvertException("Can not tile in w-dimension using c3d, use c4d");
    loArray.Fill(1);
    loArray[3] = c->m_ImageStack.size();
    }
  else
    {
    SizeType sz = c->ReadSizeVector(tileParam.c_str());
    for(int i = 0; i < VDim; i++)
      loArray[i] = sz[i];
    }

  fltTile->SetLayout(loArray);

  *c->verbose << "Tiling " << c->m_ImageStack.size() 
    << " images using layout " << loArray << endl; 

  fltTile->Update();
  
  // Put result on stack
  c->m_ImageStack.clear();
  c->m_ImageStack.push_back(fltTile->GetOutput());
}

// Invocations
template class TileImages<double, 2>;
template class TileImages<double, 3>;
template class TileImages<double, 4>;
