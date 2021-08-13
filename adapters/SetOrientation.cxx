/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SetOrientation.cxx
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

#include "SetOrientation.h"

template <class TPixel, unsigned int VDim>
void
SetOrientation<TPixel, VDim>
::operator() (std::string rai)
{
  // Only valid for 4D images or less
  if(VDim > 4)
    throw ConvertException("Orientation codes only valid for up to 4D images");
  
  // Check the RAI code validity (length must match image dimension)
  // Exception for 4D images (RAI code must be 3 characters)
  if(rai.length() != VDim)
    {
    if(VDim != 4)
      throw ConvertException("Orientation code %s is not %d characters long", rai.c_str(), VDim);
    }
  if(VDim == 4 && rai.length() != 3)
    throw ConvertException("Orientation code %s is not 3 characters long", rai.c_str());


  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Create a direction matrix
  vnl_matrix_fixed<double, VDim, VDim> eye, dm;
  eye.set_identity(); dm.set_identity();

  // RAI codes
  char codes[3][2] = { {'R', 'L'}, {'A', 'P'}, {'I', 'S'}};
  
  for(size_t i = 0; i < VDim; i++)
    {
    bool matched = false;
    for(size_t j = 0; j < VDim; j++)
      {
      // 4D Cosine Matrix is identity for 4th col/row
      if(i==3)
        {
        // Exit loop
        matched = true;
        break;
        }
      for(size_t k = 0; k < 2; k++)
        {
        if(toupper(rai[i]) == codes[j][k])
          {
          // Set the row of the direction matrix
          dm.set_column(i, (k==0 ? 1.0 : -1.0) * eye.get_row(j));

          // Clear that code (so that we catch orientation codes like RRS)
          codes[j][0] = codes[j][1] = 'X';

          // We found a code for i
          matched = true;
          }
        }
      }

    if(!matched)
      throw ConvertException("Orientation code %s is invalid", rai.c_str());
    }

  // Explain what's being done
  *c->verbose << "Setting orientation of " << c->m_ImageStack.size() << " to " << rai << endl;

  // Set the direction in the image
  img->SetDirection(itk::Matrix<double,VDim,VDim>(dm));
}

// Invocations
template class SetOrientation<double, 2>;
template class SetOrientation<double, 3>;
template class SetOrientation<double, 4>;
