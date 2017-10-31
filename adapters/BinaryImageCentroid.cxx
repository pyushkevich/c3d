/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    BinaryImageCentroid.cxx
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

#include "BinaryImageCentroid.h"

template <class TPixel, unsigned int VDim>
itk::ContinuousIndex<double, VDim> 
BinaryImageCentroid<TPixel, VDim>
::GetCentroid()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  *c->verbose << "Computing centroid of #" << c->m_ImageStack.size() << endl;

  // Find all pixels that do not match background
  itk::ContinuousIndex<double, VDim> center_idx;
  center_idx.Fill(0.0);
  size_t n = 0;
  for(Iterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(it.Get() != c->m_Background)
      {
      n++;
      for(size_t d = 0; d < VDim; d++)
        center_idx[d] += it.GetIndex()[d];
      }      
    }

  // Find the centroid in voxel units
  for(size_t d = 0; d < VDim; d++)
    center_idx[d] /= n;

  return center_idx;
}

template <class TPixel, unsigned int VDim>
void
BinaryImageCentroid<TPixel, VDim>
::operator() ()
{
  // Find the centroid in Nifti units
  itk::ContinuousIndex<double, VDim> center_idx = this->GetCentroid();
  itk::Point<double, VDim> center_pt;
  c->PeekLastImage()->TransformContinuousIndexToRASPhysicalPoint(center_idx, center_pt);

  // Print the resuts
  cout << "CENTROID_VOX " << center_idx << endl;
  cout << "CENTROID_MM " << center_pt << endl;
}

/**
 * This version of the method marks the centroid
 */
template <class TPixel, unsigned int VDim>
void
BinaryImageCentroid<TPixel, VDim>
::operator()(TPixel centroid_value)
{
  // Find the centroid
  itk::ContinuousIndex<double, VDim> centroid = this->GetCentroid();

  // Create a copy of the image on the stack - to safely work in-place
  ImageType *image = c->PopAndPushCopy();

  // Fill with zero
  image->FillBuffer(0);

  // Round the centroid to index location
  IndexType idx;
  for(int d = 0; d < VDim; d++)
    idx[d] = (long)(centroid[d] + 0.5);

  // Set the pixel
  image->SetPixel(idx, centroid_value);
  image->Modified();
}

// Invocations
template class BinaryImageCentroid<double, 2>;
template class BinaryImageCentroid<double, 3>;
template class BinaryImageCentroid<double, 4>;
