/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    TrimImage.cxx
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

#include "TrimImage.h"
#include "ExtractRegion.h"

template<unsigned int VDim>
void ExpandRegion(itk::ImageRegion<VDim> &region, const itk::Index<VDim> &idx)
{
  if(region.GetNumberOfPixels() == 0)
    {
    region.SetIndex(idx);
    for(size_t i = 0; i < VDim; i++)
      region.SetSize(i, 1);
    }
  else {
    for(size_t i = 0; i < VDim; i++)
      {
      if(region.GetIndex(i) > idx[i])
        {
        region.SetSize(i, region.GetSize(i) + (region.GetIndex(i) - idx[i]));
        region.SetIndex(i, idx[i]);
        }
      else if(region.GetIndex(i) + (long) region.GetSize(i) <= idx[i]) {
        region.SetSize(i, 1 + idx[i] - region.GetIndex(i));
        }
      }
  }
}

template <class TPixel, unsigned int VDim>
void
TrimImage<TPixel, VDim>
::operator() (const RealVector &vec, TrimMode mode)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Debugging info
  *c->verbose << "Trimming #" << c->m_ImageStack.size() << endl;
  

  // Initialize the bounding box
  RegionType bbox;

  // Find the extent of the non-background region of the image
  Iterator it(input, input->GetBufferedRegion());
  for( ; !it.IsAtEnd(); ++it)
    if(it.Value() != c->m_Background)
      ExpandRegion(bbox, it.GetIndex());

  // Pad the region by radius specified by user
  if(mode == SPECIFY_MARGIN)
    {
    *c->verbose << "  Wrapping non-background voxels with margin of " << vec << " mm." << endl;
    typename ImageType::SizeType radius;
    for(size_t i = 0; i < VDim; i++)
      radius[i] = (int) ceil(vec[i] / input->GetSpacing()[i]);
    bbox.PadByRadius(radius);
    }
  else if(mode == SPECIFY_FINALSIZE)
    {
    *c->verbose << "  Wrapping non-background voxels to create a region of size " << vec << " mm." << endl;
    // Compute the radius to pad to
    for(size_t i = 0; i < VDim; i++)
      { 
      int sznew = (int) (0.5 + vec[i] / input->GetSpacing()[i]);
      int ixnew = bbox.GetIndex()[i] + bbox.GetSize()[i]/2 - sznew/2; 
      bbox.SetIndex(i, ixnew);
      bbox.SetSize(i, sznew);
      }
    }

  // Use the extract region code for the rest
  ExtractRegion<TPixel, VDim> extract(c);
  extract(bbox);
}

// Invocations
template class TrimImage<double, 2>;
template class TrimImage<double, 3>;
template class TrimImage<double, 4>;
