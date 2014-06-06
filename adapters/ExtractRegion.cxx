/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ExtractRegion.cxx
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

#include "ExtractRegion.h"
#include "itkRegionOfInterestImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ExtractRegion<TPixel, VDim>
::operator() (RegionType bbox)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // Make sure the bounding box is within the contents of the image
  bbox.Crop(input->GetBufferedRegion());

  // Report the bounding box size
  *c->verbose << "  Extracting bounding box " << bbox.GetIndex() << " " << bbox.GetSize() << endl;

  // Chop off the region
  typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> TrimFilter;
  typename TrimFilter::Pointer fltTrim = TrimFilter::New();
  fltTrim->SetInput(input);
  fltTrim->SetRegionOfInterest(bbox);
  fltTrim->Update();

  // What happened to the origin of the image?
  ImagePointer output = fltTrim->GetOutput();

  // Update the image stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(output);
}

// Invocations
template class ExtractRegion<double, 2>;
template class ExtractRegion<double, 3>;
template class ExtractRegion<double, 4>;
