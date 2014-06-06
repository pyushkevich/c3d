/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ExtractSlice.cxx
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

#include "ExtractSlice.h"
#include "ExtractRegion.h"
#include <string>
#include <iostream>
#include "itkExtractImageFilter.h"
#include "itkImageFileWriter.h"


template <class TPixel, unsigned int VDim>
void
ExtractSlice<TPixel, VDim>
::operator() (string axis, char* position)
{
  // Check input availability
  if(c->m_ImageStack.size() < 1)
    throw ConvertException("No images on stack");

  // Get the image
  ImagePointer image = c->m_ImageStack.back();
  SizeType size = image->GetBufferedRegion().GetSize();

  // Process the first parameter
  unsigned int slicedir;
  if (!axis.compare("x") || !axis.compare("0"))
    slicedir = 0;
  else if (!axis.compare("y") || !axis.compare("1"))
    slicedir = 1;
  else if (!axis.compare("z") || !axis.compare("2"))
    slicedir = 2;
  else if (!axis.compare("w") || !axis.compare("3") || !axis.compare("t"))
    slicedir = 3;
  else
    throw ConvertException("first parameter to -slice must be x,y,z or w");

  // Process the percent parameter
  char *pos = new char[strlen(position)];
  strcpy(pos, position);
  double percent_pos;
  int slicepos;
  std::string s( pos );
  size_t ipos = s.rfind("%");
  if(ipos == s.size() - 1)
    {
    const char *tok = strtok(pos, "%");
    percent_pos = atof( tok );
    slicepos = (int)(0.5 + (percent_pos / 100.0) * (size[slicedir] -1)); 
    }
  else
    slicepos = atoi( pos );

  // Say what we are doing
  *c->verbose << "Extracting slice " << pos << " along " << axis 
    << " axis in image #" << c->m_ImageStack.size()-1 << endl;

  // Use the extractor to extract the actual region
  RegionType rslice = image->GetBufferedRegion();
  rslice.SetSize(slicedir, 1);
  rslice.SetIndex(slicedir, slicepos);

  ExtractRegion<TPixel, VDim> extractor(c);
  extractor(rslice);

  // If slicing in the last dimension, we are done
  if(slicedir == VDim - 1)
    return;

  // Now, transpose the image to make the last dimension 1 
  // (this is like MATLAB's reshape command)
  std::vector<size_t> reorder;
  for(size_t i = 0; i < VDim; i++)
    if(i != slicedir)
      reorder.push_back(i);
  reorder.push_back(slicedir);  // 0 -> 1,2,0  1 -> 0,2,1 ...

  // Create a new image for the slice
  ImagePointer imnew = ImageType::New();
  ImagePointer imext = c->m_ImageStack.back();
  imnew->CopyInformation(image);

  // Create new region, origin, spacing for the image
  RegionType reg;
  typename ImageType::PointType org;
  typename ImageType::SpacingType spc;
  typename ImageType::DirectionType dir;

  for(size_t i = 0; i < VDim; i++)
    {
    int j = reorder[i];
    reg.SetSize(i, imext->GetBufferedRegion().GetSize(j));
    reg.SetIndex(i, imext->GetBufferedRegion().GetIndex(j));
    org[i] = imext->GetOrigin()[i]; // not a bug!
    spc[i] = imext->GetSpacing()[j];
    for(size_t k = 0; k < VDim; k++)
      dir(k, i) = imext->GetDirection()(k, j);
    }

  imnew->SetRegions(reg);
  imnew->SetOrigin(org);
  imnew->SetSpacing(spc);
  imnew->SetDirection(dir);

  // Copy all the data
  imnew->Allocate();
  Iterator itrg(imnew, reg);
  ConstIterator isrc(imext, imext->GetBufferedRegion());
  for(; !isrc.IsAtEnd(); ++isrc, ++itrg)
    itrg.Set(isrc.Get());

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(imnew);

  delete pos;
}

// Invocations
template class ExtractSlice<double, 3>;
template class ExtractSlice<double, 4>;
template class ExtractSlice<double, 2>;

