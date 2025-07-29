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
#include <sstream>
#include <iostream>
#include "itkExtractImageFilter.h"
#include "itkImageFileWriter.h"

template <class TPixel, unsigned int VDim>
void
ExtractSlice<TPixel, VDim>
::ExtractOneSlice(ImageType *image, unsigned int slicedir, int slicepos)
{
  // Say what we are doing
  static const char axis[] = "XYZW";
  *c->verbose << "Extracting slice " << slicepos << " along " << axis[slicedir] 
    << " axis in image #" << c->m_ImageStack.size()-1 << endl;

  // Use the extractor to extract the actual region
  RegionType rslice = image->GetBufferedRegion();
  rslice.SetSize(slicedir, 1);
  rslice.SetIndex(slicedir, slicepos);

  // Push the big image back on the stack 
  c->m_ImageStack.push_back(image);

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
}

template <class TPixel, unsigned int VDim>
void
ExtractSlice<TPixel, VDim>
::operator() (string axis, char* positions, unsigned int n_comp)
{
  // Check input availability
  if(c->m_ImageStack.size() < n_comp)
    throw ConvertException("At least %d images are required on the stack for the -slice or -slice-all command", n_comp);

  // Check that the images are the same size
  if(n_comp > 1 && !c->CheckStackSameDimensions(n_comp))
    throw ConvertException("All images must be the same size for the -slice-all command");

  // Get the last image for the size command considerations
  SizeType size = c->PeekImage(0)->GetBufferedRegion().GetSize();

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
    throw ConvertException("first parameter to -slice or -slice-all must be one of either 0,1,2,3 or x,y,z,w/t");

  // Now determine the pattern of the second parameter. Allowed formats are
  // 1. 15      // slice 15 (0-based indexing)
  // 2. -1      // slice N-1
  // 3. 20%     // slice at 20% of the stack
  // 4. 12:-4   // range, every slice
  // 5. 12:3:-4 // range, every third slice
  // 6. 9:7.5:39 // range, equals: -slice 9 17 24 32 39
  // 7. 40%:5%:60% // range, equals: -slice 40% 45% 50% 55% 60%
  
  // Split the string on the ':'
  std::vector< std::vector<double> > pos_lists;
  std::string position;
  std::stringstream source_all(positions);
  while(std::getline(source_all, position, ','))
  {
  *c->verbose << "Processing slice string: " << position << endl;

  std::string piece;
  std::stringstream source(position);
  std::vector<double> pos_list;
  while(std::getline(source, piece, ':'))
    {
    double slicepos;

    // Process the percentage
    if(piece[piece.size()-1] == '%')
      {
      piece = piece.substr(0, piece.size()-1);
      double percent_pos = atof(piece.c_str());
      slicepos = (percent_pos / 100.0) * (size[slicedir] -1);
      *c->verbose << "Calculated position in percent: " << percent_pos << "% of " << size[slicedir] << " to be image plane position: " << slicepos << endl;
      }
    else
      {
      slicepos = atof(piece.c_str());
      }

    pos_list.push_back(slicepos);
    }
    pos_lists.push_back(pos_list);
  }

  std::vector<int> slice_extraction_list;
  for(int l= 0; l < pos_lists.size(); l++)
  {
  std::vector<double> pos_list = pos_lists[l];

  // Now we have one, two or three numbers parsed
  double pos_first, pos_step = 1, pos_last;
  if(pos_list.size() == 1)
    {
    pos_first = pos_last = pos_list[0];
    }
  else if(pos_list.size() == 2)
    {
    pos_first = pos_list[0];
    pos_last = pos_list[1];
    }
  else if(pos_list.size() == 3)
    {
    pos_first = pos_list[0];
    pos_step = pos_list[1];
    pos_last = pos_list[2];
    }

  // Deal with negative values for start/stop
  if(pos_first < 0)
    pos_first = size[slicedir] + pos_first;

  if(pos_last < 0)
    pos_last = size[slicedir] + pos_last;

  // Deal with the sign of the step
  if((pos_first < pos_last && pos_step < 0)
    || (pos_first > pos_last && pos_step > 0))
    {
    pos_step *= -1;
    }

  // Make sure all is legit
  if(pos_first < pos_last && pos_step <= 0)
    throw ConvertException(
      "Wrong slice list specification %d:%d:%d for -slice or -slice-all command! Step should be positive.",
      pos_first, pos_step, pos_last);

  if(pos_first > pos_last && pos_step >= 0)
    throw ConvertException(
      "Wrong slice list specification %d:%d:%d for -slice or -slice-all command! Step should be negative.",
      pos_first, pos_step, pos_last);

  // Store slice extration positions
  for(double pos = pos_first; pos <= pos_last; pos+=pos_step)
    {
    int i = (int)(0.5+pos);  // round to integer
    *c->verbose << "Calculated image plane position: " << pos << ", which is rounded to slice: " << i << endl;
    slice_extraction_list.push_back(i);
    }
  }

  // Take the last n_comp images on the stack
  std::vector<ImagePointer> comp_ptr = c->PopNImages(n_comp);

  // Now extract each slice
  for(int s = 0; s < slice_extraction_list.size(); s++)
  {
    int i = slice_extraction_list[s];
    for(int k = 0; k < n_comp; k++)
      this->ExtractOneSlice(comp_ptr[k], slicedir, i);
  }
}

// Invocations
template class ExtractSlice<double, 3>;
template class ExtractSlice<double, 4>;
template class ExtractSlice<double, 2>;

