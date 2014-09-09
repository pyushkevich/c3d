/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    WarpImage.cxx
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

#include "WarpImage.h"
#include "itkWarpImageFilter.h"
#include "CreateInterpolator.h"

template <class TPixel, unsigned int VDim>
void
WarpImage<TPixel, VDim>
::operator() ()
{
  // Check input availability
  if(c->m_ImageStack.size() < VDim + 1)
    {
    cerr << "Warp operation requires " << VDim+1 << " images on the stack" << endl;
    throw -1;
    }

  // Write something 
  *c->verbose << "Warping image #" << c->m_ImageStack.size() << endl;

  // Get the image to warp
  ImagePointer isrc = c->m_ImageStack[c->m_ImageStack.size() - 1];

  // Store the direction temporarily
  // itk::Matrix<double, VDim, VDim> dirin = isrc->GetDirection(), dirtemp;
  // dirtemp.SetIdentity();
  // isrc->SetDirection(dirtemp);

  // Index of the first warp image
  size_t iwarp = c->m_ImageStack.size() - (VDim + 1);

  // Create a deformation field
  typedef itk::Vector<TPixel, VDim> VectorType;
  typedef itk::OrientedRASImage<VectorType, VDim> FieldType;
  typename FieldType::Pointer field = FieldType::New();
  field->CopyInformation(c->m_ImageStack[iwarp]);
  field->SetRegions(c->m_ImageStack[iwarp]->GetBufferedRegion());
  field->Allocate();
  size_t nvox = field->GetBufferedRegion().GetNumberOfPixels();
  for(size_t d = 0; d < VDim; d++)
    {
    ImagePointer warp = c->m_ImageStack[iwarp + d];
    if(warp->GetBufferedRegion() != field->GetBufferedRegion())
      throw ConvertException("Warp field components have different dimensions");
    for(size_t i = 0; i < nvox; i++)
      field->GetBufferPointer()[i][d] = warp->GetBufferPointer()[i];
    }

  // Create the warp filter
  typedef itk::WarpImageFilter<ImageType, ImageType, FieldType> WarpType;
  typename WarpType::Pointer fltWarp = WarpType::New();
  fltWarp->SetInput(isrc);
  fltWarp->SetDisplacementField(field);

  // Create interpolator
  fltWarp->SetInterpolator(c->GetInterpolator());

  // Update the warp fileter
  fltWarp->SetOutputSpacing(field->GetSpacing());
  fltWarp->SetOutputOrigin(field->GetOrigin());
  fltWarp->SetOutputDirection(field->GetDirection());
  fltWarp->SetEdgePaddingValue(c->m_Background);
  fltWarp->Update();

  // Update the output image
  //
  ImagePointer imgout = fltWarp->GetOutput();

  // Drop image and warps from stack
  c->m_ImageStack.pop_back();
  for(size_t d = 0; d < VDim; d++)
    c->m_ImageStack.pop_back();

  // Add the warped image to the stack
  c->m_ImageStack.push_back(imgout);
}

// Invocations
template class WarpImage<double, 2>;
template class WarpImage<double, 3>;
template class WarpImage<double, 4>;
