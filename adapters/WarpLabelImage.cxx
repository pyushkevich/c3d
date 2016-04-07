/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    WarpLabelImage.cxx
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

#include "WarpLabelImage.h"
#include "itkWarpImageFilter.h"
#include "CreateInterpolator.h"
#include "ThresholdImage.h"
#include "SmoothImage.h"
#include <set>

template <class TPixel, unsigned int VDim>
void
WarpLabelImage<TPixel, VDim>
::operator() (RealVector &stdev)
{
  // Check input availability
  if(c->m_ImageStack.size() < VDim + 1)
    {
    cerr << "Warp operation requires " << VDim+1 << " images on the stack" << endl;
    throw -1;
    }

  // Write something 
  *c->verbose << "Warping image label-wise #" << c->m_ImageStack.size() << endl;

  // Index of the first warp image
  size_t iwarp = c->m_ImageStack.size() - (VDim + 1);


  // Get the image to warp
  ImagePointer isrc = c->m_ImageStack[c->m_ImageStack.size() - 1];

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
  fltWarp->SetDisplacementField(field);

  // Create interpolator
  fltWarp->SetInterpolator(c->GetInterpolator());

  // Update the warp fileter
  fltWarp->SetOutputSpacing(field->GetSpacing());
  fltWarp->SetOutputOrigin(field->GetOrigin());
  fltWarp->SetOutputDirection(field->GetDirection());
  fltWarp->SetEdgePaddingValue(c->m_Background);

  // Find all unique labels in the input image
  std::set<TPixel> labels;
  size_t n = isrc->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < n; i++)
    labels.insert(isrc->GetBufferPointer()[i]);

  // Create the output image
  ImagePointer ilabel = ImageType::New();
  ilabel->CopyInformation(field);
  ilabel->SetRegions(field->GetBufferedRegion());
  ilabel->Allocate();

  // Create the 'score' image
  ImagePointer iscore = ImageType::New();
  iscore->SetRegions(field->GetBufferedRegion());
  iscore->Allocate();
  iscore->FillBuffer(0.0);

  // For each unique label in the input image, threshold that image
  for(typename std::set<TPixel>::iterator it = labels.begin(); it!=labels.end(); ++it)
    {
    // Current label value
    TPixel val = *it;
    
    // Threshold the image (get only current label)
    ThresholdImage<TPixel, VDim> thresh(c);
    thresh(val, val, 1, 0);

    // Smooth the image
    SmoothImage<TPixel, VDim> smooth(c);
    smooth(stdev, true);

    // Apply the warp to the smoothed image
    fltWarp->SetInput(c->m_ImageStack.back());
    fltWarp->Update();

    // Iterate over voxels, update the score and label images
    ImagePointer iwarped = fltWarp->GetOutput();
    size_t nvox = iwarped->GetBufferedRegion().GetNumberOfPixels();
    TPixel *pscore = iscore->GetBufferPointer(), 
           *pwarped = iwarped->GetBufferPointer(), 
           *plabel = ilabel->GetBufferPointer();
    for(size_t i = 0; i < nvox; i++)
      {
      if(pscore[i] < pwarped[i])
        {
        pscore[i] = pwarped[i];
        plabel[i] = val;
        }
      }

    // Replace everything on the stack
    c->m_ImageStack.pop_back();
    c->m_ImageStack.push_back(isrc);
    }

  // Drop image and warps from stack
  c->m_ImageStack.pop_back();
  for(size_t d = 0; d < VDim; d++)
    c->m_ImageStack.pop_back();

  // Store the result on the stack
  c->m_ImageStack.push_back(ilabel);
}

// Invocations
template class WarpLabelImage<double, 2>;
template class WarpLabelImage<double, 3>;
template class WarpLabelImage<double, 4>;
