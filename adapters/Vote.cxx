/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    Vote.cxx
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

#include "Vote.h"
#include "itkAddImageFilter.h"
#include "itkMetaDataObject.h"


template <class TPixel, unsigned int VDim>
void
Vote<TPixel, VDim>
::operator() (bool flagUseSplitLabelSet)
{
  // Create a maximum image
  ImagePointer i0 = c->m_ImageStack[0];

  // If this is in response to a split command, retrieve the label set
  typename Converter::LabelSet lset;
  if(flagUseSplitLabelSet)
    {
    if(c->m_ImageStack.size() != c->GetSplitLabelSet().size())
      throw ConvertException(
        "Merge failed: number of images on the stack (%i) different "
        "from the number of split labels (%i)", 
        c->m_ImageStack.size(), c->GetSplitLabelSet().size());
    lset = c->GetSplitLabelSet();
    }

  // Otherwise, the label mapping is identity
  else
    {
    for(size_t i = 0; i < c->m_ImageStack.size(); i++)
      lset.push_back((double) i);
    }

  // Create a vote image
  ImagePointer ilabel = ImageType::New();
  ilabel->SetRegions(i0->GetBufferedRegion());
  ilabel->SetOrigin(i0->GetOrigin());
  ilabel->SetSpacing(i0->GetSpacing());
  ilabel->SetDirection(i0->GetDirection());
  ilabel->Allocate();
  ilabel->FillBuffer(lset[0]);

  // Say something
  *c->verbose << "Collapsing " << c->m_ImageStack.size() << 
    " images into a multi-label image via maximum voting" << endl;

  // For each of the images, update the vote
  for(size_t j = 1; j < c->m_ImageStack.size(); j++)
    {
    // Get the next image pointer
    ImagePointer ij = c->m_ImageStack[j];

    // Check the image dimensions
    if(ij->GetBufferedRegion() != ilabel->GetBufferedRegion())
      throw ConvertException("All voting images must have same dimensions");

    // Get the label corresponding to this image
    double xVoteVal = lset[j];

    // Pairwise voting
    size_t n = ilabel->GetBufferedRegion().GetNumberOfPixels();
    for(size_t k = 0; k < n; k++)
      {
      if(ij->GetBufferPointer()[k] > i0->GetBufferPointer()[k])
        {
        i0->GetBufferPointer()[k] = ij->GetBufferPointer()[k];
        ilabel->GetBufferPointer()[k] = xVoteVal;
        }
      }
    }

  // Replace the images with the product
  c->m_ImageStack.clear();
  c->m_ImageStack.push_back(ilabel);
}

// Invocations
template class Vote<double, 2>;
template class Vote<double, 3>;
template class Vote<double, 4>;
