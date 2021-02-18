/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    MorphologicalContourInterpolation.cxx
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

#include "MorphologicalContourInterpolation.h"
#include "itkMorphologicalContourInterpolator.h"

template <class TPixel, unsigned int VDim>
void
MorphologicalContourInterpolation<TPixel, VDim>
::operator() (int axis, bool heuristic_alignment, bool use_distance_transform)
{
  if (axis < -1 || axis >= (int)VDim)
    throw ConvertException("MorphologicalContourInterpolation requires that axis is in [-1, %u[, got %d", VDim, axis);

  // Get image from stack
  ImagePointer label = c->m_ImageStack.back();

  // Create a short image for the labels
  typedef itk::Image<short, VDim> LabelImageType;
  typename LabelImageType::Pointer slab = LabelImageType::New();

  // Allocate the image
  slab->CopyInformation(label);
  slab->SetRegions(label->GetBufferedRegion());
  slab->Allocate();

  // Round off doubles to create labels
  size_t nv = label->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < nv; i++)
    slab->GetBufferPointer()[i] = (short) (label->GetBufferPointer()[i] + 0.5);

  // Define the filter
  typedef itk::MorphologicalContourInterpolator<LabelImageType> MCIType;

  // Create the filter
  typename MCIType::Pointer mci = MCIType::New();
  mci->SetInput(slab);
  mci->SetAxis(axis);
  mci->SetHeuristicAlignment(heuristic_alignment);
  mci->SetUseDistanceTransform(use_distance_transform);

  // Run the filter
  *c->verbose << "Interpolating contour from #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Axis: " << mci->GetAxis() << endl;
  *c->verbose << "  HeuristicAlignment: " << mci->GetHeuristicAlignment() << endl;
  *c->verbose << "  UseDistanceTransform: " << mci->GetUseDistanceTransform() << endl;
  mci->Update();

  // Convert result from short to double
  typename LabelImageType::Pointer sres = mci->GetOutput();
  ImagePointer result = ImageType::New();
  result->CopyInformation(sres);
  result->SetRegions(sres->GetBufferedRegion());
  result->Allocate();
  nv = sres->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < nv; i++)
    result->GetBufferPointer()[i] = (TPixel) (sres->GetBufferPointer()[i]);

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

template <>
void
MorphologicalContourInterpolation<double, 2>
::operator() (int axis, bool heuristic_alignment, bool use_distance_transform)
{
  throw ConvertException("MorphologicalContourInterpolation requires images with 3 dimensions or higher");
}

// Invocations
template class MorphologicalContourInterpolation<double, 2>;
template class MorphologicalContourInterpolation<double, 3>;
template class MorphologicalContourInterpolation<double, 4>;
