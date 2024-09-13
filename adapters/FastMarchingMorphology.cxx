/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    FastMarchingMorphology.cxx
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

#include "FastMarchingMorphology.h"
#include "RetainLabels.h"
#include "ThresholdImage.h"
#include "ReplaceIntensities.h"
#include "FastMarching.h"
#include "MultiplyImages.h"
#include "AddImages.h"
#include "CompositeImages.h"

template <class TPixel, unsigned int VDim>
void
FastMarchingMorphology<TPixel, VDim>::operator()(const LabelSet & active,
                                                 const LabelSet & target,
                                                 double           newlabel,
                                                 double           radius)
{
  // Take a copy of the image
  ImagePointer img = c->m_ImageStack.back();

  // Generate a binary image of the target (paint over) labels
  c->m_ImageStack.push_back(img);
  RetainLabels<TPixel, VDim>   mask_target(c);
  mask_target(target, true, 1);

  // Generate a binary image of the active labels
  c->m_ImageStack.push_back(img);
  RetainLabels<TPixel, VDim>   mask_active(c);
  mask_active(active, true, 1);

  // Perform fast marching and threshold in range (1, radius]
  FastMarching<TPixel, VDim>       fm(c);
  ReplaceIntensities<TPixel, VDim> replace_one_fm(c);
  ThresholdImage<TPixel, VDim>     thresh_fm(c);
  fm(radius);
  std::vector<double> rep_rule = { 1., 0. };
  replace_one_fm(rep_rule);
  thresh_fm(1, radius, 1, 0);
  ImagePointer mask = c->m_ImageStack.back();

  // The composite filter does not work here because we might be using
  // background as the dilation value. Instead we can manually composite
  ThresholdImage<TPixel, VDim>     thresh_1(c);
  MultiplyImages<TPixel, VDim>     mult_1(c);
  thresh_1(1, 1, 0, 1);
  mult_1();

  c->m_ImageStack.push_back(mask);
  ThresholdImage<TPixel, VDim>     thresh_2(c);
  thresh_2(1, 1, newlabel, 0);

  AddImages<TPixel, VDim>          add(c);
  add();

  // CompositeImages<TPixel, VDim> composite(c);
  // composite();
}

// Invocations
template class FastMarchingMorphology<double, 2>;
template class FastMarchingMorphology<double, 3>;
template class FastMarchingMorphology<double, 4>;
