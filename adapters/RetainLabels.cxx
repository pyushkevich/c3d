/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    RetainLabels.cxx
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

#include "RetainLabels.h"
#include <itkUnaryFunctorImageFilter.h>

template <class TPixel, unsigned int VDim>
class RetainLabelsUnaryFunctor
{
public:
  RetainLabelsUnaryFunctor(const std::vector<double> &retain_list, TPixel bkg, bool mask, TPixel mask_value)
    : m_RetainList(retain_list), m_HaveLast(false), m_Background(bkg), m_Mask(mask), m_MaskValue(mask_value) {}

  RetainLabelsUnaryFunctor()
    : m_HaveLast(false), m_Background(0.0) {}

  TPixel operator() (const TPixel &pixel)
    {
    // Allow caching because most of our inputs are going to be piecewise constant
    if(m_HaveLast && pixel == m_LastPixel)
      return m_LastReplacement;

    // Save the last pixel value
    m_LastPixel = pixel;
    m_HaveLast = true;

    // Map the pixel to int
    int p_int = (int)(pixel + 0.5);

    // Check the pixel against the retain list
    for(int i = 0; i < m_RetainList.size(); i++)
      {
      if(p_int == m_RetainList[i])
        {
        m_LastReplacement = m_Mask ? m_MaskValue : pixel;
        return m_LastReplacement;
        }
      }

    // Failed check
    m_LastReplacement = m_Background;
    return m_Background;
    }

  bool operator != (const RetainLabelsUnaryFunctor<TPixel, VDim> &other)
    {
    if(m_RetainList != other.m_RetainList) return true;
    if(m_LastPixel != other.m_LastPixel) return true;
    if(m_LastReplacement != other.m_LastReplacement) return true;
    if(m_Background != other.m_Background) return true;
    if(m_HaveLast != other.m_HaveLast) return true;
    return false;
    }

protected:
  
  std::vector<double> m_RetainList;
  TPixel m_LastPixel, m_LastReplacement, m_Background, m_MaskValue;
  bool m_HaveLast, m_Mask;
};

template <class TPixel, unsigned int VDim>
void
RetainLabels<TPixel, VDim>
::operator() (const RetainLabels<TPixel, VDim>::LabelSet &retain, bool mask, double mask_value)
{
  // Get image from stack
  ImagePointer img = c->PopImage();

  // Create the functor
  *c->verbose << (mask ? "Masking label(s) " : "Retaining label(s) ");
  for (auto &l : retain)
    *c->verbose << l << " ";
  if(mask)
    *c->verbose << "with label " << mask_value << " ";
  *c->verbose << "in #" << c->m_ImageStack.size() << std::endl;

  typedef RetainLabelsUnaryFunctor<TPixel, VDim> FunctorType;
  FunctorType functor(retain, c->m_Background, mask, mask_value);

  // Create the filter
  typedef itk::UnaryFunctorImageFilter<ImageType, ImageType, FunctorType> FilterType;
  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(img);
  filter->SetFunctor(functor);
  filter->Update();

  // Take the output
  c->PushImage(filter->GetOutput());
}

// Invocations
template class RetainLabels<double, 2>;
template class RetainLabels<double, 3>;
template class RetainLabels<double, 4>;
