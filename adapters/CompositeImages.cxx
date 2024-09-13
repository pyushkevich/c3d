/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    CompositeImages.cxx
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

#include "CompositeImages.h"
#include "itkBinaryFunctorImageFilter.h"

template <class TInputPixel, class TOutputPixel>
class ComposeFunctor
{
public:

  ComposeFunctor(TInputPixel background = 0.) : m_Background(background) {}

  // The first parameter is the image, and the second parameter is the mask
  TOutputPixel operator() (const TInputPixel &a, const TInputPixel &b)
  {
    return b == m_Background ? a : b;
  }

  bool operator != (const ComposeFunctor<TInputPixel, TOutputPixel> &other)
  { return m_Background != other.m_Background; }
private:
  TInputPixel m_Background;
};


template <class TPixel, unsigned int VDim>
void
CompositeImages<TPixel, VDim>
::operator() ()
{
  // Get the last two images from the stack
  if(c->m_ImageStack.size() < 2)
    throw ConvertException("Composite operation requires two images on the stack");
  ImagePointer img_b = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();
  ImagePointer img_a = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();

  *c->verbose << "Compositing " << c->m_ImageStack.size()+2 << " onto " << c->m_ImageStack.size()+1 << std::endl;

  // Do some processing ...
  typedef ComposeFunctor<typename ImageType::PixelType, typename ImageType::PixelType> Functor;
  typedef itk::BinaryFunctorImageFilter<ImageType, ImageType, ImageType, Functor> Filter;
  Functor f(c->m_Background);
  typename Filter::Pointer filter = Filter::New();
  filter->SetFunctor(f);
  filter->SetInput1(img_a);
  filter->SetInput2(img_b);
  filter->Update();

  ImagePointer result = filter->GetOutput();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class CompositeImages<double, 2>;
template class CompositeImages<double, 3>;
template class CompositeImages<double, 4>;
