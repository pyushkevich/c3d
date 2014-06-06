/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ReciprocalImage.cxx
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

#include "ReciprocalImage.h"
#include "itkUnaryFunctorImageFilter.h"

template <class TIn, class TOut>
class ReciprocalFunctor
{
public:
  TOut operator() (const TIn &a)
    { return (TOut) (1.0 / a); }
};

template <class TPixel, unsigned int VDim>
void
ReciprocalImage<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Verbose
  *c->verbose << "Taking the reciprocal of #" << c->m_ImageStack.size() << endl;

  // Simply go through and divide
  typedef ReciprocalFunctor<TPixel, TPixel> FunctorType;
  typedef itk::UnaryFunctorImageFilter<ImageType, ImageType, FunctorType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput(img);
  filter->Update();

  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class ReciprocalImage<double, 2>;
template class ReciprocalImage<double, 3>;
template class ReciprocalImage<double, 4>;
