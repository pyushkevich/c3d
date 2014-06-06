/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SampleImage.cxx
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

#include "SampleImage.h"

template <class TPixel, unsigned int VDim>
void
SampleImage<TPixel, VDim>
::operator() (const RealVector &x)
{
  // Create a point
  itk::Point<double, VDim> p;
  for(size_t i = 0; i < VDim; i++)
    p[i] = x[i];

  // Map the point to a continuous index
  ImagePointer image = c->m_ImageStack.back();
  itk::ContinuousIndex<double, VDim> cidx;
  image->TransformRASPhysicalPointToContinuousIndex(p, cidx);

  // Describe what we are doing
  *c->verbose << "Probing #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Physical (RAS) Coordinates: " << p << endl;
  *c->verbose << "  Voxel Coordinates         : " << cidx << endl;

  // Use the interpolator in c
  c->GetInterpolator()->SetInputImage(image);
  m_Result = c->GetInterpolator()->EvaluateAtContinuousIndex(cidx);
  *c->verbose << "  Using " << c->m_Interpolation << " interpolation" << endl;

  // Print out the interpolated value
  cout << "Interpolated image value at " << x << " is " << m_Result << endl;
}

// Invocations
template class SampleImage<double, 2>;
template class SampleImage<double, 3>;
template class SampleImage<double, 4>;
