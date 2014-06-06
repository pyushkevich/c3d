/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ReplaceIntensities.cxx
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

#include "ReplaceIntensities.h"

template <class TPixel, unsigned int VDim>
void
ReplaceIntensities<TPixel, VDim>
::operator() (vector<double> &xRule)
{
  // Get the input image
  ImagePointer input = c->m_ImageStack.back();

  // This is a slightly sensitive procedure. We need to set an epsilon so that
  // intensities that are within machine precision are treated as equal. We define
  // this as 2 (a - b) / (a + b) < eps.
  double epsilon = 0.000001;

  // Report what we are doing
  *c->verbose << "Replacing intensities in #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Replacement Rules: ";
  for(size_t i = 0; i < xRule.size(); i+=2)
    *c->verbose << xRule[i] << " -> " << xRule[i+1] << "; ";
  *c->verbose << endl;

  // Create an iterator to process the image
  long nReps = 0;
  for(Iterator it(input, input->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    // Get the pixel value
    TPixel val = it.Get();

    // Repeat over the rules
    for(size_t k = 0; k < xRule.size(); k += 2)
      {
      double u = xRule[k], v = xRule[k+1];
      if(
        (vnl_math_isnan(val) && vnl_math_isnan(u)) ||
        (val == u) ||
        fabs(2*(val-u)/(val+u)) < epsilon)
        {
        it.Set(v);
        nReps++;
        break;
        }
      }
    }

  // Report what the filter is doing
  *c->verbose << "  Replacements Made: " << nReps << endl;

}

// Invocations
template class ReplaceIntensities<double, 2>;
template class ReplaceIntensities<double, 3>;
template class ReplaceIntensities<double, 4>;
