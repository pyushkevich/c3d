/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    StapleAlgorithm.cxx
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

#include "StapleAlgorithm.h"
#include "itkSTAPLEImageFilter.h"

template <class TPixel, unsigned int VDim>
void
StapleAlgorithm<TPixel, VDim>
::operator() (double ival)
{
  size_t i;
  
  // Create a STAPLE filter
  typedef itk::STAPLEImageFilter<ImageType, ImageType> Filter;
  typename Filter::Pointer fltStaple = Filter::New();

  // Add each of the images on the stack to the filter
  for(i = 0; i < c->m_ImageStack.size(); i++)
    fltStaple->SetInput(i, c->m_ImageStack[i]);

  // Configure the STAPLE filter
  fltStaple->SetForegroundValue(ival);

  // Describe what we are doing
  *c->verbose << "Executing STAPLE EM Algorithm on " << c->m_ImageStack.size() << " images." << endl;

  // Plug in the filter's components
  fltStaple->Update();

  // Dump sensitivity/specificity values
  *c->verbose << "  Elapsed Iterations: " << fltStaple->GetElapsedIterations() << endl;
  for(i = 0; i < c->m_ImageStack.size(); i++)
    {
    *c->verbose << "  Rater " << i 
      << ": Sensitivity = " << fltStaple->GetSensitivity(i) 
      << "; Specificity = " << fltStaple->GetSpecificity(i) << endl;
    }

  // Store the output
  c->m_ImageStack.clear();
  c->m_ImageStack.push_back(fltStaple->GetOutput());

}

// Invocations
template class StapleAlgorithm<double, 2>;
template class StapleAlgorithm<double, 3>;
template class StapleAlgorithm<double, 4>;
