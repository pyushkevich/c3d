/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConnectedComponents.cxx
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

#include "ConnectedComponents.h"
#include "ThresholdImage.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkCastImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ConnectedComponents<TPixel, VDim>
::operator() ()
{
  // The image is assumed to be binary. If background is non-zero, call binarize
  // to map the background to zero
  if(c->m_Background != 0.0)
    {
    ThresholdImage<TPixel, VDim> thresh(c);
    thresh(c->m_Background, c->m_Background, 0.0, 1.0);
    }

  // Get the last image on the stack
  ImagePointer image = c->m_ImageStack.back();

  // Integer image typedef
  typedef itk::OrientedRASImage<int, VDim> IntImageType;
  
  // Construct the connected components filter
  typedef itk::ConnectedComponentImageFilter<ImageType, IntImageType> CCFilter;

  // Relabel the components
  typedef itk::RelabelComponentImageFilter<IntImageType, IntImageType> RCFilter;

  // Describe what we are doing
  *c->verbose << "Computing connected components of #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Calling ConnectedComponentImageFilter" << endl;

  // Plug in the filter's components
  typename CCFilter::Pointer fltConnect = CCFilter::New();
  fltConnect->SetInput(image);
  fltConnect->SetFullyConnected(false);
  fltConnect->Update();

  // Describe what we are doing
  *c->verbose << "  Calling RelabelComponentImageFilter" << endl;

  // Relabel and order components
  typename RCFilter::Pointer fltRelabel = RCFilter::New();
  fltRelabel->SetInput(fltConnect->GetOutput());
  fltRelabel->Update();

  // Print the statistics about the connected components
  long szpx = fltRelabel->GetSizeOfObjectInPixels(1);
  cout << "  There are " << 
    fltRelabel->GetNumberOfObjects() << " connected components." << endl;
  cout << "  Largest component has " << szpx << " pixels." << endl;

  // We have to convert the image back to the native type
  typedef itk::CastImageFilter<IntImageType, ImageType> CastFilter;
  typename CastFilter::Pointer fltCast = CastFilter::New();
  fltCast->SetInput(fltRelabel->GetOutput());
  fltCast->Update();

  // Store the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltCast->GetOutput());
}

// Invocations
template class ConnectedComponents<double, 2>;
template class ConnectedComponents<double, 3>;
template class ConnectedComponents<double, 4>;
