/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SplitMultilabelImage.cxx
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

#include "SplitMultilabelImage.h"
#include "ThresholdImage.h"
#include "UpdateMetadataKey.h"
#include <fstream>
#include <set>
#include "vnl/vnl_math.h"

using namespace std;

template <class TPixel, unsigned int VDim>
void
SplitMultilabelImage<TPixel, VDim>
::operator() (std::vector<double> values)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // List of labels to split into
  set<double> sval;

  // Create a list of all the finite values in the image or get from supplied list
  if(values.size() > 0)
    for(auto &it : values) sval.insert(it);
  else
    for(ConstIterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
      {
      double val = it.Get();
      if(std::isfinite(val))
        sval.insert(val);
      }
  
  // The number of finite values should be reasonable
  if(sval.size() > 256)
    throw ConvertException("Number of labels passed on to -split exceeds 255");
  else if(sval.size() == 0)
    throw ConvertException("No finite labels passed on to -split");

  // A report
  *c->verbose << "Splitting #" << c->m_ImageStack.size() << " into " 
    << sval.size() << " binary images" << endl;

  // Pop the image
  c->m_ImageStack.pop_back();

  // Clear the label set
  c->GetSplitLabelSet();

  // Generate a bunch of binary copies (unfortunately, we store them as double)
  for(set<double>::iterator it = sval.begin(); it != sval.end(); ++it)
    {
    // Push our image back on the stack
    c->m_ImageStack.push_back(img);

    // Perform thresholding
    ThresholdImage<TPixel, VDim> thresh(c);
    thresh(*it, *it, 1.0, 0.0);

    // Add label to the label set
    c->GetSplitLabelSet().push_back(*it);
    }  
}

// Invocations
template class SplitMultilabelImage<double, 2>;
template class SplitMultilabelImage<double, 3>;
template class SplitMultilabelImage<double, 4>;
