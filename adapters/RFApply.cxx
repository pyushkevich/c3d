/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    RFApply.cxx
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

#include "RFApply.h"
#include "RandomForestClassifyImageFilter.h"
#include "itkVectorImage.h"
#include <iostream>

template <class TPixel, unsigned int VDim>
void
RFApply<TPixel, VDim>
::operator() (const char *train_file)
{
  // Copy all images on the stack to a temp array
  std::vector<ImagePointer> features;
  for(int k = 0; k < c->m_ImageStack.size(); k++)
    features.push_back(c->m_ImageStack[k]);
  c->m_ImageStack.clear();

  // Vector image support
  typedef itk::VectorImage<TPixel, VDim> VectorImageType;

  // Define the random forest classification filter
  typedef RandomForestClassifyImageFilter<ImageType, VectorImageType, ImageType, TPixel> FilterType;

  // Create a classifier object
  typedef typename FilterType::ClassifierType RFClassifierType;
  typename RFClassifierType::Pointer classifier = RFClassifierType::New();

  // Read the classifier object from disk
  ifstream in_file(train_file);
  classifier->Read(in_file);
  in_file.close();

  // Create the filter for this label (TODO: this is wasting computation)
  typename FilterType::Pointer filter = FilterType::New();

  // Add all the images on the stack to the filter
  for(int j = 0; j < features.size(); j++)
    filter->AddScalarImage(features[j]);

  // Pass the classifier to the filter
  filter->SetClassifier(classifier);

  // Set the filter behavior
  filter->SetGenerateClassProbabilities(true);

  // Run the filter for this set of weights
  filter->Update();

  // Append all the outputs to the stack
  for(int k = 0; k < filter->GetNumberOfIndexedOutputs(); k++)
    c->m_ImageStack.push_back(filter->GetOutput(k));

}

// Invocations
template class RFApply<double, 2>;
template class RFApply<double, 3>;
template class RFApply<double, 4>;
