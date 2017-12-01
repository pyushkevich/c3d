/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    RFTrain.cxx
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

#include "RFTrain.h"
#include "ImageCollectionToImageFilter.h"
#include "itkVectorImage.h"
#include "Library/classifier.h"
#include "Library/classification.h"
#include "RandomForestClassifier.h"
#include <iostream>

template <class TPixel, unsigned int VDim>
void
RFTrain<TPixel, VDim>
::operator() (const char * train_file, const ParametersType &param)
{
  // There should be at least two images on the stack
  if(c->m_ImageStack.size() < 2)
    throw ConvertException("At least two images (data/label) required to train RF classifier");

  // Get the label image
  ImagePointer imgLabel = c->m_ImageStack.back();

  // Iterator for grouping images into a multi-component image
  typedef itk::VectorImage<TPixel, VDim> VectorImageType;
  typedef ImageCollectionConstRegionIteratorWithIndex<
      ImageType, VectorImageType> CollectionIter;

  // Get the segmentation image - which determines the samples
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> LabelIter;

  // Shrink the buffered region by radius because we can't handle BCs
  itk::ImageRegion<VDim> reg = imgLabel->GetBufferedRegion();
  reg.ShrinkByRadius(param.patch_radius);

  // We need to iterate throught the label image once to determine the
  // number of samples to allocate.
  unsigned long nSamples = 0;
  for(LabelIter lit(imgLabel, reg); !lit.IsAtEnd(); ++lit)
    if( (int) (0.5 + lit.Value()) > 0)
      nSamples++;

  // Create an iterator for going over all the anatomical image data
  CollectionIter cit(reg);
  cit.SetRadius(param.patch_radius);

  // Add all the anatomical images to this iterator
  for(int i = 0; i < c->m_ImageStack.size() - 1; i++)
    cit.AddImage(c->m_ImageStack[i]);

  // Get the number of components
  int nComp = cit.GetTotalComponents();
  int nPatch = cit.GetNeighborhoodSize();
  int nColumns = nComp * nPatch;

  // Are we using coordinate informtion
  if(param.use_coordinate_features)
    nColumns += VDim;

  // Create a new sample
  typedef MLData<TPixel, TPixel> SampleType;
  SampleType *sample = new SampleType(nSamples, nColumns);

  // Now fill out the samples
  int iSample = 0;
  for(LabelIter lit(imgLabel, reg); !lit.IsAtEnd(); ++lit, ++cit)
    {
    int label = (int) (lit.Value() + 0.5);
    if(label > 0)
      {
      // Fill in the data
      std::vector<TPixel> &column = sample->data[iSample];
      int k = 0;
      for(int i = 0; i < nComp; i++)
        for(int j = 0; j < nPatch; j++)
          column[k++] = cit.NeighborValue(i,j);

      // Add the coordinate features if used
      if(param.use_coordinate_features)
        for(int d = 0; d < VDim; d++)
          column[k++] = lit.GetIndex()[d];

      // Fill in the label
      sample->label[iSample] = label;
      ++iSample;
      }
    }

  // Check that the sample has at least two distinct labels
  bool isValidSample = false;
  for(int iSample = 1; iSample < sample->Size(); iSample++)
    if(sample->label[iSample] != sample->label[iSample-1])
      { isValidSample = true; break; }

  // Now there is a valid sample. The text task is to train the classifier
  if(!isValidSample)
    throw ConvertException("A classifier cannot be trained because the training "
      "data contain fewer than two classes. Please label "
      "examples of two or more tissue classes in the image.");

  // Set up the classifier parameters
  TrainingParameters params;

  // TODO:
  params.treeDepth = param.tree_depth;
  params.treeNum = param.forest_size;
  params.candidateNodeClassifierNum = 10;
  params.candidateClassifierThresholdNum = 10;
  params.subSamplePercent = 0;
  params.splitIG = 0.1;
  params.leafEntropy = 0.05;
  params.verbose = true;

  // Cap the number of training voxels at some reasonable number
  if(sample->Size() > 10000)
    params.subSamplePercent = 100 * 10000.0 / sample->Size();
  else
    params.subSamplePercent = 0;

  // Create the classification engine
  typedef RandomForestClassifier<TPixel, TPixel, VDim> RFClassifierType;
  typedef typename RFClassifierType::RFAxisClassifierType RFAxisClassifierType;
  typedef Classification<TPixel, TPixel, RFAxisClassifierType> ClassificationType;

  typename RFClassifierType::Pointer classifier = RFClassifierType::New();
  ClassificationType classification;

  // Perform classifier training
  classification.Learning(
        params, *sample,
        *classifier->GetForest(),
        classifier->GetValidLabel(),
        classifier->GetClassToLabelMapping());

  // Reset the class weights to the number of classes and assign default
  int n_classes = classifier->GetClassToLabelMapping().size(), n_fore = 0, n_back = 0;
  classifier->GetClassWeights().resize(n_classes, -1.0);

  // Store the patch radius in the classifier - this remains fixed until
  // training is repeated
  classifier->SetPatchRadius(param.patch_radius);
  classifier->SetUseCoordinateFeatures(param.use_coordinate_features);

  // Dump the classifier to a file
  ofstream out_file(train_file);
  classifier->Write(out_file);
  out_file.close();
}

// Invocations
template class RFTrain<double, 2>;
template class RFTrain<double, 3>;
template class RFTrain<double, 4>;
