/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SmoothMultiLabelImage.cxx
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

#include "SmoothMultiLabelImage.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"

// This functor does pixel-wise comparison between Current Image and Max Image
// if crntPixel > maxPixel, return crntPixel as the new maxPixel
class BinaryIntensityVotingFunctor
{
public:
  typedef BinaryIntensityVotingFunctor Self;

  BinaryIntensityVotingFunctor() = default;
  ~BinaryIntensityVotingFunctor() = default;

  double operator() (const double &crntPixel, const double &maxPixel) const
  {
    return crntPixel > maxPixel ? crntPixel : maxPixel;
  }
};

// This functor does pixel-wise comparison between Current Image and Max Image
// if crntPixel > maxPixel, return current LABEL as result, otherwise 0
class BinaryLabelVotingFunctor
{
public:
  typedef BinaryLabelVotingFunctor Self;

  BinaryLabelVotingFunctor(double crntLabel)
    : m_crntLabel(crntLabel) {};

  // default constructor is needed for declaration of the functor filter
  BinaryLabelVotingFunctor() {};
  ~BinaryLabelVotingFunctor() = default;

  double operator() (const double &crntPixel, const double &maxPixel) const
  {
    return crntPixel > maxPixel ? m_crntLabel : 0;
  }

  bool operator != (const Self &other)
  {
    return other.m_crntLabel != m_crntLabel;
  }

private:
  double m_crntLabel;
};

// This functor does pixel-wise comparison between current label result and global result
// if current result is 0, return global, otherwise return current result
template<class TPixel>
class BinaryLabelDeterminationFunctor
{
public:
  typedef BinaryLabelDeterminationFunctor Self;

  BinaryLabelDeterminationFunctor() = default;
  ~BinaryLabelDeterminationFunctor() = default;

  TPixel operator() (const double &crntLabelResult, const TPixel &globalLabelResult) const
  {
    return crntLabelResult == 0 ? globalLabelResult : crntLabelResult;
  }
};

// A utility method to deep copy image
template<typename TImage>
void DeepCopy(typename TImage::Pointer input, typename TImage::Pointer output)
{
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();

  itk::ImageRegionConstIterator<TImage> inputIterator(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TImage>      outputIterator(output, output->GetLargestPossibleRegion());

  while (!inputIterator.IsAtEnd())
  {
    outputIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++outputIterator;
  }
}

template <class TPixel, unsigned int VDim>
void
SmoothMultiLabelImage<TPixel, VDim>::operator() 
(RealVector &stdev, std::vector<unsigned short> &labelsToSmooth)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();
  constexpr unsigned int sigmaDim = VDim > 3 ? 3 : VDim;

  // Process labelsToSmooth array
  // eliminate duplication and add background
  std::set<double> smoothingSet{0.0};

  set<TPixel> label_set;
  for(ConstIterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    TPixel val = it.Get();
    if(std::isfinite(val))
      label_set.insert(val);
    }

  if (labelsToSmooth.size() == 0)
  {
    // Smooth all labels
    *c->verbose << "All labels will be smoothed" << std::endl;
    for (auto cit = label_set.cbegin(); cit != label_set.cend(); ++cit)
      smoothingSet.insert(*cit);
  }
  else
  {
    for (auto cit = labelsToSmooth.cbegin(); cit != labelsToSmooth.cend(); ++cit)
      smoothingSet.insert((double)*cit);
  }

  *c->verbose << "Smoothing standard deviation (mm): (";
  for (size_t i = 0; i < sigmaDim; ++i)
  {
    *c->verbose << stdev[i];
    if (i < sigmaDim - 1)
      *c->verbose << ',';
    else
      *c->verbose << ')';
  }
    

  typedef itk::Image<double, VDim> DoubleImageType;
  typedef itk::BinaryThresholdImageFilter<ImageType, DoubleImageType> ThresholdFilterType;
  typedef itk::BinaryThresholdImageFilter<DoubleImageType, ImageType> WinningLabelImageGeneratorType;
  
  // Generate a blank image for recording maximum intensities
  typename ThresholdFilterType::Pointer miGen = ThresholdFilterType::New();
  miGen->SetInput(img);
  miGen->SetLowerThreshold(0.0);
  miGen->SetUpperThreshold(0.0);
  miGen->SetInsideValue(0.0);
  miGen->SetOutsideValue(0.0);
  miGen->Update();

  typename DoubleImageType::Pointer maxIntensityImg = miGen->GetOutput();
  typename DoubleImageType::Pointer newMaxIntensityImg = miGen->GetOutput();

  // Generate a blank image for recording winning labels
  typename WinningLabelImageGeneratorType::Pointer wlGen = WinningLabelImageGeneratorType::New();
  wlGen->SetInput(img);
  wlGen->SetLowerThreshold(0.0);
  wlGen->SetUpperThreshold(0.0);
  wlGen->SetInsideValue(0);
  wlGen->SetOutsideValue(0);
  wlGen->Update();

  typename ImageType::Pointer winningLabelsImg = wlGen->GetOutput();

  // Intensity Voter record current maximum intensity for each pixel
  typedef itk::BinaryFunctorImageFilter
    <DoubleImageType, DoubleImageType, DoubleImageType, BinaryIntensityVotingFunctor> IntensityVoterType;
  typename IntensityVoterType::Pointer intensityVoter = IntensityVoterType::New();

  // Label Voter generate a temporary result of pixels at which current label is replacing the gloabl winning label
  typedef itk::BinaryFunctorImageFilter
    <DoubleImageType, DoubleImageType, DoubleImageType, BinaryLabelVotingFunctor> LabelVoterType;
  typename LabelVoterType::Pointer labelVoter = LabelVoterType::New();
  
  // Label Determinator uses the label voting result to refresh the global winning label result
  typedef BinaryLabelDeterminationFunctor<TPixel> LabelDeterminationFunctorType;
  typedef itk::BinaryFunctorImageFilter
    <DoubleImageType, ImageType, ImageType, LabelDeterminationFunctorType> LabelDeterminatorType;
  typename LabelDeterminatorType::Pointer labelDeterminator = LabelDeterminatorType::New();

  // The smoothing filter
  typedef itk::SmoothingRecursiveGaussianImageFilter<DoubleImageType, DoubleImageType> SmoothingFilterType;

  // Iterate through all labels for smoothing process
  for (auto cit = label_set.cbegin(); cit != label_set.cend(); ++cit)
  {
    *c->verbose << "Processing Label: " << *cit << std::endl;

    // Threshold current label to 1, rest 0
    typename ThresholdFilterType::Pointer fltThreshold = ThresholdFilterType::New();

    fltThreshold->SetInput(img);
    fltThreshold->SetLowerThreshold(*cit);
    fltThreshold->SetUpperThreshold(*cit);
    fltThreshold->SetInsideValue(1.0);
    fltThreshold->SetOutsideValue(0.0);
    fltThreshold->Update();

    typename DoubleImageType::Pointer startingImg = fltThreshold->GetOutput();

    // Do the smoothing, only for selected labels
    if (smoothingSet.count(*cit))
    {
      typename SmoothingFilterType::Pointer fltSmooth = SmoothingFilterType::New();

      // Add input
      fltSmooth->SetInput(startingImg);
      
      // Populate sigma array
      typename SmoothingFilterType::SigmaArrayType sigmaArr;

      for (size_t i = 0; i < sigmaDim; ++i)
        sigmaArr[i] = stdev[i];

      fltSmooth->SetSigmaArray(sigmaArr);

      // Update
      fltSmooth->Update();

      startingImg = fltSmooth->GetOutput();

      *c->verbose << "  Smoothing Completed" << std::endl;
    }

    // Start Voting Process

    /* Intensity Voting:
      A pixel-wise intensity comparison between current label smoothing result
      and the previous highest instensity. If current intensity is greater than previous high,
      it will replace previous high in the output image. Otherwise, previous high will be
      preserved in the output image. */

    intensityVoter->SetInput1(startingImg);
    intensityVoter->SetInput2(maxIntensityImg);
    intensityVoter->Update();

    // Temporarily save the new max intensity, since old maxIntensity needs to be used again
    newMaxIntensityImg = intensityVoter->GetOutput();

    /* Label Voting:
      A pixel-wise intensity comparison between current label smoothing result
      and the previous highest intensity. If current intensity is greater than previous high,
      current label will be written to the output image. Otherwise, 0 will be written. */

    BinaryLabelVotingFunctor lvf(*cit);
    labelVoter->SetFunctor(lvf);
    labelVoter->SetInput1(startingImg);
    labelVoter->SetInput2(maxIntensityImg);
    labelVoter->Update();
    typename DoubleImageType::Pointer labelVotingResult = labelVoter->GetOutput();

    /* Label Determination:
      Pixel-wise iteration on current label voting result and the global label voting result.
      If current result is non-zero, replace global result value with the current result value. */
    
    labelDeterminator->SetInput1(labelVotingResult);
    labelDeterminator->SetInput2(winningLabelsImg);
    labelDeterminator->Update();
    winningLabelsImg = labelDeterminator->GetOutput();

    *c->verbose << "  Voting Completed" << std::endl;

    // Now copy new intensity value into maxIntensityImg
    DeepCopy<DoubleImageType>(newMaxIntensityImg, maxIntensityImg);
  }
  
  ImagePointer result = winningLabelsImg;
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class SmoothMultiLabelImage<double, 2>;
template class SmoothMultiLabelImage<double, 3>;
template class SmoothMultiLabelImage<double, 4>;
