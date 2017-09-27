/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    FillBackgroundWithNeighborhoodNoise.cxx
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

#include "FillBackgroundWithNeighborhoodNoise.h"

#include "itkBinaryDilateImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkComposeImageFilter.h"
#include "OneDimensionalInPlaceAccumulateFilter.h"
#include <vnl/vnl_random.h>

template <class TVectorImage, class TMaskImage, class TNoiseImage>
class StatisticsToGaussianNoiseImageFilter 
  : public itk::ImageToImageFilter<TNoiseImage, TNoiseImage>
{
public:
  typedef StatisticsToGaussianNoiseImageFilter<TVectorImage, TMaskImage, TNoiseImage> Self;
  typedef itk::ImageToImageFilter<TNoiseImage, TNoiseImage> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  typedef TVectorImage StatisticsImageType;
  typedef TMaskImage MaskImageType;
  typedef typename StatisticsImageType::PixelType StatisticsPixelType;
  typedef TNoiseImage OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType OutputPixelType;

  itkTypeMacro(StatisticsToGaussianNoiseImageFilter, itk::ImageToImageFilter)

  itkNewMacro(Self) 

  void SetStatisticsImage(StatisticsImageType *image) 
    { itk::ProcessObject::SetInput("statistics", image); }

  void SetMaskImage(MaskImageType *image)
    { itk::ProcessObject::SetInput("mask", image); }

  void SetGrayImage(OutputImageType *image)
    { this->SetInput(image); }

  void ThreadedGenerateData(
    const OutputImageRegionType & outputRegionForThread,
    itk::ThreadIdType threadId) ITK_OVERRIDE
    {
    StatisticsImageType *stats = 
      dynamic_cast<StatisticsImageType *>(itk::ProcessObject::GetInput("statistics"));

    MaskImageType *mask = 
      dynamic_cast<MaskImageType *>(itk::ProcessObject::GetInput("mask"));

    const OutputImageType *gray = this->GetInput();

    OutputImageType *output = this->GetOutput();

    vnl_random rnd;
    itk::ImageRegionIterator<StatisticsImageType> itStat(stats, outputRegionForThread);
    itk::ImageRegionIterator<MaskImageType> itMask(mask, outputRegionForThread);
    itk::ImageRegionConstIterator<OutputImageType> itGray(gray, outputRegionForThread);
    itk::ImageRegionIterator<OutputImageType> itOut(output, outputRegionForThread);

    while(!itOut.IsAtEnd())
      {
      OutputPixelType pix = itGray.Value();
      if(itMask.Value())
        {
        const StatisticsPixelType &stats = itStat.Get();
        OutputPixelType n = stats[0], sum_x = stats[1], sum_x2 = stats[2];
        if(n > 1)
          {
          double mean_x = sum_x / n;
          double std_x = sqrt((sum_x2 - mean_x * sum_x) / (n - 1));
          double z = rnd.normal();
          pix = z * std_x + mean_x;
          }
        }

      itOut.Set(pix);
      ++itOut; ++itMask; ++itStat; ++itGray;
      }
    }

protected:

  StatisticsToGaussianNoiseImageFilter() {}
  ~StatisticsToGaussianNoiseImageFilter() {}
};

#include "itkImageFileWriter.h"

template <class TPixel, unsigned int VDim>
void
FillBackgroundWithNeighborhoodNoise<TPixel, VDim>
::operator() (const SizeType &radius, int n_steps)
{
  // Get the mask image
  ImagePointer mask = c->PopImage();
  ImagePointer gray = c->PopImage();

  // Create dilation element 
  typedef itk::BinaryBallStructuringElement<TPixel, VDim> StructuringElementType;
  typename StructuringElementType::SizeType seRadius;
  seRadius.Fill(3);
  StructuringElementType se;
  se.SetRadius(seRadius);
  se.CreateStructuringElement();
  
  // Repeat for the given number of steps
  for(int i = 0; i < n_steps; i++)
    {
    // Dilate the mask by 1 and take the difference - this is the 'new' interface
    typedef itk::BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType> DilatorType;
    typename DilatorType::Pointer dilator = DilatorType::New();
    dilator->SetInput(mask);
    dilator->SetKernel(se);
    dilator->SetDilateValue(1);

    // Perform dilation
    dilator->Update();
    ImagePointer newMask = dilator->GetOutput();

    // Subtract original mask
    typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType> SubtractorType;
    typename SubtractorType::Pointer sub = SubtractorType::New();
    sub->SetInput1(newMask);
    sub->SetInput2(mask);
    sub->Update();

    // Generate images of mask, mask * gray, mask * gray.^2
    typedef itk::MultiplyImageFilter<ImageType,ImageType,ImageType> MultiplierType;
    typename MultiplierType::Pointer mult_mg = MultiplierType::New();
    mult_mg->SetInput1(mask);
    mult_mg->SetInput2(gray);

    typename MultiplierType::Pointer mult_mg2 = MultiplierType::New();
    mult_mg2->SetInput1(mult_mg->GetOutput());
    mult_mg2->SetInput2(gray);

    // Compose these images into a vectorimage
    typedef itk::VectorImage<TPixel,VDim> VectorImageType;
    typedef itk::ComposeImageFilter<ImageType, VectorImageType> ComposerType;
    typename ComposerType::Pointer composer = ComposerType::New();
    composer->SetInput(0, mask);
    composer->SetInput(1, mult_mg->GetOutput());
    composer->SetInput(2, mult_mg2->GetOutput());
    composer->Update();
    typename VectorImageType::Pointer stats_input = composer->GetOutput();

    // Accumulate neighborhood sums of these values
    typename VectorImageType::Pointer accum = 
      AccumulateNeighborhoodSumsInPlace(stats_input.GetPointer(), radius);

    typedef itk::ImageFileWriter<ImageType> W;
    typename W::Pointer w = W::New();
    w->SetInput(newMask);
    w->SetFileName("/tmp/debug.nii.gz");
    w->Update();

    // Iterate, generating Gaussian noise (iteration is fine because most time will be
    // spent doing random number generation anyway, which threads poorly)
    typedef StatisticsToGaussianNoiseImageFilter
      <VectorImageType, ImageType, ImageType> NoiseFilter;
    typename NoiseFilter::Pointer noise = NoiseFilter::New();
    noise->SetStatisticsImage(accum);
    noise->SetMaskImage(sub->GetOutput());
    noise->SetGrayImage(gray);
    noise->Update();

    ImagePointer result = noise->GetOutput();

    // Update the mask and gray images
    result->DisconnectPipeline();
    newMask->DisconnectPipeline();

    gray = result;
    mask = newMask;
    }

  // Put result on stack
  c->PushImage(gray);
}

// Invocations
template class FillBackgroundWithNeighborhoodNoise<double, 2>;
template class FillBackgroundWithNeighborhoodNoise<double, 3>;
template class FillBackgroundWithNeighborhoodNoise<double, 4>;
