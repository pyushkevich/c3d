/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    NormalizedCrossCorrelation.cxx
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

#include "NormalizedCrossCorrelation.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkUnaryFunctorImageFilter.h"
#include "itkBinaryFunctorImageFilter.h"
#include "itkMeanImageFilter.h"

#include "itkImageRegionSplitterDirection.h"
#include "itkTimeProbe.h"

#include "itkImageFileWriter.h"

template <class TInputImage>
class OneDimensionalInPlaceAccumulateFilter : public itk::InPlaceImageFilter<TInputImage, TInputImage>
{
public:

  typedef OneDimensionalInPlaceAccumulateFilter<TInputImage> Self;
  typedef itk::InPlaceImageFilter<TInputImage, TInputImage> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkTypeMacro(OneDimensionalInPlaceAccumulateFilter, itk::InPlaceImageFilter);

  itkNewMacro(Self);

  /** Some convenient typedefs. */
  typedef TInputImage                          InputImageType;
  typedef TInputImage                          OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;

  /** We use a custom splitter */
  typedef itk::ImageRegionSplitterDirection    SplitterType;

  /** ImageDimension constant */
  itkStaticConstMacro(OutputImageDimension, unsigned int, TInputImage::ImageDimension);

  itkGetMacro(Radius, int);
  itkSetMacro(Radius, int);

  itkGetMacro(Dimension, int);
  itkSetMacro(Dimension, int);

protected:

  OneDimensionalInPlaceAccumulateFilter() { 
    m_Radius = 0; 
    m_Dimension = 0;
    m_Splitter = SplitterType::New(); 
    this->InPlaceOn();
  }
  
  ~OneDimensionalInPlaceAccumulateFilter() {}

  virtual void ThreadedGenerateData(
    const OutputImageRegionType & outputRegionForThread,
    itk::ThreadIdType threadId);

  virtual const itk::ImageRegionSplitterBase *GetImageRegionSplitter(void) const
    {
    m_Splitter->SetDirection(m_Dimension);
    return m_Splitter;
    }

  // Dimension of accumulation
  int m_Dimension;

  // Radius of accumulation
  int m_Radius;

  // Region splitter
  typename SplitterType::Pointer m_Splitter;

};

#include <itkImageLinearIteratorWithIndex.h>
#include <queue>

template <class TInputImage>
void
OneDimensionalInPlaceAccumulateFilter<TInputImage>
::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  itk::ThreadIdType threadId)
{
  // Get the image
  InputImageType *image = const_cast<InputImageType *>(this->GetInput());

  // Set up the iterator that will go through all the lines in the
  // output region. We assume that the lines span the whole length of
  // the input, i.e., the threading direction does not interfere
  typedef itk::ImageLinearIteratorWithIndex<TInputImage> IteratorType;
  IteratorType itScan(image, outputRegionForThread);
  IteratorType itWrite(image, outputRegionForThread);
  itScan.SetDirection(m_Dimension);
  itWrite.SetDirection(m_Dimension);

  // Max size of the queue
  int que_max_size = m_Radius * 2 + 1;
  int line_length = outputRegionForThread.GetSize(m_Dimension);
  int kernel_width = 2 * m_Radius + 1;

  // Allocate an array of the length of the line
  std::vector<OutputImagePixelType> line(line_length);

  // Start iterating over lines
  itScan.GoToBegin(); itWrite.GoToBegin();
  while(!itWrite.IsAtEnd())
    {
    int i; 
    // Compute the initial sum
    OutputImagePixelType sum = itk::NumericTraits<OutputImagePixelType>::Zero;
    for(i = 0; i < m_Radius; i++)
      {
      line[i] = itScan.Get();
      sum += line[i];
      ++itScan;
      }

    // For the next Radius + 1 values, add to the sum and write
    for(; i < kernel_width; i++)
      {
      line[i] = itScan.Get();
      sum += line[i];
      itWrite.Set(sum);
      ++itScan; ++itWrite;
      }

    // Continue until we hit the end of the scanline
    for(; i < line_length; i++)
      {
      line[i] = itScan.Get();
      sum += line[i] - line[i - kernel_width];
      itWrite.Set(sum);
      ++itScan; ++itWrite;
      }

    // Fill out the last bit
    for(; i < line_length + m_Radius; i++)
      {
      sum -= line[i - kernel_width];
      itWrite.Set(sum);
      ++itWrite;
      }

    assert(itScan.IsAtEndOfLine());
    assert(itWrite.IsAtEndOfLine());

    itScan.NextLine();
    itWrite.NextLine();
    }
}


template <class TInputPixel, class TOutputPixel>
class ImageToStatsFunctor
{
public:
  TOutputPixel operator() (const TInputPixel &a, const TInputPixel &b)
    {
    TOutputPixel c;
    c[0] = a;
    c[1] = b;
    c[2] = a*b;
    c[3] = a*a;
    c[4] = b*b;
    c[5] = 1.0;
    return c;
    }

  bool operator != (const ImageToStatsFunctor<TInputPixel, TOutputPixel> &other)
    { return false; }
};

template <class TInputPixel, class TOutputPixel>
class StatsToNCCFunctor
{
public:
  StatsToNCCFunctor(int n = 0) { this->m_N = n; }

  TOutputPixel operator() (const TInputPixel &c)
    {
    double a = c[0], b = c[1], ab = c[2], a2 = c[3], b2 = c[4], n = c[5];
    double Sab = (ab - a * b / n); 
    double Saa = (a2 - a * a / n); 
    double Sbb = (b2 - b * b / n); 
    return Sab / sqrt(Saa * Sbb);
    }

  bool operator != (const StatsToNCCFunctor<TInputPixel, TOutputPixel> &other)
    { return m_N != other.m_N; }

private:
  int m_N;
};

template <class TPixel, unsigned int VDim>
class MyVector : public itk::Vector<TPixel, VDim>
{
public:
  typedef MyVector<TPixel, VDim> Self;
  Self &operator = (int value) { 
    for(int i = 0; i < VDim; i++)
      (*this)[i] = 0;
    return *this;
  }
};

template <class TPixel, unsigned int VDim>
void
NormalizedCrossCorrelation<TPixel, VDim>
::operator() (itk::Size<VDim> radius)
{
  // Get two images from stack
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size()-1];
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size()-2];

  // Iterate over both images with a neighborhood iterator
  typedef typename  itk::ConstNeighborhoodIterator<ImageType> Iter;
  Iter q1(radius, i1,i1->GetBufferedRegion());
  Iter q2(radius, i2,i2->GetBufferedRegion());

  itk::TimeProbe clock1;
  clock1.Start();
 


  typename ImageType::Pointer iout = ImageType::New();
  iout->CopyInformation(i1);
  iout->SetRegions(i1->GetBufferedRegion());
  iout->Allocate();

  Iterator it(iout, iout->GetBufferedRegion());

  for(; !q1.IsAtEnd() && !q2.IsAtEnd(); ++q1, ++q2, ++it)
    {
    // Compute correlation at this location
    double sumx=0.0, ssqx = 0.0; 
    double sumy=0.0, ssqy = 0.0; 
    double sumxy=0.0;
    size_t n = 0;

    for(size_t i = 0; i < q1.Size(); i++)
      {
      bool inb1,inb2;
      double x = q1.GetPixel(i, inb1);
      double y = q2.GetPixel(i, inb2);
      sumx += x;
      sumy += y;
      sumxy += x * y;
      ssqx += x * x;
      ssqy += y * y;
      n++;
      }

    double cc = (sumxy - sumx * sumy / n) / 
      (sqrt((ssqx - sumx * sumx / n) * (ssqy - sumy * sumy / n)));

    it.Set(cc);
    }

  clock1.Stop();
  std::cout << "Total slow code: " << clock1.GetTotal() << std::endl; 

  itk::TimeProbe clock2;
  clock2.Start();

  // Alternative approach
  typedef itk::Vector<TPixel, 6> StatsVector;
  typedef itk::Image<StatsVector, VDim> StatsImage;
  typedef ImageToStatsFunctor<TPixel,StatsVector> StatsFunctor;
  typedef itk::BinaryFunctorImageFilter<ImageType, ImageType, StatsImage, StatsFunctor> StatsFilter;
  typename StatsFilter::Pointer fltStats = StatsFilter::New();
  fltStats->SetInput1(i1);
  fltStats->SetInput2(i2);
  fltStats->Update();

  // Mean filter - use separable
  typedef OneDimensionalInPlaceAccumulateFilter<StatsImage> MeanFilter1D;

  // Create a chain of these separable filters
  typename itk::ImageSource<StatsImage>::Pointer pipeTail = fltStats.GetPointer();
  for(int d = 0; d < VDim; d++)
    {
    typename MeanFilter1D::Pointer flt = MeanFilter1D::New();
    flt->SetInput(pipeTail->GetOutput());
    flt->SetDimension(d);
    flt->SetRadius(radius[d]);
    pipeTail = flt.GetPointer();

    flt->Update();
    /*
    typedef itk::ImageFileWriter<StatsImage> W;
    typename W::Pointer w = W::New();
    w->SetInput(flt->GetOutput());
    char buffer[256];
    sprintf(buffer, "dump%02d.nii.gz", d);
    w->SetFileName(buffer);
    w->Update();
    */
    }
    
  typedef StatsToNCCFunctor<StatsVector, TPixel> NCCFunctor;
  typedef itk::UnaryFunctorImageFilter<StatsImage, ImageType, NCCFunctor> NCCFilter;
  typename NCCFilter::Pointer fltNCC = NCCFilter::New();

  int n = 1;
  for(int i = 0; i < VDim; i++)
   n *= 2 * radius[i] + 1;

  NCCFunctor funk(n);
  fltNCC->SetInput(pipeTail->GetOutput());
  fltNCC->SetFunctor(funk);
  fltNCC->Update();


  clock2.Stop();
  std::cout << "Total fast code: " << clock2.GetTotal() << std::endl; 


  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(iout);
  c->m_ImageStack.push_back(fltNCC->GetOutput());
}

// Invocations
template class NormalizedCrossCorrelation<double, 2>;
template class NormalizedCrossCorrelation<double, 3>;
template class NormalizedCrossCorrelation<double, 4>;
