/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    LevelSetSegmentation.cxx
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

#include "LevelSetSegmentation.h"
#include "itkSegmentationLevelSetImageFilter.h"
#include "itkSegmentationLevelSetFunction.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageAlgorithm.h"

template<class TImageType, class TFeatureImageType = TImageType>
class MyLevelSetFunction :
  public itk::SegmentationLevelSetFunction<TImageType, TFeatureImageType>
{
public:
  // All the regular filter stuff
  typedef MyLevelSetFunction Self;
  typedef itk::SegmentationLevelSetFunction<
    TImageType, TFeatureImageType> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  // The input/output image types
  typedef TImageType InputImageType;

  // New pointer
  itkNewMacro(Self);

  // Calculate speed image - just copy the feature image
  virtual void CalculateSpeedImage() ITK_OVERRIDE
    {
    itk::ImageAlgorithm::Copy( this->GetFeatureImage(),
      this->GetSpeedImage(),
      this->GetFeatureImage()->GetRequestedRegion(),
      this->GetFeatureImage()->GetRequestedRegion() );
    }

protected:

  MyLevelSetFunction() {};
  ~MyLevelSetFunction() {};

  virtual void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE
    { os << indent << "MyLevelSetFunction"; }

private:

};


template<class TInputImage,
         class TFeatureImage,
         class TOuputPixelType = double>
class MyLevelSetFilter :
  public itk::SegmentationLevelSetImageFilter<
    TInputImage, TFeatureImage, TOuputPixelType>
{
public:
  // All the regular filter stuff
  typedef MyLevelSetFilter Self;
  typedef itk::SegmentationLevelSetImageFilter<
    TInputImage, TFeatureImage, TOuputPixelType> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  // Some typedefs of image types
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename Superclass::FeatureImageType FeatureImageType;

  // Input pixel values
  typedef typename InputImageType::PixelType InputPixelType;

  // Standard ITK Macros
  itkTypeMacro(MyLevelSetFilter, itk::SegmentationLevelSetFilter);
  itkNewMacro(Self);

protected:
  MyLevelSetFilter() {};
  ~MyLevelSetFilter() {};

  virtual void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE
    { os << indent << "MyLevelSetFilter"; }

private:

  // Unimplemented copy stuff
  MyLevelSetFilter(const Self &s);
  void operator=(const Self &s);
};

void DumpProgress(itk::Object *object, const itk::EventObject &obj, void *client_data)
{
  itk::ProcessObject *po = (itk::ProcessObject *)object;
  cout << po->GetProgress() << endl;
}


template <class TPixel, unsigned int VDim>
void
LevelSetSegmentation<TPixel, VDim>
::operator() (int nIter, const LevelSetParameters &param)
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Level set segmentation requires two images on the stack!" << endl;
    throw -1;
    }

  // Get the last two images
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Report what the filter is doing
  *c->verbose << "Running level set segmentation (";
  *c->verbose << "#" << c->m_ImageStack.size() - 1 << " is speed, ";
  *c->verbose << "#" << c->m_ImageStack.size() << " is init)" << endl;

  // Create a segmentation filter
  typedef MyLevelSetFilter<UnorientedImageType, UnorientedImageType, TPixel> SegFilter;
  typename SegFilter::Pointer fltSegment = SegFilter::New();

  // Set up the radius
  itk::Size<VDim> rad; rad.Fill(1);

  // Create the function
  typedef MyLevelSetFunction<UnorientedImageType> SegFunction;
  typename SegFunction::Pointer fnSegment = SegFunction::New();
  fnSegment->SetCurvatureWeight(param.CurvatureWeight);
  fnSegment->SetAdvectionWeight(param.AdvectionWeight);
  fnSegment->SetPropagationWeight(1.0);
  fnSegment->Initialize(rad);
  fnSegment->SetSpeedImage(i2);

  // Set the inputs to the segmentation filter
  fltSegment->SetSegmentationFunction(fnSegment);
  fltSegment->SetInput(i1);
  fltSegment->SetFeatureImage(i2);
  fltSegment->SetNumberOfLayers(3);
  fltSegment->SetIsoSurfaceValue(0.0);
  fltSegment->SetMaximumRMSError(1.0e-4);
  fltSegment->SetNumberOfIterations(nIter);

  *c->verbose << "  NIterations:    " << nIter << endl;
  *c->verbose << "  Curv Weight:    " << param.CurvatureWeight << endl;
  *c->verbose << "  Adv Weight:     " << param.AdvectionWeight << endl;

  // Execute the filter
  fltSegment->Update();

  // Finally, map to an oriented image
  typedef itk::ShiftScaleImageFilter<UnorientedImageType, ImageType> DummyFilter;
  typename DummyFilter::Pointer fltDummy = DummyFilter::New();
  fltDummy->SetInput(fltSegment->GetOutput());
  fltDummy->SetScale(1.0);
  fltDummy->SetShift(0.0);
  fltDummy->Update();

  // Take the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltDummy->GetOutput());

  *c->verbose << "Level set done after" << fltSegment->GetElapsedIterations() << " iterations" << endl;
}

// Invocations
template class LevelSetSegmentation<double, 2>;
template class LevelSetSegmentation<double, 3>;
template class LevelSetSegmentation<double, 4>;
