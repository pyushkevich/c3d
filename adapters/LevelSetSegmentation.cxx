#include "LevelSetSegmentation.h"
#include "itkSegmentationLevelSetImageFilter.h"
#include "itkSegmentationLevelSetFunction.h"
#include "itkShiftScaleImageFilter.h"

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

protected:

  MyLevelSetFunction() {};
  ~MyLevelSetFunction() {};

  virtual void PrintSelf(std::ostream &os, itk::Indent indent) const
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

  virtual void PrintSelf(std::ostream &os, itk::Indent indent) const
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
::operator() (int nIter)
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
  fnSegment->SetCurvatureWeight(c->m_LevSetCurvature);
  fnSegment->SetAdvectionWeight(c->m_LevSetAdvection);
  fnSegment->SetPropagationWeight(1.0);
  fnSegment->Initialize(rad);
  fnSegment->SetSpeedImage(i2);

  // Set the inputs to the segmentation filter
  fltSegment->SetSegmentationFunction(fnSegment);
  fltSegment->SetInput(i1);
  fltSegment->SetFeatureImage(i2);
  fltSegment->SetNumberOfLayers(3);
  fltSegment->SetIsoSurfaceValue(0.0);
  fltSegment->SetNumberOfIterations(nIter);

  *c->verbose << "  NIterations:    " << nIter << endl;
  *c->verbose << "  Curv Weight:    " << c->m_LevSetCurvature << endl;
  *c->verbose << "  Adv Weight:     " << c->m_LevSetAdvection << endl;

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
