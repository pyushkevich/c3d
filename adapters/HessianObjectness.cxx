#include "HessianObjectness.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkHessianToObjectnessMeasureImageFilter.h"

template <class TPixel, unsigned int VDim>
void
HessianObjectness<TPixel, VDim>
::operator() (int dimension, double minscale, double maxscale)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Define the filter
  typedef typename itk::NumericTraits<TPixel>::RealType RealPixelType;
  typedef itk::SymmetricSecondRankTensor<RealPixelType, VDim> HessianPixelType;
  typedef itk::Image<HessianPixelType, VDim> HessianImageType;

  typedef itk::HessianToObjectnessMeasureImageFilter<HessianImageType, ImageType> 
    ObjectnessFilterType;

  typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType, ImageType>
    MultiScaleEnhancementFilterType;

  // Create the objectness filter
  typename ObjectnessFilterType::Pointer of = ObjectnessFilterType::New();
  of->SetScaleObjectnessMeasure(false); // why?
  of->SetBrightObject(dimension > 0);
  of->SetObjectDimension(abs(dimension));
  of->SetAlpha(0.5);
  of->SetBeta(0.5);
  of->SetGamma(5.0);

  // Create the main filter
  typename MultiScaleEnhancementFilterType::Pointer filter = MultiScaleEnhancementFilterType::New();
  filter->SetInput(img);
  filter->SetHessianToMeasureFilter(of);
  filter->SetSigmaStepMethodToLogarithmic();
  filter->SetSigmaMaximum(maxscale);
  filter->SetSigmaMinimum(minscale);
  filter->SetNumberOfSigmaSteps(minscale == maxscale ? 1 : 10);
  
  // Run the filter
  *c->verbose << "Extracting Hessian Objectness from #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Object dimension: " << of->GetObjectDimension() << endl; 
  *c->verbose << "  Object type: " << (of->GetBrightObject() ? "bright" : "dark") << endl; 
  *c->verbose << "  Sigma range: " << filter->GetSigmaMinimum() << " " << filter->GetSigmaMaximum() << endl;

  filter->Update();

  // Do some processing ...
  ImagePointer result = filter->GetOutput();
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(result);
}

// Invocations
template class HessianObjectness<double, 2>;
template class HessianObjectness<double, 3>;
