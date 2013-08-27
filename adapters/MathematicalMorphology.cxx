#include "MathematicalMorphology.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

template <class TPixel, unsigned int VDim>
void
MathematicalMorphology<TPixel, VDim>
::operator() (bool erode, TPixel value, SizeType radius)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  *c->verbose << "Applying " << (erode ? "erosion" : "dilation") << " to #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Foreground value :    " << value << endl;
  *c->verbose << "  Ball radius      :    " << radius << endl;                                                              

  // Define the structuring element
  typedef itk::BinaryBallStructuringElement<TPixel, VDim> Element;
  Element elt;
  elt.SetRadius(radius);
  elt.CreateStructuringElement();

  // Chose the right filter
  typedef itk::BinaryMorphologyImageFilter<ImageType, ImageType, Element> FilterType;
  ImagePointer output;
  if(erode)
    {
    typedef itk::BinaryErodeImageFilter<ImageType, ImageType, Element> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(img);
    filter->SetErodeValue(value);
    filter->SetKernel(elt);
    filter->Update();
    output = filter->GetOutput();
    }
  else
    {
    typedef itk::BinaryDilateImageFilter<ImageType, ImageType, Element> FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput(img);
    filter->SetDilateValue(value);
    filter->SetKernel(elt);
    filter->Update();
    output = filter->GetOutput();
    }
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(output);
}

// Invocations
template class MathematicalMorphology<double, 2>;
template class MathematicalMorphology<double, 3>;
