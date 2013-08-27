#include "BinaryMathOperation.h"

#include "itkAddImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMinimumImageFilter.h"
#include "itkMaximumImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BinaryMathOperation<TPixel, VDim>
::operator() (Operation op)
{
  // Check input availability
  if(c->m_ImageStack.size() < 2)
    {
    cerr << "Binary operations require two images on the stack" << endl;
    throw -1;
    }

  // Get the last two images
  ImagePointer i1 = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer i2 = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Write something
  *c->verbose << "Adding #" << c->m_ImageStack.size() - 1 << " and "  
    << c->m_ImageStack.size() - 2 << endl;

  // Select the operation
  typedef typename itk::ImageToImageFilter<ImageType, ImageType> BaseFilterType;
  typename BaseFilterType::Pointer filter;

  switch(op) 
    {
    case ADD:
      filter = itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case SUBTRACT:
      filter = itk::SubtractImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case MULTIPLY:
      filter = itk::MultiplyImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case DIVIDE:
      filter = itk::DivideImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case MINIMUM:
      filter = itk::MinimumImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    case MAXIMUM:
      filter = itk::MaximumImageFilter<ImageType, ImageType, ImageType>::New();
      break;
    }

  // Run the filter
  filter->SetInput(0, i1);
  filter->SetInput(1, i2);
  filter->Update();

  // Replace the images with the product
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(filter->GetOutput());
}

// Invocations
template class BinaryMathOperation<double, 2>;
template class BinaryMathOperation<double, 3>;
