#include "FlipImage.h"
#include "itkFlipImageFilter.h"

template <class TPixel, unsigned int VDim>
void
FlipImage<TPixel, VDim>
::operator() (string axis)
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Create a flip filter
  typedef itk::FlipImageFilter<ImageType> FlipType;
  typename FlipType::Pointer flipper = FlipType::New();
  typename FlipType::FlipAxesArrayType flipax;

  // Set up the axes to flip
  for(size_t i = 0; i < VDim; i++)
    if(axis.find('x'+i) != string::npos || axis.find('X'+i) != string::npos)
      flipax[i] = true;
    else
      flipax[i] = false;

  // Say what we are doing
  *c->verbose << "Flipping #" << c->m_ImageStack.size()-1 << " about " << flipax << endl;

  // Do some processing ...
  flipper->SetInput(img);
  flipper->SetFlipAxes(flipax);
  flipper->Update();

  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(flipper->GetOutput());
}

// Invocations
template class FlipImage<double, 2>;
template class FlipImage<double, 3>;
