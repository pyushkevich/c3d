#include "ConnectedComponents.h"
#include "ThresholdImage.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkCastImageFilter.h"

template <class TPixel, unsigned int VDim>
void
ConnectedComponents<TPixel, VDim>
::operator() ()
{
  // The image is assumed to be binary. If background is non-zero, call binarize
  // to map the background to zero
  if(c->m_Background != 0.0)
    {
    ThresholdImage<TPixel, VDim> thresh(c);
    thresh(c->m_Background, c->m_Background, 0.0, 1.0);
    }

  // Get the last image on the stack
  ImagePointer image = c->m_ImageStack.back();

  // Integer image typedef
  typedef itk::OrientedRASImage<int, VDim> IntImageType;
  
  // Construct the connected components filter
  typedef itk::ConnectedComponentImageFilter<ImageType, IntImageType> CCFilter;

  // Relabel the components
  typedef itk::RelabelComponentImageFilter<IntImageType, IntImageType> RCFilter;

  // Describe what we are doing
  *c->verbose << "Computing connected components of #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Calling ConnectedComponentImageFilter" << endl;

  // Plug in the filter's components
  typename CCFilter::Pointer fltConnect = CCFilter::New();
  fltConnect->SetInput(image);
  fltConnect->SetFullyConnected(false);
  fltConnect->Update();

  // Describe what we are doing
  *c->verbose << "  Calling RelabelComponentImageFilter" << endl;

  // Relabel and order components
  typename RCFilter::Pointer fltRelabel = RCFilter::New();
  fltRelabel->SetInput(fltConnect->GetOutput());
  fltRelabel->Update();

  // Print the statistics about the connected components
  long szpx = fltRelabel->GetSizeOfObjectInPixels(1);
  cout << "  There are " << 
    fltRelabel->GetNumberOfObjects() << " connected components." << endl;
  cout << "  Largest component has " << szpx << " pixels." << endl;

  // We have to convert the image back to the native type
  typedef itk::CastImageFilter<IntImageType, ImageType> CastFilter;
  typename CastFilter::Pointer fltCast = CastFilter::New();
  fltCast->SetInput(fltRelabel->GetOutput());
  fltCast->Update();

  // Store the output
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(fltCast->GetOutput());
}

// Invocations
template class ConnectedComponents<double, 2>;
template class ConnectedComponents<double, 3>;
