#include "SampleImage.h"

template <class TPixel, unsigned int VDim>
void
SampleImage<TPixel, VDim>
::operator() (const RealVector &x)
{
  // Create a point
  itk::Point<double, VDim> p;
  for(size_t i = 0; i < VDim; i++)
    p[i] = x[i];

  // Map the point to a continuous index
  ImagePointer image = c->m_ImageStack.back();
  itk::ContinuousIndex<double, VDim> cidx;
  image->TransformRASPhysicalPointToContinuousIndex(p, cidx);

  // Describe what we are doing
  *c->verbose << "Probing #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Physical (RAS) Coordinates: " << p << endl;
  *c->verbose << "  Voxel Coordinates         : " << cidx << endl;

  // Use the interpolator in c
  c->GetInterpolator()->SetInputImage(image);
  m_Result = c->GetInterpolator()->EvaluateAtContinuousIndex(cidx);
  *c->verbose << "  Using " << c->m_Interpolation << " interpolation" << endl;

  // Print out the interpolated value
  cout << "Interpolated image value at " << x << " is " << m_Result << endl;
}

// Invocations
template class SampleImage<double, 2>;
template class SampleImage<double, 3>;
