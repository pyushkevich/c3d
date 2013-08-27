#include "UnaryMathOperation.h"

template <class TPixel, unsigned int VDim>
void
UnaryMathOperation<TPixel, VDim>
::operator() (double (*func) (double))
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Print debug info
  *c->verbose << "Applying unary math operation to #" << c->m_ImageStack.size() << endl;

  // Apply the appropriate math operation (in place)
  Iterator it(img, img->GetBufferedRegion());
  for(; !it.IsAtEnd(); ++it)
    it.Set(func(it.Get()));
}

// Invocations
template class UnaryMathOperation<double, 2>;
template class UnaryMathOperation<double, 3>;
