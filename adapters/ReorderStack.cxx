#include "ReorderStack.h"

template <class TPixel, unsigned int VDim>
void
ReorderStack<TPixel, VDim>
::operator() (size_t k)
{
  // Get the number of images on the stack
  size_t n = c->m_ImageStack.size();

  // Make sure n is divisible by k
  if(n % k != 0)
    throw ConvertException("Can not reorder %d images using stride of %d;"
      " %d is not divisible by %d.", n, k, n, k);

  // Explain what we are doing
  *c->verbose << "Reordering " << n << " images with stride of " << k << endl;

  // Make a copy of the stack
  vector<ImagePointer> temp_stack = c->m_ImageStack;
  c->m_ImageStack.clear();

  // Traverse the stack
  for(size_t i = 0; i < k; i++)
    for(size_t j = i; j < n; j += k)
      c->m_ImageStack.push_back(temp_stack[j]);
}

// Invocations
template class ReorderStack<double, 2>;
template class ReorderStack<double, 3>;
