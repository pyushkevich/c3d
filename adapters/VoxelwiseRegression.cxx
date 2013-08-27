#include "VoxelwiseRegression.h"
#include "vnl/vnl_file_matrix.h"
#include "vnl/vnl_rank.h"
#include "vnl/algo/vnl_matrix_inverse.h"

template <class TPixel, unsigned int VDim>
void
VoxelwiseRegression<TPixel, VDim>
::operator() (size_t order)
{
  // Two images are needed from the stack
  ImagePointer X = c->m_ImageStack[c->m_ImageStack.size() - 1];
  ImagePointer Y = c->m_ImageStack[c->m_ImageStack.size() - 2];

  // Build the design matrix and observation matrix
  size_t n = X->GetBufferedRegion().GetNumberOfPixels();
  vnl_matrix<double> design(n, order), observ(n, 1);
  TPixel *px = X->GetBufferPointer(); TPixel *py = Y->GetBufferPointer();
  for(size_t i = 0; i < n; i++, px++, py++)
    {
    TPixel vx = *px, vxj = 1.0;
    for(size_t j = 0; j < order; j++)
      {
      design(i,j) = vxj;
      vxj *= vx;
      }
    observ(i, 0) = *py;
    }

  // Execute GLM on the two matrices
  size_t rank = vnl_rank(design, vnl_rank_row);

  // Compute A
  vnl_matrix<double> A = 
    vnl_matrix_inverse<double>(design.transpose() * design).pinverse(rank);

  // Compute bhat = Ax * Y
  vnl_matrix<double> bhat = (A * design.transpose()) * observ;

  // We should report R^2, etc.
  for(size_t j = 0; j < order; j++)
    cout << "REGCOEFF[" << j << "] = " << bhat(j,0) << endl;
}

// Invocations
template class VoxelwiseRegression<double, 2>;
template class VoxelwiseRegression<double, 3>;
