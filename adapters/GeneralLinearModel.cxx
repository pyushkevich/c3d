#include "GeneralLinearModel.h"
#include "vnl/vnl_file_matrix.h"
#include "vnl/vnl_rank.h"
#include "vnl/algo/vnl_matrix_inverse.h"

template <class TPixel, unsigned int VDim>
void
GeneralLinearModel<TPixel, VDim>
::operator() (string fn_matrix, string fn_contrast)
{
  // Read the matrix from file
  vnl_file_matrix<double> mat(fn_matrix.c_str());
  if(!mat)
    throw string("Unable to read matrix from file given");

  // Read the contrast from file
  vnl_file_matrix<double> con(fn_contrast.c_str());
  if(!con)
    throw string("Unable to read contrast from file given");

  // Check that the number of images matches
  if(c->m_ImageStack.size() != mat.rows())
    throw string("Matrix number of rows does not match stack size");

  // Check that the columns in matrix match contrast vector
  if(con.columns() != mat.columns())
    throw string("Matrix and contrast vector must have same number of columns");

  *c->verbose << "Running GLM on " << mat.rows() << " images" << endl;
  *c->verbose << "  design matrix: " << mat << endl;
  *c->verbose << "  contrast vector: " << con << endl;

  // Some matrices
  size_t rank = vnl_rank(mat, vnl_rank_row);
  vnl_matrix<double> A = 
    vnl_matrix_inverse<double>(mat.transpose() * mat).pinverse(rank);

  // Load all images into a Y matrix (can we do this)
  size_t n = c->m_ImageStack.front()->GetBufferedRegion().GetNumberOfPixels();
  vnl_matrix<double> Y(mat.rows(), n);
  for(size_t i = 0; i < mat.rows(); i++)
    {
    TPixel *pix = c->m_ImageStack[i]->GetBufferPointer();
    for(size_t j = 0; j < n; j++)
      { Y(i,j) = pix[j]; }
    }

  // Compute bhat = Ax * Y
  vnl_matrix<double> bhat = (A * mat.transpose()) * Y;

  // Compute the contrast
  vnl_matrix<double> res = con * bhat;

  // Populate the first image
  ImagePointer iout = c->m_ImageStack.front();
  for(size_t i = 0; i < n; i++)
    iout->GetBufferPointer()[i] = res(0,i);

  c->m_ImageStack.clear();
  c->m_ImageStack.push_back(iout);
}

// Invocations
template class GeneralLinearModel<double, 2>;
template class GeneralLinearModel<double, 3>;
