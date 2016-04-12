/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ComputeMoments.cxx
  Language:  C++
  Website:   itksnap.org/c3d
  Copyright (c) 2014 Paul A. Yushkevich
  
  This file is part of C3D, a command-line companion tool to ITK-SNAP

  C3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#include "ComputeMoments.h"
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_matrix.h>

template <class TPixel, unsigned int VDim>
void
ComputeMoments<TPixel, VDim>
::operator() ()
{
  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Count the number of non-zero voxels
  int n = 0;
  for(Iterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(it.Value() != c->m_Background)
      n++;
    }

  // Allocate the matrix for the PCA
  typedef vnl_matrix<double> Mat;
  typedef vnl_vector<double> Vec;
  Mat x(n, VDim);
  Vec x_row_sum(VDim, 0.0), x_mean(VDim, 0.0), x_row(VDim, 0.0);

  // Place the NIFTI coordinates of all points into the matrix
  int k = 0;
  for(Iterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
    {
    if(it.Value() != c->m_Background)
      {
      itk::Point<double, VDim> pos;
      img->TransformIndexToRASPhysicalPoint(it.GetIndex(), pos);
      for(int i = 0; i < VDim; i++)
        x(k, i) = pos[i];

      x_row_sum += x.get_row(k);

      k++;
      }
    }

  // Compute the mean of the points
  x_mean = x_row_sum / n;

  // Subtract the mean from all the points
  for(int k = 0; k < n; k++)
    x.set_row(k, x.get_row(k) - x_mean);

  // Compute A`A
  Mat xTx = x.transpose() * x;

  // Perform the eigen-analysis
  vnl_symmetric_eigensystem<double> eig(xTx);

  // Report the eigen-vectors and eigen-values
  std::cout << "Centroid: " << x_mean << std::endl;
  for(int j = 0; j < VDim; j++)
    {
    std::cout << "Mode " << j << " eigenvalue: " << eig.get_eigenvalue(j) << std::endl;
    std::cout << "Mode " << j << " eigenvector: " << eig.get_eigenvector(j) << std::endl;
    }

  // Do some processing ...
  // ImagePointer result = ...;
  
  // Put result on stack
  // c->m_ImageStack.pop_back();
  // c->m_ImageStack.push_back(result);
}

// Invocations
template class ComputeMoments<double, 2>;
template class ComputeMoments<double, 3>;
template class ComputeMoments<double, 4>;
