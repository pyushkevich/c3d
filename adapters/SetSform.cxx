/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    SetSform.cxx
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

#include "SetSform.h"
#include <string>
#include <iostream>
#include <itkAffineTransform.h>
#include <itkTransformFileReader.h>
#include <itkTransformFactory.h>
#include "itkResampleImageFilter.h"
#include "gsGSAffine3DTransform.h"

#include <itkTransformFactory.h>


template<unsigned int VDim>
void MyReadMatrix(const char *fname, itk::Matrix<double,VDim+1,VDim+1> &mat)
  {
  ifstream fin(fname);
  for(size_t i = 0; i < VDim+1; i++)
    for(size_t j = 0; j < VDim+1; j++)
      if(fin.good())
        {
        fin >> mat[i][j];
        }
      else
        {
        throw ConvertException("Unable to read matrix %s", fname);
        }
  fin.close();
  }


template <class TPixel, unsigned int VDim>
void
SetSform<TPixel, VDim>
::operator() (string fn_tran)
{

  // Check input availability
  if(c->m_ImageStack.size() < 1)
    {
    cerr << "No image to set the sform" << endl;
    throw -1;
    }

  // Get image from stack
  ImagePointer img = c->m_ImageStack.back();

  // Get the transform
  itk::Matrix<double,VDim+1,VDim+1> matrix;
  MyReadMatrix<VDim>(fn_tran.c_str(), matrix);
  vnl_matrix<double> mat(VDim+1,VDim+1,0.0);
  mat.update( matrix.GetVnlMatrix());

  // Set the matrix
  img->SetVoxelSpaceToRASPhysicalSpaceMatrix( mat );
  
  // Put result on stack
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(img);
}

// Invocations
template class SetSform<double, 2>;
template class SetSform<double, 3>;
template class SetSform<double, 4>;
