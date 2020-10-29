/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    AlignByLandmarks.cxx
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

#include "AlignByLandmarks.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkShapeLabelObject.h"
#include "itkLabelMap.h"
#include "itkShapeLabelMapFilter.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_svd.h>
#include <iostream>
#include <fstream>

template <class TPixel, unsigned int VDim>
typename AlignByLandmarks<TPixel, VDim>::CentroidMap
AlignByLandmarks<TPixel, VDim>
::ExtractCentroids(ImageType *imLabel)
{
  typedef unsigned long LabelType;
  typedef itk::ShapeLabelObject<LabelType, VDim> LabelObjectType;
  typedef itk::LabelMap<LabelObjectType> LabelMapType;

  typedef itk::LabelImageToLabelMapFilter<ImageType, LabelMapType> ConverterType;
  typename ConverterType::Pointer fltConvert = ConverterType::New();
  fltConvert->SetInput(imLabel);
  fltConvert->SetBackgroundValue(0);

  typedef itk::ShapeLabelMapFilter<LabelMapType> ShapeFilterType;
  typename ShapeFilterType::Pointer fltShape = ShapeFilterType::New();
  fltShape->SetInput(fltConvert->GetOutput());

  fltShape->Update();

  CentroidMap cm;

  LabelMapType *labelMap = fltConvert->GetOutput();
  typename LabelMapType::LabelObjectVectorType vecLabel = labelMap->GetLabelObjects();
  for(unsigned int i = 0; i < vecLabel.size(); i++)
    {
    const LabelObjectType *labelObject = vecLabel[i];
    itk::Point<double, VDim> ctr = labelObject->GetCentroid();

    // Map the point to RAS coordinate space
    itk::ContinuousIndex<double, VDim> cidx;
    imLabel->TransformPhysicalPointToContinuousIndex(ctr, cidx);
    imLabel->TransformContinuousIndexToRASPhysicalPoint(cidx, ctr);

    cm[labelObject->GetLabel()] = ctr;
    }

  return cm;
}

template <class TPixel, unsigned int VDim>
void
AlignByLandmarks<TPixel, VDim>
::operator() (int dof, std::string fn_output)
{
  // Get image from stack
  if(c->m_ImageStack.size() < 2)
    throw ConvertException("Too few images on the stack for landmark alignment");

  ImagePointer im_moving = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();
  ImagePointer im_fixed = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();

  // Check DOF
  if(dof != 6 && dof != 7 && dof != 12)
    throw ConvertException("Landmark alignment degrees of freedom parameter is wrong");

  *c->verbose << "Computing landmark alignment between images" << std::endl;
  *c->verbose << "  Degrees of freedom: " << dof << std::endl;
  *c->verbose << "  Output file: " << fn_output << std::endl;

  // Extract the centroids of the landmarks
  CentroidMap cm_fixed = ExtractCentroids(im_fixed);
  CentroidMap cm_moving = ExtractCentroids(im_moving);

  // Find common labels
  std::set<unsigned long> common;
  for(typename CentroidMap::iterator it = cm_fixed.begin(); it != cm_fixed.end(); ++it)
    {
    if(cm_moving.find(it->first) != cm_moving.end())
      common.insert(it->first);
    }

  // Construct matrices A and B containing the points
  vnl_matrix<double> A0(common.size(), VDim), B0(common.size(), VDim), A, B;
  vnl_vector<double> cA(VDim), cB(VDim); 
  cA.fill(0); cB.fill(0);
  int i = 0;
  double totaldist = 0;
  *c->verbose << "  Distances before alignment: " << std::endl;
  for(typename std::set<unsigned long>::iterator it = common.begin(); it!=common.end(); ++it)
    {
    unsigned long label = *it;
    for(int j = 0; j < VDim; j++)
      {
      A0(i,j) = cm_fixed[label][j];
      B0(i,j) = cm_moving[label][j];
      }

    double dist = (A0.get_row(i) - B0.get_row(i)).magnitude();
    totaldist += dist;
    *c->verbose << "    Label " << label << ":   " 
      << A0.get_row(i) << "       " << B0.get_row(i) << std::endl;
    *c->verbose << "    Label " << label << ":   " << dist << std::endl;
    cA += A0.get_row(i);
    cB += B0.get_row(i);
    i++;
    }
  *c->verbose << "    Total:    " << totaldist << std::endl;

  // We are solving for T(A) - B. Compute the translation and remove center from the data
  A = A0; B = B0;
  cA /= A.rows(); cB /= B.rows();
  vnl_vector<double> translation = cB - cA;
  for(i = 0; i < A0.rows(); i++)
    {
    A.set_row(i, A0.get_row(i) - cA);
    B.set_row(i, B0.get_row(i) - cB);
    }

  // The output transform components: affine component M and translation b
  vnl_matrix<double> out_M;
  vnl_vector<double> out_b;

  // Split by model
  if(dof == 6)
    {
    // Just solve for the rotation component
    vnl_svd<double> svd(A.transpose() * B);
    vnl_matrix<double> R = svd.U() * svd.V().transpose();

    // Get the transform components 
    out_M = R.transpose();
    }

  else if(dof == 7)
    {
    // Recover the uniform scaling of the data. To do this, compute the root mean 
    // squared distance in both datasets
    double rmsd_A = 0, rmsd_B = 0;
    for(i = 0; i < A.rows(); i++)
      {
      rmsd_A += (A.get_row(i)).squared_magnitude();
      rmsd_B += (B.get_row(i)).squared_magnitude();
      }
    rmsd_A = sqrt(rmsd_A / A.rows());
    rmsd_B = sqrt(rmsd_B / B.rows());

    // Remove the scale from the matrices
    A /= rmsd_A; B /= rmsd_B;

    // Solve for the rotation
    vnl_svd<double> svd(A.transpose() * B);
    vnl_matrix<double> R = svd.U() * svd.V().transpose();

    // Get the transform components 
    out_M = R.transpose() * (rmsd_B/rmsd_A);
    }

  else if(dof == 12)
    {
    // Here we just need to solve for the entire matrix M
    vnl_matrix<double> Q(VDim * VDim, VDim * VDim);
    vnl_matrix<double> Qb(VDim, VDim), P;
    Q.fill(0.0); Qb = A.transpose() * A;
    for(int i = 0; i < VDim; i++)
      Q.update(Qb, i * VDim, i * VDim);

    P = B.transpose() * A;
    vnl_vector<double> vP(P.data_block(), VDim * VDim);

    vnl_vector<double> vR = vnl_svd<double>(Q).solve(vP);
    out_M = vnl_matrix<double>(vR.data_block(), VDim, VDim);
    }

  // Put together a 4 x 4 transform
  
  // Put together into a 4x4 affine transform
  vnl_matrix<double> out_mat(VDim+1, VDim+1);
  out_mat.set_identity();
  out_mat.update(out_M, 0, 0);
  vnl_vector<double> tmat = out_mat.get_column(VDim);
  tmat.update(cB - out_M * cA, 0);
  out_mat.set_column(VDim, tmat);

  // Save it
  std::ofstream fout(fn_output.c_str());
  for(size_t i = 0; i < 4; i++)
    for(size_t j = 0; j < 4; j++)
      fout << out_mat[i][j] << (j < 3 ? " " : "\n");
  fout.close();

  // Compute the final distances
  *c->verbose << "  Distances after alignment: " << std::endl;
  totaldist = 0;
  i = 0;
  for(typename std::set<unsigned long>::iterator it = common.begin(); it!=common.end(); ++it)
    {
    unsigned long label = *it;
    vnl_vector<double> x(4), y;
    x.fill(1);
    x.update(A0.get_row(i), 0);
    y = out_mat * x;
    *c->verbose << "    Label " << label << ":   " 
      << y << "       " << B0.get_row(i) << std::endl;
    double dist = (y.extract(3) - B0.get_row(i)).magnitude();
    *c->verbose << "    Label " << label << ":   " << dist << std::endl;
    totaldist += dist;
    i++;
    }
  *c->verbose << "    Total:    " << totaldist << std::endl;
}

// Invocations
template class AlignByLandmarks<double, 2>;
template class AlignByLandmarks<double, 3>;
template class AlignByLandmarks<double, 4>;
