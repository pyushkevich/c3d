/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ExportPatches.cxx
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

#include "ExportPatches.h"
#include "itkNeighborhoodIterator.h"
#include "itkComposeImageFilter.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "FastLinearInterpolator.h"
#include <vnl/vnl_random.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>

template <unsigned int VDim>
class RandomMatrixGenerator
{
public:
  typedef vnl_matrix_fixed<double, VDim+1, VDim+1> MatrixType;
  static void Generate(double sigma_radians, vnl_random &rand, MatrixType &R)
    {
    throw ConvertException("Random rotation not implemented in 4D");
    }
};

template <>
class RandomMatrixGenerator<3>
{
public:
  typedef vnl_matrix_fixed<double, 4, 4> MatrixType;
  static void Generate(double sigma_radians, vnl_random &rand, MatrixType &R)
    {
    // Get a uniformly distributed rotation axis
    double x = rand.normal(), y = rand.normal(), z = rand.normal();
    double r = sqrt(x*x + y*y + z*z);

    // Generate a random rotation angle
    double theta = rand.normal() * sigma_radians;
    double scale = sin(theta / 2.0) / r;

    // Define the unit quaternion
    double w = cos(theta / 2.0);
    x *= scale, y *= scale, z *= scale; 

    // Generate a rotation matrix
    R(0,0) = 1.0 - 2 * y * y - 2 * z * z;
    R(1,1) = 1.0 - 2 * x * x - 2 * z * z;
    R(2,2) = 1.0 - 2 * x * x - 2 * y * y;
    R(0,1) = 2 * x * y - 2 * w * z;
    R(1,0) = 2 * x * y + 2 * w * z;
    R(2,0) = 2 * x * z - 2 * w * y;
    R(0,2) = 2 * x * z + 2 * w * y;
    R(1,2) = 2 * y * z - 2 * w * x;
    R(2,1) = 2 * y * z + 2 * w * x;
    std::cout << theta * 180 / vnl_math::pi << std::endl;
    std::cout << w << "," << x << "," << y << "," << z << std::endl;
    std::cout << R << std::endl;
    }
};


template <>
class RandomMatrixGenerator<2>
{
public:
  typedef vnl_matrix_fixed<double, 3, 3> MatrixType;
  static void Generate(double sigma_radians, vnl_random &rand, MatrixType &R)
    {
    // Get a uniformly distributed rotation axis
    double theta = rand.normal() * sigma_radians;
    double cos_theta = cos(theta), sin_theta = sin(theta);
    R(0,0) = cos_theta;
    R(1,1) = cos_theta;
    R(0,1) = -sin_theta;
    R(1,0) = sin_theta;
    }
};

template <class TPixel, unsigned int VDim> 
int ExportPatches<TPixel, VDim>::m_NumberOfAugmentations = 0;

template <class TPixel, unsigned int VDim> 
double ExportPatches<TPixel, VDim>::m_AugmentationRotationSigma = 0.0;

template <class TPixel, unsigned int VDim>
void
ExportPatches<TPixel, VDim>
::SetAugmentationParameters(int n_aug, double rot_sigma)
{
  m_NumberOfAugmentations = n_aug;
  m_AugmentationRotationSigma = rot_sigma;
}

template <class TPixel, unsigned int VDim>
void
ExportPatches<TPixel, VDim>
::operator() (const char *out_file, const SizeType &radius, double sample_frequency)
{
  // The last image is the mask - used to decide which voxels to sample
  ImagePointer mask = c->PopImage();

  // Number of channels to export
  int n_chan = c->GetStackSize();

  // Add the remaining images to the composite filter
  typedef itk::VectorImage<TPixel, VDim> VectorImageType;
  typedef itk::ComposeImageFilter<ImageType, VectorImageType> CompositeFilter;
  typename CompositeFilter::Pointer composer = CompositeFilter::New();
  for(int j = 0; j < n_chan; j++)
    composer->SetInput(j, c->PeekImage(j));
  composer->Update();
  typename VectorImageType::Pointer vecimg = composer->GetOutput();

  // Fast linear interpolator for the image
  typedef FastLinearInterpolator<VectorImageType, float, VDim> FastInterp;
  FastInterp flint(vecimg);

  // Create an image for sampling (float)
  typedef itk::VectorImage<float, VDim> FloatVectorImageType;
  typename FloatVectorImageType::Pointer sample = FloatVectorImageType::New();
  typename FloatVectorImageType::RegionType region;
  for(int i = 0; i < VDim; i++)
    {
    region.SetIndex(i, -radius[i]);
    region.SetSize(i, 2 * radius[i] + 1);
    }
  sample->SetNumberOfComponentsPerPixel(n_chan);
  sample->SetRegions(region);
  sample->Allocate();

  // Create the sample vector
  float *sample_vec = sample->GetBufferPointer();

  // Open a file for writing
  FILE *f = fopen(out_file, "wb");

  // Are we doing augmentation?
  double aug_sigma_angle_radians = m_AugmentationRotationSigma * vnl_math::pi / 180;

  // Rotation matrix for augmentation
  typedef vnl_matrix_fixed<double, VDim+1, VDim+1> Mat44;
  Mat44 R;
  Mat44 D = mask->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
  Mat44 Dinv = vnl_inverse<double>(D);

  R.set_identity();

  // Keep track of the number of patches written
  int n_patches = 0;

  // Iterate over the mask
  vnl_random rand;
  for(Iterator iter(mask, mask->GetBufferedRegion()); !iter.IsAtEnd(); ++iter)
    {
    double drand = rand.drand32(sample_frequency);
    if(iter.Value() != c->m_Background)
      {
      if(drand <= 1.0)
        {
        // This voxel is getting sampled
        for(int i = 0; i < m_NumberOfAugmentations+1; i++)
          {
          // This matrix holds the affine transform applied during sampling
          vnl_matrix_fixed<double,VDim+1,VDim+1> RV;

          // First augmentation is identity - sampling the actual image
          if(i == 0)
            {
            RV.set_identity();
            }
          else
            {
            RandomMatrixGenerator<VDim>::Generate(aug_sigma_angle_radians, rand, R);
            RV = Dinv * (R * D);
            }

          // Get the index of the current position
          IndexType idx_center = iter.GetIndex();

          // Get a pointer to the sampling region
          float *sample_buffer = sample->GetBufferPointer();
          memset(sample_buffer, 0, sizeof(float) * sample->GetPixelContainer()->Size());
          float *p = sample_buffer;

          // Perform iteration over the sample
          typedef itk::ImageLinearConstIteratorWithIndex<FloatVectorImageType> LineIter;
          for(LineIter sliter(sample, region); !sliter.IsAtEnd(); sliter.NextLine())
            {
            // Apply transformation for the beginning of the line
            // -- this is the offset of the first pixel in the line
            typename LineIter::IndexType idx = sliter.GetIndex();

            // Compute the starting index and delta for moving along X
            float cix[VDim], cix_delta[VDim];
            for(int a = 0; a < VDim; a++)
              {
              cix[a] = idx_center[a];
              for(int b = 0; b < VDim; b++)
                cix[a] += RV(a,b) * idx[b];
              cix_delta[a] = RV(a,0);
              }

            // We are ready to actually interpolate
            for(; !sliter.IsAtEndOfLine(); ++sliter)
              {
              flint.Interpolate(cix, p);
              for(int b = 0; b < VDim; b++)
                cix[b] += cix_delta[b];
              p += n_chan;
              } 
            } 

          // Save the current sample
          fwrite(sample_buffer, sizeof(float), sample->GetPixelContainer()->Size(), f);
          n_patches++;

          } // augmentation
        } // voxel selected for sampling 
      } // voxel inside mask
    } // iteration over mask image

  // Close the output
  fclose(f);

  // Report the number of patches
  std::cout << "Exported patches:" << n_patches << std::endl;

  // Restore the mask
  c->PushImage(mask);
}

// Invocations
template class ExportPatches<double, 2>;
template class ExportPatches<double, 3>;
template class ExportPatches<double, 4>;
