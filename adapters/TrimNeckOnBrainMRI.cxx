/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    TrimNeckOnBrainMRI.cxx
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

#include "TrimNeckOnBrainMRI.h"

#include "SwapDimensions.h"
#include "SmoothImage.h"
#include "ResampleImage.h"
#include "StructureTensorEigenValues.h"
#include "ScaleShiftImage.h"
#include "LandmarksToSpheres.h"
#include "RFTrain.h"
#include "RFApply.h"
#include "ThresholdImage.h"
#include "ResliceImage.h"
#include "TrimImage.h"
#include "LevelSetSegmentation.h"
#include "MathematicalMorphology.h"
#include "ExtractRegion.h"
#include "CreateInterpolator.h"

template <class TPixel>
void
TrimNeckOnBrainMRI<TPixel>::operator()(const TrimNeckOnBrainMRIParameters & param)
{
  if (c->m_ImageStack.size() == 0)
    throw ConvertException("At least one image is required on the stack for neck trim command");

  // Read the original image
  ImagePointer input = c->m_ImageStack.back();

  // Explain what we are doing
  *c->verbose << "Starting Brain MRI neck trimming algorithm on #" << c->m_ImageStack.size() << "." << std::endl;
  *c->verbose << "  Head Height:   " << param.HeadHeight << "mm" << std::endl;
  *c->verbose << "  Top Clearance: " << param.ClearanceHeight << "mm" << std::endl;

  // Get a copy of the current interpolator to replace later
  typename Converter::Interpolator::Pointer interp_orig = c->GetInterpolator();
  CreateInterpolator<TPixel, VDim>(this->c).CreateLinear();

  // Swap dimensions so that last dimension is inferior to superior
  SwapDimensions<TPixel, VDim>(this->c)(std::vector<std::string>({ "RAS" }));

  // Smooth image by one voxel
  auto stdev = c->ReadRealSize("1vox");
  SmoothImage<TPixel, VDim>(this->c)(stdev, true);

  // Resample to 2mm isotropic
  ResampleImage<TPixel, VDim> resampler(this->c);
  RealVector                  vox_resam = c->ReadRealSize("2mm");
  SizeType                    sz_resam = resampler.ComputeSizeFromTargetSpacing(vox_resam);
  resampler(sz_resam);

  // The resampled image will be needed again
  ImagePointer t1_resam = c->m_ImageStack.back();

  // Duplicate the intensity image - so its features are used by the classifier
  c->m_ImageStack.push_back(t1_resam);

  // Compute structured tensor eigenvalues - stack now contains a bunch of features
  StructureTensorEigenValues<TPixel, VDim>(this->c)(2.0, 4);

  // Create an image of zeros
  c->m_ImageStack.push_back(t1_resam);
  ScaleShiftImage<TPixel, VDim>(this->c)(0.0, 0.0);

  // Define the landmark spheres
  std::vector<std::pair<std::string, int>> lm_spec = { { "60x40x40%", 1 }, { "40x60x40%", 1 }, { "40x40x60%", 1 },
                                                       { "60x60x40%", 1 }, { "60x40x60%", 1 }, { "40x60x60%", 1 },
                                                       { "60x60x60%", 1 }, { "3x3x3%", 2 },    { "97x3x3%", 2 },
                                                       { "3x97x3%", 2 },   { "3x3x97%", 2 },   { "97x97x3%", 2 },
                                                       { "97x3x97%", 2 },  { "3x97x97%", 2 },  { "97x97x97%", 2 } };

  // Place landmark spheres. In the middle of the volume, we place foreground label
  // spheres and in the corners of the volume, background spheres. This should allow
  // us to differentiate head tissue from background
  typename LandmarksToSpheres<TPixel, VDim>::LandmarkList lms;
  for (auto & spec : lm_spec)
  {
    auto vec = c->ReadRealVector(spec.first.c_str(), true);
    typename LandmarksToSpheres<TPixel, VDim>::Landmark lm;
    for(unsigned int i = 0; i < 3; i++)
      lm.first[i] = vec[i];
    lm.second = spec.second;
    lms.push_back(lm);
  }
  LandmarksToSpheres<TPixel, VDim>(this->c).RasterizeLandmarks(lms, 15.0);

  // Retain the samples
  ImagePointer samples = c->m_ImageStack.back();

  // Train random forest classifier
  RFTrain<TPixel, VDim> rftrain(this->c);
  auto                  rfparam = *c->GetRandomForestParameters();
  rfparam.patch_radius.Fill(1);
  typename RFTrain<TPixel, VDim>::RFClassifierType::Pointer classifier = RFTrain<TPixel, VDim>::RFClassifierType::New();
  rftrain.TrainClassifier(rfparam, classifier);

  // Pop the segmentation image from the stack
  c->m_ImageStack.pop_back();

  // Apply the trained classifier
  RFApply<TPixel, VDim> rfapply(this->c);
  rfapply.ApplyClassifier(classifier);

  // Take the last image - this is the background probability map, then clear the generated images
  ImagePointer bg_prob = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.pop_back();

  // Smooth and downsample the probability image
  c->m_ImageStack.push_back(bg_prob);
  auto stdev_prob = c->ReadRealSize("1vox");
  SmoothImage<TPixel, VDim>(this->c)(stdev_prob, true);
  auto sz_res_prob = c->ReadSizeVector("50%");
  ResampleImage<TPixel, VDim>(this->c)(sz_res_prob);

  // Map range [0 1] to [-1 1] for level set, and duplicate on the stack
  ScaleShiftImage<TPixel, VDim>(this->c)(2.0, -1.0);
  c->m_ImageStack.push_back(c->m_ImageStack.back());

  // Push the samples, threshold them, and reslice to the probability image
  c->m_ImageStack.push_back(samples);
  ThresholdImage<TPixel, VDim>(this->c)(1, 1, 1, -1);
  ResliceImage<TPixel, VDim>(this->c)("identity", "");

  // Run the level set
  auto lsparam = *c->GetLevelSetParameters();
  lsparam.CurvatureWeight = 5.0;
  LevelSetSegmentation<TPixel, VDim>(this->c)(300, lsparam);

  // Resample the level set to the size of the probability image
  ImagePointer levelset = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(bg_prob);
  c->m_ImageStack.push_back(levelset);
  ResliceImage<TPixel, VDim>(this->c)("identity", "");

  // Threshold to create a binary mask of the head from the levelset
  ThresholdImage<TPixel, VDim>(this->c)(0, vnl_huge_val(0.0), 1, 0);
  ImagePointer ls_mask = c->m_ImageStack.back();

  // Compute the dimensions of the region that will be trimmed from the mask
  double z_region = (param.HeadHeight + param.ClearanceHeight) / 2.0;

  // Dilate the mask by full size in X,Y but not in Z
  MathematicalMorphology<TPixel, VDim> dilate(this->c);
  SizeType                             sz_ls = ls_mask->GetBufferedRegion().GetSize();
  SizeType                             sz_dilate = { sz_ls[0], sz_ls[1], 0 };
  dilate(MathematicalMorphology<TPixel, VDim>::DILATION, 1, sz_dilate);

  // Trim the mask by the clearance height along z axis
  RealVector sz_trim = { 0, 0, param.ClearanceHeight };
  TrimImage<TPixel, VDim>(this->c)(sz_trim, TrimImage<TPixel, VDim>::SPECIFY_MARGIN);

  // Extract region corresponding to head length, finally threshold
  RegionType bbox = c->m_ImageStack.back()->GetBufferedRegion();
  auto       sz_bbox = bbox.GetSize();
  double     z_trim_mm = param.HeadHeight + param.ClearanceHeight;
  bbox.SetSize(2, static_cast<unsigned>(0.5 + z_trim_mm / c->m_ImageStack.back()->GetSpacing()[2]));
  ExtractRegion<TPixel, VDim>(this->c)(bbox);
  ThresholdImage<TPixel, VDim>(this->c)(-vnl_huge_val(0.0), vnl_huge_val(0.0), 1, 0);

  // Resample the mask back to the original image size using nearest neighbor
  ImagePointer slab_mask = c->m_ImageStack.back();
  c->m_ImageStack.pop_back();
  c->m_ImageStack.push_back(input);
  c->m_ImageStack.push_back(slab_mask);
  CreateInterpolator<TPixel, VDim>(this->c).CreateNN();
  ResliceImage<TPixel, VDim>(this->c)("identity", "");

  *c->verbose << "Completed Brain MRI neck trimming algorithm." << std::endl;

  // Restore the original interpolator
  c->SetInterpolator(interp_orig);
}

// Invocations
template class TrimNeckOnBrainMRI<double>;
