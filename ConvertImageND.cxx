/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ConvertImageND.cxx
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

#include "ConvertImageND.h"

#include "AddImages.h"
#include "AlignByLandmarks.h"
#include "AntiAliasImage.h"
#include "ApplyMetric.h"
#include "ApplyNoise.h"
#include "BiasFieldCorrectionN4.h"
#include "BinaryHoleFill.h"
#include "BinaryImageCentroid.h"
#include "BinaryMathOperation.h"
#include "Canny.h"
#include "ClipImageIntensity.h"
#include "ComputeFFT.h"
#include "ComputeMoments.h"
#include "ComputeOverlaps.h"
#include "ConnectedComponents.h"
#include "ConvertAdapter.h"
#include "Convolution.h"
#include "CoordinateMap.h"
#include "CopyTransform.h"
#include "CreateImage.h"
#include "CreateInterpolator.h"
#include "DicomSeriesList.h"
#include "ExportPatches.h"
#include "ExtractRegion.h"
#include "ExtractSlice.h"
#include "ExtrudeSegmentation.h"
#include "FlipImage.h"
#include "FillBackgroundWithNeighborhoodNoise.h"
#include "GeneralLinearModel.h"
#include "HessianEigenValues.h"
#include "HessianObjectness.h"
#include "HistogramMatch.h"
#include "ImageERF.h"
#include "ImageGradient.h"
#include "ImageLaplacian.h"
#include "LabelOverlapMeasures.h"
#include "LabelStatistics.h"
#include "LandmarksToSpheres.h"
#include "LaplacianSharpening.h"
#include "LevelSetSegmentation.h"
#include "MatchBoundingBoxes.h"
#include "MathematicalMorphology.h"
#include "MeanFilter.h"
#include "MedianFilter.h"
#include "MixtureModel.h"
#include "MomentsFeatures.h"
#include "MRFVote.h"
#include "MultiplyImages.h"
#include "NormalizedCrossCorrelation.h"
#include "NormalizeLocalWindow.h"
#include "OtsuThreshold.h"
#include "OverlayLabelImage.h"
#include "PadImage.h"
#include "PeronaMalik.h"
#include "PrintImageInfo.h"
#include "Rank.h"
#include "ReadImage.h"
#include "ReciprocalImage.h"
#include "ReorderStack.h"
#include "ReplaceIntensities.h"
#include "ResampleImage.h"
#include "ResliceImage.h"
#include "RetainLabels.h"
#include "RFApply.h"
#include "RFTrain.h"
#include "SampleImage.h"
#include "ScalarToRGB.h"
#include "ScaleShiftImage.h"
#include "SetOrientation.h"
#include "SetSform.h"
#include "SignedDistanceTransform.h"
#include "SLICSuperVoxel.h"
#include "SmoothImage.h"
#include "SplitMultilabelImage.h"
#include "StapleAlgorithm.h"
#include "StructureTensorEigenValues.h"
#include "SwapDimensions.h"
#include "TestImage.h"
#include "ThresholdImage.h"
#include "TileImages.h"
#include "TrimImage.h"
#include "UnaryMathOperation.h"
#include "UpdateMetadataKey.h"
#include "Vote.h"
#include "VoxelwiseComponentFunction.h"
#include "VoxelwiseRegression.h"
#include "WarpImage.h"
#include "WarpLabelImage.h"
#include "WriteImage.h"
#include "WeightedSum.h"
#include "WeightedSumVoxelwise.h"
#include "WrapDimension.h"

#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

// Support for regular expressions via KWSYS in ITK
#include <itksys/RegularExpression.hxx>

// Image to image filter - for some global tolerance code
#include <itkImageToImageFilter.h>

// Documentation manual
#include "Documentation.h"

// Markdown documentation string generated at compile-time
// this looks a little odd, but works - the include file contains raw bytes
unsigned char c3d_md[] = {
  #include "markdown_docs.h"
  0x00 
};

using namespace itksys;

extern const char *ImageConverter_VERSION_STRING;

// Helper function: read a double, throw exception if unreadable
double myatof(char *str)
{
  char *end = 0;
  double d = strtod(str, &end);
  if (*end != 0)
    throw "strtod conversion failed";
  return d;
};

long myatol(char *str)
{
  char *end = 0;
  double d = strtol(str, &end, 10);
  if (*end != 0)
    throw "strtol conversion failed";
  return d;
};


std::string str_to_lower(const char *input)
{
  std::string s(input);
  std::transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
  return s;
}

// Check whether value is a valid float (no leading spaces allowed)
bool is_double(const char *input)
{
  std::istringstream iss(input);
  double f;
  iss >> std::noskipws >> f; // noskipws considers leading whitespace invalid

  // Check the entire string was consumed and if either failbit or badbit is set
  return iss.eof() && !iss.fail();
}



    /*
    out << "Command Listing: " << endl;
    out << "    -accum" << endl;
    out << "    -add" << endl;
    out << "    -align-landmarks, -alm" << endl;
    out << "    -anisotropic-diffusion, -ad" << endl;
    out << "    -antialias, -alias" << endl;
    out << "    -as, -set" << endl;
    out << "    -atan2" << endl;
    out << "    -background" << endl;
    out << "    -biascorr" << endl;
    out << "    -binarize" << endl;
     ***  out << "    -canny" << endl;
     out << "    -centroid" << endl;
     out << "    -clear" << endl;
     out << "    -clip" << endl;
     out << "    -colormap, -color-map" << endl;
     out << "    -connected-components, -connected, -comp" << endl;
     out << "    -coordinate-map-voxel, -cmv" << endl;
     out << "    -coordinate-map-physical, -cmp" << endl;
     out << "    -copy-transform, -ct" << endl;
     out << "    -cos" << endl;
     out << "    -create" << endl;
     out << "    -dilate" << endl;
     out << "    -divide" << endl;
     out << "    -dup" << endl;
     out << "    -endaccum" << endl;
     out << "    -endfor" << endl;
     out << "    -erode" << endl;
     out << "    -erf" << endl;
     out << "    -exp" << endl;
     ***  out << "    -fft" << endl;
     out << "    -flip" << endl;
     out << "    -foreach" << endl;
     out << "    -glm" << endl;
     ***  out << "    -hessobj, -hessian-objectness" << endl;
     ***  out << "    -histmatch, -histogram-match" << endl;
     out << "    -holefill, -hf" << endl;
     out << "    -info" << endl;
     out << "    -info-full" << endl;
     out << "    -insert, -ins" << endl;
     out << "    -interpolation, -interp, -int" << endl;
     out << "    -iterations" << endl;
     ***  out << "    -label-overlap" << std::endl;
     out << "    -label-statistics, -lstat" << endl;
     out << "    -landmarks-to-spheres, -lts" << endl;
     out << "    -laplacian, -laplace" << endl;
     out << "    -levelset" << endl;
     out << "    -levelset-curvature" << endl;
     out << "    -levelset-advection" << endl;
     out << "    -ln, -log" << endl;
     out << "    -log10" << endl;
     out << "    -max, -maximum" << endl;
     out << "    -mcs, -multicomponent-split" << endl;
     out << "    -mean" << endl;
     out << "    -merge" << endl;
     ***  out << "    -mf, -mean-filter" << endl;
     out << "    -mi, -mutual-info" << endl;
     out << "    -min, -minimum" << endl;
     out << "    -mixture, -mixture-model" << endl;
     out << "    -multiply, -times" << endl;
     out << "    -n4, -n4-bias-correction" << endl;
     out << "    -ncc, -normalized-cross-correlation" << endl;
     out << "    -nmi, -normalized-mutual-info" << endl;
     ***  out << "    -nlw, -normwin, -normalize-local-window" << endl;
     ***  out << "    -normpdf" << endl;
     out << "    -noround" << endl;
     out << "    -nospm" << endl;
     out << "    -o" << endl;
     out << "    -omc, -output-multicomponent" << endl;
     out << "    -oo, -output-multiple" << endl;
    out << "    -orient" << endl;
    out << "    -origin" << endl;
    out << "    -origin-voxel" << endl;
    out << "    -overlap" << endl;
    out << "    -overlay-label-image, -oli" << endl;
    out << "    -pad" << endl;
    out << "    -percent-intensity-mode, -pim" << endl;
    ***  out << "    -pixel" << endl;
    out << "    -pop" << endl;
    out << "    -popas" << endl;
    out << "    -probe" << endl;
    out << "    -push, -get" << endl;
    out << "    -rank" << endl;
    out << "    -reciprocal" << endl;
    out << "    -region" << endl;
    out << "    -reorder" << endl;
    out << "    -replace" << endl;
    out << "    -resample" << endl;
    out << "    -resample-mm" << endl;
    out << "    -reslice-itk" << endl;
    out << "    -reslice-matrix" << endl;
    out << "    -reslice-identity" << endl;
    ***  out << "    -rf-train" << endl;
    ***  out << "    -rf-apply" << endl;
    out << "    -rms" << endl;
    out << "    -round" << endl;
    out << "    -scale" << endl;
    out << "    -set-sform" << endl;
    out << "    -shift" << endl;
    out << "    -signed-distance-transform, -sdt" << endl;
    out << "    -sin" << endl;
    out << "    -slice" << endl;
    out << "    -smooth" << endl;
    out << "    -spacing" << endl;
    out << "    -split" << endl;
    out << "    -sqrt" << endl;
    out << "    -staple" << endl;
    out << "    -spm" << endl;
    out << "    -stretch" << endl;
    ***  out << "    -subtract" << endl;
    ***  out << "    -supervoxel, -sv" << endl;
    out << "    -test-image" << endl;
    out << "    -test-probe" << endl;
    out << "    -threshold, -thresh" << endl;
    out << "    -tile" << endl;
    out << "    -trim" << endl;
    out << "    -trim-to-size" << endl;
    out << "    -type" << endl;
    out << "    -verbose" << endl;
    ***  out << "    -version" << endl;
    out << "    -vote" << endl;
    out << "    -vote-label" << endl;
    out << "    -voxel-sum" << endl;
    out << "    -voxel-integral, -voxel-int" << endl;
    out << "    -voxelwise-regression, -voxreg" << endl;
    out << "    -warp" << endl;
    out << "    -wrap" << endl;
    out << "    -weighted-sum, -wsum" << endl;
    out << "    -weighted-sum-voxelwise, -wsv" << endl;
    */

/**
 * Parameters for the various algorithms. Stored in a separate structure 
 * in order to reduce number of variables declared in the header
 */
template <class TPixel, unsigned int VDim>
struct ConvertAlgorithmParameters
{
  // Root mean square error for anti-aliasing algorithm
  double m_AntiAliasRMS;

  // Level set algorithm parameters
  LevelSetParameters m_LevelSet;

  // Random forest parameters
  RFParameters<TPixel, VDim> m_RandomForest;

  ConvertAlgorithmParameters()
    {
    m_AntiAliasRMS = 0.07;
    }
};


template<class TPixel, unsigned int VDim>
ImageConverter<TPixel,VDim>
::ImageConverter()
  : verbose(&devnull)
{
  // Initialize to defaults
  m_TypeId = "float";
  m_Background = 0.0;
  m_RoundFactor = 0.5;
  m_FlagSPM = false;
  m_UseCompression = false;
  m_MultiComponentSplit = false;
  m_Iterations = 0;
  m_LoopType = LOOP_NONE;
  m_PercentIntensityMode = PIM_QUANTILE;

  // Create the parameters
  m_Param = new ParameterType();

  // Documentation initially NULL to not waste time parsing it
  m_Documentation = NULL;

  // Create an interpolator
  m_Interpolation = "linear";
  CreateInterpolator<TPixel, VDim>(this).CreateLinear();

  // Set orientation and coordinate tolerances to something more
  // reasonable than the ITK defaults
  itk::ImageToImageFilterCommon::SetGlobalDefaultCoordinateTolerance(1.0e-4);
  itk::ImageToImageFilterCommon::SetGlobalDefaultDirectionTolerance(1.0e-4);
}


template<class TPixel, unsigned int VDim>
ImageConverter<TPixel,VDim>
::~ImageConverter()
{
  delete m_Param;
  if(m_Documentation)
    delete m_Documentation;
}




template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PrintCommandListing(std::ostream &out)
{
  if(!m_Documentation)
    m_Documentation = new Documentation(c3d_md);

  // Print the automatically generated command listing
  m_Documentation->PrintCommandListing(out);

  // Print additional information on getting help
  out << "Getting help:" << std::endl;

  out << "    " 
    << std::setw(32) << std::left 
    << "-h"
    << ": List commands" << std::endl;

  out << "    " 
    << std::setw(32) << std::left 
    << "-h command"
    << ": Print help on command (e.g. -h add)" << std::endl;

  out << "    " 
    << std::setw(32) << std::left 
    << "-manual"
    << ": Print complete reference manual" << std::endl;
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PrintCommandHelp(std::ostream &out, const char *command)
{
  if(!m_Documentation)
    m_Documentation = new Documentation(c3d_md);

  if(!m_Documentation->PrintCommandHelp(out, command))
    {
    out << "No help available for command " << command << std::endl;
    }
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PrintManual(std::ostream &out)
{
  if(!m_Documentation)
    m_Documentation = new Documentation(c3d_md);

  m_Documentation->PrintManual(out);
}

template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::ProcessCommand(int argc, char *argv[])
{
  // Get the first command
  string cmd = argv[0];

  // cout << "COMMAND: " << cmd << endl;

  // Commands in alphabetical order
  if (cmd == "-accum")
    {
    if (this->m_LoopType != LOOP_NONE)
      throw ConvertException("Nested -foreach and -accum loops are not allowed");

    this->m_LoopType = LOOP_ACCUM;
    return this->AccumulateLoop(argc, argv);
    }

  else if (cmd == "-acos")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_acos);
    return 0;
    }

  else if (cmd == "-add")
    {
    BinaryMathOperation<TPixel, VDim> adapter(this);
    adapter(BinaryMathOperation<TPixel, VDim>::ADD);
    return 0;
    }

  else if (cmd == "-align-landmarks" || cmd == "-alm")
    {
    int dof = atoi(argv[1]);
    std::string fnout = argv[2];
    AlignByLandmarks<TPixel,VDim> adapter(this);
    adapter(dof, fnout);
    return 2;
    }

  else if (cmd == "-anisotropic-diffusion" || cmd == "-ad")
    {
    double cond = atof(argv[1]);
    int niter = atoi(argv[2]);
    PeronaMalik<TPixel, VDim> adapter(this);
    adapter(cond, (size_t) niter);
    return 2;
    }

  // Anti-alias a binary image, turning it into a smoother floating point image;
  // the argument is the iso-surface value
  // This command is affected by -iterations and -rms flags
  else if (cmd == "-antialias" || cmd == "-alias")
    {
    AntiAliasImage<TPixel, VDim> adapter(this);
    adapter(atof(argv[1]), m_Param->m_AntiAliasRMS);
    return 1;
    }

  // Associate variable name with image currently at the top of
  // the stack
  else if (cmd == "-as" || cmd == "-set")
    {
    string var(argv[1]);
    if (m_ImageStack.size() == 0)
      throw ConvertException("No image to assign to variable %s", var.c_str());
    m_ImageVars[var] = m_ImageStack.back();
    return 1;
    }

  else if (cmd == "-asin")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_asin);
    return 0;
    }

  else if (cmd == "-atan2")
    {
    BinaryMathOperation<TPixel, VDim> adapter(this);
    adapter(BinaryMathOperation<TPixel, VDim>::ATAN2);
    return 0;
    }

  else if (cmd == "-background")
    {
    m_Background = atof(argv[1]);
    *verbose << "Background value set to " << m_Background << endl;
    return 1;
    }

  else if (cmd == "-biascorr" || cmd == "-n4" || cmd == "-n4-bias-correction" )
    {
    BiasFieldCorrectionN4<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  // f(x) = (x == xBackground) ? 0 : 1
  else if (cmd == "-binarize")
    {
    ThresholdImage<TPixel, VDim> adapter(this);
    adapter(m_Background, m_Background, 0.0, 1.0);
    return 0;
    }

  else if (cmd == "-canny")
    {
    Canny<TPixel, VDim> adapter(this);
    RealVector sigma = ReadRealSize(argv[1]);
    double tLower = ReadIntensityValue(argv[1]);
    double tUpper = ReadIntensityValue(argv[3]);

    adapter(sigma, tLower, tUpper);
    return 3;
    }

  else if (cmd == "-ceil")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_ceil);
    return 0;
    }


  else if (cmd == "-centroid")
    {
    BinaryImageCentroid<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-centroid-mark")
    {
    double mark_value = atof(argv[1]);
    BinaryImageCentroid<TPixel, VDim> adapter(this);
    adapter(mark_value);
    return 1;
    }

  else if (cmd == "-connected-components" || cmd == "-connected" || cmd == "-comp")
    {
    ConnectedComponents<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-clear")
    {
    *verbose << "Clearing the stack" << endl;
    m_ImageStack.clear();
    return 0;
    }

  else if (cmd == "-clip")
    {
    double iMin = ReadIntensityValue(argv[1]);
    double iMax = ReadIntensityValue(argv[2]);
    ClipImageIntensity<TPixel, VDim>(this)(iMin, iMax);
    return 2;
    }

  else if (cmd == "-colormap" || cmd == "-color-map")
    {
    // Read the name of the color map
    std::string cmname = argv[1];

    // Optionally, read the range of the color map (instead of using image min/max)
    int np = 1;
    double int_min = 0.0, int_max = 0.0;

    // Try readin
    if (argc > 3)
      {
      try 
        {
        int_min = ReadIntensityValue(argv[2]);
        int_max = ReadIntensityValue(argv[3]);
        np = 3;
        }
      catch(ConvertException &exc) {}
      }

    ScalarToRGB<TPixel, VDim>(this)(cmname, int_min, int_max);
    return np;
    }

  else if (cmd == "-compress")
    {
    m_UseCompression = true;
    return 0;
    }

  else if (cmd == "-no-compress")
    {
    m_UseCompression = false;
    return 0;
    }

  else if (cmd == "-conv")
    {
    Convolution<TPixel,VDim>(this)();
    return 0;
    }

  else if (cmd == "-coordinate-map-voxel" || cmd == "-cmv")
    {
    CoordinateMap<TPixel,VDim>(this)(false);
    return 0;
    }

  else if (cmd == "-coordinate-map-physical" || cmd == "-cmp")
    {
    CoordinateMap<TPixel,VDim>(this)(true);
    return 0;
    }

  else if (cmd == "-copy-transform" || cmd == "-ct")
    {
    CopyTransform<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-cos")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_cos);
    return 0;
    }

  // Create a new image with given size and voxel size
  else if (cmd == "-create")
    {
    SizeType dims = ReadSizeVector(argv[1]);
    RealVector voxel = ReadRealSize(argv[2]);
    CreateImage<TPixel, VDim> adapter(this);
    adapter(dims, voxel);
    return 2;
    }

  else if (cmd == "-dicom-series-list")
    {
    DicomSeriesList<TPixel, VDim> adapter(this);
    adapter(argv[1]);
    return 1;
    }

  else if (cmd == "-dicom-series-read")
    {
    typedef ReadImage<TPixel, VDim> Adapter;
    typename Adapter::ImageInfo info;

    info.dicom_series_id = argv[2];

    Adapter adapter(this);
    adapter(argv[1], info);

    return 2;
    }

  else if (cmd == "-dilate")
    {
    MathematicalMorphology<TPixel,VDim> adapter(this);
    adapter(false, atof(argv[1]), ReadSizeVector(argv[2]));
    return 2;
    }

  else if (cmd == "-divide")
    {
    BinaryMathOperation<TPixel, VDim> adapter(this);
    adapter(BinaryMathOperation<TPixel, VDim>::DIVIDE);
    return 0;
    }

  else if (cmd == "-dup" || cmd == "-duplicate")
    {
    m_ImageStack.push_back(m_ImageStack.back());
    return 0;
    }

  else if (cmd == "-endaccum")
    {
    // This command ends the accum loop
    if (this->m_LoopType != LOOP_ACCUM)
      throw ConvertException("Out of place -endaccum command");
    this->m_LoopType = LOOP_NONE;
    return 0;
    }

  else if (cmd == "-endfor")
    {
    // This command ends the for loop
    if (this->m_LoopType != LOOP_FOREACH)
      throw ConvertException("Out of place -endfor command");
    this->m_LoopType = LOOP_NONE;
    return 0;
    }

  else if (cmd == "-erode")
    {
    MathematicalMorphology<TPixel,VDim> adapter(this);
    adapter(true, atof(argv[1]), ReadSizeVector(argv[2]));
    return 2;
    }

  else if (cmd == "-erf")
    {
    double thresh = atof(argv[1]);
    double scale = atof(argv[2]);
    ImageERF<TPixel, VDim> adapter(this);
    adapter(thresh, scale);
    return 2;
    }

  else if (cmd == "-exp")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_exp);
    return 0;
    }

  else if (cmd == "-export-patches" || cmd == "-xp")
    {
    ExportPatches<TPixel,VDim> adapter(this);
    std::string fn = argv[1];
    SizeType radius = this->ReadSizeVector(argv[2]);
    double freq = atof(argv[3]);
    adapter(fn.c_str(), radius, freq);
    return 3;
    }

  else if (cmd == "-export-patches-aug" || cmd == "-xpa")
    {
    int n_aug = atoi(argv[1]);
    double sigma_angle = atof(argv[2]);
    ExportPatches<TPixel,VDim>::SetAugmentationParameters(n_aug, sigma_angle);
    return 2;
    }

  else if (cmd == "-extrude-seg") 
    {
    ExtrudeSegmentation<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-fill-background-with-noise" || cmd == "-fbn")
    {
    FillBackgroundWithNeighborhoodNoise<TPixel, VDim> adapter(this);
    SizeType radius = ReadSizeVector(argv[1]);
    int steps = atoi(argv[2]);
    adapter(radius, steps);
    return 2;
    }

  else if (cmd == "-fft")
    {
    ComputeFFT<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-flip")
    {
    string flipax = argv[1];
    FlipImage<TPixel, VDim> adapter(this);
    adapter(flipax);
    return 1;
    }

  else if (cmd == "-floor")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_floor);
    return 0;
    }

  else if (cmd == "-foreach")
    {
    if (this->m_LoopType != LOOP_NONE)
      throw ConvertException("Nested loops are not allowed");
    this->m_LoopType = LOOP_FOREACH;
    return this->ForEachLoop(argc, argv);
    }

  else if (cmd == "-foreach-comp")
    {
    if (this->m_LoopType != LOOP_NONE)
      throw ConvertException("Nested loops are not allowed");
    this->m_LoopType = LOOP_FOREACH;
    int ncomp = atoi(argv[1]);
    return this->ForEachComponentLoop(ncomp, argc - 1, argv + 1) + 1;
    }

  else if (cmd == "-glm")
    {
    string mat(argv[1]);
    string con(argv[2]);
    GeneralLinearModel<TPixel, VDim> adapter(this);
    adapter(mat, con);
    return 2;
    }

  else if (cmd == "-grad" || cmd == "-gradient")
    {
    ImageGradient<TPixel,VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-h" || cmd == "-help" || cmd == "--help")
    {
    if(argc > 1 && argv[1][0] != '-')
      {
      PrintCommandHelp(std::cout, argv[1]);
      return 1;
      }
    else
      {
      PrintCommandListing(std::cout);
      return 0;
      }
    }

  else if(cmd == "-hf" || cmd == "-holefill")
    {
    double foreground = atof(argv[1]);
    bool full_conn = atoi(argv[2]);

    BinaryHoleFill<TPixel, VDim> adapter(this);
    adapter(foreground, full_conn);

    return 2;
    }

  else if (cmd == "-hesseig" || cmd == "-hessian-eigenvalues")
    {
    double scale = atof(argv[1]);
    HessianEigenValues<TPixel,VDim> adapter(this);
    adapter(scale);

    return 1;
    }

  else if (cmd == "-hessobj" || cmd == "-hessian-objectness")
    {
    int dimension = atoi(argv[1]);
    double minscale = atof(argv[2]);
    double maxscale = atof(argv[3]);
    
    HessianObjectness<TPixel, VDim> adapter(this);
    adapter(dimension, minscale, maxscale);

    return 3;
    }

  else if (cmd == "-histmatch" || cmd == "-histogram-match")
    {
    size_t nmatch = atoi(argv[1]);
    HistogramMatch<TPixel, VDim> adapter(this);
    adapter(nmatch);
    return 1;
    }

  else if (cmd == "-info")
    {
    PrintImageInfo<TPixel, VDim> adapter(this);
    adapter(false);
    return 0;
    }

  else if (cmd == "-info-full")
    {
    PrintImageInfo<TPixel, VDim> adapter(this);
    adapter(true);
    return 0;
    }

  else if (cmd == "-insert" || cmd == "-ins")
    {
    string var(argv[1]);
    size_t pos = (size_t) atoi(argv[2]);
    typename ImageVariableMap::iterator img = m_ImageVars.find(var);

    // Check if the variable exists
    if (img == m_ImageVars.end())
      throw ConvertException("No image assigned to variable %s", var.c_str());

    // Check if the position is ok
    if (m_ImageStack.size() > pos)
      throw ConvertException("Can not insert at position %i in stack of size %i", pos, m_ImageStack.size());

    // Insert at the appropriate place
    typename vector<ImagePointer>::iterator it = m_ImageStack.end();
    for(size_t i = 0; i < pos; i++) --it;
    m_ImageStack.insert(it, img->second);

    return 2;
    }

  else if (cmd == "-interpolation" || cmd == "-interp" || cmd == "-int")
    {
    // Adapter that creates interpolators
    CreateInterpolator<TPixel, VDim> adapter(this);

    // Interpret the interpolation type
    m_Interpolation = str_to_lower(argv[1]);

    if (m_Interpolation == "nearestneighbor" || m_Interpolation == "nn" || m_Interpolation == "0")
      {
      adapter.CreateNN();
      return 1;
      }
    else if (m_Interpolation == "linear" || m_Interpolation == "1")
      {
      adapter.CreateLinear();
      return 1;
      }
    else if (m_Interpolation == "cubic" || m_Interpolation == "3")
      {
      adapter.CreateCubic();
      return 1;
      }
    else if (m_Interpolation == "sinc")
      {
      adapter.CreateSinc();
      return 1;
      }
    else if (m_Interpolation == "gaussian")
      {
      RealVector sigma = ReadRealSize(argv[2]);
      adapter.CreateGaussian(sigma);
      return 2;
      }
    else if (m_Interpolation == "multilabel" || m_Interpolation == "ml")
      {
      RealVector sigma = ReadRealSize(argv[2]);
      adapter.CreateMultiLabel(sigma);
      return 2;
      }
    else
      {
      throw ConvertException("Unknown interpolation type: %s", m_Interpolation.c_str());
      }
    }

  else if (cmd == "-iterations")
    {
    m_Iterations = static_cast<size_t>(atoi(argv[1]));
    return 1;
    }

  else if (cmd == "-label-overlap")
    {
    LabelOverlapMeasures<TPixel, VDim>(this)();
    return 0;
    }

  else if (cmd == "-label-statistics" || cmd == "-lstat")
    {
    LabelStatistics<TPixel, VDim>(this)();
    return 0;
    }

  else if (cmd == "-landmarks-to-spheres" || cmd == "-lts")
    {
    char *fnland = argv[1];
    double radius = atof(argv[2]);
    LandmarksToSpheres<TPixel, VDim> (this)(fnland, radius);
    return 2;
    }

  else if (cmd == "-laplacian" || cmd == "-laplace")
    {
    ImageLaplacian<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-levelset")
    {
    int nIter = atoi(argv[1]);
    LevelSetSegmentation<TPixel, VDim> adapter(this);
    adapter(nIter, m_Param->m_LevelSet);
    return 1;
    }

  else if (cmd == "-levelset-curvature")
    {
    m_Param->m_LevelSet.CurvatureWeight = atof(argv[1]);
    return 1;
    }

  else if (cmd == "-levelset-advection")
    {
    m_Param->m_LevelSet.AdvectionWeight = atof(argv[1]);
    return 1;
    }

  else if (cmd == "-ln" || cmd == "-log")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_log);
    return 0;
    }

  else if (cmd == "-log10")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_log10);
    return 0;
    }

  // No else if here because of a windows compiler error (blocks nested too deeply)
  if (cmd == "-manual")
    {
    this->PrintManual(std::cout);
    return 0;
    }

  else if (cmd == "-match-bounding-box" || cmd == "-mbb")
    {
    MatchBoundingBoxes<TPixel,VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-maximum" || cmd == "-max")
    {
    BinaryMathOperation<TPixel, VDim> adapter(this);
    adapter(BinaryMathOperation<TPixel, VDim>::MAXIMUM);
    return 0;
    }

  else if (cmd == "-mcs" || cmd == "-multicomponent-split")
    {
    m_MultiComponentSplit = true; return 0;
    }

  else if (cmd == "-mean")
    {
    size_t n = m_ImageStack.size();
    for(size_t i = 1; i < n; i++)
      {
      AddImages<TPixel,VDim> adapter(this);
      adapter();
      }
    ScaleShiftImage<TPixel, VDim> scaler(this);
    scaler(1.0 / n, 0.0);
    return 0;
    }

  else if(cmd == "-median" || cmd == "-median-filter")
    {
    MedianFilter<TPixel, VDim> adapter(this);
    SizeType radius = this->ReadSizeVector(argv[1]);
    adapter(radius);
    return 1;
    }

  else if (cmd == "-merge")
    {
    Vote<TPixel, VDim> adapter(this);
    adapter(true);
    return 0;
    }

  else if (cmd == "-mf" || cmd == "-mean-filter")
    {
    MeanFilter<TPixel, VDim> adapter(this);
    SizeType sz = ReadSizeVector(argv[1]);
    adapter(sz);
    return 1;
    }

  else if (cmd == "-mi" || cmd == "-mutual-info")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fnf("none");
    string fnm("none");
    if (argc > 1)
      {
      fnm = argv[1];
      nret = 1;
      }
    if (argc == 3)
      {
      fnf = argv[2];
      nret = 2;
      }
    adapter("MI", fnf.c_str(), fnm.c_str());
    return nret;
    }

  else if (cmd == "-minimum" || cmd == "-min")
    {
    BinaryMathOperation<TPixel, VDim> adapter(this);
    adapter(BinaryMathOperation<TPixel, VDim>::MINIMUM);
    return 0;
    }

  else if (cmd == "-mixture" || cmd == "-mixture-model")
    {
    int ncomp = atoi(argv[1]);
    if (ncomp == 0)
      throw ConvertException("Incorrect specification of mixture model initialization");

    std::vector<double> mu, sigma;
    for(int i = 0; i < ncomp; i++)
      {
      mu.push_back(ReadIntensityValue(argv[2 + 2 * i]));
      sigma.push_back(ReadIntensityValue(argv[3 + 2 * i]));
      }

    MixtureModel<TPixel, VDim> adapter(this);
    adapter(mu, sigma);

    return 1 + 2 * ncomp;
    }

  else if (cmd == "-moments")
    {
    SizeType radius = ReadSizeVector(argv[1]);
    MomentsFeatures<TPixel, VDim>(this)(radius);
    return 1;
    }

  else if (cmd == "-mmi" || cmd == "-mattes-mutual-info")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fnf("none");
    string fnm("none");
    if (argc > 1)
      {
      fnm = argv[1];
      nret = 1;
      }
    if (argc == 3)
      {
      fnf = argv[2];
      nret = 2;
      }
    adapter("MMI", fnf.c_str(), fnm.c_str());
    return nret;
    }
  else if (cmd == "-msq" || cmd == "-mean-square")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fnf("none");
    string fnm("none");
    if (argc > 1)
      {
      fnm = argv[1];
      nret = 1;
      }
    if (argc == 3)
      {
      fnf = argv[2];
      nret = 2;
      }
    adapter("MSQ", fnf.c_str(), fnm.c_str());
    return nret;
    }
  else if (cmd == "-multiply" || cmd == "-times")
    {
    BinaryMathOperation<TPixel, VDim> adapter(this);
    adapter(BinaryMathOperation<TPixel, VDim>::MULTIPLY);
    return 0;
    }

  else if (cmd == "-ncc" || cmd == "-normalized-cross-correlation")
    {
    SizeType sz = ReadSizeVector(argv[1]);
    NormalizedCrossCorrelation<TPixel,VDim> adapter(this);
    adapter(sz);
    return 1;
    }

  else if (cmd == "-ncor" || cmd == "-normalized-correlation")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fnf("none");
    string fnm("none");
    if (argc > 1)
      {
      fnm = argv[1];
      nret = 1;
      }
    if (argc == 3)
      {
      fnf = argv[2];
      nret = 2;
      }
    adapter("NCOR", fnf.c_str(), fnm.c_str());
    return nret;
    }
  else if (cmd == "-nmi" || cmd == "-normalized-mutual-info")
    {
    ApplyMetric<TPixel, VDim> adapter(this);
    int nret = 0;
    string fnf("none");
    string fnm("none");
    if (argc > 1)
      {
      fnm = argv[1];
      nret = 1;
      }
    if (argc == 3)
      {
      fnf = argv[2];
      nret = 2;
      }
    adapter("NMI", fnf.c_str(), fnm.c_str());
    return nret;
    }

  else if (cmd == "-noise-gaussian" || cmd == "-noise")
    {
    double param = atof(argv[1]);
    ApplyNoise<TPixel, VDim> adapter(this);
    adapter(ApplyNoise<TPixel,VDim>::GAUSSIAN, param);
    return 1;
    }

  else if (cmd == "-noise-poisson" || cmd == "-noise-shot")
    {
    double param = atof(argv[1]);
    ApplyNoise<TPixel, VDim> adapter(this);
    adapter(ApplyNoise<TPixel,VDim>::POISSON, param);
    return 1;
    }

  else if (cmd == "-noise-speckle")
    {
    double param = atof(argv[1]);
    ApplyNoise<TPixel, VDim> adapter(this);
    adapter(ApplyNoise<TPixel,VDim>::SPECKLE, param);
    return 1;
    }

  else if (cmd == "-noise-salt-pepper")
    {
    double param = atof(argv[1]);
    ApplyNoise<TPixel, VDim> adapter(this);
    adapter(ApplyNoise<TPixel,VDim>::SALT_PEPPER, param);
    return 1;
    }

  else if (cmd == "-nomcs" || cmd == "-no-multicomponent-split")
    {
    m_MultiComponentSplit = false; return 0;
    }

  else if (cmd == "-nlw" || cmd == "-normwin" || cmd == "-normalize-local-window")
    {
    NormalizeLocalWindow<TPixel, VDim> adapter(this);
    SizeType radius = ReadSizeVector(argv[1]);
    adapter(radius);
    return 1;
    }

  else if (cmd == "-normpdf")
    {
    // Compute normal PDF of intensity values given sigma and mu
    double mu = atof(argv[1]);
    double s = atof(argv[2]);

    // Subtract mu
    ScaleShiftImage<TPixel, VDim> scale1(this);
    scale1(1.0, -mu);

    // Square
    m_ImageStack.push_back(m_ImageStack.back());
    MultiplyImages<TPixel, VDim> times(this);
    times();

    // Scale by -1/2s
    ScaleShiftImage<TPixel, VDim> scale2(this);
    scale2(-0.5 / s, 0.0);

    // Exponentiate
    UnaryMathOperation<TPixel, VDim> exp1(this);
    exp1(&vcl_exp);

    // Scale by factor
    ScaleShiftImage<TPixel, VDim> scale3(this);
    scale3(1.0 / sqrt(2 * vnl_math::pi * s * s), 0.0);
    return 2;
    }

  else if (cmd == "-noround")
    { m_RoundFactor = 0.0; return 0; }

  // Enable SPM extensions
  else if (cmd == "-nospm")
    { m_FlagSPM = false; return 0; }

  // Overwrite / Output command - save the image without checking if
  // it already exists.
  else if (cmd == "-o" || cmd == "-output")
    {
    WriteImage<TPixel, VDim> adapter(this);
    adapter(argv[1], true);
    return 1;
    }

  else if (cmd == "-omc" || cmd == "-output-multicomponent")
    {
    // Number of components (all by default)
    int nc, np;

    // A parameter can be optionally specified (how many components)
    RegularExpression re("^[0-9]+$");
    if (re.find(argv[1]))
      { nc = atoi(argv[1]); np=2; }
    else
      { nc = m_ImageStack.size(); np = 1; }

    // Create a writer
    WriteImage<TPixel, VDim> adapter(this);
    adapter.WriteMultiComponent(argv[np], nc);
    return np;
    }

  else if (cmd == "-oomc" || cmd == "-output-multiple-multicomponent")
    {
    // The number of components must be specified
    int nc = atoi(argv[1]);

    // Write the rest
    return 1 + this->WriteMultiple(argc-1, argv+1, nc, cmd.c_str());
    }
  
  else if (cmd == "-orient")
    {
    SetOrientation<TPixel,VDim> adapter(this);
    adapter(argv[1]);
    return 1;
    }

  // Write mulptiple images
  else if (cmd == "-oo" || cmd == "-output-multiple")
    {
    return this->WriteMultiple(argc, argv, 1, cmd.c_str());
    }

  else if (cmd == "-orient")
    {
    // Read an RAS code
    RegularExpression re("[raslpi]{3}");
    if (re.find(str_to_lower(argv[1])))
      { cout << "You supplied a RAS code" << endl; }
    else
      { cout << "I am expecting a matrix" << endl; }
    return 1;
    }


  else if (cmd == "-origin")
    {
    // The first parameter is the origin (new physical coordinate of the voxel)
    RealVector org = ReadRealVector(argv[1], true);
    m_ImageStack.back()->SetOrigin(org.data_block());
    return 1;
    }

  else if (cmd == "-origin-voxel")
    {
    // Read the physical RAS coordinate of the voxel that should be made the origin
    RealVector vec = ReadRealVector(argv[1], false);

    // This voxel should have the coordinate zero
    for (int i = 0; i <  VDim; i++)
      {
      // Here we are making the RAS/LPS switch and inverting the coordinate
      vec[i] = (i >= 2) ? -vec[i] : vec[i];
      }

    // Get physical coordinate of this voxel
    m_ImageStack.back()->SetOrigin(vec.data_block());
    return 1;
    }

  else if (cmd == "-origin-voxel-coord")
    {
    // Read the index the voxel whose coordinate we want to assign
    IndexType voxel = ReadIndexVector(argv[1]);

    // Read the new NIFTI coordinate of this position
    RealVector coord = ReadRealVector(argv[2], true);

    // Get the LPS coordinate of the voxel under the current origin
    typename ImageType::PointType lps_curr, lps_desired, lps_origin;
    m_ImageStack.back()->TransformIndexToPhysicalPoint(voxel, lps_curr);

    // Get the LPS coordinate that we want to assign to the voxel
    // The origin should be shifted by (lps_desired - lps_curr)
    for(int i = 0; i < VDim; i++)
      {
      lps_desired[i] = (i < 2) ? -coord[i] : coord[i];
      lps_origin[i] = m_ImageStack.back()->GetOrigin()[i] + lps_desired[i] - lps_curr[i];
      }

    // Assign the origin
    m_ImageStack.back()->SetOrigin(lps_origin);

    return 2;
    }

  else if (cmd == "-otsu")
    {
    OtsuThreshold<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-overlap")
    {
    double label = atof(argv[1]);
    ComputeOverlaps<TPixel, VDim> adapter(this);
    adapter(label);
    return 1;
    }

  else if (cmd == "-overlay-label-image" || cmd == "-oli")
    {
    double opacity = atof(argv[2]);
    OverlayLabelImage<TPixel, VDim> adapter(this);
    adapter(argv[1], opacity);
    return 2;
    }

  else if (cmd == "-pad")
    {
    // specify two sizes that give the padding in x,y,z
    // pad is the offset (in voxels) from the edge of the image to the
    // padExtent. Example: -pad 1x1x1vox 0x0x0vox pads on the left, posterior, inferior sides
    // by one voxel -pad 10x10x10% 10x10x10% adds 10% on all sides
    IndexType padExtentLower, padExtentUpper;

    padExtentLower = ReadIndexVector(argv[1]);
    padExtentUpper = ReadIndexVector(argv[2]);

    float padValue = atof(argv[3]);

    *verbose << "Padding image #" << m_ImageStack.size() << endl;

    PadImage<TPixel, VDim> adapter(this);
    adapter(padExtentLower, padExtentUpper, padValue);
    return 3;
    }

  else if (cmd == "-padto" || cmd == "-pad-to")
    {
    // Specify the size to which the image should be padded
    SizeType padNewSize = ReadSizeVector(argv[1]);
    float padValue = atof(argv[2]);

    // How much to add in each dimension
    IndexType padExtentLower, padExtentUpper;
    SizeType currentSize = this->PeekLastImage()->GetBufferedRegion().GetSize();
    for(int i = 0; i < VDim; i++)
      {
      // Constraint: lower + upper = desired - current
      long bilateral_pad = padNewSize[i] - currentSize[i];
      padExtentLower[i] = bilateral_pad / 2;
      padExtentUpper[i] = bilateral_pad - padExtentLower[i];
      }

    // Use the adapter
    PadImage<TPixel, VDim> adapter(this);
    adapter(padExtentLower, padExtentUpper, padValue);
    return 2;
    }

  else if (cmd == "-pca")
    {
    ComputeMoments<TPixel, VDim> adapter(this);
    adapter();
    return 0; 
    }

  else if (cmd == "-percent-intensity-mode" || cmd == "-pim")
    {
    // What does % mean when specifying intensities
    string pim = str_to_lower(argv[1]);
    if (pim == "quantile" || pim == "q")
      m_PercentIntensityMode = PIM_QUANTILE;
    else if (pim == "foregroundquantile" || pim == "fq")
      m_PercentIntensityMode = PIM_FGQUANTILE;
    else if (pim == "range" || pim == "r")
      m_PercentIntensityMode = PIM_RANGE;
    else
      throw ConvertException("Wrong -percent-intensity-mode spec %s. See help.", pim.c_str());
    return 1;
    }

  else if (cmd == "-pick")
    {
    // Create the list of images to put on the new stack
    std::vector<int> pos;
    std::vector<ImagePointer> new_stack;

    // Read all integer arguments
    for(int i = 1; i < argc; i++)
      try 
        { pos.push_back((int) myatol(argv[i])); }
      catch(...)
        { break; }

    // Retain each component
    for(int j = 0; j < pos.size(); j++)
      {
      if(pos[j] >= 0 && pos[j] < m_ImageStack.size())
        new_stack.push_back(m_ImageStack[pos[j]]);
      else if(pos[j] < 0 && -pos[j] <= m_ImageStack.size())
        new_stack.push_back(m_ImageStack[m_ImageStack.size() + pos[j]]);
      else 
        throw ConvertException("Invalid index %d in -pick command", pos[j]);
      }

    // Clear the stack
    m_ImageStack.clear();

    // Replace the stack
    for(int j = 0; j < new_stack.size(); j++)
      m_ImageStack.push_back(new_stack[j]);

    // Return number of arguments consumed
    return pos.size();
    }

  else if (cmd == "-pixel")
    {
    // Get a pixel value - no interpolation
    typename RegionType::IndexType idx = ReadIndexVector(argv[1]);
    try
      {
      double pix = m_ImageStack.back()->GetPixel(idx);
      cout << "Pixel " << idx << " has value " << pix << endl;
      }
    catch(...)
      {
      cerr << "Error: pixel " << idx << " can not be examined!" << endl;
      }
    return 1;
    }

  else if (cmd == "-pop")
    {
    *verbose << "Removing (popping) the last image from the stack" << endl;
    m_ImageStack.pop_back();
    return 0;
    }

  else if (cmd == "-popas")
    {
    string var(argv[1]);
    if (m_ImageStack.size() == 0)
      throw ConvertException("No image to assign to variable %s", var.c_str());
    m_ImageVars[var] = m_ImageStack.back();
    m_ImageStack.pop_back();
    return 1;
    }

  else if (cmd == "-probe")
    {
    // Get the probe point
    RealVector x = ReadRealVector(argv[1], true);
    SampleImage<TPixel, VDim> adapter(this);
    adapter(x);
    return 1;
    }

  else if (cmd == "-push" || cmd == "-get")
    {
    string var(argv[1]);
    typename ImageVariableMap::iterator img = m_ImageVars.find(var);
    if (img == m_ImageVars.end())
      throw ConvertException("No image assigned to variable %s", var.c_str());
    m_ImageStack.push_back(img->second);
    return 1;
    }

  else if (cmd == "-rank")
    {
    Rank<TPixel,VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-reciprocal")
    {
    // Multiply by the reciprocal (for time being at least)
    ReciprocalImage<TPixel, VDim>(this)();

    return 0;
    }

  else if (cmd == "-region")
    {
    // Get the position and index for the region
    RegionType bbox;
    bbox.SetIndex(ReadIndexVector(argv[1]));
    bbox.SetSize(ReadSizeVector(argv[2]));

    *verbose << "Extracting Subregion in #" << m_ImageStack.size() << endl;
    ExtractRegion<TPixel, VDim> adapter(this);
    adapter(bbox);
    return 2;
    }


  else if (cmd == "-reorder")
    {
    // Get the parameter, treat it as a float
    size_t k = 0;
    double k_frac = atof(argv[1]);
    if (k_frac > 0 && k_frac < 1)
      k = (size_t) (0.5 + k_frac * m_ImageStack.size());
    else if (k_frac >= 1)
      k = (size_t) atoi(argv[1]);
    else
      throw ConvertException("Parameter %s to the '-reorder' command is invalid", argv[1]);

    ReorderStack<TPixel, VDim> adapter(this);
    adapter(k);

    return 1;

    }

  else if (cmd == "-retain-labels")
    {
    // Parse the labels
    std::vector<int> vRetain;
    for(int i = 1; i < argc; i++)
      {
      try
        { vRetain.push_back((int) myatol(argv[i])); }
      catch(...)
        { break; }
      }

    // Replace the intensities with values supplie
    RetainLabels<TPixel, VDim> adapter(this);
    adapter(vRetain);

    return vRetain.size();
    }

  else if (cmd == "-rf-apply")
    {
    // Get the filename for the training output
    std::string rf_file = argv[1];

    // Get the current parameters
    RFApply<TPixel, VDim> adapter(this);
    adapter(rf_file.c_str());

    return 1;
    }

  else if (cmd == "-rf-train")
    {
    // Get the filename for the training output
    std::string rf_file = argv[1];

    // Get the current parameters
    RFTrain<TPixel, VDim> adapter(this);
    adapter(rf_file.c_str(), m_Param->m_RandomForest);

    return 1;
    }

  else if (cmd == "-rf-param-patch")
    {
    SizeType patch_radius = ReadSizeVector(argv[1]);
    m_Param->m_RandomForest.patch_radius = patch_radius;
    return 1;
    }

  else if (cmd == "-rf-param-usexyz")
    {
    m_Param->m_RandomForest.use_coordinate_features = true;
    return 0;
    }

  else if (cmd == "-rf-param-nousexyz")
    {
    m_Param->m_RandomForest.use_coordinate_features = false;
    return 0;
    }

  else if (cmd == "-rf-param-ntrees")
    {
    m_Param->m_RandomForest.forest_size = atoi(argv[1]);
    return 1;
    }

  else if (cmd == "-rf-param-treedepth")
    {
    m_Param->m_RandomForest.tree_depth = atoi(argv[1]);
    return 1;
    }

  else if (cmd == "-set-sform")
    {
    string fn_tran( argv[1] );

    *verbose << "Setting sform of image #" << m_ImageStack.size() << endl;
    SetSform<TPixel, VDim> adapter(this);
    adapter( fn_tran );

    return 1;
    }

  else if (cmd == "-replace")
    {
    vector<double> vReplace;

    // Read a list of numbers from the command line
    for(int i = 1; i < argc; i++)
      {
      try
        { vReplace.push_back(myatof(argv[i])); }
      catch(...)
        { break; }
      }

    // Make sure the number of replacement rules is even
    if (vReplace.size() % 2 == 1)
      {
      cerr << "The number of parameters to '-replace' must be even!" << endl;
      throw -1;
      }

    // Replace the intensities with values supplie
    ReplaceIntensities<TPixel, VDim> adapter(this);
    adapter(vReplace);

    return vReplace.size();
    }

  // Resample command. Retain the bounding box of the image
  // while changing voxel size
  else if (cmd == "-resample")
    {
    SizeType sz = ReadSizeVector(argv[1]);
    ResampleImage<TPixel, VDim> adapter(this);
    adapter(sz);
    return 1;
    }

  else if (cmd == "-resample-iso")
    {
    typename ImageType::SpacingType spc =  this->PeekImage(0)->GetSpacing();
    std::string mode = argv[1];
    double isospc = 0.0;
    if(mode == "min" || mode == "MIN")
      isospc = spc.GetVnlVector().min_value();
    else if(mode == "max" || mode == "MAX")
      isospc = spc.GetVnlVector().max_value();
    else
      throw ConvertException("Unsupported mode '%s' for command %s", argv[1], cmd.c_str());

    SizeType sz = this->PeekImage(0)->GetBufferedRegion().GetSize();
    for(size_t i = 0; i < VDim; i++)
      sz[i] = static_cast<size_t>((0.5 + sz[i] * spc[i] / isospc));
    ResampleImage<TPixel, VDim> adapter(this);
    adapter(sz);
    return 1;
    }

  else if (cmd == "-resample-mm")
    {
    RealVector vox = ReadRealSize(argv[1]);
    SizeType sz = m_ImageStack.back()->GetBufferedRegion().GetSize();
    for(size_t i = 0; i < VDim; i++)
      {
      sz[i] = static_cast<size_t>((0.5 + sz[i] * m_ImageStack.back()->GetSpacing()[i]) / vox[i]);
      }
    ResampleImage<TPixel, VDim> adapter(this);
    adapter(sz);
    return 1;
    }

  else if (cmd == "-reslice-itk")
    {
    string fn_tran(argv[1]);
    ResliceImage<TPixel, VDim> adapter(this);
    adapter("itk", fn_tran);
    return 1;
    }

  else if (cmd == "-reslice-matrix")
    {
    string fn_tran(argv[1]);
    ResliceImage<TPixel, VDim> adapter(this);
    adapter("matrix", fn_tran);
    return 1;
    }

  else if (cmd == "-reslice-identity")
    {
    ResliceImage<TPixel, VDim> adapter(this);
    adapter("identity", "");
    return 0;
    }

  else if (cmd == "-rgb2hsv")
    {
    VoxelwiseComponentFunction<TPixel, VDim> adapter(this);
    adapter("rgb2hsv");
    return 0;
    }

  else if (cmd == "-rms")
    { m_Param->m_AntiAliasRMS = atof(argv[1]); return 1; }

  else if (cmd == "-round")
    { m_RoundFactor = 0.5; return 0; }

  // No else if here because of a windows compiler error (blocks nested too deeply)
  if (cmd == "-scale")
    {
    double factor = atof(argv[1]);
    ScaleShiftImage<TPixel, VDim> adapter(this);
    adapter(factor, 0.0);
    return 1;
    }


  else if (cmd == "-set-sform")
    {
    string fn_tran( argv[1] );

    *verbose << "Setting sform of image #" << m_ImageStack.size() << endl;
    SetSform<TPixel, VDim> adapter(this);
    adapter( fn_tran );

    return 1;
    }

  else if (cmd == "-sin")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_sin);
    return 0;
    }

  else if (cmd == "-slice")
    {
    string axis( argv[1] );
    char * pos = argv[2];

    ExtractSlice<TPixel, VDim> adapter(this);
    adapter(axis, pos, 1);

    return 2;
    }

  else if (cmd == "-slice-all")
    {
    string axis( argv[1] );
    char * pos = argv[2];

    ExtractSlice<TPixel, VDim> adapter(this);
    adapter(axis, pos, this->GetStackSize());

    return 2;
    }


  else if (cmd == "-sharpen")
    {
    LaplacianSharpening<TPixel,VDim> adapter(this);
    adapter();

    return 0;
    }

  else if (cmd == "-shift")
    {
    double x = atof(argv[1]);
    ScaleShiftImage<TPixel, VDim> adapter(this);
    adapter(1.0, x);
    return 1;
    }

  else if (cmd == "-signed-distance-transform" || cmd == "-sdt")
    {
    SignedDistanceTransform<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  // Gaussian smoothing command
  else if (cmd == "-smooth")
    {
    RealVector stdev = ReadRealSize(argv[1]);
    SmoothImage<TPixel, VDim> adapter(this);
    adapter(stdev, false);
    return 1;
    }

  else if (cmd == "-smooth-fast")
    {
    RealVector stdev = ReadRealSize(argv[1]);
    SmoothImage<TPixel, VDim> adapter(this);
    adapter(stdev, true);
    return 1;
    }

  else if (cmd == "-spacing")
    {
    RealVector org = ReadRealSize(argv[1]);
    m_ImageStack.back()->SetSpacing(org.data_block());
    return 1;
    }

  else if (cmd == "-split")
    {
    SplitMultilabelImage<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  else if (cmd == "-sqrt")
    {
    UnaryMathOperation<TPixel, VDim> adapter(this);
    adapter(&vcl_sqrt);
    return 0;
    }

  else if (cmd == "-staple")
    {
    // Perform the STAPLE algorithm on the data
    double ival = atof(argv[1]);
    StapleAlgorithm<TPixel, VDim> adapter(this);
    adapter(ival);
    return 1;
    }

  else if (cmd == "-steig" || cmd == "-structure-tensor-eigenvalues")
    {
    double scale_grad = atof(argv[1]);
    double scale_window = atof(argv[2]);
    StructureTensorEigenValues<TPixel,VDim> adapter(this);
    adapter(scale_grad, scale_window);

    return 2;
    }

  // Enable SPM extensions
  else if (cmd == "-spm")
    { m_FlagSPM = true; return 0; }

  else if (cmd == "-subtract")
    {
    BinaryMathOperation<TPixel, VDim> adapter(this);
    adapter(BinaryMathOperation<TPixel, VDim>::SUBTRACT);
    return 0;
    }

  else if (cmd == "-supervoxel" || cmd == "-sv")
    {
    SLICSuperVoxel<TPixel,VDim> adapter(this);
    int samples = atoi(argv[1]);
    double m = atof(argv[2]);
    adapter(samples, m);
    return 2;
    }

  // Stretch the intensity range
  else if (cmd == "-stretch")
    {
    double u1 = ReadIntensityValue(argv[1]);
    double u2 = ReadIntensityValue(argv[2]);
    double v1 = ReadIntensityValue(argv[3]);
    double v2 = ReadIntensityValue(argv[4]);
    double a = (v2 - v1) / (u2 - u1);
    double b = v1 - u1 * a;
    ScaleShiftImage<TPixel, VDim> adapter(this);
    adapter(a, b);
    return 4;
    }

  else if (cmd == "-swapdim") 
    {
    // For now we only support a single string
    std::vector<std::string> code;
    code.push_back(argv[1]);
    SwapDimensions<TPixel, VDim> adapter(this);
    adapter(code);

    return 1;
    }


  // Test image equality
  else if (cmd == "-test-image")
    {
    // Check if the next argument is a tolerance value
    double tol = 1e-8;
    int np = 0;
    RegularExpression re("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)$");
    if (argc > 1 && re.find(argv[1]))
      { tol = atoi(argv[1]); np++; }

    TestImage<TPixel, VDim> adapter(this);
    adapter(true, true, tol);
    return np;
    }

  else if (cmd == "-test-probe")
    {
    // Just like the probe command
    RealVector x = ReadRealVector(argv[1], true);
    double v_test = atof(argv[2]);
    double tol = 1e-8;
    int np = 2;

    // Read the optional tolerance value
    RegularExpression re("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)$");
    if (argc > 1 && re.find(argv[1]))
      { tol = atoi(argv[1]); np++; }

    // Probe
    SampleImage<TPixel, VDim> adapter(this);
    adapter(x);

    // Compare and exit
    double abs_diff = std::fabs(v_test - adapter.GetResult());
    if(abs_diff > tol)
      {
      cout << "Probe test failed, absolute difference = " << abs_diff << endl;
      exit(1);
      }
    else
      exit(0);
    }

  // Thresholding
  else if (cmd == "-threshold" || cmd == "-thresh")
    {
    double u1 = strcmp(argv[1],"-inf") == 0 ? -vnl_huge_val(0.0) : ReadIntensityValue(argv[1]);
    double u2 = strcmp(argv[2],"inf") == 0 ? vnl_huge_val(0.0) : ReadIntensityValue(argv[2]);
    double v1 = ReadIntensityValue(argv[3]);
    double v2 = ReadIntensityValue(argv[4]);
    ThresholdImage<TPixel, VDim> adapter(this);
    adapter(u1, u2, v1, v2);
    return 4;
    }

  // Tiling
  else if (cmd == "-tile")
    {
    TileImages<TPixel, VDim> adapter(this);
    adapter(std::string(argv[1]));
    return 1;
    }

  // Trim the image (trim background values from the margin)
  else if (cmd == "-trim")
    {
    // Read the size of the wrap region
    RealVector margin = ReadRealSize(argv[1]);

    // Trim the image accordingly
    TrimImage<TPixel, VDim> adapter(this);
    adapter(margin, TrimImage<TPixel, VDim>::SPECIFY_MARGIN);

    // Return the number of arguments consumed
    return 1;
    }

  else if (cmd == "-trim-to-size")
    {
    // Read the size of the trim region
    RealVector size = ReadRealSize(argv[1]);

    // Trim the image accordingly
    TrimImage<TPixel, VDim> adapter(this);
    adapter(size, TrimImage<TPixel, VDim>::SPECIFY_FINALSIZE);

    // Return the number of arguments consumed
    return 1;
    }

  // Output type specification
  else if (cmd == "-type")
    {
    m_TypeId = str_to_lower(argv[1]);
    return 1;
    }

  // Verbose mode
  else if (cmd == "-verbose")
    { verbose = &std::cout; return 0; }

  else if (cmd == "-noverbose")
    { verbose = &devnull; return 0; }

  else if (cmd == "-version")
    {
    cout << "Version " << ImageConverter_VERSION_STRING << endl;
    return 0;
    }

  else if (cmd == "-vote")
    {
    Vote<TPixel, VDim> adapter(this);
    adapter(false);
    return 0;
    }

  else if (cmd == "-vote-mrf")
    {
    std::string mode_str = str_to_lower(argv[1]);
    double beta = atof(argv[2]);

    typedef MRFVote<TPixel, VDim> Adapter;
    typename Adapter::Mode mode;
    if(mode_str == "votes_against" || mode_str == "va")
      mode = Adapter::VOTES_AGAINST;
    else if(mode_str == "log_likelihood" || mode_str == "ll")
      mode = Adapter::LOG_LIKELIHOOD;
    else
      throw ConvertException("Unknown mode parameter %s to -vote-mrf", mode_str.c_str());

    Adapter adapter(this);
    adapter(mode, beta);
    return 2;
    }

  else if (cmd == "-vote-label")
    {
    UpdateMetadataKey<TPixel, VDim> adapter(this);
    adapter("CND:VOTE_LABEL",argv[1]);
    return 1;
    }

  else if (cmd == "-voxel-sum")
    {
    // Simply print the sum of all voxels in the image
    double sum = 0;
    size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
    for(size_t i = 0; i < n; i++)
      sum += m_ImageStack.back()->GetBufferPointer()[i];
    cout << "Voxel Sum: " << sum << endl;
    return 0;
    }

  else if (cmd == "-voxel-integral" || cmd == "-voxel-int")
    {
    // Simply print the sum of all voxels in the image
    double sum = 0;
    size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
    for(size_t i = 0; i < n; i++)
      sum += m_ImageStack.back()->GetBufferPointer()[i];
    double vol = 1.0;
    for(size_t d = 0; d < VDim; d++)
      vol *= m_ImageStack.back()->GetSpacing()[d];
    cout << "Voxel Integral: " << sum * vol << endl;
    return 0;
    }

  else if (cmd == "-voxelwise-regression" || cmd == "-voxreg")
    {
    // Get the order
    size_t order = atoi(argv[1]);
    VoxelwiseRegression<TPixel, VDim> adapter(this);
    adapter(order);
    return 1;
    }

  // Warp image
  else if (cmd == "-warp")
    {
    WarpImage<TPixel, VDim> adapter(this);
    adapter();
    return 0;
    }

  // Warp label image
  else if (cmd == "-warp-label" || cmd=="-warplabel" || cmd=="-wl")
    {
    RealVector stdev = ReadRealSize(argv[1]);
    WarpLabelImage<TPixel, VDim> adapter(this);
    adapter(stdev);
    return 1;
    }

  // Wrap image around
  else if (cmd == "-wrap")
    {
    IndexType iwrap = ReadIndexVector(argv[1]);
    WrapDimension<TPixel, VDim> adapter(this);
    adapter(iwrap);
    return 1;
    }

  else if (cmd == "-weighted-sum" || cmd == "-wsum")
    {
    std::vector<double> weights;
    for(int i = 1; i < argc; i++)
      if (is_double(argv[i]))
        weights.push_back(atof(argv[i]));
      else break;
    WeightedSum<TPixel,VDim> adapter(this);
    adapter(weights);
    return weights.size();
    }

  else if (cmd == "-weighted-sum-voxelwise" || cmd == "-wsv")
    {
    WeightedSumVoxelwise<TPixel,VDim> adapter(this);
    adapter();
    return 0;
    }

  // Unknown command
  else
    { cerr << "Unknown command " << cmd << endl; throw -1; }

  cerr << "Fell through!" << endl;
  throw -1;
}


template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::ProcessCommandLine(int argc, char *argv[])
{
  // Disable multithreading
  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);

  // The last command
  std::string lastCommand;

  // Check the number of arguments
  if (argc == 1)
    {
    cerr << "PICSL convert3d tool - from the creators of ITK-SNAP" << endl;
    cerr << "For full documentation and usage examples, see" << endl;
    cerr << "    http://www.itksnap.org/c3d" << endl;
    cerr << "To get help on available commands, call" << endl;
    cerr << "    " << argv[0] << " -h" << endl;
    return -1;
    }

  // Try processing command line
  try
    {
    // The last parameter in the command line is the output file
    string fnOutput = argv[argc-1];

    // Command line processing
    for(int i = 1; i < argc; ++i)
      {
      string cmd = argv[i];
      if (cmd[0] == '-')
        {
        // Save the last command (for exceptions, etc)
        lastCommand = argv[i];

        // A command has been supplied
        i += ProcessCommand(argc-i, argv+i);
        }
      else
        {
        lastCommand = "";

        // An image file name has been provided. If this image is followed by commands
        // read it and push in the pipeline.
        if (i != argc-1)
          {
          ReadImage<TPixel, VDim> adapter(this);
          adapter(argv[i]);
          }
        else
          {
          // Write the image, but in safe mode
          WriteImage<TPixel, VDim> adapter(this);
          adapter(argv[i], false);
          }
        }
      }
    return 0;
    }
  catch (StackAccessException &)
    {
    cerr << "Not enough images on the stack for the requested command." << endl;
    cerr << "  Requested command: " << lastCommand << endl;
    cerr << "  Note: C3D requires image operands to precede commands." << endl;
    cerr << "        message can be caused by incorrect usage, such as" << endl;
    cerr << "           c3d -command image.nii " << endl;
    cerr << "        instead of " << endl;
    cerr << "           c3d image.nii -command" << endl;
    return -1;
    }
  
  catch (std::exception &exc)
    {
    cerr << "Exception caught of type " << typeid(exc).name() << endl;
    if(lastCommand.size())
      cerr << "  When processing command: " << lastCommand << endl;
    cerr << "  Exception detail: " << exc.what() << endl;
    return -1;
    }
  catch (...)
    {
    cerr << "Unknown exception caught by convert3d" << endl;
    if(lastCommand.size())
      cerr << "  When processing command: " << lastCommand << endl;
    return -1;
    }
}

bool str_ends_with(const std::string &s, const std::string &pattern)
{
  size_t ipos = s.rfind(pattern);
  return(ipos == s.size() - pattern.size());
}

// How the specification is made
enum VecSpec { PHYSICAL, VOXELS, PERCENT, NONE };

template<unsigned int VDim>
void ReadVecSpec(const char *vec_in, vnl_vector_fixed<double,VDim> &vout, VecSpec &type)
{
  // Set up the regular expressions for numerical string parsing
  RegularExpression re1(
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)(mm|vox|%)?");
  RegularExpression re2(
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)x"
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)(mm|vox|%)?");
  RegularExpression re3(
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)x"
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)x"
    "([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)(mm|vox|%)?");

  // Lowercase string
  string vec = str_to_lower(vec_in);
  string spec;

  // Check if it's a single-component specification
  if (VDim == 2 && re2.find(vec))
    {
    vout[0] = atof(re2.match(1).c_str());
    vout[1] = atof(re2.match(3).c_str());
    spec = re2.match(5);
    }
  else if (VDim == 3 && re3.find(vec))
    {
    vout[0] = atof(re3.match(1).c_str());
    vout[1] = atof(re3.match(3).c_str());
    vout[2] = atof(re3.match(5).c_str());
    spec = re3.match(7);
    }
  else if (re1.find(vec))
    {
    vout.fill(atof(re1.match(1).c_str()));
    spec = re1.match(3);
    }
  else throw ConvertException("Invalid vector specification %s", vec_in);

  // Get the type of spec. Luckily, all suffixes have the same length
  switch(spec.length()) {
  case 0: type = NONE; break;
  case 1: type = PERCENT; break;
  case 2: type = PHYSICAL; break;
  case 3: type = VOXELS; break;
  default: throw ConvertException("Internal error in VecSpec code");
  }
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::RealVector
ImageConverter<TPixel, VDim>
::ReadRealVector(const char *vec_in, bool is_point)
{
  // Output vector
  RealVector x, scale, offset;
  VecSpec type;

  // Read the vector
  ReadVecSpec(vec_in, x, type);

  // Check the type of the vector
  if (type != VOXELS && type != PHYSICAL && type != PERCENT)
    throw ConvertException(
      "Invalid vector spec %s (must end with 'mm' or 'vox' or '%' )", vec_in);

  // If in percent, scale by voxel size
  if (type == PERCENT)
    {
    for(size_t i = 0; i < VDim; i++)
      x[i] *= m_ImageStack.back()->GetBufferedRegion().GetSize()[i] / 100.0;
    type = VOXELS;
    }

  // If the vector is in vox units, map it to physical units
  if (type == VOXELS)
    {
    // Get the matrix
    typename ImageType::TransformMatrixType MP =
      m_ImageStack.back()->GetVoxelSpaceToRASPhysicalSpaceMatrix();

    // Create the vector to multiply by
    vnl_vector_fixed<double, VDim+1> X, XP;
    for(size_t d = 0; d < VDim; d++)
      X[d] = x[d];
    X[VDim] = is_point ? 1.0 : 0.0;

    // Apply matrix
    XP = MP * X;
    for(size_t d = 0; d < VDim; d++)
      x[d] = XP[d];
    }

  return x;
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::RealVector
ImageConverter<TPixel, VDim>
::ReadRealSize(const char *vec_in)
{
  RealVector x = ReadRealVector(vec_in, false);
  for(size_t d = 0; d < VDim; d++)
    x[d] = fabs(x[d]);
  return x;
}

template<class TPixel, unsigned int VDim>
TPixel
ImageConverter<TPixel, VDim>
::ReadIntensityValue(const char *vec)
{
  // Check if the input is infinity first
  if (!strcmp(vec, "inf") || !strcmp(vec, "+inf") || !strcmp(vec, "Inf") || !strcmp(vec, "+Inf"))
    return vnl_huge_val(0.0);
  else if (!strcmp(vec, "-inf") || !strcmp(vec, "-Inf"))
    return -vnl_huge_val(0.0);

  // Read the double part
  char *endptr;

  // Read the floating point part
  TPixel val = (TPixel) strtod(vec, &endptr);

  // Check validity
  if (endptr == vec)
    throw ConvertException("Can't convert %s to an intensity spec", vec);

  // Check if there is a '%' specification
  if (*endptr == '%')
    {
    if (m_PercentIntensityMode == PIM_QUANTILE || m_PercentIntensityMode == PIM_FGQUANTILE)
      {
      // Check valid quantile
      if (val < 0.0 || val > 100.0)
        throw ConvertException("Invalid quantile spec %s, must be between 0 and 100", vec);

      // Compute the specified quantile
      double qtile = 0.01 * val;
      if (m_ImageStack.size() == 0)
        throw ConvertException("Can't use intensity quantile spec with no image on stack");

      // Allocate an array for sorting
      size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
      TPixel *asort = new TPixel[n], *p = asort, *q = m_ImageStack.back()->GetBufferPointer();

      // Copy intensity values that are legit, ignore nans
      for(size_t i = 0; i < n; i++, q++)
        {
        // We don't include nans and if FGQUANTILE, background values
        if (!vnl_math_isnan(*q))
          if (m_PercentIntensityMode == PIM_QUANTILE || *q != m_Background)
            {*p = *q; ++p;}
        }

      // Get the size of the sort array
      size_t np = p - asort;
      if (np == 0)
        {
        if (m_PercentIntensityMode == PIM_QUANTILE)
          throw ConvertException("Quantile could not be computed because the image has only NANs");
        else
          throw ConvertException("Foreground quantile could not be computed because the image has only background");
        }

      // Sort the acceptable values
      std::sort(asort, p);

      // Get the quantile
      size_t k = (size_t) (qtile * np);
      val = asort[k];
      delete[] asort;

      if (m_PercentIntensityMode == PIM_QUANTILE)
        *verbose << "Quantile " << qtile << " maps to " << val << endl;
      else
        *verbose << "Foreground quantile " << qtile <<  " (over "
        << np << " voxels) maps to " << val << endl;
      }
    else
      {
      // A range specification. We need the min and max of the image
      double imin = numeric_limits<double>::max();
      double imax = - numeric_limits<double>::max();
      size_t n = m_ImageStack.back()->GetBufferedRegion().GetNumberOfPixels();
      TPixel *q = m_ImageStack.back()->GetBufferPointer();
      for(size_t i = 0; i < n; i++, q++)
        {
        if (*q < imin) imin = *q;
        if (*q > imax) imax = *q;
        }

      double rspec = val * 0.01;
      val = imin + rspec * (imax - imin);
      *verbose << "Intensity range spec " << rspec << " maps to " << val << endl;
      }
    }

  return val;
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::SizeType
ImageConverter<TPixel, VDim>
::ReadSizeVector(const char *vec_in)
{
  // Create a copy of the input string
  char *vec = new char[strlen(vec_in) + 1];
  strcpy(vec, vec_in);

  size_t i;

  typename ImageType::SizeType sz;

  // Check if the string ends with %
  if (str_ends_with(vec, "%"))
    {
    // Read floating point size
    RealVector factor;
    char *tok = strtok(vec, "x%");
    for(i = 0; i < VDim && tok != NULL; i++)
      {
      factor[i] = atof(tok);
      if (factor[i] < 0)
        throw ConvertException("Negative percent size specification: %s", vec_in);
      tok = strtok(NULL, "x%");
      }

    if (i == 1)
      factor.fill(factor[0]);

    // Get the size of the image in voxels
    for(size_t i = 0; i < VDim; i++)
      sz[i] = (unsigned long)(m_ImageStack.back()->GetBufferedRegion().GetSize(i) * 0.01 * factor[i] + 0.5);
    }
  else
    {
    // Find all the 'x' in the string
    char *tok = strtok(vec, "x");
    for(size_t i = 0; i < VDim; i++)
      {
      if (tok == NULL)
        throw ConvertException("Invalid size specification: %s", vec_in);
      int x = atoi(tok);
      if (x < 0)
        throw ConvertException("Negative size specification: %s", vec_in);
      sz[i] = (unsigned long)(x);
      tok = strtok(NULL, "x");
      }
    }

  delete[] vec;
  return sz;
}


template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::IndexType
ImageConverter<TPixel, VDim>
::ReadIndexVector(const char *vec_in)
{
  // Create a copy of the input string
  char *vec = new char[strlen(vec_in)];
  strcpy(vec, vec_in);

  size_t i;

  typename ImageType::IndexType idx;

  // Check if the string ends with %
  if (str_ends_with(vec, "%"))
    {
    // Read floating point size
    RealVector factor;
    char *tok = strtok(vec, "x%");
    for(i = 0; i < VDim && tok != NULL; i++)
      {
      factor[i] = atof(tok);
      tok = strtok(NULL, "x%");
      }

    if (i == 1)
      factor.fill(factor[0]);

    // Get the size of the image in voxels
    for(size_t i = 0; i < VDim; i++)
      idx[i] = (long)(m_ImageStack.back()->GetBufferedRegion().GetSize(i) * 0.01 * factor[i] + 0.5);
    }
  else
    {
    // Find all the 'x' in the string
    char *tok = strtok(vec, "x");
    for(size_t i = 0; i < VDim; i++)
      {
      if (tok == NULL)
        throw ConvertException("Invalid index specification: %s", vec_in);
      int x = atoi(tok);
      idx[i] = (long)(x);
      tok = strtok(NULL, "x");
      }
    }

  delete[] vec;
  return idx;
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::CopyImage()
{
  // Get the input image
  ImagePointer input = m_ImageStack.back();

  // Simply make a copy of the input image on the stack
  ImagePointer output = ImageType::New();
  output->SetRegions(input->GetBufferedRegion());
  output->SetSpacing(input->GetSpacing());
  output->SetOrigin(input->GetOrigin());
  output->SetDirection(input->GetDirection());
  output->SetMetaDataDictionary(input->GetMetaDataDictionary());
  output->Allocate();

  // Copy everything
  size_t n = input->GetBufferedRegion().GetNumberOfPixels();
  for(size_t i = 0; i < n; i++)
    output->GetBufferPointer()[i] = input->GetBufferPointer()[i];

  // Put on the end of the stack
  m_ImageStack.pop_back();
  m_ImageStack.push_back(output);
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::GetBoundingBox(ImageType *image, RealVector &bb0, RealVector &bb1)
{
  for(size_t i = 0; i < VDim; i++)
    {
    bb0[i] = image->GetOrigin()[i];
    bb1[i] = bb0[i] + image->GetSpacing()[i] * image->GetBufferedRegion().GetSize()[i];
    }
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PrintMatrix(std::ostream &sout, vnl_matrix<double> mat, const char *fmt, const char *prefix)
{
  // Print each row and column of the matrix
  char buffer[256];
  for(size_t i = 0; i < mat.rows(); i++)
    {
    sout << prefix;
    for(size_t j = 0; j < mat.columns(); j++)
      {
      sprintf(buffer, fmt, mat(i,j));
      sout << buffer;
      }
    sout << endl;
    }
}


template<class TPixel, unsigned int VDim>
size_t
ImageConverter<TPixel, VDim>
::ForEachComponentLoop(int ncomp, int argc, char *argv[])
{
  size_t narg = 0;

  // Allow looping by components, each run of the loop applies to all the 1st, 2nd, 3rd, etc
  // components of a multicomponent image. 

  // Back up the current stack
  ImageStack<ImageType> in_stack = m_ImageStack, out_stack;

  // Print out what's going on
  *verbose << "Repeating commands for all " << in_stack.size() << " images in batches of " << ncomp << endl;

  // Make sure the stack is divisible
  if(in_stack.size() % ncomp != 0)
    throw ConvertException("Number of images on the stack (%d) is not divisible by stride (%d)",
      in_stack.size(), ncomp);

  // Loop over the component
  for(int comp = 0; comp < ncomp; comp++)
    {
    // Put the column of component images on the stack
    m_ImageStack.clear();
    for(int j = comp; j < in_stack.size(); j+=ncomp)
      {
      m_ImageStack.push_back(in_stack[j]);
      }

    // Set the in-loop flag
    m_LoopType = LOOP_FOREACH;

    // When the -endfor is encountered, the InLoop flag will be switched
    narg = 1;
    while(m_LoopType == LOOP_FOREACH)
      narg += 1 + this->ProcessCommand(argc-narg, argv+narg);

    // Place the result (if any) on the output stack
    for(int k = 0; k < m_ImageStack.size(); k++)
      out_stack.push_back(m_ImageStack.back());
    }

  // Update the stack
  m_ImageStack = out_stack;

  // Return the number of arguments to the next command
  return narg - 1;
}

template<class TPixel, unsigned int VDim>
size_t
ImageConverter<TPixel, VDim>
::ForEachLoop(int argc, char *argv[])
{
  size_t narg = 0;

  // Note: this is a rather lame implementation that uses recursion with
  // a state variable to repeat a range of command a bunch of times

  // Back up the current stack
  ImageStack<ImageType> in_stack = m_ImageStack, out_stack;

  // Print out what's going on
  *verbose << "Repeating commands for all " << in_stack.size() << " images" << endl;

  // Loop over all images
  for(size_t i = 0; i < in_stack.size(); i++)
    {
    narg = 1;

    // Set up the image stack
    m_ImageStack.clear();
    m_ImageStack.push_back(in_stack[i]);

    // Set the in-loop flag
    m_LoopType = LOOP_FOREACH;

    // When the -endfor is encountered, the InLoop flag will be switched
    while(m_LoopType == LOOP_FOREACH)
      narg += 1 + this->ProcessCommand(argc-narg, argv+narg);

    // Place the result (if any) on the output stack
    if (m_ImageStack.size() > 1)
      throw ConvertException("Commands in the -foreach clause may not produce multiple outputs");
    else if (m_ImageStack.size() == 1)
      out_stack.push_back(m_ImageStack.back());
    }

  // Update the stack
  m_ImageStack = out_stack;

  // Return the number of arguments to the next command
  return narg - 1;
}

/**
 * The -accum function allows us to apply binary operations like -add to a list of 
 * images. The commands inside the -accum/-endaccum block are repeated for consecutive
 * pairs of images, i.e., 1 2 , then result(1,2) and 3, and so on. Within each block,
 * there are always two images on the stack, and the block must produce one image
 */
template<class TPixel, unsigned int VDim>
size_t
ImageConverter<TPixel, VDim>
::AccumulateLoop(int argc, char *argv[])
{
  size_t narg = 0;

  // If there are less than two images on the stack, the accum command will not be run.
  if (m_ImageStack.size() < 2)
    {
    throw ConvertException("Too few images on the stack for the -accum command, two or more images are required!");
    }

  // Back up the current stack
  ImageStack<ImageType> in_stack = m_ImageStack;

  // Print out what's going on
  *verbose << "Accumulating result of binary operation for all " << in_stack.size() << " images" << endl;

  // Pop the last image off the stack and place it in the temporary stack
  m_ImageStack.clear();
  m_ImageStack.push_back(in_stack.back());
  in_stack.pop_back();

  // Loop until input stack is empty
  while(in_stack.size())
    {
    // Keep track of the position in the argument list
    narg = 1;

    // Push the next image on the stack
    m_ImageStack.push_back(in_stack.back());
    in_stack.pop_back();

    // Set the in-loop flag
    m_LoopType = LOOP_ACCUM;

    // When the -endfor is encountered, the InLoop flag will be switched
    while(m_LoopType == LOOP_ACCUM)
      narg += 1 + this->ProcessCommand(argc-narg, argv+narg);

    // Place the result (if any) on the output stack
    if (m_ImageStack.size() != 1)
      throw ConvertException("Commands in the -accum clause must produce exactly one output");
    }

  // Return the number of arguments to the next command
  return narg - 1;
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::LabelToRGBAMap
ImageConverter<TPixel, VDim>
::ReadLabelToRGBAMap(const char *fname)
{
  // Create a stream for reading the file
  ifstream fin(fname);
  string line;

  // Create a temporary array of color labels
  LabelToRGBAMap lmap;

  // Check that the file is readable
  if (!fin.good())
    throw ConvertException("Label file %s can not be read", fname);

  // Read each line of the file separately
  for(size_t iLine=0; !fin.eof(); iLine++)
    {
    // Read the line into a string
    std::getline(fin,line);

    // Check if the line is a comment or a blank line
    if (line[0] == '#' || line.length() == 0)
      continue;

    // Create a stream to parse that string
    istringstream iss(line);
    iss.exceptions(std::ios::badbit | std::ios::failbit);

    try
      {
      // Read in the elements of the file
      vnl_vector_fixed<double, 4> rgba;
      double idx;
      iss >> idx;
      iss >> rgba[0];
      iss >> rgba[1];
      iss >> rgba[2];
      iss >> rgba[3];

      // Add the color label to the list
      lmap[idx] = rgba;
      }
    catch( std::exception )
      {
      // create an exeption string
      throw ConvertException("Error reading file %s on line %d", fname, iLine+1);
      }
    }

  // Return the map
  return lmap;
}

template<class TPixel, unsigned int VDim>
bool
ImageConverter<TPixel, VDim>
::CheckStackSameDimensions(size_t n)
{
  if (n == 0)
    n = m_ImageStack.size();

  if (m_ImageStack.size() < n || n == 0)
    throw ConvertException("Too few images on the stack for consistency check");

  size_t top = n - 1;
  for(size_t i = 0; i < n; i++)
    {
    if (m_ImageStack[top - i]->GetBufferedRegion() != m_ImageStack[top]->GetBufferedRegion())
      {
      return false;
      }
    }

  return true;
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::SetVariable(std::string name, ImagePointer image)
{
  m_ImageVars[name] = image;
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::ImageType *
ImageConverter<TPixel, VDim>
::GetVariable(std::string name)
{
  if(m_ImageVars.find(name) != m_ImageVars.end())
    return m_ImageVars[name];
  else
    return NULL;
}


template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel, VDim>::ImagePointer
ImageConverter<TPixel, VDim>
::PopImage()
{
  if(m_ImageStack.size() == 0)
    throw ConvertException("Attempted to pop an image from empty stack");

  ImagePointer stack_end = m_ImageStack.back();
  m_ImageStack.pop_back();

  return stack_end;
}

template<class TPixel, unsigned int VDim>
std::vector<typename ImageConverter<TPixel, VDim>::ImagePointer>
ImageConverter<TPixel, VDim>
::PopNImages(unsigned int n)
{
  if(m_ImageStack.size() < n)
    throw ConvertException("Attempted to pop %d images from a stack of %d images", n, m_ImageStack.size());

  std::vector<typename ImageConverter<TPixel, VDim>::ImagePointer> retvec(n);
  for(int i = n-1; i >= 0; i--)
    retvec[i] = this->PopImage();

  return retvec;
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PushImage(ImageType *image)
{
  m_ImageStack.push_back(image);
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel,VDim>::ImageType *
ImageConverter<TPixel, VDim>
::PopAndPushCopy()
{
  // Pop the old image from the stack
  ImagePointer img_old = this->PopImage();

  // Create a new image that is a copy of the old
  ImagePointer img_new = ImageType::New();
  img_new->CopyInformation(img_old);
  img_new->SetRegions(img_old->GetBufferedRegion());
  img_new->Allocate();

  // Copy the intensities
  long nvox = img_new->GetPixelContainer()->Size();
  TPixel *p_new = img_new->GetBufferPointer(), *p_old = img_old->GetBufferPointer();
  for(long i = 0; i < nvox; i++)
    p_new[i] = p_old[i];

  // Push the copy on the stack
  this->PushImage(img_new);

  // Return the last image
  return this->PeekLastImage();
}


template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel,VDim>::ImageType *
ImageConverter<TPixel, VDim>
::PeekImage(int k)
{
  if(m_ImageStack.size() <= k || k < 0)
    throw ConvertException("Attempted to access image outside of stack range");
  return m_ImageStack[k];
}

template<class TPixel, unsigned int VDim>
typename ImageConverter<TPixel,VDim>::ImageType *
ImageConverter<TPixel, VDim>
::PeekLastImage()
{
  return this->PeekImage(this->GetStackSize() - 1);
}

template<class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::GetStackSize()
{
  return m_ImageStack.size();
}

template<class TPixel, unsigned int VDim>
void
ImageConverter<TPixel, VDim>
::PrintF(const char *fmt, ...)
{
  char buffer[4096];
  va_list args;
  va_start (args, fmt);
  vsprintf (buffer, fmt, args);
  va_end (args);

  *this->verbose << buffer;
}

template <class TPixel, unsigned int VDim>
int
ImageConverter<TPixel, VDim>
::WriteMultiple(int argc, char *argv[], int n_comp, const char *command)
{
  // Check if the argument is a printf pattern
  char buffer[1024];
  sprintf(buffer, argv[1],0);
  if (strcmp(buffer, argv[1]))
    {
    // A pattern is specified. For each image on the stack, use pattern
    for(size_t i = 0; i < m_ImageStack.size(); i+= n_comp)
      {
      WriteImage<TPixel, VDim> adapter(this);
      sprintf(buffer, argv[1], i / n_comp);
      if(n_comp == 1)
        adapter(buffer, true, i);
      else 
        adapter.WriteMultiComponent(buffer, n_comp, i);
      }
    return 1;
    }
  else
    {
    // Filenames are specified. Find out how many there are
    size_t nfiles = 0;
    for(int i = 1; i < argc; i++)
      {
      if (argv[i][0] != '-') nfiles++; else break;
      }

    // There must be files
    if (nfiles == 0)
      throw ConvertException("No files specified to %s command", command);

    if (nfiles * n_comp > m_ImageStack.size())
      throw ConvertException("Too many files specified to %s command", command);

    // Determine the starting position
    int pstart = m_ImageStack.size() - nfiles * n_comp;

    for(size_t j = pstart; j < m_ImageStack.size(); j+=n_comp)
      {
      WriteImage<TPixel, VDim> adapter(this);
      if(n_comp == 1)
        adapter(argv[j+1], true, j);
      else
        adapter.WriteMultiComponent(argv[j+1], n_comp, j);
      }

    return nfiles;
    }
}

template class ImageConverter<double, 2>;
template class ImageConverter<double, 3>;
template class ImageConverter<double, 4>;

