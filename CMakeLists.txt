PROJECT(CONVERT3D)

# Nonsense for cmake 2.6 compatibility
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# Versioning information
SET(C3D_VERSION_MAJOR 1)
SET(C3D_VERSION_MINOR 0)
SET(C3D_VERSION_PATCH 0)
SET(C3D_VERSION_FULL "${C3D_VERSION_MAJOR}.${C3D_VERSION_MINOR}.${C3D_VERSION_PATCH}")

# Do we want to build the UI
OPTION(BUILD_GUI "Do you want to build the Qt-based C3D user interfaace?" OFF)
IF(BUILD_GUI)
  SUBDIRS(gui)
ENDIF(BUILD_GUI)

# Include directories
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/adapters)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/VoxBoIO)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/PovRayIO)
INCLUDE_DIRECTORIES(${CONVERT3D_SOURCE_DIR}/itkextras/N3ITK)

# Get ITK
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
  SOURCE_GROUP("Adapter Sources" REGULAR_EXPRESSION "adapters/*cxx")
  SOURCE_GROUP("Adapter Headers" REGULAR_EXPRESSION "adapters/*h")
ENDIF(WIN32)

SET(SOURCES
  adapters/AddImages.cxx
  adapters/AlignByLandmarks.cxx
  adapters/AntiAliasImage.cxx
  adapters/ApplyMetric.cxx
  adapters/BinaryHoleFill.cxx
  adapters/BiasFieldCorrection.cxx
  adapters/BiasFieldCorrectionN4.cxx
  adapters/BinaryImageCentroid.cxx
  adapters/BinaryMathOperation.cxx
  adapters/ClipImageIntensity.cxx
  adapters/ComputeFFT.cxx
  adapters/ComputeOverlaps.cxx
  adapters/ConnectedComponents.cxx
  adapters/CoordinateMap.cxx
  adapters/CopyTransform.cxx
  adapters/CreateImage.cxx
  adapters/CreateInterpolator.cxx
  adapters/ExtractRegion.cxx
  adapters/ExtractSlice.cxx
  adapters/FlipImage.cxx
  adapters/HistogramMatch.cxx
  adapters/ImageERF.cxx
  adapters/ImageLaplacian.cxx
  adapters/GeneralLinearModel.cxx
  adapters/LabelOverlapMeasures.cxx
  adapters/LabelStatistics.cxx
  adapters/LandmarksToSpheres.cxx
  adapters/LevelSetSegmentation.cxx
  adapters/MathematicalMorphology.cxx
  adapters/MixtureModel.cxx
  adapters/MRFVote.cxx
  adapters/MultiplyImages.cxx
  adapters/NormalizeLocalWindow.cxx
  adapters/NormalizedCrossCorrelation.cxx
  adapters/OverlayLabelImage.cxx
  adapters/PadImage.cxx
  adapters/PeronaMalik.cxx
  adapters/PrintImageInfo.cxx
  adapters/Rank.cxx
  adapters/ReadImage.cxx
  adapters/ReciprocalImage.cxx
  adapters/ReorderStack.cxx
  adapters/ReplaceIntensities.cxx
  adapters/ResampleImage.cxx
  adapters/ResliceImage.cxx
  adapters/SampleImage.cxx
  adapters/ScaleShiftImage.cxx
  adapters/ScalarToRGB.cxx
  adapters/SetSform.cxx
  adapters/SetOrientation.cxx
  adapters/SignedDistanceTransform.cxx
  adapters/SmoothImage.cxx
  adapters/SplitMultilabelImage.cxx
  adapters/StapleAlgorithm.cxx
  adapters/TestImage.cxx
  adapters/ThresholdImage.cxx
  adapters/TileImages.cxx
  adapters/TrimImage.cxx
  adapters/UnaryMathOperation.cxx
  adapters/UpdateMetadataKey.cxx
  adapters/Vote.cxx
  adapters/VoxelwiseRegression.cxx
  adapters/WarpImage.cxx
  adapters/WarpLabelImage.cxx
  adapters/WeightedSum.cxx
  adapters/WeightedSumVoxelwise.cxx
  adapters/WrapDimension.cxx
  adapters/WriteImage.cxx
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx)

# Pull in FFTW
FIND_LIBRARY(FFTW_LIB fftw3f)

# Configure the version file
CONFIGURE_FILE(
  ${CONVERT3D_SOURCE_DIR}/ConvertImageVersion.cxx.in
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx @ONLY IMMEDIATE)

# Get the extra stuff compiled
SUBDIRS(itkextras)
SUBDIRS(utilities)

ADD_LIBRARY(cnd_adapters ${SOURCES})
ADD_LIBRARY(cnd_driver ConvertImageND.cxx)

LINK_LIBRARIES(
  cnd_adapters cnd_driver ${ITK_LIBRARIES} ITKVoxBoIO ITKPovRayIO)

ADD_EXECUTABLE(c3d Convert3DMain.cxx)

ADD_EXECUTABLE(c2d Convert2DMain.cxx)

INSTALL_TARGETS("/bin" c3d c2d)

INCLUDE(cmake/Package.cmake)

# CDash Support
ENABLE_TESTING()
INCLUDE(CTest)
