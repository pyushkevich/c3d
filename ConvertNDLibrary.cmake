# This is an include file to use when C3D library is to be built inside of
# another CMake project (i.e., ITK-SNAP). It is expected that the parent
# cmake file will take care of finding ITK

# Versioning information
SET(C3D_VERSION_MAJOR 1)
SET(C3D_VERSION_MINOR 2)
SET(C3D_VERSION_PATCH 0)
SET(C3D_VERSION_FULL "${C3D_VERSION_MAJOR}.${C3D_VERSION_MINOR}.${C3D_VERSION_PATCH}")

# Include directories
SET(CONVERT3D_INCLUDE_DIRS
  ${CONVERT3D_SOURCE_DIR}
  ${CONVERT3D_SOURCE_DIR}/adapters
  ${CONVERT3D_SOURCE_DIR}/external/GCv2p3
  ${CONVERT3D_SOURCE_DIR}/itkextras/
  ${CONVERT3D_SOURCE_DIR}/itkextras/VoxBoIO
  ${CONVERT3D_SOURCE_DIR}/itkextras/PovRayIO
  ${CONVERT3D_SOURCE_DIR}/itkextras/RandomForest
  ${CONVERT3D_SOURCE_DIR}/itkextras/Texture
  ${CONVERT3D_SOURCE_DIR}/utilities/doc
  ${CONVERT3D_BINARY_DIR})

IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
  ADD_DEFINITIONS(-D_SCL_SECURE_NO_WARNINGS)
  SOURCE_GROUP("Adapter Sources" REGULAR_EXPRESSION "adapters/*cxx")
  SOURCE_GROUP("Adapter Headers" REGULAR_EXPRESSION "adapters/*h")
ENDIF(WIN32)

# Markdown documentation compiled into the C code
# modified from: https://github.com/starseeker/tinyscheme-cmake/blob/master/CMakeLists.txt
# # Rather than load the init.scm file at run time,
# # with the uncertainties as to where exactly the file
# # resides, use William Ahern's hexdump to generate 
# # an embeddable version. Build our own copy of hexdump
# # to ensure consistent behavior and portability.
# # See http://25thandclement.com/~william/projects/hexdump.c.html
ADD_EXECUTABLE(markdown_to_hex utilities/hexdump.c)
set_property(TARGET markdown_to_hex APPEND PROPERTY COMPILE_DEFINITIONS "HEXDUMP_MAIN")
ADD_CUSTOM_COMMAND(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/markdown_docs.h
  COMMAND markdown_to_hex -i ${CMAKE_CURRENT_SOURCE_DIR}/doc/c3d.md > ${CMAKE_CURRENT_BINARY_DIR}/markdown_docs.h
  DEPENDS markdown_to_hex ${CMAKE_CURRENT_SOURCE_DIR}/doc/c3d.md)
ADD_CUSTOM_TARGET(markdown_docs ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/markdown_docs.h)

SET(SOURCES
  adapters/AddImages.cxx
  adapters/AlignByLandmarks.cxx
  adapters/AntiAliasImage.cxx
  adapters/ApplyMetric.cxx
  adapters/ApplyNoise.cxx
  adapters/BinaryHoleFill.cxx
  adapters/BiasFieldCorrectionN4.cxx
  adapters/BinaryImageCentroid.cxx
  adapters/BinaryMathOperation.cxx
  adapters/Canny.cxx
  adapters/ClipImageIntensity.cxx
  adapters/ComputeFFT.cxx
  adapters/ComputeMoments.cxx
  adapters/ComputeOverlaps.cxx
  adapters/ConnectedComponents.cxx
  adapters/Convolution.cxx
  adapters/CoordinateMap.cxx
  adapters/CopyTransform.cxx
  adapters/CreateImage.cxx
  adapters/CreateInterpolator.cxx
  adapters/DicomSeriesList.cxx
  adapters/ExportPatches.cxx
  adapters/ExtractRegion.cxx
  adapters/ExtractSlice.cxx
  adapters/ExtrudeSegmentation.cxx
  adapters/FillBackgroundWithNeighborhoodNoise.cxx
  adapters/FlipImage.cxx
  adapters/HaralickFeatures.cxx
  adapters/HessianEigenValues.cxx
  adapters/HessianObjectness.cxx
  adapters/HistogramMatch.cxx
  adapters/ImageERF.cxx
  adapters/ImageGradient.cxx
  adapters/ImageLaplacian.cxx
  adapters/GeneralLinearModel.cxx
  adapters/LabelOverlapMeasures.cxx
  adapters/LabelStatistics.cxx
  adapters/LandmarksToSpheres.cxx
  adapters/LaplacianSharpening.cxx
  adapters/LevelSetSegmentation.cxx
  adapters/MatchBoundingBoxes.cxx
  adapters/MathematicalMorphology.cxx
  adapters/MeanFilter.cxx
  adapters/MedianFilter.cxx
  adapters/MixtureModel.cxx
  adapters/MomentsFeatures.cxx
  adapters/MRFVote.cxx
  adapters/MultiplyImages.cxx
  adapters/NormalizeLocalWindow.cxx
  adapters/NormalizedCrossCorrelation.cxx
  adapters/OverlayLabelImage.cxx
  adapters/OtsuThreshold.cxx
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
  adapters/RetainLabels.cxx
  adapters/RFApply.cxx
  adapters/RFTrain.cxx
  adapters/SampleImage.cxx
  adapters/ScaleShiftImage.cxx
  adapters/ScalarToRGB.cxx
  adapters/SetSform.cxx
  adapters/SetOrientation.cxx
  adapters/SignedDistanceTransform.cxx
  adapters/SLICSuperVoxel.cxx
  adapters/SmoothImage.cxx
  adapters/SplitMultilabelImage.cxx
  adapters/StapleAlgorithm.cxx
  adapters/StructureTensorEigenValues.cxx
  adapters/SwapDimensions.cxx
  adapters/TestImage.cxx
  adapters/ThresholdImage.cxx
  adapters/TileImages.cxx
  adapters/TrimImage.cxx
  adapters/UnaryMathOperation.cxx
  adapters/UpdateMetadataKey.cxx
  adapters/Vote.cxx
  adapters/VoxelwiseComponentFunction.cxx
  adapters/VoxelwiseRegression.cxx
  adapters/WarpImage.cxx
  adapters/WarpLabelImage.cxx
  adapters/WeightedSum.cxx
  adapters/WeightedSumVoxelwise.cxx
  adapters/WrapDimension.cxx
  adapters/WriteImage.cxx
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx)

# Configure the version file
CONFIGURE_FILE(
  ${CONVERT3D_SOURCE_DIR}/ConvertImageVersion.cxx.in
  ${CONVERT3D_BINARY_DIR}/ConvertImageVersion.cxx @ONLY IMMEDIATE)

# Get the extra stuff compiled
SUBDIRS(${CONVERT3D_SOURCE_DIR}/itkextras)

ADD_LIBRARY(cnd_adapters ${SOURCES})

ADD_LIBRARY(cnd_maxflow 
  external/GCv2p3/GCoptimization.cpp
  external/GCv2p3/LinkedBlockList.cpp
  external/GCv2p3/graph.cpp
  external/GCv2p3/maxflow.cpp
)

ADD_LIBRARY(cnd_driver 
  ConvertImageND.cxx
  utilities/doc/Documentation.cxx)

ADD_LIBRARY(cnd_api api/ConvertAPI.cxx)

ADD_DEPENDENCIES(cnd_driver markdown_docs)

# Set up include path for the various libraries
FOREACH(target cnd_adapters cnd_driver cnd_api cnd_maxflow)
  TARGET_INCLUDE_DIRECTORIES(${target} PRIVATE ${CONVERT3D_INCLUDE_DIRS})
ENDFOREACH()

SET(C3D_LINK_LIBRARIES
  cnd_driver cnd_adapters cnd_maxflow ${ITK_LIBRARIES} ITKVoxBoIO ITKPovRayIO)
